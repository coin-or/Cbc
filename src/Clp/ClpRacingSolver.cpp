// Copyright (C) 2026, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "ClpRacingSolver.hpp"
#include "ClpEventHandler.hpp"
#include "CoinTime.hpp"
#include <atomic>
#include <memory>
#include <vector>

// On Windows (MSVC) we use std::thread/std::mutex, which is safe there.
// On Linux/Unix we use POSIX pthreads directly to avoid embedding GCC's
// std::thread weak symbols (_M_thread_deps_never_run, _State_impl vtable)
// into static binaries, which can cause at-exit crashes on certain
// glibc/libstdc++ versions.
#ifdef _WIN32
#include <mutex>
#include <thread>
#define CLP_RACING_USE_STD_THREAD 1
#else
#include <dlfcn.h>
#include <pthread.h>
#define CLP_RACING_USE_STD_THREAD 0
#endif

// Restrict OpenBLAS to 1 thread per racing thread to avoid thread explosion.
#if !defined(_WIN32)
namespace {
inline void racing_set_openblas_threads(int n)
{
  typedef void (*fn_t)(int);
  static fn_t fn = reinterpret_cast<fn_t>(dlsym(RTLD_DEFAULT, "openblas_set_num_threads"));
  if (fn)
    fn(n);
}
} // namespace
#elif defined(CLP_USE_OPENBLAS)
extern "C" void openblas_set_num_threads(int num_threads);
namespace {
inline void racing_set_openblas_threads(int n) { openblas_set_num_threads(n); }
} // namespace
#else
namespace {
inline void racing_set_openblas_threads(int) {}
} // namespace
#endif

// ─── Racing progress handler ─────────────────────────────────────────────────
// Shared across all racing threads; prints interleaved progress rows showing
// which method is currently reporting, gated by a time frequency.
namespace {

struct RacingProgressState {
#if CLP_RACING_USE_STD_THREAD
  std::mutex mu;
#else
  pthread_mutex_t mu;
#endif
  FILE *fp = nullptr;
  double startTime = 0.0;
  double lastPrintTime = 0.0;
  double timeFreq = 5.0;
  bool headerPrinted = false;
};

class RacingEventHandler : public ClpEventHandler {
public:
  RacingEventHandler(std::shared_ptr<RacingProgressState> state,
    const std::atomic<bool> *abortFlag, const char *label)
    : ClpEventHandler()
    , state_(state)
    , abortFlag_(abortFlag)
    , label_(label) {}

  RacingEventHandler(const RacingEventHandler &rhs)
    : ClpEventHandler(rhs)
    , state_(rhs.state_)
    , abortFlag_(rhs.abortFlag_)
    , label_(rhs.label_) {}

  ClpEventHandler *clone() const override
  {
    return new RacingEventHandler(*this);
  }

  int event(Event whichEvent) override
  {
    // Check abort
    if (abortFlag_ && abortFlag_->load(std::memory_order_relaxed))
      return 0; // stop

    if (whichEvent != endOfIteration || !model_)
      return -1;

    // Time-gated progress output
    if (!state_ || !state_->fp)
      return -1;

    const double now = CoinWallclockTime();
#if CLP_RACING_USE_STD_THREAD
    std::lock_guard<std::mutex> lock(state_->mu);
    if (now - state_->lastPrintTime < state_->timeFreq)
      return -1;
#else
    pthread_mutex_lock(&state_->mu);
    if (now - state_->lastPrintTime < state_->timeFreq) {
      pthread_mutex_unlock(&state_->mu);
      return -1;
    }
#endif

    if (!state_->headerPrinted) {
      fprintf(state_->fp, "  %-12s %10s %16s %12s %8s\n",
        "Method", "Iter", "Objective", "Infeas", "Time");
      state_->headerPrinted = true;
    }

    double elapsed = now - state_->startTime;
    fprintf(state_->fp, "  %-12s %10d %16.6e %12.4e %7.1fs\n",
      label_,
      model_->numberIterations(),
      model_->objectiveValue(),
      model_->sumPrimalInfeasibilities(),
      elapsed);
    fflush(state_->fp);
    state_->lastPrintTime = now;
#if !CLP_RACING_USE_STD_THREAD
    pthread_mutex_unlock(&state_->mu);
#endif
    return -1;
  }

private:
  std::shared_ptr<RacingProgressState> state_;
  const std::atomic<bool> *abortFlag_;
  const char *label_;
};

} // namespace

ClpRacingSolver::ClpRacingSolver(ClpSimplex *model, int numThreads)
  : model_(model)
  , numThreads_(numThreads)
{
}

void ClpRacingSolver::addConfig(const ClpSolve &config)
{
  configs_.push_back(config);
}

void ClpRacingSolver::addDefaultConfigs()
{
  configs_.clear();

  // Config 0: Dual simplex with presolve (default)
  {
    ClpSolve opts;
    opts.setSolveType(ClpSolve::useDual);
    opts.setPresolveType(ClpSolve::presolveOn);
    opts.setSpecialOption(2, 1); // no signal handler (thread-safe)
    configs_.push_back(opts);
  }

  // Config 1: Primal + Idiot crash (60 passes) with presolve
  {
    ClpSolve opts;
    opts.setSolveType(ClpSolve::usePrimal);
    opts.setPresolveType(ClpSolve::presolveOn);
    opts.setSpecialOption(1, 2, 60); // idiot with 60 passes
    opts.setSpecialOption(2, 1);
    configs_.push_back(opts);
  }

  // Config 2: Primal + Sprint with presolve
  {
    ClpSolve opts;
    opts.setSolveType(ClpSolve::usePrimalorSprint);
    opts.setPresolveType(ClpSolve::presolveOn);
    opts.setSpecialOption(1, 3); // sprint
    opts.setSpecialOption(2, 1);
    configs_.push_back(opts);
  }
}

int ClpRacingSolver::solve()
{
  if (configs_.empty())
    addDefaultConfigs();

  int nConfigs = static_cast<int>(configs_.size());
  int nThreads = numThreads_ > 0 ? numThreads_ : nConfigs;
  if (nThreads > nConfigs)
    nThreads = nConfigs;

  // If only one config, just solve directly (no threading overhead)
  if (nThreads <= 1) {
    model_->initialSolve(configs_[0]);
    if (model_->status() == 0) {
      winnerIndex_ = 0;
      winnerIterations_ = model_->numberIterations();
    }
    return winnerIndex_;
  }

  std::atomic<bool> abortFlag{false};
  std::atomic<int> winner{-1};

  // Clone models for each racing thread
  std::vector<ClpSimplex *> clones(nThreads, nullptr);
  for (int i = 0; i < nThreads; i++)
    clones[i] = new ClpSimplex(*model_);

  // Set up shared progress reporting
  static const char *configLabels[] = {"Dual", "Primal+Idiot", "Sprint"};
  std::shared_ptr<RacingProgressState> progressState;
  if (model_->logLevel() > 0) {
    progressState = std::make_shared<RacingProgressState>();
#if !CLP_RACING_USE_STD_THREAD
    pthread_mutex_init(&progressState->mu, nullptr);
#endif
    progressState->fp = model_->messageHandler()
      ? model_->messageHandler()->filePointer() : stdout;
    progressState->startTime = CoinGetTimeOfDay();
    progressState->lastPrintTime = progressState->startTime;
    progressState->timeFreq = 5.0;
  }

  double startTime = CoinGetTimeOfDay();

#if CLP_RACING_USE_STD_THREAD
  std::vector<std::thread> threads;
  threads.reserve(nThreads);
  for (int i = 0; i < nThreads; i++) {
    threads.emplace_back([i, &clones, &abortFlag, &winner, &progressState, this]() {
      ClpSimplex *clone = clones[i];
      clone->setLogLevel(0);
      racing_set_openblas_threads(1);
      const char *label = (i < 3) ? configLabels[i] : "Config";
      RacingEventHandler handler(progressState, &abortFlag, label);
      clone->passInEventHandler(&handler);
      clone->initialSolve(configs_[i]);
      int st = clone->status();
      if (st == 0 || st == 1 || st == 2) {
        int expected = -1;
        if (winner.compare_exchange_strong(expected, i))
          abortFlag.store(true, std::memory_order_relaxed);
      }
    });
  }
  for (auto &t : threads)
    t.join();
#else
  // Per-thread argument struct (passed to the pthread callback)
  struct ThreadArg {
    ClpSimplex *clone;
    ClpSolve *config;
    std::atomic<bool> *abortFlag;
    std::atomic<int> *winner;
    std::shared_ptr<RacingProgressState> *progressState;
    const char *label;
    int index;
  };

  std::vector<ThreadArg> args(nThreads);
  std::vector<pthread_t> threads(nThreads);

  for (int i = 0; i < nThreads; i++) {
    args[i].clone = clones[i];
    args[i].config = &configs_[i];
    args[i].abortFlag = &abortFlag;
    args[i].winner = &winner;
    args[i].progressState = &progressState;
    args[i].label = (i < 3) ? configLabels[i] : "Config";
    args[i].index = i;

    pthread_create(&threads[i], nullptr, [](void *arg) -> void * {
      ThreadArg *a = static_cast<ThreadArg *>(arg);
      ClpSimplex *clone = a->clone;
      clone->setLogLevel(0);
      racing_set_openblas_threads(1);
      RacingEventHandler handler(*a->progressState, a->abortFlag, a->label);
      clone->passInEventHandler(&handler);
      clone->initialSolve(*a->config);
      int st = clone->status();
      if (st == 0 || st == 1 || st == 2) {
        int expected = -1;
        if (a->winner->compare_exchange_strong(expected, a->index))
          a->abortFlag->store(true, std::memory_order_relaxed);
      }
      return nullptr;
    }, &args[i]);
  }

  for (int i = 0; i < nThreads; i++)
    pthread_join(threads[i], nullptr);

  if (progressState)
    pthread_mutex_destroy(&progressState->mu);
#endif

  double endTime = CoinGetTimeOfDay();
  winnerIndex_ = winner.load();

  if (winnerIndex_ >= 0) {
    ClpSimplex *w = clones[winnerIndex_];
    winnerTime_ = endTime - startTime;
    winnerIterations_ = w->numberIterations();

    // Copy solution back to original model
    int nCols = model_->numberColumns();
    int nRows = model_->numberRows();
    CoinMemcpyN(w->primalColumnSolution(), nCols,
      model_->primalColumnSolution());
    CoinMemcpyN(w->dualColumnSolution(), nCols,
      model_->dualColumnSolution());
    CoinMemcpyN(w->primalRowSolution(), nRows,
      model_->primalRowSolution());
    CoinMemcpyN(w->dualRowSolution(), nRows,
      model_->dualRowSolution());
    CoinMemcpyN(w->statusArray(), nCols + nRows,
      model_->statusArray());
    model_->setObjectiveValue(w->objectiveValue());
    model_->setProblemStatus(w->status());
    model_->setNumberIterations(winnerIterations_);
    model_->setSecondaryStatus(w->secondaryStatus());
  }

  // Clean up clones
  for (int i = 0; i < nThreads; i++)
    delete clones[i];

  return winnerIndex_;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
