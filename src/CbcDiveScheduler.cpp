// Copyright (C) 2026, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcDiveScheduler.hpp"
#include "CbcModel.hpp"
#include "CbcHeuristicDive.hpp"
#include "CbcHeuristicDiveConfigurable.hpp"
#include "CoinTime.hpp"

#include <algorithm>
#include <climits>
#include <cstdio>
#include <thread>
#include <mutex>
#include <vector>

CbcDiveScheduler::CbcDiveScheduler(CbcModel &model)
  : model_(model)
  , numThreads_(0)
  , stopOnFirst_(true)
  , bestObj_(COIN_DBL_MAX)
  , winnerIdx_(-1)
  , timeUsed_(0.0)
  , solutionFound_(false)
{
}

CbcDiveScheduler::~CbcDiveScheduler()
{
  for (auto *h : schedule_)
    delete h;
}

void CbcDiveScheduler::addHeuristic(CbcHeuristicDive *heuristic)
{
  schedule_.push_back(heuristic);
}

void CbcDiveScheduler::addDefaultSchedule()
{
  const int NL = INT_MAX >> 3;
  double timeLimit = model_.getMaximumSeconds() - model_.getCurrentSeconds();
  if (timeLimit < 1.0)
    timeLimit = 600.0;

  auto make = [&](const char *name, double pctFix, int fixCount,
                  double minFrac, double wFrac, double wLocks,
                  double randFactor, double targetFrac, int maxIter) {
    auto *h = new CbcHeuristicDiveConfigurable(model_);
    h->setHeuristicName(name);
    h->setWeightFractionality(wFrac);
    h->setWeightLocks(wLocks);
    h->setMinFractionality(minFrac);
    h->setPercentageToFix(pctFix);
    h->setFixCount(fixCount);
    h->setRandomFactor(randFactor);
    h->setTargetFractionality(targetFrac);
    h->setMaxIterations(maxIter);
    h->setMaxSimplexIterations(NL);
    h->setMaxSimplexIterationsAtRoot(NL);
    h->setMaxTime(timeLimit);
    h->setWhen(2);
    schedule_.push_back(h);
  };

  //       name          pct   fc  minF  wF   wL   rand  tgt   maxIter
  // v5 schedule: greedy-optimal 5-config from 3 rounds of experiments
  make("Blend30_r05",   0.0,  1, 0.10, 1.0, 0.0, 0.5, 0.5, 1000);
  make("W1_Blend30_f2", 0.0,  1, 0.10, 1.0, 0.0, 0.3, 0.5, 1000);
  make("base_r03",      0.0,  1, 0.10, 1.0, 0.0, 0.3, 0.5, 1000);
  make("NoFix_Blend30", 0.0,  0, 0.10, 1.0, 0.0, 0.0, 0.5, 1000);
  make("Blend30_r05_s3",0.0,  1, 0.10, 1.0, 0.0, 0.5, 0.5, 1000);

  // Set seeds
  schedule_[0]->setSeed(1);
  schedule_[1]->setSeed(1);
  schedule_[2]->setSeed(1);
  schedule_[3]->setSeed(1);
  schedule_[4]->setSeed(3);

  // Guided objective: Blend30 (origWeight=0.3) on slots 0,1,3,4
  schedule_[0]->setGuidedObjMode(1);
  schedule_[0]->setGuidedObjOrigWeight(0.3);
  schedule_[1]->setGuidedObjMode(1);
  schedule_[1]->setGuidedObjOrigWeight(0.3);
  schedule_[1]->setGuidedObjWarmup(1);
  schedule_[1]->setGuidedObjFrequency(2);
  // slot 2: no guidance (base_r03)
  schedule_[3]->setGuidedObjMode(1);
  schedule_[3]->setGuidedObjOrigWeight(0.3);
  schedule_[4]->setGuidedObjMode(1);
  schedule_[4]->setGuidedObjOrigWeight(0.3);
}

const char *CbcDiveScheduler::winnerName() const
{
  if (winnerIdx_ >= 0 && winnerIdx_ < static_cast<int>(schedule_.size()))
    return schedule_[winnerIdx_]->heuristicName();
  return "(none)";
}

int CbcDiveScheduler::run()
{
  if (schedule_.empty())
    return 0;

  // Determine thread count
  int nThreads = numThreads_;
  if (nThreads == 0) {
    // Auto: use model's thread count, default to 1
    nThreads = model_.getNumberThreads();
    if (nThreads <= 0)
      nThreads = 1;
  }

  solutionFound_.store(false);
  bestObj_ = model_.getCutoff();
  winnerIdx_ = -1;

  double t0 = CoinGetTimeOfDay();
  int result = (nThreads > 1) ? runParallel() : runSequential();
  timeUsed_ = CoinGetTimeOfDay() - t0;

  return result;
}

int CbcDiveScheduler::runSequential()
{
  int nCols = model_.getNumCols();
  std::vector<double> solution(nCols);

  for (int i = 0; i < static_cast<int>(schedule_.size()); i++) {
    if (model_.maximumSecondsReached())
      break;

    double objValue = bestObj_;
    int status = schedule_[i]->solution(objValue, solution.data());

    if (status > 0 && objValue < bestObj_) {
      bestObj_ = objValue;
      winnerIdx_ = i;
      // Store solution - don't call setBestSolution here as model may not
      // have full state. Caller is responsible for registering the solution.
      bestSolution_.assign(solution.begin(), solution.end());

      if (model_.messageHandler()->logLevel() >= 1)
        printf("  DiveScheduler: %s found solution %.6f\n",
          schedule_[i]->heuristicName(), objValue);

      if (stopOnFirst_)
        break;
    }
  }

  return (winnerIdx_ >= 0) ? 1 : 0;
}

int CbcDiveScheduler::runParallel()
{
  int nCols = model_.getNumCols();
  int nHeur = static_cast<int>(schedule_.size());

  struct Result {
    int status = 0;
    double obj = COIN_DBL_MAX;
    std::vector<double> solution;
  };

  std::vector<Result> results(nHeur);
  std::mutex mtx;

  auto worker = [&](int idx) {
    // Check if another thread already found a solution
    if (stopOnFirst_ && solutionFound_.load())
      return;

    // Set abort flag so dive exits if another thread finds solution
    schedule_[idx]->setAbortFlag(&solutionFound_);

    results[idx].solution.resize(nCols);
    results[idx].obj = bestObj_;

    int status = schedule_[idx]->solution(results[idx].obj,
      results[idx].solution.data());

    if (status > 0) {
      results[idx].status = 1;
      solutionFound_.store(true);
    }
  };

  // Launch threads
  int nThreads = numThreads_;
  if (nThreads == 0)
    nThreads = model_.getNumberThreads();
  if (nThreads <= 0)
    nThreads = static_cast<int>(schedule_.size());

  // Run in batches if more heuristics than threads
  for (int batch = 0; batch < nHeur; batch += nThreads) {
    if (stopOnFirst_ && solutionFound_.load())
      break;

    int batchEnd = std::min(batch + nThreads, nHeur);
    std::vector<std::thread> threads;

    for (int i = batch; i < batchEnd; i++)
      threads.emplace_back(worker, i);

    for (auto &t : threads)
      t.join();
  }

  // Find best result
  for (int i = 0; i < nHeur; i++) {
    if (results[i].status > 0 && results[i].obj < bestObj_) {
      bestObj_ = results[i].obj;
      winnerIdx_ = i;
      bestSolution_ = std::move(results[i].solution);
    }
  }

  if (winnerIdx_ >= 0 && model_.messageHandler()->logLevel() >= 1)
    printf("  DiveScheduler: %s found solution %.6f (parallel)\n",
      winnerName(), bestObj_);

  return (winnerIdx_ >= 0) ? 1 : 0;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
