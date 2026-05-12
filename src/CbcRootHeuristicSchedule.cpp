// Copyright (C) 2026, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcRootHeuristicSchedule.hpp"
#include "CbcModel.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicDive.hpp"
#include "CbcHeuristicDiveConfigurable.hpp"
#include "CoinTime.hpp"

#include <algorithm>
#include <mutex>
#include <thread>
#include <vector>

CbcRootHeuristicSchedule::CbcRootHeuristicSchedule(CbcModel &model)
  : model_(model)
  , stopFlag_(false)
  , maxSolutionsPhase1_(1)
  , numThreads_(0)
  , solutionsFound_(0)
{
}

CbcRootHeuristicSchedule::~CbcRootHeuristicSchedule() {}

void CbcRootHeuristicSchedule::addDefaultDivingConfigs()
{
  // v5 schedule: greedy-optimal 5-config from 3 rounds of experiments
  // Coverage: 141-144/314 instances from LP relaxation alone
  auto add = [&](const char *name, int fixCount, double rand, int seed,
                 int guideMode, double origWeight, int warmup, int freq) {
    CbcHeuristicDiveConfigurable *h = new CbcHeuristicDiveConfigurable(model_);
    h->setHeuristicName(name);
    h->setFixCount(fixCount);
    h->setRandomFactor(rand);
    h->setMinFractionality(0.10);
    h->setWeightFractionality(1.0);
    h->setTargetFractionality(0.5);
    h->setMaxIterations(500);
    // 100K simplex iters covers P99+ of successful dives in experiments
    // COIN_INT_MAX bypasses adjustHeuristics() override
    h->setMaxSimplexIterations(100000);
    h->setMaxSimplexIterationsAtRoot(100000);
    h->setSeed(seed);
    h->setGuidedObjMode(guideMode);
    h->setGuidedObjWeight(1.0);
    h->setGuidedObjOrigWeight(origWeight);
    h->setGuidedObjWarmup(warmup);
    h->setGuidedObjFrequency(freq);
    h->setWhen(2);
    model_.addHeuristic(h);
    delete h; // addHeuristic clones
  };

  // Slot 1: Blend30_r05 — best coverage (130 inst)
  add("Blend30_r05", 1, 0.5, 1, 1, 0.3, 0, 1);
  // Slot 2: W1_Blend30_f2 — best obj quality (+5 inst)
  add("W1_Blend30_f2", 1, 0.3, 1, 1, 0.3, 1, 2);
  // Slot 3: base_r03 — no guidance, complementary (+4 inst)
  add("base_r03", 1, 0.3, 1, 0, 0, 0, 1);
  // Slot 4: NoFix_Blend30 — packing/partitioning problems (+3 inst)
  add("NoFix_Blend30", 0, 0.0, 1, 1, 0.3, 0, 1);
  // Slot 5: Blend30_r05 seed=3 — diversity (+2 inst)
  add("Blend30_r05_s3", 1, 0.5, 3, 1, 0.3, 0, 1);
}

int CbcRootHeuristicSchedule::run(bool afterCuts)
{
  solutionsFound_ = 0;
  int logLevel = model_.messageHandler()->logLevel();

  // Check time/event before starting
  if (model_.maximumSecondsReached() || model_.eventHappened())
    return 0;

  // Partition heuristics into constructive vs improvement
  std::vector<CbcHeuristic *> constructive;
  std::vector<CbcHeuristic *> improvement;

  for (int i = 0; i < model_.numberHeuristics(); i++) {
    CbcHeuristic *h = model_.heuristic(i);
    if (!h)
      continue;
    // Use shouldHeurRun to respect when_ settings
    if (!h->shouldHeurRun(0))
      continue;

    switch (h->category()) {
    case HeuristicCategory::CONSTRUCTIVE:
      constructive.push_back(h);
      break;
    case HeuristicCategory::IMPROVEMENT:
    case HeuristicCategory::IMPROVEMENT_2:
      improvement.push_back(h);
      break;
    default:
      constructive.push_back(h);
      break;
    }
  }


  // Phase 1: constructive heuristics
  if (!constructive.empty()) {
    double t0 = CoinGetTimeOfDay();
    if (logLevel >= 1)
      printf("\n▶ Root heuristics%s — Phase 1: constructive (%d heuristics, %d threads)\n",
        afterCuts ? " (after cuts)" : "",
        (int)constructive.size(), numThreads_ > 0 ? numThreads_ : 1);

    int found = runPhase1(constructive);
    solutionsFound_ += found;
    double elapsed = CoinGetTimeOfDay() - t0;

    if (logLevel >= 1) {
      if (found > 0)
        printf("✔ Phase 1 — %d solution(s), best %.6g in %.3fs\n",
          found, model_.getObjValue(), elapsed);
      else
        printf("  Phase 1 — no solution (%.3fs)\n", elapsed);
    }
  }

  // Phase 2: improvement heuristics (only if we have a solution)
  if (model_.bestSolution() && !improvement.empty()
      && !model_.maximumSecondsReached() && !model_.eventHappened()) {
    double t0 = CoinGetTimeOfDay();
    double objBefore = model_.getObjValue();
    if (logLevel >= 1)
      printf("\n▶ Root heuristics%s — Phase 2: improvement (%d heuristics)\n",
        afterCuts ? " (after cuts)" : "",
        (int)improvement.size());

    int found = runPhase2(improvement);
    solutionsFound_ += found;
    double elapsed = CoinGetTimeOfDay() - t0;

    if (logLevel >= 1) {
      if (found > 0 && model_.getObjValue() < objBefore - 1.0e-5)
        printf("✔ Phase 2 — improved %.6g → %.6g in %.3fs\n",
          objBefore, model_.getObjValue(), elapsed);
      else
        printf("  Phase 2 — no improvement (%.3fs)\n", elapsed);
    }
  }

  return solutionsFound_;
}

int CbcRootHeuristicSchedule::runPhase1(
  std::vector<CbcHeuristic *> &heuristics)
{
  return runParallel(heuristics, maxSolutionsPhase1_);
}

int CbcRootHeuristicSchedule::runPhase2(
  std::vector<CbcHeuristic *> &heuristics)
{
  // For improvement phase, run all and collect best
  return runParallel(heuristics, 0); // 0 = no early stop
}

int CbcRootHeuristicSchedule::runParallel(
  std::vector<CbcHeuristic *> &heuristics, int maxSolutions)
{
  int nHeur = static_cast<int>(heuristics.size());
  if (nHeur == 0)
    return 0;

  int nThreads = numThreads_;
  if (nThreads <= 0)
    nThreads = model_.getNumberThreads();
  if (nThreads <= 0)
    nThreads = 1;

  int nCols = model_.solver()->getNumCols();
  stopFlag_.store(false);

  // Set abort flag on dive heuristics
  for (auto *h : heuristics) {
    CbcHeuristicDive *dive = dynamic_cast<CbcHeuristicDive *>(h);
    if (dive)
      dive->setAbortFlag(&stopFlag_);
  }

  // Results storage
  struct Result {
    int status = 0;
    double obj = COIN_DBL_MAX;
    std::vector<double> solution;
    double time = 0;
  };
  std::vector<Result> results(nHeur);
  std::mutex solutionMutex;
  std::atomic<int> solutionCount{0};

  // Worker function
  auto worker = [&](int idx) {
    if (maxSolutions > 0 && solutionCount.load() >= maxSolutions)
      return;
    if (model_.maximumSecondsReached() || model_.eventHappened()) {
      stopFlag_.store(true);
      return;
    }

    results[idx].solution.resize(nCols);
    results[idx].obj = model_.getCutoff();

    double t0 = CoinGetTimeOfDay();
    results[idx].status = heuristics[idx]->solution(
      results[idx].obj, results[idx].solution.data());
    results[idx].time = CoinGetTimeOfDay() - t0;

    if (results[idx].status > 0) {
      solutionCount.fetch_add(1);
      if (maxSolutions > 0 && solutionCount.load() >= maxSolutions)
        stopFlag_.store(true);
    }
  };

  // Dispatch
  if (nThreads == 1 || nHeur == 1) {
    for (int i = 0; i < nHeur; i++) {
      if (maxSolutions > 0 && solutionCount.load() >= maxSolutions)
        break;
      worker(i);
    }
  } else {
    for (int start = 0; start < nHeur; start += nThreads) {
      if (maxSolutions > 0 && solutionCount.load() >= maxSolutions)
        break;
      int end = std::min(start + nThreads, nHeur);
      std::vector<std::thread> threads;
      for (int i = start; i < end; i++)
        threads.emplace_back(worker, i);
      for (auto &t : threads)
        t.join();
    }
  }

  // Clear abort flags
  for (auto *h : heuristics) {
    CbcHeuristicDive *dive = dynamic_cast<CbcHeuristicDive *>(h);
    if (dive)
      dive->setAbortFlag(nullptr);
  }

  // Report results table (only if logLevel >= 1)
  int logLevel = model_.messageHandler()->logLevel();
  int found = 0;
  double bestObj = COIN_DBL_MAX;
  int bestIdx = -1;

  for (int i = 0; i < nHeur; i++) {
    if (results[i].status > 0) {
      if (results[i].obj < bestObj) {
        bestObj = results[i].obj;
        bestIdx = i;
      }
      found++;
    }
  }

  if (logLevel >= 1 && nHeur > 0) {
    // Print summary table — only show heuristics that took meaningful time
    bool headerPrinted = false;
    for (int i = 0; i < nHeur; i++) {
      if (results[i].time < 0.01 && results[i].status <= 0)
        continue; // skip trivially fast failures
      if (!headerPrinted) {
        printf("\n  %-24s %-22s %14s %8s\n",
          "Heuristic", "Status", "Objective", "Time(s)");
        printf("  ──────────────────────── ────────────────────── ────────────── ────────\n");
        headerPrinted = true;
      }
      char statusBuf[64];
      int statusWidth;
      if (results[i].status > 0) {
        snprintf(statusBuf, sizeof(statusBuf), "★ solution");
        statusWidth = 24; // 22 display + 2 UTF-8 bytes
      } else {
        // Show dive stats if available
        CbcHeuristicDive *dive = dynamic_cast<CbcHeuristicDive *>(heuristics[i]);
        if (dive && dive->lastDiveIterations() > 0) {
          int reason = dive->lastReasonToStop();
          const char *reasonStr = "";
          if (reason == 6) reasonStr = " (abort)";
          else if (reason / 100 > 0) reasonStr = " (prop inf)";
          else if (reason % 10 == 1) reasonStr = " (inf)";
          else if (reason % 10 == 2) reasonStr = " (simp lim)";
          else if (reason % 10 == 4) reasonStr = " (iter lim)";
          snprintf(statusBuf, sizeof(statusBuf), "%d dives, %dK lp its%s",
            dive->lastDiveIterations(), dive->lastSimplexIterations() / 1000,
            reasonStr);
        } else if (results[i].time >= 1.0) {
          snprintf(statusBuf, sizeof(statusBuf), "no solution");
        } else {
          snprintf(statusBuf, sizeof(statusBuf), "no solution");
        }
        statusWidth = 22;
      }
      if (results[i].status > 0)
        printf("  %-24s %-*s %14.4f %8.3f\n",
          heuristics[i]->heuristicName(), statusWidth, statusBuf, results[i].obj, results[i].time);
      else
        printf("  %-24s %-*s %14s %8.3f\n",
          heuristics[i]->heuristicName(), statusWidth, statusBuf, "", results[i].time);
    }
    if (headerPrinted)
      printf("\n");
  }

  // Register solutions with model (best first)
  int accepted = 0;
  if (bestIdx >= 0) {
    double before = model_.getObjValue();
    model_.setBestSolution(CBC_ROUNDING, results[bestIdx].obj,
      results[bestIdx].solution.data());
    if (model_.getObjValue() < before - 1.0e-5)
      accepted++;
  }
  // Register other solutions too
  for (int i = 0; i < nHeur; i++) {
    if (results[i].status > 0 && i != bestIdx) {
      double before = model_.getObjValue();
      model_.setBestSolution(CBC_ROUNDING, results[i].obj,
        results[i].solution.data());
      if (model_.getObjValue() < before - 1.0e-5)
        accepted++;
    }
  }

  return accepted > 0 ? accepted : found;
}
