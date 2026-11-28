// Copyright (C) COIN-OR contributors. All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcHeuristicControl.hpp"
#include <algorithm>
#include <cstdio>
#include <mutex>

CbcHeuristicControl::CbcHeuristicControl()
  : scheduleState_(new ScheduleState())
{
}

CbcHeuristicControl::CbcHeuristicControl(const CbcHeuristicControl &rhs)
  : mode(rhs.mode)
  , everyK(rhs.everyK)
  , decay(rhs.decay)
  , window(rhs.window)
  , minSuccess(rhs.minSuccess)
  , maxK(rhs.maxK)
  , maxFailures(rhs.maxFailures)
  , timeShare(rhs.timeShare)
  , runs(rhs.runs)
  , successes(rhs.successes)
  , failures(rhs.failures)
  , backoffs(rhs.backoffs)
  , skippedPaused(rhs.skippedPaused)
  , skippedTime(rhs.skippedTime)
  , skippedSchedule(rhs.skippedSchedule)
  , seconds(rhs.seconds)
  , initialK(rhs.initialK)
  , scheduleState_(new ScheduleState())
  , failuresAt(rhs.failuresAt)
{
  scheduleState_->lastRunNode.store(rhs.scheduleState_->lastRunNode.load(std::memory_order_relaxed),
    std::memory_order_relaxed);
  scheduleState_->lastIncumbentCount.store(
    rhs.scheduleState_->lastIncumbentCount.load(std::memory_order_relaxed),
    std::memory_order_relaxed);
  scheduleState_->currentEveryK.store(
    rhs.scheduleState_->currentEveryK.load(std::memory_order_relaxed),
    std::memory_order_relaxed);
  std::lock_guard< std::mutex > lock(rhs.scheduleState_->backoffMutex);
  scheduleState_->runs = rhs.scheduleState_->runs;
  scheduleState_->successes = rhs.scheduleState_->successes;
}

CbcHeuristicControl &CbcHeuristicControl::operator=(const CbcHeuristicControl &rhs)
{
  if (this != &rhs) {
    std::shared_ptr< ScheduleState > newState(new ScheduleState());
    newState->lastRunNode.store(rhs.scheduleState_->lastRunNode.load(std::memory_order_relaxed),
      std::memory_order_relaxed);
    newState->lastIncumbentCount.store(
      rhs.scheduleState_->lastIncumbentCount.load(std::memory_order_relaxed),
      std::memory_order_relaxed);
    newState->currentEveryK.store(
      rhs.scheduleState_->currentEveryK.load(std::memory_order_relaxed),
      std::memory_order_relaxed);
    {
      std::lock_guard< std::mutex > lock(rhs.scheduleState_->backoffMutex);
      newState->runs = rhs.scheduleState_->runs;
      newState->successes = rhs.scheduleState_->successes;
    }
    mode = rhs.mode;
    everyK = rhs.everyK;
    decay = rhs.decay;
    window = rhs.window;
    minSuccess = rhs.minSuccess;
    maxK = rhs.maxK;
    maxFailures = rhs.maxFailures;
    timeShare = rhs.timeShare;
    runs = rhs.runs;
    successes = rhs.successes;
    failures = rhs.failures;
    backoffs = rhs.backoffs;
    skippedPaused = rhs.skippedPaused;
    skippedTime = rhs.skippedTime;
    skippedSchedule = rhs.skippedSchedule;
    seconds = rhs.seconds;
    initialK = rhs.initialK;
    failuresAt = rhs.failuresAt;
    scheduleState_ = newState;
  }
  return *this;
}

void CbcHeuristicControl::setScheduleMode(HeuristicScheduleMode newMode, int k)
{
  mode = newMode;
  setEveryK(k);
  scheduleState_->lastRunNode.store(-1000000, std::memory_order_relaxed);
  scheduleState_->lastIncumbentCount.store(0, std::memory_order_relaxed);
  std::lock_guard< std::mutex > lock(scheduleState_->backoffMutex);
  scheduleState_->runs = 0;
  scheduleState_->successes = 0;
}

void CbcHeuristicControl::setEveryK(int k)
{
  everyK = k > 0 ? k : 1;
  initialK = everyK;
  scheduleState_->currentEveryK.store(everyK, std::memory_order_relaxed);
}

int CbcHeuristicControl::currentEveryK() const
{
  const int current = scheduleState_->currentEveryK.load(std::memory_order_relaxed);
  return current > 0 ? current : everyK;
}

void CbcHeuristicControl::shareScheduleStateWith(CbcHeuristicControl &other)
{
  scheduleState_ = other.scheduleState_;
}

const char *CbcHeuristicControl::modeName(HeuristicScheduleMode m)
{
  switch (m) {
  case HeuristicScheduleMode::EveryKDepth:
    return "depth";
  case HeuristicScheduleMode::EveryKNodes:
    return "nodes";
  case HeuristicScheduleMode::EveryKNodesNoImprove:
    return "noimprove";
  default:
    return "legacy";
  }
}

bool CbcHeuristicControl::allows(int nodeCount, double elapsed, int solutionCount, bool useGap)
{
  if (paused(solutionCount)) {
    ++skippedPaused;
    return false;
  }
  if (timeShare > 0.0 && runs > 0 && seconds > timeShare * elapsed) {
    ++skippedTime;
    return false;
  }
  if (useGap && mode != HeuristicScheduleMode::Legacy) {
    const int k = currentEveryK();
    if (nodeCount - scheduleState_->lastRunNode.load(std::memory_order_relaxed) < k) {
      ++skippedSchedule;
      return false;
    }
  }
  return true;
}

bool CbcHeuristicControl::scheduled(int nodeCount, int depth, int incumbentCount)
{
  if (mode == HeuristicScheduleMode::Legacy)
    return true;
  const int k = currentEveryK();
  bool run;
  switch (mode) {
  case HeuristicScheduleMode::EveryKDepth:
    run = (depth % k) == 0;
    break;
  case HeuristicScheduleMode::EveryKNodesNoImprove: {
    int lastIncumbent = scheduleState_->lastIncumbentCount.load(std::memory_order_relaxed);
    bool improved = false;
    while (incumbentCount > lastIncumbent) {
      if (scheduleState_->lastIncumbentCount.compare_exchange_weak(lastIncumbent, incumbentCount,
            std::memory_order_relaxed)) {
        improved = true;
        break;
      }
    }
    // a fresh incumbent is itself a trigger and restarts the K-node window
    run = improved;
    break;
  }
  case HeuristicScheduleMode::EveryKNodes:
  default:
    run = false;
    break;
  }
  if (mode != HeuristicScheduleMode::EveryKDepth) {
    if (run) {
      int lastNode = scheduleState_->lastRunNode.load(std::memory_order_relaxed);
      while (nodeCount > lastNode
        && !scheduleState_->lastRunNode.compare_exchange_weak(lastNode, nodeCount,
          std::memory_order_relaxed)) {
      }
    } else {
      int lastNode = scheduleState_->lastRunNode.load(std::memory_order_relaxed);
      while (nodeCount - lastNode >= k) {
        if (scheduleState_->lastRunNode.compare_exchange_weak(lastNode, nodeCount,
              std::memory_order_relaxed)) {
          run = true;
          break;
        }
      }
    }
  }
  if (!run)
    ++skippedSchedule;
  return run;
}

void CbcHeuristicControl::finishRun(bool success, double runSeconds, int nodeCount, int solutionCount)
{
  seconds += runSeconds;
  ++runs;
  if (success) {
    ++successes;
    failures = 0;
  } else {
    // failures only count against the incumbent they were tried on
    if (failuresAt != solutionCount)
      failures = 0;
    ++failures;
    failuresAt = solutionCount;
  }
  int lastNode = scheduleState_->lastRunNode.load(std::memory_order_relaxed);
  while (nodeCount > lastNode
    && !scheduleState_->lastRunNode.compare_exchange_weak(lastNode, nodeCount,
      std::memory_order_relaxed)) {
  }
  std::lock_guard< std::mutex > lock(scheduleState_->backoffMutex);
  ++scheduleState_->runs;
  if (success)
    ++scheduleState_->successes;
  if (window > 0 && (scheduleState_->runs % window) == 0
    && scheduleState_->successes < minSuccess * scheduleState_->runs) {
    // computed in double and clamped before truncating, so that long runs
    // cannot overflow the interval
    const int currentK = currentEveryK();
    const double grown = std::min(currentK * (1.0 + decay),
      static_cast< double >(maxK));
    const int next = std::max(static_cast< int >(grown), currentK);
    if (next != currentK) {
      scheduleState_->currentEveryK.store(next, std::memory_order_relaxed);
      ++backoffs;
    }
  }
}

std::string CbcHeuristicControl::describe() const
{
  char buf[200];
  std::snprintf(buf, sizeof(buf), "%s K=%d decay=%g window=%d minSuccess=%g maxFail=%d timeShare=%g",
    modeName(mode), currentEveryK(), decay, window, minSuccess, maxFailures, timeShare);
  return buf;
}
