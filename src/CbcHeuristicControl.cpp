// Copyright (C) COIN-OR contributors. All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcHeuristicControl.hpp"
#include <algorithm>
#include <cstdio>

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
  if (useGap && mode != HeuristicScheduleMode::Legacy && nodeCount - lastRunNode < everyK) {
    ++skippedSchedule;
    return false;
  }
  return true;
}

bool CbcHeuristicControl::scheduled(int nodeCount, int depth, int solutionCount)
{
  if (mode == HeuristicScheduleMode::Legacy)
    return true;
  bool run;
  switch (mode) {
  case HeuristicScheduleMode::EveryKDepth:
    run = (depth % everyK) == 0;
    break;
  case HeuristicScheduleMode::EveryKNodesNoImprove: {
    const bool improved = solutionCount != lastSolutionCount;
    lastSolutionCount = solutionCount;
    // a fresh incumbent is itself a trigger and restarts the K-node window
    run = improved || nodeCount - lastRunNode >= everyK;
    break;
  }
  case HeuristicScheduleMode::EveryKNodes:
  default:
    run = nodeCount - lastRunNode >= everyK;
    break;
  }
  if (run && mode != HeuristicScheduleMode::EveryKDepth)
    lastRunNode = nodeCount;
  if (!run)
    ++skippedSchedule;
  return run;
}

void CbcHeuristicControl::finishRun(bool success, double runSeconds, int nodeCount, int solutionCount)
{
  seconds += runSeconds;
  ++runs;
  lastRunNode = nodeCount;
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
  if (window > 0 && (runs % window) == 0 && successes < minSuccess * runs) {
    // computed in double and clamped before truncating, so that long runs
    // cannot overflow the interval
    const double grown = std::min(everyK * (1.0 + decay), static_cast< double >(maxK));
    const int next = std::max(static_cast< int >(grown), everyK);
    if (next != everyK) {
      everyK = next;
      ++backoffs;
    }
  }
}

std::string CbcHeuristicControl::describe() const
{
  char buf[200];
  std::snprintf(buf, sizeof(buf), "%s K=%d decay=%g window=%d minSuccess=%g maxFail=%d timeShare=%g",
    modeName(mode), everyK, decay, window, minSuccess, maxFailures, timeShare);
  return buf;
}
