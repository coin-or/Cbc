// Copyright (C) COIN-OR contributors. All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcHeuristicControl.hpp"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <thread>
#include <vector>

int main()
{
  const int workerCount = 8;
  const int nodeLimit = 10000;
  CbcHeuristicControl base;
  base.setScheduleMode(HeuristicScheduleMode::EveryKNodesNoImprove, 100);

  std::vector< CbcHeuristicControl > controls(workerCount, base);
  for (int i = 0; i < workerCount; ++i)
    controls[i].shareScheduleStateWith(base);

  std::atomic< int > runs(0);
  std::vector< std::thread > workers;
  for (int worker = 0; worker < workerCount; ++worker) {
    workers.push_back(std::thread([&, worker]() {
      for (int node = 0; node <= nodeLimit; ++node) {
        const int nodeSnapshot = std::max(0, node - 3 * worker);
        const int incumbentSnapshot = std::max(0, node / 250 - worker % 4);
        if (controls[worker].scheduled(nodeSnapshot, 0, incumbentSnapshot))
          ++runs;
      }
    }));
  }
  for (std::vector< std::thread >::iterator worker = workers.begin();
    worker != workers.end(); ++worker)
    worker->join();

  // There are at most 100 node intervals and 40 incumbent changes. A shared
  // gate should admit roughly one run per trigger, not one per worker.
  if (runs <= 0 || runs > 145) {
    std::fprintf(stderr, "FAIL threaded schedule gate admitted %d runs (expected 1..145)\n",
      runs.load());
    return 1;
  }

  CbcHeuristicControl backoffBase;
  backoffBase.decay = 1.0;
  backoffBase.window = workerCount;
  backoffBase.minSuccess = 1.0;
  backoffBase.setScheduleMode(HeuristicScheduleMode::EveryKNodesNoImprove, 100);
  std::vector< CbcHeuristicControl > backoffControls(workerCount, backoffBase);
  for (int i = 0; i < workerCount; ++i)
    backoffControls[i].shareScheduleStateWith(backoffBase);
  workers.clear();
  for (int worker = 0; worker < workerCount; ++worker) {
    workers.push_back(std::thread([&, worker]() {
      backoffControls[worker].finishRun(false, 0.0, worker, 0);
    }));
  }
  for (std::vector< std::thread >::iterator worker = workers.begin();
    worker != workers.end(); ++worker)
    worker->join();

  if (backoffControls[0].currentEveryK() != 200) {
    std::fprintf(stderr, "FAIL threaded backoff interval is %d (expected 200)\n",
      backoffControls[0].currentEveryK());
    return 1;
  }
  std::printf("PASS threaded schedule gate admitted %d runs\n", runs.load());
  std::printf("PASS shared backoff interval is %d\n", backoffControls[0].currentEveryK());
  return 0;
}
