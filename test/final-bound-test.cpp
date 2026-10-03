// Final interrupted-search bounds must reflect live nodes before cleanup,
// independently of the timed progress cadence. The binary optimum is enumerated.
// This file is licensed under the Eclipse Public License (EPL).

#include "CbcModel.hpp"
#include "CbcCompareObjective.hpp"
#include "CbcEventHandler.hpp"
#include "CbcTree.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include "OsiClpSolverInterface.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

namespace {
struct Frontier {
  bool captured = false;
  double bound = 0.0;
  double endBound = COIN_DBL_MAX;
  int endCalls = 0;
};

class RecordingTree : public CbcTree {
public:
  explicit RecordingTree(const std::shared_ptr< Frontier > &frontier)
    : frontier_(frontier)
  {
  }
  CbcTree *clone() const override
  {
    return new RecordingTree(*this);
  }
  void cleanTree(CbcModel *model, double cutoff, double &bound) override
  {
    if (cutoff == -COIN_DBL_MAX && size()) {
      frontier_->captured = true;
      frontier_->bound = std::min(getBestPossibleObjective(), model->getMinimizationObjValue());
    }
    CbcTree::cleanTree(model, cutoff, bound);
  }

private:
  std::shared_ptr< Frontier > frontier_;
};

class StopHandler : public CbcEventHandler {
public:
  StopHandler(const std::shared_ptr< Frontier > &frontier, int stopMode)
    : frontier_(frontier)
    , stopMode_(stopMode)
  {
  }
  CbcEventHandler *clone() const override
  {
    return new StopHandler(*this);
  }
  CbcAction event(CbcEvent which) override
  {
    if (which == endSearch) {
      ++frontier_->endCalls;
      frontier_->endBound = model_->getBestPossibleObjValue() * model_->solver()->getObjSense();
    }
    if (which == node && model_->getNodeCount() >= 12) {
      if (stopMode_ == 1)
        return stop;
      if (stopMode_ == 2)
        model_->setMaximumSeconds(0.0);
    }
    return noAction;
  }

private:
  std::shared_ptr< Frontier > frontier_;
  int stopMode_;
};

int run(int threads, double sense, int limit, double frequency, int stopMode = 0, int threadMode = 0)
{
  const int columns = 18, rows = 6;
  CoinPackedMatrix matrix(false, 0, 0);
  matrix.setDimensions(0, columns);
  std::vector< double > lower(columns, 0.0), upper(columns, 1.0);
  std::vector< double > objective(columns), rowLower(rows, -COIN_DBL_MAX), rowUpper(rows);
  for (int j = 0; j < columns; ++j)
    objective[j] = -sense * (13 + (j * 17 + j * j * 3) % 71);
  for (int i = 0; i < rows; ++i) {
    CoinPackedVector row;
    double total = 0.0;
    for (int j = 0; j < columns; ++j) {
      double value = 1 + (i * 19 + j * 23 + i * j * 7 + j * j) % 53;
      row.insert(j, value);
      total += value;
    }
    matrix.appendRow(row);
    rowUpper[i] = std::floor(total * 0.43);
  }
  double optimum = 0.0;
  for (unsigned int mask = 0; mask < (1U << columns); ++mask) {
    double cost = 0.0;
    bool feasible = true;
    for (int i = 0; i < rows && feasible; ++i) {
      double activity = 0.0;
      for (int j = 0; j < columns; ++j) {
        if (mask & (1U << j))
          activity += 1 + (i * 19 + j * 23 + i * j * 7 + j * j) % 53;
      }
      feasible = activity <= rowUpper[i];
    }
    if (feasible) {
      for (int j = 0; j < columns; ++j) {
        if (mask & (1U << j))
          cost += sense * objective[j];
      }
      optimum = std::min(optimum, cost);
    }
  }
  OsiClpSolverInterface solver;
  solver.loadProblem(matrix, lower.data(), upper.data(), objective.data(), rowLower.data(), rowUpper.data());
  solver.setObjSense(sense);
  for (int j = 0; j < columns; ++j)
    solver.setInteger(j);
  solver.messageHandler()->setLogLevel(0);
  solver.initialSolve();
  const double rootBound = sense * solver.getObjValue();
  CbcModel model(solver);
  model.setLogLevel(0);
  model.setNumberThreads(threads);
  model.setThreadMode(threadMode);
  model.setNumberStrong(0);
  model.setNumberBeforeTrust(0);
  CbcCompareObjective comparison;
  model.setNodeComparison(comparison);
  model.setMaximumNodes(limit);
  model.setPrintFrequency(1);
  model.setSecsPrintFrequency(frequency);
  auto frontier = std::make_shared< Frontier >();
  RecordingTree tree(frontier);
  model.passInTreeHandler(tree);
  StopHandler handler(frontier, stopMode);
  CbcBnBOutput output(stdout, false, 0);
  if (frequency == 0.0)
    handler.setOutputHandler(&output);
  model.passInEventHandler(&handler);
  model.branchAndBound();
  const double bound = sense * model.getBestPossibleObjValue();
  const bool valid = bound <= optimum + 1.0e-7;
  const bool fresh = !frontier->captured
    || (threads ? bound <= frontier->bound + 1.0e-7 && bound > rootBound + 1.0e-7
                : fabs(bound - frontier->bound) < 1.0e-7);
  const bool stopped = stopMode == 1 ? model.secondaryStatus() == 5
                                     : (stopMode == 2 ? model.isSecondsLimitReached()
                                                      : (limit <= 30 ? model.isNodeLimitReached() : model.isProvenOptimal()));
  const bool callback = limit == 10000 && !stopMode ? true : fabs(bound - frontier->endBound) < 1.0e-7;
  const bool finished = limit == 10000 && !stopMode ? fabs(bound - optimum) < 1.0e-7 : true;
  const bool frontierPresent = limit == 30 ? frontier->captured : true;
  const bool passed = valid && fresh && stopped && callback && finished && frontierPresent && frontier->endCalls == 1;
  printf("  %s: threads=%d threadMode=%d sense=%.0f limit=%d frequency=%.0f stopMode=%d nodes=%d bound=%.10g frontier=%.10g optimum=%.10g\n",
    passed ? "ok" : "FAIL", threads, threadMode, sense, limit, frequency, stopMode, model.getNodeCount(), bound, frontier->bound, optimum);
  return !passed;
}
}

int main()
{
  int failures = 0;
  {
    OsiClpSolverInterface solver;
    CbcModel model(solver);
    if (model.dealWithEventHandler(CbcEventHandler::endSearch, 0.0, NULL) != CbcEventHandler::noAction)
      ++failures;
  }
  for (double sense : { 1.0, -1.0 }) {
    for (int threads : { 0, 2 }) {
      if (threads && !CbcModel::haveMultiThreadSupport())
        continue;
      failures += run(threads, sense, 30, 1.0e9);
      failures += run(threads, sense, 30, 0.0);
      failures += run(threads, sense, 0, 1.0e9);
      failures += run(threads, sense, 10000, 1.0e9);
      if (!threads) {
        failures += run(threads, sense, 10000, 1.0e9, 1);
        failures += run(threads, sense, 10000, 1.0e9, 2);
      }
    }
    if (CbcModel::haveMultiThreadSupport()) {
      failures += run(1, sense, 30, 1.0e9);
      failures += run(4, sense, 30, 1.0e9);
      failures += run(2, sense, 30, 1.0e9, 0, 1);
    }
  }
  return failures ? 1 : 0;
}
