// Copyright (C) 2024 MIPster contributors.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
//
// Integrates the Feasibility Jump heuristic (MIT © 2022 SINTEF) into CBC/MIPster.
// Reference: Leitner, Fischetti, Toth (2023), "Feasibility Jump"
//   https://link.springer.com/article/10.1007/s12532-023-00234-8

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "CbcHeuristicFeasibilityJump.hpp"
#include "CbcModel.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"

// Suppress warnings in the third-party header.
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#endif
#include "feasibilityjump.hh"
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

// ---------------------------------------------------------------------------
// Constructors / destructor / clone
// ---------------------------------------------------------------------------

CbcHeuristicFeasibilityJump::CbcHeuristicFeasibilityJump()
  : CbcHeuristic()
{
  setHeuristicName("FeasibilityJump");
}

CbcHeuristicFeasibilityJump::CbcHeuristicFeasibilityJump(CbcModel &model)
  : CbcHeuristic(model)
{
  setHeuristicName("FeasibilityJump");
}

CbcHeuristicFeasibilityJump::CbcHeuristicFeasibilityJump(const CbcHeuristicFeasibilityJump &rhs)
  : CbcHeuristic(rhs)
  , relaxContinuous_(rhs.relaxContinuous_)
  , seed_(rhs.seed_)
  , weightUpdateDecay_(rhs.weightUpdateDecay_)
  , maxEffort_(rhs.maxEffort_)
  , effortMultiplier_(rhs.effortMultiplier_)
  , stallMultiplier_(rhs.stallMultiplier_)
  , minDepth_(rhs.minDepth_)
  , maxSolutions_(rhs.maxSolutions_)
  , feasibilityTolerance_(rhs.feasibilityTolerance_)
  , integerTolerance_(rhs.integerTolerance_)
{
}

CbcHeuristicFeasibilityJump::~CbcHeuristicFeasibilityJump() {}

CbcHeuristic *CbcHeuristicFeasibilityJump::clone() const
{
  return new CbcHeuristicFeasibilityJump(*this);
}

CbcHeuristicFeasibilityJump &CbcHeuristicFeasibilityJump::operator=(const CbcHeuristicFeasibilityJump &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    relaxContinuous_ = rhs.relaxContinuous_;
    seed_ = rhs.seed_;
    weightUpdateDecay_ = rhs.weightUpdateDecay_;
    maxEffort_ = rhs.maxEffort_;
    effortMultiplier_ = rhs.effortMultiplier_;
    stallMultiplier_ = rhs.stallMultiplier_;
    minDepth_ = rhs.minDepth_;
    maxSolutions_ = rhs.maxSolutions_;
    feasibilityTolerance_ = rhs.feasibilityTolerance_;
    integerTolerance_ = rhs.integerTolerance_;
  }
  return *this;
}

bool CbcHeuristicFeasibilityJump::shouldHeurRun(int whereFrom)
{
  if (whereFrom == 4 && minDepth_ > 0) {
    // Tree call: allow if deep enough (bypass complex distance/frequency logic)
    numCouldRun_++;
    return true;
  }
  return CbcHeuristic::shouldHeurRun(whereFrom);
}

void CbcHeuristicFeasibilityJump::resetModel(CbcModel *model)
{
  setModel(model);
}

void CbcHeuristicFeasibilityJump::setModel(CbcModel *model)
{
  model_ = model;
}

// ---------------------------------------------------------------------------
// Main heuristic entry point
// ---------------------------------------------------------------------------

int CbcHeuristicFeasibilityJump::solution(double &objectiveValue,
  double *newSolution)
{
  // Depth-based control: run at root (depth 0) always, and at tree nodes
  // only when current depth >= minDepth_ (more variables fixed by branching).
  int depth = model_->currentDepth();
  int nodeCount = model_->getNodeCount();

  if (nodeCount == 0) {
    // Root node: use standard shouldHeurRun check
    if (!shouldHeurRun(0))
      return 0;
  } else {
    // Tree node: the external shouldHeurRun(4) already passed.
    // Just check our depth requirement.
    if (minDepth_ <= 0)
      return 0; // disabled in tree
    if (depth < minDepth_)
      return 0; // not deep enough
  }

  OsiSolverInterface *solver = model_->solver();
  if (!solver)
    return 0;

  const int numCols = solver->getNumCols();
  const int numRows = solver->getNumRows();
  if (numCols == 0 || numRows == 0)
    return 0;

  // Skip if no integer variables.
  bool hasIntegers = false;
  for (int j = 0; j < numCols; ++j) {
    if (solver->isInteger(j)) {
      hasIntegers = true;
      break;
    }
  }
  if (!hasIntegers)
    return 0;

  // ------------------------------------------------------------------
  // Build FeasibilityJumpSolver from the current LP data.
  // ------------------------------------------------------------------
  FeasibilityJumpSolver fj(seed_, 0, weightUpdateDecay_,
    feasibilityTolerance_, feasibilityTolerance_);

  const double *colLower = solver->getColLower();
  const double *colUpper = solver->getColUpper();
  const double *objCoeff = solver->getObjCoefficients();

  for (int j = 0; j < numCols; ++j) {
    VarType vt = solver->isInteger(j) ? VarType::Integer : VarType::Continuous;
    fj.addVar(vt, colLower[j], colUpper[j], objCoeff[j]);
  }

  const CoinPackedMatrix *matrixByRow = solver->getMatrixByRow();
  const double *elements = matrixByRow->getElements();
  const int *colIdxs = matrixByRow->getIndices();
  const CoinBigIndex *rowStarts = matrixByRow->getVectorStarts();
  const int *rowLengths = matrixByRow->getVectorLengths();

  const char *rowSense = solver->getRowSense();
  const double *rhs = solver->getRightHandSide();
  const double *rowRange = solver->getRowRange();

  int relaxFlag = relaxContinuous_ ? 1 : 0;

  for (int i = 0; i < numRows; ++i) {
    int len = rowLengths[i];
    // OsiSolverInterface uses non-const pointers in addConstraint signature.
    std::vector< int > idxBuf(colIdxs + rowStarts[i], colIdxs + rowStarts[i] + len);
    std::vector< double > elBuf(elements + rowStarts[i], elements + rowStarts[i] + len);

    char sense = rowSense[i];

    if (sense == 'E') {
      fj.addConstraint(RowType::Equal, rhs[i], len, idxBuf.data(), elBuf.data(), relaxFlag);
    } else if (sense == 'L') {
      fj.addConstraint(RowType::Lte, rhs[i], len, idxBuf.data(), elBuf.data(), relaxFlag);
    } else if (sense == 'G') {
      fj.addConstraint(RowType::Gte, rhs[i], len, idxBuf.data(), elBuf.data(), relaxFlag);
    } else if (sense == 'R') {
      // Range row: lb <= ax <= ub, i.e. rhs[i]-rowRange[i] <= ax <= rhs[i]
      fj.addConstraint(RowType::Lte, rhs[i], len, idxBuf.data(), elBuf.data(), relaxFlag);
      fj.addConstraint(RowType::Gte, rhs[i] - rowRange[i], len, idxBuf.data(), elBuf.data(), relaxFlag);
    }
    // 'N' (free row) — ignore
  }

  // ------------------------------------------------------------------
  // Use LP relaxation solution as initial point; round integers.
  // ------------------------------------------------------------------
  const double *lpSol = solver->getColSolution();
  std::vector< double > initialValues(numCols);
  for (int j = 0; j < numCols; ++j) {
    if (solver->isInteger(j))
      initialValues[j] = std::floor(lpSol[j] + 0.5);
    else
      initialValues[j] = lpSol[j];
    // Clamp to bounds.
    initialValues[j] = std::max(colLower[j], std::min(colUpper[j], initialValues[j]));
  }

  // ------------------------------------------------------------------
  // Check CBC global time limit before starting.
  // ------------------------------------------------------------------
  double cbcMaxSecs = model_->getMaximumSeconds();
  if (cbcMaxSecs > 0.0 && model_->getCurrentSeconds() >= cbcMaxSecs)
    return 0;

  // ------------------------------------------------------------------
  // Compute effective effort budget and stall limit.
  // In the tree (depth > 0), use a fraction of the root budget to keep
  // per-node overhead low.
  // ------------------------------------------------------------------
  const int64_t nnz = matrixByRow->getNumElements();
  int64_t effectiveMaxEffort = maxEffort_;
  if (effectiveMaxEffort <= 0) {
    // NNZ-scaled: nnz * effortMultiplier_
    effectiveMaxEffort = nnz * static_cast<int64_t>(effortMultiplier_);
  }
  // In tree: use 1/4 of root effort (tighter bounds → fewer iterations needed)
  if (depth > 0)
    effectiveMaxEffort = std::max(int64_t(1000000), effectiveMaxEffort / 4);
  int64_t stallLimit = 0;
  if (stallMultiplier_ > 0) {
    stallLimit = nnz * static_cast<int64_t>(stallMultiplier_);
  }

  // ------------------------------------------------------------------
  // Run FeasibilityJump; collect best feasible solution via callback.
  // Termination conditions (whichever fires first):
  //   1. status.totalEffort >= effectiveMaxEffort (iteration budget)
  //   2. status.effortSinceLastImprovement > stallLimit (stall)
  //   3. solutionsFound >= maxSolutions_   (enough solutions found)
  //   4. CBC global time limit exceeded    (hard wall)
  // ------------------------------------------------------------------
  std::vector< double > bestSol;
  double bestObj = DBL_MAX;
  double cutoff = model_->getCutoff();
  int solutionsFound = 0;

  // Progress reporting: print a row every printInterval seconds (wall-clock,
  // used only for display — the stopping criterion remains iteration-based).
  const double printInterval = 1.0; // seconds
  double lastPrintTime = model_->getCurrentSeconds();
  double fjStartTime = lastPrintTime;
  bool headerPrinted = false;
  int logLevel = model_->messageHandler()->logLevel();
  FILE *fp = model_->messageHandler()->filePointer();
  if (!fp)
    fp = stdout;

  auto printHeader = [&]() {
    if (!headerPrinted && logLevel >= 1) {
      fprintf(fp, "\n  FeasibilityJump progress:\n");
      fprintf(fp, "  %-9s  %-9s  %-9s  %11s\n",
        "Effort(M)", "Stall(M)", "Solutions", "BestObj");
      fprintf(fp, "  %-9s  %-9s  %-9s  %11s\n",
        "---------", "---------", "---------", "-----------");
      fflush(fp);
      headerPrinted = true;
    }
  };

  auto printProgressRow = [&](FJStatus &status, bool isSolution) {
    printHeader();
    double effortM = static_cast< double >(status.totalEffort) / 1.0e6;
    double stallM = static_cast< double >(status.effortSinceLastImprovement) / 1.0e6;
    char solBuf[24];
    if (bestSol.empty())
      snprintf(solBuf, sizeof(solBuf), "%11s", "-");
    else
      snprintf(solBuf, sizeof(solBuf), "%11.2f", bestObj);
    fprintf(fp, "  %c %7.3f  %9.3f  %9d  %s\n",
      isSolution ? '*' : ' ', effortM, stallM, solutionsFound, solBuf);
    fflush(fp);
    lastPrintTime = model_->getCurrentSeconds();
  };

  int64_t callbackCount = 0;

  auto callback = [&](FJStatus status) -> CallbackControlFlow {
    double now = model_->getCurrentSeconds();

    // Hard time-limit wall: respect CBC's overall time budget.
    if (cbcMaxSecs > 0.0 && now >= cbcMaxSecs)
      return CallbackControlFlow::Terminate;

    // Deterministic iteration budget.
    if (status.totalEffort >= effectiveMaxEffort)
      return CallbackControlFlow::Terminate;

    // Stall-based early termination.
    if (stallLimit > 0 && status.effortSinceLastImprovement > stallLimit)
      return CallbackControlFlow::Terminate;

    ++callbackCount;

    if (status.solution != nullptr) {
      double obj = status.solutionObjectiveValue * solver->getObjSense();
      if (obj < bestObj && obj < cutoff) {
        bestObj = obj;
        bestSol.assign(status.solution, status.solution + status.numVars);
      }
      ++solutionsFound;
      if (logLevel >= 1)
        printProgressRow(status, true);
      if (solutionsFound >= maxSolutions_)
        return CallbackControlFlow::Terminate;
    } else if (logLevel >= 1 && (now - lastPrintTime) >= printInterval) {
      printProgressRow(status, false);
    }

    return CallbackControlFlow::Continue;
  };

  fj.solve(initialValues.data(), callback);

  // Print final summary row if header was printed.
  if (headerPrinted && logLevel >= 1) {
    fprintf(fp, "  (done — %lld callbacks, %.2fs)\n",
      (long long)callbackCount, model_->getCurrentSeconds() - fjStartTime);
    fflush(fp);
  }

  if (bestSol.empty())
    return 0;

  // Copy into CBC's output buffer and report.
  std::copy(bestSol.begin(), bestSol.end(), newSolution);
  objectiveValue = bestObj;
  return 1;
}
