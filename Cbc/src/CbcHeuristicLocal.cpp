/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcBranchActual.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"

// Default Constructor
CbcHeuristicLocal::CbcHeuristicLocal()
  : CbcHeuristic()
{
  numberSolutions_ = 0;
  swap_ = 0;
  used_ = NULL;
  lastRunDeep_ = -1000000;
  switches_ |= 16; // needs a new solution
}

// Constructor with model - assumed before cuts

CbcHeuristicLocal::CbcHeuristicLocal(CbcModel &model)
  : CbcHeuristic(model)
{
  numberSolutions_ = 0;
  swap_ = 0;
  lastRunDeep_ = -1000000;
  switches_ |= 16; // needs a new solution
  // Get a copy of original matrix
  assert(model.solver());
  if (model.solver()->getNumRows()) {
    matrix_ = *model.solver()->getMatrixByCol();
  }
  int numberColumns = model.solver()->getNumCols();
  used_ = new int[numberColumns];
  memset(used_, 0, numberColumns * sizeof(int));
}

// Destructor
CbcHeuristicLocal::~CbcHeuristicLocal()
{
  delete[] used_;
}

// Clone
CbcHeuristic *
CbcHeuristicLocal::clone() const
{
  return new CbcHeuristicLocal(*this);
}
// Create C++ lines to get to current state
void CbcHeuristicLocal::generateCpp(FILE *fp)
{
  CbcHeuristicLocal other;
  fprintf(fp, "0#include \"CbcHeuristicLocal.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicLocal heuristicLocal(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "heuristicLocal");
  if (swap_ != other.swap_)
    fprintf(fp, "3  heuristicLocal.setSearchType(%d);\n", swap_);
  else
    fprintf(fp, "4  heuristicLocal.setSearchType(%d);\n", swap_);
  fprintf(fp, "3  cbcModel->addHeuristic(&heuristicLocal);\n");
}

// Copy constructor
CbcHeuristicLocal::CbcHeuristicLocal(const CbcHeuristicLocal &rhs)
  : CbcHeuristic(rhs)
  , matrix_(rhs.matrix_)
  , numberSolutions_(rhs.numberSolutions_)
  , swap_(rhs.swap_)
{
  if (model_ && rhs.used_) {
    int numberColumns = model_->solver()->getNumCols();
    used_ = CoinCopyOfArray(rhs.used_, numberColumns);
  } else {
    used_ = NULL;
  }
}

// Assignment operator
CbcHeuristicLocal &
CbcHeuristicLocal::operator=(const CbcHeuristicLocal &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    matrix_ = rhs.matrix_;
    numberSolutions_ = rhs.numberSolutions_;
    swap_ = rhs.swap_;
    delete[] used_;
    if (model_ && rhs.used_) {
      int numberColumns = model_->solver()->getNumCols();
      used_ = CoinCopyOfArray(rhs.used_, numberColumns);
    } else {
      used_ = NULL;
    }
  }
  return *this;
}

// Resets stuff if model changes
void CbcHeuristicLocal::resetModel(CbcModel * /*model*/)
{
  //CbcHeuristic::resetModel(model);
  delete[] used_;
  if (model_ && used_) {
    int numberColumns = model_->solver()->getNumCols();
    used_ = new int[numberColumns];
    memset(used_, 0, numberColumns * sizeof(int));
  } else {
    used_ = NULL;
  }
}
/*
  Run a mini-BaB search after fixing all variables not marked as used by
  solution(). (See comments there for semantics.)

  Return values are:
    1: smallBranchAndBound found a solution
    0: everything else

  The degree of overload as return codes from smallBranchAndBound are folded
  into 0 is such that it's impossible to distinguish return codes that really
  require attention from a simple `nothing of interest'.
*/
// This version fixes stuff and does IP
int CbcHeuristicLocal::solutionFix(double &objectiveValue,
  double *newSolution,
  const int * /*keep*/)
{
  /*
  If when is set to off (0), or set to root (1) and we're not at the root,
  return. If this heuristic discovered the current solution, don't continue.
*/

  numCouldRun_++;
  // See if to do
  if (!when() || (when() == 1 && model_->phase() != 1))
    return 0; // switched off
  // Don't do if it was this heuristic which found solution!
  if (this == model_->lastHeuristic())
    return 0;
  /*
  Load up a new solver with the solution.

  Why continuousSolver(), as opposed to solver()?
*/
  OsiSolverInterface *newSolver = model_->continuousSolver()->clone();
  const double *colLower = newSolver->getColLower();
  //const double * colUpper = newSolver->getColUpper();

  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  /*
  The net effect here is that anything that hasn't moved from its lower bound
  will be fixed at lower bound.

  See comments in solution() w.r.t. asymmetric treatment of upper and lower
  bounds.
*/

  int i;
  int nFix = 0;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (!isHeuristicInteger(newSolver, iColumn))
      continue;
    const OsiObject *object = model_->object(i);
    // get original bounds
    double originalLower;
    double originalUpper;
    getIntegerInformation(object, originalLower, originalUpper);
    newSolver->setColLower(iColumn, CoinMax(colLower[iColumn], originalLower));
    if (!used_[iColumn]) {
      newSolver->setColUpper(iColumn, colLower[iColumn]);
      nFix++;
    }
  }
  /*
  Try a `small' branch-and-bound search. The notion here is that we've fixed a
  lot of variables and reduced the amount of `free' problem to a point where a
  small BaB search will suffice to fully explore the remaining problem. This
  routine will execute integer presolve, then call branchAndBound to do the
  actual search.
*/
  int returnCode = 0;
#ifdef CLP_INVESTIGATE2
  printf("Fixing %d out of %d (%d continuous)\n",
    nFix, numberIntegers, newSolver->getNumCols() - numberIntegers);
#endif
  if (nFix * 10 <= numberIntegers) {
    // see if we can fix more
    int *which = new int[2 * (numberIntegers - nFix)];
    int *sort = which + (numberIntegers - nFix);
    int n = 0;
    for (i = 0; i < numberIntegers; i++) {
      int iColumn = integerVariable[i];
      if (!isHeuristicInteger(newSolver, iColumn))
        continue;
      if (used_[iColumn]) {
        which[n] = iColumn;
        sort[n++] = used_[iColumn];
      }
    }
    CoinSort_2(sort, sort + n, which);
    // only half fixed in total
    n = CoinMin(n, numberIntegers / 2 - nFix);
    int allow = CoinMax(numberSolutions_ - 2, sort[0]);
    int nFix2 = 0;
    for (i = 0; i < n; i++) {
      int iColumn = integerVariable[i];
      if (!isHeuristicInteger(newSolver, iColumn))
        continue;
      if (used_[iColumn] <= allow) {
        newSolver->setColUpper(iColumn, colLower[iColumn]);
        nFix2++;
      } else {
        break;
      }
    }
    delete[] which;
    nFix += nFix2;
#ifdef CLP_INVESTIGATE2
    printf("Number fixed increased from %d to %d\n",
      nFix - nFix2, nFix);
#endif
  }
  if (nFix * 10 > numberIntegers) {
    returnCode = smallBranchAndBound(newSolver, numberNodes_, newSolution, objectiveValue,
      objectiveValue, "CbcHeuristicLocal");
    /*
  -2 is return due to user event, and -1 is overloaded with what look to be
  two contradictory meanings.
*/
    if (returnCode < 0) {
      returnCode = 0; // returned on size
      int numberColumns = newSolver->getNumCols();
      int numberContinuous = numberColumns - numberIntegers;
      if (numberContinuous > 2 * numberIntegers && nFix * 10 < numberColumns) {
#define LOCAL_FIX_CONTINUOUS
#ifdef LOCAL_FIX_CONTINUOUS
        //const double * colUpper = newSolver->getColUpper();
        const double *colLower = newSolver->getColLower();
        int nAtLb = 0;
        //double sumDj=0.0;
        const double *dj = newSolver->getReducedCost();
        double direction = newSolver->getObjSense();
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (!isHeuristicInteger(newSolver, iColumn)) {
            if (!used_[iColumn]) {
              //double djValue = dj[iColumn]*direction;
              nAtLb++;
              //sumDj += djValue;
            }
          }
        }
        if (nAtLb) {
          // fix some continuous
          double *sort = new double[nAtLb];
          int *which = new int[nAtLb];
          //double threshold = CoinMax((0.01*sumDj)/static_cast<double>(nAtLb),1.0e-6);
          int nFix2 = 0;
          for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (!isHeuristicInteger(newSolver, iColumn)) {
              if (!used_[iColumn]) {
                double djValue = dj[iColumn] * direction;
                if (djValue > 1.0e-6) {
                  sort[nFix2] = -djValue;
                  which[nFix2++] = iColumn;
                }
              }
            }
          }
          CoinSort_2(sort, sort + nFix2, which);
          int divisor = 2;
          nFix2 = CoinMin(nFix2, (numberColumns - nFix) / divisor);
          for (int i = 0; i < nFix2; i++) {
            int iColumn = which[i];
            newSolver->setColUpper(iColumn, colLower[iColumn]);
          }
          delete[] sort;
          delete[] which;
#ifdef CLP_INVESTIGATE2
          printf("%d integers have zero value, and %d continuous fixed at lb\n",
            nFix, nFix2);
#endif
          returnCode = smallBranchAndBound(newSolver,
            numberNodes_, newSolution,
            objectiveValue,
            objectiveValue, "CbcHeuristicLocal");
          if (returnCode < 0)
            returnCode = 0; // returned on size
        }
#endif
      }
    }
  }
  /*
  If the result is complete exploration with a solution (3) or proven
  infeasibility (2), we could generate a cut (the AI folks would call it a
  nogood) to prevent us from going down this route in the future.
*/
  if ((returnCode & 2) != 0) {
    // could add cut
    returnCode &= ~2;
  }

  delete newSolver;
  return returnCode;
}
/*
  First tries setting a variable to better value.  If feasible then
  tries setting others.  If not feasible then tries swaps
  Returns 1 if solution, 0 if not 
  The main body of this routine implements an O((q^2)/2) brute force search
  around the current solution, for q = number of integer variables. Call this
  the inc/dec heuristic.  For each integer variable x<i>, first decrement the
  value. Then, for integer variables x<i+1>, ..., x<q-1>, try increment and
  decrement. If one of these permutations produces a better solution,
  remember it.  Then repeat, with x<i> incremented. If we find a better
  solution, update our notion of current solution and continue.

  The net effect is a greedy walk: As each improving pair is found, the
  current solution is updated and the search continues from this updated
  solution.

  Way down at the end, we call solutionFix, which will create a drastically
  restricted problem based on variables marked as used, then do mini-BaC on
  the restricted problem. This can occur even if we don't try the inc/dec
  heuristic. This would be more obvious if the inc/dec heuristic were broken
  out as a separate routine and solutionFix had a name that reflected where
  it was headed.

  The return code of 0 is grossly overloaded, because it maps to a return
  code of 0 from solutionFix, which is itself grossly overloaded. See
  comments in solutionFix and in CbcHeuristic::smallBranchAndBound.
  */
int CbcHeuristicLocal::solution(double &solutionValue,
  double *betterSolution)
{
  /*
  Execute only if a new solution has been discovered since the last time we
  were called.
*/

  numCouldRun_++;
  // See if frequency kills off idea
  int swap = swap_ % 100;
  int skip = swap_ / 100;
  int nodeCount = model_->getNodeCount();
  if (nodeCount < lastRunDeep_ + skip && nodeCount != lastRunDeep_ + 1)
    return 0;
  if (numberSolutions_ == model_->getSolutionCount() && (numberSolutions_ == howOftenShallow_ || nodeCount < lastRunDeep_ + 2 * skip))
    return 0;
  howOftenShallow_ = numberSolutions_;
  numberSolutions_ = model_->getSolutionCount();
  if (nodeCount < lastRunDeep_ + skip)
    return 0;
#ifdef HEURISTIC_INFORM
  printf("Entering heuristic %s - nRuns %d numCould %d when %d\n",
    heuristicName(), numRuns_, numCouldRun_, when_);
#endif
  lastRunDeep_ = nodeCount;
  howOftenShallow_ = numberSolutions_;

  if ((swap % 10) == 2) {
    // try merge
    return solutionFix(solutionValue, betterSolution, NULL);
  }
  /*
  Exclude long (column), thin (row) systems.

  Given the n^2 nature of the search, more than 100,000 columns could get
  expensive. But I don't yet see the rationale for the second part of the
  condition (cols > 10*rows). And cost is proportional to number of integer
  variables --- shouldn't we use that?

  Why wait until we have more than one solution?
*/
  if ((model_->getNumCols() > 100000 && model_->getNumCols() > 10 * model_->getNumRows()) || numberSolutions_ <= 1)
    return 0; // probably not worth it
  // worth trying

  OsiSolverInterface *solver = model_->solver();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  const double *solution = model_->bestSolution();
  /*
  Shouldn't this test be redundant if we've already checked that
  numberSolutions_ > 1? Stronger: shouldn't this be an assertion?
*/
  if (!solution)
    return 0; // No solution found yet
  const double *objective = solver->getObjCoefficients();
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);

  int numberRows = matrix_.getNumRows();

  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();

  int i;
  double direction = solver->getObjSense();
  double newSolutionValue = model_->getObjValue() * direction;
  int returnCode = 0;
  numRuns_++;
  // Column copy
  const double *element = matrix_.getElements();
  const int *row = matrix_.getIndices();
  const CoinBigIndex *columnStart = matrix_.getVectorStarts();
  const int *columnLength = matrix_.getVectorLengths();

  // Get solution array for heuristic solution
  int numberColumns = solver->getNumCols();
  double *newSolution = new double[numberColumns];
  memcpy(newSolution, solution, numberColumns * sizeof(double));
#ifdef LOCAL_FIX_CONTINUOUS
  // mark continuous used
  const double *columnLower = solver->getColLower();
  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (!isHeuristicInteger(solver, iColumn)) {
      if (solution[iColumn] > columnLower[iColumn] + 1.0e-8)
        used_[iColumn] = numberSolutions_;
    }
  }
#endif

  // way is 1 if down possible, 2 if up possible, 3 if both possible
  char *way = new char[numberIntegers];
  // corrected costs
  double *cost = new double[numberIntegers];
  // for array to mark infeasible rows after iColumn branch
  char *mark = new char[numberRows];
  memset(mark, 0, numberRows);
  // space to save values so we don't introduce rounding errors
  double *save = new double[numberRows];
  /*
  Force variables within their original bounds, then to the nearest integer.
  Overall, we seem to be prepared to cope with noninteger bounds. Is this
  necessary? Seems like we'd be better off to force the bounds to integrality
  as part of preprocessing.  More generally, why do we need to do this? This
  solution should have been cleaned and checked when it was accepted as a
  solution!

  Once the value is set, decide whether we can move up or down.

  The only place that used_ is used is in solutionFix; if a variable is not
  flagged as used, it will be fixed (at lower bound). Why the asymmetric
  treatment? This makes some sense for binary variables (for which there are
  only two options). But for general integer variables, why not make a similar
  test against the original upper bound?
*/

  // clean solution
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (!isHeuristicInteger(solver, iColumn))
      continue;
    const OsiObject *object = model_->object(i);
    // get original bounds
    double originalLower;
    double originalUpper;
    getIntegerInformation(object, originalLower, originalUpper);
    double value = newSolution[iColumn];
    if (value < originalLower) {
      value = originalLower;
      newSolution[iColumn] = value;
    } else if (value > originalUpper) {
      value = originalUpper;
      newSolution[iColumn] = value;
    }
    double nearest = floor(value + 0.5);
    //assert(fabs(value-nearest)<10.0*primalTolerance);
    value = nearest;
    newSolution[iColumn] = nearest;
    // if away from lower bound mark that fact
    if (nearest > originalLower) {
      used_[iColumn] = numberSolutions_;
    }
    cost[i] = direction * objective[iColumn];
    /*
  Given previous computation we're checking that value is at least 1 away
  from the original bounds.
*/
    int iway = 0;

    if (value > originalLower + 0.5)
      iway = 1;
    if (value < originalUpper - 0.5)
      iway |= 2;
    way[i] = static_cast< char >(iway);
  }
  /*
  Calculate lhs of each constraint for groomed solution.
*/
  // get row activities
  double *rowActivity = new double[numberRows];
  memset(rowActivity, 0, numberRows * sizeof(double));

  for (i = 0; i < numberColumns; i++) {
    CoinBigIndex j;
    double value = newSolution[i];
    if (value) {
      for (j = columnStart[i];
           j < columnStart[i] + columnLength[i]; j++) {
        int iRow = row[j];
        rowActivity[iRow] += value * element[j];
      }
    }
  }
  /*
  Check that constraints are satisfied. For small infeasibility, force the
  activity within bound. Again, why is this necessary if the current solution
  was accepted as a valid solution?

  Why are we scanning past the first unacceptable constraint?
*/
  // check was feasible - if not adjust (cleaning may move)
  // if very infeasible then give up
  bool tryHeuristic = true;
  for (i = 0; i < numberRows; i++) {
    if (rowActivity[i] < rowLower[i]) {
      if (rowActivity[i] < rowLower[i] - 10.0 * primalTolerance)
        tryHeuristic = false;
      rowActivity[i] = rowLower[i];
    } else if (rowActivity[i] > rowUpper[i]) {
      if (rowActivity[i] < rowUpper[i] + 10.0 * primalTolerance)
        tryHeuristic = false;
      rowActivity[i] = rowUpper[i];
    }
  }
  /*
  This bit of code is not quite totally redundant: it'll bail at 10,000
  instead of 100,000. Potentially we can do a lot of work to get here, only
  to abandon it.
*/
  // Switch off if may take too long
  if (model_->getNumCols() > 10000 && model_->getNumCols() > 10 * model_->getNumRows() && swap < 10)
    tryHeuristic = false;
  /*
  Try the inc/dec heuristic?
*/
  if (tryHeuristic) {

    // total change in objective
    double totalChange = 0.0;
    // local best change in objective
    double bestChange = 0.0;
    // maybe just do 1000
    int maxIntegers = numberIntegers;
    // stop if too many goes
    int maxTries = COIN_INT_MAX;
    // integerVariable may be randomized copy!
    int *integerVariable = CoinCopyOfArray(model_->integerVariable(), numberIntegers);
#ifdef COIN_HAS_CLP
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(model_->solver());
    if (clpSolver) {
      // take out some integers
      int nn = numberIntegers;
      numberIntegers = 0;
      for (int i = 0; i < nn; i++) {
        int iColumn = integerVariable[i];
        if (clpSolver->isHeuristicInteger(iColumn))
          integerVariable[numberIntegers++] = iColumn;
      }
    }
#endif
    if (swap > 9 && numberIntegers > 500) {
      int type = swap / 10;
      if (type == 1) {
        // reduce
        maxIntegers = CoinMin(1000, numberIntegers);
      } else if (type == 2) {
        // reduce even more
        maxTries = 100000;
        maxIntegers = CoinMin(500, numberIntegers);
      } else if (type > 2) {
        assert(type < 10);
        int totals[7] = { 1000, 500, 100, 50, 50, 50, 50 };
        maxIntegers = CoinMin(totals[type - 3], numberIntegers);
        double *weight = new double[numberIntegers];
        for (int i = 0; i < numberIntegers; i++) {
          weight[i] = model_->randomNumberGenerator()->randomDouble();
        }
        CoinSort_2(weight, weight + numberIntegers, integerVariable);
        delete[] weight;
      }
    }
    /*
  Outer loop to walk integer variables. Call the current variable x<i>. At the
  end of this loop, bestChange will contain the best (negative) change in the
  objective for any single pair.

  The trouble is, we're limited to monotonically increasing improvement.
  Suppose we discover an improvement of 10 for some pair. If, later in the
  search, we discover an improvement of 9 for some other pair, we will not use
  it. That seems wasteful.
*/

    for (i = 0; i < numberIntegers; i++) {
      int iColumn = integerVariable[i];
      bestChange = 0.0;
      int endInner = CoinMin(numberIntegers, i + maxIntegers);

      double objectiveCoefficient = cost[i];
      int k;
      CoinBigIndex j;
      int goodK = -1;
      int wayK = -1, wayI = -1;
      /*
  Try decrementing x<i>.
*/
      if ((way[i] & 1) != 0) {
        int numberInfeasible = 0;
        /*
  Adjust row activities where x<i> has a nonzero coefficient. Save the old
  values for restoration. Mark any rows that become infeasible as a result
  of the decrement.
*/
        // save row activities and adjust
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          save[iRow] = rowActivity[iRow];
          rowActivity[iRow] -= element[j];
          if (rowActivity[iRow] < rowLower[iRow] - primalTolerance || rowActivity[iRow] > rowUpper[iRow] + primalTolerance) {
            // mark row
            mark[iRow] = 1;
            numberInfeasible++;
          }
        }
        /*
  Run through the remaining integer variables. Try increment and decrement on
  each one. If the potential objective change is better than anything we've
  seen so far, do a full evaluation of x<k> in that direction.  If we can
  repair all infeasibilities introduced by pushing x<i> down, we have a
  winner. Remember the best variable, and the direction for x<i> and x<k>.
*/
        // try down
        for (k = i + 1; k < endInner; k++) {
          if (!maxTries)
            break;
          maxTries--;
          if ((way[k] & 1) != 0) {
            // try down
            if (-objectiveCoefficient - cost[k] < bestChange) {
              // see if feasible down
              bool good = true;
              int numberMarked = 0;
              int kColumn = integerVariable[k];
              for (j = columnStart[kColumn];
                   j < columnStart[kColumn] + columnLength[kColumn]; j++) {
                int iRow = row[j];
                double newValue = rowActivity[iRow] - element[j];
                if (newValue < rowLower[iRow] - primalTolerance || newValue > rowUpper[iRow] + primalTolerance) {
                  good = false;
                  break;
                } else if (mark[iRow]) {
                  // made feasible
                  numberMarked++;
                }
              }
              if (good && numberMarked == numberInfeasible) {
                // better solution
                goodK = k;
                wayK = -1;
                wayI = -1;
                bestChange = -objectiveCoefficient - cost[k];
              }
            }
          }
          if ((way[k] & 2) != 0) {
            // try up
            if (-objectiveCoefficient + cost[k] < bestChange) {
              // see if feasible up
              bool good = true;
              int numberMarked = 0;
              int kColumn = integerVariable[k];
              for (j = columnStart[kColumn];
                   j < columnStart[kColumn] + columnLength[kColumn]; j++) {
                int iRow = row[j];
                double newValue = rowActivity[iRow] + element[j];
                if (newValue < rowLower[iRow] - primalTolerance || newValue > rowUpper[iRow] + primalTolerance) {
                  good = false;
                  break;
                } else if (mark[iRow]) {
                  // made feasible
                  numberMarked++;
                }
              }
              if (good && numberMarked == numberInfeasible) {
                // better solution
                goodK = k;
                wayK = 1;
                wayI = -1;
                bestChange = -objectiveCoefficient + cost[k];
              }
            }
          }
        }
        /*
  Remove effect of decrementing x<i> by restoring original lhs values.
*/
        // restore row activities
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          rowActivity[iRow] = save[iRow];
          mark[iRow] = 0;
        }
      }
      /*
  Try to increment x<i>. Actions as for decrement.
*/
      if ((way[i] & 2) != 0) {
        int numberInfeasible = 0;
        // save row activities and adjust
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          save[iRow] = rowActivity[iRow];
          rowActivity[iRow] += element[j];
          if (rowActivity[iRow] < rowLower[iRow] - primalTolerance || rowActivity[iRow] > rowUpper[iRow] + primalTolerance) {
            // mark row
            mark[iRow] = 1;
            numberInfeasible++;
          }
        }
        // try up
        for (k = i + 1; k < endInner; k++) {
          if (!maxTries)
            break;
          if ((way[k] & 1) != 0) {
            // try down
            if (objectiveCoefficient - cost[k] < bestChange) {
              // see if feasible down
              bool good = true;
              int numberMarked = 0;
              int kColumn = integerVariable[k];
              for (j = columnStart[kColumn];
                   j < columnStart[kColumn] + columnLength[kColumn]; j++) {
                int iRow = row[j];
                double newValue = rowActivity[iRow] - element[j];
                if (newValue < rowLower[iRow] - primalTolerance || newValue > rowUpper[iRow] + primalTolerance) {
                  good = false;
                  break;
                } else if (mark[iRow]) {
                  // made feasible
                  numberMarked++;
                }
              }
              if (good && numberMarked == numberInfeasible) {
                // better solution
                goodK = k;
                wayK = -1;
                wayI = 1;
                bestChange = objectiveCoefficient - cost[k];
              }
            }
          }
          if ((way[k] & 2) != 0) {
            // try up
            if (objectiveCoefficient + cost[k] < bestChange) {
              // see if feasible up
              bool good = true;
              int numberMarked = 0;
              int kColumn = integerVariable[k];
              for (j = columnStart[kColumn];
                   j < columnStart[kColumn] + columnLength[kColumn]; j++) {
                int iRow = row[j];
                double newValue = rowActivity[iRow] + element[j];
                if (newValue < rowLower[iRow] - primalTolerance || newValue > rowUpper[iRow] + primalTolerance) {
                  good = false;
                  break;
                } else if (mark[iRow]) {
                  // made feasible
                  numberMarked++;
                }
              }
              if (good && numberMarked == numberInfeasible) {
                // better solution
                goodK = k;
                wayK = 1;
                wayI = 1;
                bestChange = objectiveCoefficient + cost[k];
              }
            }
          }
        }
        // restore row activities
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          rowActivity[iRow] = save[iRow];
          mark[iRow] = 0;
        }
      }
      /*
  We've found a pair x<i> and x<k> which produce a better solution. Update our
  notion of current solution to match.

  Why does this not update newSolutionValue?
*/
      if (goodK >= 0) {
        // we found something - update solution
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          rowActivity[iRow] += wayI * element[j];
        }
        newSolution[iColumn] += wayI;
        int kColumn = integerVariable[goodK];
        for (j = columnStart[kColumn];
             j < columnStart[kColumn] + columnLength[kColumn]; j++) {
          int iRow = row[j];
          rowActivity[iRow] += wayK * element[j];
        }
        newSolution[kColumn] += wayK;
        /*
  Adjust motion range for x<k>. We may have banged up against a bound with that
  last move.
*/
        // See if k can go further ?
        const OsiObject *object = model_->object(goodK);
        // get original bounds
        double originalLower;
        double originalUpper;
        getIntegerInformation(object, originalLower, originalUpper);

        double value = newSolution[kColumn];
        int iway = 0;

        if (value > originalLower + 0.5)
          iway = 1;
        if (value < originalUpper - 0.5)
          iway |= 2;
        way[goodK] = static_cast< char >(iway);
        totalChange += bestChange;
      }
    }
    /*
  End of loop to try increment/decrement of integer variables.

  newSolutionValue does not necessarily match the current newSolution, and
  bestChange simply reflects the best single change. Still, that's sufficient
  to indicate that there's been at least one change. Check that we really do
  have a valid solution.
*/
    if (totalChange + newSolutionValue < solutionValue) {
      // paranoid check
      memset(rowActivity, 0, numberRows * sizeof(double));

      for (i = 0; i < numberColumns; i++) {
        CoinBigIndex j;
        double value = newSolution[i];
        if (value) {
          for (j = columnStart[i];
               j < columnStart[i] + columnLength[i]; j++) {
            int iRow = row[j];
            rowActivity[iRow] += value * element[j];
          }
        }
      }
      int numberBad = 0;
#ifdef COIN_DETAIL
      double sumBad = 0.0;
#endif
      // check was approximately feasible
      for (i = 0; i < numberRows; i++) {
        if (rowActivity[i] < rowLower[i]) {
#ifdef COIN_DETAIL
          sumBad += rowLower[i] - rowActivity[i];
#endif
          if (rowActivity[i] < rowLower[i] - 10.0 * primalTolerance)
            numberBad++;
        } else if (rowActivity[i] > rowUpper[i]) {
#ifdef COIN_DETAIL
          sumBad += rowUpper[i] - rowActivity[i];
#endif
          if (rowActivity[i] > rowUpper[i] + 10.0 * primalTolerance)
            numberBad++;
        }
      }
      if (!numberBad) {
        for (i = 0; i < numberIntegers; i++) {
          int iColumn = integerVariable[i];
          const OsiObject *object = model_->object(i);
          // get original bounds
          double originalLower;
          double originalUpper;
          getIntegerInformation(object, originalLower, originalUpper);

          double value = newSolution[iColumn];
          // if away from lower bound mark that fact
          if (value > originalLower) {
            used_[iColumn] = numberSolutions_;
          }
        }
        /*
  Copy the solution to the array returned to the client. Grab a basis from
  the solver (which, if it exists, is almost certainly infeasible, but it
  should be ok for a dual start). The value returned as solutionValue is
  conservative because of handling of newSolutionValue and bestChange, as
  described above.
*/
        // new solution
        memcpy(betterSolution, newSolution, numberColumns * sizeof(double));
        CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(solver->getWarmStart());
        if (basis) {
          model_->setBestSolutionBasis(*basis);
          delete basis;
        }
        returnCode = 1;
        solutionValue = newSolutionValue + bestChange;
      } else {
        // bad solution - should not happen so debug if see message
        COIN_DETAIL_PRINT(printf("Local search got bad solution with %d infeasibilities summing to %g\n",
          numberBad, sumBad));
      }
    }
    // This is just a copy!
    delete[] integerVariable;
  }
  /*
  We're done. Clean up.
*/
  delete[] newSolution;
  delete[] rowActivity;
  delete[] way;
  delete[] cost;
  delete[] save;
  delete[] mark;
  /*
  Do we want to try swapping values between solutions?
  swap_ is set elsewhere; it's not adjusted during heuristic execution.

  Again, redundant test. We shouldn't be here if numberSolutions_ = 1.
*/
  if (numberSolutions_ > 1 && (swap % 10) == 1) {
    // try merge
    int returnCode2 = solutionFix(solutionValue, betterSolution, NULL);
    if (returnCode2)
      returnCode = 1;
  }
  return returnCode;
}
// update model
void CbcHeuristicLocal::setModel(CbcModel *model)
{
  model_ = model;
  // Get a copy of original matrix
  assert(model_->solver());
  if (model_->solver()->getNumRows()) {
    matrix_ = *model_->solver()->getMatrixByCol();
  }
  delete[] used_;
  int numberColumns = model->solver()->getNumCols();
  used_ = new int[numberColumns];
  memset(used_, 0, numberColumns * sizeof(int));
}

// Default Constructor
CbcHeuristicProximity::CbcHeuristicProximity()
  : CbcHeuristic()
{
  increment_ = 0.01;
  feasibilityPump_ = NULL;
  numberSolutions_ = 0;
  used_ = NULL;
  lastRunDeep_ = -1000000;
  switches_ |= 16; // needs a new solution
}

// Constructor with model - assumed before cuts

CbcHeuristicProximity::CbcHeuristicProximity(CbcModel &model)
  : CbcHeuristic(model)
{
  increment_ = 0.01;
  feasibilityPump_ = NULL;
  numberSolutions_ = 0;
  lastRunDeep_ = -1000000;
  switches_ |= 16; // needs a new solution
  int numberColumns = model.solver()->getNumCols();
  used_ = new int[numberColumns];
  memset(used_, 0, numberColumns * sizeof(int));
}

// Destructor
CbcHeuristicProximity::~CbcHeuristicProximity()
{
  delete feasibilityPump_;
  delete[] used_;
}

// Clone
CbcHeuristic *
CbcHeuristicProximity::clone() const
{
  return new CbcHeuristicProximity(*this);
}
// Create C++ lines to get to current state
void CbcHeuristicProximity::generateCpp(FILE *fp)
{
  CbcHeuristicProximity other;
  fprintf(fp, "0#include \"CbcHeuristicProximity.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicProximity heuristicProximity(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "heuristicProximity");
  fprintf(fp, "3  cbcModel->addHeuristic(&heuristicProximity);\n");
}

// Copy constructor
CbcHeuristicProximity::CbcHeuristicProximity(const CbcHeuristicProximity &rhs)
  : CbcHeuristic(rhs)
  , numberSolutions_(rhs.numberSolutions_)
{
  increment_ = rhs.increment_;
  feasibilityPump_ = NULL;
  if (model_ && rhs.used_) {
    int numberColumns = model_->solver()->getNumCols();
    used_ = CoinCopyOfArray(rhs.used_, numberColumns);
    if (rhs.feasibilityPump_)
      feasibilityPump_ = new CbcHeuristicFPump(*rhs.feasibilityPump_);
  } else {
    used_ = NULL;
  }
}

// Assignment operator
CbcHeuristicProximity &
CbcHeuristicProximity::operator=(const CbcHeuristicProximity &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    increment_ = rhs.increment_;
    numberSolutions_ = rhs.numberSolutions_;
    delete[] used_;
    delete feasibilityPump_;
    feasibilityPump_ = NULL;
    if (model_ && rhs.used_) {
      int numberColumns = model_->solver()->getNumCols();
      used_ = CoinCopyOfArray(rhs.used_, numberColumns);
      if (rhs.feasibilityPump_)
        feasibilityPump_ = new CbcHeuristicFPump(*rhs.feasibilityPump_);
    } else {
      used_ = NULL;
    }
  }
  return *this;
}

// Resets stuff if model changes
void CbcHeuristicProximity::resetModel(CbcModel * /*model*/)
{
  //CbcHeuristic::resetModel(model);
  delete[] used_;
  if (model_ && used_) {
    int numberColumns = model_->solver()->getNumCols();
    used_ = new int[numberColumns];
    memset(used_, 0, numberColumns * sizeof(int));
  } else {
    used_ = NULL;
  }
}
/*
  Run a mini-BaB search after changing objective

  Return values are:
    1: smallBranchAndBound found a solution
    0: everything else

  The degree of overload as return codes from smallBranchAndBound are folded
  into 0 is such that it's impossible to distinguish return codes that really
  require attention from a simple `nothing of interest'.
*/
int CbcHeuristicProximity::solution(double &solutionValue,
  double *betterSolution)
{
  if (feasibilityPumpOptions_ == -3 && numCouldRun_ == 0 && !feasibilityPump_) {
    // clone feasibility pump
    for (int i = 0; i < model_->numberHeuristics(); i++) {
      const CbcHeuristicFPump *pump = dynamic_cast< const CbcHeuristicFPump * >(model_->heuristic(i));
      if (pump) {
        feasibilityPump_ = new CbcHeuristicFPump(*pump);
        break;
      }
    }
  }
  /*
  Execute only if a new solution has been discovered since the last time we
  were called.
*/

  numCouldRun_++;
  int nodeCount = model_->getNodeCount();
  if (numberSolutions_ == model_->getSolutionCount())
    return 0;
  if (!model_->bestSolution())
    return 0; // odd - because in parallel mode
  numberSolutions_ = model_->getSolutionCount();
  lastRunDeep_ = nodeCount;
  numRuns_++;
  //howOftenShallow_ = numberSolutions_;

  /*
  Load up a new solver with the solution.

  Why continuousSolver(), as opposed to solver()?
*/
  OsiSolverInterface *newSolver = model_->continuousSolver()->clone();
  int numberColumns = newSolver->getNumCols();
  double *obj = CoinCopyOfArray(newSolver->getObjCoefficients(), numberColumns);
  int *indices = new int[numberColumns];
  int n = 0;
  for (int i = 0; i < numberColumns; i++) {
    if (obj[i]) {
      indices[n] = i;
      obj[n++] = obj[i];
    }
  }
  double cutoff = model_->getCutoff();
  assert(cutoff < 1.0e20);
  if (model_->getCutoffIncrement() < 1.0e-4) {
    cutoff -= increment_;
  }
  double offset;
  newSolver->getDblParam(OsiObjOffset, offset);
  newSolver->setDblParam(OsiObjOffset, 0.0);
  newSolver->addRow(n, indices, obj, -COIN_DBL_MAX, cutoff + offset);
  delete[] indices;
  memset(obj, 0, numberColumns * sizeof(double));
  newSolver->setDblParam(OsiDualObjectiveLimit, 1.0e20);
  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  const double *solutionIn = model_->bestSolution();
  for (int i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (!isHeuristicInteger(newSolver, iColumn))
      continue;
    if (fabs(solutionIn[iColumn]) < 1.0e-5)
      obj[iColumn] = 1.0;
    else if (fabs(solutionIn[iColumn] - 1.0) < 1.0e-5)
      obj[iColumn] = -1.0;
  }
  newSolver->setObjective(obj);
  delete[] obj;
  //newSolver->writeMps("xxxx");
  int maxSolutions = model_->getMaximumSolutions();
  model_->setMaximumSolutions(1);
  bool pumpAdded = false;
  if (feasibilityPumpOptions_ == -3 && feasibilityPump_) {
    // add back feasibility pump
    pumpAdded = true;
    for (int i = 0; i < model_->numberHeuristics(); i++) {
      const CbcHeuristicFPump *pump = dynamic_cast< const CbcHeuristicFPump * >(model_->heuristic(i));
      if (pump) {
        pumpAdded = false;
        break;
      }
    }
    if (pumpAdded)
      model_->addHeuristic(feasibilityPump_);
  }
  int returnCode = smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue,
    1.0e20, "CbcHeuristicProximity");
  if (pumpAdded) {
    // take off feasibility pump
    int lastHeuristic = model_->numberHeuristics() - 1;
    model_->setNumberHeuristics(lastHeuristic);
    delete model_->heuristic(lastHeuristic);
  }
  model_->setMaximumSolutions(maxSolutions);
  /*
  -2 is return due to user event, and -1 is overloaded with what look to be
  two contradictory meanings.
*/
  if (returnCode < 0) {
    returnCode = 0;
  }
  /*
  If the result is complete exploration with a solution (3) or proven
  infeasibility (2), we could generate a cut (the AI folks would call it a
  nogood) to prevent us from going down this route in the future.
*/
  if ((returnCode & 2) != 0) {
    // could add cut
    returnCode &= ~2;
  }
  char proxPrint[200];
  if ((returnCode & 1) != 0) {
    // redo objective
    OsiSolverInterface *solver = model_->continuousSolver();
    const double *obj = solver->getObjCoefficients();
    solutionValue = -offset;
    int sumIncrease = 0.0;
    int sumDecrease = 0.0;
    int numberIncrease = 0;
    int numberDecrease = 0;
    for (int i = 0; i < numberColumns; i++) {
      solutionValue += obj[i] * betterSolution[i];
      if (isHeuristicInteger(solver, i)) {
        int change = static_cast< int >(floor(solutionIn[i] - betterSolution[i] + 0.5));
        if (change > 0) {
          numberIncrease++;
          sumIncrease += change;
        } else if (change < 0) {
          numberDecrease++;
          sumDecrease -= change;
        }
      }
    }
    sprintf(proxPrint, "Proximity search ran %d nodes (out of %d) - in new solution %d increased (%d), %d decreased (%d)",
      numberNodesDone_, numberNodes_,
      numberIncrease, sumIncrease, numberDecrease, sumDecrease);
    if (!numberIncrease && !numberDecrease) {
      // somehow tolerances are such that we can slip through
      // change for next time
      increment_ += CoinMax(increment_, fabs(solutionValue + offset) * 1.0e-10);
    }
  } else {
    sprintf(proxPrint, "Proximity search ran %d nodes - no new solution",
      numberNodesDone_);
  }
  model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
    << proxPrint
    << CoinMessageEol;

  delete newSolver;
  return returnCode;
}
// update model
void CbcHeuristicProximity::setModel(CbcModel *model)
{
  model_ = model;
  // Get a copy of original matrix
  assert(model_->solver());
  delete[] used_;
  int numberColumns = model->solver()->getNumCols();
  used_ = new int[numberColumns];
  memset(used_, 0, numberColumns * sizeof(int));
}

// Default Constructor
CbcHeuristicNaive::CbcHeuristicNaive()
  : CbcHeuristic()
{
  large_ = 1.0e6;
}

// Constructor with model - assumed before cuts

CbcHeuristicNaive::CbcHeuristicNaive(CbcModel &model)
  : CbcHeuristic(model)
{
  large_ = 1.0e6;
}

// Destructor
CbcHeuristicNaive::~CbcHeuristicNaive()
{
}

// Clone
CbcHeuristic *
CbcHeuristicNaive::clone() const
{
  return new CbcHeuristicNaive(*this);
}
// Create C++ lines to get to current state
void CbcHeuristicNaive::generateCpp(FILE *fp)
{
  CbcHeuristicNaive other;
  fprintf(fp, "0#include \"CbcHeuristicProximity.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicNaive naive(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "naive");
  if (large_ != other.large_)
    fprintf(fp, "3  naive.setLarge(%g);\n", large_);
  else
    fprintf(fp, "4  naive.setLarge(%g);\n", large_);
  fprintf(fp, "3  cbcModel->addHeuristic(&naive);\n");
}

// Copy constructor
CbcHeuristicNaive::CbcHeuristicNaive(const CbcHeuristicNaive &rhs)
  : CbcHeuristic(rhs)
  , large_(rhs.large_)
{
}

// Assignment operator
CbcHeuristicNaive &
CbcHeuristicNaive::operator=(const CbcHeuristicNaive &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    large_ = rhs.large_;
  }
  return *this;
}

// Resets stuff if model changes
void CbcHeuristicNaive::resetModel(CbcModel *model)
{
  CbcHeuristic::resetModel(model);
}
int CbcHeuristicNaive::solution(double &solutionValue,
  double *betterSolution)
{
  numCouldRun_++;
  // See if to do
  bool atRoot = model_->getNodeCount() == 0;
  int passNumber = model_->getCurrentPassNumber();
  if (!when() || (when() == 1 && model_->phase() != 1) || !atRoot || passNumber > 1)
    return 0; // switched off
  // Don't do if it was this heuristic which found solution!
  if (this == model_->lastHeuristic())
    return 0;
  numRuns_++;
  double cutoff;
  model_->solver()->getDblParam(OsiDualObjectiveLimit, cutoff);
  double direction = model_->solver()->getObjSense();
  cutoff *= direction;
  cutoff = CoinMin(cutoff, solutionValue);
  OsiSolverInterface *solver = model_->continuousSolver();
  if (!solver)
    solver = model_->solver();
  const double *colLower = solver->getColLower();
  const double *colUpper = solver->getColUpper();
  const double *objective = solver->getObjCoefficients();

  int numberColumns = model_->getNumCols();
  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();

  int i;
  bool solutionFound = false;
  CoinWarmStartBasis saveBasis;
  CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(solver->getWarmStart());
  if (basis) {
    saveBasis = *basis;
    delete basis;
  }
  // First just fix all integers as close to zero as possible
  OsiSolverInterface *newSolver = cloneBut(7); // wassolver->clone();
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (!isHeuristicInteger(newSolver, iColumn))
      continue;
    double lower = colLower[iColumn];
    double upper = colUpper[iColumn];
    double value;
    if (lower > 0.0)
      value = lower;
    else if (upper < 0.0)
      value = upper;
    else
      value = 0.0;
    newSolver->setColLower(iColumn, value);
    newSolver->setColUpper(iColumn, value);
  }
  newSolver->initialSolve();
  if (newSolver->isProvenOptimal()) {
    double solValue = newSolver->getObjValue() * direction;
    if (solValue < cutoff) {
      // we have a solution
      solutionFound = true;
      solutionValue = solValue;
      memcpy(betterSolution, newSolver->getColSolution(),
        numberColumns * sizeof(double));
      COIN_DETAIL_PRINT(printf("Naive fixing close to zero gave solution of %g\n", solutionValue));
      cutoff = solValue - model_->getCutoffIncrement();
    }
  }
  // Now fix all integers as close to zero if not zero or large cost
  int nFix = 0;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (!isHeuristicInteger(newSolver, iColumn))
      continue;
    double lower = colLower[iColumn];
    double upper = colUpper[iColumn];
    double value;
    if (fabs(objective[i]) > 0.0 && fabs(objective[i]) < large_) {
      nFix++;
      if (lower > 0.0)
        value = lower;
      else if (upper < 0.0)
        value = upper;
      else
        value = 0.0;
      newSolver->setColLower(iColumn, value);
      newSolver->setColUpper(iColumn, value);
    } else {
      // set back to original
      newSolver->setColLower(iColumn, lower);
      newSolver->setColUpper(iColumn, upper);
    }
  }
  const double *solution = solver->getColSolution();
  if (nFix) {
    newSolver->setWarmStart(&saveBasis);
    newSolver->setColSolution(solution);
    newSolver->initialSolve();
    if (newSolver->isProvenOptimal()) {
      double solValue = newSolver->getObjValue() * direction;
      if (solValue < cutoff) {
        // try branch and bound
        double *newSolution = new double[numberColumns];
        COIN_DETAIL_PRINT(printf("%d fixed after fixing costs\n", nFix));
        int returnCode = smallBranchAndBound(newSolver,
          numberNodes_, newSolution,
          solutionValue,
          solutionValue, "CbcHeuristicNaive1");
        if (returnCode < 0)
          returnCode = 0; // returned on size
        if ((returnCode & 2) != 0) {
          // could add cut
          returnCode &= ~2;
        }
        if (returnCode == 1) {
          // solution
          solutionFound = true;
          memcpy(betterSolution, newSolution,
            numberColumns * sizeof(double));
          COIN_DETAIL_PRINT(printf("Naive fixing zeros gave solution of %g\n", solutionValue));
          cutoff = solutionValue - model_->getCutoffIncrement();
        }
        delete[] newSolution;
      }
    }
  }
#if 1
  newSolver->setObjSense(-direction); // maximize
  newSolver->setWarmStart(&saveBasis);
  newSolver->setColSolution(solution);
  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    double value = solution[iColumn];
    double lower = colLower[iColumn];
    double upper = colUpper[iColumn];
    double newLower;
    double newUpper;
    if (isHeuristicInteger(newSolver, iColumn)) {
      newLower = CoinMax(lower, floor(value) - 2.0);
      newUpper = CoinMin(upper, ceil(value) + 2.0);
    } else {
      newLower = CoinMax(lower, value - 1.0e5);
      newUpper = CoinMin(upper, value + 1.0e-5);
    }
    newSolver->setColLower(iColumn, newLower);
    newSolver->setColUpper(iColumn, newUpper);
  }
  newSolver->initialSolve();
  if (newSolver->isProvenOptimal()) {
    double solValue = newSolver->getObjValue() * direction;
    if (solValue < cutoff) {
      nFix = 0;
      newSolver->setObjSense(direction); // correct direction
      //const double * thisSolution = newSolver->getColSolution();
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = solution[iColumn];
        double lower = colLower[iColumn];
        double upper = colUpper[iColumn];
        double newLower = lower;
        double newUpper = upper;
        if (isHeuristicInteger(newSolver, iColumn)) {
          if (value < lower + 1.0e-6) {
            nFix++;
            newUpper = lower;
          } else if (value > upper - 1.0e-6) {
            nFix++;
            newLower = upper;
          } else {
            newLower = CoinMax(lower, floor(value) - 2.0);
            newUpper = CoinMin(upper, ceil(value) + 2.0);
          }
        }
        newSolver->setColLower(iColumn, newLower);
        newSolver->setColUpper(iColumn, newUpper);
      }
      // try branch and bound
      double *newSolution = new double[numberColumns];
      COIN_DETAIL_PRINT(printf("%d fixed after maximizing\n", nFix));
      int returnCode = smallBranchAndBound(newSolver,
        numberNodes_, newSolution,
        solutionValue,
        solutionValue, "CbcHeuristicNaive1");
      if (returnCode < 0)
        returnCode = 0; // returned on size
      if ((returnCode & 2) != 0) {
        // could add cut
        returnCode &= ~2;
      }
      if (returnCode == 1) {
        // solution
        solutionFound = true;
        memcpy(betterSolution, newSolution,
          numberColumns * sizeof(double));
        COIN_DETAIL_PRINT(printf("Naive maximizing gave solution of %g\n", solutionValue));
        cutoff = solutionValue - model_->getCutoffIncrement();
      }
      delete[] newSolution;
    }
  }
#endif
  delete newSolver;
  return solutionFound ? 1 : 0;
}
// update model
void CbcHeuristicNaive::setModel(CbcModel *model)
{
  model_ = model;
}
// Default Constructor
CbcHeuristicCrossover::CbcHeuristicCrossover()
  : CbcHeuristic()
  , numberSolutions_(0)
  , useNumber_(3)
{
  setWhen(1);
}

// Constructor with model - assumed before cuts

CbcHeuristicCrossover::CbcHeuristicCrossover(CbcModel &model)
  : CbcHeuristic(model)
  , numberSolutions_(0)
  , useNumber_(3)
{
  setWhen(1);
  for (int i = 0; i < 10; i++)
    random_[i] = model.randomNumberGenerator()->randomDouble();
}

// Destructor
CbcHeuristicCrossover::~CbcHeuristicCrossover()
{
}

// Clone
CbcHeuristic *
CbcHeuristicCrossover::clone() const
{
  return new CbcHeuristicCrossover(*this);
}
// Create C++ lines to get to current state
void CbcHeuristicCrossover::generateCpp(FILE *fp)
{
  CbcHeuristicCrossover other;
  fprintf(fp, "0#include \"CbcHeuristicProximity.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicCrossover crossover(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "crossover");
  if (useNumber_ != other.useNumber_)
    fprintf(fp, "3  crossover.setNumberSolutions(%d);\n", useNumber_);
  else
    fprintf(fp, "4  crossover.setNumberSolutions(%d);\n", useNumber_);
  fprintf(fp, "3  cbcModel->addHeuristic(&crossover);\n");
}

// Copy constructor
CbcHeuristicCrossover::CbcHeuristicCrossover(const CbcHeuristicCrossover &rhs)
  : CbcHeuristic(rhs)
  , attempts_(rhs.attempts_)
  , numberSolutions_(rhs.numberSolutions_)
  , useNumber_(rhs.useNumber_)
{
  memcpy(random_, rhs.random_, 10 * sizeof(double));
}

// Assignment operator
CbcHeuristicCrossover &
CbcHeuristicCrossover::operator=(const CbcHeuristicCrossover &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    useNumber_ = rhs.useNumber_;
    attempts_ = rhs.attempts_;
    numberSolutions_ = rhs.numberSolutions_;
    memcpy(random_, rhs.random_, 10 * sizeof(double));
  }
  return *this;
}

// Resets stuff if model changes
void CbcHeuristicCrossover::resetModel(CbcModel *model)
{
  CbcHeuristic::resetModel(model);
}
int CbcHeuristicCrossover::solution(double &solutionValue,
  double *betterSolution)
{
  if (when_ == 0)
    return 0;
  numCouldRun_++;
  bool useBest = (numberSolutions_ != model_->getSolutionCount());
  if (!useBest && (when_ % 10) == 1)
    return 0;
  numberSolutions_ = model_->getSolutionCount();
  OsiSolverInterface *continuousSolver = model_->continuousSolver();
  int useNumber = CoinMin(model_->numberSavedSolutions(), useNumber_);
  if (useNumber < 2 || !continuousSolver)
    return 0;
  // Fix later
  if (!useBest)
    abort();
  numRuns_++;
  double cutoff;
  model_->solver()->getDblParam(OsiDualObjectiveLimit, cutoff);
  double direction = model_->solver()->getObjSense();
  cutoff *= direction;
  cutoff = CoinMin(cutoff, solutionValue);
  OsiSolverInterface *solver = cloneBut(2);
  // But reset bounds
  solver->setColLower(continuousSolver->getColLower());
  solver->setColUpper(continuousSolver->getColUpper());
  int numberColumns = solver->getNumCols();
  // Fixed
  double *fixed = new double[numberColumns];
  for (int i = 0; i < numberColumns; i++)
    fixed[i] = -COIN_DBL_MAX;
  int whichSolution[10];
  for (int i = 0; i < useNumber; i++)
    whichSolution[i] = i;
  for (int i = 0; i < useNumber; i++) {
    int k = whichSolution[i];
    const double *solution = model_->savedSolution(k);
    for (int j = 0; j < numberColumns; j++) {
      if (isHeuristicInteger(solver, j)) {
        if (fixed[j] == -COIN_DBL_MAX)
          fixed[j] = floor(solution[j] + 0.5);
        else if (fabs(fixed[j] - solution[j]) > 1.0e-7)
          fixed[j] = COIN_DBL_MAX;
      }
    }
  }
  const double *colLower = solver->getColLower();
  for (int i = 0; i < numberColumns; i++) {
    if (isHeuristicInteger(solver, i)) {
      double value = fixed[i];
      if (value != COIN_DBL_MAX) {
        if (when_ < 10) {
          solver->setColLower(i, value);
          solver->setColUpper(i, value);
        } else if (value == colLower[i]) {
          solver->setColUpper(i, value);
        }
      }
    }
  }
  int returnCode = smallBranchAndBound(solver, numberNodes_, betterSolution,
    solutionValue,
    solutionValue, "CbcHeuristicCrossover");
  if (returnCode < 0)
    returnCode = 0; // returned on size
  if ((returnCode & 2) != 0) {
    // could add cut
    returnCode &= ~2;
  }

  delete[] fixed;
  delete solver;
  return returnCode;
}
// update model
void CbcHeuristicCrossover::setModel(CbcModel *model)
{
  model_ = model;
  if (model) {
    for (int i = 0; i < 10; i++)
      random_[i] = model->randomNumberGenerator()->randomDouble();
  }
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
