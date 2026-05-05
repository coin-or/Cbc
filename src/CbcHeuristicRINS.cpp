// Copyright (C) 2006, International Business Machines
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
#include "CbcHeuristicRINS.hpp"
#include "CbcBranchActual.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"

// Default Constructor
CbcHeuristicRINS::CbcHeuristicRINS()
  : CbcHeuristic()
{
  numberSolutions_ = 0;
  numberSuccesses_ = 0;
  numberTries_ = 0;
  stateOfFixing_ = 0;
  shallowDepth_ = 0;
  lastNode_ = -999999;
  howOften_ = 100;
  decayFactor_ = 0.5;
  used_ = NULL;
  fixCloseMaxDist_ = 0.4;
  whereFrom_ = 1 + 8 + 16 + 255 * 256;
  whereFrom_ = 1 + 8 + 255 * 256;
}

// Constructor with model - assumed before cuts

CbcHeuristicRINS::CbcHeuristicRINS(CbcModel &model)
  : CbcHeuristic(model)
{
  numberSolutions_ = 0;
  numberSuccesses_ = 0;
  numberTries_ = 0;
  stateOfFixing_ = 0;
  shallowDepth_ = 0;
  lastNode_ = -999999;
  howOften_ = 100;
  decayFactor_ = 0.5;
  assert(model.solver());
  int numberColumns = model.solver()->getNumCols();
  used_ = new char[numberColumns];
  memset(used_, 0, numberColumns);
  fixCloseMaxDist_ = 0.4;
  whereFrom_ = 1 + 8 + 16 + 255 * 256;
  whereFrom_ = 1 + 8 + 255 * 256;
}

// Destructor
CbcHeuristicRINS::~CbcHeuristicRINS()
{
  delete[] used_;
}

// Clone
CbcHeuristic *
CbcHeuristicRINS::clone() const
{
  return new CbcHeuristicRINS(*this);
}

// Assignment operator
CbcHeuristicRINS &
CbcHeuristicRINS::operator=(const CbcHeuristicRINS &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    numberSolutions_ = rhs.numberSolutions_;
    howOften_ = rhs.howOften_;
    numberSuccesses_ = rhs.numberSuccesses_;
    numberTries_ = rhs.numberTries_;
    stateOfFixing_ = rhs.stateOfFixing_;
    lastNode_ = rhs.lastNode_;
    fixCloseMaxDist_ = rhs.fixCloseMaxDist_;
    delete[] used_;
    if (model_ && rhs.used_) {
      int numberColumns = model_->solver()->getNumCols();
      used_ = new char[numberColumns];
      memcpy(used_, rhs.used_, numberColumns);
    } else {
      used_ = NULL;
    }
  }
  return *this;
}

// Create C++ lines to get to current state
void CbcHeuristicRINS::generateCpp(FILE *fp)
{
  CbcHeuristicRINS other;
  fprintf(fp, "0#include \"CbcHeuristicRINS.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicRINS heuristicRINS(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "heuristicRINS");
  if (howOften_ != other.howOften_)
    fprintf(fp, "3  heuristicRINS.setHowOften(%d);\n", howOften_);
  else
    fprintf(fp, "4  heuristicRINS.setHowOften(%d);\n", howOften_);
  fprintf(fp, "3  cbcModel->addHeuristic(&heuristicRINS);\n");
}

// Copy constructor
CbcHeuristicRINS::CbcHeuristicRINS(const CbcHeuristicRINS &rhs)
  : CbcHeuristic(rhs)
  , numberSolutions_(rhs.numberSolutions_)
  , howOften_(rhs.howOften_)
  , numberSuccesses_(rhs.numberSuccesses_)
  , numberTries_(rhs.numberTries_)
  , stateOfFixing_(rhs.stateOfFixing_)
  , lastNode_(rhs.lastNode_)
  , fixCloseMaxDist_(rhs.fixCloseMaxDist_)
{
  if (model_ && rhs.used_) {
    int numberColumns = model_->solver()->getNumCols();
    used_ = new char[numberColumns];
    memcpy(used_, rhs.used_, numberColumns);
  } else {
    used_ = NULL;
  }
}
// Resets stuff if model changes
void CbcHeuristicRINS::resetModel(CbcModel * /*model*/)
{
  //CbcHeuristic::resetModel(model);
  delete[] used_;
  stateOfFixing_ = 0;
  if (model_ && used_) {
    int numberColumns = model_->solver()->getNumCols();
    used_ = new char[numberColumns];
    memset(used_, 0, numberColumns);
  } else {
    used_ = NULL;
  }
}
/*
  First tries setting a variable to better value.  If feasible then
  tries setting others.  If not feasible then tries swaps
  Returns 1 if solution, 0 if not */
int CbcHeuristicRINS::solution(double &solutionValue,
  double *betterSolution)
{
  numCouldRun_++;
  int returnCode = 0;
  const double *bestSolution = model_->bestSolution();
  if (!bestSolution)
    return 0; // No solution found yet
#ifdef HEURISTIC_INFORM
  printf("Entering heuristic %s - nRuns %d numCould %d when %d\n",
    heuristicName(), numRuns_, numCouldRun_, when_);
#endif
  if (numberSolutions_ < model_->getSolutionCount()) {
    // new solution - add info
    numberSolutions_ = model_->getSolutionCount();

    OsiSolverInterface *solver = model_->solver();
    int numberIntegers = model_->numberIntegers();
    const int *integerVariable = model_->integerVariable();

    int i;
    for (i = 0; i < numberIntegers; i++) {
      int iColumn = integerVariable[i];
      if (!isHeuristicInteger(solver, iColumn))
        continue;
      const OsiObject *object = model_->object(i);
      // get original bounds
      double originalLower;
      double originalUpper;
      getIntegerInformation(object, originalLower, originalUpper);
      double value = bestSolution[iColumn];
      if (value < originalLower) {
        value = originalLower;
      } else if (value > originalUpper) {
        value = originalUpper;
      }
      double nearest = floor(value + 0.5);
      // if away from lower bound mark that fact
      if (nearest > originalLower) {
        used_[iColumn] = 1;
      }
    }
  }
  int numberNodes = model_->getNodeCount();
  if (howOften_ == 100) {
    if (numberNodes < lastNode_ + 12)
      return 0;
    // Do at 50 and 100
    if ((numberNodes > 40 && numberNodes <= 50) || (numberNodes > 90 && numberNodes < 100))
      numberNodes = howOften_;
  }
  // Allow for infeasible nodes - so do anyway after a bit
  if (howOften_ >= 100 && numberNodes >= lastNode_ + 2 * howOften_) {
    numberNodes = howOften_;
  }
  if ((numberNodes % howOften_) == 0 && (model_->getCurrentPassNumber() <= 1 || model_->getCurrentPassNumber() == 999999)) {
    lastNode_ = model_->getNodeCount();
#ifdef RINS_CLOSE_DEBUG
    printf("RINS close-fix: entering work block at node=%d passNumber=%d\n",
      model_->getNodeCount(), model_->getCurrentPassNumber());
#endif
    OsiSolverInterface *solver = model_->solver();

    int numberIntegers = model_->numberIntegers();
    const int *integerVariable = model_->integerVariable();

    const double *currentSolution = solver->getColSolution();
    const int *used = model_->usedInSolution();
    OsiSolverInterface *newSolver = cloneBut(3); // was model_->continuousSolver()->clone();
    int numberColumns = newSolver->getNumCols();
    int numberContinuous = numberColumns - numberIntegers;

    double primalTolerance;
    solver->getDblParam(OsiPrimalTolerance, primalTolerance);

    int i;
    int nFix = 0;
    for (i = 0; i < numberIntegers; i++) {
      int iColumn = integerVariable[i];
      if (!isHeuristicInteger(solver, iColumn))
        continue;
      const OsiObject *object = model_->object(i);
      // get original bounds
      double originalLower;
      double originalUpper;
      getIntegerInformation(object, originalLower, originalUpper);
      double valueInt = bestSolution[iColumn];
      if (valueInt < originalLower) {
        valueInt = originalLower;
      } else if (valueInt > originalUpper) {
        valueInt = originalUpper;
      }
      if (fabs(currentSolution[iColumn] - valueInt) < 10.0 * primalTolerance) {
        double nearest = floor(valueInt + 0.5);
        /*
		  shallowDepth_
		  0 - normal
		  1 - only fix if at lb
		  2 - only fix if not at lb
		  3 - only fix if at lb and !used
		*/
        bool fix = false;
        switch (shallowDepth_) {
        case 0:
          fix = true;
          break;
        case 1:
          if (nearest == originalLower)
            fix = true;
          break;
        case 2:
          if (nearest != originalLower)
            fix = true;
          break;
        case 3:
          if (nearest == originalLower && !used[iColumn])
            fix = true;
          break;
        }
        if (fix) {
          newSolver->setColLower(iColumn, nearest);
          newSolver->setColUpper(iColumn, nearest);
          nFix++;
        }
      }
    }
    int divisor = 0;
    // Note: nFix counts only vars directly agreed between LP and best solution.
    // Variables that would be additionally fixed by consequence (via the fast
    // pre-processing inside smallBranchAndBound) are NOT counted here, so the
    // threshold below may be overly conservative.
#ifdef RINS_CLOSE_DEBUG
    {
      // Count candidates for close-fixing (regardless of threshold) for diagnostics
      int nExact = nFix;
      const double *cln = newSolver->getColLower();
      const double *cun = newSolver->getColUpper();
      int nCloseAvail = 0;
      for (int ii = 0; ii < numberIntegers; ii++) {
        int ic = integerVariable[ii];
        if (!isHeuristicInteger(solver, ic)) continue;
        if (cln[ic] == cun[ic]) continue;
        const OsiObject *obj = model_->object(ii);
        double lo, hi;
        getIntegerInformation(obj, lo, hi);
        double vb = bestSolution[ic];
        if (vb < lo) vb = lo; else if (vb > hi) vb = hi;
        double nr = floor(vb + 0.5);
        double dist = fabs(currentSolution[ic] - nr);
        if (dist >= 10.0 * primalTolerance && dist <= fixCloseMaxDist_) nCloseAvail++;
      }
      printf("RINS@node=%d: nExact=%d nIntegers=%d (%.1f%%) | closeAvail=%d | threshold=%s\n",
        model_->getNodeCount(), nExact, numberIntegers,
        100.0 * nExact / (numberIntegers ? numberIntegers : 1),
        nCloseAvail,
        (5 * nExact > numberIntegers) ? "MET" : "NOT MET");
    }
#endif

    // Close-fixing fallback: when the threshold is not met, additionally fix
    // integer variables whose current LP value is within fixCloseMaxDist_ of
    // the best-solution integer value.  Variables are sorted by closeness
    // (ascending) and greedily fixed until the threshold is satisfied.
    if (5 * nFix <= numberIntegers && fixCloseMaxDist_ > 0.0) {
      const double *colLowerNew = newSolver->getColLower();
      const double *colUpperNew = newSolver->getColUpper();
      double *sortDist = new double[numberIntegers];
      int *whichClose = new int[numberIntegers];
      int nClose = 0;
      for (i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        if (!isHeuristicInteger(solver, iColumn))
          continue;
        if (colLowerNew[iColumn] == colUpperNew[iColumn])
          continue; // already fixed in the exact-match pass
        const OsiObject *object = model_->object(i);
        double originalLower;
        double originalUpper;
        getIntegerInformation(object, originalLower, originalUpper);
        double valueInt = bestSolution[iColumn];
        if (valueInt < originalLower)
          valueInt = originalLower;
        else if (valueInt > originalUpper)
          valueInt = originalUpper;
        double nearest = floor(valueInt + 0.5);
        double dist = fabs(currentSolution[iColumn] - nearest);
        if (dist < 10.0 * primalTolerance || dist > fixCloseMaxDist_)
          continue;
        // Apply the same shallowDepth_ filter used in the exact-match pass
        bool include = false;
        switch (shallowDepth_) {
        case 0:
          include = true;
          break;
        case 1:
          include = (nearest == originalLower);
          break;
        case 2:
          include = (nearest != originalLower);
          break;
        case 3:
          include = (nearest == originalLower && !used[iColumn]);
          break;
        default:
          include = true;
          break;
        }
        if (include) {
          sortDist[nClose] = dist;
          whichClose[nClose++] = iColumn;
        }
      }
      // Sort closest first
      CoinSort_2(sortDist, sortDist + nClose, whichClose);
      int nFixBefore = nFix;
      // Greedily fix until threshold is met or candidates exhausted
      for (int j = 0; j < nClose && 5 * nFix <= numberIntegers; j++) {
        int iColumn = whichClose[j];
        double nearest = floor(bestSolution[iColumn] + 0.5);
        newSolver->setColLower(iColumn, nearest);
        newSolver->setColUpper(iColumn, nearest);
        nFix++;
      }
#ifdef RINS_CLOSE_DEBUG
      printf("RINS close-fix fallback: added %d vars (nFix %d->%d), candidates=%d, threshold_met=%s\n",
        nFix - nFixBefore, nFixBefore, nFix, nClose,
        (5 * nFix > numberIntegers) ? "YES" : "NO");
#endif
      delete[] sortDist;
      delete[] whichClose;
    }

    if (5 * nFix > numberIntegers) {
      if (numberContinuous > 2 * numberIntegers && ((nFix * 10 < numberColumns && !numRuns_ && numberTries_ > 2) || stateOfFixing_)) {
#define RINS_FIX_CONTINUOUS
#ifdef RINS_FIX_CONTINUOUS
        const double *colLower = newSolver->getColLower();
        //const double * colUpper = newSolver->getColUpper();
        int nAtLb = 0;
        //double sumDj=0.0;
        const double *dj = newSolver->getReducedCost();
        double direction = newSolver->getObjSense();
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (!isHeuristicInteger(newSolver, iColumn)) {
            double value = bestSolution[iColumn];
            if (value < colLower[iColumn] + 1.0e-8) {
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
          //double threshold = std::max((0.01*sumDj)/static_cast<double>(nAtLb),1.0e-6);
          int nFix2 = 0;
          for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (!isHeuristicInteger(newSolver, iColumn)) {
              double value = bestSolution[iColumn];
              if (value < colLower[iColumn] + 1.0e-8) {
                double djValue = dj[iColumn] * direction;
                if (djValue > 1.0e-6) {
                  sort[nFix2] = -djValue;
                  which[nFix2++] = iColumn;
                }
              }
            }
          }
          CoinSort_2(sort, sort + nFix2, which);
          divisor = 4;
          if (stateOfFixing_ > 0)
            divisor = stateOfFixing_;
          else if (stateOfFixing_ < -1)
            divisor = (-stateOfFixing_) - 1;
          nFix2 = std::min(nFix2, (numberColumns - nFix) / divisor);
          for (int i = 0; i < nFix2; i++) {
            int iColumn = which[i];
            newSolver->setColUpper(iColumn, colLower[iColumn]);
          }
          delete[] sort;
          delete[] which;
#ifdef CLP_INVESTIGATE2
          printf("%d integers have same value, and %d continuous fixed at lb\n",
            nFix, nFix2);
#endif
        }
#endif
      }
      if (solutionValue == -COIN_DBL_MAX) {
        // return fixings in betterSolution
        const double *colLower = newSolver->getColLower();
        const double *colUpper = newSolver->getColUpper();
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (colLower[iColumn] == colUpper[iColumn])
            betterSolution[iColumn] = colLower[iColumn];
          else
            betterSolution[iColumn] = COIN_DBL_MAX;
        }
        delete newSolver;
        return 0;
      }
      //printf("RINS %d integers have same value\n",nFix);
      if (model_->maximumSecondsReached()) {
        delete newSolver;
        return 0;
      }
      returnCode = smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue,
        model_->getCutoff(), "CbcHeuristicRINS");
      if (returnCode < 0) {
        returnCode = 0; // returned on size
        if (divisor) {
          stateOfFixing_ = -divisor; // say failed
        } else if (numberContinuous > 2 * numberIntegers && !numRuns_ && numberTries_ > 2) {
          stateOfFixing_ = -4; //start fixing
        }
      } else {
        numRuns_++;
        if (divisor)
          stateOfFixing_ = divisor; // say small enough
      }
      if ((returnCode & 1) != 0)
        numberSuccesses_++;
      //printf("return code %d",returnCode);
      if ((returnCode & 2) != 0) {
        // could add cut
        returnCode &= ~2;
        //printf("could add cut with %d elements (if all 0-1)\n",nFix);
      } else {
        //printf("\n");
      }
    }

    numberTries_++;
    if ((numberTries_ % 10) == 0 && numberSuccesses_ * 3 < numberTries_)
      howOften_ += static_cast< int >(howOften_ * decayFactor_);
    delete newSolver;
  }
  return returnCode;
}
// update model
void CbcHeuristicRINS::setModel(CbcModel *model)
{
  model_ = model;
  // Get a copy of original matrix
  assert(model_->solver());
  delete[] used_;
  int numberColumns = model->solver()->getNumCols();
  used_ = new char[numberColumns];
  memset(used_, 0, numberColumns);
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
