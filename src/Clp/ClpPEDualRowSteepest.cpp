// Corporation and others.  All Rights Reserved.
/*
   Authors

   Jeremy Omer

   Last update: april 10, 2015

 */

#include "CoinPragma.hpp"
#include "ClpSimplex.hpp"
#include "ClpPEDualRowSteepest.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpFactorization.hpp"
#include "CoinHelperFunctions.hpp"
#include "ClpMessage.hpp"
#include <cstdio>

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpPEDualRowSteepest::ClpPEDualRowSteepest(double psi, int mode)
  : ClpDualRowSteepest(mode)
  , modelPE_(NULL)
  , psi_(psi)
  , iCurrent_(0)
  , iInterval_(100)
  , updateCompatibles_(true)
  , coDegenCompatibles_(0)
  , coConsecutiveCompatibles_(0)
{
  iInterval_ = 100;
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpPEDualRowSteepest::ClpPEDualRowSteepest(const ClpPEDualRowSteepest &source)
  : ClpDualRowSteepest(source)
{
  modelPE_ = NULL;
  psi_ = source.psi_;
  iCurrent_ = source.iCurrent_;
  iInterval_ = source.iInterval_;
  updateCompatibles_ = source.updateCompatibles_;
  coDegenCompatibles_ = source.coDegenCompatibles_;
  coConsecutiveCompatibles_ = source.coConsecutiveCompatibles_;
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpPEDualRowSteepest::~ClpPEDualRowSteepest()
{
  delete modelPE_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpPEDualRowSteepest &
ClpPEDualRowSteepest::operator=(const ClpPEDualRowSteepest &rhs)
{
  if (this != &rhs) {
    ClpDualRowSteepest::operator=(rhs);
    delete modelPE_;
    modelPE_ = NULL;
  }
  return *this;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpDualRowPivot *ClpPEDualRowSteepest::clone(bool CopyData) const
{
  if (CopyData) {
    return new ClpPEDualRowSteepest(*this);
  } else {
    return new ClpPEDualRowSteepest(psi_);
  }
}
//extern int thinkDegenerate;
// Returns pivot row, -1 if none
int ClpPEDualRowSteepest::pivotRow()
{
  assert(model_);

  // Determine whether the set of compatible variables should be updated
  //
  // store the number of degenerate pivots on compatible variables and the
  // overal number of degenerate pivots
  //double progress = fabs(modelPE_->lastObjectiveValue() - model_->objectiveValue());
  //bool isLastDegenerate = progress <= 1.0e-12*fabs(model_->objectiveValue()) ? true:false;
  bool isLastDegenerate = fabs(model_->theta()) < 1.0e-7;
  bool isLastDegenerate2;
  // could fine tune a bit e.g. use 0.0 rather than tolerance
  if (model_->directionIn() > 0) {
    // coming in from lower bound
    isLastDegenerate2 = (model_->dualIn() < model_->dualTolerance());
  } else {
    // coming in from upper bound
    isLastDegenerate2 = (model_->dualIn() > -model_->dualTolerance());
  }
  isLastDegenerate = isLastDegenerate2;
  if (isLastDegenerate) {
    modelPE_->addDegeneratePivot();
    modelPE_->addDegeneratePivotConsecutive();

    if (modelPE_->isLastPivotCompatible()) {
      modelPE_->addDegenerateCompatiblePivot();
    }
  } else {
    modelPE_->resetDegeneratePivotsConsecutive();
  }

  // the compatible variables are updated when the number of degenerate pivots
  // on compatible variables is more than 20% the overall number of degenerate pivots
  if (modelPE_->isLastPivotCompatible()) {
    coConsecutiveCompatibles_++;
    if (isLastDegenerate) {
      coDegenCompatibles_++;
      if (coConsecutiveCompatibles_ >= 10 && 5 * coDegenCompatibles_ * model_->numberIterations() > modelPE_->coDegeneratePivots() * coConsecutiveCompatibles_) {
        updateCompatibles_ = true;
      }
    }
  }

  // For results comparison.
  // count the number of degenerate variables if psi==1
  // add the time spent in counting to the time limit
  //
  if (modelPE_->doStatistics()) {
    modelPE_->startTimer();
    if (psi_ >= 1 && iCurrent_ >= 100) {
      modelPE_->updateDualDegenerates();
      modelPE_->updateDualDegeneratesAvg(100);
      model_->setMaximumSeconds(36000.0 + modelPE_->timeCompatibility() - CoinCpuTime());
      iCurrent_ = 0;
    }
    modelPE_->stopTimer();
  }

  // Update the set of compatible variables if necessary
  //
  if (modelPE_->doStatistics())
    modelPE_->startTimer();
  double psiTmp = psi_;
  if ((psi_ < 1.0) && (iCurrent_ >= iInterval_) && (updateCompatibles_ || iCurrent_ >= 1000 /*|| !model_->factorization()->pivots()*/)) {
    // the compatible variables are never updated if the last pivot is non degenerate
    if (isLastDegenerate) {
      modelPE_->updateDualDegenerates();
      if (modelPE_->doStatistics()) {
        for (int i = 0; i < 4; i++)
          assert(!model_->rowArray(i)->getNumElements());
      }
      modelPE_->identifyCompatibleRows(model_->rowArray(2),
        model_->rowArray(1));
      if (modelPE_->doStatistics() > 3) {
        for (int i = 0; i < 4; i++)
          assert(!model_->rowArray(i)->getNumElements());
        char generalPrint[100];
        sprintf(generalPrint, "updating - iCurrent,iInterval %d,%d degenerate pivots %d ? %d codegen since last %d",
          iCurrent_, iInterval_,
          modelPE_->coDegeneratePivots(),
          modelPE_->coDualDegenerates(),
          coDegenCompatibles_);
        model_->messageHandler()->message(CLP_GENERAL,
          *model_->messagesPointer())
          << generalPrint << CoinMessageEol;
        // update the average number of degenerates and compatibles for statistics
        modelPE_->updateDualDegeneratesAvg(iCurrent_);
        modelPE_->updateCompatibleRowsAvg(iCurrent_);
      }
#if PE_DEBUG >= 1
      std::cout << "Number of compatible rows = " << modelPE_->coCompatibleRows() << " ; average  = " << modelPE_->coCompatibleRowsAvg() << std::endl;
      std::cout << "Number of degenerates = " << modelPE_->coDualDegenerates() << " ; average  = " << modelPE_->coDualDegeneratesAvg() << std::endl;
#endif

      // the checking frequency is modified to reflect what appears to be needed
      if (iCurrent_ == iInterval_)
        iInterval_ = std::max(50, iInterval_ - 50);
      else
        iInterval_ = std::min(300, iInterval_ + 50);

      // reset all the indicators that are used to decide whether the compatible
      // variables must be updated
      iCurrent_ = 0;
      updateCompatibles_ = false;
      coConsecutiveCompatibles_ = 0;
      coDegenCompatibles_ = 0;
    } else {
      iInterval_++;
    }
  }
  // otherwise, update the value of the priority factor depending on the number of
  // consecutive degenerate pivots
  //
  else {
    // the idea is that when a lot of consecutive pivots are degenerate, there is
    // a potentially high added value in putting a very high priroity on compatible
    // variables
    if (modelPE_->coDegeneratePivotsConsecutive() >= 10) {
      //psiTmp = 0.01;
      psiTmp = 0.25 * psi_;

#if PE_DEBUG >= 1
      std::cout << "Lower the priority factor to " << psiTmp << std::endl;
      std::cout << "Consecutive degenerate pivots=" << modelPE_->coDegeneratePivotsConsecutive() << std::endl;
#endif
    }
  }
  iCurrent_++;
  if (modelPE_->doStatistics())
    modelPE_->stopTimer();

  // Do the pricing and give priority to dual compatible variables
  //
  int i, iRow;
  double *infeas = infeasible_->denseVector();
  double largest = 0.0;
  int *index = infeasible_->getIndices();
  int number = infeasible_->getNumElements();
  const int *pivotVariable = model_->pivotVariable();
  int chosenRow = -1;
  int lastPivotRow = model_->pivotRow();
  assert(lastPivotRow < model_->numberRows());
  double tolerance = model_->currentPrimalTolerance();
  // we can't really trust infeasibilities if there is primal error
  // this coding has to mimic coding in checkPrimalSolution
  double error = std::min(1.0e-2, model_->largestPrimalError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance + error;
  // But cap
  tolerance = std::min(1000.0, tolerance);
  tolerance *= tolerance; // as we are using squares
  bool toleranceChanged = false;
  double *solution = model_->solutionRegion();
  double *lower = model_->lowerRegion();
  double *upper = model_->upperRegion();
  // do last pivot row one here
  //#define CLP_DUAL_FIXED_COLUMN_MULTIPLIER 10.0
  if (lastPivotRow >= 0 && lastPivotRow < model_->numberRows()) {
#ifdef CLP_DUAL_COLUMN_MULTIPLIER
    int numberColumns = model_->numberColumns();
#endif
    int iPivot = pivotVariable[lastPivotRow];
    double value = solution[iPivot];
    double lower = model_->lower(iPivot);
    double upper = model_->upper(iPivot);
    if (value > upper + tolerance) {
      value -= upper;
      value *= value;
#ifdef CLP_DUAL_COLUMN_MULTIPLIER
      if (iPivot < numberColumns)
        value *= CLP_DUAL_COLUMN_MULTIPLIER; // bias towards columns
#endif
      // store square in list
      if (infeas[lastPivotRow])
        infeas[lastPivotRow] = value; // already there
      else
        infeasible_->quickAdd(lastPivotRow, value);
    } else if (value < lower - tolerance) {
      value -= lower;
      value *= value;
#ifdef CLP_DUAL_COLUMN_MULTIPLIER
      if (iPivot < numberColumns)
        value *= CLP_DUAL_COLUMN_MULTIPLIER; // bias towards columns
#endif
      // store square in list
      if (infeas[lastPivotRow])
        infeas[lastPivotRow] = value; // already there
      else
        infeasible_->add(lastPivotRow, value);
    } else {
      // feasible - was it infeasible - if so set tiny
      if (infeas[lastPivotRow])
        infeas[lastPivotRow] = COIN_INDEXED_REALLY_TINY_ELEMENT;
    }
    number = infeasible_->getNumElements();
  }
  if (model_->numberIterations() < model_->lastBadIteration() + 200) {
    // we can't really trust infeasibilities if there is dual error
    if (model_->largestDualError() > model_->largestPrimalError()) {
      tolerance *= std::min(model_->largestDualError() / model_->largestPrimalError(), 1000.0);
      toleranceChanged = true;
    }
  }
  int numberWanted;
  if (mode_ < 2) {
    numberWanted = number + 1;
  } else if (mode_ == 2) {
    numberWanted = std::max(2000, number / 8);
  } else {
    int numberElements = model_->factorization()->numberElements();
    double ratio = static_cast< double >(numberElements) / static_cast< double >(model_->numberRows());
    numberWanted = std::max(2000, number / 8);
    if (ratio < 1.0) {
      numberWanted = std::max(2000, number / 20);
    } else if (ratio > 10.0) {
      ratio = number * (ratio / 80.0);
      if (ratio > number)
        numberWanted = number + 1;
      else
        numberWanted = std::max(2000, static_cast< int >(ratio));
    }
  }
  if (model_->largestPrimalError() > 1.0e-3)
    numberWanted = number + 1; // be safe
  int iPass;
  // Setup two passes
  int start[4];
  start[1] = number;
  start[2] = 0;
  double dstart = static_cast< double >(number) * model_->randomNumberGenerator()->randomDouble();
  start[0] = static_cast< int >(dstart);
  start[3] = start[0];

  int chosenRowComp = -1;
  double largestComp = 0.0;

  // only check the compatible variables when the bidimensional factor is less than 1
  // and the ratio of compatible variables is larger than 0.01
  // the percentage of compatible variables is computed as the ratio to the
  // smallest number among columns and rows
  bool checkCompatibles = true;
  double ratioCompatibles = static_cast< double >(modelPE_->coCompatibleRows()) / static_cast< double >(std::min(model_->numberRows(), model_->numberColumns()));

  if (psi_ >= 1.0 || ratioCompatibles < 0.01)
    checkCompatibles = false;

  // potentially, a partial pricing may be done
  for (iPass = 0; iPass < 2; iPass++) {
    int end = start[2 * iPass + 1];
    for (i = start[2 * iPass]; i < end; i++) {
      iRow = index[i];
      double value = infeas[iRow];
      if (value > tolerance) {
//#define OUT_EQ
#ifdef OUT_EQ
        {
          int iSequence = pivotVariable[iRow];
          if (upper[iSequence] == lower[iSequence])
            value *= 2.0;
        }
#endif
        double weight = std::min(weights_[iRow], 1.0e50);
        double largestMax = std::max(psiTmp * largest, largestComp);
        if (value > weight * largestMax) {
          // make last pivot row last resort choice
          if (iRow == lastPivotRow) {
            if (value * 1.0e-10 < largestMax * weight)
              continue;
            else
              value *= 1.0e-10;
          }
          int iSequence = pivotVariable[iRow];
          if (model_->flagged(iSequence)) {
            // just to make sure we don't exit before got something
            numberWanted++;
          } else if (checkCompatibles && modelPE_->isCompatibleRow(iRow) && value > weight * largestComp) {
            chosenRowComp = iRow;
            largestComp = value / weight;
          } else if (value > largest * weight) {
            chosenRow = iRow;
            largest = value / weight;
          }
        }
        numberWanted--;
        if (!numberWanted)
          break;
      }
    }
    if (!numberWanted)
      break;
  }
  // choose the a compatible or an incompatible row depending on the
  // largest values and on the factor of preference psiTmp
  if (modelPE_->doStatistics())
    modelPE_->startTimer();
  if (chosenRowComp >= 0 && largestComp >= psiTmp * largest) {

    if (modelPE_->doStatistics()) {
      // record the number of pivots done on compatible variables
      // that would not have been done without positive edge
      if (largestComp < largest)
        modelPE_->addPriorityPivot();
    }
    chosenRow = chosenRowComp;

#if PE_DEBUG >= 1
    modelPE_->checkCompatibilityRow(chosenRow);
#endif
  }
  if (chosenRow >= 0 && psiTmp < 1.0 && modelPE_->isCompatibleRow(chosenRow)) {
    modelPE_->isLastPivotCompatible(true);
    modelPE_->addCompatiblePivot();
  } else
    modelPE_->isLastPivotCompatible(false);

  if (modelPE_->doStatistics())
    modelPE_->stopTimer();

  if (chosenRow < 0 && toleranceChanged) {
    // won't line up with checkPrimalSolution - do again
    double saveError = model_->largestDualError();
    model_->setLargestDualError(0.0);
    // can't loop
    chosenRow = pivotRow();
    model_->setLargestDualError(saveError);
  }
  if (chosenRow < 0 && lastPivotRow < 0) {
    int nLeft = 0;
    for (int i = 0; i < number; i++) {
      int iRow = index[i];
      if (fabs(infeas[iRow]) > 1.0e-50) {
        index[nLeft++] = iRow;
      } else {
        infeas[iRow] = 0.0;
      }
    }
    infeasible_->setNumElements(nLeft);
    model_->setNumberPrimalInfeasibilities(nLeft);
  }

  modelPE_->updateLastObjectiveValue();
  return chosenRow;
}
/* Save weights - this may initialize weights as well
   This is as parent but may initialize ClpPESimplex
*/
void ClpPEDualRowSteepest::saveWeights(ClpSimplex *model, int mode)
{
  // See if we need to initialize ClpPESimplex
  if (!modelPE_ || model != modelPE_->clpModel() || !modelPE_->checkSize()) {
    delete modelPE_;
    modelPE_ = new ClpPESimplex(model);
  }
  ClpDualRowSteepest::saveWeights(model, mode);
}
/* Updates primal solution (and maybe list of candidates)
   Uses input vector which it deletes
   Computes change in objective function
   As ordinary steepest but checks for zero moves
*/
void ClpPEDualRowSteepest::updatePrimalSolution(CoinIndexedVector *input,
  double theta,
  double &changeInObjective)
{
  int iColumn = model_->sequenceIn();
  if (iColumn >= 0)
    modelPE_->updateCompatibleRows(iColumn);
  ClpDualRowSteepest::updatePrimalSolution(input, theta, changeInObjective);
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
