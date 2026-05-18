// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
/*
   Authors

   Jeremy Omer

   Last update: april 10, 2015

 */

#include "CoinPragma.hpp"
#include "ClpPEDualRowDantzig.hpp"
#include "CoinHelperFunctions.hpp"
#ifndef CLP_DUAL_COLUMN_MULTIPLIER
#define CLP_DUAL_COLUMN_MULTIPLIER 1.01
#endif

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpPEDualRowDantzig::ClpPEDualRowDantzig(double psi)
  : ClpDualRowDantzig()
  , modelPE_(NULL)
  , psi_(psi)
  , iCurrent_(0)
  , iInterval_(100)
  , updateCompatibles_(true)
  , coDegenCompatibles_(0)
  , coConsecutiveCompatibles_(0)
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpPEDualRowDantzig::ClpPEDualRowDantzig(const ClpPEDualRowDantzig &source)
  : ClpDualRowDantzig(source)
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
ClpPEDualRowDantzig::~ClpPEDualRowDantzig()
{
  delete modelPE_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpPEDualRowDantzig &
ClpPEDualRowDantzig::operator=(const ClpPEDualRowDantzig &rhs)
{
  if (this != &rhs) {
    ClpDualRowDantzig::operator=(rhs);
    delete modelPE_;
    modelPE_ = NULL;
  }
  return *this;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpDualRowPivot *ClpPEDualRowDantzig::clone(bool CopyData) const
{
  if (CopyData) {
    return new ClpPEDualRowDantzig(*this);
  } else {
    return new ClpPEDualRowDantzig(psi_);
  }
}

// Returns pivot row, -1 if none
int ClpPEDualRowDantzig::pivotRow()
{
  assert(model_);

  /* Determine whether the set of compatible variables should be updated
*/
  // store the number of degenerate pivots on compatible variables and the
  // overal number of degenerate pivots
  double progress = fabs(modelPE_->lastObjectiveValue() - model_->objectiveValue());
  bool isLastDegenerate = progress <= 1.0e-12 * fabs(model_->objectiveValue()) ? true : false;
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

  if (modelPE_->doStatistics()) {
    // For results comparison.
    // count the number of degenerate variables if psi==1
    // add the time spent in counting to the time limit

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
  if ((psi_ < 1.0) && (iCurrent_ >= iInterval_) && (updateCompatibles_ || iCurrent_ >= 1000)) {
    // the compatible variables are never updated if the last pivot is non degenerate
    // this could be counterproductive
    if (isLastDegenerate) {

      modelPE_->updateDualDegenerates();
      modelPE_->identifyCompatibleRows(model_->rowArray(2),
        model_->rowArray(1));

      if (modelPE_->doStatistics()) {
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
    } else
      iInterval_++;
  }
  // otherwise, update the value of the priority factor depending on the number of
  // consecutive degenerate pivots
  //
  else {
    // the idea is that when a lot of consecutive pivots are degenerate, there is
    // a potentially hign added value in putting a very high priroity on compatible
    // variables
    if (modelPE_->coDegeneratePivotsConsecutive() >= 10) {
      psiTmp = 0.01;

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
  int iRow;
  const int *pivotVariable = model_->pivotVariable();
  double tolerance = model_->currentPrimalTolerance();
  // we can't really trust infeasibilities if there is primal error
  if (model_->largestPrimalError() > 1.0e-8)
    tolerance *= model_->largestPrimalError() / 1.0e-8;
  double largest = 0.0;
  double largestComp = 0.0;
  int chosenRow = -1;
  int chosenRowComp = -1;
  int numberRows = model_->numberRows();
#ifdef CLP_DUAL_COLUMN_MULTIPLIER
  int numberColumns = model_->numberColumns();
#endif

  // only check the compatible variables when the bidimensional factor is less than 1
  // and the ratio of compatible variables is larger than 0.01
  // the percentage of compatible variables is computed as the ratio to the
  // smallest number among columns and rows
  bool checkCompatibles = true;
  double ratioCompatibles = static_cast< double >(modelPE_->coCompatibleRows()) / static_cast< double >(std::min(model_->numberRows(), model_->numberColumns()));

  if (psi_ >= 1.0 || ratioCompatibles < 0.01)
    checkCompatibles = false;

  // check the infeasibility of the variables (there is no partial pricing!)
  for (iRow = 0; iRow < numberRows; iRow++) {
    int iSequence = pivotVariable[iRow];
    double value = model_->solution(iSequence);
    double lower = model_->lower(iSequence);
    double upper = model_->upper(iSequence);
    double infeas = std::max(value - upper, lower - value);
    double largestMax = std::max(psi_ * largest, largestComp);
    if (infeas > tolerance) {
#ifdef CLP_DUAL_COLUMN_MULTIPLIER
      if (iSequence < numberColumns)
        infeas *= CLP_DUAL_COLUMN_MULTIPLIER;
#endif
      if (infeas > largestMax) {
        if (model_->flagged(iSequence)) {
        } else if (checkCompatibles && modelPE_->isCompatibleRow(iRow) && infeas > largestComp) {
          chosenRowComp = iRow;
          largestComp = infeas;
        } else if (infeas > largest) {
          chosenRow = iRow;
          largest = infeas;
        }
      }
    }
  }

  // choose a compatible or an incompatible row depending on the
  // largest values and on the factor of preference psi_
  if (modelPE_->doStatistics())
    modelPE_->startTimer();
  if (chosenRowComp >= 0 && largestComp >= psiTmp * largest) {
    chosenRow = chosenRowComp;

    if (modelPE_->doStatistics()) {
      // record the number of pivots done on compatible variables
      // that would not have been done without positive edge
      if (largestComp < largest)
        modelPE_->addPriorityPivot();
    }
#if PE_DEBUG >= 1
    modelPE_->checkCompatibilityRow(chosenRow);
#endif
  }
  if (psi_ < 1 && modelPE_->isCompatibleRow(chosenRow)) {
    modelPE_->isLastPivotCompatible(true);
    modelPE_->addCompatiblePivot();
  } else
    modelPE_->isLastPivotCompatible(false);
  if (modelPE_->doStatistics())
    modelPE_->stopTimer();

  modelPE_->updateLastObjectiveValue();
  return chosenRow;
}

//-------------------------------------------------------------------
// updateWeights
// Update the compatible variables and
// call the base class method to update weights
//-------------------------------------------------------------------

double ClpPEDualRowDantzig::updateWeights(CoinIndexedVector *input,
  CoinIndexedVector *spare,
  CoinIndexedVector *spare2,
  CoinIndexedVector *updatedColumn)
{
  double value = ClpDualRowDantzig::updateWeights(input, spare, spare2, updatedColumn);

  return value;
}
/* Save weights - this may initialize weights as well
   This is as parent but may initialize ClpPESimplex
*/
void ClpPEDualRowDantzig::saveWeights(ClpSimplex *model, int mode)
{
  // See if we need to initialize ClpPESimplex
  if (!modelPE_ || model != modelPE_->clpModel()) {
    delete modelPE_;
    modelPE_ = new ClpPESimplex(model);
  }
  ClpDualRowDantzig::saveWeights(model, mode);
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
