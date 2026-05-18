// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

/*
   Authors

   Jeremy Omer, Mehdi Towhidi

   Last update: april 10, 2015

 */

#include "CoinPragma.hpp"

#include <cstdio>

#include "ClpPEPrimalColumnDantzig.hpp"

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpPEPrimalColumnDantzig::ClpPEPrimalColumnDantzig(double psi)
  : ClpPrimalColumnDantzig()
  , modelPE_(NULL)
  , psi_(psi)
  , iCurrent_(0)
  , iInterval_(100)
  , coDegenCompatibles_(0)
  , coConsecutiveCompatibles_(0)
  , updateCompatibles_(true)
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpPEPrimalColumnDantzig::ClpPEPrimalColumnDantzig(const ClpPEPrimalColumnDantzig &source)
  : ClpPrimalColumnDantzig(source)
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
ClpPEPrimalColumnDantzig::~ClpPEPrimalColumnDantzig()
{
  delete modelPE_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpPEPrimalColumnDantzig &
ClpPEPrimalColumnDantzig::operator=(const ClpPEPrimalColumnDantzig &rhs)
{
  if (this != &rhs) {
    ClpPrimalColumnDantzig::operator=(rhs);
    delete modelPE_;
    modelPE_ = NULL;
  }
  return *this;
}

//-------------------------------------------------------------------
//// Clone
////-------------------------------------------------------------------
ClpPrimalColumnPivot *ClpPEPrimalColumnDantzig::clone(bool CopyData) const
{
  if (CopyData) {
    return new ClpPEPrimalColumnDantzig(*this);
  } else {
    return new ClpPEPrimalColumnDantzig(psi_);
  }
}

//-------------------------------------------------------------------
// pivotColumn
// Returns pivot column, -1 if none
// With this class, a two-dimensional pricing is applied to choose
// between the most negative reduced cost and the most negative
// recuded cost of the compatible variables.
////-------------------------------------------------------------------
int ClpPEPrimalColumnDantzig::pivotColumn(CoinIndexedVector *updates,
  CoinIndexedVector *spareRow1,
  CoinIndexedVector *spareRow2,
  CoinIndexedVector *spareColumn1,
  CoinIndexedVector *spareColumn2)
{
  assert(model_);

#if PE_DEBUG >= 1
  std::cout << "@@@@@@@@@@PEPrimalColumnDantzig::pivotColumn\n";
#endif

  int iSection, j;
  int number;
  int *index;
  double *updateBy;
  double *reducedCost;

  bool anyUpdates;

  // If the input vector updates is not empty, the dual solution and
  // the reduced cost needs to be updated
  //
  if (updates->getNumElements()) {
    anyUpdates = true;
  } else {
    // sub flip - nothing to do
    anyUpdates = false;
  }
  if (anyUpdates) {
// compute y:  y^T * A = updates^T and updates <- spareRow2 + updates
#ifdef PE_STATISTICS
    model_->factorization()->doStatistics(false);
#endif
    model_->factorization()->updateColumnTranspose(spareRow2, updates);
#ifdef PE_STATISTICS
    model_->factorization()->doStatistics(true);
#endif
    // put row of tableau in rowArray and columnArray
    model_->clpMatrix()->transposeTimes(model_, -1.0,
      updates, spareColumn2, spareColumn1);
    for (iSection = 0; iSection < 2; iSection++) {

      reducedCost = model_->djRegion(iSection);

      if (!iSection) {
        number = updates->getNumElements();
        index = updates->getIndices();
        updateBy = updates->denseVector();
      } else {
        number = spareColumn1->getNumElements();
        index = spareColumn1->getIndices();
        updateBy = spareColumn1->denseVector();
      }

      for (j = 0; j < number; j++) {
        int iSequence = index[j];
        double value = reducedCost[iSequence];
        value -= updateBy[j];
        updateBy[j] = 0.0;
        reducedCost[iSequence] = value;
      }
    }
    updates->setNumElements(0);
    spareColumn1->setNumElements(0);
  }

  // Determine whether the set of compatible variables should be updated
  //
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
    //
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
      modelPE_->updatePrimalDegenerates();
      modelPE_->identifyCompatibleCols(model_->numberRows() + model_->numberColumns(), NULL, spareRow2, spareRow1);

      if (modelPE_->doStatistics()) {
        // update the average number of degenerates and compatibles for statistics
        modelPE_->updatePrimalDegeneratesAvg(iCurrent_);
        modelPE_->updateCompatibleColsAvg(iCurrent_);
      }

#ifdef PE_TEST
      std::cout << "coDegen=" << coDegenCompatibles_ << " ; coComp = " << coConsecutiveCompatibles_
                << " ; iCurrent_ = " << iCurrent_ << " ;compatibleColumns = "
                << modelPE_->coCompatibleCols() << std::endl;
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

  // Updates of dual solutions and compatible variables finished,
  // now do the pricing
  //

  // we can't really trust infeasibilities if there is primal error
  //double largest = model_->currentPrimalTolerance();
  //if (model_->largestDualError() > 1.0e-8)
  //  largest *= model_->largestDualError() / 1.0e-8;

  // initialize the best reduced cost values
  double dualTolerance = model_->dualTolerance();
  double bestDj = 1.0e-30;
  int bestSequence = -1;
  double bestDjComp = 1.0e-30;
  int bestSequenceComp = -1;

  number = model_->numberRows() + model_->numberColumns();
  reducedCost = model_->djRegion();

  // only check the compatible variables when the bidimensional factor is less than 1
  // and the ratio of compatible variables is larger than 0.01
  bool checkCompatibles = true;
  double ratioCompatibles = static_cast< double >(modelPE_->coCompatibleCols()) / static_cast< double >((model_->numberRows() + model_->numberColumns()));

  if (psi_ >= 1.0 || ratioCompatibles < 0.01)
    checkCompatibles = false;

  int indicesToConsiderCount = 0;
  for (int iSequence = 0; iSequence < number; iSequence++) {
    // check flagged variable
    if (!model_->flagged(iSequence)) {
      double value = reducedCost[iSequence];
      double largestDj = std::max(psi_ * bestDj, bestDjComp);
      ClpSimplex::Status status = model_->getStatus(iSequence);

      // we choose the nonbasic column whose reduced cost is either
      // the most negative or the most positive one, depending on the
      // status of the variable
      // if a variable is compatible, a factor  0 < psi_ < 1 gives
      // a greater priority to it
      switch (status) {

        // basic and fixed variable cannot be chosen to enter the basis
      case ClpSimplex::basic:
      case ClpSimplex::isFixed:
        break;
        // free and superbasic are given priority using the 0.1 factor
        // since these variables never leave basis once they're in
      case ClpSimplex::isFree:
      case ClpSimplex::superBasic:
        value = fabs(value);
        if (checkCompatibles && modelPE_->isCompatibleCol(iSequence) && value > 0.1 * bestDjComp) {
          bestDjComp = 10.0 * value;
          bestSequenceComp = iSequence;
        } else if (value > 0.1 * bestDj) {
          bestDj = 10.0 * value;
          bestSequence = iSequence;
        }
        break;
        // just be careful with the sign of the reduced cost for the other
        // variables and give priority to the compatible ones
      case ClpSimplex::atUpperBound:
        if (value > largestDj) {
          if (checkCompatibles && modelPE_->isCompatibleCol(iSequence)) {
            bestDjComp = value;
            bestSequenceComp = iSequence;
          } else if (value > bestDj) {
            bestDj = value;
            bestSequence = iSequence;
          }
        }
        break;
      case ClpSimplex::atLowerBound:
        if (value < -largestDj) {
          if (checkCompatibles && modelPE_->isCompatibleCol(iSequence)) {
            bestDjComp = -value;
            bestSequenceComp = iSequence;
          } else if (value < -bestDj) {
            bestDj = -value;
            bestSequence = iSequence;
          }
        }
        break;
      }
    }
  }
  // choose a compatible or an incompatible row depending on the
  // largest values and on the factor of preference psi_
  if (modelPE_->doStatistics())
    modelPE_->startTimer();
  if (bestSequenceComp >= 0 && bestDjComp >= psiTmp * bestDj) {
    bestSequence = bestSequenceComp;

    // record the number of pivots done on compatible variables
    // that would not have been done without positive edge
    if (modelPE_->doStatistics())
      if (bestDjComp < bestDj)
        modelPE_->addPriorityPivot();
  }
  if (psi_ < 1 && modelPE_->isCompatibleCol(bestSequence)) {
    modelPE_->isLastPivotCompatible(true);
    modelPE_->addCompatiblePivot();
  } else
    modelPE_->isLastPivotCompatible(false);
  if (modelPE_->doStatistics())
    modelPE_->stopTimer();

  // save the current objective value
  modelPE_->updateLastObjectiveValue();

  return bestSequence;
}
/* Save weights - this may initialize weights as well
   This is as parent but may initialize ClpPESimplex
*/
void ClpPEPrimalColumnDantzig::saveWeights(ClpSimplex *model, int mode)
{
  // See if we need to initialize ClpPESimplex
  if (!modelPE_ || model != modelPE_->clpModel()) {
    delete modelPE_;
    modelPE_ = new ClpPESimplex(model);
  }
  ClpPrimalColumnDantzig::saveWeights(model, mode);
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
