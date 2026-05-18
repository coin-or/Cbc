// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
/*
   Authors

   Jeremy Omer, Mehdi Towhidi

   Last update: april 10, 2015

 */

#include "CoinPragma.hpp"

#include "ClpPEPrimalColumnSteepest.hpp"

#include "ClpSimplex.hpp"
#include "ClpMessage.hpp"
#include "ClpNonLinearCost.hpp"

#include <stdio.h>
//#define CLP_DEBUG
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpPEPrimalColumnSteepest::ClpPEPrimalColumnSteepest(double psi, int mode)
  : ClpPrimalColumnSteepest(mode)
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
ClpPEPrimalColumnSteepest::ClpPEPrimalColumnSteepest(const ClpPEPrimalColumnSteepest &source)
  : ClpPrimalColumnSteepest(source)
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
ClpPEPrimalColumnSteepest::~ClpPEPrimalColumnSteepest()
{
  delete modelPE_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpPEPrimalColumnSteepest &
ClpPEPrimalColumnSteepest::operator=(const ClpPEPrimalColumnSteepest &rhs)
{
  if (this != &rhs) {
    ClpPrimalColumnSteepest::operator=(rhs);
    delete modelPE_;
    modelPE_ = NULL;
  }
  return *this;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpPrimalColumnPivot *ClpPEPrimalColumnSteepest::clone(bool CopyData) const
{
  if (CopyData) {
    return new ClpPEPrimalColumnSteepest(*this);
  } else {
    return new ClpPEPrimalColumnSteepest(psi_);
  }
}

// These have to match ClpPackedMatrix version
#define TRY_NORM 1.0e-4
#define ADD_ONE 1.0

// Returns pivot column, -1 if none
/*   The Packed CoinIndexedVector updates has cost updates - for normal LP
that is just +-weight where a feasibility changed.  It also has
reduced cost from last iteration in pivot row*/
int ClpPEPrimalColumnSteepest::pivotColumn(CoinIndexedVector *updates,
  CoinIndexedVector *spareRow1,
  CoinIndexedVector *spareRow2,
  CoinIndexedVector *spareColumn1,
  CoinIndexedVector *spareColumn2)
{
  // check that the model exists and has the proper form
  assert(model_);
  assert(!model_->nonLinearCost()->lookBothWays() && model_->algorithm() != 2);

#if PE_DEBUG >= 1
  std::cout << "@@@@@@@@@@ClpPEPrimalColumnSteepest::pivotColumn\n";
#endif

  int number = 0;
  int *index;

  double tolerance = model_->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = std::min(1.0e-2, model_->largestDualError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance + error;
  int pivotRow = model_->pivotRow();
  int anyUpdates;
  double *infeas = infeasible_->denseVector();

  // Local copy of mode so can decide what to do
  int switchType;
  if (mode_ == 4)
    switchType = 5 - numberSwitched_;
  else if (mode_ >= 10)
    switchType = 3;
  else
    switchType = mode_;

  /* switchType -
	0 - all exact devex
	1 - all steepest
	2 - some exact devex
	3 - auto some exact devex
	4 - devex
	Default value is switchType = 3
	I removed the mode >= 5 out of simplicity
	*/

  // Definition in ClpGubMatrix, mode is set to 4
  model_->clpMatrix()->dualExpanded(model_, updates, NULL, 4);

  if (updates->getNumElements() > 1) {
    // would have to have two goes for devex, three for steepest
    anyUpdates = 2;
  } else if (updates->getNumElements()) {
    if (updates->getIndices()[0] == pivotRow && fabs(updates->denseVector()[0]) > 1.0e-6) {
      // reasonable size
      anyUpdates = 1;
    } else {
      // too small
      anyUpdates = 2;
    }
  } else if (pivotSequence_ >= 0) {
    // just after re-factorization
    anyUpdates = -1;
  } else {
    // sub flip - nothing to do
    anyUpdates = 0;
  }
  int sequenceOut = model_->sequenceOut();

  if (anyUpdates == 1) {
    if (switchType < 4) {
      // exact etc when can use dj
      djsAndSteepest(updates, spareRow2,
        spareColumn1, spareColumn2);
    } else {
      // devex etc when can use dj
      djsAndDevex(updates, spareRow2,
        spareColumn1, spareColumn2);
    }
  } else if (anyUpdates == -1) {
    if (switchType < 4) {
      // exact etc when djs okay
      justSteepest(updates, spareRow2,
        spareColumn1, spareColumn2);
    } else {
      // devex etc when djs okay
      justDevex(updates, spareRow2,
        spareColumn1, spareColumn2);
    }
  } else if (anyUpdates == 2) {
    if (switchType < 4) {
      // exact etc when have to use pivot
      djsAndSteepest2(updates, spareRow2,
        spareColumn1, spareColumn2);
    } else {
      // devex etc when have to use pivot
      djsAndDevex2(updates, spareRow2,
        spareColumn1, spareColumn2);
    }
  }

  // make sure outgoing from last iteration okay
  if (sequenceOut >= 0) {
    ClpSimplex::Status status = model_->getStatus(sequenceOut);
    double value = model_->reducedCost(sequenceOut);

    switch (status) {

    case ClpSimplex::basic:
    case ClpSimplex::isFixed:
      break;
    case ClpSimplex::isFree:
    case ClpSimplex::superBasic:
      if (fabs(value) > FREE_ACCEPT * tolerance) {
        // we are going to bias towards free (but only if reasonable)
        value *= FREE_BIAS;
        // store square in list
        if (infeas[sequenceOut])
          infeas[sequenceOut] = value * value; // already there
        else
          infeasible_->quickAdd(sequenceOut, value * value);
      } else {
        infeasible_->zero(sequenceOut);
      }
      break;
    case ClpSimplex::atUpperBound:
      if (value > tolerance) {
        // store square in list
        if (infeas[sequenceOut])
          infeas[sequenceOut] = value * value; // already there
        else
          infeasible_->quickAdd(sequenceOut, value * value);
      } else {
        infeasible_->zero(sequenceOut);
      }
      break;
    case ClpSimplex::atLowerBound:
      if (value < -tolerance) {
        // store square in list
        if (infeas[sequenceOut])
          infeas[sequenceOut] = value * value; // already there
        else
          infeasible_->quickAdd(sequenceOut, value * value);
      } else {
        infeasible_->zero(sequenceOut);
      }
    }
  }

  /* Determine whether the set of compatible variables should be updated
	*/
  // store the number of degenerate pivots on compatible variables and the
  // overall number of degenerate pivots
  //double progress = fabs(modelPE_->lastObjectiveValue() - model_->objectiveValue());
  //bool isLastDegenerate = progress <= 1.0e-12*fabs(model_->objectiveValue()) ? true:false;
  bool isLastDegenerate = fabs(model_->theta()) < 1.01e-7;
  bool isLastDegenerate2;
  // could fine tune a bit e.g. use 0.0 rather than tolerance
  if (model_->directionOut() < 0) {
    // going out to lower bound
    isLastDegenerate2 = (model_->valueOut() - model_->lowerOut() < model_->primalTolerance());
  } else {
    // going out to upper bound
    isLastDegenerate2 = (model_->upperOut() - model_->valueOut() < model_->primalTolerance());
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

  if (modelPE_->doStatistics()) {
    /* For results comparison.
	count the number of degenerate variables if psi==1
	add the time spent in counting to the time limit
	*/
    modelPE_->startTimer();
    if (psi_ >= 1 && iCurrent_ >= 100) {
      modelPE_->updateDualDegenerates();
      modelPE_->updateDualDegeneratesAvg(100);
      //model_->setMaximumSeconds(36000.0+modelPE_->timeCompatibility()-CoinCpuTime());
      iCurrent_ = 0;
    }
    modelPE_->stopTimer();
  }

  /* Update the set of compatible variables if necessary
	*/
  if (modelPE_->doStatistics())
    modelPE_->startTimer();
  double psiTmp = psi_;
  if ((psi_ < 1.0) && (iCurrent_ >= iInterval_) && (updateCompatibles_ || iCurrent_ >= 1000)) {
    // the compatible variables are never updated if the last pivot is non degenerate
    // this could be counterproductive
    if (isLastDegenerate) {
      modelPE_->updatePrimalDegenerates();
      for (int i = 0; i < 4; i++)
        assert(!model_->rowArray(i)->getNumElements());
      modelPE_->identifyCompatibleCols(model_->numberRows() + model_->numberColumns(), NULL, spareRow2, spareRow1);
      for (int i = 0; i < 4; i++)
        assert(!model_->rowArray(i)->getNumElements());

      if (modelPE_->doStatistics()) {
        // update the average number of degenerates and compatibles for statistics
        modelPE_->updatePrimalDegeneratesAvg(iCurrent_);
        modelPE_->updateCompatibleColsAvg(iCurrent_);
        if (modelPE_->doStatistics() > 3) {
          char generalPrint[100];

          sprintf(generalPrint, "coDegen = %d; coComp = %d; iCurrent_ = %d; compatibleColumns = %d",
            coDegenCompatibles_, coConsecutiveCompatibles_,
            iCurrent_, modelPE_->coCompatibleCols());
          model_->messageHandler()->message(CLP_GENERAL,
            *model_->messagesPointer())
            << generalPrint << CoinMessageEol;
        }
      }
      // switch off ? - maybe just until drops back
      if (modelPE_->coCompatibleCols() * 1.0 > model_->numberColumns()) {
        printf("switching off pe\n");
        psi_ = 1.0;
      }

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
  /* otherwise, update the value of the priority factor depending on the number of
	 * consecutive degenerate pivots
	 */
  else {
    // the idea is that when a lot of consecutive pivots are degenerate, there is
    // a potentially high added value in putting a very high priority on compatible
    // variables
    if (modelPE_->coDegeneratePivotsConsecutive() >= 10) {
      psiTmp = 0.01;
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

  /* Update of the dual solution and compatible variables is finished,
	** now do the pricing */

  // See what sort of pricing
  int numberWanted = 10;
  number = infeasible_->getNumElements();
  int numberColumns = model_->numberColumns();
  int numberRows = model_->numberRows();
  // ratio is done on number of rows here
  double ratio = static_cast< double >(sizeFactorization_) / static_cast< double >(numberRows);
  if (switchType == 4) {
    // Still in devex mode
    // Go to steepest if lot of iterations?
    if (ratio < 5.0) {
      numberWanted = std::max(2000, number / 10);
      numberWanted = std::max(numberWanted, numberColumns / 20);
    } else if (ratio < 7.0) {
      numberWanted = std::max(2000, number / 5);
      numberWanted = std::max(numberWanted, numberColumns / 10);
    } else {
      // we can zero out
      updates->clear();
      spareColumn1->clear();
      switchType = 3;
      // initialize
      pivotSequence_ = -1;
      pivotRow = -1;
      numberSwitched_++;
      // Make sure will re-do
      delete[] weights_;
      weights_ = NULL;
      saveWeights(model_, 4);
      COIN_DETAIL_PRINT(printf("switching to exact %d nel ratio %g\n", sizeFactorization_, ratio));
      updates->clear();
    }
    if (model_->numberIterations() % 1000 == 0)
      COIN_DETAIL_PRINT(printf("numels %d ratio %g wanted %d type x\n", sizeFactorization_, ratio, numberWanted));
  }
  if (switchType < 4) {
    if (switchType < 2) {
      numberWanted = number + 1;
    } else if (switchType == 2) {
      numberWanted = std::max(2000, number / 8);
    } else {
      if (ratio < 1.0) {
        numberWanted = std::max(2000, number / 20);
      } else if (ratio < 5.0) {
        numberWanted = std::max(2000, number / 10);
        numberWanted = std::max(numberWanted, numberColumns / 40);
      } else if (ratio < 10.0) {
        numberWanted = std::max(2000, number / 8);
        numberWanted = std::max(numberWanted, numberColumns / 20);
      } else {
        ratio = number * (ratio / 80.0);
        if (ratio > number) {
          numberWanted = number + 1;
        } else {
          numberWanted = std::max(2000, static_cast< int >(ratio));
          numberWanted = std::max(numberWanted, numberColumns / 10);
        }
      }
    }
  }

  // initialize the best reduced cost values
  double bestDj = 1.0e-30;
  int bestSequence = -1;
  double bestDjComp = 1.0e-30;
  int bestSequenceComp = -1;

  int i, iSequence;
  index = infeasible_->getIndices();
  number = infeasible_->getNumElements();
  // Re-sort infeasible every 100 pivots
  // Not a good idea
  if (0 && model_->factorization()->pivots() > 0 && (model_->factorization()->pivots() % 100) == 0) {
    int nLook = model_->numberRows() + numberColumns;
    number = 0;
    for (i = 0; i < nLook; i++) {
      if (infeas[i]) {
        if (fabs(infeas[i]) > COIN_INDEXED_TINY_ELEMENT)
          index[number++] = i;
        else
          infeas[i] = 0.0;
      }
    }
    infeasible_->setNumElements(number);
  }
  if (model_->numberIterations() < model_->lastBadIteration() + 200 && model_->factorization()->pivots() > 10) {
    // we can't really trust infeasibilities if there is dual error
    double checkTolerance = 1.0e-8;
    if (model_->largestDualError() > checkTolerance)
      tolerance *= model_->largestDualError() / checkTolerance;
    // But cap
    tolerance = std::min(1000.0, tolerance);
  }

  // stop last one coming immediately
  double saveOutInfeasibility = 0.0;
  if (sequenceOut >= 0) {
    saveOutInfeasibility = infeas[sequenceOut];
    infeas[sequenceOut] = 0.0;
  }
  if (model_->factorization()->pivots() && model_->numberPrimalInfeasibilities())
    tolerance = std::max(tolerance, 1.0e-10 * model_->infeasibilityCost());
  tolerance *= tolerance; // as we are using squares

  // only check the compatible variables when the bidimensional factor is less than 1
  // and the ratio of compatible variables is larger than 0.01
  bool checkCompatibles = true;
  double ratioCompatibles = static_cast< double >(modelPE_->coCompatibleCols()) / static_cast< double >((model_->numberRows() + model_->numberColumns()));
  double ratioCompatibles2 = static_cast< double >(modelPE_->coCompatibleCols()) / static_cast< double >(model_->numberColumns());

  if (psi_ >= 1.0 || ratioCompatibles < 0.01)
    checkCompatibles = false;
  if (ratioCompatibles2 > 0.5)
    checkCompatibles = false;
  // proceed to a partial pricing in two phases checking the weighted reduced
  // costs of the variables until numberWanted variables have been checked,
  // the number of variables explored in the first phase is randomly generated
  int iPass;
  int start[4];
  start[1] = number;
  start[2] = 0;
  double dstart = static_cast< double >(number) * model_->randomNumberGenerator()->randomDouble();
  start[0] = static_cast< int >(dstart);
  start[3] = start[0];
  for (iPass = 0; iPass < 2; iPass++) {
    int end = start[2 * iPass + 1];

    for (i = start[2 * iPass]; i < end; i++) {
      iSequence = index[i];
      double value = infeas[iSequence];
      double weight = weights_[iSequence];
      double weightedDj = weight * bestDj;
      double largestWeightedDj = std::max(psi_ * weightedDj, weight * bestDjComp);
      if (value > tolerance) {
        if (value > largestWeightedDj) {
          if (model_->flagged(iSequence)) {
            // just to make sure we don't exit before got something
            numberWanted++;
          } else if (checkCompatibles && modelPE_->isCompatibleCol(iSequence)) {
            // the condition on value is sufficient if the variable is compatible
            bestDjComp = value / weight;
            bestSequenceComp = iSequence;
          } else if (value > weightedDj) {
            // the usual condition is checked if the variable is not compatible
            bestDj = value / weight;
            bestSequence = iSequence;
          }
        }
        numberWanted--;
      }
      if (!numberWanted)
        break;
    }
    if (!numberWanted)
      break;
  }
  if (bestSequence >= 0 && model_->getStatus(bestSequence) == ClpSimplex::isFree) {
    printf("Free in %d compat %c dj %g\n",
      bestSequence, modelPE_->isCompatibleCol(bestSequence) ? 'y' : 'n', bestDj);
    bestDjComp = 0.0;
  }
  // choose a compatible or an incompatible row depending on the
  // largest values and on the factor of preference psi_
  if (bestSequenceComp >= 0 && bestDjComp >= psiTmp * bestDj) {

    // record the number of pivots done on compatible variables
    // that would not have been done without positive edge
    if (modelPE_->doStatistics()) {
      if (bestDjComp < bestDj)
        modelPE_->addPriorityPivot();
    }
    bestSequence = bestSequenceComp;
  }
  if (bestSequence >= 0 && psi_ < 1.0 && modelPE_->isCompatibleCol(bestSequence)) {
    modelPE_->isLastPivotCompatible(true);
    modelPE_->addCompatiblePivot();
  } else
    modelPE_->isLastPivotCompatible(false);

  model_->clpMatrix()->setSavedBestSequence(bestSequence);
  if (bestSequence >= 0)
    model_->clpMatrix()->setSavedBestDj(model_->djRegion()[bestSequence]);
  if (sequenceOut >= 0) {
    infeas[sequenceOut] = saveOutInfeasibility;
  }

  modelPE_->updateLastObjectiveValue();
  return bestSequence;
}
/* Save weights - this may initialize weights as well
   This is as parent but may initialize ClpPESimplex
*/
void ClpPEPrimalColumnSteepest::saveWeights(ClpSimplex *model, int mode)
{
  // See if we need to initialize ClpPESimplex
  if (!modelPE_ || model != modelPE_->clpModel() || !modelPE_->checkSize()) {
    delete modelPE_;
    modelPE_ = new ClpPESimplex(model);
  }
  ClpPrimalColumnSteepest::saveWeights(model, mode);
}
// Updates weights - as ordinary but checks for zero moves
void ClpPEPrimalColumnSteepest::updateWeights(CoinIndexedVector *input)
{
  //if (modelPE_->isLastPivotCompatible()) {
  //double theta = modelPE_->model()->theta();
  //if (theta
  //}
  ClpPrimalColumnSteepest::updateWeights(input);
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
