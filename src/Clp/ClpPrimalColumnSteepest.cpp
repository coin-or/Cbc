// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinPragma.hpp"
//#define FAKE_CILK

#include "ClpSimplex.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpFactorization.hpp"
#include "ClpNonLinearCost.hpp"
#include "ClpMessage.hpp"
#include "ClpEventHandler.hpp"
#include "ClpPackedMatrix.hpp"
#include "CoinHelperFunctions.hpp"
#include <stdio.h>
#if ALT_UPDATE_WEIGHTS == 1
extern CoinIndexedVector *altVector[3];
#endif
static void debug1(int iSequence, double value, double weight)
{
  printf("xx %d inf %.20g wt %.20g\n",
    iSequence, value, weight);
}
//#define CLP_DEBUG
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpPrimalColumnSteepest::ClpPrimalColumnSteepest(int mode)
  : ClpPrimalColumnPivot()
  , devex_(0.0)
  , weights_(NULL)
  , infeasible_(NULL)
  , alternateWeights_(NULL)
  , savedWeights_(NULL)
  , reference_(NULL)
  , state_(-1)
  , mode_(mode)
  , infeasibilitiesState_(0)
  , persistence_(normal)
  , numberSwitched_(0)
  , pivotSequence_(-1)
  , savedPivotSequence_(-1)
  , savedSequenceOut_(-1)
  , sizeFactorization_(0)
{
  type_ = 2 + 64 * mode;
}
//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpPrimalColumnSteepest::ClpPrimalColumnSteepest(const ClpPrimalColumnSteepest &rhs)
  : ClpPrimalColumnPivot(rhs)
{
  state_ = rhs.state_;
  mode_ = rhs.mode_;
  infeasibilitiesState_ = rhs.infeasibilitiesState_;
  persistence_ = rhs.persistence_;
  numberSwitched_ = rhs.numberSwitched_;
  model_ = rhs.model_;
  pivotSequence_ = rhs.pivotSequence_;
  savedPivotSequence_ = rhs.savedPivotSequence_;
  savedSequenceOut_ = rhs.savedSequenceOut_;
  sizeFactorization_ = rhs.sizeFactorization_;
  devex_ = rhs.devex_;
  if ((model_ && model_->whatsChanged() & 1) != 0) {
    if (rhs.infeasible_) {
      infeasible_ = new CoinIndexedVector(rhs.infeasible_);
    } else {
      infeasible_ = NULL;
    }
    reference_ = NULL;
    if (rhs.weights_) {
      assert(model_);
      int number = model_->numberRows() + model_->numberColumns();
      assert(number == rhs.model_->numberRows() + rhs.model_->numberColumns());
      weights_ = new double[number];
      CoinMemcpyN(rhs.weights_, number, weights_);
      savedWeights_ = new double[number];
      CoinMemcpyN(rhs.savedWeights_, number, savedWeights_);
      if (mode_ != 1) {
        reference_ = CoinCopyOfArray(rhs.reference_, (number + 31) >> 5);
      }
    } else {
      weights_ = NULL;
      savedWeights_ = NULL;
    }
    if (rhs.alternateWeights_) {
      alternateWeights_ = new CoinIndexedVector(rhs.alternateWeights_);
    } else {
      alternateWeights_ = NULL;
    }
  } else {
    infeasible_ = NULL;
    reference_ = NULL;
    weights_ = NULL;
    savedWeights_ = NULL;
    alternateWeights_ = NULL;
  }
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpPrimalColumnSteepest::~ClpPrimalColumnSteepest()
{
  delete[] weights_;
  delete infeasible_;
  delete alternateWeights_;
  delete[] savedWeights_;
  delete[] reference_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpPrimalColumnSteepest &
ClpPrimalColumnSteepest::operator=(const ClpPrimalColumnSteepest &rhs)
{
  if (this != &rhs) {
    ClpPrimalColumnPivot::operator=(rhs);
    state_ = rhs.state_;
    mode_ = rhs.mode_;
    infeasibilitiesState_ = rhs.infeasibilitiesState_;
    persistence_ = rhs.persistence_;
    numberSwitched_ = rhs.numberSwitched_;
    model_ = rhs.model_;
    pivotSequence_ = rhs.pivotSequence_;
    savedPivotSequence_ = rhs.savedPivotSequence_;
    savedSequenceOut_ = rhs.savedSequenceOut_;
    sizeFactorization_ = rhs.sizeFactorization_;
    devex_ = rhs.devex_;
    delete[] weights_;
    delete[] reference_;
    reference_ = NULL;
    delete infeasible_;
    delete alternateWeights_;
    delete[] savedWeights_;
    savedWeights_ = NULL;
    if (rhs.infeasible_ != NULL) {
      infeasible_ = new CoinIndexedVector(rhs.infeasible_);
    } else {
      infeasible_ = NULL;
    }
    if (rhs.weights_ != NULL) {
      assert(model_);
      int number = model_->numberRows() + model_->numberColumns();
      assert(number == rhs.model_->numberRows() + rhs.model_->numberColumns());
      weights_ = new double[number];
      CoinMemcpyN(rhs.weights_, number, weights_);
      savedWeights_ = new double[number];
      CoinMemcpyN(rhs.savedWeights_, number, savedWeights_);
      if (mode_ != 1) {
        reference_ = CoinCopyOfArray(rhs.reference_, (number + 31) >> 5);
      }
    } else {
      weights_ = NULL;
    }
    if (rhs.alternateWeights_ != NULL) {
      alternateWeights_ = new CoinIndexedVector(rhs.alternateWeights_);
    } else {
      alternateWeights_ = NULL;
    }
  }
  return *this;
}
// These have to match ClpPackedMatrix version
#define TRY_NORM 1.0e-4
#define ADD_ONE 1.0
static void
pivotColumnBit(clpTempInfo &info)
{
  double *COIN_RESTRICT weights = const_cast< double * >(info.lower);
  const unsigned char *COIN_RESTRICT status = info.status;
  double tolerance = info.tolerance;
  double bestDj = info.primalRatio;
  int bestSequence = -1;
  double *COIN_RESTRICT infeas = const_cast< double * >(info.infeas);
  const int *COIN_RESTRICT start = info.which;
  const int *COIN_RESTRICT index = info.index;
  for (int i = start[0]; i < start[1]; i++) {
    int iSequence = index[i];
    double value = infeas[iSequence];
    double weight = weights[iSequence];
    if (value > tolerance) {
      if (value > bestDj * weight) {
        // check flagged variable and correct dj
        if ((status[iSequence] & 64) == 0) {
          bestDj = value / weight;
          bestSequence = iSequence;
        }
      }
    }
  }
  info.primalRatio = bestDj;
  info.numberAdded = bestSequence;
}
// Returns pivot column, -1 if none
/*      The Packed CoinIndexedVector updates has cost updates - for normal LP
	that is just +-weight where a feasibility changed.  It also has
	reduced cost from last iteration in pivot row*/
int ClpPrimalColumnSteepest::pivotColumn(CoinIndexedVector *updates,
  CoinIndexedVector *spareRow1,
  CoinIndexedVector *spareRow2,
  CoinIndexedVector *spareColumn1,
  CoinIndexedVector *spareColumn2)
{
  assert(model_);
  if (model_->nonLinearCost()->lookBothWays() || model_->algorithm() == 2) {
    // Do old way
    updates->expand();
    return pivotColumnOldMethod(updates, spareRow1, spareRow2,
      spareColumn1, spareColumn2);
  }
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
        5 - dantzig
        10 - can go to mini-sprint
     */
    // Look at gub
#if 1
  model_->clpMatrix()->dualExpanded(model_, updates, NULL, 4);
#else
  updates->clear();
  model_->computeDuals(NULL);
#endif
  if (updates->getNumElements() > 1) {
    // would have to have two goes for devex, three for steepest
    anyUpdates = 2;
  } else if (updates->getNumElements()) {
    if (updates->getIndices()[0] == pivotRow && fabs(updates->denseVector()[0]) > 1.0e-6) {
      // reasonable size
      anyUpdates = 1;
      //if (fabs(model_->dualIn())<1.0e-4||fabs(fabs(model_->dualIn())-fabs(updates->denseVector()[0]))>1.0e-5)
      //printf("dualin %g pivot %g\n",model_->dualIn(),updates->denseVector()[0]);
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
  if (switchType == 5) {
    // If known matrix then we will do partial pricing
    if (model_->clpMatrix()->canDoPartialPricing()) {
      pivotSequence_ = -1;
      pivotRow = -1;
      // See if to switch
      int numberRows = model_->numberRows();
      int numberWanted = 10;
      int numberColumns = model_->numberColumns();
      int numberHiddenRows = model_->clpMatrix()->hiddenRows();
      double ratio = static_cast< double >(sizeFactorization_ + numberHiddenRows) / static_cast< double >(numberRows + 2 * numberHiddenRows);
      // Number of dual infeasibilities at last invert
      int numberDual = model_->numberDualInfeasibilities();
      int numberLook = std::min(numberDual, numberColumns / 10);
      if (ratio < 1.0) {
        numberWanted = 100;
        numberLook /= 20;
        numberWanted = std::max(numberWanted, numberLook);
      } else if (ratio < 3.0) {
        numberWanted = 500;
        numberLook /= 15;
        numberWanted = std::max(numberWanted, numberLook);
      } else if (ratio < 4.0 || mode_ == 5) {
        numberWanted = 1000;
        numberLook /= 10;
        numberWanted = std::max(numberWanted, numberLook);
      } else if (mode_ != 5) {
        switchType = 4;
        // initialize
        numberSwitched_++;
        // Make sure will re-do
        delete[] weights_;
        weights_ = NULL;
        model_->computeDuals(NULL);
        saveWeights(model_, 4);
        anyUpdates = 0;
        COIN_DETAIL_PRINT(printf("switching to devex %d nel ratio %g\n", sizeFactorization_, ratio));
      }
      if (switchType == 5) {
        numberLook *= 5; // needs tuning for gub
        if (model_->numberIterations() % 1000 == 0 && model_->logLevel() > 1) {
          COIN_DETAIL_PRINT(printf("numels %d ratio %g wanted %d look %d\n",
            sizeFactorization_, ratio, numberWanted, numberLook));
        }
        // Update duals and row djs
        // Do partial pricing
        return partialPricing(updates, spareRow2,
          numberWanted, numberLook);
      }
    }
  }
  int bestSequence = -1;
  model_->spareIntArray_[3] = -3;
  if (switchType == 5) {
    if (anyUpdates > 0) {
      justDjs(updates, spareRow2,
        spareColumn1, spareColumn2);
    }
  } else if (anyUpdates == 1) {
    if (switchType < 4) {
      // exact etc when can use dj
      djsAndSteepest(updates, spareRow2,
        spareColumn1, spareColumn2);
      if (model_->spareIntArray_[3] > -2) {
        bestSequence = model_->spareIntArray_[3];
        infeasibilitiesState_ = 2;
      } else if (model_->spareIntArray_[3] == -2) {
        redoInfeasibilities();
      }
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
  // everything may have been done if vector copy
  if (infeasibilitiesState_ == 2) {
    // all done
    infeasibilitiesState_ = 1;
    model_->clpMatrix()->setSavedBestSequence(bestSequence);
    if (bestSequence >= 0)
      model_->clpMatrix()->setSavedBestDj(model_->djRegion()[bestSequence]);
    assert(sequenceOut != bestSequence);
    return bestSequence;
  } else if (infeasibilitiesState_ == 1) {
    // need to redo
    //infeasibilitiesState_ = 0;
    redoInfeasibilities();
  }

#ifdef CLP_DEBUG
  alternateWeights_->checkClear();
#endif
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
  // update of duals finished - now do pricing
  // See what sort of pricing
  int numberWanted = 10;
  number = infeasible_->getNumElements();
  int numberColumns = model_->numberColumns();
  if (switchType == 5) {
    pivotSequence_ = -1;
    pivotRow = -1;
    // See if to switch
    int numberRows = model_->numberRows();
    // ratio is done on number of columns here
    //double ratio = static_cast<double> sizeFactorization_/static_cast<double> numberColumns;
    double ratio = static_cast< double >(sizeFactorization_) / static_cast< double >(numberRows);
    //double ratio = static_cast<double> sizeFactorization_/static_cast<double> model_->clpMatrix()->getNumElements();
    if (ratio < 1.0) {
      numberWanted = std::max(100, number / 200);
    } else if (ratio < 2.0 - 1.0) {
      numberWanted = std::max(500, number / 40);
    } else if (ratio < 4.0 - 3.0 || mode_ == 5) {
      numberWanted = std::max(2000, number / 10);
      numberWanted = std::max(numberWanted, numberColumns / 30);
    } else if (mode_ != 5) {
      switchType = 4;
      // initialize
      numberSwitched_++;
      // Make sure will re-do
      delete[] weights_;
      weights_ = NULL;
      saveWeights(model_, 4);
      COIN_DETAIL_PRINT(printf("switching to devex %d nel ratio %g\n", sizeFactorization_, ratio));
    }
    //if (model_->numberIterations()%1000==0)
    //printf("numels %d ratio %g wanted %d\n",sizeFactorization_,ratio,numberWanted);
  }
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
      numberWanted = COIN_INT_MAX - 1;
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
    //if (model_->numberIterations()%1000==0)
    //printf("numels %d ratio %g wanted %d type %d\n",sizeFactorization_,ratio,numberWanted,
    //switchType);
  }

  double bestDj = 1.0e-30;
  bestSequence = -1;
#ifdef CLP_USER_DRIVEN
  // could go parallel?
  model_->eventHandler()->event(ClpEventHandler::beforeChooseIncoming);
#endif

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
#ifdef CLP_DEBUG
  if (model_->numberDualInfeasibilities() == 1)
    printf("** %g %g %g %x %x %d\n", tolerance, model_->dualTolerance(),
      model_->largestDualError(), model_, model_->messageHandler(),
      number);
#endif
  // stop last one coming immediately
  double saveOutInfeasibility = 0.0;
  if (sequenceOut >= 0) {
    saveOutInfeasibility = infeas[sequenceOut];
    infeas[sequenceOut] = 0.0;
  }
  if (model_->factorization()->pivots() && model_->numberPrimalInfeasibilities())
    tolerance = std::max(tolerance, 1.0e-15 * model_->infeasibilityCost());
  tolerance *= tolerance; // as we are using squares

  int iPass;
  // Setup two passes (unless all)
  if (mode_ > 1) {
    int start[4];
    start[1] = number;
    start[2] = 0;
    double dstart = static_cast< double >(number) * model_->randomNumberGenerator()->randomDouble();
    start[0] = static_cast< int >(dstart);
    start[3] = start[0];
    for (iPass = 0; iPass < 2; iPass++) {
      int end = start[2 * iPass + 1];
      if (switchType < 5) {
        for (i = start[2 * iPass]; i < end; i++) {
          iSequence = index[i];
          double value = infeas[iSequence];
          double weight = weights_[iSequence];
          //assert (weight);
          if (value > tolerance) {
            //weight=1.0;
            if (value > bestDj * weight) {
              // check flagged variable and correct dj
              if (!model_->flagged(iSequence)) {
                bestDj = value / weight;
                bestSequence = iSequence;
              } else {
                // just to make sure we don't exit before got something
                numberWanted++;
              }
            }
            numberWanted--;
          }
          if (!numberWanted)
            break;
        }
      } else {
        // Dantzig
        for (i = start[2 * iPass]; i < end; i++) {
          iSequence = index[i];
          double value = infeas[iSequence];
          if (value > tolerance) {
            if (value > bestDj) {
              // check flagged variable and correct dj
              if (!model_->flagged(iSequence)) {
                bestDj = value;
                bestSequence = iSequence;
              } else {
                // just to make sure we don't exit before got something
                numberWanted++;
              }
            }
            numberWanted--;
          }
          if (!numberWanted)
            break;
        }
      }
      if (!numberWanted)
        break;
    }
  } else {
    const int numberThreads = 1;
    if (0) {
      int iSequence;
      double value;
      double weight;
      iSequence = 34841;
      value = infeas[iSequence];
      weight = weights_[iSequence];
      debug1(iSequence, value, weight);
      iSequence = 34845;
      value = infeas[iSequence];
      weight = weights_[iSequence];
      debug1(iSequence, value, weight);
    }
    clpTempInfo info[1];
    int start_lite[2];
    int chunk0 = (number + numberThreads - 1) / numberThreads;
    int n0 = 0;
    for (int i = 0; i < numberThreads; i++) {
      int *startX = start_lite + 2 * i;
      info[i].primalRatio = bestDj;
      info[i].lower = weights_;
      info[i].infeas = infeas;
      info[i].index = index;
      info[i].status = const_cast< unsigned char * >(model_->statusArray());
      info[i].which = startX;
      info[i].tolerance = tolerance;
      startX[0] = n0;
      startX[1] = std::min(n0 + chunk0, number);
      n0 += chunk0;
    }
    if (numberThreads == 1) {
      pivotColumnBit(info[0]);
    } else {
      for (int i = 0; i < numberThreads; i++) {
        cilk_spawn pivotColumnBit(info[i]);
      }
      cilk_sync;
    }
    for (int i = 0; i < numberThreads; i++) {
      double bestDjX = info[i].primalRatio;
      if (bestDjX > bestDj) {
        bestDj = bestDjX;
        bestSequence = info[i].numberAdded;
      }
    }
  }
  //double largestWeight=0.0;
  //double smallestWeight=1.0e100;
  model_->clpMatrix()->setSavedBestSequence(bestSequence);
  if (bestSequence >= 0)
    model_->clpMatrix()->setSavedBestDj(model_->djRegion()[bestSequence]);
  if (sequenceOut >= 0) {
    infeas[sequenceOut] = saveOutInfeasibility;
  }
  /*if (model_->numberIterations()%100==0)
       printf("%d best %g\n",bestSequence,bestDj);*/

#ifdef CLP_USER_DRIVEN
  // swap when working
  struct {
    int bestSequence;
    double bestDj;
  } tempInfo;
  tempInfo.bestDj = bestDj;
  tempInfo.bestSequence = bestSequence;
  model_->eventHandler()->eventWithInfo(ClpEventHandler::afterChooseIncoming,
    reinterpret_cast< void * >(&tempInfo));
#endif
#ifndef NDEBUG
  if (bestSequence >= 0) {
    if (model_->getStatus(bestSequence) == ClpSimplex::atLowerBound)
      assert(model_->reducedCost(bestSequence) < 0.0);
    if (model_->getStatus(bestSequence) == ClpSimplex::atUpperBound) {
      assert(model_->reducedCost(bestSequence) > 0.0);
    }
  }
#endif
#ifdef ALT_UPDATE_WEIGHTSz
  printf("weights");
  for (int i = 0; i < model_->numberColumns() + model_->numberRows(); i++) {
    if (model_->getColumnStatus(i) == ClpSimplex::isFixed || model_->getColumnStatus(i) == ClpSimplex::basic)
      weights_[i] = 1.0;
    printf(" %g", weights_[i]);
    if ((i % 10) == 9)
      printf("\n");
  }
  printf("\n");
#endif
  return bestSequence;
}
// Just update djs
void ClpPrimalColumnSteepest::justDjs(CoinIndexedVector *updates,
  CoinIndexedVector *spareRow2,
  CoinIndexedVector *spareColumn1,
  CoinIndexedVector *spareColumn2)
{
  int iSection, j;
  int number = 0;
  int *index;
  double *updateBy;
  double *reducedCost;
  double tolerance = model_->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = std::min(1.0e-2, model_->largestDualError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance + error;
  int pivotRow = model_->pivotRow();
  double *infeas = infeasible_->denseVector();
  //updates->scanAndPack();
  model_->factorization()->updateColumnTranspose(spareRow2, updates);

  // put row of tableau in rowArray and columnArray (packed mode)
  model_->clpMatrix()->transposeTimes(model_, -1.0,
    updates, spareColumn2, spareColumn1);
  // normal
  for (iSection = 0; iSection < 2; iSection++) {

    reducedCost = model_->djRegion(iSection);
    int addSequence;
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
    double slack_multiplier;
#endif

    if (!iSection) {
      number = updates->getNumElements();
      index = updates->getIndices();
      updateBy = updates->denseVector();
      addSequence = model_->numberColumns();
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
      slack_multiplier = CLP_PRIMAL_SLACK_MULTIPLIER;
#endif
    } else {
      number = spareColumn1->getNumElements();
      index = spareColumn1->getIndices();
      updateBy = spareColumn1->denseVector();
      addSequence = 0;
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
      slack_multiplier = 1.0;
#endif
    }

    for (j = 0; j < number; j++) {
      int iSequence = index[j];
      double value = reducedCost[iSequence];
      value -= updateBy[j];
      updateBy[j] = 0.0;
      reducedCost[iSequence] = value;
      ClpSimplex::Status status = model_->getStatus(iSequence + addSequence);

      switch (status) {

      case ClpSimplex::basic:
        infeasible_->zero(iSequence + addSequence);
      case ClpSimplex::isFixed:
        break;
      case ClpSimplex::isFree:
      case ClpSimplex::superBasic:
        if (fabs(value) > FREE_ACCEPT * tolerance) {
          // we are going to bias towards free (but only if reasonable)
          value *= FREE_BIAS;
          // store square in list
          if (infeas[iSequence + addSequence])
            infeas[iSequence + addSequence] = value * value; // already there
          else
            infeasible_->quickAdd(iSequence + addSequence, value * value);
        } else {
          infeasible_->zero(iSequence + addSequence);
        }
        break;
      case ClpSimplex::atUpperBound:
        iSequence += addSequence;
        if (value > tolerance) {
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
          value *= value * slack_multiplier;
#else
          value *= value;
#endif
          // store square in list
          if (infeas[iSequence])
            infeas[iSequence] = value; // already there
          else
            infeasible_->quickAdd(iSequence, value);
        } else {
          infeasible_->zero(iSequence);
        }
        break;
      case ClpSimplex::atLowerBound:
        iSequence += addSequence;
        if (value < -tolerance) {
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
          value *= value * slack_multiplier;
#else
          value *= value;
#endif
          // store square in list
          if (infeas[iSequence])
            infeas[iSequence] = value; // already there
          else
            infeasible_->quickAdd(iSequence, value);
        } else {
          infeasible_->zero(iSequence);
        }
      }
    }
  }
  updates->setNumElements(0);
  spareColumn1->setNumElements(0);
  if (pivotRow >= 0) {
    // make sure infeasibility on incoming is 0.0
    int sequenceIn = model_->sequenceIn();
    infeasible_->zero(sequenceIn);
  }
}
// Update djs, weights for Devex
void ClpPrimalColumnSteepest::djsAndDevex(CoinIndexedVector *updates,
  CoinIndexedVector *spareRow2,
  CoinIndexedVector *spareColumn1,
  CoinIndexedVector *spareColumn2)
{
  int j;
  int number = 0;
  int *index;
  double *updateBy;
  double *reducedCost;
  double tolerance = model_->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = std::min(1.0e-2, model_->largestDualError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance + error;
  // for weights update we use pivotSequence
  // unset in case sub flip
  assert(pivotSequence_ >= 0);
  assert(model_->pivotVariable()[pivotSequence_] == model_->sequenceIn());
  pivotSequence_ = -1;
  double *infeas = infeasible_->denseVector();
  //updates->scanAndPack();
  model_->factorization()->updateColumnTranspose(spareRow2, updates);
  // and we can see if reference
  //double referenceIn = 0.0;
  int sequenceIn = model_->sequenceIn();
  //if (mode_ != 1 && reference(sequenceIn))
  //   referenceIn = 1.0;
  // save outgoing weight round update
  double outgoingWeight = 0.0;
  int sequenceOut = model_->sequenceOut();
  if (sequenceOut >= 0)
    outgoingWeight = weights_[sequenceOut];

  double scaleFactor = 1.0 / updates->denseVector()[0]; // as formula is with 1.0
  // put row of tableau in rowArray and columnArray (packed mode)
  model_->clpMatrix()->transposeTimes(model_, -1.0,
    updates, spareColumn2, spareColumn1);
  // update weights
  double *weight;
  int numberColumns = model_->numberColumns();
  // rows
  reducedCost = model_->djRegion(0);
  int addSequence = model_->numberColumns();
  ;

  number = updates->getNumElements();
  index = updates->getIndices();
  updateBy = updates->denseVector();
  weight = weights_ + numberColumns;
  // Devex
  for (j = 0; j < number; j++) {
    double thisWeight;
    double pivot;
    double value3;
    int iSequence = index[j];
    double value = reducedCost[iSequence];
    double value2 = updateBy[j];
    updateBy[j] = 0.0;
    value -= value2;
    reducedCost[iSequence] = value;
    ClpSimplex::Status status = model_->getStatus(iSequence + addSequence);

    switch (status) {

    case ClpSimplex::basic:
      infeasible_->zero(iSequence + addSequence);
    case ClpSimplex::isFixed:
      break;
    case ClpSimplex::isFree:
    case ClpSimplex::superBasic:
      thisWeight = weight[iSequence];
      // row has -1
      pivot = value2 * scaleFactor;
      value3 = pivot * pivot * devex_;
      if (reference(iSequence + numberColumns))
        value3 += 1.0;
      weight[iSequence] = std::max(0.99 * thisWeight, value3);
      if (fabs(value) > FREE_ACCEPT * tolerance) {
        // we are going to bias towards free (but only if reasonable)
        value *= FREE_BIAS;
        // store square in list
        if (infeas[iSequence + addSequence])
          infeas[iSequence + addSequence] = value * value; // already there
        else
          infeasible_->quickAdd(iSequence + addSequence, value * value);
      } else {
        infeasible_->zero(iSequence + addSequence);
      }
      break;
    case ClpSimplex::atUpperBound:
      thisWeight = weight[iSequence];
      // row has -1
      pivot = value2 * scaleFactor;
      value3 = pivot * pivot * devex_;
      if (reference(iSequence + numberColumns))
        value3 += 1.0;
      weight[iSequence] = std::max(0.99 * thisWeight, value3);
      iSequence += addSequence;
      if (value > tolerance) {
        // store square in list
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
        value *= value * CLP_PRIMAL_SLACK_MULTIPLIER;
#else
        value *= value;
#endif
        if (infeas[iSequence])
          infeas[iSequence] = value; // already there
        else
          infeasible_->quickAdd(iSequence, value);
      } else {
        infeasible_->zero(iSequence);
      }
      break;
    case ClpSimplex::atLowerBound:
      thisWeight = weight[iSequence];
      // row has -1
      pivot = value2 * scaleFactor;
      value3 = pivot * pivot * devex_;
      if (reference(iSequence + numberColumns))
        value3 += 1.0;
      weight[iSequence] = std::max(0.99 * thisWeight, value3);
      iSequence += addSequence;
      if (value < -tolerance) {
        // store square in list
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
        value *= value * CLP_PRIMAL_SLACK_MULTIPLIER;
#else
        value *= value;
#endif
        if (infeas[iSequence])
          infeas[iSequence] = value; // already there
        else
          infeasible_->quickAdd(iSequence, value);
      } else {
        infeasible_->zero(iSequence);
      }
    }
  }

  // columns
  weight = weights_;

  scaleFactor = -scaleFactor;
  reducedCost = model_->djRegion(1);
  number = spareColumn1->getNumElements();
  index = spareColumn1->getIndices();
  updateBy = spareColumn1->denseVector();

  // Devex

  for (j = 0; j < number; j++) {
    double thisWeight;
    double pivot;
    double value3;
    int iSequence = index[j];
    double value = reducedCost[iSequence];
    double value2 = updateBy[j];
    value -= value2;
    updateBy[j] = 0.0;
    reducedCost[iSequence] = value;
    ClpSimplex::Status status = model_->getStatus(iSequence);

    switch (status) {

    case ClpSimplex::basic:
      infeasible_->zero(iSequence);
    case ClpSimplex::isFixed:
      break;
    case ClpSimplex::isFree:
    case ClpSimplex::superBasic:
      thisWeight = weight[iSequence];
      // row has -1
      pivot = value2 * scaleFactor;
      value3 = pivot * pivot * devex_;
      if (reference(iSequence))
        value3 += 1.0;
      weight[iSequence] = std::max(0.99 * thisWeight, value3);
      if (fabs(value) > FREE_ACCEPT * tolerance) {
        // we are going to bias towards free (but only if reasonable)
        value *= FREE_BIAS;
        // store square in list
        if (infeas[iSequence])
          infeas[iSequence] = value * value; // already there
        else
          infeasible_->quickAdd(iSequence, value * value);
      } else {
        infeasible_->zero(iSequence);
      }
      break;
    case ClpSimplex::atUpperBound:
      thisWeight = weight[iSequence];
      // row has -1
      pivot = value2 * scaleFactor;
      value3 = pivot * pivot * devex_;
      if (reference(iSequence))
        value3 += 1.0;
      weight[iSequence] = std::max(0.99 * thisWeight, value3);
      if (value > tolerance) {
        // store square in list
        if (infeas[iSequence])
          infeas[iSequence] = value * value; // already there
        else
          infeasible_->quickAdd(iSequence, value * value);
      } else {
        infeasible_->zero(iSequence);
      }
      break;
    case ClpSimplex::atLowerBound:
      thisWeight = weight[iSequence];
      // row has -1
      pivot = value2 * scaleFactor;
      value3 = pivot * pivot * devex_;
      if (reference(iSequence))
        value3 += 1.0;
      weight[iSequence] = std::max(0.99 * thisWeight, value3);
      if (value < -tolerance) {
        // store square in list
        if (infeas[iSequence])
          infeas[iSequence] = value * value; // already there
        else
          infeasible_->quickAdd(iSequence, value * value);
      } else {
        infeasible_->zero(iSequence);
      }
    }
  }
  // restore outgoing weight
  if (sequenceOut >= 0)
    weights_[sequenceOut] = outgoingWeight;
  // make sure infeasibility on incoming is 0.0
  infeasible_->zero(sequenceIn);
  spareRow2->setNumElements(0);
  //#define SOME_DEBUG_1
#ifdef SOME_DEBUG_1
  // check for accuracy
  int iCheck = 892;
  //printf("weight for iCheck is %g\n",weights_[iCheck]);
  int numberRows = model_->numberRows();
  //int numberColumns = model_->numberColumns();
  for (iCheck = 0; iCheck < numberRows + numberColumns; iCheck++) {
    if (model_->getStatus(iCheck) != ClpSimplex::basic && model_->getStatus(iCheck) != ClpSimplex::isFixed)
      checkAccuracy(iCheck, 1.0e-1, updates, spareRow2);
  }
#endif
  updates->setNumElements(0);
  spareColumn1->setNumElements(0);
}
// Update djs, weights for Steepest
void ClpPrimalColumnSteepest::djsAndSteepest(CoinIndexedVector *updates,
  CoinIndexedVector *spareRow2,
  CoinIndexedVector *spareColumn1,
  CoinIndexedVector *spareColumn2)
{
  int j;
  int number = 0;
  int *index;
  double *updateBy;
  double *reducedCost;
  double tolerance = model_->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = std::min(1.0e-2, model_->largestDualError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance + error;
  // for weights update we use pivotSequence
  // unset in case sub flip
  assert(pivotSequence_ >= 0);
  assert(model_->pivotVariable()[pivotSequence_] == model_->sequenceIn());
  double *infeas = infeasible_->denseVector();
  double scaleFactor = 1.0 / updates->denseVector()[0]; // as formula is with 1.0
  assert(updates->getIndices()[0] == pivotSequence_);
  pivotSequence_ = -1;
#if 1
  //updates->scanAndPack();
  model_->factorization()->updateColumnTranspose(spareRow2, updates);
  //alternateWeights_->scanAndPack();
#if ALT_UPDATE_WEIGHTS != 2
  model_->factorization()->updateColumnTranspose(spareRow2,
    alternateWeights_);
#elif ALT_UPDATE_WEIGHTS == 1
  if (altVector[1]) {
    int numberRows = model_->numberRows();
    double *work1 = altVector[1]->denseVector();
    double *worka = alternateWeights_->denseVector();
    int iRow = -1;
    double diff = 1.0e-8;
    for (int i = 0; i < numberRows; i++) {
      double dd = std::max(fabs(work1[i]), fabs(worka[i]));
      double d = fabs(work1[i] - worka[i]);
      if (dd > 1.0e-6 && d > diff * dd) {
        diff = d / dd;
        iRow = i;
      }
    }
    if (iRow >= 0)
      printf("largest difference of %g (%g,%g) on row %d\n",
        diff, work1[iRow], worka[iRow], iRow);
  }
#endif
#else
  model_->factorization()->updateTwoColumnsTranspose(spareRow2, updates, alternateWeights_);
#endif
  // and we can see if reference
  int sequenceIn = model_->sequenceIn();
  double referenceIn;
  if (mode_ != 1) {
    if (reference(sequenceIn))
      referenceIn = 1.0;
    else
      referenceIn = 0.0;
  } else {
    referenceIn = -1.0;
  }
  // save outgoing weight round update
  double outgoingWeight = 0.0;
  int sequenceOut = model_->sequenceOut();
  if (sequenceOut >= 0)
    outgoingWeight = weights_[sequenceOut];
  // update row weights here so we can scale alternateWeights_
  // update weights
  double *weight;
  double *other = alternateWeights_->denseVector();
  int numberColumns = model_->numberColumns();
  // if if (model_->clpMatrix()->canCombine(model_, pi1)) no need to index
  // rows
  reducedCost = model_->djRegion(0);
  int addSequence = model_->numberColumns();
  ;

  number = updates->getNumElements();
  index = updates->getIndices();
  updateBy = updates->denseVector();
  weight = weights_ + numberColumns;
#ifdef CLP_USER_DRIVEN
  model_->eventHandler()->eventWithInfo(ClpEventHandler::beforeChooseIncoming, updates);
#endif

  for (j = 0; j < number; j++) {
    double thisWeight;
    double pivot;
    double modification;
    double pivotSquared;
    int iSequence = index[j];
    double value2 = updateBy[j];
    ClpSimplex::Status status = model_->getStatus(iSequence + addSequence);
    double value;

    switch (status) {

    case ClpSimplex::basic:
      infeasible_->zero(iSequence + addSequence);
      reducedCost[iSequence] = 0.0;
    case ClpSimplex::isFixed:
      break;
    case ClpSimplex::isFree:
    case ClpSimplex::superBasic:
      value = reducedCost[iSequence] - value2;
      modification = other[iSequence];
      thisWeight = weight[iSequence];
      // row has -1
      pivot = value2 * scaleFactor;
      pivotSquared = pivot * pivot;

      thisWeight += pivotSquared * devex_ + pivot * modification;
      reducedCost[iSequence] = value;
      if (thisWeight < TRY_NORM) {
        if (mode_ == 1) {
          // steepest
          thisWeight = std::max(TRY_NORM, ADD_ONE + pivotSquared);
        } else {
          // exact
          thisWeight = referenceIn * pivotSquared;
          if (reference(iSequence + numberColumns))
            thisWeight += 1.0;
          thisWeight = std::max(thisWeight, TRY_NORM);
        }
      }
      weight[iSequence] = thisWeight;
      if (fabs(value) > FREE_ACCEPT * tolerance) {
        // we are going to bias towards free (but only if reasonable)
        value *= FREE_BIAS;
        // store square in list
        if (infeas[iSequence + addSequence])
          infeas[iSequence + addSequence] = value * value; // already there
        else
          infeasible_->quickAdd(iSequence + addSequence, value * value);
      } else {
        infeasible_->zero(iSequence + addSequence);
      }
      break;
    case ClpSimplex::atUpperBound:
      value = reducedCost[iSequence] - value2;
      modification = other[iSequence];
      thisWeight = weight[iSequence];
      // row has -1
      pivot = value2 * scaleFactor;
      pivotSquared = pivot * pivot;

      thisWeight += pivotSquared * devex_ + pivot * modification;
      reducedCost[iSequence] = value;
      if (thisWeight < TRY_NORM) {
        if (mode_ == 1) {
          // steepest
          thisWeight = std::max(TRY_NORM, ADD_ONE + pivotSquared);
        } else {
          // exact
          thisWeight = referenceIn * pivotSquared;
          if (reference(iSequence + numberColumns))
            thisWeight += 1.0;
          thisWeight = std::max(thisWeight, TRY_NORM);
        }
      }
      weight[iSequence] = thisWeight;
      iSequence += addSequence;
      if (value > tolerance) {
        // store square in list
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
        value *= value * CLP_PRIMAL_SLACK_MULTIPLIER;
#else
        value *= value;
#endif
        if (infeas[iSequence])
          infeas[iSequence] = value; // already there
        else
          infeasible_->quickAdd(iSequence, value);
      } else {
        infeasible_->zero(iSequence);
      }
      break;
    case ClpSimplex::atLowerBound:
      value = reducedCost[iSequence] - value2;
      modification = other[iSequence];
      thisWeight = weight[iSequence];
      // row has -1
      pivot = value2 * scaleFactor;
      pivotSquared = pivot * pivot;

      thisWeight += pivotSquared * devex_ + pivot * modification;
      reducedCost[iSequence] = value;
      if (thisWeight < TRY_NORM) {
        if (mode_ == 1) {
          // steepest
          thisWeight = std::max(TRY_NORM, ADD_ONE + pivotSquared);
        } else {
          // exact
          thisWeight = referenceIn * pivotSquared;
          if (reference(iSequence + numberColumns))
            thisWeight += 1.0;
          thisWeight = std::max(thisWeight, TRY_NORM);
        }
      }
      weight[iSequence] = thisWeight;
      iSequence += addSequence;
      if (value < -tolerance) {
        // store square in list
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
        value *= value * CLP_PRIMAL_SLACK_MULTIPLIER;
#else
        value *= value;
#endif
        if (infeas[iSequence])
          infeas[iSequence] = value; // already there
        else
          infeasible_->quickAdd(iSequence, value);
      } else {
        infeasible_->zero(iSequence);
      }
    }
  }
#ifdef CLP_USER_DRIVEN
  // could go parallel?
  //model_->eventHandler()->eventWithInfo(ClpEventHandler::beforeChooseIncoming,updates);
#endif
  // put row of tableau in rowArray and columnArray (packed)
  // get subset which have nonzero tableau elements
  int returnCode = transposeTimes2(updates, spareColumn1,
    alternateWeights_, spareColumn2, spareRow2,
    -scaleFactor);
  // zero updateBy
  CoinZeroN(updateBy, number);
  alternateWeights_->clear();
  // columns
  assert(scaleFactor);
  number = spareColumn1->getNumElements();
  index = spareColumn1->getIndices();
  updateBy = spareColumn1->denseVector();
  if (returnCode != 2 && infeasibilitiesState_) {
    //spareColumn1->clear();
    redoInfeasibilities();
  }
  if (returnCode == 1) {
    // most work already done
    for (j = 0; j < number; j++) {
      int iSequence = index[j];
      double value = updateBy[j];
      if (value) {
        updateBy[j] = 0.0;
        infeasible_->quickAdd(iSequence, value);
      } else {
        infeasible_->zero(iSequence);
      }
    }
  } else if (returnCode == 0) {
    reducedCost = model_->djRegion(1);

    for (j = 0; j < number; j++) {
      int iSequence = index[j];
      double value = reducedCost[iSequence];
      double value2 = updateBy[j];
      updateBy[j] = 0.0;
      value -= value2;
      reducedCost[iSequence] = value;
      ClpSimplex::Status status = model_->getStatus(iSequence);

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
          if (infeas[iSequence])
            infeas[iSequence] = value * value; // already there
          else
            infeasible_->quickAdd(iSequence, value * value);
        } else {
          infeasible_->zero(iSequence);
        }
        break;
      case ClpSimplex::atUpperBound:
        if (value > tolerance) {
          // store square in list
          if (infeas[iSequence])
            infeas[iSequence] = value * value; // already there
          else
            infeasible_->quickAdd(iSequence, value * value);
        } else {
          infeasible_->zero(iSequence);
        }
        break;
      case ClpSimplex::atLowerBound:
        if (value < -tolerance) {
          // store square in list
          if (infeas[iSequence])
            infeas[iSequence] = value * value; // already there
          else
            infeasible_->quickAdd(iSequence, value * value);
        } else {
          infeasible_->zero(iSequence);
        }
      }
    }
  } else {
    //assert(!number);
  }
  // restore outgoing weight
  if (sequenceOut >= 0) {
    //#define GROW_REFERENCE
#ifdef GROW_REFERENCE
    if (!reference(sequenceOut)) {
      outgoingWeight += 1.0;
      setReference(sequenceOut, true);
    }
#endif
    weights_[sequenceOut] = outgoingWeight;
    //if (model_->getStatus(sequenceOut) != ClpSimplex::basic &&
    //  model_->getStatus(sequenceOut) != ClpSimplex::isFixed)
    //checkAccuracy(sequenceOut, 1.0e-1, updates, spareRow2);
  }
  // make sure infeasibility on incoming is 0.0
  infeasible_->zero(sequenceIn);
  spareColumn2->setNumElements(0);
  //#define SOME_DEBUG_1
#ifdef SOME_DEBUG_1
  // check for accuracy
  int iCheck = 892;
  //printf("weight for iCheck is %g\n",weights_[iCheck]);
  int numberRows = model_->numberRows();
  //int numberColumns = model_->numberColumns();
  for (iCheck = 0; iCheck < numberRows + numberColumns; iCheck++) {
    if (model_->getStatus(iCheck) != ClpSimplex::basic && model_->getStatus(iCheck) != ClpSimplex::isFixed)
      checkAccuracy(iCheck, 1.0e-1, updates, spareRow2);
  }
#endif
  updates->setNumElements(0);
  spareColumn1->setNumElements(0);
}
// Update djs, weights for Devex
void ClpPrimalColumnSteepest::djsAndDevex2(CoinIndexedVector *updates,
  CoinIndexedVector *spareRow2,
  CoinIndexedVector *spareColumn1,
  CoinIndexedVector *spareColumn2)
{
  int iSection, j;
  int number = 0;
  int *index;
  double *updateBy;
  double *reducedCost;
  // dj could be very small (or even zero - take care)
  double dj = model_->dualIn();
  double tolerance = model_->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = std::min(1.0e-2, model_->largestDualError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance + error;
  int pivotRow = model_->pivotRow();
  double *infeas = infeasible_->denseVector();
  //updates->scanAndPack();
  model_->factorization()->updateColumnTranspose(spareRow2, updates);

  // put row of tableau in rowArray and columnArray
  model_->clpMatrix()->transposeTimes(model_, -1.0,
    updates, spareColumn2, spareColumn1);
  // normal
  for (iSection = 0; iSection < 2; iSection++) {

    reducedCost = model_->djRegion(iSection);
    int addSequence;
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
    double slack_multiplier;
#endif

    if (!iSection) {
      number = updates->getNumElements();
      index = updates->getIndices();
      updateBy = updates->denseVector();
      addSequence = model_->numberColumns();
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
      slack_multiplier = CLP_PRIMAL_SLACK_MULTIPLIER;
#endif
    } else {
      number = spareColumn1->getNumElements();
      index = spareColumn1->getIndices();
      updateBy = spareColumn1->denseVector();
      addSequence = 0;
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
      slack_multiplier = 1.0;
#endif
    }

    for (j = 0; j < number; j++) {
      int iSequence = index[j];
      double value = reducedCost[iSequence];
      value -= updateBy[j];
      updateBy[j] = 0.0;
      reducedCost[iSequence] = value;
      ClpSimplex::Status status = model_->getStatus(iSequence + addSequence);

      switch (status) {

      case ClpSimplex::basic:
        infeasible_->zero(iSequence + addSequence);
      case ClpSimplex::isFixed:
        break;
      case ClpSimplex::isFree:
      case ClpSimplex::superBasic:
        if (fabs(value) > FREE_ACCEPT * tolerance) {
          // we are going to bias towards free (but only if reasonable)
          value *= FREE_BIAS;
          // store square in list
          if (infeas[iSequence + addSequence])
            infeas[iSequence + addSequence] = value * value; // already there
          else
            infeasible_->quickAdd(iSequence + addSequence, value * value);
        } else {
          infeasible_->zero(iSequence + addSequence);
        }
        break;
      case ClpSimplex::atUpperBound:
        iSequence += addSequence;
        if (value > tolerance) {
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
          value *= value * slack_multiplier;
#else
          value *= value;
#endif
          // store square in list
          if (infeas[iSequence])
            infeas[iSequence] = value; // already there
          else
            infeasible_->quickAdd(iSequence, value);
        } else {
          infeasible_->zero(iSequence);
        }
        break;
      case ClpSimplex::atLowerBound:
        iSequence += addSequence;
        if (value < -tolerance) {
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
          value *= value * slack_multiplier;
#else
          value *= value;
#endif
          // store square in list
          if (infeas[iSequence])
            infeas[iSequence] = value; // already there
          else
            infeasible_->quickAdd(iSequence, value);
        } else {
          infeasible_->zero(iSequence);
        }
      }
    }
  }
  // They are empty
  updates->setNumElements(0);
  spareColumn1->setNumElements(0);
  // make sure infeasibility on incoming is 0.0
  int sequenceIn = model_->sequenceIn();
  infeasible_->zero(sequenceIn);
  // for weights update we use pivotSequence
  if (pivotSequence_ >= 0) {
    pivotRow = pivotSequence_;
    // unset in case sub flip
    pivotSequence_ = -1;
    // make sure infeasibility on incoming is 0.0
    const int *pivotVariable = model_->pivotVariable();
    sequenceIn = pivotVariable[pivotRow];
    infeasible_->zero(sequenceIn);
    // and we can see if reference
    //double referenceIn = 0.0;
    //if (mode_ != 1 && reference(sequenceIn))
    //   referenceIn = 1.0;
    // save outgoing weight round update
    double outgoingWeight = 0.0;
    int sequenceOut = model_->sequenceOut();
    if (sequenceOut >= 0)
      outgoingWeight = weights_[sequenceOut];
    // update weights
    updates->setNumElements(0);
    spareColumn1->setNumElements(0);
    // might as well set dj to 1
    dj = 1.0;
    updates->insert(pivotRow, -dj);
    model_->factorization()->updateColumnTranspose(spareRow2, updates);
    // put row of tableau in rowArray and columnArray
    model_->clpMatrix()->transposeTimes(model_, -1.0,
      updates, spareColumn2, spareColumn1);
    double *weight;
    int numberColumns = model_->numberColumns();
    // rows
    number = updates->getNumElements();
    index = updates->getIndices();
    updateBy = updates->denseVector();
    weight = weights_ + numberColumns;

    assert(devex_ > 0.0);
    for (j = 0; j < number; j++) {
      int iSequence = index[j];
      double thisWeight = weight[iSequence];
      // row has -1
      double pivot = -updateBy[iSequence];
      updateBy[iSequence] = 0.0;
      double value = pivot * pivot * devex_;
      if (reference(iSequence + numberColumns))
        value += 1.0;
      weight[iSequence] = std::max(0.99 * thisWeight, value);
    }

    // columns
    weight = weights_;

    number = spareColumn1->getNumElements();
    index = spareColumn1->getIndices();
    updateBy = spareColumn1->denseVector();
    for (j = 0; j < number; j++) {
      int iSequence = index[j];
      double thisWeight = weight[iSequence];
      // row has -1
      double pivot = updateBy[iSequence];
      updateBy[iSequence] = 0.0;
      double value = pivot * pivot * devex_;
      if (reference(iSequence))
        value += 1.0;
      weight[iSequence] = std::max(0.99 * thisWeight, value);
    }
    // restore outgoing weight
    if (sequenceOut >= 0)
      weights_[sequenceOut] = outgoingWeight;
    spareColumn2->setNumElements(0);
    //#define SOME_DEBUG_1
#ifdef SOME_DEBUG_1
    // check for accuracy
    int iCheck = 892;
    //printf("weight for iCheck is %g\n",weights_[iCheck]);
    int numberRows = model_->numberRows();
    //int numberColumns = model_->numberColumns();
    for (iCheck = 0; iCheck < numberRows + numberColumns; iCheck++) {
      if (model_->getStatus(iCheck) != ClpSimplex::basic && model_->getStatus(iCheck) != ClpSimplex::isFixed)
        checkAccuracy(iCheck, 1.0e-1, updates, spareRow2);
    }
#endif
    updates->setNumElements(0);
    spareColumn1->setNumElements(0);
  }
}
// Update djs, weights for Steepest
void ClpPrimalColumnSteepest::djsAndSteepest2(CoinIndexedVector *updates,
  CoinIndexedVector *spareRow2,
  CoinIndexedVector *spareColumn1,
  CoinIndexedVector *spareColumn2)
{
  int iSection, j;
  int number = 0;
  int *index;
  double *updateBy;
  double *reducedCost;
  // dj could be very small (or even zero - take care)
  double dj = model_->dualIn();
  double tolerance = model_->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = std::min(1.0e-2, model_->largestDualError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance + error;
  int pivotRow = model_->pivotRow();
  double *infeas = infeasible_->denseVector();
  //updates->scanAndPack();
  model_->factorization()->updateColumnTranspose(spareRow2, updates);

  // put row of tableau in rowArray and columnArray
  model_->clpMatrix()->transposeTimes(model_, -1.0,
    updates, spareColumn2, spareColumn1);
#ifdef CLP_USER_DRIVEN
  model_->eventHandler()->eventWithInfo(ClpEventHandler::beforeChooseIncoming, updates);
#endif
  // normal
  for (iSection = 0; iSection < 2; iSection++) {

    reducedCost = model_->djRegion(iSection);
    int addSequence;
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
    double slack_multiplier;
#endif

    if (!iSection) {
      number = updates->getNumElements();
      index = updates->getIndices();
      updateBy = updates->denseVector();
      addSequence = model_->numberColumns();
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
      slack_multiplier = CLP_PRIMAL_SLACK_MULTIPLIER;
#endif
    } else {
      number = spareColumn1->getNumElements();
      index = spareColumn1->getIndices();
      updateBy = spareColumn1->denseVector();
      addSequence = 0;
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
      slack_multiplier = 1.0;
#endif
    }

    for (j = 0; j < number; j++) {
      int iSequence = index[j];
      double value = reducedCost[iSequence];
      value -= updateBy[j];
      updateBy[j] = 0.0;
      reducedCost[iSequence] = value;
      ClpSimplex::Status status = model_->getStatus(iSequence + addSequence);

      switch (status) {

      case ClpSimplex::basic:
        infeasible_->zero(iSequence + addSequence);
      case ClpSimplex::isFixed:
        break;
      case ClpSimplex::isFree:
      case ClpSimplex::superBasic:
        if (fabs(value) > FREE_ACCEPT * tolerance) {
          // we are going to bias towards free (but only if reasonable)
          value *= FREE_BIAS;
          // store square in list
          if (infeas[iSequence + addSequence])
            infeas[iSequence + addSequence] = value * value; // already there
          else
            infeasible_->quickAdd(iSequence + addSequence, value * value);
        } else {
          infeasible_->zero(iSequence + addSequence);
        }
        break;
      case ClpSimplex::atUpperBound:
        iSequence += addSequence;
        if (value > tolerance) {
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
          value *= value * slack_multiplier;
#else
          value *= value;
#endif
          // store square in list
          if (infeas[iSequence])
            infeas[iSequence] = value; // already there
          else
            infeasible_->quickAdd(iSequence, value);
        } else {
          infeasible_->zero(iSequence);
        }
        break;
      case ClpSimplex::atLowerBound:
        iSequence += addSequence;
        if (value < -tolerance) {
#ifdef CLP_PRIMAL_SLACK_MULTIPLIER
          value *= value * slack_multiplier;
#else
          value *= value;
#endif
          // store square in list
          if (infeas[iSequence])
            infeas[iSequence] = value; // already there
          else
            infeasible_->quickAdd(iSequence, value);
        } else {
          infeasible_->zero(iSequence);
        }
      }
    }
  }
  // we can zero out as will have to get pivot row
  // ***** move
  updates->setNumElements(0);
  spareColumn1->setNumElements(0);
  if (pivotRow >= 0) {
    // make sure infeasibility on incoming is 0.0
    int sequenceIn = model_->sequenceIn();
    infeasible_->zero(sequenceIn);
  }
  // for weights update we use pivotSequence
  pivotRow = pivotSequence_;
  // unset in case sub flip
  pivotSequence_ = -1;
  if (pivotRow >= 0) {
    // make sure infeasibility on incoming is 0.0
    const int *pivotVariable = model_->pivotVariable();
    int sequenceIn = pivotVariable[pivotRow];
    assert(sequenceIn == model_->sequenceIn());
    infeasible_->zero(sequenceIn);
    // and we can see if reference
    double referenceIn;
    if (mode_ != 1) {
      if (reference(sequenceIn))
        referenceIn = 1.0;
      else
        referenceIn = 0.0;
    } else {
      referenceIn = -1.0;
    }
    // save outgoing weight round update
    double outgoingWeight = 0.0;
    int sequenceOut = model_->sequenceOut();
    if (sequenceOut >= 0)
      outgoingWeight = weights_[sequenceOut];
    // update weights
    updates->setNumElements(0);
    spareColumn1->setNumElements(0);
    // might as well set dj to 1
    dj = -1.0;
    updates->createPacked(1, &pivotRow, &dj);
    model_->factorization()->updateColumnTranspose(spareRow2, updates);
    bool needSubset = (mode_ < 4 || numberSwitched_ > 1 || mode_ >= 10);

    double *weight;
    double *other = alternateWeights_->denseVector();
    int numberColumns = model_->numberColumns();
    // rows
    number = updates->getNumElements();
    index = updates->getIndices();
    updateBy = updates->denseVector();
    weight = weights_ + numberColumns;
    if (needSubset) {
#if ALT_UPDATE_WEIGHTS != 2
      model_->factorization()->updateColumnTranspose(spareRow2,
        alternateWeights_);
#elif ALT_UPDATE_WEIGHTS == 1
      if (altVector[1]) {
        int numberRows = model_->numberRows();
        double *work1 = altVector[1]->denseVector();
        double *worka = alternateWeights_->denseVector();
        int iRow = -1;
        double diff = 1.0e-8;
        for (int i = 0; i < numberRows; i++) {
          double dd = std::max(fabs(work1[i]), fabs(worka[i]));
          double d = fabs(work1[i] - worka[i]);
          if (dd > 1.0e-6 && d > diff * dd) {
            diff = d / dd;
            iRow = i;
          }
        }
        if (iRow >= 0)
          printf("largest2 difference of %g (%g,%g) on row %d\n",
            diff, work1[iRow], worka[iRow], iRow);
      }
#endif
      // do alternateWeights_ here so can scale
      for (j = 0; j < number; j++) {
        int iSequence = index[j];
        assert(iSequence >= 0 && iSequence < model_->numberRows());
        double thisWeight = weight[iSequence];
        // row has -1
        double pivot = -updateBy[j];
        double modification = other[iSequence];
        double pivotSquared = pivot * pivot;

        thisWeight += pivotSquared * devex_ + pivot * modification;
        if (thisWeight < TRY_NORM) {
          if (mode_ == 1) {
            // steepest
            thisWeight = std::max(TRY_NORM, ADD_ONE + pivotSquared);
          } else {
            // exact
            thisWeight = referenceIn * pivotSquared;
            if (reference(iSequence + numberColumns))
              thisWeight += 1.0;
            thisWeight = std::max(thisWeight, TRY_NORM);
          }
        }
        weight[iSequence] = thisWeight;
      }
      transposeTimes2(updates, spareColumn1, alternateWeights_, spareColumn2, spareRow2, 0.0);
    } else {
      // put row of tableau in rowArray and columnArray
      model_->clpMatrix()->transposeTimes(model_, -1.0,
        updates, spareColumn2, spareColumn1);
    }

    if (needSubset) {
      CoinZeroN(updateBy, number);
    } else if (mode_ == 4) {
      // Devex
      assert(devex_ > 0.0);
      for (j = 0; j < number; j++) {
        int iSequence = index[j];
        double thisWeight = weight[iSequence];
        // row has -1
        double pivot = -updateBy[j];
        updateBy[j] = 0.0;
        double value = pivot * pivot * devex_;
        if (reference(iSequence + numberColumns))
          value += 1.0;
        weight[iSequence] = std::max(0.99 * thisWeight, value);
      }
    }

    // columns
    weight = weights_;

    number = spareColumn1->getNumElements();
    index = spareColumn1->getIndices();
    updateBy = spareColumn1->denseVector();
    if (needSubset) {
      // Exact - already done
    } else if (mode_ == 4) {
      // Devex
      for (j = 0; j < number; j++) {
        int iSequence = index[j];
        double thisWeight = weight[iSequence];
        // row has -1
        double pivot = updateBy[j];
        updateBy[j] = 0.0;
        double value = pivot * pivot * devex_;
        if (reference(iSequence))
          value += 1.0;
        weight[iSequence] = std::max(0.99 * thisWeight, value);
      }
    }
    // restore outgoing weight
    if (sequenceOut >= 0)
      weights_[sequenceOut] = outgoingWeight;
    alternateWeights_->clear();
    spareColumn2->setNumElements(0);
    //#define SOME_DEBUG_1
#ifdef SOME_DEBUG_1
    // check for accuracy
    int iCheck = 892;
    //printf("weight for iCheck is %g\n",weights_[iCheck]);
    int numberRows = model_->numberRows();
    //int numberColumns = model_->numberColumns();
    for (iCheck = 0; iCheck < numberRows + numberColumns; iCheck++) {
      if (model_->getStatus(iCheck) != ClpSimplex::basic && model_->getStatus(iCheck) != ClpSimplex::isFixed)
        checkAccuracy(iCheck, 1.0e-1, updates, spareRow2);
    }
#endif
  }
  updates->setNumElements(0);
  spareColumn1->setNumElements(0);
}
// Updates two arrays for steepest
int ClpPrimalColumnSteepest::transposeTimes2(const CoinIndexedVector *pi1, CoinIndexedVector *dj1,
  const CoinIndexedVector *pi2, CoinIndexedVector *dj2,
  CoinIndexedVector *spare,
  double scaleFactor)
{
  // see if reference
  int sequenceIn = model_->sequenceIn();
  double referenceIn;
  if (mode_ != 1) {
    if (reference(sequenceIn))
      referenceIn = 1.0;
    else
      referenceIn = 0.0;
  } else {
    referenceIn = -1.0;
  }
  int returnCode = 0;
  if (model_->clpMatrix()->canCombine(model_, pi1)) {
    double *infeas = scaleFactor ? infeasible_->denseVector() : NULL;
    // put row of tableau in rowArray and columnArray
    returnCode = model_->clpMatrix()->transposeTimes2(model_, pi1,
      dj1, pi2, spare,
      infeas,
      model_->djRegion(1),
      referenceIn, devex_,
      reference_,
      weights_, scaleFactor);
    if (model_->spareIntArray_[3] > -2)
      returnCode = 2;
  } else {
    // put row of tableau in rowArray and columnArray
    model_->clpMatrix()->transposeTimes(model_, -1.0,
      pi1, dj2, dj1);
    // get subset which have nonzero tableau elements
    model_->clpMatrix()->subsetTransposeTimes(model_, pi2, dj1, dj2);
    bool killDjs = (scaleFactor == 0.0);
    if (!scaleFactor)
      scaleFactor = 1.0;
    // columns
    double *weight = weights_;

    int number = dj1->getNumElements();
    const int *index = dj1->getIndices();
    double *updateBy = dj1->denseVector();
    double *updateBy2 = dj2->denseVector();

    for (int j = 0; j < number; j++) {
      double thisWeight;
      double pivot;
      double pivotSquared;
      int iSequence = index[j];
      double value2 = updateBy[j];
      if (killDjs)
        updateBy[j] = 0.0;
      double modification = updateBy2[j];
      updateBy2[j] = 0.0;
      ClpSimplex::Status status = model_->getStatus(iSequence);

      if (status != ClpSimplex::basic && status != ClpSimplex::isFixed) {
        thisWeight = weight[iSequence];
        pivot = value2 * scaleFactor;
        pivotSquared = pivot * pivot;

        thisWeight += pivotSquared * devex_ + pivot * modification;
        if (thisWeight < TRY_NORM) {
          if (referenceIn < 0.0) {
            // steepest
            thisWeight = std::max(TRY_NORM, ADD_ONE + pivotSquared);
          } else {
            // exact
            thisWeight = referenceIn * pivotSquared;
            if (reference(iSequence))
              thisWeight += 1.0;
            thisWeight = std::max(thisWeight, TRY_NORM);
          }
        }
        weight[iSequence] = thisWeight;
      }
    }
  }
  dj2->setNumElements(0);
  return returnCode;
}
// Update weights for Devex
void ClpPrimalColumnSteepest::justDevex(CoinIndexedVector *updates,
  CoinIndexedVector *spareRow2,
  CoinIndexedVector *spareColumn1,
  CoinIndexedVector *spareColumn2)
{
  int j;
  int number = 0;
  int *index;
  double *updateBy;
  // dj could be very small (or even zero - take care)
  double dj = model_->dualIn();
  double tolerance = model_->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = std::min(1.0e-2, model_->largestDualError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance + error;
  int pivotRow = model_->pivotRow();
  // for weights update we use pivotSequence
  pivotRow = pivotSequence_;
  assert(pivotRow >= 0);
  // make sure infeasibility on incoming is 0.0
  const int *pivotVariable = model_->pivotVariable();
  int sequenceIn = pivotVariable[pivotRow];
  infeasible_->zero(sequenceIn);
  // and we can see if reference
  //double referenceIn = 0.0;
  //if (mode_ != 1 && reference(sequenceIn))
  //   referenceIn = 1.0;
  // save outgoing weight round update
  double outgoingWeight = 0.0;
  int sequenceOut = model_->sequenceOut();
  if (sequenceOut >= 0)
    outgoingWeight = weights_[sequenceOut];
  assert(!updates->getNumElements());
  assert(!spareColumn1->getNumElements());
  // unset in case sub flip
  pivotSequence_ = -1;
  // might as well set dj to 1
  dj = -1.0;
  updates->createPacked(1, &pivotRow, &dj);
  model_->factorization()->updateColumnTranspose(spareRow2, updates);
  // put row of tableau in rowArray and columnArray
  model_->clpMatrix()->transposeTimes(model_, -1.0,
    updates, spareColumn2, spareColumn1);
  double *weight;
  int numberColumns = model_->numberColumns();
  // rows
  number = updates->getNumElements();
  index = updates->getIndices();
  updateBy = updates->denseVector();
  weight = weights_ + numberColumns;

  // Devex
  assert(devex_ > 0.0);
  for (j = 0; j < number; j++) {
    int iSequence = index[j];
    double thisWeight = weight[iSequence];
    // row has -1
    double pivot = -updateBy[j];
    updateBy[j] = 0.0;
    double value = pivot * pivot * devex_;
    if (reference(iSequence + numberColumns))
      value += 1.0;
    weight[iSequence] = std::max(0.99 * thisWeight, value);
  }

  // columns
  weight = weights_;

  number = spareColumn1->getNumElements();
  index = spareColumn1->getIndices();
  updateBy = spareColumn1->denseVector();
  // Devex
  for (j = 0; j < number; j++) {
    int iSequence = index[j];
    double thisWeight = weight[iSequence];
    // row has -1
    double pivot = updateBy[j];
    updateBy[j] = 0.0;
    double value = pivot * pivot * devex_;
    if (reference(iSequence))
      value += 1.0;
    weight[iSequence] = std::max(0.99 * thisWeight, value);
  }
  // restore outgoing weight
  if (sequenceOut >= 0)
    weights_[sequenceOut] = outgoingWeight;
  spareColumn2->setNumElements(0);
  //#define SOME_DEBUG_1
#ifdef SOME_DEBUG_1
  // check for accuracy
  int iCheck = 892;
  //printf("weight for iCheck is %g\n",weights_[iCheck]);
  int numberRows = model_->numberRows();
  //int numberColumns = model_->numberColumns();
  for (iCheck = 0; iCheck < numberRows + numberColumns; iCheck++) {
    if (model_->getStatus(iCheck) != ClpSimplex::basic && model_->getStatus(iCheck) != ClpSimplex::isFixed)
      checkAccuracy(iCheck, 1.0e-1, updates, spareRow2);
  }
#endif
  updates->setNumElements(0);
  spareColumn1->setNumElements(0);
}
// Update weights for Steepest
void ClpPrimalColumnSteepest::justSteepest(CoinIndexedVector *updates,
  CoinIndexedVector *spareRow2,
  CoinIndexedVector *spareColumn1,
  CoinIndexedVector *spareColumn2)
{
  int j;
  int number = 0;
  int *index;
  double *updateBy;
  // dj could be very small (or even zero - take care)
  double dj = model_->dualIn();
  double tolerance = model_->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = std::min(1.0e-2, model_->largestDualError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance + error;
  int pivotRow = model_->pivotRow();
  // for weights update we use pivotSequence
  pivotRow = pivotSequence_;
  // unset in case sub flip
  pivotSequence_ = -1;
  assert(pivotRow >= 0);
  // make sure infeasibility on incoming is 0.0
  const int *pivotVariable = model_->pivotVariable();
  int sequenceIn = pivotVariable[pivotRow];
  infeasible_->zero(sequenceIn);
  // and we can see if reference
  double referenceIn = 0.0;
  if (mode_ != 1 && reference(sequenceIn))
    referenceIn = 1.0;
  // save outgoing weight round update
  double outgoingWeight = 0.0;
  int sequenceOut = model_->sequenceOut();
  if (sequenceOut >= 0)
    outgoingWeight = weights_[sequenceOut];
  assert(!updates->getNumElements());
  assert(!spareColumn1->getNumElements());
  // update weights
  //updates->setNumElements(0);
  //spareColumn1->setNumElements(0);
  // might as well set dj to 1
  dj = -1.0;
  updates->createPacked(1, &pivotRow, &dj);
  model_->factorization()->updateColumnTranspose(spareRow2, updates);
  // put row of tableau in rowArray and columnArray
  model_->clpMatrix()->transposeTimes(model_, -1.0,
    updates, spareColumn2, spareColumn1);
  double *weight;
  double *other = alternateWeights_->denseVector();
  int numberColumns = model_->numberColumns();
  // rows
  number = updates->getNumElements();
  index = updates->getIndices();
  updateBy = updates->denseVector();
  weight = weights_ + numberColumns;

  // Exact
  // now update weight update array
  //alternateWeights_->scanAndPack();
#if ALT_UPDATE_WEIGHTS != 2
  model_->factorization()->updateColumnTranspose(spareRow2,
    alternateWeights_);
#elif ALT_UPDATE_WEIGHTS == 1
  if (altVector[1]) {
    int numberRows = model_->numberRows();
    double *work1 = altVector[1]->denseVector();
    double *worka = alternateWeights_->denseVector();
    int iRow = -1;
    double diff = 1.0e-8;
    for (int i = 0; i < numberRows; i++) {
      double dd = std::max(fabs(work1[i]), fabs(worka[i]));
      double d = fabs(work1[i] - worka[i]);
      if (dd > 1.0e-6 && d > diff * dd) {
        diff = d / dd;
        iRow = i;
      }
    }
    if (iRow >= 0)
      printf("largest3 difference of %g (%g,%g) on row %d\n",
        diff, work1[iRow], worka[iRow], iRow);
  }
#endif
  // get subset which have nonzero tableau elements
  model_->clpMatrix()->subsetTransposeTimes(model_, alternateWeights_,
    spareColumn1,
    spareColumn2);
  for (j = 0; j < number; j++) {
    int iSequence = index[j];
    double thisWeight = weight[iSequence];
    // row has -1
    double pivot = -updateBy[j];
    updateBy[j] = 0.0;
    double modification = other[iSequence];
    double pivotSquared = pivot * pivot;

    thisWeight += pivotSquared * devex_ + pivot * modification;
    if (thisWeight < TRY_NORM) {
      if (mode_ == 1) {
        // steepest
        thisWeight = std::max(TRY_NORM, ADD_ONE + pivotSquared);
      } else {
        // exact
        thisWeight = referenceIn * pivotSquared;
        if (reference(iSequence + numberColumns))
          thisWeight += 1.0;
        thisWeight = std::max(thisWeight, TRY_NORM);
      }
    }
    weight[iSequence] = thisWeight;
  }

  // columns
  weight = weights_;

  number = spareColumn1->getNumElements();
  index = spareColumn1->getIndices();
  updateBy = spareColumn1->denseVector();
  // Exact
  double *updateBy2 = spareColumn2->denseVector();
  for (j = 0; j < number; j++) {
    int iSequence = index[j];
    double thisWeight = weight[iSequence];
    double pivot = updateBy[j];
    updateBy[j] = 0.0;
    double modification = updateBy2[j];
    updateBy2[j] = 0.0;
    double pivotSquared = pivot * pivot;

    thisWeight += pivotSquared * devex_ + pivot * modification;
    if (thisWeight < TRY_NORM) {
      if (mode_ == 1) {
        // steepest
        thisWeight = std::max(TRY_NORM, ADD_ONE + pivotSquared);
      } else {
        // exact
        thisWeight = referenceIn * pivotSquared;
        if (reference(iSequence))
          thisWeight += 1.0;
        thisWeight = std::max(thisWeight, TRY_NORM);
      }
    }
    weight[iSequence] = thisWeight;
  }
  // restore outgoing weight
  if (sequenceOut >= 0)
    weights_[sequenceOut] = outgoingWeight;
  alternateWeights_->clear();
  spareColumn2->setNumElements(0);
  //#define SOME_DEBUG_1
#ifdef SOME_DEBUG_1
  // check for accuracy
  int iCheck = 892;
  //printf("weight for iCheck is %g\n",weights_[iCheck]);
  int numberRows = model_->numberRows();
  //int numberColumns = model_->numberColumns();
  for (iCheck = 0; iCheck < numberRows + numberColumns; iCheck++) {
    if (model_->getStatus(iCheck) != ClpSimplex::basic && model_->getStatus(iCheck) != ClpSimplex::isFixed)
      checkAccuracy(iCheck, 1.0e-1, updates, spareRow2);
  }
#endif
  updates->setNumElements(0);
  spareColumn1->setNumElements(0);
}
// Returns pivot column, -1 if none
int ClpPrimalColumnSteepest::pivotColumnOldMethod(CoinIndexedVector *updates,
  CoinIndexedVector *,
  CoinIndexedVector *spareRow2,
  CoinIndexedVector *spareColumn1,
  CoinIndexedVector *spareColumn2)
{
  assert(model_);
  int iSection, j;
  int number = 0;
  int *index;
  double *updateBy;
  double *reducedCost;
  // dj could be very small (or even zero - take care)
  double dj = model_->dualIn();
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
        5 - dantzig
     */
  if (updates->getNumElements()) {
    // would have to have two goes for devex, three for steepest
    anyUpdates = 2;
    // add in pivot contribution
    if (pivotRow >= 0)
      updates->add(pivotRow, -dj);
  } else if (pivotRow >= 0) {
    if (fabs(dj) > 1.0e-15) {
      // some dj
      updates->insert(pivotRow, -dj);
      if (fabs(dj) > 1.0e-6) {
        // reasonable size
        anyUpdates = 1;
      } else {
        // too small
        anyUpdates = 2;
      }
    } else {
      // zero dj
      anyUpdates = -1;
    }
  } else if (pivotSequence_ >= 0) {
    // just after re-factorization
    anyUpdates = -1;
  } else {
    // sub flip - nothing to do
    anyUpdates = 0;
  }

  if (anyUpdates > 0) {
    model_->factorization()->updateColumnTranspose(spareRow2, updates);

    // put row of tableau in rowArray and columnArray
    model_->clpMatrix()->transposeTimes(model_, -1.0,
      updates, spareColumn2, spareColumn1);
    // normal
    for (iSection = 0; iSection < 2; iSection++) {

      reducedCost = model_->djRegion(iSection);
      int addSequence;

      if (!iSection) {
        number = updates->getNumElements();
        index = updates->getIndices();
        updateBy = updates->denseVector();
        addSequence = model_->numberColumns();
      } else {
        number = spareColumn1->getNumElements();
        index = spareColumn1->getIndices();
        updateBy = spareColumn1->denseVector();
        addSequence = 0;
      }
      if (!model_->nonLinearCost()->lookBothWays()) {

        for (j = 0; j < number; j++) {
          int iSequence = index[j];
          double value = reducedCost[iSequence];
          value -= updateBy[iSequence];
          reducedCost[iSequence] = value;
          ClpSimplex::Status status = model_->getStatus(iSequence + addSequence);

          switch (status) {

          case ClpSimplex::basic:
            infeasible_->zero(iSequence + addSequence);
          case ClpSimplex::isFixed:
            break;
          case ClpSimplex::isFree:
          case ClpSimplex::superBasic:
            if (fabs(value) > FREE_ACCEPT * tolerance) {
              // we are going to bias towards free (but only if reasonable)
              value *= FREE_BIAS;
              // store square in list
              if (infeas[iSequence + addSequence])
                infeas[iSequence + addSequence] = value * value; // already there
              else
                infeasible_->quickAdd(iSequence + addSequence, value * value);
            } else {
              infeasible_->zero(iSequence + addSequence);
            }
            break;
          case ClpSimplex::atUpperBound:
            if (value > tolerance) {
              // store square in list
              if (infeas[iSequence + addSequence])
                infeas[iSequence + addSequence] = value * value; // already there
              else
                infeasible_->quickAdd(iSequence + addSequence, value * value);
            } else {
              infeasible_->zero(iSequence + addSequence);
            }
            break;
          case ClpSimplex::atLowerBound:
            if (value < -tolerance) {
              // store square in list
              if (infeas[iSequence + addSequence])
                infeas[iSequence + addSequence] = value * value; // already there
              else
                infeasible_->quickAdd(iSequence + addSequence, value * value);
            } else {
              infeasible_->zero(iSequence + addSequence);
            }
          }
        }
      } else {
        ClpNonLinearCost *nonLinear = model_->nonLinearCost();
        // We can go up OR down
        for (j = 0; j < number; j++) {
          int iSequence = index[j];
          double value = reducedCost[iSequence];
          value -= updateBy[iSequence];
          reducedCost[iSequence] = value;
          ClpSimplex::Status status = model_->getStatus(iSequence + addSequence);

          switch (status) {

          case ClpSimplex::basic:
            infeasible_->zero(iSequence + addSequence);
          case ClpSimplex::isFixed:
            break;
          case ClpSimplex::isFree:
          case ClpSimplex::superBasic:
            if (fabs(value) > FREE_ACCEPT * tolerance) {
              // we are going to bias towards free (but only if reasonable)
              value *= FREE_BIAS;
              // store square in list
              if (infeas[iSequence + addSequence])
                infeas[iSequence + addSequence] = value * value; // already there
              else
                infeasible_->quickAdd(iSequence + addSequence, value * value);
            } else {
              infeasible_->zero(iSequence + addSequence);
            }
            break;
          case ClpSimplex::atUpperBound:
            if (value > tolerance) {
              // store square in list
              if (infeas[iSequence + addSequence])
                infeas[iSequence + addSequence] = value * value; // already there
              else
                infeasible_->quickAdd(iSequence + addSequence, value * value);
            } else {
              // look other way - change up should be negative
              value -= nonLinear->changeUpInCost(iSequence + addSequence);
              if (value > -tolerance) {
                infeasible_->zero(iSequence + addSequence);
              } else {
                // store square in list
                if (infeas[iSequence + addSequence])
                  infeas[iSequence + addSequence] = value * value; // already there
                else
                  infeasible_->quickAdd(iSequence + addSequence, value * value);
              }
            }
            break;
          case ClpSimplex::atLowerBound:
            if (value < -tolerance) {
              // store square in list
              if (infeas[iSequence + addSequence])
                infeas[iSequence + addSequence] = value * value; // already there
              else
                infeasible_->quickAdd(iSequence + addSequence, value * value);
            } else {
              // look other way - change down should be positive
              value -= nonLinear->changeDownInCost(iSequence + addSequence);
              if (value < tolerance) {
                infeasible_->zero(iSequence + addSequence);
              } else {
                // store square in list
                if (infeas[iSequence + addSequence])
                  infeas[iSequence + addSequence] = value * value; // already there
                else
                  infeasible_->quickAdd(iSequence + addSequence, value * value);
              }
            }
          }
        }
      }
    }
    if (anyUpdates == 2) {
      // we can zero out as will have to get pivot row
      updates->clear();
      spareColumn1->clear();
    }
    if (pivotRow >= 0) {
      // make sure infeasibility on incoming is 0.0
      int sequenceIn = model_->sequenceIn();
      infeasible_->zero(sequenceIn);
    }
  }
  // make sure outgoing from last iteration okay
  int sequenceOut = model_->sequenceOut();
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

  // If in quadratic re-compute all
  if (model_->algorithm() == 2) {
    for (iSection = 0; iSection < 2; iSection++) {

      reducedCost = model_->djRegion(iSection);
      int addSequence;
      int iSequence;

      if (!iSection) {
        number = model_->numberRows();
        addSequence = model_->numberColumns();
      } else {
        number = model_->numberColumns();
        addSequence = 0;
      }

      if (!model_->nonLinearCost()->lookBothWays()) {
        for (iSequence = 0; iSequence < number; iSequence++) {
          double value = reducedCost[iSequence];
          ClpSimplex::Status status = model_->getStatus(iSequence + addSequence);

          switch (status) {

          case ClpSimplex::basic:
            infeasible_->zero(iSequence + addSequence);
          case ClpSimplex::isFixed:
            break;
          case ClpSimplex::isFree:
          case ClpSimplex::superBasic:
            if (fabs(value) > tolerance) {
              // we are going to bias towards free (but only if reasonable)
              value *= FREE_BIAS;
              // store square in list
              if (infeas[iSequence + addSequence])
                infeas[iSequence + addSequence] = value * value; // already there
              else
                infeasible_->quickAdd(iSequence + addSequence, value * value);
            } else {
              infeasible_->zero(iSequence + addSequence);
            }
            break;
          case ClpSimplex::atUpperBound:
            if (value > tolerance) {
              // store square in list
              if (infeas[iSequence + addSequence])
                infeas[iSequence + addSequence] = value * value; // already there
              else
                infeasible_->quickAdd(iSequence + addSequence, value * value);
            } else {
              infeasible_->zero(iSequence + addSequence);
            }
            break;
          case ClpSimplex::atLowerBound:
            if (value < -tolerance) {
              // store square in list
              if (infeas[iSequence + addSequence])
                infeas[iSequence + addSequence] = value * value; // already there
              else
                infeasible_->quickAdd(iSequence + addSequence, value * value);
            } else {
              infeasible_->zero(iSequence + addSequence);
            }
          }
        }
      } else {
        // we can go both ways
        ClpNonLinearCost *nonLinear = model_->nonLinearCost();
        for (iSequence = 0; iSequence < number; iSequence++) {
          double value = reducedCost[iSequence];
          ClpSimplex::Status status = model_->getStatus(iSequence + addSequence);

          switch (status) {

          case ClpSimplex::basic:
            infeasible_->zero(iSequence + addSequence);
          case ClpSimplex::isFixed:
            break;
          case ClpSimplex::isFree:
          case ClpSimplex::superBasic:
            if (fabs(value) > tolerance) {
              // we are going to bias towards free (but only if reasonable)
              value *= FREE_BIAS;
              // store square in list
              if (infeas[iSequence + addSequence])
                infeas[iSequence + addSequence] = value * value; // already there
              else
                infeasible_->quickAdd(iSequence + addSequence, value * value);
            } else {
              infeasible_->zero(iSequence + addSequence);
            }
            break;
          case ClpSimplex::atUpperBound:
            if (value > tolerance) {
              // store square in list
              if (infeas[iSequence + addSequence])
                infeas[iSequence + addSequence] = value * value; // already there
              else
                infeasible_->quickAdd(iSequence + addSequence, value * value);
            } else {
              // look other way - change up should be negative
              value -= nonLinear->changeUpInCost(iSequence + addSequence);
              if (value > -tolerance) {
                infeasible_->zero(iSequence + addSequence);
              } else {
                // store square in list
                if (infeas[iSequence + addSequence])
                  infeas[iSequence + addSequence] = value * value; // already there
                else
                  infeasible_->quickAdd(iSequence + addSequence, value * value);
              }
            }
            break;
          case ClpSimplex::atLowerBound:
            if (value < -tolerance) {
              // store square in list
              if (infeas[iSequence + addSequence])
                infeas[iSequence + addSequence] = value * value; // already there
              else
                infeasible_->quickAdd(iSequence + addSequence, value * value);
            } else {
              // look other way - change down should be positive
              value -= nonLinear->changeDownInCost(iSequence + addSequence);
              if (value < tolerance) {
                infeasible_->zero(iSequence + addSequence);
              } else {
                // store square in list
                if (infeas[iSequence + addSequence])
                  infeas[iSequence + addSequence] = value * value; // already there
                else
                  infeasible_->quickAdd(iSequence + addSequence, value * value);
              }
            }
          }
        }
      }
    }
  }
  // See what sort of pricing
  int numberWanted = 10;
  number = infeasible_->getNumElements();
  int numberColumns = model_->numberColumns();
  if (switchType == 5) {
    // we can zero out
    updates->clear();
    spareColumn1->clear();
    pivotSequence_ = -1;
    pivotRow = -1;
    // See if to switch
    int numberRows = model_->numberRows();
    // ratio is done on number of columns here
    //double ratio = static_cast<double> sizeFactorization_/static_cast<double> numberColumns;
    double ratio = static_cast< double >(sizeFactorization_) / static_cast< double >(numberRows);
    //double ratio = static_cast<double> sizeFactorization_/static_cast<double> model_->clpMatrix()->getNumElements();
    if (ratio < 0.1) {
      numberWanted = std::max(100, number / 200);
    } else if (ratio < 0.3) {
      numberWanted = std::max(500, number / 40);
    } else if (ratio < 0.5 || mode_ == 5) {
      numberWanted = std::max(2000, number / 10);
      numberWanted = std::max(numberWanted, numberColumns / 30);
    } else if (mode_ != 5) {
      switchType = 4;
      // initialize
      numberSwitched_++;
      // Make sure will re-do
      delete[] weights_;
      weights_ = NULL;
      saveWeights(model_, 4);
      COIN_DETAIL_PRINT(printf("switching to devex %d nel ratio %g\n", sizeFactorization_, ratio));
      updates->clear();
    }
    if (model_->numberIterations() % 1000 == 0)
      COIN_DETAIL_PRINT(printf("numels %d ratio %g wanted %d\n", sizeFactorization_, ratio, numberWanted));
  }
  if (switchType == 4) {
    // Still in devex mode
    int numberRows = model_->numberRows();
    // ratio is done on number of rows here
    double ratio = static_cast< double >(sizeFactorization_) / static_cast< double >(numberRows);
    // Go to steepest if lot of iterations?
    if (ratio < 1.0) {
      numberWanted = std::max(2000, number / 20);
    } else if (ratio < 5.0) {
      numberWanted = std::max(2000, number / 10);
      numberWanted = std::max(numberWanted, numberColumns / 20);
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
      COIN_DETAIL_PRINT(printf("numels %d ratio %g wanted %d\n", sizeFactorization_, ratio, numberWanted));
  }
  if (switchType < 4) {
    if (switchType < 2) {
      numberWanted = number + 1;
    } else if (switchType == 2) {
      numberWanted = std::max(2000, number / 8);
    } else {
      double ratio = static_cast< double >(sizeFactorization_) / static_cast< double >(model_->numberRows());
      if (ratio < 1.0) {
        numberWanted = std::max(2000, number / 20);
      } else if (ratio < 5.0) {
        numberWanted = std::max(2000, number / 10);
        numberWanted = std::max(numberWanted, numberColumns / 20);
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
  // for weights update we use pivotSequence
  pivotRow = pivotSequence_;
  // unset in case sub flip
  pivotSequence_ = -1;
  if (pivotRow >= 0) {
    // make sure infeasibility on incoming is 0.0
    const int *pivotVariable = model_->pivotVariable();
    int sequenceIn = pivotVariable[pivotRow];
    infeasible_->zero(sequenceIn);
    // and we can see if reference
    double referenceIn = 0.0;
    if (switchType != 1 && reference(sequenceIn))
      referenceIn = 1.0;
    // save outgoing weight round update
    double outgoingWeight = 0.0;
    if (sequenceOut >= 0)
      outgoingWeight = weights_[sequenceOut];
    // update weights
    if (anyUpdates != 1) {
      updates->setNumElements(0);
      spareColumn1->setNumElements(0);
      // might as well set dj to 1
      dj = 1.0;
      updates->insert(pivotRow, -dj);
      model_->factorization()->updateColumnTranspose(spareRow2, updates);
      // put row of tableau in rowArray and columnArray
      model_->clpMatrix()->transposeTimes(model_, -1.0,
        updates, spareColumn2, spareColumn1);
    }
    double *weight;
    double *other = alternateWeights_->denseVector();
    int numberColumns = model_->numberColumns();
    double scaleFactor = -1.0 / dj; // as formula is with 1.0
    // rows
    number = updates->getNumElements();
    index = updates->getIndices();
    updateBy = updates->denseVector();
    weight = weights_ + numberColumns;

    if (switchType < 4) {
      // Exact
      // now update weight update array
      model_->factorization()->updateColumnTranspose(spareRow2,
        alternateWeights_);
#ifdef ALT_UPDATE_WEIGHTS
      abort();
#endif
      for (j = 0; j < number; j++) {
        int iSequence = index[j];
        double thisWeight = weight[iSequence];
        // row has -1
        double pivot = updateBy[iSequence] * scaleFactor;
        updateBy[iSequence] = 0.0;
        double modification = other[iSequence];
        double pivotSquared = pivot * pivot;

        thisWeight += pivotSquared * devex_ + pivot * modification;
        if (thisWeight < TRY_NORM) {
          if (switchType == 1) {
            // steepest
            thisWeight = std::max(TRY_NORM, ADD_ONE + pivotSquared);
          } else {
            // exact
            thisWeight = referenceIn * pivotSquared;
            if (reference(iSequence + numberColumns))
              thisWeight += 1.0;
            thisWeight = std::max(thisWeight, TRY_NORM);
          }
        }
        weight[iSequence] = thisWeight;
      }
    } else if (switchType == 4) {
      // Devex
      assert(devex_ > 0.0);
      for (j = 0; j < number; j++) {
        int iSequence = index[j];
        double thisWeight = weight[iSequence];
        // row has -1
        double pivot = updateBy[iSequence] * scaleFactor;
        updateBy[iSequence] = 0.0;
        double value = pivot * pivot * devex_;
        if (reference(iSequence + numberColumns))
          value += 1.0;
        weight[iSequence] = std::max(0.99 * thisWeight, value);
      }
    }

    // columns
    weight = weights_;

    scaleFactor = -scaleFactor;

    number = spareColumn1->getNumElements();
    index = spareColumn1->getIndices();
    updateBy = spareColumn1->denseVector();
    if (switchType < 4) {
      // Exact
      // get subset which have nonzero tableau elements
      model_->clpMatrix()->subsetTransposeTimes(model_, alternateWeights_,
        spareColumn1,
        spareColumn2);
      double *updateBy2 = spareColumn2->denseVector();
      for (j = 0; j < number; j++) {
        int iSequence = index[j];
        double thisWeight = weight[iSequence];
        double pivot = updateBy[iSequence] * scaleFactor;
        updateBy[iSequence] = 0.0;
        double modification = updateBy2[j];
        updateBy2[j] = 0.0;
        double pivotSquared = pivot * pivot;

        thisWeight += pivotSquared * devex_ + pivot * modification;
        if (thisWeight < TRY_NORM) {
          if (switchType == 1) {
            // steepest
            thisWeight = std::max(TRY_NORM, ADD_ONE + pivotSquared);
          } else {
            // exact
            thisWeight = referenceIn * pivotSquared;
            if (reference(iSequence))
              thisWeight += 1.0;
            thisWeight = std::max(thisWeight, TRY_NORM);
          }
        }
        weight[iSequence] = thisWeight;
      }
    } else if (switchType == 4) {
      // Devex
      for (j = 0; j < number; j++) {
        int iSequence = index[j];
        double thisWeight = weight[iSequence];
        // row has -1
        double pivot = updateBy[iSequence] * scaleFactor;
        updateBy[iSequence] = 0.0;
        double value = pivot * pivot * devex_;
        if (reference(iSequence))
          value += 1.0;
        weight[iSequence] = std::max(0.99 * thisWeight, value);
      }
    }
    // restore outgoing weight
    if (sequenceOut >= 0)
      weights_[sequenceOut] = outgoingWeight;
    alternateWeights_->clear();
    spareColumn2->setNumElements(0);
    //#define SOME_DEBUG_1
#ifdef SOME_DEBUG_1
    // check for accuracy
    int iCheck = 892;
    //printf("weight for iCheck is %g\n",weights_[iCheck]);
    int numberRows = model_->numberRows();
    //int numberColumns = model_->numberColumns();
    for (iCheck = 0; iCheck < numberRows + numberColumns; iCheck++) {
      if (model_->getStatus(iCheck) != ClpSimplex::basic && model_->getStatus(iCheck) != ClpSimplex::isFixed)
        checkAccuracy(iCheck, 1.0e-1, updates, spareRow2);
    }
#endif
    updates->setNumElements(0);
    spareColumn1->setNumElements(0);
  }

  // update of duals finished - now do pricing

  double bestDj = 1.0e-30;
  int bestSequence = -1;

  int i, iSequence;
  index = infeasible_->getIndices();
  number = infeasible_->getNumElements();
  if (model_->numberIterations() < model_->lastBadIteration() + 200) {
    // we can't really trust infeasibilities if there is dual error
    double checkTolerance = 1.0e-8;
    if (!model_->factorization()->pivots())
      checkTolerance = 1.0e-6;
    if (model_->largestDualError() > checkTolerance)
      tolerance *= model_->largestDualError() / checkTolerance;
    // But cap
    tolerance = std::min(1000.0, tolerance);
  }
#ifdef CLP_DEBUG
  if (model_->numberDualInfeasibilities() == 1)
    printf("** %g %g %g %x %x %d\n", tolerance, model_->dualTolerance(),
      model_->largestDualError(), model_, model_->messageHandler(),
      number);
#endif
  // stop last one coming immediately
  double saveOutInfeasibility = 0.0;
  if (sequenceOut >= 0) {
    saveOutInfeasibility = infeas[sequenceOut];
    infeas[sequenceOut] = 0.0;
  }
  tolerance *= tolerance; // as we are using squares

  int iPass;
  // Setup two passes
  int start[4];
  start[1] = number;
  start[2] = 0;
  double dstart = static_cast< double >(number) * model_->randomNumberGenerator()->randomDouble();
  start[0] = static_cast< int >(dstart);
  start[3] = start[0];
  //double largestWeight=0.0;
  //double smallestWeight=1.0e100;
  for (iPass = 0; iPass < 2; iPass++) {
    int end = start[2 * iPass + 1];
    if (switchType < 5) {
      for (i = start[2 * iPass]; i < end; i++) {
        iSequence = index[i];
        double value = infeas[iSequence];
        if (value > tolerance) {
          double weight = weights_[iSequence];
          //weight=1.0;
          if (value > bestDj * weight) {
            // check flagged variable and correct dj
            if (!model_->flagged(iSequence)) {
              bestDj = value / weight;
              bestSequence = iSequence;
            } else {
              // just to make sure we don't exit before got something
              numberWanted++;
            }
          }
        }
        numberWanted--;
        if (!numberWanted)
          break;
      }
    } else {
      // Dantzig
      for (i = start[2 * iPass]; i < end; i++) {
        iSequence = index[i];
        double value = infeas[iSequence];
        if (value > tolerance) {
          if (value > bestDj) {
            // check flagged variable and correct dj
            if (!model_->flagged(iSequence)) {
              bestDj = value;
              bestSequence = iSequence;
            } else {
              // just to make sure we don't exit before got something
              numberWanted++;
            }
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
  if (sequenceOut >= 0) {
    infeas[sequenceOut] = saveOutInfeasibility;
  }
  /*if (model_->numberIterations()%100==0)
       printf("%d best %g\n",bestSequence,bestDj);*/
  reducedCost = model_->djRegion();
  model_->clpMatrix()->setSavedBestSequence(bestSequence);
  if (bestSequence >= 0)
    model_->clpMatrix()->setSavedBestDj(reducedCost[bestSequence]);

#ifdef CLP_DEBUG
  if (bestSequence >= 0) {
    if (model_->getStatus(bestSequence) == ClpSimplex::atLowerBound)
      assert(model_->reducedCost(bestSequence) < 0.0);
    if (model_->getStatus(bestSequence) == ClpSimplex::atUpperBound)
      assert(model_->reducedCost(bestSequence) > 0.0);
  }
#endif
  return bestSequence;
}
// Called when maximum pivots changes
void ClpPrimalColumnSteepest::maximumPivotsChanged()
{
  if (alternateWeights_ && alternateWeights_->capacity() != model_->numberRows() + model_->factorization()->maximumPivots()) {
    delete alternateWeights_;
    alternateWeights_ = new CoinIndexedVector();
    // enough space so can use it for factorization
    alternateWeights_->reserve(model_->numberRows() + model_->factorization()->maximumPivots());
  }
}
void ClpPrimalColumnSteepest::redoInfeasibilities()
{
  double *COIN_RESTRICT infeas = infeasible_->denseVector();
  int *COIN_RESTRICT index = infeasible_->getIndices();
  double tolerance = model_->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = std::min(1.0e-2, model_->largestDualError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance + error;
  // reverse sign so test is cleaner
  tolerance = -tolerance;
  int number = model_->numberRows() + model_->numberColumns();
  int numberNonZero = 0;
  const double *COIN_RESTRICT reducedCost = model_->djRegion();
  const unsigned char *COIN_RESTRICT status = model_->statusArray();
  for (int iSequence = 0; iSequence < number; iSequence++) {
    unsigned char thisStatus = status[iSequence] & 7;
    double value = reducedCost[iSequence];
    infeas[iSequence] = 0.0;
    if (thisStatus == 3) {
    } else if ((thisStatus & 1) != 0) {
      // basic or fixed
      value = 0.0;
    } else if (thisStatus == 2) {
      value = -value;
    } else {
      // free or superbasic
      if (fabs(value) > FREE_ACCEPT * -tolerance) {
        // we are going to bias towards free (but only if reasonable)
        value = -fabs(value) * FREE_BIAS;
      } else {
        value = 0.0;
      }
    }
    if (value < tolerance) {
      // store square in list
      infeas[iSequence] = value * value;
      index[numberNonZero++] = iSequence;
    }
  }
  infeasible_->setNumElements(numberNonZero);
  infeasibilitiesState_ = 0;
}
/*
   1) before factorization
   2) after factorization
   3) just redo infeasibilities
   4) restore weights
   5) at end of values pass (so need initialization)
*/
void ClpPrimalColumnSteepest::saveWeights(ClpSimplex *model, int mode)
{
  model_ = model;
  if (mode == 6) {
    // If incoming weight is 1.0 then return else as 5
    int sequenceIn = model_->sequenceIn();
    assert(sequenceIn >= 0 && sequenceIn < model_->numberRows() + model_->numberColumns());
    // possible weights_ was never set up // assert(weights_);
    if (weights_ && weights_[sequenceIn] == ((mode_ != 1) ? 1.0 : 1.0 + ADD_ONE))
      return;
    else
      mode = 5;
  }
  if (mode_ == 4 || mode_ == 5) {
    if (mode == 1 && !weights_)
      numberSwitched_ = 0; // Reset
  }
  // alternateWeights_ is defined as indexed but is treated oddly
  // at times
  int numberRows = model_->numberRows();
  int numberColumns = model_->numberColumns();
  const int *pivotVariable = model_->pivotVariable();
  bool doInfeasibilities = true;
  if (mode == 1) {
    if (!model_->numberIterations())
      pivotSequence_ = -1;
    if (weights_) {
      // Check if size has changed
      if (infeasible_->capacity() == numberRows + numberColumns && alternateWeights_->capacity() == numberRows + model_->factorization()->maximumPivots()) {
        //alternateWeights_->clear();
        if (pivotSequence_ >= 0 && pivotSequence_ < numberRows) {
#if ALT_UPDATE_WEIGHTS != 2
          // save pivot order
          CoinMemcpyN(pivotVariable,
            numberRows, alternateWeights_->getIndices());
#endif
          // pivotSequence_ stays as row index; all consumers (justSteepest,
          // justDevex, etc.) use it as a row index, not a sequence number.
        } else {
          pivotSequence_ = -1;
        }
        // indices_ now holds pivot variable sequence numbers (not sparse element
        // indices), so nElements_ must be 0.  mode=2 reads indices_ directly and
        // does not rely on nElements_, so this is safe.
        alternateWeights_->setNumElements(0);
        state_ = 1;
      } else {
        // size has changed - clear everything
        delete[] weights_;
        weights_ = NULL;
        delete infeasible_;
        infeasible_ = NULL;
        delete alternateWeights_;
        alternateWeights_ = NULL;
        delete[] savedWeights_;
        savedWeights_ = NULL;
        delete[] reference_;
        reference_ = NULL;
        state_ = -1;
        pivotSequence_ = -1;
      }
    }
  } else if (mode == 2 || mode == 4 || mode == 5) {
    // restore
    if (!weights_ || state_ == -1 || mode == 5) {
      // Partial is only allowed with certain types of matrix
      if ((mode_ != 4 && mode_ != 5) || numberSwitched_ || !model_->clpMatrix()->canDoPartialPricing()) {
        // initialize weights
        delete[] weights_;
        delete alternateWeights_;
        weights_ = new double[numberRows + numberColumns];
        alternateWeights_ = new CoinIndexedVector();
        // enough space so can use it for factorization
        alternateWeights_->reserve(numberRows + model_->factorization()->maximumPivots());
        initializeWeights();
        // create saved weights
        delete[] savedWeights_;
        savedWeights_ = CoinCopyOfArray(weights_, numberRows + numberColumns);
        // just do initialization
        mode = 3;
      } else {
        // Partial pricing
        // use region as somewhere to save non-fixed slacks
        // set up infeasibilities
        if (!infeasible_) {
          infeasible_ = new CoinIndexedVector();
          infeasible_->reserve(numberColumns + numberRows);
        }
        infeasible_->clear();
        int number = model_->numberRows() + model_->numberColumns();
        int iSequence;
        int numberLook = 0;
        int *which = infeasible_->getIndices();
        for (iSequence = model_->numberColumns(); iSequence < number; iSequence++) {
          ClpSimplex::Status status = model_->getStatus(iSequence);
          if (status != ClpSimplex::isFixed)
            which[numberLook++] = iSequence;
        }
        infeasible_->setNumElements(numberLook);
        doInfeasibilities = false;
      }
      savedPivotSequence_ = -2;
      savedSequenceOut_ = -2;
      if (pivotSequence_ < 0 || pivotSequence_ >= numberRows + numberColumns)
        pivotSequence_ = -1;

    } else {
      if (mode != 4) {
        // save
        CoinMemcpyN(weights_, (numberRows + numberColumns), savedWeights_);
        savedPivotSequence_ = pivotSequence_;
        savedSequenceOut_ = model_->sequenceOut();
      } else {
        // restore
        CoinMemcpyN(savedWeights_, (numberRows + numberColumns), weights_);
        // was - but I think should not be
        //pivotSequence_= savedPivotSequence_;
        //model_->setSequenceOut(savedSequenceOut_);
        pivotSequence_ = -1;
        model_->setSequenceOut(-1);
        // indices are wrong so clear by hand
        //alternateWeights_->clear();
        CoinZeroN(alternateWeights_->denseVector(),
          alternateWeights_->capacity());
        alternateWeights_->setNumElements(0);
      }
    }
    state_ = 0;
    // set up infeasibilities
    if (!infeasible_) {
      infeasible_ = new CoinIndexedVector();
      infeasible_->reserve(numberColumns + numberRows);
    }
  }
  if (mode >= 2 && mode != 5) {
    if (mode != 3) {
      if (pivotSequence_ >= 0) {
#if ALT_UPDATE_WEIGHTS != 2
        // restore pivot row
        int iRow;
        // permute alternateWeights
        double *temp = model_->rowArray(3)->denseVector();
        ;
        double *work = alternateWeights_->denseVector();
        int *savePivotOrder = model_->rowArray(3)->getIndices();
        int *oldPivotOrder = alternateWeights_->getIndices();
        for (iRow = 0; iRow < numberRows; iRow++) {
          int iPivot = oldPivotOrder[iRow];
          temp[iPivot] = work[iRow];
          savePivotOrder[iRow] = iPivot;
        }
        int number = 0;
        int found = -1;
        int *which = oldPivotOrder;
        // find pivot row and re-create alternateWeights
        for (iRow = 0; iRow < numberRows; iRow++) {
          int iPivot = pivotVariable[iRow];
          if (iPivot == pivotSequence_)
            found = iRow;
          work[iRow] = temp[iPivot];
          if (work[iRow])
            which[number++] = iRow;
        }
        alternateWeights_->setNumElements(number);
#ifdef CLP_DEBUG
        // Can happen but I should clean up
        assert(found >= 0);
#endif
        pivotSequence_ = found;
        for (iRow = 0; iRow < numberRows; iRow++) {
          int iPivot = savePivotOrder[iRow];
          temp[iPivot] = 0.0;
        }
#else
        for (int iRow = 0; iRow < numberRows; iRow++) {
          int iPivot = pivotVariable[iRow];
          if (iPivot == pivotSequence_) {
            pivotSequence_ = iRow;
            break;
          }
        }
#endif
      } else {
        // Just clean up
        /* If this happens when alternateWeights_ is
		       in "save" mode then alternateWeights_->clear()
		       is disastrous.
		       As will be fairly dense anyway and this
		       rarely happens just zero out */
        if (alternateWeights_ && alternateWeights_->getNumElements()) {
          //alternateWeights_->clear();
          CoinZeroN(alternateWeights_->denseVector(),
            alternateWeights_->capacity());
          alternateWeights_->setNumElements(0);
        }
      }
    }
    // Save size of factorization
    if (!model_->factorization()->pivots())
      sizeFactorization_ = model_->factorization()->numberElements();
    if (!doInfeasibilities)
      return; // don't disturb infeasibilities
    double *COIN_RESTRICT infeas = infeasible_->denseVector();
    int *COIN_RESTRICT index = infeasible_->getIndices();
    int numberNonZero = 0;
    infeasibilitiesState_ = 0;
    double tolerance = model_->currentDualTolerance();
    int number = model_->numberRows() + model_->numberColumns();
    int iSequence;

    double *reducedCost = model_->djRegion();
    const double *lower = model_->lowerRegion();
    const double *upper = model_->upperRegion();
    const double *solution = model_->solutionRegion();
    double primalTolerance = model_->currentPrimalTolerance();

    if (!model_->nonLinearCost()->lookBothWays()) {
      const unsigned char *COIN_RESTRICT status = model_->statusArray();
#ifndef CLP_PRIMAL_SLACK_MULTIPLIER
      for (iSequence = 0; iSequence < number; iSequence++) {
        double value = reducedCost[iSequence];
        infeas[iSequence] = 0.0;
        unsigned char thisStatus = status[iSequence] & 7;
        if (thisStatus == 3) {
        } else if ((thisStatus & 1) != 0) {
          // basic or fixed
          value = 0.0;
        } else if (thisStatus == 2) {
          value = -value;
        } else {
          // free or superbasic
          if (fabs(value) > FREE_ACCEPT * tolerance) {
            // we are going to bias towards free (but only if reasonable)
            value = -fabs(value) * FREE_BIAS;
          } else {
            value = 0.0;
          }
        }
        if (value < -tolerance) {
          infeas[iSequence] = value * value;
          index[numberNonZero++] = iSequence;
        }
      }
#else
      // Columns
      int numberColumns = model_->numberColumns();
      for (iSequence = 0; iSequence < numberColumns; iSequence++) {
        infeas[iSequence] = 0.0;
        double value = reducedCost[iSequence];
        unsigned char thisStatus = status[iSequence] & 7;
        if (thisStatus == 3) {
        } else if ((thisStatus & 1) != 0) {
          // basic or fixed
          value = 0.0;
        } else if (thisStatus == 2) {
          value = -value;
        } else {
          // free or superbasic
          if (fabs(value) > FREE_ACCEPT * tolerance) {
            // check hasn't slipped through
            if (solution[iSequence] < lower[iSequence] + primalTolerance) {
              model_->setStatus(iSequence, ClpSimplex::atLowerBound);
            } else if (solution[iSequence] > upper[iSequence] - primalTolerance) {
              model_->setStatus(iSequence, ClpSimplex::atUpperBound);
              value = -value;
            } else {
              // we are going to bias towards free (but only if reasonable)
              value = -fabs(value) * FREE_BIAS;
            }
          } else {
            value = 0.0;
          }
        }
        if (value < -tolerance) {
          infeas[iSequence] = value * value;
          index[numberNonZero++] = iSequence;
        }
      }
      // Rows
      for (; iSequence < number; iSequence++) {
        double value = reducedCost[iSequence];
        infeas[iSequence] = 0.0;
        unsigned char thisStatus = status[iSequence] & 7;
        if (thisStatus == 3) {
        } else if ((thisStatus & 1) != 0) {
          // basic or fixed
          value = 0.0;
        } else if (thisStatus == 2) {
          value = -value;
        } else {
          // free or superbasic
          if (fabs(value) > FREE_ACCEPT * tolerance) {
            // we are going to bias towards free (but only if reasonable)
            value = -fabs(value) * FREE_BIAS;
          } else {
            value = 0.0;
          }
        }
        if (value < -tolerance) {
          infeas[iSequence] = value * value;
          index[numberNonZero++] = iSequence;
        }
      }
#endif
      infeasible_->setNumElements(numberNonZero);
    } else {
      ClpNonLinearCost *nonLinear = model_->nonLinearCost();
      infeasible_->clear();
      // can go both ways
      for (iSequence = 0; iSequence < number; iSequence++) {
        double value = reducedCost[iSequence];
        ClpSimplex::Status status = model_->getStatus(iSequence);

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
            infeasible_->quickAdd(iSequence, value * value);
          }
          break;
        case ClpSimplex::atUpperBound:
          if (value > tolerance) {
            infeasible_->quickAdd(iSequence, value * value);
          } else {
            // look other way - change up should be negative
            value -= nonLinear->changeUpInCost(iSequence);
            if (value < -tolerance) {
              // store square in list
              infeasible_->quickAdd(iSequence, value * value);
            }
          }
          break;
        case ClpSimplex::atLowerBound:
          if (value < -tolerance) {
            infeasible_->quickAdd(iSequence, value * value);
          } else {
            // look other way - change down should be positive
            value -= nonLinear->changeDownInCost(iSequence);
            if (value > tolerance) {
              // store square in list
              infeasible_->quickAdd(iSequence, value * value);
            }
          }
        }
      }
    }
  }
}
// Gets rid of last update
void ClpPrimalColumnSteepest::unrollWeights()
{
  if ((mode_ == 4 || mode_ == 5) && !numberSwitched_)
    return;
  double *saved = alternateWeights_->denseVector();
  int number = alternateWeights_->getNumElements();
  int *which = alternateWeights_->getIndices();
  int i;
  for (i = 0; i < number; i++) {
    int iRow = which[i];
    weights_[iRow] = saved[iRow];
    saved[iRow] = 0.0;
  }
  alternateWeights_->setNumElements(0);
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpPrimalColumnPivot *ClpPrimalColumnSteepest::clone(bool CopyData) const
{
  if (CopyData) {
    return new ClpPrimalColumnSteepest(*this);
  } else {
    return new ClpPrimalColumnSteepest();
  }
}
void ClpPrimalColumnSteepest::updateWeights(CoinIndexedVector *input)
{
  // Local copy of mode so can decide what to do
  int switchType = mode_;
  if (mode_ == 4 && numberSwitched_)
    switchType = 3;
  else if (mode_ == 4 || mode_ == 5)
    return;
  int number = input->getNumElements();
  int *which = input->getIndices();
  double *work = input->denseVector();
#if ALT_UPDATE_WEIGHTS != 2
  int newNumber = 0;
  int *newWhich = alternateWeights_->getIndices();
  double *newWork = alternateWeights_->denseVector();
#endif
#ifdef ALT_UPDATE_WEIGHTSz
  {
    //int newNumber2 = 0;
    if (!altVector[0]) {
      altVector[0] = new CoinIndexedVector(2000);
      altVector[1] = new CoinIndexedVector(2000);
      altVector[2] = new CoinIndexedVector(2000);
    }
    altVector[0]->clear();
    // get updated pivot row
    int pivotRow = model_->pivotRow();
    // should it be - or what
    altVector[0]->quickAdd(pivotRow, model_->dualIn());
    model_->factorization()->updateColumnTranspose(altVector[2],
      altVector[0]);
    double *work2 = altVector[0]->denseVector();
    //altVector[1]->clear();
    int *newWhich2 = altVector[1]->getIndices();
    double *newWork2 = altVector[1]->denseVector();
    int number2 = altVector[1]->getNumElements();
    int nRow = model_->numberRows();
    int nCol = model_->numberColumns();
    int nTotal = nRow + nCol;
    double *temp = new double[2 * nTotal + nRow];
    memset(temp, 0, (2 * nTotal + nRow) * sizeof(double));
    double *pivRow = temp + nTotal;
    double *temp2 = temp + nCol;
    double *temp2P = pivRow + nCol;
    double *piU = pivRow + nTotal;
    double devex = 0.0;
    double scaleFactor = 1.0 / model_->dualIn();
    const int *pivotVariable = model_->pivotVariable();
    for (int i = 0; i < number; i++) {
      int iRow = which[i];
      int iPivot = pivotVariable[iRow];
      if (reference(iPivot)) {
        devex += work[iRow] * work[iRow];
      }
    }
    int sequenceIn = model_->sequenceIn();
    int sequenceOut = model_->sequenceOut();
    for (int i = 0; i < number2; i++) {
      int iRow = newWhich2[i];
      temp2[iRow] = newWork2[iRow];
    }
    //if (!input->packedMode()) {
    for (int i = 0; i < number; i++) {
      int iRow = which[i];
      piU[iRow] = work2[iRow];
      temp2P[iRow] = work2[iRow];
    }
    double alpha = model_->alpha();
    const int *row = model_->matrix()->getIndices();
    const CoinBigIndex *columnStart = model_->matrix()->getVectorStarts();
    const int *columnLength = model_->matrix()->getVectorLengths();
    const double *element = model_->matrix()->getElements();
    for (int i = 0; i < nCol; i++) {
      CoinBigIndex start = columnStart[i];
      CoinBigIndex end = start + columnLength[i];
      double value = 0.0;
      double value2 = 0.0;
      for (CoinBigIndex j = start; j < end; j++) {
        int iRow = row[j];
        value -= piU[iRow] * element[j];
        value2 -= newWork2[iRow] * element[j];
      }
      pivRow[i] = value;
      temp[i] = value2;
    }
    const unsigned char *COIN_RESTRICT statusArray = model_->statusArray();
    for (int i = 0; i < nTotal; i++) {
      unsigned char thisStatus = statusArray[i] & 7;
      if (thisStatus != 1 && thisStatus != 5) {
        double pivot = pivRow[i] * scaleFactor;
        double modification = temp[i];
        double thisWeight = weights_[i];
        double pivotSquared = pivot * pivot;
        double newWeight = thisWeight + pivotSquared * devexA - 2.0 * pivot * modification;
        temp[i] = newWeight;
      } else {
        temp[i] = 1.0;
      }
    }
    temp[sequenceOut] = devexA / (alpha * alpha);
    // to keep clean for debug
#ifndef NDEBUG
    {
      if (sequenceOut < nCol) {
        if (model_->columnLower()[sequenceOut] == model_->columnUpper()[sequenceOut])
          temp[sequenceOut] = 1.0;
      } else {
        if (model_->rowLower()[sequenceOut - nCol] == model_->rowUpper()[sequenceOut - nCol])
          temp[sequenceOut] = 1.0;
      }
    }
#endif
    temp[sequenceIn] = 1.0;
    for (int i = 0; i < nTotal; i++) {
      printf("%g ", temp[i]);
      if ((i % 10) == 9)
        printf("\n");
    }
    if (((nTotal - 1) % 10) != 9)
      printf("\n");
    delete[] temp;
  }
#endif
  int i;
  int sequenceIn = model_->sequenceIn();
  int sequenceOut = model_->sequenceOut();
  const int *pivotVariable = model_->pivotVariable();

  int pivotRow = model_->pivotRow();
  pivotSequence_ = pivotRow;

  devex_ = 0.0;
  // Can't create alternateWeights_ as packed as needed unpacked
  if (!input->packedMode()) {
    if (pivotRow >= 0) {
      if (switchType == 1) {
        for (i = 0; i < number; i++) {
          int iRow = which[i];
          devex_ += work[iRow] * work[iRow];
#if ALT_UPDATE_WEIGHTS != 2
          newWork[iRow] = -2.0 * work[iRow];
#endif
        }
#if ALT_UPDATE_WEIGHTS != 2
        newWork[pivotRow] = -2.0 * std::max(devex_, 0.0);
#endif
        devex_ += ADD_ONE;
        weights_[sequenceOut] = 1.0 + ADD_ONE;
#if ALT_UPDATE_WEIGHTS != 2
        CoinMemcpyN(which, number, newWhich);
        alternateWeights_->setNumElements(number);
#endif
      } else {
        if ((mode_ != 4 && mode_ != 5) || numberSwitched_ > 1) {
          for (i = 0; i < number; i++) {
            int iRow = which[i];
            int iPivot = pivotVariable[iRow];
            if (reference(iPivot)) {
              devex_ += work[iRow] * work[iRow];
#if ALT_UPDATE_WEIGHTS != 2
              newWork[iRow] = -2.0 * work[iRow];
              newWhich[newNumber++] = iRow;
#endif
            }
          }
#if ALT_UPDATE_WEIGHTS != 2
          if (!newWork[pivotRow] && devex_ > 0.0)
            newWhich[newNumber++] = pivotRow; // add if not already in
          newWork[pivotRow] = -2.0 * std::max(devex_, 0.0);
#endif
        } else {
          for (i = 0; i < number; i++) {
            int iRow = which[i];
            int iPivot = pivotVariable[iRow];
            if (reference(iPivot))
              devex_ += work[iRow] * work[iRow];
          }
        }
        if (reference(sequenceIn)) {
          devex_ += 1.0;
        } else {
        }
        if (reference(sequenceOut)) {
          weights_[sequenceOut] = 1.0 + 1.0;
        } else {
          weights_[sequenceOut] = 1.0;
        }
#if ALT_UPDATE_WEIGHTS != 2
        alternateWeights_->setNumElements(newNumber);
#endif
      }
    } else {
      if (switchType == 1) {
        for (i = 0; i < number; i++) {
          int iRow = which[i];
          devex_ += work[iRow] * work[iRow];
        }
        devex_ += ADD_ONE;
      } else {
        for (i = 0; i < number; i++) {
          int iRow = which[i];
          int iPivot = pivotVariable[iRow];
          if (reference(iPivot)) {
            devex_ += work[iRow] * work[iRow];
          }
        }
        if (reference(sequenceIn))
          devex_ += 1.0;
      }
    }
  } else {
    // packed input
    if (pivotRow >= 0) {
      if (switchType == 1) {
        for (i = 0; i < number; i++) {
#if ALT_UPDATE_WEIGHTS != 2
          int iRow = which[i];
#endif
          devex_ += work[i] * work[i];
#if ALT_UPDATE_WEIGHTS != 2
          newWork[iRow] = -2.0 * work[i];
#endif
        }
#if ALT_UPDATE_WEIGHTS != 2
        newWork[pivotRow] = -2.0 * std::max(devex_, 0.0);
#endif
        devex_ += ADD_ONE;
        weights_[sequenceOut] = 1.0 + ADD_ONE;
#if ALT_UPDATE_WEIGHTS != 2
        CoinMemcpyN(which, number, newWhich);
        alternateWeights_->setNumElements(number);
#endif
      } else {
        if ((mode_ != 4 && mode_ != 5) || numberSwitched_ > 1) {
          for (i = 0; i < number; i++) {
            int iRow = which[i];
            int iPivot = pivotVariable[iRow];
            if (reference(iPivot)) {
              devex_ += work[i] * work[i];
#if ALT_UPDATE_WEIGHTS != 2
              newWork[iRow] = -2.0 * work[i];
              newWhich[newNumber++] = iRow;
#endif
            }
          }
#if ALT_UPDATE_WEIGHTS != 2
          if (!newWork[pivotRow] && devex_ > 0.0)
            newWhich[newNumber++] = pivotRow; // add if not already in
          newWork[pivotRow] = -2.0 * std::max(devex_, 0.0);
#endif
        } else {
          for (i = 0; i < number; i++) {
            int iRow = which[i];
            int iPivot = pivotVariable[iRow];
            if (reference(iPivot))
              devex_ += work[i] * work[i];
          }
        }
        if (reference(sequenceIn)) {
          devex_ += 1.0;
        } else {
        }
        if (reference(sequenceOut)) {
          weights_[sequenceOut] = 1.0 + 1.0;
        } else {
          weights_[sequenceOut] = 1.0;
        }
#if ALT_UPDATE_WEIGHTS != 2
        alternateWeights_->setNumElements(newNumber);
#endif
      }
    } else {
      if (switchType == 1) {
        for (i = 0; i < number; i++) {
          devex_ += work[i] * work[i];
        }
        devex_ += ADD_ONE;
      } else {
        for (i = 0; i < number; i++) {
          int iRow = which[i];
          int iPivot = pivotVariable[iRow];
          if (reference(iPivot)) {
            devex_ += work[i] * work[i];
          }
        }
        if (reference(sequenceIn))
          devex_ += 1.0;
      }
    }
  }
  if (devex_ < 1.001e-30) {
    COIN_DETAIL_PRINT(printf("devex of incoming tiny %d %g\n", sequenceIn, devex_));
    devex_ = 1.0e-30;
  }
  double oldDevex = weights_[sequenceIn];
#ifdef CLP_DEBUG
  if ((model_->messageHandler()->logLevel() & 32))
    printf("old weight %g, new %g\n", oldDevex, devex_);
#endif
  double check = std::max(devex_, oldDevex) + 0.1;
  weights_[sequenceIn] = devex_;
  double testValue = 0.1;
  if (mode_ == 4 && numberSwitched_ == 1)
    testValue = 0.5;
  if (fabs(devex_ - oldDevex) > testValue * check) {
#ifdef CLP_DEBUG
    if ((model_->messageHandler()->logLevel() & 48) == 16)
      printf("old weight %g, new %g\n", oldDevex, devex_);
#endif
    //printf("old weight %g, new %g\n",oldDevex,devex_);
    testValue = 0.99;
    if (mode_ == 1)
      testValue = 1.01e1; // make unlikely to do if steepest
    else if (mode_ == 4 && numberSwitched_ == 1)
      testValue = 0.9;
    double difference = fabs(devex_ - oldDevex);
    if (difference > testValue * check) {
      // need to redo
      model_->messageHandler()->message(CLP_INITIALIZE_STEEP,
        *model_->messagesPointer())
        << oldDevex << devex_
        << CoinMessageEol;
      initializeWeights();
      // redo devex_
      if (pivotRow >= 0)
        devex_ = 1.0;
    }
  }
  if (pivotRow >= 0) {
    // set outgoing weight here
    double alpha = model_->alpha();
    if (fabs(alpha) > 1.0e15) {
      COIN_DETAIL_PRINT(printf("alpha %g for %d !!\n", alpha, model_->sequenceOut()));
      alpha = 1.0e15;
    }
    weights_[model_->sequenceOut()] = devex_ / (alpha * alpha);
  }
}
// Checks accuracy - just for debug
void ClpPrimalColumnSteepest::checkAccuracy(int sequence,
  double relativeTolerance,
  CoinIndexedVector *rowArray1,
  CoinIndexedVector *rowArray2)
{
  if ((mode_ == 4 || mode_ == 5) && !numberSwitched_)
    return;
  model_->unpack(rowArray1, sequence);
  model_->factorization()->updateColumn(rowArray2, rowArray1);
  int number = rowArray1->getNumElements();
  int *which = rowArray1->getIndices();
  double *work = rowArray1->denseVector();
  const int *pivotVariable = model_->pivotVariable();

  double devex = 0.0;
  int i;

  if (mode_ == 1) {
    for (i = 0; i < number; i++) {
      int iRow = which[i];
      devex += work[iRow] * work[iRow];
      work[iRow] = 0.0;
    }
    devex += ADD_ONE;
  } else {
    for (i = 0; i < number; i++) {
      int iRow = which[i];
      int iPivot = pivotVariable[iRow];
      if (reference(iPivot)) {
        devex += work[iRow] * work[iRow];
      }
      work[iRow] = 0.0;
    }
    if (reference(sequence))
      devex += 1.0;
  }

  double oldDevex = std::max(weights_[sequence], 1.0e-4);
  devex = std::max(devex, 1.0e-4);
  double check = std::max(devex, oldDevex);
  ;
  rowArray1->setNumElements(0);
  if (fabs(devex - oldDevex) > relativeTolerance * check) {
    //COIN_DETAIL_PRINT(printf("check %d old weight %g, new %g\n", sequence, oldDevex, devex));
    printf("check %d old weight %g, new %g\n", sequence, oldDevex, devex);
    if (mode_ == 0) {
      rowArray1->setNumElements(0);
      model_->unpack(rowArray1, sequence);
      number = rowArray1->getNumElements();
      for (i = 0; i < number; i++)
        printf("(%d,%g) ", which[i], work[which[i]]);
      printf("\n");
      model_->factorization()->updateColumn(rowArray2, rowArray1);
      number = rowArray1->getNumElements();
      for (i = 0; i < number; i++)
        printf("(%d,%g) ", which[i], work[which[i]]);
      printf("\n");
      devex = 0.0;
      for (i = 0; i < number; i++) {
        int iRow = which[i];
        int iPivot = pivotVariable[iRow];
        if (reference(iPivot)) {
          devex += work[iRow] * work[iRow];
        }
        work[iRow] = 0.0;
      }
      if (reference(sequence))
        devex += 1.0;
    }
    // update so won't print again
    weights_[sequence] = devex;
  }
}

// Initialize weights
void ClpPrimalColumnSteepest::initializeWeights()
{
  int numberRows = model_->numberRows();
  int numberColumns = model_->numberColumns();
  int number = numberRows + numberColumns;
  int iSequence;
  if (mode_ != 1) {
    // initialize to 1.0
    // and set reference framework
    if (!reference_) {
      int nWords = (number + 31) >> 5;
      reference_ = new unsigned int[nWords];
      CoinZeroN(reference_, nWords);
    }

    for (iSequence = 0; iSequence < number; iSequence++) {
      weights_[iSequence] = 1.0;
      if (model_->getStatus(iSequence) == ClpSimplex::basic) {
        setReference(iSequence, false);
      } else {
        setReference(iSequence, true);
      }
    }
  } else {
    CoinIndexedVector *temp = new CoinIndexedVector();
    temp->reserve(numberRows + model_->factorization()->maximumPivots());
    double *array = alternateWeights_->denseVector();
    int *which = alternateWeights_->getIndices();

    for (iSequence = 0; iSequence < number; iSequence++) {
      weights_[iSequence] = 1.0 + ADD_ONE;
      if (model_->getStatus(iSequence) != ClpSimplex::basic && model_->getStatus(iSequence) != ClpSimplex::isFixed) {
        model_->unpack(alternateWeights_, iSequence);
        double value = ADD_ONE;
        model_->factorization()->updateColumn(temp, alternateWeights_);
        int number = alternateWeights_->getNumElements();
        int j;
        for (j = 0; j < number; j++) {
          int iRow = which[j];
          value += array[iRow] * array[iRow];
          array[iRow] = 0.0;
        }
        alternateWeights_->setNumElements(0);
        weights_[iSequence] = value;
      }
    }
    delete temp;
  }
}
// Gets rid of all arrays
void ClpPrimalColumnSteepest::clearArrays()
{
  if (persistence_ == normal) {
    delete[] weights_;
    weights_ = NULL;
    delete infeasible_;
    infeasible_ = NULL;
    delete alternateWeights_;
    alternateWeights_ = NULL;
    delete[] savedWeights_;
    savedWeights_ = NULL;
    delete[] reference_;
    reference_ = NULL;
  }
  pivotSequence_ = -1;
  state_ = -1;
  savedPivotSequence_ = -1;
  savedSequenceOut_ = -1;
  devex_ = 0.0;
}
// Returns true if would not find any column
bool ClpPrimalColumnSteepest::looksOptimal() const
{
  if (looksOptimal_)
    return true; // user overrode
  //**** THIS MUST MATCH the action coding above
  double tolerance = model_->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = std::min(1.0e-2, model_->largestDualError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance + error;
  if (model_->numberIterations() < model_->lastBadIteration() + 200) {
    // we can't really trust infeasibilities if there is dual error
    double checkTolerance = 1.0e-8;
    if (!model_->factorization()->pivots())
      checkTolerance = 1.0e-6;
    if (model_->largestDualError() > checkTolerance)
      tolerance *= model_->largestDualError() / checkTolerance;
    // But cap
    tolerance = std::min(1000.0, tolerance);
  }
  int number = model_->numberRows() + model_->numberColumns();
  int iSequence;

  double *reducedCost = model_->djRegion();
  int numberInfeasible = 0;
  if (!model_->nonLinearCost()->lookBothWays()) {
    for (iSequence = 0; iSequence < number; iSequence++) {
      double value = reducedCost[iSequence];
      ClpSimplex::Status status = model_->getStatus(iSequence);

      switch (status) {

      case ClpSimplex::basic:
      case ClpSimplex::isFixed:
        break;
      case ClpSimplex::isFree:
      case ClpSimplex::superBasic:
        if (fabs(value) > FREE_ACCEPT * tolerance)
          numberInfeasible++;
        break;
      case ClpSimplex::atUpperBound:
        if (value > tolerance)
          numberInfeasible++;
        break;
      case ClpSimplex::atLowerBound:
        if (value < -tolerance)
          numberInfeasible++;
      }
    }
  } else {
    ClpNonLinearCost *nonLinear = model_->nonLinearCost();
    // can go both ways
    for (iSequence = 0; iSequence < number; iSequence++) {
      double value = reducedCost[iSequence];
      ClpSimplex::Status status = model_->getStatus(iSequence);

      switch (status) {

      case ClpSimplex::basic:
      case ClpSimplex::isFixed:
        break;
      case ClpSimplex::isFree:
      case ClpSimplex::superBasic:
        if (fabs(value) > FREE_ACCEPT * tolerance)
          numberInfeasible++;
        break;
      case ClpSimplex::atUpperBound:
        if (value > tolerance) {
          numberInfeasible++;
        } else {
          // look other way - change up should be negative
          value -= nonLinear->changeUpInCost(iSequence);
          if (value < -tolerance)
            numberInfeasible++;
        }
        break;
      case ClpSimplex::atLowerBound:
        if (value < -tolerance) {
          numberInfeasible++;
        } else {
          // look other way - change down should be positive
          value -= nonLinear->changeDownInCost(iSequence);
          if (value > tolerance)
            numberInfeasible++;
        }
      }
    }
  }
  return numberInfeasible == 0;
}
/* Returns number of extra columns for sprint algorithm - 0 means off.
   Also number of iterations before recompute
*/
int ClpPrimalColumnSteepest::numberSprintColumns(int &numberIterations) const
{
  numberIterations = 0;
  int numberAdd = 0;
  if (!numberSwitched_ && mode_ >= 10) {
    numberIterations = std::min(2000, model_->numberRows() / 5);
    numberIterations = std::max(numberIterations, model_->factorizationFrequency());
    numberIterations = std::max(numberIterations, 500);
    if (mode_ == 10) {
      numberAdd = std::max(300, model_->numberColumns() / 10);
      numberAdd = std::max(numberAdd, model_->numberRows() / 5);
      // fake all
      //numberAdd=1000000;
      numberAdd = std::min(numberAdd, model_->numberColumns());
    } else {
      abort();
    }
  }
  return numberAdd;
}
// Switch off sprint idea
void ClpPrimalColumnSteepest::switchOffSprint()
{
  numberSwitched_ = 10;
}
// Update djs doing partial pricing (dantzig)
int ClpPrimalColumnSteepest::partialPricing(CoinIndexedVector *updates,
  CoinIndexedVector *spareRow2,
  int numberWanted,
  int numberLook)
{
  int number = 0;
  int *index;
  double *updateBy;
  double *reducedCost;
  double saveTolerance = model_->currentDualTolerance();
  double tolerance = model_->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = std::min(1.0e-2, model_->largestDualError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance + error;
  if (model_->numberIterations() < model_->lastBadIteration() + 200) {
    // we can't really trust infeasibilities if there is dual error
    double checkTolerance = 1.0e-8;
    if (!model_->factorization()->pivots())
      checkTolerance = 1.0e-6;
    if (model_->largestDualError() > checkTolerance)
      tolerance *= model_->largestDualError() / checkTolerance;
    // But cap
    tolerance = std::min(1000.0, tolerance);
  }
  if (model_->factorization()->pivots() && model_->numberPrimalInfeasibilities())
    tolerance = std::max(tolerance, 1.0e-15 * model_->infeasibilityCost());
  // So partial pricing can use
  model_->setCurrentDualTolerance(tolerance);
  model_->factorization()->updateColumnTranspose(spareRow2, updates);
  int numberColumns = model_->numberColumns();

  // Rows
  reducedCost = model_->djRegion(0);

  number = updates->getNumElements();
  index = updates->getIndices();
  updateBy = updates->denseVector();
  int j;
  double *duals = model_->dualRowSolution();
  for (j = 0; j < number; j++) {
    int iSequence = index[j];
    double value = duals[iSequence];
    value -= updateBy[j];
    updateBy[j] = 0.0;
    duals[iSequence] = value;
  }
  //#define CLP_DEBUG
#ifdef CLP_DEBUG
  // check duals
  {
    int numberRows = model_->numberRows();
    //work space
    CoinIndexedVector arrayVector;
    arrayVector.reserve(numberRows + 1000);
    CoinIndexedVector workSpace;
    workSpace.reserve(numberRows + 1000);

    int iRow;
    double *array = arrayVector.denseVector();
    int *index = arrayVector.getIndices();
    int number = 0;
    int *pivotVariable = model_->pivotVariable();
    double *cost = model_->costRegion();
    for (iRow = 0; iRow < numberRows; iRow++) {
      int iPivot = pivotVariable[iRow];
      double value = cost[iPivot];
      if (value) {
        array[iRow] = value;
        index[number++] = iRow;
      }
    }
    arrayVector.setNumElements(number);
    // Extended duals before "updateTranspose"
    model_->clpMatrix()->dualExpanded(model_, &arrayVector, NULL, 0);

    // Btran basic costs
    model_->factorization()->updateColumnTranspose(&workSpace, &arrayVector);

    // now look at dual solution
    for (iRow = 0; iRow < numberRows; iRow++) {
      // slack
      double value = array[iRow];
      if (fabs(duals[iRow] - value) > 1.0e-3)
        printf("bad row %d old dual %g new %g\n", iRow, duals[iRow], value);
      //duals[iRow]=value;
    }
  }
#endif
#undef CLP_DEBUG
  double bestDj = tolerance;
  int bestSequence = -1;

  const double *cost = model_->costRegion(1);

  model_->clpMatrix()->setOriginalWanted(numberWanted);
  model_->clpMatrix()->setCurrentWanted(numberWanted);
  int iPassR = 0, iPassC = 0;
  // Setup two passes
  // This biases towards picking row variables
  // This probably should be fixed
  int startR[4];
  int numberRows = model_->numberRows();
  startR[1] = numberColumns+numberRows;
  startR[2] = numberColumns;
  double randomR = model_->randomNumberGenerator()->randomDouble();
  double dstart = static_cast<double> (numberRows) * randomR;
  startR[0] = numberColumns + static_cast<int> (dstart);
  startR[3] = startR[0];
  double startC[4];
  startC[1] = 1.0;
  startC[2] = 0;
  double randomC = model_->randomNumberGenerator()->randomDouble();
  startC[0] = randomC;
  startC[3] = randomC;
  reducedCost = model_->djRegion(1);
  int sequenceOut = model_->sequenceOut();
  double *duals2 = duals - numberColumns;
  int chunk = std::min(1024, (numberColumns + numberRows) / 32);
#ifdef COIN_DETAIL
  if (model_->numberIterations() % 1000 == 0 && model_->logLevel() > 1) {
    printf("%d wanted, nSlacks %d, chunk %d\n", numberWanted, nSlacks, chunk);
    int i;
    for (i = 0; i < 4; i++)
      printf("start R %d C %g ", startR[i], startC[i]);
    printf("\n");
  }
#endif
  chunk = std::max(chunk, 256);
  bool finishedR = false, finishedC = false;
  bool doingR = randomR > randomC;
  //doingR=false;
  int saveNumberWanted = numberWanted;
  while (!finishedR || !finishedC) {
    if (finishedR)
      doingR = false;
    if (doingR) {
      int saveSequence = bestSequence;
      int start = startR[iPassR];
      int end = std::min(startR[iPassR + 1], start + chunk / 2);
      int iSequence;
      for (iSequence = start; iSequence < end; iSequence++) {
	assert (iSequence>=numberColumns);
        if (iSequence != sequenceOut) {
          double value;
          ClpSimplex::Status status = model_->getStatus(iSequence);

          switch (status) {

          case ClpSimplex::basic:
          case ClpSimplex::isFixed:
            break;
          case ClpSimplex::isFree:
          case ClpSimplex::superBasic:
            value = fabs(cost[iSequence] + duals2[iSequence]);
            if (value > FREE_ACCEPT * tolerance) {
              numberWanted--;
              // we are going to bias towards free (but only if reasonable)
              value *= FREE_BIAS;
              if (value > bestDj) {
                // check flagged variable and correct dj
                if (!model_->flagged(iSequence)) {
                  bestDj = value;
                  bestSequence = iSequence;
                } else {
                  // just to make sure we don't exit before got something
                  numberWanted++;
                }
              }
            }
            break;
          case ClpSimplex::atUpperBound:
            value = cost[iSequence] + duals2[iSequence];
            if (value > tolerance) {
              numberWanted--;
              if (value > bestDj) {
                // check flagged variable and correct dj
                if (!model_->flagged(iSequence)) {
                  bestDj = value;
                  bestSequence = iSequence;
                } else {
                  // just to make sure we don't exit before got something
                  numberWanted++;
                }
              }
            }
            break;
          case ClpSimplex::atLowerBound:
            value = -(cost[iSequence] + duals2[iSequence]);
            if (value > tolerance) {
              numberWanted--;
              if (value > bestDj) {
                // check flagged variable and correct dj
                if (!model_->flagged(iSequence)) {
                  bestDj = value;
                  bestSequence = iSequence;
                } else {
                  // just to make sure we don't exit before got something
                  numberWanted++;
                }
              }
            }
            break;
          }
        }
        if (!numberWanted)
          break;
      }
      numberLook -= (end - start);
      if (numberLook < 0 && (10 * (saveNumberWanted - numberWanted) > saveNumberWanted))
        numberWanted = 0; // give up
      if (saveSequence != bestSequence) {
        // dj
        reducedCost[bestSequence] = cost[bestSequence] + duals[bestSequence - numberColumns];
        bestDj = fabs(reducedCost[bestSequence]);
        model_->clpMatrix()->setSavedBestSequence(bestSequence);
        model_->clpMatrix()->setSavedBestDj(reducedCost[bestSequence]);
      }
      model_->clpMatrix()->setCurrentWanted(numberWanted);
      if (!numberWanted)
        break;
      doingR = false;
      // update start
      startR[iPassR] = iSequence;
      if (iSequence >= startR[iPassR + 1]) {
        if (iPassR)
          finishedR = true;
        else
          iPassR = 2;
      }
    }
    if (finishedC)
      doingR = true;
    if (!doingR) {
      int saveSequence = bestSequence;
      // Columns
      double start = startC[iPassC];
      // If we put this idea back then each function needs to update endFraction **
      double end = startC[iPassC + 1]; // force end
      model_->clpMatrix()->partialPricing(model_, start, end, bestSequence, numberWanted);
      numberWanted = model_->clpMatrix()->currentWanted();
      numberLook -= static_cast< int >((end - start) * numberColumns);
      if (numberLook < 0 && (10 * (saveNumberWanted - numberWanted) > saveNumberWanted))
        numberWanted = 0; // give up
      if (saveSequence != bestSequence) {
        // dj
        bestDj = fabs(model_->clpMatrix()->reducedCost(model_, bestSequence));
      }
      if (!numberWanted)
        break;
      doingR = true;
      // update start
      startC[iPassC] = end;
      if (end >= startC[iPassC + 1] - 1.0e-8) {
        if (iPassC)
          finishedC = true;
        else
          iPassC = 2;
      }
    }
  }
  updates->setNumElements(0);

  // Restore tolerance
  model_->setCurrentDualTolerance(saveTolerance);
  // Now create variable if column generation
  model_->clpMatrix()->createVariable(model_, bestSequence);
#ifndef NDEBUG
  if (bestSequence >= 0) {
    if (model_->getStatus(bestSequence) == ClpSimplex::atLowerBound)
      assert(model_->reducedCost(bestSequence) < 0.0);
    if (model_->getStatus(bestSequence) == ClpSimplex::atUpperBound)
      assert(model_->reducedCost(bestSequence) > 0.0);
  }
#endif
  return bestSequence;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
