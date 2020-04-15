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
//#define CBC_DEBUG
//#define TRACE_ONE 19
#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchDynamic.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

// Removing magic constants.

// This is a very small number, added to something to make sure it's non-zero.
// Useful, for example in denominators of ratios to avoid any possible division by zero
#define nonZeroAmount 1.0e-30

// Increasing the size of an array when it grows to the end of its alloted space.
// In this file, only used for the history of the outcome of a branch.
// New size is   size_scale_numerator* <old value> / size_scale_denominator + additive_size_increase.

#define size_scale_numerator 3
#define size_scale_denominator 2
#define additive_size_increase 100

// Explanation of options used in this file

// TYPE2 defines a strategy for computing pseudocosts
// 0 means to just use the absolute change in objective
// 1 means use the relative change in objective
// 2 means use a convex combination of absolute and relative objective changes

// For option 2 (TYPE2 == 2), the specific combination is controlled by TYPERATIO
// Includes a TYPERATIO fraction of the absolute change and (1 - TYPERATIO) fraction of
// the relative change.  So should in general have 0 <= TYPERATIO <= 1.  But for the
// equality cases, you're better off using the other strategy (TYPE2) options.

#ifdef COIN_DEVELOP
typedef struct {
  double sumUp_;
  double upEst_; // or change in obj in update
  double sumDown_;
  double downEst_; // or movement in value in update
  int sequence_;
  int numberUp_;
  int numberUpInf_;
  int numberDown_;
  int numberDownInf_;
  char where_;
  char status_;
} History;
static History *history = NULL;
static int numberHistory = 0;
static int maxHistory = 0;
static bool getHistoryStatistics_ = true;
static void increaseHistory()
{
  if (numberHistory == maxHistory) {
    // This was originally 3 * maxHistory/2 + 100
    maxHistory = additive_size_increase + (size_scale_numerator * maxHistory) / size_scale_denominator;
    History *temp = new History[maxHistory];
    memcpy(temp, history, numberHistory * sizeof(History));
    delete[] history;
    history = temp;
  }
}
static bool addRecord(History newOne)
{
  //if (!getHistoryStatistics_)
  return false;
  bool fromCompare = false;
  int i;
  for (i = numberHistory - 1; i >= 0; i--) {
    if (newOne.sequence_ != history[i].sequence_)
      continue;
    if (newOne.where_ != history[i].where_)
      continue;
    if (newOne.numberUp_ != history[i].numberUp_)
      continue;
    if (newOne.sumUp_ != history[i].sumUp_)
      continue;
    if (newOne.numberUpInf_ != history[i].numberUpInf_)
      continue;
    if (newOne.upEst_ != history[i].upEst_)
      continue;
    if (newOne.numberDown_ != history[i].numberDown_)
      continue;
    if (newOne.sumDown_ != history[i].sumDown_)
      continue;
    if (newOne.numberDownInf_ != history[i].numberDownInf_)
      continue;
    if (newOne.downEst_ != history[i].downEst_)
      continue;
    // If B knock out previous B
    if (newOne.where_ == 'C') {
      fromCompare = true;
      if (newOne.status_ == 'B') {
        int j;
        for (j = i - 1; j >= 0; j--) {
          if (history[j].where_ == 'C') {
            if (history[j].status_ == 'I') {
              break;
            } else if (history[j].status_ == 'B') {
              history[j].status_ = ' ';
              break;
            }
          }
        }
      }
      break;
    }
  }
  if (i == -1 || fromCompare) {
    //add
    increaseHistory();
    history[numberHistory++] = newOne;
    return true;
  } else {
    return false;
  }
}
#endif

// Default Constructor
CbcBranchDynamicDecision::CbcBranchDynamicDecision()
  : CbcBranchDecision()
{
  bestCriterion_ = 0.0;
  bestChangeUp_ = 0.0;
  bestNumberUp_ = 0;
  bestChangeDown_ = 0.0;
  bestNumberDown_ = 0;
  bestObject_ = NULL;
}

// Copy constructor
CbcBranchDynamicDecision::CbcBranchDynamicDecision(
  const CbcBranchDynamicDecision &rhs)
  : CbcBranchDecision()
{
  bestCriterion_ = rhs.bestCriterion_;
  bestChangeUp_ = rhs.bestChangeUp_;
  bestNumberUp_ = rhs.bestNumberUp_;
  bestChangeDown_ = rhs.bestChangeDown_;
  bestNumberDown_ = rhs.bestNumberDown_;
  bestObject_ = rhs.bestObject_;
}

CbcBranchDynamicDecision::~CbcBranchDynamicDecision()
{
}

// Clone
CbcBranchDecision *
CbcBranchDynamicDecision::clone() const
{
  return new CbcBranchDynamicDecision(*this);
}

// Initialize i.e. before start of choosing at a node
void CbcBranchDynamicDecision::initialize(CbcModel * /*model*/)
{
  bestCriterion_ = 0.0;
  bestChangeUp_ = 0.0;
  bestNumberUp_ = 0;
  bestChangeDown_ = 0.0;
  bestNumberDown_ = 0;
  bestObject_ = NULL;
#ifdef COIN_DEVELOP
  History hist;
  hist.where_ = 'C';
  hist.status_ = 'I';
  hist.sequence_ = 55555;
  hist.numberUp_ = 0;
  hist.numberUpInf_ = 0;
  hist.sumUp_ = 0.0;
  hist.upEst_ = 0.0;
  hist.numberDown_ = 0;
  hist.numberDownInf_ = 0;
  hist.sumDown_ = 0.0;
  hist.downEst_ = 0.0;
  addRecord(hist);
#endif
}

/* Saves a clone of current branching object.  Can be used to update
      information on object causing branch - after branch */
void CbcBranchDynamicDecision::saveBranchingObject(OsiBranchingObject *object)
{
  OsiBranchingObject *obj = object->clone();
#ifndef NDEBUG
  CbcBranchingObject *obj2 = dynamic_cast< CbcBranchingObject * >(obj);
  assert(obj2);
#if COIN_DEVELOP > 1
  CbcDynamicPseudoCostBranchingObject *branchingObject = dynamic_cast< CbcDynamicPseudoCostBranchingObject * >(obj);
  if (!branchingObject)
    printf("no dynamic branching object Dynamic Decision\n");
#endif
#else
  CbcBranchingObject *obj2 = static_cast< CbcBranchingObject * >(obj);
#endif
  //object_=branchingObject;
  object_ = obj2;
}
/* Pass in information on branch just done.
   assumes object can get information from solver */
/*
  The expectation is that this method will be called after the branch has been
  imposed on the constraint system and resolve() has executed.

  Note that the CbcBranchDecision is a property of the CbcModel. Note also that
  this method is reaching right through the CbcBranchingObject to update
  information in the underlying CbcObject. That's why we delete the
  branchingObject at the end of the method --- the next time we're called,
  the CbcObject will be different.
*/
void CbcBranchDynamicDecision::updateInformation(OsiSolverInterface *solver,
  const CbcNode *node)
{
  assert(object_);
  const CbcModel *model = object_->model();
  double originalValue = node->objectiveValue();
  int originalUnsatisfied = node->numberUnsatisfied();
  double objectiveValue = solver->getObjValue() * model->getObjSense();
  int unsatisfied = 0;
  int i;
  int numberIntegers = model->numberIntegers();
  ;
  const double *solution = solver->getColSolution();
  //const double * lower = solver->getColLower();
  //const double * upper = solver->getColUpper();
  /*
	 Gain access to the associated CbcBranchingObject and its underlying
	 CbcObject.

	 Seems like we'd want to distinguish between no branching object and a
	 branching object of the wrong type. Just deleting an object of the wrong
	 type hides many sins.

	 Hmmm ... if we're using the OSI side of the hierarchy, is this indicated by a
	 null object_? Nah, then we have an assert failure off the top.
	*/

  CbcDynamicPseudoCostBranchingObject *branchingObject = dynamic_cast< CbcDynamicPseudoCostBranchingObject * >(object_);
  if (!branchingObject) {
    delete object_;
    object_ = NULL;
    return;
  }
  CbcSimpleIntegerDynamicPseudoCost *object = branchingObject->object();
  /*
	change is the change in objective due to the branch we've just imposed. It's
	possible we may have gone infeasible.
	*/
  double change = CoinMax(0.0, objectiveValue - originalValue);
  // probably should also ignore if stopped
  // FIXME. Could use enum to avoid numbers for iStatus (e.g. optimal, unknown, infeasible)
  int iStatus;
  if (solver->isProvenOptimal())
    iStatus = 0; // optimal
  else if (solver->isIterationLimitReached()
    && !solver->isDualObjectiveLimitReached())
    iStatus = 2; // unknown
  else
    iStatus = 1; // infeasible
  /*
	  If we're feasible according to the solver, evaluate integer feasibility.
	*/
  bool feasible = iStatus != 1;
  if (feasible) {
    double integerTolerance = model->getDblParam(CbcModel::CbcIntegerTolerance);
    const int *integerVariable = model->integerVariable();
    for (i = 0; i < numberIntegers; i++) {
      int j = integerVariable[i];
      double value = solution[j];
      double nearest = floor(value + 0.5);
      if (fabs(value - nearest) > integerTolerance)
        unsatisfied++;
    }
  }
  /*
	  Finally, update the object. Defaults (080104) are TYPE2 = 0, INFEAS = 1.

	  Pseudocosts are at heart the average of actual costs for a branch. We just
	  need to update the information used to calculate that average.
	*/
  int way = object_->way();
  double value = object_->value();
  //#define TYPE2 1
  //#define TYPERATIO 0.9
  if (way < 0) {
    // down
    if (feasible) {
      double movement = value - floor(value);
      movement = CoinMax(movement, MINIMUM_MOVEMENT);
      //printf("(down change %g value down %g ",change,movement);
      object->incrementNumberTimesDown();
      object->addToSumDownChange(nonZeroAmount + movement);
      object->addToSumDownDecrease(originalUnsatisfied - unsatisfied);
#if TYPE2 == 0
      object->addToSumDownCost(change / (nonZeroAmount + movement));
      object->setDownDynamicPseudoCost(object->sumDownCost() / static_cast< double >(object->numberTimesDown()));
#elif TYPE2 == 1
      object->addToSumDownCost(change);
      object->setDownDynamicPseudoCost(object->sumDownCost() / object->sumDownChange());
#elif TYPE2 == 2
      object->addToSumDownCost(change * TYPERATIO + (1.0 - TYPERATIO) * change / (nonZeroAmount + movement));
      object->setDownDynamicPseudoCost(object->sumDownCost() * (TYPERATIO / object->sumDownChange() + (1.0 - TYPERATIO) / (double)object->numberTimesDown()));
#endif
    } else {
      //printf("(down infeasible value down %g ",change,movement);
      object->incrementNumberTimesDown();
      object->incrementNumberTimesDownInfeasible();
#if INFEAS == 2
      double distanceToCutoff = 0.0;
      double objectiveValue = model->getCurrentMinimizationObjValue();
      distanceToCutoff = model->getCutoff() - originalValue;
      if (distanceToCutoff < 1.0e20)
        change = distanceToCutoff * 2.0;
      else
        change = object->downDynamicPseudoCost() * movement * 10.0;
      change = CoinMax(1.0e-12 * (1.0 + fabs(originalValue)), change);
      object->addToSumDownChange(nonZeroAmount + movement);
      object->addToSumDownDecrease(originalUnsatisfied - unsatisfied);
#if TYPE2 == 0
      object->addToSumDownCost(change / (nonZeroAmount + movement));
      object->setDownDynamicPseudoCost(object->sumDownCost() / (double)object->numberTimesDown());
#elif TYPE2 == 1
      object->addToSumDownCost(change);
      object->setDownDynamicPseudoCost(object->sumDownCost() / object->sumDownChange());
#elif TYPE2 == 2
      object->addToSumDownCost(change * TYPERATIO + (1.0 - TYPERATIO) * change / (nonZeroAmount + movement));
      object->setDownDynamicPseudoCost(object->sumDownCost() * (TYPERATIO / object->sumDownChange() + (1.0 - TYPERATIO) / (double)object->numberTimesDown()));
#endif
#endif
    }
  } else {
    // up
    if (feasible) {
      double movement = ceil(value) - value;
      movement = CoinMax(movement, MINIMUM_MOVEMENT);
      //printf("(up change %g value down %g ",change,movement);
      object->incrementNumberTimesUp();
      object->addToSumUpChange(nonZeroAmount + movement);
      object->addToSumUpDecrease(unsatisfied - originalUnsatisfied);
#if TYPE2 == 0
      object->addToSumUpCost(change / (nonZeroAmount + movement));
      object->setUpDynamicPseudoCost(object->sumUpCost() / static_cast< double >(object->numberTimesUp()));
#elif TYPE2 == 1
      object->addToSumUpCost(change);
      object->setUpDynamicPseudoCost(object->sumUpCost() / object->sumUpChange());
#elif TYPE2 == 2
      object->addToSumUpCost(change * TYPERATIO + (1.0 - TYPERATIO) * change / (nonZeroAmount + movement));
      object->setUpDynamicPseudoCost(object->sumUpCost() * (TYPERATIO / object->sumUpChange() + (1.0 - TYPERATIO) / (double)object->numberTimesUp()));
#endif
    } else {
      //printf("(up infeasible value down %g ",change,movement);
      object->incrementNumberTimesUp();
      object->incrementNumberTimesUpInfeasible();
#if INFEAS == 2
      double distanceToCutoff = 0.0;
      double objectiveValue = model->getCurrentMinimizationObjValue();
      distanceToCutoff = model->getCutoff() - originalValue;
      if (distanceToCutoff < 1.0e20)
        change = distanceToCutoff * 2.0;
      else
        change = object->upDynamicPseudoCost() * movement * 10.0;
      change = CoinMax(1.0e-12 * (1.0 + fabs(originalValue)), change);
      object->addToSumUpChange(nonZeroAmount + movement);
      object->addToSumUpDecrease(unsatisfied - originalUnsatisfied);
#if TYPE2 == 0
      object->addToSumUpCost(change / (nonZeroAmount + movement));
      object->setUpDynamicPseudoCost(object->sumUpCost() / (double)object->numberTimesUp());
#elif TYPE2 == 1
      object->addToSumUpCost(change);
      object->setUpDynamicPseudoCost(object->sumUpCost() / object->sumUpChange());
#elif TYPE2 == 2
      object->addToSumUpCost(change * TYPERATIO + (1.0 - TYPERATIO) * change / (nonZeroAmount + movement));
      object->setUpDynamicPseudoCost(object->sumUpCost() * (TYPERATIO / object->sumUpChange() + (1.0 - TYPERATIO) / (double)object->numberTimesUp()));
#endif
#endif
    }
  }
  //object->print(1,0.5);
  delete object_;
  object_ = NULL;
}

/*
  Simple dynamic decision algorithm. Compare based on infeasibility (numInfUp,
  numInfDown) until a solution is found by search, then switch to change in
  objective (changeUp, changeDown). Note that bestSoFar is remembered in
  bestObject_, so the parameter bestSoFar is unused.
*/

int CbcBranchDynamicDecision::betterBranch(CbcBranchingObject *thisOne,
  CbcBranchingObject * /*bestSoFar*/,
  double changeUp, int numInfUp,
  double changeDown, int numInfDown)
{
  CbcModel *model = thisOne->model();
  int stateOfSearch = model->stateOfSearch() % 10;
  int betterWay = 0;
  double value = 0.0;
  if (!bestObject_) {
    bestCriterion_ = -1.0e30;
    bestNumberUp_ = COIN_INT_MAX;
    bestNumberDown_ = COIN_INT_MAX;
  }
  // maybe branch up more if no solution or not many nodes done?
  if (stateOfSearch <= 2) {
    //#define TRY_STUFF 1
#ifdef TRY_STUFF
    // before solution - choose smallest number
    // could add in depth as well
    int bestNumber = CoinMin(bestNumberUp_, bestNumberDown_);
    if (numInfUp < numInfDown) {
      if (numInfUp < bestNumber) {
        betterWay = 1;
      } else if (numInfUp == bestNumber) {
        if (changeUp < bestChangeUp_)
          betterWay = 1;
      }
    } else if (numInfUp > numInfDown) {
      if (numInfDown < bestNumber) {
        betterWay = -1;
      } else if (numInfDown == bestNumber) {
        if (changeDown < bestChangeDown_)
          betterWay = -1;
      }
    } else {
      // up and down have same number
      bool better = false;
      if (numInfUp < bestNumber) {
        better = true;
      } else if (numInfUp == bestNumber) {
        if (CoinMin(changeUp, changeDown) < CoinMin(bestChangeUp_, bestChangeDown_) - 1.0e-5)
          better = true;
        ;
      }
      if (better) {
        // see which way
        if (changeUp <= changeDown)
          betterWay = 1;
        else
          betterWay = -1;
      }
    }
    if (betterWay) {
      value = CoinMin(numInfUp, numInfDown);
    }
#else
    // use pseudo shadow prices modified by locks
    // test testosi
#ifndef JJF_ONE
    double objectiveValue = model->getCurrentMinimizationObjValue();
    double distanceToCutoff = model->getCutoff() - objectiveValue;
    if (distanceToCutoff < 1.0e20)
      distanceToCutoff *= 10.0;
    else
      distanceToCutoff = 1.0e2 + fabs(objectiveValue);
    distanceToCutoff = CoinMax(distanceToCutoff, 1.0e-12 * (1.0 + fabs(objectiveValue)));
    double continuousObjective = model->getContinuousObjective();
    double distanceToCutoffC = model->getCutoff() - continuousObjective;
    if (distanceToCutoffC > 1.0e20)
      distanceToCutoffC = 1.0e2 + fabs(objectiveValue);
    distanceToCutoffC = CoinMax(distanceToCutoffC, 1.0e-12 * (1.0 + fabs(objectiveValue)));
    int numberInfC = model->getContinuousInfeasibilities();
    double perInf = distanceToCutoffC / static_cast< double >(numberInfC);
    assert(perInf > 0.0);
    //int numberIntegers = model->numberIntegers();
    changeDown += perInf * numInfDown;
    changeUp += perInf * numInfUp;
#ifdef JJF_ZERO
    if (numInfDown == 1) {
      if (numInfUp == 1) {
        changeUp += 1.0e6;
        changeDown += 1.0e6;
      } else if (changeDown <= 1.5 * changeUp) {
        changeUp += 1.0e6;
      }
    } else if (numInfUp == 1 && changeUp <= 1.5 * changeDown) {
      changeDown += 1.0e6;
    }
#endif
#endif
    double minValue = CoinMin(changeDown, changeUp);
    double maxValue = CoinMax(changeDown, changeUp);
    value = WEIGHT_BEFORE * minValue + (1.0 - WEIGHT_BEFORE) * maxValue;
    if (value > bestCriterion_ + 1.0e-8) {
      if (changeUp <= 1.5 * changeDown) {
        betterWay = 1;
      } else {
        betterWay = -1;
      }
    }
#endif
  } else {
#define TRY_STUFF 2
#if TRY_STUFF > 1
    // Get current number of infeasibilities, cutoff and current objective
    CbcNode *node = model->currentNode();
    int numberUnsatisfied = node ? node->numberUnsatisfied() : 1;
    double cutoff = model->getCutoff();
    double objectiveValue = node ? node->objectiveValue() : 0.0;
#endif
    // got a solution
    double minValue = CoinMin(changeDown, changeUp);
    double maxValue = CoinMax(changeDown, changeUp);
    // Reduce
#ifndef WEIGHT_PRODUCT
    value = WEIGHT_AFTER * minValue + (1.0 - WEIGHT_AFTER) * maxValue;
#else
    double minProductWeight = model->getDblParam(CbcModel::CbcSmallChange);
    value = CoinMax(minValue, minProductWeight) * CoinMax(maxValue, minProductWeight);
    //value += minProductWeight*minValue;
#endif
    double useValue = value;
    double useBest = bestCriterion_;
#if TRY_STUFF > 1
    if (node) {
      int thisNumber = CoinMin(numInfUp, numInfDown);
      int bestNumber = CoinMin(bestNumberUp_, bestNumberDown_);
      double distance = cutoff - objectiveValue;
      assert(distance >= 0.0);
      if (useValue + 0.1 * distance > useBest && useValue * 1.1 > useBest && useBest + 0.1 * distance > useValue && useBest * 1.1 > useValue) {
	// not much in it - look at unsatisfied
	if (thisNumber < numberUnsatisfied || bestNumber < numberUnsatisfied) {
	  double perInteger = distance / (static_cast< double >(numberUnsatisfied));
	  useValue += thisNumber * perInteger;
	  useBest += bestNumber * perInteger;
	}
      }
    }
#endif
    if (useValue > useBest + 1.0e-8) {
      if (changeUp <= 1.5 * changeDown) {
        betterWay = 1;
      } else {
        betterWay = -1;
      }
    }
  }
#ifdef COIN_DEVELOP
  History hist;
  {
    CbcDynamicPseudoCostBranchingObject *branchingObject = dynamic_cast< CbcDynamicPseudoCostBranchingObject * >(thisOne);
    if (branchingObject) {
      CbcSimpleIntegerDynamicPseudoCost *object = branchingObject->object();
      assert(object);
      hist.where_ = 'C';
      hist.status_ = ' ';
      hist.sequence_ = object->columnNumber();
      hist.numberUp_ = object->numberTimesUp();
      hist.numberUpInf_ = numInfUp;
      hist.sumUp_ = object->sumUpCost();
      hist.upEst_ = changeUp;
      hist.numberDown_ = object->numberTimesDown();
      hist.numberDownInf_ = numInfDown;
      hist.sumDown_ = object->sumDownCost();
      hist.downEst_ = changeDown;
    }
  }
#endif
  if (betterWay) {
#ifdef COIN_DEVELOP
    hist.status_ = 'B';
#endif
    // maybe change better way
    CbcDynamicPseudoCostBranchingObject *branchingObject = dynamic_cast< CbcDynamicPseudoCostBranchingObject * >(thisOne);
    if (branchingObject) {
      CbcSimpleIntegerDynamicPseudoCost *object = branchingObject->object();
      double separator = object->upDownSeparator();
      if (separator > 0.0) {
        const double *solution = thisOne->model()->testSolution();
        double valueVariable = solution[object->columnNumber()];
        betterWay = (valueVariable - floor(valueVariable) >= separator) ? 1 : -1;
      }
    }
    bestCriterion_ = value;
    bestChangeUp_ = changeUp;
    bestNumberUp_ = numInfUp;
    bestChangeDown_ = changeDown;
    bestNumberDown_ = numInfDown;
    bestObject_ = thisOne;
    // See if user is overriding way
    if (thisOne->object() && thisOne->object()->preferredWay())
      betterWay = thisOne->object()->preferredWay();
  }
#ifdef COIN_DEVELOP
  addRecord(hist);
#endif
  return betterWay;
}
/* Sets or gets best criterion so far */
void CbcBranchDynamicDecision::setBestCriterion(double value)
{
  bestCriterion_ = value;
}
double
CbcBranchDynamicDecision::getBestCriterion() const
{
  return bestCriterion_;
}
#ifdef COIN_DEVELOP
void printHistory(const char *file)
{
  if (!numberHistory)
    return;
  FILE *fp = fopen(file, "w");
  assert(fp);
  int numberIntegers = 0;
  int i;
  for (i = 0; i < numberHistory; i++) {
    if (history[i].where_ != 'C' || history[i].status_ != 'I')
      numberIntegers = CoinMax(numberIntegers, history[i].sequence_);
  }
  numberIntegers++;
  for (int iC = 0; iC < numberIntegers; iC++) {
    int n = 0;
    for (i = 0; i < numberHistory; i++) {
      if (history[i].sequence_ == iC) {
        if (!n)
          fprintf(fp, "XXX %d\n", iC);
        n++;
        fprintf(fp, "%c%c up %8d %8d %12.5f %12.5f down  %8d %8d %12.5f %12.5f\n",
          history[i].where_,
          history[i].status_,
          history[i].numberUp_,
          history[i].numberUpInf_,
          history[i].sumUp_,
          history[i].upEst_,
          history[i].numberDown_,
          history[i].numberDownInf_,
          history[i].sumDown_,
          history[i].downEst_);
      }
    }
  }
  fclose(fp);
}
#endif

// Default Constructor
CbcDynamicPseudoCostBranchingObject::CbcDynamicPseudoCostBranchingObject()
  : CbcIntegerBranchingObject()
{
  changeInGuessed_ = 1.0e-5;
  object_ = NULL;
}

// Useful constructor
CbcDynamicPseudoCostBranchingObject::CbcDynamicPseudoCostBranchingObject(CbcModel *model,
  int variable,
  int way, double value,
  CbcSimpleIntegerDynamicPseudoCost *object)
  : CbcIntegerBranchingObject(model, variable, way, value)
{
  changeInGuessed_ = 1.0e-5;
  object_ = object;
}
// Does part of work for constructor
void CbcDynamicPseudoCostBranchingObject::fillPart(int variable,
  int way, double value,
  CbcSimpleIntegerDynamicPseudoCost *object)
{
  CbcIntegerBranchingObject::fillPart(variable, way, value);
  changeInGuessed_ = 1.0e-5;
  object_ = object;
}
// Useful constructor for fixing
CbcDynamicPseudoCostBranchingObject::CbcDynamicPseudoCostBranchingObject(CbcModel *model,
  int variable, int way,
  double lowerValue,
  double /*upperValue*/)
  : CbcIntegerBranchingObject(model, variable, way, lowerValue)
{
  changeInGuessed_ = 1.0e100;
  object_ = NULL;
}

// Copy constructor
CbcDynamicPseudoCostBranchingObject::CbcDynamicPseudoCostBranchingObject(
  const CbcDynamicPseudoCostBranchingObject &rhs)
  : CbcIntegerBranchingObject(rhs)
{
  changeInGuessed_ = rhs.changeInGuessed_;
  object_ = rhs.object_;
}

// Assignment operator
CbcDynamicPseudoCostBranchingObject &
CbcDynamicPseudoCostBranchingObject::operator=(const CbcDynamicPseudoCostBranchingObject &rhs)
{
  if (this != &rhs) {
    CbcIntegerBranchingObject::operator=(rhs);
    changeInGuessed_ = rhs.changeInGuessed_;
    object_ = rhs.object_;
  }
  return *this;
}
CbcBranchingObject *
CbcDynamicPseudoCostBranchingObject::clone() const
{
  return (new CbcDynamicPseudoCostBranchingObject(*this));
}

// Destructor
CbcDynamicPseudoCostBranchingObject::~CbcDynamicPseudoCostBranchingObject()
{
}

/*
  Perform a branch by adjusting the bounds of the specified variable. Note
  that each arm of the branch advances the object to the next arm by
  advancing the value of way_.

  Providing new values for the variable's lower and upper bounds for each
  branching direction gives a little bit of additional flexibility and will
  be easily extensible to multi-way branching.
  Returns change in guessed objective on next branch
*/
double
CbcDynamicPseudoCostBranchingObject::branch()
{
  CbcIntegerBranchingObject::branch();
  return changeInGuessed_;
}
/* Some branchingObjects may claim to be able to skip
   strong branching.  If so they have to fill in CbcStrongInfo.
   The object mention in incoming CbcStrongInfo must match.
   Returns nonzero if skip is wanted */
int CbcDynamicPseudoCostBranchingObject::fillStrongInfo(CbcStrongInfo &info)
{
  assert(object_);
  assert(info.possibleBranch == this);
  info.upMovement = object_->upDynamicPseudoCost() * (ceil(value_) - value_);
  info.downMovement = object_->downDynamicPseudoCost() * (value_ - floor(value_));
  info.numIntInfeasUp -= static_cast< int >(object_->sumUpDecrease() / (1.0e-12 + static_cast< double >(object_->numberTimesUp())));
  info.numIntInfeasUp = CoinMax(info.numIntInfeasUp, 0);
  info.numObjInfeasUp = 0;
  info.finishedUp = false;
  info.numItersUp = 0;
  info.numIntInfeasDown -= static_cast< int >(object_->sumDownDecrease() / (1.0e-12 + static_cast< double >(object_->numberTimesDown())));
  info.numIntInfeasDown = CoinMax(info.numIntInfeasDown, 0);
  info.numObjInfeasDown = 0;
  info.finishedDown = false;
  info.numItersDown = 0;
  info.fix = 0;
  if (object_->numberTimesUp() < object_->numberBeforeTrust() + 2 * object_->numberTimesUpInfeasible() || object_->numberTimesDown() < object_->numberBeforeTrust() + 2 * object_->numberTimesDownInfeasible()) {
    return 0;
  } else {
    return 1;
  }
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
