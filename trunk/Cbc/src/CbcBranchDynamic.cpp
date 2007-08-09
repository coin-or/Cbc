// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchDynamic.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
CbcSimpleIntegerDynamicPseudoCost::CbcSimpleIntegerDynamicPseudoCost ()
  : CbcSimpleInteger(),
    downDynamicPseudoCost_(1.0e-5),
    upDynamicPseudoCost_(1.0e-5),
    upDownSeparator_(-1.0),
    sumDownCost_(0.0),
    sumUpCost_(0.0),
    sumDownChange_(0.0),
    sumUpChange_(0.0),
    sumDownCostSquared_(0.0),
    sumUpCostSquared_(0.0),
    sumDownDecrease_(0.0),
    sumUpDecrease_(0.0),
    lastDownCost_(0.0),
    lastUpCost_(0.0),
    lastDownDecrease_(0),
    lastUpDecrease_(0),
    numberTimesDown_(0),
    numberTimesUp_(0),
    numberTimesDownInfeasible_(0),
    numberTimesUpInfeasible_(0),
    numberBeforeTrust_(0),
    numberTimesDownLocalFixed_(0),
    numberTimesUpLocalFixed_(0),
    numberTimesDownTotalFixed_(0.0),
    numberTimesUpTotalFixed_(0.0),
    numberTimesProbingTotal_(0),
    method_(0)
{
}

/** Useful constructor

  Loads dynamic upper & lower bounds for the specified variable.
*/
CbcSimpleIntegerDynamicPseudoCost::CbcSimpleIntegerDynamicPseudoCost (CbcModel * model,
				    int iColumn, double breakEven)
  : CbcSimpleInteger(model,iColumn,breakEven),
    upDownSeparator_(-1.0),
    sumDownCost_(0.0),
    sumUpCost_(0.0),
    sumDownChange_(0.0),
    sumUpChange_(0.0),
    sumDownCostSquared_(0.0),
    sumUpCostSquared_(0.0),
    sumDownDecrease_(0.0),
    sumUpDecrease_(0.0),
    lastDownCost_(0.0),
    lastUpCost_(0.0),
    lastDownDecrease_(0),
    lastUpDecrease_(0),
    numberTimesDown_(0),
    numberTimesUp_(0),
    numberTimesDownInfeasible_(0),
    numberTimesUpInfeasible_(0),
    numberBeforeTrust_(0),
    numberTimesDownLocalFixed_(0),
    numberTimesUpLocalFixed_(0),
    numberTimesDownTotalFixed_(0.0),
    numberTimesUpTotalFixed_(0.0),
    numberTimesProbingTotal_(0),
    method_(0)
{
  const double * cost = model->getObjCoefficients();
  double costValue = CoinMax(1.0e-5,fabs(cost[iColumn]));
  // treat as if will cost what it says up
  upDynamicPseudoCost_=costValue;
  // and balance at breakeven
  downDynamicPseudoCost_=((1.0-breakEven_)*upDynamicPseudoCost_)/breakEven_;
  // so initial will have some effect
  sumUpCost_ = 2.0*upDynamicPseudoCost_;
  sumUpChange_ = 2.0;
  numberTimesUp_ = 2;
  sumDownCost_ = 2.0*downDynamicPseudoCost_;
  sumDownChange_ = 2.0;
  numberTimesDown_ = 2;
#if 0
  // No
  sumUpCost_ = 0.0;
  sumUpChange_ = 0.0;
  numberTimesUp_ = 0;
  sumDownCost_ = 0.0;
  sumDownChange_ = 0.0;
  numberTimesDown_ = 0;
#else
  sumUpCost_ = 1.0*upDynamicPseudoCost_;
  sumUpChange_ = 1.0;
  numberTimesUp_ = 1;
  sumDownCost_ = 1.0*downDynamicPseudoCost_;
  sumDownChange_ = 1.0;
  numberTimesDown_ = 1;
#endif
}

/** Useful constructor

  Loads dynamic upper & lower bounds for the specified variable.
*/
CbcSimpleIntegerDynamicPseudoCost::CbcSimpleIntegerDynamicPseudoCost (CbcModel * model,
				    int iColumn, double downDynamicPseudoCost,
							double upDynamicPseudoCost)
  : CbcSimpleInteger(model,iColumn),
    upDownSeparator_(-1.0),
    sumDownCost_(0.0),
    sumUpCost_(0.0),
    sumDownChange_(0.0),
    sumUpChange_(0.0),
    sumDownCostSquared_(0.0),
    sumUpCostSquared_(0.0),
    sumDownDecrease_(0.0),
    sumUpDecrease_(0.0),
    lastDownCost_(0.0),
    lastUpCost_(0.0),
    lastDownDecrease_(0),
    lastUpDecrease_(0),
    numberTimesDown_(0),
    numberTimesUp_(0),
    numberTimesDownInfeasible_(0),
    numberTimesUpInfeasible_(0),
    numberBeforeTrust_(0),
    numberTimesDownLocalFixed_(0),
    numberTimesUpLocalFixed_(0),
    numberTimesDownTotalFixed_(0.0),
    numberTimesUpTotalFixed_(0.0),
    numberTimesProbingTotal_(0),
    method_(0)
{
  downDynamicPseudoCost_ = downDynamicPseudoCost;
  upDynamicPseudoCost_ = upDynamicPseudoCost;
  breakEven_ = upDynamicPseudoCost_/(upDynamicPseudoCost_+downDynamicPseudoCost_);
  // so initial will have some effect
  sumUpCost_ = 2.0*upDynamicPseudoCost_;
  sumUpChange_ = 2.0;
  numberTimesUp_ = 2;
  sumDownCost_ = 2.0*downDynamicPseudoCost_;
  sumDownChange_ = 2.0;
  numberTimesDown_ = 2;
#if 0
  // No
  sumUpCost_ = 0.0;
  sumUpChange_ = 0.0;
  numberTimesUp_ = 0;
  sumDownCost_ = 0.0;
  sumDownChange_ = 0.0;
  numberTimesDown_ = 0;
#else
  sumUpCost_ = 1.0*upDynamicPseudoCost_;
  sumUpChange_ = 1.0;
  numberTimesUp_ = 1;
  sumDownCost_ = 1.0*downDynamicPseudoCost_;
  sumDownChange_ = 1.0;
  numberTimesDown_ = 1;
#endif
}
/** Useful constructor

  Loads dynamic upper & lower bounds for the specified variable.
*/
CbcSimpleIntegerDynamicPseudoCost::CbcSimpleIntegerDynamicPseudoCost (CbcModel * model,
				    int dummy, int iColumn, double downDynamicPseudoCost,
							double upDynamicPseudoCost)
{
  CbcSimpleIntegerDynamicPseudoCost(model,iColumn,downDynamicPseudoCost,upDynamicPseudoCost);
}

// Copy constructor 
CbcSimpleIntegerDynamicPseudoCost::CbcSimpleIntegerDynamicPseudoCost ( const CbcSimpleIntegerDynamicPseudoCost & rhs)
  :CbcSimpleInteger(rhs),
   downDynamicPseudoCost_(rhs.downDynamicPseudoCost_),
   upDynamicPseudoCost_(rhs.upDynamicPseudoCost_),
   upDownSeparator_(rhs.upDownSeparator_),
   sumDownCost_(rhs.sumDownCost_),
   sumUpCost_(rhs.sumUpCost_),
   sumDownChange_(rhs.sumDownChange_),
   sumUpChange_(rhs.sumUpChange_),
   sumDownCostSquared_(rhs.sumDownCostSquared_),
   sumUpCostSquared_(rhs.sumUpCostSquared_),
   sumDownDecrease_(rhs.sumDownDecrease_),
   sumUpDecrease_(rhs.sumUpDecrease_),
   lastDownCost_(rhs.lastDownCost_),
   lastUpCost_(rhs.lastUpCost_),
   lastDownDecrease_(rhs.lastDownDecrease_),
   lastUpDecrease_(rhs.lastUpDecrease_),
   numberTimesDown_(rhs.numberTimesDown_),
   numberTimesUp_(rhs.numberTimesUp_),
   numberTimesDownInfeasible_(rhs.numberTimesDownInfeasible_),
   numberTimesUpInfeasible_(rhs.numberTimesUpInfeasible_),
   numberBeforeTrust_(rhs.numberBeforeTrust_),
   numberTimesDownLocalFixed_(rhs.numberTimesDownLocalFixed_),
   numberTimesUpLocalFixed_(rhs.numberTimesUpLocalFixed_),
   numberTimesDownTotalFixed_(rhs.numberTimesDownTotalFixed_),
   numberTimesUpTotalFixed_(rhs.numberTimesUpTotalFixed_),
   numberTimesProbingTotal_(rhs.numberTimesProbingTotal_),
   method_(rhs.method_)

{
}

// Clone
CbcObject *
CbcSimpleIntegerDynamicPseudoCost::clone() const
{
  return new CbcSimpleIntegerDynamicPseudoCost(*this);
}

// Assignment operator 
CbcSimpleIntegerDynamicPseudoCost & 
CbcSimpleIntegerDynamicPseudoCost::operator=( const CbcSimpleIntegerDynamicPseudoCost& rhs)
{
  if (this!=&rhs) {
    CbcSimpleInteger::operator=(rhs);
    downDynamicPseudoCost_=rhs.downDynamicPseudoCost_;
    upDynamicPseudoCost_=rhs.upDynamicPseudoCost_;
    upDownSeparator_=rhs.upDownSeparator_;
    sumDownCost_ = rhs.sumDownCost_;
    sumUpCost_ = rhs.sumUpCost_;
    sumDownChange_ = rhs.sumDownChange_;
    sumUpChange_ = rhs.sumUpChange_;
    sumDownCostSquared_ = rhs.sumDownCostSquared_;
    sumUpCostSquared_ = rhs.sumUpCostSquared_;
    sumDownDecrease_ = rhs.sumDownDecrease_;
    sumUpDecrease_ = rhs.sumUpDecrease_;
    lastDownCost_ = rhs.lastDownCost_;
    lastUpCost_ = rhs.lastUpCost_;
    lastDownDecrease_ = rhs.lastDownDecrease_;
    lastUpDecrease_ = rhs.lastUpDecrease_;
    numberTimesDown_ = rhs.numberTimesDown_;
    numberTimesUp_ = rhs.numberTimesUp_;
    numberTimesDownInfeasible_ = rhs.numberTimesDownInfeasible_;
    numberTimesUpInfeasible_ = rhs.numberTimesUpInfeasible_;
    numberBeforeTrust_ = rhs.numberBeforeTrust_;
    numberTimesDownLocalFixed_ = rhs.numberTimesDownLocalFixed_;
    numberTimesUpLocalFixed_ = rhs.numberTimesUpLocalFixed_;
    numberTimesDownTotalFixed_ = rhs.numberTimesDownTotalFixed_;
    numberTimesUpTotalFixed_ = rhs.numberTimesUpTotalFixed_;
    numberTimesProbingTotal_ = rhs.numberTimesProbingTotal_;
    method_=rhs.method_;
  }
  return *this;
}

// Destructor 
CbcSimpleIntegerDynamicPseudoCost::~CbcSimpleIntegerDynamicPseudoCost ()
{
}
// Copy some information i.e. just variable stuff
void 
CbcSimpleIntegerDynamicPseudoCost::copySome(CbcSimpleIntegerDynamicPseudoCost * otherObject)
{
  downDynamicPseudoCost_=otherObject->downDynamicPseudoCost_;
  upDynamicPseudoCost_=otherObject->upDynamicPseudoCost_;
  sumDownCost_ = otherObject->sumDownCost_;
  sumUpCost_ = otherObject->sumUpCost_;
  sumDownChange_ = otherObject->sumDownChange_;
  sumUpChange_ = otherObject->sumUpChange_;
  sumDownCostSquared_ = otherObject->sumDownCostSquared_;
  sumUpCostSquared_ = otherObject->sumUpCostSquared_;
  sumDownDecrease_ = otherObject->sumDownDecrease_;
  sumUpDecrease_ = otherObject->sumUpDecrease_;
  lastDownCost_ = otherObject->lastDownCost_;
  lastUpCost_ = otherObject->lastUpCost_;
  lastDownDecrease_ = otherObject->lastDownDecrease_;
  lastUpDecrease_ = otherObject->lastUpDecrease_;
  numberTimesDown_ = otherObject->numberTimesDown_;
  numberTimesUp_ = otherObject->numberTimesUp_;
  numberTimesDownInfeasible_ = otherObject->numberTimesDownInfeasible_;
  numberTimesUpInfeasible_ = otherObject->numberTimesUpInfeasible_;
  numberTimesDownLocalFixed_ = otherObject->numberTimesDownLocalFixed_;
  numberTimesUpLocalFixed_ = otherObject->numberTimesUpLocalFixed_;
  numberTimesDownTotalFixed_ = otherObject->numberTimesDownTotalFixed_;
  numberTimesUpTotalFixed_ = otherObject->numberTimesUpTotalFixed_;
  numberTimesProbingTotal_ = otherObject->numberTimesProbingTotal_;
}
// Creates a branching object
CbcBranchingObject * 
CbcSimpleIntegerDynamicPseudoCost::createBranch(int way) 
{
  const double * solution = model_->testSolution();
  const double * lower = model_->getCbcColLower();
  const double * upper = model_->getCbcColUpper();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
#ifndef NDEBUG
  double nearest = floor(value+0.5);
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  assert (upper[columnNumber_]>lower[columnNumber_]);
#endif
  if (!model_->hotstartSolution()) {
    assert (fabs(value-nearest)>integerTolerance);
  } else {
    const double * hotstartSolution = model_->hotstartSolution();
    double targetValue = hotstartSolution[columnNumber_];
    if (way>0)
      value = targetValue-0.1;
    else
      value = targetValue+0.1;
  }
  CbcDynamicPseudoCostBranchingObject * newObject = 
    new CbcDynamicPseudoCostBranchingObject(model_,columnNumber_,way,
					    value,this);
  double up =  upDynamicPseudoCost_*(ceil(value)-value);
  double down =  downDynamicPseudoCost_*(value-floor(value));
  double changeInGuessed=up-down;
  if (way>0)
    changeInGuessed = - changeInGuessed;
  changeInGuessed=CoinMax(0.0,changeInGuessed);
  //if (way>0)
  //changeInGuessed += 1.0e8; // bias to stay up
  newObject->setChangeInGuessed(changeInGuessed);
  newObject->setOriginalObject(this);
  return newObject;
}
/* Create an OsiSolverBranch object
   
This returns NULL if branch not represented by bound changes
*/
OsiSolverBranch * 
CbcSimpleIntegerDynamicPseudoCost::solverBranch() const
{
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  assert (upper[columnNumber_]>lower[columnNumber_]);
#ifndef NDEBUG
  double nearest = floor(value+0.5);
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  assert (fabs(value-nearest)>integerTolerance);
#endif
  OsiSolverBranch * branch = new OsiSolverBranch();
  branch->addBranch(columnNumber_,value);
  return branch;
}
  
// Infeasibility - large is 0.5
double 
CbcSimpleIntegerDynamicPseudoCost::infeasibility(int & preferredWay) const
{
  const double * solution = model_->testSolution();
  const double * lower = model_->getCbcColLower();
  const double * upper = model_->getCbcColUpper();
  if (upper[columnNumber_]==lower[columnNumber_]) {
    // fixed
    preferredWay=1;
    return 0.0;
  }
  assert (breakEven_>0.0&&breakEven_<1.0);
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  /*printf("%d %g %g %g %g\n",columnNumber_,value,lower[columnNumber_],
    solution[columnNumber_],upper[columnNumber_]);*/
  double nearest = floor(value+0.5);
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double below = floor(value+integerTolerance);
  double above = below+1.0;
  if (above>upper[columnNumber_]) {
    above=below;
    below = above -1;
  }
#define INFEAS 1
#if INFEAS==1
  double distanceToCutoff=0.0;
  double objectiveValue = model_->getCurrentMinimizationObjValue();
  distanceToCutoff =  model_->getCutoff()  - objectiveValue;
  if (distanceToCutoff<1.0e20) 
    distanceToCutoff *= 10.0;
  else 
    distanceToCutoff = 1.0e2 + fabs(objectiveValue);
#endif
  double sum;
  int number;
  double downCost = CoinMax(value-below,0.0);
  sum = sumDownCost_;
  number = numberTimesDown_;
#if INFEAS==1
  sum += numberTimesDownInfeasible_*(distanceToCutoff/(downCost+1.0e-12));
  number += numberTimesDownInfeasible_;
#endif
  if (number>0)
    downCost *= sum / (double) number;
  else
    downCost  *=  downDynamicPseudoCost_;
  double upCost = CoinMax((above-value),0.0);
  sum = sumUpCost_;
  number = numberTimesUp_;
#if INFEAS==1
  sum += numberTimesUpInfeasible_*(distanceToCutoff/(upCost+1.0e-12));
  number += numberTimesUpInfeasible_;
#endif
  if (number>0)
    upCost *= sum / (double) number;
  else
    upCost  *=  upDynamicPseudoCost_;
  if (downCost>=upCost)
    preferredWay=1;
  else
    preferredWay=-1;
  // See if up down choice set
  if (upDownSeparator_>0.0) {
    preferredWay = (value-below>=upDownSeparator_) ? 1 : -1;
  }
  if (preferredWay_)
    preferredWay=preferredWay_;
  // weight at 1.0 is max min
#define WEIGHT_AFTER 0.7
#define WEIGHT_BEFORE 0.1
  if (fabs(value-nearest)<=integerTolerance) {
    return 0.0;
  } else {
    int stateOfSearch = model_->stateOfSearch();
    double returnValue=0.0;
    double minValue = CoinMin(downCost,upCost);
    double maxValue = CoinMax(downCost,upCost);
    if (stateOfSearch<=1||model_->currentNode()->depth()<=10/* was ||maxValue>0.2*distanceToCutoff*/) {
      // no solution
      returnValue = WEIGHT_BEFORE*minValue + (1.0-WEIGHT_BEFORE)*maxValue;
      if (0) {
	double sum;
	int number;
	double downCost2 = CoinMax(value-below,0.0);
	sum = sumDownCost_;
	number = numberTimesDown_;
	if (number>0)
	  downCost2 *= sum / (double) number;
	else
	  downCost2  *=  downDynamicPseudoCost_;
	double upCost2 = CoinMax((above-value),0.0);
	sum = sumUpCost_;
	number = numberTimesUp_;
	if (number>0)
	  upCost2 *= sum / (double) number;
	else
	  upCost2  *=  upDynamicPseudoCost_;
	double minValue2 = CoinMin(downCost2,upCost2);
	double maxValue2 = CoinMax(downCost2,upCost2);
	printf("%d value %g downC %g upC %g minV %g maxV %g downC2 %g upC2 %g minV2 %g maxV2 %g\n",
	       columnNumber_,value,downCost,upCost,minValue,maxValue,
	       downCost2,upCost2,minValue2,maxValue2);
      }
    } else {
      // some solution
      returnValue = WEIGHT_AFTER*minValue + (1.0-WEIGHT_AFTER)*maxValue;
    }
    if (numberTimesUp_<numberBeforeTrust_||
        numberTimesDown_<numberBeforeTrust_) {
      //if (returnValue<1.0e10)
      //returnValue += 1.0e12;
      //else
      returnValue *= 1.0e3;
      if (!numberTimesUp_&&!numberTimesDown_)
        returnValue=1.0e50;
    }
    //if (fabs(value-0.5)<1.0e-5) {
    //returnValue = 3.0*returnValue + 0.2;
    //} else if (value>0.9) {
    //returnValue = 2.0*returnValue + 0.1;
    //}
    if (method_==1) {
      // probing
      // average 
      double up=1.0e-15;
      double down=1.0e-15;
      if (numberTimesProbingTotal_) {
	up += numberTimesUpTotalFixed_/((double) numberTimesProbingTotal_);
	down += numberTimesDownTotalFixed_/((double) numberTimesProbingTotal_);
      }
      returnValue = 1 + 10.0*CoinMin(numberTimesDownLocalFixed_,numberTimesUpLocalFixed_) +
	CoinMin(down,up);
      returnValue *= 1.0e-3;
    }
    return CoinMax(returnValue,1.0e-15);
  }
}

double
CbcSimpleIntegerDynamicPseudoCost::infeasibility(const OsiSolverInterface * solver, const OsiBranchingInformation * info,
			 int & preferredWay) const
{
  double value = info->solution_[columnNumber_];
  value = CoinMax(value, info->lower_[columnNumber_]);
  value = CoinMin(value, info->upper_[columnNumber_]);
  if (info->upper_[columnNumber_]==info->lower_[columnNumber_]) {
    // fixed
    preferredWay=1;
    return 0.0;
  }
  assert (breakEven_>0.0&&breakEven_<1.0);
  double nearest = floor(value+0.5);
  double integerTolerance = info->integerTolerance_; 
  double below = floor(value+integerTolerance);
  double above = below+1.0;
  if (above>info->upper_[columnNumber_]) {
    above=below;
    below = above -1;
  }
#if INFEAS==1
  double objectiveValue = info->objectiveValue_;
  double distanceToCutoff =  info->cutoff_  - objectiveValue;
  if (distanceToCutoff<1.0e20) 
    distanceToCutoff *= 10.0;
  else 
    distanceToCutoff = 1.0e2 + fabs(objectiveValue);
#endif
  double sum;
  int number;
  double downCost = CoinMax(value-below,0.0);
  sum = sumDownCost_;
  number = numberTimesDown_;
#if INFEAS==1
  sum += numberTimesDownInfeasible_*(distanceToCutoff/(downCost+1.0e-12));
  number += numberTimesDownInfeasible_;
#endif
  if (number>0)
    downCost *= sum / (double) number;
  else
    downCost  *=  downDynamicPseudoCost_;
  double upCost = CoinMax((above-value),0.0);
  sum = sumUpCost_;
  number = numberTimesUp_;
#if INFEAS==1
  sum += numberTimesUpInfeasible_*(distanceToCutoff/(upCost+1.0e-12));
  number += numberTimesUpInfeasible_;
#endif
  if (number>0)
    upCost *= sum / (double) number;
  else
    upCost  *=  upDynamicPseudoCost_;
  if (downCost>=upCost)
    preferredWay=1;
  else
    preferredWay=-1;
  // See if up down choice set
  if (upDownSeparator_>0.0) {
    preferredWay = (value-below>=upDownSeparator_) ? 1 : -1;
  }
  if (preferredWay_)
    preferredWay=preferredWay_;
  // weight at 1.0 is max min
  if (fabs(value-nearest)<=integerTolerance) {
    return 0.0;
  } else {
    double returnValue=0.0;
    double minValue = CoinMin(downCost,upCost);
    double maxValue = CoinMax(downCost,upCost);
    if (!info->numberBranchingSolutions_||info->depth_<=10/* was ||maxValue>0.2*distanceToCutoff*/) {
      // no solution
      returnValue = WEIGHT_BEFORE*minValue + (1.0-WEIGHT_BEFORE)*maxValue;
    } else {
      // some solution
      returnValue = WEIGHT_AFTER*minValue + (1.0-WEIGHT_AFTER)*maxValue;
    }
    if (numberTimesUp_<numberBeforeTrust_||
        numberTimesDown_<numberBeforeTrust_) {
      //if (returnValue<1.0e10)
      //returnValue += 1.0e12;
      //else
      returnValue *= 1.0e3;
      if (!numberTimesUp_&&!numberTimesDown_)
        returnValue=1.0e50;
    }
    //if (fabs(value-0.5)<1.0e-5) {
    //returnValue = 3.0*returnValue + 0.2;
    //} else if (value>0.9) {
    //returnValue = 2.0*returnValue + 0.1;
    //}
    if (method_==1) {
      // probing
      // average 
      double up=1.0e-15;
      double down=1.0e-15;
      if (numberTimesProbingTotal_) {
	up += numberTimesUpTotalFixed_/((double) numberTimesProbingTotal_);
	down += numberTimesDownTotalFixed_/((double) numberTimesProbingTotal_);
      }
      returnValue = 1 + 10.0*CoinMin(numberTimesDownLocalFixed_,numberTimesUpLocalFixed_) +
	CoinMin(down,up);
      returnValue *= 1.0e-3;
    }
    return CoinMax(returnValue,1.0e-15);
  }
}
// Creates a branching object
CbcBranchingObject * 
CbcSimpleIntegerDynamicPseudoCost::createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) 
{
  double value = info->solution_[columnNumber_];
  value = CoinMax(value, info->lower_[columnNumber_]);
  value = CoinMin(value, info->upper_[columnNumber_]);
  assert (info->upper_[columnNumber_]>info->lower_[columnNumber_]);
  if (!info->hotstartSolution_) {
#ifndef NDEBUG
    double nearest = floor(value+0.5);
    assert (fabs(value-nearest)>info->integerTolerance_);
#endif
  } else {
    double targetValue = info->hotstartSolution_[columnNumber_];
    if (way>0)
      value = targetValue-0.1;
    else
      value = targetValue+0.1;
  }
  CbcDynamicPseudoCostBranchingObject * newObject = 
    new CbcDynamicPseudoCostBranchingObject(model_,columnNumber_,way,
					    value,this);
  double up =  upDynamicPseudoCost_*(ceil(value)-value);
  double down =  downDynamicPseudoCost_*(value-floor(value));
  double changeInGuessed=up-down;
  if (way>0)
    changeInGuessed = - changeInGuessed;
  changeInGuessed=CoinMax(0.0,changeInGuessed);
  //if (way>0)
  //changeInGuessed += 1.0e8; // bias to stay up
  newObject->setChangeInGuessed(changeInGuessed);
  newObject->setOriginalObject(this);
  return newObject;
}

// Return "up" estimate
double 
CbcSimpleIntegerDynamicPseudoCost::upEstimate() const
{
  const double * solution = model_->testSolution();
  const double * lower = model_->getCbcColLower();
  const double * upper = model_->getCbcColUpper();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  if (upper[columnNumber_]==lower[columnNumber_]) {
    // fixed
    return 0.0;
  }
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double below = floor(value+integerTolerance);
  double above = below+1.0;
  if (above>upper[columnNumber_]) {
    above=below;
    below = above -1;
  }
  double upCost = CoinMax((above-value)*upDynamicPseudoCost_,0.0);
  return upCost;
}
// Return "down" estimate
double 
CbcSimpleIntegerDynamicPseudoCost::downEstimate() const
{
  const double * solution = model_->testSolution();
  const double * lower = model_->getCbcColLower();
  const double * upper = model_->getCbcColUpper();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  if (upper[columnNumber_]==lower[columnNumber_]) {
    // fixed
    return 0.0;
  }
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double below = floor(value+integerTolerance);
  double above = below+1.0;
  if (above>upper[columnNumber_]) {
    above=below;
    below = above -1;
  }
  double downCost = CoinMax((value-below)*downDynamicPseudoCost_,0.0);
  return downCost;
}
/* Pass in information on branch just done and create CbcObjectUpdateData instance.
   If object does not need data then backward pointer will be NULL.
   Assumes can get information from solver */
CbcObjectUpdateData 
CbcSimpleIntegerDynamicPseudoCost::createUpdateInformation(const OsiSolverInterface * solver, 
							   const CbcNode * node,
							   const CbcBranchingObject * branchingObject)
{
  double originalValue=node->objectiveValue();
  int originalUnsatisfied = node->numberUnsatisfied();
  double objectiveValue = solver->getObjValue()*solver->getObjSense();
  int unsatisfied=0;
  int i;
  //might be base model - doesn't matter
  int numberIntegers = model_->numberIntegers();;
  const double * solution = solver->getColSolution();
  //const double * lower = solver->getColLower();
  //const double * upper = solver->getColUpper();
  double change = CoinMax(0.0,objectiveValue-originalValue);
  int iStatus;
  if (solver->isProvenOptimal())
    iStatus=0; // optimal
  else if (solver->isIterationLimitReached()
           &&!solver->isDualObjectiveLimitReached())
    iStatus=2; // unknown 
  else
    iStatus=1; // infeasible

  bool feasible = iStatus!=1;
  if (feasible) {
    double integerTolerance = 
      model_->getDblParam(CbcModel::CbcIntegerTolerance);
    const int * integerVariable = model_->integerVariable();
    for (i=0;i<numberIntegers;i++) {
      int j=integerVariable[i];
      double value = solution[j];
      double nearest = floor(value+0.5);
      if (fabs(value-nearest)>integerTolerance) 
        unsatisfied++;
    }
  }
  int way = branchingObject->way();
  way = - way; // because after branch so moved on
  double value = branchingObject->value();
  return CbcObjectUpdateData (this, way,
				  change, iStatus,
				  originalUnsatisfied-unsatisfied,value);
}
// Update object by CbcObjectUpdateData
void 
CbcSimpleIntegerDynamicPseudoCost::updateInformation(const CbcObjectUpdateData & data)
{
  bool feasible = data.status_!=1;
  int way = data.way_;
  double value = data.branchingValue_;
  double change = data.change_;
#define TYPE2 1
#define TYPERATIO 0.9
  if (way<0) {
    // down
    if (feasible) {
      //printf("(down change %g value down %g ",change,value-floor(value));
      incrementNumberTimesDown();
      addToSumDownChange(1.0e-30+value-floor(value));
      addToSumDownDecrease(data.intDecrease_);
#if TYPE2==0
      addToSumDownCost(change/(1.0e-30+(value-floor(value))));
      setDownDynamicPseudoCost(sumDownCost()/(double) numberTimesDown());
#elif TYPE2==1
      addToSumDownCost(change);
      setDownDynamicPseudoCost(sumDownCost()/sumDownChange());
#elif TYPE2==2
      addToSumDownCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+(value-floor(value))));
      setDownDynamicPseudoCost(sumDownCost()*(TYPERATIO/sumDownChange()+(1.0-TYPERATIO)/(double) numberTimesDown()));
#endif
    } else {
      //printf("(down infeasible value down %g ",change,value-floor(value));
      incrementNumberTimesDownInfeasible();
#if INFEAS==2
      double distanceToCutoff=0.0;
      double objectiveValue = model->getCurrentMinimizationObjValue();
      distanceToCutoff =  model->getCutoff()  - originalValue;
      if (distanceToCutoff<1.0e20) 
	change = distanceToCutoff*2.0;
      else 
	change = downDynamicPseudoCost()*(value-floor(value))*10.0; 
      change = CoinMax(1.0e-12*(1.0+fabs(originalValue)),change);
      incrementNumberTimesDown();
      addToSumDownChange(1.0e-30+value-floor(value));
      addToSumDownDecrease(data.intDecrease_);
#if TYPE2==0
      addToSumDownCost(change/(1.0e-30+(value-floor(value))));
      setDownDynamicPseudoCost(sumDownCost()/(double) numberTimesDown());
#elif TYPE2==1
      addToSumDownCost(change);
      setDownDynamicPseudoCost(sumDownCost()/sumDownChange());
#elif TYPE2==2
      addToSumDownCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+(value-floor(value))));
      setDownDynamicPseudoCost(sumDownCost()*(TYPERATIO/sumDownChange()+(1.0-TYPERATIO)/(double) numberTimesDown()));
#endif
#endif
    }
  } else {
    // up
    if (feasible) {
      //printf("(up change %g value down %g ",change,ceil(value)-value);
      incrementNumberTimesUp();
      addToSumUpChange(1.0e-30+ceil(value)-value);
      addToSumUpDecrease(data.intDecrease_);
#if TYPE2==0
      addToSumUpCost(change/(1.0e-30+(ceil(value)-value)));
      setUpDynamicPseudoCost(sumUpCost()/(double) numberTimesUp());
#elif TYPE2==1
      addToSumUpCost(change);
      setUpDynamicPseudoCost(sumUpCost()/sumUpChange());
#elif TYPE2==2
      addToSumUpCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+(ceil(value)-value)));
      setUpDynamicPseudoCost(sumUpCost()*(TYPERATIO/sumUpChange()+(1.0-TYPERATIO)/(double) numberTimesUp()));
#endif
    } else {
      //printf("(up infeasible value down %g ",change,ceil(value)-value);
      incrementNumberTimesUpInfeasible();
#if INFEAS==2
      double distanceToCutoff=0.0;
      double objectiveValue = model->getCurrentMinimizationObjValue();
      distanceToCutoff =  model->getCutoff()  - originalValue;
      if (distanceToCutoff<1.0e20) 
	change = distanceToCutoff*2.0;
      else 
	change = upDynamicPseudoCost()*(ceil(value)-value)*10.0; 
      change = CoinMax(1.0e-12*(1.0+fabs(originalValue)),change);
      incrementNumberTimesUp();
      addToSumUpChange(1.0e-30+ceil(value)-value);
      addToSumUpDecrease(data.intDecrease_);
#if TYPE2==0
      addToSumUpCost(change/(1.0e-30+(ceil(value)-value)));
      setUpDynamicPseudoCost(sumUpCost()/(double) numberTimesUp());
#elif TYPE2==1
      addToSumUpCost(change);
      setUpDynamicPseudoCost(sumUpCost()/sumUpChange());
#elif TYPE2==2
      addToSumUpCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+(ceil(value)-value)));
      setUpDynamicPseudoCost(sumUpCost()*(TYPERATIO/sumUpChange()+(1.0-TYPERATIO)/(double) numberTimesUp()));
#endif
#endif
    }
  }
  //print(1,0.5);
}
// Pass in probing information
void 
CbcSimpleIntegerDynamicPseudoCost::setProbingInformation(int fixedDown, int fixedUp)
{
  numberTimesProbingTotal_++;
  numberTimesDownLocalFixed_ = fixedDown;
  numberTimesDownTotalFixed_ += fixedDown;
  numberTimesUpLocalFixed_ = fixedUp;
  numberTimesUpTotalFixed_ += fixedUp;
}
// Print
void 
CbcSimpleIntegerDynamicPseudoCost::print(int type,double value) const
{
  if (!type) {
    double meanDown =0.0;
    double devDown =0.0;
    if (numberTimesDown_) {
      meanDown = sumDownCost_/(double) numberTimesDown_;
      devDown = meanDown*meanDown + sumDownCostSquared_ - 
        2.0*meanDown*sumDownCost_;
      if (devDown>=0.0)
        devDown = sqrt(devDown);
    }
    double meanUp =0.0;
    double devUp =0.0;
    if (numberTimesUp_) {
      meanUp = sumUpCost_/(double) numberTimesUp_;
      devUp = meanUp*meanUp + sumUpCostSquared_ - 
        2.0*meanUp*sumUpCost_;
      if (devUp>=0.0)
        devUp = sqrt(devUp);
    }
#if 0
    printf("%d down %d times (%d inf) mean %g (dev %g) up %d times (%d inf) mean %g (dev %g)\n",
           columnNumber_,
           numberTimesDown_,numberTimesDownInfeasible_,meanDown,devDown,
           numberTimesUp_,numberTimesUpInfeasible_,meanUp,devUp);
#else
    printf("%d down %d times (%d inf) mean %g  up %d times (%d inf) mean %g\n",
           columnNumber_,
           numberTimesDown_,numberTimesDownInfeasible_,meanDown,
           numberTimesUp_,numberTimesUpInfeasible_,meanUp);
#endif
  } else {
    const double * upper = model_->getCbcColUpper();
    double integerTolerance = 
      model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double below = floor(value+integerTolerance);
    double above = below+1.0;
    if (above>upper[columnNumber_]) {
      above=below;
      below = above -1;
    }
    double objectiveValue = model_->getCurrentMinimizationObjValue();
    double distanceToCutoff =  model_->getCutoff() - objectiveValue;
    if (distanceToCutoff<1.0e20) 
      distanceToCutoff *= 10.0;
    else 
      distanceToCutoff = 1.0e2 + fabs(objectiveValue);
    double sum;
    int number;
    double downCost = CoinMax(value-below,0.0);
    double downCost0 = downCost*downDynamicPseudoCost_;
    sum = sumDownCost();
    number = numberTimesDown();
    sum += numberTimesDownInfeasible()*(distanceToCutoff/(downCost+1.0e-12));
    if (number>0)
      downCost *= sum / (double) number;
    else
      downCost  *=  downDynamicPseudoCost_;
    double upCost = CoinMax((above-value),0.0);
    double upCost0 = upCost*upDynamicPseudoCost_;
    sum = sumUpCost();
    number = numberTimesUp();
    sum += numberTimesUpInfeasible()*(distanceToCutoff/(upCost+1.0e-12));
    if (number>0)
      upCost *= sum / (double) number;
    else
      upCost  *=  upDynamicPseudoCost_;
    printf("%d down %d times %g (est %g)  up %d times %g (est %g)\n",
           columnNumber_,
           numberTimesDown_,downCost,downCost0,
           numberTimesUp_,upCost,upCost0);
  }
}

// Default Constructor 
CbcDynamicPseudoCostBranchingObject::CbcDynamicPseudoCostBranchingObject()
  :CbcIntegerBranchingObject()
{
  changeInGuessed_=1.0e-5;
  object_=NULL;
}

// Useful constructor
CbcDynamicPseudoCostBranchingObject::CbcDynamicPseudoCostBranchingObject (CbcModel * model, 
                                                                          int variable, 
                                                                          int way , double value,
                                       CbcSimpleIntegerDynamicPseudoCost * object) 
  :CbcIntegerBranchingObject(model,variable,way,value)
{
  changeInGuessed_=1.0e-5;
  object_=object;
}
// Useful constructor for fixing
CbcDynamicPseudoCostBranchingObject::CbcDynamicPseudoCostBranchingObject (CbcModel * model, 
						      int variable, int way,
						      double lowerValue, 
						      double upperValue)
  :CbcIntegerBranchingObject(model,variable,way,lowerValue)
{
  changeInGuessed_=1.0e100;
  object_=NULL;
}
  

// Copy constructor 
CbcDynamicPseudoCostBranchingObject::CbcDynamicPseudoCostBranchingObject ( 
				 const CbcDynamicPseudoCostBranchingObject & rhs)
  :CbcIntegerBranchingObject(rhs)
{
  changeInGuessed_ = rhs.changeInGuessed_;
  object_=rhs.object_;
}

// Assignment operator 
CbcDynamicPseudoCostBranchingObject & 
CbcDynamicPseudoCostBranchingObject::operator=( const CbcDynamicPseudoCostBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcIntegerBranchingObject::operator=(rhs);
    changeInGuessed_ = rhs.changeInGuessed_;
    object_=rhs.object_;
  }
  return *this;
}
CbcBranchingObject * 
CbcDynamicPseudoCostBranchingObject::clone() const
{ 
  return (new CbcDynamicPseudoCostBranchingObject(*this));
}


// Destructor 
CbcDynamicPseudoCostBranchingObject::~CbcDynamicPseudoCostBranchingObject ()
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
int 
CbcDynamicPseudoCostBranchingObject::fillStrongInfo( CbcStrongInfo & info)
{
  assert (object_);
  assert (info.possibleBranch==this);
  if (object_->numberTimesUp()<object_->numberBeforeTrust()||
      object_->numberTimesDown()<object_->numberBeforeTrust()) {
    return 0;
  } else {
    info.upMovement = object_->upDynamicPseudoCost()*(ceil(value_)-value_);
    info.downMovement = object_->downDynamicPseudoCost()*(value_-floor(value_));
    info.numIntInfeasUp  -= (int) (object_->sumUpDecrease()/
                                   ((double) object_->numberTimesUp()));
    info.numObjInfeasUp = 0;
    info.finishedUp = false;
    info.numItersUp = 0;
    info.numIntInfeasDown  -= (int) (object_->sumDownDecrease()/
                                   ((double) object_->numberTimesDown()));
    info.numObjInfeasDown = 0;
    info.finishedDown = false;
    info.numItersDown = 0;
    info.fix =0;
    return 1;
  }
}
  
// Default Constructor 
CbcBranchDynamicDecision::CbcBranchDynamicDecision()
  :CbcBranchDecision()
{
  bestCriterion_ = 0.0;
  bestChangeUp_ = 0.0;
  bestNumberUp_ = 0;
  bestChangeDown_ = 0.0;
  bestNumberDown_ = 0;
  bestObject_ = NULL;
}

// Copy constructor 
CbcBranchDynamicDecision::CbcBranchDynamicDecision (
				    const CbcBranchDynamicDecision & rhs)
  :CbcBranchDecision()
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
void 
CbcBranchDynamicDecision::initialize(CbcModel * model)
{
  bestCriterion_ = 0.0;
  bestChangeUp_ = 0.0;
  bestNumberUp_ = 0;
  bestChangeDown_ = 0.0;
  bestNumberDown_ = 0;
  bestObject_ = NULL;
}

/* Saves a clone of current branching object.  Can be used to update
      information on object causing branch - after branch */
void 
CbcBranchDynamicDecision::saveBranchingObject(OsiBranchingObject * object) 
{
  OsiBranchingObject * obj = object->clone();
  CbcDynamicPseudoCostBranchingObject * branchingObject =
    dynamic_cast<CbcDynamicPseudoCostBranchingObject *>(obj);
  assert (branchingObject);
  object_=branchingObject;
}
/* Pass in information on branch just done.
   assumes object can get information from solver */
void 
CbcBranchDynamicDecision::updateInformation(OsiSolverInterface * solver,
                                            const CbcNode * node)
{
  assert (object_);
  const CbcModel * model = object_->model();
  double originalValue=node->objectiveValue();
  int originalUnsatisfied = node->numberUnsatisfied();
  double objectiveValue = solver->getObjValue()*model->getObjSense();
  int unsatisfied=0;
  int i;
  int numberIntegers = model->numberIntegers();;
  const double * solution = solver->getColSolution();
  //const double * lower = solver->getColLower();
  //const double * upper = solver->getColUpper();
  CbcDynamicPseudoCostBranchingObject * branchingObject =
    dynamic_cast<CbcDynamicPseudoCostBranchingObject *>(object_);
  if (!branchingObject) {
    delete object_;
    object_=NULL;
    return;
  }
  CbcSimpleIntegerDynamicPseudoCost *  object = branchingObject->object();
  double change = CoinMax(0.0,objectiveValue-originalValue);
  // probably should also ignore if stopped
  int iStatus;
  if (solver->isProvenOptimal())
    iStatus=0; // optimal
  else if (solver->isIterationLimitReached()
           &&!solver->isDualObjectiveLimitReached())
    iStatus=2; // unknown 
  else
    iStatus=1; // infeasible

  bool feasible = iStatus!=1;
  if (feasible) {
    double integerTolerance = 
      model->getDblParam(CbcModel::CbcIntegerTolerance);
    const int * integerVariable = model->integerVariable();
    for (i=0;i<numberIntegers;i++) {
      int j=integerVariable[i];
      double value = solution[j];
      double nearest = floor(value+0.5);
      if (fabs(value-nearest)>integerTolerance) 
        unsatisfied++;
    }
  }
  int way = object_->way();
  double value = object_->value();
  //#define TYPE2 1
  //#define TYPERATIO 0.9
  if (way<0) {
    // down
    if (feasible) {
      //printf("(down change %g value down %g ",change,value-floor(value));
      object->incrementNumberTimesDown();
      object->addToSumDownChange(1.0e-30+value-floor(value));
      object->addToSumDownDecrease(originalUnsatisfied-unsatisfied);
#if TYPE2==0
      object->addToSumDownCost(change/(1.0e-30+(value-floor(value))));
      object->setDownDynamicPseudoCost(object->sumDownCost()/(double) object->numberTimesDown());
#elif TYPE2==1
      object->addToSumDownCost(change);
      object->setDownDynamicPseudoCost(object->sumDownCost()/object->sumDownChange());
#elif TYPE2==2
      object->addToSumDownCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+(value-floor(value))));
      object->setDownDynamicPseudoCost(object->sumDownCost()*(TYPERATIO/object->sumDownChange()+(1.0-TYPERATIO)/(double) object->numberTimesDown()));
#endif
    } else {
      //printf("(down infeasible value down %g ",change,value-floor(value));
      object->incrementNumberTimesDownInfeasible();
#if INFEAS==2
      double distanceToCutoff=0.0;
      double objectiveValue = model->getCurrentMinimizationObjValue();
      distanceToCutoff =  model->getCutoff()  - originalValue;
      if (distanceToCutoff<1.0e20) 
	change = distanceToCutoff*2.0;
      else 
	change = object->downDynamicPseudoCost()*(value-floor(value))*10.0; 
      change = CoinMax(1.0e-12*(1.0+fabs(originalValue)),change);
      object->incrementNumberTimesDown();
      object->addToSumDownChange(1.0e-30+value-floor(value));
      object->addToSumDownDecrease(originalUnsatisfied-unsatisfied);
#if TYPE2==0
      object->addToSumDownCost(change/(1.0e-30+(value-floor(value))));
      object->setDownDynamicPseudoCost(object->sumDownCost()/(double) object->numberTimesDown());
#elif TYPE2==1
      object->addToSumDownCost(change);
      object->setDownDynamicPseudoCost(object->sumDownCost()/object->sumDownChange());
#elif TYPE2==2
      object->addToSumDownCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+(value-floor(value))));
      object->setDownDynamicPseudoCost(object->sumDownCost()*(TYPERATIO/object->sumDownChange()+(1.0-TYPERATIO)/(double) object->numberTimesDown()));
#endif
#endif
    }
  } else {
    // up
    if (feasible) {
      //printf("(up change %g value down %g ",change,ceil(value)-value);
      object->incrementNumberTimesUp();
      object->addToSumUpChange(1.0e-30+ceil(value)-value);
      object->addToSumUpDecrease(unsatisfied-originalUnsatisfied);
#if TYPE2==0
      object->addToSumUpCost(change/(1.0e-30+(ceil(value)-value)));
      object->setUpDynamicPseudoCost(object->sumUpCost()/(double) object->numberTimesUp());
#elif TYPE2==1
      object->addToSumUpCost(change);
      object->setUpDynamicPseudoCost(object->sumUpCost()/object->sumUpChange());
#elif TYPE2==2
      object->addToSumUpCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+(ceil(value)-value)));
      object->setUpDynamicPseudoCost(object->sumUpCost()*(TYPERATIO/object->sumUpChange()+(1.0-TYPERATIO)/(double) object->numberTimesUp()));
#endif
    } else {
      //printf("(up infeasible value down %g ",change,ceil(value)-value);
      object->incrementNumberTimesUpInfeasible();
#if INFEAS==2
      double distanceToCutoff=0.0;
      double objectiveValue = model->getCurrentMinimizationObjValue();
      distanceToCutoff =  model->getCutoff()  - originalValue;
      if (distanceToCutoff<1.0e20) 
	change = distanceToCutoff*2.0;
      else 
	change = object->upDynamicPseudoCost()*(ceil(value)-value)*10.0; 
      change = CoinMax(1.0e-12*(1.0+fabs(originalValue)),change);
      object->incrementNumberTimesUp();
      object->addToSumUpChange(1.0e-30+ceil(value)-value);
      object->addToSumUpDecrease(unsatisfied-originalUnsatisfied);
#if TYPE2==0
      object->addToSumUpCost(change/(1.0e-30+(ceil(value)-value)));
      object->setUpDynamicPseudoCost(object->sumUpCost()/(double) object->numberTimesUp());
#elif TYPE2==1
      object->addToSumUpCost(change);
      object->setUpDynamicPseudoCost(object->sumUpCost()/object->sumUpChange());
#elif TYPE2==2
      object->addToSumUpCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+(ceil(value)-value)));
      object->setUpDynamicPseudoCost(object->sumUpCost()*(TYPERATIO/object->sumUpChange()+(1.0-TYPERATIO)/(double) object->numberTimesUp()));
#endif
#endif
    }
  }
  //object->print(1,0.5);
  delete object_;
  object_=NULL;
}

/*
  Simple dynamic decision algorithm. Compare based on infeasibility (numInfUp,
  numInfDown) until a solution is found by search, then switch to change in
  objective (changeUp, changeDown). Note that bestSoFar is remembered in
  bestObject_, so the parameter bestSoFar is unused.
*/

int
CbcBranchDynamicDecision::betterBranch(CbcBranchingObject * thisOne,
			    CbcBranchingObject * bestSoFar,
			    double changeUp, int numInfUp,
			    double changeDown, int numInfDown)
{
  CbcModel * model = thisOne->model();
  int stateOfSearch = model->stateOfSearch();
  int betterWay=0;
  double value=0.0;
  if (!bestObject_) {
    bestCriterion_=-1.0;
    bestNumberUp_=COIN_INT_MAX;
    bestNumberDown_=COIN_INT_MAX;
  }
  if (stateOfSearch<=1&&thisOne->model()->currentNode()->depth()>=8) {
#define TRY_STUFF 1
#ifdef TRY_STUFF
    // before solution - choose smallest number 
    // could add in depth as well
    int bestNumber = CoinMin(bestNumberUp_,bestNumberDown_);
    if (numInfUp<numInfDown) {
      if (numInfUp<bestNumber) {
	betterWay = 1;
      } else if (numInfUp==bestNumber) {
	if (changeUp<bestChangeUp_)
	  betterWay=1;
      }
    } else if (numInfUp>numInfDown) {
      if (numInfDown<bestNumber) {
	betterWay = -1;
      } else if (numInfDown==bestNumber) {
	if (changeDown<bestChangeDown_)
	  betterWay=-1;
      }
    } else {
      // up and down have same number
      bool better=false;
      if (numInfUp<bestNumber) {
	better=true;
      } else if (numInfUp==bestNumber) {
	if (CoinMin(changeUp,changeDown)<CoinMin(bestChangeUp_,bestChangeDown_)-1.0e-5)
	  better=true;;
      }
      if (better) {
	// see which way
	if (changeUp<=changeDown)
	  betterWay=1;
	else
	  betterWay=-1;
      }
    }
    if (betterWay) {
      value = CoinMin(numInfUp,numInfDown);
    }
#else
    // got a solution
    double minValue = CoinMin(changeDown,changeUp);
    double maxValue = CoinMax(changeDown,changeUp);
    value = WEIGHT_BEFORE*minValue + (1.0-WEIGHT_BEFORE)*maxValue;
    if (value>bestCriterion_+1.0e-8) {
      if (changeUp<=changeDown) {
        betterWay=1;
      } else {
        betterWay=-1;
      }
    }
#endif
  } else {
#if TRY_STUFF > 1
    // Get current number of infeasibilities, cutoff and current objective
    CbcNode * node = model->currentNode();
    int numberUnsatisfied = node->numberUnsatisfied();
    double cutoff = model->getCutoff();
    double objectiveValue = node->objectiveValue();
#endif
    // got a solution
    double minValue = CoinMin(changeDown,changeUp);
    double maxValue = CoinMax(changeDown,changeUp);
    // Reduce
#ifdef TRY_STUFF
    //maxValue = CoinMin(maxValue,minValue*4.0);
#else
    maxValue = CoinMin(maxValue,minValue*2.0);
#endif
    value = WEIGHT_AFTER*minValue + (1.0-WEIGHT_AFTER)*maxValue;
    double useValue = value;
    double useBest = bestCriterion_;
#if TRY_STUFF>1
    int thisNumber = CoinMin(numInfUp,numInfDown);
    int bestNumber = CoinMin(bestNumberUp_,bestNumberDown_);
    double distance = cutoff-objectiveValue;
    assert (distance>=0.0);
    if (useValue+0.1*distance>useBest&&useValue*1.1>useBest&&
	useBest+0.1*distance>useValue&&useBest*1.1>useValue) {
      // not much in it - look at unsatisfied
      if (thisNumber<numberUnsatisfied||bestNumber<numberUnsatisfied) {
	double perInteger = distance/ ((double) numberUnsatisfied);
	useValue += thisNumber*perInteger;
	useBest += bestNumber*perInteger;
      }
    }
#endif
    if (useValue>useBest+1.0e-8) {
      if (changeUp<=changeDown) {
        betterWay=1;
      } else {
        betterWay=-1;
      }
    }
  }
  if (betterWay) {
    // maybe change better way
    CbcDynamicPseudoCostBranchingObject * branchingObject =
      dynamic_cast<CbcDynamicPseudoCostBranchingObject *>(thisOne);
    if (branchingObject) {
      CbcSimpleIntegerDynamicPseudoCost *  object = branchingObject->object();
      double separator = object->upDownSeparator();
      if (separator>0.0) {
	const double * solution = thisOne->model()->testSolution();
	double valueVariable = solution[object->columnNumber()];
	betterWay = (valueVariable-floor(valueVariable)>=separator) ? 1 : -1;
      }
    }
    bestCriterion_ = value;
    bestChangeUp_ = changeUp;
    bestNumberUp_ = numInfUp;
    bestChangeDown_ = changeDown;
    bestNumberDown_ = numInfDown;
    bestObject_=thisOne;
    // See if user is overriding way
    if (thisOne->object()&&thisOne->object()->preferredWay())
      betterWay = thisOne->object()->preferredWay();
  }
  return betterWay;
}
/* Sets or gets best criterion so far */
void 
CbcBranchDynamicDecision::setBestCriterion(double value)
{ 
  bestCriterion_ = value;
}
double 
CbcBranchDynamicDecision::getBestCriterion() const
{ 
  return bestCriterion_;
}
