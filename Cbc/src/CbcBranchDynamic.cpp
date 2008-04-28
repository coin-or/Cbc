// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
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
#ifdef COIN_DEVELOP
typedef struct {
  char where_;
  char status_;
  unsigned short sequence_;
  int numberUp_;
  int numberUpInf_;
  float sumUp_;
  float upEst_; // or change in obj in update
  int numberDown_;
  int numberDownInf_;
  float sumDown_;
  float downEst_; // or movement in value in update
} History;
History * history=NULL;
int numberHistory=0;
int maxHistory=0;
bool getHistoryStatistics_=true;
static void increaseHistory() {
  if (numberHistory==maxHistory) {
    maxHistory = 100+(3*maxHistory)/2;
    History * temp = new History [maxHistory];
    memcpy(temp,history,numberHistory*sizeof(History));
    delete [] history;
    history=temp;
  }
}
static bool addRecord(History newOne) {
  //if (!getHistoryStatistics_)
    return false;
  bool fromCompare=false;
  int i;
  for ( i=numberHistory-1;i>=0;i--) {
    if (newOne.sequence_!=history[i].sequence_)
      continue;
    if (newOne.where_!=history[i].where_)
      continue;
    if (newOne.numberUp_!=history[i].numberUp_)
      continue;
    if (newOne.sumUp_!=history[i].sumUp_)
      continue;
    if (newOne.numberUpInf_!=history[i].numberUpInf_)
      continue;
    if (newOne.upEst_!=history[i].upEst_)
      continue;
    if (newOne.numberDown_!=history[i].numberDown_)
      continue;
    if (newOne.sumDown_!=history[i].sumDown_)
      continue;
    if (newOne.numberDownInf_!=history[i].numberDownInf_)
      continue;
    if (newOne.downEst_!=history[i].downEst_)
      continue;
    // If B knock out previous B
    if (newOne.where_=='C') {
      fromCompare=true;
      if (newOne.status_=='B') {
	int j;
	for (j=i-1;j>=0;j--) {
	  if (history[j].where_=='C') {
	    if (history[j].status_=='I') {
	      break;
	    } else if (history[j].status_=='B') {
	      history[j].status_=' ';
	      break;
	    }
	  }
	}
      }
      break;
    }
  }
  if (i==-1||fromCompare) {
    //add
    increaseHistory();
    history[numberHistory++]=newOne;
    return true;
  } else {
    return false;
  }
}
#endif
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
#ifdef CBC_INSTRUMENT
  numberTimesInfeasible_=0;
#endif
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
#ifdef CBC_INSTRUMENT
  numberTimesInfeasible_=0;
#endif
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
#define TYPE2 0
#if TYPE2==0
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
#ifdef CBC_INSTRUMENT
  numberTimesInfeasible_=0;
#endif
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
#if TYPE2==0
  // No
  sumUpCost_ = 0.0;
  sumUpChange_ = 0.0;
  numberTimesUp_ = 0;
  sumDownCost_ = 0.0;
  sumDownChange_ = 0.0;
  numberTimesDown_ = 0;
  sumUpCost_ = 1.0e-4*upDynamicPseudoCost_;
  sumDownCost_ = 1.0e-4*downDynamicPseudoCost_;
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
#ifdef CBC_INSTRUMENT
  numberTimesInfeasible_=rhs.numberTimesInfeasible_;
#endif
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
#ifdef CBC_INSTRUMENT
    numberTimesInfeasible_=rhs.numberTimesInfeasible_;
#endif
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
CbcSimpleIntegerDynamicPseudoCost::copySome(const CbcSimpleIntegerDynamicPseudoCost * otherObject)
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
// Updates stuff like pseudocosts before threads
void 
CbcSimpleIntegerDynamicPseudoCost::updateBefore(const OsiObject * rhs) 
{
  const CbcSimpleIntegerDynamicPseudoCost * rhsObject =
    dynamic_cast <const CbcSimpleIntegerDynamicPseudoCost *>(rhs) ;
  assert (rhsObject);
  copySome(rhsObject);
}
  // was 1 - but that looks flakey
#define INFEAS 1
// Updates stuff like pseudocosts after threads finished
void 
CbcSimpleIntegerDynamicPseudoCost::updateAfter(const OsiObject * rhs, const OsiObject * baseObjectX) 
{
  const CbcSimpleIntegerDynamicPseudoCost * rhsObject =
    dynamic_cast <const CbcSimpleIntegerDynamicPseudoCost *>(rhs) ;
  assert (rhsObject);
  const CbcSimpleIntegerDynamicPseudoCost * baseObject =
    dynamic_cast <const CbcSimpleIntegerDynamicPseudoCost *>(baseObjectX) ;
  assert (baseObject);
  // compute current
  double sumDown = downDynamicPseudoCost_*(numberTimesDown_+numberTimesDownInfeasible_);
  sumDown -= baseObject->downDynamicPseudoCost_*(baseObject->numberTimesDown_+baseObject->numberTimesDownInfeasible_);
  sumDown = CoinMax(sumDown,0.0);
  sumDown += rhsObject->downDynamicPseudoCost_*(rhsObject->numberTimesDown_+rhsObject->numberTimesDownInfeasible_);
  double sumUp = upDynamicPseudoCost_*(numberTimesUp_+numberTimesUpInfeasible_);
  sumUp -= baseObject->upDynamicPseudoCost_*(baseObject->numberTimesUp_+baseObject->numberTimesUpInfeasible_);
  sumUp += rhsObject->upDynamicPseudoCost_*(rhsObject->numberTimesUp_+rhsObject->numberTimesUpInfeasible_);
  sumUp = CoinMax(sumUp,0.0);
  sumDownCost_ += rhsObject->sumDownCost_-baseObject->sumDownCost_;
  sumUpCost_ += rhsObject->sumUpCost_-baseObject->sumUpCost_;
  sumDownChange_ += rhsObject->sumDownChange_-baseObject->sumDownChange_;
  sumUpChange_ += rhsObject->sumUpChange_-baseObject->sumUpChange_;
  sumDownCostSquared_ += rhsObject->sumDownCostSquared_-baseObject->sumDownCostSquared_;
  sumUpCostSquared_ += rhsObject->sumUpCostSquared_-baseObject->sumUpCostSquared_;
  sumDownDecrease_ += rhsObject->sumDownDecrease_-baseObject->sumDownDecrease_;
  sumUpDecrease_ += rhsObject->sumUpDecrease_-baseObject->sumUpDecrease_;
  lastDownCost_ += rhsObject->lastDownCost_-baseObject->lastDownCost_;
  lastUpCost_ += rhsObject->lastUpCost_-baseObject->lastUpCost_;
  lastDownDecrease_ += rhsObject->lastDownDecrease_-baseObject->lastDownDecrease_;
  lastUpDecrease_ += rhsObject->lastUpDecrease_-baseObject->lastUpDecrease_;
  numberTimesDown_ += rhsObject->numberTimesDown_-baseObject->numberTimesDown_;
  numberTimesUp_ += rhsObject->numberTimesUp_-baseObject->numberTimesUp_;
  numberTimesDownInfeasible_ += rhsObject->numberTimesDownInfeasible_-baseObject->numberTimesDownInfeasible_;
  numberTimesUpInfeasible_ += rhsObject->numberTimesUpInfeasible_-baseObject->numberTimesUpInfeasible_;
  numberTimesDownLocalFixed_ += rhsObject->numberTimesDownLocalFixed_-baseObject->numberTimesDownLocalFixed_;
  numberTimesUpLocalFixed_ += rhsObject->numberTimesUpLocalFixed_-baseObject->numberTimesUpLocalFixed_;
  numberTimesDownTotalFixed_ += rhsObject->numberTimesDownTotalFixed_-baseObject->numberTimesDownTotalFixed_;
  numberTimesUpTotalFixed_ += rhsObject->numberTimesUpTotalFixed_-baseObject->numberTimesUpTotalFixed_;
  numberTimesProbingTotal_ += rhsObject->numberTimesProbingTotal_-baseObject->numberTimesProbingTotal_;
  if (numberTimesDown_+numberTimesDownInfeasible_>0) {
    setDownDynamicPseudoCost(sumDown/(double) (numberTimesDown_+numberTimesDownInfeasible_));
  }
  if (numberTimesUp_+numberTimesUpInfeasible_>0) {
    setUpDynamicPseudoCost(sumUp/(double) (numberTimesUp_+numberTimesUpInfeasible_));
  }
  //printf("XX %d down %d %d %g up %d %d %g\n",columnNumber_,numberTimesDown_,numberTimesDownInfeasible_,downDynamicPseudoCost_,
  // numberTimesUp_,numberTimesUpInfeasible_,upDynamicPseudoCost_);
  assert (downDynamicPseudoCost_>1.0e-40&&upDynamicPseudoCost_>1.0e-40);
}
// Same - returns true if contents match(ish)
bool 
CbcSimpleIntegerDynamicPseudoCost::same(const CbcSimpleIntegerDynamicPseudoCost * otherObject) const
{
  bool okay = true;
  if (downDynamicPseudoCost_!=otherObject->downDynamicPseudoCost_)
    okay=false;
  if (upDynamicPseudoCost_!=otherObject->upDynamicPseudoCost_)
    okay=false;
  if (sumDownCost_!= otherObject->sumDownCost_)
    okay=false;
  if (sumUpCost_!= otherObject->sumUpCost_)
    okay=false;
  if (sumDownChange_!= otherObject->sumDownChange_)
    okay=false;
  if (sumUpChange_!= otherObject->sumUpChange_)
    okay=false;
  if (sumDownCostSquared_!= otherObject->sumDownCostSquared_)
    okay=false;
  if (sumUpCostSquared_!= otherObject->sumUpCostSquared_)
    okay=false;
  if (sumDownDecrease_!= otherObject->sumDownDecrease_)
    okay=false;
  if (sumUpDecrease_!= otherObject->sumUpDecrease_)
    okay=false;
  if (lastDownCost_!= otherObject->lastDownCost_)
    okay=false;
  if (lastUpCost_!= otherObject->lastUpCost_)
    okay=false;
  if (lastDownDecrease_!= otherObject->lastDownDecrease_)
    okay=false;
  if (lastUpDecrease_!= otherObject->lastUpDecrease_)
    okay=false;
  if (numberTimesDown_!= otherObject->numberTimesDown_)
    okay=false;
  if (numberTimesUp_!= otherObject->numberTimesUp_)
    okay=false;
  if (numberTimesDownInfeasible_!= otherObject->numberTimesDownInfeasible_)
    okay=false;
  if (numberTimesUpInfeasible_!= otherObject->numberTimesUpInfeasible_)
    okay=false;
  if (numberTimesDownLocalFixed_!= otherObject->numberTimesDownLocalFixed_)
    okay=false;
  if (numberTimesUpLocalFixed_!= otherObject->numberTimesUpLocalFixed_)
    okay=false;
  if (numberTimesDownTotalFixed_!= otherObject->numberTimesDownTotalFixed_)
    okay=false;
  if (numberTimesUpTotalFixed_!= otherObject->numberTimesUpTotalFixed_)
    okay=false;
  if (numberTimesProbingTotal_!= otherObject->numberTimesProbingTotal_)
    okay=false;
  return okay;
}
// Creates a branching objecty
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
//#define FUNNY_BRANCHING  
// Infeasibility - large is 0.5
double 
CbcSimpleIntegerDynamicPseudoCost::infeasibility(int & preferredWay) const
{
  assert (downDynamicPseudoCost_>1.0e-40&&upDynamicPseudoCost_>1.0e-40);
  const double * solution = model_->testSolution();
  const double * lower = model_->getCbcColLower();
  const double * upper = model_->getCbcColUpper();
#ifdef FUNNY_BRANCHING
  const double * dj = model_->getCbcReducedCost();
  double djValue = dj[columnNumber_];
  lastDownDecrease_++;
  if (djValue>1.0e-6) {
    // wants to go down
    if (true||lower[columnNumber_]>originalLower_) {
      // Lower bound active
      lastUpDecrease_++;
      sumDownCostSquared_ += djValue;
    }
  } else if (djValue<-1.0e-6) {
    // wants to go up
    if (true||upper[columnNumber_]<originalUpper_) {
      // Upper bound active
      lastUpDecrease_++;
      sumUpCostSquared_ -= djValue;
    }
  }
#endif
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
#if INFEAS==1
  double distanceToCutoff=0.0;
  double objectiveValue = model_->getCurrentMinimizationObjValue();
  distanceToCutoff =  model_->getCutoff()  - objectiveValue;
  if (distanceToCutoff<1.0e20) 
    distanceToCutoff *= 10.0;
  else 
    distanceToCutoff = 1.0e2 + fabs(objectiveValue);
  distanceToCutoff = CoinMax(distanceToCutoff,1.0e-12*(1.0+fabs(objectiveValue)));
#endif
  double sum;
  double number;
  double downCost = CoinMax(value-below,0.0);
#if TYPE2==0
  sum = sumDownCost_;
  number = numberTimesDown_;
#if INFEAS==1
  sum += numberTimesDownInfeasible_*CoinMax(distanceToCutoff/(downCost+1.0e-12),sumDownCost_);
  number += numberTimesDownInfeasible_;
#endif
#elif TYPE2==1
  sum = sumDownCost_;
  number = sumDownChange_;
#if INFEAS==1
  sum += numberTimesDownInfeasible_*CoinMax(distanceToCutoff/(downCost+1.0e-12),sumDownCost_);
  number += numberTimesDownInfeasible_;
#endif
#elif TYPE2==2
  abort();
#if INFEAS==1
  sum += numberTimesDownInfeasible_*(distanceToCutoff/(downCost+1.0e-12));
  number += numberTimesDownInfeasible_;
#endif
#endif
  if (number>0.0)
    downCost *= sum / number;
  else
    downCost  *=  downDynamicPseudoCost_;
  double upCost = CoinMax((above-value),0.0);
#if TYPE2==0
  sum = sumUpCost_;
  number = numberTimesUp_;
#if INFEAS==1
  sum += numberTimesUpInfeasible_*CoinMax(distanceToCutoff/(upCost+1.0e-12),sumUpCost_);
  number += numberTimesUpInfeasible_;
#endif
#elif TYPE2==1
  sum = sumUpCost_;
  number = sumUpChange_;
#if INFEAS==1
  sum += numberTimesUpInfeasible_*CoinMax(distanceToCutoff/(upCost+1.0e-12),sumUpCost_);
  number += numberTimesUpInfeasible_;
#endif
#elif TYPE2==1
  abort();
#if INFEAS==1
  sum += numberTimesUpInfeasible_*(distanceToCutoff/(upCost+1.0e-12));
  number += numberTimesUpInfeasible_;
#endif
#endif
  if (number>0.0)
    upCost *= sum / number;
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
#ifdef FUNNY_BRANCHING
  if (fabs(value-nearest)>integerTolerance) {
    double ratio = (100.0+lastUpDecrease_)/(100.0+lastDownDecrease_);
    downCost *= ratio;
    upCost *= ratio;
    if ((lastUpDecrease_%100)==-1) 
      printf("col %d total %d djtimes %d down %g up %g\n",
	     columnNumber_,lastDownDecrease_,lastUpDecrease_,
	     sumDownCostSquared_,sumUpCostSquared_);
  }
#endif
  if (preferredWay_)
    preferredWay=preferredWay_;
  // weight at 1.0 is max min
#define WEIGHT_AFTER 0.7
#define WEIGHT_BEFORE 0.1
  if (fabs(value-nearest)<=integerTolerance) {
    return 0.0;
  } else {
#ifdef CBC_INSTRUMENT
    numberTimesInfeasible_++;
#endif
    int stateOfSearch = model_->stateOfSearch()%10;
    double returnValue=0.0;
    double minValue = CoinMin(downCost,upCost);
    double maxValue = CoinMax(downCost,upCost);
    char where;
    // was <= 10
    //if (stateOfSearch<=1||model_->currentNode()->depth()<=-10 /* was ||maxValue>0.2*distanceToCutoff*/) {
    if (stateOfSearch<=2) {
      // no branching solution
      where='i';
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
      where='I';
      returnValue = WEIGHT_AFTER*minValue + (1.0-WEIGHT_AFTER)*maxValue;
    }
    if (numberTimesUp_<numberBeforeTrust_||
        numberTimesDown_<numberBeforeTrust_) {
      //if (returnValue<1.0e10)
      //returnValue += 1.0e12;
      //else
      returnValue *= 1.0e3;
      if (!numberTimesUp_&&!numberTimesDown_)
        returnValue *= 1.0e10;
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
#ifdef CBC_INSTRUMENT
    int nn = numberTimesInfeasible_  - CoinMax(numberTimesUp_,numberTimesDown_);
    if (nn<0) {
      // Something to do with parallel synchronization
      numberTimesInfeasible_  = CoinMax(numberTimesUp_,numberTimesDown_);
    } else if (nn) {
      returnValue *= sqrt((double) nn);
    }
#endif
#ifdef COIN_DEVELOP
    History hist;
    hist.where_=where;
    hist.status_=' ';
    hist.sequence_=columnNumber_;
    hist.numberUp_=numberTimesUp_;
    hist.numberUpInf_=numberTimesUpInfeasible_;
    hist.sumUp_=sumUpCost_;
    hist.upEst_=upCost;
    hist.numberDown_=numberTimesDown_;
    hist.numberDownInf_=numberTimesDownInfeasible_;
    hist.sumDown_=sumDownCost_;
    hist.downEst_=downCost;
    if (stateOfSearch)
      addRecord(hist);
#endif
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
  distanceToCutoff = CoinMax(distanceToCutoff,1.0e-12*(1.0+fabs(objectiveValue)));
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
  CbcObjectUpdateData newData (this, way,
			       change, iStatus,
			       originalUnsatisfied-unsatisfied,value);
  newData.originalObjective_ = originalValue;
  // Solvers know about direction
  double direction = solver->getObjSense();
  solver->getDblParam(OsiDualObjectiveLimit,newData.cutoff_);
  newData.cutoff_ *= direction;
  return newData;
}
// Update object by CbcObjectUpdateData
void 
CbcSimpleIntegerDynamicPseudoCost::updateInformation(const CbcObjectUpdateData & data)
{
  bool feasible = data.status_!=1;
  int way = data.way_;
  double value = data.branchingValue_;
  double change = data.change_;
#define TYPERATIO 0.9
#define MINIMUM_MOVEMENT 0.1
#ifdef COIN_DEVELOP
  History hist;
  hist.where_='U'; // need to tell if hot
#endif
  double movement=0.0;
  if (way<0) {
    // down
    movement = value-floor(value);
    if (feasible) {
#ifdef COIN_DEVELOP
      hist.status_='D';
#endif
      movement = CoinMax(movement,MINIMUM_MOVEMENT);
      //printf("(down change %g value down %g ",change,movement);
      incrementNumberTimesDown();
      addToSumDownChange(1.0e-30+movement);
      addToSumDownDecrease(data.intDecrease_);
#if TYPE2==0
      addToSumDownCost(change/(1.0e-30+movement));
      setDownDynamicPseudoCost(sumDownCost()/(double) numberTimesDown());
#elif TYPE2==1
      addToSumDownCost(change);
      setDownDynamicPseudoCost(sumDownCost()/sumDownChange());
#elif TYPE2==2
      addToSumDownCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+movement));
      setDownDynamicPseudoCost(sumDownCost()*(TYPERATIO/sumDownChange()+(1.0-TYPERATIO)/(double) numberTimesDown()));
#endif
    } else {
#ifdef COIN_DEVELOP
      hist.status_='d';
#endif
      //printf("(down infeasible value down %g ",change,movement);
      incrementNumberTimesDownInfeasible();
#if INFEAS==2
      double distanceToCutoff=0.0;
      double objectiveValue = model->getCurrentMinimizationObjValue();
      distanceToCutoff =  model->getCutoff()  - originalValue;
      if (distanceToCutoff<1.0e20) 
	change = distanceToCutoff*2.0;
      else 
	change = downDynamicPseudoCost()*movement*10.0; 
      change = CoinMax(1.0e-12*(1.0+fabs(originalValue)),change);
      incrementNumberTimesDown();
      addToSumDownChange(1.0e-30+movement);
      addToSumDownDecrease(data.intDecrease_);
#if TYPE2==0
      addToSumDownCost(change/(1.0e-30+movement));
      setDownDynamicPseudoCost(sumDownCost()/(double) numberTimesDown());
#elif TYPE2==1
      addToSumDownCost(change);
      setDownDynamicPseudoCost(sumDownCost()/sumDownChange());
#elif TYPE2==2
      addToSumDownCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+movement));
      setDownDynamicPseudoCost(sumDownCost()*(TYPERATIO/sumDownChange()+(1.0-TYPERATIO)/(double) numberTimesDown()));
#endif
#endif
    }
#if INFEAS==1
    double sum = sumDownCost_;
    int number = numberTimesDown_;
    double originalValue = data.originalObjective_;
    assert (originalValue!=COIN_DBL_MAX);
    double distanceToCutoff =  data.cutoff_  - originalValue;
    if (distanceToCutoff>1.0e20) 
      distanceToCutoff=10.0+fabs(originalValue);
    sum += numberTimesDownInfeasible_*CoinMax(distanceToCutoff,1.0e-12*(1.0+fabs(originalValue)));
    number += numberTimesDownInfeasible_;
    setDownDynamicPseudoCost(sum/(double) number);
#endif
  } else {
    // up
    movement = ceil(value)-value;
    if (feasible) {
#ifdef COIN_DEVELOP
      hist.status_='U';
#endif
      movement = CoinMax(movement,MINIMUM_MOVEMENT);
      //printf("(up change %g value down %g ",change,movement);
      incrementNumberTimesUp();
      addToSumUpChange(1.0e-30+movement);
      addToSumUpDecrease(data.intDecrease_);
#if TYPE2==0
      addToSumUpCost(change/(1.0e-30+movement));
      setUpDynamicPseudoCost(sumUpCost()/(double) numberTimesUp());
#elif TYPE2==1
      addToSumUpCost(change);
      setUpDynamicPseudoCost(sumUpCost()/sumUpChange());
#elif TYPE2==2
      addToSumUpCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+movement));
      setUpDynamicPseudoCost(sumUpCost()*(TYPERATIO/sumUpChange()+(1.0-TYPERATIO)/(double) numberTimesUp()));
#endif
    } else {
#ifdef COIN_DEVELOP
      hist.status_='u';
#endif
      //printf("(up infeasible value down %g ",change,movement);
      incrementNumberTimesUpInfeasible();
#if INFEAS==2
      double distanceToCutoff=0.0;
      double objectiveValue = model->getCurrentMinimizationObjValue();
      distanceToCutoff =  model->getCutoff()  - originalValue;
      if (distanceToCutoff<1.0e20) 
	change = distanceToCutoff*2.0;
      else 
	change = upDynamicPseudoCost()*movement*10.0; 
      change = CoinMax(1.0e-12*(1.0+fabs(originalValue)),change);
      incrementNumberTimesUp();
      addToSumUpChange(1.0e-30+movement);
      addToSumUpDecrease(data.intDecrease_);
#if TYPE2==0
      addToSumUpCost(change/(1.0e-30+movement));
      setUpDynamicPseudoCost(sumUpCost()/(double) numberTimesUp());
#elif TYPE2==1
      addToSumUpCost(change);
      setUpDynamicPseudoCost(sumUpCost()/sumUpChange());
#elif TYPE2==2
      addToSumUpCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+movement));
      setUpDynamicPseudoCost(sumUpCost()*(TYPERATIO/sumUpChange()+(1.0-TYPERATIO)/(double) numberTimesUp()));
#endif
#endif
    }
#if INFEAS==1
    double sum = sumUpCost_;
    int number = numberTimesUp_;
    double originalValue = data.originalObjective_;
    assert (originalValue!=COIN_DBL_MAX);
    double distanceToCutoff =  data.cutoff_  - originalValue;
    if (distanceToCutoff>1.0e20) 
      distanceToCutoff=10.0+fabs(originalValue);
    sum += numberTimesUpInfeasible_*CoinMax(distanceToCutoff,1.0e-12*(1.0+fabs(originalValue)));
    number += numberTimesUpInfeasible_;
    setUpDynamicPseudoCost(sum/(double) number);
#endif
  }
  if (data.way_<0)
    assert (numberTimesDown_+numberTimesDownInfeasible_>0);
  else
    assert (numberTimesUp_+numberTimesUpInfeasible_>0);
  assert (downDynamicPseudoCost_>=0.0&&downDynamicPseudoCost_<1.0e100);
  downDynamicPseudoCost_ = CoinMax(1.0e-10,downDynamicPseudoCost_);
  assert (upDynamicPseudoCost_>=0.0&&upDynamicPseudoCost_<1.0e100);
  upDynamicPseudoCost_ = CoinMax(1.0e-10,upDynamicPseudoCost_);
#ifdef COIN_DEVELOP
  hist.sequence_=columnNumber_;
  hist.numberUp_=numberTimesUp_;
  hist.numberUpInf_=numberTimesUpInfeasible_;
  hist.sumUp_=sumUpCost_;
  hist.upEst_=change;
  hist.numberDown_=numberTimesDown_;
  hist.numberDownInf_=numberTimesDownInfeasible_;
  hist.sumDown_=sumDownCost_;
  hist.downEst_=movement;
  addRecord(hist);
#endif
  //print(1,0.5);
  assert (downDynamicPseudoCost_>1.0e-40&&upDynamicPseudoCost_>1.0e-40);
}
// Updates stuff like pseudocosts after mini branch and bound
void 
CbcSimpleIntegerDynamicPseudoCost::updateAfterMini(int numberDown,int numberDownInfeasible,
						   double sumDown, int numberUp,
						   int numberUpInfeasible,double sumUp)
{
#ifdef CBC_INSTRUMENT
  int difference = numberDown-numberTimesDown_;
  difference += numberUp-numberTimesUp_;
  numberTimesInfeasible_ += 2*difference;
#endif
  numberTimesDown_ = numberDown;
  numberTimesDownInfeasible_ = numberDownInfeasible;
  sumDownCost_ = sumDown;
  numberTimesUp_ = numberUp;
  numberTimesUpInfeasible_ = numberUpInfeasible;
  sumUpCost_ = sumUp;
  if (numberTimesDown_+numberTimesDownInfeasible_>0) {
    setDownDynamicPseudoCost(sumDownCost_/(double) (numberTimesDown_+numberTimesDownInfeasible_));
    assert (downDynamicPseudoCost_>0.0&&downDynamicPseudoCost_<1.0e50);
  }
  if (numberTimesUp_+numberTimesUpInfeasible_>0) {
    setUpDynamicPseudoCost(sumUpCost_/(double) (numberTimesUp_+numberTimesUpInfeasible_));
    assert (upDynamicPseudoCost_>0.0&&upDynamicPseudoCost_<1.0e50);
  }
  assert (downDynamicPseudoCost_>1.0e-40&&upDynamicPseudoCost_>1.0e-40);
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
    int n=0;
#ifdef CBC_INSTRUMENT
    n=numberTimesInfeasible_;
#endif
    printf("%d down %d times (%d inf) mean %g  up %d times (%d inf) mean %g - pseudocosts %g %g - inftimes %d\n",
           columnNumber_,
           numberTimesDown_,numberTimesDownInfeasible_,meanDown,
           numberTimesUp_,numberTimesUpInfeasible_,meanUp,downDynamicPseudoCost_,upDynamicPseudoCost_,n);
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
    distanceToCutoff = CoinMax(distanceToCutoff,1.0e-12*(1.0+fabs(objectiveValue)));
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
// Does part of work for constructor
void 
CbcDynamicPseudoCostBranchingObject::fillPart (int variable,
					   int way , double value, 
					   CbcSimpleIntegerDynamicPseudoCost * object) 
{
  CbcIntegerBranchingObject::fillPart(variable,way,value);
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
    info.upMovement = object_->upDynamicPseudoCost()*(ceil(value_)-value_);
    info.downMovement = object_->downDynamicPseudoCost()*(value_-floor(value_));
    info.numIntInfeasUp  -= (int) (object_->sumUpDecrease()/
                                   (1.0e-12+(double) object_->numberTimesUp()));
    info.numIntInfeasUp = CoinMax(info.numIntInfeasUp,0);
    info.numObjInfeasUp = 0;
    info.finishedUp = false;
    info.numItersUp = 0;
    info.numIntInfeasDown  -= (int) (object_->sumDownDecrease()/
                                   (1.0e-12+(double) object_->numberTimesDown()));
    info.numIntInfeasDown = CoinMax(info.numIntInfeasDown,0);
    info.numObjInfeasDown = 0;
    info.finishedDown = false;
    info.numItersDown = 0;
    info.fix =0;
  if (object_->numberTimesUp()<object_->numberBeforeTrust()||
      object_->numberTimesDown()<object_->numberBeforeTrust()) {
    return 0;
  } else {
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
#ifdef COIN_DEVELOP
  History hist;
  hist.where_='C';
  hist.status_='I';
  hist.sequence_=55555;
  hist.numberUp_=0;
  hist.numberUpInf_=0;
  hist.sumUp_=0.0;
  hist.upEst_=0.0;
  hist.numberDown_=0;
  hist.numberDownInf_=0;
  hist.sumDown_=0.0;
  hist.downEst_=0.0;
  addRecord(hist);
#endif
}

/* Saves a clone of current branching object.  Can be used to update
      information on object causing branch - after branch */
void 
CbcBranchDynamicDecision::saveBranchingObject(OsiBranchingObject * object) 
{
  OsiBranchingObject * obj = object->clone();
  CbcBranchingObject * obj2 =
    dynamic_cast<CbcBranchingObject *>(obj);
  assert (obj2);
  CbcDynamicPseudoCostBranchingObject * branchingObject =
    dynamic_cast<CbcDynamicPseudoCostBranchingObject *>(obj);
#if COIN_DEVELOP>1
  if (!branchingObject)
    printf("no dynamic branching object Dynamic Decision\n");
#endif
  //object_=branchingObject;
  object_ = obj2;
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
      double movement = value-floor(value);
      movement = CoinMax(movement,MINIMUM_MOVEMENT);
      //printf("(down change %g value down %g ",change,movement);
      object->incrementNumberTimesDown();
      object->addToSumDownChange(1.0e-30+movement);
      object->addToSumDownDecrease(originalUnsatisfied-unsatisfied);
#if TYPE2==0
      object->addToSumDownCost(change/(1.0e-30+movement));
      object->setDownDynamicPseudoCost(object->sumDownCost()/(double) object->numberTimesDown());
#elif TYPE2==1
      object->addToSumDownCost(change);
      object->setDownDynamicPseudoCost(object->sumDownCost()/object->sumDownChange());
#elif TYPE2==2
      object->addToSumDownCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+movement));
      object->setDownDynamicPseudoCost(object->sumDownCost()*(TYPERATIO/object->sumDownChange()+(1.0-TYPERATIO)/(double) object->numberTimesDown()));
#endif
    } else {
      //printf("(down infeasible value down %g ",change,movement);
      object->incrementNumberTimesDownInfeasible();
#if INFEAS==2
      double distanceToCutoff=0.0;
      double objectiveValue = model->getCurrentMinimizationObjValue();
      distanceToCutoff =  model->getCutoff()  - originalValue;
      if (distanceToCutoff<1.0e20) 
	change = distanceToCutoff*2.0;
      else 
	change = object->downDynamicPseudoCost()*movement*10.0; 
      change = CoinMax(1.0e-12*(1.0+fabs(originalValue)),change);
      object->incrementNumberTimesDown();
      object->addToSumDownChange(1.0e-30+movement);
      object->addToSumDownDecrease(originalUnsatisfied-unsatisfied);
#if TYPE2==0
      object->addToSumDownCost(change/(1.0e-30+movement));
      object->setDownDynamicPseudoCost(object->sumDownCost()/(double) object->numberTimesDown());
#elif TYPE2==1
      object->addToSumDownCost(change);
      object->setDownDynamicPseudoCost(object->sumDownCost()/object->sumDownChange());
#elif TYPE2==2
      object->addToSumDownCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+movement));
      object->setDownDynamicPseudoCost(object->sumDownCost()*(TYPERATIO/object->sumDownChange()+(1.0-TYPERATIO)/(double) object->numberTimesDown()));
#endif
#endif
    }
  } else {
    // up
    if (feasible) {
      double movement = ceil(value)-value;
      movement = CoinMax(movement,MINIMUM_MOVEMENT);
      //printf("(up change %g value down %g ",change,movement);
      object->incrementNumberTimesUp();
      object->addToSumUpChange(1.0e-30+movement);
      object->addToSumUpDecrease(unsatisfied-originalUnsatisfied);
#if TYPE2==0
      object->addToSumUpCost(change/(1.0e-30+movement));
      object->setUpDynamicPseudoCost(object->sumUpCost()/(double) object->numberTimesUp());
#elif TYPE2==1
      object->addToSumUpCost(change);
      object->setUpDynamicPseudoCost(object->sumUpCost()/object->sumUpChange());
#elif TYPE2==2
      object->addToSumUpCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+movement));
      object->setUpDynamicPseudoCost(object->sumUpCost()*(TYPERATIO/object->sumUpChange()+(1.0-TYPERATIO)/(double) object->numberTimesUp()));
#endif
    } else {
      //printf("(up infeasible value down %g ",change,movement);
      object->incrementNumberTimesUpInfeasible();
#if INFEAS==2
      double distanceToCutoff=0.0;
      double objectiveValue = model->getCurrentMinimizationObjValue();
      distanceToCutoff =  model->getCutoff()  - originalValue;
      if (distanceToCutoff<1.0e20) 
	change = distanceToCutoff*2.0;
      else 
	change = object->upDynamicPseudoCost()*movement*10.0; 
      change = CoinMax(1.0e-12*(1.0+fabs(originalValue)),change);
      object->incrementNumberTimesUp();
      object->addToSumUpChange(1.0e-30+movement);
      object->addToSumUpDecrease(unsatisfied-originalUnsatisfied);
#if TYPE2==0
      object->addToSumUpCost(change/(1.0e-30+movement));
      object->setUpDynamicPseudoCost(object->sumUpCost()/(double) object->numberTimesUp());
#elif TYPE2==1
      object->addToSumUpCost(change);
      object->setUpDynamicPseudoCost(object->sumUpCost()/object->sumUpChange());
#elif TYPE2==2
      object->addToSumUpCost(change*TYPERATIO+(1.0-TYPERATIO)*change/(1.0e-30+movement));
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
  int stateOfSearch = model->stateOfSearch()%10;
  int betterWay=0;
  double value=0.0;
  if (!bestObject_) {
    bestCriterion_=-1.0;
    bestNumberUp_=COIN_INT_MAX;
    bestNumberDown_=COIN_INT_MAX;
  }
  if (stateOfSearch<=2) {
    //#define TRY_STUFF 1
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
    // use pseudo shadow prices modified by locks
    // test testosi
#if 1
    double objectiveValue = model->getCurrentMinimizationObjValue();
    double distanceToCutoff =  model->getCutoff()  - objectiveValue;
    if (distanceToCutoff<1.0e20) 
      distanceToCutoff *= 10.0;
    else 
      distanceToCutoff = 1.0e2 + fabs(objectiveValue);
    distanceToCutoff = CoinMax(distanceToCutoff,1.0e-12*(1.0+fabs(objectiveValue)));
    double continuousObjective = model->getContinuousObjective();
    double distanceToCutoffC =  model->getCutoff()  - continuousObjective;
    if (distanceToCutoffC>1.0e20) 
      distanceToCutoffC = 1.0e2 + fabs(objectiveValue);
    distanceToCutoffC = CoinMax(distanceToCutoffC,1.0e-12*(1.0+fabs(objectiveValue)));
    int numberInfC = model->getContinuousInfeasibilities();
    double perInf = distanceToCutoffC/((double) numberInfC);
    assert (perInf>0.0);
    //int numberIntegers = model->numberIntegers();
    changeDown += perInf * numInfDown;
    changeUp += perInf * numInfUp;
#endif
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
    //#define TRY_STUFF 2
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
    //maxValue = CoinMin(maxValue,minValue*2.0);
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
#ifdef COIN_DEVELOP
  History hist;
  {
    CbcDynamicPseudoCostBranchingObject * branchingObject =
      dynamic_cast<CbcDynamicPseudoCostBranchingObject *>(thisOne);
    if (branchingObject) {
      CbcSimpleIntegerDynamicPseudoCost *  object = branchingObject->object();
      assert (object);
      hist.where_='C';
      hist.status_=' ';
      hist.sequence_=object->columnNumber();
      hist.numberUp_=object->numberTimesUp();
      hist.numberUpInf_=numInfUp;
      hist.sumUp_=object->sumUpCost();
      hist.upEst_=changeUp;
      hist.numberDown_=object->numberTimesDown();
      hist.numberDownInf_=numInfDown;
      hist.sumDown_=object->sumDownCost();
      hist.downEst_=changeDown;
    }
  }
#endif
  if (betterWay) {
#ifdef COIN_DEVELOP
    hist.status_='B';
#endif
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
#ifdef COIN_DEVELOP
  addRecord(hist);
#endif
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
#ifdef COIN_DEVELOP
void printHistory(const char * file)
{
  if (!numberHistory)
    return;
  FILE * fp = fopen(file,"w");
  assert(fp);
  unsigned short numberIntegers=0;
  int i;
  for (i=0;i<numberHistory;i++) {
    if (history[i].where_!='C'||history[i].status_!='I') 
      numberIntegers = CoinMax(numberIntegers,history[i].sequence_);
  }
  numberIntegers++;
  for (int iC=0;iC<numberIntegers;iC++) {
    int n=0;
    for (i=0;i<numberHistory;i++) {
      if (history[i].sequence_==iC) {
	if (!n)
	  fprintf(fp,"XXX %d\n",iC);
	n++;
	fprintf(fp,"%c%c up %8d %8d %12.5f %12.5f down  %8d %8d %12.5f %12.5f\n",
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
