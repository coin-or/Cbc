// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/17/2009 - carved out of CbcBranchDynamic

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
#include "CbcSimpleIntegerDynamicPseudoCost.hpp"
#include "OsiClpSolverInterface.hpp"
#if 0
#undef TYPE2
#undef INFEAS
#undef MOD_SHADOW
#undef WEIGHT_PRODUCT
#define TYPE2 0 // 0 can be 1 or 2 (1 or 2 errors)
//#define FUNNY_BRANCHING
//#define FUNNY_BRANCHING2
#define INFEAS 1 // 1 can be 0 (not much difference) or 2 (marginally better?)
#define INFEAS_MULTIPLIER 1.5 // default 1.5
#define MOD_SHADOW 1 // 1 can be 0 or 2
#define WEIGHT_PRODUCT // on
#endif
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
History *history = NULL;
int numberHistory = 0;
int maxHistory = 0;
bool getHistoryStatistics_ = true;
static void increaseHistory()
{
  if (numberHistory == maxHistory) {
    maxHistory = 100 + (3 * maxHistory) / 2;
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

/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
CbcSimpleIntegerDynamicPseudoCost::CbcSimpleIntegerDynamicPseudoCost()
  : CbcSimpleInteger()
  , downDynamicPseudoCost_(1.0e-5)
  , upDynamicPseudoCost_(1.0e-5)
  , upDownSeparator_(-1.0)
  , sumDownCost_(0.0)
  , sumUpCost_(0.0)
  , sumDownChange_(0.0)
  , sumUpChange_(0.0)
  , downShadowPrice_(0.0)
  , upShadowPrice_(0.0)
  , sumDownDecrease_(0.0)
  , sumUpDecrease_(0.0)
#if CBC_DYNAMIC_EXPERIMENT == 0
  , lastDownCost_(0.0)
  , lastUpCost_(0.0)
#elif CBC_DYNAMIC_EXPERIMENT > 0
  , downOver_(1.0)
  , downUnder_(1.0)
  , upOver_(1.0)
  , upUnder_(1.0)
#endif
  , lastDownDecrease_(0)
  , lastUpDecrease_(0)
  , numberTimesDown_(0)
  , numberTimesUp_(0)
  , numberTimesDownInfeasible_(0)
  , numberTimesUpInfeasible_(0)
  , numberBeforeTrust_(0)
#if CBC_DYNAMIC_EXPERIMENT == 0
  , numberTimesDownLocalFixed_(0)
  , numberTimesUpLocalFixed_(0)
  , numberTimesDownTotalFixed_(0.0)
  , numberTimesUpTotalFixed_(0.0)
  , numberTimesProbingTotal_(0)
  , method_(0)
#endif
{
}

/** Useful constructor

  Loads dynamic upper & lower bounds for the specified variable.
*/
CbcSimpleIntegerDynamicPseudoCost::CbcSimpleIntegerDynamicPseudoCost(CbcModel *model,
  int iColumn, double breakEven)
  : CbcSimpleInteger(model, iColumn, breakEven)
  , upDownSeparator_(-1.0)
  , sumDownCost_(0.0)
  , sumUpCost_(0.0)
  , sumDownChange_(0.0)
  , sumUpChange_(0.0)
  , downShadowPrice_(0.0)
  , upShadowPrice_(0.0)
  , sumDownDecrease_(0.0)
  , sumUpDecrease_(0.0)
#if CBC_DYNAMIC_EXPERIMENT == 0
  , lastDownCost_(0.0)
  , lastUpCost_(0.0)
#elif CBC_DYNAMIC_EXPERIMENT > 0
  , downOver_(1.0)
  , downUnder_(1.0)
  , upOver_(1.0)
  , upUnder_(1.0)
#endif
  , lastDownDecrease_(0)
  , lastUpDecrease_(0)
  , numberTimesDown_(0)
  , numberTimesUp_(0)
  , numberTimesDownInfeasible_(0)
  , numberTimesUpInfeasible_(0)
  , numberBeforeTrust_(0)
#if CBC_DYNAMIC_EXPERIMENT == 0
  , numberTimesDownLocalFixed_(0)
  , numberTimesUpLocalFixed_(0)
  , numberTimesDownTotalFixed_(0.0)
  , numberTimesUpTotalFixed_(0.0)
  , numberTimesProbingTotal_(0)
  , method_(0)
#endif
{
  const double *cost = model->getObjCoefficients();
  double costValue = std::max(1.0e-5, fabs(cost[iColumn]));
  // treat as if will cost what it says up
  upDynamicPseudoCost_ = costValue;
  // and balance at breakeven
  downDynamicPseudoCost_ = ((1.0 - breakEven_) * upDynamicPseudoCost_) / breakEven_;
  // so initial will have some effect
  sumUpCost_ = 2.0 * upDynamicPseudoCost_;
  sumUpChange_ = 2.0;
  numberTimesUp_ = 2;
  sumDownCost_ = 2.0 * downDynamicPseudoCost_;
  sumDownChange_ = 2.0;
  numberTimesDown_ = 2;
#if TYPE2 == 0
  // No
  sumUpCost_ = 0.0;
  sumUpChange_ = 0.0;
  numberTimesUp_ = 0;
  sumDownCost_ = 0.0;
  sumDownChange_ = 0.0;
  numberTimesDown_ = 0;
#else
  sumUpCost_ = 1.0 * upDynamicPseudoCost_;
  sumUpChange_ = 1.0;
  numberTimesUp_ = 1;
  sumDownCost_ = 1.0 * downDynamicPseudoCost_;
  sumDownChange_ = 1.0;
  numberTimesDown_ = 1;
#endif
}

/** Useful constructor

  Loads dynamic upper & lower bounds for the specified variable.
*/
CbcSimpleIntegerDynamicPseudoCost::CbcSimpleIntegerDynamicPseudoCost(CbcModel *model,
  int iColumn, double downDynamicPseudoCost,
  double upDynamicPseudoCost)
  : CbcSimpleInteger(model, iColumn)
  , upDownSeparator_(-1.0)
  , sumDownCost_(0.0)
  , sumUpCost_(0.0)
  , sumDownChange_(0.0)
  , sumUpChange_(0.0)
  , downShadowPrice_(0.0)
  , upShadowPrice_(0.0)
  , sumDownDecrease_(0.0)
  , sumUpDecrease_(0.0)
#if CBC_DYNAMIC_EXPERIMENT == 0
  , lastDownCost_(0.0)
  , lastUpCost_(0.0)
#elif CBC_DYNAMIC_EXPERIMENT > 0
  , downOver_(1.0)
  , downUnder_(1.0)
  , upOver_(1.0)
  , upUnder_(1.0)
#endif
  , lastDownDecrease_(0)
  , lastUpDecrease_(0)
  , numberTimesDown_(0)
  , numberTimesUp_(0)
  , numberTimesDownInfeasible_(0)
  , numberTimesUpInfeasible_(0)
  , numberBeforeTrust_(0)
#if CBC_DYNAMIC_EXPERIMENT == 0
  , numberTimesDownLocalFixed_(0)
  , numberTimesUpLocalFixed_(0)
  , numberTimesDownTotalFixed_(0.0)
  , numberTimesUpTotalFixed_(0.0)
  , numberTimesProbingTotal_(0)
  , method_(0)
#endif
{
  downDynamicPseudoCost_ = downDynamicPseudoCost;
  upDynamicPseudoCost_ = upDynamicPseudoCost;
  breakEven_ = upDynamicPseudoCost_ / (upDynamicPseudoCost_ + downDynamicPseudoCost_);
  // so initial will have some effect
  sumUpCost_ = 2.0 * upDynamicPseudoCost_;
  sumUpChange_ = 2.0;
  numberTimesUp_ = 2;
  sumDownCost_ = 2.0 * downDynamicPseudoCost_;
  sumDownChange_ = 2.0;
  numberTimesDown_ = 2;
#if TYPE2 == 0
  // No
  sumUpCost_ = 0.0;
  sumUpChange_ = 0.0;
  numberTimesUp_ = 0;
  sumDownCost_ = 0.0;
  sumDownChange_ = 0.0;
  numberTimesDown_ = 0;
  sumUpCost_ = 1.0e-4 * upDynamicPseudoCost_;
  sumDownCost_ = 1.0e-4 * downDynamicPseudoCost_;
#else
  sumUpCost_ = 1.0 * upDynamicPseudoCost_;
  sumUpChange_ = 1.0;
  numberTimesUp_ = 1;
  sumDownCost_ = 1.0 * downDynamicPseudoCost_;
  sumDownChange_ = 1.0;
  numberTimesDown_ = 1;
#endif
}
/** Useful constructor

  Loads dynamic upper & lower bounds for the specified variable.
*/
CbcSimpleIntegerDynamicPseudoCost::CbcSimpleIntegerDynamicPseudoCost(CbcModel *model,
  int /*dummy*/,
  int iColumn, double downDynamicPseudoCost,
  double upDynamicPseudoCost)
{
  CbcSimpleIntegerDynamicPseudoCost(model, iColumn, downDynamicPseudoCost, upDynamicPseudoCost);
}

// Copy constructor
CbcSimpleIntegerDynamicPseudoCost::CbcSimpleIntegerDynamicPseudoCost(const CbcSimpleIntegerDynamicPseudoCost &rhs)
  : CbcSimpleInteger(rhs)
  , downDynamicPseudoCost_(rhs.downDynamicPseudoCost_)
  , upDynamicPseudoCost_(rhs.upDynamicPseudoCost_)
  , upDownSeparator_(rhs.upDownSeparator_)
  , sumDownCost_(rhs.sumDownCost_)
  , sumUpCost_(rhs.sumUpCost_)
  , sumDownChange_(rhs.sumDownChange_)
  , sumUpChange_(rhs.sumUpChange_)
  , downShadowPrice_(rhs.downShadowPrice_)
  , upShadowPrice_(rhs.upShadowPrice_)
  , sumDownDecrease_(rhs.sumDownDecrease_)
  , sumUpDecrease_(rhs.sumUpDecrease_)
#if CBC_DYNAMIC_EXPERIMENT == 0
  , lastDownCost_(rhs.lastDownCost_)
  , lastUpCost_(rhs.lastUpCost_)
#elif CBC_DYNAMIC_EXPERIMENT > 0
  , downOver_(rhs.downOver_)
  , downUnder_(rhs.downUnder_)
  , upOver_(rhs.upOver_)
  , upUnder_(rhs.upUnder_)
#endif
  , lastDownDecrease_(rhs.lastDownDecrease_)
  , lastUpDecrease_(rhs.lastUpDecrease_)
  , numberTimesDown_(rhs.numberTimesDown_)
  , numberTimesUp_(rhs.numberTimesUp_)
  , numberTimesDownInfeasible_(rhs.numberTimesDownInfeasible_)
  , numberTimesUpInfeasible_(rhs.numberTimesUpInfeasible_)
  , numberBeforeTrust_(rhs.numberBeforeTrust_)
#if CBC_DYNAMIC_EXPERIMENT == 0
  , numberTimesDownLocalFixed_(rhs.numberTimesDownLocalFixed_)
  , numberTimesUpLocalFixed_(rhs.numberTimesUpLocalFixed_)
  , numberTimesDownTotalFixed_(rhs.numberTimesDownTotalFixed_)
  , numberTimesUpTotalFixed_(rhs.numberTimesUpTotalFixed_)
  , numberTimesProbingTotal_(rhs.numberTimesProbingTotal_)
  , method_(rhs.method_)
#endif
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
CbcSimpleIntegerDynamicPseudoCost::operator=(const CbcSimpleIntegerDynamicPseudoCost &rhs)
{
  if (this != &rhs) {
    CbcSimpleInteger::operator=(rhs);
    downDynamicPseudoCost_ = rhs.downDynamicPseudoCost_;
    upDynamicPseudoCost_ = rhs.upDynamicPseudoCost_;
    upDownSeparator_ = rhs.upDownSeparator_;
    sumDownCost_ = rhs.sumDownCost_;
    sumUpCost_ = rhs.sumUpCost_;
    sumDownChange_ = rhs.sumDownChange_;
    sumUpChange_ = rhs.sumUpChange_;
    downShadowPrice_ = rhs.downShadowPrice_;
    upShadowPrice_ = rhs.upShadowPrice_;
    sumDownDecrease_ = rhs.sumDownDecrease_;
    sumUpDecrease_ = rhs.sumUpDecrease_;
#if CBC_DYNAMIC_EXPERIMENT == 0
    lastDownCost_ = rhs.lastDownCost_;
    lastUpCost_ = rhs.lastUpCost_;
#elif CBC_DYNAMIC_EXPERIMENT > 0
    downOver_ = rhs.downOver_;
    downUnder_ = rhs.downUnder_;
    upOver_ = rhs.upOver_;
    upUnder_ = rhs.upUnder_;
#endif
    lastDownDecrease_ = rhs.lastDownDecrease_;
    lastUpDecrease_ = rhs.lastUpDecrease_;
    numberTimesDown_ = rhs.numberTimesDown_;
    numberTimesUp_ = rhs.numberTimesUp_;
    numberTimesDownInfeasible_ = rhs.numberTimesDownInfeasible_;
    numberTimesUpInfeasible_ = rhs.numberTimesUpInfeasible_;
    numberBeforeTrust_ = rhs.numberBeforeTrust_;
#if CBC_DYNAMIC_EXPERIMENT == 0
    numberTimesDownLocalFixed_ = rhs.numberTimesDownLocalFixed_;
    numberTimesUpLocalFixed_ = rhs.numberTimesUpLocalFixed_;
    numberTimesDownTotalFixed_ = rhs.numberTimesDownTotalFixed_;
    numberTimesUpTotalFixed_ = rhs.numberTimesUpTotalFixed_;
    numberTimesProbingTotal_ = rhs.numberTimesProbingTotal_;
    method_ = rhs.method_;
#endif
  }
  return *this;
}

// Destructor
CbcSimpleIntegerDynamicPseudoCost::~CbcSimpleIntegerDynamicPseudoCost()
{
}
// Copy some information i.e. just variable stuff
void CbcSimpleIntegerDynamicPseudoCost::copySome(const CbcSimpleIntegerDynamicPseudoCost *otherObject)
{
  downDynamicPseudoCost_ = otherObject->downDynamicPseudoCost_;
  upDynamicPseudoCost_ = otherObject->upDynamicPseudoCost_;
  sumDownCost_ = otherObject->sumDownCost_;
  sumUpCost_ = otherObject->sumUpCost_;
  sumDownChange_ = otherObject->sumDownChange_;
  sumUpChange_ = otherObject->sumUpChange_;
  downShadowPrice_ = otherObject->downShadowPrice_;
  upShadowPrice_ = otherObject->upShadowPrice_;
  sumDownDecrease_ = otherObject->sumDownDecrease_;
  sumUpDecrease_ = otherObject->sumUpDecrease_;
#if CBC_DYNAMIC_EXPERIMENT == 0
  lastDownCost_ = otherObject->lastDownCost_;
  lastUpCost_ = otherObject->lastUpCost_;
#elif CBC_DYNAMIC_EXPERIMENT > 0
  downOver_ = otherObject->downOver_;
  downUnder_ = otherObject->downUnder_;
  upOver_ = otherObject->upOver_;
  upUnder_ = otherObject->upUnder_;
#endif
  lastDownDecrease_ = otherObject->lastDownDecrease_;
  lastUpDecrease_ = otherObject->lastUpDecrease_;
  numberTimesDown_ = otherObject->numberTimesDown_;
  numberTimesUp_ = otherObject->numberTimesUp_;
  numberTimesDownInfeasible_ = otherObject->numberTimesDownInfeasible_;
  numberTimesUpInfeasible_ = otherObject->numberTimesUpInfeasible_;
#if CBC_DYNAMIC_EXPERIMENT == 0
  numberTimesDownLocalFixed_ = otherObject->numberTimesDownLocalFixed_;
  numberTimesUpLocalFixed_ = otherObject->numberTimesUpLocalFixed_;
  numberTimesDownTotalFixed_ = otherObject->numberTimesDownTotalFixed_;
  numberTimesUpTotalFixed_ = otherObject->numberTimesUpTotalFixed_;
  numberTimesProbingTotal_ = otherObject->numberTimesProbingTotal_;
#endif
}
// Updates stuff like pseudocosts before threads
void CbcSimpleIntegerDynamicPseudoCost::updateBefore(const OsiObject *rhs)
{
#ifndef NDEBUG
  const CbcSimpleIntegerDynamicPseudoCost *rhsObject = dynamic_cast< const CbcSimpleIntegerDynamicPseudoCost * >(rhs);
  assert(rhsObject);
#else
  const CbcSimpleIntegerDynamicPseudoCost *rhsObject = static_cast< const CbcSimpleIntegerDynamicPseudoCost * >(rhs);
#endif
  copySome(rhsObject);
}
// Updates stuff like pseudocosts after threads finished
void CbcSimpleIntegerDynamicPseudoCost::updateAfter(const OsiObject *rhs, const OsiObject *baseObjectX)
{
#ifndef NDEBUG
  const CbcSimpleIntegerDynamicPseudoCost *rhsObject = dynamic_cast< const CbcSimpleIntegerDynamicPseudoCost * >(rhs);
  assert(rhsObject);
  const CbcSimpleIntegerDynamicPseudoCost *baseObject = dynamic_cast< const CbcSimpleIntegerDynamicPseudoCost * >(baseObjectX);
  assert(baseObject);
#else
  const CbcSimpleIntegerDynamicPseudoCost *rhsObject = static_cast< const CbcSimpleIntegerDynamicPseudoCost * >(rhs);
  const CbcSimpleIntegerDynamicPseudoCost *baseObject = static_cast< const CbcSimpleIntegerDynamicPseudoCost * >(baseObjectX);
#endif
  // compute current
  double sumDown = downDynamicPseudoCost_ * numberTimesDown_;
  sumDown -= baseObject->downDynamicPseudoCost_ * baseObject->numberTimesDown_;
  sumDown = std::max(sumDown, 0.0);
  sumDown += rhsObject->downDynamicPseudoCost_ * rhsObject->numberTimesDown_;
  assert(rhsObject->numberTimesDown_ >= baseObject->numberTimesDown_);
  assert(rhsObject->numberTimesDownInfeasible_ >= baseObject->numberTimesDownInfeasible_);
  assert(rhsObject->sumDownCost_ >= baseObject->sumDownCost_ - 1.0e-4);
  double sumUp = upDynamicPseudoCost_ * numberTimesUp_;
  sumUp -= baseObject->upDynamicPseudoCost_ * baseObject->numberTimesUp_;
  sumUp = std::max(sumUp, 0.0);
  sumUp += rhsObject->upDynamicPseudoCost_ * rhsObject->numberTimesUp_;
  assert(rhsObject->numberTimesUp_ >= baseObject->numberTimesUp_);
  assert(rhsObject->numberTimesUpInfeasible_ >= baseObject->numberTimesUpInfeasible_);
  assert(rhsObject->sumUpCost_ >= baseObject->sumUpCost_ - 1.0e-4);
  sumDownCost_ += rhsObject->sumDownCost_ - baseObject->sumDownCost_;
  sumUpCost_ += rhsObject->sumUpCost_ - baseObject->sumUpCost_;
  sumDownChange_ += rhsObject->sumDownChange_ - baseObject->sumDownChange_;
  sumUpChange_ += rhsObject->sumUpChange_ - baseObject->sumUpChange_;
  downShadowPrice_ = 0.0;
  upShadowPrice_ = 0.0;
  sumDownDecrease_ += rhsObject->sumDownDecrease_ - baseObject->sumDownDecrease_;
  sumUpDecrease_ += rhsObject->sumUpDecrease_ - baseObject->sumUpDecrease_;
#if CBC_DYNAMIC_EXPERIMENT == 0
  lastDownCost_ += rhsObject->lastDownCost_ - baseObject->lastDownCost_;
  lastUpCost_ += rhsObject->lastUpCost_ - baseObject->lastUpCost_;
#elif CBC_DYNAMIC_EXPERIMENT > 0
  // Bad - think - hope only when threaded
  downOver_ = 0.1*rhsObject->downOver_ + 0.9*baseObject->downOver_;
  downUnder_ = 0.1*rhsObject->downUnder_ + 0.9*baseObject->downUnder_;
  upOver_ = 0.1*rhsObject->upOver_ + 0.9*baseObject->upOver_;
  upUnder_ = 0.1*rhsObject->upUnder_ + 0.9*baseObject->upUnder_;
#endif
  lastDownDecrease_ += rhsObject->lastDownDecrease_ - baseObject->lastDownDecrease_;
  lastUpDecrease_ += rhsObject->lastUpDecrease_ - baseObject->lastUpDecrease_;
  numberTimesDown_ += rhsObject->numberTimesDown_ - baseObject->numberTimesDown_;
  numberTimesUp_ += rhsObject->numberTimesUp_ - baseObject->numberTimesUp_;
  numberTimesDownInfeasible_ += rhsObject->numberTimesDownInfeasible_ - baseObject->numberTimesDownInfeasible_;
  numberTimesUpInfeasible_ += rhsObject->numberTimesUpInfeasible_ - baseObject->numberTimesUpInfeasible_;
#if CBC_DYNAMIC_EXPERIMENT == 0
  numberTimesDownLocalFixed_ += rhsObject->numberTimesDownLocalFixed_ - baseObject->numberTimesDownLocalFixed_;
  numberTimesUpLocalFixed_ += rhsObject->numberTimesUpLocalFixed_ - baseObject->numberTimesUpLocalFixed_;
  numberTimesDownTotalFixed_ += rhsObject->numberTimesDownTotalFixed_ - baseObject->numberTimesDownTotalFixed_;
  numberTimesUpTotalFixed_ += rhsObject->numberTimesUpTotalFixed_ - baseObject->numberTimesUpTotalFixed_;
  numberTimesProbingTotal_ += rhsObject->numberTimesProbingTotal_ - baseObject->numberTimesProbingTotal_;
#endif
  if (numberTimesDown_ > 0) {
    setDownDynamicPseudoCost(sumDown / static_cast< double >(numberTimesDown_));
  }
  if (numberTimesUp_ > 0) {
    setUpDynamicPseudoCost(sumUp / static_cast< double >(numberTimesUp_));
  }
  //printf("XX %d down %d %d %g up %d %d %g\n",columnNumber_,numberTimesDown_,numberTimesDownInfeasible_,downDynamicPseudoCost_,
  // numberTimesUp_,numberTimesUpInfeasible_,upDynamicPseudoCost_);
  assert(downDynamicPseudoCost_ > 1.0e-40 && upDynamicPseudoCost_ > 1.0e-40);
}
// Same - returns true if contents match(ish)
bool CbcSimpleIntegerDynamicPseudoCost::same(const CbcSimpleIntegerDynamicPseudoCost *otherObject) const
{
  bool okay = true;
  if (downDynamicPseudoCost_ != otherObject->downDynamicPseudoCost_)
    okay = false;
  if (upDynamicPseudoCost_ != otherObject->upDynamicPseudoCost_)
    okay = false;
  if (sumDownCost_ != otherObject->sumDownCost_)
    okay = false;
  if (sumUpCost_ != otherObject->sumUpCost_)
    okay = false;
  if (sumDownChange_ != otherObject->sumDownChange_)
    okay = false;
  if (sumUpChange_ != otherObject->sumUpChange_)
    okay = false;
  if (downShadowPrice_ != otherObject->downShadowPrice_)
    okay = false;
  if (upShadowPrice_ != otherObject->upShadowPrice_)
    okay = false;
  if (sumDownDecrease_ != otherObject->sumDownDecrease_)
    okay = false;
  if (sumUpDecrease_ != otherObject->sumUpDecrease_)
    okay = false;
#if CBC_DYNAMIC_EXPERIMENT == 0
  if (lastDownCost_ != otherObject->lastDownCost_)
    okay = false;
  if (lastUpCost_ != otherObject->lastUpCost_)
    okay = false;
#endif
  if (lastDownDecrease_ != otherObject->lastDownDecrease_)
    okay = false;
  if (lastUpDecrease_ != otherObject->lastUpDecrease_)
    okay = false;
  if (numberTimesDown_ != otherObject->numberTimesDown_)
    okay = false;
  if (numberTimesUp_ != otherObject->numberTimesUp_)
    okay = false;
  if (numberTimesDownInfeasible_ != otherObject->numberTimesDownInfeasible_)
    okay = false;
  if (numberTimesUpInfeasible_ != otherObject->numberTimesUpInfeasible_)
    okay = false;
#if CBC_DYNAMIC_EXPERIMENT == 0
  if (numberTimesDownLocalFixed_ != otherObject->numberTimesDownLocalFixed_)
    okay = false;
  if (numberTimesUpLocalFixed_ != otherObject->numberTimesUpLocalFixed_)
    okay = false;
  if (numberTimesDownTotalFixed_ != otherObject->numberTimesDownTotalFixed_)
    okay = false;
  if (numberTimesUpTotalFixed_ != otherObject->numberTimesUpTotalFixed_)
    okay = false;
  if (numberTimesProbingTotal_ != otherObject->numberTimesProbingTotal_)
    okay = false;
#endif
  return okay;
}
/* Create an OsiSolverBranch object

This returns NULL if branch not represented by bound changes
*/
OsiSolverBranch *
CbcSimpleIntegerDynamicPseudoCost::solverBranch() const
{
  OsiSolverInterface *solver = model_->solver();
  const double *solution = model_->testSolution();
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  double value = solution[columnNumber_];
  value = std::max(value, lower[columnNumber_]);
  value = std::min(value, upper[columnNumber_]);
  assert(upper[columnNumber_] > lower[columnNumber_]);
#ifndef NDEBUG
  double nearest = floor(value + 0.5);
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  assert(fabs(value - nearest) > integerTolerance);
#endif
  OsiSolverBranch *branch = new OsiSolverBranch();
  branch->addBranch(columnNumber_, value);
  return branch;
}
//#define FUNNY_BRANCHING
double
CbcSimpleIntegerDynamicPseudoCost::infeasibility(const OsiBranchingInformation *info,
  int &preferredWay) const
{
  assert(downDynamicPseudoCost_ > 1.0e-40 && upDynamicPseudoCost_ > 1.0e-40);
  const double *solution = model_->testSolution();
  const double *lower = model_->getCbcColLower();
  const double *upper = model_->getCbcColUpper();
#ifdef FUNNY_BRANCHING2
  const double *dj = model_->getCbcReducedCost();
  double djValue = dj[columnNumber_];
  lastDownDecrease_++;
  if (djValue > 1.0e-6) {
    // wants to go down
    if (true || lower[columnNumber_] > originalLower_) {
      // Lower bound active
      lastUpDecrease_++;
    }
  } else if (djValue < -1.0e-6) {
    // wants to go up
    if (true || upper[columnNumber_] < originalUpper_) {
      // Upper bound active
      lastUpDecrease_++;
    }
  }
#endif
  if (upper[columnNumber_] == lower[columnNumber_]) {
    // fixed
    preferredWay = 1;
    return 0.0;
  }
  assert(breakEven_ > 0.0 && breakEven_ < 1.0);
  /*
	  Find nearest integer, and integers above and below current value.

	  Given that we've already forced value within bounds, if
	  (current value)+(integer tolerance) > (upper bound)
	  shouldn't we declare this variable integer?
	*/

  double value = solution[columnNumber_];
  value = std::max(value, lower[columnNumber_]);
  value = std::min(value, upper[columnNumber_]);
  /*printf("%d %g %g %g %g\n",columnNumber_,value,lower[columnNumber_],
      solution[columnNumber_],upper[columnNumber_]);*/
  double nearest = floor(value + 0.5);
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double below = floor(value + integerTolerance);
  double above = below + 1.0;
  if (above > upper[columnNumber_]) {
    above = below;
    below = above - 1;
  }
#if INFEAS == 1
  /*
  Why do we inflate the distance to the cutoff by a factor of 10 for
  values that could be considered reachable? Why do we add 100 for values
  larger than 1e20?
*/
  double distanceToCutoff = 0.0;
  double objectiveValue = model_->getCurrentMinimizationObjValue();
  distanceToCutoff = model_->getCutoff() - objectiveValue;
  if (distanceToCutoff < 1.0e20)
    distanceToCutoff *= 10.0;
  else
    distanceToCutoff = 1.0e2 + fabs(objectiveValue);
  distanceToCutoff = std::max(distanceToCutoff, 1.0e-12 * (1.0 + fabs(objectiveValue)));
#endif
  double sum;
#ifndef INFEAS_MULTIPLIER
#define INFEAS_MULTIPLIER 1.5
#endif
  double number;
  double downCost = std::max(value - below, 0.0);
#if TYPE2 == 0
  sum = sumDownCost_;
  number = numberTimesDown_;
#if INFEAS == 1
  sum += INFEAS_MULTIPLIER * numberTimesDownInfeasible_ * std::max(distanceToCutoff / (downCost + 1.0e-12), sumDownCost_);
#endif
#elif TYPE2 == 1
  sum = sumDownCost_;
  number = sumDownChange_;
#if INFEAS == 1
  sum += INFEAS_MULTIPLIER * numberTimesDownInfeasible_ * std::max(distanceToCutoff / (downCost + 1.0e-12), sumDownCost_);
#endif
#elif TYPE2 == 2
  abort();
#if INFEAS == 1
  sum += INFEAS_MULTIPLIER * numberTimesDownInfeasible_ * (distanceToCutoff / (downCost + 1.0e-12));
#endif
#endif
#if MOD_SHADOW > 0
  if (!downShadowPrice_) {
    if (number > 0.0)
      downCost *= sum / number;
    else
      downCost *= downDynamicPseudoCost_;
  } else if (downShadowPrice_ > 0.0) {
    downCost *= downShadowPrice_;
  } else {
    downCost *= (downDynamicPseudoCost_ - downShadowPrice_);
  }
#else
  if (downShadowPrice_ <= 0.0) {
    if (number > 0.0)
      downCost *= sum / number;
    else
      downCost *= downDynamicPseudoCost_;
  } else {
    downCost *= downShadowPrice_;
  }
#endif
  double upCost = std::max((above - value), 0.0);
#if TYPE2 == 0
  sum = sumUpCost_;
  number = numberTimesUp_;
#if INFEAS == 1
  sum += INFEAS_MULTIPLIER * numberTimesUpInfeasible_ * std::max(distanceToCutoff / (upCost + 1.0e-12), sumUpCost_);
#endif
#elif TYPE2 == 1
  sum = sumUpCost_;
  number = sumUpChange_;
#if INFEAS == 1
  sum += INFEAS_MULTIPLIER * numberTimesUpInfeasible_ * std::max(distanceToCutoff / (upCost + 1.0e-12), sumUpCost_);
#endif
#elif TYPE2 == 1
  abort();
#if INFEAS == 1
  sum += INFEAS_MULTIPLIER * numberTimesUpInfeasible_ * (distanceToCutoff / (upCost + 1.0e-12));
#endif
#endif
#if MOD_SHADOW > 0
  if (!upShadowPrice_) {
    if (number > 0.0)
      upCost *= sum / number;
    else
      upCost *= upDynamicPseudoCost_;
  } else if (upShadowPrice_ > 0.0) {
    upCost *= upShadowPrice_;
  } else {
    upCost *= (upDynamicPseudoCost_ - upShadowPrice_);
  }
#else
  if (upShadowPrice_ <= 0.0) {
    if (number > 0.0)
      upCost *= sum / number;
    else
      upCost *= upDynamicPseudoCost_;
  } else {
    upCost *= upShadowPrice_;
  }
#endif
  if (downCost >= upCost)
    preferredWay = 1;
  else
    preferredWay = -1;
  // See if up down choice set
  if (upDownSeparator_ > 0.0) {
    preferredWay = (value - below >= upDownSeparator_) ? 1 : -1;
  }
#ifdef FUNNY_BRANCHING2
  if (fabs(value - nearest) > integerTolerance) {
    double ratio = (100.0 + lastUpDecrease_) / (100.0 + lastDownDecrease_);
    downCost *= ratio;
    upCost *= ratio;
    if ((lastUpDecrease_ % 100) == -1)
      printf("col %d total %d djtimes %d\n",
        columnNumber_, lastDownDecrease_, lastUpDecrease_);
  }
#endif
  if (preferredWay_)
    preferredWay = preferredWay_;
  if (info->hotstartSolution_) {
    double targetValue = info->hotstartSolution_[columnNumber_];
    if (value > targetValue)
      preferredWay = -1;
    else
      preferredWay = 1;
  }
  if (fabs(value - nearest) <= integerTolerance) {
    if (priority_ != -999)
      return 0.0;
    else
      return 1.0e-13;
  } else {
    int stateOfSearch = model_->stateOfSearch() % 10;
    double returnValue = 0.0;
    double minValue = std::min(downCost, upCost);
    double maxValue = std::max(downCost, upCost);
#ifdef COIN_DEVELOP
    char where;
#endif
    // was <= 10
    //if (stateOfSearch<=1||model_->currentNode()->depth()<=-10 /* was ||maxValue>0.2*distanceToCutoff*/) {
    if (stateOfSearch < 1) {
      // no solution
#ifdef COIN_DEVELOP
      where = 'i';
#endif
      returnValue = WEIGHT_BEFORE * minValue + (1.0 - WEIGHT_BEFORE) * maxValue;
      if (0) {
        double sum;
        int number;
        double downCost2 = std::max(value - below, 0.0);
        sum = sumDownCost_;
        number = numberTimesDown_;
        if (number > 0)
          downCost2 *= sum / static_cast< double >(number);
        else
          downCost2 *= downDynamicPseudoCost_;
        double upCost2 = std::max((above - value), 0.0);
        sum = sumUpCost_;
        number = numberTimesUp_;
        if (number > 0)
          upCost2 *= sum / static_cast< double >(number);
        else
          upCost2 *= upDynamicPseudoCost_;
        double minValue2 = std::min(downCost2, upCost2);
        double maxValue2 = std::max(downCost2, upCost2);
        printf("%d value %g downC %g upC %g minV %g maxV %g downC2 %g upC2 %g minV2 %g maxV2 %g\n",
          columnNumber_, value, downCost, upCost, minValue, maxValue,
          downCost2, upCost2, minValue2, maxValue2);
      }
    } else {
      // some solution
#ifdef COIN_DEVELOP
      where = 'I';
#endif
#ifndef WEIGHT_PRODUCT
      returnValue = WEIGHT_AFTER * minValue + (1.0 - WEIGHT_AFTER) * maxValue;
#else
      double minProductWeight = model_->getDblParam(CbcModel::CbcSmallChange);
      returnValue = std::max(minValue, minProductWeight) * std::max(maxValue, minProductWeight);
      //returnValue += minProductWeight*minValue;
#endif
    }
    if (numberTimesUp_ < numberBeforeTrust_ || numberTimesDown_ < numberBeforeTrust_) {
      //if (returnValue<1.0e10)
      //returnValue += 1.0e12;
      //else
      returnValue *= 1.0e3;
      if (!numberTimesUp_ && !numberTimesDown_)
        returnValue *= 1.0e10;
    }
    //if (fabs(value-0.5)<1.0e-5) {
    //returnValue = 3.0*returnValue + 0.2;
    //} else if (value>0.9) {
    //returnValue = 2.0*returnValue + 0.1;
    //}
#if CBC_DYNAMIC_EXPERIMENT == 0
    if (method_ == 1) {
#if 1
      double multiplier = 1.0;
      if (numberTimesProbingTotal_) {
	double useDown = numberTimesDownLocalFixed_+
	  2.0*numberTimesDownTotalFixed_/numberTimesProbingTotal_;
	double useUp = numberTimesUpLocalFixed_+
	  2.0*numberTimesUpTotalFixed_/numberTimesProbingTotal_;
#if 0
	if ((model_->specialOptions()&2048)==0&&
	    (model_->getNodeCount()%1000000)==0) {
	  int n = numberTimesProbingTotal_;
	  printf("zcol %d down (%g,%d) up (%g,%d) value %g pcost %g\n",
		 columnNumber_,
		 numberTimesDownTotalFixed_/n,numberTimesDownLocalFixed_,
		 numberTimesUpTotalFixed_/n,numberTimesUpLocalFixed_,
		 value,returnValue);
	}
#endif
      }
#else
      // probing
      // average
      double up = 1.0e-15;
      double down = 1.0e-15;
      if (numberTimesProbingTotal_) {
        up += numberTimesUpTotalFixed_ / static_cast< double >(numberTimesProbingTotal_);
        down += numberTimesDownTotalFixed_ / static_cast< double >(numberTimesProbingTotal_);
      }
      returnValue = 1 + 10.0 * std::min(numberTimesDownLocalFixed_, numberTimesUpLocalFixed_) + std::min(down, up);
      returnValue *= 1.0e-3;
#endif
    }
#endif
#ifdef COIN_DEVELOP
    History hist;
    hist.where_ = where;
    hist.status_ = ' ';
    hist.sequence_ = columnNumber_;
    hist.numberUp_ = numberTimesUp_;
    hist.numberUpInf_ = numberTimesUpInfeasible_;
    hist.sumUp_ = sumUpCost_;
    hist.upEst_ = upCost;
    hist.numberDown_ = numberTimesDown_;
    hist.numberDownInf_ = numberTimesDownInfeasible_;
    hist.sumDown_ = sumDownCost_;
    hist.downEst_ = downCost;
    if (stateOfSearch)
      addRecord(hist);
#endif
    return std::max(returnValue, 1.0e-15);
  }
}
#if CBC_DYNAMIC_EXPERIMENT == 0
// Infeasibility etc etc
usefulDynamic
CbcSimpleIntegerDynamicPseudoCost::usefulStuff(const OsiBranchingInformation *info) const
{
  usefulDynamic stuff;
  double value = model_->testSolution()[columnNumber_];
  double lower = model_->getCbcColLower()[columnNumber_];
  double upper = model_->getCbcColUpper()[columnNumber_];
  value = std::max(value, lower);
  value = std::min(value, upper);
  double nearest = floor(value + 0.5);
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double below = floor(value + integerTolerance);
  double above = below + 1.0;
  if (above > upper) {
    above = below;
    below = above - 1;
  }
  stuff.pseudoDown = above * downDynamicPseudoCost_;
  stuff.pseudoUp = below * upDynamicPseudoCost_;
  stuff.probingDown = 0.0;
  stuff.probingUp = 0.0;
  if (numberTimesProbingTotal_) {
    stuff.probingDown = numberTimesDownLocalFixed_+
      2.0*numberTimesDownTotalFixed_/numberTimesProbingTotal_;
    stuff.probingUp = numberTimesUpLocalFixed_+
      2.0*numberTimesUpTotalFixed_/numberTimesProbingTotal_;
  }
  return stuff;
}
#endif
// Creates a branching object
CbcBranchingObject *
CbcSimpleIntegerDynamicPseudoCost::createCbcBranch(OsiSolverInterface * /*solver*/,
  const OsiBranchingInformation *info, int way)
{
  double value = info->solution_[columnNumber_];
  value = std::max(value, info->lower_[columnNumber_]);
  value = std::min(value, info->upper_[columnNumber_]);
  assert(info->upper_[columnNumber_] > info->lower_[columnNumber_]);
  if (!info->hotstartSolution_ && priority_ != -999) {
#ifndef NDEBUG
#ifndef SWITCH_VARIABLES
    double nearest = floor(value + 0.5);
    assert(fabs(value - nearest) > info->integerTolerance_);
#endif
#endif
  } else if (info->hotstartSolution_) {
    double targetValue = info->hotstartSolution_[columnNumber_];
    if (way > 0)
      value = targetValue - 0.1;
    else
      value = targetValue + 0.1;
  } else {
    if (value <= info->lower_[columnNumber_])
      value += 0.1;
    else if (value >= info->upper_[columnNumber_])
      value -= 0.1;
  }
  assert(value >= info->lower_[columnNumber_] && value <= info->upper_[columnNumber_]);
  CbcDynamicPseudoCostBranchingObject *newObject = new CbcDynamicPseudoCostBranchingObject(model_, columnNumber_, way,
    value, this);
  double up = upDynamicPseudoCost_ * (ceil(value) - value);
  double down = downDynamicPseudoCost_ * (value - floor(value));
  double changeInGuessed = up - down;
  if (way > 0)
    changeInGuessed = -changeInGuessed;
  changeInGuessed = std::max(0.0, changeInGuessed);
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
  const double *solution = model_->testSolution();
  const double *lower = model_->getCbcColLower();
  const double *upper = model_->getCbcColUpper();
  double value = solution[columnNumber_];
  value = std::max(value, lower[columnNumber_]);
  value = std::min(value, upper[columnNumber_]);
  if (upper[columnNumber_] == lower[columnNumber_]) {
    // fixed
    return 0.0;
  }
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double below = floor(value + integerTolerance);
  double above = below + 1.0;
  if (above > upper[columnNumber_]) {
    above = below;
    below = above - 1;
  }
  double upCost = std::max((above - value) * upDynamicPseudoCost_, 0.0);
  return upCost;
}
// Return "down" estimate
double
CbcSimpleIntegerDynamicPseudoCost::downEstimate() const
{
  const double *solution = model_->testSolution();
  const double *lower = model_->getCbcColLower();
  const double *upper = model_->getCbcColUpper();
  double value = solution[columnNumber_];
  value = std::max(value, lower[columnNumber_]);
  value = std::min(value, upper[columnNumber_]);
  if (upper[columnNumber_] == lower[columnNumber_]) {
    // fixed
    return 0.0;
  }
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double below = floor(value + integerTolerance);
  double above = below + 1.0;
  if (above > upper[columnNumber_]) {
    above = below;
    below = above - 1;
  }
  double downCost = std::max((value - below) * downDynamicPseudoCost_, 0.0);
  return downCost;
}
// Set down pseudo cost
void CbcSimpleIntegerDynamicPseudoCost::setDownDynamicPseudoCost(double value)
{
#ifdef TRACE_ONE
  double oldDown = sumDownCost_;
#endif
  downDynamicPseudoCost_ = value;
  sumDownCost_ = std::max(sumDownCost_, value * numberTimesDown_);
#ifdef TRACE_ONE
  if (columnNumber_ == TRACE_ONE) {
    double down = downDynamicPseudoCost_ * numberTimesDown_;
    printf("For %d sumDown %g (%d), inf (%d) - pseudo %g - sumDown was %g -> %g\n",
      TRACE_ONE, down, numberTimesDown_,
      numberTimesDownInfeasible_, downDynamicPseudoCost_,
      oldDown, sumDownCost_);
  }
#endif
}
// Modify down pseudo cost in a slightly different way
void CbcSimpleIntegerDynamicPseudoCost::updateDownDynamicPseudoCost(double value)
{
  sumDownCost_ += value;
  numberTimesDown_++;
  downDynamicPseudoCost_ = sumDownCost_ / static_cast< double >(numberTimesDown_);
}
// Set up pseudo cost
void CbcSimpleIntegerDynamicPseudoCost::setUpDynamicPseudoCost(double value)
{
#ifdef TRACE_ONE
  double oldUp = sumUpCost_;
#endif
  upDynamicPseudoCost_ = value;
  sumUpCost_ = std::max(sumUpCost_, value * numberTimesUp_);
#ifdef TRACE_ONE
  if (columnNumber_ == TRACE_ONE) {
    double up = upDynamicPseudoCost_ * numberTimesUp_;
    printf("For %d sumUp %g (%d), inf (%d) - pseudo %g - sumUp was %g -> %g\n",
      TRACE_ONE, up, numberTimesUp_,
      numberTimesUpInfeasible_, upDynamicPseudoCost_,
      oldUp, sumUpCost_);
  }
#endif
}
// Modify up pseudo cost in a slightly different way
void CbcSimpleIntegerDynamicPseudoCost::updateUpDynamicPseudoCost(double value)
{
  sumUpCost_ += value;
  numberTimesUp_++;
  upDynamicPseudoCost_ = sumUpCost_ / static_cast< double >(numberTimesUp_);
}
/* Pass in information on branch just done and create CbcObjectUpdateData instance.
   If object does not need data then backward pointer will be NULL.
   Assumes can get information from solver */
CbcObjectUpdateData
CbcSimpleIntegerDynamicPseudoCost::createUpdateInformation(const OsiSolverInterface *solver,
  const CbcNode *node,
  const CbcBranchingObject *branchingObject)
{
  double originalValue = node->objectiveValue();
  int originalUnsatisfied = node->numberUnsatisfied();
  double objectiveValue = solver->getObjValue() * solver->getObjSenseInCbc();
  int unsatisfied = 0;
  int i;
  //might be base model - doesn't matter
  int numberIntegers = model_->numberIntegers();
  ;
  const double *solution = solver->getColSolution();
  //const double * lower = solver->getColLower();
  //const double * upper = solver->getColUpper();
  double change = std::max(0.0, objectiveValue - originalValue);
  int iStatus;
  if (solver->isProvenOptimal())
    iStatus = 0; // optimal
  else if (solver->isIterationLimitReached()
    && !solver->isDualObjectiveLimitReached())
    iStatus = 2; // unknown
  else
    iStatus = 1; // infeasible

  bool feasible = iStatus != 1;
  if (feasible) {
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
    const int *integerVariable = model_->integerVariable();
    for (i = 0; i < numberIntegers; i++) {
      int j = integerVariable[i];
      double value = solution[j];
      double nearest = floor(value + 0.5);
      if (fabs(value - nearest) > integerTolerance)
        unsatisfied++;
#ifdef SWITCH_VARIABLES
      const CbcSwitchingBinary *sObject = dynamic_cast< const CbcSwitchingBinary * >(this);
      if (sObject) {
        int state[3], nBadFixed;
        unsatisfied += sObject->checkAssociatedBounds(solver, solution, 0,
          state, nBadFixed);
      }
#endif
    }
  }
  int way = branchingObject->way();
  way = -way; // because after branch so moved on
  double value = branchingObject->value();
  CbcObjectUpdateData newData(this, way,
    change, iStatus,
    originalUnsatisfied - unsatisfied, value);
  newData.originalObjective_ = originalValue;
  // Solvers know about direction
  double direction = solver->getObjSenseInCbc();
  solver->getDblParam(OsiDualObjectiveLimit, newData.cutoff_);
  newData.cutoff_ *= direction;
  return newData;
}
// Just update using feasible branches and keep count of infeasible
#undef INFEAS
// Update object by CbcObjectUpdateData
void CbcSimpleIntegerDynamicPseudoCost::updateInformation(const CbcObjectUpdateData &data)
{
  bool feasible = data.status_ != 1;
  int way = data.way_;
  double value = data.branchingValue_;
  double change = data.change_;
#ifdef COIN_DEVELOP
  History hist;
  hist.where_ = 'U'; // need to tell if hot
#endif
  double movement = 0.0;
  if (way < 0) {
    // down
    movement = value - floor(value);
    if (feasible) {
#ifdef COIN_DEVELOP
      hist.status_ = 'D';
#endif
      movement = std::max(movement, MINIMUM_MOVEMENT);
      //printf("(down change %g value down %g ",change,movement);
      incrementNumberTimesDown();
      addToSumDownChange(1.0e-30 + movement);
      addToSumDownDecrease(data.intDecrease_);
#if TYPE2 == 0
      addToSumDownCost(change / (1.0e-30 + movement));
      setDownDynamicPseudoCost(sumDownCost() / static_cast< double >(numberTimesDown()));
#elif TYPE2 == 1
      addToSumDownCost(change);
      setDownDynamicPseudoCost(sumDownCost() / sumDownChange());
#elif TYPE2 == 2
      addToSumDownCost(change * TYPERATIO + (1.0 - TYPERATIO) * change / (1.0e-30 + movement));
      setDownDynamicPseudoCost(sumDownCost() * (TYPERATIO / sumDownChange() + (1.0 - TYPERATIO) / (double)numberTimesDown()));
#endif
#if CBC_DYNAMIC_EXPERIMENT > 0
      if (change>movement*downDynamicPseudoCost_) {
	double over = change/(movement*downDynamicPseudoCost_);
	downOver_ = (numberTimesDown_-1)*downOver_ + std::min(over,10.0);
	downOver_ /= numberTimesDown_;
      } else {
	double under = change/(movement*downDynamicPseudoCost_);
	downUnder_ = (numberTimesDown_-1)*downUnder_ + std::max(under,0.1);
	downUnder_ /= numberTimesDown_;
      }
#endif
    } else {
#ifdef COIN_DEVELOP
      hist.status_ = 'd';
#endif
      //printf("(down infeasible value down %g ",change,movement);
      incrementNumberTimesDown();
      incrementNumberTimesDownInfeasible();
#if INFEAS == 2
      double distanceToCutoff = 0.0;
      double objectiveValue = model->getCurrentMinimizationObjValue();
      distanceToCutoff = model->getCutoff() - originalValue;
      if (distanceToCutoff < 1.0e20)
        change = distanceToCutoff * 2.0;
      else
        change = downDynamicPseudoCost() * movement * 10.0;
      change = std::max(1.0e-12 * (1.0 + fabs(originalValue)), change);
      addToSumDownChange(1.0e-30 + movement);
      addToSumDownDecrease(data.intDecrease_);
#if TYPE2 == 0
      addToSumDownCost(change / (1.0e-30 + movement));
      setDownDynamicPseudoCost(sumDownCost() / (double)numberTimesDown());
#elif TYPE2 == 1
      addToSumDownCost(change);
      setDownDynamicPseudoCost(sumDownCost() / sumDownChange());
#elif TYPE2 == 2
      addToSumDownCost(change * TYPERATIO + (1.0 - TYPERATIO) * change / (1.0e-30 + movement));
      setDownDynamicPseudoCost(sumDownCost() * (TYPERATIO / sumDownChange() + (1.0 - TYPERATIO) / (double)numberTimesDown()));
#endif
#endif
    }
#if INFEAS == 1
    double sum = sumDownCost_;
    int number = numberTimesDown_;
    double originalValue = data.originalObjective_;
    assert(originalValue != COIN_DBL_MAX);
    double distanceToCutoff = data.cutoff_ - originalValue;
    if (distanceToCutoff > 1.0e20)
      distanceToCutoff = 10.0 + fabs(originalValue);
    sum += INFEAS_MULTIPLIER * numberTimesDownInfeasible_ * std::max(distanceToCutoff, 1.0e-12 * (1.0 + fabs(originalValue)));
    setDownDynamicPseudoCost(sum / static_cast< double >(number));
#endif
  } else {
    // up
    movement = ceil(value) - value;
    if (feasible) {
#ifdef COIN_DEVELOP
      hist.status_ = 'U';
#endif
      movement = std::max(movement, MINIMUM_MOVEMENT);
      //printf("(up change %g value down %g ",change,movement);
      incrementNumberTimesUp();
      addToSumUpChange(1.0e-30 + movement);
      addToSumUpDecrease(data.intDecrease_);
#if TYPE2 == 0
      addToSumUpCost(change / (1.0e-30 + movement));
      setUpDynamicPseudoCost(sumUpCost() / static_cast< double >(numberTimesUp()));
#elif TYPE2 == 1
      addToSumUpCost(change);
      setUpDynamicPseudoCost(sumUpCost() / sumUpChange());
#elif TYPE2 == 2
      addToSumUpCost(change * TYPERATIO + (1.0 - TYPERATIO) * change / (1.0e-30 + movement));
      setUpDynamicPseudoCost(sumUpCost() * (TYPERATIO / sumUpChange() + (1.0 - TYPERATIO) / (double)numberTimesUp()));
#endif
#if CBC_DYNAMIC_EXPERIMENT > 0
      if (change>movement*upDynamicPseudoCost_) {
	double over = change/(movement*upDynamicPseudoCost_);
	upOver_ = (numberTimesUp_-1)*upOver_ + std::min(over,10.0);
	upOver_ /= numberTimesUp_;
      } else {
	double under = change/(movement*upDynamicPseudoCost_);
	upUnder_ = (numberTimesUp_-1)*upUnder_ + std::max(under,0.1);
	upUnder_ /= numberTimesUp_;
      }
#endif
    } else {
#ifdef COIN_DEVELOP
      hist.status_ = 'u';
#endif
      //printf("(up infeasible value down %g ",change,movement);
      incrementNumberTimesUp();
      incrementNumberTimesUpInfeasible();
#if INFEAS == 2
      double distanceToCutoff = 0.0;
      double objectiveValue = model->getCurrentMinimizationObjValue();
      distanceToCutoff = model->getCutoff() - originalValue;
      if (distanceToCutoff < 1.0e20)
        change = distanceToCutoff * 2.0;
      else
        change = upDynamicPseudoCost() * movement * 10.0;
      change = std::max(1.0e-12 * (1.0 + fabs(originalValue)), change);
      addToSumUpChange(1.0e-30 + movement);
      addToSumUpDecrease(data.intDecrease_);
#if TYPE2 == 0
      addToSumUpCost(change / (1.0e-30 + movement));
      setUpDynamicPseudoCost(sumUpCost() / (double)numberTimesUp());
#elif TYPE2 == 1
      addToSumUpCost(change);
      setUpDynamicPseudoCost(sumUpCost() / sumUpChange());
#elif TYPE2 == 2
      addToSumUpCost(change * TYPERATIO + (1.0 - TYPERATIO) * change / (1.0e-30 + movement));
      setUpDynamicPseudoCost(sumUpCost() * (TYPERATIO / sumUpChange() + (1.0 - TYPERATIO) / (double)numberTimesUp()));
#endif
#endif
    }
#if INFEAS == 1
    double sum = sumUpCost_;
    int number = numberTimesUp_;
    double originalValue = data.originalObjective_;
    assert(originalValue != COIN_DBL_MAX);
    double distanceToCutoff = data.cutoff_ - originalValue;
    if (distanceToCutoff > 1.0e20)
      distanceToCutoff = 10.0 + fabs(originalValue);
    sum += INFEAS_MULTIPLIER * numberTimesUpInfeasible_ * std::max(distanceToCutoff, 1.0e-12 * (1.0 + fabs(originalValue)));
    setUpDynamicPseudoCost(sum / static_cast< double >(number));
#endif
  }
  if (data.way_ < 0)
    assert(numberTimesDown_ > 0);
  else
    assert(numberTimesUp_ > 0);
  assert(downDynamicPseudoCost_ >= 0.0 && downDynamicPseudoCost_ < 1.0e100);
  downDynamicPseudoCost_ = std::max(1.0e-10, downDynamicPseudoCost_);
  assert(upDynamicPseudoCost_ >= 0.0 && upDynamicPseudoCost_ < 1.0e100);
  upDynamicPseudoCost_ = std::max(1.0e-10, upDynamicPseudoCost_);
#ifdef COIN_DEVELOP
  hist.sequence_ = columnNumber_;
  hist.numberUp_ = numberTimesUp_;
  hist.numberUpInf_ = numberTimesUpInfeasible_;
  hist.sumUp_ = sumUpCost_;
  hist.upEst_ = change;
  hist.numberDown_ = numberTimesDown_;
  hist.numberDownInf_ = numberTimesDownInfeasible_;
  hist.sumDown_ = sumDownCost_;
  hist.downEst_ = movement;
  addRecord(hist);
#endif
  //print(1,0.5);
  assert(downDynamicPseudoCost_ > 1.0e-40 && upDynamicPseudoCost_ > 1.0e-40);
#if MOD_SHADOW > 1
  if (upShadowPrice_ > 0.0 && numberTimesDown_ >= numberBeforeTrust_
    && numberTimesUp_ >= numberBeforeTrust_) {
    // Set negative
    upShadowPrice_ = -upShadowPrice_;
    assert(downShadowPrice_ > 0.0);
    downShadowPrice_ = -downShadowPrice_;
  }
#endif
}
// Updates stuff like pseudocosts after mini branch and bound
void CbcSimpleIntegerDynamicPseudoCost::updateAfterMini(int numberDown, int numberDownInfeasible,
  double sumDown, int numberUp,
  int numberUpInfeasible, double sumUp)
{
  numberTimesDown_ = numberDown;
  numberTimesDownInfeasible_ = numberDownInfeasible;
  sumDownCost_ = sumDown;
  numberTimesUp_ = numberUp;
  numberTimesUpInfeasible_ = numberUpInfeasible;
  sumUpCost_ = sumUp;
  if (numberTimesDown_ > 0) {
    setDownDynamicPseudoCost(sumDownCost_ / static_cast< double >(numberTimesDown_));
    assert(downDynamicPseudoCost_ > 0.0 && downDynamicPseudoCost_ < 1.0e50);
  }
  if (numberTimesUp_ > 0) {
    setUpDynamicPseudoCost(sumUpCost_ / static_cast< double >(numberTimesUp_));
    assert(upDynamicPseudoCost_ > 0.0 && upDynamicPseudoCost_ < 1.0e50);
  }
  assert(downDynamicPseudoCost_ > 1.0e-40 && upDynamicPseudoCost_ > 1.0e-40);
}
#if CBC_DYNAMIC_EXPERIMENT == 0
// Pass in probing information
void CbcSimpleIntegerDynamicPseudoCost::setProbingInformation(int fixedDown, int fixedUp)
{
  numberTimesProbingTotal_++;
  numberTimesDownLocalFixed_ = fixedDown;
  numberTimesDownTotalFixed_ += fixedDown;
  numberTimesUpLocalFixed_ = fixedUp;
  numberTimesUpTotalFixed_ += fixedUp;
}
#endif
// Print
void CbcSimpleIntegerDynamicPseudoCost::print(int type, double value) const
{
  if (!type) {
    double meanDown = 0.0;
    double devDown = 0.0;
    if (numberTimesDown_) {
      meanDown = sumDownCost_ / static_cast< double >(numberTimesDown_);
      devDown = meanDown * meanDown - 2.0 * meanDown * sumDownCost_;
      if (devDown >= 0.0)
        devDown = sqrt(devDown);
    }
    double meanUp = 0.0;
    double devUp = 0.0;
    if (numberTimesUp_) {
      meanUp = sumUpCost_ / static_cast< double >(numberTimesUp_);
      devUp = meanUp * meanUp - 2.0 * meanUp * sumUpCost_;
      if (devUp >= 0.0)
        devUp = sqrt(devUp);
    }
    printf("%d down %d times (%d inf) mean %g (dev %g) up %d times (%d inf) mean %g (dev %g)\n",
      columnNumber_,
      numberTimesDown_, numberTimesDownInfeasible_, meanDown, devDown,
      numberTimesUp_, numberTimesUpInfeasible_, meanUp, devUp);
  } else {
    const double *upper = model_->getCbcColUpper();
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double below = floor(value + integerTolerance);
    double above = below + 1.0;
    if (above > upper[columnNumber_]) {
      above = below;
      below = above - 1;
    }
    double objectiveValue = model_->getCurrentMinimizationObjValue();
    double distanceToCutoff = model_->getCutoff() - objectiveValue;
    if (distanceToCutoff < 1.0e20)
      distanceToCutoff *= 10.0;
    else
      distanceToCutoff = 1.0e2 + fabs(objectiveValue);
    distanceToCutoff = std::max(distanceToCutoff, 1.0e-12 * (1.0 + fabs(objectiveValue)));
    double sum;
    int number;
    double downCost = std::max(value - below, 0.0);
    double downCost0 = downCost * downDynamicPseudoCost_;
    sum = sumDownCost();
    number = numberTimesDown();
    sum += INFEAS_MULTIPLIER * numberTimesDownInfeasible() * (distanceToCutoff / (downCost + 1.0e-12));
    if (number > 0)
      downCost *= sum / static_cast< double >(number);
    else
      downCost *= downDynamicPseudoCost_;
    double upCost = std::max((above - value), 0.0);
    double upCost0 = upCost * upDynamicPseudoCost_;
    sum = sumUpCost();
    number = numberTimesUp();
    sum += INFEAS_MULTIPLIER * numberTimesUpInfeasible() * (distanceToCutoff / (upCost + 1.0e-12));
    if (number > 0)
      upCost *= sum / static_cast< double >(number);
    else
      upCost *= upDynamicPseudoCost_;
    printf("%d down %d times %g (est %g)  up %d times %g (est %g)\n",
      columnNumber_,
      numberTimesDown_, downCost, downCost0,
      numberTimesUp_, upCost, upCost0);
  }
}

//##############################################################################

// Default Constructor
CbcIntegerPseudoCostBranchingObject::CbcIntegerPseudoCostBranchingObject()
  : CbcIntegerBranchingObject()
{
  changeInGuessed_ = 1.0e-5;
}

// Useful constructor
CbcIntegerPseudoCostBranchingObject::CbcIntegerPseudoCostBranchingObject(CbcModel *model,
  int variable, int way, double value)
  : CbcIntegerBranchingObject(model, variable, way, value)
{
}
// Useful constructor for fixing
CbcIntegerPseudoCostBranchingObject::CbcIntegerPseudoCostBranchingObject(CbcModel *model,
  int variable, int way,
  double lowerValue,
  double /*upperValue*/)
  : CbcIntegerBranchingObject(model, variable, way, lowerValue)
{
  changeInGuessed_ = 1.0e100;
}

// Copy constructor
CbcIntegerPseudoCostBranchingObject::CbcIntegerPseudoCostBranchingObject(
  const CbcIntegerPseudoCostBranchingObject &rhs)
  : CbcIntegerBranchingObject(rhs)
{
  changeInGuessed_ = rhs.changeInGuessed_;
}

// Assignment operator
CbcIntegerPseudoCostBranchingObject &
CbcIntegerPseudoCostBranchingObject::operator=(const CbcIntegerPseudoCostBranchingObject &rhs)
{
  if (this != &rhs) {
    CbcIntegerBranchingObject::operator=(rhs);
    changeInGuessed_ = rhs.changeInGuessed_;
  }
  return *this;
}
CbcBranchingObject *
CbcIntegerPseudoCostBranchingObject::clone() const
{
  return (new CbcIntegerPseudoCostBranchingObject(*this));
}

// Destructor
CbcIntegerPseudoCostBranchingObject::~CbcIntegerPseudoCostBranchingObject()
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
CbcIntegerPseudoCostBranchingObject::branch()
{
  CbcIntegerBranchingObject::branch();
  return changeInGuessed_;
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcIntegerPseudoCostBranchingObject::compareBranchingObject(const CbcBranchingObject *brObj, const bool replaceIfOverlap)
{
  const CbcIntegerPseudoCostBranchingObject *br = dynamic_cast< const CbcIntegerPseudoCostBranchingObject * >(brObj);
  assert(br);
  double *thisBd = way_ < 0 ? down_ : up_;
  const double *otherBd = br->way_ < 0 ? br->down_ : br->up_;
  return CbcCompareRanges(thisBd, otherBd, replaceIfOverlap);
}
#ifdef SWITCH_VARIABLES
/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
CbcSwitchingBinary::CbcSwitchingBinary()
  : CbcSimpleIntegerDynamicPseudoCost()
  , zeroLowerBound_(NULL)
  , oneLowerBound_(NULL)
  , zeroUpperBound_(NULL)
  , oneUpperBound_(NULL)
  , otherVariable_(NULL)
  , numberOther_(0)
  , type_(0)
{
}

/** Useful constructor
*/
CbcSwitchingBinary::CbcSwitchingBinary(CbcSimpleIntegerDynamicPseudoCost *oldObject,
  int nOdd, const int *other, const int *otherRow)
  : CbcSimpleIntegerDynamicPseudoCost(*oldObject)
  , zeroLowerBound_(NULL)
  , oneLowerBound_(NULL)
  , zeroUpperBound_(NULL)
  , oneUpperBound_(NULL)
  , otherVariable_(NULL)
  , numberOther_(0)
  , type_(0)
{
  if (nOdd)
    type_ = 2;
  const CoinPackedMatrix *rowCopy = model_->solver()->getMatrixByRow();
  const int *column = rowCopy->getIndices();
  //const int * rowLength = rowCopy->getVectorLengths();
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
  //const double * rowLower = model_->solver()->getRowLower();
  const double *rowUpper = model_->solver()->getRowUpper();
  const double *columnLower = model_->solver()->getColLower();
  const double *columnUpper = model_->solver()->getColUpper();
  const double *element = rowCopy->getElements();
  int last = other[0];
  int nPair = 0;
  int nInGroup = 1;
  for (int i = 1; i <= nOdd; i++) {
    if (other[i] == last) {
      nInGroup++;
    } else {
      if (nInGroup > 2 && model_->logLevel() > 2)
        printf("%d in group for column %d - some redundancy\n",
          nInGroup, columnNumber_);
      nPair++;
      last = other[i];
      nInGroup = 1;
    }
  }
  zeroLowerBound_ = new double[4 * nPair];
  oneLowerBound_ = zeroLowerBound_ + nPair;
  zeroUpperBound_ = oneLowerBound_ + nPair;
  oneUpperBound_ = zeroUpperBound_ + nPair;
  otherVariable_ = new int[nPair];
  numberOther_ = nPair;
  if (nPair > 1 && model_->logLevel() > 2)
    printf("%d pairs for column %d\n",
      nPair, columnNumber_);
  // Now fill
  last = other[0];
  nPair = 0;
  int rows[20];
  rows[0] = otherRow[0];
  nInGroup = 1;
  for (int i = 1; i <= nOdd; i++) {
    if (other[i] == last) {
      rows[nInGroup++] = otherRow[i];
    } else {
      double newLowerZero = 0.0;
      double newUpperZero = COIN_DBL_MAX;
      double newLowerOne = 0.0;
      double newUpperOne = COIN_DBL_MAX;
      int cColumn = -1;
      for (int j = 0; j < nInGroup; j++) {
        int iRow = rows[j];
        CoinBigIndex k = rowStart[iRow];
        double bValue, cValue;
        if (column[k] == columnNumber_) {
          bValue = element[k];
          cValue = element[k + 1];
          cColumn = column[k + 1];
        } else {
          bValue = element[k + 1];
          cValue = element[k];
          cColumn = column[k];
        }
        if (rowUpper[iRow]) {
          // G row - convert to L
          bValue = -bValue;
          cValue = -cValue;
        }
        if (bValue > 0.0) {
          // binary*abs(bValue) <= continuous*abs(cValue);
          newLowerOne = -bValue / cValue;
        } else {
          // binary*abs(bValue) >= continuous*abs(cValue);
          newUpperOne = -bValue / cValue;
          newUpperZero = 0.0;
        }
      }
      zeroLowerBound_[nPair] = newLowerZero;
      oneLowerBound_[nPair] = newLowerOne;
      zeroUpperBound_[nPair] = newUpperZero;
      oneUpperBound_[nPair] = newUpperOne;
      // make current bounds tight
      double newLower = std::min(newLowerZero, newLowerOne);
      if (newLower > columnLower[cColumn])
        model_->solver()->setColLower(cColumn, newLower);
      double newUpper = std::max(newUpperZero, newUpperOne);
      if (newUpper < columnUpper[cColumn])
        model_->solver()->setColUpper(cColumn, newUpper);
      otherVariable_[nPair++] = cColumn;
      last = other[i];
      rows[0] = otherRow[i];
      nInGroup = 1;
    }
  }
}
// Copy constructor
CbcSwitchingBinary::CbcSwitchingBinary(const CbcSwitchingBinary &rhs)
  : CbcSimpleIntegerDynamicPseudoCost(rhs)
  , numberOther_(rhs.numberOther_)
  , type_(rhs.type_)
{
  zeroLowerBound_ = CoinCopyOfArray(rhs.zeroLowerBound_, 4 * numberOther_);
  oneLowerBound_ = zeroLowerBound_ + numberOther_;
  zeroUpperBound_ = oneLowerBound_ + numberOther_;
  oneUpperBound_ = zeroUpperBound_ + numberOther_;
  otherVariable_ = CoinCopyOfArray(rhs.otherVariable_, numberOther_);
}

// Clone
CbcObject *
CbcSwitchingBinary::clone() const
{
  return new CbcSwitchingBinary(*this);
}

// Assignment operator
CbcSwitchingBinary &
CbcSwitchingBinary::operator=(const CbcSwitchingBinary &rhs)
{
  if (this != &rhs) {
    CbcSimpleIntegerDynamicPseudoCost::operator=(rhs);
    numberOther_ = rhs.numberOther_;
    type_ = rhs.type_;
    delete[] zeroLowerBound_;
    delete[] otherVariable_;
    zeroLowerBound_ = CoinCopyOfArray(rhs.zeroLowerBound_, 4 * numberOther_);
    oneLowerBound_ = zeroLowerBound_ + numberOther_;
    zeroUpperBound_ = oneLowerBound_ + numberOther_;
    oneUpperBound_ = zeroUpperBound_ + numberOther_;
    otherVariable_ = CoinCopyOfArray(rhs.otherVariable_, numberOther_);
  }
  return *this;
}

// Destructor
CbcSwitchingBinary::~CbcSwitchingBinary()
{
  delete[] zeroLowerBound_;
  delete[] otherVariable_;
}
// Add in zero switches
void CbcSwitchingBinary::addZeroSwitches(int nAdd, const int *columns)
{
  type_ |= 1;
  int nNew = numberOther_ + nAdd;
  double *bounds = new double[4 * nNew];
  int *other = new int[nNew];
  memcpy(other, otherVariable_, numberOther_ * sizeof(int));
  delete[] otherVariable_;
  otherVariable_ = other;
  memcpy(bounds, zeroLowerBound_, numberOther_ * sizeof(double));
  memcpy(bounds + nNew, oneLowerBound_, numberOther_ * sizeof(double));
  memcpy(bounds + 2 * nNew, zeroUpperBound_, numberOther_ * sizeof(double));
  memcpy(bounds + 3 * nNew, oneUpperBound_, numberOther_ * sizeof(double));
  delete[] zeroLowerBound_;
  zeroLowerBound_ = bounds;
  oneLowerBound_ = zeroLowerBound_ + nNew;
  zeroUpperBound_ = oneLowerBound_ + nNew;
  oneUpperBound_ = zeroUpperBound_ + nNew;
  for (int i = 0; i < nAdd; i++) {
    zeroLowerBound_[numberOther_] = 0.0;
    oneLowerBound_[numberOther_] = 0.0;
    zeroUpperBound_[numberOther_] = 0.0;
    oneUpperBound_[numberOther_] = COIN_DBL_MAX;
    otherVariable_[numberOther_++] = columns[i];
  }
}
// Same - returns true if contents match(ish)
bool CbcSwitchingBinary::same(const CbcSwitchingBinary *otherObject) const
{
  bool okay = CbcSimpleIntegerDynamicPseudoCost::same(otherObject);
  return okay;
}
double
CbcSwitchingBinary::infeasibility(const OsiBranchingInformation *info,
  int &preferredWay) const
{
  assert(downDynamicPseudoCost_ > 1.0e-40 && upDynamicPseudoCost_ > 1.0e-40);
  double *solution = const_cast< double * >(model_->testSolution());
  const double *lower = model_->getCbcColLower();
  const double *upper = model_->getCbcColUpper();
  double saveValue = solution[columnNumber_];
  if (!lower[columnNumber_] && upper[columnNumber_] == 1.0) {
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
    if (saveValue < integerTolerance) {
      // check others OK
      bool allGood = true;
      double tolerance;
      model_->solver()->getDblParam(OsiPrimalTolerance, tolerance);
      for (int i = 0; i < numberOther_; i++) {
        int otherColumn = otherVariable_[i];
        double value = solution[otherColumn];
        if (value < zeroLowerBound_[i] - tolerance || value > zeroUpperBound_[i] + tolerance)
          allGood = false;
      }
      if (!allGood)
        solution[columnNumber_] = 2.0 * integerTolerance;
    } else if (saveValue > 1.0 - integerTolerance) {
      // check others OK
      bool allGood = true;
      double tolerance;
      model_->solver()->getDblParam(OsiPrimalTolerance, tolerance);
      for (int i = 0; i < numberOther_; i++) {
        int otherColumn = otherVariable_[i];
        double value = solution[otherColumn];
        if (value < oneLowerBound_[i] - tolerance || value > oneUpperBound_[i] + tolerance)
          allGood = false;
      }
      if (!allGood)
        solution[columnNumber_] = 1.0 - 2.0 * integerTolerance;
    }
  }
  double inf = CbcSimpleIntegerDynamicPseudoCost::infeasibility(info, preferredWay);
  solution[columnNumber_] = saveValue;
  return inf;
}
// Set associated bounds
int CbcSwitchingBinary::setAssociatedBounds(OsiSolverInterface *solver,
  int cleanBasis) const
{
  if (!solver)
    solver = model_->solver();
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver);
  if (cleanBasis != 1)
    clpSolver = NULL;
  const double *columnLower = solver->getColLower();
  const double *columnUpper = solver->getColUpper();
  int nChanged = 0;
  if (!columnUpper[columnNumber_]) {
    if (clpSolver)
      clpSolver->setColumnStatus(columnNumber_, ClpSimplex::isFixed);
    for (int i = 0; i < numberOther_; i++) {
      int otherColumn = otherVariable_[i];
      if (zeroLowerBound_[i] > columnLower[otherColumn]) {
        solver->setColLower(otherColumn, zeroLowerBound_[i]);
        nChanged++;
      }
      if (zeroUpperBound_[i] < columnUpper[otherColumn]) {
        solver->setColUpper(otherColumn, zeroUpperBound_[i]);
#ifdef COIN_DEVELOP
        const double *solution = solver->getColSolution();
        double value = solution[otherColumn];
        if (value - zeroUpperBound_[i] > 1.0e-5 && model_->logLevel() > 1)
          printf("value for continuous %d %g - above %g - switch %d is %.12g (ub 0)\n",
            otherColumn, value, zeroUpperBound_[i], columnNumber_, solution[columnNumber_]);
#endif
        nChanged++;
      }
    }
  } else if (columnLower[columnNumber_] == 1.0) {
    if (clpSolver)
      clpSolver->setColumnStatus(columnNumber_, ClpSimplex::isFixed);
    for (int i = 0; i < numberOther_; i++) {
      int otherColumn = otherVariable_[i];
      if (oneLowerBound_[i] > columnLower[otherColumn]) {
        solver->setColLower(otherColumn, oneLowerBound_[i]);
        nChanged++;
      }
      if (oneUpperBound_[i] < columnUpper[otherColumn]) {
        solver->setColUpper(otherColumn, oneUpperBound_[i]);
        nChanged++;
      }
    }
  } else if (cleanBasis >= 2) {
    // if all OK then can fix
    int state[3];
    int nBadFixed;
    const double *solution = solver->getColSolution();
    if (!checkAssociatedBounds(solver, solution,
          0, state, nBadFixed)) {
      const double *reducedCost = solver->getReducedCost();
      double good = true;
      double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
      if (solution[columnNumber_] < integerTolerance) {
        if (cleanBasis == 2 || reducedCost[columnNumber_] > 1.0e-6)
          solver->setColUpper(columnNumber_, 0.0);
        else
          good = false;
      } else if (solution[columnNumber_] > 1.0 - integerTolerance) {
        if (cleanBasis == 2 || reducedCost[columnNumber_] < -1.0e-6)
          solver->setColLower(columnNumber_, 1.0);
        else
          good = false;
      }
      if (good)
        nChanged = setAssociatedBounds(solver, 0);
    }
  } else {
    // see if any continuous bounds force binary
    for (int i = 0; i < numberOther_; i++) {
      int otherColumn = otherVariable_[i];
      if (columnLower[otherColumn] > zeroUpperBound_[i]) {
        // can't be zero
        solver->setColLower(columnNumber_, 1.0);
        nChanged++;
      } else if (columnLower[otherColumn] > oneUpperBound_[i]) {
        // can't be one
        solver->setColUpper(columnNumber_, 0.0);
        nChanged++;
      }
      if (columnUpper[otherColumn] < zeroLowerBound_[i]) {
        // can't be zero
        solver->setColLower(columnNumber_, 1.0);
        nChanged++;
      } else if (columnUpper[otherColumn] < oneLowerBound_[i]) {
        // can't be one
        solver->setColUpper(columnNumber_, 0.0);
        nChanged++;
      }
    }
  }
  return nChanged;
}
// Check associated bounds
int CbcSwitchingBinary::checkAssociatedBounds(const OsiSolverInterface *solver,
  const double *solution, int printLevel,
  int state[3], int &nBadFixed) const
{
  state[0] = 0;
  int nBad = 0;
  if (!solver)
    solver = model_->solver();
  double tolerance;
  solver->getDblParam(OsiPrimalTolerance, tolerance);
  const double *columnLower = solver->getColLower();
  const double *columnUpper = solver->getColUpper();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  bool printIt = printLevel > 2 && model_->logLevel() > 1;
  if (solution[columnNumber_] < integerTolerance) {
    state[0] = -1;
    for (int i = 0; i < numberOther_; i++) {
      int otherColumn = otherVariable_[i];
      if (zeroLowerBound_[i] > solution[otherColumn] + tolerance * 5.0) {
        nBad++;
        if (columnUpper[columnNumber_] == 0.0) {
          nBadFixed++;
          //printIt=true;
        }
        if (printIt)
          printf("switch %d at zero, other %d at %.12g below bound of %.12g\n",
            columnNumber_, otherColumn, solution[otherColumn], zeroLowerBound_[i]);
      }
      if (zeroUpperBound_[i] < solution[otherColumn] - tolerance * 5.0) {
        nBad++;
        if (columnUpper[columnNumber_] == 0.0) {
          nBadFixed++;
          //printIt=true;
        }
        if (printIt)
          printf("switch %d at zero, other %d at %.12g above bound of %.12g\n",
            columnNumber_, otherColumn, solution[otherColumn], zeroUpperBound_[i]);
      }
    }
  } else if (solution[columnNumber_] > 1.0 - integerTolerance) {
    state[0] = 1;
    for (int i = 0; i < numberOther_; i++) {
      int otherColumn = otherVariable_[i];
      if (oneLowerBound_[i] > solution[otherColumn] + tolerance * 5.0) {
        nBad++;
        if (columnLower[columnNumber_] == 1.0) {
          nBadFixed++;
          //printIt=true;
        }
        if (printIt)
          printf("switch %d at one, other %d at %.12g below bound of %.12g\n",
            columnNumber_, otherColumn, solution[otherColumn], oneLowerBound_[i]);
      }
      if (oneUpperBound_[i] < solution[otherColumn] - tolerance * 5.0) {
        nBad++;
        if (columnLower[columnNumber_] == 1.0) {
          nBadFixed++;
          //printIt=true;
        }
        if (printIt)
          printf("switch %d at one, other %d at %.12g above bound of %.12g\n",
            columnNumber_, otherColumn, solution[otherColumn], oneUpperBound_[i]);
      }
    }
  } else {
    // in between - compute tight variables
    state[1] = 0;
    state[2] = 0;
    // for now just compute ones away from bounds
    for (int i = 0; i < numberOther_; i++) {
      int otherColumn = otherVariable_[i];
      double otherValue = solution[otherColumn];
      if (otherValue > columnLower[otherColumn] + tolerance && otherValue < columnUpper[otherColumn] - tolerance)
        state[1]++;
    }
  }
  return nBad;
}
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
