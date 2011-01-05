// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CoinTypes.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcSimpleIntegerPseudoCost.hpp"
#include "CbcSimpleIntegerDynamicPseudoCost.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

//##############################################################################

/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
CbcSimpleIntegerPseudoCost::CbcSimpleIntegerPseudoCost ()
        : CbcSimpleInteger(),
        downPseudoCost_(1.0e-5),
        upPseudoCost_(1.0e-5),
        upDownSeparator_(-1.0),
        method_(0)
{
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
CbcSimpleIntegerPseudoCost::CbcSimpleIntegerPseudoCost (CbcModel * model,
        int iColumn, double breakEven)
        : CbcSimpleInteger(model, iColumn, breakEven)
{
    const double * cost = model->getObjCoefficients();
    double costValue = CoinMax(1.0e-5, fabs(cost[iColumn]));
    // treat as if will cost what it says up
    upPseudoCost_ = costValue;
    // and balance at breakeven
    downPseudoCost_ = ((1.0 - breakEven_) * upPseudoCost_) / breakEven_;
    upDownSeparator_ = -1.0;
    method_ = 0;
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
CbcSimpleIntegerPseudoCost::CbcSimpleIntegerPseudoCost (CbcModel * model,
        int iColumn, double downPseudoCost,
        double upPseudoCost)
        : CbcSimpleInteger(model, iColumn)
{
    downPseudoCost_ = CoinMax(1.0e-10, downPseudoCost);
    upPseudoCost_ = CoinMax(1.0e-10, upPseudoCost);
    breakEven_ = upPseudoCost_ / (upPseudoCost_ + downPseudoCost_);
    upDownSeparator_ = -1.0;
    method_ = 0;
}
// Useful constructor - passed and model index and pseudo costs
CbcSimpleIntegerPseudoCost::CbcSimpleIntegerPseudoCost (CbcModel * model,
        int /*dummy*/,
        int iColumn,
        double downPseudoCost, double upPseudoCost)
{
    *this = CbcSimpleIntegerPseudoCost(model, iColumn, downPseudoCost, upPseudoCost);
    columnNumber_ = iColumn;
}

// Copy constructor
CbcSimpleIntegerPseudoCost::CbcSimpleIntegerPseudoCost ( const CbcSimpleIntegerPseudoCost & rhs)
        : CbcSimpleInteger(rhs),
        downPseudoCost_(rhs.downPseudoCost_),
        upPseudoCost_(rhs.upPseudoCost_),
        upDownSeparator_(rhs.upDownSeparator_),
        method_(rhs.method_)

{
}

// Clone
CbcObject *
CbcSimpleIntegerPseudoCost::clone() const
{
    return new CbcSimpleIntegerPseudoCost(*this);
}

// Assignment operator
CbcSimpleIntegerPseudoCost &
CbcSimpleIntegerPseudoCost::operator=( const CbcSimpleIntegerPseudoCost & rhs)
{
    if (this != &rhs) {
        CbcSimpleInteger::operator=(rhs);
        downPseudoCost_ = rhs.downPseudoCost_;
        upPseudoCost_ = rhs.upPseudoCost_;
        upDownSeparator_ = rhs.upDownSeparator_;
        method_ = rhs.method_;
    }
    return *this;
}

// Destructor
CbcSimpleIntegerPseudoCost::~CbcSimpleIntegerPseudoCost ()
{
}
CbcBranchingObject *
CbcSimpleIntegerPseudoCost::createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * /*info*/, int way)
{
    //OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    double value = solution[columnNumber_];
    value = CoinMax(value, lower[columnNumber_]);
    value = CoinMin(value, upper[columnNumber_]);
#ifndef NDEBUG
    double nearest = floor(value + 0.5);
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
    assert (upper[columnNumber_] > lower[columnNumber_]);
#endif
    if (!model_->hotstartSolution()) {
        assert (fabs(value - nearest) > integerTolerance);
    } else {
        const double * hotstartSolution = model_->hotstartSolution();
        double targetValue = hotstartSolution[columnNumber_];
        if (way > 0)
            value = targetValue - 0.1;
        else
            value = targetValue + 0.1;
    }
    CbcIntegerPseudoCostBranchingObject * newObject =
        new CbcIntegerPseudoCostBranchingObject(model_, columnNumber_, way,
                                                value);
    double up =  upPseudoCost_ * (ceil(value) - value);
    double down =  downPseudoCost_ * (value - floor(value));
    double changeInGuessed = up - down;
    if (way > 0)
        changeInGuessed = - changeInGuessed;
    changeInGuessed = CoinMax(0.0, changeInGuessed);
    //if (way>0)
    //changeInGuessed += 1.0e8; // bias to stay up
    newObject->setChangeInGuessed(changeInGuessed);
    newObject->setOriginalObject(this);
    return newObject;
}
double
CbcSimpleIntegerPseudoCost::infeasibility(const OsiBranchingInformation * /*info*/,
        int &preferredWay) const
{
    OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    if (upper[columnNumber_] == lower[columnNumber_]) {
        // fixed
        preferredWay = 1;
        return 0.0;
    }
    double value = solution[columnNumber_];
    value = CoinMax(value, lower[columnNumber_]);
    value = CoinMin(value, upper[columnNumber_]);
    /*printf("%d %g %g %g %g\n",columnNumber_,value,lower[columnNumber_],
      solution[columnNumber_],upper[columnNumber_]);*/
    double nearest = floor(value + 0.5);
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double below = floor(value + integerTolerance);
    double above = below + 1.0;
    if (above > upper[columnNumber_]) {
        above = below;
        below = above - 1;
    }
    double downCost = CoinMax((value - below) * downPseudoCost_, 0.0);
    double upCost = CoinMax((above - value) * upPseudoCost_, 0.0);
    if (downCost >= upCost)
        preferredWay = 1;
    else
        preferredWay = -1;
    // See if up down choice set
    if (upDownSeparator_ > 0.0) {
        preferredWay = (value - below >= upDownSeparator_) ? 1 : -1;
    }
    if (preferredWay_)
        preferredWay = preferredWay_;
    if (fabs(value - nearest) <= integerTolerance) {
        return 0.0;
    } else {
        // can't get at model so 1,2 don't make sense
        assert(method_ < 1 || method_ > 2);
        if (!method_)
            return CoinMin(downCost, upCost);
        else
            return CoinMax(downCost, upCost);
    }
}

// Return "up" estimate
double
CbcSimpleIntegerPseudoCost::upEstimate() const
{
    OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    double value = solution[columnNumber_];
    value = CoinMax(value, lower[columnNumber_]);
    value = CoinMin(value, upper[columnNumber_]);
    if (upper[columnNumber_] == lower[columnNumber_]) {
        // fixed
        return 0.0;
    }
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double below = floor(value + integerTolerance);
    double above = below + 1.0;
    if (above > upper[columnNumber_]) {
        above = below;
        below = above - 1;
    }
    double upCost = CoinMax((above - value) * upPseudoCost_, 0.0);
    return upCost;
}
// Return "down" estimate
double
CbcSimpleIntegerPseudoCost::downEstimate() const
{
    OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    double value = solution[columnNumber_];
    value = CoinMax(value, lower[columnNumber_]);
    value = CoinMin(value, upper[columnNumber_]);
    if (upper[columnNumber_] == lower[columnNumber_]) {
        // fixed
        return 0.0;
    }
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double below = floor(value + integerTolerance);
    double above = below + 1.0;
    if (above > upper[columnNumber_]) {
        above = below;
        below = above - 1;
    }
    double downCost = CoinMax((value - below) * downPseudoCost_, 0.0);
    return downCost;
}

