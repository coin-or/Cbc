// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/13/2009-- carved out of CbcBranchCut

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
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchCut.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"
#include "CbcBranchAllDifferent.hpp"

/** Default Constructor
*/
CbcBranchAllDifferent::CbcBranchAllDifferent ()
        : CbcBranchCut(),
        numberInSet_(0),
        which_(NULL)
{
}

/* Useful constructor - passed set of variables
*/
CbcBranchAllDifferent::CbcBranchAllDifferent (CbcModel * model, int numberInSet,
        const int * members)
        : CbcBranchCut(model)
{
    numberInSet_ = numberInSet;
    which_ = CoinCopyOfArray(members, numberInSet_);
}
// Copy constructor
CbcBranchAllDifferent::CbcBranchAllDifferent ( const CbcBranchAllDifferent & rhs)
        : CbcBranchCut(rhs)
{
    numberInSet_ = rhs.numberInSet_;
    which_ = CoinCopyOfArray(rhs.which_, numberInSet_);
}

// Clone
CbcObject *
CbcBranchAllDifferent::clone() const
{
    return new CbcBranchAllDifferent(*this);
}

// Assignment operator
CbcBranchAllDifferent &
CbcBranchAllDifferent::operator=( const CbcBranchAllDifferent & rhs)
{
    if (this != &rhs) {
        CbcBranchCut::operator=(rhs);
        delete [] which_;
        numberInSet_ = rhs.numberInSet_;
        which_ = CoinCopyOfArray(rhs.which_, numberInSet_);
    }
    return *this;
}

// Destructor
CbcBranchAllDifferent::~CbcBranchAllDifferent ()
{
    delete [] which_;
}
CbcBranchingObject *
CbcBranchAllDifferent::createCbcBranch(OsiSolverInterface * /*solver*/
                                       , const OsiBranchingInformation * /*info*/,
                                       int /*way*/)
{
    // by default way must be -1
    //assert (way==-1);
    const double * solution = model_->testSolution();
    double * values = new double[numberInSet_];
    int * which = new int[numberInSet_];
    int i;
    for (i = 0; i < numberInSet_; i++) {
        int iColumn = which_[i];
        values[i] = solution[iColumn];
        which[i] = iColumn;
    }
    CoinSort_2(values, values + numberInSet_, which);
    double last = -1.0;
    double closest = 1.0;
    int worst = -1;
    for (i = 0; i < numberInSet_; i++) {
        if (values[i] - last < closest) {
            closest = values[i] - last;
            worst = i - 1;
        }
        last = values[i];
    }
    assert (closest <= 0.99999);
    OsiRowCut down;
    down.setLb(-COIN_DBL_MAX);
    down.setUb(-1.0);
    int pair[2];
    double elements[] = {1.0, -1.0};
    pair[0] = which[worst];
    pair[1] = which[worst+1];
    delete [] values;
    delete [] which;
    down.setRow(2, pair, elements);
    // up is same - just with rhs changed
    OsiRowCut up = down;
    up.setLb(1.0);
    up.setUb(COIN_DBL_MAX);
    // Say is not a fix type branch
    CbcCutBranchingObject * newObject =
        new CbcCutBranchingObject(model_, down, up, false);
    if (model_->messageHandler()->logLevel() > 1)
        printf("creating cut in CbcBranchCut\n");
    return newObject;
}
double
CbcBranchAllDifferent::infeasibility(const OsiBranchingInformation * /*info*/,
                                     int &preferredWay) const
{
    preferredWay = -1;
    //OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    //const double * lower = solver->getColLower();
    //const double * upper = solver->getColUpper();
    double * values = new double[numberInSet_];
    int i;
    for (i = 0; i < numberInSet_; i++) {
        int iColumn = which_[i];
        values[i] = solution[iColumn];
    }
    std::sort(values, values + numberInSet_);
    double last = -1.0;
    double closest = 1.0;
    for (i = 0; i < numberInSet_; i++) {
        if (values[i] - last < closest) {
            closest = values[i] - last;
        }
        last = values[i];
    }
    delete [] values;
    if (closest > 0.99999)
        return 0.0;
    else
        return 0.5*(1.0 - closest);
}

