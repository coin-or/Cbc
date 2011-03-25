// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/9/2009-- carved out of CbcBranchActual

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
#include "CbcSimpleInteger.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"


//##############################################################################

/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
CbcSimpleInteger::CbcSimpleInteger ()
        : CbcObject(),
        originalLower_(0.0),
        originalUpper_(1.0),
        breakEven_(0.5),
        columnNumber_(-1),
        preferredWay_(0)
{
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
CbcSimpleInteger::CbcSimpleInteger ( CbcModel * model, int iColumn, double breakEven)
        : CbcObject(model)
{
    columnNumber_ = iColumn ;
    originalLower_ = model->solver()->getColLower()[columnNumber_] ;
    originalUpper_ = model->solver()->getColUpper()[columnNumber_] ;
    breakEven_ = breakEven;
    assert (breakEven_ > 0.0 && breakEven_ < 1.0);
    preferredWay_ = 0;
}


// Copy constructor
CbcSimpleInteger::CbcSimpleInteger ( const CbcSimpleInteger & rhs)
        : CbcObject(rhs)

{
    columnNumber_ = rhs.columnNumber_;
    originalLower_ = rhs.originalLower_;
    originalUpper_ = rhs.originalUpper_;
    breakEven_ = rhs.breakEven_;
    preferredWay_ = rhs.preferredWay_;
}

// Clone
CbcObject *
CbcSimpleInteger::clone() const
{
    return new CbcSimpleInteger(*this);
}

// Assignment operator
CbcSimpleInteger &
CbcSimpleInteger::operator=( const CbcSimpleInteger & rhs)
{
    if (this != &rhs) {
        CbcObject::operator=(rhs);
        columnNumber_ = rhs.columnNumber_;
        originalLower_ = rhs.originalLower_;
        originalUpper_ = rhs.originalUpper_;
        breakEven_ = rhs.breakEven_;
        preferredWay_ = rhs.preferredWay_;
    }
    return *this;
}

// Destructor
CbcSimpleInteger::~CbcSimpleInteger ()
{
}
// Construct an OsiSimpleInteger object
OsiSimpleInteger *
CbcSimpleInteger::osiObject() const
{
    OsiSimpleInteger * obj = new OsiSimpleInteger(columnNumber_,
            originalLower_, originalUpper_);
    obj->setPriority(priority());
    return obj;
}
double
CbcSimpleInteger::infeasibility(const OsiBranchingInformation * info,
                                int &preferredWay) const
{
    double value = info->solution_[columnNumber_];
    value = CoinMax(value, info->lower_[columnNumber_]);
    value = CoinMin(value, info->upper_[columnNumber_]);
    double nearest = floor(value + (1.0 - breakEven_));
    assert (breakEven_ > 0.0 && breakEven_ < 1.0);
    if (nearest > value)
        preferredWay = 1;
    else
        preferredWay = -1;
    if (preferredWay_)
        preferredWay = preferredWay_;
    double weight = fabs(value - nearest);
    // normalize so weight is 0.5 at break even
    if (nearest < value)
        weight = (0.5 / breakEven_) * weight;
    else
        weight = (0.5 / (1.0 - breakEven_)) * weight;
    if (fabs(value - nearest) <= info->integerTolerance_)
        return 0.0;
    else
        return weight;
}
double
CbcSimpleInteger::feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const
{
    double value = info->solution_[columnNumber_];
#ifdef COIN_DEVELOP
    if (fabs(value - floor(value + 0.5)) > 1.0e-5)
        printf("value for %d away from integer %g\n", columnNumber_, value);
#endif
    double newValue = CoinMax(value, info->lower_[columnNumber_]);
    newValue = CoinMin(newValue, info->upper_[columnNumber_]);
    newValue = floor(newValue + 0.5);
    solver->setColLower(columnNumber_, newValue);
    solver->setColUpper(columnNumber_, newValue);
    return fabs(value - newValue);
}

/* Create an OsiSolverBranch object

This returns NULL if branch not represented by bound changes
*/
OsiSolverBranch *
CbcSimpleInteger::solverBranch(OsiSolverInterface * /*solver*/,
                               const OsiBranchingInformation * info) const
{
    double value = info->solution_[columnNumber_];
    value = CoinMax(value, info->lower_[columnNumber_]);
    value = CoinMin(value, info->upper_[columnNumber_]);
    assert (info->upper_[columnNumber_] > info->lower_[columnNumber_]);
#ifndef NDEBUG
    double nearest = floor(value + 0.5);
    assert (fabs(value - nearest) > info->integerTolerance_);
#endif
    OsiSolverBranch * branch = new OsiSolverBranch();
    branch->addBranch(columnNumber_, value);
    return branch;
}
// Creates a branching object
CbcBranchingObject *
CbcSimpleInteger::createCbcBranch(OsiSolverInterface * /*solver*/,
                                  const OsiBranchingInformation * info, int way)
{
    CbcIntegerBranchingObject * branch = new CbcIntegerBranchingObject(model_, 0, -1, 0.5);
    fillCreateBranch(branch, info, way);
    return branch;
}
// Fills in a created branching object
void
CbcSimpleInteger::fillCreateBranch(CbcIntegerBranchingObject * branch, const OsiBranchingInformation * info, int way)
{
    branch->setOriginalObject(this);
    double value = info->solution_[columnNumber_];
    value = CoinMax(value, info->lower_[columnNumber_]);
    value = CoinMin(value, info->upper_[columnNumber_]);
    assert (info->upper_[columnNumber_] > info->lower_[columnNumber_]);
    if (!info->hotstartSolution_ && priority_ != -999) {
#ifndef NDEBUG
        double nearest = floor(value + 0.5);
        assert (fabs(value - nearest) > info->integerTolerance_);
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
    assert (value >= info->lower_[columnNumber_] &&
            value <= info->upper_[columnNumber_]);
    branch->fillPart(columnNumber_, way, value);
}
/* Column number if single column object -1 otherwise,
   so returns >= 0
   Used by heuristics
*/
int
CbcSimpleInteger::columnNumber() const
{
    return columnNumber_;
}
/* Reset variable bounds to their original values.

    Bounds may be tightened, so it may be good to be able to set this info in object.
*/
void
CbcSimpleInteger::resetBounds(const OsiSolverInterface * solver)
{
    originalLower_ = solver->getColLower()[columnNumber_] ;
    originalUpper_ = solver->getColUpper()[columnNumber_] ;
}

/*  Change column numbers after preprocessing
 */
void
CbcSimpleInteger::resetSequenceEtc(int /*numberColumns*/,
                                   const int * originalColumns)
{
    //assert (numberColumns>0);
    int iColumn;
#ifdef JJF_ZERO
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (columnNumber_ == originalColumns[iColumn])
            break;
    }
    assert (iColumn < numberColumns);
#else
    iColumn = originalColumns[columnNumber_];
    assert (iColumn >= 0);
#endif
    columnNumber_ = iColumn;
}
// This looks at solution and sets bounds to contain solution
/** More precisely: it first forces the variable within the existing
    bounds, and then tightens the bounds to fix the variable at the
    nearest integer value.
*/
void
CbcSimpleInteger::feasibleRegion()
{
    abort();
}

//##############################################################################

// Default Constructor
CbcIntegerBranchingObject::CbcIntegerBranchingObject()
        : CbcBranchingObject()
{
    down_[0] = 0.0;
    down_[1] = 0.0;
    up_[0] = 0.0;
    up_[1] = 0.0;
#ifdef FUNNY_BRANCHING
    variables_ = NULL;
    newBounds_ = NULL;
    numberExtraChangedBounds_ = 0;
#endif
}
// Useful constructor
CbcIntegerBranchingObject::CbcIntegerBranchingObject (CbcModel * model,
        int variable, int way , double value)
        : CbcBranchingObject(model, variable, way, value)
{
    int iColumn = variable;
    assert (model_->solver()->getNumCols() > 0);
    down_[0] = model_->solver()->getColLower()[iColumn];
    down_[1] = floor(value_);
    up_[0] = ceil(value_);
    up_[1] = model->getColUpper()[iColumn];
#ifdef FUNNY_BRANCHING
    variables_ = NULL;
    newBounds_ = NULL;
    numberExtraChangedBounds_ = 0;
#endif
}
// Does part of constructor
void
CbcIntegerBranchingObject::fillPart (int variable,
                                     int way , double value)
{
    //originalObject_=NULL;
    branchIndex_ = 0;
    value_ = value;
    numberBranches_ = 2;
    //model_= model;
    //originalCbcObject_=NULL;
    variable_ = variable;
    way_ = way;
    int iColumn = variable;
    down_[0] = model_->solver()->getColLower()[iColumn];
    down_[1] = floor(value_);
    up_[0] = ceil(value_);
    up_[1] = model_->getColUpper()[iColumn];
}
// Useful constructor for fixing
CbcIntegerBranchingObject::CbcIntegerBranchingObject (CbcModel * model,
        int variable, int way,
        double lowerValue,
        double upperValue)
        : CbcBranchingObject(model, variable, way, lowerValue)
{
    setNumberBranchesLeft(1);
    down_[0] = lowerValue;
    down_[1] = upperValue;
    up_[0] = lowerValue;
    up_[1] = upperValue;
#ifdef FUNNY_BRANCHING
    variables_ = NULL;
    newBounds_ = NULL;
    numberExtraChangedBounds_ = 0;
#endif
}


// Copy constructor
CbcIntegerBranchingObject::CbcIntegerBranchingObject ( const CbcIntegerBranchingObject & rhs) : CbcBranchingObject(rhs)
{
    down_[0] = rhs.down_[0];
    down_[1] = rhs.down_[1];
    up_[0] = rhs.up_[0];
    up_[1] = rhs.up_[1];
#ifdef FUNNY_BRANCHING
    numberExtraChangedBounds_ = rhs.numberExtraChangedBounds_;
    int size = numberExtraChangedBounds_ * (sizeof(double) + sizeof(int));
    char * temp = new char [size];
    newBounds_ = (double *) temp;
    variables_ = (int *) (newBounds_ + numberExtraChangedBounds_);

    int i ;
    for (i = 0; i < numberExtraChangedBounds_; i++) {
        variables_[i] = rhs.variables_[i];
        newBounds_[i] = rhs.newBounds_[i];
    }
#endif
}

// Assignment operator
CbcIntegerBranchingObject &
CbcIntegerBranchingObject::operator=( const CbcIntegerBranchingObject & rhs)
{
    if (this != &rhs) {
        CbcBranchingObject::operator=(rhs);
        down_[0] = rhs.down_[0];
        down_[1] = rhs.down_[1];
        up_[0] = rhs.up_[0];
        up_[1] = rhs.up_[1];
#ifdef FUNNY_BRANCHING
        delete [] newBounds_;
        numberExtraChangedBounds_ = rhs.numberExtraChangedBounds_;
        int size = numberExtraChangedBounds_ * (sizeof(double) + sizeof(int));
        char * temp = new char [size];
        newBounds_ = (double *) temp;
        variables_ = (int *) (newBounds_ + numberExtraChangedBounds_);

        int i ;
        for (i = 0; i < numberExtraChangedBounds_; i++) {
            variables_[i] = rhs.variables_[i];
            newBounds_[i] = rhs.newBounds_[i];
        }
#endif
    }
    return *this;
}
CbcBranchingObject *
CbcIntegerBranchingObject::clone() const
{
    return (new CbcIntegerBranchingObject(*this));
}


// Destructor
CbcIntegerBranchingObject::~CbcIntegerBranchingObject ()
{
    // for debugging threads
    way_ = -23456789;
#ifdef FUNNY_BRANCHING
    delete [] newBounds_;
#endif
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
CbcIntegerBranchingObject::branch()
{
    // for debugging threads
    if (way_ < -1 || way_ > 100000) {
        printf("way %d, left %d, iCol %d, variable %d\n",
               way_, numberBranchesLeft(),
               originalCbcObject_->columnNumber(), variable_);
        assert (way_ != -23456789);
    }
    decrementNumberBranchesLeft();
    if (down_[1] == -COIN_DBL_MAX)
        return 0.0;
    int iColumn = originalCbcObject_->columnNumber();
    assert (variable_ == iColumn);
    double olb, oub ;
    olb = model_->solver()->getColLower()[iColumn] ;
    oub = model_->solver()->getColUpper()[iColumn] ;
    //#define CBCSIMPLE_TIGHTEN_BOUNDS
#ifndef CBCSIMPLE_TIGHTEN_BOUNDS
#ifdef COIN_DEVELOP
    if (olb != down_[0] || oub != up_[1]) {
        if (way_ > 0)
            printf("branching up on var %d: [%g,%g] => [%g,%g] - other [%g,%g]\n",
                   iColumn, olb, oub, up_[0], up_[1], down_[0], down_[1]) ;
        else
            printf("branching down on var %d: [%g,%g] => [%g,%g] - other [%g,%g]\n",
                   iColumn, olb, oub, down_[0], down_[1], up_[0], up_[1]) ;
    }
#endif
#endif
    if (way_ < 0) {
#ifdef CBC_DEBUG
        { double olb, oub ;
            olb = model_->solver()->getColLower()[iColumn] ;
            oub = model_->solver()->getColUpper()[iColumn] ;
            printf("branching down on var %d: [%g,%g] => [%g,%g]\n",
                   iColumn, olb, oub, down_[0], down_[1]) ;
        }
#endif
#ifndef CBCSIMPLE_TIGHTEN_BOUNDS
        model_->solver()->setColLower(iColumn, down_[0]);
#else
        model_->solver()->setColLower(iColumn, CoinMax(down_[0], olb));
#endif
        model_->solver()->setColUpper(iColumn, down_[1]);
	//#define CBC_PRINT2
#ifdef CBC_PRINT2
        printf("%d branching down has bounds %g %g", iColumn, down_[0], down_[1]);
#endif
#ifdef FUNNY_BRANCHING
        // branch - do extra bounds
        for (int i = 0; i < numberExtraChangedBounds_; i++) {
            int variable = variables_[i];
            if ((variable&0x40000000) != 0) {
                // for going down
                int k = variable & 0x3fffffff;
                assert (k != iColumn);
                if ((variable&0x80000000) == 0) {
                    // lower bound changing
#ifdef CBC_PRINT2
                    printf(" extra for %d changes lower from %g to %g",
                           k, model_->solver()->getColLower()[k], newBounds_[i]);
#endif
                    model_->solver()->setColLower(k, newBounds_[i]);
                } else {
                    // upper bound changing
#ifdef CBC_PRINT2
                    printf(" extra for %d changes upper from %g to %g",
                           k, model_->solver()->getColUpper()[k], newBounds_[i]);
#endif
                    model_->solver()->setColUpper(k, newBounds_[i]);
                }
            }
        }
#endif
#ifdef CBC_PRINT2
        printf("\n");
#endif
        way_ = 1;
    } else {
#ifdef CBC_DEBUG
        { double olb, oub ;
            olb = model_->solver()->getColLower()[iColumn] ;
            oub = model_->solver()->getColUpper()[iColumn] ;
            printf("branching up on var %d: [%g,%g] => [%g,%g]\n",
                   iColumn, olb, oub, up_[0], up_[1]) ;
        }
#endif
        model_->solver()->setColLower(iColumn, up_[0]);
#ifndef CBCSIMPLE_TIGHTEN_BOUNDS
        model_->solver()->setColUpper(iColumn, up_[1]);
#else
        model_->solver()->setColUpper(iColumn, CoinMin(up_[1], oub));
#endif
#ifdef CBC_PRINT2
        printf("%d branching up has bounds %g %g", iColumn, up_[0], up_[1]);
#endif
#ifdef FUNNY_BRANCHING
        // branch - do extra bounds
        for (int i = 0; i < numberExtraChangedBounds_; i++) {
            int variable = variables_[i];
            if ((variable&0x40000000) == 0) {
                // for going up
                int k = variable & 0x3fffffff;
                assert (k != iColumn);
                if ((variable&0x80000000) == 0) {
                    // lower bound changing
#ifdef CBC_PRINT2
                    printf(" extra for %d changes lower from %g to %g",
                           k, model_->solver()->getColLower()[k], newBounds_[i]);
#endif
                    model_->solver()->setColLower(k, newBounds_[i]);
                } else {
                    // upper bound changing
#ifdef CBC_PRINT2
                    printf(" extra for %d changes upper from %g to %g",
                           k, model_->solver()->getColUpper()[k], newBounds_[i]);
#endif
                    model_->solver()->setColUpper(k, newBounds_[i]);
                }
            }
        }
#endif
#ifdef CBC_PRINT2
        printf("\n");
#endif
        way_ = -1;	  // Swap direction
    }
    double nlb = model_->solver()->getColLower()[iColumn];
    double nub = model_->solver()->getColUpper()[iColumn];
    if (nlb < olb) {
#ifdef CBC_PRINT2
        printf("bad lb change for column %d from %g to %g\n", iColumn, olb, nlb);
#endif
	//abort();
        model_->solver()->setColLower(iColumn, CoinMin(olb, nub));
        nlb = olb;
    }
    if (nub > oub) {
#ifdef CBC_PRINT2
        printf("bad ub change for column %d from %g to %g\n", iColumn, oub, nub);
#endif
	//abort();
        model_->solver()->setColUpper(iColumn, CoinMax(oub, nlb));
    }
#ifdef CBC_PRINT2
    if (nlb < olb + 1.0e-8 && nub > oub - 1.0e-8 && false)
        printf("bad null change for column %d - bounds %g,%g\n", iColumn, olb, oub);
#endif
    return 0.0;
}
/* Update bounds in solver as in 'branch' and update given bounds.
   branchState is -1 for 'down' +1 for 'up' */
void
CbcIntegerBranchingObject::fix(OsiSolverInterface * /*solver*/,
                               double * lower, double * upper,
                               int branchState) const
{
    int iColumn = originalCbcObject_->columnNumber();
    assert (variable_ == iColumn);
    if (branchState < 0) {
        model_->solver()->setColLower(iColumn, down_[0]);
        lower[iColumn] = down_[0];
        model_->solver()->setColUpper(iColumn, down_[1]);
        upper[iColumn] = down_[1];
    } else {
        model_->solver()->setColLower(iColumn, up_[0]);
        lower[iColumn] = up_[0];
        model_->solver()->setColUpper(iColumn, up_[1]);
        upper[iColumn] = up_[1];
    }
}
// Change (tighten) bounds in object to reflect bounds in solver.
// Return true if now fixed
bool 
CbcIntegerBranchingObject::tighten(OsiSolverInterface * solver) 
{
    double lower = solver->getColLower()[variable_];
    double upper = solver->getColUpper()[variable_];
    assert (upper>lower);
    down_[0] = CoinMax(down_[0],lower);
    up_[0] = CoinMax(up_[0],lower);
    down_[1] = CoinMin(down_[1],upper);
    up_[1] = CoinMin(up_[1],upper);
    return (down_[0]==up_[1]);
}
#ifdef FUNNY_BRANCHING
// Deactivate bounds for branching
void
CbcIntegerBranchingObject::deactivate()
{
    down_[1] = -COIN_DBL_MAX;
}
int
CbcIntegerBranchingObject::applyExtraBounds(int iColumn, double lower, double upper, int way)
{
    // branch - do bounds

    int i;
    int found = 0;
    if (variable_ == iColumn) {
        printf("odd applyExtra %d\n", iColumn);
        if (way < 0) {
            down_[0] = CoinMax(lower, down_[0]);
            down_[1] = CoinMin(upper, down_[1]);
            assert (down_[0] <= down_[1]);
        } else {
            up_[0] = CoinMax(lower, up_[0]);
            up_[1] = CoinMin(upper, up_[1]);
            assert (up_[0] <= up_[1]);
        }
        return 0;
    }
    int check = (way < 0) ? 0x40000000 : 0;
    double newLower = lower;
    double newUpper = upper;
    for (i = 0; i < numberExtraChangedBounds_; i++) {
        int variable = variables_[i];
        if ((variable&0x40000000) == check) {
            int k = variable & 0x3fffffff;
            if (k == iColumn) {
                if ((variable&0x80000000) == 0) {
                    // lower bound changing
                    found |= 1;
                    newBounds_[i] = CoinMax(lower, newBounds_[i]);
                    newLower = newBounds_[i];
                } else {
                    // upper bound changing
                    found |= 2;
                    newBounds_[i] = CoinMin(upper, newBounds_[i]);
                    newUpper = newBounds_[i];
                }
            }
        }
    }
    int nAdd = 0;
    if ((found&2) == 0) {
        // need to add new upper
        nAdd++;
    }
    if ((found&1) == 0) {
        // need to add new lower
        nAdd++;
    }
    if (nAdd) {
        int size = (numberExtraChangedBounds_ + nAdd) * (sizeof(double) + sizeof(int));
        char * temp = new char [size];
        double * newBounds = (double *) temp;
        int * variables = (int *) (newBounds + numberExtraChangedBounds_ + nAdd);

        int i ;
        for (i = 0; i < numberExtraChangedBounds_; i++) {
            variables[i] = variables_[i];
            newBounds[i] = newBounds_[i];
        }
        delete [] newBounds_;
        newBounds_ = newBounds;
        variables_ = variables;
        if ((found&2) == 0) {
            // need to add new upper
            int variable = iColumn | 0x80000000;
            variables_[numberExtraChangedBounds_] = variable;
            newBounds_[numberExtraChangedBounds_++] = newUpper;
        }
        if ((found&1) == 0) {
            // need to add new lower
            int variable = iColumn;
            variables_[numberExtraChangedBounds_] = variable;
            newBounds_[numberExtraChangedBounds_++] = newLower;
        }
    }

    return (newUpper >= newLower) ? 0 : 1;
}
#endif
// Print what would happen
void
CbcIntegerBranchingObject::print()
{
    int iColumn = originalCbcObject_->columnNumber();
    assert (variable_ == iColumn);
    if (way_ < 0) {
        {
            double olb, oub ;
            olb = model_->solver()->getColLower()[iColumn] ;
            oub = model_->solver()->getColUpper()[iColumn] ;
            printf("CbcInteger would branch down on var %d (int var %d): [%g,%g] => [%g,%g]\n",
                   iColumn, variable_, olb, oub, down_[0], down_[1]) ;
        }
    } else {
        {
            double olb, oub ;
            olb = model_->solver()->getColLower()[iColumn] ;
            oub = model_->solver()->getColUpper()[iColumn] ;
            printf("CbcInteger would branch up on var %d (int var %d): [%g,%g] => [%g,%g]\n",
                   iColumn, variable_, olb, oub, up_[0], up_[1]) ;
        }
    }
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
CbcIntegerBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
    const CbcIntegerBranchingObject* br =
        dynamic_cast<const CbcIntegerBranchingObject*>(brObj);
    assert(br);
    double* thisBd = way_ < 0 ? down_ : up_;
    const double* otherBd = br->way_ < 0 ? br->down_ : br->up_;
    return CbcCompareRanges(thisBd, otherBd, replaceIfOverlap);
}

