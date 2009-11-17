/* $Id$ */
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

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchLotsize.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"


// Default Constructor
CbcLotsizeBranchingObject::CbcLotsizeBranchingObject()
        : CbcBranchingObject()
{
    down_[0] = 0.0;
    down_[1] = 0.0;
    up_[0] = 0.0;
    up_[1] = 0.0;
}

// Useful constructor
CbcLotsizeBranchingObject::CbcLotsizeBranchingObject (CbcModel * model,
        int variable, int way , double value,
        const CbcLotsize * lotsize)
        : CbcBranchingObject(model, variable, way, value)
{
    int iColumn = lotsize->modelSequence();
    assert (variable == iColumn);
    down_[0] = model_->solver()->getColLower()[iColumn];
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
    lotsize->floorCeiling(down_[1], up_[0], value, integerTolerance);
    up_[1] = model->getColUpper()[iColumn];
}
// Useful constructor for fixing
CbcLotsizeBranchingObject::CbcLotsizeBranchingObject (CbcModel * model,
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
}


// Copy constructor
CbcLotsizeBranchingObject::CbcLotsizeBranchingObject ( const CbcLotsizeBranchingObject & rhs) : CbcBranchingObject(rhs)
{
    down_[0] = rhs.down_[0];
    down_[1] = rhs.down_[1];
    up_[0] = rhs.up_[0];
    up_[1] = rhs.up_[1];
}

// Assignment operator
CbcLotsizeBranchingObject &
CbcLotsizeBranchingObject::operator=( const CbcLotsizeBranchingObject & rhs)
{
    if (this != &rhs) {
        CbcBranchingObject::operator=(rhs);
        down_[0] = rhs.down_[0];
        down_[1] = rhs.down_[1];
        up_[0] = rhs.up_[0];
        up_[1] = rhs.up_[1];
    }
    return *this;
}
CbcBranchingObject *
CbcLotsizeBranchingObject::clone() const
{
    return (new CbcLotsizeBranchingObject(*this));
}


// Destructor
CbcLotsizeBranchingObject::~CbcLotsizeBranchingObject ()
{
}

/*
  Perform a branch by adjusting the bounds of the specified variable. Note
  that each arm of the branch advances the object to the next arm by
  advancing the value of way_.

  Providing new values for the variable's lower and upper bounds for each
  branching direction gives a little bit of additional flexibility and will
  be easily extensible to multi-way branching.
*/
double
CbcLotsizeBranchingObject::branch()
{
    decrementNumberBranchesLeft();
    int iColumn = variable_;
    if (way_ < 0) {
#ifdef CBC_DEBUG
        { double olb, oub ;
            olb = model_->solver()->getColLower()[iColumn] ;
            oub = model_->solver()->getColUpper()[iColumn] ;
            printf("branching down on var %d: [%g,%g] => [%g,%g]\n",
                   iColumn, olb, oub, down_[0], down_[1]) ;
        }
#endif
        model_->solver()->setColLower(iColumn, down_[0]);
        model_->solver()->setColUpper(iColumn, down_[1]);
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
        model_->solver()->setColUpper(iColumn, up_[1]);
        way_ = -1;	  // Swap direction
    }
    return 0.0;
}
// Print
void
CbcLotsizeBranchingObject::print()
{
    int iColumn = variable_;
    if (way_ < 0) {
        {
            double olb, oub ;
            olb = model_->solver()->getColLower()[iColumn] ;
            oub = model_->solver()->getColUpper()[iColumn] ;
            printf("branching down on var %d: [%g,%g] => [%g,%g]\n",
                   iColumn, olb, oub, down_[0], down_[1]) ;
        }
    } else {
        {
            double olb, oub ;
            olb = model_->solver()->getColLower()[iColumn] ;
            oub = model_->solver()->getColUpper()[iColumn] ;
            printf("branching up on var %d: [%g,%g] => [%g,%g]\n",
                   iColumn, olb, oub, up_[0], up_[1]) ;
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
CbcLotsizeBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
    const CbcLotsizeBranchingObject* br =
        dynamic_cast<const CbcLotsizeBranchingObject*>(brObj);
    assert(br);
    double* thisBd = way_ == -1 ? down_ : up_;
    const double* otherBd = br->way_ == -1 ? br->down_ : br->up_;
    return CbcCompareRanges(thisBd, otherBd, replaceIfOverlap);
}
