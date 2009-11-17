// Edwin 11/17/2009-- carved out of CbcBranchDynamic
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
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

// Default Constructor
CbcDynamicPseudoCostBranchingObject::CbcDynamicPseudoCostBranchingObject()
        : CbcIntegerBranchingObject()
{
    changeInGuessed_ = 1.0e-5;
    object_ = NULL;
}

// Useful constructor
CbcDynamicPseudoCostBranchingObject::CbcDynamicPseudoCostBranchingObject (CbcModel * model,
        int variable,
        int way , double value,
        CbcSimpleIntegerDynamicPseudoCost * object)
        : CbcIntegerBranchingObject(model, variable, way, value)
{
    changeInGuessed_ = 1.0e-5;
    object_ = object;
}
// Does part of work for constructor
void
CbcDynamicPseudoCostBranchingObject::fillPart (int variable,
        int way , double value,
        CbcSimpleIntegerDynamicPseudoCost * object)
{
    CbcIntegerBranchingObject::fillPart(variable, way, value);
    changeInGuessed_ = 1.0e-5;
    object_ = object;
}
// Useful constructor for fixing
CbcDynamicPseudoCostBranchingObject::CbcDynamicPseudoCostBranchingObject (CbcModel * model,
        int variable, int way,
        double lowerValue,
        double /*upperValue*/)
        : CbcIntegerBranchingObject(model, variable, way, lowerValue)
{
    changeInGuessed_ = 1.0e100;
    object_ = NULL;
}


// Copy constructor
CbcDynamicPseudoCostBranchingObject::CbcDynamicPseudoCostBranchingObject (
    const CbcDynamicPseudoCostBranchingObject & rhs)
        : CbcIntegerBranchingObject(rhs)
{
    changeInGuessed_ = rhs.changeInGuessed_;
    object_ = rhs.object_;
}

// Assignment operator
CbcDynamicPseudoCostBranchingObject &
CbcDynamicPseudoCostBranchingObject::operator=( const CbcDynamicPseudoCostBranchingObject & rhs)
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
    assert (info.possibleBranch == this);
    info.upMovement = object_->upDynamicPseudoCost() * (ceil(value_) - value_);
    info.downMovement = object_->downDynamicPseudoCost() * (value_ - floor(value_));
    info.numIntInfeasUp  -= static_cast<int> (object_->sumUpDecrease() /
                            (1.0e-12 + static_cast<double> (object_->numberTimesUp())));
    info.numIntInfeasUp = CoinMax(info.numIntInfeasUp, 0);
    info.numObjInfeasUp = 0;
    info.finishedUp = false;
    info.numItersUp = 0;
    info.numIntInfeasDown  -= static_cast<int> (object_->sumDownDecrease() /
                              (1.0e-12 + static_cast<double> (object_->numberTimesDown())));
    info.numIntInfeasDown = CoinMax(info.numIntInfeasDown, 0);
    info.numObjInfeasDown = 0;
    info.finishedDown = false;
    info.numItersDown = 0;
    info.fix = 0;
    if (object_->numberTimesUp() < object_->numberBeforeTrust() +
            2*object_->numberTimesUpInfeasible() ||
            object_->numberTimesDown() < object_->numberBeforeTrust() +
            2*object_->numberTimesDownInfeasible()) {
        return 0;
    } else {
        return 1;
    }
}
