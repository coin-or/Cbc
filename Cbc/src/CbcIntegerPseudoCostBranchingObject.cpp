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
#include "CbcIntegerPseudoCostBranchingObject.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

//##############################################################################

// Default Constructor
CbcIntegerPseudoCostBranchingObject::CbcIntegerPseudoCostBranchingObject()
        : CbcIntegerBranchingObject()
{
    changeInGuessed_ = 1.0e-5;
}

// Useful constructor
CbcIntegerPseudoCostBranchingObject::CbcIntegerPseudoCostBranchingObject (CbcModel * model,
        int variable, int way , double value)
        : CbcIntegerBranchingObject(model, variable, way, value)
{
}
// Useful constructor for fixing
CbcIntegerPseudoCostBranchingObject::CbcIntegerPseudoCostBranchingObject (CbcModel * model,
        int variable, int way,
        double lowerValue,
        double /*upperValue*/)
        : CbcIntegerBranchingObject(model, variable, way, lowerValue)
{
    changeInGuessed_ = 1.0e100;
}


// Copy constructor
CbcIntegerPseudoCostBranchingObject::CbcIntegerPseudoCostBranchingObject (
    const CbcIntegerPseudoCostBranchingObject & rhs)
        : CbcIntegerBranchingObject(rhs)
{
    changeInGuessed_ = rhs.changeInGuessed_;
}

// Assignment operator
CbcIntegerPseudoCostBranchingObject &
CbcIntegerPseudoCostBranchingObject::operator=( const CbcIntegerPseudoCostBranchingObject & rhs)
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
CbcIntegerPseudoCostBranchingObject::~CbcIntegerPseudoCostBranchingObject ()
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
CbcIntegerPseudoCostBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
    const CbcIntegerPseudoCostBranchingObject* br =
        dynamic_cast<const CbcIntegerPseudoCostBranchingObject*>(brObj);
    assert(br);
    double* thisBd = way_ < 0 ? down_ : up_;
    const double* otherBd = br->way_ < 0 ? br->down_ : br->up_;
    return CbcCompareRanges(thisBd, otherBd, replaceIfOverlap);
}
