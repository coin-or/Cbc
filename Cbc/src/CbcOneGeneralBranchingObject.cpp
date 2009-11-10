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
#include "CbcOneGeneralBranchingObject.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

#ifdef COIN_HAS_CLP
// Default Constructor
CbcOneGeneralBranchingObject::CbcOneGeneralBranchingObject()
        : CbcBranchingObject(),
        object_(NULL),
        whichOne_(-1)
{
    //printf("CbcOneGeneral %x default constructor\n",this);
}

// Useful constructor
CbcOneGeneralBranchingObject::CbcOneGeneralBranchingObject (CbcModel * model,
        CbcGeneralBranchingObject * object,
        int whichOne)
        : CbcBranchingObject(model, -1, -1, 0.5),
        object_(object),
        whichOne_(whichOne)
{
    //printf("CbcOneGeneral %x useful constructor object %x %d left\n",this,
    //	 object_,object_->numberSubLeft_);
    numberBranches_ = 1;
}

// Copy constructor
CbcOneGeneralBranchingObject::CbcOneGeneralBranchingObject ( const CbcOneGeneralBranchingObject & rhs)
        : CbcBranchingObject(rhs),
        object_(rhs.object_),
        whichOne_(rhs.whichOne_)
{
}

// Assignment operator
CbcOneGeneralBranchingObject &
CbcOneGeneralBranchingObject::operator=( const CbcOneGeneralBranchingObject & rhs)
{
    if (this != &rhs) {
        CbcBranchingObject::operator=(rhs);
        object_ = rhs.object_;
        whichOne_ = rhs.whichOne_;
    }
    return *this;
}
CbcBranchingObject *
CbcOneGeneralBranchingObject::clone() const
{
    return (new CbcOneGeneralBranchingObject(*this));
}


// Destructor
CbcOneGeneralBranchingObject::~CbcOneGeneralBranchingObject ()
{
    //printf("CbcOneGeneral %x destructor object %x %d left\n",this,
    // object_,object_->numberSubLeft_);
    assert (object_->numberSubLeft_ > 0 &&
            object_->numberSubLeft_ < 1000000);
    if (!object_->decrementNumberLeft()) {
        // printf("CbcGeneral %x yy destructor\n",object_);
        delete object_;
    }
}
double
CbcOneGeneralBranchingObject::branch()
{
    assert (numberBranchesLeft());
    decrementNumberBranchesLeft();
    assert (!numberBranchesLeft());
    object_->setWhichNode(whichOne_);
    object_->branch();
    return 0.0;
}
/* Double checks in case node can change its mind!
   Can change objective etc */
void
CbcOneGeneralBranchingObject::checkIsCutoff(double /*cutoff*/)
{
    assert (numberBranchesLeft());
}
// Print what would happen
void
CbcOneGeneralBranchingObject::print()
{
    //printf("CbcOneGeneralObject has 1 subproblem\n");
}
/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcOneGeneralBranchingObject::compareOriginalObject
(const CbcBranchingObject* /*brObj*/) const
{
    throw("must implement");
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
CbcOneGeneralBranchingObject::compareBranchingObject
(const CbcBranchingObject* /*brObj*/, const bool /*replaceIfOverlap*/)
{
    throw("must implement");
}

#endif
