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
#include "CbcNWayBranchingObject.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

// Default Constructor
CbcNWayBranchingObject::CbcNWayBranchingObject()
        : CbcBranchingObject()
{
    order_ = NULL;
    object_ = NULL;
    numberInSet_ = 0;
    way_ = 0;
}

// Useful constructor
CbcNWayBranchingObject::CbcNWayBranchingObject (CbcModel * model,
        const CbcNWay * nway,
        int number, const int * order)
        : CbcBranchingObject(model, nway->id(), -1, 0.5)
{
    numberBranches_ = number;
    order_ = new int [number];
    object_ = nway;
    numberInSet_ = number;
    memcpy(order_, order, number*sizeof(int));
}

// Copy constructor
CbcNWayBranchingObject::CbcNWayBranchingObject ( const CbcNWayBranchingObject & rhs) : CbcBranchingObject(rhs)
{
    numberInSet_ = rhs.numberInSet_;
    object_ = rhs.object_;
    if (numberInSet_) {
        order_ = new int [numberInSet_];
        memcpy(order_, rhs.order_, numberInSet_*sizeof(int));
    } else {
        order_ = NULL;
    }
}

// Assignment operator
CbcNWayBranchingObject &
CbcNWayBranchingObject::operator=( const CbcNWayBranchingObject & rhs)
{
    if (this != &rhs) {
        CbcBranchingObject::operator=(rhs);
        object_ = rhs.object_;
        delete [] order_;
        numberInSet_ = rhs.numberInSet_;
        if (numberInSet_) {
            order_ = new int [numberInSet_];
            memcpy(order_, rhs.order_, numberInSet_*sizeof(int));
        } else {
            order_ = NULL;
        }
    }
    return *this;
}
CbcBranchingObject *
CbcNWayBranchingObject::clone() const
{
    return (new CbcNWayBranchingObject(*this));
}


// Destructor
CbcNWayBranchingObject::~CbcNWayBranchingObject ()
{
    delete [] order_;
}
double
CbcNWayBranchingObject::branch()
{
    int which = branchIndex_;
    branchIndex_++;
    assert (numberBranchesLeft() >= 0);
    if (which == 0) {
        // first branch so way_ may mean something
        assert (way_ == -1 || way_ == 1);
        if (way_ == -1)
            which++;
    } else if (which == 1) {
        // second branch so way_ may mean something
        assert (way_ == -1 || way_ == 1);
        if (way_ == -1)
            which--;
        // switch way off
        way_ = 0;
    }
    const double * lower = model_->solver()->getColLower();
    const double * upper = model_->solver()->getColUpper();
    const int * members = object_->members();
    for (int j = 0; j < numberInSet_; j++) {
        int iSequence = order_[j];
        int iColumn = members[iSequence];
        if (j != which) {
            model_->solver()->setColUpper(iColumn, lower[iColumn]);
            //model_->solver()->setColLower(iColumn,lower[iColumn]);
            assert (lower[iColumn] > -1.0e20);
            // apply any consequences
            object_->applyConsequence(iSequence, -9999);
        } else {
            model_->solver()->setColLower(iColumn, upper[iColumn]);
            //model_->solver()->setColUpper(iColumn,upper[iColumn]);
#ifdef FULL_PRINT
            printf("Up Fix %d to %g\n", iColumn, upper[iColumn]);
#endif
            assert (upper[iColumn] < 1.0e20);
            // apply any consequences
            object_->applyConsequence(iSequence, 9999);
        }
    }
    return 0.0;
}
void
CbcNWayBranchingObject::print()
{
    printf("NWay - Up Fix ");
    const int * members = object_->members();
    for (int j = 0; j < way_; j++) {
        int iColumn = members[order_[j]];
        printf("%d ", iColumn);
    }
    printf("\n");
}

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcNWayBranchingObject::compareOriginalObject
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
CbcNWayBranchingObject::compareBranchingObject
(const CbcBranchingObject* /*brObj*/, const bool /*replaceIfOverlap*/)
{
    throw("must implement");
}
