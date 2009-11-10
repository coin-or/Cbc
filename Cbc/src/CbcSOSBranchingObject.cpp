// Edwin 11/10/2009-- carved out of CbcBranchActual
//##############################################################################
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
#include "CbcSOSBranchingObject.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

// Default Constructor
CbcSOSBranchingObject::CbcSOSBranchingObject()
        : CbcBranchingObject(),
        firstNonzero_(-1),
        lastNonzero_(-1)
{
    set_ = NULL;
    separator_ = 0.0;
}

// Useful constructor
CbcSOSBranchingObject::CbcSOSBranchingObject (CbcModel * model,
        const CbcSOS * set,
        int way ,
        double separator)
        : CbcBranchingObject(model, set->id(), way, 0.5)
{
    set_ = set;
    separator_ = separator;
    computeNonzeroRange();
}

// Copy constructor
CbcSOSBranchingObject::CbcSOSBranchingObject (const CbcSOSBranchingObject & rhs)
        : CbcBranchingObject(rhs),
        firstNonzero_(rhs.firstNonzero_),
        lastNonzero_(rhs.lastNonzero_)
{
    set_ = rhs.set_;
    separator_ = rhs.separator_;
}

// Assignment operator
CbcSOSBranchingObject &
CbcSOSBranchingObject::operator=( const CbcSOSBranchingObject & rhs)
{
    if (this != &rhs) {
        CbcBranchingObject::operator=(rhs);
        set_ = rhs.set_;
        separator_ = rhs.separator_;
        firstNonzero_ = rhs.firstNonzero_;
        lastNonzero_ = rhs.lastNonzero_;
    }
    return *this;
}
CbcBranchingObject *
CbcSOSBranchingObject::clone() const
{
    return (new CbcSOSBranchingObject(*this));
}


// Destructor
CbcSOSBranchingObject::~CbcSOSBranchingObject ()
{
}

void
CbcSOSBranchingObject::computeNonzeroRange()
{
    const int numberMembers = set_->numberMembers();
    const double * weights = set_->weights();
    int i = 0;
    if (way_ < 0) {
        for ( i = 0; i < numberMembers; i++) {
            if (weights[i] > separator_)
                break;
        }
        assert (i < numberMembers);
        firstNonzero_ = 0;
        lastNonzero_ = i;
    } else {
        for ( i = 0; i < numberMembers; i++) {
            if (weights[i] >= separator_)
                break;
        }
        assert (i < numberMembers);
        firstNonzero_ = i;
        lastNonzero_ = numberMembers;
    }
}

double
CbcSOSBranchingObject::branch()
{
    decrementNumberBranchesLeft();
    int numberMembers = set_->numberMembers();
    const int * which = set_->members();
    const double * weights = set_->weights();
    OsiSolverInterface * solver = model_->solver();
    //const double * lower = solver->getColLower();
    //const double * upper = solver->getColUpper();
    // *** for way - up means fix all those in down section
    if (way_ < 0) {
        int i;
        for ( i = 0; i < numberMembers; i++) {
            if (weights[i] > separator_)
                break;
        }
        assert (i < numberMembers);
        for (; i < numberMembers; i++)
            solver->setColUpper(which[i], 0.0);
        way_ = 1;	  // Swap direction
    } else {
        int i;
        for ( i = 0; i < numberMembers; i++) {
            if (weights[i] >= separator_)
                break;
            else
                solver->setColUpper(which[i], 0.0);
        }
        assert (i < numberMembers);
        way_ = -1;	  // Swap direction
    }
    computeNonzeroRange();
    return 0.0;
}
/* Update bounds in solver as in 'branch' and update given bounds.
   branchState is -1 for 'down' +1 for 'up' */
void
CbcSOSBranchingObject::fix(OsiSolverInterface * solver,
                           double * /*lower*/, double * upper,
                           int branchState) const
{
    int numberMembers = set_->numberMembers();
    const int * which = set_->members();
    const double * weights = set_->weights();
    //const double * lower = solver->getColLower();
    //const double * upper = solver->getColUpper();
    // *** for way - up means fix all those in down section
    if (branchState < 0) {
        int i;
        for ( i = 0; i < numberMembers; i++) {
            if (weights[i] > separator_)
                break;
        }
        assert (i < numberMembers);
        for (; i < numberMembers; i++) {
            solver->setColUpper(which[i], 0.0);
            upper[which[i]] = 0.0;
        }
    } else {
        int i;
        for ( i = 0; i < numberMembers; i++) {
            if (weights[i] >= separator_) {
                break;
            } else {
                solver->setColUpper(which[i], 0.0);
                upper[which[i]] = 0.0;
            }
        }
        assert (i < numberMembers);
    }
}
// Print what would happen
void
CbcSOSBranchingObject::print()
{
    int numberMembers = set_->numberMembers();
    const int * which = set_->members();
    const double * weights = set_->weights();
    OsiSolverInterface * solver = model_->solver();
    //const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    int first = numberMembers;
    int last = -1;
    int numberFixed = 0;
    int numberOther = 0;
    int i;
    for ( i = 0; i < numberMembers; i++) {
        double bound = upper[which[i]];
        if (bound) {
            first = CoinMin(first, i);
            last = CoinMax(last, i);
        }
    }
    // *** for way - up means fix all those in down section
    if (way_ < 0) {
        printf("SOS Down");
        for ( i = 0; i < numberMembers; i++) {
            double bound = upper[which[i]];
            if (weights[i] > separator_)
                break;
            else if (bound)
                numberOther++;
        }
        assert (i < numberMembers);
        for (; i < numberMembers; i++) {
            double bound = upper[which[i]];
            if (bound)
                numberFixed++;
        }
    } else {
        printf("SOS Up");
        for ( i = 0; i < numberMembers; i++) {
            double bound = upper[which[i]];
            if (weights[i] >= separator_)
                break;
            else if (bound)
                numberFixed++;
        }
        assert (i < numberMembers);
        for (; i < numberMembers; i++) {
            double bound = upper[which[i]];
            if (bound)
                numberOther++;
        }
    }
    printf(" - at %g, free range %d (%g) => %d (%g), %d would be fixed, %d other way\n",
           separator_, which[first], weights[first], which[last], weights[last], numberFixed, numberOther);
}

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcSOSBranchingObject::compareOriginalObject
(const CbcBranchingObject* brObj) const
{
    const CbcSOSBranchingObject* br =
        dynamic_cast<const CbcSOSBranchingObject*>(brObj);
    assert(br);
    const CbcSOS* s0 = set_;
    const CbcSOS* s1 = br->set_;
    if (s0->sosType() != s1->sosType()) {
        return s0->sosType() - s1->sosType();
    }
    if (s0->numberMembers() != s1->numberMembers()) {
        return s0->numberMembers() - s1->numberMembers();
    }
    const int memberCmp = memcmp(s0->members(), s1->members(),
                                 s0->numberMembers() * sizeof(int));
    if (memberCmp != 0) {
        return memberCmp;
    }
    return memcmp(s0->weights(), s1->weights(),
                  s0->numberMembers() * sizeof(double));
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
CbcSOSBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
    const CbcSOSBranchingObject* br =
        dynamic_cast<const CbcSOSBranchingObject*>(brObj);
    assert(br);
    if (firstNonzero_ < br->firstNonzero_) {
        if (lastNonzero_ >= br->lastNonzero_) {
            return CbcRangeSuperset;
        } else if (lastNonzero_ <= br->firstNonzero_) {
            return CbcRangeDisjoint;
        } else {
            // overlap
            if (replaceIfOverlap) {
                firstNonzero_ = br->firstNonzero_;
            }
            return CbcRangeOverlap;
        }
    } else if (firstNonzero_ > br->firstNonzero_) {
        if (lastNonzero_ <= br->lastNonzero_) {
            return CbcRangeSubset;
        } else if (firstNonzero_ >= br->lastNonzero_) {
            return CbcRangeDisjoint;
        } else {
            // overlap
            if (replaceIfOverlap) {
                lastNonzero_ = br->lastNonzero_;
            }
            return CbcRangeOverlap;
        }
    } else {
        if (lastNonzero_ == br->lastNonzero_) {
            return CbcRangeSame;
        }
        return lastNonzero_ < br->lastNonzero_ ? CbcRangeSubset : CbcRangeSuperset;
    }
    return CbcRangeSame; // fake return
}
