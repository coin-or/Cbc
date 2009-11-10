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
#include "CbcCliqueBranchingObject.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

// Default Constructor
CbcCliqueBranchingObject::CbcCliqueBranchingObject()
        : CbcBranchingObject()
{
    clique_ = NULL;
    downMask_[0] = 0;
    downMask_[1] = 0;
    upMask_[0] = 0;
    upMask_[1] = 0;
}

// Useful constructor
CbcCliqueBranchingObject::CbcCliqueBranchingObject (CbcModel * model,
        const CbcClique * clique,
        int way ,
        int numberOnDownSide, const int * down,
        int numberOnUpSide, const int * up)
        : CbcBranchingObject(model, clique->id(), way, 0.5)
{
    clique_ = clique;
    downMask_[0] = 0;
    downMask_[1] = 0;
    upMask_[0] = 0;
    upMask_[1] = 0;
    int i;
    for (i = 0; i < numberOnDownSide; i++) {
        int sequence = down[i];
        int iWord = sequence >> 5;
        int iBit = sequence - 32 * iWord;
        unsigned int k = 1 << iBit;
        downMask_[iWord] |= k;
    }
    for (i = 0; i < numberOnUpSide; i++) {
        int sequence = up[i];
        int iWord = sequence >> 5;
        int iBit = sequence - 32 * iWord;
        unsigned int k = 1 << iBit;
        upMask_[iWord] |= k;
    }
}

// Copy constructor
CbcCliqueBranchingObject::CbcCliqueBranchingObject ( const CbcCliqueBranchingObject & rhs) : CbcBranchingObject(rhs)
{
    clique_ = rhs.clique_;
    downMask_[0] = rhs.downMask_[0];
    downMask_[1] = rhs.downMask_[1];
    upMask_[0] = rhs.upMask_[0];
    upMask_[1] = rhs.upMask_[1];
}

// Assignment operator
CbcCliqueBranchingObject &
CbcCliqueBranchingObject::operator=( const CbcCliqueBranchingObject & rhs)
{
    if (this != &rhs) {
        CbcBranchingObject::operator=(rhs);
        clique_ = rhs.clique_;
        downMask_[0] = rhs.downMask_[0];
        downMask_[1] = rhs.downMask_[1];
        upMask_[0] = rhs.upMask_[0];
        upMask_[1] = rhs.upMask_[1];
    }
    return *this;
}
CbcBranchingObject *
CbcCliqueBranchingObject::clone() const
{
    return (new CbcCliqueBranchingObject(*this));
}


// Destructor
CbcCliqueBranchingObject::~CbcCliqueBranchingObject ()
{
}
double
CbcCliqueBranchingObject::branch()
{
    decrementNumberBranchesLeft();
    int iWord;
    int numberMembers = clique_->numberMembers();
    const int * which = clique_->members();
    const int * integerVariables = model_->integerVariable();
    int numberWords = (numberMembers + 31) >> 5;
    // *** for way - up means fix all those in down section
    if (way_ < 0) {
#ifdef FULL_PRINT
        printf("Down Fix ");
#endif
        for (iWord = 0; iWord < numberWords; iWord++) {
            int i;
            for (i = 0; i < 32; i++) {
                unsigned int k = 1 << i;
                if ((upMask_[iWord]&k) != 0) {
                    int iColumn = which[i+32*iWord];
#ifdef FULL_PRINT
                    printf("%d ", i + 32*iWord);
#endif
                    // fix weak way
                    if (clique_->type(i + 32*iWord))
                        model_->solver()->setColUpper(integerVariables[iColumn], 0.0);
                    else
                        model_->solver()->setColLower(integerVariables[iColumn], 1.0);
                }
            }
        }
        way_ = 1;	  // Swap direction
    } else {
#ifdef FULL_PRINT
        printf("Up Fix ");
#endif
        for (iWord = 0; iWord < numberWords; iWord++) {
            int i;
            for (i = 0; i < 32; i++) {
                unsigned int k = 1 << i;
                if ((downMask_[iWord]&k) != 0) {
                    int iColumn = which[i+32*iWord];
#ifdef FULL_PRINT
                    printf("%d ", i + 32*iWord);
#endif
                    // fix weak way
                    if (clique_->type(i + 32*iWord))
                        model_->solver()->setColUpper(integerVariables[iColumn], 0.0);
                    else
                        model_->solver()->setColLower(integerVariables[iColumn], 1.0);
                }
            }
        }
        way_ = -1;	  // Swap direction
    }
#ifdef FULL_PRINT
    printf("\n");
#endif
    return 0.0;
}
// Print what would happen
void
CbcCliqueBranchingObject::print()
{
    int iWord;
    int numberMembers = clique_->numberMembers();
    const int * which = clique_->members();
    const int * integerVariables = model_->integerVariable();
    int numberWords = (numberMembers + 31) >> 5;
    // *** for way - up means fix all those in down section
    if (way_ < 0) {
        printf("Clique - Down Fix ");
        for (iWord = 0; iWord < numberWords; iWord++) {
            int i;
            for (i = 0; i < 32; i++) {
                unsigned int k = 1 << i;
                if ((upMask_[iWord]&k) != 0) {
                    int iColumn = which[i+32*iWord];
                    printf("%d ", integerVariables[iColumn]);
                }
            }
        }
    } else {
        printf("Clique - Up Fix ");
        for (iWord = 0; iWord < numberWords; iWord++) {
            int i;
            for (i = 0; i < 32; i++) {
                unsigned int k = 1 << i;
                if ((downMask_[iWord]&k) != 0) {
                    int iColumn = which[i+32*iWord];
                    printf("%d ", integerVariables[iColumn]);
                }
            }
        }
    }
    printf("\n");
}

static inline int
CbcCompareCliques(const CbcClique* cl0, const CbcClique* cl1)
{
    if (cl0->cliqueType() < cl1->cliqueType()) {
        return -1;
    }
    if (cl0->cliqueType() > cl1->cliqueType()) {
        return 1;
    }
    if (cl0->numberMembers() != cl1->numberMembers()) {
        return cl0->numberMembers() - cl1->numberMembers();
    }
    if (cl0->numberNonSOSMembers() != cl1->numberNonSOSMembers()) {
        return cl0->numberNonSOSMembers() - cl1->numberNonSOSMembers();
    }
    return memcmp(cl0->members(), cl1->members(),
                  cl0->numberMembers() * sizeof(int));
}

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcCliqueBranchingObject::compareOriginalObject
(const CbcBranchingObject* brObj) const
{
    const CbcCliqueBranchingObject* br =
        dynamic_cast<const CbcCliqueBranchingObject*>(brObj);
    assert(br);
    return CbcCompareCliques(clique_, br->clique_);
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
CbcCliqueBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool /*replaceIfOverlap*/)
{
    const CbcCliqueBranchingObject* br =
        dynamic_cast<const CbcCliqueBranchingObject*>(brObj);
    assert(br);
    unsigned int* thisMask = way_ < 0 ? upMask_ : downMask_;
    const unsigned int* otherMask = br->way_ < 0 ? br->upMask_ : br->downMask_;
    const CoinUInt64 cl0 =
        (static_cast<CoinUInt64>(thisMask[0]) << 32) | thisMask[1];
    const CoinUInt64 cl1 =
        (static_cast<CoinUInt64>(otherMask[0]) << 32) | otherMask[1];
    if (cl0 == cl1) {
        return CbcRangeSame;
    }
    const CoinUInt64 cl_intersection = (cl0 & cl1);
    if (cl_intersection == cl0) {
        return CbcRangeSuperset;
    }
    if (cl_intersection == cl1) {
        return CbcRangeSubset;
    }
    const CoinUInt64 cl_xor = (cl0 ^ cl1);
    if (cl_intersection == 0 && cl_xor == 0) {
        return CbcRangeDisjoint;
    }
    const CoinUInt64 cl_union = (cl0 | cl1);
    thisMask[0] = static_cast<unsigned int>(cl_union >> 32);
    thisMask[1] = static_cast<unsigned int>(cl_union & 0xffffffff);
    return CbcRangeOverlap;
}
