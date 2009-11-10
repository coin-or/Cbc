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
#include "CbcFixingBranchingObject.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

//##############################################################################

// Default Constructor
CbcFixingBranchingObject::CbcFixingBranchingObject()
        : CbcBranchingObject()
{
    numberDown_ = 0;
    numberUp_ = 0;
    downList_ = NULL;
    upList_ = NULL;
}

// Useful constructor
CbcFixingBranchingObject::CbcFixingBranchingObject (CbcModel * model,
        int way ,
        int numberOnDownSide, const int * down,
        int numberOnUpSide, const int * up)
        : CbcBranchingObject(model, 0, way, 0.5)
{
    numberDown_ = numberOnDownSide;
    numberUp_ = numberOnUpSide;
    downList_ = CoinCopyOfArray(down, numberDown_);
    upList_ = CoinCopyOfArray(up, numberUp_);
}

// Copy constructor
CbcFixingBranchingObject::CbcFixingBranchingObject ( const CbcFixingBranchingObject & rhs) : CbcBranchingObject(rhs)
{
    numberDown_ = rhs.numberDown_;
    numberUp_ = rhs.numberUp_;
    downList_ = CoinCopyOfArray(rhs.downList_, numberDown_);
    upList_ = CoinCopyOfArray(rhs.upList_, numberUp_);
}

// Assignment operator
CbcFixingBranchingObject &
CbcFixingBranchingObject::operator=( const CbcFixingBranchingObject & rhs)
{
    if (this != &rhs) {
        CbcBranchingObject::operator=(rhs);
        delete [] downList_;
        delete [] upList_;
        numberDown_ = rhs.numberDown_;
        numberUp_ = rhs.numberUp_;
        downList_ = CoinCopyOfArray(rhs.downList_, numberDown_);
        upList_ = CoinCopyOfArray(rhs.upList_, numberUp_);
    }
    return *this;
}
CbcBranchingObject *
CbcFixingBranchingObject::clone() const
{
    return (new CbcFixingBranchingObject(*this));
}


// Destructor
CbcFixingBranchingObject::~CbcFixingBranchingObject ()
{
    delete [] downList_;
    delete [] upList_;
}
double
CbcFixingBranchingObject::branch()
{
    decrementNumberBranchesLeft();
    OsiSolverInterface * solver = model_->solver();
    const double * columnLower = solver->getColLower();
    int i;
    // *** for way - up means fix all those in up section
    if (way_ < 0) {
#ifdef FULL_PRINT
        printf("Down Fix ");
#endif
        //printf("Down Fix %d\n",numberDown_);
        for (i = 0; i < numberDown_; i++) {
            int iColumn = downList_[i];
            model_->solver()->setColUpper(iColumn, columnLower[iColumn]);
#ifdef FULL_PRINT
            printf("Setting bound on %d to lower bound\n", iColumn);
#endif
        }
        way_ = 1;	  // Swap direction
    } else {
#ifdef FULL_PRINT
        printf("Up Fix ");
#endif
        //printf("Up Fix %d\n",numberUp_);
        for (i = 0; i < numberUp_; i++) {
            int iColumn = upList_[i];
            model_->solver()->setColUpper(iColumn, columnLower[iColumn]);
#ifdef FULL_PRINT
            printf("Setting bound on %d to lower bound\n", iColumn);
#endif
        }
        way_ = -1;	  // Swap direction
    }
#ifdef FULL_PRINT
    printf("\n");
#endif
    return 0.0;
}
void
CbcFixingBranchingObject::print()
{
    int i;
    // *** for way - up means fix all those in up section
    if (way_ < 0) {
        printf("Down Fix ");
        for (i = 0; i < numberDown_; i++) {
            int iColumn = downList_[i];
            printf("%d ", iColumn);
        }
    } else {
        printf("Up Fix ");
        for (i = 0; i < numberUp_; i++) {
            int iColumn = upList_[i];
            printf("%d ", iColumn);
        }
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
CbcFixingBranchingObject::compareOriginalObject
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
CbcFixingBranchingObject::compareBranchingObject
(const CbcBranchingObject* /*brObj*/, const bool /*replaceIfOverlap*/)
{
#if 0 //ndef NDEBUG
    const CbcFixingBranchingObject* br =
        dynamic_cast<const CbcFixingBranchingObject*>(brObj);
    assert(br);
#endif
    // If two FixingBranchingObject's have the same base object then it's pretty
    // much guaranteed
    throw("must implement");
}

//##############################################################################
