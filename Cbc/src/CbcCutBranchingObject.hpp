// Edwin 11/13/2009-- carved out of CbcBranchCut
#ifndef CbcCutBranchingObject_H
#define CbcCutBranchingObject_H

#include "CbcBranchBase.hpp"
#include "OsiRowCut.hpp"
#include "CoinPackedMatrix.hpp"

/** Cut branching object

  This object can specify a two-way branch in terms of two cuts
*/

class CbcCutBranchingObject : public CbcBranchingObject {

public:

    /// Default constructor
    CbcCutBranchingObject ();

    /** Create a cut branching object

        Cut down will applied on way=-1, up on way==1
        Assumed down will be first so way_ set to -1
    */
    CbcCutBranchingObject (CbcModel * model, OsiRowCut & down, OsiRowCut &up, bool canFix);

    /// Copy constructor
    CbcCutBranchingObject ( const CbcCutBranchingObject &);

    /// Assignment operator
    CbcCutBranchingObject & operator= (const CbcCutBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    /// Destructor
    virtual ~CbcCutBranchingObject ();

    using CbcBranchingObject::branch ;
    /** \brief Sets the bounds for variables or adds a cut depending on the
               current arm of the branch and advances the object state to the next arm.
           Returns change in guessed objective on next branch
    */
    virtual double branch();

#if 0
    // No need to override. Default works fine.
    /** Reset every information so that the branching object appears to point to
        the previous child. This method does not need to modify anything in any
        solver. */
    virtual void previousBranch();
#endif

    using CbcBranchingObject::print ;
    /** \brief Print something about branch - only if log level high
    */
    virtual void print();

    /** \brief Return true if branch should fix variables
    */
    virtual bool boundBranch() const;

    /** Return the type (an integer identifier) of \c this */
    virtual int type() const {
        return 200;
    }

    /** Compare the original object of \c this with the original object of \c
        brObj. Assumes that there is an ordering of the original objects.
        This method should be invoked only if \c this and brObj are of the same
        type.
        Return negative/0/positive depending on whether \c this is
        smaller/same/larger than the argument.
    */
    virtual int compareOriginalObject(const CbcBranchingObject* brObj) const;

    /** Compare the \c this with \c brObj. \c this and \c brObj must be os the
        same type and must have the same original object, but they may have
        different feasible regions.
        Return the appropriate CbcRangeCompare value (first argument being the
        sub/superset if that's the case). In case of overlap (and if \c
        replaceIfOverlap is true) replace the current branching object with one
        whose feasible region is the overlap.
     */
    virtual CbcRangeCompare compareBranchingObject
    (const CbcBranchingObject* brObj, const bool replaceIfOverlap = false);

protected:
    /// Cut for the down arm (way_ = -1)
    OsiRowCut down_;
    /// Cut for the up arm (way_ = 1)
    OsiRowCut up_;
    /// True if one way can fix variables
    bool canFix_;
};
#endif