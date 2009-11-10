// Edwin 11/10/2009-- carved out of CbcBranchActual
#ifndef CbcFixingBranchingObject_H
#define CbcFixingBranchingObject_H

#include "CbcBranchBase.hpp"
/** General Branching Object class.
    Each way fixes some variables to lower bound
 */
class CbcFixingBranchingObject : public CbcBranchingObject {

public:

    // Default Constructor
    CbcFixingBranchingObject ();

    // Useful constructor
    CbcFixingBranchingObject (CbcModel * model,
                              int way,
                              int numberOnDownSide, const int * down,
                              int numberOnUpSide, const int * up);

    // Copy constructor
    CbcFixingBranchingObject ( const CbcFixingBranchingObject &);

    // Assignment operator
    CbcFixingBranchingObject & operator=( const CbcFixingBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    // Destructor
    virtual ~CbcFixingBranchingObject ();

    using CbcBranchingObject::branch ;
    /// Does next branch and updates state
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

    /** Return the type (an integer identifier) of \c this */
    virtual int type() const {
        return 106;
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

private:
    /// data
    /// Number on down list
    int numberDown_;
    /// Number on up list
    int numberUp_;
    /// downList - variables to fix to lb on down branch
    int * downList_;
    /// upList - variables to fix to lb on up branch
    int * upList_;
};

#endif