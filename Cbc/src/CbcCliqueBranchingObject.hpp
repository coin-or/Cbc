// Edwin 11/10/2009-- carved out of CbcBranchActual
#ifndef CbcCbcCliqueBranchingObject_H
#define CbcCbcCliqueBranchingObject_H

#include "CbcBranchBase.hpp"
#include "CbcClique.hpp"

/** Branching object for unordered cliques

    Intended for cliques which are long enough to make it worthwhile
    but <= 64 members.  There will also be ones for long cliques.

    Variable_ is the clique id number (redundant, as the object also holds a
    pointer to the clique.
 */
class CbcCliqueBranchingObject : public CbcBranchingObject {

public:

    // Default Constructor
    CbcCliqueBranchingObject ();

    // Useful constructor
    CbcCliqueBranchingObject (CbcModel * model,  const CbcClique * clique,
                              int way,
                              int numberOnDownSide, const int * down,
                              int numberOnUpSide, const int * up);

    // Copy constructor
    CbcCliqueBranchingObject ( const CbcCliqueBranchingObject &);

    // Assignment operator
    CbcCliqueBranchingObject & operator=( const CbcCliqueBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    // Destructor
    virtual ~CbcCliqueBranchingObject ();

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
        return 102;
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
    const CbcClique * clique_;
    /// downMask - bit set to fix to weak bounds, not set to leave unfixed
    unsigned int downMask_[2];
    /// upMask - bit set to fix to weak bounds, not set to leave unfixed
    unsigned int upMask_[2];
};

#endif