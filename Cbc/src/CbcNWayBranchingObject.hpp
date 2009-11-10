// Edwin 11/10/2009-- carved out of CbcBranchActual
#ifndef CbcNWayBranchingObject_H
#define CbcNWayBranchingObject_H

#include "CbcBranchBase.hpp"
#include "CbcNWay.hpp"
/** N way branching Object class.
    Variable is number of set.
 */
class CbcNWayBranchingObject : public CbcBranchingObject {

public:

    // Default Constructor
    CbcNWayBranchingObject ();

    /** Useful constructor - order had matrix indices
        way_ -1 corresponds to setting first, +1 to second, +3 etc.
        this is so -1 and +1 have similarity to normal
    */
    CbcNWayBranchingObject (CbcModel * model,  const CbcNWay * nway,
                            int numberBranches, const int * order);

    // Copy constructor
    CbcNWayBranchingObject ( const CbcNWayBranchingObject &);

    // Assignment operator
    CbcNWayBranchingObject & operator=( const CbcNWayBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    // Destructor
    virtual ~CbcNWayBranchingObject ();

    using CbcBranchingObject::branch ;
    /// Does next branch and updates state
    virtual double branch();

#if 0
    // FIXME: what do we need to do here?
    /** Reset every information so that the branching object appears to point to
        the previous child. This method does not need to modify anything in any
        solver. */
    virtual void previousBranch();
#endif

    using CbcBranchingObject::print ;
    /** \brief Print something about branch - only if log level high
    */
    virtual void print();
    /** The number of branch arms created for this branching object
    */
    virtual int numberBranches() const {
        return numberInSet_;
    }
    /// Is this a two way object (-1 down, +1 up)
    virtual bool twoWay() const {
        return false;
    }

    /** Return the type (an integer identifier) of \c this */
    virtual int type() const {
        return 105;
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
    /// order of branching - points back to CbcNWay
    int * order_;
    /// Points back to object
    const CbcNWay * object_;
    /// Number in set
    int numberInSet_;
};

#endif
