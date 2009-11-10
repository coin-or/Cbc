// Edwin 11/10/2009-- carved out of CbcBranchActual

#ifndef CbcOneGeneralBranchingObject_H
#define CbcOneGeneralBranchingObject_H

#include "CbcBranchBase.hpp"
/** Branching object for general objects - just one

 */
class CbcOneGeneralBranchingObject : public CbcBranchingObject {

public:

    // Default Constructor
    CbcOneGeneralBranchingObject ();

    // Useful constructor
    CbcOneGeneralBranchingObject (CbcModel * model,
                                  CbcGeneralBranchingObject * object,
                                  int whichOne);

    // Copy constructor
    CbcOneGeneralBranchingObject ( const CbcOneGeneralBranchingObject &);

    // Assignment operator
    CbcOneGeneralBranchingObject & operator=( const CbcOneGeneralBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    // Destructor
    virtual ~CbcOneGeneralBranchingObject ();

    using CbcBranchingObject::branch ;
    /// Does next branch and updates state
    virtual double branch();
    /** Double checks in case node can change its mind!
        Can change objective etc */
    virtual void checkIsCutoff(double cutoff);

    using CbcBranchingObject::print ;
    /** \brief Print something about branch - only if log level high
    */
    virtual void print();
    /** Return the type (an integer identifier) of \c this */
    virtual int type() const {
        return 110;
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

public:
    /// data
    /// Object
    CbcGeneralBranchingObject * object_;
    /// Which one
    int whichOne_;
};

#endif