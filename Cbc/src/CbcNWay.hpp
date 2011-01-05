// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/9/2009-- carved out of CbcBranchActual

/** Define an n-way class for variables.
    Only valid value is one at UB others at LB
    Normally 0-1
*/
#ifndef CbcNWay_H
#define CbcNWay_H

class CbcNWay : public CbcObject {

public:

    // Default Constructor
    CbcNWay ();

    /** Useful constructor (which are matrix indices)
    */
    CbcNWay (CbcModel * model, int numberMembers,
             const int * which, int identifier);

    // Copy constructor
    CbcNWay ( const CbcNWay &);

    /// Clone
    virtual CbcObject * clone() const;

    /// Assignment operator
    CbcNWay & operator=( const CbcNWay& rhs);

    /// Destructor
    virtual ~CbcNWay ();

    /// Set up a consequence for a single member
    void setConsequence(int iColumn, const CbcConsequence & consequence);

    /// Applies a consequence for a single member
    void applyConsequence(int iSequence, int state) const;

    /// Infeasibility - large is 0.5 (and 0.5 will give this)
    virtual double infeasibility(const OsiBranchingInformation * info,
                                 int &preferredWay) const;

    using CbcObject::feasibleRegion ;
    /// This looks at solution and sets bounds to contain solution
    virtual void feasibleRegion();

    /// Creates a branching object
    virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) ;

    /// Number of members
    inline int numberMembers() const {
        return numberMembers_;
    }

    /// Members (indices in range 0 ... numberColumns-1)
    inline const int * members() const {
        return members_;
    }
    /// Redoes data when sequence numbers change
    virtual void redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns);

protected:
    /// data
    /// Number of members
    int numberMembers_;

    /// Members (indices in range 0 ... numberColumns-1)
    int * members_;
    /// Consequences (normally NULL)
    CbcConsequence ** consequence_;
};
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

#ifdef JJF_ZERO
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
    virtual CbcBranchObjType type() const {
        return NWayBranchObj;
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
