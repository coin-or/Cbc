// Edwin 11/10/2009-- carved out of CbcBranchActual
#ifndef CbcSOSBranchingObject_H
#define CbcSOSBranchingObject_H

#include "CbcBranchBase.hpp"
#include "CbcSOS.hpp"

/** Branching object for Special ordered sets

    Variable_ is the set id number (redundant, as the object also holds a
    pointer to the set.
 */
class CbcSOSBranchingObject : public CbcBranchingObject {

public:

    // Default Constructor
    CbcSOSBranchingObject ();

    // Useful constructor
    CbcSOSBranchingObject (CbcModel * model,  const CbcSOS * clique,
                           int way,
                           double separator);

    // Copy constructor
    CbcSOSBranchingObject ( const CbcSOSBranchingObject &);

    // Assignment operator
    CbcSOSBranchingObject & operator=( const CbcSOSBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    // Destructor
    virtual ~CbcSOSBranchingObject ();

    using CbcBranchingObject::branch ;
    /// Does next branch and updates state
    virtual double branch();
    /** Update bounds in solver as in 'branch' and update given bounds.
        branchState is -1 for 'down' +1 for 'up' */
    virtual void fix(OsiSolverInterface * solver,
                     double * lower, double * upper,
                     int branchState) const ;

    /** Reset every information so that the branching object appears to point to
        the previous child. This method does not need to modify anything in any
        solver. */
    virtual void previousBranch() {
        CbcBranchingObject::previousBranch();
        computeNonzeroRange();
    }

    using CbcBranchingObject::print ;
    /** \brief Print something about branch - only if log level high
    */
    virtual void print();

    /** Return the type (an integer identifier) of \c this */
    virtual int type() const {
        return 104;
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

    /** Fill out the \c firstNonzero_ and \c lastNonzero_ data members */
    void computeNonzeroRange();

private:
    /// data
    const CbcSOS * set_;
    /// separator
    double separator_;
    /** The following two members describe the range in the members_ of the
        original object that whose upper bound is not fixed to 0. This is not
        necessary for Cbc to function correctly, this is there for heuristics so
        that separate branching decisions on the same object can be pooled into
        one branching object. */
    int firstNonzero_;
    int lastNonzero_;
};

#endif