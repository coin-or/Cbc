// Edwin 11/10/2009-- carved out of CbcBranchActual
#ifndef CbcGeneralBranchingObject_H
#define CbcGeneralBranchingObject_H

#include "CbcBranchBase.hpp"
#include "CbcSubProblem.hpp"

/** Branching object for general objects

 */
class CbcNode;
class CbcGeneralBranchingObject : public CbcBranchingObject {

public:

    // Default Constructor
    CbcGeneralBranchingObject ();

    // Useful constructor
    CbcGeneralBranchingObject (CbcModel * model);

    // Copy constructor
    CbcGeneralBranchingObject ( const CbcGeneralBranchingObject &);

    // Assignment operator
    CbcGeneralBranchingObject & operator=( const CbcGeneralBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    // Destructor
    virtual ~CbcGeneralBranchingObject ();

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
    /// Fill in current objective etc
    void state(double & objectiveValue, double & sumInfeasibilities,
               int & numberUnsatisfied, int which) const;
    /// Set CbcNode
    inline void setNode(CbcNode * node) {
        node_ = node;
    }
    /** Return the type (an integer identifier) of \c this */
    virtual int type() const {
        return 108;
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
    /// Number of subproblems
    inline int numberSubProblems() const {
        return numberSubProblems_;
    }
    /// Decrement number left and return number
    inline int decrementNumberLeft() {
        numberSubLeft_--;
        return numberSubLeft_;
    }
    /// Which node we want to use
    inline int whichNode() const {
        return whichNode_;
    }
    /// Set which node we want to use
    inline void setWhichNode(int value) {
        whichNode_ = value;
    }
    // Sub problem
    const CbcSubProblem * subProblem(int which) const {
        return subProblems_ + which;
    }

public:
    /// data
    // Sub problems
    CbcSubProblem * subProblems_;
    /// Node
    CbcNode * node_;
    /// Number of subproblems
    int numberSubProblems_;
    /// Number of subproblems left
    int numberSubLeft_;
    /// Which node we want to use (-1 for default)
    int whichNode_;
    /// Number of rows
    int numberRows_;
};

#endif