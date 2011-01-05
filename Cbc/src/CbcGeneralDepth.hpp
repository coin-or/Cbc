// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#ifndef CbcGeneralDepth_H
#define CbcGeneralDepth_H

#include "CbcGeneral.hpp"
#include "CbcBranchBase.hpp"
#include "CbcSubProblem.hpp"

#ifdef COIN_HAS_CLP

/** Define a catch all class.
    This will create a list of subproblems using partial evaluation
*/
#include "ClpSimplex.hpp"
#include "ClpNode.hpp"


class CbcGeneralDepth : public CbcGeneral {

public:

    // Default Constructor
    CbcGeneralDepth ();

    /** Useful constructor
        Just needs to point to model.
        Initial version does evaluation to depth N
        This is stored in CbcModel but may be
        better here
    */
    CbcGeneralDepth (CbcModel * model, int maximumDepth);

    // Copy constructor
    CbcGeneralDepth ( const CbcGeneralDepth &);

    /// Clone
    virtual CbcObject * clone() const;

    // Assignment operator
    CbcGeneralDepth & operator=( const CbcGeneralDepth& rhs);

    // Destructor
    ~CbcGeneralDepth ();

    /// Infeasibility - large is 0.5
    virtual double infeasibility(const OsiBranchingInformation * info,
                                 int &preferredWay) const;

    using CbcObject::feasibleRegion ;
    /// This looks at solution and sets bounds to contain solution
    virtual void feasibleRegion();

    /// Creates a branching object
    virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) ;
    /// Return maximum number of nodes
    inline int maximumNodes() const {
        return maximumNodes_;
    }
    /// Get maximum depth
    inline int maximumDepth() const {
        return maximumDepth_;
    }
    /// Set maximum depth
    inline void setMaximumDepth(int value) {
        maximumDepth_ = value;
    }
    /// Get which solution
    inline int whichSolution() const {
        return whichSolution_;
    }
    /// Get ClpNode info
    inline ClpNode * nodeInfo(int which) {
        return nodeInfo_->nodeInfo_[which];
    }

    /// Redoes data when sequence numbers change
    virtual void redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns);

protected:
    /// data
    /// Maximum depth
    int maximumDepth_;
    /// Maximum nodes
    int maximumNodes_;
    /// Which node has solution (or -1)
    mutable int whichSolution_;
    /// Number of valid nodes (including whichSolution_)
    mutable int numberNodes_;
    /// For solving nodes
    mutable ClpNodeStuff * nodeInfo_;
};
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
    virtual CbcBranchObjType type() const {
        return GeneralDepthBranchObj;
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
    virtual CbcBranchObjType type() const {
        return OneGeneralBranchingObj;
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
#endif //COIN_HAS_CLP
#endif

