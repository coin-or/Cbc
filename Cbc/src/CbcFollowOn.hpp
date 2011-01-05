// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#ifndef CbcFollowOn_H
#define CbcFollowOn_H

#include "CbcBranchBase.hpp"
#include "CoinPackedMatrix.hpp"

/** Define a follow on class.
    The idea of this is that in air-crew scheduling problems crew may fly in on flight A
    and out on flight B or on some other flight.  A useful branch is one which on one side
    fixes all which go out on flight B to 0, while the other branch fixes all those that do NOT
    go out on flight B to 0.

    This branching rule should be in addition to normal rules and have a high priority.
*/

class CbcFollowOn : public CbcObject {

public:

    // Default Constructor
    CbcFollowOn ();

    /** Useful constructor
    */
    CbcFollowOn (CbcModel * model);

    // Copy constructor
    CbcFollowOn ( const CbcFollowOn &);

    /// Clone
    virtual CbcObject * clone() const;

    // Assignment operator
    CbcFollowOn & operator=( const CbcFollowOn& rhs);

    // Destructor
    ~CbcFollowOn ();

    /// Infeasibility - large is 0.5
    virtual double infeasibility(const OsiBranchingInformation * info,
                                 int &preferredWay) const;

    using CbcObject::feasibleRegion ;
    /// This looks at solution and sets bounds to contain solution
    virtual void feasibleRegion();

    /// Creates a branching object
    virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) ;
    /// As some computation is needed in more than one place - returns row
    virtual int gutsOfFollowOn(int & otherRow, int & preferredWay) const;

protected:
    /// data
    /// Matrix
    CoinPackedMatrix matrix_;
    /// Matrix by row
    CoinPackedMatrix matrixByRow_;
    /// Possible rhs (if 0 then not possible)
    int * rhs_;
};

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

#ifdef JJF_ZERO
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
    virtual CbcBranchObjType type() const {
        return FollowOnBranchObj;
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

