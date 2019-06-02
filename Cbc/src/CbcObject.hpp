// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/12/2009 carved from CbcBranchBase

#ifndef CbcObject_H
#define CbcObject_H

#include <string>
#include <vector>
#include "OsiBranchingObject.hpp"
class OsiSolverInterface;
class OsiSolverBranch;

class CbcModel;
class CbcNode;
class CbcNodeInfo;
class CbcBranchingObject;
class OsiChooseVariable;
class CbcObjectUpdateData;
//#############################################################################

/** Abstract base class for `objects'.
    It now just has stuff that OsiObject does not have

  The branching model used in Cbc is based on the idea of an <i>object</i>.
  In the abstract, an object is something that has a feasible region, can be
  evaluated for infeasibility, can be branched on (<i>i.e.</i>, there's some
  constructive action to be taken to move toward feasibility), and allows
  comparison of the effect of branching.

  This class (CbcObject) is the base class for an object. To round out the
  branching model, the class CbcBranchingObject describes how to perform a
  branch, and the class CbcBranchDecision describes how to compare two
  CbcBranchingObjects.

  To create a new type of object you need to provide three methods:
  #infeasibility(), #feasibleRegion(), and #createCbcBranch(), described below.

  This base class is primarily virtual to allow for any form of structure.
  Any form of discontinuity is allowed.

  \todo The notion that all branches are binary (two arms) is wired into the
	implementation of CbcObject, CbcBranchingObject, and
	CbcBranchDecision. Changing this will require a moderate amount of
	recoding.
 */
// This can be used if object wants to skip strong branching
typedef struct {
  CbcBranchingObject *possibleBranch; // what a branch would do
  double upMovement; // cost going up (and initial away from feasible)
  double downMovement; // cost going down
  int numIntInfeasUp; // without odd ones
  int numObjInfeasUp; // just odd ones
  bool finishedUp; // true if solver finished
  int numItersUp; // number of iterations in solver
  int numIntInfeasDown; // without odd ones
  int numObjInfeasDown; // just odd ones
  bool finishedDown; // true if solver finished
  int numItersDown; // number of iterations in solver
  int objectNumber; // Which object it is
  int fix; // 0 if no fix, 1 if we can fix up, -1 if we can fix down
} CbcStrongInfo;

class CbcObject : public OsiObject {

public:
  // Default Constructor
  CbcObject();

  // Useful constructor
  CbcObject(CbcModel *model);

  // Copy constructor
  CbcObject(const CbcObject &);

  // Assignment operator
  CbcObject &operator=(const CbcObject &rhs);

  /// Clone
  virtual CbcObject *clone() const = 0;

  /// Destructor
  virtual ~CbcObject();

  /** Infeasibility of the object

        This is some measure of the infeasibility of the object. It should be
        scaled to be in the range [0.0, 0.5], with 0.0 indicating the object
        is satisfied.

        The preferred branching direction is returned in preferredWay,

        This is used to prepare for strong branching but should also think of
        case when no strong branching

        The object may also compute an estimate of cost of going "up" or "down".
        This will probably be based on pseudo-cost ideas
    */
#ifdef CBC_NEW_STYLE_BRANCH
  virtual double infeasibility(const OsiBranchingInformation *info,
    int &preferredWay) const = 0;
#else
  virtual double infeasibility(const OsiBranchingInformation * /*info*/,
    int &preferredWay) const
  {
    return infeasibility(preferredWay);
  }
  virtual double infeasibility(int & /*preferredWay*/) const
  {
    throw CoinError("Need code", "infeasibility", "CbcBranchBase");
  }
#endif

  /** For the variable(s) referenced by the object,
        look at the current solution and set bounds to match the solution.
    */
  virtual void feasibleRegion() = 0;
  /// Dummy one for compatibility
  virtual double feasibleRegion(OsiSolverInterface *solver, const OsiBranchingInformation *info) const;

  /** For the variable(s) referenced by the object,
        look at the current solution and set bounds to match the solution.
        Returns measure of how much it had to move solution to make feasible
    */
  virtual double feasibleRegion(OsiSolverInterface *solver) const;

  /** Create a branching object and indicate which way to branch first.

        The branching object has to know how to create branches (fix
        variables, etc.)
    */
#ifdef CBC_NEW_STYLE_BRANCH
  virtual CbcBranchingObject *createCbcBranch(OsiSolverInterface *solver, const OsiBranchingInformation *info, int way) = 0;
#else
  virtual CbcBranchingObject *createCbcBranch(OsiSolverInterface *
    /* solver */,
    const OsiBranchingInformation *
    /* info */,
    int /* way */)
  {
    // return createBranch(solver, info, way);
    return NULL;
  }
  virtual OsiBranchingObject *createBranch(OsiSolverInterface * /*solver*/,
    const OsiBranchingInformation * /*info*/, int /*way*/) const
  {
    throw CoinError("Need code", "createBranch", "CbcBranchBase");
  }
#endif
  /** Create an Osibranching object and indicate which way to branch first.

        The branching object has to know how to create branches (fix
        variables, etc.)
    */
  virtual OsiBranchingObject *createOsiBranch(OsiSolverInterface *solver, const OsiBranchingInformation *info, int way) const;
  /** Create an OsiSolverBranch object

        This returns NULL if branch not represented by bound changes
    */
  virtual OsiSolverBranch *solverBranch() const;

  /** \brief Given a valid solution (with reduced costs, etc.),
        return a branching object which would give a new feasible
        point in a good direction.

        If the method cannot generate a feasible point (because there aren't
        any, or because it isn't bright enough to find one), it should
        return null.
    */
  virtual CbcBranchingObject *preferredNewFeasible() const
  {
    return NULL;
  }

  /** \brief Given a valid solution (with reduced costs, etc.),
        return a branching object which would give a new feasible
        point in a bad direction.

        If the method cannot generate a feasible point (because there aren't
        any, or because it isn't bright enough to find one), it should
        return null.
    */
  virtual CbcBranchingObject *notPreferredNewFeasible() const
  {
    return NULL;
  }

  /** Reset variable bounds to their original values.

      Bounds may be tightened, so it may be good to be able to set this info in object.
     */
  virtual void resetBounds(const OsiSolverInterface *) {}

  /** Returns floor and ceiling i.e. closest valid points
    */
  virtual void floorCeiling(double &floorValue, double &ceilingValue, double value,
    double tolerance) const;

  /** Pass in information on branch just done and create CbcObjectUpdateData instance.
        If object does not need data then backward pointer will be NULL.
        Assumes can get information from solver */
  virtual CbcObjectUpdateData createUpdateInformation(const OsiSolverInterface *solver,
    const CbcNode *node,
    const CbcBranchingObject *branchingObject);

  /// Update object by CbcObjectUpdateData
  virtual void updateInformation(const CbcObjectUpdateData &) {}

  /// Identifier (normally column number in matrix)
  inline int id() const
  {
    return id_;
  }

  /** Set identifier (normally column number in matrix)
        but 1000000000 to 1100000000 means optional branching object
        i.e. code would work without it */
  inline void setId(int value)
  {
    id_ = value;
  }

  /** Return true if optional branching object
        i.e. code would work without it */
  inline bool optionalObject() const
  {
    return (id_ >= 1000000000 && id_ < 1100000000);
  }

  /// Get position in object_ list
  inline int position() const
  {
    return position_;
  }

  /// Set position in object_ list
  inline void setPosition(int position)
  {
    position_ = position;
  }

  /// update model
  inline void setModel(CbcModel *model)
  {
    model_ = model;
  }

  /// Return model
  inline CbcModel *model() const
  {
    return model_;
  }

  /// If -1 down always chosen first, +1 up always, 0 normal
  inline int preferredWay() const
  {
    return preferredWay_;
  }
  /// Set -1 down always chosen first, +1 up always, 0 normal
  inline void setPreferredWay(int value)
  {
    preferredWay_ = value;
  }
  /// Redoes data when sequence numbers change
  virtual void redoSequenceEtc(CbcModel *, int, const int *) {}
  /// Initialize for branching
  virtual void initializeForBranching(CbcModel *) {}

protected:
  /// data

  /// Model
  CbcModel *model_;
  /// Identifier (normally column number in matrix)
  int id_;
  /// Position in object list
  int position_;
  /// If -1 down always chosen first, +1 up always, 0 normal
  int preferredWay_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
