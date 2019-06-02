// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/12/2009 carved from CbcBranchBase

#ifndef CbcBranchingObject_H
#define CbcBranchingObject_H

#include <string>
#include <vector>
#include "CbcBranchBase.hpp"
#include "OsiBranchingObject.hpp"

// The types of objects that will be derived from this class.
enum CbcBranchObjType {
  SimpleIntegerBranchObj = 100,
  SimpleIntegerDynamicPseudoCostBranchObj = 101,
  CliqueBranchObj = 102,
  LongCliqueBranchObj = 103,
  SoSBranchObj = 104,
  NWayBranchObj = 105,
  FollowOnBranchObj = 106,
  DummyBranchObj = 107,
  GeneralDepthBranchObj = 108,
  OneGeneralBranchingObj = 110,
  CutBranchingObj = 200,
  LotsizeBranchObj = 300,
  DynamicPseudoCostBranchObj = 400
};

/** \brief Abstract branching object base class
    Now just difference with OsiBranchingObject

  In the abstract, an CbcBranchingObject contains instructions for how to
  branch. We want an abstract class so that we can describe how to branch on
  simple objects (<i>e.g.</i>, integers) and more exotic objects
  (<i>e.g.</i>, cliques or hyperplanes).

  The #branch() method is the crucial routine: it is expected to be able to
  step through a set of branch arms, executing the actions required to create
  each subproblem in turn. The base class is primarily virtual to allow for
  a wide range of problem modifications.

  See CbcObject for an overview of the three classes (CbcObject,
  CbcBranchingObject, and CbcBranchDecision) which make up cbc's branching
  model.
*/

class CbcBranchingObject : public OsiBranchingObject {

public:
  /// Default Constructor
  CbcBranchingObject();

  /// Constructor
  CbcBranchingObject(CbcModel *model, int variable, int way, double value);

  /// Copy constructor
  CbcBranchingObject(const CbcBranchingObject &);

  /// Assignment operator
  CbcBranchingObject &operator=(const CbcBranchingObject &rhs);

  /// Clone
  virtual CbcBranchingObject *clone() const = 0;

  /// Destructor
  virtual ~CbcBranchingObject();

  /** Some branchingObjects may claim to be able to skip
        strong branching.  If so they have to fill in CbcStrongInfo.
        The object mention in incoming CbcStrongInfo must match.
        Returns nonzero if skip is wanted */
  virtual int fillStrongInfo(CbcStrongInfo &)
  {
    return 0;
  }
  /// Reset number of branches left to original
  inline void resetNumberBranchesLeft()
  {
    branchIndex_ = 0;
  }
  /// Set number of branches to do
  inline void setNumberBranches(int value)
  {
    branchIndex_ = 0;
    numberBranches_ = value;
  }

  /** \brief Execute the actions required to branch, as specified by the
           current state of the branching object, and advance the object's
           state.  Mainly for diagnostics, whether it is true branch or
           strong branching is also passed.
           Returns change in guessed objective on next branch
    */
  virtual double branch() = 0;
  /** \brief Execute the actions required to branch, as specified by the
           current state of the branching object, and advance the object's
           state.  Mainly for diagnostics, whether it is true branch or
           strong branching is also passed.
           Returns change in guessed objective on next branch
    */
  virtual double branch(OsiSolverInterface *)
  {
    return branch();
  }
  /** Update bounds in solver as in 'branch' and update given bounds.
        branchState is -1 for 'down' +1 for 'up' */
  virtual void fix(OsiSolverInterface *,
    double *, double *,
    int) const {}

  /** Change (tighten) bounds in object to reflect bounds in solver.
	Return true if now fixed */
  virtual bool tighten(OsiSolverInterface *) { return false; }

  /** Reset every information so that the branching object appears to point to
        the previous child. This method does not need to modify anything in any
        solver. */
  virtual void previousBranch()
  {
    assert(branchIndex_ > 0);
    branchIndex_--;
    way_ = -way_;
  }

  using OsiBranchingObject::print;
  /** \brief Print something about branch - only if log level high
    */
  virtual void print() const {}

  /** \brief Index identifying the associated CbcObject within its class.

      The name is misleading, and typically the index will <i>not</i> refer
      directly to a variable.
      Rather, it identifies an CbcObject within the class of similar
      CbcObjects

      <i>E.g.</i>, for an CbcSimpleInteger, variable() is the index of the
      integer variable in the set of integer variables (<i>not</i> the index of
      the variable in the set of all variables).
    */
  inline int variable() const
  {
    return variable_;
  }

  /** Get the state of the branching object

      Returns a code indicating the active arm of the branching object.
      The precise meaning is defined in the derived class.

      \sa #way_
    */
  inline int way() const
  {
    return way_;
  }

  /** Set the state of the branching object.

      See #way()
    */
  inline void way(int way)
  {
    way_ = way;
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

  /// Return pointer back to object which created
  inline CbcObject *object() const
  {
    return originalCbcObject_;
  }
  /// Set pointer back to object which created
  inline void setOriginalObject(CbcObject *object)
  {
    originalCbcObject_ = object;
  }

  // Methods used in heuristics

  /** Return the type (an integer identifier) of \c this.
        See definition of CbcBranchObjType above for possibilities
    */

  virtual CbcBranchObjType type() const = 0;

  /** Compare the original object of \c this with the original object of \c
        brObj. Assumes that there is an ordering of the original objects.
        This method should be invoked only if \c this and brObj are of the same
        type.
        Return negative/0/positive depending on whether \c this is
        smaller/same/larger than the argument.
    */
  virtual int compareOriginalObject(const CbcBranchingObject *brObj) const
  {
    const CbcBranchingObject *br = dynamic_cast< const CbcBranchingObject * >(brObj);
    return variable() - br->variable();
  }

  /** Compare the \c this with \c brObj. \c this and \c brObj must be of the
        same type and must have the same original object, but they may have
        different feasible regions.
        Return the appropriate CbcRangeCompare value (first argument being the
        sub/superset if that's the case). In case of overlap (and if \c
        replaceIfOverlap is true) replace the current branching object with one
        whose feasible region is the overlap.
     */
  virtual CbcRangeCompare compareBranchingObject(const CbcBranchingObject *brObj, const bool replaceIfOverlap = false) = 0;

protected:
  /// The model that owns this branching object
  CbcModel *model_;
  /// Pointer back to object which created
  CbcObject *originalCbcObject_;

  /// Branching variable (0 is first integer)
  int variable_;
  // was - Way to branch - -1 down (first), 1 up, -2 down (second), 2 up (second)
  /** The state of the branching object.

      Specifies the active arm of the branching object. Coded as -1 to take
      the `down' arm, +1 for the `up' arm. `Down' and `up' are defined based on
      the natural meaning (floor and ceiling, respectively) for a simple integer.
      The precise meaning is defined in the derived class.
    */
  int way_;
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
