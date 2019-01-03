/* $Id$ */
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcBranchDynamic_H
#define CbcBranchDynamic_H

#include "CoinPackedMatrix.hpp"
#include "CbcSimpleIntegerDynamicPseudoCost.hpp"
#include "CbcBranchActual.hpp"

/** Branching decision dynamic class

  This class implements a simple algorithm
  (betterBranch()) for choosing a branching variable when dynamic pseudo costs.
*/

class CbcBranchDynamicDecision : public CbcBranchDecision {
public:
  // Default Constructor
  CbcBranchDynamicDecision();

  // Copy constructor
  CbcBranchDynamicDecision(const CbcBranchDynamicDecision &);

  virtual ~CbcBranchDynamicDecision();

  /// Clone
  virtual CbcBranchDecision *clone() const;

  /// Initialize, <i>e.g.</i> before the start of branch selection at a node
  virtual void initialize(CbcModel *model);

  /** \brief Compare two branching objects. Return nonzero if \p thisOne is
           better than \p bestSoFar.

      The routine compares branches using the values supplied in \p numInfUp and
      \p numInfDn until a solution is found by search, after which it uses the
      values supplied in \p changeUp and \p changeDn. The best branching object
      seen so far and the associated parameter values are remembered in the
      \c CbcBranchDynamicDecision object. The nonzero return value is +1 if the
      up branch is preferred, -1 if the down branch is preferred.

      As the names imply, the assumption is that the values supplied for
      \p numInfUp and \p numInfDn will be the number of infeasibilities reported
      by the branching object, and \p changeUp and \p changeDn will be the
      estimated change in objective. Other measures can be used if desired.

      Because an \c CbcBranchDynamicDecision object remembers the current best
      branching candidate (#bestObject_) as well as the values used in the
      comparison, the parameter \p bestSoFar is redundant, hence unused.
    */
  virtual int betterBranch(CbcBranchingObject *thisOne,
    CbcBranchingObject *bestSoFar,
    double changeUp, int numInfUp,
    double changeDn, int numInfDn);
  /** Sets or gets best criterion so far */
  virtual void setBestCriterion(double value);
  virtual double getBestCriterion() const;
  /** Says whether this method can handle both methods -
        1 better, 2 best, 3 both */
  virtual int whichMethod()
  {
    return 3;
  }

  /** Saves a clone of current branching object.  Can be used to update
        information on object causing branch - after branch */
  virtual void saveBranchingObject(OsiBranchingObject *object);
  /** Pass in information on branch just done.
        assumes object can get information from solver */
  virtual void updateInformation(OsiSolverInterface *solver,
    const CbcNode *node);

private:
  /// Illegal Assignment operator
  CbcBranchDynamicDecision &operator=(const CbcBranchDynamicDecision &rhs);

  /// data

  /// "best" so far
  double bestCriterion_;

  /// Change up for best
  double bestChangeUp_;

  /// Number of infeasibilities for up
  int bestNumberUp_;

  /// Change down for best
  double bestChangeDown_;

  /// Number of infeasibilities for down
  int bestNumberDown_;

  /// Pointer to best branching object
  CbcBranchingObject *bestObject_;
};
/** Simple branching object for an integer variable with pseudo costs

  This object can specify a two-way branch on an integer variable. For each
  arm of the branch, the upper and lower bounds on the variable can be
  independently specified.

  Variable_ holds the index of the integer variable in the integerVariable_
  array of the model.
*/

class CbcDynamicPseudoCostBranchingObject : public CbcIntegerBranchingObject {

public:
  /// Default constructor
  CbcDynamicPseudoCostBranchingObject();

  /** Create a standard floor/ceiling branch object

      Specifies a simple two-way branch. Let \p value = x*. One arm of the
      branch will be is lb <= x <= floor(x*), the other ceil(x*) <= x <= ub.
      Specify way = -1 to set the object state to perform the down arm first,
      way = 1 for the up arm.
    */
  CbcDynamicPseudoCostBranchingObject(CbcModel *model, int variable,
    int way, double value,
    CbcSimpleIntegerDynamicPseudoCost *object);

  /** Create a degenerate branch object

      Specifies a `one-way branch'. Calling branch() for this object will
      always result in lowerValue <= x <= upperValue. Used to fix a variable
      when lowerValue = upperValue.
    */

  CbcDynamicPseudoCostBranchingObject(CbcModel *model, int variable, int way,
    double lowerValue, double upperValue);

  /// Copy constructor
  CbcDynamicPseudoCostBranchingObject(const CbcDynamicPseudoCostBranchingObject &);

  /// Assignment operator
  CbcDynamicPseudoCostBranchingObject &operator=(const CbcDynamicPseudoCostBranchingObject &rhs);

  /// Clone
  virtual CbcBranchingObject *clone() const;

  /// Destructor
  virtual ~CbcDynamicPseudoCostBranchingObject();

  /// Does part of constructor
  void fillPart(int variable,
    int way, double value,
    CbcSimpleIntegerDynamicPseudoCost *object);

  using CbcBranchingObject::branch;
  /** \brief Sets the bounds for the variable according to the current arm
           of the branch and advances the object state to the next arm.
           This version also changes guessed objective value
    */
  virtual double branch();

  /** Some branchingObjects may claim to be able to skip
        strong branching.  If so they have to fill in CbcStrongInfo.
        The object mention in incoming CbcStrongInfo must match.
        Returns nonzero if skip is wanted */
  virtual int fillStrongInfo(CbcStrongInfo &info);

  /// Change in guessed
  inline double changeInGuessed() const
  {
    return changeInGuessed_;
  }
  /// Set change in guessed
  inline void setChangeInGuessed(double value)
  {
    changeInGuessed_ = value;
  }
  /// Return object
  inline CbcSimpleIntegerDynamicPseudoCost *object() const
  {
    return object_;
  }
  /// Set object
  inline void setObject(CbcSimpleIntegerDynamicPseudoCost *object)
  {
    object_ = object;
  }

  /** Return the type (an integer identifier) of \c this */
  virtual CbcBranchObjType type() const
  {
    return DynamicPseudoCostBranchObj;
  }

  // LL: compareOriginalObject and compareBranchingObject are inherited from
  // CbcIntegerBranchingObject thus need not be declared/defined here. After
  // all, this kind of branching object is simply using pseudocosts to make
  // decisions, but once the decisions are made they are the same kind as in
  // the underlying class.

protected:
  /// Change in guessed objective value for next branch
  double changeInGuessed_;
  /// Pointer back to object
  CbcSimpleIntegerDynamicPseudoCost *object_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
