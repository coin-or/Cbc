// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/9/2009-- carved out of CbcBranchActual

#ifndef CbcSimpleInteger_H
#define CbcSimpleInteger_H

#include "CbcBranchingObject.hpp"

/** Simple branching object for an integer variable

  This object can specify a two-way branch on an integer variable. For each
  arm of the branch, the upper and lower bounds on the variable can be
  independently specified.

  Variable_ holds the index of the integer variable in the integerVariable_
  array of the model.
*/

class CbcIntegerBranchingObject : public CbcBranchingObject {

public:
  /// Default constructor
  CbcIntegerBranchingObject();

  /** Create a standard floor/ceiling branch object

      Specifies a simple two-way branch. Let \p value = x*. One arm of the
      branch will be lb <= x <= floor(x*), the other ceil(x*) <= x <= ub.
      Specify way = -1 to set the object state to perform the down arm first,
      way = 1 for the up arm.
    */
  CbcIntegerBranchingObject(CbcModel *model, int variable,
    int way, double value);

  /** Create a degenerate branch object

      Specifies a `one-way branch'. Calling branch() for this object will
      always result in lowerValue <= x <= upperValue. Used to fix a variable
      when lowerValue = upperValue.
    */

  CbcIntegerBranchingObject(CbcModel *model, int variable, int way,
    double lowerValue, double upperValue);

  /// Copy constructor
  CbcIntegerBranchingObject(const CbcIntegerBranchingObject &);

  /// Assignment operator
  CbcIntegerBranchingObject &operator=(const CbcIntegerBranchingObject &rhs);

  /// Clone
  virtual CbcBranchingObject *clone() const;

  /// Destructor
  virtual ~CbcIntegerBranchingObject();

  /// Does part of constructor
  void fillPart(int variable, int way, double value);
  using CbcBranchingObject::branch;
  /** \brief Sets the bounds for the variable according to the current arm
           of the branch and advances the object state to the next arm.
           Returns change in guessed objective on next branch
    */
  virtual double branch();
  /** Update bounds in solver as in 'branch' and update given bounds.
        branchState is -1 for 'down' +1 for 'up' */
  virtual void fix(OsiSolverInterface *solver,
    double *lower, double *upper,
    int branchState) const;
  /** Change (tighten) bounds in object to reflect bounds in solver.
	Return true if now fixed */
  virtual bool tighten(OsiSolverInterface *);

#ifdef JJF_ZERO
  // No need to override. Default works fine.
  /** Reset every information so that the branching object appears to point to
        the previous child. This method does not need to modify anything in any
        solver. */
  virtual void previousBranch();
#endif

  using CbcBranchingObject::print;
  /** \brief Print something about branch - only if log level high
    */
  virtual void print();

  /// Lower and upper bounds for down branch
  inline const double *downBounds() const
  {
    return down_;
  }
  /// Lower and upper bounds for up branch
  inline const double *upBounds() const
  {
    return up_;
  }
  /// Set lower and upper bounds for down branch
  inline void setDownBounds(const double bounds[2])
  {
    memcpy(down_, bounds, 2 * sizeof(double));
  }
  /// Set lower and upper bounds for up branch
  inline void setUpBounds(const double bounds[2])
  {
    memcpy(up_, bounds, 2 * sizeof(double));
  }
#ifdef FUNNY_BRANCHING
  /** Which variable (top bit if upper bound changing,
        next bit if on down branch */
  inline const int *variables() const
  {
    return variables_;
  }
  // New bound
  inline const double *newBounds() const
  {
    return newBounds_;
  }
  /// Number of bound changes
  inline int numberExtraChangedBounds() const
  {
    return numberExtraChangedBounds_;
  }
  /// Just apply extra bounds to one variable - COIN_DBL_MAX ignore
  int applyExtraBounds(int iColumn, double lower, double upper, int way);
  /// Deactivate bounds for branching
  void deactivate();
  /// Are active bounds for branching
  inline bool active() const
  {
    return (down_[1] != -COIN_DBL_MAX);
  }
#endif

  /** Return the type (an integer identifier) of \c this */
  virtual CbcBranchObjType type() const
  {
    return SimpleIntegerBranchObj;
  }

  /** Compare the \c this with \c brObj. \c this and \c brObj must be os the
        same type and must have the same original object, but they may have
        different feasible regions.
        Return the appropriate CbcRangeCompare value (first argument being the
        sub/superset if that's the case). In case of overlap (and if \c
        replaceIfOverlap is true) replace the current branching object with one
        whose feasible region is the overlap.
     */
  virtual CbcRangeCompare compareBranchingObject(const CbcBranchingObject *brObj, const bool replaceIfOverlap = false);

protected:
  /// Lower [0] and upper [1] bounds for the down arm (way_ = -1)
  double down_[2];
  /// Lower [0] and upper [1] bounds for the up arm (way_ = 1)
  double up_[2];
#ifdef FUNNY_BRANCHING
  /** Which variable (top bit if upper bound changing)
        next bit if changing on down branch only */
  int *variables_;
  // New bound
  double *newBounds_;
  /// Number of Extra bound changes
  int numberExtraChangedBounds_;
#endif
};

/// Define a single integer class

class CbcSimpleInteger : public CbcObject {

public:
  // Default Constructor
  CbcSimpleInteger();

  // Useful constructor - passed model and index
  CbcSimpleInteger(CbcModel *model, int iColumn, double breakEven = 0.5);

  // Useful constructor - passed model and Osi object
  CbcSimpleInteger(CbcModel *model, const OsiSimpleInteger *object);

  // Copy constructor
  CbcSimpleInteger(const CbcSimpleInteger &);

  /// Clone
  virtual CbcObject *clone() const;

  // Assignment operator
  CbcSimpleInteger &operator=(const CbcSimpleInteger &rhs);

  // Destructor
  virtual ~CbcSimpleInteger();
  /// Construct an OsiSimpleInteger object
  OsiSimpleInteger *osiObject() const;
  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiBranchingInformation *info,
    int &preferredWay) const;

  using CbcObject::feasibleRegion;
  /** Set bounds to fix the variable at the current (integer) value.

      Given an integer value, set the lower and upper bounds to fix the
      variable. Returns amount it had to move variable.
    */
  virtual double feasibleRegion(OsiSolverInterface *solver, const OsiBranchingInformation *info) const;

  /** Create a branching object and indicate which way to branch first.

        The branching object has to know how to create branches (fix
        variables, etc.)
    */
  virtual CbcBranchingObject *createCbcBranch(OsiSolverInterface *solver, const OsiBranchingInformation *info, int way);
  /// Fills in a created branching object
  /*virtual*/ void fillCreateBranch(CbcIntegerBranchingObject *branching, const OsiBranchingInformation *info, int way);

  using CbcObject::solverBranch;
  /** Create an OsiSolverBranch object

        This returns NULL if branch not represented by bound changes
    */
  virtual OsiSolverBranch *solverBranch(OsiSolverInterface *solver, const OsiBranchingInformation *info) const;

  /** Set bounds to fix the variable at the current (integer) value.

      Given an integer value, set the lower and upper bounds to fix the
      variable. The algorithm takes a bit of care in order to compensate for
      minor numerical inaccuracy.
    */
  virtual void feasibleRegion();

  /** Column number if single column object -1 otherwise,
        so returns >= 0
        Used by heuristics
    */
  virtual int columnNumber() const;
  /// Set column number
  inline void setColumnNumber(int value)
  {
    columnNumber_ = value;
  }

  /** Reset variable bounds to their original values.

      Bounds may be tightened, so it may be good to be able to set this info in object.
     */
  virtual void resetBounds(const OsiSolverInterface *solver);

  /**  Change column numbers after preprocessing
     */
  virtual void resetSequenceEtc(int numberColumns, const int *originalColumns);
  /// Original bounds
  inline double originalLowerBound() const
  {
    return originalLower_;
  }
  inline void setOriginalLowerBound(double value)
  {
    originalLower_ = value;
  }
  inline double originalUpperBound() const
  {
    return originalUpper_;
  }
  inline void setOriginalUpperBound(double value)
  {
    originalUpper_ = value;
  }
  /// Breakeven e.g 0.7 -> >= 0.7 go up first
  inline double breakEven() const
  {
    return breakEven_;
  }
  /// Set breakeven e.g 0.7 -> >= 0.7 go up first
  inline void setBreakEven(double value)
  {
    breakEven_ = value;
  }

protected:
  /// data

  /// Original lower bound
  double originalLower_;
  /// Original upper bound
  double originalUpper_;
  /// Breakeven i.e. >= this preferred is up
  double breakEven_;
  /// Column number in model
  int columnNumber_;
  /// If -1 down always chosen first, +1 up always, 0 normal
  int preferredWay_;
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
