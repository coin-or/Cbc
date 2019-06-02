// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#ifndef CbcDummyBranchingObject_H
#define CbcDummyBranchingObject_H

#include "CbcBranchBase.hpp"
/** Dummy branching object

  This object specifies a one-way dummy branch.
  This is so one can carry on branching even when it looks feasible
*/

class CbcDummyBranchingObject : public CbcBranchingObject {

public:
  /// Default constructor
  CbcDummyBranchingObject(CbcModel *model = NULL);

  /// Copy constructor
  CbcDummyBranchingObject(const CbcDummyBranchingObject &);

  /// Assignment operator
  CbcDummyBranchingObject &operator=(const CbcDummyBranchingObject &rhs);

  /// Clone
  virtual CbcBranchingObject *clone() const;

  /// Destructor
  virtual ~CbcDummyBranchingObject();

  using CbcBranchingObject::branch;
  /** \brief Dummy branch
    */
  virtual double branch();

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

  /** Return the type (an integer identifier) of \c this */
  virtual CbcBranchObjType type() const
  {
    return DummyBranchObj;
  }

  /** Compare the original object of \c this with the original object of \c
        brObj. Assumes that there is an ordering of the original objects.
        This method should be invoked only if \c this and brObj are of the same
        type.
        Return negative/0/positive depending on whether \c this is
        smaller/same/larger than the argument.
    */
  virtual int compareOriginalObject(const CbcBranchingObject *brObj) const;

  /** Compare the \c this with \c brObj. \c this and \c brObj must be os the
        same type and must have the same original object, but they may have
        different feasible regions.
        Return the appropriate CbcRangeCompare value (first argument being the
        sub/superset if that's the case). In case of overlap (and if \c
        replaceIfOverlap is true) replace the current branching object with one
        whose feasible region is the overlap.
     */
  virtual CbcRangeCompare compareBranchingObject(const CbcBranchingObject *brObj, const bool replaceIfOverlap = false);
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
