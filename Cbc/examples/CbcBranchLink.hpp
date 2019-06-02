// $Id$
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcBranchLink_H
#define CbcBranchLink_H

#include "CbcBranchBase.hpp"

/** Define Special Linked Ordered Sets.

*/

class CbcLink : public CbcObject {

public:
  // Default Constructor
  CbcLink();

  /** Useful constructor - A valid solution is if all variables are zero
      apart from k*numberLink to (k+1)*numberLink-1 where k is 0 through
      numberInSet-1.  The length of weights array is numberInSet.
      For this constructor the variables in matrix are the numberInSet*numberLink
      starting at first. If weights null then 0,1,2..
  */
  CbcLink(CbcModel *model, int numberMembers,
    int numberLinks, int first,
    const double *weights, int setNumber);
  /** Useful constructor - A valid solution is if all variables are zero
      apart from k*numberLink to (k+1)*numberLink-1 where k is 0 through
      numberInSet-1.  The length of weights array is numberInSet.
      For this constructor the variables are given by list - grouped.
      If weights null then 0,1,2..
  */
  CbcLink(CbcModel *model, int numberMembers,
    int numberLinks, int typeSOS, const int *which,
    const double *weights, int setNumber);

  // Copy constructor
  CbcLink(const CbcLink &);

  /// Clone
  virtual CbcObject *clone() const;

  // Assignment operator
  CbcLink &operator=(const CbcLink &rhs);

  // Destructor
  ~CbcLink();

  /// Infeasibility - large is 0.5
  virtual double infeasibility(int &preferredWay) const;

  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion();
  /// Creates a branching object
  virtual CbcBranchingObject *createCbcBranch(OsiSolverInterface *solver, const OsiBranchingInformation *info, int way);

  /// Number of members
  inline int numberMembers() const
  {
    return numberMembers_;
  }

  /// Number of links for each member
  inline int numberLinks() const
  {
    return numberLinks_;
  }

  /// Which variables
  inline const int *which() const
  {
    return which_;
  }

  /** Array of weights */
  inline const double *weights() const
  {
    return weights_;
  }

private:
  /// data

  /// Weights
  double *weights_;

  /// Number of members
  int numberMembers_;
  /// Number of links
  int numberLinks_;
  /// Members
  int *which_;
  /// Type 1 or 2
  int sosType_;
};
/** Branching object for Special ordered sets

    Variable_ is the set id number (redundant, as the object also holds a
    pointer to the set.
 */
class CbcLinkBranchingObject : public CbcBranchingObject {

public:
  // Default Constructor
  CbcLinkBranchingObject();

  // Useful constructor
  CbcLinkBranchingObject(CbcModel *model, const CbcLink *set,
    int way,
    double separator);

  // Copy constructor
  CbcLinkBranchingObject(const CbcLinkBranchingObject &);

  // Assignment operator
  CbcLinkBranchingObject &operator=(const CbcLinkBranchingObject &rhs);

  /// Clone
  virtual CbcBranchingObject *clone() const;

  // Destructor
  virtual ~CbcLinkBranchingObject();

  /// Does next branch and updates state
  virtual double branch();

  /** \brief Print something about branch - only if log level high
  */
  virtual void print();
  /** Return the type (an integer identifier) of \c this */
  virtual CbcBranchObjType type() const
  {
    return CbcBranchObjType(0);
  } /*FIXME what type() should be returned here? */

  /** Compare the \c this with \c brObj. \c this and \c brObj must be os the
      same type and must have the same original object, but they may have
      different feasible regions.
      Return the appropriate CbcRangeCompare value (first argument being the
      sub/superset if that's the case). In case of overlap (and if \c
      replaceIfOverlap is true) replace the current branching object with one
      whose feasible region is the overlap.
   */
  virtual CbcRangeCompare compareBranchingObject(const CbcBranchingObject *brObj, const bool replaceIfOverlap = false);

private:
  /// data
  const CbcLink *set_;
  /// separator
  double separator_;
};
#endif
