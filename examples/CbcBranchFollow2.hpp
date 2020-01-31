// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcBranchFollowOn2_H
#define CbcBranchFollowOn2_H

#include "CbcBranchActual.hpp"
#include "CoinPackedMatrix.hpp"

/** Define a follow on class.
    The idea of this is that in air-crew scheduling problems crew may fly in on flight A
    and out on flight B or on some other flight.  A useful branch is one which on one side
    fixes all which go out on flight B to 0, while the other branch fixes all those that do NOT
    go out on flight B to 0.

    This tries to generalize so that cuts are produced with sum aij xj <= bi on each side. 
    It should be intelligent enough to fix if things can be fixed.
    We also need to make sure branch cuts work properly (i.e. persistence).

    This branching rule should be in addition to normal rules and have a high priority.
*/

class CbcFollowOn2 : public CbcObject {

public:
  // Default Constructor
  CbcFollowOn2();

  /** Useful constructor
  */
  CbcFollowOn2(CbcModel *model);

  // Copy constructor
  CbcFollowOn2(const CbcFollowOn2 &);

  /// Clone
  virtual CbcObject *clone() const;

  // Assignment operator
  CbcFollowOn2 &operator=(const CbcFollowOn2 &rhs);

  // Destructor
  ~CbcFollowOn2();

  /// Infeasibility - large is 0.5
  virtual double infeasibility(int &preferredWay) const;

  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion();
  /// Creates a branching object
  virtual CbcBranchingObject *createBranch(int way);
  /** As some computation is needed in more than one place - returns row.
      Also returns other row and effective rhs (so we can know if cut)
  */
  virtual int gutsOfFollowOn2(int &otherRow, int &preferredWay,
    int &effectiveRhs) const;

  /// get and set for maximum rhws (affects cuts as branch)
  inline int maximumRhs() const
  {
    return maximumRhs_;
  }
  inline void setMaximumRhs(int value)
  {
    maximumRhs_ = value;
  }

protected:
  /// data
  /// Matrix
  CoinPackedMatrix matrix_;
  /// Matrix by row
  CoinPackedMatrix matrixByRow_;
  /// Possible rhs (if 0 then not possible)
  int *rhs_;
  /// If >1 then allow cuts if effective rhs <= this
  int maximumRhs_;
};

#endif
