/* $Id$ */
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcStatistics_H
#define CbcStatistics_H

#include "CbcModel.hpp"

/** For gathering statistics */

class CbcStatistics {
public:
  // Default Constructor
  CbcStatistics();
  // Branch
  CbcStatistics(CbcNode *node, CbcModel *model);

  ~CbcStatistics();
  // Copy
  CbcStatistics(const CbcStatistics &rhs);
  // Assignment
  CbcStatistics &operator=(const CbcStatistics &rhs);
  // Update at end of branch
  void endOfBranch(int numberIterations, double objectiveValue);
  // Update number of infeasibilities
  void updateInfeasibility(int numberInfeasibilities);
  // Branch found to be infeasible by chooseBranch
  void sayInfeasible();
  // Just prints
  void print(const int *sequenceLookup = NULL) const;
  // Node number
  inline int node() const
  {
    return id_;
  }
  // Parent node number
  inline int parentNode() const
  {
    return parentId_;
  }
  // depth
  inline int depth() const
  {
    return depth_;
  }
  // way
  inline int way() const
  {
    return way_;
  }
  // value
  inline double value() const
  {
    return value_;
  }
  // starting objective
  inline double startingObjective() const
  {
    return startingObjective_;
  }
  // Unsatisfied at beginning
  inline int startingInfeasibility() const
  {
    return startingInfeasibility_;
  }
  // starting objective
  inline double endingObjective() const
  {
    return endingObjective_;
  }
  // Unsatisfied at end
  inline int endingInfeasibility() const
  {
    return endingInfeasibility_;
  }
  // Number iterations
  inline int numberIterations() const
  {
    return numberIterations_;
  }

protected:
  // Data
  /// Value
  double value_;
  /// Starting objective
  double startingObjective_;
  /// Ending objective
  double endingObjective_;
  /// id
  int id_;
  /// parent id
  int parentId_;
  /// way -1 or +1 is first branch -10 or +10 is second branch
  int way_;
  /// sequence number branched on
  int sequence_;
  /// depth
  int depth_;
  /// starting number of integer infeasibilities
  int startingInfeasibility_;
  /// ending number of integer infeasibilities
  int endingInfeasibility_;
  /// number of iterations
  int numberIterations_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
