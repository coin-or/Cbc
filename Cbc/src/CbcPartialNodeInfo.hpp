// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/24/09 carved from CbcNode

#ifndef CbcPartialNodeInfo_H
#define CbcPartialNodeInfo_H

#include <string>
#include <vector>

#include "CoinWarmStartBasis.hpp"
#include "CoinSearchTree.hpp"
#include "CbcBranchBase.hpp"
#include "CbcNodeInfo.hpp"

class OsiSolverInterface;
class OsiSolverBranch;

class OsiCuts;
class OsiRowCut;
class OsiRowCutDebugger;
class CoinWarmStartBasis;
class CbcCountRowCut;
class CbcModel;
class CbcNode;
class CbcSubProblem;
class CbcGeneralBranchingObject;
/** \brief Holds information for recreating a subproblem by incremental change
	   from the parent.

  A CbcPartialNodeInfo object contains changes to the bounds and basis, and
  additional cuts, required to recreate a subproblem by modifying and
  augmenting the parent subproblem.
*/

class CbcPartialNodeInfo : public CbcNodeInfo {

public:
  /** \brief Modify model according to information at node

        The routine modifies the model according to bound and basis change
        information at node and adds any cuts to the addCuts array.
    */
  virtual void applyToModel(CbcModel *model, CoinWarmStartBasis *&basis,
    CbcCountRowCut **addCuts,
    int &currentNumberCuts) const;

  /// Just apply bounds to one variable - force means overwrite by lower,upper (1=>infeasible)
  virtual int applyBounds(int iColumn, double &lower, double &upper, int force);
  /** Builds up row basis backwards (until original model).
        Returns NULL or previous one to apply .
        Depends on Free being 0 and impossible for cuts
    */
  virtual CbcNodeInfo *buildRowBasis(CoinWarmStartBasis &basis) const;
  // Default Constructor
  CbcPartialNodeInfo();

  // Constructor from current state
  CbcPartialNodeInfo(CbcNodeInfo *parent, CbcNode *owner,
    int numberChangedBounds, const int *variables,
    const double *boundChanges,
    const CoinWarmStartDiff *basisDiff);

  // Copy constructor
  CbcPartialNodeInfo(const CbcPartialNodeInfo &);

  // Destructor
  ~CbcPartialNodeInfo();

  /// Clone
  virtual CbcNodeInfo *clone() const;
  /// Basis diff information
  inline const CoinWarmStartDiff *basisDiff() const
  {
    return basisDiff_;
  }
  /// Which variable (top bit if upper bound changing)
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
  inline int numberChangedBounds() const
  {
    return numberChangedBounds_;
  }

protected:
  /* Data values */

  /// Basis diff information
  CoinWarmStartDiff *basisDiff_;
  /// Which variable (top bit if upper bound changing)
  int *variables_;
  // New bound
  double *newBounds_;
  /// Number of bound changes
  int numberChangedBounds_;

private:
  /// Illegal Assignment operator
  CbcPartialNodeInfo &operator=(const CbcPartialNodeInfo &rhs);
};

#endif //CbcPartialNodeInfo_H

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
