// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/24/09 carved from CbcNode

#ifndef CbcFullNodeInfo_H
#define CbcFullNodeInfo_H

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

//#############################################################################
/** Information required to recreate the subproblem at this node

  When a subproblem is initially created, it is represented by a CbcNode
  object and an attached CbcNodeInfo object.

  The CbcNode contains information needed while the subproblem remains live.
  The CbcNode is deleted when the last branch arm has been evaluated.

  The CbcNodeInfo contains information required to maintain the branch-and-cut
  search tree structure (links and reference counts) and to recreate the
  subproblem for this node (basis, variable bounds, cutting planes). A
  CbcNodeInfo object remains in existence until all nodes have been pruned from
  the subtree rooted at this node.

  The principle used to maintain the reference count is that the reference
  count is always the sum of all potential and actual children of the node.
  Specifically,
  <ul>
    <li> Once it's determined how the node will branch, the reference count
	 is set to the number of potential children (<i>i.e.</i>, the number
	 of arms of the branch).
    <li> As each child is created by CbcNode::branch() (converting a potential
	 child to the active subproblem), the reference count is decremented.
    <li> If the child survives and will become a node in the search tree
	 (converting the active subproblem into an actual child), increment the
	 reference count.
  </ul>
  Notice that the active subproblem lives in a sort of limbo, neither a
  potential or an actual node in the branch-and-cut tree.

  CbcNodeInfo objects come in two flavours. A CbcFullNodeInfo object contains
  a full record of the information required to recreate a subproblem.
  A CbcPartialNodeInfo object expresses this information in terms of
  differences from the parent.
*/

/** \brief Holds complete information for recreating a subproblem.

  A CbcFullNodeInfo object contains all necessary information (bounds, basis,
  and cuts) required to recreate a subproblem.

  \todo While there's no explicit statement, the code often makes the implicit
	assumption that an CbcFullNodeInfo structure will appear only at the
	root node of the search tree. Things will break if this assumption
	is violated.
*/

class CbcFullNodeInfo : public CbcNodeInfo {

public:
  /** \brief Modify model according to information at node

        The routine modifies the model according to bound information at node,
        creates a new basis according to information at node, but with the size
        passed in through basis, and adds any cuts to the addCuts array.

      \note The basis passed in via basis is solely a vehicle for passing in
        the desired basis size. It will be deleted and a new basis returned.
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
  CbcFullNodeInfo();

  /** Constructor from continuous or satisfied
    */
  CbcFullNodeInfo(CbcModel *model,
    int numberRowsAtContinuous);

  // Copy constructor
  CbcFullNodeInfo(const CbcFullNodeInfo &);

  // Destructor
  ~CbcFullNodeInfo();

  /// Clone
  virtual CbcNodeInfo *clone() const;
  /// Lower bounds
  inline const double *lower() const
  {
    return lower_;
  }
  /// Set a bound
  inline void setColLower(int sequence, double value)
  {
    lower_[sequence] = value;
  }
  /// Mutable lower bounds
  inline double *mutableLower() const
  {
    return lower_;
  }
  /// Upper bounds
  inline const double *upper() const
  {
    return upper_;
  }
  /// Set a bound
  inline void setColUpper(int sequence, double value)
  {
    upper_[sequence] = value;
  }
  /// Mutable upper bounds
  inline double *mutableUpper() const
  {
    return upper_;
  }

protected:
  // Data
  /** Full basis

      This MUST BE A POINTER to avoid cutting extra information in derived
      warm start classes.
    */
  CoinWarmStartBasis *basis_;
  int numberIntegers_;
  // Bounds stored in full
  double *lower_;
  double *upper_;

private:
  /// Illegal Assignment operator
  CbcFullNodeInfo &operator=(const CbcFullNodeInfo &rhs);
};
#endif //CbcFullNodeInfo_H

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
