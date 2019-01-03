// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/24/09 carved from CbcNode

#ifndef CbcNodeInfo_H
#define CbcNodeInfo_H

#include <string>
#include <vector>

#include "CoinWarmStartBasis.hpp"
#include "CoinSearchTree.hpp"
#include "CbcBranchBase.hpp"

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

class CbcNodeInfo {

public:
  /** \name Constructors & destructors */
  //@{
  /** Default Constructor

      Creates an empty NodeInfo object.
    */
  CbcNodeInfo();

  /// Copy constructor
  CbcNodeInfo(const CbcNodeInfo &);

#ifdef JJF_ZERO
  /** Construct with parent

      Creates a NodeInfo object which knows its parent and assumes it will
      in turn have two children.
    */
  CbcNodeInfo(CbcNodeInfo *parent);
#endif

  /** Construct with parent and owner

      As for `construct with parent', and attached to \p owner.
    */
  CbcNodeInfo(CbcNodeInfo *parent, CbcNode *owner);

  /** Destructor

      Note that the destructor will recursively delete the parent if this
      nodeInfo is the last child.
    */
  virtual ~CbcNodeInfo();
  //@}

  /** \brief Modify model according to information at node

        The routine modifies the model according to bound and basis
        information at node and adds any cuts to the addCuts array.
    */
  virtual void applyToModel(CbcModel *model, CoinWarmStartBasis *&basis,
    CbcCountRowCut **addCuts,
    int &currentNumberCuts) const = 0;
  /// Just apply bounds to one variable - force means overwrite by lower,upper (1=>infeasible)
  virtual int applyBounds(int iColumn, double &lower, double &upper, int force) = 0;

  /** Builds up row basis backwards (until original model).
        Returns NULL or previous one to apply .
        Depends on Free being 0 and impossible for cuts
    */
  virtual CbcNodeInfo *buildRowBasis(CoinWarmStartBasis &basis) const = 0;
  /// Clone
  virtual CbcNodeInfo *clone() const = 0;
  /// Called when number branches left down to zero
  virtual void allBranchesGone() {}
#ifndef JJF_ONE
  /// Increment number of references
  inline void increment(int amount = 1)
  {
    numberPointingToThis_ += amount; /*printf("CbcNodeInfo %x incremented by %d to %d\n",this,amount,numberPointingToThis_);*/
  }

  /// Decrement number of references and return number left
  inline int decrement(int amount = 1)
  {
    numberPointingToThis_ -= amount; /*printf("CbcNodeInfo %x decremented by %d to %d\n",this,amount,numberPointingToThis_);*/
    return numberPointingToThis_;
  }
#else
  /// Increment number of references
  void increment(int amount = 1);
  /// Decrement number of references and return number left
  int decrement(int amount = 1);
#endif
  /** Initialize reference counts

      Initialize the reference counts used for tree maintenance.
    */

  inline void initializeInfo(int number)
  {
    numberPointingToThis_ = number;
    numberBranchesLeft_ = number;
  }

  /// Return number of branches left in object
  inline int numberBranchesLeft() const
  {
    return numberBranchesLeft_;
  }

  /// Set number of branches left in object
  inline void setNumberBranchesLeft(int value)
  {
    numberBranchesLeft_ = value;
  }

  /// Return number of objects pointing to this
  inline int numberPointingToThis() const
  {
    return numberPointingToThis_;
  }

  /// Set number of objects pointing to this
  inline void setNumberPointingToThis(int number)
  {
    numberPointingToThis_ = number;
  }

  /// Increment number of objects pointing to this
  inline void incrementNumberPointingToThis()
  {
    numberPointingToThis_++;
  }

  /// Say one branch taken
  inline int branchedOn()
  {
    numberPointingToThis_--;
    numberBranchesLeft_--;
    return numberBranchesLeft_;
  }

  /// Say thrown away
  inline void throwAway()
  {
    numberPointingToThis_ -= numberBranchesLeft_;
    numberBranchesLeft_ = 0;
  }

  /// Parent of this
  CbcNodeInfo *parent() const
  {
    return parent_;
  }
  /// Set parent null
  inline void nullParent()
  {
    parent_ = NULL;
  }

  void addCuts(OsiCuts &cuts, int numberToBranch, //int * whichGenerator,
    int numberPointingToThis);
  void addCuts(int numberCuts, CbcCountRowCut **cuts, int numberToBranch);
  /** Delete cuts (decrements counts)
        Slow unless cuts in same order as saved
    */
  void deleteCuts(int numberToDelete, CbcCountRowCut **cuts);
  void deleteCuts(int numberToDelete, int *which);

  /// Really delete a cut
  void deleteCut(int whichOne);

  /// Decrement active cut counts
  void decrementCuts(int change = 1);

  /// Increment active cut counts
  void incrementCuts(int change = 1);

  /// Decrement all active cut counts in chain starting at parent
  void decrementParentCuts(CbcModel *model, int change = 1);

  /// Increment all active cut counts in parent chain
  void incrementParentCuts(CbcModel *model, int change = 1);

  /// Array of pointers to cuts
  inline CbcCountRowCut **cuts() const
  {
    return cuts_;
  }

  /// Number of row cuts (this node)
  inline int numberCuts() const
  {
    return numberCuts_;
  }
  inline void setNumberCuts(int value)
  {
    numberCuts_ = value;
  }

  /// Set owner null
  inline void nullOwner()
  {
    owner_ = NULL;
  }
  const inline CbcNode *owner() const
  {
    return owner_;
  }
  inline CbcNode *mutableOwner() const
  {
    return owner_;
  }
  /// The node number
  inline int nodeNumber() const
  {
    return nodeNumber_;
  }
  inline void setNodeNumber(int node)
  {
    nodeNumber_ = node;
  }
  /** Deactivate node information.
        1 - bounds
        2 - cuts
        4 - basis!
	8 - just marked
	16 - symmetry branching worked
    */
  void deactivate(int mode = 3);
  /// Say if normal
  inline bool allActivated() const
  {
    return ((active_ & 7) == 7);
  }
  /// Say if marked
  inline bool marked() const
  {
    return ((active_ & 8) != 0);
  }
  /// Mark
  inline void mark()
  {
    active_ |= 8;
  }
  /// Unmark
  inline void unmark()
  {
    active_ &= ~8;
  }
  /// Get symmetry value (true worked at this node)
  inline bool symmetryWorked() const
  {
    return (active_ & 16) != 0;
  }
  /// Say symmetry worked at this node)
  inline void setSymmetryWorked()
  {
    active_ |= 16;
  }

  /// Branching object for the parent
  inline const OsiBranchingObject *parentBranch() const
  {
    return parentBranch_;
  }
  /// If we need to take off parent based data
  void unsetParentBasedData();

protected:
  /** Number of other nodes pointing to this node.

      Number of existing and potential search tree nodes pointing to this node.
      `Existing' means referenced by #parent_ of some other CbcNodeInfo.
      `Potential' means children still to be created (#numberBranchesLeft_ of
      this CbcNodeInfo).
    */
  int numberPointingToThis_;

  /// parent
  CbcNodeInfo *parent_;

  /// Copy of the branching object of the parent when the node is created
  OsiBranchingObject *parentBranch_;

  /// Owner
  CbcNode *owner_;

  /// Number of row cuts (this node)
  int numberCuts_;

  /// The node number
  int nodeNumber_;

  /// Array of pointers to cuts
  CbcCountRowCut **cuts_;

  /** Number of rows in problem (before these cuts).  This
        means that for top of chain it must be rows at continuous */
  int numberRows_;

  /** Number of branch arms left to explore at this node

      \todo There seems to be redundancy between this field and
        CbcBranchingObject::numberBranchesLeft_. It'd be good to sort out if
        both are necessary.
    */
  int numberBranchesLeft_;
  /** Active node information.
        1 - bounds
        2 - cuts
        4 - basis!
    */
  int active_;

private:
  /// Illegal Assignment operator
  CbcNodeInfo &operator=(const CbcNodeInfo &rhs);

  /// routine common to constructors
  void setParentBasedData();
};

#endif // CbcNodeInfo_H

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
