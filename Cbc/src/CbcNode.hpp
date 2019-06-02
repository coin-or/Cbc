/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcNode_H
#define CbcNode_H

#include <string>
#include <vector>

#include "CoinWarmStartBasis.hpp"
#include "CoinSearchTree.hpp"
#include "CbcBranchBase.hpp"
#include "CbcNodeInfo.hpp"
#include "CbcFullNodeInfo.hpp"
#include "CbcPartialNodeInfo.hpp"

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

/** Information required while the node is live

  When a subproblem is initially created, it is represented by an CbcNode
  object and an attached CbcNodeInfo object.

  The CbcNode contains information (depth, branching instructions), that's
  needed while the subproblem remains `live', <i>i.e.</i>, while the
  subproblem is not fathomed and there are branch arms still be be
  evaluated.  The CbcNode is deleted when the last branch arm has been
  evaluated.

  The CbcNodeInfo object contains the information needed to maintain the
  search tree and recreate the subproblem for the node. It remains in
  existence until there are no nodes remaining in the subtree rooted at this
  node.
*/

class CbcNode : public CoinTreeNode {

public:
  /// Default Constructor
  CbcNode();

  /// Construct and increment parent reference count
  CbcNode(CbcModel *model, CbcNode *lastNode);

  /// Copy constructor
  CbcNode(const CbcNode &);

  /// Assignment operator
  CbcNode &operator=(const CbcNode &rhs);

  /// Destructor
  ~CbcNode();

  /** Create a description of the subproblem at this node

      The CbcNodeInfo structure holds the information (basis & variable bounds)
      required to recreate the subproblem for this node. It also links the node
      to its parent (via the parent's CbcNodeInfo object).

      If lastNode == NULL, a CbcFullNodeInfo object will be created. All
      parameters except \p model are unused.

      If lastNode != NULL, a CbcPartialNodeInfo object will be created. Basis and
      bounds information will be stored in the form of differences between the
      parent subproblem and this subproblem.
      (More precisely, \p lastws, \p lastUpper, \p lastLower,
      \p numberOldActiveCuts, and \p numberNewCuts are used.)
    */
  void
  createInfo(CbcModel *model,
    CbcNode *lastNode,
    const CoinWarmStartBasis *lastws,
    const double *lastLower, const double *lastUpper,
    int numberOldActiveCuts, int numberNewCuts);

  /** Create a branching object for the node

      The routine scans the object list of the model and selects a set of
      unsatisfied objects as candidates for branching. The candidates are
      evaluated, and an appropriate branch object is installed.

      The numberPassesLeft is decremented to stop fixing one variable each time
      and going on and on (e.g. for stock cutting, air crew scheduling)

      If evaluation determines that an object is monotone or infeasible,
      the routine returns immediately. In the case of a monotone object,
      the branch object has already been called to modify the model.

      Return value:
      <ul>
        <li>  0: A branching object has been installed
        <li> -1: A monotone object was discovered
        <li> -2: An infeasible object was discovered
      </ul>
    */
  int chooseBranch(CbcModel *model,
    CbcNode *lastNode,
    int numberPassesLeft);
  /** Create a branching object for the node - when dynamic pseudo costs

      The routine scans the object list of the model and selects a set of
      unsatisfied objects as candidates for branching. The candidates are
      evaluated, and an appropriate branch object is installed.
      This version gives preference in evaluation to variables which
      have not been evaluated many times.  It also uses numberStrong
      to say give up if last few tries have not changed incumbent.
      See Achterberg, Koch and Martin.

      The numberPassesLeft is decremented to stop fixing one variable each time
      and going on and on (e.g. for stock cutting, air crew scheduling)

      If evaluation determines that an object is monotone or infeasible,
      the routine returns immediately. In the case of a monotone object,
      the branch object has already been called to modify the model.

      Return value:
      <ul>
        <li>  0: A branching object has been installed
        <li> -1: A monotone object was discovered
        <li> -2: An infeasible object was discovered
        <li> >0: Number of quich branching objects (and branches will be non NULL)
      </ul>
    */
  int chooseDynamicBranch(CbcModel *model,
    CbcNode *lastNode,
    OsiSolverBranch *&branches,
    int numberPassesLeft);
  /** Create a branching object for the node

      The routine scans the object list of the model and selects a set of
      unsatisfied objects as candidates for branching. The candidates are
      evaluated, and an appropriate branch object is installed.

      The numberPassesLeft is decremented to stop fixing one variable each time
      and going on and on (e.g. for stock cutting, air crew scheduling)

      If evaluation determines that an object is monotone or infeasible,
      the routine returns immediately. In the case of a monotone object,
      the branch object has already been called to modify the model.

      Return value:
      <ul>
        <li>  0: A branching object has been installed
        <li> -1: A monotone object was discovered
        <li> -2: An infeasible object was discovered
      </ul>
      Branch state:
      <ul>
        <li> -1: start
        <li> -1: A monotone object was discovered
        <li> -2: An infeasible object was discovered
      </ul>
    */
  int chooseOsiBranch(CbcModel *model,
    CbcNode *lastNode,
    OsiBranchingInformation *usefulInfo,
    int branchState);
  /** Create a branching object for the node

      The routine scans the object list of the model and selects a set of
      unsatisfied objects as candidates for branching. It then solves a
      series of problems and a CbcGeneral branch object is installed.

      If evaluation determines that an object is infeasible,
      the routine returns immediately.

      Return value:
      <ul>
        <li>  0: A branching object has been installed
        <li> -2: An infeasible object was discovered
      </ul>
    */
  int chooseClpBranch(CbcModel *model,
    CbcNode *lastNode);
  int analyze(CbcModel *model, double *results);
  /// Decrement active cut counts
  void decrementCuts(int change = 1);

  /// Decrement all active cut counts in chain starting at parent
  void decrementParentCuts(CbcModel *model, int change = 1);

  /// Nulls out node info
  void nullNodeInfo();
  /** Initialize reference counts in attached CbcNodeInfo

      This is a convenience routine, which will initialize the reference counts
      in the attached CbcNodeInfo object based on the attached
      OsiBranchingObject.

      \sa CbcNodeInfo::initializeInfo(int).
    */
  void initializeInfo();

  /// Does next branch and updates state
  int branch(OsiSolverInterface *solver);

  /** Double checks in case node can change its mind!
        Returns objective value
        Can change objective etc */
  double checkIsCutoff(double cutoff);
  // Information to make basis and bounds
  inline CbcNodeInfo *nodeInfo() const
  {
    return nodeInfo_;
  }

  // Objective value
  inline double objectiveValue() const
  {
    return objectiveValue_;
  }
  inline void setObjectiveValue(double value)
  {
    objectiveValue_ = value;
  }
  /// Number of arms defined for the attached OsiBranchingObject.
  inline int numberBranches() const
  {
    if (branch_)
      return (branch_->numberBranches());
    else
      return (-1);
  }

  /* Active arm of the attached OsiBranchingObject.

     In the simplest instance, coded -1 for the down arm of the branch, +1 for
     the up arm. But see OsiBranchingObject::way()
       Use nodeInfo--.numberBranchesLeft_ to see how active
    */
  int way() const;
  /// Depth in branch-and-cut search tree
  inline int depth() const
  {
    return depth_;
  }
  /// Set depth in branch-and-cut search tree
  inline void setDepth(int value)
  {
    depth_ = value;
  }
  /// Get the number of objects unsatisfied at this node.
  inline int numberUnsatisfied() const
  {
    return numberUnsatisfied_;
  }
  /// Set the number of objects unsatisfied at this node.
  inline void setNumberUnsatisfied(int value)
  {
    numberUnsatisfied_ = value;
  }
  /// Get sum of "infeasibilities" reported by each object
  inline double sumInfeasibilities() const
  {
    return sumInfeasibilities_;
  }
  /// Set sum of "infeasibilities" reported by each object
  inline void setSumInfeasibilities(double value)
  {
    sumInfeasibilities_ = value;
  }
  // Guessed objective value (for solution)
  inline double guessedObjectiveValue() const
  {
    return guessedObjectiveValue_;
  }
  inline void setGuessedObjectiveValue(double value)
  {
    guessedObjectiveValue_ = value;
  }
  /// Branching object for this node
  inline const OsiBranchingObject *branchingObject() const
  {
    return branch_;
  }
  /// Modifiable branching object for this node
  inline OsiBranchingObject *modifiableBranchingObject() const
  {
    return branch_;
  }
  /// Set branching object for this node (takes ownership)
  inline void setBranchingObject(OsiBranchingObject *branchingObject)
  {
    branch_ = branchingObject;
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
  /// Returns true if on tree
  inline bool onTree() const
  {
    return (state_ & 1) != 0;
  }
  /// Sets true if on tree
  inline void setOnTree(bool yesNo)
  {
    if (yesNo)
      state_ |= 1;
    else
      state_ &= ~1;
  }
  /// Returns true if active
  inline bool active() const
  {
    return (state_ & 2) != 0;
  }
  /// Sets true if active
  inline void setActive(bool yesNo)
  {
    if (yesNo)
      state_ |= 2;
    else
      state_ &= ~2;
  }
  /// Get state (really for debug)
  inline int getState() const
  {
    return state_;
  }
  /// Set state (really for debug)
  inline void setState(int value)
  {
    state_ = value;
  }
  /// Print
  void print() const;
  /// Debug
  inline void checkInfo() const
  {
    assert(nodeInfo_->numberBranchesLeft() == branch_->numberBranchesLeft());
  }

private:
  // Data
  /// Information to make basis and bounds
  CbcNodeInfo *nodeInfo_;
  /// Objective value
  double objectiveValue_;
  /// Guessed satisfied Objective value
  double guessedObjectiveValue_;
  /// Sum of "infeasibilities" reported by each object
  double sumInfeasibilities_;
  /// Branching object for this node
  OsiBranchingObject *branch_;
  /// Depth of the node in the search tree
  int depth_;
  /// The number of objects unsatisfied at this node.
  int numberUnsatisfied_;
  /// The node number
  int nodeNumber_;
  /** State
        1 - on tree
        2 - active
    */
  int state_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
