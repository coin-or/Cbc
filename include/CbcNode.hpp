// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcNode_H
#define CbcNode_H

#include <string>
#include <vector>

#include "CoinWarmStartBasis.hpp"
#include "CbcBranchBase.hpp"

class OsiSolverInterface;

class OsiCuts;
class OsiRowCut;
class OsiRowCutDebugger;
class CoinWarmStartBasis;
class CbcCountRowCut;
class CbcModel;
class CbcNode;

//#############################################################################
/** Information required to recreate the subproblem at this node

  When a subproblem is initially created, it is represented by an CbcNode
  object and an attached CbcNodeInfo object.

  The CbcNode contains information needed while the subproblem remains live.
  The CbcNode is deleted when the last branch arm has been evaluated.

  The CbcNodeInfo contains information required to maintain the branch-and-cut
  search tree structure (links and reference counts) and to recreate the
  subproblem for this node (basis, variable bounds, cutting planes). An
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

  CbcNodeInfo objects come in two flavours. An CbcFullNodeInfo object contains
  a full record of the information required to recreate a subproblem.
  An CbcPartialNodeInfo object expresses this information in terms of
  differences from the parent.
*/

class CbcNodeInfo {

public:

/** \name Constructors & destructors */
//@{
  /** Default Constructor

    Creates an empty NodeInfo object.
  */
  CbcNodeInfo ();

  /// Copy constructor 
  CbcNodeInfo ( const CbcNodeInfo &);
   

  /** Construct with parent

    Creates a NodeInfo object which knows its parent and assumes it will
    in turn have two children.
  */
  CbcNodeInfo (CbcNodeInfo * parent);
   
  /** Construct with parent and owner

    As for `construct with parent', and attached to \p owner.
  */
  CbcNodeInfo (CbcNodeInfo * parent, CbcNode * owner);

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
  virtual void applyToModel (CbcModel *model, CoinWarmStartBasis *&basis,
			     CbcCountRowCut **addCuts,
			     int &currentNumberCuts) const = 0 ;

  /** Builds up row basis backwards (until original model).
      Returns NULL or previous one to apply .
      Depends on Free being 0 and impossible for cuts
  */
  virtual CbcNodeInfo * buildRowBasis(CoinWarmStartBasis & basis) const = 0;
  /// Clone
  virtual CbcNodeInfo * clone() const = 0;

  /// Increment number of references
  inline void increment(int amount=1)
  {numberPointingToThis_+=amount;};

  /// Decrement number of references and return number left
  inline int decrement(int amount=1)
  {numberPointingToThis_-=amount;return numberPointingToThis_;};

  /** Initialize reference counts

    Initialize the reference counts used for tree maintenance.
  */

  inline void initializeInfo(int number)
  {numberPointingToThis_=number;numberBranchesLeft_=number;};

  /// Return number of branches left in object
  inline int numberBranchesLeft() const
  {return numberBranchesLeft_;};

  /// Return number of objects pointing to this
  inline int numberPointingToThis() const
  {return numberPointingToThis_;};

  /// Say one branch taken 
  inline int branchedOn()
  {numberPointingToThis_--;numberBranchesLeft_--;return numberBranchesLeft_;};

  /// Say thrown away
  inline void throwAway()
  {numberPointingToThis_-=numberBranchesLeft_;numberBranchesLeft_=0;};

  /// Parent of this
  CbcNodeInfo * parent() const
  {return parent_;};

  void addCuts(OsiCuts & cuts,int numberToBranch, int * whichGenerator);
  void addCuts(int numberCuts, CbcCountRowCut ** cuts,int numberToBranch);
  /** Delete cuts (decrements counts)
      Slow unless cuts in same order as saved
  */
  void deleteCuts(int numberToDelete,CbcCountRowCut ** cuts);
  void deleteCuts(int numberToDelete,int * which);

  /// Really delete a cut
  void deleteCut(int whichOne);

  /// Decrement active cut counts
  void decrementCuts(int change=1);

  /// Decrement all active cut counts in chain starting at parent
  void decrementParentCuts(int change=1);

  /// Increment all active cut counts in parent chain
  void incrementParentCuts(int change=1);

  /// Array of pointers to cuts
  inline CbcCountRowCut ** cuts() const
  {return cuts_;};

  /// Number of row cuts (this node)
  inline int numberCuts() const
  {return numberCuts_;};
  inline void setNumberCuts(int value)
  {numberCuts_=value;};

  /// Set owner null
  inline void nullOwner()
  { owner_=NULL;};
protected:

  /** Number of other nodes pointing to this node.

    Number of existing and potential search tree nodes pointing to this node.
    `Existing' means referenced by #parent_ of some other CbcNodeInfo.
    `Potential' means children still to be created (#numberBranchesLeft_ of
    this CbcNodeInfo).
  */
  int numberPointingToThis_;

  /// parent
  CbcNodeInfo * parent_;

  /// Owner
  CbcNode * owner_;

  /// Number of row cuts (this node)
  int numberCuts_;

  /// Array of pointers to cuts
  CbcCountRowCut ** cuts_;

  /** Number of rows in problem (before these cuts).  This
      means that for top of chain it must be rows at continuous */
  int numberRows_;

  /** Number of branch arms left to explore at this node
  
    \todo There seems to be redundancy between this field and
	  CbcBranchingObject::numberBranchesLeft_. It'd be good to sort out if
	  both are necessary.
  */
  int numberBranchesLeft_;
      
private:
  
  /// Illegal Assignment operator 
  CbcNodeInfo & operator=(const CbcNodeInfo& rhs);
  
};

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
  virtual void applyToModel (CbcModel *model, CoinWarmStartBasis *&basis,
			     CbcCountRowCut **addCuts,
			     int &currentNumberCuts) const ;

  /** Builds up row basis backwards (until original model).
      Returns NULL or previous one to apply .
      Depends on Free being 0 and impossible for cuts
  */
  virtual CbcNodeInfo * buildRowBasis(CoinWarmStartBasis & basis) const ;
  // Default Constructor 
  CbcFullNodeInfo ();

  /** Constructor from continuous or satisfied
  */
  CbcFullNodeInfo (CbcModel * model,
		   int numberRowsAtContinuous);
  
  // Copy constructor 
  CbcFullNodeInfo ( const CbcFullNodeInfo &);
   
  // Destructor 
  ~CbcFullNodeInfo ();
  
  /// Clone
  virtual CbcNodeInfo * clone() const;
protected:
  // Data
  /** Full basis 

    This MUST BE A POINTER to avoid cutting extra information in derived
    warm start classes.
  */
  CoinWarmStartBasis *basis_;
  int numberIntegers_;
  // Bounds stored in full
  double * lower_;
  double * upper_;
private:
  /// Illegal Assignment operator 
  CbcFullNodeInfo & operator=(const CbcFullNodeInfo& rhs);
};



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
  virtual void applyToModel (CbcModel *model, CoinWarmStartBasis *&basis,
			     CbcCountRowCut **addCuts,
			     int &currentNumberCuts) const ;

  /** Builds up row basis backwards (until original model).
      Returns NULL or previous one to apply .
      Depends on Free being 0 and impossible for cuts
  */
  virtual CbcNodeInfo * buildRowBasis(CoinWarmStartBasis & basis ) const ;
  // Default Constructor 
  CbcPartialNodeInfo ();

  // Constructor from current state 
  CbcPartialNodeInfo (CbcNodeInfo * parent, CbcNode * owner,
		int numberChangedBounds,const int * variables,
		const double * boundChanges,
		const CoinWarmStartDiff *basisDiff) ;
  
  // Copy constructor 
  CbcPartialNodeInfo ( const CbcPartialNodeInfo &);
   
  // Destructor 
  ~CbcPartialNodeInfo ();
  
  /// Clone
  virtual CbcNodeInfo * clone() const;
private:
  /* Data values */

  /// Basis diff information
  CoinWarmStartDiff *basisDiff_ ;
  /// Which variable (top bit if upper bound changing)
  int * variables_;
  // New bound
  double * newBounds_;
  /// Number of bound changes
  int numberChangedBounds_;
private:
  
  /// Illegal Assignment operator 
  CbcPartialNodeInfo & operator=(const CbcPartialNodeInfo& rhs);
};



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

class CbcNode  {
 
public:
    
  /// Default Constructor 
  CbcNode ();

  /// Construct and increment parent reference count
  CbcNode (CbcModel * model, CbcNode * lastNode);

  /// Copy constructor 
  CbcNode (const CbcNode &);
   
  /// Assignment operator 
  CbcNode & operator= (const CbcNode& rhs);

  /// Destructor 
  ~CbcNode ();

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
  createInfo(CbcModel * model,
	     CbcNode * lastNode,
	     const CoinWarmStartBasis *lastws,
	     const double * lastLower, const double * lastUpper,
	     int numberOldActiveCuts,int numberNewCuts);
  
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
  int chooseBranch (CbcModel * model,
		    CbcNode * lastNode,
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
    </ul>
  */
  int chooseDynamicBranch (CbcModel * model,
		    CbcNode * lastNode,
                    int numberPassesLeft);
  
  /// Decrement active cut counts
  void decrementCuts(int change=1);

  /// Decrement all active cut counts in chain starting at parent
  void decrementParentCuts(int change=1);

  /// Nulls out node info
  void nullNodeInfo();
  /** Initialize reference counts in attached CbcNodeInfo
  
    This is a convenience routine, which will initialize the reference counts
    in the attached CbcNodeInfo object based on the attached
    CbcBranchingObject.

    \sa CbcNodeInfo::initializeInfo(int).
  */
  void initializeInfo();

  /// Does next branch and updates state
  int branch();

  // Information to make basis and bounds
  inline CbcNodeInfo * nodeInfo() const
  {return nodeInfo_;};

  // Objective value
  inline double objectiveValue() const
  { return objectiveValue_;};
  inline void setObjectiveValue(double value)
  { objectiveValue_=value;};
  /// Number of arms defined for the attached CbcBranchingObject.
  inline int numberBranches() const
  { if (branch_)
      return (branch_->numberBranches()) ;
    else
      return (-1) ; } ;

  /** Branching `variable' associated with the attached CbcBranchingObject.

    Check CbcBranchingObject::variable() for a longer explanation of
    `variable'.
  */
  inline int variable() const
  {if (branch_) return branch_->variable();else return -1;};

  /* Active arm of the attached CbcBranchingObject.
  
   In the simplest instance, coded -1 for the down arm of the branch, +1 for
   the up arm. But see CbcBranchingObject::way() 
     Use nodeInfo--.numberBranchesLeft_ to see how active
  */
  inline int way() const
  {if (branch_) return branch_->way();else return 0;};
  /// Depth in branch-and-cut search tree
  inline int depth() const
  {return depth_;};
  /// Get the number of objects unsatisfied at this node.
  inline int numberUnsatisfied() const
  {return numberUnsatisfied_;};
  /// The node number
  inline int nodeNumber() const
  { return nodeNumber_;};
  inline void setNodeNumber(int node)
  { nodeNumber_=node;};

  // Guessed objective value (for solution)
  inline double guessedObjectiveValue() const
  {return guessedObjectiveValue_;};
  inline void setGuessedObjectiveValue(double value)
  {guessedObjectiveValue_=value;};
  /// Branching object for this node
  const CbcBranchingObject * branchingObject() const
  { return branch_;};
  /// Modifiable branching object for this node
  CbcBranchingObject * modifiableBranchingObject() const
  { return branch_;};

private:
  // Data
  /// Information to make basis and bounds
  CbcNodeInfo * nodeInfo_;
  // Objective value
  double objectiveValue_;
  // Guessed satisfied Objective value
  double guessedObjectiveValue_;
  /// Branching object for this node
  CbcBranchingObject * branch_;
  /// Depth of the node in the search tree
  int depth_;
  /// The number of objects unsatisfied at this node.
  int numberUnsatisfied_;
  /// The node number
  int nodeNumber_;
};


#endif
