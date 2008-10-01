// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcNode_H
#define CbcNode_H

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
  CbcNodeInfo ();

  /// Copy constructor 
  CbcNodeInfo ( const CbcNodeInfo &);
   
#if 0
  /** Construct with parent

    Creates a NodeInfo object which knows its parent and assumes it will
    in turn have two children.
  */
  CbcNodeInfo (CbcNodeInfo * parent);
#endif
   
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
  /// Just apply bounds to one variable - force means overwrite by lower,upper (1=>infeasible)
  virtual int applyBounds(int iColumn, double & lower, double & upper,int force) = 0;

  /** Builds up row basis backwards (until original model).
      Returns NULL or previous one to apply .
      Depends on Free being 0 and impossible for cuts
  */
  virtual CbcNodeInfo * buildRowBasis(CoinWarmStartBasis & basis) const = 0;
  /// Clone
  virtual CbcNodeInfo * clone() const = 0;
  /// Called when number branches left down to zero
  virtual void allBranchesGone() {}
#if 1
  /// Increment number of references
  inline void increment(int amount=1)
  {numberPointingToThis_+=amount;/*printf("CbcNodeInfo %x incremented by %d to %d\n",this,amount,numberPointingToThis_);*/}

  /// Decrement number of references and return number left
  inline int decrement(int amount=1)
  {numberPointingToThis_-=amount;/*printf("CbcNodeInfo %x decremented by %d to %d\n",this,amount,numberPointingToThis_);*/return numberPointingToThis_;}
#else
  /// Increment number of references
  void increment(int amount=1);
  /// Decrement number of references and return number left
  int decrement(int amount=1);
#endif
  /** Initialize reference counts

    Initialize the reference counts used for tree maintenance.
  */

  inline void initializeInfo(int number)
  {numberPointingToThis_=number;numberBranchesLeft_=number;}

  /// Return number of branches left in object
  inline int numberBranchesLeft() const
  {return numberBranchesLeft_;}

  /// Set number of branches left in object
  inline void setNumberBranchesLeft(int value)
  {numberBranchesLeft_ = value;}

  /// Return number of objects pointing to this
  inline int numberPointingToThis() const
  {return numberPointingToThis_;}

  /// Set number of objects pointing to this
  inline void setNumberPointingToThis(int number)
  {numberPointingToThis_=number;}

  /// Increment number of objects pointing to this
  inline void incrementNumberPointingToThis()
  {numberPointingToThis_ ++;}

  /// Say one branch taken 
  inline int branchedOn()
  {numberPointingToThis_--;numberBranchesLeft_--;return numberBranchesLeft_;}

  /// Say thrown away
  inline void throwAway()
  {numberPointingToThis_-=numberBranchesLeft_;numberBranchesLeft_=0;}

  /// Parent of this
  CbcNodeInfo * parent() const
  {return parent_;}
  /// Set parent null
  inline void nullParent()
  { parent_=NULL;}

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

  /// Increment active cut counts
  void incrementCuts(int change=1);

  /// Decrement all active cut counts in chain starting at parent
  void decrementParentCuts(CbcModel * model, int change=1);

  /// Increment all active cut counts in parent chain
  void incrementParentCuts(CbcModel * model, int change=1);

  /// Array of pointers to cuts
  inline CbcCountRowCut ** cuts() const
  {return cuts_;}

  /// Number of row cuts (this node)
  inline int numberCuts() const
  {return numberCuts_;}
  inline void setNumberCuts(int value)
  {numberCuts_=value;}

  /// Set owner null
  inline void nullOwner()
  { owner_=NULL;}
  const inline CbcNode * owner() const
  { return owner_;}
  inline CbcNode * mutableOwner() const
  { return owner_;}
  /// The node number
  inline int nodeNumber() const
  { return nodeNumber_;}
  inline void setNodeNumber(int node)
  { nodeNumber_=node;}
  /** Deactivate node information.
      1 - bounds
      2 - cuts
      4 - basis!
  */
  void deactivate(int mode=3);
  /// Say if normal
  inline bool allActivated() const
  { return (active_==7);}
  /// Say if marked
  inline bool marked() const
  { return ((active_&8)!=0);}
  /// Mark
  inline void mark()
  { active_ |= 8;}
  /// Unmark
  inline void unmark()
  { active_ &= ~8;}

  /// Branching object for the parent
  inline const OsiBranchingObject * parentBranch() const
  { return parentBranch_;}
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
  CbcNodeInfo * parent_;

  /// Copy of the branching object of the parent when the node is created
  OsiBranchingObject * parentBranch_;
      
  /// Owner
  CbcNode * owner_;

  /// Number of row cuts (this node)
  int numberCuts_;

  /// The node number
  int nodeNumber_;

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
  /** Active node information.
      1 - bounds
      2 - cuts
      4 - basis!
  */
  int active_;

private:
  
  /// Illegal Assignment operator 
  CbcNodeInfo & operator=(const CbcNodeInfo& rhs);

  /// routine common to constructors 
  void setParentBasedData();
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

  /// Just apply bounds to one variable - force means overwrite by lower,upper (1=>infeasible)
  virtual int applyBounds(int iColumn, double & lower, double & upper,int force) ;

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
  /// Lower bounds
  inline const double * lower() const
  { return lower_;}
  /// Upper bounds
  inline const double * upper() const
  { return upper_;}
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

  /// Just apply bounds to one variable - force means overwrite by lower,upper (1=>infeasible)
  virtual int applyBounds(int iColumn, double & lower, double & upper,int force) ;
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
  /// Basis diff information
  inline const CoinWarmStartDiff *basisDiff() const
  { return basisDiff_ ;}
  /// Which variable (top bit if upper bound changing)
  inline const int * variables() const
  { return variables_;}
  // New bound
  inline const double * newBounds() const
  { return newBounds_;}
  /// Number of bound changes
  inline int numberChangedBounds() const
  { return numberChangedBounds_;}
protected:
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

class CbcNode : public CoinTreeNode {
 
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
      <li> >0: Number of quich branching objects (and branches will be non NULL)
    </ul>
  */
  int chooseDynamicBranch (CbcModel * model,
                           CbcNode * lastNode,
                           OsiSolverBranch * & branches,
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
  int chooseOsiBranch (CbcModel * model,
		       CbcNode * lastNode,
		       OsiBranchingInformation * usefulInfo,
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
  int chooseClpBranch (CbcModel * model,
		       CbcNode * lastNode);
  int analyze(CbcModel * model,double * results);
  /// Decrement active cut counts
  void decrementCuts(int change=1);

  /// Decrement all active cut counts in chain starting at parent
  void decrementParentCuts(CbcModel * model, int change=1);

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
  int branch(OsiSolverInterface * solver);

  /** Double checks in case node can change its mind!
      Returns objective value
      Can change objective etc */
  double checkIsCutoff(double cutoff);
  // Information to make basis and bounds
  inline CbcNodeInfo * nodeInfo() const
  {return nodeInfo_;}

  // Objective value
  inline double objectiveValue() const
  { return objectiveValue_;}
  inline void setObjectiveValue(double value)
  { objectiveValue_=value;}
  /// Number of arms defined for the attached OsiBranchingObject.
  inline int numberBranches() const
  { if (branch_)
      return (branch_->numberBranches()) ;
    else
      return (-1) ; } 

  /* Active arm of the attached OsiBranchingObject.
  
   In the simplest instance, coded -1 for the down arm of the branch, +1 for
   the up arm. But see OsiBranchingObject::way() 
     Use nodeInfo--.numberBranchesLeft_ to see how active
  */
  int way() const;
  /// Depth in branch-and-cut search tree
  inline int depth() const
  {return depth_;}
  /// Set depth in branch-and-cut search tree
  inline void setDepth(int value)
  {depth_ = value;}
  /// Get the number of objects unsatisfied at this node.
  inline int numberUnsatisfied() const
  { return numberUnsatisfied_;}
  /// Set the number of objects unsatisfied at this node.
  inline void setNumberUnsatisfied(int value)
  { numberUnsatisfied_ = value;}
  /// Get sum of "infeasibilities" reported by each object
  inline double sumInfeasibilities() const
  { return sumInfeasibilities_;}
  /// Set sum of "infeasibilities" reported by each object
  inline void setSumInfeasibilities(double value)
  { sumInfeasibilities_ = value;}
  // Guessed objective value (for solution)
  inline double guessedObjectiveValue() const
  {return guessedObjectiveValue_;}
  inline void setGuessedObjectiveValue(double value)
  {guessedObjectiveValue_=value;}
  /// Branching object for this node
  inline const OsiBranchingObject * branchingObject() const
  { return branch_;}
  /// Modifiable branching object for this node
  inline OsiBranchingObject * modifiableBranchingObject() const
  { return branch_;}
  /// Set branching object for this node (takes ownership)
  inline void setBranchingObject(OsiBranchingObject * branchingObject)
  { branch_ = branchingObject;}
  /// The node number  
  inline int nodeNumber() const
  { return nodeNumber_;}
  inline void setNodeNumber(int node)
  { nodeNumber_=node;}
  /// Returns true if on tree
  inline bool onTree() const
  { return (state_&1)!=0;}
  /// Sets true if on tree
  inline void setOnTree(bool yesNo)
  { if(yesNo) state_ |= 1; else state_ &= ~1; }
  /// Returns true if active
  inline bool active() const
  { return (state_&2)!=0;}
  /// Sets true if active
  inline void setActive(bool yesNo)
  { if(yesNo) state_ |= 2; else state_ &= ~2; }
  /// Print
  void print() const;
  /// Debug
  inline void checkInfo() const
  { assert (nodeInfo_->numberBranchesLeft()==
	    branch_->numberBranchesLeft());}

private:
  // Data
  /// Information to make basis and bounds
  CbcNodeInfo * nodeInfo_;
  /// Objective value
  double objectiveValue_;
  /// Guessed satisfied Objective value
  double guessedObjectiveValue_;
  /// Sum of "infeasibilities" reported by each object
  double sumInfeasibilities_;
  /// Branching object for this node
  OsiBranchingObject * branch_;
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
