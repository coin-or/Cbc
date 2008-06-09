#include "CoinTime.hpp"
#include "OsiClpSolverInterface.hpp"


// below needed for pathetic branch and bound code
#include <vector>
#include <map>

// Trivial class for Branch and Bound

class DBNodeSimple  {
public:
  enum DBNodeWay {
    DOWN_UP__NOTHING_DONE=0x11,
    DOWN_UP__DOWN_DONE=0x12,
    DOWN_UP__BOTH_DONE=0x14,
    DOWN_UP__=0x10,
    UP_DOWN__NOTHING_DONE=0x21,
    UP_DOWN__UP_DONE=0x24,
    UP_DOWN__BOTH_DONE=0x22,
    UP_DOWN__=0x20,
    DOWN_CURRENT=0x02,
    UP_CURRENT=0x04,
    WAY_UNSET
  };
  
public:
    
  // Default Constructor 
  DBNodeSimple ();

  // Constructor from current state (and list of integers)
  // Also chooses branching variable (if none set to -1)
  DBNodeSimple (OsiSolverInterface &model,
                 int numberIntegers, int * integer,
                 CoinWarmStart * basis);
  void gutsOfConstructor (OsiSolverInterface &model,
                 int numberIntegers, int * integer,
                 CoinWarmStart * basis);
  // Copy constructor 
  DBNodeSimple ( const DBNodeSimple &);
   
  // Assignment operator 
  DBNodeSimple & operator=( const DBNodeSimple& rhs);

  // Destructor 
  ~DBNodeSimple ();
  // Work of destructor
  void gutsOfDestructor();
  // Extension - true if other extension of this
  bool extension(const DBNodeSimple & other,
		 const double * originalLower,
		 const double * originalUpper) const;
  inline void incrementDescendants()
  { descendants_++;}
  // Public data
  // Basis (should use tree, but not as wasteful as bounds!)
  CoinWarmStart * basis_;
  // Objective value (COIN_DBL_MAX) if spare node
  double objectiveValue_;
  // Branching variable (0 is first integer) 
  int variable_;
  // Way to branch: see enum DBNodeWay
  int way_;
  // Number of integers (for length of arrays)
  int numberIntegers_;
  // Current value
  double value_;
  // Number of descendant nodes (so 2 is in interior)
  int descendants_;
  // Parent 
  int parent_;
  // Left child
  int child_down_;
  // Right child
  int child_up_;
  // Previous in chain
  int previous_;
  // Next in chain
  int next_;
  // Now I must use tree
  // Bounds stored in full (for integers)
  int * lower_;
  int * upper_;
};


DBNodeSimple::DBNodeSimple() :
  basis_(NULL),
  objectiveValue_(COIN_DBL_MAX),
  variable_(-100),
  way_(WAY_UNSET),
  numberIntegers_(0),
  value_(0.5),
  descendants_(-1),
  parent_(-1),
  child_down_(-1),
  child_up_(-1),
  previous_(-1),
  next_(-1),
  lower_(NULL),
  upper_(NULL)
{
}
DBNodeSimple::DBNodeSimple(OsiSolverInterface & model,
		 int numberIntegers, int * integer,CoinWarmStart * basis)
{
  gutsOfConstructor(model,numberIntegers,integer,basis);
}
void
DBNodeSimple::gutsOfConstructor(OsiSolverInterface & model,
		 int numberIntegers, int * integer,CoinWarmStart * basis)
{
  basis_ = basis;
  variable_=-1;
  way_=WAY_UNSET;
  numberIntegers_=numberIntegers;
  value_=0.0;
  descendants_ = 0;
  parent_ = -1;
  child_down_ = -1;
  child_up_ = -1;
  previous_ = -1;
  next_ = -1;
  if (model.isProvenOptimal()&&!model.isDualObjectiveLimitReached()) {
    objectiveValue_ = model.getObjSense()*model.getObjValue();
  } else {
    objectiveValue_ = 1.0e100;
    lower_ = NULL;
    upper_ = NULL;
    return; // node cutoff
  }
  lower_ = new int [numberIntegers_];
  upper_ = new int [numberIntegers_];
  assert (upper_!=NULL);
  const double * lower = model.getColLower();
  const double * upper = model.getColUpper();
  const double * solution = model.getColSolution();
  int i;
  // Hard coded integer tolerance
#define INTEGER_TOLERANCE 1.0e-6
  ///////// Start of Strong branching code - can be ignored
  // Number of strong branching candidates
#define STRONG_BRANCHING 5
#ifdef STRONG_BRANCHING
  double upMovement[STRONG_BRANCHING];
  double downMovement[STRONG_BRANCHING];
  double solutionValue[STRONG_BRANCHING];
  int chosen[STRONG_BRANCHING];
  int iSmallest=0;
  // initialize distance from integer
  for (i=0;i<STRONG_BRANCHING;i++) {
    upMovement[i]=0.0;
    chosen[i]=-1;
  }
  variable_=-1;
  // This has hard coded integer tolerance
  double mostAway=INTEGER_TOLERANCE;
  int numberAway=0;
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integer[i];
    lower_[i]=(int)lower[iColumn];
    upper_[i]=(int)upper[iColumn];
    double value = solution[iColumn];
    value = max(value,(double) lower_[i]);
    value = min(value,(double) upper_[i]);
    double nearest = floor(value+0.5);
    if (fabs(value-nearest)>INTEGER_TOLERANCE)
      numberAway++;
    if (fabs(value-nearest)>mostAway) {
      double away = fabs(value-nearest);
      if (away>upMovement[iSmallest]) {
	//add to list
	upMovement[iSmallest]=away;
	solutionValue[iSmallest]=value;
	chosen[iSmallest]=i;
	int j;
	iSmallest=-1;
	double smallest = 1.0;
	for (j=0;j<STRONG_BRANCHING;j++) {
	  if (upMovement[j]<smallest) {
	    smallest=upMovement[j];
	    iSmallest=j;
	  }
	}
      }
    }
  }
  int numberStrong=0;
  for (i=0;i<STRONG_BRANCHING;i++) {
    if (chosen[i]>=0) { 
      numberStrong ++;
      variable_ = chosen[i];
    }
  }
  // out strong branching if bit set
  OsiClpSolverInterface* clp =
    dynamic_cast<OsiClpSolverInterface*>(&model);
  if (clp&&(clp->specialOptions()&16)!=0&&numberStrong>1) {
    int j;
    int iBest=-1;
    double best = 0.0;
    for (j=0;j<STRONG_BRANCHING;j++) {
      if (upMovement[j]>best) {
        best=upMovement[j];
        iBest=j;
      }
    }
    numberStrong=1;
    variable_=chosen[iBest];
  }
  if (numberStrong==1) {
    // just one - makes it easy
    int iColumn = integer[variable_];
    double value = solution[iColumn];
    value = max(value,(double) lower_[variable_]);
    value = min(value,(double) upper_[variable_]);
    double nearest = floor(value+0.5);
    value_=value;
    if (value<=nearest)
      way_=UP_DOWN__NOTHING_DONE; // up
    else
      way_=DOWN_UP__NOTHING_DONE; // down
  } else if (numberStrong) {
    // more than one - choose
    bool chooseOne=true;
    model.markHotStart();
    for (i=0;i<STRONG_BRANCHING;i++) {
      int iInt = chosen[i];
      if (iInt>=0) {
	int iColumn = integer[iInt];
	double value = solutionValue[i]; // value of variable in original
	double objectiveChange;
	value = max(value,(double) lower_[iInt]);
	value = min(value,(double) upper_[iInt]);

	// try down

	model.setColUpper(iColumn,floor(value));
	model.solveFromHotStart();
	model.setColUpper(iColumn,upper_[iInt]);
	if (model.isProvenOptimal()&&!model.isDualObjectiveLimitReached()) {
	  objectiveChange = model.getObjSense()*model.getObjValue()
	    - objectiveValue_;
	} else {
	  objectiveChange = 1.0e100;
	}
	assert (objectiveChange>-1.0e-5);
	objectiveChange = CoinMax(objectiveChange,0.0);
	downMovement[i]=objectiveChange;

	// try up

	model.setColLower(iColumn,ceil(value));
	model.solveFromHotStart();
	model.setColLower(iColumn,lower_[iInt]);
	if (model.isProvenOptimal()&&!model.isDualObjectiveLimitReached()) {
	  objectiveChange = model.getObjSense()*model.getObjValue()
	    - objectiveValue_;
	} else {
	  objectiveChange = 1.0e100;
	}
	assert (objectiveChange>-1.0e-5);
	objectiveChange = CoinMax(objectiveChange,0.0);
	upMovement[i]=objectiveChange;
	
	/* Possibilities are:
	   Both sides feasible - store
	   Neither side feasible - set objective high and exit
	   One side feasible - change bounds and resolve
	*/
	bool solveAgain=false;
	if (upMovement[i]<1.0e100) {
	  if(downMovement[i]<1.0e100) {
	    // feasible - no action
	  } else {
	    // up feasible, down infeasible
	    solveAgain = true;
	    model.setColLower(iColumn,ceil(value));
	  }
	} else {
	  if(downMovement[i]<1.0e100) {
	    // down feasible, up infeasible
	    solveAgain = true;
	    model.setColUpper(iColumn,floor(value));
	  } else {
	    // neither side feasible
	    objectiveValue_=1.0e100;
	    chooseOne=false;
	    break;
	  }
	}
	if (solveAgain) {
	  // need to solve problem again - signal this
	  variable_ = numberIntegers;
	  chooseOne=false;
	  break;
	}
      }
    }
    if (chooseOne) {
      // choose the one that makes most difference both ways
      double best = -1.0;
      double best2 = -1.0;
      for (i=0;i<STRONG_BRANCHING;i++) {
	int iInt = chosen[i];
	if (iInt>=0) {
	  //std::cout<<"Strong branching on "
          //   <<i<<""<<iInt<<" down "<<downMovement[i]
          //   <<" up "<<upMovement[i]
          //   <<" value "<<solutionValue[i]
          //   <<std::endl;
	  bool better = false;
	  if (min(upMovement[i],downMovement[i])>best) {
	    // smaller is better
	    better=true;
	  } else if (min(upMovement[i],downMovement[i])>best-1.0e-5) {
	    if (max(upMovement[i],downMovement[i])>best2+1.0e-5) {
	      // smaller is about same, but larger is better
	      better=true;
	    }
	  }
	  if (better) {
	    best = min(upMovement[i],downMovement[i]);
	    best2 = max(upMovement[i],downMovement[i]);
	    variable_ = iInt;
	    double value = solutionValue[i];
	    value = max(value,(double) lower_[variable_]);
	    value = min(value,(double) upper_[variable_]);
	    value_=value;
	    if (upMovement[i]<=downMovement[i])
	      way_=UP_DOWN__NOTHING_DONE; // up
	    else
	      way_=DOWN_UP__NOTHING_DONE; // down
	  }
	}
      }
    }
    // Delete the snapshot
    model.unmarkHotStart();
  }
  ////// End of Strong branching
#else
  variable_=-1;
  // This has hard coded integer tolerance
  double mostAway=INTEGER_TOLERANCE;
  int numberAway=0;
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integer[i];
    lower_[i]=(int)lower[iColumn];
    upper_[i]=(int)upper[iColumn];
    double value = solution[iColumn];
    value = max(value,(double) lower_[i]);
    value = min(value,(double) upper_[i]);
    double nearest = floor(value+0.5);
    if (fabs(value-nearest)>INTEGER_TOLERANCE)
      numberAway++;
    if (fabs(value-nearest)>mostAway) {
      mostAway=fabs(value-nearest);
      variable_=i;
      value_=value;
      if (value<=nearest)
	way_=UP_DOWN__NOTHING_DONE; // up
      else
	way_=DOWN_UP__NOTHING_DONE; // down
    }
  }
#endif
}

DBNodeSimple::DBNodeSimple(const DBNodeSimple & rhs) 
{
  if (rhs.basis_)
    basis_=rhs.basis_->clone();
  else
    basis_ = NULL;
  objectiveValue_=rhs.objectiveValue_;
  variable_=rhs.variable_;
  way_=rhs.way_;
  numberIntegers_=rhs.numberIntegers_;
  value_=rhs.value_;
  descendants_ = rhs.descendants_;
  parent_ = rhs.parent_;
  child_down_ = rhs.child_down_;
  child_up_ = rhs.child_up_;
  previous_ = rhs.previous_;
  next_ = rhs.next_;
  lower_=NULL;
  upper_=NULL;
  if (rhs.lower_!=NULL) {
    lower_ = new int [numberIntegers_];
    upper_ = new int [numberIntegers_];
    assert (upper_!=NULL);
    memcpy(lower_,rhs.lower_,numberIntegers_*sizeof(int));
    memcpy(upper_,rhs.upper_,numberIntegers_*sizeof(int));
  }
}

DBNodeSimple &
DBNodeSimple::operator=(const DBNodeSimple & rhs)
{
  if (this != &rhs) {
    gutsOfDestructor();
    if (rhs.basis_)
      basis_=rhs.basis_->clone();
    objectiveValue_=rhs.objectiveValue_;
    variable_=rhs.variable_;
    way_=rhs.way_;
    numberIntegers_=rhs.numberIntegers_;
    value_=rhs.value_;
    descendants_ = rhs.descendants_;
    parent_ = rhs.parent_;
    child_down_ = rhs.child_down_;
    child_up_ = rhs.child_up_;
    previous_ = rhs.previous_;
    next_ = rhs.next_;
    if (rhs.lower_!=NULL) {
      lower_ = new int [numberIntegers_];
      upper_ = new int [numberIntegers_];
      assert (upper_!=NULL);
      memcpy(lower_,rhs.lower_,numberIntegers_*sizeof(int));
      memcpy(upper_,rhs.upper_,numberIntegers_*sizeof(int));
    }
  }
  return *this;
}


DBNodeSimple::~DBNodeSimple ()
{
  gutsOfDestructor();
}
// Work of destructor
void 
DBNodeSimple::gutsOfDestructor()
{
  delete [] lower_;
  delete [] upper_;
  delete basis_;
  lower_ = NULL;
  upper_ = NULL;
  basis_ = NULL;
  objectiveValue_ = COIN_DBL_MAX;
}
// Extension - true if other extension of this
bool 
DBNodeSimple::extension(const DBNodeSimple & other,
			 const double * originalLower,
			 const double * originalUpper) const
{
  bool ok=true;
  for (int i=0;i<numberIntegers_;i++) {
    if (upper_[i]<originalUpper[i]||
	lower_[i]>originalLower[i]) {
      if (other.upper_[i]>upper_[i]||
	  other.lower_[i]<lower_[i]) {
	ok=false;
	break;
      }
    }
  }
  return ok;
}

#include <vector>
#define FUNNY_BRANCHING 1
#define FUNNY_TREE
#ifndef FUNNY_TREE
// Vector of DBNodeSimples 
typedef std::vector<DBNodeSimple>    DBVectorNode;
#else
// Must code up by hand
class DBVectorNode  {
  
public:
    
  // Default Constructor 
  DBVectorNode ();

  // Copy constructor 
  DBVectorNode ( const DBVectorNode &);
   
  // Assignment operator 
  DBVectorNode & operator=( const DBVectorNode& rhs);

  // Destructor 
  ~DBVectorNode ();
  // Size
  inline int size() const
  { return size_-sizeDeferred_;}
  // Push
  void push_back(const DBNodeSimple & node);
  // Last one in (or other criterion)
  DBNodeSimple back() const;
  // Get rid of last one
  void pop_back();
  // Works out best one
  int best() const;
  
  // Public data
  // Maximum size
  int maximumSize_;
  // Current size
  int size_;
  // Number still hanging around
  int sizeDeferred_;
  // First spare
  int firstSpare_;
  // First 
  int first_;
  // Last 
  int last_;
  // Chosen one
  mutable int chosen_;
  // Nodes
  DBNodeSimple * nodes_;
};


DBVectorNode::DBVectorNode() :
  maximumSize_(10),
  size_(0),
  sizeDeferred_(0),
  firstSpare_(0),
  first_(-1),
  last_(-1)
{
  nodes_ = new DBNodeSimple[maximumSize_];
  for (int i=0;i<maximumSize_;i++) {
    nodes_[i].previous_=i-1;
    nodes_[i].next_=i+1;
  }
}

DBVectorNode::DBVectorNode(const DBVectorNode & rhs) 
{  
  maximumSize_ = rhs.maximumSize_;
  size_ = rhs.size_;
  sizeDeferred_ = rhs.sizeDeferred_;
  firstSpare_ = rhs.firstSpare_;
  first_ = rhs.first_;
  last_ = rhs.last_;
  nodes_ = new DBNodeSimple[maximumSize_];
  for (int i=0;i<maximumSize_;i++) {
    nodes_[i] = rhs.nodes_[i];
  }
}

DBVectorNode &
DBVectorNode::operator=(const DBVectorNode & rhs)
{
  if (this != &rhs) {
    delete [] nodes_;
    maximumSize_ = rhs.maximumSize_;
    size_ = rhs.size_;
    sizeDeferred_ = rhs.sizeDeferred_;
    firstSpare_ = rhs.firstSpare_;
    first_ = rhs.first_;
    last_ = rhs.last_;
    nodes_ = new DBNodeSimple[maximumSize_];
    for (int i=0;i<maximumSize_;i++) {
      nodes_[i] = rhs.nodes_[i];
    }
  }
  return *this;
}


DBVectorNode::~DBVectorNode ()
{
  delete [] nodes_;
}
// Push
void 
DBVectorNode::push_back(const DBNodeSimple & node)
{
  if (size_==maximumSize_) {
    assert (firstSpare_==size_);
    maximumSize_ = (maximumSize_*3)+10;
    DBNodeSimple * temp = new DBNodeSimple[maximumSize_];
    int i;
    for (i=0;i<size_;i++) {
      temp[i]=nodes_[i];
    }
    delete [] nodes_;
    nodes_ = temp;
    //firstSpare_=size_;
    int last = -1;
    for ( i=size_;i<maximumSize_;i++) {
      nodes_[i].previous_=last;
      nodes_[i].next_=i+1;
      last = i;
    }
  }
  assert (firstSpare_<maximumSize_);
  assert (nodes_[firstSpare_].previous_<0);
  int next = nodes_[firstSpare_].next_;
  nodes_[firstSpare_]=node;
  if (last_>=0) {
    assert (nodes_[last_].next_==-1);
    nodes_[last_].next_=firstSpare_;
  }
  nodes_[firstSpare_].previous_=last_;
  nodes_[firstSpare_].next_=-1;
  if (last_==-1) {
    assert (first_==-1);
    first_ = firstSpare_;
  }
  last_=firstSpare_;
  if (next>=0&&next<maximumSize_) {
    firstSpare_ = next;
    nodes_[firstSpare_].previous_=-1;
  } else {
    firstSpare_=maximumSize_;
  }
  chosen_ = -1;
  //best();
  size_++;
  assert (node.descendants_<=2);
  if (node.descendants_==2)
    sizeDeferred_++;
}
// Works out best one
int 
DBVectorNode::best() const
{
  // can modify
  chosen_=-1;
  if (chosen_<0) {
    chosen_=last_;
#if FUNNY_BRANCHING
    while (nodes_[chosen_].descendants_==2) {
      chosen_ = nodes_[chosen_].previous_;
      assert (chosen_>=0);
    }
#endif
  }
  return chosen_;
}
// Last one in (or other criterion)
DBNodeSimple 
DBVectorNode::back() const
{
  assert (last_>=0);
  return nodes_[best()];
}
// Get rid of last one
void 
DBVectorNode::pop_back()
{
  // Temporary until more sophisticated
  //assert (last_==chosen_);
  if (nodes_[chosen_].descendants_==2)
    sizeDeferred_--;
  int previous = nodes_[chosen_].previous_;
  int next = nodes_[chosen_].next_;
  nodes_[chosen_].gutsOfDestructor();
  if (previous>=0) {
    nodes_[previous].next_=next;
  } else {
    first_ = next;
  }
  if (next>=0) {
    nodes_[next].previous_ = previous;
  } else {
    last_ = previous;
  }
  nodes_[chosen_].previous_=-1;
  if (firstSpare_>=0) {
    nodes_[chosen_].next_ = firstSpare_;
  } else {
    nodes_[chosen_].next_ = -1;
  }
  firstSpare_ = chosen_;
  chosen_ = -1;
  assert (size_>0);
  size_--;
}
#endif

bool
DBNodeSimple::isGrandparentIrrelevant()
{
#if !defined(FUNNY_BRANCHING)
  return false;
#endif

  if (parent_ == -1) {
    // can't flip root higher up...
    return false;
  }
  
  if (model.isProvenDualInfeasible()) {
    // THINK: this is not going to happen, but if it does, what should we do???
    return false;
  }
  if (model.isProvenPrimalInfeasible()) {
    ...;
  }
  // Both primal and dual feasible, and in this case we don't care how we have
  // stopped (iteration limit, obj val limit, time limit, optimal solution,
  // etc.), we can just look at the reduced costs to figure out if the
  // grandparent is irrelevant. Remember, we are using dual simplex!
  const DBNodeSimple& parent = branchingTree.nodes_[parent_];
  const int iColumn = which[parent.variable_];
  double djValue = model.getReducedCost()[iColumn]*direction;
  const bool down_child = branchingTree.nodeIndex(node) == parent.child_down_;
  if (djValue>1.0e-6) {
    // wants to go down
    if (down_child) {
      return true;
    }
    const double up_lower = std::floor(parent.value_);
    if (model.getColLower()[iColumn] > up_lower) {
      return true;
    }
    return false;
  } else {
    // wants to go up
    if (!down_child) {
      return true;
    }
    const double down_upper = std::ceil(parent.value_);
    if (model.getColUpper()[iColumn] < down_upper) {
      return true;
    }
    return false;
  }
  return false;
}

void
DBVectorNode::moveNodeUp()
{
  assert(parent != -1);
  const int parent_id = node.parent_;
  const DBNodeSimple& parent = branchingTree.nodes_[parent_id];
  const int node_id = branchingTree.nodeIndex(node);
  const bool node_is_down_child = node_id == parent.child_down_;
  const int grandparent_id = parent.parent_;

  // First hang the nodes where they belong.
  node.parent_ = grandparent_id;
  parent.parent_ = node_id;
#if 1
  int& child_to_move = (node.way_ & DOWN_CURRENT) ? node.child_up_ : node.child_down_;
  if (node_is_down_child) {
    parent.child_down_ = child_to_move;
  } else {
    parent.child_up_ = child_to_move;
  }
  if (child_to_move >= 0) {
    branchingTree.nodes_[child_to_move].parent_ = parent_id;
  }
  child_to_move = parent_id;
#else
  if (node.way_ & DOWN_CURRENT) {
    if (node_is_down_child) {
      parent.child_down_ = node.child_up_;
    } else {
      parent.child_up_ = node.child_up_;
    }
    if (node.child_up_ >= 0) {
      branchingTree.nodes_[node.child_up_].parent_ = parent_id;
    }
    node.child_up_ = parent_id;
  } else { // must be UP_CURRENT
    if (node_is_down_child) {
      parent.child_down_ = node.child_down_;
    } else {
      parent.child_up_ = node.child_down_;
    }
    if (node.child_down_ >= 0) {
      branchingTree.nodes_[node.child_down_].parent_ = parent_id;
    }
    node.child_down_ = parent_id;
  }
#endif
  if (grandparent_id >= 0) {
    if (parent_id == branchingTree.nodes_[grandparent_id].child_down_) {
      branchingTree.nodes_[grandparent_id].child_down_ = node_id;
    } else {
      branchingTree.nodes_[grandparent_id].child_up_ = node_id;
    }
  }

  // Now modify bounds

  // THINK: could be avoided if always start from original bounds and go back
  // to root to apply all branching decisions. On the other hand, reduced cost
  // fixing would be lost. And so would fixings by strong branching. Actually,
  // what we'd ideally need to do is to apply flipping when strong branching
  // fixes a variable.
}





bool moveNodes(OsiSolverInterface & model,
	       DBVectorNode & branchingTree,
	       int kNode)
{

  DBNodeSimple & node = branchingTree[kNode];
  DBNodeSimple & grandParent = branchingTree[node.parent_];
  int grandParentVariable = grandParent.variable_;
  // check if branching constraint of grandparent is tight
  bool canMoveNodes = checkGrandparent();

  if(!canMoveNodes)
    return false;

  node.parent_ = grandParent.parent_;
  grandParent.parent_ = kNode;
  // change bounds of grandParent
  grandParent.lower_[


}

// Invoke solver's built-in enumeration algorithm
void 
branchAndBound(OsiSolverInterface & model) {
  double time1 = CoinCpuTime();
  // solve LP
  model.initialSolve();
  int funnyBranching=FUNNY_BRANCHING;

  if (model.isProvenOptimal()&&!model.isDualObjectiveLimitReached()) {
    // Continuous is feasible - find integers
    int numberIntegers=0;
    int numberColumns = model.getNumCols();
    int iColumn;
    int i;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if( model.isInteger(iColumn))
        numberIntegers++;
    }
    if (!numberIntegers) {
      std::cout<<"No integer variables"
               <<std::endl;
      return;
    }
    int * which = new int[numberIntegers]; // which variables are integer
    // original bounds
    int * originalLower = new int[numberIntegers];
    int * originalUpper = new int[numberIntegers];
    int * relaxedLower = new int[numberIntegers];
    int * relaxedUpper = new int[numberIntegers];
    {
      const double * lower = model.getColLower();
      const double * upper = model.getColUpper();
      numberIntegers=0;
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	if( model.isInteger(iColumn)) {
	  originalLower[numberIntegers]=(int) lower[iColumn];
	  originalUpper[numberIntegers]=(int) upper[iColumn];
	  which[numberIntegers++]=iColumn;
	}
      }
    }
    double direction = model.getObjSense();
    // empty tree
    DBVectorNode branchingTree;
    
    // Add continuous to it;
    DBNodeSimple rootNode(model,numberIntegers,which,model.getWarmStart());
    // something extra may have been fixed by strong branching
    // if so go round again
    while (rootNode.variable_==numberIntegers) {
      model.resolve();
      rootNode = DBNodeSimple(model,numberIntegers,which,model.getWarmStart());
    }
    if (rootNode.objectiveValue_<1.0e100) {
      // push on stack
      branchingTree.push_back(rootNode);
    }
    
    // For printing totals
    int numberIterations=0;
    int numberNodes =0;
    int nRedundantUp=0;
    int nRedundantDown=0;
    int nRedundantUp2=0;
    int nRedundantDown2=0;
    DBNodeSimple bestNode;
    ////// Start main while of branch and bound
    // while until nothing on stack
    while (branchingTree.size()) {
      // last node
      DBNodeSimple node = branchingTree.back();
      int kNode = branchingTree.chosen_;
      branchingTree.pop_back();
      assert (node.descendants_<2);
      numberNodes++;
      if (node.variable_>=0) {
        // branch - do bounds
        for (i=0;i<numberIntegers;i++) {
          iColumn=which[i];
          model.setColBounds( iColumn,node.lower_[i],node.upper_[i]);
        }
        // move basis
        model.setWarmStart(node.basis_);
        // do branching variable
	node.incrementDescendants();
	bool down_branch = true;
	switch (node.way_) {
	case WAY_UNSET:
	case DOWN_UP__BOTH_DONE:
	case UP_DOWN__BOTH_DONE:
	  abort();
	case DOWN_UP__NOTHING_DONE:
	  node.way_ = DOWN_UP__DOWN_DONE;
	  break;
	case DOWN_UP__DOWN_DONE:
	  node.way_ = DOWN_UP__BOTH_DONE:
	  down_branch = false;
	  break;
	case UP_DOWN__NOTHING_DONE:
	  node.way_ = UP_DOWN__UP_DONE:
	  down_branch = false;
	  break;
	case UP_DOWN__UP_DONE:
	  node.way_ = UP_DOWN__BOTH_DONE;
	  break;
	}
        if (down_branch) {
          model.setColUpper(which[node.variable_],floor(node.value_));
        } else {
          model.setColLower(which[node.variable_],ceil(node.value_));
        }
	// put back on tree anyway regardless whether any processing is left
	// to be done. We want the whole tree all the time.
	branchingTree.push_back(node);
	
        // solve
        model.resolve();
        CoinWarmStart * ws = model.getWarmStart();
        const CoinWarmStartBasis* wsb =
          dynamic_cast<const CoinWarmStartBasis*>(ws);
        assert (wsb!=NULL); // make sure not volume
        numberIterations += model.getIterationCount();
        // fix on reduced costs
        int nFixed0=0,nFixed1=0;
        double cutoff;
        model.getDblParam(OsiDualObjectiveLimit,cutoff);
        double gap=(cutoff-model.getObjValue())*direction+1.0e-4;
	//        double gap=(cutoff-modelPtr_->objectiveValue())*direction+1.0e-4;
        if (gap<1.0e10&&model.isProvenOptimal()&&!model.isDualObjectiveLimitReached()) {
          const double * dj = model.getReducedCost();
          const double * lower = model.getColLower();
          const double * upper = model.getColUpper();
          for (i=0;i<numberIntegers;i++) {
            iColumn=which[i];
            if (upper[iColumn]>lower[iColumn]) {
              double djValue = dj[iColumn]*direction;
              if (wsb->getStructStatus(iColumn)==CoinWarmStartBasis::atLowerBound&&
                  djValue>gap) {
                nFixed0++;
                model.setColUpper(iColumn,lower[iColumn]);
              } else if (wsb->getStructStatus(iColumn)==CoinWarmStartBasis::atUpperBound&&
                         -djValue>gap) {
                nFixed1++;
                model.setColLower(iColumn,upper[iColumn]);
              }
            }
          }
          //if (nFixed0+nFixed1)
          //printf("%d fixed to lower, %d fixed to upper\n",nFixed0,nFixed1);
        }
	if (model.isAbandoned()) {
	  // THINK: What the heck should we do???
	  abort();
	}
	if (node.isGrandparentIrrelevant()) {
	  branchingTree.moveNodeUp(node);
	}




	if (!model.isIterationLimitReached()) {
	  if (model.isProvenOptimal()&&!model.isDualObjectiveLimitReached()) {
#if FUNNY_BRANCHING
	    // See if branched variable off bounds
	    const double * dj = model.getReducedCost();
	    const double * lower = model.getColLower();
	    const double * upper = model.getColUpper();
	    const double * solution = model.getColSolution();
	    // Better to use "natural" value - need flag to say fixed
	    for (i=0;i<numberIntegers;i++) {
	      iColumn=which[i];
	      relaxedLower[i]=originalLower[i];
	      relaxedUpper[i]=originalUpper[i];
	      double djValue = dj[iColumn]*direction;
	      if (djValue>1.0e-6) {
		// wants to go down
		if (lower[iColumn]>originalLower[i]) {
		  // Lower bound active
		  relaxedLower[i]=(int) lower[iColumn];
		}
		if (upper[iColumn]<originalUpper[i]) {
		  // Upper bound NOT active
		}
	      } else if (djValue<-1.0e-6) {
		// wants to go up
		if (lower[iColumn]>originalLower[i]) {
		  // Lower bound NOT active
		}
		if (upper[iColumn]<originalUpper[i]) {
		  // Upper bound active
		  relaxedUpper[i]=(int) upper[iColumn];
		}
	      }
	    }
	    // See if can do anything
	    {
	      /*
		If kNode is on second branch then
		a) If other feasible could free up as well
		b) If other infeasible could do something clever.
		For now - we have to give up
	      */
	      int jNode=branchingTree.nodes_[kNode].parent_;
	      bool canDelete = (branchingTree.nodes_[kNode].descendants_<2);
	      while (jNode>=0) {
		DBNodeSimple & node = branchingTree.nodes_[jNode];
		int next = node.parent_;
		if (node.descendants_<2) {
		  int variable = node.variable_;
		  iColumn=which[variable];
		  double value = node.value_;
		  double djValue = dj[iColumn]*direction;
		  assert (node.way_==2||node.way_==-2);
		  // we don't know which branch it was - look at current bounds
		  if (upper[iColumn]<value&&node.lower_[variable]<upper[iColumn]) {
		    // must have been down branch
		    if (djValue>1.0e-3||solution[iColumn]<upper[iColumn]-1.0e-5) {
		      if (canDelete) {
			nRedundantDown++;
#if 1
			printf("%d redundant branch down with value %g current upper %g solution %g dj %g\n",
			       variable,node.value_,upper[iColumn],solution[iColumn],djValue);
#endif
			node.descendants_=2; // ignore
			branchingTree.sizeDeferred_++;
			int newUpper = originalUpper[variable];
			if (next>=0) {
			  DBNodeSimple & node2 = branchingTree.nodes_[next];
			  newUpper = node2.upper_[variable];
			}
			if (branchingTree.nodes_[jNode].parent_!=next)
			  assert (newUpper>upper[iColumn]);
			model.setColUpper(iColumn,newUpper);
			int kNode2=next;
			int jNode2=branchingTree.nodes_[kNode].parent_;
			assert (newUpper>branchingTree.nodes_[kNode].upper_[variable]);
			branchingTree.nodes_[kNode].upper_[variable]= newUpper;
			while (jNode2!=kNode2) {
			  DBNodeSimple & node2 = branchingTree.nodes_[jNode2];
			  int next = node2.parent_;
			  if (next!=kNode2)
			    assert (newUpper>node2.upper_[variable]);
			  node2.upper_[variable]= newUpper;
			  jNode2=next;
			}
		      } else {
			// can't delete but can add other way to jNode
			nRedundantDown2++;
			DBNodeSimple & node2 = branchingTree.nodes_[kNode];
			assert (node2.way_==2||node2.way_==-2);
			double value2 = node2.value_;
			int variable2 = node2.variable_;
			int iColumn2 = which[variable2];
			if (variable != variable2) {
			  if (node2.way_==2&&upper[iColumn2]<value2) {
			    // must have been down branch which was done - carry over
			    int newUpper = (int) floor(value2);
			    assert (newUpper<node.upper_[variable2]);
			    node.upper_[variable2]=newUpper;
			  } else if (node2.way_==-2&&lower[iColumn2]>value2) {
			    // must have been up branch which was done - carry over
			    int newLower = (int) ceil(value2);
			    assert (newLower>node.lower_[variable2]);
			    node.lower_[variable2]=newLower;
			  }
			  if (node.lower_[variable2]>node.upper_[variable2]) {
			    // infeasible
			    node.descendants_=2; // ignore
			    branchingTree.sizeDeferred_++;
			  }
			}
		      }
		      break;
		    } 		
		    // we don't know which branch it was - look at current bounds
		  } else if (lower[iColumn]>value&&node.upper_[variable]>lower[iColumn]) {
		    // must have been up branch
		    if (djValue<-1.0e-3||solution[iColumn]>lower[iColumn]+1.0e-5) {
		      if (canDelete) {
			nRedundantUp++;
#if 1
			printf("%d redundant branch up with value %g current lower %g solution %g dj %g\n",
			       variable,node.value_,lower[iColumn],solution[iColumn],djValue);
#endif
			node.descendants_=2; // ignore
			branchingTree.sizeDeferred_++;
			int newLower = originalLower[variable];
			if (next>=0) {
			  DBNodeSimple & node2 = branchingTree.nodes_[next];
			  newLower = node2.lower_[variable];
			}
			if (branchingTree.nodes_[jNode].parent_!=next)
			  assert (newLower<lower[iColumn]);
			model.setColLower(iColumn,newLower);
			int kNode2=next;
			int jNode2=branchingTree.nodes_[kNode].parent_;
			assert (newLower<branchingTree.nodes_[kNode].lower_[variable]);
			branchingTree.nodes_[kNode].lower_[variable]= newLower;
			while (jNode2!=kNode2) {
			  DBNodeSimple & node2 = branchingTree.nodes_[jNode2];
			  int next = node2.parent_;
			  if (next!=kNode2)
			    assert (newLower<node2.lower_[variable]);
			  node2.lower_[variable]=newLower;
			  jNode2=next;
			}
		      } else {
			// can't delete but can add other way to jNode
			nRedundantUp2++;
			DBNodeSimple & node2 = branchingTree.nodes_[kNode];
			assert (node2.way_==2||node2.way_==-2);
			double value2 = node2.value_;
			int variable2 = node2.variable_;
			int iColumn2 = which[variable2];
			if (variable != variable2) {
			  if (node2.way_==2&&upper[iColumn2]<value2) {
			    // must have been down branch which was done - carry over
			    int newUpper = (int) floor(value2);
			    assert (newUpper<node.upper_[variable2]);
			    node.upper_[variable2]=newUpper;
			  } else if (node2.way_==-2&&lower[iColumn2]>value2) {
			    // must have been up branch which was done - carry over
			    int newLower = (int) ceil(value2);
			    assert (newLower>node.lower_[variable2]);
			    node.lower_[variable2]=newLower;
			  }
			  if (node.lower_[variable2]>node.upper_[variable2]) {
			    // infeasible
			    node.descendants_=2; // ignore
			    branchingTree.sizeDeferred_++;
			  }
			}
		      }
		      break;
		    } 		
		  }
		} else {
		  break;
		}
		jNode=next;
	      }
	    }
	    // solve
	    //resolve();
	    //assert(!getIterationCount());
	    if ((numberNodes%1000)==0) 
	      printf("%d nodes, redundant down %d (%d) up %d (%d) tree size %d\n",
		     numberNodes,nRedundantDown,nRedundantDown2,nRedundantUp,nRedundantUp2,branchingTree.size());
#else
	    if ((numberNodes%1000)==0) 
	      printf("%d nodes, tree size %d\n",
		     numberNodes,branchingTree.size());
#endif
	    if (CoinCpuTime()-time1>3600.0) {
	      printf("stopping after 3600 seconds\n");
	      exit(77);
	    }
	    DBNodeSimple newNode(model,numberIntegers,which,ws);
	    // something extra may have been fixed by strong branching
	    // if so go round again
	    while (newNode.variable_==numberIntegers) {
	      model.resolve();
	      newNode = DBNodeSimple(model,numberIntegers,which,model.getWarmStart());
	    }
	    if (newNode.objectiveValue_<1.0e100) {
	      if (newNode.variable_>=0) 
		assert (fabs(newNode.value_-floor(newNode.value_+0.5))>1.0e-6);
	      newNode.parent_ = kNode;
	      // push on stack
	      branchingTree.push_back(newNode);
	      if(branchingTree.nodes_[kNode].child_down_ < 0)
		branchingTree.nodes_[kNode].child_down_ = branchingTree.last_;
	      else
		branchingTree.nodes_[kNode].child_up_ = branchingTree.last_;
#if 0
	      } else {
		// integer solution - save
		bestNode = node;
		// set cutoff (hard coded tolerance)
		model.setDblParam(OsiDualObjectiveLimit,(bestNode.objectiveValue_-1.0e-5)*direction);
		std::cout<<"Integer solution of "
			 <<bestNode.objectiveValue_
			 <<" found after "<<numberIterations
			 <<" iterations and "<<numberNodes<<" nodes"
                 <<std::endl;
	      }
#endif
	    }
	  }
        } else {
          // maximum iterations - exit
          std::cout<<"Exiting on maximum iterations"
                   <<std::endl;
	  break;
        }
      } else {
        // integer solution - save
        bestNode = node;
        // set cutoff (hard coded tolerance)
        model.setDblParam(OsiDualObjectiveLimit,(bestNode.objectiveValue_-1.0e-5)*direction);
        std::cout<<"Integer solution of "
                 <<bestNode.objectiveValue_
                 <<" found after "<<numberIterations
                 <<" iterations and "<<numberNodes<<" nodes"
                 <<std::endl;
      }
    }
    ////// End main while of branch and bound
    std::cout<<"Search took "
             <<numberIterations
             <<" iterations and "<<numberNodes<<" nodes"
             <<std::endl;
    if (bestNode.numberIntegers_) {
      // we have a solution restore
      // do bounds
      for (i=0;i<numberIntegers;i++) {
        iColumn=which[i];
        model.setColBounds( iColumn,bestNode.lower_[i],bestNode.upper_[i]);
      }
      // move basis
      model.setWarmStart(bestNode.basis_);
      // set cutoff so will be good (hard coded tolerance)
      model.setDblParam(OsiDualObjectiveLimit,(bestNode.objectiveValue_+1.0e-5)*direction);
      model.resolve();
    } else {
      OsiClpSolverInterface* clp =
	dynamic_cast<OsiClpSolverInterface*>(&model);
      if (clp) {
	ClpSimplex* modelPtr_ = clp->getModelPtr();
	modelPtr_->setProblemStatus(1);
      }
    }
    delete [] which;
    delete [] originalLower;
    delete [] originalUpper;
    delete [] relaxedLower;
    delete [] relaxedUpper;
  } else {
    std::cout<<"The LP relaxation is infeasible"
             <<std::endl;
    OsiClpSolverInterface* clp =
      dynamic_cast<OsiClpSolverInterface*>(&model);
    if (clp) {
      ClpSimplex* modelPtr_ = clp->getModelPtr();
      modelPtr_->setProblemStatus(1);
    }
    //throw CoinError("The LP relaxation is infeasible or too expensive",
    //"branchAndBound", "OsiClpSolverInterface");
  }
}



int main(int argc, char* argv[])
{
  return 0;
}
