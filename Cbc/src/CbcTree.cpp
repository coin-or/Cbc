// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CbcModel.hpp"
#include "CbcNode.hpp"
#include "CbcTree.hpp"
#include "CbcCountRowCut.hpp"

CbcTree::CbcTree()
{
}
CbcTree::~CbcTree()
{
}
// Copy constructor 
CbcTree::CbcTree ( const CbcTree & rhs)
{
  nodes_=rhs.nodes_;
}
// Assignment operator 
CbcTree &
CbcTree::operator=(const CbcTree & rhs)
{
  if (this != &rhs) {
    nodes_=rhs.nodes_;
  }
  return *this;
}
// Clone
CbcTree *
CbcTree::clone() const
{
  return new CbcTree(*this);
}
// Set comparison function and resort heap
void 
CbcTree::setComparison(CbcCompareBase  &compare)
{
  comparison_.test_ = &compare;
  make_heap(nodes_.begin(), nodes_.end(), comparison_);
}

// Return the top node of the heap 
CbcNode * 
CbcTree::top() const
{
  return nodes_.front();
}

// Add a node to the heap
void 
CbcTree::push(CbcNode * x) {
  /*printf("push obj %g, refcount %d, left %d, pointing to %d\n",
	 x->objectiveValue(),x->nodeInfo()->decrement(0),
	 x->nodeInfo()->numberBranchesLeft(),x->nodeInfo()->numberPointingToThis());*/
  assert(x->objectiveValue()!=COIN_DBL_MAX&&x->nodeInfo());
  nodes_.push_back(x);
  push_heap(nodes_.begin(), nodes_.end(), comparison_);
}

// Remove the top node from the heap
void 
CbcTree::pop() {
  pop_heap(nodes_.begin(), nodes_.end(), comparison_);
  nodes_.pop_back();
}

// Test if empty *** note may be overridden
bool 
CbcTree::empty() 
{ 
  return nodes_.empty();
}
// Gets best node and takes off heap
CbcNode * 
CbcTree::bestNode(double cutoff)
{
  CbcNode * best = NULL;
  while (!best&&nodes_.size()) {
    best = nodes_.front();
    if (best)
      assert(best->objectiveValue()!=COIN_DBL_MAX&&best->nodeInfo());
    if (best&&best->objectiveValue()!=COIN_DBL_MAX&&best->nodeInfo())
      assert (best->nodeInfo()->numberBranchesLeft());
    if (!best||best->objectiveValue()>=cutoff) {
#if 0
      // take off
      pop_heap(nodes_.begin(), nodes_.end(), comparison_);
      nodes_.pop_back();
      delete best;
      best=NULL;
#else
      // let code get rid of it
      assert (best);
#endif
    }
  }
  // switched off for now
  if (best&&comparison_.test_->fullScan()&&false) {
    CbcNode * saveBest=best;
    int n=nodes_.size();
    int iBest=-1;
    for (int i=0;i<n;i++) {
      // temp
      assert (nodes_[i]);
      assert (nodes_[i]->nodeInfo());
      if (nodes_[i]&&nodes_[i]->objectiveValue()!=COIN_DBL_MAX&&nodes_[i]->nodeInfo())
        assert (nodes_[i]->nodeInfo()->numberBranchesLeft());
      if (nodes_[i]&&nodes_[i]->objectiveValue()<cutoff
          &&comparison_.alternateTest(best,nodes_[i])) {
        best=nodes_[i];
        iBest=i;
      }
    }
    if (best==saveBest) {
      // can pop
      // take off
      pop_heap(nodes_.begin(), nodes_.end(), comparison_);
      nodes_.pop_back();
    } else {
      // make impossible
      nodes_[iBest]=NULL;
    }
  } else if (best) {
    // take off
    pop_heap(nodes_.begin(), nodes_.end(), comparison_);
    nodes_.pop_back();
  }
  return best;
}

/*! \brief Prune the tree using an objective function cutoff

  This routine removes all nodes with objective worst than the
  specified cutoff value.
*/

void 
CbcTree::cleanTree(CbcModel * model, double cutoff, double & bestPossibleObjective)
{
  int j;
  int nNodes = size();
  CbcNode ** nodeArray = new CbcNode * [nNodes];
  int * depth = new int [nNodes];
  int k=0;
  int kDelete=nNodes;
  bestPossibleObjective = 1.0e100 ;
/*
    Destructively scan the heap. Nodes to be retained go into the front of
    nodeArray, nodes to be deleted into the back. Store the depth in a
    correlated array for nodes to be deleted.
*/
  for (j=0;j<nNodes;j++) {
    CbcNode * node = top();
    pop();
    double value = node ? node->objectiveValue() : COIN_DBL_MAX;
    bestPossibleObjective = CoinMin(bestPossibleObjective,value);
    if (value >= cutoff) {
      if (node) {
        nodeArray[--kDelete] = node;
        depth[kDelete] = node->depth();
      }
    } else {
      nodeArray[k++]=node;
    }
  }
/*
  Rebuild the heap using the retained nodes.
*/
  for (j=0;j<k;j++) { push(nodeArray[j]); }
/*
  Sort the list of nodes to be deleted, nondecreasing.
*/
  CoinSort_2(depth+kDelete,depth+nNodes,nodeArray+kDelete);
/*
  Work back from deepest to shallowest. In spite of the name, addCuts1 is
  just a preparatory step. When it returns, the following will be true:
    * all cuts are removed from the solver's copy of the constraint system;
    * lastws will be a basis appropriate for the specified node;
    * variable bounds will be adjusted to be appropriate for the specified
      node;
    * addedCuts_ (returned via addedCuts()) will contain a list of cuts that
      should be added to the constraint system at this node (but they have
      not actually been added).
  Then we scan the cut list for the node. Decrement the reference count
  for the cut, and if it's gone to 0, really delete it.

  I don't yet see why the checks for status != basic and addedCuts_[i] != 0
  are necessary. When reconstructing a node, these checks are used to skip
  over loose cuts, excluding them from the reconstituted basis. But here
  we're just interested in correcting the reference count. Tight/loose should
  make no difference.

  Arguably a separate routine should be used in place of addCuts1. It's doing
  more work than needed, modifying the model to match a subproblem at a node
  that will be discarded.  Then again, we seem to need the basis.
*/
  for (j=nNodes-1;j >= kDelete;j--) {
    CbcNode * node = nodeArray[j];
    CoinWarmStartBasis *lastws = model->getEmptyBasis() ;
    
    model->addCuts1(node,lastws);
    // Decrement cut counts 
    assert (node);
    //assert (node->nodeInfo());
    int numberLeft = (node->nodeInfo()) ? node->nodeInfo()->numberBranchesLeft() : 0;
    int i;
    for (i=0;i<model->currentNumberCuts();i++) {
      // take off node
      CoinWarmStartBasis::Status status = 
	lastws->getArtifStatus(i+model->numberRowsAtContinuous());
      if (status != CoinWarmStartBasis::basic&&
	  model->addedCuts()[i]) {
	if (!model->addedCuts()[i]->decrement(numberLeft))
	  delete model->addedCuts()[i];
      }
    }
    // node should not have anything pointing to it
    if (node->nodeInfo())   
      node->nodeInfo()->throwAway();
    delete node ;
    delete lastws ;
  }
  delete [] nodeArray;
  delete [] depth;
}

// Return the best node of the heap using alternate criterion
CbcNode * 
CbcTree::bestAlternate() {
  int n=nodes_.size();
  CbcNode * best=NULL;
  if (n) {
    best = nodes_[0];
    for (int i=1;i<n;i++) {
      if (comparison_.alternateTest(best,nodes_[i])) {
        best=nodes_[i];
      }
    }
  }
  return best;
}

