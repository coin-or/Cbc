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
CbcTree::top() {
  return nodes_.front();
}

// Add a node to the heap
void 
CbcTree::push(CbcNode * x) {
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

/*! \brief Prune the tree using an objective function cutoff

  This routine removes all nodes with objective worst than the
  specified cutoff value.
*/

void 
CbcTree::cleanTree(CbcModel * model, double cutoff)
{
  int j;
  int nNodes = size();
  CbcNode ** nodeArray = new CbcNode * [nNodes];
  int * depth = new int [nNodes];
  int k=0;
  int kDelete=nNodes;
/*
    Destructively scan the heap. Nodes to be retained go into the front of
    nodeArray, nodes to be deleted into the back. Store the depth in a
    correlated array for nodes to be deleted.
*/
  for (j=0;j<nNodes;j++) {
    CbcNode * node = top();
    pop();
    if (node->objectiveValue() >= cutoff) {
      nodeArray[--kDelete] = node;
      depth[kDelete] = node->depth();
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
    int numberLeft = node->nodeInfo()->numberBranchesLeft();
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
    node->nodeInfo()->throwAway();
    delete node ;
    delete lastws ;
  }
  delete [] nodeArray;
  delete [] depth;
};


