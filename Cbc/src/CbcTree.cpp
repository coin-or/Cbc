/* $Id$ */
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcModel.hpp"
#include "CbcNode.hpp"
#include "CbcTree.hpp"
#include "CbcThread.hpp"
#include "CbcCountRowCut.hpp"
#include "CbcCompareActual.hpp"
#include "CbcBranchActual.hpp"

#if CBC_DEBUG_HEAP > 0

namespace {

/*
  The next few methods test that the heap property (parent equal or better
  than either child) is maintained in the heap. Originally created to sort out
  why `cbc -unitTest' triggered an `Invalid heap' error in a MSVS debug build.
*/
/*
  Predicate test. The parent should be better or equal to the child. Since the
  predicate comparison_(x,y) returns true if y (child) is strictly better,
  we want failure on the initial test. Clearly, success for comparison(x,y)
  and comparison(y,x) is also a failure.

  Returns true if the predicate passes, false if it fails.
*/
bool check_pred(CbcCompareBase &pred, CbcNode *parent, CbcNode *child)
{
  if (parent == 0 || child == 0)
    return (false);
  if (!pred(parent, child))
    return (true);
  else if (pred(child, parent))
    std::cout
      << " Heap predicate failure! (x<y) && (y<x)!" << std::endl;
  return (false);
}

} // end file-local namespace

/*
  Check for the heap property: the parent is better than or equal to
  either child.

  The heap is a binary tree, stored in the vector layer by layer. By advancing
  parent at half the rate of child (aka curNode), we check both children
  of a given parent.  (Draw yourself a picture; it'll help.) An empty heap
  is trivially valid. A heap with no predicate is trivially invalid.

  TODO: The heap -> vector mapping assumed here is valid for the MSVS heap
	implementation.  No guarantee it's valid elsewhere.
*/

void CbcTree::validateHeap()
{
  if (comparison_.test_ == 0) {
    std::cout
      << " Invalid heap (no predicate)!" << std::endl;
    return;
  }
  std::vector< CbcNode * >::const_iterator curNode, lastNode;
  curNode = nodes_.begin();
  lastNode = nodes_.end();
  int curNdx = 0;
  int parNdx = 0;
  if (curNode == lastNode)
    return;
  if (*curNode == 0) {
    std::cout
      << " Invalid heap[" << curNdx << "] (null entry)!" << std::endl;
  }
  std::vector< CbcNode * >::const_iterator parent;
  std::vector< CbcNode * >::const_iterator &child = curNode;
  for (parent = curNode; ++curNode != lastNode; ++parent, ++parNdx) {
    curNdx++;
    if (*parent == 0) {
      std::cout
        << " Invalid heap[" << parNdx << "] (parent null entry)!" << std::endl;
      curNode++;
      curNdx++;
      continue;
    }
    if (*curNode == 0) {
      std::cout
        << " Invalid heap[" << curNdx << "] (left child null entry)!"
        << std::endl;
    } else {
      if (!check_pred(*comparison_.test_, *parent, *child)) {
        std::cout
          << " Invalid heap (left child better)!" << std::endl;
        CbcNode *node = *parent;
        std::cout
          << "   Parent [" << parNdx << "] (" << std::hex << node << std::dec
          << ") unsat " << node->numberUnsatisfied() << ", depth "
          << node->depth() << ", obj " << node->objectiveValue() << "."
          << std::endl;
        node = *child;
        std::cout
          << "   Child [" << curNdx << "] (" << std::hex << node << std::dec
          << ") unsat " << node->numberUnsatisfied() << ", depth "
          << node->depth() << ", obj " << node->objectiveValue() << "."
          << std::endl;
      }
    }
    curNode++;
    curNdx++;
    if (curNode == lastNode)
      break;
    if (*curNode == 0) {
      std::cout
        << " Invalid heap[" << curNdx << "] (right child null entry)!"
        << std::endl;
    } else {
      if (!check_pred(*comparison_.test_, *parent, *child)) {
        std::cout
          << " Invalid heap (right child better)!" << std::endl;
        CbcNode *node = *parent;
        std::cout
          << "   Parent [" << parNdx << "] (" << std::hex << node << std::dec
          << ") unsat " << node->numberUnsatisfied() << ", depth "
          << node->depth() << ", obj " << node->objectiveValue() << "."
          << std::endl;
        node = *child;
        std::cout
          << "   Child [" << curNdx << "] (" << std::hex << node << std::dec
          << ") unsat " << node->numberUnsatisfied() << ", depth "
          << node->depth() << ", obj " << node->objectiveValue() << "."
          << std::endl;
      }
    }
  }
  return;
}

#endif // CBC_DEBUG_HEAP

CbcTree::CbcTree()
{
  maximumNodeNumber_ = 0;
  numberBranching_ = 0;
  maximumBranching_ = 0;
  branched_ = NULL;
  newBound_ = NULL;
}
CbcTree::~CbcTree()
{
  delete[] branched_;
  delete[] newBound_;
}
// Copy constructor
CbcTree::CbcTree(const CbcTree &rhs)
{
  nodes_ = rhs.nodes_;
  maximumNodeNumber_ = rhs.maximumNodeNumber_;
  numberBranching_ = rhs.numberBranching_;
  maximumBranching_ = rhs.maximumBranching_;
  if (maximumBranching_ > 0) {
    branched_ = CoinCopyOfArray(rhs.branched_, maximumBranching_);
    newBound_ = CoinCopyOfArray(rhs.newBound_, maximumBranching_);
  } else {
    branched_ = NULL;
    newBound_ = NULL;
  }
}
// Assignment operator
CbcTree &
CbcTree::operator=(const CbcTree &rhs)
{
  if (this != &rhs) {
    nodes_ = rhs.nodes_;
    maximumNodeNumber_ = rhs.maximumNodeNumber_;
    delete[] branched_;
    delete[] newBound_;
    numberBranching_ = rhs.numberBranching_;
    maximumBranching_ = rhs.maximumBranching_;
    if (maximumBranching_ > 0) {
      branched_ = CoinCopyOfArray(rhs.branched_, maximumBranching_);
      newBound_ = CoinCopyOfArray(rhs.newBound_, maximumBranching_);
    } else {
      branched_ = NULL;
      newBound_ = NULL;
    }
  }
  return *this;
}

/*
  Rebuild the heap.
*/
void CbcTree::rebuild()
{
  std::make_heap(nodes_.begin(), nodes_.end(), comparison_);
#if CBC_DEBUG_HEAP > 1
  std::cout << "  HEAP: rebuild complete." << std::endl;
#endif
#if CBC_DEBUG_HEAP > 0
  validateHeap();
#endif
}

// Adds branching information to complete state
void CbcTree::addBranchingInformation(const CbcModel *model, const CbcNodeInfo *nodeInfo,
  const double *currentLower,
  const double *currentUpper)
{
  const OsiBranchingObject *objA = nodeInfo->owner()->branchingObject();
  const CbcIntegerBranchingObject *objBranch = dynamic_cast< const CbcIntegerBranchingObject * >(objA);
  if (objBranch) {
    const CbcObject *objB = objBranch->object();
    const CbcSimpleInteger *obj = dynamic_cast< const CbcSimpleInteger * >(objB);
    assert(obj);
    int iColumn = obj->columnNumber();
    const double *down = objBranch->downBounds();
    const double *up = objBranch->upBounds();
    assert(currentLower[iColumn] == down[0]);
    assert(currentUpper[iColumn] == up[1]);
    if (dynamic_cast< const CbcPartialNodeInfo * >(nodeInfo)) {
      const CbcPartialNodeInfo *info = dynamic_cast< const CbcPartialNodeInfo * >(nodeInfo);
      const double *newBounds = info->newBounds();
      const int *variables = info->variables();
      int numberChanged = info->numberChangedBounds();
      for (int i = 0; i < numberChanged; i++) {
        int jColumn = variables[i];
        int kColumn = jColumn & (~0x80000000);
        if (iColumn == kColumn) {
          jColumn |= 0x40000000;
#ifndef NDEBUG
          double value = newBounds[i];
          if ((jColumn & 0x80000000) == 0) {
            assert(value == up[0]);
          } else {
            assert(value == down[1]);
          }
#endif
        }
        if (numberBranching_ == maximumBranching_)
          increaseSpace();
        newBound_[numberBranching_] = static_cast< int >(newBounds[i]);
        branched_[numberBranching_++] = jColumn;
      }
    } else {
      const CbcFullNodeInfo *info = dynamic_cast< const CbcFullNodeInfo * >(nodeInfo);
      int numberIntegers = model->numberIntegers();
      const int *which = model->integerVariable();
      const double *newLower = info->lower();
      const double *newUpper = info->upper();
      if (numberBranching_ == maximumBranching_)
        increaseSpace();
      assert(newLower[iColumn] == up[0] || newUpper[iColumn] == down[1]);
      int jColumn = iColumn | 0x40000000;
      if (newLower[iColumn] == up[0]) {
        newBound_[numberBranching_] = static_cast< int >(up[0]);
      } else {
        newBound_[numberBranching_] = static_cast< int >(down[1]);
        jColumn |= 0x80000000;
      }
      branched_[numberBranching_++] = jColumn;
      for (int i = 0; i < numberIntegers; i++) {
        int jColumn = which[i];
        assert(currentLower[jColumn] == newLower[jColumn] || currentUpper[jColumn] == newUpper[jColumn]);
        if (jColumn != iColumn) {
          bool changed = false;
          double value;
          if (newLower[jColumn] > currentLower[jColumn]) {
            value = newLower[jColumn];
            changed = true;
          } else if (newUpper[jColumn] < currentUpper[jColumn]) {
            value = newUpper[jColumn];
            jColumn |= 0x80000000;
            changed = true;
          }
          if (changed) {
            if (numberBranching_ == maximumBranching_)
              increaseSpace();
            newBound_[numberBranching_] = static_cast< int >(value);
            branched_[numberBranching_++] = jColumn;
          }
        }
      }
    }
  } else {
    // switch off
    delete[] branched_;
    delete[] newBound_;
    maximumBranching_ = -1;
    branched_ = NULL;
    newBound_ = NULL;
  }
}
// Increase space for data
void CbcTree::increaseSpace()
{
  assert(numberBranching_ == maximumBranching_);
  maximumBranching_ = (3 * maximumBranching_ + 10) >> 1;
  unsigned int *temp1 = CoinCopyOfArrayPartial(branched_, maximumBranching_, numberBranching_);
  delete[] branched_;
  branched_ = temp1;
  int *temp2 = CoinCopyOfArrayPartial(newBound_, maximumBranching_, numberBranching_);
  delete[] newBound_;
  newBound_ = temp2;
}
// Clone
CbcTree *
CbcTree::clone() const
{
  return new CbcTree(*this);
}

#ifndef CBC_DUBIOUS_HEAP
/*
  Set comparison predicate and re-sort the heap.

  Note that common usage is to tweak the incumbent predicate and then
  call this method to rebuild the heap. Hence we cannot check for heap
  validity at entry. rebuild() will check on the way out, if
  CBC_DEBUG_HEAP is set.

  TODO: remove the call to cleanDive and put it somewhere appropriate.
*/
void CbcTree::setComparison(CbcCompareBase &compare)
{
#if CBC_DEBUG_HEAP > 1
  std::cout << "  HEAP: resetting comparison predicate." << std::endl;
#endif
  comparison_.test_ = &compare;

  /*
  From a software engineering point of view, setComparison has no business
  knowing anything about the comparison function. Need to look for a better
  solution. Perhaps a callback comparable to newSolution, executing when
  the comparison method is set (i.e., in setComparison).
  -- lh, 100921 --
*/
  CbcCompareDefault *compareD = dynamic_cast< CbcCompareDefault * >(&compare);
  if (compareD) {
    // clean up diving
    compareD->cleanDive();
  }
  rebuild();
}

// Return the top node of the heap
CbcNode *
CbcTree::top() const
{
  return nodes_.front();
}

// Add a node to the heap
void CbcTree::push(CbcNode *x)
{
  x->setNodeNumber(maximumNodeNumber_);
  lastObjective_ = x->objectiveValue();
  lastDepth_ = x->depth();
  lastUnsatisfied_ = x->numberUnsatisfied();
  maximumNodeNumber_++;
#if CBC_DEBUG_HEAP > 2
  CbcNodeInfo *info = x->nodeInfo();
  assert(info);
  std::cout
    << "  HEAP: Pushing node " << x->nodeNumber()
    << "(" << std::hex << x << std::dec << ") obj " << x->objectiveValue()
    << ", ref " << info->decrement(0)
    << ", todo " << info->numberBranchesLeft()
    << ", refd by " << info->numberPointingToThis() << "." << std::endl;
  assert(x->objectiveValue() != COIN_DBL_MAX);
#endif
#if CBC_DEBUG_HEAP > 0
  validateHeap();
#endif
  x->setOnTree(true);
  nodes_.push_back(x);
  std::push_heap(nodes_.begin(), nodes_.end(), comparison_);
#if CBC_DEBUG_HEAP > 0
  validateHeap();
#endif
}

// Remove the top node from the heap
void CbcTree::pop()
{

#if CBC_DEBUG_HEAP > 2
  CbcNode *node = nodes_.front();
  CbcNodeInfo *info = node->nodeInfo();
  assert(info);
  std::cout
    << "  HEAP: Popping node " << node->nodeNumber()
    << "(" << std::hex << node << std::dec
    << ") obj " << node->objectiveValue()
    << ", ref " << info->decrement(0)
    << ", todo " << info->numberBranchesLeft()
    << ", refd by " << info->numberPointingToThis() << "." << std::endl;
#endif
#if CBC_DEBUG_HEAP > 0
  validateHeap();
#endif
  nodes_.front()->setOnTree(false);
  std::pop_heap(nodes_.begin(), nodes_.end(), comparison_);
  nodes_.pop_back();

#if CBC_DEBUG_HEAP > 0
  validateHeap();
#endif
}

// Test if empty *** note may be overridden
bool CbcTree::empty()
{
  return nodes_.empty();
}
/*
  Return the best node from the heap.

  Note that best is offered a chance (checkIsCutoff) to reevaluate
  itself and make arbitrary changes. A caller should be prepared
  to check that the returned node is acceptable.

  There's quite a bit of suspect code here, much of it disabled in
  some way. The net effect at present is to return the top node on
  the heap after offering the node an opportunity to reevaluate itself.
  Documentation for checkIsCutoff() puts no restrictions on allowable
  changes. -- lh, 100921 --
*/
CbcNode *
CbcTree::bestNode(double cutoff)
{
#if CBC_DEBUG_HEAP > 0
  validateHeap();
#endif
  /*
  This code is problematic. As of 100921, there's really no loop.
  If front() == null, an assert will trigger. checkIsCutoff seems to be
  work in progress; comments assert that it can make pretty much arbitrary
  changes to best. If best can change its objective, there's a good
  possibility the heap is invalid.
*/
  CbcNode *best = NULL;
  while (!best && nodes_.size()) {
    best = nodes_.front();
    if (best)
      assert(best->objectiveValue() != COIN_DBL_MAX && best->nodeInfo());
    if (best && best->objectiveValue() != COIN_DBL_MAX && best->nodeInfo())
      assert(best->nodeInfo()->numberBranchesLeft());
    if (best && best->objectiveValue() >= cutoff) {
      // double check in case node can change its mind!
      best->checkIsCutoff(cutoff);
    }
    if (!best || best->objectiveValue() >= cutoff) {
#ifdef JJF_ZERO
      // take off
      std::pop_heap(nodes_.begin(), nodes_.end(), comparison_);
      nodes_.pop_back();
      delete best;
      best = NULL;
#else
      // let code get rid of it
      assert(best);
#endif
    }
  }
  /*
  See if the comparison object wants us to do a full scan with the
  alternative criteria. The net effect is to confirm best by the
  alternative criteria, or identify a competitor and erase it.

  This code is problematic. Nulling an arbitrary node will in general
  break the heap property. Disabled some time ago, as noted in several
  places.
*/
  if (false && best && comparison_.test_->fullScan()) {
    CbcNode *saveBest = best;
    size_t n = nodes_.size();
    size_t iBest = -1;
    for (size_t i = 0; i < n; i++) {
      // temp
      assert(nodes_[i]);
      assert(nodes_[i]->nodeInfo());
      if (nodes_[i] && nodes_[i]->objectiveValue() != COIN_DBL_MAX && nodes_[i]->nodeInfo())
        assert(nodes_[i]->nodeInfo()->numberBranchesLeft());
      if (nodes_[i] && nodes_[i]->objectiveValue() < cutoff
        && comparison_.alternateTest(best, nodes_[i])) {
        best = nodes_[i];
        iBest = i;
      }
    }
    if (best == saveBest) {
      // can pop
      // take off
      std::pop_heap(nodes_.begin(), nodes_.end(), comparison_);
      nodes_.pop_back();
    } else {
      // make impossible
      nodes_[iBest] = NULL;
    }
  } else if (best) {
#if CBC_DEBUG_HEAP > 2
    CbcNode *node = nodes_.front();
    CbcNodeInfo *info = node->nodeInfo();
    assert(info);
    std::cout
      << "  bestNode: Popping node " << node->nodeNumber()
      << "(" << std::hex << node << std::dec
      << ") obj " << node->objectiveValue()
      << ", ref " << info->decrement(0)
      << ", todo " << info->numberBranchesLeft()
      << ", refd by " << info->numberPointingToThis() << "." << std::endl;
#endif
    // take off
    std::pop_heap(nodes_.begin(), nodes_.end(), comparison_);
    nodes_.pop_back();
  }
#if CBC_DEBUG_HEAP > 0
  validateHeap();
#endif
  if (best)
    best->setOnTree(false);
  return best;
}
/*! \brief Prune the tree using an objective function cutoff

  This routine removes all nodes with objective worse than the
  specified cutoff value.
*/

void CbcTree::cleanTree(CbcModel *model, double cutoff, double &bestPossibleObjective)
{
#if CBC_DEBUG_HEAP > 1
  std::cout << " cleanTree: beginning clean." << std::endl;
#endif
#if CBC_DEBUG_HEAP > 0
  validateHeap();
#endif
  int j;
  int nNodes = size();
  CbcNode **nodeArray = new CbcNode *[nNodes];
  int *depth = new int[nNodes];
  int k = 0;
  int kDelete = nNodes;
  bestPossibleObjective = 1.0e100;
  /*
        Destructively scan the heap. Nodes to be retained go into the front of
        nodeArray, nodes to be deleted into the back. Store the depth in a
        correlated array for nodes to be deleted.
    */
  for (j = 0; j < nNodes; j++) {
    CbcNode *node = top();
    pop();
    double value = node ? node->objectiveValue() : COIN_DBL_MAX;
    if (node && value >= cutoff) {
      // double check in case node can change its mind!
      value = node->checkIsCutoff(cutoff);
    }
    if (value >= cutoff || !node->active()) {
      if (node) {
        if (cutoff < -1.0e30)
          node->nodeInfo()->deactivate(7);
        nodeArray[--kDelete] = node;
        depth[kDelete] = node->depth();
      }
    } else {
      bestPossibleObjective = CoinMin(bestPossibleObjective, value);
      nodeArray[k++] = node;
    }
  }
  /*
      Rebuild the heap using the retained nodes.
    */
  for (j = 0; j < k; j++) {
    push(nodeArray[j]);
  }
#if CBC_DEBUG_HEAP > 1
  std::cout << " cleanTree: finished rebuild." << std::endl;
#endif
#if CBC_DEBUG_HEAP > 0
  validateHeap();
#endif
  /*
      Sort the list of nodes to be deleted, nondecreasing.
    */
  CoinSort_2(depth + kDelete, depth + nNodes, nodeArray + kDelete);
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
      we're just interested in correcting the reference count. Tight/loose
      should make no difference.

      Arguably a separate routine should be used in place of addCuts1. It's
      doing more work than needed, modifying the model to match a subproblem
      at a node that will be discarded.  Then again, we seem to need the basis.
    */
  for (j = nNodes - 1; j >= kDelete; j--) {
    CbcNode *node = nodeArray[j];
    CoinWarmStartBasis *lastws = (cutoff != -COIN_DBL_MAX) ? model->getEmptyBasis() : NULL;

    model->addCuts1(node, lastws);
    // Decrement cut counts
    assert(node);
    //assert (node->nodeInfo());
    int numberLeft = (node->nodeInfo()) ? node->nodeInfo()->numberBranchesLeft() : 0;
    if (cutoff != -COIN_DBL_MAX) {
      // normal
      for (int i = 0; i < model->currentNumberCuts(); i++) {
        // take off node
        CoinWarmStartBasis::Status status = lastws->getArtifStatus(i + model->numberRowsAtContinuous());
        if (status != CoinWarmStartBasis::basic && model->addedCuts()[i]) {
          if (!model->addedCuts()[i]->decrement(numberLeft))
            delete model->addedCuts()[i];
        }
      }
    } else {
      // quick
      for (int i = 0; i < model->currentNumberCuts(); i++) {
        // take off node
        if (model->addedCuts()[i]) {
          if (model->parallelMode() != 1 || true) {
            if (!model->addedCuts()[i]->decrement(numberLeft))
              delete model->addedCuts()[i];
          }
        }
      }
    }
#ifdef CBC_THREAD
    if (model->parallelMode() > 0 && model->master()) {
      // delete reference to node
      int numberThreads = model->master()->numberThreads();
      for (int i = 0; i < numberThreads; i++) {
        CbcThread *child = model->master()->child(i);
        if (child->createdNode() == node)
          child->setCreatedNode(NULL);
      }
    }
#endif
    // node should not have anything pointing to it
    if (node->nodeInfo())
      node->nodeInfo()->throwAway();
    model->deleteNode(node);
    delete lastws;
  }
  delete[] nodeArray;
  delete[] depth;
#ifdef CBC_THREAD
  if (model->parallelMode() > 0 && model->master()) {
    // need to adjust for ones not on tree
    CbcBaseModel *master = model->master();
    int numberThreads = master->numberThreads();
    for (int i = 0; i < numberThreads; i++) {
      CbcThread *child = master->child(i);
      if (child->node()) {
        double value = child->node()->objectiveValue();
        // adjust
        bestPossibleObjective = CoinMin(bestPossibleObjective, value);
      }
    }
  }
#endif
}

// Return the best node of the heap using alternate criterion
CbcNode *
CbcTree::bestAlternate()
{
  size_t n = nodes_.size();
  CbcNode *best = NULL;
  if (n) {
    best = nodes_[0];
    for (size_t i = 1; i < n; i++) {
      if (comparison_.alternateTest(best, nodes_[i])) {
        best = nodes_[i];
      }
    }
  }
  return best;
}

#ifdef JJF_ZERO // not used, reference removed in CbcModel.cpp
CbcTreeArray::CbcTreeArray()
  : CbcTree()
  , lastNode_(NULL)
  , lastNodePopped_(NULL)
  , switches_(0)
{
}
CbcTreeArray::~CbcTreeArray()
{
}
// Copy constructor
CbcTreeArray::CbcTreeArray(const CbcTreeArray &rhs)
  : CbcTree(rhs)
  , lastNode_(rhs.lastNode_)
  , lastNodePopped_(rhs.lastNodePopped_)
  , switches_(rhs.switches_)
{
}
// Assignment operator
CbcTreeArray &
CbcTreeArray::operator=(const CbcTreeArray &rhs)
{
  if (this != &rhs) {
    CbcTree::operator=(rhs);
    lastNode_ = rhs.lastNode_;
    lastNodePopped_ = rhs.lastNodePopped_;
    switches_ = rhs.switches_;
  }
  return *this;
}
// Clone
CbcTree *
CbcTreeArray::clone() const
{
  return new CbcTreeArray(*this);
}
// Set comparison function and resort heap
void CbcTreeArray::setComparison(CbcCompareBase &compare)
{
  comparison_.test_ = &compare;
  rebuild();
}

// Add a node to the heap
void CbcTreeArray::push(CbcNode *x)
{
  /*printf("push obj %g, refcount %d, left %d, pointing to %d\n",
       x->objectiveValue(),x->nodeInfo()->decrement(0),
       x->nodeInfo()->numberBranchesLeft(),x->nodeInfo()->numberPointingToThis());*/
  assert(x->objectiveValue() != COIN_DBL_MAX && x->nodeInfo());
  x->setOnTree(true);
  if (lastNode_) {
    if (lastNode_->nodeInfo()->parent() == x->nodeInfo()) {
      // x is parent of lastNode_ so put x on heap
      //#define CBCTREE_PRINT
#ifdef CBCTREE_PRINT
      printf("pushX x %x (%x at depth %d n %d) is parent of lastNode_ %x (%x at depth %d n %d)\n",
        x, x->nodeInfo(), x->depth(), x->nodeNumber(),
        lastNode_, lastNode_->nodeInfo(), lastNode_->depth(), lastNode_->nodeNumber());
#endif
      nodes_.push_back(x);
    } else {
      x->setNodeNumber(maximumNodeNumber_);
      maximumNodeNumber_++;
#ifdef CBCTREE_PRINT
      printf("pushLast x %x (%x at depth %d n %d) is parent of lastNode_ %x (%x at depth %d n %d)\n",
        x, x->nodeInfo(), x->depth(), x->nodeNumber(),
        lastNode_, lastNode_->nodeInfo(), lastNode_->depth(), lastNode_->nodeNumber());
#endif
      nodes_.push_back(lastNode_);
      lastNode_ = x;
    }
    std::push_heap(nodes_.begin(), nodes_.end(), comparison_);
  } else {
    x->setNodeNumber(maximumNodeNumber_);
    maximumNodeNumber_++;
    if (x != lastNodePopped_) {
      lastNode_ = x;
#ifdef CBCTREE_PRINT
      printf("pushNULL x %x (%x at depth %d n %d)\n",
        x, x->nodeInfo(), x->depth(), x->nodeNumber());
#endif
    } else {
      // means other way was infeasible
#ifdef CBCTREE_PRINT
      printf("push_other_infeasible x %x (%x at depth %d n %d)\n",
        x, x->nodeInfo(), x->depth(), x->nodeNumber());
#endif
      nodes_.push_back(x);
      std::push_heap(nodes_.begin(), nodes_.end(), comparison_);
    }
  }
}

// Test if empty *** note may be overridden
bool CbcTreeArray::empty()
{
  return nodes_.empty() && (lastNode_ == NULL);
}
// Gets best node and takes off heap
CbcNode *
CbcTreeArray::bestNode(double cutoff)
{
  CbcNode *best = NULL;
  // See if we want last node or best on heap
  if (lastNode_) {
#ifdef CBCTREE_PRINT
    printf("Best lastNode_ %x (%x at depth %d) - nodeNumber %d obj %g\n",
      lastNode_, lastNode_->nodeInfo(), lastNode_->depth(),
      lastNode_->nodeNumber(), lastNode_->objectiveValue());
#endif
    assert(lastNode_->onTree());
    int nodeNumber = lastNode_->nodeNumber();
    bool useLastNode = false;
    if (nodeNumber + 1 == maximumNodeNumber_) {
      // diving - look further
      CbcCompareDefault *compareDefault
        = dynamic_cast< CbcCompareDefault * >(comparison_.test_);
      assert(compareDefault);
      double bestPossible = compareDefault->getBestPossible();
      double cutoff = compareDefault->getCutoff();
      double objValue = lastNode_->objectiveValue();
      if (cutoff < 1.0e20) {
        if (objValue - bestPossible < 0.999 * (cutoff - bestPossible))
          useLastNode = true;
      } else {
        useLastNode = true;
      }
    }
    if (useLastNode) {
      lastNode_->setOnTree(false);
      best = lastNode_;
      lastNode_ = NULL;
      assert(best->objectiveValue() != COIN_DBL_MAX && best->nodeInfo());
      if (best->objectiveValue() != COIN_DBL_MAX && best->nodeInfo())
        assert(best->nodeInfo()->numberBranchesLeft());
      if (best->objectiveValue() >= cutoff) {
        // double check in case node can change its mind!
        best->checkIsCutoff(cutoff);
      }
      lastNodePopped_ = best;
      return best;
    } else {
      // put on tree
      nodes_.push_back(lastNode_);
      lastNode_->setNodeNumber(maximumNodeNumber_);
      maximumNodeNumber_++;
      lastNode_ = NULL;
      std::push_heap(nodes_.begin(), nodes_.end(), comparison_);
    }
  }
  while (!best && nodes_.size()) {
    best = nodes_.front();
    if (best)
      assert(best->objectiveValue() != COIN_DBL_MAX && best->nodeInfo());
    if (best && best->objectiveValue() != COIN_DBL_MAX && best->nodeInfo())
      assert(best->nodeInfo()->numberBranchesLeft());
    if (best && best->objectiveValue() >= cutoff) {
      // double check in case node can change its mind!
      best->checkIsCutoff(cutoff);
    }
    if (!best || best->objectiveValue() >= cutoff) {
      // let code get rid of it
      assert(best);
    }
  }
  lastNodePopped_ = best;
#ifdef CBCTREE_PRINT
  if (best)
    printf("Heap returning node %x (%x at depth %d) - nodeNumber %d - obj %g\n",
      best, best->nodeInfo(), best->depth(),
      best->nodeNumber(), best->objectiveValue());
  else
    printf("Heap returning Null\n");
#endif
  if (best) {
    // take off
    std::pop_heap(nodes_.begin(), nodes_.end(), comparison_);
    nodes_.pop_back();
  }
#ifdef DEBUG_CBC_HEAP
  if (best) {
    int n = nodes_.size();
    bool good = true;
    for (int i = 0; i < n; i++) {
      // temp
      assert(nodes_[i]);
      if (!comparison_.compareNodes(nodes_[i], best)) {
        good = false;
        CbcNode *x = nodes_[i];
        printf("i=%d x is better nun %d depth %d obj %g, best nun %d depth %d obj %g\n", i,
          x->numberUnsatisfied(), x->depth(), x->objectiveValue(),
          best->numberUnsatisfied(), best->depth(), best->objectiveValue());
      }
    }
    if (!good) {
      // compare best to all
      int i;
      for (i = 0; i < n; i++) {
        CbcNode *x = nodes_[i];
        printf("i=%d x is nun %d depth %d obj %g", i,
          x->numberUnsatisfied(), x->depth(), x->objectiveValue());
        if (!comparison_.compareNodes(x, best)) {
          printf(" - best is worse!\n");
        } else {
          printf("\n");
        }
      }
      // Now compare amongst rest
      for (i = 0; i < n; i++) {
        CbcNode *x = nodes_[i];
        printf("For i=%d ", i);
        for (int j = i + 1; j < n; j++) {
          CbcNode *y = nodes_[j];
          if (!comparison_.compareNodes(x, y)) {
            printf(" b %d", j);
          } else {
            printf(" w %d", j);
          }
        }
        printf("\n");
      }
      assert(good);
    }
  }
#endif
  if (best)
    best->setOnTree(false);
  return best;
}

double
CbcTreeArray::getBestPossibleObjective()
{
  double bestPossibleObjective = 1e100;
  for (int i = 0; i < static_cast< int >(nodes_.size()); i++) {
    if (nodes_[i] && nodes_[i]->objectiveValue() < bestPossibleObjective) {
      bestPossibleObjective = nodes_[i]->objectiveValue();
    }
  }
  if (lastNode_) {
    bestPossibleObjective = CoinMin(bestPossibleObjective, lastNode_->objectiveValue());
  }
#ifdef CBC_THREAD
  if (model->parallelMode() > 0 && model->master()) {
    // need to adjust for ones not on tree
    CbcBaseModel *master = model->master();
    int numberThreads = master->numberThreads();
    for (int i = 0; i < numberThreads; i++) {
      CbcThread *child = master->child(i);
      if (child->node()) {
        double value = child->node()->objectiveValue();
        // adjust
        bestPossibleObjective = CoinMin(bestPossibleObjective, value);
      }
    }
  }
#endif
  CbcCompareDefault *compareDefault
    = dynamic_cast< CbcCompareDefault * >(comparison_.test_);
  assert(compareDefault);
  compareDefault->setBestPossible(bestPossibleObjective);
  return bestPossibleObjective;
}
/*! \brief Prune the tree using an objective function cutoff

  This routine removes all nodes with objective worst than the
  specified cutoff value.
*/

void CbcTreeArray::cleanTree(CbcModel *model, double cutoff, double &bestPossibleObjective)
{
  int j;
  int nNodes = size();
  int lastNode = nNodes + 1;
  CbcNode **nodeArray = new CbcNode *[lastNode];
  int *depth = new int[lastNode];
  int k = 0;
  int kDelete = lastNode;
  bestPossibleObjective = 1.0e100;
  /*
        Destructively scan the heap. Nodes to be retained go into the front of
        nodeArray, nodes to be deleted into the back. Store the depth in a
        correlated array for nodes to be deleted.
    */
  for (j = 0; j < nNodes; j++) {
    CbcNode *node = nodes_.front();
    nodes_.front()->setOnTree(false);
    std::pop_heap(nodes_.begin(), nodes_.end(), comparison_);
    nodes_.pop_back();
    double value = node ? node->objectiveValue() : COIN_DBL_MAX;
    if (node && value >= cutoff) {
      // double check in case node can change its mind!
      value = node->checkIsCutoff(cutoff);
    }
    if (value >= cutoff || !node->active()) {
      if (node) {
        nodeArray[--kDelete] = node;
        depth[kDelete] = node->depth();
      }
    } else {
      bestPossibleObjective = CoinMin(bestPossibleObjective, value);
      nodeArray[k++] = node;
    }
  }
#ifdef CBC_THREAD
  if (model->parallelMode() > 0 && model->master()) {
    // need to adjust for ones not on tree
    CbcBaseModel *master = model->master();
    int numberThreads = master->numberThreads();
    for (int i = 0; i < numberThreads; i++) {
      CbcThread *child = master->child(i);
      if (child->node()) {
        double value = child->node()->objectiveValue();
        // adjust
        bestPossibleObjective = CoinMin(bestPossibleObjective, value);
      }
    }
  }
#endif
  if (lastNode_) {
    double value = lastNode_->objectiveValue();
    bestPossibleObjective = CoinMin(bestPossibleObjective, value);
    if (value >= cutoff || !lastNode_->active()) {
      nodeArray[--kDelete] = lastNode_;
      depth[kDelete] = lastNode_->depth();
      lastNode_ = NULL;
    }
  }
  CbcCompareDefault *compareDefault
    = dynamic_cast< CbcCompareDefault * >(comparison_.test_);
  assert(compareDefault);
  compareDefault->setBestPossible(bestPossibleObjective);
  compareDefault->setCutoff(cutoff);
  /*
      Rebuild the heap using the retained nodes.
    */
  for (j = 0; j < k; j++) {
    CbcNode *node = nodeArray[j];
    node->setOnTree(true);
    nodes_.push_back(node);
    std::push_heap(nodes_.begin(), nodes_.end(), comparison_);
  }
  /*
      Sort the list of nodes to be deleted, nondecreasing.
    */
  CoinSort_2(depth + kDelete, depth + lastNode, nodeArray + kDelete);
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
      we're just interested in correcting the reference count. Tight/loose
      should make no difference.

      Arguably a separate routine should be used in place of addCuts1. It's
      doing more work than needed, modifying the model to match a subproblem
      at a node that will be discarded.  Then again, we seem to need the basis.
    */
  for (j = lastNode - 1; j >= kDelete; j--) {
    CbcNode *node = nodeArray[j];
    CoinWarmStartBasis *lastws = model->getEmptyBasis();

    model->addCuts1(node, lastws);
    // Decrement cut counts
    assert(node);
    //assert (node->nodeInfo());
    int numberLeft = (node->nodeInfo()) ? node->nodeInfo()->numberBranchesLeft() : 0;
    int i;
    for (i = 0; i < model->currentNumberCuts(); i++) {
      // take off node
      CoinWarmStartBasis::Status status = lastws->getArtifStatus(i + model->numberRowsAtContinuous());
      if (status != CoinWarmStartBasis::basic && model->addedCuts()[i]) {
        if (!model->addedCuts()[i]->decrement(numberLeft))
          delete model->addedCuts()[i];
      }
    }
    // node should not have anything pointing to it
    if (node->nodeInfo())
      node->nodeInfo()->throwAway();
    delete node;
    delete lastws;
  }
  delete[] nodeArray;
  delete[] depth;
}
#endif

#else // defined(CBC_DUBIOUS_HEAP)

/*
  Unclear whether this code is useful any longer. Likely stale. See
  note in CbcCompareDefault.hpp re. CBC_DUBIOUS_HEAP.
  -- lh, 100921 --
*/

// Set comparison function and resort heap
void CbcTree::setComparison(CbcCompareBase &compare)
{
  comparison_.test_ = &compare;
  std::vector< CbcNode * > newNodes = nodes_;
  nodes_.resize(0);
  while (newNodes.size() > 0) {
    push(newNodes.back());
    newNodes.pop_back();
  }
}

// Return the top node of the heap
CbcNode *
CbcTree::top() const
{
  return nodes_.front();
}

// Add a node to the heap
void CbcTree::push(CbcNode *x)
{
  x->setNodeNumber(maximumNodeNumber_);
  maximumNodeNumber_++;
  /*printf("push obj %g, refcount %d, left %d, pointing to %d\n",
       x->objectiveValue(),x->nodeInfo()->decrement(0),
       x->nodeInfo()->numberBranchesLeft(),x->nodeInfo()->numberPointingToThis());*/
  assert(x->objectiveValue() != COIN_DBL_MAX && x->nodeInfo());
#ifdef JJF_ZERO
  nodes_.push_back(x);
  push_heap(nodes_.begin(), nodes_.end(), comparison_);
#else
  realpush(x);
#endif
}

// Remove the top node from the heap
void CbcTree::pop()
{
#ifdef JJF_ZERO
  std::pop_heap(nodes_.begin(), nodes_.end(), comparison_);
  nodes_.pop_back();
#else
  if (nodes_.size()) {
    //CbcNode* s = nodes_.front();
    realpop();
    //delete s;
  }
  assert(nodes_.size() >= 0);
#endif
}

// Test if empty *** note may be overridden
bool CbcTree::empty()
{
  return nodes_.empty();
}
// Gets best node and takes off heap
CbcNode *
CbcTree::bestNode(double cutoff)
{
  CbcNode *best = NULL;
  while (!best && nodes_.size()) {
    best = nodes_.front();
    if (best)
      assert(best->objectiveValue() != COIN_DBL_MAX && best->nodeInfo());
    if (best && best->objectiveValue() != COIN_DBL_MAX && best->nodeInfo())
      assert(best->nodeInfo()->numberBranchesLeft());
    if (best && best->objectiveValue() >= cutoff) {
      // double check in case node can change its mind!
      best->checkIsCutoff(cutoff);
    }
    if (!best || best->objectiveValue() >= cutoff) {
#ifdef JJF_ZERO
      // take off
      std::pop_heap(nodes_.begin(), nodes_.end(), comparison_);
      nodes_.pop_back();
      delete best;
      best = NULL;
#else
      // let code get rid of it
      assert(best);
#endif
    }
  }
  // switched off for now
  if (best && comparison_.test_->fullScan() && false) {
    CbcNode *saveBest = best;
    int n = nodes_.size();
    int iBest = -1;
    for (int i = 0; i < n; i++) {
      // temp
      assert(nodes_[i]);
      assert(nodes_[i]->nodeInfo());
      if (nodes_[i] && nodes_[i]->objectiveValue() != COIN_DBL_MAX && nodes_[i]->nodeInfo())
        assert(nodes_[i]->nodeInfo()->numberBranchesLeft());
      if (nodes_[i] && nodes_[i]->objectiveValue() < cutoff
        && comparison_.alternateTest(best, nodes_[i])) {
        best = nodes_[i];
        iBest = i;
      }
    }
    if (best == saveBest) {
      // can pop
      // take off
      std::pop_heap(nodes_.begin(), nodes_.end(), comparison_);
      nodes_.pop_back();
    } else {
      // make impossible
      nodes_[iBest] = NULL;
    }
  } else if (best) {
    // take off
#ifdef JJF_ZERO
    std::pop_heap(nodes_.begin(), nodes_.end(), comparison_);
    nodes_.pop_back();
#else
    realpop();
#endif
  }
#ifdef DEBUG_CBC_HEAP
  if (best) {
    int n = nodes_.size();
    bool good = true;
    for (int i = 0; i < n; i++) {
      // temp
      assert(nodes_[i]);
      if (!comparison_.compareNodes(nodes_[i], best)) {
        good = false;
        CbcNode *x = nodes_[i];
        printf("i=%d x is better nun %d depth %d obj %g, best nun %d depth %d obj %g\n", i,
          x->numberUnsatisfied(), x->depth(), x->objectiveValue(),
          best->numberUnsatisfied(), best->depth(), best->objectiveValue());
      }
    }
    if (!good) {
      // compare best to all
      int i;
      for (i = 0; i < n; i++) {
        CbcNode *x = nodes_[i];
        printf("i=%d x is nun %d depth %d obj %g", i,
          x->numberUnsatisfied(), x->depth(), x->objectiveValue());
        if (!comparison_.compareNodes(x, best)) {
          printf(" - best is worse!\n");
        } else {
          printf("\n");
        }
      }
      // Now compare amongst rest
      for (i = 0; i < n; i++) {
        CbcNode *x = nodes_[i];
        printf("For i=%d ", i);
        for (int j = i + 1; j < n; j++) {
          CbcNode *y = nodes_[j];
          if (!comparison_.compareNodes(x, y)) {
            printf(" b %d", j);
          } else {
            printf(" w %d", j);
          }
        }
        printf("\n");
      }
      assert(good);
    }
  }
#endif
  if (best)
    best->setOnTree(false);
  return best;
}

/*! \brief Prune the tree using an objective function cutoff

  This routine removes all nodes with objective worst than the
  specified cutoff value.
*/

void CbcTree::cleanTree(CbcModel *model, double cutoff, double &bestPossibleObjective)
{
  int j;
  int nNodes = nodes_.size();
  CbcNode **nodeArray = new CbcNode *[nNodes];
  int *depth = new int[nNodes];
  int k = 0;
  int kDelete = nNodes;
  bestPossibleObjective = 1.0e100;
  /*
        Destructively scan the heap. Nodes to be retained go into the front of
        nodeArray, nodes to be deleted into the back. Store the depth in a
        correlated array for nodes to be deleted.
    */
  for (j = 0; j < nNodes; j++) {
    CbcNode *node = top();
    pop();
    double value = node ? node->objectiveValue() : COIN_DBL_MAX;
    if (node && value >= cutoff) {
      // double check in case node can change its mind!
      value = node->checkIsCutoff(cutoff);
    }
    bestPossibleObjective = CoinMin(bestPossibleObjective, value);
    if (value >= cutoff) {
      if (node) {
        nodeArray[--kDelete] = node;
        depth[kDelete] = node->depth();
      }
    } else {
      nodeArray[k++] = node;
    }
  }
#ifdef CBC_THREAD
  if (model->parallelMode() > 0 && model->master()) {
    // need to adjust for ones not on tree
    CbcBaseModel *master = model->master();
    int numberThreads = master->numberThreads();
    for (int i = 0; i < numberThreads; i++) {
      CbcThread *child = master->child(i);
      if (child->node()) {
        double value = child->node()->objectiveValue();
        // adjust
        bestPossibleObjective = CoinMin(bestPossibleObjective, value);
      }
    }
  }
#endif
  /*
      Rebuild the heap using the retained nodes.
    */
  for (j = 0; j < k; j++) {
    push(nodeArray[j]);
  }
  /*
      Sort the list of nodes to be deleted, nondecreasing.
    */
  CoinSort_2(depth + kDelete, depth + nNodes, nodeArray + kDelete);
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
      we're just interested in correcting the reference count. Tight/loose
      should make no difference.

      Arguably a separate routine should be used in place of addCuts1. It's
      doing more work than needed, modifying the model to match a subproblem
      at a node that will be discarded.  Then again, we seem to need the basis.
    */
  for (j = nNodes - 1; j >= kDelete; j--) {
    CbcNode *node = nodeArray[j];
    CoinWarmStartBasis *lastws = model->getEmptyBasis();

    model->addCuts1(node, lastws);
    // Decrement cut counts
    assert(node);
    //assert (node->nodeInfo());
    int numberLeft = (node->nodeInfo()) ? node->nodeInfo()->numberBranchesLeft() : 0;
    int i;
    for (i = 0; i < model->currentNumberCuts(); i++) {
      // take off node
      CoinWarmStartBasis::Status status = lastws->getArtifStatus(i + model->numberRowsAtContinuous());
      if (status != CoinWarmStartBasis::basic && model->addedCuts()[i]) {
        if (!model->addedCuts()[i]->decrement(numberLeft))
          delete model->addedCuts()[i];
      }
    }
    // node should not have anything pointing to it
    if (node->nodeInfo())
      node->nodeInfo()->throwAway();
    delete node;
    delete lastws;
  }
  delete[] nodeArray;
  delete[] depth;
}

// Return the best node of the heap using alternate criterion
CbcNode *
CbcTree::bestAlternate()
{
  int n = nodes_.size();
  CbcNode *best = NULL;
  if (n) {
    best = nodes_[0];
    for (int i = 1; i < n; i++) {
      if (comparison_.alternateTest(best, nodes_[i])) {
        best = nodes_[i];
      }
    }
  }
  return best;
}
void CbcTree::realpop()
{
  if (nodes_.size() > 0) {
    nodes_[0] = nodes_.back();
    nodes_.pop_back();
    fixTop();
  }
  assert(nodes_.size() >= 0);
}
/* After changing data in the top node, fix the heap */
void CbcTree::fixTop()
{
  const int size = nodes_.size();
  if (size > 1) {
    CbcNode **candidates = &nodes_[0];
    CbcNode *s = candidates[0];
    --candidates;
    int pos = 1;
    int ch;
    for (ch = 2; ch < size; pos = ch, ch *= 2) {
      if (!comparison_.compareNodes(candidates[ch + 1], candidates[ch]))
        ++ch;
      if (!comparison_.compareNodes(s, candidates[ch]))
        break;
      candidates[pos] = candidates[ch];
    }
    if (ch == size) {
      if (!comparison_.compareNodes(candidates[ch], s)) {
        candidates[pos] = candidates[ch];
        pos = ch;
      }
    }
    candidates[pos] = s;
  }
}
void CbcTree::realpush(CbcNode *node)
{
  node->setOnTree(true);
  nodes_.push_back(node);
  CbcNode **candidates = &nodes_[0];
  --candidates;
  int pos = nodes_.size();
  int ch;
  for (ch = pos / 2; ch != 0; pos = ch, ch /= 2) {
    if (!comparison_.compareNodes(candidates[ch], node))
      break;
    candidates[pos] = candidates[ch];
  }
  candidates[pos] = node;
}
#endif

double
CbcTree::getBestPossibleObjective()
{
  double r_val = 1e100;
  for (int i = 0; i < static_cast< int >(nodes_.size()); i++) {
    if (nodes_[i] && nodes_[i]->objectiveValue() < r_val) {
      r_val = nodes_[i]->objectiveValue();
    }
  }
  return r_val;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
