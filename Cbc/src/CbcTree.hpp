/* $Id$ */
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcTree_H
#define CbcTree_H

#include <vector>
#include <algorithm>
#include <cmath>

#include "CoinHelperFunctions.hpp"
#include "CbcCompare.hpp"

/*! \brief Using MS heap implementation

  It's unclear if this is needed any longer, or even if it should be allowed.
  Cbc occasionally tries to do things to the tree (typically tweaking the
  comparison predicate) that can cause a violation of the heap property (parent better
  than either child). In a debug build, Microsoft's heap implementation does checks that
  detect this and fail. This symbol switched to an alternate implementation of CbcTree,
  and there are clearly differences, but no explanation as to why or what for.

  As of 100921, the code is cleaned up to make it through `cbc -unitTest' without
  triggering `Invalid heap' in an MSVS debug build. The method validateHeap() can
  be used for debugging if this turns up again.
*/
//#define CBC_DUBIOUS_HEAP
#if defined(_MSC_VER) || defined(__MNO_CYGWIN)
//#define CBC_DUBIOUS_HEAP
#endif
#if 1 //ndef CBC_DUBIOUS_HEAP

/*! \brief Controls search tree debugging

  In order to have validateHeap() available, set CBC_DEBUG_HEAP
  to 1 or higher.

  - 1 calls validateHeap() after each change to the heap
  - 2 will print a line for major operations (clean, set comparison, etc.)
  - 3 will print information about each push and pop

#define CBC_DEBUG_HEAP 1
*/

/*! \class CbcTree
    \brief Implementation of the live set as a heap.

    This class is used to hold the set of live nodes in the search tree.
*/
class CbcTree {

public:
  /*! \name Constructors and related */
  //@{
  /// Default Constructor
  CbcTree();

  /// Copy constructor
  CbcTree(const CbcTree &rhs);

  /// = operator
  CbcTree &operator=(const CbcTree &rhs);

  /// Destructor
  virtual ~CbcTree();

  /// Clone
  virtual CbcTree *clone() const;

  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *) {}
  //@}

  /*! \name Heap access and maintenance methods */
  //@{
  /// Set comparison function and resort heap
  void setComparison(CbcCompareBase &compare);

  /// Return the top node of the heap
  virtual CbcNode *top() const;

  /// Add a node to the heap
  virtual void push(CbcNode *x);

  /// Remove the top node from the heap
  virtual void pop();

  /*! \brief Gets best node and takes off heap

      Before returning the node from the top of the heap, the node
      is offered an opportunity to reevaluate itself. Callers should
      be prepared to check that the node returned is suitable for use.
    */
  virtual CbcNode *bestNode(double cutoff);

  /*! \brief Rebuild the heap */
  virtual void rebuild();
  //@}

  /*! \name Direct node access methods */
  //@{
  /// Test for an empty tree
  virtual bool empty();

  /// Return size
  virtual int size() const { return static_cast< int >(nodes_.size()); }

  /// Return a node pointer
  inline CbcNode *operator[](int i) const { return nodes_[i]; }

  /// Return a node pointer
  inline CbcNode *nodePointer(int i) const { return nodes_[i]; }
  void realpop();
  /** After changing data in the top node, fix the heap */
  void fixTop();
  void realpush(CbcNode *node);
  //@}

  /*! \name Search tree maintenance */
  //@{
  /*! \brief Prune the tree using an objective function cutoff

      This routine removes all nodes with objective worse than the
      specified cutoff value. It also sets bestPossibleObjective to
      the best objective over remaining nodes.
    */
  virtual void cleanTree(CbcModel *model, double cutoff, double &bestPossibleObjective);

  /// Get best on list using alternate method
  CbcNode *bestAlternate();

  /// We may have got an intelligent tree so give it one more chance
  virtual void endSearch() {}

  /// Get best possible objective function in the tree
  virtual double getBestPossibleObjective();

  /// Reset maximum node number
  inline void resetNodeNumbers() { maximumNodeNumber_ = 0; }

  /// Get maximum node number
  inline int maximumNodeNumber() const { return maximumNodeNumber_; }

  /// Set number of branches
  inline void setNumberBranching(int value) { numberBranching_ = value; }

  /// Get number of branches
  inline int getNumberBranching() const { return numberBranching_; }

  /// Set maximum branches
  inline void setMaximumBranching(int value) { maximumBranching_ = value; }

  /// Get maximum branches
  inline int getMaximumBranching() const { return maximumBranching_; }

  /// Get branched variables
  inline unsigned int *branched() const { return branched_; }

  /// Get bounds
  inline int *newBounds() const { return newBound_; }

  /// Last objective in branch-and-cut search tree
  inline double lastObjective() const
  {
    return lastObjective_;
  }
  /// Last depth in branch-and-cut search tree
  inline int lastDepth() const
  {
    return lastDepth_;
  }
  /// Last number of objects unsatisfied
  inline int lastUnsatisfied() const
  {
    return lastUnsatisfied_;
  }
  /// Adds branching information to complete state
  void addBranchingInformation(const CbcModel *model, const CbcNodeInfo *nodeInfo,
    const double *currentLower,
    const double *currentUpper);
  /// Increase space for data
  void increaseSpace();
  //@}

#if CBC_DEBUG_HEAP > 0
  /*! \name Debugging methods */
  //@{
  /*! \brief Check that the heap property is satisfied. */
  void validateHeap();
  //@}
#endif

protected:
  /// Storage vector for the heap
  std::vector< CbcNode * > nodes_;
  /// Sort predicate for heap ordering.
  CbcCompare comparison_;
  /// Maximum "node" number so far to split ties
  int maximumNodeNumber_;
  /// Size of variable list
  int numberBranching_;
  /// Maximum size of variable list
  int maximumBranching_;
  /// Objective of last node pushed on tree
  double lastObjective_;
  /// Depth of last node pushed on tree
  int lastDepth_;
  /// Number unsatisfied of last node pushed on tree
  int lastUnsatisfied_;
  /** Integer variables branched or bounded
        top bit set if new upper bound
        next bit set if a branch
    */
  unsigned int *branched_;
  /// New bound
  int *newBound_;
};

#ifdef JJF_ZERO // not used
/*! \brief Implementation of live set as a managed array.

    This class is used to hold the set of live nodes in the search tree.
*/
class CbcTreeArray : public CbcTree {

public:
  // Default Constructor
  CbcTreeArray();

  // Copy constructor
  CbcTreeArray(const CbcTreeArray &rhs);
  // = operator
  CbcTreeArray &operator=(const CbcTreeArray &rhs);

  virtual ~CbcTreeArray();

  /// Clone
  virtual CbcTree *clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *) {}

  /*! \name Heap access and maintenance methods */
  //@{

  /// Set comparison function and resort heap
  void setComparison(CbcCompareBase &compare);

  /// Add a node to the heap
  virtual void push(CbcNode *x);

  /// Gets best node and takes off heap
  virtual CbcNode *bestNode(double cutoff);

  //@}
  /*! \name vector methods */
  //@{

  /// Test if empty *** note may be overridden
  virtual bool empty();

  //@}

  /*! \name Search tree maintenance */
  //@{

  /*! \brief Prune the tree using an objective function cutoff

      This routine removes all nodes with objective worst than the
      specified cutoff value.
      It also sets bestPossibleObjective to best
      of all on tree before deleting.
    */

  void cleanTree(CbcModel *model, double cutoff, double &bestPossibleObjective);
  /// Get best possible objective function in the tree
  virtual double getBestPossibleObjective();
  //@}
protected:
  /// Returns
  /// Last node
  CbcNode *lastNode_;
  /// Last node popped
  CbcNode *lastNodePopped_;
  /// Not used yet
  int switches_;
};

/// New style
#include "CoinSearchTree.hpp"
/*! \class tree
    \brief Implementation of live set as a heap.

    This class is used to hold the set of live nodes in the search tree.
*/

class CbcNewTree : public CbcTree, public CoinSearchTreeManager {

public:
  // Default Constructor
  CbcNewTree();

  // Copy constructor
  CbcNewTree(const CbcNewTree &rhs);
  // = operator
  CbcNewTree &operator=(const CbcNewTree &rhs);

  virtual ~CbcNewTree();

  /// Clone
  virtual CbcNewTree *clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *) {}

  /*! \name Heap access and maintenance methods */
  //@{

  /// Set comparison function and resort heap
  void setComparison(CbcCompareBase &compare);

  /// Return the top node of the heap
  virtual CbcNode *top() const;

  /// Add a node to the heap
  virtual void push(CbcNode *x);

  /// Remove the top node from the heap
  virtual void pop();
  /// Gets best node and takes off heap
  virtual CbcNode *bestNode(double cutoff);

  //@}
  /*! \name vector methods */
  //@{

  /// Test if empty *** note may be overridden
  virtual bool empty();

  /// Return size
  inline int size() const
  {
    return nodes_.size();
  }

  /// [] operator
  inline CbcNode *operator[](int i) const
  {
    return nodes_[i];
  }

  /// Return a node pointer
  inline CbcNode *nodePointer(int i) const
  {
    return nodes_[i];
  }

  //@}

  /*! \name Search tree maintenance */
  //@{

  /*! \brief Prune the tree using an objective function cutoff

      This routine removes all nodes with objective worst than the
      specified cutoff value.
      It also sets bestPossibleObjective to best
      of all on tree before deleting.
    */

  void cleanTree(CbcModel *model, double cutoff, double &bestPossibleObjective);

  /// Get best on list using alternate method
  CbcNode *bestAlternate();

  /// We may have got an intelligent tree so give it one more chance
  virtual void endSearch() {}
  //@}
protected:
};
#endif
#else
/* CBC_DUBIOUS_HEAP is defined

  See note at top of file. This code is highly suspect.
  -- lh, 100921 --
*/
class CbcTree {

public:
  // Default Constructor
  CbcTree();

  // Copy constructor
  CbcTree(const CbcTree &rhs);
  // = operator
  CbcTree &operator=(const CbcTree &rhs);

  virtual ~CbcTree();

  /// Clone
  virtual CbcTree *clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp) {}

  /*! \name Heap access and maintenance methods */
  //@{

  /// Set comparison function and resort heap
  void setComparison(CbcCompareBase &compare);

  /// Return the top node of the heap
  virtual CbcNode *top() const;

  /// Add a node to the heap
  virtual void push(CbcNode *x);

  /// Remove the top node from the heap
  virtual void pop();
  /// Gets best node and takes off heap
  virtual CbcNode *bestNode(double cutoff);

  //@}
  /*! \name vector methods */
  //@{

  /// Test if empty *** note may be overridden
  //virtual bool empty() ;

  /// Return size
  inline int size() const
  {
    return nodes_.size();
  }

  /// [] operator
  inline CbcNode *operator[](int i) const
  {
    return nodes_[i];
  }

  /// Return a node pointer
  inline CbcNode *nodePointer(int i) const
  {
    return nodes_[i];
  }

  virtual bool empty();
  //inline int size() const { return size_; }
  void realpop();
  /** After changing data in the top node, fix the heap */
  void fixTop();
  void realpush(CbcNode *node);
  //@}

  /*! \name Search tree maintenance */
  //@{

  /*! \brief Prune the tree using an objective function cutoff

      This routine removes all nodes with objective worst than the
      specified cutoff value.
      It also sets bestPossibleObjective to best
      of all on tree before deleting.
    */

  void cleanTree(CbcModel *model, double cutoff, double &bestPossibleObjective);

  /// Get best on list using alternate method
  CbcNode *bestAlternate();

  /// We may have got an intelligent tree so give it one more chance
  virtual void endSearch() {}
  /// Reset maximum node number
  inline void resetNodeNumbers()
  {
    maximumNodeNumber_ = 0;
  }

  /// Get maximum node number
  inline int maximumNodeNumber() const { return maximumNodeNumber_; }
  //@}
protected:
  std::vector< CbcNode * > nodes_;
  CbcCompare comparison_; ///> Sort function for heap ordering.
  /// Maximum "node" number so far to split ties
  int maximumNodeNumber_;
};
#endif
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
