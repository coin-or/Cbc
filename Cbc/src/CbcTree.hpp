// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcTree_H
#define CbcTree_H
using namespace std;


#include <vector>

/*! \class tree
    \brief Implementation of live set as a heap.

    This class is used to hold the set of live nodes in the search tree.
*/

class CbcTree {

public:

  // Default Constructor 
  CbcTree ();

  // Copy constructor 
  CbcTree ( const CbcTree & rhs);
  // = operator
  CbcTree & operator=(const CbcTree & rhs);
   
  virtual ~CbcTree();

  /// Clone
  virtual CbcTree * clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) {};

/*! \name Heap access and maintenance methods */
//@{

  /// Set comparison function and resort heap
  void setComparison(CbcCompareBase  &compare);

  /// Return the top node of the heap 
  virtual CbcNode * top() const;

  /// Add a node to the heap
  virtual void push(CbcNode * x);

  /// Remove the top node from the heap
  virtual void pop() ;
  /// Gets best node and takes off heap
  virtual CbcNode * bestNode(double cutoff);

//@}
/*! \name vector methods */
//@{

  /// Test if empty *** note may be overridden
  virtual bool empty() ;

  /// Return size
  inline int size() const
  { return nodes_.size();}

  /// [] operator
  inline CbcNode * operator [] (int i) const
  { return nodes_[i];}

  /// Return a node pointer
  inline CbcNode * nodePointer (int i) const
  { return nodes_[i];}

//@}

/*! \name Search tree maintenance */
//@{

/*! \brief Prune the tree using an objective function cutoff

  This routine removes all nodes with objective worst than the
  specified cutoff value.
  It also sets bestPossibleObjective to best
  of all on tree before deleting.
*/

  void cleanTree(CbcModel * model, double cutoff, double & bestPossibleObjective);

  /// Get best on list using alternate method
  CbcNode * bestAlternate();

  /// We may have got an intelligent tree so give it one more chance
  virtual void endSearch() {}
//@}
protected:
  std::vector <CbcNode *> nodes_;
  CbcCompare comparison_;	///> Sort function for heap ordering.


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
  CbcNewTree ();

  // Copy constructor 
  CbcNewTree ( const CbcNewTree & rhs);
  // = operator
  CbcNewTree & operator=(const CbcNewTree & rhs);
   
  virtual ~CbcNewTree();

  /// Clone
  virtual CbcNewTree * clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) {};

/*! \name Heap access and maintenance methods */
//@{

  /// Set comparison function and resort heap
  void setComparison(CbcCompareBase  &compare);

  /// Return the top node of the heap 
  virtual CbcNode * top() const;

  /// Add a node to the heap
  virtual void push(CbcNode * x);

  /// Remove the top node from the heap
  virtual void pop() ;
  /// Gets best node and takes off heap
  virtual CbcNode * bestNode(double cutoff);

//@}
/*! \name vector methods */
//@{

  /// Test if empty *** note may be overridden
  virtual bool empty() ;

  /// Return size
  inline int size() const
  { return nodes_.size();}

  /// [] operator
  inline CbcNode * operator [] (int i) const
  { return nodes_[i];}

  /// Return a node pointer
  inline CbcNode * nodePointer (int i) const
  { return nodes_[i];}

//@}

/*! \name Search tree maintenance */
//@{

/*! \brief Prune the tree using an objective function cutoff

  This routine removes all nodes with objective worst than the
  specified cutoff value.
  It also sets bestPossibleObjective to best
  of all on tree before deleting.
*/

  void cleanTree(CbcModel * model, double cutoff, double & bestPossibleObjective);

  /// Get best on list using alternate method
  CbcNode * bestAlternate();

  /// We may have got an intelligent tree so give it one more chance
  virtual void endSearch() {}
//@}
protected:


};
#endif

