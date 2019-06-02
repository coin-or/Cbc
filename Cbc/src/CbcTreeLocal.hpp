/* $Id$ */
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcTreeLocal_H
#define CbcTreeLocal_H

//#############################################################################
/*  This implements (approximately) local branching as in the 2002 paper by
    Matteo Fischetti and Andrea Lodi.

    The very simple version of the algorithm for problems with
    0-1 variables and continuous is as follows:

    Obtain a feasible solution (one can be passed in).

    Add a cut which limits search to a k neighborhood of this solution.
    (At most k 0-1 variables may change value)
    Do branch and bound on this problem.

    If finished search and proven optimal then we can reverse cut so
    any solutions must be at least k+1 away from solution and we can
    add a new cut limiting search to a k neighborhood of new solution
    repeat.

    If finished search and no new solution then the simplest version
    would reverse last cut and complete search.  The version implemented
    here can use time and node limits and can widen search (increase effective k)
    .... and more

*/

#include "CbcTree.hpp"
#include "CbcNode.hpp"
#include "OsiRowCut.hpp"
class CbcModel;

class CbcTreeLocal : public CbcTree {

public:
  // Default Constructor
  CbcTreeLocal();

  /* Constructor with solution.
       If solution NULL no solution, otherwise must be integer
       range is initial upper bound (k) on difference from given solution.
       typeCuts -
                0 means just 0-1 cuts and will need to refine 0-1 solution
            1 uses weaker cuts on all integer variables
       maxDiversification is maximum number of range widenings to try
       timeLimit is seconds in subTree
       nodeLimit is nodes in subTree
       refine is whether to see if we can prove current solution is optimal
       when we fix all 0-1 (in case typeCuts==0 and there are general integer variables)
       if false then no refinement but reverse cuts weaker
    */
  CbcTreeLocal(CbcModel *model, const double *solution, int range = 10,
    int typeCuts = 0, int maxDiversification = 0,
    int timeLimit = 1000000, int nodeLimit = 1000000, bool refine = true);
  // Copy constructor
  CbcTreeLocal(const CbcTreeLocal &rhs);

  // = operator
  CbcTreeLocal &operator=(const CbcTreeLocal &rhs);

  virtual ~CbcTreeLocal();

  /// Clone
  virtual CbcTree *clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp);

  /*! \name Heap access and maintenance methods */
  //@{

  /// Return the top node of the heap
  virtual CbcNode *top() const;

  /// Add a node to the heap
  virtual void push(CbcNode *x);

  /// Remove the top node from the heap
  virtual void pop();

  //@}
  /*! \name Other stuff */
  //@{

  /// Create cut - return -1 if bad, 0 if okay and 1 if cut is everything
  int createCut(const double *solution, OsiRowCut &cut);

  /// Test if empty *** note may be overridden
  virtual bool empty();

  /// We may have got an intelligent tree so give it one more chance
  virtual void endSearch();
  /// Other side of last cut branch (if bias==rhs_ will be weakest possible)
  void reverseCut(int state, double bias = 0.0);
  /// Delete last cut branch
  void deleteCut(OsiRowCut &cut);
  /// Pass in solution (so can be used after heuristic)
  void passInSolution(const double *solution, double solutionValue);
  // range i.e. k
  inline int range() const
  {
    return range_;
  }
  // setrange i.e. k
  inline void setRange(int value)
  {
    range_ = value;
  }
  // Type of cuts - 0=just 0-1, 1=all
  inline int typeCuts() const
  {
    return typeCuts_;
  }
  // Type of cuts - 0=just 0-1, 1=all
  inline void setTypeCuts(int value)
  {
    typeCuts_ = value;
  }
  // maximum number of diversifications
  inline int maxDiversification() const
  {
    return maxDiversification_;
  }
  // maximum number of diversifications
  inline void setMaxDiversification(int value)
  {
    maxDiversification_ = value;
  }
  // time limit per subtree
  inline int timeLimit() const
  {
    return timeLimit_;
  }
  // time limit per subtree
  inline void setTimeLimit(int value)
  {
    timeLimit_ = value;
  }
  // node limit for subtree
  inline int nodeLimit() const
  {
    return nodeLimit_;
  }
  // node limit for subtree
  inline void setNodeLimit(int value)
  {
    nodeLimit_ = value;
  }
  // Whether to do refinement step
  inline bool refine() const
  {
    return refine_;
  }
  // Whether to do refinement step
  inline void setRefine(bool yesNo)
  {
    refine_ = yesNo;
  }

  //@}
private:
  // Node for local cuts
  CbcNode *localNode_;
  // best solution
  double *bestSolution_;
  // saved solution
  double *savedSolution_;
  // solution number at start of pass
  int saveNumberSolutions_;
  /* Cut.  If zero size then no solution yet.  Otherwise is left hand branch */
  OsiRowCut cut_;
  // This cut fixes all 0-1 variables
  OsiRowCut fixedCut_;
  // Model
  CbcModel *model_;
  // Original lower bounds
  double *originalLower_;
  // Original upper bounds
  double *originalUpper_;
  // range i.e. k
  int range_;
  // Type of cuts - 0=just 0-1, 1=all
  int typeCuts_;
  // maximum number of diversifications
  int maxDiversification_;
  // current diversification
  int diversification_;
  // Whether next will be strong diversification
  bool nextStrong_;
  // Current rhs
  double rhs_;
  // Save allowable gap
  double savedGap_;
  // Best solution
  double bestCutoff_;
  // time limit per subtree
  int timeLimit_;
  // time when subtree started
  int startTime_;
  // node limit for subtree
  int nodeLimit_;
  // node count when subtree started
  int startNode_;
  // -1 not started, 0 == stop on first solution, 1 don't stop on first, 2 refinement step
  int searchType_;
  // Whether to do refinement step
  bool refine_;
};

class CbcTreeVariable : public CbcTree {

public:
  // Default Constructor
  CbcTreeVariable();

  /* Constructor with solution.
       If solution NULL no solution, otherwise must be integer
       range is initial upper bound (k) on difference from given solution.
       typeCuts -
                0 means just 0-1 cuts and will need to refine 0-1 solution
            1 uses weaker cuts on all integer variables
       maxDiversification is maximum number of range widenings to try
       timeLimit is seconds in subTree
       nodeLimit is nodes in subTree
       refine is whether to see if we can prove current solution is optimal
       when we fix all 0-1 (in case typeCuts==0 and there are general integer variables)
       if false then no refinement but reverse cuts weaker
    */
  CbcTreeVariable(CbcModel *model, const double *solution, int range = 10,
    int typeCuts = 0, int maxDiversification = 0,
    int timeLimit = 1000000, int nodeLimit = 1000000, bool refine = true);
  // Copy constructor
  CbcTreeVariable(const CbcTreeVariable &rhs);

  // = operator
  CbcTreeVariable &operator=(const CbcTreeVariable &rhs);

  virtual ~CbcTreeVariable();

  /// Clone
  virtual CbcTree *clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp);

  /*! \name Heap access and maintenance methods */
  //@{

  /// Return the top node of the heap
  virtual CbcNode *top() const;

  /// Add a node to the heap
  virtual void push(CbcNode *x);

  /// Remove the top node from the heap
  virtual void pop();

  //@}
  /*! \name Other stuff */
  //@{

  /// Create cut - return -1 if bad, 0 if okay and 1 if cut is everything
  int createCut(const double *solution, OsiRowCut &cut);

  /// Test if empty *** note may be overridden
  virtual bool empty();

  /// We may have got an intelligent tree so give it one more chance
  virtual void endSearch();
  /// Other side of last cut branch (if bias==rhs_ will be weakest possible)
  void reverseCut(int state, double bias = 0.0);
  /// Delete last cut branch
  void deleteCut(OsiRowCut &cut);
  /// Pass in solution (so can be used after heuristic)
  void passInSolution(const double *solution, double solutionValue);
  // range i.e. k
  inline int range() const
  {
    return range_;
  }
  // setrange i.e. k
  inline void setRange(int value)
  {
    range_ = value;
  }
  // Type of cuts - 0=just 0-1, 1=all
  inline int typeCuts() const
  {
    return typeCuts_;
  }
  // Type of cuts - 0=just 0-1, 1=all
  inline void setTypeCuts(int value)
  {
    typeCuts_ = value;
  }
  // maximum number of diversifications
  inline int maxDiversification() const
  {
    return maxDiversification_;
  }
  // maximum number of diversifications
  inline void setMaxDiversification(int value)
  {
    maxDiversification_ = value;
  }
  // time limit per subtree
  inline int timeLimit() const
  {
    return timeLimit_;
  }
  // time limit per subtree
  inline void setTimeLimit(int value)
  {
    timeLimit_ = value;
  }
  // node limit for subtree
  inline int nodeLimit() const
  {
    return nodeLimit_;
  }
  // node limit for subtree
  inline void setNodeLimit(int value)
  {
    nodeLimit_ = value;
  }
  // Whether to do refinement step
  inline bool refine() const
  {
    return refine_;
  }
  // Whether to do refinement step
  inline void setRefine(bool yesNo)
  {
    refine_ = yesNo;
  }

  //@}
private:
  // Node for local cuts
  CbcNode *localNode_;
  // best solution
  double *bestSolution_;
  // saved solution
  double *savedSolution_;
  // solution number at start of pass
  int saveNumberSolutions_;
  /* Cut.  If zero size then no solution yet.  Otherwise is left hand branch */
  OsiRowCut cut_;
  // This cut fixes all 0-1 variables
  OsiRowCut fixedCut_;
  // Model
  CbcModel *model_;
  // Original lower bounds
  double *originalLower_;
  // Original upper bounds
  double *originalUpper_;
  // range i.e. k
  int range_;
  // Type of cuts - 0=just 0-1, 1=all
  int typeCuts_;
  // maximum number of diversifications
  int maxDiversification_;
  // current diversification
  int diversification_;
  // Whether next will be strong diversification
  bool nextStrong_;
  // Current rhs
  double rhs_;
  // Save allowable gap
  double savedGap_;
  // Best solution
  double bestCutoff_;
  // time limit per subtree
  int timeLimit_;
  // time when subtree started
  int startTime_;
  // node limit for subtree
  int nodeLimit_;
  // node count when subtree started
  int startNode_;
  // -1 not started, 0 == stop on first solution, 1 don't stop on first, 2 refinement step
  int searchType_;
  // Whether to do refinement step
  bool refine_;
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
