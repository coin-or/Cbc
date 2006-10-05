// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CoinChooseVariable_H
#define CoinChooseVariable_H

#include <string>
#include <vector>

#include "CoinWarmStartBasis.hpp"
#include "OsiBranchingObject.hpp"

class OsiSolverInterface;

/** This class chooses a variable to branch on

    The base class just chooses the variable and direction without strong branching but it 
    has information which would normally be used by strong branching e.g. to re-enter
    having fixed a variable but using same candidates for strong branching.

    The flow is :
    a) initialize the process.  This clears all information.and sets state_ to 0
    b) if we are at state 0 then decide on strong branching list and set state to 1
    c) if state is 1 then do strong branching on list.  If list is empty then just
       choose aone candidate and return without strong branching.  If not empty then
       go through list and return best.  However we may find that the node is infeasible
       or that we can fix a variable.  If so we return and it is up to user to call
       again (after fixing a variable).
    d) when finished we clear the state
*/

class CoinChooseVariable  {
 
public:
    
  /// Default Constructor 
  CoinChooseVariable ();

  /// Constructor from solver (so we can set up arrays etc)
  CoinChooseVariable (const OsiSolverInterface * solver);

  /// Copy constructor 
  CoinChooseVariable (const CoinChooseVariable &);
   
  /// Assignment operator 
  CoinChooseVariable & operator= (const CoinChooseVariable& rhs);

  /// Clone
  virtual CoinChooseVariable * clone() const;

  /// Destructor 
  virtual ~CoinChooseVariable ();

  /// Initialize
  void initialize ( OsiBranchingInformation *info=NULL);
  /** Choose a variable
      Returns - 
     -1 Node is infeasible
     0  Normal termination - we have a candidate
     1  All looks satisfied - no candidate
     2  We can change the bound on a variable - but we also have a strong branching candidate
     3  We can change the bound on a variable - but we have a non-strong branching candidate
     4  We can change the bound on a variable - no other candidates
     We can pick up branch from whichObject() and whichWay()
     We can pick up a forced branch (can change bound) from whichForcedObject() and whichForcedWay()
     If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
  */
  int chooseVariable( OsiBranchingInformation *info=NULL);
  /// Finish - deletes any solution etc
  void finalize();
  /// Objective value for feasible solution
  inline double goodObjectiveValue() const
  { return goodObjectiveValue_;};
  /// Good solution - deleted by finalize
  inline const double * goodSolution() const
  { return goodSolution_;};
  /// Get the number of objects unsatisfied at this node.
  inline int numberUnsatisfied() const
  {return numberUnsatisfied_;};
  /// Number of objects to choose for strong branching
  inline int numberStrong() const
  { return numberStrong_;};
  /// Set number of objects to choose for strong branching
  inline void setNumberStrong(int value)
  { numberStrong_ = value;};
  /// Number left on strong list
  inline int numberOnList() const
  { return numberOnList_;};
  /// Trust results from strong branching for changing bounds
  inline bool trustStrongForBound() const
  { return trustStrongForBound_;};
  /// Set trust results from strong branching for changing bounds
  inline void setTrustStrongForBound(bool yesNo)
  { trustStrongForBound_ = yesNo;};
  /// Trust results from strong branching for valid solution
  inline bool trustStrongForSolution() const
  { return trustStrongForSolution_;};
  /// Set trust results from strong branching for valid solution
  inline void setTrustStrongForSolution(bool yesNo)
  { trustStrongForSolution_ = yesNo;};
  /// Set solver
  inline void setSolver (OsiSolverInterface * solver)
  { solver_ = solver;};


protected:
  // Data
  /// Objective value for feasible solution
  double goodObjectiveValue_;
  /// Good solution - deleted by finalize
  double * goodSolution_;
  /// List of candidates
  int * list_;
  /// Useful array (for sorting etc)
  double * useful_;
  /// Pointer to solver
  const OsiSolverInterface * solver_;
  /// State 0 - no list, 1 a list
  int state_;
  /* Status -
     -1 Node is infeasible
     0  Normal termination - we have a candidate
     1  All looks satisfied - no candidate
     2  We can change the bound on a variable - but we also have a strong branching candidate
     3  We can change the bound on a variable - but we have a non-strong branching candidate
     4  We can change the bound on a variable - no other candidates
  */
  /// The number of objects unsatisfied at this node.
  int numberUnsatisfied_;
  /// Number of objects to choose for strong branching
  int numberStrong_;
  /// Number left on strong list
  int numberOnList_;
  /// List of unsatisfied objects - first numberOnList_ for strong branching
  /// Trust results from strong branching for changing bounds
  bool trustStrongForBound_;
  /// Trust results from strong branching for valid solution
  bool trustStrongForSolution_;
};


#endif
