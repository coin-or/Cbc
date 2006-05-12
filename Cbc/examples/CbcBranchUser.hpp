// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcBranchUser_H
#define CbcBranchUser_H

#include "CbcBranchBase.hpp"

/** Branching decision user class */

class CbcBranchUserDecision : public CbcBranchDecision {
public:
  // Default Constructor 
  CbcBranchUserDecision ();

  // Copy constructor 
  CbcBranchUserDecision ( const CbcBranchUserDecision &);

  virtual ~CbcBranchUserDecision();

 /// Clone
  virtual CbcBranchDecision * clone() const;

    /// Initialize i.e. before start of choosing at a node
  virtual void initialize(CbcModel * model);

  /** Returns nonzero if branching on first object is "better" than on
      second (if second NULL first wins).
      This is only used after strong branching.  The initial selection
      is done by infeasibility() for each CbcObject
      return code +1 for up branch preferred, -1 for down
      
 */
  virtual int betterBranch(CbcBranchingObject * thisOne,
			    CbcBranchingObject * bestSoFar,
			    double changeUp, int numberInfeasibilitiesUp,
			    double changeDown, int numberInfeasibilitiesDown);

  /** \brief Compare N branching objects. Return index of best
      and sets way of branching in chosen object.
    
    This routine is used only after strong branching.
    This is reccommended version as it can be more sophisticated
  */

  virtual int
  bestBranch (CbcBranchingObject ** objects, int numberObjects, int numberUnsatisfied,
	      double * changeUp, int * numberInfeasibilitiesUp,
	      double * changeDown, int * numberInfeasibilitiesDown,
	      double objectiveValue) ;
private:
  
  /// Illegal Assignment operator 
  CbcBranchUserDecision & operator=(const CbcBranchUserDecision& rhs);

};

#endif
