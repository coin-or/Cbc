// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcHeuristicDivePseudoCost_H
#define CbcHeuristicDivePseudoCost_H

#include "CbcHeuristicDive.hpp"

/** DivePseudoCost class
 */

class CbcHeuristicDivePseudoCost : public CbcHeuristicDive {
public:

  // Default Constructor 
  CbcHeuristicDivePseudoCost ();

  // Constructor with model - assumed before cuts
  CbcHeuristicDivePseudoCost (CbcModel & model);
  
  // Copy constructor 
  CbcHeuristicDivePseudoCost ( const CbcHeuristicDivePseudoCost &);
   
  // Destructor 
  ~CbcHeuristicDivePseudoCost ();

  /// Clone
  virtual CbcHeuristicDivePseudoCost * clone() const;
  
  /// Assignment operator 
  CbcHeuristicDivePseudoCost & operator=(const CbcHeuristicDivePseudoCost& rhs);

  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) ;

  /// Selects the next variable to branch on
  /** Returns true if all the fractional variables can be trivially
      rounded. Returns false, if there is at least one fractional variable
      that is not trivially roundable. In this case, the bestColumn
      returned will not be trivially roundable.
  */
  virtual bool selectVariableToBranch(OsiSolverInterface* solver,
				      const double* newSolution,
				      int& bestColumn,
				      int& bestRound);

};

#endif
