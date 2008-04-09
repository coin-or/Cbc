// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcHeuristicDiveFractional_H
#define CbcHeuristicDiveFractional_H

#include "CbcHeuristicDive.hpp"

/** DiveFractional class
 */

class CbcHeuristicDiveFractional : public CbcHeuristicDive {
public:

  // Default Constructor 
  CbcHeuristicDiveFractional ();

  // Constructor with model - assumed before cuts
  CbcHeuristicDiveFractional (CbcModel & model);
  
  // Copy constructor 
  CbcHeuristicDiveFractional ( const CbcHeuristicDiveFractional &);
   
  // Destructor 
  ~CbcHeuristicDiveFractional ();

  /// Clone
  virtual CbcHeuristicDiveFractional * clone() const;
  
  /// Assignment operator 
  CbcHeuristicDiveFractional & operator=(const CbcHeuristicDiveFractional& rhs);

  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) ;

  /// Selects the next variable to branch on
  virtual void selectVariableToBranch(OsiSolverInterface* solver,
				      const double* newSolution,
				      int& bestColumn,
				      int& bestRound);

};

#endif
