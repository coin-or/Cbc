// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcHeuristicDiveVectorLength_H
#define CbcHeuristicDiveVectorLength_H

#include "CbcHeuristicDive.hpp"

/** DiveVectorLength class
 */

class CbcHeuristicDiveVectorLength : public CbcHeuristicDive {
public:

  // Default Constructor 
  CbcHeuristicDiveVectorLength ();

  // Constructor with model - assumed before cuts
  CbcHeuristicDiveVectorLength (CbcModel & model);
  
  // Copy constructor 
  CbcHeuristicDiveVectorLength ( const CbcHeuristicDiveVectorLength &);
   
  // Destructor 
  ~CbcHeuristicDiveVectorLength ();

  /// Clone
  virtual CbcHeuristicDiveVectorLength * clone() const;
  
  /// Assignment operator 
  CbcHeuristicDiveVectorLength & operator=(const CbcHeuristicDiveVectorLength& rhs);

  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) ;

  /// Selects the next variable to branch on
  virtual void selectVariableToBranch(OsiSolverInterface* solver,
				      const double* newSolution,
				      int& bestColumn,
				      int& bestRound);

};

#endif
