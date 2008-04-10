// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcHeuristicDiveCoefficient_H
#define CbcHeuristicDiveCoefficient_H

#include "CbcHeuristicDive.hpp"

/** DiveCoefficient class
 */

class CbcHeuristicDiveCoefficient : public CbcHeuristicDive {
public:

  // Default Constructor 
  CbcHeuristicDiveCoefficient ();

  // Constructor with model - assumed before cuts
  CbcHeuristicDiveCoefficient (CbcModel & model);
  
  // Copy constructor 
  CbcHeuristicDiveCoefficient ( const CbcHeuristicDiveCoefficient &);
   
  // Destructor 
  ~CbcHeuristicDiveCoefficient ();

  /// Clone
  virtual CbcHeuristicDiveCoefficient * clone() const;
  
  /// Assignment operator 
  CbcHeuristicDiveCoefficient & operator=(const CbcHeuristicDiveCoefficient& rhs);

  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) ;

  /// Selects the next variable to branch on
  virtual void selectVariableToBranch(OsiSolverInterface* solver,
				      const double* newSolution,
				      int& bestColumn,
				      int& bestRound);

};

#endif
