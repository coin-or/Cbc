// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcHeuristicDiveGuided_H
#define CbcHeuristicDiveGuided_H

#include "CbcHeuristicDive.hpp"

/** DiveGuided class
 */

class CbcHeuristicDiveGuided : public CbcHeuristicDive {
public:

  // Default Constructor 
  CbcHeuristicDiveGuided ();

  // Constructor with model - assumed before cuts
  CbcHeuristicDiveGuided (CbcModel & model);
  
  // Copy constructor 
  CbcHeuristicDiveGuided ( const CbcHeuristicDiveGuided &);
   
  // Destructor 
  ~CbcHeuristicDiveGuided ();

  /// Clone
  virtual CbcHeuristicDiveGuided * clone() const;
  
  /// Assignment operator 
  CbcHeuristicDiveGuided & operator=(const CbcHeuristicDiveGuided& rhs);

  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) ;

  /// Tests if the heuristic can run
  virtual bool canHeuristicRun();

  /// Selects the next variable to branch on
  virtual void selectVariableToBranch(OsiSolverInterface* solver,
				      const double* newSolution,
				      int& bestColumn,
				      int& bestRound);

};

#endif
