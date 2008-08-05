// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcHeuristicDive_H
#define CbcHeuristicDive_H

#include "CbcHeuristic.hpp"

/** Dive class
 */

class CbcHeuristicDive : public CbcHeuristic {
public:

  // Default Constructor 
  CbcHeuristicDive ();

  // Constructor with model - assumed before cuts
  CbcHeuristicDive (CbcModel & model);
  
  // Copy constructor 
  CbcHeuristicDive ( const CbcHeuristicDive &);
   
  // Destructor 
  ~CbcHeuristicDive ();

  /// Clone
  virtual CbcHeuristicDive * clone() const = 0;
  
  /// Assignment operator 
  CbcHeuristicDive & operator=(const CbcHeuristicDive& rhs);

  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) {}

  /// Create C++ lines to get to current state - does work for base class
  void generateCpp( FILE * fp,const char * heuristic);

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel * model);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel * model);
  
  //  REMLOVE using CbcHeuristic::solution ;
  /** returns 0 if no solution, 1 if valid solution
      with better objective value than one passed in
      Sets solution values if good, sets objective value (only if good)
      This is called after cuts have been added - so can not add cuts
      This does Fractional Diving
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution);

  /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
  virtual void validate();

  /// Select candidate binary variables for fixing
  void selectBinaryVariables();

  /// Set percentage of integer variables to fix at bounds
  void setPercentageToFix(double value)
  { percentageToFix_ = value; }

  /// Set maximum number of iterations
  void setMaxIterations(int value)
  { maxIterations_ = value; }

  /// Set maximum number of simplex iterations
  void setMaxSimplexIterations(int value)
  { maxSimplexIterations_ = value; }

  /// Set maximum number of simplex iterations at root node
  void setMaxSimplexIterationsAtRoot(int value)
  { maxSimplexIterationsAtRoot_ = value; }

  /// Set maximum time allowed
  void setMaxTime(double value)
  { maxTime_ = value; }

  /// Tests if the heuristic can run
  virtual bool canHeuristicRun();

  /// Selects the next variable to branch on
  /** Returns true if all the fractional variables can be trivially
      rounded. Returns false, if there is at least one fractional variable
      that is not trivially roundable. In this case, the bestColumn
      returned will not be trivially roundable.
  */
  virtual bool selectVariableToBranch(OsiSolverInterface* solver,
				      const double* newSolution,
				      int& bestColumn,
				      int& bestRound) = 0;

  /// Perform reduced cost fixing on integer variables
  int reducedCostFix (OsiSolverInterface* solver);

protected:
  // Data

  // Original matrix by column
  CoinPackedMatrix matrix_;

  // Original matrix by 
  CoinPackedMatrix matrixByRow_;

  // Down locks
  unsigned short * downLocks_;

  // Up locks
  unsigned short * upLocks_;

  // Indexes of binary variables with 0 objective coefficient
  // and in variable bound constraints
  std::vector<int> binVarIndex_;

  // Indexes of variable bound rows for each binary variable
  std::vector<int> vbRowIndex_;

  // Percentage of integer variables to fix at bounds
  double percentageToFix_;

  // Maximum number of major iterations
  int maxIterations_;

  // Maximum number of simplex iterations
  int maxSimplexIterations_;

  // Maximum number of simplex iterations at root node
  int maxSimplexIterationsAtRoot_;

  // Maximum time allowed
  double maxTime_;

};
#endif
