// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcHeuristicDiveFractional_H
#define CbcHeuristicDiveFractional_H

#include "CbcHeuristic.hpp"

/** DiveFractional class
 */

class CbcHeuristicDiveFractional : public CbcHeuristic {
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

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel * model);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel * model);
  
  using CbcHeuristic::solution ;
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

  /// Set percentage of integer variables to fix at bounds
  void setPercentageToFix(double value)
  { percentageToFix_ = value; }

  /// Set maximum number of iterations
  void setMaxIterations(int value)
  { maxIterations_ = value; }

  /// Set maximum time allowed
  void setMaxTime(double value)
  { maxTime_ = value; }

protected:
  // Data

  // Original matrix by column
  CoinPackedMatrix matrix_;

  // Down locks
  unsigned short * downLocks_;

  // Up locks
  unsigned short * upLocks_;

  // Percentage of integer variables to fix at bounds
  double percentageToFix_;

  // Maximum number of iterations
  int maxIterations_;

  // Maximum time allowed
  double maxTime_;

};

#endif
