// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcHeuristicFeasibilityPump_H
#define CbcHeuristicFeasibilityPump_H

#include "CbcHeuristic.hpp"

/** Rounding class
 */

class CbcHeuristicFPump : public CbcHeuristic {
public:

  // Default Constructor 
  CbcHeuristicFPump ();

  // Constructor with model - assumed before cuts
  CbcHeuristicFPump (CbcModel & model,
		     double downValue=0.5,bool roundExpensive=false);
  
  // Copy constructor 
  CbcHeuristicFPump ( const CbcHeuristicFPump &);
   
  // Destructor 
  ~CbcHeuristicFPump ();
  
  /// Clone
  virtual CbcHeuristic * clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) ;

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel * model);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel * model);
  
  /** returns 0 if no solution, 1 if valid solution
      with better objective value than one passed in
      Sets solution values if good, sets objective value (only if good)
      This is called after cuts have been added - so can not add cuts.

      It may make sense for user to call this outside Branch and Cut to
      get solution.  Or normally is just at root node.
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution);

  /// Set maximum Time (default off) - also sets starttime to current
  void setMaximumTime(double value);
  /// Get maximum Time (default 0.0 == time limit off)
  inline double maximumTime() const
  { return maximumTime_;};
  /// Set fake cutoff (default COIN_DBL_MAX == off)
  inline void setFakeCutoff(double value)
  { fakeCutoff_ = value;};
  /// Get fake cutoff (default 0.0 == off)
  inline double fakeCutoff() const
  { return fakeCutoff_;};
  /// Set absolute increment (default 0.0 == off)
  inline void setAbsoluteIncrement(double value)
  { absoluteIncrement_ = value;};
  /// Get absolute increment (default 0.0 == off)
  inline double absoluteIncrement() const
  { return absoluteIncrement_;};
  /// Set relative increment (default 0.0 == off)
  inline void setRelativeIncrement(double value)
  { relativeIncrement_ = value;};
  /// Get relative increment (default 0.0 == off)
  inline double relativeIncrement() const
  { return relativeIncrement_;};
  /// Set maximum passes (default 100)
  inline void setMaximumPasses(int value)
  { maximumPasses_=value;};
  /// Get maximum passes (default 100)
  inline int maximumPasses() const
  { return maximumPasses_;};
  /// Set maximum retries (default 1)
  inline void setMaximumRetries(int value)
  { maximumRetries_=value;};
  /// Get maximum retries (default 1)
  inline int maximumRetries() const
  { return maximumRetries_;};

protected:
  // Data
  /// Start time
  double startTime_;
  /// Maximum Cpu seconds
  double maximumTime_;
  /// If less than this round down
  double downValue_;
  /** Fake cutoff value.
      If set then better of real cutoff and this used to add a constraint
  */
  double fakeCutoff_;
  /// If positive carry on after solution expecting gain of at least this
  double absoluteIncrement_;
  /// If positive carry on after solution expecting gain of at least this times objective
  double relativeIncrement_;
  /// Maximum number of passes
  int maximumPasses_;
  /** Maximum number of retries if we find a solution.
      If negative we clean out used array
  */
  int maximumRetries_;
  /// If true round to expensive
  bool roundExpensive_;

private:
  /// Illegal Assignment operator 
  CbcHeuristicFPump & operator=(const CbcHeuristicFPump& rhs);
  /** Rounds solution - down if < downValue
      If roundExpensive then always to more expnsive.
      returns 0 if current is solution
  */
  int rounds(double * solution, const double * objective, 
	     bool roundExpensive=false,
	     double downValue=0.5, int *flip=0);
  /* note for eagle eyed readers.
     when_ can now be exotic -
     <=10 normal
  */
};


#endif
