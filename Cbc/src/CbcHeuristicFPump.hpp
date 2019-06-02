/* $Id$ */
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcHeuristicFeasibilityPump_H
#define CbcHeuristicFeasibilityPump_H

#include "CbcHeuristic.hpp"
#include "OsiClpSolverInterface.hpp"

/** Feasibility Pump class
 */

class CbcHeuristicFPump : public CbcHeuristic {
public:
  // Default Constructor
  CbcHeuristicFPump();

  // Constructor with model - assumed before cuts
  CbcHeuristicFPump(CbcModel &model,
    double downValue = 0.5, bool roundExpensive = false);

  // Copy constructor
  CbcHeuristicFPump(const CbcHeuristicFPump &);

  // Destructor
  ~CbcHeuristicFPump();

  /// Assignment operator
  CbcHeuristicFPump &operator=(const CbcHeuristicFPump &rhs);
  /// Clone
  virtual CbcHeuristic *clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp);

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel *model);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel *model);

  using CbcHeuristic::solution;
  /** returns 0 if no solution, 1 if valid solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value (only if good)
        This is called after cuts have been added - so can not add cuts.

        It may make sense for user to call this outside Branch and Cut to
        get solution.  Or normally is just at root node.

        * new meanings for when_ - on first try then set back to 1
          11 - at end fix all integers at same bound throughout
          12 - also fix all integers staying at same internal integral value throughout
          13 - also fix all continuous variables staying at same bound throughout
          14 - also fix all continuous variables staying at same internal value throughout
          15 - as 13 but no internal integers
      And beyond that, it's apparently possible for the range to be between 21
      and 25, in which case it's reduced on entry to solution() to be between
      11 and 15 and allSlack is set to true. Then, if we're not processing
      general integers, we'll use an all-slack basis to solve ... what? Don't
      see that yet.
    */
  virtual int solution(double &objectiveValue,
    double *newSolution);
  /** If general integers then adds variables to turn into binaries round
      solution
    */
  int solutionGeneral(double &objectiveValue, double *newSolution,
    int maxAround = 1, bool fixSatisfied = false);
  /// Set maximum Time (default off) - also sets starttime to current
  void setMaximumTime(double value);
  /// Get maximum Time (default 0.0 == time limit off)
  inline double maximumTime() const
  {
    return maximumTime_;
  }
  /// Set fake cutoff (default COIN_DBL_MAX == off)
  inline void setFakeCutoff(double value)
  {
    fakeCutoff_ = value;
  }
  /// Get fake cutoff (default 0.0 == off)
  inline double fakeCutoff() const
  {
    return fakeCutoff_;
  }
  /// Set absolute increment (default 0.0 == off)
  inline void setAbsoluteIncrement(double value)
  {
    absoluteIncrement_ = value;
  }
  /// Get absolute increment (default 0.0 == off)
  inline double absoluteIncrement() const
  {
    return absoluteIncrement_;
  }
  /// Set relative increment (default 0.0 == off)
  inline void setRelativeIncrement(double value)
  {
    relativeIncrement_ = value;
  }
  /// Get relative increment (default 0.0 == off)
  inline double relativeIncrement() const
  {
    return relativeIncrement_;
  }
  /// Set default rounding (default 0.5)
  inline void setDefaultRounding(double value)
  {
    defaultRounding_ = value;
  }
  /// Get default rounding (default 0.5)
  inline double defaultRounding() const
  {
    return defaultRounding_;
  }
  /// Set initial weight (default 0.0 == off)
  inline void setInitialWeight(double value)
  {
    initialWeight_ = value;
  }
  /// Get initial weight (default 0.0 == off)
  inline double initialWeight() const
  {
    return initialWeight_;
  }
  /// Set weight factor (default 0.1)
  inline void setWeightFactor(double value)
  {
    weightFactor_ = value;
  }
  /// Get weight factor (default 0.1)
  inline double weightFactor() const
  {
    return weightFactor_;
  }
  /// Set threshold cost for using original cost - even on continuous (default infinity)
  inline void setArtificialCost(double value)
  {
    artificialCost_ = value;
  }
  /// Get threshold cost for using original cost - even on continuous (default infinity)
  inline double artificialCost() const
  {
    return artificialCost_;
  }
  /// Get iteration to size ratio
  inline double iterationRatio() const
  {
    return iterationRatio_;
  }
  /// Set iteration to size ratio
  inline void setIterationRatio(double value)
  {
    iterationRatio_ = value;
  }
  /// Set maximum passes (default 100)
  inline void setMaximumPasses(int value)
  {
    maximumPasses_ = value;
  }
  /// Get maximum passes (default 100)
  inline int maximumPasses() const
  {
    return maximumPasses_;
  }
  /// Set maximum retries (default 1)
  inline void setMaximumRetries(int value)
  {
    maximumRetries_ = value;
  }
  /// Get maximum retries (default 1)
  inline int maximumRetries() const
  {
    return maximumRetries_;
  }
  /**  Set use of multiple solutions and solves
         0 - do not reuse solves, do not accumulate integer solutions for local search
         1 - do not reuse solves, accumulate integer solutions for local search
         2 - reuse solves, do not accumulate integer solutions for local search
         3 - reuse solves, accumulate integer solutions for local search
         If we add 4 then use second form of problem (with extra rows and variables for general integers)
       At some point (date?), I added

       And then there are a few bit fields:
       4 - something about general integers
       So my (lh) guess for 4 was at least in the ballpark, but I'll have to
       rethink 8 entirely (and it may well not mean the same thing as it did
       when I added that comment.
       8 - determines whether we process general integers

       And on 090831, John added

       If we add 4 then use second form of problem (with extra rows and
       variables for general integers)
         If we add 8 then can run after initial cuts (if no solution)
    */
  inline void setAccumulate(int value)
  {
    accumulate_ = value;
  }
  /// Get accumulation option
  inline int accumulate() const
  {
    return accumulate_;
  }
  /**  Set whether to fix variables on known solution
         0 - do not fix
         1 - fix integers on reduced costs
         2 - fix integers on reduced costs but only on entry
    */
  inline void setFixOnReducedCosts(int value)
  {
    fixOnReducedCosts_ = value;
  }
  /// Get reduced cost option
  inline int fixOnReducedCosts() const
  {
    return fixOnReducedCosts_;
  }
  /**  Set reduced cost multiplier
         1.0 as normal
         <1.0 (x) - pretend gap is x* actual gap - just for fixing
    */
  inline void setReducedCostMultiplier(double value)
  {
    reducedCostMultiplier_ = value;
  }
  /// Get reduced cost multiplier
  inline double reducedCostMultiplier() const
  {
    return reducedCostMultiplier_;
  }

protected:
  // Data
  /// Start time
  double startTime_;
  /// Maximum Cpu seconds
  double maximumTime_;
  /** Fake cutoff value.
        If set then better of real cutoff and this used to add a constraint
    */
  double fakeCutoff_;
  /// If positive carry on after solution expecting gain of at least this
  double absoluteIncrement_;
  /// If positive carry on after solution expecting gain of at least this times objective
  double relativeIncrement_;
  /// Default is round up if > this
  double defaultRounding_;
  /// Initial weight for true objective
  double initialWeight_;
  /// Factor for decreasing weight
  double weightFactor_;
  /// Threshold cost for using original cost - even on continuous
  double artificialCost_;
  /** If iterationRatio >0 use instead of maximumPasses_
        test is iterations > ratio*(2*nrow+ncol) */
  double iterationRatio_;
  /**  Reduced cost multiplier
         1.0 as normal
         <1.0 (x) - pretend gap is x* actual gap - just for fixing
    */
  double reducedCostMultiplier_;
  /// Maximum number of passes
  int maximumPasses_;
  /** Maximum number of retries if we find a solution.
        If negative we clean out used array
    */
  int maximumRetries_;
  /**  Set use of multiple solutions and solves
         0 - do not reuse solves, do not accumulate integer solutions for local search
         1 - do not reuse solves, accumulate integer solutions for local search
         2 - reuse solves, do not accumulate integer solutions for local search
         3 - reuse solves, accumulate integer solutions for local search
         If we add 4 then use second form of problem (with extra rows and variables for general integers)
         If we do not accumulate solutions then no mini branch and bounds will be done
         reuse - refers to initial solve after adding in new "cut"
         If we add 8 then can run after initial cuts (if no solution)
	 16 - extra rounding
	 32,64,128 - create binaries for general 
    */
  int accumulate_;
  /**  Set whether to fix variables on known solution
         0 - do not fix
         1 - fix integers on reduced costs
         2 - fix integers on reduced costs but only on entry
    */
  int fixOnReducedCosts_;
  /// If true round to expensive
  bool roundExpensive_;

private:
  /** Rounds solution - down if < downValue
        If roundExpensive then always to more expnsive.
        returns 0 if current is solution
    */
  int rounds(OsiSolverInterface *solver, double *solution,
    /*const double * objective, */
    int numberIntegers, const int *integerVariable,
    /*char * pumpPrint,*/ int passNumber,
    /*bool roundExpensive=false,*/
    double downValue = 0.5, int *flip = 0);
  /// Does real work
  int solutionInternal(double &objectiveValue, double *newSolution);
  /* note for eagle eyed readers.
       when_ can now be exotic -
       <=10 normal
    */
};

#ifdef COIN_HAS_CLP

class CbcDisasterHandler : public OsiClpDisasterHandler {
public:
  /**@name Virtual methods that the derived classe should provide.
    */
  //@{
#ifdef JJF_ZERO
  /// Into simplex
  virtual void intoSimplex();
  /// Checks if disaster
  virtual bool check() const;
  /// saves information for next attempt
  virtual void saveInfo();
#endif
  /// Type of disaster 0 can fix, 1 abort
  virtual int typeOfDisaster();
  //@}

  /**@name Constructors, destructor */

  //@{
  /** Default constructor. */
  CbcDisasterHandler(CbcModel *model = NULL);
  /** Destructor */
  virtual ~CbcDisasterHandler();
  // Copy
  CbcDisasterHandler(const CbcDisasterHandler &);
  // Assignment
  CbcDisasterHandler &operator=(const CbcDisasterHandler &);
  /// Clone
  virtual ClpDisasterHandler *clone() const;

  //@}

  /**@name Sets/gets */

  //@{
  /** set model. */
  void setCbcModel(CbcModel *model);
  /// Get model
  inline CbcModel *cbcModel() const
  {
    return cbcModel_;
  }

  //@}

protected:
  /**@name Data members
       The data members are protected to allow access for derived classes. */
  //@{
  /// Pointer to model
  CbcModel *cbcModel_;

  //@}
};
#endif

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
