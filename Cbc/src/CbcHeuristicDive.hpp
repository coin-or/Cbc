/* $Id$ */
// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcHeuristicDive_H
#define CbcHeuristicDive_H

#include "CbcHeuristic.hpp"
class CbcSubProblem;
class OsiRowCut;
struct PseudoReducedCost {
  int var;
  double pseudoRedCost;
};

/** Dive class
 */

class CbcHeuristicDive : public CbcHeuristic {
public:
  // Default Constructor
  CbcHeuristicDive();

  // Constructor with model - assumed before cuts
  CbcHeuristicDive(CbcModel &model);

  // Copy constructor
  CbcHeuristicDive(const CbcHeuristicDive &);

  // Destructor
  ~CbcHeuristicDive();

  /// Clone
  virtual CbcHeuristicDive *clone() const = 0;

  /// Assignment operator
  CbcHeuristicDive &operator=(const CbcHeuristicDive &rhs);

  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *) {}

  /// Create C++ lines to get to current state - does work for base class
  void generateCpp(FILE *fp, const char *heuristic);

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel *model);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel *model);

  //  REMLOVE using CbcHeuristic::solution ;
  /** returns 0 if no solution, 1 if valid solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value (only if good)
        This is called after cuts have been added - so can not add cuts
        This does Fractional Diving
    */
  virtual int solution(double &objectiveValue,
    double *newSolution);
  /// inner part of dive
  int solution(double &objectiveValue, int &numberNodes,
    int &numberCuts, OsiRowCut **cuts,
    CbcSubProblem **&nodes,
    double *newSolution);
  /** returns 0 if no solution, 1 if valid solution
        with better objective value than one passed in
	also returns list of nodes
        This does Fractional Diving
    */
  int fathom(CbcModel *model, int &numberNodes, CbcSubProblem **&nodes);

  /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
  virtual void validate();

  /// Sets priorities if any
  void setPriorities();

  /// Select candidate binary variables for fixing
  void selectBinaryVariables();

  /// Set percentage of integer variables to fix at bounds
  void setPercentageToFix(double value)
  {
    percentageToFix_ = value;
  }

  /// Set maximum number of iterations
  void setMaxIterations(int value)
  {
    maxIterations_ = value;
  }

  /// Set maximum number of simplex iterations
  void setMaxSimplexIterations(int value)
  {
    maxSimplexIterations_ = value;
  }
  /// Get maximum number of simplex iterations
  inline int maxSimplexIterations() const
  {
    return maxSimplexIterations_;
  }

  /// Set maximum number of simplex iterations at root node
  void setMaxSimplexIterationsAtRoot(int value)
  {
    maxSimplexIterationsAtRoot_ = value;
  }

  /// Set maximum time allowed
  void setMaxTime(double value)
  {
    maxTime_ = value;
  }

  /// Tests if the heuristic can run
  virtual bool canHeuristicRun();

  /** Selects the next variable to branch on
        Returns true if all the fractional variables can be trivially
        rounded. Returns false, if there is at least one fractional variable
        that is not trivially roundable. In this case, the bestColumn
        returned will not be trivially roundable.
    */
  virtual bool selectVariableToBranch(OsiSolverInterface *solver,
    const double *newSolution,
    int &bestColumn,
    int &bestRound)
    = 0;
  /** Initializes any data which is going to be used repeatedly
        in selectVariableToBranch */
  virtual void initializeData() {}

  /// Perform reduced cost fixing on integer variables
  int reducedCostFix(OsiSolverInterface *solver);
  /// Fix other variables at bounds
  virtual int fixOtherVariables(OsiSolverInterface *solver,
    const double *solution,
    PseudoReducedCost *candidate,
    const double *random);

protected:
  // Data

  // Original matrix by column
  CoinPackedMatrix matrix_;

  // Original matrix by
  CoinPackedMatrix matrixByRow_;

  // Down locks
  unsigned short *downLocks_;

  // Up locks
  unsigned short *upLocks_;

  /// Extra down array (number Integers long)
  double *downArray_;

  /// Extra up array (number Integers long)
  double *upArray_;

  /// Array of priorities
  typedef struct {
    unsigned int direction : 3; //  0 bit off, 1 bit (0 down first, 1 up first) 2 bit non zero don't try other way
    unsigned int priority : 29;
  } PriorityType;
  PriorityType *priority_;
  // Indexes of binary variables with 0 objective coefficient
  // and in variable bound constraints
  std::vector< int > binVarIndex_;

  // Indexes of variable bound rows for each binary variable
  std::vector< int > vbRowIndex_;

  // Percentage of integer variables to fix at bounds
  double percentageToFix_;

  // Maximum time allowed
  double maxTime_;

  // Small objective (i.e. treat zero objective as this)
  double smallObjective_;

  // Maximum number of major iterations
  int maxIterations_;

  // Maximum number of simplex iterations
  int maxSimplexIterations_;

  // Maximum number of simplex iterations at root node
  int maxSimplexIterationsAtRoot_;
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
