/* $Id$ */
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcHeuristicRINS_H
#define CbcHeuristicRINS_H

#include "CbcHeuristic.hpp"
// for backward compatibility include 3 other headers
#include "CbcHeuristicRENS.hpp"
#include "CbcHeuristicDINS.hpp"
#include "CbcHeuristicVND.hpp"
/** LocalSearch class
 */

class CbcHeuristicRINS : public CbcHeuristic {
public:
  // Default Constructor
  CbcHeuristicRINS();

  /* Constructor with model - assumed before cuts
       Initial version does not do Lps
    */
  CbcHeuristicRINS(CbcModel &model);

  // Copy constructor
  CbcHeuristicRINS(const CbcHeuristicRINS &);

  // Destructor
  ~CbcHeuristicRINS();

  /// Clone
  virtual CbcHeuristic *clone() const;

  /// Assignment operator
  CbcHeuristicRINS &operator=(const CbcHeuristicRINS &rhs);

  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp);

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel *model);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel *model);

  using CbcHeuristic::solution;
  /** returns 0 if no solution, 1 if valid solution.
        Sets solution values if good, sets objective value (only if good)
        This does Relaxation Induced Neighborhood Search
    */
  virtual int solution(double &objectiveValue,
    double *newSolution);
  /// This version fixes stuff and does IP
  int solutionFix(double &objectiveValue,
    double *newSolution,
    const int *keep);

  /// Sets how often to do it
  inline void setHowOften(int value)
  {
    howOften_ = value;
  }
  /// Used array so we can set
  inline char *used() const
  {
    return used_;
  }
  /// Resets lastNode
  inline void setLastNode(int value)
  {
    lastNode_ = value;
  }
  /// Resets number of solutions
  inline void setSolutionCount(int value)
  {
    numberSolutions_ = value;
  }

protected:
  // Data

  /// Number of solutions so we can do something at solution
  int numberSolutions_;
  /// How often to do (code can change)
  int howOften_;
  /// Number of successes
  int numberSuccesses_;
  /// Number of tries
  int numberTries_;
  /** State of fixing continuous variables -
        0 - not tried
        +n - this divisor makes small enough
        -n - this divisor still not small enough
    */
  int stateOfFixing_;
  /// Node when last done
  int lastNode_;
  /// Whether a variable has been in a solution
  char *used_;
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
