// $Id$
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// edwin 12/5/09 carved out of CbcHeuristicRINS

#ifndef CbcHeuristicVND_H
#define CbcHeuristicVND_H

#include "CbcHeuristic.hpp"

/** LocalSearch class
 */

class CbcHeuristicVND : public CbcHeuristic {
public:
  // Default Constructor
  CbcHeuristicVND();

  /* Constructor with model - assumed before cuts
       Initial version does not do Lps
    */
  CbcHeuristicVND(CbcModel &model);

  // Copy constructor
  CbcHeuristicVND(const CbcHeuristicVND &);

  // Destructor
  ~CbcHeuristicVND();

  /// Clone
  virtual CbcHeuristic *clone() const;

  /// Assignment operator
  CbcHeuristicVND &operator=(const CbcHeuristicVND &rhs);

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
  /// base solution array so we can set
  inline double *baseSolution() const
  {
    return baseSolution_;
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
  /// Node when last done
  int lastNode_;
  /// Step size for decomposition
  int stepSize_;
  int k_;
  int kmax_;
  int nDifferent_;
  /// Base solution
  double *baseSolution_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
