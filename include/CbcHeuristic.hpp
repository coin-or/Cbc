// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcHeuristic_H
#define CbcHeuristic_H

#include <string>
#include <vector>
#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"

class OsiSolverInterface;

class CbcModel;

//#############################################################################
/** Heuristic base class */

class CbcHeuristic {
public:
  // Default Constructor 
  CbcHeuristic ();

  // Constructor with model - assumed before cuts
  CbcHeuristic (CbcModel & model);

  virtual ~CbcHeuristic();

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel * model);
  
  /// Clone
  virtual CbcHeuristic * clone() const=0;

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel * model)=0;

  /** returns 0 if no solution, 1 if valid solution
      with better objective value than one passed in
      Sets solution values if good, sets objective value 
      This is called after cuts have been added - so can not add cuts
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution)=0;

  /** returns 0 if no solution, 1 if valid solution, -1 if just
      returning an estimate of best possible solution
      with better objective value than one passed in
      Sets solution values if good, sets objective value (only if nonzero code)
      This is called at same time as cut generators - so can add cuts
      Default is do nothing
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution,
		       OsiCuts & cs) {return 0;};

  /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
  virtual void validate() {};

  /** Sets "when" flag - 0 off, 1 at root, 2 other than root, 3 always.
      If 10 added then don't worry if validate says there are funny objects
      as user knows it will be fine
  */
  inline void setWhen(int value)
  { when_=value;};
  /// Gets "when" flag - 0 off, 1 at root, 2 other than root, 3 always
  inline int when() const
  { return when_;};

protected:

  /// Model
  CbcModel * model_;
  /// When flag - 0 off, 1 at root, 2 other than root, 3 always
  int when_;
private:
  
  /// Illegal Assignment operator 
  CbcHeuristic & operator=(const CbcHeuristic& rhs);
  
};

/** Rounding class
 */

class CbcRounding : public CbcHeuristic {
public:

  // Default Constructor 
  CbcRounding ();

  // Constructor with model - assumed before cuts
  CbcRounding (CbcModel & model);
  
  // Copy constructor 
  CbcRounding ( const CbcRounding &);
   
  // Destructor 
  ~CbcRounding ();
  
  /// Clone
  virtual CbcHeuristic * clone() const;

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel * model);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel * model);
  
  /** returns 0 if no solution, 1 if valid solution
      with better objective value than one passed in
      Sets solution values if good, sets objective value (only if good)
      This is called after cuts have been added - so can not add cuts
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution);
  /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
  virtual void validate();


  /// Set seed
  void setSeed(int value)
  { seed_ = value;};

protected:
  // Data

  // Original matrix by column
  CoinPackedMatrix matrix_;

  // Original matrix by 
  CoinPackedMatrix matrixByRow_;

  // Seed for random stuff
  int seed_;

private:
  /// Illegal Assignment operator 
  CbcRounding & operator=(const CbcRounding& rhs);
};


#endif
