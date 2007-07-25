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

  // Copy constructor 
  CbcHeuristic ( const CbcHeuristic &);
   
  virtual ~CbcHeuristic();

  /// Clone
  virtual CbcHeuristic * clone() const=0;

  /// Assignment operator 
  CbcHeuristic & operator=(const CbcHeuristic& rhs);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel * model);
  
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
		       OsiCuts & cs) {return 0;}

  /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
  virtual void validate() {}

  /** Sets "when" flag - 0 off, 1 at root, 2 other than root, 3 always.
      If 10 added then don't worry if validate says there are funny objects
      as user knows it will be fine
  */
  inline void setWhen(int value)
  { when_=value;}
  /// Gets "when" flag - 0 off, 1 at root, 2 other than root, 3 always
  inline int when() const
  { return when_;}

  /// Sets number of nodes in subtree (default 200)
  inline void setNumberNodes(int value)
  { numberNodes_=value;}
  /// Gets number of nodes in a subtree (default 200)
  inline int numberNodes() const
  { return numberNodes_;}
  /// Just set model - do not do anything else
  inline void setModelOnly(CbcModel * model)
  { model_ = model;}
  

  /// Sets fraction of new(rows+columns)/old(rows+columns) before doing small branch and bound (default 1.0)
  inline void setFractionSmall(double value)
  { fractionSmall_=value;}
  /// Gets fraction of new(rows+columns)/old(rows+columns) before doing small branch and bound (default 1.0)
  inline double fractionSmall() const
  { return fractionSmall_;}

  /** Do mini branch and bound - return 
      0 not finished - no solution
      1 not finished - solution
      2 finished - no solution
      3 finished - solution
      (could add global cut if finished)
  */
  int smallBranchAndBound(OsiSolverInterface * solver,int numberNodes,
                          double * newSolution, double & newSolutionValue,
                          double cutoff , std::string name) const;
  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) {}
  /// Create C++ lines to get to current state - does work for base class
  void generateCpp( FILE * fp,const char * heuristic) ;
  /// Returns true if can deal with "odd" problems e.g. sos type 2
  virtual bool canDealWithOdd() const
  { return false;}
  /// return name of heuristic
  inline const char *heuristicName() const
  { return heuristicName_.c_str();}
  /// set name of heuristic
  inline void setHeuristicName(const char *name)
  { heuristicName_ = name;}

protected:

  /// Model
  CbcModel * model_;
  /// When flag - 0 off, 1 at root, 2 other than root, 3 always
  int when_;
  /// Number of nodes in any sub tree
  int numberNodes_;
  /// Fraction of new(rows+columns)/old(rows+columns) before doing small branch and bound
  double fractionSmall_;
  /// Name for printing
  std::string heuristicName_;
  
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
  
  /// Assignment operator 
  CbcRounding & operator=(const CbcRounding& rhs);

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
      This is called after cuts have been added - so can not add cuts
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution);
  /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
  virtual void validate();


  /// Set seed
  void setSeed(int value)
  { seed_ = value;}

protected:
  // Data

  // Original matrix by column
  CoinPackedMatrix matrix_;

  // Original matrix by 
  CoinPackedMatrix matrixByRow_;

  // Seed for random stuff
  int seed_;
};

/** heuristic - just picks up any good solution
    found by solver - see OsiBabSolver
 */

class CbcSerendipity : public CbcHeuristic {
public:

  // Default Constructor 
  CbcSerendipity ();

  /* Constructor with model
  */
  CbcSerendipity (CbcModel & model);
  
  // Copy constructor 
  CbcSerendipity ( const CbcSerendipity &);
   
  // Destructor 
  ~CbcSerendipity ();
  
  /// Assignment operator 
  CbcSerendipity & operator=(const CbcSerendipity& rhs);

  /// Clone
  virtual CbcHeuristic * clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) ;

  /// update model
  virtual void setModel(CbcModel * model);
  
  /** returns 0 if no solution, 1 if valid solution.
      Sets solution values if good, sets objective value (only if good)
      We leave all variables which are at one at this node of the
      tree to that value and will
      initially set all others to zero.  We then sort all variables in order of their cost
      divided by the number of entries in rows which are not yet covered.  We randomize that
      value a bit so that ties will be broken in different ways on different runs of the heuristic.
      We then choose the best one and set it to one and repeat the exercise.  

  */
  virtual int solution(double & objectiveValue,
		       double * newSolution);
  /// Resets stuff if model changes
  virtual void resetModel(CbcModel * model);

protected:
};

#endif
