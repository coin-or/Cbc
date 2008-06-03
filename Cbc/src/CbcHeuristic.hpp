// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcHeuristic_H
#define CbcHeuristic_H

#include <string>
#include <vector>
#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"
#include "CoinHelperFunctions.hpp"
#include "OsiBranchingObject.hpp"

class OsiSolverInterface;

class CbcModel;

//#############################################################################

class CbcHeuristicNodeList;
class CbcBranchingObject;

/** A class describing the branching decisions that were made to get
    to the node where a heuristics was invoked from */

class CbcHeuristicNode {
private:
  void gutsOfConstructor(CbcModel& model);
  CbcHeuristicNode();
  CbcHeuristicNode& operator=(const CbcHeuristicNode&);
private:
  /// The number of branching decisions made
  int numObjects_;
  /** The indices of the branching objects. Note: an index may be
      listed multiple times. E.g., a general integer variable that has
      been branched on multiple times. */
  CbcBranchingObject** brObj_;
public:
  CbcHeuristicNode(CbcModel& model);

  CbcHeuristicNode(const CbcHeuristicNode& rhs);
  ~CbcHeuristicNode();
  double distance(const CbcHeuristicNode* node) const;
  double minDistance(const CbcHeuristicNodeList& nodeList) const;
  bool minDistanceIsSmall(const CbcHeuristicNodeList& nodeList,
			  const double threshold) const;
  double avgDistance(const CbcHeuristicNodeList& nodeList) const;
};

class CbcHeuristicNodeList {
private:
  void gutsOfDelete();
  void gutsOfCopy(const CbcHeuristicNodeList& rhs);
private:
  std::vector<CbcHeuristicNode*> nodes_;
public:
  CbcHeuristicNodeList() {}
  CbcHeuristicNodeList(const CbcHeuristicNodeList& rhs);
  CbcHeuristicNodeList& operator=(const CbcHeuristicNodeList& rhs);
  ~CbcHeuristicNodeList();
  
  void append(CbcHeuristicNode*& node);
  void append(const CbcHeuristicNodeList& nodes);
  inline const CbcHeuristicNode* node(int i) const { return nodes_[i]; }
  inline int size() const { return nodes_.size(); }
};

//#############################################################################
/** Heuristic base class */

class CbcHeuristic {
private:
  void gutsOfDelete() {}
  void gutsOfCopy(const CbcHeuristic & rhs);

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
  /// Sets feasibility pump options (-1 is off)
  inline void setFeasibilityPumpOptions(int value)
  { feasibilityPumpOptions_=value;}
  /// Gets feasibility pump options (-1 is off)
  inline int feasibilityPumpOptions() const
  { return feasibilityPumpOptions_;}
  /// Just set model - do not do anything else
  inline void setModelOnly(CbcModel * model)
  { model_ = model;}
  

  /// Sets fraction of new(rows+columns)/old(rows+columns) before doing small branch and bound (default 1.0)
  inline void setFractionSmall(double value)
  { fractionSmall_=value;}
  /// Gets fraction of new(rows+columns)/old(rows+columns) before doing small branch and bound (default 1.0)
  inline double fractionSmall() const
  { return fractionSmall_;}
  /// Get how many solutions the heuristic thought it got
  inline int numberSolutionsFound() const
  { return numberSolutionsFound_;}
  /// Increment how many solutions the heuristic thought it got
  inline void incrementNumberSolutionsFound()
  { numberSolutionsFound_++;}

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
  /// Set random number generator seed
  void setSeed(int value);
  /// Set input solution
  void setInputSolution(const double * solution, double objValue);

  /** Check whether the heuristic should run */
  bool shouldHeurRun();
  bool shouldHeurRun_randomChoice();
  void debugNodes();
  void printDistanceToNodes();
  /// how many times the heuristic has actually run
  inline int numRuns() const
  { return numRuns_;}

  /// How many times the heuristic could run
  inline int numCouldRun() const
  { return numCouldRun_;}

protected:

  /// Model
  CbcModel * model_;
  /// When flag - 0 off, 1 at root, 2 other than root, 3 always
  int when_;
  /// Number of nodes in any sub tree
  int numberNodes_;
  /// Feasibility pump options (-1 is off)
  int feasibilityPumpOptions_;
  /// Fraction of new(rows+columns)/old(rows+columns) before doing small branch and bound
  double fractionSmall_;
  /// Thread specific random number generator
  CoinThreadRandom randomNumberGenerator_;
  /// Name for printing
  std::string heuristicName_;

  /// How often to do (code can change)
  int howOften_;
  /// How much to increase how often
  double decayFactor_;
  /** Upto this depth we call the tree shallow and the heuristic can be called
      multiple times. That is, the test whether the current node is far from
      the others where the jeuristic was invoked will not be done, only the
      frequency will be tested. After that depth the heuristic will can be
      invoked only once per node, right before branching. That's when it'll be
      tested whether the heur should run at all. */
  int shallowDepth_;
  /** How often to invoke the heuristics in the shallow part of the tree */
  int howOftenShallow_;
  /** How many invocations happened within the same node when in a shallow
      part of the tree. */
  int numInvocationsInShallow_;
  /** How many invocations happened when in the deep part of the tree. For
      every node we count only one invocation. */
  int numInvocationsInDeep_;
  /** After how many deep invocations was the heuristic run last time */
  int lastRunDeep_;
  /// how many times the heuristic has actually run
  int numRuns_;
  /** How "far" should this node be from every other where the heuristic was
      run in order to allow the heuristic to run in this node, too. Currently
      this is tested, but we may switch to avgDistanceToRun_ in the future. */
  int minDistanceToRun_;

  /// The description of the nodes where this heuristic has been applied
  CbcHeuristicNodeList runNodes_;

  /// How many times the heuristic could run
  int numCouldRun_;

  /// How many solutions the heuristic thought it got
  int numberSolutionsFound_;

  // Input solution - so can be used as seed
  double * inputSolution_;


#if 0
  /// Lower bounds of last node where the heuristic found a solution
  double * lowerBoundLastNode_;
  /// Upper bounds of last node where the heuristic found a solution
  double * upperBoundLastNode_;
#endif
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
  
  using CbcHeuristic::solution ;
  /** returns 0 if no solution, 1 if valid solution
      with better objective value than one passed in
      Sets solution values if good, sets objective value (only if good)
      This is called after cuts have been added - so can not add cuts
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution);
  /** returns 0 if no solution, 1 if valid solution
      with better objective value than one passed in
      Sets solution values if good, sets objective value (only if good)
      This is called after cuts have been added - so can not add cuts
      Use solutionValue rather than solvers one
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution,
		       double solutionValue);
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

  // Down locks
  unsigned short * down_;

  // Up locks
  unsigned short * up_;

  // Equality locks
  unsigned short * equal_;

  // Seed for random stuff
  int seed_;
};

/** Partial solution class
    If user knows a partial solution this tries to get an integer solution
    it uses hotstart information
 */

class CbcHeuristicPartial : public CbcHeuristic {
public:

  // Default Constructor 
  CbcHeuristicPartial ();

  /** Constructor with model - assumed before cuts
      Fixes all variables with priority <= given
      and does given number of nodes
  */
  CbcHeuristicPartial (CbcModel & model, int fixPriority=10000, int numberNodes=200);
  
  // Copy constructor 
  CbcHeuristicPartial ( const CbcHeuristicPartial &);
   
  // Destructor 
  ~CbcHeuristicPartial ();
  
  /// Assignment operator 
  CbcHeuristicPartial & operator=(const CbcHeuristicPartial& rhs);

  /// Clone
  virtual CbcHeuristic * clone() const;
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
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution);
  /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
  virtual void validate();


  /// Set priority level
  void setFixPriority(int value)
  { fixPriority_ = value;}

protected:
  // Data

  // All variables with abs priority <= this will be fixed
  int fixPriority_;
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
  
  using CbcHeuristic::solution ;
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

/** Just One class - this chooses one at random
 */

class CbcHeuristicJustOne : public CbcHeuristic {
public:

  // Default Constructor 
  CbcHeuristicJustOne ();

  // Constructor with model - assumed before cuts
  CbcHeuristicJustOne (CbcModel & model);
  
  // Copy constructor 
  CbcHeuristicJustOne ( const CbcHeuristicJustOne &);
   
  // Destructor 
  ~CbcHeuristicJustOne ();

  /// Clone
  virtual CbcHeuristicJustOne * clone() const;
  
  /// Assignment operator 
  CbcHeuristicJustOne & operator=(const CbcHeuristicJustOne& rhs);

  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) ;

  /** returns 0 if no solution, 1 if valid solution
      with better objective value than one passed in
      Sets solution values if good, sets objective value (only if good)
      This is called after cuts have been added - so can not add cuts
      This does Fractional Diving
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution);
  /// Resets stuff if model changes
  virtual void resetModel(CbcModel * model);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel * model);
  /// Selects the next variable to branch on
  /** Returns true if all the fractional variables can be trivially
      rounded. Returns false, if there is at least one fractional variable
      that is not trivially roundable. In this case, the bestColumn
      returned will not be trivially roundable.
      This is dummy as never called
  */
  virtual bool selectVariableToBranch(OsiSolverInterface* solver,
				      const double* newSolution,
				      int& bestColumn,
				      int& bestRound) 
  { return true;}
  /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
  virtual void validate();
  /// Adds an heuristic with probability
  void addHeuristic(const CbcHeuristic * heuristic, double probability);
  /// Normalize probabilities
  void normalizeProbabilities();
protected:
  // Data

  // Probability of running a heuristic
  double * probabilities_;

  // Heuristics
  CbcHeuristic ** heuristic_;

  // Number of heuristics
  int numberHeuristics_;

};

#endif
