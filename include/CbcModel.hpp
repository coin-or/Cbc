// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcModel_H
#define CbcModel_H

#include <string>
#include <vector>

#include "CoinFinite.hpp"
#include "CoinMessageHandler.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CbcCompareBase.hpp"
#include "CbcMessage.hpp"

//class OsiSolverInterface;

class CbcCutGenerator;
class OsiRowCut;
class OsiRowCutDebugger;
class CglCutGenerator;
class CbcHeuristic;
class CbcObject;
class CbcTree;
class CbcStrategy;

//#############################################################################

/** Simple Branch and bound class

  The initialSolve() method solves the initial LP relaxation of the MIP
  problem. The branchAndBound() method can then be called to finish using
  a branch and cut algorithm.

  <h3>Search Tree Traversal</h3>

  Subproblems (aka nodes) requiring additional evaluation are stored using
  the CbcNode and CbcNodeInfo objects. Ancestry linkage is maintained in the
  CbcNodeInfo object. Evaluation of a subproblem within branchAndBound()
  proceeds as follows:
  <ul>
    <li> The node representing the most promising parent subproblem is popped
	 from the heap which holds the set of subproblems requiring further
	 evaluation.
    <li> Using branching instructions stored in the node, and information in
	 its ancestors, the model and solver are adjusted to create the
	 active subproblem.
    <li> If the parent subproblem will require further evaluation
	 (<i>i.e.</i>, there are branches remaining) its node is pushed back
	 on the heap. Otherwise, the node is deleted.  This may trigger
	 recursive deletion of ancestors.
    <li> The newly created subproblem is evaluated.
    <li> If the subproblem requires further evaluation, a node is created.
	 All information needed to recreate the subproblem (branching
	 information, row and column cuts) is placed in the node and the node
	 is added to the set of subproblems awaiting further evaluation.
  </ul>
  Note that there is never a node representing the active subproblem; the model
  and solver represent the active subproblem.

  <h3>Row (Constraint) Cut Handling</h3>

  For a typical subproblem, the sequence of events is as follows:
  <ul>
    <li> The subproblem is rebuilt for further evaluation: One result of a
	 call to addCuts() is a traversal of ancestors, leaving a list of all
	 cuts used in the ancestors in #addedCuts_. This list is then scanned
	 to construct a basis that includes only tight cuts. Entries for
	 loose cuts are set to NULL.
    <li> The subproblem is evaluated: One result of a call to solveWithCuts()
         is the return of a set of newly generated cuts for the subproblem.
	 #addedCuts_ is also kept up-to-date as old cuts become loose.
    <li> The subproblem is stored for further processing: A call to
	 CbcNodeInfo::addCuts() adds the newly generated cuts to the
	 CbcNodeInfo object associated with this node.
  </ul>
  See CbcCountRowCut for details of the bookkeeping associated with cut
  management.
*/

class CbcModel  {
  
public:

enum CbcIntParam {
  /** The maximum number of nodes before terminating */
  CbcMaxNumNode=0,
  /** The maximum number of solutions before terminating */
  CbcMaxNumSol,
  /** Fathoming discipline

    Controls objective function comparisons for purposes of fathoming by bound
    or determining monotonic variables.

    If 1, action is taken only when the current objective is strictly worse
    than the target. Implementation is handled by adding a small tolerance to
    the target.
  */
  CbcFathomDiscipline,
  /** Just a marker, so that a static sized array can store parameters. */
  CbcLastIntParam
};

enum CbcDblParam {
  /** The maximum amount the value of an integer variable can vary from
      integer and still be considered feasible. */
  CbcIntegerTolerance=0,
  /** The objective is assumed to worsen by this amount for each
      integer infeasibility. */
  CbcInfeasibilityWeight,
  /** The amount by which to tighten the objective function cutoff when
      a new solution is discovered. */
  CbcCutoffIncrement,
  /** Stop when the gap between the objective value of the best known solution
    and the best bound on the objective of any solution is less than this.
  
    This is an absolute value. Conversion from a percentage is left to the
    client.
  */
  CbcAllowableGap,
  /** Stop when the gap between the objective value of the best known solution
    and the best bound on the objective of any solution is less than this
    fraction of of the absolute value of best known solution.
  
    Code stops if either this test or CbcAllowableGap test succeeds
  */
  CbcAllowableFractionGap,
  /** \brief The maximum number of seconds before terminating.
	     A double should be adequate! */
  CbcMaximumSeconds,
  /** \brief The time at start of model.
	     So that other pieces of code can access */
  CbcStartSeconds,
  /** Just a marker, so that a static sized array can store parameters. */
  CbcLastDblParam
};

  //---------------------------------------------------------------------------

public:
  ///@name Solve methods 
  //@{
    /** \brief Solve the initial LP relaxation

      Invoke the solver's %initialSolve() method.
    */
    void initialSolve(); 

    /** \brief Invoke the branch \& cut algorithm

      The method assumes that initialSolve() has been called to solve the
      LP relaxation. It processes the root node, then proceeds to explore the
      branch & cut search tree. The search ends when the tree is exhausted or
      one of several execution limits is reached.
    */
     void branchAndBound();

    /** \brief create a clean model from partially fixed problem

      The method creates a new model with given bounds and with no tree.
    */
     CbcModel *  cleanModel(const double * lower, const double * upper);
    /** \brief Invoke the branch \& cut algorithm on partially fixed problem

      The method presolves the given model and does branch and cut. The search 
      ends when the tree is exhausted or maximum nodes is reached.

      If better solution found then it is saved.

      Returns 0 if search completed and solution, 1 if not completed and solution,
      2 if completed and no solution, 3 if not completed and no solution.

      Normally okay to do cleanModel immediately followed by subBranchandBound
      (== other form of subBranchAndBound)
      but may need to get at model for advanced features.

      Deletes model2
    */
     int subBranchAndBound(CbcModel * model2,
                           CbcModel * presolvedModel,
                           int maximumNodes);
    /** \brief Invoke the branch \& cut algorithm on partially fixed problem

      The method creates a new model with given bounds, presolves it
      then proceeds to explore the branch & cut search tree. The search 
      ends when the tree is exhausted or maximum nodes is reached.

      If better solution found then it is saved.

      Returns 0 if search completed and solution, 1 if not completed and solution,
      2 if completed and no solution, 3 if not completed and no solution.

      This is just subModel immediately followed by other version of
      subBranchandBound.

    */
     int subBranchAndBound(const double * lower, const double * upper,
			    int maximumNodes);

    /** \brief Process root node and return a strengthened model

      The method assumes that initialSolve() has been called to solve the
      LP relaxation. It processes the root node and then returns a pointer
      to the strengthened model (or NULL if infeasible)
    */
     OsiSolverInterface *  strengthenedModel();
    /** \brief Evaluate a subproblem using cutting planes and heuristics

      The method invokes a main loop which generates cuts, applies heuristics,
      and reoptimises using the solver's native %resolve() method.
      It returns true if the subproblem remains feasible at the end of the
      evaluation.
    */
    bool solveWithCuts(OsiCuts & cuts, int numberTries,CbcNode * node,
		       int & numberOldActiveCuts, int & numberNewCuts,
		       int & maximumWhich, int * & whichGenerator);

    /** \brief Reoptimise an LP relaxation
    
      Invoke the solver's %resolve() method.
    */
    bool resolve();
  /// Make given rows (L or G) into global cuts and remove from lp
  void makeGlobalCuts(int numberRows,const int * which); 
  //@}

  /** \name Presolve methods */
  //@{

  /** Identify cliques and construct corresponding objects.

      Find cliques with size in the range
      [\p atLeastThisMany, \p lessThanThis] and construct corresponding
      CbcClique objects.
      If \p makeEquality is true then a new model may be returned if
      modifications had to be made, otherwise \c this is returned.
      If the problem is infeasible #numberObjects_ is set to -1.
      A client must use deleteObjects() before a second call to findCliques().
      If priorities exist, clique priority is set to the default.
  */
  CbcModel * findCliques(bool makeEquality, int atLeastThisMany,
			 int lessThanThis, int defaultValue=1000);

  /** Do integer presolve, creating a new (presolved) model.

    Returns the new model, or NULL if feasibility is lost.
    If weak is true then just does a normal presolve
  
    \todo It remains to work out the cleanest way of getting a solution to
          the original problem at the end. So this is very preliminary.
   */
  CbcModel * integerPresolve(bool weak=false);

  /** Do integer presolve, modifying the current model.

      Returns true if the model remains feasible after presolve.
  */
  bool integerPresolveThisModel(OsiSolverInterface * originalSolver,bool weak=false);


  /// Put back information into the original model after integer presolve.
  void originalModel(CbcModel * presolvedModel,bool weak);

  /** \brief For variables involved in VUB constraints, see if we can tighten
	     bounds by solving lp's

      Returns false if feasibility is lost.
      If CglProbing is available, it will be tried as well to see if it can
      tighten bounds.
      This routine is just a front end for tightenVubs(int,const int*,double).

      If <tt>type = -1</tt> all variables are processed (could be very slow).
      If <tt>type = 0</tt> only variables involved in VUBs are processed.
      If <tt>type = n > 0</tt>, only the n most expensive VUB variables
      are processed, where it is assumed that x is at its maximum so delta
      would have to go to 1 (if x not at bound).

      If \p allowMultipleBinary is true, then a VUB constraint is a row with
      one continuous variable and any number of binary variables.

      If <tt>useCutoff < 1.0e30</tt>, the original objective is installed as a
      constraint with \p useCutoff as a bound.
  */
  bool tightenVubs(int type,bool allowMultipleBinary=false,
		   double useCutoff=1.0e50);
  
  /** \brief For variables involved in VUB constraints, see if we can tighten
	     bounds by solving lp's

    This version is just handed a list of variables to be processed.
  */
  bool tightenVubs(int numberVubs, const int * which,
		   double useCutoff=1.0e50);
  /**
    Analyze problem to find a minimum change in the objective function.
  */
  void analyzeObjective();


  //@}

  /** \name Object manipulation routines
  
    See CbcObject for an explanation of `object' in the context of CbcModel.
  */
  //@{

  /// Get the number of objects
  inline int numberObjects() const { return numberObjects_;};
  /// Set the number of objects
  inline void setNumberObjects(int number) 
  {  numberObjects_=number;};

  /// Get the array of objects
  inline CbcObject ** objects() const { return object_;};

  /// Get the specified object
  const inline CbcObject * object(int which) const { return object_[which];};

  /// Delete all object information
  void deleteObjects();

  /** Add in object information.
  
    Objects are cloned; the owner can delete the originals.
  */
  void addObjects(int numberObjects, CbcObject ** objects);

  /// Ensure attached objects point to this model.
  void synchronizeModel() ;

  /** \brief Identify integer variables and create corresponding objects.
  
    Record integer variables and create an CbcSimpleInteger object for each
    one.
    If \p startAgain is true, a new scan is forced, overwriting any existing
    integer variable information.
  */

  void findIntegers(bool startAgain);
  
  //@}

  //---------------------------------------------------------------------------

  /**@name Parameter set/get methods

     The set methods return true if the parameter was set to the given value,
     false if the value of the parameter is out of range.

     The get methods return the value of the parameter.

  */
  //@{
  /// Set an integer parameter
  inline bool setIntParam(CbcIntParam key, int value) {
    intParam_[key] = value;
    return true;
  }
  /// Set a double parameter
  inline bool setDblParam(CbcDblParam key, double value) {
    dblParam_[key] = value;
    return true;
  }
  /// Get an integer parameter
  inline int getIntParam(CbcIntParam key) const {
    return intParam_[key];
  }
  /// Get a double parameter
  inline double getDblParam(CbcDblParam key) const {
    return dblParam_[key];
  }
  /*! \brief Set cutoff bound on the objective function.

    When using strict comparison, the bound is adjusted by a tolerance to
    avoid accidentally cutting off the optimal solution.
  */
  void setCutoff(double value) ;

  /// Get the cutoff bound on the objective function - always as minimize
  inline double getCutoff() const
  { double value ;
    solver_->getDblParam(OsiDualObjectiveLimit,value) ;
    return value * solver_->getObjSense() ; } ;

  /// Set the \link CbcModel::CbcMaxNumNode maximum node limit \endlink
  inline bool setMaximumNodes( int value)
  { return setIntParam(CbcMaxNumNode,value); }

  /// Get the \link CbcModel::CbcMaxNumNode maximum node limit \endlink
  inline int getMaximumNodes() const
  { return getIntParam(CbcMaxNumNode); }

  /** Set the
      \link CbcModel::CbcMaxNumSol maximum number of solutions \endlink
      desired.
  */
  inline bool setMaximumSolutions( int value) {
    return setIntParam(CbcMaxNumSol,value);
  }
  /** Get the
      \link CbcModel::CbcMaxNumSol maximum number of solutions \endlink
      desired.
  */
  inline int getMaximumSolutions() const {
    return getIntParam(CbcMaxNumSol);
  }

  /** Set the
      \link CbcModel::CbcMaximumSeconds maximum number of seconds \endlink
      desired.
  */
  inline bool setMaximumSeconds( double value) {
    return setDblParam(CbcMaximumSeconds,value);
  }
  /** Get the
      \link CbcModel::CbcMaximumSeconds maximum number of seconds \endlink
      desired.
  */
  inline double getMaximumSeconds() const {
    return getDblParam(CbcMaximumSeconds);
  }

  /** Set the
    \link CbcModel::CbcIntegerTolerance integrality tolerance \endlink
  */
  inline bool setIntegerTolerance( double value) {
    return setDblParam(CbcIntegerTolerance,value);
  }
  /** Get the
    \link CbcModel::CbcIntegerTolerance integrality tolerance \endlink
  */
  inline double getIntegerTolerance() const {
    return getDblParam(CbcIntegerTolerance);
  }

  /** Set the
      \link CbcModel::CbcInfeasibilityWeight
	    weight per integer infeasibility \endlink
  */
  inline bool setInfeasibilityWeight( double value) {
    return setDblParam(CbcInfeasibilityWeight,value);
  }
  /** Get the
      \link CbcModel::CbcInfeasibilityWeight
	    weight per integer infeasibility \endlink
  */
  inline double getInfeasibilityWeight() const {
    return getDblParam(CbcInfeasibilityWeight);
  }

  /** Set the \link CbcModel::CbcAllowableGap allowable gap \endlink
      between the best known solution and the best possible solution.
  */
  inline bool setAllowableGap( double value) {
    return setDblParam(CbcAllowableGap,value);
  }
  /** Get the \link CbcModel::CbcAllowableGap allowable gap \endlink
      between the best known solution and the best possible solution.
  */
  inline double getAllowableGap() const {
    return getDblParam(CbcAllowableGap);
  }

  /** Set the \link CbcModel::CbcAllowableFractionGap fraction allowable gap \endlink
      between the best known solution and the best possible solution.
  */
  inline bool setAllowableFractionGap( double value) {
    return setDblParam(CbcAllowableFractionGap,value);
  }
  /** Get the \link CbcModel::CbcAllowableFractionGap fraction allowable gap \endlink
      between the best known solution and the best possible solution.
  */
  inline double getAllowableFractionGap() const {
    return getDblParam(CbcAllowableFractionGap);
  }
  /** Set the \link CbcModel::CbcAllowableFractionGap percentage allowable gap \endlink
      between the best known solution and the best possible solution.
  */
  inline bool setAllowablePercentageGap( double value) {
    return setDblParam(CbcAllowableFractionGap,value*0.01);
  }
  /** Get the \link CbcModel::CbcAllowableFractionGap percentage allowable gap \endlink
      between the best known solution and the best possible solution.
  */
  inline double getAllowablePercentageGap() const {
    return 100.0*getDblParam(CbcAllowableFractionGap);
  }

  /// Set the hotstart strategy
  void setHotstartStrategy(int value) 
  { hotstartStrategy_=value;};
  /// Get the hotstart strategy 
  int getHotstartStrategy() const
  { return hotstartStrategy_;};
  
  /// Set the minimum drop to continue cuts
  inline void setMinimumDrop(double value)
  {minimumDrop_=value;};
  /// Get the minimum drop to continue cuts
  inline double getMinimumDrop() const
  { return minimumDrop_;};

  /** Set the maximum number of cut passes at root node (default 20)
      Minimum drop can also be used for fine tuning */
  inline void setMaximumCutPassesAtRoot(int value)
  {maximumCutPassesAtRoot_=value;};
  /** Get the maximum number of cut passes at root node */
  inline int getMaximumCutPassesAtRoot() const
  { return maximumCutPassesAtRoot_;};

  /** Set the maximum number of cut passes at other nodes (default 10)
      Minimum drop can also be used for fine tuning */
  inline void setMaximumCutPasses(int value)
  {maximumCutPasses_=value;};
  /** Get the maximum number of cut passes at other nodes (default 10) */
  inline int getMaximumCutPasses() const
  { return maximumCutPasses_;};
  /** Get current cut pass number in this round of cuts.
      (1 is first pass) */
  inline int getCurrentPassNumber() const
  { return currentPassNumber_;};

  /** Set the maximum number of candidates to be evaluated for strong
    branching.

    A value of 0 disables strong branching.
  */
  void setNumberStrong(int number);
  /** Get the maximum number of candidates to be evaluated for strong
    branching.
  */
  inline int numberStrong() const
  { return numberStrong_;};

  /// Set how often to scan global cuts 
  void setHowOftenGlobalScan(int number);
  /// Get how often to scan global cuts
  inline int howOftenGlobalScan() const
  { return howOftenGlobalScan_;};

  /** Set the print frequency.
  
    Controls the number of nodes evaluated between status prints.
    If <tt>number <=0</tt> the print frequency is set to 100 nodes for large
    problems, 1000 for small problems.
    Print frequency has very slight overhead if small.
  */
  void setPrintFrequency(int number)
  { printFrequency_=number;};
  /// Get the print frequency
  inline int printFrequency() const
  { return printFrequency_;};
  //@}

  //---------------------------------------------------------------------------
  ///@name Methods returning info on how the solution process terminated
  //@{
    /// Are there a numerical difficulties?
    bool isAbandoned() const;
    /// Is optimality proven?
    bool isProvenOptimal() const;
    /// Is  infeasiblity proven (or none better than cutoff)?
    bool isProvenInfeasible() const;
    /// Node limit reached?
    bool isNodeLimitReached() const;
    /// Solution limit reached?
    bool isSolutionLimitReached() const;
    /// Get how many iterations it took to solve the problem.
    int getIterationCount() const
    { return solver_->getIterationCount();};
    /// Get how many Nodes it took to solve the problem.
    int getNodeCount() const
    { return numberNodes_;};
    /** Final status of problem
    
      0 finished, 1 stopped, 2 difficulties
    */
    inline int status() const
    { return status_;};
  
  //@}

  //---------------------------------------------------------------------------
  /**@name Problem information methods 
     
     These methods call the solver's query routines to return
     information about the problem referred to by the current object.
     Querying a problem that has no data associated with it result in
     zeros for the number of rows and columns, and NULL pointers from
     the methods that return vectors.
     
     Const pointers returned from any data-query method are valid as
     long as the data is unchanged and the solver is not called.
  */
  //@{
  /// Number of rows in continuous (root) problem.
  int numberRowsAtContinuous() const
  { return numberRowsAtContinuous_;};

  /// Get number of columns
  int getNumCols() const
  { return solver_->getNumCols();};
  
  /// Get number of rows
  int getNumRows() const
  { return solver_->getNumRows();};
  
  /// Get number of nonzero elements
  int getNumElements() const
  { return solver_->getNumElements();};

  /// Number of integers in problem
  inline int numberIntegers() const
  { return numberIntegers_;};
  // Integer variables
  inline const int * integerVariable() const 
  { return integerVariable_;};
  
  /// Get pointer to array[getNumCols()] of column lower bounds
  const double * getColLower() const
  { return solver_->getColLower();};
  
  /// Get pointer to array[getNumCols()] of column upper bounds
  const double * getColUpper() const
  { return solver_->getColUpper();};
  
  /** Get pointer to array[getNumRows()] of row constraint senses.
      <ul>
      <li>'L': <= constraint
      <li>'E': =  constraint
      <li>'G': >= constraint
      <li>'R': ranged constraint
      <li>'N': free constraint
      </ul>
  */
  const char * getRowSense() const
  { return solver_->getRowSense();};
  
  /** Get pointer to array[getNumRows()] of rows right-hand sides
      <ul>
      <li> if rowsense()[i] == 'L' then rhs()[i] == rowupper()[i]
      <li> if rowsense()[i] == 'G' then rhs()[i] == rowlower()[i]
      <li> if rowsense()[i] == 'R' then rhs()[i] == rowupper()[i]
      <li> if rowsense()[i] == 'N' then rhs()[i] == 0.0
      </ul>
  */
  const double * getRightHandSide() const
  { return solver_->getRightHandSide();};
  
  /** Get pointer to array[getNumRows()] of row ranges.
      <ul>
      <li> if rowsense()[i] == 'R' then
      rowrange()[i] == rowupper()[i] - rowlower()[i]
      <li> if rowsense()[i] != 'R' then
      rowrange()[i] is 0.0
      </ul>
  */
  const double * getRowRange() const
  { return solver_->getRowRange();};
  
  /// Get pointer to array[getNumRows()] of row lower bounds
  const double * getRowLower() const
  { return solver_->getRowLower();};
  
  /// Get pointer to array[getNumRows()] of row upper bounds
  const double * getRowUpper() const
  { return solver_->getRowUpper();};
  
  /// Get pointer to array[getNumCols()] of objective function coefficients
  const double * getObjCoefficients() const
  { return solver_->getObjCoefficients();};
  
  /// Get objective function sense (1 for min (default), -1 for max)
  double getObjSense() const
  { return solver_->getObjSense();};
  
  /// Return true if variable is continuous
  bool isContinuous(int colIndex) const
  { return solver_->isContinuous(colIndex);};
  
  /// Return true if variable is binary
  bool isBinary(int colIndex) const
  { return solver_->isBinary(colIndex);};
  
  /** Return true if column is integer.
      Note: This function returns true if the the column
      is binary or a general integer.
  */
  bool isInteger(int colIndex) const
  { return solver_->isInteger(colIndex);};
  
  /// Return true if variable is general integer
  bool isIntegerNonBinary(int colIndex) const
  { return solver_->isIntegerNonBinary(colIndex);};
  
  /// Return true if variable is binary and not fixed at either bound
  bool isFreeBinary(int colIndex) const
  { return solver_->isFreeBinary(colIndex) ;};
  
  /// Get pointer to row-wise copy of matrix
  const CoinPackedMatrix * getMatrixByRow() const
  { return solver_->getMatrixByRow();};
  
  /// Get pointer to column-wise copy of matrix
  const CoinPackedMatrix * getMatrixByCol() const
  { return solver_->getMatrixByCol();};
  
  /// Get solver's value for infinity
  double getInfinity() const
  { return solver_->getInfinity();};
  //@}
  
  
  /**@name Methods related to querying the solution */
  //@{
  /// Record a new incumbent solution and update objectiveValue
  void setBestSolution(CBC_Message how,
		       double & objectiveValue, const double *solution,
		       bool fixVariables=false);
  /// Just update objectiveValue
  void setBestObjectiveValue( double objectiveValue);

  /** Call this to really test if a valid solution can be feasible
      Solution is number columns in size.
      If fixVariables true then bounds of continuous solver updated.
      Returns objective value (worse than cutoff if not feasible)
 */
  double checkSolution(double cutoff, const double * solution,
		       bool fixVariables);
  /** Test the current solution for feasiblility.

    Scan all objects for indications of infeasibility. This is broken down
    into simple integer infeasibility (\p numberIntegerInfeasibilities)
    and all other reports of infeasibility (\p numberObjectInfeasibilities).
  */
  bool feasibleSolution(int & numberIntegerInfeasibilities,
			int & numberObjectInfeasibilities) const;

  /** Solution to the most recent lp relaxation.

    The solver's solution to the most recent lp relaxation.
  */
    
  inline double * currentSolution() const
  { return currentSolution_;};
  /// Make sure region there
  void reserveCurrentSolution();

  /// Get pointer to array[getNumCols()] of primal solution vector
  inline const double * getColSolution() const
  { return solver_->getColSolution();};
  
  /// Get pointer to array[getNumRows()] of dual prices
  inline const double * getRowPrice() const
  { return solver_->getRowPrice();};
  
  /// Get a pointer to array[getNumCols()] of reduced costs
  inline const double * getReducedCost() const
  { return solver_->getReducedCost();};
  
  /// Get pointer to array[getNumRows()] of row activity levels.
  inline const double * getRowActivity() const
  { return solver_->getRowActivity();};
  
  /// Get current objective function value
  inline double getCurrentObjValue() const
  { return solver_->getObjValue();};
  
  /// Get best objective function value as minimization
  inline double getMinimizationObjValue() const
  { return bestObjective_;};
  /// Set best objective function value as minimization
  inline void setMinimizationObjValue(double value) 
  { bestObjective_=value;};
  
  /// Get best objective function value
  inline double getObjValue() const
  { return bestObjective_ * solver_->getObjSense() ; } ;
  /** Get best possible objective function value.
      This is better of best possible left on tree
      and best solution found.
      If called from within branch and cut may be optimistic.
  */
  double getBestPossibleObjValue() const;
  /// Set best objective function value
  inline void setObjValue(double value) 
  { bestObjective_=value * solver_->getObjSense() ;};
  
  /** The best solution to the integer programming problem.

    The best solution to the integer programming problem found during
    the search. If no solution is found, the method returns null.
  */

  double * bestSolution() const
  { return bestSolution_;};
  
  /// Get number of solutions
  int getSolutionCount() const
  { return numberSolutions_;};
  
  /// Set number of solutions (so heuristics will be different)
  void setSolutionCount(int value) 
  { numberSolutions_=value;};
  /** Current phase (so heuristics etc etc can find out).
      0 - initial solve
      1 - solve with cuts at root
      2 - solve with cuts
      3 - other e.g. strong branching
      4 - trying to validate a solution
      5 - at end of search
  */
  inline int phase() const
  { return phase_;};
  
  /// Get number of heuristic solutions
  int getNumberHeuristicSolutions() const { return numberHeuristicSolutions_;};

  /// Set objective function sense (1 for min (default), -1 for max,)
  void setObjSense(double s) { solver_->setObjSense(s);};

  /// Value of objective at continuous
  inline double getContinuousObjective() const
  { return originalContinuousObjective_;};
  inline void setContinuousObjective(double value)
  { originalContinuousObjective_=value;};
  /// Number of infeasibilities at continuous
  inline int getContinuousInfeasibilities() const
  { return continuousInfeasibilities_;};
  inline void setContinuousInfeasibilities(int value)
  { continuousInfeasibilities_=value;};
  /// Value of objective after root node cuts added
  inline double rootObjectiveAfterCuts() const
  { return continuousObjective_;};
  /** Number of times global cuts violated.  When global cut pool then this
      should be kept for each cut and type of cut */
  inline int numberGlobalViolations() const
  { return numberGlobalViolations_;};
  inline void clearNumberGlobalViolations()
  { numberGlobalViolations_=0;};
  /// Whether to force a resolve after takeOffCuts
  inline bool resolveAfterTakeOffCuts() const
  { return resolveAfterTakeOffCuts_;};
  inline void setResolveAfterTakeOffCuts(bool yesNo)
  { resolveAfterTakeOffCuts_=yesNo;};
  //@}

  /** \name Node selection */
  //@{
  // Comparison functions (which may be overridden by inheritance)
  inline CbcCompareBase * nodeComparison() const
  { return nodeCompare_;};
  void setNodeComparison(CbcCompareBase * compare);
  void setNodeComparison(CbcCompareBase & compare);
  //@}

  /** \name Tree methods and subtree methods */
  //@{
  /// Tree method e.g. heap (which may be overridden by inheritance)
  inline CbcTree * tree() const
  { return tree_;};
  /// For modifying tree handling (original is cloned)
  void passInTreeHandler(CbcTree & tree);
  /** For passing in an CbcModel to do a sub Tree (with derived tree handlers).
      Passed in model must exist for duration of branch and bound
  */
  void passInSubTreeModel(CbcModel & model);
  /** For retrieving a copy of subtree model with given OsiSolver.
      If no subtree model will use self (up to user to reset cutoff etc).
      If solver NULL uses current 
  */
  CbcModel * subTreeModel(OsiSolverInterface * solver=NULL) const;
  /// Returns number of times any subtree stopped on nodes, time etc
  inline int numberStoppedSubTrees() const
  { return numberStoppedSubTrees_;}
  /// Says a sub tree was stopped
  inline void incrementSubTreeStopped()
  { numberStoppedSubTrees_++;};
  /** Whether to automatically do presolve before branch and bound (subTrees).
      0 - no
      1 - ordinary presolve
      2 - integer presolve (dodgy)
  */
  inline int typePresolve() const
  { return presolve_;};
  inline void setTypePresolve(int value)
  { presolve_=value;};
  //@}

  /** \name Branching Decisions
  
    See the CbcBranchDecision class for additional information.
  */
  //@{

  /// Get the current branching decision method.
  inline CbcBranchDecision * branchingMethod() const
  { return branchingMethod_;};
  /// Set the branching decision method.
  inline void setBranchingMethod(CbcBranchDecision * method)
  { branchingMethod_ = method;};
  /** Set the branching method
  
    \overload
  */
  inline void setBranchingMethod(CbcBranchDecision & method)
  { branchingMethod_ = &method;};
  //@}

  /** \name Row (constraint) and Column (variable) cut generation */
  //@{

  /** Perform reduced cost fixing

    Fixes integer variables at their current value based on reduced cost
    penalties.
  */
  void reducedCostFix() ;

  /** Return an empty basis object of the specified size

    A useful utility when constructing a basis for a subproblem from scratch.
    The object returned will be of the requested capacity and appropriate for
    the solver attached to the model.
  */
  CoinWarmStartBasis *getEmptyBasis(int ns = 0, int na = 0) const ;

  /** Remove inactive cuts from the model

    An OsiSolverInterface is expected to maintain a valid basis, but not a
    valid solution, when loose cuts are deleted. Restoring a valid solution
    requires calling the solver to reoptimise. If it's certain the solution
    will not be required, set allowResolve to false to suppress
    reoptimisation.
  */
  void takeOffCuts(OsiCuts &cuts, int *whichGenerator,
		     int &numberOldActiveCuts, int &numberNewCuts,
		     bool allowResolve) ;

  /** Determine and install the active cuts that need to be added for
    the current subproblem

    The whole truth is a bit more complicated. The first action is a call to
    addCuts1(). addCuts() then sorts through the list, installs the tight
    cuts in the model, and does bookkeeping (adjusts reference counts).
    The basis returned from addCuts1() is adjusted accordingly.
    
    If it turns out that the node should really be fathomed by bound,
    addCuts() simply treats all the cuts as loose as it does the bookkeeping.
  */
  int addCuts(CbcNode * node, CoinWarmStartBasis *&lastws);

  /** Traverse the tree from node to root and prep the model

    addCuts1() begins the job of prepping the model to match the current
    subproblem. The model is stripped of all cuts, and the search tree is
    traversed from node to root to determine the changes required. Appropriate
    bounds changes are installed, a list of cuts is collected but not
    installed, and an appropriate basis (minus the cuts, but big enough to
    accommodate them) is constructed.

    \todo addCuts1() is called in contexts where it's known in advance that
	  all that's desired is to determine a list of cuts and do the
	  bookkeeping (adjust the reference counts). The work of installing
	  bounds and building a basis goes to waste.
  */
  void addCuts1(CbcNode * node, CoinWarmStartBasis *&lastws);

  /// Return the list of cuts initially collected for this subproblem
  CbcCountRowCut ** addedCuts() const
  { return addedCuts_;};
  /// Number of entries in the list returned by #addedCuts()
  int currentNumberCuts() const
  { return currentNumberCuts_;};
  /// Global cuts
  inline OsiCuts * globalCuts() 
  { return &globalCuts_;};
  /// Copy and set a pointer to a row cut which will be added instead of normal branching.
  void setNextRowCut(const OsiRowCut & cut);
  /// Get a pointer to current node (be careful)
  inline CbcNode * currentNode() const
  { return currentNode_;};

  /// Get the number of cut generators
  inline int numberCutGenerators() const
  { return numberCutGenerators_;};
  /// Get the list of cut generators
  inline CbcCutGenerator ** cutGenerators() const
  { return generator_;};
  ///Get the specified cut generator
  inline CbcCutGenerator * cutGenerator(int i) const
  { return generator_[i];};
  ///Get the specified cut generator before any changes
  inline CbcCutGenerator * virginCutGenerator(int i) const
  { return virginGenerator_[i];};
  /** Add one generator - up to user to delete generators.
      howoften affects how generator is used. 0 or 1 means always,
      >1 means every that number of nodes.  Negative values have same
      meaning as positive but they may be switched off (-> -100) by code if
      not many cuts generated at continuous.  -99 is just done at root.
      Name is just for printout.
      If depth >0 overrides how often generator is called (if howOften==-1 or >0).
  */
  void addCutGenerator(CglCutGenerator * generator,
		       int howOften=1, const char * name=NULL,
		       bool normal=true, bool atSolution=false, 
		       bool infeasible=false,int howOftenInSub=-100,
		       int whatDepth=-1, int whatDepthInSub=-1);
//@}
  /** \name Strategy and sub models
  
    See the CbcStrategy class for additional information.
  */
  //@{

  /// Get the current strategy
  inline CbcStrategy * strategy() const
  { return strategy_;};
  /// Set the strategy. Clones
  void setStrategy(CbcStrategy & strategy);
  /// Get the current parent model
  inline CbcModel * parentModel() const
  { return parentModel_;};
  /// Set the parent model
  inline void setParentModel(CbcModel & parentModel)
  { parentModel_ = &parentModel;};
  //@}


  /** \name Heuristics and priorities */
  //@{
  /// Add one heuristic - up to user to delete
  void addHeuristic(CbcHeuristic * generator);
  ///Get the specified heuristic
  inline CbcHeuristic * heuristic(int i) const
  { return heuristic_[i];};

  /** Pass in branching priorities.
  
      If ifClique then priorities are on cliques otherwise priorities are
      on integer variables.  
      Other type (if exists set to default)
      1 is highest priority. (well actually -INT_MAX is but that's ugly)
      If hotstart > 0 then branches are created to force
      the variable to the value given by best solution.  This enables a
      sort of hot start.  The node choice should be greatest depth
      and hotstart should normally be switched off after a solution.

      If ifNotSimpleIntegers true then appended to normal integers

      \internal Added for Kurt Spielberg.
  */
  void passInPriorities(const int * priorities, bool ifNotSimpleIntegers,
			int defaultValue=1000);

  /// Priorities
  inline const int * priority() const { return priority_;};

  /// Returns priority level for an object (or 1000 if no priorities exist)
  inline int priority(int sequence) const
  { 
    if (priority_)
      return priority_[sequence];
    else
      return 1000;
  };
  //@}
    
  /**@name Setting/Accessing application data */
  //@{
    /** Set application data.

	This is a pointer that the application can store into and
	retrieve from the solver interface.
	This field is available for the application to optionally
	define and use.
    */
    void setApplicationData (void * appData);

    /// Get application data
    void * getApplicationData() const;
  //@}
  
  //---------------------------------------------------------------------------

  /**@name Message handling */
  //@{
  /// Pass in Message handler (not deleted at end)
  void passInMessageHandler(CoinMessageHandler * handler);
  /// Set language
  void newLanguage(CoinMessages::Language language);
  inline void setLanguage(CoinMessages::Language language)
  {newLanguage(language);};
  /// Return handler
  inline CoinMessageHandler * messageHandler() const
  {return handler_;};
  /// Return messages
  inline CoinMessages messages() 
  {return messages_;};
  /// Return pointer to messages
  inline CoinMessages * messagesPointer() 
  {return &messages_;};
  /// Set log level
  inline void setLogLevel(int value)
  { handler_->setLogLevel(value);};
  /// Get log level
  inline int logLevel() const
  { return handler_->logLevel();};
  //@}
  //---------------------------------------------------------------------------


  ///@name Constructors and destructors etc
  //@{
    /// Default Constructor
    CbcModel(); 
    
    /// Constructor from solver
    CbcModel(const OsiSolverInterface &);
  
    /** Assign a solver to the model (model assumes ownership)

      On return, \p solver will be NULL.

      \note Parameter settings in the outgoing solver are not inherited by
	    the incoming solver.
    */
    void assignSolver(OsiSolverInterface *&solver);
  
    /** Copy constructor .
      If noTree is true then tree and cuts are not copied
    */  
    CbcModel(const CbcModel & rhs, bool noTree=false);
  
    /// Assignment operator 
    CbcModel & operator=(const CbcModel& rhs);
  
    /// Destructor 
     ~CbcModel ();

    /// Returns solver - has current state
    OsiSolverInterface * solver() const
    { return solver_;};

    /// Returns solver with continuous state
    OsiSolverInterface * continuousSolver() const
    { return continuousSolver_;};
  /// Clears out as much as possible (except solver)
  void gutsOfDestructor();
  //@}

//---------------------------------------------------------------------------

private:
  ///@name Private member data 
  //@{

  /// The solver associated with this model.
  OsiSolverInterface * solver_;

  /** Ownership of the solver object

    The convention is that CbcModel owns the null solver. Currently there
    is no public method to give CbcModel a solver without giving ownership,
    but the hook is here.
  */
  bool ourSolver_ ;

  /// A copy of the solver, taken at the continuous (root) node.
  OsiSolverInterface * continuousSolver_;

   /// Message handler
  CoinMessageHandler * handler_;

  /** Flag to say if handler_ is the default handler.
  
    The default handler is deleted when the model is deleted. Other
    handlers (supplied by the client) will not be deleted.
  */
  bool defaultHandler_;

  /// Cbc messages
  CoinMessages messages_;

  /// Array for integer parameters
  int intParam_[CbcLastIntParam];

  /// Array for double parameters
  double dblParam_[CbcLastDblParam];

  /** Pointer to an empty warm start object

    It turns out to be useful to have this available as a base from
    which to build custom warm start objects. This is typed as CoinWarmStart
    rather than CoinWarmStartBasis to allow for the possibility that a
    client might want to apply a solver that doesn't use a basis-based warm
    start. See getEmptyBasis for an example of how this field can be used.
  */
  mutable CoinWarmStart *emptyWarmStart_ ;

  /** Pointer to a warm start basis.  */
  CoinWarmStartBasis *basis_;

  /// Best objective
  double bestObjective_;
  /// Best possible objective
  double bestPossibleObjective_;

  /// Array holding the incumbent (best) solution.
  double * bestSolution_;

  /** Array holding the current solution.

    This array is used more as a temporary.
  */
  double * currentSolution_;

  /// Global cuts
  OsiCuts globalCuts_;

  /// Minimum degradation in objective value to continue cut generation
  double minimumDrop_;
  /// Number of solutions
  int numberSolutions_;
  /// Hotstart strategy 0 =off, 1=branch if incorrect,2=branch even if correct, ....
  int hotstartStrategy_;
  /// Number of heuristic solutions
  int numberHeuristicSolutions_;
  /// Cumulative number of nodes
  int numberNodes_;
  /// Cumulative number of iterations
  int numberIterations_;
  /// Status of problem - 0 finished, 1 stopped, 2 difficulties
  int status_;
  /// Number of integers in problem
  int numberIntegers_;
  /// Number of rows at continuous
  int numberRowsAtContinuous_;
  /// Maximum number of cuts
  int maximumNumberCuts_;
  /** Current phase (so heuristics etc etc can find out).
      0 - initial solve
      1 - solve with cuts at root
      2 - solve with cuts
      3 - other e.g. strong branching
      4 - trying to validate a solution
      5 - at end of search
  */
  int phase_;

  /// Number of entries in #addedCuts_
  int currentNumberCuts_;

  /** Current limit on search tree depth

    The allocated size of #walkback_. Increased as needed.
  */
  int maximumDepth_;
  /** Array used to assemble the path between a node and the search tree root

    The array is resized when necessary. #maximumDepth_  is the current
    allocated size.
  */
  CbcNodeInfo ** walkback_;

  /** The list of cuts initially collected for this subproblem

    When the subproblem at this node is rebuilt, a set of cuts is collected
    for inclusion in the constraint system. If any of these cuts are
    subsequently removed because they have become loose, the corresponding
    entry is set to NULL.
  */
  CbcCountRowCut ** addedCuts_;

  /** A pointer to a row cut which will be added instead of normal branching.
      After use it should be set to NULL.
  */
  OsiRowCut * nextRowCut_;

  /// Current node so can be used elsewhere
  CbcNode * currentNode_;

  /// Indices of integer variables
  int * integerVariable_;
  /// 0 bit - check if cuts valid (if on list)
  int specialOptions_;
  /// User node comparison function
  CbcCompareBase * nodeCompare_;
  /// Tree
  CbcTree * tree_;
  /// A pointer to model to be used for subtrees
  CbcModel * subTreeModel_;
  /// Number of times any subtree stopped on nodes, time etc
  int numberStoppedSubTrees_;
  /// Variable selection function
  CbcBranchDecision * branchingMethod_;
  /// Strategy
  CbcStrategy * strategy_;
  /// Parent model
  CbcModel * parentModel_;
  /** Whether to automatically do presolve before branch and bound.
      0 - no
      1 - ordinary presolve
      2 - integer presolve (dodgy)
  */
  /// Pointer to user-defined data structure
  void * appData_;
  int presolve_;
  /** Maximum number of candidates to consider for strong branching.

    To disable strong branching, set this to 0.
  */
  int numberStrong_;

  /// Print frequency
  int printFrequency_;
  /// Number of cut generators
  int numberCutGenerators_;
  // Cut generators
  CbcCutGenerator ** generator_;
  // Cut generators before any changes
  CbcCutGenerator ** virginGenerator_;
  /// Number of heuristics
  int numberHeuristics_;
  // Heuristic solvers
  CbcHeuristic ** heuristic_;

  /// Total number of objects
  int numberObjects_;

  /** \brief Integer and Clique and ... information

    \note The code assumes that the first objects on the list will be
	  SimpleInteger objects for each integer variable, followed by
	  Clique objects. Portions of the code that understand Clique objects
	  will fail if they do not immediately follow the SimpleIntegers.
	  Large chunks of the code will fail if the first objects are not
	  SimpleInteger. As of 2003.08, SimpleIntegers and Cliques are the only
	  objects.
  */
  CbcObject ** object_;

  
  /// Original columns as created by integerPresolve
  int * originalColumns_;
  /// Priorities
  int * priority_;
  /// How often to scan global cuts
  int howOftenGlobalScan_;
  /** Number of times global cuts violated.  When global cut pool then this
      should be kept for each cut and type of cut */
  int numberGlobalViolations_;
  /** Value of objective at continuous
      (Well actually after initial round of cuts)
  */
  double continuousObjective_;
  /** Value of objective before root node cuts added
  */
  double originalContinuousObjective_;
  /// Number of infeasibilities at continuous
  int continuousInfeasibilities_;
  /// Maximum number of cut passes at root
  int maximumCutPassesAtRoot_;
  /// Maximum number of cut passes
  int maximumCutPasses_;
  /// Current cut pass number
  int currentPassNumber_;
  /// Whether to force a resolve after takeOffCuts
  bool resolveAfterTakeOffCuts_;
 //@}
};

#endif
