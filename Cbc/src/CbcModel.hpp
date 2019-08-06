/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcModel_H
#define CbcModel_H
#include <string>
#include <vector>
#include "CoinMessageHandler.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiBranchingObject.hpp"
#include "OsiCuts.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CbcCompareBase.hpp"
#include "CbcCountRowCut.hpp"
#include "CbcMessage.hpp"
#include "CbcEventHandler.hpp"
#include "ClpDualRowPivot.hpp"

class CbcCutGenerator;
class CbcBaseModel;
class OsiRowCut;
class OsiBabSolver;
class OsiRowCutDebugger;
class CglCutGenerator;
class CglStored;
class CbcCutModifier;
class CglTreeProbingInfo;
class CbcHeuristic;
class OsiObject;
class CbcThread;
class CbcTree;
class CbcStrategy;
class CbcSymmetry;
class CbcFeasibilityBase;
class CbcStatistics;
class CbcFullNodeInfo;
class CbcEventHandler;
class CglPreProcess;
class OsiClpSolverInterface;
class ClpNodeStuff;

// #define CBC_CHECK_BASIS 1

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

class CbcModel {

public:
  enum CbcIntParam {
    /** The maximum number of nodes before terminating */
    CbcMaxNumNode = 0,
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
    /** Adjusts printout
            1 does different node message with number unsatisfied on last branch
        */
    CbcPrinting,
    /** Number of branches (may be more than number of nodes as may
            include strong branching) */
    CbcNumberBranches,
    /** Just a marker, so that a static sized array can store parameters. */
    CbcLastIntParam
  };

  enum CbcDblParam {
    /** The maximum amount the value of an integer variable can vary from
            integer and still be considered feasible. */
    CbcIntegerTolerance = 0,
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
    /// Cutoff - stored for speed
    CbcCurrentCutoff,
    /// Optimization direction - stored for speed
    CbcOptimizationDirection,
    /// Current objective value
    CbcCurrentObjectiveValue,
    /// Current minimization objective value
    CbcCurrentMinimizationObjectiveValue,
    /** \brief The time at start of model.
               So that other pieces of code can access */
    CbcStartSeconds,
    /** Stop doing heuristics when the gap between the objective value of the
            best known solution and the best bound on the objective of any solution
            is less than this.

          This is an absolute value. Conversion from a percentage is left to the
          client.
        */
    CbcHeuristicGap,
    /** Stop doing heuristics when the gap between the objective value of the
            best known solution and the best bound on the objective of any solution
            is less than this fraction of of the absolute value of best known
            solution.

          Code stops if either this test or CbcAllowableGap test succeeds
        */
    CbcHeuristicFractionGap,
    /// Smallest non-zero change on a branch
    CbcSmallestChange,
    /// Sum of non-zero changes on a branch
    CbcSumChange,
    /// Largest non-zero change on a branch
    CbcLargestChange,
    /// Small non-zero change on a branch to be used as guess
    CbcSmallChange,
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
      If doStatistics is 1 summary statistics are printed
      if 2 then also the path to best solution (if found by branching)
      if 3 then also one line per node
    */
  void branchAndBound(int doStatistics = 0);

private:
  /** \brief Evaluate a subproblem using cutting planes and heuristics

      The method invokes a main loop which generates cuts, applies heuristics,
      and reoptimises using the solver's native %resolve() method.
      It returns true if the subproblem remains feasible at the end of the
      evaluation.
    */
  bool solveWithCuts(OsiCuts &cuts, int numberTries, CbcNode *node);
  /** Generate one round of cuts - serial mode
      returns -
      0 - normal
      1 - must keep going
      2 - set numberTries to zero
      -1 - infeasible
    */
  int serialCuts(OsiCuts &cuts, CbcNode *node, OsiCuts &slackCuts, int lastNumberCuts);
  /** Generate one round of cuts - parallel mode
        returns -
        0 - normal
        1 - must keep going
        2 - set numberTries to zero
        -1 - infeasible
    */
  int parallelCuts(CbcBaseModel *master, OsiCuts &cuts, CbcNode *node, OsiCuts &slackCuts, int lastNumberCuts);
  /** Input one node output N nodes to put on tree and optional solution update
        This should be able to operate in parallel so is given a solver and is const(ish)
        However we will need to keep an array of solver_ and bases and more
        status is 0 for normal, 1 if solution
        Calling code should always push nodes back on tree
    */
  CbcNode **solveOneNode(int whichSolver, CbcNode *node,
    int &numberNodesOutput, int &status);
  /// Update size of whichGenerator
  void resizeWhichGenerator(int numberNow, int numberAfter);

public:
#ifdef CBC_KEEP_DEPRECATED
  // See if anyone is using these any more!!
  /** \brief create a clean model from partially fixed problem

      The method creates a new model with given bounds and with no tree.
    */
  CbcModel *cleanModel(const double *lower, const double *upper);
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
  int subBranchAndBound(CbcModel *model2,
    CbcModel *presolvedModel,
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
  int subBranchAndBound(const double *lower, const double *upper,
    int maximumNodes);

  /** \brief Process root node and return a strengthened model

      The method assumes that initialSolve() has been called to solve the
      LP relaxation. It processes the root node and then returns a pointer
      to the strengthened model (or NULL if infeasible)
    */
  OsiSolverInterface *strengthenedModel();
  /** preProcess problem - replacing solver
        If makeEquality true then <= cliques converted to ==.
        Presolve will be done numberPasses times.

        Returns NULL if infeasible

        If makeEquality is 1 add slacks to get cliques,
        if 2 add slacks to get sos (but only if looks plausible) and keep sos info
    */
  CglPreProcess *preProcess(int makeEquality = 0, int numberPasses = 5,
    int tuning = 5);
  /** Does postprocessing - original solver back.
        User has to delete process */
  void postProcess(CglPreProcess *process);
#endif
  /// Returns CglPreProcess used before branch and bound
  inline CglPreProcess *preProcess() const
  {
    return preProcess_;
  }
  /// Set CglPreProcess used before branch and bound
  inline void setPreProcess(CglPreProcess *preProcess)
  {
    preProcess_ = preProcess;
  }
  /// Adds an update information object
  void addUpdateInformation(const CbcObjectUpdateData &data);
  /** Do one node - broken out for clarity?
        also for parallel (when baseModel!=this)
        Returns 1 if solution found
        node NULL on return if no branches left
        newNode NULL if no new node created
    */
  int doOneNode(CbcModel *baseModel, CbcNode *&node, CbcNode *&newNode);

public:
  /** \brief Reoptimise an LP relaxation

      Invoke the solver's %resolve() method.
      whereFrom -
      0 - initial continuous
      1 - resolve on branch (before new cuts)
      2 - after new cuts
      3  - obsolete code or something modified problem in unexpected way
      10 - after strong branching has fixed variables at root
      11 - after strong branching has fixed variables in tree

      returns 1 feasible, 0 infeasible, -1 feasible but skip cuts
    */
  int resolve(CbcNodeInfo *parent, int whereFrom,
    double *saveSolution = NULL,
    double *saveLower = NULL,
    double *saveUpper = NULL);
  /// Make given rows (L or G) into global cuts and remove from lp
  void makeGlobalCuts(int numberRows, const int *which);
  /// Make given cut into a global cut
  int makeGlobalCut(const OsiRowCut *cut);
  /// Make given cut into a global cut
  int makeGlobalCut(const OsiRowCut &cut);
  /// Make given column cut into a global cut
  void makeGlobalCut(const OsiColCut *cut);
  /// Make given column cut into a global cut
  void makeGlobalCut(const OsiColCut &cut);
  /// Make partial cut into a global cut and save
  void makePartialCut(const OsiRowCut *cut, const OsiSolverInterface *solver = NULL);
  /// Make partial cuts into global cuts
  void makeGlobalCuts();
  /// Which cut generator generated this cut
  inline const int *whichGenerator() const
  {
    return whichGenerator_;
  }
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
  CbcModel *findCliques(bool makeEquality, int atLeastThisMany,
    int lessThanThis, int defaultValue = 1000);

  /** Do integer presolve, creating a new (presolved) model.

      Returns the new model, or NULL if feasibility is lost.
      If weak is true then just does a normal presolve

      \todo It remains to work out the cleanest way of getting a solution to
            the original problem at the end. So this is very preliminary.
     */
  CbcModel *integerPresolve(bool weak = false);

  /** Do integer presolve, modifying the current model.

        Returns true if the model remains feasible after presolve.
    */
  bool integerPresolveThisModel(OsiSolverInterface *originalSolver, bool weak = false);

  /// Put back information into the original model after integer presolve.
  void originalModel(CbcModel *presolvedModel, bool weak);

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
  bool tightenVubs(int type, bool allowMultipleBinary = false,
    double useCutoff = 1.0e50);

  /** \brief For variables involved in VUB constraints, see if we can tighten
           bounds by solving lp's

      This version is just handed a list of variables to be processed.
    */
  bool tightenVubs(int numberVubs, const int *which,
    double useCutoff = 1.0e50);
  /**
      Analyze problem to find a minimum change in the objective function.
    */
  void analyzeObjective();
  /** Returns postProcessed solution in solver(called from event handler)
     Normally used for integer solution (not really tested otherwise)
    solutionType 1 is best integer so far, 0 is current solution 
    (may not be integer) */
  const OsiSolverInterface *postProcessedSolver(int solutionType = 1);

  /**
      Add additional integers.
    */
  void AddIntegers();
  /**
      Save copy of the model.
    */
  void saveModel(OsiSolverInterface *saveSolver, double *checkCutoffForRestart, bool *feasible);
  /**
      Flip direction of optimization on all models
    */
  void flipModel();
  /**
       Clean model i.e. make SOS/integer variables exactly at bound if needed.
       Only if moreSpecialOptions2_ 15 bit set (32768) as there is a small 
       overhead (more2 in standalone cbc).
       Fine tuning can be done in configure with -DCLEAN_INTEGER_VARIABLES
       and -DZERO_ODD_TOLERANCE=1.0e-nn
       If CLEAN_INTEGER_VARIABLES not defined then cleaning is only done for 
       SOS variables.
       If ZERO_ODD_TOLERANCE not defined then 1.0e-14 used.  You can define as 
       0.0 if you are paranoid. 
       Returns number of variables forced out
       cleanVariables array will be used if exists
     */
  int cleanBounds(OsiSolverInterface *solver, char *cleanVariables);
  /// Sets up cleanVariables array (i.e. ones to be careful about)
  char *setupCleanVariables();
  //@}

  /** \name Object manipulation routines

      See OsiObject for an explanation of `object' in the context of CbcModel.
    */
  //@{

  /// Get the number of objects
  inline int numberObjects() const
  {
    return numberObjects_;
  }
  /// Set the number of objects
  inline void setNumberObjects(int number)
  {
    numberObjects_ = number;
  }

  /// Get the array of objects
  inline OsiObject **objects() const
  {
    return object_;
  }

  /// Get the specified object
  const inline OsiObject *object(int which) const
  {
    return object_[which];
  }
  /// Get the specified object
  inline OsiObject *modifiableObject(int which) const
  {
    return object_[which];
  }

  void setOptionalInteger(int index);

  /// Delete all object information (and just back to integers if true)
  void deleteObjects(bool findIntegers = true);

  /** Add in object information.

      Objects are cloned; the owner can delete the originals.
    */
  void addObjects(int numberObjects, OsiObject **objects);

  /** Add in object information.

      Objects are cloned; the owner can delete the originals.
    */
  void addObjects(int numberObjects, CbcObject **objects);

  /// Ensure attached objects point to this model.
  void synchronizeModel();

  /** \brief Identify integer variables and create corresponding objects.

      Record integer variables and create an CbcSimpleInteger object for each
      one.
      If \p startAgain is true, a new scan is forced, overwriting any existing
      integer variable information.
      If type > 0 then 1==PseudoCost, 2 new ones low priority
    */

  void findIntegers(bool startAgain, int type = 0);
  /** Add SOS info to solver -
	Overwrites SOS information in solver with information
	in CbcModel.  Has no effect with some solvers. 
	Also updates integer info. */
  void addSOSEtcToSolver();

#ifdef SWITCH_VARIABLES
  /// Convert Dynamic to Switching
  int findSwitching();
  /// Fix associated variables
  int fixAssociated(OsiSolverInterface *solver, int cleanBasis);
  /// Debug associated variables
  int checkAssociated(const OsiSolverInterface *solver,
    const double *solution, int printLevel);
#endif
  //@}

  //---------------------------------------------------------------------------

  /**@name Parameter set/get methods

       The set methods return true if the parameter was set to the given value,
       false if the value of the parameter is out of range.

       The get methods return the value of the parameter.

    */
  //@{
  /// Set an integer parameter
  inline bool setIntParam(CbcIntParam key, int value)
  {
    intParam_[key] = value;
    return true;
  }
  /// Set a double parameter
  inline bool setDblParam(CbcDblParam key, double value)
  {
    dblParam_[key] = value;
    return true;
  }
  /// Get an integer parameter
  inline int getIntParam(CbcIntParam key) const
  {
    return intParam_[key];
  }
  /// Get a double parameter
  inline double getDblParam(CbcDblParam key) const
  {
    return dblParam_[key];
  }
  /*! \brief Set cutoff bound on the objective function.

      When using strict comparison, the bound is adjusted by a tolerance to
      avoid accidentally cutting off the optimal solution.
    */
  void setCutoff(double value);

  /// Get the cutoff bound on the objective function - always as minimize
  inline double getCutoff() const
  { //double value ;
    //solver_->getDblParam(OsiDualObjectiveLimit,value) ;
    //assert( dblParam_[CbcCurrentCutoff]== value * solver_->getObjSense());
    return dblParam_[CbcCurrentCutoff];
  }

  /// Set the \link CbcModel::CbcMaxNumNode maximum node limit \endlink
  inline bool setMaximumNodes(int value)
  {
    return setIntParam(CbcMaxNumNode, value);
  }

  /// Get the \link CbcModel::CbcMaxNumNode maximum node limit \endlink
  inline int getMaximumNodes() const
  {
    return getIntParam(CbcMaxNumNode);
  }

  /** Set the
        \link CbcModel::CbcMaxNumSol maximum number of solutions \endlink
        desired.
    */
  inline bool setMaximumSolutions(int value)
  {
    return setIntParam(CbcMaxNumSol, value);
  }
  /** Get the
        \link CbcModel::CbcMaxNumSol maximum number of solutions \endlink
        desired.
    */
  inline int getMaximumSolutions() const
  {
    return getIntParam(CbcMaxNumSol);
  }
  /// Set the printing mode
  inline bool setPrintingMode(int value)
  {
    return setIntParam(CbcPrinting, value);
  }

  /// Get the printing mode
  inline int getPrintingMode() const
  {
    return getIntParam(CbcPrinting);
  }

  /** Set the
        \link CbcModel::CbcMaximumSeconds maximum number of seconds \endlink
        desired.
    */
  inline bool setMaximumSeconds(double value)
  {
    return setDblParam(CbcMaximumSeconds, value);
  }
  /** Get the
        \link CbcModel::CbcMaximumSeconds maximum number of seconds \endlink
        desired.
    */
  inline double getMaximumSeconds() const
  {
    return getDblParam(CbcMaximumSeconds);
  }
  /// Current time since start of branchAndbound
  double getCurrentSeconds() const;

  /// Return true if maximum time reached
  bool maximumSecondsReached() const;

  /** Set the
      \link CbcModel::CbcIntegerTolerance integrality tolerance \endlink
    */
  inline bool setIntegerTolerance(double value)
  {
    return setDblParam(CbcIntegerTolerance, value);
  }
  /** Get the
      \link CbcModel::CbcIntegerTolerance integrality tolerance \endlink
    */
  inline double getIntegerTolerance() const
  {
    return getDblParam(CbcIntegerTolerance);
  }

  /** Set the
        \link CbcModel::CbcInfeasibilityWeight
          weight per integer infeasibility \endlink
    */
  inline bool setInfeasibilityWeight(double value)
  {
    return setDblParam(CbcInfeasibilityWeight, value);
  }
  /** Get the
        \link CbcModel::CbcInfeasibilityWeight
          weight per integer infeasibility \endlink
    */
  inline double getInfeasibilityWeight() const
  {
    return getDblParam(CbcInfeasibilityWeight);
  }

  /** Set the \link CbcModel::CbcAllowableGap allowable gap \endlink
        between the best known solution and the best possible solution.
    */
  inline bool setAllowableGap(double value)
  {
    return setDblParam(CbcAllowableGap, value);
  }
  /** Get the \link CbcModel::CbcAllowableGap allowable gap \endlink
        between the best known solution and the best possible solution.
    */
  inline double getAllowableGap() const
  {
    return getDblParam(CbcAllowableGap);
  }

  /** Set the \link CbcModel::CbcAllowableFractionGap fraction allowable gap \endlink
        between the best known solution and the best possible solution.
    */
  inline bool setAllowableFractionGap(double value)
  {
    return setDblParam(CbcAllowableFractionGap, value);
  }
  /** Get the \link CbcModel::CbcAllowableFractionGap fraction allowable gap \endlink
        between the best known solution and the best possible solution.
    */
  inline double getAllowableFractionGap() const
  {
    return getDblParam(CbcAllowableFractionGap);
  }
  /** Set the \link CbcModel::CbcAllowableFractionGap percentage allowable gap \endlink
        between the best known solution and the best possible solution.
    */
  inline bool setAllowablePercentageGap(double value)
  {
    return setDblParam(CbcAllowableFractionGap, value * 0.01);
  }
  /** Get the \link CbcModel::CbcAllowableFractionGap percentage allowable gap \endlink
        between the best known solution and the best possible solution.
    */
  inline double getAllowablePercentageGap() const
  {
    return 100.0 * getDblParam(CbcAllowableFractionGap);
  }
  /** Set the \link CbcModel::CbcHeuristicGap heuristic gap \endlink
        between the best known solution and the best possible solution.
    */
  inline bool setHeuristicGap(double value)
  {
    return setDblParam(CbcHeuristicGap, value);
  }
  /** Get the \link CbcModel::CbcHeuristicGap heuristic gap \endlink
        between the best known solution and the best possible solution.
    */
  inline double getHeuristicGap() const
  {
    return getDblParam(CbcHeuristicGap);
  }

  /** Set the \link CbcModel::CbcHeuristicFractionGap fraction heuristic gap \endlink
        between the best known solution and the best possible solution.
    */
  inline bool setHeuristicFractionGap(double value)
  {
    return setDblParam(CbcHeuristicFractionGap, value);
  }
  /** Get the \link CbcModel::CbcHeuristicFractionGap fraction heuristic gap \endlink
        between the best known solution and the best possible solution.
    */
  inline double getHeuristicFractionGap() const
  {
    return getDblParam(CbcHeuristicFractionGap);
  }
  /** Set the
        \link CbcModel::CbcCutoffIncrement  \endlink
        desired.
    */
  inline bool setCutoffIncrement(double value)
  {
    return setDblParam(CbcCutoffIncrement, value);
  }
  /** Get the
        \link CbcModel::CbcCutoffIncrement  \endlink
        desired.
    */
  inline double getCutoffIncrement() const
  {
    return getDblParam(CbcCutoffIncrement);
  }
  /// See if can stop on gap
  bool canStopOnGap() const;

  /** Pass in target solution and optional priorities.
        If priorities then >0 means only branch if incorrect
        while <0 means branch even if correct. +1 or -1 are
        highest priority */
  void setHotstartSolution(const double *solution, const int *priorities = NULL);

  /// Set the minimum drop to continue cuts
  inline void setMinimumDrop(double value)
  {
    minimumDrop_ = value;
  }
  /// Get the minimum drop to continue cuts
  inline double getMinimumDrop() const
  {
    return minimumDrop_;
  }

  /** Set the maximum number of cut passes at root node (default 20)
        Minimum drop can also be used for fine tuning */
  inline void setMaximumCutPassesAtRoot(int value)
  {
    maximumCutPassesAtRoot_ = value;
  }
  /** Get the maximum number of cut passes at root node */
  inline int getMaximumCutPassesAtRoot() const
  {
    return maximumCutPassesAtRoot_;
  }

  /** Set the maximum number of cut passes at other nodes (default 10)
        Minimum drop can also be used for fine tuning */
  inline void setMaximumCutPasses(int value)
  {
    maximumCutPasses_ = value;
  }
  /** Get the maximum number of cut passes at other nodes (default 10) */
  inline int getMaximumCutPasses() const
  {
    return maximumCutPasses_;
  }
  /** Get current cut pass number in this round of cuts.
        (1 is first pass) */
  inline int getCurrentPassNumber() const
  {
    return currentPassNumber_;
  }
  /** Set current cut pass number in this round of cuts.
        (1 is first pass) */
  inline void setCurrentPassNumber(int value)
  {
    currentPassNumber_ = value;
  }

  /** Set the maximum number of candidates to be evaluated for strong
      branching.

      A value of 0 disables strong branching.
    */
  void setNumberStrong(int number);
  /** Get the maximum number of candidates to be evaluated for strong
      branching.
    */
  inline int numberStrong() const
  {
    return numberStrong_;
  }
  /** Set global preferred way to branch
        -1 down, +1 up, 0 no preference */
  inline void setPreferredWay(int value)
  {
    preferredWay_ = value;
  }
  /** Get the preferred way to branch (default 0) */
  inline int getPreferredWay() const
  {
    return preferredWay_;
  }
  /// Get at which depths to do cuts
  inline int whenCuts() const
  {
    return whenCuts_;
  }
  /// Set at which depths to do cuts
  inline void setWhenCuts(int value)
  {
    whenCuts_ = value;
  }
  /** Return true if we want to do cuts
        If allowForTopOfTree zero then just does on multiples of depth
        if 1 then allows for doing at top of tree
        if 2 then says if cuts allowed anywhere apart from root
    */
  bool doCutsNow(int allowForTopOfTree) const;

  /** Set the number of branches before pseudo costs believed
        in dynamic strong branching.

      A value of 0 disables dynamic strong branching.
    */
  void setNumberBeforeTrust(int number);
  /** get the number of branches before pseudo costs believed
        in dynamic strong branching. */
  inline int numberBeforeTrust() const
  {
    return numberBeforeTrust_;
  }
  /** Set the number of variables for which to compute penalties
        in dynamic strong branching.

      A value of 0 disables penalties.
    */
  void setNumberPenalties(int number);
  /** get the number of variables for which to compute penalties
        in dynamic strong branching. */
  inline int numberPenalties() const
  {
    return numberPenalties_;
  }
  /// Pointer to top of tree
  inline const CbcFullNodeInfo *topOfTree() const
  {
    return topOfTree_;
  }
  /// Number of analyze iterations to do
  inline void setNumberAnalyzeIterations(int number)
  {
    numberAnalyzeIterations_ = number;
  }
  inline int numberAnalyzeIterations() const
  {
    return numberAnalyzeIterations_;
  }
  /** Get scale factor to make penalties match strong.
        Should/will be computed */
  inline double penaltyScaleFactor() const
  {
    return penaltyScaleFactor_;
  }
  /** Set scale factor to make penalties match strong.
        Should/will be computed */
  void setPenaltyScaleFactor(double value);
  /** Problem type as set by user or found by analysis.  This will be extended
        0 - not known
        1 - Set partitioning <=
        2 - Set partitioning ==
        3 - Set covering
        4 - all +- 1 or all +1 and odd
    */
  void inline setProblemType(int number)
  {
    problemType_ = number;
  }
  inline int problemType() const
  {
    return problemType_;
  }
  /// Current depth
  inline int currentDepth() const
  {
    return currentDepth_;
  }

  /// Set how often to scan global cuts
  void setHowOftenGlobalScan(int number);
  /// Get how often to scan global cuts
  inline int howOftenGlobalScan() const
  {
    return howOftenGlobalScan_;
  }
  /// Original columns as created by integerPresolve or preprocessing
  inline int *originalColumns() const
  {
    return originalColumns_;
  }
  /// Set original columns as created by preprocessing
  void setOriginalColumns(const int *originalColumns,
    int numberGood = COIN_INT_MAX);
  /// Create conflict cut (well - most of)
  OsiRowCut *conflictCut(const OsiSolverInterface *solver, bool &localCuts);

  /** Set the print frequency.

      Controls the number of nodes evaluated between status prints.
      If <tt>number <=0</tt> the print frequency is set to 100 nodes for large
      problems, 1000 for small problems.
      Print frequency has very slight overhead if small.
    */
  inline void setPrintFrequency(int number)
  {
    printFrequency_ = number;
  }
  /// Get the print frequency
  inline int printFrequency() const
  {
    return printFrequency_;
  }
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
  /// Was continuous solution unbounded
  bool isContinuousUnbounded() const;
  /// Was continuous solution unbounded
  bool isProvenDualInfeasible() const;
  /// Node limit reached?
  bool isNodeLimitReached() const;
  /// Time limit reached?
  bool isSecondsLimitReached() const;
  /// Solution limit reached?
  bool isSolutionLimitReached() const;
  /// Get how many iterations it took to solve the problem.
  inline int getIterationCount() const
  {
    return numberIterations_;
  }
  /// Increment how many iterations it took to solve the problem.
  inline void incrementIterationCount(int value)
  {
    numberIterations_ += value;
  }
  /// Get how many Nodes it took to solve the problem (including those in complete fathoming B&B inside CLP).
  inline int getNodeCount() const
  {
    return numberNodes_;
  }
  /// Increment how many nodes it took to solve the problem.
  inline void incrementNodeCount(int value)
  {
    numberNodes_ += value;
  }
  /// Get how many Nodes were enumerated in complete fathoming B&B inside CLP
  inline int getExtraNodeCount() const
  {
    return numberExtraNodes_;
  }
  /// Get how many times complete fathoming B&B was done
  inline int getFathomCount() const
  {
    return numberFathoms_;
  }
  /** Final status of problem
        Some of these can be found out by is...... functions
        -1 before branchAndBound
        0 finished - check isProvenOptimal or isProvenInfeasible to see if solution found
        (or check value of best solution)
        1 stopped - on maxnodes, maxsols, maxtime
        2 difficulties so run was abandoned
        (5 event user programmed event occurred)
    */
  inline int status() const
  {
    return status_;
  }
  inline void setProblemStatus(int value)
  {
    status_ = value;
  }
  /** Secondary status of problem
        -1 unset (status_ will also be -1)
        0 search completed with solution
        1 linear relaxation not feasible (or worse than cutoff)
        2 stopped on gap
        3 stopped on nodes
        4 stopped on time
        5 stopped on user event
        6 stopped on solutions
        7 linear relaxation unbounded
        8 stopped on iteration limit
    */
  inline int secondaryStatus() const
  {
    return secondaryStatus_;
  }
  inline void setSecondaryStatus(int value)
  {
    secondaryStatus_ = value;
  }
  /// Are there numerical difficulties (for initialSolve) ?
  bool isInitialSolveAbandoned() const;
  /// Is optimality proven (for initialSolve) ?
  bool isInitialSolveProvenOptimal() const;
  /// Is primal infeasiblity proven (for initialSolve) ?
  bool isInitialSolveProvenPrimalInfeasible() const;
  /// Is dual infeasiblity proven (for initialSolve) ?
  bool isInitialSolveProvenDualInfeasible() const;

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
  inline int numberRowsAtContinuous() const
  {
    return numberRowsAtContinuous_;
  }

  /// Get number of columns
  inline int getNumCols() const
  {
    return solver_->getNumCols();
  }

  /// Get number of rows
  inline int getNumRows() const
  {
    return solver_->getNumRows();
  }

  /// Get number of nonzero elements
  inline CoinBigIndex getNumElements() const
  {
    return solver_->getNumElements();
  }

  /// Number of integers in problem
  inline int numberIntegers() const
  {
    return numberIntegers_;
  }
  // Integer variables
  inline const int *integerVariable() const
  {
    return integerVariable_;
  }
  /// Whether or not integer
  inline char integerType(int i) const
  {
    assert(integerInfo_);
    assert(integerInfo_[i] == 0 || integerInfo_[i] == 1);
    return integerInfo_[i];
  }
  /// Whether or not integer
  inline const char *integerType() const
  {
    return integerInfo_;
  }

  /// Get pointer to array[getNumCols()] of column lower bounds
  inline const double *getColLower() const
  {
    return solver_->getColLower();
  }

  /// Get pointer to array[getNumCols()] of column upper bounds
  inline const double *getColUpper() const
  {
    return solver_->getColUpper();
  }

  /** Get pointer to array[getNumRows()] of row constraint senses.
        <ul>
        <li>'L': <= constraint
        <li>'E': =  constraint
        <li>'G': >= constraint
        <li>'R': ranged constraint
        <li>'N': free constraint
        </ul>
    */
  inline const char *getRowSense() const
  {
    return solver_->getRowSense();
  }

  /** Get pointer to array[getNumRows()] of rows right-hand sides
        <ul>
        <li> if rowsense()[i] == 'L' then rhs()[i] == rowupper()[i]
        <li> if rowsense()[i] == 'G' then rhs()[i] == rowlower()[i]
        <li> if rowsense()[i] == 'R' then rhs()[i] == rowupper()[i]
        <li> if rowsense()[i] == 'N' then rhs()[i] == 0.0
        </ul>
    */
  inline const double *getRightHandSide() const
  {
    return solver_->getRightHandSide();
  }

  /** Get pointer to array[getNumRows()] of row ranges.
        <ul>
        <li> if rowsense()[i] == 'R' then
        rowrange()[i] == rowupper()[i] - rowlower()[i]
        <li> if rowsense()[i] != 'R' then
        rowrange()[i] is 0.0
        </ul>
    */
  inline const double *getRowRange() const
  {
    return solver_->getRowRange();
  }

  /// Get pointer to array[getNumRows()] of row lower bounds
  inline const double *getRowLower() const
  {
    return solver_->getRowLower();
  }

  /// Get pointer to array[getNumRows()] of row upper bounds
  inline const double *getRowUpper() const
  {
    return solver_->getRowUpper();
  }

  /// Get pointer to array[getNumCols()] of objective function coefficients
  inline const double *getObjCoefficients() const
  {
    return solver_->getObjCoefficients();
  }

  /// Get objective function sense (1 for min (default), -1 for max)
  inline double getObjSense() const
  {
    //assert (dblParam_[CbcOptimizationDirection]== solver_->getObjSense());
    return dblParam_[CbcOptimizationDirection];
  }

  /// Return true if variable is continuous
  inline bool isContinuous(int colIndex) const
  {
    return solver_->isContinuous(colIndex);
  }

  /// Return true if variable is binary
  inline bool isBinary(int colIndex) const
  {
    return solver_->isBinary(colIndex);
  }

  /** Return true if column is integer.
        Note: This function returns true if the the column
        is binary or a general integer.
    */
  inline bool isInteger(int colIndex) const
  {
    return solver_->isInteger(colIndex);
  }

  /// Return true if variable is general integer
  inline bool isIntegerNonBinary(int colIndex) const
  {
    return solver_->isIntegerNonBinary(colIndex);
  }

  /// Return true if variable is binary and not fixed at either bound
  inline bool isFreeBinary(int colIndex) const
  {
    return solver_->isFreeBinary(colIndex);
  }

  /// Get pointer to row-wise copy of matrix
  inline const CoinPackedMatrix *getMatrixByRow() const
  {
    return solver_->getMatrixByRow();
  }

  /// Get pointer to column-wise copy of matrix
  inline const CoinPackedMatrix *getMatrixByCol() const
  {
    return solver_->getMatrixByCol();
  }

  /// Get solver's value for infinity
  inline double getInfinity() const
  {
    return solver_->getInfinity();
  }
  /// Get pointer to array[getNumCols()] (for speed) of column lower bounds
  inline const double *getCbcColLower() const
  {
    return cbcColLower_;
  }
  /// Get pointer to array[getNumCols()] (for speed) of column upper bounds
  inline const double *getCbcColUpper() const
  {
    return cbcColUpper_;
  }
  /// Get pointer to array[getNumRows()] (for speed) of row lower bounds
  inline const double *getCbcRowLower() const
  {
    return cbcRowLower_;
  }
  /// Get pointer to array[getNumRows()] (for speed) of row upper bounds
  inline const double *getCbcRowUpper() const
  {
    return cbcRowUpper_;
  }
  /// Get pointer to array[getNumCols()] (for speed) of primal solution vector
  inline const double *getCbcColSolution() const
  {
    return cbcColSolution_;
  }
  /// Get pointer to array[getNumRows()] (for speed) of dual prices
  inline const double *getCbcRowPrice() const
  {
    return cbcRowPrice_;
  }
  /// Get a pointer to array[getNumCols()] (for speed) of reduced costs
  inline const double *getCbcReducedCost() const
  {
    return cbcReducedCost_;
  }
  /// Get pointer to array[getNumRows()] (for speed) of row activity levels.
  inline const double *getCbcRowActivity() const
  {
    return cbcRowActivity_;
  }
  //@}

  /**@name Methods related to querying the solution */
  //@{
  /// Holds solution at continuous (after cuts if branchAndBound called)
  inline double *continuousSolution() const
  {
    return continuousSolution_;
  }
  /** Array marked whenever a solution is found if non-zero.
        Code marks if heuristic returns better so heuristic
        need only mark if it wants to on solutions which
        are worse than current */
  inline int *usedInSolution() const
  {
    return usedInSolution_;
  }
  /// Increases usedInSolution for nonzeros
  void incrementUsed(const double *solution);
  /// Record a new incumbent solution and update objectiveValue
  void setBestSolution(CBC_Message how,
    double &objectiveValue, const double *solution,
    int fixVariables = 0);
  /// Just update objectiveValue
  void setBestObjectiveValue(double objectiveValue);
  /// Deals with event handler and solution
  CbcEventHandler::CbcAction dealWithEventHandler(CbcEventHandler::CbcEvent event,
    double objValue,
    const double *solution);

  /** Call this to really test if a valid solution can be feasible
        Solution is number columns in size.
        If fixVariables true then bounds of continuous solver updated.
        Returns objective value (worse than cutoff if not feasible)
        Previously computed objective value is now passed in (in case user does not do solve)
	virtual so user can override
    */
  virtual double checkSolution(double cutoff, double *solution,
    int fixVariables, double originalObjValue);
  /** Test the current solution for feasiblility.

      Scan all objects for indications of infeasibility. This is broken down
      into simple integer infeasibility (\p numberIntegerInfeasibilities)
      and all other reports of infeasibility (\p numberObjectInfeasibilities).
    */
  bool feasibleSolution(int &numberIntegerInfeasibilities,
    int &numberObjectInfeasibilities) const;

  /** Solution to the most recent lp relaxation.

      The solver's solution to the most recent lp relaxation.
    */

  inline double *currentSolution() const
  {
    return currentSolution_;
  }
  /** For testing infeasibilities - will point to
        currentSolution_ or solver-->getColSolution()
    */
  inline const double *testSolution() const
  {
    return testSolution_;
  }
  inline void setTestSolution(const double *solution)
  {
    testSolution_ = solution;
  }
  /// Make sure region there and optionally copy solution
  void reserveCurrentSolution(const double *solution = NULL);

  /// Get pointer to array[getNumCols()] of primal solution vector
  inline const double *getColSolution() const
  {
    return solver_->getColSolution();
  }

  /// Get pointer to array[getNumRows()] of dual prices
  inline const double *getRowPrice() const
  {
    return solver_->getRowPrice();
  }

  /// Get a pointer to array[getNumCols()] of reduced costs
  inline const double *getReducedCost() const
  {
    return solver_->getReducedCost();
  }

  /// Get pointer to array[getNumRows()] of row activity levels.
  inline const double *getRowActivity() const
  {
    return solver_->getRowActivity();
  }

  /// Get current objective function value
  inline double getCurrentObjValue() const
  {
    return dblParam_[CbcCurrentObjectiveValue];
  }
  /// Get current minimization objective function value
  inline double getCurrentMinimizationObjValue() const
  {
    return dblParam_[CbcCurrentMinimizationObjectiveValue];
  }

  /// Get best objective function value as minimization
  inline double getMinimizationObjValue() const
  {
    return bestObjective_;
  }
  /// Set best objective function value as minimization
  inline void setMinimizationObjValue(double value)
  {
    bestObjective_ = value;
  }

  /// Get best objective function value
  inline double getObjValue() const
  {
    return bestObjective_ * solver_->getObjSense();
  }
  /** Get best possible objective function value.
        This is better of best possible left on tree
        and best solution found.
        If called from within branch and cut may be optimistic.
    */
  double getBestPossibleObjValue() const;
  /// Set best objective function value
  inline void setObjValue(double value)
  {
    bestObjective_ = value * solver_->getObjSense();
  }
  /// Get solver objective function value (as minimization)
  inline double getSolverObjValue() const
  {
    return solver_->getObjValue() * solver_->getObjSense();
  }

  /** The best solution to the integer programming problem.

      The best solution to the integer programming problem found during
      the search. If no solution is found, the method returns null.
    */

  inline double *bestSolution() const
  {
    return bestSolution_;
  }
  /** User callable setBestSolution.
        If check false does not check valid
        If true then sees if feasible and warns if objective value
        worse than given (so just set to COIN_DBL_MAX if you don't care).
        If check true then does not save solution if not feasible
    */
  void setBestSolution(const double *solution, int numberColumns,
    double objectiveValue, bool check = false);

  /// Get number of solutions
  inline int getSolutionCount() const
  {
    return numberSolutions_;
  }

  /// Set number of solutions (so heuristics will be different)
  inline void setSolutionCount(int value)
  {
    numberSolutions_ = value;
  }
  /// Number of saved solutions (including best)
  int numberSavedSolutions() const;
  /// Maximum number of extra saved solutions
  inline int maximumSavedSolutions() const
  {
    return maximumSavedSolutions_;
  }
  /// Set maximum number of extra saved solutions
  void setMaximumSavedSolutions(int value);
  /// Return a saved solution (0==best) - NULL if off end
  const double *savedSolution(int which) const;
  /// Return a saved solution objective (0==best) - COIN_DBL_MAX if off end
  double savedSolutionObjective(int which) const;
  /// Delete a saved solution and move others up
  void deleteSavedSolution(int which);

  /** Current phase (so heuristics etc etc can find out).
        0 - initial solve
        1 - solve with cuts at root
        2 - solve with cuts
        3 - other e.g. strong branching
        4 - trying to validate a solution
        5 - at end of search
    */
  inline int phase() const
  {
    return phase_;
  }

  /// Get number of heuristic solutions
  inline int getNumberHeuristicSolutions() const
  {
    return numberHeuristicSolutions_;
  }
  /// Set number of heuristic solutions
  inline void setNumberHeuristicSolutions(int value)
  {
    numberHeuristicSolutions_ = value;
  }

  /// Set objective function sense (1 for min (default), -1 for max,)
  inline void setObjSense(double s)
  {
    dblParam_[CbcOptimizationDirection] = s;
    solver_->setObjSense(s);
  }

  /// Value of objective at continuous
  inline double getContinuousObjective() const
  {
    return originalContinuousObjective_;
  }
  inline void setContinuousObjective(double value)
  {
    originalContinuousObjective_ = value;
  }
  /// Number of infeasibilities at continuous
  inline int getContinuousInfeasibilities() const
  {
    return continuousInfeasibilities_;
  }
  inline void setContinuousInfeasibilities(int value)
  {
    continuousInfeasibilities_ = value;
  }
  /// Value of objective after root node cuts added
  inline double rootObjectiveAfterCuts() const
  {
    return continuousObjective_;
  }
  /// Sum of Changes to objective by first solve
  inline double sumChangeObjective() const
  {
    return sumChangeObjective1_;
  }
  /** Number of times global cuts violated.  When global cut pool then this
        should be kept for each cut and type of cut */
  inline int numberGlobalViolations() const
  {
    return numberGlobalViolations_;
  }
  inline void clearNumberGlobalViolations()
  {
    numberGlobalViolations_ = 0;
  }
  /// Whether to force a resolve after takeOffCuts
  inline bool resolveAfterTakeOffCuts() const
  {
    return resolveAfterTakeOffCuts_;
  }
  inline void setResolveAfterTakeOffCuts(bool yesNo)
  {
    resolveAfterTakeOffCuts_ = yesNo;
  }
  /// Maximum number of rows
  inline int maximumRows() const
  {
    return maximumRows_;
  }
  /// Work basis for temporary use
  inline CoinWarmStartBasis &workingBasis()
  {
    return workingBasis_;
  }
  /// Get number of "iterations" to stop after
  inline int getStopNumberIterations() const
  {
    return stopNumberIterations_;
  }
  /// Set number of "iterations" to stop after
  inline void setStopNumberIterations(int value)
  {
    stopNumberIterations_ = value;
  }
  /// A pointer to model from CbcHeuristic
  inline CbcModel *heuristicModel() const
  {
    return heuristicModel_;
  }
  /// Set a pointer to model from CbcHeuristic
  inline void setHeuristicModel(CbcModel *model)
  {
    heuristicModel_ = model;
  }
  //@}

  /** \name Node selection */
  //@{
  // Comparison functions (which may be overridden by inheritance)
  inline CbcCompareBase *nodeComparison() const
  {
    return nodeCompare_;
  }
  void setNodeComparison(CbcCompareBase *compare);
  void setNodeComparison(CbcCompareBase &compare);
  //@}

  /** \name Problem feasibility checking */
  //@{
  // Feasibility functions (which may be overridden by inheritance)
  inline CbcFeasibilityBase *problemFeasibility() const
  {
    return problemFeasibility_;
  }
  void setProblemFeasibility(CbcFeasibilityBase *feasibility);
  void setProblemFeasibility(CbcFeasibilityBase &feasibility);
  //@}

  /** \name Tree methods and subtree methods */
  //@{
  /// Tree method e.g. heap (which may be overridden by inheritance)
  inline CbcTree *tree() const
  {
    return tree_;
  }
  /// For modifying tree handling (original is cloned)
  void passInTreeHandler(CbcTree &tree);
  /** For passing in an CbcModel to do a sub Tree (with derived tree handlers).
        Passed in model must exist for duration of branch and bound
    */
  void passInSubTreeModel(CbcModel &model);
  /** For retrieving a copy of subtree model with given OsiSolver.
        If no subtree model will use self (up to user to reset cutoff etc).
        If solver NULL uses current
    */
  CbcModel *subTreeModel(OsiSolverInterface *solver = NULL) const;
  /// Returns number of times any subtree stopped on nodes, time etc
  inline int numberStoppedSubTrees() const
  {
    return numberStoppedSubTrees_;
  }
  /// Says a sub tree was stopped
  inline void incrementSubTreeStopped()
  {
    numberStoppedSubTrees_++;
  }
  /** Whether to automatically do presolve before branch and bound (subTrees).
        0 - no
        1 - ordinary presolve
        2 - integer presolve (dodgy)
    */
  inline int typePresolve() const
  {
    return presolve_;
  }
  inline void setTypePresolve(int value)
  {
    presolve_ = value;
  }

  //@}

  /** \name Branching Decisions

      See the CbcBranchDecision class for additional information.
    */
  //@{

  /// Get the current branching decision method.
  inline CbcBranchDecision *branchingMethod() const
  {
    return branchingMethod_;
  }
  /// Set the branching decision method.
  inline void setBranchingMethod(CbcBranchDecision *method)
  {
    delete branchingMethod_;
    branchingMethod_ = method->clone();
  }
  /** Set the branching method

      \overload
    */
  inline void setBranchingMethod(CbcBranchDecision &method)
  {
    delete branchingMethod_;
    branchingMethod_ = method.clone();
  }
  /// Get the current cut modifier method
  inline CbcCutModifier *cutModifier() const
  {
    return cutModifier_;
  }
  /// Set the cut modifier method
  void setCutModifier(CbcCutModifier *modifier);
  /** Set the cut modifier method

      \overload
    */
  void setCutModifier(CbcCutModifier &modifier);
  //@}

  /** \name Row (constraint) and Column (variable) cut generation */
  //@{

  /** State of search
        0 - no solution
        1 - only heuristic solutions
        2 - branched to a solution
        3 - no solution but many nodes
    */
  inline int stateOfSearch() const
  {
    return stateOfSearch_;
  }
  inline void setStateOfSearch(int state)
  {
    stateOfSearch_ = state;
  }
  /// Strategy worked out - mainly at root node for use by CbcNode
  inline int searchStrategy() const
  {
    return searchStrategy_;
  }
  /// Set strategy worked out - mainly at root node for use by CbcNode
  inline void setSearchStrategy(int value)
  {
    searchStrategy_ = value;
  }
  /// Stong branching strategy
  inline int strongStrategy() const
  {
    return strongStrategy_;
  }
  /// Set strong branching strategy
  inline void setStrongStrategy(int value)
  {
    strongStrategy_ = value;
  }

  /// Get the number of cut generators
  inline int numberCutGenerators() const
  {
    return numberCutGenerators_;
  }
  /// Get the list of cut generators
  inline CbcCutGenerator **cutGenerators() const
  {
    return generator_;
  }
  ///Get the specified cut generator
  inline CbcCutGenerator *cutGenerator(int i) const
  {
    return generator_[i];
  }
  ///Get the specified cut generator before any changes
  inline CbcCutGenerator *virginCutGenerator(int i) const
  {
    return virginGenerator_[i];
  }
  /** Add one generator - up to user to delete generators.
        howoften affects how generator is used. 0 or 1 means always,
        >1 means every that number of nodes.  Negative values have same
        meaning as positive but they may be switched off (-> -100) by code if
        not many cuts generated at continuous.  -99 is just done at root.
        Name is just for printout.
        If depth >0 overrides how often generator is called (if howOften==-1 or >0).
    */
  void addCutGenerator(CglCutGenerator *generator,
    int howOften = 1, const char *name = NULL,
    bool normal = true, bool atSolution = false,
    bool infeasible = false, int howOftenInSub = -100,
    int whatDepth = -1, int whatDepthInSub = -1);
  //@}
  /** \name Strategy and sub models

      See the CbcStrategy class for additional information.
    */
  //@{

  /// Get the current strategy
  inline CbcStrategy *strategy() const
  {
    return strategy_;
  }
  /// Set the strategy. Clones
  void setStrategy(CbcStrategy &strategy);
  /// Set the strategy. assigns
  inline void setStrategy(CbcStrategy *strategy)
  {
    strategy_ = strategy;
  }
  /// Get the current parent model
  inline CbcModel *parentModel() const
  {
    return parentModel_;
  }
  /// Set the parent model
  inline void setParentModel(CbcModel &parentModel)
  {
    parentModel_ = &parentModel;
  }
  //@}

  /** \name Heuristics and priorities */
  //@{
  /*! \brief Add one heuristic - up to user to delete

      The name is just used for print messages.
    */
  void addHeuristic(CbcHeuristic *generator, const char *name = NULL,
    int before = -1);
  ///Get the specified heuristic
  inline CbcHeuristic *heuristic(int i) const
  {
    return heuristic_[i];
  }
  /// Get the number of heuristics
  inline int numberHeuristics() const
  {
    return numberHeuristics_;
  }
  /// Set the number of heuristics
  inline void setNumberHeuristics(int value)
  {
    numberHeuristics_ = value;
  }
  /// Pointer to heuristic solver which found last solution (or NULL)
  inline CbcHeuristic *lastHeuristic() const
  {
    return lastHeuristic_;
  }
  /// set last heuristic which found a solution
  inline void setLastHeuristic(CbcHeuristic *last)
  {
    lastHeuristic_ = last;
  }

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

        This is now deprecated except for simple usage.  If user
        creates Cbcobjects then set priority in them

        \internal Added for Kurt Spielberg.
    */
  void passInPriorities(const int *priorities, bool ifNotSimpleIntegers);

  /// Returns priority level for an object (or 1000 if no priorities exist)
  inline int priority(int sequence) const
  {
    return object_[sequence]->priority();
  }

  /*! \brief Set an event handler

      A clone of the handler passed as a parameter is stored in CbcModel.
    */
  void passInEventHandler(const CbcEventHandler *eventHandler);

  /*! \brief Retrieve a pointer to the event handler */
  inline CbcEventHandler *getEventHandler() const
  {
    return (eventHandler_);
  }

  //@}

  /**@name Setting/Accessing application data */
  //@{
  /** Set application data.

    This is a pointer that the application can store into and
    retrieve from the solver interface.
    This field is available for the application to optionally
    define and use.
    */
  void setApplicationData(void *appData);

  /// Get application data
  void *getApplicationData() const;
  /**
        For advanced applications you may wish to modify the behavior of Cbc
        e.g. if the solver is a NLP solver then you may not have an exact
        optimum solution at each step.  Information could be built into
        OsiSolverInterface but this is an alternative so that that interface
        does not have to be changed.  If something similar is useful to
        enough solvers then it could be migrated
        You can also pass in by using solver->setAuxiliaryInfo.
        You should do that if solver is odd - if solver is normal simplex
        then use this.
        NOTE - characteristics are not cloned
    */
  void passInSolverCharacteristics(OsiBabSolver *solverCharacteristics);
  /// Get solver characteristics
  inline const OsiBabSolver *solverCharacteristics() const
  {
    return solverCharacteristics_;
  }
  //@}

  //---------------------------------------------------------------------------

  /**@name Message handling etc */
  //@{
  /// Pass in Message handler (not deleted at end)
  void passInMessageHandler(CoinMessageHandler *handler);
  /// Set language
  void newLanguage(CoinMessages::Language language);
  inline void setLanguage(CoinMessages::Language language)
  {
    newLanguage(language);
  }
  /// Return handler
  inline CoinMessageHandler *messageHandler() const
  {
    return handler_;
  }
  /// Return messages
  inline CoinMessages &messages()
  {
    return messages_;
  }
  /// Return pointer to messages
  inline CoinMessages *messagesPointer()
  {
    return &messages_;
  }
  /// Set log level
  void setLogLevel(int value);
  /// Get log level
  inline int logLevel() const
  {
    return handler_->logLevel();
  }
  /** Set flag to say if handler_ is the default handler.

      The default handler is deleted when the model is deleted. Other
      handlers (supplied by the client) will not be deleted.
    */
  inline void setDefaultHandler(bool yesNo)
  {
    defaultHandler_ = yesNo;
  }
  /// Check default handler
  inline bool defaultHandler() const
  {
    return defaultHandler_;
  }
  //@}
  //---------------------------------------------------------------------------
  ///@name Specialized
  //@{

  /**
        Set special options
        0 bit (1) - check if cuts valid (if on debugger list)
        1 bit (2) - use current basis to check integer solution (rather than all slack)
        2 bit (4) - don't check integer solution (by solving LP)
        3 bit (8) - fast analyze
        4 bit (16) - non-linear model - so no well defined CoinPackedMatrix
        5 bit (32) - keep names
        6 bit (64) - try for dominated columns
        7 bit (128) - SOS type 1 but all declared integer
        8 bit (256) - Set to say solution just found, unset by doing cuts
        9 bit (512) - Try reduced model after 100 nodes
        10 bit (1024) - Switch on some heuristics even if seems unlikely
        11 bit (2048) - Mark as in small branch and bound
        12 bit (4096) - Funny cuts so do slow way (in some places)
        13 bit (8192) - Funny cuts so do slow way (in other places)
        14 bit (16384) - Use Cplex! for fathoming
        15 bit (32768) - Try reduced model after 0 nodes
        16 bit (65536) - Original model had integer bounds
        17 bit (131072) - Perturbation switched off
        18 bit (262144) - donor CbcModel
        19 bit (524288) - recipient CbcModel
        20 bit (1048576) - waiting for sub model to return
	22 bit (4194304) - do not initialize random seed in solver (user has)
	23 bit (8388608) - leave solver_ with cuts
	24 bit (16777216) - just get feasible if no cutoff
	25 bit (33554432) - feasibility pump after root cuts
	26 bit (67108864) - child model but going for complete search
    */
  inline void setSpecialOptions(int value)
  {
    specialOptions_ = value;
  }
  /// Get special options
  inline int specialOptions() const
  {
    return specialOptions_;
  }
  /// Set random seed
  inline void setRandomSeed(int value)
  {
    randomSeed_ = value;
  }
  /// Get random seed
  inline int getRandomSeed() const
  {
    return randomSeed_;
  }
  /// Set multiple root tries
  inline void setMultipleRootTries(int value)
  {
    multipleRootTries_ = value;
  }
  /// Get multiple root tries
  inline int getMultipleRootTries() const
  {
    return multipleRootTries_;
  }
  /// Tell model to stop on event
  inline void sayEventHappened()
  {
    eventHappened_ = true;
  }
  /// Says if normal solver i.e. has well defined CoinPackedMatrix
  inline bool normalSolver() const
  {
    return (specialOptions_ & 16) == 0;
  }
  /** Says if model is sitting there waiting for mini branch and bound to finish
	This is because an event handler may only have access to parent model in
	mini branch and bound
    */
  inline bool waitingForMiniBranchAndBound() const
  {
    return (specialOptions_ & 1048576) != 0;
  }
  /** Set more special options
        at present bottom 6 bits used for shadow price mode
        1024 for experimental hotstart
        2048,4096 breaking out of cuts
        8192 slowly increase minimum drop
        16384 gomory
	32768 more heuristics in sub trees
	65536 no cuts in preprocessing
        131072 Time limits elapsed
        18 bit (262144) - Perturb fathom nodes
        19 bit (524288) - No limit on fathom nodes
        20 bit (1048576) - Reduce sum of infeasibilities before cuts
        21 bit (2097152) - Reduce sum of infeasibilities after cuts
	22 bit (4194304) - Conflict analysis
	23 bit (8388608) - Conflict analysis - temporary bit
	24 bit (16777216) - Add cutoff as LP constraint (out)
	25 bit (33554432) - diving/reordering
	26 bit (67108864) - load global cuts from file
	27 bit (134217728) - append binding global cuts to file
	28 bit (268435456) - idiot branching
        29 bit (536870912) - don't make fake objective
	30 bit (1073741824) - Funny SOS or similar - be careful
    */
  inline void setMoreSpecialOptions(int value)
  {
    moreSpecialOptions_ = value;
  }
  /// Get more special options
  inline int moreSpecialOptions() const
  {
    return moreSpecialOptions_;
  }
  /** Set more more special options
	0 bit (1) - find switching variables
	1 bit (2) - using fake objective until solution
	2 bit (4) - switching variables exist
	3 bit (8) - skip most of setBestSolution checks
	4 bit (16) - very lightweight preprocessing in smallB&B
	5 bit (32) - event handler needs to be cloned when parallel
	6 bit (64) - testing - use probing to make cliques
	7/8 bit (128) - try orbital branching (if nauty)
	9 bit (512) - branching on objective (later)
	10 bit (1024) - branching on constraints (later)
	11/12 bit 2048 - intermittent cuts
	13/14 bit 8192 - go to bitter end in strong branching (first time)
	15 bit 32768 - take care of very very small values for Integer/SOS variables
    */
  inline void setMoreSpecialOptions2(int value)
  {
    moreSpecialOptions2_ = value;
  }
  /// Get more special options2
  inline int moreSpecialOptions2() const
  {
    return moreSpecialOptions2_;
  }
  /// Set cutoff as constraint
  inline void setCutoffAsConstraint(bool yesNo)
  {
    cutoffRowNumber_ = (yesNo) ? -2 : -1;
  }
  /// Set time method
  inline void setUseElapsedTime(bool yesNo)
  {
    if (yesNo)
      moreSpecialOptions_ |= 131072;
    else
      moreSpecialOptions_ &= ~131072;
  }
  /// Get time method
  inline bool useElapsedTime() const
  {
    return (moreSpecialOptions_ & 131072) != 0;
  }
  /// Get useful temporary pointer
  inline void *temporaryPointer() const
  {
    return temporaryPointer_;
  }
  /// Set useful temporary pointer
  inline void setTemporaryPointer(void *pointer)
  {
    temporaryPointer_ = pointer;
  }
  /// Go to dantzig pivot selection if easy problem (clp only)
  void goToDantzig(int numberNodes, ClpDualRowPivot *&savePivotMethod);
  /// Now we may not own objects - just point to solver's objects
  inline bool ownObjects() const
  {
    return ownObjects_;
  }
  /// Check original model before it gets messed up
  void checkModel();
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
      If deleteSolver then current solver deleted (if model owned)

      \note Parameter settings in the outgoing solver are not inherited by
        the incoming solver.
    */
  void assignSolver(OsiSolverInterface *&solver, bool deleteSolver = true);

  /** \brief Set ownership of solver

      A parameter of false tells CbcModel it does not own the solver and
      should not delete it. Once you claim ownership of the solver, you're
      responsible for eventually deleting it. Note that CbcModel clones
      solvers with abandon.  Unless you have a deep understanding of the
      workings of CbcModel, the only time you want to claim ownership is when
      you're about to delete the CbcModel object but want the solver to
      continue to exist (as, for example, when branchAndBound has finished
      and you want to hang on to the answer).
    */
  inline void setModelOwnsSolver(bool ourSolver)
  {
    ownership_ = ourSolver ? (ownership_ | 0x80000000) : (ownership_ & (~0x80000000));
  }

  /*! \brief Get ownership of solver

      A return value of true means that CbcModel owns the solver and will
      take responsibility for deleting it when that becomes necessary.
    */
  inline bool modelOwnsSolver()
  {
    return ((ownership_ & 0x80000000) != 0);
  }

  /** Copy constructor .
      If cloneHandler is true then message handler is cloned
    */
  CbcModel(const CbcModel &rhs, bool cloneHandler = false);

  /** Clone */
  virtual CbcModel *clone(bool cloneHandler);

  /// Assignment operator
  CbcModel &operator=(const CbcModel &rhs);

  /// Destructor
  virtual ~CbcModel();

  /// Returns solver - has current state
  inline OsiSolverInterface *solver() const
  {
    return solver_;
  }

  /// Returns current solver - sets new one
  inline OsiSolverInterface *swapSolver(OsiSolverInterface *solver)
  {
    OsiSolverInterface *returnSolver = solver_;
    solver_ = solver;
    return returnSolver;
  }

  /// Returns solver with continuous state
  inline OsiSolverInterface *continuousSolver() const
  {
    return continuousSolver_;
  }

  /// Create solver with continuous state
  inline void createContinuousSolver()
  {
    continuousSolver_ = solver_->clone();
  }
  /// Clear solver with continuous state
  inline void clearContinuousSolver()
  {
    delete continuousSolver_;
    continuousSolver_ = NULL;
  }

  /// A copy of the solver, taken at constructor or by saveReferenceSolver
  inline OsiSolverInterface *referenceSolver() const
  {
    return referenceSolver_;
  }

  /// Save a copy of the current solver so can be reset to
  void saveReferenceSolver();

  /** Uses a copy of reference solver to be current solver.
        Because of possible mismatches all exotic integer information is loat
        (apart from normal information in OsiSolverInterface)
        so SOS etc and priorities will have to be redone
    */
  void resetToReferenceSolver();

  /// Clears out as much as possible (except solver)
  void gutsOfDestructor();
  /** Clears out enough to reset CbcModel as if no branch and bound done
     */
  void gutsOfDestructor2();
  /** Clears out enough to reset CbcModel cutoff etc
     */
  void resetModel();
  /** Most of copy constructor
        mode - 0 copy but don't delete before
               1 copy and delete before
           2 copy and delete before (but use virgin generators)
    */
  void gutsOfCopy(const CbcModel &rhs, int mode = 0);
  /// Move status, nodes etc etc across
  void moveInfo(const CbcModel &rhs);
  //@}

  ///@name Multithreading
  //@{
  /// Indicates whether Cbc library has been compiled with multithreading support
  static bool haveMultiThreadSupport();
  /// Get pointer to masterthread
  CbcThread *masterThread() const
  {
    return masterThread_;
  }
  /// Get pointer to walkback
  CbcNodeInfo **walkback() const
  {
    return walkback_;
  }
  /// Get number of threads
  inline int getNumberThreads() const
  {
    return numberThreads_;
  }
  /// Set number of threads
  inline void setNumberThreads(int value)
  {
    numberThreads_ = value;
  }
  /// Get thread mode
  inline int getThreadMode() const
  {
    return threadMode_;
  }
  /** Set thread mode
        always use numberThreads for branching
        1 set then deterministic
        2 set then use numberThreads for root cuts
        4 set then use numberThreads in root mini branch and bound
        8 set and numberThreads - do heuristics numberThreads at a time
        8 set and numberThreads==0 do all heuristics at once
        default is 0
    */
  inline void setThreadMode(int value)
  {
    threadMode_ = value;
  }
  /** Return
        -2 if deterministic threaded and main thread
        -1 if deterministic threaded and serial thread
        0 if serial
        1 if opportunistic threaded
    */
  inline int parallelMode() const
  {
    if (!numberThreads_) {
      if ((threadMode_ & 1) == 0)
        return 0;
      else
        return -1;
      return 0;
    } else {
      if ((threadMode_ & 1) == 0)
        return 1;
      else
        return -2;
    }
  }
  /// Thread stuff for master
  inline CbcBaseModel *master() const
  {
    return master_;
  }
  /// From here to end of section - code in CbcThread.cpp until class changed
  /// Returns true if locked
  bool isLocked() const;
  /**
       Locks a thread if parallel so that stuff like cut pool
       can be updated and/or used.
    */
  void lockThread();
  /**
       Unlocks a thread if parallel to say cut pool stuff not needed
    */
  void unlockThread();
  /** Set information in a child
        -3 pass pointer to child thread info
        -2 just stop
        -1 delete simple child stuff
        0 delete opportunistic child stuff
        1 delete deterministic child stuff
    */
  void setInfoInChild(int type, CbcThread *info);
  /** Move/copy information from one model to another
        -1 - initialization
        0 - from base model
        1 - to base model (and reset)
        2 - add in final statistics etc (and reset so can do clean destruction)
    */
  void moveToModel(CbcModel *baseModel, int mode);
  /// Split up nodes
  int splitModel(int numberModels, CbcModel **model,
    int numberNodes);
  /// Start threads
  void startSplitModel(int numberIterations);
  /// Merge models
  void mergeModels(int numberModel, CbcModel **model,
    int numberNodes);
  //@}

  ///@name semi-private i.e. users should not use
  //@{
  /// Get how many Nodes it took to solve the problem.
  int getNodeCount2() const
  {
    return numberNodes2_;
  }
  /// Set pointers for speed
  void setPointers(const OsiSolverInterface *solver);
  /** Perform reduced cost fixing

      Fixes integer variables at their current value based on reduced cost
      penalties.  Returns number fixed
    */
  int reducedCostFix();
  /** Makes all handlers same.  If makeDefault 1 then makes top level
        default and rest point to that.  If 2 then each is copy
    */
  void synchronizeHandlers(int makeDefault);
  /// Save a solution to saved list
  void saveExtraSolution(const double *solution, double objectiveValue);
  /// Save a solution to best and move current to saved
  void saveBestSolution(const double *solution, double objectiveValue);
  /// Delete best and saved solutions
  void deleteSolutions();
  /// Encapsulates solver resolve
  int resolve(OsiSolverInterface *solver);
#ifdef CLP_RESOLVE
  /// Special purpose resolve
  int resolveClp(OsiClpSolverInterface *solver, int type);
#endif

  /** Encapsulates choosing a variable -
        anyAction -2, infeasible (-1 round again), 0 done
    */
  int chooseBranch(CbcNode *&newNode, int numberPassesLeft,
    CbcNode *oldNode, OsiCuts &cuts,
    bool &resolved, CoinWarmStartBasis *lastws,
    const double *lowerBefore, const double *upperBefore,
    OsiSolverBranch *&branches);

  /** Return an empty basis object of the specified size

      A useful utility when constructing a basis for a subproblem from scratch.
      The object returned will be of the requested capacity and appropriate for
      the solver attached to the model.
    */
  CoinWarmStartBasis *getEmptyBasis(int ns = 0, int na = 0) const;

  /** Remove inactive cuts from the model

      An OsiSolverInterface is expected to maintain a valid basis, but not a
      valid solution, when loose cuts are deleted. Restoring a valid solution
      requires calling the solver to reoptimise. If it's certain the solution
      will not be required, set allowResolve to false to suppress
      reoptimisation.
      If saveCuts then slack cuts will be saved
      On input current cuts are cuts and newCuts
      on exit current cuts will be correct.  Returns number dropped
    */
  int takeOffCuts(OsiCuts &cuts,
    bool allowResolve, OsiCuts *saveCuts,
    int numberNewCuts = 0, const OsiRowCut **newCuts = NULL);

  /** Determine and install the active cuts that need to be added for
      the current subproblem

      The whole truth is a bit more complicated. The first action is a call to
      addCuts1(). addCuts() then sorts through the list, installs the tight
      cuts in the model, and does bookkeeping (adjusts reference counts).
      The basis returned from addCuts1() is adjusted accordingly.

      If it turns out that the node should really be fathomed by bound,
      addCuts() simply treats all the cuts as loose as it does the bookkeeping.

    */
  int addCuts(CbcNode *node, CoinWarmStartBasis *&lastws);

  /** Traverse the tree from node to root and prep the model

      addCuts1() begins the job of prepping the model to match the current
      subproblem. The model is stripped of all cuts, and the search tree is
      traversed from node to root to determine the changes required. Appropriate
      bounds changes are installed, a list of cuts is collected but not
      installed, and an appropriate basis (minus the cuts, but big enough to
      accommodate them) is constructed.

      Returns true if new problem similar to old

      \todo addCuts1() is called in contexts where it's known in advance that
        all that's desired is to determine a list of cuts and do the
        bookkeeping (adjust the reference counts). The work of installing
        bounds and building a basis goes to waste.
    */
  bool addCuts1(CbcNode *node, CoinWarmStartBasis *&lastws);
  /** Returns bounds just before where - initially original bounds.
        Also sets downstream nodes (lower if force 1, upper if 2)
    */
  void previousBounds(CbcNode *node, CbcNodeInfo *where, int iColumn,
    double &lower, double &upper, int force);
  /** Set objective value in a node.  This is separated out so that
       odd solvers can use.  It may look at extra information in
       solverCharacteriscs_ and will also use bound from parent node
    */
  void setObjectiveValue(CbcNode *thisNode, const CbcNode *parentNode) const;

  /** If numberBeforeTrust >0 then we are going to use CbcBranchDynamic.
        Scan and convert CbcSimpleInteger objects
    */
  void convertToDynamic();
  /// Set numberBeforeTrust in all objects
  void synchronizeNumberBeforeTrust(int type = 0);
  /// Zap integer information in problem (may leave object info)
  void zapIntegerInformation(bool leaveObjects = true);
  /// Fill in useful estimates
  void pseudoShadow(int type);
  /** Return pseudo costs
        If not all integers or not pseudo costs - returns all zero
        Length of arrays are numberIntegers() and entries
        correspond to integerVariable()[i]
        User must allocate arrays before call
    */
  void fillPseudoCosts(double *downCosts, double *upCosts,
    int *priority = NULL,
    int *numberDown = NULL, int *numberUp = NULL,
    int *numberDownInfeasible = NULL,
    int *numberUpInfeasible = NULL) const;
  /** Do heuristics at root.
        0 - don't delete
        1 - delete
        2 - just delete - don't even use
    */
  void doHeuristicsAtRoot(int deleteHeuristicsAfterwards = 0);
  /// Adjust heuristics based on model
  void adjustHeuristics();
  /// Get the hotstart solution
  inline const double *hotstartSolution() const
  {
    return hotstartSolution_;
  }
  /// Get the hotstart priorities
  inline const int *hotstartPriorities() const
  {
    return hotstartPriorities_;
  }

  /// Return the list of cuts initially collected for this subproblem
  inline CbcCountRowCut **addedCuts() const
  {
    return addedCuts_;
  }
  /// Number of entries in the list returned by #addedCuts()
  inline int currentNumberCuts() const
  {
    return currentNumberCuts_;
  }
  /// Global cuts
  inline CbcRowCuts *globalCuts()
  {
    return &globalCuts_;
  }
  /// Get rid of global cuts
  inline void zapGlobalCuts()
  {
    globalCuts_ = CbcRowCuts();
  }
  /// Copy and set a pointer to a row cut which will be added instead of normal branching.
  void setNextRowCut(const OsiRowCut &cut);
  /// Get a pointer to current node (be careful)
  inline CbcNode *currentNode() const
  {
    return currentNode_;
  }
  /// Delete a node and possibly null out currentNode_
  void deleteNode(CbcNode * node);
  /// Get a pointer to probing info
  inline CglTreeProbingInfo *probingInfo() const
  {
    return probingInfo_;
  }
  /// Thread specific random number generator
  inline CoinThreadRandom *randomNumberGenerator()
  {
    return &randomNumberGenerator_;
  }
  /// Set the number of iterations done in strong branching.
  inline void setNumberStrongIterations(int number)
  {
    numberStrongIterations_ = number;
  }
  /// Get the number of iterations done in strong branching.
  inline int numberStrongIterations() const
  {
    return numberStrongIterations_;
  }
  /// Get maximum number of iterations (designed to be used in heuristics)
  inline int maximumNumberIterations() const
  {
    return maximumNumberIterations_;
  }
  /// Set maximum number of iterations (designed to be used in heuristics)
  inline void setMaximumNumberIterations(int value)
  {
    maximumNumberIterations_ = value;
  }
#ifdef COIN_HAS_NTY
  /// Symmetry information
  inline CbcSymmetry *symmetryInfo() const
  {
    return symmetryInfo_;
  }
  /// get rid of all
  void zapSymmetry();
#endif
  /// Set depth for fast nodes
  inline void setFastNodeDepth(int value)
  {
    fastNodeDepth_ = value;
  }
  /// Get depth for fast nodes
  inline int fastNodeDepth() const
  {
    return fastNodeDepth_;
  }
  /// Get anything with priority >= this can be treated as continuous
  inline int continuousPriority() const
  {
    return continuousPriority_;
  }
  /// Set anything with priority >= this can be treated as continuous
  inline void setContinuousPriority(int value)
  {
    continuousPriority_ = value;
  }
  inline void incrementExtra(int nodes, int iterations, int fathoms = 1)
  {
    numberExtraNodes_ += nodes;
    numberExtraIterations_ += iterations;
    numberFathoms_ += fathoms;
  }
  /// Zero extra
  inline void zeroExtra()
  {
    numberExtraNodes_ = 0;
    numberExtraIterations_ = 0;
    numberFathoms_ = 0;
  }
  /// Number of extra iterations
  inline int numberExtraIterations() const
  {
    return numberExtraIterations_;
  }
  /// Increment strong info
  void incrementStrongInfo(int numberTimes, int numberIterations,
    int numberFixed, bool ifInfeasible);
  /// Return strong info
  inline const int *strongInfo() const
  {
    return strongInfo_;
  }

  /// Return mutable strong info
  inline int *mutableStrongInfo()
  {
    return strongInfo_;
  }
  /// Get stored row cuts for donor/recipient CbcModel
  CglStored *storedRowCuts() const
  {
    return storedRowCuts_;
  }
  /// Set stored row cuts for donor/recipient CbcModel
  void setStoredRowCuts(CglStored *cuts)
  {
    storedRowCuts_ = cuts;
  }
  /// Says whether all dynamic integers
  inline bool allDynamic() const
  {
    return ((ownership_ & 0x40000000) != 0);
  }
  /// Create C++ lines to get to current state
  void generateCpp(FILE *fp, int options);
  /// Generate an OsiBranchingInformation object
  OsiBranchingInformation usefulInformation() const;
  /** Warm start object produced by heuristic or strong branching

        If get a valid integer solution outside branch and bound then it can take
        a reasonable time to solve LP which produces clean solution.  If this object has
        any size then it will be used in solve.
    */
  inline void setBestSolutionBasis(const CoinWarmStartBasis &bestSolutionBasis)
  {
    bestSolutionBasis_ = bestSolutionBasis;
  }
  /// Redo walkback arrays
  void redoWalkBack();
  //@}

  void setMIPStart(const std::vector< std::pair< std::string, double > > &mipstart)
  {
    this->mipStart_ = mipstart;
  }

  /** if original column names will be preserved in preprocessed problem
   */
  void setKeepNamesPreproc( bool _keep ) 
  {
    this->keepNamesPreproc = _keep;
  }

  bool getKeepNamesPreproc() const 
  {
    return keepNamesPreproc;
  }

  /** may be safer to use this overload method: c++ string libraries 
     * implementation may not be binary compatible */
  void setMIPStart(int count, const char **colNames, const double colValues[]);


  const std::vector< std::pair< std::string, double > > &getMIPStart()
  {
    return this->mipStart_;
  }

  //---------------------------------------------------------------------------

private:
  ///@name Private member data
  //@{

  /// The solver associated with this model.
  OsiSolverInterface *solver_;

  /** Ownership of objects and other stuff

        0x80000000 model owns solver
        0x40000000 all variables CbcDynamicPseudoCost
    */
  unsigned int ownership_;

  /// A copy of the solver, taken at the continuous (root) node.
  OsiSolverInterface *continuousSolver_;

  /// A copy of the solver, taken at constructor or by saveReferenceSolver
  OsiSolverInterface *referenceSolver_;

  /// Message handler
  CoinMessageHandler *handler_;

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
  mutable CoinWarmStart *emptyWarmStart_;

  /// Best objective
  double bestObjective_;
  /// Best possible objective
  double bestPossibleObjective_;
  /// Sum of Changes to objective by first solve
  double sumChangeObjective1_;
  /// Sum of Changes to objective by subsequent solves
  double sumChangeObjective2_;

  /// Array holding the incumbent (best) solution.
  double *bestSolution_;
  /// Arrays holding other solutions.
  double **savedSolutions_;

  /** Array holding the current solution.

      This array is used more as a temporary.
    */
  double *currentSolution_;
  /** For testing infeasibilities - will point to
        currentSolution_ or solver-->getColSolution()
    */
  mutable const double *testSolution_;
  /** MIPstart values
      values for integer variables which will be converted to a complete integer initial feasible solution
    */
  std::vector< std::pair< std::string, double > > mipStart_;

  /** keepNamesPreproc
   *  if variables names will be preserved in the pre-processed problem
   *  (usefull in callbacks)
   **/
  bool keepNamesPreproc;

  /** Warm start object produced by heuristic or strong branching

        If get a valid integer solution outside branch and bound then it can take
        a reasonable time to solve LP which produces clean solution.  If this object has
        any size then it will be used in solve.
    */
  CoinWarmStartBasis bestSolutionBasis_;
  /// Global cuts
  CbcRowCuts globalCuts_;
  /// Global conflict cuts
  CbcRowCuts *globalConflictCuts_;

  /// Minimum degradation in objective value to continue cut generation
  double minimumDrop_;
  /// Number of solutions
  int numberSolutions_;
  /// Number of saved solutions
  int numberSavedSolutions_;
  /// Maximum number of saved solutions
  int maximumSavedSolutions_;
  /** State of search
        0 - no solution
        1 - only heuristic solutions
        2 - branched to a solution
        3 - no solution but many nodes
    */
  int stateOfSearch_;
  /// At which depths to do cuts
  int whenCuts_;
  /// Hotstart solution
  double *hotstartSolution_;
  /// Hotstart priorities
  int *hotstartPriorities_;
  /// Number of heuristic solutions
  int numberHeuristicSolutions_;
  /// Cumulative number of nodes
  int numberNodes_;
  /** Cumulative number of nodes for statistics.
        Must fix to match up
    */
  int numberNodes2_;
  /// Cumulative number of iterations
  int numberIterations_;
  /// Cumulative number of solves
  int numberSolves_;
  /// Status of problem - 0 finished, 1 stopped, 2 difficulties
  int status_;
  /** Secondary status of problem
        -1 unset (status_ will also be -1)
        0 search completed with solution
        1 linear relaxation not feasible (or worse than cutoff)
        2 stopped on gap
        3 stopped on nodes
        4 stopped on time
        5 stopped on user event
        6 stopped on solutions
     */
  int secondaryStatus_;
  /// Number of integers in problem
  int numberIntegers_;
  /// Number of rows at continuous
  int numberRowsAtContinuous_;
  /**
       -1 - cutoff as constraint not activated
       -2 - waiting to activate
       >=0 - activated
     */
  int cutoffRowNumber_;
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
  CbcNodeInfo **walkback_;
  /// preProcess used before branch and bound (optional)
  CglPreProcess *preProcess_;
  CbcNodeInfo **lastNodeInfo_;
  const OsiRowCut **lastCut_;
  int lastDepth_;
  int lastNumberCuts2_;
  int maximumCuts_;
  int *lastNumberCuts_;

  /** The list of cuts initially collected for this subproblem

      When the subproblem at this node is rebuilt, a set of cuts is collected
      for inclusion in the constraint system. If any of these cuts are
      subsequently removed because they have become loose, the corresponding
      entry is set to NULL.
    */
  CbcCountRowCut **addedCuts_;

  /** A pointer to a row cut which will be added instead of normal branching.
        After use it should be set to NULL.
    */
  OsiRowCut *nextRowCut_;

  /// Current node so can be used elsewhere
  CbcNode *currentNode_;

  /// Indices of integer variables
  int *integerVariable_;
  /// Whether of not integer
  char *integerInfo_;
  /// Holds solution at continuous (after cuts)
  double *continuousSolution_;
  /// Array marked whenever a solution is found if non-zero
  int *usedInSolution_;
  /**
        Special options
        0 bit (1) - check if cuts valid (if on debugger list)
        1 bit (2) - use current basis to check integer solution (rather than all slack)
        2 bit (4) - don't check integer solution (by solving LP)
        3 bit (8) - fast analyze
        4 bit (16) - non-linear model - so no well defined CoinPackedMatrix
        5 bit (32) - keep names
        6 bit (64) - try for dominated columns
        7 bit (128) - SOS type 1 but all declared integer
        8 bit (256) - Set to say solution just found, unset by doing cuts
        9 bit (512) - Try reduced model after 100 nodes
        10 bit (1024) - Switch on some heuristics even if seems unlikely
        11 bit (2048) - Mark as in small branch and bound
        12 bit (4096) - Funny cuts so do slow way (in some places)
        13 bit (8192) - Funny cuts so do slow way (in other places)
        14 bit (16384) - Use Cplex! for fathoming
        15 bit (32768) - Try reduced model after 0 nodes
        16 bit (65536) - Original model had integer bounds
        17 bit (131072) - Perturbation switched off
        18 bit (262144) - donor CbcModel
        19 bit (524288) - recipient CbcModel
        20 bit (1048576) - waiting for sub model to return
	22 bit (4194304) - do not initialize random seed in solver (user has)
	23 bit (8388608) - leave solver_ with cuts
	24 bit (16777216) - just get feasible if no cutoff
    */
  int specialOptions_;
  /** More special options
        at present bottom 6 bits used for shadow price mode
        1024 for experimental hotstart
        2048,4096 breaking out of cuts
        8192 slowly increase minimum drop
        16384 gomory
	32768 more heuristics in sub trees
	65536 no cuts in preprocessing
        131072 Time limits elapsed
        18 bit (262144) - Perturb fathom nodes
        19 bit (524288) - No limit on fathom nodes
        20 bit (1048576) - Reduce sum of infeasibilities before cuts
        21 bit (2097152) - Reduce sum of infeasibilities after cuts
    */
  int moreSpecialOptions_;
  /** More more special options
	0 bit (1) - find switching variables
	1 bit (2) - using fake objective until solution
	2 bit (4) - switching variables exist
	3 bit (8) - skip most of setBestSolution checks
	4 bit (16) - very lightweight preprocessing in smallB&B
	5 bit (32) - event handler needs to be cloned when parallel
	6 bit (64) - testing - use probing to make cliques
	7/8 bit (128) - try orbital branching (if nauty)
	9 bit (512) - branching on objective (later)
	10 bit (1024) - branching on constraints (later)
	11/12 bit 2048 - intermittent cuts
    */
  int moreSpecialOptions2_;
  /// User node comparison function
  CbcCompareBase *nodeCompare_;
  /// User feasibility function (see CbcFeasibleBase.hpp)
  CbcFeasibilityBase *problemFeasibility_;
  /// Tree
  CbcTree *tree_;
  /// Pointer to top of tree
  CbcFullNodeInfo *topOfTree_;
  /// A pointer to model to be used for subtrees
  CbcModel *subTreeModel_;
  /// A pointer to model from CbcHeuristic
  CbcModel *heuristicModel_;
  /// Number of times any subtree stopped on nodes, time etc
  int numberStoppedSubTrees_;
  /// Variable selection function
  CbcBranchDecision *branchingMethod_;
  /// Cut modifier function
  CbcCutModifier *cutModifier_;
  /// Strategy
  CbcStrategy *strategy_;
  /// Parent model
  CbcModel *parentModel_;
  /** Whether to automatically do presolve before branch and bound.
        0 - no
        1 - ordinary presolve
        2 - integer presolve (dodgy)
    */
  /// Pointer to array[getNumCols()] (for speed) of column lower bounds
  const double *cbcColLower_;
  /// Pointer to array[getNumCols()] (for speed) of column upper bounds
  const double *cbcColUpper_;
  /// Pointer to array[getNumRows()] (for speed) of row lower bounds
  const double *cbcRowLower_;
  /// Pointer to array[getNumRows()] (for speed) of row upper bounds
  const double *cbcRowUpper_;
  /// Pointer to array[getNumCols()] (for speed) of primal solution vector
  const double *cbcColSolution_;
  /// Pointer to array[getNumRows()] (for speed) of dual prices
  const double *cbcRowPrice_;
  /// Get a pointer to array[getNumCols()] (for speed) of reduced costs
  const double *cbcReducedCost_;
  /// Pointer to array[getNumRows()] (for speed) of row activity levels.
  const double *cbcRowActivity_;
  /// Pointer to user-defined data structure
  void *appData_;
  /// Presolve for CbcTreeLocal
  int presolve_;
  /** Maximum number of candidates to consider for strong branching.
      To disable strong branching, set this to 0.
    */
  int numberStrong_;
  /** \brief The number of branches before pseudo costs believed
           in dynamic strong branching.

      A value of 0 is  off.
    */
  int numberBeforeTrust_;
  /** \brief The number of variables for which to compute penalties
           in dynamic strong branching.
    */
  int numberPenalties_;
  /// For threads - stop after this many "iterations"
  int stopNumberIterations_;
  /** Scale factor to make penalties match strong.
        Should/will be computed */
  double penaltyScaleFactor_;
  /// Number of analyze iterations to do
  int numberAnalyzeIterations_;
  /// Arrays with analysis results
  double *analyzeResults_;
  /// Useful temporary pointer
  void *temporaryPointer_;
  /// Number of nodes infeasible by normal branching (before cuts)
  int numberInfeasibleNodes_;
  /** Problem type as set by user or found by analysis.  This will be extended
        0 - not known
        1 - Set partitioning <=
        2 - Set partitioning ==
        3 - Set covering
    */
  int problemType_;
  /// Print frequency
  int printFrequency_;
  /// Number of cut generators
  int numberCutGenerators_;
  // Cut generators
  CbcCutGenerator **generator_;
  // Cut generators before any changes
  CbcCutGenerator **virginGenerator_;
  /// Number of heuristics
  int numberHeuristics_;
  /// Heuristic solvers
  CbcHeuristic **heuristic_;
  /// Pointer to heuristic solver which found last solution (or NULL)
  CbcHeuristic *lastHeuristic_;
  /// Depth for fast nodes
  int fastNodeDepth_;
  /*! Pointer to the event handler */
#ifdef CBC_ONLY_CLP
  ClpEventHandler *eventHandler_;
#else
  CbcEventHandler *eventHandler_;
#endif
  /// Symmetry information
  CbcSymmetry *symmetryInfo_;
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
  OsiObject **object_;
  /// Now we may not own objects - just point to solver's objects
  bool ownObjects_;

  /// Original columns as created by integerPresolve or preprocessing
  int *originalColumns_;
  /// How often to scan global cuts
  int howOftenGlobalScan_;
  /** Number of times global cuts violated.  When global cut pool then this
        should be kept for each cut and type of cut */
  int numberGlobalViolations_;
  /// Number of extra iterations in fast lp
  int numberExtraIterations_;
  /// Number of extra nodes in fast lp
  int numberExtraNodes_;
  /// Number of times fast lp entered
  int numberFathoms_;
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
  /// Preferred way of branching
  int preferredWay_;
  /// Current cut pass number
  int currentPassNumber_;
  /// Maximum number of cuts (for whichGenerator_)
  int maximumWhich_;
  /// Maximum number of rows
  int maximumRows_;
  /// Random seed
  int randomSeed_;
  /// Multiple root tries
  int multipleRootTries_;
  /// Current depth
  int currentDepth_;
  /// Thread specific random number generator
  mutable CoinThreadRandom randomNumberGenerator_;
  /// Work basis for temporary use
  CoinWarmStartBasis workingBasis_;
  /// Which cut generator generated this cut
  int *whichGenerator_;
  /// Maximum number of statistics
  int maximumStatistics_;
  /// statistics
  CbcStatistics **statistics_;
  /// Maximum depth reached
  int maximumDepthActual_;
  /// Number of reduced cost fixings
  double numberDJFixed_;
  /// Probing info
  CglTreeProbingInfo *probingInfo_;
  /// Number of fixed by analyze at root
  int numberFixedAtRoot_;
  /// Number fixed by analyze so far
  int numberFixedNow_;
  /// Whether stopping on gap
  bool stoppedOnGap_;
  /// Whether event happened
  mutable bool eventHappened_;
  /// Number of long strong goes
  int numberLongStrong_;
  /// Number of old active cuts
  int numberOldActiveCuts_;
  /// Number of new cuts
  int numberNewCuts_;
  /// Strategy worked out - mainly at root node
  int searchStrategy_;
  /** Strategy for strong branching
	0 - normal
	when to do all fractional
	1 - root node
	2 - depth less than modifier
	4 - if objective == best possible
	6 - as 2+4 
	when to do all including satisfied
	10 - root node etc.
	If >=100 then do when depth <= strategy/100 (otherwise 5)
     */
  int strongStrategy_;
  /// Number of iterations in strong branching
  int numberStrongIterations_;
  /** 0 - number times strong branching done, 1 - number fixed, 2 - number infeasible
        Second group of three is a snapshot at node [6] */
  int strongInfo_[7];
  /**
        For advanced applications you may wish to modify the behavior of Cbc
        e.g. if the solver is a NLP solver then you may not have an exact
        optimum solution at each step.  This gives characteristics - just for one BAB.
        For actually saving/restoring a solution you need the actual solver one.
    */
  OsiBabSolver *solverCharacteristics_;
  /// Whether to force a resolve after takeOffCuts
  bool resolveAfterTakeOffCuts_;
  /// Maximum number of iterations (designed to be used in heuristics)
  int maximumNumberIterations_;
  /// Anything with priority >= this can be treated as continuous
  int continuousPriority_;
  /// Number of outstanding update information items
  int numberUpdateItems_;
  /// Maximum number of outstanding update information items
  int maximumNumberUpdateItems_;
  /// Update items
  CbcObjectUpdateData *updateItems_;
  /// Stored row cuts for donor/recipient CbcModel
  CglStored *storedRowCuts_;
  /**
       Parallel
       0 - off
       1 - testing
       2-99 threads
       other special meanings
    */
  int numberThreads_;
  /** thread mode
        always use numberThreads for branching
        1 set then deterministic
        2 set then use numberThreads for root cuts
        4 set then use numberThreads in root mini branch and bound
        default is 0
    */
  int threadMode_;
  /// Number of global cuts on entry to a node
  int numberGlobalCutsIn_;
  /// Thread stuff for master
  CbcBaseModel *master_;
  /// Pointer to masterthread
  CbcThread *masterThread_;
  //@}
};
/// So we can use osiObject or CbcObject during transition
void getIntegerInformation(const OsiObject *object, double &originalLower,
  double &originalUpper);
// So we can call from other programs
// Real main program
class OsiClpSolverInterface;
int CbcMain(int argc, const char *argv[], CbcModel &babSolver);
// four ways of calling
int callCbc(const char *input2, OsiClpSolverInterface &solver1);
int callCbc(const char *input2);
int callCbc(const std::string input2, OsiClpSolverInterface &solver1);
int callCbc(const std::string input2);
// When we want to load up CbcModel with options first
void CbcMain0(CbcModel &babSolver);
int CbcMain1(int argc, const char *argv[], CbcModel &babSolver);
// two ways of calling
int callCbc(const char *input2, CbcModel &babSolver);
int callCbc(const std::string input2, CbcModel &babSolver);
// And when CbcMain0 already called to initialize (with call back) (see CbcMain1 for whereFrom)
int callCbc1(const char *input2, CbcModel &babSolver, int(CbcModel *currentSolver, int whereFrom));
int CbcMain1(int argc, const char *argv[], CbcModel &babSolver, int(CbcModel *currentSolver, int whereFrom));
// For uniform setting of cut and heuristic options
void setCutAndHeuristicOptions(CbcModel &model);
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
