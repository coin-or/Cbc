/* $Id$ */
/* Copyright (C) 2004 International Business Machines
   Corporation and others.  All Rights Reserved.*/
#ifndef CbcModelC_H
#define CbcModelC_H

/* include all defines and ugly stuff */
#include "Coin_C_defines.h"

/** This is a first "C" interface to Cbc.
    It is mostly similar to the "C" interface to Clp and
    was contributed by Bob Entriken.
*/

#ifdef __cplusplus
  extern "C"{
#endif
  
  /**@name Constructors and destructor 
     These do not have an exact analogue in C++.
     The user does not need to know structure of Cbc_Model.
     
     For all functions outside this group there is an exact C++
     analogue created by taking the first parameter out, removing the Cbc_
     from name and applying the method to an object of type ClpSimplex.
  */
  /*@{*/

  /** Version */
  COINLIBAPI double COINLINKAGE Cbc_getVersion()
  ;
  /** Default Cbc_Model constructor */
  COINLIBAPI Cbc_Model * COINLINKAGE 
  Cbc_newModel()
  ;
  /** Cbc_Model Destructor */
  COINLIBAPI void COINLINKAGE 
  Cbc_deleteModel(Cbc_Model * model)
  ;
  /*@}*/

  /**@name Load model - loads some stuff and initializes others */
  /*@{*/
  /* Loads a problem (the constraints on the
      rows are given by lower and upper bounds). If a pointer is NULL then the
      following values are the default:
      <ul>
      <li> <code>colub</code>: all columns have upper bound infinity
      <li> <code>collb</code>: all columns have lower bound 0 
      <li> <code>rowub</code>: all rows have upper bound infinity
      <li> <code>rowlb</code>: all rows have lower bound -infinity
      <li> <code>obj</code>: all variables have 0 objective coefficient
      </ul>

   Just like the other loadProblem() method except that the matrix is
   given in a standard column major ordered format (without gaps). 
*/
  COINLIBAPI void COINLINKAGE 
  Cbc_loadProblem (Cbc_Model * model,  const int numcols, const int numrows,
          const CoinBigIndex * start, const int* index,
          const double* value,
          const double* collb, const double* colub,   
          const double* obj,
          const double* rowlb, const double* rowub)
  ;
  /** Read an mps file from the given filename */
  COINLIBAPI int COINLINKAGE 
  Cbc_readMps(Cbc_Model * model,const char *filename)
  ;
  /** Write an mps file from the given filename */
  COINLIBAPI void COINLINKAGE 
  Cbc_writeMps(Cbc_Model * model,const char *filename)
  ;
  /** Integer information */
  COINLIBAPI char * COINLINKAGE 
  Cbc_integerInformation(Cbc_Model * model)
  ;
  /** Copy in integer information */
  COINLIBAPI void COINLINKAGE 
  Cbc_copyInIntegerInformation(Cbc_Model * model,const char * information)
  ;
  /** Drop integer informations */
  COINLIBAPI void COINLINKAGE 
  Cbc_deleteIntegerInformation(Cbc_Model * model)
  ;
  /** Resizes rim part of model  */
  COINLIBAPI void COINLINKAGE 
  Cbc_resize (Cbc_Model * model, int newNumberRows, int newNumberColumns)
  ;
  /** Deletes rows */
  COINLIBAPI void COINLINKAGE 
  Cbc_deleteRows(Cbc_Model * model, int number, const int * which)
  ;
  /** Add rows */
  COINLIBAPI void COINLINKAGE 
  Cbc_addRows(Cbc_Model * model, const int number, const double * rowLower, 
         const double * rowUpper,
         const int * rowStarts, const int * columns,
         const double * elements)
         ;

  /** Deletes columns */
  COINLIBAPI void COINLINKAGE 
  Cbc_deleteColumns(Cbc_Model * model, int number, const int * which)
  ;
  /** Add columns */
  COINLIBAPI void COINLINKAGE 
  Cbc_addColumns(Cbc_Model * model, int number, const double * columnLower, 
      const double * columnUpper,
      const double * objective,
      const int * columnStarts, const int * rows,
      const double * elements);
  /** Drops names - makes lengthnames 0 and names empty */
  COINLIBAPI void COINLINKAGE 
  Cbc_dropNames(Cbc_Model * model)
  ;
  /** Copies in names */
  COINLIBAPI void COINLINKAGE 
  Cbc_copyNames(Cbc_Model * model, const char * const * rowNamesIn,
     const char * const * columnNamesIn)
  ;
  
  /*@}*/
  /**@name gets and sets - you will find some synonyms at the end of this file */
  /*@{*/ 
  /** Number of rows */
  COINLIBAPI int COINLINKAGE 
  Cbc_numberRows(Cbc_Model * model)
  ;
  /** Number of columns */
  COINLIBAPI int COINLINKAGE 
  Cbc_numberColumns(Cbc_Model * model)
  ;
  /** Primal tolerance to use */
  COINLIBAPI double COINLINKAGE 
  Cbc_primalTolerance(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setPrimalTolerance(Cbc_Model * model,  double value)
  ;
  /** Dual tolerance to use */
  COINLIBAPI double COINLINKAGE 
  Cbc_dualTolerance(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setDualTolerance(Cbc_Model * model,  double value) 
  ;
  /* Integer tolerance to use */
  COINLIBAPI double COINLINKAGE 
  Cbc_integerTolerance(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setIntegerTolerance(Cbc_Model * model,  double value)
  ;
  /** Dual objective limit */
  COINLIBAPI double COINLINKAGE 
  Cbc_dualObjectiveLimit(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setDualObjectiveLimit(Cbc_Model * model, double value)
  ;
  /** Objective offset */
  COINLIBAPI double COINLINKAGE 
  Cbc_objectiveOffset(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setObjectiveOffset(Cbc_Model * model, double value)
  ;
  /** Fills in array with problem name  */
  COINLIBAPI void COINLINKAGE 
  Cbc_problemName(Cbc_Model * model, int maxNumberCharacters, char * array)
  ;
  /** Sets problem name.  Must have \0 at end.  */
  COINLIBAPI int COINLINKAGE 
  Cbc_setProblemName(Cbc_Model * model, int maxNumberCharacters, char * array)
  ;
  /** Number of iterations */
  COINLIBAPI int COINLINKAGE 
  Cbc_numberIterations(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setNumberIterations(Cbc_Model * model, int numberIterations)
  ;
  /** Maximum number of iterations */
  COINLIBAPI int COINLINKAGE 
  Cbc_maximumIterations(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setMaximumIterations(Cbc_Model * model, int value)
  ;
  /** Maximum number of nodes */
  COINLIBAPI int COINLINKAGE 
  Cbc_maxNumNode(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setMaxNumNode(Cbc_Model * model, int value)
  ;
  /* Maximum number of solutions */
  COINLIBAPI int COINLINKAGE
  Cbc_maxNumSol(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setMaxNumSol(Cbc_Model * model, int value)
  ;
  /** Maximum time in seconds (from when set called) */
  COINLIBAPI double COINLINKAGE 
  Cbc_maximumSeconds(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setMaximumSeconds(Cbc_Model * model, double value)
  ;
  /** Returns true if hit maximum iterations (or time) */
  COINLIBAPI int COINLINKAGE 
  Cbc_hitMaximumIterations(Cbc_Model * model)
  ;
  /** Status of problem:
      0 - optimal
      1 - primal infeasible
      2 - dual infeasible
      3 - stopped on iterations etc
      4 - stopped due to errors
  */
  COINLIBAPI int COINLINKAGE 
  Cbc_status(Cbc_Model * model)
  ;
  /** Set problem status */
  COINLIBAPI void COINLINKAGE 
  Cbc_setProblemStatus(Cbc_Model * model, int problemStatus)
  ;
  /** Secondary status of problem - may get extended
      0 - none
      1 - primal infeasible because dual limit reached
      2 - scaled problem optimal - unscaled has primal infeasibilities
      3 - scaled problem optimal - unscaled has dual infeasibilities
      4 - scaled problem optimal - unscaled has both dual and primal infeasibilities
  */
  COINLIBAPI int COINLINKAGE 
  Cbc_secondaryStatus(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setSecondaryStatus(Cbc_Model * model, int status)
  ;
  /** Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore */
  COINLIBAPI double COINLINKAGE 
  Cbc_optimizationDirection(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setOptimizationDirection(Cbc_Model * model, double value)
  ;
  /** Primal row solution */
  COINLIBAPI double * COINLINKAGE 
  Cbc_primalRowSolution(Cbc_Model * model)
  ;
  /** Primal column solution */
  COINLIBAPI double * COINLINKAGE 
  Cbc_primalColumnSolution(Cbc_Model * model)
  ;
  /** Dual row solution */
  COINLIBAPI double * COINLINKAGE 
  Cbc_dualRowSolution(Cbc_Model * model)
  ;
  /** Reduced costs */
  COINLIBAPI double * COINLINKAGE 
  Cbc_dualColumnSolution(Cbc_Model * model)
  ;
  /** Row lower */
  COINLIBAPI double* COINLINKAGE 
  Cbc_rowLower(Cbc_Model * model)
  ;
  /** Row upper  */
  COINLIBAPI double* COINLINKAGE 
  Cbc_rowUpper(Cbc_Model * model)
  ;
  /** Objective */
  COINLIBAPI double * COINLINKAGE 
  Cbc_objective(Cbc_Model * model)
  ;            
  /** Column Lower */
  COINLIBAPI double * COINLINKAGE 
  Cbc_columnLower(Cbc_Model * model)
  ;
  /** Column Upper */
  COINLIBAPI double * COINLINKAGE 
  Cbc_columnUpper(Cbc_Model * model)
  ;
  /** Number of elements in matrix */
  COINLIBAPI int COINLINKAGE 
  Cbc_getNumElements(Cbc_Model * model)
  ;
  /** Column starts in matrix */
  COINLIBAPI const CoinBigIndex * COINLINKAGE 
  Cbc_getVectorStarts(Cbc_Model * model)
  ; 
  /** Row indices in matrix */
  COINLIBAPI const int * COINLINKAGE 
  Cbc_getIndices(Cbc_Model * model)
  ; 
  /** Column vector lengths in matrix */
  COINLIBAPI const int * COINLINKAGE 
  Cbc_getVectorLengths(Cbc_Model * model)
  ; 
  /** Element values in matrix */
  COINLIBAPI const double * COINLINKAGE 
  Cbc_getElements(Cbc_Model * model)
  ; 
  /** Objective value */
  COINLIBAPI double COINLINKAGE 
  Cbc_objectiveValue(Cbc_Model * model)
  ;
  /** Infeasibility/unbounded ray (NULL returned if none/wrong)
      Up to user to use delete [] on these arrays.  */
  COINLIBAPI double * COINLINKAGE 
  Cbc_infeasibilityRay(Cbc_Model * model)
  ;
  COINLIBAPI double * COINLINKAGE 
  Cbc_unboundedRay(Cbc_Model * model)
  ;
  /** See if status array exists (partly for OsiClp) */
  COINLIBAPI int COINLINKAGE 
  Cbc_statusExists(Cbc_Model * model)
  ;
  /** Return address of status array (char[numberRows+numberColumns]) */
  COINLIBAPI void  COINLINKAGE 
  Cbc_getBasisStatus(Cbc_Model * model, int * cstat, int * rstat)
  ;
  /** Copy in status vector */
  COINLIBAPI void COINLINKAGE 
  Cbc_setBasisStatus(Cbc_Model * model, int * cstat, int * rstat)
  ;
  
  /** User pointer for whatever reason */
  COINLIBAPI void COINLINKAGE 
  Cbc_setUserPointer (Cbc_Model * model, void * pointer)
  ;
  COINLIBAPI void * COINLINKAGE 
  Cbc_getUserPointer (Cbc_Model * model)
  ;
  /*@}*/
  /**@name Message handling.  Call backs are handled by ONE function */
  /*@{*/
  /** Pass in Callback function.
   Message numbers up to 1000000 are Clp, Coin ones have 1000000 added */
  COINLIBAPI void COINLINKAGE 
  Cbc_registerCallBack(Cbc_Model * model, 
               cbc_callback userCallBack)
  ;
  /** Unset Callback function */
  COINLIBAPI void COINLINKAGE 
  Cbc_clearCallBack(Cbc_Model * model)
  ;
  /** Amount of print out:
      0 - none
      1 - just final
      2 - just factorizations
      3 - as 2 plus a bit more
      4 - verbose
      above that 8,16,32 etc just for selective debug
  */
  COINLIBAPI void COINLINKAGE 
  Cbc_setLogLevel(Cbc_Model * model, int value)
  ;
  COINLIBAPI int COINLINKAGE 
  Cbc_logLevel(Cbc_Model * model)
  ;
  /** length of names (0 means no names0 */
  COINLIBAPI int COINLINKAGE 
  Cbc_lengthNames(Cbc_Model * model)
  ;
  /** Fill in array (at least lengthNames+1 long) with a row name */
  COINLIBAPI void COINLINKAGE 
  Cbc_rowName(Cbc_Model * model, int iRow, char * name)
  ;
  /** Fill in array (at least lengthNames+1 long) with a column name */
  COINLIBAPI void COINLINKAGE 
  Cbc_columnName(Cbc_Model * model, int iColumn, char * name)
  ;

  /*@}*/


  /**@name Functions most useful to user */
  /*@{*/
  /** General solve algorithm which can do presolve.
      See  ClpSolve.hpp for options
   */
  COINLIBAPI int COINLINKAGE 
  Cbc_initialSolve(Cbc_Model * model)
  ;
  /* General solve algorithm which can do presolve.
     See  CbcModel.hpp for options
  */
  COINLIBAPI int COINLINKAGE 
  Cbc_branchAndBound(Cbc_Model * model)
  ;
  /** Sets or unsets scaling, 0 -off, 1 equilibrium, 2 geometric, 3, auto, 4 dynamic(later) */
  COINLIBAPI void COINLINKAGE 
  Cbc_scaling(Cbc_Model * model, int mode)
  ;
  /** Gets scalingFlag */
  COINLIBAPI int COINLINKAGE 
  Cbc_scalingFlag(Cbc_Model * model)
  ;
  /** Crash - at present just aimed at dual, returns
      -2 if dual preferred and crash basis created
      -1 if dual preferred and all slack basis preferred
       0 if basis going in was not all slack
       1 if primal preferred and all slack basis preferred
       2 if primal preferred and crash basis created.
       
       if gap between bounds <="gap" variables can be flipped

       If "pivot" is
       0 No pivoting (so will just be choice of algorithm)
       1 Simple pivoting e.g. gub
       2 Mini iterations
  */
  COINLIBAPI int COINLINKAGE 
  Cbc_crash(Cbc_Model * model, double gap,int pivot)
  ;
  /*@}*/


  /**@name most useful gets and sets */
  /*@{*/ 
  /** If problem is primal feasible */
  COINLIBAPI int COINLINKAGE 
  Cbc_primalFeasible(Cbc_Model * model)
  ;
  /** If problem is dual feasible */
  COINLIBAPI int COINLINKAGE 
  Cbc_dualFeasible(Cbc_Model * model)
  ;
  /** Dual bound */
  COINLIBAPI double COINLINKAGE 
  Cbc_dualBound(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setDualBound(Cbc_Model * model, double value)
  ;
  /** Infeasibility cost */
  COINLIBAPI double COINLINKAGE 
  Cbc_infeasibilityCost(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setInfeasibilityCost(Cbc_Model * model, double value)
  ;
  /** Perturbation:
      50  - switch on perturbation
      100 - auto perturb if takes too long (1.0e-6 largest nonzero)
      101 - we are perturbed
      102 - don't try perturbing again
      default is 100
      others are for playing
  */
  COINLIBAPI int COINLINKAGE 
  Cbc_perturbation(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setPerturbation(Cbc_Model * model, int value)
  ;
  /** Current (or last) algorithm */
  COINLIBAPI int COINLINKAGE 
  Cbc_algorithm(Cbc_Model * model)
  ; 
  /** Set algorithm */
  COINLIBAPI void COINLINKAGE 
  Cbc_setAlgorithm(Cbc_Model * model, int value)
  ;
  /** Sum of dual infeasibilities */
  COINLIBAPI double COINLINKAGE 
  Cbc_sumDualInfeasibilities(Cbc_Model * model)
  ; 
  /** Number of dual infeasibilities */
  COINLIBAPI int COINLINKAGE 
  Cbc_numberDualInfeasibilities(Cbc_Model * model)
  ; 
  /** Sum of primal infeasibilities */
  COINLIBAPI double COINLINKAGE 
  Cbc_sumPrimalInfeasibilities(Cbc_Model * model)
  ; 
  /** Number of primal infeasibilities */
  COINLIBAPI int COINLINKAGE 
  Cbc_numberPrimalInfeasibilities(Cbc_Model * model)
  ; 
  /** Save model to file, returns 0 if success.  This is designed for
      use outside algorithms so does not save iterating arrays etc.
  It does not save any messaging information. 
  Does not save scaling values.
  It does not know about all types of virtual functions.
  */
  COINLIBAPI int COINLINKAGE 
  Cbc_saveModel(Cbc_Model * model, const char * fileName)
  ;
  /** Restore model from file, returns 0 if success,
      deletes current model */
  COINLIBAPI int COINLINKAGE 
  Cbc_restoreModel(Cbc_Model * model, const char * fileName)
  ;
  
  /** Just check solution (for external use) - sets sum of
      infeasibilities etc */
  COINLIBAPI void COINLINKAGE 
  Cbc_checkSolution(Cbc_Model * model)
  ;
  /*@}*/

  /******************** End of most useful part **************/
  /**@name gets and sets - some synonyms */
  /*@{*/ 
  /** Number of rows */
  COINLIBAPI int COINLINKAGE 
  Cbc_getNumRows(Cbc_Model * model)
  ;
  /** Number of columns */
  COINLIBAPI int COINLINKAGE 
  Cbc_getNumCols(Cbc_Model * model)
  ;
  /** Number of iterations */
  COINLIBAPI int COINLINKAGE 
  Cbc_getIterationCount(Cbc_Model * model)
  ;
  /** Are there a numerical difficulties? */
  COINLIBAPI int COINLINKAGE 
  Cbc_isAbandoned(Cbc_Model * model)
  ;
  /** Is optimality proven? */
  COINLIBAPI int COINLINKAGE 
  Cbc_isProvenOptimal(Cbc_Model * model)
  ;
  /** Is primal infeasiblity proven? */
  COINLIBAPI int COINLINKAGE 
  Cbc_isProvenPrimalInfeasible(Cbc_Model * model)
  ;
  /** Is dual infeasiblity proven? */
  COINLIBAPI int COINLINKAGE 
  Cbc_isProvenDualInfeasible(Cbc_Model * model)
  ;
  /** Is the given primal objective limit reached? */
  COINLIBAPI int COINLINKAGE 
  Cbc_isPrimalObjectiveLimitReached(Cbc_Model * model) 
  ;
  /** Is the given dual objective limit reached? */
  COINLIBAPI int COINLINKAGE 
  Cbc_isDualObjectiveLimitReached(Cbc_Model * model) 
  ;
  /** Iteration limit reached? */
  COINLIBAPI int COINLINKAGE 
  Cbc_isIterationLimitReached(Cbc_Model * model)
  ;
  /** Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore */
  COINLIBAPI double COINLINKAGE 
  Cbc_getObjSense(Cbc_Model * model)
  ;
  /** Primal row solution */
  COINLIBAPI const double * COINLINKAGE 
  Cbc_getRowActivity(Cbc_Model * model)
  ;
  /** Primal column solution */
  COINLIBAPI const double * COINLINKAGE 
  Cbc_getColSolution(Cbc_Model * model)
  ;
  COINLIBAPI void COINLINKAGE 
  Cbc_setColSolution(Cbc_Model * model, const double * input)
  ;
  /** Dual row solution */
  COINLIBAPI const double * COINLINKAGE 
  Cbc_getRowPrice(Cbc_Model * model)
  ;
  /** Reduced costs */
  COINLIBAPI const double * COINLINKAGE 
  Cbc_getReducedCost(Cbc_Model * model)
  ;
  /** Row lower */
  COINLIBAPI const double* COINLINKAGE 
  Cbc_getRowLower(Cbc_Model * model)
  ;
  /** Row upper  */
  COINLIBAPI const double* COINLINKAGE 
  Cbc_getRowUpper(Cbc_Model * model)
  ;
  /** Objective */
  COINLIBAPI const double * COINLINKAGE 
  Cbc_getObjCoefficients(Cbc_Model * model)
  ; 
  /** Column Lower */
  COINLIBAPI const double * COINLINKAGE 
  Cbc_getColLower(Cbc_Model * model)
  ;
  /** Column Upper */
  COINLIBAPI const double * COINLINKAGE 
  Cbc_getColUpper(Cbc_Model * model)
  ;
  /** Objective value */
  COINLIBAPI double COINLINKAGE 
  Cbc_getObjValue(Cbc_Model * model)
  ;
  /** Print the model */
  COINLIBAPI void COINLINKAGE 
  Cbc_printModel(Cbc_Model * model, const char * argPrefix)
  ;
  /** Determine whether the variable at location i is integer restricted */
  COINLIBAPI int COINLINKAGE 
  Cbc_isInteger(Cbc_Model * model, int i)
  ;
  /** Return CPU time */
  COINLIBAPI double COINLINKAGE 
  Cbc_cpuTime(Cbc_Model * model)
  ;
  /** Number of nodes explored in B&B tree */
  COINLIBAPI int COINLINKAGE 
  Cbc_getNodeCount(Cbc_Model * model)
  ;
  /** Return a copy of this model */
  COINLIBAPI Cbc_Model * COINLINKAGE 
  Cbc_clone(Cbc_Model * model)
  ;
  /** Set this the variable to be continuous */
  COINLIBAPI Cbc_Model * COINLINKAGE 
  Cbc_setContinuous(Cbc_Model * model, int iColumn)
  ;
  /** Add SOS constraints to the model using dense matrix */
  COINLIBAPI void  COINLINKAGE 
  Cbc_addSOS_Dense(Cbc_Model * model, int numObjects, const int * len,
             const int ** which, const double * weights, const int type)
  ;
  /** Add SOS constraints to the model using row-order matrix */
  COINLIBAPI void  COINLINKAGE 
  Cbc_addSOS_Sparse(Cbc_Model * model, const int * rowStarts,
             const int * rowIndices, const double * weights, const int type)
  ;
  /** Delete all object information */
  COINLIBAPI void  COINLINKAGE 
  Cbc_deleteObjects(Cbc_Model * model)
  ;
  /** Print the solution */
  COINLIBAPI void  COINLINKAGE 
  Cbc_printSolution(Cbc_Model * model)
  ;
  /** Dual initial solve */
  COINLIBAPI int COINLINKAGE 
  Cbc_initialDualSolve(Cbc_Model * model)
  ;
  /** Primal initial solve */
  COINLIBAPI int COINLINKAGE 
  Cbc_initialPrimalSolve(Cbc_Model * model)
  ;
  /** Dual algorithm - see ClpSimplexDual.hpp for method */
  COINLIBAPI int COINLINKAGE 
  Cbc_dual(Cbc_Model * model, int ifValuesPass)
  ;
  /** Primal algorithm - see ClpSimplexPrimal.hpp for method */
  COINLIBAPI int COINLINKAGE 
  Cbc_primal(Cbc_Model * model, int ifValuesPass)
  ;
  /*@}*/
#ifdef __cplusplus
    }
#endif
#endif
