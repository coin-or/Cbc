/*
  Copyright (C) 2002, 2003 International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).
*/
#ifndef ClpSimplexC_H
#define ClpSimplexC_H

#include "ClpConfig.h"
#include "CoinTypes.h"

#ifdef _MSC_VER
#define CLP_LINKAGE __stdcall
#define CLP_LINKAGE_CB __cdecl
#else
#define CLP_LINKAGE
#define CLP_LINKAGE_CB
#endif

struct Clp_Simplex_s;
typedef struct Clp_Simplex_s Clp_Simplex;

struct Clp_Solve_s;
typedef struct Clp_Solve_s Clp_Solve;

/** typedef for user call back.
 *
 * The cvec are constructed so don't need to be const
 */
typedef void(CLP_LINKAGE_CB *clp_callback)(Clp_Simplex *model, int msgno, int ndouble,
  const double *dvec, int nint, const CoinBigIndex *ivec,
  int nchar, char **cvec);

/** This is a first "C" interface to Clp.
    It has similarities to the OSL V3 interface
    and only has most common functions
*/

#ifdef __cplusplus
extern "C" {
#endif

/**@name Version info
      *
      * A Clp library has a version number of the form <major>.<minor>.<release>,
      * where each of major, minor, and release are nonnegative integers.
      * For a checkout of the Clp stable branch, release is 9999.
      * For a checkout of the Clp development branch, major, minor, and release are 9999.
      */
/*@{*/
/** Clp library version number as string. */
CLPLIB_EXPORT const char *CLP_LINKAGE Clp_Version(void);
/** Major number of Clp library version. */
CLPLIB_EXPORT int CLP_LINKAGE Clp_VersionMajor(void);
/** Minor number of Clp library version. */
CLPLIB_EXPORT int CLP_LINKAGE Clp_VersionMinor(void);
/** Release number of Clp library version. */
CLPLIB_EXPORT int CLP_LINKAGE Clp_VersionRelease(void);
/*@}*/

/**@name Constructors and destructor
        These do not have an exact analogue in C++.
        The user does not need to know structure of Clp_Simplex or Clp_Solve.

        For (almost) all Clp_* functions outside this group there is an exact C++
        analogue created by taking the first parameter out, removing the Clp_
        from name and applying the method to an object of type ClpSimplex.

        Similarly, for all ClpSolve_* functions there is an exact C++
        analogue created by taking the first parameter out, removing the ClpSolve_
        from name and applying the method to an object of type ClpSolve.
     */
/*@{*/

/** Default constructor */
CLPLIB_EXPORT Clp_Simplex *CLP_LINKAGE Clp_newModel(void);
/** Destructor */
CLPLIB_EXPORT void CLP_LINKAGE Clp_deleteModel(Clp_Simplex *model);
/** Default constructor */
CLPLIB_EXPORT Clp_Solve *CLP_LINKAGE ClpSolve_new(void);
/** Destructor */
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_delete(Clp_Solve *solve);
/*@}*/

/**@name Load model - loads some stuff and initializes others */
/*@{*/
/** Loads a problem (the constraints on the
         rows are given by lower and upper bounds). If a pointer is NULL then the
         following values are the default:
         <ul>
         <li> <code>colub</code>: all columns have upper bound infinity
         <li> <code>collb</code>: all columns have lower bound 0
         <li> <code>rowub</code>: all rows have upper bound infinity
         <li> <code>rowlb</code>: all rows have lower bound -infinity
         <li> <code>obj</code>: all variables have 0 objective coefficient
         </ul>
     */
/** Just like the other loadProblem() method except that the matrix is
     given in a standard column major ordered format (without gaps). */
CLPLIB_EXPORT void CLP_LINKAGE Clp_loadProblem(Clp_Simplex *model, const int numcols, const int numrows,
  const CoinBigIndex *start, const int *index,
  const double *value,
  const double *collb, const double *colub,
  const double *obj,
  const double *rowlb, const double *rowub);

/* read quadratic part of the objective (the matrix part) */
CLPLIB_EXPORT void CLP_LINKAGE
Clp_loadQuadraticObjective(Clp_Simplex *model,
  const int numberColumns,
  const CoinBigIndex *start,
  const int *column,
  const double *element);
/** Read an mps file from the given filename */
CLPLIB_EXPORT int CLP_LINKAGE Clp_readMps(Clp_Simplex *model, const char *filename,
  int keepNames,
  int ignoreErrors);
/** Write an mps file to the given filename */
/** Format type is 0 = normal, 1 = extra or 2 = hex.
    Number across is 1 or 2.
    Use objSense = -1D to flip the objective function around. */
CLPLIB_EXPORT int CLP_LINKAGE Clp_writeMps(Clp_Simplex *model, const char *filename,
  int formatType,
  int numberAcross,
  double objSense);
/** Copy in integer informations */
CLPLIB_EXPORT void CLP_LINKAGE Clp_copyInIntegerInformation(Clp_Simplex *model, const char *information);
/** Drop integer informations */
CLPLIB_EXPORT void CLP_LINKAGE Clp_deleteIntegerInformation(Clp_Simplex *model);
/** Resizes rim part of model  */
CLPLIB_EXPORT void CLP_LINKAGE Clp_resize(Clp_Simplex *model, int newNumberRows, int newNumberColumns);
/** Deletes rows */
CLPLIB_EXPORT void CLP_LINKAGE Clp_deleteRows(Clp_Simplex *model, int number, const int *which);
/** Add rows */
CLPLIB_EXPORT void CLP_LINKAGE Clp_addRows(Clp_Simplex *model, int number, const double *rowLower,
  const double *rowUpper,
  const CoinBigIndex *rowStarts, const int *columns,
  const double *elements);

/** Deletes columns */
CLPLIB_EXPORT void CLP_LINKAGE Clp_deleteColumns(Clp_Simplex *model, int number, const int *which);
/** Add columns */
CLPLIB_EXPORT void CLP_LINKAGE Clp_addColumns(Clp_Simplex *model, int number, const double *columnLower,
  const double *columnUpper,
  const double *objective,
  const CoinBigIndex *columnStarts, const int *rows,
  const double *elements);
/** Change row lower bounds */
CLPLIB_EXPORT void CLP_LINKAGE Clp_chgRowLower(Clp_Simplex *model, const double *rowLower);
/** Change row upper bounds */
CLPLIB_EXPORT void CLP_LINKAGE Clp_chgRowUpper(Clp_Simplex *model, const double *rowUpper);
/** Change column lower bounds */
CLPLIB_EXPORT void CLP_LINKAGE Clp_chgColumnLower(Clp_Simplex *model, const double *columnLower);
/** Change column upper bounds */
CLPLIB_EXPORT void CLP_LINKAGE Clp_chgColumnUpper(Clp_Simplex *model, const double *columnUpper);
/** Change objective coefficients */
CLPLIB_EXPORT void CLP_LINKAGE Clp_chgObjCoefficients(Clp_Simplex *model, const double *objIn);
/** Change matrix coefficients */
CLPLIB_EXPORT void CLP_LINKAGE Clp_modifyCoefficient(Clp_Simplex *model, int row, int column, double newElement,
  int keepZero);
/** Drops names - makes lengthnames 0 and names empty */
CLPLIB_EXPORT void CLP_LINKAGE Clp_dropNames(Clp_Simplex *model);
/** Copies in names */
CLPLIB_EXPORT void CLP_LINKAGE Clp_copyNames(Clp_Simplex *model, const char *const *rowNames,
  const char *const *columnNames);

/*@}*/
/**@name gets and sets - you will find some synonyms at the end of this file */
/*@{*/
/** The underlying ClpSimplex model */
CLPLIB_EXPORT void* CLP_LINKAGE Clp_model(Clp_Simplex *model);
/** Number of rows */
CLPLIB_EXPORT int CLP_LINKAGE Clp_numberRows(Clp_Simplex *model);
/** Number of columns */
CLPLIB_EXPORT int CLP_LINKAGE Clp_numberColumns(Clp_Simplex *model);
/** Primal tolerance to use */
CLPLIB_EXPORT double CLP_LINKAGE Clp_primalTolerance(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setPrimalTolerance(Clp_Simplex *model, double value);
/** Dual tolerance to use */
CLPLIB_EXPORT double CLP_LINKAGE Clp_dualTolerance(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setDualTolerance(Clp_Simplex *model, double value);
/** Dual objective limit */
CLPLIB_EXPORT double CLP_LINKAGE Clp_dualObjectiveLimit(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setDualObjectiveLimit(Clp_Simplex *model, double value);
/** Objective offset */
CLPLIB_EXPORT double CLP_LINKAGE Clp_objectiveOffset(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setObjectiveOffset(Clp_Simplex *model, double value);
/** Fills in array with problem name  */
CLPLIB_EXPORT void CLP_LINKAGE Clp_problemName(Clp_Simplex *model, int maxNumberCharacters, char *array);
/* Sets problem name.  Must have \0 at end.  */
CLPLIB_EXPORT int CLP_LINKAGE
Clp_setProblemName(Clp_Simplex *model, int maxNumberCharacters, char *array);
/** Number of iterations */
CLPLIB_EXPORT int CLP_LINKAGE Clp_numberIterations(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setNumberIterations(Clp_Simplex *model, int numberIterations);
/** Maximum number of iterations */
CLPLIB_EXPORT int CLP_LINKAGE Clp_maximumIterations(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setMaximumIterations(Clp_Simplex *model, int value);
/** Maximum time in seconds (from when set called) */
CLPLIB_EXPORT double CLP_LINKAGE Clp_maximumSeconds(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setMaximumSeconds(Clp_Simplex *model, double value);
/** Returns true if hit maximum iterations (or time) */
CLPLIB_EXPORT int CLP_LINKAGE Clp_hitMaximumIterations(Clp_Simplex *model);
/** Status of problem:
         0 - optimal
         1 - primal infeasible
         2 - dual infeasible
         3 - stopped on iterations etc
         4 - stopped due to errors
     */
CLPLIB_EXPORT int CLP_LINKAGE Clp_status(Clp_Simplex *model);
/** Set problem status */
CLPLIB_EXPORT void CLP_LINKAGE Clp_setProblemStatus(Clp_Simplex *model, int problemStatus);
/** Secondary status of problem - may get extended
         0 - none
         1 - primal infeasible because dual limit reached
         2 - scaled problem optimal - unscaled has primal infeasibilities
         3 - scaled problem optimal - unscaled has dual infeasibilities
         4 - scaled problem optimal - unscaled has both dual and primal infeasibilities
     */
CLPLIB_EXPORT int CLP_LINKAGE Clp_secondaryStatus(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setSecondaryStatus(Clp_Simplex *model, int status);
/** Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore */
CLPLIB_EXPORT double CLP_LINKAGE Clp_optimizationDirection(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setOptimizationDirection(Clp_Simplex *model, double value);
/** Primal row solution */
CLPLIB_EXPORT double *CLP_LINKAGE Clp_primalRowSolution(Clp_Simplex *model);
/** Primal column solution */
CLPLIB_EXPORT double *CLP_LINKAGE Clp_primalColumnSolution(Clp_Simplex *model);
/** Dual row solution */
CLPLIB_EXPORT double *CLP_LINKAGE Clp_dualRowSolution(Clp_Simplex *model);
/** Reduced costs */
CLPLIB_EXPORT double *CLP_LINKAGE Clp_dualColumnSolution(Clp_Simplex *model);
/** Row lower */
CLPLIB_EXPORT double *CLP_LINKAGE Clp_rowLower(Clp_Simplex *model);
/** Row upper  */
CLPLIB_EXPORT double *CLP_LINKAGE Clp_rowUpper(Clp_Simplex *model);
/** Objective */
CLPLIB_EXPORT double *CLP_LINKAGE Clp_objective(Clp_Simplex *model);
/** Column Lower */
CLPLIB_EXPORT double *CLP_LINKAGE Clp_columnLower(Clp_Simplex *model);
/** Column Upper */
CLPLIB_EXPORT double *CLP_LINKAGE Clp_columnUpper(Clp_Simplex *model);
/** Number of elements in matrix */
CLPLIB_EXPORT CoinBigIndex CLP_LINKAGE Clp_getNumElements(Clp_Simplex *model);
/* Column starts in matrix */
CLPLIB_EXPORT const CoinBigIndex *CLP_LINKAGE Clp_getVectorStarts(Clp_Simplex *model);
/* Row indices in matrix */
CLPLIB_EXPORT const int *CLP_LINKAGE Clp_getIndices(Clp_Simplex *model);
/* Column vector lengths in matrix */
CLPLIB_EXPORT const int *CLP_LINKAGE Clp_getVectorLengths(Clp_Simplex *model);
/* Element values in matrix */
CLPLIB_EXPORT const double *CLP_LINKAGE Clp_getElements(Clp_Simplex *model);
/** Objective value */
CLPLIB_EXPORT double CLP_LINKAGE Clp_objectiveValue(Clp_Simplex *model);
/** Integer information */
CLPLIB_EXPORT char *CLP_LINKAGE Clp_integerInformation(Clp_Simplex *model);
/** Gives Infeasibility ray.
      *
      * Use Clp_freeRay to free the returned array.
      *
      * @return infeasibility ray, or NULL returned if none/wrong.
      */
CLPLIB_EXPORT double *CLP_LINKAGE Clp_infeasibilityRay(Clp_Simplex *model);
/** Gives ray in which the problem is unbounded.
      *
      * Use Clp_freeRay to free the returned array.
      *
      * @return unbounded ray, or NULL returned if none/wrong.
      */
CLPLIB_EXPORT double *CLP_LINKAGE Clp_unboundedRay(Clp_Simplex *model);
/** Frees a infeasibility or unbounded ray. */
CLPLIB_EXPORT void CLP_LINKAGE Clp_freeRay(Clp_Simplex *model, double *ray);
/** See if status array exists (partly for OsiClp) */
CLPLIB_EXPORT int CLP_LINKAGE Clp_statusExists(Clp_Simplex *model);
/** Return address of status array (char[numberRows+numberColumns]) */
CLPLIB_EXPORT unsigned char *CLP_LINKAGE Clp_statusArray(Clp_Simplex *model);
/** Copy in status vector */
CLPLIB_EXPORT void CLP_LINKAGE Clp_copyinStatus(Clp_Simplex *model, const unsigned char *statusArray);
/* status values are as in ClpSimplex.hpp i.e. 0 - free, 1 basic, 2 at upper,
        3 at lower, 4 superbasic, (5 fixed) */
/* Get variable basis info */
CLPLIB_EXPORT int CLP_LINKAGE Clp_getColumnStatus(Clp_Simplex *model, int sequence);
/* Get row basis info */
CLPLIB_EXPORT int CLP_LINKAGE Clp_getRowStatus(Clp_Simplex *model, int sequence);
/* Set variable basis info (and value if at bound) */
CLPLIB_EXPORT void CLP_LINKAGE Clp_setColumnStatus(Clp_Simplex *model,
  int sequence, int value);
/* Set row basis info (and value if at bound) */
CLPLIB_EXPORT void CLP_LINKAGE Clp_setRowStatus(Clp_Simplex *model,
  int sequence, int value);

/** User pointer for whatever reason */
CLPLIB_EXPORT void CLP_LINKAGE Clp_setUserPointer(Clp_Simplex *model, void *pointer);
CLPLIB_EXPORT void *CLP_LINKAGE Clp_getUserPointer(Clp_Simplex *model);
/*@}*/
/**@name Message handling.  Call backs are handled by ONE function */
/*@{*/
/** Pass in Callback function.
      Message numbers up to 1000000 are Clp, Coin ones have 1000000 added */
CLPLIB_EXPORT void CLP_LINKAGE Clp_registerCallBack(Clp_Simplex *model,
  clp_callback userCallBack);
/** Unset Callback function */
CLPLIB_EXPORT void CLP_LINKAGE Clp_clearCallBack(Clp_Simplex *model);
/** Amount of print out:
         0 - none
         1 - just final
         2 - just factorizations
         3 - as 2 plus a bit more
         4 - verbose
         above that 8,16,32 etc just for selective debug
     */
CLPLIB_EXPORT void CLP_LINKAGE Clp_setLogLevel(Clp_Simplex *model, int value);
CLPLIB_EXPORT int CLP_LINKAGE Clp_logLevel(Clp_Simplex *model);
/** length of names (0 means no names0 */
CLPLIB_EXPORT int CLP_LINKAGE Clp_lengthNames(Clp_Simplex *model);
/** Fill in array (at least lengthNames+1 long) with a row name */
CLPLIB_EXPORT void CLP_LINKAGE Clp_rowName(Clp_Simplex *model, int iRow, char *name);
/** Fill in array (at least lengthNames+1 long) with a column name */
CLPLIB_EXPORT void CLP_LINKAGE Clp_columnName(Clp_Simplex *model, int iColumn, char *name);
/** Set row name - Nice if they are short - 8 chars or less I think */
CLPLIB_EXPORT void CLP_LINKAGE Clp_setRowName(Clp_Simplex *model, int iRow, char *name);
/** Set column name - Nice if they are short - 8 chars or less I think */
CLPLIB_EXPORT void CLP_LINKAGE Clp_setColumnName(Clp_Simplex *model, int iColumn, char *name);

/*@}*/

/**@name Functions most useful to user */
/*@{*/
/** General solve algorithm which can do presolve.
         See  ClpSolve.hpp for options
      */
CLPLIB_EXPORT int CLP_LINKAGE Clp_initialSolve(Clp_Simplex *model);
/** Pass solve options. (Exception to direct analogue rule) */
CLPLIB_EXPORT int CLP_LINKAGE Clp_initialSolveWithOptions(Clp_Simplex *model, Clp_Solve *);
/** Dual initial solve */
CLPLIB_EXPORT int CLP_LINKAGE Clp_initialDualSolve(Clp_Simplex *model);
/** Primal initial solve */
CLPLIB_EXPORT int CLP_LINKAGE Clp_initialPrimalSolve(Clp_Simplex *model);
/** Barrier initial solve */
CLPLIB_EXPORT int CLP_LINKAGE Clp_initialBarrierSolve(Clp_Simplex *model);
/** Barrier initial solve, no crossover */
CLPLIB_EXPORT int CLP_LINKAGE Clp_initialBarrierNoCrossSolve(Clp_Simplex *model);
/** Dual algorithm - see ClpSimplexDual.hpp for method */
CLPLIB_EXPORT int CLP_LINKAGE Clp_dual(Clp_Simplex *model, int ifValuesPass);
/** Primal algorithm - see ClpSimplexPrimal.hpp for method */
CLPLIB_EXPORT int CLP_LINKAGE Clp_primal(Clp_Simplex *model, int ifValuesPass);
#ifndef SLIM_CLP
/** Solve the problem with the idiot code */
CLPLIB_EXPORT void CLP_LINKAGE Clp_idiot(Clp_Simplex *model, int tryhard);
#endif
/** Sets or unsets scaling, 0 -off, 1 equilibrium, 2 geometric, 3, auto, 4 dynamic(later) */
CLPLIB_EXPORT void CLP_LINKAGE Clp_scaling(Clp_Simplex *model, int mode);
/** Gets scalingFlag */
CLPLIB_EXPORT int CLP_LINKAGE Clp_scalingFlag(Clp_Simplex *model);
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
CLPLIB_EXPORT int CLP_LINKAGE Clp_crash(Clp_Simplex *model, double gap, int pivot);
/*@}*/

/**@name most useful gets and sets */
/*@{*/
/** If problem is primal feasible */
CLPLIB_EXPORT int CLP_LINKAGE Clp_primalFeasible(Clp_Simplex *model);
/** If problem is dual feasible */
CLPLIB_EXPORT int CLP_LINKAGE Clp_dualFeasible(Clp_Simplex *model);
/** Dual bound */
CLPLIB_EXPORT double CLP_LINKAGE Clp_dualBound(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setDualBound(Clp_Simplex *model, double value);
/** Infeasibility cost */
CLPLIB_EXPORT double CLP_LINKAGE Clp_infeasibilityCost(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setInfeasibilityCost(Clp_Simplex *model, double value);
/** Perturbation:
         50  - switch on perturbation
         100 - auto perturb if takes too long (1.0e-6 largest nonzero)
         101 - we are perturbed
         102 - don't try perturbing again
         default is 100
         others are for playing
     */
CLPLIB_EXPORT int CLP_LINKAGE Clp_perturbation(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setPerturbation(Clp_Simplex *model, int value);
/** Current (or last) algorithm */
CLPLIB_EXPORT int CLP_LINKAGE Clp_algorithm(Clp_Simplex *model);
/** Set algorithm */
CLPLIB_EXPORT void CLP_LINKAGE Clp_setAlgorithm(Clp_Simplex *model, int value);
/** Sum of dual infeasibilities */
CLPLIB_EXPORT double CLP_LINKAGE Clp_sumDualInfeasibilities(Clp_Simplex *model);
/** Number of dual infeasibilities */
CLPLIB_EXPORT int CLP_LINKAGE Clp_numberDualInfeasibilities(Clp_Simplex *model);
/** Sum of primal infeasibilities */
CLPLIB_EXPORT double CLP_LINKAGE Clp_sumPrimalInfeasibilities(Clp_Simplex *model);
/** Number of primal infeasibilities */
CLPLIB_EXPORT int CLP_LINKAGE Clp_numberPrimalInfeasibilities(Clp_Simplex *model);
/** Save model to file, returns 0 if success.  This is designed for
         use outside algorithms so does not save iterating arrays etc.
     It does not save any messaging information.
     Does not save scaling values.
     It does not know about all types of virtual functions.
     */
CLPLIB_EXPORT int CLP_LINKAGE Clp_saveModel(Clp_Simplex *model, const char *fileName);
/** Restore model from file, returns 0 if success,
         deletes current model */
CLPLIB_EXPORT int CLP_LINKAGE Clp_restoreModel(Clp_Simplex *model, const char *fileName);

/** Just check solution (for external use) - sets sum of
         infeasibilities etc */
CLPLIB_EXPORT void CLP_LINKAGE Clp_checkSolution(Clp_Simplex *model);
/*@}*/

/******************** End of most useful part **************/
/**@name gets and sets - some synonyms */
/*@{*/
/** Number of rows */
CLPLIB_EXPORT int CLP_LINKAGE Clp_getNumRows(Clp_Simplex *model);
/** Number of columns */
CLPLIB_EXPORT int CLP_LINKAGE Clp_getNumCols(Clp_Simplex *model);
/** Number of iterations */
CLPLIB_EXPORT int CLP_LINKAGE Clp_getIterationCount(Clp_Simplex *model);
/** Are there a numerical difficulties? */
CLPLIB_EXPORT int CLP_LINKAGE Clp_isAbandoned(Clp_Simplex *model);
/** Is optimality proven? */
CLPLIB_EXPORT int CLP_LINKAGE Clp_isProvenOptimal(Clp_Simplex *model);
/** Is primal infeasiblity proven? */
CLPLIB_EXPORT int CLP_LINKAGE Clp_isProvenPrimalInfeasible(Clp_Simplex *model);
/** Is dual infeasiblity proven? */
CLPLIB_EXPORT int CLP_LINKAGE Clp_isProvenDualInfeasible(Clp_Simplex *model);
/** Is the given primal objective limit reached? */
CLPLIB_EXPORT int CLP_LINKAGE Clp_isPrimalObjectiveLimitReached(Clp_Simplex *model);
/** Is the given dual objective limit reached? */
CLPLIB_EXPORT int CLP_LINKAGE Clp_isDualObjectiveLimitReached(Clp_Simplex *model);
/** Iteration limit reached? */
CLPLIB_EXPORT int CLP_LINKAGE Clp_isIterationLimitReached(Clp_Simplex *model);
/** Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore */
CLPLIB_EXPORT double CLP_LINKAGE Clp_getObjSense(Clp_Simplex *model);
/** Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore */
CLPLIB_EXPORT void CLP_LINKAGE Clp_setObjSense(Clp_Simplex *model, double objsen);
/** Primal row solution */
CLPLIB_EXPORT const double *CLP_LINKAGE Clp_getRowActivity(Clp_Simplex *model);
/** Primal column solution */
CLPLIB_EXPORT const double *CLP_LINKAGE Clp_getColSolution(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setColSolution(Clp_Simplex *model, const double *input);
/** Dual row solution */
CLPLIB_EXPORT const double *CLP_LINKAGE Clp_getRowPrice(Clp_Simplex *model);
/** Reduced costs */
CLPLIB_EXPORT const double *CLP_LINKAGE Clp_getReducedCost(Clp_Simplex *model);
/** Row lower */
CLPLIB_EXPORT const double *CLP_LINKAGE Clp_getRowLower(Clp_Simplex *model);
/** Row upper  */
CLPLIB_EXPORT const double *CLP_LINKAGE Clp_getRowUpper(Clp_Simplex *model);
/** Objective */
CLPLIB_EXPORT const double *CLP_LINKAGE Clp_getObjCoefficients(Clp_Simplex *model);
/** Column Lower */
CLPLIB_EXPORT const double *CLP_LINKAGE Clp_getColLower(Clp_Simplex *model);
/** Column Upper */
CLPLIB_EXPORT const double *CLP_LINKAGE Clp_getColUpper(Clp_Simplex *model);
/** Objective value */
CLPLIB_EXPORT double CLP_LINKAGE Clp_getObjValue(Clp_Simplex *model);
/** Set random seed*/
CLPLIB_EXPORT void CLP_LINKAGE Clp_setRandomSeed(Clp_Simplex *model, int seed);
/** Print model for debugging purposes */
CLPLIB_EXPORT void CLP_LINKAGE Clp_printModel(Clp_Simplex *model, const char *prefix);
/* Small element value - elements less than this set to zero,
        default is 1.0e-20 */
CLPLIB_EXPORT double CLP_LINKAGE Clp_getSmallElementValue(Clp_Simplex *model);
CLPLIB_EXPORT void CLP_LINKAGE Clp_setSmallElementValue(Clp_Simplex *model, double value);
/*@}*/

/**@name Get and set ClpSolve options
     */
/*@{*/
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setSpecialOption(Clp_Solve *, int which, int value, int extraInfo);
CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_getSpecialOption(Clp_Solve *, int which);

/** method: (see ClpSolve::SolveType)
         0 - dual simplex
         1 - primal simplex
         2 - primal or sprint
         3 - barrier
         4 - barrier no crossover
         5 - automatic
         6 - not implemented
       -- pass extraInfo == -1 for default behavior */
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setSolveType(Clp_Solve *, int method, int extraInfo);
CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_getSolveType(Clp_Solve *);

/** amount: (see ClpSolve::PresolveType)
         0 - presolve on
         1 - presolve off
         2 - presolve number
         3 - presolve number cost
       -- pass extraInfo == -1 for default behavior */
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setPresolveType(Clp_Solve *, int amount, int extraInfo);
CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_getPresolveType(Clp_Solve *);

CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_getPresolvePasses(Clp_Solve *);
CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_getExtraInfo(Clp_Solve *, int which);
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setInfeasibleReturn(Clp_Solve *, int trueFalse);
CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_infeasibleReturn(Clp_Solve *);

CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_doDual(Clp_Solve *);
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setDoDual(Clp_Solve *, int doDual);

CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_doSingleton(Clp_Solve *);
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setDoSingleton(Clp_Solve *, int doSingleton);

CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_doDoubleton(Clp_Solve *);
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setDoDoubleton(Clp_Solve *, int doDoubleton);

CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_doTripleton(Clp_Solve *);
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setDoTripleton(Clp_Solve *, int doTripleton);

CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_doTighten(Clp_Solve *);
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setDoTighten(Clp_Solve *, int doTighten);

CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_doForcing(Clp_Solve *);
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setDoForcing(Clp_Solve *, int doForcing);

CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_doImpliedFree(Clp_Solve *);
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setDoImpliedFree(Clp_Solve *, int doImpliedFree);

CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_doDupcol(Clp_Solve *);
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setDoDupcol(Clp_Solve *, int doDupcol);

CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_doDuprow(Clp_Solve *);
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setDoDuprow(Clp_Solve *, int doDuprow);

CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_doSingletonColumn(Clp_Solve *);
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setDoSingletonColumn(Clp_Solve *, int doSingleton);

CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_presolveActions(Clp_Solve *);
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setPresolveActions(Clp_Solve *, int action);

CLPLIB_EXPORT int CLP_LINKAGE ClpSolve_substitution(Clp_Solve *);
CLPLIB_EXPORT void CLP_LINKAGE ClpSolve_setSubstitution(Clp_Solve *, int value);
/*@}*/

/**@name Functions for expert users */
/*@{*/
/** gives pointer to ClpSimplex object (C++ class), return should be cast to ClpSimplex* */
CLPLIB_EXPORT void* CLP_LINKAGE Clp_getClpSimplex(Clp_Simplex *model);
/*@}*/

#ifdef __cplusplus
}
#endif
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
