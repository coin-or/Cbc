/* $Id$ */
/*
  Copyright (C) 2004 International Business Machines Corporation and others.
  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).
*/
#ifndef CbcModelC_H
#define CbcModelC_H

/* include all defines and ugly stuff */
#include "Coin_C_defines.h"
#include <stddef.h>

/*
 * Original version contributed by Bob Entriken,
 * significantly updated by Miles Lubin.
 * 2018: several updates by Haroldo 
 */

#ifdef __cplusplus
extern "C" {
#endif

/** Current version of Cbc */
COINLIBAPI const char *COINLINKAGE Cbc_getVersion(void);

/** \name Problem creation and modification routines */
//@{

/** @brief Creates an empty problem */
COINLIBAPI Cbc_Model *COINLINKAGE
Cbc_newModel(void);

/** @brief Sets problem name.
     *
     * @param model problem object
     * @param array string with problem name
     **/
COINLIBAPI int COINLINKAGE
Cbc_setProblemName(Cbc_Model *model, const char *array);

/** @brief Creates a new column
     *
     * Creates a new column (variable)
     *
     * @param model problem object
     * @param name variable name
     * @param lb column lower bound
     * @param ub column upper bound
     * @param obj objective function coefficient
     * @param isInteger 1 if variable is integral, 0 otherwise
     * @param nz number of rows (constraints) where this column appears, can be 0 if constraints will be added later
     * @param rows index of rows where this column appears, NULL if rows will be added later
     * @param coefs coefficients that this column appears in its rows, NULL if rows will be added later
     ***/
COINLIBAPI void COINLINKAGE
Cbc_addCol(Cbc_Model *model, const char *name, double lb,
  double ub, double obj, char isInteger,
  int nz, int *rows, double *coefs);

/** @brief Adds a new row 
     *
     *  Adds a new row (linear constraint) to the problem
     *
     *  @param model problem object
     *  @param name constraint name
     *  @param nz number of variables with non-zero coefficients in this row
     *  @param cols index of variables that appear in this row
     *  @param coefs cofficients that that variables appear
     *  @param sense constraint sense: L if <=, G if >=, E if =, R if ranged and N if free
     *  @param rhs right hand size
     * */
COINLIBAPI void COINLINKAGE
Cbc_addRow(Cbc_Model *model, const char *name, int nz,
  const int *cols, const double *coefs, char sense, double rhs);

/** @brief Add SOS constraints to the model using row-order matrix */
COINLIBAPI void COINLINKAGE
Cbc_addSOS(Cbc_Model *model, int numRows, const int *rowStarts,
  const int *colIndices, const double *weights, const int type);

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

     The constraint matrix is
     given in standard compressed sparse column (without gaps).
     <ul>
     <li> <code>start[i]</code> stores the starting index of the ith column
     <li> <code>index[k]</code> stores the row index of the kth nonzero element
     <li> <code>value[k]</code> stores the coefficient of the kth nonzero element
     </ul>
    */
COINLIBAPI void COINLINKAGE
Cbc_loadProblem(Cbc_Model *model, const int numcols, const int numrows,
  const CoinBigIndex *start, const int *index,
  const double *value,
  const double *collb, const double *colub,
  const double *obj,
  const double *rowlb, const double *rowub);

/** @brief Set the name of a column 
     *
     * @param model problem object 
     * @param iColumn column index
     * @param column name
     **/
COINLIBAPI void COINLINKAGE
Cbc_setColName(Cbc_Model *model, int iColumn, const char *name);

/** @brief Set the name of a row 
     *
     * @param model problem object 
     * @param iRow row index
     * @param name row name
     **/
COINLIBAPI void COINLINKAGE
Cbc_setRowName(Cbc_Model *model, int iRow, const char *name);

/** @brief Sets optimization direction
    *
    * @param model problem object 
    * @param sense: direction of optimization (1 - minimize, -1 - maximize, 0 - ignore)
    **/
COINLIBAPI void COINLINKAGE
Cbc_setObjSense(Cbc_Model *model, double sense);

/** @brief Set the lower bound of a single constraint 
     *
     * @param model problem object 
     * @param index row index
     * @param value new row lower bound
     **/
COINLIBAPI void COINLINKAGE
Cbc_setRowLower(Cbc_Model *model, int index, double value);

/** @brief  Set the upper bound of a single constraint 
     *
     * @param model problem object 
     * @param index row index
     * @param value new row upper bound
     **/
COINLIBAPI void COINLINKAGE
Cbc_setRowUpper(Cbc_Model *model, int index, double value);

/** @brief Set the objective coefficient of a single variable 
     *
     * @param model problem object 
     * @param index variable index
     * @param value new objective function coefficient for this variable
     **/
COINLIBAPI void COINLINKAGE
Cbc_setObjCoeff(Cbc_Model *model, int index, double value);

/** @brief Set the lower bound of a single variable 
     *
     * @param model problem object 
     * @param index variable index
     * @param value variable lower bound
     **/
COINLIBAPI void COINLINKAGE
Cbc_setColLower(Cbc_Model *model, int index, double value);

/** @brief Set the upper bound of a single variable 
     *
     * @param model problem object 
     * @param index variable index
     * @param value new variable upper bound
     **/
COINLIBAPI void COINLINKAGE
Cbc_setColUpper(Cbc_Model *model, int index, double value);

/** @brief Set this variable to be continuous 
     *
     * @param model problem object 
     * @param iColumn column index
     **/
COINLIBAPI void COINLINKAGE
Cbc_setContinuous(Cbc_Model *model, int iColumn);

/** @brief Set this variable to be integer 
     *
     * @param model problem object 
     * @param iColumn column index
     **/
COINLIBAPI void COINLINKAGE
Cbc_setInteger(Cbc_Model *model, int iColumn);

/** @brief Cbc_Model destructor */
COINLIBAPI void COINLINKAGE
Cbc_deleteModel(Cbc_Model *model);

/** @brief Enter initial feasible solution 
     *
     * Enter an initial feasible solution. Only the non-zero main 
     * binary/integer decision variables need to be informed. 
     * Auxiliary and/or continuous variables are computed 
     * automatically.
     * 
     * @param model problem object 
     * @param count number of variables
     * @param colNames names of variables
     * @param colValues variable values
     *
     **/
COINLIBAPI void COINLINKAGE
Cbc_setMIPStart(Cbc_Model *model, int count, const char **colNames, const double colValues[]);

/** @brief Enter initial feasible solution 
     *
     * Enter an initial feasible solution. Only the non-zero main 
     * binary/integer decision variables need to be informed. 
     * Auxiliary and/or continuous variables are computed 
     * automatically. Same as setMIPStart but using variable indexes.
     * 
     * @param model problem object 
     * @param count number of variables
     * @param colIdxs indexes of variables
     * @param colValues variable values
     *
     **/
COINLIBAPI void COINLINKAGE
Cbc_setMIPStartI(Cbc_Model *model, int count, const int colIdxs[], const double colValues[]);

/** @brief Creates a copy of the current model 
     *
     * @param model problem object 
     * @return model copy
     **/
COINLIBAPI Cbc_Model *COINLINKAGE
Cbc_clone(Cbc_Model *model);

//@}

/** \name Routines to query problem contents
*/
//@{

/** @brief Queries problem name 
     *
     * @param model problem object
     * @param maxNumberCharacters space in string array
     * @param array string where problem name will be saved
     **/
COINLIBAPI void COINLINKAGE
Cbc_problemName(Cbc_Model *model, int maxNumberCharacters, char *array);

/** @brief Number of nonzero elements in constraint matrix 
     *
     * @param model problem object
     * @return number of non-zero entries in constraint matrix
     **/
COINLIBAPI int COINLINKAGE
Cbc_getNumElements(Cbc_Model *model);

/** @brief Number of variables in the model 
     * @param model problem object
     * @return number of columns (variables)
     **/
COINLIBAPI int COINLINKAGE
Cbc_getNumCols(Cbc_Model *model);

/** @brief Number of integer variables in the model 
     *
     * @param model problem object
     * @return number of integer variables in this model
     **/
COINLIBAPI int COINLINKAGE
Cbc_getNumIntegers(Cbc_Model *model);

/** Number of constraints in the model 
     * @param model problem object
     * @return number of rows (constraints) in the model
     **/
COINLIBAPI int COINLINKAGE
Cbc_getNumRows(Cbc_Model *model);

/** @brief Queries row name 
     *
     * @param model problem object 
     * @param row index
     * @param name string where row name will be stored
     * @param string where row name will be stored
     **/
COINLIBAPI void COINLINKAGE
Cbc_getRowName(Cbc_Model *model, int iRow, char *name, size_t maxLength);

/** Queries column name
     *
     * @param model problem object 
     * @param iColumn column index
     * @param name where name will be stored
     * @param maxLength maximum length of name string
     **/
COINLIBAPI void COINLINKAGE
Cbc_getColName(Cbc_Model *model, int iColumn, char *name, size_t maxLength);

/** @brief Number of non-zero entries in a row 
     *
     * @param model problem object 
     * @param row row index
     * @return number of non-zero entries in row
     **/
COINLIBAPI int COINLINKAGE
Cbc_getRowNz(Cbc_Model *model, int row);

/** @brief Indices of variables that appear on a row 
     *
     * @param model problem object 
     * @param row row index
     * @return vector with indexes of columns that appear on this row
     **/
COINLIBAPI const int *COINLINKAGE
Cbc_getRowIndices(Cbc_Model *model, int row);

/** @brief Coefficients of variables that appear on this row 
     *
     * @param model problem object 
     * @param row row index
     * @return coefficients of variables that appear on this row
     **/
COINLIBAPI const double *COINLINKAGE
Cbc_getRowCoeffs(Cbc_Model *model, int row);

/** @brief Number of non-zero entries in a column 
     *
     * @param model problem object 
     * @param col column index
     * @return numbef of rows that this column appears
     **/
COINLIBAPI int COINLINKAGE
Cbc_getColNz(Cbc_Model *model, int col);

/** @brief Indices of rows that a column appears 
     *
     * @param model problem object 
     * @param col column index
     * @return indices of rows that this column appears
     **/
COINLIBAPI const int *COINLINKAGE
Cbc_getColIndices(Cbc_Model *model, int col);

/** @brief Coefficients that a column appear in rows 
     *
     * @param model problem object 
     * @param col column index
     * @return coefficients of this column in rows
     **/
COINLIBAPI const double *COINLINKAGE
Cbc_getColCoeffs(Cbc_Model *model, int col);

/** @brief Right hand side of a row 
     *
     * @param model problem object 
     * @param row row index
     * @return row right hand side
     **/
COINLIBAPI double COINLINKAGE
Cbc_getRowRHS(Cbc_Model *model, int row);

/** @brief Sense a row 
     * @param model problem object 
     * @param row row index
     * @return row sense: E for =, L for <=, G for >= and R for ranged row
     **/
COINLIBAPI char COINLINKAGE
Cbc_getRowSense(Cbc_Model *model, int row);

/** @brief Direction of optimization 
     *
     * @param model problem object 
     * @return Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore) 
     **/
COINLIBAPI double COINLINKAGE
Cbc_getObjSense(Cbc_Model *model);

/** @brief Constraint lower bounds 
     *
     * @param model problem object 
     * @return vector with lower bounds of constraints
     **/
COINLIBAPI const double *COINLINKAGE
Cbc_getRowLower(Cbc_Model *model);

/** @brief Constraint upper bounds 
     *
     * @param model problem object 
     * @return constraint upper bounds
     **/
COINLIBAPI const double *COINLINKAGE
Cbc_getRowUpper(Cbc_Model *model);

/** @brief Objective vector 
     *
     * @param model problem object 
     * @return vector with coefficients of variables in the objective function
     **/
COINLIBAPI const double *COINLINKAGE
Cbc_getObjCoefficients(Cbc_Model *model);

/** @brief Variable lower bounds 
     *
     * @param model problem object 
     * @return vector with lower bounds of variables
     **/
COINLIBAPI const double *COINLINKAGE
Cbc_getColLower(Cbc_Model *model);

/** @brief Variable upper bounds 
     *
     * @param model problem object 
     * @return vector with column upper bounds
     **/
COINLIBAPI const double *COINLINKAGE
Cbc_getColUpper(Cbc_Model *model);

/** @brief Determine whether the ith variable is integer restricted 
     * 
     * @param model problem object 
     * @param i variable index
     * @return 1 if variable is integer, 0 otherwise
     **/
COINLIBAPI int COINLINKAGE
Cbc_isInteger(Cbc_Model *model, int i);


//@}

/** \name Routines to load and save problems from disk
*/
//@{

/** @brief Read an mps file from the given filename 
    * 
    * @param model problem object
    * @param fileName file name 
    **/
COINLIBAPI int COINLINKAGE
Cbc_readMps(Cbc_Model *model, const char *filename);

/** @brief Read an lp file from the given filename 
     *
     * @param model problem object
     * @param fileName file name 
     **/
COINLIBAPI int COINLINKAGE
Cbc_readLp(Cbc_Model *model, const char *filename);

/** @brief Write an mps file from the given filename 
     *
     * @param model problem object
     * @param fileName file name 
     **/
COINLIBAPI void COINLINKAGE
Cbc_writeMps(Cbc_Model *model, const char *filename);

/** @brief Write an lp file from the given filename 
     *
     * @param model problem object
     * @param fileName file name 
     **/
COINLIBAPI void COINLINKAGE
Cbc_writeLp(Cbc_Model *model, const char *filename);

//@}

/**@name Getting and setting model data
     Note that problem access and modification methods,
       such as getColLower and setColLower,
       are *not valid* after calling Cbc_solve().
       Therefore it is not recommended to reuse a Cbc_Model
       object for multiple solves. A workaround is to call Cbc_clone()
       before solving.
     * */
/*@{*/

/** Provide an initial feasible solution to accelerate branch-and-bound 
     Note that feasibility of the solution is *not* verified.
    */
COINLIBAPI void COINLINKAGE
Cbc_setInitialSolution(Cbc_Model *model, const double *sol);
/** "Column start" vector of constraint matrix. Same format as Cbc_loadProblem() */
COINLIBAPI const CoinBigIndex *COINLINKAGE
Cbc_getVectorStarts(Cbc_Model *model);
/** "Row index" vector of constraint matrix */
COINLIBAPI const int *COINLINKAGE
Cbc_getIndices(Cbc_Model *model);
/** Coefficient vector of constraint matrix */
COINLIBAPI const double *COINLINKAGE
Cbc_getElements(Cbc_Model *model);

/** Maximum lenght of a row or column name */
COINLIBAPI size_t COINLINKAGE
Cbc_maxNameLength(Cbc_Model *model);
/** Print the model */
COINLIBAPI void COINLINKAGE
Cbc_printModel(Cbc_Model *model, const char *argPrefix);
/*@}*/

/**@name Solver parameters */
/*@{*/
/** Set parameter "name" to value "value". Note that this
     * translates directly to using "-name value" as a 
     * command-line argument to Cbc.*/
COINLIBAPI void COINLINKAGE
Cbc_setParameter(Cbc_Model *model, const char *name, const char *value);


/** returns the allowable gap
 */
COINLIBAPI double COINLINKAGE
Cbc_getAllowableGap(Cbc_Model *model);

/** sets the allowable gap
 */
COINLIBAPI void COINLINKAGE
Cbc_setAllowableGap(Cbc_Model *model, double allowedGap);

/** returns the allowable fraction gap
 */
COINLIBAPI double COINLINKAGE
Cbc_getAllowableFractionGap(Cbc_Model *model);

/** sets the allowable fraction gap
 */
COINLIBAPI void COINLINKAGE
Cbc_setAllowableFractionGap(Cbc_Model *model, double allowedFracionGap);

/** returns the allowable percentage gap
 */
COINLIBAPI double COINLINKAGE
Cbc_getAllowablePercentageGap(Cbc_Model *model);

/** sets the allowable percentage gap
 */
COINLIBAPI void COINLINKAGE
Cbc_setAllowablePercentageGap(Cbc_Model *model, double allowedPercentageGap);

/** returns the time limit for the search process
 */
COINLIBAPI double COINLINKAGE
Cbc_getMaximumSeconds(Cbc_Model *model);

/** sets the time limit for the search process
 */
COINLIBAPI void COINLINKAGE
Cbc_setMaximumSeconds(Cbc_Model *model, double maxSeconds);


/** returns the maximum number of nodes that can be explored in the search tree
 */
COINLIBAPI int COINLINKAGE
Cbc_getMaximumNodes(Cbc_Model *model);

/** sets the maximum number of nodes that can be explored in the search tree
 */
COINLIBAPI void COINLINKAGE
Cbc_setMaximumNodes(Cbc_Model *model, int maxNodes);

/** returns solution limit for the search process
 */
COINLIBAPI int COINLINKAGE
Cbc_getMaximumSolutions(Cbc_Model *model);

/** sets a solution limit as a stopping criterion 
 */
COINLIBAPI void COINLINKAGE
Cbc_setMaximumSolutions(Cbc_Model *model, int maxSolutions);

/** returns the current log leven
 */
COINLIBAPI int COINLINKAGE
Cbc_getLogLevel(Cbc_Model *model);

/** sets the log level
 */
COINLIBAPI void COINLINKAGE
Cbc_setLogLevel(Cbc_Model *model, int logLevel);


/** returns the cutoff
 */
COINLIBAPI double COINLINKAGE
Cbc_getCutoff(Cbc_Model *model);

/** sets the cutoff
 */
COINLIBAPI void COINLINKAGE
Cbc_setCutoff(Cbc_Model *model, double cutoff);



/*@}*/
/**@name Message handling.  Call backs are handled by ONE function */
/*@{*/
/** Pass in Callback function.
     Message numbers up to 1000000 are Clp, Coin ones have 1000000 added */
COINLIBAPI void COINLINKAGE
Cbc_registerCallBack(Cbc_Model *model,
  cbc_callback userCallBack);

/** Unset Callback function */
COINLIBAPI void COINLINKAGE
Cbc_clearCallBack(Cbc_Model *model);

COINLIBAPI void COINLINKAGE Cbc_addCutCallback( 
    Cbc_Model *model, cbc_cut_callback cutcb, 
    const char *name, void *appData );

/*@}*/

/**@name Solving the model */
/*@{*/
/* Solve the model with Cbc (using CbcMain1).
    */
COINLIBAPI int COINLINKAGE
Cbc_solve(Cbc_Model *model);
/*@}*/

/**@name Accessing the solution and optimization status */
/*@{*/

/** @brief Best feasible solution vector 
     *
     * @param model problem object
     * @return vector with best solution found
     **/
COINLIBAPI const double *COINLINKAGE
Cbc_getColSolution(Cbc_Model *model);


/** @brief Best known bound on the optimal objective value 
     *
     * @param model problem object
     * @return best possible cost (lower bound)
     **/
COINLIBAPI double COINLINKAGE
Cbc_getBestPossibleObjValue(Cbc_Model *model);

/** @brief Best integer feasible solution 
     *
     * Best integer feasible solution or NULL if no integer feas sol found 
     *
     * @param model problem object
     * @return vector with the best solution found or NULL if no feasible solution was found
     **/
COINLIBAPI double *COINLINKAGE
Cbc_bestSolution(Cbc_Model *model);

/** @brief number of integer feasible solution saved
     *
     * @param model problem object 
     * @return number of saved solutions
     **/
COINLIBAPI int COINLINKAGE
Cbc_numberSavedSolutions(Cbc_Model *model);

/** @brief Vector with the i-th saved solution
     * 
     * @param model problem object 
     * @param whichSol index of the solution to be retrieved
     * @return vector with integer feasible solution
     **/
COINLIBAPI const double *COINLINKAGE
Cbc_savedSolution(Cbc_Model *model, int whichSol);

/** @brief Cost of the whichSol solution
     *
     * @param model problem object 
     * @param whichSol solution index
     * @return solution cost
     **/
COINLIBAPI double COINLINKAGE
Cbc_savedSolutionObj(Cbc_Model *model, int whichSol);

/** @brief Queries vector of reduced costs
     *
     * @param model problem object
     * @return reduced cost vector
     **/
COINLIBAPI const double *COINLINKAGE
Cbc_getReducedCost(Cbc_Model *model);

/** If optimization was abandoned due to numerical difficulties
     *
     * @param model problem object 
     * @return 1 if numerical difficulties interrupted the optimization, 0 otherwise
     * */
COINLIBAPI int COINLINKAGE
Cbc_isAbandoned(Cbc_Model *model);

/** @brief If the optimal solution was found 
     *
     * @param model problem object 
     * @return 1 if optimal solution was found, 0 otherwise
     **/
COINLIBAPI int COINLINKAGE
Cbc_isProvenOptimal(Cbc_Model *model);

/** @brief If infeasibility was proven
     *
     * If model is infeasible, please note that infeasibility can also be declared 
     * if cutoff is informed and no solution better than the cutoff exists.
     *
     * @param model problem object 
     * @return 1 if model is infeasible, 0 otherwise
     **/
COINLIBAPI int COINLINKAGE
Cbc_isProvenInfeasible(Cbc_Model *model);

/** @brief Is continuous model unbounded ?
    *
    * @param model problem object 
    * @return 1 if model is unbounded, 0 otherwise
    * */
COINLIBAPI int COINLINKAGE
Cbc_isContinuousUnbounded(Cbc_Model *model);

/** Objective value of best feasible solution 
     *
     * @param model problem object
     * @return cost of the best solution found
     * */
COINLIBAPI double COINLINKAGE
Cbc_getObjValue(Cbc_Model *model);

/** @brief Final optimization status
     *
     * Returns the optimization status. For more info check function
     * isProvenOptimal, isProvenInfeasible, etc. Check also secondary status.
     * Possible status are:
     *
     * -1 before branchAndBound
     * 0 finished - check isProvenOptimal or isProvenInfeasible to see if solution found (or check value of best solution)
     * 1 stopped - on maxnodes, maxsols, maxtime
     * 2 execution abandoned due to numerical dificulties
     * 5 user programmed interruption
     *
     * @param model problem object 
     * @return problem status
    */
COINLIBAPI int COINLINKAGE Cbc_status(Cbc_Model *model);

/** @brief Secondary status of problem
     *
     * Returns additional information regarding the optimization status
     *
     * -1 unset (status_ will also be -1)
     *  0 search completed with solution
     *  1 linear relaxation not feasible (or worse than cutoff)
     *  2 stopped on gap
     *  3 stopped on nodes
     *  4 stopped on time
     *  5 stopped on user event
     *  6 stopped on solutions
     *  7 linear relaxation unbounded
     *  8 stopped on iteration limit
     *
     *  @model problem object 
     *  @return optimization status
    */
COINLIBAPI int COINLINKAGE
Cbc_secondaryStatus(Cbc_Model *model);

/** Sum of primal infeasibilities */
COINLIBAPI double COINLINKAGE
Cbc_sumPrimalInfeasibilities(Cbc_Model *model);

/** Number of primal infeasibilities */
COINLIBAPI int COINLINKAGE
Cbc_numberPrimalInfeasibilities(Cbc_Model *model);

/** Just check solution (for external use) - sets sum of
        infeasibilities etc */
COINLIBAPI void COINLINKAGE
Cbc_checkSolution(Cbc_Model *model);

/** Number of iterations */
COINLIBAPI int COINLINKAGE
Cbc_getIterationCount(Cbc_Model *model);

/** Node limit reached? */
COINLIBAPI int COINLINKAGE
Cbc_isNodeLimitReached(Cbc_Model *model);
/** Time limit reached? */
COINLIBAPI int COINLINKAGE
Cbc_isSecondsLimitReached(Cbc_Model *model);
/** Solution limit reached? */
COINLIBAPI int COINLINKAGE
Cbc_isSolutionLimitReached(Cbc_Model *model);
/** Are there numerical difficulties (for initialSolve) ? */
COINLIBAPI int COINLINKAGE
Cbc_isInitialSolveAbandoned(Cbc_Model *model);
/** Is optimality proven (for initialSolve) ? */
COINLIBAPI int COINLINKAGE
Cbc_isInitialSolveProvenOptimal(Cbc_Model *model);
/** Is primal infeasiblity proven (for initialSolve) ? */
COINLIBAPI int COINLINKAGE
Cbc_isInitialSolveProvenPrimalInfeasible(Cbc_Model *model);
/** "row" solution
     *  This is the vector A*x, where A is the constraint matrix
     *  and x is the current solution. */
COINLIBAPI const double *COINLINKAGE
Cbc_getRowActivity(Cbc_Model *model);
/** Number of nodes explored in B&B tree */
COINLIBAPI int COINLINKAGE
Cbc_getNodeCount(Cbc_Model *model);
/** Print the solution */
COINLIBAPI void COINLINKAGE
Cbc_printSolution(Cbc_Model *model);

/*@}*/

/** \name OsiSolverInterface related routines (used in callbacks) */
//@{

/** @brief Returns number of cols in OsiSolverInterface object */
COINLIBAPI int COINLINKAGE
Osi_getNumCols( void *osi );

/** @brief Returns column name in OsiSolverInterface object */
COINLIBAPI void COINLINKAGE
Osi_getColName( void *osi, int i, char *name, int maxLen );

/** @brief Returns column lower bounds in OsiSolverInterface object */
COINLIBAPI const double * COINLINKAGE
Osi_getColLower( void *osi );

/** @brief Returns column upper bounds in OsiSolverInterface object */
COINLIBAPI const double * COINLINKAGE
Osi_getColUpper( void *osi );

/** @brief Returns integrality information for columns in OsiSolverInterface object */
COINLIBAPI int COINLINKAGE
Osi_isInteger( void *osi, int col );

/** @brief Returns number of rows in OsiSolverInterface object */
COINLIBAPI int COINLINKAGE
Osi_getNumRows( void *osi );

COINLIBAPI int COINLINKAGE
Osi_getRowNz(void *osi, int row);

/** @brief Indices of variables that appear on a row */
COINLIBAPI const int *COINLINKAGE
Osi_getRowIndices(void *osi, int row);

/** @brief Coefficients of variables that appear on this row 
     *
     * @param model problem object 
     * @param row row index
     * @return coefficients of variables that appear on this row
     **/
COINLIBAPI const double *COINLINKAGE
Osi_getRowCoeffs(void *osi, int row);

/** @brief Right hand side of a row 
     *
     * @param model problem object 
     * @param row row index
     * @return row right hand side
     **/
COINLIBAPI double COINLINKAGE
Osi_getRowRHS(void *osi, int row);

/** @brief Sense a row 
     * @param model problem object 
     * @param row row index
     * @return row sense: E for =, L for <=, G for >= and R for ranged row
     **/
COINLIBAPI char COINLINKAGE
Osi_getRowSense(void *osi, int row);

/** @brief Returns solution vector in OsiSolverInterface object */
COINLIBAPI const double * COINLINKAGE
Osi_getColSolution( void *osi );


/*@}*/

/** \name OsiCuts related routines (used in callbacks) */
//@{

/** adds a row cut (used in callback) */
COINLIBAPI void COINLINKAGE 
OsiCuts_addRowCut( void *osiCuts, int nz, const int *idx, const double *coef, char sense, double rhs );

/*@}*/

#ifdef __cplusplus
}
#endif
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
