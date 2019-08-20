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

/**
 * @file Cbc_C_Interface.h
 * @author COIN-OR CBC Development team
 * @date 15 Aug 2019
 *
 * The C API for the COIN-OR Branch-and-Cut solver
 * 
 */

/*
 * Original version contributed by Bob Entriken,
 * significantly updated by Miles Lubin.
 * 2018: several updates by Haroldo 
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Cbc_Model Cbc_Model;

/** typedef for cbc callback to monitor the progress of the search
 * in terms of improved upper and lower bounds */
typedef int(COINLINKAGE_CB *cbc_progress_callback)(void *model,
                                                   int phase,
                                                   int step,
                                                   const char *phaseName,
                                                   double seconds,
                                                   double lb,
                                                   double ub,
                                                   int nint,
                                                   int *vecint,
                                                   void *cbData
                                                   );


/** typedef for cbc callback to monitor the discovery
 * of new integer feasible solutions */
typedef int (COINLINKAGE_CB *cbc_incumbent_callback)(void *cbcModel, double obj, int nz, char **vnames, double *x, void *appData);

typedef void(COINLINKAGE_CB *cbc_callback)(Cbc_Model *model, int msgno, int ndouble,
  const double *dvec, int nint, const int *ivec,
  int nchar, char **cvec);

/** typedef for cbc cut callback osiSolver needs to be an OsiSolverInterface object,
 * osiCuts is an OsiCuts object and appdata is a pointer that will be passed to the cut
 * generation, you can use it to point to a data structure with information about the original problem,
 * for instance
 **/
typedef void(COINLINKAGE_CB *cbc_cut_callback)(void *osiSolver, void *osiCuts, void *appdata);

/** Current version of Cbc */
COINLIBAPI const char *COINLINKAGE Cbc_getVersion(void);

/** \name Problem creation and modification routines */

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

/** @brief activates/deactivates name indexes
 *
 * When name indexes are active column/row indexes can be queried fast. 
 *
 * @param model problem object
 * @param store: 1 maintain indexes of column and constraints names for searching indexes, 0 not
 **/
COINLIBAPI void COINLINKAGE
Cbc_storeNameIndexes(Cbc_Model *model, char _store);

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


/** @brief Deletes some columns
  *
  *  Deletes some columns (variables)
  *
  *  @param model problem object
  *  @param numCols number of columns that will be deleted
  *  @param cols Vector with indexes of columns that will be deleted 
  * */
COINLIBAPI void COINLINKAGE
Cbc_deleteCols(Cbc_Model *model, int numCols, const int cols[]);

/** @brief Adds a new row 
  *
  *  Adds a new row (linear constraint) to the problem
  *
  *  @param model problem object
  *  @param name constraint name
  *  @param nz number of variables with non-zero coefficients in this row
  *  @param cols index of variables that appear in this row
  *  @param coefs coefficients that that variables appear
  *  @param sense constraint sense: L if <=, G if >=, E if =, R if ranged and N if free
  *  @param rhs right hand size
  * */
COINLIBAPI void COINLINKAGE
Cbc_addRow(Cbc_Model *model, const char *name, int nz,
  const int *cols, const double *coefs, char sense, double rhs);

/** @brief Deletes some rows
 *
 *  Deletes some rows from the model
 *
 *  @param model problem object
 *  @param numRows number of rows
 *  @param rows rows to be deleted
 * */
COINLIBAPI void COINLINKAGE
Cbc_deleteRows(Cbc_Model *model, int numRows, const int rows[]);

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

     The constraint matrix is given in standard compressed sparse column
     (without gaps).
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

/** @brief Frees memory of model object 
  *
  * @param model problem object */
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

/** \name Routines to query problem contents
*/

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
 *
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

/** @brief Number of constraints in the model 
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

/** @brief Queries column name
  *
  * @param model problem object 
  * @param iColumn column index
  * @param name where name will be stored
  * @param maxLength maximum length of name string
  **/
COINLIBAPI void COINLINKAGE
Cbc_getColName(Cbc_Model *model, int iColumn, char *name, size_t maxLength);

/** @brief searches columns by name and returns its index
 *
 * call Cbc_storeNameIndexes to enable search by name
 *
 * @param model problem object
 * @param name column (variable) name
 * @return column index or -1 if not found
 **/
COINLIBAPI int COINLINKAGE
Cbc_getColNameIndex(Cbc_Model *model, const char *name);

/** @brief searches rows by name and returns its index
 *
 * call Cbc_storeNameIndexes to enable search by name
 *
 * @param model problem object
 * @param name row (constraint) name
 * @return row index or -1 if not found
 **/
COINLIBAPI int COINLINKAGE
Cbc_getRowNameIndex(Cbc_Model *model, const char *name);

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

/** @brief Sense of a row 
  * 
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

/** @brief Determine whether the i-th variable is restricted to be integral
  * 
  * @param model problem object 
  * @param i variable index
  * @return 1 if variable is integer, 0 otherwise
  **/
COINLIBAPI int COINLINKAGE
Cbc_isInteger(Cbc_Model *model, int i);

/** \name Routines to load and save problems from disk
*/

/** @brief Read an MPS file from the given filename 
  * 
  * @param model problem object
  * @param fileName file name 
  **/
COINLIBAPI int COINLINKAGE
Cbc_readMps(Cbc_Model *model, const char *filename);

/** @brief Read a LP file from the given filename 
  *
  * @param model problem object
  * @param fileName file name 
  **/
COINLIBAPI int COINLINKAGE
Cbc_readLp(Cbc_Model *model, const char *filename);

/** @brief Write an MPS file from the given filename 
  *
  * @param model problem object
  * @param fileName file name 
  **/
COINLIBAPI void COINLINKAGE
Cbc_writeMps(Cbc_Model *model, const char *filename);

/** @brief Write an LP file from the given filename 
  *
  * @param model problem object
  * @param fileName file name 
  **/
COINLIBAPI void COINLINKAGE
Cbc_writeLp(Cbc_Model *model, const char *filename);

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

/** Sets a parameter
 *
 * @param model problem object
 * @param name parameter name, e.g. cuts
 * @param name parameter value, e.g. off
 * 
 **/
COINLIBAPI void COINLINKAGE
Cbc_setParameter(Cbc_Model *model, const char *name, const char *value);

/** @brief returns the allowable gap
 *
 * @param model model object
 * @return the maximum allowable gap between the lower bound and the upper bound, when 
 *         the gap decrease to a smaller value the search is concluded
 */
COINLIBAPI double COINLINKAGE
Cbc_getAllowableGap(Cbc_Model *model);

/** @brief sets the allowable gap
 *
 * @param model model object
 * @param allowedGap the maximum allowable gap between the lower bound and the upper bound, when 
 *         the gap decrease to a smaller value the search is concluded
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

/** gets the tolerance for infeasibility in the LP solver
 */
COINLIBAPI double COINLINKAGE
Cbc_getPrimalTolerance(Cbc_Model *model);

/** sets the tolerance for infeasibility in the LP solver
 */
COINLIBAPI void COINLINKAGE
Cbc_setPrimalTolerance(Cbc_Model *model, double tol);

/** gets the tolerance for optimality in the LP solver
 */
COINLIBAPI double COINLINKAGE
Cbc_getDualTolerance(Cbc_Model *model);

/** sets the tolerance for optimality in the LP solver
 */
COINLIBAPI void COINLINKAGE
Cbc_setDualTolerance(Cbc_Model *model, double tol);

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
/**@name Message handling.  */
/*@{*/
/** Pass in Callback function.
     Message numbers up to 1000000 are Clp, Coin ones have 1000000 added */
COINLIBAPI void COINLINKAGE
Cbc_registerCallBack(Cbc_Model *model,
  cbc_callback userCallBack);

/** Unset Callback function */
COINLIBAPI void COINLINKAGE
Cbc_clearCallBack(Cbc_Model *model);

/** @brief adds a callback to generate cutting planes
 *
 * @param model mip model
 * @param cutcb cut callback function
 * @param name cut generator name
 * @param appData optional pointer to some additional information that the cut generator may need
 * @param howOften 1 if the cut generator should be called at every node, > 1 at every howOften nodes negative
 *        values have the same meaning but in this case the cut generator may be disable if not bound improvement
 *        was obtained with these cuts. -99 for cut generators that will be called only at the root node
 * @param atSolution if the cut generator must to be called also when an integer solution if found (=1) or zero otherwise
 **/
COINLIBAPI void COINLINKAGE Cbc_addCutCallback( 
    Cbc_Model *model, 
    cbc_cut_callback cutcb, 
    const char *name, 
    void *appData, 
    int howOften,
    char atSolution );

/** callback to monitor new incumbent solutions **/
COINLIBAPI void COINLINKAGE Cbc_addIncCallback(
    Cbc_Model *model, cbc_incumbent_callback inccb, 
    void *appData );

/** callback to monitor improvements in lower or upper
 * bounds */
COINLIBAPI void COINLINKAGE Cbc_addProgrCallback(
  Cbc_Model *model, cbc_progress_callback prgcbc,
  void *appData);

/*@}*/

/**@name Solving the model */
/*@{*/

/** @brief Solves the model with CBC
   *
   * @param model problem object
   * @return execution status, for MIPs: 
   *   -1 before branchAndBound 
   *   0  finished - check isProvenOptimal or isProvenInfeasible to see if solution found (or check value of best solution) 
   *   1  stopped - on maxnodes, maxsols, maxtime 
   *   2  difficulties so run was abandoned 
   *   5  event user programmed event occurred
   **/
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

/** @brief Queries vector of row prices (values for dual variables)
  *
  * @param model problem object
  * @return reduced cost vector
   */
COINLIBAPI const double *COINLINKAGE
Cbc_getRowPrice(Cbc_Model *model);

/** If optimization was abandoned due to numerical difficulties
  *
  * @param model problem object 
  * @returns 1 if numerical difficulties interrupted the optimization, 0 otherwise
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

/** @brief Creates a new OsiClpSolverInterface and returns a pointer to an OsiSolverInterface object */
COINLIBAPI void * COINLINKAGE
Osi_newSolver();

/** @brief Solves initial LP relaxation */
COINLIBAPI void COINLINKAGE
Osi_initialSolve(void *osi);

/** @brief Reoptimizes linear program  */
COINLIBAPI void COINLINKAGE
Osi_resolve(void *osi);

/** @brief Performs branch and bound */
COINLIBAPI void COINLINKAGE
Osi_branchAndBound(void *osi);


/** @brief Checks if optimization was abandoned */
COINLIBAPI char COINLINKAGE
Osi_isAbandoned(void *osi);

/** @brief Checks if optimal solution was found */
COINLIBAPI char COINLINKAGE
Osi_isProvenOptimal(void *osi);

/** @brief Checks if problem is primal infeasible */
COINLIBAPI char COINLINKAGE
Osi_isProvenPrimalInfeasible(void *osi);

/** @brief Checks if problem is dual infeasible */
COINLIBAPI char COINLINKAGE
Osi_isProvenDualInfeasible(void *osi);

/** @brief Checks if primal objective limit was reached */
COINLIBAPI char COINLINKAGE
Osi_isPrimalObjectiveLimitReached(void *osi);

/** @brief Checks if dual objective limit was reached */
COINLIBAPI char COINLINKAGE
Osi_isDualObjectiveLimitReached(void *osi);

/** @brief Checks if iteration limit was reached */
COINLIBAPI char COINLINKAGE
Osi_isIterationLimitReached(void *osi);

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

/** @brief Returns number non-zeros in the constraint matrix */
COINLIBAPI int COINLINKAGE
Osi_getNumNz( void *osi );

/** @brief Returns number integer/binary variables */
COINLIBAPI int COINLINKAGE
Osi_getNumIntegers( void *osi );

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

/** @brief Returns vector with objective function coefficients */
COINLIBAPI const double * COINLINKAGE
Osi_getObjCoefficients();

/** @brief Returns the objective sense: 1 for MIN and -1 for MAX */
COINLIBAPI double COINLINKAGE
Osi_getObjSense();

/** @brief Returns solution vector in OsiSolverInterface object */
COINLIBAPI const double * COINLINKAGE
Osi_getColSolution( void *osi );

/** @brief Returns vector of reduced costs */
COINLIBAPI const double * COINLINKAGE
Osi_getReducedCost( void *osi );

/** @brief Returns of dual variables */
COINLIBAPI const double *COINLINKAGE
Osi_getRowPrice( void *osi );

/** @brief Returns the objective function value */
COINLIBAPI double COINLINKAGE
Osi_getObjValue( void *osi );

/** @brief Sets column upper bound */
COINLIBAPI void COINLINKAGE
Osi_setColUpper (void *osi, int elementIndex, double ub);

/** @brief Sets column upper bound */
COINLIBAPI void COINLINKAGE
Osi_setColLower(void *osi, int elementIndex, double lb);

/** @brief Sets one objective function coefficient */
COINLIBAPI void COINLINKAGE
Osi_setObjCoef(void *osi, int elementIndex, double obj);

/** @brief Sets optimization direction
    *
    * @param osi OsiSolverInterface object
    * @param sense: direction of optimization (1 - minimize, -1 - maximize, 0 - ignore)
    **/
COINLIBAPI void COINLINKAGE
Osi_setObjSense(void *osi, double sense);

/** @brief Sets a variable to integer */
COINLIBAPI void COINLINKAGE
Osi_setInteger(void *osi, int index);

/** @brief Sets a variable to continuous */
COINLIBAPI void COINLINKAGE
Osi_setContinuous(void *osi, int index);

/** @brief Number of non-zero entries in a column 
     *
     * @param model problem object 
     * @param col column index
     * @return numbef of rows that this column appears
     **/
COINLIBAPI int COINLINKAGE
Osi_getColNz(void *model, int col);

/** @brief Indices of rows that a column appears 
     *
     * @param model problem object 
     * @param col column index
     * @return indices of rows that this column appears
     **/
COINLIBAPI const int *COINLINKAGE
Osi_getColIndices(void *model, int col);

/** @brief Coefficients that a column appear in rows 
     *
     * @param model problem object 
     * @param col column index
     * @return coefficients of this column in rows
     **/
COINLIBAPI const double *COINLINKAGE
Osi_getColCoeffs(void *model, int col);

/** @brief Creates a new column
     *
     * Creates a new column (variable)
     *
     * @param osi OsiSolverInterface object
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
Osi_addCol(void *osi, const char *name, double lb,
  double ub, double obj, char isInteger,
  int nz, int *rows, double *coefs);

/** @brief Adds a new row 
     *
     *  Adds a new row (linear constraint) to the problem
     *
     *  @param osi OsiSolverInterface object
     *  @param name constraint name
     *  @param nz number of variables with non-zero coefficients in this row
     *  @param cols index of variables that appear in this row
     *  @param coefs cofficients that that variables appear
     *  @param sense constraint sense: L if <=, G if >=, E if =, R if ranged and N if free
     *  @param rhs right hand size
     * */
COINLIBAPI void COINLINKAGE
Osi_addRow(void *osi, const char *name, int nz,
  const int *cols, const double *coefs, char sense, double rhs);

/** @brief Deletes an OsiSolverInterface object */
COINLIBAPI void COINLINKAGE
Osi_deleteSolver( void *osi );

/*@}*/

/** \name OsiCuts related routines (used in callbacks) */

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
