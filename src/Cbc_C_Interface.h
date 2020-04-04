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
#include "CbcSolverConfig.h"

#ifdef _MSC_VER
#define CBC_LINKAGE __stdcall
#define CBC_LINKAGE_CB __cdecl
#else
#define CBC_LINKAGE
#define CBC_LINKAGE_CB
#endif

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

/*! Which method should be used to solve the linear programming problem.
 *  If a problem with integer variables, this affects only the root node.
 * */
enum LPMethod {
  LPM_Auto    = 0,  /*! Solver will decide automatically which method to use */
  LPM_Dual    = 1,  /*! Dual simplex */
  LPM_Primal  = 2,  /*! Primal simplex */
  LPM_Barrier = 3   /*! The barrier algorithm. */
};

/*! Selects the pivot selection strategy to be used
 * in the dual simplex algorithm.
 * */
enum DualPivot {
  DP_Auto       = 0,  /*! Solver will decide automatically which method to use */
  DP_Dantzig    = 1,  /*! Simple strategy, implemented as example. */
  DP_Steepest   = 2,  /*! Default strategy */
  DP_Partial    = 3,  /*! Same as steepest, but examines a subset of choices. */
  DP_PESteepest = 4  /*! Positive edge criterion, tries to avoid degenerate moves. Influenced by the psi parameter */
};

/*! Type of cutting plane */
enum CutType {
  CT_Gomory         = 0,  /*! Gomory cuts obtained from the tableau */
  CT_MIR            = 1,  /*! Mixed integer rounding cuts */
  CT_ZeroHalf       = 2,  /*! Zero-half cuts */
  CT_Clique         = 3,  /*! Clique cuts */
  CT_KnapsackCover  = 4,  /*! Knapsack cover cuts */
  CT_LiftAndProject = 5   /*! Lift and project cuts */
};

/*! Double parameters
 * */
enum DblParam {
  DBL_PARAM_PRIMAL_TOL    = 0,  /*! Tollerance to consider a solution feasible in the linear programming solver. */
  DBL_PARAM_DUAL_TOL      = 1,  /*! Tollerance for a solution to be considered optimal in the linear programming solver. */
  DBL_PARAM_ZERO_TOL      = 2,  /*! Coefficients less that this value will be ignored when reading instances */
  DBL_PARAM_INT_TOL       = 3,  /*! Maximum allowed distance from integer value for a variable to be considered integral */
  DBL_PARAM_PRESOLVE_TOL  = 4,  /*! Tollerance used in the presolver, should be increased if the pre-solver is declaring infeasible a feasible problem */
  DBL_PARAM_TIME_LIMIT    = 5,  /*! Time limit in seconds */
  DBL_PARAM_PSI           = 6,  /*! Two dimensional princing factor in the Positive Edge pivot strategy. */
  DBL_PARAM_CUTOFF        = 7,  /*! Only search for solutions with cost less-or-equal to this value. */
  DBL_PARAM_ALLOWABLE_GAP = 8,  /*! Allowable gap between the lower and upper bound to conclude the search */
  DBL_PARAM_GAP_RATIO     = 9   /*! Stops the search when the difference between the upper and lower bound is less than this fraction of the larger value */
};
#define N_DBL_PARAMS 10

/*! Integer parameters */
enum IntParam {
  INT_PARAM_PERT_VALUE          = 0,  /*! Method of perturbation, -5000 to 102, default 50 */
  INT_PARAM_IDIOT               = 1,  /*! Parameter of the "idiot" method to try to produce an initial feasible basis. -1 let the solver decide if this should be applied; 0 deactivates it and >0 sets number of passes. */
  INT_PARAM_STRONG_BRANCHING    = 2,  /*! Number of variables to be evaluated in strong branching. */
  INT_PARAM_CUT_DEPTH           = 3,  /*! Sets the application of cuts to every depth multiple of this value. -1, the default value, let the solve decide. */
  INT_PARAM_MAX_NODES           = 4,  /*! Maximum number of nodes to be explored in the search tree */
  INT_PARAM_NUMBER_BEFORE       = 5,  /*! Number of branches before trusting pseudocodes computed in strong branching. */
  INT_PARAM_FPUMP_ITS           = 6,  /*! Maximum number of iterations in the feasibility pump method. */
  INT_PARAM_MAX_SOLS            = 7,  /*! Maximum number of solutions generated during the search. Stops the search when this number of solutions is found. */
  INT_PARAM_CUT_PASS_IN_TREE    = 8,  /*! Maximum number of cuts passes in the search tree (with the exception of the root node). Default 1. */
  INT_PARAM_THREADS             = 9,  /*! Number of threads that can be used in the branch-and-bound method.*/
  INT_PARAM_CUT_PASS            = 10, /*! Number of cut passes in the root node. Default -1, solver decides */
  INT_PARAM_LOG_LEVEL           = 11, /*! Verbosity level, from 0 to 2 */
  INT_PARAM_MAX_SAVED_SOLS      = 12, /*! Size of the pool to save the best solutions found during the search. */
  INT_PARAM_MULTIPLE_ROOTS      = 13, /*! Multiple root passes to get additional cuts and solutions. */
  INT_PARAM_ROUND_INT_VARS      = 14, /*! If integer variables should be round to remove small infeasibilities. This can increase the overall amount of infeasibilities in problems with both continuous and integer variables */
  INT_PARAM_RANDOM_SEED         = 15  /*! When solving LP and MIP, randomization is used to break ties in some decisions. This changes the random seed so that multiple executions can produce different results */
};
#define N_INT_PARAMS 16
  
/** typedef for cbc callback to monitor the progress of the search
 * in terms of improved upper and lower bounds */
typedef int(CBC_LINKAGE_CB *cbc_progress_callback)(void *model,
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
typedef int (CBC_LINKAGE_CB *cbc_incumbent_callback)(void *cbcModel, double obj, int nz, char **vnames, double *x, void *appData);

typedef void(CBC_LINKAGE_CB *cbc_callback)(Cbc_Model *model, int msgno, int ndouble,
  const double *dvec, int nint, const int *ivec,
  int nchar, char **cvec);

/** 
 * \brief CBC cut callback
 *  
 * The CBC cut generation callback
 *
 * \param osiSolver an OsiSolverInterface object with the problem and the fractional solution
 *        please note that if pre-processing is on then the number of variables and constraints 
 *        in this problem will be smaller than the dimensions of the initial problem. if keepNames
 *        is true than the original variable names will be preseved to ease the mapping of the
 *        pre-processed variables and the original ones.
 *
 * \param osiCuts an OsiCuts object, where cuts should be added.
 *
 * \param appData optional pointer to an object contained some original problem informatio that
 *        can be filled by the user.
 *
 **/
typedef void(CBC_LINKAGE_CB *cbc_cut_callback)(void *osiSolver, void *osiCuts, void *appdata);

/** Current version of Cbc */
CBCSOLVERLIB_EXPORT const char *CBC_LINKAGE Cbc_getVersion(void);

/** \name Problem creation and modification routines */

/** @brief Creates an empty problem */
CBCSOLVERLIB_EXPORT Cbc_Model *CBC_LINKAGE
Cbc_newModel(void);

/** @brief Sets problem name.
   *
   * @param model problem object
   * @param array string with problem name
   **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_setProblemName(Cbc_Model *model, const char *array);

/** @brief activates/deactivates name indexes
 *
 * When name indexes are active column/row indexes can be queried fast. 
 *
 * @param model problem object
 * @param store: 1 maintain indexes of column and constraints names for searching indexes, 0 not
 **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
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
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
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
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
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
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_addRow(Cbc_Model *model, const char *name, int nz,
  const int *cols, const double *coefs, char sense, double rhs);


/** @brief Adds a lazy constraint
 *
 *  This method adds a lazy constraint, i.e. a constraint
 *  that will be included in the model only after the first
 *  integer solution violating it is generated. 
 *
 **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_addLazyConstraint(Cbc_Model *model, int nz,
  int *cols, double *coefs, char sense, double rhs);

/** @brief Deletes some rows
 *
 *  Deletes some rows from the model
 *
 *  @param model problem object
 *  @param numRows number of rows
 *  @param rows rows to be deleted
 * */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_deleteRows(Cbc_Model *model, int numRows, const int rows[]);

/** @brief Add SOS constraints to the model using row-order matrix */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
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
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
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
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setColName(Cbc_Model *model, int iColumn, const char *name);

/** @brief Set the name of a row 
  *
  * @param model problem object 
  * @param iRow row index
  * @param name row name
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setRowName(Cbc_Model *model, int iRow, const char *name);

/** @brief Sets optimization direction
  *
  * @param model problem object 
  * @param sense: direction of optimization (1 - minimize, -1 - maximize, 0 - ignore)
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setObjSense(Cbc_Model *model, double sense);

/** @brief Set the lower bound of a single constraint 
  *
  * @param model problem object 
  * @param index row index
  * @param value new row lower bound
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setRowLower(Cbc_Model *model, int index, double value);

/** @brief  Set the upper bound of a single constraint 
  *
  * @param model problem object 
  * @param index row index
  * @param value new row upper bound
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setRowUpper(Cbc_Model *model, int index, double value);

/** @brief  Sets the RHS of a constraint
  *
  * @param model problem object 
  * @param row row index
  * @param rhs value of the new RHS
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setRowRHS(Cbc_Model *model, int row, double rhs);

/** @brief Set the objective coefficient of a single variable 
  *
  * @param model problem object 
  * @param index variable index
  * @param value new objective function coefficient for this variable
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setObjCoeff(Cbc_Model *model, int index, double value);

/** @brief Set the lower bound of a single variable 
  *
  * @param model problem object 
  * @param index variable index
  * @param value variable lower bound
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setColLower(Cbc_Model *model, int index, double value);

/** @brief Set the upper bound of a single variable 
  *
  * @param model problem object 
  * @param index variable index
  * @param value new variable upper bound
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setColUpper(Cbc_Model *model, int index, double value);

/** @brief Set this variable to be continuous 
  *
  * @param model problem object 
  * @param iColumn column index
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setContinuous(Cbc_Model *model, int iColumn);

/** @brief Set this variable to be integer 
  *
  * @param model problem object 
  * @param iColumn column index
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setInteger(Cbc_Model *model, int iColumn);

/** @brief Frees memory of model object 
  *
  * @param model problem object */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
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
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
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
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setMIPStartI(Cbc_Model *model, int count, const int colIdxs[], const double colValues[]);

/** @brief Reads an initial feasible solution from a file
  *
  * Reads an initial feasible solution from a file. The file format
  * is the same used as output by CBC. In the case of a Mixed-Integer
  * Linear Program only the non-zero integer/binary variables need to 
  * be informed.
  *
  * @param model problem object 
  * @param fileName problem object 
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_readMIPStart(Cbc_Model *model, const char fileName[]);

/** @brief Creates a copy of the current model 
  *
  * @param model problem object 
  * @return model copy
  **/
CBCSOLVERLIB_EXPORT Cbc_Model *CBC_LINKAGE
Cbc_clone(Cbc_Model *model);

/** \name Routines to query problem contents
*/

/** @brief Queries problem name 
  *
  * @param model problem object
  * @param maxNumberCharacters space in string array
  * @param array string where problem name will be saved
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_problemName(Cbc_Model *model, int maxNumberCharacters, char *array);

/** @brief Number of nonzero elements in constraint matrix 
  *
  * @param model problem object
  * @return number of non-zero entries in constraint matrix
  **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_getNumElements(Cbc_Model *model);

/** @brief Number of variables in the model 
 *
  * @param model problem object
  * @return number of columns (variables)
  **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_getNumCols(Cbc_Model *model);

/** @brief Number of integer variables in the model 
  *
  * @param model problem object
  * @return number of integer variables in this model
  **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_getNumIntegers(Cbc_Model *model);

/** @brief Number of constraints in the model 
  * @param model problem object
  * @return number of rows (constraints) in the model
  **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_getNumRows(Cbc_Model *model);

/** @brief Queries row name 
  *
  * @param model problem object 
  * @param row index
  * @param name string where row name will be stored
  * @param string where row name will be stored
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_getRowName(Cbc_Model *model, int iRow, char *name, size_t maxLength);

/** @brief Queries column name
  *
  * @param model problem object 
  * @param iColumn column index
  * @param name where name will be stored
  * @param maxLength maximum length of name string
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_getColName(Cbc_Model *model, int iColumn, char *name, size_t maxLength);

/** @brief searches columns by name and returns its index
 *
 * call Cbc_storeNameIndexes to enable search by name
 *
 * @param model problem object
 * @param name column (variable) name
 * @return column index or -1 if not found
 **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_getColNameIndex(Cbc_Model *model, const char *name);

/** @brief searches rows by name and returns its index
 *
 * call Cbc_storeNameIndexes to enable search by name
 *
 * @param model problem object
 * @param name row (constraint) name
 * @return row index or -1 if not found
 **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_getRowNameIndex(Cbc_Model *model, const char *name);

/** @brief Number of non-zero entries in a row 
  *
  * @param model problem object 
  * @param row row index
  * @return number of non-zero entries in row
  **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_getRowNz(Cbc_Model *model, int row);

/** @brief Indices of variables that appear on a row 
  *
  * @param model problem object 
  * @param row row index
  * @return vector with indexes of columns that appear on this row
  **/
CBCSOLVERLIB_EXPORT const int *CBC_LINKAGE
Cbc_getRowIndices(Cbc_Model *model, int row);

/** @brief Coefficients of variables that appear on this row 
  *
  * @param model problem object 
  * @param row row index
  * @return coefficients of variables that appear on this row
  **/
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_getRowCoeffs(Cbc_Model *model, int row);

/** @brief Number of non-zero entries in a column 
  *
  * @param model problem object 
  * @param col column index
  * @return numbef of rows that this column appears
  **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_getColNz(Cbc_Model *model, int col);

/** @brief Indices of rows that a column appears 
  *
  * @param model problem object 
  * @param col column index
  * @return indices of rows that this column appears
  **/
CBCSOLVERLIB_EXPORT const int *CBC_LINKAGE
Cbc_getColIndices(Cbc_Model *model, int col);

/** @brief Coefficients that a column appear in rows 
  *
  * @param model problem object 
  * @param col column index
  * @return coefficients of this column in rows
  **/
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_getColCoeffs(Cbc_Model *model, int col);

/** @brief Right hand side of a row 
  *
  * @param model problem object 
  * @param row row index
  * @return row right hand side
  **/
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Cbc_getRowRHS(Cbc_Model *model, int row);

/** @brief Sense of a row 
  * 
  * @param model problem object 
  * @param row row index
  * @return row sense: E for =, L for <=, G for >= and R for ranged row
  **/
CBCSOLVERLIB_EXPORT char CBC_LINKAGE
Cbc_getRowSense(Cbc_Model *model, int row);

/** @brief Direction of optimization 
  *
  * @param model problem object 
  * @return Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore) 
  **/
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Cbc_getObjSense(Cbc_Model *model);

/** @brief Constraint lower bounds 
  *
  * @param model problem object 
  * @return vector with lower bounds of constraints
  **/
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_getRowLower(Cbc_Model *model);

/** @brief Constraint upper bounds 
 *
 * @param model problem object 
 * @return constraint upper bounds
 **/
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_getRowUpper(Cbc_Model *model);

/** @brief Objective vector 
  *
  * @param model problem object 
  * @return vector with coefficients of variables in the objective function
  **/
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_getObjCoefficients(Cbc_Model *model);

/** @brief Variable lower bounds 
  *
  * @param model problem object 
  * @return vector with lower bounds of variables
  **/
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_getColLower(Cbc_Model *model);

/** @brief Variable upper bounds 
  *
  * @param model problem object 
  * @return vector with column upper bounds
  **/
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_getColUpper(Cbc_Model *model);

/** @brief Determine whether the i-th variable is restricted to be integral
  * 
  * @param model problem object 
  * @param i variable index
  * @return 1 if variable is integer, 0 otherwise
  **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_isInteger(Cbc_Model *model, int i);

/** \name Routines to load and save problems from disk
*/

/** @brief Read an MPS file from the given filename 
  * 
  * @param model problem object
  * @param fileName file name 
  **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_readMps(Cbc_Model *model, const char *filename);

/** @brief Read a LP file from the given filename 
  *
  * @param model problem object
  * @param fileName file name 
  **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_readLp(Cbc_Model *model, const char *filename);

/** @brief Write an MPS file from the given filename 
  *
  * @param model problem object
  * @param fileName file name 
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_writeMps(Cbc_Model *model, const char *filename);

/** @brief Write an LP file from the given filename 
  *
  * @param model problem object
  * @param fileName file name 
  **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_writeLp(Cbc_Model *model, const char *filename);

/** @brief If Cbc was built with gzip compressed files support
  *
  * @return 1 if yes, 0 otherwise
  **/
CBCSOLVERLIB_EXPORT char CBC_LINKAGE
Cbc_supportsGzip();

/** @brief If Cbc was built with bzip2 compressed files support
  *
  * @return 1 if yes, 0 otherwise
  **/
CBCSOLVERLIB_EXPORT char CBC_LINKAGE
Cbc_supportsBzip2();

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
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setInitialSolution(Cbc_Model *model, const double *sol);

/** "Column start" vector of constraint matrix. Same format as Cbc_loadProblem() */
CBCSOLVERLIB_EXPORT const CoinBigIndex *CBC_LINKAGE
Cbc_getVectorStarts(Cbc_Model *model);
/** "Row index" vector of constraint matrix */
CBCSOLVERLIB_EXPORT const int *CBC_LINKAGE

Cbc_getIndices(Cbc_Model *model);

/** Coefficient vector of constraint matrix */
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_getElements(Cbc_Model *model);

/** Maximum lenght of a row or column name */
CBCSOLVERLIB_EXPORT size_t CBC_LINKAGE
Cbc_maxNameLength(Cbc_Model *model);

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
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setParameter(Cbc_Model *model, const char *name, const char *value);

/** Sets an integer parameter
 *
 * @param model problem object
 * @param which which integer parameter
 * @param val  value
 * 
 **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setIntParam(Cbc_Model *model, enum IntParam which, const int val);

/** Sets a double parameter
 *
 * @param model problem object
 * @param which which integer parameter
 * @param val  value
 * 
 **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setDblParam(Cbc_Model *model, enum DblParam which, const double val);





/** @brief returns the allowable gap
 *
 * @param model model object
 * @return the maximum allowable gap between the lower bound and the upper bound, when 
 *         the gap decrease to a smaller value the search is concluded
 */
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Cbc_getAllowableGap(Cbc_Model *model);

/** @brief sets the allowable gap
 *
 * @param model model object
 * @param allowedGap the maximum allowable gap between the lower bound and the upper bound, when 
 *         the gap decrease to a smaller value the search is concluded
 */

CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setAllowableGap(Cbc_Model *model, double allowedGap);

/** returns the allowable fraction gap
 */
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Cbc_getAllowableFractionGap(Cbc_Model *model);

/** sets the allowable fraction gap
 */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setAllowableFractionGap(Cbc_Model *model, double allowedFracionGap);

/** gets the tolerance for infeasibility in the LP solver
 */
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Cbc_getPrimalTolerance(Cbc_Model *model);

/** sets the tolerance for infeasibility in the LP solver
 */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setPrimalTolerance(Cbc_Model *model, double tol);

/** gets the tolerance for optimality in the LP solver
 */
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Cbc_getDualTolerance(Cbc_Model *model);

/** sets the tolerance for optimality in the LP solver
 */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setDualTolerance(Cbc_Model *model, double tol);

/** returns the time limit for the search process
 */
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Cbc_getMaximumSeconds(Cbc_Model *model);

/** sets the time limit for the search process
 */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setMaximumSeconds(Cbc_Model *model, double maxSeconds);

/** returns the maximum number of nodes that can be explored in the search tree
 */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_getMaximumNodes(Cbc_Model *model);

/** sets the maximum number of nodes that can be explored in the search tree
 */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setMaximumNodes(Cbc_Model *model, int maxNodes);

/** returns solution limit for the search process
 */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_getMaximumSolutions(Cbc_Model *model);

/** sets a solution limit as a stopping criterion 
 */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setMaximumSolutions(Cbc_Model *model, int maxSolutions);

/** returns the current log leven
 */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_getLogLevel(Cbc_Model *model);

/** sets the log level
 */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setLogLevel(Cbc_Model *model, int logLevel);

/** returns the cutoff
 */
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Cbc_getCutoff(Cbc_Model *model);

/** sets the cutoff
 */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setCutoff(Cbc_Model *model, double cutoff);

/** sets which method will be used to solve the linear programming problem
 */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setLPmethod(Cbc_Model *model, enum LPMethod lpm );

/** Returns a pointer to the OsiClpSolverInterface object 
 * containing the problem
 */
CBCSOLVERLIB_EXPORT void * CBC_LINKAGE
Cbc_getSolverPtr(Cbc_Model *model);


/** sets which pivotting method should be used in the dual simplex
 */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_setDualPivot(Cbc_Model *model, enum DualPivot dp );

/*@}*/
/**@name Message handling.  */
/*@{*/
/** Pass in Callback function.
     Message numbers up to 1000000 are Clp, Coin ones have 1000000 added */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Cbc_registerCallBack(Cbc_Model *model,
  cbc_callback userCallBack);

/** Unset Callback function */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
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
CBCSOLVERLIB_EXPORT void CBC_LINKAGE Cbc_addCutCallback( 
    Cbc_Model *model, 
    cbc_cut_callback cutcb, 
    const char *name, 
    void *appData, 
    int howOften,
    char atSolution );

/** callback to monitor new incumbent solutions **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE Cbc_addIncCallback(
    Cbc_Model *model, cbc_incumbent_callback inccb, 
    void *appData );

/** callback to monitor improvements in lower or upper
 * bounds */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE Cbc_addProgrCallback(
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
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_solve(Cbc_Model *model);

/** @brief Solves only the linear programming relaxation
  *
  * @param model problem object
  * @return execution status
  *   0  optimal 
  *   1  incomplete search (stopped on time, iterations)
  *   2  unfeasible
  *   3  unbounded
  **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_solveLinearProgram(Cbc_Model *model);


/*@}*/

/**@name Accessing the solution and optimization status */
/*@{*/

/** @brief Best feasible solution vector 
  *
  * @param model problem object
  * @return vector with best solution found
  **/
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_getColSolution(Cbc_Model *model);

/** @brief Best known bound on the optimal objective value 
  *
  * @param model problem object
  * @return best possible cost (lower bound)
  **/
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Cbc_getBestPossibleObjValue(Cbc_Model *model);

/** @brief Best integer feasible solution 
  *
  * Best integer feasible solution or NULL if no integer feas sol found 
  *
  * @param model problem object
  * @return vector with the best solution found or NULL if no feasible solution was found
  **/
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_bestSolution(Cbc_Model *model);

/** @brief number of integer feasible solution saved
  *
  * @param model problem object 
  * @return number of saved solutions
  **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_numberSavedSolutions(Cbc_Model *model);

/** @brief Vector with the i-th saved solution
  * 
  * @param model problem object 
  * @param whichSol index of the solution to be retrieved
  * @return vector with integer feasible solution
  **/
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_savedSolution(Cbc_Model *model, int whichSol);

/** @brief Cost of the whichSol solution
  *
  * @param model problem object 
  * @param whichSol solution index
  * @return solution cost
  **/
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Cbc_savedSolutionObj(Cbc_Model *model, int whichSol);

/** @brief Queries vector of reduced costs
  *
  * @param model problem object
  * @return reduced cost vector
  **/
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_getReducedCost(Cbc_Model *model);

/** @brief Queries vector of row prices (values for dual variables)
  *
  * @param model problem object
  * @return reduced cost vector
   */
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_getRowPrice(Cbc_Model *model);

/** If optimization was abandoned due to numerical difficulties
  *
  * @param model problem object 
  * @returns 1 if numerical difficulties interrupted the optimization, 0 otherwise
  * */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_isAbandoned(Cbc_Model *model);

/** @brief If the optimal solution was found 
  *
  * @param model problem object 
  * @return 1 if optimal solution was found, 0 otherwise
  **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_isProvenOptimal(Cbc_Model *model);

/** @brief If infeasibility was proven
  *
  * If model is infeasible, please note that infeasibility can also be declared 
  * if cutoff is informed and no solution better than the cutoff exists.
  *
  * @param model problem object 
  * @return 1 if model is infeasible, 0 otherwise
  **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_isProvenInfeasible(Cbc_Model *model);

/** @brief Is continuous model unbounded ?
    *
    * @param model problem object 
    * @return 1 if model is unbounded, 0 otherwise
    * */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_isContinuousUnbounded(Cbc_Model *model);

/** Objective value of best feasible solution 
  *
  * @param model problem object
  * @return cost of the best solution found
  * */
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Cbc_getObjValue(Cbc_Model *model);

/** @brief Final optimization status
  *
  * Returns the optimization status of MIP models. For more info check function
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
CBCSOLVERLIB_EXPORT int CBC_LINKAGE Cbc_status(Cbc_Model *model);

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
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_secondaryStatus(Cbc_Model *model);

/** Number of iterations */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_getIterationCount(Cbc_Model *model);

/** Node limit reached? */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_isNodeLimitReached(Cbc_Model *model);

/** Time limit reached? */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_isSecondsLimitReached(Cbc_Model *model);

/** Solution limit reached? */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_isSolutionLimitReached(Cbc_Model *model);

/** Are there numerical difficulties (for initialSolve) ? */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_isInitialSolveAbandoned(Cbc_Model *model);

/** Is optimality proven (for initialSolve) ? */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_isInitialSolveProvenOptimal(Cbc_Model *model);

/** Is primal infeasiblity proven (for initialSolve) ? */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_isInitialSolveProvenPrimalInfeasible(Cbc_Model *model);

/** "row" solution
  *  This is the vector A*x, where A is the constraint matrix
  *  and x is the current solution. */
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Cbc_getRowActivity(Cbc_Model *model);

/** Number of nodes explored in B&B tree */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Cbc_getNodeCount(Cbc_Model *model);

/*@}*/

/** \name OsiSolverInterface related routines (used in callbacks) */

/** @brief Creates a new OsiClpSolverInterface and returns a pointer to an OsiSolverInterface object */
CBCSOLVERLIB_EXPORT void * CBC_LINKAGE
Osi_newSolver();

/** @brief Solves initial LP relaxation */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Osi_initialSolve(void *osi);

/** @brief Reoptimizes linear program  */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Osi_resolve(void *osi);

/** @brief Performs branch and bound */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Osi_branchAndBound(void *osi);


/** @brief Checks if optimization was abandoned */
CBCSOLVERLIB_EXPORT char CBC_LINKAGE
Osi_isAbandoned(void *osi);

/** @brief Checks if optimal solution was found */
CBCSOLVERLIB_EXPORT char CBC_LINKAGE
Osi_isProvenOptimal(void *osi);

/** @brief Checks if problem is primal infeasible */
CBCSOLVERLIB_EXPORT char CBC_LINKAGE
Osi_isProvenPrimalInfeasible(void *osi);

/** @brief Checks if problem is dual infeasible */
CBCSOLVERLIB_EXPORT char CBC_LINKAGE
Osi_isProvenDualInfeasible(void *osi);

/** @brief Checks if primal objective limit was reached */
CBCSOLVERLIB_EXPORT char CBC_LINKAGE
Osi_isPrimalObjectiveLimitReached(void *osi);

/** @brief Checks if dual objective limit was reached */
CBCSOLVERLIB_EXPORT char CBC_LINKAGE
Osi_isDualObjectiveLimitReached(void *osi);

/** @brief Checks if iteration limit was reached */
CBCSOLVERLIB_EXPORT char CBC_LINKAGE
Osi_isIterationLimitReached(void *osi);

/** @brief Returns number of cols in OsiSolverInterface object */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Osi_getNumCols( void *osi );

/** @brief Returns column name in OsiSolverInterface object */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Osi_getColName( void *osi, int i, char *name, int maxLen );

/** @brief Returns column lower bounds in OsiSolverInterface object */
CBCSOLVERLIB_EXPORT const double * CBC_LINKAGE
Osi_getColLower( void *osi );

/** @brief Returns column upper bounds in OsiSolverInterface object */
CBCSOLVERLIB_EXPORT const double * CBC_LINKAGE
Osi_getColUpper( void *osi );

/** @brief Returns integrality information for columns in OsiSolverInterface object */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Osi_isInteger( void *osi, int col );

/** @brief Returns number of rows in OsiSolverInterface object */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Osi_getNumRows( void *osi );

/** @brief Returns number non-zeros in the constraint matrix */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Osi_getNumNz( void *osi );

/** @brief Returns number integer/binary variables */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Osi_getNumIntegers( void *osi );

CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Osi_getRowNz(void *osi, int row);

/** @brief Indices of variables that appear on a row */
CBCSOLVERLIB_EXPORT const int *CBC_LINKAGE
Osi_getRowIndices(void *osi, int row);

/** @brief Coefficients of variables that appear on this row 
     *
     * @param model problem object 
     * @param row row index
     * @return coefficients of variables that appear on this row
     **/
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Osi_getRowCoeffs(void *osi, int row);

/** @brief Right hand side of a row 
     *
     * @param model problem object 
     * @param row row index
     * @return row right hand side
     **/
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Osi_getRowRHS(void *osi, int row);

/** @brief Sense a row 
     * @param model problem object 
     * @param row row index
     * @return row sense: E for =, L for <=, G for >= and R for ranged row
     **/
CBCSOLVERLIB_EXPORT char CBC_LINKAGE
Osi_getRowSense(void *osi, int row);

/** @brief Returns vector with objective function coefficients */
CBCSOLVERLIB_EXPORT const double * CBC_LINKAGE
Osi_getObjCoefficients();

/** @brief Returns the objective sense: 1 for MIN and -1 for MAX */
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Osi_getObjSense();

/** @brief Returns solution vector in OsiSolverInterface object */
CBCSOLVERLIB_EXPORT const double * CBC_LINKAGE
Osi_getColSolution( void *osi );

/** @brief Returns vector of reduced costs */
CBCSOLVERLIB_EXPORT const double * CBC_LINKAGE
Osi_getReducedCost( void *osi );

/** @brief Returns of dual variables */
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
Osi_getRowPrice( void *osi );

/** @brief Returns the objective function value */
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Osi_getObjValue( void *osi );

/** @brief Sets column upper bound */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Osi_setColUpper (void *osi, int elementIndex, double ub);

/** @brief Sets column upper bound */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Osi_setColLower(void *osi, int elementIndex, double lb);

/** @brief Sets one objective function coefficient */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Osi_setObjCoef(void *osi, int elementIndex, double obj);

/** @brief Sets optimization direction
    *
    * @param osi OsiSolverInterface object
    * @param sense: direction of optimization (1 - minimize, -1 - maximize, 0 - ignore)
    **/
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Osi_setObjSense(void *osi, double sense);

/** @brief Sets a variable to integer */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Osi_setInteger(void *osi, int index);

/** @brief Sets a variable to continuous */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Osi_setContinuous(void *osi, int index);

/** @brief Number of non-zero entries in a column 
     *
     * @param model problem object 
     * @param col column index
     * @return numbef of rows that this column appears
     **/
CBCSOLVERLIB_EXPORT int CBC_LINKAGE
Osi_getColNz(void *model, int col);

/** @brief Indices of rows that a column appears 
     *
     * @param model problem object 
     * @param col column index
     * @return indices of rows that this column appears
     **/
CBCSOLVERLIB_EXPORT const int *CBC_LINKAGE
Osi_getColIndices(void *model, int col);

/** @brief Coefficients that a column appear in rows 
     *
     * @param model problem object 
     * @param col column index
     * @return coefficients of this column in rows
     **/
CBCSOLVERLIB_EXPORT const double *CBC_LINKAGE
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
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
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
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Osi_addRow(void *osi, const char *name, int nz,
  const int *cols, const double *coefs, char sense, double rhs);

/** @brief Returns the integer tolerance
 **/ 
CBCSOLVERLIB_EXPORT double CBC_LINKAGE
Osi_getIntegerTolerance(void *osi);

/** @brief Deletes an OsiSolverInterface object */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE
Osi_deleteSolver( void *osi );

/*@}*/

/** \name Cgl related routines */

CBCSOLVERLIB_EXPORT void CBC_LINKAGE Cgl_generateCuts( void *osiClpSolver, enum CutType ct, void *oc, int strength );

/*@}*/


/** \name OsiCuts related routines (used in callbacks) */

/** Creates a new cut pool and returns its pointer */
CBCSOLVERLIB_EXPORT void * CBC_LINKAGE 
OsiCuts_new();

/** adds a row cut (used in callback) */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE 
OsiCuts_addRowCut( void *osiCuts, int nz, const int *idx, const double *coef, char sense, double rhs );

/** adds a row cut (used in callback), stating that this is a globally valid cut */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE 
OsiCuts_addGlobalRowCut( void *osiCuts, int nz, const int *idx, const double *coef, char sense, double rhs );

/** Returns the number of row cuts stored */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE 
OsiCuts_sizeRowCuts( void *osiCuts );

/** Returns the number of row cuts stored */
CBCSOLVERLIB_EXPORT int CBC_LINKAGE 
OsiCuts_nzRowCut( void *osiCuts, int iRowCut );

/** Returns the variable indexes in a row cut */
CBCSOLVERLIB_EXPORT const int * CBC_LINKAGE 
OsiCuts_idxRowCut( void *osiCuts, int iRowCut );

/** Returns the variable coefficients in a row cut */
CBCSOLVERLIB_EXPORT const double * CBC_LINKAGE 
OsiCuts_coefRowCut( void *osiCuts, int iRowCut );

/** Returns the variable coefficients in a row cut */
CBCSOLVERLIB_EXPORT double CBC_LINKAGE 
OsiCuts_rhsRowCut( void *osiCuts, int iRowCut );

/** Returns the sense of a row cut */
CBCSOLVERLIB_EXPORT char CBC_LINKAGE 
OsiCuts_senseRowCut( void *osiCuts, int iRowCut );

/** Deletes a cut pool */
CBCSOLVERLIB_EXPORT void CBC_LINKAGE 
OsiCuts_delete( void *osiCuts );



/*@}*/

#ifdef __cplusplus
}
#endif
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
