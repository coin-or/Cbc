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
 * Original verison contributed by Bob Entriken,
 * significantly updated by Miles Lubin.
*/


#ifdef __cplusplus
extern "C" {
#endif

    /**@name Constructors and destructor
      This is a "C" interface to Cbc.
      The user does not need to know structure of Cbc_Model.
    */
    /*@{*/

    /** Default Cbc_Model constructor */
    COINLIBAPI Cbc_Model * COINLINKAGE
    Cbc_newModel(void)
    ;
    /** Cbc_Model Destructor */
    COINLIBAPI void COINLINKAGE
    Cbc_deleteModel(Cbc_Model * model)
    ;
    /** Current version of Cbc */
    COINLIBAPI const char* COINLINKAGE Cbc_getVersion(void)
    ;
    /*@}*/

    /**@name Getting and setting model data
     Note that problem access and modification methods,
       such as getColLower and setColLower,
       are *not valid* after calling Cbc_solve().
       Therefore it is not recommended to reuse a Cbc_Model
       object for multiple solves. A workaround is to call Cbc_clone()
       before solving.
     * */
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

     The constraint matrix is
     given in standard compressed sparse column (without gaps).
     <ul>
     <li> <code>start[i]</code> stores the starting index of the ith column
     <li> <code>index[k]</code> stores the row index of the kth nonzero element
     <li> <code>value[k]</code> stores the coefficient of the kth nonzero element
     </ul>
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
    Cbc_readMps(Cbc_Model * model, const char *filename)
    ;
    /** Read an lp file from the given filename */
    COINLIBAPI int COINLINKAGE
    Cbc_readLp(Cbc_Model * model, const char *filename)
    ;
    /** Write an mps file from the given filename */
    COINLIBAPI void COINLINKAGE
    Cbc_writeMps(Cbc_Model * model, const char *filename)
    ;
    /** Write an lp file from the given filename */
    COINLIBAPI void COINLINKAGE
    Cbc_writeLp(Cbc_Model * model, const char *filename)
    ;
    /** Provide an initial feasible solution to accelerate branch-and-bound 
     Note that feasibility of the solution is *not* verified.
    */
    COINLIBAPI void COINLINKAGE
    Cbc_setInitialSolution(Cbc_Model *model, const double * sol)
    ;
    /** Fills in array with problem name  */
    COINLIBAPI void COINLINKAGE
    Cbc_problemName(Cbc_Model * model, int maxNumberCharacters, char * array)
    ;
    /** Sets problem name.
    
      \p array must be a null-terminated string.
    */
    COINLIBAPI int COINLINKAGE
    Cbc_setProblemName(Cbc_Model * model, const char * array)
    ;

    /** Number of nonzero elements in constraint matrix */
    COINLIBAPI int COINLINKAGE
    Cbc_getNumElements(Cbc_Model * model)
    ;
    /** "Column start" vector of constraint matrix. Same format as Cbc_loadProblem() */
    COINLIBAPI const CoinBigIndex * COINLINKAGE
    Cbc_getVectorStarts(Cbc_Model * model)
    ;
    /** "Row index" vector of constraint matrix */
    COINLIBAPI const int * COINLINKAGE
    Cbc_getIndices(Cbc_Model * model)
    ;
    /** Coefficient vector of constraint matrix */
    COINLIBAPI const double * COINLINKAGE
    Cbc_getElements(Cbc_Model * model)
    ;

    /** Maximum lenght of a row or column name */
    COINLIBAPI size_t COINLINKAGE
    Cbc_maxNameLength(Cbc_Model * model)
    ;
    /** Fill in first maxLength bytes of name array with a row name */
    COINLIBAPI void COINLINKAGE
    Cbc_getRowName(Cbc_Model * model, int iRow, char * name, size_t maxLength)
    ;
    /** Fill in first maxLength bytes of name array with a column name */
    COINLIBAPI void COINLINKAGE
    Cbc_getColName(Cbc_Model * model, int iColumn, char * name, size_t maxLength)
    ;
    /** Set the name of a column */
    COINLIBAPI void COINLINKAGE
    Cbc_setColName(Cbc_Model * model, int iColumn, const char * name)
    ;
    /** Set the name of a row */
    COINLIBAPI void COINLINKAGE
    Cbc_setRowName(Cbc_Model * model, int iRow, const char * name)
    ;
    /** Number of constraints in the model */
    COINLIBAPI int COINLINKAGE
    Cbc_getNumRows(Cbc_Model * model)
    ;
    /** Number of non-zero entries in a row */
    COINLIBAPI int COINLINKAGE
    Cbc_getRowNz(Cbc_Model * model, int row)
    ;
    /** Indices of variables that appear on this row */
    COINLIBAPI const int * COINLINKAGE
    Cbc_getRowIndices(Cbc_Model * model, int row)
    ;
    /** Coefficients of variables that appear on this row */
    COINLIBAPI const double * COINLINKAGE
    Cbc_getRowCoeffs(Cbc_Model * model, int row)
    ;
    /** Number of non-zero entries in a column */
    COINLIBAPI int COINLINKAGE
    Cbc_getColNz(Cbc_Model * model, int col)
    ;
    /** Indices of rows that a column appears */
    COINLIBAPI const int * COINLINKAGE
    Cbc_getColIndices(Cbc_Model * model, int col)
    ;
    /** Coefficients that a column appear in rows */
    COINLIBAPI const double * COINLINKAGE
    Cbc_getColCoeffs(Cbc_Model * model, int col)
    ;
    /** Right hand side of a row */
    COINLIBAPI double COINLINKAGE
    Cbc_getRowRHS(Cbc_Model * model, int row)
    ;
    /** Sense a row */
    COINLIBAPI char COINLINKAGE
    Cbc_getRowSense(Cbc_Model * model, int row)
    ;
    /** Number of variables in the model */
    COINLIBAPI int COINLINKAGE
    Cbc_getNumCols(Cbc_Model * model)
    ;
    /** Number of integer variables in the model */
    COINLIBAPI int COINLINKAGE
    Cbc_getNumIntegers(Cbc_Model * model)
    ;
    /** Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore) */
    COINLIBAPI void COINLINKAGE
    Cbc_setObjSense(Cbc_Model * model, double sense)
    ;
    /** Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore) */
    COINLIBAPI double COINLINKAGE
    Cbc_getObjSense(Cbc_Model * model)
    ;
    /** Constraint lower bounds */
    COINLIBAPI const double* COINLINKAGE
    Cbc_getRowLower(Cbc_Model * model)
    ;
    /** Set the lower bound of a single constraint */
    COINLIBAPI void COINLINKAGE
    Cbc_setRowLower(Cbc_Model * model, int index, double value)
    ;
    /** Constraint upper bounds */
    COINLIBAPI const double* COINLINKAGE
    Cbc_getRowUpper(Cbc_Model * model)
    ;
    /** Set the upper bound of a single constraint */
    COINLIBAPI void COINLINKAGE
    Cbc_setRowUpper(Cbc_Model * model, int index, double value)
    ;
    /** Objective vector */
    COINLIBAPI const double * COINLINKAGE
    Cbc_getObjCoefficients(Cbc_Model * model)
    ;
    /** Set the objective coefficient of a single variable */
    COINLIBAPI void COINLINKAGE
    Cbc_setObjCoeff(Cbc_Model * model, int index, double value)
    ;
    /** Variable lower bounds */
    COINLIBAPI const double * COINLINKAGE
    Cbc_getColLower(Cbc_Model * model)
    ;
    /** Set the lower bound of a single variable */
    COINLIBAPI void COINLINKAGE
    Cbc_setColLower(Cbc_Model * model, int index, double value)
    ;
    /** Variable upper bounds */
    COINLIBAPI const double * COINLINKAGE
    Cbc_getColUpper(Cbc_Model * model)
    ;
    /** Set the upper bound of a single variable */
    COINLIBAPI void COINLINKAGE
    Cbc_setColUpper(Cbc_Model * model, int index, double value)
    ;
    /** Determine whether the ith variable is integer restricted */
    COINLIBAPI int COINLINKAGE
    Cbc_isInteger(Cbc_Model * model, int i)
    ;
    /** Set this variable to be continuous */
    COINLIBAPI void COINLINKAGE
    Cbc_setContinuous(Cbc_Model * model, int iColumn)
    ;
    /** Set this variable to be integer */
    COINLIBAPI void COINLINKAGE
    Cbc_setInteger(Cbc_Model * model, int iColumn)
    ;
    /** Adds a new column */
    COINLIBAPI void COINLINKAGE
    Cbc_addCol( Cbc_Model *model, const char *name, double lb, 
            double ub, double obj, char isInteger,
            int nz, int *rows, double *coefs )
    ;
    /** Adds a new row */
    COINLIBAPI void COINLINKAGE
    Cbc_addRow( Cbc_Model *model, const char *name, int nz,
            const int *cols, const double *coefs, char sense, double rhs )
    ;
    /** Add SOS constraints to the model using row-order matrix */
    COINLIBAPI void  COINLINKAGE
    Cbc_addSOS(Cbc_Model * model, int numRows, const int * rowStarts,
               const int * colIndices, const double * weights, const int type)
    ;
    /** Print the model */
    COINLIBAPI void COINLINKAGE
    Cbc_printModel(Cbc_Model * model, const char * argPrefix)
    ;
    /** Return a copy of this model */
    COINLIBAPI Cbc_Model * COINLINKAGE
    Cbc_clone(Cbc_Model * model)
    ;
    /*@}*/
    /**@name Solver parameters */
    /*@{*/
    /** Set parameter "name" to value "value". Note that this
     * translates directly to using "-name value" as a 
     * command-line argument to Cbc.*/
    COINLIBAPI void COINLINKAGE
    Cbc_setParameter(Cbc_Model * model, const char * name, const char * value)
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

    /*@}*/


    /**@name Solving the model */
    /*@{*/
    /* Solve the model with Cbc (using CbcMain1).
    */
    COINLIBAPI int COINLINKAGE
    Cbc_solve(Cbc_Model * model)
    ;
    /*@}*/


    /**@name Accessing the solution and solution status */
    /*@{*/

    /** Sum of primal infeasibilities */
    COINLIBAPI double COINLINKAGE
    Cbc_sumPrimalInfeasibilities(Cbc_Model * model)
    ;
    /** Number of primal infeasibilities */
    COINLIBAPI int COINLINKAGE
    Cbc_numberPrimalInfeasibilities(Cbc_Model * model)
    ;

    /** Just check solution (for external use) - sets sum of
        infeasibilities etc */
    COINLIBAPI void COINLINKAGE
    Cbc_checkSolution(Cbc_Model * model)
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
    /** Is infeasiblity proven (or none better than cutoff)? */
    COINLIBAPI int COINLINKAGE
    Cbc_isProvenInfeasible(Cbc_Model * model)
    ;
    /** Was continuous solution unbounded? */
    COINLIBAPI int COINLINKAGE
    Cbc_isContinuousUnbounded(Cbc_Model * model)
    ;
    /** Node limit reached? */
    COINLIBAPI int COINLINKAGE
    Cbc_isNodeLimitReached(Cbc_Model * model)
    ;
    /** Time limit reached? */
    COINLIBAPI int COINLINKAGE
    Cbc_isSecondsLimitReached(Cbc_Model * model)
    ;
    /** Solution limit reached? */
    COINLIBAPI int COINLINKAGE
    Cbc_isSolutionLimitReached(Cbc_Model * model)
    ;
    /** Are there numerical difficulties (for initialSolve) ? */
    COINLIBAPI int COINLINKAGE
    Cbc_isInitialSolveAbandoned(Cbc_Model * model)
    ;
    /** Is optimality proven (for initialSolve) ? */
    COINLIBAPI int COINLINKAGE
    Cbc_isInitialSolveProvenOptimal(Cbc_Model * model)
    ;
    /** Is primal infeasiblity proven (for initialSolve) ? */
    COINLIBAPI int COINLINKAGE
    Cbc_isInitialSolveProvenPrimalInfeasible(Cbc_Model * model)
    ;
    /** "row" solution
     *  This is the vector A*x, where A is the constraint matrix
     *  and x is the current solution. */
    COINLIBAPI const double * COINLINKAGE
    Cbc_getRowActivity(Cbc_Model * model)
    ;
    /** Best feasible solution vector */
    COINLIBAPI const double * COINLINKAGE
    Cbc_getColSolution(Cbc_Model * model)
    ;
    /** Objective value of best feasible solution */
    COINLIBAPI double COINLINKAGE
    Cbc_getObjValue(Cbc_Model * model)
    ;
    /** Best known bound on the optimal objective value */
    COINLIBAPI double COINLINKAGE
    Cbc_getBestPossibleObjValue(Cbc_Model * model)
    ;
    /** Best integer feasible solution or NULL if no integer feas sol found */
    COINLIBAPI double*  COINLINKAGE
    Cbc_bestSolution(Cbc_Model * model)
    ;
    /** Number of nodes explored in B&B tree */
    COINLIBAPI int COINLINKAGE
    Cbc_getNodeCount(Cbc_Model * model)
    ;
    /** Print the solution */
    COINLIBAPI void  COINLINKAGE
    Cbc_printSolution(Cbc_Model * model)
    ;
    /** Final status of problem
        Some of these can be found out by is...... functions
        -1 before branchAndBound
        0 finished - check isProvenOptimal or isProvenInfeasible to see if solution found
        (or check value of best solution)
        1 stopped - on maxnodes, maxsols, maxtime
        2 difficulties so run was abandoned
        (5 event user programmed event occurred)
    */
    COINLIBAPI int COINLINKAGE
    Cbc_status(Cbc_Model * model)
    ;
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
    COINLIBAPI int COINLINKAGE
    Cbc_secondaryStatus(Cbc_Model * model)
    ;
    /*@}*/
#ifdef __cplusplus
}
#endif
#endif
