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

/** This is a first "C" interface to Cbc.
    It is mostly similar to the "C" interface to Clp and
    was contributed by Bob Entriken.
*/

#ifdef __cplusplus
extern "C" {
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
    Cbc_readMps(Cbc_Model * model, const char *filename)
    ;
    /** Write an mps file from the given filename */
    COINLIBAPI void COINLINKAGE
    Cbc_writeMps(Cbc_Model * model, const char *filename)
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
    /** Set parameter "name" to value "value". Note that this
     * translates directly to using "-name value" as a 
     * command-line argument to Cbc.*/
    COINLIBAPI void COINLINKAGE
    Cbc_setParameter(Cbc_Model * model, const char * name, const char * value)
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
    /* Solve using CbcMain1. This is the recommended default solve function.
    */
    COINLIBAPI int COINLINKAGE
    Cbc_solve(Cbc_Model * model)
    ;
    /*@}*/


    /**@name most useful gets and sets */
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
    /** Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore) */
    COINLIBAPI void COINLINKAGE
    Cbc_setObjSense(Cbc_Model * model, double sense)
    ;
    /** Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore) */
    COINLIBAPI double COINLINKAGE
    Cbc_getObjSense(Cbc_Model * model)
    ;
    /** Primal row solution
     *  This is the vector A*x, where A is the constraint matrix
     *  and x is the current solution. */
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
    /** Set this the variable to be integer */
    COINLIBAPI Cbc_Model * COINLINKAGE
    Cbc_setInteger(Cbc_Model * model, int iColumn)
    ;
    /** Add SOS constraints to the model using dense matrix */
    COINLIBAPI void  COINLINKAGE
    Cbc_addSOS_Dense(Cbc_Model * model, int numObjects, const int * len,
                     const int * const * which, const double * weights, const int type)
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
    /*@}*/
#ifdef __cplusplus
}
#endif
#endif
