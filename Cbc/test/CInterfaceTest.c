/* $Id$ */
/* Copyright (C) 2014, International Business Machines
   Corporation and others.  All Rights Reserved.
   This code is licensed under the terms of the Eclipse Public License (EPL). */

#undef NDEBUG /* force asserts to work */
#include "Cbc_C_Interface.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef INFINITY /* workaround for non-C99 compilers */
#define INFINITY (HUGE_VAL * 2)
#endif


static int callback_called = 0;

void (COINLINKAGE_CB test_callback)(Cbc_Model * model,int  msgno, int ndouble,
                            const double * dvec, int nint, const int * ivec,
                            int nchar, char ** cvec) {

    callback_called = 1;
    printf("In callback: message %d\n", msgno);

}


void testKnapsack() {

    Cbc_Model *model = Cbc_newModel();

    /* Simple knapsack problem
       Maximize  5x[1] + 3x[2] + 2x[3] + 7x[4] + 4x[5]
       s.t.      2x[1] + 8x[2] + 4x[3] + 2x[4] + 5x[5] <= 10
       All x binary
       */
    
    CoinBigIndex start[] = {0, 1, 2, 3, 4, 5, 6};
    int rowindex[] = {0, 0, 0, 0, 0};
    double value[] = {2, 8, 4, 2, 5};
    double collb[] = {0,0,0,0,0};
    double colub[] = {1,1,1,1,1};
    double obj[] = {5, 3, 2, 7, 4};
    double feasible[] = {1,1,0,0,0};
    double rowlb[] = {-INFINITY};
    double rowub[] = {10};
    const double *sol;
    const char* setname = "test model";
    char *getname = malloc(20);
    int i;

    printf("Interface reports Cbc version %s\n", Cbc_getVersion());

    Cbc_loadProblem(model, 5, 1, start, rowindex, value, collb, colub, obj, rowlb, rowub);

    Cbc_setColName(model, 2, "var2");
    Cbc_setRowName(model, 0, "constr0");


    assert(Cbc_getNumCols(model) == 5);
    assert(Cbc_getNumRows(model) == 1);

    for (i = 0; i < 5; i++) {
        Cbc_setInteger(model, i);
        assert(Cbc_isInteger(model,i));
    }

    Cbc_setObjSense(model, -1);
    assert(Cbc_getObjSense(model) == -1);

    Cbc_setProblemName(model, setname);

    Cbc_registerCallBack(model, test_callback);

    Cbc_setInitialSolution(model, feasible);

    Cbc_solve(model);

    assert(Cbc_isProvenOptimal(model));
    assert(!Cbc_isAbandoned(model));
    assert(!Cbc_isProvenInfeasible(model));
    assert(!Cbc_isContinuousUnbounded(model));
    assert(!Cbc_isNodeLimitReached(model));
    assert(!Cbc_isSecondsLimitReached(model));
    assert(!Cbc_isSolutionLimitReached(model));
    assert(fabs( Cbc_getObjValue(model)- (16.0) < 1e-6));
    assert(fabs( Cbc_getBestPossibleObjValue(model)- (16.0) < 1e-6));

    assert(callback_called == 1);
    
    sol = Cbc_getColSolution(model);
    
    assert(fabs(sol[0] - 1.0) < 1e-6);
    assert(fabs(sol[1] - 0.0) < 1e-6);
    assert(fabs(sol[2] - 0.0) < 1e-6);
    assert(fabs(sol[3] - 1.0) < 1e-6);
    assert(fabs(sol[4] - 1.0) < 1e-6);

    Cbc_problemName(model, 20, getname);
    i = strcmp(getname,setname);
    assert( (i == 0) );

    Cbc_getColName(model, 2, getname, 20);
    i = strcmp(getname, "var2");
    assert( (i == 0) );
    Cbc_getRowName(model, 0, getname, 20);
    i = strcmp(getname, "constr0");
    assert( (i == 0) );
    assert( Cbc_maxNameLength(model) >= 7 );
    
    Cbc_deleteModel(model);
    free(getname);

}

/*
void testProblemModification() {

    Cbc_Model *model = Cbc_newModel();

    / * Simple knapsack problem
       Maximize  5x[1] + 3x[2] + 2x[3] + 7x[4] + 4x[5]
       s.t.      2x[1] + 8x[2] + 4x[3] + 2x[4] + 5x[5] <= 10
       All x binary
       * /
    
    CoinBigIndex start[] = {0, 1, 2, 3, 4, 5, 6};
    int rowindex[] = {0, 0, 0, 0, 0};
    double value[] = {2, 8, 4, 2, 5};
    double collb[] = {0,0,0,0,0};
    double colub[] = {1,1,1,1,1};
    double obj[] = {5, 3, 2, 7, 4};
    double rowlb[] = {-INFINITY};
    double rowub[] = {10};
    const double *sol;
    int i;

    printf("Interface reports Cbc version %s\n", Cbc_getVersion());

    Cbc_loadProblem(model, 5, 1, start, rowindex, value, collb, colub, obj, rowlb, rowub);

    for (i = 0; i < 5; i++) {
        Cbc_setInteger(model, i);
        assert(Cbc_isInteger(model,i));
    }

    Cbc_setObjSense(model, -1);
    assert(Cbc_getObjSense(model) == -1);

    Cbc_solve(model);

    assert(Cbc_isProvenOptimal(model));
    assert(fabs( Cbc_getObjValue(model)- (16.0) < 1e-6));

    sol = Cbc_getColSolution(model);
    
    assert(fabs(sol[0] - 1.0) < 1e-6);
    assert(fabs(sol[1] - 0.0) < 1e-6);
    assert(fabs(sol[2] - 0.0) < 1e-6);
    assert(fabs(sol[3] - 1.0) < 1e-6);
    assert(fabs(sol[4] - 1.0) < 1e-6);

    Cbc_setColUpper(model, 0, 0.0);
    Cbc_solve(model);

    assert(Cbc_isProvenOptimal(model));
    assert(fabs( Cbc_getObjValue(model)- (11.0) < 1e-6));

    sol = Cbc_getColSolution(model);
    
    assert(fabs(sol[0] - 0.0) < 1e-6);
    assert(fabs(sol[1] - 0.0) < 1e-6);
    assert(fabs(sol[2] - 0.0) < 1e-6);
    assert(fabs(sol[3] - 1.0) < 1e-6);
    assert(fabs(sol[4] - 1.0) < 1e-6);


    Cbc_setColLower(model, 1, 1.0);

    assert(Cbc_isProvenOptimal(model));
    assert(fabs( Cbc_getObjValue(model)- (10.0) < 1e-6));

    sol = Cbc_getColSolution(model);
    
    assert(fabs(sol[0] - 0.0) < 1e-6);
    assert(fabs(sol[1] - 1.0) < 1e-6);
    assert(fabs(sol[2] - 0.0) < 1e-6);
    assert(fabs(sol[3] - 1.0) < 1e-6);
    assert(fabs(sol[4] - 0.0) < 1e-6);

    
    Cbc_deleteModel(model);

}
*/


void testSOS() {

    Cbc_Model *model = Cbc_newModel();

    /*
       Maximize  5x[1] + 3x[2] + 2x[3] + 7x[4] + 4x[5]
       s.t.       x[1] +  x[2] +  x[3] +  x[4] +  x[5] == 1
       All x binary
       */
    
    CoinBigIndex start[] = {0, 0, 0, 0, 0, 0, 0};
    double collb[] = {0,0,0,0,0};
    double colub[] = {1,1,1,1,1};
    double obj[] = {5, 3, 2, 7, 4};
    int sosrowstart[] = {0,5};
    int soscolindex[] = {0,1,2,3,4}; 
    const double *sol;
    int i;

    Cbc_loadProblem(model, 5, 0, start, NULL, NULL, collb, colub, obj, NULL, NULL);

    assert(Cbc_getNumCols(model) == 5);
    assert(Cbc_getNumRows(model) == 0);

    for (i = 0; i < 5; i++) {
        Cbc_setInteger(model, i);
        assert(Cbc_isInteger(model,i));
    }

    Cbc_setObjSense(model, -1);
    assert(Cbc_getObjSense(model) == -1);
    
    Cbc_addSOS(model,1,sosrowstart,soscolindex,obj,1); 

    Cbc_solve(model);

    assert(Cbc_isProvenOptimal(model));
    assert(!Cbc_isAbandoned(model));
    assert(!Cbc_isProvenInfeasible(model));
    assert(!Cbc_isContinuousUnbounded(model));
    assert(!Cbc_isNodeLimitReached(model));
    assert(!Cbc_isSecondsLimitReached(model));
    assert(!Cbc_isSolutionLimitReached(model));
    assert(fabs( Cbc_getObjValue(model)- (7.0) < 1e-6));
    assert(fabs( Cbc_getBestPossibleObjValue(model)- (7.0) < 1e-6));

    sol = Cbc_getColSolution(model);
    
    assert(fabs(sol[0] - 0.0) < 1e-6);
    assert(fabs(sol[1] - 0.0) < 1e-6);
    assert(fabs(sol[2] - 0.0) < 1e-6);
    assert(fabs(sol[3] - 1.0) < 1e-6);
    assert(fabs(sol[4] - 0.0) < 1e-6);

    Cbc_deleteModel(model);

}


void testIntegerInfeasible() {

    Cbc_Model *model = Cbc_newModel();

    /* Minimize x
     * s.t.     x <= -10
     * x binary */

    CoinBigIndex start[] = {0, 1};
    int rowindex[] = {0};
    double value[] = {1.0};
    double rowlb[] = {-INFINITY};
    double rowub[] = {-10};

    double collb[] = {0.0};
    double colub[] = {1.0};
    double obj[] = {1.0};

    Cbc_loadProblem(model, 1, 1, start, rowindex, value, collb, colub, obj, rowlb, rowub);

    Cbc_setInteger(model, 0);

    assert(Cbc_getNumCols(model) == 1);
    assert(Cbc_getNumRows(model) == 1);

    Cbc_solve(model);
    
    assert(!Cbc_isProvenOptimal(model));
    assert(Cbc_isProvenInfeasible(model));
    
    Cbc_deleteModel(model);

}

void testIntegerUnbounded() {

    Cbc_Model *model = Cbc_newModel();

    /* http://list.coin-or.org/pipermail/cbc/2014-March/001276.html
     * Minimize x
     * s.t. x + y <= 3
     *      x - y == 0
     *      x,y Free
     *      x integer */

    CoinBigIndex start[] = {0,2,4};
    int rowindex[] = {0, 1, 0, 1};
    double value[] = {1, 1, 1, -1};
    double rowlb[] = {-INFINITY, 0.0};
    double rowub[] = {3.0,0.0};
    double collb[] = {-INFINITY, -INFINITY};
    double colub[] = {INFINITY, INFINITY};
    double obj[] = {1.0,0.0};

    Cbc_loadProblem(model, 2, 2, start, rowindex, value, collb, colub, obj, rowlb, rowub);

    Cbc_setInteger(model, 0);

    Cbc_setParameter(model, "log", "0");
    
    printf("About to solve problem silently. You should see no output except \"Done\".\n");
    Cbc_solve(model);
    printf("Done\n");
    
    assert(!Cbc_isProvenOptimal(model));
    assert(!Cbc_isProvenInfeasible(model));
    assert(Cbc_isContinuousUnbounded(model));

    Cbc_deleteModel(model);


}

void testIntegerBounds() {
    /* max 1.1x + 100.0z
       st     x +      z <= 3
         0 <= x <= 3
         0 <= z <= 1.5, Integer
     x* = 2, z* = 1 */
    
    Cbc_Model *model = Cbc_newModel();

    CoinBigIndex start[] = {0,1,2};
    int rowindex[] = {0, 0};
    double value[] = {1, 1};
    double rowlb[] = {-INFINITY};
    double rowub[] = {3.0};
    double collb[] = {0.0, 0.0};
    double colub[] = {3.0, 1.5};
    double obj[] = {1.1,100.0};
    const double *sol;

    Cbc_loadProblem(model, 2, 1, start, rowindex, value, collb, colub, obj, rowlb, rowub);

    Cbc_setInteger(model, 1);
    Cbc_setObjSense(model, -1);

    Cbc_solve(model);
    
    assert(Cbc_isProvenOptimal(model));

    sol = Cbc_getColSolution(model);
    
    assert(fabs(sol[0] - 2.0) < 1e-6);
    assert(fabs(sol[1] - 1.0) < 1e-6);

    Cbc_deleteModel(model);

}

void testQueens(int n) {
    int **x = malloc( sizeof(int*)*n );
    int i, j, k, p;
    x[0] = malloc( sizeof(int)*n*n );
    for ( i=1 ; (i<n) ; ++i )
        x[i] = x[i-1] + n;

    Cbc_Model *m = Cbc_newModel();
    
    // adding variables
    k = 0;
    for ( i=0 ; (i<n) ; ++i )
    {
        for ( j=0 ; (j<n) ; ++j )
        {
            x[i][j] = k++;
            char name[256];
            sprintf(name, "x(%d,%d)", i, j);
            Cbc_addCol(m, name, 0.0, 1.0, 0.0, 1, 0, NULL, NULL);
        }
    }

    int *idx = malloc(sizeof(int)*n);
    double *coef = malloc(sizeof(double)*n);

    // constraint one per row
    for ( i=0 ; (i<n) ; ++i )
    {
        for ( j=0 ; j<n ; ++j )
        {
            idx[j] = x[i][j];
            coef[j] = 1.0;
        }
        char name[256];
        sprintf(name, "row(%d)", i);
        Cbc_addRow(m, name, n, idx, coef, 'E', 1.0);
    }

    // constraint one per column
    for ( j=0 ; (j<n) ; ++j )
    {
        for ( i=0 ; i<n ; ++i )
        {
            idx[i] = x[i][j];
            coef[i] = 1.0;
        }
        char name[256];
        sprintf(name, "col(%d)", j);
        Cbc_addRow(m, name, n, idx, coef, 'E', 1.0);
    }

    // diagonal 
    p = 0;
    for ( k=2-n ; k<(n-1) ; ++k, ++p )
    {
        int nz = 0;
        for ( i=0 ; (i<n) ; ++i )
        {
            for ( j=0 ; (j<n) ; ++j )
            {
                if (i-j==k)
                {
                    idx[nz] = x[i][j];
                    coef[nz] = 1.0;
                    ++nz;
                }
            }
        }
        char name[256];
        sprintf(name, "diag1(%d)", k);
        Cbc_addRow(m, name, nz, idx, coef, 'L', 1.0);
    }

    // diagonal /
    p = 0;
    for ( k=3 ; k<(n+n) ; ++k, ++p )
    {
        int nz = 0;
        for ( i=0 ; (i<n) ; ++i )
        {
            for ( j=0 ; (j<n) ; ++j )
            {
                if (i+j==k)
                {
                    idx[nz] = x[i][j];
                    coef[nz] = 1.0;
                    ++nz;
                }
            }
        }
        char name[256];
        sprintf(name, "diag2(%d)", k);
        Cbc_addRow(m, name, nz, idx, coef, 'L', 1.0);
    } 

    Cbc_setMaximumSeconds(m, 60);
    int status = Cbc_solve(m);
    const double *xs = Cbc_getColSolution(m);
    if (n<=75)
    {
        // should find the optimal for small problems
        assert(Cbc_isProvenOptimal(m));
        assert(xs);
    }
    if (xs) {
        // solution check

        // total number of queens
        int nq = 0;
        for ( i=0 ; (i<n) ; ++i )
            for ( j=0 ; (j<n) ; ++j )
                if ((fabs(xs[x[i][j]]-1.0))<1e-5)
                    nq++;
        assert(nq == n);
        // one per row
        for ( i=0 ; (i<n) ; ++i )
        {
            nq = 0;
            for ( j=0 ; (j<n) ; ++j )
                if ((fabs(xs[x[i][j]]-1.0))<1e-5)
                    nq++;
            assert( nq == 1);
        }
        // one per column
        for ( j=0 ; (j<n) ; ++j )
        {
            nq = 0;
            for ( i=0 ; (i<n) ; ++i )
                if ((fabs(xs[x[i][j]]-1.0))<1e-5)
                    nq++;
            assert( nq == 1);
        }
    }

    free(idx);
    free(coef);
    free(x[0]);
    free(x);

    Cbc_deleteModel(m);
}

int main() {

    printf("Knapsack test\n");
    testKnapsack();
    printf("SOS test\n");
    testSOS();
    printf("Infeasible test\n");
    testIntegerInfeasible();
    printf("Unbounded test\n");
    testIntegerUnbounded();
    /*printf("Problem modification test\n");
    testProblemModification();*/
    printf("Integer bounds test\n");
    testIntegerBounds();
    printf("n-Queens test\n");
    testQueens(10);
    testQueens(25);
    testQueens(50);
    testQueens(75);

    return 0;
}
