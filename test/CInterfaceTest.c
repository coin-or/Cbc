/* $Id$ */
/* Copyright (C) 2014, International Business Machines
   Corporation and others.  All Rights Reserved.
   This code is licensed under the terms of the Eclipse Public License (EPL). */

#undef NDEBUG /* force asserts to work */
#include "Cbc_C_Interface.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#ifndef INFINITY /* workaround for non-C99 compilers */
#define INFINITY (HUGE_VAL * 2)
#endif


static int callback_called = 0;

void (CBC_LINKAGE_CB test_callback)(Cbc_Model * model,int  msgno, int ndouble,
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
    assert(fabs(Cbc_getObjValue(model)- (16.0)) < 1e-6);
    assert(fabs(Cbc_getBestPossibleObjValue(model)- 16.0) < 1e-6);

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
    assert(fabs(Cbc_getObjValue(model)- 7.0) < 1e-6);
    assert(fabs(Cbc_getBestPossibleObjValue(model)-7.0) < 1e-6);

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
    Cbc_Model *model;
    int *idx;
    double *coef;
    const double *xs;
    int **x = malloc( sizeof(int*)*n );
    int i, j, k, p;

    x[0] = malloc( sizeof(int)*n*n );
    for ( i=1 ; (i<n) ; ++i )
        x[i] = x[i-1] + n;

    model = Cbc_newModel();
    
    /* adding variables */
    k = 0;
    for ( i=0 ; (i<n) ; ++i )
    {
        for ( j=0 ; (j<n) ; ++j )
        {
            char name[256];
            x[i][j] = k++;
            sprintf(name, "x(%d,%d)", i, j);
            Cbc_addCol(model, name, 0.0, 1.0, 0.0, 1, 0, NULL, NULL);
        }
    }

    idx = malloc(sizeof(int)*n);
    coef = malloc(sizeof(double)*n);

    /* constraint one per row */
    for ( i=0 ; (i<n) ; ++i )
    {
        char name[256];
        for ( j=0 ; j<n ; ++j )
        {
            idx[j] = x[i][j];
            coef[j] = 1.0;
        }
        sprintf(name, "row(%d)", i);
        Cbc_addRow(model, name, n, idx, coef, 'E', 1.0);
    }

    /* constraint one per column */
    for ( j=0 ; (j<n) ; ++j )
    {
        char name[256];
        for ( i=0 ; i<n ; ++i )
        {
            idx[i] = x[i][j];
            coef[i] = 1.0;
        }
        sprintf(name, "col(%d)", j);
        Cbc_addRow(model, name, n, idx, coef, 'E', 1.0);
    }

    /* diagonal  */
    p = 0;
    for ( k=2-n ; k<(n-1) ; ++k, ++p )
    {
        char name[256];
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
        sprintf(name, "diag1(%d)", k);
        char *s = name;
        while (*s != '\0') {
            if (*s == '-')
                *s = 'm';
            ++s;
        }
        Cbc_addRow(model, name, nz, idx, coef, 'L', 1.0);
    }

    /* diagonal */
    p = 0;
    for ( k=3 ; k<(n+n) ; ++k, ++p )
    {
        char name[256];
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
        sprintf(name, "diag2(%d)", k);
        char *s = name;
        while (*s != '\0') {
            if (*s == '-')
                *s = 'm';
            ++s;
        }
 
        Cbc_addRow(model, name, nz, idx, coef, 'L', 1.0);
    } 

    Cbc_setMaximumSeconds(model, 100);
    (void) Cbc_solve(model);
    xs = Cbc_getColSolution(model);
    if (n<=30)
    {
        /* should find the optimal for small problems */
        assert(Cbc_isProvenOptimal(model));
        assert(xs);
    }
    if (xs) {
        /* solution check 

         total number of queens */
        int nq = 0;
        for ( i=0 ; (i<n) ; ++i )
            for ( j=0 ; (j<n) ; ++j )
                if ((fabs(xs[x[i][j]]-1.0))<1e-5)
                    nq++;
        assert(nq == n);
        for ( i=0 ; (i<n) ; ++i )
        {
            nq = 0;
            for ( j=0 ; (j<n) ; ++j )
                if ((fabs(xs[x[i][j]]-1.0))<1e-5)
                    nq++;
            assert( nq == 1);
        }
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

    Cbc_setProblemName(model, "CrazyQueens");
    Cbc_writeMps(model, "q");
    Cbc_writeLp(model, "q");

    Cbc_deleteModel(model);
}

/* asMIP 0 solves only the LP
 *       1 solves as MIP */
void testTSP(char asMIP) {
    #define N 7
    int oo = INT_MAX;

    /* distance matrix */
    const double d[][N] = 
            {      /* a   b   c   d   e   f   g*/
          /* a */  { oo, 49, 80, 56, oo, oo, 47},
          /* b */  { 50, oo, oo, 37, 21, oo, 25},
          /* c */  { 99, oo, oo, 52, oo, 35, oo},
          /* d */  { 67, 39, 37, oo, 15, oo, oo},
          /* e */  { oo, 30, oo, 20, oo, 20, 49},
          /* f */  { oo, oo, 35, oo, 20, oo, 32},
          /* g */  { 68, 35, oo, oo, 38, 37, oo},
            };

    /* variable indexes */
    int x[N][N];
    int y[N];
    int i, j, k;
    int idx[N];
    double coef[N];
    char name[256];
    int ia;
    int nz = 0;
    int newConstraints;
    double sum;
    double opt;
    int arcs[6][2];
    int nArcs;

    Cbc_Model *m = Cbc_newModel();

    /* x variables */
    for ( i=0 ; (i<N) ; ++i ) {
        for ( j=0 ; j<N ; ++j ) {
            if (d[i][j] == oo) {
                x[i][j] = -1;
            } else {
                snprintf(name, 256, "x(%d,%d)", i, j);
                x[i][j] = Cbc_getNumCols(m);
                Cbc_addCol(m, name, 0.0, 1.0, d[i][j], asMIP, 0, NULL, NULL);
            }
        }
    }
    for ( i=0 ; (i<N) ; ++i ) {
        snprintf(name, 256, "y(%d)", i);
        y[i] = Cbc_getNumCols(m);
        Cbc_addCol(m, name, 0.0, N, 0.0, asMIP, 0, NULL, NULL);
    }

    /* outbound arc selection */
    for ( i=0 ; (i<N) ; ++i ) {
        nz = 0;
        for ( j=0 ; (j<N) ; ++j ) {
            if (d[i][j] == oo)
                continue;
            coef[nz] = 1.0;
            idx[nz++] = x[i][j];
        }
        snprintf(name, 256, "out(%d)", i);
        Cbc_addRow(m, name, nz, idx, coef, 'E', 1.0);
    }

    /* inbound arc selection */
    for ( j=0 ; (j<N) ; ++j ) {
        nz = 0;
        for ( i=0 ; (i<N) ; ++i ) {
            if (d[i][j] == oo)
                continue;
            coef[nz] = 1.0;
            idx[nz++] = x[i][j];
        }
        snprintf(name, 256, "in(%d)", j);
        Cbc_addRow(m, name, nz, idx, coef, 'E', 1.0);
    }

    /* weak sub-tour elimination constraints */
    for ( i=1 ; (i<N) ; ++i ) {
        for ( j=1 ; (j<N) ; ++j ) {
            nz = 0;

            if (d[i][j]==oo)
                continue;

            coef[nz] = 1.0;
            idx[nz++] = y[i];

            coef[nz] = -1.0;
            idx[nz++] = y[j];

            coef[nz] = -(N+1);
            idx[nz++] = x[i][j];

            snprintf(name, 256, "from(%d)to(%d)", i, j);
            Cbc_addRow(m, name, nz, idx, coef, 'G', -N);
        }
    }

    Cbc_solve(m);
    assert(Cbc_isProvenOptimal(m));
    opt = asMIP ? 262 : 238.75;
    assert( fabs(Cbc_getObjValue(m)-opt) <= 1e-4 );

    const double *s = Cbc_getColSolution(m);

    if (!asMIP) {

        do {
            newConstraints = 0;
            /* eliminating subtours of size 2 and 3 and reoptimize */
            for ( i=0 ; (i<N) ; ++i ) {
                for ( j=i+1 ; j<N ; ++j ) {
                    nArcs = 0;
                    if (d[i][j] == oo)
                        continue;
                    arcs[nArcs][0] = i;
                    arcs[nArcs++][1] = j;
                    sum = s[x[i][j]];
                    if (d[j][i] != oo) {
                        arcs[nArcs][0] = j;
                        arcs[nArcs++][1] = i;
                        sum += s[x[j][i]];
                        if (sum > 1.01) {
                            idx[0] = x[i][j];
                            idx[1] = x[j][i];
                            coef[0] = 1.0;
                            coef[1] = 1.0;
                            snprintf(name, 256, "noSub(%d,%d)", i, j);
                            Cbc_addRow(m, name, 2, idx, coef, 'L', 1.0);
                            ++newConstraints;
                        }
                    }
                    for ( k=j+1 ; k<N ; ++k ) {
                        int pNArcs = nArcs;
                        double pSum = sum;
                        if (d[i][k] != oo) {
                            arcs[nArcs][0]   = i;
                            arcs[nArcs++][1] = k;
                            sum += s[x[i][k]];
                        }
                        if (d[k][i] != oo) {
                            arcs[nArcs][0]   = k;
                            arcs[nArcs++][1] = i;
                            sum += s[x[k][i]];
                        }
                        if (d[j][k] != oo) {
                            arcs[nArcs][0]   = j;
                            arcs[nArcs++][1] = k;
                            sum += s[x[j][k]];
                        }
                        if (d[k][j] != oo) {
                            arcs[nArcs][0]   = k;
                            arcs[nArcs++][1] = j;
                            sum += s[x[k][j]];
                        }

                        if (sum >= 2.01)
                        {
                            coef[0] = coef[1] = coef[2] = 1.0;
                            coef[3] = coef[4] = coef[5] = 1.0;

                            for ( ia=0 ; (ia<nArcs) ; ++ia )
                                idx[ia] = x[arcs[ia][0]][arcs[ia][1]];

                            snprintf(name, 256, "noSub(%d,%d,%d)", i, j, k);
                            Cbc_addRow(m, name, nArcs, idx, coef, 'L', 2.0);
                            ++newConstraints;
                        }
                        nArcs = pNArcs;
                        sum = pSum;
                    }
                }
            }
            
            printf("Model strengthened with %d sub-tour elimination constraints, reoptimizing it.", newConstraints);
            Cbc_solve(m);
            assert(Cbc_isProvenOptimal(m));
            printf("New bound now %g\n", Cbc_getObjValue(m));
        } while (newConstraints); /* reoptimize loop */

        assert( fabs(Cbc_getObjValue(m)-261) <= 1e-4 );
    } /* initially not as MIP */

    Cbc_deleteModel(m);

    #undef N
}

static char arc_var(const char *vname, int *u, int *v) {
    // gets  arc in the name of a x variable in the format x(u,v)
    char str[256];
    strcpy(str, vname);

    char *s = strstr(str, "x(");
    if (!s)
        return 0;
    s += 2;
    char *s2 = strstr(s, ",");
    if (!s2)
        return 0;
    *s2 = '\0';
    ++s2;

    char *s3 = strstr(s2, ")");
    if (!s3)
        return 0;
    *s3 = '\0';

    *u = atoi(s);
    *v = atoi(s2);

    return 1;
}

static double rad(double x)
{     /* convert input coordinate to longitude/latitude, in radians */
      double pi = 3.141592, deg, min;
      deg = (int)x;
      min = x - deg;
      return pi * (deg + 5.0 * min / 3.0) / 180.0;
}

double dist(double x1, double y1, double x2, double y2) {
    double rrr = 6378.388;
    double latitude_i = rad(x1);
    double latitude_j = rad(x2);
    double longitude_i = rad(y1);
    double longitude_j = rad(y2);
    double q1 = cos(longitude_i - longitude_j);
    double q2 = cos(latitude_i - latitude_j);
    double q3 = cos(latitude_i + latitude_j);
    double dij = (int)(rrr * acos(0.5 * ((1.0 + q1) * q2 -
       (1.0 - q1) *q3)) + 1.0);
    return dij;
}

// checks a solution for the TSP, assuming that x variables 
// are the first ones the the graph is complete
// starting in st computes all other points in the route and stores in 
// el - can be used to identify subtours in integer solutions
static int tspRouteStarting(const double **xsol, int st, int n, int *el) {
    int size = 0;
    int next = st;
    do {
        el[size++] = next;
        char found = 0;
        for ( int i=0 ; (i<n) ; ++i ) {
            if (xsol[next][i] >= 0.99) {
                next = i;
                found = 1;
                break;
            }
        }
        if (!found) {
            fprintf(stderr, "TSP solution does not satisfy degree constraints.");
            abort();
        }
    } while ( next != st );

    return size;
}

void subTourSep(void *osiSolver, void *osiCuts, void *appdata) {
    if (!Osi_isProvenOptimal(osiSolver))
        return;

    int n = *((int *)appdata);
    
    const double *x = Osi_getColSolution(osiSolver);
    double **xs = malloc( sizeof(double*)*n );
    xs[0] = malloc( sizeof(double)*n*n );
    for ( int i=1 ; (i<n) ; ++i )
        xs[i] = xs[i-1] + n;
    for ( int i=0 ; (i<(n*n)) ; ++i )
        xs[0][i] = 0.0;

    int nFrac = 0;
    for ( int i=0 ; (i<Osi_getNumCols(osiSolver)) ; ++i ) {
        if (fabs(x[i]) <= 1e-4)
            continue;
        char vname[256] = "";
        Osi_getColName(osiSolver, i, vname, 256);
        int u = -1, v = -1;

        if (!arc_var(vname, &u, &v))
            continue;

        xs[u][v] = x[i];

        if (fabs(x[i] - floor(x[i]+0.5)) > 1e-5)
            nFrac++;
    }
    if (nFrac == 0) {
        // integer sol, search for 
        // disconnected sub-routes with DFS
        //
        int *el = malloc(sizeof(int)*n);
        char *iv = malloc(sizeof(char)*n);
        int *idx = malloc(sizeof(int)*n*n);
        double *coef = malloc(sizeof(double)*n*n);
        for ( int i=0 ; (i<n*n) ; ++i )
            coef[i] = 1.0;

        for ( int st=0 ; (st<n) ; ++st ) {
            int nz = 0;
            memset(iv, 0, sizeof(char)*n);
            int nel = tspRouteStarting((const double **) xs, st, n, el);
            if (nel == n)
                break; // no sub-tour
            else {
                if (nel <= (n/2)) {
                    // only for the smaller subset
                    for ( int j=0 ; (j<nel) ; ++j )
                        iv[el[j]] = 1;
                    for ( int ic=0 ; (ic<Osi_getNumCols(osiSolver)) ; ++ic ) {
                        char vname[256] = "";
                        Osi_getColName(osiSolver, ic, vname, 256);
                        int u = -1, v = -1;

                        if (!arc_var(vname, &u, &v))
                            continue;

                        if (iv[u] + iv[v] != 2)
                            continue;
                        
                        // both in subset
                        idx[nz++] = ic;
                    }
                    OsiCuts_addGlobalRowCut( osiCuts, nz, idx, coef, 'L', nel-1 );
                }
            } // found subroute
        } // checking for subroutes
        free(iv);
        free(el);
        free(idx);
        free(coef);
    } // integer sol

    free(xs[0]); free(xs);
}

 int newTSPSol(void *cbcModel, double obj, int nz, char **vnames, double *x, void *appData) {
     /*printf("Found TSP Solution with Cost: %g\n", obj);
     
     int col = 0;
     for ( int i=0 ; (i<nz) ; ++i ) {
         printf("\t%s %g", vnames[i], x[i]);
         if (++col % 5)
             printf("\t");
         else 
             printf("\n");
     } 
     printf("\n");*/
     return 0;  // FIXME what is the correct return value?
 }

void testTSPUlysses22( char lazy ) {
    if (lazy) {
        printf("TSP Test with 22 cities, with lazy constraints\n");
        printf("==============================================\n");
    }
    else {
        printf("TSP Test with 22 cities, without lazy constraints\n");
        printf("=================================================\n");
    }
    const int n = 22;
    double coord[22][2] = {
         {38.24, 20.42}, {39.57, 26.15}, {40.56, 25.32}, {36.26, 23.12},
         {33.48, 10.54}, {37.56, 12.19}, {38.42, 13.11}, {37.52, 20.44},
         {41.23,  9.10}, {41.17, 13.05}, {36.08, -5.21}, {38.47, 15.13},
         {38.15, 15.35}, {37.51, 15.17}, {35.49, 14.32}, {39.36, 19.56},
         {38.09, 24.36}, {36.09, 23.00}, {40.44, 13.57}, {40.33, 14.15},
         {40.37, 14.23}, {37.57, 22.56}
    };

    double c[22][22];
    for ( int i=0 ; i<22 ; ++i ) {
        for ( int j=0 ; j<22 ; ++j ) {
            if (i==0)
                c[i][j] = 0;
            else
                c[i][j] = dist(coord[i][0], coord[i][1], coord[j][0], coord[j][1]);
        }
    }

    int idx[22];
    double coef[22] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    Cbc_Model *m = Cbc_newModel();
    Cbc_storeNameIndexes(m, 1);
                              
    int x[22][22];
    for ( int i=0 ; (i<22) ; ++i ) {
        for ( int j=0 ; (j<22) ; ++j ) {
            x[i][j] = Cbc_getNumCols(m);
            char vname[256];
            sprintf(vname, "x(%d,%d)", i, j);
            Cbc_addCol(m, vname, 0.0, 1.0, c[i][j], 1, 0, NULL, NULL);
        }
    }

    // out degree
    for ( int i=0 ; (i<n) ; ++i ) {
        int nz = 0;
        for ( int j=0 ; (j<n) ; ++j )  {
            if (i==j)
                continue;
            idx[nz++] = x[i][j];
        }
        char rname[256];
        sprintf(rname, "dout(%d)", i);
        Cbc_addRow(m, rname, nz, idx, coef, 'E', 1.0);
    }

    // in degree
    for ( int i=0 ; (i<n) ; ++i ) {
        int nz = 0;
        for ( int j=0 ; (j<n) ; ++j )  {
            if (i==j)
                continue;
            idx[nz++] = x[j][i];
        }
        char rname[256];
        sprintf(rname, "din(%d)", i);
        Cbc_addRow(m, rname, nz, idx, coef, 'E', 1.0);
    }

    // subtours of size 2
    for ( int i=0 ; (i<n) ; ++i ) {
        for ( int j=i+1 ; (j<n) ; ++j )  {
            char rname[256];
            sprintf(rname, "no2sub(%d)", i);
            idx[0] = x[i][j];
            idx[1] = x[j][i];
            Cbc_addRow(m, rname, 2, idx, coef, 'L', 1.0);
        }
    }

    int y[22];

    if (!lazy) {
        // y vars
        for ( int i=0 ; (i<n) ; ++i ) {
            char vname[256];
            sprintf(vname, "y(%d)", i);
            y[i] = Cbc_getNumCols(m);
            Cbc_addCol(m, vname, 0.0, DBL_MAX, 0, 1, 0, NULL, NULL);
        }

        // weak sub-tour elimination constraints
        coef[1] = -(n+1);
        coef[2] = -1.0;
        for ( int i=1 ; (i<n) ; ++i ) {
            for ( int j=1 ; (j<n) ; ++j )  {
                if (i==j)
                    continue;
                char rname[256];
                sprintf(rname, "noSub(%d,%d)", i, j);
                idx[0] = y[i];
                idx[1] = x[i][j];
                idx[2] = y[j];

                Cbc_addRow(m, rname, 3, idx, coef, 'G', -n);
            }
        }
    }

    
    if (lazy)
        Cbc_addCutCallback(m, subTourSep, "Sub-tour separator", (void *) &n, 1, 1);
    
    Cbc_addIncCallback(m, newTSPSol, NULL);
    
    // adding initial solution
    /*char msColNames[][64] = { "x(0,10)", "x(1,16)", "x(2,1)", "x(3,21)", "x(4,14)", "x(5,4)",
                              "x(6,5)", "x(7,0)", "x(8,9)", "x(9,18)", "x(10,8)", "x(11,15)",
                              "x(12,11)", "x(13,12)", "x(14,13)", "x(15,7)", "x(16,2)", "x(17,3)",
                              "x(18,20)", "x(19,6)", "x(20,19)", "x(21,17)"}; */

    char msColNames[][64] = { "x(0,10)", "x(1,16)", "x(2,1)", "x(3,17)", "x(4,13)", "x(5,14)", "x(6,5)",
                              "x(7,0)", "x(8,9)", "x(9,18)", "x(10,8)", "x(11,15)", "x(12,11)", "x(13,12)",
                              "x(14,4)", "x(15,2)", "x(16,21)", "x(17,7)", "x(18,20)", "x(19,6)", 
                              "x(20,19)", "x(21,3)" };
    char **mscn = malloc(sizeof(char*)*n);
    for ( int i=0 ; (i<n) ; ++i )
        mscn[i] = msColNames[i];
    double msColValues[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                            1.0, 1.0 };
                            
    if (!lazy)
        Cbc_setMIPStart(m, n, (const char **) mscn, msColValues);
    free(mscn);
    Cbc_setMaximumNodes(m, 5000);
    Cbc_solve(m);
    assert(!Cbc_isProvenInfeasible(m));
    assert(!Cbc_isContinuousUnbounded(m));
    assert(!Cbc_isAbandoned(m));
    assert(Cbc_getBestPossibleObjValue(m) >= 5216-1e-4);
    double *sol = NULL;
    if (Cbc_isProvenOptimal(m)) {
        printf("Optimal solution with cost %g found.\n", Cbc_getObjValue(m));
        assert( Cbc_numberSavedSolutions(m) >= 1 );
        assert(fabs(Cbc_getObjValue(m)-5423) < 1e-4);
        assert(fabs(Cbc_getBestPossibleObjValue(m) - Cbc_getObjValue(m))<1e-4);
        assert(Cbc_bestSolution(m) != NULL);
        sol = (double *) Cbc_bestSolution(m);
    }
    else {
        if (Cbc_bestSolution(m)) {
            assert( Cbc_numberSavedSolutions(m) >= 1 );
            assert( Cbc_getObjValue(m) >= 5423-1e-4);
            sol = (double *) Cbc_bestSolution(m);
        }
    }

    if (sol) {
        double **xs = malloc( sizeof(double*)*n );
        xs[0] = malloc( sizeof(double)*n*n );
        for ( int i=1 ; (i<n) ; ++i )
            xs[i] = xs[i-1] + n;
        for ( int i=0 ; i<n ; ++i )
            for ( int j=0 ; j<n ; ++j )
                xs[i][j] = sol[x[i][j]];

        double tot_d = 0.0;
        int *el = malloc(sizeof(int)*(n+1));
        int nel = tspRouteStarting((const double **) xs, 0, n, el);
        assert( nel == n );
        printf("route : ");
        el[n] = 0;
        nel++;
        for ( int i=0 ; (i<nel) ; ++i ) {
            if (i) {
                printf(" -> ");
                tot_d += c[el[i-1]][el[i]];
            }
            printf("%d", el[i]);
        }
        printf("\n");
        printf("Total distance: %g\n", tot_d);
        assert(fabs(tot_d - Cbc_getObjValue(m)) < 1e-4);

        free(el);
        free(xs[0]); free(xs);
    }

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
    
    /* TSP test with 7 cities */
    testTSP(0);  /* only the LP */
    testTSP(1);  /* solving as MIP */

    /* TSP test with 22 cities, with and without
     * lazy constraints */
    testTSPUlysses22( 1 );
    testTSPUlysses22( 0 );

    return 0;
}
