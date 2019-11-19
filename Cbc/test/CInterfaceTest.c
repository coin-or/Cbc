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
#include <limits.h>

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
    //assert(!Cbc_isProvenInfeasible(model));
    //assert(Cbc_isContinuousUnbounded(model));

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

    testTSP(0);
    testTSP(1);

    return 0;
}
