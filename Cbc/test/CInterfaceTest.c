/* $Id$ */
/* Copyright (C) 2014, International Business Machines
   Corporation and others.  All Rights Reserved.
   This code is licensed under the terms of the Eclipse Public License (EPL). */

#undef NDEBUG /* force asserts to work */
#include "Cbc_C_Interface.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>


void testKnapsack() {

    Cbc_Model *model = Cbc_newModel();

    /* Simple knapsack problem
       Minimize -5x[1] - 3x[2] - 2x[3] - 7x[4] - 4x[5]
       s.t.      2x[1] + 8x[2] + 4x[3] + 2x[4] + 5x[5] <= 10
       All x binary
       */
    
    CoinBigIndex start[] = {0, 1, 2, 3, 4, 5, 6};
    int rowindex[] = {0, 0, 0, 0, 0};
    double value[] = {2, 8, 4, 2, 5};
    double collb[] = {0,0,0,0,0};
    double colub[] = {1,1,1,1,1};
    double obj[] = {-5, -3, -2, -7, -4};
    double rowlb[] = {-INFINITY};
    double rowub[] = {10};
    char integer[] = {1,1,1,1,1};
    char *information;
    const double *sol;
    int i;

    Cbc_loadProblem(model, 5, 1, start, rowindex, value, collb, colub, obj, rowlb, rowub);

    Cbc_copyInIntegerInformation(model, integer);

    assert(Cbc_getNumCols(model) == 5);
    assert(Cbc_getNumRows(model) == 1);

    information = Cbc_integerInformation(model);
    for (i = 0; i < 5; i++) {
        assert(information[i] == 1);
    }

    assert(Cbc_optimizationDirection(model) == 1);

    Cbc_solve(model);

    assert(Cbc_isProvenOptimal(model));
    assert(fabs( Cbc_objectiveValue(model)- (-16.0) < 1e-6));
    
    sol = Cbc_getColSolution(model);
    
    assert(fabs(sol[0] - 1.0) < 1e-6);
    assert(fabs(sol[1] - 0.0) < 1e-6);
    assert(fabs(sol[2] - 0.0) < 1e-6);
    assert(fabs(sol[3] - 1.0) < 1e-6);
    assert(fabs(sol[4] - 1.0) < 1e-6);

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
    char integer[] = {1};

    Cbc_loadProblem(model, 1, 1, start, rowindex, value, collb, colub, obj, rowlb, rowub);

    Cbc_copyInIntegerInformation(model, integer);

    assert(Cbc_getNumCols(model) == 1);
    assert(Cbc_getNumRows(model) == 1);

    Cbc_solve(model);
    
    assert(!Cbc_isProvenOptimal(model));
    assert(Cbc_isProvenPrimalInfeasible(model));

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
    double obj[] = {1.0};
    char integer[] = {1,0};

    Cbc_loadProblem(model, 2, 2, start, rowindex, value, collb, colub, obj, rowlb, rowub);

    Cbc_copyInIntegerInformation(model, integer);

    Cbc_solve(model);
    
    assert(!Cbc_isProvenOptimal(model));
    assert(!Cbc_isProvenPrimalInfeasible(model));
    assert(Cbc_isProvenDualInfeasible(model));



}


int main() {

    testKnapsack();
    testIntegerInfeasible();
    testIntegerUnbounded();

    return 0;
}
