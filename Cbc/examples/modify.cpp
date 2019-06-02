// $Id$
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*
  This example shows the creation of a model from arrays, solution
  and then changes to objective and adding a row
*/

#include "CbcModel.hpp"

// Using as solver
#include "OsiClpSolverInterface.hpp"

int main(int argc, const char *argv[])
{
  // model is as exmip1.mps from Data/samples
  int numberRows = 5;
  int numberColumns = 8;
  int numberElements = 14;
  // matrix data - column ordered
  CoinBigIndex start[9] = { 0, 2, 4, 6, 8, 10, 11, 12, 14 };
  int length[8] = { 2, 2, 2, 2, 2, 1, 1, 2 };
  int rows[14] = { 0, 4, 0, 1, 1, 2, 0, 3, 0, 4, 2, 3, 0, 4 };
  double elements[14] = { 3, 5.6, 1, 2, 1.1, 1, -2, 2.8, -1, 1, 1, -1.2, -1, 1.9 };
  CoinPackedMatrix matrix(true, numberRows, numberColumns, numberElements, elements, rows, start, length);

  // rim data
  double objective[8] = { 1, 0, 0, 0, 2, 0, 0, -1 };
  double rowLower[5] = { 2.5, -COIN_DBL_MAX, 4, 1.8, 3 };
  double rowUpper[5] = { COIN_DBL_MAX, 2.1, 4, 5, 15 };
  double colLower[8] = { 2.5, 0, 0, 0, 0.5, 0, 0, 0 };
  double colUpper[8] = { COIN_DBL_MAX, 4.1, 1, 1, 4, COIN_DBL_MAX, COIN_DBL_MAX, 4.3 };
  OsiClpSolverInterface solver1;
  // load problem
  solver1.loadProblem(matrix, colLower, colUpper, objective,
    rowLower, rowUpper);
  // mark integer
  solver1.setInteger(2);
  solver1.setInteger(3);

  // Solve
  solver1.initialSolve();

  // Pass data and solver to CbcModel
  CbcModel model(solver1);

  // reduce printout
  model.setLogLevel(1);
  model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  // Do complete search
  model.branchAndBound();
  /* Print solution.  CbcModel clones solver so we
     need to get current copy */

  const double *solution = model.solver()->getColSolution();
  int i;

  for (i = 0; i < numberColumns; i++) {
    double value = solution[i];
    if (fabs(value) > 1.0e-7 && model.solver()->isInteger(i))
      printf("i %d has value %g\n", i, value);
    else if (fabs(value) > 1.0e-7)
      printf("c %d has value %g\n", i, value);
  }

  // Change objective
  solver1.setObjCoeff(0, -100.0);

  // Now model has too much information e.g. best solution
  // simplest is to start again
  // Pass data and solver to CbcModel
  model = CbcModel(solver1);

  // reduce printout
  model.setLogLevel(1);
  model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  // Do complete search
  model.branchAndBound();

  solution = model.solver()->getColSolution();

  for (i = 0; i < numberColumns; i++) {
    double value = solution[i];
    if (fabs(value) > 1.0e-7 && model.solver()->isInteger(i))
      printf("i %d has value %g\n", i, value);
    else if (fabs(value) > 1.0e-7)
      printf("c %d has value %g\n", i, value);
  }

  // Add constraint
  int column[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  double element2[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
  solver1.addRow(8, column, element2, 7.8, COIN_DBL_MAX);

  // Now model has too much information e.g. best solution
  // simplest is to start again
  // Pass data and solver to CbcModel
  model = CbcModel(solver1);

  // reduce printout
  model.setLogLevel(1);
  model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  // Do complete search
  model.branchAndBound();

  solution = model.solver()->getColSolution();

  for (i = 0; i < numberColumns; i++) {
    double value = solution[i];
    if (fabs(value) > 1.0e-7 && model.solver()->isInteger(i))
      printf("i %d has value %g\n", i, value);
    else if (fabs(value) > 1.0e-7)
      printf("c %d has value %g\n", i, value);
  }
  return 0;
}
