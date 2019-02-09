// $Id$
// Copyright (C) 2008, Stefan Vigerske, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif
#include <cassert>
#include <iostream>
using namespace std;
#include "CoinHelperFunctions.hpp"
#include "CoinError.hpp"
#include "CbcSolver.hpp"
#include "CbcModel.hpp"
#include "CbcBranchActual.hpp" //for CbcSOS
#include "CbcBranchLotsize.hpp" //for CbcLotsize
#include "OsiClpSolverInterface.hpp"
#define testtol 1e-6
/** model sos1a from the GAMS test library
 * http://www.gams.com/testlib/libhtml/sos1a.htm */
void sos1a(int &error_count, int &warning_count);
/** model sos2a from the GAMS test library
 * http://www.gams.com/testlib/libhtml/sos2a.htm */
void sos2a(int &error_count, int &warning_count);
/** model semicon1 from the GAMS test library
 * http://www.gams.com/testlib/libhtml/semicon1.htm */
void semicon1(int &error_count, int &warning_count);
/** model semiint1 from the GAMS test library
 * http://www.gams.com/testlib/libhtml/semiint1.htm */
void semiint1(int &error_count, int &warning_count);
int main(int argc, const char *argv[])
{
  WindowsErrorPopupBlocker();
  // only in CoinUtils/trunk: WindowsErrorPopupBlocker();
  int error_count = 0;
  int warning_count = 0;

  sos1a(error_count, warning_count);
  cout << "\n***********************\n"
       << endl;
  sos2a(error_count, warning_count);
  cout << "\n***********************\n"
       << endl;
  semicon1(error_count, warning_count);
  cout << "\n***********************\n"
       << endl;
  semiint1(error_count, warning_count);

  cout << endl
       << "Finished - there have been " << error_count << " errors and " << warning_count << " warnings." << endl;
  return error_count;
}
void sos1a(int &error_count, int &warning_count)
{
  OsiClpSolverInterface solver1;

  int numcols = 3;
  int numrows = 1;
  int nnz = 3;
  CoinBigIndex *start = new int[numcols + 1];
  int *index = new int[nnz];
  double *value = new double[nnz];
  double *collb = new double[numcols];
  double *colub = new double[numcols];
  double *obj = new double[numcols];
  double *rowlb = new double[numrows];
  double *rowub = new double[numrows];
  // objective
  obj[0] = .9;
  obj[1] = 1.;
  obj[2] = 1.1;

  // column bounds
  collb[0] = 0.;
  colub[0] = .8;
  collb[1] = 0.;
  colub[1] = .6;
  collb[2] = 0.;
  colub[2] = .6;
  // matrix
  start[0] = 0;
  index[0] = 0;
  value[0] = 1.;
  start[1] = 1;
  index[1] = 0;
  value[1] = 1.;
  start[2] = 2;
  index[2] = 0;
  value[2] = 1.;
  start[3] = 3;

  // row bounds
  rowlb[0] = -solver1.getInfinity();
  rowub[0] = 1.;
  solver1.loadProblem(numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub);
  solver1.setObjSense(-1);

  CbcSolverUsefulData data;
  CbcModel model(solver1);
  CbcMain0(model, data);

  int which[3] = { 0, 1, 2 };
  CbcObject *sosobject = new CbcSOS(&model, 3, which, NULL, 0, 1);
  model.addObjects(1, &sosobject);
  delete sosobject;

  const char *argv2[] = { "gamstest_sos1a", "-solve", "-quit" };
  CbcMain1(3, argv2, model, NULL, data);
  cout << endl;
  if (!model.isProvenOptimal()) {
    cerr << "Error: Model sos1a not solved to optimality." << endl;
    ++error_count;
    return; // other tests make no sense ---- memory leak here
  }

  OsiSolverInterface *solver = model.solver();
  assert(solver);

  cout << "Objective value model: " << model.getObjValue()
       << "\t solver: " << solver->getObjValue()
       << "\t expected: 0.72" << endl;
  if (CoinAbs(model.getObjValue() - 0.72) > testtol || CoinAbs(solver->getObjValue() - 0.72) > testtol) {
    cerr << "Error: Objective value incorrect." << endl;
    ++error_count;
  }

  cout << "Primal value variable 0 in model: " << model.bestSolution()[0]
       << "\t in solver: " << solver->getColSolution()[0]
       << "\t expected: 0.8" << endl;
  if (CoinAbs(model.bestSolution()[0] - 0.8) > testtol || CoinAbs(solver->getColSolution()[0] - 0.8) > testtol) {
    cerr << "Error: Primal value incorrect." << endl;
    ++error_count;
  }
  cout << "Primal value variable 1 in model: " << model.bestSolution()[1]
       << "\t in solver: " << solver->getColSolution()[1]
       << "\t expected: 0.0" << endl;
  if (CoinAbs(model.bestSolution()[1]) > testtol || CoinAbs(solver->getColSolution()[1]) > testtol) {
    cerr << "Error: Primal value incorrect." << endl;
    ++error_count;
  }
  cout << "Primal value variable 2 in model: " << model.bestSolution()[2]
       << "\t in solver: " << solver->getColSolution()[2]
       << "\t expected: 0.0" << endl;
  if (CoinAbs(model.bestSolution()[2]) > testtol || CoinAbs(solver->getColSolution()[2]) > testtol) {
    cerr << "Error: Primal value incorrect." << endl;
    ++error_count;
  }
  delete[] start;
  delete[] index;
  delete[] value;
  delete[] collb;
  delete[] colub;
  delete[] obj;
  delete[] rowlb;
  delete[] rowub;
}
void sos2a(int &error_count, int &warning_count)
{
  OsiClpSolverInterface solver1;

  int numcols = 7; // w1, w2, w3, x, fx, fplus, fminus
  int numrows = 5; // wsum, xdef, fxdef, gapplus, gapminus
  int nnz = 15;
  CoinBigIndex *start = new int[numcols + 1];
  int *index = new int[nnz];
  double *value = new double[nnz];
  double *collb = new double[numcols];
  double *colub = new double[numcols];
  double *obj = new double[numcols];
  double *rowlb = new double[numrows];
  double *rowub = new double[numrows];
  // objective
  obj[0] = 0.;
  obj[1] = 0.;
  obj[2] = 0.;
  obj[3] = 0.;
  obj[4] = 0.;
  obj[5] = 1.;
  obj[6] = 1.;

  // column bounds
  collb[0] = 0.;
  colub[0] = solver1.getInfinity();
  collb[1] = 0.;
  colub[1] = solver1.getInfinity();
  collb[2] = 0.;
  colub[2] = solver1.getInfinity();
  collb[3] = -solver1.getInfinity();
  colub[3] = solver1.getInfinity();
  collb[4] = -solver1.getInfinity();
  colub[4] = solver1.getInfinity();
  collb[5] = 0.;
  colub[5] = solver1.getInfinity();
  collb[6] = 0.;
  colub[6] = solver1.getInfinity();
  // matrix
  start[0] = 0;
  index[0] = 0;
  value[0] = 1.;
  index[1] = 1;
  value[1] = 1.;
  index[2] = 2;
  value[2] = 1.;
  start[1] = 3;
  index[3] = 0;
  value[3] = 1.;
  index[4] = 1;
  value[4] = 2.;
  index[5] = 2;
  value[5] = 2.;
  start[2] = 6;
  index[6] = 0;
  value[6] = 1.;
  index[7] = 1;
  value[7] = 3.;
  index[8] = 2;
  value[8] = 3.;
  start[3] = 9;
  index[9] = 1;
  value[9] = -1.;
  start[4] = 10;
  index[10] = 2;
  value[10] = -1.;
  index[11] = 3;
  value[11] = -1.;
  index[12] = 4;
  value[12] = 1.;
  start[5] = 13;
  index[13] = 3;
  value[13] = 1.;
  start[6] = 14;
  index[14] = 4;
  value[14] = 1.;
  start[7] = 15;

  // row bounds
  rowlb[0] = 1.;
  rowub[0] = 1.;
  rowlb[1] = 0.;
  rowub[1] = 0.;
  rowlb[2] = 0.;
  rowub[2] = 0.;
  rowlb[3] = -1.3;
  rowub[3] = solver1.getInfinity();
  rowlb[4] = 1.3;
  rowub[4] = solver1.getInfinity();
  solver1.loadProblem(numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub);
  double *primalval = new double[numcols];
  double *redcost = new double[numcols];
  double optvalue = solver1.getInfinity();
  for (int testcase = 0; testcase < 2; ++testcase) {
    switch (testcase) {
    case 0:
      solver1.setColLower(0, 0.);
      optvalue = 0.;
      primalval[0] = .7;
      redcost[0] = 0.;
      primalval[1] = .3;
      redcost[1] = 0.;
      primalval[2] = 0.;
      redcost[2] = 0.;
      primalval[3] = 1.3;
      redcost[3] = 0.;
      primalval[4] = 1.3;
      redcost[4] = 0.;
      primalval[5] = 0.;
      redcost[5] = 1.;
      primalval[6] = 0.;
      redcost[6] = 1.;
      break;
    case 1:
      solver1.setColLower(0, .8);
      optvalue = 0.1;
      primalval[0] = .8;
      redcost[0] = 1.;
      primalval[1] = .2;
      redcost[1] = 0.;
      primalval[2] = 0.;
      redcost[2] = -1.;
      primalval[3] = 1.2;
      redcost[3] = 0.;
      primalval[4] = 1.2;
      redcost[4] = 0.;
      primalval[5] = 0.;
      redcost[5] = 1.;
      primalval[6] = 0.1;
      redcost[6] = 0.;
      break;
    }
    CbcModel model(solver1);
    CbcSolverUsefulData data;
    CbcMain0(model, data);
    int which[3] = { 0, 1, 2 };
    CbcObject *sosobject = new CbcSOS(&model, 3, which, NULL, 0, 2);
    model.addObjects(1, &sosobject);
    delete sosobject;
    const char *argv2[] = { "gamstest_sos2a", "-solve", "-quit" };
    cout << "\nSolving sos2a model with w1 having lower bound " << solver1.getColLower()[0] << endl;
    CbcMain1(3, argv2, model, NULL, data);
    cout << endl;
    if (!model.isProvenOptimal()) {
      cerr << "Error: Model sos2a not solved to optimality." << endl;
      ++error_count;
      continue; // other tests make no sense
    }
    OsiSolverInterface *solver = model.solver();
    assert(solver);
    cout << "Objective value model: " << model.getObjValue()
         << "\t solver: " << solver->getObjValue()
         << "\t expected: " << optvalue << endl;
    if (CoinAbs(model.getObjValue() - optvalue) > testtol || CoinAbs(solver->getObjValue() - optvalue) > testtol) {
      cerr << "Error: Objective value incorrect." << endl;
      ++error_count;
    }
    for (int i = 0; i < numcols; ++i) {
      cout << "Primal value variable " << i << " in model: " << model.bestSolution()[i]
           << "\t in solver: " << solver->getColSolution()[i]
           << "\t expected: " << primalval[i]
           << endl;
      if (CoinAbs(model.bestSolution()[i] - primalval[i]) > testtol || CoinAbs(solver->getColSolution()[i] - primalval[i]) > testtol) {
        cerr << "Error: Primal value incorrect." << endl;
        ++error_count;
      }
    }
    for (int i = 0; i < numcols; ++i) {
      cout << "Reduced cost variable " << i << " in model: " << model.getReducedCost()[i]
           << "\t in solver: " << solver->getReducedCost()[i]
           << "\t expected: " << redcost[i]
           << endl;
      if (CoinAbs(model.getReducedCost()[i] - redcost[i]) > testtol || CoinAbs(solver->getReducedCost()[i] - redcost[i]) > testtol) {
        cerr << "Warning: Reduced cost incorrect." << endl;
        ++warning_count;
      }
    }
  }
  delete[] start;
  delete[] index;
  delete[] value;
  delete[] collb;
  delete[] colub;
  delete[] obj;
  delete[] rowlb;
  delete[] rowub;
  delete[] primalval;
  delete[] redcost;
}
void semicon1(int &error_count, int &warning_count)
{
  OsiClpSolverInterface solver1;

  int numcols = 4; // s, pup, plo, x
  int numrows = 3; // bigx, smallx, f
  int nnz = 6;
  CoinBigIndex *start = new int[numcols + 1];
  int *index = new int[nnz];
  double *value = new double[nnz];
  double *collb = new double[numcols];
  double *colub = new double[numcols];
  double *obj = new double[numcols];
  double *rowlb = new double[numrows];
  double *rowub = new double[numrows];
  // objective
  obj[0] = 0;
  obj[1] = 1.;
  obj[2] = 1;
  obj[3] = 0;

  // column bounds
  collb[0] = 0.;
  colub[0] = 10.;
  collb[1] = 0.;
  colub[1] = solver1.getInfinity();
  collb[2] = 0.;
  colub[2] = solver1.getInfinity();
  collb[3] = 0.;
  colub[3] = solver1.getInfinity();
  // matrix
  start[0] = 0;
  index[0] = 2;
  value[0] = 1.;
  start[1] = 1;
  index[1] = 0;
  value[1] = -1.;
  start[2] = 2;
  index[2] = 1;
  value[2] = 1.;
  start[3] = 3;
  index[3] = 0;
  value[3] = 1.;
  index[4] = 1;
  value[4] = 1.;
  index[5] = 2;
  value[5] = 1.;
  start[4] = nnz;

  // row bounds
  rowlb[0] = -solver1.getInfinity();
  rowub[0] = 8.9;
  rowlb[1] = 8.9;
  rowub[1] = solver1.getInfinity();
  rowlb[2] = 10.;
  rowub[2] = 10.;
  solver1.loadProblem(numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub);
  for (int testcase = 0; testcase < 5; ++testcase) {
    CbcModel model(solver1);
    CbcSolverUsefulData data;
    CbcMain0(model, data);

    double points[4] = { 0., 0., 0., 10. };
    double objval;
    double primalval[4];
    double redcost[4];
    double row2marg;
    redcost[1] = 1.0;
    redcost[2] = 1.0;
    redcost[3] = 0.0;
    switch (testcase) {
    case 0:
      points[2] = 0.;
      objval = 0.;
      primalval[0] = 1.1;
      primalval[1] = 0.0;
      primalval[2] = 0.0;
      primalval[3] = 8.9;
      redcost[0] = 0.0;
      row2marg = 0.0;
      break;
    case 1:
      points[2] = 1.;
      objval = 0.;
      primalval[0] = 1.1;
      primalval[1] = 0.0;
      primalval[2] = 0.0;
      primalval[3] = 8.9;
      redcost[0] = 0.0;
      row2marg = 0.0;
      break;
    case 2:
      points[2] = 1.5;
      objval = 0.4;
      primalval[0] = 1.5;
      primalval[1] = 0.0;
      primalval[2] = 0.4;
      primalval[3] = 8.5;
      redcost[0] = 1.0;
      row2marg = -1.0;
      break;
    case 3:
      points[2] = 2.1;
      objval = 1.0;
      primalval[0] = 2.1;
      primalval[1] = 0.0;
      primalval[2] = 1.0;
      primalval[3] = 7.9;
      redcost[0] = 1.0;
      row2marg = -1.0;
      break;
    case 4:
      points[2] = 2.8;
      objval = 1.1;
      primalval[0] = 0.0;
      primalval[1] = 1.1;
      primalval[2] = 0.0;
      primalval[3] = 10.0;
      redcost[0] = -1.0;
      row2marg = 1.0;
      break;
    default: // to please the compile
      redcost[0] = 0.;
      row2marg = 0.;
      objval = 0.;
    }

    CbcObject *semiconobject = new CbcLotsize(&model, 0, 2, points, true);
    model.addObjects(1, &semiconobject);
    delete semiconobject;

    cout << "\nSolving semicon1 model for lotsize variable being either 0 or between " << points[2] << " and 10.\n"
         << endl;
    const char *argv2[] = { "gamstest_semicon1", "-solve", "-quit" };
    CbcMain1(3, argv2, model, NULL, data);
    cout << endl;
    if (!model.isProvenOptimal()) {
      cerr << "Error: Model semicon1 not solved to optimality." << endl;
      ++error_count;
      continue; // other tests make no sense
    }
    OsiSolverInterface *solver = model.solver();
    assert(solver);
    cout << "Objective value in model: " << model.getObjValue()
         << "\t in solver: " << solver->getObjValue()
         << "\t expected: " << objval << endl;
    if (CoinAbs(model.getObjValue() - objval) > testtol || CoinAbs(solver->getObjValue() - objval) > testtol) {
      cerr << "Error: Objective value incorrect." << endl;
      ++error_count;
    }
    for (int i = 0; i < numcols; ++i) {
      cout << "Primal value variable " << i << " in model: " << model.bestSolution()[i]
           << "\t in solver: " << solver->getColSolution()[i]
           << "\t expected: " << primalval[i]
           << endl;
      if (CoinAbs(model.bestSolution()[i] - primalval[i]) > testtol || CoinAbs(solver->getColSolution()[i] - primalval[i]) > testtol) {
        cerr << "Error: Primal value incorrect." << endl;
        ++error_count;
      }
    }
    cout << "Reduced cost variable " << 0 << " in model: " << model.getReducedCost()[0]
         << "\t in solver: " << solver->getReducedCost()[0]
         << "\t expected: " << redcost[0]
         << endl;
    if (CoinAbs(model.getReducedCost()[0] - redcost[0]) > testtol || CoinAbs(solver->getReducedCost()[0] - redcost[0]) > testtol) {
      cerr << "Warning: Reduced cost incorrect." << endl;
      ++warning_count;
    }
    cout << "Reduced cost variable " << 3 << " in model: " << model.getReducedCost()[3]
         << "\t in solver: " << solver->getReducedCost()[3]
         << "\t expected: " << redcost[3]
         << endl;
    if (CoinAbs(model.getReducedCost()[3] - redcost[3]) > testtol || CoinAbs(solver->getReducedCost()[3] - redcost[3]) > testtol) {
      cerr << "Warning: Reduced cost incorrect." << endl;
      ++warning_count;
    }
    cout << "Reduced cost variable 1 plus - dual of row 0 in model: " << model.getReducedCost()[1] - model.getRowPrice()[0]
         << "\t expected: " << redcost[1]
         << endl;
    if (CoinAbs(model.getReducedCost()[1] - model.getRowPrice()[0] - redcost[1]) > testtol) {
      cerr << "Warning: Reduced cost or row margin incorrect." << endl;
      ++warning_count;
    }
    cout << "Reduced cost variable 2 plus + dual of row 1 in model: " << model.getReducedCost()[2] + model.getRowPrice()[1]
         << "\t expected: " << redcost[2]
         << endl;
    if (CoinAbs(model.getReducedCost()[2] + model.getRowPrice()[1] - redcost[2]) > testtol) {
      cerr << "Warning: Reduced cost or row margin incorrect." << endl;
      ++warning_count;
    }

    cout << "Row 2 marginal (price) in model: " << model.getRowPrice()[2]
         << "\t in solver: " << solver->getRowPrice()[2]
         << "\t expected: " << row2marg << endl;
    if (CoinAbs(model.getRowPrice()[2] - row2marg) > testtol || CoinAbs(solver->getRowPrice()[2] - row2marg) > testtol) {
      cerr << "Warning: Row price incorrect." << endl;
      ++warning_count;
    }
  }

  delete[] start;
  delete[] index;
  delete[] value;
  delete[] collb;
  delete[] colub;
  delete[] obj;
  delete[] rowlb;
  delete[] rowub;
}
void semiint1(int &error_count, int &warning_count)
{
  OsiClpSolverInterface solver1;

  int numcols = 4; // s, pup, plo, x
  int numrows = 3; // bigx, smallx, f
  int nnz = 6;
  CoinBigIndex *start = new int[numcols + 1];
  int *index = new int[nnz];
  double *value = new double[nnz];
  double *collb = new double[numcols];
  double *colub = new double[numcols];
  double *obj = new double[numcols];
  double *rowlb = new double[numrows];
  double *rowub = new double[numrows];
  // objective
  obj[0] = 0;
  obj[1] = 1.;
  obj[2] = 1;
  obj[3] = 0;

  // column bounds
  collb[0] = 0.;
  colub[0] = 10.;
  collb[1] = 0.;
  colub[1] = solver1.getInfinity();
  collb[2] = 0.;
  colub[2] = solver1.getInfinity();
  collb[3] = 0.;
  colub[3] = solver1.getInfinity();
  // matrix
  start[0] = 0;
  index[0] = 2;
  value[0] = 1.;
  start[1] = 1;
  index[1] = 0;
  value[1] = -1.;
  start[2] = 2;
  index[2] = 1;
  value[2] = 1.;
  start[3] = 3;
  index[3] = 0;
  value[3] = 1.;
  index[4] = 1;
  value[4] = 1.;
  index[5] = 2;
  value[5] = 1.;
  start[4] = nnz;

  // row bounds
  rowlb[0] = -solver1.getInfinity();
  rowub[0] = 7.9;
  rowlb[1] = 7.9;
  rowub[1] = solver1.getInfinity();
  rowlb[2] = 10.;
  rowub[2] = 10.;
  solver1.loadProblem(numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub);
  solver1.setInteger(0);

  for (int testcase = 0; testcase < 6; ++testcase) {
    CbcModel model(solver1);
    CbcSolverUsefulData data;
    CbcMain0(model, data);

    double points[10];
    points[0] = 0.;
    int nrpoints = 0;
    double objval;
    double primalval[4];
    double redcost[4];
    double row2marg;
    redcost[2] = 1.0;
    redcost[3] = 0.0;
    switch (testcase) {
    case 0:
      nrpoints = 0; // pure integer case
      objval = 0.1;
      primalval[0] = 2.0;
      primalval[1] = 0.1;
      primalval[2] = 0.0;
      primalval[3] = 8;
      redcost[0] = -1.0;
      redcost[1] = 0.0;
      row2marg = 1.0;
      break;
    case 1:
      nrpoints = 0; // pure integer case too
      objval = 0.1;
      primalval[0] = 2.0;
      primalval[1] = 0.1;
      primalval[2] = 0.0;
      primalval[3] = 8.0;
      redcost[0] = -1.0;
      redcost[1] = 0.0;
      row2marg = 1.0;
      break;
    case 2:
      for (nrpoints = 1; nrpoints < 10; ++nrpoints)
        points[nrpoints] = nrpoints + 1;
      objval = 0.1;
      primalval[0] = 2.0;
      primalval[1] = 0.1;
      primalval[2] = 0.0;
      primalval[3] = 8.0;
      redcost[0] = -1.0;
      redcost[1] = 0.0;
      row2marg = 1.0;
      break;
    case 3:
      for (nrpoints = 1; nrpoints < 9; ++nrpoints)
        points[nrpoints] = nrpoints + 2;
      objval = 0.9;
      primalval[0] = 3.0;
      primalval[1] = 0.0;
      primalval[2] = 0.9;
      primalval[3] = 7.0;
      redcost[0] = 1.0;
      redcost[1] = 1.0;
      row2marg = -1.0;
      break;
    case 4:
      for (nrpoints = 1; nrpoints < 8; ++nrpoints)
        points[nrpoints] = nrpoints + 3;
      objval = 1.9;
      primalval[0] = 4.0;
      primalval[1] = 0.0;
      primalval[2] = 1.9;
      primalval[3] = 6.0;
      redcost[0] = 1.0;
      redcost[1] = 1.0;
      row2marg = -1.0;
      break;
    case 5:
      for (nrpoints = 1; nrpoints < 7; ++nrpoints)
        points[nrpoints] = nrpoints + 4;
      objval = 2.1;
      primalval[0] = 0.0;
      primalval[1] = 2.1;
      primalval[2] = 0.0;
      primalval[3] = 10.0;
      redcost[0] = -1.0;
      redcost[1] = 0.0;
      row2marg = 1.0;
      break;
    default: // to please the compile
      redcost[0] = 0.;
      redcost[1] = 0.;
      row2marg = 0.;
      objval = 0.;
    }
    if (nrpoints) {
      CbcObject *semiintobject = new CbcLotsize(&model, 0, nrpoints, points);
      model.addObjects(1, &semiintobject);
      delete semiintobject;
    }

    cout << "\nSolving semiint1 model for integer lotsize variable being either 0 or between " << points[2] << " and 10.\n"
         << endl;
    const char *argv2[] = { "gamstest_semiint1", "-solve", "-quit" };
    CbcMain1(3, argv2, model, NULL, data);
    cout << endl;
    if (!model.isProvenOptimal()) {
      cerr << "Error: Model semiint1 not solved to optimality." << endl;
      ++error_count;
      continue; // other tests make no sense
    }
    OsiSolverInterface *solver = model.solver();
    assert(solver);
    cout << "Objective value in model: " << model.getObjValue()
         << "\t in solver: " << solver->getObjValue()
         << "\t expected: " << objval << endl;
    if (CoinAbs(model.getObjValue() - objval) > testtol || CoinAbs(solver->getObjValue() - objval) > testtol) {
      cerr << "Error: Objective value incorrect." << endl;
      ++error_count;
    }
    for (int i = 0; i < numcols; ++i) {
      cout << "Primal value variable " << i << " in model: " << model.bestSolution()[i]
           << "\t in solver: " << solver->getColSolution()[i]
           << "\t expected: " << primalval[i]
           << endl;
      if (CoinAbs(model.bestSolution()[i] - primalval[i]) > testtol || CoinAbs(solver->getColSolution()[i] - primalval[i]) > testtol) {
        cerr << "Error: Primal value incorrect." << endl;
        ++error_count;
      }
    }
    cout << "Reduced cost variable " << 0 << " in model: " << model.getReducedCost()[0]
         << "\t in solver: " << solver->getReducedCost()[0]
         << "\t expected: " << redcost[0]
         << endl;
    if (CoinAbs(model.getReducedCost()[0] - redcost[0]) > testtol || CoinAbs(solver->getReducedCost()[0] - redcost[0]) > testtol) {
      cerr << "Warning: Reduced cost incorrect." << endl;
      ++warning_count;
    }
    cout << "Reduced cost variable " << 3 << " in model: " << model.getReducedCost()[3]
         << "\t in solver: " << solver->getReducedCost()[3]
         << "\t expected: " << redcost[3]
         << endl;
    if (CoinAbs(model.getReducedCost()[3] - redcost[3]) > testtol || CoinAbs(solver->getReducedCost()[3] - redcost[3]) > testtol) {
      cerr << "Warning: Reduced cost incorrect." << endl;
      ++warning_count;
    }
    cout << "Row 2 marginal (price) in model: " << model.getRowPrice()[2]
         << "\t in solver: " << solver->getRowPrice()[2]
         << "\t expected: " << row2marg << endl;
    if (CoinAbs(model.getRowPrice()[2] - row2marg) > testtol || CoinAbs(solver->getRowPrice()[2] - row2marg) > testtol) {
      cerr << "Warning: Row price incorrect." << endl;
      ++warning_count;
    }

    cout << "Row 2 marginal (price) in model: " << model.getRowPrice()[2]
         << "\t in solver: " << solver->getRowPrice()[2]
         << "\t expected: " << row2marg << endl;
    if (CoinAbs(model.getRowPrice()[2] - row2marg) > testtol || CoinAbs(solver->getRowPrice()[2] - row2marg) > testtol) {
      cerr << "Warning: Row price incorrect." << endl;
      ++warning_count;
    }
  }

  delete[] start;
  delete[] index;
  delete[] value;
  delete[] collb;
  delete[] colub;
  delete[] obj;
  delete[] rowlb;
  delete[] rowub;
}
