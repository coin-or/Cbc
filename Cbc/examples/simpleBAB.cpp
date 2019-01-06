// $Id$
// Copyright (C) 2009, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinPackedVector.hpp"
//#define USE_CBC
#ifdef USE_CBC
#include "CbcModel.hpp"
#endif

int main(int argc, const char *argv[])
{

  OsiClpSolverInterface model;

  CoinBigIndex start[] = { 0, 1, 2 };
  int index[] = { 0, 0 };
  double values[] = { 1.0, 2.0 };
  double collb[] = { 0.0, 0.0 };
  double colub[] = { 10.0, 10.0 };
  double obj[] = { 1.0, 1.0 };
  double rowlb[] = { 0.0 };
  double rowub[] = { 3.9 };

  // obj: Max x0 + x1
  //  st. x0 + 2 x1 <= 3.9
  //          0 <= x0 <= 10 and integer
  //          0 <= x1 <= 10
  model.loadProblem(2, 1, start, index, values, collb, colub, obj, rowlb, rowub);
  model.setInteger(0);
  model.setObjSense(-1.0);
  //bool optimal;

#ifndef USE_CBC
  // Save bounds - that and dual limit should be all that is needed
  // For this simple example we could just re-use collb and colub
  double saveLower[2];
  double saveUpper[2];
  int numberColumns = model.getNumCols();
  CoinCopyN(model.getColLower(), numberColumns, saveLower);
  CoinCopyN(model.getColUpper(), numberColumns, saveUpper);
  double objLimit;
  model.getDblParam(OsiDualObjectiveLimit, objLimit);
  model.branchAndBound();
  //optimal = model.isProvenOptimal();
  const double *val = model.getColSolution(); // x0 = 3, x1 = 0.45
  printf("Solution %g %g\n", val[0], val[1]);
  // Restore bounds and dual limit
  model.setColLower(saveLower);
  model.setColUpper(saveUpper);
  model.setDblParam(OsiDualObjectiveLimit, objLimit);
#else
  {
    CbcModel model2(model);
    model2.branchAndBound();
    //optimal = model2.isProvenOptimal();
    const double *val = model2.getColSolution(); // x0 = 3, x1 = 0.45
    printf("Solution %g %g\n", val[0], val[1]);
  }
#endif

  const int rowCols[] = { 0 };
  const double rowElements = { 1.0 };

  // add x0 <= 2, and solve once again.
  CoinPackedVector v(1, rowCols, rowElements);
  model.addRow(v, 0.0, 2.0);
#ifndef USE_CBC
  model.branchAndBound();
  //optimal = model.isProvenOptimal(); // should be x0 = 2, x1 = 0.95
  // Address of solution will be same as only adding rows - but be safe
  val = model.getColSolution();
  printf("Solution %g %g\n", val[0], val[1]);
#else
  {
    CbcModel model2(model);
    model2.branchAndBound();
    //optimal = model2.isProvenOptimal(); // should be x0 = 2, x1 = 0.95
    const double *val = model2.getColSolution();
    printf("Solution %g %g\n", val[0], val[1]);
  }
#endif
  return 0;
}
