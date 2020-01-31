// $Id$
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcModel.hpp"

// Using as solver
#include "OsiClpSolverInterface.hpp"

int main(int argc, const char *argv[])
{
  OsiClpSolverInterface solver1;
  // Read in example model
  // and assert that it is a clean model
#if defined(SAMPLEDIR)
  int numMpsReadErrors = solver1.readMps(SAMPLEDIR "/p0033.mps", "");
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }
#else
  fprintf(stderr, "Do not know where to find sample MPS files.\n");
  exit(1);
#endif

  // Pass data and solver to CbcModel
  CbcModel model(solver1);

  // uncomment to reduce printout
  //model.setLogLevel(1);
  //model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  // Do complete search
  model.branchAndBound();
  /* Print solution.  CbcModel clones solver so we
     need to get current copy */
  int numberColumns = model.solver()->getNumCols();

  const double *solution = model.solver()->getColSolution();

  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    double value = solution[iColumn];
    if (fabs(value) > 1.0e-7 && model.solver()->isInteger(iColumn))
      printf("%d has value %g\n", iColumn, value);
  }
  return 0;
}
