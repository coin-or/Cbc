// $Id$
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
// For Branch and bound
#include "CbcModel.hpp"
#include "CbcBranchLotsize.hpp"
#include "OsiClpSolverInterface.hpp"

// Time
#include "CoinTime.hpp"

/************************************************************************

This main program reads in an integer model from an mps file.

It then replaces all 0-1 variables by lotsizing variables
which can take values 0.0,0.45-0.55 or 1.0

*************************************************************************/
int main(int argc, const char *argv[])
{

  // Define your favorite OsiSolver

  OsiClpSolverInterface solver1;

  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName;
#if defined(MIPLIB3DIR)
  mpsFileName = MIPLIB3DIR "/10teams";
#else
  if (argc < 2) {
    fprintf(stderr, "Do not know where to find miplib3 MPS files.\n");
    exit(1);
  }
#endif
  if (argc >= 2)
    mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(), "");
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }

  int iColumn;
  int numberColumns = solver1.getNumCols();
  int numberLot = 0;
  char *mark = new char[numberColumns];
  // take off integers but find where they are
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (solver1.isBinary(iColumn)) {
      solver1.setContinuous(iColumn);
      mark[iColumn] = 1;
      numberLot++;
    } else {
      mark[iColumn] = 0;
    }
  }
  CbcModel model(solver1);
  // Do lotsizing
  CbcObject **objects = new CbcObject *[numberLot];
  numberLot = 0;
  /* For semi-continuous variables numberRanges is 2
     and ranges[]={0.0,0.0,K,COIN_DBL_MAX};
  */
  // valid ranges are 0.0 to 0.0, 0.45 to 0.55, 1.0 to 1.0
  double ranges[] = { 0.0, 0.0, 0.45, 0.55, 1.0, 1.0 };
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (mark[iColumn])
      objects[numberLot++] = new CbcLotsize(&model, iColumn, 3, ranges, true);
  }
  delete[] mark;
  model.addObjects(numberLot, objects);
  for (iColumn = 0; iColumn < numberLot; iColumn++)
    delete objects[iColumn];
  delete[] objects;

  // If time is given then stop after that number of minutes
  if (argc > 2) {
    int minutes = atoi(argv[2]);
    std::cout << "Stopping after " << minutes << " minutes" << std::endl;
    assert(minutes >= 0);
    model.setDblParam(CbcModel::CbcMaximumSeconds, 60.0 * minutes);
  }
  // Switch off most output
  model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  if (model.getNumCols() < 3000) {
    model.messageHandler()->setLogLevel(1);
    //model.solver()->messageHandler()->setLogLevel(0);
  } else {
    model.messageHandler()->setLogLevel(2);
    model.solver()->messageHandler()->setLogLevel(1);
  }
  model.messageHandler()->setLogLevel(1);

  double time1 = CoinCpuTime();

  // Do complete search

  model.branchAndBound();

  std::cout << mpsFileName << " took " << CoinCpuTime() - time1 << " seconds, "
            << model.getNodeCount() << " nodes with objective "
            << model.getObjValue()
            << (!model.status() ? " Finished" : " Not finished")
            << std::endl;

  // Print solution - we can't get names from Osi!

  if (model.getMinimizationObjValue() < 1.0e50) {
    int numberColumns = model.solver()->getNumCols();

    const double *solution = model.solver()->getColSolution();

    int iColumn;
    std::cout << std::setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14);

    std::cout << "--------------------------------------" << std::endl;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      if (fabs(value) > 1.0e-7)
        std::cout << std::setw(6) << iColumn << " " << value << std::endl;
    }
    std::cout << "--------------------------------------" << std::endl;

    std::cout << std::resetiosflags(std::ios::fixed | std::ios::showpoint | std::ios::scientific);
  }
  return 0;
}
