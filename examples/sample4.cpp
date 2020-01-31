// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcConfig.h"

#include "CoinPragma.hpp"

#include <cassert>
#include <iomanip>

// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcBranchUser.hpp"
#include "CbcBranchActual.hpp"
#include "CbcCompareUser.hpp"
#include "CoinTime.hpp"
#include "OsiClpSolverInterface.hpp"

//#############################################################################

/************************************************************************

This main program reads in an SOS model (rgn) from an mps file.

It then solves it three ways :-

a) As normal
b) SOS 1
c) SOS 2 (so answer will be different)

************************************************************************/

int main(int argc, const char *argv[])
{

  // Define your favorite OsiSolver

  OsiClpSolverInterface solver1;
  //solver1.messageHandler()->setLogLevel(0);
  CbcModel model(solver1);
  model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);

  // Read in rgn.mps
  std::string mpsFileName;
#if defined(MIPLIB3DIR)
  mpsFileName = MIPLIB3DIR "/rgn";
#else
  if (argc < 2) {
    fprintf(stderr, "Do not know where to find miplib3 MPS files.\n");
    exit(1);
  }
#endif
  if (argc >= 2)
    mpsFileName = argv[1];
  int numMpsReadErrors = model.solver()->readMps(mpsFileName.c_str(), "");
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }

  // Definition of node choice
  CbcCompareUser compare;
  compare.setWeight(0.0);
  model.setNodeComparison(compare);
  // Reduce output
  model.messageHandler()->setLogLevel(1);
  // Get branching messages
  model.messageHandler()->setLogLevel(3);

  // Do initial solve to continuous
  model.initialSolve();

  // Save model
  CbcModel model2 = model;
  int numberColumns = model.getNumCols();
  int numberIntegers = 0;
  int *integerVariable = new int[numberColumns];
  int i;
  for (i = 0; i < numberColumns; i++) {
    if (model.isInteger(i)) {
      integerVariable[numberIntegers++] = i;
    }
  }

  if (numberColumns != 180 || numberIntegers != 100) {
    printf("Incorrect model for example\n");
    exit(1);
  }

  double time1 = CoinCpuTime();

  model.branchAndBound();

  std::cout << "rgn.mps"
            << " took " << CoinCpuTime() - time1 << " seconds, "
            << model.getNodeCount() << " nodes with objective "
            << model.getObjValue()
            << (!model.status() ? " Finished" : " Not finished")
            << std::endl;

  const double *solution = model.solver()->getColSolution();

  std::cout << std::setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14);

  std::cout << "--------------------------------------" << std::endl;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    double value = solution[iColumn];
    if (fabs(value) > 1.0e-7)
      std::cout << std::setw(6) << iColumn << " " << value << std::endl;
  }
  std::cout << "--------------------------------------" << std::endl;

  std::cout << std::resetiosflags(std::ios::fixed | std::ios::showpoint | std::ios::scientific);

  // Restore model
  model = model2;

  // Convert slacks to variables
  CoinBigIndex start[5] = { 0, 1, 2, 3, 4 };
  int row[4] = { 0, 1, 2, 3 };
  double element[4] = { 1.0, 1.0, 1.0, 1.0 };
  double up[4] = { 1.0, 1.0, 1.0, 1.0 };
  model.solver()->addCols(4, start, row, element, NULL, up, NULL);
  // Now use SOS1
  int numberSets = 4;
  int which[104];
  double weights[104];
  int starts[5];
  // load
  int n = 0;
  starts[0] = 0;
  for (int iSet = 0; iSet < 4; iSet++) {
    for (int i = 0; i < 25; i++) {
      weights[n] = i + 1.0;
      which[n] = iSet * 25 + i;
      n++;
    }
    // slack - make sure first branch is on slack
    weights[n] = 1000.0;
    which[n] = 180 + iSet;
    n++;
    starts[iSet + 1] = n;
  }
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    // Stop being integer
    model.solver()->setContinuous(iColumn);
  }
  // save model in this state
  CbcModel modelSOS = model;
  CbcObject **objects = new CbcObject *[numberSets];
  for (i = 0; i < numberSets; i++) {
    objects[i] = new CbcSOS(&model, starts[i + 1] - starts[i], which + starts[i],
      weights, i);
  }
  model.addObjects(numberSets, objects);
  for (i = 0; i < numberSets; i++)
    delete objects[i];
  delete[] objects;

  time1 = CoinCpuTime();

  model.branchAndBound();

  std::cout << "rgn.mps"
            << " took " << CoinCpuTime() - time1 << " seconds, "
            << model.getNodeCount() << " nodes with objective "
            << model.getObjValue()
            << (!model.status() ? " Finished" : " Not finished")
            << std::endl;

  solution = model.solver()->getColSolution();

  std::cout << std::setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14);

  std::cout << "--------------------------------------" << std::endl;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    double value = solution[iColumn];
    if (fabs(value) > 1.0e-7)
      std::cout << std::setw(6) << iColumn << " " << value << std::endl;
  }
  std::cout << "--------------------------------------" << std::endl;

  std::cout << std::resetiosflags(std::ios::fixed | std::ios::showpoint | std::ios::scientific);

  // Restore SOS model
  model = modelSOS;

  // Now use SOS2
  objects = new CbcObject *[numberSets];
  for (i = 0; i < numberSets; i++) {
    objects[i] = new CbcSOS(&model, starts[i + 1] - starts[i], which + starts[i],
      weights, i, 2);
  }
  model.addObjects(numberSets, objects);
  for (i = 0; i < numberSets; i++)
    delete objects[i];
  delete[] objects;

  time1 = CoinCpuTime();

  model.branchAndBound();

  std::cout << "rgn.mps"
            << " took " << CoinCpuTime() - time1 << " seconds, "
            << model.getNodeCount() << " nodes with objective "
            << model.getObjValue()
            << (!model.status() ? " Finished" : " Not finished")
            << std::endl;

  solution = model.solver()->getColSolution();

  std::cout << std::setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14);

  std::cout << "--------------------------------------" << std::endl;

  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    double value = solution[iColumn];
    if (fabs(value) > 1.0e-7)
      std::cout << std::setw(6) << iColumn << " " << value
                << std::endl;
  }
  std::cout << "--------------------------------------" << std::endl;

  std::cout << std::resetiosflags(std::ios::fixed | std::ios::showpoint | std::ios::scientific);

  delete[] integerVariable;
  return 0;
}
