// $Id$
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
// For Branch and bound
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcBranchActual.hpp"
#include "CbcCompareUser.hpp"

#include "CoinTime.hpp"

/************************************************************************

This main program reads in an integer model from an mps file.
It expects it to be unit coefficients and unit rhs.

Branching is follow-on branching plus simple binary branching on integer variables.

*/
int main(int argc, const char *argv[])
{

  // Define a Solver for long thin problems

  OsiClpSolverInterface solver1;

  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName;
#if defined(SAMPLEDIR)
  mpsFileName = SAMPLEDIR "/p0033.mps";
#else
  if (argc < 2) {
    fprintf(stderr, "Do not know where to find sample MPS files.\n");
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
  double time1 = CoinCpuTime();

  solver1.initialSolve();
  // Reduce printout
  solver1.setHintParam(OsiDoReducePrint, true, OsiHintTry);

  OsiSolverInterface *solver2 = &solver1;
  CbcModel model(*solver2);
  // Point to solver
  OsiSolverInterface *solver3 = model.solver();
  assert(dynamic_cast< OsiClpSolverInterface * >(solver3));

  // Definition of node choice
  CbcCompareUser compare;
  model.setNodeComparison(compare);

  int iColumn;
  int numberColumns = solver3->getNumCols();
  /* We are going to add a single follow on object but we
     want to give low priority to existing integers.
     As the default priority is 1000 we don't actually need to give
     integer priorities but it is here to show how.
  */
  // Normal integer priorities
  int *priority = new int[numberColumns];
  int numberIntegers = 0;
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (solver3->isInteger(iColumn)) {
      priority[numberIntegers++] = 100; // low priority
    }
  }
  /* Second parameter is true if we are adding objects,
     false if integers.  So this does integers */
  model.passInPriorities(priority, false);
  delete[] priority;
  /* Add in objects before we can give priority.
     In this case just one - but this shows general method
  */
  CbcObject **objects = new CbcObject *[1];
  objects[0] = new CbcFollowOn(&model);
  model.addObjects(1, objects);
  delete objects[0];
  delete[] objects;
  // High priority
  int followPriority = 1;
  model.passInPriorities(&followPriority, true);

  // Do initial solve to continuous
  model.initialSolve();

  // Do more strong branching if small
  // Switch off strong branching if wanted
  model.setNumberStrong(0);

  model.solver()->setIntParam(OsiMaxNumIterationHotStart, 100);

  // Switch off most output
  if (model.getNumCols() < 3000) {
    model.messageHandler()->setLogLevel(1);
    //model.solver()->messageHandler()->setLogLevel(0);
  } else {
    model.messageHandler()->setLogLevel(2);
    model.solver()->messageHandler()->setLogLevel(1);
  }
  //model.setPrintFrequency(50);

  // Do complete search
  try {
    model.branchAndBound();
  } catch (CoinError e) {
    e.print();
    if (e.lineNumber() >= 0)
      std::cout << "This was from a CoinAssert" << std::endl;
    exit(0);
  }
  //void printHowMany();
  //printHowMany();
  std::cout << mpsFileName << " took " << CoinCpuTime() - time1 << " seconds, "
            << model.getNodeCount() << " nodes with objective "
            << model.getObjValue()
            << (!model.status() ? " Finished" : " Not finished")
            << std::endl;

  // Print solution if finished - we can't get names from Osi!

  if (model.getMinimizationObjValue() < 1.0e50) {
    int numberColumns = model.solver()->getNumCols();

    const double *solution = model.solver()->getColSolution();

    int iColumn;
    std::cout << std::setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14);

    std::cout << "--------------------------------------" << std::endl;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      if (fabs(value) > 1.0e-7 && model.solver()->isInteger(iColumn))
        std::cout << std::setw(6) << iColumn << " " << value << std::endl;
    }
    std::cout << "--------------------------------------" << std::endl;

    std::cout << std::resetiosflags(std::ios::fixed | std::ios::showpoint | std::ios::scientific);
  }
  return 0;
}
