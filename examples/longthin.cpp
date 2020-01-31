// $Id$
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcBranchActual.hpp"
#include "CbcBranchUser.hpp"
#include "CbcCompareUser.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcSolver2.hpp"
#include "CoinModel.hpp"

// Cuts

#include "CglProbing.hpp"

#include "CoinTime.hpp"

/************************************************************************

This main program reads in an integer model from an mps file.
It expects it to be unit coefficients and unit rhs and long and thin

Branching is simple binary branching on integer variables.

*/
int main(int argc, const char *argv[])
{

  // Define a Solver for long thin problems

  CbcSolver2 solver1;

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
  CbcSolver2 *osiclp = dynamic_cast< CbcSolver2 * >(solver3);
  assert(osiclp);
  osiclp->initialize(&model, NULL);
  osiclp->setAlgorithm(2);
  osiclp->setMemory(1000);
  // Set up some cut generators and defaults
  // Probing first as gets tight bounds on continuous

  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(3);
  // Number of unsatisfied variables to look at
  generator1.setMaxProbe(10);
  // How far to follow the consequences
  generator1.setMaxLook(50);
  // Only look at rows with fewer than this number of elements
  generator1.setMaxElements(200);
  generator1.setRowCuts(3);

  // Add in generators
  // Experiment with -1 and -99 etc
  model.addCutGenerator(&generator1, -99, "Probing");
  // Allow rounding heuristic

  CbcRounding heuristic1(model);
  model.addHeuristic(&heuristic1);

  // And Greedy heuristic

  CbcHeuristicGreedyCover heuristic2(model);
  // Use original upper and perturb more
  heuristic2.setAlgorithm(11);
  model.addHeuristic(&heuristic2);

  // Redundant definition of default branching (as Default == User)
  CbcBranchUserDecision branch;
  model.setBranchingMethod(&branch);

  // Definition of node choice
  CbcCompareUser compare;
  model.setNodeComparison(compare);

  int iColumn;
  int numberColumns = solver3->getNumCols();
  // do pseudo costs
  CbcObject **objects = new CbcObject *[numberColumns];
  const CoinPackedMatrix *matrix = solver3->getMatrixByCol();
  // Column copy
  const int *columnLength = matrix->getVectorLengths();
  const double *objective = model.getObjCoefficients();
  int numberIntegers = 0;
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (solver3->isInteger(iColumn)) {
      /*  Branching up gets us much closer to an integer solution so we want
          to encourage up - so we will branch up if variable value > 0.333333.
          The expected cost of going up obviously depends on the cost of the
          variable so we just choose pseudo costs to reflect that.  We could also
          decide to try and use the pseudo costs to make it more likely to branch
          on a variable with many coefficients.  This leads to the computation below.
      */
      double cost = objective[iColumn] * (1.0 + 0.2 * ((double)columnLength[iColumn]));
      CbcSimpleIntegerPseudoCost *newObject = new CbcSimpleIntegerPseudoCost(&model, iColumn,
        2.0 * cost, cost);
      newObject->setMethod(3);
      objects[numberIntegers++] = newObject;
    }
  }
  model.addObjects(numberIntegers, objects);
  for (iColumn = 0; iColumn < numberIntegers; iColumn++)
    delete objects[iColumn];
  delete[] objects;

  // Do initial solve to continuous
  model.initialSolve();

  // Do more strong branching if small
  // Switch off strong branching if wanted
  model.setNumberStrong(5);

  // say use resolve for strong branching
  osiclp->setSpecialOptions(16);
  // We had better allow a lot
  model.solver()->setIntParam(OsiMaxNumIterationHotStart, 10000);
  // So use strategy to keep rows
  osiclp->setStrategy(1);

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
