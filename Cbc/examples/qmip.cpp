// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"

// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcBranchUser.hpp"
#include "CbcCompareUser.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristicLocal.hpp"
#include "ClpQuadInterface.hpp"

// Cuts - some may work but be very very careful

#include "CglProbing.hpp"

// Heuristics would need adapting

#include "CbcHeuristic.hpp"

// Time
#include "CoinTime.hpp"

/************************************************************************

This main program reads in a quadratic integer model from an mps file.

It then sets up some Cgl cut generators and calls branch and cut.

Branching is simple binary branching on integer variables.

Node selection is depth first until first solution is found and then
based on objective and number of unsatisfied integer variables.
In this example the functionality is the same as default but it is
a user comparison function.

Variable branching selection is on maximum minimum-of-up-down change
after strong branching on 5 variables closest to 0.5.

A simple rounding heuristic could be used.but is switched off as needs work for quadratic

This is NOT meant to be a serious MIQP code;  it is to show how you can use
ClpQuadInterface.hpp to solve a QP at each node.  You could pick up data in that interface 
and use any QP solver (e.g. one that works)

************************************************************************/

int main(int argc, const char *argv[])
{

  // Define a Solver which inherits from OsiClpsolverInterface -> OsiSolverInterface

  ClpQuadInterface solver1;

  // Read in model using argv[1]
  if (argc <= 1) {
    printf("using %s <modelfile>\n", argv[0]);
    return 1;
  }
  // must use clp to get a quadratic model
  ClpSimplex *clp = solver1.getModelPtr();
  int numMpsReadErrors = clp->readMps(argv[1]);
  // and assert that it is a clean model
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }

  // This clones solver
  CbcModel model(solver1);
  // But now model doesn't know about integers!
  const char *integerInfo = clp->integerInformation();
  int numberColumns = clp->numberColumns();
  // and point to solver
  OsiSolverInterface *solver2 = model.solver();
  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (integerInfo[iColumn])
      solver2->setInteger(iColumn);
  }
  // Okay - now we have a good MIQP solver
  // Within branch and cut it is better (at present) to switch off all objectives

  ClpQuadInterface *osiclp = dynamic_cast< ClpQuadInterface * >(solver2);
  assert(osiclp);
  // Set fake objective so Cbc not confused
  osiclp->initialize();
  solver2->setHintParam(OsiDoReducePrint, true, OsiHintTry);

  // Set up some cut generators and defaults
  // Probing first as gets tight bounds on continuous

  CglProbing generator1;
  // can not use objective
  generator1.setUsingObjective(false);
  generator1.setMaxPass(3);
  generator1.setMaxProbe(100);
  generator1.setMaxLook(50);
  generator1.setRowCuts(3);

  // Add in generators
  // Only some generators work (and even then try without first)
  model.addCutGenerator(&generator1, 1, "Probing");
  // Allow rounding heuristic

  CbcRounding heuristic1(model);
  // do not add yet as don't know how to deal with quadratic objective
  //model.addHeuristic(&heuristic1);

  // Redundant definition of default branching (as Default == User)
  CbcBranchUserDecision branch;
  model.setBranchingMethod(&branch);

  // Definition of node choice
  CbcCompareUser compare;
  // breadth first
  //compare.setWeight(0.0);
  model.setNodeComparison(compare);

  // Do initial solve to continuous
  model.initialSolve();

  // Could tune more
  model.setMinimumDrop(CoinMin(1.0,
    fabs(model.getMinimizationObjValue()) * 1.0e-3 + 1.0e-4));

  model.setMaximumCutPassesAtRoot(0);
  model.setMaximumCutPasses(0);

  // Switch off strong branching if wanted
  //model.setNumberStrong(5);

  model.solver()->setIntParam(OsiMaxNumIterationHotStart, 10000);

  // If time is given then stop after that number of minutes
  if (argc > 2) {
    int minutes = atoi(argv[2]);
    std::cout << "Stopping after " << minutes << " minutes" << std::endl;
    assert(minutes >= 0);
    model.setDblParam(CbcModel::CbcMaximumSeconds, 60.0 * minutes);
  }
  // Switch off most output
  if (model.getNumCols() < 3000) {
    model.messageHandler()->setLogLevel(1);
    //model.solver()->messageHandler()->setLogLevel(0);
  } else {
    model.messageHandler()->setLogLevel(2);
    model.solver()->messageHandler()->setLogLevel(1);
  }
  model.setPrintFrequency(50);

  double time1 = CoinCpuTime();

  // Do complete search

  model.branchAndBound();

  std::cout << argv[1] << " took " << CoinCpuTime() - time1 << " seconds, "
            << model.getNodeCount() << " nodes with objective "
            << model.getObjValue()
            << (!model.status() ? " Finished" : " Not finished")
            << std::endl;

  // Print more statistics
  std::cout << "Cuts at root node changed objective from " << model.getContinuousObjective()
            << " to " << model.rootObjectiveAfterCuts() << std::endl;

  int numberGenerators = model.numberCutGenerators();
  for (int iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CbcCutGenerator *generator = model.cutGenerator(iGenerator);
    std::cout << generator->cutGeneratorName() << " was tried "
              << generator->numberTimesEntered() << " times and created "
              << generator->numberCutsInTotal() << " cuts of which "
              << generator->numberCutsActive() << " were active after adding rounds of cuts"
              << std::endl;
  }
  // Print solution if finished - we can't get names from Osi!

  if (!model.status() && model.getMinimizationObjValue() < 1.0e50) {
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
