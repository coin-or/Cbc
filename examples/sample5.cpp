// $Id$
// Copyright (C) 2002, International Business Machines
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
#include "CbcCompareUser.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristicLocal.hpp"
#include "OsiClpSolverInterface.hpp"

// Cuts

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding.hpp"

// Heuristics

#include "CbcHeuristic.hpp"

// Methods of building

#include "CoinBuild.hpp"
#include "CoinModel.hpp"

#include "CoinTime.hpp"

/************************************************************************

This main program creates an integer model and then solves it

It then sets up some Cgl cut generators and calls branch and cut.

Branching is simple binary branching on integer variables.

Node selection is depth first until first solution is found and then
based on objective and number of unsatisfied integer variables.
In this example the functionality is the same as default but it is
a user comparison function.

Variable branching selection is on maximum minimum-of-up-down change
after strong branching on 5 variables closest to 0.5.

A simple rounding heuristic is used.


************************************************************************/

int main(int argc, const char *argv[])
{

  /* Define your favorite OsiSolver.

     CbcModel clones the solver so use solver1 up to the time you pass it
     to CbcModel then use a pointer to cloned solver (model.solver())
  */

  OsiClpSolverInterface solver1;
  /* From now on we can build model in a solver independent way.
     You can add rows one at a time but for large problems this is slow so
     this example uses CoinBuild or CoinModel
  */
  OsiSolverInterface *solver = &solver1;
  // Data (is exmip1.mps in Mps/Sample
  // Objective
  double objValue[] = { 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0 };
  // Lower bounds for columns
  double columnLower[] = { 2.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0 };
  // Upper bounds for columns
  double columnUpper[] = { COIN_DBL_MAX, 4.1, 1.0, 1.0, 4.0,
    COIN_DBL_MAX, COIN_DBL_MAX, 4.3 };
  // Lower bounds for row activities
  double rowLower[] = { 2.5, -COIN_DBL_MAX, -COIN_DBL_MAX, 1.8, 3.0 };
  // Upper bounds for row activities
  double rowUpper[] = { COIN_DBL_MAX, 2.1, 4.0, 5.0, 15.0 };
  // Matrix stored packed
  int column[] = { 0, 1, 3, 4, 7,
    1, 2,
    2, 5,
    3, 6,
    4, 7 };
  double element[] = { 3.0, 1.0, -2.0, -1.0, -1.0,
    2.0, 1.1,
    1.0, 1.0,
    2.8, -1.2,
    1.0, 1.9 };
  int starts[] = { 0, 5, 7, 9, 11, 13 };
  // Integer variables (note upper bound already 1.0)
  int whichInt[] = { 2, 3 };
  int numberRows = (int)(sizeof(rowLower) / sizeof(double));
  int numberColumns = (int)(sizeof(columnLower) / sizeof(double));
#define BUILD 2
#if BUILD == 1
  // Using CoinBuild
  // First do columns (objective and bounds)
  int i;
  // We are not adding elements
  for (i = 0; i < numberColumns; i++) {
    solver->addCol(0, NULL, NULL, columnLower[i], columnUpper[i],
      objValue[i]);
  }
  // mark as integer
  for (i = 0; i < (int)(sizeof(whichInt) / sizeof(int)); i++)
    solver->setInteger(whichInt[i]);
  // Now build rows
  CoinBuild build;
  for (i = 0; i < numberRows; i++) {
    int startRow = starts[i];
    int numberInRow = starts[i + 1] - starts[i];
    build.addRow(numberInRow, column + startRow, element + startRow,
      rowLower[i], rowUpper[i]);
  }
  // add rows into solver
  solver->addRows(build);
#else
  /* using CoinModel - more flexible but still beta.
     Can do exactly same way but can mix and match much more.
     Also all operations are on building object
  */
  CoinModel build;
  // First do columns (objective and bounds)
  int i;
  for (i = 0; i < numberColumns; i++) {
    build.setColumnBounds(i, columnLower[i], columnUpper[i]);
    build.setObjective(i, objValue[i]);
  }
  // mark as integer
  for (i = 0; i < (int)(sizeof(whichInt) / sizeof(int)); i++)
    build.setInteger(whichInt[i]);
  // Now build rows
  for (i = 0; i < numberRows; i++) {
    int startRow = starts[i];
    int numberInRow = starts[i + 1] - starts[i];
    build.addRow(numberInRow, column + startRow, element + startRow,
      rowLower[i], rowUpper[i]);
  }
  // add rows into solver
  solver->loadFromCoinModel(build);
#endif

  // Pass to solver
  CbcModel model(*solver);
  model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);

  // Set up some cut generators and defaults
  // Probing first as gets tight bounds on continuous

  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(3);
  generator1.setMaxProbe(100);
  generator1.setMaxLook(50);
  generator1.setRowCuts(3);
  //  generator1.snapshot(*model.solver());
  //generator1.createCliques(*model.solver(),2,1000,true);
  //generator1.setMode(0);

  CglGomory generator2;
  // try larger limit
  generator2.setLimit(300);

  CglKnapsackCover generator3;

  CglOddHole generator4;
  generator4.setMinimumViolation(0.005);
  generator4.setMinimumViolationPer(0.00002);
  // try larger limit
  generator4.setMaximumEntries(200);

  CglClique generator5;
  generator5.setStarCliqueReport(false);
  generator5.setRowCliqueReport(false);

  CglMixedIntegerRounding mixedGen;
  CglFlowCover flowGen;

  // Add in generators
  model.addCutGenerator(&generator1, -1, "Probing");
  model.addCutGenerator(&generator2, -1, "Gomory");
  model.addCutGenerator(&generator3, -1, "Knapsack");
  model.addCutGenerator(&generator4, -1, "OddHole");
  model.addCutGenerator(&generator5, -1, "Clique");
  model.addCutGenerator(&flowGen, -1, "FlowCover");
  model.addCutGenerator(&mixedGen, -1, "MixedIntegerRounding");

  OsiClpSolverInterface *osiclp = dynamic_cast< OsiClpSolverInterface * >(model.solver());
  // go faster stripes
  if (osiclp->getNumRows() < 300 && osiclp->getNumCols() < 500) {
    osiclp->setupForRepeatedUse(2, 0);
    printf("trying slightly less reliable but faster version (? Gomory cuts okay?)\n");
    printf("may not be safe if doing cuts in tree which need accuracy (level 2 anyway)\n");
  }

  // Allow rounding heuristic

  CbcRounding heuristic1(model);
  model.addHeuristic(&heuristic1);

  // And local search when new solution found

  CbcHeuristicLocal heuristic2(model);
  model.addHeuristic(&heuristic2);

  // Redundant definition of default branching (as Default == User)
  CbcBranchUserDecision branch;
  model.setBranchingMethod(&branch);

  // Definition of node choice
  CbcCompareUser compare;
  model.setNodeComparison(compare);

  // Do initial solve to continuous
  model.initialSolve();

  // Could tune more
  model.setMinimumDrop(CoinMin(1.0,
    fabs(model.getMinimizationObjValue()) * 1.0e-3 + 1.0e-4));

  if (model.getNumCols() < 500)
    model.setMaximumCutPassesAtRoot(-100); // always do 100 if possible
  else if (model.getNumCols() < 5000)
    model.setMaximumCutPassesAtRoot(100); // use minimum drop
  else
    model.setMaximumCutPassesAtRoot(20);
  //model.setMaximumCutPasses(5);

  // Switch off strong branching if wanted
  // model.setNumberStrong(0);
  // Do more strong branching if small
  if (model.getNumCols() < 5000)
    model.setNumberStrong(10);

  model.solver()->setIntParam(OsiMaxNumIterationHotStart, 100);

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
  double time1 = CoinCpuTime();

  // Do complete search

  model.branchAndBound();

  std::cout << " Branch and cut took " << CoinCpuTime() - time1 << " seconds, "
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
  // Print solution if any - we can't get names from Osi!

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
