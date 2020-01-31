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
#include "CbcBranchCut.hpp"
#include "CbcBranchToFixLots.hpp"
#include "CbcCompareUser.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcSolver3.hpp"

// Cuts

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding.hpp"
#include "CglMixedIntegerRounding2.hpp"
// Preprocessing
#include "CglPreProcess.hpp"

#include "CoinTime.hpp"

//#############################################################################

/************************************************************************

This main program reads in an integer model from an mps file.

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

// ****** define comparison to choose best next node

int main(int argc, const char *argv[])
{

  CbcSolver3 solver1;

  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName;
  if (argc >= 2)
    mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(), "");
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }
  double time1 = CoinCpuTime();

  /* Options are:
     preprocess to do preprocessing
     time in minutes
     if 2 parameters and numeric taken as time
  */
  bool preProcess = false;
  double minutes = -1.0;
  int nGoodParam = 0;
  for (int iParam = 2; iParam < argc; iParam++) {
    if (!strcmp(argv[iParam], "preprocess")) {
      preProcess = true;
      nGoodParam++;
    } else if (!strcmp(argv[iParam], "time")) {
      if (iParam + 1 < argc && isdigit(argv[iParam + 1][0])) {
        minutes = atof(argv[iParam + 1]);
        if (minutes >= 0.0) {
          nGoodParam += 2;
          iParam++; // skip time
        }
      }
    }
  }
  if (nGoodParam == 0 && argc == 3 && isdigit(argv[2][0])) {
    // If time is given then stop after that number of minutes
    minutes = atof(argv[2]);
    if (minutes >= 0.0)
      nGoodParam = 1;
  }
  if (nGoodParam != argc - 2) {
    printf("Usage <file> [preprocess] [time <minutes>] or <file> <minutes>\n");
    exit(1);
  }
  solver1.initialSolve();
  // Reduce printout
  solver1.setHintParam(OsiDoReducePrint, true, OsiHintTry);
  // Say we want scaling
  //solver1.setHintParam(OsiDoScale,true,OsiHintTry);
  //solver1.setCleanupScaling(1);
  // See if we want preprocessing
  OsiSolverInterface *solver2 = &solver1;
  CglPreProcess process;
  if (preProcess) {
    /* Do not try and produce equality cliques and
       do up to 5 passes */
    solver2 = process.preProcess(solver1, false, 5);
    if (!solver2) {
      printf("Pre-processing says infeasible\n");
      exit(2);
    }
    solver2->resolve();
  }
  CbcModel model(*solver2);
  // Point to solver
  OsiSolverInterface *solver3 = model.solver();
  CbcSolver3 *osiclp = dynamic_cast< CbcSolver3 * >(solver3);
  assert(osiclp);
  const double fractionFix = 0.985;
  osiclp->initialize(&model, NULL);
  osiclp->setAlgorithm(2);
  osiclp->setMemory(1000);
  osiclp->setNested(fractionFix);
  //osiclp->setNested(1.0); //off
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
  /* This is same as default constructor - 
     (1,true,1)
      I presume if maxAggregate larger then slower but maybe better
      criterion can be  1 through 3
      Reference: 
      Hugues Marchand and Laurence A. Wolsey
      Aggregation and Mixed Integer Rounding to Solve MIPs
      Operations Research, 49(3), May-June 2001.
   */
  int maxAggregate = 1;
  bool multiply = true;
  int criterion = 1;
  CglMixedIntegerRounding2 mixedGen2(maxAggregate, multiply, criterion);
  CglFlowCover flowGen;

  // Add in generators
  // Experiment with -1 and -99 etc
  model.addCutGenerator(&generator1, -99, "Probing");
  //model.addCutGenerator(&generator2,-1,"Gomory");
  //model.addCutGenerator(&generator3,-1,"Knapsack");
  //model.addCutGenerator(&generator4,-1,"OddHole");
  //model.addCutGenerator(&generator5,-1,"Clique");
  //model.addCutGenerator(&flowGen,-1,"FlowCover");
  //model.addCutGenerator(&mixedGen,-1,"MixedIntegerRounding");
  //model.addCutGenerator(&mixedGen2,-1,"MixedIntegerRounding2");
  // Say we want timings
  int numberGenerators = model.numberCutGenerators();
  int iGenerator;
  for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CbcCutGenerator *generator = model.cutGenerator(iGenerator);
    generator->setTiming(true);
  }

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
  CbcObject **objects = new CbcObject *[numberColumns + 1];
  const CoinPackedMatrix *matrix = solver3->getMatrixByCol();
  // Column copy
  const int *columnLength = matrix->getVectorLengths();
  const double *objective = model.getObjCoefficients();
  int n = 0;
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (solver3->isInteger(iColumn)) {
      double costPer = objective[iColumn] / ((double)columnLength[iColumn]);
      CbcSimpleIntegerPseudoCost *newObject = new CbcSimpleIntegerPseudoCost(&model, n, iColumn,
        costPer, costPer);
      newObject->setMethod(3);
      objects[n++] = newObject;
    }
  }
  // and special fix lots branch
  objects[n++] = new CbcBranchToFixLots(&model, -1.0e-6, fractionFix + 0.01, 1, 0, NULL);
  model.addObjects(n, objects);
  for (iColumn = 0; iColumn < n; iColumn++)
    delete objects[iColumn];
  delete[] objects;
  // High priority for odd object
  int followPriority = 1;
  model.passInPriorities(&followPriority, true);

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
  //model.setMaximumCutPasses(1);

  // Do more strong branching if small
  //if (model.getNumCols()<5000)
  //model.setNumberStrong(10);
  // Switch off strong branching if wanted
  model.setNumberStrong(0);

  model.solver()->setIntParam(OsiMaxNumIterationHotStart, 100);

  // If time is given then stop after that number of minutes
  if (minutes >= 0.0) {
    std::cout << "Stopping after " << minutes << " minutes" << std::endl;
    model.setDblParam(CbcModel::CbcMaximumSeconds, 60.0 * minutes);
  }
  // Switch off most output
  if (model.getNumCols() < 300000) {
    model.messageHandler()->setLogLevel(1);
    //model.solver()->messageHandler()->setLogLevel(0);
  } else {
    model.messageHandler()->setLogLevel(2);
    model.solver()->messageHandler()->setLogLevel(1);
  }
  //model.messageHandler()->setLogLevel(2);
  //model.solver()->messageHandler()->setLogLevel(2);
  //model.setPrintFrequency(50);
#define DEBUG_CUTS
#ifdef DEBUG_CUTS
  // Set up debugger by name (only if no preprocesing)
  if (!preProcess) {
    std::string problemName;
    //model.solver()->getStrParam(OsiProbName,problemName) ;
    //model.solver()->activateRowCutDebugger(problemName.c_str()) ;
    model.solver()->activateRowCutDebugger("cap6000a");
  }
#endif

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

  // Print more statistics
  std::cout << "Cuts at root node changed objective from " << model.getContinuousObjective()
            << " to " << model.rootObjectiveAfterCuts() << std::endl;

  for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CbcCutGenerator *generator = model.cutGenerator(iGenerator);
    std::cout << generator->cutGeneratorName() << " was tried "
              << generator->numberTimesEntered() << " times and created "
              << generator->numberCutsInTotal() << " cuts of which "
              << generator->numberCutsActive() << " were active after adding rounds of cuts";
    if (generator->timing())
      std::cout << " ( " << generator->timeInCutGenerator() << " seconds)" << std::endl;
    else
      std::cout << std::endl;
  }
  // Print solution if finished - we can't get names from Osi!

  if (model.getMinimizationObjValue() < 1.0e50) {
    // post process
    if (preProcess)
      process.postProcess(*model.solver());
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
