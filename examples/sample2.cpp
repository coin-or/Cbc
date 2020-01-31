// $Id$
// Copyright (C) 2002, International Business Machines
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
#include "CbcStrategy.hpp"
#include "CbcHeuristicLocal.hpp"
#include "OsiClpSolverInterface.hpp"

// Cuts

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglRedSplit.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
// Preprocessing
#include "CglPreProcess.hpp"

// Heuristics

#include "CbcHeuristic.hpp"

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

Preprocessing may be selected in two ways - the second is preferred now

************************************************************************/
#define PREPROCESS 2

int main(int argc, const char *argv[])
{

  // Define your favorite OsiSolver

  OsiClpSolverInterface solver1;

  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string dirsep(1, CoinFindDirSeparator());
  std::string mpsFileName;
#if defined(SAMPLEDIR)
  mpsFileName = SAMPLEDIR;
  mpsFileName += dirsep + "p0033.mps";
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
  if (nGoodParam != argc - 2 && argc >= 2) {
    printf("Usage <file> [preprocess] [time <minutes>] or <file> <minutes>\n");
    exit(1);
  }
  solver1.initialSolve();
  // Reduce printout
  solver1.setHintParam(OsiDoReducePrint, true, OsiHintTry);
  // See if we want preprocessing
  OsiSolverInterface *solver2 = &solver1;
#if PREPROCESS == 1
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
#endif
  CbcModel model(*solver2);
  model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  // Set up some cut generators and defaults
  // Probing first as gets tight bounds on continuous

  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(1);
  generator1.setMaxPassRoot(5);
  // Number of unsatisfied variables to look at
  generator1.setMaxProbe(10);
  generator1.setMaxProbeRoot(1000);
  // How far to follow the consequences
  generator1.setMaxLook(50);
  generator1.setMaxLookRoot(500);
  // Only look at rows with fewer than this number of elements
  generator1.setMaxElements(200);
  generator1.setRowCuts(3);

  CglGomory generator2;
  // try larger limit
  generator2.setLimit(300);

  CglKnapsackCover generator3;

  CglRedSplit generator4;
  // try larger limit
  generator4.setLimit(200);

  CglClique generator5;
  generator5.setStarCliqueReport(false);
  generator5.setRowCliqueReport(false);

  CglMixedIntegerRounding2 mixedGen;
  CglFlowCover flowGen;

  // Add in generators
  // Experiment with -1 and -99 etc
  model.addCutGenerator(&generator1, -1, "Probing");
  model.addCutGenerator(&generator2, -1, "Gomory");
  model.addCutGenerator(&generator3, -1, "Knapsack");
  // model.addCutGenerator(&generator4,-1,"RedSplit");
  model.addCutGenerator(&generator5, -1, "Clique");
  model.addCutGenerator(&flowGen, -1, "FlowCover");
  model.addCutGenerator(&mixedGen, -1, "MixedIntegerRounding");
  // Say we want timings
  int numberGenerators = model.numberCutGenerators();
  int iGenerator;
  for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CbcCutGenerator *generator = model.cutGenerator(iGenerator);
    generator->setTiming(true);
  }
  OsiClpSolverInterface *osiclp = dynamic_cast< OsiClpSolverInterface * >(model.solver());
  // go faster stripes
  if (osiclp) {
    // Turn this off if you get problems
    // Used to be automatically set
    osiclp->setSpecialOptions(128);
    if (osiclp->getNumRows() < 300 && osiclp->getNumCols() < 500) {
      //osiclp->setupForRepeatedUse(2,0);
      osiclp->setupForRepeatedUse(0, 0);
    }
  }
  // Uncommenting this should switch off all CBC messages
  // model.messagesPointer()->setDetailMessages(10,10000,NULL);
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
  double objValue = model.solver()->getObjSense() * model.solver()->getObjValue();
  double minimumDropA = CoinMin(1.0, fabs(objValue) * 1.0e-3 + 1.0e-4);
  double minimumDrop = fabs(objValue) * 1.0e-4 + 1.0e-4;
  printf("min drop %g (A %g)\n", minimumDrop, minimumDropA);
  model.setMinimumDrop(minimumDrop);

  if (model.getNumCols() < 500)
    model.setMaximumCutPassesAtRoot(-100); // always do 100 if possible
  else if (model.getNumCols() < 5000)
    model.setMaximumCutPassesAtRoot(100); // use minimum drop
  else
    model.setMaximumCutPassesAtRoot(20);
  model.setMaximumCutPasses(10);
  //model.setMaximumCutPasses(2);

  // Switch off strong branching if wanted
  // model.setNumberStrong(0);
  // Do more strong branching if small
  if (model.getNumCols() < 5000)
    model.setNumberStrong(10);
  model.setNumberStrong(20);
  //model.setNumberStrong(5);
  model.setNumberBeforeTrust(5);

  model.solver()->setIntParam(OsiMaxNumIterationHotStart, 100);

  // If time is given then stop after that number of minutes
  if (minutes >= 0.0) {
    std::cout << "Stopping after " << minutes << " minutes" << std::endl;
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
  //model.messageHandler()->setLogLevel(2);
  //model.solver()->messageHandler()->setLogLevel(2);
  //model.setPrintFrequency(50);
  //#define DEBUG_CUTS
#ifdef DEBUG_CUTS
  // Set up debugger by name (only if no preprocesing)
  if (!preProcess) {
    std::string problemName;
    model.solver()->getStrParam(OsiProbName, problemName);
    model.solver()->activateRowCutDebugger(problemName.c_str());
  }
#endif
#if PREPROCESS == 2
  // Default strategy will leave cut generators as they exist already
  // so cutsOnlyAtRoot (1) ignored
  // numberStrong (2) is 5 (default)
  // numberBeforeTrust (3) is 5 (default is 0)
  // printLevel (4) defaults (0)
  CbcStrategyDefault strategy(true, 5, 5);
  // Set up pre-processing to find sos if wanted
  if (preProcess)
    strategy.setupPreProcessing(2);
  model.setStrategy(strategy);
#endif
  // Do complete search

  model.branchAndBound();

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
  // Print solution if finished - we can't get names from Osi! - so get from OsiClp

  if (model.getMinimizationObjValue() < 1.0e50) {
#if PREPROCESS == 1
    // post process
    OsiSolverInterface *solver;
    if (preProcess) {
      process.postProcess(*model.solver());
      // Solution now back in solver1
      solver = &solver1;
    } else {
      solver = model.solver();
    }
#else
    OsiSolverInterface *solver = model.solver();
#endif
    int numberColumns = solver->getNumCols();

    const double *solution = solver->getColSolution();

    // Get names from solver1 (as OsiSolverInterface may lose)
    std::vector< std::string > columnNames = *solver1.getModelPtr()->columnNames();

    int iColumn;
    std::cout << std::setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14);

    std::cout << "--------------------------------------" << std::endl;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      if (fabs(value) > 1.0e-7 && solver->isInteger(iColumn))
        std::cout << std::setw(6) << iColumn << " "
                  << columnNames[iColumn] << " "
                  << value << std::endl;
    }
    std::cout << "--------------------------------------" << std::endl;

    std::cout << std::resetiosflags(std::ios::fixed | std::ios::showpoint | std::ios::scientific);
  }
  return 0;
}
