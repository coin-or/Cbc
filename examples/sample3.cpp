// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcConfig.h"

#include "CoinPragma.hpp"

#include <assert.h>
#include <iomanip>

// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcBranchUser.hpp"
#include "CbcCompareUser.hpp"
#include "CbcTreeLocal.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristicLocal.hpp"
#include "OsiClpSolverInterface.hpp"

// Cuts

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"

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


************************************************************************/
// See if option there and return integer value - or -123456
bool option(int argc, const char *argv[], const char *check, int &value)
{
  value = -123456;
  for (int i = 1; i < argc; i++) {
    if (const char *x = strstr(argv[i], check)) {
      // see if at beginning
      if (x == argv[i] || (x == argv[i] + 1 && argv[i][0] == '-')) {
        // see if =
        x = strchr(argv[i], '=');
        if (x) {
          value = atoi(x + 1);
        } else if (i < argc - 1) {
          value = atoi(argv[i + 1]);
        }
        return true;
      }
    }
  }
  return false;
}

// ****** define comparison to choose best next node

int main(int argc, const char *argv[])
{

  // Define your favorite OsiSolver

  OsiClpSolverInterface solver1;
  //solver1.messageHandler()->setLogLevel(0);
  solver1.getModelPtr()->setDualBound(1.0e10);

  if (argc <= 1) {
    printf("using %s <modelfile>\n", argv[0]);
    return 1;
  }

  // Read in model using argv[1]
  // and assert that it is a clean model
  int numMpsReadErrors = solver1.readMps(argv[1], "");
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }
  // do here so integers correct
  CbcModel model(solver1);
  model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);

  // Set up some cut generators and defaults
  // Probing first as gets tight bounds on continuous

  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(3);
  generator1.setMaxProbe(100);
  generator1.setMaxLook(50);
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

  // Add in generators
  model.addCutGenerator(&generator1, -1, "Probing");
  // You may want to switch these off for local search on some problems
  //model.addCutGenerator(&generator2,-1,"Gomory");
  model.addCutGenerator(&generator2, -99, "Gomory");
  //model.addCutGenerator(&generator3,-1,"Knapsack");
  model.addCutGenerator(&generator4, -99, "OddHole");
  model.addCutGenerator(&generator5, -99, "Clique");

  // Allow rounding heuristic

  CbcRounding heuristic1(model);
  model.addHeuristic(&heuristic1);

  // Redundant definition of default branching (as Default == User)
  CbcBranchUserDecision branch;
  model.setBranchingMethod(&branch);

  // Definition of node choice
  CbcCompareUser compare;
  model.setNodeComparison(compare);

  // Do initial solve to continuous
  model.initialSolve();
  //model.solver()->setHintParam(OsiDoScale,false,OsiHintTry);

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

  // Do more strong branching if small
  if (model.getNumCols() < 5000)
    model.setNumberStrong(10);
  // Switch off strong branching if wanted
  //model.setNumberStrong(0);

  model.solver()->setIntParam(OsiMaxNumIterationHotStart, 100);

  // If time is given then stop after that number of minutes
  if (argc == 3) {
    int minutes = atoi(argv[2]);
    if (minutes > 0 && minutes < 1000000) {
      std::cout << "Stopping after " << minutes << " minutes" << std::endl;
      model.setDblParam(CbcModel::CbcMaximumSeconds, 60.0 * minutes);
      argc = 0; // stop other options
    }
  }
  int optionValue;
  if (option(argc, argv, "minute", optionValue)) {
    int minutes = optionValue;
    if (minutes > 0 && minutes < 1000000) {
      std::cout << "Stopping after " << minutes << " minutes" << std::endl;
      model.setDblParam(CbcModel::CbcMaximumSeconds, 60.0 * minutes);
    }
  }

  // Switch off most output
  if (model.getNumCols() < 3000) {
    model.messageHandler()->setLogLevel(1);
    //model.solver()->messageHandler()->setLogLevel(0);
  } else {
    model.messageHandler()->setLogLevel(2);
    model.solver()->messageHandler()->setLogLevel(1);
  }
  if (option(argc, argv, "print", optionValue)) {
    if (optionValue == 3) {
      model.messageHandler()->setLogLevel(3);
      model.solver()->messageHandler()->setLogLevel(3);
    } else {
      model.messageHandler()->setLogLevel(2);
      model.solver()->messageHandler()->setLogLevel(1);
    }
  }
  //model.setPrintFrequency(50);

  double time1 = CoinCpuTime();

  if (option(argc, argv, "clique", optionValue)) {
    printf("Finding cliques\n");
    // Find cliques
    int numberColumns = model.getNumCols();
    int numberIntegers = 0;
    int *integerVariable = new int[numberColumns];
    int i;
    for (i = 0; i < numberColumns; i++) {
      if (model.isInteger(i)) {
        integerVariable[numberIntegers++] = i;
        //model.solver()->setContinuous(i);
      }
    }
    model.findCliques(false, 2, 1000);
    // give cliques high priority but leave integers to improve time
    int numberCliques = model.numberObjects() - numberIntegers;
    int *priority = new int[numberIntegers + numberCliques];
    for (i = 0; i < numberIntegers; i++)
      priority[i] = 10000 + i;
    for (; i < numberIntegers + numberCliques; i++)
      priority[i] = i;
    model.passInPriorities(priority, false);
    delete[] priority;
    delete[] integerVariable;
    model.messageHandler()->setLogLevel(2);
  }
  if (option(argc, argv, "presolvelocal", optionValue)) {
    // integer presolve
    CbcModel *model2 = model.integerPresolve();
    if (model2) {
      // Do complete search
      model2->branchAndBound();
      // get back solution
      model.originalModel(model2, false);
    } else {
      // infeasible
      exit(1);
    }
  } else if (option(argc, argv, "presolve", optionValue)) {
    // integer presolve
    CbcModel *model2 = model.integerPresolve();
    if (model2) {
      int maxNodes = 2000;
      int k = 10;
      CbcTreeLocal localTree(model2, NULL, k, 0, 1, 10000, maxNodes);
      model2->passInTreeHandler(localTree);
      // Do complete search
      model2->branchAndBound();
      // get back solution
      model.originalModel(model2, false);
    } else {
      // infeasible
      exit(1);
    }
  } else if (option(argc, argv, "localprove", optionValue)) {
    int maxNodes = 2000;
    int k = 10;
    CbcTreeLocal localTree(&model, NULL, k, 0, 1, 10000, maxNodes);
    model.passInTreeHandler(localTree);
    model.branchAndBound();
  } else if (option(argc, argv, "local", optionValue)) {
    int maxNodes = 2000;
    int k = 10;
    double *solution = NULL;
    // trying to get best solution to fast0507
    // (so we can debug why localTree (prove) cut off solution!
    int numberColumns = model.getNumCols();
    if (numberColumns == 63009) {
      solution = new double[numberColumns];
      memset(solution, 0, numberColumns * sizeof(double));
      int seq[] = {
        245, 246, 395, 582, 813, 1432, 2009, 2063, 2455, 2879,
        3027, 3116, 3154, 3258, 3675, 5548, 5890, 6229, 6318, 6595,
        6852, 6859, 7279, 7565, 7604, 7855, 9160, 9259, 9261, 9296,
        9318, 10080, 10339, 10349, 10652, 10874, 11046, 11396, 11812, 11968,
        12129, 12321, 12495, 12650, 13182, 13350, 14451, 15198, 15586, 17876,
        18537, 19742, 19939, 19969, 20646, 20667, 20668, 23810, 25973, 27238,
        27510, 27847, 30708, 31164, 33189, 36506, 36881, 38025, 38070, 38216,
        38381, 39076, 39397, 39708, 40770, 40831, 41017, 41164, 41195, 41227,
        41384, 41743, 42451, 42776, 42833, 43991, 44904, 45045, 45718, 47860,
        47896, 48180, 48257, 48461, 48961, 49227, 51278, 52534, 53275, 54856,
        55672, 55674, 55760, 55822, 56641, 56703, 56892, 56964, 57117, 57556,
        57905, 59936, 61649, 62319, 62387
      };
      double intSolnV[] = {
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1.
      };
      int n = sizeof(seq) / sizeof(int);
      for (int i = 0; i < n; i++)
        solution[seq[i]] = intSolnV[i];
      k = 50;
    }
    CbcTreeLocal localTree(&model, solution, k, 0, 50, 10000, maxNodes);
    delete[] solution;
    model.passInTreeHandler(localTree);
    model.branchAndBound();
  } else if (option(argc, argv, "1local", optionValue)) {
    // test giving cbc a solution
    CbcModel model2 = model;
    model2.setMaximumSolutions(1);
    model2.branchAndBound();
    if (model2.status() && model2.getMinimizationObjValue() < 1.0e50) {
      // Definition of local search
      const double *solution = model2.solver()->getColSolution();
      CbcTreeLocal localTree(&model, solution, 2);
      model.passInTreeHandler(localTree);
      model.branchAndBound();
    } else {
      model = model2;
    }
  } else {
    // Do complete search

    model.branchAndBound();
  }

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
  if (model.numberGlobalViolations())
    std::cout << "Global cuts were active " << model.numberGlobalViolations() << " times" << std::endl;
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
    int n = 0;
    int i;
    bool comma = false;
    bool newLine = false;
    for (i = 0; i < numberColumns; i++) {
      if (solution[i] > 0.5 && model.solver()->isInteger(i)) {
        if (comma)
          printf(",");
        if (newLine)
          printf("\n");
        printf("%d ", i);
        comma = true;
        newLine = false;
        n++;
        if (n == 10) {
          n = 0;
          newLine = true;
        }
      }
    }
    printf("};\n");
    n = 0;
    comma = false;
    newLine = false;
    printf("\tdouble intSolnV[]={\n");
    for (i = 0; i < numberColumns; i++) {
      if (solution[i] > 0.5 && model.solver()->isInteger(i)) {
        if (comma)
          printf(",");
        if (newLine)
          printf("\n");
        int value = (int)(solution[i] + 0.5);
        printf("%d. ", value);
        comma = true;
        newLine = false;
        n++;
        if (n == 10) {
          n = 0;
          newLine = true;
        }
      }
    }
    printf("};\n");
  }
  return 0;
}
