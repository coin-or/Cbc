// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iomanip>


// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcBranchUser.hpp"
#include "CbcCompareUser.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristicUser.hpp"
#ifdef COIN_USE_CLP
#include "OsiClpSolverInterface.hpp"
#endif
#ifdef COIN_USE_OSL
#include "OsiOslSolverInterface.hpp"
#endif

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

#include  "CoinTime.hpp"

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

int main (int argc, const char *argv[])
{

  // Define your favorite OsiSolver
  
#ifdef COIN_USE_CLP
  OsiClpSolverInterface solver1;
  //solver1.messageHandler()->setLogLevel(0);
  CbcModel model(solver1);
#elif COIN_USE_OSL
  OsiOslSolverInterface solver1;
  //solver1.messageHandler()->setLogLevel(0);
  CbcModel model(solver1);
#endif
  model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);

  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName = "../../Mps/Sample/p0033.mps";
  if (argc>=2) mpsFileName = argv[1];
  int numMpsReadErrors = model.solver()->readMps(mpsFileName.c_str(),"");
  assert(numMpsReadErrors==0);

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
  model.addCutGenerator(&generator1,-1,"Probing");
  model.addCutGenerator(&generator2,-1,"Gomory");
  model.addCutGenerator(&generator3,-1,"Knapsack");
  model.addCutGenerator(&generator4,-1,"OddHole");
  model.addCutGenerator(&generator5,-1,"Clique");
  model.addCutGenerator(&flowGen,-1,"FlowCover");
  model.addCutGenerator(&mixedGen,-1,"MixedIntegerRounding");

#ifdef COIN_USE_CLP
  OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (model.solver());
  // go faster stripes
  if (osiclp->getNumRows()<300&&osiclp->getNumCols()<500) {
    osiclp->setupForRepeatedUse(2,0);
    printf("trying slightly less reliable but faster version (? Gomory cuts okay?)\n");
    printf("may not be safe if doing cuts in tree which need accuracy (level 2 anyway)\n");
  }
#endif

  // Allow rounding heuristic

  CbcRounding heuristic1(model);
  model.addHeuristic(&heuristic1);

  // And local search when new solution found

  CbcLocalSearch heuristic2(model);
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
  model.setMinimumDrop(min(1.0,
			     fabs(model.getMinimizationObjValue())*1.0e-3+1.0e-4));

  if (model.getNumCols()<500)
    model.setMaximumCutPassesAtRoot(-100); // always do 100 if possible
  else if (model.getNumCols()<5000)
    model.setMaximumCutPassesAtRoot(100); // use minimum drop
  else
    model.setMaximumCutPassesAtRoot(20);
  //model.setMaximumCutPasses(5);

  // Switch off strong branching if wanted
  // model.setNumberStrong(0);
  // Do more strong branching if small
  if (model.getNumCols()<5000)
    model.setNumberStrong(10);

  model.solver()->setIntParam(OsiMaxNumIterationHotStart,100);

  // If time is given then stop after that number of minutes
  if (argc>2) {
    int minutes = atoi(argv[2]);
    std::cout<<"Stopping after "<<minutes<<" minutes"<<std::endl;
    assert (minutes>=0);
    model.setDblParam(CbcModel::CbcMaximumSeconds,60.0*minutes);
  }
  // Switch off most output
  if (model.getNumCols()<3000) {
    model.messageHandler()->setLogLevel(1);
    //model.solver()->messageHandler()->setLogLevel(0);
  } else {
    model.messageHandler()->setLogLevel(2);
    model.solver()->messageHandler()->setLogLevel(1);
  }
  //model.messageHandler()->setLogLevel(2);
  //model.solver()->messageHandler()->setLogLevel(2);
  //model.setPrintFrequency(50);
  //#define DEBUG_CUTS 2
#if DEBUG_CUTS==1
  // Set up debugger by value
  {
    int numberColumns=model.getNumCols();
    double * solution = new double[numberColumns];
    memset(solution, 0,numberColumns*sizeof(double));
    // Set nonzero integers
    solution[1]=1.0;
    solution[18]=1.0;
    solution[33]=1.0;
    solution[59]=1.0;
    model.solver()->activateRowCutDebugger(solution);
    delete [] solution;
  }
#elif DEBUG_CUTS==2
  // Set up debugger by name
  {
    std::string problemName ;
    model.solver()->getStrParam(OsiProbName,problemName) ;
    model.solver()->activateRowCutDebugger(problemName.c_str()) ;
  }
#endif
  double time1 = CoinCpuTime();

  if (0) {
    // integer presolve
    CbcModel * model2 = model.integerPresolve();
    if (model2) {
      // Do complete search
  
      model2->branchAndBound();
      // get back solution
      model.originalModel(model2,false);
    } else {
      // infeasible
      exit(1);
    }
  } else {
    // Do complete search
  
    model.branchAndBound();
  }

  std::cout<<mpsFileName<<" took "<<CoinCpuTime()-time1<<" seconds, "
	   <<model.getNodeCount()<<" nodes with objective "
	   <<model.getObjValue()
	   <<(!model.status() ? " Finished" : " Not finished")
	   <<std::endl;

  printf("Exact obj %18.18g\n",model.getObjValue());
  // Print more statistics
  std::cout<<"Cuts at root node changed objective from "<<model.getContinuousObjective()
	   <<" to "<<model.rootObjectiveAfterCuts()<<std::endl;

  int numberGenerators = model.numberCutGenerators();
  for (int iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
    CbcCutGenerator * generator = model.cutGenerator(iGenerator);
    std::cout<<generator->cutGeneratorName()<<" was tried "
	     <<generator->numberTimesEntered()<<" times and created "
	     <<generator->numberCutsInTotal()<<" cuts of which "
	     <<generator->numberCutsActive()<<" were active after adding rounds of cuts"
	     <<std::endl;
  }
  // Print solution if finished - we can't get names from Osi!

  if (model.getMinimizationObjValue()<1.0e50) {
    int numberColumns = model.solver()->getNumCols();
    
    const double * solution = model.solver()->getColSolution();
    
    int iColumn;
    std::cout<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);
    
    std::cout<<"--------------------------------------"<<std::endl;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double value=solution[iColumn];
      if (fabs(value)>1.0e-7&&model.solver()->isInteger(iColumn)) 
	std::cout<<std::setw(6)<<iColumn<<" "<<value<<std::endl;
    }
    std::cout<<"--------------------------------------"<<std::endl;
  
    std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);
  }
  return 0;
}    
