// Copyright (C) 2006, International Business Machines
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
#include "ClpDynamicInterface.hpp"
// Time
#include "CoinTime.hpp"


/************************************************************************

This main program reads in mkc from an mps file.

It then decomposes problem - subproblems are multi-color knapsack

This is to try and let Cbc use a dynamic model

************************************************************************/

int main (int argc, const char *argv[])
{

  // Define a Solver which inherits from OsiClpsolverInterface -> OsiSolverInterface
  
  ClpDynamicInterface solver1;

  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName = "../../Data/miplib3/mkc";
  if (argc>=2) mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(),"");
  assert(numMpsReadErrors==0);

  // This clones solver 
  CbcModel model(solver1);
  OsiSolverInterface * solver2 = model.solver();
  ClpDynamicInterface * osiclp = dynamic_cast< ClpDynamicInterface*> (solver2);
  assert (osiclp);
  // Setup data
  int numberRows = osiclp->getNumRows();
  int * block = new int[numberRows];
  int i;
  for (i=0;i<numberRows;i++)
    block[i]=-2;
  if (numberRows==67) {
    // ltw
    for (i=0;i<58;i++)
      block[i]=-1;
    block[66]=-1;
  } else if (numberRows==3411) {
    // mkc
    block[0]=-1;
    block[1]=-1;
    for (i=26;i<465;i++)
      block[i]=-1;
  } else {
    printf("unknown model\n");
    exit(1);;
  }

  // Now that we know master rows let solver decompose it
  osiclp->initialize(block,&model);
  delete [] block;

  // Allow fake heuristic to move serendipity solutions across
  CbcHeuristicDynamic heuristic2(model);
  model.addHeuristic(&heuristic2);

  // Do initial solve to continuous
  model.initialSolve();

  // Switch off strong branching and trust
  model.setNumberStrong(0);
  model.setNumberBeforeTrust(0);
  model.messageHandler()->setLogLevel(2);
  model.solver()->messageHandler()->setLogLevel(1);
  
  double time1 = CoinCpuTime();

  // Do complete search
  
  model.branchAndBound();

  std::cout<<argv[1]<<" took "<<CoinCpuTime()-time1<<" seconds, "
	   <<model.getNodeCount()<<" nodes with objective "
	   <<model.getObjValue()
	   <<(!model.status() ? " Finished" : " Not finished")
	   <<std::endl;
  // Print solution if finished - we can't get names from Osi!

  if (!model.status()&&model.getMinimizationObjValue()<1.0e50) {
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
