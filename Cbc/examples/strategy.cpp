// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iomanip>


// For Branch and bound
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#include "CbcCutGenerator.hpp"
#include "CoinTime.hpp"

//#############################################################################


/************************************************************************

This main program reads in an integer model from an mps file.

It then sets up a strategy and solves a problem

Sensible users may wish to look at CbcStrategy.hpp to see what the default means and then tweak it
using - say a subset of miplib to tune.  

CbcStrategyDefault was there as an example - the more adventurous may wish to look at 

strategy2.cpp and try adding in other possibilities 
************************************************************************/

int main (int argc, const char *argv[])
{

  // Define your favorite OsiSolver
  
  OsiClpSolverInterface solver1;

  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName = "../../Data/Sample/p0033.mps";
  if (argc>=2) mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(),"");
  assert(numMpsReadErrors==0);
  double time1 = CoinCpuTime();

  solver1.initialSolve();
  // Reduce printout
  solver1.setHintParam(OsiDoReducePrint,true,OsiHintTry);
  // See if we want preprocessing
  OsiSolverInterface * solver2=&solver1;
  CbcModel model(solver1);
  // Stop after 20 minutes
  int minutes=20;
  std::cout<<"Stopping after "<<minutes<<" minutes"<<std::endl;
  model.setDblParam(CbcModel::CbcMaximumSeconds,60.0*minutes);
  ////////////////////////////////////////////////////////////

  CbcStrategyDefault strategy;
  // Default strategy will leave cut generators as they exist already
  // so cutsOnlyAtRoot (1) true
  // numberStrong (2) is 5 (default)
  // numberBeforeTrust (3) is 0 (default)
  // printLevel (4) defaults (0)
  // So same as CbcStrategyDefault strategy(true,5,0,0);
  // Set up pre-processing if wanted
  //strategy.setupPreProcessing(1);
  model.setStrategy(strategy);

  ////////////////////////////////////////////////////////////
  // Do complete search
  
  model.branchAndBound();

  std::cout<<mpsFileName<<" took "<<CoinCpuTime()-time1<<" seconds, "
	   <<model.getNodeCount()<<" nodes with objective "
	   <<model.getObjValue()
	   <<(!model.status() ? " Finished" : " Not finished")
	   <<std::endl;

  // Print more statistics
  std::cout<<"Cuts at root node changed objective from "<<model.getContinuousObjective()
	   <<" to "<<model.rootObjectiveAfterCuts()<<std::endl;

  for (int iGenerator=0;iGenerator<model.numberCutGenerators();iGenerator++) {
    CbcCutGenerator * generator = model.cutGenerator(iGenerator);
    std::cout<<generator->cutGeneratorName()<<" was tried "
	     <<generator->numberTimesEntered()<<" times and created "
	     <<generator->numberCutsInTotal()<<" cuts of which "
	     <<generator->numberCutsActive()<<" were active after adding rounds of cuts";
    if (generator->timing())
      std::cout<<" ( "<<generator->timeInCutGenerator()<<" seconds)"<<std::endl;
    else
      std::cout<<std::endl;
  }
  // Print solution if finished - we can't get names from Osi! - so get from OsiClp

  if (model.getMinimizationObjValue()<1.0e50) {
    OsiSolverInterface * solver = model.solver();
    int numberColumns = solver->getNumCols();
    
    const double * solution = solver->getColSolution();

    // Get names from solver1 (as OsiSolverInterface may lose)
    std::vector<std::string> columnNames = *solver1.getModelPtr()->columnNames();
    
    int iColumn;
    std::cout<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);
    
    std::cout<<"--------------------------------------"<<std::endl;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double value=solution[iColumn];
      if (fabs(value)>1.0e-7&&solver->isInteger(iColumn)) 
	std::cout<<std::setw(6)<<iColumn<<" "
                 <<columnNames[iColumn]<<" "
                 <<value<<std::endl;
    }
    std::cout<<"--------------------------------------"<<std::endl;
  
    std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);
  }
  return 0;
}    
