// $Id$
// Copyright (C) 2004, International Business Machines
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
#include "CbcBranchActual.hpp"
#include "CbcCompareUser.hpp"
#include "CoinTime.hpp"
#include "OsiClpSolverInterface.hpp"

//#############################################################################


/************************************************************************

This main program reads in an SOSr model (ltw) from an mps file.

It then solves it three ways :-

a) As normal
b) SOS 1
c) SOS 2(so answer will be different)

************************************************************************/

int main (int argc, const char *argv[])
{

  // Define your favorite OsiSolver
  
  OsiClpSolverInterface solver1;
  //solver1.messageHandler()->setLogLevel(0);
  CbcModel model(solver1);
  model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);

  // Read in ltw.mps
  // and assert that it is a clean model
  int numMpsReadErrors = model.solver()->readMps("./ltw.mps","");
  assert(numMpsReadErrors==0);

  // Definition of node choice
  CbcCompareUser compare;
  compare.setWeight(0.0);
  model.setNodeComparison(compare);
  // Reduce output
  model.messageHandler()->setLogLevel(1);
  // Get branching messages
  model.messageHandler()->setLogLevel(3);

  // Do initial solve to continuous
  model.initialSolve();

  // Save model
  CbcModel model2 = model;
  int numberColumns=model.getNumCols();
  int numberIntegers = 0;
  int * integerVariable = new int[numberColumns];
  int i;
  for ( i=0;i<numberColumns;i++) {
    if (model.isInteger(i)) {
      integerVariable[numberIntegers++]=i;
    }
  }

  
  double time1 = CoinCpuTime() ;

  model.branchAndBound();

  std::cout<<"ltw.mps"<<" took "<<CoinCpuTime()-time1<<" seconds, "
	   <<model.getNodeCount()<<" nodes with objective "
	   <<model.getObjValue()
	   <<(!model.status() ? " Finished" : " Not finished")
	   <<std::endl;

  const double * solution = model.solver()->getColSolution();
  
  std::cout<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);
  
  std::cout<<"--------------------------------------"<<std::endl;
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    double value=solution[iColumn];
    if (fabs(value)>1.0e-7) 
      std::cout<<std::setw(6)<<iColumn<<" "<<value<<std::endl;
  }
  std::cout<<"--------------------------------------"<<std::endl;
  
  std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);

  // Restore model
  model = model2;

  // Now use SOS1
  int numberSets=8;
  int which[28]={20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,
		 39,40,41,42,43,44,45,46,47};
  double weights[]={1.0,2.0,3.0,4.0,5.0};
  int starts[]={0,2,4,6,8,13,18,23,28};
  CbcObject ** objects = new CbcObject * [numberSets];
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    // Stop being integer
    model.solver()->setContinuous(iColumn);
  }
  for (i=0;i<numberSets;i++) {
    objects[i]= new CbcSOS(&model,starts[i+1]-starts[i],which+starts[i],
			   weights,i);
  }
  model.addObjects(numberSets,objects);
  for (i=0;i<numberSets;i++)
    delete objects[i];
  delete [] objects;

  time1 = CoinCpuTime() ;

  model.branchAndBound();

  std::cout<<"ltw.mps"<<" took "<<CoinCpuTime()-time1<<" seconds, "
	   <<model.getNodeCount()<<" nodes with objective "
	   <<model.getObjValue()
	   <<(!model.status() ? " Finished" : " Not finished")
	   <<std::endl;

  solution = model.solver()->getColSolution();
  
  std::cout<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);
  
  std::cout<<"--------------------------------------"<<std::endl;
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    double value=solution[iColumn];
    if (fabs(value)>1.0e-7) 
      std::cout<<std::setw(6)<<iColumn<<" "<<value<<std::endl;
  }
  std::cout<<"--------------------------------------"<<std::endl;
  
  std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);


  // Restore model
  model = model2;

// Now use SOS2
  objects = new CbcObject * [numberSets];
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    // Stop being integer
    model.solver()->setContinuous(iColumn);
  }
  for (i=0;i<numberSets;i++) {
    objects[i]= new CbcSOS(&model,starts[i+1]-starts[i],which+starts[i],
			   weights,i,2);
  }
  model.addObjects(numberSets,objects);
  for (i=0;i<numberSets;i++)
    delete objects[i];
  delete [] objects;

  time1 = CoinCpuTime() ;

  model.branchAndBound();

  std::cout<<"ltw.mps"<<" took "<<CoinCpuTime()-time1<<" seconds, "
	   <<model.getNodeCount()<<" nodes with objective "
	   <<model.getObjValue()
	   <<(!model.status() ? " Finished" : " Not finished")
	   <<std::endl;

  solution = model.solver()->getColSolution();
  
  std::cout<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);
  
  std::cout<<"--------------------------------------"<<std::endl;

  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    double value=solution[iColumn];
    if (fabs(value)>1.0e-7) 
      std::cout<<std::setw(6)<<iColumn<<" "<<value
	       <<std::endl;
  }
  std::cout<<"--------------------------------------"<<std::endl;
  
  std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);


  delete [] integerVariable;
  return 0;
}    
