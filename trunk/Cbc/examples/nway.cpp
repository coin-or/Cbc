// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iomanip>


// For Branch and bound
#include "CbcModel.hpp"
#include "CbcBranchActual.hpp"
#include "OsiClpSolverInterface.hpp"

// Time
#include "CoinTime.hpp"


/************************************************************************

This main program reads in an integer model from an mps file.

It then tries to find SOS structure and solves using N-way variables

*************************************************************************/
int main (int argc, const char *argv[])
{
  
  // Define your favorite OsiSolver
  
  OsiClpSolverInterface solver1;
  
  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName = "../../Data/miplib3/10teams";
  if (argc>=2) mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(),"");
  assert(numMpsReadErrors==0);
  
  int iRow, iColumn;
  int numberColumns = solver1.getNumCols();
  int numberRows = solver1.getNumRows();
  // get row copy
  const CoinPackedMatrix * matrix = solver1.getMatrixByRow();
  const double * element = matrix->getElements();
  const int * column = matrix->getIndices();
  const CoinBigIndex * rowStart = matrix->getVectorStarts();
  const int * rowLength = matrix->getVectorLengths();
  const double * rowLower = solver1.getRowLower();
  const double * rowUpper = solver1.getRowUpper();
  const double * columnLower = solver1.getColLower();
  
  // Look for possible SOS
  int numberSOS=0;
  int * mark = new int[numberColumns];
  CoinFillN(mark,numberColumns,-1);
  for (iRow=0;iRow<numberRows;iRow++) {
    if (rowLower[iRow]==1.0&&rowUpper[iRow]==1.0) {
      bool goodRow=true;
      for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
        int iColumn = column[j];
        if (element[j]!=1.0||!solver1.isInteger(iColumn)||mark[iColumn]>=0||columnLower[iColumn]) {
          goodRow=false;
          break;
        }
      }
      if (goodRow) {
        // mark all
        for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
          int iColumn = column[j];
          mark[iColumn]=numberSOS;
        }
        numberSOS++;
      }
    }
  }
  std::cout<<numberSOS<<" SOS"<<std::endl;
  if (!numberSOS)
    return 0;
  CbcModel model(solver1);
  // Do sets and priorities
  CbcObject ** objects = new CbcObject * [numberSOS];
  int numberIntegers = model.numberIntegers();
  /* model may not have created objects
     If none then create
  */
  if (!numberIntegers||!model.numberObjects()) {
    model.findIntegers(true);
    numberIntegers = model.numberIntegers();
  }
  int * priority = new int[numberSOS];
  // Set SOS priorities high
  CoinFillN(priority,numberSOS,1);
  // Set up SOS
  int * which = new int[numberColumns];
  for (int iSOS =0;iSOS<numberSOS;iSOS++) {
    int n=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (mark[iColumn]==iSOS)
        which[n++]=iColumn;
    }
    // NULL uses 0,1,2 .. as weights
    objects[iSOS]= new CbcNWay(&model,n,which,iSOS);
  }
  delete [] mark;
  delete [] which;
  model.addObjects(numberSOS,objects);
  for (iColumn=0;iColumn<numberSOS;iColumn++)
    delete objects[iColumn];
  delete [] objects;
  model.passInPriorities(priority,true);
  delete [] priority;

  // If time is given then stop after that number of minutes
  if (argc>2) {
    int minutes = atoi(argv[2]);
    std::cout<<"Stopping after "<<minutes<<" minutes"<<std::endl;
    assert (minutes>=0);
    model.setDblParam(CbcModel::CbcMaximumSeconds,60.0*minutes);
  }
  // Switch off most output
  model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  if (model.getNumCols()<3000) {
    model.messageHandler()->setLogLevel(1);
    //model.solver()->messageHandler()->setLogLevel(0);
  } else {
    model.messageHandler()->setLogLevel(2);
    //model.solver()->messageHandler()->setLogLevel(1);
  }
  //model.messageHandler()->setLogLevel(1);
  
  double time1 = CoinCpuTime();

  // Do complete search
  
  model.branchAndBound();

  std::cout<<mpsFileName<<" took "<<CoinCpuTime()-time1<<" seconds, "
	   <<model.getNodeCount()<<" nodes with objective "
	   <<model.getObjValue()
	   <<(!model.status() ? " Finished" : " Not finished")
	   <<std::endl;

  // Print solution - we can't get names from Osi!

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
