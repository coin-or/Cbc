// Copyright (C) 2005, International Business Machines
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
#include "CbcBranchActual.hpp"
#include "OsiChooseVariable.hpp"
#include "CoinModel.hpp"
// For Linked Ordered Sets
#include "CbcLinked.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CglStored.hpp"

#include  "CoinTime.hpp"
#include  "CoinHelperFunctions.hpp"


/************************************************************************
This shows how to get a global solution to a non-convex quadratically 
constrained problem.

Having variables integer would make easier
************************************************************************/

int main (int argc, const char *argv[])
{


  // Create model 
  CoinModel build;
  // size
  int numberX=5;
  int numberConstraints=20;
  // variables are x sub j and w sub ij which is x sub j * x sub i
  char name[20];
  char product[20];
  int i;
  // pointers to W (-1 if not existing)
  int ** whichW = new int * [numberX];
  int numberW=0;
  // build columns
  for (i=0;i<numberX;i++) {
    sprintf(name,"X_%d",i);
    build.addColumn(0,NULL,NULL,0.0,1.0,0.0,name);
  }
  // now W
  for (i=0;i<numberX;i++) {
    // and do pointers to W
    int * sequenceW = new int[numberX];
    int j;
    for (j=0;j<i;j++)
      sequenceW[j]=-1;
    for (;j<numberX;j++) {
      sequenceW[j]=numberX+numberW;
      assert (sequenceW[j]==build.numberColumns());
      sprintf(name,"W_%d_%d",i,j);
      // minimize sums of squares
      double value = (i==j) ? 1.0 : 0.0;
      build.addColumn(0,NULL,NULL,0.0,1.0,value,name);
      numberW++;
    }
    whichW[i]=sequenceW;
  }
  int numberColumns = build.numberColumns();
  int * column = new int [numberColumns];
  double * element = new double [numberColumns];
  // first row says sum of x >= 1.0
  for (i=0;i<numberX;i++) {
    column[i]=i;
    element[i]=1.0;
  }
  build.addRow(numberX,column,element,1.0,COIN_DBL_MAX);
  // Now define W
  for (i=0;i<numberX;i++) {
    int * sequenceW = whichW[i];
    int j;
    for (j=0;j<numberX;j++) {
      int wColumn = sequenceW[j];
      if (wColumn>=0) {
	// add in w part of row
	int numberRows=build.numberRows();
	element[0]=-1.0;
	build.addRow(1,&wColumn,element,0.0,0.0);
	// add in product
	sprintf(product,"X_%d",j);
	build.setElement(numberRows,i,product);
      }
    }
  }
#if 1
  // Optional linear constraints 
  for (i=0;i<numberX;i++) {
    // We want to add in sum over j W_i_j >= X_i
    int * sequenceW = whichW[i];
    int j;
    int n=0;
    for (j=0;j<numberX;j++) {
      int wColumn = sequenceW[j];
      if (wColumn<0) {
	int * sequenceW2 = whichW[j];
	wColumn = sequenceW2[i];
	assert (wColumn>=0);
      }
      column[n]=wColumn;
      element[n++]=1.0;
    }
    // add in X part of row
    column[n]=i;
    element[n++]=-1.0;
    build.addRow(n,column,element,0.0,COIN_DBL_MAX);
  }
#endif
  // Now add in random quadratics - no linear term but that would be easy (X variables)
  double scaleFactor=0.5;
#define USE_CUTS
#ifndef USE_CUTS
  // Put constraints in matrix
  for (i=0;i<numberConstraints;i++) {
    int n=0;
    for (int j=0;j<numberW;j++) {
      // make random values clean i.e. multiple of 1.0e-5
      double value = CoinDrand48()*1.0e5;
      value = floor(value)*1.0e-5;
      if (value) {
	column[n]=j+numberX;
	element[n++]=value;
      }
    }
    build.addRow(n,column,element,-COIN_DBL_MAX,scaleFactor);
  }
#else
  // Put constraints in as cuts
  // note that if constraints may tighten bounds on X
  // So copy build and capture tight bounds
  CglStored storedCuts;
  CoinModel build2=build;
  for (i=0;i<numberConstraints;i++) {
    int n=0;
    for (int j=0;j<numberW;j++) {
      // make random values clean i.e. multiple of 1.0e-5
      double value = CoinDrand48()*1.0e5;
      value = floor(value)*1.0e-5;
      if (value) {
	column[n]=j+numberX;
	element[n++]=value;
      }
    }
    storedCuts.addCut(-COIN_DBL_MAX,scaleFactor,n,column,element);
    build2.addRow(n,column,element,-COIN_DBL_MAX,scaleFactor);
  }
  // load from coin model
  OsiSolverLink solver3;
  solver3.setDefaultMeshSize(0.001);
  solver3.setIntegerPriority(1000);
  solver3.setBiLinearPriority(10000);
  solver3.load(build2,true,2);
  // move bounds
  for (i=0;i<numberX;i++) {
    build.setColumnUpper(i,solver3.getColUpper()[i]);
  }
#endif
  delete [] column;
  delete [] element;

  // load from coin model
  OsiSolverLink solver1;
  //OsiSolverInterface * solver2 = solver1.clone();
  ///CbcModel model;
  //model.assignSolver(solver2,true);
  //OsiSolverLink * si =
  //dynamic_cast<OsiSolverLink *>(model.solver()) ;
  //assert (si != NULL);
  //solver1.setSpecialOptions2(16);
  solver1.setDefaultMeshSize(0.001);
  // need some relative granularity
  solver1.setDefaultMeshSize(0.01);
  solver1.setIntegerPriority(1000);
  solver1.setBiLinearPriority(10000);
  solver1.load(build,true,2);
  CbcModel model(solver1);
  // set more stuff
  OsiSolverLink * clpSolver = dynamic_cast< OsiSolverLink*> (model.solver());
  //ClpSimplex * lpSolver = clpSolver->getModelPtr();
  clpSolver->messageHandler()->setLogLevel(0) ;
  
  double time1 = CoinCpuTime();
  model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  model.solver()->initialSolve();
  model.solver()->writeMps("link");
  // Do complete search
  
  CbcBranchDefaultDecision decision;
  OsiChooseStrong choose(model.solver());
  choose.setNumberBeforeTrusted(model.numberBeforeTrust());
  choose.setNumberStrong(model.numberStrong());
  //choose.setShadowPriceMode(testOsiOptions);
  decision.setChooseMethod(choose);
  model.setBranchingMethod(decision);
  model.setLogLevel(1);
#ifdef USE_CUTS
  model.addCutGenerator(&storedCuts,1,"Stored cuts",false,true);
  model.setMaximumCutPasses(4);
#endif
  model.branchAndBound();

  std::cout<<"took "<<CoinCpuTime()-time1<<" seconds, "
	   <<model.getNodeCount()<<" nodes with objective "
	   <<model.getObjValue()
	   <<(!model.status() ? " Finished" : " Not finished")
	   <<std::endl;


  if (model.getMinimizationObjValue()<1.0e50) {
    
    const double * solution = model.bestSolution();
    for (i=0;i<numberColumns;i++)
      printf("%d %g\n",i,solution[i]);
  }
  // deletes
  for (i=0;i<numberX;i++) {
    delete [] whichW[i];
  }
  delete [] whichW;
  return 0;
}    
