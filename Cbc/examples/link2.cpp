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
#include "CbcCutGenerator.hpp"
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

Having variables integer would probably make easier
************************************************************************/

int main (int argc, const char *argv[])
{


  // Create model 
  CoinModel build;
  // size
  int numberX=20;
  int numberConstraints=1000;
  // variables are x sub j and w sub ij which is x sub j * x sub i
  char name[20];
  char product[20];
  int i;
  // pointers to W (-1 if not existing)
  int ** whichW = new int * [numberX];
  for (i=0;i<numberX;i++) {
    int * sequenceW = new int[numberX];
    int j;
    for (j=0;j<numberX;j++)
      sequenceW[j]=-1;
    whichW[i]=sequenceW;
  }
  /* now decide on W.
     When both x and y are continuous then the default method for x*y will branch on x
     until it is fixed and x*y can be modeled exactly.  In this case it is better not to
     use a triangular version for W but to randomize so each variable is the "x" one
     about the same number of times.
  */
  int neededW=(numberX*(numberX+1))/2;
  int * wI = new int[neededW];
  int * wJ = new int[neededW];
  while (neededW) {
    bool gotColumn=false;
    int iColumn=-1;
    int jColumn=-1;
    while(!gotColumn) {
      iColumn=(int) floor(CoinDrand48()*numberX);
      jColumn=(int) floor(CoinDrand48()*numberX);
      if (iColumn==numberX||jColumn==numberX)
	continue;
      int * sequenceW = whichW[iColumn];
      if (sequenceW[jColumn]==-1) {
	// available
	gotColumn=true;
	break;
      }
    }
    int * sequenceW = whichW[iColumn];
    sequenceW[jColumn]=1;
    neededW--;
    if (iColumn!=jColumn) {
      int * sequenceW = whichW[jColumn];
      sequenceW[iColumn]=-2;
    }
  }
  // build columns
  for (i=0;i<numberX;i++) {
    sprintf(name,"X_%d",i);
    // If we know tighter bounds then better
    build.addColumn(0,NULL,NULL,0.0,1.0,0.0,name);
  }
  // now do W
  int numberW=0;
  for (i=0;i<numberX;i++) {
    int * sequenceW = whichW[i];
    for (int j=0;j<numberX;j++) {
      if (sequenceW[j]>0) {
	sequenceW[j]=numberX+numberW;
	assert (sequenceW[j]==build.numberColumns());
	sprintf(name,"W_%d_%d",i,j);
	// minimize sums of squares of x so ..
	double value = (i==j) ? 1.0 : 0.0;
	build.addColumn(0,NULL,NULL,0.0,1.0,value,name);
	wI[numberW]=i;
	wJ[numberW]=j;
	numberW++;
      }
    }
  }
  int numberColumns = build.numberColumns();
  int * column = new int [numberColumns];
  double * element = new double [numberColumns];
  // first row says sum of x >= 1.0 *** no make equality?
  for (i=0;i<numberX;i++) {
    column[i]=i;
    element[i]=1.0;
  }
  // Should be stronger if we make sum x == 1 rather >= 1
#define EQUALITY
#ifndef EQUALITY
  build.addRow(numberX,column,element,1.0,COIN_DBL_MAX);
#else
  build.addRow(numberX,column,element,1.0,1.0);
#endif
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
  /* Optional linear constraints 
     These are theoretically redundant but may tighten relaxation
  */
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
#ifndef EQUALITY
    build.addRow(n,column,element,0.0,COIN_DBL_MAX);
#else
    build.addRow(n,column,element,0.0,0.0);
#endif
  }
#endif
  // Now add in random quadratics - no linear term but that would be easy (X variables)
  double rhs=0.5; // So constraints do something
  if (numberX==8)
    rhs=0.366;
  else if (numberX==20)
    rhs=0.25; // so will be quickish!
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
    // make sure last one feasible at 1.0
    if (element[n-1]>rhs)
      element[n-1]=0.9*rhs;
    build.addRow(n,column,element,-COIN_DBL_MAX,rhs);
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
    // make sure last one feasible at 1.0
    if (element[n-1]>rhs)
      element[n-1]=0.9*rhs;
    storedCuts.addCut(-COIN_DBL_MAX,rhs,n,column,element);
    build2.addRow(n,column,element,-COIN_DBL_MAX,rhs);
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
  solver1.setDefaultMeshSize(0.001);
  // Use coarse grid so we can solve quickly
  solver1.setDefaultMeshSize(0.01);
  solver1.setDefaultMeshSize(0.001);
  //solver1.setDefaultMeshSize(0.0001);
  solver1.setIntegerPriority(1000);
  solver1.setBiLinearPriority(10000);
  solver1.load(build,true,2);
  // Add more objects with even coarser grid
  solver1.setBiLinearPriorities(10,0.2);
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
  model.addCutGenerator(&storedCuts,1,"Stored cuts",true,false);
  model.cutGenerator(0)->setMustCallAgain(true);;
  //model.addCutGenerator(&storedCuts,1,"Stored cuts",false,true);
  //model.setMaximumCutPasses(4);
#endif
  //model.setNumberThreads(4);
  model.branchAndBound();

  std::cout<<"took "<<CoinCpuTime()-time1<<" seconds, "
	   <<model.getNodeCount()<<" nodes with objective "
	   <<model.getObjValue()
	   <<(!model.status() ? " Finished" : " Not finished")
	   <<std::endl;

  int numberGenerators = model.numberCutGenerators();
  for (int iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
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

  if (model.getMinimizationObjValue()<1.0e50) {
    
    const double * solution = model.bestSolution();
    double obj=0.0;
    for (i=0;i<numberX;i++) {
      printf("%d %s %g\n",i,build.columnName(i),solution[i]);
      obj += solution[i]*solution[i];
    }
    printf("true objective %g\n",obj);
    for (;i<numberColumns;i++) {
      int iX = wI[i-numberX];
      int jX = wJ[i-numberX];
      printf("%d %s %g - should be %g\n",i,build.columnName(i),solution[i],solution[iX]*solution[jX]);
    }
#ifdef USE_CUTS
    double worstViolation=0.0;
    double sumViolation=0.0;
    for (i=0;i<numberConstraints;i++) {
      const OsiRowCut * cut = storedCuts.rowCutPointer(i);
      const CoinPackedVector row = cut->row();
      int n = row.getNumElements();
      const int * column = row.getIndices();
      const double * element = row.getElements();
      double sum=0.0;
      double sumApprox=0.0;
      for (int j=0;j<n;j++) {
	int iColumn = column[j];
	sumApprox += solution[iColumn]*element[j];
	int iX = wI[iColumn-numberX];
	int jX = wJ[iColumn-numberX];
	sum += solution[iX]*solution[jX]*element[j];
      }
      if (sum>cut->ub()+1.0e-6) 
	printf("*** %d %g (approx %g) > %g\n",i,sum,sumApprox,cut->ub());
      else if (sum>cut->ub()-1.0e-3) 
	printf("    %d %g (approx %g) =~ %g\n",i,sum,sumApprox,cut->ub());
      if (sum>cut->ub()) {
	double violation = sum-cut->ub();
	sumViolation += violation;
	worstViolation = CoinMax(worstViolation,violation);
      }
    }
    printf("Total violation %g, worst was %g\n",sumViolation,worstViolation);
#endif
  }
  // deletes
  for (i=0;i<numberX;i++) {
    delete [] whichW[i];
  }
  delete [] whichW;
  delete [] wI;
  delete [] wJ;
  return 0;
}    
