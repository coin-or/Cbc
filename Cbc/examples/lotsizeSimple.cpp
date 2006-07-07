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
#include "CbcBranchLotsizeSimple.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinHelperFunctions.hpp"

// Time
#include "CoinTime.hpp"


/************************************************************************

This main program reads in an integer model from an mps file.

It then replaces all integer variables by lotsizing variables
which can take integral value

OR

it uses 0,2,4... or 0,1,3,

*************************************************************************/
int main (int argc, const char *argv[])
{
  
  // Define your favorite OsiSolver
  
  OsiClpSolverInterface solver1;
  
  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName = "miplib3/10teams";
  if (argc>=2) mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(),"");
  assert(numMpsReadErrors==0);
  
  int iColumn;
  int numberColumns = solver1.getNumCols();
  int numberLot=0;
  char * mark = new char[numberColumns];
  int largestBound=0;
  // take off integers but find where they are
  const double * upper = solver1.getColUpper();
  const double * lower = solver1.getColLower();
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (solver1.isInteger(iColumn)) {
      solver1.setContinuous(iColumn);
      mark[iColumn]=1;
      numberLot++;
      assert (lower[iColumn]>=0.0);
      largestBound = CoinMax(largestBound,(int) upper[iColumn]);
    } else {
      mark[iColumn]=0;
    }
  }
  CbcModel model(solver1);
  // Do lotsizing
  CbcObject ** objects = new CbcObject * [numberLot];
  double * points = new double [largestBound+1];
  bool ordinary = argc==2;
  CoinIotaN(points,largestBound+1,0.0);
  
  numberLot=0;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (mark[iColumn]) {
      int iLo = (int) lower[iColumn];
      int iUp = (int) upper[iColumn];
      if (ordinary) {
        objects[numberLot++]= new CbcLotsizeSimple(&model,iColumn,iUp-iLo+1,
                                                   points+iLo);
      } else {
        // recreate points
        int i;
        int n=0;
        if (((iUp-iLo)&1)==0) {
          // every second
          for (i=iLo;i<=iUp;i+=2) 
            points[n++]=i;
        } else {
          // iLo,iLo+1 then every second
          points[n++]=iLo;
          for (i=iLo+1;i<=iUp;i+=2) 
            points[n++]=i;
        }
        objects[numberLot++]= new CbcLotsizeSimple(&model,iColumn,n,
                                                   points);
      }
    }
  }
  delete [] points;
  delete [] mark;
  model.addObjects(numberLot,objects);
  for (iColumn=0;iColumn<numberLot;iColumn++)
    delete objects[iColumn];
  delete [] objects;

  // Switch off most output
  model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  model.messageHandler()->setLogLevel(1);
  
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
      if (fabs(value)>1.0e-7) 
	std::cout<<std::setw(6)<<iColumn<<" "<<value<<std::endl;
    }
    std::cout<<"--------------------------------------"<<std::endl;
  
    std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);
  }
  return 0;
}    
