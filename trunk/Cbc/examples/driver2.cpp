// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iomanip>


// For Branch and bound
//#include "CbcStrategy.hpp"
#include "OsiCbcSolverInterface.hpp"

#include  "CoinTime.hpp"

//#############################################################################


/************************************************************************

This main program reads in an integer model from an mps file.
It then uses default strategy - just cuts at root node

************************************************************************/

int main (int argc, const char *argv[])
{

  // This would just do cuts at root
  // OsiCbcSolverInterface solver1;
  // This does cuts in tree and uses Clp
  CbcStrategyDefault strategy(false);
  OsiCbcSolverInterface solver1(NULL,&strategy);
  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName = "../../Data/Sample/p0033.mps";
  if (argc>=2) mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(),"");
  assert(numMpsReadErrors==0);
  double time1 = CoinCpuTime();
  // Do complete search
  
  solver1.branchAndBound();

  std::cout<<mpsFileName<<" took "<<CoinCpuTime()-time1<<" seconds, "
	   <<solver1.getNodeCount()<<" nodes with objective "
	   <<solver1.getObjValue()
	   <<(!solver1.status() ? " Finished" : " Not finished")
	   <<std::endl;

  // Print solution if finished - we can't get names from Osi!

  if (solver1.getObjValue()*solver1.getObjSense()<1.0e50) {
    int numberColumns = solver1.getNumCols();
    
    const double * solution = solver1.getColSolution();
    
    int iColumn;
    std::cout<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);
    
    std::cout<<"--------------------------------------"<<std::endl;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double value=solution[iColumn];
      if (fabs(value)>1.0e-7&&solver1.isInteger(iColumn)) 
	std::cout<<std::setw(6)<<iColumn<<" "<<value<<std::endl;
    }
    std::cout<<"--------------------------------------"<<std::endl;
  
    std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);
  }
  return 0;
}    
