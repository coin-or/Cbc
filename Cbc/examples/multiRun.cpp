// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iomanip>

#include "CbcModel.hpp"
#include "OsiClpSolverInterface.hpp"

#include  "CoinTime.hpp"

//#############################################################################


/************************************************************************

This main program reads in an integer model from an mps file.
It then tries four variants of brach and bound - this is to find leaks in callCbc

************************************************************************/

int main (int argc, const char *argv[])
{

  OsiClpSolverInterface solver1;

  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName = "../../Data/Sample/p0033.mps";
  if (argc>=2) mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(),"");
  assert(numMpsReadErrors==0);
  double time1 = CoinCpuTime();

  std::string test[4]= {
    "-solve",
    "-rins on -solve",
    "-preprocess off -solve",
    "-probing on -solve"};
  // do 8 runs
  for (int i=0;i<8;i++) {
    OsiClpSolverInterface solver2 = solver1;
    int type = i/4;
    int j=i-4*type;
    if (type==0) {
      callCbc(test[j].c_str(),solver2);
      printf ("Doing char * %s - took %g seconds\n",test[j].c_str(),CoinCpuTime()-time1);
    } else {
      callCbc(test[j],solver2);
      printf ("Doing string %s - took %g seconds\n",test[j].c_str(),CoinCpuTime()-time1);
    }
    int numberColumns = solver2.getNumCols();
    
    const double * solution = solver2.getColSolution();
    
    int iColumn;
    std::cout<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);
    
    std::cout<<"--------------------------------------"<<std::endl;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double value=solution[iColumn];
      if (fabs(value)>1.0e-7&&solver2.isInteger(iColumn)) 
	std::cout<<std::setw(6)<<iColumn<<" "<<value<<std::endl;
    }
    std::cout<<"--------------------------------------"<<std::endl;
  
    std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);
    time1 = CoinCpuTime();
  }
  return 0;
}    
