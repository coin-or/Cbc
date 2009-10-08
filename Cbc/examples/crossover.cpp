// Copyright (C) 2009, International Business Machines
// Corporation and others.  All Rights Reserved.

// Using Clp as solver
#include "OsiClpSolverInterface.hpp"
#include "ClpInterior.hpp"

/* This is to show use of crossover from a close to optimal solution
   to a basic optimal solution.

   It should not really be in Cbc/examples but it can't go in Osi/examples
   
   options defaults to 0 and basis to 1.  If file name given then 
   argv[2] can override options and argv[3] basis.

   options - 0 no presolve (use primal and dual)
             1 presolve (just use primal)
	     2 no presolve (just use primal)
   basis -   0 use all slack basis
             1 try and put some in basis
*/

int main (int argc, const char *argv[])
{
  OsiClpSolverInterface solver1;
  // Read in example model
  // and assert that it is a clean model
  std::string mpsFileName = "../../Data/Sample/p0033.mps";
  if (argc>=2) mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(),"");
  assert(numMpsReadErrors==0);

  ClpSimplex * simplex = solver1.getModelPtr();
  ClpInterior barrier;
  barrier.borrowModel(*simplex);
  barrier.primalDual();
  barrier.returnModel(*simplex);

  // Do crossover
  int options=0;
  int basis=1;
  if (argc>=3) {
    options=atoi(argv[2]);
    if (argc>=4) 
      basis=atoi(argv[3]);
    printf("options %d basis %d\n",options,basis);
  }
  //simplex->setLogLevel(4);
  solver1.crossover(options,basis);
  printf("Following should be zero iterations\n");

  solver1.resolve(); 

  return 0;
}
   
 
