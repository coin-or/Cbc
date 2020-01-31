// $Id: presolveBefore.cpp
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>

#include "CbcConfig.h"
#include "CoinPragma.hpp"

#include "CbcModel.hpp"

#include "OsiClpSolverInterface.hpp"
#include "ClpPresolve.hpp"
#include "CoinTime.hpp"

//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif

/************************************************************************

This main program reads in an integer model from an mps or lp file.

It presolves problem (to file), passes that to standalone solver and then
postsolves.

************************************************************************/

int main(int argc, const char *argv[])
{

  double time0 = CoinCpuTime(), time1, time2;
  // Define your favorite OsiSolver

  ClpSimplex * simplex = new ClpSimplex();;
  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName;
  if (argc < 2) {
    fprintf(stderr, "Do not know where to find sample MPS files.\n");
    exit(1);
  }
  if (argc >= 2)
    mpsFileName = argv[1];
  int numMpsReadErrors;
  if (strstr(mpsFileName.c_str(), ".mps"))
    numMpsReadErrors = simplex->readMps(mpsFileName.c_str(), "");
  else
    numMpsReadErrors = simplex->readLp(mpsFileName.c_str());
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }
  ClpPresolve pinfo;
  // dont allow some things
  pinfo.setDoDupcol(false);
  ClpSimplex *simplexA = pinfo.presolvedModel(*simplex, 1.0e-8);
  if (!simplexA) {
    std::cout << "Problem is not feasible" << std::endl;
    exit(77);
  }
  //#define SAVE_MEMORY
#ifdef SAVE_MEMORY
  delete simplex;
#endif
  OsiClpSolverInterface solver1(simplexA);
  // Do initial solve to continuous
  solver1.initialSolve();
  time1 = CoinCpuTime();
  std::cout << "Initialization " << time1 - time0 << " seconds" << std::endl;
  CbcModel model(solver1);
  // initialize
  CbcMain0(model);
  /* Now go into code for standalone solver
     Could copy arguments and add -quit at end to be safe
     but this will do
  */
  if (argc > 2) {
    CbcMain1(argc - 1, argv + 1, model);
  } else {
    const char *argv2[] = { "presolveBefore", "-solve", "-quit" };
    CbcMain1(3, argv2, model);
  }

  if (!simplexA->problemStatus()) {
    std::cout << "Objective value " << simplexA->objectiveValue() << std::endl;
  } else {
    std::cout << "Infeasible!" << std::endl;
  }
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(model.solver());
  assert(clpSolver);
  ClpSimplex *simplex2 = clpSolver->getModelPtr();
  time1 = CoinCpuTime();

#ifdef SAVE_MEMORY
  simplex = new ClpSimplex();
  if (strstr(mpsFileName.c_str(), ".mps"))
    simplex->readMps(mpsFileName.c_str(), "");
  else
    simplex->readLp(mpsFileName.c_str());
  pinfo.setOriginalModel(simplex);
#endif
  pinfo.postsolve(true);
  // Fix all integers (if there are any)
  const char *info2 = simplex2->integerInformation();
  if (info2) {
    const int *original = pinfo.originalColumns();
    double *lower2 = simplex2->columnLower();
    double *upper2 = simplex2->columnUpper();
    double *lower = simplex->columnLower();
    double *upper = simplex->columnUpper();
    int i;
    for (i = 0; i < simplex2->numberColumns(); i++) {
      if (info2[i]) {
	int iSeq = original[i];
	upper[iSeq] = upper2[i];
      lower[iSeq] = lower2[i];
      }
    }
  }
  delete simplexA;
  simplex->initialSolve();
  time2 = CoinCpuTime();
  std::cout << "Cleanup took " << time2 - time1 << " seconds" << std::endl;
  std::cout << "Total time " << time2 - time0 << " seconds" << std::endl;
  delete simplex;
  return 0;
}
