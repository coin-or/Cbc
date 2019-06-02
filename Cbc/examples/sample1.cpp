// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>

#include "CbcConfig.h"
#include "CoinPragma.hpp"

// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"

#include "OsiClpSolverInterface.hpp"
#include "ClpPresolve.hpp"
#include "CbcCompareUser.hpp"
#include "CglProbing.hpp"

//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif
// Time

#include <time.h>
#if !defined(_MSC_VER)
#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>
#endif
static double cpuTime()
{
  double cpu_temp;
#if defined(_MSC_VER)
  unsigned int ticksnow; /* clock_t is same as int */

  ticksnow = (unsigned int)clock();

  cpu_temp = (double)((double)ticksnow / CLOCKS_PER_SEC);
#else
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  cpu_temp = (double)usage.ru_utime.tv_sec;
  cpu_temp += 1.0e-6 * ((double)usage.ru_utime.tv_usec);
#endif
  return cpu_temp;
}

/************************************************************************

This main program reads in an integer model from an mps file.

It then sets up some Cgl cut generators and calls branch and cut.

Branching is simple binary branching on integer variables.

Node selection is depth first until first solution is found and then
based on objective and number of unsatisfied integer variables.

Variable branching selection is on maximum minimum-of-up-down change
after strong branching on 5 variables closest to 0.5.

A simple rounding heuristic is used.

Any cut generators based on Cgl can be added in same way

You may also wish to look at CbcModel.hpp


************************************************************************/

int main(int argc, const char *argv[])
{

  // Define your favorite OsiSolver

  ClpSimplex simplex;
  double time0 = cpuTime();
  double time1 = time0;
  double time2;

  // Read in model using argv[1]
  // and assert that it is a clean model

  if (argc <= 1) {
    printf("using %s <modelfile>\n", argv[0]);
    return 1;
  }
  int numMpsReadErrors = simplex.readMps(argv[1], "");
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }
  time2 = cpuTime();
  std::cout << "Input took " << time2 - time1 << " seconds" << std::endl;
  ;
  time1 = time2;
  // Should work with OsiPresolve but not sure - so this is complicated
  ClpPresolve pinfo;
  ClpSimplex *simplex2 = pinfo.presolvedModel(simplex, 1.0e-8);
  time2 = cpuTime();
  std::cout << "Presolve took " << time2 - time1 << " seconds" << std::endl;
  ;
  time1 = time2;
  if (!simplex2 || !simplex2->integerInformation()) {
    std::cout << "Please use a feasible problem which has integers after presolve" << std::endl;
    exit(77);
  }
  OsiClpSolverInterface solver1(simplex2);
  solver1.writeMps("bad2");
  // Do initial solve to continuous
  solver1.initialSolve();
  time2 = cpuTime();
  std::cout << "Continuous solve took " << time2 - time1 << " seconds" << std::endl;
  ;
  time1 = time2;
  solver1.messageHandler()->setLogLevel(0);
  CbcModel model(solver1);
  // Definition of node choice
  CbcCompareUser compare;
  model.setNodeComparison(compare);

  // Maybe probing due to very large coefficients

  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(3);
  generator1.setMaxProbe(100);
  generator1.setMaxLook(50);
  generator1.setRowCuts(3);

  // Add in generators
  // model.addCutGenerator(&generator1,-1,"Probing");
  // Switch off strong branching if wanted
  model.setNumberStrong(0);
  model.solver()->setIntParam(OsiMaxNumIterationHotStart, 100);
  //model.solver()->setHintParam(OsiDoScale,false,OsiHintTry);
  // Switch off most output
  if (model.getNumCols() < 3000) {
    model.messageHandler()->setLogLevel(1);
    model.solver()->messageHandler()->setLogLevel(0);
  } else {
    model.messageHandler()->setLogLevel(2);
    model.solver()->messageHandler()->setLogLevel(1);
  }

  // Do complete search

  model.branchAndBound();
  time2 = cpuTime();
  std::cout << "Search took " << time2 - time1 << " seconds" << std::endl;
  // as we made such a mess of presolve lets be safe
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(model.solver());
  assert(clpSolver);
  ClpSimplex *clp = clpSolver->getModelPtr();
  *simplex2 = *clp;
  pinfo.postsolve(true);
  time1 = time2;
  // Fix all integers
  const int *original = pinfo.originalColumns();
  double *lower2 = simplex2->columnLower();
  double *upper2 = simplex2->columnUpper();
  const char *info2 = simplex2->integerInformation();
  double *lower = simplex.columnLower();
  double *upper = simplex.columnUpper();
  int i;
  for (i = 0; i < simplex2->numberColumns(); i++) {
    if (info2[i]) {
      int iSeq = original[i];
      upper[iSeq] = upper2[i];
      lower[iSeq] = lower2[i];
    }
  }

  simplex.primal();
  time2 = cpuTime();
  std::cout << "Cleanup took " << time2 - time1 << " seconds" << std::endl;
  ;
  std::cout << "Total time " << time2 - time0 << " seconds" << std::endl;
  ;
  return 0;
}
