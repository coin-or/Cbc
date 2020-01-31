// $Id: driver4.cpp 1898 2013-04-09 18:06:04Z stefan $
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
#include "CbcModel.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcSolver.hpp"

#include "CoinTime.hpp"

//#############################################################################

/************************************************************************

This main program shows how to take advantage of the standalone cbc in your program,
while still making major modifications.
First it reads in an integer model from an mps file
Then it initializes the integer model with cbc defaults
Then it calls CbcMain1 passing all parameters apart from first but with callBack to modify stuff
Finally it prints solution
*/
int main(int argc, const char *argv[])
{

  OsiClpSolverInterface solver1;
  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName;
#if defined(SAMPLEDIR)
  mpsFileName = SAMPLEDIR "/p0033.mps";
#else
  if (argc < 2) {
    fprintf(stderr, "Do not know where to find sample MPS files.\n");
    exit(1);
  }
#endif
  if (argc >= 2)
    mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(), "");
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }
  // Tell solver to return fast if presolve or initial solve infeasible
  solver1.getModelPtr()->setMoreSpecialOptions(3);

  // Pass to Cbc initialize defaults
  CbcModel modelA(solver1);
  CbcModel *model = &modelA;
  CbcMain0(modelA);
  modelA.setMaximumSavedSolutions(5);
  /* Now go into code for standalone solver
     Could copy arguments and add -quit at end to be safe
     but this will do
  */
  if (argc > 2) {
    CbcMain1(argc - 1, argv + 1, modelA);
  } else {
    const char *argv2[] = { "driver4", "-solve", "-quit" };
    CbcMain1(3, argv2, modelA);
  }
  // Solver was cloned so get current copy
  OsiSolverInterface *solver = model->solver();
  // Print solution if finished (could get from model->bestSolution() as well

  if (model->bestSolution()) {

    int numberColumns = solver->getNumCols();
    std::cout << std::setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14);

    // do solutions
    for (int solutionNumber = 0; solutionNumber < modelA.numberSavedSolutions();
         solutionNumber++) {
      const double *solution = modelA.savedSolution(solutionNumber);
      std::cout << "-------------------------------------- solution "
                << solutionNumber + 1 << " objective " << modelA.savedSolutionObjective(solutionNumber)
                << std::endl;
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = solution[iColumn];
        if (fabs(value) > 1.0e-7 && solver1.isInteger(iColumn))
          std::cout << std::setw(6) << iColumn << " " << std::setw(8) << setiosflags(std::ios::left) << solver1.getModelPtr()->columnName(iColumn)
                    << resetiosflags(std::ios::adjustfield) << std::setw(14) << " " << value << std::endl;
      }
      std::cout << "--------------------------------------" << std::endl;
    }
    std::cout << std::resetiosflags(std::ios::fixed | std::ios::showpoint | std::ios::scientific);
  } else {
    std::cout << " No solution!" << std::endl;
  }
  return 0;
}
