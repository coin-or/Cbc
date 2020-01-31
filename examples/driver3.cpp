// $Id$
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
#include "CbcModel.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcBranchDynamic.hpp"
#include "CbcSolver.hpp"

#include "CoinTime.hpp"

//#############################################################################

/************************************************************************

This main program shows how to take advantage of the standalone cbc in your program.
It should perform very nearly the same as cbc  
First it reads in an integer model from an mps file and saves and strips off integer information.
Then it initializes the integer model with cbc defaults
Then it puts back integers - here you could do anything and also set parameters
Then it calls CbcMain1 passing all parameters apart from first
Finally it prints solution

************************************************************************/

static int dummyCallBack(CbcModel * /*model*/, int /*whereFrom*/)
{
  return 0;
}

int main(int argc, const char *argv[])
{

  OsiClpSolverInterface solver1;
  //#define USE_OSI_NAMES
#ifdef USE_OSI_NAMES
  // Say we are keeping names (a bit slower this way)
  solver1.setIntParam(OsiNameDiscipline, 1);
#endif
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

  // Strip off integer information and save
  int numberColumns = solver1.getNumCols();
  char *integer = new char[numberColumns];
  int i;
  for (i = 0; i < numberColumns; i++) {
    if (solver1.isInteger(i)) {
      integer[i] = 1;
      solver1.setContinuous(i);
    } else {
      integer[i] = 0;
    }
  }
  // Pass to Cbc initialize defaults
  CbcSolverUsefulData cbcData;
  CbcModel model(solver1);
  CbcMain0(model, cbcData);

  // Solve just to show there are no integers
  model.branchAndBound();
  // Set cutoff etc back in model and solver
  model.resetModel();
  // Solver was cloned so get it
  OsiSolverInterface *solver = model.solver();
  // Put back integers.  Here the user could do anything really
#define ADD_DIRECTLY
#ifndef ADD_DIRECTLY
  for (i = 0; i < numberColumns; i++) {
    if (integer[i])
      solver->setInteger(i);
  }
#else
  CbcObject **objects = new CbcObject *[numberColumns];
  int n = 0;
  for (i = 0; i < numberColumns; i++) {
    if (integer[i]) {
      CbcSimpleIntegerDynamicPseudoCost *newObject = new CbcSimpleIntegerDynamicPseudoCost(&model, i);
      objects[n++] = newObject;
    }
  }
  model.addObjects(n, objects);
  for (i = 0; i < n; i++)
    delete objects[i];
  delete[] objects;
#endif
  delete[] integer;
  /* Now go into code for standalone solver
     Could copy arguments and add -quit at end to be safe
     but this will do
  */
  if (argc > 2) {
    CbcMain1(argc - 1, argv + 1, model, dummyCallBack, cbcData);
  } else {
    const char *argv2[] = { "driver3", "-solve", "-quit" };
    CbcMain1(3, argv2, model, dummyCallBack, cbcData);
  }

  // Print solution if finished (could get from model.bestSolution() as well

  if (solver->getObjValue() * solver->getObjSense() < 1.0e50) {

    const double *solution = solver->getColSolution();

    int iColumn;
    std::cout << std::setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14);

    std::cout << "--------------------------------------" << std::endl;
#ifdef USE_OSI_NAMES

    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      if (fabs(value) > 1.0e-7 && solver->isInteger(iColumn))
        std::cout << std::setw(6) << iColumn << " " << std::setw(8) << setiosflags(std::ios::left) << solver->getColName(iColumn)
                  << resetiosflags(std::ios::adjustfield) << std::setw(14) << " " << value << std::endl;
    }
#else
    // names may not be in current solver - use original

    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      if (fabs(value) > 1.0e-7 && solver->isInteger(iColumn))
        std::cout << std::setw(6) << iColumn << " " << std::setw(8) << setiosflags(std::ios::left) << solver1.getModelPtr()->columnName(iColumn)
                  << resetiosflags(std::ios::adjustfield) << std::setw(14) << " " << value << std::endl;
    }
#endif
    std::cout << "--------------------------------------" << std::endl;

    std::cout << std::resetiosflags(std::ios::fixed | std::ios::showpoint | std::ios::scientific);
  }
  return 0;
}
