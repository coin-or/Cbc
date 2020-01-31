// $Id: parallel.cpp 1902 2013-04-10 16:58:16Z stefan $
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
This is like driver4 but executes in parallel
First it reads in a model from an mps file
Then it initializes three integer models with cbc defaults
Then it calls CbcMain0/1 using threads passing parameters
Finally it prints solution (just for first model)

All models have same parameters unless "-switch" is found so --

miplib/p0033 -solve -switch -heuristic off -solve

would solve first model with heuristics and subsequent ones without

This could be used to try different ideas OR on different models

NOTE -
The minimum has been done to make thread safe so
no interrupts
-quit is added to make sure just reads from argv
*/
/*************************************************************************/
#define USE_PTHREAD
#ifdef USE_PTHREAD
#include <pthread.h>
#endif
/* Return non-zero to return quickly */
static int callBack(CbcModel *model, int whereFrom)
{
  int returnCode = 0;
  switch (whereFrom) {
  case 1:
  case 2:
    if (!model->status() && model->secondaryStatus())
      returnCode = 1;
    break;
  case 3: {
    //CbcCompareUser compare;
    //model->setNodeComparison(compare);
  } break;
  case 4:
    // If not good enough could skip postprocessing
    break;
  case 5:
    break;
  default:
    abort();
  }
  return returnCode;
}
// For threads
typedef struct {
  CbcModel *model;
  CbcSolverUsefulData *data;
  int argc;
  char **argv;
} threadStuff;
static void *doThread(void *voidInfo)
{
  threadStuff *stuff = reinterpret_cast< threadStuff * >(voidInfo);
  CbcModel *model = stuff->model;
  CbcMain0(*model, *stuff->data);
  // Now go into code for standalone solver
  CbcMain1(stuff->argc, const_cast< const char ** >(stuff->argv),
    *model, callBack, *stuff->data);
  return NULL;
}
int main(int argc, const char *argv[])
{
  // Number of models to do at once
  int numberModels = 3;
  // Stuff for each model
  CbcModel *allModels = new CbcModel[numberModels];
  OsiClpSolverInterface *allSolvers = new OsiClpSolverInterface[numberModels];
  CbcSolverUsefulData *data = new CbcSolverUsefulData[numberModels];
  threadStuff *allData = new threadStuff[numberModels];
  // First populate first model
  /* in a real application models would be different */
  OsiClpSolverInterface &solver1 = allSolvers[0];
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
  // create models
  for (int iModel = 1; iModel < numberModels; iModel++) {
    allSolvers[iModel] = solver1;
  }
  // now create CbcModels and copy parameters (up to end or -switch)
  // Pass to Cbc initialize defaults
  CbcModel modelA(solver1);
  int lastArgPosition = 2;
  int lastArgc = 0;
  for (int iModel = 0; iModel < numberModels; iModel++) {
    allModels[iModel] = CbcModel(allSolvers[iModel]);
    // default does NOT allow interrupts
    // allow printing
    data[iModel].noPrinting_ = false;
    allData[iModel].model = allModels + iModel;
    allData[iModel].data = data + iModel;
    // See how many parameters
    int endArgc = lastArgPosition;
    int thisArgc = lastArgc;
    for (; endArgc < argc; endArgc++) {
      if (!strcmp(argv[endArgc], "-switch"))
        break;
    }
    thisArgc = endArgc - lastArgPosition;
    if (thisArgc > 0)
      lastArgc = thisArgc;
    else
      thisArgc = lastArgc;
    // allow extra for -quit
    char **thisArgv = new char *[thisArgc + 2];
    thisArgv[0] = strdup(argv[0]);
    int put = 1;
    for (int iArgc = lastArgPosition; iArgc < lastArgPosition + thisArgc; iArgc++)
      thisArgv[put++] = strdup(argv[iArgc]);
    // add -quit
    thisArgv[put++] = strdup("-quit");
    allData[iModel].argc = put;
    allData[iModel].argv = thisArgv;
    if (endArgc < argc)
      lastArgPosition = endArgc + 1;
  }
#ifdef USE_PTHREAD
  pthread_t *threadId = new pthread_t[numberModels];
  // solve
  for (int iModel = 0; iModel < numberModels; iModel++) {
    pthread_create(threadId + iModel, NULL,
      doThread,
      allData + iModel);
  }
  // wait
  for (int iModel = 0; iModel < numberModels; iModel++) {
    pthread_join(threadId[iModel], NULL);
  }
#else
  for (int iModel = 0; iModel < numberModels; iModel++) {
    doThread(allData + iModel);
  }
#endif
  // Just print first one
  CbcModel *model = allModels;
  // Solver was cloned so get current copy
  OsiSolverInterface *solver = model->solver();
  // Print solution if finished (could get from model->bestSolution() as well

  if (model->bestSolution()) {

    const double *solution = solver->getColSolution();

    int iColumn;
    int numberColumns = solver->getNumCols();
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
  } else {
    std::cout << " No solution!" << std::endl;
  }
  delete[] allModels;
  delete[] allSolvers;
  delete[] data;
  for (int iModel = 0; iModel < numberModels; iModel++) {
    char **argv = allData[iModel].argv;
    int argc = allData[iModel].argc;
    for (int i = 0; i < argc; i++)
      free(argv[i]);
    delete[] argv;
  }
  delete[] allData;
#ifdef USE_PTHREAD
  delete[] threadId;
#endif
  return 0;
}
