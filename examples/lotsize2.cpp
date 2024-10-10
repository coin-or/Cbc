// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
#include "CbcModel.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcSolver.hpp"
#include "CbcBranchLotsize.hpp"

#include "CoinTime.hpp"

//#############################################################################

/************************************************************************

This main program shows how to take advantage of the standalone cbc in your program,
while still making major modifications.
First it reads in an integer model
Then it initializes the integer model with cbc defaults
Then it calls CbcMain1 passing all parameters apart from first two but with callBack to modify stuff
The callBack adds lotsizing fom second parameter
Initially switch off preprocessing
Finally it prints solution

************************************************************************/
/* Meaning of whereFrom:
   1 after initial solve by dualsimplex etc
   2 after preprocessing
   3 just before branchAndBound (so user can override)
   4 just after branchAndBound (before postprocessing)
   5 after postprocessing
*/
/* Meaning of model status is as normal
   status
      -1 before branchAndBound
      0 finished - check isProvenOptimal or isProvenInfeasible to see if solution found
      (or check value of best solution)
      1 stopped - on maxnodes, maxsols, maxtime
      2 difficulties so run was abandoned
      (5 event user programmed event occurred) 

      cbc secondary status of problem
        -1 unset (status_ will also be -1)
	0 search completed with solution
	1 linear relaxation not feasible (or worse than cutoff)
	2 stopped on gap
	3 stopped on nodes
	4 stopped on time
	5 stopped on user event
	6 stopped on solutions
	7 linear relaxation unbounded

   but initially check if status is 0 and secondary status is 1 -> infeasible
   or you can check solver status.
*/
static const char * lotsizingFile=NULL;
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
  case 3:
    {
      // Add in lot sizing
      /* format of file is
	 variable name,a,b,c....
	 where each entry is either a value or value~value pair.
      */
      FILE * fp = fopen(lotsizingFile,"r");
      if (!fp) {
	printf("bad filename %s\n",lotsizingFile);
	exit(2);
      }
      OsiSolverInterface * solver = model->solver();
      int numberColumns = solver->getNumCols();
      const double * lower = solver->getColLower();
      const double * upper = solver->getColUpper();
      CbcObject **objects = new CbcObject *[numberColumns];
      int numberLot = 0;
#define MAXIMUM_ENTRIES 100
      double pairs[2*MAXIMUM_ENTRIES];
      double loValue[MAXIMUM_ENTRIES];
      double upValue[MAXIMUM_ENTRIES];
      double sort[MAXIMUM_ENTRIES];
      int which[MAXIMUM_ENTRIES];
      char line[500];
      while (fgets(line,500,fp)) {
	// add terminating string
	char * newLine = strchr(line,'\n');
	if (newLine)
	  *newLine='\0';
	strcat(line,",zz,");
	// name
	char * comma = strchr(line,',');
	assert (comma);
	*comma='\0';
	std::string name=line;
	int iColumn;
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  if (name==solver->getColName(iColumn))
	    break;
	}
	if (iColumn==numberColumns) {
	  printf("variable %s not found\n",line);
	  exit(4);
	}
	char * pos = comma+1;
	int nPairs = 0;
	while (*pos!='z') {
	  comma = strchr(pos,',');
	  assert (comma);
	  *comma = '\0';
	  double lo,up;
	  if (strchr(pos,'~')) {
	    // pair of values
	    char * sep = strchr(pos,'~');
	    *sep = '\0';
	    lo = atof(pos);
	    up = atof(sep+1);
	    assert (lo<=up);
	  } else {
	    // single value
	    lo = atof(pos);
	    up = lo;
	  }
	  sort[nPairs] = lo;
	  loValue[nPairs] = lo;
	  upValue[nPairs] = up;
	  which[nPairs] = nPairs;
	  nPairs++;
	  pos = comma+1;
	}
	// sort and truncate
	CoinSort_2(sort,sort+nPairs,which);
	int n = nPairs;
	nPairs = 0;
	double lowerBound = lower[iColumn];
	double upperBound = upper[iColumn];
	for (int i=0;i<n;i++) {
	  if (upValue[i]<lowerBound)
	    continue;
	  if (loValue[i]>upperBound)
	    break;
	  loValue[i] = std::max(lowerBound,loValue[i]);
	  upValue[i] = std::min(upperBound,upValue[i]);
	  pairs[2*nPairs] = loValue[i];
	  pairs[2*nPairs+1] = upValue[i];
	  nPairs++;
	}
	// update bounds in solver
	if (pairs[0]>lower[iColumn])
	  solver->setColLower(iColumn,pairs[0]);
	if (pairs[2*nPairs-1]<upper[iColumn])
	  solver->setColUpper(iColumn,pairs[2*nPairs-1]);
	// create object
	objects[numberLot++] = new CbcLotsize(model, iColumn, nPairs, pairs, true);
      }
      int numberObjects = model->numberObjects();
      model->addObjects(numberLot, objects);
      for (int iColumn = 0; iColumn < numberLot; iColumn++)
	delete objects[iColumn];
      delete[] objects;
      fclose(fp);
      // Set priorities
      OsiObject ** objectsNow = model->objects();
      for (int i=0;i<numberObjects;i++)
	objectsNow[i]->setPriority(2000);
      for (int i=numberObjects;i<numberObjects+numberLot;i++)
	objectsNow[i]->setPriority(1000);
    }
    break;
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
  const char * fileName;
  if (argc < 3) {
    fprintf(stderr, "Do not know where to find files.\n");
    exit(1);
  }
  fileName = argv[1];
  lotsizingFile = argv[2];
  int numReadErrors;
  if (!strstr(fileName,".lp"))
    numReadErrors = solver1.readMps(fileName, "");
  else
    numReadErrors = solver1.readLp(fileName);
  if (numReadErrors != 0) {
    printf("%d errors reading file\n", numReadErrors);
    return numReadErrors;
  }
  // Tell solver to return fast if presolve or initial solve infeasible
  solver1.getModelPtr()->setMoreSpecialOptions(3);

  /* Two ways of doing this depending on whether NEW_STYLE_SOLVER defined.
     So we need pointer to model.  Old way could use modelA. rather than model->
   */
  // Messy code below copied from CbcSolver.cpp
  // Pass to Cbc initialize defaults
  CbcModel modelA(solver1);
  CbcSolverUsefulData solverData;
  CbcModel *model = &modelA;
  CbcMain0(modelA,solverData);
  /* Now go into code for standalone solver
     Could copy arguments and add -quit at end to be safe
     but this will do
  */
  if (argc > 3) {
    CbcMain1(argc - 2, argv + 2, modelA, callBack,solverData);
  } else {
    const char *argv2[] = { "lotsizing", "-solve", "-quit" };
    CbcMain1(3, argv2, modelA, callBack,solverData);
  }
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
  return 0;
}
