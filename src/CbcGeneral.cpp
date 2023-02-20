// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CoinTypes.h"
#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcGeneral.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

// Default Constructor
CbcGeneral::CbcGeneral()
  : CbcObject()
{
}

// Constructor from model
CbcGeneral::CbcGeneral(CbcModel *model)
  : CbcObject(model)
{
}

// Destructor
CbcGeneral::~CbcGeneral()
{
}

// Copy constructor
CbcGeneral::CbcGeneral(const CbcGeneral &rhs)
  : CbcObject(rhs)
{
}
#ifdef CBC_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"
#include "ClpNode.hpp"
#include "CbcBranchDynamic.hpp"
// Assignment operator
CbcGeneral &
CbcGeneral::operator=(const CbcGeneral &rhs)
{
  if (this != &rhs) {
    CbcObject::operator=(rhs);
  }
  return *this;
}
// Infeasibility - large is 0.5
double
CbcGeneral::infeasibility(const OsiBranchingInformation * /*info*/,
  int & /*preferredWay*/) const
{
  abort();
  return 0.0;
}
CbcBranchingObject *
CbcGeneral::createCbcBranch(OsiSolverInterface * /*solver*/, const OsiBranchingInformation * /*info*/, int /*way*/)
{
  abort();
  return NULL;
}
#endif
#ifdef CBC_TRY_SCIP
// I have put Scip "interface" here to encourage modifying it
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
/* type -
   0 - use continuous solver
   1 - use current solver
*/
int tryScip(CbcModel * model, int type)
{
  double time1 = CoinCpuTime();
  OsiSolverInterface * solver2;
  if (type==1)
    solver2 = model->solver();
  else
    solver2 = model->continuousSolver();
  const double *objective = solver2->getObjCoefficients();
  const double *lower = solver2->getColLower();
  const double *upper = solver2->getColUpper();
  const double *rowLower = solver2->getRowLower();
  const double *rowUpper = solver2->getRowUpper();
  int numberRows = solver2->getNumRows();
  int numberColumns = solver2->getNumCols();
  // Row copy
  CoinPackedMatrix matrixByRow(*solver2->getMatrixByRow());
  double * elementByRow = matrixByRow.getMutableElements();
  const int *column = matrixByRow.getIndices();
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();
  
  // Scip is under Apache 2 license
  //  http://www.apache.org/licenses/LICENSE-2.0 
  SCIP* scip = NULL;
  SCIP_CALL(SCIPcreate(&scip));
  SCIP_CALL(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL(SCIPcreateProbBasic(scip,"fromCbc"));
  // TODO - By now objective has been flipped - later redo
  if (model->modelFlipped()) 
    model->messageHandler()->message(CBC_GENERAL, model->messages())
    << "Objective already flipped - will fix" << CoinMessageEol;
  SCIP_CALL(SCIPsetObjsense(scip,SCIP_OBJSENSE_MINIMIZE));
  // Columns
  std::vector<SCIP_VAR *> variable;
  variable.reserve(numberColumns);
  for (int i=0;i<numberColumns;i++) {
    SCIP_VAR* var = NULL;
    SCIP_VARTYPE columnType = SCIP_VARTYPE_CONTINUOUS;
    if (solver2->isInteger(i)) {
      if (!lower[i] && upper[i] == 1.0)
	columnType = SCIP_VARTYPE_BINARY;
      else
	columnType = SCIP_VARTYPE_INTEGER;
    }
    SCIP_CALL(SCIPcreateVarBasic(scip,&var,solver2->getColName(i).c_str(), 
				 lower[i],upper[i],objective[i],
				 columnType));
    SCIP_CALL(SCIPaddVar(scip,var));
    variable[i]=var;
  }
  // Constraints
  // Might be able to add all at once ??? - for speed
  std::vector< SCIP_CONS* > constraint;
  SCIP_VAR** tempvars = new SCIP_VAR * [numberColumns];;
  constraint.reserve(numberRows);
  for (int i=0;i<numberRows;i++) {
    SCIP_CONS* cons = NULL;
    CoinBigIndex startrow = rowStart[i];
    int length = rowLength[i];
    int endrow = startrow+length;
    for (CoinBigIndex j=startrow;j<endrow;j++)
      tempvars[j-startrow] = variable[column[j]];
    SCIP_CALL(SCIPcreateConsBasicLinear(scip,&cons,
					solver2->getRowName(i).c_str(),length,tempvars,
					elementByRow+startrow,rowLower[i],rowUpper[i]));
    SCIP_CALL(SCIPaddCons(scip, cons));
    constraint[i]=cons;
  }
  delete [] tempvars;
  double offset;
  solver2->getDblParam(OsiObjOffset,offset);
  SCIP_CALL(SCIPaddOrigObjoffset(scip,offset));
  // TODO later add code for SOS
#if 0
  // Important TODO s
  double * bestSolution = model->bestSolution();
  if (bestSolution) {
    SCIP_CALL(SCIPincludeHeurTrySol(scip));
    SCIP_HEUR * heuristic;
    SCIP_SOL * solution;
    SCIP_CALL(SCIPheurPassSolTrySol(scip,heuristic,solution));
  }
  if (model->getCutoff() < 1.0e50) {
    // how do I pass in ???
    //scip->primal->cutoffbound = model->getCutoff();
  }
#endif
  int logLevel = 3;
  if (model->logLevel() ==0)
    logLevel = 0;
  else if (model->logLevel() > 3)
    logLevel = 5;
  SCIP_CALL(SCIPsetIntParam(scip, "display/verblevel", logLevel));
  SCIP_CALL(SCIPsolve(scip));
  SCIP_STATUS solutionStatus = SCIPgetStatus(scip);
  long numberNodes = SCIPgetNTotalNodes(scip);
  long numberIterations = SCIPgetNLPIterations(scip);
  int status;
  if (solutionStatus==SCIP_STATUS_OPTIMAL) {
    status = 0;
  } else if (solutionStatus==SCIP_STATUS_USERINTERRUPT) {
    status = 5;
  } else if (solutionStatus==SCIP_STATUS_UNKNOWN) {
    status = 2;
  } else if (solutionStatus==SCIP_STATUS_INFEASIBLE) {
    status = 0;
  } else if (solutionStatus==SCIP_STATUS_UNBOUNDED) {
    status = 0;
  } else {
    status = 1;
  } 
  double timeTaken = CoinCpuTime() - time1;
  char printBuffer[100];
  SCIP_SOL* sol = SCIPgetBestSol(scip);
  if(sol) {
    double objectiveValue = SCIPgetSolOrigObj(scip,sol);
    sprintf(printBuffer,
	    "Solution value of %g found by Scip taking %.2f seconds",
	    objectiveValue, timeTaken);
  } else {
    sprintf(printBuffer,"Scip took %.2f seconds", timeTaken);
  }
  model->messageHandler()->message(CBC_GENERAL, model->messages())
    << printBuffer << CoinMessageEol;
  model->incrementNodeCount(numberNodes);
  model->incrementIterationCount(numberIterations);
  if(sol) {
    double * newSolution = new double [numberColumns];
    double objectiveValue = SCIPgetSolOrigObj(scip,sol);
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = SCIPgetSolVal(scip, sol, variable[iColumn]);
      newSolution[iColumn] = value; 
    }
    model->setBestSolution(CBC_ROUNDING,objectiveValue,newSolution);
    delete [] newSolution;
  }
  for (int i=0;i<numberColumns;i++) {
    SCIP_CALL(SCIPreleaseVar(scip, &variable[i]));
  }
  variable.clear();
    
  for (int i=0;i<numberRows;i++) {
    SCIP_CALL(SCIPreleaseCons(scip, &constraint[i]));
  }
  constraint.clear();
  
  SCIP_CALL(SCIPfree(&scip));
  return status;
}
#endif
