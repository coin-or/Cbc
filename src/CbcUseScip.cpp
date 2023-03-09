// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"
#if CBC_TRY_SCIP
// I have put Scip "interface" here to encourage modifying it
#include "CbcUseScip.hpp"
// Default Constructor
CbcUseScip::CbcUseScip()
  : model_(NULL)
  , originalSolver_(NULL)
  , modSolver_(NULL)
  , scip_(NULL)
{
}

// Constructor from model
CbcUseScip::CbcUseScip(CbcModel *model)
  : model_(model)
  , originalSolver_(NULL)
  , modSolver_(NULL)
  , scip_(NULL)
{
}

// Destructor
CbcUseScip::~CbcUseScip()
{
}

// Copy constructor
CbcUseScip::CbcUseScip(const CbcUseScip &rhs)
{
  model_ = rhs.model_;
  originalSolver_ = rhs.originalSolver_;
  modSolver_ = rhs.modSolver_;
  scip_ = rhs.scip_;
}
// Assignment operator
CbcUseScip &
CbcUseScip::operator=(const CbcUseScip &rhs)
{
  if (this != &rhs) {
    model_ = rhs.model_;
    originalSolver_ = rhs.originalSolver_;
    modSolver_ = rhs.modSolver_;
    scip_ = rhs.scip_;
  }
  return *this;
}
/* type -
   0 - use continuous solver
   1 - use current solver
*/
int CbcUseScip::tryScip(int type)
{
  double time1 = CoinCpuTime();
  OsiSolverInterface * solver2;
  if (type==1)
    solver2 = model_->solver();
  else
    solver2 = model_->continuousSolver();
  int numberRows = solver2->getNumRows();
  int numberColumns = solver2->getNumCols();
  // in case we need to flip
  double *objective =
    CoinCopyOfArray(solver2->getObjCoefficients(),numberColumns);
  const double *lower = solver2->getColLower();
  const double *upper = solver2->getColUpper();
  const double *rowLower = solver2->getRowLower();
  const double *rowUpper = solver2->getRowUpper();
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
  double offset;
  double flipFactor = 1.0;
  solver2->getDblParam(OsiObjOffset,offset);
  offset = -offset;
  if (model_->modelFlipped()) {
    model_->messageHandler()->message(CBC_GENERAL, model_->messages())
      << "Objective flipped - will be flipped back" << CoinMessageEol;
    flipFactor = -1.0;
    SCIP_CALL(SCIPsetObjsense(scip,SCIP_OBJSENSE_MAXIMIZE));
    for (int i=0;i<numberColumns;i++)
      objective[i] = -objective[i];
    offset = -offset;
  } else {
    SCIP_CALL(SCIPsetObjsense(scip,SCIP_OBJSENSE_MINIMIZE));
  }
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
  SCIP_VAR** tempvars = new SCIP_VAR * [numberColumns];
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
  // Set offset
  SCIP_CALL(SCIPaddOrigObjoffset(scip,offset));
  delete [] objective;
  // code for SOS
  OsiClpSolverInterface *clpSolver =
    dynamic_cast<OsiClpSolverInterface *>(model_->solver());
  if (clpSolver && clpSolver->numberSOS()) {
    double * tempvals = new double [numberColumns];
    int numberSets = clpSolver->numberSOS();
    const CoinSet * setInfo = clpSolver->setInfo();
    char name[20];
    for (int i=0;i<numberSets;i++) {
      int n = setInfo[i].numberEntries();
      const int *which = setInfo[i].which();
      const double *weights = setInfo[i].weights();
      sprintf(name,"SOS_%d",i);
      SCIP_CONS* cons = NULL;
      for (int j = 0; j < n; j++) {
	int iVar = which[j];
	tempvars[j] = variable[iVar];
      }
      memcpy(tempvals,weights,n*sizeof(double));
      int type = setInfo[i].setType();
      if (type==1) {
	// SOS 1
	SCIP_CALL(SCIPcreateConsBasicSOS1(scip,&cons,name,
					  n,tempvars,tempvals));
      } else {
	// SOS 2
	SCIP_CALL(SCIPcreateConsBasicSOS2(scip,&cons,name,
					  n,tempvars,tempvals));
      }
      SCIP_CALL(SCIPaddCons(scip, cons));
    }
    delete [] tempvals;
  }
  double * bestSolution = model_->bestSolution();
  if (bestSolution) {
    SCIP_SOL* sol;
    SCIP_Real* vals;
    SCIP_Bool stored;
    SCIP_CALL(SCIPcreateOrigSol(scip,&sol,NULL));
    SCIP_CALL(SCIPallocBufferArray(scip,&vals,numberColumns));
    for (int i=0;i<numberColumns;i++) {
      tempvars[i] = variable[i];
      vals[i] = bestSolution[i];
    }
    SCIP_CALL(SCIPsetSolVals(scip,sol,numberColumns,tempvars, vals) );
    SCIP_CALL(SCIPaddSolFree(scip,&sol,&stored));
    assert(stored);
    SCIPfreeBufferArray(scip,&vals);
  } 
  if (model_->getCutoff() < 1.0e50) {
    SCIP_CALL(SCIPsetObjlimit(scip,model_->getCutoff()*flipFactor));
  }
  delete [] tempvars;
  int logLevel = 3;
  if (model_->logLevel() ==0)
    logLevel = 0;
  else if (model_->logLevel() > 3)
    logLevel = 5;
  if (0) {
    SCIP_CALL(SCIPsetIntParam(scip,"display/verblevel",logLevel));
    SCIP_CALL(SCIPpresolve(scip));
  }
  SCIP_CALL(SCIPsetIntParam(scip,"display/verblevel",logLevel));
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
  model_->messageHandler()->message(CBC_GENERAL, model_->messages())
    << printBuffer << CoinMessageEol;
  model_->incrementNodeCount(numberNodes);
  model_->incrementIterationCount(numberIterations);
  if(sol) {
    double * newSolution = new double [numberColumns];
    double objectiveValue = SCIPgetSolOrigObj(scip,sol);
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = SCIPgetSolVal(scip, sol, variable[iColumn]);
      newSolution[iColumn] = value; 
    }
    model_->setBestSolution(CBC_ROUNDING,objectiveValue,newSolution);
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
#if CBC_TRY_SCIP > 1
/** transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to vars array to get active variables for */
   SCIP_Real**           scalars,            /**< pointer to scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c  */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(scalars != NULL);
   assert(*vars != NULL);
   assert(*scalars != NULL);
   assert(nvars != NULL);
   assert(constant != NULL);

   if( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

      if( requiredsize > *nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, vars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, scalars, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, requiredsize, constant, &requiredsize, TRUE) );
         assert( requiredsize <= *nvars );
      }
   }
   else
   {
      for( v = 0; v < *nvars; ++v )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&(*vars)[v], &(*scalars)[v], constant) );

         /* negated variables with an original counterpart may also be returned by SCIPvarGetOrigvarSum();
          * make sure we get the original variable in that case
          */
         if( SCIPvarGetStatus((*vars)[v]) == SCIP_VARSTATUS_NEGATED )
         {
            (*vars)[v] = SCIPvarGetNegatedVar((*vars)[v]);
            (*scalars)[v] *= -1.0;
            *constant += 1.0;
         }
      }
   }
   return SCIP_OKAY;
}
#define	SCIPsetIsInfinity(set,val)        ( (val) >= (set)->num_infinity) 
#define SCIPsetInfinity(set)               ( (set)->num_infinity )
// create presolved model
int
CbcUseScip::presolveModel(OsiSolverInterface * solver)
{
  double time1 = CoinCpuTime();
  originalSolver_ = solver;
  int numberRows = solver->getNumRows();
  int numberColumns = solver->getNumCols();
  const double *objective = solver->getObjCoefficients();
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  // Row copy
  CoinPackedMatrix matrixByRow(*solver->getMatrixByRow());
  double * elementByRow = matrixByRow.getMutableElements();
  const int *column = matrixByRow.getIndices();
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();
  
  // Scip is under Apache 2 license
  //  http://www.apache.org/licenses/LICENSE-2.0 
  scip_ = NULL;
  SCIP_CALL(SCIPcreate(&scip_));
  SCIP_CALL(SCIPincludeDefaultPlugins(scip_));
  SCIP_CALL(SCIPcreateProbBasic(scip_,"fromCbc"));
  double offset;
  solver->getDblParam(OsiObjOffset,offset);
  if (solver->getObjSense()<0)
    SCIP_CALL(SCIPsetObjsense(scip_,SCIP_OBJSENSE_MAXIMIZE));
   else 
    SCIP_CALL(SCIPsetObjsense(scip_,SCIP_OBJSENSE_MINIMIZE));
  // Columns
  std::vector<SCIP_VAR *> variable;
  variable.reserve(numberColumns);
  for (int i=0;i<numberColumns;i++) {
    SCIP_VAR* var = NULL;
    SCIP_VARTYPE columnType = SCIP_VARTYPE_CONTINUOUS;
    if (solver->isInteger(i)) {
      if (!lower[i] && upper[i] == 1.0)
	columnType = SCIP_VARTYPE_BINARY;
      else
	columnType = SCIP_VARTYPE_INTEGER;
    }
    SCIP_CALL(SCIPcreateVarBasic(scip_,&var,solver->getColName(i).c_str(), 
				 lower[i],upper[i],objective[i],
				 columnType));
    SCIP_CALL(SCIPaddVar(scip_,var));
    variable[i]=var;
  }
  // Constraints
  // Might be able to add all at once ??? - for speed
  std::vector< SCIP_CONS* > constraint;
  SCIP_VAR** tempvars = new SCIP_VAR * [numberColumns];
  double * tempvals = new double [numberColumns];
  constraint.reserve(numberRows);
  for (int i=0;i<numberRows;i++) {
    SCIP_CONS* cons = NULL;
    CoinBigIndex startrow = rowStart[i];
    int length = rowLength[i];
    int endrow = startrow+length;
    for (CoinBigIndex j=startrow;j<endrow;j++)
      tempvars[j-startrow] = variable[column[j]];
    SCIP_CALL(SCIPcreateConsBasicLinear(scip_,&cons,
					solver->getRowName(i).c_str(),length,tempvars,
					elementByRow+startrow,rowLower[i],rowUpper[i]));
    SCIP_CALL(SCIPaddCons(scip_, cons));
    constraint[i]=cons;
  }
  if (offset)
    printf("OOOOO offset %g\n",offset);// Set offset
  // code for SOS - BUT something wrong
  OsiClpSolverInterface *clpSolver =
      dynamic_cast<OsiClpSolverInterface *>(solver);
  if (clpSolver && clpSolver->numberSOS()) {
    int numberSets = clpSolver->numberSOS();
    const CoinSet * setInfo = clpSolver->setInfo();
    char name[20];
    for (int i=0;i<numberSets;i++) {
      int n = setInfo[i].numberEntries();
      const int *which = setInfo[i].which();
      const double *weights = setInfo[i].weights();
      sprintf(name,"SOS_%d",i);
      SCIP_CONS* cons = NULL;
      for (int j = 0; j < n; j++) {
	int iVar = which[j];
	tempvars[j] = variable[iVar];
      }
      memcpy(tempvals,weights,n*sizeof(double));
      int type = setInfo[i].setType();
      if (type==1) {
	// SOS 1
	SCIP_CALL(SCIPcreateConsBasicSOS1(scip_,&cons,name,
					  n,tempvars,tempvals));
      } else {
	// SOS 2
	SCIP_CALL(SCIPcreateConsBasicSOS2(scip_,&cons,name,
					  n,tempvars,tempvals));
      }
      SCIP_CALL(SCIPaddCons(scip_, cons));
    }
  }
  // Presolve
  int logLevel = 3;
  SCIP_CALL(SCIPsetIntParam(scip_,"display/verblevel",logLevel));
  // ? should problem be transformed before - and not after
  //SCIP_CALL(SCIPtransformProb(scip_));
  SCIP_CALL(SCIPpresolve(scip_));
  int nsols = scip_->primal->nexistingsols;
  for (int i=0;i<nsols;i++) {
    if (scip_->primal->existingsols[i]->vals->valssize)
      printf("SOLSOL %d obj %g valssize %d\n",i,scip_->primal->existingsols[i]->obj,
	   scip_->primal->existingsols[i]->vals->valssize);
  }
  int presolveStatus = 0; // need 1 to say ignore
  if (scip_->set->stage==SCIP_STAGE_SOLVED)
    presolveStatus = 1;
  char printBuffer[100];
  sprintf(printBuffer,"Return code from presolve %d %g seconds - %d rows, %d columns",
	  presolveStatus,CoinCpuTime()-time1,scip_->transprob->nconss,
	  scip_->transprob->nvars);
  solver->messageHandler()->message(COIN_GENERAL_INFO, solver->messages())
    << printBuffer << CoinMessageEol;
  if (presolveStatus)
    return presolveStatus;
  modSolver_ = solver->clone();
  CoinBigIndex guessElements = numberColumns+
    (11*modSolver_->getNumElements())/10;
  // get empty model
  int * whichRow = new int [numberRows+numberColumns];
  int * whichColumn = whichRow+numberRows;
  for (int i=0;i<numberRows;i++)
    whichRow[i] = i;
  for (int i=0;i<numberColumns;i++)
    whichColumn[i] = i;
  modSolver_->deleteRows(numberRows,whichRow);
  modSolver_->deleteCols(numberColumns,whichColumn);
  delete [] whichRow;
  double objscale2 = scip_->transprob->objscale;
  SCIPprobAddObjoffset(scip_->transprob, -offset/objscale2);
#ifdef CBC_WRITESTUFF
  SCIPwriteTransProblem(scip_,"/tmp/realscip.mps","mps",false);
#endif
  SCIP_PROB* prob = scip_->transprob; //SCIPgetProbData(scip_);
  bool transformed=true;
  // code from reader_lp.c
  char *valuestr=NULL;
  SCIP_SET * set = scip_->set;
  SCIP_VAR** vars;
  SCIP_VAR** fixedvars;
  SCIP_CONS** conss;
  SCIP_CONS* cons;
  SCIP_Real objscale;
  char* name="none";
  int nfixedvars;
  int nconss;
  int nvars;
  int i;
  
  vars = prob->vars;
  nvars = prob->nvars;
  fixedvars = prob->fixedvars;
  nfixedvars = prob->nfixedvars;
  
  /* case of the transformed problem, we want to write currently valid problem */
  SCIP_CONSHDLR** conshdlrs;
  int nconshdlrs;
  
  conshdlrs = set->conshdlrs;
  nconshdlrs = set->nconshdlrs;
  
  /* collect number of constraints which have to be enforced; these are the constraints which currency (locally)
   * enabled; these also includes the local constraints
   */
  nconss = 0;
  for( i = 0; i < nconshdlrs; ++i )
    {
      nconss += SCIPconshdlrGetNEnfoConss(conshdlrs[i]);
    }
  
  SCIP_CALL( SCIPsetAllocBufferArray(set, &conss, nconss) );
  
  /* copy the constraints */
  nconss = 0;
  for( i = 0; i < nconshdlrs; ++i )
    {
      SCIP_CONS** conshdlrconss;
      int nconshdlrconss;
      int c;
      
      {
	conshdlrconss = SCIPconshdlrGetEnfoConss(conshdlrs[i]);
	nconshdlrconss = SCIPconshdlrGetNEnfoConss(conshdlrs[i]);
      }
      for( c = 0; c < nconshdlrconss; ++c )
	{
	  conss[nconss] = conshdlrconss[c];
	  nconss++;
	}
    }
  // Osi arrays
  double * columnLower2 = new double[nvars];
  double * columnUpper2 = new double[nvars];
  double * objective2 = new double[nvars];
  int * columnType2 = new int[nvars];
  double * rowLower2 = new double[nconss];
  double * rowUpper2 = new double[nconss];
  CoinBigIndex * rowStart2 = new CoinBigIndex [nconss+1];
  int * column2 = new int [guessElements];
  double * element2 = new double[guessElements];
  
  /* adapt objective scale for transformed problem (for the original no change is necessary) */
  objscale = prob->objscale;
  if( prob->objsense == SCIP_OBJSENSE_MAXIMIZE )
    objscale *= -1.0;
  
  int c;
  int v;
  int k;
  char* namestr;
  
  //char* consname;
  const char** consnames;
  
  SCIP_CONSHDLR* conshdlr;
  const char* conshdlrname;
  
  SCIP_Real lhs;
  SCIP_Real rhs;
  SCIP_Real* rhss;
  SCIP_Real value;
  
  SCIP_VAR* var = NULL;
  
  SCIP_CONS** consIndicator;
  SCIP_CONS** consSOS1;
  SCIP_CONS** consSOS2;
  SCIP_CONS** consQuadratic;
  int nConsIndicator;
  int nConsSOS1;
  int nConsSOS2;
  int nConsQuadratic;
  SCIP_VAR** aggvars;
  int naggvars = 0;
  int saggvars;
  SCIP_HASHTABLE* varFixedHash;
  SCIP_HASHTABLE* indicatorSlackHash;
  
  SCIP_VAR** fixvars = NULL;
  int nfixvars = 0;
  
  SCIP_VAR** consvars;
  int nconsvars;
  SCIP_Real* vals;
  SCIP_Longint* weights;
  
  SCIP_Bool needRANGES;
  unsigned int maxnamelen;
  
  SCIP_Bool error=true;
  
  assert(scip_ != NULL);
  
  needRANGES = FALSE;
  maxnamelen = 0;
  nConsSOS1 = 0;
  nConsSOS2 = 0;
  nConsQuadratic = 0;
  nConsIndicator = 0;
  int nbinvars=prob->nbinvars;
  int nintvars=prob->nintvars;
  int nimplvars=prob->nimplvars;
  int ncontvars=prob->ncontvars;
  double objoffset =prob->objoffset;
  printf("objoffset %g old %g from scip %g\n",objoffset,offset,prob->objoffset);
  SCIP_Bool linearizeands;
  SCIP_Bool aggrlinearizationands;
  
  SCIP_CONSHDLR* conshdlrInd;
  SCIP_CONS** consExpr;
  
  SCIP_VAR** tmpvars;
  int tmpvarssize;
  SCIP_HASHTABLE* varAggregated;
  SCIP_HASHMAP* consHidden;
  
  SCIP_Real* consvals;
  SCIP_Real lb;
  SCIP_Real ub;
  
  
  assert(scip_ != NULL);
  
  for( v = 0; v < nvars; ++v )
    {
      var = vars[v];
      objective2[v] = SCIPvarGetObj(var)*objscale; // ? min/max objscale
      columnLower2[v] = SCIPvarGetLbLocal(var);
      columnUpper2[v] = SCIPvarGetUbLocal(var);
    }
  
  /* collect SOS, quadratic, and SOC constraints in array for later output */
  SCIP_CALL( SCIPallocBufferArray(scip_, &consSOS1, nconss) );
  SCIP_CALL( SCIPallocBufferArray(scip_, &consSOS2, nconss) );
  SCIP_CALL( SCIPallocBufferArray(scip_, &consExpr, nconss) );
  SCIP_CALL( SCIPallocBufferArray(scip_, &consIndicator, nconss) );
  
  tmpvarssize = SCIPgetNTotalVars(scip_);
  SCIP_CALL( SCIPallocBufferArray(scip_, &tmpvars, tmpvarssize) );

  rowStart2[0] = 0;
  CoinBigIndex nels = 0;
  CoinBigIndex increaseAfter = guessElements-numberColumns;
  for( c = 0; c < nconss; ++c )
    {
      cons = conss[c];
      assert( cons != NULL);
      assert(SCIPconsIsEnabled(cons));
      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );
      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );
      if( strcmp(conshdlrname, "linear") == 0 )
	{
	  consvars = SCIPgetVarsLinear(scip_, cons);
	  consvals = SCIPgetValsLinear(scip_, cons);
	  nconsvars = SCIPgetNVarsLinear(scip_, cons);
	  memcpy(tempvars,consvars,nconsvars*sizeof(SCIP_VAR *));
	  memcpy(tempvals,consvals,nconsvars*sizeof(double));
	  lb = SCIPgetLhsLinear(scip_, cons);
	  ub = SCIPgetRhsLinear(scip_, cons);
	  if (lb<-1.999999e19)
	    lb = -COIN_DBL_MAX;
	  if(ub>1.999999e19)
	    ub = COIN_DBL_MAX;
	}
      else if( strcmp(conshdlrname, "setppc") == 0 )
	{
	  consvars = SCIPgetVarsSetppc(scip_, cons);
	  nconsvars = SCIPgetNVarsSetppc(scip_, cons);
	  for (int i=0;i<nconsvars;i++) {
	    double scalar = 1.0;
	    double activeconstant = 0.0;
	    SCIP_VAR *tmpVar = consvars[i];
	    SCIP_CALL(SCIPvarGetProbvarSum(&tmpVar,scip_->set,
	    			   &scalar,&activeconstant));
	    tempvars[i] =consvars[i];// ??tmpVar;
	    tempvals[i] = 1.0;
	  }
	  switch( SCIPgetTypeSetppc(scip_, cons) )
	    {
	    case SCIP_SETPPCTYPE_PARTITIONING :
	      lb = 1.0;
	      ub = 1.0;
	      break;
	    case SCIP_SETPPCTYPE_PACKING :
	      lb = -COIN_DBL_MAX;
	      ub = 1.0;
	      break;
	    case SCIP_SETPPCTYPE_COVERING :
	      lb = 1.0;
	      ub = COIN_DBL_MAX;
	      break;
	    }
	}
      else if( strcmp(conshdlrname, "logicor") == 0 )
	{
	  consvars = SCIPgetVarsLogicor(scip_, cons);
	  nconsvars = SCIPgetNVarsLogicor(scip_, cons);
	  memcpy(tempvars,consvars,nconsvars*sizeof(SCIP_VAR *));
	  
	  for( v = 0; v < nconsvars; ++v )
            tempvals[v] = 1.0;
	  lb = 1.0;
	  ub = COIN_DBL_MAX;
	}
      else if( strcmp(conshdlrname, "knapsack") == 0 )
	{
	  SCIP_Longint* weights;
	  
	  consvars = SCIPgetVarsKnapsack(scip_, cons);
	  nconsvars = SCIPgetNVarsKnapsack(scip_, cons);
	  memcpy(tempvars,consvars,nconsvars*sizeof(SCIP_VAR *));
	  
	  weights = SCIPgetWeightsKnapsack(scip_, cons);
	  for( v = 0; v < nconsvars; ++v )
            tempvals[v] = (SCIP_Real)weights[v];
	  lb = -COIN_DBL_MAX;
	  ub = SCIPgetCapacityKnapsack(scip_, cons);
	}
      else if( strcmp(conshdlrname, "varbound") == 0 )
	{
	  tempvars[0] = SCIPgetVarVarbound(scip_, cons);
	  tempvars[1] = SCIPgetVbdvarVarbound(scip_, cons);
	  
	  tempvals[0] = 1.0;
	  tempvals[1] = SCIPgetVbdcoefVarbound(scip_, cons);
	  nconsvars = 2;
	  lb = SCIPgetLhsVarbound(scip_, cons);
	  ub = SCIPgetRhsVarbound(scip_, cons);
	  if (lb<-1.999999e19)
	    lb = -COIN_DBL_MAX;
	  if(ub>1.999999e19)
	    ub = COIN_DBL_MAX;
	}
      else
	{
	  abort(); // no handler
	}
      {
	double activeconstant = 0.0;
	/* retransform given variables to active variables */
      	SCIP_CALL( getActiveVariables(scip_, &tempvars, &tempvals, &nconsvars, &activeconstant, true) );
	if (lb>-1.999999e19)
	   lb -= activeconstant;
	if (ub<1.999999e19)
	   ub -= activeconstant;
	rowLower2[c] = lb;
	rowUpper2[c] = ub;
	for (int i=0;i<nconsvars;i++) {
	  int iColumn = tempvars[i]->probindex;
	  assert (iColumn>=0&&iColumn<numberColumns);
	  double value = tempvals[i];
	  if (value) {
	    column2[nels] = iColumn;
	    element2[nels++] = value;
	  }
	}
	rowStart2[c+1] = nels;
      }
      if (nels > increaseAfter) {
	guessElements += numberColumns+guessElements/5;
	increaseAfter = guessElements-numberColumns;
	int * column2a = new int [guessElements];
	memcpy(column2a,column2,nels*sizeof(int));
	delete [] column2;
	column2 = column2a;
	double * element2a = new double[guessElements];
	memcpy(element2a,element2,nels*sizeof(double));
	delete [] element2; 
	element2 = element2a;
      }
    }
  {
    printf("Cgl0004 scip presolve %d rows and  %d columns\n",nconss,nvars);
    CoinBigIndex * dummyStart = new CoinBigIndex [nvars+1];
    memset(dummyStart,0,(nvars+1)*sizeof(CoinBigIndex));
    modSolver_->addCols(nvars,dummyStart,NULL,NULL,
			   columnLower2,columnUpper2,objective2);
    delete [] dummyStart;
    for (int i=0;i<nvars;i++) {
      SCIP_VARTYPE type = SCIPvarGetType(vars[i]);
      if (type == SCIP_VARTYPE_BINARY || type == SCIP_VARTYPE_INTEGER) {
	modSolver_->setInteger(i);
      } else {
	assert (type == SCIP_VARTYPE_IMPLINT
		|| type == SCIP_VARTYPE_CONTINUOUS);
      }
    }
    modSolver_->addRows(nconss,rowStart2,column2,element2,
			rowLower2,rowUpper2);
    modSolver_->setDblParam(OsiObjOffset,-objoffset*objscale);
#ifdef CBC_WRITESTUFF
    modSolver_->writeMps("/tmp/modcbc");
#endif
  }
  delete [] tempvars;
  delete [] tempvals;
  delete [] columnLower2;
  delete [] columnUpper2;
  delete [] objective2;
  delete [] columnType2;
  delete [] rowLower2;
  delete [] rowUpper2;
  delete [] rowStart2;
  delete [] column2;
  delete [] element2;
  
  /* free space */
  //SCIPfreeBlockMemoryArray(scip_, &aggvars, saggvars);
  //SCIPhashtableFree(&varAggregated);
  //if( conshdlrInd != NULL )
  // SCIPhashmapFree(&consHidden);
  
  
  /* free space */
  SCIPfreeBufferArray(scip_, &tmpvars);
  SCIPfreeBufferArray(scip_, &consIndicator);
  SCIPfreeBufferArray(scip_, &consExpr);
  SCIPfreeBufferArray(scip_, &consSOS2);
  SCIPfreeBufferArray(scip_, &consSOS1);
  return presolveStatus;
}
// postprocess model
// saveOriginalSolver has been saved
int CbcUseScip::afterSolve(const double * solin, double objValue)
{
  // only call if solution
  if (!solin) {
    // try and release memory
    return 1;
  }
  int numberColumns = modSolver_->getNumCols();
  int maxColumns = 0;
  SCIP_VAR** vars = scip_->transprob->vars;
  SCIP_VAR** vars2 = scip_->origprob->vars;
  int nvars2 = scip_->origprob->nvars;
  for (int i=0;i<numberColumns;i++) {
    maxColumns = CoinMax(maxColumns,vars[i]->index);
  }
  maxColumns++;
  SCIP_SOL *sol;
  BMS_BLKMEM *blkmem = scip_->mem->setmem;
  SCIP_HEUR *dummyHeuristic = scip_->set->heurs[58];
  // change if things change
  assert (!strcmp(SCIPheurGetName(dummyHeuristic),"trysol"));
  SCIP_CALL(SCIPsolCreate(&sol,blkmem,scip_->set,scip_->stat,
			  scip_->primal,NULL,dummyHeuristic));
  sol->obj = objValue;
  SCIP_CALL(SCIPrealarrayExtend(sol->vals,0,1.2,0,maxColumns));
  SCIP_CALL(SCIPboolarrayExtend(sol->valid,0,1.2,0,maxColumns));
  sol->valid->minusedidx = 0;
  sol->valid->maxusedidx = maxColumns;
  sol->vals->minusedidx = 0;
  sol->vals->maxusedidx = maxColumns;
  for (int i=0;i<maxColumns;i++) {
    sol->valid->vals[i] = false;
    sol->vals->vals[i] = 0.0;
  }
  for (int i=0;i<numberColumns;i++) {
    int index = vars[i]->index;
    sol->valid->vals[index] = true;
    sol->vals->vals[index] = solin[i];
  }
  SCIP_Bool hasinfval;
  
  /* retransform solution into the original problem space */
  SCIP_CALL(SCIPsolRetransform(sol,scip_->set,
			       scip_->stat, scip_->origprob,
			       scip_->transprob, &hasinfval) );
  int n = originalSolver_->getNumCols();
  for (int i=0;i<n;i++) {
    double value = sol->vals->vals[i];
    if (originalSolver_->isInteger(i)) {
      originalSolver_->setColLower(i,value);
      originalSolver_->setColUpper(i,value);
    }
  }
  originalSolver_->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
  originalSolver_->resolve();
  printf("Objective value now %g status %s\n", 
	 originalSolver_->getObjValue(),
	 originalSolver_->isProvenOptimal() ? "good" : "bad");
  // When I try and free scip I get lots of errors - can someone help
  // try and release memory
  return 0;
}
#endif
