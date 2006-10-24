// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <cassert>

#include "CoinTime.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinMpsIO.hpp"
#include "CoinModel.hpp"
#include "ClpSimplex.hpp"
#include "CbcSolverLink.hpp"
#include "CbcModel.hpp"
#include "ClpPackedMatrix.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"
#include "CoinTime.hpp"
#include "CglProbing.hpp"
#include "CglProbing.hpp"
#include "CglGomory.hpp"
#include "CglKnapsackCover.hpp"
#include "CglRedSplit.hpp"
#include "CglClique.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglFlowCover.hpp"
#include "CglTwomir.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcCompareActual.hpp"
#include "CbcBranchLink.hpp"
//#############################################################################
// Solve methods
//#############################################################################
void CbcSolverLink::initialSolve()
{
  specialOptions_ =0;
  modelPtr_->setWhatsChanged(0);
  ClpMatrixBase * save = NULL;
  if (numberVariables_) {
    CoinPackedMatrix * temp = new CoinPackedMatrix(*matrix_);
    // update all bounds before coefficients
    for (int i=0;i<numberVariables_;i++ ) {
      info_[i].updateBounds(modelPtr_);
    }
    for (int i=0;i<numberVariables_;i++ ) {
      info_[i].updateCoefficients(modelPtr_,temp);
    }
    temp->removeGaps(1.0e-14);
    save = modelPtr_->clpMatrix();
    ClpPackedMatrix * clpMatrix = dynamic_cast<ClpPackedMatrix *> (save);
    assert (clpMatrix);
    if (save->getNumRows()>temp->getNumRows()) {
      // add in cuts
      int numberRows = temp->getNumRows();
      int * which = new int[numberRows];
      for (int i=0;i<numberRows;i++)
	which[i]=i;
      save->deleteRows(numberRows,which);
      delete [] which;
      temp->bottomAppendPackedMatrix(*clpMatrix->matrix());
    }
    modelPtr_->replaceMatrix(temp);
  }
  OsiClpSolverInterface::initialSolve();
  if (save) {
    delete save;
  }
}

//-----------------------------------------------------------------------------
void CbcSolverLink::resolve()
{
  specialOptions_ =0;
  modelPtr_->setWhatsChanged(0);
  ClpMatrixBase * save = NULL;
  bool allFixed=false;
  bool feasible=true;
  if (numberVariables_) {
    CoinPackedMatrix * temp = new CoinPackedMatrix(*matrix_);
    allFixed=true;
    //bool best=true;
    const double * lower = modelPtr_->columnLower();
    const double * upper = modelPtr_->columnUpper();
    // update all bounds before coefficients
    for (int i=0;i<numberVariables_;i++ ) {
      info_[i].updateBounds(modelPtr_);
      int iColumn = info_[i].variable();
      double lo = lower[iColumn];
      double up = upper[iColumn];
      if (up>lo)
	allFixed=false;
      else if (up<lo)
	feasible=false;
    }
    for (int i=0;i<numberVariables_;i++ ) {
      info_[i].updateCoefficients(modelPtr_,temp);
    }
    temp->removeGaps(1.0e-14);
    save = modelPtr_->clpMatrix();
    ClpPackedMatrix * clpMatrix = dynamic_cast<ClpPackedMatrix *> (save);
    assert (clpMatrix);
    if (save->getNumRows()>temp->getNumRows()) {
      // add in cuts
      int numberRows = temp->getNumRows();
      int * which = new int[numberRows];
      for (int i=0;i<numberRows;i++)
	which[i]=i;
      CoinPackedMatrix * mat = clpMatrix->matrix();
      // for debug
      //mat = new CoinPackedMatrix(*mat);
      mat->deleteRows(numberRows,which);
      delete [] which;
      temp->bottomAppendPackedMatrix(*mat);
      temp->removeGaps(1.0e-14);
    }
    if (0) {
      const CoinPackedMatrix * matrix = modelPtr_->matrix();
      int numberColumns = matrix->getNumCols();
      assert (numberColumns==temp->getNumCols());
      const double * element1 = temp->getMutableElements();
      const int * row1 = temp->getIndices();
      const CoinBigIndex * columnStart1 = temp->getVectorStarts();
      const int * columnLength1 = temp->getVectorLengths();
      const double * element2 = matrix->getMutableElements();
      const int * row2 = matrix->getIndices();
      const CoinBigIndex * columnStart2 = matrix->getVectorStarts();
      const int * columnLength2 = matrix->getVectorLengths();
      for (int i=0;i<numberColumns;i++) {
	assert (columnLength2[i]==columnLength1[i]);
	int offset = columnStart2[i]-columnStart1[i];
	for (int j=columnStart1[i];j<columnStart1[i]+columnLength1[i];j++) {
	  assert (row1[j]==row2[j+offset]);
	  assert (element1[j]==element2[j+offset]);
	}
      }
    }
    modelPtr_->replaceMatrix(temp);
  }
  if (!feasible)
    allFixed=false;
  if ((specialOptions2_&1)!=0)
    allFixed=false;
  int returnCode=-1;
  if (feasible) {
    // may do lots of work
    returnCode=fathom(allFixed);
  }
  if (returnCode>=0) {
    if (returnCode==0)
      OsiClpSolverInterface::resolve();
    if (!allFixed&&(specialOptions2_&1)==0) {
      const double * solution = getColSolution();
      bool satisfied=true;
      for (int i=0;i<numberVariables_;i++) {
	int iColumn = info_[i].variable();
	double value = solution[iColumn];
	if (fabs(value-floor(value+0.5))>0.0001)
	  satisfied=false;
      }
      //if (satisfied)
      //printf("satisfied but not fixed\n");
    }
  } else {
    modelPtr_->setProblemStatus(1);
    modelPtr_->setObjectiveValue(COIN_DBL_MAX);
  }
  if (save) {
    delete save;
  }
}

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CbcSolverLink::CbcSolverLink ()
  : OsiClpSolverInterface()
{
  gutsOfDestructor(true);
}
/* returns
   sequence of nonlinear or
   -1 numeric
   -2 not found
   -3 too many terms
*/
static int getVariable(const CoinModel & model, char * expression,
		       int & linear)
{
  int non=-1;
  linear=-1;
  if (strcmp(expression,"Numeric")) {
    // function
    char * first = strchr(expression,'*');
    int numberColumns = model.numberColumns();
    int j;
    if (first) {
      *first='\0';
      for (j=0;j<numberColumns;j++) {
	if (!strcmp(expression,model.columnName(j))) {
	  linear=j;
	  memmove(expression,first+1,strlen(first+1)+1);
	  break;
	}
      }
    }
    // find nonlinear
    for (j=0;j<numberColumns;j++) {
      const char * name = model.columnName(j);
      first = strstr(expression,name);
      if (first) {
	if (first!=expression&&isalnum(*(first-1)))
	  continue; // not real match
	first += strlen(name);
	if (!isalnum(*first)) {
	  // match
	  non=j;
	  // but check no others
	  j++;
	  for (;j<numberColumns;j++) {
	    const char * name = model.columnName(j);
	    first = strstr(expression,name);
	    if (first) {
	      if (isalnum(*(first-1)))
		continue; // not real match
	      first += strlen(name);
	      if (!isalnum(*first)) {
		// match - ouch
		non=-3;
		break;
	      }
	    }
	  }
	  break;
	}
      }
    }
    if (non==-1)
      non=-2;
  } 
  return non;
}
/* This creates from a coinModel object 

   if errors.then number of sets is -1
      
   This creates linked ordered sets information.  It assumes -

   for product terms syntax is yy*f(zz)
   also just f(zz) is allowed 
   and even a constant

   modelObject not const as may be changed as part of process.
*/
CbcSolverLink::CbcSolverLink ( CoinModel & coinModel)
  : OsiClpSolverInterface()
{
  gutsOfDestructor(true);
  gdb(coinModel);
}
void CbcSolverLink::gdb ( CoinModel & coinModel)
{
  // first check and set up arrays
  int numberColumns = coinModel.numberColumns();
  int numberRows = coinModel.numberRows();
  // List of nonlinear entries
  int * which = new int[numberColumns];
  numberSets_=0;
  numberNonlinear_=0;
  int iColumn;
  int numberErrors=0;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinModelLink triple=coinModel.firstInColumn(iColumn);
    bool linear=true;
    int n=0;
    while (triple.row()>=0) {
      int iRow = triple.row();
      const char * expr = coinModel.getElementAsString(iRow,iColumn);
      if (strcmp(expr,"Numeric")) {
	linear=false;
	// associate to avoid errors - much better to use current solution values
	coinModel.associateElement(expr,0.0);
      }
      triple=coinModel.next(triple);
      n++;
    }
    if (!linear) {
      which[numberSets_++]=iColumn;
      numberNonlinear_ += n;
    }
  }
  // return if nothing
  simplexModel_ = new ClpSimplex(*modelPtr_);
  simplexModel_->loadProblem(coinModel);
  if (!numberSets_) {
    delete [] which;
    return;
  } else {
    simplexModel_->deleteColumns(numberSets_,which);
    coinModel_ = coinModel;
    which_ = CoinCopyOfArray(which,numberSets_);
    delete [] which;
    if (!modelPtr_->numberColumns()) {
      delete modelPtr_;
      modelPtr_ = new ClpSimplex(*simplexModel_);
      // synchronize integers
      int numberColumns=simplexModel_->numberColumns();
      for (int i=0;i<numberColumns;i++) {
	if (simplexModel_->isInteger(i))
	  this->setInteger(i);
      }
    }
  }
  // set starts
  startSet_ = new int[numberSets_+1];
  definitionRow_ = new int [numberNonlinear_];
  definitionRowL_ = new int [numberNonlinear_];
  functionColumn_ = new int [numberSets_];
  definitionColumnL_ = new int [numberNonlinear_];
  usageRow_ = new int [numberNonlinear_];
  convexityRow_ = new int[numberNonlinear_];

  startSet_[0]=0;
  int nonLinear=0;
  int nRow=numberRows;
  int iSet;
  char * mentioned=new char[numberColumns];
  memset(mentioned,0,numberColumns);
  for (iSet=0;iSet<numberSets_;iSet++) {
    iColumn = which_[iSet];
    // first get constants
    CoinModelLink triple=coinModel.firstInColumn(iColumn);
    int functionColumn=-1;
    while (triple.row()>=0) {
      int iRow = triple.row();
      const char * expr = coinModel.getElementAsString(iRow,iColumn);
      int linear;
      assert (strlen(expr)<1000);
      char expr2 [1000];
      strcpy(expr2,expr);
      int non = getVariable(coinModel,expr2,linear);
      if (non==-1) {
	// add constant
	definitionRowL_[nonLinear]=-2;
	definitionRow_[nonLinear]=-2;
	definitionColumnL_[nonLinear]=-2;
	usageRow_[nonLinear]=iRow;
      } else if (non<0) {
	// error
	printf("bad expression %s for column %s and row %s\n",
	       expr,coinModel.getColumnName(iColumn),coinModel.getRowName(iRow));
	numberErrors++;
      } else {
	if (functionColumn==-1) {
	  functionColumn=non;
	} else if (functionColumn!=non) {
	  printf("second function of %s when already got %s - expression %s\n",
		 coinModel.getColumnName(non),coinModel.getColumnName(functionColumn),
		 expr);
	  numberErrors++;
	}
      }
      triple=coinModel.next(triple);
    }
    if (!numberErrors) {
      triple=coinModel.firstInColumn(iColumn);
      functionColumn_[iSet]=functionColumn;
      while (triple.row()>=0) {
	int iRow = triple.row();
	const char * expr = coinModel.getElementAsString(iRow,iColumn);
	int linear;
	assert (strlen(expr)<1000);
	char expr2 [1000];
	strcpy(expr2,expr);
	int non = getVariable(coinModel,expr2,linear);
	if (non>=0) {
	  mentioned[non]=1;
	  if(linear>=0)
	    mentioned[linear]=1;
	  // definition row
	  double elementValue;
	  definitionRow_[nonLinear]=nRow;
	  elementValue = -1.0;
	  simplexModel_->addRow(1,&functionColumn,&elementValue,0.0,0.0);
	  nRow++;
	  // convexity row
	  convexityRow_[nonLinear]=nRow;
 	  simplexModel_->addRow(0,NULL,NULL,1.0,1.0);
	  nRow++;
	  // other definition row
	  if (linear>=0) {
	    definitionRowL_[nonLinear]=nRow;
	    definitionColumnL_[nonLinear]=linear;
	    elementValue = -1.0;
 	    simplexModel_->addRow(1,&linear,&elementValue,0.0,0.0);
	    nRow++;
	  } else {
	    definitionRowL_[nonLinear]=-1;
	    definitionColumnL_[nonLinear]=-1;
	  }
	  usageRow_[nonLinear]=iRow;
	  coinModel_.setElement(iRow,functionColumn,expr2);
	  nonLinear++;
	}
	triple=coinModel.next(triple);
      }
    }
    startSet_[iSet+1]=nonLinear;
  }
  if (numberErrors) {
    // errors
    gutsOfDestructor();
    numberSets_=-1;
  } else {
    int numberRows=coinModel_.numberRows();
    for (int i=0;i<numberColumns;i++) {
      if (mentioned[i])
	coinModel_.setElement(numberRows,i,coinModel_.getColumnName(i));
    }
  }
  delete [] mentioned;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface * 
CbcSolverLink::clone(bool copyData) const
{
  assert (copyData);
  return new CbcSolverLink(*this);
}


//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CbcSolverLink::CbcSolverLink (
                  const CbcSolverLink & rhs)
  : OsiClpSolverInterface(rhs)
{
  gutsOfDestructor(true);
  gutsOfCopy(rhs);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CbcSolverLink::~CbcSolverLink ()
{
  gutsOfDestructor();
}

//-------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CbcSolverLink &
CbcSolverLink::operator=(const CbcSolverLink& rhs)
{
  if (this != &rhs) { 
    gutsOfDestructor();
    OsiClpSolverInterface::operator=(rhs);
    gutsOfCopy(rhs);
  }
  return *this;
}
void 
CbcSolverLink::gutsOfDestructor(bool justNullify)
{
  if (!justNullify) {
    delete matrix_;
    delete [] info_;
    delete [] bestSolution_;
    delete [] startSet_;
    delete [] definitionRow_;
    delete [] definitionRowL_;
    delete [] functionColumn_;
    delete [] definitionColumnL_;
    delete [] usageRow_;
    delete [] convexityRow_;
    delete [] granularityTop_;
    delete [] granularityBottom_;
    delete [] which_;
    
    delete simplexModel_;
    assert (savedModel_==NULL);
  } 
  matrix_ = NULL;
  info_ = NULL;
  bestSolution_ = NULL;
  startSet_ = NULL;
  definitionRow_ = NULL;
  definitionRowL_ = NULL;
  definitionColumnL_ = NULL;
  functionColumn_ = NULL; 
  usageRow_ = NULL;
  convexityRow_ = NULL;
  simplexModel_ = NULL;
  savedModel_=NULL;
  granularityTop_ = NULL;
  granularityBottom_ = NULL; 
  which_ = NULL;

  model_=NULL;
  numberVariables_ = 0;
  bestObjectiveValue_ =1.0e100;
  specialOptions2_ = 0;
  numberSets_ = 0;
  numberNonlinear_ = 0;
}
void 
CbcSolverLink::gutsOfCopy(const CbcSolverLink & rhs)
{
  model_ = rhs.model_;
  coinModel_ = rhs.coinModel_;
  numberVariables_ = rhs.numberVariables_;
  bestObjectiveValue_ = rhs.bestObjectiveValue_;
  specialOptions2_ = rhs.specialOptions2_;
  if (numberVariables_) { 
    matrix_ = new CoinPackedMatrix(*rhs.matrix_);
    info_ = new CbcLinkedBound [numberVariables_];
    for (int i=0;i<numberVariables_;i++) {
      info_[i] = CbcLinkedBound(rhs.info_[i]);
    }
    if (rhs.bestSolution_) {
      bestSolution_ = CoinCopyOfArray(rhs.bestSolution_,modelPtr_->getNumCols());
    } else {
      bestSolution_=NULL;
    }
  }
  numberSets_ = rhs.numberSets_;
  numberNonlinear_ = rhs.numberNonlinear_;
  if (numberSets_>0) {
    startSet_ = CoinCopyOfArray(rhs.startSet_,numberSets_+1);
    definitionRow_ = CoinCopyOfArray(rhs.definitionRow_,numberNonlinear_);
    definitionRowL_ = CoinCopyOfArray(rhs.definitionRowL_,numberNonlinear_);
    functionColumn_ = CoinCopyOfArray(rhs.functionColumn_,numberSets_);
    definitionColumnL_ = CoinCopyOfArray(rhs.definitionColumnL_,numberNonlinear_);
    convexityRow_ = CoinCopyOfArray(rhs.convexityRow_,numberNonlinear_);
    usageRow_ = CoinCopyOfArray(rhs.usageRow_,numberNonlinear_);
    simplexModel_ = new ClpSimplex(*rhs.simplexModel_);
    assert (rhs.savedModel_==NULL);
    granularityTop_ = CoinCopyOfArray(rhs.granularityTop_,numberSets_);
    granularityBottom_ = CoinCopyOfArray(rhs.granularityBottom_,numberSets_);
    which_ = CoinCopyOfArray(rhs.which_,numberSets_);
  }
}
//-------------------------------------------------------------------
// Real initializer
//-------------------------------------------------------------------
void
CbcSolverLink::initialize(CbcModel * model, const CoinPackedMatrix * matrix,
			  int numberVariables, const int * which,
			  const int *starts, const int * positionL, const int * positionU)
{
  model_=model;
  numberVariables_ = numberVariables;
  matrix_ = new CoinPackedMatrix(*matrix);
  info_ = new CbcLinkedBound [numberVariables_];
  // get multipliers from matrix
  int n = starts[numberVariables_];
  double * multipliers = new double[n];
  int i;
  const double * elements = matrix_->getElements();
  for ( i=0;i<n;i++) {
    int position = positionL[i];
    double value = elements[position];
    multipliers[i]=value;
    position = positionU[i];
    value = elements[position];
    assert(multipliers[i]==value);
  }
  for ( i=0;i<numberVariables_;i++) {
    info_[i] = CbcLinkedBound(model,which[i],starts[i+1]-starts[i],
			      positionL+starts[i],
			      positionU+starts[i],
			      multipliers+starts[i]);
  }
  delete [] multipliers;
}
void
CbcSolverLink::initialize ( CbcModel * model, bool top)
{
  model_=model;
  int numberExtraColumns=0;
  int numberExtraElements=0;
  int iSet;
  int iColumn;
  double * granularity = CoinCopyOfArray(top ? granularityTop_ : granularityBottom_, numberSets_);
  if (!granularity) {
    printf("no granularity given\n");
    abort();
  }
  int * numberInSet = new int[numberSets_];
  int * setType = new int[numberSets_];
  assert (simplexModel_);
  // Use Current bounds
  const double * lower = modelPtr_->columnLower();
  const double * upper = modelPtr_->columnUpper();
  // To point to interesting variables
  int numberColumns = simplexModel_->numberColumns();
  int * back = new int[numberColumns];
  CoinZeroN(back,numberColumns);
  int largestSet=0;
  for (iSet=0;iSet<numberSets_;iSet++) {
    int zColumn = functionColumn_[iSet];
    double lo = lower[zColumn];
    double up = upper[zColumn];
    assert (up-lo<1.0e10);
    double g = granularity[iSet];
    assert (g);
    int number;
    if (g>0) {
      number = (int) floor((up-lo)/g+0.5);
      if (number<3)
	number = (int) ceil((up-lo)/g);
      setType[iSet]=2;
    } else {
      g = -g;
      // must be exact
      number = (int) floor((up-lo)/g+1.0e-6*g);
      if (fabs(up-lo+number*g)>1.0e-6) {
	printf("granularity of %g for set %d does not match bounds %g and %g\n",
	       -g,iSet,lo,up);
	abort();
      }
      setType[iSet]=1;
    }
    granularity[iSet] = (up-lo) / ((double) number);
    // and end point
    number++;
    numberInSet[iSet]=number;
    int nElements=0;
    int nLink=0;
    int i;
    // allow zeros to stay
    for (i=startSet_[iSet];i<startSet_[iSet+1];i++) {
      if (definitionRow_[i]>=0) {
	nElements += 2*number;
	nLink++;
      }
      if (definitionRowL_[i]>=0)
	nElements += 2*number;
      // output
      nElements += 2*number;
      // convexity
      nElements += 2*number;
      int iColumn = definitionColumnL_[i];
      if (iColumn>=0)
	back[iColumn]=1;
    }
    int functionColumn = functionColumn_[iSet];
    back[functionColumn]=1;
    number *= 2*nLink;
    numberExtraColumns+=number;
    largestSet = CoinMax(largestSet,number);
    numberExtraElements+=number*nElements;
  }
  CbcObject ** objects = new CbcObject * [numberSets_];
  double * newLower = new double[numberExtraColumns];
  double * newUpper = new double[numberExtraColumns];
  double * newObjective = new double[numberExtraColumns];
  for (iColumn=0;iColumn<numberExtraColumns;iColumn++) {
    newLower[iColumn]=0.0;
    newUpper[iColumn]=COIN_DBL_MAX;
    newObjective[iColumn]=0.0;
  }
  int * start = new int[numberExtraColumns+1];
  int * row = new int[numberExtraElements];
  double * element = new double[numberExtraElements];
  // For objects
  double * weight = new double[largestSet];
  int * sequence = new int [largestSet];
  // get copy 
  CbcSolverLink * newSolver = new CbcSolverLink(*this);
  // replace model
  newSolver->freeCachedResults();
  delete newSolver->modelPtr_;
  newSolver->modelPtr_=new ClpSimplex(*simplexModel_);
  ClpSimplex * newClp = newSolver->modelPtr_;
  // correct bounds
  memcpy(newClp->columnLower(),lower,numberColumns*sizeof(double));
  memcpy(newClp->columnUpper(),upper,numberColumns*sizeof(double));
  // space for bound information
  numberVariables_ = 0;
  delete [] info_;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (back[iColumn]) {
      back[iColumn]=numberVariables_;
      numberVariables_++;
    } else {
      back[iColumn]=-1;
    }
  }
  info_ = new CbcLinkedBound [numberVariables_];
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    int iVar = back[iColumn];
    if (iVar>=0) 
      info_[iVar] = CbcLinkedBound(NULL,iColumn,0,NULL,NULL,NULL);
  }
  // First run just puts in junk and does objects
  start[0]=0;
  int nel=0;
  int n=0;
  for (iSet=0;iSet<numberSets_;iSet++) {
    int functionColumn = functionColumn_[iSet];
    int iSeq=0;
    int iNon;
    for (iNon =startSet_[iSet];iNon<startSet_[iSet+1];iNon++) {
      if (definitionRowL_[iNon]>=0)
	break;
    }
    int startNon=iNon;
    int number = numberInSet[iSet];
    double lo = lower[functionColumn];
    double up = upper[functionColumn];
    double g = granularity[iSet];
    double valueZ=lo;
    for (int i=0;i<number;i++) {
      for (int iNon =startNon;iNon<startSet_[iSet+1];iNon++) {
	// for Link
	weight[iSeq]=valueZ;
	sequence[iSeq++]=n+numberColumns;
	for (int i=0;i<2;i++) {
	  if (iNon==startNon) {
	    // linear terms
	    for (int iNon =startSet_[iSet];iNon<startNon;iNon++) {
	      row[nel]=usageRow_[iNon];
	      element[nel++]=1.0;
	    }
	  }
	  // convexity
	  row[nel]=convexityRow_[iNon];
	  element[nel++]=1.0;
	  // definition of z
	  assert (definitionRow_[iNon]>=0);
	  row[nel]=definitionRow_[iNon];
	  element[nel++]=1.0;
	  if (definitionRowL_[iNon]) {
	    row[nel]=definitionRowL_[iNon];
	    element[nel++]=1.0;
	  }
	  row[nel]=usageRow_[iNon];
	  element[nel++]=1.0;
	  n++;
	  start[n]=nel;
	}
      }
      valueZ += g;
      if (fabs(valueZ-up)<1.0e-7)
	valueZ=up;
    }
    int nLink = iSeq/number;
    objects[iSet]=new CbcLink(NULL,number,nLink,setType[iSet],sequence,weight,iSet);
    objects[iSet]->setPriority(1000);
    objects[iSet]->setPriority(10000);
    
  }
  newClp->addColumns(n,newLower,newUpper,newObjective,
		     start,row,element);
  delete [] newLower;
  delete [] newUpper;
  delete [] newObjective;
  delete [] start;
  delete [] row;
  delete [] element;
  lower = newClp->columnLower();
  upper = newClp->columnUpper();
  // pack down
  CoinPackedMatrix * matrix = newClp->matrix();
  matrix->removeGaps();
  row = matrix->getMutableIndices();
  start = matrix->getMutableVectorStarts();
  element = matrix->getMutableElements();
  // position
  n = numberColumns;
  nel = start[n];
  // Now do real work
  for (iSet=0;iSet<numberSets_;iSet++) {
    int functionColumn = functionColumn_[iSet];
    const char * name = coinModel_.getColumnName(functionColumn);
    int iSeq=0;
    int iNon;
    for (iNon =startSet_[iSet];iNon<startSet_[iSet+1];iNon++) {
      if (definitionRowL_[iNon]>=0)
	break;
    }
    int startNon=iNon;
    int number = numberInSet[iSet];
    double lo = lower[functionColumn];
    double up = upper[functionColumn];
    double g = granularity[iSet];
    double valueZ=lo;
    for (int i=0;i<number;i++) {
      for (int iNon =startNon;iNon<startSet_[iSet+1];iNon++) {
	// for Link
	weight[iSeq]=valueZ;
	coinModel_.associateElement(name,valueZ);
	sequence[iSeq++]=n+numberColumns;
	for (int i=0;i<2;i++) {
	  if (iNon==startNon) {
	    // linear terms
	    for (int iNon =startSet_[iSet];iNon<startNon;iNon++) {
	      row[nel]=usageRow_[iNon];
	      element[nel++]=coinModel_.getElement(usageRow_[iNon],functionColumn)*valueZ;
	    }
	  }
	  // convexity
	  row[nel]=convexityRow_[iNon];
	  element[nel++]=1.0;
	  // definition of z
	  assert (definitionRow_[iNon]>=0);
	  row[nel]=definitionRow_[iNon];
	  element[nel++]=valueZ;
	  if (definitionRowL_[iNon]) {
	    row[nel]=definitionRowL_[iNon];
	    int jColumn = definitionColumnL_[iNon];
	    int jVar = back[jColumn];
	    assert (jVar>=0);
	    double bound = (!i) ? lower[jColumn] :upper[jColumn];
	    info_[jVar].addCoefficientModifier((i!=0),nel,1.0);
	    element[nel++]=bound;
	    row[nel]=usageRow_[iNon];
#if 0
	    CoinModel dummy;
	    dummy.setElement(0,0,coinModel_.getColumnName(jColumn));
	    dummy.setElement(0,1,coinModel_.getElementAsString(row[nel],functionColumn));
	    double associated[2];
	    associated[0]=valueZ;
	    associated[1]=dummy.unsetValue();
	    dummy.computeAssociated(associated);
	    double value = associated[1];
#endif
	    double value = getFunctionValueFromString(coinModel_.getElementAsString(row[nel],functionColumn),
					       name,valueZ);
	    info_[jVar].addCoefficientModifier((i!=0),nel,value);
	    element[nel++]=bound*value;
	  } else {
	    // just z
	    row[nel]=usageRow_[iNon];
	    CoinModel dummy;
	    dummy.setElement(0,0,name);
	    dummy.setElement(0,1,coinModel_.getElementAsString(row[nel],functionColumn));
	    double associated[2];
	    associated[0]=valueZ;
	    associated[1]=dummy.unsetValue();
	    dummy.computeAssociated(associated);
	    double value = associated[1];
	    element[nel++]=value;
	    element[nel++]=1.0;
	  }
	  n++;
	  start[n]=nel;
	}
      }
      valueZ += g;
      if (fabs(valueZ-up)<1.0e-7)
	valueZ=up;
    }
  }
  matrix_ = new CoinPackedMatrix(*matrix);
  delete [] weight;
  delete [] sequence;
  delete [] numberInSet;
  delete [] setType;
  delete [] granularity;
  // now load up model
  assert (!savedModel_);
  ClpSimplex * saved = modelPtr_;
  modelPtr_=newClp;
  OsiSolverInterface * solver = this->clone();
  savedModel_ = saved;
  model_->assignSolver(solver);
  for (iSet=0;iSet<numberSets_;iSet++) {
    objects[iSet]->setModel(model_);
  }
  model_->addObjects(numberSets_,objects);
  for (iSet=0;iSet<numberSets_;iSet++) {
    delete objects[iSet];
  }
  delete [] objects;
}
// revert to previous solver
void 
CbcSolverLink::revert()
{
  delete modelPtr_;
  modelPtr_=savedModel_;
  savedModel_=NULL;
}
// Add a bound modifier
void 
CbcSolverLink::addBoundModifier(bool upperBoundAffected, bool useUpperBound, int whichVariable, int whichVariableAffected, 
				double multiplier)
{
  bool found=false;
  int i;
  for ( i=0;i<numberVariables_;i++) {
    if (info_[i].variable()==whichVariable) {
      found=true;
      break;
    }
  }
  if (!found) {
    // add in
    CbcLinkedBound * temp = new CbcLinkedBound [numberVariables_+1];
    for (int i=0;i<numberVariables_;i++) 
      temp[i]= info_[i];
    delete [] info_;
    info_=temp;
    info_[numberVariables_++] = CbcLinkedBound(model_,whichVariable,0,NULL,NULL,NULL);
  }
  info_[i].addBoundModifier(upperBoundAffected, useUpperBound,whichVariableAffected,multiplier);
}
// Sets granularity
void 
CbcSolverLink::setGranularity(double top, double bottom, const double * topG, const double * bottomG)
{
  delete [] granularityTop_;
  delete [] granularityBottom_;
  granularityTop_ = NULL;
  granularityBottom_ = NULL;
  if (numberSets_>0) {
    if (topG) {
      granularityTop_ = CoinCopyOfArray(topG,numberSets_);
    } else {
      granularityTop_ = new double [numberSets_];
      for (int i=0;i<numberSets_;i++) {
	granularityTop_[i]=top;
      }
    }
    if (bottomG) {
      granularityBottom_ = CoinCopyOfArray(bottomG,numberSets_);
    } else {
      granularityBottom_ = new double [numberSets_];
      for (int i=0;i<numberSets_;i++) {
	granularityBottom_[i]=bottom;
      }
    }
  }
}
//#############################################################################
// Constructors, destructors  and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CbcLinkedBound::CbcLinkedBound ()
{
  model_ = NULL;
  variable_ = -1;
  numberAffected_ = 0;
  maximumAffected_ = numberAffected_;
  affected_ = NULL;
  useObject_=-1;
  numberSets_=0;
  whichSets_ = NULL;
}
// Useful Constructor
CbcLinkedBound::CbcLinkedBound(CbcModel * model, int variable,
			       int numberAffected, const int * positionL, 
			       const int * positionU, const double * multiplier)
{
  model_ = model;
  variable_ = variable;
  numberAffected_ = 2*numberAffected;
  maximumAffected_ = numberAffected_;
  numberSets_ = 0;
  whichSets_ = NULL;
  useObject_=-1;
  if (numberAffected_) { 
    affected_ = new boundElementAction[numberAffected_];
    int n=0;
    for (int i=0;i<numberAffected;i++) {
      // LB
      boundElementAction action;
      action.affect=2;
      action.ubUsed=0;
      action.type=0;
      action.affected=positionL[i];
      action.multiplier=multiplier[i];
      affected_[n++]=action;
      // UB
      action.affect=2;
      action.ubUsed=1;
      action.type=0;
      action.affected=positionU[i];
      action.multiplier=multiplier[i];
      affected_[n++]=action;
    }
  } else {
    affected_ = NULL;
  }
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CbcLinkedBound::CbcLinkedBound (
                  const CbcLinkedBound & rhs)
{
  model_ = rhs.model_;
  variable_ = rhs.variable_;
  numberAffected_ = rhs.numberAffected_;
  maximumAffected_ = rhs.maximumAffected_;
  numberSets_ = rhs.numberSets_;
  useObject_=rhs.useObject_;
  if (numberAffected_) { 
    affected_ = new boundElementAction[maximumAffected_];
    memcpy(affected_,rhs.affected_,numberAffected_*sizeof(boundElementAction));
    whichSets_ = CoinCopyOfArray(rhs.whichSets_,numberSets_);
  } else {
    affected_ = NULL;
    whichSets_ = NULL;
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CbcLinkedBound::~CbcLinkedBound ()
{
  delete [] affected_;
  delete [] whichSets_;
}

//-------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CbcLinkedBound &
CbcLinkedBound::operator=(const CbcLinkedBound& rhs)
{
  if (this != &rhs) { 
    delete [] affected_;
    delete [] whichSets_;
    model_ = rhs.model_;
    variable_ = rhs.variable_;
    numberAffected_ = rhs.numberAffected_;
    maximumAffected_ = rhs.maximumAffected_;
    numberSets_ = rhs.numberSets_;
    useObject_=rhs.useObject_;
    if (numberAffected_) { 
      affected_ = new boundElementAction[maximumAffected_];
      memcpy(affected_,rhs.affected_,numberAffected_*sizeof(boundElementAction));
      whichSets_ = CoinCopyOfArray(rhs.whichSets_,numberSets_);
    } else {
      affected_ = NULL;
      whichSets_ = NULL;
    }
  }
  return *this;
}
// Add a bound modifier
void 
CbcLinkedBound::addBoundModifier(bool upperBoundAffected, bool useUpperBound, int whichVariable, 
				 double multiplier)
{
  if (numberAffected_==maximumAffected_) {
    maximumAffected_ = maximumAffected_+10+maximumAffected_/4;
    boundElementAction * temp = new boundElementAction[maximumAffected_];
    memcpy(temp,affected_,numberAffected_*sizeof(boundElementAction));
    delete [] affected_;
    affected_ = temp;
  }
  boundElementAction action;
  action.affect=upperBoundAffected ? 1 : 0;
  action.ubUsed=useUpperBound ? 1 : 0;
  action.type=2;
  action.affected=whichVariable;
  action.multiplier=multiplier;
  affected_[numberAffected_++]=action;
  
}
// Update other bounds
void 
CbcLinkedBound::updateBounds(ClpSimplex * solver)
{
  double * lower = solver->columnLower();
  double * upper = solver->columnUpper();
  double lo = lower[variable_];
  double up = upper[variable_];
  // printf("bounds for %d are %g and %g\n",variable_,lo,up);
  for (int j=0;j<numberAffected_;j++) {
    if (affected_[j].affect<2) {
      double multiplier = affected_[j].multiplier;
      assert (affected_[j].type==2);
      int iColumn = affected_[j].affected;
      double useValue = (affected_[j].ubUsed) ? up : lo;
      if (affected_[j].affect==0) 
	lower[iColumn] = CoinMin(upper[iColumn],CoinMax(lower[iColumn],multiplier*useValue));
      else
	upper[iColumn] = CoinMax(lower[iColumn],CoinMin(upper[iColumn],multiplier*useValue));
    }
  }
}
// Add an element modifier
void 
CbcLinkedBound::addCoefficientModifier(bool useUpperBound, int position, 
				 double multiplier)
{
  if (numberAffected_==maximumAffected_) {
    maximumAffected_ = maximumAffected_+10+maximumAffected_/4;
    boundElementAction * temp = new boundElementAction[maximumAffected_];
    memcpy(temp,affected_,numberAffected_*sizeof(boundElementAction));
    delete [] affected_;
    affected_ = temp;
  }
  boundElementAction action;
  action.affect=2;
  action.ubUsed=useUpperBound ? 1 : 0;
  action.type=0;
  action.affected=position;
  action.multiplier=multiplier;
  affected_[numberAffected_++]=action;
  
}
// Update coefficients
void 
CbcLinkedBound::updateCoefficients(ClpSimplex * solver, CoinPackedMatrix * matrix)
{
  double * lower = solver->columnLower();
  double * upper = solver->columnUpper();
  double * element = matrix->getMutableElements();
  double lo = lower[variable_];
  double up = upper[variable_];
  // printf("bounds for %d are %g and %g\n",variable_,lo,up);
  for (int j=0;j<numberAffected_;j++) {
    if (affected_[j].affect==2) {
      double multiplier = affected_[j].multiplier;
      assert (affected_[j].type==0);
      int position = affected_[j].affected;
      //double old = element[position];
      if (affected_[j].ubUsed)
	element[position] = multiplier*up;
      else
	element[position] = multiplier*lo;
      //if ( old != element[position])
      //printf("change at %d from %g to %g\n",position,old,element[position]);
    }
  }
}
// Fix variables in LOS
void 
CbcLinkedBound::fixLOS(ClpSimplex * solver)
{
  abort();
}
// Update this variables bounds from reference row
void 
CbcLinkedBound::updateBoundsForThis(ClpSimplex * solver)
{
  if (useObject_>=0) {
    const OsiObject * object = model_->object(useObject_);
    const CbcLink * link = dynamic_cast<const CbcLink *> (object);
    assert (link);
    double * lower = solver->columnLower();
    double * upper = solver->columnUpper();
    double lo = lower[variable_];
    double up = upper[variable_];
    int numberMembers = link->numberMembers();
    int numberLinks = link->numberLinks();
    const int * which = link->which();
    const double * weights = link->weights();
    int i;
    double loWeight=-COIN_DBL_MAX;
    double upWeight=COIN_DBL_MAX;
    int base=0;
    for (i=0;i<numberMembers;i++) {
      int nFixed=0;
      for (int j=0;j<numberLinks;j++) {
	int iColumn = which[j+base];
	if (!upper[iColumn])
	  nFixed++;
      }
      if (nFixed) {
	assert (nFixed==numberLinks);
	if (loWeight!=-COIN_DBL_MAX) {
	  assert(i>0);
	  upWeight=weights[i-1];
	  break;
	}
      } else if (loWeight==-COIN_DBL_MAX) {
	loWeight=weights[i];
      }
    }
    upper[variable_] = CoinMin(up,upWeight);
    lower[variable_] = CoinMax(lo,loWeight);
    assert (upper[variable_]>=lower[variable_]);
  }
}
// Default Constructor
CbcHeuristicDynamic2::CbcHeuristicDynamic2() 
  :CbcHeuristic()
{
}

// Constructor from model
CbcHeuristicDynamic2::CbcHeuristicDynamic2(CbcModel & model)
  :CbcHeuristic(model)
{
}

// Destructor 
CbcHeuristicDynamic2::~CbcHeuristicDynamic2 ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicDynamic2::clone() const
{
  return new CbcHeuristicDynamic2(*this);
}

// Copy constructor 
CbcHeuristicDynamic2::CbcHeuristicDynamic2(const CbcHeuristicDynamic2 & rhs)
:
  CbcHeuristic(rhs)
{
}

// Returns 1 if solution, 0 if not
int
CbcHeuristicDynamic2::solution(double & solutionValue,
			 double * betterSolution)
{
  if (!model_)
    return 0;
  CbcSolverLink * clpSolver 
    = dynamic_cast<CbcSolverLink *> (model_->solver());
  assert (clpSolver); 
  double newSolutionValue = clpSolver->bestObjectiveValue();
  const double * solution = clpSolver->bestSolution();
  if (newSolutionValue<solutionValue&&solution) {
    int numberColumns = clpSolver->getNumCols();
    // new solution
    memcpy(betterSolution,solution,numberColumns*sizeof(double));
    solutionValue = newSolutionValue;
    return 1;
  } else {
    return 0;
  }
}
// update model
void CbcHeuristicDynamic2::setModel(CbcModel * model)
{
  model_ = model;
}
// Resets stuff if model changes
void 
CbcHeuristicDynamic2::resetModel(CbcModel * model)
{
  model_ = model;
}
