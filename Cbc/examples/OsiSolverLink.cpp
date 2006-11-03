// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <cassert>

#include "CoinTime.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinMpsIO.hpp"
#include "CoinModel.hpp"
#include "ClpSimplex.hpp"
#include "OsiSolverLink.hpp"
#include "OsiBranchLink.hpp"
#include "ClpPackedMatrix.hpp"
#include "CoinTime.hpp"
//#############################################################################
// Solve methods
//#############################################################################
void OsiSolverLink::initialSolve()
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
    updateCoefficients(modelPtr_,temp);
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
//#define WRITE_MATRIX
#ifdef WRITE_MATRIX
static int xxxxxx=0;
#endif
//-----------------------------------------------------------------------------
void OsiSolverLink::resolve()
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
    updateCoefficients(modelPtr_,temp);
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
#ifdef WRITE_MATRIX
  {
    xxxxxx++;
    char temp[50];
    sprintf(temp,"bb%d",xxxxxx);
    writeMps(temp);
  }
#endif
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
OsiSolverLink::OsiSolverLink ()
  : OsiClpSolverInterface()
{
  gutsOfDestructor(true);
}
#if 0
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
#endif
/* This creates from a coinModel object 

   if errors.then number of sets is -1
      
   This creates linked ordered sets information.  It assumes -

   for product terms syntax is yy*f(zz)
   also just f(zz) is allowed 
   and even a constant

   modelObject not const as may be changed as part of process.
*/
OsiSolverLink::OsiSolverLink ( CoinModel & coinModel)
  : OsiClpSolverInterface()
{
  gutsOfDestructor(true);
  gdb(coinModel);
}
void OsiSolverLink::gdb ( CoinModel & coinModel)
{
  // first check and set up arrays
  int numberColumns = coinModel.numberColumns();
  int numberRows = coinModel.numberRows();
  // List of nonlinear entries
  int * which = new int[numberColumns];
  numberVariables_=0;
  specialOptions2_=0;
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
      }
      triple=coinModel.next(triple);
      n++;
    }
    if (!linear) {
      which[numberVariables_++]=iColumn;
    }
  }
  // return if nothing
  if (!numberVariables_) {
    delete [] which;
    return;
  } else {
    coinModel_ = coinModel;
    int nBi=0;
    int iRow;
    for (iRow=0;iRow<numberRows;iRow++) {   
      CoinModelLink triple=coinModel_.firstInRow(iRow);
      while (triple.column()>=0) {
	int iColumn = triple.column();
	const char *  el = coinModel_.getElementAsString(iRow,iColumn);
	if (strcmp("Numeric",el)) {
	  // check if value*x*y
	  char temp[100];
	  strcpy(temp,el);
	  char * ast = strchr(temp,'*');
	  if (ast) {
	    char * pos = temp;
	    while (pos!=ast) {
	      char x = *pos;
	      pos++;
	      assert ((x>='0'&&x<='9')||x=='.'||x=='+'||x=='-');
	    }
	    // must be column
	    assert (coinModel_.column(ast+1)>=0);
	  } else {
	    // must be column
	    assert (coinModel_.column(temp)>=0);
	  }
	  coinModel.setElement(iRow,iColumn,0.0);
	  nBi++;
	}
	triple=coinModel_.next(triple);
      }
    }
    if (!nBi)
      exit(1);
    int nInt=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {   
      if (coinModel_.isInteger(iColumn))
	nInt++;
    }
    printf("There are %d bilinear and %d integers\n",nBi,nInt);
    loadFromCoinModel(coinModel,true);
    OsiObject ** objects = new OsiObject * [nBi+nInt];
    char * marked = new char [numberColumns];
    memset(marked,0,numberColumns);
    nBi=nInt;
    for (iRow=0;iRow<numberRows;iRow++) {   
      CoinModelLink triple=coinModel_.firstInRow(iRow);
      while (triple.column()>=0) {
	int iColumn = triple.column();
	const char *  el = coinModel_.getElementAsString(iRow,iColumn);
	if (strcmp("Numeric",el)) {
	  // value*x*y
	  char temp[100];
	  strcpy(temp,el);
	  char * ast = strchr(temp,'*');
	  double value=1.0;
	  if (ast) {
	    *ast='\0';
	    value = atof(temp);
	    ast++;
	  } else {
	    ast=temp;
	  }
	  // other column
	  int jColumn = coinModel_.column(ast);
	  double meshI = coinModel_.isInteger(iColumn) ? 1.0 : 0.0;
	  if (meshI)
	    marked[iColumn]=1;
	  double meshJ = coinModel_.isInteger(jColumn) ? 1.0 : 0.0;
	  if (meshJ)
	    marked[jColumn]=1;
	  if (!meshJ&&!meshI) {
	    meshI=5.0e-1;
	    meshJ=5.0e-1;
	  }
	  objects[nBi] = new OsiBiLinear(this,iColumn,jColumn,iRow,value,meshI,meshJ,
					 nBi-nInt,(const OsiObject **) (objects+nInt));
	  objects[nBi]->setPriority(10000);
	  nBi++;
	}
	triple=coinModel_.next(triple);
      }
    }
    nInt=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {   
      if (coinModel_.isInteger(iColumn)) {
	objects[nInt] = new OsiSimpleInteger(this,iColumn);
	if (marked[iColumn])
	  objects[nInt]->setPriority(100000);
	else
	  objects[nInt]->setPriority(100);
	nInt++;
      }
    }
    nInt=nBi;
    delete [] marked;
    if (numberErrors) {
      // errors
      gutsOfDestructor();
      numberVariables_=-1;
    } else {
      addObjects(nInt,objects);
      int i;
      for (i=0;i<nInt;i++)
	delete objects[i];
      delete [] objects;
      // Now do dummy bound stuff
      matrix_ = new CoinPackedMatrix(*getMatrixByCol());
      info_ = new OsiLinkedBound [numberVariables_];
      for ( i=0;i<numberVariables_;i++) {
	info_[i] = OsiLinkedBound(this,which[i],0,NULL,NULL,NULL);
      }
    }
  }
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface * 
OsiSolverLink::clone(bool copyData) const
{
  assert (copyData);
  return new OsiSolverLink(*this);
}


//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiSolverLink::OsiSolverLink (
                  const OsiSolverLink & rhs)
  : OsiClpSolverInterface(rhs)
{
  gutsOfDestructor(true);
  gutsOfCopy(rhs);
  // something odd happens - try this
  OsiSolverInterface::operator=(rhs);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiSolverLink::~OsiSolverLink ()
{
  gutsOfDestructor();
}

//-------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiSolverLink &
OsiSolverLink::operator=(const OsiSolverLink& rhs)
{
  if (this != &rhs) { 
    gutsOfDestructor();
    OsiClpSolverInterface::operator=(rhs);
    gutsOfCopy(rhs);
  }
  return *this;
}
void 
OsiSolverLink::gutsOfDestructor(bool justNullify)
{
  if (!justNullify) {
    delete matrix_;
    delete [] info_;
  } 
  matrix_ = NULL;
  info_ = NULL;
  numberVariables_ = 0;
  specialOptions2_ = 0;
}
void 
OsiSolverLink::gutsOfCopy(const OsiSolverLink & rhs)
{
  coinModel_ = rhs.coinModel_;
  numberVariables_ = rhs.numberVariables_;
  specialOptions2_ = rhs.specialOptions2_;
  if (numberVariables_) { 
    matrix_ = new CoinPackedMatrix(*rhs.matrix_);
    info_ = new OsiLinkedBound [numberVariables_];
    for (int i=0;i<numberVariables_;i++) {
      info_[i] = OsiLinkedBound(rhs.info_[i]);
    }
  }
}
// Add a bound modifier
void 
OsiSolverLink::addBoundModifier(bool upperBoundAffected, bool useUpperBound, int whichVariable, int whichVariableAffected, 
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
    OsiLinkedBound * temp = new OsiLinkedBound [numberVariables_+1];
    for (int i=0;i<numberVariables_;i++) 
      temp[i]= info_[i];
    delete [] info_;
    info_=temp;
    info_[numberVariables_++] = OsiLinkedBound(this,whichVariable,0,NULL,NULL,NULL);
  }
  info_[i].addBoundModifier(upperBoundAffected, useUpperBound,whichVariableAffected,multiplier);
}
// Update coefficients
void 
OsiSolverLink::updateCoefficients(ClpSimplex * solver, CoinPackedMatrix * matrix)
{
  double * lower = solver->columnLower();
  double * upper = solver->columnUpper();
  for (int iObject =0;iObject<numberObjects_;iObject++) {
    OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[iObject]);
    if (obj) {
      obj->updateCoefficients(lower,upper,matrix,&basis_);
    }
  }
}
//#############################################################################
// Constructors, destructors  and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiLinkedBound::OsiLinkedBound ()
{
  model_ = NULL;
  variable_ = -1;
  numberAffected_ = 0;
  maximumAffected_ = numberAffected_;
  affected_ = NULL;
}
// Useful Constructor
OsiLinkedBound::OsiLinkedBound(OsiSolverInterface * model, int variable,
			       int numberAffected, const int * positionL, 
			       const int * positionU, const double * multiplier)
{
  model_ = model;
  variable_ = variable;
  numberAffected_ = 2*numberAffected;
  maximumAffected_ = numberAffected_;
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
OsiLinkedBound::OsiLinkedBound (
                  const OsiLinkedBound & rhs)
{
  model_ = rhs.model_;
  variable_ = rhs.variable_;
  numberAffected_ = rhs.numberAffected_;
  maximumAffected_ = rhs.maximumAffected_;
  if (numberAffected_) { 
    affected_ = new boundElementAction[maximumAffected_];
    memcpy(affected_,rhs.affected_,numberAffected_*sizeof(boundElementAction));
  } else {
    affected_ = NULL;
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiLinkedBound::~OsiLinkedBound ()
{
  delete [] affected_;
}

//-------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiLinkedBound &
OsiLinkedBound::operator=(const OsiLinkedBound& rhs)
{
  if (this != &rhs) { 
    delete [] affected_;
    model_ = rhs.model_;
    variable_ = rhs.variable_;
    numberAffected_ = rhs.numberAffected_;
    maximumAffected_ = rhs.maximumAffected_;
    if (numberAffected_) { 
      affected_ = new boundElementAction[maximumAffected_];
      memcpy(affected_,rhs.affected_,numberAffected_*sizeof(boundElementAction));
    } else {
      affected_ = NULL;
    }
  }
  return *this;
}
// Add a bound modifier
void 
OsiLinkedBound::addBoundModifier(bool upperBoundAffected, bool useUpperBound, int whichVariable, 
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
OsiLinkedBound::updateBounds(ClpSimplex * solver)
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
#if 0
// Add an element modifier
void 
OsiLinkedBound::addCoefficientModifier(bool useUpperBound, int position, 
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
OsiLinkedBound::updateCoefficients(ClpSimplex * solver, CoinPackedMatrix * matrix)
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
#endif
