// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

// Default Constructor 
CbcClique::CbcClique ()
  : CbcObject(),
    numberMembers_(0),
    numberNonSOSMembers_(0),
    members_(NULL),
    type_(NULL),
    cliqueType_(-1),
    slack_(-1)
{
}

// Useful constructor (which are integer indices)
CbcClique::CbcClique (CbcModel * model, int cliqueType, int numberMembers,
	   const int * which, const char * type, int identifier,int slack)
  : CbcObject(model)
{
  id_=identifier;
  numberMembers_=numberMembers;
  if (numberMembers_) {
    members_ = new int[numberMembers_];
    memcpy(members_,which,numberMembers_*sizeof(int));
    type_ = new char[numberMembers_];
    if (type) {
      memcpy(type_,type,numberMembers_*sizeof(char));
    } else {
      for (int i=0;i<numberMembers_;i++)
	type_[i]=1;
    }
  } else {
    members_ = NULL;
    type_ = NULL;
  }
  // Find out how many non sos
  int i;
  numberNonSOSMembers_=0;
  for (i=0;i<numberMembers_;i++)
    if (!type_[i])
      numberNonSOSMembers_++;
  cliqueType_ = cliqueType;
  slack_ = slack;
}

// Copy constructor 
CbcClique::CbcClique ( const CbcClique & rhs)
  :CbcObject(rhs)
{
  numberMembers_ = rhs.numberMembers_;
  numberNonSOSMembers_ = rhs.numberNonSOSMembers_;
  if (numberMembers_) {
    members_ = new int[numberMembers_];
    memcpy(members_,rhs.members_,numberMembers_*sizeof(int));
    type_ = new char[numberMembers_];
    memcpy(type_,rhs.type_,numberMembers_*sizeof(char));
  } else {
    members_ = NULL;
    type_ = NULL;
  }
  cliqueType_ = rhs.cliqueType_;
  slack_ = rhs.slack_;
}

// Clone
CbcObject *
CbcClique::clone() const
{
  return new CbcClique(*this);
}

// Assignment operator 
CbcClique & 
CbcClique::operator=( const CbcClique& rhs)
{
  if (this!=&rhs) {
    CbcObject::operator=(rhs);
    delete [] members_;
    delete [] type_;
    numberMembers_ = rhs.numberMembers_;
    numberNonSOSMembers_ = rhs.numberNonSOSMembers_;
    if (numberMembers_) {
      members_ = new int[numberMembers_];
      memcpy(members_,rhs.members_,numberMembers_*sizeof(int));
      type_ = new char[numberMembers_];
      memcpy(type_,rhs.type_,numberMembers_*sizeof(char));
    } else {
      members_ = NULL;
      type_ = NULL;
    }
    cliqueType_ = rhs.cliqueType_;
    slack_ = rhs.slack_;
  }
  return *this;
}

// Destructor 
CbcClique::~CbcClique ()
{
  delete [] members_;
  delete [] type_;
}

// Infeasibility - large is 0.5
double 
CbcClique::infeasibility(int & preferredWay) const
{
  int numberUnsatis=0, numberFree=0;
  int j;
  const int * integer = model_->integerVariable();
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  double largestValue=0.0;
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double * sort = new double[numberMembers_];

  double slackValue=0.0;
  for (j=0;j<numberMembers_;j++) {
    int sequence = members_[j];
    int iColumn = integer[sequence];
    double value = solution[iColumn];
    value = CoinMax(value, lower[iColumn]);
    value = CoinMin(value, upper[iColumn]);
    double nearest = floor(value+0.5);
    double distance = fabs(value-nearest);
    if (distance>integerTolerance) {
      if (!type_[j])
	value = 1.0-value; // non SOS
      // if slack then choose that
      if (j==slack_&&value>0.05)
	slackValue = value;
      largestValue = CoinMax(value,largestValue);
      sort[numberUnsatis++]=-value;
    } else if (upper[iColumn]>lower[iColumn]) {
      numberFree++;
    }
  }
  preferredWay=1;
  double otherWay = 0.0;
  if (numberUnsatis) {
    // sort
    std::sort(sort,sort+numberUnsatis);
    for (j=0;j<numberUnsatis;j++) {
      if ((j&1)!=0)
	otherWay += -sort[j];
    }
    // Need to think more
    double value = 0.2*numberUnsatis+0.01*(numberMembers_-numberFree);
    if (fabs(largestValue-0.5)<0.1) {
      // close to half
      value +=0.1;
    }
    if (slackValue) {
      // branching on slack
      value += slackValue;
    }
    // scale other way
    otherWay *= value/(1.0-otherWay);
    delete [] sort;
    return value;
  } else {
    delete [] sort;
    return 0.0; // satisfied
  }
}

// This looks at solution and sets bounds to contain solution
void 
CbcClique::feasibleRegion()
{
  int j;
  const int * integer = model_->integerVariable();
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
#ifndef NDEBUG
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
#endif  
  for (j=0;j<numberMembers_;j++) {
    int sequence = members_[j];
    int iColumn = integer[sequence];
    double value = solution[iColumn];
    value = CoinMax(value, lower[iColumn]);
    value = CoinMin(value, upper[iColumn]);
    double nearest = floor(value+0.5);
#ifndef NDEBUG
    double distance = fabs(value-nearest);
    assert(distance<=integerTolerance);
#endif
    solver->setColLower(iColumn,nearest);
    solver->setColUpper(iColumn,nearest);
  }
}
// Redoes data when sequence numbers change
void 
CbcClique::redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns)
{
  model_=model;
  int n2=0;
  for (int j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    int i;
    for (i=0;i<numberColumns;i++) {
      if (originalColumns[i]==iColumn)
        break;
    }
    if (i<numberColumns) {
      members_[n2]=i;
      type_[n2++]=type_[j];
    }
  }
  if (n2<numberMembers_) {
    //printf("** SOS number of members reduced from %d to %d!\n",numberMembers_,n2);
    numberMembers_=n2;
  }
  // Find out how many non sos
  int i;
  numberNonSOSMembers_=0;
  for (i=0;i<numberMembers_;i++)
    if (!type_[i])
      numberNonSOSMembers_++;
}


// Creates a branching object
CbcBranchingObject * 
CbcClique::createBranch(int way) 
{
  int numberUnsatis=0;
  int j;
  int nUp=0;
  int nDown=0;
  int numberFree=numberMembers_;
  const int * integer = model_->integerVariable();
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  int * upList = new int[numberMembers_];
  int * downList = new int[numberMembers_];
  double * sort = new double[numberMembers_];
  double integerTolerance = 
      model_->getDblParam(CbcModel::CbcIntegerTolerance);

  double slackValue=0.0;
  for (j=0;j<numberMembers_;j++) {
    int sequence = members_[j];
    int iColumn = integer[sequence];
    double value = solution[iColumn];
    value = CoinMax(value, lower[iColumn]);
    value = CoinMin(value, upper[iColumn]);
    double nearest = floor(value+0.5);
    double distance = fabs(value-nearest);
    if (distance>integerTolerance) {
      if (!type_[j])
	value = 1.0-value; // non SOS
      // if slack then choose that
      if (j==slack_&&value>0.05)
	slackValue = value;
      value = -value; // for sort
      upList[numberUnsatis]=j;
      sort[numberUnsatis++]=value;
    } else if (upper[iColumn]>lower[iColumn]) {
      upList[--numberFree]=j;
    }
  }
  assert (numberUnsatis);
  if (!slackValue) {
    // sort
    CoinSort_2(sort,sort+numberUnsatis,upList);
    // put first in up etc
    int kWay=1;
    for (j=0;j<numberUnsatis;j++) {
      if (kWay>0) 
	upList[nUp++]=upList[j];
      else
	downList[nDown++]=upList[j];
      kWay = -kWay;
    }
    for (j=numberFree;j<numberMembers_;j++) {
      if (kWay>0)
	upList[nUp++]=upList[j];
      else
	downList[nDown++]=upList[j];
      kWay = -kWay;
    }
  } else {
    // put slack to 0 in first way
    nUp = 1;
    upList[0]=slack_;
    for (j=0;j<numberUnsatis;j++) {
      downList[nDown++]=upList[j];
    }
    for (j=numberFree;j<numberMembers_;j++) {
      downList[nDown++]=upList[j];
    }
  }
  // create object
  CbcBranchingObject * branch;
  if (numberMembers_ <=64)
     branch = new CbcCliqueBranchingObject(model_,this,way,
					 nDown,downList,nUp,upList);
  else
    branch = new CbcLongCliqueBranchingObject(model_,this,way,
					    nDown,downList,nUp,upList);
  delete [] upList;
  delete [] downList;
  delete [] sort;
  return branch;
}

// Default Constructor 
CbcSOS::CbcSOS ()
  : CbcObject(),
    members_(NULL),
    weights_(NULL),
    numberMembers_(0),
    sosType_(-1),
    integerValued_(false)
{
}

// Useful constructor (which are indices)
CbcSOS::CbcSOS (CbcModel * model,  int numberMembers,
	   const int * which, const double * weights, int identifier,int type)
  : CbcObject(model),
    numberMembers_(numberMembers),
    sosType_(type)
{
  id_=identifier;
  integerValued_ = type==1;
  if (numberMembers_) {
    members_ = new int[numberMembers_];
    weights_ = new double[numberMembers_];
    memcpy(members_,which,numberMembers_*sizeof(int));
    if (weights) {
      memcpy(weights_,weights,numberMembers_*sizeof(double));
    } else {
      for (int i=0;i<numberMembers_;i++)
        weights_[i]=i;
    }
    // sort so weights increasing
    CoinSort_2(weights_,weights_+numberMembers_,members_);
    double last = -COIN_DBL_MAX;
    int i;
    for (i=0;i<numberMembers_;i++) {
      double possible = CoinMax(last+1.0e-10,weights_[i]);
      weights_[i] = possible;
      last=possible;
    }
  } else {
    members_ = NULL;
    weights_ = NULL;
  }
  assert (sosType_>0&&sosType_<3);
}

// Copy constructor 
CbcSOS::CbcSOS ( const CbcSOS & rhs)
  :CbcObject(rhs)
{
  numberMembers_ = rhs.numberMembers_;
  sosType_ = rhs.sosType_;
  integerValued_ = rhs.integerValued_;
  if (numberMembers_) {
    members_ = new int[numberMembers_];
    weights_ = new double[numberMembers_];
    memcpy(members_,rhs.members_,numberMembers_*sizeof(int));
    memcpy(weights_,rhs.weights_,numberMembers_*sizeof(double));
  } else {
    members_ = NULL;
    weights_ = NULL;
  }
}

// Clone
CbcObject *
CbcSOS::clone() const
{
  return new CbcSOS(*this);
}

// Assignment operator 
CbcSOS & 
CbcSOS::operator=( const CbcSOS& rhs)
{
  if (this!=&rhs) {
    CbcObject::operator=(rhs);
    delete [] members_;
    delete [] weights_;
    numberMembers_ = rhs.numberMembers_;
    sosType_ = rhs.sosType_;
    integerValued_ = rhs.integerValued_;
    if (numberMembers_) {
      members_ = new int[numberMembers_];
      weights_ = new double[numberMembers_];
      memcpy(members_,rhs.members_,numberMembers_*sizeof(int));
      memcpy(weights_,rhs.weights_,numberMembers_*sizeof(double));
    } else {
      members_ = NULL;
      weights_ = NULL;
    }
  }
  return *this;
}

// Destructor 
CbcSOS::~CbcSOS ()
{
  delete [] members_;
  delete [] weights_;
}

// Infeasibility - large is 0.5
double 
CbcSOS::infeasibility(int & preferredWay) const
{
  int j;
  int firstNonZero=-1;
  int lastNonZero = -1;
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  //const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  //double largestValue=0.0;
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double weight = 0.0;
  double sum =0.0;

  // check bounds etc
  double lastWeight=-1.0e100;
  for (j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    if (lastWeight>=weights_[j]-1.0e-7)
      throw CoinError("Weights too close together in SOS","infeasibility","CbcSOS");
    double value = CoinMax(0.0,solution[iColumn]);
    sum += value;
    if (value>integerTolerance&&upper[iColumn]) {
      // Possibly due to scaling a fixed variable might slip through
      if (value>upper[iColumn]) {
        value=upper[iColumn];
	// Could change to #ifdef CBC_DEBUG
#ifndef NDEBUG
	if (model_->messageHandler()->logLevel()>2)
	  printf("** Variable %d (%d) has value %g and upper bound of %g\n",
		 iColumn,j,value,upper[iColumn]);
#endif
      } 
      weight += weights_[j]*value;
      if (firstNonZero<0)
        firstNonZero=j;
      lastNonZero=j;
    }
  }
  preferredWay=1;
  if (lastNonZero-firstNonZero>=sosType_) {
    // find where to branch
    assert (sum>0.0);
    weight /= sum;
    //int iWhere;
    //for (iWhere=firstNonZero;iWhere<lastNonZero;iWhere++) 
    //if (weight<weights_[iWhere+1])
    //break;
    // probably best to use pseudo duals
    double value = lastNonZero-firstNonZero+1;
    value *= 0.5/((double) numberMembers_);
    return value;
  } else {
    return 0.0; // satisfied
  }
}

// This looks at solution and sets bounds to contain solution
void 
CbcSOS::feasibleRegion()
{
  int j;
  int firstNonZero=-1;
  int lastNonZero = -1;
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  //const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double weight = 0.0;
  double sum =0.0;

  for (j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    double value = CoinMax(0.0,solution[iColumn]);
    sum += value;
    if (value>integerTolerance&&upper[iColumn]) {
      weight += weights_[j]*value;
      if (firstNonZero<0)
	firstNonZero=j;
      lastNonZero=j;
    }
  }
  assert (lastNonZero-firstNonZero<sosType_) ;
  for (j=0;j<firstNonZero;j++) {
    int iColumn = members_[j];
    solver->setColUpper(iColumn,0.0);
  }
  for (j=lastNonZero+1;j<numberMembers_;j++) {
    int iColumn = members_[j];
    solver->setColUpper(iColumn,0.0);
  }
}
// Redoes data when sequence numbers change
void 
CbcSOS::redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns)
{
  model_=model;
  int n2=0;
  for (int j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    int i;
    for (i=0;i<numberColumns;i++) {
      if (originalColumns[i]==iColumn)
        break;
    }
    if (i<numberColumns) {
      members_[n2]=i;
      weights_[n2++]=weights_[j];
    }
  }
  if (n2<numberMembers_) {
    //printf("** SOS number of members reduced from %d to %d!\n",numberMembers_,n2);
    numberMembers_=n2;
  }
}


// Creates a branching object
CbcBranchingObject * 
CbcSOS::createBranch(int way) 
{
  int j;
  const double * solution = model_->testSolution();
  double integerTolerance = 
      model_->getDblParam(CbcModel::CbcIntegerTolerance);
  OsiSolverInterface * solver = model_->solver();
  const double * upper = solver->getColUpper();
  int firstNonFixed=-1;
  int lastNonFixed=-1;
  int firstNonZero=-1;
  int lastNonZero = -1;
  double weight = 0.0;
  double sum =0.0;
  for (j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    if (upper[iColumn]) {
      double value = CoinMax(0.0,solution[iColumn]);
      sum += value;
      if (firstNonFixed<0)
	firstNonFixed=j;
      lastNonFixed=j;
      if (value>integerTolerance) {
	weight += weights_[j]*value;
	if (firstNonZero<0)
	  firstNonZero=j;
	lastNonZero=j;
      }
    }
  }
  assert (lastNonZero-firstNonZero>=sosType_) ;
  // find where to branch
  assert (sum>0.0);
  weight /= sum;
  int iWhere;
  double separator=0.0;
  for (iWhere=firstNonZero;iWhere<lastNonZero;iWhere++) 
    if (weight<weights_[iWhere+1])
      break;
  if (sosType_==1) {
    // SOS 1
    separator = 0.5 *(weights_[iWhere]+weights_[iWhere+1]);
  } else {
    // SOS 2
    if (iWhere==firstNonFixed)
      iWhere++;;
    if (iWhere==lastNonFixed-1)
      iWhere = lastNonFixed-2;
    separator = weights_[iWhere+1];
  }
  // create object
  CbcBranchingObject * branch;
  branch = new CbcSOSBranchingObject(model_,this,way,separator);
  return branch;
}

/* Create an OsiSolverBranch object
   
This returns NULL if branch not represented by bound changes
*/
OsiSolverBranch * 
CbcSOS::solverBranch() const
{
  int j;
  const double * solution = model_->testSolution();
  double integerTolerance = 
      model_->getDblParam(CbcModel::CbcIntegerTolerance);
  OsiSolverInterface * solver = model_->solver();
  const double * upper = solver->getColUpper();
  int firstNonFixed=-1;
  int lastNonFixed=-1;
  int firstNonZero=-1;
  int lastNonZero = -1;
  double weight = 0.0;
  double sum =0.0;
  double * fix = new double[numberMembers_];
  int * which = new int[numberMembers_];
  for (j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    // fix all on one side or other (even if fixed)
    fix[j]=0.0;
    which[j]=iColumn;
    if (upper[iColumn]) {
      double value = CoinMax(0.0,solution[iColumn]);
      sum += value;
      if (firstNonFixed<0)
	firstNonFixed=j;
      lastNonFixed=j;
      if (value>integerTolerance) {
	weight += weights_[j]*value;
	if (firstNonZero<0)
	  firstNonZero=j;
	lastNonZero=j;
      }
    }
  }
  assert (lastNonZero-firstNonZero>=sosType_) ;
  // find where to branch
  assert (sum>0.0);
  weight /= sum;
  // down branch fixes ones above weight to 0
  int iWhere;
  int iDownStart=0;
  int iUpEnd=0;
  for (iWhere=firstNonZero;iWhere<lastNonZero;iWhere++) 
    if (weight<weights_[iWhere+1])
      break;
  if (sosType_==1) {
    // SOS 1
    iUpEnd=iWhere+1;
    iDownStart=iUpEnd;
  } else {
    // SOS 2
    if (iWhere==firstNonFixed)
      iWhere++;;
    if (iWhere==lastNonFixed-1)
      iWhere = lastNonFixed-2;
    iUpEnd=iWhere+1;
    iDownStart=iUpEnd+1;
  }
  // 
  OsiSolverBranch * branch = new OsiSolverBranch();
  branch->addBranch(-1,0,NULL,NULL,numberMembers_-iDownStart,which+iDownStart,fix);
  branch->addBranch(1,0,NULL,NULL,iUpEnd,which,fix);
  delete [] fix;
  delete [] which;
  return branch;
}
// Construct an OsiSOS object
OsiSOS * 
CbcSOS::osiObject(const OsiSolverInterface * solver) const
{
  OsiSOS * obj = new OsiSOS(solver,numberMembers_,members_,weights_,sosType_);
  obj->setPriority(priority());
  return obj;
}

/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
CbcSimpleInteger::CbcSimpleInteger ()
  : CbcObject(),
    originalLower_(0.0),
    originalUpper_(1.0),
    breakEven_(0.5),
    columnNumber_(-1),
    preferredWay_(0)
{
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
CbcSimpleInteger::CbcSimpleInteger ( CbcModel * model, int iColumn, double breakEven)
  : CbcObject(model)
{
  columnNumber_ = iColumn ;
  originalLower_ = model->solver()->getColLower()[columnNumber_] ;
  originalUpper_ = model->solver()->getColUpper()[columnNumber_] ;
  breakEven_ = breakEven;
  assert (breakEven_>0.0&&breakEven_<1.0);
  preferredWay_ = 0;
}


// Copy constructor 
CbcSimpleInteger::CbcSimpleInteger ( const CbcSimpleInteger & rhs)
  :CbcObject(rhs)

{
  columnNumber_ = rhs.columnNumber_;
  originalLower_ = rhs.originalLower_;
  originalUpper_ = rhs.originalUpper_;
  breakEven_ = rhs.breakEven_;
  preferredWay_ = rhs.preferredWay_;
}

// Clone
CbcObject *
CbcSimpleInteger::clone() const
{
  return new CbcSimpleInteger(*this);
}

// Assignment operator 
CbcSimpleInteger & 
CbcSimpleInteger::operator=( const CbcSimpleInteger& rhs)
{
  if (this!=&rhs) {
    CbcObject::operator=(rhs);
    columnNumber_ = rhs.columnNumber_;
    originalLower_ = rhs.originalLower_;
    originalUpper_ = rhs.originalUpper_;
    breakEven_ = rhs.breakEven_;
    preferredWay_ = rhs.preferredWay_;
  }
  return *this;
}

// Destructor 
CbcSimpleInteger::~CbcSimpleInteger ()
{
}
// Construct an OsiSimpleInteger object
OsiSimpleInteger * 
CbcSimpleInteger::osiObject() const
{
  OsiSimpleInteger * obj = new OsiSimpleInteger(columnNumber_,
						originalLower_,originalUpper_);
  obj->setPriority(priority());
  return obj;
}

double
CbcSimpleInteger::infeasibility(const OsiSolverInterface * solver, const OsiBranchingInformation * info,
			 int & preferredWay) const
{
  double value = info->solution_[columnNumber_];
  value = CoinMax(value, info->lower_[columnNumber_]);
  value = CoinMin(value, info->upper_[columnNumber_]);
  double nearest = floor(value+(1.0-breakEven_));
  assert (breakEven_>0.0&&breakEven_<1.0);
  if (nearest>value) 
    preferredWay=1;
  else
    preferredWay=-1;
  if (preferredWay_)
    preferredWay=preferredWay_;
  double weight = fabs(value-nearest);
  // normalize so weight is 0.5 at break even
  if (nearest<value)
    weight = (0.5/breakEven_)*weight;
  else
    weight = (0.5/(1.0-breakEven_))*weight;
  if (fabs(value-nearest)<=info->integerTolerance_) 
    return 0.0;
  else
    return weight;
}
double 
CbcSimpleInteger::feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const 
{
  double value = info->solution_[columnNumber_];
  double newValue = CoinMax(value, info->lower_[columnNumber_]);
  newValue = CoinMin(newValue, info->upper_[columnNumber_]);
  newValue = floor(newValue+0.5);
  solver->setColLower(columnNumber_,newValue);
  solver->setColUpper(columnNumber_,newValue);
  return fabs(value-newValue);
}

/* Create an OsiSolverBranch object
   
This returns NULL if branch not represented by bound changes
*/
OsiSolverBranch * 
CbcSimpleInteger::solverBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info) const
{
  double value = info->solution_[columnNumber_];
  value = CoinMax(value, info->lower_[columnNumber_]);
  value = CoinMin(value, info->upper_[columnNumber_]);
  assert (info->upper_[columnNumber_]>info->lower_[columnNumber_]);
#ifndef NDEBUG
  double nearest = floor(value+0.5);
  assert (fabs(value-nearest)>info->integerTolerance_);
#endif
  OsiSolverBranch * branch = new OsiSolverBranch();
  branch->addBranch(columnNumber_,value);
  return branch;
}
// Creates a branching object
CbcBranchingObject * 
CbcSimpleInteger::createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) 
{
  double value = info->solution_[columnNumber_];
  value = CoinMax(value, info->lower_[columnNumber_]);
  value = CoinMin(value, info->upper_[columnNumber_]);
  assert (info->upper_[columnNumber_]>info->lower_[columnNumber_]);
  if (!info->hotstartSolution_) {
#ifndef NDEBUG
    double nearest = floor(value+0.5);
    assert (fabs(value-nearest)>info->integerTolerance_);
#endif
  } else {
    double targetValue = info->hotstartSolution_[columnNumber_];
    if (way>0)
      value = targetValue-0.1;
    else
      value = targetValue+0.1;
  }
  CbcBranchingObject * branch = new CbcIntegerBranchingObject(model_,columnNumber_,way,
					     value);
  branch->setOriginalObject(this);
  return branch;
}
/* Column number if single column object -1 otherwise,
   so returns >= 0
   Used by heuristics
*/
int 
CbcSimpleInteger::columnNumber() const
{
  return columnNumber_;
}
/* Reset variable bounds to their original values.
  
    Bounds may be tightened, so it may be good to be able to set this info in object.
*/
void 
CbcSimpleInteger::resetBounds(const OsiSolverInterface * solver) 
{
  originalLower_ = solver->getColLower()[columnNumber_] ;
  originalUpper_ = solver->getColUpper()[columnNumber_] ;
}

/*  Change column numbers after preprocessing
 */
void 
CbcSimpleInteger::resetSequenceEtc(int numberColumns, const int * originalColumns) 
{
  int iColumn;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (columnNumber_==originalColumns[iColumn])
      break;
  }
  assert (iColumn<numberColumns);
  columnNumber_ = iColumn;
}

// Infeasibility - large is 0.5
double 
CbcSimpleInteger::infeasibility(int & preferredWay) const
{
  OsiBranchingInformation info(model_->solver(),model_->normalSolver(),false);
  return infeasibility(model_->solver(),&info,preferredWay);
}

// This looks at solution and sets bounds to contain solution
/** More precisely: it first forces the variable within the existing
    bounds, and then tightens the bounds to fix the variable at the
    nearest integer value.
*/
void 
CbcSimpleInteger::feasibleRegion()
{
  abort();
}
CbcBranchingObject * 
CbcSimpleInteger::createBranch( int way) 
{
  abort();
  return NULL;
}
// Default Constructor 
CbcIntegerBranchingObject::CbcIntegerBranchingObject()
  :CbcBranchingObject()
{
  down_[0] = 0.0;
  down_[1] = 0.0;
  up_[0] = 0.0;
  up_[1] = 0.0;
}

// Useful constructor
CbcIntegerBranchingObject::CbcIntegerBranchingObject (CbcModel * model, 
						      int variable, int way , double value)
  :CbcBranchingObject(model,variable,way,value)
{
  int iColumn = variable;
  down_[0] = model_->solver()->getColLower()[iColumn];
  down_[1] = floor(value_);
  up_[0] = ceil(value_);
  up_[1] = model->getColUpper()[iColumn];
}
// Useful constructor for fixing
CbcIntegerBranchingObject::CbcIntegerBranchingObject (CbcModel * model, 
						      int variable, int way,
						      double lowerValue, 
						      double upperValue)
  :CbcBranchingObject(model,variable,way,lowerValue)
{
  setNumberBranchesLeft(1);
  down_[0] = lowerValue;
  down_[1] = upperValue;
  up_[0] = lowerValue;
  up_[1] = upperValue;
}
  

// Copy constructor 
CbcIntegerBranchingObject::CbcIntegerBranchingObject ( const CbcIntegerBranchingObject & rhs) :CbcBranchingObject(rhs)
{
  down_[0] = rhs.down_[0];
  down_[1] = rhs.down_[1];
  up_[0] = rhs.up_[0];
  up_[1] = rhs.up_[1];
}

// Assignment operator 
CbcIntegerBranchingObject & 
CbcIntegerBranchingObject::operator=( const CbcIntegerBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    down_[0] = rhs.down_[0];
    down_[1] = rhs.down_[1];
    up_[0] = rhs.up_[0];
    up_[1] = rhs.up_[1];
  }
  return *this;
}
CbcBranchingObject * 
CbcIntegerBranchingObject::clone() const
{ 
  return (new CbcIntegerBranchingObject(*this));
}


// Destructor 
CbcIntegerBranchingObject::~CbcIntegerBranchingObject ()
{
  // for debugging threads
  way_=-23456789;
}

/*
  Perform a branch by adjusting the bounds of the specified variable. Note
  that each arm of the branch advances the object to the next arm by
  advancing the value of way_.

  Providing new values for the variable's lower and upper bounds for each
  branching direction gives a little bit of additional flexibility and will
  be easily extensible to multi-way branching.
  Returns change in guessed objective on next branch
*/
double
CbcIntegerBranchingObject::branch()
{
  // for debugging threads
  if (way_<-1||way_>100000) {
    printf("way %d, left %d, iCol %d, variable %d\n",
	   way_,numberBranchesLeft(),
	   originalCbcObject_->columnNumber(),variable_);
    assert (way_!=-23456789);
  }
  decrementNumberBranchesLeft();
  int iColumn = originalCbcObject_->columnNumber();
  assert (variable_==iColumn);
  double olb,oub ;
  olb = model_->solver()->getColLower()[iColumn] ;
  oub = model_->solver()->getColUpper()[iColumn] ;
#ifdef COIN_DEVELOP
  if (olb!=down_[0]||oub!=up_[1]) {
    if (way_>0)
      printf("branching up on var %d: [%g,%g] => [%g,%g] - other [%g,%g]\n",
	     iColumn,olb,oub,up_[0],up_[1],down_[0],down_[1]) ; 
    else
      printf("branching down on var %d: [%g,%g] => [%g,%g] - other [%g,%g]\n",
	     iColumn,olb,oub,down_[0],down_[1],up_[0],up_[1]) ; 
  }
#endif
  if (way_<0) {
#ifdef CBC_DEBUG
  { double olb,oub ;
    olb = model_->solver()->getColLower()[iColumn] ;
    oub = model_->solver()->getColUpper()[iColumn] ;
    printf("branching down on var %d: [%g,%g] => [%g,%g]\n",
	   iColumn,olb,oub,down_[0],down_[1]) ; }
#endif
    model_->solver()->setColLower(iColumn,down_[0]);
    model_->solver()->setColUpper(iColumn,down_[1]);
    way_=1;
  } else {
#ifdef CBC_DEBUG
  { double olb,oub ;
    olb = model_->solver()->getColLower()[iColumn] ;
    oub = model_->solver()->getColUpper()[iColumn] ;
    printf("branching up on var %d: [%g,%g] => [%g,%g]\n",
	   iColumn,olb,oub,up_[0],up_[1]) ; }
#endif
    model_->solver()->setColLower(iColumn,up_[0]);
    model_->solver()->setColUpper(iColumn,up_[1]);
    way_=-1;	  // Swap direction
  }
  double nlb = model_->solver()->getColLower()[iColumn];
  double nub = model_->solver()->getColUpper()[iColumn];
  if (nlb<olb) {
#ifndef NDEBUG
    printf("bad lb change for column %d from %g to %g\n",iColumn,olb,nlb);
#endif
    model_->solver()->setColLower(iColumn,CoinMin(olb,nub));
    nlb=olb;
  }
  if (nub>oub) {
#ifndef NDEBUG
    printf("bad ub change for column %d from %g to %g\n",iColumn,oub,nub);
#endif
    model_->solver()->setColUpper(iColumn,CoinMax(oub,nlb));
  }
#ifndef NDEBUG
  if (nlb<olb+1.0e-8&&nub>oub-1.0e-8)
    printf("bad null change for column %d - bounds %g,%g\n",iColumn,olb,oub);
#endif
  return 0.0;
}
// Print what would happen  
void
CbcIntegerBranchingObject::print()
{
  int iColumn = originalCbcObject_->columnNumber();
  assert (variable_==iColumn);
  if (way_<0) {
  { double olb,oub ;
    olb = model_->solver()->getColLower()[iColumn] ;
    oub = model_->solver()->getColUpper()[iColumn] ;
    printf("CbcInteger would branch down on var %d (int var %d): [%g,%g] => [%g,%g]\n",
	   iColumn,variable_,olb,oub,down_[0],down_[1]) ; }
  } else {
  { double olb,oub ;
    olb = model_->solver()->getColLower()[iColumn] ;
    oub = model_->solver()->getColUpper()[iColumn] ;
    printf("CbcInteger would branch up on var %d (int var %d): [%g,%g] => [%g,%g]\n",
	   iColumn,variable_,olb,oub,up_[0],up_[1]) ; }
  }
}


/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
CbcSimpleIntegerPseudoCost::CbcSimpleIntegerPseudoCost ()
  : CbcSimpleInteger(),
    downPseudoCost_(1.0e-5),
    upPseudoCost_(1.0e-5),
    upDownSeparator_(-1.0),
    method_(0)
{
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
CbcSimpleIntegerPseudoCost::CbcSimpleIntegerPseudoCost (CbcModel * model,
				    int iColumn, double breakEven)
  : CbcSimpleInteger(model,iColumn,breakEven)
{
  const double * cost = model->getObjCoefficients();
  double costValue = CoinMax(1.0e-5,fabs(cost[iColumn]));
  // treat as if will cost what it says up
  upPseudoCost_=costValue;
  // and balance at breakeven
  downPseudoCost_=((1.0-breakEven_)*upPseudoCost_)/breakEven_;
  upDownSeparator_ = -1.0;
  method_=0;
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
CbcSimpleIntegerPseudoCost::CbcSimpleIntegerPseudoCost (CbcModel * model,
				    int iColumn, double downPseudoCost,
							double upPseudoCost)
  : CbcSimpleInteger(model,iColumn)
{
  downPseudoCost_ = CoinMax(1.0e-10,downPseudoCost);
  upPseudoCost_ = CoinMax(1.0e-10,upPseudoCost);
  breakEven_ = upPseudoCost_/(upPseudoCost_+downPseudoCost_);
  upDownSeparator_ = -1.0;
  method_=0;
}
// Useful constructor - passed and model index and pseudo costs
CbcSimpleIntegerPseudoCost::CbcSimpleIntegerPseudoCost (CbcModel * model, int dummy,int iColumn, 
							double downPseudoCost, double upPseudoCost)
{
  CbcSimpleIntegerPseudoCost(model,iColumn,downPseudoCost,upPseudoCost);
  columnNumber_=iColumn;
}

// Copy constructor 
CbcSimpleIntegerPseudoCost::CbcSimpleIntegerPseudoCost ( const CbcSimpleIntegerPseudoCost & rhs)
  :CbcSimpleInteger(rhs),
   downPseudoCost_(rhs.downPseudoCost_),
   upPseudoCost_(rhs.upPseudoCost_),
   upDownSeparator_(rhs.upDownSeparator_),
   method_(rhs.method_)

{
}

// Clone
CbcObject *
CbcSimpleIntegerPseudoCost::clone() const
{
  return new CbcSimpleIntegerPseudoCost(*this);
}

// Assignment operator 
CbcSimpleIntegerPseudoCost & 
CbcSimpleIntegerPseudoCost::operator=( const CbcSimpleIntegerPseudoCost& rhs)
{
  if (this!=&rhs) {
    CbcSimpleInteger::operator=(rhs);
    downPseudoCost_=rhs.downPseudoCost_;
    upPseudoCost_=rhs.upPseudoCost_;
    upDownSeparator_=rhs.upDownSeparator_;
    method_=rhs.method_;
  }
  return *this;
}

// Destructor 
CbcSimpleIntegerPseudoCost::~CbcSimpleIntegerPseudoCost ()
{
}
// Creates a branching object
CbcBranchingObject * 
CbcSimpleIntegerPseudoCost::createBranch(int way) 
{
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
#ifndef NDEBUG
  double nearest = floor(value+0.5);
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  assert (upper[columnNumber_]>lower[columnNumber_]);
#endif
  if (!model_->hotstartSolution()) {
    assert (fabs(value-nearest)>integerTolerance);
  } else {
    const double * hotstartSolution = model_->hotstartSolution();
    double targetValue = hotstartSolution[columnNumber_];
    if (way>0)
      value = targetValue-0.1;
    else
      value = targetValue+0.1;
  }
  CbcIntegerPseudoCostBranchingObject * newObject = 
    new CbcIntegerPseudoCostBranchingObject(model_,columnNumber_,way,
					    value);
  double up =  upPseudoCost_*(ceil(value)-value);
  double down =  downPseudoCost_*(value-floor(value));
  double changeInGuessed=up-down;
  if (way>0)
    changeInGuessed = - changeInGuessed;
  changeInGuessed=CoinMax(0.0,changeInGuessed);
  //if (way>0)
  //changeInGuessed += 1.0e8; // bias to stay up
  newObject->setChangeInGuessed(changeInGuessed);
  newObject->setOriginalObject(this);
  return newObject;
}
// Infeasibility - large is 0.5
double 
CbcSimpleIntegerPseudoCost::infeasibility(int & preferredWay) const
{
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  if (upper[columnNumber_]==lower[columnNumber_]) {
    // fixed
    preferredWay=1;
    return 0.0;
  }
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  /*printf("%d %g %g %g %g\n",columnNumber_,value,lower[columnNumber_],
    solution[columnNumber_],upper[columnNumber_]);*/
  double nearest = floor(value+0.5);
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double below = floor(value+integerTolerance);
  double above = below+1.0;
  if (above>upper[columnNumber_]) {
    above=below;
    below = above -1;
  }
  double downCost = CoinMax((value-below)*downPseudoCost_,0.0);
  double upCost = CoinMax((above-value)*upPseudoCost_,0.0);
  if (downCost>=upCost)
    preferredWay=1;
  else
    preferredWay=-1;
  // See if up down choice set
  if (upDownSeparator_>0.0) {
    preferredWay = (value-below>=upDownSeparator_) ? 1 : -1;
  }
  if (preferredWay_)
    preferredWay=preferredWay_;
  if (fabs(value-nearest)<=integerTolerance) {
    return 0.0;
  } else {
    // can't get at model so 1,2 don't make sense
    assert(method_<1||method_>2);
    if (!method_)
      return CoinMin(downCost,upCost);
    else
      return CoinMax(downCost,upCost);
  }
}

// Return "up" estimate
double 
CbcSimpleIntegerPseudoCost::upEstimate() const
{
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  if (upper[columnNumber_]==lower[columnNumber_]) {
    // fixed
    return 0.0;
  }
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double below = floor(value+integerTolerance);
  double above = below+1.0;
  if (above>upper[columnNumber_]) {
    above=below;
    below = above -1;
  }
  double upCost = CoinMax((above-value)*upPseudoCost_,0.0);
  return upCost;
}
// Return "down" estimate
double 
CbcSimpleIntegerPseudoCost::downEstimate() const
{
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  if (upper[columnNumber_]==lower[columnNumber_]) {
    // fixed
    return 0.0;
  }
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double below = floor(value+integerTolerance);
  double above = below+1.0;
  if (above>upper[columnNumber_]) {
    above=below;
    below = above -1;
  }
  double downCost = CoinMax((value-below)*downPseudoCost_,0.0);
  return downCost;
}

// Default Constructor 
CbcIntegerPseudoCostBranchingObject::CbcIntegerPseudoCostBranchingObject()
  :CbcIntegerBranchingObject()
{
  changeInGuessed_=1.0e-5;
}

// Useful constructor
CbcIntegerPseudoCostBranchingObject::CbcIntegerPseudoCostBranchingObject (CbcModel * model, 
						      int variable, int way , double value)
  :CbcIntegerBranchingObject(model,variable,way,value)
{
}
// Useful constructor for fixing
CbcIntegerPseudoCostBranchingObject::CbcIntegerPseudoCostBranchingObject (CbcModel * model, 
						      int variable, int way,
						      double lowerValue, 
						      double upperValue)
  :CbcIntegerBranchingObject(model,variable,way,lowerValue)
{
  changeInGuessed_=1.0e100;
}
  

// Copy constructor 
CbcIntegerPseudoCostBranchingObject::CbcIntegerPseudoCostBranchingObject ( 
				 const CbcIntegerPseudoCostBranchingObject & rhs)
  :CbcIntegerBranchingObject(rhs)
{
  changeInGuessed_ = rhs.changeInGuessed_;
}

// Assignment operator 
CbcIntegerPseudoCostBranchingObject & 
CbcIntegerPseudoCostBranchingObject::operator=( const CbcIntegerPseudoCostBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcIntegerBranchingObject::operator=(rhs);
    changeInGuessed_ = rhs.changeInGuessed_;
  }
  return *this;
}
CbcBranchingObject * 
CbcIntegerPseudoCostBranchingObject::clone() const
{ 
  return (new CbcIntegerPseudoCostBranchingObject(*this));
}


// Destructor 
CbcIntegerPseudoCostBranchingObject::~CbcIntegerPseudoCostBranchingObject ()
{
}

/*
  Perform a branch by adjusting the bounds of the specified variable. Note
  that each arm of the branch advances the object to the next arm by
  advancing the value of way_.

  Providing new values for the variable's lower and upper bounds for each
  branching direction gives a little bit of additional flexibility and will
  be easily extensible to multi-way branching.
  Returns change in guessed objective on next branch
*/
double
CbcIntegerPseudoCostBranchingObject::branch()
{
  CbcIntegerBranchingObject::branch();
  return changeInGuessed_;
}


// Default Constructor 
CbcCliqueBranchingObject::CbcCliqueBranchingObject()
  :CbcBranchingObject()
{
  clique_ = NULL;
  downMask_[0]=0;
  downMask_[1]=0;
  upMask_[0]=0;
  upMask_[1]=0;
}

// Useful constructor
CbcCliqueBranchingObject::CbcCliqueBranchingObject (CbcModel * model,
						    const CbcClique * clique,
						    int way ,
						    int numberOnDownSide, const int * down,
						    int numberOnUpSide, const int * up)
  :CbcBranchingObject(model,clique->id(),way,0.5)
{
  clique_ = clique;
  downMask_[0]=0;
  downMask_[1]=0;
  upMask_[0]=0;
  upMask_[1]=0;
  int i;
  for (i=0;i<numberOnDownSide;i++) {
    int sequence = down[i];
    int iWord = sequence>>5;
    int iBit = sequence - 32*iWord;
    unsigned int k = 1<<iBit;
    downMask_[iWord] |= k;
  }
  for (i=0;i<numberOnUpSide;i++) {
    int sequence = up[i];
    int iWord = sequence>>5;
    int iBit = sequence - 32*iWord;
    unsigned int k = 1<<iBit;
    upMask_[iWord] |= k;
  }
}

// Copy constructor 
CbcCliqueBranchingObject::CbcCliqueBranchingObject ( const CbcCliqueBranchingObject & rhs) :CbcBranchingObject(rhs)
{
  clique_=rhs.clique_;
  downMask_[0]=rhs.downMask_[0];
  downMask_[1]=rhs.downMask_[1];
  upMask_[0]=rhs.upMask_[0];
  upMask_[1]=rhs.upMask_[1];
}

// Assignment operator 
CbcCliqueBranchingObject & 
CbcCliqueBranchingObject::operator=( const CbcCliqueBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    clique_=rhs.clique_;
    downMask_[0]=rhs.downMask_[0];
    downMask_[1]=rhs.downMask_[1];
    upMask_[0]=rhs.upMask_[0];
    upMask_[1]=rhs.upMask_[1];
  }
  return *this;
}
CbcBranchingObject * 
CbcCliqueBranchingObject::clone() const
{ 
  return (new CbcCliqueBranchingObject(*this));
}


// Destructor 
CbcCliqueBranchingObject::~CbcCliqueBranchingObject ()
{
}
double
CbcCliqueBranchingObject::branch()
{
  decrementNumberBranchesLeft();
  int iWord;
  int numberMembers = clique_->numberMembers();
  const int * which = clique_->members();
  const int * integerVariables = model_->integerVariable();
  int numberWords=(numberMembers+31)>>5;
  // *** for way - up means fix all those in down section
  if (way_<0) {
#ifdef FULL_PRINT
    printf("Down Fix ");
#endif
    for (iWord=0;iWord<numberWords;iWord++) {
      int i;
      for (i=0;i<32;i++) {
	unsigned int k = 1<<i;
	if ((upMask_[iWord]&k)!=0) {
	  int iColumn = which[i+32*iWord];
#ifdef FULL_PRINT
	  printf("%d ",i+32*iWord);
#endif
	  // fix weak way
	  if (clique_->type(i+32*iWord))
	    model_->solver()->setColUpper(integerVariables[iColumn],0.0);
	  else
	    model_->solver()->setColLower(integerVariables[iColumn],1.0);
	}
      }
    }
    way_=1;	  // Swap direction
  } else {
#ifdef FULL_PRINT
    printf("Up Fix ");
#endif
    for (iWord=0;iWord<numberWords;iWord++) {
      int i;
      for (i=0;i<32;i++) {
	unsigned int k = 1<<i;
	if ((downMask_[iWord]&k)!=0) {
	  int iColumn = which[i+32*iWord];
#ifdef FULL_PRINT
	  printf("%d ",i+32*iWord);
#endif
	  // fix weak way
	  if (clique_->type(i+32*iWord))
	    model_->solver()->setColUpper(integerVariables[iColumn],0.0);
	  else
	    model_->solver()->setColLower(integerVariables[iColumn],1.0);
	}
      }
    }
    way_=-1;	  // Swap direction
  }
#ifdef FULL_PRINT
  printf("\n");
#endif
  return 0.0;
}
// Print what would happen  
void
CbcCliqueBranchingObject::print()
{
  int iWord;
  int numberMembers = clique_->numberMembers();
  const int * which = clique_->members();
  const int * integerVariables = model_->integerVariable();
  int numberWords=(numberMembers+31)>>5;
  // *** for way - up means fix all those in down section
  if (way_<0) {
    printf("Clique - Down Fix ");
    for (iWord=0;iWord<numberWords;iWord++) {
      int i;
      for (i=0;i<32;i++) {
	unsigned int k = 1<<i;
	if ((upMask_[iWord]&k)!=0) {
	  int iColumn = which[i+32*iWord];
	  printf("%d ",integerVariables[iColumn]);
	}
      }
    }
  } else {
    printf("Clique - Up Fix ");
    for (iWord=0;iWord<numberWords;iWord++) {
      int i;
      for (i=0;i<32;i++) {
	unsigned int k = 1<<i;
	if ((downMask_[iWord]&k)!=0) {
	  int iColumn = which[i+32*iWord];
	  printf("%d ",integerVariables[iColumn]);
	}
      }
    }
  }
  printf("\n");
}
  
// Default Constructor 
CbcLongCliqueBranchingObject::CbcLongCliqueBranchingObject()
  :CbcBranchingObject()
{
  clique_=NULL;
  downMask_=NULL;
  upMask_=NULL;
}

// Useful constructor
CbcLongCliqueBranchingObject::CbcLongCliqueBranchingObject (CbcModel * model,
							    const CbcClique * clique, 
							    int way ,
						    int numberOnDownSide, const int * down,
						    int numberOnUpSide, const int * up)
  :CbcBranchingObject(model,clique->id(),way,0.5)
{
  clique_ = clique;
  int numberMembers = clique_->numberMembers();
  int numberWords=(numberMembers+31)>>5;
  downMask_ = new unsigned int [numberWords];
  upMask_ = new unsigned int [numberWords];
  memset(downMask_,0,numberWords*sizeof(unsigned int));
  memset(upMask_,0,numberWords*sizeof(unsigned int));
  int i;
  for (i=0;i<numberOnDownSide;i++) {
    int sequence = down[i];
    int iWord = sequence>>5;
    int iBit = sequence - 32*iWord;
    unsigned int k = 1<<iBit;
    downMask_[iWord] |= k;
  }
  for (i=0;i<numberOnUpSide;i++) {
    int sequence = up[i];
    int iWord = sequence>>5;
    int iBit = sequence - 32*iWord;
    unsigned int k = 1<<iBit;
    upMask_[iWord] |= k;
  }
}

// Copy constructor 
CbcLongCliqueBranchingObject::CbcLongCliqueBranchingObject ( const CbcLongCliqueBranchingObject & rhs) :CbcBranchingObject(rhs)
{
  clique_=rhs.clique_;
  if (rhs.downMask_) {
    int numberMembers = clique_->numberMembers();
    int numberWords=(numberMembers+31)>>5;
    downMask_ = new unsigned int [numberWords];
    memcpy(downMask_,rhs.downMask_,numberWords*sizeof(unsigned int));
    upMask_ = new unsigned int [numberWords];
    memcpy(upMask_,rhs.upMask_,numberWords*sizeof(unsigned int));
  } else {
    downMask_=NULL;
    upMask_=NULL;
  }    
}

// Assignment operator 
CbcLongCliqueBranchingObject & 
CbcLongCliqueBranchingObject::operator=( const CbcLongCliqueBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    clique_=rhs.clique_;
    delete [] downMask_;
    delete [] upMask_;
    if (rhs.downMask_) {
      int numberMembers = clique_->numberMembers();
      int numberWords=(numberMembers+31)>>5;
      downMask_ = new unsigned int [numberWords];
      memcpy(downMask_,rhs.downMask_,numberWords*sizeof(unsigned int));
      upMask_ = new unsigned int [numberWords];
      memcpy(upMask_,rhs.upMask_,numberWords*sizeof(unsigned int));
    } else {
      downMask_=NULL;
      upMask_=NULL;
    }    
  }
  return *this;
}
CbcBranchingObject * 
CbcLongCliqueBranchingObject::clone() const
{ 
  return (new CbcLongCliqueBranchingObject(*this));
}


// Destructor 
CbcLongCliqueBranchingObject::~CbcLongCliqueBranchingObject ()
{
  delete [] downMask_;
  delete [] upMask_;
}
double
CbcLongCliqueBranchingObject::branch()
{
  decrementNumberBranchesLeft();
  int iWord;
  int numberMembers = clique_->numberMembers();
  const int * which = clique_->members();
  const int * integerVariables = model_->integerVariable();
  int numberWords=(numberMembers+31)>>5;
  // *** for way - up means fix all those in down section
  if (way_<0) {
#ifdef FULL_PRINT
    printf("Down Fix ");
#endif
    for (iWord=0;iWord<numberWords;iWord++) {
      int i;
      for (i=0;i<32;i++) {
	unsigned int k = 1<<i;
	if ((upMask_[iWord]&k)!=0) {
	  int iColumn = which[i+32*iWord];
#ifdef FULL_PRINT
	  printf("%d ",i+32*iWord);
#endif
	  // fix weak way
	  if (clique_->type(i+32*iWord))
	    model_->solver()->setColUpper(integerVariables[iColumn],0.0);
	  else
	    model_->solver()->setColLower(integerVariables[iColumn],1.0);
	}
      }
    }
    way_=1;	  // Swap direction
  } else {
#ifdef FULL_PRINT
    printf("Up Fix ");
#endif
    for (iWord=0;iWord<numberWords;iWord++) {
      int i;
      for (i=0;i<32;i++) {
	unsigned int k = 1<<i;
	if ((downMask_[iWord]&k)!=0) {
	  int iColumn = which[i+32*iWord];
#ifdef FULL_PRINT
	  printf("%d ",i+32*iWord);
#endif
	  // fix weak way
	  if (clique_->type(i+32*iWord))
	    model_->solver()->setColUpper(integerVariables[iColumn],0.0);
	  else
	    model_->solver()->setColLower(integerVariables[iColumn],1.0);
	}
      }
    }
    way_=-1;	  // Swap direction
  }
#ifdef FULL_PRINT
  printf("\n");
#endif
  return 0.0;
}
void
CbcLongCliqueBranchingObject::print()
{
  int iWord;
  int numberMembers = clique_->numberMembers();
  const int * which = clique_->members();
  const int * integerVariables = model_->integerVariable();
  int numberWords=(numberMembers+31)>>5;
  // *** for way - up means fix all those in down section
  if (way_<0) {
    printf("Clique - Down Fix ");
    for (iWord=0;iWord<numberWords;iWord++) {
      int i;
      for (i=0;i<32;i++) {
	unsigned int k = 1<<i;
	if ((upMask_[iWord]&k)!=0) {
	  int iColumn = which[i+32*iWord];
	  printf("%d ",integerVariables[iColumn]);
	}
      }
    }
  } else {
    printf("Clique - Up Fix ");
    for (iWord=0;iWord<numberWords;iWord++) {
      int i;
      for (i=0;i<32;i++) {
	unsigned int k = 1<<i;
	if ((downMask_[iWord]&k)!=0) {
	  int iColumn = which[i+32*iWord];
	  printf("%d ",integerVariables[iColumn]);
	}
      }
    }
  }
  printf("\n");
}
// Default Constructor 
CbcSOSBranchingObject::CbcSOSBranchingObject()
  :CbcBranchingObject()
{
  set_ = NULL;
  separator_=0.0;
}

// Useful constructor
CbcSOSBranchingObject::CbcSOSBranchingObject (CbcModel * model,
					      const CbcSOS * set,
					      int way ,
					      double separator)
  :CbcBranchingObject(model,set->id(),way,0.5)
{
  set_ = set;
  separator_ = separator;
}

// Copy constructor 
CbcSOSBranchingObject::CbcSOSBranchingObject ( const CbcSOSBranchingObject & rhs) :CbcBranchingObject(rhs)
{
  set_=rhs.set_;
  separator_ = rhs.separator_;
}

// Assignment operator 
CbcSOSBranchingObject & 
CbcSOSBranchingObject::operator=( const CbcSOSBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    set_=rhs.set_;
    separator_ = rhs.separator_;
  }
  return *this;
}
CbcBranchingObject * 
CbcSOSBranchingObject::clone() const
{ 
  return (new CbcSOSBranchingObject(*this));
}


// Destructor 
CbcSOSBranchingObject::~CbcSOSBranchingObject ()
{
}
double
CbcSOSBranchingObject::branch()
{
  decrementNumberBranchesLeft();
  int numberMembers = set_->numberMembers();
  const int * which = set_->members();
  const double * weights = set_->weights();
  OsiSolverInterface * solver = model_->solver();
  //const double * lower = solver->getColLower();
  //const double * upper = solver->getColUpper();
  // *** for way - up means fix all those in down section
  if (way_<0) {
    int i;
    for ( i=0;i<numberMembers;i++) {
      if (weights[i] > separator_)
	break;
    }
    assert (i<numberMembers);
    for (;i<numberMembers;i++) 
      solver->setColUpper(which[i],0.0);
    way_=1;	  // Swap direction
  } else {
    int i;
    for ( i=0;i<numberMembers;i++) {
      if (weights[i] >= separator_)
	break;
      else
	solver->setColUpper(which[i],0.0);
    }
    assert (i<numberMembers);
    way_=-1;	  // Swap direction
  }
  return 0.0;
}
// Print what would happen  
void
CbcSOSBranchingObject::print()
{
  int numberMembers = set_->numberMembers();
  const int * which = set_->members();
  const double * weights = set_->weights();
  OsiSolverInterface * solver = model_->solver();
  //const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  int first=numberMembers;
  int last=-1;
  int numberFixed=0;
  int numberOther=0;
  int i;
  for ( i=0;i<numberMembers;i++) {
    double bound = upper[which[i]];
    if (bound) {
      first = CoinMin(first,i);
      last = CoinMax(last,i);
    }
  }
  // *** for way - up means fix all those in down section
  if (way_<0) {
    printf("SOS Down");
    for ( i=0;i<numberMembers;i++) {
      double bound = upper[which[i]];
      if (weights[i] > separator_)
	break;
      else if (bound)
	numberOther++;
    }
    assert (i<numberMembers);
    for (;i<numberMembers;i++) {
      double bound = upper[which[i]];
      if (bound)
	numberFixed++;
    }
  } else {
    printf("SOS Up");
    for ( i=0;i<numberMembers;i++) {
      double bound = upper[which[i]];
      if (weights[i] >= separator_)
	break;
      else if (bound)
	numberFixed++;
    }
    assert (i<numberMembers);
    for (;i<numberMembers;i++) {
      double bound = upper[which[i]];
      if (bound)
	numberOther++;
    }
  }
  printf(" - at %g, free range %d (%g) => %d (%g), %d would be fixed, %d other way\n",
	 separator_,which[first],weights[first],which[last],weights[last],numberFixed,numberOther);
}
  
// Default Constructor 
CbcBranchDefaultDecision::CbcBranchDefaultDecision()
  :CbcBranchDecision()
{
  bestCriterion_ = 0.0;
  bestChangeUp_ = 0.0;
  bestNumberUp_ = 0;
  bestChangeDown_ = 0.0;
  bestObject_ = NULL;
  bestNumberDown_ = 0;
}

// Copy constructor 
CbcBranchDefaultDecision::CbcBranchDefaultDecision (
				    const CbcBranchDefaultDecision & rhs)
  :CbcBranchDecision(rhs)
{
  bestCriterion_ = rhs.bestCriterion_;
  bestChangeUp_ = rhs.bestChangeUp_;
  bestNumberUp_ = rhs.bestNumberUp_;
  bestChangeDown_ = rhs.bestChangeDown_;
  bestNumberDown_ = rhs.bestNumberDown_;
  bestObject_ = rhs.bestObject_;
  model_ = rhs.model_;
}

CbcBranchDefaultDecision::~CbcBranchDefaultDecision()
{
}

// Clone
CbcBranchDecision * 
CbcBranchDefaultDecision::clone() const
{
  return new CbcBranchDefaultDecision(*this);
}

// Initialize i.e. before start of choosing at a node
void 
CbcBranchDefaultDecision::initialize(CbcModel * model)
{
  bestCriterion_ = 0.0;
  bestChangeUp_ = 0.0;
  bestNumberUp_ = 0;
  bestChangeDown_ = 0.0;
  bestNumberDown_ = 0;
  bestObject_ = NULL;
  model_ = model;
}


/*
  Simple default decision algorithm. Compare based on infeasibility (numInfUp,
  numInfDn) until a solution is found by search, then switch to change in
  objective (changeUp, changeDn). Note that bestSoFar is remembered in
  bestObject_, so the parameter bestSoFar is unused.
*/

int
CbcBranchDefaultDecision::betterBranch(CbcBranchingObject * thisOne,
			    CbcBranchingObject * bestSoFar,
			    double changeUp, int numInfUp,
			    double changeDn, int numInfDn)
{
  bool beforeSolution = cbcModel()->getSolutionCount()==
    cbcModel()->getNumberHeuristicSolutions();;
  int betterWay=0;
  if (beforeSolution) {
    if (!bestObject_) {
      bestNumberUp_=COIN_INT_MAX;
      bestNumberDown_=COIN_INT_MAX;
    }
    // before solution - choose smallest number 
    // could add in depth as well
    int bestNumber = CoinMin(bestNumberUp_,bestNumberDown_);
    if (numInfUp<numInfDn) {
      if (numInfUp<bestNumber) {
	betterWay = 1;
      } else if (numInfUp==bestNumber) {
	if (changeUp<bestCriterion_)
	  betterWay=1;
      }
    } else if (numInfUp>numInfDn) {
      if (numInfDn<bestNumber) {
	betterWay = -1;
      } else if (numInfDn==bestNumber) {
	if (changeDn<bestCriterion_)
	  betterWay=-1;
      }
    } else {
      // up and down have same number
      bool better=false;
      if (numInfUp<bestNumber) {
	better=true;
      } else if (numInfUp==bestNumber) {
	if (min(changeUp,changeDn)<bestCriterion_)
	  better=true;;
      }
      if (better) {
	// see which way
	if (changeUp<=changeDn)
	  betterWay=1;
	else
	  betterWay=-1;
      }
    }
  } else {
    if (!bestObject_) {
      bestCriterion_=-1.0;
    }
    // got a solution
    if (changeUp<=changeDn) {
      if (changeUp>bestCriterion_)
	betterWay=1;
    } else {
      if (changeDn>bestCriterion_)
	betterWay=-1;
    }
  }
  if (betterWay) {
    bestCriterion_ = CoinMin(changeUp,changeDn);
    bestChangeUp_ = changeUp;
    bestNumberUp_ = numInfUp;
    bestChangeDown_ = changeDn;
    bestNumberDown_ = numInfDn;
    bestObject_=thisOne;
    // See if user is overriding way
    if (thisOne->object()&&thisOne->object()->preferredWay())
      betterWay = thisOne->object()->preferredWay();
  }
  return betterWay;
}
/* Sets or gets best criterion so far */
void 
CbcBranchDefaultDecision::setBestCriterion(double value)
{ 
  bestCriterion_ = value;
}
double 
CbcBranchDefaultDecision::getBestCriterion() const
{ 
  return bestCriterion_;
}

/* Compare N branching objects. Return index of best
   and sets way of branching in chosen object.
   
   This routine is used only after strong branching.
*/

int
CbcBranchDefaultDecision::bestBranch (CbcBranchingObject ** objects, int numberObjects,
				   int numberUnsatisfied,
				   double * changeUp, int * numberInfeasibilitiesUp,
				   double * changeDown, int * numberInfeasibilitiesDown,
				   double objectiveValue) 
{

  int bestWay=0;
  int whichObject = -1;
  if (numberObjects) {
    CbcModel * model = cbcModel();
    // at continuous
    //double continuousObjective = model->getContinuousObjective();
    //int continuousInfeasibilities = model->getContinuousInfeasibilities();
    
    // average cost to get rid of infeasibility
    //double averageCostPerInfeasibility = 
    //(objectiveValue-continuousObjective)/
    //(double) (abs(continuousInfeasibilities-numberUnsatisfied)+1);
    /* beforeSolution is :
       0 - before any solution
       n - n heuristic solutions but no branched one
       -1 - branched solution found
    */
    int numberSolutions = model->getSolutionCount();
    double cutoff = model->getCutoff();
    int method=0;
    int i;
    if (numberSolutions) {
      int numberHeuristic = model->getNumberHeuristicSolutions();
      if (numberHeuristic<numberSolutions) {
	method = 1;
      } else {
	method = 2;
	// look further
	for ( i = 0 ; i < numberObjects ; i++) {
	  int numberNext = numberInfeasibilitiesUp[i];
	  
	  if (numberNext<numberUnsatisfied) {
	    int numberUp = numberUnsatisfied - numberInfeasibilitiesUp[i];
	    double perUnsatisfied = changeUp[i]/(double) numberUp;
	    double estimatedObjective = objectiveValue + numberUnsatisfied * perUnsatisfied;
	    if (estimatedObjective<cutoff) 
	      method=3;
	  }
	  numberNext = numberInfeasibilitiesDown[i];
	  if (numberNext<numberUnsatisfied) {
	    int numberDown = numberUnsatisfied - numberInfeasibilitiesDown[i];
	    double perUnsatisfied = changeDown[i]/(double) numberDown;
	    double estimatedObjective = objectiveValue + numberUnsatisfied * perUnsatisfied;
	    if (estimatedObjective<cutoff) 
	      method=3;
	  }
	}
      }
      method=2;
    } else {
      method = 0;
    }
    // Uncomment next to force method 4
    //method=4;
    /* Methods :
       0 - fewest infeasibilities
       1 - largest min change in objective
       2 - as 1 but use sum of changes if min close
       3 - predicted best solution
       4 - take cheapest up branch if infeasibilities same
    */
    int bestNumber=COIN_INT_MAX;
    double bestCriterion=-1.0e50;
    double alternativeCriterion = -1.0;
    double bestEstimate = 1.0e100;
    switch (method) {
    case 0:
      // could add in depth as well
      for ( i = 0 ; i < numberObjects ; i++) {
	int thisNumber = min(numberInfeasibilitiesUp[i],numberInfeasibilitiesDown[i]);
	if (thisNumber<=bestNumber) {
	  int betterWay=0;
	  if (numberInfeasibilitiesUp[i]<numberInfeasibilitiesDown[i]) {
	    if (numberInfeasibilitiesUp[i]<bestNumber) {
	      betterWay = 1;
	    } else {
	      if (changeUp[i]<bestCriterion)
		betterWay=1;
	    }
	  } else if (numberInfeasibilitiesUp[i]>numberInfeasibilitiesDown[i]) {
	    if (numberInfeasibilitiesDown[i]<bestNumber) {
	      betterWay = -1;
	    } else {
	      if (changeDown[i]<bestCriterion)
		betterWay=-1;
	    }
	  } else {
	    // up and down have same number
	    bool better=false;
	    if (numberInfeasibilitiesUp[i]<bestNumber) {
	      better=true;
	    } else if (numberInfeasibilitiesUp[i]==bestNumber) {
	      if (min(changeUp[i],changeDown[i])<bestCriterion)
		better=true;;
	    }
	    if (better) {
	      // see which way
	      if (changeUp[i]<=changeDown[i])
		betterWay=1;
	      else
		betterWay=-1;
	    }
	  }
	  if (betterWay) {
	    bestCriterion = min(changeUp[i],changeDown[i]);
	    bestNumber = thisNumber;
	    whichObject = i;
	    bestWay = betterWay;
	  }
	}
      }
      break;
    case 1:
      for ( i = 0 ; i < numberObjects ; i++) {
	int betterWay=0;
	if (changeUp[i]<=changeDown[i]) {
	  if (changeUp[i]>bestCriterion)
	    betterWay=1;
	} else {
	  if (changeDown[i]>bestCriterion)
	    betterWay=-1;
	}
	if (betterWay) {
	  bestCriterion = min(changeUp[i],changeDown[i]);
	  whichObject = i;
	  bestWay = betterWay;
	}
      }
      break;
    case 2:
      for ( i = 0 ; i < numberObjects ; i++) {
	double change = min(changeUp[i],changeDown[i]);
	double sum = changeUp[i] + changeDown[i];
	bool take=false;
	if (change>1.1*bestCriterion) 
	  take=true;
	else if (change>0.9*bestCriterion&&sum+change>bestCriterion+alternativeCriterion) 
	  take=true;
	if (take) {
	  if (changeUp[i]<=changeDown[i]) {
	    if (changeUp[i]>bestCriterion)
	      bestWay=1;
	  } else {
	    if (changeDown[i]>bestCriterion)
	      bestWay=-1;
	  }
	  bestCriterion = change;
	  alternativeCriterion = sum;
	  whichObject = i;
	}
      }
      break;
    case 3:
      for ( i = 0 ; i < numberObjects ; i++) {
	int numberNext = numberInfeasibilitiesUp[i];
	
	if (numberNext<numberUnsatisfied) {
	  int numberUp = numberUnsatisfied - numberInfeasibilitiesUp[i];
	  double perUnsatisfied = changeUp[i]/(double) numberUp;
	  double estimatedObjective = objectiveValue + numberUnsatisfied * perUnsatisfied;
	  if (estimatedObjective<bestEstimate) {
	    bestEstimate = estimatedObjective;
	    bestWay=1;
	    whichObject=i;
	  }
	}
	numberNext = numberInfeasibilitiesDown[i];
	if (numberNext<numberUnsatisfied) {
	  int numberDown = numberUnsatisfied - numberInfeasibilitiesDown[i];
	  double perUnsatisfied = changeDown[i]/(double) numberDown;
	  double estimatedObjective = objectiveValue + numberUnsatisfied * perUnsatisfied;
	  if (estimatedObjective<bestEstimate) {
	    bestEstimate = estimatedObjective;
	    bestWay=-1;
	    whichObject=i;
	  }
	}
      }
      break;
    case 4:
      // if number infeas same then cheapest up
      // first get best number or when going down
      // now choose smallest change up amongst equal number infeas
      for ( i = 0 ; i < numberObjects ; i++) {
	int thisNumber = min(numberInfeasibilitiesUp[i],numberInfeasibilitiesDown[i]);
	if (thisNumber<=bestNumber) {
	  int betterWay=0;
	  if (numberInfeasibilitiesUp[i]<numberInfeasibilitiesDown[i]) {
	    if (numberInfeasibilitiesUp[i]<bestNumber) {
	      betterWay = 1;
	    } else {
	      if (changeUp[i]<bestCriterion)
		betterWay=1;
	    }
	  } else if (numberInfeasibilitiesUp[i]>numberInfeasibilitiesDown[i]) {
	    if (numberInfeasibilitiesDown[i]<bestNumber) {
	      betterWay = -1;
	    } else {
	      if (changeDown[i]<bestCriterion)
		betterWay=-1;
	    }
	  } else {
	    // up and down have same number
	    bool better=false;
	    if (numberInfeasibilitiesUp[i]<bestNumber) {
	      better=true;
	    } else if (numberInfeasibilitiesUp[i]==bestNumber) {
	      if (min(changeUp[i],changeDown[i])<bestCriterion)
		better=true;;
	    }
	    if (better) {
	      // see which way
	      if (changeUp[i]<=changeDown[i])
		betterWay=1;
	      else
		betterWay=-1;
	    }
	  }
	  if (betterWay) {
	    bestCriterion = min(changeUp[i],changeDown[i]);
	    bestNumber = thisNumber;
	    whichObject = i;
	    bestWay = betterWay;
	  }
	}
      }
      bestCriterion=1.0e50;
      for ( i = 0 ; i < numberObjects ; i++) {
	int thisNumber = numberInfeasibilitiesUp[i];
	if (thisNumber==bestNumber&&changeUp) {
	  if (changeUp[i]<bestCriterion) {
	    bestCriterion = changeUp[i];
	    whichObject = i;
	    bestWay = 1;
	  }
	}
      }
      break;
    }
    // set way in best
    if (whichObject>=0) {
      CbcBranchingObject * bestObject = objects[whichObject];
      if (bestObject->object()&&bestObject->object()->preferredWay()) 
        bestWay = bestObject->object()->preferredWay();
      bestObject->way(bestWay);
    } else {
      printf("debug\n");
    }
  }
  return whichObject;
}

// Default Constructor 
CbcFollowOn::CbcFollowOn ()
  : CbcObject(),
    rhs_(NULL)
{
}

// Useful constructor
CbcFollowOn::CbcFollowOn (CbcModel * model)
  : CbcObject(model)
{
  assert (model);
  OsiSolverInterface * solver = model_->solver();
  matrix_ = *solver->getMatrixByCol();
  matrix_.removeGaps();
  matrix_.setExtraGap(0.0);
  matrixByRow_ = *solver->getMatrixByRow();
  int numberRows = matrix_.getNumRows();
  
  rhs_ = new int[numberRows];
  int i;
  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();
  // Row copy
  const double * elementByRow = matrixByRow_.getElements();
  const int * column = matrixByRow_.getIndices();
  const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
  const int * rowLength = matrixByRow_.getVectorLengths();
  for (i=0;i<numberRows;i++) {
    rhs_[i]=0;
    double value = rowLower[i];
    if (value==rowUpper[i]) {
      if (floor(value)==value&&value>=1.0&&value<10.0) {
	// check elements
	bool good=true;
	for (int j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	  int iColumn = column[j];
	  if (!solver->isBinary(iColumn))
	    good=false;
	  double elValue = elementByRow[j];
	  if (floor(elValue)!=elValue||value<1.0)
	    good=false;
	}
	if (good)
	  rhs_[i]=(int) value;
      }
    }
  }
}

// Copy constructor 
CbcFollowOn::CbcFollowOn ( const CbcFollowOn & rhs)
  :CbcObject(rhs),
   matrix_(rhs.matrix_),
   matrixByRow_(rhs.matrixByRow_)
{
  int numberRows = matrix_.getNumRows();
  rhs_= CoinCopyOfArray(rhs.rhs_,numberRows);
}

// Clone
CbcObject *
CbcFollowOn::clone() const
{
  return new CbcFollowOn(*this);
}

// Assignment operator 
CbcFollowOn & 
CbcFollowOn::operator=( const CbcFollowOn& rhs)
{
  if (this!=&rhs) {
    CbcObject::operator=(rhs);
    delete [] rhs_;
    matrix_ = rhs.matrix_;
    matrixByRow_ = rhs.matrixByRow_;
    int numberRows = matrix_.getNumRows();
    rhs_= CoinCopyOfArray(rhs.rhs_,numberRows);
  }
  return *this;
}

// Destructor 
CbcFollowOn::~CbcFollowOn ()
{
  delete [] rhs_;
}
// As some computation is needed in more than one place - returns row
int 
CbcFollowOn::gutsOfFollowOn(int & otherRow, int & preferredWay) const
{
  int whichRow=-1;
  otherRow=-1;
  int numberRows = matrix_.getNumRows();
  
  int i;
  // For sorting
  int * sort = new int [numberRows];
  int * isort = new int [numberRows];
  // Column copy
  //const double * element = matrix_.getElements();
  const int * row = matrix_.getIndices();
  const CoinBigIndex * columnStart = matrix_.getVectorStarts();
  const int * columnLength = matrix_.getVectorLengths();
  // Row copy
  const double * elementByRow = matrixByRow_.getElements();
  const int * column = matrixByRow_.getIndices();
  const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
  const int * rowLength = matrixByRow_.getVectorLengths();
  OsiSolverInterface * solver = model_->solver();
  const double * columnLower = solver->getColLower();
  const double * columnUpper = solver->getColUpper();
  const double * solution = solver->getColSolution();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  int nSort=0;
  for (i=0;i<numberRows;i++) {
    if (rhs_[i]) {
      // check elements
      double smallest=1.0e10;
      double largest=0.0;
      int rhsValue=rhs_[i];
      int number1=0;
      int numberUnsatisfied=0;
      for (int j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	int iColumn = column[j];
	double value = elementByRow[j];
	double solValue = solution[iColumn];
	if (columnLower[iColumn]!=columnUpper[iColumn]) {
	  smallest = CoinMin(smallest,value);
	  largest = CoinMax(largest,value);
	  if (value==1.0)
	    number1++;
	  if (solValue<1.0-integerTolerance&&solValue>integerTolerance)
	    numberUnsatisfied++;
	} else {
	  rhsValue -= (int)(value*floor(solValue+0.5));
	}
      }
      if (numberUnsatisfied>1) {
	if (smallest<largest) {
	  // probably no good but check a few things
	  assert (largest<=rhsValue);
	  if (number1==1&&largest==rhsValue)
	    printf("could fix\n");
	} else if (largest==rhsValue) {
	  sort[nSort]=i;
	  isort[nSort++]=-numberUnsatisfied;
	}
      }
    }
  }
  if (nSort>1) {
    CoinSort_2(isort,isort+nSort,sort);
    CoinZeroN(isort,numberRows);
    double * other = new double[numberRows];
    CoinZeroN(other,numberRows);
    int * which = new int[numberRows];
    //#define COUNT
#ifndef COUNT
    bool beforeSolution = model_->getSolutionCount()==0;
#endif
    for (int k=0;k<nSort-1;k++) {
      i=sort[k];
      int numberUnsatisfied = 0;
      int n=0;
      int j;
      for (j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	int iColumn = column[j];
	if (columnLower[iColumn]!=columnUpper[iColumn]) {
	  double solValue = solution[iColumn]-columnLower[iColumn];
	  if (solValue<1.0-integerTolerance&&solValue>integerTolerance) {
	    numberUnsatisfied++;
	    for (int jj=columnStart[iColumn];jj<columnStart[iColumn]+columnLength[iColumn];jj++) {
	      int iRow = row[jj];
	      if (rhs_[iRow]) {
		other[iRow]+=solValue;
		if (isort[iRow]) {
		  isort[iRow]++;
		} else {
		  isort[iRow]=1;
		  which[n++]=iRow;
		}
	      }
	    }
	  }
	}
      }
      double total=0.0;
      // Take out row
      double sumThis=other[i];
      other[i]=0.0;
      assert (numberUnsatisfied==isort[i]);
      // find one nearest half if solution, one if before solution
      int iBest=-1;
      double dtarget=0.5*total;
#ifdef COUNT
      int target = (numberUnsatisfied+1)>>1;
      int best=numberUnsatisfied;
#else
      double best;
      if (beforeSolution)
	best=dtarget;
      else
	best=1.0e30;
#endif
      for (j=0;j<n;j++) {
	int iRow = which[j];
	double dvalue=other[iRow];
	other[iRow]=0.0;
#ifdef COUNT
	int value = isort[iRow];
#endif
	isort[iRow]=0;
        if (fabs(dvalue)<1.0e-8||fabs(sumThis-dvalue)<1.0e-8)
          continue;
	if (dvalue<integerTolerance||dvalue>1.0-integerTolerance)
	  continue;
#ifdef COUNT
	if (abs(value-target)<best&&value!=numberUnsatisfied) {
	  best=abs(value-target);
	  iBest=iRow;
	  if (dvalue<dtarget)
	    preferredWay=1;
	  else
	    preferredWay=-1;
	}
#else
	if (beforeSolution) {
	  if (fabs(dvalue-dtarget)>best) {
	    best = fabs(dvalue-dtarget);
	    iBest=iRow;
	    if (dvalue<dtarget)
	      preferredWay=1;
	    else
	      preferredWay=-1;
	  }
	} else {
	  if (fabs(dvalue-dtarget)<best) {
	    best = fabs(dvalue-dtarget);
	    iBest=iRow;
	    if (dvalue<dtarget)
	      preferredWay=1;
	    else
	      preferredWay=-1;
	  }
	}
#endif
      }
      if (iBest>=0) {
	whichRow=i;
	otherRow=iBest;
	break;
      }
    }
    delete [] which;
    delete [] other;
  }
  delete [] sort;
  delete [] isort;
  return whichRow;
}

// Infeasibility - large is 0.5
double 
CbcFollowOn::infeasibility(int & preferredWay) const
{
  int otherRow=0;
  int whichRow = gutsOfFollowOn(otherRow,preferredWay);
  if (whichRow<0)
    return 0.0;
  else
  return 2.0* model_->getDblParam(CbcModel::CbcIntegerTolerance);
}

// This looks at solution and sets bounds to contain solution
void 
CbcFollowOn::feasibleRegion()
{
}


// Creates a branching object
CbcBranchingObject * 
CbcFollowOn::createBranch(int way) 
{
  int otherRow=0;
  int preferredWay;
  int whichRow = gutsOfFollowOn(otherRow,preferredWay);
  assert(way==preferredWay);
  assert (whichRow>=0);
  int numberColumns = matrix_.getNumCols();
  
  // Column copy
  //const double * element = matrix_.getElements();
  const int * row = matrix_.getIndices();
  const CoinBigIndex * columnStart = matrix_.getVectorStarts();
  const int * columnLength = matrix_.getVectorLengths();
  // Row copy
  //const double * elementByRow = matrixByRow_.getElements();
  const int * column = matrixByRow_.getIndices();
  const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
  const int * rowLength = matrixByRow_.getVectorLengths();
  OsiSolverInterface * solver = model_->solver();
  const double * columnLower = solver->getColLower();
  const double * columnUpper = solver->getColUpper();
  //const double * solution = solver->getColSolution();
  int nUp=0;
  int nDown=0;
  int * upList = new int[numberColumns];
  int * downList = new int[numberColumns];
  int j;
  for (j=rowStart[whichRow];j<rowStart[whichRow]+rowLength[whichRow];j++) {
    int iColumn = column[j];
    if (columnLower[iColumn]!=columnUpper[iColumn]) {
      bool up=true;
      for (int jj=columnStart[iColumn];jj<columnStart[iColumn]+columnLength[iColumn];jj++) {
	int iRow = row[jj];
	if (iRow==otherRow) {
	  up=false;
	  break;
	}
      }
      if (up)
	upList[nUp++]=iColumn;
      else
	downList[nDown++]=iColumn;
    }
  }
  //printf("way %d\n",way);
  // create object
  //printf("would fix %d down and %d up\n",nDown,nUp);
  CbcBranchingObject * branch
     = new CbcFixingBranchingObject(model_,way,
					 nDown,downList,nUp,upList);
  delete [] upList;
  delete [] downList;
  return branch;
}
// Default Constructor 
CbcFixingBranchingObject::CbcFixingBranchingObject()
  :CbcBranchingObject()
{
  numberDown_=0;
  numberUp_=0;
  downList_=NULL;
  upList_=NULL;
}

// Useful constructor
CbcFixingBranchingObject::CbcFixingBranchingObject (CbcModel * model,
						    int way ,
						    int numberOnDownSide, const int * down,
						    int numberOnUpSide, const int * up)
  :CbcBranchingObject(model,0,way,0.5)
{
  numberDown_=numberOnDownSide;
  numberUp_=numberOnUpSide;
  downList_ = CoinCopyOfArray(down,numberDown_);
  upList_ = CoinCopyOfArray(up,numberUp_);
}

// Copy constructor 
CbcFixingBranchingObject::CbcFixingBranchingObject ( const CbcFixingBranchingObject & rhs) :CbcBranchingObject(rhs)
{
  numberDown_=rhs.numberDown_;
  numberUp_=rhs.numberUp_;
  downList_ = CoinCopyOfArray(rhs.downList_,numberDown_);
  upList_ = CoinCopyOfArray(rhs.upList_,numberUp_);
}

// Assignment operator 
CbcFixingBranchingObject & 
CbcFixingBranchingObject::operator=( const CbcFixingBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    delete [] downList_;
    delete [] upList_;
    numberDown_=rhs.numberDown_;
    numberUp_=rhs.numberUp_;
    downList_ = CoinCopyOfArray(rhs.downList_,numberDown_);
    upList_ = CoinCopyOfArray(rhs.upList_,numberUp_);
  }
  return *this;
}
CbcBranchingObject * 
CbcFixingBranchingObject::clone() const
{ 
  return (new CbcFixingBranchingObject(*this));
}


// Destructor 
CbcFixingBranchingObject::~CbcFixingBranchingObject ()
{
  delete [] downList_;
  delete [] upList_;
}
double
CbcFixingBranchingObject::branch()
{
  decrementNumberBranchesLeft();
  OsiSolverInterface * solver = model_->solver();
  const double * columnLower = solver->getColLower();
  int i;
  // *** for way - up means fix all those in up section
  if (way_<0) {
#ifdef FULL_PRINT
    printf("Down Fix ");
#endif
    //printf("Down Fix %d\n",numberDown_);
    for (i=0;i<numberDown_;i++) {
      int iColumn = downList_[i];
      model_->solver()->setColUpper(iColumn,columnLower[iColumn]);
#ifdef FULL_PRINT
      printf("Setting bound on %d to lower bound\n",iColumn);
#endif
    }
    way_=1;	  // Swap direction
  } else {
#ifdef FULL_PRINT
    printf("Up Fix ");
#endif
    //printf("Up Fix %d\n",numberUp_);
    for (i=0;i<numberUp_;i++) {
      int iColumn = upList_[i];
      model_->solver()->setColUpper(iColumn,columnLower[iColumn]);
#ifdef FULL_PRINT
      printf("Setting bound on %d to lower bound\n",iColumn);
#endif
    }
    way_=-1;	  // Swap direction
  }
#ifdef FULL_PRINT
  printf("\n");
#endif
  return 0.0;
}
void
CbcFixingBranchingObject::print()
{
  int i;
  // *** for way - up means fix all those in up section
  if (way_<0) {
    printf("Down Fix ");
    for (i=0;i<numberDown_;i++) {
      int iColumn = downList_[i];
      printf("%d ",iColumn);
    }
  } else {
    printf("Up Fix ");
    for (i=0;i<numberUp_;i++) {
      int iColumn = upList_[i];
      printf("%d ",iColumn);
    }
  }
  printf("\n");
}
// Default Constructor 
CbcNWay::CbcNWay ()
  : CbcObject(),
    numberMembers_(0),
    members_(NULL),
    consequence_(NULL)
{
}

// Useful constructor (which are integer indices)
CbcNWay::CbcNWay (CbcModel * model, int numberMembers,
                  const int * which, int identifier)
  : CbcObject(model)
{
  id_=identifier;
  numberMembers_=numberMembers;
  if (numberMembers_) {
    members_ = new int[numberMembers_];
    memcpy(members_,which,numberMembers_*sizeof(int));
  } else {
    members_ = NULL;
  }
  consequence_ = NULL;
}

// Copy constructor 
CbcNWay::CbcNWay ( const CbcNWay & rhs)
  :CbcObject(rhs)
{
  numberMembers_ = rhs.numberMembers_;
  consequence_ = NULL;
  if (numberMembers_) {
    members_ = new int[numberMembers_];
    memcpy(members_,rhs.members_,numberMembers_*sizeof(int));
    if (rhs.consequence_) {
      consequence_ = new CbcConsequence * [numberMembers_];
      for (int i=0;i<numberMembers_;i++) {
        if (rhs.consequence_[i])
          consequence_[i]= rhs.consequence_[i]->clone();
        else
          consequence_[i]=NULL;
      }
    }
  } else {
    members_ = NULL;
  }
}

// Clone
CbcObject *
CbcNWay::clone() const
{
  return new CbcNWay(*this);
}

// Assignment operator 
CbcNWay & 
CbcNWay::operator=( const CbcNWay& rhs)
{
  if (this!=&rhs) {
    CbcObject::operator=(rhs);
    delete [] members_;
    numberMembers_ = rhs.numberMembers_;
    if (consequence_) {
      for (int i=0;i<numberMembers_;i++) 
        delete consequence_[i];
      delete [] consequence_;
      consequence_=NULL;
    }
    if (numberMembers_) {
      members_ = new int[numberMembers_];
      memcpy(members_,rhs.members_,numberMembers_*sizeof(int));
    } else {
      members_ = NULL;
    }
    if (rhs.consequence_) {
      consequence_ = new CbcConsequence * [numberMembers_];
      for (int i=0;i<numberMembers_;i++) {
        if (rhs.consequence_[i])
          consequence_[i]= rhs.consequence_[i]->clone();
        else
          consequence_[i]=NULL;
      }
    }
  }
  return *this;
}

// Destructor 
CbcNWay::~CbcNWay ()
{
  delete [] members_;
  if (consequence_) {
    for (int i=0;i<numberMembers_;i++) 
      delete consequence_[i];
    delete [] consequence_;
  }
}
// Set up a consequence for a single member
void 
CbcNWay::setConsequence(int iColumn, const CbcConsequence & consequence)
{
  if (!consequence_) {
    consequence_ = new CbcConsequence * [numberMembers_];
    for (int i=0;i<numberMembers_;i++) 
      consequence_[i]=NULL;
  }
  for (int i=0;i<numberMembers_;i++) {
    if (members_[i]==iColumn) {
      consequence_[i]=consequence.clone();
      break;
    }
  }
}

// Applies a consequence for a single member
void 
CbcNWay::applyConsequence(int iSequence, int state) const
{
  assert (state==-9999||state==9999);
  if (consequence_) {
    CbcConsequence * consequence = consequence_[iSequence];
    if (consequence) 
      consequence->applyToSolver(model_->solver(),state);
  }
}
  
// Infeasibility - large is 0.5
double 
CbcNWay::infeasibility(int & preferredWay) const
{
  int numberUnsatis=0;
  int j;
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  double largestValue=0.0;
  
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);

  for (j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    double value = solution[iColumn];
    value = CoinMax(value, lower[iColumn]);
    value = CoinMin(value, upper[iColumn]);
    double distance = CoinMin(value-lower[iColumn],upper[iColumn]-value);
    if (distance>integerTolerance) {
      numberUnsatis++;
      largestValue = CoinMax(distance,largestValue);
    }
  }
  preferredWay=1;
  if (numberUnsatis) {
    return largestValue;
  } else {
    return 0.0; // satisfied
  }
}

// This looks at solution and sets bounds to contain solution
void 
CbcNWay::feasibleRegion()
{
  int j;
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  for (j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    double value = solution[iColumn];
    value = CoinMax(value, lower[iColumn]);
    value = CoinMin(value, upper[iColumn]);
    if (value>=upper[iColumn]-integerTolerance) {
      solver->setColLower(iColumn,upper[iColumn]);
    } else {
      assert (value<=lower[iColumn]+integerTolerance);
      solver->setColUpper(iColumn,lower[iColumn]);
    }
  }
}
// Redoes data when sequence numbers change
void 
CbcNWay::redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns)
{
  model_=model;
  int n2=0;
  for (int j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    int i;
    for (i=0;i<numberColumns;i++) {
      if (originalColumns[i]==iColumn)
        break;
    }
    if (i<numberColumns) {
      members_[n2]=i;
      consequence_[n2++]=consequence_[j];
    } else {
      delete consequence_[j];
    }
  }
  if (n2<numberMembers_) {
    printf("** NWay number of members reduced from %d to %d!\n",numberMembers_,n2);
    numberMembers_=n2;
  }
}


// Creates a branching object
CbcBranchingObject * 
CbcNWay::createBranch(int way) 
{
  int numberFree=0;
  int j;

  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  int * list = new int[numberMembers_];
  double * sort = new double[numberMembers_];

  for (j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    double value = solution[iColumn];
    value = CoinMax(value, lower[iColumn]);
    value = CoinMin(value, upper[iColumn]);
    if (upper[iColumn]>lower[iColumn]) {
      double distance = upper[iColumn]-value;
      list[numberFree]=j;
      sort[numberFree++]=distance;
    }
  }
  assert (numberFree);
  // sort
  CoinSort_2(sort,sort+numberFree,list);
  // create object
  CbcBranchingObject * branch;
  branch = new CbcNWayBranchingObject(model_,this,numberFree,list);
  branch->setOriginalObject(this);
  delete [] list;
  delete [] sort;
  return branch;
}
  
// Default Constructor 
CbcNWayBranchingObject::CbcNWayBranchingObject()
  :CbcBranchingObject()
{
  order_=NULL;
  object_=NULL;
  numberInSet_=0;
  way_=0;
}

// Useful constructor
CbcNWayBranchingObject::CbcNWayBranchingObject (CbcModel * model,
                                                const CbcNWay * nway, 
                                                int number, const int * order)
  :CbcBranchingObject(model,nway->id(),-1,0.5)
{
  numberBranches_ = number;
  order_ = new int [number];
  object_=nway;
  numberInSet_=number;
  memcpy(order_,order,number*sizeof(int));
}

// Copy constructor 
CbcNWayBranchingObject::CbcNWayBranchingObject ( const CbcNWayBranchingObject & rhs) :CbcBranchingObject(rhs)
{
  numberInSet_=rhs.numberInSet_;
  object_=rhs.object_;
  if (numberInSet_) {
    order_ = new int [numberInSet_];
    memcpy(order_,rhs.order_,numberInSet_*sizeof(int));
  } else {
    order_=NULL;
  }    
}

// Assignment operator 
CbcNWayBranchingObject & 
CbcNWayBranchingObject::operator=( const CbcNWayBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    object_=rhs.object_;
    delete [] order_;
    numberInSet_=rhs.numberInSet_;
    if (numberInSet_) {
      order_ = new int [numberInSet_];
      memcpy(order_,rhs.order_,numberInSet_*sizeof(int));
    } else {
      order_=NULL;
    }    
  }
  return *this;
}
CbcBranchingObject * 
CbcNWayBranchingObject::clone() const
{ 
  return (new CbcNWayBranchingObject(*this));
}


// Destructor 
CbcNWayBranchingObject::~CbcNWayBranchingObject ()
{
  delete [] order_;
}
double
CbcNWayBranchingObject::branch()
{
  int which = branchIndex_;
  branchIndex_++;
  assert (numberBranchesLeft()>=0);
  if (which==0) {
    // first branch so way_ may mean something
    assert (way_==-1||way_==1);
    if (way_==-1)
      which++;
  } else if (which==1) {
    // second branch so way_ may mean something
    assert (way_==-1||way_==1);
    if (way_==-1)
      which--;
    // switch way off
    way_=0;
  }
  const double * lower = model_->solver()->getColLower();
  const double * upper = model_->solver()->getColUpper();
  const int * members = object_->members();
  for (int j=0;j<numberInSet_;j++) {
    int iSequence = order_[j];
    int iColumn = members[iSequence];
    if (j!=which) {
      model_->solver()->setColUpper(iColumn,lower[iColumn]);
      //model_->solver()->setColLower(iColumn,lower[iColumn]);
      assert (lower[iColumn]>-1.0e20);
      // apply any consequences
      object_->applyConsequence(iSequence,-9999);
    } else {
      model_->solver()->setColLower(iColumn,upper[iColumn]);
      //model_->solver()->setColUpper(iColumn,upper[iColumn]);
#ifdef FULL_PRINT
      printf("Up Fix %d to %g\n",iColumn,upper[iColumn]);
#endif
      assert (upper[iColumn]<1.0e20);
      // apply any consequences
      object_->applyConsequence(iSequence,9999);
    }
  }
  return 0.0;
}
void
CbcNWayBranchingObject::print()
{
  printf("NWay - Up Fix ");
  const int * members = object_->members();
  for (int j=0;j<way_;j++) {
    int iColumn = members[order_[j]];
    printf("%d ",iColumn);
  }
  printf("\n");
}

// Default Constructor 
CbcFixVariable::CbcFixVariable ()
  : CbcConsequence(),
    numberStates_(0),
    states_(NULL),
    startLower_(NULL),
    startUpper_(NULL),
    newBound_(NULL),
    variable_(NULL)
{
}

// One useful Constructor 
CbcFixVariable::CbcFixVariable (int numberStates,const int * states, const int * numberNewLower, 
                                const int ** newLowerValue,
                                const int ** lowerColumn,
                                const int * numberNewUpper, const int ** newUpperValue,
                                const int ** upperColumn)
  : CbcConsequence(),
    states_(NULL),
    startLower_(NULL),
    startUpper_(NULL),
    newBound_(NULL),
    variable_(NULL)
{
  // How much space
  numberStates_ = numberStates;
  if (numberStates_) {
    states_ = new int[numberStates_];
    memcpy(states_,states,numberStates_*sizeof(int));
    int i;
    int n=0;
    startLower_ = new int[numberStates_+1];
    startUpper_ = new int[numberStates_+1];
    startLower_[0]=0;
    //count
    for (i=0;i<numberStates_;i++) {
      n += numberNewLower[i];
      startUpper_[i]=n;
      n += numberNewUpper[i];
      startLower_[i+1]=n;
    }
    newBound_ = new double [n];
    variable_ = new int [n];
    n=0;
    for (i=0;i<numberStates_;i++) {
      int j;
      int k;
      const int * bound;
      const int * variable;
      k=numberNewLower[i];
      bound = newLowerValue[i];
      variable = lowerColumn[i];
      for (j=0;j<k;j++) {
        newBound_[n]=bound[j];
        variable_[n++]=variable[j];
      }
      k=numberNewUpper[i];
      bound = newUpperValue[i];
      variable = upperColumn[i];
      for (j=0;j<k;j++) {
        newBound_[n]=bound[j];
        variable_[n++]=variable[j];
      }
    }
  }
}

// Copy constructor 
CbcFixVariable::CbcFixVariable ( const CbcFixVariable & rhs)
  :CbcConsequence(rhs)
{
  numberStates_ = rhs.numberStates_;
  states_ = NULL;
  startLower_ = NULL;
  startUpper_ = NULL;
  newBound_ = NULL;
  variable_ = NULL;
  if (numberStates_) {
    states_ = CoinCopyOfArray(rhs.states_,numberStates_);
    startLower_ = CoinCopyOfArray(rhs.startLower_,numberStates_+1);
    startUpper_ = CoinCopyOfArray(rhs.startUpper_,numberStates_+1);
    int n=startLower_[numberStates_];
    newBound_ = CoinCopyOfArray(rhs.newBound_,n);
    variable_ = CoinCopyOfArray(rhs.variable_,n);
  }
}

// Clone
CbcConsequence *
CbcFixVariable::clone() const
{
  return new CbcFixVariable(*this);
}

// Assignment operator 
CbcFixVariable & 
CbcFixVariable::operator=( const CbcFixVariable& rhs)
{
  if (this!=&rhs) {
    CbcConsequence::operator=(rhs);
    delete [] states_;
    delete [] startLower_;
    delete [] startUpper_;
    delete [] newBound_;
    delete [] variable_;
    states_ = NULL;
    startLower_ = NULL;
    startUpper_ = NULL;
    newBound_ = NULL;
    variable_ = NULL;
    numberStates_ = rhs.numberStates_;
    if (numberStates_) {
      states_ = CoinCopyOfArray(rhs.states_,numberStates_);
      startLower_ = CoinCopyOfArray(rhs.startLower_,numberStates_+1);
      startUpper_ = CoinCopyOfArray(rhs.startUpper_,numberStates_+1);
      int n=startLower_[numberStates_];
      newBound_ = CoinCopyOfArray(rhs.newBound_,n);
      variable_ = CoinCopyOfArray(rhs.variable_,n);
    }
  }
  return *this;
}

// Destructor 
CbcFixVariable::~CbcFixVariable ()
{
  delete [] states_;
  delete [] startLower_;
  delete [] startUpper_;
  delete [] newBound_;
  delete [] variable_;
}
// Set up a startLower for a single member
void 
CbcFixVariable::applyToSolver(OsiSolverInterface * solver, int state) const
{
  assert (state==-9999||state==9999);
  // Find state
  int find;
  for (find=0;find<numberStates_;find++) 
    if (states_[find]==state)
      break;
  if (find==numberStates_)
    return;
  int i;
  // Set new lower bounds
  for (i=startLower_[find];i<startUpper_[find];i++) {
    int iColumn = variable_[i];
    double value = newBound_[i];
    double oldValue = solver->getColLower()[iColumn];
    //printf("for %d old lower bound %g, new %g",iColumn,oldValue,value);
    solver->setColLower(iColumn,CoinMax(value,oldValue));
    //printf(" => %g\n",solver->getColLower()[iColumn]);
  }
  // Set new upper bounds
  for (i=startUpper_[find];i<startLower_[find+1];i++) {
    int iColumn = variable_[i];
    double value = newBound_[i];
    double oldValue = solver->getColUpper()[iColumn];
    //printf("for %d old upper bound %g, new %g",iColumn,oldValue,value);
    solver->setColUpper(iColumn,CoinMin(value,oldValue));
    //printf(" => %g\n",solver->getColUpper()[iColumn]);
  }
}

// Default Constructor 
CbcDummyBranchingObject::CbcDummyBranchingObject(CbcModel * model)
  :CbcBranchingObject(model,0,0,0.5)
{
  setNumberBranchesLeft(1);
}


// Copy constructor 
CbcDummyBranchingObject::CbcDummyBranchingObject ( const CbcDummyBranchingObject & rhs) :CbcBranchingObject(rhs)
{
}

// Assignment operator 
CbcDummyBranchingObject & 
CbcDummyBranchingObject::operator=( const CbcDummyBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
  }
  return *this;
}
CbcBranchingObject * 
CbcDummyBranchingObject::clone() const
{ 
  return (new CbcDummyBranchingObject(*this));
}


// Destructor 
CbcDummyBranchingObject::~CbcDummyBranchingObject ()
{
}

/*
  Perform a dummy branch
*/
double
CbcDummyBranchingObject::branch()
{
  decrementNumberBranchesLeft();
  return 0.0;
}
// Print what would happen  
void
CbcDummyBranchingObject::print()
{
  printf("Dummy branch\n");
}
