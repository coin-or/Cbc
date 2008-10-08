// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CoinTypes.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

//##############################################################################

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

//##############################################################################

// Default Constructor 
CbcSOS::CbcSOS ()
  : CbcObject(),
    members_(NULL),
    weights_(NULL),
    shadowEstimateDown_(1.0),
    shadowEstimateUp_(1.0),
    downDynamicPseudoRatio_(0.0),
    upDynamicPseudoRatio_(0.0),
    numberTimesDown_(0),
    numberTimesUp_(0),
    numberMembers_(0),
    sosType_(-1),
    integerValued_(false)
{
}

// Useful constructor (which are indices)
CbcSOS::CbcSOS (CbcModel * model,  int numberMembers,
	   const int * which, const double * weights, int identifier,int type)
  : CbcObject(model),
    shadowEstimateDown_(1.0),
    shadowEstimateUp_(1.0),
    downDynamicPseudoRatio_(0.0),
    upDynamicPseudoRatio_(0.0),
    numberTimesDown_(0),
    numberTimesUp_(0),
    numberMembers_(numberMembers),
    sosType_(type)
{
  id_=identifier;
  integerValued_ = type==1;
  if (integerValued_) {
    // check all members integer
    OsiSolverInterface * solver = model->solver();
    if (solver) {
      for (int i=0;i<numberMembers_;i++) {
	if (!solver->isInteger(which[i]))
	  integerValued_=false;
      }
    } else {
      // can't tell
      integerValued_=false;
    }
  }
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
  shadowEstimateDown_ = rhs.shadowEstimateDown_;
  shadowEstimateUp_ = rhs.shadowEstimateUp_;
  downDynamicPseudoRatio_ = rhs.downDynamicPseudoRatio_;
  upDynamicPseudoRatio_ = rhs.upDynamicPseudoRatio_;
  numberTimesDown_ = rhs.numberTimesDown_;
  numberTimesUp_ = rhs.numberTimesUp_;
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
    shadowEstimateDown_ = rhs.shadowEstimateDown_;
    shadowEstimateUp_ = rhs.shadowEstimateUp_;
    downDynamicPseudoRatio_ = rhs.downDynamicPseudoRatio_;
    upDynamicPseudoRatio_ = rhs.upDynamicPseudoRatio_;
    numberTimesDown_ = rhs.numberTimesDown_;
    numberTimesUp_ = rhs.numberTimesUp_;
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
    // adjust??
    return value;
  } else {
    return 0.0; // satisfied
  }
}
// Infeasibility - large is 0.5
double 
CbcSOS::infeasibility(const OsiBranchingInformation * info, 
		      int & preferredWay) const
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
    if (info->defaultDual_>=0.0&&info->usefulRegion_&&info->columnStart_) {
      assert (sosType_==1);
      int iWhere;
      for (iWhere=firstNonZero;iWhere<lastNonZero-1;iWhere++) {
	if (weight<weights_[iWhere+1]) {
	  break;
	}
      }
      int jColumnDown = members_[iWhere];
      int jColumnUp = members_[iWhere+1];
      int n=0;
      CoinBigIndex j;
      double objMove = info->objective_[jColumnDown];
      for (j=info->columnStart_[jColumnDown];
	   j<info->columnStart_[jColumnDown]+info->columnLength_[jColumnDown];j++) {
	double value = info->elementByColumn_[j];
	int iRow = info->row_[j];
	info->indexRegion_[n++]=iRow;
	info->usefulRegion_[iRow]=value;
      }
      for (iWhere=firstNonZero;iWhere<lastNonZero;iWhere++) {
	int jColumn = members_[iWhere];
	double solValue = info->solution_[jColumn];
	if (!solValue)
	  continue;
	objMove -= info->objective_[jColumn]*solValue;
	for (j=info->columnStart_[jColumn];
	     j<info->columnStart_[jColumn]+info->columnLength_[jColumn];j++) {
	  double value = -info->elementByColumn_[j]*solValue;
	  int iRow = info->row_[j];
	  double oldValue = info->usefulRegion_[iRow];
	  if (!oldValue) {
	    info->indexRegion_[n++]=iRow;
	  } else {
	    value += oldValue;
	    if (!value)
	      value = 1.0e-100;
	  }
	  info->usefulRegion_[iRow]=value;
	}
      }
      const double * pi = info->pi_;
      const double * activity = info->rowActivity_;
      const double * lower = info->rowLower_;
      const double * upper = info->rowUpper_;
      double tolerance = info->primalTolerance_;
      double direction = info->direction_;
      shadowEstimateDown_ = objMove*direction;
      bool infeasible=false;
      for (int k=0;k<n;k++) {
	int iRow = info->indexRegion_[k];
	double movement=info->usefulRegion_[iRow];
	// not this time info->usefulRegion_[iRow]=0.0;
	if (lower[iRow]<-1.0e20) 
	  assert (pi[iRow]<=1.0e-3);
	if (upper[iRow]>1.0e20) 
	  assert (pi[iRow]>=-1.0e-3);
	double valueP = pi[iRow]*direction;
	// if move makes infeasible then make at least default
	double newValue = activity[iRow] + movement;
	if (newValue>upper[iRow]+tolerance||newValue<lower[iRow]-tolerance) {
	  shadowEstimateDown_ += fabs(movement)*CoinMax(fabs(valueP),info->defaultDual_);
	  infeasible=true;
	}
      }
      if (shadowEstimateDown_<info->integerTolerance_) {
	if (!infeasible) {
	  shadowEstimateDown_=1.0e-10;
#ifdef COIN_DEVELOP
	  printf("zero pseudoShadowPrice\n");
#endif
	} else
	  shadowEstimateDown_ = info->integerTolerance_;
      }
      // And other way
      // take off
      objMove -= info->objective_[jColumnDown];
      for (j=info->columnStart_[jColumnDown];
	   j<info->columnStart_[jColumnDown]+info->columnLength_[jColumnDown];j++) {
	double value = -info->elementByColumn_[j];
	int iRow = info->row_[j];
	double oldValue = info->usefulRegion_[iRow];
	if (!oldValue) {
	  info->indexRegion_[n++]=iRow;
	  } else {
	  value += oldValue;
	  if (!value)
	    value = 1.0e-100;
	}
	info->usefulRegion_[iRow]=value;
      }
      // add on
      objMove += info->objective_[jColumnUp];
      for (j=info->columnStart_[jColumnUp];
	   j<info->columnStart_[jColumnUp]+info->columnLength_[jColumnUp];j++) {
	double value = info->elementByColumn_[j];
	int iRow = info->row_[j];
	double oldValue = info->usefulRegion_[iRow];
	if (!oldValue) {
	  info->indexRegion_[n++]=iRow;
	  } else {
	  value += oldValue;
	  if (!value)
	    value = 1.0e-100;
	}
	info->usefulRegion_[iRow]=value;
      }
      shadowEstimateUp_ = objMove*direction;
      infeasible=false;
      for (int k=0;k<n;k++) {
	int iRow = info->indexRegion_[k];
	double movement=info->usefulRegion_[iRow];
	info->usefulRegion_[iRow]=0.0;
	if (lower[iRow]<-1.0e20) 
	  assert (pi[iRow]<=1.0e-3);
	if (upper[iRow]>1.0e20) 
	  assert (pi[iRow]>=-1.0e-3);
	double valueP = pi[iRow]*direction;
	// if move makes infeasible then make at least default
	double newValue = activity[iRow] + movement;
	if (newValue>upper[iRow]+tolerance||newValue<lower[iRow]-tolerance) {
	  shadowEstimateUp_ += fabs(movement)*CoinMax(fabs(valueP),info->defaultDual_);
	  infeasible=true;
	}
      }
      if (shadowEstimateUp_<info->integerTolerance_) {
	if (!infeasible) {
	  shadowEstimateUp_=1.0e-10;
#ifdef COIN_DEVELOP
	  printf("zero pseudoShadowPrice\n");
#endif
	} else
	  shadowEstimateUp_ = info->integerTolerance_;
      }
      // adjust
      double downCost = shadowEstimateDown_;
      double upCost = shadowEstimateUp_;
      if (numberTimesDown_) 
	downCost *= downDynamicPseudoRatio_/
	  ((double) numberTimesDown_);
      if (numberTimesUp_) 
	upCost *= upDynamicPseudoRatio_/
	  ((double) numberTimesUp_);
#define WEIGHT_AFTER 0.7
#define WEIGHT_BEFORE 0.1
      int stateOfSearch = model_->stateOfSearch()%10;
      double returnValue=0.0;
      double minValue = CoinMin(downCost,upCost);
      double maxValue = CoinMax(downCost,upCost);
      if (stateOfSearch<=2) {
	// no branching solution
	returnValue = WEIGHT_BEFORE*minValue + (1.0-WEIGHT_BEFORE)*maxValue;
      } else {
	returnValue = WEIGHT_AFTER*minValue + (1.0-WEIGHT_AFTER)*maxValue;
      }
#ifdef PRINT_SHADOW
      printf("%d id - down %d %g up %d %g shadow %g, %g returned %g\n",
	     id_,numberTimesDown_,downDynamicPseudoRatio_,
	     numberTimesUp_,upDynamicPseudoRatio_,shadowEstimateDown_,
	     shadowEstimateUp_,returnValue);
#endif
      return returnValue;
    } else {
      double value = lastNonZero-firstNonZero+1;
      value *= 0.5/((double) numberMembers_);
      return value;
    }
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
  branch->setOriginalObject(this);
  return branch;
}
/* Pass in information on branch just done and create CbcObjectUpdateData instance.
   If object does not need data then backward pointer will be NULL.
   Assumes can get information from solver */
CbcObjectUpdateData 
CbcSOS::createUpdateInformation(const OsiSolverInterface * solver, 
				const CbcNode * node,
				const CbcBranchingObject * branchingObject)
{
  double originalValue=node->objectiveValue();
  int originalUnsatisfied = node->numberUnsatisfied();
  double objectiveValue = solver->getObjValue()*solver->getObjSense();
  int unsatisfied=0;
  int i;
  //might be base model - doesn't matter
  int numberIntegers = model_->numberIntegers();;
  const double * solution = solver->getColSolution();
  //const double * lower = solver->getColLower();
  //const double * upper = solver->getColUpper();
  double change = CoinMax(0.0,objectiveValue-originalValue);
  int iStatus;
  if (solver->isProvenOptimal())
    iStatus=0; // optimal
  else if (solver->isIterationLimitReached()
           &&!solver->isDualObjectiveLimitReached())
    iStatus=2; // unknown 
  else
    iStatus=1; // infeasible

  bool feasible = iStatus!=1;
  if (feasible) {
    double integerTolerance = 
      model_->getDblParam(CbcModel::CbcIntegerTolerance);
    const int * integerVariable = model_->integerVariable();
    for (i=0;i<numberIntegers;i++) {
      int j=integerVariable[i];
      double value = solution[j];
      double nearest = floor(value+0.5);
      if (fabs(value-nearest)>integerTolerance) 
        unsatisfied++;
    }
  }
  int way = branchingObject->way();
  way = - way; // because after branch so moved on
  double value = branchingObject->value();
  CbcObjectUpdateData newData (this, way,
			       change, iStatus,
			       originalUnsatisfied-unsatisfied,value);
  newData.originalObjective_ = originalValue;
  // Solvers know about direction
  double direction = solver->getObjSense();
  solver->getDblParam(OsiDualObjectiveLimit,newData.cutoff_);
  newData.cutoff_ *= direction;
  return newData;
}
// Update object by CbcObjectUpdateData
void 
CbcSOS::updateInformation(const CbcObjectUpdateData & data)
{
  bool feasible = data.status_!=1;
  int way = data.way_;
  //double value = data.branchingValue_;
  double originalValue = data.originalObjective_;
  double change = data.change_;
  if (way<0) {
    // down
    if (!feasible) {
      double distanceToCutoff=0.0;
      //double objectiveValue = model_->getCurrentMinimizationObjValue();
      distanceToCutoff =  model_->getCutoff()  - originalValue;
      if (distanceToCutoff<1.0e20) 
	change = distanceToCutoff*2.0;
      else 
	change = (downDynamicPseudoRatio_*shadowEstimateDown_+1.0e-3)*10.0;
    } 
    change = CoinMax(1.0e-12*(1.0+fabs(originalValue)),change);
#ifdef PRINT_SHADOW
    if (numberTimesDown_)
      printf("Updating id %d - down change %g (true %g) - ndown %d estimated change %g - raw shadow estimate %g\n",
	     id_,change,data.change_,numberTimesDown_,shadowEstimateDown_*
	     (downDynamicPseudoRatio_/((double) numberTimesDown_)),
	     shadowEstimateDown_);
    else
      printf("Updating id %d - down change %g (true %g) - shadow estimate %g\n",
	     id_,change,data.change_,shadowEstimateDown_);
#endif
    numberTimesDown_++;
    downDynamicPseudoRatio_ += change/shadowEstimateDown_;
  } else {
    // up
    if (!feasible) {
      double distanceToCutoff=0.0;
      //double objectiveValue = model_->getCurrentMinimizationObjValue();
      distanceToCutoff =  model_->getCutoff()  - originalValue;
      if (distanceToCutoff<1.0e20) 
	change = distanceToCutoff*2.0;
      else 
	change = (upDynamicPseudoRatio_*shadowEstimateUp_+1.0e-3)*10.0;
    } 
    change = CoinMax(1.0e-12*(1.0+fabs(originalValue)),change);
#ifdef PRINT_SHADOW
    if (numberTimesUp_)
      printf("Updating id %d - up change %g (true %g) - nup %d estimated change %g - raw shadow estimate %g\n",
	     id_,change,data.change_,numberTimesUp_,shadowEstimateUp_*
	     (upDynamicPseudoRatio_/((double) numberTimesUp_)),
	     shadowEstimateUp_);
    else
      printf("Updating id %d - up change %g (true %g) - shadow estimate %g\n",
	     id_,change,data.change_,shadowEstimateUp_);
#endif
    numberTimesUp_++;
    upDynamicPseudoRatio_ += change/shadowEstimateUp_;
  }
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

//##############################################################################

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
#ifdef COIN_DEVELOP
  if (fabs(value-floor(value+0.5))>1.0e-5)
    printf("value for %d away from integer %g\n",columnNumber_,value);
#endif
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
  CbcIntegerBranchingObject * branch = new CbcIntegerBranchingObject(model_,0,-1,0.5);
  fillCreateBranch(branch,info,way);
  return branch;
}
// Fills in a created branching object
void 
CbcSimpleInteger::fillCreateBranch(CbcIntegerBranchingObject * branch, const OsiBranchingInformation * info, int way) 
{
  branch->setOriginalObject(this);
  double value = info->solution_[columnNumber_];
  value = CoinMax(value, info->lower_[columnNumber_]);
  value = CoinMin(value, info->upper_[columnNumber_]);
  assert (info->upper_[columnNumber_]>info->lower_[columnNumber_]);
  if (!info->hotstartSolution_&&priority_!=-999) {
#ifndef NDEBUG
    double nearest = floor(value+0.5);
    assert (fabs(value-nearest)>info->integerTolerance_);
#endif
  } else if (info->hotstartSolution_) {
    double targetValue = info->hotstartSolution_[columnNumber_];
    if (way>0)
      value = targetValue-0.1;
    else
      value = targetValue+0.1;
  } else {
    if (value<=info->lower_[columnNumber_])
      value += 0.1;
    else if (value>=info->upper_[columnNumber_])
      value -= 0.1;
  }
  assert (value>=info->lower_[columnNumber_]&&
	  value<=info->upper_[columnNumber_]);
  branch->fillPart(columnNumber_,way,value);
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
  assert (numberColumns>0);
  int iColumn;
#if 0
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (columnNumber_==originalColumns[iColumn])
      break;
  }
  assert (iColumn<numberColumns);
#else
  iColumn=originalColumns[columnNumber_];
  assert (iColumn>=0);
#endif
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

//##############################################################################

// Default Constructor 
CbcIntegerBranchingObject::CbcIntegerBranchingObject()
  :CbcBranchingObject()
{
  down_[0] = 0.0;
  down_[1] = 0.0;
  up_[0] = 0.0;
  up_[1] = 0.0;
#ifdef FUNNY_BRANCHING
  variables_ = NULL;
  newBounds_ = NULL;
  numberExtraChangedBounds_ = 0;
#endif
}
// Useful constructor
CbcIntegerBranchingObject::CbcIntegerBranchingObject (CbcModel * model, 
						      int variable, int way , double value)
  :CbcBranchingObject(model,variable,way,value)
{
  int iColumn = variable;
  assert (model_->solver()->getNumCols()>0);
  down_[0] = model_->solver()->getColLower()[iColumn];
  down_[1] = floor(value_);
  up_[0] = ceil(value_);
  up_[1] = model->getColUpper()[iColumn];
#ifdef FUNNY_BRANCHING
  variables_ = NULL;
  newBounds_ = NULL;
  numberExtraChangedBounds_ = 0;
#endif
}
// Does part of constructor
void 
CbcIntegerBranchingObject::fillPart (int variable,
				 int way , double value) 
{
  //originalObject_=NULL;
  branchIndex_=0;
  value_=value;
  numberBranches_=2;
  //model_= model;
  //originalCbcObject_=NULL;
  variable_=variable;
  way_=way;
  int iColumn = variable;
  down_[0] = model_->solver()->getColLower()[iColumn];
  down_[1] = floor(value_);
  up_[0] = ceil(value_);
  up_[1] = model_->getColUpper()[iColumn];
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
#ifdef FUNNY_BRANCHING
  variables_ = NULL;
  newBounds_ = NULL;
  numberExtraChangedBounds_ = 0;
#endif
}
  

// Copy constructor 
CbcIntegerBranchingObject::CbcIntegerBranchingObject ( const CbcIntegerBranchingObject & rhs) :CbcBranchingObject(rhs)
{
  down_[0] = rhs.down_[0];
  down_[1] = rhs.down_[1];
  up_[0] = rhs.up_[0];
  up_[1] = rhs.up_[1];
#ifdef FUNNY_BRANCHING
  numberExtraChangedBounds_ = rhs.numberExtraChangedBounds_;
  int size = numberExtraChangedBounds_*(sizeof(double)+sizeof(int));
  char * temp = new char [size];
  newBounds_ = (double *) temp;
  variables_ = (int *) (newBounds_+numberExtraChangedBounds_);

  int i ;
  for (i=0;i<numberExtraChangedBounds_;i++) {
    variables_[i]=rhs.variables_[i];
    newBounds_[i]=rhs.newBounds_[i];
  }
#endif
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
#ifdef FUNNY_BRANCHING
    delete [] newBounds_;
    numberExtraChangedBounds_ = rhs.numberExtraChangedBounds_;
    int size = numberExtraChangedBounds_*(sizeof(double)+sizeof(int));
    char * temp = new char [size];
    newBounds_ = (double *) temp;
    variables_ = (int *) (newBounds_+numberExtraChangedBounds_);
    
    int i ;
    for (i=0;i<numberExtraChangedBounds_;i++) {
      variables_[i]=rhs.variables_[i];
      newBounds_[i]=rhs.newBounds_[i];
    }
#endif
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
#ifdef FUNNY_BRANCHING
  delete [] newBounds_;
#endif
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
  if (down_[1]==-COIN_DBL_MAX)
    return 0.0;
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
    //#define CBC_PRINT2
#ifdef CBC_PRINT2
    printf("%d branching down has bounds %g %g",iColumn,down_[0],down_[1]);
#endif
#ifdef FUNNY_BRANCHING
    // branch - do extra bounds
    for (int i=0;i<numberExtraChangedBounds_;i++) {
      int variable = variables_[i];
      if ((variable&0x40000000)!=0) {
	// for going down
	int k = variable&0x3fffffff;
	assert (k!=iColumn);
	if ((variable&0x80000000)==0) {
	  // lower bound changing
#ifdef CBC_PRINT2
	  printf(" extra for %d changes lower from %g to %g",
		 k,model_->solver()->getColLower()[k],newBounds_[i]);
#endif
	  model_->solver()->setColLower(k,newBounds_[i]);
	} else {
	  // upper bound changing
#ifdef CBC_PRINT2
	  printf(" extra for %d changes upper from %g to %g",
		 k,model_->solver()->getColUpper()[k],newBounds_[i]);
#endif
	  model_->solver()->setColUpper(k,newBounds_[i]);
	}
      }
    }
#endif
#ifdef CBC_PRINT2
    printf("\n");
#endif
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
#ifdef CBC_PRINT2
    printf("%d branching up has bounds %g %g",iColumn,up_[0],up_[1]);
#endif
#ifdef FUNNY_BRANCHING
    // branch - do extra bounds
    for (int i=0;i<numberExtraChangedBounds_;i++) {
      int variable = variables_[i];
      if ((variable&0x40000000)==0) {
	// for going up
	int k = variable&0x3fffffff;
	assert (k!=iColumn);
	if ((variable&0x80000000)==0) {
	  // lower bound changing
#ifdef CBC_PRINT2
	  printf(" extra for %d changes lower from %g to %g",
		 k,model_->solver()->getColLower()[k],newBounds_[i]);
#endif
	  model_->solver()->setColLower(k,newBounds_[i]);
	} else {
	  // upper bound changing
#ifdef CBC_PRINT2
	  printf(" extra for %d changes upper from %g to %g",
		 k,model_->solver()->getColUpper()[k],newBounds_[i]);
#endif
	  model_->solver()->setColUpper(k,newBounds_[i]);
	}
      }
    }
#endif
#ifdef CBC_PRINT2
    printf("\n");
#endif
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
  if (nlb<olb+1.0e-8&&nub>oub-1.0e-8&&false)
    printf("bad null change for column %d - bounds %g,%g\n",iColumn,olb,oub);
#endif
  return 0.0;
}
/* Update bounds in solver as in 'branch' and update given bounds.
   branchState is -1 for 'down' +1 for 'up' */
void 
CbcIntegerBranchingObject::fix(OsiSolverInterface * solver,
			       double * lower, double * upper,
			       int branchState) const 
{
  int iColumn = originalCbcObject_->columnNumber();
  assert (variable_==iColumn);
  if (branchState<0) {
    model_->solver()->setColLower(iColumn,down_[0]);
    lower[iColumn]=down_[0];
    model_->solver()->setColUpper(iColumn,down_[1]);
    upper[iColumn]=down_[1];
  } else {
    model_->solver()->setColLower(iColumn,up_[0]);
    lower[iColumn]=up_[0];
    model_->solver()->setColUpper(iColumn,up_[1]);
    upper[iColumn]=up_[1];
  }
}
#ifdef FUNNY_BRANCHING
// Deactivate bounds for branching
void 
CbcIntegerBranchingObject::deactivate()
{
  down_[1]=-COIN_DBL_MAX;
}
int
CbcIntegerBranchingObject::applyExtraBounds(int iColumn, double lower, double upper, int way)
{
  // branch - do bounds

  int i;
  int found=0;
  if (variable_==iColumn) {
    printf("odd applyExtra %d\n",iColumn);
    if (way<0) {
      down_[0]=CoinMax(lower,down_[0]);
      down_[1]=CoinMin(upper,down_[1]);
      assert (down_[0]<=down_[1]);
    } else {
      up_[0]=CoinMax(lower,up_[0]);
      up_[1]=CoinMin(upper,up_[1]);
      assert (up_[0]<=up_[1]);
    }
    return 0;
  }
  int check = (way<0) ? 0x40000000 : 0;
  double newLower=lower;
  double newUpper=upper;
  for (i=0;i<numberExtraChangedBounds_;i++) {
    int variable = variables_[i];
    if ((variable&0x40000000)==check) {
      int k = variable&0x3fffffff;
      if (k==iColumn) {
	if ((variable&0x80000000)==0) {
	  // lower bound changing
	  found |= 1;
	  newBounds_[i] = CoinMax(lower,newBounds_[i]);
	  newLower = newBounds_[i];
	} else {
	  // upper bound changing
	  found |= 2;
	  newBounds_[i] = CoinMin(upper,newBounds_[i]);
	  newUpper = newBounds_[i];
	}
      }
    }
  }
  int nAdd=0;
  if ((found&2)==0) {
    // need to add new upper
    nAdd++;
  }
  if ((found&1)==0) {
    // need to add new lower
    nAdd++;
  }
  if (nAdd) { 
    int size = (numberExtraChangedBounds_+nAdd)*(sizeof(double)+sizeof(int));
    char * temp = new char [size];
    double * newBounds = (double *) temp;
    int * variables = (int *) (newBounds+numberExtraChangedBounds_+nAdd);

    int i ;
    for (i=0;i<numberExtraChangedBounds_;i++) {
      variables[i]=variables_[i];
      newBounds[i]=newBounds_[i];
    }
    delete [] newBounds_;
    newBounds_ = newBounds;
    variables_ = variables;
    if ((found&2)==0) {
      // need to add new upper
      int variable = iColumn | 0x80000000;
      variables_[numberExtraChangedBounds_]=variable;
      newBounds_[numberExtraChangedBounds_++]=newUpper;
    }
    if ((found&1)==0) {
      // need to add new lower
      int variable = iColumn;
      variables_[numberExtraChangedBounds_]=variable;
      newBounds_[numberExtraChangedBounds_++]=newLower;
    }
  }
  
  return (newUpper>=newLower) ? 0 : 1;
}
#endif
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

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
   */
CbcRangeCompare
CbcIntegerBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
  const CbcIntegerBranchingObject* br =
    dynamic_cast<const CbcIntegerBranchingObject*>(brObj);
  assert(br);
  double* thisBd = way_ < 0 ? down_ : up_;
  const double* otherBd = br->way_ < 0 ? br->down_ : br->up_;
  return CbcCompareRanges(thisBd, otherBd, replaceIfOverlap);
}

//##############################################################################

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
  *this=CbcSimpleIntegerPseudoCost(model,iColumn,downPseudoCost,upPseudoCost);
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

//##############################################################################

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

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcIntegerPseudoCostBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
  const CbcIntegerPseudoCostBranchingObject* br =
    dynamic_cast<const CbcIntegerPseudoCostBranchingObject*>(brObj);
  assert(br);
  double* thisBd = way_ < 0 ? down_ : up_;
  const double* otherBd = br->way_ < 0 ? br->down_ : br->up_;
  return CbcCompareRanges(thisBd, otherBd, replaceIfOverlap);
}


//##############################################################################

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

static inline int
CbcCompareCliques(const CbcClique* cl0, const CbcClique* cl1)
{
  if (cl0->cliqueType() < cl1->cliqueType()) {
    return -1;
  }
  if (cl0->cliqueType() > cl1->cliqueType()) {
    return 1;
  }
  if (cl0->numberMembers() != cl1->numberMembers()) {
    return cl0->numberMembers() - cl1->numberMembers();
  }
  if (cl0->numberNonSOSMembers() != cl1->numberNonSOSMembers()) {
    return cl0->numberNonSOSMembers() - cl1->numberNonSOSMembers();
  }
  return memcmp(cl0->members(), cl1->members(),
		cl0->numberMembers() * sizeof(int));
}

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type. 
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcCliqueBranchingObject::compareOriginalObject
(const CbcBranchingObject* brObj) const
{
  const CbcCliqueBranchingObject* br =
    dynamic_cast<const CbcCliqueBranchingObject*>(brObj);
  assert(br);
  return CbcCompareCliques(clique_, br->clique_);
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcCliqueBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
  const CbcCliqueBranchingObject* br =
    dynamic_cast<const CbcCliqueBranchingObject*>(brObj);
  assert(br);
  unsigned int* thisMask = way_ < 0 ? upMask_ : downMask_;
  const unsigned int* otherMask = br->way_ < 0 ? br->upMask_ : br->downMask_;
  const CoinUInt64 cl0 =
    (static_cast<CoinUInt64>(thisMask[0]) << 32) | thisMask[1];
  const CoinUInt64 cl1 =
    (static_cast<CoinUInt64>(otherMask[0]) << 32) | otherMask[1];
  if (cl0 == cl1) {
    return CbcRangeSame;
  }
  const CoinUInt64 cl_intersection = (cl0 & cl1);
  if (cl_intersection == cl0) {
    return CbcRangeSuperset;
  }
  if (cl_intersection == cl1) {
    return CbcRangeSubset;
  }
  const CoinUInt64 cl_xor = (cl0 ^ cl1);
  if (cl_intersection == 0 && cl_xor == 0) {
    return CbcRangeDisjoint;
  }
  const CoinUInt64 cl_union = (cl0 | cl1);
  thisMask[0] = static_cast<unsigned int>(cl_union >> 32);
  thisMask[1] = static_cast<unsigned int>(cl_union & 0xffffffff);
  return CbcRangeOverlap;
}

//##############################################################################

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

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type. 
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcLongCliqueBranchingObject::compareOriginalObject
(const CbcBranchingObject* brObj) const
{
  const CbcLongCliqueBranchingObject* br =
    dynamic_cast<const CbcLongCliqueBranchingObject*>(brObj);
  assert(br);
  return CbcCompareCliques(clique_, br->clique_);
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcLongCliqueBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
  const CbcLongCliqueBranchingObject* br =
    dynamic_cast<const CbcLongCliqueBranchingObject*>(brObj);
  assert(br);
  const int numberMembers = clique_->numberMembers();
  const int numberWords=(numberMembers+31)>>5;
  unsigned int* thisMask = way_ < 0 ? upMask_ : downMask_;
  const unsigned int* otherMask = br->way_ < 0 ? br->upMask_ : br->downMask_;

  if (memcmp(thisMask, otherMask, numberWords * sizeof(unsigned int)) == 0) {
    return CbcRangeSame;
  }
  bool canBeSuperset = true;
  bool canBeSubset = true;
  int i;
  for (i = numberWords-1; i >= 0 && (canBeSuperset || canBeSubset); --i) {
    const unsigned int both = (thisMask[i] & otherMask[i]);
    canBeSuperset &= (both == thisMask[i]);
    canBeSubset &= (both == otherMask[i]);
  }
  if (canBeSuperset) {
    return CbcRangeSuperset;
  }
  if (canBeSubset) {
    return CbcRangeSubset;
  }

  for (i = numberWords-1; i >= 0; --i) {
    if ((thisMask[i] ^ otherMask[i]) != 0) {
      break;
    }
  }
  if (i == -1) { // complement
    return CbcRangeDisjoint;
  }
  // must be overlap
  for (i = numberWords-1; i >= 0; --i) {
    thisMask[i] |= otherMask[i];
  }
  return CbcRangeOverlap;
}

//##############################################################################

// Default Constructor 
CbcSOSBranchingObject::CbcSOSBranchingObject()
  :CbcBranchingObject(),
   firstNonzero_(-1),
   lastNonzero_(-1)
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
  computeNonzeroRange();
}

// Copy constructor 
CbcSOSBranchingObject::CbcSOSBranchingObject (const CbcSOSBranchingObject & rhs)
 :CbcBranchingObject(rhs),
  firstNonzero_(rhs.firstNonzero_),
  lastNonzero_(rhs.lastNonzero_)
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
    firstNonzero_ = rhs.firstNonzero_;
    lastNonzero_ = rhs.lastNonzero_;
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

void
CbcSOSBranchingObject::computeNonzeroRange()
{
  const int numberMembers = set_->numberMembers();
  const double * weights = set_->weights();
  int i = 0;
  if (way_ < 0) {
    for ( i=0;i<numberMembers;i++) {
      if (weights[i] > separator_)
	break;
    }
    assert (i<numberMembers);
    firstNonzero_ = 0;
    lastNonzero_ = i;
  } else {
    for ( i=0;i<numberMembers;i++) {
      if (weights[i] >= separator_)
	break;
    }
    assert (i<numberMembers);
    firstNonzero_ = i;
    lastNonzero_ = numberMembers;
  }
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
  computeNonzeroRange();
  return 0.0;
}
/* Update bounds in solver as in 'branch' and update given bounds.
   branchState is -1 for 'down' +1 for 'up' */
void 
CbcSOSBranchingObject::fix(OsiSolverInterface * solver,
			       double * lower, double * upper,
			       int branchState) const 
{
  int numberMembers = set_->numberMembers();
  const int * which = set_->members();
  const double * weights = set_->weights();
  //const double * lower = solver->getColLower();
  //const double * upper = solver->getColUpper();
  // *** for way - up means fix all those in down section
  if (branchState<0) {
    int i;
    for ( i=0;i<numberMembers;i++) {
      if (weights[i] > separator_)
	break;
    }
    assert (i<numberMembers);
    for (;i<numberMembers;i++) {
      solver->setColUpper(which[i],0.0);
      upper[which[i]]=0.0;
    }
  } else {
    int i;
    for ( i=0;i<numberMembers;i++) {
      if (weights[i] >= separator_) {
	break;
      } else {
	solver->setColUpper(which[i],0.0);
	upper[which[i]]=0.0;
      }
    }
    assert (i<numberMembers);
  }
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
  
/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type. 
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcSOSBranchingObject::compareOriginalObject
(const CbcBranchingObject* brObj) const
{
  const CbcSOSBranchingObject* br =
    dynamic_cast<const CbcSOSBranchingObject*>(brObj);
  assert(br);
  const CbcSOS* s0 = set_;
  const CbcSOS* s1 = br->set_;
  if (s0->sosType() != s1->sosType()) {
    return s0->sosType() - s1->sosType();
  }
  if (s0->numberMembers() != s1->numberMembers()) {
    return s0->numberMembers() - s1->numberMembers();
  }
  const int memberCmp = memcmp(s0->members(), s1->members(),
			       s0->numberMembers() * sizeof(int));
  if (memberCmp != 0) {
    return memberCmp;
  }
  return memcmp(s0->weights(), s1->weights(),
		s0->numberMembers() * sizeof(double));
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcSOSBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
  const CbcSOSBranchingObject* br =
    dynamic_cast<const CbcSOSBranchingObject*>(brObj);
  assert(br);
  if (firstNonzero_ < br->firstNonzero_) {
    if (lastNonzero_ >= br->lastNonzero_) {
      return CbcRangeSuperset;
    } else if (lastNonzero_ <= br->firstNonzero_) {
      return CbcRangeDisjoint;
    } else {
      // overlap
      if (replaceIfOverlap) {
	firstNonzero_ = br->firstNonzero_;
      }
      return CbcRangeOverlap;
    }
  } else if (firstNonzero_ > br->firstNonzero_) {
    if (lastNonzero_ <= br->lastNonzero_) {
      return CbcRangeSubset;
    } else if (firstNonzero_ >= br->lastNonzero_) {
      return CbcRangeDisjoint;
    } else {
      // overlap
      if (replaceIfOverlap) {
	lastNonzero_ = br->lastNonzero_;
      }
      return CbcRangeOverlap;
    }
  } else {
    if (lastNonzero_ == br->lastNonzero_) {
      return CbcRangeSame;
    }
    return lastNonzero_ < br->lastNonzero_ ? CbcRangeSubset : CbcRangeSuperset;
  }
  return CbcRangeSame; // fake return
}

//##############################################################################

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
	if (CoinMin(changeUp,changeDn)<bestCriterion_)
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
	int thisNumber = CoinMin(numberInfeasibilitiesUp[i],numberInfeasibilitiesDown[i]);
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
	      if (CoinMin(changeUp[i],changeDown[i])<bestCriterion)
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
	    bestCriterion = CoinMin(changeUp[i],changeDown[i]);
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
	  bestCriterion = CoinMin(changeUp[i],changeDown[i]);
	  whichObject = i;
	  bestWay = betterWay;
	}
      }
      break;
    case 2:
      for ( i = 0 ; i < numberObjects ; i++) {
	double change = CoinMin(changeUp[i],changeDown[i]);
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
	int thisNumber = CoinMin(numberInfeasibilitiesUp[i],numberInfeasibilitiesDown[i]);
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
	      if (CoinMin(changeUp[i],changeDown[i])<bestCriterion)
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
	    bestCriterion = CoinMin(changeUp[i],changeDown[i]);
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

//##############################################################################

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

//##############################################################################

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

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type. 
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcFixingBranchingObject::compareOriginalObject
(const CbcBranchingObject* brObj) const
{
  throw("must implement");
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
   */
CbcRangeCompare
CbcFixingBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
#ifndef NDEBUG
  const CbcFixingBranchingObject* br =
    dynamic_cast<const CbcFixingBranchingObject*>(brObj);
  assert(br);
#endif
  // If two FixingBranchingObject's have the same base object then it's pretty
  // much guaranteed
  throw("must implement");
}

//##############################################################################

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
  
//##############################################################################

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

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type. 
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcNWayBranchingObject::compareOriginalObject
(const CbcBranchingObject* brObj) const
{
  throw("must implement");
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcNWayBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
  throw("must implement");
}

//##############################################################################

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

//##############################################################################

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

// Default Constructor
CbcGeneral::CbcGeneral() 
  : CbcObject()
{
}

// Constructor from model
CbcGeneral::CbcGeneral(CbcModel * model)
  : CbcObject(model)
{
}


// Destructor 
CbcGeneral::~CbcGeneral ()
{
}

// Copy constructor 
CbcGeneral::CbcGeneral ( const CbcGeneral & rhs)
  : CbcObject(rhs)
{
}
#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"
#include "ClpNode.hpp"
#include "CbcBranchDynamic.hpp"
// Assignment operator 
CbcGeneral & 
CbcGeneral::operator=( const CbcGeneral& rhs)
{
  if (this!=&rhs) {
    CbcObject::operator=(rhs);
  }
  return *this;
}

// Default Constructor 
CbcGeneralDepth::CbcGeneralDepth ()
  : CbcGeneral(),
    maximumDepth_(0),
    whichSolution_(-1),
    numberNodes_(0),
    nodeInfo_(NULL)
{
}

// Useful constructor (which are integer indices)
CbcGeneralDepth::CbcGeneralDepth (CbcModel * model, int maximumDepth)
  : CbcGeneral(model),
    maximumDepth_(maximumDepth),
    whichSolution_(-1),
    numberNodes_(0),
    nodeInfo_(NULL)
{
  int nNodes = maximumNodes();
  if (nNodes) {
    nodeInfo_ = new ClpNodeStuff();
    ClpNodeStuff * info = nodeInfo_;
    // for reduced costs and duals
    info->solverOptions_ |= 7;
    if (maximumDepth_>0) {
      info->nDepth_ = maximumDepth_;
    } else {
      info->nDepth_ = - maximumDepth_;
      info->solverOptions_ |= 32;
    }
    ClpNode ** nodeInfo = new ClpNode * [nNodes];
    for (int i=0;i<nNodes;i++) 
      nodeInfo[i]=NULL;
    info->nodeInfo_ = nodeInfo;
  } else {
    nodeInfo_ = NULL;
  }
}

// Copy constructor 
CbcGeneralDepth::CbcGeneralDepth ( const CbcGeneralDepth & rhs)
  :CbcGeneral(rhs)
{
  maximumDepth_ = rhs.maximumDepth_;
  whichSolution_ = -1;
  numberNodes_ = 0;
  int nNodes = maximumNodes();
  if (nNodes) {
    assert (rhs.nodeInfo_);
    nodeInfo_ = new ClpNodeStuff(*rhs.nodeInfo_);
    ClpNodeStuff * info = nodeInfo_;
    if (maximumDepth_>0) {
      info->nDepth_ = maximumDepth_;
    } else {
      info->nDepth_ = - maximumDepth_;
      info->solverOptions_ |= 32;
    }
    if (!info->nodeInfo_) {
      ClpNode ** nodeInfo = new ClpNode * [nNodes];
      for (int i=0;i<nNodes;i++) 
	nodeInfo[i]=NULL;
      info->nodeInfo_ = nodeInfo;
    }
  } else {
    nodeInfo_ = NULL;
  }
}

// Clone
CbcObject *
CbcGeneralDepth::clone() const
{
  return new CbcGeneralDepth(*this);
}

// Assignment operator 
CbcGeneralDepth & 
CbcGeneralDepth::operator=( const CbcGeneralDepth& rhs)
{
  if (this!=&rhs) {
    CbcGeneral::operator=(rhs);
    delete nodeInfo_;
    maximumDepth_ = rhs.maximumDepth_;
    whichSolution_ = -1;
    numberNodes_ = 0;
    if (maximumDepth_) {
      assert (rhs.nodeInfo_);
      nodeInfo_ = new ClpNodeStuff(*rhs.nodeInfo_);
    } else {
      nodeInfo_ = NULL;
    }
  }
  return *this;
}

// Destructor 
CbcGeneralDepth::~CbcGeneralDepth ()
{
  delete nodeInfo_;
}
// Return maximum number of nodes
int 
CbcGeneralDepth::maximumNodes() const
{
  int n;
  if (maximumDepth_>0) 
    n = (1<<maximumDepth_)+1+maximumDepth_;
  else if (maximumDepth_<0)
    n = 1+1-maximumDepth_;
  else
    n = 0;
  return n;
}

// Infeasibility - large is 0.5
double 
CbcGeneralDepth::infeasibility(int & preferredWay) const
{
  whichSolution_ = -1;
  // should use genuine OsiBranchingInformation usefulInfo = model_->usefulInformation();
  // for now assume only called when correct
  //if (usefulInfo.depth_>=4&&!model_->parentModel()
  //     &&(usefulInfo.depth_%2)==0) {
  if (true) {
    OsiSolverInterface * solver = model_->solver();
    OsiClpSolverInterface * clpSolver 
      = dynamic_cast<OsiClpSolverInterface *> (solver);
    if (clpSolver) {
      ClpNodeStuff * info = nodeInfo_;
      info->integerTolerance_=model_->getIntegerTolerance();
      info->integerIncrement_=model_->getCutoffIncrement();
      int numberIntegers = model_->numberIntegers();
      double * down = new double[numberIntegers];
      double * up = new double[numberIntegers];
      int * numberDown = new int[numberIntegers];
      int * numberUp = new int[numberIntegers];
      int * numberDownInfeasible = new int[numberIntegers];
      int * numberUpInfeasible = new int[numberIntegers];
      model_->fillPseudoCosts(down,up,numberDown,numberUp,
		      numberDownInfeasible,numberUpInfeasible);
      info->fillPseudoCosts(down,up,numberDown,numberUp,
			    numberDownInfeasible,
			    numberUpInfeasible,numberIntegers);
      info->presolveType_= 1; 
      delete [] down;
      delete [] up;
      delete [] numberDown;
      delete [] numberUp;
      delete [] numberDownInfeasible;
      delete [] numberUpInfeasible;
      bool takeHint;
      OsiHintStrength strength;
      solver->getHintParam(OsiDoReducePrint,takeHint,strength);
      ClpSimplex * simplex = clpSolver->getModelPtr();
      int saveLevel = simplex->logLevel();
      if (strength!=OsiHintIgnore&&takeHint&&saveLevel==1)
	simplex->setLogLevel(0);
      clpSolver->setBasis();
      whichSolution_ = simplex->fathomMany(info);
      model_->incrementExtra(info->numberNodesExplored_,
			     info->numberIterations_);
      // update pseudo costs
      double smallest=1.0e50;
      double largest=-1.0;
      OsiObject ** objects = model_->objects();
#ifndef NDEBUG
      const int * integerVariable = model_->integerVariable();
#endif
      for (int i=0;i<numberIntegers;i++) {
	CbcSimpleIntegerDynamicPseudoCost * obj =
	  dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(objects[i]) ;
	assert (obj&&obj->columnNumber()==integerVariable[i]);
	if (info->numberUp_[i]>0) {
	  if (info->downPseudo_[i]>largest)
	    largest=info->downPseudo_[i];
	  if (info->downPseudo_[i]<smallest) 
	    smallest=info->downPseudo_[i];
	  if (info->upPseudo_[i]>largest) 
	    largest=info->upPseudo_[i];
	  if (info->upPseudo_[i]<smallest)
	    smallest=info->upPseudo_[i];
	  obj->updateAfterMini(info->numberDown_[i],
			       info->numberDownInfeasible_[i],
			       info->downPseudo_[i],
			       info->numberUp_[i],
			       info->numberUpInfeasible_[i],
			       info->upPseudo_[i]);
	}
      }
      //printf("range of costs %g to %g\n",smallest,largest);
      simplex->setLogLevel(saveLevel);
      numberNodes_ = info->nNodes_;
      int numberDo = numberNodes_;
      if (whichSolution_>=0)
	numberDo--;
      if (numberDo>0) {
	return 0.5;
      } else {
	// no solution
	return COIN_DBL_MAX; // say infeasible
      }
    } else {
      return -1.0;
    }
  } else {
    return -1.0;
  }
}

// This looks at solution and sets bounds to contain solution
void 
CbcGeneralDepth::feasibleRegion()
{
  // Other stuff should have done this
}
// Redoes data when sequence numbers change
void 
CbcGeneralDepth::redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns)
{
}

//#define CHECK_PATH
#ifdef CHECK_PATH
extern const double * debuggerSolution_Z;
extern int numberColumns_Z;
extern int gotGoodNode_Z;
#endif

// Creates a branching object
CbcBranchingObject * 
CbcGeneralDepth::createBranch(int way) 
{
  int numberDo = numberNodes_;
  if (whichSolution_>=0)
    numberDo--;
  assert (numberDo>0);
  // create object
  CbcGeneralBranchingObject * branch = new CbcGeneralBranchingObject(model_);
  // skip solution
  branch->numberSubProblems_ = numberDo;
  // If parentBranch_ back in then will have to be 2*
  branch->numberSubLeft_ = numberDo;
  branch->setNumberBranches(numberDo);
  CbcSubProblem * sub = new CbcSubProblem[numberDo];
  int iProb=0;
  branch->subProblems_ = sub;
  branch->numberRows_ = model_->solver()->getNumRows();
  int iNode;
  OsiSolverInterface * solver = model_->solver();
  OsiClpSolverInterface * clpSolver 
    = dynamic_cast<OsiClpSolverInterface *> (solver);
  assert (clpSolver);
  ClpSimplex * simplex = clpSolver->getModelPtr();
  int numberColumns = simplex->numberColumns();
  double * lowerBefore=CoinCopyOfArray(simplex->getColLower(),
				       numberColumns);
  double * upperBefore=CoinCopyOfArray(simplex->getColUpper(),
				       numberColumns);
  ClpNodeStuff * info = nodeInfo_;
  double * weight = new double[numberNodes_];
  int * whichNode = new int [numberNodes_];
  // Sort 
  for (iNode=0;iNode<numberNodes_;iNode++) {
    if (iNode!=whichSolution_) {
      double objectiveValue = info->nodeInfo_[iNode]->objectiveValue();
      double sumInfeasibilities = info->nodeInfo_[iNode]->sumInfeasibilities();
      int numberInfeasibilities = info->nodeInfo_[iNode]->numberInfeasibilities();
      double thisWeight = 0.0;
#if 1
      // just closest
      thisWeight = 1.0e9*numberInfeasibilities;
      thisWeight += sumInfeasibilities;
      thisWeight += 1.0e-7*objectiveValue;
      // Try estimate
      thisWeight = info->nodeInfo_[iNode]->estimatedSolution();
#else
      thisWeight = 1.0e-3*numberInfeasibilities;
      thisWeight += 1.0e-5*sumInfeasibilities;
      thisWeight += objectiveValue;
#endif
      whichNode[iProb]=iNode;
      weight[iProb++]=thisWeight;
    }
  }
  assert (iProb==numberDo);
  CoinSort_2(weight,weight+numberDo,whichNode);
  for (iProb=0;iProb<numberDo;iProb++) {
    iNode = whichNode[iProb];
    ClpNode * node = info->nodeInfo_[iNode];
    // move bounds
    node->applyNode(simplex,3);
    // create subproblem
    sub[iProb]=CbcSubProblem(clpSolver,lowerBefore,upperBefore,
			     node->statusArray(),node->depth());
    sub[iProb].objectiveValue_ = node->objectiveValue();
    sub[iProb].sumInfeasibilities_ = node->sumInfeasibilities();
    sub[iProb].numberInfeasibilities_ = node->numberInfeasibilities();
#ifdef CHECK_PATH
    if (simplex->numberColumns()==numberColumns_Z) {
      bool onOptimal=true;
      const double * columnLower = simplex->columnLower();
      const double * columnUpper = simplex->columnUpper();
      for (int i=0;i<numberColumns_Z;i++) {
	if (iNode==gotGoodNode_Z)
	  printf("good %d %d %g %g\n",iNode,i,columnLower[i],columnUpper[i]);
	if (columnUpper[i]<debuggerSolution_Z[i]||columnLower[i]>debuggerSolution_Z[i]&&simplex->isInteger(i)) {
	  onOptimal=false;
	  break;
	}
      }
      if (onOptimal) {
	printf("adding to node %x as %d - objs\n",this,iProb);
	for (int j=0;j<=iProb;j++)
	  printf("%d %g\n",j,sub[j].objectiveValue_);
      }
    }
#endif
  }
  delete [] weight;
  delete [] whichNode;
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  // restore bounds
  for ( int j=0;j<numberColumns;j++) {
    if (lowerBefore[j] != lower[j])
      solver->setColLower(j,lowerBefore[j]);
    if (upperBefore[j] != upper[j])
      solver->setColUpper(j,upperBefore[j]);
  }
  delete [] upperBefore;
  delete [] lowerBefore;
  return branch;
}

// Default Constructor 
CbcGeneralBranchingObject::CbcGeneralBranchingObject()
  :CbcBranchingObject(),
   subProblems_(NULL),
   node_(NULL),
   numberSubProblems_(0),
   numberSubLeft_(0),
   whichNode_(-1),
   numberRows_(0)
{
  //  printf("CbcGeneral %x default constructor\n",this);
}

// Useful constructor
CbcGeneralBranchingObject::CbcGeneralBranchingObject (CbcModel * model)
  :CbcBranchingObject(model,-1,-1,0.5),
   subProblems_(NULL),
   node_(NULL),
   numberSubProblems_(0),
   numberSubLeft_(0),
   whichNode_(-1),
   numberRows_(0)
{
  //printf("CbcGeneral %x useful constructor\n",this);
}

// Copy constructor 
CbcGeneralBranchingObject::CbcGeneralBranchingObject ( const CbcGeneralBranchingObject & rhs) 
  :CbcBranchingObject(rhs),
   subProblems_(NULL),
   node_(rhs.node_),
   numberSubProblems_(rhs.numberSubProblems_),
   numberSubLeft_(rhs.numberSubLeft_),
   whichNode_(rhs.whichNode_),
   numberRows_(rhs.numberRows_)
{
  abort();
  if (numberSubProblems_) {
    subProblems_ = new CbcSubProblem[numberSubProblems_];
    for (int i=0;i<numberSubProblems_;i++)
      subProblems_[i]=rhs.subProblems_[i];
  }
}

// Assignment operator 
CbcGeneralBranchingObject & 
CbcGeneralBranchingObject::operator=( const CbcGeneralBranchingObject& rhs)
{
  if (this != &rhs) {
    abort();
    CbcBranchingObject::operator=(rhs);
    delete [] subProblems_;
    numberSubProblems_ = rhs.numberSubProblems_;
    numberSubLeft_ = rhs.numberSubLeft_;
    whichNode_ = rhs.whichNode_;
    numberRows_ = rhs.numberRows_;
    if (numberSubProblems_) {
      subProblems_ = new CbcSubProblem[numberSubProblems_];
      for (int i=0;i<numberSubProblems_;i++)
	subProblems_[i]=rhs.subProblems_[i];
    } else {
      subProblems_ = NULL;
    }
    node_ = rhs.node_;
  }
  return *this;
}
CbcBranchingObject * 
CbcGeneralBranchingObject::clone() const
{ 
  return (new CbcGeneralBranchingObject(*this));
}


// Destructor 
CbcGeneralBranchingObject::~CbcGeneralBranchingObject ()
{
  //printf("CbcGeneral %x destructor\n",this);
  delete [] subProblems_;
}
bool doingDoneBranch=false;
double
CbcGeneralBranchingObject::branch()
{
  double cutoff=model_->getCutoff();
  if (whichNode_<0) {
    assert (node_);
    bool applied=false;
    while (numberBranchesLeft()) {
      int which = branchIndex();
      decrementNumberBranchesLeft();
      CbcSubProblem * thisProb = subProblems_+which;
      if (thisProb->objectiveValue_<cutoff) {
	//printf("branch %x (sub %x) which now %d\n",this,
	//     subProblems_,which);
	OsiSolverInterface * solver = model_->solver();
	thisProb->apply(solver);
	OsiClpSolverInterface * clpSolver 
	  = dynamic_cast<OsiClpSolverInterface *> (solver);
	assert (clpSolver);
	// Move status to basis
	clpSolver->setWarmStart(NULL);
	//ClpSimplex * simplex = clpSolver->getModelPtr();
	node_->setObjectiveValue(thisProb->objectiveValue_);
	node_->setSumInfeasibilities(thisProb->sumInfeasibilities_);
	node_->setNumberUnsatisfied(thisProb->numberInfeasibilities_);
	applied=true;
	doingDoneBranch=true;
	break;
      } else if (numberBranchesLeft()) {
	node_->nodeInfo()->branchedOn() ;
      }
    }
    if (!applied) {
      // no good one
      node_->setObjectiveValue(cutoff+1.0e20);
      node_->setSumInfeasibilities(1.0);
      node_->setNumberUnsatisfied(1);
      assert (whichNode_<0);
    }
  } else {
    decrementNumberBranchesLeft();
    CbcSubProblem * thisProb = subProblems_+whichNode_;
    assert (thisProb->objectiveValue_<cutoff);
    OsiSolverInterface * solver = model_->solver();
    thisProb->apply(solver);
    OsiClpSolverInterface * clpSolver 
      = dynamic_cast<OsiClpSolverInterface *> (solver);
    assert (clpSolver);
    // Move status to basis
    clpSolver->setWarmStart(NULL);
  }
  return 0.0;
}
/* Double checks in case node can change its mind!
   Can change objective etc */
void 
CbcGeneralBranchingObject::checkIsCutoff(double cutoff)
{
  assert (node_);
  int first = branchIndex();
  int last = first + numberBranchesLeft();
  for (int which=first;which<last;which++) {
    CbcSubProblem * thisProb = subProblems_+which;
    if (thisProb->objectiveValue_<cutoff) {
      node_->setObjectiveValue(thisProb->objectiveValue_);
      node_->setSumInfeasibilities(thisProb->sumInfeasibilities_);
      node_->setNumberUnsatisfied(thisProb->numberInfeasibilities_);
      break;
    }
  }
}
// Print what would happen  
void
CbcGeneralBranchingObject::print()
{
  //printf("CbcGeneralObject has %d subproblems\n",numberSubProblems_);
}
// Fill in current objective etc
void 
CbcGeneralBranchingObject::state(double & objectiveValue,
				 double & sumInfeasibilities,
				 int & numberUnsatisfied,int which) const
{
  assert (which>=0&&which<numberSubProblems_);
  const CbcSubProblem * thisProb = subProblems_+which;
  objectiveValue = thisProb->objectiveValue_;
  sumInfeasibilities = thisProb->sumInfeasibilities_;
  numberUnsatisfied = thisProb->numberInfeasibilities_;
}
/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type. 
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcGeneralBranchingObject::compareOriginalObject
(const CbcBranchingObject* brObj) const
{
  throw("must implement");
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcGeneralBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
  throw("must implement");
}

// Default Constructor 
CbcOneGeneralBranchingObject::CbcOneGeneralBranchingObject()
  :CbcBranchingObject(),
   object_(NULL),
   whichOne_(-1)
{
  //printf("CbcOneGeneral %x default constructor\n",this);
}

// Useful constructor
CbcOneGeneralBranchingObject::CbcOneGeneralBranchingObject (CbcModel * model,
				 CbcGeneralBranchingObject * object,
				 int whichOne)
  :CbcBranchingObject(model,-1,-1,0.5),
   object_(object),
   whichOne_(whichOne)
{
  //printf("CbcOneGeneral %x useful constructor object %x %d left\n",this,
  //	 object_,object_->numberSubLeft_);
  numberBranches_=1;
}

// Copy constructor 
CbcOneGeneralBranchingObject::CbcOneGeneralBranchingObject ( const CbcOneGeneralBranchingObject & rhs) 
  :CbcBranchingObject(rhs),
   object_(rhs.object_),
   whichOne_(rhs.whichOne_)
{
}

// Assignment operator 
CbcOneGeneralBranchingObject & 
CbcOneGeneralBranchingObject::operator=( const CbcOneGeneralBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    object_ = rhs.object_;
    whichOne_ = rhs.whichOne_;
  }
  return *this;
}
CbcBranchingObject * 
CbcOneGeneralBranchingObject::clone() const
{ 
  return (new CbcOneGeneralBranchingObject(*this));
}


// Destructor 
CbcOneGeneralBranchingObject::~CbcOneGeneralBranchingObject ()
{
  //printf("CbcOneGeneral %x destructor object %x %d left\n",this,
  // object_,object_->numberSubLeft_);
  assert (object_->numberSubLeft_>0&&
	  object_->numberSubLeft_<1000000);
  if (!object_->decrementNumberLeft()) {
    // printf("CbcGeneral %x yy destructor\n",object_);
    delete object_;
  }
}
double
CbcOneGeneralBranchingObject::branch()
{
  assert (numberBranchesLeft());
  decrementNumberBranchesLeft();
  assert (!numberBranchesLeft());
  object_->setWhichNode(whichOne_);
  object_->branch();
  return 0.0;
}
/* Double checks in case node can change its mind!
   Can change objective etc */
void 
CbcOneGeneralBranchingObject::checkIsCutoff(double cutoff)
{
  assert (numberBranchesLeft());
}
// Print what would happen  
void
CbcOneGeneralBranchingObject::print()
{
  //printf("CbcOneGeneralObject has 1 subproblem\n");
}
/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type. 
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcOneGeneralBranchingObject::compareOriginalObject
(const CbcBranchingObject* brObj) const
{
  throw("must implement");
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcOneGeneralBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
  throw("must implement");
}
// Default Constructor 
CbcSubProblem::CbcSubProblem()
  : objectiveValue_(0.0),
    sumInfeasibilities_(0.0),
    variables_(NULL),
    newBounds_(NULL),
    status_(NULL),
    depth_(0),
    numberChangedBounds_(0),
    numberInfeasibilities_(0)
{
}

// Useful constructor
CbcSubProblem::CbcSubProblem (const OsiSolverInterface * solver,
			      const double * lastLower,
			      const double * lastUpper,
			      const unsigned char * status,
			      int depth)
  : objectiveValue_(0.0),
    sumInfeasibilities_(0.0),
    variables_(NULL),
    newBounds_(NULL),
    status_(NULL),
    depth_(depth),
    numberChangedBounds_(0),
    numberInfeasibilities_(0)
{
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();

  numberChangedBounds_=0;
  int numberColumns = solver->getNumCols();
  int i;
  for (i=0;i<numberColumns;i++) {
    if (lower[i]!=lastLower[i]) 
      numberChangedBounds_++;
    if (upper[i]!=lastUpper[i]) 
      numberChangedBounds_++;
  }
  if (numberChangedBounds_) {
    newBounds_ = new double [numberChangedBounds_] ;
    variables_ = new int [numberChangedBounds_] ;
    numberChangedBounds_=0;
    for (i=0;i<numberColumns;i++) {
      if (lower[i]!=lastLower[i]) {
	variables_[numberChangedBounds_]=i;
	newBounds_[numberChangedBounds_++]=lower[i];
      }
      if (upper[i]!=lastUpper[i]) {
	variables_[numberChangedBounds_]=i|0x80000000;
	newBounds_[numberChangedBounds_++]=upper[i];
      }
#ifdef CBC_DEBUG
      if (lower[i] != lastLower[i]) {
	std::cout
	  << "lower on " << i << " changed from "
	  << lastLower[i] << " to " << lower[i] << std::endl ;
      }
      if (upper[i] != lastUpper[i]) {
	std::cout
	  << "upper on " << i << " changed from "
	  << lastUpper[i] << " to " << upper[i] << std::endl ;
      }
#endif
    }
#ifdef CBC_DEBUG
    std::cout << numberChangedBounds_ << " changed bounds." << std::endl ;
#endif
  }
  const OsiClpSolverInterface * clpSolver 
    = dynamic_cast<const OsiClpSolverInterface *> (solver);
  assert (clpSolver);
  // Do difference
  // Current basis
  status_ = clpSolver->getBasis(status);
}

// Copy constructor 
CbcSubProblem::CbcSubProblem ( const CbcSubProblem & rhs) 
  : objectiveValue_(rhs.objectiveValue_),
    sumInfeasibilities_(rhs.sumInfeasibilities_),
    variables_(NULL),
    newBounds_(NULL),
    status_(NULL),
    depth_(rhs.depth_),
    numberChangedBounds_(rhs.numberChangedBounds_),
    numberInfeasibilities_(rhs.numberInfeasibilities_)
{
  if (numberChangedBounds_) {
    variables_ = CoinCopyOfArray(rhs.variables_,numberChangedBounds_);
    newBounds_ = CoinCopyOfArray(rhs.newBounds_,numberChangedBounds_);
  }
  if (rhs.status_) {
    status_ = new CoinWarmStartBasis(*rhs.status_);
  }
}

// Assignment operator 
CbcSubProblem & 
CbcSubProblem::operator=( const CbcSubProblem& rhs)
{
  if (this != &rhs) {
    delete [] variables_;
    delete [] newBounds_;
    delete status_;
    objectiveValue_ = rhs.objectiveValue_;
    sumInfeasibilities_ = rhs.sumInfeasibilities_;
    depth_ = rhs.depth_;
    numberChangedBounds_ = rhs.numberChangedBounds_;
    numberInfeasibilities_ = rhs.numberInfeasibilities_;
    if (numberChangedBounds_) {
      variables_ = CoinCopyOfArray(rhs.variables_,numberChangedBounds_);
      newBounds_ = CoinCopyOfArray(rhs.newBounds_,numberChangedBounds_);
    } else {
      variables_ = NULL;
      newBounds_ = NULL;
    }
    if (rhs.status_) {
      status_ = new CoinWarmStartBasis(*rhs.status_);
    } else {
      status_ = NULL;
    }
  }
  return *this;
}

// Destructor 
CbcSubProblem::~CbcSubProblem ()
{
  delete [] variables_;
  delete [] newBounds_;
  delete status_;
}
// Apply subproblem
void 
CbcSubProblem::apply(OsiSolverInterface * solver, int what) const
{
  int i;
  if ((what&1)!=0) {
    int nSame=0;
    for (i=0;i<numberChangedBounds_;i++) {
      int variable = variables_[i];
      int k = variable&0x3fffffff;
      if ((variable&0x80000000)==0) {
	// lower bound changing
	//#define CBC_PRINT2
#ifdef CBC_PRINT2
	if(solver->getColLower()[k]!=newBounds_[i])
	  printf("lower change for column %d - from %g to %g\n",k,solver->getColLower()[k],newBounds_[i]);
#endif
#ifndef NDEBUG
	if ((variable&0x40000000)==0&&true) {
	  double oldValue = solver->getColLower()[k];
	  assert (newBounds_[i]>oldValue-1.0e-8);
	  if (newBounds_[i]<oldValue+1.0e-8) {
#ifdef CBC_PRINT2
	    printf("bad null lower change for column %d - bound %g\n",k,oldValue);
#endif
	    if (newBounds_[i]==oldValue) 
	      nSame++;
	  }
	}
#endif
	solver->setColLower(k,newBounds_[i]);
      } else {
	// upper bound changing
#ifdef CBC_PRINT2
	if(solver->getColUpper()[k]!=newBounds_[i])
	  printf("upper change for column %d - from %g to %g\n",k,solver->getColUpper()[k],newBounds_[i]);
#endif
#ifndef NDEBUG
	if ((variable&0x40000000)==0&&true) {
	  double oldValue = solver->getColUpper()[k];
	  assert (newBounds_[i]<oldValue+1.0e-8);
	  if (newBounds_[i]>oldValue-1.0e-8) {
#ifdef CBC_PRINT2
	    printf("bad null upper change for column %d - bound %g\n",k,oldValue);
#endif
	    if (newBounds_[i]==oldValue) 
	      nSame++;
	  }
	}
#endif
	solver->setColUpper(k,newBounds_[i]);
      }
    }
    if (nSame&&(nSame<numberChangedBounds_||(what&3)!=3))
      printf("%d changes out of %d redundant %d\n",
	     nSame,numberChangedBounds_,what);
    else if (numberChangedBounds_&&what==7&&!nSame)
      printf("%d good changes %d\n",
	     numberChangedBounds_,what);
  }
#if 0
  if ((what&2)!=0) {
    OsiClpSolverInterface * clpSolver 
      = dynamic_cast<OsiClpSolverInterface *> (solver);
    assert (clpSolver);
    //assert (clpSolver->getNumRows()==numberRows_);
    //clpSolver->setBasis(*status_);
    // Current basis
    CoinWarmStartBasis * basis=clpSolver->getPointerToWarmStart();
    printf("BBBB\n");
    basis->print();
    assert (basis->fullBasis());
    basis->applyDiff(status_);
    printf("diff applied %x\n",status_);
    printf("CCCC\n");
    basis->print();
    assert (basis->fullBasis());
#ifndef NDEBUG
    if (!basis->fullBasis())
      printf("Debug this basis!!\n");
#endif
    clpSolver->setBasis(*basis);
  }
#endif
  if ((what&8)!=0) {
    OsiClpSolverInterface * clpSolver 
      = dynamic_cast<OsiClpSolverInterface *> (solver);
    assert (clpSolver);
    clpSolver->setBasis(*status_);
    delete status_;
    status_=NULL;
  }
}
#endif

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type. 
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcDummyBranchingObject::compareOriginalObject
(const CbcBranchingObject* brObj) const
{
  throw("must implement");
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcDummyBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool replaceIfOverlap)
{
  throw("must implement");
}
