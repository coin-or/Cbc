// Copyright (C) 2005, International Business Machines
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
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchLink.hpp"
#include "CoinError.hpp"

// Default Constructor 
CbcLink::CbcLink ()
  : CbcObject(),
    weights_(NULL),
    numberMembers_(0),
    numberLinks_(0),
    first_(-1)
{
}

// Useful constructor (which are indices)
CbcLink::CbcLink (CbcModel * model,  int numberMembers,
	   int numberLinks, int first , const double * weights, int identifier)
  : CbcObject(model),
    numberMembers_(numberMembers),
    numberLinks_(numberLinks),
    first_(first)
{
  id_=identifier;
  if (numberMembers_) {
    weights_ = new double[numberMembers_];
    if (weights) {
      memcpy(weights_,weights,numberMembers_*sizeof(double));
    } else {
      for (int i=0;i<numberMembers_;i++)
        weights_[i]=i;
    }
    // weights must be increasing
    int i;
    double last=-COIN_DBL_MAX;
    for (i=0;i<numberMembers_;i++) {
      assert (weights_[i]>last+1.0e-12);
      last=weights_[i];
    }
  } else {
    weights_ = NULL;
  }
}

// Copy constructor 
CbcLink::CbcLink ( const CbcLink & rhs)
  :CbcObject(rhs)
{
  numberMembers_ = rhs.numberMembers_;
  numberLinks_ = rhs.numberLinks_;
  first_ = rhs.first_;
  if (numberMembers_) {
    weights_ = new double[numberMembers_];
    memcpy(weights_,rhs.weights_,numberMembers_*sizeof(double));
  } else {
    weights_ = NULL;
  }
}

// Clone
CbcObject *
CbcLink::clone() const
{
  return new CbcLink(*this);
}

// Assignment operator 
CbcLink & 
CbcLink::operator=( const CbcLink& rhs)
{
  if (this!=&rhs) {
    CbcObject::operator=(rhs);
    delete [] weights_;
    numberMembers_ = rhs.numberMembers_;
    numberLinks_ = rhs.numberLinks_;
    first_ = rhs.first_;
    if (numberMembers_) {
      weights_ = new double[numberMembers_];
      memcpy(weights_,rhs.weights_,numberMembers_*sizeof(double));
    } else {
      weights_ = NULL;
    }
  }
  return *this;
}

// Destructor 
CbcLink::~CbcLink ()
{
  delete [] weights_;
}

// Infeasibility - large is 0.5
double 
CbcLink::infeasibility(int & preferredWay) const
{
  int j;
  int firstNonZero=-1;
  int lastNonZero = -1;
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double weight = 0.0;
  double sum =0.0;

  // check bounds etc
  double lastWeight=-1.0e100;
  int base=first_;
  for (j=0;j<numberMembers_;j++) {
    for (int k=0;k<numberLinks_;k++) {
      int iColumn = base+k;
      if (lower[iColumn])
        throw CoinError("Non zero lower bound in CBCLink","infeasibility","CbcLink");
      if (lastWeight>=weights_[j]-1.0e-7)
        throw CoinError("Weights too close together in CBCLink","infeasibility","CbcLink");
      double value = CoinMax(0.0,solution[iColumn]);
      sum += value;
      if (value>integerTolerance&&upper[iColumn]) {
        // Possibly due to scaling a fixed variable might slip through
        if (value>upper[iColumn]+1.0e-8) {
          // Could change to #ifdef CBC_DEBUG
#ifndef NDEBUG
          if (model_->messageHandler()->logLevel()>1)
            printf("** Variable %d (%d) has value %g and upper bound of %g\n",
                   iColumn,j,value,upper[iColumn]);
#endif
        } 
        weight += weights_[j]*CoinMin(value,upper[iColumn]);
        if (firstNonZero<0)
          firstNonZero=j;
        lastNonZero=j;
      }
    }
    base += numberLinks_;
  }
  preferredWay=1;
  if (lastNonZero-firstNonZero>=1) {
    // find where to branch
    assert (sum>0.0);
    weight /= sum;
    double value = lastNonZero-firstNonZero+1;
    value *= 0.5/((double) numberMembers_);
    return value;
  } else {
    return 0.0; // satisfied
  }
}

// This looks at solution and sets bounds to contain solution
void 
CbcLink::feasibleRegion()
{
  int j;
  int firstNonZero=-1;
  int lastNonZero = -1;
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * upper = solver->getColUpper();
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double weight = 0.0;
  double sum =0.0;

  int base=first_;
  for (j=0;j<numberMembers_;j++) {
    for (int k=0;k<numberLinks_;k++) {
      int iColumn = base+k;
      double value = CoinMax(0.0,solution[iColumn]);
      sum += value;
      if (value>integerTolerance&&upper[iColumn]) {
        weight += weights_[j]*value;
        if (firstNonZero<0)
          firstNonZero=j;
        lastNonZero=j;
      }
    }
    base += numberLinks_;
  }
  assert (lastNonZero-firstNonZero==0) ;
  base=first_;
  for (j=0;j<firstNonZero;j++) {
    for (int k=0;k<numberLinks_;k++) {
      int iColumn = base+k;
      solver->setColUpper(iColumn,0.0);
    }
    base += numberLinks_;
  }
  // skip
  base += numberLinks_;
  for (j=lastNonZero+1;j<numberMembers_;j++) {
    for (int k=0;k<numberLinks_;k++) {
      int iColumn = base+k;
      solver->setColUpper(iColumn,0.0);
    }
    base += numberLinks_;
  }
}


// Creates a branching object
CbcBranchingObject * 
CbcLink::createBranch(int way) 
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
  int base=first_;
  for (j=0;j<numberMembers_;j++) {
    for (int k=0;k<numberLinks_;k++) {
      int iColumn = base+k;
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
    base += numberLinks_;
  }
  assert (lastNonZero-firstNonZero>=1) ;
  // find where to branch
  assert (sum>0.0);
  weight /= sum;
  int iWhere;
  double separator=0.0;
  for (iWhere=firstNonZero;iWhere<lastNonZero;iWhere++) 
    if (weight<weights_[iWhere+1])
      break;
  separator = 0.5 *(weights_[iWhere]+weights_[iWhere+1]);
  // create object
  CbcBranchingObject * branch;
  branch = new CbcLinkBranchingObject(model_,this,way,separator);
  return branch;
}
// Useful constructor
CbcLinkBranchingObject::CbcLinkBranchingObject (CbcModel * model,
					      const CbcLink * set,
					      int way ,
					      double separator)
  :CbcBranchingObject(model,set->id(),way,0.5)
{
  set_ = set;
  separator_ = separator;
}

// Copy constructor 
CbcLinkBranchingObject::CbcLinkBranchingObject ( const CbcLinkBranchingObject & rhs) :CbcBranchingObject(rhs)
{
  set_=rhs.set_;
  separator_ = rhs.separator_;
}

// Assignment operator 
CbcLinkBranchingObject & 
CbcLinkBranchingObject::operator=( const CbcLinkBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    set_=rhs.set_;
    separator_ = rhs.separator_;
  }
  return *this;
}
CbcBranchingObject * 
CbcLinkBranchingObject::clone() const
{ 
  return (new CbcLinkBranchingObject(*this));
}


// Destructor 
CbcLinkBranchingObject::~CbcLinkBranchingObject ()
{
}
double
CbcLinkBranchingObject::branch(bool normalBranch)
{
  if (model_->messageHandler()->logLevel()>2&&normalBranch)
    print(normalBranch);
  numberBranchesLeft_--;
  int numberMembers = set_->numberMembers();
  int numberLinks = set_->numberLinks();
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
    int base=set_->first()+i*numberLinks;;
    for (;i<numberMembers;i++) {
      for (int k=0;k<numberLinks;k++) {
        int iColumn = base+k;
        solver->setColUpper(iColumn,0.0);
      }
      base += numberLinks;
    }
    way_=1;	  // Swap direction
  } else {
    int i;
    int base=set_->first();
    for ( i=0;i<numberMembers;i++) { 
      if (weights[i] >= separator_) {
	break;
      } else {
        for (int k=0;k<numberLinks;k++) {
          int iColumn = base+k;
          solver->setColUpper(iColumn,0.0);
        }
        base += numberLinks;
      }
    }
    assert (i<numberMembers);
    way_=-1;	  // Swap direction
  }
  return 0.0;
}
// Print what would happen  
void
CbcLinkBranchingObject::print(bool normalBranch)
{
  int numberMembers = set_->numberMembers();
  int numberLinks = set_->numberLinks();
  const double * weights = set_->weights();
  OsiSolverInterface * solver = model_->solver();
  const double * upper = solver->getColUpper();
  int first=numberMembers;
  int last=-1;
  int numberFixed=0;
  int numberOther=0;
  int i;
  int base=set_->first();
  for ( i=0;i<numberMembers;i++) {
    for (int k=0;k<numberLinks;k++) {
      int iColumn = base+k;
      double bound = upper[iColumn];
      if (bound) {
        first = CoinMin(first,i);
        last = CoinMax(last,i);
      }
    }
    base += numberLinks;
  }
  // *** for way - up means fix all those in down section
  base=set_->first();
  if (way_<0) {
    printf("SOS Down");
    for ( i=0;i<numberMembers;i++) {
      for (int k=0;k<numberLinks;k++) {
        int iColumn = base+k;
        double bound = upper[iColumn];
        if (weights[i] > separator_)
          break;
        else if (bound)
          numberOther++;
      }
      base += numberLinks;
    }
    assert (i<numberMembers);
    for (;i<numberMembers;i++) {
      for (int k=0;k<numberLinks;k++) {
        int iColumn = base+k;
        double bound = upper[iColumn];
        if (bound)
          numberFixed++;
      }
      base += numberLinks;
    }
  } else {
    printf("SOS Up");
    for ( i=0;i<numberMembers;i++) {
      for (int k=0;k<numberLinks;k++) {
        int iColumn = base+k;
        double bound = upper[iColumn];
        if (weights[i] >= separator_)
          break;
        else if (bound)
          numberFixed++;
      }
      base += numberLinks;
    }
    assert (i<numberMembers);
    for (;i<numberMembers;i++) {
      for (int k=0;k<numberLinks;k++) {
        int iColumn = base+k;
        double bound = upper[iColumn];
        if (bound)
          numberOther++;
      }
      base += numberLinks;
    }
  }
  assert ((numberFixed%numberLinks)==0);
  assert ((numberFixed%numberOther)==0);
  printf(" - at %g, free range %d (%g) => %d (%g), %d would be fixed, %d other way\n",
	 separator_,first,weights[first],last,weights[last],numberFixed/numberLinks,
         numberOther/numberLinks);
}
