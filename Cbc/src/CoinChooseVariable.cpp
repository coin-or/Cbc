// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CbcConfig.h"

#include <string>
#include <cassert>
#include <cfloat>
#include "OsiSolverInterface.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiSolverBranch.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinTime.hpp"
#include "CoinChooseVariable.hpp"
using namespace std;

CoinChooseVariable::CoinChooseVariable() :
  goodObjectiveValue_(COIN_DBL_MAX),
  goodSolution_(NULL),
  list_(NULL),
  useful_(NULL),
  solver_(NULL),
  state_(-1),
  numberUnsatisfied_(0),
  numberStrong_(0),
  numberOnList_(0),
  trustStrongForBound_(true),
  trustStrongForSolution_(true)
{
}

CoinChooseVariable::CoinChooseVariable(const OsiSolverInterface * solver) :
  goodObjectiveValue_(COIN_DBL_MAX),
  goodSolution_(NULL),
  solver_(solver),
  state_(-1),
  numberUnsatisfied_(0),
  numberStrong_(0),
  numberOnList_(0),
  trustStrongForBound_(true),
  trustStrongForSolution_(true)
{
  // create useful arrays
  int numberObjects = solver_->numberObjects();
  list_ = new int [numberObjects];
  useful_ = new double [numberObjects];
}

CoinChooseVariable::CoinChooseVariable(const CoinChooseVariable & rhs) 
{  
  goodObjectiveValue_ = rhs.goodObjectiveValue_;
  state_ = rhs.state_;
  numberUnsatisfied_ = rhs.numberUnsatisfied_;
  numberStrong_ = rhs.numberStrong_;
  numberOnList_ = rhs.numberOnList_;
  trustStrongForBound_ = rhs.trustStrongForBound_;
  trustStrongForSolution_ = rhs.trustStrongForSolution_;
  solver_ = rhs.solver_;
  if (solver_) {
    int numberObjects = solver_->numberObjects();
    int numberColumns = solver_->getNumCols();
    if (rhs.goodSolution_) {
      goodSolution_ = CoinCopyOfArray(rhs.goodSolution_,numberColumns);
    } else {
      goodSolution_ = NULL;
    }
    list_ = CoinCopyOfArray(rhs.list_,numberObjects);
    useful_ = CoinCopyOfArray(rhs.useful_,numberObjects);
  } else {
    goodSolution_ = NULL;
    list_ = NULL;
    useful_ = NULL;
  }
}

CoinChooseVariable &
CoinChooseVariable::operator=(const CoinChooseVariable & rhs)
{
  if (this != &rhs) {
    delete [] goodSolution_;
    delete [] list_;
    delete [] useful_;
    goodObjectiveValue_ = rhs.goodObjectiveValue_;
    state_ = rhs.state_;
    numberUnsatisfied_ = rhs.numberUnsatisfied_;
    numberStrong_ = rhs.numberStrong_;
    numberOnList_ = rhs.numberOnList_;
    trustStrongForBound_ = rhs.trustStrongForBound_;
    trustStrongForSolution_ = rhs.trustStrongForSolution_;
    solver_ = rhs.solver_;
    if (solver_) {
      int numberObjects = solver_->numberObjects();
      int numberColumns = solver_->getNumCols();
      if (rhs.goodSolution_) {
	goodSolution_ = CoinCopyOfArray(rhs.goodSolution_,numberColumns);
      } else {
	goodSolution_ = NULL;
      }
      list_ = CoinCopyOfArray(rhs.list_,numberObjects);
      useful_ = CoinCopyOfArray(rhs.useful_,numberObjects);
    } else {
      goodSolution_ = NULL;
      list_ = NULL;
      useful_ = NULL;
    }
  }
  return *this;
}


CoinChooseVariable::~CoinChooseVariable ()
{
  delete [] goodSolution_;
  delete [] list_;
  delete [] useful_;
}

// Clone
CoinChooseVariable *
CoinChooseVariable::clone() const
{
  return new CoinChooseVariable(*this);
}


// Initialize
void 
CoinChooseVariable::initialize ( OsiBranchingInformation *info)
{
  state_=0;
  assert (!goodSolution_);
  goodObjectiveValue_ = COIN_DBL_MAX;
}
/* Choose a variable
   Returns - 
   -1 Node is infeasible
   0  Normal termination - we have a candidate
   1  All looks satisfied - no candidate
   2  We can change the bound on a variable - but we also have a strong branching candidate
   3  We can change the bound on a variable - but we have a non-strong branching candidate
   4  We can change the bound on a variable - no other candidates
   We can pick up branch from whichObject() and whichWay()
   We can pick up a forced branch (can change bound) from whichForcedObject() and whichForcedWay()
   If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
*/
int 
CoinChooseVariable::chooseVariable( OsiBranchingInformation *info)
{
}
// Finish - deletes any solution etc
void 
CoinChooseVariable::finalize()
{
  state_=0;
  delete [] goodSolution_;
}
