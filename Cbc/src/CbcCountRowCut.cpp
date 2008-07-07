// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>

#include "OsiRowCut.hpp"
#include "CbcCountRowCut.hpp"
#include "CbcNode.hpp"
//#define CHECK_CUT_COUNTS
// Default Constructor 
CbcCountRowCut::CbcCountRowCut ()
  :
  OsiRowCut(),
  owner_(NULL),
  ownerCut_(-1),
  numberPointingToThis_(0),
  whichCutGenerator_(-1)
{
#ifdef CBC_DETERMINISTIC_THREAD
  numberPointingToThis_=10;
#endif
#ifdef CHECK_CUT_COUNTS
  printf("CbcCountRowCut default constructor %x\n",this);
#endif
}
   
// Copy Constructor 
CbcCountRowCut::CbcCountRowCut (const OsiRowCut & rhs)
  : OsiRowCut(rhs),
    owner_(NULL),
    ownerCut_(-1),
    numberPointingToThis_(0),
    whichCutGenerator_(-1)
{
#ifdef CBC_DETERMINISTIC_THREAD
  numberPointingToThis_=10;
#endif
#ifdef CHECK_CUT_COUNTS
  printf("CbcCountRowCut constructor %x from RowCut\n",this);
#endif
}
// Copy Constructor 
CbcCountRowCut::CbcCountRowCut (const OsiRowCut & rhs,
				CbcNodeInfo * info, int whichOne,
				int whichGenerator)
  : OsiRowCut(rhs),
    owner_(info),
    ownerCut_(whichOne),
    numberPointingToThis_(0),
    whichCutGenerator_(whichGenerator)
{
#ifdef CBC_DETERMINISTIC_THREAD
  numberPointingToThis_=10;
#endif
#ifdef CHECK_CUT_COUNTS
  printf("CbcCountRowCut constructor %x from RowCut and info\n",this);
#endif
}
CbcCountRowCut::~CbcCountRowCut()
{
#ifdef CHECK_CUT_COUNTS
  printf("CbcCountRowCut destructor %x - references %d\n",this,
	 numberPointingToThis_);
#endif
  // Look at owner and delete
  owner_->deleteCut(ownerCut_);
  ownerCut_=-1234567;
}
// Increment number of references
void 
CbcCountRowCut::increment(int change)
{
  assert(ownerCut_!=-1234567);
#ifndef CBC_DETERMINISTIC_THREAD
  numberPointingToThis_+=change;
#endif
}

// Decrement number of references and return number left
int 
CbcCountRowCut::decrement(int change)
{
  assert(ownerCut_!=-1234567);
#ifndef CBC_DETERMINISTIC_THREAD
  //assert(numberPointingToThis_>=change);
  assert(numberPointingToThis_>=0);
  if(numberPointingToThis_<change) {
    assert(numberPointingToThis_>0);
    printf("negative cut count %d - %d\n",numberPointingToThis_, change);
    change = numberPointingToThis_;
  }
  numberPointingToThis_-=change;
#else
  assert (numberPointingToThis_==10);
#endif
  return numberPointingToThis_;
}

// Set information
void 
CbcCountRowCut::setInfo(CbcNodeInfo * info, int whichOne)
{
  owner_=info;
  ownerCut_=whichOne;
}
// Returns true if can drop cut if slack basic
bool 
CbcCountRowCut::canDropCut(const OsiSolverInterface * solver, int iRow) const
{
  // keep if COIN_DBL_MAX otherwise keep if slack zero
  if (effectiveness()<1.0e20) {
    return true;
  } else if (effectiveness()!=COIN_DBL_MAX) {
    if (iRow>=solver->getNumRows())
      return true;
    const double * rowActivity = solver->getRowActivity();
    const double * rowLower = solver->getRowLower();
    const double * rowUpper = solver->getRowUpper();
    double tolerance;
    solver->getDblParam(OsiPrimalTolerance,tolerance) ;
    double value = rowActivity[iRow];
    if (value<rowLower[iRow]+tolerance||
	value>rowUpper[iRow]-tolerance)
      return false;
    else
      return true;
  } else {
    return false;
  }
}

