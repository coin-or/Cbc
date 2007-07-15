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
// Default Constructor 
CbcCountRowCut::CbcCountRowCut ()
  :
  OsiRowCut(),
  owner_(NULL),
  ownerCut_(-1),
  numberPointingToThis_(0),
  whichCutGenerator_(-1)
{
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
  numberPointingToThis_+=change;
}

// Decrement number of references and return number left
int 
CbcCountRowCut::decrement(int change)
{
  assert(ownerCut_!=-1234567);
  //assert(numberPointingToThis_>=change);
  assert(numberPointingToThis_>=0);
  if(numberPointingToThis_<change) {
    assert(numberPointingToThis_>0);
    printf("negative cut count %d - %d\n",numberPointingToThis_, change);
    change = numberPointingToThis_;
  }
  numberPointingToThis_-=change;
  return numberPointingToThis_;
}

// Set information
void 
CbcCountRowCut::setInfo(CbcNodeInfo * info, int whichOne)
{
  owner_=info;
  ownerCut_=whichOne;
}

