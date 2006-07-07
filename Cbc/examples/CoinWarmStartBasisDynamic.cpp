// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>

#include "CoinWarmStartBasisDynamic.hpp"
#include <cmath>
#include <iostream>

CoinWarmStartBasisDynamic::CoinWarmStartBasisDynamic()
  : CoinWarmStartBasis(),
    numberCommonVariables_(0),
    numberCommonRows_(0),
    numberDynamicVariables_(0),
    dynamicVariables_(NULL)
{
}
CoinWarmStartBasisDynamic::CoinWarmStartBasisDynamic(int ns, int na, 
                                                     const char* sStat, const char* aStat,
                                                     int numberCommon, int numberDynamicVariables,
                                                     const int * dynamicVariables) 
  : CoinWarmStartBasis(ns, na, sStat, aStat),
    numberCommonVariables_(numberCommon),
    numberCommonRows_(0),
    numberDynamicVariables_(numberDynamicVariables)
{
  dynamicVariables_ = CoinCopyOfArray(dynamicVariables,numberDynamicVariables);
}

CoinWarmStartBasisDynamic::CoinWarmStartBasisDynamic(const CoinWarmStartBasisDynamic& rhs) 
  : CoinWarmStartBasis(rhs),
    numberCommonVariables_(rhs.numberCommonVariables_),
    numberCommonRows_(rhs.numberCommonRows_),
    numberDynamicVariables_(rhs.numberDynamicVariables_)
{
  dynamicVariables_ = CoinCopyOfArray(rhs.dynamicVariables_,numberDynamicVariables_);
}

CoinWarmStartBasisDynamic& 
CoinWarmStartBasisDynamic::operator=(const CoinWarmStartBasisDynamic& rhs)
{
  if (this != &rhs) {
    CoinWarmStartBasis::operator=(rhs);
    delete [] dynamicVariables_;
    numberCommonVariables_ = rhs.numberCommonVariables_;
    numberCommonRows_ = rhs.numberCommonRows_;
    numberDynamicVariables_ = rhs.numberDynamicVariables_;
    dynamicVariables_ = CoinCopyOfArray(rhs.dynamicVariables_,numberDynamicVariables_);
  }
  return *this;
}
CoinWarmStartBasisDynamic::~CoinWarmStartBasisDynamic()
{
  delete[] dynamicVariables_;
}
/*  Save list of dynamic variables */

void 
CoinWarmStartBasisDynamic::setDynamicVariables(int numberDynamicVariables, const int * dynamicVariables)
{
  delete[] dynamicVariables_;
  dynamicVariables_ = CoinCopyOfArray(dynamicVariables,numberDynamicVariables);
  numberDynamicVariables_=numberDynamicVariables;
}
void 
CoinWarmStartBasisDynamic::assignBasisStatus(int ns, int na, char*& sStat, 
                                             char*& aStat,
                                             int numberCommon, int numberDynamicVariables,
                                             int *& dynamicVariables) 
{
  CoinWarmStartBasis::assignBasisStatus(ns, na, sStat, aStat);
  delete[] dynamicVariables_;
  dynamicVariables_ = dynamicVariables;
  dynamicVariables=NULL;
  numberCommonVariables_=numberCommon;
  numberDynamicVariables_=numberDynamicVariables;
  assert (numberCommonVariables_+numberDynamicVariables_<=ns);
  numberCommonRows_=0;
}
// Deletes columns
void 
CoinWarmStartBasisDynamic::deleteColumns(int number, const int * which)
{
  int i ;
  int numStructural = getNumStructural();
  char * deleted = new char[numStructural];
  int numberDeleted=0;
  memset(deleted,0,numStructural*sizeof(char));
  for (i=0;i<number;i++) {
    int j = which[i];
    if (j>=0&&j<numStructural&&!deleted[j]) {
      numberDeleted++;
      deleted[j]=1;
    }
  }
  int put=0;
  for (i=0;i<numberDynamicVariables_;i++) {
    if (!deleted[i+numberCommonVariables_]) 
      dynamicVariables_[put++]=dynamicVariables_[i];
  }
  delete [] deleted;
  CoinWarmStartBasis::deleteColumns(number, which);
}
// Prints in readable format (for debug)
void 
CoinWarmStartBasisDynamic::print() const
{
  CoinWarmStartBasis::print();
  std::cout<<"There are "<<numberCommonVariables_<<" permanent and "
	   <<numberDynamicVariables_<<" dynamic variables"<<std::endl;
  int i;
  bool first=true;
  for (i=0;i<numberDynamicVariables_;i++) {
    if ((i%5)==0) {
      if (!first) 
        std::cout<<std::endl;
      else
        first=false;
      std::cout<<i;
    }
    std::cout<<" "<<dynamicVariables_[i];
  }
  std::cout<<std::endl;
}

CoinWarmStartBasisDiffDynamic &
CoinWarmStartBasisDiffDynamic::operator= (const CoinWarmStartBasisDiffDynamic &rhs) 
{
  if (this != &rhs) {
    CoinWarmStartBasisDiff::operator=(rhs);
    delete [] dynamic_;
    numberDynamic_ = rhs.numberDynamic_;
    dynamic_ = CoinCopyOfArray(rhs.dynamic_,numberDynamic_);
  }
  return *this;
}
CoinWarmStartBasisDiffDynamic::CoinWarmStartBasisDiffDynamic (const CoinWarmStartBasisDiffDynamic &rhs) 
  : CoinWarmStartBasisDiff(rhs)
{
  numberDynamic_ = rhs.numberDynamic_;
  dynamic_ = CoinCopyOfArray(rhs.dynamic_,numberDynamic_);
}
CoinWarmStartBasisDiffDynamic::CoinWarmStartBasisDiffDynamic (int sze, const unsigned int *const diffNdxs,
                                                              const unsigned int *const diffVals,
                                                              int numberDynamic, const int * dynamic) 
  : CoinWarmStartBasisDiff(sze, diffNdxs, diffVals)
{
  numberDynamic_ = numberDynamic;
  dynamic_ = CoinCopyOfArray(dynamic,numberDynamic_);
}

/*
  Generate a diff that'll convert oldCWS into the basis pointed to by this.

  This routine is a bit of a hack, for efficiency's sake. Rather than work
  with individual status vector entries, we're going to treat the vectors as
  int's --- in effect, we create one diff entry for each block of 16 status
  entries. Diffs for logicals are tagged with 0x80000000.
*/

CoinWarmStartDiff*
CoinWarmStartBasisDynamic::generateDiff (const CoinWarmStart *const oldCWS) const
{ 
/*
  Make sure the parameter is CoinWarmStartBasis or derived class.
*/
  const CoinWarmStartBasis *oldBasis =
      dynamic_cast<const CoinWarmStartBasis *>(oldCWS) ;
  if (!oldBasis)
  { throw CoinError("Old basis not derived from CoinWarmStartBasis.",
		    "generateDiff","CoinWarmStartBasis") ; }
  const CoinWarmStartBasis *newBasis = this ;
/*
  Make sure newBasis is equal or bigger than oldBasis. Calculate the worst case
  number of diffs and allocate vectors to hold them.
*/
  const int oldArtifCnt = oldBasis->getNumArtificial() ;
  const int oldStructCnt = oldBasis->getNumStructural() ;
  const int newArtifCnt = newBasis->getNumArtificial() ;
  const int newStructCnt = newBasis->getNumStructural() ;

  assert(newArtifCnt >= oldArtifCnt) ;
  assert(newStructCnt >= oldStructCnt) ;

  int sizeOldArtif = (oldArtifCnt+15)>>4 ;
  int sizeNewArtif = (newArtifCnt+15)>>4 ;
  int sizeOldStruct = (oldStructCnt+15)>>4 ;
  int sizeNewStruct = (newStructCnt+15)>>4 ;
  int maxBasisLength = sizeNewArtif+sizeNewStruct ;

  unsigned int *diffNdx = new unsigned int [maxBasisLength]; 
  unsigned int *diffVal = new unsigned int [maxBasisLength]; 
/*
  Ok, setup's over. Now scan the logicals (aka artificials, standing in for
  constraints). For the portion of the status arrays which overlap, create
  diffs. Then add any additional status from newBasis.

  I removed the following bit of code & comment:

    if (sizeNew == sizeOld) sizeOld--; // make sure all taken

  I assume this is meant to trap cases where oldBasis does not occupy all of
  the final int, but I can't see where it's necessary.
*/
  const unsigned int *oldStatus =
      reinterpret_cast<const unsigned int *>(oldBasis->getArtificialStatus()) ;
  const unsigned int *newStatus = 
      reinterpret_cast<const unsigned int *>(newBasis->getArtificialStatus()) ;
  int numberChanged = 0 ;
  int i ;
  for (i = 0 ; i < sizeOldArtif ; i++)
  { if (oldStatus[i] != newStatus[i])
    { diffNdx[numberChanged] = i|0x80000000 ;
      diffVal[numberChanged++] = newStatus[i] ; } }
  for ( ; i < sizeNewArtif ; i++)
  { diffNdx[numberChanged] = i|0x80000000 ;
    diffVal[numberChanged++] = newStatus[i] ; }
/*
  Repeat for structural variables.
*/
  oldStatus =
      reinterpret_cast<const unsigned int *>(oldBasis->getStructuralStatus()) ;
  newStatus =
      reinterpret_cast<const unsigned int *>(newBasis->getStructuralStatus()) ;
  for (i = 0 ; i < sizeOldStruct ; i++)
  { if (oldStatus[i] != newStatus[i])
    { diffNdx[numberChanged] = i ;
      diffVal[numberChanged++] = newStatus[i] ; } }
  for ( ; i < sizeNewStruct ; i++)
  { diffNdx[numberChanged] = i ;
    diffVal[numberChanged++] = newStatus[i] ; }
/*
  Create the object of our desire.
*/
  CoinWarmStartBasisDiffDynamic *diff =
    new CoinWarmStartBasisDiffDynamic(numberChanged,diffNdx,diffVal,numberDynamicVariables_,
                                      dynamicVariables_) ;
/*
  Clean up and return.
*/
  delete[] diffNdx ;
  delete[] diffVal ;

  return (dynamic_cast<CoinWarmStartDiff *>(diff)) ; }


