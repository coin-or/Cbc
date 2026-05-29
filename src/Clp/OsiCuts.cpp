// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include <algorithm>
#include <cassert>

#include "OsiCuts.hpp"
#include "CoinHelperFunctions.hpp"

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
OsiCuts::OsiCuts()
  : rowCutPtrs_()
  , colCutPtrs_()
{
  // nothing to do here
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
OsiCuts::OsiCuts(const OsiCuts &source)
  : rowCutPtrs_()
  , colCutPtrs_()
{
  gutsOfCopy(source);
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
OsiCuts::~OsiCuts()
{
  gutsOfDestructor();
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
OsiCuts &
OsiCuts::operator=(const OsiCuts &rhs)
{
  if (this != &rhs) {
    gutsOfDestructor();
    gutsOfCopy(rhs);
  }
  return *this;
}

//-------------------------------------------------------------------
void OsiCuts::gutsOfCopy(const OsiCuts &source)
{
  assert(sizeRowCuts() == 0);
  assert(sizeColCuts() == 0);
  assert(sizeCuts() == 0);
  int i;
  int ne = source.sizeRowCuts();
  for (i = 0; i < ne; i++)
    insert(source.rowCut(i));
  ne = source.sizeColCuts();
  for (i = 0; i < ne; i++)
    insert(source.colCut(i));
}

//-------------------------------------------------------------------
void OsiCuts::gutsOfDestructor()
{
  int i;

  int ne = static_cast< int >(rowCutPtrs_.size());
  for (i = 0; i < ne; i++) {
    if (rowCutPtrs_[i]->globallyValidAsInteger() != 2)
      delete rowCutPtrs_[i];
  }
  rowCutPtrs_.clear();

  ne = static_cast< int >(colCutPtrs_.size());
  for (i = 0; i < ne; i++) {
    if (colCutPtrs_[i]->globallyValidAsInteger() != 2)
      delete colCutPtrs_[i];
  }
  colCutPtrs_.clear();

  assert(sizeRowCuts() == 0);
  assert(sizeColCuts() == 0);
  assert(sizeCuts() == 0);
}

//------------------------------------------------------------
//
// Embedded iterator class implementation
//
//------------------------------------------------------------
OsiCuts::iterator::iterator(OsiCuts &cuts)
  : cuts_(cuts)
  , rowCutIndex_(-1)
  , colCutIndex_(-1)
  , cutP_(NULL)
{
  this->operator++();
}

OsiCuts::iterator::iterator(const OsiCuts::iterator &src)
  : cuts_(src.cuts_)
  , rowCutIndex_(src.rowCutIndex_)
  , colCutIndex_(src.colCutIndex_)
  , cutP_(src.cutP_)
{
  // nothing to do here
}

OsiCuts::iterator &OsiCuts::iterator::operator=(const OsiCuts::iterator &rhs)
{
  if (this != &rhs) {
    cuts_ = rhs.cuts_;
    rowCutIndex_ = rhs.rowCutIndex_;
    colCutIndex_ = rhs.colCutIndex_;
    cutP_ = rhs.cutP_;
  }
  return *this;
}

OsiCuts::iterator::~iterator()
{
  //nothing to do
}

OsiCuts::iterator OsiCuts::iterator::begin()
{
  rowCutIndex_ = -1;
  colCutIndex_ = -1;
  this->operator++();
  return *this;
}

OsiCuts::iterator OsiCuts::iterator::end()
{
  rowCutIndex_ = cuts_.sizeRowCuts();
  colCutIndex_ = cuts_.sizeColCuts() - 1;
  cutP_ = NULL;
  return *this;
}

OsiCuts::iterator OsiCuts::iterator::operator++()
{
  cutP_ = NULL;

  // Are there any more row cuts to consider?
  if ((rowCutIndex_ + 1) >= cuts_.sizeRowCuts()) {
    // Only column cuts left.
    colCutIndex_++;
    // Only update cutP if there is a column cut.
    // This is necessary for the iterator to work for
    // OsiCuts that don't have any cuts.
    if (cuts_.sizeColCuts() > 0 && colCutIndex_ < cuts_.sizeColCuts())
      cutP_ = cuts_.colCutPtr(colCutIndex_);
  }

  // Are there any more col cuts to consider?
  else if ((colCutIndex_ + 1) >= cuts_.sizeColCuts()) {
    // Only row cuts left
    rowCutIndex_++;
    if (rowCutIndex_ < cuts_.sizeRowCuts())
      cutP_ = cuts_.rowCutPtr(rowCutIndex_);
  }

  // There are still Row & column cuts left to consider
  else {
    double nextColCutE = cuts_.colCut(colCutIndex_ + 1).effectiveness();
    double nextRowCutE = cuts_.rowCut(rowCutIndex_ + 1).effectiveness();
    if (nextColCutE > nextRowCutE) {
      colCutIndex_++;
      cutP_ = cuts_.colCutPtr(colCutIndex_);
    } else {
      rowCutIndex_++;
      cutP_ = cuts_.rowCutPtr(rowCutIndex_);
    }
  }
  return *this;
}

//------------------------------------------------------------
//
// Embedded const_iterator class implementation
//
//------------------------------------------------------------
OsiCuts::const_iterator::const_iterator(const OsiCuts &cuts)
  : cutsPtr_(&cuts)
  , rowCutIndex_(-1)
  , colCutIndex_(-1)
  , cutP_(NULL)
{
  this->operator++();
}

OsiCuts::const_iterator::const_iterator(const OsiCuts::const_iterator &src)
  : cutsPtr_(src.cutsPtr_)
  , rowCutIndex_(src.rowCutIndex_)
  , colCutIndex_(src.colCutIndex_)
  , cutP_(src.cutP_)
{
  // nothing to do here
}

OsiCuts::const_iterator &
OsiCuts::const_iterator::operator=(const OsiCuts::const_iterator &rhs)
{
  if (this != &rhs) {
    cutsPtr_ = rhs.cutsPtr_;
    rowCutIndex_ = rhs.rowCutIndex_;
    colCutIndex_ = rhs.colCutIndex_;
    cutP_ = rhs.cutP_;
  }
  return *this;
}

OsiCuts::const_iterator::~const_iterator()
{
  //nothing to do
}

OsiCuts::const_iterator OsiCuts::const_iterator::begin()
{
  rowCutIndex_ = -1;
  colCutIndex_ = -1;
  this->operator++();
  return *this;
}

OsiCuts::const_iterator OsiCuts::const_iterator::end()
{
  rowCutIndex_ = cutsPtr_->sizeRowCuts();
  colCutIndex_ = cutsPtr_->sizeColCuts() - 1;
  cutP_ = NULL;
  return *this;
}

OsiCuts::const_iterator OsiCuts::const_iterator::operator++()
{
  cutP_ = NULL;

  // Are there any more row cuts to consider?
  if ((rowCutIndex_ + 1) >= cutsPtr_->sizeRowCuts()) {
    // Only column cuts left.
    colCutIndex_++;
    // Only update cutP if there is a column cut.
    // This is necessary for the iterator to work for
    // OsiCuts that don't have any cuts.
    if (cutsPtr_->sizeRowCuts() > 0 && colCutIndex_ < cutsPtr_->sizeColCuts())
      cutP_ = cutsPtr_->colCutPtr(colCutIndex_);
  }

  // Are there any more col cuts to consider?
  else if ((colCutIndex_ + 1) >= cutsPtr_->sizeColCuts()) {
    // Only row cuts left
    rowCutIndex_++;
    if (rowCutIndex_ < cutsPtr_->sizeRowCuts())
      cutP_ = cutsPtr_->rowCutPtr(rowCutIndex_);
  }

  // There are still Row & column cuts left to consider
  else {
    double nextColCutE = cutsPtr_->colCut(colCutIndex_ + 1).effectiveness();
    double nextRowCutE = cutsPtr_->rowCut(rowCutIndex_ + 1).effectiveness();
    if (nextColCutE > nextRowCutE) {
      colCutIndex_++;
      cutP_ = cutsPtr_->colCutPtr(colCutIndex_);
    } else {
      rowCutIndex_++;
      cutP_ = cutsPtr_->rowCutPtr(rowCutIndex_);
    }
  }
  return *this;
}
static bool scaleCutIntegral(double* cutElem, int* cutIndex, int cutNz,
			     double& cutRhs, double maxdelta);
/* Insert a row cut unless it is a duplicate (CoinAbsFltEq)
       returns true if inserted */
bool OsiCuts::insertIfNotDuplicate(OsiRowCut &rc, CoinAbsFltEq treatAsSame)
{
  double newLb = rc.lb();
  double newUb = rc.ub();
  CoinPackedVector vector = rc.row();
  int numberElements = vector.getNumElements();
  int *newIndices = vector.getIndices();
  double *newElements = vector.getElements();
  CoinSort_2(newIndices, newIndices + numberElements, newElements);
  bool notDuplicate = true;
  int numberRowCuts = sizeRowCuts();
  for (int i = 0; i < numberRowCuts; i++) {
    const OsiRowCut *cutPtr = rowCutPtr(i);
    if (cutPtr->row().getNumElements() != numberElements)
      continue;
    if (!treatAsSame(cutPtr->lb(), newLb))
      continue;
    if (!treatAsSame(cutPtr->ub(), newUb))
      continue;
    const CoinPackedVector *thisVector = &(cutPtr->row());
    const int *indices = thisVector->getIndices();
    const double *elements = thisVector->getElements();
    int j;
    for (j = 0; j < numberElements; j++) {
      if (indices[j] != newIndices[j])
        break;
      if (!treatAsSame(elements[j], newElements[j]))
        break;
    }
    if (j == numberElements) {
      notDuplicate = false;
      break;
    }
  }
  if (notDuplicate) {
    OsiRowCut *newCutPtr = new OsiRowCut();
    newCutPtr->setRow(vector);
    newCutPtr->setLb(newLb);
    newCutPtr->setUb(newUb);
    newCutPtr->setGloballyValid(rc.globallyValid());
    newCutPtr->setEffectiveness(rc.effectiveness());
    rowCutPtrs_.push_back(newCutPtr);
  }
  return notDuplicate;
}

/* Insert a row cut unless it is a duplicate (CoinRelFltEq)*/
void OsiCuts::insertIfNotDuplicate(OsiRowCut &rc, CoinRelFltEq treatAsSame)
{
  double newLb = rc.lb();
  double newUb = rc.ub();
  CoinPackedVector vector = rc.row();
  int numberElements = vector.getNumElements();
  int *newIndices = vector.getIndices();
  double *newElements = vector.getElements();
  CoinSort_2(newIndices, newIndices + numberElements, newElements);
  bool notDuplicate = true;
  int numberRowCuts = sizeRowCuts();
  for (int i = 0; i < numberRowCuts; i++) {
    const OsiRowCut *cutPtr = rowCutPtr(i);
    if (cutPtr->row().getNumElements() != numberElements)
      continue;
    if (!treatAsSame(cutPtr->lb(), newLb))
      continue;
    if (!treatAsSame(cutPtr->ub(), newUb))
      continue;
    const CoinPackedVector *thisVector = &(cutPtr->row());
    const int *indices = thisVector->getIndices();
    const double *elements = thisVector->getElements();
    int j;
    for (j = 0; j < numberElements; j++) {
      if (indices[j] != newIndices[j])
        break;
      if (!treatAsSame(elements[j], newElements[j]))
        break;
    }
    if (j == numberElements) {
      notDuplicate = false;
      break;
    }
  }
  if (notDuplicate) {
    OsiRowCut *newCutPtr = new OsiRowCut();
    newCutPtr->setLb(newLb);
    newCutPtr->setUb(newUb);
    newCutPtr->setRow(vector);
    newCutPtr->setGloballyValid(rc.globallyValid());
    newCutPtr->setEffectiveness(rc.effectiveness());
    rowCutPtrs_.push_back(newCutPtr);
  }
}
#define USE__RATIONAL 60
#include "CoinRational.hpp"
#define SMALL_VALUE1 1.0e-14
static long computeGcd(long a, long b) {
  // This is the standard Euclidean algorithm for gcd
  long remainder = 1;
  // Make sure a<=b (will always remain so)
  if (a > b) {
    // Swap a and b
    long temp = a;
    a = b;
    b = temp;
  }
  // If zero then gcd is nonzero
  if (!a) {
    if (b) {
      return b;
    }
    else {
      return 0;
    }
  }
  while (remainder) {
    remainder = b % a;
    b = a;
    a = remainder;
  }
  return b;
} /* computeGcd */
static bool scaleCutIntegral(double* cutElem, int* cutIndex, int cutNz,
			     double& cutRhs, double maxdelta) {
  long gcd, lcm;
  double maxscale = 1000;
  long maxdnom = USE__RATIONAL;
  //long numerator = 0, denominator = 0;
  // Initialize gcd and lcm
  CoinRational r = CoinRational(cutRhs, maxdelta, maxdnom);
  if (r.getNumerator() != 0){
     gcd = labs(r.getNumerator());
     lcm = r.getDenominator();
  }
  else{
    return false;
  }
  for (int i = 0; i < cutNz; ++i) {
    CoinRational r = CoinRational(cutElem[i], maxdelta, maxdnom);
    if (r.getNumerator() != 0){
       gcd = computeGcd(gcd, r.getNumerator());
       lcm *= r.getDenominator()/(computeGcd(lcm,r.getDenominator()));
    }
    else{
      return false;
    }
  }
  double scale = ((double)lcm)/((double)gcd);
  if (fabs(scale) > maxscale) {
      return false;
  }
  scale = fabs(scale);
  // Pre-check that every scaled value rounds to an integer before modifying
  // anything, so the caller's arrays are never left in a partially-scaled state.
  for (int i = 0; i < cutNz; ++i) {
    double value = cutElem[i]*scale;
    if (fabs(floor(value+0.5)-value) >= 1.0e-9) return false;
  }
  {
    double value = cutRhs*scale;
    if (fabs(floor(value+0.5)-value) >= 1.0e-9) return false;
  }
  // Looks like we have a good scaling factor; scale and return.
  for (int i = 0; i < cutNz; ++i) {
    cutElem[i] = floor(cutElem[i]*scale+0.5);
  }
  cutRhs = floor(cutRhs*scale+0.5);
  return true;
} /* scaleCutIntegral */
// Returns value - floor but allowing for small errors
inline double above_integer(double value) {
  double value2=floor(value);
  double value3=floor(value+0.5);
  if (fabs(value3-value)<1.0e-9*(fabs(value3)+1.0))
    return 0.0;
  return value-value2;
}
/* Insert a row cut unless it is a duplicate - cut may get sorted.
   Duplicate is defined as CoinAbsFltEq says same
   returns true if inserted.
   Also tries to "integerize" cut and checks for accuracy */
bool OsiCuts::insertIfNotDuplicateAndClean(OsiRowCut &rc,
					   int typeCut ,
					   CoinAbsFltEq treatAsSame)
{
  double newLb = rc.lb();
  double newUb = rc.ub();
  CoinPackedVector vector = rc.row();
  int numberElements = vector.getNumElements();
  int *newIndices = vector.getIndices();
  double *newElements = vector.getElements();
  CoinSort_2(newIndices, newIndices + numberElements, newElements);
  bool notDuplicate = true;
  int numberRowCuts = sizeRowCuts();
  for (int i = 0; i < numberRowCuts; i++) {
    const OsiRowCut *cutPtr = rowCutPtr(i);
    if (cutPtr->row().getNumElements() != numberElements)
      continue;
    if (!treatAsSame(cutPtr->lb(), newLb))
      continue;
    if (!treatAsSame(cutPtr->ub(), newUb))
      continue;
    const CoinPackedVector *thisVector = &(cutPtr->row());
    const int *indices = thisVector->getIndices();
    const double *elements = thisVector->getElements();
    int j;
    for (j = 0; j < numberElements; j++) {
      if (indices[j] != newIndices[j])
        break;
      if (!treatAsSame(elements[j], newElements[j]))
        break;
    }
    if (j == numberElements) {
      notDuplicate = false;
      break;
    }
  }
  if (notDuplicate) {
    OsiRowCut *newCutPtr = new OsiRowCut();
    // scale
    if (newLb<-1.0e30 || newUb > 1.0e30) {
      double maxdelta = 1.0e-12; //param.getEPS();
      double rhs = newLb<-1.0e30 ? newUb : newLb;
      bool goodScale = scaleCutIntegral(newElements,newIndices,numberElements,
					rhs,maxdelta);
      if (goodScale) {
	if (newLb<-1.0e30)
	  newUb = rhs;
	else
	  newLb = rhs;
      } else {
	// scale
	double rhs = fabs((newLb<-1.0e30) ? newUb : newLb);
	double largest = 1.0e-100;
	double smallest = 1.0e100;
	for (int i=0;i<numberElements;i++) {
	  double value = fabs(newElements[i]);
	  largest = std::max(largest,value);
	  smallest = std::min(smallest,value);
	}
	if (largest>smallest*1.0e8||rhs>smallest*1.0e8) {
	  //printf ("badly scaled cut - rhs %g els %g -> %g - type %d\n",rhs,smallest,largest,
	  //	  typeCut);
	  //if (rhs>1.0e5 && smallest < 1.0e-4)
	  if (typeCut==61)  {// CglTwoMir
	    delete newCutPtr;
	    return false;
	  }
	}
	// TwoMIR cut derivation can accumulate floating-point errors of up to ~1e-6.
	// When integral scaling fails, relax the RHS slightly so that numerically
	// tight cuts do not incorrectly exclude the optimal integer solution.
	if (typeCut >= 61 && typeCut <= 63) {
	  if (newLb < -1.0e30)
	    newUb += 2.0e-6;
	  else
	    newLb -= 2.0e-6;
	}
      }
      newCutPtr->setRow(numberElements,newIndices,newElements);
    } else {
      newCutPtr->setRow(vector);
    }
    newCutPtr->setLb(newLb);
    newCutPtr->setUb(newUb);
    newCutPtr->setGloballyValid(rc.globallyValid());
    newCutPtr->setEffectiveness(rc.effectiveness());
    rowCutPtrs_.push_back(newCutPtr);
  }
  return notDuplicate;
}

