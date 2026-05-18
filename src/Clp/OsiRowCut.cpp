// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cfloat>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>

#include "CoinPragma.hpp"
#include "CoinFinite.hpp"
#include "OsiRowCut.hpp"

#ifndef OSI_INLINE_ROWCUT_METHODS

//  //-------------------------------------------------------------------
//  // Set/Get lower & upper bounds
//  //-------------------------------------------------------------------
double OsiRowCut::lb() const { return lb_; }

void OsiRowCut::setLb(double lb) { lb_ = lb; }

double OsiRowCut::ub() const { return ub_; }

void OsiRowCut::setUb(double ub) { ub_ = ub; }

//-------------------------------------------------------------------
// Set row elements
//-------------------------------------------------------------------
void OsiRowCut::setRow(int size,
  const int *colIndices, const double *elements,
  bool testForDuplicateIndex)
{
  row_.setVector(size, colIndices, elements, testForDuplicateIndex);
}

void OsiRowCut::setRow(const CoinPackedVector &v)
{
  row_ = v;
}

void OsiRowCut::setRow(CoinPackedVector &&v)
{
  row_ = std::move(v);
}

//-------------------------------------------------------------------
// Get the row
//-------------------------------------------------------------------
const CoinPackedVector &OsiRowCut::row() const
{
  return row_;
}

//-------------------------------------------------------------------
// Get the row for changing
//-------------------------------------------------------------------
CoinPackedVector &OsiRowCut::mutableRow()
{
  return row_;
}

//----------------------------------------------------------------
// == operator
//-------------------------------------------------------------------
bool OsiRowCut::operator==(const OsiRowCut &rhs) const
{
  if (this->OsiCut::operator!=(rhs))
    return false;
  if (row() != rhs.row())
    return false;
  if (lb() != rhs.lb())
    return false;
  if (ub() != rhs.ub())
    return false;

  return true;
}

bool OsiRowCut::operator!=(const OsiRowCut &rhs) const
{
  return !((*this) == rhs);
}

//----------------------------------------------------------------
// consistent & infeasible
//-------------------------------------------------------------------
bool OsiRowCut::consistent() const
{
  const CoinPackedVector &r = row();
  r.duplicateIndex("consistent", "OsiRowCut");
  if (r.getMinIndex() < 0)
    return false;
  return true;
}

bool OsiRowCut::consistent(const OsiSolverInterface &im) const
{
  const CoinPackedVector &r = row();
  if (r.getMaxIndex() >= im.getNumCols())
    return false;

  return true;
}

bool OsiRowCut::infeasible(const OsiSolverInterface &) const
{
  if (lb() > ub())
    return true;

  return false;
}

#endif
/* Returns infeasibility of the cut with respect to solution
    passed in i.e. is positive if cuts off that solution.
    solution is getNumCols() long..
*/
double
OsiRowCut::violated(const double *solution) const
{
  int i;
  double sum = 0.0;
  const int *column = row_.getIndices();
  int number = row_.getNumElements();
  const double *element = row_.getElements();
  for (i = 0; i < number; i++) {
    int colIndx = column[i];
    sum += solution[colIndx] * element[i];
  }
  if (sum > ub_)
    return sum - ub_;
  else if (sum < lb_)
    return lb_ - sum;
  else
    return 0.0;
}

//-------------------------------------------------------------------
// Row sense, rhs, range
//-------------------------------------------------------------------
char OsiRowCut::sense() const
{
  if (lb_ == ub_)
    return 'E';
  else if (lb_ == -COIN_DBL_MAX && ub_ == COIN_DBL_MAX)
    return 'N';
  else if (lb_ == -COIN_DBL_MAX)
    return 'L';
  else if (ub_ == COIN_DBL_MAX)
    return 'G';
  else
    return 'R';
}

double OsiRowCut::rhs() const
{
  if (lb_ == ub_)
    return ub_;
  else if (lb_ == -COIN_DBL_MAX && ub_ == COIN_DBL_MAX)
    return 0.0;
  else if (lb_ == -COIN_DBL_MAX)
    return ub_;
  else if (ub_ == COIN_DBL_MAX)
    return lb_;
  else
    return ub_;
}

double OsiRowCut::range() const
{
  if (lb_ == ub_)
    return 0.0;
  else if (lb_ == -COIN_DBL_MAX && ub_ == COIN_DBL_MAX)
    return 0.0;
  else if (lb_ == -COIN_DBL_MAX)
    return 0.0;
  else if (ub_ == COIN_DBL_MAX)
    return 0.0;
  else
    return ub_ - lb_;
}

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
OsiRowCut::OsiRowCut()
  : OsiCut()
  , row_()
  , lb_(-COIN_DBL_MAX)
  , ub_(COIN_DBL_MAX)
{
  //#ifdef NDEBUG
  //row_.setTestForDuplicateIndex(false);
  //#endif
}

//-------------------------------------------------------------------
// Ownership constructor
//-------------------------------------------------------------------

OsiRowCut::OsiRowCut(double cutlb, double cutub,
  int capacity, int size,
  int *&colIndices, double *&elements)
  : OsiCut()
  , row_(capacity, size, colIndices, elements)
  , lb_(cutlb)
  , ub_(cutub)
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
OsiRowCut::OsiRowCut(const OsiRowCut &source)
  : OsiCut(source)
  , row_(source.row_)
  , lb_(source.lb_)
  , ub_(source.ub_)
{
  // Nothing to do here
}

//----------------------------------------------------------------
// Move constructor
//----------------------------------------------------------------
OsiRowCut::OsiRowCut(OsiRowCut &&rhs) noexcept
  : OsiCut(rhs)
  , row_(std::move(rhs.row_))
  , lb_(rhs.lb_)
  , ub_(rhs.ub_)
{
}

//----------------------------------------------------------------
// Move assignment
//----------------------------------------------------------------
OsiRowCut &
OsiRowCut::operator=(OsiRowCut &&rhs) noexcept
{
  if (this != &rhs) {
    OsiCut::operator=(rhs);
    row_ = std::move(rhs.row_);
    lb_ = rhs.lb_;
    ub_ = rhs.ub_;
  }
  return *this;
}

//----------------------------------------------------------------
// Clone
//----------------------------------------------------------------
OsiRowCut *OsiRowCut::clone() const
{
  return (new OsiRowCut(*this));
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
OsiRowCut::~OsiRowCut()
{
  // Nothing to do here
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
OsiRowCut &
OsiRowCut::operator=(const OsiRowCut &rhs)
{
  if (this != &rhs) {
    OsiCut::operator=(rhs);
    row_ = rhs.row_;
    lb_ = rhs.lb_;
    ub_ = rhs.ub_;
  }
  return *this;
}

//----------------------------------------------------------------
// Print
//-------------------------------------------------------------------

void OsiRowCut::print() const
{
  int i;
  char temp[80];
  std::cout << "Row cut has " << row_.getNumElements()
            << " elements";
#define MORE_PREC
#ifndef MORE_PREC
  if (lb_ < -1.0e20 && ub_ < 1.0e20)
    std::cout << " with upper rhs of " << ub_;
  else if (lb_ > -1.0e20 && ub_ > 1.0e20)
    std::cout << " with lower rhs of " << lb_;
  else
    std::cout << " !!! with lower, upper rhs of " << lb_ << " and " << ub_;
#else
  sprintf(temp," lower %.10g, upper %.10g",lb_,ub_);
  std::cout << temp;
#endif
  std::cout << std::endl;
  for (i = 0; i < row_.getNumElements(); i++) {
    int colIndx = row_.getIndices()[i];
    double element = row_.getElements()[i];
    if (i > 0 && element > 0)
      std::cout << " +";
    std::cout << element << " * x" << colIndx << " ";
  }
  std::cout << std::endl;
}

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
OsiRowCut2::OsiRowCut2(int row)
  : OsiRowCut()
  , whichRow_(row)
{
  // nothing to do here
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
OsiRowCut2::OsiRowCut2(const OsiRowCut2 &source)
  : OsiRowCut(source)
  , whichRow_(source.whichRow_)
{
  // Nothing to do here
}

//----------------------------------------------------------------
// Clone
//----------------------------------------------------------------
OsiRowCut *OsiRowCut2::clone() const
{
  return (new OsiRowCut2(*this));
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
OsiRowCut2::~OsiRowCut2()
{
  // Nothing to do here
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
OsiRowCut2 &
OsiRowCut2::operator=(const OsiRowCut2 &rhs)
{
  if (this != &rhs) {
    OsiRowCut::operator=(rhs);
    whichRow_ = rhs.whichRow_;
  }
  return *this;
}

static const double multiplier[] = {1.23456789e2,-9.87654321};
int hashCut (const OsiRowCut &x, int size)
{
    int xN =x.row().getNumElements();
    double xLb = x.lb();
    double xUb = x.ub();
    const int * xIndices = x.row().getIndices();
    const double * xElements = x.row().getElements();
    unsigned int hashValue;
    double value=1.0;
    if (xLb>-1.0e10)
        value += xLb*multiplier[0];
    if (xUb<1.0e10)
        value += xUb*multiplier[1];
    for( int j=0;j<xN;j++) {
        int xColumn = xIndices[j];
        double xValue = xElements[j];
        int k=(j&1);
        value += (j+1)*multiplier[k]*(xColumn+1)*xValue;
    }
    // should be compile time but too lazy for now
    union { double d; unsigned int i[2]; } xx;
    if (sizeof(value)>sizeof(hashValue)) {
        assert (sizeof(value)==2*sizeof(hashValue));
        xx.d = value;
        hashValue = (xx.i[0] + xx.i[1]);
    } else {
        assert (sizeof(value)==sizeof(hashValue));
        xx.d = value;
        hashValue = xx.i[0];
    }
    return hashValue%(size);
}

bool same(const OsiRowCut &x, const OsiRowCut &y)
{
    int xN = x.row().getNumElements();
    int yN = y.row().getNumElements();
    bool identical = false;
    if (xN == yN) {
        double xLb = x.lb();
        double xUb = x.ub();
        double yLb = y.lb();
        double yUb = y.ub();
        if (std::fabs(xLb - yLb) < 1.0e-8 && fabs(xUb - yUb) < 1.0e-8) {
            const int *xIndices = x.row().getIndices();
            const double *xElements = x.row().getElements();
            const int *yIndices = y.row().getIndices();
            const double *yElements = y.row().getElements();
            int j;
            for (j = 0; j < xN; j++) {
                if (xIndices[j] != yIndices[j])
                    break;
                if (fabs(xElements[j] - yElements[j]) > 1.0e-12)
                    break;
            }
            identical = (j == xN);
        }
    }
    return identical;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
