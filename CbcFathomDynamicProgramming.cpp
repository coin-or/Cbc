// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcFathomDynamicProgramming.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinSort.hpp"
// Default Constructor
CbcFathomDynamicProgramming::CbcFathomDynamicProgramming() 
  :CbcFathom(),
   size_(0),
   type_(-1),
   cost_(NULL),
   back_(NULL),
   id_(NULL),
   lookup_(NULL),
   numberActive_(0),
   maximumSizeAllowed_(1000000),
   startBit_(NULL),
   numberBits_(NULL),
   rhs_(NULL),
   numberNonOne_(0)
{

}

// Constructor from model
CbcFathomDynamicProgramming::CbcFathomDynamicProgramming(CbcModel & model)
  :CbcFathom(model),
   cost_(NULL),
   back_(NULL),
   id_(NULL),
   lookup_(NULL),
   numberActive_(0),
   maximumSizeAllowed_(1000000),
   startBit_(NULL),
   numberBits_(NULL),
   rhs_(NULL),
   numberNonOne_(0)
{
  type_=gutsOfCheckPossible();
}

// Destructor 
CbcFathomDynamicProgramming::~CbcFathomDynamicProgramming ()
{
  gutsOfDelete();
}
// Does deleteions
void 
CbcFathomDynamicProgramming::gutsOfDelete()
{
  delete [] cost_;
  delete [] back_;
  delete [] id_;
  delete [] lookup_;
  delete [] startBit_;
  delete [] numberBits_;
  delete [] rhs_;
  cost_ = NULL;
  back_ = NULL;
  id_ = NULL;
  lookup_ = NULL;
  startBit_ = NULL;
  numberBits_ = NULL;
  rhs_ = NULL;
}
// Clone
CbcFathom *
CbcFathomDynamicProgramming::clone() const
{
  return new CbcFathomDynamicProgramming(*this);
}

// Copy constructor 
CbcFathomDynamicProgramming::CbcFathomDynamicProgramming(const CbcFathomDynamicProgramming & rhs)
:
  CbcFathom(rhs),
  size_(rhs.size_),
  type_(rhs.type_),
  cost_(NULL),
  back_(NULL),
  id_(NULL),
  lookup_(NULL),
  numberActive_(rhs.numberActive_),
  maximumSizeAllowed_(rhs.maximumSizeAllowed_),
  startBit_(NULL),
  numberBits_(NULL),
  rhs_(NULL),
  numberNonOne_(rhs.numberNonOne_)
{
  if (size_) {
    cost_=CoinCopyOfArray(rhs.cost_,size_);
    back_=CoinCopyOfArray(rhs.back_,size_);
    id_=CoinCopyOfArray(rhs.id_,size_);
    int numberRows=model_->getNumRows();
    lookup_=CoinCopyOfArray(rhs.lookup_,numberRows);
    startBit_=CoinCopyOfArray(rhs.startBit_,numberActive_);
    numberBits_=CoinCopyOfArray(rhs.numberBits_,numberActive_);
    rhs_=CoinCopyOfArray(rhs.rhs_,numberActive_);
  }
}
// Returns type
int 
CbcFathomDynamicProgramming::gutsOfCheckPossible(int allowableSize)
{
  assert(model_->solver());
  OsiSolverInterface * solver = model_->solver();
  const CoinPackedMatrix * matrix = solver->getMatrixByCol();

  int numberIntegers = model_->numberIntegers();
  int numberColumns = solver->getNumCols();
  size_=0;
  if (numberIntegers!=numberColumns)
    return -1; // can't do dynamic programming

  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  const double * rowUpper = solver->getRowUpper();

  int numberRows = model_->getNumRows();
  int i;

  // First check columns to see if possible
  double * rhs = new double [numberRows];
  CoinCopyN(rowUpper,numberRows,rhs);

  // Column copy
  const double * element = matrix->getElements();
  const int * row = matrix->getIndices();
  const CoinBigIndex * columnStart = matrix->getVectorStarts();
  const int * columnLength = matrix->getVectorLengths();
  bool bad=false;
  /* It is just possible that we could say okay as
     variables may get fixed but seems unlikely */
  for (i=0;i<numberColumns;i++) {
    int j;
    double lowerValue = lower[i];
    assert (lowerValue==floor(lowerValue));
    for (j=columnStart[i];
         j<columnStart[i]+columnLength[i];j++) {
      int iRow=row[j];
      double value = element[j];
      if (upper[i]>lowerValue&&(value<=0.0||value!=floor(value)))
        bad=true;
      if (lowerValue)
	rhs[iRow] -= lowerValue*value;
    }
  }
  // check possible (at present do not allow covering)
  int numberActive=0;
  for (i=0;i<numberRows;i++) {
    if (rhs[i]>1.0e5||fabs(rhs[i]-floor(rhs[i]+0.5))>1.0e-7)
      bad=true;
    else if (rhs[i]>0.0)
      numberActive++;
  }
  if (bad) {
    delete [] rhs;
    return -1;
  }
  // check size of array needed
  double size=1.0;
  double check = INT_MAX;
  for (i=0;i<numberRows;i++) {
    int n= (int) floor(rhs[i]+0.5);
    if (n) {
      n++; // allow for 0,1... n
      if (numberActive!=1) {
        // power of 2
        int iBit=0;
        int k=n;
        k &= ~1;
        while (k) {
          iBit++;
          k &= ~(1<<iBit);
        }
        // See if exact power
        if (n!=(1<<iBit)) {
          // round up to next power of 2
          n= 1<<(iBit+1);
        }
        size *= n;
        if (size>=check)
          break;
      } else {
        size = n; // just one constraint
      }
    }
  }
  // set size needed
  if (size>=check)
    size_=INT_MAX;
  else
    size_=(int) size;
        
  int n01=0;
  int nbadcoeff=0;
  // See if we can tighten bounds
  for (i=0;i<numberColumns;i++) {
    int j;
    double lowerValue = lower[i];
    double gap = upper[i]-lowerValue;
    for (j=columnStart[i];
         j<columnStart[i]+columnLength[i];j++) {
      int iRow=row[j];
      double value = element[j];
      if (value!=1.0)
        nbadcoeff++;
      if (gap*value>rhs[iRow]+1.0e-8)
        gap = rhs[iRow]/value;
    }
    gap=lowerValue+floor(gap+1.0e-7);
    if (gap<upper[i])
      solver->setColUpper(i,gap);
    if (gap<=1.0)
      n01++;
  }
  if (allowableSize&&size_<=allowableSize) {
    numberActive_=numberActive;
    cost_ = new double [size_];
    CoinFillN(cost_,size_,COIN_DBL_MAX);
    // but do nothing is okay
    cost_[0]=0.0;
    back_ = new int[size_];
    CoinFillN(back_,size_,-1);
    id_ = new int[size_];
    CoinFillN(id_,size_,-1);
    startBit_=new int[numberActive_];
    numberBits_=new int[numberActive_];
    lookup_ = new int [numberRows];
    rhs_ = new int [numberActive_];
    numberActive=0;
    int kBit=0;
    for (i=0;i<numberRows;i++) {
      int n= (int) floor(rhs[i]+0.5);
      if (n) {
        lookup_[i]=numberActive;
        rhs_[numberActive]=n;
        startBit_[numberActive]=kBit;
        n++; // allow for 0,1... n
        int iBit=0;
        // power of 2
        int k=n;
        k &= ~1;
        while (k) {
          iBit++;
          k &= ~(1<<iBit);
        }
        // See if exact power
        if (n!=(1<<iBit)) {
          // round up to next power of 2
          iBit++;
        }
        if (numberActive!=1) {
          n= 1<<iBit;
          size *= n;
          if (size>=check)
            break;
        } else {
          size = n; // just one constraint
        }
        numberBits_[numberActive++]=iBit;
        kBit += iBit;
      } else {
        lookup_[i]=-1;
      }
    }
  }
  delete [] rhs;
  if (allowableSize&&size_>allowableSize) {
    printf("Too large - need %d entries x 16 bytes\n",size_);
    return -1; // too big
  }
  if (n01==numberColumns&&!nbadcoeff)
    return 0; // easiest
  else
    return 1;
}

// Resets stuff if model changes
void 
CbcFathomDynamicProgramming::resetModel(CbcModel * model)
{
  model_=model;
  type_=gutsOfCheckPossible();
}
int
CbcFathomDynamicProgramming::fathom(double * & betterSolution)
{
  int returnCode=0;
  int type=gutsOfCheckPossible(maximumSizeAllowed_);
  if (type>=0) {
    bool gotSolution=false;
    OsiSolverInterface * solver = model_->solver();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    const double * objective = solver->getObjCoefficients();
    double direction = solver->getObjSense();
    const CoinPackedMatrix * matrix = solver->getMatrixByCol();
    // Column copy
    const double * element = matrix->getElements();
    const int * row = matrix->getIndices();
    const CoinBigIndex * columnStart = matrix->getVectorStarts();
    const int * columnLength = matrix->getVectorLengths();
    const double * rowLower = solver->getRowLower();
    const double * rowUpper = solver->getRowUpper();
    int numberRows = model_->getNumRows();
    
    int numberColumns = solver->getNumCols();
    int * indices = new int [numberActive_];
    double offset;
    solver->getDblParam(OsiObjOffset,offset);
    double fixedObj=-offset;
    int i;
    // may be possible
    if (type==0) {
      // rhs 1 and coefficients 1
      for (i=0;i<numberColumns;i++) {
        int j;
        double lowerValue = lower[i];
        assert (lowerValue==floor(lowerValue));
        double cost = direction * objective[i];
        fixedObj += lowerValue*cost;
        if (lowerValue==upper[i])
          continue;
        int n=0;
        for (j=columnStart[i];
             j<columnStart[i]+columnLength[i];j++) {
          int iRow=row[j];
          double value = element[j];
          int newRow = lookup_[iRow];
          if (newRow<0||value>rhs_[newRow]) {
            n=0;
            break; //can't use
          } else {
            indices[n++]=newRow;
          }
        }
        if (n) {
          addOneColumn0(i,n,indices,cost);
        }
      }
      returnCode=1;
    } else {
      // more complex
      int * coefficients = new int[numberActive_];
      // If not too many general rhs then we can be more efficient
      numberNonOne_=0;
      for (i=0;i<numberActive_;i++) {
        if (rhs_[i]!=1)
          numberNonOne_++;
      }
      if (numberNonOne_*2<numberActive_) {
        // put rhs >1 every second
        int * permute = new int[numberActive_];
        int * temp = new int[numberActive_];
        // try different ways
        int k=0;
        for (i=0;i<numberRows;i++) {
          int newRow = lookup_[i];
          if (newRow>=0&&rhs_[newRow]>1) {
            permute[newRow]=k;
            k +=2;
          }
        }
        // adjust so k points to last
        k -= 2;
        // and now rest
        int k1=1;
        for (i=0;i<numberRows;i++) {
          int newRow = lookup_[i];
          if (newRow>=0&&rhs_[newRow]==1) {
            permute[newRow]=k1;
            k1++;
            if (k1<=k)
              k1++;
          }
        }
        for (i=0;i<numberActive_;i++) {
          int put = permute[i];
          temp[put]=rhs_[i];
        }
        memcpy(rhs_,temp,numberActive_*sizeof(int));
        for (i=0;i<numberActive_;i++) {
          int put = permute[i];
          temp[put]=numberBits_[i];
        }
        memcpy(numberBits_,temp,numberActive_*sizeof(int));
        k=0;
        for (i=0;i<numberActive_;i++) {
          startBit_[i]=k;
          k+= numberBits_[i];
        }
        for (i=0;i<numberRows;i++) {
          int newRow = lookup_[i];
          if (newRow>=0)
            lookup_[i]=permute[newRow];
        }
        delete [] permute;
        delete [] temp;
        // mark new method
        type=2;
      }
      for (i=0;i<numberColumns;i++) {
        int j;
        double lowerValue = lower[i];
        assert (lowerValue==floor(lowerValue));
        double cost = direction * objective[i];
        fixedObj += lowerValue*cost;
        if (lowerValue==upper[i])
          continue;
        int n=0;
        double gap=upper[i]-lowerValue;
        for (j=columnStart[i];
             j<columnStart[i]+columnLength[i];j++) {
          int iRow=row[j];
          double value = element[j];
          int iValue = (int) value;
          int newRow = lookup_[iRow];
          if (newRow<0||iValue>rhs_[newRow]) {
            n=0;
            break; //can't use
          } else {
            coefficients[n]=iValue;
            indices[n++]=newRow;
            if (gap*iValue>rhs_[newRow]) {
              gap = rhs_[newRow]/iValue;
            }
          }
        }
        int nTotal = (int) gap;
        if (n) {
          if (type==1) {
            for (int k=1;k<=nTotal;k++) {
              addOneColumn1(i,n,indices,coefficients,cost);
            }
          } else {
            // may be able to take sort out with new method
            CoinSort_2(indices,indices+n,coefficients);
            for (int k=1;k<=nTotal;k++) {
              addOneColumn1A(i,n,indices,coefficients,cost);
            }
          }
        }
      }
      delete [] coefficients;
      returnCode=1;
    }
    int needed=0;
    double bestValue=COIN_DBL_MAX;
    int iBest=-1;
    if (type==0) {
      int numberActive=0;
      for (i=0;i<numberRows;i++) {
        int newRow = lookup_[i];
        if (newRow>=0) {
          if (rowLower[i]==rowUpper[i]) {
            needed += 1<<numberActive;
            numberActive++;
          }
        }
      }
      for (i=0;i<size_;i++) {
        if ((i&needed)==needed) {
          // this one will do
          if (cost_[i]<bestValue) {
            bestValue=cost_[i];
            iBest=i;
          }
        }
      }
    } else {
      int * lower = new int[numberActive_];
      for (i=0;i<numberRows;i++) {
        int newRow = lookup_[i];
        if (newRow>=0) {
          int gap=(int) (rowUpper[i]-CoinMax(0.0,rowLower[i]));
          lower[newRow]=rhs_[newRow]-gap;
          int numberBits = numberBits_[newRow];
          int startBit = startBit_[newRow];
          if (numberBits==1&&!gap) {
            needed |= 1<<startBit;
          }
        }
      }
      for (i=0;i<size_;i++) {
        if ((i&needed)==needed) {
        // this one may do
          bool good = true;
          for (int kk=0;kk<numberActive_;kk++) {
            int numberBits = numberBits_[kk];
            int startBit = startBit_[kk];
            int size = 1<<numberBits;
            int start = 1<<startBit;
            int mask = start*(size-1);
            int level=(i&mask)>>startBit;
            if (level<lower[kk]) {
              good=false;
              break;
            }
          }
          if (good&&cost_[i]<bestValue) {
            bestValue=cost_[i];
            iBest=i;
          }
        }
      }
      delete [] lower;
    }
    if (bestValue<COIN_DBL_MAX) {
      bestValue += fixedObj;
      printf("Can get solution of %g\n",bestValue);
      if (bestValue<model_->getMinimizationObjValue()) {
        // set up solution
        betterSolution = new double[numberColumns];
        memcpy(betterSolution,lower,numberColumns*sizeof(double));
        while (iBest>0) {
          int iColumn = id_[iBest];
          assert (iColumn>=0);
          betterSolution[iColumn]++;
          assert (betterSolution[iColumn]<=upper[iColumn]);
          iBest = back_[iBest];
        }
      }
    }
    delete [] indices;
    gutsOfDelete();
    if (gotSolution) {
      int i;
      // paranoid check
      double * rowActivity = new double [numberRows];
      memset(rowActivity,0,numberRows*sizeof(double));
      for (i=0;i<numberColumns;i++) {
        int j;
        double value = betterSolution[i];
        if (value) {
          for (j=columnStart[i];
               j<columnStart[i]+columnLength[i];j++) {
            int iRow=row[j];
            rowActivity[iRow] += value*element[j];
          }
        }
      }
      // check was feasible
      bool feasible=true;
      for (i=0;i<numberRows;i++) {
        if(rowActivity[i]<rowLower[i]) {
          if (rowActivity[i]<rowLower[i]-1.0e-8)
            feasible = false;
        } else if(rowActivity[i]>rowUpper[i]) {
          if (rowActivity[i]>rowUpper[i]+1.0e-8)
            feasible = false;
        }
      }
      if (feasible) {
        if (model_->messageHandler()->logLevel()>1)
          printf("** good solution by dynamic programming\n");
      }
      delete [] rowActivity;
    }
  }
  return returnCode;
}
/* Adds one column if type 0,
   returns true if was used in making any changes
*/
bool 
CbcFathomDynamicProgramming::addOneColumn0(int id,int numberElements, const int * rows,
                     double cost)
{
  // build up mask
  int mask=0;
  int i;
  for (i=0;i<numberElements;i++) {
    int iRow=rows[i];
    mask |= 1<<iRow;
  }
  i=size_-1-mask;
  bool touched = false;
  while (i>=0) {
    int kMask = i&mask;
    if (kMask==0) {
      double thisCost = cost_[i];
      if (thisCost!=COIN_DBL_MAX) {
        // possible
        double newCost=thisCost+cost;
        int next = i + mask;
        if (cost_[next]>newCost) {
          cost_[next]=newCost;
          back_[next]=i;
          id_[next]=id;
          touched=true;
        }
      }
      i--;
    } else {
      // we can skip some
      int k=(i&~mask)-1;
#ifdef CBC_DEBUG
      for (int j=i-1;j>k;j--) {
        int jMask = j&mask;
        assert (jMask!=0);
      }
#endif
      i=k;
    }
  }
  return touched;
}
/* Adds one attempt of one column of type 1,
   returns true if was used in making any changes.
   At present the user has to call it once for each possible value
*/
bool 
CbcFathomDynamicProgramming::addOneColumn1(int id,int numberElements, const int * rows,
                                           const int * coefficients, double cost)
{
  /* build up masks.
     a) mask for 1 rhs
     b) mask for addition
     c) mask so adding will overflow
     d) individual masks
  */
  int mask1=0;
  int maskAdd=0;
  int mask2=0;
  int i;
  int n2=0;
  int mask[40];
  int nextbit[40];
  int adjust[40];
  assert (numberElements<=40);
  for (i=0;i<numberElements;i++) {
    int iRow=rows[i];
    int numberBits = numberBits_[iRow];
    int startBit = startBit_[iRow];
    if (numberBits==1) {
      mask1 |= 1<<startBit;
      maskAdd |= 1<<startBit;
      mask2 |= 1<<startBit;
    } else {
      int value=coefficients[i];
      int size = 1<<numberBits;
      int start = 1<<startBit;
      assert (value<size);
      maskAdd |= start*value;
      int gap = size-rhs_[iRow]-1;
      assert (gap>=0);
      int hi2=rhs_[iRow]-value;
      if (hi2<size-1)
        hi2++;
      adjust[n2] = start*hi2;
      mask2 += start*gap;
      nextbit[n2]=startBit+numberBits;
      mask[n2++] = start*(size-1);
    }
  }
  i=size_-1-maskAdd;
  bool touched = false;
  while (i>=0) {
    int kMask = i&mask1;
    if (kMask==0) {
      bool good=true;
      for (int kk=n2-1;kk>=0;kk--) {
        int iMask = mask[kk];
        int jMask = iMask&mask2;
        int kkMask = iMask&i;
        kkMask += jMask;
        if (kkMask>iMask) {
          // we can skip some
          int k=(i&~iMask);
          k |= adjust[kk];  
#ifdef CBC_DEBUG
          for (int j=i-1;j>k;j--) {
            int jMask = j&mask1;
            if (jMask==0) {
              bool good=true;
              for (int kk=n2-1;kk>=0;kk--) {
                int iMask = mask[kk];
                int jMask = iMask&mask2;
                int kkMask = iMask&i;
                kkMask += jMask;
                if (kkMask>iMask) {
                  good=false;
                  break;
                }
              }
              assert (!good);
            }
          }
#endif
          i=k;
          good=false;
          break;
        }
      }
      if (good) {
        double thisCost = cost_[i];
        if (thisCost!=COIN_DBL_MAX) {
          // possible
          double newCost=thisCost+cost;
          int next = i + maskAdd;
          if (cost_[next]>newCost) {
            cost_[next]=newCost;
            back_[next]=i;
            id_[next]=id;
            touched=true;
          }
        }
      }
      i--;
    } else {
      // we can skip some
      // we can skip some
      int k=(i&~mask1)-1;
#ifdef CBC_DEBUG
      for (int j=i-1;j>k;j--) {
        int jMask = j&mask;
        assert (jMask!=0);
      }
#endif
      i=k;
    }
  }
  return touched;
}
/* Adds one attempt of one column of type 1,
   returns true if was used in making any changes.
   At present the user has to call it once for each possible value
   This version is when there are enough 1 rhs to do faster
*/
bool 
CbcFathomDynamicProgramming::addOneColumn1A(int id,int numberElements, const int * rows,
                                           const int * coefficients, double cost)
{
  /* build up masks.
     a) mask for 1 rhs
     b) mask for addition
     c) mask so adding will overflow
     d) mask for non 1 rhs
  */
  int maskA=0;
  int maskAdd=0;
  int maskC=0;
  int maskD=0;
  int i;
  for (i=0;i<numberElements;i++) {
    int iRow=rows[i];
    int numberBits = numberBits_[iRow];
    int startBit = startBit_[iRow];
    if (numberBits==1) {
      maskA |= 1<<startBit;
      maskAdd |= 1<<startBit;
    } else {
      int value=coefficients[i];
      int size = 1<<numberBits;
      int start = 1<<startBit;
      assert (value<size);
      maskAdd |= start*value;
      int gap = size-rhs_[iRow]-1;
      assert (gap>=0);
      maskC |= start*gap;
      maskD |= start*(size-1);
    }
  }
  int maskDiff = maskD-maskC;
  i=size_-1-maskAdd;
  bool touched = false;
  while (i>=0) {
    int kMask = i&maskA;
    if (kMask==0) {
      int added = i & maskD; // just bits belonging to non 1 rhs
      added += maskC; // will overflow mask if bad
      added &= (~maskD);
      if (added == 0) {
        double thisCost = cost_[i];
        if (thisCost!=COIN_DBL_MAX) {
          // possible
          double newCost=thisCost+cost;
          int next = i + maskAdd;
          if (cost_[next]>newCost) {
            cost_[next]=newCost;
            back_[next]=i;
            id_[next]=id;
            touched=true;
          }
        }
        i--;
      } else {
        // we can skip some 
        int k = i & ~ maskD; // clear all
        // Put back enough - but only below where we are
        int kk=(numberNonOne_<<1)-2;
        assert (rhs_[kk]>1);
        int iMask=0;
        for(;kk>=0;kk-=2) {
          iMask = 1<<startBit_[kk+1];
          if ((added&iMask)!=0) {
            iMask--;
            break;
          }
        }
        assert (kk>=0);
        iMask &= maskDiff;
        k |= iMask;
        assert (k<i);
        i=k;
      }
    } else {
      // we can skip some
      // we can skip some
      int k=(i&~maskA)-1;
#ifdef CBC_DEBUG
      for (int j=i-1;j>k;j--) {
        int jMask = j&mask;
        assert (jMask!=0);
      }
#endif
      i=k;
    }
  }
  return touched;
}
// update model
void CbcFathomDynamicProgramming::setModel(CbcModel * model)
{
  model_ = model;
  type_=gutsOfCheckPossible();
}

  
