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
   startBit_(NULL),
   numberBits_(NULL),
   rhs_(NULL)
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
   startBit_(NULL),
   numberBits_(NULL),
   rhs_(NULL)
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
  numberActive_ = 0;
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
  startBit_(NULL),
  numberBits_(NULL),
  rhs_(NULL)
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
  int type=gutsOfCheckPossible(1000000);
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
    // may be possible
    // for test - just if all 1
    if (type==0) {
      int i;

      int * indices = new int [numberActive_];
      double offset;
      solver->getDblParam(OsiObjOffset,offset);
      double fixedObj=-offset;
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
      int needed=0;
      int numberActive=0;
      for (i=0;i<numberRows;i++) {
        int newRow = lookup_[i];
        if (newRow>=0) {
          if (rowLower[i]==rowUpper[i]) {
            needed += 1<<numberActive;
          }
          numberActive++;
        }
      }
      double bestValue=COIN_DBL_MAX;
      int iBest=-1;
      for (i=0;i<size_;i++) {
        if ((i&needed)==needed) {
          // this one will do
          if (cost_[i]<bestValue) {
            bestValue=cost_[i];
            iBest=i;
          }
        }
      }
      returnCode=1;
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
        } else {
        }
      }
      delete [] indices;
    }
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
  i=0;
  bool touched = false;
  while (i<size_) {
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
      i++;
    } else {
      // we can skip some
      int k=i;
      int iBit=0;
      k &= ~1;
      while ((k&kMask)!=0) {
        iBit++;
        k &= ~(1<<iBit);
      }
      // onto next
      k += 1<<(iBit+1);
#ifdef CBC_DEBUG
      for (int j=i+1;j<k;j++) {
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

  
