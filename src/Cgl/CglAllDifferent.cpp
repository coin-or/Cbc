// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>
//#define PRINT_DEBUG
//#define CGL_DEBUG 1
//#undef NDEBUG
#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinFinite.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CglAllDifferent.hpp"

#ifdef CGL_DEBUG
// A declaration is required somewhere, eh? I'm assuming static so the value
// carries over between calls to generateCuts.
namespace { int nPath = 0 ; }
#endif
//-------------------------------------------------------------------
// Generate cuts
//------------------------------------------------------------------- 
void CglAllDifferent::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			      const CglTreeInfo )
{
#ifndef NDEBUG
  int nCols=si.getNumCols();
#endif
  int i;
  const double * lower = si.getColLower();
  const double * upper = si.getColUpper();
#ifdef CGL_DEBUG
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
  if (debugger&&debugger->onOptimalPath(si)) {
    printf("On optimal path %d\n",nPath);
    nPath++;
    int nCols=si.getNumCols();
    const double * solution = si.getColSolution();
    const double * optimal = debugger->optimalSolution();
    const double * objective = si.getObjCoefficients();
    double objval1=0.0,objval2=0.0;
    for (i=0;i<nCols;i++) {
#if CGL_DEBUG>1
      printf("%d %g %g %g %g\n",i,lower[i],solution[i],upper[i],optimal[i]);
#endif
      objval1 += solution[i]*objective[i];
      objval2 += optimal[i]*objective[i];
      assert(optimal[i]>=lower[i]&&optimal[i]<=upper[i]);
    }
    printf("current obj %g, integer %g\n",objval1,objval2);
  }
#endif
  int * lo = new int[numberDifferent_];
  int * up = new int[numberDifferent_];
  for (i=0;i<numberDifferent_;i++) {
    int iColumn = originalWhich_[i];
    assert (iColumn<nCols);
    lo[i]  = static_cast<int> (lower[iColumn]);
    assert (floor(lower[iColumn]+0.5)==lower[iColumn]);
    up[i]  = static_cast<int> (upper[iColumn]);
    assert (floor(upper[iColumn]+0.5)==upper[iColumn]);
    assert (up[i]>=lo[i]);
  }
  // We are going to assume we can just have one big 2d array!
  // Could save by going to bits
  // also could skip sets where all are fixed
  // could do some of above by separate first pass
  // once a variable fixed - can take out of list
  // so need to redo complete stuff (including temp which_) every big pass
  int offset = COIN_INT_MAX;
  int maxValue = -COIN_INT_MAX;
  int numberLook=0;
  // copies
  //int * which = new int [numberTotal];
  //int * start = new int [numberSets_+1];
  for (i=0;i<numberSets_;i++) {
    for (int j=start_[i];j<start_[i+1];j++) {
      int k=which_[j];
      offset = std::min(offset,lo[k]);
      maxValue = std::max(maxValue,up[k]);
    }
    numberLook++;
    int gap = maxValue-offset+1;
    double size = static_cast<double> (gap) * numberDifferent_;
    if (size>1.0e7) {
      if (logLevel_)
        printf("Only looking at %d sets\n",numberLook);
      break;
    }
  }
  // Which sets a variable is in
  int * back = new int [start_[numberSets_]];
  int * backStart = new int[numberDifferent_+1];
  memset(backStart,0,(numberDifferent_+1)*sizeof(int));
  int numberTotal = start_[numberLook];
  for (i=0;i<numberTotal;i++) {
    int k=which_[i];
    // note +1 
    backStart[k+1]++;
  }
  int n=0;
  for (i=0;i<numberDifferent_;i++) {
    int nThis = backStart[i+1];
    backStart[i+1]=n;
    n+= nThis;
  }
  // at end all backStart correct!
  for (i=0;i<numberLook;i++) {
    for (int j=start_[i];j<start_[i+1];j++) {
      int k=which_[j];
      // note +1 
      int iPut = backStart[k+1];
      back[iPut]=i;
      backStart[k+1]=iPut+1;
    }
  }
  // value is possible for variable k if possible[k*gap+value] is nonzero
  int gap = maxValue-offset+1;
  char * possible = new char[gap*numberDifferent_];
  memset(possible,0,gap*numberDifferent_);
  // initialize
  int numberFixed=0;
  int * alreadyFixed = new int[numberDifferent_];
  for (i=0;i<numberDifferent_;i++) {
    alreadyFixed[i]=-1;
    int startV = i*gap + lo[i] - offset;
    int n = up[i]-lo[i]+1;
    memset(possible+startV,1,n);
  }
  for (i=0;i<numberDifferent_;i++) {
    int n = up[i]-lo[i]+1;
    if (n==1) {
      int fixedAt = lo[i]-offset;
      numberFixed++;
      alreadyFixed[i]=fixedAt;
      // take out of all others
      for (int j=backStart[i];j<backStart[i+1];j++) {
        int iSet = back[j];
        for (int jj=start_[iSet];jj<start_[iSet+1];jj++) {
          int k=which_[jj];
          if (k!=i) {
            // impossible
            possible[k*gap+fixedAt]=0;
          }
        }
      }
    }
  }
  bool finished=false;
  //int numberTightened=0;
  bool infeasible=false;
  // space to see which values possible
  int * check = new int[gap];
  unsigned int * bitmap = new unsigned int[numberDifferent_];
  int * stack = new int[numberDifferent_+1];
  int * first = new int[numberDifferent_+1];
  // just for valgrind etc
  memset(stack,0,(numberDifferent_+1)*sizeof(int));
  memset(first,0,(numberDifferent_+1)*sizeof(int));
  // do one set at a time
  while (!finished) {
    finished=true;
    int fixed=numberFixed;
    for (i=0;i<numberLook;i++) {
      memset(check,0,gap*sizeof(int));
      for (int j=start_[i];j<start_[i+1];j++) {
        int k=which_[j];
        if (alreadyFixed[k]>=0) {
          if (check[alreadyFixed[k]]==0) {
            check[alreadyFixed[k]]=1;
            continue;
          } else {
            // infeasible
            infeasible=true;
            i=numberLook;
            break;
          }
        }
        char * allowed = possible + k*gap;
        int n=0;
        for (int jj=0;jj<gap;jj++) {
          if (allowed[jj]) {
            n++;
            check[jj]++;
          }
        }
        if (n<2) {
          if (n==1) {
            // fix
            int fixedAt = -1;
            for (int jj=0;jj<gap;jj++) {
              if (allowed[jj]) {
                fixedAt=jj;
                break;
              }
            }
            numberFixed++;
            alreadyFixed[k]=fixedAt;
            check[fixedAt]=1;
            // take out of all others
            for (int j=backStart[k];j<backStart[k+1];j++) {
              int iSet = back[j];
              for (int jj=start_[iSet];jj<start_[iSet+1];jj++) {
                int kk=which_[jj];
                if (kk!=k) {
                  // impossible
                  possible[kk*gap+fixedAt]=0;
                }
              }
            }
          } else {
            // infeasible
            infeasible=true;
            j=numberTotal;
            i=numberLook;
            break;
          }
        }
      }
      // now check set
      // If number covered < number in set infeasible
      if (gap<30&&!infeasible) {
        int n=start_[i+1]-start_[i];
        memset(bitmap,0,n*sizeof(unsigned int));
        int j;
        int * which = which_+start_[i];
        unsigned int covered=0;
        bool good=true;
        for (j=0;j<n;j++) {
          int k=which[j];
          char * allowed = possible + k*gap;
          int jj;
          for (jj=0;jj<gap;jj++) 
            if (allowed[jj]) 
              break;
          assert (jj<gap);
          first[j]=jj;
          unsigned int iBit = 1<<jj;
          if ((covered&iBit)==0) {
            stack[j]=jj;
            covered |= iBit;
          } else {
            // can't
            jj++;
            for (;jj<gap;jj++) {
              iBit  = iBit << 1;
              if (allowed[jj]&&(covered&iBit)==0) 
                break;
            }
            if (jj<gap) {
              stack[j]=jj;
              covered |= iBit;
            } else {
              good = false;
              break;
            }
          }
        }
        int nStack=j;
        // just do first for rest
        for (;j<n;j++) {
          int k=which[j];
          char * allowed = possible + k*gap;
          int jj;
          for (jj=0;jj<gap;jj++) 
            if (allowed[jj]) 
              break;
          assert (jj<gap);
          first[j]=jj;
        }
        int kLook=0;
        while (nStack) {
          nStack--;
          if (good) {
#if 0
            printf("con %d = ",i);
            for (j=0;j<n;j++) 
              printf("%d ",stack[j]+1);
            printf("\n");
#endif
            // bug - kLook >= 0
            kLook=0;
            for (j=kLook;j<n;j++) {
              int iBit = 1 << stack[j];
              bitmap[j] |= iBit;
            }
          }
          kLook=nStack;
          int jj=stack[nStack];
          unsigned int iBit = 1<<jj;
          covered &= ~iBit;
          {
#ifndef NDEBUG
            unsigned int kBit=0;
            for (int k=0;k<nStack;k++) {
              int kk=stack[k];
              kBit |= 1<<kk;
            }
            assert (covered==kBit);
#endif
          }
          jj++;
          stack[nStack]=jj;
          while (nStack<n) {
            int k=which[nStack];
            char * allowed = possible + k*gap;
            for (;jj<gap;jj++) {
              iBit  = 1 << jj;
              if (allowed[jj]&&(covered&iBit)==0) 
                break;
            }
            if (jj<gap) {
              stack[nStack]=jj;
              covered |= iBit;
              nStack++;
              stack[nStack]=first[nStack];
              jj = first[nStack];
              good=true;
            } else {
              good = false;
              break;
            }
          }
        }
        int nnFix=0;
        // Now see if we can fix any
        for (j=0;j<n;j++) {
          int k=which[j];
          unsigned int mapped = bitmap[j];
          char * allowed = possible + k*gap;
          unsigned int iBit=1;
          for (int jj=0;jj<gap;jj++) {
            if ((mapped&iBit)==0) {
              if (allowed[jj]) {
                if (!nnFix)
                  printf("for con %d x ",i);
                nnFix++;
                printf("%d not %d ",j,jj+1);
                allowed[jj]=0;
                finished=false;
              }
            }
            iBit  = iBit << 1;
          }
        }
        if (nnFix)
          printf("\n");
      }
    }
    if (numberFixed>fixed)
      finished=false; // try again
  }
  // Could try two sets
  if (infeasible) {
    // create infeasible cut
    OsiRowCut rc;
    rc.setLb(COIN_DBL_MAX);
    rc.setUb(0.0);   
    cs.insertIfNotDuplicate(rc);
  } else {
    // check to see if can tighten bounds
    CoinPackedVector lbs;
    CoinPackedVector ubs;
    int nTightened=0;
    for (i=0;i<numberDifferent_;i++) {
      int iColumn = originalWhich_[i];
      char * allowed = possible+i*gap;
      int firstLo=-1;
      int lastUp=-1;
      for (int jj=0;jj<gap;jj++) {
        if (allowed[jj]) {
          if (firstLo<0)
            firstLo=jj;
          lastUp = jj;
        }
      }
      if (firstLo+offset>lo[i]) {
        lbs.insert(iColumn,static_cast<double> (firstLo+offset));
        nTightened++;
      }
      if (lastUp+offset<up[i]) {
        ubs.insert(iColumn,static_cast<double> (lastUp+offset));
        nTightened++;
      }
    }
    if (nTightened) {
      OsiColCut cc;
      cc.setUbs(ubs);
      cc.setLbs(lbs);
      cc.setEffectiveness(100.0);
      cs.insert(cc);
    }
  }
  //delete [] which;
  //delete [] start;
  delete [] first;
  delete [] stack;
  delete [] bitmap;
  delete [] check;
  delete [] alreadyFixed;
  delete [] back;
  delete [] backStart;
  delete [] possible;
  delete [] lo;
  delete [] up;
}

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglAllDifferent::CglAllDifferent ()
:
CglCutGenerator(),
numberSets_(0),
numberDifferent_(0),
maxLook_(2),
logLevel_(0),
start_(NULL),
which_(NULL),
originalWhich_(NULL)
{
}

//-------------------------------------------------------------------
// Useful Constructor 
//-------------------------------------------------------------------
CglAllDifferent::CglAllDifferent (int numberSets,
                                  const int * starts, const int * which)
:
CglCutGenerator(),
numberSets_(numberSets),
maxLook_(2),
logLevel_(0),
start_(NULL),
which_(NULL),
originalWhich_(NULL)
{
  if (numberSets_>0) {
    int n = starts[numberSets_];
    start_ = CoinCopyOfArray(starts,numberSets_+1);
    originalWhich_ = CoinCopyOfArray(which,n);
    which_ = new int[n];
    int i;
    int maxValue=-1;
    for (i=0;i<n;i++) {
      int iColumn = which[i];
      assert (iColumn>=0);
      maxValue = std::max(iColumn,maxValue);
    }
    maxValue++;
    int * translate = new int[maxValue];
    for (i=0;i<maxValue;i++)
      translate[i]=-1;
    for (i=0;i<n;i++) {
      int iColumn = which[i];
      translate[iColumn]=0;
    }
    numberDifferent_=0;
    for (i=0;i<maxValue;i++) {
      if (!translate[i]) 
        translate[i]=numberDifferent_++;
    }
    // Now translate
    for (i=0;i<n;i++) {
      int iColumn = which[i];
      iColumn = translate[iColumn];
      assert (iColumn>=0);
      which_[i]=iColumn;
    }
    delete [] translate;
  }
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglAllDifferent::CglAllDifferent (  const CglAllDifferent & rhs)
                                                              :
  CglCutGenerator(rhs),
  numberSets_(rhs.numberSets_),
  numberDifferent_(rhs.numberDifferent_),
  maxLook_(rhs.maxLook_),
  logLevel_(rhs.logLevel_)
{  
  if (numberSets_) {
    int n = rhs.start_[numberSets_];
    start_ = CoinCopyOfArray(rhs.start_,numberSets_+1);
    which_ = CoinCopyOfArray(rhs.which_,n);
    originalWhich_ = CoinCopyOfArray(rhs.originalWhich_,n);
  } else {
    start_=NULL;
    which_=NULL;
    originalWhich_=NULL;
  }
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglAllDifferent::clone() const
{
  return new CglAllDifferent(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglAllDifferent::~CglAllDifferent ()
{
  // free memory
  delete [] start_;
  delete [] which_;
  delete [] originalWhich_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglAllDifferent &
CglAllDifferent::operator=(
                                         const CglAllDifferent& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    // free memory
    delete [] start_;
    delete [] which_;
    delete [] originalWhich_;
    numberSets_ = rhs.numberSets_;
    numberDifferent_ = rhs.numberDifferent_;
    maxLook_ = rhs.maxLook_;
    logLevel_ = rhs.logLevel_;
    if (numberSets_) {
      int n = rhs.start_[numberSets_];
      start_ = CoinCopyOfArray(rhs.start_,numberSets_+1);
      which_ = CoinCopyOfArray(rhs.which_,n);
      originalWhich_ = CoinCopyOfArray(rhs.originalWhich_,n);
    } else {
      start_=NULL;
      which_=NULL;
      originalWhich_=NULL;
    }
  }
  return *this;
}

/// This can be used to refresh any inforamtion
void 
CglAllDifferent::refreshSolver(OsiSolverInterface * solver)
{
  // Get integer information
  solver->getColType(true);
}
// Create C++ lines to get to current state
std::string
CglAllDifferent::generateCpp( FILE * fp) 
{
  CglAllDifferent other;
  fprintf(fp,"0#include \"CglAllDifferent.hpp\"\n");
  fprintf(fp,"3  CglAllDifferent allDifferent;\n");
  if (logLevel_!=other.logLevel_)
    fprintf(fp,"3  allDifferent.setLogLevel(%d);\n",logLevel_);
  else
    fprintf(fp,"4  allDifferent.setLogLevel(%d);\n",logLevel_);
  if (maxLook_!=other.maxLook_)
    fprintf(fp,"3  allDifferent.setMaxLook(%d);\n",maxLook_);
  else
    fprintf(fp,"4  allDifferent.setMaxLook(%d);\n",maxLook_);
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  allDifferent.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  allDifferent.setAggressiveness(%d);\n",getAggressiveness());
  return "allDifferent";
}
