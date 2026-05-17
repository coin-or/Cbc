// Copyright (C) 2010, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cfloat> 
#include <cassert>

#include "CoinPragma.hpp"
#include "CglZeroHalf.hpp" 
#include "CoinPackedVector.hpp"
#include "CoinSort.hpp"
#include "CoinPackedMatrix.hpp"

//-------------------------------------------------------------
void
CglZeroHalf::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
				const CglTreeInfo info)
{
  si.getColType(true);
  if (mnz_) {
    int cnum=0,cnzcnt=0;
    int *cbeg=NULL, *ccnt=NULL,*cind=NULL,*cval=NULL,*crhs=NULL;
    char *csense=NULL;
    const double * solution = si.getColSolution();
    if ((flags_&1)==0) {
      // redo bounds
      const double * columnLower = si.getColLower();
      const double * columnUpper = si.getColUpper();
      int numberColumns = si.getNumCols();
      for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	if (vlb_[iColumn]!=COIN_INT_MAX) {
	  int ilo,iup;
	  double lo = columnLower[iColumn];
	  if (lo<-COIN_INT_MAX)
	    lo=-COIN_INT_MAX;
	  ilo= static_cast<int> (ceil(lo));
	  double up = columnUpper[iColumn];
	  if (up>COIN_INT_MAX)
	    up=COIN_INT_MAX;
	  iup= static_cast<int> (floor(up));
	  vlb_[iColumn]=ilo;
	  vub_[iColumn]=iup;
	}
      }
    }
    if (true) {
    cutInfo_.setMaxSeconds(getMaxSeconds());
    cutInfo_.sep_012_cut(mr_,mc_,mnz_,
				 mtbeg_,mtcnt_, mtind_, mtval_,
				 vlb_, vub_,
				 mrhs_, msense_,
				 solution,
				 info.inTree ? false : true,
				 &cnum,&cnzcnt,
				 &cbeg,&ccnt,&cind,&cval,&crhs,&csense);
    } else {
      int k = 4*mr_+2*mnz_;
      int * temp = new int[k];
      int * mtbeg = temp;
      int * mtcnt = mtbeg + mr_;
      int * mtind = mtcnt+mr_;
      int * mtval = mtind+mnz_;
      int * mrhs = mtval+mnz_;
      char * msense = reinterpret_cast<char*> (mrhs+mr_);
      int i;
      k=0;
      int kel=0;
      for (i=0;i<mr_;i++) {
	int kel2=kel;
	int rhs = mrhs_[i];
	for (int j=mtbeg_[i];j<mtbeg_[i]+mtcnt_[i];j++) {
	  int iColumn=mtind_[j];
	  int value=mtval_[j];
	  if (vlb_[iColumn]<vub_[iColumn]) {
	    mtind[kel]=mtind_[j];
	    mtval[kel++]=mtval_[j];
	  } else {
	    rhs -= vlb_[iColumn]*value;
	  }
	}
	if (kel>kel2) {
	  mtcnt[k]=kel-kel2;
	  mtbeg[k]=kel2;
	  mrhs[k]=rhs;
	  msense[k++]=msense_[i];
	}
      }
      if (kel) {
	cutInfo_.sep_012_cut(k,mc_,kel,
				 mtbeg,mtcnt, mtind, mtval,
				 vlb_, vub_,
				 mrhs, msense,
				 solution,
				 info.inTree ? false : true,
				 &cnum,&cnzcnt,
				 &cbeg,&ccnt,&cind,&cval,&crhs,&csense);
      }
      delete [] temp;
    }
    if (cnum) {
      // add cuts
      double * element = new double[mc_];
      for (int i=0;i<cnum;i++) {
	int n = ccnt[i];
	int start = cbeg[i];
	for (int j=0;j<n;j++) 
	  element[j]=cval[start+j];
	OsiRowCut rc;
	if (csense[i]=='L') {
	  rc.setLb(-COIN_DBL_MAX);
	  rc.setUb(crhs[i]);
	} else if (csense[i]=='G') {
	  rc.setLb(crhs[i]);
	  rc.setUb(COIN_DBL_MAX);
	} else {
	  abort();
	}
	rc.setRow(n,cind+start,element,false);
	if ((flags_&1)!=0)
	  rc.setGloballyValid();
	//double violation = rc.violated(solution);
	//if (violation>1.0e-6)
	  cs.insertIfNotDuplicate(rc);
	  //else
	  //printf("violation of %g\n",violation);
      }
      delete [] element;
      free(cbeg); 
      free(ccnt);
      free(cind);
      free(cval);
      free(crhs);
      free(csense);
    }
  }
}

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglZeroHalf::CglZeroHalf ()
:
CglCutGenerator(),
  mr_(0),
  mc_(0),
  mnz_(0),
  mtbeg_(NULL),
  mtcnt_(NULL),
  mtind_(NULL),
  mtval_(NULL),
  vlb_(NULL),
  vub_(NULL),
  mrhs_(NULL),
  msense_(NULL),
  flags_(0)
{
  cutInfo_=Cgl012Cut();
}
//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglZeroHalf::CglZeroHalf (
                  const CglZeroHalf & source)
:
  CglCutGenerator(source),
  mtbeg_(NULL),
  mtcnt_(NULL),
  mtind_(NULL),
  mtval_(NULL),
  vlb_(NULL),
  vub_(NULL),
  mrhs_(NULL),
  msense_(NULL),
  flags_(source.flags_)
{  
  mr_ = source.mr_;
  mc_ = source.mc_;
  mnz_ = source.mnz_;
  if (mr_) {
    mtbeg_ = CoinCopyOfArray(source.mtbeg_,mr_);
    mtcnt_ = CoinCopyOfArray(source.mtcnt_,mr_);
    mtind_ = CoinCopyOfArray(source.mtind_,mnz_);
    mtval_ = CoinCopyOfArray(source.mtval_,mnz_);
    vlb_ = CoinCopyOfArray(source.vlb_,mc_);
    vub_ = CoinCopyOfArray(source.vub_,mc_);
    mrhs_ = CoinCopyOfArray(source.mrhs_,mr_);
    msense_ = CoinCopyOfArray(source.msense_,mr_);
  }
  //cutInfo_ = Cgl012Cut(source.cutInfo_);
  cutInfo_ = Cgl012Cut();
  cutInfo_.setSepGraphSparseThreshold(source.cutInfo_.getSepGraphSparseThreshold());
  cutInfo_.setRowMaxPairCount(source.cutInfo_.getRowMaxPairCount());
  cutInfo_.setRowMaxFractionalCount(source.cutInfo_.getRowMaxFractionalCount());
}


//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglZeroHalf::clone() const
{
  return new CglZeroHalf(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglZeroHalf::~CglZeroHalf ()
{
  delete []  mtbeg_;
  delete []  mtcnt_;
  delete []  mtind_;
  delete []  mtval_;
  delete []  vlb_;
  delete []  vub_;
  delete []  mrhs_;
  delete []  msense_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglZeroHalf &
CglZeroHalf::operator=(
                   const CglZeroHalf& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    delete []  mtbeg_;
    delete []  mtcnt_;
    delete []  mtind_;
    delete []  mtval_;
    delete []  vlb_;
    delete []  vub_;
    delete []  mrhs_;
    delete []  msense_;
    mr_ = rhs.mr_;
    mc_ = rhs.mc_;
    mnz_ = rhs.mnz_;
    flags_ = rhs.flags_;
    if (mr_) {
      mtbeg_ = CoinCopyOfArray(rhs.mtbeg_,mr_);
      mtcnt_ = CoinCopyOfArray(rhs.mtcnt_,mr_);
      mtind_ = CoinCopyOfArray(rhs.mtind_,mnz_);
      mtval_ = CoinCopyOfArray(rhs.mtval_,mnz_);
      vlb_ = CoinCopyOfArray(rhs.vlb_,mc_);
      vub_ = CoinCopyOfArray(rhs.vub_,mc_);
      mrhs_ = CoinCopyOfArray(rhs.mrhs_,mr_);
      msense_ = CoinCopyOfArray(rhs.msense_,mr_);
    } else {
      mtbeg_ = NULL;
      mtcnt_ = NULL;
      mtind_ = NULL;
      mtval_ = NULL;
      vlb_ = NULL;
      vub_ = NULL;
      mrhs_ = NULL;
      msense_ = NULL;
    }
    //cutInfo_=Cgl012Cut(rhs.cutInfo_);
    cutInfo_=Cgl012Cut();
    cutInfo_.setSepGraphSparseThreshold(rhs.cutInfo_.getSepGraphSparseThreshold());
    cutInfo_.setRowMaxPairCount(rhs.cutInfo_.getRowMaxPairCount());
    cutInfo_.setRowMaxFractionalCount(rhs.cutInfo_.getRowMaxFractionalCount());
  }
  return *this;
}
void 
CglZeroHalf::refreshSolver(OsiSolverInterface * solver)
{
  if (!solver||!solver->getNumRows())
    return; // no solver
  delete []  mtbeg_;
  delete []  mtcnt_;
  delete []  mtind_;
  delete []  mtval_;
  delete []  vlb_;
  delete []  vub_;
  delete []  mrhs_;
  delete []  msense_;
  mr_ = 0;
  mc_ = 0;
  mnz_ = 0;
  mtbeg_ = NULL;
  mtcnt_ = NULL;
  mtind_ = NULL;
  mtval_ = NULL;
  vlb_ = NULL;
  vub_ = NULL;
  mrhs_ = NULL;
  msense_ = NULL;
  cutInfo_.free_log_var();
  cutInfo_.free_parity_ilp();
  cutInfo_.free_ilp();
  CoinPackedMatrix rowCopy(*solver->getMatrixByRow());
  const int * column = rowCopy.getIndices();
  const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
  const int * rowLength = rowCopy.getVectorLengths(); 
  const double * rowElements = rowCopy.getElements();
  const double * columnLower = solver->getColLower();
  const double * columnUpper = solver->getColUpper();
  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();
  // Get integer information
  const char * intVar = solver->getColType();
  int iColumn,iRow;
  // count number of possible
  int numberColumns = solver->getNumCols();
  int numberRows = solver->getNumRows();
  vlb_ = new int [numberColumns];
  vub_ = new int [numberColumns];
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    int ilo,iup;
    if (intVar[iColumn]) {
      double lo = columnLower[iColumn];
      if (lo<-COIN_INT_MAX)
	lo=-COIN_INT_MAX;
      ilo= static_cast<int> (ceil(lo));
      double up = columnUpper[iColumn];
      if (up>COIN_INT_MAX)
	up=COIN_INT_MAX;
      iup= static_cast<int> (floor(up));
    } else {
      ilo=COIN_INT_MAX;
      iup=-COIN_INT_MAX;
    }
    vlb_[iColumn]=ilo;
    vub_[iColumn]=iup;
  }
  for (iRow=0;iRow<numberRows;iRow++) {
    int n = rowLength[iRow];
    bool good=(n>0);
    for (CoinBigIndex j=rowStart[iRow];
	 j<rowStart[iRow]+n;j++) {
      int jColumn = column[j];
      if (vlb_[jColumn]==COIN_INT_MAX) {
	// continuous
	good=false;
	break;
      } else {
	double value = rowElements[j];
	if (fabs(value-floor(value+0.5))>1.0e-15||fabs(value)>=COIN_INT_MAX) {
	  // not integer coefficient
	  good=false;
	  break;
	}
      }
    }
    double lo = rowLower[iRow];
    double up = rowUpper[iRow];
    int iType=1;
    double rhs=1.0e20;
    if (lo>-1.0e20) {
      if (fabs(lo-floor(lo+0.5))>1.0e-15) {
	// not integer coefficient
	good=false;
      }
      rhs=fabs(lo);
      if (up<1.0e20) {
	rhs=std::max(fabs(lo),fabs(up));
	if (lo!=up)
	  iType=2; // ranged so make copy
	if (fabs(up-floor(up+0.5))>1.0e-12) {
	  // not integer coefficient
	  good=false;
	}
      }
    } else if (up<1.0e20) {
      rhs=fabs(up);
      if (up<1.0e20) {
	if (fabs(up-floor(up+0.5))>1.0e-12) {
	  // not integer coefficient
	  good=false;
	}
      }
    }
    if (good&&rhs<COIN_INT_MAX) {
      mr_+=iType;
      mnz_ += iType*n;
    }
  }
  int saveMr=mr_;
  int saveMnz=mnz_;
  if (mnz_) {
    mc_ = numberColumns;
    mtbeg_ = new int [mr_];
    mtcnt_ = new int [mr_];
    mtind_ = new int [mnz_];
    mtval_ = new int [mnz_];
    mrhs_ = new int [mr_];
    msense_ = new char [mr_];
    mr_=0;
    mnz_=0;
    for (iRow=0;iRow<numberRows;iRow++) {
      int n = rowLength[iRow];
      bool good=(n>0);
      for (CoinBigIndex j=rowStart[iRow];
	   j<rowStart[iRow]+n;j++) {
	int jColumn = column[j];
	if (vlb_[jColumn]==COIN_INT_MAX) {
	  // continuous
	  good=false;
	  break;
	} else {
	  double value = rowElements[j];
	  if (fabs(value-floor(value+0.5))>1.0e-15||fabs(value)>=COIN_INT_MAX) {
	    // not integer coefficient
	    good=false;
	    break;
	  }
	}
      }
      double lo = rowLower[iRow];
      double up = rowUpper[iRow];
      int iType=1;
      double rhs=1.0e20;
      if (lo>-1.0e20) {
	if (fabs(lo-floor(lo+0.5))>1.0e-15) {
	  // not integer coefficient
	  good=false;
	}
	rhs=fabs(lo);
	if (up<1.0e20) {
	  rhs=std::max(fabs(lo),fabs(up));
	  if (lo!=up)
	    iType=2; // ranged so make copy
	  if (fabs(up-floor(up+0.5))>1.0e-12) {
	    // not integer coefficient
	    good=false;
	  }
	}
      } else if (up<1.0e20) {
	rhs=fabs(up);
	if (up<1.0e20) {
	  if (fabs(up-floor(up+0.5))>1.0e-12) {
	    // not integer coefficient
	    good=false;
	  }
	}
      }
      if (good&&rhs<COIN_INT_MAX) {
	mtbeg_[mr_]=mnz_;
	for (CoinBigIndex j=rowStart[iRow];
	     j<rowStart[iRow]+n;j++) {
	  int jColumn = column[j];
	  double value = rowElements[j];
	  assert (fabs(value)<COIN_INT_MAX);
	  int iValue = static_cast<int> (floor(value+0.5));
	  if (iValue) {
	    mtind_[mnz_]=jColumn;
	    mtval_[mnz_++]=iValue;
	  }
	}
	mtcnt_[mr_]=mnz_-mtbeg_[mr_];
	if (iType==1) {
	  if (lo>-1.0e20) {
	    mrhs_[mr_]=static_cast<int> (floor(lo+0.5));
	    msense_[mr_]='G';
	  } else {
	    mrhs_[mr_]=static_cast<int> (floor(up+0.5));
	    msense_[mr_]='L';
	  }
	  mr_++;
	} else {
	  // ranged!
	  mrhs_[mr_]=static_cast<int> (floor(lo+0.5));
	  msense_[mr_]='G';
	  int k = mnz_-mtbeg_[mr_]; 
	  mr_++;
	  mtbeg_[mr_]=mnz_;
	  memcpy(mtind_+mnz_,mtind_+mnz_-k,k*sizeof(int));
	  memcpy(mtval_+mnz_,mtval_+mnz_-k,k*sizeof(int));
	  mnz_+= mtcnt_[mr_-1];
	  mtcnt_[mr_]=mnz_-mtbeg_[mr_];
	  mrhs_[mr_]=static_cast<int> (floor(up+0.5));
	  msense_[mr_]='L';
	  mr_++;
	}
      }
    }
    assert(saveMr==mr_);
    assert(saveMnz==mnz_);
    cutInfo_.ilp_load(mr_,mc_,mnz_,mtbeg_,mtcnt_,mtind_,mtval_,
	     vlb_,vub_,mrhs_,msense_);
    cutInfo_.alloc_parity_ilp(mr_,mc_,mnz_);
    cutInfo_.initialize_log_var();
  } else {
    // no good
    delete [] vlb_;
    delete [] vub_;
    vlb_ = NULL;
    vub_ = NULL;
    mr_=0;
    mnz_=0;
  }
}
// Create C++ lines to get to current state
std::string
CglZeroHalf::generateCpp( FILE * fp) 
{
  CglZeroHalf other;
  fprintf(fp,"0#include \"CglZeroHalf.hpp\"\n");
  fprintf(fp,"3  CglZeroHalf zeroHalf;\n");
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  zeroHalf.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  zeroHalf.setAggressiveness(%d);\n",getAggressiveness());
  if (getSepGraphSparseThreshold()!=other.getSepGraphSparseThreshold())
    fprintf(fp,"3  zeroHalf.setSepGraphSparseThreshold(%d);\n",getSepGraphSparseThreshold());
  if (getRowMaxPairCount()!=other.getRowMaxPairCount())
    fprintf(fp,"3  zeroHalf.setRowMaxPairCount(%d);\n",getRowMaxPairCount());
  if (getRowMaxFractionalCount()!=other.getRowMaxFractionalCount())
    fprintf(fp,"3  zeroHalf.setRowMaxFractionalCount(%d);\n",getRowMaxFractionalCount());
  return "zeroHalf";
}
#include <vector>
#include <algorithm>
//bool operator() (cgl_node * x, cgl_node * y) {
bool best(cgl_node * x, cgl_node * y) {
  return (x->distanceBack>y->distanceBack);
}
#ifndef CGL_NEW_SHORT
void cglShortestPath(cgl_graph * graph, int source, int maximumLength)
#else
void cglShortestPath(auxiliary_graph * graph, int source, int maximumLength)
#endif
{
  int numberNodes=graph->nnodes;
#define HEAP
#ifndef HEAP
  int * candidate = new int [numberNodes];
#endif
  cgl_node * nodes = graph->nodes;
  int i;
  for ( i=0;i<numberNodes;i++) {
    nodes[i].parentNode=-1;
    nodes[i].distanceBack=COIN_INT_MAX;
#ifndef HEAP
    candidate[i]=i;
#endif
  }
  nodes[source].distanceBack=0;
#ifdef HEAP
  std::vector <cgl_node *> nodes_;
  for ( i=0;i<numberNodes;i++) 
    nodes_.push_back(nodes+i);
  // create heap
  std::make_heap(nodes_.begin(), nodes_.end(), best);
#endif
  int numberCandidates = numberNodes;
  while (numberCandidates>0) {
#ifdef HEAP
  cgl_node * bestNode = nodes_.front();
  int iNode = bestNode->index;
  std::pop_heap(nodes_.begin(), nodes_.end(), best);
  nodes_.pop_back();
  if (nodes[iNode].distanceBack==COIN_INT_MAX)
    break;
  // Prune: no useful odd cycle can have weight > maximumLength
  if (nodes[iNode].distanceBack >= maximumLength)
    break;
  numberCandidates--;
#else
    int best=-1;
    int bestDistance=COIN_INT_MAX;
    for (i=0;i<numberCandidates;i++) {
      int iNode = candidate[i];
      if (nodes[iNode].distanceBack<bestDistance) {
	best=i;
	bestDistance=nodes[iNode].distanceBack;
      }
    }
    if (best<0)
      break; // disconnected graph?
    int iNode = candidate[best];
    numberCandidates--;
    candidate[best]=candidate[numberCandidates];
#endif
    cgl_arc * nextNode = nodes[iNode+1].firstArc;
    int thisDistance = nodes[iNode].distanceBack;
    for (cgl_arc * arc=nodes[iNode].firstArc;arc!=nextNode;arc++) {
      int toNode = arc->to;
      int dist = arc->length;
      if (thisDistance+dist<nodes[toNode].distanceBack) {
	nodes[toNode].distanceBack=thisDistance+dist;
	nodes[toNode].parentNode=iNode;
#ifdef HEAP
	nodes_.push_back(nodes+toNode);
	//printf("size %d\n",nodes_.size());
#endif
      }
    }
  }
}
