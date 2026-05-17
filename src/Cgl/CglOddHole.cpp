// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <iostream>

#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CglOddHole.hpp"
//#define CGL_DEBUG
// We may want to sort cut
typedef struct {double dj;double element; int sequence;} 
double_double_int_triple;
class double_double_int_triple_compare {
public:
  bool operator() (double_double_int_triple x , double_double_int_triple y) const
  {
    return ( x.dj < y.dj);
  }
}; 
//-------------------------------------------------------------------------------
// Generate three cycle cuts
//------------------------------------------------------------------- 
void CglOddHole::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			      const CglTreeInfo info)
{
  // Get basic problem information
  int nRows=si.getNumRows(); 
  int nCols=si.getNumCols(); 
  
  const CoinPackedMatrix * rowCopy = si.getMatrixByRow();

  // Could do cliques and extra OSL cliques
  // For moment just easy ones
  
  // If no information exists then get a list of suitable rows
  // If it does then suitable rows are subset of information
  
  CglOddHole temp;
  int * checkRow = new int[nRows];
  int i;
  if (!suitableRows_) {
    for (i=0;i<nRows;i++) {
      checkRow[i]=1;
    }
  } else {
    // initialize and extend rows to current size
    memset(checkRow,0,nRows*sizeof(int));
    memcpy(checkRow,suitableRows_,std::min(nRows,numberRows_)*sizeof(int));
  }
  temp.createRowList(si,checkRow);
  // now cut down further by only allowing rows with fractional solution
  double * solution = new double[nCols];
  memcpy(solution,si.getColSolution(),nCols*sizeof(double));
  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * collower = si.getColLower();
  const double * colupper = si.getColUpper();
  int * suitable = temp.suitableRows_;

  // At present I am using new and delete as easier to see arrays in debugger
  int * fixed = new int[nCols]; // mark fixed columns 
  const char * intVar = si.getColType();
  for (i=0;i<nCols;i++) {
    if (intVar[i]==1) {
      fixed[i]=0;
      if (colupper[i]-collower[i]<epsilon_) {
	solution[i]=0.0;
	fixed[i]=2;
      } else if (solution[i]<epsilon_) {
	solution[i]=0.0;
	fixed[i]=-1;
      } else if (solution[i]>onetol_) {
	solution[i]=1.0;
	fixed[i]=+1;
      }
    } else {
      //mark as fixed even if not (can not intersect any interesting rows)
      solution[i]=0.0;
      fixed[i]=3;
    }
  }
  // first do packed
  const double * rowlower = si.getRowLower();
  const double * rowupper = si.getRowUpper();
  for (i=0;i<nRows;i++) {
    if (suitable[i]) {
      CoinBigIndex k;
      double sum=0.0;
      if (rowupper[i]>1.001) suitable[i]=-1;
      for (k=rowStart[i]; k<rowStart[i]+rowLength[i];k++) {
	int icol=column[k];
	if (!fixed[icol]) sum += solution[icol];
      }
      if (sum<0.9) suitable[i]=-1; //say no good
    }
  }
#ifdef CGL_DEBUG
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
  if (debugger&&!debugger->onOptimalPath(si))
    debugger = NULL;
#else
  const OsiRowCutDebugger * debugger = NULL;
#endif
  temp.generateCuts(debugger, *rowCopy,solution,
		    si.getReducedCost(),cs,suitable,fixed,info,true);
  // now cover
  //if no >= then skip
  bool doCover=false;
  int nsuitable=0;
  for (i=0;i<nRows;i++) {
    suitable[i]=abs(suitable[i]);
    if (suitable[i]) {
      CoinBigIndex k;
      double sum=0.0;
      if (rowlower[i]<0.999) sum=2.0;
      if (rowupper[i]>1.001) doCover=true;
      for (k=rowStart[i]; k<rowStart[i]+rowLength[i];k++) {
	int icol=column[k];
	if (!fixed[icol]) sum += solution[icol];
	if (fixed[icol]==1) sum=2.0; //don't use if any at 1
      }
      if (sum>1.1) {
	suitable[i]=-1; //say no good
      } else {
	nsuitable++;
      }
    }
  }
  if (doCover&&nsuitable) 
    temp.generateCuts(debugger, *rowCopy,solution,si.getReducedCost(),
		      cs,suitable,fixed,info,false);
  delete [] checkRow;
  delete [] solution;
  delete [] fixed;
    
}
void CglOddHole::generateCuts(const OsiRowCutDebugger * /*debugger*/,
			      const CoinPackedMatrix & rowCopy, 
				 const double * solution, 
			      const double * dj, OsiCuts & cs,
				 const int * suitableRow,
			      const int * fixedColumn,
			      const CglTreeInfo info,
			      bool packed)
{
  CoinPackedMatrix columnCopy = rowCopy;
  columnCopy.reverseOrdering();

  // Get basic problem information
  int nRows=columnCopy.getNumRows(); 
  int nCols=columnCopy.getNumCols(); 
  
  const int * column = rowCopy.getIndices();
  const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
  const int * rowLength = rowCopy.getVectorLengths(); 
  
  const int * row = columnCopy.getIndices();
  const CoinBigIndex * columnStart = columnCopy.getVectorStarts();
  const int * columnLength = columnCopy.getVectorLengths(); 

  // we need only look at suitable rows and variables with unsatisfied 0-1
  // lookup from true row to compressed matrix
  int * mrow = new int[nRows];
  // lookup from true column to compressed
  int * lookup = new int[nCols];
  // number of columns in compressed matrix
  int nSmall=0;
  int i;
  //do lookup from true sequence to compressed
  int n=0;
  for (i=0;i<nRows;i++) {
    if (suitableRow[i]>0) {
      mrow[i]=n++;
    } else {
      mrow[i]=-1;
    }
  }
  for (i=0;i<nCols;i++) {
    if (!fixedColumn[i]) {
      lookup[i]=nSmall++;
    } else {
      lookup[i]=-1;
    }
  }
  int nSmall2=2*nSmall;
  // we don't know how big matrix will be
#define MAXELS 50000
  int maxels=MAXELS;
  //How do I do reallocs in C++?
  // 1.0 - value x(i) - value x(j) for each node pair (or reverse if cover) 
  double * cost = reinterpret_cast<double *> (malloc(maxels*sizeof(double)));
  // arc i.e. j which can be reached from i
  int * to= reinterpret_cast<int *> (malloc(maxels*sizeof(int)));
  //original row for each arc
  int * rowfound=reinterpret_cast<int *> (malloc(maxels*sizeof(int)));
  // start of each column
  int * starts=new int[2*nSmall+1];
  starts[0]=0;
  // useful array for marking if already connected
  int * mark =new int[nSmall2];
  memset(mark,0,nSmall2*sizeof(int));
  n=0; //number of elements in matrix
  for (i=0;i<nCols;i++) {
    int icol=lookup[i];
    if (icol>=0) {
      // column in compressed matrix
      CoinBigIndex k;
      double dd=1.0000001-solution[i];
      mark[icol]=1;
      // reallocate if matrix reached size limit
      if (n+nCols>maxels) {
	maxels*=2;
	cost=reinterpret_cast<double *> (realloc(cost,maxels*sizeof(double)));
	to=reinterpret_cast<int *> (realloc(to,maxels*sizeof(int)));
	rowfound=reinterpret_cast<int *> (realloc(rowfound,maxels*sizeof(int)));
      }
      // get all other connected variables
      for (k=columnStart[i];k<columnStart[i]+columnLength[i];k++) {
	int irow=row[k];
	int jrow=mrow[irow];
	// but only if row in compressed matrix
	if (jrow>=0) {
	  CoinBigIndex j;
	  for (j=rowStart[irow];j<rowStart[irow]+rowLength[irow];j++) {
	    int jcol=column[j];
	    int kcol=lookup[jcol];
	    if (kcol>=0&&!mark[kcol]) {
	      cost[n]=dd-solution[jcol];
	      to[n]=kcol;
	      rowfound[n++]=irow;//original row
	      mark[kcol]=1;
	    }
	  }
	}
      }
      starts[icol+1]=n;
      // zero out markers for next column
      mark[icol]=0;
      for (k=starts[icol];k<starts[icol+1];k++) {
	int ito=to[k];
	if (ito<0||ito>=nSmall) abort();
	mark[to[k]]=0;
      }
    }
  }
  //if cover then change sign - otherwise make sure positive
  if (packed) {
    for (i=0;i<n;i++) {
      if (cost[i]<1.0e-10) {
	cost[i]=1.0e-10;
      }
    }
  } else {
    for (i=0;i<n;i++) {
      cost[i]=-cost[i];
      if (cost[i]<1.0e-10) {
	cost[i]=1.0e-10;
      }
    }
  }
  // we are going to double size 

  if (2*n>maxels) {
    maxels=2*n;
    cost=reinterpret_cast<double *> (realloc(cost,maxels*sizeof(double)));
    to=reinterpret_cast<int *> (realloc(to,maxels*sizeof(int)));
    rowfound=reinterpret_cast<int *> (realloc(rowfound,maxels*sizeof(int)));
  }
  /* copy and make bipartite*/

  for (i=0;i<nSmall;i++) {
    int k,j=i+nSmall;
    for (k=starts[i];k<starts[i+1];k++) {
      int ito=to[k];
      to[n]=ito;
      to[k]=ito+nSmall;
      cost[n]=cost[k];
      rowfound[n++]=rowfound[k];;
    }
    starts[j+1]=n;
  }
  //random numbers to winnow out duplicate cuts
  double * check = new double[nCols];
  if (info.randomNumberGenerator) {
    const CoinThreadRandom * randomGenerator = info.randomNumberGenerator;
    for (i=0;i<nCols;i++) {
      check[i]=randomGenerator->randomDouble();
    }
  } else {
    CoinSeedRandom(13579);
    for (i=0;i<nCols;i++) {
      check[i]=CoinDrand48(); // NOT on a thread by thread basis
    }
  }

  // Shortest path algorithm from Dijkstra - is there a better one?

  typedef struct {
    double cost; //cost to starting node
    int back; //previous node
  } Path;
  typedef struct {
    double cost; //cost to starting node
    int node; //node
  } Item;
  Item * stack = new Item [nSmall2];
  Path * path = new Path [nSmall2];
  // arrays below are used only if looks promising
  // allocate here
  // we don't know how many cuts will be generated
  int ncuts=0;
  int maxcuts=1000;
  double * hash = reinterpret_cast<double *> (malloc(maxcuts*sizeof(double)));
  // to clean (should not be needed)
  int * clean = new int[nSmall2];
  int * candidate = new int[std::max(nSmall2,nCols)];
  double * element = new double[nCols];
  // in case we want to sort
  double_double_int_triple * sortit = 
    new double_double_int_triple [nCols];
  memset(mark,0,nSmall2*sizeof(int));
  int * countcol = new int[nCols];
  memset(countcol,0,nCols*sizeof(int));
  int bias = packed ? 0 : 1; //amount to add before halving
  // If nSmall large then should do a randomized subset
  // Improvement 1
  int icol;
  for (icol=0;icol<nSmall;icol++) {
    int j;
    int jcol=icol+nSmall;
    int istack=1;
    for (j=0;j<nSmall2;j++) {
      path[j].cost=1.0e70;
      path[j].back=nSmall2+1;
    }
    path[icol].cost=0.0;
    path[icol].back=-1;
    stack[0].cost=0.0;
    stack[0].node=icol;
    mark[icol]=1;
    while(istack) {
      Item thisItem=stack[--istack];
      double thisCost=thisItem.cost;
      int inode=thisItem.node;
      int k;
      mark[inode]=0; //say available for further work
      // See if sorting every so many would help (and which way)?
      // Improvement 2
      for (k=starts[inode];k<starts[inode+1];k++) {
	int jnode=to[k];
	if (!mark[jnode]&&thisCost+cost[k]<path[jnode].cost-1.0e-12) {
	  path[jnode].cost=thisCost+cost[k];
	  path[jnode].back=inode;
	  // add to stack
	  stack[istack].cost=path[jnode].cost;
	  stack[istack++].node=jnode;
	  mark[jnode]=1;
#ifdef CGL_DEBUG
	  assert (istack<=nSmall2);
#endif
	}
      }
    }
    bool good=(path[jcol].cost<0.9999);

    if (good)  { /* try */
      int ii;
      int nrow2=0;
      int nclean=0;
      double sum=0;
#ifdef CGL_DEBUG
      printf("** %d ",jcol-nSmall);
#endif
      ii=1;
      candidate[0]=jcol;
      while(jcol!=icol) {
	int jjcol;
	jcol=path[jcol].back;
	if (jcol>=nSmall) {
	  jjcol=jcol-nSmall;
	} else {
	  jjcol=jcol;
	}
#ifdef CGL_DEBUG
	printf(" %d",jjcol);
#endif
	if (mark[jjcol]) {
	  // good=false;
	  // probably means this is from another cycle (will have been found)
	  // one of cycles must be zero cost
	  // printf("variable already on chain!\n");
	} else {
	  mark[jjcol]=1;
	  clean[nclean++]=jjcol;
	  candidate[ii++]=jcol;
#ifdef CGL_DEBUG
	  assert (ii<=nSmall2);
#endif
	}
      }
#ifdef CGL_DEBUG
      printf("\n");
#endif
      for (j=0;j<nclean;j++) {
	int k=clean[j];
	mark[k]=0;
      }
      if (good) {
	int k;
	for (k=ii-1;k>0;k--) {
	  int jk,kk=candidate[k];
	  int ix=0;
	  for (jk=starts[kk];jk<starts[kk+1];jk++) {
	    int ito=to[jk];
	    if (ito==candidate[k-1]) {
	      ix=1;
	      // back to original row
	      mrow[nrow2++]=rowfound[jk];
	      break;
	    }
	  }
	  if (!ix) {
	    good=false;
	  }
	}
	if ((nrow2&1)!=1) {
	  good=false;
	}
	if (good) {
	  int nincut=0;
	  for (k=0;k<nrow2;k++) {
	    CoinBigIndex j;
	    int irow=mrow[k];
	    for (j=rowStart[irow];j<rowStart[irow]+rowLength[irow];j++) {
	      int icol=column[j];
	      if (!countcol[icol]) candidate[nincut++]=icol;
	      countcol[icol]++;
	    }
	  }
#ifdef CGL_DEBUG
	  printf("true constraint %d",nrow2);
#endif
	  nrow2=nrow2>>1;
	  double rhs=nrow2; 
	  if (!packed) rhs++; // +1 for cover
	  ii=0;
	  for (k=0;k<nincut;k++) {
	    int jcol=candidate[k];
	    if (countcol[jcol]) {
#ifdef CGL_DEBUG
	      printf(" %d %d",jcol,countcol[jcol]);
#endif
	      int ihalf=(countcol[jcol]+bias)>>1;
	      if (ihalf) {
		element[ii]=ihalf;
		sum+=solution[jcol]*element[ii];
		/*printf("%d %g %g\n",jcol,element[ii],sumall[jcol]);*/
		candidate[ii++]=jcol;
	      }
	      countcol[jcol]=0;
	    }
	  }
#ifdef CGL_DEBUG
          printf("\n");
#endif
	  OsiRowCut rc;
	  double violation=0.0;
	  if (packed) {
	    violation = sum-rhs;
	    rc.setLb(-COIN_DBL_MAX);
	    rc.setUb(rhs);   
	  } else {
	    // other way for cover
	    violation = rhs-sum;
	    rc.setUb(COIN_DBL_MAX);
	    rc.setLb(rhs);   
	  }
	  if (violation<minimumViolation_) {
#ifdef CGL_DEBUG
	    printf("why no cut\n");
#endif
	    good=false;
	  } else {
	    if (static_cast<double> (ii) * minimumViolationPer_>violation||
		ii>maximumEntries_) {
#ifdef CGL_DEBUG
	      printf("why no cut\n");
#endif
	      if (packed) {
		// sort and see if we can get down to length
		// relax by taking out ones with solution 0.0
		nincut=ii;
		for (k=0;k<nincut;k++) {
		  int jcol=candidate[k];
		  double value = fabs(dj[jcol]);
		  if (solution[jcol])
		  value = -solution[jcol];
		  sortit[k].dj=value;
		  sortit[k].element=element[k];
		  sortit[k].sequence=jcol;
		}
		// sort 
		std::sort(sortit,sortit+nincut,double_double_int_triple_compare());
		nincut = std::min(nincut,maximumEntries_);
		sum=0.0;
		for (k=0;k<nincut;k++) {
		  int jcol=sortit[k].sequence;
		  candidate[k]=jcol;
		  element[k]=sortit[k].element;
		  sum+=solution[jcol]*element[k];
		}
		violation = sum-rhs;
		ii=nincut;
		if (violation<minimumViolation_) {
		  good=false;
		}
	      } else { 
		good=false;
	      }
	    }
	  }
	  if (good) {
	    //this assumes not many cuts
	    int j;
#if 0
	    double value=0.0;
	    for (j=0;j<ii;j++) {
	      int icol=candidate[j];
	      value += check[icol]*element[j];
	    }
#else
            CoinPackedVector candidatePv(ii,candidate,element);
            candidatePv.sortIncrIndex();
            double value = candidatePv.dotProduct(check);
#endif

	    for (j=0;j<ncuts;j++) {
	      if (value==hash[j]) {
		//could check equality - quicker just to assume
		break;
	      }
	    }
	    if (j==ncuts) {
	      //new
	      if (ncuts==maxcuts) {
		maxcuts *= 2;
		hash = reinterpret_cast<double *> (realloc(hash,maxcuts*sizeof(double)));
	      }
	      hash[ncuts++]=value;
	      rc.setRow(ii,candidate,element);
#ifdef CGL_DEBUG
	      printf("sum %g rhs %g %d\n",sum,rhs,ii);
	      if (debugger) 
		assert(!debugger->invalidCut(rc)); 
#endif
	      cs.insertIfNotDuplicate(rc);
	    }
	  }
	}
	/* end of adding cut */
      }
    }
  }
  delete [] countcol;
  delete [] element;
  delete [] candidate;
  delete [] sortit;
  delete [] clean;
  delete [] path;
  delete [] stack;
  free(hash);
  delete [] check;
  delete [] mark;
  delete [] starts;
  delete [] lookup;
  delete [] mrow;
  free(rowfound);
  free(to);
  free(cost);
}

// Create a list of rows which might yield cuts
// The possible parameter is a list to cut down search
void CglOddHole::createRowList( const OsiSolverInterface & si,
		      const int * possible)
{
  // Get basic problem information
  int nRows=si.getNumRows(); 
  
  const CoinPackedMatrix * rowCopy = si.getMatrixByRow();

  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  
  int rowIndex;
  delete [] suitableRows_;
  numberRows_=nRows;

  const double * rowElements = rowCopy->getElements();
  const double * rowupper = si.getRowUpper();
  const double * rowlower = si.getRowLower();
  const double * collower = si.getColLower();
  const double * colupper = si.getColUpper();

  suitableRows_=new int[nRows];
  if (possible) {
    memcpy(suitableRows_,possible,nRows*sizeof(int));
  } else {
    int i;
    for (i=0;i<nRows;i++) {
      suitableRows_[i]=1;
    }
  }
  const char * intVar = si.getColType();
  for (rowIndex=0; rowIndex<nRows; rowIndex++){
    double rhs1=rowupper[rowIndex];
    double rhs2=rowlower[rowIndex];
    if (suitableRows_[rowIndex]) {
      CoinBigIndex i;
      bool goodRow=true;
      for (i=rowStart[rowIndex];
	   i<rowStart[rowIndex]+rowLength[rowIndex];i++) {
	int thisCol=column[i];
	if (colupper[thisCol]-collower[thisCol]>epsilon_) {
	  // could allow general integer variables but unlikely
	  if (intVar[thisCol] !=1) {
	    goodRow=false;
	    break;
	  }
	  if (fabs(rowElements[i]-1.0)>epsilon_) {
	    goodRow=false;
	    break;
	  }
	} else {
	  rhs1 -= collower[thisCol]*rowElements[i];
	  rhs2 -= collower[thisCol]*rowElements[i];
	}
      }
      if (fabs(rhs1-1.0)>epsilon_&&fabs(rhs2-1.0)>epsilon_) {
	goodRow=false;
      }
      if (goodRow) {
	suitableRows_[rowIndex]=1;
      } else {
	suitableRows_[rowIndex]=0;
      }
    }
  }
}
  /// This version passes in a list - 1 marks possible
void CglOddHole::createRowList(int numberRows, const int * whichRow)
{
  suitableRows_=new int [numberRows];
  numberRows_=numberRows;
  memcpy(suitableRows_,whichRow,numberRows*sizeof(int));
}

// Create a list of extra row cliques which may not be in matrix
// At present these are classical cliques
void CglOddHole::createCliqueList(int numberCliques, const int * cliqueStart,
		     const int * cliqueMember)
{
  numberCliques_=numberCliques;
  startClique_=new int[numberCliques_+1];
  memcpy(startClique_,cliqueStart,(numberCliques_+1)*sizeof(int));
  int length=startClique_[numberCliques_];
  member_=new int[length];
  memcpy(member_,cliqueMember,length*sizeof(int));
}
// Returns how many rows might give three cycle cuts
int CglOddHole::numberPossible()
{
  int i,n=0;
  for (i=0;i<numberRows_;i++) {
    if (suitableRows_[i]) n++;
  }
  return n;
}


//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglOddHole::CglOddHole ()
:
CglCutGenerator(),
epsilon_(1.0e-08),
onetol_(1-epsilon_)
{
  // null copy of suitable rows
  numberRows_=0;
  suitableRows_=NULL;
  startClique_=NULL;
  numberCliques_=0;
  member_=NULL;
  minimumViolation_=0.001;
  minimumViolationPer_=0.0003;
  maximumEntries_=100;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglOddHole::CglOddHole (
                                                              const CglOddHole & source)
                                                              :
CglCutGenerator(source),
epsilon_(source.epsilon_),
onetol_(source.onetol_)
{  
  // copy list of suitable rows
  numberRows_=source.numberRows_;
  if (numberRows_) {
    suitableRows_=new int[numberRows_];
    memcpy(suitableRows_,source.suitableRows_,numberRows_*sizeof(int));
  } else {
    suitableRows_=NULL;
  }
  // copy list of cliques
  numberCliques_=source.numberCliques_;
  if (numberCliques_) {
    startClique_=new int[numberCliques_+1];
    memcpy(startClique_,source.startClique_,(numberCliques_+1)*sizeof(int));
    int length=startClique_[numberCliques_];
    member_=new int[length];
    memcpy(member_,source.member_,length*sizeof(int));
  } else {
    startClique_=NULL;
    member_=NULL;
  }
  minimumViolation_=source.minimumViolation_;
  minimumViolationPer_=source.minimumViolationPer_;
  maximumEntries_=source.maximumEntries_;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglOddHole::clone() const
{
  return new CglOddHole(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglOddHole::~CglOddHole ()
{
  // free memory
  delete [] suitableRows_;
  delete [] startClique_;
  delete [] member_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglOddHole &
CglOddHole::operator=(
                                         const CglOddHole& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    epsilon_=rhs.epsilon_;
    onetol_=rhs.onetol_;
    delete [] suitableRows_;
    // copy list of suitable rows
    numberRows_=rhs.numberRows_;
    suitableRows_=new int[numberRows_];
    memcpy(suitableRows_,rhs.suitableRows_,numberRows_*sizeof(int));
    delete [] startClique_;
    delete [] member_;
    // copy list of cliques
    numberCliques_=rhs.numberCliques_;
    if (numberCliques_) {
      startClique_=new int[numberCliques_+1];
      memcpy(startClique_,rhs.startClique_,(numberCliques_+1)*sizeof(int));
      int length=startClique_[numberCliques_];
      member_=new int[length];
      memcpy(member_,rhs.member_,length*sizeof(int));
    } else {
      startClique_=NULL;
      member_=NULL;
    }
    minimumViolation_=rhs.minimumViolation_;
    minimumViolationPer_=rhs.minimumViolationPer_;
    maximumEntries_=rhs.maximumEntries_;
  }
  return *this;
}
// Minimum violation
double 
CglOddHole::getMinimumViolation() const
{
  return minimumViolation_;
}
void 
CglOddHole::setMinimumViolation(double value)
{
  if (value>1.0e-8&&value<=0.5)
    minimumViolation_=value;
}
// Minimum violation per entry
double 
CglOddHole::getMinimumViolationPer() const
{
  return minimumViolationPer_;
}
void 
CglOddHole::setMinimumViolationPer(double value)
{
  if (value>1.0e-8&&value<=0.25)
    minimumViolationPer_=value;
}
// Maximum number of entries in a cut
int 
CglOddHole::getMaximumEntries() const
{
  return maximumEntries_;
}
void 
CglOddHole::setMaximumEntries(int value)
{
  if (value>2)
    maximumEntries_=value;
}

// This can be used to refresh any inforamtion
void 
CglOddHole::refreshSolver(OsiSolverInterface * solver)
{
  // Get integer information
  solver->getColType(true);
}
