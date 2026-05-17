// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>
#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CoinFactorization.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinFinite.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CglTwomir.hpp"
class CoinWarmStartBasis;
#ifdef CBC_HAS_CLP
#define CGL_HAS_CLP_TWOMIR
#endif
#ifdef CGL_HAS_CLP_TWOMIR
#include "OsiClpSolverInterface.hpp"
#endif
#undef DGG_DEBUG_DGG

//#define DGG_DEBUG_DGG 1
//#define CGL_DEBUG
//#define CGL_DEBUG_ZERO
#define  q_max  data->cparams.q_max
#define  q_min  data->cparams.q_min
#define  t_max  data->cparams.t_max
#define  t_min  data->cparams.t_min
#define  a_max  data->cparams.a_max
#define  max_elements  data->cparams.max_elements
//#define CGL_DEBUG
#ifdef CGL_DEBUG
// Declarations and defines for debug build.

#define talk true

namespace {
  const OsiSolverInterface *six ;
}

void write_cut( DGG_constraint_t *cut){ //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       printf("2mir_test: cut: !!!!!!!!!!!!!!!!!!!!!!!***********************************\n");
      for (int i=0; i<cut->nz; i++)
	{ printf(" %12.10f x[%d] ", cut->coeff[i],cut->index[i]);}
      printf(" >= %12.10f \n", cut->rhs); 
}

void testus( DGG_constraint_t *cut){ //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //  const OsiRowCutDebugger debugg(*six, "flugpl.mps"); 
  const OsiRowCutDebugger debugg(*six, "egout.mps"); 
  const OsiRowCutDebugger *debugger = &debugg;
    if (debugger&&!debugger->onOptimalPath(*six))
      return;
   
    OsiRowCut rowcut;

    rowcut.setRow(cut->nz, cut->index, cut->coeff);
    rowcut.setUb(DBL_MAX);
    rowcut.setLb(cut->rhs);
   
    if(debugger->invalidCut(rowcut)){
      write_cut(cut);
      //assert(0); 
    }
}

#else	// CGL_DEBUG

#define talk false

#endif	// CGL_DEBUG

//-------------------------------------------------------------------
// Generate  cuts
//------------------------------------------------------------------- 
void CglTwomir::generateCuts(const OsiSolverInterface & si, OsiCuts & cs, 
			     const CglTreeInfo info )
{
# ifdef CGL_DEBUG
  //!!!!!!!!!!!!!!!!!!
  six = &si;
# endif
  const double * colUpper = si.getColUpper();
  const double * colLower = si.getColLower();
  const OsiSolverInterface * useSolver;
#ifdef CGL_HAS_CLP_TWOMIR
  double * objective = NULL;
  OsiClpSolverInterface * clpSolver = originalSolver_ ? getClpSolver(originalSolver_) : NULL;
  int numberOriginalRows;
  int numberColumns=si.getNumCols();
  int twomirType=0;
  if (!clpSolver) {
#endif
    useSolver=&si;
    // Temp - check if free variables
    int ncol = useSolver->getNumCols();
    int numberFree=0;
    for (int i=0;i<ncol;i++) {
      if (colLower[i]<-1.0e20&&colUpper[i]>1.0e20)
        numberFree++;
    }
    if (numberFree) {
#ifdef COIN_DEVELOP
      if (!info.inTree&&!info.pass)
        printf("CglTwoMir - %d free variables - returning\n",numberFree);
#endif
      return;
    }
#ifdef CGL_HAS_CLP_TWOMIR
  } else {
    useSolver = originalSolver_;
    assert (twomirType_);
    // check simplex is plausible
    if (!clpSolver->getNumRows()||numberColumns!=clpSolver->getNumCols()) {
      delete originalSolver_;
      originalSolver_=si.clone();
      clpSolver = getClpSolver(originalSolver_);
      assert (clpSolver);
      useSolver = originalSolver_;
    }
    ClpSimplex * simplex = clpSolver->getModelPtr();
    numberOriginalRows = simplex->numberRows();
    int numberRows = si.getNumRows();
    assert (numberOriginalRows<=numberRows);
    // only do if different (unless type 2x)
    twomirType = twomirType_%10;
    int whenToDo = twomirType_/10;
    if (whenToDo==2 ||(numberRows>numberOriginalRows && whenToDo==1
		       && (info.options&512)==0) ||
	((info.options&1024)!=0 && (info.options&512)==0)) {
      // bounds
      const double * solution = si.getColSolution();
      memcpy(simplex->columnLower(),colLower,numberColumns*sizeof(double));
      memcpy(simplex->columnUpper(),colUpper,numberColumns*sizeof(double));
      for (int i=0;i<numberColumns;i++) {
	if (colLower[i]<-1.0e20&&colUpper[i]>1.0e20) {
	  double lower=-si.getInfinity();
	  double upper=si.getInfinity();
	  if (solution[i]>0.0)
	    lower=-1.0e10;
	  else
	    upper=1.0e10;
	  originalSolver_->setColLower(i,lower);
	  originalSolver_->setColUpper(i,upper);
	}
      }
      double * obj = simplex->objective();
      objective = CoinCopyOfArray(obj,numberColumns);
      const double * pi = si.getRowPrice();
      const CoinPackedMatrix * rowCopy = si.getMatrixByRow();
      const int * column = rowCopy->getIndices();
      const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
      const int * rowLength = rowCopy->getVectorLengths(); 
      const double * rowElements = rowCopy->getElements();
      const double * rowLower = si.getRowLower();
      const double * rowUpper = si.getRowUpper();
      int numberCopy;
      int numberAdd;
      double * rowLower2 = NULL;
      double * rowUpper2 = NULL;
      int * column2 = NULL;
      CoinBigIndex * rowStart2 = NULL;
      double * rowElements2 = NULL;
      char * copy = new char [numberRows-numberOriginalRows];
      memset(copy,0,numberRows-numberOriginalRows);
      if (twomirType==2) {
	numberCopy=0;
	numberAdd=0;
	for (int iRow=numberOriginalRows;iRow<numberRows;iRow++) {
	  bool simple = true;
	  for (CoinBigIndex k=rowStart[iRow];
	       k<rowStart[iRow]+rowLength[iRow];k++) {
	    double value = rowElements[k];
	    if (value!=floor(value+0.5)) {
	      simple=false;
	      break;
	    }
	  }
	  if (simple) {
	    numberCopy++;
	    numberAdd+=rowLength[iRow];
	    copy[iRow-numberOriginalRows]=1;
	  }
	}
	if (numberCopy) {
	  //printf("Using %d rows out of %d\n",numberCopy,numberRows-numberOriginalRows);
	  rowLower2 = new double [numberCopy];
	  rowUpper2 = new double [numberCopy];
	  rowStart2 = new CoinBigIndex [numberCopy+1];
	  rowStart2[0]=0;
	  column2 = new int [numberAdd];
	  rowElements2 = new double [numberAdd];
	}
      }
      numberCopy=0;
      numberAdd=0;
      //const double * rowSolution = si.getRowActivity();
      //double offset=0.0;
      for (int iRow=numberOriginalRows;iRow<numberRows;iRow++) {
	if (!copy[iRow-numberOriginalRows]) {
	  double value = pi[iRow];
	  //offset += rowSolution[iRow]*value;
	  for (CoinBigIndex k=rowStart[iRow];
	       k<rowStart[iRow]+rowLength[iRow];k++) {
	    int iColumn=column[k];
	    obj[iColumn] -= value*rowElements[k];
	  }
	} else {
	  rowLower2[numberCopy]=rowLower[iRow];
	  rowUpper2[numberCopy]=rowUpper[iRow];
	  for (CoinBigIndex k=rowStart[iRow];
	       k<rowStart[iRow]+rowLength[iRow];k++) {
	    column2[numberAdd]=column[k];
	    rowElements2[numberAdd++]=rowElements[k];
	  }
	  numberCopy++;
	  rowStart2[numberCopy]=numberAdd;
	}
      }
#if 0
      CoinThreadRandom randomNumberGenerator;
      const double * solution = si.getColSolution();
      for (int i=0;i<numberColumns;i++) {
	if (intVar[i]==1) {
	  double randomNumber = randomNumberGenerator.randomDouble();
	  //randomNumber = 0.001*floor(randomNumber*1000.0);
	  if (solution[i]>0.5)
	    obj[i] -= randomNumber*0.001*fabs(obj[i]);
	  else
	    obj[i] += randomNumber*0.001*fabs(obj[i]);
	}
      }
#endif
      if (numberCopy) {
	clpSolver->addRows(numberCopy,
			   rowStart2,column2,rowElements2,
			   rowLower2,rowUpper2);
	delete [] rowLower2 ;
	delete [] rowUpper2 ;
	delete [] column2 ;
	delete [] rowStart2 ;
	delete [] rowElements2 ;
      }
      delete [] copy;
      memcpy(simplex->primalColumnSolution(),si.getColSolution(),
	     numberColumns*sizeof(double));
      CoinWarmStart * warmstart = si.getWarmStart();
      CoinWarmStartBasis* warm =
	dynamic_cast<CoinWarmStartBasis*>(warmstart);
      warm->resize(simplex->numberRows(),numberColumns);
      clpSolver->setBasis(*warm);
      delete warm;
      simplex->setDualObjectiveLimit(COIN_DBL_MAX);  //FIXME should this be si->getInfinity()?
      simplex->setLogLevel(0);
      clpSolver->resolve();
      //printf("Trying - %d its status %d objs %g %g - with offset %g\n",
      //     simplex->numberIterations(),simplex->status(),
      //     simplex->objectiveValue(),si.getObjValue(),simplex->objectiveValue()+offset);
      //simplex->setLogLevel(0);
      if (simplex->status()) {
	//printf("BAD status %d\n",simplex->status());
	//clpSolver->writeMps("clp");
	//si.writeMps("si");
	delete [] objective;
	objective=NULL;
	useSolver=&si;
      }
    }
  }
#endif

  useSolver->getStrParam(OsiProbName,probname_) ;
  int numberRowCutsBefore = cs.sizeRowCuts();
  
  DGG_list_t cut_list;
  DGG_list_init (&cut_list);

  DGG_data_t* data = DGG_getData(reinterpret_cast<const void *> (useSolver));
  // Note that the lhs variables are hash defines to data->cparams.*
  q_max = q_max_;
  q_min = q_min_;
  t_max = t_max_;
  t_min = t_min_;
  a_max = a_max_;
  max_elements = info.inTree ? max_elements_ : max_elements_root_;
  data->gomory_threshold = info.inTree ? away_ : awayAtRoot_;
  if (!info.inTree) {
    //const CoinPackedMatrix * columnCopy = useSolver->getMatrixByCol();
    //int numberColumns=columnCopy->getNumCols(); 
    if (!info.pass||(info.options&32)!=0) {
      max_elements=useSolver->getNumCols();
      //} else {
      //int numberRows=columnCopy.getNumRows();
      //int numberElements=columnCopy->getNumElements();
      //if (max_elements>500&&numberElements>10*numberColumns)
      //max_elements=numberColumns;
    }
  }

  if (!do_mir_) t_max = t_min - 1;
  if (!do_2mir_) q_max = q_min - 1;

  if (do_tab_ && info.level < 1 && info.pass < 6)
    DGG_generateTabRowCuts( &cut_list, data, reinterpret_cast<const void *> (useSolver) );
  
  if (do_form_)
    DGG_generateFormulationCuts( &cut_list, data, reinterpret_cast<const void *> (useSolver),
				 info.formulation_rows,
				 randomNumberGenerator_);
  
#ifdef CGL_DEBUG
  const OsiRowCutDebugger debugg(si,probname_.c_str()) ;
  const OsiRowCutDebugger *debugger = &debugg;
  if (debugger&&!debugger->onOptimalPath(si))
    debugger = NULL;
  else
    {if(talk) printf ("2mir_test: debug success\n");}
#endif
  const char * intVar = si.getColType();
  int i;
  for ( i=0; i<cut_list.n; i++){
    DGG_constraint_t *cut = cut_list.c[i];
    OsiRowCut rowcut;
    if (cut->nz<max_elements) {
      // See if any zero coefficients!!!!!!!
      int nZero=0;
      for( int i=0; i < cut->nz; i++) {
	if (!cut->coeff[i])
	  nZero++;
      }
#ifdef CGL_DEBUG_ZERO
      if (nZero) {
	printf("Cut ");
	for( int i=0; i < cut->nz; i++) {
	  printf("%d %g ",cut->index[i],cut->coeff[i]);
	}
	printf("\n");
      }
#endif
      if (nZero) {
#ifdef CGL_DEBUG_ZERO
	printf("TwoMir cut had %d zero coefficients!\n",nZero);
#endif
      } else {
	//#define CBC_CHECK_CUT
#ifdef CBC_CHECK_CUT
	double rhs = cut->rhs;
	int * cutIndex = cut->index;
	double * packed = cut->coeff;
	int i,number2=cut->nz;
	int number=0;
	double largest=0.0;
	double smallest=1.0e30;
	const double *colUpper = useSolver->getColUpper();
	const double *colLower = useSolver->getColLower();
	bool goodCut=true;
	for (i=0;i<number2;i++) {
	  double value=fabs(packed[i]);
#if 1
	  // add more stringent test (i.e. if close to integer tolerance)
	  if (value<1.0e-6) { 
	    // throw away
	    goodCut=false;
	    break;
#else
	  if (value<5.0e-7) { // was 1.0e-9
	    int iColumn = cutIndex[i];
	    if (colUpper[iColumn]-colLower[iColumn]<100.0) {
	      // weaken cut
	      if (packed[i]>0.0) 
		rhs -= value*colUpper[iColumn];
	      else
	 	rhs += value*colLower[iColumn];
	    } else {
	      // throw away
	      goodCut=false;
	      break;
	    }
#endif
	  } else {
	    int iColumn = cutIndex[i];
	    if (colUpper[iColumn]!=colLower[iColumn]) {
	      largest=std::max(largest,value);
	      smallest=std::min(smallest,value);
	      cutIndex[number]=cutIndex[i];
	      packed[number++]=packed[i];
	    } else {
	      // fixed so subtract out
	      rhs -= packed[i]*colLower[iColumn];
	    }
	  }
	}
	if (largest<1.0e8*smallest&&smallest>1.0e-5&&goodCut) {
	  if (number > 1) {
	    rowcut.setRow(number, cutIndex, packed);
	    rowcut.setUb(si.getInfinity());
	    rowcut.setLb(rhs);
	    cs.insertIfNotDuplicateAndClean(rowcut,61);
	  } else {
	    // singleton row cut!
	    double lb = rhs;
	    double ub = COIN_DBL_MAX;
	    double value = packed[0];
	    int iColumn = cutIndex[0];
	    double lbCol = colLower[iColumn];
	    double ubCol = colUpper[iColumn];
	    bool isInteger =  intVar[iColumn]!=0;
	    // turn lb,ub into new bounds on column
	    if (value>0) {
	      lb = lb/value;
	      if (isInteger)
		lb = ceil(lb-1.0e-4);
	      ub = ubCol;
	    } else {
	      ub = lb/value;
	      if (isInteger)
		ub = floor(ub+1.0e-4);
	      lb = lbCol;
	    }
	    if (lb>ub+1.0e-4) {
	      // infeasible
	      //printf("CUTinf\n");
	      OsiRowCut rc;
	      rc.setRow(0,cutIndex,packed,false);
	      rc.setLb(1.0);
	      rc.setUb(0.0);
	      cs.insertIfNotDuplicateAndClean(rc,62);
	    } else if (lb>lbCol || ub<ubCol) {
	      if (!isInteger) {
		// think
		//printf("CUTnotint\n");
	      } else {
		OsiColCut cc;
		if (lb>lbCol)
		  cc.setLbs(1,&iColumn,&lb);
		if (ub<ubCol)
		  cc.setUbs(1,&iColumn,&ub);
		cs.insert(cc);
		//printf("CUT %g<=%g -> %g<= %g\n",
		//	   colLower[iColumn],colUpper[iColumn],
		//  lb,ub);
	      }
	    } else {
	    //printf("CUT whynone\n");
	    }
	  }
	}
#else
	rowcut.setRow(cut->nz, cut->index, cut->coeff);
	rowcut.setUb(si->getInfinity());
	rowcut.setLb(cut->rhs);
	cs.insertIfNotDuplicateAndClean(rowcut,63);
#endif
      }
    
#ifdef CGL_DEBUG
      if (debugger) {
        if (debugger->invalidCut(rowcut)) {
          write_cut(cut);
          printf ("2mir_test: debug failed, mayday, mayday **********************************\n");} 
	//assert(0);
      }
      //assert(!debugger->invalidCut(rowcut));
#endif
    }
  }
  for ( i=0; i<cut_list.n; i++)
    DGG_freeConstraint (cut_list.c[i]);
  DGG_list_free (&cut_list);
  DGG_freeData (data);
  if (!info.inTree&&((info.options&4)==4||((info.options&8)&&!info.pass))) {
    int numberRowCutsAfter = cs.sizeRowCuts();
    for (int i=numberRowCutsBefore;i<numberRowCutsAfter;i++) {
      int length = cs.rowCutPtr(i)->row().getNumElements();
      if (length<=max_elements_)
	cs.rowCutPtr(i)->setGloballyValid();
    }
  }
#ifdef CGL_HAS_CLP_TWOMIR
  if (objective) {
    int numberRowCutsAfter = cs.sizeRowCuts();
    ClpSimplex * simplex = clpSolver->getModelPtr();
    memcpy(simplex->objective(),objective,numberColumns*sizeof(double));
    delete [] objective;
    // take out locally useless cuts
    const double * solution = si.getColSolution();
    double primalTolerance = 1.0e-7;
    for (int k = numberRowCutsAfter - 1; k >= numberRowCutsBefore; k--) {
      const OsiRowCut * thisCut = cs.rowCutPtr(k) ;
      double sum = 0.0;
      int n = thisCut->row().getNumElements();
      const int * column = thisCut->row().getIndices();
      const double * element = thisCut->row().getElements();
      assert (n);
      for (int i = 0; i < n; i++) {
	double value = element[i];
	sum += value * solution[column[i]];
      }
      if (sum > thisCut->ub() + primalTolerance) {
	sum = sum - thisCut->ub();
      } else if (sum < thisCut->lb() - primalTolerance) {
	sum = thisCut->lb() - sum;
      } else {
	sum = 0.0;
      }
      if (!sum) {
	// take out
	cs.eraseRowCut(k);
      }
    }
#ifdef CLP_INVESTIGATE2
    printf("OR %p pass %d inTree %c - %d cuts (but %d deleted)\n",
       originalSolver_,info.pass,info.inTree?'Y':'N',
       numberRowCutsAfter-numberRowCutsBefore,
       numberRowCutsAfter-cs.sizeRowCuts());
#endif
  }

  int numberRowCutsAfter = cs.sizeRowCuts();
  if (!info.inTree) {
    for (int i=numberRowCutsBefore;i<numberRowCutsAfter;i++) {
      cs.rowCutPtr(i)->setGloballyValid();
    }
  }
  if (twomirType==2) {
    // back to original
    int numberRows = clpSolver->getNumRows();
    if (numberRows>numberOriginalRows) {
      int numberDelete = numberRows-numberOriginalRows;
      int * delRow = new int [numberDelete];
      for (int i=0;i<numberDelete;i++)
	delRow[i]=i+numberOriginalRows;
      clpSolver->deleteRows(numberDelete,delRow);
      delete [] delRow;
    }
  }
#endif
}


//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglTwomir::CglTwomir () :
  CglCutGenerator(),
  probname_(),
  randomNumberGenerator_(987654321),
  away_(0.0005),awayAtRoot_(0.0005),twomirType_(0),
  do_mir_(true), do_2mir_(true), do_tab_(true), do_form_(true),
  t_min_(1), t_max_(1), q_min_(1), q_max_(1), a_max_(2),max_elements_(50000),
  max_elements_root_(50000),form_nrows_(0) {}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglTwomir::CglTwomir (const CglTwomir & source) :
  CglCutGenerator(source),
  randomNumberGenerator_(source.randomNumberGenerator_),
  away_(source.away_),
  awayAtRoot_(source.awayAtRoot_),
  twomirType_(source.twomirType_),
  do_mir_(source.do_mir_),
  do_2mir_(source.do_2mir_),
  do_tab_(source.do_tab_), 
  do_form_(source.do_form_),
  t_min_(source.t_min_),
  t_max_(source.t_max_),
  q_min_(source.q_min_),
  q_max_(source.q_max_),
  a_max_(source.a_max_),
  max_elements_(source.max_elements_),
  max_elements_root_(source.max_elements_root_),
  form_nrows_(source.form_nrows_)
{
  probname_ = source.probname_ ;
  // already done by CglCutGenerator if (source.originalSolver_)
  //originalSolver_ = source.originalSolver_->clone();
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglTwomir::clone() const
{
  return new CglTwomir(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglTwomir::~CglTwomir ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglTwomir &
CglTwomir::operator=(const CglTwomir& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    randomNumberGenerator_ = rhs.randomNumberGenerator_;
    away_=rhs.away_;
    awayAtRoot_=rhs.awayAtRoot_;
    twomirType_ = rhs.twomirType_;
    do_mir_=rhs.do_mir_;
    do_2mir_=rhs.do_2mir_;
    do_tab_=rhs.do_tab_; 
    do_form_=rhs.do_form_;
    t_min_=rhs.t_min_;
    t_max_=rhs.t_max_;
    q_min_=rhs.q_min_;
    q_max_=rhs.q_max_;
    a_max_=rhs.a_max_;
    max_elements_=rhs.max_elements_;
    max_elements_root_ = rhs.max_elements_root_;
    form_nrows_=rhs.form_nrows_;
    delete originalSolver_;
    if (rhs.originalSolver_)
      originalSolver_ = rhs.originalSolver_->clone();
    else
      originalSolver_=NULL;
  }
  return *this;
}
// Pass in a copy of original solver (clone it)
void 
CglTwomir::passInOriginalSolver(OsiSolverInterface * solver)
{
  delete originalSolver_;
  if (solver) {
    if (!twomirType_)
      twomirType_=1;
    originalSolver_ = solver->clone();
    originalSolver_->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
    // Temp - check if free variables
    const double *colUpper = originalSolver_->getColUpper();
    const double *colLower = originalSolver_->getColLower();
    int ncol = originalSolver_->getNumCols();
    int numberFree=0;
    for (int i=0;i<ncol;i++) {
      if (colLower[i]<-1.0e20&&colUpper[i]>1.0e20) 
        numberFree++;
    }
    if (numberFree)
      printf("CglTwoMir - %d free variables - take care\n",numberFree);
  } else {
    twomirType_=0;
    originalSolver_=NULL;
  }
}

int DGG_freeData( DGG_data_t *data )
{
  free(data->info);
  free(data->lb);
  free(data->ub);
  free(data->x);
  free(data->rc);
#ifdef TWOMIR_LESS_MALLOC
  free(data->xtemp_);
  free(data->rctemp_);
  free(data->skalatemp_);
  DGG_freeConstraint(data->constrainttemp0_);
  DGG_freeConstraint(data->constrainttemp1_);
  free(data->pitemp_);
  free(data->rowIsBasictemp_);
  // is row+nrowfree(data->colIsBasictemp_);
  delete data->vector0_;
  delete data->vector1_;
  //free (data->spareArray0_);
#endif

  free(data);
  return 0;
}

DGG_data_t* DGG_getData(const void *osi_ptr )
{
  DGG_data_t *data = NULL;
  const OsiSolverInterface *si = reinterpret_cast<const OsiSolverInterface *> (osi_ptr);

  data = reinterpret_cast<DGG_data_t*> (malloc( sizeof(DGG_data_t)) );

  /* retrieve basis information */
  CoinWarmStart *startbasis = si->getWarmStart();
  const CoinWarmStartBasis *basis = dynamic_cast<const CoinWarmStartBasis*>(startbasis);

  /* retrieve bounds information */
  const double *colUpper = si->getColUpper();
  const double *colLower = si->getColLower();
  const double *rowUpper = si->getRowUpper();
  const double *rowLower = si->getRowLower();
  const double *redCost  = si->getReducedCost();
  const double *dualVal  = si->getRowPrice();

  /* retrieve current optimal solution */
  const double *colSolut = si->getColSolution();

  /* retrieve the matrix in row format */
  const CoinPackedMatrix *rowMatrixPtr = si->getMatrixByRow();
  const int *rowCnt , *rowInd;
  const CoinBigIndex *rowBeg = 0;
  const double *rowMat;
    
  rowBeg = rowMatrixPtr->getVectorStarts();
  rowCnt = rowMatrixPtr->getVectorLengths();
  rowMat = rowMatrixPtr->getElements();
  rowInd = rowMatrixPtr->getIndices();

  /* set number of columns and number of rows */
  data->ncol = si->getNumCols();
  data->nrow = si->getNumRows();
  int total = data->ncol+data->nrow;
  /* set ninteger */
  data->ninteger = 0;

  /* allocate memory for the arrays in 'data' */
  data->info = reinterpret_cast<int*> (malloc( sizeof(int)*total) );
  data->lb = reinterpret_cast<double*> (malloc( sizeof(double)*total) );
  data->ub = reinterpret_cast<double*> (malloc( sizeof(double)*total) );
  data->x  = reinterpret_cast<double*> (malloc( sizeof(double)*total) );
  data->rc = reinterpret_cast<double*> (malloc( sizeof(double)*total) );
#ifdef TWOMIR_LESS_MALLOC
  // get space for scaling etc
  data->xtemp_ = reinterpret_cast<double*> (malloc( sizeof(double)*total));
  data->rctemp_ = reinterpret_cast<double*> (malloc( sizeof(double)*total));
  data->pitemp_ = reinterpret_cast<char*>   (malloc( sizeof(char)*total));
  data->skalatemp_ = reinterpret_cast<int*>   (malloc( sizeof(int)*total));
  // need to work on this one i.e. copy
  data->rowIsBasictemp_ = reinterpret_cast<int*>   (malloc( sizeof(int)*total));
  data->colIsBasictemp_ = data->rowIsBasictemp_+data->nrow; 
  data->constrainttemp0_ = DGG_newConstraint(total);
  data->constrainttemp1_ = DGG_newConstraint(total);
  data->vector0_ = new CoinIndexedVector(data->nrow);
  data->vector1_ = new CoinIndexedVector(data->nrow);
  const char * intVar = si->getColType();
  //data->spareArray0_ = reinterpret_cast<double*> (malloc( sizeof(double)*total));
  //memset(data->spareArray0_,0,sizeof(double)*total);
#endif

  memset(data->info, 0, sizeof(int)*total);

  /* set parameters for column variables */
  data->nbasic_col = 0;

  for(int i=0; i < data->ncol; i++){
    /* is variable basic */
    if ( basis->getStructStatus(i) == CoinWarmStartBasis::basic ){
      data->nbasic_col++;
      DGG_setIsBasic(data,i);
    }

#if DGG_DEBUG_DGG
    {
      int error = 0;

      if ( basis->getStructStatus(i) != CoinWarmStartBasis::basic )
	if ( fabs(colSolut[i] - colUpper[i]) > DGG_BOUND_THRESH )
	  if ( fabs(colSolut[i] - colLower[i]) > DGG_BOUND_THRESH ){
	    fprintf(stdout, "WARNING!!!! : ");
	    fprintf(stdout, "variable %d non-basic, lb = %f, ub = %f, x = %f\n",
		    i, colLower[i], colUpper[i], colSolut[i]);
	    error+=1;
	  }
	
      if (error)
	fprintf(stdout, "\nFOUND %d errors. BYE.\n", error);
    }
#endif

    /* set variable bounds*/
    data->lb[i] = colLower[i];
    data->ub[i] = colUpper[i];


    /* is variable integer */
    if ( intVar[i] ){
      data->ninteger++;
      DGG_setIsInteger(data,i);
      /* tighten variable bounds*/
      data->lb[i] = ceil(colLower[i]);
      data->ub[i] = floor(colUpper[i]);
    }
 
    /* set x value */
    data->x[i] = colSolut[i];

    /* WARNING: remember to set rc!! Its not set!! */
    data->rc[i] = redCost[i];
  }

  /* set parameters for row variables */

  /* slack variables (row variables) work as follows:
     for a ranged constraint, b_dw < ax < b_up, define a variable s so that
     1) if b_up is not infinity:
        ax + s = b_up,   0 < s < b_up - b_dw
     2) if b_up is infinity:
        ax - s = b_dw,   0 < s < b_up - b_dw
  */
  {
    int i,j;
    double activity;

    data->nbasic_row = 0;
    for(i=0, j=data->ncol; i < data->nrow; i++, j++){

      /* check if the row is an equality constraint */ 
      if ( fabs( rowUpper[i] - rowLower[i] ) <= DGG_BOUND_THRESH )
        DGG_setEqualityConstraint(data,j);
      
      /* check if the row is bounded above/below and define variable bounds */
      if ( rowUpper[i] < si->getInfinity() )
        DGG_setIsConstraintBoundedAbove(data,j);
      if ( rowLower[i] > -si->getInfinity() )
        DGG_setIsConstraintBoundedBelow(data,j);

      data->lb[j] = 0.0;
      if (DGG_isConstraintBoundedAbove(data,j) && DGG_isConstraintBoundedBelow(data,j)) 
          data->ub[j] = rowUpper[i] - rowLower[i];
      else
          data->ub[j] = si->getInfinity();

      /* compute row activity. for this we need to go to the row in question,
         and multiply all the coefficients times their respective variables.
         For the moment, we will store the inverse of this value in 
         the 'x' field (since it is in fact a partial computation of it)  */
      {
        CoinBigIndex k;
 
        activity = 0.0;
	for(k=rowBeg[i]; k < rowBeg[i]+rowCnt[i]; k++)
          activity += rowMat[k]*colSolut[rowInd[k]];
      }
  
      /* compute x value */
      if ( DGG_isConstraintBoundedAbove(data,j) )
        data->x[j] = rowUpper[i] - activity;
      else
        data->x[j] = activity - rowLower[i];

      if ( data->x[j] < -DGG_NULL_SLACK ){
#if DGG_DEBUG_DGG
	int k;
	double norm = 0.0, min = DBL_MAX, amin = DBL_MAX, max = DBL_MIN;
        printf("** warning: row %d has negative slack!\n", i);
        for(k=rowBeg[i]; k < rowBeg[i]+rowCnt[i]; k++){
	  norm += rowMat[k]*rowMat[k];
	  if ( fabs(rowMat[k]) < amin ) amin = fabs(rowMat[k]);
	  if ( rowMat[k] < min ) min = rowMat[k];
	  if ( rowMat[k] > max ) max = rowMat[k];
	}
	norm = sqrt(norm);
	printf("min = %f  amin = %f  max = %f\n", min, amin, max);
	printf("rlower = %f activity = %f\n", rowLower[i], activity);
	printf("norm = %f   (b-ax) = %f\n", norm, (rowLower[i] - activity));
	printf("steepn = %f\n", (rowLower[i] - activity)/norm);
#endif
      }

      data->rc[j] = dualVal[i];
#if DGG_DEBUG_SOLVER
      DGG_IF_EXIT( !DGG_isConstraintBoundedAbove(data,j) && !DGG_isConstraintBoundedBelow(data,j),
                   1, "some row is not bounded above or below");
#endif

      /* is variable basic */
      if ( basis->getArtifStatus(i) == CoinWarmStartBasis::basic ){
        data->nbasic_row++;
        DGG_setIsBasic(data,j);
      }

      /* is variable integer. For this we need to go to the row in question,
         and check that the rhs is integer, and that all of the coefficients
         and variables participating in the constraint are also integer.    */
      {
        CoinBigIndex k;
     
        if( DGG_isConstraintBoundedAbove(data,j)) {
          if ( frac_part(rowUpper[i]) > DGG_INTEGRALITY_THRESH )
            goto DONE_ROW; 
	  // need to check this as well
	  if( DGG_isConstraintBoundedBelow(data,j)) {
	    if ( frac_part(rowLower[i]) > DGG_INTEGRALITY_THRESH )
	      goto DONE_ROW;
	  }
        }
        else
          if ( frac_part(rowLower[i]) > DGG_INTEGRALITY_THRESH )
            goto DONE_ROW;
 
        for(k=rowBeg[i]; k < rowBeg[i]+rowCnt[i]; k++)
          if ( frac_part(rowMat[k]) > DGG_INTEGRALITY_THRESH || !DGG_isInteger(data, rowInd[k]))
            goto DONE_ROW;
        
        DGG_setIsInteger(data, j); 
        data->ninteger++;
      }

    DONE_ROW:;
 
      /* set variable bounds: careful!! Later, remember to adjust 
         the INFINITY to a DGG standard (to deal with neq solvers). */
      /* WARNING: remember to set rc!! Its not set!! */
    }
  }

  /* CLEANUP */
  delete basis;
  return data;
}

DGG_constraint_t*
DGG_getSlackExpression(const void *osi_ptr, DGG_data_t* data, int row_index)
{
  DGG_constraint_t *row = 0;
  CoinBigIndex i;
  int j;

  /* retrieve the matrix in row format */
  const OsiSolverInterface *si = reinterpret_cast<const OsiSolverInterface *> (osi_ptr);
  const CoinPackedMatrix *rowMatrixPtr = si->getMatrixByRow();
  const int *rowCnt , *rowInd ;
  const CoinBigIndex *rowBeg;
  const double *rowMat;
  const double *rowUpper;
  const double *rowLower;
    
#ifndef TWOMIR_LESS_MALLOC
  row = DGG_newConstraint(data->ncol);
#else
  row = DGG_newConstraint(data,1);
#endif

  rowBeg = rowMatrixPtr->getVectorStarts();
  rowCnt = rowMatrixPtr->getVectorLengths();
  rowMat = rowMatrixPtr->getElements();
  rowInd = rowMatrixPtr->getIndices();

  rowUpper = si->getRowUpper();
  rowLower = si->getRowLower();

#if DGG_DEBUG_DGG
  if ( row_index < 0 || row_index > data->nrow )
    DGG_THROW(0, "bad row index"); 
#endif

  /* copy the information into the row ADT */
  row->nz = rowCnt[row_index];
  for(j=0, i=rowBeg[row_index]; i < rowBeg[row_index]+rowCnt[row_index]; i++, j++){
    row->coeff[j] = rowMat[i];
    row->index[j] = rowInd[i];
    if (DGG_isConstraintBoundedAbove (data, data->ncol + row_index))
      row->coeff[j] = -row->coeff[j];
  }
  
  row->sense = '?';
  if ( DGG_isConstraintBoundedAbove(data, data->ncol + row_index) )
    row->rhs = rowUpper[row_index];
  else 
    row->rhs = -rowLower[row_index];

  return row;
}

int
DGG_getTableauConstraint( int index,  const void *osi_ptr, DGG_data_t *data,
                          DGG_constraint_t* tabrow, 
                          const int * colIsBasic,
                          const int * /*rowIsBasic*/,
                          CoinFactorization & factorization,
                          int mode )
{

#if DGG_DEBUG_DGG
  /* ensure that the index corresponds to a basic variable */
  if ( !DGG_isBasic(data, index) )
     DGG_THROW(1, "index is non-basic");

  /* ensure that the index corresponds to a column variable */
  if ( index < 0 || index > (data->ncol - 1) )
    DGG_THROW(1, "index not a column variable");
#endif

  /* obtain pointer to solver interface */
  const OsiSolverInterface *si = reinterpret_cast<const OsiSolverInterface *> (osi_ptr);
  DGG_TEST(!si, 1, "null OsiSolverInterfave");

  /* obtain address of the LP matrix */
  const CoinPackedMatrix *colMatrixPtr = si->getMatrixByCol();
  const CoinBigIndex* colBeg = colMatrixPtr->getVectorStarts();
  const int* colCnt = colMatrixPtr->getVectorLengths();
  const int* colInd = colMatrixPtr->getIndices();
  const double* colMat = colMatrixPtr->getElements();

  /* obtain row right-hand-sides */
  const double *rowUpper = si->getRowUpper();
  const double *rowLower = si->getRowLower();

  /* allocate memory for constraint in non-sparse form */
  double *value = NULL, rhs = 0.0;
#ifndef TWOMIR_LESS_MALLOC
  value = reinterpret_cast<double*>(malloc(sizeof(double)*(data->nrow+data->ncol)));
  memset(value, 0, sizeof(double)*(data->nrow+data->ncol));
#endif


  /* obtain the tableau row coefficients for all variables */
  int nz = 0;
  int nz1;
  {
    int i, j, cnt = 0;
    double one = 1.0;
    CoinIndexedVector *work = data->vector0_;
    CoinIndexedVector *array = data->vector1_;

    array->setVector(1,&colIsBasic[index],&one);
 
    factorization.updateColumnTranspose ( work, array );

    int * arrayRows = array->getIndices();
    double *arrayElements = array->denseVector();
    cnt = array->getNumElements();
    double * coeff = tabrow->coeff;
    int * indices = tabrow->index;
    // pivot
    //coeff[0] = 1.0;
    //indices[0] = index;
    int ncol = data->ncol;
    /* compute column (structural) variable coefficients */
    for(j = 0; j < ncol; j++) {
      if ( !DGG_isBasic(data, j) || j==index) {
	double value = 0.0;
	for(CoinBigIndex i=colBeg[j]; i < colBeg[j]+colCnt[j]; i++)
	  value += colMat[i]*arrayElements[ colInd[i] ];
	if ( fabs(value) > DGG_MIN_TABLEAU_COEFFICIENT ) {
	  coeff[nz] = value;
	  indices[nz] = j;
	  nz++;
	}
      }
    }
    /* compute row variable (slack/logical) variable coefficients */
    /* and compute rhs */
    rhs = 0.0;
    nz1 = nz;
    for(j = 0; j < cnt; j++){
      int iRow = arrayRows[j];
      int iSeq = ncol+iRow;
      double value = arrayElements[iRow];
      arrayElements[iRow] = 0.0;
      if ( DGG_isEqualityConstraint(data,iSeq) && !mode ) {
        //value[ iSeq ] = 0.0;
      } else if ( DGG_isConstraintBoundedAbove(data, iSeq) ) {
        //value[ iSeq ] = arrayElements[ iRow ];
	coeff[nz] = value;
	indices[nz++] = iSeq;
      } else {
        //value[ iSeq ] = -1*arrayElements[ iRow ];
	coeff[nz] = -value;
	indices[nz++] = iSeq;
      }
      if ( DGG_isConstraintBoundedAbove(data,iSeq) )
        rhs += value*rowUpper[iRow];
      else
        rhs += value*rowLower[iRow];
    }
    // say array clear
    array->setNumElements(0);
  }


  tabrow->nz = nz;
  CoinSort_2(tabrow->index+nz1,tabrow->index+nz,tabrow->coeff+nz1);
#if 0
  printf("NN index %d nz %d rhs %g\n",index,nz,rhs);
  for (int i=0;i<nz;i++)
    printf("%d %g\n",tabrow->index[i],tabrow->coeff[i]);
#endif
  tabrow->sense = 'E';
  tabrow->rhs = rhs;

#ifndef TWOMIR_LESS_MALLOC
  /* CLEANUP */
  free(value);
#endif

  return 0;
}

int
DGG_getFormulaConstraint( int da_row,  
                                  const void *osi_ptr,   
				  DGG_data_t *data, 
                                  DGG_constraint_t* form_row )
{
  /* ensure that the da_row corresponds to a row */
  if ( data->nrow <= da_row || 0> da_row)      DGG_THROW(1, "row out of range...");
  
  /* obtain pointer to solver interface */
  const OsiSolverInterface *si = reinterpret_cast<const OsiSolverInterface *> (osi_ptr);
  //DGG_TEST(!si, 1, "null OsiSolverInterfave");
  
  /* obtain address of the LP matrix */
  const CoinPackedMatrix *rowMatrixPtr = si->getMatrixByRow();
  const CoinBigIndex* rowBeg = rowMatrixPtr->getVectorStarts();
  const int* rowCnt = rowMatrixPtr->getVectorLengths();
  const int* rowInd = rowMatrixPtr->getIndices();
  const double* rowMat = rowMatrixPtr->getElements();

  /* obtain row right-hand-sides */
  const double *rowUpper = si->getRowUpper();
  const double *rowLower = si->getRowLower();
 
  int nz = rowCnt[da_row]; 

  form_row->nz = nz; 
  form_row->max_nz = nz+1;
 
  int i;
  for( i=0; i < nz; i++) form_row->coeff[i] = rowMat[rowBeg[da_row]+i];
  for( i=0; i < nz; i++) form_row->index[i] = rowInd[rowBeg[da_row]+i];

  if ( DGG_isConstraintBoundedAbove(data,data->ncol + da_row) ){
    form_row->rhs = rowUpper[da_row];
    form_row->sense = 'L';
  }
  else{
    form_row->rhs = rowLower[da_row];
    form_row->sense = 'G';
  }
  if ( DGG_isEqualityConstraint(data,data->ncol + da_row)  )
    form_row->sense = 'E';

  /* now add slack/surplus if there is one */
  if ( DGG_isEqualityConstraint(data,data->ncol + da_row) == 0 ){
    form_row->index[nz] =  data->ncol + da_row; 
    if ( DGG_isConstraintBoundedAbove(data, data->ncol + da_row) )
      form_row->coeff[nz] = 1; 
    else
      form_row->coeff[nz] = -1; 
    form_row->nz +=1;
  }

  return 0;

}
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------

/******************** CONSTRAINT ADTs *****************************************/
DGG_constraint_t* DGG_newConstraint(int max_arrays)
{
   DGG_constraint_t *c = NULL;
   
   if (max_arrays <= 0) return NULL;
   c = reinterpret_cast<DGG_constraint_t*> (malloc(sizeof(DGG_constraint_t)));
   c->nz = 0;
   c->max_nz = max_arrays;
   c->rhs = 0.0;
   c->sense = '?';

   c->coeff = NULL;
   c->index = NULL;
   c->coeff = reinterpret_cast<double*>(malloc(sizeof(double)*max_arrays));
   c->index = reinterpret_cast<int*>(malloc(sizeof(int)*max_arrays));
   return c;
}

void DGG_freeConstraint(DGG_constraint_t *c)
{
  if (c == NULL) return;
  if (c->coeff) free(c->coeff);
  if (c->index) free(c->index);
  free(c);
}

DGG_constraint_t *DGG_copyConstraint(DGG_constraint_t* c)
{
  DGG_constraint_t *nc = NULL;

  if (!c || c->max_nz <= 0) return nc;
  nc = DGG_newConstraint(c->max_nz);
  if (nc == NULL) return nc;

  nc->nz = c->nz;
  nc->rhs = c->rhs;
  nc->sense = c->sense;

  memcpy(nc->coeff, c->coeff, sizeof(double)*nc->nz);
  memcpy(nc->index, c->index, sizeof(int)*nc->nz);

  return nc;
}
#ifdef TWOMIR_LESS_MALLOC
DGG_constraint_t* DGG_newConstraint(DGG_data_t *data, int which) 
{
  DGG_constraint_t *c = which ? data->constrainttemp1_ : data->constrainttemp0_;
  c->nz = 0;
  // c->max_nz = max_arrays;
  c->rhs = 0.0;
  c->sense = '?';

  return c;
}

DGG_constraint_t *DGG_copyConstraint(DGG_constraint_t* c,
				     DGG_data_t *data, int which)
{

  if (!c)
    return NULL; 
  DGG_constraint_t *nc = which ? data->constrainttemp1_ : data->constrainttemp0_;
  nc->nz = c->nz;
  nc->rhs = c->rhs;
  nc->sense = c->sense;

  memcpy(nc->coeff, c->coeff, sizeof(double)*nc->nz);
  memcpy(nc->index, c->index, sizeof(int)*nc->nz);

  return nc;
}
#endif
void DGG_scaleConstraint(DGG_constraint_t *c, int t)
{
  int i;

  c->rhs *= t;
  if (t < 0){
    if (c->sense == 'G') c->sense = 'L';
    else if (c->sense == 'L') c->sense = 'G';
  }
  for(i=0; i<c->nz; i++)  c->coeff[i] *= t;
}

void DGG_list_init (DGG_list_t *l)
{
  l->n = 0;
  l->c = NULL;
  l->ctype = NULL;
  l->alpha = NULL;
}

void DGG_list_free(DGG_list_t *l)
{
  if (l->c != NULL) free (l->c);
  if (l->ctype != NULL) free (l->ctype);
  if (l->alpha != NULL) free (l->alpha);
}

int DGG_list_addcut (DGG_list_t *l, DGG_constraint_t *cut, int ctype, double alpha)
{
  l->n ++;
  l->c = reinterpret_cast<DGG_constraint_t **>(realloc (l->c, l->n * sizeof(DGG_constraint_t *)));
  l->ctype = reinterpret_cast<int *>(realloc (l->ctype, l->n * sizeof (int)));
  l->alpha = reinterpret_cast<double *>(realloc (l->alpha, l->n * sizeof (double)));

  if (l->c == NULL || l->ctype == NULL || l->alpha == NULL){
    printf ("No memory, bailing out\n");
    return -1;
  }

  l->c[l->n - 1] = cut;
  l->ctype[l->n - 1] = ctype;
  l->alpha[l->n - 1] = alpha;
  return 0;
}

void DGG_list_delcut (DGG_list_t *l, int i)
{
  if (i >= l->n && i < 0) return;

  DGG_freeConstraint (l->c[i]);
  l->c[i] = l->c[l->n - 1];
  l->ctype[i] = l->ctype[l->n - 1];
  l->alpha[i] = l->alpha[l->n - 1];
  l->n --;
}

/******************* CONSTRAINT MANIPULATION **********************************/
/* VARIABLES CLOSE TO UPPER BOUNDS:
 we will substitute: x' = (u - x); hence the constraint will change from ax ~ b
 to -ax' ~ b - au note: the new bounds of x' will be, 0 <= x' <= u - l
 VARIABLES  CLOSE TO LOWER BOUNDS:
 we will substitute: x' = (x - l); hence, the constraint will change from
 ax ~ b to ax' ~ b - al. note: some variable lower bounds may have changed
 when doing the complement  in the previous stage - this must be taken into
 account. note: the new bounds  of x' will be, 0 <= x'  <= u - l */

int DGG_transformConstraint( DGG_data_t *data,
                             double **x_out,
			     double **rc_out,
                             char **isint_out,
                             DGG_constraint_t *constraint )
{
  double half;
#ifndef TWOMIR_LESS_MALLOC
  double *px = reinterpret_cast<double*> (malloc( sizeof(double)*constraint->max_nz ));
  double *rc = reinterpret_cast<double*> (malloc( sizeof(double)*constraint->max_nz ));
  char   *pi = reinterpret_cast<char*>   (malloc( sizeof(char)  *constraint->max_nz ));
#else
  double *px = data->xtemp_;
  double *rc = data->rctemp_;
  char   *pi = data->pitemp_;
#endif

  {
    int i, idx;

    for(i=0; i < constraint->nz; i++){
      idx = constraint->index[i];

      px[i] = data->x[idx];
      rc[i] = data->rc[idx];
      pi[i] = static_cast<char>(DGG_isInteger(data, idx)); 
      half = (data->ub[idx] - data->lb[idx]) / 2;

      if ( data->ub[idx] - data->x[idx] < half ){
	px[i] = data->ub[idx] - data->x[idx];
	if (fabs(px[i]) <= DGG_BOUND_THRESH)
	  px[i] = 0.0;
	constraint->rhs -= constraint->coeff[i]*data->ub[idx];
	constraint->coeff[i] *= -1;
      }
      else {
	px[i] = data->x[idx] - data->lb[idx];
	if (fabs(px[i]) <= DGG_BOUND_THRESH)
	  px[i] = 0.0;
	constraint->rhs -= constraint->coeff[i]*data->lb[idx];
      }
    }
  }

  *x_out = px;
  *rc_out = rc;
  *isint_out = pi;

#if DGG_DEBUG_DGG
  DGG_TEST(DGG_isConstraintViolated(data, constraint), 1, "bad transformation");
#endif
 
  return 0;
}

int DGG_unTransformConstraint( DGG_data_t *data, 
                               DGG_constraint_t *constraint )
{
  int i, idx;
  double half;
  
  for(i=0; i < constraint->nz; i++){
    idx = constraint->index[i];
    half = (data->ub[idx] - data->lb[idx]) / 2;
    
    if ( data->ub[idx] - data->x[idx] < half ){
      constraint->rhs -= constraint->coeff[i]*data->ub[idx];	
      constraint->coeff[i] *= -1;
    }
    else 
      constraint->rhs += constraint->coeff[i]*data->lb[idx];
  }
  return 0;
}

int
DGG_substituteSlacks( const void *solver_ptr, 
                          DGG_data_t *data, 
                          DGG_constraint_t *cut )
{
  int i,j, lnz;
  double *lcut, lrhs;
  DGG_constraint_t *row=NULL;
 
  /* lcut will store all the column coefficients. allocate space and init. */
  lcut = reinterpret_cast<double*>(malloc(sizeof(double)*data->ncol)); 
  memset(lcut, 0, sizeof(double)*data->ncol);
 
  /* initialize lrhs */
  lrhs = cut->rhs;

  /* set coefficients in lcut */
  /* technical: we could speed this up by re-using allocated memory 
     for row->coeff and row->index                                  */
  for(i=0; i < cut->nz; i++){
    if ( cut->index[i] < data->ncol )
      lcut[ cut->index[i] ] += cut->coeff[i];
    else{
      row = DGG_getSlackExpression(solver_ptr, data, (cut->index[i] - data->ncol));
      
      for(j=0; j < row->nz; j++)
	lcut[ row->index[j] ] += row->coeff[j]*cut->coeff[i];
      lrhs -= row->rhs*cut->coeff[i];
#ifndef TWOMIR_LESS_MALLOC
      DGG_freeConstraint(row);
#endif
    }
  }

  /* count nz in new constraint */
  lnz = 0;
  for(i=0; i < data->ncol; i++)
    if ( fabs(lcut[i]) > DGG_MIN_TABLEAU_COEFFICIENT )
      lnz += 1;

  /* free row->coeff and row->index, and re-allocate */
  free(cut->coeff); cut->coeff = 0;
  free(cut->index); cut->index = 0;

  cut->nz = lnz;
  cut->max_nz = lnz;
  if (lnz)
    {
      cut->coeff = reinterpret_cast<double*> (malloc( sizeof(double)*lnz ));
      cut->index = reinterpret_cast<int*> (malloc( sizeof(int)*lnz ));
    }

  /* set new constraint */
  lnz = 0;
  for(i=0; i < data->ncol; i++){
    if ( fabs(lcut[i]) > DGG_MIN_TABLEAU_COEFFICIENT ){
      cut->coeff[lnz] = lcut[i];
      cut->index[lnz] = i;
      lnz += 1;
    }
  }
  cut->rhs = lrhs;

  free(lcut);
  return 0; 
}

int DGG_nicefyConstraint( const void * /*solver_ptr*/, 
                          DGG_data_t *data,
			  DGG_constraint_t *cut)
													
{
  
  double min_coef = COIN_DBL_MAX, max_coef = COIN_DBL_MIN;
  
  DGG_TEST(cut->sense == 'L', 1, "can't nicefy an L constraint");
  
  int i;
  for( i=0; i<cut->nz; i++) // first clean out noise
    if( fabs(cut->coeff[i]) < DGG_NICEFY_MIN_ABSVALUE)
      cut->coeff[i] = 0;

  for( i=0; i<cut->nz; i++){
    
    if( DGG_isInteger(data, cut->index[i])){// look at integral vars.

      double aht = ABOV(cut->coeff[i]);
      double ub  = data->ub[ cut->index[i]];

      if(aht <  DGG_NICEFY_MIN_FIX){// coefficient = integer + epsylon
	
	cut->coeff[i] = floor( cut->coeff[i]);
	double ahtu = aht * ub;
	
	if(ahtu<DGG_NICEFY_MAX_PADDING)
	  cut->rhs -= ahtu;// safely remove the fractional part
       	else 
	  cut->coeff[i] += DGG_NICEFY_MIN_FIX; // inflate the fractional part
      }
      else 
	if (1-aht <  DGG_NICEFY_MIN_FIX) // coefficient = integer - epsylon
	  cut->coeff[i] = ceil( cut->coeff[i]);
      
    }// done with integers
    else // now look at continuous variables
      if ( cut->coeff[i] < DGG_NICEFY_MIN_ABSVALUE) // delete all negative and noise
	cut->coeff[i] = 0.0;
      else 
	if(cut->coeff[i] <  DGG_NICEFY_MIN_FIX) {// coefficient = epsylon
	  double au = cut->coeff[i] * data->ub[ cut->index[i]];
	
	  if(au<DGG_NICEFY_MAX_PADDING){ // safely remove the variable
	    cut->coeff[i] = 0.0;
	    cut->rhs -= au;
	  }
	  else 
	    cut->coeff[i] = DGG_NICEFY_MIN_FIX; // inflate the coefficient
	}// done with continuous variables too

    double abs_coef = fabs(cut->coeff[i]);
    min_coef = DGG_MIN(min_coef, abs_coef);
    max_coef = DGG_MAX(max_coef, abs_coef);
  }

  cut->sense = 'G';
  /*
  if ( max_coef > DGG_NICEFY_MAX_RATIO*min_coef ) // kill the cut if numbers are all over the place
    cut->nz = 0;
  */
  return 0;

}

/******************* CUT GENERATION *******************************************/
int
DGG_generateTabRowCuts( DGG_list_t *cut_list,
			    DGG_data_t *data,
			    const void *solver_ptr )

{
  int k, rval = 0;
  DGG_constraint_t *base = NULL;
  int nc = cut_list->n;

  base = DGG_newConstraint(data->ncol + data->nrow);

  if(talk) printf ("2mir_test: generating tab row cuts\n");
  /* allocate memory for basic column/row indicators */
  int *rowIsBasic = 0, *colIsBasic = 0;
#ifndef TWOMIR_LESS_MALLOC
  rowIsBasic = reinterpret_cast<int*>(malloc(sizeof(int)*data->nrow));
  colIsBasic = reinterpret_cast<int*>(malloc(sizeof(int)*data->ncol));
#else
  rowIsBasic = data->rowIsBasictemp_;
  colIsBasic = data->colIsBasictemp_;
#endif
    
  /* initialize the IsBasic arrays with -1 / 1 values indicating 
     where the basic rows and columns are. NOTE: WE could do this 
     only once and keep it in osi_data at the expense of space!! */

  int i;
  for( i=0; i<data->ncol; i++){
    if ( DGG_isBasic(data,i) ) colIsBasic[i] = 1;
    else                       colIsBasic[i] = -1;
  }
  for( i=0; i<data->nrow; i++){
    if ( DGG_isBasic(data,i+data->ncol) ) rowIsBasic[i] = 1;
    else                                  rowIsBasic[i] = -1;
  }

  /* obtain factorization */
  CoinFactorization factorization;
  /* obtain address of the LP matrix */
  const OsiSolverInterface *si = reinterpret_cast<const OsiSolverInterface *> (solver_ptr);
  const CoinPackedMatrix *colMatrixPtr = si->getMatrixByCol();
  rval = factorization.factorize(*colMatrixPtr, rowIsBasic, colIsBasic); 
  /* 0 = okay. -1 = singular. -2 = too many in basis. -99 = memory. */
  DGG_TEST2(rval, 1, "factorization error = %d", rval);

  for(k=0; k<data->ncol; k++){
    if (!(DGG_isBasic(data, k) && DGG_isInteger(data,k))) continue;

    double frac = frac_part (data->x[k]);
    if (frac < data->gomory_threshold || frac > 1-data->gomory_threshold) continue;

    base->nz = 0;
    rval = DGG_getTableauConstraint(k, solver_ptr, data, base, 
                                    colIsBasic,rowIsBasic,factorization,0);
    DGG_CHECKRVAL(rval, rval);

    if (base->nz == 0){
      printf ("2mir_test: why does constraint not exist ?\n");
      continue;
    }

    if (base->nz > 500) continue;
    rval = DGG_generateCutsFromBase(base, cut_list, data, solver_ptr);
    DGG_CHECKRVAL(rval, rval);
  }

#ifndef TWOMIR_LESS_MALLOC
   free(rowIsBasic);
  free(colIsBasic);
#endif

  if(talk)
    printf ("2mir_test: generated %d tab cuts\n", cut_list->n - nc);
  DGG_freeConstraint(base);
  return rval;
}

int DGG_generateFormulationCuts( DGG_list_t *cut_list,
				 DGG_data_t *data,
				 const void *solver_ptr,
				 int nrows,
				 CoinThreadRandom & generator)
{
  int k, rval = 0;
  DGG_constraint_t *base = NULL;
  int num_rows = (data->nrow < nrows) ? data->nrow : nrows;
  int nc = cut_list->n;

  base = DGG_newConstraint(data->ncol + data->nrow);

  if(talk)
    printf ("2mir_test: generating form row cuts %d\n", num_rows);
  for(k=0; k<num_rows; k++) {
    base->nz = 0;

    rval = DGG_getFormulaConstraint(k, solver_ptr, data, base);
    DGG_CHECKRVAL1(rval, rval);

    //printf ("generating formulation for row %d\n", k);
    rval = DGG_generateFormulationCutsFromBase(base, data->x[data->ncol+k],
					       cut_list, data, solver_ptr,
					       generator);
    DGG_CHECKRVAL1(rval, rval);
    if (base->nz == 0){
#ifdef COIN_DEVELOP
      printf ("why does constraint not exist ?\n");
#endif
      continue;
    }
  }

 CLEANUP:
  if(talk)
    printf ("2mir_test: generated %d form cuts\n", cut_list->n - nc);
  DGG_freeConstraint(base);
  return rval;
}


int DGG_generateFormulationCutsFromBase( DGG_constraint_t *base,
					 double slack,
					 DGG_list_t *cut_list,
					 DGG_data_t *data,
					 const void *solver_ptr,
					 CoinThreadRandom & generator)
{
  int i, p, rval;
  int int_skala;
  double skala;
  int num_inlist = 0;
#ifndef TWOMIR_LESS_MALLOC
  int* skala_list = reinterpret_cast<int*> (malloc( sizeof(int)*base->nz ));
#else
  int* skala_list = data->skalatemp_;
#endif
  char *isint = NULL;
  double *xout = NULL, *rcout = NULL;
  DGG_constraint_t *scaled_base = NULL;
  int tot_int = 0;
  double prob_choose = 0.0;
  rval = DGG_transformConstraint(data, &xout, &rcout, &isint, base);
  DGG_CHECKRVAL1(rval, rval);

  for(p = 0; p < base->nz; p++)  if(isint[p]) tot_int ++;
  if (tot_int == 0) goto CLEANUP;

  prob_choose = 5.0/tot_int;

  for(p = 0; p < base->nz; p++) {
    if(isint[p]) if(generator.randomDouble() < prob_choose){
      if(xout[p]<0.01) continue;

      skala =fabs(base->coeff[p]);
      if(skala<0.01)  continue;

      // check if slack is too large
      if (fabs(slack/skala) > 0.5) continue;

#ifndef TWOMIR_LESS_MALLOC
      scaled_base = DGG_copyConstraint(base);
#else
      // need to think and write new copy
      //scaled_base = DGG_copyConstraint(base);
      scaled_base = DGG_copyConstraint(base,data,0);
#endif
      DGG_CHECKRVAL1((scaled_base == NULL),-1);

      if(base->sense == 'L') {
	skala = -skala; 
	scaled_base->sense = 'G';
      }
      
      int_skala = int(100*skala);

      for(i = 0; i< num_inlist; i++) 
	if(int_skala == skala_list[i]) 
	  goto END_LOOP;

      skala_list[num_inlist++] = int_skala;

      scaled_base->rhs = base->rhs/skala;
      for(i = 0; i<base->nz; i++) 
	scaled_base->coeff[i] = base->coeff[i] / skala;
 
      rval = DGG_unTransformConstraint(data, scaled_base);
      DGG_CHECKRVAL1(rval, rval);

      rval = DGG_generateCutsFromBase(scaled_base, cut_list,
				      data, solver_ptr);
      DGG_CHECKRVAL1(rval, rval);
      
    END_LOOP:
#ifndef TWOMIR_LESS_MALLOC
      DGG_freeConstraint(scaled_base);
#endif
      scaled_base = NULL;
    }      
  }

 CLEANUP:
#ifndef TWOMIR_LESS_MALLOC
  if (isint) free(isint);    
  if (xout)  free(xout);   
  if (rcout)  free(rcout);
  if (skala_list) free(skala_list);
  if (scaled_base != NULL) DGG_freeConstraint (scaled_base);
#endif
  return rval;
}

int
DGG_generateCutsFromBase( DGG_constraint_t *orig_base,
			  DGG_list_t *cut_list,
			  DGG_data_t *data,
			  const void *solver_ptr )
{
  int rval = 0;
  int t;
  double *x = NULL, *rc = NULL;
  char *isint = NULL;
  DGG_constraint_t *base = NULL;
  bool   not_nicefied = true;
  int new_pos = cut_list->n;

  //  DGG_constraint_t *keep_origbase = DGG_copyConstraint(orig_base); //for debug only ------

  if (orig_base->sense == 'L') return 0;
  if (orig_base->nz == 0) return 0;

#define CGL_TWOMIR_LARGE_RHS 1.0e10
#ifdef CGL_TWOMIR_LARGE_RHS 
  if (fabs(orig_base->rhs) > CGL_TWOMIR_LARGE_RHS )
    return 0;
#endif
  rval = DGG_transformConstraint(data, &x, &rc, &isint, orig_base);
  double frac = frac_part(orig_base->rhs);
  //printf ("frac = %.7f, r %.7f, fr %.7f\n", frac, orig_base->rhs, floor(orig_base->rhs));
  if (rval || frac < data->gomory_threshold || frac > 1-data->gomory_threshold){
#ifndef TWOMIR_LESS_MALLOC
    free (x); free (rc); free (isint);
#endif
    return 0;
  }

  int min_t = t_min;
  int min_q = q_min;
  if (orig_base->sense == 'G' && min_t < 1) min_t = 1;
  if (orig_base->sense == 'G' && min_q < 1) min_q = 1;
  
  if (min_q > 0 &&  min_t > 0 ) {
    not_nicefied = false;
    rval = DGG_nicefyConstraint(solver_ptr, data, orig_base);
    DGG_CHECKRVAL(rval, rval);

    if (orig_base->nz == 0){
      if(talk)
	printf ("2mir_test: Nicefy returns empty constraint\n");
      rval = 0; goto CLEANUP;
    }
  }

  for(t = min_t; t <= t_max ; t++){
    if (t == 0) continue;

    base = DGG_copyConstraint(orig_base);
    DGG_TEST(!base, 1, "error making copy of base");

    DGG_scaleConstraint (base, t);

    if(not_nicefied){
      rval = DGG_nicefyConstraint(solver_ptr, data, base);
      DGG_CHECKRVAL(rval, rval);
      if (base->nz == 0){
	 if(talk)
	   printf ("2mir_test: Nicefy returns empty constraint\n");
	 goto MIR_DONE;
      }
    }
    
    if ( DGG_isBaseTrivial(data, base) ) goto MIR_DONE;

    rval = DGG_addMirToList(base, isint, x, cut_list, data, orig_base);
    DGG_CHECKRVAL(rval, rval);

  MIR_DONE:
    DGG_freeConstraint(base);
  }

  for( t = min_q; t <= q_max; t++ ){
    if (t == 0) continue;

    base = DGG_copyConstraint(orig_base);
    DGG_TEST(!base, 1, "error making copy of base");

    DGG_scaleConstraint (base, t);

    if(not_nicefied){
      rval = DGG_nicefyConstraint(solver_ptr, data, base);
      DGG_CHECKRVAL(rval, rval);
      if (base->nz == 0){
	 if(talk)
	   printf ("2mir_test: Nicefy returns empty constraint\n");
	 goto TWOMIR_DONE;
      }
    }
    
    if ( DGG_isBaseTrivial(data, base) ) goto TWOMIR_DONE;

    rval = DGG_add2stepToList(base, isint, x, rc, cut_list, data, orig_base);
    DGG_CHECKRVAL(rval, rval);
    
  TWOMIR_DONE:
    DGG_freeConstraint(base);
  }

  int i;
  for ( i = cut_list->n-1; i>=new_pos; i--){
    DGG_constraint_t *lcut = cut_list->c[i];

    rval = DGG_unTransformConstraint(data, lcut);
    DGG_CHECKRVAL(rval, rval);
    
    rval = DGG_substituteSlacks(solver_ptr, data, lcut);
    DGG_CHECKRVAL(rval, rval);

 

    if ( !DGG_isCutDesirable(lcut, data) ){
      DGG_list_delcut (cut_list, i);
      continue;
    }
    //else  testus(lcut);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    /*
    if ( data->opt_x && DGG_cutsOffPoint(data->opt_x, lcut) ){
      //DGG_cutDisplay_sparse(data, data->opt_x, lcut, stdout);
      DGG_TEST(1,1, "new cut is infeasible for optimal solution\n");
    }
    */
  }

 CLEANUP:
#ifndef TWOMIR_LESS_MALLOC
  if (x) free(x);
  if (rc) free (rc);
  if (isint) free(isint);
#endif
  return 0;
}

int
DGG_addMirToList ( DGG_constraint_t *base, char *isint, double * /*x*/,
		   DGG_list_t *list, DGG_data_t * /*data*/,
		   DGG_constraint_t * /*orig_base*/ )
{
  int rval = 0;
  DGG_constraint_t *cut = NULL;

  rval = DGG_buildMir(isint, base, &cut); 
  DGG_CHECKRVAL(rval, rval);

  DGG_list_addcut(list, cut, DGG_TMIR_CUT, 0.0);

  return 0;
}

int
DGG_add2stepToList ( DGG_constraint_t *base, char *isint, double * /*x*/,
				double *rc, DGG_list_t *list, DGG_data_t *data,
		     DGG_constraint_t * /*orig_base*/ )
{
  int rval;
  DGG_constraint_t *cut = NULL;
  int i;
  double norm_val, best_norm_val, best_norm_alpha=-1.0;
  double rc_val, best_rc_val,  best_rc_alpha=-1.0;
  double vht, bht, alpha;
  
  best_rc_val = best_norm_val = COIN_DBL_MAX;
  
  bht = ABOV(base->rhs);

  double best_rc = 0;
  for(i=0; i<base->nz; i++) if (isint[i]) best_rc = std::max(best_rc, fabs(rc[i]));
  double  rc_cutoff = best_rc / 10;

  for(i=0; i<base->nz; i++){
    if (!isint[i]) continue;
    if (fabs(rc[i]) <= rc_cutoff) continue; //too unimportant

    vht = ABOV(base->coeff[i]);
    if(vht >= bht)  continue;  // too big
    if(vht < bht/a_max) continue; // too small
    alpha = vht;
    int kk = 1;
    while ( !DGG_is2stepValid(alpha, bht) &&  bht/alpha <= a_max) {
      alpha = vht/kk; 
      kk++;
      if (kk>1000)
        break;
    }
    if ( !DGG_is2stepValid(alpha, bht) )    continue;
      
    rval = DGG_build2step(alpha, isint, base, &cut);
    DGG_CHECKRVAL(rval, rval);

    rc_val = COIN_DBL_MAX; // this gives a lower bound on obj. fn. improvement

    for(i=0; i<cut->nz; i++) if(cut->coeff[i]> 1E-6){
      rc_val = std::min(rc_val, fabs(rc[i])/cut->coeff[i]);
    }
    rc_val *= cut->rhs;

    norm_val = 0; // this is the square of the L2 norm
 
    for(i=0; i<cut->nz; i++) if(cut->coeff[i]> 1E-6){
      norm_val += (cut->coeff[i]*cut->coeff[i]);
    }

    norm_val /= (cut->rhs * cut->rhs + 1.0);
         
    if (rc_val < best_rc_val )  {	
      best_rc_val = rc_val; best_rc_alpha = alpha;  }

    if (norm_val < best_norm_val ) {	
      best_norm_val = norm_val;  best_norm_alpha = alpha;  }

    DGG_freeConstraint(cut);
  }
 
  if( best_rc_val> 1E-6 && best_rc_alpha != -1.0){
    rval = DGG_build2step(best_rc_alpha, isint, base, &cut);
    DGG_CHECKRVAL(rval, rval);
    DGG_list_addcut(list, cut, DGG_2STEP_CUT, best_rc_alpha);
  }
  else if (best_norm_alpha != -1.0){
    rval = DGG_build2step(best_norm_alpha, isint, base, &cut);
    DGG_CHECKRVAL(rval, rval);
    DGG_list_addcut(list, cut, DGG_2STEP_CUT, best_norm_alpha);
  }

  return 0;
}

int DGG_buildMir( char *isint,
		  DGG_constraint_t *base,
		  DGG_constraint_t **cut_out )
{
  int i, lnz = 0;
  double b   = (base->rhs);
  double bht = ABOV(b);
  double bup = ceil(b);
  DGG_constraint_t *tmir = NULL;

  DGG_TEST( base->sense == 'L', 1, "this form not valid for L");
  DGG_TEST( base->nz == 0, 1, "base must have some coefficients\n");

  tmir = DGG_newConstraint( base->nz );
  tmir->sense = 'G';
  tmir->rhs = bht * bup;

  for(i=0; i<base->nz; i++){
    double v   = base->coeff[i];

    if (!isint[i]) {
      if (v > 0.0) tmir->coeff[lnz] = v;
      else         tmir->coeff[lnz] = 0.0;
    }
    else {
      double vht = ABOV(v); 
      DGG_IF_EXIT( vht<0, 1, "negative vht");
      tmir->coeff[lnz]  = bht * floor(v) + DGG_MIN(bht,vht);
    }

    tmir->index[lnz] = base->index[i];
    lnz += 1;
  }

  tmir->nz = lnz;
  *cut_out = tmir;

  return 0;
}

int DGG_build2step( double alpha,
		    char *isint,
		    DGG_constraint_t *base,
		    DGG_constraint_t **cut_out )

{
  DGG_constraint_t *tmir = 0;

  int i,  lnz = 0;
  double vht, bht, bup, rho, tau, k;
  double b = (base->rhs);

  DGG_TEST( base->sense == 'L', 1, "this form not valid for L");
  DGG_TEST( base->nz == 0, 1, "base must have some coefficients\n");

  bht = ABOV(b);  
  bup = ceil(b); 
  tau = ceil(bht/alpha); 
  rho = bht - alpha*floor(bht/alpha);

  /* ensure bht > alpha > 0 */
  DGG_TEST3( (bht <= alpha) || (alpha <= 0.0), 1, "bad alpha (%f) / bht (%f) pair", alpha, bht);
  /* ensure that we are not in a limiting case */
  DGG_TEST( DGG_is_a_multiple_of_b(alpha, bht), 1, "can't generate simple 2mir in limiting case");
  /* ensure that rho is not zero */
  DGG_TEST2( rho < DGG_MIN_RHO, 1, "rho (%f) too small", rho);

  /* initialize constraint */
  tmir = DGG_newConstraint( base->nz );

  tmir->rhs = bup*tau*rho;
  tmir->sense = 'G';

  /* compute cut coefficients */
  for(i=0; i<base->nz; i++){
    double v   = base->coeff[i];

    if (!isint[i]) {
      if (v > 0.0) tmir->coeff[lnz] = v;
      else         tmir->coeff[lnz] = 0.0;
    }
    else {
      vht = v - floor(v);
      DGG_IF_EXIT( vht < 0.0, 1, "negative vht");
      k   = DGG_MIN(tau-1,floor(vht/alpha));
      tmir->coeff[lnz]  = floor(v)*tau*rho +  k*rho + DGG_MIN(rho,vht-k*alpha);
    }

    tmir->index[lnz] = base->index[i];
    lnz += 1;
  }

  tmir->nz = lnz;
  *cut_out = tmir;

  return 0;
}

/******************* TEST / DEBUGGING ROUTINES ********************************/

/* DGG_is2stepValid:
   checks that:

   bht > alpha > 0
   (1/alpha) >= tau > (bht/alpha)
*/
int DGG_is2stepValid(double alpha, double bht)
{

  /* d */
  double tau;

  /* ensure that alpha is not null or negative */
  if ( alpha < DGG_MIN_ALPHA )
    return 0;

  /* compute tau and tau_lim */
  tau = ceil( bht / alpha );

  /* make sure alpha is not a divisor of bht */
  if ( DGG_is_a_multiple_of_b(alpha, bht) )
    return 0;

  /* page 15, definition 12 */
  /* check if alpha is admissible for simple-2-step-tmir */

  if ( (bht > alpha) && (alpha > 0.0) )
    if ( (1/alpha) >= tau )
      return 1;

  /* not admissible */
  return 0;
}

/* checks that its worth doing a 1MIR on the constraint. More precisely,

- Is the RHS null? 
- Are there any integer variables set at fractional values?           */

int DGG_isBaseTrivial(DGG_data_t *d, DGG_constraint_t* c)
{

  /* is rhs sufficiently fractional */
  if ( frac_part(ABOV(c->rhs)) < d->gomory_threshold )
    return 1;

  if ( (1.0 - frac_part(ABOV(c->rhs))) < d->gomory_threshold )
    return 1;

  return 0;
}


/* tests lhs vs rhs of a constraint */
int DGG_isConstraintViolated(DGG_data_t *d, DGG_constraint_t *c)
{
  double lhs = DGG_cutLHS(c, d->x);
  double rhs = c->rhs;

  /* compare LHS and RHS */
  if (c->sense == 'G')
    if ( lhs > (rhs - DGG_NULL_SLACK) )
      return 0;
  if (c->sense == 'L')
    if ( lhs < (rhs + DGG_NULL_SLACK) )
      return 0;
  if (c->sense == 'E')
    if ( fabs(lhs - rhs) < DGG_NULL_SLACK )
      return 0;

  return 0;

}

double DGG_cutLHS(DGG_constraint_t *c, double *x)
{

  int i;
  double lhs = 0.0;

  for(i=0; i < c->nz; i++)
    lhs += c->coeff[i]*x[c->index[i]];

  return lhs;
}

int DGG_isCutDesirable(DGG_constraint_t *c, DGG_data_t *d)
{
  double lhs, rhs;

  lhs = DGG_cutLHS(c, d->x);
  rhs = c->rhs;

  if (c->nz > 500) return 0;

  /* if the cut is not violated, return 0 */
  if (c->sense == 'G')
    if ( lhs > (rhs - DGG_NULL_SLACK) )
      return 0;
  if (c->sense == 'L')
    if ( lhs < (rhs + DGG_NULL_SLACK) )
      return 0;
  if (c->sense == 'E')
    if ( fabs(lhs - rhs) < DGG_NULL_SLACK )
      return 0;
  return 1;
}

/******************** SIMPLE MACROS AND FUNCTIONS *****************************/

int DGG_is_even(double vht, double bht, int tau, int q)
{

  double v2 = V2I(bht, tau, q);

  if ( vht > v2 )
    return 1;

  return 0;
}

int DGG_is_a_multiple_of_b(double a, double b)
{
  double c = b/a;
  
  if ( (b - a*floor(c)) < DGG_MIN_RHO )
    return 1;

  return 0;
}

int DGG_cutsOffPoint(double *x, DGG_constraint_t *cut)
{

  int i;
  double LHS = 0.0;

	for(i=0; i < cut->nz; i++)
	  LHS += cut->coeff[i]*(x[ cut->index[i] ]);

  //fprintf(stdout, "LHS = %f, SENSE = %c, RHS = %f\n", LHS, cut->sense, cut->rhs);
  if ( cut->sense == 'E' )
    if ( fabs(LHS - cut->rhs) > DGG_NULL_SLACK )
      goto BAD;
  if (cut->sense == 'G' )
    if ( (cut->rhs - LHS) > DGG_NULL_SLACK )
      goto BAD;
  if (cut->sense == 'L' )
    if ( (LHS - cut->rhs) > DGG_NULL_SLACK )
      goto BAD;

  return 0;

  BAD:

  fprintf(stdout, "LHS = %f, SENSE = %c, RHS = %f\n", LHS, cut->sense, cut->rhs);
  DGG_TEST(1, 1, "found a bad cut!");
  return 0;
}
// Returns true if needs optimal basis to do cuts
bool 
CglTwomir::needsOptimalBasis() const
{
  return true;
}

// Away stuff
void CglTwomir::setAway(double value)
{
  if (value>0.0&&value<=0.5)
    away_=value;
}
double CglTwomir::getAway() const
{
  return away_;
}

// Away stuff at root
void CglTwomir::setAwayAtRoot(double value)
{
  if (value>0.0&&value<=0.5)
    awayAtRoot_=value;
}
double CglTwomir::getAwayAtRoot() const
{
  return awayAtRoot_;
}

// This can be used to refresh any information
void 
CglTwomir::refreshSolver(OsiSolverInterface * solver)
{
  if (originalSolver_) {
    delete originalSolver_;
    originalSolver_ = solver->clone();
  }
  // Get integer information
  solver->getColType(true);
}
// Create C++ lines to get to current state
std::string
CglTwomir::generateCpp( FILE * fp) 
{
  CglTwomir other;
  fprintf(fp,"0#include \"CglTwomir.hpp\"\n");
  fprintf(fp,"3  CglTwomir twomir;\n");
  if (t_min_!=other.t_min_||t_max_!=other.t_max_)
    fprintf(fp,"3  twomir.setMirScale(%d,%d);\n",t_min_,t_max_);
  else
    fprintf(fp,"4  twomir.setMirScale(%d,%d);\n",t_min_,t_max_);
  if (q_min_!=other.q_min_||q_max_!=other.q_max_)
    fprintf(fp,"3  twomir.setTwomirScale(%d,%d);\n",q_min_,q_max_);
  else
    fprintf(fp,"4  twomir.setTwomirScale(%d,%d);\n",q_min_,q_max_);
  if (do_mir_!=other.do_mir_||do_2mir_!=other.do_2mir_||
      do_tab_!=other.do_tab_||do_form_!=other.do_form_)
    fprintf(fp,"3  twomir.setCutTypes(%s,%s,%s,%s);\n",
	    do_mir_ ? "true" : "false",
	    do_2mir_ ? "true" : "false",
	    do_tab_ ? "true" : "false",
	    do_form_ ? "true" : "false");
  else
    fprintf(fp,"4  twomir.setCutTypes(%s,%s,%s,%s);\n",
	    do_mir_ ? "true" : "false",
	    do_2mir_ ? "true" : "false",
	    do_tab_ ? "true" : "false",
	    do_form_ ? "true" : "false");
  if (a_max_!=other.a_max_)
    fprintf(fp,"3  twomir.setAMax(%d);\n",a_max_);
  else
    fprintf(fp,"4  twomir.setAMax(%d);\n",a_max_);
  if (max_elements_!=other.max_elements_)
    fprintf(fp,"3  twomir.setMaxElements(%d);\n",max_elements_);
  else
    fprintf(fp,"4  twomir.setMaxElements(%d);\n",max_elements_);
  if (max_elements_root_!=other.max_elements_root_)
    fprintf(fp,"3  twomir.setMaxElementsRoot(%d);\n",max_elements_root_);
  else
    fprintf(fp,"4  twomir.setMaxElementsRoot(%d);\n",max_elements_root_);
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  twomir.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  twomir.setAggressiveness(%d);\n",getAggressiveness());
  return "twomir";
}
