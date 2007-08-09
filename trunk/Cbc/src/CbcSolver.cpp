// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
   
#include "CbcConfig.h"
#include "CoinPragma.hpp"

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <iostream>


#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
// Version
#define CBCVERSION "1.04.00"

#include "CoinMpsIO.hpp"
#include "CoinModel.hpp"

#include "ClpFactorization.hpp"
#include "ClpQuadraticObjective.hpp"
#include "CoinTime.hpp"
#include "ClpSimplex.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpSolve.hpp"
#include "ClpMessage.hpp"
#include "ClpPackedMatrix.hpp"
#include "ClpPlusMinusOneMatrix.hpp"
#include "ClpNetworkMatrix.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpLinearObjective.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "ClpPrimalColumnDantzig.hpp"
#include "ClpPresolve.hpp"
#include "CbcOrClpParam.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiChooseVariable.hpp"
//#define USER_HAS_FAKE_CLP
//#define USER_HAS_FAKE_CBC
//#define CLP_MALLOC_STATISTICS
#ifdef CLP_MALLOC_STATISTICS
#include <malloc.h>
#include <exception>
#include <new>
static double malloc_times=0.0;
static double malloc_total=0.0;
static int malloc_amount[]={0,32,128,256,1024,4096,16384,65536,262144,INT_MAX};
static int malloc_n=10;
double malloc_counts[10]={0,0,0,0,0,0,0,0,0,0};
void * operator new (size_t size) throw (std::bad_alloc)
{
  malloc_times ++;
  malloc_total += size;
  int i;
  for (i=0;i<malloc_n;i++) {
    if ((int) size<=malloc_amount[i]) {
      malloc_counts[i]++;
      break;
    }
  }
  void * p =malloc(size);
  return p;
}
void operator delete (void *p) throw()
{
  free(p);
}
static void malloc_stats2()
{
  double average = malloc_total/malloc_times;
  printf("count %g bytes %g - average %g\n",malloc_times,malloc_total,average);
  for (int i=0;i<malloc_n;i++) 
    printf("%g ",malloc_counts[i]);
  printf("\n");
  malloc_times=0.0;
  malloc_total=0.0;
  memset(malloc_counts,0,sizeof(malloc_counts));
}
#endif
#ifdef DMALLOC
#include "dmalloc.h"
#endif
#ifdef WSSMP_BARRIER
#define FOREIGN_BARRIER
#endif
#ifdef UFL_BARRIER
#define FOREIGN_BARRIER
#endif
#ifdef TAUCS_BARRIER
#define FOREIGN_BARRIER
#endif
#include "CoinWarmStartBasis.hpp"

#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#ifndef COIN_HAS_LINK
//#ifdef COIN_HAS_ASL
#define COIN_HAS_LINK
//#endif
#endif
#ifdef COIN_HAS_LINK
#include "CbcLinked.hpp"
#endif
#include "CglPreProcess.hpp"
#include "CglCutGenerator.hpp"
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglRedSplit.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglTwomir.hpp"
#include "CglDuplicateRow.hpp"
#include "CglStored.hpp"
#include "CglLandP.hpp"
#include "CglResidualCapacity.hpp"

#include "CbcModel.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicRINS.hpp"
#include "CbcTreeLocal.hpp"
#include "CbcCompareActual.hpp"
#include "CbcBranchActual.hpp"
#include  "CbcOrClpParam.hpp"
#include  "CbcCutGenerator.hpp"
#include  "CbcStrategy.hpp"
#include "CbcBranchLotsize.hpp"

#include "OsiClpSolverInterface.hpp"
#ifdef COIN_HAS_ASL
#include "Cbc_ampl.h"
static bool usingAmpl=false;
#endif
static double totalTime=0.0;
static void statistics(ClpSimplex * originalModel, ClpSimplex * model);
static bool maskMatches(const int * starts, char ** masks,
			std::string & check);
static void generateCode(CbcModel * model, const char * fileName,int type,int preProcess);
//#ifdef NDEBUG
//#undef NDEBUG
//#endif
// define (probably dummy fake main programs for UserClp and UserCbc
void fakeMain (ClpSimplex & model,OsiSolverInterface & osiSolver, CbcModel & babSolver);
void fakeMain2 (ClpSimplex & model,OsiClpSolverInterface & osiSolver,int options);

// Allow for interrupts
// But is this threadsafe ? (so switched off by option)

#include "CoinSignal.hpp"
static CbcModel * currentBranchModel = NULL;

extern "C" {
   static void signal_handler(int whichSignal)
   {
     if (currentBranchModel!=NULL) {
       currentBranchModel->setMaximumNodes(0); // stop at next node
       currentBranchModel->setMaximumSeconds(0.0); // stop 
     }
     return;
   }
}

int CbcOrClpRead_mode=1;
FILE * CbcOrClpReadCommand=stdin;
static bool noPrinting=false;
static int * analyze(OsiClpSolverInterface * solverMod, int & numberChanged, double & increment,
                     bool changeInt,  CoinMessageHandler * generalMessageHandler)
{
  OsiSolverInterface * solver = solverMod->clone();
  char generalPrint[200];
  if (0) {
    // just get increment
    CbcModel model(*solver);
    model.analyzeObjective();
    double increment2=model.getCutoffIncrement();
    printf("initial cutoff increment %g\n",increment2);
  }
  const double *objective = solver->getObjCoefficients() ;
  const double *lower = solver->getColLower() ;
  const double *upper = solver->getColUpper() ;
  int numberColumns = solver->getNumCols() ;
  int numberRows = solver->getNumRows();
  double direction = solver->getObjSense();
  int iRow,iColumn;

  // Row copy
  CoinPackedMatrix matrixByRow(*solver->getMatrixByRow());
  const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();

  // Column copy
  CoinPackedMatrix  matrixByCol(*solver->getMatrixByCol());
  const double * element = matrixByCol.getElements();
  const int * row = matrixByCol.getIndices();
  const CoinBigIndex * columnStart = matrixByCol.getVectorStarts();
  const int * columnLength = matrixByCol.getVectorLengths();

  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();

  char * ignore = new char [numberRows];
  int * changed = new int[numberColumns];
  int * which = new int[numberRows];
  double * changeRhs = new double[numberRows];
  memset(changeRhs,0,numberRows*sizeof(double));
  memset(ignore,0,numberRows);
  numberChanged=0;
  int numberInteger=0;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (upper[iColumn] > lower[iColumn]+1.0e-8&&solver->isInteger(iColumn)) 
      numberInteger++;
  }
  bool finished=false;
  while (!finished) {
    int saveNumberChanged = numberChanged;
    for (iRow=0;iRow<numberRows;iRow++) {
      int numberContinuous=0;
      double value1=0.0,value2=0.0;
      bool allIntegerCoeff=true;
      double sumFixed=0.0;
      int jColumn1=-1,jColumn2=-1;
      for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
        int jColumn = column[j];
        double value = elementByRow[j];
        if (upper[jColumn] > lower[jColumn]+1.0e-8) {
          if (!solver->isInteger(jColumn)) {
            if (numberContinuous==0) {
              jColumn1=jColumn;
              value1=value;
            } else {
              jColumn2=jColumn;
              value2=value;
            }
            numberContinuous++;
          } else {
            if (fabs(value-floor(value+0.5))>1.0e-12)
              allIntegerCoeff=false;
          }
        } else {
          sumFixed += lower[jColumn]*value;
        }
      }
      double low = rowLower[iRow];
      if (low>-1.0e20) {
        low -= sumFixed;
        if (fabs(low-floor(low+0.5))>1.0e-12)
          allIntegerCoeff=false;
      }
      double up = rowUpper[iRow];
      if (up<1.0e20) {
        up -= sumFixed;
        if (fabs(up-floor(up+0.5))>1.0e-12)
          allIntegerCoeff=false;
      }
      if (!allIntegerCoeff)
        continue; // can't do
      if (numberContinuous==1) {
        // see if really integer
        // This does not allow for complicated cases
        if (low==up) {
          if (fabs(value1)>1.0e-3) {
            value1 = 1.0/value1;
            if (fabs(value1-floor(value1+0.5))<1.0e-12) {
              // integer
              changed[numberChanged++]=jColumn1;
              solver->setInteger(jColumn1);
              if (upper[jColumn1]>1.0e20)
                solver->setColUpper(jColumn1,1.0e20);
              if (lower[jColumn1]<-1.0e20)
                solver->setColLower(jColumn1,-1.0e20);
            }
          }
        } else {
          if (fabs(value1)>1.0e-3) {
            value1 = 1.0/value1;
            if (fabs(value1-floor(value1+0.5))<1.0e-12) {
              // This constraint will not stop it being integer
              ignore[iRow]=1;
            }
          }
        }
      } else if (numberContinuous==2) {
        if (low==up) {
          /* need general theory - for now just look at 2 cases -
             1 - +- 1 one in column and just costs i.e. matching objective
             2 - +- 1 two in column but feeds into G/L row which will try and minimize
          */
          if (fabs(value1)==1.0&&value1*value2==-1.0&&!lower[jColumn1]
              &&!lower[jColumn2]) {
            int n=0;
            int i;
            double objChange=direction*(objective[jColumn1]+objective[jColumn2]);
            double bound = CoinMin(upper[jColumn1],upper[jColumn2]);
            bound = CoinMin(bound,1.0e20);
            for ( i=columnStart[jColumn1];i<columnStart[jColumn1]+columnLength[jColumn1];i++) {
              int jRow = row[i];
              double value = element[i];
              if (jRow!=iRow) {
                which[n++]=jRow;
                changeRhs[jRow]=value;
              }
            }
            for ( i=columnStart[jColumn1];i<columnStart[jColumn1]+columnLength[jColumn1];i++) {
              int jRow = row[i];
              double value = element[i];
              if (jRow!=iRow) {
                if (!changeRhs[jRow]) {
                  which[n++]=jRow;
                  changeRhs[jRow]=value;
                } else {
                  changeRhs[jRow]+=value;
                }
              }
            }
            if (objChange>=0.0) {
              // see if all rows OK
              bool good=true;
              for (i=0;i<n;i++) {
                int jRow = which[i];
                double value = changeRhs[jRow];
                if (value) {
                  value *= bound;
                  if (rowLength[jRow]==1) {
                    if (value>0.0) {
                      double rhs = rowLower[jRow];
                      if (rhs>0.0) {
                        double ratio =rhs/value;
                        if (fabs(ratio-floor(ratio+0.5))>1.0e-12)
                          good=false;
                      }
                    } else {
                      double rhs = rowUpper[jRow];
                      if (rhs<0.0) {
                        double ratio =rhs/value;
                        if (fabs(ratio-floor(ratio+0.5))>1.0e-12)
                          good=false;
                      }
                    }
                  } else if (rowLength[jRow]==2) {
                    if (value>0.0) {
                      if (rowLower[jRow]>-1.0e20)
                        good=false;
                    } else {
                      if (rowUpper[jRow]<1.0e20)
                        good=false;
                    }
                  } else {
                    good=false;
                  }
                }
              }
              if (good) {
                // both can be integer
                changed[numberChanged++]=jColumn1;
                solver->setInteger(jColumn1);
                if (upper[jColumn1]>1.0e20)
                  solver->setColUpper(jColumn1,1.0e20);
                if (lower[jColumn1]<-1.0e20)
                  solver->setColLower(jColumn1,-1.0e20);
                changed[numberChanged++]=jColumn2;
                solver->setInteger(jColumn2);
                if (upper[jColumn2]>1.0e20)
                  solver->setColUpper(jColumn2,1.0e20);
                if (lower[jColumn2]<-1.0e20)
                  solver->setColLower(jColumn2,-1.0e20);
              }
            }
            // clear
            for (i=0;i<n;i++) {
              changeRhs[which[i]]=0.0;
            }
          }
        }
      }
    }
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (upper[iColumn] > lower[iColumn]+1.0e-8&&!solver->isInteger(iColumn)) {
        double value;
        value = upper[iColumn];
        if (value<1.0e20&&fabs(value-floor(value+0.5))>1.0e-12) 
          continue;
        value = lower[iColumn];
        if (value>-1.0e20&&fabs(value-floor(value+0.5))>1.0e-12) 
          continue;
        bool integer=true;
        for (CoinBigIndex j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
          int iRow = row[j];
          if (!ignore[iRow]) {
            integer=false;
            break;
          }
        }
        if (integer) {
          // integer
          changed[numberChanged++]=iColumn;
          solver->setInteger(iColumn);
          if (upper[iColumn]>1.0e20)
            solver->setColUpper(iColumn,1.0e20);
          if (lower[iColumn]<-1.0e20)
            solver->setColLower(iColumn,-1.0e20);
        }
      }
    }
    finished = numberChanged==saveNumberChanged;
  }
  delete [] which;
  delete [] changeRhs;
  delete [] ignore;
  //if (numberInteger&&!noPrinting)
  //printf("%d integer variables",numberInteger);
  if (changeInt) {
    //if (!noPrinting) {
    //if (numberChanged)
    //  printf(" and %d variables made integer\n",numberChanged);
    //else
    //  printf("\n");
    //}
    delete [] ignore;
    //increment=0.0;
    if (!numberChanged) {
      delete [] changed;
      delete solver;
      return NULL;
    } else {
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
        if (solver->isInteger(iColumn))
          solverMod->setInteger(iColumn);
      }
      delete solver;
      return changed;
    }
  } else {
    //if (!noPrinting) {
    //if (numberChanged)
    //  printf(" and %d variables could be made integer\n",numberChanged);
    //else
    //  printf("\n");
    //}
    // just get increment
    int logLevel=generalMessageHandler->logLevel();
    CbcModel model(*solver);
    model.passInMessageHandler(generalMessageHandler);
    if (noPrinting)
      model.setLogLevel(0);
    model.analyzeObjective();
    generalMessageHandler->setLogLevel(logLevel);
    double increment2=model.getCutoffIncrement();
    if (increment2>increment&&increment2>0.0) {
      if (!noPrinting) {
	sprintf(generalPrint,"Cutoff increment increased from %g to %g",increment,increment2);
	CoinMessages generalMessages = solverMod->getModelPtr()->messages();
	generalMessageHandler->message(CLP_GENERAL,generalMessages)
	  << generalPrint
	  <<CoinMessageEol;
      }
      increment=increment2;
    }
    delete solver;
    numberChanged=0;
    delete [] changed;
    return NULL;
  }
}
#ifdef COIN_HAS_ASL
/*  Returns OsiSolverInterface (User should delete)
    On entry numberKnapsack is maximum number of Total entries
*/
static OsiSolverInterface *  
expandKnapsack(CoinModel & model, int * whichColumn, int * knapsackStart, 
	       int * knapsackRow, int &numberKnapsack,
	       CglStored & stored, int logLevel,
	       int fixedPriority, int SOSPriority,CoinModel & tightenedModel)
{
  int maxTotal = numberKnapsack;
  // load from coin model
  OsiSolverLink *si = new OsiSolverLink();
  OsiSolverInterface * finalModel=NULL;
  si->setDefaultMeshSize(0.001);
  // need some relative granularity
  si->setDefaultBound(100.0);
  si->setDefaultMeshSize(0.01);
  si->setDefaultBound(100000.0);
  si->setIntegerPriority(1000);
  si->setBiLinearPriority(10000);
  si->load(model,true,logLevel);
  // get priorities
  const int * priorities=model.priorities();
  int numberColumns = model.numberColumns();
  if (priorities) {
    OsiObject ** objects = si->objects();
    int numberObjects = si->numberObjects();
    for (int iObj = 0;iObj<numberObjects;iObj++) {
      int iColumn = objects[iObj]->columnNumber();
      if (iColumn>=0&&iColumn<numberColumns) {
#ifndef NDEBUG
	OsiSimpleInteger * obj =
	  dynamic_cast <OsiSimpleInteger *>(objects[iObj]) ;
#endif
	assert (obj);
	int iPriority = priorities[iColumn];
	if (iPriority>0)
	  objects[iObj]->setPriority(iPriority);
      }
    }
    if (fixedPriority>0) {
      si->setFixedPriority(fixedPriority);
    }
    if (SOSPriority<0)
      SOSPriority=100000;
  } 
  CoinModel coinModel=*si->coinModel();
  assert(coinModel.numberRows()>0);
  tightenedModel = coinModel;
  int numberRows = coinModel.numberRows();
  // Mark variables
  int * whichKnapsack = new int [numberColumns];
  int iRow,iColumn;
  for (iColumn=0;iColumn<numberColumns;iColumn++) 
    whichKnapsack[iColumn]=-1;
  int kRow;
  bool badModel=false;
  // analyze
  if (logLevel>1) {
    for (iRow=0;iRow<numberRows;iRow++) {
      /* Just obvious one at first
	 positive non unit coefficients 
	 all integer
	 positive rowUpper
	 for now - linear (but further down in code may use nonlinear)
	 column bounds should be tight
      */
      //double lower = coinModel.getRowLower(iRow);
      double upper = coinModel.getRowUpper(iRow);
      if (upper<1.0e10) {
	CoinModelLink triple=coinModel.firstInRow(iRow);
	bool possible=true;
	int n=0;
	int n1=0;
	while (triple.column()>=0) {
	  int iColumn = triple.column();
	  const char *  el = coinModel.getElementAsString(iRow,iColumn);
	  if (!strcmp("Numeric",el)) {
	    if (coinModel.columnLower(iColumn)==coinModel.columnUpper(iColumn)) {
	      triple=coinModel.next(triple);
	      continue; // fixed
	    }
	    double value=coinModel.getElement(iRow,iColumn);
	    if (value<0.0) {
	      possible=false;
	    } else {
	      n++;
	      if (value==1.0)
		n1++;
	      if (coinModel.columnLower(iColumn)<0.0)
		possible=false;
	      if (!coinModel.isInteger(iColumn))
		possible=false;
	      if (whichKnapsack[iColumn]>=0)
		possible=false;
	    }
	  } else {
	    possible=false; // non linear
	  } 
	  triple=coinModel.next(triple);
	}
	if (n-n1>1&&possible) {
	  double lower = coinModel.getRowLower(iRow);
	  double upper = coinModel.getRowUpper(iRow);
	  CoinModelLink triple=coinModel.firstInRow(iRow);
	  while (triple.column()>=0) {
	    int iColumn = triple.column();
	    lower -= coinModel.columnLower(iColumn)*triple.value();
	    upper -= coinModel.columnLower(iColumn)*triple.value();
	    triple=coinModel.next(triple);
	  }
	  printf("%d is possible %g <=",iRow,lower);
	  // print
	  triple=coinModel.firstInRow(iRow);
	  while (triple.column()>=0) {
	    int iColumn = triple.column();
	    if (coinModel.columnLower(iColumn)!=coinModel.columnUpper(iColumn))
	      printf(" (%d,el %g up %g)",iColumn,triple.value(),
		     coinModel.columnUpper(iColumn)-coinModel.columnLower(iColumn));
	    triple=coinModel.next(triple);
	  }
	  printf(" <= %g\n",upper);
	}
      }
    }
  }
  numberKnapsack=0;
  for (kRow=0;kRow<numberRows;kRow++) {
    iRow=kRow;
    /* Just obvious one at first
       positive non unit coefficients 
       all integer
       positive rowUpper
       for now - linear (but further down in code may use nonlinear)
       column bounds should be tight
    */
    //double lower = coinModel.getRowLower(iRow);
    double upper = coinModel.getRowUpper(iRow);
    if (upper<1.0e10) {
      CoinModelLink triple=coinModel.firstInRow(iRow);
      bool possible=true;
      int n=0;
      int n1=0;
      while (triple.column()>=0) {
	int iColumn = triple.column();
	const char *  el = coinModel.getElementAsString(iRow,iColumn);
	if (!strcmp("Numeric",el)) {
	  if (coinModel.columnLower(iColumn)==coinModel.columnUpper(iColumn)) {
	    triple=coinModel.next(triple);
	    continue; // fixed
	  }
	  double value=coinModel.getElement(iRow,iColumn);
	  if (value<0.0) {
	    possible=false;
	  } else {
	    n++;
	    if (value==1.0)
	      n1++;
	    if (coinModel.columnLower(iColumn)<0.0)
	      possible=false;
	    if (!coinModel.isInteger(iColumn))
	      possible=false;
	    if (whichKnapsack[iColumn]>=0)
	      possible=false;
	  }
	} else {
	  possible=false; // non linear
	} 
	triple=coinModel.next(triple);
      }
      if (n-n1>1&&possible) {
	// try
	CoinModelLink triple=coinModel.firstInRow(iRow);
	while (triple.column()>=0) {
	  int iColumn = triple.column();
	  if (coinModel.columnLower(iColumn)!=coinModel.columnUpper(iColumn))
	    whichKnapsack[iColumn]=numberKnapsack;
	  triple=coinModel.next(triple);
	}
	knapsackRow[numberKnapsack++]=iRow;
      }
    }
  }
  if (logLevel>0)
    printf("%d out of %d candidate rows are possible\n",numberKnapsack,numberRows);
  // Check whether we can get rid of nonlinearities
  /* mark rows
     -2 in knapsack and other variables
     -1 not involved
     n only in knapsack n
  */
  int * markRow = new int [numberRows];
  for (iRow=0;iRow<numberRows;iRow++) 
    markRow[iRow]=-1;
  int canDo=1; // OK and linear
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinModelLink triple=coinModel.firstInColumn(iColumn);
    int iKnapsack = whichKnapsack[iColumn];
    bool linear=true;
    // See if quadratic objective
    const char * expr = coinModel.getColumnObjectiveAsString(iColumn);
    if (strcmp(expr,"Numeric")) {
      linear=false;
    }
    while (triple.row()>=0) {
      int iRow = triple.row();
      if (iKnapsack>=0) {
	if (markRow[iRow]==-1) {
	  markRow[iRow]=iKnapsack;
	} else if (markRow[iRow]!=iKnapsack) {
	  markRow[iRow]=-2;
	}
      }
      const char * expr = coinModel.getElementAsString(iRow,iColumn);
      if (strcmp(expr,"Numeric")) {
	linear=false;
      }
      triple=coinModel.next(triple);
    }
    if (!linear) {
      if (whichKnapsack[iColumn]<0) {
	canDo=0;
	break;
      } else {
	canDo=2;
      }
    }
  }
  int * markKnapsack = NULL;
  double * coefficient = NULL;
  double * linear = NULL;
  int * whichRow = NULL;
  int * lookupRow = NULL;
  badModel=(canDo==0);
  if (numberKnapsack&&canDo) {
    /* double check - OK if
       no nonlinear
       nonlinear only on columns in knapsack
       nonlinear only on columns in knapsack * ONE other - same for all in knapsack
       AND that is only row connected to knapsack
       (theoretically could split knapsack if two other and small numbers)
       also ONE could be ONE expression - not just a variable
    */
    int iKnapsack;
    markKnapsack = new int [numberKnapsack];
    coefficient = new double [numberKnapsack];
    linear = new double [numberColumns];
    for (iKnapsack=0;iKnapsack<numberKnapsack;iKnapsack++) 
      markKnapsack[iKnapsack]=-1;
    if (canDo==2) {
      for (iRow=-1;iRow<numberRows;iRow++) {
	int numberOdd;
	CoinPackedMatrix * row = coinModel.quadraticRow(iRow,linear,numberOdd);
	if (row) {
	  // see if valid
	  const double * element = row->getElements();
	  const int * column = row->getIndices();
	  const CoinBigIndex * columnStart = row->getVectorStarts();
	  const int * columnLength = row->getVectorLengths();
	  int numberLook = row->getNumCols();
	  for (int i=0;i<numberLook;i++) {
	    int iKnapsack=whichKnapsack[i];
	    if (iKnapsack<0) {
	      // might be able to swap - but for now can't have knapsack in
	      for (int j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
		int iColumn = column[j];
		if (whichKnapsack[iColumn]>=0) {
		  canDo=0; // no good
		  badModel=true;
		  break;
		}
	      }
	    } else {
	      // OK if in same knapsack - or maybe just one
	      int marked=markKnapsack[iKnapsack];
	      for (int j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
		int iColumn = column[j];
		if (whichKnapsack[iColumn]!=iKnapsack&&whichKnapsack[iColumn]>=0) {
		  canDo=0; // no good
		  badModel=true;
		  break;
		} else if (marked==-1) {
		  markKnapsack[iKnapsack]=iColumn;
		  marked=iColumn;
		  coefficient[iKnapsack]=element[j];
		  coinModel.associateElement(coinModel.columnName(iColumn),1.0);
		} else if (marked!=iColumn) {
		  badModel=true;
		  canDo=0; // no good
		  break;
		} else {
		  // could manage with different coefficients - but for now ...
		  assert(coefficient[iKnapsack]==element[j]);
		}
	      }
	    }
	  }
	  delete row;
	}
      }
    }
    if (canDo) {
      // for any rows which are cuts
      whichRow = new int [numberRows];
      lookupRow = new int [numberRows];
      bool someNonlinear=false;
      double maxCoefficient=1.0;
      for (iKnapsack=0;iKnapsack<numberKnapsack;iKnapsack++) {
	if (markKnapsack[iKnapsack]>=0) {
	  someNonlinear=true;
	  int iColumn = markKnapsack[iKnapsack];
	  maxCoefficient = CoinMax(maxCoefficient,fabs(coefficient[iKnapsack]*coinModel.columnUpper(iColumn)));
	}
      }
      if (someNonlinear) {
	// associate all columns to stop possible error messages
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  coinModel.associateElement(coinModel.columnName(iColumn),1.0);
	}
      }
      ClpSimplex tempModel;
      tempModel.loadProblem(coinModel);
      // Create final model - first without knapsacks
      int nCol=0;
      int nRow=0;
      for (iRow=0;iRow<numberRows;iRow++) {
	if (markRow[iRow]<0) {
	  lookupRow[iRow]=nRow;
	  whichRow[nRow++]=iRow;
	} else {
	  lookupRow[iRow]=-1;
	}
      }
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	if (whichKnapsack[iColumn]<0)
	  whichColumn[nCol++]=iColumn;
      }
      ClpSimplex finalModelX(&tempModel,nRow,whichRow,nCol,whichColumn,false,false,false);
      OsiClpSolverInterface finalModelY(&finalModelX,true);
      finalModel = finalModelY.clone();
      finalModelY.releaseClp();
      // Put back priorities
      const int * priorities=model.priorities();
      if (priorities) {
	finalModel->findIntegers(false);
	OsiObject ** objects = finalModel->objects();
	int numberObjects = finalModel->numberObjects();
	for (int iObj = 0;iObj<numberObjects;iObj++) {
	  int iColumn = objects[iObj]->columnNumber();
	  if (iColumn>=0&&iColumn<nCol) {
#ifndef NDEBUG 
	    OsiSimpleInteger * obj =
	      dynamic_cast <OsiSimpleInteger *>(objects[iObj]) ;
#endif
	    assert (obj);
	    int iPriority = priorities[whichColumn[iColumn]];
	    if (iPriority>0)
	      objects[iObj]->setPriority(iPriority);
	  }
	}
      }
      for (iRow=0;iRow<numberRows;iRow++) {
	whichRow[iRow]=iRow;
      }
      int numberOther=finalModel->getNumCols();
      int nLargest=0;
      int nelLargest=0;
      int nTotal=0;
      for (iKnapsack=0;iKnapsack<numberKnapsack;iKnapsack++) {
	iRow = knapsackRow[iKnapsack];
	int nCreate = maxTotal;
	int nelCreate=coinModel.expandKnapsack(iRow,nCreate,NULL,NULL,NULL,NULL);
	if (nelCreate<0)
	  badModel=true;
	nTotal+=nCreate;
	nLargest = CoinMax(nLargest,nCreate);
	nelLargest = CoinMax(nelLargest,nelCreate);
      }
      if (nTotal>maxTotal) 
	badModel=true;
      if (!badModel) {
	// Now arrays for building
	nelLargest = CoinMax(nelLargest,nLargest)+1;
	double * buildObj = new double [nLargest];
	double * buildElement = new double [nelLargest];
	int * buildStart = new int[nLargest+1];
	int * buildRow = new int[nelLargest];
	// alow for integers in knapsacks
	OsiObject ** object = new OsiObject * [numberKnapsack+nTotal];
	int nSOS=0;
	int nObj=numberKnapsack;
	for (iKnapsack=0;iKnapsack<numberKnapsack;iKnapsack++) {
	  knapsackStart[iKnapsack]=finalModel->getNumCols();
	  iRow = knapsackRow[iKnapsack];
	  int nCreate = 10000;
	  coinModel.expandKnapsack(iRow,nCreate,buildObj,buildStart,buildRow,buildElement);
	  // Redo row numbers
	  for (iColumn=0;iColumn<nCreate;iColumn++) {
	    for (int j=buildStart[iColumn];j<buildStart[iColumn+1];j++) {
	      int jRow=buildRow[j];
	      jRow=lookupRow[jRow];
	      assert (jRow>=0&&jRow<nRow);
	      buildRow[j]=jRow;
	    }
	  }
	  finalModel->addCols(nCreate,buildStart,buildRow,buildElement,NULL,NULL,buildObj);
	  int numberFinal=finalModel->getNumCols();
	  for (iColumn=numberOther;iColumn<numberFinal;iColumn++) {
	    if (markKnapsack[iKnapsack]<0) {
	      finalModel->setColUpper(iColumn,maxCoefficient);
	      finalModel->setInteger(iColumn);
	    } else {
	      finalModel->setColUpper(iColumn,maxCoefficient+1.0);
	      finalModel->setInteger(iColumn);
	    }
	    OsiSimpleInteger * sosObject = new OsiSimpleInteger(finalModel,iColumn);
	    sosObject->setPriority(1000000);
	    object[nObj++]=sosObject;
	    buildRow[iColumn-numberOther]=iColumn;
	    buildElement[iColumn-numberOther]=1.0;
	  }
	  if (markKnapsack[iKnapsack]<0) {
	    // convexity row
	    finalModel->addRow(numberFinal-numberOther,buildRow,buildElement,1.0,1.0);
	  } else {
	    int iColumn = markKnapsack[iKnapsack];
	    int n=numberFinal-numberOther;
	    buildRow[n]=iColumn;
	    buildElement[n++]=-fabs(coefficient[iKnapsack]);
	    // convexity row (sort of)
	    finalModel->addRow(n,buildRow,buildElement,0.0,0.0);
	    OsiSOS * sosObject = new OsiSOS(finalModel,n-1,buildRow,NULL,1);
	    sosObject->setPriority(iKnapsack+SOSPriority);
	    // Say not integral even if is (switch off heuristics)
	    sosObject->setIntegerValued(false);
	    object[nSOS++]=sosObject;
	  }
	  numberOther=numberFinal;
	}
	finalModel->addObjects(nObj,object);
	for (iKnapsack=0;iKnapsack<nObj;iKnapsack++) 
	  delete object[iKnapsack];
	delete [] object;
	// Can we move any rows to cuts
	const int * cutMarker = coinModel.cutMarker();
	if (cutMarker&&0) {
	  printf("AMPL CUTS OFF until global cuts fixed\n");
	  cutMarker=NULL;
	}
	if (cutMarker) {
	  // Row copy
	  const CoinPackedMatrix * matrixByRow = finalModel->getMatrixByRow();
	  const double * elementByRow = matrixByRow->getElements();
	  const int * column = matrixByRow->getIndices();
	  const CoinBigIndex * rowStart = matrixByRow->getVectorStarts();
	  const int * rowLength = matrixByRow->getVectorLengths();
	  
	  const double * rowLower = finalModel->getRowLower();
	  const double * rowUpper = finalModel->getRowUpper();
	  int nDelete=0;
	  for (iRow=0;iRow<numberRows;iRow++) {
	    if (cutMarker[iRow]&&lookupRow[iRow]>=0) {
	      int jRow=lookupRow[iRow];
	      whichRow[nDelete++]=jRow;
	      int start = rowStart[jRow];
	      stored.addCut(rowLower[jRow],rowUpper[jRow],
			    rowLength[jRow],column+start,elementByRow+start);
	    }
	  }
	  finalModel->deleteRows(nDelete,whichRow);
	}
	knapsackStart[numberKnapsack]=finalModel->getNumCols();
	delete [] buildObj;
	delete [] buildElement;
	delete [] buildStart;
	delete [] buildRow;
	finalModel->writeMps("full");
      }
    }
  }
  delete [] whichKnapsack;
  delete [] markRow;
  delete [] markKnapsack;
  delete [] coefficient;
  delete [] linear;
  delete [] whichRow;
  delete [] lookupRow;
  delete si;
  si=NULL;
  if (!badModel&&finalModel) {
    finalModel->setDblParam(OsiObjOffset,coinModel.objectiveOffset());
    return finalModel;
  } else {
    delete finalModel;
    printf("can't make knapsacks - did you set fixedPriority (extra1)\n");
    return NULL;
  }
}
// Fills in original solution (coinModel length)
static void 
afterKnapsack(const CoinModel & coinModel2, const int * whichColumn, const int * knapsackStart, 
		   const int * knapsackRow, int numberKnapsack,
		   const double * knapsackSolution, double * solution, int logLevel)
{
  CoinModel coinModel = coinModel2;
  int numberColumns = coinModel.numberColumns();
  int iColumn;
  // associate all columns to stop possible error messages
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    coinModel.associateElement(coinModel.columnName(iColumn),1.0);
  }
  CoinZeroN(solution,numberColumns);
  int nCol=knapsackStart[0];
  for (iColumn=0;iColumn<nCol;iColumn++) {
    int jColumn = whichColumn[iColumn];
    solution[jColumn]=knapsackSolution[iColumn];
  }
  int * buildRow = new int [numberColumns]; // wild overkill
  double * buildElement = new double [numberColumns];
  int iKnapsack;
  for (iKnapsack=0;iKnapsack<numberKnapsack;iKnapsack++) {
    int k=-1;
    double value=0.0;
    for (iColumn=knapsackStart[iKnapsack];iColumn<knapsackStart[iKnapsack+1];iColumn++) {
      if (knapsackSolution[iColumn]>1.0e-5) {
	if (k>=0) {
	  printf("Two nonzero values for knapsack %d at (%d,%g) and (%d,%g)\n",iKnapsack,
		 k,knapsackSolution[k],iColumn,knapsackSolution[iColumn]);
	  abort();
	}
	k=iColumn;
	value=floor(knapsackSolution[iColumn]+0.5);
	assert (fabs(value-knapsackSolution[iColumn])<1.0e-5);
      }
    }
    if (k>=0) {
      int iRow = knapsackRow[iKnapsack];
      int nCreate = 10000;
      int nel=coinModel.expandKnapsack(iRow,nCreate,NULL,NULL,buildRow,buildElement,k-knapsackStart[iKnapsack]);
      assert (nel);
      if (logLevel>0) 
	printf("expanded column %d in knapsack %d has %d nonzero entries:\n",
	       k-knapsackStart[iKnapsack],iKnapsack,nel);
      for (int i=0;i<nel;i++) {
	int jColumn = buildRow[i];
	double value = buildElement[i];
	if (logLevel>0) 
	  printf("%d - original %d has value %g\n",i,jColumn,value);
	solution[jColumn]=value;
      }
    }
  }
  delete [] buildRow;
  delete [] buildElement;
#if 0
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (solution[iColumn]>1.0e-5&&coinModel.isInteger(iColumn))
      printf("%d %g\n",iColumn,solution[iColumn]);
  }
#endif
}
#endif
#if 0
static int outDupRow(OsiSolverInterface * solver) 
{
  CglDuplicateRow dupCuts(solver);
  CglTreeInfo info;
  info.level = 0;
  info.pass = 0;
  int numberRows = solver->getNumRows();
  info.formulation_rows = numberRows;
  info.inTree = false;
  info.strengthenRow= NULL;
  info.pass = 0;
  OsiCuts cs;
  dupCuts.generateCuts(*solver,cs,info);
  const int * duplicate = dupCuts.duplicate();
  // Get rid of duplicate rows
  int * which = new int[numberRows]; 
  int numberDrop=0;
  for (int iRow=0;iRow<numberRows;iRow++) {
    if (duplicate[iRow]==-2||duplicate[iRow]>=0) 
      which[numberDrop++]=iRow;
  }
  if (numberDrop) {
    solver->deleteRows(numberDrop,which);
  }
  delete [] which;
  // see if we have any column cuts
  int numberColumnCuts = cs.sizeColCuts() ;
  const double * columnLower = solver->getColLower();
  const double * columnUpper = solver->getColUpper();
  for (int k = 0;k<numberColumnCuts;k++) {
    OsiColCut * thisCut = cs.colCutPtr(k) ;
    const CoinPackedVector & lbs = thisCut->lbs() ;
    const CoinPackedVector & ubs = thisCut->ubs() ;
    int j ;
    int n ;
    const int * which ;
    const double * values ;
    n = lbs.getNumElements() ;
    which = lbs.getIndices() ;
    values = lbs.getElements() ;
    for (j = 0;j<n;j++) {
      int iColumn = which[j] ;
      if (values[j]>columnLower[iColumn]) 
        solver->setColLower(iColumn,values[j]) ;
    }
    n = ubs.getNumElements() ;
    which = ubs.getIndices() ;
    values = ubs.getElements() ;
    for (j = 0;j<n;j++) {
      int iColumn = which[j] ;
      if (values[j]<columnUpper[iColumn]) 
        solver->setColUpper(iColumn,values[j]) ;
    }
  }
  return numberDrop;
}
#endif
void checkSOS(CbcModel * babModel, const OsiSolverInterface * solver)
{
#ifdef COIN_DEVELOP
  if (!babModel->ownObjects())
    return;
  //const double *objective = solver->getObjCoefficients() ;
  const double *columnLower = solver->getColLower() ;
  const double * columnUpper = solver->getColUpper() ;
  const double * solution = solver->getColSolution();
  //int numberColumns = solver->getNumCols() ;
  //int numberRows = solver->getNumRows();
  //double direction = solver->getObjSense();
  //int iRow,iColumn;

  // Row copy
  CoinPackedMatrix matrixByRow(*solver->getMatrixByRow());
  //const double * elementByRow = matrixByRow.getElements();
  //const int * column = matrixByRow.getIndices();
  //const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();

  // Column copy
  CoinPackedMatrix  matrixByCol(*solver->getMatrixByCol());
  const double * element = matrixByCol.getElements();
  const int * row = matrixByCol.getIndices();
  const CoinBigIndex * columnStart = matrixByCol.getVectorStarts();
  const int * columnLength = matrixByCol.getVectorLengths();

  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();
  OsiObject ** objects = babModel->objects();
  int numberObjects = babModel->numberObjects();
  for (int iObj = 0;iObj<numberObjects;iObj++) {
    CbcSOS * objSOS =
      dynamic_cast <CbcSOS *>(objects[iObj]) ;
    if (objSOS) {
      int n=objSOS->numberMembers();
      const int * which = objSOS->members();
      const double * weight = objSOS->weights();
      int type = objSOS->sosType();
      // convexity row?
      int iColumn;
      iColumn=which[0];
      int j;
      int convex=-1;
      for (j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
	int iRow = row[j];
	double value = element[j];
	if (rowLower[iRow]==1.0&&rowUpper[iRow]==1.0&&
	    value==1.0) {
	  // possible
	  if (rowLength[iRow]==n) {
	    if (convex==-1)
	      convex=iRow;
	    else
	      convex=-2;
	  }
	}
      }
      printf ("set %d of type %d has %d members - possible convexity row %d\n",
	      iObj,type,n,convex);
      for (int i=0;i<n;i++) {
	iColumn = which[i];
	int convex2=-1;
	for (j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  if (iRow==convex) {
	    double value = element[j];
	    if (value==1.0) {
	      convex2=iRow;
	    }
	  }
	}
	if (convex2<0&&convex>=0) {
	  printf("odd convexity row\n");
	  convex=-2;
	}
	printf("col %d has weight %g and value %g, bounds %g %g\n",
	       iColumn,weight[i],solution[iColumn],columnLower[iColumn],
	       columnUpper[iColumn]);
      }
    }
  }
#endif
}
int callCbc1(const char * input2, CbcModel & model)
{
  char * input = strdup(input2);
  int length = strlen(input);
  bool blank = input[0]=='0';
  int n=blank ? 0 : 1;
  for (int i=0;i<length;i++) {
    if (blank) {
      // look for next non blank
      if (input[i]==' ') {
	continue;
      } else {
	n++;
	blank=false;
      }
    } else {
      // look for next blank
      if (input[i]!=' ') {
	continue;
      } else {
	blank=true;
      }
    }
  }
  char ** argv = new char * [n+2];
  argv[0]=strdup("cbc");
  int i=0;
  while(input[i]==' ')
    i++;
  for (int j=0;j<n;j++) {
    int saveI=i;
    for (;i<length;i++) {
      // look for next blank
      if (input[i]!=' ') {
	continue;
      } else {
	break;
      }
    }
    input[i]='\0';
    argv[j+1]=strdup(input+saveI);
    while(input[i]==' ')
      i++;
  }
  argv[n+1]=strdup("-quit");
  free(input);
  totalTime=0.0;
  currentBranchModel = NULL;
  CbcOrClpRead_mode=1;
  CbcOrClpReadCommand=stdin;
  noPrinting=false;
  int returnCode = CbcMain1(n+2,const_cast<const char **>(argv),model);
  for (int k=0;k<n+2;k++)
    free(argv[k]);
  delete [] argv;
  return returnCode;
}
int callCbc1(const std::string input2, CbcModel & babSolver)
{
  char * input3 = strdup(input2.c_str());
  int returnCode=callCbc1(input3,babSolver);
  free(input3);
  return returnCode;
}
int callCbc(const char * input2, CbcModel & babSolver)
{
  CbcMain0(babSolver);
  return callCbc1(input2,babSolver);
}
int callCbc(const std::string input2, CbcModel & babSolver)
{
  char * input3 = strdup(input2.c_str());
  CbcMain0(babSolver);
  int returnCode=callCbc1(input3,babSolver);
  free(input3);
  return returnCode;
}
int callCbc(const char * input2, OsiClpSolverInterface& solver1) 
{
  CbcModel model(solver1);
  return callCbc(input2,model);
}
int callCbc(const char * input2)
{
  {
    OsiClpSolverInterface solver1;
    return callCbc(input2,solver1);
  }
}
int callCbc(const std::string input2, OsiClpSolverInterface& solver1) 
{
  char * input3 = strdup(input2.c_str());
  int returnCode=callCbc(input3,solver1);
  free(input3);
  return returnCode;
}
int callCbc(const std::string input2) 
{
  char * input3 = strdup(input2.c_str());
  OsiClpSolverInterface solver1;
  int returnCode=callCbc(input3,solver1);
  free(input3);
  return returnCode;
}

int CbcMain (int argc, const char *argv[],
	     CbcModel  & model)
{
  CbcMain0(model);
  return CbcMain1(argc,argv,model);
}
#define CBCMAXPARAMETERS 200
static CbcOrClpParam parameters[CBCMAXPARAMETERS];
static int numberParameters=0 ;
void CbcMain0 (CbcModel  & model)
{
  OsiClpSolverInterface * originalSolver = dynamic_cast<OsiClpSolverInterface *> (model.solver());
  assert (originalSolver);
  CoinMessageHandler * generalMessageHandler = originalSolver->messageHandler();
  generalMessageHandler->setPrefix(true);
  OsiSolverInterface * solver = model.solver();
  OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
  ClpSimplex * lpSolver = clpSolver->getModelPtr();
  lpSolver->setPerturbation(50);
  lpSolver->messageHandler()->setPrefix(false);
  establishParams(numberParameters,parameters) ;
  const char dirsep =  CoinFindDirSeparator();
  std::string directory = (dirsep == '/' ? "./" : ".\\");
  std::string defaultDirectory = directory;
  std::string importFile ="";
  std::string exportFile ="default.mps";
  std::string importBasisFile ="";
  std::string importPriorityFile ="";
  std::string debugFile="";
  std::string printMask="";
  std::string exportBasisFile ="default.bas";
  std::string saveFile ="default.prob";
  std::string restoreFile ="default.prob";
  std::string solutionFile ="stdout";
  std::string solutionSaveFile ="solution.file";
  int doIdiot=-1;
  int outputFormat=2;
  int substitution=3;
  int dualize=0;
  int preSolve=5;
  int doSprint=-1;
  int testOsiParameters=-1;
  parameters[whichParam(BASISIN,numberParameters,parameters)].setStringValue(importBasisFile);
  parameters[whichParam(PRIORITYIN,numberParameters,parameters)].setStringValue(importPriorityFile);
  parameters[whichParam(BASISOUT,numberParameters,parameters)].setStringValue(exportBasisFile);
  parameters[whichParam(DEBUG,numberParameters,parameters)].setStringValue(debugFile);
  parameters[whichParam(PRINTMASK,numberParameters,parameters)].setStringValue(printMask);
  parameters[whichParam(DIRECTORY,numberParameters,parameters)].setStringValue(directory);
  parameters[whichParam(DUALBOUND,numberParameters,parameters)].setDoubleValue(lpSolver->dualBound());
  parameters[whichParam(DUALTOLERANCE,numberParameters,parameters)].setDoubleValue(lpSolver->dualTolerance());
  parameters[whichParam(EXPORT,numberParameters,parameters)].setStringValue(exportFile);
  parameters[whichParam(IDIOT,numberParameters,parameters)].setIntValue(doIdiot);
  parameters[whichParam(IMPORT,numberParameters,parameters)].setStringValue(importFile);
  parameters[whichParam(PRESOLVETOLERANCE,numberParameters,parameters)].setDoubleValue(1.0e-8);
  int slog = whichParam(SOLVERLOGLEVEL,numberParameters,parameters);
  int log = whichParam(LOGLEVEL,numberParameters,parameters);
  parameters[slog].setIntValue(1);
  clpSolver->messageHandler()->setLogLevel(1) ;
  model.messageHandler()->setLogLevel(1);
  lpSolver->setLogLevel(1);
  parameters[log].setIntValue(1);
  parameters[whichParam(MAXFACTOR,numberParameters,parameters)].setIntValue(lpSolver->factorizationFrequency());
  parameters[whichParam(MAXITERATION,numberParameters,parameters)].setIntValue(lpSolver->maximumIterations());
  parameters[whichParam(OUTPUTFORMAT,numberParameters,parameters)].setIntValue(outputFormat);
  parameters[whichParam(PRESOLVEPASS,numberParameters,parameters)].setIntValue(preSolve);
  parameters[whichParam(PERTVALUE,numberParameters,parameters)].setIntValue(lpSolver->perturbation());
  parameters[whichParam(PRIMALTOLERANCE,numberParameters,parameters)].setDoubleValue(lpSolver->primalTolerance());
  parameters[whichParam(PRIMALWEIGHT,numberParameters,parameters)].setDoubleValue(lpSolver->infeasibilityCost());
  parameters[whichParam(RESTORE,numberParameters,parameters)].setStringValue(restoreFile);
  parameters[whichParam(SAVE,numberParameters,parameters)].setStringValue(saveFile);
  //parameters[whichParam(TIMELIMIT,numberParameters,parameters)].setDoubleValue(1.0e8);
  parameters[whichParam(TIMELIMIT_BAB,numberParameters,parameters)].setDoubleValue(1.0e8);
  parameters[whichParam(SOLUTION,numberParameters,parameters)].setStringValue(solutionFile);
  parameters[whichParam(SAVESOL,numberParameters,parameters)].setStringValue(solutionSaveFile);
  parameters[whichParam(SPRINT,numberParameters,parameters)].setIntValue(doSprint);
  parameters[whichParam(SUBSTITUTION,numberParameters,parameters)].setIntValue(substitution);
  parameters[whichParam(DUALIZE,numberParameters,parameters)].setIntValue(dualize);
  model.setNumberBeforeTrust(10);
  parameters[whichParam(NUMBERBEFORE,numberParameters,parameters)].setIntValue(5);
  parameters[whichParam(MAXNODES,numberParameters,parameters)].setIntValue(model.getMaximumNodes());
  model.setNumberStrong(5);
  parameters[whichParam(STRONGBRANCHING,numberParameters,parameters)].setIntValue(model.numberStrong());
  parameters[whichParam(INFEASIBILITYWEIGHT,numberParameters,parameters)].setDoubleValue(model.getDblParam(CbcModel::CbcInfeasibilityWeight));
  parameters[whichParam(INTEGERTOLERANCE,numberParameters,parameters)].setDoubleValue(model.getDblParam(CbcModel::CbcIntegerTolerance));
  parameters[whichParam(INCREMENT,numberParameters,parameters)].setDoubleValue(model.getDblParam(CbcModel::CbcCutoffIncrement));
  parameters[whichParam(TESTOSI,numberParameters,parameters)].setIntValue(testOsiParameters);
  parameters[whichParam(FPUMPTUNE,numberParameters,parameters)].setIntValue(1003);
#ifdef CBC_THREAD
  parameters[whichParam(THREADS,numberParameters,parameters)].setIntValue(0);
#endif
  // Set up likely cut generators and defaults
  parameters[whichParam(PREPROCESS,numberParameters,parameters)].setCurrentOption("on");
  parameters[whichParam(MIPOPTIONS,numberParameters,parameters)].setIntValue(128|64|1);
  parameters[whichParam(MIPOPTIONS,numberParameters,parameters)].setIntValue(1);
  parameters[whichParam(CUTPASSINTREE,numberParameters,parameters)].setIntValue(1);
  parameters[whichParam(MOREMIPOPTIONS,numberParameters,parameters)].setIntValue(-1);
  parameters[whichParam(MAXHOTITS,numberParameters,parameters)].setIntValue(100);
  parameters[whichParam(CUTSSTRATEGY,numberParameters,parameters)].setCurrentOption("on");
  parameters[whichParam(HEURISTICSTRATEGY,numberParameters,parameters)].setCurrentOption("on");
  parameters[whichParam(NODESTRATEGY,numberParameters,parameters)].setCurrentOption("fewest");
  parameters[whichParam(GOMORYCUTS,numberParameters,parameters)].setCurrentOption("ifmove");
  parameters[whichParam(PROBINGCUTS,numberParameters,parameters)].setCurrentOption("ifmove");
  parameters[whichParam(KNAPSACKCUTS,numberParameters,parameters)].setCurrentOption("ifmove");
  parameters[whichParam(REDSPLITCUTS,numberParameters,parameters)].setCurrentOption("off");
  parameters[whichParam(CLIQUECUTS,numberParameters,parameters)].setCurrentOption("ifmove");
  parameters[whichParam(MIXEDCUTS,numberParameters,parameters)].setCurrentOption("ifmove");
  parameters[whichParam(FLOWCUTS,numberParameters,parameters)].setCurrentOption("ifmove");
  parameters[whichParam(TWOMIRCUTS,numberParameters,parameters)].setCurrentOption("root");
  parameters[whichParam(LANDPCUTS,numberParameters,parameters)].setCurrentOption("off");
  parameters[whichParam(RESIDCUTS,numberParameters,parameters)].setCurrentOption("off");
  parameters[whichParam(ROUNDING,numberParameters,parameters)].setCurrentOption("on");
  parameters[whichParam(FPUMP,numberParameters,parameters)].setCurrentOption("on");
  parameters[whichParam(GREEDY,numberParameters,parameters)].setCurrentOption("on");
  parameters[whichParam(COMBINE,numberParameters,parameters)].setCurrentOption("on");
  parameters[whichParam(RINS,numberParameters,parameters)].setCurrentOption("off");
  parameters[whichParam(COSTSTRATEGY,numberParameters,parameters)].setCurrentOption("off");
}
int CbcMain1 (int argc, const char *argv[],
	     CbcModel  & model)
{
  /* Note
     This is meant as a stand-alone executable to do as much of coin as possible. 
     It should only have one solver known to it.
  */
  OsiClpSolverInterface * originalSolver = dynamic_cast<OsiClpSolverInterface *> (model.solver());
  assert (originalSolver);
  CoinMessageHandler * generalMessageHandler = model.messageHandler();
  generalMessageHandler->setPrefix(false);
  // Move handler across if not default
  if (!originalSolver->defaultHandler()&&originalSolver->getModelPtr()->defaultHandler())
    originalSolver->getModelPtr()->passInMessageHandler(originalSolver->messageHandler());
  CoinMessages generalMessages = originalSolver->getModelPtr()->messages();
  char generalPrint[10000];
  if (originalSolver->getModelPtr()->logLevel()==0)
    noPrinting=true;
  // see if log in list
  for (int i=1;i<argc;i++) {
    if (!strncmp(argv[i],"log",3)) {
      const char * equals = strchr(argv[i],'=');
      if (equals&&atoi(equals+1)>0) 
	noPrinting=false;
      else
	noPrinting=true;
      break;
    } else if (!strncmp(argv[i],"-log",4)&&i<argc-1) {
      if (atoi(argv[i+1])>0) 
	noPrinting=false;
      else
	noPrinting=true;
      break;
    }
  }
  CbcModel * babModel = NULL;
  {
    double time1 = CoinCpuTime(),time2;
    bool goodModel=(originalSolver->getNumCols()) ? true : false;

    CoinSighandler_t saveSignal=SIG_DFL;
    // register signal handler
    saveSignal = signal(SIGINT,signal_handler);
    // Set up all non-standard stuff
    int cutPass=-1234567;
    int cutPassInTree=-1234567;
    int tunePreProcess=5;
    int testOsiParameters=-1;
    // 0 normal, 1 from ampl or MIQP etc (2 allows cuts)
    int complicatedInteger=0;
    OsiSolverInterface * solver = model.solver();
    OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
    ClpSimplex * lpSolver = clpSolver->getModelPtr();
    if (noPrinting) {
      setCbcOrClpPrinting(false);
      lpSolver->setLogLevel(0);
    }
    // For priorities etc
    int * priorities=NULL;
    int * branchDirection=NULL;
    double * pseudoDown=NULL;
    double * pseudoUp=NULL;
    double * solutionIn = NULL;
    int * prioritiesIn = NULL;
    int numberSOS = 0;
    int * sosStart = NULL;
    int * sosIndices = NULL;
    char * sosType = NULL;
    double * sosReference = NULL;
    int * cut=NULL;
    int * sosPriority=NULL;
    CglStored storedAmpl;
    CoinModel * coinModel = NULL;
#ifdef COIN_HAS_ASL
    ampl_info info;
    CoinModel saveCoinModel;
    CoinModel saveTightenedModel;
    int * whichColumn = NULL;
    int * knapsackStart=NULL;
    int * knapsackRow=NULL;
    int numberKnapsack=0;
    memset(&info,0,sizeof(info));
    if (argc>2&&!strcmp(argv[2],"-AMPL")) {
      usingAmpl=true;
      // see if log in list
      noPrinting=true;
      for (int i=1;i<argc;i++) {
        if (!strncmp(argv[i],"log",3)) {
	  const char * equals = strchr(argv[i],'=');
          if (equals&&atoi(equals+1)>0) {
            noPrinting=false;
	    info.logLevel=atoi(equals+1);
	    int log = whichParam(LOGLEVEL,numberParameters,parameters);
	    parameters[log].setIntValue(info.logLevel);
	    // mark so won't be overWritten
	    info.numberRows=-1234567;
	    break;
	  }
        }
      }
      int returnCode = readAmpl(&info,argc,const_cast<char **>(argv),(void **) (& coinModel));
      if (returnCode)
        return returnCode;
      CbcOrClpRead_mode=2; // so will start with parameters
      // see if log in list (including environment)
      for (int i=1;i<info.numberArguments;i++) {
        if (!strcmp(info.arguments[i],"log")) {
          if (i<info.numberArguments-1&&atoi(info.arguments[i+1])>0)
            noPrinting=false;
          break;
        }
      }
      if (noPrinting) {
        model.messageHandler()->setLogLevel(0);
        setCbcOrClpPrinting(false);
      }
      if (!noPrinting)
        printf("%d rows, %d columns and %d elements\n",
               info.numberRows,info.numberColumns,info.numberElements);
#ifdef COIN_HAS_LINK
      if (!coinModel) {
#endif
      solver->loadProblem(info.numberColumns,info.numberRows,info.starts,
                          info.rows,info.elements,
                          info.columnLower,info.columnUpper,info.objective,
                          info.rowLower,info.rowUpper);
      // take off cuts if ampl wants that
      if (info.cut&&0) {
	printf("AMPL CUTS OFF until global cuts fixed\n");
	info.cut=NULL;
      }
      if (info.cut) {
	int numberRows = info.numberRows;
	int * whichRow = new int [numberRows];
	// Row copy
	const CoinPackedMatrix * matrixByRow = solver->getMatrixByRow();
	const double * elementByRow = matrixByRow->getElements();
	const int * column = matrixByRow->getIndices();
	const CoinBigIndex * rowStart = matrixByRow->getVectorStarts();
	const int * rowLength = matrixByRow->getVectorLengths();
	
	const double * rowLower = solver->getRowLower();
	const double * rowUpper = solver->getRowUpper();
	int nDelete=0;
	for (int iRow=0;iRow<numberRows;iRow++) {
	  if (info.cut[iRow]) {
	    whichRow[nDelete++]=iRow;
	    int start = rowStart[iRow];
	    storedAmpl.addCut(rowLower[iRow],rowUpper[iRow],
			  rowLength[iRow],column+start,elementByRow+start);
	  }
	}
	solver->deleteRows(nDelete,whichRow);
	delete [] whichRow;
      }
#ifdef COIN_HAS_LINK
      } else {
	// save
	saveCoinModel = *coinModel;
	// load from coin model
	OsiSolverLink solver1;
	OsiSolverInterface * solver2 = solver1.clone();
	model.assignSolver(solver2,false);
	OsiSolverLink * si =
	  dynamic_cast<OsiSolverLink *>(model.solver()) ;
	assert (si != NULL);
	si->setDefaultMeshSize(0.001);
	// need some relative granularity
	si->setDefaultBound(100.0);
	si->setDefaultMeshSize(0.001);
	si->setDefaultBound(100000.0);
	si->setIntegerPriority(1000);
	si->setBiLinearPriority(10000);
	CoinModel * model2 = (CoinModel *) coinModel;
	int logLevel = whichParam(LOGLEVEL,numberParameters,parameters);
	si->load(*model2,true,logLevel);
	// redo
	solver = model.solver();
	clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
	lpSolver = clpSolver->getModelPtr();
	clpSolver->messageHandler()->setLogLevel(0) ;
	testOsiParameters=0;
	parameters[whichParam(TESTOSI,numberParameters,parameters)].setIntValue(0);
	complicatedInteger=1;
	if (info.cut) {
	  printf("Sorry - can't do cuts with LOS as ruins delicate row order\n");
	  abort();
	  int numberRows = info.numberRows;
	  int * whichRow = new int [numberRows];
	  // Row copy
	  const CoinPackedMatrix * matrixByRow = solver->getMatrixByRow();
	  const double * elementByRow = matrixByRow->getElements();
	  const int * column = matrixByRow->getIndices();
	  const CoinBigIndex * rowStart = matrixByRow->getVectorStarts();
	  const int * rowLength = matrixByRow->getVectorLengths();
	  
	  const double * rowLower = solver->getRowLower();
	  const double * rowUpper = solver->getRowUpper();
	  int nDelete=0;
	  for (int iRow=0;iRow<numberRows;iRow++) {
	    if (info.cut[iRow]) {
	      whichRow[nDelete++]=iRow;
	      int start = rowStart[iRow];
	      storedAmpl.addCut(rowLower[iRow],rowUpper[iRow],
				rowLength[iRow],column+start,elementByRow+start);
	    }
	  }
	  solver->deleteRows(nDelete,whichRow);
	  // and special matrix
	  si->cleanMatrix()->deleteRows(nDelete,whichRow);
	  delete [] whichRow;
	}
      }
#endif
      // If we had a solution use it
      if (info.primalSolution) {
        solver->setColSolution(info.primalSolution);
      }
      // status
      if (info.rowStatus) {
        unsigned char * statusArray = lpSolver->statusArray();
        int i;
        for (i=0;i<info.numberColumns;i++)
          statusArray[i]=(char)info.columnStatus[i];
        statusArray+=info.numberColumns;
        for (i=0;i<info.numberRows;i++)
          statusArray[i]=(char)info.rowStatus[i];
        CoinWarmStartBasis * basis = lpSolver->getBasis();
        solver->setWarmStart(basis);
        delete basis;
      }
      freeArrays1(&info);
      // modify objective if necessary
      solver->setObjSense(info.direction);
      solver->setDblParam(OsiObjOffset,info.offset);
      if (info.offset) {
	sprintf(generalPrint,"Ampl objective offset is %g",
		info.offset);
	generalMessageHandler->message(CLP_GENERAL,generalMessages)
	  << generalPrint
	  <<CoinMessageEol;
      }
      // Set integer variables (unless nonlinear when set)
      if (!info.nonLinear) {
	for (int i=info.numberColumns-info.numberIntegers;
	     i<info.numberColumns;i++)
	  solver->setInteger(i);
      }
      goodModel=true;
      // change argc etc
      argc = info.numberArguments;
      argv = const_cast<const char **>(info.arguments);
    }
#endif    
    // default action on import
    int allowImportErrors=0;
    int keepImportNames=1;
    int doIdiot=-1;
    int outputFormat=2;
    int slpValue=-1;
    int cppValue=-1;
    int printOptions=0;
    int printMode=0;
    int presolveOptions=0;
    int substitution=3;
    int dualize=0;
    int doCrash=0;
    int doVector=0;
    int doSprint=-1;
    int doScaling=1;
    // set reasonable defaults
    int preSolve=5;
    int preProcess=1;
    bool useStrategy=false;
    bool preSolveFile=false;
    bool strongChanged=false;
   
    double djFix=1.0e100;
    double gapRatio=1.0e100;
    double tightenFactor=0.0;
    const char dirsep =  CoinFindDirSeparator();
    std::string directory = (dirsep == '/' ? "./" : ".\\");
    std::string defaultDirectory = directory;
    std::string importFile ="";
    std::string exportFile ="default.mps";
    std::string importBasisFile ="";
    std::string importPriorityFile ="";
    std::string debugFile="";
    std::string printMask="";
    double * debugValues = NULL;
    int numberDebugValues = -1;
    int basisHasValues=0;
    std::string exportBasisFile ="default.bas";
    std::string saveFile ="default.prob";
    std::string restoreFile ="default.prob";
    std::string solutionFile ="stdout";
    std::string solutionSaveFile ="solution.file";
    int slog = whichParam(SOLVERLOGLEVEL,numberParameters,parameters);
    int log = whichParam(LOGLEVEL,numberParameters,parameters);
    double normalIncrement=model.getCutoffIncrement();;
    if (testOsiParameters>=0) {
      // trying nonlinear - switch off some stuff
      preProcess=0;
    }
    // Set up likely cut generators and defaults
    int nodeStrategy=0;
    int doSOS=1;
    int verbose=0;
    CglGomory gomoryGen;
    // try larger limit
    gomoryGen.setLimitAtRoot(512);
    gomoryGen.setLimit(50);
    // set default action (0=off,1=on,2=root)
    int gomoryAction=3;

    CglProbing probingGen;
    probingGen.setUsingObjective(1);
    probingGen.setMaxPass(3);
    probingGen.setMaxPassRoot(3);
    // Number of unsatisfied variables to look at
    probingGen.setMaxProbe(10);
    probingGen.setMaxProbeRoot(50);
    // How far to follow the consequences
    probingGen.setMaxLook(10);
    probingGen.setMaxLookRoot(50);
    probingGen.setMaxLookRoot(10);
    // Only look at rows with fewer than this number of elements
    probingGen.setMaxElements(200);
    probingGen.setRowCuts(3);
    // set default action (0=off,1=on,2=root)
    int probingAction=1;

    CglKnapsackCover knapsackGen;
    //knapsackGen.switchOnExpensive();
    // set default action (0=off,1=on,2=root)
    int knapsackAction=3;

    CglRedSplit redsplitGen;
    //redsplitGen.setLimit(100);
    // set default action (0=off,1=on,2=root)
    // Off as seems to give some bad cuts
    int redsplitAction=0;

    CglClique cliqueGen(false,true);
    cliqueGen.setStarCliqueReport(false);
    cliqueGen.setRowCliqueReport(false);
    cliqueGen.setMinViolation(0.1);
    // set default action (0=off,1=on,2=root)
    int cliqueAction=3;

    CglMixedIntegerRounding2 mixedGen;
    // set default action (0=off,1=on,2=root)
    int mixedAction=3;

    CglFlowCover flowGen;
    // set default action (0=off,1=on,2=root)
    int flowAction=3;

    CglTwomir twomirGen;
    twomirGen.setMaxElements(250);
    // set default action (0=off,1=on,2=root)
    int twomirAction=2;
    CglLandP landpGen;
    // set default action (0=off,1=on,2=root)
    int landpAction=0;
    CglResidualCapacity residualCapacityGen;
    // set default action (0=off,1=on,2=root)
    int residualCapacityAction=0;
    // Stored cuts
    bool storedCuts = false;

    bool useRounding=true;
    bool useFpump=true;
    bool useGreedy=true;
    bool useCombine=true;
    bool useLocalTree=false;
    bool useRINS=false;
    int useCosts=0;
    // don't use input solution
    int useSolution=0;
    
    // total number of commands read
    int numberGoodCommands=0;
    // Set false if user does anything advanced
    bool defaultSettings=true;

    // Hidden stuff for barrier
    int choleskyType = 0;
    int gamma=0;
    int scaleBarrier=0;
    int doKKT=0;
    int crossover=2; // do crossover unless quadratic
    // For names
    int lengthName = 0;
    std::vector<std::string> rowNames;
    std::vector<std::string> columnNames;
    
    std::string field;
    if (!noPrinting) {
      sprintf(generalPrint,"Coin Cbc and Clp Solver version %s, build %s",
	      CBCVERSION,__DATE__);
      generalMessageHandler->message(CLP_GENERAL,generalMessages)
	<< generalPrint
	<<CoinMessageEol;
      // Print command line
      if (argc>1) {
        sprintf(generalPrint,"command line - ");
        for (int i=0;i<argc;i++)
          sprintf(generalPrint+strlen(generalPrint),"%s ",argv[i]);
	generalMessageHandler->message(CLP_GENERAL,generalMessages)
	  << generalPrint
	  <<CoinMessageEol;
      }
    }
    while (1) {
      // next command
      field=CoinReadGetCommand(argc,argv);
      // adjust field if has odd trailing characters
      char temp [200];
      strcpy(temp,field.c_str());
      int length = strlen(temp);
      for (int k=length-1;k>=0;k--) {
	if (temp[k]<' ')
	  length--;
	else
	  break;
      }
      temp[length]='\0';
      field=temp;
      // exit if null or similar
      if (!field.length()) {
	if (numberGoodCommands==1&&goodModel) {
	  // we just had file name - do branch and bound
	  field="branch";
	} else if (!numberGoodCommands) {
	  // let's give the sucker a hint
	  std::cout
	    <<"CoinSolver takes input from arguments ( - switches to stdin)"
	    <<std::endl
	    <<"Enter ? for list of commands or help"<<std::endl;
	  field="-";
	} else {
	  break;
	}
      }
      
      // see if ? at end
      int numberQuery=0;
      if (field!="?"&&field!="???") {
	int length = field.length();
	int i;
	for (i=length-1;i>0;i--) {
	  if (field[i]=='?') 
	    numberQuery++;
	  else
	    break;
	}
	field=field.substr(0,length-numberQuery);
      }
      // find out if valid command
      int iParam;
      int numberMatches=0;
      int firstMatch=-1;
      for ( iParam=0; iParam<numberParameters; iParam++ ) {
	int match = parameters[iParam].matches(field);
	if (match==1) {
	  numberMatches = 1;
	  firstMatch=iParam;
	  break;
	} else {
	  if (match&&firstMatch<0)
	    firstMatch=iParam;
	  numberMatches += match>>1;
	}
      }
      if (iParam<numberParameters&&!numberQuery) {
	// found
	CbcOrClpParam found = parameters[iParam];
	CbcOrClpParameterType type = found.type();
	int valid;
	numberGoodCommands++;
	if (type==BAB&&goodModel) {
	  // check if any integers
#ifdef COIN_HAS_ASL
	  if (info.numberSos&&doSOS&&usingAmpl) {
	    // SOS
	    numberSOS = info.numberSos;
	  }
#endif
	  if (!lpSolver->integerInformation()&&!numberSOS&&
	      !clpSolver->numberSOS()&&!model.numberObjects()&&!clpSolver->numberObjects())
	    type=DUALSIMPLEX;
	}
	if (type==GENERALQUERY) {
	  bool evenHidden=false;
	  if ((verbose&8)!=0) {
	    // even hidden
	    evenHidden = true;
	    verbose &= ~8;
	  }
#ifdef COIN_HAS_ASL
          if (verbose<4&&usingAmpl)
            verbose +=4;
#endif
          if (verbose<4) {
            std::cout<<"In argument list keywords have leading - "
              ", -stdin or just - switches to stdin"<<std::endl;
            std::cout<<"One command per line (and no -)"<<std::endl;
            std::cout<<"abcd? gives list of possibilities, if only one + explanation"<<std::endl;
            std::cout<<"abcd?? adds explanation, if only one fuller help"<<std::endl;
            std::cout<<"abcd without value (where expected) gives current value"<<std::endl;
            std::cout<<"abcd value sets value"<<std::endl;
            std::cout<<"Commands are:"<<std::endl;
          } else {
            std::cout<<"Cbc options are set within AMPL with commands like:"<<std::endl<<std::endl;
            std::cout<<"         option cbc_options \"cuts=root log=2 feas=on slog=1\""<<std::endl<<std::endl;
            std::cout<<"only maximize, dual, primal, help and quit are recognized without ="<<std::endl;
          }
	  int maxAcross=5;
          if ((verbose%4)!=0)
            maxAcross=1;
	  int limits[]={1,51,101,151,201,251,301,351,401};
	  std::vector<std::string> types;
	  types.push_back("Double parameters:");
	  types.push_back("Branch and Cut double parameters:");
	  types.push_back("Integer parameters:");
	  types.push_back("Branch and Cut integer parameters:");
	  types.push_back("Keyword parameters:");
	  types.push_back("Branch and Cut keyword parameters:");
	  types.push_back("Actions or string parameters:");
	  types.push_back("Branch and Cut actions:");
	  int iType;
	  for (iType=0;iType<8;iType++) {
	    int across=0;
            if ((verbose%4)!=0)
              std::cout<<std::endl;
	    std::cout<<types[iType]<<std::endl;
            if ((verbose&2)!=0)
              std::cout<<std::endl;
	    for ( iParam=0; iParam<numberParameters; iParam++ ) {
	      int type = parameters[iParam].type();
	      if ((parameters[iParam].displayThis()||evenHidden)&&
		  type>=limits[iType]
		  &&type<limits[iType+1]) {
                // but skip if not useful for ampl (and in ampl mode)
                if (verbose>=4&&(parameters[iParam].whereUsed()&4)==0)
                  continue;
		if (!across) {
                  if ((verbose&2)==0) 
                    std::cout<<"  ";
                  else
                    std::cout<<"Command ";
                }
                std::cout<<parameters[iParam].matchName()<<"  ";
		across++;
		if (across==maxAcross) {
		  across=0;
                  if ((verbose%4)!=0) {
                    // put out description as well
                    if ((verbose&1)!=0) 
                      std::cout<<parameters[iParam].shortHelp();
                    std::cout<<std::endl;
                    if ((verbose&2)!=0) {
                      std::cout<<"---- description"<<std::endl;
                      parameters[iParam].printLongHelp();
                      std::cout<<"----"<<std::endl<<std::endl;
                    }
                  } else {
                    std::cout<<std::endl;
                  }
		}
	      }
	    }
	    if (across)
	      std::cout<<std::endl;
	  }
	} else if (type==FULLGENERALQUERY) {
	  std::cout<<"Full list of commands is:"<<std::endl;
	  int maxAcross=5;
	  int limits[]={1,51,101,151,201,251,301,351,401};
	  std::vector<std::string> types;
	  types.push_back("Double parameters:");
	  types.push_back("Branch and Cut double parameters:");
	  types.push_back("Integer parameters:");
	  types.push_back("Branch and Cut integer parameters:");
	  types.push_back("Keyword parameters:");
	  types.push_back("Branch and Cut keyword parameters:");
	  types.push_back("Actions or string parameters:");
	  types.push_back("Branch and Cut actions:");
	  int iType;
	  for (iType=0;iType<8;iType++) {
	    int across=0;
	    std::cout<<types[iType]<<"  ";
	    for ( iParam=0; iParam<numberParameters; iParam++ ) {
	      int type = parameters[iParam].type();
	      if (type>=limits[iType]
		  &&type<limits[iType+1]) {
		if (!across)
		  std::cout<<"  ";
		std::cout<<parameters[iParam].matchName()<<"  ";
		across++;
		if (across==maxAcross) {
		  std::cout<<std::endl;
		  across=0;
		}
	      }
	    }
	    if (across)
	      std::cout<<std::endl;
	  }
	} else if (type<101) {
	  // get next field as double
	  double value = CoinReadGetDoubleField(argc,argv,&valid);
	  if (!valid) {
	    if (type<51) {
	      parameters[iParam].setDoubleParameter(lpSolver,value);
	    } else if (type<81) {
	      parameters[iParam].setDoubleParameter(model,value);
	    } else {
	      parameters[iParam].setDoubleParameter(lpSolver,value);
	      switch(type) {
	      case DJFIX:
		djFix=value;
		preSolve=5;
                defaultSettings=false; // user knows what she is doing
		break;
	      case GAPRATIO:
		gapRatio=value;
		break;
	      case TIGHTENFACTOR:
		tightenFactor=value;
		if(!complicatedInteger)
		  defaultSettings=false; // user knows what she is doing
		break;
	      default:
		break;
	      }
	    }
	  } else if (valid==1) {
	    abort();
	  } else {
	    std::cout<<parameters[iParam].name()<<" has value "<<
	      parameters[iParam].doubleValue()<<std::endl;
	  }
	} else if (type<201) {
	  // get next field as int
	  int value = CoinReadGetIntField(argc,argv,&valid);
	  if (!valid) {
	    if (type<151) {
	      if (parameters[iParam].type()==PRESOLVEPASS)
		preSolve = value;
	      else if (parameters[iParam].type()==IDIOT)
		doIdiot = value;
	      else if (parameters[iParam].type()==SPRINT)
		doSprint = value;
	      else if (parameters[iParam].type()==OUTPUTFORMAT)
		outputFormat = value;
	      else if (parameters[iParam].type()==SLPVALUE)
		slpValue = value;
              else if (parameters[iParam].type()==CPP)
	        cppValue = value;
	      else if (parameters[iParam].type()==PRESOLVEOPTIONS)
		presolveOptions = value;
	      else if (parameters[iParam].type()==PRINTOPTIONS)
		printOptions = value;
              else if (parameters[iParam].type()==SUBSTITUTION)
                substitution = value;
              else if (parameters[iParam].type()==DUALIZE)
                dualize = value;
	      else if (parameters[iParam].type()==PROCESSTUNE)
		tunePreProcess = value;
	      else if (parameters[iParam].type()==VERBOSE)
		verbose = value;
              parameters[iParam].setIntParameter(lpSolver,value);
	    } else {
	      if (parameters[iParam].type()==CUTPASS)
		cutPass = value;
	      else if (parameters[iParam].type()==CUTPASSINTREE)
		cutPassInTree = value;
	      else if (parameters[iParam].type()==FPUMPITS)
		{ useFpump = true;parameters[iParam].setIntValue(value);}
	      else if (parameters[iParam].type()==STRONGBRANCHING||
		       parameters[iParam].type()==NUMBERBEFORE)
		strongChanged=true;
	      parameters[iParam].setIntParameter(model,value);
	    }
	  } else if (valid==1) {
	    abort();
	  } else {
	    std::cout<<parameters[iParam].name()<<" has value "<<
	      parameters[iParam].intValue()<<std::endl;
	  }
	} else if (type<301) {
	  // one of several strings
	  std::string value = CoinReadGetString(argc,argv);
	  int action = parameters[iParam].parameterOption(value);
	  if (action<0) {
	    if (value!="EOL") {
	      // no match
	      parameters[iParam].printOptions();
	    } else {
	      // print current value
	      std::cout<<parameters[iParam].name()<<" has value "<<
		parameters[iParam].currentOption()<<std::endl;
	    }
	  } else {
	    parameters[iParam].setCurrentOption(action,!noPrinting);
	    // for now hard wired
	    switch (type) {
	    case DIRECTION:
	      if (action==0)
		lpSolver->setOptimizationDirection(1);
	      else if (action==1)
		lpSolver->setOptimizationDirection(-1);
	      else
		lpSolver->setOptimizationDirection(0);
	      break;
	    case DUALPIVOT:
	      if (action==0) {
		ClpDualRowSteepest steep(3);
		lpSolver->setDualRowPivotAlgorithm(steep);
	      } else if (action==1) {
		//ClpDualRowDantzig dantzig;
		ClpDualRowSteepest dantzig(5);
		lpSolver->setDualRowPivotAlgorithm(dantzig);
	      } else if (action==2) {
		// partial steep
		ClpDualRowSteepest steep(2);
		lpSolver->setDualRowPivotAlgorithm(steep);
	      } else {
		ClpDualRowSteepest steep;
		lpSolver->setDualRowPivotAlgorithm(steep);
	      }
	      break;
	    case PRIMALPIVOT:
	      if (action==0) {
		ClpPrimalColumnSteepest steep(3);
		lpSolver->setPrimalColumnPivotAlgorithm(steep);
	      } else if (action==1) {
		ClpPrimalColumnSteepest steep(0);
		lpSolver->setPrimalColumnPivotAlgorithm(steep);
	      } else if (action==2) {
		ClpPrimalColumnDantzig dantzig;
		lpSolver->setPrimalColumnPivotAlgorithm(dantzig);
	      } else if (action==3) {
		ClpPrimalColumnSteepest steep(2);
		lpSolver->setPrimalColumnPivotAlgorithm(steep);
	      } else if (action==4) {
		ClpPrimalColumnSteepest steep(1);
		lpSolver->setPrimalColumnPivotAlgorithm(steep);
	      } else if (action==5) {
		ClpPrimalColumnSteepest steep(4);
		lpSolver->setPrimalColumnPivotAlgorithm(steep);
	      } else if (action==6) {
		ClpPrimalColumnSteepest steep(10);
		lpSolver->setPrimalColumnPivotAlgorithm(steep);
	      }
	      break;
	    case SCALING:
	      lpSolver->scaling(action);
	      solver->setHintParam(OsiDoScale,action!=0,OsiHintTry);
	      doScaling = 1-action;
	      break;
	    case AUTOSCALE:
	      lpSolver->setAutomaticScaling(action!=0);
	      break;
	    case SPARSEFACTOR:
	      lpSolver->setSparseFactorization((1-action)!=0);
	      break;
	    case BIASLU:
	      lpSolver->factorization()->setBiasLU(action);
	      break;
	    case PERTURBATION:
	      if (action==0)
		lpSolver->setPerturbation(50);
	      else
		lpSolver->setPerturbation(100);
	      break;
	    case ERRORSALLOWED:
	      allowImportErrors = action;
	      break;
            case INTPRINT:
              printMode=action;
              break;
              //case ALGORITHM:
	      //algorithm  = action;
              //defaultSettings=false; // user knows what she is doing
              //abort();
	      //break;
	    case KEEPNAMES:
	      keepImportNames = 1-action;
	      break;
	    case PRESOLVE:
	      if (action==0)
		preSolve = 5;
	      else if (action==1)
		preSolve=0;
	      else if (action==2)
		preSolve=10;
	      else
		preSolveFile=true;
	      break;
	    case PFI:
	      lpSolver->factorization()->setForrestTomlin(action==0);
	      break;
	    case CRASH:
	      doCrash=action;
	      break;
	    case VECTOR:
	      doVector=action;
	      break;
	    case MESSAGES:
	      lpSolver->messageHandler()->setPrefix(action!=0);
	      break;
	    case CHOLESKY:
	      choleskyType = action;
	      break;
	    case GAMMA:
	      gamma=action;
	      break;
	    case BARRIERSCALE:
	      scaleBarrier=action;
	      break;
	    case KKT:
	      doKKT=action;
	      break;
	    case CROSSOVER:
	      crossover=action;
	      break;
	    case SOS:
	      doSOS=action;
	      break;
	    case GOMORYCUTS:
              defaultSettings=false; // user knows what she is doing
	      gomoryAction = action;
	      break;
	    case PROBINGCUTS:
              defaultSettings=false; // user knows what she is doing
	      probingAction = action;
	      break;
	    case KNAPSACKCUTS:
              defaultSettings=false; // user knows what she is doing
	      knapsackAction = action;
	      break;
	    case REDSPLITCUTS:
              defaultSettings=false; // user knows what she is doing
	      redsplitAction = action;
	      break;
	    case CLIQUECUTS:
              defaultSettings=false; // user knows what she is doing
	      cliqueAction = action;
	      break;
	    case FLOWCUTS:
              defaultSettings=false; // user knows what she is doing
	      flowAction = action;
	      break;
	    case MIXEDCUTS:
              defaultSettings=false; // user knows what she is doing
	      mixedAction = action;
	      break;
	    case TWOMIRCUTS:
              defaultSettings=false; // user knows what she is doing
	      twomirAction = action;
	      break;
	    case LANDPCUTS:
              defaultSettings=false; // user knows what she is doing
	      landpAction = action;
	      break;
	    case RESIDCUTS:
              defaultSettings=false; // user knows what she is doing
	      residualCapacityAction = action;
	      break;
	    case ROUNDING:
              defaultSettings=false; // user knows what she is doing
	      useRounding = action;
	      break;
	    case FPUMP:
              defaultSettings=false; // user knows what she is doing
              useFpump=action;
	      break;
	    case RINS:
              useRINS=action;
	      break;
            case CUTSSTRATEGY:
	      gomoryAction = action;
	      probingAction = action;
	      knapsackAction = action;
	      cliqueAction = action;
	      flowAction = action;
	      mixedAction = action;
	      twomirAction = action;
	      //landpAction = action;
              parameters[whichParam(GOMORYCUTS,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(PROBINGCUTS,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(KNAPSACKCUTS,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(CLIQUECUTS,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(FLOWCUTS,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(MIXEDCUTS,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(TWOMIRCUTS,numberParameters,parameters)].setCurrentOption(action);
              if (!action) {
                redsplitAction = action;
                parameters[whichParam(REDSPLITCUTS,numberParameters,parameters)].setCurrentOption(action);
                landpAction = action;
                parameters[whichParam(LANDPCUTS,numberParameters,parameters)].setCurrentOption(action);
                residualCapacityAction = action;
                parameters[whichParam(RESIDCUTS,numberParameters,parameters)].setCurrentOption(action);
              }
              break;
            case HEURISTICSTRATEGY:
              useRounding = action;
              useGreedy = action;
              useCombine = action;
              //useLocalTree = action;
              useFpump=action;
              parameters[whichParam(ROUNDING,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(GREEDY,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(COMBINE,numberParameters,parameters)].setCurrentOption(action);
              //parameters[whichParam(LOCALTREE,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(FPUMP,numberParameters,parameters)].setCurrentOption(action);
              break;
	    case GREEDY:
              defaultSettings=false; // user knows what she is doing
	      useGreedy = action;
	      break;
	    case COMBINE:
              defaultSettings=false; // user knows what she is doing
	      useCombine = action;
	      break;
	    case LOCALTREE:
              defaultSettings=false; // user knows what she is doing
	      useLocalTree = action;
	      break;
	    case COSTSTRATEGY:
	      useCosts=action;
	      break;
	    case NODESTRATEGY:
	      nodeStrategy=action;
	      break;
	    case PREPROCESS:
	      preProcess = action;
	      break;
	    case USESOLUTION:
	      useSolution = action;
	      break;
	    default:
	      abort();
	    }
	  }
	} else {
	  // action
	  if (type==EXIT) {
#ifdef COIN_HAS_ASL
            if(usingAmpl) {
              if (info.numberIntegers||info.numberBinary) {
                // integer
              } else {
                // linear
              }
              writeAmpl(&info);
              freeArrays2(&info);
              freeArgs(&info);
            }
#endif
	    break; // stop all
          }
	  switch (type) {
	  case DUALSIMPLEX:
	  case PRIMALSIMPLEX:
	  case SOLVECONTINUOUS:
	  case BARRIER:
	    if (goodModel) {
              double objScale = 
                parameters[whichParam(OBJSCALE2,numberParameters,parameters)].doubleValue();
              if (objScale!=1.0) {
                int iColumn;
                int numberColumns=lpSolver->numberColumns();
                double * dualColumnSolution = 
                  lpSolver->dualColumnSolution();
                ClpObjective * obj = lpSolver->objectiveAsObject();
                assert(dynamic_cast<ClpLinearObjective *> (obj));
                double offset;
                double * objective = obj->gradient(NULL,NULL,offset,true);
                for (iColumn=0;iColumn<numberColumns;iColumn++) {
                  dualColumnSolution[iColumn] *= objScale;
                  objective[iColumn] *= objScale;;
                }
                int iRow;
                int numberRows=lpSolver->numberRows();
                double * dualRowSolution = 
                  lpSolver->dualRowSolution();
                for (iRow=0;iRow<numberRows;iRow++) 
                  dualRowSolution[iRow] *= objScale;
                lpSolver->setObjectiveOffset(objScale*lpSolver->objectiveOffset());
              }
	      ClpSolve::SolveType method;
	      ClpSolve::PresolveType presolveType;
	      ClpSimplex * model2 = lpSolver;
              if (dualize) {
		//printf("dualize %d\n",dualize);
                model2 = ((ClpSimplexOther *) model2)->dualOfModel();
                sprintf(generalPrint,"Dual of model has %d rows and %d columns\n",
                       model2->numberRows(),model2->numberColumns());
		generalMessageHandler->message(CLP_GENERAL,generalMessages)
		  << generalPrint
		  <<CoinMessageEol;
                model2->setOptimizationDirection(1.0);
              }
              if (noPrinting)
                lpSolver->setLogLevel(0);
	      ClpSolve solveOptions;
              solveOptions.setPresolveActions(presolveOptions);
              solveOptions.setSubstitution(substitution);
	      if (preSolve!=5&&preSolve) {
		presolveType=ClpSolve::presolveNumber;
                if (preSolve<0) {
                  preSolve = - preSolve;
                  if (preSolve<=100) {
                    presolveType=ClpSolve::presolveNumber;
                    sprintf(generalPrint,"Doing %d presolve passes - picking up non-costed slacks\n",
                           preSolve);
		    generalMessageHandler->message(CLP_GENERAL,generalMessages)
		      << generalPrint
		      <<CoinMessageEol;
                    solveOptions.setDoSingletonColumn(true);
                  } else {
                    preSolve -=100;
                    presolveType=ClpSolve::presolveNumberCost;
                    sprintf(generalPrint,"Doing %d presolve passes - picking up costed slacks\n",
                           preSolve);
		    generalMessageHandler->message(CLP_GENERAL,generalMessages)
		      << generalPrint
		      <<CoinMessageEol;
                  }
                } 
	      } else if (preSolve) {
		presolveType=ClpSolve::presolveOn;
	      } else {
		presolveType=ClpSolve::presolveOff;
              }
	      solveOptions.setPresolveType(presolveType,preSolve);
	      if (type==DUALSIMPLEX||type==SOLVECONTINUOUS) {
		method=ClpSolve::useDual;
	      } else if (type==PRIMALSIMPLEX) {
		method=ClpSolve::usePrimalorSprint;
	      } else {
		method = ClpSolve::useBarrier;
		if (crossover==1) {
		  method=ClpSolve::useBarrierNoCross;
		} else if (crossover==2) {
		  ClpObjective * obj = lpSolver->objectiveAsObject();
		  if (obj->type()>1) {
		    method=ClpSolve::useBarrierNoCross;
		    presolveType=ClpSolve::presolveOff;
		    solveOptions.setPresolveType(presolveType,preSolve);
		  } 
		}
	      }
	      solveOptions.setSolveType(method);
	      if(preSolveFile)
		presolveOptions |= 0x40000000;
	      solveOptions.setSpecialOption(4,presolveOptions);
	      solveOptions.setSpecialOption(5,printOptions);
	      if (doVector) {
		ClpMatrixBase * matrix = lpSolver->clpMatrix();
		if (dynamic_cast< ClpPackedMatrix*>(matrix)) {
		  ClpPackedMatrix * clpMatrix = dynamic_cast< ClpPackedMatrix*>(matrix);
		  clpMatrix->makeSpecialColumnCopy();
		}
	      }
	      if (method==ClpSolve::useDual) {
		// dual
		if (doCrash)
		  solveOptions.setSpecialOption(0,1,doCrash); // crash
		else if (doIdiot)
		  solveOptions.setSpecialOption(0,2,doIdiot); // possible idiot
	      } else if (method==ClpSolve::usePrimalorSprint) {
		// primal
		// if slp turn everything off
		if (slpValue>0) {
		  doCrash=false;
		  doSprint=0;
		  doIdiot=-1;
		  solveOptions.setSpecialOption(1,10,slpValue); // slp
		  method=ClpSolve::usePrimal;
		}
		if (doCrash) {
		  solveOptions.setSpecialOption(1,1,doCrash); // crash
		} else if (doSprint>0) {
		  // sprint overrides idiot
		  solveOptions.setSpecialOption(1,3,doSprint); // sprint
		} else if (doIdiot>0) {
		  solveOptions.setSpecialOption(1,2,doIdiot); // idiot
		} else if (slpValue<=0) {
		  if (doIdiot==0) {
		    if (doSprint==0)
		      solveOptions.setSpecialOption(1,4); // all slack
		    else
		      solveOptions.setSpecialOption(1,9); // all slack or sprint
		  } else {
		    if (doSprint==0)
		      solveOptions.setSpecialOption(1,8); // all slack or idiot
		    else
		      solveOptions.setSpecialOption(1,7); // initiative
		  }
		}
		if (basisHasValues==-1)
		  solveOptions.setSpecialOption(1,11); // switch off values
	      } else if (method==ClpSolve::useBarrier||method==ClpSolve::useBarrierNoCross) {
		int barrierOptions = choleskyType;
		if (scaleBarrier)
		  barrierOptions |= 8;
		if (doKKT)
		  barrierOptions |= 16;
		if (gamma)
		  barrierOptions |= 32*gamma;
		if (crossover==3) 
		  barrierOptions |= 256; // try presolve in crossover
		solveOptions.setSpecialOption(4,barrierOptions);
	      }
	      model2->setMaximumSeconds(model.getMaximumSeconds());
#ifdef COIN_HAS_LINK
	      OsiSolverInterface * coinSolver = model.solver();
	      OsiSolverLink * linkSolver = dynamic_cast< OsiSolverLink*> (coinSolver);
	      if (!linkSolver) {
		model2->initialSolve(solveOptions);
	      } else {
		// special solver
		int testOsiOptions = parameters[whichParam(TESTOSI,numberParameters,parameters)].intValue();
		double * solution = NULL;
		if (testOsiOptions<10) {
		  solution = linkSolver->nonlinearSLP(slpValue>0 ? slpValue : 20 ,1.0e-5);
		} else if (testOsiOptions>=10) {
		  CoinModel coinModel = *linkSolver->coinModel();
		  ClpSimplex * tempModel = approximateSolution(coinModel,slpValue>0 ? slpValue : 50 ,1.0e-5,0);
		  assert (tempModel);
		  solution = CoinCopyOfArray(tempModel->primalColumnSolution(),coinModel.numberColumns());
		  delete tempModel;
		}
		if (solution) {
		  memcpy(model2->primalColumnSolution(),solution,
			 CoinMin(model2->numberColumns(),linkSolver->coinModel()->numberColumns())*sizeof(double));
		  delete [] solution;
		} else {
		  printf("No nonlinear solution\n");
		}
	      }
#else
	      model2->initialSolve(solveOptions);
#endif
	      {
		// map states
		/* clp status
		   -1 - unknown e.g. before solve or if postSolve says not optimal
		   0 - optimal
		   1 - primal infeasible
		   2 - dual infeasible
		   3 - stopped on iterations or time
		   4 - stopped due to errors
		   5 - stopped by event handler (virtual int ClpEventHandler::event()) */
		/* cbc status
		   -1 before branchAndBound
		   0 finished - check isProvenOptimal or isProvenInfeasible to see if solution found
		   (or check value of best solution)
		   1 stopped - on maxnodes, maxsols, maxtime
		   2 difficulties so run was abandoned
		   (5 event user programmed event occurred) */
		/* clp secondary status of problem - may get extended
		   0 - none
		   1 - primal infeasible because dual limit reached OR probably primal
		   infeasible but can't prove it (main status 4)
		   2 - scaled problem optimal - unscaled problem has primal infeasibilities
		   3 - scaled problem optimal - unscaled problem has dual infeasibilities
		   4 - scaled problem optimal - unscaled problem has primal and dual infeasibilities
		   5 - giving up in primal with flagged variables
		   6 - failed due to empty problem check 
		   7 - postSolve says not optimal
		   8 - failed due to bad element check
		   9 - status was 3 and stopped on time
		   100 up - translation of enum from ClpEventHandler
		*/
		/* cbc secondary status of problem
		   -1 unset (status_ will also be -1)
		   0 search completed with solution
		   1 linear relaxation not feasible (or worse than cutoff)
		   2 stopped on gap
		   3 stopped on nodes
		   4 stopped on time
		   5 stopped on user event
		   6 stopped on solutions
		   7 linear relaxation unbounded
		*/
		int iStatus = model2->status();
		int iStatus2 = model2->secondaryStatus();
		if (iStatus==0) {
		  iStatus2=0;
		} else if (iStatus==1) {
		  iStatus=0;
		  iStatus2=1; // say infeasible
		} else if (iStatus==2) {
		  iStatus=0;
		  iStatus2=7; // say unbounded
		} else if (iStatus==3) {
		  iStatus=1;
		  if (iStatus2==9)
		    iStatus2=4;
		  else
		    iStatus2=3; // Use nodes - as closer than solutions
		} else if (iStatus==4) {
		  iStatus=2; // difficulties
		  iStatus2=0; 
		}
		model.setProblemStatus(iStatus);
		model.setSecondaryStatus(iStatus2);
		assert (lpSolver==clpSolver->getModelPtr());
		assert (clpSolver==model.solver());
		clpSolver->setWarmStart(NULL);
		// and in babModel if exists
		if (babModel) {
		  babModel->setProblemStatus(iStatus);
		  babModel->setSecondaryStatus(iStatus2);
		}
	      }
	      basisHasValues=1;
              if (dualize) {
                int returnCode=((ClpSimplexOther *) lpSolver)->restoreFromDual(model2);
                delete model2;
		if (returnCode&&dualize!=2)
		  lpSolver->primal(1);
                model2=lpSolver;
              }
#ifdef COIN_HAS_ASL
              if (usingAmpl) {
                double value = model2->getObjValue()*model2->getObjSense();
                char buf[300];
                int pos=0;
                int iStat = model2->status();
                if (iStat==0) {
                  pos += sprintf(buf+pos,"optimal," );
                } else if (iStat==1) {
                  // infeasible
                  pos += sprintf(buf+pos,"infeasible,");
                } else if (iStat==2) {
                  // unbounded
                  pos += sprintf(buf+pos,"unbounded,");
                } else if (iStat==3) {
                  pos += sprintf(buf+pos,"stopped on iterations or time,");
                } else if (iStat==4) {
                  iStat = 7;
                  pos += sprintf(buf+pos,"stopped on difficulties,");
                } else if (iStat==5) {
                  iStat = 3;
                  pos += sprintf(buf+pos,"stopped on ctrl-c,");
                } else {
                  pos += sprintf(buf+pos,"status unknown,");
                  iStat=6;
                }
                info.problemStatus=iStat;
                info.objValue = value;
                pos += sprintf(buf+pos," objective %.*g",ampl_obj_prec(),
                               value);
                sprintf(buf+pos,"\n%d iterations",
                        model2->getIterationCount());
                free(info.primalSolution);
                int numberColumns=model2->numberColumns();
                info.primalSolution = (double *) malloc(numberColumns*sizeof(double));
                CoinCopyN(model2->primalColumnSolution(),numberColumns,info.primalSolution);
                int numberRows = model2->numberRows();
                free(info.dualSolution);
                info.dualSolution = (double *) malloc(numberRows*sizeof(double));
                CoinCopyN(model2->dualRowSolution(),numberRows,info.dualSolution);
                CoinWarmStartBasis * basis = model2->getBasis();
                free(info.rowStatus);
                info.rowStatus = (int *) malloc(numberRows*sizeof(int));
                free(info.columnStatus);
                info.columnStatus = (int *) malloc(numberColumns*sizeof(int));
                // Put basis in 
                int i;
                // free,basic,ub,lb are 0,1,2,3
                for (i=0;i<numberRows;i++) {
                  CoinWarmStartBasis::Status status = basis->getArtifStatus(i);
                  info.rowStatus[i]=status;
                }
                for (i=0;i<numberColumns;i++) {
                  CoinWarmStartBasis::Status status = basis->getStructStatus(i);
                  info.columnStatus[i]=status;
                }
                // put buffer into info
                strcpy(info.buffer,buf);
                delete basis;
              }
#endif
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
          case STATISTICS:
	    if (goodModel) {
              // If presolve on look at presolved
              bool deleteModel2=false;
              ClpSimplex * model2 = lpSolver;
              if (preSolve) {
                ClpPresolve pinfo;
                int presolveOptions2 = presolveOptions&~0x40000000;
                if ((presolveOptions2&0xffff)!=0)
                  pinfo.setPresolveActions(presolveOptions2);
                pinfo.setSubstitution(substitution);
                if ((printOptions&1)!=0)
                  pinfo.statistics();
                double presolveTolerance = 
                  parameters[whichParam(PRESOLVETOLERANCE,numberParameters,parameters)].doubleValue();
                model2 = 
                  pinfo.presolvedModel(*lpSolver,presolveTolerance,
                                       true,preSolve);
                if (model2) {
                  printf("Statistics for presolved model\n");
                  deleteModel2=true;
                } else {
                  printf("Presolved model looks infeasible - will use unpresolved\n");
                  model2 = lpSolver;
                }
              } else {
                printf("Statistics for unpresolved model\n");
                model2 =  lpSolver;
              }
              statistics(lpSolver,model2);
              if (deleteModel2)
                delete model2;
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case TIGHTEN:
	    if (goodModel) {
     	      int numberInfeasibilities = lpSolver->tightenPrimalBounds();
	      if (numberInfeasibilities)
		std::cout<<"** Analysis indicates model infeasible"<<std::endl;
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case PLUSMINUS:
	    if (goodModel) {
	      ClpMatrixBase * saveMatrix = lpSolver->clpMatrix();
	      ClpPackedMatrix* clpMatrix =
		dynamic_cast< ClpPackedMatrix*>(saveMatrix);
	      if (clpMatrix) {
		ClpPlusMinusOneMatrix * newMatrix = new ClpPlusMinusOneMatrix(*(clpMatrix->matrix()));
		if (newMatrix->getIndices()) {
		  lpSolver->replaceMatrix(newMatrix);
		  delete saveMatrix;
		  std::cout<<"Matrix converted to +- one matrix"<<std::endl;
		} else {
		  std::cout<<"Matrix can not be converted to +- 1 matrix"<<std::endl;
		}
	      } else {
		std::cout<<"Matrix not a ClpPackedMatrix"<<std::endl;
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case OUTDUPROWS:
	    if (goodModel) {
	      int numberRows = clpSolver->getNumRows();
              //int nOut = outDupRow(clpSolver);
	      CglDuplicateRow dupcuts(clpSolver);
	      storedCuts = dupcuts.outDuplicates(clpSolver)!=0;
	      int nOut = numberRows-clpSolver->getNumRows();
              if (nOut&&!noPrinting)
                sprintf(generalPrint,"%d rows eliminated\n",nOut);
	      generalMessageHandler->message(CLP_GENERAL,generalMessages)
		<< generalPrint
		<<CoinMessageEol;
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case NETWORK:
	    if (goodModel) {
	      ClpMatrixBase * saveMatrix = lpSolver->clpMatrix();
	      ClpPackedMatrix* clpMatrix =
		dynamic_cast< ClpPackedMatrix*>(saveMatrix);
	      if (clpMatrix) {
		ClpNetworkMatrix * newMatrix = new ClpNetworkMatrix(*(clpMatrix->matrix()));
		if (newMatrix->getIndices()) {
		  lpSolver->replaceMatrix(newMatrix);
		  delete saveMatrix;
		  std::cout<<"Matrix converted to network matrix"<<std::endl;
		} else {
		  std::cout<<"Matrix can not be converted to network matrix"<<std::endl;
		}
	      } else {
		std::cout<<"Matrix not a ClpPackedMatrix"<<std::endl;
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
          case MIPLIB:
            // User can set options - main difference is lack of model and CglPreProcess
            goodModel=true;
/*
  Run branch-and-cut. First set a few options -- node comparison, scaling. If
  the solver is Clp, consider running some presolve code (not yet converted
  this to generic OSI) with branch-and-cut. If presolve is disabled, or the
  solver is not Clp, simply run branch-and-cut. Print elapsed time at the end.
*/
	  case BAB: // branchAndBound
          case STRENGTHEN:
            if (goodModel) {
              bool miplib = type==MIPLIB;
              int logLevel = parameters[slog].intValue();
              // Reduce printout
              if (logLevel<=1)
                model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
              {
                OsiSolverInterface * solver = model.solver();
                OsiClpSolverInterface * si =
                  dynamic_cast<OsiClpSolverInterface *>(solver) ;
                assert (si != NULL);
		// See if quadratic
#ifdef COIN_HAS_LINK
		if (!complicatedInteger) {
		  ClpSimplex * lpSolver = si->getModelPtr();
		  ClpQuadraticObjective * obj = (dynamic_cast< ClpQuadraticObjective*>(lpSolver->objectiveAsObject()));
		  if (obj) {
		    preProcess=0;
		    int testOsiOptions = parameters[whichParam(TESTOSI,numberParameters,parameters)].intValue();
		    parameters[whichParam(TESTOSI,numberParameters,parameters)].setIntValue(CoinMax(0,testOsiOptions));
		    // create coin model
		    coinModel = lpSolver->createCoinModel();
		    assert (coinModel);
		    // load from coin model
		    OsiSolverLink solver1;
		    OsiSolverInterface * solver2 = solver1.clone();
		    model.assignSolver(solver2,false);
		    OsiSolverLink * si =
		      dynamic_cast<OsiSolverLink *>(model.solver()) ;
		    assert (si != NULL);
		    si->setDefaultMeshSize(0.001);
		    // need some relative granularity
		    si->setDefaultBound(100.0);
		    si->setDefaultMeshSize(0.001);
		    si->setDefaultBound(1000.0);
		    si->setIntegerPriority(1000);
		    si->setBiLinearPriority(10000);
		    si->setSpecialOptions2(2+4+8);
		    CoinModel * model2 = (CoinModel *) coinModel;
		    si->load(*model2,true, parameters[log].intValue());
		    // redo
		    solver = model.solver();
		    clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
		    lpSolver = clpSolver->getModelPtr();
		    clpSolver->messageHandler()->setLogLevel(0) ;
		    testOsiParameters=0;
		    complicatedInteger=2;  // allow cuts
		    OsiSolverInterface * coinSolver = model.solver();
		    OsiSolverLink * linkSolver = dynamic_cast< OsiSolverLink*> (coinSolver);
		    if (linkSolver->quadraticModel()) {
		      ClpSimplex * qp = linkSolver->quadraticModel();
		      //linkSolver->nonlinearSLP(CoinMax(slpValue,10),1.0e-5);
		      qp->nonlinearSLP(CoinMax(slpValue,40),1.0e-5);
		      qp->primal(1);
		      OsiSolverLinearizedQuadratic solver2(qp);
		      const double * solution=NULL;
		      // Reduce printout
		      solver2.setHintParam(OsiDoReducePrint,true,OsiHintTry);
		      CbcModel model2(solver2);
		      // Now do requested saves and modifications
		      CbcModel * cbcModel = & model2;
		      OsiSolverInterface * osiModel = model2.solver();
		      OsiClpSolverInterface * osiclpModel = dynamic_cast< OsiClpSolverInterface*> (osiModel);
		      ClpSimplex * clpModel = osiclpModel->getModelPtr();
		      
		      // Set changed values
		      
		      CglProbing probing;
		      probing.setMaxProbe(10);
		      probing.setMaxLook(10);
		      probing.setMaxElements(200);
		      probing.setMaxProbeRoot(50);
		      probing.setMaxLookRoot(10);
		      probing.setRowCuts(3);
		      probing.setUsingObjective(true);
		      cbcModel->addCutGenerator(&probing,-1,"Probing",true,false,false,-100,-1,-1);
		      cbcModel->cutGenerator(0)->setTiming(true);
		      
		      CglGomory gomory;
		      gomory.setLimitAtRoot(512);
		      cbcModel->addCutGenerator(&gomory,-98,"Gomory",true,false,false,-100,-1,-1);
		      cbcModel->cutGenerator(1)->setTiming(true);
		      
		      CglKnapsackCover knapsackCover;
		      cbcModel->addCutGenerator(&knapsackCover,-98,"KnapsackCover",true,false,false,-100,-1,-1);
		      cbcModel->cutGenerator(2)->setTiming(true);
		      
		      CglRedSplit redSplit;
		      cbcModel->addCutGenerator(&redSplit,-99,"RedSplit",true,false,false,-100,-1,-1);
		      cbcModel->cutGenerator(3)->setTiming(true);
		      
		      CglClique clique;
		      clique.setStarCliqueReport(false);
		      clique.setRowCliqueReport(false);
		      clique.setMinViolation(0.1);
		      cbcModel->addCutGenerator(&clique,-98,"Clique",true,false,false,-100,-1,-1);
		      cbcModel->cutGenerator(4)->setTiming(true);
		      
		      CglMixedIntegerRounding2 mixedIntegerRounding2;
		      cbcModel->addCutGenerator(&mixedIntegerRounding2,-98,"MixedIntegerRounding2",true,false,false,-100,-1,-1);
		      cbcModel->cutGenerator(5)->setTiming(true);
		      
		      CglFlowCover flowCover;
		      cbcModel->addCutGenerator(&flowCover,-98,"FlowCover",true,false,false,-100,-1,-1);
		      cbcModel->cutGenerator(6)->setTiming(true);
		      
		      CglTwomir twomir;
		      twomir.setMaxElements(250);
		      cbcModel->addCutGenerator(&twomir,-99,"Twomir",true,false,false,-100,-1,-1);
		      cbcModel->cutGenerator(7)->setTiming(true);
		      
		      CbcHeuristicFPump heuristicFPump(*cbcModel);
		      heuristicFPump.setWhen(13);
		      heuristicFPump.setMaximumPasses(20);
		      heuristicFPump.setMaximumRetries(7);
		      heuristicFPump.setHeuristicName("feasibility pump");
		      heuristicFPump.setInitialWeight(1);
		      heuristicFPump.setFractionSmall(0.6);
		      cbcModel->addHeuristic(&heuristicFPump);
		      
		      CbcRounding rounding(*cbcModel);
		      rounding.setHeuristicName("rounding");
		      cbcModel->addHeuristic(&rounding);
		      
		      CbcHeuristicLocal heuristicLocal(*cbcModel);
		      heuristicLocal.setHeuristicName("combine solutions");
		      heuristicLocal.setSearchType(1);
		      heuristicLocal.setFractionSmall(0.6);
		      cbcModel->addHeuristic(&heuristicLocal);
		      
		      CbcHeuristicGreedyCover heuristicGreedyCover(*cbcModel);
		      heuristicGreedyCover.setHeuristicName("greedy cover");
		      cbcModel->addHeuristic(&heuristicGreedyCover);
		      
		      CbcHeuristicGreedyEquality heuristicGreedyEquality(*cbcModel);
		      heuristicGreedyEquality.setHeuristicName("greedy equality");
		      cbcModel->addHeuristic(&heuristicGreedyEquality);
		      
		      CbcCompareDefault compare;
		      cbcModel->setNodeComparison(compare);
		      cbcModel->setNumberBeforeTrust(5);
		      cbcModel->setSpecialOptions(2);
		      cbcModel->messageHandler()->setLogLevel(1);
		      cbcModel->setMaximumCutPassesAtRoot(-100);
		      cbcModel->setMaximumCutPasses(1);
		      cbcModel->setMinimumDrop(0.05);
		      // For branchAndBound this may help
		      clpModel->defaultFactorizationFrequency();
		      clpModel->setDualBound(1.0001e+08);
		      clpModel->setPerturbation(50);
		      osiclpModel->setSpecialOptions(193);
		      osiclpModel->messageHandler()->setLogLevel(0);
		      osiclpModel->setIntParam(OsiMaxNumIterationHotStart,100);
		      osiclpModel->setHintParam(OsiDoReducePrint,true,OsiHintTry);
		      // You can save some time by switching off message building
		      // clpModel->messagesPointer()->setDetailMessages(100,10000,(int *) NULL);
		      
		      // Solve
		      
		      cbcModel->initialSolve();
		      if (clpModel->tightenPrimalBounds()!=0) {
			std::cout<<"Problem is infeasible - tightenPrimalBounds!"<<std::endl;
			exit(1);
		      }
		      clpModel->dual();  // clean up
		      cbcModel->initialSolve();
#ifdef CBC_THREAD
		      int numberThreads =parameters[whichParam(THREADS,numberParameters,parameters)].intValue();
		      cbcModel->setNumberThreads(numberThreads%100);
		      cbcModel->setThreadMode(numberThreads/100);
#endif
		      cbcModel->branchAndBound();
		      OsiSolverLinearizedQuadratic * solver3 = dynamic_cast<OsiSolverLinearizedQuadratic *> (model2.solver());
		      assert (solver3);
		      solution = solver3->bestSolution();
		      double bestObjectiveValue = solver3->bestObjectiveValue();
		      linkSolver->setBestObjectiveValue(bestObjectiveValue);
		      linkSolver->setBestSolution(solution,solver3->getNumCols());
		      CbcHeuristicDynamic3 dynamic(model);
		      dynamic.setHeuristicName("dynamic pass thru");
		      model.addHeuristic(&dynamic);
		      // if convex
		      if ((linkSolver->specialOptions2()&4)!=0) {
			int numberColumns = coinModel->numberColumns();
			assert (linkSolver->objectiveVariable()==numberColumns);
			// add OA cut
			double offset;
			double * gradient = new double [numberColumns+1];
			memcpy(gradient,qp->objectiveAsObject()->gradient(qp,solution,offset,true,2),
			       numberColumns*sizeof(double));
			double rhs = 0.0;
			int * column = new int[numberColumns+1];
			int n=0;
			for (int i=0;i<numberColumns;i++) {
			  double value = gradient[i];
			  if (fabs(value)>1.0e-12) {
			    gradient[n]=value;
			    rhs += value*solution[i];
			    column[n++]=i;
			  }
			}
			gradient[n]=-1.0;
			column[n++]=numberColumns;
			storedAmpl.addCut(-COIN_DBL_MAX,offset+1.0e-7,n,column,gradient);
			delete [] gradient;
			delete [] column;
		      }
		      // could do three way branching round a) continuous b) best solution
		      printf("obj %g\n",bestObjectiveValue);
		      linkSolver->initialSolve();
		    }
		  }
		}
#endif
                si->setSpecialOptions(0x40000000);
              }
              if (!miplib) {
		if (!preSolve) {
		  model.solver()->setHintParam(OsiDoPresolveInInitial,false,OsiHintTry);
		  model.solver()->setHintParam(OsiDoPresolveInResolve,false,OsiHintTry);
		}
		double time1a = CoinCpuTime();
                model.initialSolve();
                OsiSolverInterface * solver = model.solver();
                OsiClpSolverInterface * si =
                  dynamic_cast<OsiClpSolverInterface *>(solver) ;
		ClpSimplex * clpSolver = si->getModelPtr();
		clpSolver->setSpecialOptions(clpSolver->specialOptions()|0x01000000); // say is Cbc (and in branch and bound)
		if (!noPrinting) {
		  sprintf(generalPrint,"Continuous objective value is %g - %.2f seconds",
			  solver->getObjValue(),CoinCpuTime()-time1a);
		  generalMessageHandler->message(CLP_GENERAL,generalMessages)
		    << generalPrint
		    <<CoinMessageEol;
		}
		if (!complicatedInteger&&clpSolver->tightenPrimalBounds()!=0) {
		  std::cout<<"Problem is infeasible - tightenPrimalBounds!"<<std::endl;
		  exit(1);
		}
		if (clpSolver->dualBound()==1.0e10) {
		  // user did not set - so modify
		  // get largest scaled away from bound
		  double largest=1.0e-12;
		  int numberRows = clpSolver->numberRows();
		  const double * rowPrimal = clpSolver->primalRowSolution();
		  const double * rowLower = clpSolver->rowLower();
		  const double * rowUpper = clpSolver->rowUpper();
		  const double * rowScale = clpSolver->rowScale();
		  int iRow;
		  for (iRow=0;iRow<numberRows;iRow++) {
		    double value = rowPrimal[iRow];
		    double above = value-rowLower[iRow];
		    double below = rowUpper[iRow]-value;
		    if (rowScale) {
		      double multiplier = rowScale[iRow];
		      above *= multiplier;
		      below *= multiplier;
		    }
		    if (above<1.0e12)
		      largest = CoinMax(largest,above);
		    if (below<1.0e12)
		      largest = CoinMax(largest,below);
		  }
		  
		  int numberColumns = clpSolver->numberColumns();
		  const double * columnPrimal = clpSolver->primalColumnSolution();
		  const double * columnLower = clpSolver->columnLower();
		  const double * columnUpper = clpSolver->columnUpper();
		  const double * columnScale = clpSolver->columnScale();
		  int iColumn;
		  for (iColumn=0;iColumn<numberColumns;iColumn++) {
		    double value = columnPrimal[iColumn];
		    double above = value-columnLower[iColumn];
		    double below = columnUpper[iColumn]-value;
		    if (columnScale) {
		      double multiplier = 1.0/columnScale[iColumn];
		      above *= multiplier;
		      below *= multiplier;
		    }
		    if (above<1.0e12)
		      largest = CoinMax(largest,above);
		    if (below<1.0e12)
		      largest = CoinMax(largest,below);
		  }
		  //if (!noPrinting)
		  //std::cout<<"Largest (scaled) away from bound "<<largest<<std::endl;
		  clpSolver->setDualBound(CoinMax(1.0001e8,CoinMin(1000.0*largest,1.00001e10)));
		}
		si->resolve();  // clean up
	      }
              // If user made settings then use them
              if (!defaultSettings) {
                OsiSolverInterface * solver = model.solver();
                if (!doScaling)
                  solver->setHintParam(OsiDoScale,false,OsiHintTry);
                OsiClpSolverInterface * si =
                  dynamic_cast<OsiClpSolverInterface *>(solver) ;
                assert (si != NULL);
                // get clp itself
                ClpSimplex * modelC = si->getModelPtr();
                //if (modelC->tightenPrimalBounds()!=0) {
                //std::cout<<"Problem is infeasible!"<<std::endl;
                //break;
                //}
                // bounds based on continuous
                if (tightenFactor&&!complicatedInteger) {
                  if (modelC->tightenPrimalBounds(tightenFactor)!=0) {
                    std::cout<<"Problem is infeasible!"<<std::endl;
                    break;
                  }
                }
                if (djFix<1.0e20) {
                  // do some fixing
                  int numberColumns = modelC->numberColumns();
                  int i;
                  const char * type = modelC->integerInformation();
                  double * lower = modelC->columnLower();
                  double * upper = modelC->columnUpper();
                  double * solution = modelC->primalColumnSolution();
                  double * dj = modelC->dualColumnSolution();
                  int numberFixed=0;
                  for (i=0;i<numberColumns;i++) {
                    if (type[i]) {
                      double value = solution[i];
                      if (value<lower[i]+1.0e-5&&dj[i]>djFix) {
                        solution[i]=lower[i];
                        upper[i]=lower[i];
                        numberFixed++;
                      } else if (value>upper[i]-1.0e-5&&dj[i]<-djFix) {
                        solution[i]=upper[i];
                        lower[i]=upper[i];
                        numberFixed++;
                      }
                    }
                  }
                  sprintf(generalPrint,"%d columns fixed\n",numberFixed);
		  generalMessageHandler->message(CLP_GENERAL,generalMessages)
		    << generalPrint
		    <<CoinMessageEol;
                }
              }
              // See if we want preprocessing
              OsiSolverInterface * saveSolver=NULL;
              CglPreProcess process;
	      // Say integers in sync 
	      bool integersOK=true;
              delete babModel;
              babModel = new CbcModel(model);
              OsiSolverInterface * solver3 = clpSolver->clone();
              babModel->assignSolver(solver3);
              OsiClpSolverInterface * clpSolver2 = dynamic_cast< OsiClpSolverInterface*> (babModel->solver());
              int numberChanged=0;
              if (clpSolver2->messageHandler()->logLevel())
                clpSolver2->messageHandler()->setLogLevel(1);
              if (logLevel>-1)
                clpSolver2->messageHandler()->setLogLevel(logLevel);
              lpSolver = clpSolver2->getModelPtr();
              if (lpSolver->factorizationFrequency()==200&&!miplib) {
                // User did not touch preset
                int numberRows = lpSolver->numberRows();
                const int cutoff1=10000;
                const int cutoff2=100000;
                const int base=75;
                const int freq0 = 50;
                const int freq1=200;
                const int freq2=400;
                const int maximum=1000;
                int frequency;
                if (numberRows<cutoff1)
                  frequency=base+numberRows/freq0;
                else if (numberRows<cutoff2)
                  frequency=base+cutoff1/freq0 + (numberRows-cutoff1)/freq1;
                else
                  frequency=base+cutoff1/freq0 + (cutoff2-cutoff1)/freq1 + (numberRows-cutoff2)/freq2;
                lpSolver->setFactorizationFrequency(CoinMin(maximum,frequency));
              }
              time2 = CoinCpuTime();
              totalTime += time2-time1;
              time1 = time2;
              double timeLeft = babModel->getMaximumSeconds();
              int numberOriginalColumns = babModel->solver()->getNumCols();
              if (preProcess==7) {
		// use strategy instead
		preProcess=0;
		useStrategy=true;
#ifdef COIN_HAS_LINK
		// empty out any cuts
		if (storedAmpl.sizeRowCuts()) {
		  printf("Emptying ampl stored cuts as internal preprocessing\n");
		  CglStored temp;
		  storedAmpl=temp;
		}
#endif
	      }
              if (preProcess&&type==BAB) {
		// See if sos from mps file
		if (numberSOS==0&&clpSolver->numberSOS()&&doSOS) {
		  // SOS
		  numberSOS = clpSolver->numberSOS();
		  const CoinSet * setInfo = clpSolver->setInfo();
		  sosStart = new int [numberSOS+1];
		  sosType = new char [numberSOS];
		  int i;
		  int nTotal=0;
		  sosStart[0]=0;
		  for ( i=0;i<numberSOS;i++) {
		    int type = setInfo[i].setType();
		    int n=setInfo[i].numberEntries();
		    sosType[i]=type;
		    nTotal += n;
		    sosStart[i+1] = nTotal;
		  }
		  sosIndices = new int[nTotal];
		  sosReference = new double [nTotal];
		  for (i=0;i<numberSOS;i++) {
		    int n=setInfo[i].numberEntries();
		    const int * which = setInfo[i].which();
		    const double * weights = setInfo[i].weights();
		    int base = sosStart[i];
		    for (int j=0;j<n;j++) {
		      int k=which[j];
		      sosIndices[j+base]=k;
		      sosReference[j+base] = weights ? weights[j] : (double) j;
		    }
                  }
		}
                saveSolver=babModel->solver()->clone();
                /* Do not try and produce equality cliques and
                   do up to 10 passes */
                OsiSolverInterface * solver2;
                {
                  // Tell solver we are in Branch and Cut
                  saveSolver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo) ;
                  // Default set of cut generators
                  CglProbing generator1;
                  generator1.setUsingObjective(1);
                  generator1.setMaxPass(3);
                  generator1.setMaxProbeRoot(saveSolver->getNumCols());
                  generator1.setMaxElements(100);
                  generator1.setMaxLookRoot(50);
                  generator1.setRowCuts(3);
                  // Add in generators
                  process.addCutGenerator(&generator1);
                  int translate[]={9999,0,0,-1,2,3,-2};
		  process.passInMessageHandler(babModel->messageHandler());
                  //process.messageHandler()->setLogLevel(babModel->logLevel());
#ifdef COIN_HAS_ASL
                  if (info.numberSos&&doSOS&&usingAmpl) {
                    // SOS
                    numberSOS = info.numberSos;
                    sosStart = info.sosStart;
                    sosIndices = info.sosIndices;
                  }
#endif
                  if (numberSOS&&doSOS) {
                    // SOS
                    int numberColumns = saveSolver->getNumCols();
                    char * prohibited = new char[numberColumns];
                    memset(prohibited,0,numberColumns);
                    int n=sosStart[numberSOS];
                    for (int i=0;i<n;i++) {
                      int iColumn = sosIndices[i];
                      prohibited[iColumn]=1;
                    }
                    process.passInProhibited(prohibited,numberColumns);
                    delete [] prohibited;
                  }
		  if (model.numberObjects()) {
		    OsiObject ** oldObjects = babModel->objects();
		    int numberOldObjects = babModel->numberObjects();
                    // SOS
                    int numberColumns = saveSolver->getNumCols();
                    char * prohibited = new char[numberColumns];
                    memset(prohibited,0,numberColumns);
		    for (int iObj = 0;iObj<numberOldObjects;iObj++) {
		      CbcSOS * obj =
			dynamic_cast <CbcSOS *>(oldObjects[iObj]) ;
		      if (obj) {
			int n=obj->numberMembers();
			const int * which = obj->members();
			for (int i=0;i<n;i++) {
			  int iColumn = which[i];
			  prohibited[iColumn]=1;
			}
		      }
                    }
                    process.passInProhibited(prohibited,numberColumns);
                    delete [] prohibited;
		  }
		  int numberPasses = 10;
		  if (tunePreProcess>=1000000) {
		    numberPasses = (tunePreProcess/1000000)-1;
		    tunePreProcess = tunePreProcess % 1000000;
		  } else if (tunePreProcess>=1000) {
		    numberPasses = (tunePreProcess/1000)-1;
		    tunePreProcess = tunePreProcess % 1000;
		  }
		  if (doSprint>0) {
		    // Sprint for primal solves
		    ClpSolve::SolveType method = ClpSolve::usePrimalorSprint;
		    ClpSolve::PresolveType presolveType = ClpSolve::presolveOff;
		    int numberPasses = 5;
		    int options[] = {0,3,0,0,0,0};
		    int extraInfo[] = {-1,20,-1,-1,-1,-1};
		    extraInfo[1]=doSprint;
		    int independentOptions[] = {0,0,3};
		    ClpSolve clpSolve(method,presolveType,numberPasses,
				      options,extraInfo,independentOptions);
		    // say use in OsiClp
		    clpSolve.setSpecialOption(6,1);
		    OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (saveSolver);
		    osiclp->setSolveOptions(clpSolve);
		    osiclp->setHintParam(OsiDoDualInResolve,false);
		    // switch off row copy
		    osiclp->getModelPtr()->setSpecialOptions(osiclp->getModelPtr()->specialOptions()|256);
		    osiclp->getModelPtr()->setInfeasibilityCost(1.0e11);
		  }
                  solver2 = process.preProcessNonDefault(*saveSolver,translate[preProcess],numberPasses,
							 tunePreProcess);
		  integersOK=false; // We need to redo if CbcObjects exist
                  // Tell solver we are not in Branch and Cut
                  saveSolver->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo) ;
                  if (solver2)
                    solver2->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo) ;
                }
#ifdef COIN_HAS_ASL
                if (!solver2&&usingAmpl) {
                  // infeasible
                  info.problemStatus=1;
                  info.objValue = 1.0e100;
                  sprintf(info.buffer,"infeasible/unbounded by pre-processing");
                  info.primalSolution=NULL;
                  info.dualSolution=NULL;
                  break;
                }
#endif
                if (!noPrinting) {
                  if (!solver2) {
                    sprintf(generalPrint,"Pre-processing says infeasible or unbounded\n");
		    generalMessageHandler->message(CLP_GENERAL,generalMessages)
		      << generalPrint
		      <<CoinMessageEol;
                    break;
                  } else {
                    //printf("processed model has %d rows, %d columns and %d elements\n",
		    //     solver2->getNumRows(),solver2->getNumCols(),solver2->getNumElements());
                  }
                }
                //solver2->resolve();
                if (preProcess==2) {
                  OsiClpSolverInterface * clpSolver2 = dynamic_cast< OsiClpSolverInterface*> (solver2);
                  ClpSimplex * lpSolver = clpSolver2->getModelPtr();
                  lpSolver->writeMps("presolved.mps",0,1,lpSolver->optimizationDirection());
                  printf("Preprocessed model (minimization) on presolved.mps\n");
                }
                // we have to keep solver2 so pass clone
                solver2 = solver2->clone();
                babModel->assignSolver(solver2);
                babModel->initialSolve();
                babModel->setMaximumSeconds(timeLeft-(CoinCpuTime()-time1));
              }
              // now tighten bounds
              if (!miplib) {
                OsiClpSolverInterface * si =
                  dynamic_cast<OsiClpSolverInterface *>(babModel->solver()) ;
                assert (si != NULL);
                // get clp itself
                ClpSimplex * modelC = si->getModelPtr();
                //if (noPrinting)
		//modelC->setLogLevel(0);
                if (!complicatedInteger&&modelC->tightenPrimalBounds()!=0) {
                  std::cout<<"Problem is infeasible!"<<std::endl;
                  break;
                }
                si->resolve();
              }
#if 0
	      numberDebugValues=599;
	      debugValues = new double[numberDebugValues];
	      CoinZeroN(debugValues,numberDebugValues);
	      debugValues[3]=1.0;
	      debugValues[6]=25.0;
	      debugValues[9]=4.0;
	      debugValues[26]=4.0;
	      debugValues[27]=6.0;
	      debugValues[35]=8.0;
	      debugValues[53]=21.0;
	      debugValues[56]=4.0;
#endif
              if (debugValues) {
                // for debug
                std::string problemName ;
                babModel->solver()->getStrParam(OsiProbName,problemName) ;
                babModel->solver()->activateRowCutDebugger(problemName.c_str()) ;
                twomirGen.probname_=strdup(problemName.c_str());
                // checking seems odd
                //redsplitGen.set_given_optsol(babModel->solver()->getRowCutDebuggerAlways()->optimalSolution(),
                //                         babModel->getNumCols());
              }
	      int testOsiOptions = parameters[whichParam(TESTOSI,numberParameters,parameters)].intValue();
#ifdef COIN_HAS_ASL
	      // If linked then see if expansion wanted
	      {
		OsiSolverLink * solver3 = dynamic_cast<OsiSolverLink *> (babModel->solver());
		int options = parameters[whichParam(MIPOPTIONS,numberParameters,parameters)].intValue()/10000;
		if (solver3||(options&16)!=0) {
		  if (options) {
		    /*
		      1 - force mini branch and bound
		      2 - set priorities high on continuous
		      4 - try adding OA cuts
		      8 - try doing quadratic linearization
		      16 - try expanding knapsacks
		    */
		    if ((options&16)) {
		      int numberColumns = saveCoinModel.numberColumns();
		      int numberRows = saveCoinModel.numberRows();
		      whichColumn = new int[numberColumns];
		      knapsackStart=new int[numberRows+1];
		      knapsackRow=new int[numberRows];
		      numberKnapsack=10000;
		      int extra1 = parameters[whichParam(EXTRA1,numberParameters,parameters)].intValue();
		      int extra2 = parameters[whichParam(EXTRA2,numberParameters,parameters)].intValue();
		      int logLevel = parameters[log].intValue();
		      OsiSolverInterface * solver = expandKnapsack(saveCoinModel,whichColumn,knapsackStart,
								   knapsackRow,numberKnapsack,
								   storedAmpl,logLevel,extra1,extra2,
								   saveTightenedModel);
		      if (solver) {
			clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
			assert (clpSolver);
			lpSolver = clpSolver->getModelPtr();
			babModel->assignSolver(solver);
			testOsiOptions=0;
			// allow gomory
			complicatedInteger=0;
			// Priorities already done
			free(info.priorities);
			info.priorities=NULL;
		      } else {
			numberKnapsack=0;
			delete [] whichColumn;
			delete [] knapsackStart;
			delete [] knapsackRow;
			whichColumn=NULL;
			knapsackStart=NULL;
			knapsackRow=NULL;
		      }
		    }
		  }
		}
	      }
#endif
	      if (useCosts&&testOsiOptions<0) {
		int numberColumns = babModel->getNumCols();
		int * sort = new int[numberColumns];
		double * dsort = new double[numberColumns];
		int * priority = new int [numberColumns];
		const double * objective = babModel->getObjCoefficients();
		const double * lower = babModel->getColLower() ;
		const double * upper = babModel->getColUpper() ;
		const CoinPackedMatrix * matrix = babModel->solver()->getMatrixByCol();
		const int * columnLength = matrix->getVectorLengths();
		int iColumn;
		int n=0;
		for (iColumn=0;iColumn<numberColumns;iColumn++) {
		  if (babModel->isInteger(iColumn)) {
		    sort[n]=n;
		    if (useCosts==1)
		      dsort[n++]=-fabs(objective[iColumn]);
		    else if (useCosts==2)
		      dsort[n++]=iColumn;
		    else if (useCosts==3)
		      dsort[n++]=upper[iColumn]-lower[iColumn];
		    else if (useCosts==4)
		      dsort[n++]=-(upper[iColumn]-lower[iColumn]);
		    else if (useCosts==5)
		      dsort[n++]=-columnLength[iColumn];
		  }
		}
		CoinSort_2(dsort,dsort+n,sort);
		int level=0;
		double last = -1.0e100;
		for (int i=0;i<n;i++) {
		  int iPut=sort[i];
		  if (dsort[i]!=last) {
		    level++;
		    last=dsort[i];
		  }
		  priority[iPut]=level;
		}
		babModel->passInPriorities( priority,false);
		integersOK=true;
		delete [] priority;
		delete [] sort;
		delete [] dsort;
	      }
              // FPump done first as it only works if no solution
              CbcHeuristicFPump heuristic4(*babModel);
	      heuristic4.setFractionSmall(0.6);
              if (useFpump) {
                heuristic4.setMaximumPasses(parameters[whichParam(FPUMPITS,numberParameters,parameters)].intValue());
                int pumpTune=parameters[whichParam(FPUMPTUNE,numberParameters,parameters)].intValue();
		if (pumpTune>0) {
		  /*
		    >=10000000 for using obj
		    >=1000000 use as accumulate switch
		    >=1000 use index+1 as number of large loops
		    >=100 use 0.05 objvalue as increment
		    >=10 use +0.1 objvalue for cutoff (add)
		    1 == fix ints at bounds, 2 fix all integral ints, 3 and continuous at bounds
		    4 and static continuous, 5 as 3 but no internal integers
		    6 as 3 but all slack basis!
		  */
		  int logLevel = parameters[log].intValue();
		  double value = babModel->solver()->getObjSense()*babModel->solver()->getObjValue();
		  int w = pumpTune/10;
		  int c = w % 10;
		  w /= 10;
		  int i = w % 10;
		  w /= 10;
		  int r = w;
		  int accumulate = r/1000;
		  r -= 1000*accumulate;
		  if (accumulate>=10) {
		    int which = accumulate/10;
		    accumulate -= 10*which;
		    which--;
		    // weights and factors
		    double weight[]={0.1,0.1,0.5,0.5,1.0,1.0,5.0,5.0};
		    double factor[] = {0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5};
		    heuristic4.setInitialWeight(weight[which]);
		    heuristic4.setWeightFactor(factor[which]);
		  }
		  // fake cutoff
		  if (logLevel>1)
		    printf("Setting ");
		  if (c) {
		    double cutoff;
		    babModel->solver()->getDblParam(OsiDualObjectiveLimit,cutoff);
		    cutoff = CoinMin(cutoff,value + 0.1*fabs(value)*c);
		    double dextra1 = parameters[whichParam(DEXTRA1,numberParameters,parameters)].doubleValue();
		    if (dextra1)
		      cutoff=dextra1;
		    heuristic4.setFakeCutoff(cutoff);
		    if (logLevel>1)
		      printf("fake cutoff of %g ",cutoff);
		  }
		  if (i||r) {
		    // also set increment
		    double increment = (0.01*i+0.005)*(fabs(value)+1.0e-12);
		    double dextra2 = parameters[whichParam(DEXTRA2,numberParameters,parameters)].doubleValue();
		    if (dextra2)
		      increment = dextra2;
		    heuristic4.setAbsoluteIncrement(increment);
		    heuristic4.setAccumulate(accumulate);
		    heuristic4.setMaximumRetries(r+1);
		    if (logLevel>1) {
		      if (i) 
			printf("increment of %g ",heuristic4.absoluteIncrement());
		      if (accumulate) 
			printf("accumulate of %d ",accumulate);
		      printf("%d retries ",r+2);
		    }
		  }
		  pumpTune = pumpTune%100;
		  if (logLevel>1)
		    printf("and setting when to %d\n",pumpTune+10);
		  if (pumpTune==6)
		    pumpTune =13;
		  heuristic4.setWhen(pumpTune+10);
		}
		heuristic4.setHeuristicName("feasibility pump");
                babModel->addHeuristic(&heuristic4);
              }
              if (!miplib) {
                CbcRounding heuristic1(*babModel);
		heuristic1.setHeuristicName("rounding");
                if (useRounding)
                  babModel->addHeuristic(&heuristic1) ;
                CbcHeuristicLocal heuristic2(*babModel);
		heuristic2.setHeuristicName("combine solutions");
		heuristic2.setFractionSmall(0.6);
                heuristic2.setSearchType(1);
                if (useCombine)
                  babModel->addHeuristic(&heuristic2);
                CbcHeuristicGreedyCover heuristic3(*babModel);
		heuristic3.setHeuristicName("greedy cover");
                CbcHeuristicGreedyEquality heuristic3a(*babModel);
		heuristic3a.setHeuristicName("greedy equality");
                if (useGreedy) {
                  babModel->addHeuristic(&heuristic3);
                  babModel->addHeuristic(&heuristic3a);
                }
                if (useLocalTree) {
                  CbcTreeLocal localTree(babModel,NULL,10,0,0,10000,2000);
                  babModel->passInTreeHandler(localTree);
                }
              }
	      CbcHeuristicRINS heuristic5(*babModel);
	      heuristic5.setHeuristicName("RINS");
	      heuristic5.setFractionSmall(0.6);
	      if (useRINS)
		babModel->addHeuristic(&heuristic5) ;
	      if (type==MIPLIB) {
		if (babModel->numberStrong()==5&&babModel->numberBeforeTrust()==5) 
		  babModel->setNumberBeforeTrust(10);
	      }
              // add cut generators if wanted
              int switches[20];
              int numberGenerators=0;
	      int translate[]={-100,-1,-99,-98,1,1,1,1};
              if (probingAction) {
		if (probingAction==5||probingAction==7)
		  probingGen.setRowCuts(-3); // strengthening etc just at root
		if (probingAction==6||probingAction==7) {
		  // Number of unsatisfied variables to look at
		  probingGen.setMaxProbe(1000);
		  probingGen.setMaxProbeRoot(1000);
		  // How far to follow the consequences
		  probingGen.setMaxLook(50);
		  probingGen.setMaxLookRoot(50);
		}
                babModel->addCutGenerator(&probingGen,translate[probingAction],"Probing");
                switches[numberGenerators++]=0;
              }
              if (gomoryAction&&(complicatedInteger!=1||
				 (gomoryAction==1||gomoryAction==4))) {
                babModel->addCutGenerator(&gomoryGen,translate[gomoryAction],"Gomory");
                switches[numberGenerators++]=-1;
              }
              if (knapsackAction) {
                babModel->addCutGenerator(&knapsackGen,translate[knapsackAction],"Knapsack");
                switches[numberGenerators++]=0;
              }
              if (redsplitAction&&!complicatedInteger) {
                babModel->addCutGenerator(&redsplitGen,translate[redsplitAction],"Reduce-and-split");
                switches[numberGenerators++]=1;
              }
              if (cliqueAction) {
                babModel->addCutGenerator(&cliqueGen,translate[cliqueAction],"Clique");
                switches[numberGenerators++]=0;
              }
              if (mixedAction) {
                babModel->addCutGenerator(&mixedGen,translate[mixedAction],"MixedIntegerRounding2");
                switches[numberGenerators++]=-1;
              }
              if (flowAction) {
                babModel->addCutGenerator(&flowGen,translate[flowAction],"FlowCover");
                switches[numberGenerators++]=1;
              }
              if (twomirAction&&!complicatedInteger) {
                babModel->addCutGenerator(&twomirGen,translate[twomirAction],"TwoMirCuts");
                switches[numberGenerators++]=1;
              }
              if (landpAction) {
                babModel->addCutGenerator(&landpGen,translate[landpAction],"LiftAndProject");
                switches[numberGenerators++]=1;
              }
              if (residualCapacityAction) {
                babModel->addCutGenerator(&residualCapacityGen,translate[residualCapacityAction],"ResidualCapacity");
                switches[numberGenerators++]=1;
              }
	      if (storedCuts) 
		babModel->setSpecialOptions(babModel->specialOptions()|64);
              // Say we want timings
              numberGenerators = babModel->numberCutGenerators();
              int iGenerator;
              int cutDepth=
                parameters[whichParam(CUTDEPTH,numberParameters,parameters)].intValue();
              for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
                CbcCutGenerator * generator = babModel->cutGenerator(iGenerator);
                int howOften = generator->howOften();
                if (howOften==-98||howOften==-99) 
                  generator->setSwitchOffIfLessThan(switches[iGenerator]);
                generator->setTiming(true);
                if (cutDepth>=0)
                  generator->setWhatDepth(cutDepth) ;
              }
              // Could tune more
              if (!miplib) {
                babModel->setMinimumDrop(min(5.0e-2,
                                             fabs(babModel->getMinimizationObjValue())*1.0e-3+1.0e-4));
                if (cutPass==-1234567) {
                  if (babModel->getNumCols()<500)
                    babModel->setMaximumCutPassesAtRoot(-100); // always do 100 if possible
                  else if (babModel->getNumCols()<5000)
                    babModel->setMaximumCutPassesAtRoot(100); // use minimum drop
                  else
                    babModel->setMaximumCutPassesAtRoot(20);
                } else {
                  babModel->setMaximumCutPassesAtRoot(cutPass);
                }
                if (cutPassInTree==-1234567) 
		  babModel->setMaximumCutPasses(1);
		else
		  babModel->setMaximumCutPasses(cutPassInTree);
              }
              // Do more strong branching if small
              //if (babModel->getNumCols()<5000)
              //babModel->setNumberStrong(20);
              // Switch off strong branching if wanted
              //if (babModel->getNumCols()>10*babModel->getNumRows())
              //babModel->setNumberStrong(0);
              if (!noPrinting) {
		int iLevel = parameters[log].intValue();
		if (iLevel<0) {
		  babModel->setPrintingMode(1);
		  iLevel = -iLevel;
		}
                babModel->messageHandler()->setLogLevel(iLevel);
                if (babModel->getNumCols()>2000||babModel->getNumRows()>1500||
                    babModel->messageHandler()->logLevel()>1)
                  babModel->setPrintFrequency(100);
              }
              
              babModel->solver()->setIntParam(OsiMaxNumIterationHotStart,
                    parameters[whichParam(MAXHOTITS,numberParameters,parameters)].intValue());
              OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel->solver());
              // go faster stripes
              if (osiclp->getNumRows()<300&&osiclp->getNumCols()<500) {
                osiclp->setupForRepeatedUse(2,parameters[slog].intValue());
              } else {
                osiclp->setupForRepeatedUse(0,parameters[slog].intValue());
              }
              double increment=babModel->getCutoffIncrement();;
              int * changed = NULL;
              if (!miplib&&increment==normalIncrement)
                changed=analyze( osiclp,numberChanged,increment,false,generalMessageHandler);
              if (debugValues) {
                if (numberDebugValues==babModel->getNumCols()) {
                  // for debug
                  babModel->solver()->activateRowCutDebugger(debugValues) ;
                } else {
                  printf("debug file has incorrect number of columns\n");
                }
              }
              babModel->setCutoffIncrement(CoinMax(babModel->getCutoffIncrement(),increment));
              // Turn this off if you get problems
              // Used to be automatically set
              int mipOptions = parameters[whichParam(MIPOPTIONS,numberParameters,parameters)].intValue()%10000;
              if (mipOptions!=(1)) {
                sprintf(generalPrint,"mip options %d\n",mipOptions);
		generalMessageHandler->message(CLP_GENERAL,generalMessages)
		  << generalPrint
		  <<CoinMessageEol;
	      }
              osiclp->setSpecialOptions(mipOptions);
              if (gapRatio < 1.0e100) {
                double value = babModel->solver()->getObjValue() ;
                double value2 = gapRatio*(1.0e-5+fabs(value)) ;
                babModel->setAllowableGap(value2) ;
                std::cout << "Continuous " << value
                          << ", so allowable gap set to "
                          << value2 << std::endl ;
              }
              // probably faster to use a basis to get integer solutions
              babModel->setSpecialOptions(babModel->specialOptions()|2);
              currentBranchModel = babModel;
              OsiSolverInterface * strengthenedModel=NULL;
              if (type==BAB||type==MIPLIB) {
		int moreMipOptions = parameters[whichParam(MOREMIPOPTIONS,numberParameters,parameters)].intValue();
                if (moreMipOptions>=0) {
                  sprintf(generalPrint,"more mip options %d\n",moreMipOptions);
		  generalMessageHandler->message(CLP_GENERAL,generalMessages)
		    << generalPrint
		    <<CoinMessageEol;
		  OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel->solver());
		  if (moreMipOptions==10000) {
		    // test memory saving
		    moreMipOptions -= 10000;
		    ClpSimplex * lpSolver = osiclp->getModelPtr();
                    lpSolver->setPersistenceFlag(1);
		    // switch off row copy if few rows
		    if (lpSolver->numberRows()<150)
		      lpSolver->setSpecialOptions(lpSolver->specialOptions()|256);
		  }
		  if (((moreMipOptions+1)%1000000)!=0)
		    babModel->setSearchStrategy(moreMipOptions%1000000);
		  // go faster stripes
		  if( moreMipOptions >=999999) {
		    if (osiclp) {
		      int save = osiclp->specialOptions();
		      osiclp->setupForRepeatedUse(2,0);
		      osiclp->setSpecialOptions(save|osiclp->specialOptions());
		    }
		  } 
                }
              }
              if (type==BAB) {
#ifdef COIN_HAS_ASL
                if (usingAmpl) {
                  priorities=info.priorities;
                  branchDirection=info.branchDirection;
                  pseudoDown=info.pseudoDown;
                  pseudoUp=info.pseudoUp;
                  solutionIn=info.primalSolution;
                  prioritiesIn = info.priorities;
                  if (info.numberSos&&doSOS) {
                    // SOS
                    numberSOS = info.numberSos;
                    sosStart = info.sosStart;
                    sosIndices = info.sosIndices;
                    sosType = info.sosType;
                    sosReference = info.sosReference;
                    sosPriority = info.sosPriority;
                  }
                }
#endif                
                const int * originalColumns = preProcess ? process.originalColumns() : NULL;
                if (solutionIn&&useSolution) {
                  if (preProcess) {
                    int numberColumns = babModel->getNumCols();
                    // extend arrays in case SOS
                    int n = originalColumns[numberColumns-1]+1;
                    int nSmaller = CoinMin(n,numberOriginalColumns);
                    double * solutionIn2 = new double [n];
                    int * prioritiesIn2 = new int[n];
                    int i;
                    for (i=0;i<nSmaller;i++) {
                      solutionIn2[i]=solutionIn[i];
                      prioritiesIn2[i]=prioritiesIn[i];
                    }
                    for (;i<n;i++) {
                      solutionIn2[i]=0.0;
                      prioritiesIn2[i]=1000000;
                    }
                    int iLast=-1;
                    for (i=0;i<numberColumns;i++) {
                      int iColumn = originalColumns[i];
                      assert (iColumn>iLast);
                      iLast=iColumn;
                      solutionIn2[i]=solutionIn2[iColumn];
                      if (prioritiesIn)
                        prioritiesIn2[i]=prioritiesIn2[iColumn];
                    }
                    babModel->setHotstartSolution(solutionIn2,prioritiesIn2);
                    delete [] solutionIn2;
                    delete [] prioritiesIn2;
                  } else {
                    babModel->setHotstartSolution(solutionIn,prioritiesIn);
                  }
                }
		OsiSolverInterface * testOsiSolver= (testOsiOptions>=0) ? babModel->solver() : NULL;
		if (!testOsiSolver) {
		  // *************************************************************
		  // CbcObjects
		  if (preProcess&&(process.numberSOS()||babModel->numberObjects())) {
		    int numberSOS = process.numberSOS();
		    int numberIntegers = babModel->numberIntegers();
		    /* model may not have created objects
		       If none then create
		    */
		    if (!numberIntegers||!babModel->numberObjects()) {
		      int type = (pseudoUp) ? 1 : 0;
		      babModel->findIntegers(true,type);
		      numberIntegers = babModel->numberIntegers();
		      integersOK=true;
		    }
		    OsiObject ** oldObjects = babModel->objects();
		    // Do sets and priorities
		    OsiObject ** objects = new OsiObject * [numberSOS];
		    // set old objects to have low priority
		    int numberOldObjects = babModel->numberObjects();
		    int numberColumns = babModel->getNumCols();
		    // backward pointer to new variables
		    int * newColumn = new int[numberOriginalColumns];
		    int i;
		    for (i=0;i<numberOriginalColumns;i++)
		      newColumn[i]=-1;
		    assert (originalColumns);
		    for (i=0;i<numberColumns;i++)
		      newColumn[originalColumns[i]]=i;
		    if (!integersOK) {
		      // Change column numbers etc
		      int n=0;
		      for (int iObj = 0;iObj<numberOldObjects;iObj++) {
			int iColumn = oldObjects[iObj]->columnNumber();
			if (iColumn<0||iColumn>=numberOriginalColumns) {
			  oldObjects[n++]=oldObjects[iObj];
			} else {
			  iColumn = newColumn[iColumn];
			  if (iColumn>=0) {
			    CbcSimpleInteger * obj =
			      dynamic_cast <CbcSimpleInteger *>(oldObjects[iObj]) ;
			    assert (obj);
			    obj->setColumnNumber(iColumn);
			    oldObjects[n++]=oldObjects[iObj];
			  } else {
			    delete oldObjects[iObj];
			  }
			}
		      }
		      babModel->setNumberObjects(n);
		    }
		    int nMissing=0;
		    for (int iObj = 0;iObj<numberOldObjects;iObj++) {
		      if (process.numberSOS())
			oldObjects[iObj]->setPriority(numberColumns+1);
		      int iColumn = oldObjects[iObj]->columnNumber();
		      if (iColumn<0||iColumn>=numberOriginalColumns) {
			CbcSOS * obj =
			  dynamic_cast <CbcSOS *>(oldObjects[iObj]) ;
			if (obj) {
			  int n=obj->numberMembers();
			  int * which = obj->mutableMembers();
			  double * weights = obj->mutableWeights();
			  int nn=0;
			  for (i=0;i<n;i++) {
			    int iColumn = which[i];
			    int jColumn = newColumn[iColumn];
			    if (jColumn>=0) { 
			      which[nn] = jColumn;
			      weights[nn++]=weights[i];
			    } else {
			      nMissing++;
			    }
			  }
			  obj->setNumberMembers(nn);
			}
			continue;
		      }
		      if (originalColumns)
			iColumn = originalColumns[iColumn];
		      if (branchDirection) {
			CbcSimpleInteger * obj =
			  dynamic_cast <CbcSimpleInteger *>(oldObjects[iObj]) ;
			if (obj) { 
			  obj->setPreferredWay(branchDirection[iColumn]);
			} else {
			  CbcObject * obj =
			    dynamic_cast <CbcObject *>(oldObjects[iObj]) ;
			  assert (obj);
			  obj->setPreferredWay(branchDirection[iColumn]);
			}
		      }
		      if (pseudoUp) {
			CbcSimpleIntegerPseudoCost * obj1a =
			  dynamic_cast <CbcSimpleIntegerPseudoCost *>(oldObjects[iObj]) ;
			assert (obj1a);
			if (pseudoDown[iColumn]>0.0)
			  obj1a->setDownPseudoCost(pseudoDown[iColumn]);
			if (pseudoUp[iColumn]>0.0)
			  obj1a->setUpPseudoCost(pseudoUp[iColumn]);
		      }
		    }
		    if (nMissing) {
		      sprintf(generalPrint,"%d SOS variables vanished due to pre processing? - check validity?\n",nMissing);
		      generalMessageHandler->message(CLP_GENERAL,generalMessages)
			<< generalPrint
			<<CoinMessageEol;
		    }
		    delete [] newColumn;
		    const int * starts = process.startSOS();
		    const int * which = process.whichSOS();
		    const int * type = process.typeSOS();
		    const double * weight = process.weightSOS();
		    int iSOS;
		    for (iSOS =0;iSOS<numberSOS;iSOS++) {
		      int iStart = starts[iSOS];
		      int n=starts[iSOS+1]-iStart;
		      objects[iSOS] = new CbcSOS(babModel,n,which+iStart,weight+iStart,
						 iSOS,type[iSOS]);
		      // branch on long sets first
		      objects[iSOS]->setPriority(numberColumns-n);
		    }
		    if (numberSOS)
		      babModel->addObjects(numberSOS,objects);
		    for (iSOS=0;iSOS<numberSOS;iSOS++)
		      delete objects[iSOS];
		    delete [] objects;
		  } else if (priorities||branchDirection||pseudoDown||pseudoUp||numberSOS) {
		    // do anyway for priorities etc
		    int numberIntegers = babModel->numberIntegers();
		    /* model may not have created objects
		       If none then create
		    */
		    if (!numberIntegers||!babModel->numberObjects()) {
		      int type = (pseudoUp) ? 1 : 0;
		      babModel->findIntegers(true,type);
		    }
		    if (numberSOS) {
		      // Do sets and priorities
		      OsiObject ** objects = new OsiObject * [numberSOS];
		      int iSOS;
		      if (originalColumns) {
			// redo sequence numbers
			int numberColumns = babModel->getNumCols();
			int nOld = originalColumns[numberColumns-1]+1;
			int * back = new int[nOld];
			int i;
			for (i=0;i<nOld;i++)
			  back[i]=-1;
			for (i=0;i<numberColumns;i++)
			  back[originalColumns[i]]=i;
			// Really need better checks
			int nMissing=0;
			int n=sosStart[numberSOS];
			for (i=0;i<n;i++) {
			  int iColumn = sosIndices[i];
			  int jColumn = back[iColumn];
			  if (jColumn>=0) 
			    sosIndices[i] = jColumn;
			  else 
			    nMissing++;
			}
			delete [] back;
			if (nMissing) {
			  sprintf(generalPrint,"%d SOS variables vanished due to pre processing? - check validity?\n",nMissing);
			  generalMessageHandler->message(CLP_GENERAL,generalMessages)
			    << generalPrint
			    <<CoinMessageEol;
			}
		      }
		      for (iSOS =0;iSOS<numberSOS;iSOS++) {
			int iStart = sosStart[iSOS];
			int n=sosStart[iSOS+1]-iStart;
			objects[iSOS] = new CbcSOS(babModel,n,sosIndices+iStart,sosReference+iStart,
						   iSOS,sosType[iSOS]);
			if (sosPriority)
			  objects[iSOS]->setPriority(sosPriority[iSOS]);
			else if (!prioritiesIn)
			  objects[iSOS]->setPriority(10);  // rather than 1000 
		      }
		      // delete any existing SOS objects
		      int numberObjects=babModel->numberObjects();
		      OsiObject ** oldObjects=babModel->objects();
		      int nNew=0;
		      for (int i=0;i<numberObjects;i++) {
			OsiObject * objThis = oldObjects[i];
			CbcSOS * obj1 =
			  dynamic_cast <CbcSOS *>(objThis) ;
			OsiSOS * obj2 =
			  dynamic_cast <OsiSOS *>(objThis) ;
			if (!obj1&&!obj2) {
			  oldObjects[nNew++]=objThis;
			} else {
			  delete objThis;
			}
		      }
		      babModel->setNumberObjects(nNew);
		      babModel->addObjects(numberSOS,objects);
		      for (iSOS=0;iSOS<numberSOS;iSOS++)
			delete objects[iSOS];
		      delete [] objects;
		    }
                  }
                  OsiObject ** objects = babModel->objects();
                  int numberObjects = babModel->numberObjects();
                  for (int iObj = 0;iObj<numberObjects;iObj++) {
                    // skip sos
                    CbcSOS * objSOS =
                      dynamic_cast <CbcSOS *>(objects[iObj]) ;
                    if (objSOS)
                      continue;
                    int iColumn = objects[iObj]->columnNumber();
                    assert (iColumn>=0);
                    if (originalColumns)
                      iColumn = originalColumns[iColumn];
                    if (branchDirection) {
		      CbcSimpleInteger * obj =
			dynamic_cast <CbcSimpleInteger *>(objects[iObj]) ;
		      if (obj) { 
			obj->setPreferredWay(branchDirection[iColumn]);
		      } else {
			CbcObject * obj =
			  dynamic_cast <CbcObject *>(objects[iObj]) ;
			assert (obj);
			obj->setPreferredWay(branchDirection[iColumn]);
		      }
		    }
                    if (priorities) {
                      int iPriority = priorities[iColumn];
                      if (iPriority>0)
                        objects[iObj]->setPriority(iPriority);
                    }
                    if (pseudoUp&&pseudoUp[iColumn]) {
                      CbcSimpleIntegerPseudoCost * obj1a =
                        dynamic_cast <CbcSimpleIntegerPseudoCost *>(objects[iObj]) ;
                      assert (obj1a);
                      if (pseudoDown[iColumn]>0.0)
                        obj1a->setDownPseudoCost(pseudoDown[iColumn]);
                      if (pseudoUp[iColumn]>0.0)
                        obj1a->setUpPseudoCost(pseudoUp[iColumn]);
                    }
                  }
		  // *************************************************************
		} else {
		  // *************************************************************
		  // OsiObjects
		  // Find if none
		  int numberIntegers = testOsiSolver->getNumIntegers();
		  /* model may not have created objects
		     If none then create
		  */
		  if (!numberIntegers||!testOsiSolver->numberObjects()) {
		    //int type = (pseudoUp) ? 1 : 0;
		    testOsiSolver->findIntegers(false);
		    numberIntegers = testOsiSolver->getNumIntegers();
		  }
		  if (preProcess&&process.numberSOS()) {
		    int numberSOS = process.numberSOS();
		    OsiObject ** oldObjects = testOsiSolver->objects();
		    // Do sets and priorities
		    OsiObject ** objects = new OsiObject * [numberSOS];
		    // set old objects to have low priority
		    int numberOldObjects = testOsiSolver->numberObjects();
		    int numberColumns = testOsiSolver->getNumCols();
		    for (int iObj = 0;iObj<numberOldObjects;iObj++) {
		      oldObjects[iObj]->setPriority(numberColumns+1);
		      int iColumn = oldObjects[iObj]->columnNumber();
		      assert (iColumn>=0);
		      if (iColumn>=numberOriginalColumns)
			continue;
		      if (originalColumns)
			iColumn = originalColumns[iColumn];
		      if (branchDirection) {
			OsiSimpleInteger * obj =
			  dynamic_cast <OsiSimpleInteger *>(oldObjects[iObj]) ;
			if (obj) { 
			  obj->setPreferredWay(branchDirection[iColumn]);
			} else {
			  OsiObject2 * obj =
			    dynamic_cast <OsiObject2 *>(oldObjects[iObj]) ;
			  if (obj)
			    obj->setPreferredWay(branchDirection[iColumn]);
			}
		      }
		      if (pseudoUp) {
			abort();
		      }
		    }
		    const int * starts = process.startSOS();
		    const int * which = process.whichSOS();
		    const int * type = process.typeSOS();
		    const double * weight = process.weightSOS();
		    int iSOS;
		    for (iSOS =0;iSOS<numberSOS;iSOS++) {
		      int iStart = starts[iSOS];
		      int n=starts[iSOS+1]-iStart;
		      objects[iSOS] = new OsiSOS(testOsiSolver,n,which+iStart,weight+iStart,
						 type[iSOS]);
		      // branch on long sets first
		      objects[iSOS]->setPriority(numberColumns-n);
		    }
		    testOsiSolver->addObjects(numberSOS,objects);
		    for (iSOS=0;iSOS<numberSOS;iSOS++)
		      delete objects[iSOS];
		    delete [] objects;
		  } else if (priorities||branchDirection||pseudoDown||pseudoUp||numberSOS) {
		    if (numberSOS) {
		      // Do sets and priorities
		      OsiObject ** objects = new OsiObject * [numberSOS];
		      int iSOS;
		      if (originalColumns) {
			// redo sequence numbers
			int numberColumns = testOsiSolver->getNumCols();
			int nOld = originalColumns[numberColumns-1]+1;
			int * back = new int[nOld];
			int i;
			for (i=0;i<nOld;i++)
			  back[i]=-1;
			for (i=0;i<numberColumns;i++)
			  back[originalColumns[i]]=i;
			// Really need better checks
			int nMissing=0;
			int n=sosStart[numberSOS];
			for (i=0;i<n;i++) {
			  int iColumn = sosIndices[i];
			  int jColumn = back[iColumn];
			  if (jColumn>=0) 
			    sosIndices[i] = jColumn;
			  else 
			    nMissing++;
			}
			delete [] back;
			if (nMissing) {
			  sprintf(generalPrint,"%d SOS variables vanished due to pre processing? - check validity?\n",nMissing);
			  generalMessageHandler->message(CLP_GENERAL,generalMessages)
			    << generalPrint
			    <<CoinMessageEol;
			}
		      }
		      for (iSOS =0;iSOS<numberSOS;iSOS++) {
			int iStart = sosStart[iSOS];
			int n=sosStart[iSOS+1]-iStart;
			objects[iSOS] = new OsiSOS(testOsiSolver,n,sosIndices+iStart,sosReference+iStart,
						   sosType[iSOS]);
			if (sosPriority)
			  objects[iSOS]->setPriority(sosPriority[iSOS]);
			else if (!prioritiesIn)
			  objects[iSOS]->setPriority(10);  // rather than 1000 
		      }
		      // delete any existing SOS objects
		      int numberObjects=testOsiSolver->numberObjects();
		      OsiObject ** oldObjects=testOsiSolver->objects();
		      int nNew=0;
		      for (int i=0;i<numberObjects;i++) {
			OsiObject * objThis = oldObjects[i];
			OsiSOS * obj1 =
			  dynamic_cast <OsiSOS *>(objThis) ;
			OsiSOS * obj2 =
			  dynamic_cast <OsiSOS *>(objThis) ;
			if (!obj1&&!obj2) {
			  oldObjects[nNew++]=objThis;
			} else {
			  delete objThis;
			}
		      }
		      testOsiSolver->setNumberObjects(nNew);
		      testOsiSolver->addObjects(numberSOS,objects);
		      for (iSOS=0;iSOS<numberSOS;iSOS++)
			delete objects[iSOS];
		      delete [] objects;
		    }
                  }
                  OsiObject ** objects = testOsiSolver->objects();
                  int numberObjects = testOsiSolver->numberObjects();
		  int logLevel = parameters[log].intValue();
                  for (int iObj = 0;iObj<numberObjects;iObj++) {
                    // skip sos
                    OsiSOS * objSOS =
                      dynamic_cast <OsiSOS *>(objects[iObj]) ;
                    if (objSOS) {
		      if (logLevel>2)
			printf("Set %d is SOS - priority %d\n",iObj,objSOS->priority());
                      continue;
		    }
                    int iColumn = objects[iObj]->columnNumber();
                    if (iColumn>=0) {
		      if (originalColumns)
			iColumn = originalColumns[iColumn];
		      if (branchDirection) {
			OsiSimpleInteger * obj =
			  dynamic_cast <OsiSimpleInteger *>(objects[iObj]) ;
			if (obj) { 
			  obj->setPreferredWay(branchDirection[iColumn]);
			} else {
			  OsiObject2 * obj =
			    dynamic_cast <OsiObject2 *>(objects[iObj]) ;
			  if (obj)
			    obj->setPreferredWay(branchDirection[iColumn]);
			}
		      }
		      if (priorities) {
			int iPriority = priorities[iColumn];
			if (iPriority>0)
			  objects[iObj]->setPriority(iPriority);
		      }
		      if (logLevel>2)
			printf("Obj %d is int? - priority %d\n",iObj,objects[iObj]->priority());
		      if (pseudoUp&&pseudoUp[iColumn]) {
			abort();
		      }
                    }
                  }
		  // *************************************************************
                }
                int statistics = (printOptions>0) ? printOptions: 0;
#ifdef COIN_HAS_ASL
                if (!usingAmpl) {
#endif
                  free(priorities);
                  priorities=NULL;
                  free(branchDirection);
                  branchDirection=NULL;
                  free(pseudoDown);
                  pseudoDown=NULL;
                  free(pseudoUp);
                  pseudoUp=NULL;
                  free(solutionIn);
                  solutionIn=NULL;
                  free(prioritiesIn);
                  prioritiesIn=NULL;
                  free(sosStart);
                  sosStart=NULL;
                  free(sosIndices);
                  sosIndices=NULL;
                  free(sosType);
                  sosType=NULL;
                  free(sosReference);
                  sosReference=NULL;
		  free(cut);
		  cut=NULL;
                  free(sosPriority);
                  sosPriority=NULL;
#ifdef COIN_HAS_ASL
                }
#endif                
		if (nodeStrategy) {
		  // change default
		  if (nodeStrategy>2) {
		    // up or down
		    int way = (((nodeStrategy-1)%1)==1) ? -1 : +1;
		    babModel->setPreferredWay(way);
#if 0
		    OsiObject ** objects = babModel->objects();
		    int numberObjects = babModel->numberObjects();
		    for (int iObj = 0;iObj<numberObjects;iObj++) {
		      CbcObject * obj =
			dynamic_cast <CbcObject *>(objects[iObj]) ;
		      assert (obj);
		      obj->setPreferredWay(way);
		    }
#endif
		  }
		  if (nodeStrategy==2||nodeStrategy>4) {
		    // depth
		    CbcCompareDefault compare;
		    compare.setWeight(-3.0);
		    babModel->setNodeComparison(compare);
		  } else if (nodeStrategy==0) {
		    // hybrid was default i.e. mixture of low depth and infeasibility
		  } else if (nodeStrategy==1) {
		    // real fewest
		    CbcCompareDefault compare;
		    compare.setWeight(-2.0);
		    babModel->setNodeComparison(compare);
		  }
		}
	        if (cppValue>=0) {
		  int prepro = useStrategy ? -1 : preProcess;
                  // generate code
                  FILE * fp = fopen("user_driver.cpp","w");
	          if (fp) {
	            // generate enough to do BAB
		    babModel->generateCpp(fp,1);
                    OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel->solver());
	            // Make general so do factorization
                    int factor = osiclp->getModelPtr()->factorizationFrequency();
                    osiclp->getModelPtr()->setFactorizationFrequency(200);
                    osiclp->generateCpp(fp);
                    osiclp->getModelPtr()->setFactorizationFrequency(factor);
                    //solveOptions.generateCpp(fp);
                    fclose(fp);
                    // now call generate code
                    generateCode(babModel,"user_driver.cpp",cppValue,prepro);
                  } else {
                    std::cout<<"Unable to open file user_driver.cpp"<<std::endl;
                  }
                }
		if (!babModel->numberStrong())
		  babModel->setNumberBeforeTrust(0);
		if (useStrategy) {
		  CbcStrategyDefault strategy(true,babModel->numberStrong(),babModel->numberBeforeTrust());
                  strategy.setupPreProcessing(1);
		  babModel->setStrategy(strategy);
		}
                if (testOsiOptions>=0) {
                  sprintf(generalPrint,"Testing OsiObject options %d\n",testOsiOptions);
		  generalMessageHandler->message(CLP_GENERAL,generalMessages)
		    << generalPrint
		    <<CoinMessageEol;
		  if (!numberSOS) {
		    babModel->solver()->findIntegersAndSOS(false);
#ifdef COIN_HAS_LINK
		    // If linked then pass in model
		    OsiSolverLink * solver3 = dynamic_cast<OsiSolverLink *> (babModel->solver());
		    if (solver3) {
		      CbcHeuristicDynamic3 serendipity(*babModel);
		      serendipity.setHeuristicName("linked");
		      babModel->addHeuristic(&serendipity);
		      int options = parameters[whichParam(MIPOPTIONS,numberParameters,parameters)].intValue()/10000;
		      CglStored stored;
		      if (options) {
			printf("nlp options %d\n",options);
			/*
			  1 - force mini branch and bound
			  2 - set priorities high on continuous
			  4 - try adding OA cuts
			  8 - try doing quadratic linearization
			  16 - try expanding knapsacks
                          32 - OA cuts strictly concave
			  64 - no branching at all on bilinear x-x!
			*/
			if ((options&2)) {
			  solver3->setBiLinearPriorities(10,tightenFactor > 0.0 ? tightenFactor : 1.0);
			} else if (tightenFactor>0.0) {
			  // set grid size for all continuous bi-linear
			  solver3->setMeshSizes(tightenFactor);
			}
			if ((options&4)) {
			  solver3->setSpecialOptions2(solver3->specialOptions2()|(8+4));
			  // say convex
			  solver3->sayConvex((options&32)==0);
			}
			int extra1 = parameters[whichParam(EXTRA1,numberParameters,parameters)].intValue();
			if ((options&1)!=0&&extra1>0)
			  solver3->setFixedPriority(extra1);
			double cutoff=COIN_DBL_MAX;
			if ((options&8))
			  cutoff=solver3->linearizedBAB(&stored);
			if (cutoff<babModel->getCutoff()) {
			  babModel->setCutoff(cutoff);
			  // and solution
			  //babModel->setBestObjectiveValue(solver3->bestObjectiveValue());
			  babModel->setBestSolution(solver3->bestSolution(),solver3->getNumCols(),
						    solver3->bestObjectiveValue());
			}
			if ((options&64))
			  solver3->setBranchingStrategyOnVariables(16,-1,4);
		      }
		      solver3->setCbcModel(babModel);
		      if (stored.sizeRowCuts()) 
			babModel->addCutGenerator(&stored,1,"Stored");
		      CglTemporary temp;
		      babModel->addCutGenerator(&temp,1,"OnceOnly");
		      //choose.setNumberBeforeTrusted(2000);
		      //choose.setNumberStrong(20);
		    }
		    // For temporary testing of heuristics
		    //int testOsiOptions = parameters[whichParam(TESTOSI,numberParameters,parameters)].intValue();
		    if (testOsiOptions>=10) {
		      if (testOsiOptions>=20)
			testOsiOptions -= 10;
		      printf("*** Temp heuristic with mode %d\n",testOsiOptions-10);
		      OsiSolverLink * solver3 = dynamic_cast<OsiSolverLink *> (babModel->solver());
		      assert (solver3) ;
		      int extra1 = parameters[whichParam(EXTRA1,numberParameters,parameters)].intValue();
		      solver3->setBiLinearPriority(extra1);
		      printf("bilinear priority now %d\n",extra1);
		      int extra2 = parameters[whichParam(EXTRA2,numberParameters,parameters)].intValue();
		      double saveDefault = solver3->defaultBound();
		      solver3->setDefaultBound((double) extra2);
		      double * solution = solver3->heuristicSolution(slpValue>0 ? slpValue : 40 ,1.0e-5,testOsiOptions-10);
		      solver3->setDefaultBound(saveDefault);
		      if (!solution)
			printf("Heuristic failed\n");
		    }
#endif
		  } else {
		    // move across
		    babModel->deleteObjects(false);
		    //babModel->addObjects(babModel->solver()->numberObjects(),babModel->solver()->objects());
		  }
		  CbcBranchDefaultDecision decision;
		  if (babModel->numberStrong()) {
		    OsiChooseStrong choose(babModel->solver());
		    choose.setNumberBeforeTrusted(babModel->numberBeforeTrust());
		    choose.setNumberStrong(babModel->numberStrong());
		    choose.setShadowPriceMode(testOsiOptions);
		    decision.setChooseMethod(choose);
		  } else {
		    OsiChooseVariable choose(babModel->solver());
		    decision.setChooseMethod(choose);
		  }
		  babModel->setBranchingMethod(decision);
		  if (useCosts&&testOsiOptions>=0) {
		    int numberColumns = babModel->getNumCols();
		    int * sort = new int[numberColumns];
		    double * dsort = new double[numberColumns];
		    int * priority = new int [numberColumns];
		    const double * objective = babModel->getObjCoefficients();
		    const double * lower = babModel->getColLower() ;
		    const double * upper = babModel->getColUpper() ;
		    const CoinPackedMatrix * matrix = babModel->solver()->getMatrixByCol();
		    const int * columnLength = matrix->getVectorLengths();
		    int iColumn;
		    for (iColumn=0;iColumn<numberColumns;iColumn++) {
		      sort[iColumn]=iColumn;
		      if (useCosts==1)
			dsort[iColumn]=-fabs(objective[iColumn]);
		      else if (useCosts==2)
			dsort[iColumn]=iColumn;
		      else if (useCosts==3)
			dsort[iColumn]=upper[iColumn]-lower[iColumn];
		      else if (useCosts==4)
			dsort[iColumn]=-(upper[iColumn]-lower[iColumn]);
		      else if (useCosts==5)
			dsort[iColumn]=-columnLength[iColumn];
		    }
		    CoinSort_2(dsort,dsort+numberColumns,sort);
		    int level=0;
		    double last = -1.0e100;
		    for (int i=0;i<numberColumns;i++) {
		      int iPut=sort[i];
		      if (dsort[i]!=last) {
			level++;
			last=dsort[i];
		      }
		      priority[iPut]=level;
		    }
		    OsiObject ** objects = babModel->objects();
		    int numberObjects = babModel->numberObjects();
		    for (int iObj = 0;iObj<numberObjects;iObj++) {
		      OsiObject * obj = objects[iObj] ;
		      int iColumn = obj->columnNumber();
		      if (iColumn>=0)
			obj->setPriority(priority[iColumn]);
		    }
		    delete [] priority;
		    delete [] sort;
		    delete [] dsort;
		  }
		}
		checkSOS(babModel, babModel->solver());
		if (doSprint>0) {
		  // Sprint for primal solves
		  ClpSolve::SolveType method = ClpSolve::usePrimalorSprint;
		  ClpSolve::PresolveType presolveType = ClpSolve::presolveOff;
		  int numberPasses = 5;
		  int options[] = {0,3,0,0,0,0};
		  int extraInfo[] = {-1,20,-1,-1,-1,-1};
		  extraInfo[1]=doSprint;
		  int independentOptions[] = {0,0,3};
		  ClpSolve clpSolve(method,presolveType,numberPasses,
				    options,extraInfo,independentOptions);
		  // say use in OsiClp
		  clpSolve.setSpecialOption(6,1);
		  OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel->solver());
		  osiclp->setSolveOptions(clpSolve);
		  osiclp->setHintParam(OsiDoDualInResolve,false);
		  // switch off row copy
		  osiclp->getModelPtr()->setSpecialOptions(osiclp->getModelPtr()->specialOptions()|256);
		  osiclp->getModelPtr()->setInfeasibilityCost(1.0e11);
		}
#ifdef COIN_HAS_LINK
		if (storedAmpl.sizeRowCuts()) {
		  if (preProcess) {
		    const int * originalColumns = process.originalColumns();
                    int numberColumns = babModel->getNumCols();
                    int * newColumn = new int[numberOriginalColumns];
                    int i;
                    for (i=0;i<numberOriginalColumns;i++) 
		      newColumn[i]=-1;
                    for (i=0;i<numberColumns;i++) {
                      int iColumn = originalColumns[i];
		      newColumn[iColumn]=i;
                    }
		    int * buildColumn = new int[numberColumns];
		    // Build up valid cuts
		    int nBad=0;
		    int nCuts = storedAmpl.sizeRowCuts();
		    CglStored newCuts;
		    for (i=0;i<nCuts;i++) {
		      const OsiRowCut * cut = storedAmpl.rowCutPointer(i);
		      double lb = cut->lb();
		      double ub = cut->ub();
		      int n=cut->row().getNumElements();
		      const int * column = cut->row().getIndices();
		      const double * element = cut->row().getElements();
		      bool bad=false;
		      for (int i=0;i<n;i++) {
			int iColumn = column[i];
			iColumn = newColumn[iColumn];
			if (iColumn>=0) {
			  buildColumn[i]=iColumn;
			} else {
			  bad=true;
			  break;
			}
		      }
		      if (!bad) {
			newCuts.addCut(lb,ub,n,buildColumn,element);
		      } else {
			nBad++;
		      }
		    }
		    storedAmpl=newCuts;
		    if (nBad)
		      printf("%d cuts dropped\n",nBad);
		    delete [] newColumn;
		    delete [] buildColumn;
		  }
		}
#endif
#ifdef CLP_MALLOC_STATISTICS
		malloc_stats();
		malloc_stats2();
#endif
		if (outputFormat==5) {
		  osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel->solver());
		  lpSolver = osiclp->getModelPtr();
		  lpSolver->setPersistenceFlag(1);
		}
#ifdef COIN_HAS_ASL
		// add in lotsizing
		if (usingAmpl&&info.special) {
		  int numberColumns = babModel->getNumCols();
		  int i;
		  int n=0;
		  if (preProcess) {
		    const int * originalColumns = process.originalColumns();
                    for (i=0;i<numberColumns;i++) {
                      int iColumn = originalColumns[i];
		      assert (iColumn>=i);
		      int iType = info.special[iColumn];
		      if (iType) {
			assert (iType==1);
			n++;
		      }
		      info.special[i]=iType;
                    }
		  }
		  if (n) {
		    int numberIntegers=0;
		    int numberOldObjects=0;
		    OsiObject ** oldObjects=NULL;
		    const double * lower = babModel->solver()->getColLower();
		    const double * upper = babModel->solver()->getColUpper();
		    if (testOsiOptions<0) {
		      // *************************************************************
		      // CbcObjects
		      numberIntegers = babModel->numberIntegers();
		      /* model may not have created objects
			 If none then create
		      */
		      if (!numberIntegers||!babModel->numberObjects()) {
			int type = (pseudoUp) ? 1 : 0;
			babModel->findIntegers(true,type);
			numberIntegers = babModel->numberIntegers();
		      }
		      oldObjects = babModel->objects();
		      numberOldObjects = babModel->numberObjects();
		    } else {
		      numberIntegers = testOsiSolver->getNumIntegers();
		      if (!numberIntegers||!testOsiSolver->numberObjects()) {
			/* model may not have created objects
			   If none then create
			*/
			testOsiSolver->findIntegers(false);
			numberIntegers = testOsiSolver->getNumIntegers();
		      }
		      oldObjects = testOsiSolver->objects();
		      numberOldObjects = testOsiSolver->numberObjects();
		    }
		    OsiObject ** objects = new OsiObject * [n];
		    n=0;
		    // set new objects to have one lower priority
		    double ranges[] = {-COIN_DBL_MAX,-1.0,1.0,COIN_DBL_MAX};
		    for (int iObj = 0;iObj<numberOldObjects;iObj++) {
		      int iColumn = oldObjects[iObj]->columnNumber();
		      if (iColumn>=0&&info.special[iColumn]) {
			if (lower[iColumn]<=-1.0&&upper[iColumn]>=0.0) {
			  ranges[0]=lower[iColumn];
			  ranges[3]=upper[iColumn];
			  int priority = oldObjects[iObj]->priority();
			  if (testOsiOptions<0) {
			    objects[n] = new CbcLotsize(babModel,iColumn,2,ranges,true);
			  } else {
			    objects[n] = new OsiLotsize(testOsiSolver,iColumn,2,ranges,true);
			  }
			  objects[n++]->setPriority (priority-1);
			}
		      }
		    }
		    if (testOsiOptions<0) {
		      babModel->addObjects(n,objects);
		    } else {
		      testOsiSolver->addObjects(n,objects);
		    }
		    for (i=0;i<n;i++)
		      delete objects[i];
		    delete [] objects;
		  }
		}
#endif
		if (storedAmpl.sizeRowCuts()) {
		  //babModel->addCutGenerator(&storedAmpl,1,"AmplStored");
		  int numberRowCuts = storedAmpl.sizeRowCuts();
		  for (int i=0;i<numberRowCuts;i++) {
		    const OsiRowCut * rowCutPointer = storedAmpl.rowCutPointer(i);
		    babModel->makeGlobalCut(rowCutPointer);
		  }
		}
		// If defaults then increase trust for small models
		if (!strongChanged) {
		  int numberColumns = babModel->getNumCols();
		  if (numberColumns<=50)
		    babModel->setNumberBeforeTrust(1000);
		  else if (numberColumns<=100)
		    babModel->setNumberBeforeTrust(100);
		  else if (numberColumns<=300)
		    babModel->setNumberBeforeTrust(50);
		}
#ifdef CBC_THREAD
                int numberThreads =parameters[whichParam(THREADS,numberParameters,parameters)].intValue();
		babModel->setNumberThreads(numberThreads%100);
		babModel->setThreadMode(numberThreads/100);
#endif
                babModel->branchAndBound(statistics);
#ifdef CLP_MALLOC_STATISTICS
		malloc_stats();
		malloc_stats2();
#endif
		checkSOS(babModel, babModel->solver());
              } else if (type==MIPLIB) {
		CbcStrategyDefault strategy(true,babModel->numberStrong(),babModel->numberBeforeTrust());
                // Set up pre-processing 
		int translate2[]={9999,1,1,3,2,4,5};
                if (preProcess)
                  strategy.setupPreProcessing(translate2[preProcess]);
                babModel->setStrategy(strategy);
#ifdef CBC_THREAD
                int numberThreads =parameters[whichParam(THREADS,numberParameters,parameters)].intValue();
		babModel->setNumberThreads(numberThreads%100);
		babModel->setThreadMode(numberThreads/100);
#endif
		if (outputFormat==5) {
		  osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel->solver());
		  lpSolver = osiclp->getModelPtr();
		  lpSolver->setPersistenceFlag(1);
		}
                if (testOsiOptions>=0) {
                  printf("Testing OsiObject options %d\n",testOsiOptions);
		  CbcBranchDefaultDecision decision;
		  OsiChooseStrong choose(babModel->solver());
		  choose.setNumberBeforeTrusted(babModel->numberBeforeTrust());
		  choose.setNumberStrong(babModel->numberStrong());
		  choose.setShadowPriceMode(testOsiOptions);
		  //babModel->deleteObjects(false);
		  decision.setChooseMethod(choose);
		  babModel->setBranchingMethod(decision);
		}
		model = *babModel;
		return 777;
              } else {
                strengthenedModel = babModel->strengthenedModel();
              }
              currentBranchModel = NULL;
              osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel->solver());
              if (debugFile=="createAfterPre"&&babModel->bestSolution()) {
                lpSolver = osiclp->getModelPtr();
                //move best solution (should be there -- but ..)
                int n = lpSolver->getNumCols();
                memcpy(lpSolver->primalColumnSolution(),babModel->bestSolution(),n*sizeof(double));
                saveSolution(osiclp->getModelPtr(),"debug.file");
              }
              if (!noPrinting) {
                // Print more statistics
		sprintf(generalPrint,"Cuts at root node changed objective from %g to %g",
			babModel->getContinuousObjective(),babModel->rootObjectiveAfterCuts());
		generalMessageHandler->message(CLP_GENERAL,generalMessages)
		  << generalPrint
		  <<CoinMessageEol;
                
		numberGenerators = babModel->numberCutGenerators();
		char timing[30];
                for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
                  CbcCutGenerator * generator = babModel->cutGenerator(iGenerator);
		  sprintf(generalPrint,"%s was tried %d times and created %d cuts of which %d were active after adding rounds of cuts",
			  generator->cutGeneratorName(),
			  generator->numberTimesEntered(),
			  generator->numberCutsInTotal(),
			  generator->numberCutsActive());
                  if (generator->timing()) {
		    sprintf(timing," (%.3f seconds)",generator->timeInCutGenerator());
		    strcat(generalPrint,timing);
		  }
		  generalMessageHandler->message(CLP_GENERAL,generalMessages)
		    << generalPrint
		    <<CoinMessageEol;
                }
              }
	      // adjust time to allow for children on some systems
              time2 = CoinCpuTime() + CoinCpuTimeJustChildren();
              totalTime += time2-time1;
              // For best solution
              double * bestSolution = NULL;
              if (babModel->getMinimizationObjValue()<1.0e50&&type==BAB) {
                // post process
		int n;
                if (preProcess) {
                  n = saveSolver->getNumCols();
                  bestSolution = new double [n];
		  OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (babModel->solver());
		  ClpSimplex * lpSolver = clpSolver->getModelPtr();
		  lpSolver->setSpecialOptions(lpSolver->specialOptions()|0x01000000); // say is Cbc (and in branch and bound)
                  process.postProcess(*babModel->solver());
                  // Solution now back in saveSolver
                  babModel->assignSolver(saveSolver);
                  memcpy(bestSolution,babModel->solver()->getColSolution(),n*sizeof(double));
                } else {
                  n = babModel->solver()->getNumCols();
                  bestSolution = new double [n];
                  memcpy(bestSolution,babModel->solver()->getColSolution(),n*sizeof(double));
                }
		model.setBestSolution(bestSolution,n,babModel->getObjValue());
		// and put back in very original solver
		{
		  ClpSimplex * original = originalSolver->getModelPtr();
		  double * lower = original->columnLower();
		  double * upper = original->columnUpper();
		  double * solution = original->primalColumnSolution();
		  int n = original->numberColumns();
		  //assert (!n||n==babModel->solver()->getNumCols());
		  for (int i=0;i<n;i++) {
		    solution[i]=bestSolution[i];
		    if (originalSolver->isInteger(i)) {
		      lower[i]=solution[i];
		      upper[i]=solution[i];
		    }
		  }
		}
		checkSOS(babModel, babModel->solver());
              }
              if (type==STRENGTHEN&&strengthenedModel)
                clpSolver = dynamic_cast< OsiClpSolverInterface*> (strengthenedModel);
#ifdef COIN_HAS_ASL
	      else if (usingAmpl) 
		clpSolver = dynamic_cast< OsiClpSolverInterface*> (babModel->solver());
#endif
              lpSolver = clpSolver->getModelPtr();
              if (numberChanged) {
                for (int i=0;i<numberChanged;i++) {
                  int iColumn=changed[i];
                  clpSolver->setContinuous(iColumn);
                }
                delete [] changed;
              }
              if (type==BAB) {
                //move best solution (should be there -- but ..)
                int n = lpSolver->getNumCols();
                if (bestSolution) {
                  memcpy(lpSolver->primalColumnSolution(),bestSolution,n*sizeof(double));
		  // now see what that does to row solution
		  int numberRows=lpSolver->numberRows();
		  double * rowSolution = lpSolver->primalRowSolution();
		  memset (rowSolution,0,numberRows*sizeof(double));
		  lpSolver->clpMatrix()->times(1.0,bestSolution,rowSolution);
		  lpSolver->setObjectiveValue(babModel->getObjValue());
		}
                if (debugFile=="create"&&bestSolution) {
                  saveSolution(lpSolver,"debug.file");
                }
                delete [] bestSolution;
                std::string statusName[]={"Finished","Stopped on ","Difficulties",
                                          "","","User ctrl-c"};
                std::string minor[]={"","","gap","nodes","time","","solutions","user ctrl-c"};
                int iStat = babModel->status();
                int iStat2 = babModel->secondaryStatus();
                if (!noPrinting) {
		  sprintf(generalPrint,"Result - %s%s objective %.16g after %d nodes and %d iterations - took %.2f seconds",
			  statusName[iStat].c_str(),minor[iStat2].c_str(),
                          babModel->getObjValue(),babModel->getNodeCount(),
                          babModel->getIterationCount(),time2-time1);
		  generalMessageHandler->message(CLP_GENERAL,generalMessages)
		    << generalPrint
		    <<CoinMessageEol;
		}
#ifdef COIN_HAS_ASL
                if (usingAmpl) {
		  clpSolver = dynamic_cast< OsiClpSolverInterface*> (babModel->solver());
		  lpSolver = clpSolver->getModelPtr();
                  double value = babModel->getObjValue()*lpSolver->getObjSense();
                  char buf[300];
                  int pos=0;
                  if (iStat==0) {
                    if (babModel->getObjValue()<1.0e40) {
                      pos += sprintf(buf+pos,"optimal," );
                    } else {
                      // infeasible
                      iStat=1;
                      pos += sprintf(buf+pos,"infeasible,");
                    }
                  } else if (iStat==1) {
                    if (iStat2!=6)
                      iStat=3;
                    else
                      iStat=4;
                    pos += sprintf(buf+pos,"stopped on %s,",minor[iStat2].c_str());
                  } else if (iStat==2) {
                    iStat = 7;
                    pos += sprintf(buf+pos,"stopped on difficulties,");
                  } else if (iStat==5) {
                    iStat = 3;
                    pos += sprintf(buf+pos,"stopped on ctrl-c,");
                  } else {
                    pos += sprintf(buf+pos,"status unknown,");
                    iStat=6;
                  }
                  info.problemStatus=iStat;
                  info.objValue = value;
                  if (babModel->getObjValue()<1.0e40) {
		    int precision = ampl_obj_prec();
		    if (precision>0)
		      pos += sprintf(buf+pos," objective %.*g",precision,
				     value);
		    else
		      pos += sprintf(buf+pos," objective %g",value);
		  }
                  sprintf(buf+pos,"\n%d nodes, %d iterations, %g seconds",
                          babModel->getNodeCount(),
                          babModel->getIterationCount(),
			  totalTime);
                  if (bestSolution) {
                    free(info.primalSolution);
		    if (!numberKnapsack) {
		      info.primalSolution = (double *) malloc(n*sizeof(double));
		      CoinCopyN(lpSolver->primalColumnSolution(),n,info.primalSolution);
		      int numberRows = lpSolver->numberRows();
		      free(info.dualSolution);
		      info.dualSolution = (double *) malloc(numberRows*sizeof(double));
		      CoinCopyN(lpSolver->dualRowSolution(),numberRows,info.dualSolution);
		    } else {
		      // expanded knapsack
		      info.dualSolution=NULL;
		      int numberColumns = saveCoinModel.numberColumns();
		      info.primalSolution = (double *) malloc(numberColumns*sizeof(double));
		      // Fills in original solution (coinModel length)
		      afterKnapsack(saveTightenedModel,  whichColumn,  knapsackStart, 
				    knapsackRow,  numberKnapsack,
				    lpSolver->primalColumnSolution(), info.primalSolution,1);
		    }
                  } else {
                    info.primalSolution=NULL;
                    info.dualSolution=NULL;
                  }
                  // put buffer into info
                  strcpy(info.buffer,buf);
                }
#endif
              } else {
                std::cout<<"Model strengthened - now has "<<clpSolver->getNumRows()
                         <<" rows"<<std::endl;
              }
              time1 = time2;
#ifdef COIN_HAS_ASL
              if (usingAmpl) {
		// keep if going to be destroyed
		OsiSolverInterface * solver = babModel->solver();
		OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
		ClpSimplex * lpSolver2 = clpSolver->getModelPtr();
		if (lpSolver==lpSolver2)
		  babModel->setModelOwnsSolver(false);
	      }
#endif
              //delete babModel;
              //babModel=NULL;
            } else {
              std::cout << "** Current model not valid" << std::endl ; 
            }
            break ;
	  case IMPORT:
	    {
#ifdef COIN_HAS_ASL
              if (!usingAmpl) {
#endif
                free(priorities);
                priorities=NULL;
                free(branchDirection);
                branchDirection=NULL;
                free(pseudoDown);
                pseudoDown=NULL;
                free(pseudoUp);
                pseudoUp=NULL;
                free(solutionIn);
                solutionIn=NULL;
                free(prioritiesIn);
                prioritiesIn=NULL;
                free(sosStart);
                sosStart=NULL;
                free(sosIndices);
                sosIndices=NULL;
                free(sosType);
                sosType=NULL;
                free(sosReference);
                sosReference=NULL;
		free(cut);
		cut=NULL;
                free(sosPriority);
                sosPriority=NULL;
#ifdef COIN_HAS_ASL
              }
#endif                
              //delete babModel;
              //babModel=NULL;
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      bool canOpen=false;
              // See if gmpl file
              int gmpl=0;
              std::string gmplData;
	      if (field=="-") {
		// stdin
		canOpen=true;
		fileName = "-";
	      } else {
		// See if .lp
		{
		  const char * c_name = field.c_str();
		  int length = strlen(c_name);
		  if (length>3&&!strncmp(c_name+length-3,".lp",3))
		    gmpl=-1; // .lp
		}
                bool absolutePath;
                if (dirsep=='/') {
                  // non Windows (or cygwin)
                  absolutePath=(field[0]=='/');
                } else {
                  //Windows (non cycgwin)
                  absolutePath=(field[0]=='\\');
                  // but allow for :
                  if (strchr(field.c_str(),':'))
                    absolutePath=true;
                }
		if (absolutePath) {
		  fileName = field;
		} else if (field[0]=='~') {
		  char * environVar = getenv("HOME");
		  if (environVar) {
		    std::string home(environVar);
		    field=field.erase(0,1);
		    fileName = home+field;
		  } else {
		    fileName=field;
		  }
		} else {
		  fileName = directory+field;
                  // See if gmpl (model & data) - or even lp file
                  int length = field.size();
                  int percent = field.find('%');
                  if (percent<length&&percent>0) {
                    gmpl=1;
                    fileName = directory+field.substr(0,percent);
                    gmplData = directory+field.substr(percent+1);
                    if (percent<length-1)
                      gmpl=2; // two files
                    printf("GMPL model file %s and data file %s\n",
                           fileName.c_str(),gmplData.c_str());
		  }
		}
                std::string name=fileName;
                if (fileCoinReadable(name)) {
		  // can open - lets go for it
		  canOpen=true;
                  if (gmpl==2) {
                    FILE *fp;
                    fp=fopen(gmplData.c_str(),"r");
                    if (fp) {
                      fclose(fp);
                    } else {
                      canOpen=false;
                      std::cout<<"Unable to open file "<<gmplData<<std::endl;
                    }
                  }
		} else {
		  std::cout<<"Unable to open file "<<fileName<<std::endl;
		}
	      }
	      if (canOpen) {
                int status;
		ClpSimplex * lpSolver = clpSolver->getModelPtr();
                if (!gmpl) {
                  status =clpSolver->readMps(fileName.c_str(),
                                                 keepImportNames!=0,
                                                 allowImportErrors!=0);
		} else if (gmpl>0) {
                  status= lpSolver->readGMPL(fileName.c_str(),
					     (gmpl==2) ? gmplData.c_str() : NULL,
					     keepImportNames!=0);
		} else {
                  status= lpSolver->readLp(fileName.c_str(),1.0e-12);
		}
		if (!status||(status>0&&allowImportErrors)) {
		  if (keepImportNames&&gmpl<=0) {
		    lengthName = lpSolver->lengthNames();
		    rowNames = *(lpSolver->rowNames());
		    columnNames = *(lpSolver->columnNames());
		  } else {
		    lengthName=0;
		  }
		  goodModel=true;
		  // sets to all slack (not necessary?)
		  lpSolver->createStatus();
		  // make sure integer
		  int numberColumns = lpSolver->numberColumns();
		  for (int i=0;i<numberColumns;i++) {
		    if (lpSolver->isInteger(i))
		      clpSolver->setInteger(i);
		  }
		  time2 = CoinCpuTime();
		  totalTime += time2-time1;
		  time1=time2;
		  // Go to canned file if just input file
		  if (CbcOrClpRead_mode==2&&argc==2) {
		    // only if ends .mps
		    char * find = (char *)strstr(fileName.c_str(),".mps");
		    if (find&&find[4]=='\0') {
		      find[1]='p'; find[2]='a';find[3]='r';
		      FILE *fp=fopen(fileName.c_str(),"r");
		      if (fp) {
			CbcOrClpReadCommand=fp; // Read from that file
			CbcOrClpRead_mode=-1;
		      }
		    }
		  }
		} else {
		  // errors
		  std::cout<<"There were "<<status<<
		    " errors on input"<<std::endl;
		}
	      }
	    }
	    break;
	  case MODELIN:
#ifdef COIN_HAS_LINK
	    {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      bool canOpen=false;
	      if (field=="-") {
		// stdin
		canOpen=true;
		fileName = "-";
	      } else {
                bool absolutePath;
                if (dirsep=='/') {
                  // non Windows (or cygwin)
                  absolutePath=(field[0]=='/');
                } else {
                  //Windows (non cycgwin)
                  absolutePath=(field[0]=='\\');
                  // but allow for :
                  if (strchr(field.c_str(),':'))
                    absolutePath=true;
                }
		if (absolutePath) {
		  fileName = field;
		} else if (field[0]=='~') {
		  char * environVar = getenv("HOME");
		  if (environVar) {
		    std::string home(environVar);
		    field=field.erase(0,1);
		    fileName = home+field;
		  } else {
		    fileName=field;
		  }
		} else {
		  fileName = directory+field;
		}
		FILE *fp=fopen(fileName.c_str(),"r");
		if (fp) {
		  // can open - lets go for it
		  fclose(fp);
		  canOpen=true;
		} else {
		  std::cout<<"Unable to open file "<<fileName<<std::endl;
		}
	      }
	      if (canOpen) {
		CoinModel coinModel(fileName.c_str(),2);
		// load from coin model
		OsiSolverLink solver1;
		OsiSolverInterface * solver2 = solver1.clone();
		model.assignSolver(solver2,false);
		OsiSolverLink * si =
		  dynamic_cast<OsiSolverLink *>(model.solver()) ;
		assert (si != NULL);
		si->setDefaultMeshSize(0.001);
		// need some relative granularity
		si->setDefaultBound(100.0);
		si->setDefaultMeshSize(0.001);
		si->setDefaultBound(100.0);
		si->setIntegerPriority(1000);
		si->setBiLinearPriority(10000);
		CoinModel * model2 = (CoinModel *) &coinModel;
		si->load(*model2);
		// redo
		solver = model.solver();
		clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
		lpSolver = clpSolver->getModelPtr();
		clpSolver->messageHandler()->setLogLevel(0) ;
		testOsiParameters=0;
		complicatedInteger=2;
	      }
	    }
#endif
	    break;
	  case EXPORT:
	    if (goodModel) {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      bool canOpen=false;
	      if (field[0]=='/'||field[0]=='\\') {
		fileName = field;
	      } else if (field[0]=='~') {
		char * environVar = getenv("HOME");
		if (environVar) {
		  std::string home(environVar);
		  field=field.erase(0,1);
		  fileName = home+field;
		} else {
		  fileName=field;
		}
	      } else {
		fileName = directory+field;
	      }
	      FILE *fp=fopen(fileName.c_str(),"w");
	      if (fp) {
		// can open - lets go for it
		fclose(fp);
		canOpen=true;
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	      if (canOpen) {
		// If presolve on then save presolved
		bool deleteModel2=false;
		ClpSimplex * model2 = lpSolver;
		if (dualize) {
		  model2 = ((ClpSimplexOther *) model2)->dualOfModel();
		  sprintf(generalPrint,"Dual of model has %d rows and %d columns\n",
			 model2->numberRows(),model2->numberColumns());
		  generalMessageHandler->message(CLP_GENERAL,generalMessages)
		    << generalPrint
		    <<CoinMessageEol;
		  model2->setOptimizationDirection(1.0);
		}
#ifdef COIN_HAS_ASL
		if (info.numberSos&&doSOS&&usingAmpl) {
		  // SOS
		  numberSOS = info.numberSos;
		  sosStart = info.sosStart;
		  sosIndices = info.sosIndices;
		  sosReference = info.sosReference;
		  preSolve=false;
  		  clpSolver->setSOSData(numberSOS,info.sosType,sosStart,sosIndices,sosReference);
		}
#endif
		if (preSolve) {
		  ClpPresolve pinfo;
		  int presolveOptions2 = presolveOptions&~0x40000000;
		  if ((presolveOptions2&0xffff)!=0)
		    pinfo.setPresolveActions(presolveOptions2);
		  if ((printOptions&1)!=0)
		    pinfo.statistics();
                  double presolveTolerance = 
                    parameters[whichParam(PRESOLVETOLERANCE,numberParameters,parameters)].doubleValue();
                  model2 = 
		    pinfo.presolvedModel(*lpSolver,presolveTolerance,
					 true,preSolve);
		  if (model2) {
		    printf("Saving presolved model on %s\n",
			   fileName.c_str());
		    deleteModel2=true;
		  } else {
		    printf("Presolved model looks infeasible - saving original on %s\n",
			   fileName.c_str());
		    deleteModel2=false;
		    model2 = lpSolver;

		  }
		  model2->writeMps(fileName.c_str(),(outputFormat-1)/2,1+((outputFormat-1)&1));
		  if (deleteModel2)
		    delete model2;
		} else {
		  printf("Saving model on %s\n",
			   fileName.c_str());
		  if (numberSOS) {
		    // Convert names
		    int iRow;
		    int numberRows=model2->numberRows();
		    int iColumn;
		    int numberColumns=model2->numberColumns();
		    
		    char ** rowNames = NULL;
		    char ** columnNames = NULL;
		    if (model2->lengthNames()) {
		      rowNames = new char * [numberRows];
		      for (iRow=0;iRow<numberRows;iRow++) {
			rowNames[iRow] = 
			  strdup(model2->rowName(iRow).c_str());
		      }
		      
		      columnNames = new char * [numberColumns];
		      for (iColumn=0;iColumn<numberColumns;iColumn++) {
			columnNames[iColumn] = 
			  strdup(model2->columnName(iColumn).c_str());
		      }
		    }
		    clpSolver->writeMpsNative(fileName.c_str(),const_cast<const char **> (rowNames),const_cast<const char **> (columnNames),
					      (outputFormat-1)/2,1+((outputFormat-1)&1));
		    if (rowNames) {
		      for (iRow=0;iRow<numberRows;iRow++) {
			free(rowNames[iRow]);
		      }
		      delete [] rowNames;
		      for (iColumn=0;iColumn<numberColumns;iColumn++) {
			free(columnNames[iColumn]);
		      }
		      delete [] columnNames;
		    }
		  } else {
#ifdef COIN_HAS_LINK
		    OsiSolverLink * linkSolver = dynamic_cast< OsiSolverLink*> (clpSolver);
		    if (!linkSolver||!linkSolver->quadraticModel()) 
		      model2->writeMps(fileName.c_str(),(outputFormat-1)/2,1+((outputFormat-1)&1));
		    else
		      linkSolver->quadraticModel()->writeMps(fileName.c_str(),(outputFormat-1)/2,1+((outputFormat-1)&1));
#endif
		  }
		}
		time2 = CoinCpuTime();
		totalTime += time2-time1;
		time1=time2;
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case BASISIN:
	    if (goodModel) {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      bool canOpen=false;
	      if (field=="-") {
		// stdin
		canOpen=true;
		fileName = "-";
	      } else {
		if (field[0]=='/'||field[0]=='\\') {
		  fileName = field;
		} else if (field[0]=='~') {
		  char * environVar = getenv("HOME");
		  if (environVar) {
		    std::string home(environVar);
		    field=field.erase(0,1);
		    fileName = home+field;
		  } else {
		    fileName=field;
		  }
		} else {
		  fileName = directory+field;
		}
		FILE *fp=fopen(fileName.c_str(),"r");
		if (fp) {
		  // can open - lets go for it
		  fclose(fp);
		  canOpen=true;
		} else {
		  std::cout<<"Unable to open file "<<fileName<<std::endl;
		}
	      }
	      if (canOpen) {
		int values = lpSolver->readBasis(fileName.c_str());
		if (values==0)
		  basisHasValues=-1;
		else
		  basisHasValues=1;
		assert (lpSolver==clpSolver->getModelPtr());
		clpSolver->setWarmStart(NULL);
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case PRIORITYIN:
	    if (goodModel) {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
              if (field[0]=='/'||field[0]=='\\') {
                fileName = field;
              } else if (field[0]=='~') {
                char * environVar = getenv("HOME");
                if (environVar) {
                  std::string home(environVar);
                  field=field.erase(0,1);
                  fileName = home+field;
                } else {
                  fileName=field;
                }
              } else {
                fileName = directory+field;
              }
              FILE *fp=fopen(fileName.c_str(),"r");
              if (fp) {
                // can open - lets go for it
                std::string headings[]={"name","number","direction","priority","up","down",
                                        "solution","priin"};
                int got[]={-1,-1,-1,-1,-1,-1,-1,-1};
                int order[8];
                assert(sizeof(got)==sizeof(order));
                int nAcross=0;
                char line[1000];
                int numberColumns = lpSolver->numberColumns();
                if (!fgets(line,1000,fp)) {
                  std::cout<<"Odd file "<<fileName<<std::endl;
                } else {
                  char * pos = line;
                  char * put = line;
                  while (*pos>=' '&&*pos!='\n') {
                    if (*pos!=' '&&*pos!='\t') {
                      *put=tolower(*pos);
                      put++;
                    }
                    pos++;
                  }
                  *put='\0';
                  pos=line;
                  int i;
                  bool good=true;
                  while (pos) {
                    char * comma = strchr(pos,',');
                    if (comma)
                      *comma='\0';
                    for (i=0;i<(int) (sizeof(got)/sizeof(int));i++) {
                      if (headings[i]==pos) {
                        if (got[i]<0) {
                          order[nAcross]=i;
                          got[i]=nAcross++;
                        } else {
                          // duplicate
                          good=false;
                        }
                        break;
                      }
                    }
                    if (i==(int) (sizeof(got)/sizeof(int)))
                      good=false;
                    if (comma) {
                      *comma=',';
                      pos=comma+1;
                    } else {
                      break;
                    }
                  }
                  if (got[0]<0&&got[1]<0)
                    good=false;
                  if (got[0]>=0&&got[1]>=0)
                    good=false;
                  if (got[0]>=0&&!lpSolver->lengthNames())
                    good=false;
                  if (good) {
                    char ** columnNames = new char * [numberColumns];
                    pseudoDown= (double *) malloc(numberColumns*sizeof(double));
                    pseudoUp = (double *) malloc(numberColumns*sizeof(double));
                    branchDirection = (int *) malloc(numberColumns*sizeof(int));
                    priorities= (int *) malloc(numberColumns*sizeof(int));
                    free(solutionIn);
                    solutionIn=NULL;
                    free(prioritiesIn);
                    prioritiesIn=NULL;
                    int iColumn;
                    if (got[6]>=0) {
                      solutionIn = (double *) malloc(numberColumns*sizeof(double));
                      CoinZeroN(solutionIn,numberColumns);
                    }
                    if (got[7]>=0) {
                      prioritiesIn = (int *) malloc(numberColumns*sizeof(int));
                      for (iColumn=0;iColumn<numberColumns;iColumn++) 
                        prioritiesIn[iColumn]=10000;
                    }
                    for (iColumn=0;iColumn<numberColumns;iColumn++) {
                      columnNames[iColumn] = 
                        strdup(lpSolver->columnName(iColumn).c_str());
                      pseudoDown[iColumn]=0.0;
                      pseudoUp[iColumn]=0.0;
                      branchDirection[iColumn]=0;
                      priorities[iColumn]=0;
                    }
                    int nBadPseudo=0;
                    int nBadDir=0;
                    int nBadPri=0;
                    int nBadName=0;
                    int nBadLine=0;
                    int nLine=0;
                    while (fgets(line,1000,fp)) {
                      nLine++;
                      iColumn = -1;
                      double up =0.0;
                      double down=0.0;
                      int pri=0;
                      int dir=0;
                      double solValue=COIN_DBL_MAX;
                      int priValue=1000000;
                      char * pos = line;
                      char * put = line;
                      while (*pos>=' '&&*pos!='\n') {
                        if (*pos!=' '&&*pos!='\t') {
                          *put=tolower(*pos);
                          put++;
                        }
                        pos++;
                      }
                      *put='\0';
                      pos=line;
                      for (int i=0;i<nAcross;i++) {
                        char * comma = strchr(pos,',');
                        if (comma) {
                          *comma='\0';
                        } else if (i<nAcross-1) {
                          nBadLine++;
                          break;
                        }
                        switch (order[i]) {
                          // name
                        case 0:
                          for (iColumn=0;iColumn<numberColumns;iColumn++) {
                            if (!strcmp(columnNames[iColumn],pos))
                              break;
                          }
                          if (iColumn==numberColumns)
                            iColumn=-1;
                          break;
                          // number
                        case 1:
                          iColumn = atoi(pos);
                          if (iColumn<0||iColumn>=numberColumns)
                            iColumn=-1;
                          break;
                          // direction
                        case 2:
                          if (*pos=='D')
                            dir=-1;
                          else if (*pos=='U')
                            dir=1;
                          else if (*pos=='N')
                            dir=0;
                          else if (*pos=='1'&&*(pos+1)=='\0')
                            dir=1;
                          else if (*pos=='0'&&*(pos+1)=='\0')
                            dir=0;
                          else if (*pos=='1'&&*(pos+1)=='1'&&*(pos+2)=='\0')
                            dir=-1;
                          else
                            dir=-2; // bad
                          break;
                          // priority
                        case 3:
                          pri=atoi(pos);
                          break;
                          // up
                        case 4:
                          up = atof(pos);
                          break;
                          // down
                        case 5:
                          down = atof(pos);
                          break;
                          // sol value
                        case 6:
                          solValue = atof(pos);
                          break;
                          // priority in value
                        case 7:
                          priValue = atoi(pos);
                          break;
                        }
                        if (comma) {
                          *comma=',';
                          pos=comma+1;
                        }
                      }
                      if (iColumn>=0) {
                        if (down<0.0) {
                          nBadPseudo++;
                          down=0.0;
                        }
                        if (up<0.0) {
                          nBadPseudo++;
                          up=0.0;
                        }
                        if (!up)
                          up=down;
                        if (!down)
                          down=up;
                        if (dir<-1||dir>1) {
                          nBadDir++;
                          dir=0;
                        }
                        if (pri<0) {
                          nBadPri++;
                          pri=0;
                        }
                        pseudoDown[iColumn]=down;
                        pseudoUp[iColumn]=up;
                        branchDirection[iColumn]=dir;
                        priorities[iColumn]=pri;
                        if (solValue!=COIN_DBL_MAX) {
                          assert (solutionIn);
                          solutionIn[iColumn]=solValue;
                        }
                        if (priValue!=1000000) {
                          assert (prioritiesIn);
                          prioritiesIn[iColumn]=priValue;
                        }
                      } else {
                        nBadName++;
                      }
                    }
                    if (!noPrinting) {
                      printf("%d fields and %d records",nAcross,nLine);
                      if (nBadPseudo)
                        printf(" %d bad pseudo costs",nBadPseudo);
                      if (nBadDir)
                        printf(" %d bad directions",nBadDir);
                      if (nBadPri)
                        printf(" %d bad priorities",nBadPri);
                      if (nBadName)
                        printf(" ** %d records did not match on name/sequence",nBadName);
                      printf("\n");
                    }
                    for (iColumn=0;iColumn<numberColumns;iColumn++) {
                      free(columnNames[iColumn]);
                    }
                    delete [] columnNames;
                  } else {
                    std::cout<<"Duplicate or unknown keyword - or name/number fields wrong"<<line<<std::endl;
                  }
                }
                fclose(fp);
              } else {
                std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case DEBUG:
	    if (goodModel) {
              delete [] debugValues;
              debugValues=NULL;
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
                debugFile=field;
                if (debugFile=="create"||
                    debugFile=="createAfterPre") {
                  printf("Will create a debug file so this run should be a good one\n");
                  break;
                }
	      }
	      std::string fileName;
              if (field[0]=='/'||field[0]=='\\') {
                fileName = field;
              } else if (field[0]=='~') {
                char * environVar = getenv("HOME");
                if (environVar) {
                  std::string home(environVar);
                  field=field.erase(0,1);
                  fileName = home+field;
                } else {
                  fileName=field;
                }
              } else {
                fileName = directory+field;
              }
              FILE *fp=fopen(fileName.c_str(),"rb");
              if (fp) {
                // can open - lets go for it
                int numRows;
                double obj;
                fread(&numRows,sizeof(int),1,fp);
                fread(&numberDebugValues,sizeof(int),1,fp);
                fread(&obj,sizeof(double),1,fp);
                debugValues = new double[numberDebugValues+numRows];
                fread(debugValues,sizeof(double),numRows,fp);
                fread(debugValues,sizeof(double),numRows,fp);
                fread(debugValues,sizeof(double),numberDebugValues,fp);
                printf("%d doubles read into debugValues\n",numberDebugValues);
                fclose(fp);
              } else {
                std::cout<<"Unable to open file "<<fileName<<std::endl;
              }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case PRINTMASK:
            // get next field
	    {
	      std::string name = CoinReadGetString(argc,argv);
	      if (name!="EOL") {
		parameters[iParam].setStringValue(name);
                printMask = name;
	      } else {
		parameters[iParam].printString();
	      }
	    }
	    break;
	  case BASISOUT:
	    if (goodModel) {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      bool canOpen=false;
	      if (field[0]=='/'||field[0]=='\\') {
		fileName = field;
	      } else if (field[0]=='~') {
		char * environVar = getenv("HOME");
		if (environVar) {
		  std::string home(environVar);
		  field=field.erase(0,1);
		  fileName = home+field;
		} else {
		  fileName=field;
		}
	      } else {
		fileName = directory+field;
	      }
	      FILE *fp=fopen(fileName.c_str(),"w");
	      if (fp) {
		// can open - lets go for it
		fclose(fp);
		canOpen=true;
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	      if (canOpen) {
		ClpSimplex * model2 = lpSolver;
		model2->writeBasis(fileName.c_str(),outputFormat>1,outputFormat-2);
		time2 = CoinCpuTime();
		totalTime += time2-time1;
		time1=time2;
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case SAVE:
	    {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      bool canOpen=false;
	      if (field[0]=='/'||field[0]=='\\') {
		fileName = field;
	      } else if (field[0]=='~') {
		char * environVar = getenv("HOME");
		if (environVar) {
		  std::string home(environVar);
		  field=field.erase(0,1);
		  fileName = home+field;
		} else {
		  fileName=field;
		}
	      } else {
		fileName = directory+field;
	      }
	      FILE *fp=fopen(fileName.c_str(),"wb");
	      if (fp) {
		// can open - lets go for it
		fclose(fp);
		canOpen=true;
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	      if (canOpen) {
		int status;
		// If presolve on then save presolved
		bool deleteModel2=false;
		ClpSimplex * model2 = lpSolver;
		if (preSolve) {
		  ClpPresolve pinfo;
                  double presolveTolerance = 
                    parameters[whichParam(PRESOLVETOLERANCE,numberParameters,parameters)].doubleValue();
		  model2 = 
		    pinfo.presolvedModel(*lpSolver,presolveTolerance,
					 false,preSolve);
		  if (model2) {
		    printf("Saving presolved model on %s\n",
			   fileName.c_str());
		    deleteModel2=true;
		  } else {
		    printf("Presolved model looks infeasible - saving original on %s\n",
			   fileName.c_str());
		    deleteModel2=false;
		    model2 = lpSolver;

		  }
		} else {
		  printf("Saving model on %s\n",
			   fileName.c_str());
		}
		status =model2->saveModel(fileName.c_str());
		if (deleteModel2)
		  delete model2;
		if (!status) {
		  goodModel=true;
		  time2 = CoinCpuTime();
		  totalTime += time2-time1;
		  time1=time2;
		} else {
		  // errors
		  std::cout<<"There were errors on output"<<std::endl;
		}
	      }
	    }
	    break;
	  case RESTORE:
	    {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      bool canOpen=false;
	      if (field[0]=='/'||field[0]=='\\') {
		fileName = field;
	      } else if (field[0]=='~') {
		char * environVar = getenv("HOME");
		if (environVar) {
		  std::string home(environVar);
		  field=field.erase(0,1);
		  fileName = home+field;
		} else {
		  fileName=field;
		}
	      } else {
		fileName = directory+field;
	      }
	      FILE *fp=fopen(fileName.c_str(),"rb");
	      if (fp) {
		// can open - lets go for it
		fclose(fp);
		canOpen=true;
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	      if (canOpen) {
		int status =lpSolver->restoreModel(fileName.c_str());
		if (!status) {
		  goodModel=true;
		  time2 = CoinCpuTime();
		  totalTime += time2-time1;
		  time1=time2;
		} else {
		  // errors
		  std::cout<<"There were errors on input"<<std::endl;
		}
	      }
	    }
	    break;
	  case MAXIMIZE:
	    lpSolver->setOptimizationDirection(-1);
	    break;
	  case MINIMIZE:
	    lpSolver->setOptimizationDirection(1);
	    break;
	  case ALLSLACK:
	    lpSolver->allSlackBasis(true);
	    break;
	  case REVERSE:
	    if (goodModel) {
	      int iColumn;
	      int numberColumns=lpSolver->numberColumns();
	      double * dualColumnSolution = 
		lpSolver->dualColumnSolution();
	      ClpObjective * obj = lpSolver->objectiveAsObject();
	      assert(dynamic_cast<ClpLinearObjective *> (obj));
	      double offset;
	      double * objective = obj->gradient(NULL,NULL,offset,true);
	      for (iColumn=0;iColumn<numberColumns;iColumn++) {
		dualColumnSolution[iColumn] = dualColumnSolution[iColumn];
		objective[iColumn] = -objective[iColumn];
	      }
	      int iRow;
	      int numberRows=lpSolver->numberRows();
	      double * dualRowSolution = 
		lpSolver->dualRowSolution();
	      for (iRow=0;iRow<numberRows;iRow++) 
		dualRowSolution[iRow] = dualRowSolution[iRow];
	    }
	    break;
	  case DIRECTORY:
	    {
	      std::string name = CoinReadGetString(argc,argv);
	      if (name!="EOL") {
		int length=name.length();
		if (name[length-1]=='/'||name[length-1]=='\\')
		  directory=name;
		else
		  directory = name+"/";
		parameters[iParam].setStringValue(directory);
	      } else {
		parameters[iParam].printString();
	      }
	    }
	    break;
	  case STDIN:
	    CbcOrClpRead_mode=-1;
	    break;
	  case NETLIB_DUAL:
	  case NETLIB_EITHER:
	  case NETLIB_BARRIER:
	  case NETLIB_PRIMAL:
	  case NETLIB_TUNE:
	    {
	      printf("unit test is now only from clp - does same thing\n");
	      //return(22);
	    }
	    break;
	  case UNITTEST:
	    {
	      printf("unit test is now only from clp - does same thing\n");
	      //return(22);
	    }
	    break;
	  case FAKEBOUND:
	    if (goodModel) {
	      // get bound
	      double value = CoinReadGetDoubleField(argc,argv,&valid);
	      if (!valid) {
		std::cout<<"Setting "<<parameters[iParam].name()<<
		  " to DEBUG "<<value<<std::endl;
		int iRow;
		int numberRows=lpSolver->numberRows();
		double * rowLower = lpSolver->rowLower();
		double * rowUpper = lpSolver->rowUpper();
		for (iRow=0;iRow<numberRows;iRow++) {
		  // leave free ones for now
		  if (rowLower[iRow]>-1.0e20||rowUpper[iRow]<1.0e20) {
		    rowLower[iRow]=CoinMax(rowLower[iRow],-value);
		    rowUpper[iRow]=CoinMin(rowUpper[iRow],value);
		  }
		}
		int iColumn;
		int numberColumns=lpSolver->numberColumns();
		double * columnLower = lpSolver->columnLower();
		double * columnUpper = lpSolver->columnUpper();
		for (iColumn=0;iColumn<numberColumns;iColumn++) {
		  // leave free ones for now
		  if (columnLower[iColumn]>-1.0e20||
		      columnUpper[iColumn]<1.0e20) {
		    columnLower[iColumn]=CoinMax(columnLower[iColumn],-value);
		    columnUpper[iColumn]=CoinMin(columnUpper[iColumn],value);
		  }
		}
	      } else if (valid==1) {
		abort();
	      } else {
		std::cout<<"enter value for "<<parameters[iParam].name()<<
		  std::endl;
	      }
	    }
	    break;
	  case REALLY_SCALE:
	    if (goodModel) {
	      ClpSimplex newModel(*lpSolver,
				  lpSolver->scalingFlag());
	      printf("model really really scaled\n");
	      *lpSolver=newModel;
	    }
	    break;
	  case USERCLP:
#ifdef USER_HAS_FAKE_CLP
            // Replace the sample code by whatever you want
	    if (goodModel) {
              // Way of using an existing piece of code
              OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (model.solver());
              ClpSimplex * lpSolver = clpSolver->getModelPtr();
              // set time from integer model
              double timeToGo = model.getMaximumSeconds();
              lpSolver->setMaximumSeconds(timeToGo);
	      int extra1 = parameters[whichParam(EXTRA1,numberParameters,parameters)].intValue();
              fakeMain2(*lpSolver,*clpSolver,extra1);
              lpSolver = clpSolver->getModelPtr();
#ifdef COIN_HAS_ASL
	      // My actual usage has objective only in clpSolver
	      //double objectiveValue=clpSolver->getObjValue();
	      //int iStat = lpSolver->status();
	      //int iStat2 = lpSolver->secondaryStatus();
#endif
	    }
#endif
	    break;
	  case USERCBC:
#ifdef USER_HAS_FAKE_CBC
            // Replace the sample code by whatever you want
	    if (goodModel) {
              // Way of using an existing piece of code
              OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (model.solver());
              ClpSimplex * lpSolver = clpSolver->getModelPtr();
              // set time from integer model
              double timeToGo = model.getMaximumSeconds();
              lpSolver->setMaximumSeconds(timeToGo);
              fakeMain(*lpSolver,*clpSolver,model);
#ifdef COIN_HAS_ASL
	      // My actual usage has objective only in clpSolver
	      double objectiveValue=clpSolver->getObjValue();
	      int iStat = lpSolver->status();
	      int iStat2 = lpSolver->secondaryStatus();
#endif
	      // make sure solution back in correct place
	      clpSolver = dynamic_cast< OsiClpSolverInterface*> (model.solver());
	      lpSolver = clpSolver->getModelPtr();
#ifdef COIN_HAS_ASL
	      if (usingAmpl) {
		int n = clpSolver->getNumCols();
		double value = objectiveValue*lpSolver->getObjSense();
		char buf[300];
		int pos=0;
                std::string minor[]={"","","gap","nodes","time","","solutions","user ctrl-c"};
		if (iStat==0) {
		  if (objectiveValue<1.0e40) {
		    pos += sprintf(buf+pos,"optimal," );
		  } else {
		    // infeasible
		    iStat=1;
		    pos += sprintf(buf+pos,"infeasible,");
		  }
		} else if (iStat==1) {
		  if (iStat2!=6)
		    iStat=3;
		  else
		    iStat=4;
		  pos += sprintf(buf+pos,"stopped on %s,",minor[iStat2].c_str());
		} else if (iStat==2) {
		  iStat = 7;
		  pos += sprintf(buf+pos,"stopped on difficulties,");
		} else if (iStat==5) {
		  iStat = 3;
		  pos += sprintf(buf+pos,"stopped on ctrl-c,");
		} else {
		  pos += sprintf(buf+pos,"status unknown,");
		  iStat=6;
		}
		info.problemStatus=iStat;
		info.objValue = value;
		if (objectiveValue<1.0e40) 
		  pos += sprintf(buf+pos," objective %.*g",ampl_obj_prec(),
				 value);
		sprintf(buf+pos,"\n%d nodes, %d iterations",
			model.getNodeCount(),
			model.getIterationCount());
		if (objectiveValue<1.0e50) {
		  free(info.primalSolution);
		  info.primalSolution = (double *) malloc(n*sizeof(double));
		  CoinCopyN(lpSolver->primalColumnSolution(),n,info.primalSolution);
		  int numberRows = lpSolver->numberRows();
		  free(info.dualSolution);
		  info.dualSolution = (double *) malloc(numberRows*sizeof(double));
		  CoinCopyN(lpSolver->dualRowSolution(),numberRows,info.dualSolution);
		} else {
		  info.primalSolution=NULL;
		  info.dualSolution=NULL;
		}
		// put buffer into info
		strcpy(info.buffer,buf);
              }
#endif
	    }
#endif
	    break;
	  case HELP:
	    std::cout<<"Coin Solver version "<<CBCVERSION
		     <<", build "<<__DATE__<<std::endl;
	    std::cout<<"Non default values:-"<<std::endl;
	    std::cout<<"Perturbation "<<lpSolver->perturbation()<<" (default 100)"
		     <<std::endl;
	    CoinReadPrintit(
		    "Presolve being done with 5 passes\n\
Dual steepest edge steep/partial on matrix shape and factorization density\n\
Clpnnnn taken out of messages\n\
If Factorization frequency default then done on size of matrix\n\n\
(-)unitTest, (-)netlib or (-)netlibp will do standard tests\n\n\
You can switch to interactive mode at any time so\n\
clp watson.mps -scaling off -primalsimplex\nis the same as\n\
clp watson.mps -\nscaling off\nprimalsimplex"
		    );
  	    break;
	  case SOLUTION:
	    if (goodModel) {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      FILE *fp=NULL;
	      if (field=="-"||field=="EOL"||field=="stdout") {
		// stdout
		fp=stdout;
	      } else if (field=="stderr") {
		// stderr
		fp=stderr;
	      } else {
		if (field[0]=='/'||field[0]=='\\') {
		  fileName = field;
		} else if (field[0]=='~') {
		  char * environVar = getenv("HOME");
		  if (environVar) {
		    std::string home(environVar);
		    field=field.erase(0,1);
		    fileName = home+field;
		  } else {
		    fileName=field;
		  }
		} else {
		  fileName = directory+field;
		}
		fp=fopen(fileName.c_str(),"w");
	      }
	      if (fp) {
		// make fancy later on
		int iRow;
		int numberRows=lpSolver->numberRows();
		double * dualRowSolution = lpSolver->dualRowSolution();
		double * primalRowSolution = 
		  lpSolver->primalRowSolution();
		double * rowLower = lpSolver->rowLower();
		double * rowUpper = lpSolver->rowUpper();
		double primalTolerance = lpSolver->primalTolerance();
		char format[6];
		sprintf(format,"%%-%ds",CoinMax(lengthName,8));
                bool doMask = (printMask!=""&&lengthName);
		int * maskStarts=NULL;
		int maxMasks=0;
		char ** masks =NULL;
		if (doMask) {
		  int nAst =0;
		  const char * pMask2 = printMask.c_str();
		  char pMask[100];
		  int iChar;
		  int lengthMask = strlen(pMask2);
		  assert (lengthMask<100);
		  if (*pMask2=='"') {
		    if (pMask2[lengthMask-1]!='"') {
		      printf("mismatched \" in mask %s\n",pMask2);
		      break;
		    } else {
		      strcpy(pMask,pMask2+1);
		      *strchr(pMask,'"')='\0';
		    }
		  } else if (*pMask2=='\'') {
		    if (pMask2[lengthMask-1]!='\'') {
		      printf("mismatched ' in mask %s\n",pMask2);
		      break;
		    } else {
		      strcpy(pMask,pMask2+1);
		      *strchr(pMask,'\'')='\0';
		    }
		  } else {
		    strcpy(pMask,pMask2);
		  }
		  if (lengthMask>lengthName) {
		    printf("mask %s too long - skipping\n",pMask);
		    break;
		  }
		  maxMasks = 1;
		  for (iChar=0;iChar<lengthMask;iChar++) {
		    if (pMask[iChar]=='*') {
		      nAst++;
		      maxMasks *= (lengthName+1);
		    }
		  }
		  int nEntries = 1;
		  maskStarts = new int[lengthName+2];
		  masks = new char * [maxMasks];
		  char ** newMasks = new char * [maxMasks];
		  int i;
		  for (i=0;i<maxMasks;i++) {
		    masks[i] = new char[lengthName+1];
		    newMasks[i] = new char[lengthName+1];
		  }
		  strcpy(masks[0],pMask);
		  for (int iAst=0;iAst<nAst;iAst++) {
		    int nOldEntries = nEntries;
		    nEntries=0;
		    for (int iEntry = 0;iEntry<nOldEntries;iEntry++) {
		      char * oldMask = masks[iEntry];
		      char * ast = strchr(oldMask,'*');
		      assert (ast);
		      int length = strlen(oldMask)-1;
		      int nBefore = ast-oldMask;
		      int nAfter = length-nBefore;
		      // and add null
		      nAfter++;
		      for (int i=0;i<=lengthName-length;i++) {
			char * maskOut = newMasks[nEntries];
			memcpy(maskOut,oldMask,nBefore);
			for (int k=0;k<i;k++) 
			  maskOut[k+nBefore]='?';
			memcpy(maskOut+nBefore+i,ast+1,nAfter);
			nEntries++;
			assert (nEntries<=maxMasks);
		      }
		    }
		    char ** temp = masks;
		    masks = newMasks;
		    newMasks = temp;
		  }
		  // Now extend and sort
		  int * sort = new int[nEntries];
		  for (i=0;i<nEntries;i++) {
		    char * maskThis = masks[i];
		    int length = strlen(maskThis);
		    while (maskThis[length-1]==' ')
		      length--;
		    maskThis[length]='\0';
		    sort[i]=length;
		  }
		  CoinSort_2(sort,sort+nEntries,masks);
		  int lastLength=-1;
		  for (i=0;i<nEntries;i++) {
		    int length = sort[i];
		    while (length>lastLength) 
		      maskStarts[++lastLength] = i;
		  }
		  maskStarts[++lastLength]=nEntries;
		  delete [] sort;
		  for (i=0;i<maxMasks;i++)
		    delete [] newMasks[i];
		  delete [] newMasks;
		}
                if (printMode>2) {
                  for (iRow=0;iRow<numberRows;iRow++) {
                    int type=printMode-3;
                    if (primalRowSolution[iRow]>rowUpper[iRow]+primalTolerance||
                        primalRowSolution[iRow]<rowLower[iRow]-primalTolerance) {
                      fprintf(fp,"** ");
                      type=2;
                    } else if (fabs(primalRowSolution[iRow])>1.0e-8) {
                      type=1;
                    } else if (numberRows<50) {
                      type=3;
                    }
                    if (doMask&&!maskMatches(maskStarts,masks,rowNames[iRow]))
                      type =0;
                    if (type) {
                      fprintf(fp,"%7d ",iRow);
                      if (lengthName)
                        fprintf(fp,format,rowNames[iRow].c_str());
                      fprintf(fp,"%15.8g        %15.8g\n",primalRowSolution[iRow],
                              dualRowSolution[iRow]);
                    }
                  }
                }
		int iColumn;
		int numberColumns=lpSolver->numberColumns();
		double * dualColumnSolution = 
		  lpSolver->dualColumnSolution();
		double * primalColumnSolution = 
		  lpSolver->primalColumnSolution();
		double * columnLower = lpSolver->columnLower();
		double * columnUpper = lpSolver->columnUpper();
                if (printMode!=2) {
                  for (iColumn=0;iColumn<numberColumns;iColumn++) {
                    int type=(printMode>3) ? 1 :0;
                    if (primalColumnSolution[iColumn]>columnUpper[iColumn]+primalTolerance||
                        primalColumnSolution[iColumn]<columnLower[iColumn]-primalTolerance) {
                      fprintf(fp,"** ");
                      type=2;
                    } else if (fabs(primalColumnSolution[iColumn])>1.0e-8) {
                      type=1;
                    } else if (numberColumns<50) {
                      type=3;
                    }
                    // see if integer
                    if ((!lpSolver->isInteger(iColumn)||fabs(primalColumnSolution[iColumn])<1.0e-8)
                         &&printMode==1)
                      type=0;
		    if (doMask&&!maskMatches(maskStarts,masks,
					     columnNames[iColumn]))
                      type =0;
                    if (type) {
                      fprintf(fp,"%7d ",iColumn);
                      if (lengthName)
                        fprintf(fp,format,columnNames[iColumn].c_str());
                      fprintf(fp,"%15.8g        %15.8g\n",
                              primalColumnSolution[iColumn],
                              dualColumnSolution[iColumn]);
                    }
                  }
                } else {
                  // special format suitable for OsiRowCutDebugger
                  int n=0;
                  bool comma=false;
                  bool newLine=false;
                  fprintf(fp,"\tint intIndicesV[]={\n");
                  for (iColumn=0;iColumn<numberColumns;iColumn++) {
                    if(primalColumnSolution[iColumn]>0.5&&model.solver()->isInteger(iColumn)) {
                      if (comma)
                        fprintf(fp,",");
                      if (newLine)
                        fprintf(fp,"\n");
                      fprintf(fp,"%d ",iColumn);
                      comma=true;
                      newLine=false;
                      n++;
                      if (n==10) {
                        n=0;
                        newLine=true;
                      }
                    }
                  }
                  fprintf(fp,"};\n");
                  n=0;
                  comma=false;
                  newLine=false;
                  fprintf(fp,"\tdouble intSolnV[]={\n");
                  for ( iColumn=0;iColumn<numberColumns;iColumn++) {
                    if(primalColumnSolution[iColumn]>0.5&&model.solver()->isInteger(iColumn)) {
                      if (comma)
                        fprintf(fp,",");
                      if (newLine)
                        fprintf(fp,"\n");
                      int value = (int) (primalColumnSolution[iColumn]+0.5);
                      fprintf(fp,"%d. ",value);
                      comma=true;
                      newLine=false;
                      n++;
                      if (n==10) {
                        n=0;
                        newLine=true;
                      }
                    }
                  }
                  fprintf(fp,"};\n");
                }
		if (fp!=stdout)
		  fclose(fp);
		if (masks) {
		  delete [] maskStarts;
		  for (int i=0;i<maxMasks;i++)
		    delete [] masks[i];
		  delete [] masks;
		}
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	      
	    }
	    break;
	  case SAVESOL:
	    if (goodModel) {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
              if (field[0]=='/'||field[0]=='\\') {
                fileName = field;
              } else if (field[0]=='~') {
                char * environVar = getenv("HOME");
                if (environVar) {
                  std::string home(environVar);
                  field=field.erase(0,1);
                  fileName = home+field;
                } else {
                  fileName=field;
                }
              } else {
                fileName = directory+field;
              }
              saveSolution(lpSolver,fileName);
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	      
	    }
	    break;
          case DUMMY:
            break;
	  default:
	    abort();
	  }
	} 
      } else if (!numberMatches) {
	std::cout<<"No match for "<<field<<" - ? for list of commands"
		 <<std::endl;
      } else if (numberMatches==1) {
	if (!numberQuery) {
	  std::cout<<"Short match for "<<field<<" - completion: ";
	  std::cout<<parameters[firstMatch].matchName()<<std::endl;
	} else if (numberQuery) {
	  std::cout<<parameters[firstMatch].matchName()<<" : ";
	  std::cout<<parameters[firstMatch].shortHelp()<<std::endl;
	  if (numberQuery>=2) 
	    parameters[firstMatch].printLongHelp();
	}
      } else {
	if (!numberQuery) 
	  std::cout<<"Multiple matches for "<<field<<" - possible completions:"
		   <<std::endl;
	else
	  std::cout<<"Completions of "<<field<<":"<<std::endl;
	for ( iParam=0; iParam<numberParameters; iParam++ ) {
	  int match = parameters[iParam].matches(field);
	  if (match&&parameters[iParam].displayThis()) {
	    std::cout<<parameters[iParam].matchName();
	    if (numberQuery>=2) 
	      std::cout<<" : "<<parameters[iParam].shortHelp();
	    std::cout<<std::endl;
	  }
	}
      }
    }
  }
  // By now all memory should be freed
#ifdef DMALLOC
  dmalloc_log_unfreed();
  dmalloc_shutdown();
#endif
  if (babModel)
    model.moveInfo(*babModel);
  delete babModel;
  model.solver()->setWarmStart(NULL);
  return 0;
}    
static void breakdown(const char * name, int numberLook, const double * region)
{
  double range[] = {
    -COIN_DBL_MAX,
    -1.0e15,-1.0e11,-1.0e8,-1.0e5,-1.0e4,-1.0e3,-1.0e2,-1.0e1,
    -1.0,
    -1.0e-1,-1.0e-2,-1.0e-3,-1.0e-4,-1.0e-5,-1.0e-8,-1.0e-11,-1.0e-15,
    0.0,
    1.0e-15,1.0e-11,1.0e-8,1.0e-5,1.0e-4,1.0e-3,1.0e-2,1.0e-1,
    1.0,
    1.0e1,1.0e2,1.0e3,1.0e4,1.0e5,1.0e8,1.0e11,1.0e15,
    COIN_DBL_MAX};
  int nRanges = (int) (sizeof(range)/sizeof(double));
  int * number = new int[nRanges];
  memset(number,0,nRanges*sizeof(int));
  int * numberExact = new int[nRanges];
  memset(numberExact,0,nRanges*sizeof(int));
  int i;
  for ( i=0;i<numberLook;i++) {
    double value = region[i];
    for (int j=0;j<nRanges;j++) {
      if (value==range[j]) {
        numberExact[j]++;
        break;
      } else if (value<range[j]) {
        number[j]++;
        break;
      }
    }
  }
  printf("\n%s has %d entries\n",name,numberLook);
  for (i=0;i<nRanges;i++) {
    if (number[i]) 
      printf("%d between %g and %g",number[i],range[i-1],range[i]);
    if (numberExact[i]) {
      if (number[i])
        printf(", ");
      printf("%d exactly at %g",numberExact[i],range[i]);
    }
    if (number[i]+numberExact[i])
      printf("\n");
  }
  delete [] number;
  delete [] numberExact;
}
static void statistics(ClpSimplex * originalModel, ClpSimplex * model)
{
  int numberColumns = originalModel->numberColumns();
  const char * integerInformation  = originalModel->integerInformation(); 
  const double * columnLower = originalModel->columnLower();
  const double * columnUpper = originalModel->columnUpper();
  int numberIntegers=0;
  int numberBinary=0;
  int iRow,iColumn;
  if (integerInformation) {
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (integerInformation[iColumn]) {
        if (columnUpper[iColumn]>columnLower[iColumn]) {
          numberIntegers++;
          if (columnUpper[iColumn]==0.0&&columnLower[iColumn]==1) 
            numberBinary++;
        }
      }
    }
  }
  numberColumns = model->numberColumns();
  int numberRows = model->numberRows();
  columnLower = model->columnLower();
  columnUpper = model->columnUpper();
  const double * rowLower = model->rowLower();
  const double * rowUpper = model->rowUpper();
  const double * objective = model->objective();
  CoinPackedMatrix * matrix = model->matrix();
  CoinBigIndex numberElements = matrix->getNumElements();
  const int * columnLength = matrix->getVectorLengths();
  //const CoinBigIndex * columnStart = matrix->getVectorStarts();
  const double * elementByColumn = matrix->getElements();
  int * number = new int[numberRows+1];
  memset(number,0,(numberRows+1)*sizeof(int));
  int numberObjSingletons=0;
  /* cType
     0 0/inf, 1 0/up, 2 lo/inf, 3 lo/up, 4 free, 5 fix, 6 -inf/0, 7 -inf/up,
     8 0/1
  */ 
  int cType[9];
  std::string cName[]={"0.0->inf,","0.0->up,","lo->inf,","lo->up,","free,","fixed,","-inf->0.0,",
                       "-inf->up,","0.0->1.0"};
  int nObjective=0;
  memset(cType,0,sizeof(cType));
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    int length=columnLength[iColumn];
    if (length==1&&objective[iColumn])
      numberObjSingletons++;
    number[length]++;
    if (objective[iColumn])
      nObjective++;
    if (columnLower[iColumn]>-1.0e20) {
      if (columnLower[iColumn]==0.0) {
        if (columnUpper[iColumn]>1.0e20)
          cType[0]++;
        else if (columnUpper[iColumn]==1.0)
          cType[8]++;
        else if (columnUpper[iColumn]==0.0)
          cType[5]++;
        else
          cType[1]++;
      } else {
        if (columnUpper[iColumn]>1.0e20) 
          cType[2]++;
        else if (columnUpper[iColumn]==columnLower[iColumn])
          cType[5]++;
        else
          cType[3]++;
      }
    } else {
      if (columnUpper[iColumn]>1.0e20) 
        cType[4]++;
      else if (columnUpper[iColumn]==0.0) 
        cType[6]++;
      else
        cType[7]++;
    }
  }
  /* rType
     0 E 0, 1 E 1, 2 E -1, 3 E other, 4 G 0, 5 G 1, 6 G other, 
     7 L 0,  8 L 1, 9 L other, 10 Range 0/1, 11 Range other, 12 free 
  */ 
  int rType[13];
  std::string rName[]={"E 0.0,","E 1.0,","E -1.0,","E other,","G 0.0,","G 1.0,","G other,",
                       "L 0.0,","L 1.0,","L other,","Range 0.0->1.0,","Range other,","Free"};
  memset(rType,0,sizeof(rType));
  for (iRow=0;iRow<numberRows;iRow++) {
    if (rowLower[iRow]>-1.0e20) {
      if (rowLower[iRow]==0.0) {
        if (rowUpper[iRow]>1.0e20)
          rType[4]++;
        else if (rowUpper[iRow]==1.0)
          rType[10]++;
        else if (rowUpper[iRow]==0.0)
          rType[0]++;
        else
          rType[11]++;
      } else if (rowLower[iRow]==1.0) {
        if (rowUpper[iRow]>1.0e20) 
          rType[5]++;
        else if (rowUpper[iRow]==rowLower[iRow])
          rType[1]++;
        else
          rType[11]++;
      } else if (rowLower[iRow]==-1.0) {
        if (rowUpper[iRow]>1.0e20) 
          rType[6]++;
        else if (rowUpper[iRow]==rowLower[iRow])
          rType[2]++;
        else
          rType[11]++;
      } else {
        if (rowUpper[iRow]>1.0e20) 
          rType[6]++;
        else if (rowUpper[iRow]==rowLower[iRow])
          rType[3]++;
        else
          rType[11]++;
      }
    } else {
      if (rowUpper[iRow]>1.0e20) 
        rType[12]++;
      else if (rowUpper[iRow]==0.0) 
        rType[7]++;
      else if (rowUpper[iRow]==1.0) 
        rType[8]++;
      else
        rType[9]++;
    }
  }
  // Basic statistics
  printf("\n\nProblem has %d rows, %d columns (%d with objective) and %d elements\n",
         numberRows,numberColumns,nObjective,numberElements);
  if (number[0]+number[1]) {
    printf("There are ");
    if (numberObjSingletons)
      printf("%d singletons with objective ",numberObjSingletons);
    int numberNoObj = number[1]-numberObjSingletons;
    if (numberNoObj)
      printf("%d singletons with no objective ",numberNoObj);
    if (number[0])
      printf("** %d columns have no entries",number[0]);
    printf("\n");
  }
  printf("Column breakdown:\n");
  int k;
  for (k=0;k<(int) (sizeof(cType)/sizeof(int));k++) {
    printf("%d of type %s ",cType[k],cName[k].c_str());
    if (((k+1)%3)==0)
      printf("\n");
  }
  if ((k%3)!=0)
    printf("\n");
  printf("Row breakdown:\n");
  for (k=0;k<(int) (sizeof(rType)/sizeof(int));k++) {
    printf("%d of type %s ",rType[k],rName[k].c_str());
    if (((k+1)%3)==0)
      printf("\n");
  }
  if ((k%3)!=0)
    printf("\n");
  if (model->logLevel()<2)
    return ;
  int kMax = model->logLevel()>3 ? 1000000 : 10;
  k=0;
  for (iRow=1;iRow<=numberRows;iRow++) {
    if (number[iRow]) {
      k++;
      printf("%d columns have %d entries\n",number[iRow],iRow);
      if (k==kMax)
        break;
    }
  }
  if (k<numberRows) {
    int kk=k;
    k=0;
    for (iRow=numberRows;iRow>=1;iRow--) {
      if (number[iRow]) {
        k++;
        if (k==kMax)
          break;
      }
    }
    if (k>kk) {
      printf("\n    .........\n\n");
      iRow=k;
      k=0;
      for (;iRow<numberRows;iRow++) {
        if (number[iRow]) {
          k++;
          printf("%d columns have %d entries\n",number[iRow],iRow);
          if (k==kMax)
            break;
        }
      }
    }
  }
  delete [] number;
  printf("\n\n");
  // get row copy
  CoinPackedMatrix rowCopy = *matrix;
  rowCopy.reverseOrdering();
  //const int * column = rowCopy.getIndices();
  const int * rowLength = rowCopy.getVectorLengths();
  //const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
  //const double * element = rowCopy.getElements();
  number = new int[numberColumns+1];
  memset(number,0,(numberColumns+1)*sizeof(int));
  for (iRow=0;iRow<numberRows;iRow++) {
    int length=rowLength[iRow];
    number[length]++;
  }
  if (number[0])
    printf("** %d rows have no entries\n",number[0]);
  k=0;
  for (iColumn=1;iColumn<=numberColumns;iColumn++) {
    if (number[iColumn]) {
      k++;
      printf("%d rows have %d entries\n",number[iColumn],iColumn);
      if (k==kMax)
        break;
    }
  }
  if (k<numberColumns) {
    int kk=k;
    k=0;
    for (iColumn=numberColumns;iColumn>=1;iColumn--) {
      if (number[iColumn]) {
        k++;
        if (k==kMax)
          break;
      }
    }
    if (k>kk) {
      printf("\n    .........\n\n");
      iColumn=k;
      k=0;
      for (;iColumn<numberColumns;iColumn++) {
        if (number[iColumn]) {
          k++;
          printf("%d rows have %d entries\n",number[iColumn],iColumn);
          if (k==kMax)
            break;
        }
      }
    }
  }
  delete [] number;
  // Now do breakdown of ranges
  breakdown("Elements",numberElements,elementByColumn);
  breakdown("RowLower",numberRows,rowLower);
  breakdown("RowUpper",numberRows,rowUpper);
  breakdown("ColumnLower",numberColumns,columnLower);
  breakdown("ColumnUpper",numberColumns,columnUpper);
  breakdown("Objective",numberColumns,objective);
}
static bool maskMatches(const int * starts, char ** masks,
			std::string & check)
{
  // back to char as I am old fashioned
  const char * checkC = check.c_str();
  int length = strlen(checkC);
  while (checkC[length-1]==' ')
    length--;
  for (int i=starts[length];i<starts[length+1];i++) {
    char * thisMask = masks[i];
    int k;
    for ( k=0;k<length;k++) {
      if (thisMask[k]!='?'&&thisMask[k]!=checkC[k]) 
	break;
    }
    if (k==length)
      return true;
  }
  return false;
}
static void clean(char * temp)
{
  char * put = temp;
  while (*put>=' ')
    put++;
  *put='\0';
}
static void generateCode(CbcModel * model, const char * fileName,int type,int preProcess)
{
  // options on code generation
  bool sizecode = (type&4)!=0;
  type &= 3;
  FILE * fp = fopen(fileName,"r");
  assert (fp);
  int numberLines=0;
#define MAXLINES 5000
#define MAXONELINE 200
  char line[MAXLINES][MAXONELINE];
  strcpy(line[numberLines++],"0#if defined(_MSC_VER)");
  strcpy(line[numberLines++],"0// Turn off compiler warning about long names");
  strcpy(line[numberLines++],"0#  pragma warning(disable:4786)");
  strcpy(line[numberLines++],"0#endif\n");
  strcpy(line[numberLines++],"0#include <cassert>");
  strcpy(line[numberLines++],"0#include <iomanip>");
  strcpy(line[numberLines++],"0#include \"OsiClpSolverInterface.hpp\"");
  strcpy(line[numberLines++],"0#include \"CbcModel.hpp\"");
  strcpy(line[numberLines++],"0#include \"CbcCutGenerator.hpp\"");
  strcpy(line[numberLines++],"0#include \"CbcStrategy.hpp\"");
  strcpy(line[numberLines++],"0#include \"CglPreProcess.hpp\"");
  strcpy(line[numberLines++],"0#include \"CoinTime.hpp\"");
  if (preProcess>0) 
    strcpy(line[numberLines++],"0#include \"CglProbing.hpp\""); // possibly redundant
  // To allow generated 5's to be just before branchAndBound - do rest here
  strcpy(line[numberLines++],"5  cbcModel->initialSolve();");
  strcpy(line[numberLines++],"5  if (clpModel->tightenPrimalBounds()!=0) {");
  strcpy(line[numberLines++],"5    std::cout<<\"Problem is infeasible - tightenPrimalBounds!\"<<std::endl;");
  strcpy(line[numberLines++],"5    exit(1);");
  strcpy(line[numberLines++],"5  }");
  strcpy(line[numberLines++],"5  clpModel->dual();  // clean up");
  if (sizecode) {
    // override some settings
    strcpy(line[numberLines++],"5  // compute some things using problem size");
    strcpy(line[numberLines++],"5  cbcModel->setMinimumDrop(min(5.0e-2,");
    strcpy(line[numberLines++],"5       fabs(cbcModel->getMinimizationObjValue())*1.0e-3+1.0e-4));");
    strcpy(line[numberLines++],"5  if (cbcModel->getNumCols()<500)");
    strcpy(line[numberLines++],"5    cbcModel->setMaximumCutPassesAtRoot(-100); // always do 100 if possible");
    strcpy(line[numberLines++],"5  else if (cbcModel->getNumCols()<5000)");
    strcpy(line[numberLines++],"5    cbcModel->setMaximumCutPassesAtRoot(100); // use minimum drop");
    strcpy(line[numberLines++],"5  else");
    strcpy(line[numberLines++],"5    cbcModel->setMaximumCutPassesAtRoot(20);");
    strcpy(line[numberLines++],"5  cbcModel->setMaximumCutPasses(1);");
  }
  if (preProcess<=0) {
    // no preprocessing or strategy
    if (preProcess) {
      strcpy(line[numberLines++],"5  // Preprocessing using CbcStrategy");
      strcpy(line[numberLines++],"5  CbcStrategyDefault strategy(true,5,5);");
      strcpy(line[numberLines++],"5  strategy.setupPreProcessing(1);");
      strcpy(line[numberLines++],"5  cbcModel->setStrategy(strategy);");
    }
  } else {
    int translate[]={9999,0,0,-1,2,3,-2};
    strcpy(line[numberLines++],"5  // Hand coded preprocessing");
    strcpy(line[numberLines++],"5  CglPreProcess process;");
    strcpy(line[numberLines++],"5  OsiSolverInterface * saveSolver=cbcModel->solver()->clone();");
    strcpy(line[numberLines++],"5  // Tell solver we are in Branch and Cut");
    strcpy(line[numberLines++],"5  saveSolver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo) ;");
    strcpy(line[numberLines++],"5  // Default set of cut generators");
    strcpy(line[numberLines++],"5  CglProbing generator1;");
    strcpy(line[numberLines++],"5  generator1.setUsingObjective(1);");
    strcpy(line[numberLines++],"5  generator1.setMaxPass(3);");
    strcpy(line[numberLines++],"5  generator1.setMaxProbeRoot(saveSolver->getNumCols());");
    strcpy(line[numberLines++],"5  generator1.setMaxElements(100);");
    strcpy(line[numberLines++],"5  generator1.setMaxLookRoot(50);");
    strcpy(line[numberLines++],"5  generator1.setRowCuts(3);");
    strcpy(line[numberLines++],"5  // Add in generators");
    strcpy(line[numberLines++],"5  process.addCutGenerator(&generator1);");
    strcpy(line[numberLines++],"5  process.messageHandler()->setLogLevel(cbcModel->logLevel());");
    strcpy(line[numberLines++],"5  OsiSolverInterface * solver2 = ");
    sprintf(line[numberLines++],"5    process.preProcessNonDefault(*saveSolver,%d,10);",translate[preProcess]);
    strcpy(line[numberLines++],"5  // Tell solver we are not in Branch and Cut");
    strcpy(line[numberLines++],"5  saveSolver->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo) ;");
    strcpy(line[numberLines++],"5  if (solver2)");
    strcpy(line[numberLines++],"5    solver2->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo) ;");
    strcpy(line[numberLines++],"5  if (!solver2) {");
    strcpy(line[numberLines++],"5    std::cout<<\"Pre-processing says infeasible!\"<<std::endl;");
    strcpy(line[numberLines++],"5    exit(1);");
    strcpy(line[numberLines++],"5  } else {");
    strcpy(line[numberLines++],"5    std::cout<<\"processed model has \"<<solver2->getNumRows()");
    strcpy(line[numberLines++],"5	     <<\" rows, \"<<solver2->getNumCols()");
    strcpy(line[numberLines++],"5	     <<\" columns and \"<<solver2->getNumElements()");
    strcpy(line[numberLines++],"5	     <<\" elements\"<<solver2->getNumElements()<<std::endl;");
    strcpy(line[numberLines++],"5  }");
    strcpy(line[numberLines++],"5  // we have to keep solver2 so pass clone");
    strcpy(line[numberLines++],"5  solver2 = solver2->clone();");
    strcpy(line[numberLines++],"5  cbcModel->assignSolver(solver2);");
    strcpy(line[numberLines++],"5  cbcModel->initialSolve();");
  }
  while (fgets(line[numberLines],MAXONELINE,fp)) {
    assert (numberLines<MAXLINES);
    clean(line[numberLines]);
    numberLines++;
  }
  fclose(fp);
  strcpy(line[numberLines++],"0\nint main (int argc, const char *argv[])\n{");
  strcpy(line[numberLines++],"0  OsiClpSolverInterface solver1;");
  strcpy(line[numberLines++],"0  int status=1;");
  strcpy(line[numberLines++],"0  if (argc<2)");
  strcpy(line[numberLines++],"0    std::cout<<\"Please give file name\"<<std::endl;");
  strcpy(line[numberLines++],"0  else");
  strcpy(line[numberLines++],"0    status=solver1.readMps(argv[1],\"\");");
  strcpy(line[numberLines++],"0  if (status) {");
  strcpy(line[numberLines++],"0    std::cout<<\"Bad readMps \"<<argv[1]<<std::endl;");
  strcpy(line[numberLines++],"0    exit(1);");
  strcpy(line[numberLines++],"0  }\n");
  strcpy(line[numberLines++],"0  double time1 = CoinCpuTime();");
  strcpy(line[numberLines++],"0  CbcModel model(solver1);");
  strcpy(line[numberLines++],"0  // Now do requested saves and modifications");
  strcpy(line[numberLines++],"0  CbcModel * cbcModel = & model;");
  strcpy(line[numberLines++],"0  OsiSolverInterface * osiModel = model.solver();");
  strcpy(line[numberLines++],"0  OsiClpSolverInterface * osiclpModel = dynamic_cast< OsiClpSolverInterface*> (osiModel);");
  strcpy(line[numberLines++],"0  ClpSimplex * clpModel = osiclpModel->getModelPtr();");
  // add in comments about messages
  strcpy(line[numberLines++],"3  // You can save some time by switching off message building");
  strcpy(line[numberLines++],"3  // clpModel->messagesPointer()->setDetailMessages(100,10000,(int *) NULL);");
  // add in actual solve
  strcpy(line[numberLines++],"5  cbcModel->branchAndBound();");
  strcpy(line[numberLines++],"8  std::cout<<argv[1]<<\" took \"<<CoinCpuTime()-time1<<\" seconds, \"");
  strcpy(line[numberLines++],"8	   <<cbcModel->getNodeCount()<<\" nodes with objective \"");
  strcpy(line[numberLines++],"8	   <<cbcModel->getObjValue()");
  strcpy(line[numberLines++],"8	   <<(!cbcModel->status() ? \" Finished\" : \" Not finished\")");
  strcpy(line[numberLines++],"8	   <<std::endl;");
  strcpy(line[numberLines++],"5  // For best solution");
  strcpy(line[numberLines++],"5  int numberColumns = solver1.getNumCols();");
  strcpy(line[numberLines++],"5  if (cbcModel->getMinimizationObjValue()<1.0e50) {");
  if (preProcess>0) {
    strcpy(line[numberLines++],"5    // post process");
    strcpy(line[numberLines++],"5    process.postProcess(*cbcModel->solver());");
    strcpy(line[numberLines++],"5    // Solution now back in saveSolver");
    strcpy(line[numberLines++],"5    cbcModel->assignSolver(saveSolver);");
    strcpy(line[numberLines++],"5    memcpy(cbcModel->bestSolution(),cbcModel->solver()->getColSolution(),");
    strcpy(line[numberLines++],"5	   numberColumns*sizeof(double));");
  }
  strcpy(line[numberLines++],"5    // put back in original solver");
  strcpy(line[numberLines++],"5    solver1.setColSolution(cbcModel->bestSolution());");
  strcpy(line[numberLines++],"5    const double * solution = solver1.getColSolution();");
  strcpy(line[numberLines++],"8  \n  // Now you would use solution etc etc\n");
  strcpy(line[numberLines++],"5");
  strcpy(line[numberLines++],"5    // Get names from solver1 (as OsiSolverInterface may lose)");
  strcpy(line[numberLines++],"5    std::vector<std::string> columnNames = *solver1.getModelPtr()->columnNames();");
  strcpy(line[numberLines++],"5    ");
  strcpy(line[numberLines++],"5    int iColumn;");
  strcpy(line[numberLines++],"5    std::cout<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);");
  strcpy(line[numberLines++],"5    ");
  strcpy(line[numberLines++],"5    std::cout<<\"--------------------------------------\"<<std::endl;");
  strcpy(line[numberLines++],"5    for (iColumn=0;iColumn<numberColumns;iColumn++) {");
  strcpy(line[numberLines++],"5      double value=solution[iColumn];");
  strcpy(line[numberLines++],"5      if (fabs(value)>1.0e-7&&solver1.isInteger(iColumn)) ");
  strcpy(line[numberLines++],"5	std::cout<<std::setw(6)<<iColumn<<\" \"");
  strcpy(line[numberLines++],"5                 <<columnNames[iColumn]<<\" \"");
  strcpy(line[numberLines++],"5                 <<value<<std::endl;");
  strcpy(line[numberLines++],"5    }");
  strcpy(line[numberLines++],"5    std::cout<<\"--------------------------------------\"<<std::endl;");
  strcpy(line[numberLines++],"5  ");
  strcpy(line[numberLines++],"5    std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);");
  strcpy(line[numberLines++],"5  }");
  strcpy(line[numberLines++],"8  return 0;\n}");
  fp = fopen(fileName,"w");
  assert (fp);

  int wanted[9];
  memset(wanted,0,sizeof(wanted));
  wanted[0]=wanted[3]=wanted[5]=wanted[8]=1;
  if (type>0) 
    wanted[1]=wanted[6]=1;
  if (type>1) 
    wanted[2]=wanted[4]=wanted[7]=1;
  std::string header[9]=
  { "","Save values","Redundant save of default values","Set changed values",
    "Redundant set default values","Solve","Restore values","Redundant restore values","Finish up"};
  for (int iType=0;iType<9;iType++) {
    if (!wanted[iType])
      continue;
    int n=0;
    int iLine;
    for (iLine=0;iLine<numberLines;iLine++) {
      if (line[iLine][0]=='0'+iType) {
        if (!n&&header[iType]!="")
          fprintf(fp,"\n  // %s\n\n",header[iType].c_str());
        n++;
	// skip save and clp as cloned
	if (!strstr(line[iLine],"save")||(!strstr(line[iLine],"clpMo")&&
					  !strstr(line[iLine],"_Osi")))
	  fprintf(fp,"%s\n",line[iLine]+1);
      }
    }
  }
  fclose(fp);
  printf("C++ file written to %s\n",fileName);
}
/*
  Version 1.00.00 November 16 2005.
  This is to stop me (JJF) messing about too much.
  Tuning changes should be noted here.
  The testing next version may be activated by CBC_NEXT_VERSION
  This applies to OsiClp, Clp etc
  Version 1.00.01 November 24 2005
  Added several classes for advanced users.  This can't affect code (if you don't use it)
  Made some tiny changes (for N way branching) which should not change anything.
  CbcNWay object class - for N way branching this also allows use of CbcConsequence class.
  CbcBranchAllDifferent object class - for branching on general integer variables
  to stop them having same value so branches are x >= y+1 and x <= y-1.
  Added two new Cgl classes - CglAllDifferent which does column fixing (too slowly)
  and CglStored which just has a list of cuts which can be activated.
  Modified preprocess option to SOS
  Version 1.00.02 December 9 2005
  Added use of CbcStrategy to do clean preprocessing
  Added use of referenceSolver for cleaner repetition of Cbc
  Version 1.01.00 February 2 2006
  Added first try at Ampl interface
  Version 1.04 June 2007
  Goes parallel
*/
