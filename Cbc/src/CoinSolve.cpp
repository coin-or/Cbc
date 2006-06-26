// Copyright (C) 2004, International Business Machines
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
// Same version as CBC
#define CBCVERSION "1.01.00"

#include "CoinMpsIO.hpp"

#include "ClpFactorization.hpp"
#include "CoinTime.hpp"
#include "ClpSimplex.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpSolve.hpp"
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

#include "CbcModel.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcTreeLocal.hpp"
#include "CbcCompareActual.hpp"
#include "CbcBranchActual.hpp"
#include  "CbcOrClpParam.hpp"
#include  "CbcCutGenerator.hpp"
#include  "CbcStrategy.hpp"

#include "OsiClpSolverInterface.hpp"
#ifdef COIN_HAS_ASL
#include "Cbc_ampl.h"
static bool usingAmpl=false;
#endif
static double totalTime=0.0;
static void statistics(ClpSimplex * originalModel, ClpSimplex * model);
static bool maskMatches(const int * starts, char ** masks,
			std::string & check);
static void generateCode(const char * fileName,int type,int preProcess);
//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif
// Allow for interrupts
// But is this threadsafe ? (so switched off by option)

#include "CoinSignal.hpp"
static CbcModel * currentBranchModel = NULL;

extern "C" {
   static void signal_handler(int whichSignal)
   {
      if (currentBranchModel!=NULL) 
	 currentBranchModel->setMaximumNodes(0); // stop at next node
      return;
   }
}

int mainTest (int argc, const char *argv[],int algorithm,
	      ClpSimplex empty, bool doPresolve,int switchOff);
void CbcClpUnitTest (const CbcModel & saveModel);
int CbcOrClpRead_mode=1;
FILE * CbcOrClpReadCommand=stdin;
static bool noPrinting=false;
static int * analyze(OsiClpSolverInterface * solverMod, int & numberChanged, double & increment,
                     bool changeInt)
{
  OsiSolverInterface * solver = solverMod->clone();
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
  if (numberInteger&&!noPrinting)
    printf("%d integer variables",numberInteger);
  if (changeInt) {
    if (!noPrinting) {
      if (numberChanged)
        printf(" and %d variables made integer\n",numberChanged);
      else
        printf("\n");
    }
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
    if (!noPrinting) {
      if (numberChanged)
        printf(" and %d variables could be made integer\n",numberChanged);
      else
        printf("\n");
    }
    // just get increment
    CbcModel model(*solver);
    if (noPrinting)
      model.setLogLevel(0);
    model.analyzeObjective();
    double increment2=model.getCutoffIncrement();
    if (increment2>increment) {
      if (!noPrinting)
        printf("cutoff increment increased from %g to %g\n",increment,increment2);
      increment=increment2;
    }
    delete solver;
    numberChanged=0;
    delete [] changed;
    return NULL;
  }
}
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
int main (int argc, const char *argv[])
{
  /* Note
     This is meant as a stand-alone executable to do as much of coin as possible. 
     It should only have one solver known to it.
  */
  {
    double time1 = CoinCpuTime(),time2;
    bool goodModel=false;
    CoinSighandler_t saveSignal=SIG_DFL;
    // register signal handler
    saveSignal = signal(SIGINT,signal_handler);
    // Set up all non-standard stuff
    OsiClpSolverInterface solver1;
    CbcModel model(solver1);
    CbcModel * babModel = NULL;
    model.setNumberBeforeTrust(21);
    int cutPass=-1234567;
    OsiSolverInterface * solver = model.solver();
    OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
    ClpSimplex * lpSolver = clpSolver->getModelPtr();
    clpSolver->messageHandler()->setLogLevel(0) ;
    model.messageHandler()->setLogLevel(1);
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
    int * sosPriority=NULL;
#ifdef COIN_HAS_ASL
    ampl_info info;
    if (argc>2&&!strcmp(argv[2],"-AMPL")) {
      usingAmpl=true;
      int returnCode = readAmpl(&info,argc,const_cast<char **>(argv));
      if (returnCode)
        return returnCode;
      CbcOrClpRead_mode=2; // so will start with parameters
      // see if log in list
      noPrinting=true;
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
      solver->loadProblem(info.numberColumns,info.numberRows,info.starts,
                          info.rows,info.elements,
                          info.columnLower,info.columnUpper,info.objective,
                          info.rowLower,info.rowUpper);
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
      // Set integer variables
      for (int i=info.numberColumns-info.numberBinary-info.numberIntegers;
           i<info.numberColumns;i++)
        solver->setInteger(i);
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
    int doSprint=-1;
    int doScaling=1;
    // set reasonable defaults
    int preSolve=5;
    int preProcess=4;
    bool useStrategy=false;
    bool preSolveFile=false;
   
    double djFix=1.0e100;
    double gapRatio=1.0e100;
    double tightenFactor=0.0;
    lpSolver->setPerturbation(50);
    lpSolver->messageHandler()->setPrefix(false);
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
#define CBCMAXPARAMETERS 200
    CbcOrClpParam parameters[CBCMAXPARAMETERS];
    int numberParameters ;
    establishParams(numberParameters,parameters) ;
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
    parameters[slog].setIntValue(0);
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
    model.setNumberBeforeTrust(5);
    parameters[whichParam(NUMBERBEFORE,numberParameters,parameters)].setIntValue(5);
    parameters[whichParam(MAXNODES,numberParameters,parameters)].setIntValue(model.getMaximumNodes());
    model.setNumberStrong(5);
    parameters[whichParam(STRONGBRANCHING,numberParameters,parameters)].setIntValue(model.numberStrong());
    parameters[whichParam(INFEASIBILITYWEIGHT,numberParameters,parameters)].setDoubleValue(model.getDblParam(CbcModel::CbcInfeasibilityWeight));
    parameters[whichParam(INTEGERTOLERANCE,numberParameters,parameters)].setDoubleValue(model.getDblParam(CbcModel::CbcIntegerTolerance));
    parameters[whichParam(INCREMENT,numberParameters,parameters)].setDoubleValue(model.getDblParam(CbcModel::CbcCutoffIncrement));
    // Set up likely cut generators and defaults
    parameters[whichParam(PREPROCESS,numberParameters,parameters)].setCurrentOption("sos");
    parameters[whichParam(MIPOPTIONS,numberParameters,parameters)].setIntValue(128|64|1);
    parameters[whichParam(MOREMIPOPTIONS,numberParameters,parameters)].setIntValue(-1);
    parameters[whichParam(MAXHOTITS,numberParameters,parameters)].setIntValue(100);
    parameters[whichParam(CUTSSTRATEGY,numberParameters,parameters)].setCurrentOption("on");
    parameters[whichParam(HEURISTICSTRATEGY,numberParameters,parameters)].setCurrentOption("on");
    int doSOS=1;
    int verbose=0;
    CglGomory gomoryGen;
    // try larger limit
    gomoryGen.setLimitAtRoot(512);
    gomoryGen.setLimit(50);
    // set default action (0=off,1=on,2=root)
    int gomoryAction=3;
    parameters[whichParam(GOMORYCUTS,numberParameters,parameters)].setCurrentOption("ifmove");

    CglProbing probingGen;
    probingGen.setUsingObjective(true);
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
    int probingAction=3;
    parameters[whichParam(PROBINGCUTS,numberParameters,parameters)].setCurrentOption("ifmove");

    CglKnapsackCover knapsackGen;
    //knapsackGen.switchOnExpensive();
    // set default action (0=off,1=on,2=root)
    int knapsackAction=3;
    parameters[whichParam(KNAPSACKCUTS,numberParameters,parameters)].setCurrentOption("ifmove");

    CglRedSplit redsplitGen;
    //redsplitGen.setLimit(100);
    // set default action (0=off,1=on,2=root)
    // Off as seems to give some bad cuts
    int redsplitAction=2;
    parameters[whichParam(REDSPLITCUTS,numberParameters,parameters)].setCurrentOption("root");

    CglClique cliqueGen(false,true);
    cliqueGen.setStarCliqueReport(false);
    cliqueGen.setRowCliqueReport(false);
    cliqueGen.setMinViolation(0.1);
    // set default action (0=off,1=on,2=root)
    int cliqueAction=3;
    parameters[whichParam(CLIQUECUTS,numberParameters,parameters)].setCurrentOption("ifmove");

    CglMixedIntegerRounding2 mixedGen;
    // set default action (0=off,1=on,2=root)
    int mixedAction=3;
    parameters[whichParam(MIXEDCUTS,numberParameters,parameters)].setCurrentOption("ifmove");

    CglFlowCover flowGen;
    // set default action (0=off,1=on,2=root)
    int flowAction=3;
    parameters[whichParam(FLOWCUTS,numberParameters,parameters)].setCurrentOption("ifmove");

    CglTwomir twomirGen;
    twomirGen.setMaxElements(250);
    // set default action (0=off,1=on,2=root)
    int twomirAction=2;
    parameters[whichParam(TWOMIRCUTS,numberParameters,parameters)].setCurrentOption("root");

    bool useRounding=true;
    parameters[whichParam(ROUNDING,numberParameters,parameters)].setCurrentOption("on");
    bool useFpump=true;
    parameters[whichParam(FPUMP,numberParameters,parameters)].setCurrentOption("on");
    bool useGreedy=true;
    parameters[whichParam(GREEDY,numberParameters,parameters)].setCurrentOption("on");
    bool useCombine=true;
    parameters[whichParam(COMBINE,numberParameters,parameters)].setCurrentOption("on");
    bool useLocalTree=false;
    parameters[whichParam(COSTSTRATEGY,numberParameters,parameters)].setCurrentOption("off");
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
      std::cout<<"Coin Cbc and Clp Solver version "<<CBCVERSION
               <<", build "<<__DATE__<<std::endl;
      // Print command line
      if (argc>1) {
        printf("command line - ");
        for (int i=0;i<argc;i++)
          printf("%s ",argv[i]);
        printf("\n");
      }
    }
    while (1) {
      // next command
      field=CoinReadGetCommand(argc,argv);
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
	  if (!lpSolver->integerInformation())
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
                defaultSettings=false; // user knows what she is doing
		break;
	      default:
		abort();
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
	      else if (parameters[iParam].type()==CUTPASS)
		cutPass = value;
	      else if (parameters[iParam].type()==VERBOSE)
		verbose = value;
	      else if (parameters[iParam].type()==FPUMPITS)
		{ useFpump = true;parameters[iParam].setIntValue(value);}
              parameters[iParam].setIntParameter(lpSolver,value);
	    } else {
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
	    case ROUNDING:
              defaultSettings=false; // user knows what she is doing
	      useRounding = action;
	      break;
	    case FPUMP:
              defaultSettings=false; // user knows what she is doing
              useFpump=action;
	      break;
            case CUTSSTRATEGY:
	      gomoryAction = action;
	      probingAction = action;
	      knapsackAction = action;
	      cliqueAction = action;
	      flowAction = action;
	      mixedAction = action;
	      twomirAction = action;
              parameters[whichParam(GOMORYCUTS,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(PROBINGCUTS,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(KNAPSACKCUTS,numberParameters,parameters)].setCurrentOption(action);
              if (!action) {
                redsplitAction = action;
                parameters[whichParam(REDSPLITCUTS,numberParameters,parameters)].setCurrentOption(action);
              }
              parameters[whichParam(CLIQUECUTS,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(FLOWCUTS,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(MIXEDCUTS,numberParameters,parameters)].setCurrentOption(action);
              parameters[whichParam(TWOMIRCUTS,numberParameters,parameters)].setCurrentOption(action);
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
                model2 = ((ClpSimplexOther *) model2)->dualOfModel();
                printf("Dual of model has %d rows and %d columns\n",
                       model2->numberRows(),model2->numberColumns());
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
                    printf("Doing %d presolve passes - picking up non-costed slacks\n",
                           preSolve);
                    solveOptions.setDoSingletonColumn(true);
                  } else {
                    preSolve -=100;
                    presolveType=ClpSolve::presolveNumberCost;
                    printf("Doing %d presolve passes - picking up costed slacks\n",
                           preSolve);
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
	      model2->initialSolve(solveOptions);
	      basisHasValues=1;
              if (dualize) {
                ((ClpSimplexOther *) lpSolver)->restoreFromDual(model2);
                delete model2;
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
              int nOut = outDupRow(clpSolver);
              if (nOut&&!noPrinting)
                printf("%d rows eliminated\n",nOut);
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
            // User can set options - main differenec is lack of model and CglPreProcess
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
              // Don't switch off all output
              {
                OsiSolverInterface * solver = model.solver();
                OsiClpSolverInterface * si =
                  dynamic_cast<OsiClpSolverInterface *>(solver) ;
                assert (si != NULL);
                si->setSpecialOptions(0x40000000);
              }
              if (!miplib) {
                model.initialSolve();
                OsiSolverInterface * solver = model.solver();
                OsiClpSolverInterface * si =
                  dynamic_cast<OsiClpSolverInterface *>(solver) ;
		if (si->getModelPtr()->tightenPrimalBounds()!=0) {
		  std::cout<<"Problem is infeasible - tightenPrimalBounds!"<<std::endl;
		  exit(1);
		}
		si->getModelPtr()->dual();  // clean up
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
                if (tightenFactor) {
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
                  printf("%d columns fixed\n",numberFixed);
                }
              }
              // See if we want preprocessing
              OsiSolverInterface * saveSolver=NULL;
              CglPreProcess process;
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
              if (preProcess==6) {
		// use strategy instead
		preProcess=0;
		useStrategy=true;
	      }
              if (preProcess&&type==BAB) {
                saveSolver=babModel->solver()->clone();
                /* Do not try and produce equality cliques and
                   do up to 10 passes */
                OsiSolverInterface * solver2;
                {
                  // Tell solver we are in Branch and Cut
                  saveSolver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo) ;
                  // Default set of cut generators
                  CglProbing generator1;
                  generator1.setUsingObjective(true);
                  generator1.setMaxPass(3);
                  generator1.setMaxProbeRoot(saveSolver->getNumCols());
                  generator1.setMaxElements(100);
                  generator1.setMaxLookRoot(50);
                  generator1.setRowCuts(3);
                  // Add in generators
                  process.addCutGenerator(&generator1);
                  int translate[]={9999,0,0,-1,2,3};
                  process.messageHandler()->setLogLevel(babModel->logLevel());
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
                  solver2 = process.preProcessNonDefault(*saveSolver,translate[preProcess],10);
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
                  sprintf(info.buffer,"infeasible by pre-processing");
                  info.primalSolution=NULL;
                  info.dualSolution=NULL;
                  break;
                }
#endif
                if (!noPrinting) {
                  if (!solver2) {
                    printf("Pre-processing says infeasible\n");
                    break;
                  } else {
                    printf("processed model has %d rows, %d columns and %d elements\n",
                           solver2->getNumRows(),solver2->getNumCols(),solver2->getNumElements());
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
                if (noPrinting)
                  modelC->setLogLevel(0);
                if (modelC->tightenPrimalBounds()!=0) {
                  std::cout<<"Problem is infeasible!"<<std::endl;
                  break;
                }
                modelC->dual();
              }
              if (debugValues) {
                // for debug
                std::string problemName ;
                babModel->solver()->getStrParam(OsiProbName,problemName) ;
                //babModel->solver()->activateRowCutDebugger(problemName.c_str()) ;
                twomirGen.probname_=strdup(problemName.c_str());
                // checking seems odd
                //redsplitGen.set_given_optsol(babModel->solver()->getRowCutDebuggerAlways()->optimalSolution(),
                //                         babModel->getNumCols());
              }
	      if (useCosts) {
		int numberColumns = babModel->getNumCols();
		int * sort = new int[numberColumns];
		double * dsort = new double[numberColumns];
		int * priority = new int [numberColumns];
		const double * objective = babModel->getObjCoefficients();
		int iColumn;
		int n=0;
		for (iColumn=0;iColumn<numberColumns;iColumn++) {
		  if (babModel->isInteger(iColumn)) {
		    sort[n]=n;
		    if (useCosts==1)
		      dsort[n++]=-objective[iColumn];
		    else
		      dsort[n++]=iColumn;
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
		delete [] priority;
		delete [] sort;
		delete [] dsort;
	      }
              // FPump done first as it only works if no solution
              CbcHeuristicFPump heuristic4(*babModel);
              if (useFpump) {
                heuristic4.setMaximumPasses(parameters[whichParam(FPUMPITS,numberParameters,parameters)].intValue());
                babModel->addHeuristic(&heuristic4);
              }
              if (!miplib) {
                CbcRounding heuristic1(*babModel);
                if (useRounding)
                  babModel->addHeuristic(&heuristic1) ;
                CbcHeuristicLocal heuristic2(*babModel);
                heuristic2.setSearchType(1);
                if (useCombine)
                  babModel->addHeuristic(&heuristic2);
                CbcHeuristicGreedyCover heuristic3(*babModel);
                CbcHeuristicGreedyEquality heuristic3a(*babModel);
                if (useGreedy) {
                  babModel->addHeuristic(&heuristic3);
                  babModel->addHeuristic(&heuristic3a);
                }
                if (useLocalTree) {
                  CbcTreeLocal localTree(babModel,NULL,10,0,0,10000,2000);
                  babModel->passInTreeHandler(localTree);
                }
              }
              // add cut generators if wanted
              int switches[20];
              int numberGenerators=0;
              if (probingAction==1) {
                babModel->addCutGenerator(&probingGen,-1,"Probing");
                switches[numberGenerators++]=0;
              } else if (probingAction>=2) {
                babModel->addCutGenerator(&probingGen,-101+probingAction,"Probing");
                switches[numberGenerators++]=0;
              }
              if (gomoryAction==1) {
                babModel->addCutGenerator(&gomoryGen,-1,"Gomory");
                switches[numberGenerators++]=1;
              } else if (gomoryAction>=2) {
                babModel->addCutGenerator(&gomoryGen,-101+gomoryAction,"Gomory");
                switches[numberGenerators++]=-1;
              }
              if (knapsackAction==1) {
                babModel->addCutGenerator(&knapsackGen,-1,"Knapsack");
                switches[numberGenerators++]=0;
              } else if (knapsackAction>=2) {
                babModel->addCutGenerator(&knapsackGen,-101+knapsackAction,"Knapsack");
                switches[numberGenerators++]=0;
              }
              if (redsplitAction==1) {
                babModel->addCutGenerator(&redsplitGen,-1,"Reduce-and-split");
                switches[numberGenerators++]=1;
              } else if (redsplitAction>=2) {
                babModel->addCutGenerator(&redsplitGen,-101+redsplitAction,"Reduce-and-split");
                switches[numberGenerators++]=1;
              }
              if (cliqueAction==1) {
                babModel->addCutGenerator(&cliqueGen,-1,"Clique");
                switches[numberGenerators++]=1;
              } else if (cliqueAction>=2) {
                babModel->addCutGenerator(&cliqueGen,-101+cliqueAction,"Clique");
                switches[numberGenerators++]=-1;
              }
              if (mixedAction==1) {
                babModel->addCutGenerator(&mixedGen,-1,"MixedIntegerRounding2");
                switches[numberGenerators++]=1;
              } else if (mixedAction>=2) {
                babModel->addCutGenerator(&mixedGen,-101+mixedAction,"MixedIntegerRounding2");
                switches[numberGenerators++]=-1;
              }
              if (flowAction==1) {
                babModel->addCutGenerator(&flowGen,-1,"FlowCover");
                switches[numberGenerators++]=1;
              } else if (flowAction>=2) {
                babModel->addCutGenerator(&flowGen,-101+flowAction,"FlowCover");
                switches[numberGenerators++]=1;
              }
              if (twomirAction==1) {
                babModel->addCutGenerator(&twomirGen,-1,"TwoMirCuts");
                switches[numberGenerators++]=1;
              } else if (twomirAction>=2) {
                babModel->addCutGenerator(&twomirGen,-101+twomirAction,"TwoMirCuts");
                switches[numberGenerators++]=1;
              }
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
                babModel->setMaximumCutPasses(1);
              }
              // Do more strong branching if small
              //if (babModel->getNumCols()<5000)
              //babModel->setNumberStrong(20);
              // Switch off strong branching if wanted
              //if (babModel->getNumCols()>10*babModel->getNumRows())
              //babModel->setNumberStrong(0);
              if (!noPrinting) {
                babModel->messageHandler()->setLogLevel(parameters[log].intValue());
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
              if (!miplib)
                changed=analyze( osiclp,numberChanged,increment,false);
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
              int mipOptions = parameters[whichParam(MIPOPTIONS,numberParameters,parameters)].intValue();
              if (mipOptions!=(128|64|1))
                printf("mip options %d\n",mipOptions);
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
              babModel->setSpecialOptions(2);
              currentBranchModel = babModel;
              OsiSolverInterface * strengthenedModel=NULL;
              if (type==BAB||type==MIPLIB) {
                int moreMipOptions = parameters[whichParam(MOREMIPOPTIONS,numberParameters,parameters)].intValue();
                if (moreMipOptions>=0) {
                  printf("more mip options %d\n",moreMipOptions);
                  babModel->setSearchStrategy(moreMipOptions);
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
                if (preProcess&&process.numberSOS()) {
                  int numberSOS = process.numberSOS();
                  int numberIntegers = babModel->numberIntegers();
                  /* model may not have created objects
                     If none then create
                  */
                  if (!numberIntegers||!babModel->numberObjects()) {
                    int type = (pseudoUp) ? 1 : 0;
                    babModel->findIntegers(true,type);
                    numberIntegers = babModel->numberIntegers();
                  }
                  CbcObject ** oldObjects = babModel->objects();
                  // Do sets and priorities
                  CbcObject ** objects = new CbcObject * [numberSOS];
                  // set old objects to have low priority
                  int numberOldObjects = babModel->numberObjects();
                  int numberColumns = babModel->getNumCols();
                  for (int iObj = 0;iObj<numberOldObjects;iObj++) {
                    oldObjects[iObj]->setPriority(numberColumns+1);
                    int iColumn = oldObjects[iObj]->columnNumber();
                    assert (iColumn>=0);
                    if (iColumn>=numberOriginalColumns)
                      continue;
                    if (originalColumns)
                      iColumn = originalColumns[iColumn];
                    if (branchDirection)
                      oldObjects[iObj]->setPreferredWay(branchDirection[iColumn]);
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
                    CbcObject ** objects = new CbcObject * [numberSOS];
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
                      if (nMissing)
                        printf("%d SOS variables vanished due to pre processing? - check validity?\n",nMissing);
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
                    babModel->addObjects(numberSOS,objects);
                    for (iSOS=0;iSOS<numberSOS;iSOS++)
                      delete objects[iSOS];
                    delete [] objects;
                  }
                  CbcObject ** objects = babModel->objects();
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
                    if (branchDirection)
                      objects[iObj]->setPreferredWay(branchDirection[iColumn]);
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
                  free(sosPriority);
                  sosPriority=NULL;
#ifdef COIN_HAS_ASL
                }
#endif                
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
                    generateCode("user_driver.cpp",cppValue,prepro);
                  } else {
                    std::cout<<"Unable to open file user_driver.cpp"<<std::endl;
                  }
                }
		if (useStrategy) {
		  CbcStrategyDefault strategy(true,5,5);
                  strategy.setupPreProcessing(1);
		  babModel->setStrategy(strategy);
		}
                babModel->branchAndBound(statistics);
              } else if (type==MIPLIB) {
                CbcStrategyDefault strategy(true,5,5);
                // Set up pre-processing to find sos if wanted
                if (preProcess)
                  strategy.setupPreProcessing(2);
                babModel->setStrategy(strategy);
                CbcClpUnitTest(*babModel);
                goodModel=false;
                break;
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
                std::cout<<"Cuts at root node changed objective from "<<babModel->getContinuousObjective()
                         <<" to "<<babModel->rootObjectiveAfterCuts()<<std::endl;
                
                for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
                  CbcCutGenerator * generator = babModel->cutGenerator(iGenerator);
                  std::cout<<generator->cutGeneratorName()<<" was tried "
                           <<generator->numberTimesEntered()<<" times and created "
                           <<generator->numberCutsInTotal()<<" cuts of which "
                           <<generator->numberCutsActive()<<" were active after adding rounds of cuts";
                  if (generator->timing())
                    std::cout<<" ( "<<generator->timeInCutGenerator()<<" seconds)"<<std::endl;
                  else
                    std::cout<<std::endl;
                }
              }
              time2 = CoinCpuTime();
              totalTime += time2-time1;
              // For best solution
              double * bestSolution = NULL;
              if (babModel->getMinimizationObjValue()<1.0e50&&type==BAB) {
                // post process
                if (preProcess) {
                  int n = saveSolver->getNumCols();
                  bestSolution = new double [n];
                  process.postProcess(*babModel->solver());
                  // Solution now back in saveSolver
                  babModel->assignSolver(saveSolver);
                  memcpy(bestSolution,babModel->solver()->getColSolution(),n*sizeof(double));
                } else {
                  int n = babModel->solver()->getNumCols();
                  bestSolution = new double [n];
                  memcpy(bestSolution,babModel->solver()->getColSolution(),n*sizeof(double));
                }
              }
              if (type==STRENGTHEN&&strengthenedModel)
                clpSolver = dynamic_cast< OsiClpSolverInterface*> (strengthenedModel);
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
                if (bestSolution)
                  memcpy(lpSolver->primalColumnSolution(),bestSolution,n*sizeof(double));
                if (debugFile=="create"&&bestSolution) {
                  saveSolution(lpSolver,"debug.file");
                }
                delete [] bestSolution;
                std::string statusName[]={"Finished","Stopped on ","Difficulties",
                                          "","","User ctrl-c"};
                std::string minor[]={"","","gap","nodes","time","","solutions","user ctrl-c"};
                int iStat = babModel->status();
                int iStat2 = babModel->secondaryStatus();
                if (!noPrinting)
                  std::cout<<"Result - "<<statusName[iStat]<<minor[iStat2]
                           <<" objective "<<babModel->getObjValue()<<
                    " after "<<babModel->getNodeCount()<<" nodes and "
                           <<babModel->getIterationCount()<<
                    " iterations - took "<<time2-time1<<" seconds"<<std::endl;
#ifdef COIN_HAS_ASL
                if (usingAmpl) {
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
                  if (babModel->getObjValue()<1.0e40) 
                    pos += sprintf(buf+pos," objective %.*g",ampl_obj_prec(),
                                   value);
                  sprintf(buf+pos,"\n%d nodes, %d iterations",
                          babModel->getNodeCount(),
                          babModel->getIterationCount());
                  if (bestSolution) {
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
              } else {
                std::cout<<"Model strengthened - now has "<<clpSolver->getNumRows()
                         <<" rows"<<std::endl;
              }
              time1 = time2;
              delete babModel;
              babModel=NULL;
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
                free(sosPriority);
                sosPriority=NULL;
#ifdef COIN_HAS_ASL
              }
#endif                
              delete babModel;
              babModel=NULL;
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
		int status =lpSolver->readMps(fileName.c_str(),
						   keepImportNames!=0,
						   allowImportErrors!=0);
		if (!status||(status>0&&allowImportErrors)) {
                  if (keepImportNames) {
                    lengthName = lpSolver->lengthNames();
                    rowNames = *(lpSolver->rowNames());
                    columnNames = *(lpSolver->columnNames());
                  } else {
                    lengthName=0;
                  }
		  goodModel=true;
		  //Set integers in clpsolver
		  const char * info = lpSolver->integerInformation();
		  if (info) {
		    int numberColumns = lpSolver->numberColumns();
		    int i;
		    for (i=0;i<numberColumns;i++) {
		      if (info[i]) 
			clpSolver->setInteger(i);
		    }
		  }
		  // sets to all slack (not necessary?)
		  lpSolver->createStatus();
		  time2 = CoinCpuTime();
		  totalTime += time2-time1;
		  time1=time2;
		  // Go to canned file if just input file
		  if (CbcOrClpRead_mode==2&&argc==2) {
		    // only if ends .mps
		    std::string::size_type loc = fileName.find(".mps") ;
		    if (loc != std::string::npos &&
			fileName.length() == loc+3)
		    { fileName.replace(loc+1,3,"par") ;
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
		} else {
		  printf("Saving model on %s\n",
			   fileName.c_str());
		}
#if 0
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
#ifdef STRIPBLANKS
		    char * xx = rowNames[iRow];
		    int i;
		    int length = strlen(xx);
		    int n=0;
		    for (i=0;i<length;i++) {
		      if (xx[i]!=' ')
			xx[n++]=xx[i];
		    }
		    xx[n]='\0';
#endif
		  }
		  
		  columnNames = new char * [numberColumns];
		  for (iColumn=0;iColumn<numberColumns;iColumn++) {
		    columnNames[iColumn] = 
		      strdup(model2->columnName(iColumn).c_str());
#ifdef STRIPBLANKS
		    char * xx = columnNames[iColumn];
		    int i;
		    int length = strlen(xx);
		    int n=0;
		    for (i=0;i<length;i++) {
		      if (xx[i]!=' ')
			xx[n++]=xx[i];
		    }
		    xx[n]='\0';
#endif
		  }
		}
		CoinMpsIO writer;
		writer.setMpsData(*model2->matrix(), COIN_DBL_MAX,
				  model2->getColLower(), model2->getColUpper(),
				  model2->getObjCoefficients(),
				  (const char*) 0 /*integrality*/,
				  model2->getRowLower(), model2->getRowUpper(),
				  columnNames, rowNames);
		// Pass in array saying if each variable integer
		writer.copyInIntegerInformation(model2->integerInformation());
		writer.setObjectiveOffset(model2->objectiveOffset());
		writer.writeMps(fileName.c_str(),0,1,1);
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
#else
		model2->writeMps(fileName.c_str(),(outputFormat-1)/2,1+((outputFormat-1)&1));
#endif
		if (deleteModel2)
		  delete model2;
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
                    char ** columnNames = columnNames = new char * [numberColumns];
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
	      // create fields for unitTest
	      const char * fields[4];
	      int nFields=2;
	      fields[0]="fake main from unitTest";
	      fields[1]="-netlib";
	      if (directory!=defaultDirectory) {
		fields[2]="-netlibDir";
		fields[3]=directory.c_str();
		nFields=4;
	      }
	      int algorithm;
	      if (type==NETLIB_DUAL) {
		std::cerr<<"Doing netlib with dual algorithm"<<std::endl;
		algorithm =0;
	      } else if (type==NETLIB_BARRIER) {
		std::cerr<<"Doing netlib with barrier algorithm"<<std::endl;
		algorithm =2;
	      } else if (type==NETLIB_EITHER) {
		std::cerr<<"Doing netlib with dual or primal algorithm"<<std::endl;
		algorithm =3;
	      } else if (type==NETLIB_TUNE) {
		std::cerr<<"Doing netlib with best algorithm!"<<std::endl;
		algorithm =5;
                // uncomment next to get active tuning
                // algorithm=6;
	      } else {
		std::cerr<<"Doing netlib with primal agorithm"<<std::endl;
		algorithm=1;
	      }
              int specialOptions = lpSolver->specialOptions();
              lpSolver->setSpecialOptions(0);
	      mainTest(nFields,fields,algorithm,*lpSolver,
		       (preSolve!=0),specialOptions);
	    }
	    break;
	  case UNITTEST:
	    {
	      // create fields for unitTest
	      const char * fields[3];
	      int nFields=1;
	      fields[0]="fake main from unitTest";
	      if (directory!=defaultDirectory) {
		fields[1]="-mpsDir";
		fields[2]=directory.c_str();
		nFields=3;
	      }
	      mainTest(nFields,fields,false,*lpSolver,(preSolve!=0),
		       false);
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
            // Replace the sample code by whatever you want
	    if (goodModel) {
              printf("Dummy user clp code - model has %d rows and %d columns\n",
                     lpSolver->numberRows(),lpSolver->numberColumns());
	    }
	    break;
	  case USERCBC:
            // Replace the sample code by whatever you want
	    if (goodModel) {
              printf("Dummy user cbc code - model has %d rows and %d columns\n",
                     model.getNumRows(),model.getNumCols());
  // Reduce printout
  //solver1.setHintParam(OsiDoReducePrint,true,OsiHintTry);
  OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (model.solver());
  // go faster stripes
  if (osiclp&&0) {
    // Turn this off if you get problems
    // Used to be automatically set
    osiclp->setSpecialOptions(128);
    if(osiclp->getNumRows()<300&&osiclp->getNumCols()<500) {
      //osiclp->setupForRepeatedUse(2,0);
      osiclp->setupForRepeatedUse(0,0);
    }
  } 
  // Allow rounding heuristic

  CbcRounding heuristic1(model);
  model.addHeuristic(&heuristic1);

  // Do initial solve to continuous
  ClpPrimalColumnSteepest steepest(5);
  osiclp->getModelPtr()->setPrimalColumnPivotAlgorithm(steepest);
  osiclp->getModelPtr()->setPerturbation(50);
  osiclp->getModelPtr()->setInfeasibilityCost(1.0e9);
  osiclp->setHintParam(OsiDoPresolveInInitial,false,OsiHintTry);
  osiclp->setHintParam(OsiDoDualInInitial,false,OsiHintTry);
  //osiclp->setHintParam(OsiDoDualInResolve,false,OsiHintTry);
  model.setSpecialOptions(model.specialOptions()|4);
  osiclp->getModelPtr()->defaultFactorizationFrequency();
  {
    ClpSimplex * clp = osiclp->getModelPtr();
    // fix integers to 1
    int numberColumns = clp->numberColumns();
    double * lower = clp->columnLower();
    int i;
    for (i=0;i<numberColumns;i++) {
      if (osiclp->isInteger(i))
        lower[i]=1.0;
    }
    clp->primal();
    double objValue = clp->objectiveValue();
    osiclp->setDblParam(OsiDualObjectiveLimit,objValue+1.0e-4);
    // unfix integers 
    for (i=0;i<numberColumns;i++) {
      if (osiclp->isInteger(i))
        lower[i]=0.0;
    }
    clp->primal();
    //clp->dual();
    int nArt=0;
    int nFixed0=0,nFixed1=0;
    double gap=objValue-clp->objectiveValue();
    // for those at one fix anyway
    double gap2=1.0;
    const double * solution = clp->primalColumnSolution();
    const double * dj = clp->dualColumnSolution();
    const double * objective = clp->objective();
    double * upper = clp->columnUpper();
    for (i=0;i<numberColumns;i++) {
      if (objective[i]>1.0e5&&solution[i]>1.0e-8)
        nArt++;
      if (osiclp->isInteger(i)) {
        if(dj[i]>gap+1.0e-4) {
          nFixed0++;
          upper[i]=0.0;
        }
        if(-dj[i]>gap2+1.0e-4) {
          nFixed1++;
        lower[i]=1.0;
        }
      }
    }
    printf("%d artificials, %d fixed to 0, %d fixed to 1\n",nArt,nFixed0,nFixed1);
    //osiclp->getModelPtr()->setPerturbation(100);
    osiclp->setWarmStart(NULL); // set basis in osiclp
  }
  osiclp->initialSolve();

  // Switch off strong branching if wanted
  // model.setNumberStrong(0);
  // Do more strong branching if small
  model.setNumberStrong(0);
  model.setNumberBeforeTrust(0);

  // TEMP - set gap - better to see if no improvement in last few nodes
  model.setAllowableGap(600.0);
  // message levels
  model.messageHandler()->setLogLevel(2);
  model.solver()->messageHandler()->setLogLevel(2);
  // Do complete search
  
  model.branchAndBound();
	    }
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
static void generateCode(const char * fileName,int type,int preProcess)
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
    int translate[]={9999,0,0,-1,2,3};
    strcpy(line[numberLines++],"5  // Hand coded preprocessing");
    strcpy(line[numberLines++],"5  CglPreProcess process;");
    strcpy(line[numberLines++],"5  OsiSolverInterface * saveSolver=cbcModel->solver()->clone();");
    strcpy(line[numberLines++],"5  // Tell solver we are in Branch and Cut");
    strcpy(line[numberLines++],"5  saveSolver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo) ;");
    strcpy(line[numberLines++],"5  // Default set of cut generators");
    strcpy(line[numberLines++],"5  CglProbing generator1;");
    strcpy(line[numberLines++],"5  generator1.setUsingObjective(true);");
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
    strcpy(line[numberLines++],"5	     <<\" and \"<<solver2->getNumElements()<<std::endl;");
    strcpy(line[numberLines++],"5  }");
    strcpy(line[numberLines++],"5  // we have to keep solver2 so pass clone");
    strcpy(line[numberLines++],"5  solver2 = solver2->clone();");
    strcpy(line[numberLines++],"5  cbcModel->assignSolver(solver2);");
    strcpy(line[numberLines++],"5  cbcModel->initialSolve();");
  }
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
    "Redundant set default values","Solve","Restore values","Redundant restore values","Add to model"};
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
*/
