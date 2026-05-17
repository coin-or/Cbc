// Auto-extracted from CbcSolver::run() — preprocessing block
// Part of CbcSolver refactoring

#include "CbcConfig.h"
#include "CoinPragma.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>

#include "CoinTime.hpp"
#include "CoinWarmStartBasis.hpp"

#include "OsiClpSolverInterface.hpp"

#include "ClpSimplex.hpp"

#include "CglPreProcess.hpp"
#include "CglProbing.hpp"

#include "CbcModel.hpp"
#include "CbcSolver.hpp"
#include "CbcSolverStatistics.hpp"
#include "CbcMessage.hpp"
#include "CbcOutput.hpp"
#include "CbcBranchLotsize.hpp"
#include "CbcSOS.hpp"

#include "CoinTable.hpp"
#include "OsiRowCutDebugger.hpp"

extern void printGeneralMessage(CbcModel &model, std::string message, int type);
extern OsiClpSolverInterface *getClpSolver(OsiSolverInterface *solver);

int CbcSolver::preprocess(int preProcess, int cbcParamCode,
  OsiClpSolverInterface *&clpSolver, ClpSimplex *&lpSolver,
  CglPreProcess &process, CbcSolverStatistics &statistics,
  int &returnCode, ampl_info *info)
{
  CbcParameters &parameters = parameters_;
  std::ostringstream buffer;
  // lotStruct and LotStruct have identical layout
  typedef struct {
    double low;
    double high;
    int column;
  } lotStruct;
  lotStruct *lotsize = reinterpret_cast<lotStruct *>(lotsize_);
  typedef struct {
    lotStruct *lotsize;
    int numberLotSizing;
  } lotStruct2;
  lotStruct2 passSCtopreprocess;

  // Variables that cross the block boundary — declared here as locals.
  // The caller is responsible for reading them back if needed.
  CbcPreprocHandler *preprocHandler = nullptr;
  double preprocStart = 0.0;
  bool redoSOS = false;
  bool integersOK = true;
  int truncateColumns = COIN_INT_MAX;
  int truncateRows = -1;
  double *truncatedRhsLower = NULL;
  double *truncatedRhsUpper = NULL;
  int *newPriorities = NULL;

  saveSolver_ = babModel_->solver()->clone();
  OsiSolverInterface *saveSolver = saveSolver_;
  /* Do not try and produce equality cliques and
     do up to 10 passes */
  OsiSolverInterface *solver2;
  {
    // Tell solver we are in Branch and Cut
    saveSolver->setHintParam(OsiDoInBranchAndCut, true,
      OsiHintDo);
    // Default set of cut generators
    CglProbing generator1;
    generator1.setUsingObjective(1);
    generator1.setMaxPass(1);
    generator1.setMaxPassRoot(1);
    generator1.setMaxProbeRoot(
      std::min(3000, saveSolver->getNumCols()));
    generator1.setMaxElements(100);
    generator1.setMaxElementsRoot(200);
    generator1.setMaxLookRoot(50);
    if (saveSolver->getNumCols() > 3000)
      generator1.setMaxProbeRoot(123);
    generator1.setRowCuts(3);
    // switch off duplicate columns if we have a solution
    if (model_.bestSolution() /*||debugValues*/)
      tunePreProcess_ |= 4096;
    // take off top
    int tune2 = tunePreProcess_ % 10000;
    if ((tune2 & (1 | 512)) != 0) {
      // heavy probing
      generator1.setMaxPassRoot(2);
#ifndef CBC_EXPERIMENT_JJF
      generator1.setMaxElements(1000);
#else
      if ((tune2 & 512) != 0)
        generator1.setMaxElementsRoot(saveSolver->getNumCols());
      else
        generator1.setMaxElements(1000);
#endif
      generator1.setMaxProbeRoot(saveSolver->getNumCols());
      if ((tune2 & 512) != 0)
        generator1.setMaxLookRoot(std::min(saveSolver->getNumCols(), 1000));
      else
        generator1.setMaxLookRoot(std::min(saveSolver->getNumCols(), 400));
    }
    if ((babModel_->specialOptions() & 65536) != 0)
      process.setOptions(1);
    // Add in generators
    if ((model_.moreSpecialOptions() & 65536) == 0)
      process.addCutGenerator(&generator1);
    int translate[] = { 9999, 0, 0, -3, 2, 3, -2, 9999, 4, 5, 0, -2 };
    process.passInMessageHandler(babModel_->messageHandler());
    // process.messageHandler()->setLogLevel(babModel_->logLevel());
    if (info && info->numberSos && doSOS_ && statusUserFunction_[0]) {
      // SOS
      numberSOS_ = info->numberSos;
      sosStart_ = info->sosStart;
      sosIndices_ = info->sosIndices;
    }
    if (numberSOS_ && doSOS_) {
      // SOS
      int numberColumns = saveSolver->getNumCols();
      char *prohibited = new char[numberColumns];
      memset(prohibited, 0, numberColumns);
      // worth looking to see if any members can be made integer

      int numberRows = saveSolver->getNumRows();
      const CoinPackedMatrix *matrixByCol = saveSolver->getMatrixByCol();
      const double *element = matrixByCol->getElements();
      const int *row = matrixByCol->getIndices();
      const CoinBigIndex *columnStart = matrixByCol->getVectorStarts();
      const int *columnLength = matrixByCol->getVectorLengths();
      const double *columnLower = saveSolver->getColLower();
      const double *columnUpper = saveSolver->getColUpper();
      const double *rowLower = saveSolver->getRowLower();
      const double *rowUpper = saveSolver->getRowUpper();
      double *sameElement = new double[numberRows];
      int *rowCount = new int[2 * numberRows];
      int *rowUsed = rowCount + numberRows;
      memset(sameElement, 0, numberRows * sizeof(double));
      memset(rowCount, 0, numberRows * sizeof(int));
      int numberInteresting1 = 0;
      int numberInteresting2 = 0;
      int numberChanged = 0;
      for (int iSet = 0; iSet < numberSOS_; iSet++) {
        if (sosType_[iSet] != 1) {
          for (int i = sosStart_[iSet]; i < sosStart_[iSet + 1];
            i++) {
            numberInteresting2++;
            int iColumn = sosIndices_[i];
            prohibited[iColumn] = 1;
          }
        } else {
          int nUsed = 0;
          for (int i = sosStart_[iSet]; i < sosStart_[iSet + 1];
            i++) {
            int iColumn = sosIndices_[i];
            for (CoinBigIndex j = columnStart[iColumn];
              j < columnStart[iColumn] + columnLength[iColumn];
              j++) {
              int iRow = row[j];
              double el = element[j];
              if (rowCount[iRow]) {
                if (el != sameElement[iRow])
                  sameElement[iRow] = 0.0;
              } else {
                sameElement[iRow] = el;
                rowUsed[nUsed++] = iRow;
              }
              rowCount[iRow]++;
            }
          }
          int nInSet = sosStart_[iSet + 1] - sosStart_[iSet];
          double nonzeroValue = COIN_DBL_MAX;
          for (int iUsed = 0; iUsed < nUsed; iUsed++) {
            int iRow = rowUsed[iUsed];
            if (rowCount[iRow] == nInSet && sameElement[iRow] && rowLower[iRow] == rowUpper[iRow]) {
              // all entries must be 0.0 or xx
              nonzeroValue = rowLower[iRow] / sameElement[iRow];
            }
            rowCount[iRow] = 0;
            sameElement[iRow] = 0.0;
          }
          if (nonzeroValue != COIN_DBL_MAX) {
            // could do scaling otherwise
            if (fabs(nonzeroValue - 1.0) < 1.0e-8) {
              for (int i = sosStart_[iSet]; i < sosStart_[iSet + 1];
                i++) {
                int iColumn = sosIndices_[i];
                if (columnUpper[iColumn] < 0.0 || columnLower[iColumn] > 1.0) {
                  printf("sos says infeasible\n");
                }
                if (!saveSolver->isInteger(iColumn)) {
                  numberChanged++;
                  saveSolver->setInteger(iColumn);
                }
                if (columnUpper[iColumn] < 1.0)
                  saveSolver->setColUpper(iColumn, 0.0);
                else
                  saveSolver->setColUpper(iColumn, 1.0);
                if (columnLower[iColumn] > 0.0)
                  saveSolver->setColLower(iColumn, 1.0);
                else
                  saveSolver->setColLower(iColumn, 0.0);
#ifndef DO_LESS_PROHIBITED
                prohibited[iColumn] = 1;
#endif
              }
            } else {
              for (int i = sosStart_[iSet]; i < sosStart_[iSet + 1];
                i++) {
                int iColumn = sosIndices_[i];
#ifndef DO_LESS_PROHIBITED
                prohibited[iColumn] = 1;
#endif
              }
            }
          } else {
            for (int i = sosStart_[iSet]; i < sosStart_[iSet + 1];
              i++) {
              int iColumn = sosIndices_[i];
              if (!saveSolver->isInteger(iColumn))
                numberInteresting1++;
#ifdef DO_LESS_PROHIBITED
              if (!saveSolver->isInteger(iColumn))
#endif
                prohibited[iColumn] = 1;
            }
          }
        }
      }
      // SOS info shown in Problem loading section; skip redundant message here.
      delete[] sameElement;
      delete[] rowCount;
      process.passInProhibited(prohibited, numberColumns);
      delete[] prohibited;
    }
    if (0) {

      // Special integers
      int numberColumns = saveSolver->getNumCols();
      char *prohibited = new char[numberColumns];
      memset(prohibited, 0, numberColumns);
      const CoinPackedMatrix *matrix = saveSolver->getMatrixByCol();
      const int *columnLength = matrix->getVectorLengths();
      int numberProhibited = 0;
      for (int iColumn = numberColumns - 1; iColumn >= 0;
        iColumn--) {
        if (!saveSolver->isInteger(iColumn) || columnLength[iColumn] > 1)
          break;
        numberProhibited++;
        prohibited[iColumn] = 1;
      }
      if (numberProhibited) {
        process.passInProhibited(prohibited, numberColumns);
        printf("**** Treating last %d integers as special - give "
               "high priority?\n",
          numberProhibited);
      }
      delete[] prohibited;
    }
    if (!model_.numberObjects() && true) {
      /* model may not have created objects
         If none then create
      */
      model_.findIntegers(true);
    }
    // Lotsizing
    if (numberLotSizing_) {
      int numberColumns = saveSolver->getNumCols();
      int numberRows = saveSolver->getNumRows();
#if 0
      // Create model which uses 0-1 variables
      {
        double * els = new double[8*numberLotSizing_];
        int * cols = new int[4*numberLotSizing_];
        CoinBigIndex * starts = new CoinBigIndex [2*numberLotSizing_+1];
        // add 0-1 variables
        double * tup = els+4*numberLotSizing_;
        double * tlo = tup+2*numberLotSizing_;
        for (int i=0;i<numberLotSizing_;i++) {
          tlo[i]=0.0;
          tup[i]=1.0;
          els[i]=0.0; // objective
        }
        OsiClpSolverInterface *si = getClpSolver(saveSolver);
        assert(si != NULL);
        // get clp itself
        ClpSimplex *model = si->getModelPtr();
        model->tightenPrimalBounds(0.0,0);
        model->addColumns(numberLotSizing_,tlo,tup,els,NULL,NULL,NULL);
        for (int i=0;i<numberLotSizing_;i++)
          model->setInteger(i+numberColumns);
        starts[0]=0;
        CoinBigIndex n=0;
        const double * columnUpper = model->columnUpper();
        for (int i=0;i<numberLotSizing_;i++) {
          int iColumn = lotsize[i].column;
          tlo[n] = -COIN_DBL_MAX;
          tup[n]=0.0;
          els[2*n] = 1.0;
          // should be able to get correct value
          els[2*n+1] = -std::min(10000.0,columnUpper[iColumn]);
          cols[2*n] = iColumn;
          cols[2*n+1] = i+numberColumns;
          n++;
          starts[n]=2*n;
          tlo[n] = 0.0;
          tup[n]=COIN_DBL_MAX;
          els[2*n] = 1.0;
          els[2*n+1] = -lotsize[i].low;
          cols[2*n] = iColumn;
          cols[2*n+1] = i+numberColumns;
          n++;
          starts[n]=2*n;
        }
        model->addRows(n,tlo,tup,starts,cols,els);
        model->writeMps("equivalentInteger.mps");
        exit(1);
      }
#endif
      char *prohibited = new char[numberColumns + numberRows];
      char *prohibitedRow = prohibited + numberColumns;
      memset(prohibited, 0, numberColumns + numberRows);
      const CoinPackedMatrix *matrix = saveSolver->getMatrixByCol();
      // const double *element = matrix->getElements();
      const int *row = matrix->getIndices();
      const CoinBigIndex *columnStart = matrix->getVectorStarts();
      const int *columnLength = matrix->getVectorLengths();
      for (int i = 0; i < numberLotSizing_; i++) {
        int iColumn = lotsize[i].column;
        prohibited[iColumn] = 2;
#if 0
        // also leave all rows
        for (CoinBigIndex j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          prohibitedRow[iRow] = -2;
        }
#endif
      }
      process.passInProhibited(prohibited, numberColumns);
      process.passInRowTypes(prohibitedRow, numberRows);
      passSCtopreprocess.lotsize = lotsize;
      passSCtopreprocess.numberLotSizing = numberLotSizing_;
      process.setApplicationData(&passSCtopreprocess);
      process.setOptions(128 | 256 | process.options());
      delete[] prohibited;
    }
    if (model_.numberObjects()) {
      OsiObject **oldObjects = babModel_->objects();
      int numberOldObjects = babModel_->numberObjects();
      if (!numberOldObjects) {
        oldObjects = model_.objects();
        numberOldObjects = model_.numberObjects();
      }
      // SOS
      int numberColumns = saveSolver->getNumCols();
      char *prohibited = new char[numberColumns];
      memset(prohibited, 0, numberColumns);
      int numberProhibited = 0;
      for (int iObj = 0; iObj < numberOldObjects; iObj++) {
        CbcSOS *obj = dynamic_cast< CbcSOS * >(oldObjects[iObj]);
        if (obj) {
          int n = obj->numberMembers();
          const int *which = obj->members();
          for (int i = 0; i < n; i++) {
            int iColumn = which[i];
#ifdef DO_LESS_PROHIBITED
            if (!saveSolver->isInteger(iColumn))
#endif
              prohibited[iColumn] = 1;
            numberProhibited++;
          }
        }
        CbcLotsize *obj2 = dynamic_cast< CbcLotsize * >(oldObjects[iObj]);
        if (obj2) {
          int iColumn = obj2->columnNumber();
          prohibited[iColumn] = 1;
          numberProhibited++;
        }
      }
      if (numberProhibited)
        process.passInProhibited(prohibited, numberColumns);
      delete[] prohibited;
    }
    int numberPasses = 10;
    if (doSprint_ > 0) {
      // Sprint for primal solves
      ClpSolve::SolveType method = ClpSolve::usePrimalorSprint;
      ClpSolve::PresolveType presolveType = ClpSolve::presolveOff;
      int numberPasses = 5;
      int options[] = { 0, 3, 0, 0, 0, 0 };
      int extraInfo[] = { -1, 20, -1, -1, -1, -1 };
      extraInfo[1] = doSprint_;
      int independentOptions[] = { 0, 0, 3 };
      ClpSolve clpSolve(method, presolveType, numberPasses,
        options, extraInfo, independentOptions);
      // say use in OsiClp
      clpSolve.setSpecialOption(6, 1);
      OsiClpSolverInterface *osiclp = getClpSolver(saveSolver);
      osiclp->setSolveOptions(clpSolve);
      osiclp->setHintParam(OsiDoDualInResolve, false);
      // switch off row copy
      osiclp->getModelPtr()->setSpecialOptions(
        osiclp->getModelPtr()->specialOptions() | 256);
      osiclp->getModelPtr()->setInfeasibilityCost(1.0e11);
    }
    {
      OsiClpSolverInterface *osiclp = getClpSolver(saveSolver);
      osiclp->setSpecialOptions(osiclp->specialOptions() | 1024);
      osiclp->getModelPtr()->setClpScaledMatrix(NULL); // safer
      int savePerturbation = osiclp->getModelPtr()->perturbation();
      // #define CBC_TEMP1
#ifdef CBC_TEMP1
      if (savePerturbation == 50)
        osiclp->getModelPtr()->setPerturbation(52); // try less
#endif
      if ((model_.moreSpecialOptions() & 65536) != 0)
        process.setOptions(2 + 4 + 8); // no cuts
      // cbcPreProcessPointer removed (dead code)
      redoSOS = true;
      int saveOptions = osiclp->getModelPtr()->moreSpecialOptions();
      if ((model_.specialOptions() & 16777216) != 0 && model_.getCutoff() > 1.0e30) {
        osiclp->getModelPtr()->setMoreSpecialOptions(saveOptions | 262144);
      }
#if DEBUG_PREPROCESS > 1
      if (debugValues_) {
        process.setApplicationData(
          const_cast< double * >(debugValues_));
      }
#endif
      if (parameters[CbcParam::DEBUGFILE]->fileName() == "unitTest") {
        // This is probably wrong, will need to debug
        babModel_->solver()->activateRowCutDebugger(saveInputQueue_[0].c_str());
        OsiRowCutDebugger *debugger = babModel_->solver()->getRowCutDebuggerAlways();
        numberDebugValues_ = babModel_->getNumCols();
        debugValues_ = CoinCopyOfArray(debugger->optimalSolution(),
          numberDebugValues_);
      }
      redoSOS = true;
      bool keepPPN = parameters[CbcParam::PREPROCNAMES]->modeVal();
#ifdef SAVE_NAUTY
      keepPPN = 1;
#endif
      process.setKeepColumnNames(keepPPN);
      process.setTimeLimit(babModel_->getMaximumSeconds() - babModel_->getCurrentSeconds(),
        babModel_->useElapsedTime());
      if (model_.getKeepNamesPreproc())
        process.setKeepColumnNames(true);
      if (keepPPN)
        babModel_->setKeepNamesPreproc(1);
      setPreProcessingMode(saveSolver, 1);

      // Install preprocessing output handler
      {
        int ll = babModel_->messageHandler()->logLevel();
        if (ll >= 1) {
          FILE *fp = babModel_->messageHandler()->filePointer();
          bool u8 = CbcOutput::useUtf8();
          preprocHandler = new CbcPreprocHandler(fp, u8, ll);
          process.passInMessageHandler(preprocHandler);
          fprintf(fp, "\n%s\n\n",
            CoinTable::phaseStart("Preprocessing", u8).c_str());
          preprocStart = CoinWallclockTime();
          fflush(fp);
        }
      }
#if CBC_USE_PAPILO
      extern void zapPapilo(int pOptions, CglPreProcess *process);
      int pOptions = 0;
      int tune2 = preProcess;
      // Convert to minimize if papilo
      bool maximize = false;
      if (tune2 > 11) {
        OsiClpSolverInterface *clpSolver = getClpSolver(saveSolver);
        if (clpSolver->getObjSense() == -1.0) {
          maximize = true;
          clpSolver->setObjSense(1.0);
          double objOffset;
          clpSolver->getDblParam(OsiObjOffset, objOffset);
          int numberColumns = clpSolver->getNumCols();
          double *objective = clpSolver->getModelPtr()->objective();
          for (int i = 0; i < numberColumns; i++)
            objective[i] = -objective[i];
          clpSolver->setDblParam(OsiObjOffset, -objOffset);
        }
        bool stopAfter = false;
        if (tune2 > 15) {
          preProcess = 10; // say stop
          tune2 -= 4;
        } else {
          preProcess = 1;
        }
#ifdef CBC_THREAD
        pOptions = (tune2 & 1) != 0 ? 2 : 1;
#endif
        if ((tune2 & 2) == 0)
          pOptions |= 8; // at beginning
        else
          pOptions |= 16; // at end
      }
      zapPapilo(pOptions, &process);
#endif
      solver2 = process.preProcessNonDefault(*saveSolver, translate[preProcess], numberPasses,
        tunePreProcess_);
      if (!solver2) {
        // Case A: preprocessing itself detected infeasibility —
        // retry with simpler settings (no double check possible)
        process.clean();
        solver2 = process.preProcessNonDefault(*saveSolver,
          0, 99, 0);
      } else if (!solver2->isProvenOptimal()) {
        /* Infeasible - but most real problems are not
           infeasible - so try simpler preprocessing which
           is less affected by tolerance issues */
        // Case B: LP relaxation infeasible — double check with resolve()
        solver2->resolve();
        if (!solver2->isProvenOptimal()) {
          process.clean();
          solver2 = process.preProcessNonDefault(*saveSolver,
            0, 99, 0);
        }
      }
      setPreProcessingMode(saveSolver, 0);

      // Print preprocessing table summary
      if (preprocHandler)
        preprocHandler->printTableEnd();
#if CBC_USE_PAPILO
      // Convert back
      if (maximize) {
        OsiClpSolverInterface *clpSolver = getClpSolver(saveSolver);
        double objOffset;
        clpSolver->setObjSense(-1.0);
        clpSolver->getDblParam(OsiObjOffset, objOffset);
        int numberColumns = clpSolver->getNumCols();
        double *objective = clpSolver->getModelPtr()->objective();
        for (int i = 0; i < numberColumns; i++)
          objective[i] = -objective[i];
        clpSolver->setDblParam(OsiObjOffset, -objOffset);
        if (solver2) {
          OsiClpSolverInterface *clpSolver = getClpSolver(solver2);
          double objOffset;
          clpSolver->setObjSense(-1.0);
          clpSolver->getDblParam(OsiObjOffset, objOffset);
          int numberColumns = clpSolver->getNumCols();
          double *objective = clpSolver->getModelPtr()->objective();
          for (int i = 0; i < numberColumns; i++)
            objective[i] = -objective[i];
          clpSolver->setDblParam(OsiObjOffset, -objOffset);
        }
      }
#endif
      if (solver2) {
        setPreProcessingMode(solver2, 0);
        model_.setOriginalColumns(process.originalColumns(),
          solver2->getNumCols());

        osiclp->getModelPtr()->setPerturbation(savePerturbation);
        osiclp->getModelPtr()->setMoreSpecialOptions(saveOptions);
        /* clean solvers - should be done in preProcess but
           that doesn't know about Clp */
        OsiClpSolverInterface *solver;
        solver = dynamic_cast< OsiClpSolverInterface * >(solver2);
        solver->getModelPtr()->cleanScalingEtc();
        solver = getClpSolver(process.originalModel());
        solver->getModelPtr()->cleanScalingEtc();
        solver = getClpSolver(process.startModel());
        if (solver)
          solver->getModelPtr()->cleanScalingEtc();
        int numberSolvers = process.numberSolvers();
        if (numberSolvers == 99)
          numberSolvers = 1; // really just 1
        // some of these may be same
        for (int i = 0; i < numberSolvers; i++) {
          solver = getClpSolver(process.modelAtPass(i));
          if (solver)
            solver->getModelPtr()->cleanScalingEtc();
          solver = getClpSolver(process.modifiedModel(i));
          if (solver)
            solver->getModelPtr()->cleanScalingEtc();
        }
      }
    }
    integersOK = false; // We need to redo if CbcObjects exist
    // Tell solver we are not in Branch and Cut
    saveSolver->setHintParam(OsiDoInBranchAndCut, false,
      OsiHintDo);
    if (solver2)
      solver2->setHintParam(OsiDoInBranchAndCut, false,
        OsiHintDo);
  }
  if (info && !solver2 && statusUserFunction_[0]) {
    // infeasible
    info->problemStatus = 1;
    info->objValue = 1.0e100;
    sprintf(info->buffer,
      "infeasible/unbounded by pre-processing");
    info->primalSolution = NULL;
    info->dualSolution = NULL;
    if (preprocHandler) {
      preprocHandler->markInfeasible("infeasible or unbounded");
      process.passInMessageHandler(model_.messageHandler());
      delete preprocHandler;
      preprocHandler = nullptr;
    }
  }
  if (!solver2) {
    printGeneralMessage(model_,
      "Pre-processing says infeasible or unbounded");
    if (preprocHandler)
      preprocHandler->markInfeasible("infeasible or unbounded");
    // say infeasible for solution
    integerStatus_ = 6;
    delete saveSolver_;
    saveSolver_ = NULL;
    model_.setProblemStatus(0);
    model_.setSecondaryStatus(1);
    babModel_->setProblemStatus(0);
    babModel_->setSecondaryStatus(1);
  } else {
    statistics.nprocessedrows = solver2->getNumRows();
    statistics.nprocessedcols = solver2->getNumCols();
    model_.setProblemStatus(-1);
    babModel_->setProblemStatus(-1);
  }
  returnCode = 0;
  if (callBack_ != NULL)
    returnCode = callBack_->callBack(babModel_, 2);
  if (returnCode) {
    // exit if user wants — close preprocessing section first
    if (preprocHandler) {
      process.passInMessageHandler(model_.messageHandler());
      delete preprocHandler;
      preprocHandler = nullptr;
    }
    delete babModel_;
    babModel_ = NULL;
    return 3;
  }
  if (!solver2) {
    // preprocessing says infeasible — section closed by destructor below
    process.passInMessageHandler(model_.messageHandler());
    delete preprocHandler;
    preprocHandler = nullptr;
    return 1;
  }
  if (model_.bestSolution()) {
    // need to redo - in case no better found in BAB
    // just get integer part right
    const int *originalColumns = process.originalColumns();
    int numberColumns = std::min(solver2->getNumCols(), babModel_->getNumCols());
#if 0
    double *bestSolution = babModel_->bestSolution();
    const double *oldBestSolution = model_.bestSolution();
    for (int i = 0; i < numberColumns; i++) {
      int jColumn = originalColumns[i];
      bestSolution[i] = oldBestSolution[jColumn];
    }
#else
    int numberColumnsB = babModel_->getNumCols();
    int numberColumns2 = std::max(solver2->getNumCols(), numberColumnsB);
    double *bestSolution = new double[numberColumns2];
    memset(bestSolution, 0, numberColumns2 * sizeof(double));
    const double *oldBestSolution = model_.bestSolution();
    for (int i = 0; i < numberColumns; i++) {
      int jColumn = originalColumns[i];
      if (jColumn < numberColumnsB)
        bestSolution[i] = oldBestSolution[jColumn];
    }
    double obj = model_.getObjValue();
    double newCutoff = std::min(model_.getCutoff(), obj + 1.0e-4);
    babModel_->setBestSolution(bestSolution, numberColumns, 1.0e10, false);
    babModel_->setCutoff(newCutoff);
    delete[] bestSolution;
#endif
  }
  // solver2->resolve();
#ifdef CBC_NAMES_FOR_COMPARE
  {
    OsiClpSolverInterface *solver = getClpSolver(solver2);
    ClpSimplex *simplex = solver->getModelPtr();
    int numberRows = simplex->numberRows();
    int numberColumns = simplex->numberColumns();
    int numberIntegers = 0;
    for (int i = 0; i < numberColumns; i++) {
      if (simplex->isInteger(i))
        numberIntegers++;
    }
    printf("CbC %s after preprocessing %d rows, %d columns, %d integers\n",
      simplex->problemName().c_str(), numberRows,
      numberColumns, numberIntegers);
  }
#endif
  if (preProcess == 2 || preProcess >= 10) {
    // names are wrong - redo
    const int *originalColumns = process.originalColumns();
    int numberColumns = solver2->getNumCols();
    OsiSolverInterface *originalSolver = model_.solver();
    int numberOriginalColumns = originalSolver->getNumCols();
    int oddColumn = 1;
    for (int i = 0; i < numberColumns; i++) {
      int iColumn = originalColumns[i];
      if (iColumn < numberOriginalColumns) {
        std::string columnName = originalSolver->getColName(iColumn);
        // check if odd name
        if (columnName[0] == 'N') {
          try {
            int sequence = std::atoi(columnName.substr(1).c_str());
            oddColumn = sequence + 1;
          } catch (...) {
          }
        }
        solver2->setColName(i, columnName);
      } else {
        char name[15];
        sprintf(name, "N%.7d", oddColumn);
        oddColumn++;
        solver2->setColName(i, name);
      }
    }
    OsiClpSolverInterface *clpSolver2 = getClpSolver(solver2);
    ClpSimplex *lpSolver2 = clpSolver2->getModelPtr();
    char name[100];
    if (preProcess == 2) {
      strcpy(name, "presolved.mps");
    } else {
      strcpy(name,
        parameters[CbcParam::IMPORTFILE]->fileName().c_str());
      char *dot = strstr(name, ".mps");
      if (!dot)
        dot = strstr(name, ".lp");
      if (preProcess >= 10 && !dot) {
        // get some sort of name
        dot = name + strlen(name);
      }
      if (dot) {
        *dot = '\0';
        int n = static_cast< int >(dot - name);
        int i;
        for (i = n - 1; i >= 0; i--) {
          if (name[i] == '/')
            break;
        }
        if (i >= 0)
          memmove(name, name + i + 1, n);
        strcat(name, "_preprocessed.mps");
      } else {
        strcpy(name, "preprocessed.mps");
      }
    }
    lpSolver2->writeMps(name, 0, 1, lpSolver2->optimizationDirection());
    printf("Preprocessed model (minimization) on %s - size %d %d \n",
      name, lpSolver2->getNumRows(), lpSolver2->getNumCols());
    if (preProcess >= 10) {
      printf("user wanted to stop\n");
      delete saveSolver_;
      saveSolver_ = NULL;
      delete babModel_;
      babModel_ = NULL;
      returnCode = 0;
      return 3;
    }
  }
  {
    // look at new integers
    int numberOriginalColumns = process.originalModel()->getNumCols();
    const int *originalColumns = process.originalColumns();
    OsiClpSolverInterface *osiclp2 = getClpSolver(solver2);
    int numberColumns = osiclp2->getNumCols();
    OsiClpSolverInterface *osiclp = getClpSolver(saveSolver_);
    for (int i = 0; i < numberColumns; i++) {
      int iColumn = originalColumns[i];
      if (iColumn < numberOriginalColumns) {
        if (osiclp2->isInteger(i) && !osiclp->isInteger(iColumn))
          osiclp2->setOptionalInteger(i); // say optional
      }
    }
    // do lotsizing
    if (numberLotSizing_) {
      const double *columnLower = solver2->getColLower();
      const double *columnUpper = solver2->getColUpper();
      CbcObject **objects = new CbcObject *[numberLotSizing_];
      double points[] = { 0.0, 0.0, 0.0, 0.0 };
      int *back = new int[numberOriginalColumns];
      for (int i = 0; i < numberOriginalColumns; i++)
        back[i] = -1;
      for (int i = 0; i < numberColumns; i++) {
        int iColumn = originalColumns[i];
        back[iColumn] = i;
      }
      int n = numberLotSizing_;
      numberLotSizing_ = 0;
      for (int i = 0; i < n; i++) {
        int iColumn = lotsize[i].column;
        iColumn = back[iColumn];
        if (iColumn >= 0) {
          if (columnLower[iColumn] < lotsize[i].low && columnUpper[iColumn] > lotsize[i].low - 1.0e-8) {
            points[2] = lotsize[i].low;
            points[3] = columnUpper[iColumn];
            objects[numberLotSizing_++] = new CbcLotsize(babModel_, iColumn, 2,
              points, true);
          } else {
            // printf("SC %d new bounds %g,%g was %d low %g high %g\n",iColumn,
            //	 columnLower[iColumn],columnUpper[iColumn],
            //	 lotsize[i].column,lotsize[i].low,lotsize[i].high);
            if (columnUpper[iColumn] < lotsize[i].low) {
              if (columnLower[iColumn]) {
                printf("Infeasible due to SC variables?!\n");
              } else {
                solver2->setColUpper(iColumn, 0.0);
              }
            }
          }
        }
      }
      delete[] back;
      babModel_->addObjects(numberLotSizing_, objects);
      for (int i = 0; i < numberLotSizing_; i++)
        delete objects[i];
      delete[] objects;
      if (numberLotSizing_ < n) {
        buffer.str("");
        buffer << "Pre-processing reduced number of SC"
               << " variables from " << n << " to "
               << numberLotSizing_;
        printGeneralMessage(model_, buffer.str());
      }
    }
    // redo existing SOS
    if (osiclp->numberSOS()) {
      redoSOS = false;
      int *back = new int[numberOriginalColumns];
      for (int i = 0; i < numberOriginalColumns; i++)
        back[i] = -1;
      for (int i = 0; i < numberColumns; i++) {
        int iColumn = originalColumns[i];
        back[iColumn] = i;
      }
      int numberSOSOld = osiclp->numberSOS();
      int numberSOS = osiclp2->numberSOS();
      assert(numberSOS == numberSOSOld);
      CoinSet *setInfo = const_cast< CoinSet * >(osiclp2->setInfo());
      for (int i = 0; i < numberSOS; i++) {
        // int type = setInfo[i].setType();
        int n = setInfo[i].numberEntries();
        int *which = setInfo[i].modifiableWhich();
#ifndef DO_LESS_PROHIBITED
        for (int j = 0; j < n; j++) {
          int iColumn = which[j];
          iColumn = back[iColumn];
          assert(iColumn >= 0);
          which[j] = iColumn;
        }
#else
        double *weights = setInfo[i].modifiableWeights();
        int n2 = 0;
        for (int j = 0; j < n; j++) {
          int iColumn = which[j];
          iColumn = back[iColumn];
          if (iColumn >= 0) {
            which[n2] = iColumn;
            weights[n2++] = weights[j];
          }
        }
        setInfo[i].setNumberEntries(n2);
#endif
      }
      delete[] back;
    }
  }
  // we have to keep solver2 so pass clone
  solver2 = solver2->clone();
  // see if extra variables wanted
  int threshold = parameters[CbcParam::EXTRAVARIABLES]->intVal();
  int more2 = parameters[CbcParam::MOREMOREMIPOPTIONS]->intVal();
  if (threshold || (more2 & (512 | 1024)) != 0) {
    int numberColumns = solver2->getNumCols();
    truncateRows = solver2->getNumRows();
    bool modifiedModel = false;
    int highPriority = 0;
    /*
      normal - no priorities
      >10000 equal high priority
      >20000 higher priority for higher cost
    */
    if (threshold > 10000) {
      highPriority = threshold / 10000;
      threshold -= 10000 * highPriority;
    }
    // If 1000 set then don't put obj on ne variables
    bool moveObjective = true;
    if (threshold > 1000) {
      moveObjective = false;
      threshold -= 1000;
    }
    const double *columnLower = solver2->getColLower();
    const double *columnUpper = solver2->getColUpper();
    const double *objective = solver2->getObjCoefficients();
    int numberIntegers = 0;
    int numberBinary = 0;
    int numberTotalIntegers = 0;
    double *obj = new double[numberColumns];
    int *which = new int[numberColumns];
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (solver2->isInteger(iColumn)) {
        numberTotalIntegers++;
        if (columnUpper[iColumn] > columnLower[iColumn]) {
          numberIntegers++;
          if (columnLower[iColumn] == 0.0 && columnUpper[iColumn] == 1)
            numberBinary++;
        }
      }
    }
    int numberSort = 0;
    int numberZero = 0;
    int numberZeroContinuous = 0;
    int numberDifferentObj = 0;
    int numberContinuous = 0;
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (columnUpper[iColumn] > columnLower[iColumn]) {
        if (solver2->isInteger(iColumn)) {
          if (!objective[iColumn]) {
            numberZero++;
          } else {
            obj[numberSort] = fabs(objective[iColumn]);
            which[numberSort++] = iColumn;
          }
        } else if (objective[iColumn]) {
          numberContinuous++;
        } else {
          numberZeroContinuous++;
        }
      }
    }
    CoinSort_2(obj, obj + numberSort, which);
    double last = obj[0];
    for (int jColumn = 1; jColumn < numberSort; jColumn++) {
      if (fabs(obj[jColumn] - last) > 1.0e-12) {
        numberDifferentObj++;
        last = obj[jColumn];
      }
    }
    numberDifferentObj++;
    if (numberDifferentObj <= threshold + (numberZero)
        ? 1
        : 0 && numberDifferentObj) {
      int *backward = NULL;
      if (highPriority) {
        newPriorities = new int[numberTotalIntegers + numberDifferentObj + numberColumns];
        backward = newPriorities + numberTotalIntegers + numberDifferentObj;
        numberTotalIntegers = 0;
        for (int iColumn = 0; iColumn < numberColumns;
          iColumn++) {
          if (solver2->isInteger(iColumn)) {
            backward[iColumn] = numberTotalIntegers;
            newPriorities[numberTotalIntegers++] = 10000;
          }
        }
      }
      int iLast = 0;
      double last = obj[0];
      for (int jColumn = 1; jColumn < numberSort; jColumn++) {
        if (fabs(obj[jColumn] - last) > 1.0e-12) {
          buffer.str("");
          buffer << jColumn - iLast << " variables have objective of "
                 << last;
          printGeneralMessage(model_, buffer.str());
          iLast = jColumn;
          last = obj[jColumn];
        }
      }
      buffer.str("");
      buffer << numberSort - iLast << " variables have objective of "
             << last;
      printGeneralMessage(model_, buffer.str());
      int spaceNeeded = numberSort + numberDifferentObj;
      CoinBigIndex *columnAddDummy = new CoinBigIndex[numberDifferentObj + 1];
      int *columnAdd = new int[spaceNeeded];
      double *elementAdd = new double[spaceNeeded];
      CoinBigIndex *rowAdd = new CoinBigIndex[numberDifferentObj + 1];
      double *objectiveNew = new double[3 * numberDifferentObj];
      double *lowerNew = objectiveNew + numberDifferentObj;
      double *upperNew = lowerNew + numberDifferentObj;
      memset(columnAddDummy, 0,
        (numberDifferentObj + 1) * sizeof(CoinBigIndex));
      iLast = 0;
      last = obj[0];
      numberDifferentObj = 0;
      int priorityLevel = 9999;
      int numberElements = 0;
      rowAdd[0] = 0;
      for (int jColumn = 1; jColumn < numberSort + 1; jColumn++) {
        if (jColumn == numberSort || fabs(obj[jColumn] - last) > 1.0e-12) {
          // not if just one
          if (jColumn - iLast > 1) {
            // do priority
            if (highPriority == 1) {
              newPriorities[numberTotalIntegers + numberDifferentObj] = 500;
            } else if (highPriority == 2) {
              newPriorities[numberTotalIntegers + numberDifferentObj] = priorityLevel;
              priorityLevel--;
            }
            int iColumn = which[iLast];
            objectiveNew[numberDifferentObj] = objective[iColumn];
            double lower = 0.0;
            double upper = 0.0;
            for (int kColumn = iLast; kColumn < jColumn;
              kColumn++) {
              iColumn = which[kColumn];
              if (moveObjective)
                solver2->setObjCoeff(iColumn, 0.0);
              double lowerValue = columnLower[iColumn];
              double upperValue = columnUpper[iColumn];
              double elementValue = -1.0;
              if (objectiveNew[numberDifferentObj] * objective[iColumn] < 0.0) {
                lowerValue = -columnUpper[iColumn];
                upperValue = -columnLower[iColumn];
                elementValue = 1.0;
              }
              if (!moveObjective)
                objectiveNew[numberDifferentObj] = 0.0;
              columnAdd[numberElements] = iColumn;
              elementAdd[numberElements++] = elementValue;
              if (lower != -COIN_DBL_MAX) {
                if (lowerValue != -COIN_DBL_MAX)
                  lower += lowerValue;
                else
                  lower = -COIN_DBL_MAX;
              }
              if (upper != COIN_DBL_MAX) {
                if (upperValue != COIN_DBL_MAX)
                  upper += upperValue;
                else
                  upper = COIN_DBL_MAX;
              }
            }
            columnAdd[numberElements] = numberColumns + numberDifferentObj;
            elementAdd[numberElements++] = 1.0;
            lowerNew[numberDifferentObj] = lower;
            upperNew[numberDifferentObj] = upper;
            numberDifferentObj++;
            rowAdd[numberDifferentObj] = numberElements;
          } else if (highPriority) {
            // just one
            // do priority
            int iColumn = which[iLast];
            int iInt = backward[iColumn];
            if (highPriority == 1) {
              newPriorities[iInt] = 500;
            } else {
              newPriorities[iInt] = priorityLevel;
              priorityLevel--;
            }
          }
          if (jColumn < numberSort) {
            iLast = jColumn;
            last = obj[jColumn];
          }
        }
      }
      if (numberDifferentObj) {
        // add columns
        solver2->addCols(numberDifferentObj, columnAddDummy, NULL,
          NULL, lowerNew, upperNew, objectiveNew);
        // add constraints and make integer if all integer in
        // group
        OsiClpSolverInterface *clpSolver2 = getClpSolver(solver2);
        for (int iObj = 0; iObj < numberDifferentObj; iObj++) {
          lowerNew[iObj] = 0.0;
          upperNew[iObj] = 0.0;
          solver2->setInteger(numberColumns + iObj);
          if (clpSolver2)
            clpSolver2->setOptionalInteger(numberColumns + iObj);
        }
        solver2->addRows(numberDifferentObj, rowAdd, columnAdd,
          elementAdd, lowerNew, upperNew);
        buffer.str("");
        buffer << "Replacing model - "
               << numberDifferentObj << " new variables";
        modifiedModel = true;
      }
      delete[] columnAdd;
      delete[] columnAddDummy;
      delete[] elementAdd;
      delete[] rowAdd;
      delete[] objectiveNew;
    }
    delete[] which;
    delete[] obj;
    if ((more2 & (512 | 1024)) != 0) {
      // try for row slacks etc
      // later do row branching
      int iRow, iColumn;
      int numberColumns = solver2->getNumCols();
      int numberRows = solver2->getNumRows();
      int fudgeObjective = more2 & 512;
      int addSlacks = more2 & 1024;
      if (fudgeObjective) {
        bool moveObj = false;
        fudgeObjective = 0;
        const double *objective = solver2->getObjCoefficients();
        const double *columnLower = solver2->getColLower();
        const double *columnUpper = solver2->getColUpper();
        double *newValues = new double[numberColumns + 1];
        int *newColumn = new int[numberColumns + 1];
        bool allInteger = true;
        int n = 0;
        double newLower = 0.0;
        double newUpper = 0.0;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (objective[iColumn]) {
            if (!solver2->isInteger(iColumn)) {
              allInteger = false;
              break;
            } else {
              double value = objective[iColumn];
              double nearest = floor(value + 0.5);
              if (fabs(value - nearest) > 1.0e-8) {
                allInteger = false;
                break;
              } else {
                newValues[n] = nearest;
                newColumn[n++] = iColumn;
                if (nearest > 0.0) {
                  newLower += std::max(columnLower[iColumn], -1.0e20) * nearest;
                  newUpper += std::min(columnUpper[iColumn], 1.0e20) * nearest;
                } else {
                  newUpper += std::max(columnLower[iColumn], -1.0e20) * nearest;
                  newLower += std::min(columnUpper[iColumn], 1.0e20) * nearest;
                }
              }
            }
          }
        }
        if (allInteger && n) {
          fudgeObjective = n;
          solver2->addCol(0, NULL, NULL, newLower, newUpper, 0.0,
            "obj_col");
          solver2->setInteger(numberColumns);
          newValues[n] = -1.0;
          newColumn[n++] = numberColumns;
          solver2->addRow(n, newColumn, newValues, 0.0, 0.0);
          if (moveObj) {
            memset(newValues, 0, numberColumns * sizeof(double));
            newValues[numberColumns] = 1.0;
            solver2->setObjective(newValues);
          }
          numberRows++;
          numberColumns++;
        }
        delete[] newValues;
        delete[] newColumn;
      }
      if (addSlacks) {
        bool moveObj = false;
        addSlacks = 0;
        // get row copy
        const CoinPackedMatrix *matrix = solver2->getMatrixByRow();
        const double *element = matrix->getElements();
        const int *column = matrix->getIndices();
        const CoinBigIndex *rowStart = matrix->getVectorStarts();
        const int *rowLength = matrix->getVectorLengths();
        const double *rowLower = solver2->getRowLower();
        const double *rowUpper = solver2->getRowUpper();
        const double *columnLower = solver2->getColLower();
        const double *columnUpper = solver2->getColUpper();

        // maximum space for additional columns
        CoinBigIndex *newColumnStart = new CoinBigIndex[numberRows + 1];
        newColumnStart[0] = 0;
        int *newRow = new int[numberRows];
        double *newElement = new double[numberRows];
        double *newObjective = new double[numberRows];
        double *newColumnLower = new double[numberRows];
        double *newColumnUpper = new double[numberRows];
        double *oldObjective = CoinCopyOfArray(
          solver2->getObjCoefficients(), numberColumns);
        for (iRow = 0; iRow < numberRows; iRow++) {
          if (rowLower[iRow] != rowUpper[iRow]) {
            bool allInteger = true;
            double newLower = 0.0;
            double newUpper = 0.0;
            double constantObjective = 0.0;
            for (CoinBigIndex j = rowStart[iRow];
              j < rowStart[iRow] + rowLength[iRow]; j++) {
              int iColumn = column[j];
              if (!solver2->isInteger(iColumn)) {
                allInteger = false;
                break;
              } else {
                double value = element[j];
                double nearest = floor(value + 0.5);
                if (fabs(value - nearest) > 1.0e-8) {
                  allInteger = false;
                  break;
                } else {
                  if (!oldObjective[iColumn])
                    constantObjective = COIN_DBL_MAX;
                  if (!constantObjective) {
                    constantObjective = oldObjective[iColumn] / nearest;
                  } else if (constantObjective != COIN_DBL_MAX) {
                    double newConstant = oldObjective[iColumn] / nearest;
                    if (constantObjective > 0.0) {
                      if (newConstant <= 0.0)
                        constantObjective = COIN_DBL_MAX;
                      else
                        constantObjective = std::min(
                          constantObjective, newConstant);
                    } else {
                      if (newConstant >= 0.0)
                        constantObjective = COIN_DBL_MAX;
                      else
                        constantObjective = std::max(
                          constantObjective, newConstant);
                    }
                  }
                  if (nearest > 0.0) {
                    newLower += std::max(columnLower[iColumn], -1.0e20) * nearest;
                    newUpper += std::min(columnUpper[iColumn], 1.0e20) * nearest;
                  } else {
                    newUpper += std::max(columnLower[iColumn], -1.0e20) * nearest;
                    newLower += std::min(columnUpper[iColumn], 1.0e20) * nearest;
                  }
                }
              }
            }
            if (allInteger) {
              newColumnStart[addSlacks + 1] = addSlacks + 1;
              newRow[addSlacks] = iRow;
              newElement[addSlacks] = -1.0;
              newObjective[addSlacks] = 0.0;
              if (moveObj && constantObjective != COIN_DBL_MAX) {
                // move some of objective here if looks constant
                newObjective[addSlacks] = constantObjective;
                for (CoinBigIndex j = rowStart[iRow];
                  j < rowStart[iRow] + rowLength[iRow]; j++) {
                  int iColumn = column[j];
                  double value = element[j];
                  double nearest = floor(value + 0.5);
                  oldObjective[iColumn] -= nearest * constantObjective;
                }
              }
              newColumnLower[addSlacks] = std::max(newLower, ceil(rowLower[iRow]));
              ;
              newColumnUpper[addSlacks] = std::min(newUpper, floor(rowUpper[iRow]));
              addSlacks++;
            }
          }
        }
        if (addSlacks) {
          solver2->setObjective(oldObjective);
          solver2->addCols(addSlacks, newColumnStart, newRow,
            newElement, newColumnLower,
            newColumnUpper, newObjective);
          truncatedRhsLower = CoinCopyOfArray(solver2->getRowLower(), numberRows);
          truncatedRhsUpper = CoinCopyOfArray(solver2->getRowUpper(), numberRows);
          for (int j = 0; j < addSlacks; j++) {
            int iRow = newRow[j];
            solver2->setRowLower(iRow, 0.0);
            solver2->setRowUpper(iRow, 0.0);
            int iColumn = j + numberColumns;
            solver2->setInteger(iColumn);
            std::string name = solver2->getRowName(iRow);
            name += "_int";
            solver2->setColName(iColumn, name);
          }
        }
      }
      if (fudgeObjective || addSlacks) {
        modifiedModel = true;
        if (fudgeObjective && addSlacks) {
          buffer.str("");
          buffer << "Objective integer added with "
                 << fudgeObjective << " elements and "
                 << addSlacks << " Integer slacks added";
        } else if (fudgeObjective) {
          // just objective
          buffer.str("");
          buffer << "Objective integer added with "
                 << fudgeObjective << " elements",
            more2 &= ~1024;
        } else {
          // just slacks
          buffer.str("");
          buffer << addSlacks << " Integer slacks added",
            more2 &= ~512;
        }
      } else {
        more2 &= ~(512 | 1024);
      }
      parameters[CbcParam::MOREMOREMIPOPTIONS]->setVal(more2);
    }
    if (modifiedModel) {
      printGeneralMessage(model_, buffer.str());
      truncateColumns = numberColumns;
    }
  }
  babModel_->assignSolver(solver2);
  babModel_->setOriginalColumns(process.originalColumns(),
    truncateColumns);
  babModel_->initialSolve();
#if CBC_USE_INITIAL_TIME == 2
  // time starts from here?
  // NOTE: time1Elapsed and time1 were locals in run() — this block
  // is normally dead code (CBC_USE_INITIAL_TIME defaults to 1).
  if (babModel_->useElapsedTime())
    babModel_->setDblParam(CbcModel::CbcStartSeconds,
      CoinGetTimeOfDay());
  else
    babModel_->setDblParam(CbcModel::CbcStartSeconds,
      CoinCpuTime());
#endif

  return 0; // ok
}
