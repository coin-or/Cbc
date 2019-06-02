// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>

#include "CoinTime.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpFactorization.hpp"
#include "ClpObjective.hpp"
#include "ClpSimplex.hpp"
#include "CbcSolverLongThin.hpp"
#include "CbcModel.hpp"
#include "ClpPresolve.hpp"
#include "CbcBranchActual.hpp"
#include "CbcBranchFollow2.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcCompareUser.hpp"
// Cuts

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding.hpp"
#include "CglTwomir.hpp"
#include "CglDuplicateRow.hpp"
#include "CbcFathomDynamicProgramming.hpp"

static int timesBad_ = 0;
//#############################################################################
// Solve methods
//#############################################################################
static CglDuplicateRow *tryCut = NULL;
void CbcSolverLongThin::initialSolve()
{
  modelPtr_->scaling(0);
  setBasis(basis_, modelPtr_);
  modelPtr_->dual();
  basis_ = getBasis(modelPtr_);
  assert(!modelPtr_->specialOptions());
  modelPtr_->setLogLevel(0);
  if (!tryCut) {
    tryCut = new CglDuplicateRow(this);
    tryCut->setLogLevel(2);
  }
}

//-----------------------------------------------------------------------------
void CbcSolverLongThin::resolve()
{
  int *whichRow = NULL;
  int *whichColumn = NULL;
  // problem may be small enough to do nested search
  const double *colLower = modelPtr_->columnLower();
  const double *colUpper = modelPtr_->columnUpper();

  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  int numberRows = modelPtr_->numberRows();
  int numberColumns = modelPtr_->numberColumns();

  int i;
  int nFix = 0;
  int nNewRow = 0;
  int nNewCol = 0;
  int sizeDynamic = COIN_INT_MAX;
  int smallOriginalNumberRows = 0;
  if (algorithm_ == 0) {
    for (i = 0; i < numberIntegers; i++) {
      int iColumn = integerVariable[i];
      if (colLower[iColumn] == colUpper[iColumn])
        nFix++;
    }
  } else {
    whichRow = new int[numberRows];
    whichColumn = new int[numberColumns];
    // more sophisticated
    OsiCuts cs;
    tryCut->generateCuts(*this, cs);
    int numberCuts = cs.sizeColCuts();
    if (numberCuts) {
      for (i = 0; i < numberCuts; i++) {
        const OsiColCut *thisCut = cs.colCutPtr(i);
        const CoinPackedVector &ubs = thisCut->ubs();
        int n = ubs.getNumElements();
        const int *which = ubs.getIndices();
        const double *values = ubs.getElements();
        for (int j = 0; j < n; j++) {
          int iColumn = which[j];
          this->setColUpper(iColumn, values[j]);
        }
      }
    }
#if 1
    const int *duplicate = tryCut->duplicate();
    sizeDynamic = tryCut->sizeDynamic();
    int nOrig = tryCut->numberOriginalRows();
    for (i = 0; i < nOrig; i++) {
      if (duplicate[i] == -1)
        whichRow[nNewRow++] = i;
      else
        modelPtr_->setRowStatus(i, ClpSimplex::basic);
    }
    smallOriginalNumberRows = nNewRow;
    for (; i < numberRows; i++) {
      whichRow[nNewRow++] = i;
    }
#else
    for (i = 0; i < numberRows; i++)
      whichRow[i] = i;
    nNewRow = numberRows;
#endif
    for (i = 0; i < numberIntegers; i++) {
      int iColumn = integerVariable[i];
      if (colLower[iColumn] == colUpper[iColumn])
        nFix++;
      bool choose;
      if (algorithm_ == 1)
        choose = true;
      else
        choose = (node_[i] > count_ - memory_ && node_[i] > 0);
      if ((choose && colUpper[i])
        || (modelPtr_->getStatus(i) != ClpSimplex::atLowerBound && modelPtr_->getStatus(i) != ClpSimplex::isFixed)
        || colLower[i] > 0.0)
        whichColumn[nNewCol++] = i;
    }
  }
  if (nestedSearch_ < 1.0 && model_ && model_->phase() == 2) {
    if (((double)sizeDynamic) * ((double)nNewCol) < 1000000000 && sizeDynamic < 10000000) {
      // could do Dynamic Programming
      // back to original number of rows
      nNewRow = smallOriginalNumberRows;
      // and get rid of any basics
      int nNewCol = 0;
      for (i = 0; i < numberColumns; i++) {
        if (colUpper[i] || colLower[i] > 0.0)
          whichColumn[nNewCol++] = i;
      }
      ClpSimplex temp(modelPtr_, nNewRow, whichRow, nNewCol, whichColumn);
      int returnCode;
      double *rowLower2 = temp.rowLower();
      double *rowUpper2 = temp.rowUpper();
      int numberColumns2 = temp.numberColumns();
      double *colLower2 = temp.columnLower();
      double *colUpper2 = temp.columnUpper();
      const CoinPackedMatrix *matrix = temp.matrix();
      const double *element = matrix->getElements();
      const int *row = matrix->getIndices();
      const CoinBigIndex *columnStart = matrix->getVectorStarts();
      const int *columnLength = matrix->getVectorLengths();
      double offset = 0.0;
      const double *objective = temp.objective();
      bool feasible = true;
      for (i = 0; i < numberColumns2; i++) {
        double value = colLower2[i];
        if (value) {
          offset += value * objective[i];
          colLower2[i] = 0.0;
          colUpper2[i] -= value;
          for (int j = columnStart[i];
               j < columnStart[i] + columnLength[i]; j++) {
            int iRow = row[j];
            rowLower2[iRow] -= value * element[j];
            rowUpper2[iRow] -= value * element[j];
            if (rowUpper2[iRow] < -1.0e-8) {
              feasible = false;
              printf("odd - problem is infeasible\n");
            }
          }
        }
      }
      temp.setObjectiveOffset(-offset);
      OsiClpSolverInterface temp2(&temp);
      double *solutionDP = NULL;
      if (feasible) {
        for (i = 0; i < numberColumns2; i++)
          temp2.setInteger(i);
        CbcModel modelSmall(temp2);
        modelSmall.messageHandler()->setLogLevel(0);
        CbcFathomDynamicProgramming fathom1(modelSmall);
        // Set maximum space allowed
        fathom1.setMaximumSize(100000000);
        temp2.writeMps("small");
        returnCode = fathom1.fathom(solutionDP);
        if (returnCode != 1) {
          printf("probably not enough memory\n");
          abort();
        }
      }
      if (solutionDP) {
        double objValue = 0.0;
        double *solution = modelPtr_->primalColumnSolution();
        const double *objective = modelPtr_->objective();
        for (i = 0; i < numberColumns; i++)
          solution[i] = colLower[i];
        for (i = 0; i < nNewCol; i++) {
          int iColumn = whichColumn[i];
          solution[iColumn] += solutionDP[i];
        }
        for (i = 0; i < numberColumns; i++)
          objValue += solution[i] * objective[i];
        if (objValue < model_->getCutoff()) {
          printf("good solution %g by dynamic programming\n", objValue);
          returnCode = 0;
          // paranoid check
          double *rowLower = modelPtr_->rowLower();
          double *rowUpper = modelPtr_->rowUpper();
          // Column copy
          const CoinPackedMatrix *matrix2 = modelPtr_->matrix();
          element = matrix2->getElements();
          row = matrix2->getIndices();
          columnStart = matrix2->getVectorStarts();
          columnLength = matrix2->getVectorLengths();
          double *rowActivity = new double[numberRows];
          memset(rowActivity, 0, numberRows * sizeof(double));
          for (i = 0; i < numberColumns; i++) {
            int j;
            double value = solution[i];
            assert(value >= colLower[i] && value <= colUpper[i]);
            if (value) {
              printf("%d has value %g\n", i, value);
              for (j = columnStart[i];
                   j < columnStart[i] + columnLength[i]; j++) {
                int iRow = row[j];
                rowActivity[iRow] += value * element[j];
              }
            }
          }
          // check was feasible
          bool feasible = true;
          for (i = 0; i < numberRows; i++) {
            if (rowActivity[i] < rowLower[i]) {
              if (rowActivity[i] < rowLower[i] - 1.0e-8)
                feasible = false;
            } else if (rowActivity[i] > rowUpper[i]) {
              if (rowActivity[i] > rowUpper[i] + 1.0e-8)
                feasible = false;
            }
          }
          if (!feasible) {
            printf("** Bad solution by dynamic programming\n");
            abort();
          }
          delete[] rowActivity;
          model_->setBestSolution(CBC_TREE_SOL, objValue, solution);
        } else {
          returnCode = 2;
        }
      } else {
        returnCode = 2;
      }
      temp2.releaseClp();
      modelPtr_->setProblemStatus(1);
      delete[] whichRow;
      delete[] whichColumn;
      return;
    }
    if (nFix > nestedSearch_ * numberIntegers) {
      // Do nested search
      // back to original number of rows
      nNewRow = smallOriginalNumberRows;
      // and get rid of any basics
      int nNewCol = 0;
      for (i = 0; i < numberColumns; i++) {
        if (colUpper[i] || colLower[i] > 0.0)
          whichColumn[nNewCol++] = i;
      }
#if 0
      // We clone from continuous solver so set some stuff
      OsiSolverInterface * solver = model_->continuousSolver();
      CbcSolverLongThin * osiclp = dynamic_cast< CbcSolverLongThin*> (solver);
      assert (osiclp);
      // up special options
      if (osiclp->specialOptions()==3)
	osiclp->setSpecialOptions(7);
      double saveNested = osiclp->getNested();
      int saveAlgorithm = osiclp->getAlgorithm();
      osiclp->setNested(1.0);
      osiclp->setAlgorithm(0);
      int numberObjects = model_->numberObjects();
      if (numberObjects>model_->numberIntegers()) {
	// for now only integers
	//assert (numberObjects == model_->numberIntegers()+1);
	model_->setNumberObjects(model_->numberIntegers());
        // try follow on
	//model_->setNumberObjects(model_->numberIntegers()+1);
      }
      double saveMaxTime = model_->getDblParam(CbcModel::CbcMaximumSeconds);
      model_->setDblParam(CbcModel::CbcMaximumSeconds,1.0e5);
      // up special options
#if 1
      int returnCode= model_->subBranchAndBound(colLower,colUpper,2000);
#else
      CbcModel * model3 = model_->cleanModel(colLower,colUpper);
      // integer presolve
      int returnCode=0;
      CbcModel * model2 = model3->integerPresolve(false);
      if (!model2||!model2->getNumRows()) {
        delete model2;
        delete model3;
        returnCode= 2;
      } else {
        if (handler_->logLevel()>1)
          printf("Reduced model has %d rows and %d columns\n",
                 model2->getNumRows(),model2->getNumCols());
        if (true) {
          OsiSolverInterface * solver = model2->solver();
          OsiSolverInterface * osiclp = dynamic_cast< OsiSolverInterface*> (solver);
          assert (osiclp);
          int * priority = new int [numberColumns+1];
          int n=0;
          int iColumn;
          for ( iColumn=0;iColumn<numberColumns;iColumn++) {
            if (solver->isInteger(iColumn)) {
              priority[n++]=10000;
            }
          }
          priority[n]=1;
          CbcObject * newObject =new CbcFollowOn2(model2);
          model2->addObjects(1,&newObject);
          delete newObject;
          model2->passInPriorities(priority,false);
          delete [] priority;
        }
        returnCode= model_->subBranchAndBound(model3,model2,4000);
      }
#endif
      model_->setDblParam(CbcModel::CbcMaximumSeconds,saveMaxTime);
      model_->setNumberObjects(numberObjects);
      osiclp->setNested(saveNested);
      osiclp->setAlgorithm(saveAlgorithm);
#else
      // start again very simply
      ClpSimplex temp(modelPtr_, nNewRow, whichRow, nNewCol, whichColumn);
      int returnCode;
      OsiClpSolverInterface temp2(&temp);
      temp2.setupForRepeatedUse(2);
      int numberColumns2 = temp.numberColumns();
      const double *colUpper2 = temp2.getColUpper();
      const double *colLower2 = temp2.getColLower();
      const double *solution2 = temp.getColSolution();
      double *cleanSolution2 = new double[numberColumns2];
      for (i = 0; i < numberColumns2; i++) {
        temp2.setInteger(i);
        double value = solution2[i];
        value = CoinMin(CoinMax(value, colLower2[i]), colUpper2[i]);
        cleanSolution2[i] = value;
      }
      temp2.setColSolution(cleanSolution2);
      delete[] cleanSolution2;
      CbcModel modelSmall(temp2);
      modelSmall.setNumberStrong(0);
      CglProbing generator1;
      generator1.setUsingObjective(true);
      generator1.setMaxPass(3);
      generator1.setMaxProbe(100);
      generator1.setMaxLook(50);
      generator1.setRowCuts(3);

      CglGomory generator2;
      // try larger limit
      generator2.setLimit(300);

      CglKnapsackCover generator3;

      CglOddHole generator4;
      generator4.setMinimumViolation(0.005);
      generator4.setMinimumViolationPer(0.00002);
      // try larger limit
      generator4.setMaximumEntries(200);

      CglClique generator5;
      generator5.setStarCliqueReport(false);
      generator5.setRowCliqueReport(false);

      CglMixedIntegerRounding mixedGen;
      CglFlowCover flowGen;

      // Add in generators
      modelSmall.addCutGenerator(&generator1, -1, "Probing", true, false, false, -1);
      modelSmall.addCutGenerator(&generator2, -99, "Gomory", true, false, false, -99);
      modelSmall.addCutGenerator(&generator3, -99, "Knapsack", true, false, false, -99);
      modelSmall.addCutGenerator(&generator4, -99, "OddHole", true, false, false, -99);
      modelSmall.addCutGenerator(&generator5, -99, "Clique", true, false, false, -99);
      modelSmall.addCutGenerator(&flowGen, -99, "FlowCover", true, false, false, -99);
      modelSmall.addCutGenerator(&mixedGen, -99, "MixedIntegerRounding", true, false, false, -100);
#if 1
      const CoinPackedMatrix *matrix = temp2.getMatrixByCol();
      const int *columnLength = matrix->getVectorLengths();
      int *priority = new int[numberColumns2 + 1];
      // do pseudo costs and priorities - take a reasonable guess
      CbcObject **objects = new CbcObject *[numberColumns2 + 1];
      int n = 0;
      const double *objective = modelSmall.getObjCoefficients();
      for (i = 0; i < numberColumns2; i++) {
        CbcSimpleIntegerPseudoCost *newObject = new CbcSimpleIntegerPseudoCost(&modelSmall, n, i, objective[i], 0.5 * objective[i]);
        newObject->setMethod(3);
        objects[n] = newObject;
        priority[n++] = 10000 - columnLength[i];
      }
      priority[n] = 1;
      objects[n++] = new CbcFollowOn2(&modelSmall);
      modelSmall.addObjects(n, objects);
      for (i = 0; i < n; i++)
        delete objects[i];
      delete[] objects;
      modelSmall.passInPriorities(priority, false);
      delete[] priority;
#endif
      modelSmall.setCutoff(model_->getCutoff());
      //if (!onPathX&&modelSmall.getCutoff()>480.5)
      //modelSmall.setCutoff(480.5);
      //printf("cutoff %g\n",model_->getCutoff());
      modelSmall.messageHandler()->setLogLevel(1);
      modelSmall.solver()->messageHandler()->setLogLevel(0);
      modelSmall.messagesPointer()->setDetailMessage(3, 9);
      modelSmall.messagesPointer()->setDetailMessage(3, 6);
      modelSmall.messagesPointer()->setDetailMessage(3, 4);
      modelSmall.messagesPointer()->setDetailMessage(3, 13);
      modelSmall.messagesPointer()->setDetailMessage(3, 14);
      modelSmall.messagesPointer()->setDetailMessage(3, 1);
      modelSmall.messagesPointer()->setDetailMessage(3, 3007);
      modelSmall.branchAndBound();
      temp2.releaseClp();
      if (modelSmall.bestSolution()) {
        double objValue = 0.0;
        const double *solution2 = modelSmall.bestSolution();
        double *solution = modelPtr_->primalColumnSolution();
        const double *objective = modelPtr_->objective();
        for (i = 0; i < numberColumns; i++)
          solution[i] = colLower[i];
        for (i = 0; i < nNewCol; i++) {
          int iColumn = whichColumn[i];
          solution[iColumn] = solution2[i];
        }
        for (i = 0; i < numberColumns; i++)
          objValue += solution[i] * objective[i];
        assert(objValue < model_->getCutoff());
        if (objValue < model_->getCutoff()) {
          //printf("good solution \n");
          model_->setBestSolution(CBC_TREE_SOL, objValue, solution);
          returnCode = 0;
        } else {
          returnCode = 2;
        }
      } else {
        returnCode = 2;
      }
#endif
      if (returnCode != 0 && returnCode != 2) {
        printf("pretending entire search done\n");
        returnCode = 0;
      }
      if (returnCode == 0 || returnCode == 2) {
        modelPtr_->setProblemStatus(1);
        delete[] whichRow;
        delete[] whichColumn;
        return;
      }
    }
  }
  if ((count_ < 100 && algorithm_ == 2) || !algorithm_) {
    delete[] whichRow;
    delete[] whichColumn;
    assert(!modelPtr_->specialOptions());
    int saveOptions = modelPtr_->specialOptions();
    bool takeHint;
    OsiHintStrength strength;
    getHintParam(OsiDoInBranchAndCut, takeHint, strength);
    if (strength != OsiHintIgnore && takeHint) {
      // could do something - think about it
      //printf("thin hint %d %c\n",strength,takeHint ? 'T' :'F');
    }
    if ((specialOptions_ & 1) == 0) {
      modelPtr_->setSpecialOptions(saveOptions | (64 | 1024));
    } else {
      if ((specialOptions_ & 4) == 0)
        modelPtr_->setSpecialOptions(saveOptions | (64 | 128 | 512 | 1024 | 4096));
      else
        modelPtr_->setSpecialOptions(saveOptions | (64 | 128 | 512 | 1024 | 2048 | 4096));
    }
    //printf("thin options %d size %d\n",modelPtr_->specialOptions(),modelPtr_->numberColumns());
    setBasis(basis_, modelPtr_);
    //modelPtr_->setLogLevel(1);
    modelPtr_->dual(0, 0);
    basis_ = getBasis(modelPtr_);
    modelPtr_->setSpecialOptions(saveOptions);
    if (modelPtr_->status() == 0) {
      count_++;
      double *solution = modelPtr_->primalColumnSolution();
      int i;
      for (i = 0; i < numberColumns; i++) {
        if (solution[i] > 1.0e-6 || modelPtr_->getStatus(i) == ClpSimplex::basic) {
          node_[i] = CoinMax(count_, node_[i]);
          howMany_[i]++;
        }
      }
    } else {
      if (!algorithm_ == 2)
        printf("infeasible early on\n");
    }
  } else {
    // use counts
    int i;
    const double *lower = modelPtr_->columnLower();
    const double *upper = modelPtr_->columnUpper();
    setBasis(basis_, modelPtr_);
    ClpSimplex *temp = new ClpSimplex(modelPtr_, nNewRow, whichRow, nNewCol, whichColumn);
    //temp->setLogLevel(2);
    //printf("small has %d rows and %d columns\n",nNewRow,nNewCol);
    temp->setSpecialOptions(128 + 512);
    temp->setDualObjectiveLimit(1.0e50);
    temp->dual();
    if (temp->status()) {
      // In some cases we know that it must be infeasible
      if (believeInfeasible_ || algorithm_ == 1) {
        modelPtr_->setProblemStatus(1);
        printf("assuming infeasible!\n");
        //modelPtr_->writeMps("infeas.mps");
        //temp->writeMps("infeas2.mps");
        //abort();
        delete temp;
        delete[] whichRow;
        delete[] whichColumn;
        return;
      }
    }
    double *solution = modelPtr_->primalColumnSolution();
    if (!temp->status()) {
      const double *solution2 = temp->primalColumnSolution();
      memset(solution, 0, numberColumns * sizeof(double));
      for (i = 0; i < nNewCol; i++) {
        int iColumn = whichColumn[i];
        solution[iColumn] = solution2[i];
        modelPtr_->setStatus(iColumn, temp->getStatus(i));
      }
      double *rowSolution = modelPtr_->primalRowSolution();
      const double *rowSolution2 = temp->primalRowSolution();
      double *dual = modelPtr_->dualRowSolution();
      const double *dual2 = temp->dualRowSolution();
      memset(dual, 0, numberRows * sizeof(double));
      for (i = 0; i < nNewRow; i++) {
        int iRow = whichRow[i];
        modelPtr_->setRowStatus(iRow, temp->getRowStatus(i));
        rowSolution[iRow] = rowSolution2[i];
        dual[iRow] = dual2[i];
      }
      // See if optimal
      double *dj = modelPtr_->dualColumnSolution();
      // get reduced cost for large problem
      // this assumes minimization
      memcpy(dj, modelPtr_->objective(), numberColumns * sizeof(double));
      modelPtr_->transposeTimes(-1.0, dual, dj);
      modelPtr_->setObjectiveValue(temp->objectiveValue());
      modelPtr_->setProblemStatus(0);
      int nBad = 0;

      for (i = 0; i < numberColumns; i++) {
        if (modelPtr_->getStatus(i) == ClpSimplex::atLowerBound
          && upper[i] > lower[i] && dj[i] < -1.0e-5)
          nBad++;
      }
      //modelPtr_->writeMps("bada.mps");
      //temp->writeMps("badb.mps");
      if (nBad) {
        assert(algorithm_ == 2);
        //printf("%d bad\n",nBad);
        timesBad_++;
        modelPtr_->primal();
      }
    } else {
      // infeasible - do all
      modelPtr_->setSpecialOptions(64 + 128 + 512);
      setBasis(basis_, modelPtr_);
      //modelPtr_->setLogLevel(1);
      modelPtr_->dual(0, 0);
      basis_ = getBasis(modelPtr_);
      modelPtr_->setSpecialOptions(0);
      if (modelPtr_->status()) {
        printf("really infeasible!\n");
        delete temp;
        delete[] whichRow;
        delete[] whichColumn;
        return;
      } else {
        printf("initially infeasible\n");
      }
    }
    delete temp;
    delete[] whichRow;
    delete[] whichColumn;
    basis_ = getBasis(modelPtr_);
    modelPtr_->setSpecialOptions(0);
    count_++;
    if ((count_ % 100) == 0 && algorithm_ == 2)
      printf("count %d, bad %d\n", count_, timesBad_);
    for (i = 0; i < numberColumns; i++) {
      if (solution[i] > 1.0e-6 || modelPtr_->getStatus(i) == ClpSimplex::basic) {
        node_[i] = CoinMax(count_, node_[i]);
        howMany_[i]++;
      }
    }
    if (modelPtr_->objectiveValue() >= modelPtr_->dualObjectiveLimit())
      modelPtr_->setProblemStatus(1);
  }
}

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CbcSolverLongThin::CbcSolverLongThin()
  : OsiClpSolverInterface()
{
  node_ = NULL;
  howMany_ = NULL;
  count_ = 0;
  model_ = NULL;
  memory_ = 300;
  believeInfeasible_ = false;
  nestedSearch_ = 1.0;
  algorithm_ = 0;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface *
CbcSolverLongThin::clone(bool CopyData) const
{
  if (CopyData) {
    return new CbcSolverLongThin(*this);
  } else {
    printf("warning CbcSolveUser clone with copyData false\n");
    return new CbcSolverLongThin();
  }
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CbcSolverLongThin::CbcSolverLongThin(
  const CbcSolverLongThin &rhs)
  : OsiClpSolverInterface(rhs)
{
  model_ = rhs.model_;
  int numberColumns = modelPtr_->numberColumns();
  node_ = CoinCopyOfArray(rhs.node_, numberColumns);
  howMany_ = CoinCopyOfArray(rhs.howMany_, numberColumns);
  count_ = rhs.count_;
  memory_ = rhs.memory_;
  believeInfeasible_ = rhs.believeInfeasible_;
  nestedSearch_ = rhs.nestedSearch_;
  algorithm_ = rhs.algorithm_;
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CbcSolverLongThin::~CbcSolverLongThin()
{
  delete[] node_;
  delete[] howMany_;
}

//-------------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CbcSolverLongThin &
CbcSolverLongThin::operator=(const CbcSolverLongThin &rhs)
{
  if (this != &rhs) {
    OsiClpSolverInterface::operator=(rhs);
    delete[] node_;
    delete[] howMany_;
    model_ = rhs.model_;
    int numberColumns = modelPtr_->numberColumns();
    node_ = CoinCopyOfArray(rhs.node_, numberColumns);
    howMany_ = CoinCopyOfArray(rhs.howMany_, numberColumns);
    count_ = rhs.count_;
    memory_ = rhs.memory_;
    believeInfeasible_ = rhs.believeInfeasible_;
    nestedSearch_ = rhs.nestedSearch_;
    algorithm_ = rhs.algorithm_;
  }
  return *this;
}
//-------------------------------------------------------------------
// Real initializer
//-------------------------------------------------------------------
void CbcSolverLongThin::initialize(CbcModel *model, const char *keep)
{
  model_ = model;
  int numberColumns = modelPtr_->numberColumns();
  if (numberColumns) {
    node_ = new int[numberColumns];
    howMany_ = new int[numberColumns];
    for (int i = 0; i < numberColumns; i++) {
      if (keep[i])
        node_[i] = COIN_INT_MAX;
      else
        node_[i] = 0;
      howMany_[i] = 0;
    }
  } else {
    node_ = NULL;
    howMany_ = NULL;
  }
}
