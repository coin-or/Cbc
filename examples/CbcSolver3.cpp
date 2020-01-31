// $Id$
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>

#include "CoinTime.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpFactorization.hpp"
#include "ClpObjective.hpp"
#include "ClpSimplex.hpp"
#include "CbcSolver3.hpp"
#include "CbcModel.hpp"
#include "ClpPresolve.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcBranchActual.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcCompareUser.hpp"
#include "CbcCompareActual.hpp"
#include "CbcCompareObjective.hpp"
// Cuts

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding.hpp"
#include "CglTwomir.hpp"

static int timesBad_ = 0;
//#############################################################################
// Solve methods
//#############################################################################
void CbcSolver3::initialSolve()
{
  modelPtr_->scaling(0);
  setBasis(basis_, modelPtr_);
  // Do long thin by sprint
  ClpSolve options;
  options.setSolveType(ClpSolve::usePrimalorSprint);
  options.setPresolveType(ClpSolve::presolveOff);
  options.setSpecialOption(1, 3, 30);
  modelPtr_->initialSolve(options);
  basis_ = getBasis(modelPtr_);
  modelPtr_->setLogLevel(0);
}

//-----------------------------------------------------------------------------
void CbcSolver3::resolve()
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
    for (i = 0; i < numberRows; i++)
      whichRow[i] = i;
    nNewRow = numberRows;
    smallOriginalNumberRows = nNewRow;
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
      generator2.setLimit(3000);

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
      // Add in generators (just at root)
      modelSmall.addCutGenerator(&generator1, -99, "Probing", true, false, false, -1);
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
      modelSmall.addObjects(n, objects);
      for (i = 0; i < n; i++)
        delete objects[i];
      delete[] objects;
      modelSmall.passInPriorities(priority, false);
      delete[] priority;
#endif
      modelSmall.setCutoff(model_->getCutoff());
      modelSmall.messageHandler()->setLogLevel(1);
      modelSmall.solver()->messageHandler()->setLogLevel(0);
      modelSmall.messagesPointer()->setDetailMessage(3, 9);
      modelSmall.messagesPointer()->setDetailMessage(3, 6);
      modelSmall.messagesPointer()->setDetailMessage(3, 4);
      modelSmall.messagesPointer()->setDetailMessage(3, 13);
      modelSmall.messagesPointer()->setDetailMessage(3, 14);
      modelSmall.messagesPointer()->setDetailMessage(3, 1);
      modelSmall.messagesPointer()->setDetailMessage(3, 3007);
      // Definition of node choice
      CbcCompareObjective compare;
      modelSmall.setNodeComparison(compare);
      // And Greedy heuristic
      CbcHeuristicGreedyCover heuristic2(modelSmall);
      // Use original upper and perturb more
      heuristic2.setAlgorithm(11);
      modelSmall.addHeuristic(&heuristic2);
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
    printf("model has %d rows\n", modelPtr_->numberRows());
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
    printf("model has %d rows\n", nNewRow);
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
CbcSolver3::CbcSolver3()
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
CbcSolver3::clone(bool CopyData) const
{
  if (CopyData) {
    return new CbcSolver3(*this);
  } else {
    printf("warning CbcSolveUser clone with copyData false\n");
    return new CbcSolver3();
  }
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CbcSolver3::CbcSolver3(
  const CbcSolver3 &rhs)
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
CbcSolver3::~CbcSolver3()
{
  delete[] node_;
  delete[] howMany_;
}

//-------------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CbcSolver3 &
CbcSolver3::operator=(const CbcSolver3 &rhs)
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
void CbcSolver3::initialize(CbcModel *model, const char *keep)
{
  model_ = model;
  int numberColumns = modelPtr_->numberColumns();
  if (numberColumns) {
    node_ = new int[numberColumns];
    howMany_ = new int[numberColumns];
    for (int i = 0; i < numberColumns; i++) {
      if (keep && keep[i])
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
