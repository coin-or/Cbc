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
#include "ClpSolve.hpp"
#include "CbcSolver2.hpp"
#include "CbcModel.hpp"
#include "CoinModel.hpp"

static int timesBad_ = 0;
static int iterationsBad_ = 0;
//#############################################################################
// Solve methods
//#############################################################################
void CbcSolver2::initialSolve()
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
void CbcSolver2::resolve()
{
  int numberColumns = modelPtr_->numberColumns();
  if ((count_ < 10 && algorithm_ == 2) || !algorithm_) {
    OsiClpSolverInterface::resolve();
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
    int numberRows = modelPtr_->numberRows();
    int *whichRow = new int[numberRows];
    int *whichColumn = new int[numberColumns];
    int i;
    const double *lower = modelPtr_->columnLower();
    const double *upper = modelPtr_->columnUpper();
    const double *rowUpper = modelPtr_->rowUpper();
    bool equality = false;
    for (i = 0; i < numberRows; i++) {
      if (rowUpper[i] == 1.0) {
        equality = true;
        break;
      }
    }
    setBasis(basis_, modelPtr_);
    int nNewCol = 0;
    // Column copy
    //const double * element = modelPtr_->matrix()->getElements();
    const int *row = modelPtr_->matrix()->getIndices();
    const CoinBigIndex *columnStart = modelPtr_->matrix()->getVectorStarts();
    const int *columnLength = modelPtr_->matrix()->getVectorLengths();

    int *rowActivity = new int[numberRows];
    memset(rowActivity, 0, numberRows * sizeof(int));
    int *rowActivity2 = new int[numberRows];
    memset(rowActivity2, 0, numberRows * sizeof(int));
    char *mark = new char[numberColumns];
    memset(mark, 0, numberColumns);
    // Get rows which are satisfied
    for (i = 0; i < numberColumns; i++) {
      if (lower[i] > 0.0) {
        CoinBigIndex j;
        for (j = columnStart[i];
             j < columnStart[i] + columnLength[i]; j++) {
          int iRow = row[j];
          rowActivity2[iRow]++;
        }
      } else if (!upper[i]) {
        mark[i] = 2; // no good
      }
    }
    // If equality - check not infeasible
    if (equality) {
      bool feasible = true;
      for (i = 0; i < numberRows; i++) {
        if (rowActivity2[i] > 1) {
          feasible = false;
          break;
        }
      }
      if (!feasible) {
        delete[] rowActivity;
        delete[] rowActivity2;
        modelPtr_->setProblemStatus(1);
        delete[] whichRow;
        delete[] whichColumn;
        delete[] mark;
        printf("infeasible by inspection (over)\n");
        return;
      }
    }
    int nNoGood = 0;
    for (i = 0; i < numberColumns; i++) {
      if (mark[i] == 2) {
        nNoGood++;
        continue;
      }
      bool choose;
      if (algorithm_ == 1)
        choose = true;
      else
        choose = (node_[i] > count_ - memory_ && node_[i] > 0);
      bool any;
      if (equality) {
        // See if forced to be zero
        CoinBigIndex j;
        any = true;
        for (j = columnStart[i];
             j < columnStart[i] + columnLength[i]; j++) {
          int iRow = row[j];
          if (rowActivity2[iRow])
            any = false; // can't be in
        }
      } else {
        // See if not useful
        CoinBigIndex j;
        any = false;
        for (j = columnStart[i];
             j < columnStart[i] + columnLength[i]; j++) {
          int iRow = row[j];
          if (!rowActivity2[iRow])
            any = true; // useful
        }
      }
      if (!any && !lower[i]) {
        choose = false;
        // and say can't be useful
        mark[i] = 2;
        nNoGood++;
      }
      if (strategy_ && modelPtr_->getColumnStatus(i) == ClpSimplex::basic)
        choose = true;
      if (choose || lower[i] > 0.0) {
        mark[i] = 1;
        whichColumn[nNewCol++] = i;
        CoinBigIndex j;
        double value = upper[i];
        if (value) {
          for (j = columnStart[i];
               j < columnStart[i] + columnLength[i]; j++) {
            int iRow = row[j];
            rowActivity[iRow]++;
          }
        }
      }
    }
    // If equality add in slacks
    CoinModel build;
    if (equality) {
      int row = 0;
      for (i = 0; i < numberRows; i++) {
        // put in all rows if wanted
        if (strategy_)
          rowActivity2[i] = 0;
        if (!rowActivity2[i]) {
          double element = 1.0;
          build.addColumn(1, &row, &element, 0.0, 1.0, 1.0e8); // large cost
          row++;
        }
      }
    }
    int nOK = 0;
    int nNewRow = 0;
    for (i = 0; i < numberRows; i++) {
      if (rowActivity[i])
        nOK++;
      if (!rowActivity2[i])
        whichRow[nNewRow++] = i; // not satisfied
      else
        modelPtr_->setRowStatus(i, ClpSimplex::basic); // make slack basic
    }
    if (nOK < numberRows) {
      for (i = 0; i < numberColumns; i++) {
        if (!mark[i]) {
          CoinBigIndex j;
          int good = 0;
          for (j = columnStart[i];
               j < columnStart[i] + columnLength[i]; j++) {
            int iRow = row[j];
            if (!rowActivity[iRow]) {
              rowActivity[iRow]++;
              good++;
            }
          }
          if (good) {
            nOK += good;
            whichColumn[nNewCol++] = i;
          }
        }
      }
    }
    delete[] rowActivity;
    delete[] rowActivity2;
    if (nOK < numberRows) {
      modelPtr_->setProblemStatus(1);
      delete[] whichRow;
      delete[] whichColumn;
      delete[] mark;
      printf("infeasible by inspection\n");
      return;
    }
    bool allIn = false;
    if (nNewCol + nNoGood + numberRows > numberColumns) {
      // add in all
      allIn = true;
      for (i = 0; i < numberColumns; i++) {
        if (!mark[i]) {
          whichColumn[nNewCol++] = i;
        }
      }
    }
    ClpSimplex *temp = new ClpSimplex(modelPtr_, nNewRow, whichRow, nNewCol, whichColumn);
    if (equality)
      temp->addColumns(build);
    temp->setLogLevel(1);
    printf("small has %d rows and %d columns (%d impossible to help) %s\n",
      nNewRow, nNewCol, nNoGood, allIn ? "all in" : "");
    temp->setSpecialOptions(128 + 512);
    temp->setDualObjectiveLimit(1.0e50);
    temp->dual();
    assert(!temp->status());
    double *solution = modelPtr_->primalColumnSolution();
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
    delete temp;
    if (nBad && !allIn) {
      assert(algorithm_ == 2);
      //printf("%d bad\n",nBad);
      timesBad_++;
      // just non mark==2
      int nAdded = 0;
      for (i = 0; i < numberColumns; i++) {
        if (!mark[i]) {
          whichColumn[nNewCol++] = i;
          nAdded++;
        }
      }
      assert(nAdded);
      {
        temp = new ClpSimplex(modelPtr_, nNewRow, whichRow, nNewCol, whichColumn);
        if (equality)
          temp->addColumns(build);
        temp->setLogLevel(2);
        temp->setSpecialOptions(128 + 512);
        temp->setDualObjectiveLimit(1.0e50);
        temp->primal(1);
        assert(!temp->status());
        double *solution = modelPtr_->primalColumnSolution();
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
        modelPtr_->setObjectiveValue(temp->objectiveValue());
        modelPtr_->setProblemStatus(0);
        iterationsBad_ += temp->numberIterations();
        printf("clean %d\n", temp->numberIterations());
        delete temp;
      }
    }
    delete[] mark;
    delete[] whichRow;
    delete[] whichColumn;
    basis_ = getBasis(modelPtr_);
    modelPtr_->setSpecialOptions(0);
    count_++;
    if ((count_ % 100) == 0 && algorithm_ == 2)
      printf("count %d, bad %d - iterations %d\n", count_, timesBad_, iterationsBad_);
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
CbcSolver2::CbcSolver2()
  : OsiClpSolverInterface()
{
  node_ = NULL;
  howMany_ = NULL;
  count_ = 0;
  model_ = NULL;
  memory_ = 300;
  algorithm_ = 0;
  strategy_ = 0;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface *
CbcSolver2::clone(bool CopyData) const
{
  if (CopyData) {
    return new CbcSolver2(*this);
  } else {
    printf("warning CbcSolveUser clone with copyData false\n");
    return new CbcSolver2();
  }
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CbcSolver2::CbcSolver2(
  const CbcSolver2 &rhs)
  : OsiClpSolverInterface(rhs)
{
  model_ = rhs.model_;
  int numberColumns = modelPtr_->numberColumns();
  node_ = CoinCopyOfArray(rhs.node_, numberColumns);
  howMany_ = CoinCopyOfArray(rhs.howMany_, numberColumns);
  count_ = rhs.count_;
  memory_ = rhs.memory_;
  algorithm_ = rhs.algorithm_;
  strategy_ = rhs.strategy_;
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CbcSolver2::~CbcSolver2()
{
  delete[] node_;
  delete[] howMany_;
}

//-------------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CbcSolver2 &
CbcSolver2::operator=(const CbcSolver2 &rhs)
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
    algorithm_ = rhs.algorithm_;
    strategy_ = rhs.strategy_;
  }
  return *this;
}
//-------------------------------------------------------------------
// Real initializer
//-------------------------------------------------------------------
void CbcSolver2::initialize(CbcModel *model, const char *keep)
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
