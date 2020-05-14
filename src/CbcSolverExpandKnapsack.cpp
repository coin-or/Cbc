// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file CbcSolverExpandKnapsack.cpp

  Returns OsiSolverInterface (User should delete)
  On entry numberKnapsack is maximum number of Total entries

  Expanding possibilities of x*y, where x*y are both integers, constructing
  a knapsack constraint. Results in a tighter model.

*/

#include "CbcConfig.h"
#include "CoinPragma.hpp"

#include "OsiSolverInterface.hpp"

#include "CglStored.hpp"

#ifndef COIN_HAS_LINK
#define COIN_HAS_LINK
#endif
#ifdef COIN_HAS_LINK
#include "CbcLinked.hpp"
#endif

#ifdef COIN_HAS_LINK


/* Expands out all possible combinations for a knapsack
   If buildObj NULL then just computes space needed - returns number elements
   On entry numberOutput is maximum allowed, on exit it is number needed or
   -1 (as will be number elements) if maximum exceeded.  numberOutput will have at
   least space to return values which reconstruct input.
   Rows returned will be original rows but no entries will be returned for
   any rows all of whose entries are in knapsack.  So up to user to allow for this.
   If reConstruct >=0 then returns number of entrie which make up item "reConstruct"
   in expanded knapsack.  Values in buildRow and buildElement;
*/
int expandKnapsack(CoinModel& cm, int knapsackRow, int &numberOutput, double *buildObj, CoinBigIndex *buildStart,
  int *buildRow, double *buildElement, int reConstruct = -1)
{
  /* mark rows
       -2 in knapsack and other variables
       -1 not involved
       0 only in knapsack
    */
  int *markRow = new int[cm.numberRows_];
  int iRow;
  int iColumn;
  int *whichColumn = new int[cm.numberColumns_];
  for (iColumn = 0; iColumn < cm.numberColumns_; iColumn++)
    whichColumn[iColumn] = -1;
  int numJ = 0;
  for (iRow = 0; iRow < cm.numberRows_; iRow++)
    markRow[iRow] = -1;
  CoinModelLink triple;
  triple = cm.firstInRow(knapsackRow);
  while (triple.column() >= 0) {
    int iColumn = triple.column();
#ifndef NDEBUG
    const char *el = cm.getElementAsString(knapsackRow, iColumn);
    assert(!strcmp("Numeric", el));
#endif
    whichColumn[iColumn] = numJ;
    numJ++;
    triple = cm.next(triple);
  }
  for (iRow = 0; iRow < cm.numberRows_; iRow++) {
    triple = cm.firstInRow(iRow);
    int type = -3;
    while (triple.column() >= 0) {
      int iColumn = triple.column();
      if (whichColumn[iColumn] >= 0) {
        if (type == -3)
          type = 0;
        else if (type != 0)
          type = -2;
      } else {
        if (type == -3)
          type = -1;
        else if (type == 0)
          type = -2;
      }
      triple = cm.next(triple);
    }
    if (type == -3)
      type = -1;
    markRow[iRow] = type;
  }
  int *bound = new int[cm.numberColumns_ + 1];
  int *whichRow = new int[cm.numberRows_];
  ClpSimplex tempModel;
  CoinModel tempModel2(cm);
  tempModel.loadProblem(tempModel2);
  int *stack = new int[cm.numberColumns_ + 1];
  double *size = new double[cm.numberColumns_ + 1];
  double *rhsOffset = new double[cm.numberRows_];
  int *build = new int[cm.numberColumns_];
  int maxNumber = numberOutput;
  numJ = 0;
  double minSize = cm.getRowLower(knapsackRow);
  double maxSize = cm.getRowUpper(knapsackRow);
  double offset = 0.0;
  triple = cm.firstInRow(knapsackRow);
  while (triple.column() >= 0) {
    iColumn = triple.column();
    double lowerColumn = cm.columnLower(iColumn);
    double upperColumn = cm.columnUpper(iColumn);
    double gap = upperColumn - lowerColumn;
    if (gap > 1.0e8)
      gap = 1.0e8;
    assert(fabs(floor(gap + 0.5) - gap) < 1.0e-5);
    whichColumn[numJ] = iColumn;
    bound[numJ] = static_cast< int >(gap);
    size[numJ++] = triple.value();
    offset += triple.value() * lowerColumn;
    triple = cm.next(triple);
  }
  int jRow;
  for (iRow = 0; iRow < cm.numberRows_; iRow++)
    whichRow[iRow] = iRow;
  ClpSimplex smallModel(&tempModel, cm.numberRows_, whichRow, numJ, whichColumn, true, true, true);
  // modify rhs to allow for nonzero lower bounds
  double *rowLower = smallModel.rowLower();
  double *rowUpper = smallModel.rowUpper();
  const double *columnLower = smallModel.columnLower();
  //const double * columnUpper = smallModel.columnUpper();
  const CoinPackedMatrix *matrix = smallModel.matrix();
  const double *element = matrix->getElements();
  const int *row = matrix->getIndices();
  const CoinBigIndex *columnStart = matrix->getVectorStarts();
  const int *columnLength = matrix->getVectorLengths();
  const double *objective = smallModel.objective();
  double objectiveOffset = 0.0;
  CoinZeroN(rhsOffset, cm.numberRows_);
  for (iColumn = 0; iColumn < numJ; iColumn++) {
    double lower = columnLower[iColumn];
    if (lower) {
      objectiveOffset += objective[iColumn];
      for (CoinBigIndex j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        double value = element[j] * lower;
        int kRow = row[j];
        rhsOffset[kRow] += value;
        if (rowLower[kRow] > -1.0e20)
          rowLower[kRow] -= value;
        if (rowUpper[kRow] < 1.0e20)
          rowUpper[kRow] -= value;
      }
    }
  }
  // relax
  for (jRow = 0; jRow < cm.numberRows_; jRow++) {
    if (markRow[jRow] == 0 && knapsackRow != jRow) {
      if (rowLower[jRow] > -1.0e20)
        rowLower[jRow] -= 1.0e-7;
      if (rowUpper[jRow] < 1.0e20)
        rowUpper[jRow] += 1.0e-7;
    } else {
      rowLower[jRow] = -COIN_DBL_MAX;
      rowUpper[jRow] = COIN_DBL_MAX;
    }
  }
  double *rowActivity = smallModel.primalRowSolution();
  CoinZeroN(rowActivity, cm.numberRows_);
  maxSize -= offset;
  minSize -= offset;
  // now generate
  int i;
  int iStack = numJ;
  for (i = 0; i < numJ; i++) {
    stack[i] = 0;
  }
  double tooMuch = 10.0 * maxSize;
  stack[numJ] = 1;
  size[numJ] = tooMuch;
  bound[numJ] = 0;
  double sum = tooMuch;
  numberOutput = 0;
  int nelCreate = 0;
  /* typeRun is - 0 for initial sizes
                    1 for build
      	  2 for reconstruct
    */
  int typeRun = buildObj ? 1 : 0;
  if (reConstruct >= 0) {
    assert(buildRow && buildElement);
    typeRun = 2;
  }
  if (typeRun == 1)
    buildStart[0] = 0;
  while (iStack >= 0) {
    if (sum >= minSize && sum <= maxSize) {
      double checkSize = 0.0;
      bool good = true;
      int nRow = 0;
      double obj = objectiveOffset;
      // nRow is zero? CoinZeroN(rowActivity,nRow);
      for (iColumn = 0; iColumn < numJ; iColumn++) {
        int iValue = stack[iColumn];
        if (iValue > bound[iColumn]) {
          good = false;
          break;
        } else if (iValue) {
          obj += objective[iColumn] * iValue;
          for (CoinBigIndex j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            double value = element[j] * iValue;
            int kRow = row[j];
            if (rowActivity[kRow]) {
              rowActivity[kRow] += value;
              if (!rowActivity[kRow])
                rowActivity[kRow] = 1.0e-100;
            } else {
              build[nRow++] = kRow;
              rowActivity[kRow] = value;
            }
          }
        }
      }
      if (good) {
        for (jRow = 0; jRow < nRow; jRow++) {
          int kRow = build[jRow];
          double value = rowActivity[kRow];
          if (value > rowUpper[kRow] || value < rowLower[kRow]) {
            good = false;
            break;
          }
        }
      }
      if (good) {
        if (typeRun == 1) {
          buildObj[numberOutput] = obj;
          for (jRow = 0; jRow < nRow; jRow++) {
            int kRow = build[jRow];
            double value = rowActivity[kRow];
            if (markRow[kRow] < 0 && fabs(value) > 1.0e-13) {
              buildElement[nelCreate] = value;
              buildRow[nelCreate++] = kRow;
            }
          }
          buildStart[numberOutput + 1] = nelCreate;
        } else if (!typeRun) {
          for (jRow = 0; jRow < nRow; jRow++) {
            int kRow = build[jRow];
            double value = rowActivity[kRow];
            if (markRow[kRow] < 0 && fabs(value) > 1.0e-13) {
              nelCreate++;
            }
          }
        }
        if (typeRun == 2 && reConstruct == numberOutput) {
          // build and exit
          nelCreate = 0;
          for (iColumn = 0; iColumn < numJ; iColumn++) {
            int iValue = stack[iColumn];
            if (iValue) {
              buildRow[nelCreate] = whichColumn[iColumn];
              buildElement[nelCreate++] = iValue;
            }
          }
          numberOutput = 1;
          for (i = 0; i < numJ; i++) {
            bound[i] = 0;
          }
          break;
        }
        numberOutput++;
        if (numberOutput > maxNumber) {
          nelCreate = -1;
          numberOutput = -1;
          for (i = 0; i < numJ; i++) {
            bound[i] = 0;
          }
          break;
        } else if (typeRun == 1 && numberOutput == maxNumber) {
          // On second run
          for (i = 0; i < numJ; i++) {
            bound[i] = 0;
          }
          break;
        }
        for (int j = 0; j < numJ; j++) {
          checkSize += stack[j] * size[j];
        }
        assert(fabs(sum - checkSize) < 1.0e-3);
      }
      for (jRow = 0; jRow < nRow; jRow++) {
        int kRow = build[jRow];
        rowActivity[kRow] = 0.0;
      }
    }
    if (sum > maxSize || stack[iStack] > bound[iStack]) {
      sum -= size[iStack] * stack[iStack];
      stack[iStack--] = 0;
      if (iStack >= 0) {
        stack[iStack]++;
        sum += size[iStack];
      }
    } else {
      // must be less
      // add to last possible
      iStack = numJ - 1;
      sum += size[iStack];
      stack[iStack]++;
    }
  }
  //printf("%d will be created\n",numberOutput);
  delete[] whichColumn;
  delete[] whichRow;
  delete[] bound;
  delete[] stack;
  delete[] size;
  delete[] rhsOffset;
  delete[] build;
  delete[] markRow;
  return nelCreate;
}

OsiSolverInterface *
expandKnapsack(CoinModel &model, int *whichColumn, int *knapsackStart,
  int *knapsackRow, int &numberKnapsack,
  CglStored &stored, int logLevel,
  int fixedPriority, int SOSPriority, CoinModel &tightenedModel)
{
  int maxTotal = numberKnapsack;
  // load from coin model
  OsiSolverLink *si = new OsiSolverLink();
  OsiSolverInterface *finalModel = NULL;
  si->setDefaultMeshSize(0.001);
  // need some relative granularity
  si->setDefaultBound(100.0);
  si->setDefaultMeshSize(0.01);
  si->setDefaultBound(100000.0);
  si->setIntegerPriority(1000);
  si->setBiLinearPriority(10000);
  si->load(model, true, logLevel);
  // get priorities
  const int *priorities = model.priorities();
  int numberColumns = model.numberColumns();
  if (priorities) {
    OsiObject **objects = si->objects();
    int numberObjects = si->numberObjects();
    for (int iObj = 0; iObj < numberObjects; iObj++) {
      int iColumn = objects[iObj]->columnNumber();
      if (iColumn >= 0 && iColumn < numberColumns) {
#ifndef NDEBUG
        OsiSimpleInteger *obj = dynamic_cast< OsiSimpleInteger * >(objects[iObj]);
#endif
        assert(obj);
        int iPriority = priorities[iColumn];
        if (iPriority > 0)
          objects[iObj]->setPriority(iPriority);
      }
    }
    if (fixedPriority > 0) {
      si->setFixedPriority(fixedPriority);
    }
    if (SOSPriority < 0)
      SOSPriority = 100000;
  }
  CoinModel coinModel = *si->coinModel();
  assert(coinModel.numberRows() > 0);
  tightenedModel = coinModel;
  int numberRows = coinModel.numberRows();
  // Mark variables
  int *whichKnapsack = new int[numberColumns];
  int iRow, iColumn;
  for (iColumn = 0; iColumn < numberColumns; iColumn++)
    whichKnapsack[iColumn] = -1;
  int kRow;
  bool badModel = false;
  // analyze
  if (logLevel > 1) {
    for (iRow = 0; iRow < numberRows; iRow++) {
      /* Just obvious one at first
            positive non unit coefficients
            all integer
            positive rowUpper
            for now - linear (but further down in code may use nonlinear)
            column bounds should be tight
            */
      //double lower = coinModel.getRowLower(iRow);
      double upper = coinModel.getRowUpper(iRow);
      if (upper < 1.0e10) {
        CoinModelLink triple = coinModel.firstInRow(iRow);
        bool possible = true;
        int n = 0;
        int n1 = 0;
        while (triple.column() >= 0) {
          int iColumn = triple.column();
          const char *el = coinModel.getElementAsString(iRow, iColumn);
          if (!strcmp("Numeric", el)) {
            if (coinModel.columnLower(iColumn) == coinModel.columnUpper(iColumn)) {
              triple = coinModel.next(triple);
              continue; // fixed
            }
            double value = coinModel.getElement(iRow, iColumn);
            if (value < 0.0) {
              possible = false;
            } else {
              n++;
              if (value == 1.0)
                n1++;
              if (coinModel.columnLower(iColumn) < 0.0)
                possible = false;
              if (!coinModel.isInteger(iColumn))
                possible = false;
              if (whichKnapsack[iColumn] >= 0)
                possible = false;
            }
          } else {
            possible = false; // non linear
          }
          triple = coinModel.next(triple);
        }
        if (n - n1 > 1 && possible) {
          double lower = coinModel.getRowLower(iRow);
          double upper = coinModel.getRowUpper(iRow);
          CoinModelLink triple = coinModel.firstInRow(iRow);
          while (triple.column() >= 0) {
            int iColumn = triple.column();
            lower -= coinModel.columnLower(iColumn) * triple.value();
            upper -= coinModel.columnLower(iColumn) * triple.value();
            triple = coinModel.next(triple);
          }
          printf("%d is possible %g <=", iRow, lower);
          // print
          triple = coinModel.firstInRow(iRow);
          while (triple.column() >= 0) {
            int iColumn = triple.column();
            if (coinModel.columnLower(iColumn) != coinModel.columnUpper(iColumn))
              printf(" (%d,el %g up %g)", iColumn, triple.value(),
                coinModel.columnUpper(iColumn) - coinModel.columnLower(iColumn));
            triple = coinModel.next(triple);
          }
          printf(" <= %g\n", upper);
        }
      }
    }
  }
  numberKnapsack = 0;
  for (kRow = 0; kRow < numberRows; kRow++) {
    iRow = kRow;
    /* Just obvious one at first
           positive non unit coefficients
           all integer
           positive rowUpper
           for now - linear (but further down in code may use nonlinear)
           column bounds should be tight
        */
    //double lower = coinModel.getRowLower(iRow);
    double upper = coinModel.getRowUpper(iRow);
    if (upper < 1.0e10) {
      CoinModelLink triple = coinModel.firstInRow(iRow);
      bool possible = true;
      int n = 0;
      int n1 = 0;
      while (triple.column() >= 0) {
        int iColumn = triple.column();
        const char *el = coinModel.getElementAsString(iRow, iColumn);
        if (!strcmp("Numeric", el)) {
          if (coinModel.columnLower(iColumn) == coinModel.columnUpper(iColumn)) {
            triple = coinModel.next(triple);
            continue; // fixed
          }
          double value = coinModel.getElement(iRow, iColumn);
          if (value < 0.0) {
            possible = false;
          } else {
            n++;
            if (value == 1.0)
              n1++;
            if (coinModel.columnLower(iColumn) < 0.0)
              possible = false;
            if (!coinModel.isInteger(iColumn))
              possible = false;
            if (whichKnapsack[iColumn] >= 0)
              possible = false;
          }
        } else {
          possible = false; // non linear
        }
        triple = coinModel.next(triple);
      }
      if (n - n1 > 1 && possible) {
        // try
        CoinModelLink triple = coinModel.firstInRow(iRow);
        while (triple.column() >= 0) {
          int iColumn = triple.column();
          if (coinModel.columnLower(iColumn) != coinModel.columnUpper(iColumn))
            whichKnapsack[iColumn] = numberKnapsack;
          triple = coinModel.next(triple);
        }
        knapsackRow[numberKnapsack++] = iRow;
      }
    }
  }
  if (logLevel > 0)
    printf("%d out of %d candidate rows are possible\n", numberKnapsack, numberRows);
  // Check whether we can get rid of nonlinearities
  /* mark rows
       -2 in knapsack and other variables
       -1 not involved
       n only in knapsack n
    */
  int *markRow = new int[numberRows];
  for (iRow = 0; iRow < numberRows; iRow++)
    markRow[iRow] = -1;
  int canDo = 1; // OK and linear
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    CoinModelLink triple = coinModel.firstInColumn(iColumn);
    int iKnapsack = whichKnapsack[iColumn];
    bool linear = true;
    // See if quadratic objective
    const char *expr = coinModel.getColumnObjectiveAsString(iColumn);
    if (strcmp(expr, "Numeric")) {
      linear = false;
    }
    while (triple.row() >= 0) {
      int iRow = triple.row();
      if (iKnapsack >= 0) {
        if (markRow[iRow] == -1) {
          markRow[iRow] = iKnapsack;
        } else if (markRow[iRow] != iKnapsack) {
          markRow[iRow] = -2;
        }
      }
      const char *expr = coinModel.getElementAsString(iRow, iColumn);
      if (strcmp(expr, "Numeric")) {
        linear = false;
      }
      triple = coinModel.next(triple);
    }
    if (!linear) {
      if (whichKnapsack[iColumn] < 0) {
        canDo = 0;
        break;
      } else {
        canDo = 2;
      }
    }
  }
  int *markKnapsack = NULL;
  double *coefficient = NULL;
  double *linear = NULL;
  int *whichRow = NULL;
  int *lookupRow = NULL;
  badModel = (canDo == 0);
  if (numberKnapsack && canDo) {
    /* double check - OK if
           no nonlinear
           nonlinear only on columns in knapsack
           nonlinear only on columns in knapsack * ONE other - same for all in knapsack
           AND that is only row connected to knapsack
           (theoretically could split knapsack if two other and small numbers)
           also ONE could be ONE expression - not just a variable
        */
    int iKnapsack;
    markKnapsack = new int[numberKnapsack];
    coefficient = new double[numberKnapsack];
    linear = new double[numberColumns];
    for (iKnapsack = 0; iKnapsack < numberKnapsack; iKnapsack++)
      markKnapsack[iKnapsack] = -1;
    if (canDo == 2) {
      for (iRow = -1; iRow < numberRows; iRow++) {
        int numberOdd;
        CoinPackedMatrix *row = coinModel.quadraticRow(iRow, linear, numberOdd);
        if (row) {
          // see if valid
          const double *element = row->getElements();
          const int *column = row->getIndices();
          const CoinBigIndex *columnStart = row->getVectorStarts();
          const int *columnLength = row->getVectorLengths();
          int numberLook = row->getNumCols();
          for (int i = 0; i < numberLook; i++) {
            int iKnapsack = whichKnapsack[i];
            if (iKnapsack < 0) {
              // might be able to swap - but for now can't have knapsack in
              for (CoinBigIndex j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
                int iColumn = column[j];
                if (whichKnapsack[iColumn] >= 0) {
                  canDo = 0; // no good
                  badModel = true;
                  break;
                }
              }
            } else {
              // OK if in same knapsack - or maybe just one
              int marked = markKnapsack[iKnapsack];
              for (CoinBigIndex j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
                int iColumn = column[j];
                if (whichKnapsack[iColumn] != iKnapsack && whichKnapsack[iColumn] >= 0) {
                  canDo = 0; // no good
                  badModel = true;
                  break;
                } else if (marked == -1) {
                  markKnapsack[iKnapsack] = iColumn;
                  marked = iColumn;
                  coefficient[iKnapsack] = element[j];
                  coinModel.associateElement(coinModel.columnName(iColumn), 1.0);
                } else if (marked != iColumn) {
                  badModel = true;
                  canDo = 0; // no good
                  break;
                } else {
                  // could manage with different coefficients - but for now ...
                  assert(coefficient[iKnapsack] == element[j]);
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
      whichRow = new int[numberRows];
      lookupRow = new int[numberRows];
      bool someNonlinear = false;
      double maxCoefficient = 1.0;
      for (iKnapsack = 0; iKnapsack < numberKnapsack; iKnapsack++) {
        if (markKnapsack[iKnapsack] >= 0) {
          someNonlinear = true;
          int iColumn = markKnapsack[iKnapsack];
          maxCoefficient = CoinMax(maxCoefficient, fabs(coefficient[iKnapsack] * coinModel.columnUpper(iColumn)));
        }
      }
      if (someNonlinear) {
        // associate all columns to stop possible error messages
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          coinModel.associateElement(coinModel.columnName(iColumn), 1.0);
        }
      }
      ClpSimplex tempModel;
      tempModel.loadProblem(coinModel);
      // Create final model - first without knapsacks
      int nCol = 0;
      int nRow = 0;
      for (iRow = 0; iRow < numberRows; iRow++) {
        if (markRow[iRow] < 0) {
          lookupRow[iRow] = nRow;
          whichRow[nRow++] = iRow;
        } else {
          lookupRow[iRow] = -1;
        }
      }
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (whichKnapsack[iColumn] < 0)
          whichColumn[nCol++] = iColumn;
      }
      ClpSimplex finalModelX(&tempModel, nRow, whichRow, nCol, whichColumn, false, false, false);
      OsiClpSolverInterface finalModelY(&finalModelX, true);
      finalModel = finalModelY.clone();
      finalModelY.releaseClp();
      // Put back priorities
      const int *priorities = model.priorities();
      if (priorities) {
        finalModel->findIntegers(false);
        OsiObject **objects = finalModel->objects();
        int numberObjects = finalModel->numberObjects();
        for (int iObj = 0; iObj < numberObjects; iObj++) {
          int iColumn = objects[iObj]->columnNumber();
          if (iColumn >= 0 && iColumn < nCol) {
#ifndef NDEBUG
            OsiSimpleInteger *obj = dynamic_cast< OsiSimpleInteger * >(objects[iObj]);
#endif
            assert(obj);
            int iPriority = priorities[whichColumn[iColumn]];
            if (iPriority > 0)
              objects[iObj]->setPriority(iPriority);
          }
        }
      }
      for (iRow = 0; iRow < numberRows; iRow++) {
        whichRow[iRow] = iRow;
      }
      int numberOther = finalModel->getNumCols();
      int nLargest = 0;
      int nelLargest = 0;
      int nTotal = 0;
      for (iKnapsack = 0; iKnapsack < numberKnapsack; iKnapsack++) {
        iRow = knapsackRow[iKnapsack];
        int nCreate = maxTotal;
        int nelCreate = expandKnapsack(coinModel, iRow, nCreate, NULL, NULL, NULL, NULL);
        if (nelCreate < 0)
          badModel = true;
        nTotal += nCreate;
        nLargest = CoinMax(nLargest, nCreate);
        nelLargest = CoinMax(nelLargest, nelCreate);
      }
      if (nTotal > maxTotal)
        badModel = true;
      if (!badModel) {
        // Now arrays for building
        nelLargest = CoinMax(nelLargest, nLargest) + 1;
        double *buildObj = new double[nLargest];
        double *buildElement = new double[nelLargest];
        CoinBigIndex *buildStart = new CoinBigIndex[nLargest + 1];
        int *buildRow = new int[nelLargest];
        // alow for integers in knapsacks
        OsiObject **object = new OsiObject *[numberKnapsack + nTotal];
        int nSOS = 0;
        int nObj = numberKnapsack;
        for (iKnapsack = 0; iKnapsack < numberKnapsack; iKnapsack++) {
          knapsackStart[iKnapsack] = finalModel->getNumCols();
          iRow = knapsackRow[iKnapsack];
          int nCreate = 10000;
          expandKnapsack(coinModel, iRow, nCreate, buildObj, buildStart, buildRow, buildElement);
          // Redo row numbers
          for (iColumn = 0; iColumn < nCreate; iColumn++) {
            for (CoinBigIndex j = buildStart[iColumn]; j < buildStart[iColumn + 1]; j++) {
              int jRow = buildRow[j];
              jRow = lookupRow[jRow];
              assert(jRow >= 0 && jRow < nRow);
              buildRow[j] = jRow;
            }
          }
          finalModel->addCols(nCreate, buildStart, buildRow, buildElement, NULL, NULL, buildObj);
          int numberFinal = finalModel->getNumCols();
          for (iColumn = numberOther; iColumn < numberFinal; iColumn++) {
            if (markKnapsack[iKnapsack] < 0) {
              finalModel->setColUpper(iColumn, maxCoefficient);
              finalModel->setInteger(iColumn);
            } else {
              finalModel->setColUpper(iColumn, maxCoefficient + 1.0);
              finalModel->setInteger(iColumn);
            }
            OsiSimpleInteger *sosObject = new OsiSimpleInteger(finalModel, iColumn);
            sosObject->setPriority(1000000);
            object[nObj++] = sosObject;
            buildRow[iColumn - numberOther] = iColumn;
            buildElement[iColumn - numberOther] = 1.0;
          }
          if (markKnapsack[iKnapsack] < 0) {
            // convexity row
            finalModel->addRow(numberFinal - numberOther, buildRow, buildElement, 1.0, 1.0);
          } else {
            int iColumn = markKnapsack[iKnapsack];
            int n = numberFinal - numberOther;
            buildRow[n] = iColumn;
            buildElement[n++] = -fabs(coefficient[iKnapsack]);
            // convexity row (sort of)
            finalModel->addRow(n, buildRow, buildElement, 0.0, 0.0);
            OsiSOS *sosObject = new OsiSOS(finalModel, n - 1, buildRow, NULL, 1);
            sosObject->setPriority(iKnapsack + SOSPriority);
            // Say not integral even if is (switch off heuristics)
            sosObject->setIntegerValued(false);
            object[nSOS++] = sosObject;
          }
          numberOther = numberFinal;
        }
        finalModel->addObjects(nObj, object);
        for (iKnapsack = 0; iKnapsack < nObj; iKnapsack++)
          delete object[iKnapsack];
        delete[] object;
        // Can we move any rows to cuts
        const int *cutMarker = coinModel.cutMarker();
        if (cutMarker && 0) {
          printf("AMPL CUTS OFF until global cuts fixed\n");
          cutMarker = NULL;
        }
        if (cutMarker) {
          // Row copy
          const CoinPackedMatrix *matrixByRow = finalModel->getMatrixByRow();
          const double *elementByRow = matrixByRow->getElements();
          const int *column = matrixByRow->getIndices();
          const CoinBigIndex *rowStart = matrixByRow->getVectorStarts();
          const int *rowLength = matrixByRow->getVectorLengths();

          const double *rowLower = finalModel->getRowLower();
          const double *rowUpper = finalModel->getRowUpper();
          int nDelete = 0;
          for (iRow = 0; iRow < numberRows; iRow++) {
            if (cutMarker[iRow] && lookupRow[iRow] >= 0) {
              int jRow = lookupRow[iRow];
              whichRow[nDelete++] = jRow;
              CoinBigIndex start = rowStart[jRow];
              stored.addCut(rowLower[jRow], rowUpper[jRow],
                rowLength[jRow], column + start, elementByRow + start);
            }
          }
          finalModel->deleteRows(nDelete, whichRow);
        }
        knapsackStart[numberKnapsack] = finalModel->getNumCols();
        delete[] buildObj;
        delete[] buildElement;
        delete[] buildStart;
        delete[] buildRow;
        finalModel->writeMps("full");
      }
    }
  }
  delete[] whichKnapsack;
  delete[] markRow;
  delete[] markKnapsack;
  delete[] coefficient;
  delete[] linear;
  delete[] whichRow;
  delete[] lookupRow;
  delete si;
  si = NULL;
  if (!badModel && finalModel) {
    finalModel->setDblParam(OsiObjOffset, coinModel.objectiveOffset());
    return finalModel;
  } else {
    delete finalModel;
    printf("can't make knapsacks - did you set fixedPriority (extra1)\n");
    return NULL;
  }
}
#endif //COIN_HAS_LINK

// Fills in original solution (coinModel length)
void afterKnapsack(const CoinModel &coinModel2, const int *whichColumn, const int *knapsackStart,
  const int *knapsackRow, int numberKnapsack,
  const double *knapsackSolution, double *solution, int logLevel)
{
  CoinModel coinModel = coinModel2;
  int numberColumns = coinModel.numberColumns();
  int iColumn;
  // associate all columns to stop possible error messages
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    coinModel.associateElement(coinModel.columnName(iColumn), 1.0);
  }
  CoinZeroN(solution, numberColumns);
  int nCol = knapsackStart[0];
  for (iColumn = 0; iColumn < nCol; iColumn++) {
    int jColumn = whichColumn[iColumn];
    solution[jColumn] = knapsackSolution[iColumn];
  }
  int *buildRow = new int[numberColumns]; // wild overkill
  double *buildElement = new double[numberColumns];
  int iKnapsack;
  for (iKnapsack = 0; iKnapsack < numberKnapsack; iKnapsack++) {
    int k = -1;
    for (iColumn = knapsackStart[iKnapsack]; iColumn < knapsackStart[iKnapsack + 1]; iColumn++) {
      if (knapsackSolution[iColumn] > 1.0e-5) {
        if (k >= 0) {
          printf("Two nonzero values for knapsack %d at (%d,%g) and (%d,%g)\n", iKnapsack,
            k, knapsackSolution[k], iColumn, knapsackSolution[iColumn]);
          abort();
        }
        k = iColumn;
        assert(fabs(floor(knapsackSolution[iColumn] + 0.5) - knapsackSolution[iColumn]) < 1.0e-5);
      }
    }
    if (k >= 0) {
      int iRow = knapsackRow[iKnapsack];
      int nCreate = 10000;
      int nel = expandKnapsack(coinModel, iRow, nCreate, NULL, NULL, buildRow, buildElement, k - knapsackStart[iKnapsack]);
      assert(nel);
      if (logLevel > 0)
        printf("expanded column %d in knapsack %d has %d nonzero entries:\n",
          k - knapsackStart[iKnapsack], iKnapsack, nel);
      for (int i = 0; i < nel; i++) {
        int jColumn = buildRow[i];
        double value = buildElement[i];
        if (logLevel > 0)
          printf("%d - original %d has value %g\n", i, jColumn, value);
        solution[jColumn] = value;
      }
    }
  }
  delete[] buildRow;
  delete[] buildElement;
#if 0
   for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (solution[iColumn]>1.0e-5&&coinModel.isInteger(iColumn))
	 printf("%d %g\n",iColumn,solution[iColumn]);
   }
#endif
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
