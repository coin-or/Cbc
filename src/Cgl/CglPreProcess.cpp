// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <string>
#include <cassert>
#include <cmath> 
#include <vector>
#include <algorithm>
#include <cfloat>
#include <climits>
#include "CoinPragma.hpp" 
#include "CglPreProcess.hpp"
#include "CglMessage.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CglStored.hpp"
#include "CglCutGenerator.hpp"
#include "CoinTime.hpp"
#include "CoinSort.hpp"
#include "CoinDenseFactorization.hpp"
#include "CoinBuild.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinWarmStartBasis.hpp"
#ifdef CBC_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#else
#undef CBC_USE_PAPILO
#endif
#include "CglProbing.hpp"
#include "CglDuplicateRow.hpp"
#include "CoinPresolveDupcol.hpp"
#include "CglClique.hpp"
#include "CglBKClique.hpp"
//#define PRINT_DEBUG 1
//#define COIN_DEVELOP 1
#ifdef COIN_DEVELOP
static int whichMps = 0;
char nameMps[50];
#endif
#if CBC_USE_PAPILO
#include "papilo/core/Problem.hpp"
#include "papilo/core/Presolve.hpp"
static papilo::PresolveResult<double> presolveResult;
typedef struct {
  CglPreProcess * preProcess;  
  void * papiloProb; 
  int * mapping;
  ClpSimplex * presolvedModel;
  ClpSimplex * originalModel;
  OsiSolverInterface *beforePapiloEnd;
  int returnCode; // 0 unchanged, 1 reduced, -1 infeasible or unbounded
  int presolveType; // 1 beginning, 2 end ? 3 both //
  int numberThreads;
} papiloStruct;
static papiloStruct initialTry={0};
static papiloStruct papiloPresolve(ClpSimplex * inModel,
				   bool treatAsContinuous,
				   double timeLimit);
static ClpSimplex *postSolvedModel(papiloStruct papilo);
#endif 

OsiSolverInterface * 
CglPreProcess::preProcess(OsiSolverInterface &model,
  bool makeEquality, int numberPasses)
{
  // Tell solver we are in Branch and Cut
  model.setHintParam(OsiDoInBranchAndCut, true, OsiHintDo);
  // Default set of cut generators
  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(3);
  generator1.setMaxProbeRoot(model.getNumCols());
  generator1.setMaxElements(100);
  generator1.setMaxLookRoot(50);
  generator1.setRowCuts(3);
  // Add in generators
  addCutGenerator(&generator1);
  OsiSolverInterface *newSolver = preProcessNonDefault(model, makeEquality ? 1 : 0, numberPasses);
  // Tell solver we are not in Branch and Cut
  model.setHintParam(OsiDoInBranchAndCut, false, OsiHintDo);
  if (newSolver)
    newSolver->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo);
  return newSolver;
}
static void outSingletons(int &nCol, int &nRow,
  int *startCol, int *row, double *element,
  int *startRow, int *column)
{
  int iRow, iCol;
  bool singletons = false;
  int *countRow = new int[nRow];
  int *countCol = new int[nCol];
  int *temp = new int[nRow];
  // make row copy
  memset(countRow, 0, nRow * sizeof(int));
  memset(countCol, 0, nCol * sizeof(int));
  for (iCol = 0; iCol < nCol; iCol++) {
    for (int j = startCol[iCol]; j < startCol[iCol + 1]; j++) {
      int iRow = row[j];
      countRow[iRow]++;
      countCol[iCol]++;
    }
  }
  startRow[0] = 0;
  for (iRow = 0; iRow < nRow; iRow++) {
    int k = countRow[iRow] + startRow[iRow];
    temp[iRow] = startRow[iRow];
    startRow[iRow + 1] = k;
  }
  for (iCol = 0; iCol < nCol; iCol++) {
    for (int j = startCol[iCol]; j < startCol[iCol + 1]; j++) {
      int iRow = row[j];
      int k = temp[iRow];
      temp[iRow]++;
      column[k] = iCol;
    }
  }
  for (iRow = 0; iRow < nRow; iRow++) {
    if (countRow[iRow] <= 1)
      singletons = true;
  }
  for (iCol = 0; iCol < nCol; iCol++) {
    if (countCol[iCol] <= 1)
      singletons = true;
  }
  if (singletons) {
    while (singletons) {
      singletons = false;
      for (iCol = 0; iCol < nCol; iCol++) {
        if (countCol[iCol] == 1) {
          singletons = true;
          countCol[iCol] = 0;
          int iRow = row[startCol[iCol]];
          int start = startRow[iRow];
          int end = start + countRow[iRow];
          countRow[iRow]--;
          int j;
          for (j = start; j < end; j++) {
            if (column[j] == iCol) {
              column[j] = column[end - 1];
              break;
            }
          }
          assert(j < end);
        }
      }
      for (iRow = 0; iRow < nRow; iRow++) {
        if (countRow[iRow] == 1) {
          singletons = true;
          countRow[iRow] = 0;
          int iCol = column[startRow[iRow]];
          int start = startCol[iCol];
          int end = start + countCol[iCol];
          countCol[iCol]--;
          int j;
          for (j = start; j < end; j++) {
            if (row[j] == iRow) {
              row[j] = row[end - 1];
              if (element)
                element[j] = element[end - 1];
              break;
            }
          }
          assert(j < end);
        }
      }
    }
    // Pack down
    int newNrow = 0;
    for (iRow = 0; iRow < nRow; iRow++) {
      if (countRow[iRow] == 0) {
        temp[iRow] = -1;
      } else {
        assert(countRow[iRow] > 1);
        temp[iRow] = newNrow;
        newNrow++;
      }
    }
    int newNcol = 0;
    int nEl = 0;
    int iNext = 0;
    for (iCol = 0; iCol < nCol; iCol++) {
      int start = iNext;
      iNext = startCol[iCol + 1];
      if (countCol[iCol] == 0) {
        countCol[iCol] = -1;
      } else {
        assert(countCol[iCol] > 1);
        int end = start + countCol[iCol];
        countCol[iCol] = newNcol;
        int j;
        for (j = start; j < end; j++) {
          int iRow = row[j];
          iRow = temp[iRow];
          assert(iRow >= 0);
          row[nEl] = iRow;
          if (element)
            element[nEl] = element[j];
          nEl++;
        }
        newNcol++;
        startCol[newNcol] = nEl;
      }
    }
    newNrow = 0;
    nEl = 0;
    iNext = 0;
    for (iRow = 0; iRow < nRow; iRow++) {
      int start = iNext;
      iNext = startRow[iRow + 1];
      if (countRow[iRow] > 1) {
        int end = start + countRow[iRow];
        int j;
        for (j = start; j < end; j++) {
          int iCol = column[j];
          iCol = countCol[iCol];
          assert(iCol >= 0);
          column[nEl++] = iCol;
        }
        newNrow++;
        startRow[newNrow] = nEl;
      }
    }
    nRow = newNrow;
    nCol = newNcol;
  }
  delete[] countCol;
  delete[] countRow;
  delete[] temp;
}
static int makeIntegers2(OsiSolverInterface *model, int mode)
{
  // See whether we should make variables integer
  const double *objective = model->getObjCoefficients();
  const double *lower = model->getColLower();
  const double *upper = model->getColUpper();
  const double *rowLower = model->getRowLower();
  const double *rowUpper = model->getRowUpper();
  int numberRows = model->getNumRows();
  double *rhs = new double[numberRows];
  int *count = new int[numberRows];
  int iColumn;
  bool makeAll = (mode > 1);
  int numberColumns = model->getNumCols();
  // Column copy of matrix
  const double *element = model->getMatrixByCol()->getElements();
  const int *row = model->getMatrixByCol()->getIndices();
  const CoinBigIndex *columnStart = model->getMatrixByCol()->getVectorStarts();
  const int *columnLength = model->getMatrixByCol()->getVectorLengths();
  // Row copy
  CoinPackedMatrix matrixByRow(*model->getMatrixByRow());
  //const double * elementByRow = matrixByRow.getElements();
  const int *column = matrixByRow.getIndices();
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();
  int numberIntegers = 1;
  int totalNumberIntegers = 0;
  while (numberIntegers) {
    memset(rhs, 0, numberRows * sizeof(double));
    memset(count, 0, numberRows * sizeof(int));
    int currentNumber = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = start + columnLength[iColumn];
      if (upper[iColumn] == lower[iColumn]) {
        for (CoinBigIndex j = start; j < end; j++) {
          int iRow = row[j];
          rhs[iRow] += lower[iColumn] * element[j];
        }
      } else if (model->isInteger(iColumn)) {
        currentNumber++;
        for (CoinBigIndex j = start; j < end; j++) {
          int iRow = row[j];
          if (fabs(element[j] - floor(element[j] + 0.5)) > 1.0e-10)
            rhs[iRow] = COIN_DBL_MAX;
        }
      } else {
        for (CoinBigIndex j = start; j < end; j++) {
          int iRow = row[j];
          count[iRow]++;
          if (fabs(element[j]) != 1.0)
            rhs[iRow] = COIN_DBL_MAX;
        }
      }
    }
    //#define CBC_USEFUL_PRINTING 1
#if CBC_USEFUL_PRINTING > 1
    printf("Current number of integers is %d\n", currentNumber);
#endif
    // now look at continuous
    bool allGood = true;
    double direction = model->getObjSense();
    int numberObj = 0;
#if CBC_USEFUL_PRINTING > 1
    int numberEq = 0;
    int numberEqI = 0;
#endif
    int numberZero = 0;
    int numberNonZero = 0;
    if (false) {
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (upper[iColumn] > lower[iColumn]) {
          if (!model->isInteger(iColumn)) {
            CoinBigIndex start = columnStart[iColumn];
            CoinBigIndex end = start + columnLength[iColumn];
            int nC = 0;
            for (CoinBigIndex j = start; j < end; j++) {
              int iRow = row[j];
              if (count[iRow] > 1) {
                nC++;
              }
            }
            if (nC > 2) {
              for (CoinBigIndex j = start; j < end; j++) {
                int iRow = row[j];
                if (count[iRow] > 1)
                  count[iRow] = 999999;
              }
            }
          }
        }
      }
    }
    int *newInts = new int[numberColumns];
    // Columns to zap
    int nColumnZap = 0;
    int *columnZap = new int[numberColumns];
    char *noGoodColumn = new char[numberColumns];
    memset(noGoodColumn, 0, numberColumns);
    int nNew = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (upper[iColumn] > lower[iColumn]) {
        double objValue = objective[iColumn] * direction;
        bool thisGood = true;
        if ((objValue || makeAll) && !model->isInteger(iColumn)) {
          if (objValue) {
            numberObj++;
          } else if (columnLength[iColumn] == 1) {
            continue; // don't bother with singletons
          }
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = start + columnLength[iColumn];
          if (objValue >= 0.0) {
            // wants to be as low as possible
            if (lower[iColumn] < -1.0e10 || fabs(lower[iColumn] - floor(lower[iColumn] + 0.5)) > 1.0e-10) {
              thisGood = false;
            } else if (upper[iColumn] < 1.0e10 && fabs(upper[iColumn] - floor(upper[iColumn] + 0.5)) > 1.0e-10) {
              thisGood = false;
            }
            bool singletonRow = true;
            bool equality = false;
            int xxxx = 0;
#if CBC_USEFUL_PRINTING
            int nC = 0;
            bool badCount = false;
#endif
            for (CoinBigIndex j = start; j < end; j++) {
              int iRow = row[j];
              if (count[iRow] > 1) {
                singletonRow = false;
                //printf("col %d row%d element %g - row count %d\n",iColumn,iRow,element[j],count[iRow]);
#if CBC_USEFUL_PRINTING
                if (count[iRow] == 999999)
                  badCount = true;
#endif
                if (element[j] == 1.0) {
                  if ((xxxx & 1) == 0)
                    xxxx |= 1;
                  else
                    xxxx = 15;
                } else {
                  if ((xxxx & 2) == 0)
                    xxxx |= 2;
                  else
                    xxxx = 15;
                }
#if CBC_USEFUL_PRINTING
                nC++;
#endif
              } else if (rowLower[iRow] == rowUpper[iRow]) {
                equality = true;
              }
              double rhsValue = rhs[iRow];
              double lowerValue = rowLower[iRow];
              double upperValue = rowUpper[iRow];
              if (rhsValue < 1.0e20) {
                if (lowerValue > -1.0e20)
                  lowerValue -= rhsValue;
                if (upperValue < 1.0e20)
                  upperValue -= rhsValue;
              }
              if (fabs(rhsValue) > 1.0e20 || fabs(rhsValue - floor(rhsValue + 0.5)) > 1.0e-10
                || fabs(element[j]) != 1.0) {
                // no good
                thisGood = false;
                break;
              }
              if (element[j] > 0.0) {
                if (lowerValue > -1.0e20 && fabs(lowerValue - floor(lowerValue + 0.5)) > 1.0e-10) {
                  // no good
                  thisGood = false;
                  break;
                }
              } else {
                if (upperValue < 1.0e20 && fabs(upperValue - floor(upperValue + 0.5)) > 1.0e-10) {
                  // no good
                  thisGood = false;
                  break;
                }
              }
            }
#if CBC_USEFUL_PRINTING
            if (!model->isInteger(iColumn) && false)
              printf("%d has %d rows with >1 - state network %s interaction %s\n", iColumn, nC, xxxx > 3 ? "bad" : "good",
                badCount ? "too much" : "ok");
#endif
            // If not good here then mark rows
            if (!thisGood) {
              for (CoinBigIndex j = start; j < end; j++) {
                int iRow = row[j];
                count[iRow] = 999999;
              }
            }
            if (!singletonRow && end > start + 1 && !equality)
              thisGood = false;
            // Can we make equality
            if (end == start + 1 && !equality && false) {
#if CBC_USEFUL_PRINTING > 1
              numberEq++;
#endif
              int iRow = row[start];
              if (element[start] > 0.0)
                model->setRowUpper(iRow, rowLower[iRow]);
              else
                model->setRowLower(iRow, rowUpper[iRow]);
            }
          } else {
            // wants to be as high as possible
            if (upper[iColumn] > 1.0e10 || fabs(upper[iColumn] - floor(upper[iColumn] + 0.5)) > 1.0e-10) {
              thisGood = false;
            } else if (lower[iColumn] > -1.0e10 && fabs(lower[iColumn] - floor(lower[iColumn] + 0.5)) > 1.0e-10) {
              thisGood = false;
            }
            bool singletonRow = true;
            bool equality = false;
            for (CoinBigIndex j = start; j < end; j++) {
              int iRow = row[j];
              if (count[iRow] > 1) {
                singletonRow = false;
                thisGood = false;
              } else if (rowLower[iRow] == rowUpper[iRow]) {
                equality = true;
              }
              double rhsValue = rhs[iRow];
              double lowerValue = rowLower[iRow];
              double upperValue = rowUpper[iRow];
              if (rhsValue < 1.0e20) {
                if (lowerValue > -1.0e20)
                  lowerValue -= rhsValue;
                if (upperValue < 1.0e20)
                  upperValue -= rhsValue;
              }
              if (fabs(rhsValue) > 1.0e20 || fabs(rhsValue - floor(rhsValue + 0.5)) > 1.0e-10
                || fabs(element[j]) != 1.0) {
                // no good
                thisGood = false;
                break;
              }
              if (element[j] < 0.0) {
                if (lowerValue > -1.0e20 && fabs(lowerValue - floor(lowerValue + 0.5)) > 1.0e-10) {
                  // no good
                  thisGood = false;
                  break;
                }
              } else {
                if (upperValue < 1.0e20 && fabs(upperValue - floor(upperValue + 0.5)) > 1.0e-10) {
                  // no good
                  thisGood = false;
                  break;
                }
              }
            }
            if (!singletonRow && end > start + 1 && !equality)
              thisGood = false;
            // If not good here then mark rows
            if (!thisGood) {
              for (CoinBigIndex j = start; j < end; j++) {
                int iRow = row[j];
                count[iRow] = 999999;
              }
            }
            // Can we make equality
            if (end == start + 1 && !equality && false) {
#if CBC_USEFUL_PRINTING > 1
              numberEq++;
#endif
              int iRow = row[start];
              if (element[start] < 0.0)
                model->setRowUpper(iRow, rowLower[iRow]);
              else
                model->setRowLower(iRow, rowUpper[iRow]);
            }
          }
        } else if (objValue) {
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = start + columnLength[iColumn];
          if (end == start + 1) {
            int iRow = row[start];
            if (rowUpper[iRow] > rowLower[iRow] && !count[iRow]) {
              if (fabs(rhs[iRow]) > 1.0e20 || fabs(rhs[iRow] - floor(rhs[iRow] + 0.5)) > 1.0e-10
                || fabs(element[start]) != 1.0) {
                // no good
              } else if (false) {
#if CBC_USEFUL_PRINTING > 1
                numberEqI++;
#endif
                if (element[start] * objValue > 0.0)
                  model->setRowUpper(iRow, rowLower[iRow]);
                else
                  model->setRowLower(iRow, rowUpper[iRow]);
              }
            }
          }
        }
        if (!thisGood) {
          if (objValue)
            allGood = false;
          // look at again
          columnZap[nColumnZap++] = iColumn;
        } else if (makeAll && !model->isInteger(iColumn) && upper[iColumn] - lower[iColumn] < 10) {
          newInts[nNew++] = iColumn;
        }
      }
    }
    // Rows to look at
    int *rowLook = new int[numberRows];
    while (nColumnZap) {
      int nRowLook = 0;
      for (int i = 0; i < nColumnZap; i++) {
        int iColumn = columnZap[i];
        noGoodColumn[iColumn] = 1;
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = start + columnLength[iColumn];
        for (CoinBigIndex j = start; j < end; j++) {
          int iRow = row[j];
          if (count[iRow] != 999999) {
            count[iRow] = 999999;
            rowLook[nRowLook++] = iRow;
          }
        }
      }
      nColumnZap = 0;
      if (nRowLook) {
        for (int i = 0; i < nRowLook; i++) {
          int iRow = rowLook[i];
          CoinBigIndex start = rowStart[iRow];
          CoinBigIndex end = start + rowLength[iRow];
          for (CoinBigIndex j = start; j < end; j++) {
            int iColumn = column[j];
            if (upper[iColumn] > lower[iColumn] && !model->isInteger(iColumn)) {
              if (!noGoodColumn[iColumn]) {
                noGoodColumn[iColumn] = 1;
                columnZap[nColumnZap++] = iColumn;
              }
            }
          }
        }
      }
    }
    delete[] rowLook;
    delete[] noGoodColumn;
    delete[] columnZap;
    // Final look
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (upper[iColumn] > lower[iColumn] && !model->isInteger(iColumn) && objective[iColumn]) {
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = start + columnLength[iColumn];
        for (CoinBigIndex j = start; j < end; j++) {
          int iRow = row[j];
          if (count[iRow] == 999999)
            allGood = false;
        }
      }
    }
    // do if some but not too many
    if (nNew && nNew < currentNumber) {
      double tolerance;
      model->getDblParam(OsiPrimalTolerance, tolerance);
      for (int i = 0; i < nNew; i++) {
        int iColumn = newInts[i];
        double objValue = objective[iColumn];
        bool thisGood = true;
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = start + columnLength[iColumn];
        for (CoinBigIndex j = start; j < end; j++) {
          int iRow = row[j];
          if (count[iRow] == 999999) {
            thisGood = false;
            break;
          }
        }
        if (thisGood && upper[iColumn] < lower[iColumn] + 10.0) {
          model->setInteger(iColumn);
          // clean up bounds
          model->setColLower(iColumn, ceil(lower[iColumn] - tolerance));
          model->setColUpper(iColumn, floor(upper[iColumn] + tolerance));
          if (objValue)
            numberNonZero++;
          else
            numberZero++;
        } else if (objValue) {
          // unable to fix all with obj
          allGood = false;
        }
      }
    }
    delete[] newInts;
    // Can we look at remainder and make any integer
    if (makeAll && false) {
      int nLook = 0;
      int nEl = 0;
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (upper[iColumn] > lower[iColumn] && !model->isInteger(iColumn)) {
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = start + columnLength[iColumn];
          bool possible = true;
          int n = 0;
          for (CoinBigIndex j = start; j < end; j++) {
            int iRow = row[j];
            if (count[iRow] > 1) {
              if (count[iRow] == 999999) {
                possible = false;
                break;
              } else {
                n++;
              }
            }
          }
          if (possible) {
            nLook++;
            nEl += n;
          }
        }
      }
      if (nLook) {
        int *startC = new int[nLook + 1];
        int *back = new int[nLook];
        int *row2 = new int[nEl];
        double *element2 = new double[nEl];
        int *backRow = new int[numberRows];
        int jRow;
        for (jRow = 0; jRow < numberRows; jRow++) {
          backRow[jRow] = -1;
        }
        int nCol = nLook;
        nLook = 0;
        nEl = 0;
        startC[0] = 0;
        int nRow = 0;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (upper[iColumn] > lower[iColumn] && !model->isInteger(iColumn)) {
            CoinBigIndex start = columnStart[iColumn];
            CoinBigIndex end = start + columnLength[iColumn];
            bool possible = true;
            int n = 0;
            for (CoinBigIndex j = start; j < end; j++) {
              int iRow = row[j];
              if (count[iRow] > 1) {
                if (count[iRow] == 999999) {
                  possible = false;
                  break;
                } else {
                  n++;
                }
              }
            }
            if (!n)
              possible = false; // may be done later
            if (possible) {
              back[nLook] = iColumn;
              for (CoinBigIndex j = start; j < end; j++) {
                int iRow = row[j];
                if (count[iRow] > 1) {
                  int jRow = backRow[iRow];
                  if (jRow < 0) {
                    // new row
                    backRow[iRow] = nRow;
                    jRow = nRow;
                    nRow++;
                  }
                  element2[nEl] = element[j];
                  row2[nEl++] = jRow;
                }
              }
              nLook++;
              startC[nLook] = nEl;
            }
          }
        }
        // Redo nCol
        nCol = nLook;
        delete[] backRow;
        int *startRow = new int[nRow + 1];
        int *column2 = new int[nEl];
        // take out singletons and do row copy
        outSingletons(nCol, nRow,
          startC, row2, element2,
          startRow, column2);
        // Decompose
        int *rowBlock = new int[nRow];
        int *stack = new int[nRow];
        for (int iRow = 0; iRow < nRow; iRow++)
          rowBlock[iRow] = -2;
        int numberBlocks = 0;
        // to say if column looked at
        int *columnBlock = new int[nCol];
        int iColumn;
        for (iColumn = 0; iColumn < nCol; iColumn++)
          columnBlock[iColumn] = -2;
        for (iColumn = 0; iColumn < nCol; iColumn++) {
          int kstart = startC[iColumn];
          int kend = startC[iColumn + 1];
          if (columnBlock[iColumn] == -2) {
            // column not allocated
            int j;
            int nstack = 0;
            for (j = kstart; j < kend; j++) {
              int iRow = row2[j];
              if (rowBlock[iRow] != -1) {
                assert(rowBlock[iRow] == -2);
                rowBlock[iRow] = numberBlocks; // mark
                stack[nstack++] = iRow;
              }
            }
            if (nstack) {
              // new block - put all connected in
              numberBlocks++;
              columnBlock[iColumn] = numberBlocks - 1;
              while (nstack) {
                int iRow = stack[--nstack];
                int k;
                for (k = startRow[iRow]; k < startRow[iRow + 1]; k++) {
                  int iColumn = column2[k];
                  int kkstart = startC[iColumn];
                  int kkend = startC[iColumn + 1];
                  if (columnBlock[iColumn] == -2) {
                    columnBlock[iColumn] = numberBlocks - 1; // mark
                    // column not allocated
                    int jj;
                    for (jj = kkstart; jj < kkend; jj++) {
                      int jRow = row2[jj];
                      if (rowBlock[jRow] == -2) {
                        rowBlock[jRow] = numberBlocks - 1;
                        stack[nstack++] = jRow;
                      }
                    }
                  } else {
                    assert(columnBlock[iColumn] == numberBlocks - 1);
                  }
                }
              }
            } else {
              // Only in master
              columnBlock[iColumn] = -1;
              // empty - should already be integer
              abort();
            }
          }
        }
        // See if each block OK
        for (int iBlock = 0; iBlock < numberBlocks; iBlock++) {
          // Get block
          int *startCB = new int[nCol + 1];
          int *row2B = new int[nEl];
          int *startCC = new int[nCol + 1];
          int *row2C = new int[nEl];
          int *startRowC = new int[nRow + 1];
          int *column2C = new int[nEl];
          int *whichRow = new int[nRow];
          int *whichCol = new int[nCol];
          int i;
          int nRowB = 0;
          int nColB = 0;
          int nElB = 0;
          for (i = 0; i < nRow; i++) {
            if (rowBlock[i] == iBlock) {
              whichRow[i] = nRowB;
              nRowB++;
            } else {
              whichRow[i] = -1;
            }
          }
          bool network = true;
          // even if not network - take out network columns NO
          startCB[0] = 0;
          for (i = 0; i < nCol; i++) {
            if (columnBlock[i] == iBlock) {
              int type = 0;
              whichCol[i] = nColB;
              for (int j = startC[i]; j < startC[i + 1]; j++) {
                int iRow = row2[j];
                iRow = whichRow[iRow];
                if (iRow >= 0) {
                  if (element2[j] == 1.0) {
                    if ((type & 1) == 0)
                      type |= 1;
                    else
                      type = 7;
                  } else {
                    assert(element2[j] == -1.0);
                    if ((type & 2) == 0)
                      type |= 2;
                    else
                      type = 7;
                  }
                  row2B[nElB++] = iRow;
                }
              }
              if (type != 3)
                network = false;
              nColB++;
              startCB[nColB] = nElB;
              assert(startCB[nColB] > startCB[nColB - 1] + 1);
            } else {
              whichCol[i] = -1;
            }
          }
          // See if network
          bool goodInteger = false;
          if (!network) {
            // take out singletons
            outSingletons(nColB, nRowB,
              startCB, row2B, NULL,
              startRowC, column2C);
            // See if totally balanced;
            int *split = new int[nRowB];
            int *best = new int[nRowB];
            int *current = new int[nRowB];
            int *size = new int[nRowB];
            {
              memset(size, 0, nRowB * sizeof(int));
              for (i = 0; i < nColB; i++) {
                int j;
                for (j = startCB[i]; j < startCB[i + 1]; j++) {
                  int iRow = row2B[j];
                  size[iRow]++;
                }
              }
#if CBC_USEFUL_PRINTING
              for (i = 0; i < nRowB; i++)
                if (size[i] < 2)
                  printf("%d entries in row %d\n", size[i], i);
#endif
            }
            for (i = 0; i < nColB; i++)
              whichCol[i] = i;
            for (i = 0; i < nRowB; i++)
              whichRow[i] = 0;
            int nLeft = nColB;
            int nSet = 1;
            size[0] = nRowB;
            while (nLeft) {
              // find best column
              int iBest = -1;
              memset(best, 0, nSet * sizeof(int));
              memset(current, 0, nSet * sizeof(int));
              for (i = 0; i < nColB; i++) {
                if (whichCol[i] < nLeft) {
                  int j;
                  for (j = startCB[i]; j < startCB[i + 1]; j++) {
                    int iRow = row2B[j];
                    int iSet = whichRow[iRow];
                    current[iSet]++;
                  }
                  // See if better - could this be done faster
                  bool better = false;
                  for (j = nSet - 1; j >= 0; j--) {
                    if (current[j] > best[j]) {
                      better = true;
                      break;
                    } else if (current[j] < best[j]) {
                      break;
                    }
                  }
                  if (better) {
                    iBest = i;
                    memcpy(best, current, nSet * sizeof(int));
                  }
                  for (j = startCB[i]; j < startCB[i + 1]; j++) {
                    int iRow = row2B[j];
                    int iSet = whichRow[iRow];
                    current[iSet] = 0;
                  }
                }
              }
              assert(iBest >= 0);
              // swap
              for (i = 0; i < nColB; i++) {
                if (whichCol[i] == nLeft - 1) {
                  whichCol[i] = whichCol[iBest];
                  whichCol[iBest] = nLeft - 1;
                  break;
                }
              }
              // See which ones will have to split
              int nMore = 0;
              for (i = 0; i < nSet; i++) {
                current[i] = i + nMore;
                if (best[i] > 0 && best[i] < size[i]) {
                  split[i] = i + nMore;
                  nMore++;
                } else {
                  split[i] = -1;
                }
              }
              if (nMore) {
                int j;
                for (j = startCB[iBest]; j < startCB[iBest + 1]; j++) {
                  int iRow = row2B[j];
                  int iSet = whichRow[iRow];
                  int newSet = split[iSet];
                  if (newSet >= 0) {
                    whichRow[iRow] = newSet + 1 + nRowB;
                  }
                }
                nSet += nMore;
                memset(size, 0, nSet * sizeof(int));
                for (i = 0; i < nRowB; i++) {
                  int iSet = whichRow[i];
                  if (iSet >= nRowB) {
                    // has 1 - correct it
                    iSet -= nRowB;
                  } else {
                    // 0 part of split set or not split
                    iSet = current[iSet];
                  }
                  whichRow[i] = iSet;
                  size[iSet]++;
                }
              }
              nLeft--;
            }
            if (nSet < nRowB) {
              // ties - need to spread out whichRow
              memset(split, 0, nRowB * sizeof(int));
              for (i = 0; i < nRowB; i++) {
                int iSet = whichRow[i];
                split[iSet]++;
              }
              current[0] = 0;
              for (i = 0; i < nSet; i++) {
                current[i + 1] = current[i] + split[i];
                split[i] = current[i];
              }
              for (i = 0; i < nRowB; i++) {
                int iSet = whichRow[i];
                int k = split[iSet];
                split[iSet] = k;
                whichRow[i] = k;
              }
            }
            // Get inverse of whichCol
            for (i = 0; i < nColB; i++) {
              int iColumn = whichCol[i];
              startCC[iColumn] = i;
            }
            memcpy(whichCol, startCC, nColB * sizeof(int));
            // Permute matrix
            startCC[0] = 0;
            int nelB = 0;
            memset(split, 0, nRowB * sizeof(int));
            for (i = 0; i < nColB; i++) {
              int iColumn = whichCol[i];
              int j;
              for (j = startCB[iColumn]; j < startCB[iColumn + 1]; j++) {
                int iRow = row2B[j];
                int iSet = whichRow[iRow];
                row2C[nelB++] = iSet;
                split[iSet]++;
              }
              startCC[i + 1] = nelB;
            }
            startRowC[0] = 0;
            for (i = 0; i < nRowB; i++) {
              startRowC[i + 1] = startRowC[i] + split[i];
              split[i] = 0;
            }
            for (i = 0; i < nColB; i++) {
              int j;
              for (j = startCC[i]; j < startCC[i + 1]; j++) {
                int iRow = row2C[j];
                int k = split[iRow] + startRowC[iRow];
                split[iRow]++;
                column2C[k] = i;
              }
            }
            for (i = 0; i < nRowB; i++)
              split[i] = 0;
            goodInteger = true;
            for (i = nColB - 1; i > 0; i--) {
              int j;
              for (j = startCC[i]; j < startCC[i + 1]; j++) {
                int iRow = row2C[j];
                split[iRow] = 1;
              }
              for (j = startCC[i]; j < startCC[i + 1]; j++) {
                int iRow = row2C[j];
                for (int k = startRowC[iRow]; k < startRowC[iRow + 1]; k++) {
                  int iColumn = column2C[k];
                  if (iColumn < i) {
                    for (int jj = startCC[iColumn]; jj < startCC[iColumn + 1]; jj++) {
                      int jRow = row2C[jj];
                      if (jRow > iRow && !split[jRow]) {
                        // bad
                        goodInteger = false;
                        break;
                      }
                    }
                  }
                }
              }
              if (!goodInteger)
                break;
              for (j = startCC[i]; j < startCC[i + 1]; j++) {
                int iRow = row2C[j];
                split[iRow] = 0;
              }
            }
            delete[] split;
            delete[] best;
            delete[] current;
            delete[] size;
          } else {
            // was network
            goodInteger = true;
          }
          if (goodInteger) {
#if CBC_USEFUL_PRINTING
            printf("Block %d can be integer\n", iBlock);
#endif
            for (i = 0; i < nCol; i++) {
              if (columnBlock[i] == iBlock) {
                int iBack = back[i];
                model->setInteger(iBack);
              }
            }
          }
          delete[] startRowC;
          delete[] column2C;
          delete[] startCB;
          delete[] row2B;
          delete[] startCC;
          delete[] row2C;
          delete[] whichRow;
          delete[] whichCol;
        }
        delete[] startRow;
        delete[] column2;
        delete[] element2;
        delete[] startC;
        delete[] row2;
        delete[] back;
      }
    }
    numberIntegers = numberNonZero;
    if (allGood && numberObj) {
#if CBC_USEFUL_PRINTING > 1
      int numberPossible = 0;
#endif
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (upper[iColumn] > lower[iColumn] && objective[iColumn] && !model->isInteger(iColumn)) {
#if CBC_USEFUL_PRINTING > 1
          numberPossible++;
#endif
          if (upper[iColumn] <= lower[iColumn] + 10) {
            model->setInteger(iColumn);
            numberIntegers++;
          }
        }
      }
#if CBC_USEFUL_PRINTING > 1
      printf("ZZZZYY CglPreProcess analysis says all (%d) continuous with costs could be made integer - %d were\n", numberPossible, numberIntegers - numberNonZero);
#endif
    }
#if CBC_USEFUL_PRINTING > 1
    if (numberZero)
      printf("ZZZZYY %d continuous with zero cost were made integer\n", numberZero);
#endif
    numberIntegers += numberZero;
#if CBC_USEFUL_PRINTING > 1
    if (numberEq || numberEqI)
      printf("ZZZZYY %d rows made equality from continuous, %d from integer\n", numberEq, numberEqI);
#endif
    totalNumberIntegers += numberIntegers;
    if (!makeAll)
      numberIntegers = 0;
  }
  delete[] rhs;
  delete[] count;
  return (totalNumberIntegers);
}
//#define DEBUG_PREPROCESS 1
#if DEBUG_PREPROCESS
#if DEBUG_PREPROCESS > 1
extern double *debugSolution;
extern int debugNumberColumns;
const OsiRowCutDebugger * debugger = NULL;
#endif
static int mpsNumber = 0;
static void writeDebugMps(const OsiSolverInterface *solver,
  const char *where,
  OsiPresolve *pinfo)
{
  mpsNumber++; 
  char name[30];  
  sprintf(name, "/tmp/presolve%2.2d.mps", mpsNumber);
  printf("saving %s from %s - %d row, %d columns\n", 
    name, where, solver->getNumRows(), solver->getNumCols());
  if (mpsNumber>20) {
    printf("Not saving\n");
  } else {
    solver->writeMpsNative(name, NULL, NULL, 0, 1, 0);
    sprintf(name, "/tmp/presolve%2.2d.bas", mpsNumber);
    solver->writeBasisNative(name);
  }
#if DEBUG_PREPROCESS > 1
  if (pinfo && debugSolution) {
    int n = solver->getNumCols();
    if (n < debugNumberColumns) {
      const int *original = pinfo->originalColumns();
      if (!original) {
        printf("No original columns\n");
        abort();
      }
      for (int i = 0; i < n; i++)
        debugSolution[i] = debugSolution[original[i]];
      debugNumberColumns = n;
    }
  }
  if (debugSolution) {
    OsiSolverInterface *newSolver = solver->clone();
    const double *lower = newSolver->getColLower();
    const double *upper = newSolver->getColUpper();
    for (int i = 0; i < debugNumberColumns; i++) {
      if (newSolver->isInteger(i)) {
        double value = floor(debugSolution[i] + 0.5);
        if (value < lower[i] || value > upper[i]) {
          printf("Bad value %d - %g %g %g\n", i, lower[i], debugSolution[i],
            upper[i]);
        } else {
          newSolver->setColLower(i, value);
          newSolver->setColUpper(i, value);
        }
      }
    }
    printf("Starting solve %d\n", mpsNumber);
    newSolver->resolve();
    printf("Ending solve %d - status %s obj %g\n", mpsNumber,
      newSolver->isProvenOptimal() ? "ok" : "bad",
      newSolver->getObjValue());
    if (!newSolver->isProvenOptimal()) {
      newSolver->writeMpsNative("badmodel.mps", NULL, NULL, 0, 1, 0);
      exit(77);
    } else {
      // update
      int n = newSolver->getNumCols();
      if (n>debugNumberColumns) {
	double * temp = CoinCopyOfArray(newSolver->getColSolution(),n);
	delete [] debugSolution;
	debugSolution = temp;
	debugNumberColumns = n;
      }
    }
    delete newSolver;
  }
#endif
}
#else
#define writeDebugMps(x, y, z)
#endif
#define USE_CGL_RATIONAL 1
#if USE_CGL_RATIONAL>0
#include "CoinRational.hpp"
static int64_t computeGcd(int64_t a, int64_t b) {
  // This is the standard Euclidean algorithm for gcd
  int64_t remainder = 1;
  // Make sure a<=b (will always remain so)
  if (a > b) {
    // Swap a and b
    int64_t temp = a;
    a = b;
    b = temp;
  }
  // If zero then gcd is nonzero
  if (!a) {
    if (b) {
      return b;
    } 
    else {
      printf("### WARNING: CglGMI::computeGcd() given two zeroes!\n");
      exit(1);
    }
  }
  while (remainder) {
    remainder = b % a;
    b = a;
    a = remainder;
  }
  return b;
} /* computeGcd */
static bool scaleRowIntegral(double* rowElem, int rowNz)
{
  int64_t gcd, lcm;
  double maxdelta = 1.0e-13;
  double maxscale = 1000; 
  int64_t maxdnom = 1000;
  //int64_t numerator = 0, denominator = 0;
  // Initialize gcd and lcm
  CoinRational r = CoinRational(rowElem[0], maxdelta, maxdnom);
  if (r.getNumerator() != 0){
    gcd = llabs(r.getNumerator());
    lcm = r.getDenominator();
  } else {
    return false;
  } 
  for (int i = 1; i < rowNz; ++i) {
    if (rowElem[i]) {
      r = CoinRational(rowElem[i], maxdelta, maxdnom);
      if (r.getNumerator() != 0){
	gcd = computeGcd(gcd, r.getNumerator());
	lcm *= r.getDenominator()/(computeGcd(lcm,r.getDenominator()));
      } else {
	return false;
      }
    }
  }
  double scale = ((double)lcm)/((double)gcd);
  if (fabs(scale) > maxscale) {
    return false;
  }
  scale = fabs(scale);
  // Looks like we have a good scaling factor; scale and return;
  for (int i = 0; i < rowNz; ++i) {
    double value = rowElem[i]*scale;
    rowElem[i] = floor(value+0.5);
    if (fabs(rowElem[i]-value)>1.0e-10) 
      return false;
  }
  return true;
} /* scaleRowIntegral */
#endif
// returns -1 if infeasible, +n made integer
static int analyze(OsiSolverInterface * solver)
{
  const double *objective = solver->getObjCoefficients() ;
  int numberColumns = solver->getNumCols() ;
  char * intVar = new char[numberColumns];
  for (int i=0;i<numberColumns;i++) {
    if (solver->isInteger(i))
      intVar[i] = 1;
    else
      intVar[i] = 0;
  }
  double * lower = CoinCopyOfArray(solver->getColLower(),numberColumns);
  double * upper = CoinCopyOfArray(solver->getColUpper(),numberColumns);
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
  int * which = new int[numberRows];
  double * changeRhs = new double[numberRows];
  memset(changeRhs,0,numberRows*sizeof(double));
  memset(ignore,0,numberRows);
  int numberChanged=0;
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
          if (!intVar[jColumn]) {
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
              numberChanged++;
              intVar[jColumn1]=77;
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
             (take out 2 for now - until fixed)
          */
          if (fabs(value1)==1.0&&value1*value2==-1.0&&!lower[jColumn1]
              &&!lower[jColumn2]&&columnLength[jColumn1]==1&&columnLength[jColumn2]==1) {
            int n=0;
            CoinBigIndex i;
            double objChange=direction*(objective[jColumn1]+objective[jColumn2]);
            double bound = std::min(upper[jColumn1],upper[jColumn2]);
            bound = std::min(bound,1.0e20);
            for ( i=columnStart[jColumn1];i<columnStart[jColumn1]+columnLength[jColumn1];i++) {
              int jRow = row[i];
              double value = element[i];
              if (jRow!=iRow) {
                which[n++]=jRow;
                changeRhs[jRow]=value;
              }
            }
            for ( i=columnStart[jColumn2];i<columnStart[jColumn2]+columnLength[jColumn2];i++) {
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
                numberChanged++;
                intVar[jColumn1]=77;
                numberChanged++;
                intVar[jColumn2]=77;
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
      if (upper[iColumn] > lower[iColumn]+1.0e-8&&!intVar[iColumn]) {
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
          numberChanged++;
          intVar[iColumn]=77;
        }
      }
    }
    finished = numberChanged==saveNumberChanged;
  }
  bool feasible=true;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (intVar[iColumn]==77) {
      if (upper[iColumn]>1.0e20) {
        upper[iColumn] = 1.0e20;
      } else {
        upper[iColumn] = floor(upper[iColumn]+1.0e-5);
      }
      if (lower[iColumn]<-1.0e20) {
        lower[iColumn] = -1.0e20;
      } else {
        lower[iColumn] = ceil(lower[iColumn]-1.0e-5);
        if (lower[iColumn]>upper[iColumn])
          feasible=false;
      }
      if (lower[iColumn]==0.0&&upper[iColumn]==1.0)
        intVar[iColumn]=1;
      else if (lower[iColumn]==upper[iColumn])
        intVar[iColumn]=0;
      else
        intVar[iColumn]=2;
    }
  }
  delete [] which;
  delete [] changeRhs;
  delete [] ignore;
  if (numberChanged) {
    //printf("%d variables made integer\n",numberChanged);
    for (int i=0;i<numberColumns;i++) {
      if (intVar[i])
	solver->setInteger(i);
    }
  }
  delete [] intVar;
  delete [] lower;
  delete [] upper;
  if (feasible)
    return numberChanged;
  else
    return -1;
}
OsiSolverInterface *
CglPreProcess::preProcessNonDefault(OsiSolverInterface &model,
  int makeEquality, int numberPasses,
  int tuning)
{
  double ppstart = getCurrentCPUTime();
  OsiSolverInterface * modelIn = &model;
#if CBC_USE_PAPILO
  papiloStruct keepPapilo = initialTry;
  OsiClpSolverInterface * clpSolverIn
    = getClpSolverInterface(&model);
  if ((initialTry.presolveType&3)!=0)
    makeEquality = -2;
  bool doPapiloPresolve = (initialTry.presolveType&1)!=0;
  if (clpSolverIn && doPapiloPresolve) {
    initialTry = papiloPresolve(clpSolverIn->getModelPtr(),false,
				timeLimit_-ppstart);
    initialTry.presolveType &= ~1;
    char generalPrint[100];
    sprintf(generalPrint,
	    "papilo presolve took %g seconds",getCurrentCPUTime()-ppstart);
    handler_->message(CGL_GENERAL, messages_)
      << generalPrint
      << CoinMessageEol;
    if (initialTry.returnCode<0) {
      // infeasible
      sprintf(generalPrint,"papilo says infeasiblen");
      handler_->message(CGL_GENERAL, messages_)
	<< generalPrint
	<< CoinMessageEol;
      return NULL;
    } else if (initialTry.returnCode==0) {
      // no change - look out for memory leaks
      doPapiloPresolve = false;
      memset(&initialTry,0,sizeof(papiloStruct));
    } else {
      initialTry.preProcess = this;
      initialTry.presolveType |= 64; // say type 1
      // swap
      ClpSimplex * presolved = initialTry.presolvedModel;
      int numberColumns2 = presolved->numberColumns();
      assert (numberColumns2<=clpSolverIn->getNumCols());
      for (int i=0;i<numberColumns2;i++) {
	if (presolved->isInteger(i))
	  clpSolverIn->setIntegerType(i,1);
	else
	  clpSolverIn->setIntegerType(i,0);
      }      
      initialTry.originalModel = clpSolverIn->swapModelPtr(initialTry.presolvedModel);
      if (makeEquality==-2)
	makeEquality = 1; // don't do below
    }
  }
#endif
#define CGL_TRY_MINI_DUAL_STUFF
#ifdef CGL_TRY_MINI_DUAL_STUFF
  if (makeEquality==-2) {
    OsiPresolve dummy;
    // Just to do dual stuff using existing coding
    //printf("start mini\n");
    OsiSolverInterface * returnedModel = dummy.miniPresolvedModel(*modelIn);
    delete returnedModel; // throw it away
    //printf("end mini\n");
  }
#endif
#if DEBUG_PREPROCESS > 1
  bool rcdActive = true;
  std::string modelName;
  modelIn->getStrParam(OsiProbName, modelName);
  writeDebugMps(&model, "IPP:preProcessNonDefault", 0);
  std::cout
    << "  Attempting to activate row cut debugger for "
    << modelName << " ... ";
  if (!appData_) {
    modelIn->activateRowCutDebugger(modelName.c_str());
  } else {
    // see if passed in
    double *solution = CoinCopyOfArray(reinterpret_cast< double * >(appData_),
      modelIn->getNumCols());
    modelIn->activateRowCutDebugger(solution);
    delete[] solution;
  }
  debugger = modelIn->getRowCutDebugger();
  if (debugger)
    std::cout << "on optimal path." << std::endl;
  else if (modelIn->getRowCutDebuggerAlways())
    std::cout << "not on optimal path." << std::endl;
  else {
    std::cout << "failure." << std::endl;
    rcdActive = false;
  }
  if (rcdActive) {
    const OsiRowCutDebugger *debugger = modelIn->getRowCutDebuggerAlways();
    std::cout << "  Optimal solution is:" << std::endl;
    debugger->printOptimalSolution(*modelIn);
  }
  debugger = NULL;
#else
  const OsiRowCutDebugger *debugger = NULL;
#endif
  originalModel_ = modelIn;
  // See if SC variables
  double * scBound = NULL;
  if ((options_&128)!=0) {
    assert (appData_);
    typedef struct {
      double low;
      double high;
      int column;
    } lotStruct;
    typedef struct {lotStruct * lotsize;int numberLotSizing;} templot;
    templot * temp = reinterpret_cast<templot *>(appData_);
    lotStruct * lotsize=temp->lotsize;
    int numberLotSizing=temp->numberLotSizing;
    int numberColumns = modelIn->getNumCols();
    scBound = new double [numberColumns];
    const double * lower = modelIn->getColLower();
    const double * upper = modelIn->getColUpper();
    for (int i=0;i<numberColumns;i++)
      scBound[i]=-COIN_DBL_MAX;
    for (int i=0;i<numberLotSizing;i++) {
      int iColumn = lotsize[i].column;
      scBound[iColumn] = lotsize[i].low;
      assert (lower[iColumn]==0.0);
      assert (upper[iColumn]==lotsize[i].high);
    }
  }
  if (tuning >= 1000000) {
    numberPasses = tuning / 1000000;
    tuning %= 1000000;
    //minimumLength = tuning;
  }
  if (numberPasses != 99) {
    numberSolvers_ = numberPasses;
  } else {
    numberSolvers_ = 1;
  }
  model_ = new OsiSolverInterface *[numberSolvers_];
  modifiedModel_ = new OsiSolverInterface *[numberSolvers_];
  presolve_ = new OsiPresolve *[numberSolvers_];
  for (int i = 0; i < numberSolvers_; i++) {
    model_[i] = NULL;
    modifiedModel_[i] = NULL;
    presolve_[i] = NULL;
  }
  // clear original
  delete[] originalColumn_;
  delete[] originalRow_;
  originalColumn_ = NULL;
  originalRow_ = NULL;
  // Should not be hardwired tolerance
#ifndef CGL_PREPROCESS_TOLERANCE
  // So standalone version can switch off
  double feasibilityTolerance;
  if (((tuning%10000) & 1024) == 0)
    modelIn->getDblParam(OsiPrimalTolerance, feasibilityTolerance);
   else
     feasibilityTolerance = 1.0e-4;
#else
  // So standalone version can switch off
  double feasibilityTolerance = ((tuning & 1024) == 0)
    ? CGL_PREPROCESS_TOLERANCE : 1.0e-4;
#endif
  if (numberPasses==99) {
    // keep very simple
    OsiSolverInterface *presolvedModel;
    OsiPresolve *pinfo = new OsiPresolve();
    int presolveActions = 0;
    // Allow dual stuff on integers
    // Allow stuff which may not unroll cleanly - unless told not to
    if ((tuning & 4096) == 0)
      presolveActions = 1 + 16;
    else
      presolveActions = 16; // actually just switch off duplicate columns for ints
    if ((tuning & 32) != 0)
      presolveActions |= 32;
    presolveActions |= 8;
    pinfo->setPresolveActions(presolveActions);
    if (prohibited_)
      assert(numberProhibited_ == originalModel_->getNumCols());
    presolvedModel =
      pinfo->presolvedModel(*originalModel_, feasibilityTolerance, true,
			    5, prohibited_, true, rowType_,scBound);
    startModel_ = originalModel_;
    if (presolvedModel) {
      // update prohibited and rowType
      update(pinfo, presolvedModel,scBound);
      model_[0] = presolvedModel;
      presolve_[0] = pinfo;
      modifiedModel_[0] = presolvedModel->clone();
      createOriginalIndices();
      numberSolvers_ = 99; // mark as odd
      return modifiedModel_[0];
    } else {
      numberSolvers_ = 1;
      return NULL;
    }
  }
  //startModel_=&model;
  // make clone
  delete startModel_;
  startModel_ = originalModel_->clone();
  int numberRows = originalModel_->getNumRows();
  int numberColumns = originalModel_->getNumCols();
  int *rows = new int[2*numberRows];
  double *element = new double[numberRows];
  // probably okay with rowType_ but for now ..
  CoinPackedMatrix matrixByRow(*originalModel_->getMutableMatrixByRow());
  if (!rowType_ && true) {
    writeDebugMps(originalModel_, "before_scaled", 0);
    // clean matrix
    double *elementByRow = matrixByRow.getMutableElements();
    const int *column = matrixByRow.getIndices();
    const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
    const int *rowLength = matrixByRow.getVectorLengths();
    
    const double *rowLower = originalModel_->getRowLower();
    const double *rowUpper = originalModel_->getRowUpper();
    double * temp = new double [2*numberColumns+4];
    double * tempSave = temp+numberColumns+2;
    int nChanged = 0;
    int nScaled = 0;
    for (int iRow=0;iRow<numberRows;iRow++) {
      int n = 0;
      // make majority positive (unless upper rhs +1)
      if (!rowLength[iRow])
	continue;
      CoinBigIndex start = rowStart[iRow];
      CoinBigIndex end = start+rowLength[iRow];
      double multiplier = 1.0;
      for (CoinBigIndex j=start;j<end;j++) {
	if (elementByRow[j]<0)
	  n++;
      }
      if (2*n>end-start && rowUpper[iRow] != 1.0)
	multiplier = -1.0;
      int nInRow = end-start;
      n = nInRow;
      double smallest = COIN_DBL_MAX;
      if (multiplier == 1.0) {
	for (int i=0;i<n;i++) {
	  double value = elementByRow[i+start];
	  if (fabs(value) < smallest)
	    smallest = fabs(value);
	  temp[i] = value;
	}
      } else {
	for (int i=0;i<n;i++) {
	  double value = elementByRow[i+start];
	  if (fabs(value) < smallest)
	    smallest = fabs(value);
	  temp[i] = -value;
	}
      }
      double lower = rowLower[iRow];
      double upper = rowUpper[iRow];
      if (smallest != 1.0) {
	if (lower>-1.0e20) 
	  temp[n++] = multiplier*lower;
	if (upper<1.0e20) 
	  temp[n++] = multiplier*upper;
	memcpy(tempSave,temp,n*sizeof(double));
	if (scaleRowIntegral(temp, n)) {
	  // double check
	  double largestError = 0.0;
	  double mult = temp[0]/elementByRow[start];
	  if (fabs(mult-floor(mult+0.1))<1.0e-12)
	    mult = floor(mult+0.1);
	  for (int i=0;i<n;i++) {
	    double value = mult*tempSave[i];
	    if (value) {
	      double vint = floor(value+0.01);
	      largestError = std::max(largestError,fabs(value-vint));
	    }
	  }
	  if (largestError<1.0e-9) {
	    multiplier = mult;
	  }
	}
      }
      element[iRow] = multiplier;
      if (multiplier!=1.0) {
	nChanged++;
	if (fabs(multiplier)!=1.0)
	  nScaled++;
	memcpy(elementByRow+start,temp,(end-start)*sizeof(double));
	if (multiplier<0.0) {
	  double tempV = lower;
	  lower = -upper;
	  upper = -tempV;
	}
	if (lower>-1.0e20)
	  lower *= fabs(multiplier);
	if (upper<1.0e20) 
	  upper *= fabs(multiplier);
	originalModel_->setRowLower(iRow,lower);
	originalModel_->setRowUpper(iRow,upper);
      }
    }
    if (nChanged) {
      CoinPackedMatrix * matrixByColumn =
	originalModel_->getMutableMatrixByCol();
      const int *row = matrixByColumn->getIndices();
      const CoinBigIndex *columnStart = matrixByColumn->getVectorStarts();
      const int *columnLength = matrixByColumn->getVectorLengths();
      double *columnElements = matrixByColumn->getMutableElements();
      double largestDelta = 0.0;
      for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	for (CoinBigIndex j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  double value = columnElements[j];
	  if (fabs(element[iRow])!=1.0) {
	    value *= element[iRow];
	    double vint = floor(value+0.01);
	    largestDelta = std::max(largestDelta,fabs(value-vint));
	    assert (largestDelta<1.0e-9);
	    columnElements[j] = vint;
	    assert (fabs(vint)>0.9);
	  } else if (element[iRow]==-1.0) {
	    columnElements[j] = -value;
	  }
	}
      }
#if CBC_USEFUL_PRINTING > 1
      printf("%d rows changed %d scaled - largest error %g\n",
	     nChanged,nScaled,largestDelta);
#endif
      writeDebugMps(originalModel_, "scaled", 0);
      delete startModel_;
      startModel_ = originalModel_->clone();
    }
    delete [] temp;
  }
  if (rowType_)
    assert(numberRowType_ == numberRows);
  //int originalNumberColumns=numberColumns;
  int minimumLength = 5;
  int numberModifiedPasses = 10;
  if (numberPasses <= 1)
    numberModifiedPasses = 1; // lightweight preprocessing
  else if (numberPasses <= 2)
    numberModifiedPasses = 2; // fairly lightweight preprocessing
  if (tuning >= 10000) {
    numberModifiedPasses = tuning / 10000;
    tuning %= 10000;
    //minimumLength = tuning;
  }
  if ((tuning & 1) != 0)
    options_ |= (16|512); // heavy stuff
  //bool heavyProbing = (tuning&1)!=0;
  int makeIntegers = (tuning & 6) >> 1;
  // See if we want to do initial presolve
  int doInitialPresolve = 1;
  if (numberSolvers_ < 2)
    doInitialPresolve = 0;
  // We want to add columns
  int numberSlacks = 0;

  int iRow;

  int numberCliques = 0;
  int *which = new int[numberColumns];

  // Statistics
  int totalP1 = 0, totalM1 = 0;
  int numberFixed = 0;
  // May just find it is infeasible
  bool feasible = true;

  // Row copy
  const double *elementByRow = matrixByRow.getElements();
  const int *column = matrixByRow.getIndices();
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();

  const double *lower = originalModel_->getColLower();
  const double *upper = originalModel_->getColUpper();
  const double *rowLower = originalModel_->getRowLower();
  const double *rowUpper = originalModel_->getRowUpper();
  // Clean bounds
  int iColumn;
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (originalModel_->isInteger(iColumn)) {
      double lo = std::max(lower[iColumn], ceil(lower[iColumn] - 1.0e-6));
      if (lo > lower[iColumn])
        originalModel_->setColLower(iColumn, lo);
      double up = std::min(upper[iColumn], floor(upper[iColumn] + 1.0e-6));
      if (up < upper[iColumn])
        originalModel_->setColUpper(iColumn, up);
      if (lo > up)
        feasible = false;
    }
  }
  bool allToGub = makeEquality == 5;
  if (allToGub)
    makeEquality = 3;
  // Initialize random seed
  CoinThreadRandom randomGenerator(987654321);
  bool justOnesWithObj = false;
  if (makeEquality == 2 || makeEquality == 3 || makeEquality == 4) {
    int iRow, iColumn;
    int numberIntegers = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (originalModel_->isInteger(iColumn))
        numberIntegers++;
    }
    // Look for possible SOS
    int numberSOS = 0;
    int *mark = new int[numberColumns];
    CoinFillN(mark, numberColumns, -1);
    int numberOverlap = 0;
    int numberInSOS = 0;
    // See if worthwhile creating accumulation variables
    int firstOther = numberRows;
    int *whichRow = new int[numberRows];
    for (iRow = 0; iRow < numberRows; iRow++) {
      if (rowUpper[iRow] == 1.0) {
        if (rowLength[iRow] < 5)
          continue;
        bool goodRow = true;
        bool overlap = false;
        for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
          int iColumn = column[j];
          if (elementByRow[j] != 1.0 || !originalModel_->isInteger(iColumn) || lower[iColumn]) {
            goodRow = false;
            break;
          }
          if (mark[iColumn] >= 0) {
            overlap = true;
            numberOverlap++;
          }
        }
        if (goodRow) {
          if (!overlap) {
            // mark all
            for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
              int iColumn = column[j];
              mark[iColumn] = numberSOS;
            }
            numberSOS++;
            numberInSOS += rowLength[iRow];
          }
          // May still be interesting even if overlap
          if (rowLength[iRow] >= 5) {
            firstOther--;
            whichRow[firstOther] = iRow;
          }
        }
      }
    }
    if (makeEquality == 2 && false) {
      if (numberOverlap || numberIntegers > numberInSOS + 1) {
        // try just ones with costs
        CoinFillN(mark, numberColumns, -1);
        numberOverlap = 0;
        numberInSOS = 0;
        bool allCostsInSOS = true;
        const double *objective = originalModel_->getObjCoefficients();
        for (iRow = 0; iRow < numberRows; iRow++) {
          if (rowUpper[iRow] == 1.0 && rowLength[iRow] >= 5) {
            bool goodRow = true;
            bool overlap = false;
            int nObj = 0;
            for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
              int iColumn = column[j];
              if (elementByRow[j] != 1.0 || !originalModel_->isInteger(iColumn) || lower[iColumn]) {
                goodRow = false;
              }
              if (objective[iColumn])
                nObj++;
              if (mark[iColumn] >= 0) {
                overlap = true;
                numberOverlap++;
              }
            }
            if (nObj && nObj >= rowLength[iRow] - 1) {
              if (goodRow) {
                if (!overlap) {
                  // mark all
                  for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
                    int iColumn = column[j];
                    mark[iColumn] = numberSOS;
                  }
                  numberSOS++;
                  numberInSOS += rowLength[iRow];
                }
              } else {
                // no good
                allCostsInSOS = false;
              }
            }
          }
        }
        if (numberInSOS && allCostsInSOS) {
          int nGoodObj = 0;
          int nBadObj = 0;
          for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (objective[iColumn]) {
              if (mark[iColumn] >= 0)
                nGoodObj++;
              else
                nBadObj++;
            }
          }
          if (nBadObj * 10 < nGoodObj) {
            justOnesWithObj = true;
            makeEquality = 3;
#if CBC_USEFUL_PRINTING > 1
            printf("trying SOS as all costs there\n");
#endif
          }
        }
      }
    }
    if (firstOther < numberRows && makeEquality == 4) {
      CoinPackedMatrix *matrixByColumn = const_cast< CoinPackedMatrix * >(startModel_->getMatrixByCol());
      // Column copy
      const int *row = matrixByColumn->getIndices();
      const CoinBigIndex *columnStart = matrixByColumn->getVectorStarts();
      const int *columnLength = matrixByColumn->getVectorLengths();
      double *columnElements = matrixByColumn->getMutableElements();
      int *rowCount = new int[numberRows];
      memset(rowCount, 0, numberRows * sizeof(int));
      double *rowValue = new double[numberRows];
      int numberY = 0;
      CoinBigIndex numberElements = 0;
      int numberSOS = 0;
      for (int kRow = firstOther; kRow < numberRows; kRow++) {
        int iRow = whichRow[kRow];
        int n = 0;
        for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
          int iColumn = column[j];
          for (CoinBigIndex k = columnStart[iColumn]; k < columnStart[iColumn] + columnLength[iColumn];
               k++) {
            int jRow = row[k];
            double value = columnElements[k];
            if (jRow != iRow) {
              if (rowCount[jRow] > 0) {
                if (value != rowValue[jRow])
                  rowCount[jRow] = -1; // no good
                else
                  rowCount[jRow]++;
              } else if (!rowCount[jRow]) {
                whichRow[n++] = jRow;
                rowCount[jRow] = 1;
                rowValue[jRow] = value;
              }
            }
          }
        }
        int bestRow = -1;
        int bestCount = 4;
        for (int j = 0; j < n; j++) {
          int jRow = whichRow[j];
          int count = rowCount[jRow];
          rowCount[jRow] = 0;
          if (count >= 5) {
            numberY++;
            numberElements += count;
          }
          if (count > bestCount) {
            // possible
            bestRow = jRow;
            bestCount = count;
          }
        }
        if (bestRow >= 0) {
          numberSOS++;
          numberY++;
          numberElements += bestCount;
        }
      }
      if (numberY) {
        // Some may be duplicates
        // make sure ordered
        matrixByRow.orderMatrix();
        elementByRow = matrixByRow.getElements();
        column = matrixByRow.getIndices();
        rowStart = matrixByRow.getVectorStarts();
        rowLength = matrixByRow.getVectorLengths();
        CoinBigIndex *newStart = new CoinBigIndex[numberY + 1];
        int *newColumn = new int[numberElements];
        double *newValue = new double[numberElements];
        double *hash = new double[numberY];
        double *hashColumn = new double[numberColumns];
        int i;
        for (i = 0; i < numberColumns; i++)
          hashColumn[i] = randomGenerator.randomDouble();
        double *valueY = new double[3 * numberY];
        int *rowY = new int[3 * numberY];
        int *columnY = new int[3 * numberY];
        // For new solution
        double *newSolution = new double[numberColumns + numberY];
        memcpy(newSolution, startModel_->getColSolution(), numberColumns * sizeof(double));
        memset(rowCount, 0, numberRows * sizeof(int));
        // List of SOS entries to zero out
        CoinBigIndex *where = new CoinBigIndex[numberColumns];
        numberY = 0;
        numberElements = 0;
        int numberElementsY = 0;
        newStart[0] = 0;
        for (int kRow = firstOther; kRow < numberRows; kRow++) {
          int iRow = whichRow[kRow];
          int n = 0;
          int saveNumberY = numberY;
          for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
            int iColumn = column[j];
            for (CoinBigIndex k = columnStart[iColumn]; k < columnStart[iColumn] + columnLength[iColumn];
                 k++) {
              int jRow = row[k];
              double value = columnElements[k];
              if (jRow != iRow) {
                if (rowCount[jRow] > 0) {
                  if (value != rowValue[jRow])
                    rowCount[jRow] = -1; // no good
                  else
                    rowCount[jRow]++;
                } else if (!rowCount[jRow]) {
                  whichRow[n++] = jRow;
                  rowCount[jRow] = 1;
                  rowValue[jRow] = value;
                  assert(value);
                }
              }
            }
          }
          for (i = 0; i < n; i++) {
            // Sort so fewest first
            std::sort(whichRow, whichRow + n);
            int jRow = whichRow[i];
            int count = rowCount[jRow];
            rowCount[jRow] = 0;
            if (count >= 5) {
              //assert (count<rowLength[jRow]); // not error - just need to think
              // mark so not looked at again
              rowCount[jRow] = -count;
              // form new row
              double value = 0.0;
              double hashValue = 0.0;
              int nInSOS = 0;
              double valueOfY = 0.0;
              for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
                int iColumn = column[j];
                for (CoinBigIndex k = columnStart[iColumn]; k < columnStart[iColumn] + columnLength[iColumn];
                     k++) {
                  if (row[k] == jRow) {
                    value = columnElements[k];
                    newColumn[numberElements] = iColumn;
                    newValue[numberElements++] = 1.0;
                    hashValue += hashColumn[iColumn];
                    columnElements[k] = 0.0;
                    valueOfY += newSolution[iColumn];
                  } else if (row[k] == iRow) {
                    if (columnElements[k])
                      where[nInSOS++] = k;
                  }
                }
              }
              // See if already exists
              int n = static_cast< int >(numberElements - newStart[numberY]);
              int j;
              for (j = 0; j < numberY; j++) {
                if (hashValue == hash[j]) {
                  // Double check
                  CoinBigIndex offset = newStart[numberY] - newStart[j];
                  if (n == newStart[j + 1] - newStart[j]) {
                    CoinBigIndex k;
                    for (k = newStart[j]; k < newStart[j] + n; k++) {
                      if (newColumn[k] != newColumn[k + offset])
                        break;
                    }
                    if (k == newStart[j + 1])
                      break;
                  }
                }
              }
              if (j == numberY) {
                // not duplicate
                newSolution[numberY + numberColumns] = valueOfY;
                numberY++;
                newStart[numberY] = numberElements;
                hash[j] = hashValue;
                // Now do -1
                rowY[numberElementsY] = j + numberRows;
                columnY[numberElementsY] = j;
                valueY[numberElementsY++] = -1;
                if (n == nInSOS) {
                  // SOS entry
                  rowY[numberElementsY] = iRow;
                  columnY[numberElementsY] = j;
                  valueY[numberElementsY++] = 1;
                  for (int i = 0; i < n; i++) {
                    CoinBigIndex iEl = where[i];
                    columnElements[iEl] = 0.0;
                  }
                }
              } else {
                // duplicate
                numberElements = newStart[numberY];
              }
              // Now do
              rowY[numberElementsY] = jRow;
              columnY[numberElementsY] = j;
              valueY[numberElementsY++] = value;
            }
          }
          if (numberY > saveNumberY)
            rowCount[iRow] = -1000;
        }
        delete[] hash;
        delete[] hashColumn;
        matrixByColumn->cleanMatrix();
        // Now add rows
        double *rhs = new double[numberY];
        memset(rhs, 0, numberY * sizeof(double));
        startModel_->addRows(numberY, newStart, newColumn, newValue, rhs, rhs);
        delete[] rhs;
        delete[] newStart;
        delete[] newColumn;
        delete[] newValue;
        delete[] where;
        // Redo matrix
        CoinPackedMatrix add(true, rowY, columnY, valueY, numberElementsY);
        delete[] valueY;
        delete[] rowY;
        delete[] columnY;
        const int *row = add.getIndices();
        const CoinBigIndex *columnStart = add.getVectorStarts();
        //const int * columnLength = add.getVectorLengths();
        double *columnElements = add.getMutableElements();
        double *lo = new double[numberY];
        double *up = new double[numberY];
        for (i = 0; i < numberY; i++) {
          lo[i] = 0.0;
          up[i] = 1.0;
        }
        startModel_->addCols(numberY, columnStart, row, columnElements, lo, up, NULL);
        delete[] lo;
        delete[] up;
        for (i = 0; i < numberY; i++)
          startModel_->setInteger(i + numberColumns);
        CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(startModel_->getWarmStart());
        if (basis) {
          for (i = 0; i < numberY; i++) {
            basis->setArtifStatus(i + numberRows, CoinWarmStartBasis::atLowerBound);
            basis->setStructStatus(i + numberColumns, CoinWarmStartBasis::basic);
          }
          startModel_->setWarmStart(basis);
          delete basis;
        }
        startModel_->setColSolution(newSolution);
        delete[] newSolution;
        writeDebugMps(startModel_, "start", NULL);
        if (numberElements < 10 * std::min(numberColumns, 100 * numberY)) {
          handler_->message(CGL_ADDED_INTEGERS, messages_)
            << numberY << numberSOS << numberElements
            << CoinMessageEol;
          numberColumns += numberY;
          delete[] which;
          which = new int[numberColumns];
          bool saveTakeHint;
          OsiHintStrength saveStrength;
          startModel_->getHintParam(OsiDoDualInResolve,
            saveTakeHint, saveStrength);
          startModel_->setHintParam(OsiDoDualInResolve, false, OsiHintTry);
          startModel_->resolve();
          numberIterationsPre_ += startModel_->getIterationCount();
          startModel_->setHintParam(OsiDoDualInResolve, saveTakeHint, saveStrength);
        } else {
          // not such a good idea?
          delete startModel_;
          startModel_ = NULL;
        }
      }
      delete[] rowValue;
      delete[] rowCount;
    }
    if (makeEquality == 4) {
      makeEquality = 0;
#if 1
      // Try and make continuous variables integer
      // make clone
      if (!startModel_)
        startModel_ = originalModel_->clone();
      makeInteger();
#endif
    }
    delete[] whichRow;
    delete[] mark;
    if (numberSOS) {
      if (makeEquality == 2) {
        if (numberOverlap || numberIntegers * 4 > numberInSOS * 5 + 1) {
          handler_->message(CGL_PROCESS_SOS2, messages_)
            << numberSOS << numberInSOS << numberIntegers << numberOverlap
            << CoinMessageEol;
          makeEquality = 0;
        }
      }
    } else {
      // no sos
      makeEquality = 0;
    }
  }
  if (startModel_) {
    lower = originalModel_->getColLower();
    upper = originalModel_->getColUpper();
    rowLower = originalModel_->getRowLower();
    rowUpper = originalModel_->getRowUpper();
  }
  // See if all + 1
  bool allPlusOnes = true;
  int nPossible = 0;
  for (iRow = 0; iRow < numberRows; iRow++) {
    int numberP1 = 0, numberM1 = 0;
    int numberTotal = 0;
    CoinBigIndex j;
    double upperValue = rowUpper[iRow];
    double lowerValue = rowLower[iRow];
    bool good = true;
    bool possibleSlack = true;
    bool allPlus = true;
    for (j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
      int iColumn = column[j];
      double value = elementByRow[j];
      if (upper[iColumn] - lower[iColumn] < 1.0e-8) {
        // fixed
        upperValue -= lower[iColumn] * value;
        lowerValue -= lower[iColumn] * value;
        continue;
      } else if (!originalModel_->isBinary(iColumn)) {
        good = false;
        possibleSlack = false;
        //break;
      } else {
        numberTotal++;
      }

      if (fabs(value - floor(value + 0.5)) > 1.0e-12)
        possibleSlack = false;
      ;
      if (fabs(value) != 1.0) {
        good = false;
        allPlus = false;
      } else if (value > 0.0) {
        which[numberP1++] = iColumn;
      } else {
        numberM1++;
        which[numberColumns - numberM1] = iColumn;
        allPlus = false;
      }
    }
    if (possibleSlack) {
      if (upperValue > 1.0e20 && lowerValue > -1.0e12) {
        possibleSlack = (fabs(lowerValue - floor(lowerValue + 0.5)) < 1.0e-12);
      } else if (lowerValue < -1.0e20 && upperValue < 1.0e12) {
        possibleSlack = (fabs(upperValue - floor(upperValue + 0.5)) < 1.0e-12);
      } else {
        possibleSlack = false;
      }
    }
    if (allPlus)
      nPossible++;

    // avoid that we overflow in the int conversions below
    if (upperValue >= 1.0e6 || lowerValue <= -1.0e6)
      continue;

    int iUpper = static_cast<int>(floor(upperValue + 1.0e-5));
    int iLower = static_cast<int>(ceil(lowerValue - 1.0e-5));
    int state = 0;

    if (iUpper == 1 - numberM1)
      state = 1;
    else if (iUpper == -numberM1)
      state = 2;
    else if (iUpper < -numberM1)
      state = 3;
    if (fabs((static_cast< double >(iUpper)) - upperValue) > 1.0e-9)
      state = -1;

    if (!state) {
      if (-iLower == 1 - numberP1)
        state = -1;
      else if (-iLower == -numberP1)
        state = -2;
      else if (-iLower < -numberP1)
        state = -3;
      if (fabs((static_cast< double >(iLower)) - lowerValue) > 1.0e-9)
        state = -1;
    }

    if (good && state > 0) {
      if (abs(state) == 3) {
        // infeasible
        feasible = false;
        break;
      } else if (abs(state) == 2) {
        // we can fix all
        numberFixed += numberP1 + numberM1;
        int i;
        if (state > 0) {
          // fix all +1 at 0, -1 at 1
          for (i = 0; i < numberP1; i++)
            originalModel_->setColUpper(which[i], 0.0);
          for (i = 0; i < numberM1; i++)
            originalModel_->setColLower(which[numberColumns - i - 1], 1.0);
        } else {
          // fix all +1 at 1, -1 at 0
          for (i = 0; i < numberP1; i++)
            originalModel_->setColLower(which[i], 1.0);
          for (i = 0; i < numberM1; i++)
            originalModel_->setColUpper(which[numberColumns - i - 1], 0.0);
        }
      } else {
        if (!makeEquality || (makeEquality == -1 && numberM1 + numberP1 < minimumLength))
          continue;
        if (makeEquality == 2 || makeEquality == 3) {
          if (numberM1 || numberP1 < minimumLength)
            continue;
        }
        numberCliques++;
        if (iLower != iUpper) {
          element[numberSlacks] = state;
          rows[numberSlacks++] = iRow;
        }
        if (state > 0) {
          totalP1 += numberP1;
          totalM1 += numberM1;
        } else {
          totalP1 += numberM1;
          totalM1 += numberP1;
        }
      }
    }
  }
  // allow if some +1's
  allPlusOnes = 10 * nPossible > numberRows;
  delete[] which;
  if (!feasible) {
    handler_->message(CGL_INFEASIBLE, messages_)
      << CoinMessageEol;
    delete[] rows;
    delete[] element;
    return NULL;
  } else {
    if (numberCliques) {
      handler_->message(CGL_CLIQUES, messages_)
        << numberCliques
        << (static_cast< double >(totalP1 + totalM1)) / (static_cast< double >(numberCliques))
        << CoinMessageEol;
      //printf("%d of these could be converted to equality constraints\n",
      //     numberSlacks);
    }
    if (numberFixed)
      handler_->message(CGL_FIXED, messages_)
        << numberFixed
        << CoinMessageEol;
  }
  if (makeEquality == -2) {
    // do all
    numberSlacks = 0;
    // get row copy
    const CoinPackedMatrix * rowCopy = originalModel_->getMatrixByRow();
    const int * column = rowCopy->getIndices();
    const int * rowLength = rowCopy->getVectorLengths();
    const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
    const double * elementByRow = rowCopy->getElements();
    for (int iRow = 0; iRow < numberRows; iRow++) {
      // skip free rows - unlikely to have got here but ... !
      if (rowLower[iRow]<rowUpper[iRow] &&
	  (rowLower[iRow] > -1.0e30 || rowUpper[iRow] < 1.0e30)) {
	bool allInteger = true;
	if (rowLower[iRow]>-1.0e20&&
	    floor(rowLower[iRow]) != rowLower[iRow])
	  allInteger = false;
	if (rowUpper[iRow]<1.0e20&&
	    floor(rowUpper[iRow]) != rowUpper[iRow])
	  allInteger = false;
	CoinBigIndex kstart = rowStart[iRow];
	CoinBigIndex kend = rowStart[iRow] + rowLength[iRow];
	for (CoinBigIndex j = kstart; j < kend; j++) {
	  int iColumn = column[j];
	  if (allInteger) {
	    if (modelIn->isInteger(iColumn)) {
	      if (floor(elementByRow[j])!=elementByRow[j])
		allInteger = false;
	    } else {
	      allInteger = false;
	    }
	  }
	}
	rows[numberSlacks++] = allInteger ? iRow|0x80000000 : iRow;
      }
    }
  }
  if (numberSlacks && ((makeEquality == 2 && !justOnesWithObj)
		       || makeEquality == -2)) {
#ifdef CBC_KEEP_OLD_MESSAGES
    handler_->message(CGL_SLACKS, messages_)
      << numberSlacks
      << CoinMessageEol;
#endif
    // add variables to make equality rows
    // Get new model
    if (!startModel_) {
      startModel_ = originalModel_->clone();
    }
    CoinBigIndex *starts = new CoinBigIndex[numberSlacks + 1];
    double *loweretc = new double[3 * numberSlacks];
    double *upperetc = loweretc + numberSlacks;
    double *objetc = upperetc + numberSlacks;
    for (int i = 0; i < numberSlacks; i++) {
      int iRow = rows[i];
      double upperValue = 1.0;
      if (makeEquality == -2) {
	rows[numberSlacks+i] = iRow;
	int jRow = iRow & 0x7fffffff;
	upperValue = rowUpper[jRow] - rowLower[jRow];
	if (fabs(rowLower[jRow]) < fabs(rowUpper[jRow])) {
	  element[i] = -1.0;
	} else {
	  element[i] = 1.0;
	}
	iRow = jRow;
	if (upperValue>1.0e20)
	  upperValue = COIN_DBL_MAX;
      } else if (iRow >= numberRows) {
        // just a slack not a clique
        upperValue = COIN_DBL_MAX;
        iRow -= numberRows;
      }
      rows[i] = iRow;
      loweretc[i] = 0.0;
      upperetc[i] = upperValue;
      objetc[i] = 0.0;
      starts[i] = i;
    }
    starts[numberSlacks] = numberSlacks;
    startModel_->addCols(numberSlacks, starts, rows, element,
      loweretc, upperetc, objetc);
    delete[] starts;
    delete[] loweretc;
    CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(startModel_->getWarmStart());
    if (basis) {
      const double *rowActivity = startModel_->getRowActivity();
      double *solution = CoinCopyOfArray(startModel_->getColSolution(),
        numberColumns + numberSlacks);
      for (int i = 0; i < numberSlacks; i++) {
        int iRow = rows[i];
        if (basis->getArtifStatus(iRow) == CoinWarmStartBasis::basic) {
          basis->setArtifStatus(iRow, CoinWarmStartBasis::atLowerBound);
          basis->setStructStatus(i + numberColumns, CoinWarmStartBasis::basic);
        }
        double value = element[i];
        double slackvalue = rowActivity[iRow];
        if (value > 0) {
          slackvalue = rowUpper[iRow] - slackvalue;
        } else {
          slackvalue = slackvalue - rowLower[iRow];
        }
        solution[i + numberColumns] = slackvalue;
      }
      startModel_->setWarmStart(basis);
      startModel_->setColSolution(solution);
      delete basis;
      delete[] solution;
    }
    for (int i = 0; i < numberSlacks; i++) {
      int iRow = rows[i];
      double value = element[i];
      // set integer
      if (makeEquality != -2 || (rows[i+numberSlacks] &0x80000000) !=0)
	startModel_->setInteger(numberColumns + i);
      if (value > 0)
        startModel_->setRowLower(iRow, rowUpper[iRow]);
      else
        startModel_->setRowUpper(iRow, rowLower[iRow]);
    }
    if (prohibited_) {
      // extend
      int numberColumns = startModel_->getNumCols();
      char *temp = new char[numberColumns];
      memcpy(temp, prohibited_, numberProhibited_);
      memset(temp + numberProhibited_, 0, numberColumns - numberProhibited_);
      numberProhibited_ = numberColumns;
      delete[] prohibited_;
      prohibited_ = temp;
    }
  } else if (!startModel_) {
    // make clone anyway so can tighten bounds
    startModel_ = originalModel_->clone();
  }
  writeDebugMps(startModel_, "b1", NULL);
  // move objective to integers or to aggregated
  lower = startModel_->getColLower();
  upper = startModel_->getColUpper();
  rowLower = startModel_->getRowLower();
  rowUpper = startModel_->getRowUpper();
  matrixByRow = CoinPackedMatrix(*startModel_->getMatrixByRow());
  elementByRow = matrixByRow.getElements();
  column = matrixByRow.getIndices();
  rowStart = matrixByRow.getVectorStarts();
  rowLength = matrixByRow.getVectorLengths();
  numberRows = startModel_->getNumRows();
  numberColumns = startModel_->getNumCols();
  char *marked = new char[numberColumns];
  memset(marked, 0, numberColumns);
  // Column copy
#ifdef CBC_PREPROCESS_EXPERIMENT
  //CoinPackedMatrix * matrixByColumn = const_cast<CoinPackedMatrix *>(startModel_->getMatrixByCol());
  CoinPackedMatrix * matrixByColumn = startModel_->getMutableMatrixByCol();
  const int * row = matrixByColumn->getIndices();
  const CoinBigIndex * columnStart = matrixByColumn->getVectorStarts();
  const int * columnLength = startModel_->getMatrixByCol()->getVectorLengths();
  const double * columnElements = matrixByColumn->getElements();
  // look for costed slacks on equality rows first
  int nCosted = 0;
  int nCostedI = 0;
#define CAN_MOVE_LARGE_OBJ
#ifndef CAN_MOVE_LARGE_OBJ
  double goodRatio=0.001;
#else
  double goodRatio=0.0;
#endif
#endif
  double *obj = CoinCopyOfArray(startModel_->getObjCoefficients(), numberColumns);
#ifdef CBC_PREPROCESS_EXPERIMENT
#ifdef PRINT_STUFF
  {
    int * counts = new int[numberRows];
    memset(counts,0,numberRows*sizeof(int));
    int nn=0;
    int nz=0;
    int neq=0;
    for (int i=0;i<numberColumns;i++) {
      if (columnLength[i]==1) {
	int iRow = row[columnStart[i]];
	counts[iRow]++;
	if (rowLower[iRow]==rowUpper[iRow])
	  neq++;
	if (obj[i])
	  nn++;
	else
	  nz++;
      }
    }
    int cc[10];
    memset(cc,0,sizeof(cc));
    for (int i=0;i<numberRows;i++) {
      if (rowLower[i]==rowUpper[i]) {
	int k = std::min(counts[i],9);
	cc[k]++;
      }
    }
    for (int i=1;i<10;i++)
      if (cc[i])
	printf("(%d E rows have %d singletons) ",cc[i],i);
    printf("\n");
    memset(cc,0,sizeof(cc));
    for (int i=0;i<numberRows;i++) {
      if (rowLower[i]!=rowUpper[i]) {
	int k = std::min(counts[i],9);
	cc[k]++;
      }
    }
    for (int i=1;i<10;i++)
      if (cc[i])
	printf("(%d L/G rows have %d singletons) ",cc[i],i);
    printf("\n");
    delete [] counts;
    printf("%d single with cost %d zero cost %d equality rows\n",
	   nn,nz,neq);
  }
#endif
#endif
  double offset;
#if CBC_USEFUL_PRINTING > 1
  int numberMoved = 0;
#endif
  startModel_->getDblParam(OsiObjOffset, offset);
#if CBC_PREPROCESS_EXPERIMENT > 1
  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
    if (columnLength[iColumn]==1 && obj[iColumn]) {
      double testObj = fabs(obj[iColumn])*goodRatio;
      int iRow = row[columnStart[iColumn]];
      if (rowLower[iRow]==rowUpper[iRow]) {
	double rhs = rowLower[iRow];
	bool canDo = false;
        if (startModel_->isInteger(iColumn)) {
	  // are all integer
	  bool allInteger=true;
	  double element = columnElements[columnStart[iColumn]];
	  if ((element-floor(element)))
	    allInteger=false;
	  for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
	    int jColumn = column[j];
	    if (!startModel_->isInteger(jColumn)) {
	      allInteger=false;
	      break;
	    } else if (fabs(obj[jColumn])<testObj) {
	      // scaling
	      allInteger=false;
	      break;
	    }
	  }
	  if (allInteger) {
	    //printf("can take off cost on integer %d and make slack\n",iColumn);
	    nCostedI++;
	    canDo = true;
	    marked[iColumn]=1;
	    //printf("zz %d\n",iColumn);
	  }
	} else {
	  bool goodScaling=true;
	  for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
	    int jColumn = column[j];
	    if (fabs(obj[jColumn])<testObj) {
	      // scaling
	      goodScaling=false;
	      break;
	    }
	  }
	  if (goodScaling) {
	    //printf("can take off cost on %d and make slack\n",iColumn);
	    nCosted++;
	    marked[iColumn]=1;
	    canDo = true;
	    //printf("zz %d\n",iColumn);
	  }
	}
	if (canDo) {
	  // make cost zero
	  double valueSlack = columnElements[columnStart[iColumn]];
	  double objValue = obj[iColumn];
	  double subtract = objValue/valueSlack;
	  offset -= subtract * rhs;
	  startModel_->setContinuous(iColumn); // to allow change 
	  for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
	    int jColumn = column[j];
	    if (iColumn != jColumn) {
	      double value = elementByRow[j];
	      obj[jColumn] -= subtract*value;
	    } else {
	      obj[iColumn] = 0.0;
	    }
	  }
	}
      }
    }
  }
  if (nCosted||nCostedI) {
    startModel_->setDblParam(OsiObjOffset, offset);
    startModel_->setObjective(obj);
#ifdef PRINT_STUFF
    printf("%d integer costed slacks and %d continuous\n",nCostedI,nCosted);
    int nn=0;
    int nz=0;
    for (int i=0;i<numberColumns;i++) {
      if (columnLength[i]==1) {
	if (obj[i])
	  nn++;
	else
	  nz++;
      }
    }
    printf("%d single with cost %d zero cost\n",nn,nz);
#endif
  }
#endif
  writeDebugMps(startModel_, "b2", NULL);
  // This is not a vital loop so be careful
#ifndef SKIP_MOVE_COSTS_TO_INTEGERS
  for (iRow = 0; iRow < numberRows; iRow++) {
    int nPlus = 0;
    int nMinus = 0;
    int iPlus = -1;
    int iMinus = -1;
    double valuePlus = 0;
    double valueMinus = 0;
    //bool allInteger=true;
    double rhs = rowLower[iRow];
    if (rhs != rowUpper[iRow])
      continue;
    int numberContinuous = 0;
    for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
      int iColumn = column[j];
      assert (iColumn<numberColumns);
#ifdef CBC_PREPROCESS_EXPERIMENT
      if (marked[iColumn]) {
	nPlus = -numberColumns;
	nMinus = -numberColumns;
	break;
      }
#endif
      double value = elementByRow[j];
      if (upper[iColumn] > lower[iColumn]) {
        if (startModel_->isInteger(iColumn)) {
#if 0
	  if (columnLength[iColumn]==1) {
	    if (value==1.0) {
	    }
	  }
	  if (value!=floor(value+0.5))
	    allInteger=false;
	  if (allInteger&&fabs(value)<1.0e8) {
	    if (!multiple)
	      multiple = static_cast<int> (fabs(value));
	    else if (multiple>0)
	      multiple = gcd(multiple,static_cast<int> (fabs(value)));
	  } else {
	    allInteger=false;
	  }
#endif
        } else {
          numberContinuous++;
        }
        if (value > 0.0) {
          if (nPlus > 0 && value != valuePlus) {
            nPlus = -numberColumns;
          } else if (!nPlus) {
            nPlus = 1;
            iPlus = iColumn;
            valuePlus = value;
          } else {
            nPlus++;
          }
        } else {
          if (nMinus > 0 && value != valueMinus) {
            nMinus = -numberColumns;
          } else if (!nMinus) {
            nMinus = 1;
            iMinus = iColumn;
            valueMinus = value;
          } else {
            nMinus++;
          }
        }
      }
    }
    if (((nPlus == 1 && startModel_->isInteger(iPlus) && nMinus > 0 && !obj[iPlus]) ||
	 (nMinus == 1 && startModel_->isInteger(iMinus) && nPlus > 0 && !obj[iMinus]))
	&& numberContinuous) {
      int jColumn;
      double multiplier;
      double inverseValueOther;
      if (nPlus == 1) {
        jColumn = iPlus;
        multiplier = fabs(valuePlus / valueMinus);
        rhs /= -valueMinus;
	inverseValueOther = 1.0/valueMinus;
      } else {
        jColumn = iMinus;
        multiplier = fabs(valueMinus / valuePlus);
        rhs /= valuePlus;
	inverseValueOther = 1.0/valuePlus;
      }
      double smallestPos = COIN_DBL_MAX;
      double smallestNeg = -COIN_DBL_MAX;
      for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
        int iColumn = column[j];
        if (iColumn != jColumn) {
          double objValue = obj[iColumn];
          if (upper[iColumn] > lower[iColumn]) {
            if (objValue >= 0.0)
              smallestPos = std::min(smallestPos, objValue);
            else
              smallestNeg = std::max(smallestNeg, objValue);
          }
        }
      }
      if (smallestPos > 0.0) {
        double move = 0.0;
        double multiply = (smallestNeg == -COIN_DBL_MAX) ? -1.0 : 1.0;
        if (smallestNeg == -COIN_DBL_MAX)
          move = smallestPos;
        else if (smallestPos == COIN_DBL_MAX)
          move = smallestNeg;
        if (move) {
          // can move objective
#if CBC_USEFUL_PRINTING > 1
          numberMoved++;
          if (rhs)
            printf("ZZZ on col %d move %g offset %g\n",
              jColumn, move, move * rhs);
#endif
          offset += move * multiply * rhs;
          for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
            int iColumn = column[j];
            if (iColumn != jColumn) {
              if (upper[iColumn] > lower[iColumn]) {
                obj[iColumn] -= move;
	      } else {
		double value = elementByRow[j];
		obj[iColumn] -= move*value*inverseValueOther;
	      }
            } else {
              obj[jColumn] += move * multiplier;
            }
          }
        }
      }
    }
  }
#endif
#if CBC_USEFUL_PRINTING > 1
  if (numberMoved)
    printf("ZZZ %d costs moved\n", numberMoved);
#endif
  startModel_->setDblParam(OsiObjOffset, offset);
  startModel_->setObjective(obj);
  delete[] obj;
  delete[] marked;
  delete[] rows;
  delete[] element;
  writeDebugMps(startModel_, "b3", NULL);
  if (makeIntegers) {
    makeIntegers2(startModel_, makeIntegers);
  }
  writeDebugMps(startModel_, "b4", NULL);
  int infeas = 0;
  OsiSolverInterface *startModel2 = startModel_;
  // Do we want initial presolve
  if (doInitialPresolve) {
    assert(doInitialPresolve == 1);
    OsiSolverInterface *presolvedModel;
    OsiSolverInterface *oldModel = startModel2;
    OsiPresolve *pinfo = new OsiPresolve();
    int presolveActions = 0;
    // Allow dual stuff on integers
    // Allow stuff which may not unroll cleanly - unless told not to
    if ((tuning & 4096) == 0)
      presolveActions = 1 + 16;
    else
      presolveActions = 16; // actually just switch off duplicate columns for ints
    if ((tuning & 32) != 0)
      presolveActions |= 32;
    if ((tuning & 512) != 0)
      presolveActions |= 0x4000;
    // Do not allow all +1 to be tampered with
    //if (allPlusOnes)
    //presolveActions |= 2;
    // allow transfer of costs
    presolveActions |= 4; // can be slow
#ifdef CBC_PREPROCESS_EXPERIMENT
    if ((tuning&8192)!=0)
      presolveActions |= 0x40 | 0x80; // more for dupcol and doubleton
#endif
    // If trying for SOS don't allow some transfers
    if (makeEquality == 2 || makeEquality == 3)
      presolveActions |= 8;
    pinfo->setPresolveActions(presolveActions);
    if (prohibited_)
      assert(numberProhibited_ == oldModel->getNumCols());
    int saveLogLevel = oldModel->messageHandler()->logLevel();
    if (saveLogLevel == 1)
      oldModel->messageHandler()->setLogLevel(0);
    std::string solverName;
    oldModel->getStrParam(OsiSolverName, solverName);
    // Extend if you want other solvers to keep solution
    bool keepSolution = solverName == "clp";
    //double feasibilityTolerance = ((tuning & 1024) == 0) ? CGL_PREPROCESS_TOLERANCE : 1.0e-4;
    presolvedModel = pinfo->presolvedModel(*oldModel, feasibilityTolerance, true, 5, prohibited_, keepSolution, rowType_,scBound);
    oldModel->messageHandler()->setLogLevel(saveLogLevel);
    if (presolvedModel) {
      //#define MAKE_LESS_THAN
#ifdef MAKE_LESS_THAN
      {
        int numberRows = presolvedModel->getNumRows();
        int numberColumns = presolvedModel->getNumCols();
        CoinPackedMatrix *matrix = presolvedModel->getMutableMatrixByCol();
        const int *row = matrix->getIndices();
        const CoinBigIndex *columnStart = matrix->getVectorStarts();
        const int *columnLength = matrix->getVectorLengths();
        double *element = const_cast< double * >(matrix->getElements());
        const double *rowLower = presolvedModel->getRowLower();
        const double *rowUpper = presolvedModel->getRowUpper();
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          for (CoinBigIndex j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            if (rowUpper[iRow] == COIN_DBL_MAX)
              element[j] = -element[j];
          }
        }
        for (int iRow = 0; iRow < numberRows; iRow++) {
          if (rowUpper[iRow] == COIN_DBL_MAX) {
            presolvedModel->setRowUpper(iRow, -rowLower[iRow]);
            presolvedModel->setRowLower(iRow, -COIN_DBL_MAX);
          }
        }
      }
#endif
      presolvedModel->messageHandler()->setLogLevel(saveLogLevel);
      //presolvedModel->writeMps("new");
      writeDebugMps(presolvedModel, "ordinary", pinfo);
      // update prohibited and rowType
      update(pinfo, presolvedModel,scBound);
      if (!presolvedModel->getNumRows()) {
        doInitialPresolve = 0;
        delete presolvedModel;
        delete pinfo;
      } else {
        model_[0] = presolvedModel;
        presolve_[0] = pinfo;
        modifiedModel_[0] = presolvedModel->clone();
        startModel2 = modifiedModel_[0];
      }
    } else {
      infeas = 1;
      doInitialPresolve = 0;
      delete presolvedModel;
      delete pinfo;
    }
  }
  // tighten bounds
  /*

  Virtuous solvers may require a refresh via initialSolve if this
  call is ever changed to give a nonzero value to the (default) second
  parameter. Previous actions may have made significant changes to the
  constraint system. Safe as long as tightenPrimalBounds doesn't ask for
  the current solution.
*/
  if (!infeas) {
    // may be better to just do at end
    writeDebugMps(startModel2, "before", NULL);
    infeas = tightenPrimalBounds(*startModel2,false,scBound);
    infeas = (infeas<0) ? 1 : 0;
    writeDebugMps(startModel2, "after", NULL);
    // make as many integer as possible
    int numberChanged = analyze(startModel2);
    if (numberChanged<0)
      infeas = 1;
    else if (numberChanged)
      handler_->message(CGL_MADE_INTEGER, messages_)
	<< numberChanged
	<< CoinMessageEol;
  }
  if (infeas) {
    handler_->message(CGL_INFEASIBLE, messages_)
      << CoinMessageEol;
    return NULL;
  } else {
    //printf("skipping tightenPrimalBounds\n");
  }
  OsiSolverInterface *returnModel = NULL;
  int numberChanges;
  if ((tuning & (128 + 1)) != 0) {
    OsiSolverInterface *newSolver = NULL;
    if ((tuning & 2048) != 0) {
      CglProbing generator1;
      generator1.setUsingObjective(false);
      generator1.setMaxPass(1);
      generator1.setMaxPassRoot(1);
      generator1.setMaxLook(100);
      generator1.setRowCuts(3);
      generator1.setMaxElements(300);
      generator1.setMaxProbeRoot(startModel2->getNumCols());
      CoinThreadRandom randomGenerator;
      CglTreeProbingInfo info(startModel2);
      info.level = 0;
      info.formulation_rows = startModel2->getNumRows();
      info.inTree = false;
      info.options = !numberProhibited_ ? 0 : 2;
      info.randomNumberGenerator = &randomGenerator;
      info.pass = 4;
      generator1.setMode(8);
      OsiCuts cs;
      generator1.generateCutsAndModify(*startModel2, cs, &info);
      OsiSolverInterface *temp = generator1.cliqueModel(startModel2, 2);
      OsiSolverInterface *temp2 = cliqueIt(*temp, 0.0001);
      delete temp;
      if (temp2) {
        if (temp2->getNumRows() < startModel2->getNumRows()) {
          delete newSolver;
          newSolver = temp2;
        } else {
          delete temp2;
        }
      }
    } else {
      // take out cliques
      newSolver = cliqueIt(*startModel2, 0.0001);
    }
    if (newSolver) {
      if (startModel2 == modifiedModel_[0])
        modifiedModel_[0] = newSolver;
      delete startModel2;
      startModel2 = newSolver;
      newSolver->setHintParam(OsiDoDualInInitial, true, OsiHintTry);
      newSolver->initialSolve();
      assert(newSolver->isProvenOptimal());
      //printf("new size %d rows, %d columns\n",
      //     newSolver->getNumRows(),newSolver->getNumCols());
    }
  }
  {
    // Give a hint to do dual
    bool saveTakeHint;
    OsiHintStrength saveStrength;
    startModel2->getHintParam(OsiDoDualInInitial,
      saveTakeHint, saveStrength);
    startModel2->setHintParam(OsiDoDualInInitial, true, OsiHintTry);
    startModel2->initialSolve();
    numberIterationsPre_ += startModel2->getIterationCount();
    // double check
    if (!startModel2->isProvenOptimal()) {
      if (!startModel2->isProvenDualInfeasible()) {
	// relax any fixed (or almost fixed) continuous variables
	if (numberColumns==startModel2->getNumCols()) {
	  const double * originalLower = model.getColLower();
	  const double * originalUpper = model.getColUpper();
	  const double * lower = startModel2->getColLower();
	  const double * upper = startModel2->getColUpper();
	  int numberColumns = model.getNumCols();
	  // will need different coding on later infeas?
	  for (int i=0;i<numberColumns;i++) {
	    if (!startModel2->isInteger(i)&&
		upper[i]-lower[i]<feasibilityTolerance) {
	      // relax bounds
	      double lo = std::max(originalLower[i],lower[i]-feasibilityTolerance);
	      double up = std::min(originalUpper[i],upper[i]+feasibilityTolerance);
	      startModel2->setColLower(i,lo);
	      startModel2->setColUpper(i,up);
	    }
	  }
	}
        // Do presolves
        bool saveHint;
        OsiHintStrength saveStrength;
        startModel2->getHintParam(OsiDoPresolveInInitial, saveHint, saveStrength);
        startModel2->setHintParam(OsiDoPresolveInInitial, true, OsiHintTry);
        startModel2->setHintParam(OsiDoDualInInitial, false, OsiHintTry);
        startModel2->initialSolve();
        numberIterationsPre_ += startModel2->getIterationCount();
        if (!startModel2->isProvenDualInfeasible()) {
          CoinWarmStart *empty = startModel2->getEmptyWarmStart();
          startModel2->setWarmStart(empty);
          delete empty;
          startModel2->setHintParam(OsiDoDualInInitial, true, OsiHintTry);
          startModel2->initialSolve();
          numberIterationsPre_ += startModel2->getIterationCount();
        }
        startModel2->setHintParam(OsiDoPresolveInInitial, saveHint, saveStrength);
      }
    }
    startModel2->setHintParam(OsiDoDualInInitial, saveTakeHint, saveStrength);
  }
  if (!startModel2->isProvenOptimal()) {
    if (!startModel2->isProvenDualInfeasible()) {
      handler_->message(CGL_INFEASIBLE, messages_) << CoinMessageEol;
#if CBC_USEFUL_PRINTING > 1
      startModel2->writeMps("infeas");
#endif
    } else {
      handler_->message(CGL_UNBOUNDED, messages_) << CoinMessageEol;
    }
    return NULL;
  }
  reducedCostFix(*startModel2);
  if (!numberSolvers_) {
    // just fix
    OsiSolverInterface *newModel = modified(startModel2, false, numberChanges, 0, numberModifiedPasses);
    if (startModel_ != originalModel_)
      delete startModel_;
    if (startModel2 != startModel_)
      delete startModel2;
    startModel_ = newModel;
    returnModel = startModel_;
  } else {
    OsiSolverInterface *presolvedModel;
    OsiSolverInterface *oldModel = startModel2;
    if (doInitialPresolve)
      oldModel = modifiedModel_[0];
      //CglDuplicateRow dupCuts(oldModel);
      //dupCuts.setLogLevel(1);
      // If +1 try duplicate rows
#define USECGLCLIQUE 512
    if ((options_ & 8) != 0)
      tuning &= ~USECGLCLIQUE;
    if ((options_ & 4) != 0)
      allPlusOnes = false;
    if ((allPlusOnes && (options_ & 8) == 0) || (tuning & USECGLCLIQUE) != 0) {
#if 1
      // put at beginning
      int nAdd = ((tuning & (64 + USECGLCLIQUE)) == 64 + USECGLCLIQUE && allPlusOnes) ? 2 : 1;
      CglCutGenerator **temp = generator_;
      generator_ = new CglCutGenerator *[numberCutGenerators_ + nAdd];
      if (temp)
        memcpy(generator_ + nAdd, temp, numberCutGenerators_ * sizeof(CglCutGenerator *));
      delete[] temp;
      numberCutGenerators_ += nAdd;
      if (nAdd == 2 || (tuning & USECGLCLIQUE) != 0) {
        if ((options_ & 1024) == 0) {
          // Modern cgraph-based clique separator (default)
          CglBKClique *bkGen = new CglBKClique();
          generator_[0] = bkGen;
        } else {
          // Legacy CglClique (options_ & 1024 to force old behavior)
          CglClique *cliqueGen = new CglClique(false, true);
          cliqueGen->setStarCliqueReport(false);
          cliqueGen->setRowCliqueReport(false);
          if ((tuning & USECGLCLIQUE) == 0)
            cliqueGen->setMinViolation(-2.0);
          else
            cliqueGen->setMinViolation(-3.0);
          generator_[0] = cliqueGen;
        }
      }
      if ((allPlusOnes || (tuning & 256) != 0) && (options_ & 4) == 0) {
        CglDuplicateRow *dupCuts = new CglDuplicateRow(oldModel);
        if ((tuning & 256) != 0)
          dupCuts->setMaximumDominated(numberColumns);
        generator_[nAdd - 1] = dupCuts;
      }
#else
      CglDuplicateRow dupCuts(oldModel);
      addCutGenerator(&dupCuts);
#endif
    }
#if 1
    // Dominated columns
    if ((options_&512)!=0) {
    // Column copy
      const CoinPackedMatrix * matrixByCol = oldModel->getMatrixByCol();
      const double * element = matrixByCol->getElements();
      const int * row = matrixByCol->getIndices();
      const CoinBigIndex * columnStart = matrixByCol->getVectorStarts();
      const int * columnLength = matrixByCol->getVectorLengths();
      
      const double * rowLower = oldModel->getRowLower();
      const double * rowUpper = oldModel->getRowUpper();
      const double * columnLower = oldModel->getColLower();
      const double * columnUpper = oldModel->getColUpper();
      const double * cost = oldModel->getObjCoefficients();
      // Need to do something about implied upper bounds
      // Row copy
      const CoinPackedMatrix * matrixByRow = oldModel->getMatrixByRow();
      const double * elementByRow = matrixByRow->getElements();
      const int * column = matrixByRow->getIndices();
      const CoinBigIndex * rowStart = matrixByRow->getVectorStarts();
      const int * rowLength = matrixByRow->getVectorLengths();
      int numberRows = oldModel->getNumRows();
      int numberColumns = oldModel->getNumCols();
      typedef struct {
	unsigned int isInteger : 1;
	unsigned int columnType : 2;
	unsigned int column : 29;
      } variableType;
      // candidate columns - just ones in L or G rows
      variableType * possible = new variableType [numberColumns];
      int * rowSet = new int [2*numberRows];
      /* 0 - no good
	 > 0 no need to flip if 2 then implied bounds OK
	 < 0 flip and -2 */
      int * flip = rowSet+numberRows;
      for (int iRow=0;iRow < numberRows;iRow++) {
	flip[iRow]=1;
	if (rowUpper[iRow]<COIN_DBL_MAX && rowLower[iRow]>-COIN_DBL_MAX)
	  flip[iRow]=0;
	else if (rowUpper[iRow]==COIN_DBL_MAX)
	  flip[iRow]=-1;
	// sure we could do much better
	bool allOK = flip[iRow] != 0;
	for (CoinBigIndex i = rowStart[iRow];
	     i < rowStart[iRow] + rowLength[iRow]; i++) {
	  int iColumn = column[i];
	  if (!oldModel->isInteger(iColumn)) {
	    allOK = false;
	    break;
	  }
	}
	if (allOK) {
	  if (flip[iRow]>0) {
	    // see if implied bounds OK
	    for (CoinBigIndex i = rowStart[iRow];
		 i < rowStart[iRow] + rowLength[iRow]; i++) {
	      int iColumn = column[i];
	      double value = elementByRow[i];
	      if (value>0.0) {
		allOK = false;
		break;
	      } else if (value*columnUpper[iColumn]>rowUpper[iRow]) {
		allOK = false;
		break;
	      }
	    }
	    if (allOK)
	      flip[iRow] = 2;
	  } else {
	    // see if implied bounds OK
	    for (CoinBigIndex i = rowStart[iRow];
		 i < rowStart[iRow] + rowLength[iRow]; i++) {
	      int iColumn = column[i];
	      double value = elementByRow[i];
	      if (value<0.0) {
		allOK = false;
		break;
	      } else if (value*columnUpper[iColumn]<rowLower[iRow]) {
		allOK = false;
		break;
	      }
	    }
	    if (allOK)
	      flip[iRow] = -2;
	  }
	}
      }
      int nPossible = 0;
      for (int iColumn=0;iColumn < numberColumns;iColumn++) {
	if (columnLower[iColumn]==0.0 &&
	    (!prohibited_ || !prohibited_[iColumn])) {
	  bool ok=true;
	  double upper = columnUpper[iColumn];
	  // 0 - all good, 1 - cost likes zero, 2 cost does not like 0
	  int columnType = 0;
	  bool checkFurther = false;
	  for (CoinBigIndex i = columnStart[iColumn];
	       i < columnStart[iColumn] + columnLength[iColumn]; i++) {
	    int iRow = row[i];
	    double value = element[i];
	    if (!flip[iRow]) {
	      ok=false;
	      break;
	    } else if (flip[iRow]<0) {
	      value = -value;
	    }
	    if (abs(flip[iRow])==1) {
	      checkFurther = true;
	    }
	    // clean this up
	    // may need min/max on rows
	    if (value < 0.0 && abs(flip[iRow])==1) {
	      if (upper<1.0e20) {
		ok=false;
		break;
	      }
	    }
	  }
	  if (ok && checkFurther) 
	    columnType = (cost[iColumn]<0) ? 2 : 1;
	  if (ok) {
	    variableType thisColumn = {0};
	    thisColumn.column = iColumn;
	    if (oldModel->isInteger(iColumn))
	      thisColumn.isInteger = 1;
	    thisColumn.columnType = columnType;
	    possible[nPossible++] = thisColumn;
	  }
	}
      }
      unsigned int * test = new unsigned int [nPossible];
      double * compEls = new double [numberRows];
      memset(compEls,0,numberRows*sizeof(double));
      int nDominated = 0;
      for (int j=0;j<nPossible;j++) {
	int iColumn = possible[j].column;
	unsigned int thisTest = 0;
	int nEls = columnLength[iColumn];
	int nNeg = 0;
	for (int k = 0;k<nEls;k++) {
	  CoinBigIndex i = columnStart[iColumn]+k;
	  int iRow = row[i];
	  double value = element[i];
	  if (flip[iRow]<0)
	    value = -value;
	  int iShift = iRow%32;
	  thisTest |= 1<<iShift;
	  compEls[iRow]=value;
	  if (value<0.0)
	    nNeg++;
	  rowSet[k]=iRow;
	}
	int nPos = nEls-nNeg;
	test[j] = thisTest;
	// may be able to do better if sort?
	// test all before
	int isInteger = possible[j].isInteger;
	for (int j2=0;j2<j;j2++) {
	  unsigned int test2 = test[j2];
	  unsigned int testBoth = test2|thisTest;
	  int jColumn = possible[j2].column;
	  if (possible[j].columnType==3)
	    continue;
	  if ((testBoth==test2 || testBoth==thisTest)
	      && isInteger == possible[j2].isInteger) {
	    // possible - maybe see which one is possible
	    int choose = 0;
	    if (cost[iColumn]<cost[jColumn])
	      choose = 1;
	    else if (cost[iColumn]>cost[jColumn])
	      choose = -1;
	    int nNegUnmatched = nNeg;
	    int nPosUnmatched = nPos;
	    for (CoinBigIndex k = columnStart[jColumn];
		 k < columnStart[jColumn] + columnLength[jColumn]; k++) {
	      int iRow = row[k];
	      double value = element[k];
	      if (flip[iRow]<0)
		value = - value;
	      if (compEls[iRow]<0.0)
		nNegUnmatched--;
	      else if (compEls[iRow]>0.0)
		nPosUnmatched--;
	      if (value<compEls[iRow]) {
		if (!choose) {
		  choose = -1;
		} else if (choose>0) {
		  choose = 2;
		  break;
		}
	      } else if (value>compEls[iRow]) {
		if (!choose) {
		  choose = 1;
		} else if (choose<0) {
		  choose = 2;
		  break;
		}
	      }
	    }
	    if (nNegUnmatched&&nPosUnmatched)
	      choose=2;
	    if (choose!=2) {
	      if (choose==0) {
		if (nNegUnmatched)
		  choose=1;
		else if (nPosUnmatched)
		  choose=-1;
	      } else if (choose==1) {
		if (nPosUnmatched)
		  choose=2;
	      } else {
		if (nNegUnmatched)
		  choose=2;
	      }
	      if (choose==1 && possible[j2].columnType<2) {
		// fix jColumn to zero
		possible[j2].columnType=3;
		oldModel->setColUpper(jColumn,0.0);
		nDominated++;
	      } else if (choose==-1 && possible[j].columnType<2) {
		// fix iColumn to zero
		possible[j].columnType=3;
		oldModel->setColUpper(iColumn,0.0);
		nDominated++;
		break;
	      } else if (choose==0 && possible[j2].columnType<2) {
		// fix jColumn to zero
		possible[j2].columnType=3;
		oldModel->setColUpper(jColumn,0.0);
		nDominated++;
	      }
	    }
	  }
	}
	for (int k = 0;k<nEls;k++) {
	  int iRow = rowSet[k];
	  compEls[iRow]=0.0;
	}
      }
#ifdef CBC_KEEP_OLD_MESSAGES
      if (nDominated) {
	char generalPrint[100];
	sprintf(generalPrint, "%d variables fixed as dominated",
		nDominated);
	handler_->message(CGL_GENERAL, messages_)
	  << generalPrint
	  << CoinMessageEol;
      }
#endif
      delete [] compEls;
      delete [] test;
      delete [] rowSet;
      delete [] possible;
    }
#endif
    for (int iPass = doInitialPresolve; ((iPass < numberSolvers_) && ((getCurrentCPUTime() - ppstart) < timeLimit_)); iPass++) {
      if (inspect_) {
        int nR = oldModel->getNumRows();
        int nC = oldModel->getNumCols();
        int nNZ = oldModel->getNumElements();
        int nInts = 0;
        for (int i = 0; i < nC; i++)
          if (oldModel->isInteger(i))
            nInts++;
        char buf[256];
        snprintf(buf, sizeof(buf),
          "[Preproc major pass %d/%d] %d rows, %d cols (%d integer), %d NZ  (t=%.2fs)",
          iPass - doInitialPresolve + 1, numberSolvers_ - doInitialPresolve,
          nR, nC, nInts, nNZ, getCurrentCPUTime() - ppstart);
        FILE *fp = handler_->filePointer();
        if (fp) { fprintf(fp, "  %s\n", buf); fflush(fp); }
      }
      // Look at Vubs
      {
        const double *columnLower = oldModel->getColLower();
        const double *columnUpper = oldModel->getColUpper();
        const CoinPackedMatrix *rowCopy = oldModel->getMatrixByRow();
        const int *column = rowCopy->getIndices();
        const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
        const int *rowLength = rowCopy->getVectorLengths();
        const double *rowElements = rowCopy->getElements();
        const CoinPackedMatrix *columnCopy = oldModel->getMatrixByCol();
        //const int * row = columnCopy->getIndices();
        //const CoinBigIndex * columnStart = columnCopy->getVectorStarts();
        const int *columnLength = columnCopy->getVectorLengths();
        //const double * columnElements = columnCopy->getElements();
        const double *rowLower = oldModel->getRowLower();
        const double *rowUpper = oldModel->getRowUpper();
        const double *objective = oldModel->getObjCoefficients();
        double direction = oldModel->getObjSense();
        int numberRows = oldModel->getNumRows();
        for (int iRow = 0; iRow < numberRows; iRow++) {
          if (rowLength[iRow] == 2 && (rowLower[iRow] < -1.0e20 || rowUpper[iRow] > 1.0e20)) {
            CoinBigIndex start = rowStart[iRow];
            int iColumn1 = column[start];
            int iColumn2 = column[start + 1];
            double value1 = rowElements[start];
            double value2 = rowElements[start + 1];
            double upper;
            if (rowLower[iRow] < -1.0e20) {
              if (rowUpper[iRow] < 1.0e20)
                upper = rowUpper[iRow];
              else
                continue; // free row
            } else {
              upper = -rowLower[iRow];
              value1 = -value1;
              value2 = -value2;
            }
            //for now just singletons
            bool integer1 = oldModel->isInteger(iColumn1);
            bool integer2 = oldModel->isInteger(iColumn2);
            int debug = 0;
            if (columnLength[iColumn1] == 1) {
              if (integer1) {
                debug = 0; // no good
              } else if (integer2) {
                // possible
                debug = 1;
              }
            } else if (columnLength[iColumn2] == 1) {
              if (integer2) {
                debug = -1; // print and skip
              } else if (integer1) {
                // possible
                debug = 1;
                double valueD = value1;
                value1 = value2;
                value2 = valueD;
                int valueI = iColumn1;
                iColumn1 = iColumn2;
                iColumn2 = valueI;
                bool valueB = integer1;
                integer1 = integer2;
                integer2 = valueB;
              }
            }
#if CBC_USEFUL_PRINTING
            if (debug && 0) {
              printf("%d %d elements%selement %g and %d %d elements%selement %g <= %g\n",
                iColumn1, columnLength[iColumn1], integer1 ? " (integer) " : " ", value1,
                iColumn2, columnLength[iColumn2], integer2 ? " (integer) " : " ", value2,
                upper);
            }
#endif
            if (debug > 0) {
              if (value1 > 0.0 && objective[iColumn1] * direction < 0.0) {
                // will push as high as possible so make ==
                // highest effective rhs
                if (value2 > 0)
                  upper -= value2 * columnLower[iColumn2];
                else
                  upper -= value2 * columnUpper[iColumn2];
                if (columnUpper[iColumn1] > 1.0e20 || columnUpper[iColumn1] * value1 >= upper) {
                  //printf("looks possible\n");
                  // make equality
                  if (rowLower[iRow] < -1.0e20)
                    oldModel->setRowLower(iRow, rowUpper[iRow]);
                  else
                    oldModel->setRowUpper(iRow, rowLower[iRow]);
                } else {
                  // may be able to make integer
                  // may just be better to use to see objective integral
                  if (upper == floor(upper) && value2 == floor(value2) && value1 == floor(value1) && objective[iColumn1] == floor(objective[iColumn1]))
                    oldModel->setInteger(iColumn1);
                  //printf("odd3\n");
                }
              } else if (value1 < 0.0 && objective[iColumn1] * direction > 0.0) {
                //printf("odd4\n");
              } else {
                //printf("odd2\n");
              }
            } else if (debug < 0) {
              //printf("odd1\n");
            }
          }
        }
      }
      OsiPresolve *pinfo = new OsiPresolve();
      int presolveActions = 0;
      // Allow dual stuff on integers
      // Allow stuff which may not unroll cleanly
      if ((tuning & 4096) == 0)
        presolveActions = 1 + 16;
      else
        presolveActions = 16; // actually just switch off duplicate columns for ints
      if ((tuning & 32) != 0)
	presolveActions |= 32;
      // Do not allow all +1 to be tampered with
      //if (allPlusOnes)
      //presolveActions |= 2;
      // allow transfer of costs
      presolveActions |= 4;
      // If trying for SOS don't allow some transfers
      if (makeEquality == 2 || makeEquality == 3)
        presolveActions |= 8;
      pinfo->setPresolveActions(presolveActions);
      if (prohibited_)
        assert(numberProhibited_ == oldModel->getNumCols());
      /*
  VIRTUOUS but possible bad for performance 
  
  At this point, the solution is most likely stale: we may have added cuts as
  we left the previous call to modified(), or we may have changed row bounds
  in VUB analysis just above. Continuous presolve doesn't need a solution
  unless we want it to transform the current solution to match the presolved
  model.
*/
      int saveLogLevel = oldModel->messageHandler()->logLevel();
      if (saveLogLevel == 1)
        oldModel->messageHandler()->setLogLevel(0);
      std::string solverName;
      oldModel->getStrParam(OsiSolverName, solverName);
      // Extend if you want other solvers to keep solution
      bool keepSolution = solverName == "clp";
      presolvedModel = pinfo->presolvedModel(*oldModel, feasibilityTolerance, true, 5,
        prohibited_, keepSolution, rowType_,scBound);
      oldModel->messageHandler()->setLogLevel(saveLogLevel);
      if (!presolvedModel) {
        returnModel = NULL;
        delete pinfo;
        break;
      }
      presolvedModel->messageHandler()->setLogLevel(saveLogLevel);
      // update prohibited and rowType
      update(pinfo, presolvedModel,scBound);
      writeDebugMps(presolvedModel, "ordinary2", pinfo);
      model_[iPass] = presolvedModel;
      presolve_[iPass] = pinfo;
      if (!presolvedModel->getNumRows()) {
        // was returnModel=oldModel;
        returnModel = presolvedModel;
        numberSolvers_ = iPass + 1;
        break; // model totally solved
      }
      bool constraints = iPass < numberPasses - 1;
      // Give a hint to do primal
      bool saveTakeHint;
      OsiHintStrength saveStrength;
      presolvedModel->getHintParam(OsiDoDualInInitial,
        saveTakeHint, saveStrength);
      //if (iPass)
      presolvedModel->setHintParam(OsiDoDualInInitial, false, OsiHintTry);
      double inspectLpStart1 = inspect_ ? CoinWallclockTime() : 0.0;
      presolvedModel->initialSolve();
      numberIterationsPre_ += presolvedModel->getIterationCount();
      if (inspect_) {
        FILE *fp = handler_->filePointer();
        if (fp) {
          const char *stat = presolvedModel->isProvenOptimal() ? "optimal"
            : presolvedModel->isProvenPrimalInfeasible() ? "infeasible" : "not optimal";
          fprintf(fp, "  [Preproc LP pass %d] initialSolve: %.2fs  iters=%d  obj=%.6g  (%s)\n",
            iPass - doInitialPresolve + 1,
            CoinWallclockTime() - inspectLpStart1,
            presolvedModel->getIterationCount(),
            presolvedModel->getObjValue(), stat);
          fflush(fp);
        }
      }
      presolvedModel->setHintParam(OsiDoDualInInitial, saveTakeHint, saveStrength);
      if (!presolvedModel->isProvenOptimal()) {
        writeDebugMps(presolvedModel, "bad2", NULL);
        CoinWarmStartBasis *slack = dynamic_cast< CoinWarmStartBasis * >(presolvedModel->getEmptyWarmStart());
        presolvedModel->setWarmStart(slack);
        delete slack;
        presolvedModel->resolve();
        if (!presolvedModel->isProvenOptimal()) {
          returnModel = NULL;
          //printf("infeasible\n");
          break;
        } else {
          //printf("feasible on second try\n");
        }
      }
      // maybe we can fix some
      int numberFixed = reducedCostFix(*presolvedModel);
#if CBC_USEFUL_PRINTING > 1
      if (numberFixed)
        printf("%d variables fixed on reduced cost\n", numberFixed);
#endif
#if DEBUG_PREPROCESS > 1
      const OsiRowCutDebugger *debugger = presolvedModel->getRowCutDebugger();
      if (debugger)
        printf("Contains optimal before modified\n");
#endif
      OsiSolverInterface *newModel = modified(presolvedModel, constraints, numberChanges, iPass - doInitialPresolve, numberModifiedPasses);
#if DEBUG_PREPROCESS > 1
      if (debugger)
        assert(newModel->getRowCutDebugger());
#endif
      returnModel = newModel;
      if (!newModel) {
        break;
      }
      modifiedModel_[iPass] = newModel;
      oldModel = newModel;
      writeDebugMps(newModel, "ordinary3", NULL);
      // tighten bounds
 #if DEBUG_PREPROCESS > 1
      const OsiRowCutDebugger *debugger = newModel->getRowCutDebugger();
      if (debugger)
 	printf("Contains optimal before tightenA\n");
 #endif
      if (!numberChanges && !numberFixed) {
	int change = tightenPrimalBounds(*newModel,true,scBound);
	if (change > 0) {
	  if ((change&0x40000000)!=0) {
	    /** some free rows found!
		always go round again as more likely to have
		bugs if free rows in model */
	    change -= 0x40000000;
	    numberChanges+=change;
	    if (iPass==numberSolvers_-1) {
	      OsiSolverInterface ** modelOld = model_;
	      OsiSolverInterface ** modifiedModelOld = modifiedModel_;
	      OsiPresolve ** presolveOld = presolve_;
	      model_ = new OsiSolverInterface *[numberSolvers_+1];
	      modifiedModel_ = new OsiSolverInterface *[numberSolvers_+1];
	      presolve_ = new OsiPresolve *[numberSolvers_+1];
	      for (int i = 0; i < numberSolvers_; i++) {
		model_[i] = modelOld[i];
		modifiedModel_[i] = modifiedModelOld[i];
		presolve_[i] = presolveOld[i];
	      }
	      delete [] modelOld;
	      delete [] modifiedModelOld;
	      delete [] presolveOld;
	      model_[numberSolvers_] = NULL;
	      modifiedModel_[numberSolvers_] = NULL;
	      presolve_[numberSolvers_] = NULL;
	      numberSolvers_++;
	    }
	    // round again
	  }
	  numberChanges += change;
	} else if (change<0) {
	  // infeasible
	  newModel = NULL;
	  returnModel=NULL;
	  break;
	}
#ifdef CBC_HAS_CLP
	OsiClpSolverInterface * clpSolver = getClpSolver(newModel);
	if (clpSolver) {
	  change = clpSolver->tightenBounds();
	  if (change > 0) {
	    numberChanges += change;
	    //printf("Clp tightenBounds tightened %d\n",change);
	  }
	}
#endif
      }
 #if DEBUG_PREPROCESS > 1
       if (debugger)
         assert(newModel->getRowCutDebugger());
#endif
       {
	 bool saveTakeHint;
	 OsiHintStrength saveStrength;
	 newModel->getHintParam(OsiDoDualInResolve,
				saveTakeHint, saveStrength);
	 newModel->setHintParam(OsiDoDualInResolve, true, OsiHintTry);
	 double inspectLpStart2 = inspect_ ? CoinWallclockTime() : 0.0;
	 newModel->resolve();
	 if (inspect_) {
	   FILE *fp = handler_->filePointer();
	   if (fp) {
	     const char *stat = newModel->isProvenOptimal() ? "optimal"
	       : newModel->isProvenPrimalInfeasible() ? "infeasible" : "not optimal";
	     fprintf(fp, "  [Preproc LP pass %d] resolve after cuts: %.2fs  iters=%d  obj=%.6g  (%s)\n",
	       iPass - doInitialPresolve + 1,
	       CoinWallclockTime() - inspectLpStart2,
	       newModel->getIterationCount(),
	       newModel->getObjValue(), stat);
	     fflush(fp);
	   }
	 }
	 newModel->setHintParam(OsiDoDualInResolve, saveTakeHint, saveStrength);
       }
      if (!newModel->isProvenOptimal()) {
	numberSolvers_ = iPass + 1;
	break;
      }
      if ((!numberChanges && !numberFixed) ||
	  (makeEquality==-2 && iPass == numberSolvers_-1)) {
	if (makeEquality != -2) {
#if CBC_USEFUL_PRINTING > 1
	  printf("exiting after pass %d of %d\n", iPass, numberSolvers_);
#endif
	  numberSolvers_ = iPass + 1;
	  break;
	} else {
	  bool goRoundAgain = false;
	  makeEquality = 0;
	  const CoinPackedMatrix * matrixByColumn = newModel->getMatrixByCol();
	  const int * row = matrixByColumn->getIndices();
	  const CoinBigIndex * columnStart = matrixByColumn->getVectorStarts();
	  const int * columnLength = newModel->getMatrixByCol()->getVectorLengths();
	  const double * columnElements = matrixByColumn->getElements();
	  const double * rowLower = newModel->getRowLower();
	  const double * rowUpper = newModel->getRowUpper();
	  const double * columnLower = newModel->getColLower();
	  const double * columnUpper = newModel->getColUpper();
	  const double * cost = newModel->getObjCoefficients();
	  // Row copy
	  const CoinPackedMatrix * matrixByRow = newModel->getMatrixByRow();
	  const double * elementByRow = matrixByRow->getElements();
	  const int * column = matrixByRow->getIndices();
	  const CoinBigIndex * rowStart = matrixByRow->getVectorStarts();
	  const int * rowLength = matrixByRow->getVectorLengths();
	  int numberRows = newModel->getNumRows();
	  int numberColumns = newModel->getNumCols();
	  int nMadeContinuous = 0;
	  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	    if (columnLength[iColumn]==1 && !cost[iColumn]
		&& newModel->isInteger(iColumn)) {
	      CoinBigIndex start = columnStart[iColumn];
	      int iRow = row[start];
	      if (rowLower[iRow]==rowUpper[iRow] &&
		  fabs(columnElements[start]) == 1.0) {
		bool canDo = false;
		// are all integer
		bool allInteger=true;
		double element = columnElements[columnStart[iColumn]];
		if ((element-floor(element)))
		  allInteger=false;
		for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
		  int jColumn = column[j];
		  if (!newModel->isInteger(jColumn)) {
		    allInteger=false;
		    break;
		  }
		}
		if (allInteger) {
		  nMadeContinuous++;
		  newModel->setContinuous(iColumn);
		}
	      }
	    }
	  }
	  goRoundAgain = nMadeContinuous > 0;
#ifdef CBC_KEEP_OLD_MESSAGES
	  if (nMadeContinuous) {
	    char generalPrint[100];
	    sprintf(generalPrint,"%d integer variables made continuous",nMadeContinuous);
	    handler_->message(CGL_GENERAL, messages_)
	      << generalPrint
	      << CoinMessageEol;
	  }
#endif
	  if (!goRoundAgain) {
	    numberSolvers_ = iPass + 1;
	    break;
	  } else if (iPass==numberSolvers_-1) {
	    OsiSolverInterface ** modelOld = model_;
	    OsiSolverInterface ** modifiedModelOld = modifiedModel_;
	    OsiPresolve ** presolveOld = presolve_;
	    model_ = new OsiSolverInterface *[numberSolvers_+1];
	    modifiedModel_ = new OsiSolverInterface *[numberSolvers_+1];
	    presolve_ = new OsiPresolve *[numberSolvers_+1];
	    for (int i = 0; i < numberSolvers_; i++) {
	      model_[i] = modelOld[i];
	      modifiedModel_[i] = modifiedModelOld[i];
	      presolve_[i] = presolveOld[i];
	    }
	    delete [] modelOld;
	    delete [] modifiedModelOld;
	    delete [] presolveOld;
	    model_[numberSolvers_] = NULL;
	    modifiedModel_[numberSolvers_] = NULL;
	    presolve_[numberSolvers_] = NULL;
	    numberSolvers_++;
	  }
	}
      }
    }
  }
  if (returnModel) {
#if 0
    if (returnModel->getNumRows()) {
      // tighten bounds
#if DEBUG_PREPROCESS > 1
      const OsiRowCutDebugger *debugger = returnModel->getRowCutDebugger();
      if (debugger)
	printf("Contains optimal before tighten\n");
#endif
      int infeas = tightenPrimalBounds(*returnModel,true,scBound);
      infeas = (infeas <0) ? -infeas : 0;
#if DEBUG_PREPROCESS > 1
      if (debugger)
        assert(returnModel->getRowCutDebugger());
      writeDebugMps(returnModel, "afterTighten", NULL);
#endif
      if (infeas) {
	handler_->message(CGL_INFEASIBLE, messages_)
	  << CoinMessageEol;
        delete returnModel;
        for (int iPass = 0; iPass < numberSolvers_; iPass++) {
          if (returnModel == modifiedModel_[iPass])
            modifiedModel_[iPass] = NULL;
        }
        //printf("startModel_ %p startModel2 %p originalModel_ %p returnModel %p\n",
        //     startModel_,startModel2,originalModel_,returnModel);
        if (returnModel == startModel_ && startModel_ != originalModel_)
          startModel_ = NULL;
        returnModel = NULL;
      }
    }
#endif
  } else {
    double timeDone = getCurrentCPUTime()-ppstart;
    if (timeDone < timeLimit_) {
      handler_->message(CGL_INFEASIBLE, messages_)
	<< CoinMessageEol;
    } else {
      char generalPrint[100];
      sprintf(generalPrint,
	      "Preprocessing exiting on time after %g seconds - ignore infeasibility message",timeDone);
      handler_->message(CGL_GENERAL, messages_)
	<< generalPrint
	<< CoinMessageEol;
    }
  }
  int numberIntegers = 0;
  if (returnModel) {
    int iColumn;
    int numberColumns = returnModel->getNumCols();
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (returnModel->isInteger(iColumn))
        numberIntegers++;
    }
  }
  if ((makeEquality == 2 || makeEquality == 3) && numberCliques && returnModel) {
    int iRow, iColumn;
    int numberColumns = returnModel->getNumCols();
    int numberRows = returnModel->getNumRows();
    const double *objective = returnModel->getObjCoefficients();
    // get row copy
    const CoinPackedMatrix *matrix = returnModel->getMatrixByRow();
    const double *element = matrix->getElements();
    const int *column = matrix->getIndices();
    const CoinBigIndex *rowStart = matrix->getVectorStarts();
    const int *rowLength = matrix->getVectorLengths();
    const double *rowLower = returnModel->getRowLower();
    const double *rowUpper = returnModel->getRowUpper();
    const double *columnLower = returnModel->getColLower();

    // Look for possible SOS
    int numberSOS = 0;
    int *mark = new int[numberColumns];
    int *sosRow = new int[numberRows];
    CoinZeroN(sosRow, numberRows);
    CoinFillN(mark, numberColumns, -1);
    int numberOverlap = 0;
    int numberInSOS = 0;
    for (iRow = 0; iRow < numberRows; iRow++) {
      if (rowLower[iRow] == 1.0 && rowUpper[iRow] == 1.0) {
        if ((rowLength[iRow] < 5 && !justOnesWithObj) || (rowLength[iRow] < 20 && allToGub))
          continue;
        bool goodRow = true;
        int nObj = 0;
        for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
          iColumn = column[j];
          if (element[j] != 1.0 || !returnModel->isInteger(iColumn) || columnLower[iColumn]) {
            goodRow = false;
            break;
          }
          if (mark[iColumn] >= 0 && !allToGub) {
            goodRow = false;
            numberOverlap++;
          }
          if (objective[iColumn])
            nObj++;
        }
        if (goodRow && justOnesWithObj) {
          if (!nObj || nObj < rowLength[iRow] - 1)
            goodRow = false;
        }
        if (goodRow) {
          // mark all
          for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
            int iColumn = column[j];
            mark[iColumn] = numberSOS;
          }
          sosRow[numberSOS++] = iRow;
          numberInSOS += rowLength[iRow];
        }
      }
    }
    if (numberSOS) {
      if (makeEquality == 2 && (numberOverlap || numberIntegers > numberInSOS + 1)) {
        handler_->message(CGL_PROCESS_SOS2, messages_)
          << numberSOS << numberInSOS << numberIntegers << numberOverlap
          << CoinMessageEol;
      } else {
        handler_->message(CGL_PROCESS_SOS1, messages_)
          << numberSOS << numberInSOS
          << CoinMessageEol;
        numberSOS_ = numberSOS;
        typeSOS_ = new int[numberSOS_];
        startSOS_ = new int[numberSOS_ + 1];
        whichSOS_ = new int[numberInSOS];
        weightSOS_ = new double[numberInSOS];
        numberInSOS = 0;
        startSOS_[0] = 0;
        const CoinPackedMatrix *columnCopy = returnModel->getMatrixByCol();
        const int *row = columnCopy->getIndices();
        const CoinBigIndex *columnStart = columnCopy->getVectorStarts();
        const int *columnLength = columnCopy->getVectorLengths();
        const double *columnElements = columnCopy->getElements();
        const double *objective = returnModel->getObjCoefficients();
        int *numberInRow = new int[numberRows];
        double *sort = new double[numberColumns];
        int *which = new int[numberColumns];
        for (int iSOS = 0; iSOS < numberSOS_; iSOS++) {
          int n = 0;
          int numberObj = 0;
          CoinZeroN(numberInRow, numberRows);
          int kRow = sosRow[iSOS];
          for (CoinBigIndex j = rowStart[kRow]; j < rowStart[kRow] + rowLength[kRow]; j++) {
            int iColumn = column[j];
            whichSOS_[numberInSOS] = iColumn;
            weightSOS_[numberInSOS] = n;
            numberInSOS++;
            n++;
            if (objective[iColumn])
              numberObj++;
            for (CoinBigIndex j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              if (!sosRow[iRow])
                numberInRow[iRow]++;
            }
          }
          // See if any rows look good
          int bestRow = -1;
          int numberDifferent = 1;
          int start = startSOS_[iSOS];
          for (int iRow = 0; iRow < numberRows; iRow++) {
            if (numberInRow[iRow] >= n - 1) {
              // See how many different
              int i;
              for (i = 0; i < n; i++) {
                int iColumn = whichSOS_[i + start];
                sort[i] = 0.0;
                which[i] = iColumn;
                for (CoinBigIndex j = columnStart[iColumn];
                     j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                  int jRow = row[j];
                  if (jRow == iRow) {
                    sort[i] = columnElements[j];
                    break;
                  }
                }
              }
              // sort
              CoinSort_2(sort, sort + n, which);
              double last = sort[0];
              int nDiff = 1;
              for (i = 1; i < n; i++) {
                if (sort[i] > last + std::max(fabs(last) * 1.0e-8, 1.0e-5)) {
                  nDiff++;
                }
                last = sort[i];
              }
              if (nDiff > numberDifferent) {
                numberDifferent = nDiff;
                bestRow = iRow;
              }
            }
          }
          if (numberObj >= n - 1 || bestRow < 0) {
            int i;
            for (i = 0; i < n; i++) {
              int iColumn = whichSOS_[i + start];
              sort[i] = objective[iColumn];
              which[i] = iColumn;
            }
            // sort
            CoinSort_2(sort, sort + n, which);
            double last = sort[0];
            int nDiff = 1;
            for (i = 1; i < n; i++) {
              if (sort[i] > last + std::max(fabs(last) * 1.0e-8, 1.0e-5)) {
                nDiff++;
              }
              last = sort[i];
            }
            if (nDiff > numberDifferent) {
              numberDifferent = nDiff;
              bestRow = numberRows;
            }
          }
          if (bestRow >= 0) {
            // if not objective - recreate
            if (bestRow < numberRows) {
              int i;
              for (i = 0; i < n; i++) {
                int iColumn = whichSOS_[i + start];
                sort[i] = 0.0;
                which[i] = iColumn;
                for (CoinBigIndex j = columnStart[iColumn];
                     j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                  int jRow = row[j];
                  if (jRow == bestRow) {
                    sort[i] = columnElements[j];
                    break;
                  }
                }
              }
              // sort
              CoinSort_2(sort, sort + n, which);
            }
            // make sure gaps OK
            double last = sort[0];
            for (int i = 1; i < n; i++) {
              double next = last + std::max(fabs(last) * 1.0e-8, 1.0e-5);
              sort[i] = std::max(sort[i], next);
              last = sort[i];
            }
            //CoinCopyN(sort,n,weightSOS_+start);
            //CoinCopyN(which,n,whichSOS_+start);
          }
          typeSOS_[iSOS] = 1;
          startSOS_[iSOS + 1] = numberInSOS;
        }
        delete[] numberInRow;
        delete[] sort;
        delete[] which;
      }
    }
    delete[] mark;
    delete[] sosRow;
  }
#if DEBUG_PREPROCESS > 1
  if (debugger)
    assert(returnModel->getRowCutDebugger());
#endif
  if (returnModel) {
    if (makeIntegers)
      makeIntegers2(returnModel, makeIntegers);
#if DEBUG_PREPROCESS > 1
    if (debugger)
      assert(returnModel->getRowCutDebugger());
#endif
#ifdef SOS_SUB_STUFF
    // start of sos sub stuff
    if (returnModel && !rowType_ && (tuning&6) !=0) {
      int numberColumns = returnModel->getNumCols();
      int numberRows = returnModel->getNumRows();
      const CoinPackedMatrix * rowCopy = returnModel->getMatrixByRow();
      numberRows = rowCopy->getNumRows();
      const int * column = rowCopy->getIndices();
      const int * rowLength = rowCopy->getVectorLengths();
      const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
      const double * rowLower = returnModel->getRowLower();
      const double * rowUpper = returnModel->getRowUpper();
      const double * element = rowCopy->getElements();
      int *sort = new int[2*numberRows];
      int *delrows = sort+numberRows;
      int numberLook = 0;
      int k = numberRows;
      for (int i = 0; i < numberRows; i++) {
	if (rowLength[i] == 0)
	  continue;
	bool possible = true;
	if (rowLower[i]<rowUpper[i])
	  possible = false;
	for (CoinBigIndex j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	  if (element[j]!=1.0) {
	    possible = false;
	    break;
	  }
	}
	if (possible) 
	  sort[numberLook++] = i;
	else
	  sort[--k] = i;
      }
      int numberTotal = numberLook;
      if (numberLook) {
	// move others
	for (int i=k;i<numberRows;i++)
	  sort[numberTotal++] = sort[i];
	
	double *workrow = new double[numberTotal + 1];
	
	double *workcol = new double[numberColumns + 1];
	coin_init_random_vec(workcol, numberColumns);
	for (int i=0;i<numberTotal;i++) {
	  double sum = 0.0;
	  int iRow = sort[i];
	  for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
	    int iColumn = column[j];
	    sum += workcol[iColumn];
	  }
	  workrow[i] = sum;
	}
	CoinSort_2(workrow, workrow + numberLook, sort);
	CoinSort_2(workrow+numberLook, workrow + numberTotal, sort+numberLook);
	
	int numberOut = 0;
	
	int iLook = numberLook;
	for (int kk = 0; kk < numberLook; kk++) {
	  double dval = workrow[kk];
	  int ilast = sort[kk];
	  double req = rowLower[ilast];
	  for (int jj = iLook; jj < numberTotal; jj++) {
	    if (workrow[jj] == dval) {
	      int ithis = sort[jj];
	      CoinBigIndex krs = rowStart[ithis];
	      CoinBigIndex kre = krs + rowLength[ithis];
	      if (rowLength[ithis] == rowLength[ilast]) {
		CoinBigIndex ishift = rowStart[ilast] - krs;
		CoinBigIndex k;
		for (k = krs; k < kre; k++) {
		  if (column[k] != column[k + ishift]) {
		    break;
		  }
		}
		if (k == kre) {
		  /* now check to see if row satisfied */
		  double rlo2 = rowLower[ithis];
		  double rup2 = rowUpper[ithis];
		  for (k = krs; k < kre; k++) {
		    double value = req*element[k];
		    if (value<rlo2-feasibilityTolerance ||
			value>rup2+feasibilityTolerance)
		      break;
		  }
		}
		if (k == kre) {
		  delrows[numberOut++] = ithis;
		}
	      }
	    } else if (workrow[jj] < dval) {
	      iLook++;
	    } else {
	      break;
	    }
	  }
	}
	
	delete[] workrow;
	delete[] workcol;
	
	if (numberOut) {
#ifdef CBC_KEEP_OLD_MESSAGES
	  char generalPrint[100];
	  sprintf(generalPrint, "%d more redundant rows found",
		  numberOut);
	  handler_->message(CGL_GENERAL, messages_)
	    << generalPrint
	    << CoinMessageEol;
#endif
	  returnModel->deleteRows(numberOut,delrows);
	  returnModel->resolve();
	}
      }
      delete[] sort;
    }
#endif // end of sos sub stuff
#if DEBUG_PREPROCESS > 1
    if (debugger)
      assert(returnModel->getRowCutDebugger());
#endif
    // If can make some cuts then do so
    if (rowType_) {
      int numberRows = returnModel->getNumRows();
      int numberCuts = 0;
      for (int i = 0; i < numberRows; i++) {
        if (rowType_[i] > 0)
          numberCuts++;
      }
      if (numberCuts) {
        CglStored stored;

        int *whichRow = new int[numberRows];
        // get row copy
        const CoinPackedMatrix *rowCopy = returnModel->getMatrixByRow();
        const int *column = rowCopy->getIndices();
        const int *rowLength = rowCopy->getVectorLengths();
        const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
        const double *rowLower = returnModel->getRowLower();
        const double *rowUpper = returnModel->getRowUpper();
        const double *element = rowCopy->getElements();
        int iRow, nDelete = 0;
        for (iRow = 0; iRow < numberRows; iRow++) {
          if (rowType_[iRow] == 1) {
            // take out
            whichRow[nDelete++] = iRow;
          }
        }
        for (int jRow = 0; jRow < nDelete; jRow++) {
          iRow = whichRow[jRow];
          CoinBigIndex start = rowStart[iRow];
          stored.addCut(rowLower[iRow], rowUpper[iRow], rowLength[iRow],
            column + start, element + start);
        }
        returnModel->deleteRows(nDelete, whichRow);
        delete[] whichRow;
        cuts_ = stored;
      }
    }
  }
#if DEBUG_PREPROCESS > 1
  if (debugger)
    assert(returnModel->getRowCutDebugger());
#endif
#if 0
  if (returnModel) {
    int numberColumns = returnModel->getNumCols();
    int numberRows = returnModel->getNumRows();
    int * del = new int [std::max(numberColumns,numberRows)];
    int * original = new int [numberColumns];
    int nDel=0;
    for (int i=0;i<numberColumns;i++) {
      original[i]=i;
      if (returnModel->isInteger(i))
	del[nDel++]=i;
    }
    int nExtra=0;
    if (nDel&&nDel!=numberColumns&&(options_&1)!=0&&false) {
      OsiSolverInterface * yyyy = returnModel->clone();
      int nPass=0;
      while (nDel&&nPass<10) {
	nPass++;
	OsiSolverInterface * xxxx = yyyy->clone();
	int nLeft=0;
	for (int i=0;i<nDel;i++) 
	  original[del[i]]=-1;
	for (int i=0;i<numberColumns;i++) {
	  int kOrig=original[i];
	  if (kOrig>=0)
	    original[nLeft++]=kOrig;
	}
	assert (nLeft==numberColumns-nDel);
	xxxx->deleteCols(nDel,del);
	numberColumns = xxxx->getNumCols();
	const CoinPackedMatrix * rowCopy = xxxx->getMatrixByRow();
	numberRows = rowCopy->getNumRows();
	const int * column = rowCopy->getIndices();
	const int * rowLength = rowCopy->getVectorLengths();
	const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
	const double * rowLower = xxxx->getRowLower();
	const double * rowUpper = xxxx->getRowUpper();
	const double * element = rowCopy->getElements();
        const CoinPackedMatrix * columnCopy = xxxx->getMatrixByCol();
        const int * columnLength = columnCopy->getVectorLengths(); 
	nDel=0;
	// Could do gcd stuff on ones with costs
	for (int i=0;i<numberRows;i++) {
	  if (!rowLength[i]) {
	    del[nDel++]=i;
	  } else if (rowLength[i]==1) {
	    int k=rowStart[i];
	    int iColumn = column[k];
	    if (!xxxx->isInteger(iColumn)) {
	      double mult =1.0/fabs(element[k]);
	      if (rowLower[i]<-1.0e20) {
		double value = rowUpper[i]*mult;
		if (fabs(value-floor(value+0.5))<1.0e-8) {
		  del[nDel++]=i;
		  if (columnLength[iColumn]==1) {
		    xxxx->setInteger(iColumn);
		    int kOrig=original[iColumn];
		    returnModel->setInteger(kOrig);
		  }
		}
	      } else if (rowUpper[i]>1.0e20) {
		double value = rowLower[i]*mult;
		if (fabs(value-floor(value+0.5))<1.0e-8) {
		  del[nDel++]=i;
		  if (columnLength[iColumn]==1) {
		    xxxx->setInteger(iColumn);
		    int kOrig=original[iColumn];
		    returnModel->setInteger(kOrig);
		  }
		}
	      } else {
		double value = rowUpper[i]*mult;
		if (rowLower[i]==rowUpper[i]&&
		    fabs(value-floor(value+0.5))<1.0e-8) {
		  del[nDel++]=i;
		  xxxx->setInteger(iColumn);
		  int kOrig=original[iColumn];
		  returnModel->setInteger(kOrig);
		}
	      }
	    }
	  } else {
	    // only if all singletons
	    bool possible=false;
	    if (rowLower[i]<-1.0e20) {
	      double value = rowUpper[i];
	      if (fabs(value-floor(value+0.5))<1.0e-8) 
		possible=true;
	    } else if (rowUpper[i]>1.0e20) {
	      double value = rowLower[i];
	      if (fabs(value-floor(value+0.5))<1.0e-8) 
		possible=true;
	    } else {
	      double value = rowUpper[i];
	      if (rowLower[i]==rowUpper[i]&&
		  fabs(value-floor(value+0.5))<1.0e-8)
		possible=true;
	    }
	    if (possible) {
	      for (CoinBigIndex j=rowStart[i];
		   j<rowStart[i]+rowLength[i];j++) {
		int iColumn = column[j];
		if (columnLength[iColumn]!=1||fabs(element[j])!=1.0) {
		  possible=false;
		  break;
		}
	      }
	      if (possible) {
		for (CoinBigIndex j=rowStart[i];
		     j<rowStart[i]+rowLength[i];j++) {
		  int iColumn = column[j];
		  if (!xxxx->isInteger(iColumn)) {
		    xxxx->setInteger(iColumn);
		    int kOrig=original[iColumn];
		    returnModel->setInteger(kOrig);
		  }
		}
		del[nDel++]=i;
	      }
	    }
	  }
	}
	if (nDel) {
	  xxxx->deleteRows(nDel,del);
	}
	if (nDel!=numberRows) {
	  nDel=0;
	  for (int i=0;i<numberColumns;i++) {
	    if (xxxx->isInteger(i)) {
	      del[nDel++]=i;
	      nExtra++;
	    }
	  }
	} 
	delete yyyy;
	yyyy=xxxx->clone();
      }
      numberColumns = yyyy->getNumCols();
      numberRows = yyyy->getNumRows();
      if (!numberColumns||!numberRows) {
	printf("All gone\n");
	int numberColumns = returnModel->getNumCols();
	for (int i=0;i<numberColumns;i++)
	  assert(returnModel->isInteger(i));
      }
      // Would need to check if original bounds integer
      //yyyy->writeMps("noints");
      delete yyyy;
      printf("Creating simplified model with %d rows and %d columns - %d extra integers\n",
	     numberRows,numberColumns,nExtra);
    }
    delete [] del;
    delete [] original;
    //exit(2);
  }
#endif
  //writeDebugMps(returnModel, "returnModel", NULL);
#if DEBUG_PREPROCESS > 1
  if (debugger)
    assert(returnModel->getRowCutDebugger());
#endif


  if (returnModel && returnModel != &model && keepColumnNames_)
  {
    returnModel->setIntParam( OsiNameDiscipline, 1 );
    for ( int i=0 ; (i<returnModel->getNumCols()) ; i++ ) {
      int iColumn = originalColumns()[i];
      if (iColumn<modelIn->getNumCols())
	returnModel->setColName( i, modelIn->getColName(iColumn) );
    }
  }
  // clean model
  if (returnModel) {
    int numberRows = returnModel->getNumRows();
    int numberColumns = returnModel->getNumCols();
    CoinPackedMatrix matrixByRow(*returnModel->getMutableMatrixByRow());
    // clean rhs
    double *elementByRow = matrixByRow.getMutableElements();
    const int *column = matrixByRow.getIndices();
    const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
    const int *rowLength = matrixByRow.getVectorLengths();
    
    const double *rowLower = returnModel->getRowLower();
    const double *rowUpper = returnModel->getRowUpper();
    for (int iRow=0;iRow<numberRows;iRow++) {
      CoinBigIndex start = rowStart[iRow];
      CoinBigIndex end = start+rowLength[iRow];
      bool allInteger = true;
      for (CoinBigIndex j=start;j<end;j++) {
	if (!returnModel->isInteger(column[j])) {
	  allInteger = false;
	  break;
	}
	double value = fabs(elementByRow[j]);
	if (value-floor(value+0.5)>1.0e-13) {
	  allInteger = false;
	  break;
	}
      }
      if (allInteger) {
	double lower = rowLower[iRow];
	double upper = rowUpper[iRow];
	double lowerNew = lower;
	double upperNew = upper;
	if (lower>-1.0e15&&lower!=floor(lower+0.5)) {
	  if (fabs(lower-floor(lower+0.5))<1.0e-13) {
	    lowerNew = floor(lower+0.5);
	  }
	}
	if (upper<1.0e15&&upper!=floor(upper+0.5)) {
	  if (fabs(upper-floor(upper+0.5))<1.0e-13) {
	    upperNew = floor(upper+0.5);
	  }
	}
	if (lowerNew<=upperNew) {
	  if (lowerNew!=lower)
	    returnModel->setRowLower(iRow,lowerNew);
	  if (upperNew!=upper)
	    returnModel->setRowUpper(iRow,upperNew);
	} else {
	  printf("LOOKS infeas\n");
	}
      }
    }
  }
  // end clean model
  delete [] scBound;
#if CBC_USE_PAPILO
  if (doPapiloPresolve) {
    // put back
    clpSolverIn->swapModelPtr(initialTry.originalModel);
  }
  // may be better to do at end ??
  if (clpSolverIn && (keepPapilo.presolveType&2) !=0 && returnModel) {
    double ppstart2 = getCurrentCPUTime();
    OsiClpSolverInterface * clpSolverOut =
      getClpSolver(returnModel);
    ClpSimplex * simplexOut = clpSolverOut->getModelPtr();
    initialTry = papiloPresolve(simplexOut,false,timeLimit_-ppstart2);
    initialTry.presolveType &= ~2;
    initialTry.presolveType |= 8;
    initialTry.beforePapiloEnd = clpSolverOut;
    char generalPrint[100];
    sprintf(generalPrint,
	    "papilo presolve took %g seconds - total %g seconds",
	    getCurrentCPUTime()-ppstart2,getCurrentCPUTime()-ppstart);
    handler_->message(CGL_GENERAL, messages_)
      << generalPrint
      << CoinMessageEol;
    if (initialTry.returnCode<0 && getCurrentCPUTime()<timeLimit_) {
      // infeasible
      sprintf(generalPrint,"papilo says infeasible");
      handler_->message(CGL_GENERAL, messages_)
	<< generalPrint
	<< CoinMessageEol;
      return NULL;
    } else if (initialTry.returnCode==0) {
      memset(&initialTry,0,sizeof(papiloStruct));// no change
      handler_->message(CGL_GENERAL, messages_)
	<< "papilo did not reduce size of problem"
	<< CoinMessageEol;
    } else {
      initialTry.preProcess = this;
      initialTry.presolveType |= 128; // say type 2
      // swap
      ClpSimplex * presolved = initialTry.presolvedModel;
      int numberColumns2 = presolved->numberColumns();
      int nRowsWon = simplexOut->numberRows()-presolved->numberRows();
      int nColsWon = simplexOut->numberColumns()-numberColumns2;
      char rwon[10],cwon[10];
      if (nRowsWon>0)
	sprintf(rwon,"(-%d)",nRowsWon);
      else if (nRowsWon<0)
	sprintf(rwon,"(+%d!!)",-nRowsWon);
      else
	rwon[0]='\0';
      if (nColsWon>0)
	sprintf(cwon,"(-%d)",nColsWon);
      else if (nColsWon<0)
	sprintf(cwon,"(+%d!!)",-nColsWon);
      else
	cwon[0]='\0';
      bool onlyBounds = nRowsWon==0&&nColsWon==0;
      sprintf(generalPrint,"After papilo - problem has %d%s rows, %d%s columns %s",
	      presolved->numberRows(),rwon,
	      numberColumns2,cwon,onlyBounds ? "only bounds tightened" : "");
      handler_->message(CGL_GENERAL, messages_)
	<< generalPrint
	<< CoinMessageEol;
      assert (numberColumns2<=clpSolverOut->getNumCols());
      if (clpSolverOut->integerInformation()) {
	for (int i=0;i<numberColumns2;i++) {
	  if (presolved->isInteger(i))
	    clpSolverOut->setIntegerType(i,1);
	  else
	    clpSolverOut->setIntegerType(i,0);
	}
      }
      initialTry.originalModel = clpSolverOut->swapModelPtr(initialTry.presolvedModel);
      initialTry.beforePapiloEnd = returnModel;
    }
  }
#endif
  if (returnModel) {
    int numberIntegers = 0;
    int numberBinary = 0;
    int numberColumns = returnModel->getNumCols();
    for (int i = 0; i < numberColumns; i++) {
      if (returnModel->isInteger(i)) {
        numberIntegers++;
        if (returnModel->isBinary(i)) {
          numberBinary++;
        }
      }
    }
    handler_->message(CGL_PROCESS_STATS2, messages_)
      << returnModel->getNumRows() << numberColumns
      << numberIntegers << numberBinary << returnModel->getNumElements()
      << CoinMessageEol;
  }
  return returnModel;
}
static int 
tighten(double *colLower, double * colUpper,
	const int *column, const double *rowElements, 
	const CoinBigIndex *rowStart, 
	const CoinBigIndex * rowStartPos,const int * rowLength,
	double *rowLower, double *rowUpper, 
	int nRows,int nCols,char * intVar,int maxpass,
	double tolerance);

/* Tightens primal bounds to make dual and branch and cutfaster.  Unless
   fixed, bounds are slightly looser than they could be.
   Returns negative if problem infeasible, number tightened if feasible
*/
int
CglPreProcess::tightenPrimalBounds(OsiSolverInterface &model,
				   bool tightenRowBounds,
				   double * scBound)
{
  // temporary
  if (scBound)
    return 0;
  // Get a row copy in standard format
  CoinPackedMatrix copy = *model.getMatrixByRow();
  // get matrix data pointers
  int *column = copy.getMutableIndices();
  const CoinBigIndex *rowStart = copy.getVectorStarts();
  const int *rowLength = copy.getVectorLengths();
  double *element = copy.getMutableElements();
  int numberChanged = 1, iPass = 0;
  double large = 1.0e100; // treat bounds > this as infinite
  int numberInfeasible = 0;
  int totalTightened = 0;

  double tolerance;
  model.getDblParam(OsiPrimalTolerance, tolerance);

  int numberColumns = model.getNumCols();
  const double *colLower = model.getColLower();
  const double *colUpper = model.getColUpper();
#ifdef CBC_PREPROCESS_EXPERIMENT
  const double *objective = model.getObjCoefficients();
#endif
  double direction = model.getObjSense();
  const int *columnLength = model.getMatrixByCol()->getVectorLengths();
  // New and saved column bounds
  double *newLower = new double[numberColumns];
  memcpy(newLower, colLower, numberColumns * sizeof(double));
  double *newUpper = new double[numberColumns];
  memcpy(newUpper, colUpper, numberColumns * sizeof(double));
  double *columnLower = new double[numberColumns];
  memcpy(columnLower, colLower, numberColumns * sizeof(double));
  double *columnUpper = new double[numberColumns];
  memcpy(columnUpper, colUpper, numberColumns * sizeof(double));
  int iRow, iColumn;

  int numberRows = model.getNumRows();
  const double *rowLower = model.getRowLower();
  const double *rowUpper = model.getRowUpper();
  int nFreed = 0;
#if 1
  char * intVar = new char [numberColumns];
  for (int i=0;i<numberColumns;i++) {
    if (model.isInteger(i))
      intVar[i]=1;
    else // wasif (!prohibited_[i])
      intVar[i]=0;
    //else
    //intVar[i]=-1;
  }
  // TEMP COPIES
  {
    double * cLower = new double [2*numberRows+2*numberColumns];
    double * cUpper = cLower+numberColumns;
    double * rLower = cUpper+numberColumns;
    double * rUpper = rLower+numberRows;
    memcpy(cLower, colLower, numberColumns * sizeof(double));
    memcpy(cUpper, colUpper, numberColumns * sizeof(double));
    memcpy(rLower, rowLower, numberRows * sizeof(double));
    memcpy(rUpper, rowUpper, numberRows * sizeof(double));
    CoinBigIndex * rowStartPos = new CoinBigIndex[numberRows];
    // sort rows and set up rowStartPos
    for (int iRow=0;iRow<numberRows;iRow++) {
      CoinBigIndex rStart = rowStart[iRow];
      CoinBigIndex rEnd = rStart + rowLength[iRow];
      CoinSort_2(element+rStart,element+rEnd,column+rStart);
      CoinBigIndex j;
      for (j=rStart;j<rEnd;j++) {
	if (element[j]>0)
	  break;
      }
      rowStartPos[iRow]=j;
    }
    int nInfeas = tighten(cLower, cUpper, column, element, 
			  rowStart, rowStartPos, rowLength, rLower, rUpper, 
			  numberRows, numberColumns, intVar,20,tolerance);
    // see if any free rows
    for (int iRow=0;iRow<numberRows;iRow++) {
      if (rLower[iRow]<-1.0e20&&rUpper[iRow]>1.0e20) {
	model.setRowLower(iRow,-COIN_DBL_MAX);
	model.setRowUpper(iRow,COIN_DBL_MAX);
      }
    }
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      if (cLower[iColumn]>newLower[iColumn]) {
	newLower[iColumn] = cLower[iColumn];
	columnLower[iColumn] = cLower[iColumn];
	model.setColLower(iColumn,cLower[iColumn]);
      }
      if (cUpper[iColumn]<newUpper[iColumn]) {
	newUpper[iColumn] = cUpper[iColumn];
	columnUpper[iColumn] = cUpper[iColumn];
	model.setColUpper(iColumn,cUpper[iColumn]);
      }
    }
    delete [] rowStartPos;
    delete [] cLower;
  }
#endif
#ifndef NDEBUG
  double large2 = 1.0e10 * large;
#endif
#define MAXPASS 10

  // Loop round seeing if we can tighten bounds
  // Would be faster to have a stack of possible rows
  // and we put altered rows back on stack
  int numberCheck = -1;
  while (numberChanged > numberCheck) {

    numberChanged = 0; // Bounds tightened this pass

    if (iPass == MAXPASS)
      break;
    iPass++;

    for (iRow = 0; iRow < numberRows; iRow++) {

      if (rowLower[iRow] > -large || rowUpper[iRow] < large) {

        // possible row
        int infiniteUpper = 0;
        int infiniteLower = 0;
        double maximumUp = 0.0;
        double maximumDown = 0.0;
        double newBound;
        CoinBigIndex rStart = rowStart[iRow];
        CoinBigIndex rEnd = rowStart[iRow] + rowLength[iRow];
        CoinBigIndex j;
#ifdef CBC_PREPROCESS_EXPERIMENT
	bool allInteger = true;
#endif
        // Compute possible lower and upper ranges

        for (j = rStart; j < rEnd; ++j) {
          double value = element[j];
          iColumn = column[j];
#ifdef CBC_PREPROCESS_EXPERIMENT
	  if (!model.isInteger(iColumn) || value != floor(value+0.5))
	    allInteger = false;
#endif
          if (value > 0.0) {
            if (newUpper[iColumn] >= large) {
              ++infiniteUpper;
            } else {
              maximumUp += newUpper[iColumn] * value;
            }
            if (newLower[iColumn] <= -large) {
              ++infiniteLower;
            } else {
              maximumDown += newLower[iColumn] * value;
            }
          } else if (value < 0.0) {
            if (newUpper[iColumn] >= large) {
              ++infiniteLower;
            } else {
              maximumDown += newUpper[iColumn] * value;
            }
            if (newLower[iColumn] <= -large) {
              ++infiniteUpper;
            } else {
              maximumUp += newLower[iColumn] * value;
            }
          }
        }
        // Build in a margin of error
        maximumUp += 1.0e-8 * fabs(maximumUp);
        maximumDown -= 1.0e-8 * fabs(maximumDown);
        double maxUp = (!infiniteUpper) ? maximumUp : COIN_DBL_MAX;
        double maxDown = (!infiniteLower) ? maximumDown : -COIN_DBL_MAX;
        if (maxUp <= rowUpper[iRow] + tolerance && maxDown >= rowLower[iRow] - tolerance) {

          // Row is redundant - make totally free
	  model.setRowLower(iRow,-COIN_DBL_MAX);
	  model.setRowUpper(iRow,COIN_DBL_MAX);
        } else {
          if (maxUp < rowLower[iRow] - 100.0 * tolerance || maxDown > rowUpper[iRow] + 100.0 * tolerance) {
            // problem is infeasible - exit at once
            numberInfeasible++;
            break;
          }
          double lower = rowLower[iRow];
          double upper = rowUpper[iRow];
#ifdef CBC_PREPROCESS_EXPERIMENT
	  if (tightenRowBounds && upper>lower && allInteger) {
	    double lowerNew = lower;
	    double upperNew = upper;
	    maxUp = std::max(maxUp,lower);
	    maxDown = std::min(maxDown,upper);
	    if (!infiniteLower && maxDown > lower + 1.0e-6) 
	      lowerNew = std::max(maxDown-1.0e-6,lower);
	    if (!infiniteUpper && maxUp < upper - 1.0e-6)
	      upperNew = std::min(maxUp+1.0e-6,upper);
	    if (lowerNew > upperNew) {
	      printf("BAD bounds of %g (%g) and %g (%g) on row %d\n",
		     lowerNew,lower,upperNew,upper,iRow);
	      lowerNew = upperNew;
	    }
	    if (lower!=lowerNew||upper!=upperNew) {
	      if (allInteger) {
#ifdef LOTS_OF_PRINTING
		printf("On row %d bounds -> %g,%g\n",
		       iRow,
		       lower,upper);
#endif
		lower = std::max(lower,ceil(lowerNew-1.0e-1));
		upper = std::min(upper,floor(upperNew+1.0e-1));
	      } else {
		lower = lowerNew;
		upper = upperNew;
	      }
	      if (lower!=rowLower[iRow]||upper!=rowUpper[iRow])
#ifdef LOTS_OF_PRINTING
		printf("On row %d bounds %g,%g -> %g,%g\n",
		       iRow,rowLower[iRow],rowUpper[iRow],
		       lower,upper);
#endif
	      model.setRowLower(iRow,lower);
	      model.setRowUpper(iRow,upper);
	    }
	  }
#endif
          for (j = rStart; j < rEnd; ++j) {
            double value = element[j];
            iColumn = column[j];
            double nowLower = newLower[iColumn];
            double nowUpper = newUpper[iColumn];
	    bool anyChange = false;
            if (value > 0.0) {
              // positive value
              if (lower > -large) {
                if (!infiniteUpper) {
                  assert(nowUpper < large2);
                  newBound = nowUpper + (lower - maximumUp) / value;
                  // relax if original was large
                  if (fabs(maximumUp) > 1.0e8)
                    newBound -= 1.0e-12 * fabs(maximumUp);
                } else if (infiniteUpper == 1 && nowUpper > large) {
                  newBound = (lower - maximumUp) / value;
                  // relax if original was large
                  if (fabs(maximumUp) > 1.0e8)
                    newBound -= 1.0e-12 * fabs(maximumUp);
                } else {
                  newBound = -COIN_DBL_MAX;
                }
                if (newBound > nowLower + 1.0e-12 && newBound > -large) {
                  // Tighten the lower bound
                  newLower[iColumn] = newBound;
		  anyChange = true;
                  // check infeasible (relaxed)
                  if (nowUpper - newBound < -100.0 * tolerance) {
                    numberInfeasible++;
                  }
                  // adjust
                  double now;
                  if (nowLower < -large) {
                    now = 0.0;
                    infiniteLower--;
                  } else {
                    now = nowLower;
                  }
                  maximumDown += (newBound - now) * value;
                  nowLower = newBound;
                }
              }
              if (upper < large) {
                if (!infiniteLower) {
                  assert(nowLower > -large2);
                  newBound = nowLower + (upper - maximumDown) / value;
                  // relax if original was large
                  if (fabs(maximumDown) > 1.0e8)
                    newBound += 1.0e-12 * fabs(maximumDown);
                } else if (infiniteLower == 1 && nowLower < -large) {
                  newBound = (upper - maximumDown) / value;
                  // relax if original was large
                  if (fabs(maximumDown) > 1.0e8)
                    newBound += 1.0e-12 * fabs(maximumDown);
                } else {
                  newBound = COIN_DBL_MAX;
                }
                if (newBound < nowUpper - 1.0e-12 && newBound < large) {
                  // Tighten the upper bound
                  newUpper[iColumn] = newBound;
		  anyChange = true;
                  // check infeasible (relaxed)
                  if (newBound - nowLower < -100.0 * tolerance) {
                    numberInfeasible++;
                  }
                  // adjust
                  double now;
                  if (nowUpper > large) {
                    now = 0.0;
                    infiniteUpper--;
                  } else {
                    now = nowUpper;
                  }
                  maximumUp += (newBound - now) * value;
                  nowUpper = newBound;
                }
              }
            } else {
              // negative value
              if (lower > -large) {
                if (!infiniteUpper) {
                  assert(nowLower < large2);
                  newBound = nowLower + (lower - maximumUp) / value;
                  // relax if original was large
                  if (fabs(maximumUp) > 1.0e8)
                    newBound += 1.0e-12 * fabs(maximumUp);
                } else if (infiniteUpper == 1 && nowLower < -large) {
                  newBound = (lower - maximumUp) / value;
                  // relax if original was large
                  if (fabs(maximumUp) > 1.0e8)
                    newBound += 1.0e-12 * fabs(maximumUp);
                } else {
                  newBound = COIN_DBL_MAX;
                }
                if (newBound < nowUpper - 1.0e-12 && newBound < large) {
                  // Tighten the upper bound
                  newUpper[iColumn] = newBound;
		  anyChange = true;
                  // check infeasible (relaxed)
                  if (newBound - nowLower < -100.0 * tolerance) {
                    numberInfeasible++;
                  }
                  // adjust
                  double now;
                  if (nowUpper > large) {
                    now = 0.0;
                    infiniteLower--;
                  } else {
                    now = nowUpper;
                  }
                  maximumDown += (newBound - now) * value;
                  nowUpper = newBound;
                }
              }
              if (upper < large) {
                if (!infiniteLower) {
                  assert(nowUpper < large2);
                  newBound = nowUpper + (upper - maximumDown) / value;
                  // relax if original was large
                  if (fabs(maximumDown) > 1.0e8)
                    newBound -= 1.0e-12 * fabs(maximumDown);
                } else if (infiniteLower == 1 && nowUpper > large) {
                  newBound = (upper - maximumDown) / value;
                  // relax if original was large
                  if (fabs(maximumDown) > 1.0e8)
                    newBound -= 1.0e-12 * fabs(maximumDown);
                } else {
                  newBound = -COIN_DBL_MAX;
                }
                if (newBound > nowLower + 1.0e-12 && newBound > -large) {
                  // Tighten the lower bound
                  newLower[iColumn] = newBound;
		  anyChange = true;
                  // check infeasible (relaxed)
                  if (nowUpper - newBound < -100.0 * tolerance) {
                    numberInfeasible++;
                  }
                  // adjust
                  double now;
                  if (nowLower < -large) {
                    now = 0.0;
                    infiniteUpper--;
                  } else {
                    now = nowLower;
                  }
                  maximumUp += (newBound - now) * value;
                  nowLower = newBound;
                }
              }
            }
	    if (anyChange) {
	      numberChanged++;
	    } else if (columnLength[iColumn] == 1) {
#if 0 //def CBC_PREPROCESS_EXPERIMENT
	      // may be able to do better
	      // should be picked up elsewhere if no objective
	      if (objective[iColumn]) {
		double newBound;
		if (direction*objective[iColumn]>0.0) {
		  // want objective as low as possible so reduce upper bound
		  if (value < 0.0) {
		    double gap = maxUp-upper;
		    newBound = newLower[iColumn] - gap/value;
		  } else {
		    double gap = lower-maxDown;
		    newBound = newLower[iColumn] + gap/value;
		  }
		  if (newBound>1.0e50)
		    newBound = COIN_DBL_MAX;
		  if (newBound<newUpper[iColumn]-1.0e-7) {
		    if (model.isInteger(iColumn)) {
		      newBound = ceil(newBound);
		      newBound = std::max(newLower[iColumn],newBound);
		    } else {
		      newBound = std::max(newLower[iColumn],newBound+1.0e-7);
		    }
		  }
		  if (newBound<newUpper[iColumn]-1.0e-7) {
		    numberChanged++;
#ifdef LOTS_OF_PRINTING
		    printf("singleton %d obj %g %g <= %g rlo %g rup %g maxd %g maxu %g el %g\n",
			   iColumn,objective[iColumn],newLower[iColumn],newUpper[iColumn],
			   lower,upper,maxDown,maxUp,value);
		    printf("upperbound changed to %g\n",newBound);
#endif
		    newUpper[iColumn] = newBound;
		    anyChange = true;
		    // check infeasible (relaxed)
		    if (newBound - nowLower < -100.0 * tolerance) {
		      numberInfeasible++;
		    }
		    // adjust
		    double now;
		    if (nowUpper > large) {
		      now = 0.0;
		      infiniteUpper--;
		    } else {
		      now = nowUpper;
		    }
		    maximumUp += (newBound - now) * value;
		    maxUp += (newBound - now) * value;
		    nowUpper = newBound;
		  }
		} else {
		  // want objective as low as possible so increase lower bound
		  if (value > 0.0) { // ?
		    double gap = maxUp-upper;
		    newBound = newUpper[iColumn] - gap/value;
		  } else {
		    double gap = lower-maxDown;
		    newBound = newUpper[iColumn] + gap/value;
		  }
		  if (newBound<-1.0e50)
		    newBound = -COIN_DBL_MAX;
		  if (newBound>newLower[iColumn]+1.0e-7) {
		    if (model.isInteger(iColumn)) {
		      newBound = ceil(newBound);
		      newBound = std::min(newUpper[iColumn],newBound);
		    } else {
		      newBound = std::min(newUpper[iColumn],newBound+1.0e-7);
		    }
		  }
		  if (newBound>newLower[iColumn]+1.0e-7) {
		    numberChanged++;
#ifdef LOTS_OF_PRINTING
		    printf("singleton %d obj %g %g <= %g rlo %g rup %g maxd %g maxu %g el %g\n",
			   iColumn,objective[iColumn],newLower[iColumn],newUpper[iColumn],
			   lower,upper,maxDown,maxUp,value);
		    printf("lowerbound changed to %g\n",newBound);
#endif
		    newLower[iColumn] = newBound;
		    anyChange = true;
		    // check infeasible (relaxed)
		    if (nowUpper - newBound < -100.0 * tolerance) {
		      numberInfeasible++;
		    }
		    // adjust
		    double now;
		    if (nowLower < -large) {
		      now = 0.0;
		      infiniteLower--;
		    } else {
		      now = nowLower;
		    }
		    maximumDown += (newBound - now) * value;
		    maxDown += (newBound - now) * value;
		    nowLower = newBound;
		  }
		}
	      }
#endif
	    }
          }
        }
      }
    }
    totalTightened += numberChanged;
    if (iPass == 1)
      numberCheck = numberChanged >> 4;
    if (numberInfeasible)
      break;
  }
  delete [] intVar;
  if (!numberInfeasible) {
    // Set bounds slightly loose unless integral - now tighter
    double useTolerance = 1.0e-5;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (columnUpper[iColumn] > columnLower[iColumn]) {
        double lower = newLower[iColumn];
        double upper = newUpper[iColumn];
        if (model.isInteger(iColumn)) {
          if (fabs(lower - floor(lower + 0.5)) < 1.0e-5)
            lower = floor(lower + 0.5);
          else
            lower = ceil(lower);
          if (fabs(upper - floor(upper + 0.5)) < 1.0e-5)
            upper = floor(upper + 0.5);
          else
            upper = floor(upper);
          if (lower > upper)
            numberInfeasible++;
        } else {
          if (fabs(upper) < 1.0e-8 && fabs(lower) < 1.0e-8) {
            lower = 0.0;
            upper = 0.0;
	  } else if (scBound &&scBound[iColumn]!=-COIN_DBL_MAX) {
	    double lowerSc =scBound[iColumn];
	    if (upper<lowerSc-1.0e-5) {
	      upper = 0.0;
	      scBound[iColumn] = -COIN_DBL_MAX;
	      prohibited_[iColumn]=0;
	    } else if (lower > 1.0e-5) {
	      lower = lowerSc;
	      scBound[iColumn] = -COIN_DBL_MAX;
	      prohibited_[iColumn]=0;
	    }
	    if (lower > upper+1.0e-5)
	      numberInfeasible++;
          } else {
            // Relax unless integral
            if (fabs(lower - floor(lower + 0.5)) > 1.0e-9)
              lower -= useTolerance;
            else
              lower = floor(lower + 0.5);
            lower = std::max(columnLower[iColumn], lower);
            if (fabs(upper - floor(upper + 0.5)) > 1.0e-9)
              upper += useTolerance;
            else
              upper = floor(upper + 0.5);
            upper = std::min(columnUpper[iColumn], upper);
          }
        }
        model.setColLower(iColumn, lower);
        model.setColUpper(iColumn, upper);
        newLower[iColumn] = lower;
        newUpper[iColumn] = upper;
      } else if (columnUpper[iColumn] < columnLower[iColumn]-1.0e-7) {
	numberInfeasible++;
      }
    }
    if (!numberInfeasible) {
      // check common bad formulations
      int numberChanges = 0;
      for (iRow = 0; iRow < numberRows; iRow++) {
        if (rowLower[iRow] > -large || rowUpper[iRow] < large) {
          // possible row
          double sumFixed = 0.0;
          int infiniteUpper = 0;
          int infiniteLower = 0;
          double maximumUp = 0.0;
          double maximumDown = 0.0;
          double largest = 0.0;
          CoinBigIndex rStart = rowStart[iRow];
          CoinBigIndex rEnd = rowStart[iRow] + rowLength[iRow];
          CoinBigIndex j;
          //int numberInteger = 0;
          //int whichInteger=-1;
          // Compute possible lower and upper ranges
          for (j = rStart; j < rEnd; ++j) {
            double value = element[j];
            iColumn = column[j];
            if (newUpper[iColumn] > newLower[iColumn]) {
              if (model.isInteger(iColumn)) {
                //numberInteger++;
                //whichInteger=iColumn;
              }
              largest = std::max(largest, fabs(value));
              if (value > 0.0) {
                if (newUpper[iColumn] >= large) {
                  ++infiniteUpper;
                } else {
                  maximumUp += newUpper[iColumn] * value;
                }
                if (newLower[iColumn] <= -large) {
                  ++infiniteLower;
                } else {
                  maximumDown += newLower[iColumn] * value;
                }
              } else if (value < 0.0) {
                if (newUpper[iColumn] >= large) {
                  ++infiniteLower;
                } else {
                  maximumDown += newUpper[iColumn] * value;
                }
                if (newLower[iColumn] <= -large) {
                  ++infiniteUpper;
                } else {
                  maximumUp += newLower[iColumn] * value;
                }
              }
            } else {
              // fixed
              sumFixed += newLower[iColumn] * value;
            }
          }
          // Adjust
          maximumUp += sumFixed;
          maximumDown += sumFixed;
          // For moment just when all one sign and ints
          //maximumUp += 1.0e-8*fabs(maximumUp);
          //maximumDown -= 1.0e-8*fabs(maximumDown);
          double gap = 0.0;
          if ((rowLower[iRow] > maximumDown && largest > rowLower[iRow] - maximumDown) && ((maximumUp <= rowUpper[iRow] && !infiniteUpper) || rowUpper[iRow] >= 1.0e30)) {
            gap = rowLower[iRow] - maximumDown;
            if (infiniteLower)
              gap = 0.0; // switch off
          } else if ((maximumUp > rowUpper[iRow] && largest > maximumUp - rowUpper[iRow]) && ((maximumDown >= rowLower[iRow] && !infiniteLower) || rowLower[iRow] <= -1.0e30)) {
            gap = -(maximumUp - rowUpper[iRow]);
            if (infiniteUpper)
              gap = 0.0; // switch off
          }
          if (fabs(gap) > 1.0e-8) {
            for (j = rStart; j < rEnd; ++j) {
              double value = element[j];
              iColumn = column[j];
              double difference = newUpper[iColumn] - newLower[iColumn];
              if (difference > 0.0 && difference <= 1.0) {
                double newValue = value;
                if (value * gap > 0.0 && model.isInteger(iColumn)) {
                  if (fabs(value * difference) > fabs(gap)) {
                    // No need for it to be larger than
                    newValue = gap / difference;
                  }
                  if (fabs(value - newValue) > 1.0e-12) {
                    numberChanges++;
                    // BUT variable may have bound
                    double rhsAdjust = 0.0;
                    if (gap > 0.0) {
                      // rowLower
                      if (value > 0.0) {
                        // new value is based on going up from lower bound
                        if (colLower[iColumn])
                          rhsAdjust = colLower[iColumn] * (value - newValue);
                      } else {
                        // new value is based on going down from upper bound
                        if (colUpper[iColumn])
                          rhsAdjust = colUpper[iColumn] * (value - newValue);
                      }
                    } else {
                      // rowUpper
                      if (value < 0.0) {
                        // new value is based on going up from lower bound
                        if (colLower[iColumn])
                          rhsAdjust = colLower[iColumn] * (value - newValue);
                      } else {
                        // new value is based on going down from upper bound
                        if (colUpper[iColumn])
                          rhsAdjust = colUpper[iColumn] * (value - newValue);
                      }
                    }
                    if (rhsAdjust) {
#if CBC_USEFUL_PRINTING > 1
                      printf("FFor column %d bounds %g, %g on row %d bounds %g, %g coefficient was changed from %g to %g with rhs adjustment of %g\n",
                        iColumn, colLower[iColumn], colUpper[iColumn],
                        iRow, rowLower[iRow], rowUpper[iRow],
                        value, newValue, rhsAdjust);
#endif
                      if (rowLower[iRow] > -1.0e20)
                        model.setRowLower(iRow, rowLower[iRow] - rhsAdjust);
                      if (rowUpper[iRow] < 1.0e20)
                        model.setRowUpper(iRow, rowUpper[iRow] - rhsAdjust);
#if CBC_USEFUL_PRINTING > 1
                      printf("FFor column %d bounds %g, %g on row %d bounds %g, %g coefficient was changed from %g to %g with rhs adjustment of %g\n",
                        iColumn, colLower[iColumn], colUpper[iColumn],
                        iRow, rowLower[iRow], rowUpper[iRow],
                        value, newValue, rhsAdjust);
#endif
                    }
                    element[j] = newValue;
                    handler_->message(CGL_ELEMENTS_CHANGED2, messages_)
                      << iRow << iColumn << value << newValue
                      << CoinMessageEol;
#if DEBUG_PREPROCESS > 1
                    const OsiRowCutDebugger *debugger = model.getRowCutDebugger();
                    if (debugger && debugger->numberColumns() == numberColumns) {
                      const double *optimal = debugger->optimalSolution();
                      double sum = 0.0;
                      for (int jj = rStart; jj < rEnd; ++jj) {
                        double value = element[j];
                        int jColumn = column[jj];
                        sum += value * optimal[jColumn];
                      }
                      assert(sum >= rowLower[iRow] - 1.0e7 && sum <= rowUpper[iRow] + 1.0e-7);
                    }
#endif
                  }
                }
              }
            }
          }
	} else {
	  // free row!
	  nFreed++;
        }
      }
      if (numberChanges) {
        handler_->message(CGL_ELEMENTS_CHANGED1, messages_)
          << numberChanges << CoinMessageEol;
        model.replaceMatrixOptional(copy);
      }
    }
  }
  if (numberInfeasible) {
    // restore column bounds
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      model.setColLower(iColumn, columnLower[iColumn]);
      model.setColUpper(iColumn, columnUpper[iColumn]);
    }
  } else if ((options_&128)!=0) {
    // check SC
  }
  delete[] newLower;
  delete[] newUpper;
  delete[] columnLower;
  delete[] columnUpper;
  if (numberInfeasible) {
    return -numberInfeasible;
  } else {
    if (nFreed) 
      totalTightened += 0x40000000+nFreed;
    return (totalTightened);
  }
}
#define CGL_REASONABLE_INTEGER_BOUND 1.23456789e10
// This tightens column bounds (and can declare infeasibility)
// It may also declare rows to be redundant
static int 
tighten(double *colLower, double * colUpper,
	const int *column, const double *rowElements, 
	const CoinBigIndex *rowStart, 
	const CoinBigIndex * rowStartPos,const int * rowLength,
	double *rowLower, double *rowUpper, 
	int nRows,int nCols,char * intVar,int maxpass,
	double tolerance)
{
  int i, j;
  CoinBigIndex k, krs, kre;
  int dolrows;
  int iflagu, iflagl;
  int ntotal=0,nchange=1,jpass=0;
  double dmaxup, dmaxdown, dbound;
  int ninfeas=0;
  // Later try with cliques????
  assert (rowStartPos);
  while(nchange) {
    nchange = 0; 
    if (jpass==maxpass) break;
    jpass++;
    dolrows = (jpass & 1) == 1;
    
    for (i = 0; i < nRows; ++i) {
      if (rowLower[i]>-1.0e20||rowUpper[i]<1.0e20) {
	int iflagu = 0;
	int iflagl = 0;
	double dmaxup = 0.0;
	double dmaxdown = 0.0;
	CoinBigIndex krs = rowStart[i];
	CoinBigIndex krs2 = rowStartPos[i];
	CoinBigIndex kre = rowStart[i]+rowLength[i];
	
	/* ------------------------------------------------------------*/
	/* Compute L(i) and U(i) */
	/* ------------------------------------------------------------*/
	for (k = krs; k < krs2; ++k) {
	  double value=rowElements[k];
	  int j = column[k];
	  if (colUpper[j] < 1.0e12) 
	    dmaxdown += colUpper[j] * value;
	  else
	    ++iflagl;
	  if (colLower[j] > -1.0e12) 
	    dmaxup += colLower[j] * value;
	  else
	    ++iflagu;
	}
	for (k = krs2; k < kre; ++k) {
	  double value=rowElements[k];
	  int j = column[k];
	  if (colUpper[j] < 1.0e12) 
	    dmaxup += colUpper[j] * value;
	  else
	    ++iflagu;
	  if (colLower[j] > -1.0e12) 
	    dmaxdown += colLower[j] * value;
	  else
	    ++iflagl;
	}
	if (iflagu)
	  dmaxup=1.0e31;
	if (iflagl)
	  dmaxdown=-1.0e31;
	if (dmaxup <= rowUpper[i] + tolerance && dmaxdown >= rowLower[i] - tolerance) {
	  /*
	   * The sum of the column maxs is at most the row ub, and
	   * the sum of the column mins is at least the row lb;
	   * this row says nothing at all.
	   * I suspect that this corresponds to
	   * an implied column singleton in the paper (viii, on p. 325),
	   * where the singleton in question is the row slack.
	   */
	  ++nchange;
	  rowLower[i]=-COIN_DBL_MAX;
	  rowUpper[i]=COIN_DBL_MAX;
	} else {
	  if (dmaxup < rowLower[i] -tolerance || dmaxdown > rowUpper[i]+tolerance) {
	    ninfeas++;
	    break;
	  }
	  /*        Finite U(i) */
	  /* -------------------------------------------------------------*/
	  /* below is deliberate mistake (previously was by chance) */
	  /*        never do both */
	  if (iflagu == 0 && rowLower[i] > 0.0 && iflagl == 0 && rowUpper[i] < 1e15) {
	    if (dolrows) {
	      iflagu = 1;
	    } else {
	      iflagl = 1;
	    }
	  }
	  if (iflagu == 0 && rowLower[i] > -1e15) {
	    for (k = krs; k < kre; ++k) {
	      double value=rowElements[k];
	      j = column[k];
	      if (value > 0.0) {
		if (colUpper[j] < 1.0e12) {
		  dbound = colUpper[j] + (rowLower[i] - dmaxup) / value;
		  if (dbound > colLower[j] + 1.0e-8) {
		    /* we can tighten the lower bound */
		    /* the paper mentions this as a possibility on p. 227 */
		    colLower[j] = dbound;
		    ++nchange;
		    
		    /* this may have fixed the variable */
		    /* I believe that this roughly corresponds to a
		     * forcing constraint in the paper (p. 226).
		     * If there is a forcing constraint (with respect
		     * to the original, untightened bounds), then in this 
		     * loop we will go through all the columns and fix
		     * each of them to their implied bound, rather than
		     * determining that the row as a whole is forced
		     * and just fixing them without doing computation for
		     * each column (as in the paper).
		     * By doing it this way, we can tighten bounds and
		     * get forcing constraints for free.
		     */
		    if (colUpper[j] - colLower[j] <= tolerance) {
		      /* --------------------------------------------------*/
		      /*                check if infeasible !!!!! */
		      /* --------------------------------------------------*/
		      if (colUpper[j] - colLower[j] < -100.0*tolerance) {
			ninfeas++;
		      }
		    }
		  }
		}
	      } else {
		if (colLower[j] > -1.0e12) {
		  dbound = colLower[j] + (rowLower[i] - dmaxup) / value;
		  if (dbound < colUpper[j] - 1.0e-8) {
		    colUpper[j] = dbound;
		    ++nchange;
		    if (colUpper[j] - colLower[j] <= tolerance) {
		      /* --------------------------------------------------*/
		      /*                check if infeasible !!!!! */
		      /* --------------------------------------------------*/
		      if (colUpper[j] - colLower[j] < -100.0*tolerance) {
			ninfeas++;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  
	  /* ----------------------------------------------------------------*/
	  /*        Finite L(i) */
	  /* ----------------------------------------------------------------*/
	  if (iflagl == 0 && rowUpper[i] < 1e15) {
	    for (k = krs; k < kre; ++k) {
	      double value=rowElements[k];
	      j = column[k];
	      if (value < 0.0) {
		if (colUpper[j] < 1.0e12) {
		  dbound = colUpper[j] + (rowUpper[i] - dmaxdown) / value;
		  if (dbound > colLower[j] + 1.0e-8) {
		    colLower[j] = dbound;
		    ++nchange;
		    if (! (colUpper[j] - colLower[j] > tolerance)) {
		      /* --------------------------------------------------*/
		      /*                check if infeasible !!!!! */
		      /* --------------------------------------------------*/
		      if (colUpper[j] - colLower[j] < -100.0*tolerance) {
			ninfeas++;
		      }
		    }
		  }
		} 
	      } else {
		if (colLower[j] > -1.0e12) {
		  dbound = colLower[j] + (rowUpper[i] - dmaxdown) / value;
		  if (dbound < colUpper[j] - 1.0e-8) {
		    colUpper[j] = dbound;
		    ++nchange;
		    if (! (colUpper[j] - colLower[j] > tolerance)) {
		      /* --------------------------------------------------*/
		      /*                check if infeasible !!!!! */
		      /* --------------------------------------------------*/
		      if (colUpper[j] - colLower[j] < -100.0*tolerance) {
			ninfeas++;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    for (j = 0; j < nCols; ++j) {
      if (intVar[j]) {
	if (colUpper[j]-colLower[j]>1.0e-8) {
	  if (floor(colUpper[j]+1.0e-4)<colUpper[j]) 
	    nchange++;
	  // clean up anyway
	  colUpper[j]=floor(colUpper[j]+1.0e-4);
	  if (ceil(colLower[j]-1.0e-4)>colLower[j]) 
	    nchange++;
	  // clean up anyway
	  colLower[j]=ceil(colLower[j]-1.0e-4);
	  if (colUpper[j]<colLower[j]) {
	    /*printf("infeasible\n");*/
	    ninfeas++;
	  }
	} else {
	  // clean
	  colUpper[j]=floor(colUpper[j]+1.0e-4);
	  colLower[j]=ceil(colLower[j]-1.0e-4);
	  if (colUpper[j]<colLower[j]) {
	    /*printf("infeasible\n");*/
	    ninfeas++;
	  }
	}
      }
    }
    if (ninfeas) break;
  }
  return (ninfeas);
}

/* Creates solution in original model
   deleteStuff 0 - don't, 1 do (but not if infeasible), 2 always */
void CglPreProcess::postProcess(OsiSolverInterface &modelIn, int deleteStuff)
{
  // Do presolves
  bool saveHint;
  bool solveWithDual = false;
  OsiHintStrength saveStrength;
  originalModel_->getHintParam(OsiDoPresolveInInitial, saveHint, saveStrength);
  bool saveHint2;
  OsiSolverInterface *modelM = &modelIn;
  OsiHintStrength saveStrength2;
  originalModel_->getHintParam(OsiDoDualInInitial,
    saveHint2, saveStrength2);
  OsiSolverInterface *clonedCopy = NULL;
  double saveObjectiveValue = modelM->getObjValue();
  // Make sure cutoff is ignored
  modelM->setDblParam(OsiDualObjectiveLimit, 1.0e50);
  if (!modelM->isProvenOptimal()) {
    CoinWarmStartBasis *slack = dynamic_cast< CoinWarmStartBasis * >(modelM->getEmptyWarmStart());
    modelM->setWarmStart(slack);
    delete slack;
    modelM->resolve();
  }
  double * scBound = NULL;
  if (modelM->isProvenOptimal()) {
#if COIN_DEVELOP > 1
    whichMps++;
    sprintf(nameMps, "start_%d", whichMps);
    modelM->writeMps(nameMps);
    printf("Mps file %s saved in %s at line %d\n",
      nameMps, __FILE__, __LINE__);
#endif
#if CBC_USE_PAPILO
      if (this==initialTry.preProcess && (initialTry.presolveType&128)!=0) {
      initialTry.preProcess = NULL;
      OsiClpSolverInterface * presolvedSolver =
	getClpSolver(modelM);
      delete initialTry.presolvedModel;
      initialTry.presolvedModel = presolvedSolver->getModelPtr();
      ClpSimplex * original = postSolvedModel(initialTry);
      OsiClpSolverInterface * originalSolver = 
	getClpSolver(initialTry.beforePapiloEnd);
      modelM= originalSolver;
      if (originalSolver->integerInformation()) {
	int numberColumns = original->numberColumns();
	for (int i=0;i<numberColumns;i++) {
	  if (original->isInteger(i))
	    originalSolver->setIntegerType(i,1);
	  else
	    originalSolver->setIntegerType(i,0);
	}
      }
      ClpSimplex * simplex = originalSolver->swapModelPtr(original);
    }
#endif
    // See if SC variables
    if ((options_&128)!=0) {
      assert (appData_);
      typedef struct {
	double low;
	double high;
	int column;
      } lotStruct;
      typedef struct {lotStruct * lotsize;int numberLotSizing;} templot;
      templot * temp = reinterpret_cast<templot *>(appData_);
      lotStruct * lotsize=temp->lotsize;
      int numberLotSizing=temp->numberLotSizing;
      int numberColumns = originalModel_->getNumCols();
      scBound = new double [numberColumns];
      for (int i=0;i<numberColumns;i++)
	scBound[i]=-COIN_DBL_MAX;
      for (int i=0;i<numberLotSizing;i++) {
	int iColumn = lotsize[i].column;
	scBound[iColumn] = lotsize[i].low;
      }
      numberColumns = modelM->getNumCols();
      const double *solution = modelM->getColSolution();
      for (int i=0;i<numberColumns;i++) {
	int iColumn = originalColumn_[i];
	if (scBound[iColumn]!=-COIN_DBL_MAX) {
	  double lower =scBound[iColumn];
	  double value = solution[i];
	  if (value<1.0e-5) {
	    modelM->setColUpper(i,0.0);
	  } else {
	    assert (value>lower-1.0e-3);
	    modelM->setColLower(i,lower);
	  }
	}
      }
    }
    if ((options_&256)!=0&&false) {
      int numberColumns = modelM->getNumCols();
      const double *solution = modelM->getColSolution();
      const double *columnLower = modelM->getColLower();
      const double *columnUpper = modelM->getColUpper();
      //OsiSolverInterface * originalModel = originalModel_->clone();
      OsiSolverInterface * originalModel = originalModel_;
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	int jColumn = originalColumn_[iColumn];
	if (modelM->isInteger(iColumn)) {
	  double value = solution[iColumn];
	  double value2 = floor(value + 0.5);
	  // if test fails then empty integer
	  if (fabs(value - value2) < 1.0e-3) {
	    originalModel->setColLower(jColumn, value2);
	    originalModel->setColUpper(jColumn, value2);
	  }
	} else if (columnUpper[iColumn] == columnLower[iColumn]) {
	  originalModel->setColUpper(jColumn, columnLower[iColumn]);
	  originalModel->setColLower(jColumn, columnLower[iColumn]);
	} else if (scBound) {
	  if (scBound[jColumn]!=-COIN_DBL_MAX) {
	    double lower =scBound[jColumn];
	    originalModel->setColLower(jColumn, lower);
          }
        }
      }
      //originalModel->setHintParam(OsiDoReducePrint, false, OsiHintTry);
      originalModel->initialSolve();
      if (deleteStuff) {
	for (int iPass = numberSolvers_ - 1; iPass >= 0; iPass--) {
	  delete modifiedModel_[iPass];
	  ;
	  delete model_[iPass];
	  ;
	  delete presolve_[iPass];
	  modifiedModel_[iPass] = NULL;
	  model_[iPass] = NULL;
	  presolve_[iPass] = NULL;
	}
      }
      delete [] scBound;
      return;
    }
    // If some cuts add back rows
    if (cuts_.sizeRowCuts()) {
      clonedCopy = modelM->clone();
      modelM = clonedCopy;
      int numberRowCuts = cuts_.sizeRowCuts();
      // add in
      CoinBuild build;
      for (int k = 0; k < numberRowCuts; k++) {
        const OsiRowCut *thisCut = cuts_.rowCutPointer(k);
        int n = thisCut->row().getNumElements();
        const int *columnCut = thisCut->row().getIndices();
        const double *elementCut = thisCut->row().getElements();
        double lower = thisCut->lb();
        double upper = thisCut->ub();
        build.addRow(n, columnCut, elementCut, lower, upper);
      }
      modelM->addRows(build);
      // basis is wrong
      CoinWarmStartBasis empty;
      modelM->setWarmStart(&empty);
    }
    if (numberSolvers_==99) {
      // was simple presolve
      numberSolvers_ = 1;
    }
    for (int iPass = numberSolvers_ - 1; iPass >= 0; iPass--) {
      OsiSolverInterface *model = model_[iPass];
#ifdef CBC_HAS_CLP
      OsiClpSolverInterface * postsolvedSolver =
	getClpSolver(model);
      if (postsolvedSolver) // make sure can't stop
	postsolvedSolver->getModelPtr()->setMaximumSeconds(-1.0);
#endif
      int * original = NULL;
      if (model->getNumCols()) {
        CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(modelM->getWarmStart());
        if (basis) {
          model->setWarmStart(basis);
          delete basis;
        }
        int numberColumns = modelM->getNumCols();
        const double *solutionM = modelM->getColSolution();
        const double *columnLower2 = model->getColLower();
        const double *columnUpper2 = model->getColUpper();
        const double *columnLower = modelM->getColLower();
        const double *columnUpper = modelM->getColUpper();
	if (scBound) {
	  // looks odd
	  int numberColumns0 =
	    std::max(originalModel_->getNumCols(),model_[0]->getNumCols());
	  original = new int [numberColumns0];
	  for (int i=0;i<numberColumns0;i++)
	    original[i]=i;
	  for (int jPass=0;jPass<=iPass;jPass++) {
	    OsiPresolve * pinfo = presolve_[jPass];
	    const int *forward = pinfo->originalColumns();
	    OsiSolverInterface * solver = model_[jPass];
	    int numberColumns = solver->getNumCols();
	    for (int i=0;i<numberColumns;i++) {
	      int iColumn = forward[i];
	      original[i] = original[iColumn];
	    }
	  }
	}
        int iColumn;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (modelM->isInteger(iColumn)) {
            double value = solutionM[iColumn];
            double value2 = floor(value + 0.5);
            // if test fails then empty integer
            if (fabs(value - value2) < 1.0e-3) {
              model->setColLower(iColumn, value2);
              model->setColUpper(iColumn, value2);
            }
          } else if (columnUpper[iColumn] == columnLower[iColumn] && !scBound) {
            if (columnUpper2[iColumn] > columnLower2[iColumn] && !model->isInteger(iColumn)) {
              model->setColUpper(iColumn, columnLower[iColumn]);
              model->setColLower(iColumn, columnLower[iColumn]);
            }
	  } else if (scBound) {
            double value = solutionM[iColumn];
	    int jColumn = original[iColumn];
	    if (scBound[jColumn]==-COIN_DBL_MAX) {
	      if (columnUpper[iColumn] == columnLower[iColumn]) {
		if (columnUpper2[iColumn] > columnLower2[iColumn] && !model->isInteger(iColumn)) {
		  model->setColUpper(iColumn, columnLower[iColumn]);
		  model->setColLower(iColumn, columnLower[iColumn]);
		}
	      }
	    } else {
	      double lower =scBound[jColumn];
	      assert (value<1.0e-5||value>lower-1.0e-5);
	      if (value<1.0e-5) 
		model->setColUpper(iColumn,0.0);
	      else
		model->setColLower(iColumn, lower);
	    }
          }
        }
        // IMPORTANT: Copy solution from modelM to model.
        // The loop above only fixed bounds for integers, but model still has the old
        // LP relaxation solution. We must copy the values from modelM (which has the
        // correct integer solution) to model.
        {
          int nc = model->getNumCols();
          double *newSol = new double[nc];
          const double *oldSol = model->getColSolution();
          const double *solM = modelM->getColSolution();
          const double *lo = model->getColLower();
          const double *hi = model->getColUpper();
          for (int i = 0; i < nc; i++) {
            if (i < numberColumns) {
              // This column exists in modelM, use its value
              newSol[i] = solM[i];
            } else {
              // This column doesn't exist in modelM (modelM has fewer columns),
              // use the old value but clip to bounds
              newSol[i] = oldSol[i];
            }
            // Clip to bounds
            newSol[i] = std::max(newSol[i], lo[i]);
            newSol[i] = std::min(newSol[i], hi[i]);
          }
          model->setColSolution(newSol);
          delete[] newSol;
        }
      }
      // After the fix loop has tightened bounds (fixing integers), the warm-start
      // basis stored on model may be stale: variables that are now fixed (lo==hi)
      // could be BASIC in the warm start at fractional LP-relaxation values. CLP's
      // dual simplex would then declare the basis "optimal" in 0 iterations without
      // checking primal feasibility, returning the old LP relaxation solution.
      // Fix: mark all fixed (lo==hi) variables as non-basic (atLowerBound) in the
      // warm start so CLP will re-solve from a primal-feasible point.
      {
        CoinWarmStartBasis *wsb = dynamic_cast<CoinWarmStartBasis*>(model->getWarmStart());
        if (wsb) {
          const double *lo2 = model->getColLower(), *hi2 = model->getColUpper();
          int nc = model->getNumCols();
          for (int k = 0; k < nc; k++) {
            if (lo2[k] >= hi2[k] - 1e-10) {
              wsb->setStructStatus(k, CoinWarmStartBasis::atLowerBound);
            }
          }
          model->setWarmStart(wsb);
          delete wsb;
        }
      }
      int numberColumns = modelM->getNumCols();
      const double *solutionM = modelM->getColSolution();
      int iColumn;
      // Give a hint to do primal
      //model->setHintParam(OsiDoPresolveInInitial,true,OsiHintTry);
      model->setHintParam(OsiDoDualInInitial, false, OsiHintTry);
      // clean
      /*
  VIRTUOUS - I am not happy here (JJF) - This was put in for Special Ordered Sets of type 2

  Previous loop has likely made nontrivial bound changes, hence invalidated
  solution. Why do we need this? We're about to do an initialSolve, which
  will overwrite solution. Perhaps belongs in same guarded block with
  following feasibility check? If this is necessary for clp, solution should
  be acquired before bounds changes.
*/
      if (0) {
        int numberColumns = model->getNumCols();
        const double *lower = model->getColLower();
        const double *upper = model->getColUpper();
        double *solution = CoinCopyOfArray(model->getColSolution(), numberColumns);
        int i;
        for (i = 0; i < numberColumns; i++) {
          double value = solution[i];
          value = std::min(value, upper[i]);
          value = std::max(value, lower[i]);
          solution[i] = value;
        }
        model->setColSolution(solution);
        delete[] solution;
      }
#if CBC_USEFUL_PRINTING > 1
      {
        int numberColumns = model->getNumCols();
        int numberRows = model->getNumRows();
        const double *solution = model->getColSolution();
        const double *lower = model->getColLower();
        const double *upper = model->getColUpper();
        const double *rowLower = model->getRowLower();
        const double *rowUpper = model->getRowUpper();
#ifndef NDEBUG
        double primalTolerance = 1.0e-8;
#endif
        // Column copy
        const CoinPackedMatrix *matrix = model->getMatrixByCol();
        const double *element = matrix->getElements();
        const int *row = matrix->getIndices();
        const CoinBigIndex *columnStart = matrix->getVectorStarts();
        const int *columnLength = matrix->getVectorLengths();
        double *rowActivity = new double[numberRows];
        memset(rowActivity, 0, numberRows * sizeof(double));
        int i;
        for (i = 0; i < numberColumns; i++) {
          int j;
          double value = solution[i];
          if (value < lower[i]) {
            value = lower[i];
          } else if (value > upper[i]) {
            value = upper[i];
          }
          assert(upper[i] >= lower[i]);
          assert((fabs(value) < 1.0e20));
          if (value) {
            for (j = columnStart[i];
                 j < columnStart[i] + columnLength[i]; j++) {
              int iRow = row[j];
              rowActivity[iRow] += value * element[j];
            }
          }
        }
        // check was feasible - if not adjust (cleaning may move)
        bool feasible = true;
        for (i = 0; i < numberRows; i++) {
          if (rowActivity[i] < rowLower[i]) {
            if (rowActivity[i] < rowLower[i] - 1000.0 * primalTolerance) {
              feasible = false;
#if CBC_USEFUL_PRINTING
              printf("Bad row %d %g <= %g <= %g\n",
                i, rowLower[i], rowActivity[i], rowUpper[i]);
#endif
            }
            rowActivity[i] = rowLower[i];
          } else if (rowActivity[i] > rowUpper[i]) {
            if (rowActivity[i] > rowUpper[i] + 1000.0 * primalTolerance) {
              feasible = false;
#if CBC_USEFUL_PRINTING
              printf("Bad row %d %g <= %g <= %g\n",
                i, rowLower[i], rowActivity[i], rowUpper[i]);
#endif
            }
            rowActivity[i] = rowUpper[i];
          }
        }
        if (!feasible)
	  printf("ZZ not feasible??\n");
      }
#endif
      {
        int numberFixed = 0;
        int numberColumns = model->getNumCols();
        const double *columnLower = model->getColLower();
        const double *columnUpper = model->getColUpper();
        int iColumn;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (columnLower[iColumn] == columnUpper[iColumn])
            numberFixed++;
        }
        if (numberColumns > 2000 && numberFixed < numberColumns && numberFixed * 5 > numberColumns) {
          model->setHintParam(OsiDoPresolveInInitial, true, OsiHintTry);
        }
      }
      model->setHintParam(OsiDoDualInInitial, true, OsiHintTry);
      model->setHintParam(OsiDoPresolveInInitial, false, OsiHintTry);
      model->initialSolve();
      numberIterationsPost_ += model->getIterationCount();
      if (!model->isProvenOptimal()) {
        // try without basis
        CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(model->getEmptyWarmStart());
        model->setWarmStart(basis);
        delete basis;
        model->initialSolve();
      }
      if (!model->isProvenOptimal()) {
#if COIN_DEVELOP
        whichMps++;
        sprintf(nameMps, "bad2_%d", whichMps);
        model->writeMps(nameMps);
        printf("Mps file %s saved in %s at line %d\n",
          nameMps, __FILE__, __LINE__);
        printf("bad unwind in postprocess\n");
        OsiSolverInterface *temp = model->clone();
        temp->setDblParam(OsiDualObjectiveLimit, 1.0e30);
        temp->setHintParam(OsiDoReducePrint, false, OsiHintTry);
        temp->initialSolve();
        if (temp->isProvenOptimal()) {
          printf("Was infeasible on objective limit\n");
        }
        delete temp;
#endif
      } else {
#if COIN_DEVELOP > 1
        whichMps++;
        sprintf(nameMps, "good2_%d", whichMps);
        model->writeMps(nameMps);
        printf("Mps file %s saved in %s at line %d\n",
          nameMps, __FILE__, __LINE__);
#endif
      }
      const int *originalColumns = presolve_[iPass]->originalColumns();
      const double *columnLower = modelM->getColLower();
      const double *columnUpper = modelM->getColUpper();
      OsiSolverInterface *modelM2;
      if (iPass)
        modelM2 = modifiedModel_[iPass - 1];
      else
        modelM2 = startModel_;
      const double *solutionM2 = modelM2->getColSolution();
      const double *columnLower2 = modelM2->getColLower();
      const double *columnUpper2 = modelM2->getColUpper();
      double primalTolerance;
      modelM->getDblParam(OsiPrimalTolerance, primalTolerance);
      /* clean up status for any bound alterations made by preprocess which 
	 postsolve won't understand.
	 Could move inside OsiPresolve but some people might object */
      CoinWarmStartBasis *presolvedBasis = dynamic_cast< CoinWarmStartBasis * >(model->getWarmStart());
      assert(presolvedBasis);
      int numberChanged = 0;
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        int jColumn = originalColumns[iColumn];
        switch (presolvedBasis->getStructStatus(iColumn)) {
        case CoinWarmStartBasis::basic:
        case CoinWarmStartBasis::superBasic:
        case CoinWarmStartBasis::isFree:
          break;
        case CoinWarmStartBasis::atLowerBound:
          if (solutionM[iColumn] > columnLower2[jColumn] + primalTolerance) {
	    if (columnLower2[jColumn]<-1.0e50&&
		columnUpper2[iColumn]>1.0e50)
	      presolvedBasis->setStructStatus(iColumn, CoinWarmStartBasis::isFree);
	    else
	      presolvedBasis->setStructStatus(iColumn, CoinWarmStartBasis::superBasic);
            numberChanged++;
          }
          break;
        case CoinWarmStartBasis::atUpperBound:
          if (solutionM[iColumn] < columnUpper2[jColumn] - primalTolerance) {
	    if (columnLower2[jColumn]<-1.0e50&&
		columnUpper2[iColumn]>1.0e50)
	      presolvedBasis->setStructStatus(iColumn, CoinWarmStartBasis::isFree);
	    else
	      presolvedBasis->setStructStatus(iColumn, CoinWarmStartBasis::superBasic);
            numberChanged++;
          }
          break;
        }
      }
      if (numberChanged) {
	// say do primal
	solveWithDual = false;
        model->setWarmStart(presolvedBasis);
      }
      delete presolvedBasis;
      presolve_[iPass]->postsolve(true);
      // and fix values
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        int jColumn = originalColumns[iColumn];
        if (!modelM2->isInteger(jColumn)) {
          if (columnUpper[iColumn] == columnLower[iColumn]) {
            if (columnUpper2[jColumn] > columnLower2[jColumn] && !modelM2->isInteger(jColumn)) {
              double value = solutionM[iColumn];
              value = std::max(value, columnLower[iColumn]);
              value = std::min(value, columnUpper[iColumn]);
#if CBC_USEFUL_PRINTING > 1
              printf("assuming %d fixed to solution of %g (was %g) - bounds %g %g, old bounds and sol %g %g\n",
                jColumn, value, solutionM2[jColumn], columnLower2[jColumn], columnUpper2[jColumn],
                columnLower[iColumn], columnUpper[iColumn]);
#endif
              modelM2->setColLower(jColumn, value);
              modelM2->setColUpper(jColumn, value);
            }
          } else {
#if CBC_USEFUL_PRINTING
            if (columnUpper2[jColumn] > columnLower2[jColumn] && !modelM2->isInteger(jColumn)) {
              double value = solutionM[iColumn];
              value = std::max(value, columnLower[iColumn]);
              value = std::min(value, columnUpper[iColumn]);
              printf("assuming %d not fixed to solution of %g (was %g) - bounds %g %g, old bounds and sol %g %g\n",
                jColumn, value, solutionM2[jColumn], columnLower2[jColumn], columnUpper2[jColumn],
                columnLower[iColumn], columnUpper[iColumn]);
            }
#endif
	    if (scBound) {
              double value = solutionM[iColumn];
	      int jColumn2 = original[iColumn];
	      if (scBound[jColumn2]!=-COIN_DBL_MAX) {
		double lower =scBound[jColumn2];
		assert (value<1.0e-5||value>lower-1.0e-5);
		if (value<1.0e-5) { 
		  modelM2->setColUpper(jColumn,0.0);
		} else {
		  modelM2->setColLower(jColumn, lower);
		  //modelM2->setColUpper(iColumn, value);
		}
	      }
	    }
          }
        } else {
          // integer - dupcol bounds may be odd so use solution
          double value = floor(solutionM2[jColumn] + 0.5);
          if (value < columnLower2[jColumn]) {
            //printf("changing lower bound for %d from %g to %g to allow feasibility\n",
            //	   jColumn,columnLower2[jColumn],value);
            modelM2->setColLower(jColumn, value);
          } else if (value > columnUpper2[jColumn]) {
            //printf("changing upper bound for %d from %g to %g to allow feasibility\n",
            //	   jColumn,columnUpper2[jColumn],value);
            modelM2->setColUpper(jColumn, value);
          }
        }
      }
      if (deleteStuff) {
        delete modifiedModel_[iPass];
        ;
        delete model_[iPass];
        ;
        delete presolve_[iPass];
        modifiedModel_[iPass] = NULL;
        model_[iPass] = NULL;
        presolve_[iPass] = NULL;
      }
      // After postsolve, presolve back-substitution may have computed fractional
      // values for eliminated integer variables via equations like
      //   x_elim = (rhs - a*x_surv) / b  with non-unit b.
      // These fractional values cascade: the next pass sees them as "unknown"
      // integers, the LP returns fractional solutions, and the bad values
      // propagate all the way to startModel_.
      //
      // Fix: for each fractional integer in modelM2 try rounding down then
      // rounding up, checking row activities to ensure feasibility.  Only
      // commit a rounding if at least one direction is row-feasible; update
      // row activities after each committed rounding so subsequent roundings
      // see the correct state (same greedy strategy as numberBadValues code).
      // Variables for which neither direction is feasible are left fractional
      // and handled by the outer numberBadValues repair pass.
      {
        int nM2Cols = modelM2->getNumCols();
        int nM2Rows = modelM2->getNumRows();
        const double *sM2 = modelM2->getColSolution();
        const double *loM2 = modelM2->getColLower();
        const double *hiM2 = modelM2->getColUpper();
        bool anyFrac = false;
        for (int k = 0; k < nM2Cols && !anyFrac; k++) {
          if (modelM2->isInteger(k) && fabs(sM2[k] - floor(sM2[k] + 0.5)) > 1e-5)
            anyFrac = true;
        }
        if (anyFrac) {
          const double *rowLoM2 = modelM2->getRowLower();
          const double *rowHiM2 = modelM2->getRowUpper();
          const CoinPackedMatrix *colMtx2 = modelM2->getMatrixByCol();
          const double *elems2 = colMtx2->getElements();
          const int *rowIdx2 = colMtx2->getIndices();
          const CoinBigIndex *colStarts2 = colMtx2->getVectorStarts();
          const int *colLens2 = colMtx2->getVectorLengths();
          // Compute initial row activities.
          std::vector<double> rowAct(nM2Rows, 0.0);
          for (int k = 0; k < nM2Cols; k++) {
            double val = sM2[k];
            if (val == 0.0) continue;
            for (CoinBigIndex j = colStarts2[k]; j < colStarts2[k] + colLens2[k]; j++)
              rowAct[rowIdx2[j]] += val * elems2[j];
          }
          std::vector<double> newSol(sM2, sM2 + nM2Cols);
          const double rowTol = 1e-5;
          for (int k = 0; k < nM2Cols; k++) {
            if (!modelM2->isInteger(k)) continue;
            double v = newSol[k];
            double vr = floor(v + 0.5);
            if (fabs(v - vr) <= 1e-5) continue; // already close to integer
            double vDown = std::max(loM2[k], floor(v));
            double vUp   = std::min(hiM2[k], ceil(v));
            // Check round-down feasibility against current row activities.
            bool canDown = (vDown >= loM2[k] - 1e-10);
            if (canDown) {
              double delta = vDown - v;
              for (CoinBigIndex j = colStarts2[k]; j < colStarts2[k] + colLens2[k]; j++) {
                int r = rowIdx2[j];
                double newA = rowAct[r] + delta * elems2[j];
                if (newA < rowLoM2[r] - rowTol || newA > rowHiM2[r] + rowTol) {
                  canDown = false;
                  break;
                }
              }
            }
            // Check round-up feasibility against current row activities.
            bool canUp = (vUp <= hiM2[k] + 1e-10);
            if (canUp) {
              double delta = vUp - v;
              for (CoinBigIndex j = colStarts2[k]; j < colStarts2[k] + colLens2[k]; j++) {
                int r = rowIdx2[j];
                double newA = rowAct[r] + delta * elems2[j];
                if (newA < rowLoM2[r] - rowTol || newA > rowHiM2[r] + rowTol) {
                  canUp = false;
                  break;
                }
              }
            }
            double chosen;
            if (canDown && canUp)
              chosen = (fabs(v - vDown) <= fabs(v - vUp)) ? vDown : vUp;
            else if (canDown)
              chosen = vDown;
            else if (canUp)
              chosen = vUp;
            else
              continue; // neither direction feasible: leave fractional
            // Commit the rounding: update row activities.
            double delta = chosen - v;
            for (CoinBigIndex j = colStarts2[k]; j < colStarts2[k] + colLens2[k]; j++)
              rowAct[rowIdx2[j]] += delta * elems2[j];
            newSol[k] = chosen;
          }
          modelM2->setColSolution(newSol.data());
        }
      }
      modelM = modelM2;
      delete [] original;
    }
#if CBC_USE_PAPILO
      if (this==initialTry.preProcess && (initialTry.presolveType&64)!=0) {
      initialTry.preProcess = NULL;
      OsiClpSolverInterface * presolvedSolver =
	getClpSolver(modelM);
      delete initialTry.presolvedModel;
      initialTry.presolvedModel = presolvedSolver->getModelPtr();
      ClpSimplex * original = postSolvedModel(initialTry);
      OsiClpSolverInterface * originalSolver =
	getClpSolver(originalModel_);
      modelM= originalSolver;
      int numberColumns = original->numberColumns();
      for (int i=0;i<numberColumns;i++) {
	if (original->isInteger(i))
	  originalSolver->setIntegerType(i,1);
	else
	  originalSolver->setIntegerType(i,0);
      }      
      ClpSimplex * simplex = originalSolver->swapModelPtr(original);
      //delete simplex;
    }
#endif
    // should be back to startModel_;
    OsiSolverInterface *model = originalModel_;
    // Use number of columns in original
    int numberColumns = model->getNumCols();
    const double *solutionM = modelM->getColSolution();
    int iColumn;
    const double *columnLower2 = model->getColLower();
    const double *columnUpper2 = model->getColUpper();
    const double *columnLower = modelM->getColLower();
    const double *columnUpper = modelM->getColUpper();
    int numberBadValues = 0;
    CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(modelM->getWarmStart());
    if (basis) {
      model->setWarmStart(basis);
      delete basis;
    }
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (modelM->isInteger(iColumn)) {
        double value = solutionM[iColumn];
        double value2 = floor(value + 0.5);
        // if test fails then empty integer
        if (fabs(value - value2) < 1.0e-3) {
          value2 = std::max(std::min(value2, columnUpper[iColumn]), columnLower[iColumn]);
          model->setColLower(iColumn, value2);
          model->setColUpper(iColumn, value2);
        } else {
#if CBC_USEFUL_PRINTING > 1
          printf("NPASS=%d, ipass end var %d values %g %g %g\n",
            numberSolvers_, iColumn, model->getColLower()[iColumn],
            value, model->getColUpper()[iColumn]);
#endif
          numberBadValues++;
        }
      } else if (columnUpper[iColumn] == columnLower[iColumn]) {
        if (columnUpper2[iColumn] > columnLower2[iColumn] && !model->isInteger(iColumn)) {
          model->setColUpper(iColumn, columnLower[iColumn]);
          model->setColLower(iColumn, columnLower[iColumn]);
        }
      }
    }
    if (numberBadValues) {
      const CoinPackedMatrix *columnCopy = model->getMatrixByCol();
      const int *row = columnCopy->getIndices();
      const CoinBigIndex *columnStart = columnCopy->getVectorStarts();
      const int *columnLength = columnCopy->getVectorLengths();
      const double *element = columnCopy->getElements();
      int numberRows = model->getNumRows();
      double *rowActivity = new double[numberRows];
      memset(rowActivity, 0, numberRows * sizeof(double));
      double *solution = CoinCopyOfArray(solutionM, numberColumns);
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = solutionM[iColumn];
        if (modelM->isInteger(iColumn)) {
          double value2 = floor(value + 0.5);
          // if test fails then empty integer
          if (fabs(value - value2) < 1.0e-3)
            value = value2;
        }
        solution[iColumn] = value;
        for (CoinBigIndex j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          rowActivity[iRow] += value * element[j];
        }
      }
      const double *rowLower = model->getRowLower();
      const double *rowUpper = model->getRowUpper();
      //const double * columnLower = model->getColLower();
      //const double * columnUpper = model->getColUpper();
      const double *objective = model->getObjCoefficients();
      double direction = model->getObjSense();
#ifndef NDEBUG
      int numberCheck = 0;
#endif
      double tolerance;
      model->getDblParam(OsiPrimalTolerance, tolerance);
      tolerance *= 10.0;
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = solution[iColumn];
        if (model->isInteger(iColumn)) {
          double value2 = floor(value);
          // See if empty integer
          if (value != value2) {
#ifndef NDEBUG
            numberCheck++;
#endif
            int allowed = 0;
            // can we go up
            double movement = value2 + 1.0 - value;
            CoinBigIndex j;
            bool good = true;
            for (j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
#if CBC_USEFUL_PRINTING > 1
              if (rowLower[iRow] > -1.0e20 && rowUpper[iRow] < 1.0e20)
                printf("odd row with both bounds %d %g %g - element %g\n",
                  iRow, rowLower[iRow], rowUpper[iRow], element[j]);
#endif
              double newActivity = rowActivity[iRow] + movement * element[j];
              if (newActivity > rowUpper[iRow] + tolerance || newActivity < rowLower[iRow] - tolerance)
                good = false;
            }
            if (good)
              allowed = 1;
            // can we go down
            movement = value2 - value;
            good = true;
            for (j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              double newActivity = rowActivity[iRow] + movement * element[j];
              if (newActivity > rowUpper[iRow] + tolerance || newActivity < rowLower[iRow] - tolerance)
                good = false;
            }
            if (good)
              allowed |= 2;
            if (allowed) {
              if (allowed == 3) {
                if (direction * objective[iColumn] > 0.0)
                  allowed = 2;
                else
                  allowed = 1;
              }
              if (allowed == 1)
                value2++;
              movement = value2 - value;
              solution[iColumn] = value2;
              model->setColLower(iColumn, value2);
              model->setColUpper(iColumn, value2);
              for (j = columnStart[iColumn];
                   j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                int iRow = row[j];
                rowActivity[iRow] += movement * element[j];
              }
            } else {
#if CBC_USEFUL_PRINTING > 1
              printf("Unable to move %d\n", iColumn);
#endif
            }
          }
        }
      }
      // numberBadValues counts fractional integers in startModel_ (which can
      // include CGL-discovered implicit-integer variables not present in
      // originalModel_).  numberCheck counts only those in originalModel_.
      // Since startModel_ is a superset, numberCheck <= numberBadValues always.
      assert(numberCheck <= numberBadValues);
      model->setColSolution(solution);
      delete[] rowActivity;
      delete[] solution;
    }

  } else {
    // infeasible
    // Back to startModel_;
    OsiSolverInterface *model = originalModel_;
    // Use number of columns in original
    int numberColumns = model->getNumCols();
    const double *columnLower = model->getColLower();
    int iColumn;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (model->isInteger(iColumn))
        model->setColUpper(iColumn, columnLower[iColumn]);
    }
  }
  delete clonedCopy;
  originalModel_->setHintParam(OsiDoPresolveInInitial, false, OsiHintTry);
  originalModel_->setHintParam(OsiDoDualInInitial, false, OsiHintTry);
  {
    int numberFixed = 0;
    int numberColumns = originalModel_->getNumCols();
    const double *columnLower = originalModel_->getColLower();
    const double *columnUpper = originalModel_->getColUpper();
    int iColumn;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (columnLower[iColumn] == columnUpper[iColumn]) {
        numberFixed++;
      } else if (scBound) {
	if (scBound[iColumn]!=-COIN_DBL_MAX) {
	  double lower =scBound[iColumn];
	  originalModel_->setColLower(iColumn, lower);
	}
      }
    }
    if (numberColumns < 10000 || numberFixed == numberColumns) {
      CoinWarmStart *empty = originalModel_->getEmptyWarmStart();
      originalModel_->setWarmStart(empty);
      delete empty;
    }
  }
  delete [] scBound;
  //double time1 = CoinCpuTime();
#ifdef CBC_HAS_CLP
  OsiClpSolverInterface * originalSolver =
    getClpSolver(originalModel_);
#ifndef CBC_SKIP_CLP_TEST
  if (originalSolver) // make sure can't stop
#endif
    originalSolver->getModelPtr()->setMaximumSeconds(-1.0);
#endif
  originalModel_->initialSolve();
  numberIterationsPost_ += originalModel_->getIterationCount();
  double objectiveValue = originalModel_->getObjValue();
  double testObj = 1.0e-8 * std::max(fabs(saveObjectiveValue), fabs(objectiveValue)) + 1.0e-4;
  if (!originalModel_->isProvenOptimal()) {
#if COIN_DEVELOP
    whichMps++;
    sprintf(nameMps, "bad3_%d", whichMps);
    originalModel_->writeMps(nameMps);
    printf("Mps file %s saved in %s at line %d\n",
      nameMps, __FILE__, __LINE__);
    printf("bad end unwind in postprocess\n");
#endif
    handler_->message(CGL_POST_INFEASIBLE, messages_)
      << CoinMessageEol;
    if (deleteStuff) {
      for (int iPass = numberSolvers_ - 1; iPass >= 0; iPass--) {
        delete modifiedModel_[iPass];
        ;
        delete model_[iPass];
        ;
        delete presolve_[iPass];
        modifiedModel_[iPass] = NULL;
        model_[iPass] = NULL;
        presolve_[iPass] = NULL;
      }
    }
  } else if (fabs(saveObjectiveValue - objectiveValue) > testObj
    && deleteStuff) {
    // Only warn when postprocessing makes the solution worse; an improvement
    // is expected when the MIP solution was found by a heuristic (e.g.
    // Feasibility Jump) that does not re-optimise continuous variables.
    double direction = originalModel_->getObjSense(); // +1 min, -1 max
    if (direction * (objectiveValue - saveObjectiveValue) > testObj) {
      handler_->message(CGL_POST_CHANGED, messages_)
        << saveObjectiveValue << objectiveValue
        << CoinMessageEol;
    }
  }
  originalModel_->setHintParam(OsiDoDualInInitial, saveHint2, saveStrength2);
  originalModel_->setHintParam(OsiDoPresolveInInitial, saveHint, saveStrength);
}
//-------------------------------------------------------------------
// Returns the greatest common denominator of two
// positive integers, a and b, found using Euclid's algorithm
//-------------------------------------------------------------------
static int gcd(int a, int b)
{
  int remainder = -1;
  // make sure a<=b (will always remain so)
  if (a > b) {
    // Swap a and b
    int temp = a;
    a = b;
    b = temp;
  }
  // if zero then gcd is nonzero (zero may occur in rhs of packed)
  if (!a) {
    if (b) {
      return b;
    } else {
      printf("**** gcd given two zeros!!\n");
      abort();
    }
  }
  while (remainder) {
    remainder = b % a;
    b = a;
    a = remainder;
  }
  return b;
}
#define CGL_PREPROCESS_DENSE_CODE
#define F77_FUNC(x, y) x##_
/* Type of Fortran integer translated into C */
#ifndef ipfint
//typedef ipfint FORTRAN_INTEGER_TYPE ;
typedef int ipfint;
typedef const int cipfint;
#endif
//#define COIN_HAS_LAPACK
//#include "CoinFactorization.hpp"
#ifdef CGL_PREPROCESS_DENSE_CODE
// using simple lapack interface
extern "C" {
/** LAPACK Fortran subroutine DGETRF. */
void F77_FUNC(dgetrf, DGETRF)(ipfint *m, ipfint *n,
  double *A, ipfint *ldA,
  ipfint *ipiv, ipfint *info);
/** LAPACK Fortran subroutine DGETRS. */
void F77_FUNC(dgetrs, DGETRS)(char *trans, cipfint *n,
  cipfint *nrhs, const double *A, cipfint *ldA,
  cipfint *ipiv, double *B, cipfint *ldB, ipfint *info,
  int trans_len);
}
#endif
/* Return model with useful modifications.  
   If constraints true then adds any x+y=1 or x-y=0 constraints
   If NULL infeasible
*/
OsiSolverInterface *
CglPreProcess::modified(OsiSolverInterface *model,
  bool constraints,
  int &numberChanges,
  int iBigPass,
  int numberPasses)
{
  OsiSolverInterface *newModel = model->clone();
  int numberRows = newModel->getNumRows();
  CglUniqueRowCuts twoCuts(numberRows);
  int numberColumns = newModel->getNumCols();
  bool solveWithDual = true;
  //int number01Integers = 0;
  //int iColumn;
  //for (iColumn = 0; iColumn < numberColumns; iColumn++) {
  //  if (newModel->isBinary(iColumn))
  //    number01Integers++;
  //}
#if CBC_USEFUL_PRINTING > 0
  printf("At beginning of modified %d rows and %d columns (pass %d - doing %d minor passes)\n",
    numberRows, numberColumns, iBigPass, numberPasses);
#endif
  if (inspect_) {
    int nInts = 0;
    for (int i = 0; i < numberColumns; i++)
      if (newModel->isInteger(i))
        nInts++;
    char buf[256];
    snprintf(buf, sizeof(buf),
      "[Preproc major pass %d] modified(): %d rows, %d cols (%d integer), %d NZ, %d minor passes",
      iBigPass, numberRows, numberColumns, nInts, newModel->getNumElements(), numberPasses);
    FILE *fp = handler_->filePointer();
    if (fp) { fprintf(fp, "  %s\n", buf); fflush(fp); }
  }
  OsiRowCut **whichCut = new OsiRowCut *[numberRows + 1];
  memset(whichCut, 0, (numberRows + 1) * sizeof(OsiRowCut *));
  numberChanges = 0;
  CoinThreadRandom randomGenerator;
  CglTreeProbingInfo info(model);
  info.level = 0;
  info.pass = 0;
  info.formulation_rows = numberRows;
  info.inTree = false;
  info.options = !numberProhibited_ ? 0 : 2;
  info.randomNumberGenerator = &randomGenerator;
  info.strengthenRow = whichCut;
#ifdef HEAVY_PROBING
  // See if user asked for heavy probing
  bool heavyProbing = false;
  for (int iGenerator = 0; iGenerator < numberCutGenerators_; iGenerator++) {
    CglProbing *probingCut = dynamic_cast< CglProbing * >(generator_[iGenerator]);
    if (probingCut && probingCut->getMaxPassRoot() > 1) {
      heavyProbing = true;
      break;
    }
  }
#endif
  bool feasible = true;
  int firstGenerator = 0;
  int lastGenerator = numberCutGenerators_;
  bool useSolution = getCutoff() < 1.0e20;
#if 1
  // Do triple stuff
  if (iBigPass == 0) {
    // Row copy
    CoinPackedMatrix matrixByRow(*newModel->getMatrixByRow());
    const double *elementByRow = matrixByRow.getElements();
    const int *column = matrixByRow.getIndices();
    const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
    const int *rowLength = matrixByRow.getVectorLengths();

    // Column copy
    CoinPackedMatrix matrixByCol(*newModel->getMatrixByCol());
    //const double * element = matrixByCol.getElements();
    const int *row = matrixByCol.getIndices();
    const CoinBigIndex *columnStart = matrixByCol.getVectorStarts();
    const int *columnLength = matrixByCol.getVectorLengths();

    const double *rowLower = newModel->getRowLower();
    const double *rowUpper = newModel->getRowUpper();
    const double *columnLower = newModel->getColLower();
    const double *columnUpper = newModel->getColUpper();
#define TRIPLE_ROWS 4
#define TRIPLE_COLS 3
    // just allow TRIPLE_ROWS rows
    double el[TRIPLE_COLS][2 * TRIPLE_ROWS];
    double rhs[2 * TRIPLE_ROWS], lower[TRIPLE_ROWS], upper[TRIPLE_ROWS];
    double modifiedRhs[2 * TRIPLE_ROWS], scaleFactor[2 * TRIPLE_ROWS];
    memset(modifiedRhs, 0, sizeof(modifiedRhs));
    memset(scaleFactor, 0, sizeof(scaleFactor));
    int rowNumber[2 * TRIPLE_ROWS];
    double colLower[TRIPLE_COLS], colUpper[TRIPLE_COLS];
#if CBC_USEFUL_PRINTING > 0
    int binary[TRIPLE_COLS];
    int colNumber[TRIPLE_COLS];
#endif
    unsigned char result[2 * TRIPLE_ROWS][2 * TRIPLE_ROWS];
    int rowType[2 * TRIPLE_ROWS];
#define MAX_ELS 20
#define MAX_LOOK_ROWS 40
#define MAX_LOOK_COLS 32 // for more go to long
    unsigned int bitMask[MAX_LOOK_ROWS];
    int *which = new int[2 * numberColumns + 3 * numberRows + MAX_LOOK_ROWS + MAX_LOOK_COLS];
    int *marked = which + MAX_LOOK_COLS;
    int *whichRow = marked + numberColumns;
    int *markRow = whichRow + MAX_LOOK_ROWS;
    int *backwardRow = markRow + numberRows;
    int *backwardColumn = backwardRow + numberRows;
    int *rowTypeAll = backwardColumn + numberColumns;
    memset(marked, 0, numberColumns * sizeof(int));
    memset(markRow, 0, numberRows * sizeof(int));
    for (int iRow = 0; iRow < numberRows; iRow++) {
      int type = 0;
      if (rowLower[iRow] != rowUpper[iRow]) {
        if (rowLower[iRow] > -1.0e30) {
          if (rowUpper[iRow] < 1.0e30) {
            type = 2;
          } else {
            type = 1;
          }
        } else {
          if (rowUpper[iRow] < 1.0e30) {
            type = 2;
          } else {
            type = 0;
          }
        }
      } else {
        type = 3;
      }
      if (rowType_ && rowType_[iRow])
        type = 0; // not allowed if may be cut
      rowTypeAll[iRow] = type;
    }
    // clean model
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (newModel->isInteger(iColumn)) {
        double lower = columnLower[iColumn];
        double upper = columnUpper[iColumn];
        double saveLower = lower;
        double saveUpper = upper;
        if (lower != ceil(lower)) {
          if (fabs(lower - floor(lower + 0.5)) < 1.0e-6)
            lower = floor(lower + 0.5);
          else
            lower = ceil(lower);
          newModel->setColLower(iColumn, lower);
        }
        if (upper != floor(upper)) {
          if (fabs(upper - floor(upper + 0.5)) < 1.0e-6)
            upper = floor(upper + 0.5);
          else
            upper = floor(upper);
          newModel->setColUpper(iColumn, upper);
        }
        if (lower > upper) {
          char generalPrint[100];
          sprintf(generalPrint, "%d input bounds %g %g",
            iColumn, saveLower, saveUpper);
          handler_->message(CGL_GENERAL, messages_)
            << generalPrint
            << CoinMessageEol;
          handler_->message(CGL_INFEASIBLE, messages_)
            << CoinMessageEol;
          feasible = false;
        }
        if (lower < 0.0) {
          // take out for now
          for (CoinBigIndex j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            rowTypeAll[iRow] = 0;
          }
        }
      } else {
	// continuous - take out to be on safe side
	for (CoinBigIndex j = columnStart[iColumn];
	     j < columnStart[iColumn] + columnLength[iColumn]; j++) {
	  int iRow = row[j];
	  rowTypeAll[iRow] = 0;
	}
      }
    }
#if CBC_USEFUL_PRINTING > 0
    int nFreed = 0;
#endif
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (!feasible)
        break;
      if (newModel->isInteger(iColumn) && columnLength[iColumn] <= MAX_ELS
        && columnLength[iColumn] > 1 && columnLower[iColumn] == 0.0 && columnUpper[iColumn] == 1.0) {
        int nMarked = 0;
        backwardColumn[iColumn] = 0;
        marked[iColumn] = 1;
        which[nMarked++] = iColumn;
        int nMarkRow = 0;
        for (CoinBigIndex j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          markRow[iRow] = nMarkRow + 1;
          backwardRow[iRow] = nMarkRow;
          bitMask[nMarkRow] = 0;
          if (!rowTypeAll[iRow]) {
            // no good
            bitMask[nMarkRow] = 0xffffffff;
          }
          whichRow[nMarkRow++] = iRow;
        }
        for (CoinBigIndex j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          if (rowLength[iRow] <= 3) {
            for (CoinBigIndex k = rowStart[iRow];
                 k < rowStart[iRow] + rowLength[iRow]; k++) {
              int kColumn = column[k];
              if (!marked[kColumn]) {
                int iMask = backwardRow[iRow];
                if (nMarked == MAX_LOOK_COLS) {
                  // say row no good
                  bitMask[iMask] = 0xffffffff;
                  continue;
                } else if (bitMask[iMask] == 0xffffffff) {
                  continue;
                }
                backwardColumn[kColumn] = nMarked;
                marked[kColumn] = nMarked + 1;
                which[nMarked++] = kColumn;
              }
              int iMask = backwardRow[iRow];
              int iShift = backwardColumn[kColumn];
              bitMask[iMask] |= (1 << iShift);
            }
          }
        }
        // See if any match
        for (int iLook = 0; iLook < nMarkRow - 1; iLook++) {
          unsigned int mask = bitMask[iLook];
          if (mask && mask != 0xffffffff && markRow) {
            int nMatch = 0;
            // count bits
            int nBits = 0;
            for (int i = 0; i < nMarked; i++) {
              if ((mask & 1) != 0)
                nBits++;
              mask = mask >> 1;
            }
            mask = bitMask[iLook];
            for (int jLook = iLook + 1; jLook < nMarkRow; jLook++) {
              if (bitMask[jLook] == mask)
                nMatch++;
            }
            if (nMatch) {
#if CBC_USEFUL_PRINTING > 0
              printf("For column %d %d matches on size %d\n",
                iColumn, nMatch, nBits);
              if (nMatch > 1) {
                printf("WHAT NOW\n");
              }
#endif
              // Pad out here
              int jColumn1 = -1, jColumn2 = -1;
              int iRow1 = whichRow[iLook];
              int iRow2 = -1;
              for (int jLook = iLook + 1; jLook < nMarkRow; jLook++) {
                if (bitMask[jLook] == mask) {
                  iRow2 = whichRow[jLook];
                  break;
                }
              }
              rowNumber[0] = iRow1;
              rowNumber[1] = iRow2;
              colLower[0] = 0.0;
              colUpper[0] = 1.0;
#if CBC_USEFUL_PRINTING > 0
              binary[0] = 1;
              int nCol = nBits;
#endif
              assert(rowLength[iRow1] <= 3);
              int nRow = 2;
              for (int j = 0; j < 2; j++) {
                int iRow = rowNumber[j];
                lower[j] = rowLower[iRow];
                upper[j] = rowUpper[iRow];
                for (CoinBigIndex k = rowStart[iRow];
                     k < rowStart[iRow] + rowLength[iRow]; k++) {
                  int kColumn = column[k];
                  if (kColumn == iColumn) {
                    el[0][j] = elementByRow[k];
                  } else if (kColumn == jColumn1) {
                    el[1][j] = elementByRow[k];
                  } else if (kColumn == jColumn2) {
                    el[2][j] = elementByRow[k];
                  } else if (jColumn1 < 0) {
                    jColumn1 = kColumn;
                    el[1][j] = elementByRow[k];
                    colLower[1] = columnLower[kColumn];
                    colUpper[1] = columnUpper[kColumn];
#if CBC_USEFUL_PRINTING > 0
                    binary[1] = 0;
                    if (!colLower[1] && colUpper[1] == 1.0 && newModel->isInteger(kColumn))
                      binary[1] = 1;
#endif
                  } else if (jColumn2 < 0) {
                    jColumn2 = kColumn;
                    el[2][j] = elementByRow[k];
                    colLower[2] = columnLower[kColumn];
                    colUpper[2] = columnUpper[kColumn];
#if CBC_USEFUL_PRINTING > 0
                    binary[2] = 0;
                    if (!colLower[2] && colUpper[2] == 1.0 && newModel->isInteger(kColumn))
                      binary[2] = 1;
#endif
                  } else {
                    abort();
                  }
                }
              }
              if (jColumn1 < 0) {
#if CBC_USEFUL_PRINTING > 0
                printf("**Why jColumn1 negative\n");
#endif
                continue;
              }
#if CBC_USEFUL_PRINTING > 0
              colNumber[0] = iColumn;
              colNumber[1] = jColumn1;
              colNumber[2] = jColumn2;
#endif
              // do something
              assert(columnLower[iColumn] == 0.0);
              assert(columnUpper[iColumn] == 1.0);
              if (nBits == 2) {
#if CBC_USEFUL_PRINTING > 0
                printf("iColumn %d %g %g binary\n", iColumn, colLower[0],
                  colUpper[0]);
                printf("jColumn1 %d %g %g %s\n", jColumn1, colLower[1],
                  colUpper[1], newModel->isInteger(jColumn1) ? "integer" : "continuous");
                for (int i = 0; i < 2; i++) {
                  printf("%d row %d %g <= %g*x0%s%g*x1 <= %g\n",
                    i, rowNumber[i], lower[i], el[0][i],
                    el[1][i] > 0.0 ? " +" : " -", fabs(el[1][i]), upper[i]);
                }
#endif
                double lower1[2], upper1[2];
                lower1[0] = lower1[1] = colLower[1];
                upper1[0] = upper1[1] = colUpper[1];
                double l, u;
                // binary at 0
                l = lower[0];
                u = upper[0];
                if (el[1][0] > 0.0) {
                  lower1[0] = std::max(lower1[0], l / el[1][0]);
                  upper1[0] = std::min(upper1[0], u / el[1][0]);
                } else {
                  lower1[0] = std::max(lower1[0], u / el[1][0]);
                  upper1[0] = std::min(upper1[0], l / el[1][0]);
                }
                l = lower[1];
                u = upper[1];
                if (el[1][1] > 0.0) {
                  lower1[0] = std::max(lower1[0], l / el[1][1]);
                  upper1[0] = std::min(upper1[0], u / el[1][1]);
                } else {
                  lower1[0] = std::max(lower1[0], u / el[1][1]);
                  upper1[0] = std::min(upper1[0], l / el[1][1]);
                }
                l = lower[0] - el[0][0];
                u = upper[0] - el[0][0];
                if (el[1][0] > 0.0) {
                  lower1[1] = std::max(lower1[1], l / el[1][0]);
                  upper1[1] = std::min(upper1[1], u / el[1][0]);
                } else {
                  lower1[1] = std::max(lower1[1], u / el[1][0]);
                  upper1[1] = std::min(upper1[1], l / el[1][0]);
                }
                l = lower[1] - el[0][1];
                u = upper[1] - el[0][1];
                if (el[1][1] > 0.0) {
                  lower1[1] = std::max(lower1[1], l / el[1][1]);
                  upper1[1] = std::min(upper1[1], u / el[1][1]);
                } else {
                  lower1[1] = std::max(lower1[1], u / el[1][1]);
                  upper1[1] = std::min(upper1[1], l / el[1][1]);
                }
                if (std::min(lower1[0], lower1[1]) > colLower[1] + 1.0e-6) {
#if CBC_USEFUL_PRINTING > 0
                  printf("for jColumn1 0-bounds %g,%g 1-bounds %g,%g\n",
                    lower1[0], upper1[0], lower1[1], upper1[1]);
#endif
                  double value = std::min(lower1[0], lower1[1]);
                  if (newModel->isInteger(jColumn1))
                    value = ceil(value - 1.0e-5);
#if CBC_USEFUL_PRINTING > 0
                  printf("increasing lb on %d from %g to %g\n",
                    jColumn1, colLower[1], value);
#endif
                  colLower[1] = value;
                  newModel->setColLower(jColumn1, value);
                }
                if (std::max(upper1[0], upper1[1]) < colUpper[1] - 1.0e-6) {
#if CBC_USEFUL_PRINTING > 0
                  printf("for jColumn1 0-bounds %g,%g 1-bounds %g,%g\n",
                    lower1[0], upper1[0], lower1[1], upper1[1]);
#endif
                  double value = std::max(upper1[0], upper1[1]);
                  if (newModel->isInteger(jColumn1))
                    value = floor(value + 1.0e-5);
#if CBC_USEFUL_PRINTING > 0
                  printf("decreasing ub on %d from %g to %g\n",
                    jColumn1, colUpper[1], value);
#endif
                  colUpper[1] = value;
                  newModel->setColUpper(jColumn1, value);
                }
                if (lower1[0] > colUpper[1] + 1.0e-6 || upper1[0] < colLower[1] - 1.0e-6) {
#if CBC_USEFUL_PRINTING > 0
                  printf("for jColumn1 0-bounds %g,%g 1-bounds %g,%g\n",
                    lower1[0], upper1[0], lower1[1], upper1[1]);
                  printf("fixing %d to 1\n", iColumn);
#endif
                  colLower[0] = 1.0;
                  newModel->setColLower(iColumn, 1.0);
                  nMarkRow = 0; // stop looking
                }
                if (lower1[1] > colUpper[1] + 1.0e-6 || upper1[1] < colLower[1] - 1.0e-6) {
#if CBC_USEFUL_PRINTING > 0
                  printf("for jColumn1 0-bounds %g,%g 1-bounds %g,%g\n",
                    lower1[0], upper1[0], lower1[1], upper1[1]);
                  printf("fixing %d to 0\n", iColumn);
#endif
                  colUpper[0] = 0.0;
                  newModel->setColLower(iColumn, 0.0);
                  nMarkRow = 0; // stop looking
                }
                if (colLower[0] > colUpper[0] + 1.0e-6 || colLower[1] > colUpper[1] + 1.0e-6) {
#if CBC_USEFUL_PRINTING > 0
                  printf("** infeasible\n");
#endif
                  feasible = false;
                }
              } else {
#if CBC_USEFUL_PRINTING > 0
                printf("iColumn %d %g %g binary\n", iColumn, colLower[0],
                  colUpper[0]);
                printf("jColumn1 %d %g %g %s\n", jColumn1, colLower[1],
                  colUpper[1], newModel->isInteger(jColumn1) ? "integer" : "continuous");
                printf("jColumn2 %d %g %g %s\n", jColumn2, colLower[2],
                  colUpper[2], newModel->isInteger(jColumn2) ? "integer" : "continuous");
                for (int i = 0; i < 2; i++) {
                  printf("%d row %d %g <= %g*x0%s%g*x1%s%g*x2 <= %g\n",
                    i, rowNumber[i], lower[i], el[0][i],
                    el[1][i] > 0.0 ? " +" : " -", fabs(el[1][i]),
                    el[2][i] > 0.0 ? " +" : " -", fabs(el[2][i]), upper[i]);
                }
#endif
                // Find other doubleton rows
                for (CoinBigIndex j1 = columnStart[jColumn1];
                     j1 < columnStart[jColumn1] + columnLength[jColumn1]; j1++) {
                  int iRow = row[j1];
                  if (rowLength[iRow] == 2 && rowTypeAll[iRow]) {
                    CoinBigIndex start = rowStart[iRow];
                    bool good = false;
                    if (jColumn1 == column[start] && jColumn2 == column[start + 1]) {
                      good = true;
                      el[1][nRow] = elementByRow[start];
                      el[2][nRow] = elementByRow[start + 1];
                    } else if (jColumn1 == column[start + 1] && jColumn2 == column[start]) {
                      good = true;
                      el[1][nRow] = elementByRow[start + 1];
                      el[2][nRow] = elementByRow[start];
                    }
                    if (good) {
                      rowNumber[nRow] = iRow;
                      lower[nRow] = rowLower[iRow];
                      upper[nRow] = rowUpper[iRow];
                      el[0][nRow] = 0.0;
#if CBC_USEFUL_PRINTING > 0
                      printf("%d row %d %g <= %g*x0%s%g*x1%s%g*x2 <= %g\n",
                        nRow, iRow, lower[nRow], el[0][nRow],
                        el[1][nRow] > 0.0 ? " +" : " -", fabs(el[1][nRow]),
                        el[2][nRow] > 0.0 ? " +" : " -", fabs(el[2][nRow]), upper[nRow]);
#endif
                      nRow++;
                      if (nRow == TRIPLE_ROWS)
                        break;
                    }
                  }
                }
              }
              // Now put in <= form
              int nLook = nRow;
              // types 1 <=, -1 >= 2 first ==, -2 second ==, 3 <= from range, -3 >= from range
              for (int i = 0; i < nRow; i++) {
                rowType[i] = 1;
                if (lower[i] < -1.0e20) {
                  rhs[i] = upper[i];
                } else if (upper[i] > 1.0e20) {
                  rowType[i] = -1;
                  rhs[i] = -lower[i];
                  for (int j = 0; j < nBits; j++)
                    el[j][i] = -el[j][i];
                } else if (lower[i] == upper[i]) {
                  rhs[i] = upper[i];
                  rowType[i] = 2;
                  rowType[nLook] = -2;
                  rowNumber[nLook] = rowNumber[i];
                  rhs[nLook] = -upper[i];
                  for (int j = 0; j < nBits; j++)
                    el[j][nLook] = -el[j][i];
                  nLook++;
                } else {
                  rhs[i] = upper[i];
                  rowType[i] = 3;
                  rowType[nLook] = -3;
                  rowNumber[nLook] = rowNumber[i];
                  rhs[nLook] = -lower[i];
                  for (int j = 0; j < nBits; j++)
                    el[j][nLook] = -el[j][i];
                  nLook++;
                }
              }
              // scale
              for (int i = 0; i < nRow; i++) {
                double largest = 0.0;
                double smallest = COIN_DBL_MAX;
                for (int k = 1; k < nBits; k++) {
                  double value = fabs(el[k][i]);
                  if (value) {
                    largest = std::max(largest, value);
                    smallest = std::min(smallest, value);
                  }
                }
                scaleFactor[i] = 1.0 / sqrt(largest * smallest);
              }
              // Look at all possible combinations
              // For now just look at first
              // also ignore bounds (should subtract out lbs from rhs)
              assert(nBits <= 3);
              memset(result, 0, sizeof(result));
              // bottom 4 bits for 0, next 4 for 1
              // 0 not stronger, 1 ==, 2 stronger (els)
              // 0 not stronger, 4 ==, 8 stronger (rhs)
              for (int j = 0; j < 2; j++) {
                for (int k = 0; k < nLook; k++)
                  modifiedRhs[k] = scaleFactor[k] * (rhs[k] - j * el[0][k]);
                for (int k1 = 0; k1 < nLook; k1++) {
                  for (int k2 = 0; k2 < nLook; k2++) {
                    if (k1 != k2 && rowNumber[k1] != rowNumber[k2]) {
                      int stronger = 8;
                      if (modifiedRhs[k1] > modifiedRhs[k2] + 1.0e-7) {
                        stronger = 0;
                      } else if (modifiedRhs[k1] > modifiedRhs[k2] - 1.0e-7) {
                        stronger = 4;
                      }
                      if (stronger) {
                        int strongerEl = 2;
                        for (int k3 = 1; k3 < nBits; k3++) {
                          if (scaleFactor[k1] * el[k3][k1] < scaleFactor[k2] * el[k3][k2] - 1.0e-9) {
                            stronger = 0;
                            break;
                          } else if (scaleFactor[k1] * el[k3][k1] < scaleFactor[k2] * el[k3][k2] + 1.0e-9) {
                            strongerEl = 1;
                          }
                        }
                        if (stronger) {
                          stronger |= strongerEl;
                        }
                        result[k1][k2] |= static_cast< unsigned char >(stronger << (4 * j));
                      }
                    }
                  }
                }
              }
              int dropped = -1;
              for (int k1 = 0; k1 < nLook; k1++) {
                for (int k2 = 0; k2 < nLook; k2++) {
                  if (k1 != k2 && rowNumber[k1] != rowNumber[k2]) {
                    int state0 = result[k1][k2] & 15;
                    int state1 = result[k1][k2] >> 4;
                    if (state0 && state1) {
                      if (state0 == 5) {
                        if (state1 == 5) {
                          // same
                          if (abs(rowType[k1]) == 1)
                            dropped = k1;
                          else
                            dropped = k2;
#if CBC_USEFUL_PRINTING > 0
                          printf("ZZZsame ");
#endif
                        } else if (state1 == 9) {
                          // drop second
                          dropped = k2;
#if CBC_USEFUL_PRINTING > 0
                          printf("ZZZfirst ");
#endif
                        } else {
#if CBC_USEFUL_PRINTING > 0
                          printf("ZZYY ");
#endif
                        }
                      } else if (state0 == 9) {
                        // drop second
                        dropped = k2;
#if CBC_USEFUL_PRINTING > 0
                        printf("ZZZsecond ");
#endif
                      } else {
#if CBC_USEFUL_PRINTING > 0
                        printf("ZZYY ");
#endif
                      }
#if CBC_USEFUL_PRINTING > 0
                      printf("row %d (%d) and row %d (%d) status at 0 is %d status at 1 is %d\n",
                        k1, rowNumber[k1], k2, rowNumber[k2],
                        state0, state1);
#endif
                      if (dropped >= 0)
                        break;
                    }
                  }
                }
                if (dropped >= 0)
                  break;
              }
              if (dropped >= 0) {
                int iRow = rowNumber[dropped];
                // Think if ranged or ==
                if (rowLower[iRow] < -1.0e30 || rowUpper[iRow] > 1.0e30) {
                  newModel->setRowLower(iRow, -COIN_DBL_MAX);
                  newModel->setRowUpper(iRow, COIN_DBL_MAX);
#if CBC_USEFUL_PRINTING > 0
                  nFreed++;
#endif
                } else {
#if CBC_USEFUL_PRINTING > 0
                  printf("XXXYYY - shouldn't drop ranged/equality row????\n");
#endif
                }
              }
              // stop rest
              for (int jLook = iLook; jLook < nMarkRow; jLook++) {
                if (bitMask[jLook] == mask)
                  bitMask[jLook] = 0;
              }
            }
          }
        }
        for (int i = 0; i < nMarked; i++) {
          marked[which[i]] = 0;
        }
        for (int i = 0; i < nMarkRow; i++) {
          markRow[whichRow[i]] = 0;
        }
      }
    }
#if CBC_USEFUL_PRINTING > 0
    if (nFreed)
      printf("%d rows freed up\n", nFreed);
#endif
    delete[] which;
  }
#endif
#if 0
  // Do domination stuff
  if (iBigPass==0) {
    // Row copy
    CoinPackedMatrix matrixByRow(*newModel->getMatrixByRow());
    const double * elementByRow = matrixByRow.getElements();
    const int * column = matrixByRow.getIndices();
    const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
    const int * rowLength = matrixByRow.getVectorLengths();

    // Column copy
    CoinPackedMatrix  matrixByCol(*newModel->getMatrixByCol());
    //const double * element = matrixByCol.getElements();
    const int * row = matrixByCol.getIndices();
    const CoinBigIndex * columnStart = matrixByCol.getVectorStarts();
    const int * columnLength = matrixByCol.getVectorLengths();

    const double * rowLower = newModel->getRowLower();
    const double * rowUpper = newModel->getRowUpper();
    const double * columnLower = newModel->getColLower();
    const double * columnUpper = newModel->getColUpper();
    // get sizes for canonical form (overestimate if free columns)
    int nRows=numberRows;
    CoinBigIndex nEls=matrixByRow.getNumElements();
    for (int iRow=0;iRow<numberRows;iRow++) {
      if (rowLower[iRow]>-1.0e30&&rowUpper[iRow]<1.0e30) {
	nRows++;
	nEls += rowLength[iRow];
      }
    }
    int * rowNumber = new int[3*nRows+numberColumns];
    int * rowBinary = rowNumber+nRows;
    int * rowPos = rowBinary+nRows;
    int * whichColumn = rowPos+nRows;
    double * elementByRow2 = new double [nEls+numberColumns+nRows];
    double * columnValue = elementByRow2+nEls;
    double * rhs = columnValue+numberColumns;
    int * column2 = new int [nEls];
    CoinBigIndex * rowStart2 = new CoinBigIndex[nRows+1];
    char * marked = new char [numberColumns+nRows];
    char * markedRow = marked + numberColumns;
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      columnValue[iColumn]=0.0;
      if (columnLower[iColumn]<-1.0e10&&columnUpper[iColumn]>1.0e10) {
	marked[iColumn]=-1;
      } else if (fabs(columnUpper[iColumn])<fabs(columnLower[iColumn])) {
	// flip
	marked[iColumn]=2;
	columnValue[iColumn]=-columnUpper[iColumn];
	if (newModel->isInteger(iColumn) &&
	    columnUpper[iColumn]==columnLower[iColumn]+1) 
	  marked[iColumn]=3;
      } else {
	marked[iColumn]=0;
	columnValue[iColumn]=columnLower[iColumn];
	if (newModel->isInteger(iColumn) &&
	    columnUpper[iColumn]==columnLower[iColumn]+1) 
	  marked[iColumn]=1;
      }
    }
    nRows=0;
    nEls=0;
    rowStart2[0]=0;
    for (int iRow=0;iRow<numberRows;iRow++) {
      CoinBigIndex start = rowStart[iRow];
      CoinBigIndex end = start + rowLength[iRow];
      for (int iTry=0;iTry<2;iTry++) {
	double multiplier;
	double rhsValue;
	if (!iTry) {
	  multiplier=1.0;
	  rhsValue = rowUpper[iRow];
	} else {
	  multiplier=-1.0;
	  rhsValue = -rowLower[iRow];
	}
	if (rhsValue<1.0e30) {
	  char typeRow=iTry;
	  int nPos=0;
	  int nInt=0;
	  double largest=0.0;
	  double smallest=COIN_DBL_MAX;
	  for (CoinBigIndex k=start;k<end;k++) {
	    int kColumn = column[k];
	    int type = marked[kColumn];
	    double value = multiplier*elementByRow[k];
	    if (type<0) {
	      nEls=rowStart2[nRows];
	      typeRow=-1;;
	      break;
	    } else if ((type&2)!=0) {
	      value = -value;
	    }
	    if ((type&1)!=0)
	      nInt++;
	    rhsValue -= value*columnValue[kColumn];
	    elementByRow2[nEls]=value;
	    if (value>0.0)
	      nPos++;
	    largest=std::max(fabs(value),largest);
	    smallest=std::min(fabs(value),smallest);
	    column2[nEls++]=kColumn;
	  }
	  if (typeRow>=0 && smallest*1.0e7>largest) {
	    double scale = sqrt(largest*smallest);
	    if (fabs(rhsValue)>1.0e6*scale||
		(rhsValue&&fabs(rhsValue)<1.0e-6*scale)) {
	      scale=0.0;
	    } else if (rhsValue) {
	      scale=1.0/fabs(rhsValue);
	    }
	    if (scale) {
	      rhs[nRows]=scale*rhsValue;
	      for (CoinBigIndex k=rowStart2[nRows];k<nEls;k++) 
		elementByRow2[k] *= scale;
	      rowPos[nRows]=nPos;
	      markedRow[nRows]=typeRow;
	      rowBinary[nRows]=nInt;
	      rowNumber[nRows++]=iRow;
	      rowStart2[nRows]=nEls;
	    } else {
	      nEls=rowStart2[nRows];
	    }
	  }
	}
      }
    }
    memset(columnValue,0,numberColumns*sizeof(double));
    double tolerance = 1.0e-9;
    for (int iRow=0;iRow<nRows;iRow++) {
      CoinBigIndex start = rowStart2[iRow];
      CoinBigIndex end = rowStart2[iRow+1];
      int n=0;
      int nInt=rowBinary[iRow];
      for (CoinBigIndex k=start;k<end;k++) {
	int kColumn = column2[k];
	double value = elementByRow2[k];
	columnValue[kColumn]=value;
	whichColumn[n++]=kColumn;
      }
      double rhsValue=rhs[iRow];
      int nPos=rowPos[iRow];
      int nNeg=n-nPos;
      // initially only short integer rows
      if (n>3)
	nInt=0;
      // for first try ignore integers!
      nInt=0;
      if (nInt) {
      } else {
	for (int jRow=iRow+1;jRow<nRows;jRow++) {
	  CoinBigIndex start2 = rowStart2[jRow];
	  CoinBigIndex end2 = rowStart2[jRow+1];
	  int n2=end2-start2;
	  int nPos2=rowPos[jRow];
	  int nNeg2=n2-nPos2;
	  int nInt2=rowBinary[jRow];
	  double rhsValue2=rhs[jRow];
	  // initially only short integer rows
	  if (n2>3)
	    nInt2=0;
	  // for first try ignore integers!
	  nInt2=0;
	  if (nInt2) {
	  } else {
	    // continuous tests
	    // -1 iRow may be stronger, +1 jRow may be stronger, 0 continue, 2 == els
	    int way=2;
	    if (rhsValue>rhsValue2+tolerance)
	      way=1;
	    else if (rhsValue2>rhsValue+tolerance)
	      way=-1;
	    if (nNeg2>nNeg) {
	      // iRow can be stronger
	      if (nPos2>nPos || way == 1) 
		way=0;
	      else
		way=-1;
	    } else if (nNeg2==nNeg) {
	      // iRow can be either way
	      if (nPos2>nPos) {
		if (way!=-1)
		  way=1;
		else
		  way=0;
	      } else if (nPos2<nPos) { 
		if (way!=1)
		  way=-1;
		else
		  way=0;
	      }
	    } else {
	      // jRow can be stronger
	      if (nPos2<nPos || way==-1) 
		way=0;
	      else
		way=1;
	    }
	    int nHitPos=0;
	    int nHitNeg=0;
	    if (way==-1) {
	      // iRow may be stronger
	      for (CoinBigIndex k2=start2;k2<end2;k2++) {
		int kColumn = column2[k2];
		double value = elementByRow2[k2];
		double valueI=columnValue[kColumn];
		if (value>valueI+1.0e-12) {
		  way=0;
		  break;
		} else {
		  if (valueI<0.0)
		    nHitNeg++;
		}
	      }
	      if (nHitNeg<nNeg2)
		way=0;
	    } else if (way==1) {
	      // jRow may be stronger
	      for (CoinBigIndex k2=start2;k2<end2;k2++) {
		int kColumn = column2[k2];
		double value = elementByRow2[k2];
		double valueI=columnValue[kColumn];
		if (value<valueI-1.0e-12) {
		  way=0;
		  break;
		} else {
		  if (valueI>0.0)
		    nHitPos++;
		}
	      }
	      if (nHitPos<nPos)
		way=0;
	    } else if (way==2) {
	      // same number and rhs - could go either way
	      CoinBigIndex k2;
	      for (k2=start2;k2<end2;k2++) {
		int kColumn = column2[k2];
		double value = elementByRow2[k2];
		if (value<columnValue[kColumn]-1.0e-12) {
		  way=-1;
		  break;
		} else if (value>columnValue[kColumn]+1.0e-12) {
		  way=1;
		  break;
		}
	      }
	      k2++;
	      if (way==1) {
		for (;k2<end2;k2++) {
		  int kColumn = column2[k2];
		  double value = elementByRow2[k2];
		  double valueI=columnValue[kColumn];
		  if (value<valueI-1.0e-12) {
		    way=0;
		    break;
		  } else {
		    if (valueI>0.0)
		      nHitPos++;
		  }
		}
		if (nHitPos<nPos)
		  way=0;
	      } else if (way==-1) {
		for (;k2<end2;k2++) {
		  int kColumn = column2[k2];
		  double value = elementByRow2[k2];
		  double valueI=columnValue[kColumn];
		  if (value>valueI+1.0e-12) {
		    way=0;
		    break;
		  } else {
		    if (valueI<0.0)
		      nHitNeg++;
		  }
		}
		if (nHitNeg<nNeg2)
		  way=0;
	      }
	    }
#if CBC_USEFUL_PRINTING
	    if (way) {
	      int iRowX=rowNumber[iRow];
	      int jRowX=rowNumber[jRow];
	      CoinBigIndex startI = rowStart[iRowX];
	      CoinBigIndex endI = startI + rowLength[iRowX];
	      CoinBigIndex startJ = rowStart[jRowX];
	      CoinBigIndex endJ = startJ + rowLength[jRowX];
	      printf("way %d for row %d (%d - %d els) and %d (%d - %d els)\n",
		     way,iRow,iRowX,endI-startI,jRow,jRowX,endJ-startJ);
	      printf("%g <= ",rowLower[iRowX]);
	      if (endI-startI<100&&endJ-startJ<10) {
		for (CoinBigIndex k=startI;k<endI;k++) 
		  printf("(%d,%g) ",column[k],elementByRow[k]);
	      } else {
		printf("something ");
	      }
	      printf("<= %g\n",rowUpper[iRowX]);
	      printf("%g <= ",rowLower[jRowX]);
	      if (endI-startI<100&&endJ-startJ<10) {
		for (CoinBigIndex k=startJ;k<endJ;k++) 
		  printf("(%d,%g) ",column[k],elementByRow[k]);
	      } else {
		printf("something ");
	      }
	      printf("<= %g\n",rowUpper[jRowX]);
	    }
#endif
	  }
	}
      }
      for (int j=0;j<n;j++) {
	int kColumn = whichColumn[j];
	columnValue[kColumn]=0.0;
      }
    }
    delete [] rowNumber;
    delete [] elementByRow2;
    delete [] column2;
    delete [] rowStart2;
    delete [] marked;
  }
#endif
  bool noStrengthening = false;
  for (int iPass = 0; iPass < numberPasses; iPass++) {
    // Statistics
    int numberFixed = 0;
    int numberTwo = twoCuts.sizeRowCuts();
    int numberStrengthened = 0;
    info.pass = iPass;
    info.options = 0;
    int numberChangedThisPass = 0;
#if 1
    // look at cliques every time
    if ((options_ & 32) != 0) {
      OsiSolverInterface *temp = cliqueIt(*newModel, 0.0001);
      if (temp) {
#if CBC_USEFUL_PRINTING
        printf("bigpass %d pass %d after cliques %d rows, before %d\n",
          iBigPass, iPass, temp->getNumRows(), newModel->getNumRows());
#endif
        if (temp->getNumRows() < newModel->getNumRows()) {
          numberChangedThisPass += newModel->getNumRows() - temp->getNumRows();
          delete newModel;
          newModel = temp;
          numberRows = newModel->getNumRows();
        } else {
          delete temp;
        }
      } else {
#if CBC_USEFUL_PRINTING
        printf("bigpass %d pass %d no more cliques %d rows\n",
          iBigPass, iPass, newModel->getNumRows());
#endif
      }
    }
#endif
    /*
      needResolve    solution is stale
      rebuilt   constraint system deleted and recreated (implies initialSolve)
    */
    for (int iGenerator = firstGenerator; iGenerator < lastGenerator; iGenerator++) {
      bool needResolve = false;
      bool rebuilt = false;
      OsiCuts cs;
      CoinZeroN(whichCut, numberRows);
      CglProbing *probingCut = NULL;
      int numberFromCglDuplicate = 0;
      const int *duplicate = NULL;
      CglDuplicateRow *dupRow = NULL;
      CglClique *cliqueGen = NULL;
      CglBKClique *bkCliqueGen = NULL;
      if (iGenerator >= 0) {
        //char name[20];
        //sprintf(name,"prex%2.2d.mps",iGenerator);
        //newModel->writeMpsNative(name, NULL, NULL,0,1,0);
        // refresh as model may have changed
        generator_[iGenerator]->refreshSolver(newModel);
        // skip duplicate rows except once
        dupRow = dynamic_cast< CglDuplicateRow * >(generator_[iGenerator]);
        cliqueGen = dynamic_cast< CglClique * >(generator_[iGenerator]);
        bkCliqueGen = dynamic_cast< CglBKClique * >(generator_[iGenerator]);
        if ((cliqueGen || bkCliqueGen) && iPass)
          continue;
        if (dupRow && (iPass || iBigPass))
          continue;
        probingCut = dynamic_cast< CglProbing * >(generator_[iGenerator]);
#if CBC_USEFUL_PRINTING > 0
        double time1 = CoinCpuTime();
#endif
        double inspectTime1 = inspect_ ? CoinWallclockTime() : 0.0;
        if (inspect_) {
          char buf[256];
          FILE *fp = handler_->filePointer();
          if (fp) {
            if (probingCut) {
              snprintf(buf, sizeof(buf),
                "[Preproc pass %d.%d] Probing: maxProbe=%d, maxPass=%d, maxLook=%d, maxElem=%d (minor pass %d/%d)",
                iBigPass, iPass,
                probingCut->getMaxProbeRoot(), probingCut->getMaxPassRoot(),
                probingCut->getMaxLookRoot(), probingCut->getMaxElementsRoot(),
                iPass + 1, numberPasses);
            } else if (cliqueGen || bkCliqueGen) {
              snprintf(buf, sizeof(buf),
                "[Preproc pass %d.%d] Clique separator (minor pass %d/%d)",
                iBigPass, iPass, iPass + 1, numberPasses);
            } else if (dupRow) {
              snprintf(buf, sizeof(buf),
                "[Preproc pass %d.%d] DuplicateRow (minor pass %d/%d)",
                iBigPass, iPass, iPass + 1, numberPasses);
            } else {
              snprintf(buf, sizeof(buf),
                "[Preproc pass %d.%d] Generator %d (minor pass %d/%d)",
                iBigPass, iPass, iGenerator, iPass + 1, numberPasses);
            }
            fprintf(fp, "  %s\n", buf);
            fflush(fp);
          }
        }
        if (!probingCut) {
          generator_[iGenerator]->generateCuts(*newModel, cs, info);
        } else {
          info.options = 64 | 2048;
          probingCut->setMode(useSolution ? 4 : 4 | 64);
          int saveMaxElements = probingCut->getMaxElementsRoot();
          int saveMaxProbe = probingCut->getMaxProbeRoot();
          int saveMaxLook = probingCut->getMaxLookRoot();
	  if ((!iBigPass||(options_&64)!=0)&&!iPass&&(options_&(16|64))!=0) {
	    //if (/*!iBigPass &&*/ !iPass /*&&(options_&(16|64))!=0*/) {
            noStrengthening = true;
            numberPasses = 1;
            probingCut->setMaxProbeRoot(std::max(saveMaxProbe, 1000));
            probingCut->setMaxElementsRoot(std::max(saveMaxElements, 2000));
	    int maxLook = std::min(numberColumns, numberRows)/2;
	    if ((options_&16)!=0) {
	      maxLook = numberColumns;
	      info.options |= 32768;
	    }
	    maxLook = std::min(maxLook,2000);
            probingCut->setMaxLookRoot(std::max(saveMaxLook, maxLook));
            options_ &= ~16;
          } else if (iPass || (options_ & 64) == 0) {
            // cut back
            probingCut->setMaxElementsRoot(probingCut->getMaxElements());
            probingCut->setMaxProbeRoot(probingCut->getMaxProbe());
            probingCut->setMaxLookRoot(probingCut->getMaxLook());
          }
          if (inspect_) {
            char buf[256];
            FILE *fp = handler_->filePointer();
            if (fp) {
              snprintf(buf, sizeof(buf),
                "[Preproc pass %d.%d] Probing actual: maxProbe=%d, maxPass=%d, maxLook=%d, maxElem=%d  cols=%d, rows=%d",
                iBigPass, iPass,
                probingCut->getMaxProbeRoot(), probingCut->getMaxPassRoot(),
                probingCut->getMaxLookRoot(), probingCut->getMaxElementsRoot(),
                numberColumns, numberRows);
              fprintf(fp, "  %s\n", buf);
              fflush(fp);
            }
          }
          probingCut->generateCutsAndModify(*newModel, cs, &info);
          probingCut->setMaxElementsRoot(saveMaxElements);
          probingCut->setMaxProbeRoot(saveMaxProbe);
          probingCut->setMaxLookRoot(saveMaxLook);
          if (!iPass && (!cs.sizeColCuts() || iBigPass > 2))
            options_ &= ~64; // switch off heavy
        }
#if CBC_USEFUL_PRINTING > 0
        printf("Generator %d took %g seconds\n",
          iGenerator, CoinCpuTime() - time1);
        printf("After probing1 %d row cuts and %d column cuts\n",
          cs.sizeRowCuts(), cs.sizeColCuts());
#endif
        if (inspect_) {
          FILE *fp = handler_->filePointer();
          if (fp) {
            double elapsed = CoinWallclockTime() - inspectTime1;
            int nR = newModel->getNumRows();
            int nC = newModel->getNumCols();
            int nNZ = newModel->getNumElements();
            // cut pool NZ stats
            int nRowCuts = cs.sizeRowCuts();
            int cutNZtotal = 0, cutNZmax = 0;
            for (int ci = 0; ci < nRowCuts; ci++) {
              int n = cs.rowCut(ci).row().getNumElements();
              cutNZtotal += n;
              if (n > cutNZmax) cutNZmax = n;
            }
            double cutNZavg = nRowCuts > 0 ? (double)cutNZtotal / nRowCuts : 0.0;
            if (nRowCuts > 0)
              fprintf(fp, "  [Preproc pass %d.%d] Generator done: %.2fs  "
                "row cuts=%d (NZ: total=%d avg=%.1f max=%d)  col cuts=%d  "
                "model now: %d rows, %d cols, %d NZ\n",
                iBigPass, iPass, elapsed,
                nRowCuts, cutNZtotal, cutNZavg, cutNZmax,
                cs.sizeColCuts(), nR, nC, nNZ);
            else
              fprintf(fp, "  [Preproc pass %d.%d] Generator done: %.2fs  "
                "row cuts=0  col cuts=%d  model now: %d rows, %d cols, %d NZ\n",
                iBigPass, iPass, elapsed, cs.sizeColCuts(), nR, nC, nNZ);
            fflush(fp);
          }
        }
        if (cs.sizeColCuts() && iPass < numberPasses - 100 && !iBigPass) {
          // delete all row cuts for now????
          int n = cs.sizeRowCuts();
          for (int i = n - 1; i >= 0; i--)
            cs.eraseRowCut(i);
          for (int i = 0; i < numberRows; i++)
            whichCut[i] = 0;
        }
#if 1 //def CLIQUE_ANALYSIS
        if (probingCut) {
          //printf("ordinary probing\n");
          info.analyze(*newModel);
        }
#endif
        // If CglDuplicate may give us useless rows
        if (dupRow) {
          numberFromCglDuplicate = dupRow->numberOriginalRows();
          duplicate = dupRow->duplicate();
          if (cs.sizeRowCuts()) {
            // add to twoCuts (may be more, but ....)
            int numberRowCuts = cs.sizeRowCuts();
            for (int k = 0; k < numberRowCuts; k++) {
              OsiRowCut *thisCut = cs.rowCutPtr(k);
              twoCuts.insert(*thisCut);
            }
          }
        }
        if (cliqueGen && cs.sizeRowCuts()) {
          int n = cs.sizeRowCuts();
#if CBC_USEFUL_PRINTING
          printf("%d clique cuts\n", n);
#endif
          OsiSolverInterface *copySolver = newModel->clone();
          numberRows = copySolver->getNumRows();
          copySolver->applyCuts(cs);
          //static int kk=0;
          //char name[20];
          //kk++;
          //sprintf(name,"matrix%d",kk);
          //printf("writing matrix %s\n",name);
          //copySolver->writeMps(name);
          CglDuplicateRow dupCuts(copySolver);
          dupCuts.setMode(8);
          OsiCuts cs2;
          dupCuts.generateCuts(*copySolver, cs2, info);
#if CBC_USEFUL_PRINTING > 0
          printf("After probing dupCuts %d row cuts and %d column cuts\n",
            cs2.sizeRowCuts(), cs2.sizeColCuts());
#endif
          // -1 not used, -2 delete, -3 not clique
          const int *duplicate = dupCuts.duplicate();
          // -1 not used, >=0 earliest row affected
          const int *used = duplicate + numberRows + n;
          int numberDrop = 0;
          int *drop = new int[numberRows];
          for (int iRow = 0; iRow < numberRows; iRow++) {
            if (duplicate[iRow] == -2)
              drop[numberDrop++] = iRow;
          }
          int nOther = 0;
          for (int iRow = numberRows + n - 1; iRow >= numberRows; iRow--) {
#if 1
            int earliest = used[iRow];
            while (earliest >= numberRows) {
              if (duplicate[earliest] == -2)
                earliest = used[earliest];
              else
                break;
            }
#else
            int earliest = 0;
#endif
            if (duplicate[iRow] == -2 || earliest == -1 || earliest >= numberRows) {
              cs.eraseRowCut(iRow - numberRows);
              nOther++;
            }
          }
          n -= nOther;
          int newNumberRows = numberRows - numberDrop + n;
          bool special = (cliqueGen->getMinViolation() == -3.0);
#if CBC_USEFUL_PRINTING
          printf("could drop %d rows - current nrows %d other %d - new nrows %d\n",
            numberDrop, numberRows, nOther, newNumberRows);
#endif
          if (n <= numberDrop || special) {
#if CBC_USEFUL_PRINTING
            printf("Dropping rows current nrows %d - new nrows %d\n",
              numberRows, newNumberRows);
#endif
            if (newNumberRows > numberRows) {
              // need new array
              delete[] whichCut;
              whichCut = new OsiRowCut *[newNumberRows + 1];
              CoinZeroN(whichCut, newNumberRows);
              info.strengthenRow = whichCut;
            }
            newModel->deleteRows(numberDrop, drop);
            // may be able to delete some added cliques
            newModel->applyCuts(cs);
            numberRows = newModel->getNumRows();
	    {
	      bool saveTakeHint;
	      OsiHintStrength saveStrength;
	      newModel->getHintParam(OsiDoDualInResolve,
				     saveTakeHint, saveStrength);
	      newModel->setHintParam(OsiDoDualInResolve, solveWithDual, OsiHintTry);
	      double inspectLpStartM = inspect_ ? CoinWallclockTime() : 0.0;
	      newModel->resolve();
	      if (inspect_) {
		FILE *fp = handler_->filePointer();
		if (fp) {
		  const char *stat = newModel->isProvenOptimal() ? "optimal"
		    : newModel->isProvenPrimalInfeasible() ? "infeasible" : "not optimal";
		  fprintf(fp, "  [Preproc LP pass %d.%d] resolve (row drop): %.2fs  iters=%d  obj=%.6g  (%s)\n",
		    iBigPass, iPass,
		    CoinWallclockTime() - inspectLpStartM,
		    newModel->getIterationCount(),
		    newModel->getObjValue(), stat);
		  fflush(fp);
		}
	      }
	      newModel->setHintParam(OsiDoDualInResolve, saveTakeHint, saveStrength);
	      solveWithDual = true;
	    }
#if 0
	    int numberRows2=copySolver->getNumRows();
	    const double * rowLower = copySolver->getRowLower();
	    const double * rowUpper = copySolver->getRowUpper();
	    const CoinPackedMatrix * matrixByRow = copySolver->getMatrixByRow();
	    // Row copy
	    const double * elementByRow = matrixByRow->getElements();
	    const int * column = matrixByRow->getIndices();
	    const CoinBigIndex * rowStart = matrixByRow->getVectorStarts();
	    const int * rowLength = matrixByRow->getVectorLengths();
	    const double * solution = newModel->getColSolution();
	    for (int iRow=0;iRow<numberRows2;iRow++) {
	      double sum=0.0;
	      for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
		int iColumn = column[j];
		double value = elementByRow[j];
		sum += value*solution[iColumn];
	      }
	      assert (sum>rowLower[iRow]-1.0e-4&&sum<rowUpper[iRow]+1.0e-4);
	    }
#endif
          }
          delete copySolver;
          delete[] drop;
          continue;
          //for (int i=0;i<n;i++) {
          //OsiRowCut & thisCut = cs.rowCut(i);
          //thisCut.print();
          //}
        }
      } else {
#ifdef HEAVY_PROBING
        // special probing
        CglProbing generator1;
        probingCut = &generator1;
        generator1.setUsingObjective(false);
        generator1.setMaxPass(1);
        generator1.setMaxPassRoot(1);
        generator1.setMaxProbeRoot(100);
        generator1.setMaxLook(100);
        generator1.setRowCuts(3);
        if (heavyProbing) {
          generator1.setMaxElements(400);
          //generator1.setMaxLook(10000);
          generator1.setMaxProbeRoot(model->getNumCols());
        }
        // out for now - think about cliques
        if (!generator1.snapshot(*newModel, NULL, false)) {
          generator1.createCliques(*newModel, 2, 1000);
          // To get special stuff
          info.pass = 4;
          CoinZeroN(whichCut, numberRows);
          generator1.setMode(16 + 4);
          generator1.generateCutsAndModify(*newModel, cs, &info);
#if CBC_USEFUL_PRINTING > 0
          printf("After probing clique stuff %d row cuts and %d column cuts\n",
            cs.sizeRowCuts(), cs.sizeColCuts());
#endif
          // can we extend cliques?
          // make fake model
          OsiSolverInterface *fakeModel = generator1.cliqueModel(newModel, 1);
          // if above added rows then take out duplicates
          OsiSolverInterface *fakeModel2 = cliqueIt(*fakeModel, 0.0);
          delete fakeModel;
          //delete fakeModel2;
          delete newModel;
          newModel = fakeModel2;
#ifdef CLIQUE_ANALYSIS
          printf("special probing\n");
          info.analyze(*newModel);
#endif
        } else {
          feasible = false;
        }
#endif
      }
      // check changes
      // first are any rows strengthened by cuts
      int iRow;
#ifdef MAX_ADD_ELEMENTS_PREPROCESS
      const CoinPackedMatrix *tempRowCopy = newModel->getMatrixByRow();
      const int *tempRowLength = tempRowCopy->getVectorLengths();
#endif
      for (iRow = 0; iRow < numberRows; iRow++) {
        if (whichCut[iRow]) {
#ifdef MAX_ADD_ELEMENTS_PREPROCESS
          OsiRowCut *thisCut = whichCut[iRow];
          CoinPackedVector row = thisCut->row();
          if (row.getNumElements() <= tempRowLength[iRow]
              + MAX_ADD_ELEMENTS_PREPROCESS) {
            numberStrengthened++;
          } else {
            delete thisCut;
            whichCut[iRow] = NULL;
          }
#else
          numberStrengthened++;
#endif
        }
      }
      // Also can we get rid of duplicate rows
      int numberDrop = 0;
      for (iRow = 0; iRow < numberFromCglDuplicate; iRow++) {
        if (duplicate[iRow] == -2 || duplicate[iRow] >= 0) {
          numberDrop++;
          newModel->setRowBounds(iRow, -COIN_DBL_MAX, COIN_DBL_MAX);
        }
      }
      const double *columnLower = newModel->getColLower();
      const double *columnUpper = newModel->getColUpper();
      if ((numberStrengthened || numberDrop) && feasible) {
        /*
	  
	Deleting all rows and rebuilding invalidates everything, initialSolve will
	be required.
	*/
        if (useSolution)
          needResolve = true;
        rebuilt = true;
        // Easier to recreate entire matrix
        const CoinPackedMatrix *rowCopy = newModel->getMatrixByRow();
        const int *column = rowCopy->getIndices();
        const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
        const int *rowLength = rowCopy->getVectorLengths();
        const double *rowElements = rowCopy->getElements();
        const double *rowLower = newModel->getRowLower();
        const double *rowUpper = newModel->getRowUpper();
        CoinBuild build;
        // For basis
        char *keepRow = new char[numberRows];
        for (iRow = 0; iRow < numberRows; iRow++) {
          keepRow[iRow] = 0;
          OsiRowCut *thisCut = whichCut[iRow];
          //whichCut[iRow]=NULL;
          if (rowLower[iRow] > -1.0e20 || rowUpper[iRow] < 1.0e20) {
#if 0
	    if (thisCut) {
	      double * allColumns = new double[numberColumns];
	      int which[]={0,8,11,19,21,29,30,38,42,61,77,90,104,105,7,37};
	      memset(allColumns,0,numberColumns*sizeof(double));
	      for (int k=0;k<sizeof(which)/sizeof(int);k++) {
		allColumns[which[k]]=1.0;
	      }
	      double lb = thisCut->lb();
	      double ub = thisCut->ub();
	      CoinPackedVector row = thisCut->row();
	      printf("Cut on row %d - %g <= ",iRow,lb);
	      bool feas1=true;
	      double sum1=0.0;
	      for (int k = 0; k < row.getNumElements(); ++k) {
		int j = row.getIndices()[k];
		double value = row.getElements()[k];
		printf("(%d,%g) ",j,value);
		sum1 += value*allColumns[j];
	      }
	      if (sum1<lb-1.0e-3||sum1>ub+1.0e-3) {
		printf(" ******** ");
		feas1 = false;
	      }
	      printf("<= %g\n",ub);
	      printf("Old row %g <= ",rowLower[iRow]);
	      bool feas2=true;
	      double sum2=0.0;
              int start=rowStart[iRow];
	      int end = start + rowLength[iRow];
	      for (int k = start; k < end; ++k) {
		CoinBigIndex j = column[k];
		double value = rowElements[k];
		printf("(%d,%g) ",j,value);
		sum2 += value*allColumns[j];
	      }
	      if (sum2<rowLower[iRow]-1.0e-3||sum2>rowUpper[iRow]+1.0e-3) {
		printf(" ******** ");
		feas2 = false;
	      }
	      printf("<= %g\n",rowUpper[iRow]);
	      if (feas1 && !feas2)
		abort();
	      delete [] allColumns;
	    }
#endif
            if (!thisCut) {
              // put in old row
              CoinBigIndex start = rowStart[iRow];
              int kInt = -1;
              double newValue = 0.0;
              if (!iPass && !iBigPass) {
                // worthwhile seeing if odd gcd
                CoinBigIndex end = start + rowLength[iRow];
                double rhsAdjustment = 0.0;
                int nPosInt = 0;
                int nNegInt = 0;
                // Find largest integer coefficient
                CoinBigIndex k;
                for (k = start; k < end; ++k) {
                  int j = column[k];
                  if (columnUpper[j] > columnLower[j]) {
                    if (newModel->isInteger(j)) {
                      if (rowElements[k] > 0.0)
                        nPosInt++;
                      else
                        nNegInt++;
                    } else {
                      break; // no good
                    }
                  } else {
                    rhsAdjustment += columnLower[j] * rowElements[k];
                  }
                }
                if (k == end) {
                  // see if singleton coefficient can be strengthened
                  if ((nPosInt == 1 && nNegInt > 1) || (nNegInt == 1 && nPosInt > 1)) {
                    double lo;
                    double up;
                    if (rowLower[iRow] > -1.0e20)
                      lo = rowLower[iRow] - rhsAdjustment;
                    else
                      lo = -COIN_DBL_MAX;
                    if (rowUpper[iRow] < 1.0e20)
                      up = rowUpper[iRow] - rhsAdjustment;
                    else
                      up = COIN_DBL_MAX;
                    double multiplier = 1.0;
                    if (nNegInt == 1) {
                      // swap signs
                      multiplier = lo;
                      lo = -up;
                      up = -multiplier;
                      multiplier = -1.0;
                    }
                    bool possible = true;
                    double singletonValue = 0;
                    double scale = 4.0 * 5.0 * 6.0;
                    int kGcd = -1;
                    double smallestSum = 0.0;
                    //double largestSum = 0.0;
                    for (k = start; k < end; ++k) {
                      int j = column[k];
                      double value = multiplier * rowElements[k];
                      if (columnUpper[j] > columnLower[j]) {
                        if (value > 0.0) {
                          // singleton
                          kInt = j;
                          if (columnUpper[j] - columnLower[j] != 1.0) {
                            possible = false;
                            break;
                          } else {
                            singletonValue = value;
                          }
                        } else {
                          if (columnLower[j] > -1.0e10)
                            smallestSum += value * columnLower[j];
                          else
                            smallestSum = -COIN_DBL_MAX;
                          //if (columnUpper[j] < -1.0e10)
			  //  largestSum += value * columnUpper[j];
                          //else
			  //  largestSum = COIN_DBL_MAX;
                          value *= -scale;
                          if (fabs(value - floor(value + 0.5)) > 1.0e-12) {
                            possible = false;
                            break;
                          } else {
                            int kVal = static_cast< int >(floor(value + 0.5));
                            if (kGcd > 0)
                              kGcd = gcd(kGcd, kVal);
                            else
                              kGcd = kVal;
                          }
                        }
                      }
                    }
                    if (possible) {
                      double multiple = (static_cast< double >(kGcd)) / scale;
                      int interesting = 0;
                      double saveLo = lo;
                      double saveUp = up;
#if CBC_USEFUL_PRINTING > 1
                      double nearestLo0 = lo;
                      double nearestLo1 = lo;
#endif
                      double nearestUp0 = up;
                      double nearestUp1 = up;
                      // adjust rhs for singleton
                      if (lo != -COIN_DBL_MAX) {
                        // singleton at lb
                        lo -= columnLower[kInt] * singletonValue;
                        double exact = lo / multiple;
                        if (fabs(exact - floor(exact + 0.5)) > 1.0e-4) {
                          interesting += 1;
#if CBC_USEFUL_PRINTING > 1
                          nearestLo0 = ceil(exact) * multiple;
#endif
                        }
                        // singleton at ub
                        lo -= singletonValue;
                        exact = lo / multiple;
                        if (fabs(exact - floor(exact + 0.5)) > 1.0e-4) {
                          interesting += 2;
#if CBC_USEFUL_PRINTING > 1
                          nearestLo1 = ceil(exact) * multiple;
#endif
                        }
                      }
                      if (up != COIN_DBL_MAX) {
                        // singleton at lb
                        up -= columnLower[kInt] * singletonValue;
                        double exact = up / multiple;
                        if (fabs(exact - floor(exact + 0.5)) > 1.0e-4) {
                          interesting += 4;
                          nearestUp0 = floor(exact) * multiple;
                        }
                        // singleton at ub
                        up -= singletonValue;
                        exact = up / multiple;
                        if (fabs(exact - floor(exact + 0.5)) > 1.0e-4) {
                          interesting += 8;
                          nearestUp1 = floor(exact) * multiple;
                        }
                      }
                      if (interesting) {
#if CBC_USEFUL_PRINTING > 1
                        printf("Row %d interesting %d lo,up %g,%g singleton %d value %g bounds %g,%g - gcd %g\n", iRow, interesting, saveLo, saveUp, kInt, singletonValue,
                          columnLower[kInt], columnUpper[kInt], multiple);
                        printf("Orig lo,up %g,%g %d\n", rowLower[iRow], rowUpper[iRow], end - start);
                        for (k = start; k < end; ++k) {
                          CoinBigIndex j = column[k];
                          double value = multiplier * rowElements[k];
                          printf(" (%d, %g - bds %g, %g)", j, value,
                            columnLower[j], columnUpper[j]);
                        }
                        printf("\n");
#endif
                        if (columnLower[kInt]) {
#if CBC_USEFUL_PRINTING > 1
                          printf("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ\n"); //think
#endif
                          interesting = 0;
                        }
                        newValue = singletonValue;
#if CBC_USEFUL_PRINTING > 1
                        double newLoRhs = rowLower[iRow];
                        double newUpRhs = rowUpper[iRow];
                        if ((interesting & 3) != 0) {
                          newLoRhs = nearestLo0;
                          newValue = nearestLo0 - nearestLo1;
                        }
#endif
                        if (saveLo == saveUp && ((interesting & 5) == 5 || (interesting & 10) == 10)) {
#if CBC_USEFUL_PRINTING > 1
                          printf("INFEAS? ");
#endif
                          interesting = 0; //ninfeas++;
                        }
                        if ((interesting & 12)) {
#if CBC_USEFUL_PRINTING > 1
                          double value2 = newValue;
                          newUpRhs = nearestUp0;
#endif
                          newValue = nearestUp0 - nearestUp1;
#if CBC_USEFUL_PRINTING > 1
                          if (newValue != value2) {
                            printf("??? old newvalue %g ", newValue);
                          }
#endif
                        }
#if CBC_USEFUL_PRINTING > 1
                        printf("guess is new lo %g, new up %g, new value %g\n",
                          newLoRhs, newUpRhs, newValue);
#endif
                      }
                      // Just do mzzv11 case
                      double exact = singletonValue / multiple;
                      if (fabs(exact - floor(exact + 0.5)) < 1.0e-5)
                        interesting &= ~2;
                      if (!smallestSum && interesting == 2 && !saveLo && saveUp > 1.0e20) {
                        newValue = multiple * floor(exact);
                        newValue *= multiplier;
#if CBC_USEFUL_PRINTING > 1
                        printf("New coefficient for %d will be %g\n", kInt, newValue);
#endif
                      } else {
                        // don't do
                        kInt = -1;
                      }
                    } else {
                      kInt = -1;
                    }
                  }
                }
                // endgcd
              }
              int length = rowLength[iRow];
              if (kInt < 0) {
                build.addRow(length, column + start, rowElements + start,
                  rowLower[iRow], rowUpper[iRow]);
              } else {
                double *els = CoinCopyOfArray(rowElements + start, length);
                if (fabs(newValue) > 1.0e-13) {
                  for (int k = 0; k < length; ++k) {
                    int j = column[k + start];
                    if (j == kInt) {
                      els[k] = newValue;
                      break;
                    }
                  }
                  build.addRow(length, column + start, els,
                    rowLower[iRow], rowUpper[iRow]);
                } else {
                  // strengthened to zero!
#if CBC_USEFUL_PRINTING > 1
                  printf("CglPreProcess - element strenthened to zero!\n");
#endif
                  int *cols = CoinCopyOfArray(column + start, length);
                  int n = 0;
                  for (int k = 0; k < length; ++k) {
                    int j = cols[k];
                    if (j != kInt) {
                      els[n] = els[k];
                      cols[n++] = j;
                    }
                  }
                  build.addRow(n, cols, els,
                    rowLower[iRow], rowUpper[iRow]);
                  delete[] cols;
                }
                delete[] els;
              }
              keepRow[iRow] = 1;
            } else {
              // strengthens this row (should we check?)
              // may be worth going round again
              numberChangedThisPass++;
              int n = thisCut->row().getNumElements();
              const int *columnCut = thisCut->row().getIndices();
              const double *elementCut = thisCut->row().getElements();
              double lower = thisCut->lb();
              double upper = thisCut->ub();
              if (probingCut) {
                int i;
                int nFree = 0;
                for (i = 0; i < n; i++) {
                  int iColumn = columnCut[i];
                  if (columnUpper[iColumn] > columnLower[iColumn] + 1.0e-12)
                    nFree++;
                }
                bool good = (n == nFree);
                nFree = 0;
                int n1 = rowLength[iRow];
                CoinBigIndex start = rowStart[iRow];
                for (i = 0; i < n1; i++) {
                  int iColumn = column[start + i];
                  if (columnUpper[iColumn] > columnLower[iColumn] + 1.0e-12)
                    nFree++;
                }
                if (n1 != nFree)
                  good = false;
                if (noStrengthening && n > n1) {
                  good = false;
                  //printf("Skipping row %d strengthened row has %d (%d)\n",
                  //	 iRow,n,n1);
                }
                if (good) {
#if PRINT_DEBUG > 1
                  printf("Original row %.8d %g <= ", iRow, rowLower[iRow]);
                  for (i = 0; i < n1; i++)
                    printf("%g * x%d ", rowElements[start + i], column[start + i]);
                  printf("<= %g\n", rowUpper[iRow]);
                  printf("New                   %g <= ", lower);
                  for (i = 0; i < n; i++)
                    printf("%g * x%d ", elementCut[i], columnCut[i]);
                  printf("<= %g\n", upper);
#endif
                } else {
                  // can't use
                  n = -1;
                  numberStrengthened--;
                  // put in old row
                  CoinBigIndex start = rowStart[iRow];
                  build.addRow(rowLength[iRow], column + start, rowElements + start,
                    rowLower[iRow], rowUpper[iRow]);
                  keepRow[iRow] = 1;
                }
              }
              if (n > 0) {
                build.addRow(n, columnCut, elementCut, lower, upper);
                keepRow[iRow] = 1;
              } else if (!n) {
                // Either infeasible or redundant
                if (lower <= 0.0 && upper >= 0.0) {
                  // redundant - row will go
                } else {
                  // infeasible!
                  feasible = false;
                  break;
                }
              }
            }
          }
        }
        // Get rid of all whichCut not in cs
        int nCuts = cs.sizeRowCuts();
        if (nCuts) {
          std::sort(whichCut, whichCut + numberRows);
          OsiRowCut **tempCuts = new OsiRowCut *[nCuts + 1];
          for (int i = 0; i < nCuts; i++)
            tempCuts[i] = cs.rowCutPtr(i);
          std::sort(tempCuts, tempCuts + nCuts);
          tempCuts[nCuts] = std::max(whichCut[numberRows - 1], tempCuts[nCuts - 1]) + 1;
          int iCut = 0;
          void *cut = tempCuts[0];
          for (int i = 0; i < numberRows; i++) {
            OsiRowCut *thisCut = whichCut[i];
            if (!thisCut)
              continue;
            while (cut < thisCut)
              cut = tempCuts[++iCut];
            if (cut > thisCut)
              delete thisCut;
            whichCut[i] = NULL;
          }
          delete[] tempCuts;
        } else {
          for (int i = 0; i < numberRows; i++) {
            delete whichCut[i];
            whichCut[i] = NULL;
          }
        }
        if (rowType_) {
          assert(numberRowType_ == numberRows);
          int numberRowType_ = 0;
          for (iRow = 0; iRow < numberRows; iRow++) {
            if (keepRow[iRow]) {
              rowType_[numberRowType_++] = rowType_[iRow];
            }
          }
        }
        // recreate
        int *del = new int[numberRows];
        // save and modify basis
        CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(newModel->getWarmStart());
        if (basis) {
          int n = 0;
          for (iRow = 0; iRow < numberRows; iRow++) {
            if (!keepRow[iRow]) {
              del[n++] = iRow;
            }
          }
          basis->deleteRows(n, del);
        }
        for (iRow = 0; iRow < numberRows; iRow++)
          del[iRow] = iRow;
        newModel->deleteRows(numberRows, del);
        newModel->addRows(build);
        numberRows = newModel->getNumRows();
#if CBC_USEFUL_PRINTING > 0
        printf("After build %d rows and %d columns (%d rows strengthened)\n",
          numberRows, newModel->getNumCols(), numberStrengthened);
#endif
        if (basis) {
          assert(numberRows == basis->getNumArtificial());
          newModel->setWarmStart(basis);
          delete basis;
        }
        if (dupRow && cs.sizeRowCuts())
          newModel->applyCuts(cs);
        delete[] keepRow;
        delete[] del;
        /*
  VIRTUOUS

  A solver is not required to hold these pointers constant across complete
  row deletion and rebuild.
*/
        columnLower = newModel->getColLower();
        columnUpper = newModel->getColUpper();
      }
      if (!feasible)
        break;
      // now see if we have any x=y x+y=1
      if (constraints) {
        int numberRowCuts = cs.sizeRowCuts();
        for (int k = 0; k < numberRowCuts; k++) {
          OsiRowCut *thisCut = cs.rowCutPtr(k);
          int n = thisCut->row().getNumElements();
          double lower = thisCut->lb();
          double upper = thisCut->ub();
          if (n == 2 && lower == upper) {
            twoCuts.insert(*thisCut);
          }
        }
      }
      if (probingCut) {
        // we could look harder for infeasibilities
        assert(info.numberVariables() == numberColumns);
        int number01 = info.numberIntegers();
        const CliqueEntry *entry = info.fixEntries();
        const int *toZero = info.toZero();
        const int *toOne = info.toOne();
        const int *which = info.integerVariable();
        int numberBounds = 0;
        char *markLB = new char[number01];
        memset(markLB, 0, number01);
        char *markUB = new char[number01];
        memset(markUB, 0, number01);
        for (int k = 0; k < number01; k++) {
          int start = toZero[k];
          int end = toOne[k];
          // to zero
          int j;
          for (j = start; j < end; j++) {
            int goingToOne = oneFixesInCliqueEntry(entry[j]);
            int v = sequenceInCliqueEntry(entry[j]);
            if (v >= number01)
              continue;
            if (goingToOne) {
              // If v -> 1 means k -> 0 we have k+v==1
              int startV = toOne[v];
              int endV = toZero[v + 1];
              for (int jv = startV; jv < endV; jv++) {
                if (k == static_cast< int >(sequenceInCliqueEntry(entry[jv]))) {
                  int goingToOneV = oneFixesInCliqueEntry(entry[jv]);
                  double lo, up;
                  if (!goingToOneV) {
                    lo = 1.0;
                    up = 1.0;
                    OsiRowCut thisCut;
                    thisCut.setLb(lo);
                    thisCut.setUb(up);
                    double values[] = { 1.0, 1.0 };
                    int columns[2];
                    columns[0] = which[k];
                    columns[1] = which[v];
                    thisCut.setRow(2, columns, values, false);
                    twoCuts.insertIfNotDuplicate(thisCut);
                  } else {
                    // means infeasible for k to go to 0
                    markLB[k] = 1;
                    numberBounds++;
                  }
                  break;
                }
              }
            } else {
              // If v -> 0 means k -> 0 we have k==v
              int startV = toZero[v];
              int endV = toOne[v];
              for (int jv = startV; jv < endV; jv++) {
                if (k == static_cast< int >(sequenceInCliqueEntry(entry[jv]))) {
                  int goingToOneV = oneFixesInCliqueEntry(entry[jv]);
                  double lo, up;
                  if (!goingToOneV) {
                    lo = 0.0;
                    up = 0.0;
                    OsiRowCut thisCut;
                    thisCut.setLb(lo);
                    thisCut.setUb(up);
                    double values[] = { 1.0, -1.0 };
                    int columns[2];
                    columns[0] = which[k];
                    columns[1] = which[v];
                    thisCut.setRow(2, columns, values, false);
                    twoCuts.insertIfNotDuplicate(thisCut);
                  } else {
                    // means infeasible for k to go to 0
                    markLB[k] = 1;
                    numberBounds++;
                  }
                  break;
                }
              }
            }
          }
          start = toOne[k];
          end = toZero[k + 1];
          // to one
          for (j = start; j < end; j++) {
            int goingToOne = oneFixesInCliqueEntry(entry[j]);
            int v = sequenceInCliqueEntry(entry[j]);
            if (v >= number01)
              continue;
            if (goingToOne) {
              // If v -> 1 means k -> 1 we have k==v
              int startV = toOne[v];
              int endV = toZero[v + 1];
              for (int jv = startV; jv < endV; jv++) {
                if (k == static_cast< int >(sequenceInCliqueEntry(entry[jv]))) {
                  int goingToOneV = oneFixesInCliqueEntry(entry[jv]);
                  double lo, up;
                  if (goingToOneV) {
                    lo = 0.0;
                    up = 0.0;
                    OsiRowCut thisCut;
                    thisCut.setLb(lo);
                    thisCut.setUb(up);
                    double values[] = { 1.0, -1.0 };
                    int columns[2];
                    columns[0] = which[k];
                    columns[1] = which[v];
                    thisCut.setRow(2, columns, values, false);
                    twoCuts.insertIfNotDuplicate(thisCut);
                  } else {
                    // means infeasible for k to go to 1
                    markUB[k] = 1;
                    numberBounds++;
                  }
                  break;
                }
              }
            } else {
              // If v -> 0 means k -> 1 we have k+v==1
              int startV = toZero[v];
              int endV = toOne[v];
              for (int jv = startV; jv < endV; jv++) {
                if (k == static_cast< int >(sequenceInCliqueEntry(entry[jv]))) {
                  int goingToOneV = oneFixesInCliqueEntry(entry[jv]);
                  double lo, up;
                  if (goingToOneV) {
                    lo = 1.0;
                    up = 1.0;
                    OsiRowCut thisCut;
                    thisCut.setLb(lo);
                    thisCut.setUb(up);
                    double values[] = { 1.0, 1.0 };
                    int columns[2];
                    columns[0] = which[k];
                    columns[1] = which[v];
                    thisCut.setRow(2, columns, values, false);
                    twoCuts.insertIfNotDuplicate(thisCut);
                  } else {
                    // means infeasible for k to go to 1
                    markUB[k] = 1;
                    numberBounds++;
                  }
                  break;
                }
              }
            }
          }
        }
        if (numberBounds) {
          CoinPackedVector lbs;
          CoinPackedVector ubs;
          for (int k = 0; k < number01; k++) {
            if (markLB[k] && markUB[k]) {
              // infeasible
              feasible = false;
              break;
            } else if (markLB[k]) {
              lbs.insert(which[k], 1.0);
            } else if (markUB[k]) {
              ubs.insert(which[k], 0.0);
            }
          }
          OsiColCut cc;
          cc.setUbs(ubs);
          cc.setLbs(lbs);
          cc.setEffectiveness(1.0e-5);
          cs.insert(cc);
        }
        delete[] markLB;
        delete[] markUB;
      }
      // see if we have any column cuts
      int numberColumnCuts = cs.sizeColCuts();
      int numberBounds = 0;
      writeDebugMps(newModel, "beforecolcut", NULL);
      for (int k = 0; k < numberColumnCuts; k++) {
        OsiColCut *thisCut = cs.colCutPtr(k);
        /*
	  Nontrivial bound changes will invalidate current solution.
	*/
        if (thisCut->effectiveness() > 1.0 && useSolution) {
          needResolve = true;
        }
        const CoinPackedVector &lbs = thisCut->lbs();
        const CoinPackedVector &ubs = thisCut->ubs();
        int j;
        int n;
        const int *which;
        const double *values;
        n = lbs.getNumElements();
        which = lbs.getIndices();
        values = lbs.getElements();
        for (j = 0; j < n; j++) {
          int iColumn = which[j];
          if (values[j] > columnLower[iColumn] && values[j] > -1.0e20) {
            //printf("%d lower from %g to %g\n",iColumn,columnLower[iColumn],values[j]);
            newModel->setColLower(iColumn, values[j]);
            if (false) {
              OsiSolverInterface *xx = newModel->clone();
              xx->initialSolve();
              assert(xx->isProvenOptimal());
              delete xx;
            }
            numberChangedThisPass++;
            if (columnLower[iColumn] == columnUpper[iColumn]) {
              numberFixed++;
            } else {
              numberBounds++;
              if (columnLower[iColumn] > columnUpper[iColumn] + 1.0e-12)
                feasible = false;
            }
          }
        }
        n = ubs.getNumElements();
        which = ubs.getIndices();
        values = ubs.getElements();
        for (j = 0; j < n; j++) {
          int iColumn = which[j];
          if (values[j] < columnUpper[iColumn] && values[j] < 1.0e20) {
            //printf("%d upper from %g to %g\n",iColumn,columnUpper[iColumn],values[j]);
            newModel->setColUpper(iColumn, values[j]);
            if (false) {
              OsiSolverInterface *xx = newModel->clone();
              xx->initialSolve();
              assert(xx->isProvenOptimal());
              delete xx;
            }
            numberChangedThisPass++;
            if (columnLower[iColumn] == columnUpper[iColumn]) {
              numberFixed++;
            } else {
              numberBounds++;
              if (columnLower[iColumn] > columnUpper[iColumn]+1.0e-12)
                feasible = false;
            }
          }
        }
      }
      writeDebugMps(newModel, "aftercolcut", NULL);
      numberTwo = twoCuts.sizeRowCuts() - numberTwo;
      numberChanges += numberTwo + numberStrengthened / 10;
      if (numberFixed || numberTwo || numberStrengthened || numberBounds)
        handler_->message(CGL_PROCESS_STATS, messages_)
          << numberFixed << numberBounds << numberStrengthened << numberTwo
          << CoinMessageEol;
      if (!feasible)
        break;
      /*
	If solution needs to be refreshed, do resolve or initialSolve as appropriate.
      */
      if (needResolve) {
        if (rebuilt) {
          // basis shot to bits?
          //CoinWarmStartBasis *slack =
          //dynamic_cast<CoinWarmStartBasis *>(newModel->getEmptyWarmStart()) ;
          //newModel->setWarmStart(slack);
          //delete slack ;
          bool saveHint;
          OsiHintStrength saveStrength;
          newModel->getHintParam(OsiDoPresolveInInitial, saveHint, saveStrength);
          // Do presolves
          if ((numberFixed + numberTwo) * 4 > numberColumns)
            newModel->setHintParam(OsiDoPresolveInInitial, true, OsiHintTry);
          newModel->setHintParam(OsiDoDualInInitial, true, OsiHintTry);
          newModel->initialSolve();
          newModel->setHintParam(OsiDoPresolveInInitial, saveHint, saveStrength);
        } else {
	  bool saveTakeHint;
	  OsiHintStrength saveStrength;
	  newModel->getHintParam(OsiDoDualInResolve,
				 saveTakeHint, saveStrength);
	  newModel->setHintParam(OsiDoDualInResolve, solveWithDual, OsiHintTry);
	  newModel->resolve();
	  newModel->setHintParam(OsiDoDualInResolve, saveTakeHint, saveStrength);
	  solveWithDual = true;
        }
        numberIterationsPre_ += newModel->getIterationCount();
        feasible = newModel->isProvenOptimal();
        if (!feasible) {
          // Better double check
          CoinWarmStartBasis *slack = dynamic_cast< CoinWarmStartBasis * >(newModel->getEmptyWarmStart());
          newModel->setWarmStart(slack);
          delete slack;
          newModel->resolve();
          numberIterationsPre_ += newModel->getIterationCount();
          feasible = newModel->isProvenOptimal();
          //if (!feasible)
          //newModel->writeMpsNative("infeas.mps",NULL,NULL,2,1);
        }
      }
      if (!feasible) {
        writeDebugMps(newModel, "infeasible", NULL);
        break;
      }
    }
    if (!feasible)
      break;
    numberChanges += numberChangedThisPass;
    if (iPass < numberPasses - 1) {
      int multiplier = (numberPasses > 10) ? numberRows + 1000 : 1;
      if ((!numberFixed && numberChangedThisPass * multiplier < numberRows + numberColumns) || iPass == numberPasses - 2) {
        // do special probing at end - but not if very last pass
        if (iBigPass < numberSolvers_ - 1 && numberPasses > 4) {
          firstGenerator = -1;
          lastGenerator = 0;
        }
        iPass = numberPasses - 2;
      }
    }
  }
  delete[] whichCut;
  int numberRowCuts = twoCuts.sizeRowCuts();
  if (numberRowCuts) {
    writeDebugMps(newModel, "beforeaddtwo", NULL);
    // add in x=y etc
    CoinBuild build;
    for (int k = 0; k < numberRowCuts; k++) {
      OsiRowCut *thisCut = twoCuts.rowCutPtr(k);
      int n = thisCut->row().getNumElements();
      const int *columnCut = thisCut->row().getIndices();
      const double *elementCut = thisCut->row().getElements();
      double lower = thisCut->lb();
      double upper = thisCut->ub();
      build.addRow(n, columnCut, elementCut, lower, upper);
#if DEBUG_PREPROCESS > 1
      if (debugSolution) {
	double high = 0.0;
	double low = 0.0;
	const double *lowerM = newModel->getColLower();
	const double *upperM = newModel->getColUpper();
	for (int i=0;i<n;i++) {
	  int iColumn = columnCut[i];
	  double value = debugSolution[iColumn];
	  double el = elementCut[i];
	  if (newModel->isInteger(iColumn)) {
	    value = floor(value+0.5);
	    high += value*el;
	    low += value*el;
	  } else if (el > 0.0) {
	    high += upperM[iColumn]*el;
	    low += lowerM[iColumn]*el;
	  } else {
	    low += upperM[iColumn]*el;
	    high += lowerM[iColumn]*el;
	  }
	}
	if (high < lower-1.0e-4 || low > upper+1.0e-4) {
	  printf("bad cut\n");
	}
      }
#endif
    }
    newModel->addRows(build);
    writeDebugMps(newModel, "afteraddtwo", NULL);
    if (rowType_) {
      // adjust
      int numberRows = newModel->getNumRows();
      char *temp = CoinCopyOfArrayPartial(rowType_, numberRows, numberRowType_);
      delete[] rowType_;
      rowType_ = temp;
      for (int iRow = numberRowType_; iRow < numberRows; iRow++)
        rowType_[iRow] = -1;
      numberRowType_ = numberRows;
    }
    writeDebugMps(newModel, "afterBaddtwo", NULL);
  }
  if (!feasible) {
    delete newModel;
    newModel = NULL;
  }
  return newModel;
}

/* Default Constructor

*/
CglPreProcess::CglPreProcess()
  : originalModel_(NULL)
  , startModel_(NULL)
  , numberSolvers_(0)
  , model_(NULL)
  , modifiedModel_(NULL)
  , presolve_(NULL)
  , handler_(NULL)
  , defaultHandler_(true)
  , appData_(NULL)
  , originalColumn_(NULL)
  , originalRow_(NULL)
  , numberCutGenerators_(0)
  , generator_(NULL)
  , numberSOS_(0)
  , typeSOS_(NULL)
  , startSOS_(NULL)
  , whichSOS_(NULL)
  , weightSOS_(NULL)
  , numberProhibited_(0)
  , numberIterationsPre_(0)
  , numberIterationsPost_(0)
  , prohibited_(NULL)
  , numberRowType_(0)
  , options_(0)
  , rowType_(NULL)
  , useElapsedTime_(true)
  , timeLimit_(COIN_DBL_MAX)
  , keepColumnNames_(false)
  , inspect_(false)
{
  handler_ = new CoinMessageHandler();
  handler_->setLogLevel(2);
  messages_ = CglMessage();
}

CglPreProcess::CglPreProcess(const CglPreProcess &rhs)
  : numberSolvers_(rhs.numberSolvers_)
  , defaultHandler_(rhs.defaultHandler_)
  , appData_(rhs.appData_)
  , originalColumn_(NULL)
  , originalRow_(NULL)
  , numberCutGenerators_(rhs.numberCutGenerators_)
  , numberProhibited_(rhs.numberProhibited_)
  , numberIterationsPre_(rhs.numberIterationsPre_)
  , numberIterationsPost_(rhs.numberIterationsPost_)
  , numberRowType_(rhs.numberRowType_)
  , options_(rhs.options_)
  , useElapsedTime_(true)
  , timeLimit_(COIN_DBL_MAX)
  , keepColumnNames_(false)
  , inspect_(rhs.inspect_)
{
  if (defaultHandler_) {
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(rhs.handler_->logLevel());
  } else {
    handler_ = rhs.handler_;
  }
  messages_ = rhs.messages_;
  if (numberCutGenerators_) {
    generator_ = new CglCutGenerator *[numberCutGenerators_];
    for (int i = 0; i < numberCutGenerators_; i++) {
      generator_[i] = rhs.generator_[i]->clone();
    }
  } else {
    generator_ = NULL;
  }
  if (rhs.originalModel_) {
    originalModel_ = rhs.originalModel_;
    // If no make equality then solvers are same
    if (rhs.originalModel_ != rhs.startModel_) {
      startModel_ = rhs.startModel_->clone();
    } else {
      startModel_ = originalModel_;
    }
  } else {
    originalModel_ = NULL;
    startModel_ = NULL;
  }
  if (numberSolvers_) {
    model_ = new OsiSolverInterface *[numberSolvers_];
    modifiedModel_ = new OsiSolverInterface *[numberSolvers_];
    presolve_ = new OsiPresolve *[numberSolvers_];
    for (int i = 0; i < numberSolvers_; i++) {
      model_[i] = rhs.model_[i]->clone();
      modifiedModel_[i] = rhs.modifiedModel_[i]->clone();
      presolve_[i] = new OsiPresolve(*rhs.presolve_[i]);
    }
  } else {
    model_ = NULL;
    presolve_ = NULL;
  }
  numberSOS_ = rhs.numberSOS_;
  if (numberSOS_) {
    int numberTotal = rhs.startSOS_[numberSOS_];
    typeSOS_ = CoinCopyOfArray(rhs.typeSOS_, numberSOS_);
    startSOS_ = CoinCopyOfArray(rhs.startSOS_, numberSOS_ + 1);
    whichSOS_ = CoinCopyOfArray(rhs.whichSOS_, numberTotal);
    weightSOS_ = CoinCopyOfArray(rhs.weightSOS_, numberTotal);
  } else {
    typeSOS_ = NULL;
    startSOS_ = NULL;
    whichSOS_ = NULL;
    weightSOS_ = NULL;
  }
  prohibited_ = CoinCopyOfArray(rhs.prohibited_, numberProhibited_);
  rowType_ = CoinCopyOfArray(rhs.rowType_, numberRowType_);
  cuts_ = rhs.cuts_;
}

// Assignment operator
CglPreProcess &
CglPreProcess::operator=(const CglPreProcess &rhs)
{
  if (this != &rhs) {
    gutsOfDestructor();
    numberSolvers_ = rhs.numberSolvers_;
    defaultHandler_ = rhs.defaultHandler_;
    appData_ = rhs.appData_;
    numberCutGenerators_ = rhs.numberCutGenerators_;
    numberProhibited_ = rhs.numberProhibited_;
    numberIterationsPre_ = rhs.numberIterationsPre_;
    numberIterationsPost_ = rhs.numberIterationsPost_;
    numberRowType_ = rhs.numberRowType_;
    options_ = rhs.options_;
    inspect_ = rhs.inspect_;
    if (defaultHandler_) {
      handler_ = new CoinMessageHandler();
      handler_->setLogLevel(rhs.handler_->logLevel());
    } else {
      handler_ = rhs.handler_;
    }
    messages_ = rhs.messages_;
    if (numberCutGenerators_) {
      generator_ = new CglCutGenerator *[numberCutGenerators_];
      for (int i = 0; i < numberCutGenerators_; i++) {
        generator_[i] = rhs.generator_[i]->clone();
      }
    } else {
      generator_ = NULL;
    }
    if (rhs.originalModel_) {
      originalModel_ = rhs.originalModel_;
      // If no make equality then solvers are same
      if (rhs.originalModel_ != rhs.startModel_) {
        startModel_ = rhs.startModel_->clone();
      } else {
        startModel_ = originalModel_;
      }
    } else {
      originalModel_ = NULL;
      startModel_ = NULL;
    }
    if (numberSolvers_) {
      model_ = new OsiSolverInterface *[numberSolvers_];
      modifiedModel_ = new OsiSolverInterface *[numberSolvers_];
      presolve_ = new OsiPresolve *[numberSolvers_];
      for (int i = 0; i < numberSolvers_; i++) {
        model_[i] = rhs.model_[i]->clone();
        modifiedModel_[i] = rhs.modifiedModel_[i]->clone();
        presolve_[i] = new OsiPresolve(*rhs.presolve_[i]);
      }
    } else {
      model_ = NULL;
      presolve_ = NULL;
    }
    numberSOS_ = rhs.numberSOS_;
    if (numberSOS_) {
      int numberTotal = rhs.startSOS_[numberSOS_];
      typeSOS_ = CoinCopyOfArray(rhs.typeSOS_, numberSOS_);
      startSOS_ = CoinCopyOfArray(rhs.startSOS_, numberSOS_ + 1);
      whichSOS_ = CoinCopyOfArray(rhs.whichSOS_, numberTotal);
      weightSOS_ = CoinCopyOfArray(rhs.weightSOS_, numberTotal);
    } else {
      typeSOS_ = NULL;
      startSOS_ = NULL;
      whichSOS_ = NULL;
      weightSOS_ = NULL;
    }
    prohibited_ = CoinCopyOfArray(rhs.prohibited_, numberProhibited_);
    rowType_ = CoinCopyOfArray(rhs.rowType_, numberRowType_);
    cuts_ = rhs.cuts_;
    timeLimit_ = rhs.timeLimit_;
    keepColumnNames_ = rhs.keepColumnNames_;
  }
  return *this;
}

// Destructor
CglPreProcess::~CglPreProcess()
{
  gutsOfDestructor();
}
// Clears out as much as possible (except solver)
void CglPreProcess::gutsOfDestructor()
{
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
  if (startModel_ != originalModel_)
    delete startModel_;
  startModel_ = NULL;
  //delete originalModel_;
  originalModel_ = NULL;
  int i;
  for (i = 0; i < numberCutGenerators_; i++) {
    delete generator_[i];
  }
  delete[] generator_;
  generator_ = NULL;
  if (numberSolvers_==99) 
    numberSolvers_ = 1;
  for (i = 0; i < numberSolvers_; i++) {
    delete model_[i];
    delete modifiedModel_[i];
    delete presolve_[i];
  }
  delete[] model_;
  delete[] modifiedModel_;
  delete[] presolve_;
  model_ = NULL;
  presolve_ = NULL;
  delete[] originalColumn_;
  delete[] originalRow_;
  originalColumn_ = NULL;
  originalRow_ = NULL;
  delete[] typeSOS_;
  delete[] startSOS_;
  delete[] whichSOS_;
  delete[] weightSOS_;
  typeSOS_ = NULL;
  startSOS_ = NULL;
  whichSOS_ = NULL;
  weightSOS_ = NULL;
  delete[] prohibited_;
  prohibited_ = NULL;
  numberProhibited_ = 0;
  numberIterationsPre_ = 0;
  numberIterationsPost_ = 0;
  delete[] rowType_;
  rowType_ = NULL;
  numberRowType_ = 0;
}
// Clears models
void CglPreProcess::clean()
{
  for (int i = 0; i < numberSolvers_; i++) {
    delete model_[i];
    delete modifiedModel_[i];
    delete presolve_[i];
  }
  delete[] model_;
  delete[] modifiedModel_;
  delete[] presolve_;
  model_ = NULL;
  modifiedModel_ = NULL;
  presolve_ = NULL;
  delete startModel_;
  startModel_ = NULL;
 }
// Add one generator
void CglPreProcess::addCutGenerator(CglCutGenerator *generator)
{
  CglCutGenerator **temp = generator_;
  generator_ = new CglCutGenerator *[numberCutGenerators_ + 1];
  if( numberCutGenerators_ > 0 )
     memcpy(generator_, temp, numberCutGenerators_ * sizeof(CglCutGenerator *));
  delete[] temp;
  generator_[numberCutGenerators_++] = generator->clone();
}
//#############################################################################
// Set/Get Application Data
// This is a pointer that the application can store into and retrieve
// This field is the application to optionally define and use.
//#############################################################################

void CglPreProcess::setApplicationData(void *appData)
{
  appData_ = appData;
}
//-----------------------------------------------------------------------------
void *CglPreProcess::getApplicationData() const
{
  return appData_;
}
/* Set cutoff bound on the objective function.
   
When using strict comparison, the bound is adjusted by a tolerance to
avoid accidentally cutting off the optimal solution.
*/
void CglPreProcess::setCutoff(double value)
{
  // Solvers know about direction
  double direction = originalModel_->getObjSense();
  originalModel_->setDblParam(OsiDualObjectiveLimit, value * direction);
}

// Get the cutoff bound on the objective function - always as minimize
double
CglPreProcess::getCutoff() const
{
  double value;
  originalModel_->getDblParam(OsiDualObjectiveLimit, value);
  return value * originalModel_->getObjSense();
}
// Pass in Message handler (not deleted at end)
void CglPreProcess::passInMessageHandler(CoinMessageHandler *handler)
{
  if (defaultHandler_)
    delete handler_;
  defaultHandler_ = false;
  handler_ = handler;
}
// Set language
void CglPreProcess::newLanguage(CoinMessages::Language language)
{
  messages_ = CglMessage(language);
}
// Return a pointer to the original columns (without clique slacks)
const int *
CglPreProcess::originalColumns()
{
  if (!originalColumn_)
    createOriginalIndices();
  return originalColumn_;
}
// Return a pointer to the original rows
const int *
CglPreProcess::originalRows()
{
  if (!originalRow_)
    createOriginalIndices();
  return originalRow_;
}
// create original columns and rows
void CglPreProcess::createOriginalIndices()
{
  // Find last model and presolve
  int iPass;
  for (iPass = numberSolvers_ - 1; iPass >= 0; iPass--) {
    if (presolve_[iPass])
      break;
  }
  int nRows, nColumns;
  if (iPass >= 0) {
    nRows = model_[iPass]->getNumRows();
    nColumns = model_[iPass]->getNumCols();
  } else {
    nRows = originalModel_->getNumRows();
    nColumns = originalModel_->getNumCols();
  }
  delete[] originalColumn_;
  originalColumn_ = new int[nColumns];
  delete[] originalRow_;
  originalRow_ = new int[nRows];
  if (iPass >= 0) {
    memcpy(originalColumn_, presolve_[iPass]->originalColumns(),
      nColumns * sizeof(int));
    memcpy(originalRow_, presolve_[iPass]->originalRows(),
      nRows * sizeof(int));
    iPass--;
    for (; iPass >= 0; iPass--) {
      const int *originalColumns = presolve_[iPass]->originalColumns();
      int i;
      for (i = 0; i < nColumns; i++)
        originalColumn_[i] = originalColumns[originalColumn_[i]];
      const int *originalRows = presolve_[iPass]->originalRows();
      int nRowsNow = model_[iPass]->getNumRows();
      for (i = 0; i < nRows; i++) {
        int iRow = originalRow_[i];
        if (iRow >= 0 && iRow < nRowsNow)
          originalRow_[i] = originalRows[iRow];
        else
          originalRow_[i] = -1;
      }
    }
#if CBC_USE_PAPILO
    if (this==initialTry.preProcess) {
      int *mapping = initialTry.mapping;
      int n2=initialTry.presolvedModel->numberColumns();
      if ((initialTry.presolveType&64)!=0) {
	for (int i=0;i<nColumns;i++) {
	  int j = originalColumn_[i];
	  int iColumn = mapping[j];
	  originalColumn_[i] = iColumn;
	}
      } else {
	assert (n2<=nColumns);
	int * orig2 = new int[n2];
	for (int i=0;i<n2;i++) {
	  int j = originalColumn_[i];
	  int iColumn = mapping[i];
	  orig2[i] = originalColumn_[iColumn];
	}
	// move to originalColumn_ (and look at leaks)
	memcpy(originalColumn_,orig2,n2*sizeof(int));
	delete [] orig2;
      }
    }
    if (this!=initialTry.preProcess)
#endif
    std::sort(originalColumn_, originalColumn_ + nColumns);
  } else {
    int i;
    for (i = 0; i < nColumns; i++)
      originalColumn_[i] = i;
    for (i = 0; i < nRows; i++)
      originalRow_[i] = i;
  }
}
// Update prohibited and rowType
void CglPreProcess::update(const OsiPresolve *pinfo,
			   const OsiSolverInterface *solver, double * scBound)
{
  if (prohibited_) {
    const int *original = pinfo->originalColumns();
    int numberColumns = solver->getNumCols();
    // number prohibited must stay constant
    int i;
#ifndef NDEBUG
    int n = 0;
    for (i = 0; i < numberProhibited_; i++) {
      if (prohibited_[i])
        n++;
    }
    int n2 = 0;
#endif
    for (i = 0; i < numberColumns; i++) {
      int iColumn = original[i];
      assert(i == 0 || iColumn > original[i - 1]);
      if (scBound)
	scBound[i] = scBound[iColumn];
      char p = prohibited_[iColumn];
#ifndef NDEBUG
      if (p)
        n2++;
#endif
      prohibited_[i] = p;
    }
    assert(n == n2);
    numberProhibited_ = numberColumns;
  }
  if (rowType_) {
    const int *original = pinfo->originalRows();
    int numberRows = solver->getNumRows();
#if CBC_USEFUL_PRINTING > 1
    int nMarked1 = 0;
    for (int i = 0; i < pinfo->getNumRows(); i++) {
      if (rowType_[i])
        nMarked1++;
    }
    int nMarked2 = 0;
    int k = -1;
    for (int i = 0; i < numberRows; i++) {
      int iRow = original[i];
      if (iRow < i)
        abort();
      if (iRow <= k)
        abort();
      k = iRow;
      if (rowType_[iRow])
        nMarked2++;
    }
    if (nMarked1 > nMarked2)
      printf("Marked rows reduced from %d to %d\n",
        nMarked1, nMarked2);
#endif
    for (int i = 0; i < numberRows; i++) {
      int iRow = original[i];
      rowType_[i] = rowType_[iRow];
    }
    numberRowType_ = numberRows;
  }
}
/* Fix some of problem - returning new problem.
   Uses reduced costs.
   Optional signed character array
   1 always keep, -1 always discard, 0 use djs
   
*/
OsiSolverInterface *
CglPreProcess::someFixed(OsiSolverInterface &model,
  double fractionToKeep,
  bool fixContinuousAsWell,
  char *keep) const
{
  model.resolve();
  int numberColumns = model.getNumCols();
  OsiSolverInterface *newModel = model.clone();
  int i;
  const double *lower = model.getColLower();
  const double *upper = model.getColUpper();
  const double *solution = model.getColSolution();
  double *dj = CoinCopyOfArray(model.getReducedCost(), numberColumns);
  int *sort = new int[numberColumns];
  int number = 0;
  int numberThrow = 0;
  int numberContinuous = 0;
  for (i = 0; i < numberColumns; i++) {
    if (!model.isInteger(i) && upper[i] > lower[i])
      numberContinuous++;
    if (model.isInteger(i) || fixContinuousAsWell) {
      if (keep) {
        if (keep[i] == 1) {
          continue; // always keep
        } else if (keep[i] == -1) {
          // always fix
          dj[number] = -1.0e30;
          numberThrow++;
          sort[number++] = i;
          continue;
        }
      }
      double value = solution[i];
      if (value < lower[i] + 1.0e-8) {
        dj[number] = -dj[i];
        sort[number++] = i;
      } else if (value > upper[number] - 1.0e-8) {
        dj[number] = -dj[i];
        sort[number++] = i;
      }
    }
  }
  CoinSort_2(dj, dj + number, sort);
  int numberToFix = static_cast< int >(numberColumns * (1.0 - fractionToKeep));
  if (!fixContinuousAsWell)
    numberToFix = static_cast< int >((numberColumns - numberContinuous) * (1.0 - fractionToKeep));
  numberToFix = std::max(numberToFix, numberThrow);
  numberToFix = std::min(number, numberToFix);
#if CBC_USEFUL_PRINTING
  printf("%d columns fixed\n", numberToFix);
#endif
  for (i = 0; i < numberToFix; i++) {
    int iColumn = sort[i];
    double value = solution[iColumn];
    if (value < lower[iColumn] + 1.0e-8) {
      newModel->setColUpper(iColumn, lower[iColumn]);
    } else if (value > upper[number] - 1.0e-8) {
      newModel->setColLower(iColumn, lower[iColumn]);
    } else {
      // must be a throw away on - go to lower
      newModel->setColUpper(iColumn, lower[iColumn]);
    }
  }
  delete[] sort;
  delete[] dj;
  return newModel;
}
// If we have a cutoff - fix variables
int CglPreProcess::reducedCostFix(OsiSolverInterface &model)
{
  double cutoff;
  model.getDblParam(OsiDualObjectiveLimit, cutoff);
  double direction = model.getObjSense();
  cutoff *= direction;
  double gap = cutoff - model.getObjValue() * direction;
  double tolerance;
  model.getDblParam(OsiDualTolerance, tolerance);
  if (gap <= 0.0 || fabs(cutoff) > 1.0e20)
    return 0;
  gap += 100.0 * tolerance;
  // not really but thats all we can get
  double integerTolerance;
  model.getDblParam(OsiPrimalTolerance, integerTolerance);
  int numberColumns = model.getNumCols();

  const double *lower = model.getColLower();
  const double *upper = model.getColUpper();
  const double *solution = model.getColSolution();
  const double *reducedCost = model.getReducedCost();

  int numberFixed = 0;
  int iColumn;
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (model.isInteger(iColumn) && upper[iColumn] > lower[iColumn]) {
      double djValue = direction * reducedCost[iColumn];
      if (solution[iColumn] < lower[iColumn] + integerTolerance && djValue > gap) {
        model.setColUpper(iColumn, lower[iColumn]);
        numberFixed++;
      } else if (solution[iColumn] > upper[iColumn] - integerTolerance && -djValue > gap) {
        model.setColLower(iColumn, upper[iColumn]);
        numberFixed++;
      }
    }
  }
  return numberFixed;
}
// Pass in prohibited columns
void CglPreProcess::passInProhibited(const char *prohibited, int numberColumns)
{
  char *temp = prohibited_;
  prohibited_ = CoinCopyOfArray(prohibited, numberColumns);
  if (temp && numberProhibited_ == numberColumns) {
    // merge
    for (int i = 0; i < numberColumns; i++)
      prohibited_[i] |= temp[i];
  }
  numberProhibited_ = numberColumns;
  delete[] temp;
}
/* Pass in row types
   0 normal
   1 cut rows - will be dropped if remain in
   At end of preprocess cut rows will be dropped
   and put into cuts
*/
void CglPreProcess::passInRowTypes(const char *rowTypes, int numberRows)
{
  delete[] rowType_;
  rowType_ = CoinCopyOfArray(rowTypes, numberRows);
  numberRowType_ = numberRows;
  cuts_ = CglStored();
}
// Make continuous variables integer
void CglPreProcess::makeInteger()
{
  // First check if we need to
  //int numberInteger = 0;
  {
    const double *lower = startModel_->getColLower();
    const double *upper = startModel_->getColUpper();
    int numberColumns = startModel_->getNumCols();
    int iColumn;
    int numberContinuous = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (upper[iColumn] > lower[iColumn] + 1.0e-8) {
        if (startModel_->isInteger(iColumn))
          ;//numberInteger++;
        else
          numberContinuous++;
      }
    }
    if (!numberContinuous)
      return;
  }
  OsiSolverInterface *solver = startModel_->clone();
  const double *objective = solver->getObjCoefficients();
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  int numberColumns = solver->getNumCols();
  int numberRows = solver->getNumRows();
  double direction = solver->getObjSense();
  int iRow, iColumn;

  // Row copy
  CoinPackedMatrix matrixByRow(*solver->getMatrixByRow());
  const double *elementByRow = matrixByRow.getElements();
  const int *column = matrixByRow.getIndices();
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();

  // Column copy
  CoinPackedMatrix matrixByCol(*solver->getMatrixByCol());
  const double *element = matrixByCol.getElements();
  const int *row = matrixByCol.getIndices();
  const CoinBigIndex *columnStart = matrixByCol.getVectorStarts();
  const int *columnLength = matrixByCol.getVectorLengths();

  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();

  char *ignore = new char[numberRows];
  int *changed = new int[numberColumns];
  int *which = new int[numberRows];
  double *changeRhs = new double[numberRows];
  memset(changeRhs, 0, numberRows * sizeof(double));
  memset(ignore, 0, numberRows);
  int numberChanged = 0;
  bool finished = false;
  while (!finished) {
    int saveNumberChanged = numberChanged;
    for (iRow = 0; iRow < numberRows; iRow++) {
      int numberContinuous = 0;
      double value1 = 0.0, value2 = 0.0;
      bool allIntegerCoeff = true;
      double sumFixed = 0.0;
      int jColumn1 = -1, jColumn2 = -1;
      for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
        int jColumn = column[j];
        double value = elementByRow[j];
        if (upper[jColumn] > lower[jColumn] + 1.0e-8) {
          if (!solver->isInteger(jColumn)) {
            if (numberContinuous == 0) {
              jColumn1 = jColumn;
              value1 = value;
            } else {
              jColumn2 = jColumn;
              value2 = value;
            }
            numberContinuous++;
          } else {
            if (fabs(value - floor(value + 0.5)) > 1.0e-12)
              allIntegerCoeff = false;
          }
        } else {
          sumFixed += lower[jColumn] * value;
        }
      }
      double low = rowLower[iRow];
      if (low > -1.0e20) {
        low -= sumFixed;
        if (fabs(low - floor(low + 0.5)) > 1.0e-12)
          allIntegerCoeff = false;
      }
      double up = rowUpper[iRow];
      if (up < 1.0e20) {
        up -= sumFixed;
        if (fabs(up - floor(up + 0.5)) > 1.0e-12)
          allIntegerCoeff = false;
      }
      if (!allIntegerCoeff)
        continue; // can't do
      if (numberContinuous == 1) {
        // see if really integer
        // This does not allow for complicated cases
        if (low == up) {
          if (fabs(value1) > 1.0e-3) {
            value1 = 1.0 / value1;
            if (fabs(value1 - floor(value1 + 0.5)) < 1.0e-12) {
              // integer
              changed[numberChanged++] = jColumn1;
              solver->setInteger(jColumn1);
              if (upper[jColumn1] > 1.0e20)
                solver->setColUpper(jColumn1, 1.0e20);
              if (lower[jColumn1] < -1.0e20)
                solver->setColLower(jColumn1, -1.0e20);
            }
          }
        } else {
          if (fabs(value1) > 1.0e-3) {
            value1 = 1.0 / value1;
            if (fabs(value1 - floor(value1 + 0.5)) < 1.0e-12) {
              // This constraint will not stop it being integer
              ignore[iRow] = 1;
            }
          }
        }
      } else if (numberContinuous == 2) {
        if (low == up) {
          /* need general theory - for now just look at 2 cases -
             1 - +- 1 one in column and just costs i.e. matching objective
             2 - +- 1 two in column but feeds into G/L row which will try and minimize
          */
          if (fabs(value1) == 1.0 && value1 * value2 == -1.0 && !lower[jColumn1]
            && !lower[jColumn2]) {
            int n = 0;
            CoinBigIndex i;
            double objChange = direction * (objective[jColumn1] + objective[jColumn2]);
            double bound = std::min(upper[jColumn1], upper[jColumn2]);
            bound = std::min(bound, 1.0e20);
            for (i = columnStart[jColumn1]; i < columnStart[jColumn1] + columnLength[jColumn1]; i++) {
              int jRow = row[i];
              double value = element[i];
              if (jRow != iRow) {
                which[n++] = jRow;
                changeRhs[jRow] = value;
              }
            }
            for (i = columnStart[jColumn1]; i < columnStart[jColumn1] + columnLength[jColumn1]; i++) {
              int jRow = row[i];
              double value = element[i];
              if (jRow != iRow) {
                if (!changeRhs[jRow]) {
                  which[n++] = jRow;
                  changeRhs[jRow] = value;
                } else {
                  changeRhs[jRow] += value;
                }
              }
            }
            if (objChange >= 0.0) {
              // see if all rows OK
              bool good = true;
              for (i = 0; i < n; i++) {
                int jRow = which[i];
                double value = changeRhs[jRow];
                if (value) {
                  value *= bound;
                  if (rowLength[jRow] == 1) {
                    if (value > 0.0) {
                      double rhs = rowLower[jRow];
                      if (rhs > 0.0) {
                        double ratio = rhs / value;
                        if (fabs(ratio - floor(ratio + 0.5)) > 1.0e-12)
                          good = false;
                      }
                    } else {
                      double rhs = rowUpper[jRow];
                      if (rhs < 0.0) {
                        double ratio = rhs / value;
                        if (fabs(ratio - floor(ratio + 0.5)) > 1.0e-12)
                          good = false;
                      }
                    }
                  } else if (rowLength[jRow] == 2) {
                    if (value > 0.0) {
                      if (rowLower[jRow] > -1.0e20)
                        good = false;
                    } else {
                      if (rowUpper[jRow] < 1.0e20)
                        good = false;
                    }
                  } else {
                    good = false;
                  }
                }
              }
              if (good) {
                // both can be integer
                changed[numberChanged++] = jColumn1;
                solver->setInteger(jColumn1);
                if (upper[jColumn1] > 1.0e20)
                  solver->setColUpper(jColumn1, 1.0e20);
                if (lower[jColumn1] < -1.0e20)
                  solver->setColLower(jColumn1, -1.0e20);
                changed[numberChanged++] = jColumn2;
                solver->setInteger(jColumn2);
                if (upper[jColumn2] > 1.0e20)
                  solver->setColUpper(jColumn2, 1.0e20);
                if (lower[jColumn2] < -1.0e20)
                  solver->setColLower(jColumn2, -1.0e20);
              }
            }
            // clear
            for (i = 0; i < n; i++) {
              changeRhs[which[i]] = 0.0;
            }
          }
        }
      }
    }
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (upper[iColumn] > lower[iColumn] + 1.0e-8 && !solver->isInteger(iColumn)) {
        double value;
        value = upper[iColumn];
        if (value < 1.0e20 && fabs(value - floor(value + 0.5)) > 1.0e-12)
          continue;
        value = lower[iColumn];
        if (value > -1.0e20 && fabs(value - floor(value + 0.5)) > 1.0e-12)
          continue;
        bool integer = true;
        for (CoinBigIndex j = columnStart[iColumn]; j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          if (!ignore[iRow]) {
            integer = false;
            break;
          }
        }
        if (integer) {
          // integer
          changed[numberChanged++] = iColumn;
          solver->setInteger(iColumn);
          if (upper[iColumn] > 1.0e20)
            solver->setColUpper(iColumn, 1.0e20);
          if (lower[iColumn] < -1.0e20)
            solver->setColLower(iColumn, -1.0e20);
        }
      }
    }
    finished = numberChanged == saveNumberChanged;
  }
  delete[] which;
  delete[] changeRhs;
  delete[] ignore; 
  //increment=0.0;
  if (numberChanged) {
    handler_->message(CGL_MADE_INTEGER, messages_)
      << numberChanged
      << CoinMessageEol;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (solver->isInteger(iColumn) && objective[iColumn])
        startModel_->setInteger(iColumn);
    }
  }
  delete solver;
  delete[] changed;
}
//#define BRON_TIMES
#ifdef BRON_TIMES
static int numberTimesX = 0;
#endif
/* Replace cliques by more maximal cliques
   Returns NULL if rows not reduced by greater than cliquesNeeded*rows
   
*/
OsiSolverInterface *
CglPreProcess::cliqueIt(OsiSolverInterface &model,
  double cliquesNeeded) const
{
  /*
    Initial arrays
    * Candidate nodes (columns)
    First nIn already in
    Next nCandidate are candidates
    numberColumns-1 back to nNot are Nots
    * Starts
    * Other node
    * Original row (paired with other node)
    * CliqueIn expanded array with 1 in, 2 not, 3 out, 0 possible, -1 never
    * Type (for original row)
    */
  const double *lower = model.getColLower();
  const double *upper = model.getColUpper();
  const double *rowLower = model.getRowLower();
  const double *rowUpper = model.getRowUpper();
  int numberRows = model.getNumRows();
  int numberColumns = model.getNumCols();
  // Column copy of matrix
  //const double * element = model.getMatrixByCol()->getElements();
  //const int * row = model.getMatrixByCol()->getIndices();
  //const CoinBigIndex * columnStart = model.getMatrixByCol()->getVectorStarts();
  //const int * columnLength = model.getMatrixByCol()->getVectorLengths();
  // Row copy
  CoinPackedMatrix matrixByRow(*model.getMatrixByRow());
  const double *elementByRow = matrixByRow.getElements();
  const int *column = matrixByRow.getIndices();
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();
  char *type = new char[numberRows + 3 * numberColumns];
  char *numberInColumn = type + numberRows;
  char *negativeInColumn = numberInColumn + numberColumns;
  char *positiveInColumn = negativeInColumn + numberColumns;
  memset(numberInColumn, 0, 3 * numberColumns);
  // First pass to mark columns
  int numberCliques = 0;
  //int numberOddCliques=0;
  for (int i = 0; i < numberRows; i++) {
    type[i] = -1;
    double rupper = rowUpper[i];
    double rlower = rowLower[i];
    if (rupper == 1.0 && (rlower <= 0.0 || rlower == 1.0)) {
      bool possible = true;
      CoinBigIndex start = rowStart[i];
      CoinBigIndex end = start + rowLength[i];
      int n = 0;
      for (CoinBigIndex j = start; j < end; j++) {
        int iColumn = column[j];
        if (upper[iColumn] == 1.0 && lower[iColumn] == 0.0 && model.isInteger(iColumn) && elementByRow[j] == 1.0) {
          n++;
        } else {
          possible = false;
          break;
        }
      }
      if (n > 1000)
        possible = false;
      if (possible) {
        for (CoinBigIndex j = start; j < end; j++) {
          int iColumn = column[j];
          if (numberInColumn[iColumn] < 100)
            numberInColumn[iColumn]++;
        }
        numberCliques++;
        if (rowLower[i] > 0.0)
          type[i] = 1;
        else
          type[i] = 0;
      }
    } else if ((rupper == 0.0 || rlower == 0.0) && rowLength[i] == 2) {
      int multiplier;
      if (rupper == 0.0 && rlower < -1.0e20)
        multiplier = 1;
      else if (rlower == 0.0 && rupper > 1.0e20)
        multiplier = -1;
      else
        multiplier = 0;
      if (multiplier) {
        CoinBigIndex start = rowStart[i];
        if (fabs(elementByRow[start]) == 1.0 && fabs(elementByRow[start + 1]) == 1.0) {
          if (elementByRow[start] * elementByRow[start + 1] == -1.0) {
            int iPColumn, iNColumn;
            if (multiplier * elementByRow[start] == 1.0) {
              iPColumn = column[start];
              iNColumn = column[start + 1];
            } else {
              iNColumn = column[start];
              iPColumn = column[start + 1];
            }
            if (upper[iPColumn] == 1.0 && lower[iPColumn] == 0.0 && model.isInteger(iPColumn) && upper[iNColumn] == 1.0 && lower[iNColumn] == 0.0 && model.isInteger(iNColumn)) {
              type[i] = -2;
              if (positiveInColumn[iPColumn] < 100)
                positiveInColumn[iPColumn]++;
              if (negativeInColumn[iNColumn] < 100)
                negativeInColumn[iNColumn]++;
            }
          }
        }
      }
    }
  }
#if CBC_USEFUL_PRINTING > 0
  // look at odd cliques
  int nOdd = 0;
  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (!numberInColumn[iColumn] && negativeInColumn[iColumn] > 1) {
      nOdd++;
    }
  }
  if (nOdd)
    printf("%d possible odd cliques?\n", nOdd);
#endif
  double numberElements = 0;
  if (numberCliques > std::max(1, static_cast< int >(cliquesNeeded * numberRows))) {
    numberCliques = 0;
    for (int i = 0; i < numberRows; i++) {
      if (type[i] >= 0) {
        bool possible = true;
        int n = 0;
        CoinBigIndex start = rowStart[i];
        CoinBigIndex end = start + rowLength[i];
        for (CoinBigIndex j = start; j < end; j++) {
          int iColumn = column[j];
          if (numberInColumn[iColumn] < 2) {
            possible = false;
            break;
          } else {
            n++;
          }
        }
        if (possible) {
          numberElements += n * (n - 1);
          numberCliques++;
        } else {
          type[i] = -1;
        }
      }
    }
  }
  OsiSolverInterface *newSolver = NULL;
  if (numberCliques > std::max(1, static_cast< int >(cliquesNeeded * numberRows))) {
    if (numberElements < 5.0e7 && numberElements < numberCliques * 100) {
#if CBC_USEFUL_PRINTING > 0
      printf("%d cliques needing 2 * %g ints\n",
        numberCliques, numberElements);
#endif
#ifdef BRON_TIMES
      double time1 = CoinCpuTime();
#endif
      CglBK bk(model, type, static_cast< int >(numberElements));
      bk.bronKerbosch();
      newSolver = bk.newSolver(model);
#ifdef BRON_TIMES
      printf("Time %g - bron called %d times\n", CoinCpuTime() - time1, numberTimesX);
#endif
    } else {
#if CBC_USEFUL_PRINTING > 0
      printf("*** %d cliques needing 2 * %g ints\n",
        numberCliques, numberElements);
#endif
    }
  }
  delete[] type;
  return newSolver;
}
// Default constructor
CglBK::CglBK()
{
  candidates_ = NULL;
  mark_ = NULL;
  start_ = NULL;
  otherColumn_ = NULL;
  originalRow_ = NULL;
  dominated_ = NULL;
  cliqueMatrix_ = NULL;
  rowType_ = NULL;
  numberColumns_ = 0;
  numberRows_ = 0;
  numberPossible_ = 0;
  numberCandidates_ = 0;
  firstNot_ = 0;
  numberIn_ = 0;
  left_ = 0;
  lastColumn_ = 0;
}

// Useful constructor
CglBK::CglBK(const OsiSolverInterface &model, const char *rowType,
  int numberElements)
{
  const double *lower = model.getColLower();
  const double *upper = model.getColUpper();
  const double *rowLower = model.getRowLower();
  const double *rowUpper = model.getRowUpper();
  numberRows_ = model.getNumRows();
  numberColumns_ = model.getNumCols();
  // Column copy of matrix
#ifndef NDEBUG
  const double *element = model.getMatrixByCol()->getElements();
#endif
  const int *row = model.getMatrixByCol()->getIndices();
  const CoinBigIndex *columnStart = model.getMatrixByCol()->getVectorStarts();
  const int *columnLength = model.getMatrixByCol()->getVectorLengths();
  start_ = new CoinBigIndex[numberColumns_ + 1];
  otherColumn_ = new int[numberElements];
  candidates_ = new int[2 * numberColumns_];
  CoinZeroN(candidates_, 2 * numberColumns_); // for valgrind
  originalRow_ = new int[numberElements];
  dominated_ = new int[numberRows_];
  CoinZeroN(dominated_, numberRows_);
  //int inputNumber=numberElements;
  numberElements = 0;
  numberPossible_ = 0;
  rowType_ = rowType;
  // Row copy
  CoinPackedMatrix matrixByRow(*model.getMatrixByRow());
  const double *elementByRow = matrixByRow.getElements();
  const int *column = matrixByRow.getIndices();
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();
#if 1
  // take out duplicate doubleton rows
  double *sort = new double[numberRows_];
  int *which = new int[numberRows_];
  double *randomValues = new double[numberColumns_];
  // Initialize random seed
  CoinThreadRandom randomGenerator(987654321);
  for (int i = 0; i < numberColumns_; i++)
    randomValues[i] = randomGenerator.randomDouble();
  int nSort = 0;
  for (int i = 0; i < numberRows_; i++) {
    if (rowLength[i] == 2 && rowUpper[i] == 1.0) {
      CoinBigIndex first = rowStart[i];
      CoinBigIndex last = first + 1;
      if (column[first] > column[last]) {
        first = last;
        last = rowStart[i];
      }
      int iColumn1 = column[first];
      int iColumn2 = column[last];
      double value = elementByRow[first] * randomValues[iColumn1] + elementByRow[last] * randomValues[iColumn2];
      sort[nSort] = value;
      which[nSort++] = i;
    }
  }
  CoinSort_2(sort, sort + nSort, which);
  double value = sort[0];
  //int nDominated = 0;
  for (int i = 1; i < nSort; i++) {
    if (sort[i] == value) {
      int i1 = which[i - 1];
      int i2 = which[i];
      if (rowLower[i1] == rowLower[i2]) {
        CoinBigIndex first1 = rowStart[i1];
        CoinBigIndex last1 = first1 + 1;
        if (column[first1] > column[last1]) {
          first1 = last1;
          last1 = rowStart[i1];
        }
        int iColumn11 = column[first1];
        int iColumn12 = column[last1];
        CoinBigIndex first2 = rowStart[i2];
        CoinBigIndex last2 = first2 + 1;
        if (column[first2] > column[last2]) {
          first2 = last2;
          last2 = rowStart[i2];
        }
        int iColumn21 = column[first2];
        int iColumn22 = column[last2];
        if (iColumn11 == iColumn21 && iColumn12 == iColumn22 && elementByRow[first1] == elementByRow[first2] && elementByRow[last1] == elementByRow[last2]) {
          dominated_[i2] = 1;
          //nDominated++;
        }
      }
    }
    value = sort[i];
  }
  //if (nDominated)
  //printf("%d duplicate doubleton rows!\n",nDominated);
  delete[] randomValues;
  delete[] sort;
  delete[] which;
#endif
  for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
    start_[iColumn] = numberElements;
    CoinBigIndex start = columnStart[iColumn];
    CoinBigIndex end = start + columnLength[iColumn];
    if (upper[iColumn] == 1.0 && lower[iColumn] == 0.0 && model.isInteger(iColumn)) {
      for (CoinBigIndex j = start; j < end; j++) {
        int iRow = row[j];
        if (rowType[iRow] >= 0 && !dominated_[iRow]) {
          assert(element[j] == 1.0);
#if 0
	  CoinBigIndex r=rowStart[iRow];
	  assert (rowLength[iRow]==2);
	  int kColumn = column[r];
	  if (kColumn==iColumn)
	    kColumn=column[r+1];
	  originalRow_[numberElements]=iRow;
	  otherColumn_[numberElements++]=kColumn;
#else
          for (CoinBigIndex r = rowStart[iRow]; r < rowStart[iRow] + rowLength[iRow]; r++) {
            int kColumn = column[r];
            if (kColumn != iColumn) {
              originalRow_[numberElements] = iRow;
              otherColumn_[numberElements++] = kColumn;
            }
          }

#endif
        }
      }
      if (numberElements > start_[iColumn]) {
        candidates_[numberPossible_++] = iColumn;
      }
    }
  }
  //printf("input number %d, computed number %d\n",inputNumber,numberElements);
  start_[numberColumns_] = numberElements;
  numberCandidates_ = numberPossible_;
  numberIn_ = 0;
  firstNot_ = numberPossible_;
  left_ = numberPossible_;
  lastColumn_ = -1;
  mark_ = new char[numberColumns_];
  memset(mark_, 0, numberColumns_);
  cliqueMatrix_ = new CoinPackedMatrix(false, 0.5, 0.0);
  int n = 0;
  for (int i = 0; i < numberRows_; i++) {
    if (rowType[i] >= 0)
      n++;
  }
  cliqueMatrix_->reserve(std::min(100, n), 5 * numberPossible_);
}

// Copy constructor .
CglBK::CglBK(const CglBK &rhs)
{
  // This only copies data in candidates_
  // rest just points
  candidates_ = CoinCopyOfArray(rhs.candidates_, 2 * rhs.numberPossible_);
  mark_ = rhs.mark_;
  start_ = rhs.start_;
  otherColumn_ = rhs.otherColumn_;
  originalRow_ = rhs.originalRow_;
  dominated_ = rhs.dominated_;
  cliqueMatrix_ = rhs.cliqueMatrix_;
  rowType_ = rhs.rowType_;
  numberColumns_ = rhs.numberColumns_;
  numberRows_ = rhs.numberRows_;
  numberPossible_ = rhs.numberPossible_;
  numberCandidates_ = rhs.numberCandidates_;
  firstNot_ = rhs.firstNot_;
  numberIn_ = rhs.numberIn_;
  left_ = rhs.left_;
  lastColumn_ = rhs.lastColumn_;
}

// Assignment operator
CglBK &CglBK::operator=(const CglBK &rhs)
{
  if (this != &rhs) {
    delete[] candidates_;
    // This only copies data in candidates_
    // rest just points
    candidates_ = CoinCopyOfArray(rhs.candidates_, 2 * numberPossible_);
    mark_ = rhs.mark_;
    start_ = rhs.start_;
    otherColumn_ = rhs.otherColumn_;
    originalRow_ = rhs.originalRow_;
    dominated_ = rhs.dominated_;
    cliqueMatrix_ = rhs.cliqueMatrix_;
    rowType_ = rhs.rowType_;
    numberColumns_ = rhs.numberColumns_;
    numberRows_ = rhs.numberRows_;
    numberPossible_ = rhs.numberPossible_;
    numberCandidates_ = rhs.numberCandidates_;
    firstNot_ = rhs.firstNot_;
    numberIn_ = rhs.numberIn_;
    left_ = rhs.left_;
    lastColumn_ = rhs.lastColumn_;
  }
  return *this;
}

// Destructor
CglBK::~CglBK()
{
  delete[] candidates_;
  // only deletes if left_==-1
  if (left_ == -1) {
    delete[] mark_;
    delete[] start_;
    delete[] otherColumn_;
    delete[] originalRow_;
    delete[] dominated_;
    delete cliqueMatrix_;
  }
}
// For Bron-Kerbosch
void CglBK::bronKerbosch()
{
#ifdef BRON_TIMES
  numberTimesX++;
  if ((numberTimesX % 1000) == 0)
    printf("times %d - %d candidates left\n", numberTimesX, numberCandidates_);
#endif
  if (!numberCandidates_ && firstNot_ == numberPossible_) {
    // mark original rows which are dominated
    // save if clique size >2
    if (numberIn_ > 2) {
      double *elements = new double[numberIn_];
      int *column = candidates_ + numberPossible_;
      // mark those in clique
      for (int i = 0; i < numberIn_; i++) {
        int iColumn = column[i];
        mark_[iColumn] = 1;
      }
      for (int i = 0; i < numberIn_; i++) {
        elements[i] = 1.0;
        int iColumn = column[i];
        for (CoinBigIndex j = start_[iColumn]; j < start_[iColumn + 1]; j++) {
          int jColumn = otherColumn_[j];
          if (mark_[jColumn]) {
            int iRow = originalRow_[j];
            /* only get rid of <= cliques
	       can we find dominating clique??
	       should really look further if dominating clique also == 
	       or could make that == */
            if (rowType_[iRow] == 0)
              dominated_[iRow]++;
          }
        }
      }
      for (int i = 0; i < numberIn_; i++) {
        int iColumn = column[i];
        mark_[iColumn] = 0;
      }
      cliqueMatrix_->appendRow(numberIn_, column, elements);
      delete[] elements;
    }
  } else {
#if 0
    int nCplusN=numberCandidates_+(numberPossible_-firstNot_);
    int iChoose = CoinDrand48()*nCplusN;
    iChoose=std::min(0,nCplusN-1);
    if (iChoose>=numberCandidates_) {
      iChoose -= numberCandidates_;
      iChoose += firstNot_;
    }
#else
    for (int i = 0; i < numberCandidates_; i++) {
      int jColumn = candidates_[i];
      mark_[jColumn] = 1;
    }
    int nMax = 0;
    int iChoose = 0;
    for (int i = numberPossible_ - 1; i >= firstNot_; i--) {
      int iColumn = candidates_[i];
      int n = 0;
      for (CoinBigIndex j = start_[iColumn]; j < start_[iColumn + 1]; j++) {
        int jColumn = otherColumn_[j];
        n += mark_[jColumn];
      }
      if (n > nMax) {
        nMax = n;
        iChoose = i;
      }
    }
    if (nMax < numberCandidates_ - 1 || !nMax) {
      for (int i = 0; i < numberCandidates_; i++) {
        int iColumn = candidates_[i];
        int n = 0;
        for (CoinBigIndex j = start_[iColumn]; j < start_[iColumn + 1]; j++) {
          int jColumn = otherColumn_[j];
          n += mark_[jColumn];
        }
        if (n > nMax) {
          nMax = n;
          iChoose = i;
        }
      }
    }
    for (int i = 0; i < numberCandidates_; i++) {
      int jColumn = candidates_[i];
      mark_[jColumn] = 0;
    }
#endif
    iChoose = candidates_[iChoose];
    int *temp = candidates_ + numberPossible_ + numberIn_;
    int nTemp = 0;
    if (nMax < numberCandidates_) {
      // Neighborhood of iColumn
      for (CoinBigIndex j = start_[iChoose]; j < start_[iChoose + 1]; j++) {
        int jColumn = otherColumn_[j];
        mark_[jColumn] = 1;
      }
      for (int i = 0; i < numberCandidates_; i++) {
        int jColumn = candidates_[i];
        if (!mark_[jColumn])
          temp[nTemp++] = jColumn;
      }
      for (CoinBigIndex j = start_[iChoose]; j < start_[iChoose + 1]; j++) {
        int jColumn = otherColumn_[j];
        mark_[jColumn] = 0;
      }
    }
    //if (nMax==numberCandidates_)
    //assert (!nTemp);
    for (int kk = 0; kk < nTemp; kk++) {
      int iColumn = temp[kk];
      // move up
      int put = 0;
      for (int i = 0; i < numberCandidates_; i++) {
        if (candidates_[i] != iColumn)
          candidates_[put++] = candidates_[i];
      }
      numberCandidates_--;
      CglBK bk2(*this);
      int *newCandidates = bk2.candidates_;
#if 0
      printf("%p (next %p) iColumn %d, %d candidates %d not %d in\n",
	     this,&bk2,iColumn,numberCandidates_,
	     numberPossible_-firstNot_,numberIn_);
      for (int i=0;i<numberCandidates_;i++) {
	printf(" %d",candidates_[i]);
      }
      printf("\n");
#endif
      newCandidates[numberPossible_ + numberIn_] = iColumn;
      bk2.numberIn_ = numberIn_ + 1;
      // Neighborhood of iColumn
      for (CoinBigIndex j = start_[iColumn]; j < start_[iColumn + 1]; j++) {
        int jColumn = otherColumn_[j];
        mark_[jColumn] = 1;
      }
      // Intersection of candidates with neighborhood
      int numberCandidates = 0;
      for (int i = 0; i < bk2.numberCandidates_; i++) {
        int jColumn = newCandidates[i];
        if (mark_[jColumn])
          newCandidates[numberCandidates++] = jColumn;
      }
      bk2.numberCandidates_ = numberCandidates;
      // Intersection of not with neighborhood
      int firstNot = numberPossible_;
      for (int i = numberPossible_ - 1; i >= bk2.firstNot_; i--) {
        int jColumn = newCandidates[i];
        if (mark_[jColumn])
          newCandidates[--firstNot] = jColumn;
      }
      bk2.firstNot_ = firstNot;
      for (CoinBigIndex j = start_[iColumn]; j < start_[iColumn + 1]; j++) {
        int jColumn = otherColumn_[j];
        mark_[jColumn] = 0;
      }
      //for (int i=0;i<numberColumns_;i++)
      //assert (!mark_[i]);
      bk2.bronKerbosch();
      candidates_[--firstNot_] = iColumn;
    }
  }
}
// Creates strengthened smaller model
OsiSolverInterface *
CglBK::newSolver(const OsiSolverInterface &model)
{
  // See how many rows can be deleted
  int nDelete = 0;
  int *deleted = new int[numberRows_];
  for (int i = 0; i < numberRows_; i++) {
    if (dominated_[i]) {
      deleted[nDelete++] = i;
    }
  }
  int nAdd = cliqueMatrix_->getNumRows();
#if CBC_USEFUL_PRINTING > 0
  printf("%d rows can be deleted with %d new cliques\n",
    nDelete, nAdd);
#endif

  OsiSolverInterface *newSolver = NULL;
  if (nDelete > nAdd) {
    newSolver = model.clone();
    newSolver->deleteRows(nDelete, deleted);
    double *lower = new double[nAdd];
    double *upper = new double[nAdd];
    for (int i = 0; i < nAdd; i++) {
      lower[i] = -COIN_DBL_MAX;
      upper[i] = 1.0;
    }
    const double *elementByRow = cliqueMatrix_->getElements();
    const int *column = cliqueMatrix_->getIndices();
    const CoinBigIndex *rowStart = cliqueMatrix_->getVectorStarts();
    //const int * rowLength = cliqueMatrix_->getVectorLengths();
    assert(cliqueMatrix_->getNumElements() == rowStart[nAdd]);
    newSolver->addRows(nAdd, rowStart, column, elementByRow, lower, upper);
#if PRINT_DEBUG
    for (int i = 0; i < nAdd; i++) {
      if (rowStart[i + 1] - rowStart[i] > 10)
        printf("Clique %d has %d entries\n", i, rowStart[i + 1] - rowStart[i]);
    }
#endif
    delete[] lower;
    delete[] upper;
  }
  delete[] deleted;
  // mark so everything will be deleted
  left_ = -1;
  return newSolver;
}
CglUniqueRowCuts::CglUniqueRowCuts(int initialMaxSize, int hashMultiplier)
{
  numberCuts_ = 0;
  size_ = initialMaxSize;
  hashMultiplier_ = hashMultiplier;
  int hashSize = hashMultiplier_ * size_;
  if (size_) {
    rowCut_ = new OsiRowCut *[size_];
    hash_ = new CglHashLink[hashSize];
  } else {
    rowCut_ = NULL;
    hash_ = NULL;
  }
  for (int i = 0; i < hashSize; i++) {
    hash_[i].index = -1;
    hash_[i].next = -1;
  }
  lastHash_ = -1;
}
CglUniqueRowCuts::~CglUniqueRowCuts()
{
  for (int i = 0; i < numberCuts_; i++)
    delete rowCut_[i];
  delete[] rowCut_;
  delete[] hash_;
}
CglUniqueRowCuts::CglUniqueRowCuts(const CglUniqueRowCuts &rhs)
{
  numberCuts_ = rhs.numberCuts_;
  hashMultiplier_ = rhs.hashMultiplier_;
  size_ = rhs.size_;
  int hashSize = size_ * hashMultiplier_;
  lastHash_ = rhs.lastHash_;
  if (size_) {
    rowCut_ = new OsiRowCut *[size_];
    hash_ = new CglHashLink[hashSize];
    for (int i = 0; i < hashSize; i++) {
      hash_[i] = rhs.hash_[i];
    }
    for (int i = 0; i < size_; i++) {
      if (rhs.rowCut_[i])
        rowCut_[i] = new OsiRowCut(*rhs.rowCut_[i]);
      else
        rowCut_[i] = NULL;
    }
  } else {
    rowCut_ = NULL;
    hash_ = NULL;
  }
}
CglUniqueRowCuts &
CglUniqueRowCuts::operator=(const CglUniqueRowCuts &rhs)
{
  if (this != &rhs) {
    for (int i = 0; i < numberCuts_; i++)
      delete rowCut_[i];
    delete[] rowCut_;
    delete[] hash_;
    numberCuts_ = rhs.numberCuts_;
    hashMultiplier_ = rhs.hashMultiplier_;
    size_ = rhs.size_;
    lastHash_ = rhs.lastHash_;
    if (size_) {
      rowCut_ = new OsiRowCut *[size_];
      int hashSize = size_ * hashMultiplier_;
      hash_ = new CglHashLink[hashSize];
      for (int i = 0; i < hashSize; i++) {
        hash_[i] = rhs.hash_[i];
      }
      for (int i = 0; i < size_; i++) {
        if (rhs.rowCut_[i])
          rowCut_[i] = new OsiRowCut(*rhs.rowCut_[i]);
        else
          rowCut_[i] = NULL;
      }
    } else {
      rowCut_ = NULL;
      hash_ = NULL;
    }
  }
  return *this;
}
void CglUniqueRowCuts::eraseRowCut(int sequence)
{
  // find
  assert(sequence >= 0 && sequence < numberCuts_);
  OsiRowCut *cut = rowCut_[sequence];
  int hashSize = size_ * hashMultiplier_;
  int ipos = hashCut(*cut, hashSize);
  int found = -1;
  while (true) {
    int j1 = hash_[ipos].index;
    if (j1 >= 0) {
      if (j1 != sequence) {
        int k = hash_[ipos].next;
        if (k != -1)
          ipos = k;
        else
          break;
      } else {
        found = j1;
        break;
      }
    } else {
      break;
    }
  }
  assert(found >= 0);
  assert(hash_[ipos].index == sequence);
  // shuffle up
  while (hash_[ipos].next >= 0) {
    int k = hash_[ipos].next;
    hash_[ipos] = hash_[k];
    ipos = k;
  }
  delete cut;
  // move last to found
  numberCuts_--;
  if (numberCuts_) {
    ipos = hashCut(*rowCut_[numberCuts_], hashSize);
    while (true) {
      int j1 = hash_[ipos].index;
      if (j1 != numberCuts_) {
        int k = hash_[ipos].next;
        ipos = k;
      } else {
        // change
        hash_[ipos].index = found;
        rowCut_[found] = rowCut_[numberCuts_];
        rowCut_[numberCuts_] = NULL;
        break;
      }
    }
  }
  assert(!rowCut_[numberCuts_]);
}
// Return 0 if added, 1 if not
int CglUniqueRowCuts::insertIfNotDuplicate(const OsiRowCut &cut)
{
  int hashSize = size_ * hashMultiplier_;
  if (numberCuts_ == size_) {
    size_ = 2 * size_ + 100;
    hashSize = hashMultiplier_ * size_;
    OsiRowCut **temp = new OsiRowCut *[size_];
    delete[] hash_;
    hash_ = new CglHashLink[hashSize];
    for (int i = 0; i < hashSize; i++) {
      hash_[i].index = -1;
      hash_[i].next = -1;
    }
    for (int i = 0; i < numberCuts_; i++) {
      temp[i] = rowCut_[i];
      int ipos = hashCut(*temp[i], hashSize);
      int found = -1;
      int jpos = ipos;
      while (true) {
        int j1 = hash_[ipos].index;

        if (j1 >= 0) {
          if (!same(*temp[i], *temp[j1])) {
            int k = hash_[ipos].next;
            if (k != -1)
              ipos = k;
            else
              break;
          } else {
            found = j1;
            break;
          }
        } else {
          break;
        }
      }
      if (found < 0) {
        assert(hash_[ipos].next == -1);
        if (ipos == jpos) {
          // first
          hash_[ipos].index = i;
        } else {
          // find next space
          while (true) {
            ++lastHash_;
            assert(lastHash_ < hashSize);
            if (hash_[lastHash_].index == -1)
              break;
          }
          hash_[ipos].next = lastHash_;
          hash_[lastHash_].index = i;
        }
      }
    }
    delete[] rowCut_;
    rowCut_ = temp;
  }
  if (numberCuts_ < size_) {
    double newLb = cut.lb();
    double newUb = cut.ub();
    CoinPackedVector vector = cut.row();
    int numberElements = vector.getNumElements();
    int *newIndices = vector.getIndices();
    double *newElements = vector.getElements();
    CoinSort_2(newIndices, newIndices + numberElements, newElements);
    int i;
    bool bad = false;
    for (i = 0; i < numberElements; i++) {
      double value = fabs(newElements[i]);
      if (value < 1.0e-12 || value > 1.0e12)
        bad = true;
    }
    if (bad)
      return 1;
    OsiRowCut newCut;
    newCut.setLb(newLb);
    newCut.setUb(newUb);
    newCut.setRow(vector);
    int ipos = hashCut(newCut, hashSize);
    int found = -1;
    int jpos = ipos;
    while (true) {
      int j1 = hash_[ipos].index;

      if (j1 >= 0) {
        if (!same(newCut, *rowCut_[j1])) {
          int k = hash_[ipos].next;
          if (k != -1)
            ipos = k;
          else
            break;
        } else {
          found = j1;
          break;
        }
      } else {
        break;
      }
    }
    if (found < 0) {
      assert(hash_[ipos].next == -1);
      if (ipos == jpos) {
        // first
        hash_[ipos].index = numberCuts_;
      } else {
        // find next space
        while (true) {
          ++lastHash_;
          assert(lastHash_ < hashSize);
          if (hash_[lastHash_].index == -1)
            break;
        }
        hash_[ipos].next = lastHash_;
        hash_[lastHash_].index = numberCuts_;
      }
      OsiRowCut *newCutPtr = new OsiRowCut();
      newCutPtr->setLb(newLb);
      newCutPtr->setUb(newUb);
      newCutPtr->setRow(vector);
      rowCut_[numberCuts_++] = newCutPtr;
      return 0;
    } else {
      return 1;
    }
  } else {
    return -1;
  }
}
// Add in cuts as normal cuts and delete
void CglUniqueRowCuts::addCuts(OsiCuts &cs)
{
  for (int i = 0; i < numberCuts_; i++) {
    cs.insertIfNotDuplicate(*rowCut_[i]);
    delete rowCut_[i];
    rowCut_[i] = NULL;
  }
  numberCuts_ = 0;
}

void CglPreProcess::setTimeLimit(const double timeLimit, const bool useElapsedTime)
{
  this->timeLimit_ = timeLimit;
  this->useElapsedTime_ = useElapsedTime;
}

void CglPreProcess::setKeepColumnNames(const bool keep)
{
  this->keepColumnNames_ = keep;
}

double CglPreProcess::getCurrentCPUTime() const
{
  return CoinGetTimeOfDay();
}
#if CBC_USE_PAPILO
static papiloStruct papiloPresolve(ClpSimplex * inModel,
				   bool treatAsContinuous,
				   double timeLimit)
{
  // move to papilo
  papilo::Problem<double> *papiloProb = new papilo::Problem<double>;
  papiloStruct pproblem = initialTry;
  papilo::Presolve<double> presolve;
  {
    pproblem.originalModel = inModel;
    int numberColumns = inModel->numberColumns();
    int numberRows = inModel->numberRows();
    const double * objective = inModel->objective();
    double offset = inModel->objectiveOffset();
    papiloProb->setObjective(
			     std::vector <double> (objective,objective+numberColumns),-offset);
    CoinPackedMatrix rowCopy;
    rowCopy.reverseOrderedCopyOf(*inModel->matrix());
    int *column = rowCopy.getMutableIndices();
    CoinBigIndex *rowStart = rowCopy.getMutableVectorStarts();
    double *element = rowCopy.getMutableElements();
    int numberElements = rowCopy.getNumElements();
    double spareRatio = 2.0; 
    int gap = 4; 
    if (numberElements>10000) {
      spareRatio = 1.2;
      gap = 2;
    }
    papilo::SparseStorage<double> copy(numberRows,
				       numberColumns,numberElements,
				       spareRatio,gap);
    // build storage
    int * columns = copy.getColumns();
    //papilo::Vec<int>  rowstart = copy.getRowStarts();
    papilo::IndexRange * rowranges = copy.getRowRanges();;
    double * values = copy.getValues();
    int shift = 0;
    for( int r = 0; r < numberRows; r++ ) {
      rowranges[r].start = rowStart[r] + shift;
      
      for( int j = rowStart[r]; j < rowStart[r + 1]; j++ ) {
	assert( j + shift >= 0 );
	values[j + shift] = element[j];
	columns[j + shift] = column[j];
      }
      rowranges[r].end = rowStart[r + 1] + shift;
      const int rowsize = rowranges[r].end - rowranges[r].start;
      const int rowalloc = copy.computeRowAlloc( rowsize );
      shift += rowalloc - rowsize;
    }
    rowranges[numberRows].start = rowStart[numberRows] + shift;
    rowranges[numberRows].end = rowranges[numberRows].start;
    const double * rLower = inModel->rowLower();
    const double * rUpper = inModel->rowUpper();
    const double * cLower = inModel->columnLower();
    const double * cUpper = inModel->columnUpper();
    std::vector <double> rlower(rLower,rLower+numberRows);
    std::vector <double> rupper(rUpper,rUpper+numberRows);
    std::vector <papilo::RowFlags> flags(numberRows);
    for (int i=0;i<numberRows;i++) {
      if (rLower[i]<-1.0e30) {
	if (rUpper[i]<1.0e30)
	  flags[i] = papilo::RowFlag::kLhsInf;
	else
	  flags[i] = papilo::RowFlag::kRedundant;
      } else if (rUpper[i]>1.0e30) {
	flags[i] = papilo::RowFlag::kRhsInf;
      } else if (rUpper[i]==rLower[i]) {
	flags[i] = papilo::RowFlag::kEquation;
      }
    }
    papiloProb->setConstraintMatrix(copy,rlower,rupper,flags);
    std::vector <papilo::ColFlags> cflags(numberColumns);
    // may wish to treat as continuous (parameter)
    bool treatAsContinuous = false;
    int nInt = 0;
    int nBinary = 0;
    for (int i=0;i<numberColumns;i++) {
      unsigned int flag = static_cast<unsigned int>(papilo::ColFlag::kNone);
      if (inModel->isInteger(i) && !treatAsContinuous) { 
	flag |= static_cast<unsigned int>(papilo::ColFlag::kIntegral);
	nInt++;
	if (!cLower[i] && cUpper[i] ==1.0)
	  nBinary++;
      }
      if (cLower[i]<-1.0e100)
	flag |= static_cast<unsigned int>(papilo::ColFlag::kLbInf);
      else if (cLower[i]<-1.0e20)
	flag = static_cast<unsigned int>(papilo::ColFlag::kLbHuge);
      if (cUpper[i]>1.0e100)
	flag |= static_cast<unsigned int>(papilo::ColFlag::kUbInf);
      else if (cUpper[i]>1.0e20)
	flag = static_cast<unsigned int>(papilo::ColFlag::kUbHuge);
      //if (cLower[i] == cUpper[i])
      //flag |= static_cast<unsigned int>(papilo::ColFlag::kFixed);
      cflags[i] = static_cast<papilo::ColFlag>(flag);
    }
    unsigned int pflag = 0;
    if (!nInt) {
      pflag |= static_cast<unsigned int>(papilo::ProblemFlag::kLinear);
    } else {
      if (nInt<numberColumns) {
	pflag |= static_cast<unsigned int>(papilo::ProblemFlag::kMixedInteger);
      } else if (nInt==nBinary) {
	pflag |= static_cast<unsigned int>(papilo::ProblemFlag::kBinary);
      } else {
	pflag |= static_cast<unsigned int>(papilo::ProblemFlag::kInteger);
      } 
    }
    papiloProb->set_problem_type(static_cast<papilo::ProblemFlag>(pflag));
    std::vector <double> clower(cLower,cLower+numberColumns); 
    std::vector <double> cupper(cUpper,cUpper+numberColumns);
    papiloProb->setVariableDomains(clower,cupper,cflags);
    presolve.addDefaultPresolvers();
    presolve.getPresolveOptions().threads = initialTry.numberThreads;
    presolve.setVerbosityLevel(papilo::VerbosityLevel::kWarning);
    if (inModel->logLevel())
      presolve.setVerbosityLevel(papilo::VerbosityLevel::kInfo);
    // keep sparser even if more rows
    presolve.getPresolveOptions().maxfillinpersubstitution = 3;
    presolve.getPresolveOptions().tlim = timeLimit;
  }
  presolveResult = presolve.apply( *papiloProb, false );
  switch (presolveResult.status) {
  case papilo::PresolveStatus::kUnchanged:
    delete papiloProb;
    pproblem.papiloProb = NULL;
    pproblem.presolvedModel = NULL;
    pproblem.returnCode = 0;
    return pproblem;
  case papilo::PresolveStatus::kReduced:
    pproblem.papiloProb = papiloProb;
    pproblem.presolvedModel = new ClpSimplex(*inModel);
    pproblem.returnCode = 1;
    break;
  case papilo::PresolveStatus::kUnbndOrInfeas:
  case papilo::PresolveStatus::kUnbounded:
  case papilo::PresolveStatus::kInfeasible:
    delete papiloProb;
    pproblem.papiloProb = NULL;
    pproblem.presolvedModel = NULL;
    pproblem.returnCode = -1;
    return pproblem;
  }
  ClpSimplex *model2 = pproblem.presolvedModel;
  {
    // create model2
    int numberColumns2 = papiloProb->getNCols();
    int numberRows2 = papiloProb->getNRows();
    papilo::ConstraintMatrix<double> matrix2 =
      papiloProb->getConstraintMatrix();
    papilo::Objective<double> obj2 = papiloProb->getObjective();
    papilo::Vec<double> lower2 = papiloProb->getLowerBounds();
    papilo::Vec<double> upper2 = papiloProb->getUpperBounds();
    // what about symmetry? - getSymmetries()
    CoinBigIndex nnz = matrix2.getNnz();
    CoinBigIndex * starts = new CoinBigIndex [numberColumns2+1];
    int * row = new int [nnz];
    double * element = new double [nnz];
    nnz = 0;
    starts[0] = 0;
    for (int i=0;i<numberColumns2;i++) {
      papilo::SparseVectorView<double> column =
	matrix2.getColumnCoefficients(i);
      int n = column.getLength();
      memcpy(element+nnz,column.getValues(),n*sizeof(double));
      memcpy(row+nnz,column.getIndices(),n*sizeof(int));
      nnz += n;
      starts[i+1] = nnz;
    }
    double * cLower = new double[3*numberColumns2];
    double * cUpper = cLower+numberColumns2;
    double * obj = cUpper+numberColumns2;
    const papilo::Vec<double> obj2p = obj2.coefficients;
    for (int i=0;i<numberColumns2;i++) {
      cLower[i] = lower2[i];
      cUpper[i] = upper2[i];
      assert (!papiloProb->getColFlags()[i].test( papilo::ColFlag::kFixed ));
      obj[i] = obj2p[i];
    }
    double offset = papiloProb->getObjective().offset;
    double * rLower = new double[2*numberRows2];
    double * rUpper = rLower+numberRows2;
    lower2 = matrix2.getLeftHandSides();
    upper2 = matrix2.getRightHandSides();
    for (int i=0;i<numberRows2;i++) {
      rLower[i] = lower2[i];
      rUpper[i] = upper2[i];
      if (papiloProb->getRowFlags()[i].test( papilo::RowFlag::kLhsInf ))
	rLower[i] = -COIN_DBL_MAX;
      else if (papiloProb->getRowFlags()[i].test( papilo::RowFlag::kRhsInf ))
	rUpper[i] = COIN_DBL_MAX;
      else if (papiloProb->getRowFlags()[i].test( papilo::RowFlag::kEquation ))
	assert (rLower[i] == rUpper[i]);
      // what do kIntegral (and kRedundant) mean
    }
    model2->loadProblem(numberColumns2,numberRows2,starts,row,element,
		       cLower,cUpper,obj,rLower,rUpper);
    delete [] starts;
    delete [] rLower;
    delete [] row;
    delete [] element;
    delete [] cLower;
    model2->setObjectiveOffset(-papiloProb->getObjective().offset);
    std::vector <int> mapping = presolveResult.postsolve.origcol_mapping;
    std::string columnName;
    int *mapping2 = new int [numberColumns2];
    pproblem.mapping = mapping2;
    for (int i=0;i<numberColumns2;i++) {
      int iColumn = mapping[i]; 
      mapping2[i] = iColumn;
      columnName = inModel->columnName(iColumn);
      model2->setColumnName(i,columnName);
      if (inModel->isInteger(iColumn))
	model2->setInteger(i);
      else
	model2->setContinuous(i);
    }
  }
  return pproblem;
}
static ClpSimplex *postSolvedModel(papiloStruct papilo)
{
  ClpSimplex * simplex = papilo.presolvedModel;
  papilo::Problem<double> *papiloProb =
    static_cast<papilo::Problem<double> *>(papilo.papiloProb);
  int numberRows = simplex->numberRows();
  int numberColumns = simplex->numberColumns();
  double * rowSol = simplex->primalRowSolution();
  double * colSol = simplex->primalColumnSolution();
  double * duals = simplex->dualRowSolution();
  double * djs = simplex->dualColumnSolution();
  papilo::Solution<double> solution;
  for (int i=0;i<numberColumns;i++) {
    solution.primal.push_back(colSol[i]);
    solution.reducedCosts.push_back(djs[i]);
    papilo::VarBasisStatus statusP = papilo::VarBasisStatus::UNDEFINED;
    ClpSimplex::Status status = simplex->getStatus(i);
    if (status==ClpSimplex::Status::basic)
      statusP = papilo::VarBasisStatus::BASIC;
    else if (status==ClpSimplex::Status::atLowerBound)
      statusP = papilo::VarBasisStatus::ON_LOWER;
    else if (status==ClpSimplex::Status::atUpperBound)
      statusP = papilo::VarBasisStatus::ON_UPPER;
    solution.varBasisStatus.push_back(statusP);
  }
  for (int i=0;i<numberRows;i++) {
    solution.dual.push_back(duals[i]);
    solution.slack.push_back(rowSol[i]); // probably not right
    papilo::VarBasisStatus statusP = papilo::VarBasisStatus::UNDEFINED;
    ClpSimplex::Status status = simplex->getRowStatus(i);
    if (status==ClpSimplex::Status::basic)
      statusP = papilo::VarBasisStatus::BASIC;
    else if (status==ClpSimplex::Status::atLowerBound)
      statusP = papilo::VarBasisStatus::ON_LOWER; // maybe wrong way round
    else if (status==ClpSimplex::Status::atUpperBound)
      statusP = papilo::VarBasisStatus::ON_UPPER; // maybe wrong way round
    solution.rowBasisStatus.push_back(statusP);
  }
  solution.basisAvailabe = true;
  solution.type = papilo::SolutionType::kPrimal;
  ClpSimplex *inModel = papilo.originalModel;
  int numberColumns2 = inModel->numberColumns();
  int numberRows2 = inModel->numberRows();
  const papilo::Message msg{};
  papilo::Postsolve postsolve(msg,presolveResult.postsolve.getNum());
  papilo::Solution<double> solutionOut;
  solutionOut.type = papilo::SolutionType::kPrimal;
  solutionOut.basisAvailabe = true;
  postsolve.undo(solution,solutionOut,presolveResult.postsolve);
  // just fix
  double * cLower = inModel->columnLower();
  double * cUpper = inModel->columnUpper();
  double * solutionClp = inModel->primalColumnSolution();
  for (int i=0;i<numberColumns2;i++) {
    double value = solutionOut.primal[i];
    if (inModel->isInteger(i)) {
      double value2 = floor(value+0.5);
      assert (fabs(value-value2)<1.0e-5);
      value = value2;
      solutionClp[i] = value;
      cLower[i] = value;
      cUpper[i] = value;
    } else {
      solutionClp[i] = value;
    }
  }
  inModel->primal();
  delete [] initialTry.mapping;
  delete papiloProb;
  return inModel;
}
void zapPapilo(int pOptions,CglPreProcess * process)
{
  memset(&initialTry,0,sizeof(papiloStruct));
  initialTry.numberThreads=pOptions&7;
  initialTry.presolveType = pOptions/8;
  if (initialTry.presolveType)
    initialTry.preProcess = process;
}
#endif
/*
  PAPILO.
  This is a free presolver.  Not much work has been done recently in Coin on
  pre-processing for integer problems and this does well on some problems.
  It is a bit slow on some problems, but it may be worth trying.
  To use - download papilo.  You may wish to download the whole SCIP
  optimization suite.  SCIP has some very good features.  One I like and
  would be pleased if a student could investigate is that they have developed
  a reasonable estimate of how long a branch and bound run is going to take.
  If after some time it looks as if it is going to be very long, it restarts
  (which may result in a smaller problem) and then uses all the information
  it had gathered to do a better job on organizing the tree.

  The instructions below are if you have downloaded to scipoptsuite-9.1.0.
  Simplest is to build everything - look at README.md.  The instructions create 
  a directory "build".  If you have done this then build Cbc and add 
  -DCBC_USE_PAPILO=1
  and 
  "-I ~/scipoptsuite-9.1.0/papilo/src/ -I ~/scipoptsuite-9.1.0/build/papilo" 

  You may need to point to extra libraries although I don't need to and I am
  not at all sure I need clusol at all so e.g.
  LDFLAGS="-lopenblas -ltbb -lclusol -lpthread -lz -lrt" 
  You probably can build without some of this but I have not bothered.

  To use in standalone solver (or one that uses that code)  
  do -preprocess papilo. For full list -preprocess??

  The default is to do normal Cgl preprocessing and then do papilo. You
  can also do papiloBegin which does papilo and then rest of Cgl
  preprocessing.  In the majority of cases you won't get much from
  using papilo but it is worth seeing what happens and it may prompt some
  further work.
*/
