/* $Id: CbcSolverExpandKnapsack.cpp 1240 2009-10-02 18:41:44Z forrest $ */
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

OsiSolverInterface *
expandKnapsack(CoinModel & model, int * whichColumn, int * knapsackStart,
               int * knapsackRow, int &numberKnapsack,
               CglStored & stored, int logLevel,
               int fixedPriority, int SOSPriority, CoinModel & tightenedModel)
{
    int maxTotal = numberKnapsack;
    // load from coin model
    OsiSolverLink *si = new OsiSolverLink();
    OsiSolverInterface * finalModel = NULL;
    si->setDefaultMeshSize(0.001);
    // need some relative granularity
    si->setDefaultBound(100.0);
    si->setDefaultMeshSize(0.01);
    si->setDefaultBound(100000.0);
    si->setIntegerPriority(1000);
    si->setBiLinearPriority(10000);
    si->load(model, true, logLevel);
    // get priorities
    const int * priorities = model.priorities();
    int numberColumns = model.numberColumns();
    if (priorities) {
        OsiObject ** objects = si->objects();
        int numberObjects = si->numberObjects();
        for (int iObj = 0; iObj < numberObjects; iObj++) {
            int iColumn = objects[iObj]->columnNumber();
            if (iColumn >= 0 && iColumn < numberColumns) {
#ifndef NDEBUG
                OsiSimpleInteger * obj =
                    dynamic_cast <OsiSimpleInteger *>(objects[iObj]) ;
#endif
                assert (obj);
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
    int * whichKnapsack = new int [numberColumns];
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
                    const char *  el = coinModel.getElementAsString(iRow, iColumn);
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
                const char *  el = coinModel.getElementAsString(iRow, iColumn);
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
    int * markRow = new int [numberRows];
    for (iRow = 0; iRow < numberRows; iRow++)
        markRow[iRow] = -1;
    int canDo = 1; // OK and linear
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        CoinModelLink triple = coinModel.firstInColumn(iColumn);
        int iKnapsack = whichKnapsack[iColumn];
        bool linear = true;
        // See if quadratic objective
        const char * expr = coinModel.getColumnObjectiveAsString(iColumn);
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
            const char * expr = coinModel.getElementAsString(iRow, iColumn);
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
    int * markKnapsack = NULL;
    double * coefficient = NULL;
    double * linear = NULL;
    int * whichRow = NULL;
    int * lookupRow = NULL;
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
        markKnapsack = new int [numberKnapsack];
        coefficient = new double [numberKnapsack];
        linear = new double [numberColumns];
        for (iKnapsack = 0; iKnapsack < numberKnapsack; iKnapsack++)
            markKnapsack[iKnapsack] = -1;
        if (canDo == 2) {
            for (iRow = -1; iRow < numberRows; iRow++) {
                int numberOdd;
                CoinPackedMatrix * row = coinModel.quadraticRow(iRow, linear, numberOdd);
                if (row) {
                    // see if valid
                    const double * element = row->getElements();
                    const int * column = row->getIndices();
                    const CoinBigIndex * columnStart = row->getVectorStarts();
                    const int * columnLength = row->getVectorLengths();
                    int numberLook = row->getNumCols();
                    for (int i = 0; i < numberLook; i++) {
                        int iKnapsack = whichKnapsack[i];
                        if (iKnapsack < 0) {
                            // might be able to swap - but for now can't have knapsack in
                            for (int j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
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
                            for (int j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
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
            whichRow = new int [numberRows];
            lookupRow = new int [numberRows];
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
            const int * priorities = model.priorities();
            if (priorities) {
                finalModel->findIntegers(false);
                OsiObject ** objects = finalModel->objects();
                int numberObjects = finalModel->numberObjects();
                for (int iObj = 0; iObj < numberObjects; iObj++) {
                    int iColumn = objects[iObj]->columnNumber();
                    if (iColumn >= 0 && iColumn < nCol) {
#ifndef NDEBUG
                        OsiSimpleInteger * obj =
                            dynamic_cast <OsiSimpleInteger *>(objects[iObj]) ;
#endif
                        assert (obj);
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
                int nelCreate = coinModel.expandKnapsack(iRow, nCreate, NULL, NULL, NULL, NULL);
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
                double * buildObj = new double [nLargest];
                double * buildElement = new double [nelLargest];
                int * buildStart = new int[nLargest+1];
                int * buildRow = new int[nelLargest];
                // alow for integers in knapsacks
                OsiObject ** object = new OsiObject * [numberKnapsack+nTotal];
                int nSOS = 0;
                int nObj = numberKnapsack;
                for (iKnapsack = 0; iKnapsack < numberKnapsack; iKnapsack++) {
                    knapsackStart[iKnapsack] = finalModel->getNumCols();
                    iRow = knapsackRow[iKnapsack];
                    int nCreate = 10000;
                    coinModel.expandKnapsack(iRow, nCreate, buildObj, buildStart, buildRow, buildElement);
                    // Redo row numbers
                    for (iColumn = 0; iColumn < nCreate; iColumn++) {
                        for (int j = buildStart[iColumn]; j < buildStart[iColumn+1]; j++) {
                            int jRow = buildRow[j];
                            jRow = lookupRow[jRow];
                            assert (jRow >= 0 && jRow < nRow);
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
                        OsiSimpleInteger * sosObject = new OsiSimpleInteger(finalModel, iColumn);
                        sosObject->setPriority(1000000);
                        object[nObj++] = sosObject;
                        buildRow[iColumn-numberOther] = iColumn;
                        buildElement[iColumn-numberOther] = 1.0;
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
                        OsiSOS * sosObject = new OsiSOS(finalModel, n - 1, buildRow, NULL, 1);
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
                delete [] object;
                // Can we move any rows to cuts
                const int * cutMarker = coinModel.cutMarker();
                if (cutMarker && 0) {
                    printf("AMPL CUTS OFF until global cuts fixed\n");
                    cutMarker = NULL;
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
                    int nDelete = 0;
                    for (iRow = 0; iRow < numberRows; iRow++) {
                        if (cutMarker[iRow] && lookupRow[iRow] >= 0) {
                            int jRow = lookupRow[iRow];
                            whichRow[nDelete++] = jRow;
                            int start = rowStart[jRow];
                            stored.addCut(rowLower[jRow], rowUpper[jRow],
                                          rowLength[jRow], column + start, elementByRow + start);
                        }
                    }
                    finalModel->deleteRows(nDelete, whichRow);
                }
                knapsackStart[numberKnapsack] = finalModel->getNumCols();
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
#endif	//COIN_HAS_LINK


