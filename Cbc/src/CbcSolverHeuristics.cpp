/* $Id$ */
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).


/*! \file CbcSolverHeuristics.cpp
    \brief Second level routines for the cbc stand-alone solver.
*/

#include "CbcConfig.h"
#include "CoinPragma.hpp"


#include "CoinTime.hpp"

#include "OsiClpSolverInterface.hpp"

#include "ClpPresolve.hpp"

#include "CbcOrClpParam.hpp"

#include "CbcModel.hpp"

#include "CbcHeuristicLocal.hpp"
#include "CbcHeuristicPivotAndFix.hpp"
//#include "CbcHeuristicPivotAndComplement.hpp"
#include "CbcHeuristicRandRound.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicRINS.hpp"

#include "CbcHeuristicDiveCoefficient.hpp"
#include "CbcHeuristicDiveFractional.hpp"
#include "CbcHeuristicDiveGuided.hpp"
#include "CbcHeuristicDiveVectorLength.hpp"
#include "CbcHeuristicDivePseudoCost.hpp"
#include "CbcHeuristicDiveLineSearch.hpp"

#include "CbcStrategy.hpp"
#include "OsiAuxInfo.hpp"

#include "ClpSimplexOther.hpp"

// Crunch down model
void
crunchIt(ClpSimplex * model)
{
#ifdef JJF_ZERO
    model->dual();
#else
    int numberColumns = model->numberColumns();
    int numberRows = model->numberRows();
    // Use dual region
    double * rhs = model->dualRowSolution();
    int * whichRow = new int[3*numberRows];
    int * whichColumn = new int[2*numberColumns];
    int nBound;
    ClpSimplex * small = static_cast<ClpSimplexOther *> (model)->crunch(rhs, whichRow, whichColumn,
                         nBound, false, false);
    if (small) {
        small->dual();
        if (small->problemStatus() == 0) {
            model->setProblemStatus(0);
            static_cast<ClpSimplexOther *> (model)->afterCrunch(*small, whichRow, whichColumn, nBound);
        } else if (small->problemStatus() != 3) {
            model->setProblemStatus(1);
        } else {
            if (small->problemStatus() == 3) {
                // may be problems
                small->computeObjectiveValue();
                model->setObjectiveValue(small->objectiveValue());
                model->setProblemStatus(3);
            } else {
                model->setProblemStatus(3);
            }
        }
        delete small;
    } else {
        model->setProblemStatus(1);
    }
    delete [] whichRow;
    delete [] whichColumn;
#endif
}
/*
  On input
  doAction - 0 just fix in original and return NULL
             1 return fixed non-presolved solver
             2 as one but use presolve Inside this
	     3 use presolve and fix ones with large cost
             ? do heuristics and set best solution
	     ? do BAB and just set best solution
	     10+ then use lastSolution and relax a few
             -2 cleanup afterwards if using 2
  On output - number fixed
*/
OsiClpSolverInterface *
fixVubs(CbcModel & model, int skipZero2,
        int & doAction,
        CoinMessageHandler * /*generalMessageHandler*/,
        const double * lastSolution, double dextra[6],
        int extra[5])
{
    if (doAction == 11 && !lastSolution)
        lastSolution = model.bestSolution();
    assert (((doAction >= 0 && doAction <= 3) && !lastSolution) || (doAction == 11 && lastSolution));
    double fractionIntFixed = dextra[3];
    double fractionFixed = dextra[4];
    double fixAbove = dextra[2];
    double fixAboveValue = (dextra[5] > 0.0) ? dextra[5] : 1.0;
#ifdef COIN_DETAIL
    double time1 = CoinCpuTime();
#endif
    int leaveIntFree = extra[1];
    OsiSolverInterface * originalSolver = model.solver();
    OsiClpSolverInterface * originalClpSolver = dynamic_cast< OsiClpSolverInterface*> (originalSolver);
    ClpSimplex * originalLpSolver = originalClpSolver->getModelPtr();
    int * originalColumns = NULL;
    OsiClpSolverInterface * clpSolver;
    ClpSimplex * lpSolver;
    ClpPresolve pinfo;
    assert(originalSolver->getObjSense() > 0);
    if (doAction == 2 || doAction == 3) {
        double * saveLB = NULL;
        double * saveUB = NULL;
        int numberColumns = originalLpSolver->numberColumns();
        if (fixAbove > 0.0) {
#ifdef COIN_DETAIL
            double time1 = CoinCpuTime();
#endif
            originalClpSolver->initialSolve();
            COIN_DETAIL_PRINT(printf("first solve took %g seconds\n", CoinCpuTime() - time1));
            double * columnLower = originalLpSolver->columnLower() ;
            double * columnUpper = originalLpSolver->columnUpper() ;
            const double * solution = originalLpSolver->primalColumnSolution();
            saveLB = CoinCopyOfArray(columnLower, numberColumns);
            saveUB = CoinCopyOfArray(columnUpper, numberColumns);
            const double * objective = originalLpSolver->getObjCoefficients() ;
            int iColumn;
            int nFix = 0;
            int nArt = 0;
            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                if (objective[iColumn] > fixAbove) {
                    if (solution[iColumn] < columnLower[iColumn] + 1.0e-8) {
                        columnUpper[iColumn] = columnLower[iColumn];
                        nFix++;
                    } else {
                        nArt++;
                    }
                } else if (objective[iColumn] < -fixAbove) {
                    if (solution[iColumn] > columnUpper[iColumn] - 1.0e-8) {
                        columnLower[iColumn] = columnUpper[iColumn];
                        nFix++;
                    } else {
                        nArt++;
                    }
                }
            }
            COIN_DETAIL_PRINT(printf("%d artificials fixed, %d left as in solution\n", nFix, nArt));
            lpSolver = pinfo.presolvedModel(*originalLpSolver, 1.0e-8, true, 10);
            if (!lpSolver || doAction == 2) {
                // take off fixing in original
                memcpy(columnLower, saveLB, numberColumns*sizeof(double));
                memcpy(columnUpper, saveUB, numberColumns*sizeof(double));
            }
            delete [] saveLB;
            delete [] saveUB;
            if (!lpSolver) {
                // try again
                pinfo.destroyPresolve();
                lpSolver = pinfo.presolvedModel(*originalLpSolver, 1.0e-8, true, 10);
                assert (lpSolver);
            }
        } else {
            lpSolver = pinfo.presolvedModel(*originalLpSolver, 1.0e-8, true, 10);
            assert (lpSolver);
        }
        clpSolver = new OsiClpSolverInterface(lpSolver, true);
        assert(lpSolver == clpSolver->getModelPtr());
        numberColumns = lpSolver->numberColumns();
        originalColumns = CoinCopyOfArray(pinfo.originalColumns(), numberColumns);
        doAction = 1;
    } else {
        OsiSolverInterface * solver = originalSolver->clone();
        clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
        lpSolver = clpSolver->getModelPtr();
    }
    // Tighten bounds
    lpSolver->tightenPrimalBounds(0.0, 11, true);
    int numberColumns = clpSolver->getNumCols() ;
    double * saveColumnLower = CoinCopyOfArray(lpSolver->columnLower(), numberColumns);
    double * saveColumnUpper = CoinCopyOfArray(lpSolver->columnUpper(), numberColumns);
    //char generalPrint[200];
    const double *objective = lpSolver->getObjCoefficients() ;
    double *columnLower = lpSolver->columnLower() ;
    double *columnUpper = lpSolver->columnUpper() ;
    int numberRows = clpSolver->getNumRows();
    int iRow, iColumn;

    // Row copy
    CoinPackedMatrix matrixByRow(*clpSolver->getMatrixByRow());
    const double * elementByRow = matrixByRow.getElements();
    const int * column = matrixByRow.getIndices();
    const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
    const int * rowLength = matrixByRow.getVectorLengths();

    // Column copy
    CoinPackedMatrix  matrixByCol(*clpSolver->getMatrixByCol());
    //const double * element = matrixByCol.getElements();
    const int * row = matrixByCol.getIndices();
    const CoinBigIndex * columnStart = matrixByCol.getVectorStarts();
    const int * columnLength = matrixByCol.getVectorLengths();

    const double * rowLower = clpSolver->getRowLower();
    const double * rowUpper = clpSolver->getRowUpper();

    // Get maximum size of VUB tree
    // otherColumn is one fixed to 0 if this one zero
    int nEl = matrixByCol.getNumElements();
    CoinBigIndex * fixColumn = new CoinBigIndex [numberColumns+1];
    int * otherColumn = new int [nEl];
    int * fix = new int[numberColumns];
    char * mark = new char [numberColumns];
    memset(mark, 0, numberColumns);
    int numberInteger = 0;
    int numberOther = 0;
    fixColumn[0] = 0;
    double large = lpSolver->largeValue(); // treat bounds > this as infinite
#ifndef NDEBUG
    double large2 = 1.0e10 * large;
#endif
    double tolerance = lpSolver->primalTolerance();
    int * check = new int[numberRows];
    for (iRow = 0; iRow < numberRows; iRow++) {
        check[iRow] = -2; // don't check
        if (rowLower[iRow] < -1.0e6 && rowUpper[iRow] > 1.0e6)
            continue;// unlikely
        // possible row
        int numberPositive = 0;
        int iPositive = -1;
        int numberNegative = 0;
        int iNegative = -1;
        CoinBigIndex rStart = rowStart[iRow];
        CoinBigIndex rEnd = rowStart[iRow] + rowLength[iRow];
        CoinBigIndex j;
        int kColumn;
        for (j = rStart; j < rEnd; ++j) {
            double value = elementByRow[j];
            kColumn = column[j];
            if (columnUpper[kColumn] > columnLower[kColumn]) {
                if (value > 0.0) {
                    numberPositive++;
                    iPositive = kColumn;
                } else {
                    numberNegative++;
                    iNegative = kColumn;
                }
            }
        }
        if (numberPositive == 1 && numberNegative == 1)
            check[iRow] = -1; // try both
        if (numberPositive == 1 && rowLower[iRow] > -1.0e20)
            check[iRow] = iPositive;
        else if (numberNegative == 1 && rowUpper[iRow] < 1.0e20)
            check[iRow] = iNegative;
    }
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        fix[iColumn] = -1;
        if (columnUpper[iColumn] > columnLower[iColumn] + 1.0e-8) {
            if (clpSolver->isInteger(iColumn))
                numberInteger++;
            if (columnLower[iColumn] == 0.0) {
                bool infeasible = false;
                fix[iColumn] = 0;
                // fake upper bound
                double saveUpper = columnUpper[iColumn];
                columnUpper[iColumn] = 0.0;
                for (CoinBigIndex i = columnStart[iColumn];
                        i < columnStart[iColumn] + columnLength[iColumn]; i++) {
                    iRow = row[i];
                    if (check[iRow] != -1 && check[iRow] != iColumn)
                        continue; // unlikely
                    // possible row
                    int infiniteUpper = 0;
                    int infiniteLower = 0;
                    double maximumUp = 0.0;
                    double maximumDown = 0.0;
                    double newBound;
                    CoinBigIndex rStart = rowStart[iRow];
                    CoinBigIndex rEnd = rowStart[iRow] + rowLength[iRow];
                    CoinBigIndex j;
                    int kColumn;
                    // Compute possible lower and upper ranges
                    for (j = rStart; j < rEnd; ++j) {
                        double value = elementByRow[j];
                        kColumn = column[j];
                        if (value > 0.0) {
                            if (columnUpper[kColumn] >= large) {
                                ++infiniteUpper;
                            } else {
                                maximumUp += columnUpper[kColumn] * value;
                            }
                            if (columnLower[kColumn] <= -large) {
                                ++infiniteLower;
                            } else {
                                maximumDown += columnLower[kColumn] * value;
                            }
                        } else if (value < 0.0) {
                            if (columnUpper[kColumn] >= large) {
                                ++infiniteLower;
                            } else {
                                maximumDown += columnUpper[kColumn] * value;
                            }
                            if (columnLower[kColumn] <= -large) {
                                ++infiniteUpper;
                            } else {
                                maximumUp += columnLower[kColumn] * value;
                            }
                        }
                    }
                    // Build in a margin of error
                    maximumUp += 1.0e-8 * fabs(maximumUp);
                    maximumDown -= 1.0e-8 * fabs(maximumDown);
                    double maxUp = maximumUp + infiniteUpper * 1.0e31;
                    double maxDown = maximumDown - infiniteLower * 1.0e31;
                    if (maxUp <= rowUpper[iRow] + tolerance &&
                            maxDown >= rowLower[iRow] - tolerance) {
                        //printf("Redundant row in vubs %d\n",iRow);
                    } else {
                        if (maxUp < rowLower[iRow] - 100.0*tolerance ||
                                maxDown > rowUpper[iRow] + 100.0*tolerance) {
                            infeasible = true;
                            break;
                        }
                        double lower = rowLower[iRow];
                        double upper = rowUpper[iRow];
                        for (j = rStart; j < rEnd; ++j) {
                            double value = elementByRow[j];
                            kColumn = column[j];
                            double nowLower = columnLower[kColumn];
                            double nowUpper = columnUpper[kColumn];
                            if (value > 0.0) {
                                // positive value
                                if (lower > -large) {
                                    if (!infiniteUpper) {
                                        assert(nowUpper < large2);
                                        newBound = nowUpper +
                                                   (lower - maximumUp) / value;
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
                                        // check infeasible (relaxed)
                                        if (nowUpper < newBound) {
                                            if (nowUpper - newBound <
                                                    -100.0*tolerance) {
                                                infeasible = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                                if (upper < large) {
                                    if (!infiniteLower) {
                                        assert(nowLower > - large2);
                                        newBound = nowLower +
                                                   (upper - maximumDown) / value;
                                        // relax if original was large
                                        if (fabs(maximumDown) > 1.0e8)
                                            newBound += 1.0e-12 * fabs(maximumDown);
                                    } else if (infiniteLower == 1 && nowLower < -large) {
                                        newBound =   (upper - maximumDown) / value;
                                        // relax if original was large
                                        if (fabs(maximumDown) > 1.0e8)
                                            newBound += 1.0e-12 * fabs(maximumDown);
                                    } else {
                                        newBound = COIN_DBL_MAX;
                                    }
                                    if (newBound < nowUpper - 1.0e-12 && newBound < large) {
                                        // Tighten the upper bound
                                        // check infeasible (relaxed)
                                        if (nowLower > newBound) {
                                            if (newBound - nowLower <
                                                    -100.0*tolerance) {
                                                infeasible = true;
                                                break;
                                            } else {
                                                newBound = nowLower;
                                            }
                                        }
                                        if (!newBound || (clpSolver->isInteger(kColumn) && newBound < 0.999)) {
                                            // fix to zero
                                            if (!mark[kColumn]) {
                                                otherColumn[numberOther++] = kColumn;
                                                mark[kColumn] = 1;
                                                if (check[iRow] == -1)
                                                    check[iRow] = iColumn;
                                                else
                                                    assert(check[iRow] == iColumn);
                                            }
                                        }
                                    }
                                }
                            } else {
                                // negative value
                                if (lower > -large) {
                                    if (!infiniteUpper) {
                                        assert(nowLower < large2);
                                        newBound = nowLower +
                                                   (lower - maximumUp) / value;
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
                                        // check infeasible (relaxed)
                                        if (nowLower > newBound) {
                                            if (newBound - nowLower <
                                                    -100.0*tolerance) {
                                                infeasible = true;
                                                break;
                                            } else {
                                                newBound = nowLower;
                                            }
                                        }
                                        if (!newBound || (clpSolver->isInteger(kColumn) && newBound < 0.999)) {
                                            // fix to zero
                                            if (!mark[kColumn]) {
                                                otherColumn[numberOther++] = kColumn;
                                                mark[kColumn] = 1;
                                                if (check[iRow] == -1)
                                                    check[iRow] = iColumn;
                                                else
                                                    assert(check[iRow] == iColumn);
                                            }
                                        }
                                    }
                                }
                                if (upper < large) {
                                    if (!infiniteLower) {
                                        assert(nowUpper < large2);
                                        newBound = nowUpper +
                                                   (upper - maximumDown) / value;
                                        // relax if original was large
                                        if (fabs(maximumDown) > 1.0e8)
                                            newBound -= 1.0e-12 * fabs(maximumDown);
                                    } else if (infiniteLower == 1 && nowUpper > large) {
                                        newBound =   (upper - maximumDown) / value;
                                        // relax if original was large
                                        if (fabs(maximumDown) > 1.0e8)
                                            newBound -= 1.0e-12 * fabs(maximumDown);
                                    } else {
                                        newBound = -COIN_DBL_MAX;
                                    }
                                    if (newBound > nowLower + 1.0e-12 && newBound > -large) {
                                        // Tighten the lower bound
                                        // check infeasible (relaxed)
                                        if (nowUpper < newBound) {
                                            if (nowUpper - newBound <
                                                    -100.0*tolerance) {
                                                infeasible = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                for (int i = fixColumn[iColumn]; i < numberOther; i++)
                    mark[otherColumn[i]] = 0;
                // reset bound unless infeasible
                if (!infeasible || !clpSolver->isInteger(iColumn))
                    columnUpper[iColumn] = saveUpper;
                else if (clpSolver->isInteger(iColumn))
                    columnLower[iColumn] = 1.0;
            }
        }
        fixColumn[iColumn+1] = numberOther;
    }
    delete [] check;
    delete [] mark;
    // Now do reverse way
    int * counts = new int [numberColumns];
    CoinZeroN(counts, numberColumns);
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        for (int i = fixColumn[iColumn]; i < fixColumn[iColumn+1]; i++)
            counts[otherColumn[i]]++;
    }
    numberOther = 0;
    CoinBigIndex * fixColumn2 = new CoinBigIndex [numberColumns+1];
    int * otherColumn2 = new int [fixColumn[numberColumns]];
    fixColumn2[0] = 0;
    for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
        numberOther += counts[iColumn];
        counts[iColumn] = 0;
        fixColumn2[iColumn+1] = numberOther;
    }
    // Create other way
    for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
        for (int i = fixColumn[iColumn]; i < fixColumn[iColumn+1]; i++) {
            int jColumn = otherColumn[i];
            CoinBigIndex put = fixColumn2[jColumn] + counts[jColumn];
            counts[jColumn]++;
            otherColumn2[put] = iColumn;
        }
    }
    // get top layer i.e. those which are not fixed by any other
    int kLayer = 0;
    while (true) {
        int numberLayered = 0;
        for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (fix[iColumn] == kLayer) {
                for (int i = fixColumn2[iColumn]; i < fixColumn2[iColumn+1]; i++) {
                    int jColumn = otherColumn2[i];
                    if (fix[jColumn] == kLayer) {
                        fix[iColumn] = kLayer + 100;
                    }
                }
            }
            if (fix[iColumn] == kLayer) {
                numberLayered++;
            }
        }
        if (numberLayered) {
            kLayer += 100;
        } else {
            break;
        }
    }
    for (int iPass = 0; iPass < 2; iPass++) {
        for (int jLayer = 0; jLayer < kLayer; jLayer++) {
            int check[] = { -1, 0, 1, 2, 3, 4, 5, 10, 50, 100, 500, 1000, 5000, 10000, COIN_INT_MAX};
            int nCheck = static_cast<int> (sizeof(check) / sizeof(int));
            int countsI[20];
            int countsC[20];
            assert (nCheck <= 20);
            memset(countsI, 0, nCheck*sizeof(int));
            memset(countsC, 0, nCheck*sizeof(int));
            check[nCheck-1] = numberColumns;
            int numberLayered = 0;
            int numberInteger = 0;
            for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
                if (fix[iColumn] == jLayer) {
                    numberLayered++;
                    int nFix = fixColumn[iColumn+1] - fixColumn[iColumn];
                    if (iPass) {
                        // just integers
                        nFix = 0;
                        for (int i = fixColumn[iColumn]; i < fixColumn[iColumn+1]; i++) {
                            int jColumn = otherColumn[i];
                            if (clpSolver->isInteger(jColumn))
                                nFix++;
                        }
                    }
                    int iFix;
                    for (iFix = 0; iFix < nCheck; iFix++) {
                        if (nFix <= check[iFix])
                            break;
                    }
                    assert (iFix < nCheck);
                    if (clpSolver->isInteger(iColumn)) {
                        numberInteger++;
                        countsI[iFix]++;
                    } else {
                        countsC[iFix]++;
                    }
                }
            }
#ifdef COIN_DETAIL
            if (numberLayered) {
	        printf("%d (%d integer) at priority %d\n", numberLayered, numberInteger, 1 + (jLayer / 100));
                char buffer[50];
                for (int i = 1; i < nCheck; i++) {
                    if (countsI[i] || countsC[i]) {
                        if (i == 1)
                            sprintf(buffer, " ==    zero            ");
                        else if (i < nCheck - 1)
                            sprintf(buffer, "> %6d and <= %6d ", check[i-1], check[i]);
                        else
                            sprintf(buffer, "> %6d                ", check[i-1]);
                        printf("%s %8d integers and %8d continuous\n", buffer, countsI[i], countsC[i]);
                    }
                }
            }
#endif
        }
    }
    delete [] counts;
    // Now do fixing
    {
        // switch off presolve and up weight
        ClpSolve solveOptions;
        //solveOptions.setPresolveType(ClpSolve::presolveOff,0);
        solveOptions.setSolveType(ClpSolve::usePrimalorSprint);
        //solveOptions.setSpecialOption(1,3,30); // sprint
        int numberColumns = lpSolver->numberColumns();
        int iColumn;
        bool allSlack = true;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (lpSolver->getColumnStatus(iColumn) == ClpSimplex::basic) {
                allSlack = false;
                break;
            }
        }
        if (allSlack)
            solveOptions.setSpecialOption(1, 2, 50); // idiot
        lpSolver->setInfeasibilityCost(1.0e11);
        lpSolver->defaultFactorizationFrequency();
        if (doAction != 11)
            lpSolver->initialSolve(solveOptions);
        double * columnLower = lpSolver->columnLower();
        double * columnUpper = lpSolver->columnUpper();
        double * fullSolution = lpSolver->primalColumnSolution();
        const double * dj = lpSolver->dualColumnSolution();
        int iPass = 0;
#define MAXPROB 2
        ClpSimplex models[MAXPROB];
        int pass[MAXPROB];
        int kPass = -1;
        int kLayer = 0;
        int skipZero = 0;
        if (skipZero2 == -1)
            skipZero2 = 40; //-1;
        /* 0 fixed to 0 by choice
           1 lb of 1 by choice
           2 fixed to 0 by another
           3 as 2 but this go
           -1 free
        */
        char * state = new char [numberColumns];
        for (iColumn = 0; iColumn < numberColumns; iColumn++)
            state[iColumn] = -1;
        while (true) {
            double largest = -0.1;
            double smallest = 1.1;
            int iLargest = -1;
            int iSmallest = -1;
            int atZero = 0;
            int atOne = 0;
            int toZero = 0;
            int toOne = 0;
            int numberFree = 0;
            int numberGreater = 0;
            columnLower = lpSolver->columnLower();
            columnUpper = lpSolver->columnUpper();
            fullSolution = lpSolver->primalColumnSolution();
            if (doAction == 11) {
                {
                    double * columnLower = lpSolver->columnLower();
                    double * columnUpper = lpSolver->columnUpper();
                    //	  lpSolver->dual();
                    memcpy(columnLower, saveColumnLower, numberColumns*sizeof(double));
                    memcpy(columnUpper, saveColumnUpper, numberColumns*sizeof(double));
                    //	  lpSolver->dual();
                    int iColumn;
                    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                        if (columnUpper[iColumn] > columnLower[iColumn] + 1.0e-8) {
                            if (clpSolver->isInteger(iColumn)) {
                                double value = lastSolution[iColumn];
                                int iValue = static_cast<int> (value + 0.5);
                                assert (fabs(value - static_cast<double> (iValue)) < 1.0e-3);
                                assert (iValue >= columnLower[iColumn] &&
                                        iValue <= columnUpper[iColumn]);
                                columnLower[iColumn] = iValue;
                                columnUpper[iColumn] = iValue;
                            }
                        }
                    }
                    lpSolver->initialSolve(solveOptions);
                    memcpy(columnLower, saveColumnLower, numberColumns*sizeof(double));
                    memcpy(columnUpper, saveColumnUpper, numberColumns*sizeof(double));
                }
                for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (columnUpper[iColumn] > columnLower[iColumn] + 1.0e-8) {
                        if (clpSolver->isInteger(iColumn)) {
                            double value = lastSolution[iColumn];
                            int iValue = static_cast<int> (value + 0.5);
                            assert (fabs(value - static_cast<double> (iValue)) < 1.0e-3);
                            assert (iValue >= columnLower[iColumn] &&
                                    iValue <= columnUpper[iColumn]);
                            if (!fix[iColumn]) {
                                if (iValue == 0) {
                                    state[iColumn] = 0;
                                    assert (!columnLower[iColumn]);
                                    columnUpper[iColumn] = 0.0;
                                } else if (iValue == 1) {
                                    state[iColumn] = 1;
                                    columnLower[iColumn] = 1.0;
                                } else {
                                    // leave fixed
                                    columnLower[iColumn] = iValue;
                                    columnUpper[iColumn] = iValue;
                                }
                            } else if (iValue == 0) {
                                state[iColumn] = 10;
                                columnUpper[iColumn] = 0.0;
                            } else {
                                // leave fixed
                                columnLower[iColumn] = iValue;
                                columnUpper[iColumn] = iValue;
                            }
                        }
                    }
                }
                int jLayer = 0;
                int nFixed = -1;
                int nTotalFixed = 0;
                while (nFixed) {
                    nFixed = 0;
                    for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
                        if (columnUpper[iColumn] == 0.0 && fix[iColumn] == jLayer) {
                            for (int i = fixColumn[iColumn]; i < fixColumn[iColumn+1]; i++) {
                                int jColumn = otherColumn[i];
                                if (columnUpper[jColumn]) {
                                    bool canFix = true;
                                    for (int k = fixColumn2[jColumn]; k < fixColumn2[jColumn+1]; k++) {
                                        int kColumn = otherColumn2[k];
                                        if (state[kColumn] == 1) {
                                            canFix = false;
                                            break;
                                        }
                                    }
                                    if (canFix) {
                                        columnUpper[jColumn] = 0.0;
                                        nFixed++;
                                    }
                                }
                            }
                        }
                    }
                    nTotalFixed += nFixed;
                    jLayer += 100;
                }
                COIN_DETAIL_PRINT(printf("This fixes %d variables in lower priorities\n", nTotalFixed));
                break;
            }
            for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
                if (!clpSolver->isInteger(iColumn) || fix[iColumn] > kLayer)
                    continue;
                // skip if fixes nothing
                if (fixColumn[iColumn+1] - fixColumn[iColumn] <= skipZero2)
                    continue;
                double value = fullSolution[iColumn];
                if (value > 1.00001) {
                    numberGreater++;
                    continue;
                }
                double lower = columnLower[iColumn];
                double upper = columnUpper[iColumn];
                if (lower == upper) {
                    if (lower)
                        atOne++;
                    else
                        atZero++;
                    continue;
                }
                if (value < 1.0e-7) {
                    toZero++;
                    columnUpper[iColumn] = 0.0;
                    state[iColumn] = 10;
                    continue;
                }
                if (value > 1.0 - 1.0e-7) {
                    toOne++;
                    columnLower[iColumn] = 1.0;
                    state[iColumn] = 1;
                    continue;
                }
                numberFree++;
                // skip if fixes nothing
                if (fixColumn[iColumn+1] - fixColumn[iColumn] <= skipZero)
                    continue;
                if (value < smallest) {
                    smallest = value;
                    iSmallest = iColumn;
                }
                if (value > largest) {
                    largest = value;
                    iLargest = iColumn;
                }
            }
            if (toZero || toOne)
	      COIN_DETAIL_PRINT(printf("%d at 0 fixed and %d at one fixed\n", toZero, toOne));
            COIN_DETAIL_PRINT(printf("%d variables free, %d fixed to 0, %d to 1 - smallest %g, largest %g\n",
				     numberFree, atZero, atOne, smallest, largest));
            if (numberGreater && !iPass)
	      COIN_DETAIL_PRINT(printf("%d variables have value > 1.0\n", numberGreater));
            //skipZero2=0; // leave 0 fixing
            int jLayer = 0;
            int nFixed = -1;
            int nTotalFixed = 0;
            while (nFixed) {
                nFixed = 0;
                for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (columnUpper[iColumn] == 0.0 && fix[iColumn] == jLayer) {
                        for (int i = fixColumn[iColumn]; i < fixColumn[iColumn+1]; i++) {
                            int jColumn = otherColumn[i];
                            if (columnUpper[jColumn]) {
                                bool canFix = true;
                                for (int k = fixColumn2[jColumn]; k < fixColumn2[jColumn+1]; k++) {
                                    int kColumn = otherColumn2[k];
                                    if (state[kColumn] == 1) {
                                        canFix = false;
                                        break;
                                    }
                                }
                                if (canFix) {
                                    columnUpper[jColumn] = 0.0;
                                    nFixed++;
                                }
                            }
                        }
                    }
                }
                nTotalFixed += nFixed;
                jLayer += 100;
            }
            COIN_DETAIL_PRINT(printf("This fixes %d variables in lower priorities\n", nTotalFixed));
            if (iLargest < 0 || numberFree <= leaveIntFree)
                break;
            double movement;
            int way;
            if (smallest <= 1.0 - largest && smallest < 0.2 && largest < fixAboveValue) {
                columnUpper[iSmallest] = 0.0;
                state[iSmallest] = 0;
                movement = smallest;
                way = -1;
            } else {
                columnLower[iLargest] = 1.0;
                state[iLargest] = 1;
                movement = 1.0 - largest;
                way = 1;
            }
            double saveObj = lpSolver->objectiveValue();
            iPass++;
            kPass = iPass % MAXPROB;
            models[kPass] = *lpSolver;
            if (way == -1) {
                // fix others
                for (int i = fixColumn[iSmallest]; i < fixColumn[iSmallest+1]; i++) {
                    int jColumn = otherColumn[i];
                    if (state[jColumn] == -1) {
                        columnUpper[jColumn] = 0.0;
                        state[jColumn] = 3;
                    }
                }
            }
            pass[kPass] = iPass;
            double maxCostUp = COIN_DBL_MAX;
            objective = lpSolver->getObjCoefficients() ;
            if (way == -1)
                maxCostUp = (1.0 - movement) * objective[iSmallest];
            lpSolver->setDualObjectiveLimit(saveObj + maxCostUp);
            crunchIt(lpSolver);
            double moveObj = lpSolver->objectiveValue() - saveObj;
            COIN_DETAIL_PRINT(printf("movement %s was %g costing %g\n",
				     (way == -1) ? "down" : "up", movement, moveObj));
            if (way == -1 && (moveObj >= maxCostUp || lpSolver->status())) {
                // go up
                columnLower = models[kPass].columnLower();
                columnUpper = models[kPass].columnUpper();
                columnLower[iSmallest] = 1.0;
                columnUpper[iSmallest] = saveColumnUpper[iSmallest];
                *lpSolver = models[kPass];
                columnLower = lpSolver->columnLower();
                columnUpper = lpSolver->columnUpper();
                fullSolution = lpSolver->primalColumnSolution();
                dj = lpSolver->dualColumnSolution();
                columnLower[iSmallest] = 1.0;
                columnUpper[iSmallest] = saveColumnUpper[iSmallest];
                state[iSmallest] = 1;
                // unfix others
                for (int i = fixColumn[iSmallest]; i < fixColumn[iSmallest+1]; i++) {
                    int jColumn = otherColumn[i];
                    if (state[jColumn] == 3) {
                        columnUpper[jColumn] = saveColumnUpper[jColumn];
                        state[jColumn] = -1;
                    }
                }
                crunchIt(lpSolver);
            }
            models[kPass] = *lpSolver;
        }
        lpSolver->dual();
        COIN_DETAIL_PRINT(printf("Fixing took %g seconds\n", CoinCpuTime() - time1));
        columnLower = lpSolver->columnLower();
        columnUpper = lpSolver->columnUpper();
        fullSolution = lpSolver->primalColumnSolution();
        dj = lpSolver->dualColumnSolution();
        int * sort = new int[numberColumns];
        double * dsort = new double[numberColumns];
        int chunk = 20;
        int iRelax = 0;
        //double fractionFixed=6.0/8.0;
        // relax while lots fixed
        while (true) {
            if (skipZero2 > 10 && doAction < 10)
                break;
            iRelax++;
            int n = 0;
            double sum0 = 0.0;
            double sum00 = 0.0;
            double sum1 = 0.0;
            for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
                if (!clpSolver->isInteger(iColumn) || fix[iColumn] > kLayer)
                    continue;
                // skip if fixes nothing
                if (fixColumn[iColumn+1] - fixColumn[iColumn] == 0 && doAction < 10)
                    continue;
                double djValue = dj[iColumn];
                if (state[iColumn] == 1) {
                    assert (columnLower[iColumn]);
                    assert (fullSolution[iColumn] > 0.1);
                    if (djValue > 0.0) {
                        //printf("YY dj of %d at %g is %g\n",iColumn,value,djValue);
                        sum1 += djValue;
                        sort[n] = iColumn;
                        dsort[n++] = -djValue;
                    } else {
                        //printf("dj of %d at %g is %g\n",iColumn,value,djValue);
                    }
                } else if (state[iColumn] == 0 || state[iColumn] == 10) {
                    assert (fullSolution[iColumn] < 0.1);
                    assert (!columnUpper[iColumn]);
                    double otherValue = 0.0;
                    int nn = 0;
                    for (int i = fixColumn[iColumn]; i < fixColumn[iColumn+1]; i++) {
                        int jColumn = otherColumn[i];
                        if (columnUpper[jColumn] == 0.0) {
                            if (dj[jColumn] < -1.0e-5) {
                                nn++;
                                otherValue += dj[jColumn]; // really need to look at rest
                            }
                        }
                    }
                    if (djValue < -1.0e-2 || otherValue < -1.0e-2) {
                        //printf("XX dj of %d at %g is %g - %d out of %d contribute %g\n",iColumn,value,djValue,
                        // nn,fixColumn[iColumn+1]-fixColumn[iColumn],otherValue);
                        if (djValue < 1.0e-8) {
                            sum0 -= djValue;
                            sum00 -= otherValue;
                            sort[n] = iColumn;
                            if (djValue < -1.0e-2)
                                dsort[n++] = djValue + otherValue;
                            else
                                dsort[n++] = djValue + 0.001 * otherValue;
                        }
                    } else {
                        //printf("dj of %d at %g is %g - no contribution from %d\n",iColumn,value,djValue,
                        //   fixColumn[iColumn+1]-fixColumn[iColumn]);
                    }
                }
            }
            CoinSort_2(dsort, dsort + n, sort);
            double * originalColumnLower = saveColumnLower;
            double * originalColumnUpper = saveColumnUpper;
            double * lo = CoinCopyOfArray(columnLower, numberColumns);
            double * up = CoinCopyOfArray(columnUpper, numberColumns);
            for (int k = 0; k < CoinMin(chunk, n); k++) {
                iColumn = sort[k];
                state[iColumn] = -2;
            }
            memcpy(columnLower, originalColumnLower, numberColumns*sizeof(double));
            memcpy(columnUpper, originalColumnUpper, numberColumns*sizeof(double));
            int nFixed = 0;
            int nFixed0 = 0;
            int nFixed1 = 0;
            for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
                if (state[iColumn] == 0 || state[iColumn] == 10) {
                    columnUpper[iColumn] = 0.0;
                    assert (lo[iColumn] == 0.0);
                    nFixed++;
                    nFixed0++;
                    for (int i = fixColumn[iColumn]; i < fixColumn[iColumn+1]; i++) {
                        int jColumn = otherColumn[i];
                        if (columnUpper[jColumn]) {
                            bool canFix = true;
                            for (int k = fixColumn2[jColumn]; k < fixColumn2[jColumn+1]; k++) {
                                int kColumn = otherColumn2[k];
                                if (state[kColumn] == 1 || state[kColumn] == -2) {
                                    canFix = false;
                                    break;
                                }
                            }
                            if (canFix) {
                                columnUpper[jColumn] = 0.0;
                                assert (lo[jColumn] == 0.0);
                                nFixed++;
                            }
                        }
                    }
                } else if (state[iColumn] == 1) {
                    columnLower[iColumn] = 1.0;
                    nFixed1++;
                }
            }
            COIN_DETAIL_PRINT(printf("%d fixed %d orig 0 %d 1\n", nFixed, nFixed0, nFixed1));
            int jLayer = 0;
            nFixed = -1;
            int nTotalFixed = 0;
            while (nFixed) {
                nFixed = 0;
                for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (columnUpper[iColumn] == 0.0 && fix[iColumn] == jLayer) {
                        for (int i = fixColumn[iColumn]; i < fixColumn[iColumn+1]; i++) {
                            int jColumn = otherColumn[i];
                            if (columnUpper[jColumn]) {
                                bool canFix = true;
                                for (int k = fixColumn2[jColumn]; k < fixColumn2[jColumn+1]; k++) {
                                    int kColumn = otherColumn2[k];
                                    if (state[kColumn] == 1 || state[kColumn] == -2) {
                                        canFix = false;
                                        break;
                                    }
                                }
                                if (canFix) {
                                    columnUpper[jColumn] = 0.0;
                                    assert (lo[jColumn] == 0.0);
                                    nFixed++;
                                }
                            }
                        }
                    }
                }
                nTotalFixed += nFixed;
                jLayer += 100;
            }
            nFixed = 0;
            int nFixedI = 0;
            for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
                if (columnLower[iColumn] == columnUpper[iColumn]) {
                    if (clpSolver->isInteger(iColumn))
                        nFixedI++;
                    nFixed++;
                }
            }
            COIN_DETAIL_PRINT(printf("This fixes %d variables in lower priorities - total %d (%d integer) - all target %d, int target %d\n",
				     nTotalFixed, nFixed, nFixedI, static_cast<int>(fractionFixed*numberColumns), static_cast<int> (fractionIntFixed*numberInteger)));
            int nBad = 0;
            int nRelax = 0;
            for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
                if (lo[iColumn] < columnLower[iColumn] ||
                        up[iColumn] > columnUpper[iColumn]) {
                    COIN_DETAIL_PRINT(printf("bad %d old %g %g, new %g %g\n", iColumn, lo[iColumn], up[iColumn],
					     columnLower[iColumn], columnUpper[iColumn]));
                    nBad++;
                }
                if (lo[iColumn] > columnLower[iColumn] ||
                        up[iColumn] < columnUpper[iColumn]) {
                    nRelax++;
                }
            }
            COIN_DETAIL_PRINT(printf("%d relaxed\n", nRelax));
            if (iRelax > 20 && nRelax == chunk)
                nRelax = 0;
            if (iRelax > 50)
                nRelax = 0;
            assert (!nBad);
            delete [] lo;
            delete [] up;
            lpSolver->primal(1);
            if (nFixed < fractionFixed*numberColumns || nFixedI < fractionIntFixed*numberInteger || !nRelax)
                break;
        }
        delete [] state;
        delete [] sort;
        delete [] dsort;
    }
    delete [] fix;
    delete [] fixColumn;
    delete [] otherColumn;
    delete [] otherColumn2;
    delete [] fixColumn2;
    // See if was presolved
    if (originalColumns) {
        columnLower = lpSolver->columnLower();
        columnUpper = lpSolver->columnUpper();
        for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
            saveColumnLower[iColumn] = columnLower[iColumn];
            saveColumnUpper[iColumn] = columnUpper[iColumn];
        }
        pinfo.postsolve(true);
        columnLower = originalLpSolver->columnLower();
        columnUpper = originalLpSolver->columnUpper();
        double * newColumnLower = lpSolver->columnLower();
        double * newColumnUpper = lpSolver->columnUpper();
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            int jColumn = originalColumns[iColumn];
            columnLower[jColumn] = CoinMax(columnLower[jColumn], newColumnLower[iColumn]);
            columnUpper[jColumn] = CoinMin(columnUpper[jColumn], newColumnUpper[iColumn]);
        }
        numberColumns = originalLpSolver->numberColumns();
        delete [] originalColumns;
    }
    delete [] saveColumnLower;
    delete [] saveColumnUpper;
    if (!originalColumns) {
        // Basis
        memcpy(originalLpSolver->statusArray(), lpSolver->statusArray(), numberRows + numberColumns);
        memcpy(originalLpSolver->primalColumnSolution(), lpSolver->primalColumnSolution(), numberColumns*sizeof(double));
        memcpy(originalLpSolver->primalRowSolution(), lpSolver->primalRowSolution(), numberRows*sizeof(double));
        // Fix in solver
        columnLower = lpSolver->columnLower();
        columnUpper = lpSolver->columnUpper();
    }
    double * originalColumnLower = originalLpSolver->columnLower();
    double * originalColumnUpper = originalLpSolver->columnUpper();
    // number fixed
    doAction = 0;
    for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
        originalColumnLower[iColumn] = columnLower[iColumn];
        originalColumnUpper[iColumn] = columnUpper[iColumn];
        if (columnLower[iColumn] == columnUpper[iColumn])
            doAction++;
    }
    COIN_DETAIL_PRINT(printf("%d fixed by vub preprocessing\n", doAction));
    if (originalColumns) {
        originalLpSolver->initialSolve();
    }
    delete clpSolver;
    return NULL;
}

int doHeuristics(CbcModel * model, int type, CbcOrClpParam* parameters_,
		 int numberParameters_,int noPrinting_,int initialPumpTune)
{
#ifdef JJF_ZERO //NEW_STYLE_SOLVER==0
    CbcOrClpParam * parameters_ = parameters;
    int numberParameters_ = numberParameters;
    bool noPrinting_ = noPrinting_;
#endif
    char generalPrint[10000];
    CoinMessages generalMessages = model->messages();
    CoinMessageHandler * generalMessageHandler = model->messageHandler();
    //generalMessageHandler->setPrefix(false);
    bool anyToDo = false;
    int logLevel = parameters_[whichParam(CLP_PARAM_INT_LOGLEVEL, numberParameters_, parameters_)].intValue();
    int useFpump = parameters_[whichParam(CBC_PARAM_STR_FPUMP, numberParameters_, parameters_)].currentOptionAsInteger();
    int useRounding = parameters_[whichParam(CBC_PARAM_STR_ROUNDING, numberParameters_, parameters_)].currentOptionAsInteger();
    int useGreedy = parameters_[whichParam(CBC_PARAM_STR_GREEDY, numberParameters_, parameters_)].currentOptionAsInteger();
    int useCombine = parameters_[whichParam(CBC_PARAM_STR_COMBINE, numberParameters_, parameters_)].currentOptionAsInteger();
    int useCrossover = parameters_[whichParam(CBC_PARAM_STR_CROSSOVER2, numberParameters_, parameters_)].currentOptionAsInteger();
    //int usePivotC = parameters_[whichParam(CBC_PARAM_STR_PIVOTANDCOMPLEMENT, numberParameters_, parameters_)].currentOptionAsInteger();
    int usePivotF = parameters_[whichParam(CBC_PARAM_STR_PIVOTANDFIX, numberParameters_, parameters_)].currentOptionAsInteger();
    int useRand = parameters_[whichParam(CBC_PARAM_STR_RANDROUND, numberParameters_, parameters_)].currentOptionAsInteger();
    int useRINS = parameters_[whichParam(CBC_PARAM_STR_RINS, numberParameters_, parameters_)].currentOptionAsInteger();
    int useRENS = parameters_[whichParam(CBC_PARAM_STR_RENS, numberParameters_, parameters_)].currentOptionAsInteger();
    int useDINS = parameters_[whichParam(CBC_PARAM_STR_DINS, numberParameters_, parameters_)].currentOptionAsInteger();
    int useDIVING2 = parameters_[whichParam(CBC_PARAM_STR_DIVINGS, numberParameters_, parameters_)].currentOptionAsInteger();
    int useNaive = parameters_[whichParam(CBC_PARAM_STR_NAIVE, numberParameters_, parameters_)].currentOptionAsInteger();
    int kType = (type < 10) ? type : 1;
    assert (kType == 1 || kType == 2);
    // FPump done first as it only works if no solution
    if (useFpump >= kType && useFpump <= kType + 1) {
        anyToDo = true;
        CbcHeuristicFPump heuristic4(*model);
        double dextra3 = parameters_[whichParam(CBC_PARAM_DBL_SMALLBAB, numberParameters_, parameters_)].doubleValue();
        heuristic4.setFractionSmall(dextra3);
        double dextra1 = parameters_[whichParam(CBC_PARAM_DBL_ARTIFICIALCOST, numberParameters_, parameters_)].doubleValue();
        if (dextra1)
            heuristic4.setArtificialCost(dextra1);
        heuristic4.setMaximumPasses(parameters_[whichParam(CBC_PARAM_INT_FPUMPITS, numberParameters_, parameters_)].intValue());
        if (parameters_[whichParam(CBC_PARAM_INT_FPUMPITS, numberParameters_, parameters_)].intValue() == 21)
            heuristic4.setIterationRatio(1.0);
        int pumpTune = parameters_[whichParam(CBC_PARAM_INT_FPUMPTUNE, numberParameters_, parameters_)].intValue();
        int pumpTune2 = parameters_[whichParam(CBC_PARAM_INT_FPUMPTUNE2, numberParameters_, parameters_)].intValue();
        if (pumpTune > 0) {
            bool printStuff = (pumpTune != initialPumpTune || logLevel > 1 || pumpTune2 > 0)
                              && !noPrinting_;
            if (printStuff) {
                generalMessageHandler->message(CBC_GENERAL, generalMessages)
                << "Options for feasibility pump - "
                << CoinMessageEol;
            }
            /*
            >=10000000 for using obj
            >=1000000 use as accumulate switch
            >=1000 use index+1 as number of large loops
            >=100 use dextra1 as cutoff
            %100 == 10,20 etc for experimentation
            1 == fix ints at bounds, 2 fix all integral ints, 3 and continuous at bounds
            4 and static continuous, 5 as 3 but no internal integers
            6 as 3 but all slack basis!
            */
            double value = model->solver()->getObjSense() * model->solver()->getObjValue();
            int w = pumpTune / 10;
            int i = w % 10;
            w /= 10;
            int c = w % 10;
            w /= 10;
            int r = w;
            int accumulate = r / 1000;
            r -= 1000 * accumulate;
            if (accumulate >= 10) {
                int which = accumulate / 10;
                accumulate -= 10 * which;
                which--;
                // weights and factors
                double weight[] = {0.01, 0.01, 0.1, 0.1, 0.5, 0.5, 1.0, 1.0, 5.0, 5.0};
                double factor[] = {0.1, 0.5, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5};
                heuristic4.setInitialWeight(weight[which]);
                heuristic4.setWeightFactor(factor[which]);
                if (printStuff) {
                    sprintf(generalPrint, "Initial weight for objective %g, decay factor %g",
                            weight[which], factor[which]);
                    generalMessageHandler->message(CBC_GENERAL, generalMessages)
                    << generalPrint
                    << CoinMessageEol;
                }

            }
            // fake cutoff
            if (c) {
                double cutoff;
                model->solver()->getDblParam(OsiDualObjectiveLimit, cutoff);
                cutoff = CoinMin(cutoff, value + 0.05 * fabs(value) * c);
                double fakeCutoff = parameters_[whichParam(CBC_PARAM_DBL_FAKECUTOFF, numberParameters_, parameters_)].doubleValue();
                if (fakeCutoff)
                    cutoff = fakeCutoff;
                heuristic4.setFakeCutoff(cutoff);
                if (printStuff) {
                    sprintf(generalPrint, "Fake cutoff of %g", cutoff);
                    generalMessageHandler->message(CBC_GENERAL, generalMessages)
                    << generalPrint
                    << CoinMessageEol;
                }
            }
            int offRandomEtc = 0;
            if (pumpTune2) {
                if ((pumpTune2 / 1000) != 0) {
                    offRandomEtc = 1000000 * (pumpTune2 / 1000);
                    if (printStuff) {
                        generalMessageHandler->message(CBC_GENERAL, generalMessages)
                        << "Feasibility pump may run twice"
                        << CoinMessageEol;
                    }
                    pumpTune2 = pumpTune2 % 1000;
                }
                if ((pumpTune2 / 100) != 0) {
                    offRandomEtc += 100 * (pumpTune2 / 100);
                    if (printStuff) {
                        generalMessageHandler->message(CBC_GENERAL, generalMessages)
                        << "Not using randomized objective"
                        << CoinMessageEol;
                    }
                }
                int maxAllowed = pumpTune2 % 100;
                if (maxAllowed) {
                    offRandomEtc += 1000 * maxAllowed;
                    if (printStuff) {
                        sprintf(generalPrint, "Fixing if same for %d passes",
                                maxAllowed);
                        generalMessageHandler->message(CBC_GENERAL, generalMessages)
                        << generalPrint
                        << CoinMessageEol;
                    }
                }
            }
            if (accumulate) {
                heuristic4.setAccumulate(accumulate);
                if (printStuff) {
                    if (accumulate) {
                        sprintf(generalPrint, "Accumulate of %d", accumulate);
                        generalMessageHandler->message(CBC_GENERAL, generalMessages)
                        << generalPrint
                        << CoinMessageEol;
                    }
                }
            }
            if (r) {
                // also set increment
                //double increment = (0.01*i+0.005)*(fabs(value)+1.0e-12);
                double increment = 0.0;
                double fakeIncrement = parameters_[whichParam(CBC_PARAM_DBL_FAKEINCREMENT, numberParameters_, parameters_)].doubleValue();
                if (fakeIncrement)
                    increment = fakeIncrement;
                heuristic4.setAbsoluteIncrement(increment);
                heuristic4.setMaximumRetries(r + 1);
                if (printStuff) {
                    if (increment) {
                        sprintf(generalPrint, "Increment of %g", increment);
                        generalMessageHandler->message(CBC_GENERAL, generalMessages)
                        << generalPrint
                        << CoinMessageEol;
                    }
                    sprintf(generalPrint, "%d retries", r + 1);
                    generalMessageHandler->message(CBC_GENERAL, generalMessages)
                    << generalPrint
                    << CoinMessageEol;
                }
            }
            if (i + offRandomEtc) {
                heuristic4.setFeasibilityPumpOptions(i*10 + offRandomEtc);
                if (printStuff) {
                    sprintf(generalPrint, "Feasibility pump options of %d",
                            i*10 + offRandomEtc);
                    generalMessageHandler->message(CBC_GENERAL, generalMessages)
                    << generalPrint
                    << CoinMessageEol;
                }
            }
            pumpTune = pumpTune % 100;
            if (pumpTune == 6)
                pumpTune = 13;
            heuristic4.setWhen((pumpTune % 10) + 10);
            if (printStuff) {
                sprintf(generalPrint, "Tuning (fixing) %d", pumpTune % 10);
                generalMessageHandler->message(CBC_GENERAL, generalMessages)
                << generalPrint
                << CoinMessageEol;
            }
        }
        heuristic4.setHeuristicName("feasibility pump");
        //#define ROLF
#ifdef ROLF
        CbcHeuristicFPump pump(*model);
        pump.setMaximumTime(60);
        pump.setMaximumPasses(100);
        pump.setMaximumRetries(1);
        pump.setFixOnReducedCosts(0);
        pump.setHeuristicName("Feasibility pump");
        pump.setFractionSmall(1.0);
        pump.setWhen(13);
        model->addHeuristic(&pump);
#else
        model->addHeuristic(&heuristic4);
#endif
    }
    if (useRounding >= type && useRounding >= kType && useRounding <= kType + 1) {
        CbcRounding heuristic1(*model);
        heuristic1.setHeuristicName("rounding");
        model->addHeuristic(&heuristic1) ;
        anyToDo = true;
    }
    if (useGreedy >= type && useGreedy >= kType && useGreedy <= kType + 1) {
        CbcHeuristicGreedyCover heuristic3(*model);
        heuristic3.setHeuristicName("greedy cover");
        CbcHeuristicGreedyEquality heuristic3a(*model);
        heuristic3a.setHeuristicName("greedy equality");
        model->addHeuristic(&heuristic3);
        model->addHeuristic(&heuristic3a);
        anyToDo = true;
    }
    if ((useRENS==7 && kType==1) || (useRENS==8 && kType==2)) {
        useRENS=1+2*(useRENS-7);
        CbcHeuristicRENS heuristic6a(*model);
        heuristic6a.setHeuristicName("RENSdj");
        heuristic6a.setFractionSmall(0.6/*3.4*/);
        heuristic6a.setFeasibilityPumpOptions(3);
        heuristic6a.setNumberNodes(10);
	heuristic6a.setWhereFrom(4*256+4*1);
	heuristic6a.setWhen(2);
	heuristic6a.setRensType(1+16);
        model->addHeuristic(&heuristic6a) ;
        heuristic6a.setHeuristicName("RENSub");
        heuristic6a.setFractionSmall(0.4);
        heuristic6a.setFeasibilityPumpOptions(1008003);
        heuristic6a.setNumberNodes(50);
	heuristic6a.setWhereFrom(4*256+4*1);
	heuristic6a.setWhen(2);
	heuristic6a.setRensType(2+16);
        model->addHeuristic(&heuristic6a) ;
    }
    if (useRENS >= kType && useRENS <= kType + 1) {
#ifndef JJF_ONE
        CbcHeuristicRENS heuristic6(*model);
        heuristic6.setHeuristicName("RENS");
        heuristic6.setFractionSmall(0.4);
        heuristic6.setFeasibilityPumpOptions(1008003);
        int nodes [] = { -2, 50, 50, 50, 200, 1000, 10000};
        heuristic6.setNumberNodes(nodes[useRENS]);
#else
        CbcHeuristicVND heuristic6(*model);
        heuristic6.setHeuristicName("VND");
        heuristic5.setFractionSmall(0.5);
        heuristic5.setDecayFactor(5.0);
#endif
        model->addHeuristic(&heuristic6) ;
        anyToDo = true;
    }
    if (useNaive >= kType && useNaive <= kType + 1) {
        CbcHeuristicNaive heuristic5b(*model);
        heuristic5b.setHeuristicName("Naive");
        heuristic5b.setFractionSmall(0.4);
        heuristic5b.setNumberNodes(50);
        model->addHeuristic(&heuristic5b) ;
        anyToDo = true;
    }
    int useDIVING = 0;
    {
        int useD;
        useD = parameters_[whichParam(CBC_PARAM_STR_DIVINGV, numberParameters_, parameters_)].currentOptionAsInteger();
        useDIVING |= 1 * ((useD >= kType) ? 1 : 0);
        useD = parameters_[whichParam(CBC_PARAM_STR_DIVINGG, numberParameters_, parameters_)].currentOptionAsInteger();
        useDIVING |= 2 * ((useD >= kType) ? 1 : 0);
        useD = parameters_[whichParam(CBC_PARAM_STR_DIVINGF, numberParameters_, parameters_)].currentOptionAsInteger();
        useDIVING |= 4 * ((useD >= kType) ? 1 : 0);
        useD = parameters_[whichParam(CBC_PARAM_STR_DIVINGC, numberParameters_, parameters_)].currentOptionAsInteger();
        useDIVING |= 8 * ((useD >= kType) ? 1 : 0);
        useD = parameters_[whichParam(CBC_PARAM_STR_DIVINGL, numberParameters_, parameters_)].currentOptionAsInteger();
        useDIVING |= 16 * ((useD >= kType) ? 1 : 0);
        useD = parameters_[whichParam(CBC_PARAM_STR_DIVINGP, numberParameters_, parameters_)].currentOptionAsInteger();
        useDIVING |= 32 * ((useD >= kType) ? 1 : 0);
    }
    if (useDIVING2 >= kType && useDIVING2 <= kType + 1) {
        int diveOptions = parameters_[whichParam(CBC_PARAM_INT_DIVEOPT, numberParameters_, parameters_)].intValue();
        if (diveOptions < 0 || diveOptions > 10)
            diveOptions = 2;
        CbcHeuristicJustOne heuristicJustOne(*model);
        heuristicJustOne.setHeuristicName("DiveAny");
        heuristicJustOne.setWhen(diveOptions);
        // add in others
        CbcHeuristicDiveCoefficient heuristicDC(*model);
        heuristicDC.setHeuristicName("DiveCoefficient");
        heuristicJustOne.addHeuristic(&heuristicDC, 1.0) ;
        CbcHeuristicDiveFractional heuristicDF(*model);
        heuristicDF.setHeuristicName("DiveFractional");
        heuristicJustOne.addHeuristic(&heuristicDF, 1.0) ;
        CbcHeuristicDiveGuided heuristicDG(*model);
        heuristicDG.setHeuristicName("DiveGuided");
        heuristicJustOne.addHeuristic(&heuristicDG, 1.0) ;
        CbcHeuristicDiveLineSearch heuristicDL(*model);
        heuristicDL.setHeuristicName("DiveLineSearch");
        heuristicJustOne.addHeuristic(&heuristicDL, 1.0) ;
        CbcHeuristicDivePseudoCost heuristicDP(*model);
        heuristicDP.setHeuristicName("DivePseudoCost");
        heuristicJustOne.addHeuristic(&heuristicDP, 1.0) ;
        CbcHeuristicDiveVectorLength heuristicDV(*model);
        heuristicDV.setHeuristicName("DiveVectorLength");
        heuristicJustOne.addHeuristic(&heuristicDV, 1.0) ;
        // Now normalize probabilities
        heuristicJustOne.normalizeProbabilities();
        model->addHeuristic(&heuristicJustOne) ;
    }

    if (useDIVING > 0) {
        int diveOptions2 = parameters_[whichParam(CBC_PARAM_INT_DIVEOPT, numberParameters_, parameters_)].intValue();
        int diveOptions;
        if (diveOptions2 > 99) {
            // switch on various active set stuff
            diveOptions = diveOptions2 % 100;
            diveOptions2 -= diveOptions;
        } else {
            diveOptions = diveOptions2;
            diveOptions2 = 0;
        }
        if (diveOptions < 0 || diveOptions > 9)
            diveOptions = 2;
        if ((useDIVING&1) != 0) {
            CbcHeuristicDiveVectorLength heuristicDV(*model);
            heuristicDV.setHeuristicName("DiveVectorLength");
            heuristicDV.setWhen(diveOptions);
            model->addHeuristic(&heuristicDV) ;
        }
        if ((useDIVING&2) != 0) {
            CbcHeuristicDiveGuided heuristicDG(*model);
            heuristicDG.setHeuristicName("DiveGuided");
            heuristicDG.setWhen(diveOptions);
            model->addHeuristic(&heuristicDG) ;
        }
        if ((useDIVING&4) != 0) {
            CbcHeuristicDiveFractional heuristicDF(*model);
            heuristicDF.setHeuristicName("DiveFractional");
            heuristicDF.setWhen(diveOptions);
            model->addHeuristic(&heuristicDF) ;
        }
        if ((useDIVING&8) != 0) {
            CbcHeuristicDiveCoefficient heuristicDC(*model);
            heuristicDC.setHeuristicName("DiveCoefficient");
            heuristicDC.setWhen(diveOptions);
            model->addHeuristic(&heuristicDC) ;
        }
        if ((useDIVING&16) != 0) {
            CbcHeuristicDiveLineSearch heuristicDL(*model);
            heuristicDL.setHeuristicName("DiveLineSearch");
            heuristicDL.setWhen(diveOptions);
            model->addHeuristic(&heuristicDL) ;
        }
        if ((useDIVING&32) != 0) {
            CbcHeuristicDivePseudoCost heuristicDP(*model);
            heuristicDP.setHeuristicName("DivePseudoCost");
            heuristicDP.setWhen(diveOptions + diveOptions2);
            model->addHeuristic(&heuristicDP) ;
        }
        anyToDo = true;
    }
#ifdef JJF_ZERO
    if (usePivotC >= type && usePivotC <= kType + 1) {
        CbcHeuristicPivotAndComplement heuristic(*model);
        heuristic.setHeuristicName("pivot and complement");
        heuristic.setFractionSmall(10.0); // normally 0.5
        model->addHeuristic(&heuristic);
        anyToDo = true;
    }
#endif
    if (usePivotF >= type && usePivotF <= kType + 1) {
        CbcHeuristicPivotAndFix heuristic(*model);
        heuristic.setHeuristicName("pivot and fix");
        heuristic.setFractionSmall(10.0); // normally 0.5
        model->addHeuristic(&heuristic);
        anyToDo = true;
    }
    if (useRand >= type && useRand <= kType + 1) {
        CbcHeuristicRandRound heuristic(*model);
        heuristic.setHeuristicName("randomized rounding");
        heuristic.setFractionSmall(10.0); // normally 0.5
        model->addHeuristic(&heuristic);
        anyToDo = true;
    }
    if (useDINS >= kType && useDINS <= kType + 1) {
        CbcHeuristicDINS heuristic5a(*model);
        heuristic5a.setHeuristicName("DINS");
        heuristic5a.setFractionSmall(0.6);
        if (useDINS < 4)
            heuristic5a.setDecayFactor(5.0);
        else
            heuristic5a.setDecayFactor(1.5);
        heuristic5a.setNumberNodes(1000);
        model->addHeuristic(&heuristic5a) ;
        anyToDo = true;
    }
    if (useRINS >= kType && useRINS <= kType + 1) {
        CbcHeuristicRINS heuristic5(*model);
        heuristic5.setHeuristicName("RINS");
        if (useRINS < 4) {
            heuristic5.setFractionSmall(0.5);
            heuristic5.setDecayFactor(5.0);
        } else {
            heuristic5.setFractionSmall(0.6);
            heuristic5.setDecayFactor(1.5);
        }
        model->addHeuristic(&heuristic5) ;
        anyToDo = true;
    }
    if (useCombine >= kType && useCombine <= kType + 1) {
        CbcHeuristicLocal heuristic2(*model);
        heuristic2.setHeuristicName("combine solutions");
        heuristic2.setFractionSmall(0.5);
        heuristic2.setSearchType(1);
        model->addHeuristic(&heuristic2);
        anyToDo = true;
    }
    if (useCrossover >= kType && useCrossover <= kType + 1) {
        CbcHeuristicCrossover heuristic2a(*model);
        heuristic2a.setHeuristicName("crossover");
        heuristic2a.setFractionSmall(0.3);
        // just fix at lower
        heuristic2a.setWhen(11);
        model->addHeuristic(&heuristic2a);
        model->setMaximumSavedSolutions(5);
        anyToDo = true;
    }
    int heurSwitches = parameters_[whichParam(CBC_PARAM_INT_HOPTIONS, numberParameters_, parameters_)].intValue() % 100;
    if (heurSwitches) {
        for (int iHeur = 0; iHeur < model->numberHeuristics(); iHeur++) {
            CbcHeuristic * heuristic = model->heuristic(iHeur);
            heuristic->setSwitches(heurSwitches);
        }
    }
    if (type == 2 && anyToDo) {
        // Do heuristics
#ifndef JJF_ONE
        // clean copy
        CbcModel model2(*model);
        // But get rid of heuristics in model
        model->doHeuristicsAtRoot(2);
        if (logLevel <= 1)
            model2.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
        OsiBabSolver defaultC;
        //solver_->setAuxiliaryInfo(&defaultC);
        model2.passInSolverCharacteristics(&defaultC);
        // Save bounds
        int numberColumns = model2.solver()->getNumCols();
        model2.createContinuousSolver();
        bool cleanModel = !model2.numberIntegers() && !model2.numberObjects();
        model2.findIntegers(false);
        int heurOptions = (parameters_[whichParam(CBC_PARAM_INT_HOPTIONS, numberParameters_, parameters_)].intValue() / 100) % 100;
        if (heurOptions == 0 || heurOptions == 2) {
            model2.doHeuristicsAtRoot(1);
        } else if (heurOptions == 1 || heurOptions == 3) {
            model2.setMaximumNodes(-1);
            CbcStrategyDefault strategy(0, 5, 5);
            strategy.setupPreProcessing(1, 0);
            model2.setStrategy(strategy);
            model2.branchAndBound();
        }
        if (cleanModel)
            model2.zapIntegerInformation(false);
        if (model2.bestSolution()) {
            double value = model2.getMinimizationObjValue();
            model->setCutoff(value);
            model->setBestSolution(model2.bestSolution(), numberColumns, value);
            model->setSolutionCount(1);
            model->setNumberHeuristicSolutions(1);
        }
#else
        if (logLevel <= 1)
            model->solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
        OsiBabSolver defaultC;
        //solver_->setAuxiliaryInfo(&defaultC);
        model->passInSolverCharacteristics(&defaultC);
        // Save bounds
        int numberColumns = model->solver()->getNumCols();
        model->createContinuousSolver();
        bool cleanModel = !model->numberIntegers() && !model->numberObjects();
        model->findIntegers(false);
        model->doHeuristicsAtRoot(1);
        if (cleanModel)
            model->zapIntegerInformation(false);
#endif
        return 0;
    } else {
        return 0;
    }
}

