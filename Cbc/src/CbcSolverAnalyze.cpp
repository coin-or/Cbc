/* $Id: CbcSolverAnalyze.cpp 1240 2009-10-02 18:41:44Z forrest $ */
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file CbcSolverAnalyze.cpp

  Look to see if a constraint is all-integer (variables & coeffs), or could be
  all integer. Consider whether one or two continuous variables can be declared
  integer. John's comment is that integer preprocessing might do a better job,
  so we should consider whether this routine should stay.

  No hurry to get rid of it, in my opinion  -- lh, 091210 --

*/

#include "CbcConfig.h"
#include "CoinPragma.hpp"

#include "OsiClpSolverInterface.hpp"

#include "ClpMessage.hpp"

#include "CbcModel.hpp"


#ifndef CBC_OTHER_SOLVER

int * analyze(OsiClpSolverInterface * solverMod, int & numberChanged,
		     double & increment, bool changeInt,
		     CoinMessageHandler * generalMessageHandler, bool noPrinting)
{
    bool noPrinting_ = noPrinting;
    OsiSolverInterface * solver = solverMod->clone();
    char generalPrint[200];
    if (0) {
        // just get increment
        CbcModel model(*solver);
        model.analyzeObjective();
        double increment2 = model.getCutoffIncrement();
        printf("initial cutoff increment %g\n", increment2);
    }
    const double *objective = solver->getObjCoefficients() ;
    const double *lower = solver->getColLower() ;
    const double *upper = solver->getColUpper() ;
    int numberColumns = solver->getNumCols() ;
    int numberRows = solver->getNumRows();
    double direction = solver->getObjSense();
    int iRow, iColumn;

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
    memset(changeRhs, 0, numberRows*sizeof(double));
    memset(ignore, 0, numberRows);
    numberChanged = 0;
    int numberInteger = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (upper[iColumn] > lower[iColumn] + 1.0e-8 && solver->isInteger(iColumn))
            numberInteger++;
    }
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
                    if (fabs(value1) == 1.0 && value1*value2 == -1.0 && !lower[jColumn1]
                            && !lower[jColumn2]) {
                        int n = 0;
                        int i;
                        double objChange = direction * (objective[jColumn1] + objective[jColumn2]);
                        double bound = CoinMin(upper[jColumn1], upper[jColumn2]);
                        bound = CoinMin(bound, 1.0e20);
                        for ( i = columnStart[jColumn1]; i < columnStart[jColumn1] + columnLength[jColumn1]; i++) {
                            int jRow = row[i];
                            double value = element[i];
                            if (jRow != iRow) {
                                which[n++] = jRow;
                                changeRhs[jRow] = value;
                            }
                        }
                        for ( i = columnStart[jColumn1]; i < columnStart[jColumn1] + columnLength[jColumn1]; i++) {
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
    delete [] which;
    delete [] changeRhs;
    delete [] ignore;
    //if (numberInteger&&!noPrinting_)
    //printf("%d integer variables",numberInteger);
    if (changeInt) {
        //if (!noPrinting_) {
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
            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                if (solver->isInteger(iColumn))
                    solverMod->setInteger(iColumn);
            }
            delete solver;
            return changed;
        }
    } else {
        //if (!noPrinting_) {
        //if (numberChanged)
        //  printf(" and %d variables could be made integer\n",numberChanged);
        //else
        //  printf("\n");
        //}
        // just get increment
        int logLevel = generalMessageHandler->logLevel();
        CbcModel model(*solver);
        model.passInMessageHandler(generalMessageHandler);
        if (noPrinting_)
            model.setLogLevel(0);
        model.analyzeObjective();
        generalMessageHandler->setLogLevel(logLevel);
        double increment2 = model.getCutoffIncrement();
        if (increment2 > increment && increment2 > 0.0) {
            if (!noPrinting_) {
                sprintf(generalPrint, "Cutoff increment increased from %g to %g", increment, increment2);
                CoinMessages generalMessages = solverMod->getModelPtr()->messages();
                generalMessageHandler->message(CLP_GENERAL, generalMessages)
                << generalPrint
                << CoinMessageEol;
            }
            increment = increment2;
        }
        delete solver;
        numberChanged = 0;
        delete [] changed;
        return NULL;
    }
}
#endif	// ifndef CBC_OTHER_SOLVER

