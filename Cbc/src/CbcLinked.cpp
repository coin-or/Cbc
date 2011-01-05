/* $Id: CbcLinked.cpp 1200 2009-07-25 08:44:13Z forrest $ */
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcConfig.h"

#include "CoinTime.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinModel.hpp"
#include "ClpSimplex.hpp"
// returns jColumn (-2 if linear term, -1 if unknown) and coefficient
static
int decodeBit(char * phrase, char * & nextPhrase, double & coefficient, bool ifFirst, const CoinModel & model)
{
    char * pos = phrase;
    // may be leading - (or +)
    char * pos2 = pos;
    double value = 1.0;
    if (*pos2 == '-' || *pos2 == '+')
        pos2++;
    // next terminator * or + or -
    while (*pos2) {
        if (*pos2 == '*') {
            break;
        } else if (*pos2 == '-' || *pos2 == '+') {
            if (pos2 == pos || *(pos2 - 1) != 'e')
                break;
        }
        pos2++;
    }
    // if * must be number otherwise must be name
    if (*pos2 == '*') {
        char * pos3 = pos;
        while (pos3 != pos2) {
            pos3++;
#ifndef NDEBUG
            char x = *(pos3 - 1);
            assert ((x >= '0' && x <= '9') || x == '.' || x == '+' || x == '-' || x == 'e');
#endif
        }
        char saved = *pos2;
        *pos2 = '\0';
        value = atof(pos);
        *pos2 = saved;
        // and down to next
        pos2++;
        pos = pos2;
        while (*pos2) {
            if (*pos2 == '-' || *pos2 == '+')
                break;
            pos2++;
        }
    }
    char saved = *pos2;
    *pos2 = '\0';
    // now name
    // might have + or -
    if (*pos == '+') {
        pos++;
    } else if (*pos == '-') {
        pos++;
        assert (value == 1.0);
        value = - value;
    }
    int jColumn = model.column(pos);
    // must be column unless first when may be linear term
    if (jColumn < 0) {
        if (ifFirst) {
            char * pos3 = pos;
            while (pos3 != pos2) {
                pos3++;
#ifndef NDEBUG
                char x = *(pos3 - 1);
                assert ((x >= '0' && x <= '9') || x == '.' || x == '+' || x == '-' || x == 'e');
#endif
            }
            assert(*pos2 == '\0');
            // keep possible -
            value = value * atof(pos);
            jColumn = -2;
        } else {
            // bad
            *pos2 = saved;
            printf("bad nonlinear term %s\n", phrase);
            abort();
        }
    }
    *pos2 = saved;
    pos = pos2;
    coefficient = value;
    nextPhrase = pos;
    return jColumn;
}
#include "ClpQuadraticObjective.hpp"
#include <cassert>
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include "CbcLinked.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinMpsIO.hpp"
//#include "OsiSolverLink.hpp"
//#include "OsiBranchLink.hpp"
#include "ClpPackedMatrix.hpp"
#include "CoinTime.hpp"
#include "CbcModel.hpp"
#include "CbcCutGenerator.hpp"
#include "CglStored.hpp"
#include "CglPreProcess.hpp"
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglRedSplit.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglTwomir.hpp"
#include "CglDuplicateRow.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "ClpLinearObjective.hpp"
#include "CbcBranchActual.hpp"
#include "CbcCompareActual.hpp"
//#############################################################################
// Solve methods
//#############################################################################
void OsiSolverLink::initialSolve()
{
    //writeMps("yy");
    //exit(77);
    specialOptions_ = 0;
    modelPtr_->setWhatsChanged(0);
    if (numberVariables_) {
        CoinPackedMatrix * temp = new CoinPackedMatrix(*matrix_);
        // update all bounds before coefficients
        for (int i = 0; i < numberVariables_; i++ ) {
            info_[i].updateBounds(modelPtr_);
        }
        int updated = updateCoefficients(modelPtr_, temp);
        if (updated || 1) {
            temp->removeGaps(1.0e-14);
            ClpMatrixBase * save = modelPtr_->clpMatrix();
            ClpPackedMatrix * clpMatrix = dynamic_cast<ClpPackedMatrix *> (save);
            assert (clpMatrix);
            if (save->getNumRows() > temp->getNumRows()) {
                // add in cuts
                int numberRows = temp->getNumRows();
                int * which = new int[numberRows];
                for (int i = 0; i < numberRows; i++)
                    which[i] = i;
                save->deleteRows(numberRows, which);
                delete [] which;
                temp->bottomAppendPackedMatrix(*clpMatrix->matrix());
            }
            modelPtr_->replaceMatrix(temp, true);
        } else {
            delete temp;
        }
    }
    if (0) {
        const double * lower = getColLower();
        const double * upper = getColUpper();
        int n = 0;
        for (int i = 84; i < 84 + 16; i++) {
            if (lower[i] + 0.01 < upper[i]) {
                n++;
            }
        }
        if (!n)
            writeMps("sol_query");

    }
    //static int iPass=0;
    //char temp[50];
    //iPass++;
    //sprintf(temp,"cc%d",iPass);
    //writeMps(temp);
    //writeMps("tight");
    //exit(33);
    //printf("wrote cc%d\n",iPass);
    OsiClpSolverInterface::initialSolve();
    int secondaryStatus = modelPtr_->secondaryStatus();
    if (modelPtr_->status() == 0 && (secondaryStatus == 2 || secondaryStatus == 4))
        modelPtr_->cleanup(1);
    //if (!isProvenOptimal())
    //writeMps("yy");
    if (isProvenOptimal() && quadraticModel_ && modelPtr_->numberColumns() == quadraticModel_->numberColumns()) {
        // see if qp can get better solution
        const double * solution = modelPtr_->primalColumnSolution();
        int numberColumns = modelPtr_->numberColumns();
        bool satisfied = true;
        for (int i = 0; i < numberColumns; i++) {
            if (isInteger(i)) {
                double value = solution[i];
                if (fabs(value - floor(value + 0.5)) > 1.0e-6) {
                    satisfied = false;
                    break;
                }
            }
        }
        if (satisfied) {
            ClpSimplex qpTemp(*quadraticModel_);
            double * lower = qpTemp.columnLower();
            double * upper = qpTemp.columnUpper();
            double * lower2 = modelPtr_->columnLower();
            double * upper2 = modelPtr_->columnUpper();
            for (int i = 0; i < numberColumns; i++) {
                if (isInteger(i)) {
                    double value = floor(solution[i] + 0.5);
                    lower[i] = value;
                    upper[i] = value;
                } else {
                    lower[i] = lower2[i];
                    upper[i] = upper2[i];
                }
            }
            //qpTemp.writeMps("bad.mps");
            //modelPtr_->writeMps("bad2.mps");
            //qpTemp.objectiveAsObject()->setActivated(0);
            //qpTemp.primal();
            //qpTemp.objectiveAsObject()->setActivated(1);
            qpTemp.primal();
            //assert (!qpTemp.problemStatus());
            if (qpTemp.objectiveValue() < bestObjectiveValue_ - 1.0e-3 && !qpTemp.problemStatus()) {
                delete [] bestSolution_;
                bestSolution_ = CoinCopyOfArray(qpTemp.primalColumnSolution(), numberColumns);
                bestObjectiveValue_ = qpTemp.objectiveValue();
                printf("better qp objective of %g\n", bestObjectiveValue_);
                // If model has stored then add cut (if convex)
                if (cbcModel_ && (specialOptions2_&4) != 0) {
                    int numberGenerators = cbcModel_->numberCutGenerators();
                    int iGenerator;
                    cbcModel_->lockThread();
                    for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
                        CbcCutGenerator * generator = cbcModel_->cutGenerator(iGenerator);
                        CglCutGenerator * gen = generator->generator();
                        CglStored * gen2 = dynamic_cast<CglStored *> (gen);
                        if (gen2) {
                            // add OA cut
                            double offset;
                            double * gradient = new double [numberColumns+1];
                            memcpy(gradient, qpTemp.objectiveAsObject()->gradient(&qpTemp, bestSolution_, offset, true, 2),
                                   numberColumns*sizeof(double));
                            // assume convex
                            double rhs = 0.0;
                            int * column = new int[numberColumns+1];
                            int n = 0;
                            for (int i = 0; i < numberColumns; i++) {
                                double value = gradient[i];
                                if (fabs(value) > 1.0e-12) {
                                    gradient[n] = value;
                                    rhs += value * solution[i];
                                    column[n++] = i;
                                }
                            }
                            gradient[n] = -1.0;
                            column[n++] = objectiveVariable_;
                            gen2->addCut(-COIN_DBL_MAX, offset + 1.0e-7, n, column, gradient);
                            delete [] gradient;
                            delete [] column;
                            break;
                        }
                    }
                    cbcModel_->unlockThread();
                }
            }
        }
    }
}
//#define WRITE_MATRIX
#ifdef WRITE_MATRIX
static int xxxxxx = 0;
#endif
//-----------------------------------------------------------------------------
void OsiSolverLink::resolve()
{
    if (false) {
        bool takeHint;
        OsiHintStrength strength;
        // Switch off printing if asked to
        getHintParam(OsiDoReducePrint, takeHint, strength);
        if (strength != OsiHintIgnore && takeHint) {
            printf("no printing\n");
        } else {
            printf("printing\n");
        }
    }
    specialOptions_ = 0;
    modelPtr_->setWhatsChanged(0);
    bool allFixed = numberFix_ > 0;
    bool feasible = true;
    if (numberVariables_) {
        CoinPackedMatrix * temp = new CoinPackedMatrix(*matrix_);
        //bool best=true;
        const double * lower = modelPtr_->columnLower();
        const double * upper = modelPtr_->columnUpper();
        // update all bounds before coefficients
        for (int i = 0; i < numberVariables_; i++ ) {
            info_[i].updateBounds(modelPtr_);
            int iColumn = info_[i].variable();
            double lo = lower[iColumn];
            double up = upper[iColumn];
            if (up < lo)
                feasible = false;
        }
        int updated = updateCoefficients(modelPtr_, temp);
        if (updated) {
            temp->removeGaps(1.0e-14);
            ClpMatrixBase * save = modelPtr_->clpMatrix();
            ClpPackedMatrix * clpMatrix = dynamic_cast<ClpPackedMatrix *> (save);
            assert (clpMatrix);
            if (save->getNumRows() > temp->getNumRows()) {
                // add in cuts
                int numberRows = temp->getNumRows();
                int * which = new int[numberRows];
                for (int i = 0; i < numberRows; i++)
                    which[i] = i;
                CoinPackedMatrix * mat = clpMatrix->matrix();
                // for debug
                //mat = new CoinPackedMatrix(*mat);
                mat->deleteRows(numberRows, which);
                delete [] which;
                temp->bottomAppendPackedMatrix(*mat);
                temp->removeGaps(1.0e-14);
            }
            modelPtr_->replaceMatrix(temp, true);
        } else {
            delete temp;
        }
    }
#ifdef WRITE_MATRIX
    {
        xxxxxx++;
        char temp[50];
        sprintf(temp, "bb%d", xxxxxx);
        writeMps(temp);
        printf("wrote bb%d\n", xxxxxx);
    }
#endif
    if (0) {
        const double * lower = getColLower();
        const double * upper = getColUpper();
        int n = 0;
        for (int i = 60; i < 64; i++) {
            if (lower[i]) {
                printf("%d bounds %g %g\n", i, lower[i], upper[i]);
                n++;
            }
        }
        if (n == 1)
            printf("just one?\n");
    }
    // check feasible
    {
        const double * lower = getColLower();
        const double * upper = getColUpper();
        int numberColumns = getNumCols();
        for (int i = 0; i < numberColumns; i++) {
            if (lower[i] > upper[i] + 1.0e-12) {
                feasible = false;
                break;
            }
        }
    }
    if (!feasible)
        allFixed = false;
    if ((specialOptions2_&1) == 0)
        allFixed = false;
    int returnCode = -1;
    // See if in strong branching
    int maxIts = modelPtr_->maximumIterations();
    if (feasible) {
        if (maxIts > 10000) {
            // may do lots of work
            if ((specialOptions2_&1) != 0) {
                // see if fixed
                const double * lower = modelPtr_->columnLower();
                const double * upper = modelPtr_->columnUpper();
                for (int i = 0; i < numberFix_; i++ ) {
                    int iColumn = fixVariables_[i];
                    double lo = lower[iColumn];
                    double up = upper[iColumn];
                    if (up > lo) {
                        allFixed = false;
                        break;
                    }
                }
                returnCode = allFixed ? fathom(allFixed) : 0;
            } else {
                returnCode = 0;
            }
        } else {
            returnCode = 0;
        }
    }
    if (returnCode >= 0) {
        if (returnCode == 0)
            OsiClpSolverInterface::resolve();
        if (!allFixed && (specialOptions2_&1) != 0) {
            const double * solution = getColSolution();
            bool satisfied = true;
            for (int i = 0; i < numberVariables_; i++) {
                int iColumn = info_[i].variable();
                double value = solution[iColumn];
                if (fabs(value - floor(value + 0.5)) > 0.0001)
                    satisfied = false;
            }
            //if (satisfied)
            //printf("satisfied but not fixed\n");
        }
        int satisfied = 2;
        const double * solution = getColSolution();
        const double * lower = getColLower();
        const double * upper = getColUpper();
        int numberColumns2 = coinModel_.numberColumns();
        for (int i = 0; i < numberColumns2; i++) {
            if (isInteger(i)) {
                double value = solution[i];
                if (fabs(value - floor(value + 0.5)) > 1.0e-6) {
                    satisfied = 0;
                    break;
                } else if (upper[i] > lower[i]) {
                    satisfied = 1;
                }
            }
        }
        if (isProvenOptimal()) {
            //if (satisfied==2)
            //printf("satisfied %d\n",satisfied);
            if (satisfied && (specialOptions2_&2) != 0) {
                assert (quadraticModel_);
                // look at true objective
#ifndef NDEBUG
                double direction = modelPtr_->optimizationDirection();
                assert (direction == 1.0);
#endif
                double value = - quadraticModel_->objectiveOffset();
                const double * objective = quadraticModel_->objective();
                int i;
                for ( i = 0; i < numberColumns2; i++)
                    value += solution[i] * objective[i];
                // and now rest
                for ( i = 0; i < numberObjects_; i++) {
                    OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[i]);
                    if (obj) {
                        value += obj->xyCoefficient(solution);
                    }
                }
                if (value < bestObjectiveValue_ - 1.0e-3) {
                    printf("obj of %g\n", value);
                    //modelPtr_->setDualObjectiveLimit(value);
                    delete [] bestSolution_;
                    bestSolution_ = CoinCopyOfArray(modelPtr_->getColSolution(), modelPtr_->getNumCols());
                    bestObjectiveValue_ = value;
                    if (maxIts <= 10000 && cbcModel_) {
                        OsiSolverLink * solver2 = dynamic_cast<OsiSolverLink *> (cbcModel_->solver());
                        assert (solver2);
                        if (solver2 != this) {
                            // in strong branching - need to store in original solver
                            if (value < solver2->bestObjectiveValue_ - 1.0e-3) {
                                delete [] solver2->bestSolution_;
                                solver2->bestSolution_ = CoinCopyOfArray(bestSolution_, modelPtr_->getNumCols());
                                solver2->bestObjectiveValue_ = value;
                            }
                        }
                    }
                    // If model has stored then add cut (if convex)
                    if (cbcModel_ && (specialOptions2_&4) != 0 && quadraticModel_) {
                        int numberGenerators = cbcModel_->numberCutGenerators();
                        int iGenerator;
                        for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
                            CbcCutGenerator * generator = cbcModel_->cutGenerator(iGenerator);
                            CglCutGenerator * gen = generator->generator();
                            CglStored * gen2 = dynamic_cast<CglStored *> (gen);
                            if (gen2) {
                                cbcModel_->lockThread();
                                // add OA cut
                                double offset = 0.0;
                                int numberColumns = quadraticModel_->numberColumns();
                                double * gradient = new double [numberColumns+1];
                                // gradient from bilinear
                                int i;
                                CoinZeroN(gradient, numberColumns + 1);
                                //const double * objective = modelPtr_->objective();
                                assert (objectiveRow_ >= 0);
                                const double * element = originalRowCopy_->getElements();
                                const int * column2 = originalRowCopy_->getIndices();
                                const CoinBigIndex * rowStart = originalRowCopy_->getVectorStarts();
                                //const int * rowLength = originalRowCopy_->getVectorLengths();
                                //int numberColumns2 = coinModel_.numberColumns();
                                for ( i = rowStart[objectiveRow_]; i < rowStart[objectiveRow_+1]; i++)
                                    gradient[column2[i]] = element[i];
                                for ( i = 0; i < numberObjects_; i++) {
                                    OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[i]);
                                    if (obj) {
                                        int xColumn = obj->xColumn();
                                        int yColumn = obj->yColumn();
                                        if (xColumn != yColumn) {
                                            double coefficient = /* 2.0* */obj->coefficient();
                                            gradient[xColumn] += coefficient * solution[yColumn];
                                            gradient[yColumn] += coefficient * solution[xColumn];
                                            offset += coefficient * solution[xColumn] * solution[yColumn];
                                        } else {
                                            double coefficient = obj->coefficient();
                                            gradient[xColumn] += 2.0 * coefficient * solution[yColumn];
                                            offset += coefficient * solution[xColumn] * solution[yColumn];
                                        }
                                    }
                                }
                                // assume convex
                                double rhs = 0.0;
                                int * column = new int[numberColumns+1];
                                int n = 0;
                                for (int i = 0; i < numberColumns; i++) {
                                    double value = gradient[i];
                                    if (fabs(value) > 1.0e-12) {
                                        gradient[n] = value;
                                        rhs += value * solution[i];
                                        column[n++] = i;
                                    }
                                }
                                gradient[n] = -1.0;
                                column[n++] = objectiveVariable_;
                                gen2->addCut(-COIN_DBL_MAX, offset + 1.0e-4, n, column, gradient);
                                delete [] gradient;
                                delete [] column;
                                cbcModel_->unlockThread();
                                break;
                            }
                        }
                    }
                }
            } else if (satisfied == 2) {
                // is there anything left to do?
                int i;
                int numberContinuous = 0;
                double gap = 0.0;
                for ( i = 0; i < numberObjects_; i++) {
                    OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[i]);
                    if (obj) {
                        if (obj->xMeshSize() < 1.0 && obj->yMeshSize() < 1.0) {
                            numberContinuous++;
                            int xColumn = obj->xColumn();
                            double gapX = upper[xColumn] - lower[xColumn];
                            int yColumn = obj->yColumn();
                            double gapY = upper[yColumn] - lower[yColumn];
                            gap = CoinMax(gap, CoinMax(gapX, gapY));
                        }
                    }
                }
                if (numberContinuous && 0) {
                    // iterate to get solution and fathom node
                    int numberColumns2 = coinModel_.numberColumns();
                    double * lower2 = CoinCopyOfArray(getColLower(), numberColumns2);
                    double * upper2 = CoinCopyOfArray(getColUpper(), numberColumns2);
                    while (gap > defaultMeshSize_) {
                        gap *= 0.9;
                        const double * solution = getColSolution();
                        double newGap = 0.0;
                        for ( i = 0; i < numberObjects_; i++) {
                            OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[i]);
                            if (obj && (obj->branchingStrategy()&8) == 0) {
                                if (obj->xMeshSize() < 1.0 && obj->yMeshSize() < 1.0) {
                                    numberContinuous++;
                                    // need to make sure new xy value in range
                                    double xB[3];
                                    double yB[3];
                                    double xybar[4];
                                    obj->getCoefficients(this, xB, yB, xybar);
                                    //double xyTrue = obj->xyCoefficient(solution);
                                    double xyLambda = 0.0;
                                    int firstLambda = obj->firstLambda();
                                    for (int j = 0; j < 4; j++) {
                                        xyLambda += solution[firstLambda+j] * xybar[j];
                                    }
                                    //printf ("x %d, y %d - true %g lambda %g\n",obj->xColumn(),obj->yColumn(),
                                    //  xyTrue,xyLambda);
                                    int xColumn = obj->xColumn();
                                    double gapX = upper[xColumn] - lower[xColumn];
                                    int yColumn = obj->yColumn();
                                    if (gapX > gap) {
                                        double value = solution[xColumn];
                                        double newLower = CoinMax(lower2[xColumn], value - 0.5 * gap);
                                        double newUpper = CoinMin(upper2[xColumn], value + 0.5 * gap);
                                        if (newUpper - newLower < 0.99*gap) {
                                            if (newLower == lower2[xColumn])
                                                newUpper = CoinMin(upper2[xColumn], newLower + gap);
                                            else if (newUpper == upper2[xColumn])
                                                newLower = CoinMax(lower2[xColumn], newUpper - gap);
                                        }
                                        // see if problem
#ifndef NDEBUG
                                        double lambda[4];
#endif
                                        xB[0] = newLower;
                                        xB[1] = newUpper;
                                        xB[2] = value;
                                        yB[2] = solution[yColumn];
                                        xybar[0] = xB[0] * yB[0];
                                        xybar[1] = xB[0] * yB[1];
                                        xybar[2] = xB[1] * yB[0];
                                        xybar[3] = xB[1] * yB[1];
#ifndef NDEBUG
                                        double infeasibility = obj->computeLambdas(xB, yB, xybar, lambda);
                                        assert (infeasibility < 1.0e-9);
#endif
                                        setColLower(xColumn, newLower);
                                        setColUpper(xColumn, newUpper);
                                    }
                                    double gapY = upper[yColumn] - lower[yColumn];
                                    if (gapY > gap) {
                                        double value = solution[yColumn];
                                        double newLower = CoinMax(lower2[yColumn], value - 0.5 * gap);
                                        double newUpper = CoinMin(upper2[yColumn], value + 0.5 * gap);
                                        if (newUpper - newLower < 0.99*gap) {
                                            if (newLower == lower2[yColumn])
                                                newUpper = CoinMin(upper2[yColumn], newLower + gap);
                                            else if (newUpper == upper2[yColumn])
                                                newLower = CoinMax(lower2[yColumn], newUpper - gap);
                                        }
                                        // see if problem
#ifndef NDEBUG
                                        double lambda[4];
#endif
                                        yB[0] = newLower;
                                        yB[1] = newUpper;
                                        xybar[0] = xB[0] * yB[0];
                                        xybar[1] = xB[0] * yB[1];
                                        xybar[2] = xB[1] * yB[0];
                                        xybar[3] = xB[1] * yB[1];
#ifndef NDEBUG
                                        double infeasibility = obj->computeLambdas(xB, yB, xybar, lambda);
                                        assert (infeasibility < 1.0e-9);
#endif
                                        setColLower(yColumn, newLower);
                                        setColUpper(yColumn, newUpper);
                                    }
                                    newGap = CoinMax(newGap, CoinMax(gapX, gapY));
                                }
                            }
                        }
                        printf("solving with gap of %g\n", gap);
                        //OsiClpSolverInterface::resolve();
                        initialSolve();
                        if (!isProvenOptimal())
                            break;
                    }
                    delete [] lower2;
                    delete [] upper2;
                    //if (isProvenOptimal())
                    //writeMps("zz");
                }
            }
            // ???  - try
            // But skip if strong branching
            CbcModel * cbcModel = (modelPtr_->maximumIterations() > 10000) ? cbcModel_ : NULL;
            if ((specialOptions2_&2) != 0) {
                // If model has stored then add cut (if convex)
                // off until I work out problem with ibell3a
                if (cbcModel && (specialOptions2_&4) != 0 && quadraticModel_) {
                    int numberGenerators = cbcModel_->numberCutGenerators();
                    int iGenerator;
                    for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
                        CbcCutGenerator * generator = cbcModel_->cutGenerator(iGenerator);
                        CglCutGenerator * gen = generator->generator();
                        CglTemporary * gen2 = dynamic_cast<CglTemporary *> (gen);
                        if (gen2) {
                            double * solution2 = NULL;
                            int numberColumns = quadraticModel_->numberColumns();
                            int depth = cbcModel_->currentNode() ? cbcModel_->currentNode()->depth() : 0;
                            if (depth < 5) {
                                ClpSimplex qpTemp(*quadraticModel_);
                                double * lower = qpTemp.columnLower();
                                double * upper = qpTemp.columnUpper();
                                double * lower2 = modelPtr_->columnLower();
                                double * upper2 = modelPtr_->columnUpper();
                                for (int i = 0; i < numberColumns; i++) {
                                    lower[i] = lower2[i];
                                    upper[i] = upper2[i];
                                }
                                qpTemp.setLogLevel(modelPtr_->logLevel());
                                qpTemp.primal();
                                assert (!qpTemp.problemStatus());
                                if (qpTemp.objectiveValue() < bestObjectiveValue_ - 1.0e-3 && !qpTemp.problemStatus()) {
                                    solution2 = CoinCopyOfArray(qpTemp.primalColumnSolution(), numberColumns);
                                } else {
                                    printf("QP says expensive - kill\n");
                                    modelPtr_->setProblemStatus(1);
                                    modelPtr_->setObjectiveValue(COIN_DBL_MAX);
                                    break;
                                }
                            }
                            const double * solution = getColSolution();
                            // add OA cut
                            doAOCuts(gen2, solution, solution);
                            if (solution2) {
                                doAOCuts(gen2, solution, solution2);
                                delete [] solution2;
                            }
                            break;
                        }
                    }
                }
            } else if (cbcModel && (specialOptions2_&8) == 8) {
                // convex and nonlinear in constraints
                int numberGenerators = cbcModel_->numberCutGenerators();
                int iGenerator;
                for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
                    CbcCutGenerator * generator = cbcModel_->cutGenerator(iGenerator);
                    CglCutGenerator * gen = generator->generator();
                    CglTemporary * gen2 = dynamic_cast<CglTemporary *> (gen);
                    if (gen2) {
                        const double * solution = getColSolution();
                        const double * rowUpper = getRowUpper();
                        const double * rowLower = getRowLower();
                        const double * element = originalRowCopy_->getElements();
                        const int * column2 = originalRowCopy_->getIndices();
                        const CoinBigIndex * rowStart = originalRowCopy_->getVectorStarts();
                        //const int * rowLength = originalRowCopy_->getVectorLengths();
                        int numberColumns2 = CoinMax(coinModel_.numberColumns(), objectiveVariable_ + 1);
                        double * gradient = new double [numberColumns2];
                        int * column = new int[numberColumns2];
                        //const double * columnLower = modelPtr_->columnLower();
                        //const double * columnUpper = modelPtr_->columnUpper();
                        cbcModel_->lockThread();
                        for (int iNon = 0; iNon < numberNonLinearRows_; iNon++) {
                            int iRow = rowNonLinear_[iNon];
                            bool convex = convex_[iNon] > 0;
                            if (!convex_[iNon])
                                continue; // can't use this row
                            // add OA cuts
                            double offset = 0.0;
                            // gradient from bilinear
                            int i;
                            CoinZeroN(gradient, numberColumns2);
                            //const double * objective = modelPtr_->objective();
                            for ( i = rowStart[iRow]; i < rowStart[iRow+1]; i++)
                                gradient[column2[i]] = element[i];
                            for ( i = startNonLinear_[iNon]; i < startNonLinear_[iNon+1]; i++) {
                                OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[whichNonLinear_[i]]);
                                assert (obj);
                                int xColumn = obj->xColumn();
                                int yColumn = obj->yColumn();
                                if (xColumn != yColumn) {
                                    double coefficient = /* 2.0* */obj->coefficient();
                                    gradient[xColumn] += coefficient * solution[yColumn];
                                    gradient[yColumn] += coefficient * solution[xColumn];
                                    offset += coefficient * solution[xColumn] * solution[yColumn];
                                } else {
                                    double coefficient = obj->coefficient();
                                    gradient[xColumn] += 2.0 * coefficient * solution[yColumn];
                                    offset += coefficient * solution[xColumn] * solution[yColumn];
                                }
                            }
                            // assume convex
                            double rhs = 0.0;
                            int n = 0;
                            for (int i = 0; i < numberColumns2; i++) {
                                double value = gradient[i];
                                if (fabs(value) > 1.0e-12) {
                                    gradient[n] = value;
                                    rhs += value * solution[i];
                                    column[n++] = i;
                                }
                            }
                            if (iRow == objectiveRow_) {
                                gradient[n] = -1.0;
                                assert (objectiveVariable_ >= 0);
                                rhs -= solution[objectiveVariable_];
                                column[n++] = objectiveVariable_;
                                assert (convex);
                            } else if (convex) {
                                offset += rowUpper[iRow];
                            } else if (!convex) {
                                offset += rowLower[iRow];
                            }
                            if (convex && rhs > offset + 1.0e-5)
                                gen2->addCut(-COIN_DBL_MAX, offset + 1.0e-7, n, column, gradient);
                            else if (!convex && rhs < offset - 1.0e-5)
                                gen2->addCut(offset - 1.0e-7, COIN_DBL_MAX, n, column, gradient);
                        }
                        cbcModel_->unlockThread();
                        delete [] gradient;
                        delete [] column;
                        break;
                    }
                }
            }
        }
    } else {
        modelPtr_->setProblemStatus(1);
        modelPtr_->setObjectiveValue(COIN_DBL_MAX);
    }
}
// Do OA cuts
int
OsiSolverLink::doAOCuts(CglTemporary * cutGen, const double * solution, const double * solution2)
{
    cbcModel_->lockThread();
    // add OA cut
    double offset = 0.0;
    int numberColumns = quadraticModel_->numberColumns();
    double * gradient = new double [numberColumns+1];
    // gradient from bilinear
    int i;
    CoinZeroN(gradient, numberColumns + 1);
    //const double * objective = modelPtr_->objective();
    assert (objectiveRow_ >= 0);
    const double * element = originalRowCopy_->getElements();
    const int * column2 = originalRowCopy_->getIndices();
    const CoinBigIndex * rowStart = originalRowCopy_->getVectorStarts();
    //const int * rowLength = originalRowCopy_->getVectorLengths();
    //int numberColumns2 = coinModel_.numberColumns();
    for ( i = rowStart[objectiveRow_]; i < rowStart[objectiveRow_+1]; i++)
        gradient[column2[i]] = element[i];
    //const double * columnLower = modelPtr_->columnLower();
    //const double * columnUpper = modelPtr_->columnUpper();
    for ( i = 0; i < numberObjects_; i++) {
        OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[i]);
        if (obj) {
            int xColumn = obj->xColumn();
            int yColumn = obj->yColumn();
            if (xColumn != yColumn) {
                double coefficient = /* 2.0* */obj->coefficient();
                gradient[xColumn] += coefficient * solution2[yColumn];
                gradient[yColumn] += coefficient * solution2[xColumn];
                offset += coefficient * solution2[xColumn] * solution2[yColumn];
            } else {
                double coefficient = obj->coefficient();
                gradient[xColumn] += 2.0 * coefficient * solution2[yColumn];
                offset += coefficient * solution2[xColumn] * solution2[yColumn];
            }
        }
    }
    // assume convex
    double rhs = 0.0;
    int * column = new int[numberColumns+1];
    int n = 0;
    for (int i = 0; i < numberColumns; i++) {
        double value = gradient[i];
        if (fabs(value) > 1.0e-12) {
            gradient[n] = value;
            rhs += value * solution[i];
            column[n++] = i;
        }
    }
    gradient[n] = -1.0;
    assert (objectiveVariable_ >= 0);
    rhs -= solution[objectiveVariable_];
    column[n++] = objectiveVariable_;
    int returnCode = 0;
    if (rhs > offset + 1.0e-5) {
        cutGen->addCut(-COIN_DBL_MAX, offset + 1.0e-7, n, column, gradient);
        //printf("added cut with %d elements\n",n);
        returnCode = 1;
    }
    delete [] gradient;
    delete [] column;
    cbcModel_->unlockThread();
    return returnCode;
}

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
OsiSolverLink::OsiSolverLink ()
        : CbcOsiSolver()
{
    gutsOfDestructor(true);
}
#ifdef JJF_ZERO
/* returns
   sequence of nonlinear or
   -1 numeric
   -2 not found
   -3 too many terms
*/
static int getVariable(const CoinModel & model, char * expression,
                       int & linear)
{
    int non = -1;
    linear = -1;
    if (strcmp(expression, "Numeric")) {
        // function
        char * first = strchr(expression, '*');
        int numberColumns = model.numberColumns();
        int j;
        if (first) {
            *first = '\0';
            for (j = 0; j < numberColumns; j++) {
                if (!strcmp(expression, model.columnName(j))) {
                    linear = j;
                    memmove(expression, first + 1, strlen(first + 1) + 1);
                    break;
                }
            }
        }
        // find nonlinear
        for (j = 0; j < numberColumns; j++) {
            const char * name = model.columnName(j);
            first = strstr(expression, name);
            if (first) {
                if (first != expression && isalnum(*(first - 1)))
                    continue; // not real match
                first += strlen(name);
                if (!isalnum(*first)) {
                    // match
                    non = j;
                    // but check no others
                    j++;
                    for (; j < numberColumns; j++) {
                        const char * name = model.columnName(j);
                        first = strstr(expression, name);
                        if (first) {
                            if (isalnum(*(first - 1)))
                                continue; // not real match
                            first += strlen(name);
                            if (!isalnum(*first)) {
                                // match - ouch
                                non = -3;
                                break;
                            }
                        }
                    }
                    break;
                }
            }
        }
        if (non == -1)
            non = -2;
    }
    return non;
}
#endif
/* This creates from a coinModel object

   if errors.then number of sets is -1

   This creates linked ordered sets information.  It assumes -

   for product terms syntax is yy*f(zz)
   also just f(zz) is allowed
   and even a constant

   modelObject not const as may be changed as part of process.
*/
OsiSolverLink::OsiSolverLink ( CoinModel & coinModel)
        : CbcOsiSolver()
{
    gutsOfDestructor(true);
    load(coinModel);
}
// need bounds
static void fakeBounds(OsiSolverInterface * solver, int column, double maximumValue,
                       CoinModel * model1, CoinModel * model2)
{
    double lo = solver->getColLower()[column];
    if (lo < -maximumValue) {
        solver->setColLower(column, -maximumValue);
        if (model1)
            model1->setColLower(column, -maximumValue);
        if (model2)
            model2->setColLower(column, -maximumValue);
    }
    double up = solver->getColUpper()[column];
    if (up > maximumValue) {
        solver->setColUpper(column, maximumValue);
        if (model1)
            model1->setColUpper(column, maximumValue);
        if (model2)
            model2->setColUpper(column, maximumValue);
    }
}
void OsiSolverLink::load ( CoinModel & coinModelOriginal, bool tightenBounds, int logLevel)
{
    // first check and set up arrays
    int numberColumns = coinModelOriginal.numberColumns();
    int numberRows = coinModelOriginal.numberRows();
    // List of nonlinear entries
    int * which = new int[numberColumns];
    numberVariables_ = 0;
    //specialOptions2_=0;
    int iColumn;
    int numberErrors = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        CoinModelLink triple = coinModelOriginal.firstInColumn(iColumn);
        bool linear = true;
        int n = 0;
        // See if quadratic objective
        const char * expr = coinModelOriginal.getColumnObjectiveAsString(iColumn);
        if (strcmp(expr, "Numeric")) {
            linear = false;
        }
        while (triple.row() >= 0) {
            int iRow = triple.row();
            const char * expr = coinModelOriginal.getElementAsString(iRow, iColumn);
            if (strcmp(expr, "Numeric")) {
                linear = false;
            }
            triple = coinModelOriginal.next(triple);
            n++;
        }
        if (!linear) {
            which[numberVariables_++] = iColumn;
        }
    }
    // return if nothing
    if (!numberVariables_) {
        delete [] which;
        coinModel_ = coinModelOriginal;
        int nInt = 0;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (coinModel_.isInteger(iColumn))
                nInt++;
        }
        printf("There are %d integers\n", nInt);
        loadFromCoinModel(coinModelOriginal, true);
        OsiObject ** objects = new OsiObject * [nInt];
        nInt = 0;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (coinModel_.isInteger(iColumn)) {
                objects[nInt] = new OsiSimpleInteger(this, iColumn);
                objects[nInt]->setPriority(integerPriority_);
                nInt++;
            }
        }
        addObjects(nInt, objects);
        int i;
        for (i = 0; i < nInt; i++)
            delete objects[i];
        delete [] objects;
        return;
    } else {
        coinModel_ = coinModelOriginal;
        // arrays for tightening bounds
        int * freeRow = new int [numberRows];
        CoinZeroN(freeRow, numberRows);
        int * tryColumn = new int [numberColumns];
        CoinZeroN(tryColumn, numberColumns);
        int nBi = 0;
        int numberQuadratic = 0;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            // See if quadratic objective
            const char * expr = coinModel_.getColumnObjectiveAsString(iColumn);
            if (strcmp(expr, "Numeric")) {
                // check if value*x+-value*y....
                assert (strlen(expr) < 20000);
                tryColumn[iColumn] = 1;
                char temp[20000];
                strcpy(temp, expr);
                char * pos = temp;
                bool ifFirst = true;
                double linearTerm = 0.0;
                while (*pos) {
                    double value;
                    int jColumn = decodeBit(pos, pos, value, ifFirst, coinModel_);
                    // must be column unless first when may be linear term
                    if (jColumn >= 0) {
                        tryColumn[jColumn] = 1;
                        numberQuadratic++;
                        nBi++;
                    } else if (jColumn == -2) {
                        linearTerm = value;
                    } else {
                        printf("bad nonlinear term %s\n", temp);
                        abort();
                    }
                    ifFirst = false;
                }
                coinModelOriginal.setObjective(iColumn, linearTerm);
            }
        }
        int iRow;
        int saveNBi = nBi;
        for (iRow = 0; iRow < numberRows; iRow++) {
            CoinModelLink triple = coinModel_.firstInRow(iRow);
            while (triple.column() >= 0) {
                int iColumn = triple.column();
                const char *  el = coinModel_.getElementAsString(iRow, iColumn);
                if (strcmp("Numeric", el)) {
                    // check if value*x+-value*y....
                    assert (strlen(el) < 20000);
                    char temp[20000];
                    strcpy(temp, el);
                    char * pos = temp;
                    bool ifFirst = true;
                    double linearTerm = 0.0;
                    tryColumn[iColumn] = 1;
                    freeRow[iRow] = 1;
                    while (*pos) {
                        double value;
                        int jColumn = decodeBit(pos, pos, value, ifFirst, coinModel_);
                        // must be column unless first when may be linear term
                        if (jColumn >= 0) {
                            tryColumn[jColumn] = 1;
                            nBi++;
                        } else if (jColumn == -2) {
                            linearTerm = value;
                        } else {
                            printf("bad nonlinear term %s\n", temp);
                            abort();
                        }
                        ifFirst = false;
                    }
                    coinModelOriginal.setElement(iRow, iColumn, linearTerm);
                }
                triple = coinModel_.next(triple);
            }
        }
        if (!nBi)
            exit(1);
        bool quadraticObjective = false;
        int nInt = 0;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (coinModel_.isInteger(iColumn))
                nInt++;
        }
        printf("There are %d bilinear and %d integers\n", nBi, nInt);
        loadFromCoinModel(coinModelOriginal, true);
        CoinModel coinModel = coinModelOriginal;
        if (tightenBounds && numberColumns < 100) {
            // first fake bounds
            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                if (tryColumn[iColumn]) {
                    fakeBounds(this, iColumn, defaultBound_, &coinModel, &coinModel_);
                }
            }
            ClpSimplex tempModel(*modelPtr_);
            int nDelete = 0;
            for (iRow = 0; iRow < numberRows; iRow++) {
                if (freeRow[iRow])
                    freeRow[nDelete++] = iRow;
            }
            tempModel.deleteRows(nDelete, freeRow);
            tempModel.setOptimizationDirection(1.0);
            if (logLevel < 3) {
                tempModel.setLogLevel(0);
                tempModel.setSpecialOptions(32768);
            }
            double * objective = tempModel.objective();
            CoinZeroN(objective, numberColumns);
            // now up and down
            double * columnLower = modelPtr_->columnLower();
            double * columnUpper = modelPtr_->columnUpper();
            const double * solution = tempModel.primalColumnSolution();
            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                if (tryColumn[iColumn]) {
                    objective[iColumn] = 1.0;
                    tempModel.primal(1);
                    if (solution[iColumn] > columnLower[iColumn] + 1.0e-3) {
                        double value = solution[iColumn];
                        if (coinModel_.isInteger(iColumn))
                            value = ceil(value - 0.9e-3);
                        if (logLevel > 1)
                            printf("lower bound on %d changed from %g to %g\n", iColumn, columnLower[iColumn], value);
                        columnLower[iColumn] = value;
                        coinModel_.setColumnLower(iColumn, value);
                        coinModel.setColumnLower(iColumn, value);
                    }
                    objective[iColumn] = -1.0;
                    tempModel.primal(1);
                    if (solution[iColumn] < columnUpper[iColumn] - 1.0e-3) {
                        double value = solution[iColumn];
                        if (coinModel_.isInteger(iColumn))
                            value = floor(value + 0.9e-3);
                        if (logLevel > 1)
                            printf("upper bound on %d changed from %g to %g\n", iColumn, columnUpper[iColumn], value);
                        columnUpper[iColumn] = value;
                        coinModel_.setColumnUpper(iColumn, value);
                        coinModel.setColumnUpper(iColumn, value);
                    }
                    objective[iColumn] = 0.0;
                }
            }
        }
        delete [] freeRow;
        delete [] tryColumn;
        CoinBigIndex * startQuadratic = NULL;
        int * columnQuadratic = NULL;
        double * elementQuadratic = NULL;
        if ( saveNBi == nBi) {
            printf("all bilinearity in objective\n");
            specialOptions2_ |= 2;
            quadraticObjective = true;
            // save copy as quadratic model
            quadraticModel_ = new ClpSimplex(*modelPtr_);
            startQuadratic = new CoinBigIndex [numberColumns+1];
            columnQuadratic = new int [numberQuadratic];
            elementQuadratic = new double [numberQuadratic];
            numberQuadratic = 0;
        }
        //if (quadraticObjective||((specialOptions2_&8)!=0&&saveNBi)) {
        if (saveNBi) {
            // add in objective as constraint
            objectiveVariable_ = numberColumns;
            objectiveRow_ = coinModel.numberRows();
            coinModel.addColumn(0, NULL, NULL, -COIN_DBL_MAX, COIN_DBL_MAX, 1.0);
            int * column = new int[numberColumns+1];
            double * element = new double[numberColumns+1];
            double * objective = coinModel.objectiveArray();
            int n = 0;
            for (int i = 0; i < numberColumns; i++) {
                if (objective[i]) {
                    column[n] = i;
                    element[n++] = objective[i];
                    objective[i] = 0.0;
                }
            }
            column[n] = objectiveVariable_;
            element[n++] = -1.0;
            double offset = - coinModel.objectiveOffset();
            //assert (!offset); // get sign right if happens
            printf("***** offset %g\n", offset);
            coinModel.setObjectiveOffset(0.0);
            double lowerBound = -COIN_DBL_MAX;
            coinModel.addRow(n, column, element, lowerBound, offset);
            delete [] column;
            delete [] element;
        }
        OsiObject ** objects = new OsiObject * [nBi+nInt];
        char * marked = new char [numberColumns];
        memset(marked, 0, numberColumns);
        // statistics I-I I-x x-x
        int stats[3] = {0, 0, 0};
        double * sort = new double [nBi];
        nBi = nInt;
        const OsiObject ** justBi = const_cast<const OsiObject **> (objects + nInt);
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (quadraticObjective)
                startQuadratic[iColumn] = numberQuadratic;
            // See if quadratic objective
            const char * expr = coinModel_.getColumnObjectiveAsString(iColumn);
            if (strcmp(expr, "Numeric")) {
                // need bounds
                fakeBounds(this, iColumn, defaultBound_, &coinModel, &coinModel_);
                // value*x*y
                char temp[20000];
                strcpy(temp, expr);
                char * pos = temp;
                bool ifFirst = true;
                while (*pos) {
                    double value;
                    int jColumn = decodeBit(pos, pos, value, ifFirst, coinModel_);
                    // must be column unless first when may be linear term
                    if (jColumn >= 0) {
                        if (quadraticObjective) {
                            columnQuadratic[numberQuadratic] = jColumn;
                            if (jColumn == iColumn)
                                elementQuadratic[numberQuadratic++] = 2.0 * value; // convention
                            else
                                elementQuadratic[numberQuadratic++] = 1.0 * value; // convention
                        }
                        // need bounds
                        fakeBounds(this, jColumn, defaultBound_, &coinModel, &coinModel_);
                        double meshI = coinModel_.isInteger(iColumn) ? 1.0 : 0.0;
                        if (meshI)
                            marked[iColumn] = 1;
                        double meshJ = coinModel_.isInteger(jColumn) ? 1.0 : 0.0;
                        if (meshJ)
                            marked[jColumn] = 1;
                        // stats etc
                        if (meshI) {
                            if (meshJ)
                                stats[0]++;
                            else
                                stats[1]++;
                        } else {
                            if (meshJ)
                                stats[1]++;
                            else
                                stats[2]++;
                        }
                        if (iColumn <= jColumn)
                            sort[nBi-nInt] = iColumn + numberColumns * jColumn;
                        else
                            sort[nBi-nInt] = jColumn + numberColumns * iColumn;
                        if (!meshJ && !meshI) {
                            meshI = defaultMeshSize_;
                            meshJ = 0.0;
                        }
                        OsiBiLinear * newObj = new OsiBiLinear(&coinModel, iColumn, jColumn, objectiveRow_, value, meshI, meshJ,
                                                               nBi - nInt, justBi);
                        newObj->setPriority(biLinearPriority_);
                        objects[nBi++] = newObj;
                    } else if (jColumn == -2) {
                    } else {
                        printf("bad nonlinear term %s\n", temp);
                        abort();
                    }
                    ifFirst = false;
                }
            }
        }
        // stats
        printf("There were %d I-I, %d I-x and %d x-x bilinear in objective\n",
               stats[0], stats[1], stats[2]);
        if (quadraticObjective) {
            startQuadratic[numberColumns] = numberQuadratic;
            quadraticModel_->loadQuadraticObjective(numberColumns, startQuadratic,
                                                    columnQuadratic, elementQuadratic);
            delete [] startQuadratic;
            delete [] columnQuadratic;
            delete [] elementQuadratic;
        }
        for (iRow = 0; iRow < numberRows; iRow++) {
            CoinModelLink triple = coinModel_.firstInRow(iRow);
            while (triple.column() >= 0) {
                int iColumn = triple.column();
                const char *  el = coinModel_.getElementAsString(iRow, iColumn);
                if (strcmp("Numeric", el)) {
                    // need bounds
                    fakeBounds(this, iColumn, defaultBound_, &coinModel, &coinModel_);
                    // value*x*y
                    char temp[20000];
                    strcpy(temp, el);
                    char * pos = temp;
                    bool ifFirst = true;
                    while (*pos) {
                        double value;
                        int jColumn = decodeBit(pos, pos, value, ifFirst, coinModel_);
                        // must be column unless first when may be linear term
                        if (jColumn >= 0) {
                            // need bounds
                            fakeBounds(this, jColumn, defaultBound_, &coinModel, &coinModel_);
                            double meshI = coinModel_.isInteger(iColumn) ? 1.0 : 0.0;
                            if (meshI)
                                marked[iColumn] = 1;
                            double meshJ = coinModel_.isInteger(jColumn) ? 1.0 : 0.0;
                            if (meshJ)
                                marked[jColumn] = 1;
                            // stats etc
                            if (meshI) {
                                if (meshJ)
                                    stats[0]++;
                                else
                                    stats[1]++;
                            } else {
                                if (meshJ)
                                    stats[1]++;
                                else
                                    stats[2]++;
                            }
                            if (iColumn <= jColumn)
                                sort[nBi-nInt] = iColumn + numberColumns * jColumn;
                            else
                                sort[nBi-nInt] = jColumn + numberColumns * iColumn;
                            if (!meshJ && !meshI) {
                                meshI = defaultMeshSize_;
                                meshJ = 0.0;
                            }
                            OsiBiLinear * newObj = new OsiBiLinear(&coinModel, iColumn, jColumn, iRow, value, meshI, meshJ,
                                                                   nBi - nInt, justBi);
                            newObj->setPriority(biLinearPriority_);
                            objects[nBi++] = newObj;
                        } else if (jColumn == -2) {
                        } else {
                            printf("bad nonlinear term %s\n", temp);
                            abort();
                        }
                        ifFirst = false;
                    }
                }
                triple = coinModel_.next(triple);
            }
        }
        {
            // stats
            std::sort(sort, sort + nBi - nInt);
            int nDiff = 0;
            double last = -1.0;
            for (int i = 0; i < nBi - nInt; i++) {
                if (sort[i] != last)
                    nDiff++;
                last = sort[i];
            }
            delete [] sort;
            printf("There were %d I-I, %d I-x and %d x-x bilinear in total of which %d were duplicates\n",
                   stats[0], stats[1], stats[2], nBi - nInt - nDiff);
        }
        // reload with all bilinear stuff
        loadFromCoinModel(coinModel, true);
        //exit(77);
        nInt = 0;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (coinModel_.isInteger(iColumn)) {
                objects[nInt] = new OsiSimpleInteger(this, iColumn);
                if (marked[iColumn])
                    objects[nInt]->setPriority(integerPriority_);
                else
                    objects[nInt]->setPriority(integerPriority_);
                nInt++;
            }
        }
        nInt = nBi;
        delete [] marked;
        if (numberErrors) {
            // errors
            gutsOfDestructor();
            numberVariables_ = -1;
        } else {
            addObjects(nInt, objects);
            int i;
            for (i = 0; i < nInt; i++)
                delete objects[i];
            delete [] objects;
            // Now do dummy bound stuff
            matrix_ = new CoinPackedMatrix(*getMatrixByCol());
            info_ = new OsiLinkedBound [numberVariables_];
            for ( i = 0; i < numberVariables_; i++) {
                info_[i] = OsiLinkedBound(this, which[i], 0, NULL, NULL, NULL);
            }
            // Do row copy but just part
            int numberRows2 = objectiveRow_ >= 0 ? numberRows + 1 : numberRows;
            int * whichRows = new int [numberRows2];
            int * whichColumns = new int [numberColumns];
            CoinIotaN(whichRows, numberRows2, 0);
            CoinIotaN(whichColumns, numberColumns, 0);
            originalRowCopy_ = new CoinPackedMatrix(*getMatrixByRow(),
                                                    numberRows2, whichRows,
                                                    numberColumns, whichColumns);
            delete [] whichColumns;
            numberNonLinearRows_ = 0;
            CoinZeroN(whichRows, numberRows2);
            for ( i = 0; i < numberObjects_; i++) {
                OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[i]);
                if (obj) {
                    int xyRow = obj->xyRow();
                    assert (xyRow >= 0 && xyRow < numberRows2); // even if obj we should move
                    whichRows[xyRow]++;
                }
            }
            int * pos = new int [numberRows2];
            int n = 0;
            for (i = 0; i < numberRows2; i++) {
                if (whichRows[i]) {
                    pos[numberNonLinearRows_] = n;
                    n += whichRows[i];
                    whichRows[i] = numberNonLinearRows_;
                    numberNonLinearRows_++;
                } else {
                    whichRows[i] = -1;
                }
            }
            startNonLinear_ = new int [numberNonLinearRows_+1];
            memcpy(startNonLinear_, pos, numberNonLinearRows_*sizeof(int));
            startNonLinear_[numberNonLinearRows_] = n;
            rowNonLinear_ = new int [numberNonLinearRows_];
            convex_ = new int [numberNonLinearRows_];
            // do row numbers now
            numberNonLinearRows_ = 0;
            for (i = 0; i < numberRows2; i++) {
                if (whichRows[i] >= 0) {
                    rowNonLinear_[numberNonLinearRows_++] = i;
                }
            }
            whichNonLinear_ = new int [n];
            for ( i = 0; i < numberObjects_; i++) {
                OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[i]);
                if (obj) {
                    int xyRow = obj->xyRow();
                    int k = whichRows[xyRow];
                    int put = pos[k];
                    pos[k]++;
                    whichNonLinear_[put] = i;
                }
            }
            delete [] pos;
            delete [] whichRows;
            analyzeObjects();
        }
    }
    // See if there are any quadratic bounds
    int nQ = 0;
    const CoinPackedMatrix * rowCopy = getMatrixByRow();
    //const double * element = rowCopy->getElements();
    //const int * column = rowCopy->getIndices();
    //const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
    const int * rowLength = rowCopy->getVectorLengths();
    const double * rowLower = getRowLower();
    const double * rowUpper = getRowUpper();
    for (int iObject = 0; iObject < numberObjects_; iObject++) {
        OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[iObject]);
        if (obj) {
            int xyRow = obj->xyRow();
            if (rowLength[xyRow] == 4 && false) {
                // we have simple bound
                nQ++;
                double coefficient = obj->coefficient();
                double lo = rowLower[xyRow];
                double up = rowUpper[xyRow];
                if (coefficient != 1.0) {
                    printf("*** double check code here\n");
                    if (coefficient < 0.0) {
                        double temp = lo;
                        lo = - up;
                        up = - temp;
                        coefficient = - coefficient;
                    }
                    if (lo > -1.0e20)
                        lo /= coefficient;
                    if (up < 1.0e20)
                        up /= coefficient;
                    setRowLower(xyRow, lo);
                    setRowUpper(xyRow, up);
                    // we also need to change elements in matrix_
                }
                int type = 0;
                if (lo == up) {
                    // good news
                    type = 3;
                    coefficient = lo;
                } else if (lo < -1.0e20) {
                    assert (up < 1.0e20);
                    coefficient = up;
                    type = 1;
                    // can we make equality?
                } else if (up > 1.0e20) {
                    coefficient = lo;
                    type = 2;
                    // can we make equality?
                } else {
                    // we would need extra code
                    abort();
                }
                obj->setBoundType(type);
                obj->setCoefficient(coefficient);
                // can do better if integer?
                assert (!isInteger(obj->xColumn()));
                assert (!isInteger(obj->yColumn()));
            }
        }
    }
    delete [] which;
    if ((specialOptions2_&16) != 0)
        addTighterConstraints();
}
// Add reformulated bilinear constraints
void
OsiSolverLink::addTighterConstraints()
{
    // This is first attempt - for now get working on trimloss
    int numberW = 0;
    int * xW = new int[numberObjects_];
    int * yW = new int[numberObjects_];
    // Points to firstlambda
    int * wW = new int[numberObjects_];
    // Coefficient
    double * alphaW = new double[numberObjects_];
    // Objects
    OsiBiLinear ** objW = new OsiBiLinear * [numberObjects_];
    int numberColumns = getNumCols();
    int firstLambda = numberColumns;
    // set up list (better to rethink and do properly as column ordered)
    int * list = new int[numberColumns];
    memset(list, 0, numberColumns*sizeof(int));
    int i;
    for ( i = 0; i < numberObjects_; i++) {
        OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[i]);
        if (obj) {
            //obj->setBranchingStrategy(4); // ***** temp
            objW[numberW] = obj;
            xW[numberW] = obj->xColumn();
            yW[numberW] = obj->yColumn();
            list[xW[numberW]] = 1;
            list[yW[numberW]] = 1;
            wW[numberW] = obj->firstLambda();
            firstLambda = CoinMin(firstLambda, obj->firstLambda());
            alphaW[numberW] = obj->coefficient();
            //assert (alphaW[numberW]==1.0); // fix when occurs
            numberW++;
        }
    }
    int nList = 0;
    for (i = 0; i < numberColumns; i++) {
        if (list[i])
            list[nList++] = i;
    }
    // set up mark array
    char * mark = new char [firstLambda*firstLambda];
    memset(mark, 0, firstLambda*firstLambda);
    for (i = 0; i < numberW; i++) {
        int x = xW[i];
        int y = yW[i];
        mark[x*firstLambda+y] = 1;
        mark[y*firstLambda+x] = 1;
    }
    int numberRows2 = originalRowCopy_->getNumRows();
    int * addColumn = new int [numberColumns];
    double * addElement = new double [numberColumns];
    int * addW = new int [numberColumns];
    assert (objectiveRow_ < 0); // fix when occurs
    for (int iRow = 0; iRow < numberRows2; iRow++) {
        for (int iList = 0; iList < nList; iList++) {
            int kColumn = list[iList];
#ifndef NDEBUG
            const double * columnLower = getColLower();
#endif
            //const double * columnUpper = getColUpper();
            const double * rowLower = getRowLower();
            const double * rowUpper = getRowUpper();
            const CoinPackedMatrix * rowCopy = getMatrixByRow();
            const double * element = rowCopy->getElements();
            const int * column = rowCopy->getIndices();
            const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
            const int * rowLength = rowCopy->getVectorLengths();
            CoinBigIndex j;
            int numberElements = rowLength[iRow];
            int n = 0;
            for (j = rowStart[iRow]; j < rowStart[iRow] + numberElements; j++) {
                int iColumn = column[j];
                if (iColumn >= firstLambda) {
                    // no good
                    n = -1;
                    break;
                }
                if (mark[iColumn*firstLambda+kColumn])
                    n++;
            }
            if (n == numberElements) {
                printf("can add row %d\n", iRow);
                assert (columnLower[kColumn] >= 0); // might be able to fix
                n = 0;
                for (j = rowStart[iRow]; j < rowStart[iRow] + numberElements; j++) {
                    int xColumn = kColumn;
                    int yColumn = column[j];
                    int k;
                    for (k = 0; k < numberW; k++) {
                        if ((xW[k] == yColumn && yW[k] == xColumn) ||
                                (yW[k] == yColumn && xW[k] == xColumn))
                            break;
                    }
                    assert (k < numberW);
                    if (xW[k] != xColumn) {
                        int temp = xColumn;
                        xColumn = yColumn;
                        yColumn = temp;
                    }
                    addW[n/4] = k;
                    int start = wW[k];
                    double value = element[j];
                    for (int kk = 0; kk < 4; kk++) {
                        // Dummy value
                        addElement[n] = value;
                        addColumn[n++] = start + kk;
                    }
                }
                addColumn[n++] = kColumn;
                double lo = rowLower[iRow];
                double up = rowUpper[iRow];
                if (lo > -1.0e20) {
                    // and tell object
                    for (j = 0; j < n - 1; j += 4) {
                        int iObject = addW[j/4];
                        objW[iObject]->addExtraRow(matrix_->getNumRows(), addElement[j]);
                    }
                    addElement[n-1] = -lo;
                    if (lo == up)
                        addRow(n, addColumn, addElement, 0.0, 0.0);
                    else
                        addRow(n, addColumn, addElement, 0.0, COIN_DBL_MAX);
                    matrix_->appendRow(n, addColumn, addElement);
                }
                if (up<1.0e20 && up>lo) {
                    // and tell object
                    for (j = 0; j < n - 1; j += 4) {
                        int iObject = addW[j/4];
                        objW[iObject]->addExtraRow(matrix_->getNumRows(), addElement[j]);
                    }
                    addElement[n-1] = -up;
                    addRow(n, addColumn, addElement, -COIN_DBL_MAX, 0.0);
                    matrix_->appendRow(n, addColumn, addElement);
                }
            }
        }
    }
#ifdef JJF_ZERO
    // possibly do bounds
    for (int iColumn = 0; iColumn < firstLambda; iColumn++) {
        for (int iList = 0; iList < nList; iList++) {
            int kColumn = list[iList];
            const double * columnLower = getColLower();
            const double * columnUpper = getColUpper();
            if (mark[iColumn*firstLambda+kColumn]) {
                printf("can add column %d\n", iColumn);
                assert (columnLower[kColumn] >= 0); // might be able to fix
                int xColumn = kColumn;
                int yColumn = iColumn;
                int k;
                for (k = 0; k < numberW; k++) {
                    if ((xW[k] == yColumn && yW[k] == xColumn) ||
                            (yW[k] == yColumn && xW[k] == xColumn))
                        break;
                }
                assert (k < numberW);
                if (xW[k] != xColumn) {
                    int temp = xColumn;
                    xColumn = yColumn;
                    yColumn = temp;
                }
                int start = wW[k];
                int n = 0;
                for (int kk = 0; kk < 4; kk++) {
                    // Dummy value
                    addElement[n] = 1.0e-19;
                    addColumn[n++] = start + kk;
                }
                // Tell object about this
                objW[k]->addExtraRow(matrix_->getNumRows(), 1.0);
                addColumn[n++] = kColumn;
                double lo = columnLower[iColumn];
                double up = columnUpper[iColumn];
                if (lo > -1.0e20) {
                    addElement[n-1] = -lo;
                    if (lo == up)
                        addRow(n, addColumn, addElement, 0.0, 0.0);
                    else
                        addRow(n, addColumn, addElement, 0.0, COIN_DBL_MAX);
                    matrix_->appendRow(n, addColumn, addElement);
                }
                if (up<1.0e20 && up>lo) {
                    addElement[n-1] = -up;
                    addRow(n, addColumn, addElement, -COIN_DBL_MAX, 0.0);
                    matrix_->appendRow(n, addColumn, addElement);
                }
            }
        }
    }
#endif
    delete [] xW;
    delete [] yW;
    delete [] wW;
    delete [] alphaW;
    delete [] addColumn;
    delete [] addElement;
    delete [] addW;
    delete [] mark;
    delete [] list;
    delete [] objW;
}
// Set all biLinear priorities on x-x variables
void
OsiSolverLink::setBiLinearPriorities(int value, double meshSize)
{
    OsiObject ** newObject = new OsiObject * [numberObjects_];
    int numberOdd = 0;
    int i;
    for ( i = 0; i < numberObjects_; i++) {
        OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[i]);
        if (obj) {
            if (obj->xMeshSize() < 1.0 && obj->yMeshSize() < 1.0) {
                double oldSatisfied = CoinMax(obj->xSatisfied(),
                                              obj->ySatisfied());
                OsiBiLinear * objNew = new OsiBiLinear(*obj);
                newObject[numberOdd++] = objNew;
                objNew->setXSatisfied(0.5*meshSize);
                obj->setXOtherSatisfied(0.5*meshSize);
                objNew->setXOtherSatisfied(oldSatisfied);
                objNew->setXMeshSize(meshSize);
                objNew->setYSatisfied(0.5*meshSize);
                obj->setYOtherSatisfied(0.5*meshSize);
                objNew->setYOtherSatisfied(oldSatisfied);
                objNew->setYMeshSize(meshSize);
                objNew->setXYSatisfied(0.25*meshSize);
                objNew->setPriority(value);
                objNew->setBranchingStrategy(8);
            }
        }
    }
    addObjects(numberOdd, newObject);
    for (i = 0; i < numberOdd; i++)
        delete newObject[i];
    delete [] newObject;
}
/* Set options and priority on all or some biLinear variables
   1 - on I-I
   2 - on I-x
   4 - on x-x
      or combinations.
      -1 means leave (for priority value and strategy value)
*/
void
OsiSolverLink::setBranchingStrategyOnVariables(int strategyValue, int priorityValue,
        int mode)
{
    int i;
    for ( i = 0; i < numberObjects_; i++) {
        OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[i]);
        if (obj) {
            bool change = false;
            if (obj->xMeshSize() < 1.0 && obj->yMeshSize() < 1.0 && (mode&4) != 0)
                change = true;
            else if (((obj->xMeshSize() == 1.0 && obj->yMeshSize() < 1.0) ||
                      (obj->xMeshSize() < 1.0 && obj->yMeshSize() == 1.0)) && (mode&2) != 0)
                change = true;
            else if (obj->xMeshSize() == 1.0 && obj->yMeshSize() == 1.0 && (mode&1) != 0)
                change = true;
            else if (obj->xMeshSize() > 1.0 || obj->yMeshSize() > 1.0)
                abort();
            if (change) {
                if (strategyValue >= 0)
                    obj->setBranchingStrategy(strategyValue);
                if (priorityValue >= 0)
                    obj->setPriority(priorityValue);
            }
        }
    }
}

// Say convex (should work it out)
void
OsiSolverLink::sayConvex(bool convex)
{
    specialOptions2_ |= 4;
    if (convex_) {
        for (int iNon = 0; iNon < numberNonLinearRows_; iNon++) {
            convex_[iNon] = convex ? 1 : -1;
        }
    }
}
// Set all mesh sizes on x-x variables
void
OsiSolverLink::setMeshSizes(double value)
{
    int i;
    for ( i = 0; i < numberObjects_; i++) {
        OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[i]);
        if (obj) {
            if (obj->xMeshSize() < 1.0 && obj->yMeshSize() < 1.0) {
#ifdef JJF_ZERO
                numberContinuous++;
                int xColumn = obj->xColumn();
                double gapX = upper[xColumn] - lower[xColumn];
                int yColumn = obj->yColumn();
                double gapY = upper[yColumn] - lower[yColumn];
                gap = CoinMax(gap, CoinMax(gapX, gapY));
#endif
                obj->setMeshSizes(this, value, value);
            }
        }
    }
}
/* Solves nonlinear problem from CoinModel using SLP - may be used as crash
   for other algorithms when number of iterations small.
   Also exits if all problematical variables are changing
   less than deltaTolerance
   Returns solution array
*/
double *
OsiSolverLink::nonlinearSLP(int numberPasses, double deltaTolerance)
{
    if (!coinModel_.numberRows()) {
        printf("Model not set up or nonlinear arrays not created!\n");
        return NULL;
    }
    // first check and set up arrays
    int numberColumns = coinModel_.numberColumns();
    int numberRows = coinModel_.numberRows();
    char * markNonlinear = new char [numberColumns+numberRows];
    CoinZeroN(markNonlinear, numberColumns + numberRows);
    // List of nonlinear entries
    int * listNonLinearColumn = new int[numberColumns];
    // List of nonlinear constraints
    int * whichRow = new int [numberRows];
    CoinZeroN(whichRow, numberRows);
    int numberNonLinearColumns = 0;
    int iColumn;
    CoinModel coinModel = coinModel_;
    //const CoinModelHash * stringArray = coinModel.stringArray();
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        CoinModelLink triple = coinModel.firstInColumn(iColumn);
        bool linear = true;
        int n = 0;
        // See if nonlinear objective
        const char * expr = coinModel.getColumnObjectiveAsString(iColumn);
        if (strcmp(expr, "Numeric")) {
            linear = false;
            // try and see which columns
            assert (strlen(expr) < 20000);
            char temp[20000];
            strcpy(temp, expr);
            char * pos = temp;
            bool ifFirst = true;
            double linearTerm = 0.0;
            while (*pos) {
                double value;
                int jColumn = decodeBit(pos, pos, value, ifFirst, coinModel_);
                // must be column unless first when may be linear term
                if (jColumn >= 0) {
                    markNonlinear[jColumn] = 1;
                } else if (jColumn == -2) {
                    linearTerm = value;
                } else {
                    printf("bad nonlinear term %s\n", temp);
                    abort();
                }
                ifFirst = false;
            }
        }
        while (triple.row() >= 0) {
            int iRow = triple.row();
            const char * expr = coinModel.getElementAsString(iRow, iColumn);
            if (strcmp(expr, "Numeric")) {
                linear = false;
                whichRow[iRow]++;
                // try and see which columns
                assert (strlen(expr) < 20000);
                char temp[20000];
                strcpy(temp, expr);
                char * pos = temp;
                bool ifFirst = true;
                double linearTerm = 0.0;
                while (*pos) {
                    double value;
                    int jColumn = decodeBit(pos, pos, value, ifFirst, coinModel_);
                    // must be column unless first when may be linear term
                    if (jColumn >= 0) {
                        markNonlinear[jColumn] = 1;
                    } else if (jColumn == -2) {
                        linearTerm = value;
                    } else {
                        printf("bad nonlinear term %s\n", temp);
                        abort();
                    }
                    ifFirst = false;
                }
            }
            triple = coinModel.next(triple);
            n++;
        }
        if (!linear) {
            markNonlinear[iColumn] = 1;
        }
    }
    //int xxxx[]={3,2,0,4,3,0};
    //double initialSolution[6];
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (markNonlinear[iColumn]) {
            // put in something
            double lower = coinModel.columnLower(iColumn);
            double upper = CoinMin(coinModel.columnUpper(iColumn), lower + 1000.0);
            coinModel.associateElement(coinModel.columnName(iColumn), 0.5*(lower + upper));
            //coinModel.associateElement(coinModel.columnName(iColumn),xxxx[iColumn]);
            listNonLinearColumn[numberNonLinearColumns++] = iColumn;
            //initialSolution[iColumn]=xxxx[iColumn];
        }
    }
    // if nothing just solve
    if (!numberNonLinearColumns) {
        delete [] listNonLinearColumn;
        delete [] whichRow;
        delete [] markNonlinear;
        ClpSimplex tempModel;
        tempModel.loadProblem(coinModel, true);
        tempModel.initialSolve();
        double * solution = CoinCopyOfArray(tempModel.getColSolution(), numberColumns);
        return solution;
    }
    // Create artificials
    ClpSimplex tempModel;
    tempModel.loadProblem(coinModel, true);
    const double * rowLower = tempModel.rowLower();
    const double * rowUpper = tempModel.rowUpper();
    bool takeAll = false;
    int iRow;
    int numberArtificials = 0;
    for (iRow = 0; iRow < numberRows; iRow++) {
        if (whichRow[iRow] || takeAll) {
            if (rowLower[iRow] > -1.0e30)
                numberArtificials++;
            if (rowUpper[iRow] < 1.0e30)
                numberArtificials++;
        }
    }
    CoinBigIndex * startArtificial = new CoinBigIndex [numberArtificials+1];
    int * rowArtificial = new int [numberArtificials];
    double * elementArtificial = new double [numberArtificials];
    double * objectiveArtificial = new double [numberArtificials];
    numberArtificials = 0;
    startArtificial[0] = 0;
    double artificialCost = 1.0e9;
    for (iRow = 0; iRow < numberRows; iRow++) {
        if (whichRow[iRow] || takeAll) {
            if (rowLower[iRow] > -1.0e30) {
                rowArtificial[numberArtificials] = iRow;
                elementArtificial[numberArtificials] = 1.0;
                objectiveArtificial[numberArtificials] = artificialCost;
                numberArtificials++;
                startArtificial[numberArtificials] = numberArtificials;
            }
            if (rowUpper[iRow] < 1.0e30) {
                rowArtificial[numberArtificials] = iRow;
                elementArtificial[numberArtificials] = -1.0;
                objectiveArtificial[numberArtificials] = artificialCost;
                numberArtificials++;
                startArtificial[numberArtificials] = numberArtificials;
            }
        }
    }
    // Get first solution
    int numberColumnsSmall = numberColumns;
    ClpSimplex model;
    model.loadProblem(coinModel, true);
    model.addColumns(numberArtificials, NULL, NULL, objectiveArtificial,
                     startArtificial, rowArtificial, elementArtificial);
    double * columnLower = model.columnLower();
    double * columnUpper = model.columnUpper();
    double * trueLower = new double[numberNonLinearColumns];
    double * trueUpper = new double[numberNonLinearColumns];
    int jNon;
    for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
        iColumn = listNonLinearColumn[jNon];
        trueLower[jNon] = columnLower[iColumn];
        trueUpper[jNon] = columnUpper[iColumn];
        //columnLower[iColumn]=initialSolution[iColumn];
        //columnUpper[iColumn]=initialSolution[iColumn];
    }
    model.initialSolve();
    //model.writeMps("bad.mps");
    // redo number of columns
    numberColumns = model.numberColumns();
    int * last[3];
    double * solution = model.primalColumnSolution();

    double * trust = new double[numberNonLinearColumns];
    for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
        iColumn = listNonLinearColumn[jNon];
        trust[jNon] = 0.5;
        if (solution[iColumn] < trueLower[jNon])
            solution[iColumn] = trueLower[jNon];
        else if (solution[iColumn] > trueUpper[jNon])
            solution[iColumn] = trueUpper[jNon];
    }
    int iPass;
    double lastObjective = 1.0e31;
    double * saveSolution = new double [numberColumns];
    double * saveRowSolution = new double [numberRows];
    memset(saveRowSolution, 0, numberRows*sizeof(double));
    double * savePi = new double [numberRows];
    double * safeSolution = new double [numberColumns];
    unsigned char * saveStatus = new unsigned char[numberRows+numberColumns];
    double targetDrop = 1.0e31;
    //double objectiveOffset;
    //model.getDblParam(ClpObjOffset,objectiveOffset);
    // 1 bound up, 2 up, -1 bound down, -2 down, 0 no change
    for (iPass = 0; iPass < 3; iPass++) {
        last[iPass] = new int[numberNonLinearColumns];
        for (jNon = 0; jNon < numberNonLinearColumns; jNon++)
            last[iPass][jNon] = 0;
    }
    // goodMove +1 yes, 0 no, -1 last was bad - just halve gaps, -2 do nothing
    int goodMove = -2;
    char * statusCheck = new char[numberColumns];
    double * changeRegion = new double [numberColumns];
    int logLevel = 63;
    double dualTolerance = model.dualTolerance();
    double primalTolerance = model.primalTolerance();
    int lastGoodMove = 1;
    for (iPass = 0; iPass < numberPasses; iPass++) {
        lastGoodMove = goodMove;
        columnLower = model.columnLower();
        columnUpper = model.columnUpper();
        solution = model.primalColumnSolution();
        double * rowActivity = model.primalRowSolution();
        // redo objective
        ClpSimplex tempModel;
        // load new values
        for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
            iColumn = listNonLinearColumn[jNon];
            coinModel.associateElement(coinModel.columnName(iColumn), solution[iColumn]);
        }
        tempModel.loadProblem(coinModel);
        double objectiveOffset;
        tempModel.getDblParam(ClpObjOffset, objectiveOffset);
        double objValue = -objectiveOffset;
        const double * objective = tempModel.objective();
        for (iColumn = 0; iColumn < numberColumnsSmall; iColumn++)
            objValue += solution[iColumn] * objective[iColumn];
        double * rowActivity2 = tempModel.primalRowSolution();
        const double * rowLower2 = tempModel.rowLower();
        const double * rowUpper2 = tempModel.rowUpper();
        memset(rowActivity2, 0, numberRows*sizeof(double));
        tempModel.times(1.0, solution, rowActivity2);
        for (iRow = 0; iRow < numberRows; iRow++) {
            if (rowActivity2[iRow] < rowLower2[iRow] - primalTolerance)
                objValue += (rowLower2[iRow] - rowActivity2[iRow] - primalTolerance) * artificialCost;
            else if (rowActivity2[iRow] > rowUpper2[iRow] + primalTolerance)
                objValue -= (rowUpper2[iRow] - rowActivity2[iRow] + primalTolerance) * artificialCost;
        }
        double theta = -1.0;
        double maxTheta = COIN_DBL_MAX;
        if (objValue <= lastObjective + 1.0e-15*fabs(lastObjective) || !iPass)
            goodMove = 1;
        else
            goodMove = -1;
        //maxTheta=1.0;
        if (iPass) {
            int jNon = 0;
            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                changeRegion[iColumn] = solution[iColumn] - saveSolution[iColumn];
                double alpha = changeRegion[iColumn];
                double oldValue = saveSolution[iColumn];
                if (markNonlinear[iColumn] == 0) {
                    // linear
                    if (alpha < -1.0e-15) {
                        // variable going towards lower bound
                        double bound = columnLower[iColumn];
                        oldValue -= bound;
                        if (oldValue + maxTheta*alpha < 0.0) {
                            maxTheta = CoinMax(0.0, oldValue / (-alpha));
                        }
                    } else if (alpha > 1.0e-15) {
                        // variable going towards upper bound
                        double bound = columnUpper[iColumn];
                        oldValue = bound - oldValue;
                        if (oldValue - maxTheta*alpha < 0.0) {
                            maxTheta = CoinMax(0.0, oldValue / alpha);
                        }
                    }
                } else {
                    // nonlinear
                    if (alpha < -1.0e-15) {
                        // variable going towards lower bound
                        double bound = trueLower[jNon];
                        oldValue -= bound;
                        if (oldValue + maxTheta*alpha < 0.0) {
                            maxTheta = CoinMax(0.0, oldValue / (-alpha));
                        }
                    } else if (alpha > 1.0e-15) {
                        // variable going towards upper bound
                        double bound = trueUpper[jNon];
                        oldValue = bound - oldValue;
                        if (oldValue - maxTheta*alpha < 0.0) {
                            maxTheta = CoinMax(0.0, oldValue / alpha);
                        }
                    }
                    jNon++;
                }
            }
            // make sure both accurate
            memset(rowActivity, 0, numberRows*sizeof(double));
            model.times(1.0, solution, rowActivity);
            memset(saveRowSolution, 0, numberRows*sizeof(double));
            model.times(1.0, saveSolution, saveRowSolution);
            for (int iRow = 0; iRow < numberRows; iRow++) {
                double alpha = rowActivity[iRow] - saveRowSolution[iRow];
                double oldValue = saveRowSolution[iRow];
                if (alpha < -1.0e-15) {
                    // variable going towards lower bound
                    double bound = rowLower[iRow];
                    oldValue -= bound;
                    if (oldValue + maxTheta*alpha < 0.0) {
                        maxTheta = CoinMax(0.0, oldValue / (-alpha));
                    }
                } else if (alpha > 1.0e-15) {
                    // variable going towards upper bound
                    double bound = rowUpper[iRow];
                    oldValue = bound - oldValue;
                    if (oldValue - maxTheta*alpha < 0.0) {
                        maxTheta = CoinMax(0.0, oldValue / alpha);
                    }
                }
            }
        } else {
            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                changeRegion[iColumn] = 0.0;
                saveSolution[iColumn] = solution[iColumn];
            }
            memcpy(saveRowSolution, rowActivity, numberRows*sizeof(double));
        }
        if (goodMove >= 0) {
            //theta = CoinMin(theta2,maxTheta);
            theta = maxTheta;
            if (theta > 0.0 && theta <= 1.0) {
                // update solution
                double lambda = 1.0 - theta;
                for (iColumn = 0; iColumn < numberColumns; iColumn++)
                    solution[iColumn] = lambda * saveSolution[iColumn]
                                        + theta * solution[iColumn];
                memset(rowActivity, 0, numberRows*sizeof(double));
                model.times(1.0, solution, rowActivity);
                if (lambda > 0.999) {
                    memcpy(model.dualRowSolution(), savePi, numberRows*sizeof(double));
                    memcpy(model.statusArray(), saveStatus, numberRows + numberColumns);
                }
                // redo rowActivity
                memset(rowActivity, 0, numberRows*sizeof(double));
                model.times(1.0, solution, rowActivity);
            }
        }
        // load new values
        for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
            iColumn = listNonLinearColumn[jNon];
            coinModel.associateElement(coinModel.columnName(iColumn), solution[iColumn]);
        }
        double * sol2 = CoinCopyOfArray(model.primalColumnSolution(), numberColumns);
        unsigned char * status2 = CoinCopyOfArray(model.statusArray(), numberColumns);
        model.loadProblem(coinModel);
        model.addColumns(numberArtificials, NULL, NULL, objectiveArtificial,
                         startArtificial, rowArtificial, elementArtificial);
        memcpy(model.primalColumnSolution(), sol2, numberColumns*sizeof(double));
        memcpy(model.statusArray(), status2, numberColumns);
        delete [] sol2;
        delete [] status2;
        columnLower = model.columnLower();
        columnUpper = model.columnUpper();
        solution = model.primalColumnSolution();
        rowActivity = model.primalRowSolution();
        int * temp = last[2];
        last[2] = last[1];
        last[1] = last[0];
        last[0] = temp;
        for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
            iColumn = listNonLinearColumn[jNon];
            double change = solution[iColumn] - saveSolution[iColumn];
            if (change < -1.0e-5) {
                if (fabs(change + trust[jNon]) < 1.0e-5)
                    temp[jNon] = -1;
                else
                    temp[jNon] = -2;
            } else if (change > 1.0e-5) {
                if (fabs(change - trust[jNon]) < 1.0e-5)
                    temp[jNon] = 1;
                else
                    temp[jNon] = 2;
            } else {
                temp[jNon] = 0;
            }
        }
        // goodMove +1 yes, 0 no, -1 last was bad - just halve gaps, -2 do nothing
        double maxDelta = 0.0;
        if (goodMove >= 0) {
            if (objValue <= lastObjective + 1.0e-15*fabs(lastObjective))
                goodMove = 1;
            else
                goodMove = 0;
        } else {
            maxDelta = 1.0e10;
        }
        double maxGap = 0.0;
        int numberSmaller = 0;
        int numberSmaller2 = 0;
        int numberLarger = 0;
        for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
            iColumn = listNonLinearColumn[jNon];
            maxDelta = CoinMax(maxDelta,
                               fabs(solution[iColumn] - saveSolution[iColumn]));
            if (goodMove > 0) {
                if (last[0][jNon]*last[1][jNon] < 0) {
                    // halve
                    trust[jNon] *= 0.5;
                    numberSmaller2++;
                } else {
                    if (last[0][jNon] == last[1][jNon] &&
                            last[0][jNon] == last[2][jNon])
                        trust[jNon] = CoinMin(1.5 * trust[jNon], 1.0e6);
                    numberLarger++;
                }
            } else if (goodMove != -2 && trust[jNon] > 10.0*deltaTolerance) {
                trust[jNon] *= 0.2;
                numberSmaller++;
            }
            maxGap = CoinMax(maxGap, trust[jNon]);
        }
#ifdef CLP_DEBUG
        if (logLevel&32)
            std::cout << "largest gap is " << maxGap << " "
                      << numberSmaller + numberSmaller2 << " reduced ("
                      << numberSmaller << " badMove ), "
                      << numberLarger << " increased" << std::endl;
#endif
        if (iPass > 10000) {
            for (jNon = 0; jNon < numberNonLinearColumns; jNon++)
                trust[jNon] *= 0.0001;
        }
        printf("last good %d goodMove %d\n", lastGoodMove, goodMove);
        if (goodMove > 0) {
            double drop = lastObjective - objValue;
            printf("Pass %d, objective %g - drop %g maxDelta %g\n", iPass, objValue, drop, maxDelta);
            if (iPass > 20 && drop < 1.0e-12*fabs(objValue) && lastGoodMove > 0)
                drop = 0.999e-4; // so will exit
            if (maxDelta < deltaTolerance && drop < 1.0e-4 && goodMove && theta<0.99999 && lastGoodMove>0) {
                if (logLevel > 1)
                    std::cout << "Exiting as maxDelta < tolerance and small drop" << std::endl;
                break;
            }
        } else if (!numberSmaller && iPass > 1) {
            if (logLevel > 1)
                std::cout << "Exiting as all gaps small" << std::endl;
            break;
        }
        if (!iPass)
            goodMove = 1;
        targetDrop = 0.0;
        double * r = model.dualColumnSolution();
        for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
            iColumn = listNonLinearColumn[jNon];
            columnLower[iColumn] = CoinMax(solution[iColumn]
                                           - trust[jNon],
                                           trueLower[jNon]);
            columnUpper[iColumn] = CoinMin(solution[iColumn]
                                           + trust[jNon],
                                           trueUpper[jNon]);
        }
        if (iPass) {
            // get reduced costs
            model.matrix()->transposeTimes(savePi,
                                           model.dualColumnSolution());
            const double * objective = model.objective();
            for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
                iColumn = listNonLinearColumn[jNon];
                double dj = objective[iColumn] - r[iColumn];
                r[iColumn] = dj;
                if (dj < -dualTolerance)
                    targetDrop -= dj * (columnUpper[iColumn] - solution[iColumn]);
                else if (dj > dualTolerance)
                    targetDrop -= dj * (columnLower[iColumn] - solution[iColumn]);
            }
        } else {
            memset(r, 0, numberColumns*sizeof(double));
        }
#ifdef JJF_ZERO
        for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
            iColumn = listNonLinearColumn[jNon];
            if (statusCheck[iColumn] == 'L' && r[iColumn] < -1.0e-4) {
                columnLower[iColumn] = CoinMax(solution[iColumn],
                                               trueLower[jNon]);
                columnUpper[iColumn] = CoinMin(solution[iColumn]
                                               + trust[jNon],
                                               trueUpper[jNon]);
            } else if (statusCheck[iColumn] == 'U' && r[iColumn] > 1.0e-4) {
                columnLower[iColumn] = CoinMax(solution[iColumn]
                                               - trust[jNon],
                                               trueLower[jNon]);
                columnUpper[iColumn] = CoinMin(solution[iColumn],
                                               trueUpper[jNon]);
            } else {
                columnLower[iColumn] = CoinMax(solution[iColumn]
                                               - trust[jNon],
                                               trueLower[jNon]);
                columnUpper[iColumn] = CoinMin(solution[iColumn]
                                               + trust[jNon],
                                               trueUpper[jNon]);
            }
        }
#endif
        if (goodMove > 0) {
            memcpy(saveSolution, solution, numberColumns*sizeof(double));
            memcpy(saveRowSolution, rowActivity, numberRows*sizeof(double));
            memcpy(savePi, model.dualRowSolution(), numberRows*sizeof(double));
            memcpy(saveStatus, model.statusArray(), numberRows + numberColumns);

#ifdef CLP_DEBUG
            if (logLevel&32)
                std::cout << "Pass - " << iPass
                          << ", target drop is " << targetDrop
                          << std::endl;
#endif
            lastObjective = objValue;
            if (targetDrop < CoinMax(1.0e-8, CoinMin(1.0e-6, 1.0e-6*fabs(objValue))) && lastGoodMove && iPass > 3) {
                if (logLevel > 1)
                    printf("Exiting on target drop %g\n", targetDrop);
                break;
            }
#ifdef CLP_DEBUG
            {
                double * r = model.dualColumnSolution();
                for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
                    iColumn = listNonLinearColumn[jNon];
                    if (logLevel&32)
                        printf("Trust %d %g - solution %d %g obj %g dj %g state %c - bounds %g %g\n",
                               jNon, trust[jNon], iColumn, solution[iColumn], objective[iColumn],
                               r[iColumn], statusCheck[iColumn], columnLower[iColumn],
                               columnUpper[iColumn]);
                }
            }
#endif
            model.scaling(false);
            model.primal(1);
            for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
                iColumn = listNonLinearColumn[jNon];
                printf("%d bounds etc %g %g %g\n", iColumn, columnLower[iColumn], solution[iColumn], columnUpper[iColumn]);
            }
            char temp[20];
            sprintf(temp, "pass%d.mps", iPass);
            //model.writeMps(temp);
#ifdef CLP_DEBUG
            if (model.status()) {
                model.writeMps("xx.mps");
            }
#endif
            if (model.status() == 1) {
                // not feasible ! - backtrack and exit
                // use safe solution
                memcpy(solution, safeSolution, numberColumns*sizeof(double));
                memcpy(saveSolution, solution, numberColumns*sizeof(double));
                memset(rowActivity, 0, numberRows*sizeof(double));
                model.times(1.0, solution, rowActivity);
                memcpy(saveRowSolution, rowActivity, numberRows*sizeof(double));
                memcpy(model.dualRowSolution(), savePi, numberRows*sizeof(double));
                memcpy(model.statusArray(), saveStatus, numberRows + numberColumns);
                for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
                    iColumn = listNonLinearColumn[jNon];
                    columnLower[iColumn] = CoinMax(solution[iColumn]
                                                   - trust[jNon],
                                                   trueLower[jNon]);
                    columnUpper[iColumn] = CoinMin(solution[iColumn]
                                                   + trust[jNon],
                                                   trueUpper[jNon]);
                }
                break;
            } else {
                // save in case problems
                memcpy(safeSolution, solution, numberColumns*sizeof(double));
            }
            goodMove = 1;
        } else {
            // bad pass - restore solution
#ifdef CLP_DEBUG
            if (logLevel&32)
                printf("Backtracking\n");
#endif
            // load old values
            for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
                iColumn = listNonLinearColumn[jNon];
                coinModel.associateElement(coinModel.columnName(iColumn), saveSolution[iColumn]);
            }
            model.loadProblem(coinModel);
            model.addColumns(numberArtificials, NULL, NULL, objectiveArtificial,
                             startArtificial, rowArtificial, elementArtificial);
            solution = model.primalColumnSolution();
            rowActivity = model.primalRowSolution();
            memcpy(solution, saveSolution, numberColumns*sizeof(double));
            memcpy(rowActivity, saveRowSolution, numberRows*sizeof(double));
            memcpy(model.dualRowSolution(), savePi, numberRows*sizeof(double));
            memcpy(model.statusArray(), saveStatus, numberRows + numberColumns);
            columnLower = model.columnLower();
            columnUpper = model.columnUpper();
            for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
                iColumn = listNonLinearColumn[jNon];
                columnLower[iColumn] = solution[iColumn];
                columnUpper[iColumn] = solution[iColumn];
            }
            model.primal(1);
            //model.writeMps("xx.mps");
            iPass--;
            goodMove = -1;
        }
    }
    // restore solution
    memcpy(solution, saveSolution, numberColumns*sizeof(double));
    delete [] statusCheck;
    delete [] savePi;
    delete [] saveStatus;
    // load new values
    for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
        iColumn = listNonLinearColumn[jNon];
        coinModel.associateElement(coinModel.columnName(iColumn), solution[iColumn]);
    }
    double * sol2 = CoinCopyOfArray(model.primalColumnSolution(), numberColumns);
    unsigned char * status2 = CoinCopyOfArray(model.statusArray(), numberColumns);
    model.loadProblem(coinModel);
    model.addColumns(numberArtificials, NULL, NULL, objectiveArtificial,
                     startArtificial, rowArtificial, elementArtificial);
    memcpy(model.primalColumnSolution(), sol2, numberColumns*sizeof(double));
    memcpy(model.statusArray(), status2, numberColumns);
    delete [] sol2;
    delete [] status2;
    columnLower = model.columnLower();
    columnUpper = model.columnUpper();
    solution = model.primalColumnSolution();
    for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
        iColumn = listNonLinearColumn[jNon];
        columnLower[iColumn] = CoinMax(solution[iColumn],
                                       trueLower[jNon]);
        columnUpper[iColumn] = CoinMin(solution[iColumn],
                                       trueUpper[jNon]);
    }
    model.primal(1);
    for (jNon = 0; jNon < numberNonLinearColumns; jNon++) {
        iColumn = listNonLinearColumn[jNon];
        columnLower[iColumn] = trueLower[jNon];
        columnUpper[iColumn] = trueUpper[jNon];
    }
    delete [] saveSolution;
    delete [] safeSolution;
    delete [] saveRowSolution;
    for (iPass = 0; iPass < 3; iPass++)
        delete [] last[iPass];
    delete [] trust;
    delete [] trueUpper;
    delete [] trueLower;
    delete [] changeRegion;
    delete [] startArtificial;
    delete [] rowArtificial;
    delete [] elementArtificial;
    delete [] objectiveArtificial;
    delete [] listNonLinearColumn;
    delete [] whichRow;
    delete [] markNonlinear;
    return CoinCopyOfArray(solution, coinModel.numberColumns());
}
/* Solve linearized quadratic objective branch and bound.
   Return cutoff and OA cut
*/
double
OsiSolverLink::linearizedBAB(CglStored * cut)
{
    double bestObjectiveValue = COIN_DBL_MAX;
    if (quadraticModel_) {
        ClpSimplex * qp = new ClpSimplex(*quadraticModel_);
        // bounds
        int numberColumns = qp->numberColumns();
        double * lower = qp->columnLower();
        double * upper = qp->columnUpper();
        const double * lower2 = getColLower();
        const double * upper2 = getColUpper();
        for (int i = 0; i < numberColumns; i++) {
            lower[i] = CoinMax(lower[i], lower2[i]);
            upper[i] = CoinMin(upper[i], upper2[i]);
        }
        qp->nonlinearSLP(20, 1.0e-5);
        qp->primal();
        OsiSolverLinearizedQuadratic solver2(qp);
        const double * solution = NULL;
        // Reduce printout
        solver2.setHintParam(OsiDoReducePrint, true, OsiHintTry);
        CbcModel model2(solver2);
        // Now do requested saves and modifications
        CbcModel * cbcModel = & model2;
        OsiSolverInterface * osiModel = model2.solver();
        OsiClpSolverInterface * osiclpModel = dynamic_cast< OsiClpSolverInterface*> (osiModel);
        ClpSimplex * clpModel = osiclpModel->getModelPtr();

        // Set changed values

        CglProbing probing;
        probing.setMaxProbe(10);
        probing.setMaxLook(10);
        probing.setMaxElements(200);
        probing.setMaxProbeRoot(50);
        probing.setMaxLookRoot(10);
        probing.setRowCuts(3);
        probing.setUsingObjective(true);
        cbcModel->addCutGenerator(&probing, -1, "Probing", true, false, false, -100, -1, -1);
        cbcModel->cutGenerator(0)->setTiming(true);

        CglGomory gomory;
        gomory.setLimitAtRoot(512);
        cbcModel->addCutGenerator(&gomory, -98, "Gomory", true, false, false, -100, -1, -1);
        cbcModel->cutGenerator(1)->setTiming(true);

        CglKnapsackCover knapsackCover;
        cbcModel->addCutGenerator(&knapsackCover, -98, "KnapsackCover", true, false, false, -100, -1, -1);
        cbcModel->cutGenerator(2)->setTiming(true);

        CglClique clique;
        clique.setStarCliqueReport(false);
        clique.setRowCliqueReport(false);
        clique.setMinViolation(0.1);
        cbcModel->addCutGenerator(&clique, -98, "Clique", true, false, false, -100, -1, -1);
        cbcModel->cutGenerator(3)->setTiming(true);

        CglMixedIntegerRounding2 mixedIntegerRounding2;
        cbcModel->addCutGenerator(&mixedIntegerRounding2, -98, "MixedIntegerRounding2", true, false, false, -100, -1, -1);
        cbcModel->cutGenerator(4)->setTiming(true);

        CglFlowCover flowCover;
        cbcModel->addCutGenerator(&flowCover, -98, "FlowCover", true, false, false, -100, -1, -1);
        cbcModel->cutGenerator(5)->setTiming(true);

        CglTwomir twomir;
        twomir.setMaxElements(250);
        cbcModel->addCutGenerator(&twomir, -99, "Twomir", true, false, false, -100, -1, -1);
        cbcModel->cutGenerator(6)->setTiming(true);
        // For now - switch off most heuristics (because CglPreProcess is bad with QP)
#ifndef JJF_ONE
        CbcHeuristicFPump heuristicFPump(*cbcModel);
        heuristicFPump.setWhen(13);
        heuristicFPump.setMaximumPasses(20);
        heuristicFPump.setMaximumRetries(7);
        heuristicFPump.setAbsoluteIncrement(4332.64);
        cbcModel->addHeuristic(&heuristicFPump);
        heuristicFPump.setInitialWeight(1);

        CbcHeuristicLocal heuristicLocal(*cbcModel);
        heuristicLocal.setSearchType(1);
        cbcModel->addHeuristic(&heuristicLocal);

        CbcHeuristicGreedyCover heuristicGreedyCover(*cbcModel);
        cbcModel->addHeuristic(&heuristicGreedyCover);

        CbcHeuristicGreedyEquality heuristicGreedyEquality(*cbcModel);
        cbcModel->addHeuristic(&heuristicGreedyEquality);
#endif

        CbcRounding rounding(*cbcModel);
        rounding.setHeuristicName("rounding");
        cbcModel->addHeuristic(&rounding);

        cbcModel->setNumberBeforeTrust(5);
        cbcModel->setSpecialOptions(2);
        cbcModel->messageHandler()->setLogLevel(1);
        cbcModel->setMaximumCutPassesAtRoot(-100);
        cbcModel->setMaximumCutPasses(1);
        cbcModel->setMinimumDrop(0.05);
        // For branchAndBound this may help
        clpModel->defaultFactorizationFrequency();
        clpModel->setDualBound(1.0001e+08);
        clpModel->setPerturbation(50);
        osiclpModel->setSpecialOptions(193);
        osiclpModel->messageHandler()->setLogLevel(0);
        osiclpModel->setIntParam(OsiMaxNumIterationHotStart, 100);
        osiclpModel->setHintParam(OsiDoReducePrint, true, OsiHintTry);
        // You can save some time by switching off message building
        // clpModel->messagesPointer()->setDetailMessages(100,10000,(int *) NULL);

        // Solve

        cbcModel->initialSolve();
        if (clpModel->tightenPrimalBounds() != 0) {
            std::cout << "Problem is infeasible - tightenPrimalBounds!" << std::endl;
            delete qp;
            return COIN_DBL_MAX;
        }
        clpModel->dual();  // clean up
        cbcModel->initialSolve();
        cbcModel->branchAndBound();
        OsiSolverLinearizedQuadratic * solver3 = dynamic_cast<OsiSolverLinearizedQuadratic *> (model2.solver());
        assert (solver3);
        solution = solver3->bestSolution();
        bestObjectiveValue = solver3->bestObjectiveValue();
        setBestObjectiveValue(bestObjectiveValue);
        setBestSolution(solution, solver3->getNumCols());
        // if convex
        if ((specialOptions2()&4) != 0) {
            if (cbcModel_)
                cbcModel_->lockThread();
            // add OA cut
            double offset;
            double * gradient = new double [numberColumns+1];
            memcpy(gradient, qp->objectiveAsObject()->gradient(qp, solution, offset, true, 2),
                   numberColumns*sizeof(double));
            double rhs = 0.0;
            int * column = new int[numberColumns+1];
            int n = 0;
            for (int i = 0; i < numberColumns; i++) {
                double value = gradient[i];
                if (fabs(value) > 1.0e-12) {
                    gradient[n] = value;
                    rhs += value * solution[i];
                    column[n++] = i;
                }
            }
            gradient[n] = -1.0;
            column[n++] = numberColumns;
            cut->addCut(-COIN_DBL_MAX, offset + 1.0e-7, n, column, gradient);
            delete [] gradient;
            delete [] column;
            if (cbcModel_)
                cbcModel_->unlockThread();
        }
        delete qp;
        printf("obj %g\n", bestObjectiveValue);
    }
    return bestObjectiveValue;
}
/* Solves nonlinear problem from CoinModel using SLP - and then tries to get
   heuristic solution
   Returns solution array
*/
double *
OsiSolverLink::heuristicSolution(int numberPasses, double deltaTolerance, int mode)
{
    // get a solution
    CoinModel tempModel = coinModel_;
    ClpSimplex * temp = approximateSolution(tempModel, numberPasses, deltaTolerance);
    int numberColumns = coinModel_.numberColumns();
    double * solution = CoinCopyOfArray(temp->primalColumnSolution(), numberColumns);
    delete temp;
    if (mode == 0) {
        return solution;
    } else if (mode == 2) {
        const double * lower = getColLower();
        const double * upper = getColUpper();
        for (int iObject = 0; iObject < numberObjects_; iObject++) {
            OsiSimpleInteger * obj = dynamic_cast<OsiSimpleInteger *> (object_[iObject]);
            if (obj && (obj->priority() < biLinearPriority_ || biLinearPriority_ <= 0)) {
                int iColumn = obj->columnNumber();
                double value = solution[iColumn];
                value = floor(value + 0.5);
                if (fabs(value - solution[iColumn]) > 0.01) {
                    setColLower(iColumn, CoinMax(lower[iColumn], value - CoinMax(defaultBound_, 0.0)));
                    setColUpper(iColumn, CoinMin(upper[iColumn], value + CoinMax(defaultBound_, 1.0)));
                } else {
                    // could fix to integer
                    setColLower(iColumn, CoinMax(lower[iColumn], value - CoinMax(defaultBound_, 0.0)));
                    setColUpper(iColumn, CoinMin(upper[iColumn], value + CoinMax(defaultBound_, 0.0)));
                }
            }
        }
        return solution;
    }
    OsiClpSolverInterface newSolver;
    if (mode == 1) {
        // round all with priority < biLinearPriority_
        setFixedPriority(biLinearPriority_);
        // ? should we save and restore coin model
        tempModel = coinModel_;
        // solve modified problem
        char * mark = new char[numberColumns];
        memset(mark, 0, numberColumns);
        for (int iObject = 0; iObject < numberObjects_; iObject++) {
            OsiSimpleInteger * obj = dynamic_cast<OsiSimpleInteger *> (object_[iObject]);
            if (obj && obj->priority() < biLinearPriority_) {
                int iColumn = obj->columnNumber();
                double value = solution[iColumn];
                value = ceil(value - 1.0e-7);
                tempModel.associateElement(coinModel_.columnName(iColumn), value);
                mark[iColumn] = 1;
            }
            OsiBiLinear * objB = dynamic_cast<OsiBiLinear *> (object_[iObject]);
            if (objB) {
                // if one or both continuous then fix one
                if (objB->xMeshSize() < 1.0) {
                    int xColumn = objB->xColumn();
                    double value = solution[xColumn];
                    tempModel.associateElement(coinModel_.columnName(xColumn), value);
                    mark[xColumn] = 1;
                } else if (objB->yMeshSize() < 1.0) {
                    int yColumn = objB->yColumn();
                    double value = solution[yColumn];
                    tempModel.associateElement(coinModel_.columnName(yColumn), value);
                    mark[yColumn] = 1;
                }
            }
        }
        CoinModel * reOrdered = tempModel.reorder(mark);
        assert (reOrdered);
        tempModel = *reOrdered;
        delete reOrdered;
        delete [] mark;
        newSolver.loadFromCoinModel(tempModel, true);
        for (int iObject = 0; iObject < numberObjects_; iObject++) {
            OsiSimpleInteger * obj = dynamic_cast<OsiSimpleInteger *> (object_[iObject]);
            if (obj && obj->priority() < biLinearPriority_) {
                int iColumn = obj->columnNumber();
                double value = solution[iColumn];
                value = ceil(value - 1.0e-7);
                newSolver.setColLower(iColumn, value);
                newSolver.setColUpper(iColumn, value);
            }
            OsiBiLinear * objB = dynamic_cast<OsiBiLinear *> (object_[iObject]);
            if (objB) {
                // if one or both continuous then fix one
                if (objB->xMeshSize() < 1.0) {
                    int xColumn = objB->xColumn();
                    double value = solution[xColumn];
                    newSolver.setColLower(xColumn, value);
                    newSolver.setColUpper(xColumn, value);
                } else if (objB->yMeshSize() < 1.0) {
                    int yColumn = objB->yColumn();
                    double value = solution[yColumn];
                    newSolver.setColLower(yColumn, value);
                    newSolver.setColUpper(yColumn, value);
                }
            }
        }
    }
    CbcModel model(newSolver);
    // Now do requested saves and modifications
    CbcModel * cbcModel = & model;
    OsiSolverInterface * osiModel = model.solver();
    OsiClpSolverInterface * osiclpModel = dynamic_cast< OsiClpSolverInterface*> (osiModel);
    ClpSimplex * clpModel = osiclpModel->getModelPtr();
    CglProbing probing;
    probing.setMaxProbe(10);
    probing.setMaxLook(10);
    probing.setMaxElements(200);
    probing.setMaxProbeRoot(50);
    probing.setMaxLookRoot(10);
    probing.setRowCuts(3);
    probing.setRowCuts(0);
    probing.setUsingObjective(true);
    cbcModel->addCutGenerator(&probing, -1, "Probing", true, false, false, -100, -1, -1);

    CglGomory gomory;
    gomory.setLimitAtRoot(512);
    cbcModel->addCutGenerator(&gomory, -98, "Gomory", true, false, false, -100, -1, -1);

    CglKnapsackCover knapsackCover;
    cbcModel->addCutGenerator(&knapsackCover, -98, "KnapsackCover", true, false, false, -100, -1, -1);

    CglClique clique;
    clique.setStarCliqueReport(false);
    clique.setRowCliqueReport(false);
    clique.setMinViolation(0.1);
    cbcModel->addCutGenerator(&clique, -98, "Clique", true, false, false, -100, -1, -1);
    CglMixedIntegerRounding2 mixedIntegerRounding2;
    cbcModel->addCutGenerator(&mixedIntegerRounding2, -98, "MixedIntegerRounding2", true, false, false, -100, -1, -1);

    CglFlowCover flowCover;
    cbcModel->addCutGenerator(&flowCover, -98, "FlowCover", true, false, false, -100, -1, -1);

    CglTwomir twomir;
    twomir.setMaxElements(250);
    cbcModel->addCutGenerator(&twomir, -99, "Twomir", true, false, false, -100, -1, -1);
    cbcModel->cutGenerator(6)->setTiming(true);

    CbcHeuristicFPump heuristicFPump(*cbcModel);
    heuristicFPump.setWhen(1);
    heuristicFPump.setMaximumPasses(20);
    heuristicFPump.setDefaultRounding(0.5);
    cbcModel->addHeuristic(&heuristicFPump);

    CbcRounding rounding(*cbcModel);
    cbcModel->addHeuristic(&rounding);

    CbcHeuristicLocal heuristicLocal(*cbcModel);
    heuristicLocal.setSearchType(1);
    cbcModel->addHeuristic(&heuristicLocal);

    CbcHeuristicGreedyCover heuristicGreedyCover(*cbcModel);
    cbcModel->addHeuristic(&heuristicGreedyCover);

    CbcHeuristicGreedyEquality heuristicGreedyEquality(*cbcModel);
    cbcModel->addHeuristic(&heuristicGreedyEquality);

    CbcCompareDefault compare;
    cbcModel->setNodeComparison(compare);
    cbcModel->setNumberBeforeTrust(5);
    cbcModel->setSpecialOptions(2);
    cbcModel->messageHandler()->setLogLevel(1);
    cbcModel->setMaximumCutPassesAtRoot(-100);
    cbcModel->setMaximumCutPasses(1);
    cbcModel->setMinimumDrop(0.05);
    clpModel->setNumberIterations(1);
    // For branchAndBound this may help
    clpModel->defaultFactorizationFrequency();
    clpModel->setDualBound(6.71523e+07);
    clpModel->setPerturbation(50);
    osiclpModel->setSpecialOptions(193);
    osiclpModel->messageHandler()->setLogLevel(0);
    osiclpModel->setIntParam(OsiMaxNumIterationHotStart, 100);
    osiclpModel->setHintParam(OsiDoReducePrint, true, OsiHintTry);
    // You can save some time by switching off message building
    // clpModel->messagesPointer()->setDetailMessages(100,10000,(int *) NULL);
    // Solve

    cbcModel->initialSolve();
    //double cutoff = model_->getCutoff();
    if (!cbcModel_)
        cbcModel->setCutoff(1.0e50);
    else
        cbcModel->setCutoff(cbcModel_->getCutoff());
    // to change exits
    bool isFeasible = false;
    int saveLogLevel = clpModel->logLevel();
    clpModel->setLogLevel(0);
    int returnCode = 0;
    if (clpModel->tightenPrimalBounds() != 0) {
        clpModel->setLogLevel(saveLogLevel);
        returnCode = -1; // infeasible//std::cout<<"Problem is infeasible - tightenPrimalBounds!"<<std::endl;
        //clpModel->writeMps("infeas2.mps");
    } else {
        clpModel->setLogLevel(saveLogLevel);
        clpModel->dual();  // clean up
        // compute some things using problem size
        cbcModel->setMinimumDrop(CoinMin(5.0e-2,
                                         fabs(cbcModel->getMinimizationObjValue())*1.0e-3 + 1.0e-4));
        if (cbcModel->getNumCols() < 500)
            cbcModel->setMaximumCutPassesAtRoot(-100); // always do 100 if possible
        else if (cbcModel->getNumCols() < 5000)
            cbcModel->setMaximumCutPassesAtRoot(100); // use minimum drop
        else
            cbcModel->setMaximumCutPassesAtRoot(20);
        cbcModel->setMaximumCutPasses(1);
        // Hand coded preprocessing
        CglPreProcess process;
        OsiSolverInterface * saveSolver = cbcModel->solver()->clone();
        // Tell solver we are in Branch and Cut
        saveSolver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo) ;
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
        process.messageHandler()->setLogLevel(cbcModel->logLevel());
        OsiSolverInterface * solver2 =
            process.preProcessNonDefault(*saveSolver, 0, 10);
        // Tell solver we are not in Branch and Cut
        saveSolver->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo) ;
        if (solver2)
            solver2->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo) ;
        if (!solver2) {
            std::cout << "Pre-processing says infeasible!" << std::endl;
            delete saveSolver;
            returnCode = -1;
        } else {
            std::cout << "processed model has " << solver2->getNumRows()
                      << " rows, " << solver2->getNumCols()
                      << " and " << solver2->getNumElements() << std::endl;
            // we have to keep solver2 so pass clone
            solver2 = solver2->clone();
            //solver2->writeMps("intmodel");
            cbcModel->assignSolver(solver2);
            cbcModel->initialSolve();
            cbcModel->branchAndBound();
            // For best solution
            int numberColumns = newSolver.getNumCols();
            if (cbcModel->getMinimizationObjValue() < 1.0e50) {
                // post process
                process.postProcess(*cbcModel->solver());
                // Solution now back in saveSolver
                cbcModel->assignSolver(saveSolver);
                memcpy(cbcModel->bestSolution(), cbcModel->solver()->getColSolution(),
                       numberColumns*sizeof(double));
                // put back in original solver
                newSolver.setColSolution(cbcModel->bestSolution());
                isFeasible = true;
            } else {
                delete saveSolver;
            }
        }
    }
    assert (!returnCode);
    abort();
    return solution;
}
// Analyze constraints to see which are convex (quadratic)
void
OsiSolverLink::analyzeObjects()
{
    // space for starts
    int numberColumns = coinModel_.numberColumns();
    int * start = new int [numberColumns+1];
    const double * rowLower = getRowLower();
    const double * rowUpper = getRowUpper();
    for (int iNon = 0; iNon < numberNonLinearRows_; iNon++) {
        int iRow = rowNonLinear_[iNon];
        int numberElements = startNonLinear_[iNon+1] - startNonLinear_[iNon];
        // triplet arrays
        int * iColumn = new int [2*numberElements+1];
        int * jColumn = new int [2*numberElements];
        double * element = new double [2*numberElements];
        int i;
        int n = 0;
        for ( i = startNonLinear_[iNon]; i < startNonLinear_[iNon+1]; i++) {
            OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[whichNonLinear_[i]]);
            assert (obj);
            int xColumn = obj->xColumn();
            int yColumn = obj->yColumn();
            double coefficient = obj->coefficient();
            if (xColumn != yColumn) {
                iColumn[n] = xColumn;
                jColumn[n] = yColumn;
                element[n++] = coefficient;
                iColumn[n] = yColumn;
                jColumn[n] = xColumn;
                element[n++] = coefficient;
            } else {
                iColumn[n] = xColumn;
                jColumn[n] = xColumn;
                element[n++] = coefficient;
            }
        }
        // First sort in column order
        CoinSort_3(iColumn, iColumn + n, jColumn, element);
        // marker at end
        iColumn[n] = numberColumns;
        int lastI = iColumn[0];
        // compute starts
        start[0] = 0;
        for (i = 1; i < n + 1; i++) {
            if (iColumn[i] != lastI) {
                while (lastI < iColumn[i]) {
                    start[lastI+1] = i;
                    lastI++;
                }
                lastI = iColumn[i];
            }
        }
        // -1 unknown, 0 convex, 1 nonconvex
        int status = -1;
        int statusNegative = -1;
        int numberLong = 0; // number with >2 elements
        for (int k = 0; k < numberColumns; k++) {
            int first = start[k];
            int last = start[k+1];
            if (last > first) {
                int j;
                double diagonal = 0.0;
                int whichK = -1;
                for (j = first; j < last; j++) {
                    if (jColumn[j] == k) {
                        diagonal = element[j];
                        status = diagonal > 0 ? 0 : 1;
                        statusNegative = diagonal < 0 ? 0 : 1;
                        whichK = (j == first) ? j + 1 : j - 1;
                        break;
                    }
                }
                if (last == first + 1) {
                    // just one entry
                    if (!diagonal) {
                        // one off diagonal - not positive semi definite
                        status = 1;
                        statusNegative = 1;
                    }
                } else if (diagonal) {
                    if (last == first + 2) {
                        // other column and element
                        double otherElement = element[whichK];;
                        int otherColumn = jColumn[whichK];
                        double otherDiagonal = 0.0;
                        // check 2x2 determinant - unless past and 2 long
                        if (otherColumn > i || start[otherColumn+1] > start[otherColumn] + 2) {
                            for (j = start[otherColumn]; j < start[otherColumn+1]; j++) {
                                if (jColumn[j] == otherColumn) {
                                    otherDiagonal = element[j];
                                    break;
                                }
                            }
                            // determinant
                            double determinant = diagonal * otherDiagonal - otherElement * otherElement;
                            if (determinant < -1.0e-12) {
                                // not positive semi definite
                                status = 1;
                                statusNegative = 1;
                            } else if (start[otherColumn+1] > start[otherColumn] + 2 && determinant < 1.0e-12) {
                                // not positive semi definite
                                status = 1;
                                statusNegative = 1;
                            }
                        }
                    } else {
                        numberLong++;
                    }
                }
            }
        }
        if ((status == 0 || statusNegative == 0) && numberLong) {
            // need to do more work
            //printf("Needs more work\n");
        }
        assert (status > 0 || statusNegative > 0);
        if (!status) {
            convex_[iNon] = 1;
            // equality may be ok
            if (rowUpper[iRow] < 1.0e20)
                specialOptions2_ |= 8;
            else
                convex_[iNon] = 0;
        } else if (!statusNegative) {
            convex_[iNon] = -1;
            // equality may be ok
            if (rowLower[iRow] > -1.0e20)
                specialOptions2_ |= 8;
            else
                convex_[iNon] = 0;
        } else {
            convex_[iNon] = 0;
        }
        //printf("Convexity of row %d is %d\n",iRow,convex_[iNon]);
        delete [] iColumn;
        delete [] jColumn;
        delete [] element;
    }
    delete [] start;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface *
OsiSolverLink::clone(bool /*copyData*/) const
{
    //assert (copyData);
    OsiSolverLink * newModel = new OsiSolverLink(*this);
    return newModel;
}


//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
OsiSolverLink::OsiSolverLink (
    const OsiSolverLink & rhs)
        : OsiSolverInterface(rhs),
        CbcOsiSolver(rhs)
{
    gutsOfDestructor(true);
    gutsOfCopy(rhs);
    // something odd happens - try this
    OsiSolverInterface::operator=(rhs);
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
OsiSolverLink::~OsiSolverLink ()
{
    gutsOfDestructor();
}

//-------------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
OsiSolverLink &
OsiSolverLink::operator=(const OsiSolverLink & rhs)
{
    if (this != &rhs) {
        gutsOfDestructor();
        CbcOsiSolver::operator=(rhs);
        gutsOfCopy(rhs);
    }
    return *this;
}
void
OsiSolverLink::gutsOfDestructor(bool justNullify)
{
    if (!justNullify) {
        delete matrix_;
        delete originalRowCopy_;
        delete [] info_;
        delete [] bestSolution_;
        delete quadraticModel_;
        delete [] startNonLinear_;
        delete [] rowNonLinear_;
        delete [] convex_;
        delete [] whichNonLinear_;
        delete [] fixVariables_;
    }
    matrix_ = NULL;
    originalRowCopy_ = NULL;
    quadraticModel_ = NULL;
    numberNonLinearRows_ = 0;
    startNonLinear_ = NULL;
    rowNonLinear_ = NULL;
    convex_ = NULL;
    whichNonLinear_ = NULL;
    info_ = NULL;
    fixVariables_ = NULL;
    numberVariables_ = 0;
    specialOptions2_ = 0;
    objectiveRow_ = -1;
    objectiveVariable_ = -1;
    bestSolution_ = NULL;
    bestObjectiveValue_ = 1.0e100;
    defaultMeshSize_ = 1.0e-4;
    defaultBound_ = 1.0e5;
    integerPriority_ = 1000;
    biLinearPriority_ = 10000;
    numberFix_ = 0;
}
void
OsiSolverLink::gutsOfCopy(const OsiSolverLink & rhs)
{
    coinModel_ = rhs.coinModel_;
    numberVariables_ = rhs.numberVariables_;
    numberNonLinearRows_ = rhs.numberNonLinearRows_;
    specialOptions2_ = rhs.specialOptions2_;
    objectiveRow_ = rhs.objectiveRow_;
    objectiveVariable_ = rhs.objectiveVariable_;
    bestObjectiveValue_ = rhs.bestObjectiveValue_;
    defaultMeshSize_ = rhs.defaultMeshSize_;
    defaultBound_ = rhs.defaultBound_;
    integerPriority_ = rhs.integerPriority_;
    biLinearPriority_ = rhs.biLinearPriority_;
    numberFix_ = rhs.numberFix_;
    if (numberVariables_) {
        if (rhs.matrix_)
            matrix_ = new CoinPackedMatrix(*rhs.matrix_);
        else
            matrix_ = NULL;
        if (rhs.originalRowCopy_)
            originalRowCopy_ = new CoinPackedMatrix(*rhs.originalRowCopy_);
        else
            originalRowCopy_ = NULL;
        info_ = new OsiLinkedBound [numberVariables_];
        for (int i = 0; i < numberVariables_; i++) {
            info_[i] = OsiLinkedBound(rhs.info_[i]);
        }
        if (rhs.bestSolution_) {
            bestSolution_ = CoinCopyOfArray(rhs.bestSolution_, modelPtr_->getNumCols());
        } else {
            bestSolution_ = NULL;
        }
    }
    if (numberNonLinearRows_) {
        startNonLinear_ = CoinCopyOfArray(rhs.startNonLinear_, numberNonLinearRows_ + 1);
        rowNonLinear_ = CoinCopyOfArray(rhs.rowNonLinear_, numberNonLinearRows_);
        convex_ = CoinCopyOfArray(rhs.convex_, numberNonLinearRows_);
        int numberEntries = startNonLinear_[numberNonLinearRows_];
        whichNonLinear_ = CoinCopyOfArray(rhs.whichNonLinear_, numberEntries);
    }
    if (rhs.quadraticModel_) {
        quadraticModel_ = new ClpSimplex(*rhs.quadraticModel_);
    } else {
        quadraticModel_ = NULL;
    }
    fixVariables_ = CoinCopyOfArray(rhs.fixVariables_, numberFix_);
}
// Add a bound modifier
void
OsiSolverLink::addBoundModifier(bool upperBoundAffected, bool useUpperBound, int whichVariable, int whichVariableAffected,
                                double multiplier)
{
    bool found = false;
    int i;
    for ( i = 0; i < numberVariables_; i++) {
        if (info_[i].variable() == whichVariable) {
            found = true;
            break;
        }
    }
    if (!found) {
        // add in
        OsiLinkedBound * temp = new OsiLinkedBound [numberVariables_+1];
        for (int i = 0; i < numberVariables_; i++)
            temp[i] = info_[i];
        delete [] info_;
        info_ = temp;
        info_[numberVariables_++] = OsiLinkedBound(this, whichVariable, 0, NULL, NULL, NULL);
    }
    info_[i].addBoundModifier(upperBoundAffected, useUpperBound, whichVariableAffected, multiplier);
}
// Update coefficients
int
OsiSolverLink::updateCoefficients(ClpSimplex * solver, CoinPackedMatrix * matrix)
{
    double * lower = solver->columnLower();
    double * upper = solver->columnUpper();
    double * objective = solver->objective();
    int numberChanged = 0;
    for (int iObject = 0; iObject < numberObjects_; iObject++) {
        OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (object_[iObject]);
        if (obj) {
            numberChanged += obj->updateCoefficients(lower, upper, objective, matrix, &basis_);
        }
    }
    return numberChanged;
}
// Set best solution found internally
void
OsiSolverLink::setBestSolution(const double * solution, int numberColumns)
{
    delete [] bestSolution_;
    int numberColumnsThis = modelPtr_->numberColumns();
    bestSolution_ = new double [numberColumnsThis];
    CoinZeroN(bestSolution_, numberColumnsThis);
    memcpy(bestSolution_, solution, CoinMin(numberColumns, numberColumnsThis)*sizeof(double));
}
/* Two tier integer problem where when set of variables with priority
   less than this are fixed the problem becomes an easier integer problem
*/
void
OsiSolverLink::setFixedPriority(int priorityValue)
{
    delete [] fixVariables_;
    fixVariables_ = NULL;
    numberFix_ = 0;
    int i;
    for ( i = 0; i < numberObjects_; i++) {
        OsiSimpleInteger * obj = dynamic_cast<OsiSimpleInteger *> (object_[i]);
        if (obj) {
#ifndef NDEBUG
            int iColumn = obj->columnNumber();
            assert (iColumn >= 0);
#endif
            if (obj->priority() < priorityValue)
                numberFix_++;
        }
    }
    if (numberFix_) {
        specialOptions2_ |= 1;
        fixVariables_ = new int [numberFix_];
        numberFix_ = 0;
        // need to make sure coinModel_ is correct
        int numberColumns = coinModel_.numberColumns();
        char * highPriority = new char [numberColumns];
        CoinZeroN(highPriority, numberColumns);
        for ( i = 0; i < numberObjects_; i++) {
            OsiSimpleInteger * obj = dynamic_cast<OsiSimpleInteger *> (object_[i]);
            if (obj) {
                int iColumn = obj->columnNumber();
                assert (iColumn >= 0);
                if (iColumn < numberColumns) {
                    if (obj->priority() < priorityValue) {
                        object_[i] = new OsiSimpleFixedInteger(*obj);
                        delete obj;
                        fixVariables_[numberFix_++] = iColumn;
                        highPriority[iColumn] = 1;
                    }
                }
            }
        }
        CoinModel * newModel = coinModel_.reorder(highPriority);
        if (newModel) {
            coinModel_ = * newModel;
        } else {
            printf("Unable to use priorities\n");
            delete [] fixVariables_;
            fixVariables_ = NULL;
            numberFix_ = 0;
        }
        delete newModel;
        delete [] highPriority;
    }
}
// Gets correct form for a quadratic row - user to delete
CoinPackedMatrix *
OsiSolverLink::quadraticRow(int rowNumber, double * linearRow) const
{
    int numberColumns = coinModel_.numberColumns();
    CoinZeroN(linearRow, numberColumns);
    int numberElements = 0;
#ifndef NDEBUG
    int numberRows = coinModel_.numberRows();
    assert (rowNumber >= 0 && rowNumber < numberRows);
#endif
    CoinModelLink triple = coinModel_.firstInRow(rowNumber);
    while (triple.column() >= 0) {
        int iColumn = triple.column();
        const char * expr = coinModel_.getElementAsString(rowNumber, iColumn);
        if (strcmp(expr, "Numeric")) {
            // try and see which columns
            assert (strlen(expr) < 20000);
            char temp[20000];
            strcpy(temp, expr);
            char * pos = temp;
            bool ifFirst = true;
            while (*pos) {
                double value;
                int jColumn = decodeBit(pos, pos, value, ifFirst, coinModel_);
                // must be column unless first when may be linear term
                if (jColumn >= 0) {
                    numberElements++;
                } else if (jColumn == -2) {
                    linearRow[iColumn] = value;
                } else {
                    printf("bad nonlinear term %s\n", temp);
                    abort();
                }
                ifFirst = false;
            }
        } else {
            linearRow[iColumn] = coinModel_.getElement(rowNumber, iColumn);
        }
        triple = coinModel_.next(triple);
    }
    if (!numberElements) {
        return NULL;
    } else {
        int * column = new int[numberElements];
        int * column2 = new int[numberElements];
        double * element = new double[numberElements];
        numberElements = 0;
        CoinModelLink triple = coinModel_.firstInRow(rowNumber);
        while (triple.column() >= 0) {
            int iColumn = triple.column();
            const char * expr = coinModel_.getElementAsString(rowNumber, iColumn);
            if (strcmp(expr, "Numeric")) {
                // try and see which columns
                assert (strlen(expr) < 20000);
                char temp[20000];
                strcpy(temp, expr);
                char * pos = temp;
                bool ifFirst = true;
                while (*pos) {
                    double value;
                    int jColumn = decodeBit(pos, pos, value, ifFirst, coinModel_);
                    // must be column unless first when may be linear term
                    if (jColumn >= 0) {
                        column[numberElements] = iColumn;
                        column2[numberElements] = jColumn;
                        element[numberElements++] = value;
                    } else if (jColumn != -2) {
                        printf("bad nonlinear term %s\n", temp);
                        abort();
                    }
                    ifFirst = false;
                }
            }
            triple = coinModel_.next(triple);
        }
        return new CoinPackedMatrix(true, column2, column, element, numberElements);
    }
}
/*
  Problem specific
  Returns -1 if node fathomed and no solution
  0 if did nothing
  1 if node fathomed and solution
  allFixed is true if all LinkedBound variables are fixed
*/
int
OsiSolverLink::fathom(bool allFixed)
{
    int returnCode = 0;
    if (allFixed) {
        // solve anyway
        OsiClpSolverInterface::resolve();
        if (!isProvenOptimal()) {
            printf("cutoff before fathoming\n");
            return -1;
        }
        // all fixed so we can reformulate
        OsiClpSolverInterface newSolver;
        // set values
        const double * lower = modelPtr_->columnLower();
        const double * upper = modelPtr_->columnUpper();
        int i;
        for (i = 0; i < numberFix_; i++ ) {
            int iColumn = fixVariables_[i];
            double lo = lower[iColumn];
#ifndef NDEBUG
            double up = upper[iColumn];
            assert (lo == up);
#endif
            //printf("column %d fixed to %g\n",iColumn,lo);
            coinModel_.associateElement(coinModel_.columnName(iColumn), lo);
        }
        newSolver.loadFromCoinModel(coinModel_, true);
        for (i = 0; i < numberFix_; i++ ) {
            int iColumn = fixVariables_[i];
            newSolver.setColLower(iColumn, lower[iColumn]);
            newSolver.setColUpper(iColumn, lower[iColumn]);
        }
        // see if everything with objective fixed
        const double * objective = modelPtr_->objective();
        int numberColumns = newSolver.getNumCols();
        bool zeroObjective = true;
        double sum = 0.0;
        for (i = 0; i < numberColumns; i++) {
            if (upper[i] > lower[i] && objective[i]) {
                zeroObjective = false;
                break;
            } else {
                sum += lower[i] * objective[i];
            }
        }
        int fake[] = {5, 4, 3, 2, 0, 0, 0};
        bool onOptimalPath = true;
        for (i = 0; i < 7; i++) {
            if (static_cast<int> (upper[i]) != fake[i])
                onOptimalPath = false;
        }
        if (onOptimalPath)
            printf("possible\n");
        if (zeroObjective) {
            // randomize objective
            ClpSimplex * clpModel = newSolver.getModelPtr();
            const double * element = clpModel->matrix()->getMutableElements();
            //const int * row = clpModel->matrix()->getIndices();
            const CoinBigIndex * columnStart = clpModel->matrix()->getVectorStarts();
            const int * columnLength = clpModel->matrix()->getVectorLengths();
            double * objective = clpModel->objective();
            for (i = 0; i < numberColumns; i++) {
                if (clpModel->isInteger(i)) {
                    double value = 0.0;
                    for (int j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
                        value += fabs(element[j]);
                    }
                    objective[i] = value;
                }
            }
        }
        //newSolver.writeMps("xx");
        CbcModel model(newSolver);
        // Now do requested saves and modifications
        CbcModel * cbcModel = & model;
        OsiSolverInterface * osiModel = model.solver();
        OsiClpSolverInterface * osiclpModel = dynamic_cast< OsiClpSolverInterface*> (osiModel);
        ClpSimplex * clpModel = osiclpModel->getModelPtr();
        CglProbing probing;
        probing.setMaxProbe(10);
        probing.setMaxLook(10);
        probing.setMaxElements(200);
        probing.setMaxProbeRoot(50);
        probing.setMaxLookRoot(10);
        probing.setRowCuts(3);
        probing.setRowCuts(0);
        probing.setUsingObjective(true);
        cbcModel->addCutGenerator(&probing, -1, "Probing", true, false, false, -100, -1, -1);

        CglGomory gomory;
        gomory.setLimitAtRoot(512);
        cbcModel->addCutGenerator(&gomory, -98, "Gomory", true, false, false, -100, -1, -1);

        CglKnapsackCover knapsackCover;
        cbcModel->addCutGenerator(&knapsackCover, -98, "KnapsackCover", true, false, false, -100, -1, -1);

        CglClique clique;
        clique.setStarCliqueReport(false);
        clique.setRowCliqueReport(false);
        clique.setMinViolation(0.1);
        cbcModel->addCutGenerator(&clique, -98, "Clique", true, false, false, -100, -1, -1);
        CglMixedIntegerRounding2 mixedIntegerRounding2;
        cbcModel->addCutGenerator(&mixedIntegerRounding2, -98, "MixedIntegerRounding2", true, false, false, -100, -1, -1);

        CglFlowCover flowCover;
        cbcModel->addCutGenerator(&flowCover, -98, "FlowCover", true, false, false, -100, -1, -1);

        CglTwomir twomir;
        twomir.setMaxElements(250);
        cbcModel->addCutGenerator(&twomir, -99, "Twomir", true, false, false, -100, -1, -1);
        cbcModel->cutGenerator(6)->setTiming(true);

        CbcHeuristicFPump heuristicFPump(*cbcModel);
        heuristicFPump.setWhen(1);
        heuristicFPump.setMaximumPasses(20);
        heuristicFPump.setDefaultRounding(0.5);
        cbcModel->addHeuristic(&heuristicFPump);

        CbcRounding rounding(*cbcModel);
        cbcModel->addHeuristic(&rounding);

        CbcHeuristicLocal heuristicLocal(*cbcModel);
        heuristicLocal.setSearchType(1);
        cbcModel->addHeuristic(&heuristicLocal);

        CbcHeuristicGreedyCover heuristicGreedyCover(*cbcModel);
        cbcModel->addHeuristic(&heuristicGreedyCover);

        CbcHeuristicGreedyEquality heuristicGreedyEquality(*cbcModel);
        cbcModel->addHeuristic(&heuristicGreedyEquality);

        CbcCompareDefault compare;
        cbcModel->setNodeComparison(compare);
        cbcModel->setNumberBeforeTrust(5);
        cbcModel->setSpecialOptions(2);
        cbcModel->messageHandler()->setLogLevel(1);
        cbcModel->setMaximumCutPassesAtRoot(-100);
        cbcModel->setMaximumCutPasses(1);
        cbcModel->setMinimumDrop(0.05);
        clpModel->setNumberIterations(1);
        // For branchAndBound this may help
        clpModel->defaultFactorizationFrequency();
        clpModel->setDualBound(6.71523e+07);
        clpModel->setPerturbation(50);
        osiclpModel->setSpecialOptions(193);
        osiclpModel->messageHandler()->setLogLevel(0);
        osiclpModel->setIntParam(OsiMaxNumIterationHotStart, 100);
        osiclpModel->setHintParam(OsiDoReducePrint, true, OsiHintTry);
        // You can save some time by switching off message building
        // clpModel->messagesPointer()->setDetailMessages(100,10000,(int *) NULL);
        // Solve

        cbcModel->initialSolve();
        //double cutoff = model_->getCutoff();
        if (zeroObjective || !cbcModel_)
            cbcModel->setCutoff(1.0e50);
        else
            cbcModel->setCutoff(cbcModel_->getCutoff());
        // to change exits
        bool isFeasible = false;
        int saveLogLevel = clpModel->logLevel();
        clpModel->setLogLevel(0);
        if (clpModel->tightenPrimalBounds() != 0) {
            clpModel->setLogLevel(saveLogLevel);
            returnCode = -1; // infeasible//std::cout<<"Problem is infeasible - tightenPrimalBounds!"<<std::endl;
        } else {
            clpModel->setLogLevel(saveLogLevel);
            clpModel->dual();  // clean up
            // compute some things using problem size
            cbcModel->setMinimumDrop(CoinMin(5.0e-2,
                                             fabs(cbcModel->getMinimizationObjValue())*1.0e-3 + 1.0e-4));
            if (cbcModel->getNumCols() < 500)
                cbcModel->setMaximumCutPassesAtRoot(-100); // always do 100 if possible
            else if (cbcModel->getNumCols() < 5000)
                cbcModel->setMaximumCutPassesAtRoot(100); // use minimum drop
            else
                cbcModel->setMaximumCutPassesAtRoot(20);
            cbcModel->setMaximumCutPasses(1);
            // Hand coded preprocessing
            CglPreProcess process;
            OsiSolverInterface * saveSolver = cbcModel->solver()->clone();
            // Tell solver we are in Branch and Cut
            saveSolver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo) ;
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
            process.messageHandler()->setLogLevel(cbcModel->logLevel());
            OsiSolverInterface * solver2 =
                process.preProcessNonDefault(*saveSolver, 0, 10);
            // Tell solver we are not in Branch and Cut
            saveSolver->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo) ;
            if (solver2)
                solver2->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo) ;
            if (!solver2) {
                std::cout << "Pre-processing says infeasible!" << std::endl;
                delete saveSolver;
                returnCode = -1;
            } else {
                std::cout << "processed model has " << solver2->getNumRows()
                          << " rows, " << solver2->getNumCols()
                          << " and " << solver2->getNumElements() << std::endl;
                // we have to keep solver2 so pass clone
                solver2 = solver2->clone();
                //solver2->writeMps("intmodel");
                cbcModel->assignSolver(solver2);
                cbcModel->initialSolve();
                if (zeroObjective) {
                    cbcModel->setMaximumSolutions(1); // just getting a solution
#ifdef JJF_ZERO
                    OsiClpSolverInterface * osiclpModel = dynamic_cast< OsiClpSolverInterface*> (cbcModel->solver());
                    ClpSimplex * clpModel = osiclpModel->getModelPtr();
                    const double * element = clpModel->matrix()->getMutableElements();
                    //const int * row = clpModel->matrix()->getIndices();
                    const CoinBigIndex * columnStart = clpModel->matrix()->getVectorStarts();
                    const int * columnLength = clpModel->matrix()->getVectorLengths();
                    int n = clpModel->numberColumns();
                    int * sort2 = new int[n];
                    int * pri = new int[n];
                    double * sort = new double[n];
                    int i;
                    int nint = 0;
                    for (i = 0; i < n; i++) {
                        if (clpModel->isInteger(i)) {
                            double largest = 0.0;
                            for (int j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
                                largest = CoinMax(largest, fabs(element[j]));
                            }
                            sort2[nint] = nint;
                            sort[nint++] = -largest;
                        }
                    }
                    CoinSort_2(sort, sort + nint, sort2);
                    int kpri = 1;
                    double last = sort[0];
                    for (i = 0; i < nint; i++) {
                        if (sort[i] != last) {
                            kpri++;
                            last = sort[i];
                        }
                        pri[sort2[i]] = kpri;
                    }
                    cbcModel->passInPriorities(pri, false);
                    delete [] sort;
                    delete [] sort2;
                    delete [] pri;
#endif
                }
                cbcModel->branchAndBound();
                // For best solution
                int numberColumns = newSolver.getNumCols();
                if (cbcModel->getMinimizationObjValue() < 1.0e50) {
                    // post process
                    process.postProcess(*cbcModel->solver());
                    // Solution now back in saveSolver
                    cbcModel->assignSolver(saveSolver);
                    memcpy(cbcModel->bestSolution(), cbcModel->solver()->getColSolution(),
                           numberColumns*sizeof(double));
                    // put back in original solver
                    newSolver.setColSolution(cbcModel->bestSolution());
                    isFeasible = true;
                } else {
                    delete saveSolver;
                }
            }
            //const double * solution = newSolver.getColSolution();
            if (isFeasible && cbcModel->getMinimizationObjValue() < 1.0e50) {
                int numberColumns = this->getNumCols();
                int i;
                const double * solution = cbcModel->bestSolution();
                int numberColumns2 = newSolver.getNumCols();
                for (i = 0; i < numberColumns2; i++) {
                    double value = solution[i];
                    assert (fabs(value - floor(value + 0.5)) < 0.0001);
                    value = floor(value + 0.5);
                    this->setColLower(i, value);
                    this->setColUpper(i, value);
                }
                for (; i < numberColumns; i++) {
                    this->setColLower(i, 0.0);
                    this->setColUpper(i, 1.1);
                }
                // but take off cuts
                int numberRows = getNumRows();
                int numberRows2 = cbcModel_->continuousSolver()->getNumRows();

                for (i = numberRows2; i < numberRows; i++)
                    setRowBounds(i, -COIN_DBL_MAX, COIN_DBL_MAX);
                initialSolve();
                //if (!isProvenOptimal())
                //getModelPtr()->writeMps("bad.mps");
                if (isProvenOptimal()) {
                    delete [] bestSolution_;
                    bestSolution_ = CoinCopyOfArray(modelPtr_->getColSolution(), modelPtr_->getNumCols());
                    bestObjectiveValue_ = modelPtr_->objectiveValue();
                    printf("BB best value %g\n", bestObjectiveValue_);
                    returnCode = 1;
                } else {
                    printf("*** WHY BAD SOL\n");
                    returnCode = -1;
                }
            } else {
                modelPtr_->setProblemStatus(1);
                modelPtr_->setObjectiveValue(COIN_DBL_MAX);
                returnCode = -1;
            }
        }
    }
    return returnCode;
}
//#############################################################################
// Constructors, destructors  and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
OsiLinkedBound::OsiLinkedBound ()
{
    model_ = NULL;
    variable_ = -1;
    numberAffected_ = 0;
    maximumAffected_ = numberAffected_;
    affected_ = NULL;
}
// Useful Constructor
OsiLinkedBound::OsiLinkedBound(OsiSolverInterface * model, int variable,
                               int numberAffected, const int * positionL,
                               const int * positionU, const double * multiplier)
{
    model_ = model;
    variable_ = variable;
    numberAffected_ = 2 * numberAffected;
    maximumAffected_ = numberAffected_;
    if (numberAffected_) {
        affected_ = new boundElementAction[numberAffected_];
        int n = 0;
        for (int i = 0; i < numberAffected; i++) {
            // LB
            boundElementAction action;
            action.affect = 2;
            action.ubUsed = 0;
            action.type = 0;
            action.affected = positionL[i];
            action.multiplier = multiplier[i];
            affected_[n++] = action;
            // UB
            action.affect = 2;
            action.ubUsed = 1;
            action.type = 0;
            action.affected = positionU[i];
            action.multiplier = multiplier[i];
            affected_[n++] = action;
        }
    } else {
        affected_ = NULL;
    }
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
OsiLinkedBound::OsiLinkedBound (
    const OsiLinkedBound & rhs)
{
    model_ = rhs.model_;
    variable_ = rhs.variable_;
    numberAffected_ = rhs.numberAffected_;
    maximumAffected_ = rhs.maximumAffected_;
    if (numberAffected_) {
        affected_ = new boundElementAction[maximumAffected_];
        memcpy(affected_, rhs.affected_, numberAffected_*sizeof(boundElementAction));
    } else {
        affected_ = NULL;
    }
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
OsiLinkedBound::~OsiLinkedBound ()
{
    delete [] affected_;
}

//-------------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
OsiLinkedBound &
OsiLinkedBound::operator=(const OsiLinkedBound & rhs)
{
    if (this != &rhs) {
        delete [] affected_;
        model_ = rhs.model_;
        variable_ = rhs.variable_;
        numberAffected_ = rhs.numberAffected_;
        maximumAffected_ = rhs.maximumAffected_;
        if (numberAffected_) {
            affected_ = new boundElementAction[maximumAffected_];
            memcpy(affected_, rhs.affected_, numberAffected_*sizeof(boundElementAction));
        } else {
            affected_ = NULL;
        }
    }
    return *this;
}
// Add a bound modifier
void
OsiLinkedBound::addBoundModifier(bool upperBoundAffected, bool useUpperBound, int whichVariable,
                                 double multiplier)
{
    if (numberAffected_ == maximumAffected_) {
        maximumAffected_ = maximumAffected_ + 10 + maximumAffected_ / 4;
        boundElementAction * temp = new boundElementAction[maximumAffected_];
        memcpy(temp, affected_, numberAffected_*sizeof(boundElementAction));
        delete [] affected_;
        affected_ = temp;
    }
    boundElementAction action;
    action.affect = static_cast<unsigned char>(upperBoundAffected ? 1 : 0);
    action.ubUsed = static_cast<unsigned char>(useUpperBound ? 1 : 0);
    action.type = 2;
    action.affected = static_cast<short int>(whichVariable);
    action.multiplier = multiplier;
    affected_[numberAffected_++] = action;

}
// Update other bounds
void
OsiLinkedBound::updateBounds(ClpSimplex * solver)
{
    double * lower = solver->columnLower();
    double * upper = solver->columnUpper();
    double lo = lower[variable_];
    double up = upper[variable_];
    // printf("bounds for %d are %g and %g\n",variable_,lo,up);
    for (int j = 0; j < numberAffected_; j++) {
        if (affected_[j].affect < 2) {
            double multiplier = affected_[j].multiplier;
            assert (affected_[j].type == 2);
            int iColumn = affected_[j].affected;
            double useValue = (affected_[j].ubUsed) ? up : lo;
            if (affected_[j].affect == 0)
                lower[iColumn] = CoinMin(upper[iColumn], CoinMax(lower[iColumn], multiplier * useValue));
            else
                upper[iColumn] = CoinMax(lower[iColumn], CoinMin(upper[iColumn], multiplier * useValue));
        }
    }
}
#ifdef JJF_ZERO
// Add an element modifier
void
OsiLinkedBound::addCoefficientModifier(bool useUpperBound, int position,
                                       double multiplier)
{
    if (numberAffected_ == maximumAffected_) {
        maximumAffected_ = maximumAffected_ + 10 + maximumAffected_ / 4;
        boundElementAction * temp = new boundElementAction[maximumAffected_];
        memcpy(temp, affected_, numberAffected_*sizeof(boundElementAction));
        delete [] affected_;
        affected_ = temp;
    }
    boundElementAction action;
    action.affect = 2;
    action.ubUsed = useUpperBound ? 1 : 0;
    action.type = 0;
    action.affected = position;
    action.multiplier = multiplier;
    affected_[numberAffected_++] = action;

}
// Update coefficients
void
OsiLinkedBound::updateCoefficients(ClpSimplex * solver, CoinPackedMatrix * matrix)
{
    double * lower = solver->columnLower();
    double * upper = solver->columnUpper();
    double * element = matrix->getMutableElements();
    double lo = lower[variable_];
    double up = upper[variable_];
    // printf("bounds for %d are %g and %g\n",variable_,lo,up);
    for (int j = 0; j < numberAffected_; j++) {
        if (affected_[j].affect == 2) {
            double multiplier = affected_[j].multiplier;
            assert (affected_[j].type == 0);
            int position = affected_[j].affected;
            //double old = element[position];
            if (affected_[j].ubUsed)
                element[position] = multiplier * up;
            else
                element[position] = multiplier * lo;
            //if ( old != element[position])
            //printf("change at %d from %g to %g\n",position,old,element[position]);
        }
    }
}
#endif
// Default Constructor
CbcHeuristicDynamic3::CbcHeuristicDynamic3()
        : CbcHeuristic()
{
}

// Constructor from model
CbcHeuristicDynamic3::CbcHeuristicDynamic3(CbcModel & model)
        : CbcHeuristic(model)
{
}

// Destructor
CbcHeuristicDynamic3::~CbcHeuristicDynamic3 ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicDynamic3::clone() const
{
    return new CbcHeuristicDynamic3(*this);
}

// Copy constructor
CbcHeuristicDynamic3::CbcHeuristicDynamic3(const CbcHeuristicDynamic3 & rhs)
        :
        CbcHeuristic(rhs)
{
}

// Returns 1 if solution, 0 if not
int
CbcHeuristicDynamic3::solution(double & solutionValue,
                               double * betterSolution)
{
    if (!model_)
        return 0;
    OsiSolverLink * clpSolver
    = dynamic_cast<OsiSolverLink *> (model_->solver());
    assert (clpSolver);
    double newSolutionValue = clpSolver->bestObjectiveValue();
    const double * solution = clpSolver->bestSolution();
    if (newSolutionValue < solutionValue && solution) {
        int numberColumns = clpSolver->getNumCols();
        // new solution
        memcpy(betterSolution, solution, numberColumns*sizeof(double));
        solutionValue = newSolutionValue;
        return 1;
    } else {
        return 0;
    }
}
// update model
void CbcHeuristicDynamic3::setModel(CbcModel * model)
{
    model_ = model;
}
// Resets stuff if model changes
void
CbcHeuristicDynamic3::resetModel(CbcModel * model)
{
    model_ = model;
}
#include <cassert>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "OsiSolverInterface.hpp"
//#include "OsiBranchLink.hpp"
#include "CoinError.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinWarmStartBasis.hpp"

// Default Constructor
OsiOldLink::OsiOldLink ()
        : OsiSOS(),
        numberLinks_(0)
{
}

// Useful constructor (which are indices)
OsiOldLink::OsiOldLink (const OsiSolverInterface * /*solver*/,  int numberMembers,
                        int numberLinks, int first , const double * weights, int /*identifier*/)
        : OsiSOS(),
        numberLinks_(numberLinks)
{
    numberMembers_ = numberMembers;
    members_ = NULL;
    sosType_ = 1;
    if (numberMembers_) {
        weights_ = new double[numberMembers_];
        members_ = new int[numberMembers_*numberLinks_];
        if (weights) {
            memcpy(weights_, weights, numberMembers_*sizeof(double));
        } else {
            for (int i = 0; i < numberMembers_; i++)
                weights_[i] = i;
        }
        // weights must be increasing
        int i;
        double last = -COIN_DBL_MAX;
        for (i = 0; i < numberMembers_; i++) {
            assert (weights_[i] > last + 1.0e-12);
            last = weights_[i];
        }
        for (i = 0; i < numberMembers_*numberLinks_; i++) {
            members_[i] = first + i;
        }
    } else {
        weights_ = NULL;
    }
}

// Useful constructor (which are indices)
OsiOldLink::OsiOldLink (const OsiSolverInterface * /*solver*/,  int numberMembers,
                        int numberLinks, int /*sosType*/, const int * which ,
                        const double * weights, int /*identifier*/)
        : OsiSOS(),
        numberLinks_(numberLinks)
{
    numberMembers_ = numberMembers;
    members_ = NULL;
    sosType_ = 1;
    if (numberMembers_) {
        weights_ = new double[numberMembers_];
        members_ = new int[numberMembers_*numberLinks_];
        if (weights) {
            memcpy(weights_, weights, numberMembers_*sizeof(double));
        } else {
            for (int i = 0; i < numberMembers_; i++)
                weights_[i] = i;
        }
        // weights must be increasing
        int i;
        double last = -COIN_DBL_MAX;
        for (i = 0; i < numberMembers_; i++) {
            assert (weights_[i] > last + 1.0e-12);
            last = weights_[i];
        }
        for (i = 0; i < numberMembers_*numberLinks_; i++) {
            members_[i] = which[i];
        }
    } else {
        weights_ = NULL;
    }
}

// Copy constructor
OsiOldLink::OsiOldLink ( const OsiOldLink & rhs)
        : OsiSOS(rhs)
{
    numberLinks_ = rhs.numberLinks_;
    if (numberMembers_) {
        delete [] members_;
        members_ = CoinCopyOfArray(rhs.members_, numberMembers_ * numberLinks_);
    }
}

// Clone
OsiObject *
OsiOldLink::clone() const
{
    return new OsiOldLink(*this);
}

// Assignment operator
OsiOldLink &
OsiOldLink::operator=( const OsiOldLink & rhs)
{
    if (this != &rhs) {
        OsiSOS::operator=(rhs);
        delete [] members_;
        numberLinks_ = rhs.numberLinks_;
        if (numberMembers_) {
            members_ = CoinCopyOfArray(rhs.members_, numberMembers_ * numberLinks_);
        } else {
            members_ = NULL;
        }
    }
    return *this;
}

// Destructor
OsiOldLink::~OsiOldLink ()
{
}

// Infeasibility - large is 0.5
double
OsiOldLink::infeasibility(const OsiBranchingInformation * info, int & whichWay) const
{
    int j;
    int firstNonZero = -1;
    int lastNonZero = -1;
    const double * solution = info->solution_;
    //const double * lower = info->lower_;
    const double * upper = info->upper_;
    double integerTolerance = info->integerTolerance_;
    double weight = 0.0;
    double sum = 0.0;

    // check bounds etc
    double lastWeight = -1.0e100;
    int base = 0;
    for (j = 0; j < numberMembers_; j++) {
        for (int k = 0; k < numberLinks_; k++) {
            int iColumn = members_[base+k];
            if (lastWeight >= weights_[j] - 1.0e-7)
                throw CoinError("Weights too close together in OsiLink", "infeasibility", "OsiLink");
            lastWeight = weights_[j];
            double value = CoinMax(0.0, solution[iColumn]);
            sum += value;
            if (value > integerTolerance && upper[iColumn]) {
                // Possibly due to scaling a fixed variable might slip through
                if (value > upper[iColumn] + 1.0e-8) {
#ifdef OSI_DEBUG
                    printf("** Variable %d (%d) has value %g and upper bound of %g\n",
                           iColumn, j, value, upper[iColumn]);
#endif
                }
                value = CoinMin(value, upper[iColumn]);
                weight += weights_[j] * value;
                if (firstNonZero < 0)
                    firstNonZero = j;
                lastNonZero = j;
            }
        }
        base += numberLinks_;
    }
    double valueInfeasibility;
    whichWay = 1;
    whichWay_ = 1;
    if (lastNonZero - firstNonZero >= sosType_) {
        // find where to branch
        assert (sum > 0.0);
        weight /= sum;
        valueInfeasibility = lastNonZero - firstNonZero + 1;
        valueInfeasibility *= 0.5 / static_cast<double> (numberMembers_);
        //#define DISTANCE
#ifdef DISTANCE
        assert (sosType_ == 1); // code up
        /* may still be satisfied.
           For LOS type 2 we might wish to move coding around
           and keep initial info in model_ for speed
        */
        int iWhere;
        bool possible = false;
        for (iWhere = firstNonZero; iWhere <= lastNonZero; iWhere++) {
            if (fabs(weight - weights_[iWhere]) < 1.0e-8) {
                possible = true;
                break;
            }
        }
        if (possible) {
            // One could move some of this (+ arrays) into model_
            const CoinPackedMatrix * matrix = solver->getMatrixByCol();
            const double * element = matrix->getMutableElements();
            const int * row = matrix->getIndices();
            const CoinBigIndex * columnStart = matrix->getVectorStarts();
            const int * columnLength = matrix->getVectorLengths();
            const double * rowSolution = solver->getRowActivity();
            const double * rowLower = solver->getRowLower();
            const double * rowUpper = solver->getRowUpper();
            int numberRows = matrix->getNumRows();
            double * array = new double [numberRows];
            CoinZeroN(array, numberRows);
            int * which = new int [numberRows];
            int n = 0;
            int base = numberLinks_ * firstNonZero;
            for (j = firstNonZero; j <= lastNonZero; j++) {
                for (int k = 0; k < numberLinks_; k++) {
                    int iColumn = members_[base+k];
                    double value = CoinMax(0.0, solution[iColumn]);
                    if (value > integerTolerance && upper[iColumn]) {
                        value = CoinMin(value, upper[iColumn]);
                        for (int j = columnStart[iColumn]; j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                            int iRow = row[j];
                            double a = array[iRow];
                            if (a) {
                                a += value * element[j];
                                if (!a)
                                    a = 1.0e-100;
                            } else {
                                which[n++] = iRow;
                                a = value * element[j];
                                assert (a);
                            }
                            array[iRow] = a;
                        }
                    }
                }
                base += numberLinks_;
            }
            base = numberLinks_ * iWhere;
            for (int k = 0; k < numberLinks_; k++) {
                int iColumn = members_[base+k];
                const double value = 1.0;
                for (int j = columnStart[iColumn]; j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                    int iRow = row[j];
                    double a = array[iRow];
                    if (a) {
                        a -= value * element[j];
                        if (!a)
                            a = 1.0e-100;
                    } else {
                        which[n++] = iRow;
                        a = -value * element[j];
                        assert (a);
                    }
                    array[iRow] = a;
                }
            }
            for (j = 0; j < n; j++) {
                int iRow = which[j];
                // moving to point will increase row solution by this
                double distance = array[iRow];
                if (distance > 1.0e-8) {
                    if (distance + rowSolution[iRow] > rowUpper[iRow] + 1.0e-8) {
                        possible = false;
                        break;
                    }
                } else if (distance < -1.0e-8) {
                    if (distance + rowSolution[iRow] < rowLower[iRow] - 1.0e-8) {
                        possible = false;
                        break;
                    }
                }
            }
            for (j = 0; j < n; j++)
                array[which[j]] = 0.0;
            delete [] array;
            delete [] which;
            if (possible) {
                valueInfeasibility = 0.0;
                printf("possible %d %d %d\n", firstNonZero, lastNonZero, iWhere);
            }
        }
#endif
    } else {
        valueInfeasibility = 0.0; // satisfied
    }
    infeasibility_ = valueInfeasibility;
    otherInfeasibility_ = 1.0 - valueInfeasibility;
    return valueInfeasibility;
}

// This looks at solution and sets bounds to contain solution
double
OsiOldLink::feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const
{
    int j;
    int firstNonZero = -1;
    int lastNonZero = -1;
    const double * solution = info->solution_;
    const double * upper = info->upper_;
    double integerTolerance = info->integerTolerance_;
    double weight = 0.0;
    double sum = 0.0;

    int base = 0;
    for (j = 0; j < numberMembers_; j++) {
        for (int k = 0; k < numberLinks_; k++) {
            int iColumn = members_[base+k];
            double value = CoinMax(0.0, solution[iColumn]);
            sum += value;
            if (value > integerTolerance && upper[iColumn]) {
                weight += weights_[j] * value;
                if (firstNonZero < 0)
                    firstNonZero = j;
                lastNonZero = j;
            }
        }
        base += numberLinks_;
    }
#ifdef DISTANCE
    if (lastNonZero - firstNonZero > sosType_ - 1) {
        /* may still be satisfied.
           For LOS type 2 we might wish to move coding around
           and keep initial info in model_ for speed
        */
        int iWhere;
        bool possible = false;
        for (iWhere = firstNonZero; iWhere <= lastNonZero; iWhere++) {
            if (fabs(weight - weights_[iWhere]) < 1.0e-8) {
                possible = true;
                break;
            }
        }
        if (possible) {
            // One could move some of this (+ arrays) into model_
            const CoinPackedMatrix * matrix = solver->getMatrixByCol();
            const double * element = matrix->getMutableElements();
            const int * row = matrix->getIndices();
            const CoinBigIndex * columnStart = matrix->getVectorStarts();
            const int * columnLength = matrix->getVectorLengths();
            const double * rowSolution = solver->getRowActivity();
            const double * rowLower = solver->getRowLower();
            const double * rowUpper = solver->getRowUpper();
            int numberRows = matrix->getNumRows();
            double * array = new double [numberRows];
            CoinZeroN(array, numberRows);
            int * which = new int [numberRows];
            int n = 0;
            int base = numberLinks_ * firstNonZero;
            for (j = firstNonZero; j <= lastNonZero; j++) {
                for (int k = 0; k < numberLinks_; k++) {
                    int iColumn = members_[base+k];
                    double value = CoinMax(0.0, solution[iColumn]);
                    if (value > integerTolerance && upper[iColumn]) {
                        value = CoinMin(value, upper[iColumn]);
                        for (int j = columnStart[iColumn]; j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                            int iRow = row[j];
                            double a = array[iRow];
                            if (a) {
                                a += value * element[j];
                                if (!a)
                                    a = 1.0e-100;
                            } else {
                                which[n++] = iRow;
                                a = value * element[j];
                                assert (a);
                            }
                            array[iRow] = a;
                        }
                    }
                }
                base += numberLinks_;
            }
            base = numberLinks_ * iWhere;
            for (int k = 0; k < numberLinks_; k++) {
                int iColumn = members_[base+k];
                const double value = 1.0;
                for (int j = columnStart[iColumn]; j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                    int iRow = row[j];
                    double a = array[iRow];
                    if (a) {
                        a -= value * element[j];
                        if (!a)
                            a = 1.0e-100;
                    } else {
                        which[n++] = iRow;
                        a = -value * element[j];
                        assert (a);
                    }
                    array[iRow] = a;
                }
            }
            for (j = 0; j < n; j++) {
                int iRow = which[j];
                // moving to point will increase row solution by this
                double distance = array[iRow];
                if (distance > 1.0e-8) {
                    if (distance + rowSolution[iRow] > rowUpper[iRow] + 1.0e-8) {
                        possible = false;
                        break;
                    }
                } else if (distance < -1.0e-8) {
                    if (distance + rowSolution[iRow] < rowLower[iRow] - 1.0e-8) {
                        possible = false;
                        break;
                    }
                }
            }
            for (j = 0; j < n; j++)
                array[which[j]] = 0.0;
            delete [] array;
            delete [] which;
            if (possible) {
                printf("possible feas region %d %d %d\n", firstNonZero, lastNonZero, iWhere);
                firstNonZero = iWhere;
                lastNonZero = iWhere;
            }
        }
    }
#else
    assert (lastNonZero - firstNonZero < sosType_) ;
#endif
    base = 0;
    for (j = 0; j < firstNonZero; j++) {
        for (int k = 0; k < numberLinks_; k++) {
            int iColumn = members_[base+k];
            solver->setColUpper(iColumn, 0.0);
        }
        base += numberLinks_;
    }
    // skip
    base += numberLinks_;
    for (j = lastNonZero + 1; j < numberMembers_; j++) {
        for (int k = 0; k < numberLinks_; k++) {
            int iColumn = members_[base+k];
            solver->setColUpper(iColumn, 0.0);
        }
        base += numberLinks_;
    }
    // go to coding as in OsiSOS
    abort();
    return -1.0;
}

// Redoes data when sequence numbers change
void
OsiOldLink::resetSequenceEtc(int numberColumns, const int * originalColumns)
{
    int n2 = 0;
    for (int j = 0; j < numberMembers_*numberLinks_; j++) {
        int iColumn = members_[j];
        int i;
#ifdef JJF_ZERO
        for (i = 0; i < numberColumns; i++) {
            if (originalColumns[i] == iColumn)
                break;
        }
#else
        i = originalColumns[iColumn];
#endif
        if (i >= 0 && i < numberColumns) {
            members_[n2] = i;
            weights_[n2++] = weights_[j];
        }
    }
    if (n2 < numberMembers_) {
        printf("** SOS number of members reduced from %d to %d!\n", numberMembers_, n2 / numberLinks_);
        numberMembers_ = n2 / numberLinks_;
    }
}

// Creates a branching object
OsiBranchingObject *
OsiOldLink::createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const
{
    int j;
    const double * solution = info->solution_;
    double tolerance = info->primalTolerance_;
    const double * upper = info->upper_;
    int firstNonFixed = -1;
    int lastNonFixed = -1;
    int firstNonZero = -1;
    int lastNonZero = -1;
    double weight = 0.0;
    double sum = 0.0;
    int base = 0;
    for (j = 0; j < numberMembers_; j++) {
        for (int k = 0; k < numberLinks_; k++) {
            int iColumn = members_[base+k];
            if (upper[iColumn]) {
                double value = CoinMax(0.0, solution[iColumn]);
                sum += value;
                if (firstNonFixed < 0)
                    firstNonFixed = j;
                lastNonFixed = j;
                if (value > tolerance) {
                    weight += weights_[j] * value;
                    if (firstNonZero < 0)
                        firstNonZero = j;
                    lastNonZero = j;
                }
            }
        }
        base += numberLinks_;
    }
    assert (lastNonZero - firstNonZero >= sosType_) ;
    // find where to branch
    assert (sum > 0.0);
    weight /= sum;
    int iWhere;
    double separator = 0.0;
    for (iWhere = firstNonZero; iWhere < lastNonZero; iWhere++)
        if (weight < weights_[iWhere+1])
            break;
    if (sosType_ == 1) {
        // SOS 1
        separator = 0.5 * (weights_[iWhere] + weights_[iWhere+1]);
    } else {
        // SOS 2
        if (iWhere == firstNonFixed)
            iWhere++;;
        if (iWhere == lastNonFixed - 1)
            iWhere = lastNonFixed - 2;
        separator = weights_[iWhere+1];
    }
    // create object
    OsiBranchingObject * branch;
    branch = new OsiOldLinkBranchingObject(solver, this, way, separator);
    return branch;
}
OsiOldLinkBranchingObject::OsiOldLinkBranchingObject()
        : OsiSOSBranchingObject()
{
}

// Useful constructor
OsiOldLinkBranchingObject::OsiOldLinkBranchingObject (OsiSolverInterface * solver,
        const OsiOldLink * set,
        int way ,
        double separator)
        : OsiSOSBranchingObject(solver, set, way, separator)
{
}

// Copy constructor
OsiOldLinkBranchingObject::OsiOldLinkBranchingObject ( const OsiOldLinkBranchingObject & rhs) : OsiSOSBranchingObject(rhs)
{
}

// Assignment operator
OsiOldLinkBranchingObject &
OsiOldLinkBranchingObject::operator=( const OsiOldLinkBranchingObject & rhs)
{
    if (this != &rhs) {
        OsiSOSBranchingObject::operator=(rhs);
    }
    return *this;
}
OsiBranchingObject *
OsiOldLinkBranchingObject::clone() const
{
    return (new OsiOldLinkBranchingObject(*this));
}


// Destructor
OsiOldLinkBranchingObject::~OsiOldLinkBranchingObject ()
{
}
double
OsiOldLinkBranchingObject::branch(OsiSolverInterface * solver)
{
    const OsiOldLink * set =
        dynamic_cast <const OsiOldLink *>(originalObject_) ;
    assert (set);
    int way = (!branchIndex_) ? (2 * firstBranch_ - 1) : -(2 * firstBranch_ - 1);
    branchIndex_++;
    int numberMembers = set->numberMembers();
    const int * which = set->members();
    const double * weights = set->weights();
    int numberLinks = set->numberLinks();
    //const double * lower = info->lower_;
    //const double * upper = solver->getColUpper();
    // *** for way - up means fix all those in down section
    if (way < 0) {
        int i;
        for ( i = 0; i < numberMembers; i++) {
            if (weights[i] > value_)
                break;
        }
        assert (i < numberMembers);
        int base = i * numberLinks;;
        for (; i < numberMembers; i++) {
            for (int k = 0; k < numberLinks; k++) {
                int iColumn = which[base+k];
                solver->setColUpper(iColumn, 0.0);
            }
            base += numberLinks;
        }
    } else {
        int i;
        int base = 0;
        for ( i = 0; i < numberMembers; i++) {
            if (weights[i] >= value_) {
                break;
            } else {
                for (int k = 0; k < numberLinks; k++) {
                    int iColumn = which[base+k];
                    solver->setColUpper(iColumn, 0.0);
                }
                base += numberLinks;
            }
        }
        assert (i < numberMembers);
    }
    return 0.0;
}
// Print what would happen
void
OsiOldLinkBranchingObject::print(const OsiSolverInterface * solver)
{
    const OsiOldLink * set =
        dynamic_cast <const OsiOldLink *>(originalObject_) ;
    assert (set);
    int way = (!branchIndex_) ? (2 * firstBranch_ - 1) : -(2 * firstBranch_ - 1);
    int numberMembers = set->numberMembers();
    int numberLinks = set->numberLinks();
    const double * weights = set->weights();
    const int * which = set->members();
    const double * upper = solver->getColUpper();
    int first = numberMembers;
    int last = -1;
    int numberFixed = 0;
    int numberOther = 0;
    int i;
    int base = 0;
    for ( i = 0; i < numberMembers; i++) {
        for (int k = 0; k < numberLinks; k++) {
            int iColumn = which[base+k];
            double bound = upper[iColumn];
            if (bound) {
                first = CoinMin(first, i);
                last = CoinMax(last, i);
            }
        }
        base += numberLinks;
    }
    // *** for way - up means fix all those in down section
    base = 0;
    if (way < 0) {
        printf("SOS Down");
        for ( i = 0; i < numberMembers; i++) {
            if (weights[i] > value_)
                break;
            for (int k = 0; k < numberLinks; k++) {
                int iColumn = which[base+k];
                double bound = upper[iColumn];
                if (bound)
                    numberOther++;
            }
            base += numberLinks;
        }
        assert (i < numberMembers);
        for (; i < numberMembers; i++) {
            for (int k = 0; k < numberLinks; k++) {
                int iColumn = which[base+k];
                double bound = upper[iColumn];
                if (bound)
                    numberFixed++;
            }
            base += numberLinks;
        }
    } else {
        printf("SOS Up");
        for ( i = 0; i < numberMembers; i++) {
            if (weights[i] >= value_)
                break;
            for (int k = 0; k < numberLinks; k++) {
                int iColumn = which[base+k];
                double bound = upper[iColumn];
                if (bound)
                    numberFixed++;
            }
            base += numberLinks;
        }
        assert (i < numberMembers);
        for (; i < numberMembers; i++) {
            for (int k = 0; k < numberLinks; k++) {
                int iColumn = which[base+k];
                double bound = upper[iColumn];
                if (bound)
                    numberOther++;
            }
            base += numberLinks;
        }
    }
    assert ((numberFixed % numberLinks) == 0);
    assert ((numberOther % numberLinks) == 0);
    printf(" - at %g, free range %d (%g) => %d (%g), %d would be fixed, %d other way\n",
           value_, first, weights[first], last, weights[last], numberFixed / numberLinks,
           numberOther / numberLinks);
}
// Default Constructor
OsiBiLinear::OsiBiLinear ()
        : OsiObject2(),
        coefficient_(0.0),
        xMeshSize_(0.0),
        yMeshSize_(0.0),
        xSatisfied_(1.0e-6),
        ySatisfied_(1.0e-6),
        xOtherSatisfied_(0.0),
        yOtherSatisfied_(0.0),
        xySatisfied_(1.0e-6),
        xyBranchValue_(0.0),
        xColumn_(-1),
        yColumn_(-1),
        firstLambda_(-1),
        branchingStrategy_(0),
        boundType_(0),
        xRow_(-1),
        yRow_(-1),
        xyRow_(-1),
        convexity_(-1),
        numberExtraRows_(0),
        multiplier_(NULL),
        extraRow_(NULL),
        chosen_(-1)
{
}

// Useful constructor
OsiBiLinear::OsiBiLinear (OsiSolverInterface * solver, int xColumn,
                          int yColumn, int xyRow, double coefficient,
                          double xMesh, double yMesh,
                          int numberExistingObjects, const OsiObject ** objects )
        : OsiObject2(),
        coefficient_(coefficient),
        xMeshSize_(xMesh),
        yMeshSize_(yMesh),
        xSatisfied_(1.0e-6),
        ySatisfied_(1.0e-6),
        xOtherSatisfied_(0.0),
        yOtherSatisfied_(0.0),
        xySatisfied_(1.0e-6),
        xyBranchValue_(0.0),
        xColumn_(xColumn),
        yColumn_(yColumn),
        firstLambda_(-1),
        branchingStrategy_(0),
        boundType_(0),
        xRow_(-1),
        yRow_(-1),
        xyRow_(xyRow),
        convexity_(-1),
        numberExtraRows_(0),
        multiplier_(NULL),
        extraRow_(NULL),
        chosen_(-1)
{
    double columnLower[4];
    double columnUpper[4];
    double objective[4];
    double rowLower[3];
    double rowUpper[3];
    CoinBigIndex starts[5];
    int index[16];
    double element[16];
    int i;
    starts[0] = 0;
    // rows
    int numberRows = solver->getNumRows();
    // convexity
    rowLower[0] = 1.0;
    rowUpper[0] = 1.0;
    convexity_ = numberRows;
    starts[1] = 0;
    // x
    rowLower[1] = 0.0;
    rowUpper[1] = 0.0;
    index[0] = xColumn_;
    element[0] = -1.0;
    xRow_ = numberRows + 1;
    starts[2] = 1;
    int nAdd = 2;
    if (xColumn_ != yColumn_) {
        rowLower[2] = 0.0;
        rowUpper[2] = 0.0;
        index[1] = yColumn;
        element[1] = -1.0;
        nAdd = 3;
        yRow_ = numberRows + 2;
        starts[3] = 2;
    } else {
        yRow_ = -1;
        branchingStrategy_ = 1;
    }
    // may be objective
    assert (xyRow_ >= -1);
    solver->addRows(nAdd, starts, index, element, rowLower, rowUpper);
    int n = 0;
    // order is LxLy, LxUy, UxLy and UxUy
    firstLambda_ = solver->getNumCols();
    // bit sloppy as theoretically could be infeasible but otherwise need to do more work
    double xB[2];
    double yB[2];
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    xB[0] = lower[xColumn_];
    xB[1] = upper[xColumn_];
    yB[0] = lower[yColumn_];
    yB[1] = upper[yColumn_];
    if (xMeshSize_ != floor(xMeshSize_)) {
        // not integral
        xSatisfied_ = CoinMax(xSatisfied_, 0.51 * xMeshSize_);
        if (!yMeshSize_) {
            xySatisfied_ = CoinMax(xySatisfied_, xSatisfied_ * CoinMax(fabs(yB[0]), fabs(yB[1])));
        }
    }
    if (yMeshSize_ != floor(yMeshSize_)) {
        // not integral
        ySatisfied_ = CoinMax(ySatisfied_, 0.51 * yMeshSize_);
        if (!xMeshSize_) {
            xySatisfied_ = CoinMax(xySatisfied_, ySatisfied_ * CoinMax(fabs(xB[0]), fabs(xB[1])));
        }
    }
    // adjust
    double distance;
    double steps;
    if (xMeshSize_) {
        distance = xB[1] - xB[0];
        steps = floor ((distance + 0.5 * xMeshSize_) / xMeshSize_);
        distance = xB[0] + xMeshSize_ * steps;
        if (fabs(xB[1] - distance) > xSatisfied_) {
            printf("bad x mesh %g %g %g -> %g\n", xB[0], xMeshSize_, xB[1], distance);
            //double newValue = CoinMax(fabs(xB[1]-distance),xMeshSize_);
            //printf("xSatisfied increased to %g\n",newValue);
            //xSatisfied_ = newValue;
            //xB[1]=distance;
            //solver->setColUpper(xColumn_,distance);
        }
    }
    if (yMeshSize_) {
        distance = yB[1] - yB[0];
        steps = floor ((distance + 0.5 * yMeshSize_) / yMeshSize_);
        distance = yB[0] + yMeshSize_ * steps;
        if (fabs(yB[1] - distance) > ySatisfied_) {
            printf("bad y mesh %g %g %g -> %g\n", yB[0], yMeshSize_, yB[1], distance);
            //double newValue = CoinMax(fabs(yB[1]-distance),yMeshSize_);
            //printf("ySatisfied increased to %g\n",newValue);
            //ySatisfied_ = newValue;
            //yB[1]=distance;
            //solver->setColUpper(yColumn_,distance);
        }
    }
    for (i = 0; i < 4; i++) {
        double x = (i < 2) ? xB[0] : xB[1];
        double y = ((i & 1) == 0) ? yB[0] : yB[1];
        columnLower[i] = 0.0;
        columnUpper[i] = 2.0;
        objective[i] = 0.0;
        double value;
        // xy
        value = coefficient_ * x * y;
        if (xyRow_ >= 0) {
            if (fabs(value) < 1.0e-19)
                value = 1.0e-19;
            element[n] = value;
            index[n++] = xyRow_;
        } else {
            objective[i] = value;
        }
        // convexity
        value = 1.0;
        element[n] = value;
        index[n++] = 0 + numberRows;
        // x
        value = x;
        if (fabs(value) < 1.0e-19)
            value = 1.0e-19;
        element[n] = value;
        index[n++] = 1 + numberRows;
        if (xColumn_ != yColumn_) {
            // y
            value = y;
            if (fabs(value) < 1.0e-19)
                value = 1.0e-19;
            element[n] = value;
            index[n++] = 2 + numberRows;
        }
        starts[i+1] = n;
    }
    solver->addCols(4, starts, index, element, columnLower, columnUpper, objective);
    // At least one has to have a mesh
    if (!xMeshSize_ && (!yMeshSize_ || yRow_ < 0)) {
        printf("one of x and y must have a mesh size\n");
        abort();
    } else if (yRow_ >= 0) {
        if (!xMeshSize_)
            branchingStrategy_ = 2;
        else if (!yMeshSize_)
            branchingStrategy_ = 1;
    }
    // Now add constraints to link in x and or y to existing ones.
    bool xDone = false;
    bool yDone = false;
    // order is LxLy, LxUy, UxLy and UxUy
    for (i = numberExistingObjects - 1; i >= 0; i--) {
        const OsiObject * obj = objects[i];
        const OsiBiLinear * obj2 =
            dynamic_cast <const OsiBiLinear *>(obj) ;
        if (obj2) {
            if (xColumn_ == obj2->xColumn_ && !xDone) {
                // make sure y equal
                double rhs = 0.0;
                CoinBigIndex starts[2];
                int index[4];
                double element[4] = {1.0, 1.0, -1.0, -1.0};
                starts[0] = 0;
                starts[1] = 4;
                index[0] = firstLambda_ + 0;
                index[1] = firstLambda_ + 1;
                index[2] = obj2->firstLambda_ + 0;
                index[3] = obj2->firstLambda_ + 1;
                solver->addRows(1, starts, index, element, &rhs, &rhs);
                xDone = true;
            }
            if (yColumn_ == obj2->yColumn_ && yRow_ >= 0 && !yDone) {
                // make sure x equal
                double rhs = 0.0;
                CoinBigIndex starts[2];
                int index[4];
                double element[4] = {1.0, 1.0, -1.0, -1.0};
                starts[0] = 0;
                starts[1] = 4;
                index[0] = firstLambda_ + 0;
                index[1] = firstLambda_ + 2;
                index[2] = obj2->firstLambda_ + 0;
                index[3] = obj2->firstLambda_ + 2;
                solver->addRows(1, starts, index, element, &rhs, &rhs);
                yDone = true;
            }
        }
    }
}
// Set sizes and other stuff
void
OsiBiLinear::setMeshSizes(const OsiSolverInterface * solver, double x, double y)
{
    xMeshSize_ = x;
    yMeshSize_ = y;
    double xB[2];
    double yB[2];
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    xB[0] = lower[xColumn_];
    xB[1] = upper[xColumn_];
    yB[0] = lower[yColumn_];
    yB[1] = upper[yColumn_];
    if (xMeshSize_ != floor(xMeshSize_)) {
        // not integral
        xSatisfied_ = CoinMax(xSatisfied_, 0.51 * xMeshSize_);
        if (!yMeshSize_) {
            xySatisfied_ = CoinMax(xySatisfied_, xSatisfied_ * CoinMax(fabs(yB[0]), fabs(yB[1])));
        }
    }
    if (yMeshSize_ != floor(yMeshSize_)) {
        // not integral
        ySatisfied_ = CoinMax(ySatisfied_, 0.51 * yMeshSize_);
        if (!xMeshSize_) {
            xySatisfied_ = CoinMax(xySatisfied_, ySatisfied_ * CoinMax(fabs(xB[0]), fabs(xB[1])));
        }
    }
}
// Useful constructor
OsiBiLinear::OsiBiLinear (CoinModel * coinModel, int xColumn,
                          int yColumn, int xyRow, double coefficient,
                          double xMesh, double yMesh,
                          int numberExistingObjects, const OsiObject ** objects )
        : OsiObject2(),
        coefficient_(coefficient),
        xMeshSize_(xMesh),
        yMeshSize_(yMesh),
        xSatisfied_(1.0e-6),
        ySatisfied_(1.0e-6),
        xOtherSatisfied_(0.0),
        yOtherSatisfied_(0.0),
        xySatisfied_(1.0e-6),
        xyBranchValue_(0.0),
        xColumn_(xColumn),
        yColumn_(yColumn),
        firstLambda_(-1),
        branchingStrategy_(0),
        boundType_(0),
        xRow_(-1),
        yRow_(-1),
        xyRow_(xyRow),
        convexity_(-1),
        numberExtraRows_(0),
        multiplier_(NULL),
        extraRow_(NULL),
        chosen_(-1)
{
    double columnLower[4];
    double columnUpper[4];
    double objective[4];
    double rowLower[3];
    double rowUpper[3];
    CoinBigIndex starts[5];
    int index[16];
    double element[16];
    int i;
    starts[0] = 0;
    // rows
    int numberRows = coinModel->numberRows();
    // convexity
    rowLower[0] = 1.0;
    rowUpper[0] = 1.0;
    convexity_ = numberRows;
    starts[1] = 0;
    // x
    rowLower[1] = 0.0;
    rowUpper[1] = 0.0;
    index[0] = xColumn_;
    element[0] = -1.0;
    xRow_ = numberRows + 1;
    starts[2] = 1;
    int nAdd = 2;
    if (xColumn_ != yColumn_) {
        rowLower[2] = 0.0;
        rowUpper[2] = 0.0;
        index[1] = yColumn;
        element[1] = -1.0;
        nAdd = 3;
        yRow_ = numberRows + 2;
        starts[3] = 2;
    } else {
        yRow_ = -1;
        branchingStrategy_ = 1;
    }
    // may be objective
    assert (xyRow_ >= -1);
    for (i = 0; i < nAdd; i++) {
        CoinBigIndex iStart = starts[i];
        coinModel->addRow(starts[i+1] - iStart, index + iStart, element + iStart, rowLower[i], rowUpper[i]);
    }
    int n = 0;
    // order is LxLy, LxUy, UxLy and UxUy
    firstLambda_ = coinModel->numberColumns();
    // bit sloppy as theoretically could be infeasible but otherwise need to do more work
    double xB[2];
    double yB[2];
    const double * lower = coinModel->columnLowerArray();
    const double * upper = coinModel->columnUpperArray();
    xB[0] = lower[xColumn_];
    xB[1] = upper[xColumn_];
    yB[0] = lower[yColumn_];
    yB[1] = upper[yColumn_];
    if (xMeshSize_ != floor(xMeshSize_)) {
        // not integral
        xSatisfied_ = CoinMax(xSatisfied_, 0.51 * xMeshSize_);
        if (!yMeshSize_) {
            xySatisfied_ = CoinMax(xySatisfied_, xSatisfied_ * CoinMax(fabs(yB[0]), fabs(yB[1])));
        }
    }
    if (yMeshSize_ != floor(yMeshSize_)) {
        // not integral
        ySatisfied_ = CoinMax(ySatisfied_, 0.51 * yMeshSize_);
        if (!xMeshSize_) {
            xySatisfied_ = CoinMax(xySatisfied_, ySatisfied_ * CoinMax(fabs(xB[0]), fabs(xB[1])));
        }
    }
    // adjust
    double distance;
    double steps;
    if (xMeshSize_) {
        distance = xB[1] - xB[0];
        steps = floor ((distance + 0.5 * xMeshSize_) / xMeshSize_);
        distance = xB[0] + xMeshSize_ * steps;
        if (fabs(xB[1] - distance) > xSatisfied_) {
            printf("bad x mesh %g %g %g -> %g\n", xB[0], xMeshSize_, xB[1], distance);
            //double newValue = CoinMax(fabs(xB[1]-distance),xMeshSize_);
            //printf("xSatisfied increased to %g\n",newValue);
            //xSatisfied_ = newValue;
            //xB[1]=distance;
            //coinModel->setColUpper(xColumn_,distance);
        }
    }
    if (yMeshSize_) {
        distance = yB[1] - yB[0];
        steps = floor ((distance + 0.5 * yMeshSize_) / yMeshSize_);
        distance = yB[0] + yMeshSize_ * steps;
        if (fabs(yB[1] - distance) > ySatisfied_) {
            printf("bad y mesh %g %g %g -> %g\n", yB[0], yMeshSize_, yB[1], distance);
            //double newValue = CoinMax(fabs(yB[1]-distance),yMeshSize_);
            //printf("ySatisfied increased to %g\n",newValue);
            //ySatisfied_ = newValue;
            //yB[1]=distance;
            //coinModel->setColUpper(yColumn_,distance);
        }
    }
    for (i = 0; i < 4; i++) {
        double x = (i < 2) ? xB[0] : xB[1];
        double y = ((i & 1) == 0) ? yB[0] : yB[1];
        columnLower[i] = 0.0;
        columnUpper[i] = 2.0;
        objective[i] = 0.0;
        double value;
        // xy
        value = coefficient_ * x * y;
        if (xyRow_ >= 0) {
            if (fabs(value) < 1.0e-19)
                value = 1.0e-19;
            element[n] = value;
            index[n++] = xyRow_;
        } else {
            objective[i] = value;
        }
        // convexity
        value = 1.0;
        element[n] = value;
        index[n++] = 0 + numberRows;
        // x
        value = x;
        if (fabs(value) < 1.0e-19)
            value = 1.0e-19;
        element[n] = value;
        index[n++] = 1 + numberRows;
        if (xColumn_ != yColumn_) {
            // y
            value = y;
            if (fabs(value) < 1.0e-19)
                value = 1.0e-19;
            element[n] = value;
            index[n++] = 2 + numberRows;
        }
        starts[i+1] = n;
    }
    for (i = 0; i < 4; i++) {
        CoinBigIndex iStart = starts[i];
        coinModel->addColumn(starts[i+1] - iStart, index + iStart, element + iStart, columnLower[i],
                             columnUpper[i], objective[i]);
    }
    // At least one has to have a mesh
    if (!xMeshSize_ && (!yMeshSize_ || yRow_ < 0)) {
        printf("one of x and y must have a mesh size\n");
        abort();
    } else if (yRow_ >= 0) {
        if (!xMeshSize_)
            branchingStrategy_ = 2;
        else if (!yMeshSize_)
            branchingStrategy_ = 1;
    }
    // Now add constraints to link in x and or y to existing ones.
    bool xDone = false;
    bool yDone = false;
    // order is LxLy, LxUy, UxLy and UxUy
    for (i = numberExistingObjects - 1; i >= 0; i--) {
        const OsiObject * obj = objects[i];
        const OsiBiLinear * obj2 =
            dynamic_cast <const OsiBiLinear *>(obj) ;
        if (obj2) {
            if (xColumn_ == obj2->xColumn_ && !xDone) {
                // make sure y equal
                double rhs = 0.0;
                CoinBigIndex starts[2];
                int index[4];
                double element[4] = {1.0, 1.0, -1.0, -1.0};
                starts[0] = 0;
                starts[1] = 4;
                index[0] = firstLambda_ + 0;
                index[1] = firstLambda_ + 1;
                index[2] = obj2->firstLambda_ + 0;
                index[3] = obj2->firstLambda_ + 1;
                coinModel->addRow(4, index, element, rhs, rhs);
                xDone = true;
            }
            if (yColumn_ == obj2->yColumn_ && yRow_ >= 0 && !yDone) {
                // make sure x equal
                double rhs = 0.0;
                CoinBigIndex starts[2];
                int index[4];
                double element[4] = {1.0, 1.0, -1.0, -1.0};
                starts[0] = 0;
                starts[1] = 4;
                index[0] = firstLambda_ + 0;
                index[1] = firstLambda_ + 2;
                index[2] = obj2->firstLambda_ + 0;
                index[3] = obj2->firstLambda_ + 2;
                coinModel->addRow(4, index, element, rhs, rhs);
                yDone = true;
            }
        }
    }
}

// Copy constructor
OsiBiLinear::OsiBiLinear ( const OsiBiLinear & rhs)
        : OsiObject2(rhs),
        coefficient_(rhs.coefficient_),
        xMeshSize_(rhs.xMeshSize_),
        yMeshSize_(rhs.yMeshSize_),
        xSatisfied_(rhs.xSatisfied_),
        ySatisfied_(rhs.ySatisfied_),
        xOtherSatisfied_(rhs.xOtherSatisfied_),
        yOtherSatisfied_(rhs.yOtherSatisfied_),
        xySatisfied_(rhs.xySatisfied_),
        xyBranchValue_(rhs.xyBranchValue_),
        xColumn_(rhs.xColumn_),
        yColumn_(rhs.yColumn_),
        firstLambda_(rhs.firstLambda_),
        branchingStrategy_(rhs.branchingStrategy_),
        boundType_(rhs.boundType_),
        xRow_(rhs.xRow_),
        yRow_(rhs.yRow_),
        xyRow_(rhs.xyRow_),
        convexity_(rhs.convexity_),
        numberExtraRows_(rhs.numberExtraRows_),
        multiplier_(NULL),
        extraRow_(NULL),
        chosen_(rhs.chosen_)
{
    if (numberExtraRows_) {
        multiplier_ = CoinCopyOfArray(rhs.multiplier_, numberExtraRows_);
        extraRow_ = CoinCopyOfArray(rhs.extraRow_, numberExtraRows_);
    }
}

// Clone
OsiObject *
OsiBiLinear::clone() const
{
    return new OsiBiLinear(*this);
}

// Assignment operator
OsiBiLinear &
OsiBiLinear::operator=( const OsiBiLinear & rhs)
{
    if (this != &rhs) {
        OsiObject2::operator=(rhs);
        coefficient_ = rhs.coefficient_;
        xMeshSize_ = rhs.xMeshSize_;
        yMeshSize_ = rhs.yMeshSize_;
        xSatisfied_ = rhs.xSatisfied_;
        ySatisfied_ = rhs.ySatisfied_;
        xOtherSatisfied_ = rhs.xOtherSatisfied_;
        yOtherSatisfied_ = rhs.yOtherSatisfied_;
        xySatisfied_ = rhs.xySatisfied_;
        xyBranchValue_ = rhs.xyBranchValue_;
        xColumn_ = rhs.xColumn_;
        yColumn_ = rhs.yColumn_;
        firstLambda_ = rhs.firstLambda_;
        branchingStrategy_ = rhs.branchingStrategy_;
        boundType_ = rhs.boundType_;
        xRow_ = rhs.xRow_;
        yRow_ = rhs.yRow_;
        xyRow_ = rhs.xyRow_;
        convexity_ = rhs.convexity_;
        numberExtraRows_ = rhs.numberExtraRows_;
        delete [] multiplier_;
        delete [] extraRow_;
        if (numberExtraRows_) {
            multiplier_ = CoinCopyOfArray(rhs.multiplier_, numberExtraRows_);
            extraRow_ = CoinCopyOfArray(rhs.extraRow_, numberExtraRows_);
        } else {
            multiplier_ = NULL;
            extraRow_ = NULL;
        }
        chosen_ = rhs.chosen_;
    }
    return *this;
}

// Destructor
OsiBiLinear::~OsiBiLinear ()
{
    delete [] multiplier_;
    delete [] extraRow_;
}
// Adds in data for extra row with variable coefficients
void
OsiBiLinear::addExtraRow(int row, double multiplier)
{
    int * tempI = new int [numberExtraRows_+1];
    double * tempD = new double [numberExtraRows_+1];
    memcpy(tempI, extraRow_, numberExtraRows_*sizeof(int));
    memcpy(tempD, multiplier_, numberExtraRows_*sizeof(double));
    tempI[numberExtraRows_] = row;
    tempD[numberExtraRows_] = multiplier;
    if (numberExtraRows_)
        assert (row > tempI[numberExtraRows_-1]);
    numberExtraRows_++;
    delete [] extraRow_;
    extraRow_ = tempI;
    delete [] multiplier_;
    multiplier_ = tempD;
}
static bool testCoarse = true;
// Infeasibility - large is 0.5
double
OsiBiLinear::infeasibility(const OsiBranchingInformation * info, int & whichWay) const
{
    // order is LxLy, LxUy, UxLy and UxUy
    double xB[2];
    double yB[2];
    xB[0] = info->lower_[xColumn_];
    xB[1] = info->upper_[xColumn_];
    yB[0] = info->lower_[yColumn_];
    yB[1] = info->upper_[yColumn_];
#ifdef JJF_ZERO
    if (info->lower_[1] <= 43.0 && info->upper_[1] >= 43.0) {
        if (info->lower_[4] <= 49.0 && info->upper_[4] >= 49.0) {
            if (info->lower_[2] <= 16.0 && info->upper_[2] >= 16.0) {
                if (info->lower_[3] <= 19.0 && info->upper_[3] >= 19.0) {
                    printf("feas %g %g %g %g  p %g t %g\n",
                           info->solution_[1],
                           info->solution_[2],
                           info->solution_[3],
                           info->solution_[4],
                           info->solution_[0],
                           info->solution_[5]);
                }
            }
        }
    }
#endif
    double x = info->solution_[xColumn_];
    x = CoinMax(x, xB[0]);
    x = CoinMin(x, xB[1]);
    double y = info->solution_[yColumn_];
    y = CoinMax(y, yB[0]);
    y = CoinMin(y, yB[1]);
    int j;
#ifndef NDEBUG
    double xLambda = 0.0;
    double yLambda = 0.0;
    if ((branchingStrategy_&4) == 0) {
        for (j = 0; j < 4; j++) {
            int iX = j >> 1;
            int iY = j & 1;
            xLambda += xB[iX] * info->solution_[firstLambda_+j];
            if (yRow_ >= 0)
                yLambda += yB[iY] * info->solution_[firstLambda_+j];
        }
    } else {
        const double * element = info->elementByColumn_;
        const int * row = info->row_;
        const CoinBigIndex * columnStart = info->columnStart_;
        const int * columnLength = info->columnLength_;
        for (j = 0; j < 4; j++) {
            int iColumn = firstLambda_ + j;
            int iStart = columnStart[iColumn];
            int iEnd = iStart + columnLength[iColumn];
            int k = iStart;
            double sol = info->solution_[iColumn];
            for (; k < iEnd; k++) {
                if (xRow_ == row[k])
                    xLambda += element[k] * sol;
                if (yRow_ == row[k])
                    yLambda += element[k] * sol;
            }
        }
    }
    assert (fabs(x - xLambda) < 1.0e-1);
    if (yRow_ >= 0)
        assert (fabs(y - yLambda) < 1.0e-1);
#endif
    // If x or y not satisfied then branch on that
    double distance;
    double steps;
    bool xSatisfied;
    double xNew = xB[0];
    if (xMeshSize_) {
        if (x < 0.5*(xB[0] + xB[1])) {
            distance = x - xB[0];
            steps = floor ((distance + 0.5 * xMeshSize_) / xMeshSize_);
            xNew = xB[0] + steps * xMeshSize_;
            assert (xNew <= xB[1] + xSatisfied_);
            xSatisfied =  (fabs(xNew - x) < xSatisfied_);
        } else {
            distance = xB[1] - x;
            steps = floor ((distance + 0.5 * xMeshSize_) / xMeshSize_);
            xNew = xB[1] - steps * xMeshSize_;
            assert (xNew >= xB[0] - xSatisfied_);
            xSatisfied =  (fabs(xNew - x) < xSatisfied_);
        }
        // but if first coarse grid then only if gap small
        if (testCoarse && (branchingStrategy_&8) != 0 && xSatisfied &&
                xB[1] - xB[0] >= xMeshSize_) {
            // but allow if fine grid would allow
            if (fabs(xNew - x) >= xOtherSatisfied_ && fabs(yB[0] - y) > yOtherSatisfied_
                    && fabs(yB[1] - y) > yOtherSatisfied_) {
                xNew = 0.5 * (xB[0] + xB[1]);
                x = xNew;
                xSatisfied = false;
            }
        }
    } else {
        xSatisfied = true;
    }
    bool ySatisfied;
    double yNew = yB[0];
    if (yMeshSize_) {
        if (y < 0.5*(yB[0] + yB[1])) {
            distance = y - yB[0];
            steps = floor ((distance + 0.5 * yMeshSize_) / yMeshSize_);
            yNew = yB[0] + steps * yMeshSize_;
            assert (yNew <= yB[1] + ySatisfied_);
            ySatisfied =  (fabs(yNew - y) < ySatisfied_);
        } else {
            distance = yB[1] - y;
            steps = floor ((distance + 0.5 * yMeshSize_) / yMeshSize_);
            yNew = yB[1] - steps * yMeshSize_;
            assert (yNew >= yB[0] - ySatisfied_);
            ySatisfied =  (fabs(yNew - y) < ySatisfied_);
        }
        // but if first coarse grid then only if gap small
        if (testCoarse && (branchingStrategy_&8) != 0 && ySatisfied &&
                yB[1] - yB[0] >= yMeshSize_) {
            // but allow if fine grid would allow
            if (fabs(yNew - y) >= yOtherSatisfied_ && fabs(xB[0] - x) > xOtherSatisfied_
                    && fabs(xB[1] - x) > xOtherSatisfied_) {
                yNew = 0.5 * (yB[0] + yB[1]);
                y = yNew;
                ySatisfied = false;
            }
        }
    } else {
        ySatisfied = true;
    }
    /* There are several possibilities
       1 - one or both are unsatisfied and branching strategy tells us what to do
       2 - both are unsatisfied and branching strategy is 0
       3 - both are satisfied but xy is not
           3a one has bounds within satisfied_ - other does not
       (or neither have but branching strategy tells us what to do)
       3b neither do - and branching strategy does not tell us
       3c both do - treat as feasible knowing another copy of object will fix
       4 - both are satisfied and xy is satisfied - as 3c
    */
    chosen_ = -1;
    xyBranchValue_ = COIN_DBL_MAX;
    whichWay_ = 0;
    double xyTrue = x * y;
    double xyLambda = 0.0;
    if ((branchingStrategy_&4) == 0) {
        for (j = 0; j < 4; j++) {
            int iX = j >> 1;
            int iY = j & 1;
            xyLambda += xB[iX] * yB[iY] * info->solution_[firstLambda_+j];
        }
    } else {
        if (xyRow_ >= 0) {
            const double * element = info->elementByColumn_;
            const int * row = info->row_;
            const CoinBigIndex * columnStart = info->columnStart_;
            const int * columnLength = info->columnLength_;
            for (j = 0; j < 4; j++) {
                int iColumn = firstLambda_ + j;
                int iStart = columnStart[iColumn];
                int iEnd = iStart + columnLength[iColumn];
                int k = iStart;
                double sol = info->solution_[iColumn];
                for (; k < iEnd; k++) {
                    if (xyRow_ == row[k])
                        xyLambda += element[k] * sol;
                }
            }
        } else {
            // objective
            const double * objective = info->objective_;
            for (j = 0; j < 4; j++) {
                int iColumn = firstLambda_ + j;
                double sol = info->solution_[iColumn];
                xyLambda += objective[iColumn] * sol;
            }
        }
        xyLambda /= coefficient_;
    }
    if (0) {
        // only true with positive values
        // see if all convexification constraints OK with true
        assert (xyTrue + 1.0e-5 > xB[0]*y + yB[0]*x - xB[0]*yB[0]);
        assert (xyTrue + 1.0e-5 > xB[1]*y + yB[1]*x - xB[1]*yB[1]);
        assert (xyTrue - 1.0e-5 < xB[1]*y + yB[0]*x - xB[1]*yB[0]);
        assert (xyTrue - 1.0e-5 < xB[0]*y + yB[1]*x - xB[0]*yB[1]);
        // see if all convexification constraints OK with lambda version
#ifndef JJF_ONE
        assert (xyLambda + 1.0e-5 > xB[0]*y + yB[0]*x - xB[0]*yB[0]);
        assert (xyLambda + 1.0e-5 > xB[1]*y + yB[1]*x - xB[1]*yB[1]);
        assert (xyLambda - 1.0e-5 < xB[1]*y + yB[0]*x - xB[1]*yB[0]);
        assert (xyLambda - 1.0e-5 < xB[0]*y + yB[1]*x - xB[0]*yB[1]);
#endif
        // see if other bound stuff true
        assert (xyLambda + 1.0e-5 > xB[0]*y);
        assert (xyLambda + 1.0e-5 > yB[0]*x);
        assert (xyLambda - 1.0e-5 < xB[1]*y);
        assert (xyLambda - 1.0e-5 < yB[1]*x);
#define SIZE 2
        if (yColumn_ == xColumn_ + SIZE) {
#if SIZE==6
            double bMax = 2200.0;
            double bMin = bMax - 100.0;
            double b[] = {330.0, 360.0, 380.0, 430.0, 490.0, 530.0};
#elif SIZE==2
            double bMax = 1900.0;
            double bMin = bMax - 200.0;
            double b[] = {460.0, 570.0};
#else
            abort();
#endif
            double sum = 0.0;
            double sum2 = 0.0;
            int m = xColumn_;
            double x = info->solution_[m];
            double xB[2];
            double yB[2];
            xB[0] = info->lower_[m];
            xB[1] = info->upper_[m];
            for (int i = 0; i < SIZE*SIZE; i += SIZE) {
                int n = i + SIZE + m;
                double y = info->solution_[n];
                yB[0] = info->lower_[n];
                yB[1] = info->upper_[n];
                int firstLambda = SIZE * SIZE + 2 * SIZE + 4 * i + 4 * m;
                double xyLambda = 0.0;
                for (int j = 0; j < 4; j++) {
                    int iX = j >> 1;
                    int iY = j & 1;
                    xyLambda += xB[iX] * yB[iY] * info->solution_[firstLambda+j];
                }
                sum += xyLambda * b[i/SIZE];
                double xyTrue = x * y;
                sum2 += xyTrue * b[i/SIZE];
            }
            if (sum > bMax*x + 1.0e-5 || sum < bMin*x - 1.0e-5) {
                //if (sum<bMax*x+1.0e-5&&sum>bMin*x-1.0e-5) {
                printf("bmin*x %g b*w %g bmax*x %g (true) %g\n", bMin*x, sum, bMax*x, sum2);
                printf("m %d lb %g value %g up %g\n",
                       m, xB[0], x, xB[1]);
                sum = 0.0;
                for (int i = 0; i < SIZE*SIZE; i += SIZE) {
                    int n = i + SIZE + m;
                    double y = info->solution_[n];
                    yB[0] = info->lower_[n];
                    yB[1] = info->upper_[n];
                    printf("n %d lb %g value %g up %g\n",
                           n, yB[0], y, yB[1]);
                    int firstLambda = SIZE * SIZE + 2 * SIZE + 4 * i + m * 4;
                    double xyLambda = 0.0;
                    for (int j = 0; j < 4; j++) {
                        int iX = j >> 1;
                        int iY = j & 1;
                        xyLambda += xB[iX] * yB[iY] * info->solution_[firstLambda+j];
                        printf("j %d l %d new xylambda %g ", j, firstLambda + j, xyLambda);
                    }
                    sum += xyLambda * b[i/SIZE];
                    printf(" - sum now %g\n", sum);
                }
            }
            if (sum2 > bMax*x + 1.0e-5 || sum2 < bMin*x - 1.0e-5) {
                printf("bmin*x %g b*x*y %g bmax*x %g (estimate) %g\n", bMin*x, sum2, bMax*x, sum);
                printf("m %d lb %g value %g up %g\n",
                       m, xB[0], x, xB[1]);
                sum2 = 0.0;
                for (int i = 0; i < SIZE*SIZE; i += SIZE) {
                    int n = i + SIZE + m;
                    double y = info->solution_[n];
                    yB[0] = info->lower_[n];
                    yB[1] = info->upper_[n];
                    printf("n %d lb %g value %g up %g\n",
                           n, yB[0], y, yB[1]);
                    double xyTrue = x * y;
                    sum2 += xyTrue * b[i/SIZE];
                    printf("xyTrue %g - sum now %g\n", xyTrue, sum2);
                }
            }
        }
    }
    // If pseudo shadow prices then see what would happen
    //double pseudoEstimate = 0.0;
    if (info->defaultDual_ >= 0.0) {
        // If we move to xy then we move by coefficient * (xyTrue-xyLambda) on row xyRow_
        double move = xyTrue - xyLambda;
        assert (xyRow_ >= 0);
        if (boundType_ == 0) {
            move *= coefficient_;
            move *= info->pi_[xyRow_];
            move = CoinMax(move, 0.0);
        } else if (boundType_ == 1) {
            // if OK then say satisfied
        } else if (boundType_ == 2) {
        } else {
            // == row so move x and y not xy
        }
    }
    if ((branchingStrategy_&16) != 0) {
        // always treat as satisfied!!
        xSatisfied = true;
        ySatisfied = true;
        xyTrue = xyLambda;
    }
    if ( !xSatisfied) {
        if (!ySatisfied) {
            if ((branchingStrategy_&3) == 0) {
                // If pseudo shadow prices then see what would happen
                if (info->defaultDual_ >= 0.0) {
                    // need coding here
                    if (fabs(x - xNew) > fabs(y - yNew)) {
                        chosen_ = 0;
                        xyBranchValue_ = x;
                    } else {
                        chosen_ = 1;
                        xyBranchValue_ = y;
                    }
                } else {
                    if (fabs(x - xNew) > fabs(y - yNew)) {
                        chosen_ = 0;
                        xyBranchValue_ = x;
                    } else {
                        chosen_ = 1;
                        xyBranchValue_ = y;
                    }
                }
            } else if ((branchingStrategy_&3) == 1) {
                chosen_ = 0;
                xyBranchValue_ = x;
            } else {
                chosen_ = 1;
                xyBranchValue_ = y;
            }
        } else {
            // y satisfied
            chosen_ = 0;
            xyBranchValue_ = x;
        }
    } else {
        // x satisfied
        if (!ySatisfied) {
            chosen_ = 1;
            xyBranchValue_ = y;
        } else {
            /*
            3 - both are satisfied but xy is not
            3a one has bounds within satisfied_ - other does not
            (or neither have but branching strategy tells us what to do)
            3b neither do - and branching strategy does not tell us
            3c both do - treat as feasible knowing another copy of object will fix
              4 - both are satisfied and xy is satisfied - as 3c
            */
            if (fabs(xyLambda - xyTrue) < xySatisfied_ || (xB[0] == xB[1] && yB[0] == yB[1])) {
                // satisfied
#ifdef JJF_ZERO
                printf("all satisfied true %g lambda %g\n",
                       xyTrue, xyLambda);
                printf("x %d (%g,%g,%g) y %d (%g,%g,%g)\n",
                       xColumn_, xB[0], x, xB[1],
                       yColumn_, yB[0], y, yB[1]);
#endif
            } else {
                // May be infeasible - check
                bool feasible = true;
                if (xB[0] == xB[1] && yB[0] == yB[1]) {
                    double lambda[4];
                    computeLambdas(info->solver_, lambda);
                    for (int j = 0; j < 4; j++) {
                        int iColumn = firstLambda_ + j;
                        if (info->lower_[iColumn] > lambda[j] + 1.0e-5 ||
                                info->upper_[iColumn] < lambda[j] - 1.0e-5)
                            feasible = false;
                    }
                }
                if (testCoarse && (branchingStrategy_&8) != 0 && xB[1] - xB[0] < 1.0001*xSatisfied_ &&
                        yB[1] - yB[0] < 1.0001*ySatisfied_)
                    feasible = true;
                if (feasible) {
                    if (xB[1] - xB[0] >= xSatisfied_ && xMeshSize_) {
                        if (yB[1] - yB[0] >= ySatisfied_ && yMeshSize_) {
                            if ((branchingStrategy_&3) == 0) {
                                // If pseudo shadow prices then see what would happen
                                if (info->defaultDual_ >= 0.0) {
                                    // need coding here
                                    if (xB[1] - xB[0] > yB[1] - yB[0]) {
                                        chosen_ = 0;
                                        xyBranchValue_ = 0.5 * (xB[0] + xB[1]);
                                    } else {
                                        chosen_ = 1;
                                        xyBranchValue_ = 0.5 * (yB[0] + yB[1]);
                                    }
                                } else {
                                    if (xB[1] - xB[0] > yB[1] - yB[0]) {
                                        chosen_ = 0;
                                        xyBranchValue_ = 0.5 * (xB[0] + xB[1]);
                                    } else {
                                        chosen_ = 1;
                                        xyBranchValue_ = 0.5 * (yB[0] + yB[1]);
                                    }
                                }
                            } else if ((branchingStrategy_&3) == 1) {
                                chosen_ = 0;
                                xyBranchValue_ = 0.5 * (xB[0] + xB[1]);
                            } else {
                                chosen_ = 1;
                                xyBranchValue_ = 0.5 * (yB[0] + yB[1]);
                            }
                        } else {
                            // y satisfied
                            chosen_ = 0;
                            xyBranchValue_ = 0.5 * (xB[0] + xB[1]);
                        }
                    } else if (yB[1] - yB[0] >= ySatisfied_ && yMeshSize_) {
                        chosen_ = 1;
                        xyBranchValue_ = 0.5 * (yB[0] + yB[1]);
                    } else {
                        // treat as satisfied unless no coefficient tightening
                        if ((branchingStrategy_&4) != 0) {
                            chosen_ = 0; // fix up in branch
                            xyBranchValue_ = x;
                        }
                    }
                } else {
                    // node not feasible!!!
                    chosen_ = 0;
                    infeasibility_ = COIN_DBL_MAX;
                    otherInfeasibility_ = COIN_DBL_MAX;
                    whichWay = whichWay_;
                    return infeasibility_;
                }
            }
        }
    }
    if (chosen_ == -1) {
        infeasibility_ = 0.0;
    } else if (chosen_ == 0) {
        infeasibility_ = CoinMax(fabs(xyBranchValue_ - x), 1.0e-12);
        //assert (xyBranchValue_>=info->lower_[xColumn_]&&xyBranchValue_<=info->upper_[xColumn_]);
    } else {
        infeasibility_ = CoinMax(fabs(xyBranchValue_ - y), 1.0e-12);
        //assert (xyBranchValue_>=info->lower_[yColumn_]&&xyBranchValue_<=info->upper_[yColumn_]);
    }
    if (info->defaultDual_ < 0.0) {
        // not using pseudo shadow prices
        otherInfeasibility_ = 1.0 - infeasibility_;
    } else {
        abort();
    }
    if (infeasibility_) {
        bool fixed = true;
        for (int j = 0; j < 4; j++) {
            int iColumn = firstLambda_ + j;
            //if (info->lower_[iColumn]) printf("lower of %g on %d\n",info->lower_[iColumn],iColumn);
            if (info->lower_[iColumn] < info->upper_[iColumn])
                fixed = false;
        }
        if (fixed) {
            //printf("must be tolerance problem - xy true %g lambda %g\n",xyTrue,xyLambda);
            chosen_ = -1;
            infeasibility_ = 0.0;
        }
    }
    whichWay = whichWay_;
    //if (infeasibility_&&priority_==10)
    //printf("x %d %g %g %g, y %d %g %g %g\n",xColumn_,xB[0],x,xB[1],yColumn_,yB[0],y,yB[1]);
    return infeasibility_;
}
// Sets infeasibility and other when pseudo shadow prices
void
OsiBiLinear::getPseudoShadow(const OsiBranchingInformation * info)
{
    // order is LxLy, LxUy, UxLy and UxUy
    double xB[2];
    double yB[2];
    xB[0] = info->lower_[xColumn_];
    xB[1] = info->upper_[xColumn_];
    yB[0] = info->lower_[yColumn_];
    yB[1] = info->upper_[yColumn_];
    double x = info->solution_[xColumn_];
    x = CoinMax(x, xB[0]);
    x = CoinMin(x, xB[1]);
    double y = info->solution_[yColumn_];
    y = CoinMax(y, yB[0]);
    y = CoinMin(y, yB[1]);
    int j;
    double xyTrue = x * y;
    double xyLambda = 0.0;
    if ((branchingStrategy_&4) == 0) {
        for (j = 0; j < 4; j++) {
            int iX = j >> 1;
            int iY = j & 1;
            xyLambda += xB[iX] * yB[iY] * info->solution_[firstLambda_+j];
        }
    } else {
        if (xyRow_ >= 0) {
            const double * element = info->elementByColumn_;
            const int * row = info->row_;
            const CoinBigIndex * columnStart = info->columnStart_;
            const int * columnLength = info->columnLength_;
            for (j = 0; j < 4; j++) {
                int iColumn = firstLambda_ + j;
                int iStart = columnStart[iColumn];
                int iEnd = iStart + columnLength[iColumn];
                int k = iStart;
                double sol = info->solution_[iColumn];
                for (; k < iEnd; k++) {
                    if (xyRow_ == row[k])
                        xyLambda += element[k] * sol;
                }
            }
        } else {
            // objective
            const double * objective = info->objective_;
            for (j = 0; j < 4; j++) {
                int iColumn = firstLambda_ + j;
                double sol = info->solution_[iColumn];
                xyLambda += objective[iColumn] * sol;
            }
        }
        xyLambda /= coefficient_;
    }
    assert (info->defaultDual_ >= 0.0);
    // If we move to xy then we move by coefficient * (xyTrue-xyLambda) on row xyRow_
    double movement = xyTrue - xyLambda;
    infeasibility_ = 0.0;
    const double * pi = info->pi_;
    const double * activity = info->rowActivity_;
    const double * lower = info->rowLower_;
    const double * upper = info->rowUpper_;
    double tolerance = info->primalTolerance_;
    double direction = info->direction_;
    bool infeasible = false;
    if (xyRow_ >= 0) {
        assert (!boundType_);
        if (lower[xyRow_] < -1.0e20)
            assert (pi[xyRow_] <= 1.0e-3);
        if (upper[xyRow_] > 1.0e20)
            assert (pi[xyRow_] >= -1.0e-3);
        double valueP = pi[xyRow_] * direction;
        // if move makes infeasible then make at least default
        double newValue = activity[xyRow_] + movement * coefficient_;
        if (newValue > upper[xyRow_] + tolerance || newValue < lower[xyRow_] - tolerance) {
            infeasibility_ += fabs(movement * coefficient_) * CoinMax(fabs(valueP), info->defaultDual_);
            infeasible = true;
        }
    } else {
        // objective
        assert (movement > -1.0e-7);
        infeasibility_ += movement;
    }
    for (int i = 0; i < numberExtraRows_; i++) {
        int iRow = extraRow_[i];
        if (lower[iRow] < -1.0e20)
            assert (pi[iRow] <= 1.0e-3);
        if (upper[iRow] > 1.0e20)
            assert (pi[iRow] >= -1.0e-3);
        double valueP = pi[iRow] * direction;
        // if move makes infeasible then make at least default
        double newValue = activity[iRow] + movement * multiplier_[i];
        if (newValue > upper[iRow] + tolerance || newValue < lower[iRow] - tolerance) {
            infeasibility_ += fabs(movement * multiplier_[i]) * CoinMax(fabs(valueP), info->defaultDual_);
            infeasible = true;
        }
    }
    if (infeasibility_ < info->integerTolerance_) {
        if (!infeasible)
            infeasibility_ = 0.0;
        else
            infeasibility_ = info->integerTolerance_;
    }
    otherInfeasibility_ = CoinMax(1.0e-12, infeasibility_ * 10.0);
}
// Gets sum of movements to correct value
double
OsiBiLinear::getMovement(const OsiBranchingInformation * info)
{
    // order is LxLy, LxUy, UxLy and UxUy
    double xB[2];
    double yB[2];
    xB[0] = info->lower_[xColumn_];
    xB[1] = info->upper_[xColumn_];
    yB[0] = info->lower_[yColumn_];
    yB[1] = info->upper_[yColumn_];
    double x = info->solution_[xColumn_];
    x = CoinMax(x, xB[0]);
    x = CoinMin(x, xB[1]);
    double y = info->solution_[yColumn_];
    y = CoinMax(y, yB[0]);
    y = CoinMin(y, yB[1]);
    int j;
    double xyTrue = x * y;
    double xyLambda = 0.0;
    if ((branchingStrategy_&4) == 0) {
        for (j = 0; j < 4; j++) {
            int iX = j >> 1;
            int iY = j & 1;
            xyLambda += xB[iX] * yB[iY] * info->solution_[firstLambda_+j];
        }
    } else {
        if (xyRow_ >= 0) {
            const double * element = info->elementByColumn_;
            const int * row = info->row_;
            const CoinBigIndex * columnStart = info->columnStart_;
            const int * columnLength = info->columnLength_;
            for (j = 0; j < 4; j++) {
                int iColumn = firstLambda_ + j;
                int iStart = columnStart[iColumn];
                int iEnd = iStart + columnLength[iColumn];
                int k = iStart;
                double sol = info->solution_[iColumn];
                for (; k < iEnd; k++) {
                    if (xyRow_ == row[k])
                        xyLambda += element[k] * sol;
                }
            }
        } else {
            // objective
            const double * objective = info->objective_;
            for (j = 0; j < 4; j++) {
                int iColumn = firstLambda_ + j;
                double sol = info->solution_[iColumn];
                xyLambda += objective[iColumn] * sol;
            }
        }
        xyLambda /= coefficient_;
    }
    // If we move to xy then we move by coefficient * (xyTrue-xyLambda) on row xyRow_
    double movement = xyTrue - xyLambda;
    double mesh = CoinMax(xMeshSize_, yMeshSize_);
    if (fabs(movement) < xySatisfied_ && (xB[1] - xB[0] < mesh || yB[1] - yB[0] < mesh))
        return 0.0; // say feasible
    const double * activity = info->rowActivity_;
    const double * lower = info->rowLower_;
    const double * upper = info->rowUpper_;
    double tolerance = info->primalTolerance_;
    double  infeasibility = 0.0;
    if (xyRow_ >= 0) {
        assert (!boundType_);
        // if move makes infeasible
        double newValue = activity[xyRow_] + movement * coefficient_;
        if (newValue > upper[xyRow_] + tolerance)
            infeasibility += newValue - upper[xyRow_];
        else if (newValue < lower[xyRow_] - tolerance)
            infeasibility += lower[xyRow_] - newValue;
    } else {
        // objective
        assert (movement > -1.0e-7);
        infeasibility += movement;
    }
    for (int i = 0; i < numberExtraRows_; i++) {
        int iRow = extraRow_[i];
        // if move makes infeasible
        double newValue = activity[iRow] + movement * multiplier_[i];
        if (newValue > upper[iRow] + tolerance)
            infeasibility += newValue - upper[iRow];
        else if (newValue < lower[iRow] - tolerance)
            infeasibility += lower[iRow] - newValue;
    }
    return infeasibility;
}
// This looks at solution and sets bounds to contain solution
double
OsiBiLinear::feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const
{
    // If another object has finer mesh ignore this
    if ((branchingStrategy_&8) != 0)
        return 0.0;
    // order is LxLy, LxUy, UxLy and UxUy
    double xB[2];
    double yB[2];
    xB[0] = info->lower_[xColumn_];
    xB[1] = info->upper_[xColumn_];
    yB[0] = info->lower_[yColumn_];
    yB[1] = info->upper_[yColumn_];
    double x = info->solution_[xColumn_];
    double y = info->solution_[yColumn_];
    int j;
#ifndef NDEBUG
    double xLambda = 0.0;
    double yLambda = 0.0;
    if ((branchingStrategy_&4) == 0) {
        for (j = 0; j < 4; j++) {
            int iX = j >> 1;
            int iY = j & 1;
            xLambda += xB[iX] * info->solution_[firstLambda_+j];
            if (yRow_ >= 0)
                yLambda += yB[iY] * info->solution_[firstLambda_+j];
        }
    } else {
        const double * element = info->elementByColumn_;
        const int * row = info->row_;
        const CoinBigIndex * columnStart = info->columnStart_;
        const int * columnLength = info->columnLength_;
        for (j = 0; j < 4; j++) {
            int iColumn = firstLambda_ + j;
            int iStart = columnStart[iColumn];
            int iEnd = iStart + columnLength[iColumn];
            int k = iStart;
            double sol = info->solution_[iColumn];
            for (; k < iEnd; k++) {
                if (xRow_ == row[k])
                    xLambda += element[k] * sol;
                if (yRow_ == row[k])
                    yLambda += element[k] * sol;
            }
        }
    }
    if (yRow_ < 0)
        yLambda = xLambda;
#ifdef JJF_ZERO
    if (fabs(x - xLambda) > 1.0e-4 ||
            fabs(y - yLambda) > 1.0e-4)
        printf("feasibleregion x %d given %g lambda %g y %d given %g lambda %g\n",
               xColumn_, x, xLambda,
               yColumn_, y, yLambda);
#endif
#endif
    double infeasibility = 0.0;
    double distance;
    double steps;
    double xNew = x;
    if (xMeshSize_) {
        distance = x - xB[0];
        if (x < 0.5*(xB[0] + xB[1])) {
            distance = x - xB[0];
            steps = floor ((distance + 0.5 * xMeshSize_) / xMeshSize_);
            xNew = xB[0] + steps * xMeshSize_;
            assert (xNew <= xB[1] + xSatisfied_);
        } else {
            distance = xB[1] - x;
            steps = floor ((distance + 0.5 * xMeshSize_) / xMeshSize_);
            xNew = xB[1] - steps * xMeshSize_;
            assert (xNew >= xB[0] - xSatisfied_);
        }
        if (xMeshSize_ < 1.0 && fabs(xNew - x) <= xSatisfied_) {
            double lo = CoinMax(xB[0], x - 0.5 * xSatisfied_);
            double up = CoinMin(xB[1], x + 0.5 * xSatisfied_);
            solver->setColLower(xColumn_, lo);
            solver->setColUpper(xColumn_, up);
        } else {
            infeasibility +=  fabs(xNew - x);
            solver->setColLower(xColumn_, xNew);
            solver->setColUpper(xColumn_, xNew);
        }
    }
    double yNew = y;
    if (yMeshSize_) {
        distance = y - yB[0];
        if (y < 0.5*(yB[0] + yB[1])) {
            distance = y - yB[0];
            steps = floor ((distance + 0.5 * yMeshSize_) / yMeshSize_);
            yNew = yB[0] + steps * yMeshSize_;
            assert (yNew <= yB[1] + ySatisfied_);
        } else {
            distance = yB[1] - y;
            steps = floor ((distance + 0.5 * yMeshSize_) / yMeshSize_);
            yNew = yB[1] - steps * yMeshSize_;
            assert (yNew >= yB[0] - ySatisfied_);
        }
        if (yMeshSize_ < 1.0 && fabs(yNew - y) <= ySatisfied_) {
            double lo = CoinMax(yB[0], y - 0.5 * ySatisfied_);
            double up = CoinMin(yB[1], y + 0.5 * ySatisfied_);
            solver->setColLower(yColumn_, lo);
            solver->setColUpper(yColumn_, up);
        } else {
            infeasibility +=  fabs(yNew - y);
            solver->setColLower(yColumn_, yNew);
            solver->setColUpper(yColumn_, yNew);
        }
    }
    if (0) {
        // temp
        solver->setColLower(xColumn_, x);
        solver->setColUpper(xColumn_, x);
        solver->setColLower(yColumn_, y);
        solver->setColUpper(yColumn_, y);
    }
    if ((branchingStrategy_&4)) {
        // fake to make correct
        double lambda[4];
        computeLambdas(solver, lambda);
        for (int j = 0; j < 4; j++) {
            int iColumn = firstLambda_ + j;
            double value = lambda[j];
            solver->setColLower(iColumn, value);
            solver->setColUpper(iColumn, value);
        }
    }
    double xyTrue = xNew * yNew;
    double xyLambda = 0.0;
    for (j = 0; j < 4; j++) {
        int iX = j >> 1;
        int iY = j & 1;
        xyLambda += xB[iX] * yB[iY] * info->solution_[firstLambda_+j];
    }
    infeasibility += fabs(xyTrue - xyLambda);
    return infeasibility;
}
// Returns true value of single xyRow coefficient
double
OsiBiLinear::xyCoefficient(const double * solution) const
{
    // If another object has finer mesh ignore this
    if ((branchingStrategy_&8) != 0)
        return 0.0;
    double x = solution[xColumn_];
    double y = solution[yColumn_];
    //printf("x (%d,%g) y (%d,%g) x*y*coefficient %g\n",
    // xColumn_,x,yColumn_,y,x*y*coefficient_);
    return x*y*coefficient_;
}

// Redoes data when sequence numbers change
void
OsiBiLinear::resetSequenceEtc(int numberColumns, const int * originalColumns)
{
    int i;
#ifdef JJF_ZERO
    for (i = 0; i < numberColumns; i++) {
        if (originalColumns[i] == firstLambda_)
            break;
    }
#else
    i = originalColumns[firstLambda_];
#endif
    if (i >= 0 && i < numberColumns) {
        firstLambda_ = i;
        for (int j = 0; j < 4; j++) {
            assert (originalColumns[j+i] - firstLambda_ == j);
        }
    } else {
        printf("lost set\n");
        abort();
    }
    // rows will be out anyway
    abort();
}

// Creates a branching object
OsiBranchingObject *
OsiBiLinear::createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * /*info*/, int way) const
{
    // create object
    OsiBranchingObject * branch;
    assert (chosen_ == 0 || chosen_ == 1);
    //if (chosen_==0)
    //assert (xyBranchValue_>=info->lower_[xColumn_]&&xyBranchValue_<=info->upper_[xColumn_]);
    //else
    //assert (xyBranchValue_>=info->lower_[yColumn_]&&xyBranchValue_<=info->upper_[yColumn_]);
    branch = new OsiBiLinearBranchingObject(solver, this, way, xyBranchValue_, chosen_);
    return branch;
}
// Does work of branching
void
OsiBiLinear::newBounds(OsiSolverInterface * solver, int way, short xOrY, double separator) const
{
    int iColumn;
    double mesh;
    double satisfied;
    if (xOrY == 0) {
        iColumn = xColumn_;
        mesh = xMeshSize_;
        satisfied = xSatisfied_;
    } else {
        iColumn = yColumn_;
        mesh = yMeshSize_;
        satisfied = ySatisfied_;
    }
    const double * columnLower = solver->getColLower();
    const double * columnUpper = solver->getColUpper();
    double lower = columnLower[iColumn];
    double distance;
    double steps;
    double zNew = separator;
    distance = separator - lower;
    assert (mesh);
    steps = floor ((distance + 0.5 * mesh) / mesh);
    if (mesh < 1.0)
        zNew = lower + steps * mesh;
    if (zNew > columnUpper[iColumn] - satisfied)
        zNew = 0.5 * (columnUpper[iColumn] - lower);
    double oldUpper = columnUpper[iColumn] ;
    double oldLower = columnLower[iColumn] ;
#ifndef NDEBUG
    int nullChange = 0;
#endif
    if (way < 0) {
        if (zNew > separator && mesh < 1.0)
            zNew -= mesh;
        double oldUpper = columnUpper[iColumn] ;
        if (zNew + satisfied >= oldUpper)
            zNew = 0.5 * (oldUpper + oldLower);
        if (mesh == 1.0)
            zNew = floor(separator);
#ifndef NDEBUG
        if (oldUpper < zNew + 1.0e-8)
            nullChange = -1;
#endif
        solver->setColUpper(iColumn, zNew);
    } else {
        if (zNew < separator && mesh < 1.0)
            zNew += mesh;
        if (zNew - satisfied <= oldLower)
            zNew = 0.5 * (oldUpper + oldLower);
        if (mesh == 1.0)
            zNew = ceil(separator);
#ifndef NDEBUG
        if (oldLower > zNew - 1.0e-8)
            nullChange = 1;
#endif
        solver->setColLower(iColumn, zNew);
    }
    if ((branchingStrategy_&4) != 0
            && columnLower[xColumn_] == columnUpper[xColumn_]
            && columnLower[yColumn_] == columnUpper[yColumn_]) {
        // fake to make correct
        double lambda[4];
        computeLambdas(solver, lambda);
        for (int j = 0; j < 4; j++) {
            int iColumn = firstLambda_ + j;
            double value = lambda[j];
#ifndef NDEBUG
            if (fabs(value - columnLower[iColumn]) > 1.0e-5 ||
                    fabs(value - columnUpper[iColumn]) > 1.0e-5)
                nullChange = 0;
#endif
            solver->setColLower(iColumn, value);
            solver->setColUpper(iColumn, value);
        }
    }
#ifndef NDEBUG
    if (nullChange)
        printf("null change on column%s %d - bounds %g,%g\n", nullChange > 0 ? "Lower" : "Upper",
               iColumn, oldLower, oldUpper);
#endif
#ifdef JJF_ZERO
    // always free up lambda
    for (int i = firstLambda_; i < firstLambda_ + 4; i++) {
        solver->setColLower(i, 0.0);
        solver->setColUpper(i, 2.0);
    }
#endif
    double xB[3];
    xB[0] = columnLower[xColumn_];
    xB[1] = columnUpper[xColumn_];
    double yB[3];
    yB[0] = columnLower[yColumn_];
    yB[1] = columnUpper[yColumn_];
    if (false && (branchingStrategy_&4) != 0 && yRow_ >= 0 &&
            xMeshSize_ == 1.0 && yMeshSize_ == 1.0) {
        if ((xB[1] - xB[0])*(yB[1] - yB[0]) < 40) {
            // try looking at all solutions
            double lower[4];
            double upper[4];
            double lambda[4];
            int i;
            double lowerLambda[4];
            double upperLambda[4];
            for (i = 0; i < 4; i++) {
                lower[i] = CoinMax(0.0, columnLower[firstLambda_+i]);
                upper[i] = CoinMin(1.0, columnUpper[firstLambda_+i]);
                lowerLambda[i] = 1.0;
                upperLambda[i] = 0.0;
            }
            // get coefficients
            double xybar[4];
            getCoefficients(solver, xB, yB, xybar);
            double x, y;
            for (x = xB[0]; x <= xB[1]; x++) {
                xB[2] = x;
                for (y = yB[0]; y <= yB[1]; y++) {
                    yB[2] = y;
                    computeLambdas(xB, yB, xybar, lambda);
                    for (i = 0; i < 4; i++) {
                        lowerLambda[i] = CoinMin(lowerLambda[i], lambda[i]);
                        upperLambda[i] = CoinMax(upperLambda[i], lambda[i]);
                    }
                }
            }
            double change = 0.0;;
            for (i = 0; i < 4; i++) {
                if (lowerLambda[i] > lower[i] + 1.0e-12) {
                    solver->setColLower(firstLambda_ + i, lowerLambda[i]);
                    change += lowerLambda[i] - lower[i];
                }
                if (upperLambda[i] < upper[i] - 1.0e-12) {
                    solver->setColUpper(firstLambda_ + i, upperLambda[i]);
                    change -= upperLambda[i] - upper[i];
                }
            }
            if (change > 1.0e-5)
                printf("change of %g\n", change);
        }
    }
    if (boundType_) {
        assert (!xMeshSize_ || !yMeshSize_);
        if (xMeshSize_) {
            // can tighten bounds on y
            if ((boundType_&1) != 0) {
                if (xB[0]*yB[1] > coefficient_) {
                    // tighten upper bound on y
                    solver->setColUpper(yColumn_, coefficient_ / xB[0]);
                }
            }
            if ((boundType_&2) != 0) {
                if (xB[1]*yB[0] < coefficient_) {
                    // tighten lower bound on y
                    solver->setColLower(yColumn_, coefficient_ / xB[1]);
                }
            }
        } else {
            // can tighten bounds on x
            if ((boundType_&1) != 0) {
                if (yB[0]*xB[1] > coefficient_) {
                    // tighten upper bound on x
                    solver->setColUpper(xColumn_, coefficient_ / yB[0]);
                }
            }
            if ((boundType_&2) != 0) {
                if (yB[1]*xB[0] < coefficient_) {
                    // tighten lower bound on x
                    solver->setColLower(xColumn_, coefficient_ / yB[1]);
                }
            }
        }
    }
}
// Compute lambdas if coefficients not changing
void
OsiBiLinear::computeLambdas(const OsiSolverInterface * solver, double lambda[4]) const
{
    // fix so correct
    double xB[3], yB[3];
    double xybar[4];
    getCoefficients(solver, xB, yB, xybar);
    double x, y;
    x = solver->getColLower()[xColumn_];
    assert(x == solver->getColUpper()[xColumn_]);
    xB[2] = x;
    y = solver->getColLower()[yColumn_];
    assert(y == solver->getColUpper()[yColumn_]);
    yB[2] = y;
    computeLambdas(xB, yB, xybar, lambda);
    assert (xyRow_ >= 0);
}
// Get LU coefficients from matrix
void
OsiBiLinear::getCoefficients(const OsiSolverInterface * solver, double xB[2], double yB[2],
                             double xybar[4]) const
{
    const CoinPackedMatrix * matrix = solver->getMatrixByCol();
    const double * element = matrix->getElements();
    const double * objective = solver->getObjCoefficients();
    const int * row = matrix->getIndices();
    const CoinBigIndex * columnStart = matrix->getVectorStarts();
    const int * columnLength = matrix->getVectorLengths();
    // order is LxLy, LxUy, UxLy and UxUy
    int j;
    double multiplier = (boundType_ == 0) ? 1.0 / coefficient_ : 1.0;
    if (yRow_ >= 0) {
        for (j = 0; j < 4; j++) {
            int iColumn = firstLambda_ + j;
            int iStart = columnStart[iColumn];
            int iEnd = iStart + columnLength[iColumn];
            int k = iStart;
            double x = 0.0;
            double y = 0.0;
            xybar[j] = 0.0;
            for (; k < iEnd; k++) {
                if (xRow_ == row[k])
                    x = element[k];
                if (yRow_ == row[k])
                    y = element[k];
                if (xyRow_ == row[k])
                    xybar[j] = element[k] * multiplier;
            }
            if (xyRow_ < 0)
                xybar[j] = objective[iColumn] * multiplier;
            if (j == 0)
                xB[0] = x;
            else if (j == 1)
                yB[1] = y;
            else if (j == 2)
                yB[0] = y;
            else if (j == 3)
                xB[1] = x;
            assert (fabs(xybar[j] - x*y) < 1.0e-4);
        }
    } else {
        // x==y
        for (j = 0; j < 4; j++) {
            int iColumn = firstLambda_ + j;
            int iStart = columnStart[iColumn];
            int iEnd = iStart + columnLength[iColumn];
            int k = iStart;
            double x = 0.0;
            xybar[j] = 0.0;
            for (; k < iEnd; k++) {
                if (xRow_ == row[k])
                    x = element[k];
                if (xyRow_ == row[k])
                    xybar[j] = element[k] * multiplier;
            }
            if (xyRow_ < 0)
                xybar[j] = objective[iColumn] * multiplier;
            if (j == 0) {
                xB[0] = x;
                yB[0] = x;
            } else if (j == 2) {
                xB[1] = x;
                yB[1] = x;
            }
        }
        assert (fabs(xybar[0] - xB[0]*yB[0]) < 1.0e-4);
        assert (fabs(xybar[1] - xB[0]*yB[1]) < 1.0e-4);
        assert (fabs(xybar[2] - xB[1]*yB[0]) < 1.0e-4);
        assert (fabs(xybar[3] - xB[1]*yB[1]) < 1.0e-4);
    }
}
// Compute lambdas (third entry in each .B is current value)
double
OsiBiLinear::computeLambdas(const double xB[3], const double yB[3], const double xybar[4], double lambda[4]) const
{
    // fake to make correct
    double x = xB[2];
    double y = yB[2];
    // order is LxLy, LxUy, UxLy and UxUy
    // l0 + l1 = this
    double rhs1 = (xB[1] - x) / (xB[1] - xB[0]);
    // l0 + l2 = this
    double rhs2 = (yB[1] - y) / (yB[1] - yB[0]);
    // For xy (taking out l3)
    double rhs3 = xB[1] * yB[1] - x * y;
    double a0 = xB[1] * yB[1] - xB[0] * yB[0];
    double a1 = xB[1] * yB[1] - xB[0] * yB[1];
    double a2 = xB[1] * yB[1] - xB[1] * yB[0];
    // divide through to get l0 coefficient
    rhs3 /= a0;
    a1 /= a0;
    a2 /= a0;
    // subtract out l0
    double b[2][2];
    double rhs[2];
    // first for l1 and l2
    b[0][0] = 1.0 - a1;
    b[0][1] = -a2;
    rhs[0] = rhs1 - rhs3;
    // second
    b[1][0] = -a1;
    b[1][1] = 1.0 - a2;
    rhs[1] = rhs2 - rhs3;
    if (fabs(b[0][0]) > fabs(b[0][1])) {
        double sub = b[1][0] / b[0][0];
        b[1][1] -= sub * b[0][1];
        rhs[1] -= sub * rhs[0];
        assert (fabs(b[1][1]) > 1.0e-12);
        lambda[2] = rhs[1] / b[1][1];
        lambda[0] = rhs2 - lambda[2];
        lambda[1] = rhs1 - lambda[0];
    } else {
        double sub = b[1][1] / b[0][1];
        b[1][0] -= sub * b[0][0];
        rhs[1] -= sub * rhs[0];
        assert (fabs(b[1][0]) > 1.0e-12);
        lambda[1] = rhs[1] / b[1][0];
        lambda[0] = rhs1 - lambda[1];
        lambda[2] = rhs2 - lambda[0];
    }
    lambda[3] = 1.0 - (lambda[0] + lambda[1] + lambda[2]);
    double infeasibility = 0.0;
    double xy = 0.0;
    for (int j = 0; j < 4; j++) {
        double value = lambda[j];
        if (value > 1.0) {
            infeasibility += value - 1.0;
            value = 1.0;
        }
        if (value < 0.0) {
            infeasibility -= value;
            value = 0.0;
        }
        lambda[j] = value;
        xy += xybar[j] * value;
    }
    assert (fabs(xy - x*y) < 1.0e-4);
    return infeasibility;
}
// Updates coefficients
int
OsiBiLinear::updateCoefficients(const double * lower, const double * upper, double * objective,
                                CoinPackedMatrix * matrix, CoinWarmStartBasis * basis) const
{
    // Return if no updates
    if ((branchingStrategy_&4) != 0)
        return 0;
    int numberUpdated = 0;
    double * element = matrix->getMutableElements();
    const int * row = matrix->getIndices();
    const CoinBigIndex * columnStart = matrix->getVectorStarts();
    const int * columnLength = matrix->getVectorLengths();
    // order is LxLy, LxUy, UxLy and UxUy
    double xB[2];
    double yB[2];
    xB[0] = lower[xColumn_];
    xB[1] = upper[xColumn_];
    yB[0] = lower[yColumn_];
    yB[1] = upper[yColumn_];
    //printf("x %d (%g,%g) y %d (%g,%g)\n",
    // xColumn_,xB[0],xB[1],
    // yColumn_,yB[0],yB[1]);
    CoinWarmStartBasis::Status status[4];
    int numStruct = basis ? basis->getNumStructural() - firstLambda_ : 0;
    double coefficient = (boundType_ == 0) ? coefficient_ : 1.0;
    for (int j = 0; j < 4; j++) {
        status[j] = (j < numStruct) ? basis->getStructStatus(j + firstLambda_) : CoinWarmStartBasis::atLowerBound;
        int iX = j >> 1;
        double x = xB[iX];
        int iY = j & 1;
        double y = yB[iY];
        CoinBigIndex k = columnStart[j+firstLambda_];
        CoinBigIndex last = k + columnLength[j+firstLambda_];
        double value;
        // xy
        value = coefficient * x * y;
        if (xyRow_ >= 0) {
            assert (row[k] == xyRow_);
#if BI_PRINT > 1
            printf("j %d xy (%d,%d) coeff from %g to %g\n", j, xColumn_, yColumn_, element[k], value);
#endif
            element[k++] = value;
        } else {
            // objective
            objective[j+firstLambda_] = value;
        }
        numberUpdated++;
        // convexity
        assert (row[k] == convexity_);
        k++;
        // x
        value = x;
#if BI_PRINT > 1
        printf("j %d x (%d) coeff from %g to %g\n", j, xColumn_, element[k], value);
#endif
        assert (row[k] == xRow_);
        element[k++] = value;
        numberUpdated++;
        if (yRow_ >= 0) {
            // y
            value = y;
#if BI_PRINT > 1
            printf("j %d y (%d) coeff from %g to %g\n", j, yColumn_, element[k], value);
#endif
            assert (row[k] == yRow_);
            element[k++] = value;
            numberUpdated++;
        }
        // Do extra rows
        for (int i = 0; i < numberExtraRows_; i++) {
            int iRow = extraRow_[i];
            for (; k < last; k++) {
                if (row[k] == iRow)
                    break;
            }
            assert (k < last);
            element[k++] = x * y * multiplier_[i];
        }
    }

    if (xB[0] == xB[1]) {
        if (yB[0] == yB[1]) {
            // only one basic
            bool first = true;
            for (int j = 0; j < 4; j++) {
                if (status[j] == CoinWarmStartBasis::basic) {
                    if (first) {
                        first = false;
                    } else {
                        basis->setStructStatus(j + firstLambda_, CoinWarmStartBasis::atLowerBound);
#if BI_PRINT
                        printf("zapping %d (x=%d,y=%d)\n", j, xColumn_, yColumn_);
#endif
                    }
                }
            }
        } else {
            if (status[0] == CoinWarmStartBasis::basic &&
                    status[2] == CoinWarmStartBasis::basic) {
                basis->setStructStatus(2 + firstLambda_, CoinWarmStartBasis::atLowerBound);
#if BI_PRINT
                printf("zapping %d (x=%d,y=%d)\n", 2, xColumn_, yColumn_);
#endif
            }
            if (status[1] == CoinWarmStartBasis::basic &&
                    status[3] == CoinWarmStartBasis::basic) {
                basis->setStructStatus(3 + firstLambda_, CoinWarmStartBasis::atLowerBound);
#if BI_PRINT
                printf("zapping %d (x=%d,y=%d)\n", 3, xColumn_, yColumn_);
#endif
            }
        }
    } else if (yB[0] == yB[1]) {
        if (status[0] == CoinWarmStartBasis::basic &&
                status[1] == CoinWarmStartBasis::basic) {
            basis->setStructStatus(1 + firstLambda_, CoinWarmStartBasis::atLowerBound);
#if BI_PRINT
            printf("zapping %d (x=%d,y=%d)\n", 1, xColumn_, yColumn_);
#endif
        }
        if (status[2] == CoinWarmStartBasis::basic &&
                status[3] == CoinWarmStartBasis::basic) {
            basis->setStructStatus(3 + firstLambda_, CoinWarmStartBasis::atLowerBound);
#if BI_PRINT
            printf("zapping %d (x=%d,y=%d)\n", 3, xColumn_, yColumn_);
#endif
        }
    }
    return numberUpdated;
}
// This does NOT set mutable stuff
double
OsiBiLinear::checkInfeasibility(const OsiBranchingInformation * info) const
{
    // If another object has finer mesh ignore this
    if ((branchingStrategy_&8) != 0)
        return 0.0;
    int way;
    double saveInfeasibility = infeasibility_;
    short int saveWhichWay = whichWay_;
    double saveXyBranchValue = xyBranchValue_;
    short saveChosen = chosen_;
    double value = infeasibility(info, way);
    infeasibility_ = saveInfeasibility;
    whichWay_ = saveWhichWay;
    xyBranchValue_ = saveXyBranchValue;
    chosen_ = saveChosen;
    return value;
}
OsiBiLinearBranchingObject::OsiBiLinearBranchingObject()
        : OsiTwoWayBranchingObject(),
        chosen_(0)
{
}

// Useful constructor
OsiBiLinearBranchingObject::OsiBiLinearBranchingObject (OsiSolverInterface * solver,
        const OsiBiLinear * set,
        int way ,
        double separator,
        int chosen)
        : OsiTwoWayBranchingObject(solver, set, way, separator),
        chosen_(static_cast<short int>(chosen))
{
    assert (chosen_ >= 0 && chosen_ < 2);
}

// Copy constructor
OsiBiLinearBranchingObject::OsiBiLinearBranchingObject ( const OsiBiLinearBranchingObject & rhs)
        : OsiTwoWayBranchingObject(rhs),
        chosen_(rhs.chosen_)
{
}

// Assignment operator
OsiBiLinearBranchingObject &
OsiBiLinearBranchingObject::operator=( const OsiBiLinearBranchingObject & rhs)
{
    if (this != &rhs) {
        OsiTwoWayBranchingObject::operator=(rhs);
        chosen_ = rhs.chosen_;
    }
    return *this;
}
OsiBranchingObject *
OsiBiLinearBranchingObject::clone() const
{
    return (new OsiBiLinearBranchingObject(*this));
}


// Destructor
OsiBiLinearBranchingObject::~OsiBiLinearBranchingObject ()
{
}
double
OsiBiLinearBranchingObject::branch(OsiSolverInterface * solver)
{
    const OsiBiLinear * set =
        dynamic_cast <const OsiBiLinear *>(originalObject_) ;
    assert (set);
    int way = (!branchIndex_) ? (2 * firstBranch_ - 1) : -(2 * firstBranch_ - 1);
    branchIndex_++;
    set->newBounds(solver, way, chosen_, value_);
    return 0.0;
}
/* Return true if branch should only bound variables
 */
bool
OsiBiLinearBranchingObject::boundBranch() const
{
    const OsiBiLinear * set =
        dynamic_cast <const OsiBiLinear *>(originalObject_) ;
    assert (set);
    return (set->branchingStrategy()&4) != 0;
}
// Print what would happen
void
OsiBiLinearBranchingObject::print(const OsiSolverInterface * /*solver*/)
{
    const OsiBiLinear * set =
        dynamic_cast <const OsiBiLinear *>(originalObject_) ;
    assert (set);
    int way = (!branchIndex_) ? (2 * firstBranch_ - 1) : -(2 * firstBranch_ - 1);
    int iColumn = (chosen_ == 1) ? set->xColumn() : set->yColumn();
    printf("OsiBiLinear would branch %s on %c variable %d from value %g\n",
           (way < 0) ? "down" : "up",
           (chosen_ == 0) ? 'X' : 'Y', iColumn, value_);
}
// Default Constructor
OsiBiLinearEquality::OsiBiLinearEquality ()
        : OsiBiLinear(),
        numberPoints_(0)
{
}

// Useful constructor
OsiBiLinearEquality::OsiBiLinearEquality (OsiSolverInterface * solver, int xColumn,
        int yColumn, int xyRow, double rhs,
        double xMesh)
        : OsiBiLinear(),
        numberPoints_(0)
{
    double xB[2];
    double yB[2];
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    xColumn_ = xColumn;
    yColumn_ = yColumn;
    xyRow_ = xyRow;
    coefficient_ = rhs;
    xB[0] = lower[xColumn_];
    xB[1] = upper[xColumn_];
    yB[0] = lower[yColumn_];
    yB[1] = upper[yColumn_];
    if (xB[1]*yB[1] < coefficient_ + 1.0e-12 ||
            xB[0]*yB[0] > coefficient_ - 1.0e-12) {
        printf("infeasible row - reformulate\n");
        abort();
    }
    // reduce range of x if possible
    if (yB[0]*xB[1] > coefficient_ + 1.0e12) {
        xB[1] = coefficient_ / yB[0];
        solver->setColUpper(xColumn_, xB[1]);
    }
    if (yB[1]*xB[0] < coefficient_ - 1.0e12) {
        xB[0] = coefficient_ / yB[1];
        solver->setColLower(xColumn_, xB[0]);
    }
    // See how many points
    numberPoints_ = static_cast<int> ((xB[1] - xB[0] + 0.5 * xMesh) / xMesh);
    // redo exactly
    xMeshSize_ = (xB[1] - xB[0]) / static_cast<double> (numberPoints_);
    numberPoints_++;
    //#define KEEPXY
#ifndef KEEPXY
    // Take out xyRow
    solver->setRowLower(xyRow_, 0.0);
    solver->setRowUpper(xyRow_, 0.0);
#else
    // make >=
    solver->setRowLower(xyRow_, coefficient_ - 0.05);
    solver->setRowUpper(xyRow_, COIN_DBL_MAX);
#endif
    double rowLower[3];
    double rowUpper[3];
#ifndef KEEPXY
    double * columnLower = new double [numberPoints_];
    double * columnUpper = new double [numberPoints_];
    double * objective = new double [numberPoints_];
    CoinBigIndex *starts = new CoinBigIndex[numberPoints_+1];
    int * index = new int[3*numberPoints_];
    double * element = new double [3*numberPoints_];
#else
    double * columnLower = new double [numberPoints_+2];
    double * columnUpper = new double [numberPoints_+2];
    double * objective = new double [numberPoints_+2];
    CoinBigIndex *starts = new CoinBigIndex[numberPoints_+3];
    int * index = new int[4*numberPoints_+2];
    double * element = new double [4*numberPoints_+2];
#endif
    int i;
    starts[0] = 0;
    // rows
    int numberRows = solver->getNumRows();
    // convexity
    rowLower[0] = 1.0;
    rowUpper[0] = 1.0;
    convexity_ = numberRows;
    starts[1] = 0;
    // x
    rowLower[1] = 0.0;
    rowUpper[1] = 0.0;
    index[0] = xColumn_;
    element[0] = -1.0;
    xRow_ = numberRows + 1;
    starts[2] = 1;
    rowLower[2] = 0.0;
    rowUpper[2] = 0.0;
    index[1] = yColumn;
    element[1] = -1.0;
    yRow_ = numberRows + 2;
    starts[3] = 2;
    solver->addRows(3, starts, index, element, rowLower, rowUpper);
    int n = 0;
    firstLambda_ = solver->getNumCols();
    double x = xB[0];
    assert(xColumn_ != yColumn_);
    for (i = 0; i < numberPoints_; i++) {
        double y = coefficient_ / x;
        columnLower[i] = 0.0;
        columnUpper[i] = 2.0;
        objective[i] = 0.0;
        double value;
#ifdef KEEPXY
        // xy
        value = coefficient_;
        element[n] = value;
        index[n++] = xyRow_;
#endif
        // convexity
        value = 1.0;
        element[n] = value;
        index[n++] = 0 + numberRows;
        // x
        value = x;
        if (fabs(value) < 1.0e-19)
            value = 1.0e-19;
        element[n] = value;
        index[n++] = 1 + numberRows;
        // y
        value = y;
        if (fabs(value) < 1.0e-19)
            value = 1.0e-19;
        element[n] = value;
        index[n++] = 2 + numberRows;
        starts[i+1] = n;
        x += xMeshSize_;
    }
#ifdef KEEPXY
    // costed slacks
    columnLower[numberPoints_] = 0.0;
    columnUpper[numberPoints_] = xMeshSize_;
    objective[numberPoints_] = 1.0e3;;
    // convexity
    element[n] = 1.0;
    index[n++] = 0 + numberRows;
    starts[numberPoints_+1] = n;
    columnLower[numberPoints_+1] = 0.0;
    columnUpper[numberPoints_+1] = xMeshSize_;
    objective[numberPoints_+1] = 1.0e3;;
    // convexity
    element[n] = -1.0;
    index[n++] = 0 + numberRows;
    starts[numberPoints_+2] = n;
    solver->addCols(numberPoints_ + 2, starts, index, element, columnLower, columnUpper, objective);
#else
    solver->addCols(numberPoints_, starts, index, element, columnLower, columnUpper, objective);
#endif
    delete [] columnLower;
    delete [] columnUpper;
    delete [] objective;
    delete [] starts;
    delete [] index;
    delete [] element;
}

// Copy constructor
OsiBiLinearEquality::OsiBiLinearEquality ( const OsiBiLinearEquality & rhs)
        : OsiBiLinear(rhs),
        numberPoints_(rhs.numberPoints_)
{
}

// Clone
OsiObject *
OsiBiLinearEquality::clone() const
{
    return new OsiBiLinearEquality(*this);
}

// Assignment operator
OsiBiLinearEquality &
OsiBiLinearEquality::operator=( const OsiBiLinearEquality & rhs)
{
    if (this != &rhs) {
        OsiBiLinear::operator=(rhs);
        numberPoints_ = rhs.numberPoints_;
    }
    return *this;
}

// Destructor
OsiBiLinearEquality::~OsiBiLinearEquality ()
{
}
// Possible improvement
double
OsiBiLinearEquality::improvement(const OsiSolverInterface * solver) const
{
    const double * pi = solver->getRowPrice();
    int i;
    const double * solution = solver->getColSolution();
    printf(" for x %d y %d - pi %g %g\n", xColumn_, yColumn_, pi[xRow_], pi[yRow_]);
    for (i = 0; i < numberPoints_; i++) {
        if (fabs(solution[i+firstLambda_]) > 1.0e-7)
            printf("(%d %g) ", i, solution[i+firstLambda_]);
    }
    printf("\n");
    return 0.0;
}
/* change grid
   if type 0 then use solution and make finer
   if 1 then back to original
*/
double
OsiBiLinearEquality::newGrid(OsiSolverInterface * solver, int type) const
{
    CoinPackedMatrix * matrix = solver->getMutableMatrixByCol();
    if (!matrix) {
        printf("Unable to modify matrix\n");
        abort();
    }
    double * element = matrix->getMutableElements();
#ifndef NDEBUG
    const int * row = matrix->getIndices();
#endif
    const CoinBigIndex * columnStart = matrix->getVectorStarts();
    //const int * columnLength = matrix->getVectorLengths();
    // get original bounds
    double xB[2];
    double yB[2];
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    xB[0] = lower[xColumn_];
    xB[1] = upper[xColumn_];
    assert (fabs((xB[1] - xB[0]) - xMeshSize_*(numberPoints_ - 1)) < 1.0e-7);
    yB[0] = lower[yColumn_];
    yB[1] = upper[yColumn_];
    double mesh = 0.0;
    int i;
    if (type == 0) {
        const double * solution = solver->getColSolution();
        int first = -1;
        int last = -1;
        double xValue = 0.0;
        double step = 0.0;
        for (i = 0; i < numberPoints_; i++) {
            int iColumn = i + firstLambda_;
            if (fabs(solution[iColumn]) > 1.0e-7) {
                int k = columnStart[iColumn] + 1;
                xValue += element[k] * solution[iColumn];
                if (first == -1) {
                    first = i;
                    step = -element[k];
                } else {
                    step += element[k];
                }
                last = i;
            }
        }
        if (last > first + 1) {
            printf("not adjacent - presuming small djs\n");
        }
        // new step size
        assert (numberPoints_ > 2);
        step = CoinMax((1.5 * step) / static_cast<double> (numberPoints_ - 1), 0.5 * step);
        xB[0] = CoinMax(xB[0], xValue - 0.5 * step);
        xB[1] = CoinMin(xB[1], xValue + 0.5 * step);
        // and now divide these
        mesh = (xB[1] - xB[0]) / static_cast<double> (numberPoints_ - 1);
    } else {
        // back to original
        mesh = xMeshSize_;
    }
    double x = xB[0];
    for (i = 0; i < numberPoints_; i++) {
        int iColumn = i + firstLambda_;
        double y = coefficient_ / x;
        //assert (columnLength[iColumn]==3); - could have cuts
        int k = columnStart[iColumn];
#ifdef KEEPXY
        // xy
        assert (row[k] == xyRow_);
        k++;
#endif
        assert (row[k] == convexity_);
        k++;
        double value;
        // x
        value = x;
        assert (row[k] == xRow_);
        assert (fabs(value) > 1.0e-10);
        element[k++] = value;
        // y
        value = y;
        assert (row[k] == yRow_);
        assert (fabs(value) > 1.0e-10);
        element[k++] = value;
        x += mesh;
    }
    return mesh;
}
/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
OsiSimpleFixedInteger::OsiSimpleFixedInteger ()
        : OsiSimpleInteger()
{
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
OsiSimpleFixedInteger::OsiSimpleFixedInteger (const OsiSolverInterface * solver, int iColumn)
        : OsiSimpleInteger(solver, iColumn)
{
}


// Useful constructor - passed solver index and original bounds
OsiSimpleFixedInteger::OsiSimpleFixedInteger ( int iColumn, double lower, double upper)
        : OsiSimpleInteger(iColumn, lower, upper)
{
}

// Useful constructor - passed simple integer
OsiSimpleFixedInteger::OsiSimpleFixedInteger ( const OsiSimpleInteger &rhs)
        : OsiSimpleInteger(rhs)
{
}

// Copy constructor
OsiSimpleFixedInteger::OsiSimpleFixedInteger ( const OsiSimpleFixedInteger & rhs)
        : OsiSimpleInteger(rhs)

{
}

// Clone
OsiObject *
OsiSimpleFixedInteger::clone() const
{
    return new OsiSimpleFixedInteger(*this);
}

// Assignment operator
OsiSimpleFixedInteger &
OsiSimpleFixedInteger::operator=( const OsiSimpleFixedInteger & rhs)
{
    if (this != &rhs) {
        OsiSimpleInteger::operator=(rhs);
    }
    return *this;
}

// Destructor
OsiSimpleFixedInteger::~OsiSimpleFixedInteger ()
{
}
// Infeasibility - large is 0.5
double
OsiSimpleFixedInteger::infeasibility(const OsiBranchingInformation * info, int & whichWay) const
{
    double value = info->solution_[columnNumber_];
    value = CoinMax(value, info->lower_[columnNumber_]);
    value = CoinMin(value, info->upper_[columnNumber_]);
    double nearest = floor(value + (1.0 - 0.5));
    if (nearest > value) {
        whichWay = 1;
    } else {
        whichWay = 0;
    }
    infeasibility_ = fabs(value - nearest);
    bool satisfied = false;
    if (infeasibility_ <= info->integerTolerance_) {
        otherInfeasibility_ = 1.0;
        satisfied = true;
        if (info->lower_[columnNumber_] != info->upper_[columnNumber_])
            infeasibility_ = 1.0e-5;
        else
            infeasibility_ = 0.0;
    } else if (info->defaultDual_ < 0.0) {
        otherInfeasibility_ = 1.0 - infeasibility_;
    } else {
        const double * pi = info->pi_;
        const double * activity = info->rowActivity_;
        const double * lower = info->rowLower_;
        const double * upper = info->rowUpper_;
        const double * element = info->elementByColumn_;
        const int * row = info->row_;
        const CoinBigIndex * columnStart = info->columnStart_;
        const int * columnLength = info->columnLength_;
        double direction = info->direction_;
        double downMovement = value - floor(value);
        double upMovement = 1.0 - downMovement;
        double valueP = info->objective_[columnNumber_] * direction;
        CoinBigIndex start = columnStart[columnNumber_];
        CoinBigIndex end = start + columnLength[columnNumber_];
        double upEstimate = 0.0;
        double downEstimate = 0.0;
        if (valueP > 0.0)
            upEstimate = valueP * upMovement;
        else
            downEstimate -= valueP * downMovement;
        double tolerance = info->primalTolerance_;
        for (CoinBigIndex j = start; j < end; j++) {
            int iRow = row[j];
            if (lower[iRow] < -1.0e20)
                assert (pi[iRow] <= 1.0e-3);
            if (upper[iRow] > 1.0e20)
                assert (pi[iRow] >= -1.0e-3);
            valueP = pi[iRow] * direction;
            double el2 = element[j];
            double value2 = valueP * el2;
            double u = 0.0;
            double d = 0.0;
            if (value2 > 0.0)
                u = value2;
            else
                d = -value2;
            // if up makes infeasible then make at least default
            double newUp = activity[iRow] + upMovement * el2;
            if (newUp > upper[iRow] + tolerance || newUp < lower[iRow] - tolerance)
                u = CoinMax(u, info->defaultDual_);
            upEstimate += u * upMovement * fabs(el2);
            // if down makes infeasible then make at least default
            double newDown = activity[iRow] - downMovement * el2;
            if (newDown > upper[iRow] + tolerance || newDown < lower[iRow] - tolerance)
                d = CoinMax(d, info->defaultDual_);
            downEstimate += d * downMovement * fabs(el2);
        }
        if (downEstimate >= upEstimate) {
            infeasibility_ = CoinMax(1.0e-12, upEstimate);
            otherInfeasibility_ = CoinMax(1.0e-12, downEstimate);
            whichWay = 1;
        } else {
            infeasibility_ = CoinMax(1.0e-12, downEstimate);
            otherInfeasibility_ = CoinMax(1.0e-12, upEstimate);
            whichWay = 0;
        }
    }
    if (preferredWay_ >= 0 && !satisfied)
        whichWay = preferredWay_;
    whichWay_ = static_cast<short int>(whichWay);
    return infeasibility_;
}
// Creates a branching object
OsiBranchingObject *
OsiSimpleFixedInteger::createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const
{
    double value = info->solution_[columnNumber_];
    value = CoinMax(value, info->lower_[columnNumber_]);
    value = CoinMin(value, info->upper_[columnNumber_]);
    assert (info->upper_[columnNumber_] > info->lower_[columnNumber_]);
    double nearest = floor(value + 0.5);
    double integerTolerance = info->integerTolerance_;
    if (fabs(value - nearest) < integerTolerance) {
        // adjust value
        if (nearest != info->upper_[columnNumber_])
            value = nearest + 2.0 * integerTolerance;
        else
            value = nearest - 2.0 * integerTolerance;
    }
    OsiBranchingObject * branch = new OsiIntegerBranchingObject(solver, this, way,
            value);
    return branch;
}

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>
//#define CGL_DEBUG 2
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CoinWarmStartBasis.hpp"
//#include "CglTemporary.hpp"
#include "CoinFinite.hpp"
//-------------------------------------------------------------------
// Generate Stored cuts
//-------------------------------------------------------------------
void
CglTemporary::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
                           const CglTreeInfo /*info*/) const
{
    // Get basic problem information
    const double * solution = si.getColSolution();
    int numberRowCuts = cuts_.sizeRowCuts();
    for (int i = 0; i < numberRowCuts; i++) {
        const OsiRowCut * rowCutPointer = cuts_.rowCutPtr(i);
        double violation = rowCutPointer->violated(solution);
        if (violation >= requiredViolation_)
            cs.insert(*rowCutPointer);
    }
    // delete
    cuts_ = OsiCuts();
}

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CglTemporary::CglTemporary ()
        :
        CglStored()
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CglTemporary::CglTemporary (const CglTemporary & source) :
        CglStored(source)
{
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglTemporary::clone() const
{
    return new CglTemporary(*this);
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CglTemporary::~CglTemporary ()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CglTemporary &
CglTemporary::operator=(const CglTemporary & rhs)
{
    if (this != &rhs) {
        CglStored::operator=(rhs);
    }
    return *this;
}
void checkQP(ClpSimplex * /*model*/)
{
#ifdef JJF_ZERO
    printf("Checking quadratic model %x\n", model);
    if (model) {
        ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(model->objectiveAsObject()));
        assert (quadraticObj);
        CoinPackedMatrix * quadraticObjective = quadraticObj->quadraticObjective();
        int numberColumns = quadraticObj->numberColumns();
        const int * columnQuadratic = quadraticObjective->getIndices();
        const CoinBigIndex * columnQuadraticStart = quadraticObjective->getVectorStarts();
        const int * columnQuadraticLength = quadraticObjective->getVectorLengths();
        //const double * quadraticElement = quadraticObjective->getElements();
        for (int i = 0; i < numberColumns; i++) {
            for (int j = columnQuadraticStart[i]; j < columnQuadraticStart[i] + columnQuadraticLength[i]; j++)
                assert (columnQuadratic[j] >= 0 && columnQuadratic[j] < 1000);
        }
    }
#endif
}
//#############################################################################
// Solve methods
//#############################################################################
void OsiSolverLinearizedQuadratic::initialSolve()
{
    OsiClpSolverInterface::initialSolve();
    int secondaryStatus = modelPtr_->secondaryStatus();
    if (modelPtr_->status() == 0 && (secondaryStatus == 2 || secondaryStatus == 4))
        modelPtr_->cleanup(1);
    if (isProvenOptimal() && modelPtr_->numberColumns() == quadraticModel_->numberColumns()) {
        // see if qp can get better solution
        const double * solution = modelPtr_->primalColumnSolution();
        int numberColumns = modelPtr_->numberColumns();
        bool satisfied = true;
        for (int i = 0; i < numberColumns; i++) {
            if (isInteger(i)) {
                double value = solution[i];
                if (fabs(value - floor(value + 0.5)) > 1.0e-6) {
                    satisfied = false;
                    break;
                }
            }
        }
        if (satisfied) {
            checkQP(quadraticModel_);
            ClpSimplex qpTemp(*quadraticModel_);
            checkQP(&qpTemp);
            double * lower = qpTemp.columnLower();
            double * upper = qpTemp.columnUpper();
            double * lower2 = modelPtr_->columnLower();
            double * upper2 = modelPtr_->columnUpper();
            for (int i = 0; i < numberColumns; i++) {
                if (isInteger(i)) {
                    double value = floor(solution[i] + 0.5);
                    lower[i] = value;
                    upper[i] = value;
                } else {
                    lower[i] = lower2[i];
                    upper[i] = upper2[i];
                }
            }
            //qpTemp.writeMps("bad.mps");
            //modelPtr_->writeMps("bad2.mps");
            //qpTemp.objectiveAsObject()->setActivated(0);
            //qpTemp.primal();
            //qpTemp.objectiveAsObject()->setActivated(1);
            qpTemp.primal();
            //assert (!qpTemp.problemStatus());
            if (qpTemp.objectiveValue() < bestObjectiveValue_ && !qpTemp.problemStatus()) {
                delete [] bestSolution_;
                bestSolution_ = CoinCopyOfArray(qpTemp.primalColumnSolution(), numberColumns);
                bestObjectiveValue_ = qpTemp.objectiveValue();
                printf("better qp objective of %g\n", bestObjectiveValue_);
            }
        }
    }
}
//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################
//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
OsiSolverLinearizedQuadratic::OsiSolverLinearizedQuadratic ()
        : OsiClpSolverInterface()
{
    bestObjectiveValue_ = COIN_DBL_MAX;
    bestSolution_ = NULL;
    specialOptions3_ = 0;
    quadraticModel_ = NULL;
}
OsiSolverLinearizedQuadratic::OsiSolverLinearizedQuadratic ( ClpSimplex * quadraticModel)
        : OsiClpSolverInterface(new ClpSimplex(*quadraticModel), true)
{
    bestObjectiveValue_ = COIN_DBL_MAX;
    bestSolution_ = NULL;
    specialOptions3_ = 0;
    quadraticModel_ = new ClpSimplex(*quadraticModel);
    // linearize
    int numberColumns = modelPtr_->numberColumns();
    const double * solution = modelPtr_->primalColumnSolution();
    // Replace objective
    ClpObjective * trueObjective = modelPtr_->objectiveAsObject();
    ClpObjective * objective = new ClpLinearObjective(NULL, numberColumns);
    modelPtr_->setObjectivePointer(objective);
    double offset;
    double saveOffset = modelPtr_->objectiveOffset();
    memcpy(modelPtr_->objective(), trueObjective->gradient(modelPtr_, solution, offset, true, 2),
           numberColumns*sizeof(double));
    modelPtr_->setObjectiveOffset(saveOffset + offset);
    delete trueObjective;
    checkQP(quadraticModel_);
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface *
OsiSolverLinearizedQuadratic::clone(bool /*copyData*/) const
{
    //assert (copyData);
    return new OsiSolverLinearizedQuadratic(*this);
}


//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
OsiSolverLinearizedQuadratic::OsiSolverLinearizedQuadratic (
    const OsiSolverLinearizedQuadratic & rhs)
        : OsiSolverInterface(rhs)
        , OsiClpSolverInterface(rhs)
{
    bestObjectiveValue_ = rhs.bestObjectiveValue_;
    if (rhs.bestSolution_) {
        bestSolution_ = CoinCopyOfArray(rhs.bestSolution_, modelPtr_->numberColumns());
    } else {
        bestSolution_ = NULL;
    }
    specialOptions3_ = rhs.specialOptions3_;
    if (rhs.quadraticModel_) {
        quadraticModel_ = new ClpSimplex(*rhs.quadraticModel_);
    } else {
        quadraticModel_ = NULL;
    }
    checkQP(rhs.quadraticModel_);
    checkQP(quadraticModel_);
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
OsiSolverLinearizedQuadratic::~OsiSolverLinearizedQuadratic ()
{
    delete [] bestSolution_;
    delete quadraticModel_;
}

//-------------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
OsiSolverLinearizedQuadratic &
OsiSolverLinearizedQuadratic::operator=(const OsiSolverLinearizedQuadratic & rhs)
{
    if (this != &rhs) {
        delete [] bestSolution_;
        delete quadraticModel_;
        OsiClpSolverInterface::operator=(rhs);
        bestObjectiveValue_ = rhs.bestObjectiveValue_;
        if (rhs.bestSolution_) {
            bestSolution_ = CoinCopyOfArray(rhs.bestSolution_, modelPtr_->numberColumns());
        } else {
            bestSolution_ = NULL;
        }
        specialOptions3_ = rhs.specialOptions3_;
        if (rhs.quadraticModel_) {
            quadraticModel_ = new ClpSimplex(*rhs.quadraticModel_);
        } else {
            quadraticModel_ = NULL;
        }
        checkQP(rhs.quadraticModel_);
        checkQP(quadraticModel_);
    }
    return *this;
}
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
int
CoinModel::expandKnapsack(int knapsackRow, int & numberOutput, double * buildObj, CoinBigIndex * buildStart,
                          int * buildRow, double * buildElement, int reConstruct) const
{
    /* mark rows
       -2 in knapsack and other variables
       -1 not involved
       0 only in knapsack
    */
    int * markRow = new int [numberRows_];
    int iRow;
    int iColumn;
    int * whichColumn = new int [numberColumns_];
    for (iColumn = 0; iColumn < numberColumns_; iColumn++)
        whichColumn[iColumn] = -1;
    int numJ = 0;
    for (iRow = 0; iRow < numberRows_; iRow++)
        markRow[iRow] = -1;
    CoinModelLink triple;
    triple = firstInRow(knapsackRow);
    while (triple.column() >= 0) {
        int iColumn = triple.column();
#ifndef NDEBUG
        const char *  el = getElementAsString(knapsackRow, iColumn);
        assert (!strcmp("Numeric", el));
#endif
        whichColumn[iColumn] = numJ;
        numJ++;
        triple = next(triple);
    }
    for (iRow = 0; iRow < numberRows_; iRow++) {
        triple = firstInRow(iRow);
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
            triple = next(triple);
        }
        if (type == -3)
            type = -1;
        markRow[iRow] = type;
    }
    int * bound = new int [numberColumns_+1];
    int * whichRow = new int [numberRows_];
    ClpSimplex tempModel;
    CoinModel tempModel2(*this);
    tempModel.loadProblem(tempModel2);
    int * stack = new int [numberColumns_+1];
    double * size = new double [numberColumns_+1];
    double * rhsOffset = new double[numberRows_];
    int * build = new int[numberColumns_];
    int maxNumber = numberOutput;
    numJ = 0;
    double minSize = getRowLower(knapsackRow);
    double maxSize = getRowUpper(knapsackRow);
    double offset = 0.0;
    triple = firstInRow(knapsackRow);
    while (triple.column() >= 0) {
        iColumn = triple.column();
        double lowerColumn = columnLower(iColumn);
        double upperColumn = columnUpper(iColumn);
        double gap = upperColumn - lowerColumn;
        if (gap > 1.0e8)
            gap = 1.0e8;
        assert (fabs(floor(gap + 0.5) - gap) < 1.0e-5);
        whichColumn[numJ] = iColumn;
        bound[numJ] = static_cast<int> (gap);
        size[numJ++] = triple.value();
        offset += triple.value() * lowerColumn;
        triple = next(triple);
    }
    int jRow;
    for (iRow = 0; iRow < numberRows_; iRow++)
        whichRow[iRow] = iRow;
    ClpSimplex smallModel(&tempModel, numberRows_, whichRow, numJ, whichColumn, true, true, true);
    // modify rhs to allow for nonzero lower bounds
    double * rowLower = smallModel.rowLower();
    double * rowUpper = smallModel.rowUpper();
    const double * columnLower = smallModel.columnLower();
    //const double * columnUpper = smallModel.columnUpper();
    const CoinPackedMatrix * matrix = smallModel.matrix();
    const double * element = matrix->getElements();
    const int * row = matrix->getIndices();
    const CoinBigIndex * columnStart = matrix->getVectorStarts();
    const int * columnLength = matrix->getVectorLengths();
    const double * objective = smallModel.objective();
    double objectiveOffset = 0.0;
    CoinZeroN(rhsOffset, numberRows_);
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
    for (jRow = 0; jRow < numberRows_; jRow++) {
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
    double * rowActivity = smallModel.primalRowSolution();
    CoinZeroN(rowActivity, numberRows_);
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
        assert (buildRow && buildElement);
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
                    buildStart[numberOutput+1] = nelCreate;
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
                assert (fabs(sum - checkSize) < 1.0e-3);
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
                stack[iStack] ++;
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
    delete [] whichColumn;
    delete [] whichRow;
    delete [] bound;
    delete [] stack;
    delete [] size;
    delete [] rhsOffset;
    delete [] build;
    delete [] markRow;
    return nelCreate;
}
#include "ClpConstraint.hpp"
#include "ClpConstraintLinear.hpp"
#include "ClpConstraintQuadratic.hpp"
#ifdef COIN_HAS_ASL
//#include "ClpAmplObjective.hpp"
#endif
/* Return an approximate solution to a CoinModel.
    Lots of bounds may be odd to force a solution.
    mode = 0 just tries to get a continuous solution
*/
ClpSimplex *
approximateSolution(CoinModel & coinModel,
                    int numberPasses, double deltaTolerance,
                    int /*mode*/)
{
#ifndef JJF_ONE
    //#ifdef COIN_HAS_ASL
    // matrix etc will be changed
    CoinModel coinModel2 = coinModel;
    if (coinModel2.moreInfo()) {
        // for now just ampl objective
        ClpSimplex * model = new ClpSimplex();
        model->loadProblem(coinModel2);
        int numberConstraints;
        ClpConstraint ** constraints = NULL;
        int type = model->loadNonLinear(coinModel2.moreInfo(),
                                        numberConstraints, constraints);
        if (type == 1 || type == 3) {
            model->nonlinearSLP(numberPasses, deltaTolerance);
        } else if (type == 2 || type == 4) {
            model->nonlinearSLP(numberConstraints, constraints,
                                numberPasses, deltaTolerance);
        } else {
            printf("error or linear - fix %d\n", type);
        }
        //exit(66);
        return model;
    }
    // first check and set up arrays
    int numberColumns = coinModel.numberColumns();
    int numberRows = coinModel.numberRows();
    // List of nonlinear rows
    int * which = new int[numberRows];
    bool testLinear = false;
    int numberConstraints = 0;
    int iColumn;
    bool linearObjective = true;
    int maximumQuadraticElements = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        // See if quadratic objective
        const char * expr = coinModel.getColumnObjectiveAsString(iColumn);
        if (strcmp(expr, "Numeric")) {
            linearObjective = false;
            // check if value*x+-value*y....
            assert (strlen(expr) < 20000);
            char temp[20000];
            strcpy(temp, expr);
            char * pos = temp;
            bool ifFirst = true;
            double linearTerm = 0.0;
            while (*pos) {
                double value;
                int jColumn = decodeBit(pos, pos, value, ifFirst, coinModel);
                // must be column unless first when may be linear term
                if (jColumn >= 0) {
                    maximumQuadraticElements++;
                } else if (jColumn == -2) {
                    linearTerm = value;
                } else {
                    printf("bad nonlinear term %s\n", temp);
                    abort();
                }
                ifFirst = false;
            }
        }
    }
    if (!linearObjective) {
        // zero objective
        for (iColumn = 0; iColumn < numberColumns; iColumn++)
            coinModel2.setObjective(iColumn, 0.0);
    }
    int iRow;
    for (iRow = 0; iRow < numberRows; iRow++) {
        int numberQuadratic = 0;
        bool linear = true;
        CoinModelLink triple = coinModel.firstInRow(iRow);
        while (triple.column() >= 0) {
            int iColumn = triple.column();
            const char *  expr = coinModel.getElementAsString(iRow, iColumn);
            if (strcmp("Numeric", expr)) {
                linear = false;
                // check if value*x+-value*y....
                assert (strlen(expr) < 20000);
                char temp[20000];
                strcpy(temp, expr);
                char * pos = temp;
                bool ifFirst = true;
                double linearTerm = 0.0;
                while (*pos) {
                    double value;
                    int jColumn = decodeBit(pos, pos, value, ifFirst, coinModel);
                    // must be column unless first when may be linear term
                    if (jColumn >= 0) {
                        numberQuadratic++;
                    } else if (jColumn == -2) {
                        linearTerm = value;
                    } else {
                        printf("bad nonlinear term %s\n", temp);
                        abort();
                    }
                    ifFirst = false;
                }
            }
            triple = coinModel.next(triple);
        }
        if (!linear || testLinear) {
            CoinModelLink triple = coinModel.firstInRow(iRow);
            while (triple.column() >= 0) {
                int iColumn = triple.column();
                coinModel2.setElement(iRow, iColumn, 0.0);
                triple = coinModel.next(triple);
            }
            which[numberConstraints++] = iRow;
            maximumQuadraticElements = CoinMax(maximumQuadraticElements, numberQuadratic);
        }
    }
    ClpSimplex * model = new ClpSimplex();
    // return if nothing
    if (!numberConstraints && linearObjective) {
        delete [] which;
        model->loadProblem(coinModel);
        model->dual();
        return model;
    }
    // space for quadratic
    // allow for linear term
    maximumQuadraticElements += numberColumns;
    CoinBigIndex * startQuadratic = new CoinBigIndex [numberColumns+1];
    int * columnQuadratic = new int [maximumQuadraticElements];
    double * elementQuadratic = new double [maximumQuadraticElements];
    ClpConstraint ** constraints = new ClpConstraint * [numberConstraints];
    double * linearTerm = new double [numberColumns];
    int saveNumber = numberConstraints;
    numberConstraints = 0;
    ClpQuadraticObjective * quadObj = NULL;
    if (!linearObjective) {
        int numberQuadratic = 0;
        CoinZeroN(linearTerm, numberColumns);
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            startQuadratic[iColumn] = numberQuadratic;
            // See if quadratic objective
            const char * expr = coinModel.getColumnObjectiveAsString(iColumn);
            if (strcmp(expr, "Numeric")) {
                // value*x*y
                char temp[20000];
                strcpy(temp, expr);
                char * pos = temp;
                bool ifFirst = true;
                while (*pos) {
                    double value;
                    int jColumn = decodeBit(pos, pos, value, ifFirst, coinModel);
                    // must be column unless first when may be linear term
                    if (jColumn >= 0) {
                        columnQuadratic[numberQuadratic] = jColumn;
                        if (jColumn != iColumn)
                            elementQuadratic[numberQuadratic++] = 1.0 * value; // convention
                        else if (jColumn == iColumn)
                            elementQuadratic[numberQuadratic++] = 2.0 * value; // convention
                    } else if (jColumn == -2) {
                        linearTerm[iColumn] = value;
                    } else {
                        printf("bad nonlinear term %s\n", temp);
                        abort();
                    }
                    ifFirst = false;
                }
            } else {
                // linear part
                linearTerm[iColumn] = coinModel.getColumnObjective(iColumn);
            }
        }
        startQuadratic[numberColumns] = numberQuadratic;
        quadObj = new ClpQuadraticObjective(linearTerm, numberColumns,
                                            startQuadratic, columnQuadratic, elementQuadratic);
    }
    int iConstraint;
    for (iConstraint = 0; iConstraint < saveNumber; iConstraint++) {
        iRow = which[iConstraint];
        if (iRow >= 0) {
            int numberQuadratic = 0;
            int lastColumn = -1;
            int largestColumn = -1;
            CoinZeroN(linearTerm, numberColumns);
            CoinModelLink triple = coinModel.firstInRow(iRow);
            while (triple.column() >= 0) {
                int iColumn = triple.column();
                while (lastColumn < iColumn) {
                    startQuadratic[lastColumn+1] = numberQuadratic;
                    lastColumn++;
                }
                const char *  expr = coinModel.getElementAsString(iRow, iColumn);
                if (strcmp("Numeric", expr)) {
                    largestColumn = CoinMax(largestColumn, iColumn);
                    // value*x*y
                    char temp[20000];
                    strcpy(temp, expr);
                    char * pos = temp;
                    bool ifFirst = true;
                    while (*pos) {
                        double value;
                        int jColumn = decodeBit(pos, pos, value, ifFirst, coinModel);
                        // must be column unless first when may be linear term
                        if (jColumn >= 0) {
                            columnQuadratic[numberQuadratic] = jColumn;
                            if (jColumn == iColumn)
                                elementQuadratic[numberQuadratic++] = 2.0 * value; // convention
                            else
                                elementQuadratic[numberQuadratic++] = 1.0 * value; // convention
                            largestColumn = CoinMax(largestColumn, jColumn);
                        } else if (jColumn == -2) {
                            linearTerm[iColumn] = value;
                            // and put in as row -1
                            columnQuadratic[numberQuadratic] = -1;
                            if (jColumn == iColumn)
                                elementQuadratic[numberQuadratic++] = 2.0 * value; // convention
                            else
                                elementQuadratic[numberQuadratic++] = 1.0 * value; // convention
                            largestColumn = CoinMax(largestColumn, iColumn);
                        } else {
                            printf("bad nonlinear term %s\n", temp);
                            abort();
                        }
                        ifFirst = false;
                    }
                } else {
                    // linear part
                    linearTerm[iColumn] = coinModel.getElement(iRow, iColumn);
                    // and put in as row -1
                    columnQuadratic[numberQuadratic] = -1;
                    elementQuadratic[numberQuadratic++] = linearTerm[iColumn];
                    if (linearTerm[iColumn])
                        largestColumn = CoinMax(largestColumn, iColumn);
                }
                triple = coinModel.next(triple);
            }
            while (lastColumn < numberColumns) {
                startQuadratic[lastColumn+1] = numberQuadratic;
                lastColumn++;
            }
            // here we create ClpConstraint
            if (testLinear) {
                int n = 0;
                int * indices = new int[numberColumns];
                for (int j = 0; j < numberColumns; j++) {
                    if (linearTerm[j]) {
                        linearTerm[n] = linearTerm[j];
                        indices[n++] = j;
                    }
                }
                /// Constructor from constraint
                constraints[numberConstraints++] = new ClpConstraintLinear(iRow, n, numberColumns,
                        indices, linearTerm);
                delete [] indices;
            } else {
                constraints[numberConstraints++] = new ClpConstraintQuadratic(iRow, largestColumn + 1, numberColumns,
                        startQuadratic, columnQuadratic, elementQuadratic);
            }
        }
    }
    delete [] startQuadratic;
    delete [] columnQuadratic;
    delete [] elementQuadratic;
    delete [] linearTerm;
    delete [] which;
    model->loadProblem(coinModel2);
    if (quadObj)
        model->setObjective(quadObj);
    delete quadObj;
    int returnCode;
    if (numberConstraints) {
        returnCode = model->nonlinearSLP(numberConstraints, constraints,
                                         numberPasses, deltaTolerance);
        for (iConstraint = 0; iConstraint < saveNumber; iConstraint++)
            delete constraints[iConstraint];
    } else {
        returnCode = model->nonlinearSLP(numberPasses, deltaTolerance);
    }
    delete [] constraints;
    assert (!returnCode);
    return model;
#else
    printf("loadNonLinear needs ampl\n");
    abort();
    return NULL;
#endif
}
OsiChooseStrongSubset::OsiChooseStrongSubset() :
        OsiChooseStrong(),
        numberObjectsToUse_(0)
{
}

OsiChooseStrongSubset::OsiChooseStrongSubset(const OsiSolverInterface * solver) :
        OsiChooseStrong(solver),
        numberObjectsToUse_(-1)
{
}

OsiChooseStrongSubset::OsiChooseStrongSubset(const OsiChooseStrongSubset & rhs)
        : OsiChooseStrong(rhs)
{
    numberObjectsToUse_ = -1;
}

OsiChooseStrongSubset &
OsiChooseStrongSubset::operator=(const OsiChooseStrongSubset & rhs)
{
    if (this != &rhs) {
        OsiChooseStrong::operator=(rhs);
        numberObjectsToUse_ = -1;
    }
    return *this;
}


OsiChooseStrongSubset::~OsiChooseStrongSubset ()
{
}

// Clone
OsiChooseVariable *
OsiChooseStrongSubset::clone() const
{
    return new OsiChooseStrongSubset(*this);
}
// Initialize
int
OsiChooseStrongSubset::setupList ( OsiBranchingInformation *info, bool initialize)
{
    assert (solver_ == info->solver_);
    // Only has to work with Clp
    OsiSolverInterface * solverA = const_cast<OsiSolverInterface *> (solver_);
    OsiSolverLink * solver = dynamic_cast<OsiSolverLink *> (solverA);
    assert (solver);
    int numberObjects = solver->numberObjects();
    if (numberObjects > pseudoCosts_.numberObjects()) {
        // redo useful arrays
        pseudoCosts_.initialize(numberObjects);
    }
    int numObj = numberObjects;
    if (numberObjectsToUse_ < 0) {
        // Sort objects so bilinear at end
        OsiObject ** sorted = new OsiObject * [numberObjects];
        OsiObject ** objects = solver->objects();
        numObj = 0;
        int numberBiLinear = 0;
        int i;
        for (i = 0; i < numberObjects; i++) {
            OsiObject * obj = objects[i];
            OsiBiLinear * objB = dynamic_cast<OsiBiLinear *> (obj);
            if (!objB)
                objects[numObj++] = obj;
            else
                sorted[numberBiLinear++] = obj;
        }
        numberObjectsToUse_ = numObj;
        for (i = 0; i < numberBiLinear; i++)
            objects[numObj++] = sorted[i];
        delete [] sorted;
        // See if any master objects
        for (i = 0; i < numberObjectsToUse_; i++) {
            OsiUsesBiLinear * obj = dynamic_cast<OsiUsesBiLinear *> (objects[i]);
            if (obj)
                obj->addBiLinearObjects(solver);
        }
    }
    solver->setNumberObjects(numberObjectsToUse_);
    numObj = numberObjectsToUse_;
    // Use shadow prices
    //info->defaultDual_=0.0;
    int numberUnsatisfied = OsiChooseStrong::setupList ( info, initialize);
    solver->setNumberObjects(numberObjects);
    numObj = numberObjects;
    return numberUnsatisfied;
}
/* Choose a variable
   Returns -
   -1 Node is infeasible
   0  Normal termination - we have a candidate
   1  All looks satisfied - no candidate
   2  We can change the bound on a variable - but we also have a strong branching candidate
   3  We can change the bound on a variable - but we have a non-strong branching candidate
   4  We can change the bound on a variable - no other candidates
   We can pick up branch from whichObject() and whichWay()
   We can pick up a forced branch (can change bound) from whichForcedObject() and whichForcedWay()
   If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
*/
int
OsiChooseStrongSubset::chooseVariable( OsiSolverInterface * solver, OsiBranchingInformation *info, bool fixVariables)
{
    //int numberObjects = solver->numberObjects();
    //solver->setNumberObjects(numberObjectsToUse_);
    //numberObjects_=numberObjectsToUse_;
    // Use shadow prices
    //info->defaultDual_=0.0;
    int returnCode = OsiChooseStrong::chooseVariable(solver, info, fixVariables);
    //solver->setNumberObjects(numberObjects);
    //numberObjects_=numberObjects;
    return returnCode;
}
/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
OsiUsesBiLinear::OsiUsesBiLinear ()
        : OsiSimpleInteger(),
        numberBiLinear_(0),
        type_(0),
        objects_(NULL)
{
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
OsiUsesBiLinear::OsiUsesBiLinear (const OsiSolverInterface * solver, int iColumn, int type)
        : OsiSimpleInteger(solver, iColumn),
        numberBiLinear_(0),
        type_(type),
        objects_(NULL)
{
    if (type_) {
        assert(originalLower_ == floor(originalLower_ + 0.5));
        assert(originalUpper_ == floor(originalUpper_ + 0.5));
    }
}


// Useful constructor - passed solver index and original bounds
OsiUsesBiLinear::OsiUsesBiLinear ( int iColumn, double lower, double upper, int type)
        : OsiSimpleInteger(iColumn, lower, upper),
        numberBiLinear_(0),
        type_(type),
        objects_(NULL)
{
    if (type_) {
        assert(originalLower_ == floor(originalLower_ + 0.5));
        assert(originalUpper_ == floor(originalUpper_ + 0.5));
    }
}

// Useful constructor - passed simple integer
OsiUsesBiLinear::OsiUsesBiLinear ( const OsiSimpleInteger &rhs, int type)
        : OsiSimpleInteger(rhs),
        numberBiLinear_(0),
        type_(type),
        objects_(NULL)
{
    if (type_) {
        assert(originalLower_ == floor(originalLower_ + 0.5));
        assert(originalUpper_ == floor(originalUpper_ + 0.5));
    }
}

// Copy constructor
OsiUsesBiLinear::OsiUsesBiLinear ( const OsiUsesBiLinear & rhs)
        : OsiSimpleInteger(rhs),
        numberBiLinear_(0),
        type_(rhs.type_),
        objects_(NULL)

{
}

// Clone
OsiObject *
OsiUsesBiLinear::clone() const
{
    return new OsiUsesBiLinear(*this);
}

// Assignment operator
OsiUsesBiLinear &
OsiUsesBiLinear::operator=( const OsiUsesBiLinear & rhs)
{
    if (this != &rhs) {
        OsiSimpleInteger::operator=(rhs);
        delete [] objects_;
        numberBiLinear_ = 0;
        type_ = rhs.type_;
        objects_ = NULL;
    }
    return *this;
}

// Destructor
OsiUsesBiLinear::~OsiUsesBiLinear ()
{
    delete [] objects_;
}
// Infeasibility - large is 0.5
double
OsiUsesBiLinear::infeasibility(const OsiBranchingInformation * info, int & whichWay) const
{
    assert (type_ == 0); // just continuous for now
    double value = info->solution_[columnNumber_];
    value = CoinMax(value, info->lower_[columnNumber_]);
    value = CoinMin(value, info->upper_[columnNumber_]);
    infeasibility_ = 0.0;
    for (int i = 0; i < numberBiLinear_; i++) {
        OsiBiLinear * obj = dynamic_cast<OsiBiLinear *> (objects_[i]);
        assert (obj);
        //obj->getPseudoShadow(info);
        //infeasibility_ += objects_[i]->infeasibility(info,whichWay);
        infeasibility_ += obj->getMovement(info);
    }
    bool satisfied = false;
    whichWay = -1;
    if (!infeasibility_) {
        otherInfeasibility_ = 1.0;
        satisfied = true;
        infeasibility_ = 0.0;
    } else {
        otherInfeasibility_ = 10.0 * infeasibility_;
        if (value - info->lower_[columnNumber_] >
                info->upper_[columnNumber_] - value)
            whichWay = 1;
        else
            whichWay = -1;
    }
    if (preferredWay_ >= 0 && !satisfied)
        whichWay = preferredWay_;
    whichWay_ = static_cast<short int>(whichWay);
    return infeasibility_;
}
// Creates a branching object
OsiBranchingObject *
OsiUsesBiLinear::createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const
{
    double value = info->solution_[columnNumber_];
    value = CoinMax(value, info->lower_[columnNumber_]);
    value = CoinMin(value, info->upper_[columnNumber_]);
    assert (info->upper_[columnNumber_] > info->lower_[columnNumber_]);
    double nearest = floor(value + 0.5);
    double integerTolerance = info->integerTolerance_;
    if (fabs(value - nearest) < integerTolerance) {
        // adjust value
        if (nearest != info->upper_[columnNumber_])
            value = nearest + 2.0 * integerTolerance;
        else
            value = nearest - 2.0 * integerTolerance;
    }
    OsiBranchingObject * branch = new OsiIntegerBranchingObject(solver, this, way,
            value, value, value);
    return branch;
}
// This looks at solution and sets bounds to contain solution
/** More precisely: it first forces the variable within the existing
    bounds, and then tightens the bounds to fix the variable at the
    nearest integer value.
*/
double
OsiUsesBiLinear::feasibleRegion(OsiSolverInterface * solver,
                                const OsiBranchingInformation * info) const
{
    double value = info->solution_[columnNumber_];
    double newValue = CoinMax(value, info->lower_[columnNumber_]);
    newValue = CoinMin(newValue, info->upper_[columnNumber_]);
    solver->setColLower(columnNumber_, newValue);
    solver->setColUpper(columnNumber_, newValue);
    return fabs(value - newValue);
}
// Add all bi-linear objects
void
OsiUsesBiLinear::addBiLinearObjects(OsiSolverLink * solver)
{
    delete [] objects_;
    numberBiLinear_ = 0;
    OsiObject ** objects = solver->objects();
    int i;
    int numberObjects = solver->numberObjects();
    for (i = 0; i < numberObjects; i++) {
        OsiObject * obj = objects[i];
        OsiBiLinear * objB = dynamic_cast<OsiBiLinear *> (obj);
        if (objB) {
            if (objB->xColumn() == columnNumber_ || objB->yColumn() == columnNumber_)
                numberBiLinear_++;
        }
    }
    if (numberBiLinear_) {
        objects_ = new OsiObject * [numberBiLinear_];
        numberBiLinear_ = 0;
        for (i = 0; i < numberObjects; i++) {
            OsiObject * obj = objects[i];
            OsiBiLinear * objB = dynamic_cast<OsiBiLinear *> (obj);
            if (objB) {
                if (objB->xColumn() == columnNumber_ || objB->yColumn() == columnNumber_)
                    objects_[numberBiLinear_++] = obj;;
            }
        }
    } else {
        objects_ = NULL;
    }
}

