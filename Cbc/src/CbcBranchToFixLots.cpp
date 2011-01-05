// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/13/2009-- carved out of CbcBranchCut

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchCut.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"
#include "CbcBranchToFixLots.hpp"

/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
CbcBranchToFixLots::CbcBranchToFixLots ()
        : CbcBranchCut(),
        djTolerance_(COIN_DBL_MAX),
        fractionFixed_(1.0),
        mark_(NULL),
        depth_(-1),
        numberClean_(0),
        alwaysCreate_(false)
{
}

/* Useful constructor - passed reduced cost tolerance and fraction we would like fixed.
   Also depth level to do at.
   Also passed number of 1 rows which when clean triggers fix
   Always does if all 1 rows cleaned up and number>0 or if fraction columns reached
   Also whether to create branch if can't reach fraction.
*/
CbcBranchToFixLots::CbcBranchToFixLots (CbcModel * model, double djTolerance,
                                        double fractionFixed, int depth,
                                        int numberClean,
                                        const char * mark, bool alwaysCreate)
        : CbcBranchCut(model)
{
    djTolerance_ = djTolerance;
    fractionFixed_ = fractionFixed;
    if (mark) {
        int numberColumns = model->getNumCols();
        mark_ = new char[numberColumns];
        memcpy(mark_, mark, numberColumns);
    } else {
        mark_ = NULL;
    }
    depth_ = depth;
    assert (model);
    OsiSolverInterface * solver = model_->solver();
    matrixByRow_ = *solver->getMatrixByRow();
    numberClean_ = numberClean;
    alwaysCreate_ = alwaysCreate;
}
// Copy constructor
CbcBranchToFixLots::CbcBranchToFixLots ( const CbcBranchToFixLots & rhs)
        : CbcBranchCut(rhs)
{
    djTolerance_ = rhs.djTolerance_;
    fractionFixed_ = rhs.fractionFixed_;
    int numberColumns = model_->getNumCols();
    mark_ = CoinCopyOfArray(rhs.mark_, numberColumns);
    matrixByRow_ = rhs.matrixByRow_;
    depth_ = rhs.depth_;
    numberClean_ = rhs.numberClean_;
    alwaysCreate_ = rhs.alwaysCreate_;
}

// Clone
CbcObject *
CbcBranchToFixLots::clone() const
{
    return new CbcBranchToFixLots(*this);
}

// Assignment operator
CbcBranchToFixLots &
CbcBranchToFixLots::operator=( const CbcBranchToFixLots & rhs)
{
    if (this != &rhs) {
        CbcBranchCut::operator=(rhs);
        djTolerance_ = rhs.djTolerance_;
        fractionFixed_ = rhs.fractionFixed_;
        int numberColumns = model_->getNumCols();
        delete [] mark_;
        mark_ = CoinCopyOfArray(rhs.mark_, numberColumns);
        matrixByRow_ = rhs.matrixByRow_;
        depth_ = rhs.depth_;
        numberClean_ = rhs.numberClean_;
        alwaysCreate_ = rhs.alwaysCreate_;
    }
    return *this;
}

// Destructor
CbcBranchToFixLots::~CbcBranchToFixLots ()
{
    delete [] mark_;
}
CbcBranchingObject *
CbcBranchToFixLots::createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * /*info*/, int /*way*/)
{
    // by default way must be -1
    //assert (way==-1);
    //OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    const double * dj = solver->getReducedCost();
    int i;
    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
    // make smaller ?
    double tolerance = CoinMin(1.0e-8, integerTolerance);
    // How many fixed are we aiming at
    int wantedFixed = static_cast<int> (static_cast<double>(numberIntegers) * fractionFixed_);
    int nSort = 0;
    int numberFixed = 0;
    int numberColumns = solver->getNumCols();
    int * sort = new int[numberColumns];
    double * dsort = new double[numberColumns];
    if (djTolerance_ != -1.234567) {
        int type = shallWe();
        assert (type);
        // Take clean first
        if (type == 1) {
            for (i = 0; i < numberIntegers; i++) {
                int iColumn = integerVariable[i];
                if (upper[iColumn] > lower[iColumn]) {
                    if (!mark_ || !mark_[iColumn]) {
                        if (solution[iColumn] < lower[iColumn] + tolerance) {
                            if (dj[iColumn] > djTolerance_) {
                                dsort[nSort] = -dj[iColumn];
                                sort[nSort++] = iColumn;
                            }
                        } else if (solution[iColumn] > upper[iColumn] - tolerance) {
                            if (dj[iColumn] < -djTolerance_) {
                                dsort[nSort] = dj[iColumn];
                                sort[nSort++] = iColumn;
                            }
                        }
                    }
                } else {
                    numberFixed++;
                }
            }
            // sort
            CoinSort_2(dsort, dsort + nSort, sort);
            nSort = CoinMin(nSort, wantedFixed - numberFixed);
        } else if (type < 10) {
            int i;
            //const double * rowLower = solver->getRowLower();
            const double * rowUpper = solver->getRowUpper();
            // Row copy
            const double * elementByRow = matrixByRow_.getElements();
            const int * column = matrixByRow_.getIndices();
            const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
            const int * rowLength = matrixByRow_.getVectorLengths();
            const double * columnLower = solver->getColLower();
            const double * columnUpper = solver->getColUpper();
            const double * solution = solver->getColSolution();
            int numberColumns = solver->getNumCols();
            int numberRows = solver->getNumRows();
            for (i = 0; i < numberColumns; i++) {
                sort[i] = i;
                if (columnLower[i] != columnUpper[i]) {
                    dsort[i] = 1.0e100;
                } else {
                    dsort[i] = 1.0e50;
                    numberFixed++;
                }
            }
            for (i = 0; i < numberRows; i++) {
                double rhsValue = rowUpper[i];
                bool oneRow = true;
                // check elements
                int numberUnsatisfied = 0;
                for (int j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
                    int iColumn = column[j];
                    double value = elementByRow[j];
                    double solValue = solution[iColumn];
                    if (columnLower[iColumn] != columnUpper[iColumn]) {
                        if (solValue < 1.0 - integerTolerance && solValue > integerTolerance)
                            numberUnsatisfied++;
                        if (value != 1.0) {
                            oneRow = false;
                            break;
                        }
                    } else {
                        rhsValue -= value * floor(solValue + 0.5);
                    }
                }
                if (oneRow && rhsValue <= 1.0 + tolerance) {
                    if (!numberUnsatisfied) {
                        for (int j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
                            int iColumn = column[j];
                            if (dsort[iColumn] > 1.0e50) {
                                dsort[iColumn] = 0;
                                nSort++;
                            }
                        }
                    }
                }
            }
            // sort
            CoinSort_2(dsort, dsort + numberColumns, sort);
        } else {
            // new way
            for (i = 0; i < numberIntegers; i++) {
                int iColumn = integerVariable[i];
                if (upper[iColumn] > lower[iColumn]) {
                    if (!mark_ || !mark_[iColumn]) {
                        double distanceDown = solution[iColumn] - lower[iColumn];
                        double distanceUp = upper[iColumn] - solution[iColumn];
                        double distance = CoinMin(distanceDown, distanceUp);
                        if (distance > 0.001 && distance < 0.5) {
                            dsort[nSort] = distance;
                            sort[nSort++] = iColumn;
                        }
                    }
                }
            }
            // sort
            CoinSort_2(dsort, dsort + nSort, sort);
            int n = 0;
            double sum = 0.0;
            for (int k = 0; k < nSort; k++) {
                sum += dsort[k];
                if (sum <= djTolerance_)
                    n = k;
                else
                    break;
            }
            nSort = CoinMin(n, numberClean_ / 1000000);
        }
    } else {
#define FIX_IF_LESS -0.1
        // 3 in same row and sum <FIX_IF_LESS?
        int numberRows = matrixByRow_.getNumRows();
        const double * solution = model_->testSolution();
        const int * column = matrixByRow_.getIndices();
        const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
        const int * rowLength = matrixByRow_.getVectorLengths();
        double bestSum = 1.0;
        int nBest = -1;
        int kRow = -1;
        OsiSolverInterface * solver = model_->solver();
        for (int i = 0; i < numberRows; i++) {
            int numberUnsatisfied = 0;
            double sum = 0.0;
            for (int j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
                int iColumn = column[j];
                if (solver->isInteger(iColumn)) {
                    double solValue = solution[iColumn];
                    if (solValue > 1.0e-5 && solValue < FIX_IF_LESS) {
                        numberUnsatisfied++;
                        sum += solValue;
                    }
                }
            }
            if (numberUnsatisfied >= 3 && sum < FIX_IF_LESS) {
                // possible
                if (numberUnsatisfied > nBest ||
                        (numberUnsatisfied == nBest && sum < bestSum)) {
                    nBest = numberUnsatisfied;
                    bestSum = sum;
                    kRow = i;
                }
            }
        }
        assert (nBest > 0);
        for (int j = rowStart[kRow]; j < rowStart[kRow] + rowLength[kRow]; j++) {
            int iColumn = column[j];
            if (solver->isInteger(iColumn)) {
                double solValue = solution[iColumn];
                if (solValue > 1.0e-5 && solValue < FIX_IF_LESS) {
                    sort[nSort++] = iColumn;
                }
            }
        }
    }
    OsiRowCut down;
    down.setLb(-COIN_DBL_MAX);
    double rhs = 0.0;
    for (i = 0; i < nSort; i++) {
        int iColumn = sort[i];
        double distanceDown = solution[iColumn] - lower[iColumn];
        double distanceUp = upper[iColumn] - solution[iColumn];
        if (distanceDown < distanceUp) {
            rhs += lower[iColumn];
            dsort[i] = 1.0;
        } else {
            rhs -= upper[iColumn];
            dsort[i] = -1.0;
        }
    }
    down.setUb(rhs);
    down.setRow(nSort, sort, dsort);
    down.setEffectiveness(COIN_DBL_MAX); // so will persist
    delete [] sort;
    delete [] dsort;
    // up is same - just with rhs changed
    OsiRowCut up = down;
    up.setLb(rhs + 1.0);
    up.setUb(COIN_DBL_MAX);
    // Say can fix one way
    CbcCutBranchingObject * newObject =
        new CbcCutBranchingObject(model_, down, up, true);
    if (model_->messageHandler()->logLevel() > 1)
        printf("creating cut in CbcBranchCut\n");
    return newObject;
}
/* Does a lot of the work,
   Returns 0 if no good, 1 if dj, 2 if clean, 3 if both
   10 if branching on ones away from bound
*/
int
CbcBranchToFixLots::shallWe() const
{
    int returnCode = 0;
    OsiSolverInterface * solver = model_->solver();
    int numberRows = matrixByRow_.getNumRows();
    //if (numberRows!=solver->getNumRows())
    //return 0;
    const double * solution = model_->testSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    const double * dj = solver->getReducedCost();
    int i;
    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();
    if (numberClean_ > 1000000) {
        int wanted = numberClean_ % 1000000;
        int * sort = new int[numberIntegers];
        double * dsort = new double[numberIntegers];
        int nSort = 0;
        for (i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            if (upper[iColumn] > lower[iColumn]) {
                if (!mark_ || !mark_[iColumn]) {
                    double distanceDown = solution[iColumn] - lower[iColumn];
                    double distanceUp = upper[iColumn] - solution[iColumn];
                    double distance = CoinMin(distanceDown, distanceUp);
                    if (distance > 0.001 && distance < 0.5) {
                        dsort[nSort] = distance;
                        sort[nSort++] = iColumn;
                    }
                }
            }
        }
        // sort
        CoinSort_2(dsort, dsort + nSort, sort);
        int n = 0;
        double sum = 0.0;
        for (int k = 0; k < nSort; k++) {
            sum += dsort[k];
            if (sum <= djTolerance_)
                n = k;
            else
                break;
        }
        delete [] sort;
        delete [] dsort;
        return (n >= wanted) ? 10 : 0;
    }
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
    // make smaller ?
    double tolerance = CoinMin(1.0e-8, integerTolerance);
    // How many fixed are we aiming at
    int wantedFixed = static_cast<int> (static_cast<double>(numberIntegers) * fractionFixed_);
    if (djTolerance_ < 1.0e10) {
        int nSort = 0;
        int numberFixed = 0;
        for (i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            if (upper[iColumn] > lower[iColumn]) {
                if (!mark_ || !mark_[iColumn]) {
                    if (solution[iColumn] < lower[iColumn] + tolerance) {
                        if (dj[iColumn] > djTolerance_) {
                            nSort++;
                        }
                    } else if (solution[iColumn] > upper[iColumn] - tolerance) {
                        if (dj[iColumn] < -djTolerance_) {
                            nSort++;
                        }
                    }
                }
            } else {
                numberFixed++;
            }
        }
        if (numberFixed + nSort < wantedFixed && !alwaysCreate_) {
            returnCode = 0;
        } else if (numberFixed < wantedFixed) {
            returnCode = 1;
        } else {
            returnCode = 0;
        }
    }
    if (numberClean_) {
        // see how many rows clean
        int i;
        //const double * rowLower = solver->getRowLower();
        const double * rowUpper = solver->getRowUpper();
        // Row copy
        const double * elementByRow = matrixByRow_.getElements();
        const int * column = matrixByRow_.getIndices();
        const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
        const int * rowLength = matrixByRow_.getVectorLengths();
        const double * columnLower = solver->getColLower();
        const double * columnUpper = solver->getColUpper();
        const double * solution = solver->getColSolution();
        int numberClean = 0;
        bool someToDoYet = false;
        int numberColumns = solver->getNumCols();
        char * mark = new char[numberColumns];
        int numberFixed = 0;
        for (i = 0; i < numberColumns; i++) {
            if (columnLower[i] != columnUpper[i]) {
                mark[i] = 0;
            } else {
                mark[i] = 1;
                numberFixed++;
            }
        }
        int numberNewFixed = 0;
        for (i = 0; i < numberRows; i++) {
            double rhsValue = rowUpper[i];
            bool oneRow = true;
            // check elements
            int numberUnsatisfied = 0;
            for (int j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
                int iColumn = column[j];
                double value = elementByRow[j];
                double solValue = solution[iColumn];
                if (columnLower[iColumn] != columnUpper[iColumn]) {
                    if (solValue < 1.0 - integerTolerance && solValue > integerTolerance)
                        numberUnsatisfied++;
                    if (value != 1.0) {
                        oneRow = false;
                        break;
                    }
                } else {
                    rhsValue -= value * floor(solValue + 0.5);
                }
            }
            if (oneRow && rhsValue <= 1.0 + tolerance) {
                if (numberUnsatisfied) {
                    someToDoYet = true;
                } else {
                    numberClean++;
                    for (int j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
                        int iColumn = column[j];
                        if (columnLower[iColumn] != columnUpper[iColumn] && !mark[iColumn]) {
                            mark[iColumn] = 1;
                            numberNewFixed++;
                        }
                    }
                }
            }
        }
        delete [] mark;
        //printf("%d clean, %d old fixed, %d new fixed\n",
        //   numberClean,numberFixed,numberNewFixed);
        if (someToDoYet && numberClean < numberClean_
                && numberNewFixed + numberFixed < wantedFixed) {
        } else if (numberFixed < wantedFixed) {
            returnCode |= 2;
        } else {
        }
    }
    return returnCode;
}
double
CbcBranchToFixLots::infeasibility(const OsiBranchingInformation * /*info*/,
                                  int &preferredWay) const
{
    preferredWay = -1;
    CbcNode * node = model_->currentNode();
    int depth;
    if (node)
        depth = CoinMax(node->depth(), 0);
    else
        return 0.0;
    if (depth_ < 0) {
        return 0.0;
    } else if (depth_ > 0) {
        if ((depth % depth_) != 0)
            return 0.0;
    }
    if (djTolerance_ != -1.234567) {
        if (!shallWe())
            return 0.0;
        else
            return 1.0e20;
    } else {
        // See if 3 in same row and sum <FIX_IF_LESS?
        int numberRows = matrixByRow_.getNumRows();
        const double * solution = model_->testSolution();
        const int * column = matrixByRow_.getIndices();
        const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
        const int * rowLength = matrixByRow_.getVectorLengths();
        double bestSum = 1.0;
        int nBest = -1;
        OsiSolverInterface * solver = model_->solver();
        for (int i = 0; i < numberRows; i++) {
            int numberUnsatisfied = 0;
            double sum = 0.0;
            for (int j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
                int iColumn = column[j];
                if (solver->isInteger(iColumn)) {
                    double solValue = solution[iColumn];
                    if (solValue > 1.0e-5 && solValue < FIX_IF_LESS) {
                        numberUnsatisfied++;
                        sum += solValue;
                    }
                }
            }
            if (numberUnsatisfied >= 3 && sum < FIX_IF_LESS) {
                // possible
                if (numberUnsatisfied > nBest ||
                        (numberUnsatisfied == nBest && sum < bestSum)) {
                    nBest = numberUnsatisfied;
                    bestSum = sum;
                }
            }
        }
        if (nBest > 0)
            return 1.0e20;
        else
            return 0.0;
    }
}
// Redoes data when sequence numbers change
void
CbcBranchToFixLots::redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns)
{
    model_ = model;
    if (mark_) {
        OsiSolverInterface * solver = model_->solver();
        int numberColumnsNow = solver->getNumCols();
        char * temp = new char[numberColumnsNow];
        memset(temp, 0, numberColumnsNow);
        for (int i = 0; i < numberColumns; i++) {
            int j = originalColumns[i];
            temp[i] = mark_[j];
        }
        delete [] mark_;
        mark_ = temp;
    }
    OsiSolverInterface * solver = model_->solver();
    matrixByRow_ = *solver->getMatrixByRow();
}

