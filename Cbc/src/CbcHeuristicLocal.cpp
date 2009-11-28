/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcBranchActual.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"

// Default Constructor
CbcHeuristicLocal::CbcHeuristicLocal()
        : CbcHeuristic()
{
    numberSolutions_ = 0;
    swap_ = 0;
    used_ = NULL;
}

// Constructor with model - assumed before cuts

CbcHeuristicLocal::CbcHeuristicLocal(CbcModel & model)
        : CbcHeuristic(model)
{
    numberSolutions_ = 0;
    swap_ = 0;
    // Get a copy of original matrix
    assert(model.solver());
    if (model.solver()->getNumRows()) {
        matrix_ = *model.solver()->getMatrixByCol();
    }
    int numberColumns = model.solver()->getNumCols();
    used_ = new int[numberColumns];
    memset(used_, 0, numberColumns*sizeof(int));
}

// Destructor
CbcHeuristicLocal::~CbcHeuristicLocal ()
{
    delete [] used_;
}

// Clone
CbcHeuristic *
CbcHeuristicLocal::clone() const
{
    return new CbcHeuristicLocal(*this);
}
// Create C++ lines to get to current state
void
CbcHeuristicLocal::generateCpp( FILE * fp)
{
    CbcHeuristicLocal other;
    fprintf(fp, "0#include \"CbcHeuristicLocal.hpp\"\n");
    fprintf(fp, "3  CbcHeuristicLocal heuristicLocal(*cbcModel);\n");
    CbcHeuristic::generateCpp(fp, "heuristicLocal");
    if (swap_ != other.swap_)
        fprintf(fp, "3  heuristicLocal.setSearchType(%d);\n", swap_);
    else
        fprintf(fp, "4  heuristicLocal.setSearchType(%d);\n", swap_);
    fprintf(fp, "3  cbcModel->addHeuristic(&heuristicLocal);\n");
}

// Copy constructor
CbcHeuristicLocal::CbcHeuristicLocal(const CbcHeuristicLocal & rhs)
        :
        CbcHeuristic(rhs),
        matrix_(rhs.matrix_),
        numberSolutions_(rhs.numberSolutions_),
        swap_(rhs.swap_)
{
    if (model_ && rhs.used_) {
        int numberColumns = model_->solver()->getNumCols();
        used_ = CoinCopyOfArray(rhs.used_, numberColumns);
    } else {
        used_ = NULL;
    }
}

// Assignment operator
CbcHeuristicLocal &
CbcHeuristicLocal::operator=( const CbcHeuristicLocal & rhs)
{
    if (this != &rhs) {
        CbcHeuristic::operator=(rhs);
        matrix_ = rhs.matrix_;
        numberSolutions_ = rhs.numberSolutions_;
        swap_ = rhs.swap_;
        delete [] used_;
        if (model_ && rhs.used_) {
            int numberColumns = model_->solver()->getNumCols();
            used_ = CoinCopyOfArray(rhs.used_, numberColumns);
        } else {
            used_ = NULL;
        }
    }
    return *this;
}

// Resets stuff if model changes
void
CbcHeuristicLocal::resetModel(CbcModel * /*model*/)
{
    //CbcHeuristic::resetModel(model);
    delete [] used_;
    if (model_ && used_) {
        int numberColumns = model_->solver()->getNumCols();
        used_ = new int[numberColumns];
        memset(used_, 0, numberColumns*sizeof(int));
    } else {
        used_ = NULL;
    }
}
// This version fixes stuff and does IP
int
CbcHeuristicLocal::solutionFix(double & objectiveValue,
                               double * newSolution,
                               const int * /*keep*/)
{
    numCouldRun_++;
    // See if to do
    if (!when() || (when() == 1 && model_->phase() != 1))
        return 0; // switched off
    // Don't do if it was this heuristic which found solution!
    if (this == model_->lastHeuristic())
        return 0;
    OsiSolverInterface * newSolver = model_->continuousSolver()->clone();
    const double * colLower = newSolver->getColLower();
    //const double * colUpper = newSolver->getColUpper();

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();

    int i;
    int nFix = 0;
    for (i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        const OsiObject * object = model_->object(i);
        // get original bounds
        double originalLower;
        double originalUpper;
        getIntegerInformation( object, originalLower, originalUpper);
        newSolver->setColLower(iColumn, CoinMax(colLower[iColumn], originalLower));
        if (!used_[iColumn]) {
            newSolver->setColUpper(iColumn, colLower[iColumn]);
            nFix++;
        }
    }
    int returnCode = 0;
#ifdef CLP_INVESTIGATE2
    printf("Fixing %d out of %d (%d continuous)\n",
           nFix, numberIntegers, newSolver->getNumCols() - numberIntegers);
#endif
    if (nFix*10 <= numberIntegers) {
        // see if we can fix more
        int * which = new int [2*(numberIntegers-nFix)];
        int * sort = which + (numberIntegers - nFix);
        int n = 0;
        for (i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            if (used_[iColumn]) {
                which[n] = iColumn;
                sort[n++] = used_[iColumn];
            }
        }
        CoinSort_2(sort, sort + n, which);
        // only half fixed in total
        n = CoinMin(n, numberIntegers / 2 - nFix);
        int allow = CoinMax(numberSolutions_ - 2, sort[0]);
        int nFix2 = 0;
        for (i = 0; i < n; i++) {
            int iColumn = integerVariable[i];
            if (used_[iColumn] <= allow) {
                newSolver->setColUpper(iColumn, colLower[iColumn]);
                nFix2++;
            } else {
                break;
            }
        }
        delete [] which;
        nFix += nFix2;
        printf("Number fixed increased from %d to %d\n",
               nFix - nFix2, nFix);
    }
    if (nFix*10 > numberIntegers) {
        returnCode = smallBranchAndBound(newSolver, numberNodes_, newSolution, objectiveValue,
                                         objectiveValue, "CbcHeuristicLocal");
        if (returnCode < 0) {
            returnCode = 0; // returned on size
            int numberColumns = newSolver->getNumCols();
            int numberContinuous = numberColumns - numberIntegers;
            if (numberContinuous > 2*numberIntegers &&
                    nFix*10 < numberColumns) {
#define LOCAL_FIX_CONTINUOUS
#ifdef LOCAL_FIX_CONTINUOUS
                //const double * colUpper = newSolver->getColUpper();
                const double * colLower = newSolver->getColLower();
                int nAtLb = 0;
                //double sumDj=0.0;
                const double * dj = newSolver->getReducedCost();
                double direction = newSolver->getObjSense();
                for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (!newSolver->isInteger(iColumn)) {
                        if (!used_[iColumn]) {
                            //double djValue = dj[iColumn]*direction;
                            nAtLb++;
                            //sumDj += djValue;
                        }
                    }
                }
                if (nAtLb) {
                    // fix some continuous
                    double * sort = new double[nAtLb];
                    int * which = new int [nAtLb];
                    //double threshold = CoinMax((0.01*sumDj)/static_cast<double>(nAtLb),1.0e-6);
                    int nFix2 = 0;
                    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                        if (!newSolver->isInteger(iColumn)) {
                            if (!used_[iColumn]) {
                                double djValue = dj[iColumn] * direction;
                                if (djValue > 1.0e-6) {
                                    sort[nFix2] = -djValue;
                                    which[nFix2++] = iColumn;
                                }
                            }
                        }
                    }
                    CoinSort_2(sort, sort + nFix2, which);
                    int divisor = 2;
                    nFix2 = CoinMin(nFix2, (numberColumns - nFix) / divisor);
                    for (int i = 0; i < nFix2; i++) {
                        int iColumn = which[i];
                        newSolver->setColUpper(iColumn, colLower[iColumn]);
                    }
                    delete [] sort;
                    delete [] which;
#ifdef CLP_INVESTIGATE2
                    printf("%d integers have zero value, and %d continuous fixed at lb\n",
                           nFix, nFix2);
#endif
                    returnCode = smallBranchAndBound(newSolver,
                                                     numberNodes_, newSolution,
                                                     objectiveValue,
                                                     objectiveValue, "CbcHeuristicLocal");
                    if (returnCode < 0)
                        returnCode = 0; // returned on size
                }
#endif
            }
        }
    }
    if ((returnCode&2) != 0) {
        // could add cut
        returnCode &= ~2;
    }

    delete newSolver;
    return returnCode;
}
/*
  First tries setting a variable to better value.  If feasible then
  tries setting others.  If not feasible then tries swaps
  Returns 1 if solution, 0 if not */
int
CbcHeuristicLocal::solution(double & solutionValue,
                            double * betterSolution)
{

    numCouldRun_++;
    if (numberSolutions_ == model_->getSolutionCount())
        return 0;
    numberSolutions_ = model_->getSolutionCount();
    if ((model_->getNumCols() > 100000 && model_->getNumCols() >
            10*model_->getNumRows()) || numberSolutions_ <= 1)
        return 0; // probably not worth it
    // worth trying

    OsiSolverInterface * solver = model_->solver();
    const double * rowLower = solver->getRowLower();
    const double * rowUpper = solver->getRowUpper();
    const double * solution = model_->bestSolution();
    if (!solution)
        return 0; // No solution found yet
    const double * objective = solver->getObjCoefficients();
    double primalTolerance;
    solver->getDblParam(OsiPrimalTolerance, primalTolerance);

    int numberRows = matrix_.getNumRows();

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();

    int i;
    double direction = solver->getObjSense();
    double newSolutionValue = model_->getObjValue() * direction;
    int returnCode = 0;
    numRuns_++;
    // Column copy
    const double * element = matrix_.getElements();
    const int * row = matrix_.getIndices();
    const CoinBigIndex * columnStart = matrix_.getVectorStarts();
    const int * columnLength = matrix_.getVectorLengths();

    // Get solution array for heuristic solution
    int numberColumns = solver->getNumCols();
    double * newSolution = new double [numberColumns];
    memcpy(newSolution, solution, numberColumns*sizeof(double));
#ifdef LOCAL_FIX_CONTINUOUS
    // mark continuous used
    const double * columnLower = solver->getColLower();
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (!solver->isInteger(iColumn)) {
            if (solution[iColumn] > columnLower[iColumn] + 1.0e-8)
                used_[iColumn] = numberSolutions_;
        }
    }
#endif

    // way is 1 if down possible, 2 if up possible, 3 if both possible
    char * way = new char[numberIntegers];
    // corrected costs
    double * cost = new double[numberIntegers];
    // for array to mark infeasible rows after iColumn branch
    char * mark = new char[numberRows];
    memset(mark, 0, numberRows);
    // space to save values so we don't introduce rounding errors
    double * save = new double[numberRows];

    // clean solution
    for (i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        const OsiObject * object = model_->object(i);
        // get original bounds
        double originalLower;
        double originalUpper;
        getIntegerInformation( object, originalLower, originalUpper);
        double value = newSolution[iColumn];
        if (value < originalLower) {
            value = originalLower;
            newSolution[iColumn] = value;
        } else if (value > originalUpper) {
            value = originalUpper;
            newSolution[iColumn] = value;
        }
        double nearest = floor(value + 0.5);
        //assert(fabs(value-nearest)<10.0*primalTolerance);
        value = nearest;
        newSolution[iColumn] = nearest;
        // if away from lower bound mark that fact
        if (nearest > originalLower) {
            used_[iColumn] = numberSolutions_;
        }
        cost[i] = direction * objective[iColumn];
        int iway = 0;

        if (value > originalLower + 0.5)
            iway = 1;
        if (value < originalUpper - 0.5)
            iway |= 2;
        way[i] = static_cast<char>(iway);
    }
    // get row activities
    double * rowActivity = new double[numberRows];
    memset(rowActivity, 0, numberRows*sizeof(double));

    for (i = 0; i < numberColumns; i++) {
        int j;
        double value = newSolution[i];
        if (value) {
            for (j = columnStart[i];
                    j < columnStart[i] + columnLength[i]; j++) {
                int iRow = row[j];
                rowActivity[iRow] += value * element[j];
            }
        }
    }
    // check was feasible - if not adjust (cleaning may move)
    // if very infeasible then give up
    bool tryHeuristic = true;
    for (i = 0; i < numberRows; i++) {
        if (rowActivity[i] < rowLower[i]) {
            if (rowActivity[i] < rowLower[i] - 10.0*primalTolerance)
                tryHeuristic = false;
            rowActivity[i] = rowLower[i];
        } else if (rowActivity[i] > rowUpper[i]) {
            if (rowActivity[i] < rowUpper[i] + 10.0*primalTolerance)
                tryHeuristic = false;
            rowActivity[i] = rowUpper[i];
        }
    }
    // Switch off if may take too long
    if (model_->getNumCols() > 10000 && model_->getNumCols() >
            10*model_->getNumRows())
        tryHeuristic = false;
    if (tryHeuristic) {

        // best change in objective
        double bestChange = 0.0;

        for (i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];

            double objectiveCoefficient = cost[i];
            int k;
            int j;
            int goodK = -1;
            int wayK = -1, wayI = -1;
            if ((way[i]&1) != 0) {
                int numberInfeasible = 0;
                // save row activities and adjust
                for (j = columnStart[iColumn];
                        j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                    int iRow = row[j];
                    save[iRow] = rowActivity[iRow];
                    rowActivity[iRow] -= element[j];
                    if (rowActivity[iRow] < rowLower[iRow] - primalTolerance ||
                            rowActivity[iRow] > rowUpper[iRow] + primalTolerance) {
                        // mark row
                        mark[iRow] = 1;
                        numberInfeasible++;
                    }
                }
                // try down
                for (k = i + 1; k < numberIntegers; k++) {
                    if ((way[k]&1) != 0) {
                        // try down
                        if (-objectiveCoefficient - cost[k] < bestChange) {
                            // see if feasible down
                            bool good = true;
                            int numberMarked = 0;
                            int kColumn = integerVariable[k];
                            for (j = columnStart[kColumn];
                                    j < columnStart[kColumn] + columnLength[kColumn]; j++) {
                                int iRow = row[j];
                                double newValue = rowActivity[iRow] - element[j];
                                if (newValue < rowLower[iRow] - primalTolerance ||
                                        newValue > rowUpper[iRow] + primalTolerance) {
                                    good = false;
                                    break;
                                } else if (mark[iRow]) {
                                    // made feasible
                                    numberMarked++;
                                }
                            }
                            if (good && numberMarked == numberInfeasible) {
                                // better solution
                                goodK = k;
                                wayK = -1;
                                wayI = -1;
                                bestChange = -objectiveCoefficient - cost[k];
                            }
                        }
                    }
                    if ((way[k]&2) != 0) {
                        // try up
                        if (-objectiveCoefficient + cost[k] < bestChange) {
                            // see if feasible up
                            bool good = true;
                            int numberMarked = 0;
                            int kColumn = integerVariable[k];
                            for (j = columnStart[kColumn];
                                    j < columnStart[kColumn] + columnLength[kColumn]; j++) {
                                int iRow = row[j];
                                double newValue = rowActivity[iRow] + element[j];
                                if (newValue < rowLower[iRow] - primalTolerance ||
                                        newValue > rowUpper[iRow] + primalTolerance) {
                                    good = false;
                                    break;
                                } else if (mark[iRow]) {
                                    // made feasible
                                    numberMarked++;
                                }
                            }
                            if (good && numberMarked == numberInfeasible) {
                                // better solution
                                goodK = k;
                                wayK = 1;
                                wayI = -1;
                                bestChange = -objectiveCoefficient + cost[k];
                            }
                        }
                    }
                }
                // restore row activities
                for (j = columnStart[iColumn];
                        j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                    int iRow = row[j];
                    rowActivity[iRow] = save[iRow];
                    mark[iRow] = 0;
                }
            }
            if ((way[i]&2) != 0) {
                int numberInfeasible = 0;
                // save row activities and adjust
                for (j = columnStart[iColumn];
                        j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                    int iRow = row[j];
                    save[iRow] = rowActivity[iRow];
                    rowActivity[iRow] += element[j];
                    if (rowActivity[iRow] < rowLower[iRow] - primalTolerance ||
                            rowActivity[iRow] > rowUpper[iRow] + primalTolerance) {
                        // mark row
                        mark[iRow] = 1;
                        numberInfeasible++;
                    }
                }
                // try up
                for (k = i + 1; k < numberIntegers; k++) {
                    if ((way[k]&1) != 0) {
                        // try down
                        if (objectiveCoefficient - cost[k] < bestChange) {
                            // see if feasible down
                            bool good = true;
                            int numberMarked = 0;
                            int kColumn = integerVariable[k];
                            for (j = columnStart[kColumn];
                                    j < columnStart[kColumn] + columnLength[kColumn]; j++) {
                                int iRow = row[j];
                                double newValue = rowActivity[iRow] - element[j];
                                if (newValue < rowLower[iRow] - primalTolerance ||
                                        newValue > rowUpper[iRow] + primalTolerance) {
                                    good = false;
                                    break;
                                } else if (mark[iRow]) {
                                    // made feasible
                                    numberMarked++;
                                }
                            }
                            if (good && numberMarked == numberInfeasible) {
                                // better solution
                                goodK = k;
                                wayK = -1;
                                wayI = 1;
                                bestChange = objectiveCoefficient - cost[k];
                            }
                        }
                    }
                    if ((way[k]&2) != 0) {
                        // try up
                        if (objectiveCoefficient + cost[k] < bestChange) {
                            // see if feasible up
                            bool good = true;
                            int numberMarked = 0;
                            int kColumn = integerVariable[k];
                            for (j = columnStart[kColumn];
                                    j < columnStart[kColumn] + columnLength[kColumn]; j++) {
                                int iRow = row[j];
                                double newValue = rowActivity[iRow] + element[j];
                                if (newValue < rowLower[iRow] - primalTolerance ||
                                        newValue > rowUpper[iRow] + primalTolerance) {
                                    good = false;
                                    break;
                                } else if (mark[iRow]) {
                                    // made feasible
                                    numberMarked++;
                                }
                            }
                            if (good && numberMarked == numberInfeasible) {
                                // better solution
                                goodK = k;
                                wayK = 1;
                                wayI = 1;
                                bestChange = objectiveCoefficient + cost[k];
                            }
                        }
                    }
                }
                // restore row activities
                for (j = columnStart[iColumn];
                        j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                    int iRow = row[j];
                    rowActivity[iRow] = save[iRow];
                    mark[iRow] = 0;
                }
            }
            if (goodK >= 0) {
                // we found something - update solution
                for (j = columnStart[iColumn];
                        j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                    int iRow = row[j];
                    rowActivity[iRow]  += wayI * element[j];
                }
                newSolution[iColumn] += wayI;
                int kColumn = integerVariable[goodK];
                for (j = columnStart[kColumn];
                        j < columnStart[kColumn] + columnLength[kColumn]; j++) {
                    int iRow = row[j];
                    rowActivity[iRow]  += wayK * element[j];
                }
                newSolution[kColumn] += wayK;
                // See if k can go further ?
                const OsiObject * object = model_->object(goodK);
                // get original bounds
                double originalLower;
                double originalUpper;
                getIntegerInformation( object, originalLower, originalUpper);

                double value = newSolution[kColumn];
                int iway = 0;

                if (value > originalLower + 0.5)
                    iway = 1;
                if (value < originalUpper - 0.5)
                    iway |= 2;
                way[goodK] = static_cast<char>(iway);
            }
        }
        if (bestChange + newSolutionValue < solutionValue) {
            // paranoid check
            memset(rowActivity, 0, numberRows*sizeof(double));

            for (i = 0; i < numberColumns; i++) {
                int j;
                double value = newSolution[i];
                if (value) {
                    for (j = columnStart[i];
                            j < columnStart[i] + columnLength[i]; j++) {
                        int iRow = row[j];
                        rowActivity[iRow] += value * element[j];
                    }
                }
            }
            int numberBad = 0;
            double sumBad = 0.0;
            // check was approximately feasible
            for (i = 0; i < numberRows; i++) {
                if (rowActivity[i] < rowLower[i]) {
                    sumBad += rowLower[i] - rowActivity[i];
                    if (rowActivity[i] < rowLower[i] - 10.0*primalTolerance)
                        numberBad++;
                } else if (rowActivity[i] > rowUpper[i]) {
                    sumBad += rowUpper[i] - rowActivity[i];
                    if (rowActivity[i] > rowUpper[i] + 10.0*primalTolerance)
                        numberBad++;
                }
            }
            if (!numberBad) {
                for (i = 0; i < numberIntegers; i++) {
                    int iColumn = integerVariable[i];
                    const OsiObject * object = model_->object(i);
                    // get original bounds
                    double originalLower;
                    double originalUpper;
                    getIntegerInformation( object, originalLower, originalUpper);

                    double value = newSolution[iColumn];
                    // if away from lower bound mark that fact
                    if (value > originalLower) {
                        used_[iColumn] = numberSolutions_;
                    }
                }
                // new solution
                memcpy(betterSolution, newSolution, numberColumns*sizeof(double));
                CoinWarmStartBasis * basis =
                    dynamic_cast<CoinWarmStartBasis *>(solver->getWarmStart()) ;
                if (basis) {
                    model_->setBestSolutionBasis(* basis);
                    delete basis;
                }
                returnCode = 1;
                solutionValue = newSolutionValue + bestChange;
            } else {
                // bad solution - should not happen so debug if see message
                printf("Local search got bad solution with %d infeasibilities summing to %g\n",
                       numberBad, sumBad);
            }
        }
    }
    delete [] newSolution;
    delete [] rowActivity;
    delete [] way;
    delete [] cost;
    delete [] save;
    delete [] mark;
    if (numberSolutions_ > 1 && swap_ == 1) {
        // try merge
        int returnCode2 = solutionFix( solutionValue, betterSolution, NULL);
        if (returnCode2)
            returnCode = 1;
    }
    return returnCode;
}
// update model
void CbcHeuristicLocal::setModel(CbcModel * model)
{
    model_ = model;
    // Get a copy of original matrix
    assert(model_->solver());
    if (model_->solver()->getNumRows()) {
        matrix_ = *model_->solver()->getMatrixByCol();
    }
    delete [] used_;
    int numberColumns = model->solver()->getNumCols();
    used_ = new int[numberColumns];
    memset(used_, 0, numberColumns*sizeof(int));
}
// Default Constructor
CbcHeuristicNaive::CbcHeuristicNaive()
        : CbcHeuristic()
{
    large_ = 1.0e6;
}

// Constructor with model - assumed before cuts

CbcHeuristicNaive::CbcHeuristicNaive(CbcModel & model)
        : CbcHeuristic(model)
{
    large_ = 1.0e6;
}

// Destructor
CbcHeuristicNaive::~CbcHeuristicNaive ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicNaive::clone() const
{
    return new CbcHeuristicNaive(*this);
}
// Create C++ lines to get to current state
void
CbcHeuristicNaive::generateCpp( FILE * fp)
{
    CbcHeuristicNaive other;
    fprintf(fp, "0#include \"CbcHeuristicLocal.hpp\"\n");
    fprintf(fp, "3  CbcHeuristicNaive naive(*cbcModel);\n");
    CbcHeuristic::generateCpp(fp, "naive");
    if (large_ != other.large_)
        fprintf(fp, "3  naive.setLarge(%g);\n", large_);
    else
        fprintf(fp, "4  naive.setLarge(%g);\n", large_);
    fprintf(fp, "3  cbcModel->addHeuristic(&naive);\n");
}

// Copy constructor
CbcHeuristicNaive::CbcHeuristicNaive(const CbcHeuristicNaive & rhs)
        :
        CbcHeuristic(rhs),
        large_(rhs.large_)
{
}

// Assignment operator
CbcHeuristicNaive &
CbcHeuristicNaive::operator=( const CbcHeuristicNaive & rhs)
{
    if (this != &rhs) {
        CbcHeuristic::operator=(rhs);
        large_ = rhs.large_;
    }
    return *this;
}

// Resets stuff if model changes
void
CbcHeuristicNaive::resetModel(CbcModel * model)
{
    CbcHeuristic::resetModel(model);
}
int
CbcHeuristicNaive::solution(double & solutionValue,
                            double * betterSolution)
{
    numCouldRun_++;
    // See if to do
    bool atRoot = model_->getNodeCount() == 0;
    int passNumber = model_->getCurrentPassNumber();
    if (!when() || (when() == 1 && model_->phase() != 1) || !atRoot || passNumber != 1)
        return 0; // switched off
    // Don't do if it was this heuristic which found solution!
    if (this == model_->lastHeuristic())
        return 0;
    numRuns_++;
    double cutoff;
    model_->solver()->getDblParam(OsiDualObjectiveLimit, cutoff);
    double direction = model_->solver()->getObjSense();
    cutoff *= direction;
    cutoff = CoinMin(cutoff, solutionValue);
    OsiSolverInterface * solver = model_->continuousSolver();
    if (!solver)
        solver = model_->solver();
    const double * colLower = solver->getColLower();
    const double * colUpper = solver->getColUpper();
    const double * objective = solver->getObjCoefficients();

    int numberColumns = model_->getNumCols();
    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();

    int i;
    bool solutionFound = false;
    CoinWarmStartBasis saveBasis;
    CoinWarmStartBasis * basis =
        dynamic_cast<CoinWarmStartBasis *>(solver->getWarmStart()) ;
    if (basis) {
        saveBasis = * basis;
        delete basis;
    }
    // First just fix all integers as close to zero as possible
    OsiSolverInterface * newSolver = cloneBut(7); // wassolver->clone();
    for (i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        double lower = colLower[iColumn];
        double upper = colUpper[iColumn];
        double value;
        if (lower > 0.0)
            value = lower;
        else if (upper < 0.0)
            value = upper;
        else
            value = 0.0;
        newSolver->setColLower(iColumn, value);
        newSolver->setColUpper(iColumn, value);
    }
    newSolver->initialSolve();
    if (newSolver->isProvenOptimal()) {
        double solValue = newSolver->getObjValue() * direction ;
        if (solValue < cutoff) {
            // we have a solution
            solutionFound = true;
            solutionValue = solValue;
            memcpy(betterSolution, newSolver->getColSolution(),
                   numberColumns*sizeof(double));
            printf("Naive fixing close to zero gave solution of %g\n", solutionValue);
            cutoff = solValue - model_->getCutoffIncrement();
        }
    }
    // Now fix all integers as close to zero if zero or large cost
    int nFix = 0;
    for (i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        double lower = colLower[iColumn];
        double upper = colUpper[iColumn];
        double value;
        if (fabs(objective[i]) > 0.0 && fabs(objective[i]) < large_) {
            nFix++;
            if (lower > 0.0)
                value = lower;
            else if (upper < 0.0)
                value = upper;
            else
                value = 0.0;
            newSolver->setColLower(iColumn, value);
            newSolver->setColUpper(iColumn, value);
        } else {
            // set back to original
            newSolver->setColLower(iColumn, lower);
            newSolver->setColUpper(iColumn, upper);
        }
    }
    const double * solution = solver->getColSolution();
    if (nFix) {
        newSolver->setWarmStart(&saveBasis);
        newSolver->setColSolution(solution);
        newSolver->initialSolve();
        if (newSolver->isProvenOptimal()) {
            double solValue = newSolver->getObjValue() * direction ;
            if (solValue < cutoff) {
                // try branch and bound
                double * newSolution = new double [numberColumns];
                printf("%d fixed after fixing costs\n", nFix);
                int returnCode = smallBranchAndBound(newSolver,
                                                     numberNodes_, newSolution,
                                                     solutionValue,
                                                     solutionValue, "CbcHeuristicNaive1");
                if (returnCode < 0)
                    returnCode = 0; // returned on size
                if ((returnCode&2) != 0) {
                    // could add cut
                    returnCode &= ~2;
                }
                if (returnCode == 1) {
                    // solution
                    solutionFound = true;
                    memcpy(betterSolution, newSolution,
                           numberColumns*sizeof(double));
                    printf("Naive fixing zeros gave solution of %g\n", solutionValue);
                    cutoff = solutionValue - model_->getCutoffIncrement();
                }
                delete [] newSolution;
            }
        }
    }
#if 1
    newSolver->setObjSense(-direction); // maximize
    newSolver->setWarmStart(&saveBasis);
    newSolver->setColSolution(solution);
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = solution[iColumn];
        double lower = colLower[iColumn];
        double upper = colUpper[iColumn];
        double newLower;
        double newUpper;
        if (newSolver->isInteger(iColumn)) {
            newLower = CoinMax(lower, floor(value) - 2.0);
            newUpper = CoinMin(upper, ceil(value) + 2.0);
        } else {
            newLower = CoinMax(lower, value - 1.0e5);
            newUpper = CoinMin(upper, value + 1.0e-5);
        }
        newSolver->setColLower(iColumn, newLower);
        newSolver->setColUpper(iColumn, newUpper);
    }
    newSolver->initialSolve();
    if (newSolver->isProvenOptimal()) {
        double solValue = newSolver->getObjValue() * direction ;
        if (solValue < cutoff) {
            nFix = 0;
            newSolver->setObjSense(direction); // correct direction
            //const double * thisSolution = newSolver->getColSolution();
            for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                double value = solution[iColumn];
                double lower = colLower[iColumn];
                double upper = colUpper[iColumn];
                double newLower = lower;
                double newUpper = upper;
                if (newSolver->isInteger(iColumn)) {
                    if (value < lower + 1.0e-6) {
                        nFix++;
                        newUpper = lower;
                    } else if (value > upper - 1.0e-6) {
                        nFix++;
                        newLower = upper;
                    } else {
                        newLower = CoinMax(lower, floor(value) - 2.0);
                        newUpper = CoinMin(upper, ceil(value) + 2.0);
                    }
                }
                newSolver->setColLower(iColumn, newLower);
                newSolver->setColUpper(iColumn, newUpper);
            }
            // try branch and bound
            double * newSolution = new double [numberColumns];
            printf("%d fixed after maximizing\n", nFix);
            int returnCode = smallBranchAndBound(newSolver,
                                                 numberNodes_, newSolution,
                                                 solutionValue,
                                                 solutionValue, "CbcHeuristicNaive1");
            if (returnCode < 0)
                returnCode = 0; // returned on size
            if ((returnCode&2) != 0) {
                // could add cut
                returnCode &= ~2;
            }
            if (returnCode == 1) {
                // solution
                solutionFound = true;
                memcpy(betterSolution, newSolution,
                       numberColumns*sizeof(double));
                printf("Naive maximizing gave solution of %g\n", solutionValue);
                cutoff = solutionValue - model_->getCutoffIncrement();
            }
            delete [] newSolution;
        }
    }
#endif
    delete newSolver;
    return solutionFound ? 1 : 0;
}
// update model
void CbcHeuristicNaive::setModel(CbcModel * model)
{
    model_ = model;
}
// Default Constructor
CbcHeuristicCrossover::CbcHeuristicCrossover()
        : CbcHeuristic(),
        numberSolutions_(0),
        useNumber_(3)
{
    setWhen(1);
}

// Constructor with model - assumed before cuts

CbcHeuristicCrossover::CbcHeuristicCrossover(CbcModel & model)
        : CbcHeuristic(model),
        numberSolutions_(0),
        useNumber_(3)
{
    setWhen(1);
    for (int i = 0; i < 10; i++)
        random_[i] = model.randomNumberGenerator()->randomDouble();
}

// Destructor
CbcHeuristicCrossover::~CbcHeuristicCrossover ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicCrossover::clone() const
{
    return new CbcHeuristicCrossover(*this);
}
// Create C++ lines to get to current state
void
CbcHeuristicCrossover::generateCpp( FILE * fp)
{
    CbcHeuristicCrossover other;
    fprintf(fp, "0#include \"CbcHeuristicLocal.hpp\"\n");
    fprintf(fp, "3  CbcHeuristicCrossover crossover(*cbcModel);\n");
    CbcHeuristic::generateCpp(fp, "crossover");
    if (useNumber_ != other.useNumber_)
        fprintf(fp, "3  crossover.setNumberSolutions(%d);\n", useNumber_);
    else
        fprintf(fp, "4  crossover.setNumberSolutions(%d);\n", useNumber_);
    fprintf(fp, "3  cbcModel->addHeuristic(&crossover);\n");
}

// Copy constructor
CbcHeuristicCrossover::CbcHeuristicCrossover(const CbcHeuristicCrossover & rhs)
        :
        CbcHeuristic(rhs),
        attempts_(rhs.attempts_),
        numberSolutions_(rhs.numberSolutions_),
        useNumber_(rhs.useNumber_)
{
    memcpy(random_, rhs.random_, 10*sizeof(double));
}

// Assignment operator
CbcHeuristicCrossover &
CbcHeuristicCrossover::operator=( const CbcHeuristicCrossover & rhs)
{
    if (this != &rhs) {
        CbcHeuristic::operator=(rhs);
        useNumber_ = rhs.useNumber_;
        attempts_ = rhs.attempts_;
        numberSolutions_ = rhs.numberSolutions_;
        memcpy(random_, rhs.random_, 10*sizeof(double));
    }
    return *this;
}

// Resets stuff if model changes
void
CbcHeuristicCrossover::resetModel(CbcModel * model)
{
    CbcHeuristic::resetModel(model);
}
int
CbcHeuristicCrossover::solution(double & solutionValue,
                                double * betterSolution)
{
    if (when_ == 0)
        return 0;
    numCouldRun_++;
    bool useBest = (numberSolutions_ != model_->getSolutionCount());
    if (!useBest && (when_ % 10) == 1)
        return 0;
    numberSolutions_ = model_->getSolutionCount();
    OsiSolverInterface * continuousSolver = model_->continuousSolver();
    int useNumber = CoinMin(model_->numberSavedSolutions(), useNumber_);
    if (useNumber < 2 || !continuousSolver)
        return 0;
    // Fix later
    if (!useBest)
        abort();
    numRuns_++;
    double cutoff;
    model_->solver()->getDblParam(OsiDualObjectiveLimit, cutoff);
    double direction = model_->solver()->getObjSense();
    cutoff *= direction;
    cutoff = CoinMin(cutoff, solutionValue);
    OsiSolverInterface * solver = cloneBut(2);
    // But reset bounds
    solver->setColLower(continuousSolver->getColLower());
    solver->setColUpper(continuousSolver->getColUpper());
    int numberColumns = solver->getNumCols();
    // Fixed
    double * fixed = new double [numberColumns];
    for (int i = 0; i < numberColumns; i++)
        fixed[i] = -COIN_DBL_MAX;
    int whichSolution[10];
    for (int i = 0; i < useNumber; i++)
        whichSolution[i] = i;
    for (int i = 0; i < useNumber; i++) {
        int k = whichSolution[i];
        const double * solution = model_->savedSolution(k);
        for (int j = 0; j < numberColumns; j++) {
            if (solver->isInteger(j)) {
                if (fixed[j] == -COIN_DBL_MAX)
                    fixed[j] = floor(solution[j] + 0.5);
                else if (fabs(fixed[j] - solution[j]) > 1.0e-7)
                    fixed[j] = COIN_DBL_MAX;
            }
        }
    }
    const double * colLower = solver->getColLower();
    for (int i = 0; i < numberColumns; i++) {
        if (solver->isInteger(i)) {
            double value = fixed[i];
            if (value != COIN_DBL_MAX) {
                if (when_ < 10) {
                    solver->setColLower(i, value);
                    solver->setColUpper(i, value);
                } else if (value == colLower[i]) {
                    solver->setColUpper(i, value);
                }
            }
        }
    }
    int returnCode = smallBranchAndBound(solver, numberNodes_, betterSolution,
                                         solutionValue,
                                         solutionValue, "CbcHeuristicCrossover");
    if (returnCode < 0)
        returnCode = 0; // returned on size
    if ((returnCode&2) != 0) {
        // could add cut
        returnCode &= ~2;
    }

    delete solver;
    return returnCode;
}
// update model
void CbcHeuristicCrossover::setModel(CbcModel * model)
{
    model_ = model;
    if (model) {
        for (int i = 0; i < 10; i++)
            random_[i] = model->randomNumberGenerator()->randomDouble();
    }
}


