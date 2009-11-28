/* $Id: CbcHeuristicRINS.cpp 1240 2009-10-02 18:41:44Z forrest $ */
// Copyright (C) 2006, International Business Machines
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
#include "CbcHeuristicRINS.hpp"
#include "CbcBranchActual.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"

// Default Constructor
CbcHeuristicRINS::CbcHeuristicRINS()
        : CbcHeuristic()
{
    numberSolutions_ = 0;
    numberSuccesses_ = 0;
    numberTries_ = 0;
    stateOfFixing_ = 0;
    lastNode_ = -999999;
    howOften_ = 100;
    decayFactor_ = 0.5;
    used_ = NULL;
    whereFrom_ = 1 + 8 + 16 + 255 * 256;
    whereFrom_ = 1 + 8 + 255 * 256;
}

// Constructor with model - assumed before cuts

CbcHeuristicRINS::CbcHeuristicRINS(CbcModel & model)
        : CbcHeuristic(model)
{
    numberSolutions_ = 0;
    numberSuccesses_ = 0;
    numberTries_ = 0;
    stateOfFixing_ = 0;
    lastNode_ = -999999;
    howOften_ = 100;
    decayFactor_ = 0.5;
    assert(model.solver());
    int numberColumns = model.solver()->getNumCols();
    used_ = new char[numberColumns];
    memset(used_, 0, numberColumns);
    whereFrom_ = 1 + 8 + 16 + 255 * 256;
    whereFrom_ = 1 + 8 + 255 * 256;
}

// Destructor
CbcHeuristicRINS::~CbcHeuristicRINS ()
{
    delete [] used_;
}

// Clone
CbcHeuristic *
CbcHeuristicRINS::clone() const
{
    return new CbcHeuristicRINS(*this);
}

// Assignment operator
CbcHeuristicRINS &
CbcHeuristicRINS::operator=( const CbcHeuristicRINS & rhs)
{
    if (this != &rhs) {
        CbcHeuristic::operator=(rhs);
        numberSolutions_ = rhs.numberSolutions_;
        howOften_ = rhs.howOften_;
        numberSuccesses_ = rhs.numberSuccesses_;
        numberTries_ = rhs.numberTries_;
        stateOfFixing_ = rhs.stateOfFixing_;
        lastNode_ = rhs.lastNode_;
        delete [] used_;
        if (model_ && rhs.used_) {
            int numberColumns = model_->solver()->getNumCols();
            used_ = new char[numberColumns];
            memcpy(used_, rhs.used_, numberColumns);
        } else {
            used_ = NULL;
        }
    }
    return *this;
}

// Create C++ lines to get to current state
void
CbcHeuristicRINS::generateCpp( FILE * fp)
{
    CbcHeuristicRINS other;
    fprintf(fp, "0#include \"CbcHeuristicRINS.hpp\"\n");
    fprintf(fp, "3  CbcHeuristicRINS heuristicRINS(*cbcModel);\n");
    CbcHeuristic::generateCpp(fp, "heuristicRINS");
    if (howOften_ != other.howOften_)
        fprintf(fp, "3  heuristicRINS.setHowOften(%d);\n", howOften_);
    else
        fprintf(fp, "4  heuristicRINS.setHowOften(%d);\n", howOften_);
    fprintf(fp, "3  cbcModel->addHeuristic(&heuristicRINS);\n");
}

// Copy constructor
CbcHeuristicRINS::CbcHeuristicRINS(const CbcHeuristicRINS & rhs)
        :
        CbcHeuristic(rhs),
        numberSolutions_(rhs.numberSolutions_),
        howOften_(rhs.howOften_),
        numberSuccesses_(rhs.numberSuccesses_),
        numberTries_(rhs.numberTries_),
        stateOfFixing_(rhs.stateOfFixing_),
        lastNode_(rhs.lastNode_)
{
    if (model_ && rhs.used_) {
        int numberColumns = model_->solver()->getNumCols();
        used_ = new char[numberColumns];
        memcpy(used_, rhs.used_, numberColumns);
    } else {
        used_ = NULL;
    }
}
// Resets stuff if model changes
void
CbcHeuristicRINS::resetModel(CbcModel * /*model*/)
{
    //CbcHeuristic::resetModel(model);
    delete [] used_;
    stateOfFixing_ = 0;
    if (model_ && used_) {
        int numberColumns = model_->solver()->getNumCols();
        used_ = new char[numberColumns];
        memset(used_, 0, numberColumns);
    } else {
        used_ = NULL;
    }
}
/*
  First tries setting a variable to better value.  If feasible then
  tries setting others.  If not feasible then tries swaps
  Returns 1 if solution, 0 if not */
int
CbcHeuristicRINS::solution(double & solutionValue,
                           double * betterSolution)
{
    numCouldRun_++;
    int returnCode = 0;
    const double * bestSolution = model_->bestSolution();
    if (!bestSolution)
        return 0; // No solution found yet
    if (numberSolutions_ < model_->getSolutionCount()) {
        // new solution - add info
        numberSolutions_ = model_->getSolutionCount();

        int numberIntegers = model_->numberIntegers();
        const int * integerVariable = model_->integerVariable();

        int i;
        for (i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            const OsiObject * object = model_->object(i);
            // get original bounds
            double originalLower;
            double originalUpper;
            getIntegerInformation( object, originalLower, originalUpper);
            double value = bestSolution[iColumn];
            if (value < originalLower) {
                value = originalLower;
            } else if (value > originalUpper) {
                value = originalUpper;
            }
            double nearest = floor(value + 0.5);
            // if away from lower bound mark that fact
            if (nearest > originalLower) {
                used_[iColumn] = 1;
            }
        }
    }
    int numberNodes = model_->getNodeCount();
    if (howOften_ == 100) {
        if (numberNodes < lastNode_ + 12)
            return 0;
        // Do at 50 and 100
        if ((numberNodes > 40 && numberNodes <= 50) || (numberNodes > 90 && numberNodes < 100))
            numberNodes = howOften_;
    }
    // Allow for infeasible nodes - so do anyway after a bit
    if (howOften_ >= 100 && numberNodes >= lastNode_ + 2*howOften_) {
        numberNodes = howOften_;
    }
    if ((numberNodes % howOften_) == 0 && (model_->getCurrentPassNumber() == 1 ||
                                           model_->getCurrentPassNumber() == 999999)) {
        lastNode_ = model_->getNodeCount();
        OsiSolverInterface * solver = model_->solver();

        int numberIntegers = model_->numberIntegers();
        const int * integerVariable = model_->integerVariable();

        const double * currentSolution = solver->getColSolution();
        OsiSolverInterface * newSolver = cloneBut(3); // was model_->continuousSolver()->clone();
        int numberColumns = newSolver->getNumCols();
        int numberContinuous = numberColumns - numberIntegers;

        double primalTolerance;
        solver->getDblParam(OsiPrimalTolerance, primalTolerance);

        int i;
        int nFix = 0;
        for (i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            const OsiObject * object = model_->object(i);
            // get original bounds
            double originalLower;
            double originalUpper;
            getIntegerInformation( object, originalLower, originalUpper);
            double valueInt = bestSolution[iColumn];
            if (valueInt < originalLower) {
                valueInt = originalLower;
            } else if (valueInt > originalUpper) {
                valueInt = originalUpper;
            }
            if (fabs(currentSolution[iColumn] - valueInt) < 10.0*primalTolerance) {
                double nearest = floor(valueInt + 0.5);
                newSolver->setColLower(iColumn, nearest);
                newSolver->setColUpper(iColumn, nearest);
                nFix++;
            }
        }
        int divisor = 0;
        if (5*nFix > numberIntegers) {
            if (numberContinuous > 2*numberIntegers && ((nFix*10 < numberColumns &&
                    !numRuns_ && numberTries_ > 2) ||
                    stateOfFixing_)) {
#define RINS_FIX_CONTINUOUS
#ifdef RINS_FIX_CONTINUOUS
                const double * colLower = newSolver->getColLower();
                //const double * colUpper = newSolver->getColUpper();
                int nAtLb = 0;
                //double sumDj=0.0;
                const double * dj = newSolver->getReducedCost();
                double direction = newSolver->getObjSense();
                for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (!newSolver->isInteger(iColumn)) {
                        double value = bestSolution[iColumn];
                        if (value < colLower[iColumn] + 1.0e-8) {
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
                            double value = bestSolution[iColumn];
                            if (value < colLower[iColumn] + 1.0e-8) {
                                double djValue = dj[iColumn] * direction;
                                if (djValue > 1.0e-6) {
                                    sort[nFix2] = -djValue;
                                    which[nFix2++] = iColumn;
                                }
                            }
                        }
                    }
                    CoinSort_2(sort, sort + nFix2, which);
                    divisor = 4;
                    if (stateOfFixing_ > 0)
                        divisor = stateOfFixing_;
                    else if (stateOfFixing_ < -1)
                        divisor = (-stateOfFixing_) - 1;
                    nFix2 = CoinMin(nFix2, (numberColumns - nFix) / divisor);
                    for (int i = 0; i < nFix2; i++) {
                        int iColumn = which[i];
                        newSolver->setColUpper(iColumn, colLower[iColumn]);
                    }
                    delete [] sort;
                    delete [] which;
#ifdef CLP_INVESTIGATE2
                    printf("%d integers have same value, and %d continuous fixed at lb\n",
                           nFix, nFix2);
#endif
                }
#endif
            }
            //printf("%d integers have same value\n",nFix);
            returnCode = smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue,
                                             model_->getCutoff(), "CbcHeuristicRINS");
            if (returnCode < 0) {
                returnCode = 0; // returned on size
                if (divisor) {
                    stateOfFixing_ = - divisor; // say failed
                } else if (numberContinuous > 2*numberIntegers &&
                           !numRuns_ && numberTries_ > 2) {
                    stateOfFixing_ = -4; //start fixing
                }
            } else {
                numRuns_++;
                if (divisor)
                    stateOfFixing_ =  divisor; // say small enough
            }
            if ((returnCode&1) != 0)
                numberSuccesses_++;
            //printf("return code %d",returnCode);
            if ((returnCode&2) != 0) {
                // could add cut
                returnCode &= ~2;
                //printf("could add cut with %d elements (if all 0-1)\n",nFix);
            } else {
                //printf("\n");
            }
        }

        numberTries_++;
        if ((numberTries_ % 10) == 0 && numberSuccesses_*3 < numberTries_)
            howOften_ += static_cast<int> (howOften_ * decayFactor_);
        delete newSolver;
    }
    return returnCode;
}
// update model
void CbcHeuristicRINS::setModel(CbcModel * model)
{
    model_ = model;
    // Get a copy of original matrix
    assert(model_->solver());
    delete [] used_;
    int numberColumns = model->solver()->getNumCols();
    used_ = new char[numberColumns];
    memset(used_, 0, numberColumns);
}
// Default Constructor
CbcHeuristicRENS::CbcHeuristicRENS()
        : CbcHeuristic()
{
    numberTries_ = 0;
    whereFrom_ = 256 + 1;
}

// Constructor with model - assumed before cuts

CbcHeuristicRENS::CbcHeuristicRENS(CbcModel & model)
        : CbcHeuristic(model)
{
    numberTries_ = 0;
    whereFrom_ = 256 + 1;
}

// Destructor
CbcHeuristicRENS::~CbcHeuristicRENS ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicRENS::clone() const
{
    return new CbcHeuristicRENS(*this);
}

// Assignment operator
CbcHeuristicRENS &
CbcHeuristicRENS::operator=( const CbcHeuristicRENS & rhs)
{
    if (this != &rhs) {
        CbcHeuristic::operator=(rhs);
        numberTries_ = rhs.numberTries_;
    }
    return *this;
}

// Copy constructor
CbcHeuristicRENS::CbcHeuristicRENS(const CbcHeuristicRENS & rhs)
        :
        CbcHeuristic(rhs),
        numberTries_(rhs.numberTries_)
{
}
// Resets stuff if model changes
void
CbcHeuristicRENS::resetModel(CbcModel * )
{
}
int
CbcHeuristicRENS::solution(double & solutionValue,
                           double * betterSolution)
{
    int returnCode = 0;
    const double * bestSolution = model_->bestSolution();
    if (numberTries_ || (when() < 2 && bestSolution))
        return 0;
    numberTries_++;
    OsiSolverInterface * solver = model_->solver();

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();

    const double * currentSolution = solver->getColSolution();
    OsiSolverInterface * newSolver = cloneBut(3); // was model_->continuousSolver()->clone();
    const double * colLower = newSolver->getColLower();
    const double * colUpper = newSolver->getColUpper();

    double primalTolerance;
    solver->getDblParam(OsiPrimalTolerance, primalTolerance);

    int i;
    int numberFixed = 0;
    int numberTightened = 0;
    int numberAtBound = 0;
    int numberColumns = newSolver->getNumCols();
    int numberContinuous = numberColumns - numberIntegers;

    for (i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        double value = currentSolution[iColumn];
        double lower = colLower[iColumn];
        double upper = colUpper[iColumn];
        value = CoinMax(value, lower);
        value = CoinMin(value, upper);
#define RENS_FIX_ONLY_LOWER
#ifndef RENS_FIX_ONLY_LOWER
        if (fabs(value - floor(value + 0.5)) < 1.0e-8) {
            value = floor(value + 0.5);
            if (value == lower || value == upper)
                numberAtBound++;
            newSolver->setColLower(iColumn, value);
            newSolver->setColUpper(iColumn, value);
            numberFixed++;
        } else if (colUpper[iColumn] - colLower[iColumn] >= 2.0) {
            numberTightened++;
            newSolver->setColLower(iColumn, floor(value));
            newSolver->setColUpper(iColumn, ceil(value));
        }
#else
        if (fabs(value - floor(value + 0.5)) < 1.0e-8 &&
                floor(value + 0.5) == lower) {
            value = floor(value + 0.5);
            numberAtBound++;
            newSolver->setColLower(iColumn, value);
            newSolver->setColUpper(iColumn, value);
            numberFixed++;
        } else if (colUpper[iColumn] - colLower[iColumn] >= 2.0) {
            numberTightened++;
            if (fabs(value - floor(value + 0.5)) < 1.0e-8) {
                value = floor(value + 0.5);
                if (value < upper) {
                    newSolver->setColLower(iColumn, CoinMax(value - 1.0, lower));
                    newSolver->setColUpper(iColumn, CoinMin(value + 1.0, upper));
                } else {
                    newSolver->setColLower(iColumn, upper - 1.0);
                }
            } else {
                newSolver->setColLower(iColumn, floor(value));
                newSolver->setColUpper(iColumn, ceil(value));
            }
        }
#endif
    }
    if (numberFixed > numberIntegers / 5) {
        if (numberContinuous > numberIntegers && numberFixed < numberColumns / 5) {
#define RENS_FIX_CONTINUOUS
#ifdef RENS_FIX_CONTINUOUS
            const double * colLower = newSolver->getColLower();
            //const double * colUpper = newSolver->getColUpper();
            int nAtLb = 0;
            double sumDj = 0.0;
            const double * dj = newSolver->getReducedCost();
            double direction = newSolver->getObjSense();
            for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                if (!newSolver->isInteger(iColumn)) {
                    double value = currentSolution[iColumn];
                    if (value < colLower[iColumn] + 1.0e-8) {
                        double djValue = dj[iColumn] * direction;
                        nAtLb++;
                        sumDj += djValue;
                    }
                }
            }
            if (nAtLb) {
                // fix some continuous
                double * sort = new double[nAtLb];
                int * which = new int [nAtLb];
                double threshold = CoinMax((0.01 * sumDj) / static_cast<double>(nAtLb), 1.0e-6);
                int nFix2 = 0;
                for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (!newSolver->isInteger(iColumn)) {
                        double value = currentSolution[iColumn];
                        if (value < colLower[iColumn] + 1.0e-8) {
                            double djValue = dj[iColumn] * direction;
                            if (djValue > threshold) {
                                sort[nFix2] = -djValue;
                                which[nFix2++] = iColumn;
                            }
                        }
                    }
                }
                CoinSort_2(sort, sort + nFix2, which);
                nFix2 = CoinMin(nFix2, (numberColumns - numberFixed) / 2);
                for (int i = 0; i < nFix2; i++) {
                    int iColumn = which[i];
                    newSolver->setColUpper(iColumn, colLower[iColumn]);
                }
                delete [] sort;
                delete [] which;
#ifdef CLP_INVESTIGATE2
                printf("%d integers fixed (%d tightened) (%d at bound), and %d continuous fixed at lb\n",
                       numberFixed, numberTightened, numberAtBound, nFix2);
#endif
            }
#endif
        }
#ifdef COIN_DEVELOP
        printf("%d integers fixed and %d tightened\n", numberFixed, numberTightened);
#endif
        returnCode = smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue,
                                         model_->getCutoff(), "CbcHeuristicRENS");
        if (returnCode < 0) {
            returnCode = 0; // returned on size
#ifdef RENS_FIX_CONTINUOUS
            if (numberContinuous > numberIntegers && numberFixed >= numberColumns / 5) {
                const double * colLower = newSolver->getColLower();
                //const double * colUpper = newSolver->getColUpper();
                int nAtLb = 0;
                double sumDj = 0.0;
                const double * dj = newSolver->getReducedCost();
                double direction = newSolver->getObjSense();
                for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (!newSolver->isInteger(iColumn)) {
                        double value = currentSolution[iColumn];
                        if (value < colLower[iColumn] + 1.0e-8) {
                            double djValue = dj[iColumn] * direction;
                            nAtLb++;
                            sumDj += djValue;
                        }
                    }
                }
                if (nAtLb) {
                    // fix some continuous
                    double * sort = new double[nAtLb];
                    int * which = new int [nAtLb];
                    double threshold = CoinMax((0.01 * sumDj) / static_cast<double>(nAtLb), 1.0e-6);
                    int nFix2 = 0;
                    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                        if (!newSolver->isInteger(iColumn)) {
                            double value = currentSolution[iColumn];
                            if (value < colLower[iColumn] + 1.0e-8) {
                                double djValue = dj[iColumn] * direction;
                                if (djValue > threshold) {
                                    sort[nFix2] = -djValue;
                                    which[nFix2++] = iColumn;
                                }
                            }
                        }
                    }
                    CoinSort_2(sort, sort + nFix2, which);
                    nFix2 = CoinMin(nFix2, (numberColumns - numberFixed) / 2);
                    for (int i = 0; i < nFix2; i++) {
                        int iColumn = which[i];
                        newSolver->setColUpper(iColumn, colLower[iColumn]);
                    }
                    delete [] sort;
                    delete [] which;
#ifdef CLP_INVESTIGATE2
                    printf("%d integers fixed (%d tightened) (%d at bound), and %d continuous fixed at lb\n",
                           numberFixed, numberTightened, numberAtBound, nFix2);
#endif
                }
                returnCode = smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue,
                                                 model_->getCutoff(), "CbcHeuristicRENS");
#endif
            }
        }
        //printf("return code %d",returnCode);
        if ((returnCode&2) != 0) {
            // could add cut
            returnCode &= ~2;
#ifdef COIN_DEVELOP
            if (!numberTightened && numberFixed == numberAtBound)
                printf("could add cut with %d elements\n", numberFixed);
#endif
        } else {
            //printf("\n");
        }
    }

    delete newSolver;
    return returnCode;
}
// update model
void CbcHeuristicRENS::setModel(CbcModel * model)
{
    model_ = model;
}

// Default Constructor
CbcHeuristicDINS::CbcHeuristicDINS()
        : CbcHeuristic()
{
    numberSolutions_ = 0;
    numberSuccesses_ = 0;
    numberTries_ = 0;
    howOften_ = 100;
    decayFactor_ = 0.5;
    maximumKeepSolutions_ = 5;
    numberKeptSolutions_ = 0;
    numberIntegers_ = -1;
    localSpace_ = 10;
    values_ = NULL;
}

// Constructor with model - assumed before cuts

CbcHeuristicDINS::CbcHeuristicDINS(CbcModel & model)
        : CbcHeuristic(model)
{
    numberSolutions_ = 0;
    numberSuccesses_ = 0;
    numberTries_ = 0;
    howOften_ = 100;
    decayFactor_ = 0.5;
    assert(model.solver());
    maximumKeepSolutions_ = 5;
    numberKeptSolutions_ = 0;
    numberIntegers_ = -1;
    localSpace_ = 10;
    values_ = NULL;
}

// Destructor
CbcHeuristicDINS::~CbcHeuristicDINS ()
{
    for (int i = 0; i < numberKeptSolutions_; i++)
        delete [] values_[i];
    delete [] values_;
}

// Clone
CbcHeuristic *
CbcHeuristicDINS::clone() const
{
    return new CbcHeuristicDINS(*this);
}

// Assignment operator
CbcHeuristicDINS &
CbcHeuristicDINS::operator=( const CbcHeuristicDINS & rhs)
{
    if (this != &rhs) {
        CbcHeuristic::operator=(rhs);
        numberSolutions_ = rhs.numberSolutions_;
        howOften_ = rhs.howOften_;
        numberSuccesses_ = rhs.numberSuccesses_;
        numberTries_ = rhs.numberTries_;
        for (int i = 0; i < numberKeptSolutions_; i++)
            delete [] values_[i];
        delete [] values_;
        maximumKeepSolutions_ = rhs.maximumKeepSolutions_;
        numberKeptSolutions_ = rhs.numberKeptSolutions_;
        numberIntegers_ = rhs.numberIntegers_;
        localSpace_ = rhs.localSpace_;
        if (model_ && rhs.values_) {
            assert (numberIntegers_ >= 0);
            values_ = new int * [maximumKeepSolutions_];
            for (int i = 0; i < maximumKeepSolutions_; i++)
                values_[i] = CoinCopyOfArray(rhs.values_[i], numberIntegers_);
        } else {
            values_ = NULL;
        }
    }
    return *this;
}

// Create C++ lines to get to current state
void
CbcHeuristicDINS::generateCpp( FILE * fp)
{
    CbcHeuristicDINS other;
    fprintf(fp, "0#include \"CbcHeuristicDINS.hpp\"\n");
    fprintf(fp, "3  CbcHeuristicDINS heuristicDINS(*cbcModel);\n");
    CbcHeuristic::generateCpp(fp, "heuristicDINS");
    if (howOften_ != other.howOften_)
        fprintf(fp, "3  heuristicDINS.setHowOften(%d);\n", howOften_);
    else
        fprintf(fp, "4  heuristicDINS.setHowOften(%d);\n", howOften_);
    if (maximumKeepSolutions_ != other.maximumKeepSolutions_)
        fprintf(fp, "3  heuristicDINS.setMaximumKeep(%d);\n", maximumKeepSolutions_);
    else
        fprintf(fp, "4  heuristicDINS.setMaximumKeep(%d);\n", maximumKeepSolutions_);
    fprintf(fp, "3  cbcModel->addHeuristic(&heuristicDINS);\n");
}

// Copy constructor
CbcHeuristicDINS::CbcHeuristicDINS(const CbcHeuristicDINS & rhs)
        :
        CbcHeuristic(rhs),
        numberSolutions_(rhs.numberSolutions_),
        howOften_(rhs.howOften_),
        numberSuccesses_(rhs.numberSuccesses_),
        numberTries_(rhs.numberTries_),
        maximumKeepSolutions_(rhs.maximumKeepSolutions_),
        numberKeptSolutions_(rhs.numberKeptSolutions_),
        numberIntegers_(rhs.numberIntegers_),
        localSpace_(rhs.localSpace_)
{
    if (model_ && rhs.values_) {
        assert (numberIntegers_ >= 0);
        values_ = new int * [maximumKeepSolutions_];
        for (int i = 0; i < maximumKeepSolutions_; i++)
            values_[i] = CoinCopyOfArray(rhs.values_[i], numberIntegers_);
    } else {
        values_ = NULL;
    }
}
// Resets stuff if model changes
void
CbcHeuristicDINS::resetModel(CbcModel * )
{
    //CbcHeuristic::resetModel(model);
    for (int i = 0; i < numberKeptSolutions_; i++)
        delete [] values_[i];
    delete [] values_;
    numberKeptSolutions_ = 0;
    numberIntegers_ = -1;
    numberSolutions_ = 0;
    values_ = NULL;
}
/*
  First tries setting a variable to better value.  If feasible then
  tries setting others.  If not feasible then tries swaps
  Returns 1 if solution, 0 if not */
int
CbcHeuristicDINS::solution(double & solutionValue,
                           double * betterSolution)
{
    numCouldRun_++;
    int returnCode = 0;
    const double * bestSolution = model_->bestSolution();
    if (!bestSolution)
        return 0; // No solution found yet
    if (numberSolutions_ < model_->getSolutionCount()) {
        // new solution - add info
        numberSolutions_ = model_->getSolutionCount();

        int numberIntegers = model_->numberIntegers();
        const int * integerVariable = model_->integerVariable();
        if (numberIntegers_ < 0) {
            numberIntegers_ = numberIntegers;
            assert (!values_);
            values_ = new int * [maximumKeepSolutions_];
            for (int i = 0; i < maximumKeepSolutions_; i++)
                values_[i] = NULL;
        } else {
            assert (numberIntegers == numberIntegers_);
        }
        // move solutions (0 will be most recent)
        {
            int * temp = values_[maximumKeepSolutions_-1];
            for (int i = maximumKeepSolutions_ - 1; i > 0; i--)
                values_[i] = values_[i-1];
            if (!temp)
                temp = new int [numberIntegers_];
            values_[0] = temp;
        }
        int i;
        for (i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            double value = bestSolution[iColumn];
            double nearest = floor(value + 0.5);
            values_[0][i] = static_cast<int> (nearest);
        }
        numberKeptSolutions_ = CoinMin(numberKeptSolutions_ + 1, maximumKeepSolutions_);
    }
    int finalReturnCode = 0;
    if (((model_->getNodeCount() % howOften_) == howOften_ / 2 || !model_->getNodeCount()) && (model_->getCurrentPassNumber() == 1 || model_->getCurrentPassNumber() == 999999)) {
        OsiSolverInterface * solver = model_->solver();

        int numberIntegers = model_->numberIntegers();
        const int * integerVariable = model_->integerVariable();

        const double * currentSolution = solver->getColSolution();
        int localSpace = localSpace_;
        // 0 means finished but no solution, 1 solution, 2 node limit
        int status = -1;
        double cutoff = model_->getCutoff();
        while (status) {
            status = 0;
            OsiSolverInterface * newSolver = cloneBut(3); // was model_->continuousSolver()->clone();
            const double * colLower = solver->getColLower();
            const double * colUpper = solver->getColUpper();

            double primalTolerance;
            solver->getDblParam(OsiPrimalTolerance, primalTolerance);
            const double * continuousSolution = newSolver->getColSolution();
            // Space for added constraint
            double * element = new double [numberIntegers];
            int * column = new int [numberIntegers];
            int i;
            int nFix = 0;
            int nCouldFix = 0;
            int nCouldFix2 = 0;
            int nBound = 0;
            int nEl = 0;
            double bias = localSpace;
            int okSame = numberKeptSolutions_ - 1;
            for (i = 0; i < numberIntegers; i++) {
                int iColumn = integerVariable[i];
                const OsiObject * object = model_->object(i);
                // get original bounds
                double originalLower;
                double originalUpper;
                getIntegerInformation( object, originalLower, originalUpper);
                double valueInt = bestSolution[iColumn];
                if (valueInt < originalLower) {
                    valueInt = originalLower;
                } else if (valueInt > originalUpper) {
                    valueInt = originalUpper;
                }
                int intValue = static_cast<int> (floor(valueInt + 0.5));
                double currentValue = currentSolution[iColumn];
                double currentLower = colLower[iColumn];
                double currentUpper = colUpper[iColumn];
                if (fabs(valueInt - currentValue) >= 0.5) {
                    // Re-bound
                    nBound++;
                    if (intValue >= currentValue) {
                        currentLower = CoinMax(currentLower, ceil(2 * currentValue - intValue));
                        currentUpper = intValue;
                    } else {
                        currentLower = intValue;
                        currentUpper = CoinMin(currentUpper, floor(2 * currentValue - intValue));
                    }
                    newSolver->setColLower(iColumn, currentLower);
                    newSolver->setColUpper(iColumn, currentUpper);
                } else {
                    // See if can fix
                    bool canFix = false;
                    double continuousValue = continuousSolution[iColumn];
                    if (fabs(currentValue - valueInt) < 10.0*primalTolerance) {
                        if (currentUpper - currentLower > 1.0) {
                            // General integer variable
                            canFix = true;
                        } else if (fabs(continuousValue - valueInt) < 10.0*primalTolerance) {
                            int nSame = 1;
                            //assert (intValue==values_[0][i]);
                            for (int k = 1; k < numberKeptSolutions_; k++) {
                                if (intValue == values_[k][i])
                                    nSame++;
                            }
                            if (nSame >= okSame) {
                                // can fix
                                canFix = true;
                            } else {
                                nCouldFix++;
                            }
                        } else {
                            nCouldFix2++;
                        }
                    }
                    if (canFix) {
                        newSolver->setColLower(iColumn, intValue);
                        newSolver->setColUpper(iColumn, intValue);
                        nFix++;
                    } else {
                        if (currentUpper - currentLower > 1.0) {
                            // General integer variable
                            currentLower = floor(currentValue);
                            if (intValue >= currentLower && intValue <= currentLower + 1) {
                                newSolver->setColLower(iColumn, currentLower);
                                newSolver->setColUpper(iColumn, currentLower + 1.0);
                            } else {
                                // fix
                                double value;
                                if (intValue < currentLower)
                                    value = currentLower;
                                else
                                    value = currentLower + 1;
                                newSolver->setColLower(iColumn, value);
                                newSolver->setColUpper(iColumn, value);
                                nFix++;
                            }
                        } else {
                            // 0-1 (ish)
                            column[nEl] = iColumn;
                            if (intValue == currentLower) {
                                bias += currentLower;
                                element[nEl++] = 1.0;
                            } else if (intValue == currentUpper) {
                                bias += currentUpper;
                                element[nEl++] = -1.0;
                            } else {
                                printf("bad DINS logic\n");
                                abort();
                            }
                        }
                    }
                }
            }
            char generalPrint[200];
            sprintf(generalPrint,
                    "%d fixed, %d same as cont/int, %d same as int - %d bounded %d in cut\n",
                    nFix, nCouldFix, nCouldFix2, nBound, nEl);
            model_->messageHandler()->message(CBC_FPUMP2, model_->messages())
            << generalPrint
            << CoinMessageEol;
            if (nFix > numberIntegers / 10) {
#if 0
                newSolver->initialSolve();
                printf("obj %g\n", newSolver->getObjValue());
                for (i = 0; i < numberIntegers; i++) {
                    int iColumn = integerVariable[i];
                    printf("%d new bounds %g %g - solutions %g %g\n",
                           iColumn, newSolver->getColLower()[iColumn],
                           newSolver->getColUpper()[iColumn],
                           bestSolution[iColumn],
                           currentSolution[iColumn]);
                }
#endif
                if (nEl > 0)
                    newSolver->addRow(nEl, column, element, -COIN_DBL_MAX, bias);
                //printf("%d integers have same value\n",nFix);
                returnCode = smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue,
                                                 cutoff, "CbcHeuristicDINS");
                if (returnCode < 0) {
                    returnCode = 0; // returned on size
                    status = 0;
                } else {
                    numRuns_++;
                    if ((returnCode&2) != 0) {
                        // could add cut as complete search
                        returnCode &= ~2;
                        if ((returnCode&1) != 0) {
                            numberSuccesses_++;
                            status = 1;
                        } else {
                            // no solution
                            status = 0;
                        }
                    } else {
                        if ((returnCode&1) != 0) {
                            numberSuccesses_++;
                            status = 1;
                        } else {
                            // no solution but node limit
                            status = 2;
                            if (nEl)
                                localSpace -= 5;
                            else
                                localSpace = -1;
                            if (localSpace < 0)
                                status = 0;
                        }
                    }
                    if ((returnCode&1) != 0) {
                        cutoff = CoinMin(cutoff, solutionValue - model_->getCutoffIncrement());
                        finalReturnCode = 1;
                    }
                }
            }
            delete [] element;
            delete [] column;
            delete newSolver;
        }
        numberTries_++;
        if ((numberTries_ % 10) == 0 && numberSuccesses_*3 < numberTries_)
            howOften_ += static_cast<int> (howOften_ * decayFactor_);
    }
    return finalReturnCode;
}
// update model
void CbcHeuristicDINS::setModel(CbcModel * model)
{
    model_ = model;
    // Get a copy of original matrix
    assert(model_->solver());
    for (int i = 0; i < numberKeptSolutions_; i++)
        delete [] values_[i];
    delete [] values_;
    numberKeptSolutions_ = 0;
    numberIntegers_ = -1;
    numberSolutions_ = 0;
    values_ = NULL;
}

// Default Constructor
CbcHeuristicVND::CbcHeuristicVND()
        : CbcHeuristic()
{
    numberSolutions_ = 0;
    numberSuccesses_ = 0;
    numberTries_ = 0;
    lastNode_ = -999999;
    howOften_ = 100;
    decayFactor_ = 0.5;
    baseSolution_ = NULL;
    whereFrom_ = 1 + 8 + 255 * 256;
    stepSize_ = 0;
    k_ = 0;
    kmax_ = 0;
    nDifferent_ = 0;
}

// Constructor with model - assumed before cuts

CbcHeuristicVND::CbcHeuristicVND(CbcModel & model)
        : CbcHeuristic(model)
{
    numberSolutions_ = 0;
    numberSuccesses_ = 0;
    numberTries_ = 0;
    lastNode_ = -999999;
    howOften_ = 100;
    decayFactor_ = 0.5;
    assert(model.solver());
    int numberColumns = model.solver()->getNumCols();
    baseSolution_ = new double [numberColumns];
    memset(baseSolution_, 0, numberColumns*sizeof(double));
    whereFrom_ = 1 + 8 + 255 * 256;
    stepSize_ = 0;
    k_ = 0;
    kmax_ = 0;
    nDifferent_ = 0;
}

// Destructor
CbcHeuristicVND::~CbcHeuristicVND ()
{
    delete [] baseSolution_;
}

// Clone
CbcHeuristic *
CbcHeuristicVND::clone() const
{
    return new CbcHeuristicVND(*this);
}

// Assignment operator
CbcHeuristicVND &
CbcHeuristicVND::operator=( const CbcHeuristicVND & rhs)
{
    if (this != &rhs) {
        CbcHeuristic::operator=(rhs);
        numberSolutions_ = rhs.numberSolutions_;
        howOften_ = rhs.howOften_;
        numberSuccesses_ = rhs.numberSuccesses_;
        numberTries_ = rhs.numberTries_;
        lastNode_ = rhs.lastNode_;
        delete [] baseSolution_;
        if (model_ && rhs.baseSolution_) {
            int numberColumns = model_->solver()->getNumCols();
            baseSolution_ = new double [numberColumns];
            memcpy(baseSolution_, rhs.baseSolution_, numberColumns*sizeof(double));
        } else {
            baseSolution_ = NULL;
        }
        stepSize_ = rhs.stepSize_;
        k_ = rhs.k_;
        kmax_ = rhs.kmax_;
        nDifferent_ = rhs.nDifferent_;
    }
    return *this;
}

// Create C++ lines to get to current state
void
CbcHeuristicVND::generateCpp( FILE * fp)
{
    CbcHeuristicVND other;
    fprintf(fp, "0#include \"CbcHeuristicVND.hpp\"\n");
    fprintf(fp, "3  CbcHeuristicVND heuristicVND(*cbcModel);\n");
    CbcHeuristic::generateCpp(fp, "heuristicVND");
    if (howOften_ != other.howOften_)
        fprintf(fp, "3  heuristicVND.setHowOften(%d);\n", howOften_);
    else
        fprintf(fp, "4  heuristicVND.setHowOften(%d);\n", howOften_);
    fprintf(fp, "3  cbcModel->addHeuristic(&heuristicVND);\n");
}

// Copy constructor
CbcHeuristicVND::CbcHeuristicVND(const CbcHeuristicVND & rhs)
        :
        CbcHeuristic(rhs),
        numberSolutions_(rhs.numberSolutions_),
        howOften_(rhs.howOften_),
        numberSuccesses_(rhs.numberSuccesses_),
        numberTries_(rhs.numberTries_),
        lastNode_(rhs.lastNode_)
{
    if (model_ && rhs.baseSolution_) {
        int numberColumns = model_->solver()->getNumCols();
        baseSolution_ = new double [numberColumns];
        memcpy(baseSolution_, rhs.baseSolution_, numberColumns*sizeof(double));
    } else {
        baseSolution_ = NULL;
    }
    stepSize_ = rhs.stepSize_;
    k_ = rhs.k_;
    kmax_ = rhs.kmax_;
    nDifferent_ = rhs.nDifferent_;
}
// Resets stuff if model changes
void
CbcHeuristicVND::resetModel(CbcModel * /*model*/)
{
    //CbcHeuristic::resetModel(model);
    delete [] baseSolution_;
    if (model_ && baseSolution_) {
        int numberColumns = model_->solver()->getNumCols();
        baseSolution_ = new double [numberColumns];
        memset(baseSolution_, 0, numberColumns*sizeof(double));
    } else {
        baseSolution_ = NULL;
    }
}
/*
  First tries setting a variable to better value.  If feasible then
  tries setting others.  If not feasible then tries swaps
  Returns 1 if solution, 0 if not */
int
CbcHeuristicVND::solution(double & solutionValue,
                          double * betterSolution)
{
    numCouldRun_++;
    int returnCode = 0;
    const double * bestSolution = model_->bestSolution();
    if (!bestSolution)
        return 0; // No solution found yet
    if (numberSolutions_ < model_->getSolutionCount()) {
        // new solution - add info
        numberSolutions_ = model_->getSolutionCount();

        int numberIntegers = model_->numberIntegers();
        const int * integerVariable = model_->integerVariable();

        int i;
        for (i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            const OsiObject * object = model_->object(i);
            // get original bounds
            double originalLower;
            double originalUpper;
            getIntegerInformation( object, originalLower, originalUpper);
            double value = bestSolution[iColumn];
            if (value < originalLower) {
                value = originalLower;
            } else if (value > originalUpper) {
                value = originalUpper;
            }
        }
    }
    int numberNodes = model_->getNodeCount();
    if (howOften_ == 100) {
        if (numberNodes < lastNode_ + 12)
            return 0;
        // Do at 50 and 100
        if ((numberNodes > 40 && numberNodes <= 50) || (numberNodes > 90 && numberNodes < 100))
            numberNodes = howOften_;
    }
    if ((numberNodes % howOften_) == 0 && (model_->getCurrentPassNumber() == 1 ||
                                           model_->getCurrentPassNumber() == 999999)) {
        lastNode_ = model_->getNodeCount();
        OsiSolverInterface * solver = model_->solver();

        int numberIntegers = model_->numberIntegers();
        const int * integerVariable = model_->integerVariable();

        const double * currentSolution = solver->getColSolution();
        OsiSolverInterface * newSolver = cloneBut(3); // was model_->continuousSolver()->clone();
        //const double * colLower = newSolver->getColLower();
        //const double * colUpper = newSolver->getColUpper();

        double primalTolerance;
        solver->getDblParam(OsiPrimalTolerance, primalTolerance);

        // Sort on distance
        double * distance = new double [numberIntegers];
        int * which = new int [numberIntegers];

        int i;
        int nFix = 0;
        double tolerance = 10.0 * primalTolerance;
        for (i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            const OsiObject * object = model_->object(i);
            // get original bounds
            double originalLower;
            double originalUpper;
            getIntegerInformation( object, originalLower, originalUpper);
            double valueInt = bestSolution[iColumn];
            if (valueInt < originalLower) {
                valueInt = originalLower;
            } else if (valueInt > originalUpper) {
                valueInt = originalUpper;
            }
            baseSolution_[iColumn] = currentSolution[iColumn];
            distance[i] = fabs(currentSolution[iColumn] - valueInt);
            which[i] = i;
            if (fabs(currentSolution[iColumn] - valueInt) < tolerance)
                nFix++;
        }
        CoinSort_2(distance, distance + numberIntegers, which);
        nDifferent_ = numberIntegers - nFix;
        stepSize_ = nDifferent_ / 10;
        k_ = stepSize_;
        //nFix = numberIntegers-stepSize_;
        for (i = 0; i < nFix; i++) {
            int j = which[i];
            int iColumn = integerVariable[j];
            const OsiObject * object = model_->object(i);
            // get original bounds
            double originalLower;
            double originalUpper;
            getIntegerInformation( object, originalLower, originalUpper);
            double valueInt = bestSolution[iColumn];
            if (valueInt < originalLower) {
                valueInt = originalLower;
            } else if (valueInt > originalUpper) {
                valueInt = originalUpper;
            }
            double nearest = floor(valueInt + 0.5);
            newSolver->setColLower(iColumn, nearest);
            newSolver->setColUpper(iColumn, nearest);
        }
        delete [] distance;
        delete [] which;
        if (nFix > numberIntegers / 5) {
            //printf("%d integers have samish value\n",nFix);
            returnCode = smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue,
                                             model_->getCutoff(), "CbcHeuristicVND");
            if (returnCode < 0)
                returnCode = 0; // returned on size
            else
                numRuns_++;
            if ((returnCode&1) != 0)
                numberSuccesses_++;
            //printf("return code %d",returnCode);
            if ((returnCode&2) != 0) {
                // could add cut
                returnCode &= ~2;
                //printf("could add cut with %d elements (if all 0-1)\n",nFix);
            } else {
                //printf("\n");
            }
            numberTries_++;
            if ((numberTries_ % 10) == 0 && numberSuccesses_*3 < numberTries_)
                howOften_ += static_cast<int> (howOften_ * decayFactor_);
        }

        delete newSolver;
    }
    return returnCode;
}
// update model
void CbcHeuristicVND::setModel(CbcModel * model)
{
    model_ = model;
    // Get a copy of original matrix
    assert(model_->solver());
    delete [] baseSolution_;
    int numberColumns = model->solver()->getNumCols();
    baseSolution_ = new double [numberColumns];
    memset(baseSolution_, 0, numberColumns*sizeof(double));
}


