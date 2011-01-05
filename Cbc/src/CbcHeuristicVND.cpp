// $Id$
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// edwin 12/5/09 carved out of CbcHeuristicRINS

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
#include "CbcHeuristicVND.hpp"
#include "CbcBranchActual.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"


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

