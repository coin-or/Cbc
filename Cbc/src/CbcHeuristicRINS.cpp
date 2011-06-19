/* $Id: CbcHeuristicRINS.cpp 1240 2009-10-02 18:41:44Z forrest $ */
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

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
    shallowDepth_ = 0;
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
    shallowDepth_ = 0;
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
	const int * used = model_->usedInSolution();
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
		/*
		  shallowDepth_
		  0 - normal
		  1 - only fix if at lb
		  2 - only fix if not at lb
		  3 - only fix if at lb and !used
		*/
		bool fix=false;
		switch (shallowDepth_) {
		case 0:
		  fix = true;
		  break;
		case 1:
		if (nearest==originalLower) 
		  fix = true;
		  break;
		case 2:
		if (nearest!=originalLower) 
		  fix = true;
		  break;
		case 3:
		if (nearest==originalLower && !used[iColumn]) 
		  fix = true;
		  break;
		}
		if (fix) {
		  newSolver->setColLower(iColumn, nearest);
		  newSolver->setColUpper(iColumn, nearest);
		  nFix++;
		}
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



