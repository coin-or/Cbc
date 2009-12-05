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
#include "CbcHeuristicRENS.hpp"
#include "CbcBranchActual.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"

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
