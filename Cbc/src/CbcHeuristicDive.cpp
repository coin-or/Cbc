/* $Id: CbcHeuristicDive.cpp 1240 2009-10-02 18:41:44Z forrest $ */
// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CbcHeuristicDive.hpp"
#include "CbcStrategy.hpp"
#include "OsiAuxInfo.hpp"
#include  "CoinTime.hpp"

#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif

//#define DIVE_FIX_BINARY_VARIABLES
//#define DIVE_DEBUG

// Default Constructor
CbcHeuristicDive::CbcHeuristicDive()
        : CbcHeuristic()
{
    // matrix and row copy will automatically be empty
    downLocks_ = NULL;
    upLocks_ = NULL;
    downArray_ = NULL;
    upArray_ = NULL;
    percentageToFix_ = 0.2;
    maxIterations_ = 100;
    maxSimplexIterations_ = 10000;
    maxSimplexIterationsAtRoot_ = 1000000;
    maxTime_ = 600;
    whereFrom_ = 255 - 2 - 16 + 256;
    decayFactor_ = 1.0;
}

// Constructor from model
CbcHeuristicDive::CbcHeuristicDive(CbcModel & model)
        : CbcHeuristic(model)
{
    downLocks_ = NULL;
    upLocks_ = NULL;
    downArray_ = NULL;
    upArray_ = NULL;
    // Get a copy of original matrix
    assert(model.solver());
    // model may have empty matrix - wait until setModel
    const CoinPackedMatrix * matrix = model.solver()->getMatrixByCol();
    if (matrix) {
        matrix_ = *matrix;
        matrixByRow_ = *model.solver()->getMatrixByRow();
        validate();
    }
    percentageToFix_ = 0.2;
    maxIterations_ = 100;
    maxSimplexIterations_ = 10000;
    maxSimplexIterationsAtRoot_ = 1000000;
    maxTime_ = 600;
    whereFrom_ = 255 - 2 - 16 + 256;
    decayFactor_ = 1.0;
}

// Destructor
CbcHeuristicDive::~CbcHeuristicDive ()
{
    delete [] downLocks_;
    delete [] upLocks_;
    assert (!downArray_);
}

// Create C++ lines to get to current state
void
CbcHeuristicDive::generateCpp( FILE * fp, const char * heuristic)
{
    // hard coded as CbcHeuristic virtual
    CbcHeuristic::generateCpp(fp, heuristic);
    if (percentageToFix_ != 0.2)
        fprintf(fp, "3  %s.setPercentageToFix(%.f);\n", heuristic, percentageToFix_);
    else
        fprintf(fp, "4  %s.setPercentageToFix(%.f);\n", heuristic, percentageToFix_);
    if (maxIterations_ != 100)
        fprintf(fp, "3  %s.setMaxIterations(%d);\n", heuristic, maxIterations_);
    else
        fprintf(fp, "4  %s.setMaxIterations(%d);\n", heuristic, maxIterations_);
    if (maxSimplexIterations_ != 10000)
        fprintf(fp, "3  %s.setMaxSimplexIterations(%d);\n", heuristic, maxSimplexIterations_);
    else
        fprintf(fp, "4  %s.setMaxSimplexIterations(%d);\n", heuristic, maxSimplexIterations_);
    if (maxTime_ != 600)
        fprintf(fp, "3  %s.setMaxTime(%.2f);\n", heuristic, maxTime_);
    else
        fprintf(fp, "4  %s.setMaxTime(%.2f);\n", heuristic, maxTime_);
}

// Copy constructor
CbcHeuristicDive::CbcHeuristicDive(const CbcHeuristicDive & rhs)
        :
        CbcHeuristic(rhs),
        matrix_(rhs.matrix_),
        matrixByRow_(rhs.matrixByRow_),
        percentageToFix_(rhs.percentageToFix_),
        maxIterations_(rhs.maxIterations_),
        maxSimplexIterations_(rhs.maxSimplexIterations_),
        maxSimplexIterationsAtRoot_(rhs.maxSimplexIterationsAtRoot_),
        maxTime_(rhs.maxTime_)
{
    downArray_ = NULL;
    upArray_ = NULL;
    if (rhs.downLocks_) {
        int numberIntegers = model_->numberIntegers();
        downLocks_ = CoinCopyOfArray(rhs.downLocks_, numberIntegers);
        upLocks_ = CoinCopyOfArray(rhs.upLocks_, numberIntegers);
    } else {
        downLocks_ = NULL;
        upLocks_ = NULL;
    }
}

// Assignment operator
CbcHeuristicDive &
CbcHeuristicDive::operator=( const CbcHeuristicDive & rhs)
{
    if (this != &rhs) {
        CbcHeuristic::operator=(rhs);
        matrix_ = rhs.matrix_;
        matrixByRow_ = rhs.matrixByRow_;
        percentageToFix_ = rhs.percentageToFix_;
        maxIterations_ = rhs.maxIterations_;
        maxSimplexIterations_ = rhs.maxSimplexIterations_;
        maxSimplexIterationsAtRoot_ = rhs.maxSimplexIterationsAtRoot_;
        maxTime_ = rhs.maxTime_;
        delete [] downLocks_;
        delete [] upLocks_;
        if (rhs.downLocks_) {
            int numberIntegers = model_->numberIntegers();
            downLocks_ = CoinCopyOfArray(rhs.downLocks_, numberIntegers);
            upLocks_ = CoinCopyOfArray(rhs.upLocks_, numberIntegers);
        } else {
            downLocks_ = NULL;
            upLocks_ = NULL;
        }
    }
    return *this;
}

// Resets stuff if model changes
void
CbcHeuristicDive::resetModel(CbcModel * model)
{
    model_ = model;
    assert(model_->solver());
    // Get a copy of original matrix
    const CoinPackedMatrix * matrix = model_->solver()->getMatrixByCol();
    // model may have empty matrix - wait until setModel
    if (matrix) {
        matrix_ = *matrix;
        matrixByRow_ = *model->solver()->getMatrixByRow();
        validate();
    }
}

// update model
void CbcHeuristicDive::setModel(CbcModel * model)
{
    model_ = model;
    assert(model_->solver());
    // Get a copy of original matrix
    const CoinPackedMatrix * matrix = model_->solver()->getMatrixByCol();
    if (matrix) {
        matrix_ = *matrix;
        matrixByRow_ = *model->solver()->getMatrixByRow();
        // make sure model okay for heuristic
        validate();
    }
}

bool CbcHeuristicDive::canHeuristicRun()
{
    return shouldHeurRun_randomChoice();
}

inline bool compareBinaryVars(const PseudoReducedCost obj1,
                              const PseudoReducedCost obj2)
{
    return obj1.pseudoRedCost > obj2.pseudoRedCost;
}

// See if diving will give better solution
// Sets value of solution
// Returns 1 if solution, 0 if not
int
CbcHeuristicDive::solution(double & solutionValue,
                           double * betterSolution)
{
    int nodeCount = model_->getNodeCount();
    if (feasibilityPumpOptions_>0 && (nodeCount % feasibilityPumpOptions_) != 0)
        return 0;
#ifdef DIVE_DEBUG
    std::cout << "solutionValue = " << solutionValue << std::endl;
#endif
    ++numCouldRun_;

    // test if the heuristic can run
    if (!canHeuristicRun())
        return 0;

#ifdef JJF_ZERO
    // See if to do
    if (!when() || (when() % 10 == 1 && model_->phase() != 1) ||
            (when() % 10 == 2 && (model_->phase() != 2 && model_->phase() != 3)))
        return 0; // switched off
#endif

#ifdef DIVE_DEBUG
    int nRoundInfeasible = 0;
    int nRoundFeasible = 0;
    int reasonToStop = 0;
#endif
    double time1 = CoinCpuTime();
    int numberSimplexIterations = 0;
    int maxSimplexIterations = (model_->getNodeCount()) ? maxSimplexIterations_
                               : maxSimplexIterationsAtRoot_;

    OsiSolverInterface * solver = cloneBut(6); // was model_->solver()->clone();
# ifdef COIN_HAS_CLP
    OsiClpSolverInterface * clpSolver
    = dynamic_cast<OsiClpSolverInterface *> (solver);
    if (clpSolver) {
        // say give up easily
        ClpSimplex * clpSimplex = clpSolver->getModelPtr();
        clpSimplex->setMoreSpecialOptions(clpSimplex->moreSpecialOptions() | 64);
    }
# endif
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    const double * rowLower = solver->getRowLower();
    const double * rowUpper = solver->getRowUpper();
    const double * solution = solver->getColSolution();
    const double * objective = solver->getObjCoefficients();
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double primalTolerance;
    solver->getDblParam(OsiPrimalTolerance, primalTolerance);

    int numberRows = matrix_.getNumRows();
    assert (numberRows <= solver->getNumRows());
    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();
    double direction = solver->getObjSense(); // 1 for min, -1 for max
    double newSolutionValue = direction * solver->getObjValue();
    int returnCode = 0;
    // Column copy
    const double * element = matrix_.getElements();
    const int * row = matrix_.getIndices();
    const CoinBigIndex * columnStart = matrix_.getVectorStarts();
    const int * columnLength = matrix_.getVectorLengths();
#ifdef DIVE_FIX_BINARY_VARIABLES
    // Row copy
    const double * elementByRow = matrixByRow_.getElements();
    const int * column = matrixByRow_.getIndices();
    const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
    const int * rowLength = matrixByRow_.getVectorLengths();
#endif

    // Get solution array for heuristic solution
    int numberColumns = solver->getNumCols();
    double * newSolution = new double [numberColumns];
    memcpy(newSolution, solution, numberColumns*sizeof(double));

    // vectors to store the latest variables fixed at their bounds
    int* columnFixed = new int [numberIntegers];
    double* originalBound = new double [numberIntegers];
    bool * fixedAtLowerBound = new bool [numberIntegers];
    PseudoReducedCost * candidate = new PseudoReducedCost [numberIntegers];
    double * random = new double [numberIntegers];

    int maxNumberAtBoundToFix = static_cast<int> (floor(percentageToFix_ * numberIntegers));

    // count how many fractional variables
    int numberFractionalVariables = 0;
    for (int i = 0; i < numberIntegers; i++) {
        random[i] = randomNumberGenerator_.randomDouble() + 0.3;
        int iColumn = integerVariable[i];
        double value = newSolution[iColumn];
        if (fabs(floor(value + 0.5) - value) > integerTolerance) {
            numberFractionalVariables++;
        }
    }

    const double* reducedCost = NULL;
    // See if not NLP
    if (model_->solverCharacteristics()->reducedCostsAccurate())
        reducedCost = solver->getReducedCost();

    int iteration = 0;
    while (numberFractionalVariables) {
        iteration++;

        // initialize any data
        initializeData();

        // select a fractional variable to bound
        int bestColumn = -1;
        int bestRound; // -1 rounds down, +1 rounds up
        bool canRound = selectVariableToBranch(solver, newSolution,
                                               bestColumn, bestRound);

        // if the solution is not trivially roundable, we don't try to round;
        // if the solution is trivially roundable, we try to round. However,
        // if the rounded solution is worse than the current incumbent,
        // then we don't round and proceed normally. In this case, the
        // bestColumn will be a trivially roundable variable
        if (canRound) {
            // check if by rounding all fractional variables
            // we get a solution with an objective value
            // better than the current best integer solution
            double delta = 0.0;
            for (int i = 0; i < numberIntegers; i++) {
                int iColumn = integerVariable[i];
                double value = newSolution[iColumn];
                if (fabs(floor(value + 0.5) - value) > integerTolerance) {
                    assert(downLocks_[i] == 0 || upLocks_[i] == 0);
                    double obj = objective[iColumn];
                    if (downLocks_[i] == 0 && upLocks_[i] == 0) {
                        if (direction * obj >= 0.0)
                            delta += (floor(value) - value) * obj;
                        else
                            delta += (ceil(value) - value) * obj;
                    } else if (downLocks_[i] == 0)
                        delta += (floor(value) - value) * obj;
                    else
                        delta += (ceil(value) - value) * obj;
                }
            }
            if (direction*(solver->getObjValue() + delta) < solutionValue) {
#ifdef DIVE_DEBUG
                nRoundFeasible++;
#endif
                // Round all the fractional variables
                for (int i = 0; i < numberIntegers; i++) {
                    int iColumn = integerVariable[i];
                    double value = newSolution[iColumn];
                    if (fabs(floor(value + 0.5) - value) > integerTolerance) {
                        assert(downLocks_[i] == 0 || upLocks_[i] == 0);
                        if (downLocks_[i] == 0 && upLocks_[i] == 0) {
                            if (direction * objective[iColumn] >= 0.0)
                                newSolution[iColumn] = floor(value);
                            else
                                newSolution[iColumn] = ceil(value);
                        } else if (downLocks_[i] == 0)
                            newSolution[iColumn] = floor(value);
                        else
                            newSolution[iColumn] = ceil(value);
                    }
                }
                break;
            }
#ifdef DIVE_DEBUG
            else
                nRoundInfeasible++;
#endif
        }

        // do reduced cost fixing
#ifdef DIVE_DEBUG
        int numberFixed = reducedCostFix(solver);
        std::cout << "numberReducedCostFixed = " << numberFixed << std::endl;
#else
        reducedCostFix(solver);
#endif

        int numberAtBoundFixed = 0;
#ifdef DIVE_FIX_BINARY_VARIABLES
        // fix binary variables based on pseudo reduced cost
        if (binVarIndex_.size()) {
            int cnt = 0;
            int n = static_cast<int>(binVarIndex_.size());
            for (int j = 0; j < n; j++) {
                int iColumn1 = binVarIndex_[j];
                double value = newSolution[iColumn1];
                if (fabs(value) <= integerTolerance &&
                        lower[iColumn1] != upper[iColumn1]) {
                    double maxPseudoReducedCost = 0.0;
#ifdef DIVE_DEBUG
                    std::cout << "iColumn1 = " << iColumn1 << ", value = " << value << std::endl;
#endif
                    int iRow = vbRowIndex_[j];
                    double chosenValue = 0.0;
                    for (int k = rowStart[iRow]; k < rowStart[iRow] + rowLength[iRow]; k++) {
                        int iColumn2 = column[k];
#ifdef DIVE_DEBUG
                        std::cout << "iColumn2 = " << iColumn2 << std::endl;
#endif
                        if (iColumn1 != iColumn2) {
                            double pseudoReducedCost = fabs(reducedCost[iColumn2] *
                                                            elementByRow[k]);
#ifdef DIVE_DEBUG
                            int k2;
                            for (k2 = rowStart[iRow]; k2 < rowStart[iRow] + rowLength[iRow]; k2++) {
                                if (column[k2] == iColumn1)
                                    break;
                            }
                            std::cout << "reducedCost[" << iColumn2 << "] = "
                                      << reducedCost[iColumn2]
                                      << ", elementByRow[" << iColumn2 << "] = " << elementByRow[k]
                                      << ", elementByRow[" << iColumn1 << "] = " << elementByRow[k2]
                                      << ", pseudoRedCost = " << pseudoReducedCost
                                      << std::endl;
#endif
                            if (pseudoReducedCost > maxPseudoReducedCost)
                                maxPseudoReducedCost = pseudoReducedCost;
                        } else {
                            // save value
                            chosenValue = fabs(elementByRow[k]);
                        }
                    }
                    assert (chosenValue);
                    maxPseudoReducedCost /= chosenValue;
#ifdef DIVE_DEBUG
                    std::cout << ", maxPseudoRedCost = " << maxPseudoReducedCost << std::endl;
#endif
                    candidate[cnt].var = iColumn1;
                    candidate[cnt++].pseudoRedCost = maxPseudoReducedCost;
                }
            }
#ifdef DIVE_DEBUG
            std::cout << "candidates for rounding = " << cnt << std::endl;
#endif
            std::sort(candidate, candidate + cnt, compareBinaryVars);
            for (int i = 0; i < cnt; i++) {
                int iColumn = candidate[i].var;
                if (numberAtBoundFixed < maxNumberAtBoundToFix) {
                    columnFixed[numberAtBoundFixed] = iColumn;
                    originalBound[numberAtBoundFixed] = upper[iColumn];
                    fixedAtLowerBound[numberAtBoundFixed] = true;
                    solver->setColUpper(iColumn, lower[iColumn]);
                    numberAtBoundFixed++;
                    if (numberAtBoundFixed == maxNumberAtBoundToFix)
                        break;
                }
            }
        }
#endif

        // fix other integer variables that are at their bounds
        int cnt = 0;
#ifdef GAP
        double gap = 1.0e30;
#endif
        if (reducedCost && true) {
#ifndef JJF_ONE
            cnt = fixOtherVariables(solver, solution, candidate, random);
#else
#ifdef GAP
            double cutoff = model_->getCutoff() ;
            if (cutoff < 1.0e20 && false) {
                double direction = solver->getObjSense() ;
                gap = cutoff - solver->getObjValue() * direction ;
                gap *= 0.1; // Fix more if plausible
                double tolerance;
                solver->getDblParam(OsiDualTolerance, tolerance) ;
                if (gap <= 0.0)
                    gap = tolerance;
                gap += 100.0 * tolerance;
            }
            int nOverGap = 0;
#endif
            int numberFree = 0;
            int numberFixed = 0;
            for (int i = 0; i < numberIntegers; i++) {
                int iColumn = integerVariable[i];
                if (upper[iColumn] > lower[iColumn]) {
                    numberFree++;
                    double value = newSolution[iColumn];
                    if (fabs(floor(value + 0.5) - value) <= integerTolerance) {
                        candidate[cnt].var = iColumn;
                        candidate[cnt++].pseudoRedCost =
                            fabs(reducedCost[iColumn] * random[i]);
#ifdef GAP
                        if (fabs(reducedCost[iColumn]) > gap)
                            nOverGap++;
#endif
                    }
                } else {
                    numberFixed++;
                }
            }
#ifdef GAP
            int nLeft = maxNumberAtBoundToFix - numberAtBoundFixed;
#ifdef CLP_INVESTIGATE4
            printf("cutoff %g obj %g nover %d - %d free, %d fixed\n",
                   cutoff, solver->getObjValue(), nOverGap, numberFree, numberFixed);
#endif
            if (nOverGap > nLeft && true) {
                nOverGap = CoinMin(nOverGap, nLeft + maxNumberAtBoundToFix / 2);
                maxNumberAtBoundToFix += nOverGap - nLeft;
            }
#else
#ifdef CLP_INVESTIGATE4
            printf("cutoff %g obj %g - %d free, %d fixed\n",
                   model_->getCutoff(), solver->getObjValue(), numberFree, numberFixed);
#endif
#endif
#endif
        } else {
            for (int i = 0; i < numberIntegers; i++) {
                int iColumn = integerVariable[i];
                if (upper[iColumn] > lower[iColumn]) {
                    double value = newSolution[iColumn];
                    if (fabs(floor(value + 0.5) - value) <= integerTolerance) {
                        candidate[cnt].var = iColumn;
                        candidate[cnt++].pseudoRedCost = numberIntegers - i;
                    }
                }
            }
        }
        std::sort(candidate, candidate + cnt, compareBinaryVars);
        for (int i = 0; i < cnt; i++) {
            int iColumn = candidate[i].var;
            if (upper[iColumn] > lower[iColumn]) {
                double value = newSolution[iColumn];
                if (fabs(floor(value + 0.5) - value) <= integerTolerance &&
                        numberAtBoundFixed < maxNumberAtBoundToFix) {
                    // fix the variable at one of its bounds
                    if (fabs(lower[iColumn] - value) <= integerTolerance) {
                        columnFixed[numberAtBoundFixed] = iColumn;
                        originalBound[numberAtBoundFixed] = upper[iColumn];
                        fixedAtLowerBound[numberAtBoundFixed] = true;
                        solver->setColUpper(iColumn, lower[iColumn]);
                        numberAtBoundFixed++;
                    } else if (fabs(upper[iColumn] - value) <= integerTolerance) {
                        columnFixed[numberAtBoundFixed] = iColumn;
                        originalBound[numberAtBoundFixed] = lower[iColumn];
                        fixedAtLowerBound[numberAtBoundFixed] = false;
                        solver->setColLower(iColumn, upper[iColumn]);
                        numberAtBoundFixed++;
                    }
                    if (numberAtBoundFixed == maxNumberAtBoundToFix)
                        break;
                }
            }
        }
#ifdef DIVE_DEBUG
        std::cout << "numberAtBoundFixed = " << numberAtBoundFixed << std::endl;
#endif

        double originalBoundBestColumn;
        if (bestColumn >= 0) {
            if (bestRound < 0) {
                originalBoundBestColumn = upper[bestColumn];
                solver->setColUpper(bestColumn, floor(newSolution[bestColumn]));
            } else {
                originalBoundBestColumn = lower[bestColumn];
                solver->setColLower(bestColumn, ceil(newSolution[bestColumn]));
            }
        } else {
            break;
        }
        int originalBestRound = bestRound;
        int saveModelOptions = model_->specialOptions();
        while (1) {

            model_->setSpecialOptions(saveModelOptions | 2048);
            solver->resolve();
            model_->setSpecialOptions(saveModelOptions);
            if (!solver->isAbandoned()) {
                numberSimplexIterations += solver->getIterationCount();
            } else {
                numberSimplexIterations = maxSimplexIterations + 1;
                break;
            }

            if (!solver->isProvenOptimal()) {
                if (numberAtBoundFixed > 0) {
                    // Remove the bound fix for variables that were at bounds
                    for (int i = 0; i < numberAtBoundFixed; i++) {
                        int iColFixed = columnFixed[i];
                        if (fixedAtLowerBound[i])
                            solver->setColUpper(iColFixed, originalBound[i]);
                        else
                            solver->setColLower(iColFixed, originalBound[i]);
                    }
                    numberAtBoundFixed = 0;
                } else if (bestRound == originalBestRound) {
                    bestRound *= (-1);
                    if (bestRound < 0) {
                        solver->setColLower(bestColumn, originalBoundBestColumn);
                        solver->setColUpper(bestColumn, floor(newSolution[bestColumn]));
                    } else {
                        solver->setColLower(bestColumn, ceil(newSolution[bestColumn]));
                        solver->setColUpper(bestColumn, originalBoundBestColumn);
                    }
                } else
                    break;
            } else
                break;
        }

        if (!solver->isProvenOptimal() ||
                direction*solver->getObjValue() >= solutionValue) {
#ifdef DIVE_DEBUG
            reasonToStop = 1;
#endif
            break;
        }

        if (iteration > maxIterations_) {
#ifdef DIVE_DEBUG
            reasonToStop = 2;
#endif
            break;
        }

        if (CoinCpuTime() - time1 > maxTime_) {
#ifdef DIVE_DEBUG
            reasonToStop = 3;
#endif
            break;
        }

        if (numberSimplexIterations > maxSimplexIterations) {
#ifdef DIVE_DEBUG
            reasonToStop = 4;
#endif
            // also switch off
#ifdef CLP_INVESTIGATE
            printf("switching off diving as too many iterations %d, %d allowed\n",
                   numberSimplexIterations, maxSimplexIterations);
#endif
            when_ = 0;
            break;
        }

        if (solver->getIterationCount() > 1000 && iteration > 3) {
#ifdef DIVE_DEBUG
            reasonToStop = 5;
#endif
            // also switch off
#ifdef CLP_INVESTIGATE
            printf("switching off diving one iteration took %d iterations (total %d)\n",
                   solver->getIterationCount(), numberSimplexIterations);
#endif
            when_ = 0;
            break;
        }

        memcpy(newSolution, solution, numberColumns*sizeof(double));
        numberFractionalVariables = 0;
        for (int i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            double value = newSolution[iColumn];
            if (fabs(floor(value + 0.5) - value) > integerTolerance) {
                numberFractionalVariables++;
            }
        }

    }


    double * rowActivity = new double[numberRows];
    memset(rowActivity, 0, numberRows*sizeof(double));

    // re-compute new solution value
    double objOffset = 0.0;
    solver->getDblParam(OsiObjOffset, objOffset);
    newSolutionValue = -objOffset;
    for (int i = 0 ; i < numberColumns ; i++ )
        newSolutionValue += objective[i] * newSolution[i];
    newSolutionValue *= direction;
    //printf("new solution value %g %g\n",newSolutionValue,solutionValue);
    if (newSolutionValue < solutionValue) {
        // paranoid check
        memset(rowActivity, 0, numberRows*sizeof(double));
        for (int i = 0; i < numberColumns; i++) {
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
        // check was approximately feasible
        bool feasible = true;
        for (int i = 0; i < numberRows; i++) {
            if (rowActivity[i] < rowLower[i]) {
                if (rowActivity[i] < rowLower[i] - 1000.0*primalTolerance)
                    feasible = false;
            } else if (rowActivity[i] > rowUpper[i]) {
                if (rowActivity[i] > rowUpper[i] + 1000.0*primalTolerance)
                    feasible = false;
            }
        }
        for (int i = 0; i < numberIntegers; i++) {
            int iColumn = integerVariable[i];
            double value = newSolution[iColumn];
            if (fabs(floor(value + 0.5) - value) > integerTolerance) {
                feasible = false;
                break;
            }
        }
        if (feasible) {
            // new solution
            memcpy(betterSolution, newSolution, numberColumns*sizeof(double));
            solutionValue = newSolutionValue;
            //printf("** Solution of %g found by CbcHeuristicDive\n",newSolutionValue);
            returnCode = 1;
        } else {
            // Can easily happen
            //printf("Debug CbcHeuristicDive giving bad solution\n");
        }
    }

#ifdef DIVE_DEBUG
    std::cout << "nRoundInfeasible = " << nRoundInfeasible
              << ", nRoundFeasible = " << nRoundFeasible
              << ", returnCode = " << returnCode
              << ", reasonToStop = " << reasonToStop
              << ", simplexIts = " << numberSimplexIterations
              << ", iterations = " << iteration << std::endl;
#endif

    delete [] newSolution;
    delete [] columnFixed;
    delete [] originalBound;
    delete [] fixedAtLowerBound;
    delete [] candidate;
    delete [] rowActivity;
    delete [] random;
    delete [] downArray_;
    downArray_ = NULL;
    delete [] upArray_;
    upArray_ = NULL;
    delete solver;
    return returnCode;
}

// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void
CbcHeuristicDive::validate()
{
    if (model_ && (when() % 100) < 10) {
        if (model_->numberIntegers() !=
                model_->numberObjects() && (model_->numberObjects() ||
                                            (model_->specialOptions()&1024) == 0)) {
            int numberOdd = 0;
            for (int i = 0; i < model_->numberObjects(); i++) {
                if (!model_->object(i)->canDoHeuristics())
                    numberOdd++;
            }
            if (numberOdd)
                setWhen(0);
        }
    }

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();
    delete [] downLocks_;
    delete [] upLocks_;
    downLocks_ = new unsigned short [numberIntegers];
    upLocks_ = new unsigned short [numberIntegers];
    // Column copy
    const double * element = matrix_.getElements();
    const int * row = matrix_.getIndices();
    const CoinBigIndex * columnStart = matrix_.getVectorStarts();
    const int * columnLength = matrix_.getVectorLengths();
    const double * rowLower = model_->solver()->getRowLower();
    const double * rowUpper = model_->solver()->getRowUpper();
    for (int i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        int down = 0;
        int up = 0;
        if (columnLength[iColumn] > 65535) {
            setWhen(0);
            break; // unlikely to work
        }
        for (CoinBigIndex j = columnStart[iColumn];
                j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            if (rowLower[iRow] > -1.0e20 && rowUpper[iRow] < 1.0e20) {
                up++;
                down++;
            } else if (element[j] > 0.0) {
                if (rowUpper[iRow] < 1.0e20)
                    up++;
                else
                    down++;
            } else {
                if (rowLower[iRow] > -1.0e20)
                    up++;
                else
                    down++;
            }
        }
        downLocks_[i] = static_cast<unsigned short> (down);
        upLocks_[i] = static_cast<unsigned short> (up);
    }

#ifdef DIVE_FIX_BINARY_VARIABLES
    selectBinaryVariables();
#endif
}

// Select candidate binary variables for fixing
void
CbcHeuristicDive::selectBinaryVariables()
{
    // Row copy
    const double * elementByRow = matrixByRow_.getElements();
    const int * column = matrixByRow_.getIndices();
    const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
    const int * rowLength = matrixByRow_.getVectorLengths();

    const int numberRows = matrixByRow_.getNumRows();
    const int numberCols = matrixByRow_.getNumCols();

    const double * lower = model_->solver()->getColLower();
    const double * upper = model_->solver()->getColUpper();
    const double * rowLower = model_->solver()->getRowLower();
    const double * rowUpper = model_->solver()->getRowUpper();

    //  const char * integerType = model_->integerType();


    //  const int numberIntegers = model_->numberIntegers();
    //  const int * integerVariable = model_->integerVariable();
    const double * objective = model_->solver()->getObjCoefficients();

    // vector to store the row number of variable bound rows
    int* rowIndexes = new int [numberCols];
    memset(rowIndexes, -1, numberCols*sizeof(int));

    for (int i = 0; i < numberRows; i++) {
        int positiveBinary = -1;
        int negativeBinary = -1;
        int nPositiveOther = 0;
        int nNegativeOther = 0;
        for (int k = rowStart[i]; k < rowStart[i] + rowLength[i]; k++) {
            int iColumn = column[k];
            if (model_->solver()->isInteger(iColumn) &&
                    lower[iColumn] == 0.0 && upper[iColumn] == 1.0 &&
                    objective[iColumn] == 0.0 &&
                    elementByRow[k] > 0.0 &&
                    positiveBinary < 0)
                positiveBinary = iColumn;
            else if (model_->solver()->isInteger(iColumn) &&
                     lower[iColumn] == 0.0 && upper[iColumn] == 1.0 &&
                     objective[iColumn] == 0.0 &&
                     elementByRow[k] < 0.0 &&
                     negativeBinary < 0)
                negativeBinary = iColumn;
            else if ((elementByRow[k] > 0.0 &&
                      lower[iColumn] >= 0.0) ||
                     (elementByRow[k] < 0.0 &&
                      upper[iColumn] <= 0.0))
                nPositiveOther++;
            else if ((elementByRow[k] > 0.0 &&
                      lower[iColumn] <= 0.0) ||
                     (elementByRow[k] < 0.0 &&
                      upper[iColumn] >= 0.0))
                nNegativeOther++;
            if (nPositiveOther > 0 && nNegativeOther > 0)
                break;
        }
        int binVar = -1;
        if (positiveBinary >= 0 &&
                (negativeBinary >= 0 || nNegativeOther > 0) &&
                nPositiveOther == 0 &&
                rowLower[i] == 0.0 &&
                rowUpper[i] > 0.0)
            binVar = positiveBinary;
        else if (negativeBinary >= 0 &&
                 (positiveBinary >= 0 || nPositiveOther > 0) &&
                 nNegativeOther == 0 &&
                 rowLower[i] < 0.0 &&
                 rowUpper[i] == 0.0)
            binVar = negativeBinary;
        if (binVar >= 0) {
            if (rowIndexes[binVar] == -1)
                rowIndexes[binVar] = i;
            else if (rowIndexes[binVar] >= 0)
                rowIndexes[binVar] = -2;
        }
    }

    for (int j = 0; j < numberCols; j++) {
        if (rowIndexes[j] >= 0) {
            binVarIndex_.push_back(j);
            vbRowIndex_.push_back(rowIndexes[j]);
        }
    }

#ifdef DIVE_DEBUG
    std::cout << "number vub Binary = " << binVarIndex_.size() << std::endl;
#endif

    delete [] rowIndexes;

}

/*
  Perform reduced cost fixing on integer variables.

  The variables in question are already nonbasic at bound. We're just nailing
  down the current situation.
*/

int CbcHeuristicDive::reducedCostFix (OsiSolverInterface* solver)

{
    //return 0; // temp
#ifndef JJF_ONE
    if (!model_->solverCharacteristics()->reducedCostsAccurate())
        return 0; //NLP
#endif
    double cutoff = model_->getCutoff() ;
    if (cutoff > 1.0e20)
        return 0;
#ifdef DIVE_DEBUG
    std::cout << "cutoff = " << cutoff << std::endl;
#endif
    double direction = solver->getObjSense() ;
    double gap = cutoff - solver->getObjValue() * direction ;
    gap *= 0.5; // Fix more
    double tolerance;
    solver->getDblParam(OsiDualTolerance, tolerance) ;
    if (gap <= 0.0)
        gap = tolerance; //return 0;
    gap += 100.0 * tolerance;
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);

    const double *lower = solver->getColLower() ;
    const double *upper = solver->getColUpper() ;
    const double *solution = solver->getColSolution() ;
    const double *reducedCost = solver->getReducedCost() ;

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();

    int numberFixed = 0 ;

# ifdef COIN_HAS_CLP
    OsiClpSolverInterface * clpSolver
    = dynamic_cast<OsiClpSolverInterface *> (solver);
    ClpSimplex * clpSimplex = NULL;
    if (clpSolver)
        clpSimplex = clpSolver->getModelPtr();
# endif
    for (int i = 0 ; i < numberIntegers ; i++) {
        int iColumn = integerVariable[i] ;
        double djValue = direction * reducedCost[iColumn] ;
        if (upper[iColumn] - lower[iColumn] > integerTolerance) {
            if (solution[iColumn] < lower[iColumn] + integerTolerance && djValue > gap) {
#ifdef COIN_HAS_CLP
                // may just have been fixed before
                if (clpSimplex) {
                    if (clpSimplex->getColumnStatus(iColumn) == ClpSimplex::basic) {
#ifdef COIN_DEVELOP
                        printf("DJfix %d has status of %d, dj of %g gap %g, bounds %g %g\n",
                               iColumn, clpSimplex->getColumnStatus(iColumn),
                               djValue, gap, lower[iColumn], upper[iColumn]);
#endif
                    } else {
                        assert(clpSimplex->getColumnStatus(iColumn) == ClpSimplex::atLowerBound ||
                               clpSimplex->getColumnStatus(iColumn) == ClpSimplex::isFixed);
                    }
                }
#endif
                solver->setColUpper(iColumn, lower[iColumn]) ;
                numberFixed++ ;
            } else if (solution[iColumn] > upper[iColumn] - integerTolerance && -djValue > gap) {
#ifdef COIN_HAS_CLP
                // may just have been fixed before
                if (clpSimplex) {
                    if (clpSimplex->getColumnStatus(iColumn) == ClpSimplex::basic) {
#ifdef COIN_DEVELOP
                        printf("DJfix %d has status of %d, dj of %g gap %g, bounds %g %g\n",
                               iColumn, clpSimplex->getColumnStatus(iColumn),
                               djValue, gap, lower[iColumn], upper[iColumn]);
#endif
                    } else {
                        assert(clpSimplex->getColumnStatus(iColumn) == ClpSimplex::atUpperBound ||
                               clpSimplex->getColumnStatus(iColumn) == ClpSimplex::isFixed);
                    }
                }
#endif
                solver->setColLower(iColumn, upper[iColumn]) ;
                numberFixed++ ;
            }
        }
    }
    return numberFixed;
}
// Fix other variables at bounds
int
CbcHeuristicDive::fixOtherVariables(OsiSolverInterface * solver,
                                    const double * solution,
                                    PseudoReducedCost * candidate,
                                    const double * random)
{
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double primalTolerance;
    solver->getDblParam(OsiPrimalTolerance, primalTolerance);

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();
    const double* reducedCost = solver->getReducedCost();
    // fix other integer variables that are at their bounds
    int cnt = 0;
#ifdef GAP
    double direction = solver->getObjSense(); // 1 for min, -1 for max
    double gap = 1.0e30;
#endif
#ifdef GAP
    double cutoff = model_->getCutoff() ;
    if (cutoff < 1.0e20 && false) {
        double direction = solver->getObjSense() ;
        gap = cutoff - solver->getObjValue() * direction ;
        gap *= 0.1; // Fix more if plausible
        double tolerance;
        solver->getDblParam(OsiDualTolerance, tolerance) ;
        if (gap <= 0.0)
            gap = tolerance;
        gap += 100.0 * tolerance;
    }
    int nOverGap = 0;
#endif
    int numberFree = 0;
    int numberFixedAlready = 0;
    for (int i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        if (upper[iColumn] > lower[iColumn]) {
            numberFree++;
            double value = solution[iColumn];
            if (fabs(floor(value + 0.5) - value) <= integerTolerance) {
                candidate[cnt].var = iColumn;
                candidate[cnt++].pseudoRedCost =
                    fabs(reducedCost[iColumn] * random[i]);
#ifdef GAP
                if (fabs(reducedCost[iColumn]) > gap)
                    nOverGap++;
#endif
            }
        } else {
            numberFixedAlready++;
        }
    }
#ifdef GAP
    int nLeft = maxNumberToFix - numberFixedAlready;
#ifdef CLP_INVESTIGATE4
    printf("cutoff %g obj %g nover %d - %d free, %d fixed\n",
           cutoff, solver->getObjValue(), nOverGap, numberFree,
           numberFixedAlready);
#endif
    if (nOverGap > nLeft && true) {
        nOverGap = CoinMin(nOverGap, nLeft + maxNumberToFix / 2);
        maxNumberToFix += nOverGap - nLeft;
    }
#else
#ifdef CLP_INVESTIGATE4
    printf("cutoff %g obj %g - %d free, %d fixed\n",
           model_->getCutoff(), solver->getObjValue(), numberFree,
           numberFixedAlready);
#endif
#endif
    return cnt;
}

