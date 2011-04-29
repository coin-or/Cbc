/* $Id: CbcHeuristicDiveGuided.cpp 1173 2009-06-04 09:44:10Z forrest $ */
// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CbcHeuristicDiveGuided.hpp"
#include "CbcStrategy.hpp"

// Default Constructor
CbcHeuristicDiveGuided::CbcHeuristicDiveGuided()
        : CbcHeuristicDive()
{
}

// Constructor from model
CbcHeuristicDiveGuided::CbcHeuristicDiveGuided(CbcModel & model)
        : CbcHeuristicDive(model)
{
}

// Destructor
CbcHeuristicDiveGuided::~CbcHeuristicDiveGuided ()
{
}

// Clone
CbcHeuristicDiveGuided *
CbcHeuristicDiveGuided::clone() const
{
    return new CbcHeuristicDiveGuided(*this);
}

// Create C++ lines to get to current state
void
CbcHeuristicDiveGuided::generateCpp( FILE * fp)
{
    CbcHeuristicDiveGuided other;
    fprintf(fp, "0#include \"CbcHeuristicDiveGuided.hpp\"\n");
    fprintf(fp, "3  CbcHeuristicDiveGuided heuristicDiveGuided(*cbcModel);\n");
    CbcHeuristic::generateCpp(fp, "heuristicDiveGuided");
    fprintf(fp, "3  cbcModel->addHeuristic(&heuristicDiveGuided);\n");
}

// Copy constructor
CbcHeuristicDiveGuided::CbcHeuristicDiveGuided(const CbcHeuristicDiveGuided & rhs)
        :
        CbcHeuristicDive(rhs)
{
}

// Assignment operator
CbcHeuristicDiveGuided &
CbcHeuristicDiveGuided::operator=( const CbcHeuristicDiveGuided & rhs)
{
    if (this != &rhs) {
        CbcHeuristicDive::operator=(rhs);
    }
    return *this;
}

bool
CbcHeuristicDiveGuided::canHeuristicRun()
{
    double* bestIntegerSolution = model_->bestSolution();
    if (bestIntegerSolution == NULL)
        return false; // no integer solution available. Switch off heuristic

    return CbcHeuristicDive::canHeuristicRun();
}

bool
CbcHeuristicDiveGuided::selectVariableToBranch(OsiSolverInterface* solver,
        const double* newSolution,
        int& bestColumn,
        int& bestRound)
{
    double* bestIntegerSolution = model_->bestSolution();

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);

    bestColumn = -1;
    bestRound = -1; // -1 rounds down, +1 rounds up
    double bestFraction = COIN_DBL_MAX;
    bool allTriviallyRoundableSoFar = true;
    for (int i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        double value = newSolution[iColumn];
        double fraction = value - floor(value);
        int round = 0;
        if (fabs(floor(value + 0.5) - value) > integerTolerance) {
            if (allTriviallyRoundableSoFar || (downLocks_[i] > 0 && upLocks_[i] > 0)) {

                if (allTriviallyRoundableSoFar && downLocks_[i] > 0 && upLocks_[i] > 0) {
                    allTriviallyRoundableSoFar = false;
                    bestFraction = COIN_DBL_MAX;
                }

                if (value >= bestIntegerSolution[iColumn])
                    round = -1;
                else {
                    round = 1;
                    fraction = 1.0 - fraction;
                }

                // if variable is not binary, penalize it
                if (!solver->isBinary(iColumn))
                    fraction *= 1000.0;

                if (fraction < bestFraction) {
                    bestColumn = iColumn;
                    bestFraction = fraction;
                    bestRound = round;
                }
            }
        }
    }
    return allTriviallyRoundableSoFar;
}

