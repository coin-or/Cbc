// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CoinTypes.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcFixVariable.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"


//##############################################################################

// Default Constructor
CbcFixVariable::CbcFixVariable ()
        : CbcConsequence(),
        numberStates_(0),
        states_(NULL),
        startLower_(NULL),
        startUpper_(NULL),
        newBound_(NULL),
        variable_(NULL)
{
}

// One useful Constructor
CbcFixVariable::CbcFixVariable (int numberStates, const int * states, const int * numberNewLower,
                                const int ** newLowerValue,
                                const int ** lowerColumn,
                                const int * numberNewUpper, const int ** newUpperValue,
                                const int ** upperColumn)
        : CbcConsequence(),
        states_(NULL),
        startLower_(NULL),
        startUpper_(NULL),
        newBound_(NULL),
        variable_(NULL)
{
    // How much space
    numberStates_ = numberStates;
    if (numberStates_) {
        states_ = new int[numberStates_];
        memcpy(states_, states, numberStates_*sizeof(int));
        int i;
        int n = 0;
        startLower_ = new int[numberStates_+1];
        startUpper_ = new int[numberStates_+1];
        startLower_[0] = 0;
        //count
        for (i = 0; i < numberStates_; i++) {
            n += numberNewLower[i];
            startUpper_[i] = n;
            n += numberNewUpper[i];
            startLower_[i+1] = n;
        }
        newBound_ = new double [n];
        variable_ = new int [n];
        n = 0;
        for (i = 0; i < numberStates_; i++) {
            int j;
            int k;
            const int * bound;
            const int * variable;
            k = numberNewLower[i];
            bound = newLowerValue[i];
            variable = lowerColumn[i];
            for (j = 0; j < k; j++) {
                newBound_[n] = bound[j];
                variable_[n++] = variable[j];
            }
            k = numberNewUpper[i];
            bound = newUpperValue[i];
            variable = upperColumn[i];
            for (j = 0; j < k; j++) {
                newBound_[n] = bound[j];
                variable_[n++] = variable[j];
            }
        }
    }
}

// Copy constructor
CbcFixVariable::CbcFixVariable ( const CbcFixVariable & rhs)
        : CbcConsequence(rhs)
{
    numberStates_ = rhs.numberStates_;
    states_ = NULL;
    startLower_ = NULL;
    startUpper_ = NULL;
    newBound_ = NULL;
    variable_ = NULL;
    if (numberStates_) {
        states_ = CoinCopyOfArray(rhs.states_, numberStates_);
        startLower_ = CoinCopyOfArray(rhs.startLower_, numberStates_ + 1);
        startUpper_ = CoinCopyOfArray(rhs.startUpper_, numberStates_ + 1);
        int n = startLower_[numberStates_];
        newBound_ = CoinCopyOfArray(rhs.newBound_, n);
        variable_ = CoinCopyOfArray(rhs.variable_, n);
    }
}

// Clone
CbcConsequence *
CbcFixVariable::clone() const
{
    return new CbcFixVariable(*this);
}

// Assignment operator
CbcFixVariable &
CbcFixVariable::operator=( const CbcFixVariable & rhs)
{
    if (this != &rhs) {
        CbcConsequence::operator=(rhs);
        delete [] states_;
        delete [] startLower_;
        delete [] startUpper_;
        delete [] newBound_;
        delete [] variable_;
        states_ = NULL;
        startLower_ = NULL;
        startUpper_ = NULL;
        newBound_ = NULL;
        variable_ = NULL;
        numberStates_ = rhs.numberStates_;
        if (numberStates_) {
            states_ = CoinCopyOfArray(rhs.states_, numberStates_);
            startLower_ = CoinCopyOfArray(rhs.startLower_, numberStates_ + 1);
            startUpper_ = CoinCopyOfArray(rhs.startUpper_, numberStates_ + 1);
            int n = startLower_[numberStates_];
            newBound_ = CoinCopyOfArray(rhs.newBound_, n);
            variable_ = CoinCopyOfArray(rhs.variable_, n);
        }
    }
    return *this;
}

// Destructor
CbcFixVariable::~CbcFixVariable ()
{
    delete [] states_;
    delete [] startLower_;
    delete [] startUpper_;
    delete [] newBound_;
    delete [] variable_;
}
// Set up a startLower for a single member
void
CbcFixVariable::applyToSolver(OsiSolverInterface * solver, int state) const
{
    assert (state == -9999 || state == 9999);
    // Find state
    int find;
    for (find = 0; find < numberStates_; find++)
        if (states_[find] == state)
            break;
    if (find == numberStates_)
        return;
    int i;
    // Set new lower bounds
    for (i = startLower_[find]; i < startUpper_[find]; i++) {
        int iColumn = variable_[i];
        double value = newBound_[i];
        double oldValue = solver->getColLower()[iColumn];
        //printf("for %d old lower bound %g, new %g",iColumn,oldValue,value);
        solver->setColLower(iColumn, CoinMax(value, oldValue));
        //printf(" => %g\n",solver->getColLower()[iColumn]);
    }
    // Set new upper bounds
    for (i = startUpper_[find]; i < startLower_[find+1]; i++) {
        int iColumn = variable_[i];
        double value = newBound_[i];
        double oldValue = solver->getColUpper()[iColumn];
        //printf("for %d old upper bound %g, new %g",iColumn,oldValue,value);
        solver->setColUpper(iColumn, CoinMin(value, oldValue));
        //printf(" => %g\n",solver->getColUpper()[iColumn]);
    }
}

