// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CoinPragma.hpp"
#include "CoinTypes.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcSubProblem.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif

// Default Constructor
CbcSubProblem::CbcSubProblem()
        : objectiveValue_(0.0),
        sumInfeasibilities_(0.0),
        variables_(NULL),
        newBounds_(NULL),
        status_(NULL),
        depth_(0),
        numberChangedBounds_(0),
        numberInfeasibilities_(0)
{
}

// Useful constructor
CbcSubProblem::CbcSubProblem (const OsiSolverInterface * solver,
                              const double * lastLower,
                              const double * lastUpper,
                              const unsigned char * status,
                              int depth)
        : objectiveValue_(0.0),
        sumInfeasibilities_(0.0),
        variables_(NULL),
        newBounds_(NULL),
        status_(NULL),
        depth_(depth),
        numberChangedBounds_(0),
        numberInfeasibilities_(0)
{
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();

    numberChangedBounds_ = 0;
    int numberColumns = solver->getNumCols();
    int i;
    for (i = 0; i < numberColumns; i++) {
        if (lower[i] != lastLower[i])
            numberChangedBounds_++;
        if (upper[i] != lastUpper[i])
            numberChangedBounds_++;
    }
    if (numberChangedBounds_) {
        newBounds_ = new double [numberChangedBounds_] ;
        variables_ = new int [numberChangedBounds_] ;
        numberChangedBounds_ = 0;
        for (i = 0; i < numberColumns; i++) {
            if (lower[i] != lastLower[i]) {
                variables_[numberChangedBounds_] = i;
                newBounds_[numberChangedBounds_++] = lower[i];
            }
            if (upper[i] != lastUpper[i]) {
                variables_[numberChangedBounds_] = i | 0x80000000;
                newBounds_[numberChangedBounds_++] = upper[i];
            }
#ifdef CBC_DEBUG
            if (lower[i] != lastLower[i]) {
                std::cout
                    << "lower on " << i << " changed from "
                    << lastLower[i] << " to " << lower[i] << std::endl ;
            }
            if (upper[i] != lastUpper[i]) {
                std::cout
                    << "upper on " << i << " changed from "
                    << lastUpper[i] << " to " << upper[i] << std::endl ;
            }
#endif
        }
#ifdef CBC_DEBUG
        std::cout << numberChangedBounds_ << " changed bounds." << std::endl ;
#endif
    }
    const OsiClpSolverInterface * clpSolver
    = dynamic_cast<const OsiClpSolverInterface *> (solver);
    assert (clpSolver);
    // Do difference
    // Current basis
    status_ = clpSolver->getBasis(status);
    assert (status_->fullBasis());
    //status_->print();
}

// Copy constructor
CbcSubProblem::CbcSubProblem ( const CbcSubProblem & rhs)
        : objectiveValue_(rhs.objectiveValue_),
        sumInfeasibilities_(rhs.sumInfeasibilities_),
        variables_(NULL),
        newBounds_(NULL),
        status_(NULL),
        depth_(rhs.depth_),
        numberChangedBounds_(rhs.numberChangedBounds_),
        numberInfeasibilities_(rhs.numberInfeasibilities_)
{
    if (numberChangedBounds_) {
        variables_ = CoinCopyOfArray(rhs.variables_, numberChangedBounds_);
        newBounds_ = CoinCopyOfArray(rhs.newBounds_, numberChangedBounds_);
    }
    if (rhs.status_) {
        status_ = new CoinWarmStartBasis(*rhs.status_);
    }
}

// Assignment operator
CbcSubProblem &
CbcSubProblem::operator=( const CbcSubProblem & rhs)
{
    if (this != &rhs) {
        delete [] variables_;
        delete [] newBounds_;
        delete status_;
        objectiveValue_ = rhs.objectiveValue_;
        sumInfeasibilities_ = rhs.sumInfeasibilities_;
        depth_ = rhs.depth_;
        numberChangedBounds_ = rhs.numberChangedBounds_;
        numberInfeasibilities_ = rhs.numberInfeasibilities_;
        if (numberChangedBounds_) {
            variables_ = CoinCopyOfArray(rhs.variables_, numberChangedBounds_);
            newBounds_ = CoinCopyOfArray(rhs.newBounds_, numberChangedBounds_);
        } else {
            variables_ = NULL;
            newBounds_ = NULL;
        }
        if (rhs.status_) {
            status_ = new CoinWarmStartBasis(*rhs.status_);
        } else {
            status_ = NULL;
        }
    }
    return *this;
}

// Destructor
CbcSubProblem::~CbcSubProblem ()
{
    delete [] variables_;
    delete [] newBounds_;
    delete status_;
}
// Apply subproblem
void
CbcSubProblem::apply(OsiSolverInterface * solver, int what) const
{
    int i;
    if ((what&1) != 0) {
#ifndef NDEBUG
        int nSame = 0;
#endif
        for (i = 0; i < numberChangedBounds_; i++) {
            int variable = variables_[i];
            int k = variable & 0x3fffffff;
            if ((variable&0x80000000) == 0) {
                // lower bound changing
                //#define CBC_PRINT2
#ifdef CBC_PRINT2
                if (solver->getColLower()[k] != newBounds_[i])
                    printf("lower change for column %d - from %g to %g\n", k, solver->getColLower()[k], newBounds_[i]);
#endif
#ifndef NDEBUG
                if ((variable&0x40000000) == 0 && true) {
                    double oldValue = solver->getColLower()[k];
                    assert (newBounds_[i] > oldValue - 1.0e-8);
                    if (newBounds_[i] < oldValue + 1.0e-8) {
#ifdef CBC_PRINT2
                        printf("bad null lower change for column %d - bound %g\n", k, oldValue);
#endif
                        if (newBounds_[i] == oldValue)
                            nSame++;
                    }
                }
#endif
                solver->setColLower(k, newBounds_[i]);
            } else {
                // upper bound changing
#ifdef CBC_PRINT2
                if (solver->getColUpper()[k] != newBounds_[i])
                    printf("upper change for column %d - from %g to %g\n", k, solver->getColUpper()[k], newBounds_[i]);
#endif
#ifndef NDEBUG
                if ((variable&0x40000000) == 0 && true) {
                    double oldValue = solver->getColUpper()[k];
                    assert (newBounds_[i] < oldValue + 1.0e-8);
                    if (newBounds_[i] > oldValue - 1.0e-8) {
#ifdef CBC_PRINT2
                        printf("bad null upper change for column %d - bound %g\n", k, oldValue);
#endif
                        if (newBounds_[i] == oldValue)
                            nSame++;
                    }
                }
#endif
                solver->setColUpper(k, newBounds_[i]);
            }
        }
#ifndef NDEBUG
#ifdef CBC_PRINT2
        if (nSame && (nSame < numberChangedBounds_ || (what&3) != 3))
            printf("%d changes out of %d redundant %d\n",
                   nSame, numberChangedBounds_, what);
        else if (numberChangedBounds_ && what == 7 && !nSame)
            printf("%d good changes %d\n",
                   numberChangedBounds_, what);
#endif
#endif
    }
#ifdef JJF_ZERO
    if ((what&2) != 0) {
        OsiClpSolverInterface * clpSolver
        = dynamic_cast<OsiClpSolverInterface *> (solver);
        assert (clpSolver);
        //assert (clpSolver->getNumRows()==numberRows_);
        //clpSolver->setBasis(*status_);
        // Current basis
        CoinWarmStartBasis * basis = clpSolver->getPointerToWarmStart();
        printf("BBBB\n");
        basis->print();
        assert (basis->fullBasis());
        basis->applyDiff(status_);
        printf("diff applied %x\n", status_);
        printf("CCCC\n");
        basis->print();
        assert (basis->fullBasis());
#ifndef NDEBUG
        if (!basis->fullBasis())
            printf("Debug this basis!!\n");
#endif
        clpSolver->setBasis(*basis);
    }
#endif
    if ((what&8) != 0) {
        OsiClpSolverInterface * clpSolver
        = dynamic_cast<OsiClpSolverInterface *> (solver);
        assert (clpSolver);
        clpSolver->setBasis(*status_);
        delete status_;
        status_ = NULL;
    }
}

