// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/24/09 carved from CbcNode

#include "CbcConfig.h"

#include <string>
//#define CBC_DEBUG 1
//#define CHECK_CUT_COUNTS
//#define CHECK_NODE
//#define CBC_CHECK_BASIS
#include <cassert>
#include <cfloat>
#define CUTS
#include "CoinPragma.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiChooseVariable.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiSolverBranch.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinTime.hpp"
#include "CbcModel.hpp"
#include "CbcNode.hpp"
#include "CbcStatistics.hpp"
#include "CbcStrategy.hpp"
#include "CbcBranchActual.hpp"
#include "CbcBranchDynamic.hpp"
#include "OsiRowCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiCuts.hpp"
#include "CbcCountRowCut.hpp"
#include "CbcFeasibilityBase.hpp"
#include "CbcMessage.hpp"
#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#include "ClpSimplexOther.hpp"
#endif
using namespace std;
#include "CglCutGenerator.hpp"

// Default constructor
CbcPartialNodeInfo::CbcPartialNodeInfo()

        : CbcNodeInfo(),
        basisDiff_(NULL),
        variables_(NULL),
        newBounds_(NULL),
        numberChangedBounds_(0)

{ /* this space intentionally left blank */ }

// Constructor from current state
CbcPartialNodeInfo::CbcPartialNodeInfo (CbcNodeInfo *parent, CbcNode *owner,
                                        int numberChangedBounds,
                                        const int *variables,
                                        const double *boundChanges,
                                        const CoinWarmStartDiff *basisDiff)
        : CbcNodeInfo(parent, owner)
{
    basisDiff_ = basisDiff->clone() ;
#ifdef CBC_CHECK_BASIS
    std::cout << "Constructor (" << this << ") " << std::endl ;
#endif

    numberChangedBounds_ = numberChangedBounds;
    size_t size = numberChangedBounds_ * (sizeof(double) + sizeof(int));
    char * temp = new char [size];
    newBounds_ = reinterpret_cast<double *> (temp);
    variables_ = reinterpret_cast<int *> (newBounds_ + numberChangedBounds_);

    int i ;
    for (i = 0; i < numberChangedBounds_; i++) {
        variables_[i] = variables[i];
        newBounds_[i] = boundChanges[i];
    }
}

CbcPartialNodeInfo::CbcPartialNodeInfo (const CbcPartialNodeInfo & rhs)

        : CbcNodeInfo(rhs)

{
    basisDiff_ = rhs.basisDiff_->clone() ;

#ifdef CBC_CHECK_BASIS
    std::cout << "Copy constructor (" << this << ") from " << this << std::endl ;
#endif
    numberChangedBounds_ = rhs.numberChangedBounds_;
    size_t size = numberChangedBounds_ * (sizeof(double) + sizeof(int));
    char * temp = new char [size];
    newBounds_ = reinterpret_cast<double *> (temp);
    variables_ = reinterpret_cast<int *> (newBounds_ + numberChangedBounds_);

    int i ;
    for (i = 0; i < numberChangedBounds_; i++) {
        variables_[i] = rhs.variables_[i];
        newBounds_[i] = rhs.newBounds_[i];
    }
}

CbcNodeInfo *
CbcPartialNodeInfo::clone() const
{
    return (new CbcPartialNodeInfo(*this));
}


CbcPartialNodeInfo::~CbcPartialNodeInfo ()
{
    delete basisDiff_ ;
    delete [] newBounds_;
}


/**
   The basis supplied as a parameter is incrementally modified, and lower and
   upper bounds on variables in the model are incrementally modified. Any
   cuts associated with this node are added to the list in addCuts.
*/

void CbcPartialNodeInfo::applyToModel (CbcModel *model,
                                       CoinWarmStartBasis *&basis,
                                       CbcCountRowCut **addCuts,
                                       int &currentNumberCuts) const

{
    OsiSolverInterface *solver = model->solver();
    if ((active_&4) != 0) {
        basis->applyDiff(basisDiff_) ;
#ifdef CBC_CHECK_BASIS
        std::cout << "Basis (after applying " << this << ") " << std::endl ;
        basis->print() ;
#endif
    }

    // branch - do bounds
    int i;
    if ((active_&1) != 0) {
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
                if ((variable&0x40000000) == 0 && false) {
                    double oldValue = solver->getColLower()[k];
                    assert (newBounds_[i] > oldValue - 1.0e-8);
                    if (newBounds_[i] < oldValue + 1.0e-8)
                        printf("bad null lower change for column %d - bound %g\n", k, oldValue);
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
                if ((variable&0x40000000) == 0 && false) {
                    double oldValue = solver->getColUpper()[k];
                    assert (newBounds_[i] < oldValue + 1.0e-8);
                    if (newBounds_[i] > oldValue - 1.0e-8)
                        printf("bad null upper change for column %d - bound %g\n", k, oldValue);
                }
#endif
                solver->setColUpper(k, newBounds_[i]);
            }
        }
    }
    if ((active_&2) != 0) {
        for (i = 0; i < numberCuts_; i++) {
            addCuts[currentNumberCuts+i] = cuts_[i];
            if (cuts_[i] && model->messageHandler()->logLevel() > 4) {
                cuts_[i]->print();
            }
        }

        currentNumberCuts += numberCuts_;
    }
    return ;
}
// Just apply bounds to one variable (1=>infeasible)
int
CbcPartialNodeInfo::applyBounds(int iColumn, double & lower, double & upper, int force)
{
    // branch - do bounds
    int i;
    int found = 0;
    double newLower = -COIN_DBL_MAX;
    double newUpper = COIN_DBL_MAX;
    for (i = 0; i < numberChangedBounds_; i++) {
        int variable = variables_[i];
        int k = variable & 0x3fffffff;
        if (k == iColumn) {
            if ((variable&0x80000000) == 0) {
                // lower bound changing
                found |= 1;
                newLower = CoinMax(newLower, newBounds_[i]);
                if ((force&1) == 0) {
                    if (lower > newBounds_[i])
		      COIN_DETAIL_PRINT(printf("%d odd lower going from %g to %g\n", iColumn, lower, newBounds_[i]));
                    lower = newBounds_[i];
                } else {
                    newBounds_[i] = lower;
                    variables_[i] |= 0x40000000; // say can go odd way
                }
            } else {
                // upper bound changing
                found |= 2;
                newUpper = CoinMin(newUpper, newBounds_[i]);
                if ((force&2) == 0) {
                    if (upper < newBounds_[i])
		      COIN_DETAIL_PRINT(printf("%d odd upper going from %g to %g\n", iColumn, upper, newBounds_[i]));
                    upper = newBounds_[i];
                } else {
                    newBounds_[i] = upper;
                    variables_[i] |= 0x40000000; // say can go odd way
                }
            }
        }
    }
    newLower = CoinMax(newLower, lower);
    newUpper = CoinMin(newUpper, upper);
    int nAdd = 0;
    if ((force&2) != 0 && (found&2) == 0) {
        // need to add new upper
        nAdd++;
    }
    if ((force&1) != 0 && (found&1) == 0) {
        // need to add new lower
        nAdd++;
    }
    if (nAdd) {
        size_t size = (numberChangedBounds_ + nAdd) * (sizeof(double) + sizeof(int));
        char * temp = new char [size];
        double * newBounds = reinterpret_cast<double *> (temp);
        int * variables = reinterpret_cast<int *> (newBounds + numberChangedBounds_ + nAdd);

        int i ;
        for (i = 0; i < numberChangedBounds_; i++) {
            variables[i] = variables_[i];
            newBounds[i] = newBounds_[i];
        }
        delete [] newBounds_;
        newBounds_ = newBounds;
        variables_ = variables;
        if ((force&2) != 0 && (found&2) == 0) {
            // need to add new upper
            int variable = iColumn | 0x80000000;
            variables_[numberChangedBounds_] = variable;
            newBounds_[numberChangedBounds_++] = newUpper;
        }
        if ((force&1) != 0 && (found&1) == 0) {
            // need to add new lower
            int variable = iColumn;
            variables_[numberChangedBounds_] = variable;
            newBounds_[numberChangedBounds_++] = newLower;
        }
    }

    return (newUpper >= newLower) ? 0 : 1;
}

/* Builds up row basis backwards (until original model).
   Returns NULL or previous one to apply .
   Depends on Free being 0 and impossible for cuts
*/

CbcNodeInfo *
CbcPartialNodeInfo::buildRowBasis(CoinWarmStartBasis & basis ) const

{
    basis.applyDiff(basisDiff_) ;

    return parent_ ;
}

