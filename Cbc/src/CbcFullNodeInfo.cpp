// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/24/09 carved from CbcNode

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CbcConfig.h"

#include <string>
//#define CBC_DEBUG 1
//#define CHECK_CUT_COUNTS
//#define CHECK_NODE
//#define CBC_CHECK_BASIS
#include <cassert>
#include <cfloat>
#define CUTS
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

CbcFullNodeInfo::CbcFullNodeInfo() :
        CbcNodeInfo(),
        basis_(),
        numberIntegers_(0),
        lower_(NULL),
        upper_(NULL)
{
}
CbcFullNodeInfo::CbcFullNodeInfo(CbcModel * model,
                                 int numberRowsAtContinuous) :
        CbcNodeInfo(NULL, model->currentNode())
{
    OsiSolverInterface * solver = model->solver();
    numberRows_ = numberRowsAtContinuous;
    numberIntegers_ = model->numberIntegers();
    int numberColumns = model->getNumCols();
    lower_ = new double [numberColumns];
    upper_ = new double [numberColumns];
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    int i;

    for (i = 0; i < numberColumns; i++) {
        lower_[i] = lower[i];
        upper_[i] = upper[i];
    }

    basis_ =  dynamic_cast<CoinWarmStartBasis*>(solver->getWarmStart());
}

CbcFullNodeInfo::CbcFullNodeInfo(const CbcFullNodeInfo & rhs) :
        CbcNodeInfo(rhs)
{
    basis_ = dynamic_cast<CoinWarmStartBasis *>(rhs.basis_->clone()) ;
    numberIntegers_ = rhs.numberIntegers_;
    lower_ = NULL;
    upper_ = NULL;
    if (rhs.lower_ != NULL) {
        int numberColumns = basis_->getNumStructural();
        lower_ = new double [numberColumns];
        upper_ = new double [numberColumns];
        assert (upper_ != NULL);
        memcpy(lower_, rhs.lower_, numberColumns*sizeof(double));
        memcpy(upper_, rhs.upper_, numberColumns*sizeof(double));
    }
}

CbcNodeInfo *
CbcFullNodeInfo::clone() const
{
    return (new CbcFullNodeInfo(*this));
}

CbcFullNodeInfo::~CbcFullNodeInfo ()
{
    delete basis_ ;
    delete [] lower_;
    delete [] upper_;
}

/*
  The basis supplied as a parameter is deleted and replaced with a new basis
  appropriate for the node, and lower and upper bounds on variables are
  reset according to the stored bounds arrays. Any cuts associated with this
  node are added to the list in addCuts, but not actually added to the
  constraint system in the model.

  Why pass in a basis at all? The short answer is ``We need the parameter to
  pass out a basis, so might as well use it to pass in the size.''

  A longer answer is that in practice we take a memory allocation hit up in
  addCuts1 (the only place applyToModel is called) when we setSize() the
  basis that's passed in. It's immediately tossed here in favour of a clone
  of the basis attached to this nodeInfo. This can probably be fixed, given
  a bit of thought.
*/

void CbcFullNodeInfo::applyToModel (CbcModel *model,
                                    CoinWarmStartBasis *&basis,
                                    CbcCountRowCut **addCuts,
                                    int &currentNumberCuts) const

{
    OsiSolverInterface *solver = model->solver() ;

    // branch - do bounds
    assert (active_ == 7 || active_ == 15);
    int i;
    solver->setColLower(lower_);
    solver->setColUpper(upper_);
    int numberColumns = model->getNumCols();
    // move basis - but make sure size stays
    // for bon-min - should not be needed int numberRows = model->getNumRows();
    int numberRows = basis->getNumArtificial();
    delete basis ;
    if (basis_) {
        basis = dynamic_cast<CoinWarmStartBasis *>(basis_->clone()) ;
        basis->resize(numberRows, numberColumns);
#ifdef CBC_CHECK_BASIS
        std::cout << "Basis (after applying root " << this << ") " << std::endl ;
        basis->print() ;
#endif
    } else {
        // We have a solver without a basis
        basis = NULL;
    }
    for (i = 0; i < numberCuts_; i++)
        addCuts[currentNumberCuts+i] = cuts_[i];
    currentNumberCuts += numberCuts_;
    assert(!parent_);
    return ;
}
// Just apply bounds to one variable (1=>infeasible)
int
CbcFullNodeInfo::applyBounds(int iColumn, double & lower, double & upper, int force)
{
    if ((force && 1) == 0) {
      if (lower > lower_[iColumn])
	COIN_DETAIL_PRINT(printf("%d odd lower going from %g to %g\n", iColumn, lower, lower_[iColumn]));
        lower = lower_[iColumn];
    } else {
        lower_[iColumn] = lower;
    }
    if ((force && 2) == 0) {
      if (upper < upper_[iColumn])
	COIN_DETAIL_PRINT(printf("%d odd upper going from %g to %g\n", iColumn, upper, upper_[iColumn]));
        upper = upper_[iColumn];
    } else {
        upper_[iColumn] = upper;
    }
    return (upper_[iColumn] >= lower_[iColumn]) ? 0 : 1;
}

/* Builds up row basis backwards (until original model).
   Returns NULL or previous one to apply .
   Depends on Free being 0 and impossible for cuts
*/
CbcNodeInfo *
CbcFullNodeInfo::buildRowBasis(CoinWarmStartBasis & basis ) const
{
    const unsigned int * saved =
        reinterpret_cast<const unsigned int *> (basis_->getArtificialStatus());
    unsigned int * now =
        reinterpret_cast<unsigned int *> (basis.getArtificialStatus());
    int number = basis_->getNumArtificial() >> 4;;
    int i;
    for (i = 0; i < number; i++) {
        if (!now[i])
            now[i] = saved[i];
    }
    return NULL;
}

