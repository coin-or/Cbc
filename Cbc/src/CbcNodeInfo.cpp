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
#include "CbcNodeInfo.hpp"

// Default Constructor
CbcNodeInfo::CbcNodeInfo ()
        :
        numberPointingToThis_(0),
        parent_(NULL),
        parentBranch_(NULL),
        owner_(NULL),
        numberCuts_(0),
        nodeNumber_(0),
        cuts_(NULL),
        numberRows_(0),
        numberBranchesLeft_(0),
        active_(7)
{
#ifdef CHECK_NODE
    printf("CbcNodeInfo %x Constructor\n", this);
#endif
}

void
CbcNodeInfo::setParentBasedData()
{
    if (parent_) {
        numberRows_ = parent_->numberRows_ + parent_->numberCuts_;
        //parent_->increment();
        if (parent_->owner()) {
            const OsiBranchingObject* br = parent_->owner()->branchingObject();
            assert(br);
            parentBranch_ = br->clone();
        }
    }
}

void
CbcNodeInfo::unsetParentBasedData()
{
    if (parent_) {
        numberRows_ = 0;
        if (parent_->owner()) {
            delete parentBranch_;
            parentBranch_ = NULL;
        }
    }
}

#ifdef JJF_ZERO
// Constructor given parent
CbcNodeInfo::CbcNodeInfo (CbcNodeInfo * parent)
        :
        numberPointingToThis_(2),
        parent_(parent),
        parentBranch_(NULL),
        owner_(NULL),
        numberCuts_(0),
        nodeNumber_(0),
        cuts_(NULL),
        numberRows_(0),
        numberBranchesLeft_(2),
        active_(7)
{
#ifdef CHECK_NODE
    printf("CbcNodeInfo %x Constructor from parent %x\n", this, parent_);
#endif
    //setParentBasedData();
}
#endif

// Copy Constructor
CbcNodeInfo::CbcNodeInfo (const CbcNodeInfo & rhs)
        :
        numberPointingToThis_(rhs.numberPointingToThis_),
        parent_(rhs.parent_),
        parentBranch_(NULL),
        owner_(rhs.owner_),
        numberCuts_(rhs.numberCuts_),
        nodeNumber_(rhs.nodeNumber_),
        cuts_(NULL),
        numberRows_(rhs.numberRows_),
        numberBranchesLeft_(rhs.numberBranchesLeft_),
        active_(rhs.active_)
{
#ifdef CHECK_NODE
    printf("CbcNodeInfo %x Copy constructor\n", this);
#endif
    if (numberCuts_) {
        cuts_ = new CbcCountRowCut * [numberCuts_];
        int n = 0;
        for (int i = 0; i < numberCuts_; i++) {
            CbcCountRowCut * thisCut = rhs.cuts_[i];
            if (thisCut) {
                // I think this is correct - new one should take priority
                thisCut->setInfo(this, n);
                thisCut->increment(numberBranchesLeft_);
                cuts_[n++] = thisCut;
            }
        }
        numberCuts_ = n;
    }
    if (rhs.parentBranch_) {
        parentBranch_ = rhs.parentBranch_->clone();
    }
}
// Constructor given parent and owner
CbcNodeInfo::CbcNodeInfo (CbcNodeInfo * parent, CbcNode * owner)
        :
        numberPointingToThis_(2),
        parent_(parent),
        parentBranch_(NULL),
        owner_(owner),
        numberCuts_(0),
        nodeNumber_(0),
        cuts_(NULL),
        numberRows_(0),
        numberBranchesLeft_(2),
        active_(7)
{
#ifdef CHECK_NODE
    printf("CbcNodeInfo %x Constructor from parent %x\n", this, parent_);
#endif
    //setParentBasedData();
}

/**
   Take care to detach from the owning CbcNode and decrement the reference
   count in the parent.  If this is the last nodeInfo object pointing to the
   parent, make a recursive call to delete the parent.
*/
CbcNodeInfo::~CbcNodeInfo()
{
#ifdef CHECK_NODE
    printf("CbcNodeInfo %x Destructor parent %x\n", this, parent_);
#endif

    assert(!numberPointingToThis_);
    // But there may be some left (max nodes?)
    for (int i = 0; i < numberCuts_; i++) {
        if (cuts_[i]) {
#ifndef GLOBAL_CUTS_JUST_POINTERS
            delete cuts_[i];
#else
            if (cuts_[i]->globallyValidAsInteger() != 2)
                delete cuts_[i];
#endif
        }
    }
    delete [] cuts_;
    if (owner_)
        owner_->nullNodeInfo();
    if (parent_) {
        int numberLinks = parent_->decrement();
        if (!numberLinks) delete parent_;
    }
    delete parentBranch_;
}


//#define ALLCUTS
void
CbcNodeInfo::decrementCuts(int change)
{
    int i;
    // get rid of all remaining if negative
    int changeThis;
    if (change < 0)
        changeThis = numberBranchesLeft_;
    else
        changeThis = change;
    // decrement cut counts
    for (i = 0; i < numberCuts_; i++) {
        if (cuts_[i]) {
            int number = cuts_[i]->decrement(changeThis);
            if (!number) {
                //printf("info %x del cut %d %x\n",this,i,cuts_[i]);
#ifndef GLOBAL_CUTS_JUST_POINTERS
                delete cuts_[i];
#else
                if (cuts_[i]->globallyValidAsInteger() != 2)
                    delete cuts_[i];
#endif
                cuts_[i] = NULL;
            }
        }
    }
}
void
CbcNodeInfo::incrementCuts(int change)
{
    int i;
    assert (change > 0);
    // increment cut counts
    for (i = 0; i < numberCuts_; i++) {
        if (cuts_[i])
            cuts_[i]->increment(change);
    }
}
void
CbcNodeInfo::decrementParentCuts(CbcModel * model, int change)
{
    if (parent_) {
        // get rid of all remaining if negative
        int changeThis;
        if (change < 0)
            changeThis = numberBranchesLeft_;
        else
            changeThis = change;
        int i;
        // Get over-estimate of space needed for basis
        CoinWarmStartBasis & dummy = model->workingBasis();
        dummy.setSize(0, numberRows_ + numberCuts_);
        buildRowBasis(dummy);
        /* everything is zero (i.e. free) so we can use to see
           if latest basis */
        CbcNodeInfo * thisInfo = parent_;
        while (thisInfo)
            thisInfo = thisInfo->buildRowBasis(dummy);
        // decrement cut counts
        thisInfo = parent_;
        int numberRows = numberRows_;
        while (thisInfo) {
            for (i = thisInfo->numberCuts_ - 1; i >= 0; i--) {
                CoinWarmStartBasis::Status status = dummy.getArtifStatus(--numberRows);
#ifdef ALLCUTS
                status = CoinWarmStartBasis::isFree;
#endif
                if (thisInfo->cuts_[i]) {
                    int number = 1;
                    if (status != CoinWarmStartBasis::basic) {
                        // tight - drop 1 or 2
                        if (change < 0)
                            number = thisInfo->cuts_[i]->decrement(changeThis);
                        else
                            number = thisInfo->cuts_[i]->decrement(change);
                    }
                    if (!number) {
#ifndef GLOBAL_CUTS_JUST_POINTERS
                        delete thisInfo->cuts_[i];
#else
                        if (thisInfo->cuts_[i]->globallyValidAsInteger() != 2)
                            delete thisInfo->cuts_[i];
#endif
                        thisInfo->cuts_[i] = NULL;
                    }
                }
            }
            thisInfo = thisInfo->parent_;
        }
    }
}
#ifdef JJF_ZERO
void
CbcNodeInfo::incrementParentCuts(CbcModel * model, int change)
{
    if (parent_) {
        int i;
        // Get over-estimate of space needed for basis
        CoinWarmStartBasis & dummy = model->workingBasis();
        dummy.setSize(0, numberRows_ + numberCuts_);
        /* everything is zero (i.e. free) so we can use to see
           if latest basis */
        buildRowBasis(dummy);
        CbcNodeInfo * thisInfo = parent_;
        while (thisInfo)
            thisInfo = thisInfo->buildRowBasis(dummy);
        // increment cut counts
        thisInfo = parent_;
        int numberRows = numberRows_;
        while (thisInfo) {
            for (i = thisInfo->numberCuts_ - 1; i >= 0; i--) {
                CoinWarmStartBasis::Status status = dummy.getArtifStatus(--numberRows);
#ifdef ALLCUTS
                status = CoinWarmStartBasis::isFree;
#endif
                if (thisInfo->cuts_[i] && status != CoinWarmStartBasis::basic) {
                    thisInfo->cuts_[i]->increment(change);
                }
            }
            thisInfo = thisInfo->parent_;
        }
    }
}
#endif
/*
  Append cuts to the cuts_ array in a nodeInfo. The initial reference count
  is set to numberToBranchOn, which will normally be the number of arms
  defined for the CbcBranchingObject attached to the CbcNode that owns this
  CbcNodeInfo.
*/
void
CbcNodeInfo::addCuts (OsiCuts & cuts, int numberToBranchOn,
                      /*int * whichGenerator,*/int numberPointingToThis)
{
    int numberCuts = cuts.sizeRowCuts();
    if (numberCuts) {
        int i;
        if (!numberCuts_) {
            cuts_ = new CbcCountRowCut * [numberCuts];
        } else {
            CbcCountRowCut ** temp = new CbcCountRowCut * [numberCuts+numberCuts_];
            memcpy(temp, cuts_, numberCuts_*sizeof(CbcCountRowCut *));
            delete [] cuts_;
            cuts_ = temp;
        }
        for (i = 0; i < numberCuts; i++) {
            CbcCountRowCut * thisCut = new CbcCountRowCut(*cuts.rowCutPtr(i),
                    this, numberCuts_,
                    -1, numberPointingToThis);
            thisCut->increment(numberToBranchOn);
            cuts_[numberCuts_++] = thisCut;
#ifdef CBC_DEBUG
#if CBC_DEBUG>1
            int n = thisCut->row().getNumElements();
            printf("Cut %d has %d entries, rhs %g %g =>", i, n, thisCut->lb(),
                   thisCut->ub());
            int j;
            const int * index = thisCut->row().getIndices();
            const double * element = thisCut->row().getElements();
            for (j = 0; j < n; j++) {
                printf(" (%d,%g)", index[j], element[j]);
                assert(fabs(element[j]) > 1.00e-12);
            }
            printf("\n");
#else
            int n = thisCut->row().getNumElements();
            int j;
            const double * element = thisCut->row().getElements();
            for (j = 0; j < n; j++) {
                assert(fabs(element[j]) > 1.00e-12);
            }
#endif
#endif
        }
    }
}

void
CbcNodeInfo::addCuts(int numberCuts, CbcCountRowCut ** cut,
                     int numberToBranchOn)
{
    if (numberCuts) {
        int i;
        if (!numberCuts_) {
            cuts_ = new CbcCountRowCut * [numberCuts];
        } else {
            CbcCountRowCut ** temp = new CbcCountRowCut * [numberCuts+numberCuts_];
            memcpy(temp, cuts_, numberCuts_*sizeof(CbcCountRowCut *));
            delete [] cuts_;
            cuts_ = temp;
        }
        for (i = 0; i < numberCuts; i++) {
            CbcCountRowCut * thisCut = cut[i];
            thisCut->setInfo(this, numberCuts_);
            //printf("info %x cut %d %x\n",this,i,thisCut);
            thisCut->increment(numberToBranchOn);
            cuts_[numberCuts_++] = thisCut;
#ifdef CBC_DEBUG
            int n = thisCut->row().getNumElements();
#if CBC_DEBUG>1
            printf("Cut %d has %d entries, rhs %g %g =>", i, n, thisCut->lb(),
                   thisCut->ub());
#endif
            int j;
#if CBC_DEBUG>1
            const int * index = thisCut->row().getIndices();
#endif
            const double * element = thisCut->row().getElements();
            for (j = 0; j < n; j++) {
#if CBC_DEBUG>1
                printf(" (%d,%g)", index[j], element[j]);
#endif
                assert(fabs(element[j]) > 1.00e-12);
            }
            printf("\n");
#endif
        }
    }
}

// delete cuts
void
CbcNodeInfo::deleteCuts(int numberToDelete, CbcCountRowCut ** cuts)
{
    int i;
    int j;
    int last = -1;
    for (i = 0; i < numberToDelete; i++) {
        CbcCountRowCut * next = cuts[i];
        for (j = last + 1; j < numberCuts_; j++) {
            if (next == cuts_[j])
                break;
        }
        if (j == numberCuts_) {
            // start from beginning
            for (j = 0; j < last; j++) {
                if (next == cuts_[j])
                    break;
            }
            assert(j < last);
        }
        last = j;
        int number = cuts_[j]->decrement();
        if (!number) {
#ifndef GLOBAL_CUTS_JUST_POINTERS
            delete cuts_[j];
#else
            if (cuts_[j]->globallyValidAsInteger() != 2)
                delete cuts_[j];
#endif
        }
        cuts_[j] = NULL;
    }
    j = 0;
    for (i = 0; i < numberCuts_; i++) {
        if (cuts_[i])
            cuts_[j++] = cuts_[i];
    }
    numberCuts_ = j;
}

// delete cuts
void
CbcNodeInfo::deleteCuts(int numberToDelete, int * which)
{
    int i;
    for (i = 0; i < numberToDelete; i++) {
        int iCut = which[i];
        int number = cuts_[iCut]->decrement();
        if (!number) {
#ifndef GLOBAL_CUTS_JUST_POINTERS
            delete cuts_[iCut];
#else
            if (cuts_[iCut]->globallyValidAsInteger() != 2)
                delete cuts_[iCut];
#endif
        }
        cuts_[iCut] = NULL;
    }
    int j = 0;
    for (i = 0; i < numberCuts_; i++) {
        if (cuts_[i])
            cuts_[j++] = cuts_[i];
    }
    numberCuts_ = j;
}

// Really delete a cut
void
CbcNodeInfo::deleteCut(int whichOne)
{
    assert(whichOne < numberCuts_);
    cuts_[whichOne] = NULL;
}
/* Deactivate node information.
   1 - bounds
   2 - cuts
   4 - basis!
*/
void
CbcNodeInfo::deactivate(int mode)
{
    active_ &= (~mode);
}

