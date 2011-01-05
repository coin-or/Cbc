// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/12/2009 carved from CbcBranchBase

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "OsiChooseVariable.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchBase.hpp"
#include "CbcBranchDecision.hpp"

// Default Constructor
CbcBranchDecision::CbcBranchDecision ()
        : object_(NULL), model_(NULL), chooseMethod_(NULL)
{
}

// Copy Constructor
CbcBranchDecision::CbcBranchDecision (const CbcBranchDecision &rhs)
        : object_(NULL), model_(rhs.model_), chooseMethod_(NULL)
{
    if (rhs.chooseMethod_)
        chooseMethod_ = rhs.chooseMethod_->clone();
}

CbcBranchDecision::~CbcBranchDecision()
{
    delete object_;
    delete chooseMethod_;
}
/* Compare N branching objects. Return index of best
   and sets way of branching in chosen object.

   This routine is used only after strong branching.
   This is reccommended version as it can be more sophisticated
*/

int
CbcBranchDecision::bestBranch (CbcBranchingObject ** objects, int numberObjects,
                               int /*numberUnsatisfied*/,
                               double * changeUp, int * numberInfeasibilitiesUp,
                               double * changeDown, int * numberInfeasibilitiesDown,
                               double /*objectiveValue*/)
{
    int bestWay = 0;
    int whichObject = -1;
    if (numberObjects) {
        initialize(objects[0]->model());
        CbcBranchingObject * bestObject = NULL;
        for (int i = 0 ; i < numberObjects ; i++) {
            int betterWay = betterBranch(objects[i],
                                         bestObject,
                                         changeUp[i],
                                         numberInfeasibilitiesUp [i],
                                         changeDown[i],
                                         numberInfeasibilitiesDown[i] );
            if (betterWay) {
                bestObject = objects[i];
                bestWay = betterWay;
                whichObject = i;
            }
        }
        // set way in best
        if (whichObject >= 0)
            objects[whichObject]->way(bestWay);
    }
    return whichObject;
}
// Set (clone) chooseMethod
void
CbcBranchDecision::setChooseMethod(const OsiChooseVariable & method)
{
    delete chooseMethod_;
    chooseMethod_ = method.clone();
}

