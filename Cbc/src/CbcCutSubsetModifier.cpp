// $Id$
// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

//Edwin 11/25/09 carved out of CbcCutGenerator

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include "CbcConfig.h"
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#else
#include "OsiSolverInterface.hpp"
#endif
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcBranchDynamic.hpp"
#include "CglProbing.hpp"
#include "CoinTime.hpp"
#include "CbcCutSubsetModifier.hpp"

// Default Constructor
CbcCutSubsetModifier::CbcCutSubsetModifier ()
        : CbcCutModifier(),
        firstOdd_(COIN_INT_MAX)
{
}

// Useful constructor
CbcCutSubsetModifier::CbcCutSubsetModifier (int firstOdd)
        : CbcCutModifier()
{
    firstOdd_ = firstOdd;
}

// Copy constructor
CbcCutSubsetModifier::CbcCutSubsetModifier ( const CbcCutSubsetModifier & rhs)
        : CbcCutModifier(rhs)
{
    firstOdd_ = rhs.firstOdd_;
}

// Clone
CbcCutModifier *
CbcCutSubsetModifier::clone() const
{
    return new CbcCutSubsetModifier(*this);
}

// Assignment operator
CbcCutSubsetModifier &
CbcCutSubsetModifier::operator=( const CbcCutSubsetModifier & rhs)
{
    if (this != &rhs) {
        CbcCutModifier::operator=(rhs);
        firstOdd_ = rhs.firstOdd_;
    }
    return *this;
}

// Destructor
CbcCutSubsetModifier::~CbcCutSubsetModifier ()
{
}
/* Returns
   0 unchanged
   1 strengthened
   2 weakened
   3 deleted
*/
int
CbcCutSubsetModifier::modify(const OsiSolverInterface * /*solver*/,
                             OsiRowCut & cut)
{
    int n = cut.row().getNumElements();
    if (!n)
        return 0;
    const int * column = cut.row().getIndices();
    //const double * element = cut.row().getElements();
    int returnCode = 0;
    for (int i = 0; i < n; i++) {
        if (column[i] >= firstOdd_) {
            returnCode = 3;
            break;
        }
    }
#ifdef COIN_DETAIL
    if (!returnCode) {
        const double * element = cut.row().getElements();
        printf("%g <= ", cut.lb());
        for (int i = 0; i < n; i++) {
            printf("%g*x%d ", element[i], column[i]);
        }
        printf("<= %g\n", cut.ub());
    }
#endif
    //return 3;
    return returnCode;
}

