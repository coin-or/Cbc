// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

//Edwin 11/25/09 carved out of CbcCompareActual

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CbcMessage.hpp"
#include "CbcModel.hpp"
#include "CbcTree.hpp"
#include "CbcCompareActual.hpp"
#include "CoinError.hpp"
#include "CbcCompareEstimate.hpp"
/** Default Constructor

*/
CbcCompareEstimate::CbcCompareEstimate ()
        : CbcCompareBase()
{
    test_ = this;
}

// Copy constructor
CbcCompareEstimate::CbcCompareEstimate ( const CbcCompareEstimate & rhs)
        : CbcCompareBase(rhs)

{
}

// Clone
CbcCompareBase *
CbcCompareEstimate::clone() const
{
    return new CbcCompareEstimate(*this);
}

// Assignment operator
CbcCompareEstimate &
CbcCompareEstimate::operator=( const CbcCompareEstimate & rhs)
{
    if (this != &rhs) {
        CbcCompareBase::operator=(rhs);
    }
    return *this;
}

// Destructor
CbcCompareEstimate::~CbcCompareEstimate ()
{
}

// Returns true if y better than x
bool
CbcCompareEstimate::test (CbcNode * x, CbcNode * y)
{
    double testX = x->guessedObjectiveValue();
    double testY = y->guessedObjectiveValue();
    if (testX != testY)
        return testX > testY;
    else
        return equalityTest(x, y); // so ties will be broken in consistent manner
}

// Create C++ lines to get to current state
void
CbcCompareEstimate::generateCpp( FILE * fp)
{
    fprintf(fp, "0#include \"CbcCompareActual.hpp\"\n");
    fprintf(fp, "3  CbcCompareEstimate compare;\n");
    fprintf(fp, "3  cbcModel->setNodeComparison(compare);\n");
}

