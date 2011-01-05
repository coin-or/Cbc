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
#include "CbcCompareObjective.hpp"
/** Default Constructor

*/
CbcCompareObjective::CbcCompareObjective ()
        : CbcCompareBase()
{
    test_ = this;
}

// Copy constructor
CbcCompareObjective::CbcCompareObjective ( const CbcCompareObjective & rhs)
        : CbcCompareBase(rhs)

{
}

// Clone
CbcCompareBase *
CbcCompareObjective::clone() const
{
    return new CbcCompareObjective(*this);
}

// Assignment operator
CbcCompareObjective &
CbcCompareObjective::operator=( const CbcCompareObjective & rhs)
{
    if (this != &rhs) {
        CbcCompareBase::operator=(rhs);
    }
    return *this;
}

// Destructor
CbcCompareObjective::~CbcCompareObjective ()
{
}

// Returns true if y better than x
bool
CbcCompareObjective::test (CbcNode * x, CbcNode * y)
{
    double testX = x->objectiveValue();
    double testY = y->objectiveValue();
    if (testX != testY)
        return testX > testY;
    else
        return equalityTest(x, y); // so ties will be broken in consistent manner
}
// Create C++ lines to get to current state
void
CbcCompareObjective::generateCpp( FILE * fp)
{
    fprintf(fp, "0#include \"CbcCompareActual.hpp\"\n");
    fprintf(fp, "3  CbcCompareObjective compare;\n");
    fprintf(fp, "3  cbcModel->setNodeComparison(compare);\n");
}
