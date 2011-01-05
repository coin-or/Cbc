// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

//Edwin 11/24/09 carved out of CbcCompareActual

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
#include "CbcCompareDepth.hpp"
/** Default Constructor

*/
CbcCompareDepth::CbcCompareDepth ()
        : CbcCompareBase()
{
    test_ = this;
}

// Copy constructor
CbcCompareDepth::CbcCompareDepth ( const CbcCompareDepth & rhs)
        : CbcCompareBase(rhs)

{
}

// Clone
CbcCompareBase *
CbcCompareDepth::clone() const
{
    return new CbcCompareDepth(*this);
}

// Assignment operator
CbcCompareDepth &
CbcCompareDepth::operator=( const CbcCompareDepth & rhs)
{
    if (this != &rhs) {
        CbcCompareBase::operator=(rhs);
    }
    return *this;
}

// Destructor
CbcCompareDepth::~CbcCompareDepth ()
{
}

// Returns true if y better than x
bool
CbcCompareDepth::test (CbcNode * x, CbcNode * y)
{
    int testX = x->depth();
    int testY = y->depth();
    if (testX != testY)
        return testX < testY;
    else
        return equalityTest(x, y); // so ties will be broken in consistent manner
}
// Create C++ lines to get to current state
void
CbcCompareDepth::generateCpp( FILE * fp)
{
    fprintf(fp, "0#include \"CbcCompareActual.hpp\"\n");
    fprintf(fp, "3  CbcCompareDepth compare;\n");
    fprintf(fp, "3  cbcModel->setNodeComparison(compare);\n");
}

