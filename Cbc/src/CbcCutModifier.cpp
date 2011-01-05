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
#include "CbcCutModifier.hpp"

// Default Constructor
CbcCutModifier::CbcCutModifier()
{
}


// Destructor
CbcCutModifier::~CbcCutModifier ()
{
}

// Copy constructor
CbcCutModifier::CbcCutModifier ( const CbcCutModifier & /*rhs*/)
{
}

// Assignment operator
CbcCutModifier &
CbcCutModifier::operator=( const CbcCutModifier & rhs)
{
    if (this != &rhs) {
    }
    return *this;
}

