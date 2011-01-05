// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CoinTypes.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcGeneral.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

// Default Constructor
CbcGeneral::CbcGeneral()
        : CbcObject()
{
}

// Constructor from model
CbcGeneral::CbcGeneral(CbcModel * model)
        : CbcObject(model)
{
}


// Destructor
CbcGeneral::~CbcGeneral ()
{
}

// Copy constructor
CbcGeneral::CbcGeneral ( const CbcGeneral & rhs)
        : CbcObject(rhs)
{
}
#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"
#include "ClpNode.hpp"
#include "CbcBranchDynamic.hpp"
// Assignment operator
CbcGeneral &
CbcGeneral::operator=( const CbcGeneral & rhs)
{
    if (this != &rhs) {
        CbcObject::operator=(rhs);
    }
    return *this;
}
// Infeasibility - large is 0.5
double
CbcGeneral::infeasibility(const OsiBranchingInformation * /*info*/,
                          int &/*preferredWay*/) const
{
    abort();
    return 0.0;
}
CbcBranchingObject *
CbcGeneral::createCbcBranch(OsiSolverInterface * /*solver*/, const OsiBranchingInformation * /*info*/, int /*way*/)
{
    abort();
    return NULL;
}
#endif

