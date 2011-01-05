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


// Default Constructor
CbcBranchingObject::CbcBranchingObject()
        : OsiBranchingObject()
{
    model_ = NULL;
    originalCbcObject_ = NULL;
    variable_ = -1;
    way_ = 0;
}

// Useful constructor
CbcBranchingObject::CbcBranchingObject (CbcModel * model, int variable, int way , double value)
        : OsiBranchingObject(model->solver(), value)
{
    model_ = model;
    originalCbcObject_ = NULL;
    variable_ = variable;
    way_ = way;
}

// Copy constructor
CbcBranchingObject::CbcBranchingObject ( const CbcBranchingObject & rhs)
        : OsiBranchingObject(rhs)
{
    model_ = rhs.model_;
    originalCbcObject_ = rhs.originalCbcObject_;
    variable_ = rhs.variable_;
    way_ = rhs.way_;
    value_ = rhs.value_;
}

// Assignment operator
CbcBranchingObject &
CbcBranchingObject::operator=( const CbcBranchingObject & rhs)
{
    if (this != &rhs) {
        OsiBranchingObject::operator=(rhs);
        model_ = rhs.model_;
        originalCbcObject_ = rhs.originalCbcObject_;
        variable_ = rhs.variable_;
        way_ = rhs.way_;
    }
    return *this;
}

// Destructor
CbcBranchingObject::~CbcBranchingObject ()
{
}

