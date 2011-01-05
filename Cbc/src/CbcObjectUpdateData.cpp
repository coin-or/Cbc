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
#include "CbcObjectUpdateData.hpp"

// Default constructor
CbcObjectUpdateData::CbcObjectUpdateData()
        : object_(NULL),
        way_(0),
        objectNumber_(-1),
        change_(0.0),
        status_(0),
        intDecrease_(0),
        branchingValue_(0.0),
        originalObjective_(COIN_DBL_MAX),
        cutoff_(COIN_DBL_MAX)
{
}

// Useful constructor
CbcObjectUpdateData::CbcObjectUpdateData (CbcObject * object,
        int way,
        double change,
        int status,
        int intDecrease,
        double branchingValue)
        : object_(object),
        way_(way),
        objectNumber_(-1),
        change_(change),
        status_(status),
        intDecrease_(intDecrease),
        branchingValue_(branchingValue),
        originalObjective_(COIN_DBL_MAX),
        cutoff_(COIN_DBL_MAX)
{
}

// Destructor
CbcObjectUpdateData::~CbcObjectUpdateData ()
{
}

// Copy constructor
CbcObjectUpdateData::CbcObjectUpdateData ( const CbcObjectUpdateData & rhs)
        : object_(rhs.object_),
        way_(rhs.way_),
        objectNumber_(rhs.objectNumber_),
        change_(rhs.change_),
        status_(rhs.status_),
        intDecrease_(rhs.intDecrease_),
        branchingValue_(rhs.branchingValue_),
        originalObjective_(rhs.originalObjective_),
        cutoff_(rhs.cutoff_)
{
}

// Assignment operator
CbcObjectUpdateData &
CbcObjectUpdateData::operator=( const CbcObjectUpdateData & rhs)
{
    if (this != &rhs) {
        object_ = rhs.object_;
        way_ = rhs.way_;
        objectNumber_ = rhs.objectNumber_;
        change_ = rhs.change_;
        status_ = rhs.status_;
        intDecrease_ = rhs.intDecrease_;
        branchingValue_ = rhs.branchingValue_;
        originalObjective_ = rhs.originalObjective_;
        cutoff_ = rhs.cutoff_;
    }
    return *this;
}

