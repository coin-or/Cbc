// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/12/2009 carved from CbcBranchBase

#ifndef CbcObjectUpdateData_H
#define CbcObjectUpdateData_H

#include "CbcObject.hpp"
/*  This stores data so an object can be updated
 */
class CbcObjectUpdateData {

public:

    /// Default Constructor
    CbcObjectUpdateData ();

    /// Useful constructor
    CbcObjectUpdateData (CbcObject * object,
                         int way,
                         double change,
                         int status,
                         int intDecrease_,
                         double branchingValue);

    /// Copy constructor
    CbcObjectUpdateData ( const CbcObjectUpdateData &);

    /// Assignment operator
    CbcObjectUpdateData & operator=( const CbcObjectUpdateData& rhs);

    /// Destructor
    virtual ~CbcObjectUpdateData ();


public:
    /// data

    /// Object
    CbcObject * object_;
    /// Branch as defined by instance of CbcObject
    int way_;
    /// Object number
    int objectNumber_;
    /// Change in objective
    double change_;
    /// Status 0 Optimal, 1 infeasible, 2 unknown
    int status_;
    /// Decrease in number unsatisfied
    int intDecrease_;
    /// Branching value
    double branchingValue_;
    /// Objective value before branching
    double originalObjective_;
    /// Current cutoff
    double cutoff_;

};

#endif

