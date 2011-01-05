// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

//Edwin 11/25/09 carved out of CbcCompareActual

#ifndef CbcCompareObjective_H
#define CbcCompareObjective_H


//#############################################################################
/*  These are alternative strategies for node traversal.
    They can take data etc for fine tuning

    At present the node list is stored as a heap and the "test"
    comparison function returns true if node y is better than node x.

*/
#include "CbcNode.hpp"
#include "CbcCompareBase.hpp"
#include "CbcCompare.hpp"

class CbcModel;

class CbcCompareObjective  : public CbcCompareBase {
public:
    // Default Constructor
    CbcCompareObjective ();

    virtual ~CbcCompareObjective();
    // Copy constructor
    CbcCompareObjective ( const CbcCompareObjective &rhs);

    // Assignment operator
    CbcCompareObjective & operator=( const CbcCompareObjective& rhs);

    /// Clone
    virtual CbcCompareBase * clone() const;
    /// Create C++ lines to get to current state
    virtual void generateCpp( FILE * fp);

    /* This returns true if objective value of node y is less than
       objective value of node x */
    virtual bool test (CbcNode * x, CbcNode * y);
};

#endif //CbcCompareObjective_H

