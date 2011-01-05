// $Id$
// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

//Edwin 11/25/09 carved out of CbcCutGenerator

#ifndef CbcCutModifier_H
#define CbcCutModifier_H

#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "CglCutGenerator.hpp"

class CbcModel;
class OsiRowCut;
class OsiRowCutDebugger;
/** Abstract cut modifier base class

    In exotic circumstances - cuts may need to be modified
    a) strengthened - changed
    b) weakened - changed
    c) deleted - set to NULL
    d) unchanged
*/

class CbcCutModifier {
public:
    /// Default Constructor
    CbcCutModifier ();

    // Copy constructor
    CbcCutModifier ( const CbcCutModifier &);

    /// Destructor
    virtual ~CbcCutModifier();

    /// Assignment
    CbcCutModifier & operator=(const CbcCutModifier& rhs);
/// Clone
    virtual CbcCutModifier * clone() const = 0;

    /** Returns
        0 unchanged
        1 strengthened
        2 weakened
        3 deleted
    */
    virtual int modify(const OsiSolverInterface * solver, OsiRowCut & cut) = 0;
    /// Create C++ lines to get to current state
    virtual void generateCpp( FILE * ) {}
protected:

};

#endif //CbcCutModifier_H

