// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/13/2009-- carved out of CbcBranchCut

#ifndef CbcBranchAllDifferent_H
#define CbcBranchAllDifferent_H

#include "CbcBranchBase.hpp"
#include "OsiRowCut.hpp"
#include "CoinPackedMatrix.hpp"
#include "CbcBranchCut.hpp"

/** Define a branch class that branches so that it is only satsified if all
    members have different values
    So cut is x <= y-1 or x >= y+1
*/


class CbcBranchAllDifferent : public CbcBranchCut {

public:

    // Default Constructor
    CbcBranchAllDifferent ();

    /** Useful constructor - passed set of integer variables which must all be different
    */
    CbcBranchAllDifferent (CbcModel * model, int number, const int * which);

    // Copy constructor
    CbcBranchAllDifferent ( const CbcBranchAllDifferent &);

    /// Clone
    virtual CbcObject * clone() const;

    // Assignment operator
    CbcBranchAllDifferent & operator=( const CbcBranchAllDifferent& rhs);

    // Destructor
    ~CbcBranchAllDifferent ();

    /// Infeasibility - large is 0.5
    virtual double infeasibility(const OsiBranchingInformation * info,
                                 int &preferredWay) const;

    /// Creates a branching object
    virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) ;


protected:
    /// data

    /// Number of entries
    int numberInSet_;
    /// Which variables
    int * which_;
};
#endif

