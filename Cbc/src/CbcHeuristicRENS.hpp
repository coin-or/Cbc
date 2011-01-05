// $Id$
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// edwin 12/5/09 carved out of CbcHeuristicRINS

#ifndef CbcHeuristicRENS_H
#define CbcHeuristicRENS_H

#include "CbcHeuristic.hpp"

/** LocalSearch class
 */

class CbcHeuristicRENS : public CbcHeuristic {
public:

    // Default Constructor
    CbcHeuristicRENS ();

    /* Constructor with model - assumed before cuts
       Initial version does not do Lps
    */
    CbcHeuristicRENS (CbcModel & model);

    // Copy constructor
    CbcHeuristicRENS ( const CbcHeuristicRENS &);

    // Destructor
    ~CbcHeuristicRENS ();

    /// Clone
    virtual CbcHeuristic * clone() const;


    /// Assignment operator
    CbcHeuristicRENS & operator=(const CbcHeuristicRENS& rhs);

    /// Resets stuff if model changes
    virtual void resetModel(CbcModel * model);

    /// update model (This is needed if cliques update matrix etc)
    virtual void setModel(CbcModel * model);

    using CbcHeuristic::solution ;
    /** returns 0 if no solution, 1 if valid solution.
        Sets solution values if good, sets objective value (only if good)
        This does Relaxation Extension Neighborhood Search
        Does not run if when_<2 and a solution exists
    */
    virtual int solution(double & objectiveValue,
                         double * newSolution);

    /// Set type
    inline void setRensType(int value)
    { rensType_ = value;}

protected:
    // Data
    /// Number of tries
    int numberTries_;
    /** Type
        0 - fix at LB
        1 - fix on dj
        2 - fix at UB as well
	3 - fix on 0.01*average dj
	add 16 to allow two tries
    */
    int rensType_;
};

#endif

