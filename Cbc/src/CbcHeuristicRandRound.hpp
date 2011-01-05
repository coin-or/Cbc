/* $Id: CbcHeuristicRandRound.hpp 1173 2009-06-04 09:44:10Z forrest $ */
// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcHeuristicRandRound_H
#define CbcHeuristicRandRound_H

#include "CbcHeuristic.hpp"
/** LocalSearch class
 */

class CbcHeuristicRandRound : public CbcHeuristic {
public:

    // Default Constructor
    CbcHeuristicRandRound ();

    /* Constructor with model - assumed before cuts
       Initial version does not do Lps
    */
    CbcHeuristicRandRound (CbcModel & model);

    // Copy constructor
    CbcHeuristicRandRound ( const CbcHeuristicRandRound &);

    // Destructor
    ~CbcHeuristicRandRound ();

    /// Clone
    virtual CbcHeuristic * clone() const;

    /// Assignment operator
    CbcHeuristicRandRound & operator=(const CbcHeuristicRandRound& rhs);

    /// Create C++ lines to get to current state
    virtual void generateCpp( FILE * fp) ;

    /// Resets stuff if model changes
    virtual void resetModel(CbcModel * model);

    /// update model (This is needed if cliques update matrix etc)
    virtual void setModel(CbcModel * model);

    using CbcHeuristic::solution ;
    /** returns 0 if no solution, 1 if valid solution.
        Sets solution values if good, sets objective value (only if good)
        needs comments
    */
    virtual int solution(double & objectiveValue,
                         double * newSolution);

protected:
};


#endif

