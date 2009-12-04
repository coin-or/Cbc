/* $Id: CbcHeuristicPivotAndComplement.hpp 1173 2009-06-04 09:44:10Z forrest $ */
// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcHeuristicPivotAndComplement_H
#define CbcHeuristicPivotAndComplement_H

#include "CbcHeuristic.hpp"
/** LocalSearch class
 */

class CbcHeuristicPivotAndComplement : public CbcHeuristic {
public:

    // Default Constructor
    CbcHeuristicPivotAndComplement ();

    /* Constructor with model - assumed before cuts
       Initial version does not do Lps
    */
    CbcHeuristicPivotAndComplement (CbcModel & model);

    // Copy constructor
    CbcHeuristicPivotAndComplement ( const CbcHeuristicPivotAndComplement &);

    // Destructor
    ~CbcHeuristicPivotAndComplement ();

    /// Clone
    virtual CbcHeuristic * clone() const;

    /// Assignment operator
    CbcHeuristicPivotAndComplement & operator=(const CbcHeuristicPivotAndComplement& rhs);

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
