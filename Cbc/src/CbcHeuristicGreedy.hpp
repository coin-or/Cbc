/* $Id$ */
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcHeuristicGreedy_H
#define CbcHeuristicGreedy_H

#include "CbcHeuristic.hpp"
/** Greedy heuristic classes
 */

class CbcHeuristicGreedyCover : public CbcHeuristic {
public:

    // Default Constructor
    CbcHeuristicGreedyCover ();

    /* Constructor with model - assumed before cuts
       Initial version does not do Lps
    */
    CbcHeuristicGreedyCover (CbcModel & model);

    // Copy constructor
    CbcHeuristicGreedyCover ( const CbcHeuristicGreedyCover &);

    // Destructor
    ~CbcHeuristicGreedyCover ();

    /// Clone
    virtual CbcHeuristic * clone() const;
    /// Assignment operator
    CbcHeuristicGreedyCover & operator=(const CbcHeuristicGreedyCover& rhs);
    /// Create C++ lines to get to current state
    virtual void generateCpp( FILE * fp) ;

    /// update model (This is needed if cliques update matrix etc)
    virtual void setModel(CbcModel * model);

    using CbcHeuristic::solution ;
    /** returns 0 if no solution, 1 if valid solution.
        Sets solution values if good, sets objective value (only if good)
        We leave all variables which are at one at this node of the
        tree to that value and will
        initially set all others to zero.  We then sort all variables in order of their cost
        divided by the number of entries in rows which are not yet covered.  We randomize that
        value a bit so that ties will be broken in different ways on different runs of the heuristic.
        We then choose the best one and set it to one and repeat the exercise.

    */
    virtual int solution(double & objectiveValue,
                         double * newSolution);
    /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
    virtual void validate() ;
    /// Resets stuff if model changes
    virtual void resetModel(CbcModel * model);
    /* Algorithm
       0 - use current upper bounds
       1 - use original upper bounds
       If 10 added perturb ratios more
       if 100 added round up all >=0.5
    */
    inline int algorithm() const {
        return algorithm_;
    }
    inline void setAlgorithm(int value) {
        algorithm_ = value;
    }
    // Only do this many times
    inline int numberTimes() const {
        return numberTimes_;
    }
    inline void setNumberTimes(int value) {
        numberTimes_ = value;
    }

protected:
    /// Guts of constructor from a CbcModel
    void gutsOfConstructor(CbcModel * model);
    // Data

    // Original matrix by column
    CoinPackedMatrix matrix_;
    // original number of rows
    int originalNumberRows_;
    /* Algorithm
       0 - use current upper bounds
       1 - use original upper bounds
       If 10 added perturb ratios more
    */
    int algorithm_;
    /// Do this many times
    int numberTimes_;

};


class CbcHeuristicGreedyEquality : public CbcHeuristic {
public:

    // Default Constructor
    CbcHeuristicGreedyEquality ();

    /* Constructor with model - assumed before cuts
       Initial version does not do Lps
    */
    CbcHeuristicGreedyEquality (CbcModel & model);

    // Copy constructor
    CbcHeuristicGreedyEquality ( const CbcHeuristicGreedyEquality &);

    // Destructor
    ~CbcHeuristicGreedyEquality ();

    /// Clone
    virtual CbcHeuristic * clone() const;
    /// Assignment operator
    CbcHeuristicGreedyEquality & operator=(const CbcHeuristicGreedyEquality& rhs);
    /// Create C++ lines to get to current state
    virtual void generateCpp( FILE * fp) ;

    /// update model (This is needed if cliques update matrix etc)
    virtual void setModel(CbcModel * model);

    using CbcHeuristic::solution ;
    /** returns 0 if no solution, 1 if valid solution.
        Sets solution values if good, sets objective value (only if good)
        We leave all variables which are at one at this node of the
        tree to that value and will
        initially set all others to zero.  We then sort all variables in order of their cost
        divided by the number of entries in rows which are not yet covered.  We randomize that
        value a bit so that ties will be broken in different ways on different runs of the heuristic.
        We then choose the best one and set it to one and repeat the exercise.

    */
    virtual int solution(double & objectiveValue,
                         double * newSolution);
    /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
    virtual void validate() ;
    /// Resets stuff if model changes
    virtual void resetModel(CbcModel * model);
    /* Algorithm
       0 - use current upper bounds
       1 - use original upper bounds
       If 10 added perturb ratios more
       if 100 added round up all >=0.5
    */
    inline int algorithm() const {
        return algorithm_;
    }
    inline void setAlgorithm(int value) {
        algorithm_ = value;
    }
    // Fraction of rhs to cover before branch and cut
    inline void setFraction(double value) {
        fraction_ = value;
    }
    inline double fraction() const {
        return fraction_;
    }
    // Only do this many times
    inline int numberTimes() const {
        return numberTimes_;
    }
    inline void setNumberTimes(int value) {
        numberTimes_ = value;
    }
protected:
    /// Guts of constructor from a CbcModel
    void gutsOfConstructor(CbcModel * model);
    // Data

    // Original matrix by column
    CoinPackedMatrix matrix_;
    // Fraction of rhs to cover before branch and cut
    double fraction_;
    // original number of rows
    int originalNumberRows_;
    /* Algorithm
       0 - use current upper bounds
       1 - use original upper bounds
       If 10 added perturb ratios more
    */
    int algorithm_;
    /// Do this many times
    int numberTimes_;

};


#endif

