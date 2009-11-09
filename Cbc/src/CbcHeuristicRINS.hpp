/* $Id: CbcHeuristicRINS.hpp 1261 2009-10-30 12:45:20Z forrest $ */
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcHeuristicRINS_H
#define CbcHeuristicRINS_H

#include "CbcHeuristic.hpp"
/** LocalSearch class
 */

class CbcHeuristicRINS : public CbcHeuristic {
public:

    // Default Constructor
    CbcHeuristicRINS ();

    /* Constructor with model - assumed before cuts
       Initial version does not do Lps
    */
    CbcHeuristicRINS (CbcModel & model);

    // Copy constructor
    CbcHeuristicRINS ( const CbcHeuristicRINS &);

    // Destructor
    ~CbcHeuristicRINS ();

    /// Clone
    virtual CbcHeuristic * clone() const;


    /// Assignment operator
    CbcHeuristicRINS & operator=(const CbcHeuristicRINS& rhs);

    /// Create C++ lines to get to current state
    virtual void generateCpp( FILE * fp) ;

    /// Resets stuff if model changes
    virtual void resetModel(CbcModel * model);

    /// update model (This is needed if cliques update matrix etc)
    virtual void setModel(CbcModel * model);

    using CbcHeuristic::solution ;
    /** returns 0 if no solution, 1 if valid solution.
        Sets solution values if good, sets objective value (only if good)
        This does Relaxation Induced Neighborhood Search
    */
    virtual int solution(double & objectiveValue,
                         double * newSolution);
    /// This version fixes stuff and does IP
    int solutionFix(double & objectiveValue,
                    double * newSolution,
                    const int * keep);

    /// Sets how often to do it
    inline void setHowOften(int value) {
        howOften_ = value;
    }
    /// Used array so we can set
    inline char * used() const {
        return used_;
    }

protected:
    // Data

    /// Number of solutions so we can do something at solution
    int numberSolutions_;
    /// How often to do (code can change)
    int howOften_;
    /// Number of successes
    int numberSuccesses_;
    /// Number of tries
    int numberTries_;
    /** State of fixing continuous variables -
        0 - not tried
        +n - this divisor makes small enough
        -n - this divisor still not small enough
    */
    int stateOfFixing_;
    /// Node when last done
    int lastNode_;
    /// Whether a variable has been in a solution
    char * used_;
};

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

protected:
    // Data
    /// Number of tries
    int numberTries_;
};

class CbcHeuristicDINS : public CbcHeuristic {
public:

    // Default Constructor
    CbcHeuristicDINS ();

    /* Constructor with model - assumed before cuts
       Initial version does not do Lps
    */
    CbcHeuristicDINS (CbcModel & model);

    // Copy constructor
    CbcHeuristicDINS ( const CbcHeuristicDINS &);

    // Destructor
    ~CbcHeuristicDINS ();

    /// Clone
    virtual CbcHeuristic * clone() const;


    /// Assignment operator
    CbcHeuristicDINS & operator=(const CbcHeuristicDINS& rhs);

    /// Create C++ lines to get to current state
    virtual void generateCpp( FILE * fp) ;

    /// Resets stuff if model changes
    virtual void resetModel(CbcModel * model);

    /// update model (This is needed if cliques update matrix etc)
    virtual void setModel(CbcModel * model);

    using CbcHeuristic::solution ;
    /** returns 0 if no solution, 1 if valid solution.
        Sets solution values if good, sets objective value (only if good)
        This does Relaxation Induced Neighborhood Search
    */
    virtual int solution(double & objectiveValue,
                         double * newSolution);
    /// This version fixes stuff and does IP
    int solutionFix(double & objectiveValue,
                    double * newSolution,
                    const int * keep);

    /// Sets how often to do it
    inline void setHowOften(int value) {
        howOften_ = value;
    }
    /// Sets maximum number of solutions kept
    inline void setMaximumKeep(int value) {
        maximumKeepSolutions_ = value;
    }
    /// Sets tightness of extra constraint
    inline void setConstraint(int value) {
        localSpace_ = value;
    }

protected:
    // Data

    /// Number of solutions so we can do something at solution
    int numberSolutions_;
    /// How often to do (code can change)
    int howOften_;
    /// Number of successes
    int numberSuccesses_;
    /// Number of tries
    int numberTries_;
    /// Maximum number of solutions to keep
    int maximumKeepSolutions_;
    /// Number of solutions kept
    int numberKeptSolutions_;
    /// Number of integer variables
    int numberIntegers_;
    /// Local parameter
    int localSpace_;
    /// Values of integer variables
    int ** values_;
};

/** LocalSearch class
 */

class CbcHeuristicVND : public CbcHeuristic {
public:

    // Default Constructor
    CbcHeuristicVND ();

    /* Constructor with model - assumed before cuts
       Initial version does not do Lps
    */
    CbcHeuristicVND (CbcModel & model);

    // Copy constructor
    CbcHeuristicVND ( const CbcHeuristicVND &);

    // Destructor
    ~CbcHeuristicVND ();

    /// Clone
    virtual CbcHeuristic * clone() const;


    /// Assignment operator
    CbcHeuristicVND & operator=(const CbcHeuristicVND& rhs);

    /// Create C++ lines to get to current state
    virtual void generateCpp( FILE * fp) ;

    /// Resets stuff if model changes
    virtual void resetModel(CbcModel * model);

    /// update model (This is needed if cliques update matrix etc)
    virtual void setModel(CbcModel * model);

    using CbcHeuristic::solution ;
    /** returns 0 if no solution, 1 if valid solution.
        Sets solution values if good, sets objective value (only if good)
        This does Relaxation Induced Neighborhood Search
    */
    virtual int solution(double & objectiveValue,
                         double * newSolution);
    /// This version fixes stuff and does IP
    int solutionFix(double & objectiveValue,
                    double * newSolution,
                    const int * keep);

    /// Sets how often to do it
    inline void setHowOften(int value) {
        howOften_ = value;
    }
    /// base solution array so we can set
    inline double * baseSolution() const {
        return baseSolution_;
    }

protected:
    // Data

    /// Number of solutions so we can do something at solution
    int numberSolutions_;
    /// How often to do (code can change)
    int howOften_;
    /// Number of successes
    int numberSuccesses_;
    /// Number of tries
    int numberTries_;
    /// Node when last done
    int lastNode_;
    /// Step size for decomposition
    int stepSize_;
    int k_;
    int kmax_;
    int nDifferent_;
    /// Base solution
    double * baseSolution_;
};

#endif
