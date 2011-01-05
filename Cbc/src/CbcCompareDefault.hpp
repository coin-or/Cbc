// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

//Edwin 11/25/09 carved out of CbcCompareActual

#ifndef CbcCompareDefault_H
#define CbcCompareDefault_H


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

/* This is an example of a more complex rule with data
   It is default after first solution
   If weight is 0.0 then it is computed to hit first solution
   less 5%
*/
class CbcCompareDefault  : public CbcCompareBase {
public:
    /// Default Constructor
    CbcCompareDefault () ;
    /// Constructor with weight
    CbcCompareDefault (double weight);

    /// Copy constructor
    CbcCompareDefault ( const CbcCompareDefault &rhs);

    /// Assignment operator
    CbcCompareDefault & operator=( const CbcCompareDefault& rhs);

    /// Clone
    virtual CbcCompareBase * clone() const;
    /// Create C++ lines to get to current state
    virtual void generateCpp( FILE * fp);

    ~CbcCompareDefault() ;
    /* This returns true if weighted value of node y is less than
       weighted value of node x */
    virtual bool test (CbcNode * x, CbcNode * y) ;

    using CbcCompareBase::newSolution ;
    /// This allows method to change behavior as it is called
    /// after each solution
    virtual bool newSolution(CbcModel * model,
                             double objectiveAtContinuous,
                             int numberInfeasibilitiesAtContinuous) ;
    /// This allows method to change behavior
    /// Return true if want tree re-sorted
    virtual bool every1000Nodes(CbcModel * model, int numberNodes);

    /* if weight == -1.0 then fewest infeasibilities (before solution)
       if -2.0 then do breadth first just for first 1000 nodes
       if -3.0 then depth first before solution
    */
    inline double getWeight() const {
        return weight_;
    }
    inline void setWeight(double weight) {
        weight_ = weight;
    }
    /// Cutoff
    inline double getCutoff() const {
        return cutoff_;
    }
    inline void setCutoff(double cutoff) {
        cutoff_ = cutoff;
    }
    /// Best possible solution
    inline double getBestPossible() const {
        return bestPossible_;
    }
    inline void setBestPossible(double bestPossible) {
        bestPossible_ = bestPossible;
    }
    /// Depth above which want to explore first
    inline void setBreadthDepth(int value) {
        breadthDepth_ = value;
    }
    /// Start dive
    void startDive(CbcModel * model);
    /// Clean up diving (i.e. switch off or prepare)
    void cleanDive();
protected:
    /// Weight for each infeasibility
    double weight_;
    /// Weight for each infeasibility - computed from solution
    double saveWeight_;
    /// Cutoff
    double cutoff_;
    /// Best possible solution
    double bestPossible_;
    /// Number of solutions
    int numberSolutions_;
    /// Tree size (at last check)
    int treeSize_;
    /// Depth above which want to explore first
    int breadthDepth_;
    /// Chosen node from estimated (-1 is off)
    int startNodeNumber_;
    /// Node number when dive started
    int afterNodeNumber_;
    /// Indicates doing setup for diving
    bool setupForDiving_ ;
};

#endif //CbcCompareDefault_H

