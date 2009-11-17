//Edwin 11/17/2009 - carved out of CbcBranchLotsize
#ifndef CbcLotsize_H
#define CbcLotsize_H

#include "CbcBranchBase.hpp"
/** Lotsize class */


class CbcLotsize : public CbcObject {

public:

    // Default Constructor
    CbcLotsize ();

    /* Useful constructor - passed model index.
       Also passed valid values - if range then pairs
    */
    CbcLotsize (CbcModel * model, int iColumn,
                int numberPoints, const double * points, bool range = false);

    // Copy constructor
    CbcLotsize ( const CbcLotsize &);

    /// Clone
    virtual CbcObject * clone() const;

    // Assignment operator
    CbcLotsize & operator=( const CbcLotsize& rhs);

    // Destructor
    ~CbcLotsize ();

    /// Infeasibility - large is 0.5
    virtual double infeasibility(const OsiBranchingInformation * info,
                                 int &preferredWay) const;

    using CbcObject::feasibleRegion ;
    /** Set bounds to contain the current solution.

      More precisely, for the variable associated with this object, take the
      value given in the current solution, force it within the current bounds
      if required, then set the bounds to fix the variable at the integer
      nearest the solution value.
    */
    virtual void feasibleRegion();

    /// Creates a branching object
    virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) ;

    /** \brief Given a valid solution (with reduced costs, etc.),
        return a branching object which would give a new feasible
        point in the good direction.

      The preferred branching object will force the variable to be +/-1 from
      its current value, depending on the reduced cost and objective sense.  If
      movement in the direction which improves the objective is impossible due
      to bounds on the variable, the branching object will move in the other
      direction.  If no movement is possible, the method returns NULL.

      Only the bounds on this variable are considered when determining if the new
      point is feasible.
    */
    virtual CbcBranchingObject * preferredNewFeasible() const;

    /** \brief Given a valid solution (with reduced costs, etc.),
        return a branching object which would give a new feasible
        point in a bad direction.

      As for preferredNewFeasible(), but the preferred branching object will
      force movement in a direction that degrades the objective.
    */
    virtual CbcBranchingObject * notPreferredNewFeasible() const ;

    /** Reset original upper and lower bound values from the solver.

      Handy for updating bounds held in this object after bounds held in the
      solver have been tightened.
     */
    virtual void resetBounds(const OsiSolverInterface * solver);

    /** Finds range of interest so value is feasible in range range_ or infeasible
        between hi[range_] and lo[range_+1].  Returns true if feasible.
    */
    bool findRange(double value) const;

    /** Returns floor and ceiling
    */
    virtual void floorCeiling(double & floorLotsize, double & ceilingLotsize, double value,
                              double tolerance) const;

    /// Model column number
    inline int modelSequence() const {
        return columnNumber_;
    }
    /// Set model column number
    inline void setModelSequence(int value) {
        columnNumber_ = value;
    }

    /** Column number if single column object -1 otherwise,
        so returns >= 0
        Used by heuristics
    */
    virtual int columnNumber() const;
    /// Original bounds
    inline double originalLowerBound() const {
        return bound_[0];
    }
    inline double originalUpperBound() const {
        return bound_[rangeType_*numberRanges_-1];
    }
    /// Type - 1 points, 2 ranges
    inline int rangeType() const {
        return rangeType_;
    }
    /// Number of points
    inline int numberRanges() const {
        return numberRanges_;
    }
    /// Ranges
    inline double * bound() const {
        return bound_;
    }
    /** \brief Return true if object can take part in normal heuristics
    */
    virtual bool canDoHeuristics() const {
        return false;
    }

private:
    /// Just for debug (CBC_PRINT defined in CbcBranchLotsize.cpp)
    void printLotsize(double value, bool condition, int type) const;

private:
    /// data

    /// Column number in model
    int columnNumber_;
    /// Type - 1 points, 2 ranges
    int rangeType_;
    /// Number of points
    int numberRanges_;
    // largest gap
    double largestGap_;
    /// Ranges
    double * bound_;
    /// Current range
    mutable int range_;
};

#endif