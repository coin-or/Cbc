/* $Id$ */
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcBranchCut_H
#define CbcBranchCut_H

#include "CbcBranchBase.hpp"
#include "OsiRowCut.hpp"
#include "CoinPackedMatrix.hpp"
#include "CbcCutBranchingObject.hpp"

/** Define a cut branching class.
    At present empty - all stuff in descendants
*/

class CbcBranchCut : public CbcObject {

public:

    // Default Constructor
    CbcBranchCut ();

    /** In to maintain normal methods
    */
    CbcBranchCut (CbcModel * model);
    // Copy constructor
    CbcBranchCut ( const CbcBranchCut &);

    /// Clone
    virtual CbcObject * clone() const;

    // Assignment operator
    CbcBranchCut & operator=( const CbcBranchCut& rhs);

    // Destructor
    ~CbcBranchCut ();

    /// Infeasibility
    virtual double infeasibility(const OsiBranchingInformation * info,
                                 int &preferredWay) const;

    using CbcObject::feasibleRegion ;
    /** Set bounds to contain the current solution.

      More precisely, for the variable associated with this object, take the
      value given in the current solution, force it within the current bounds
      if required, then set the bounds to fix the variable at the integer
      nearest the solution value.

      At present this will do nothing
    */
    virtual void feasibleRegion();

    /** \brief Return true if branch created by object should fix variables
    */
    virtual bool boundBranch() const ;

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

      At present this does nothing
    */
    virtual CbcBranchingObject * preferredNewFeasible() const;

    /** \brief Given a valid solution (with reduced costs, etc.),
        return a branching object which would give a new feasible
        point in a bad direction.

      As for preferredNewFeasible(), but the preferred branching object will
      force movement in a direction that degrades the objective.

      At present this does nothing
    */
    virtual CbcBranchingObject * notPreferredNewFeasible() const ;

    using CbcObject::resetBounds ;
    /** Reset original upper and lower bound values from the solver.

      Handy for updating bounds held in this object after bounds held in the
      solver have been tightened.
     */
    virtual void resetBounds();


protected:
    /// data

};
#endif
