/* $Id$ */
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcBranchDynamic_H
#define CbcBranchDynamic_H

#include "CbcBranchActual.hpp"
#include "CoinPackedMatrix.hpp"
#include "CbcSimpleIntegerDynamicPseudoCost.hpp"
#include "CbcDynamicPseudoCostBranchingObject.hpp"

/** Branching decision dynamic class

  This class implements a simple algorithm
  (betterBranch()) for choosing a branching variable when dynamic pseudo costs.
*/

class CbcBranchDynamicDecision : public CbcBranchDecision {
public:
    // Default Constructor
    CbcBranchDynamicDecision ();

    // Copy constructor
    CbcBranchDynamicDecision ( const CbcBranchDynamicDecision &);

    virtual ~CbcBranchDynamicDecision();

    /// Clone
    virtual CbcBranchDecision * clone() const;

    /// Initialize, <i>e.g.</i> before the start of branch selection at a node
    virtual void initialize(CbcModel * model);

    /** \brief Compare two branching objects. Return nonzero if \p thisOne is
           better than \p bestSoFar.

      The routine compares branches using the values supplied in \p numInfUp and
      \p numInfDn until a solution is found by search, after which it uses the
      values supplied in \p changeUp and \p changeDn. The best branching object
      seen so far and the associated parameter values are remembered in the
      \c CbcBranchDynamicDecision object. The nonzero return value is +1 if the
      up branch is preferred, -1 if the down branch is preferred.

      As the names imply, the assumption is that the values supplied for
      \p numInfUp and \p numInfDn will be the number of infeasibilities reported
      by the branching object, and \p changeUp and \p changeDn will be the
      estimated change in objective. Other measures can be used if desired.

      Because an \c CbcBranchDynamicDecision object remembers the current best
      branching candidate (#bestObject_) as well as the values used in the
      comparison, the parameter \p bestSoFar is redundant, hence unused.
    */
    virtual int betterBranch(CbcBranchingObject * thisOne,
                             CbcBranchingObject * bestSoFar,
                             double changeUp, int numInfUp,
                             double changeDn, int numInfDn);
    /** Sets or gets best criterion so far */
    virtual void setBestCriterion(double value);
    virtual double getBestCriterion() const;
    /** Says whether this method can handle both methods -
        1 better, 2 best, 3 both */
    virtual int whichMethod() {
        return 3;
    }

    /** Saves a clone of current branching object.  Can be used to update
        information on object causing branch - after branch */
    virtual void saveBranchingObject(OsiBranchingObject * object) ;
    /** Pass in information on branch just done.
        assumes object can get information from solver */
    virtual void updateInformation(OsiSolverInterface * solver,
                                   const CbcNode * node);


private:

    /// Illegal Assignment operator
    CbcBranchDynamicDecision & operator=(const CbcBranchDynamicDecision& rhs);

    /// data

    /// "best" so far
    double bestCriterion_;

    /// Change up for best
    double bestChangeUp_;

    /// Number of infeasibilities for up
    int bestNumberUp_;

    /// Change down for best
    double bestChangeDown_;

    /// Number of infeasibilities for down
    int bestNumberDown_;

    /// Pointer to best branching object
    CbcBranchingObject * bestObject_;
};
#endif
