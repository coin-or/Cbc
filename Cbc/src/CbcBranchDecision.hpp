// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/12/2009 carved from CbcBranchBase

#ifndef CbcBranchDecision_H
#define CbcBranchDecision_H

#include "CbcBranchBase.hpp"

/** Abstract branching decision base class

  In the abstract, an CbcBranchDecision object is expected to be able to
  compare two possible branching choices.

  The #betterBranch() method is the crucial routine. It is expected to be able
  to compare two \link CbcBranchingObject CbcBranchingObjects \endlink.

  See CbcObject for an overview of the three classes (CbcObject,
  CbcBranchingObject, and CbcBranchDecision) which make up cbc's branching
  model.
*/
class CbcModel;
class OsiChooseVariable;

class CbcBranchDecision {
public:
    /// Default Constructor
    CbcBranchDecision ();

    // Copy constructor
    CbcBranchDecision ( const CbcBranchDecision &);

    /// Destructor
    virtual ~CbcBranchDecision();

/// Clone
    virtual CbcBranchDecision * clone() const = 0;

    /// Initialize <i>e.g.</i> before starting to choose a branch at a node
    virtual void initialize(CbcModel * model) = 0;

    /** \brief Compare two branching objects. Return nonzero if branching
           using \p thisOne is better than branching using \p bestSoFar.

      If \p bestSoFar is NULL, the routine should return a nonzero value.
      This routine is used only after strong branching.
      Either this or bestBranch is used depending which user wants.

    */

    virtual int
    betterBranch (CbcBranchingObject * thisOne,
                  CbcBranchingObject * bestSoFar,
                  double changeUp, int numberInfeasibilitiesUp,
                  double changeDown, int numberInfeasibilitiesDown) = 0 ;

    /** \brief Compare N branching objects. Return index of best
        and sets way of branching in chosen object.

      Either this or betterBranch is used depending which user wants.
    */

    virtual int
    bestBranch (CbcBranchingObject ** objects, int numberObjects, int numberUnsatisfied,
                double * changeUp, int * numberInfeasibilitiesUp,
                double * changeDown, int * numberInfeasibilitiesDown,
                double objectiveValue) ;

    /** Says whether this method can handle both methods -
        1 better, 2 best, 3 both */
    virtual int whichMethod() {
        return 2;
    }

    /** Saves a clone of current branching object.  Can be used to update
        information on object causing branch - after branch */
    virtual void saveBranchingObject(OsiBranchingObject * ) {}
    /** Pass in information on branch just done.
        assumes object can get information from solver */
    virtual void updateInformation(OsiSolverInterface * ,
                                   const CbcNode * ) {}
    /** Sets or gets best criterion so far */
    virtual void setBestCriterion(double ) {}
    virtual double getBestCriterion() const {
        return 0.0;
    }
    /// Create C++ lines to get to current state
    virtual void generateCpp( FILE * ) {}
    /// Model
    inline CbcModel * cbcModel() const {
        return model_;
    }
    /* If chooseMethod_ id non-null then the rest is fairly pointless
       as choosemethod_ will be doing all work
     This comment makes more sense if you realise that there's a conversion in
     process from the Cbc branching classes to Osi branching classes. The test
     for use of the Osi branching classes is CbcModel::branchingMethod_
     non-null (i.e., it points to one of these CbcBranchDecision objects) and
     that branch decision object has an OsiChooseVariable method set. In which
     case, we'll use it, rather than the choose[*]Variable methods defined in
     CbcNode.
	*/

    OsiChooseVariable * chooseMethod() const {
        return chooseMethod_;
    }
    /// Set (clone) chooseMethod
    void setChooseMethod(const OsiChooseVariable & method);

protected:

    // Clone of branching object
    CbcBranchingObject * object_;
    /// Pointer to model
    CbcModel * model_;
    /* If chooseMethod_ id non-null then the rest is fairly pointless
       as choosemethod_ will be doing all work
    */
    OsiChooseVariable * chooseMethod_;
private:
    /// Assignment is illegal
    CbcBranchDecision & operator=(const CbcBranchDecision& rhs);

};
#endif

