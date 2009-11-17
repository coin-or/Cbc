/* $Id$ */
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcBranchLotsize_H
#define CbcBranchLotsize_H

#include "CbcBranchBase.hpp"
#include "CbcLotsize.hpp"

/** Lotsize branching object

  This object can specify a two-way branch on an integer variable. For each
  arm of the branch, the upper and lower bounds on the variable can be
  independently specified.

  Variable_ holds the index of the integer variable in the integerVariable_
  array of the model.
*/

class CbcLotsizeBranchingObject : public CbcBranchingObject {

public:

    /// Default constructor
    CbcLotsizeBranchingObject ();

    /** Create a lotsize floor/ceiling branch object

      Specifies a simple two-way branch. Let \p value = x*. One arm of the
      branch will be is lb <= x <= valid range below(x*), the other valid range above(x*) <= x <= ub.
      Specify way = -1 to set the object state to perform the down arm first,
      way = 1 for the up arm.
    */
    CbcLotsizeBranchingObject (CbcModel *model, int variable,
                               int way , double value, const CbcLotsize * lotsize) ;

    /** Create a degenerate branch object

      Specifies a `one-way branch'. Calling branch() for this object will
      always result in lowerValue <= x <= upperValue. Used to fix in valid range
    */

    CbcLotsizeBranchingObject (CbcModel *model, int variable, int way,
                               double lowerValue, double upperValue) ;

    /// Copy constructor
    CbcLotsizeBranchingObject ( const CbcLotsizeBranchingObject &);

    /// Assignment operator
    CbcLotsizeBranchingObject & operator= (const CbcLotsizeBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    /// Destructor
    virtual ~CbcLotsizeBranchingObject ();

    using CbcBranchingObject::branch ;
    /** \brief Sets the bounds for the variable according to the current arm
           of the branch and advances the object state to the next arm.
    */
    virtual double branch();

#if 0
    // No need to override. Default works fine.
    /** Reset every information so that the branching object appears to point to
        the previous child. This method does not need to modify anything in any
        solver. */
    virtual void previousBranch();
#endif

    using CbcBranchingObject::print ;
    /** \brief Print something about branch - only if log level high
    */
    virtual void print();

    /** Return the type (an integer identifier) of \c this */
    virtual int type() const {
        return 300;
    }

    // LL: compareOriginalObject can be inherited from the CbcBranchingObject
    // since variable_ uniquely defines the lot sizing object.

    /** Compare the \c this with \c brObj. \c this and \c brObj must be os the
        same type and must have the same original object, but they may have
        different feasible regions.
        Return the appropriate CbcRangeCompare value (first argument being the
        sub/superset if that's the case). In case of overlap (and if \c
        replaceIfOverlap is true) replace the current branching object with one
        whose feasible region is the overlap.
     */
    virtual CbcRangeCompare compareBranchingObject
    (const CbcBranchingObject* brObj, const bool replaceIfOverlap = false);

protected:
    /// Lower [0] and upper [1] bounds for the down arm (way_ = -1)
    double down_[2];
    /// Lower [0] and upper [1] bounds for the up arm (way_ = 1)
    double up_[2];
};


#endif
