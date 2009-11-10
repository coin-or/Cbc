// Edwin 11/10/2009-- carved out of CbcBranchActual
#ifndef CbcCbcIntegerPseudoCostBranchingObject_H
#define CbcCbcIntegerPseudoCostBranchingObject_H

#include "CbcIntegerBranchingObject.hpp"
/** Simple branching object for an integer variable with pseudo costs

  This object can specify a two-way branch on an integer variable. For each
  arm of the branch, the upper and lower bounds on the variable can be
  independently specified.

  Variable_ holds the index of the integer variable in the integerVariable_
  array of the model.
*/

class CbcIntegerPseudoCostBranchingObject : public CbcIntegerBranchingObject {

public:

    /// Default constructor
    CbcIntegerPseudoCostBranchingObject ();

    /** Create a standard floor/ceiling branch object

      Specifies a simple two-way branch. Let \p value = x*. One arm of the
      branch will be is lb <= x <= floor(x*), the other ceil(x*) <= x <= ub.
      Specify way = -1 to set the object state to perform the down arm first,
      way = 1 for the up arm.
    */
    CbcIntegerPseudoCostBranchingObject (CbcModel *model, int variable,
                                         int way , double value) ;

    /** Create a degenerate branch object

      Specifies a `one-way branch'. Calling branch() for this object will
      always result in lowerValue <= x <= upperValue. Used to fix a variable
      when lowerValue = upperValue.
    */

    CbcIntegerPseudoCostBranchingObject (CbcModel *model, int variable, int way,
                                         double lowerValue, double upperValue) ;

    /// Copy constructor
    CbcIntegerPseudoCostBranchingObject ( const CbcIntegerPseudoCostBranchingObject &);

    /// Assignment operator
    CbcIntegerPseudoCostBranchingObject & operator= (const CbcIntegerPseudoCostBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    /// Destructor
    virtual ~CbcIntegerPseudoCostBranchingObject ();

    using CbcBranchingObject::branch ;
    /** \brief Sets the bounds for the variable according to the current arm
           of the branch and advances the object state to the next arm.
           This version also changes guessed objective value
    */
    virtual double branch();

    /// Change in guessed
    inline double changeInGuessed() const {
        return changeInGuessed_;
    }
    /// Set change in guessed
    inline void setChangeInGuessed(double value) {
        changeInGuessed_ = value;
    }

    /** Return the type (an integer identifier) of \c this */
    virtual int type() const {
        return 101;
    }

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
    /// Change in guessed objective value for next branch
    double changeInGuessed_;
};

#endif