// Edwin 11/17/2009-- carved out of CbcBranchDynamic
#ifndef CbcDynamicPseudoCostBranchingObject_H
#define CbcDynamicPseudoCostBranchingObject_H

//#include "CbcSimpleIntegerDynamicPseudoCost.hpp"
#include "CbcBranchDynamic.hpp"
/** Simple branching object for an integer variable with pseudo costs

  This object can specify a two-way branch on an integer variable. For each
  arm of the branch, the upper and lower bounds on the variable can be
  independently specified.

  Variable_ holds the index of the integer variable in the integerVariable_
  array of the model.
*/

class CbcDynamicPseudoCostBranchingObject : public CbcIntegerBranchingObject {

public:

    /// Default constructor
    CbcDynamicPseudoCostBranchingObject ();

    /** Create a standard floor/ceiling branch object

      Specifies a simple two-way branch. Let \p value = x*. One arm of the
      branch will be is lb <= x <= floor(x*), the other ceil(x*) <= x <= ub.
      Specify way = -1 to set the object state to perform the down arm first,
      way = 1 for the up arm.
    */
    CbcDynamicPseudoCostBranchingObject (CbcModel *model, int variable,
                                         int way , double value,
                                         CbcSimpleIntegerDynamicPseudoCost * object) ;

    /** Create a degenerate branch object

      Specifies a `one-way branch'. Calling branch() for this object will
      always result in lowerValue <= x <= upperValue. Used to fix a variable
      when lowerValue = upperValue.
    */

    CbcDynamicPseudoCostBranchingObject (CbcModel *model, int variable, int way,
                                         double lowerValue, double upperValue) ;

    /// Copy constructor
    CbcDynamicPseudoCostBranchingObject ( const CbcDynamicPseudoCostBranchingObject &);

    /// Assignment operator
    CbcDynamicPseudoCostBranchingObject & operator= (const CbcDynamicPseudoCostBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    /// Destructor
    virtual ~CbcDynamicPseudoCostBranchingObject ();

    /// Does part of constructor
    void fillPart (int variable,
                   int way , double value,
                   CbcSimpleIntegerDynamicPseudoCost * object) ;

    using CbcBranchingObject::branch ;
    /** \brief Sets the bounds for the variable according to the current arm
           of the branch and advances the object state to the next arm.
           This version also changes guessed objective value
    */
    virtual double branch();

    /** Some branchingObjects may claim to be able to skip
        strong branching.  If so they have to fill in CbcStrongInfo.
        The object mention in incoming CbcStrongInfo must match.
        Returns nonzero if skip is wanted */
    virtual int fillStrongInfo( CbcStrongInfo & info);

    /// Change in guessed
    inline double changeInGuessed() const {
        return changeInGuessed_;
    }
    /// Set change in guessed
    inline void setChangeInGuessed(double value) {
        changeInGuessed_ = value;
    }
    /// Return object
    inline CbcSimpleIntegerDynamicPseudoCost * object() const {
        return object_;
    }
    /// Set object
    inline void setObject(CbcSimpleIntegerDynamicPseudoCost * object) {
        object_ = object;
    }

    /** Return the type (an integer identifier) of \c this */
    virtual int type() const {
        return 400;
    }

    // LL: compareOriginalObject and compareBranchingObject are inherited from
    // CbcIntegerBranchingObject thus need not be declared/defined here. After
    // all, this kind of branching object is simply using pseudocosts to make
    // decisions, but once the decisions are made they are the same kind as in
    // the underlying class.

protected:
    /// Change in guessed objective value for next branch
    double changeInGuessed_;
    /// Pointer back to object
    CbcSimpleIntegerDynamicPseudoCost * object_;

};

#endif