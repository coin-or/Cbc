// Edwin 11/10/2009-- carved out of CbcBranchActual
#ifndef CbcIntegerBranchingObject_H
#define CbcIntegerBranchingObject_H

/** Simple branching object for an integer variable

  This object can specify a two-way branch on an integer variable. For each
  arm of the branch, the upper and lower bounds on the variable can be
  independently specified.

  Variable_ holds the index of the integer variable in the integerVariable_
  array of the model.
*/

class CbcIntegerBranchingObject : public CbcBranchingObject {

public:

    /// Default constructor
    CbcIntegerBranchingObject ();

    /** Create a standard floor/ceiling branch object

      Specifies a simple two-way branch. Let \p value = x*. One arm of the
      branch will be lb <= x <= floor(x*), the other ceil(x*) <= x <= ub.
      Specify way = -1 to set the object state to perform the down arm first,
      way = 1 for the up arm.
    */
    CbcIntegerBranchingObject (CbcModel *model, int variable,
                               int way , double value) ;

    /** Create a degenerate branch object

      Specifies a `one-way branch'. Calling branch() for this object will
      always result in lowerValue <= x <= upperValue. Used to fix a variable
      when lowerValue = upperValue.
    */

    CbcIntegerBranchingObject (CbcModel *model, int variable, int way,
                               double lowerValue, double upperValue) ;

    /// Copy constructor
    CbcIntegerBranchingObject ( const CbcIntegerBranchingObject &);

    /// Assignment operator
    CbcIntegerBranchingObject & operator= (const CbcIntegerBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    /// Destructor
    virtual ~CbcIntegerBranchingObject ();

    /// Does part of constructor
    void fillPart ( int variable, int way , double value) ;
    using CbcBranchingObject::branch ;
    /** \brief Sets the bounds for the variable according to the current arm
           of the branch and advances the object state to the next arm.
           Returns change in guessed objective on next branch
    */
    virtual double branch();
    /** Update bounds in solver as in 'branch' and update given bounds.
        branchState is -1 for 'down' +1 for 'up' */
    virtual void fix(OsiSolverInterface * solver,
                     double * lower, double * upper,
                     int branchState) const ;

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

    /// Lower and upper bounds for down branch
    inline const double * downBounds() const {
        return down_;
    }
    /// Lower and upper bounds for up branch
    inline const double * upBounds() const {
        return up_;
    }
    /// Set lower and upper bounds for down branch
    inline void setDownBounds(const double bounds[2]) {
        memcpy(down_, bounds, 2*sizeof(double));
    }
    /// Set lower and upper bounds for up branch
    inline void setUpBounds(const double bounds[2]) {
        memcpy(up_, bounds, 2*sizeof(double));
    }
#ifdef FUNNY_BRANCHING
    /** Which variable (top bit if upper bound changing,
        next bit if on down branch */
    inline const int * variables() const {
        return variables_;
    }
    // New bound
    inline const double * newBounds() const {
        return newBounds_;
    }
    /// Number of bound changes
    inline int numberExtraChangedBounds() const {
        return numberExtraChangedBounds_;
    }
    /// Just apply extra bounds to one variable - COIN_DBL_MAX ignore
    int applyExtraBounds(int iColumn, double lower, double upper, int way) ;
    /// Deactivate bounds for branching
    void deactivate();
    /// Are active bounds for branching
    inline bool active() const {
        return (down_[1] != -COIN_DBL_MAX);
    }
#endif

    /** Return the type (an integer identifier) of \c this */
    virtual int type() const {
        return 100;
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
    /// Lower [0] and upper [1] bounds for the down arm (way_ = -1)
    double down_[2];
    /// Lower [0] and upper [1] bounds for the up arm (way_ = 1)
    double up_[2];
#ifdef FUNNY_BRANCHING
    /** Which variable (top bit if upper bound changing)
        next bit if changing on down branch only */
    int * variables_;
    // New bound
    double * newBounds_;
    /// Number of Extra bound changes
    int numberExtraChangedBounds_;
#endif
};
#endif