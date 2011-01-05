// $Id$
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/17/2009 - carved out of CbcBranchDynamic

#ifndef CbcSimpleIntegerDynamicPseudoCost_H
#define CbcSimpleIntegerDynamicPseudoCost_H

#include "CbcSimpleInteger.hpp"

#define TYPERATIO 0.9
#define MINIMUM_MOVEMENT 0.1
#define TYPE2 0
// was 1 - but that looks flakey
#define INFEAS 1
#define MOD_SHADOW 1
// weight at 1.0 is max min
#define WEIGHT_AFTER 0.8
#define WEIGHT_BEFORE 0.1
//Stolen from Constraint Integer Programming book (with epsilon change)
#define WEIGHT_PRODUCT


/** Define a single integer class but with dynamic pseudo costs.
    Based on work by Achterberg, Koch and Martin.

    It is wild overkill but to keep design all twiddly things are in each.
    This could be used for fine tuning.

 */


class CbcSimpleIntegerDynamicPseudoCost : public CbcSimpleInteger {

public:

    // Default Constructor
    CbcSimpleIntegerDynamicPseudoCost ();

    // Useful constructor - passed  model index
    CbcSimpleIntegerDynamicPseudoCost (CbcModel * model,  int iColumn, double breakEven = 0.5);

    // Useful constructor - passed  model index and pseudo costs
    CbcSimpleIntegerDynamicPseudoCost (CbcModel * model, int iColumn,
                                       double downDynamicPseudoCost, double upDynamicPseudoCost);

    // Useful constructor - passed  model index and pseudo costs
    CbcSimpleIntegerDynamicPseudoCost (CbcModel * model, int dummy, int iColumn,
                                       double downDynamicPseudoCost, double upDynamicPseudoCost);

    // Copy constructor
    CbcSimpleIntegerDynamicPseudoCost ( const CbcSimpleIntegerDynamicPseudoCost &);

    /// Clone
    virtual CbcObject * clone() const;

    // Assignment operator
    CbcSimpleIntegerDynamicPseudoCost & operator=( const CbcSimpleIntegerDynamicPseudoCost& rhs);

    // Destructor
    virtual ~CbcSimpleIntegerDynamicPseudoCost ();

    /// Infeasibility - large is 0.5
    virtual double infeasibility(const OsiBranchingInformation * info,
                                 int &preferredWay) const;

    /// Creates a branching object
    virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) ;


    /// Fills in a created branching object
    void fillCreateBranch(CbcIntegerBranchingObject * branching, const OsiBranchingInformation * info, int way) ;


    /** Pass in information on branch just done and create CbcObjectUpdateData instance.
        If object does not need data then backward pointer will be NULL.
        Assumes can get information from solver */
    virtual CbcObjectUpdateData createUpdateInformation(const OsiSolverInterface * solver,
            const CbcNode * node,
            const CbcBranchingObject * branchingObject);
    /// Update object by CbcObjectUpdateData
    virtual void updateInformation(const CbcObjectUpdateData & data) ;
    /// Copy some information i.e. just variable stuff
    void copySome(const CbcSimpleIntegerDynamicPseudoCost * otherObject);
    /// Updates stuff like pseudocosts before threads
    virtual void updateBefore(const OsiObject * rhs) ;
    /// Updates stuff like pseudocosts after threads finished
    virtual void updateAfter(const OsiObject * rhs, const OsiObject * baseObject) ;
    /// Updates stuff like pseudocosts after mini branch and bound
    void updateAfterMini(int numberDown, int numberDownInfeasible, double sumDown,
                         int numberUp, int numberUpInfeasible, double sumUp);

    using CbcSimpleInteger::solverBranch ;
    /** Create an OsiSolverBranch object

        This returns NULL if branch not represented by bound changes
    */
    virtual OsiSolverBranch * solverBranch() const;

    /// Down pseudo cost
    inline double downDynamicPseudoCost() const {
        return downDynamicPseudoCost_;
    }
    /// Set down pseudo cost
    void setDownDynamicPseudoCost(double value) ;
    /// Modify down pseudo cost in a slightly different way
    void updateDownDynamicPseudoCost(double value);

    /// Up pseudo cost
    inline double upDynamicPseudoCost() const {
        return upDynamicPseudoCost_;
    }
    /// Set up pseudo cost
    void setUpDynamicPseudoCost(double value);
    /// Modify up pseudo cost in a slightly different way
    void updateUpDynamicPseudoCost(double value);

    /// Down pseudo shadow price cost
    inline double downShadowPrice() const {
        return downShadowPrice_;
    }
    /// Set down pseudo shadow price cost
    inline void setDownShadowPrice(double value) {
        downShadowPrice_ = value;
    }
    /// Up pseudo shadow price cost
    inline double upShadowPrice() const {
        return upShadowPrice_;
    }
    /// Set up pseudo shadow price cost
    inline void setUpShadowPrice(double value) {
        upShadowPrice_ = value;
    }

    /// Up down separator
    inline double upDownSeparator() const {
        return upDownSeparator_;
    }
    /// Set up down separator
    inline void setUpDownSeparator(double value) {
        upDownSeparator_ = value;
    }

    /// Down sum cost
    inline double sumDownCost() const {
        return sumDownCost_;
    }
    /// Set down sum cost
    inline void setSumDownCost(double value) {
        sumDownCost_ = value;
    }
    /// Add to down sum cost and set last and square
    inline void addToSumDownCost(double value) {
        sumDownCost_ += value;
        lastDownCost_ = value;
    }

    /// Up sum cost
    inline double sumUpCost() const {
        return sumUpCost_;
    }
    /// Set up sum cost
    inline void setSumUpCost(double value) {
        sumUpCost_ = value;
    }
    /// Add to up sum cost and set last and square
    inline void addToSumUpCost(double value) {
        sumUpCost_ += value;
        lastUpCost_ = value;
    }

    /// Down sum change
    inline double sumDownChange() const {
        return sumDownChange_;
    }
    /// Set down sum change
    inline void setSumDownChange(double value) {
        sumDownChange_ = value;
    }
    /// Add to down sum change
    inline void addToSumDownChange(double value) {
        sumDownChange_ += value;
    }

    /// Up sum change
    inline double sumUpChange() const {
        return sumUpChange_;
    }
    /// Set up sum change
    inline void setSumUpChange(double value) {
        sumUpChange_ = value;
    }
    /// Add to up sum change and set last and square
    inline void addToSumUpChange(double value) {
        sumUpChange_ += value;
    }

    /// Sum down decrease number infeasibilities from strong or actual
    inline double sumDownDecrease() const {
        return sumDownDecrease_;
    }
    /// Set sum down decrease number infeasibilities from strong or actual
    inline void setSumDownDecrease(double value) {
        sumDownDecrease_ = value;
    }
    /// Add to sum down decrease number infeasibilities from strong or actual
    inline void addToSumDownDecrease(double value) {
        sumDownDecrease_ += value;/*lastDownDecrease_ = (int) value;*/
    }

    /// Sum up decrease number infeasibilities from strong or actual
    inline double sumUpDecrease() const {
        return sumUpDecrease_;
    }
    /// Set sum up decrease number infeasibilities from strong or actual
    inline void setSumUpDecrease(double value) {
        sumUpDecrease_ = value;
    }
    /// Add to sum up decrease number infeasibilities from strong or actual
    inline void addToSumUpDecrease(double value) {
        sumUpDecrease_ += value;/*lastUpDecrease_ = (int) value;*/
    }

    /// Down number times
    inline int numberTimesDown() const {
        return numberTimesDown_;
    }
    /// Set down number times
    inline void setNumberTimesDown(int value) {
        numberTimesDown_ = value;
    }
    /// Increment down number times
    inline void incrementNumberTimesDown() {
        numberTimesDown_++;
    }

    /// Up number times
    inline int numberTimesUp() const {
        return numberTimesUp_;
    }
    /// Set up number times
    inline void setNumberTimesUp(int value) {
        numberTimesUp_ = value;
    }
    /// Increment up number times
    inline void incrementNumberTimesUp() {
        numberTimesUp_++;
    }

    /// Down number times infeasible
    inline int numberTimesDownInfeasible() const {
        return numberTimesDownInfeasible_;
    }
    /// Set down number times infeasible
    inline void setNumberTimesDownInfeasible(int value) {
        numberTimesDownInfeasible_ = value;
    }
    /// Increment down number times infeasible
    inline void incrementNumberTimesDownInfeasible() {
        numberTimesDownInfeasible_++;
    }

    /// Up number times infeasible
    inline int numberTimesUpInfeasible() const {
        return numberTimesUpInfeasible_;
    }
    /// Set up number times infeasible
    inline void setNumberTimesUpInfeasible(int value) {
        numberTimesUpInfeasible_ = value;
    }
    /// Increment up number times infeasible
    inline void incrementNumberTimesUpInfeasible() {
        numberTimesUpInfeasible_++;
    }

    /// Number of times before trusted
    inline int numberBeforeTrust() const {
        return numberBeforeTrust_;
    }
    /// Set number of times before trusted
    inline void setNumberBeforeTrust(int value) {
        numberBeforeTrust_ = value;
    }
    /// Increment number of times before trusted
    inline void incrementNumberBeforeTrust() {
        numberBeforeTrust_++;
    }

    /// Return "up" estimate
    virtual double upEstimate() const;
    /// Return "down" estimate (default 1.0e-5)
    virtual double downEstimate() const;

    /// method - see below for details
    inline int method() const {
        return method_;
    }
    /// Set method
    inline void setMethod(int value) {
        method_ = value;
    }

    /// Pass in information on a down branch
    void setDownInformation(double changeObjectiveDown, int changeInfeasibilityDown);
    /// Pass in information on a up branch
    void setUpInformation(double changeObjectiveUp, int changeInfeasibilityUp);
    /// Pass in probing information
    void setProbingInformation(int fixedDown, int fixedUp);

    /// Print - 0 -summary, 1 just before strong
    void print(int type = 0, double value = 0.0) const;
    /// Same - returns true if contents match(ish)
    bool same(const CbcSimpleIntegerDynamicPseudoCost * obj) const;
protected:
    /// data

    /// Down pseudo cost
    double downDynamicPseudoCost_;
    /// Up pseudo cost
    double upDynamicPseudoCost_;
    /** Up/down separator
        If >0.0 then do first branch up if value-floor(value)
        >= this value
    */
    double upDownSeparator_;
    /// Sum down cost from strong or actual
    double sumDownCost_;
    /// Sum up cost from strong or actual
    double sumUpCost_;
    /// Sum of all changes to x when going down
    double sumDownChange_;
    /// Sum of all changes to x when going up
    double sumUpChange_;
    /// Current pseudo-shadow price estimate down
    mutable double downShadowPrice_;
    /// Current pseudo-shadow price estimate up
    mutable double upShadowPrice_;
    /// Sum down decrease number infeasibilities from strong or actual
    double sumDownDecrease_;
    /// Sum up decrease number infeasibilities from strong or actual
    double sumUpDecrease_;
    /// Last down cost from strong (i.e. as computed by last strong)
    double lastDownCost_;
    /// Last up cost from strong (i.e. as computed by last strong)
    double lastUpCost_;
    /// Last down decrease number infeasibilities from strong (i.e. as computed by last strong)
    mutable int lastDownDecrease_;
    /// Last up decrease number infeasibilities from strong (i.e. as computed by last strong)
    mutable int lastUpDecrease_;
    /// Number of times we have gone down
    int numberTimesDown_;
    /// Number of times we have gone up
    int numberTimesUp_;
    /// Number of times we have been infeasible going down
    int numberTimesDownInfeasible_;
    /// Number of times we have been infeasible going up
    int numberTimesUpInfeasible_;
    /// Number of branches before we trust
    int numberBeforeTrust_;
    /// Number of local probing fixings going down
    int numberTimesDownLocalFixed_;
    /// Number of local probing fixings going up
    int numberTimesUpLocalFixed_;
    /// Number of total probing fixings going down
    double numberTimesDownTotalFixed_;
    /// Number of total probing fixings going up
    double numberTimesUpTotalFixed_;
    /// Number of times probing done
    int numberTimesProbingTotal_;
    /// Number of times infeasible when tested
    /** Method -
        0 - pseudo costs
        1 - probing
    */
    int method_;
};
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
    virtual CbcBranchObjType type() const {
        return SimpleIntegerDynamicPseudoCostBranchObj;
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

