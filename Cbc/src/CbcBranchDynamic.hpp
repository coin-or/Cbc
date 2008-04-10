// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcBranchDynamic_H
#define CbcBranchDynamic_H

#include "CbcBranchActual.hpp"
#include "CoinPackedMatrix.hpp"


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
  CbcSimpleIntegerDynamicPseudoCost (CbcModel * model,  int iColumn, double breakEven=0.5);
  
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
  ~CbcSimpleIntegerDynamicPseudoCost ();
  
  using CbcObject::infeasibility ;
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const;

  using CbcObject::createBranch ;
  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way) ;

  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiSolverInterface * solver, 
			       const OsiBranchingInformation * info, int & preferredWay) const;


  /** Create a branching object and indicate which way to branch first.
      
      The branching object has to know how to create branches (fix
      variables, etc.)
  */
  virtual CbcBranchingObject * createBranch(OsiSolverInterface * solver,
					    const OsiBranchingInformation * info, int way) ;
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

  using CbcSimpleInteger::solverBranch ;
  /** Create an OsiSolverBranch object

      This returns NULL if branch not represented by bound changes
  */
  virtual OsiSolverBranch * solverBranch() const;
  
  /// Down pseudo cost
  inline double downDynamicPseudoCost() const
  { return downDynamicPseudoCost_;}
  /// Set down pseudo cost
  inline void setDownDynamicPseudoCost(double value)
  { downDynamicPseudoCost_=value;}

  /// Up pseudo cost
  inline double upDynamicPseudoCost() const
  { return upDynamicPseudoCost_;}
  /// Set up pseudo cost
  inline void setUpDynamicPseudoCost(double value)
  { upDynamicPseudoCost_=value;}

  /// Up down separator
  inline double upDownSeparator() const
  { return upDownSeparator_;}
  /// Set up down separator
  inline void setUpDownSeparator(double value)
  { upDownSeparator_=value;}

  /// Down sum cost
  inline double sumDownCost() const
  { return sumDownCost_;}
  /// Set down sum cost
  inline void setSumDownCost(double value)
  { sumDownCost_=value;}
  /// Add to down sum cost and set last and square
  inline void addToSumDownCost(double value)
  { sumDownCost_+=value;lastDownCost_=value;sumDownCostSquared_ += value*value;}

  /// Up sum cost
  inline double sumUpCost() const
  { return sumUpCost_;}
  /// Set up sum cost
  inline void setSumUpCost(double value)
  { sumUpCost_=value;}
  /// Add to up sum cost and set last and square
  inline void addToSumUpCost(double value)
  { sumUpCost_+=value;lastUpCost_=value;sumUpCostSquared_ += value*value;}

  /// Down sum change
  inline double sumDownChange() const
  { return sumDownChange_;}
  /// Set down sum change
  inline void setSumDownChange(double value)
  { sumDownChange_=value;}
  /// Add to down sum change
  inline void addToSumDownChange(double value)
  { sumDownChange_+=value;}

  /// Up sum change
  inline double sumUpChange() const
  { return sumUpChange_;}
  /// Set up sum change
  inline void setSumUpChange(double value)
  { sumUpChange_=value;}
  /// Add to up sum change and set last and square
  inline void addToSumUpChange(double value)
  { sumUpChange_+=value;}

  /// Sum down decrease number infeasibilities from strong or actual
  inline double sumDownDecrease() const
  { return sumDownDecrease_;}
  /// Set sum down decrease number infeasibilities from strong or actual
  inline void setSumDownDecrease(double value)
  { sumDownDecrease_=value;}
  /// Add to sum down decrease number infeasibilities from strong or actual
  inline void addToSumDownDecrease(double value)
  { sumDownDecrease_+=value;/*lastDownDecrease_ = (int) value;*/}

  /// Sum up decrease number infeasibilities from strong or actual
  inline double sumUpDecrease() const
  { return sumUpDecrease_;}
  /// Set sum up decrease number infeasibilities from strong or actual
  inline void setSumUpDecrease(double value)
  { sumUpDecrease_=value;}
  /// Add to sum up decrease number infeasibilities from strong or actual
  inline void addToSumUpDecrease(double value)
  { sumUpDecrease_+=value;/*lastUpDecrease_ = (int) value;*/}

  /// Down number times
  inline int numberTimesDown() const
  { return numberTimesDown_;}
  /// Set down number times
  inline void setNumberTimesDown(int value)
  { numberTimesDown_=value;}
  /// Increment down number times
  inline void incrementNumberTimesDown()
  { numberTimesDown_++;}

  /// Up number times
  inline int numberTimesUp() const
  { return numberTimesUp_;}
  /// Set up number times
  inline void setNumberTimesUp(int value)
  { numberTimesUp_=value;}
  /// Increment up number times
  inline void incrementNumberTimesUp()
  { numberTimesUp_++;}

  /// Down number times infeasible
  inline int numberTimesDownInfeasible() const
  { return numberTimesDownInfeasible_;}
  /// Set down number times infeasible
  inline void setNumberTimesDownInfeasible(int value)
  { numberTimesDownInfeasible_=value;}
  /// Increment down number times infeasible
  inline void incrementNumberTimesDownInfeasible()
  { numberTimesDownInfeasible_++;}

  /// Up number times infeasible
  inline int numberTimesUpInfeasible() const
  { return numberTimesUpInfeasible_;}
  /// Set up number times infeasible
  inline void setNumberTimesUpInfeasible(int value)
  { numberTimesUpInfeasible_=value;}
  /// Increment up number times infeasible
  inline void incrementNumberTimesUpInfeasible()
  { numberTimesUpInfeasible_++;}

  /// Number of times before trusted
  inline int numberBeforeTrust() const
  { return numberBeforeTrust_;}
  /// Set number of times before trusted
  inline void setNumberBeforeTrust(int value)
  { numberBeforeTrust_=value;}

  /// Return "up" estimate
  virtual double upEstimate() const;
  /// Return "down" estimate (default 1.0e-5)
  virtual double downEstimate() const;
  
  /// method - see below for details
  inline int method() const
  { return method_;}
  /// Set method
  inline void setMethod(int value)
  { method_=value;}

  /// Pass in information on a down branch
  void setDownInformation(double changeObjectiveDown, int changeInfeasibilityDown);
  /// Pass in information on a up branch
  void setUpInformation(double changeObjectiveUp, int changeInfeasibilityUp);
  /// Pass in probing information
  void setProbingInformation(int fixedDown, int fixedUp);

  /// Print - 0 -summary, 1 just before strong
  void print(int type=0, double value=0.0) const;
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
  /// Sum down cost from strong or actual squared
  mutable double sumDownCostSquared_;
  /// Sum up cost from strong or actual squared
  mutable double sumUpCostSquared_;
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
#define CBC_INSTRUMENT
#ifdef CBC_INSTRUMENT
  mutable int numberTimesInfeasible_;
#endif
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
  inline double changeInGuessed() const
  { return changeInGuessed_;}
  /// Set change in guessed
  inline void setChangeInGuessed(double value)
  { changeInGuessed_=value;}
  /// Return object
  inline CbcSimpleIntegerDynamicPseudoCost * object() const
  { return object_;}
  /// Set object
  inline void setObject(CbcSimpleIntegerDynamicPseudoCost * object)
  { object_=object;}

  /** Return the type (an integer identifier) of \c this */
  virtual int type() const { return 400; }

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
  virtual int whichMethod() {return 3;}

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
