// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcBranchLotsizeSimple_H
#define CbcBranchLotsizeSimple_H

#include "CbcBranchBase.hpp"

/** LotsizeSimple class
    This is just the lotsize class restricted to points not ranges
    and with debug etc taken out
 */


class CbcLotsizeSimple : public CbcObject {

public:

  // Default Constructor 
  CbcLotsizeSimple ();

  /* Useful constructor - passed model index.
     Also passed valid values 
  */
  CbcLotsizeSimple (CbcModel * model, int iColumn,
	      int numberPoints, const double * points);
  
  // Copy constructor 
  CbcLotsizeSimple ( const CbcLotsizeSimple & rhs);
   
  /// Clone
  virtual CbcObject * clone() const;

  // Assignment operator 
  CbcLotsizeSimple & operator=( const CbcLotsizeSimple& rhs);

  // Destructor 
  ~CbcLotsizeSimple ();
  
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const;

  /** Set bounds to contain the current solution.

    More precisely, for the variable associated with this object, take the
    value given in the current solution, force it within the current bounds
    if required, then set the bounds to fix the variable at the integer
    nearest the solution value.
  */
  virtual void feasibleRegion();

  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way) ;
  
  /// Model column number
  inline int modelSequence() const
  {return columnNumber_;};

  /** Column number if single column object -1 otherwise,
      so returns >= 0
      Used by heuristics
  */
  virtual int columnNumber() const
  {return columnNumber_;};
  /// Original bounds
  inline double originalLowerBound() const
  { return bound_[0];};
  inline double originalUpperBound() const
  { return bound_[numberRanges_-1];};
  /// Number of points
  inline int numberRanges() const
  { return numberRanges_;};
  /// Ranges
  inline double * bound() const
  { return bound_;};

  /** Finds range of interest so value is feasible in range range_ or infeasible 
      between bound_[range_] and bound_[range_+1].  Returns true if feasible.
  */
  bool findRange(double value) const;
  
  /** Returns floor and ceiling
  */
  virtual void floorCeiling(double & floorLotsizeSimple, double & ceilingLotsizeSimple, double value,
			    double tolerance) const;
  /// Simple function to return a solution value strictly feasible
  double cleanValue() const;
protected:
  /// data

  /// Column number in model
  int columnNumber_;
  /// Number of points
  int numberRanges_;
  // largest gap (used for normalizing infeasibility reporting)
  double largestGap_;
  /// Ranges
  double * bound_;
  /// Current range
  mutable int range_;
};


/** LotsizeSimple branching object

  This object can specify a two-way branch on an integer variable. For each
  arm of the branch, the upper and lower bounds on the variable can be
  independently specified.
  
  Variable_ holds the index of the integer variable in the integerVariable_
  array of the model.
*/

class CbcLotsizeSimpleBranchingObject : public CbcBranchingObject {

public:

  /// Default constructor 
  CbcLotsizeSimpleBranchingObject ();

  /** Create a lotsize floor/ceiling branch object

    Specifies a simple two-way branch. Let \p value = x*. One arm of the
    branch will be is lb <= x <= valid range below(x*), the other valid range above(x*) <= x <= ub.
    Specify way = -1 to set the object state to perform the down arm first,
    way = 1 for the up arm.
  */
  CbcLotsizeSimpleBranchingObject (CbcModel *model, int variable,
			     int way , double value,const CbcLotsizeSimple * lotsize) ;
  
  /// Copy constructor 
  CbcLotsizeSimpleBranchingObject ( const CbcLotsizeSimpleBranchingObject &);
   
  /// Assignment operator 
  CbcLotsizeSimpleBranchingObject & operator= (const CbcLotsizeSimpleBranchingObject& rhs);

  /// Clone
  virtual CbcBranchingObject * clone() const;

  /// Destructor 
  virtual ~CbcLotsizeSimpleBranchingObject ();
  
  /** \brief Sets the bounds for the variable according to the current arm
	     of the branch and advances the object state to the next arm.
  */
  virtual double branch();

protected:
  /// Lower [0] and upper [1] bounds for the down arm (way_ = -1)
  double down_[2];
  /// Lower [0] and upper [1] bounds for the up arm (way_ = 1)
  double up_[2];
};


#endif
