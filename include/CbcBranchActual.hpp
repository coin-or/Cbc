// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcBranchActual_H
#define CbcBranchActual_H

#include "CbcBranchBase.hpp"
#include "CoinPackedMatrix.hpp"

/// Define a clique class


class CbcClique : public CbcObject {

public:

  // Default Constructor 
  CbcClique ();

  /** Useful constructor (which are integer indices)
      slack can denote a slack in set.
      If type == NULL then as if 1
  */
  CbcClique (CbcModel * model, int cliqueType, int numberMembers,
	     const int * which, const char * type,
	     int identifier,int slack=-1);
  
  // Copy constructor 
  CbcClique ( const CbcClique &);
   
  /// Clone
  virtual CbcObject * clone() const;

  // Assignment operator 
  CbcClique & operator=( const CbcClique& rhs);

  // Destructor 
  ~CbcClique ();
  
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const;

  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion();
  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way) const;
  /// Number of members
  inline int numberMembers() const
  {return numberMembers_;};

  /// Number of Non SOS members i.e. fixing to zero is strong
  inline int numberNonSOSMembers() const
  {return numberNonSOSMembers_;};

  /// Members (indices in range 0 ... numberIntegers_-1)
  inline const int * members() const
  {return members_;};

  /** Type of each member i.e. which way is strong 0=non SOS, 1 =SOS,
      index is 0 ... numberMembers_-1 */
  inline const char type(int index) const
  {if (type_) return type_[index]; else return 1;};

  /// Clique type - 0 <=, 1 == 
  inline int cliqueType() const
  {return cliqueType_;};

protected:
  /// data
  /// Number of members
  int numberMembers_;

  /// Number of Non SOS members i.e. fixing to zero is strong
  int numberNonSOSMembers_;

  /// Members (indices in range 0 ... numberIntegers_-1)
  int * members_;

  /// Type of each member 0=SOS, 1 =clique
  char * type_;

  /// Clique type - 0 <=, 1 ==
   int cliqueType_;

  /// Which one is slack (if any) sequence within this set
  int slack_;
};

/** Define Special Ordered Sets of type 1 and 2.  These do not have to be
    integer - so do not appear in lists of integers.
    
    which_ points directly to columns of matrix
*/


class CbcSOS : public CbcObject {

public:

  // Default Constructor 
  CbcSOS ();

  /** Useful constructor - which are integer indices
      and  weights are also given.  If null then 0,1,2..
      type is SOS type
  */
  CbcSOS (CbcModel * model, int numberMembers,
	   const int * which, const double * weights, int identifier,
	  int type=1);
  
  // Copy constructor 
  CbcSOS ( const CbcSOS &);
   
  /// Clone
  virtual CbcObject * clone() const;

  // Assignment operator 
  CbcSOS & operator=( const CbcSOS& rhs);

  // Destructor 
  ~CbcSOS ();
  
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const;

  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion();
  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way) const;

  /// Number of members
  inline int numberMembers() const
  {return numberMembers_;};

  /// Members (indices in range 0 ... numberColumns-1)
  inline const int * members() const
  {return members_;};

  /// SOS type
  inline int sosType() const
  {return sosType_;};

  /** Array of weights */
  inline const double * weights() const
  { return weights_;};

private:
  /// data

  /// Members (indices in range 0 ... numberColumns-1)
  int * members_;
  /// Weights
  double * weights_;

  /// Number of members
  int numberMembers_;
  /// SOS type
   int sosType_;
};

/// Define a single integer class


class CbcSimpleInteger : public CbcObject {

public:

  // Default Constructor 
  CbcSimpleInteger ();

  // Useful constructor - passed integer index and model index
  CbcSimpleInteger (CbcModel * model, int sequence, int iColumn, double breakEven=0.5);
  
  // Copy constructor 
  CbcSimpleInteger ( const CbcSimpleInteger &);
   
  /// Clone
  virtual CbcObject * clone() const;

  // Assignment operator 
  CbcSimpleInteger & operator=( const CbcSimpleInteger& rhs);

  // Destructor 
  ~CbcSimpleInteger ();
  
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
  virtual CbcBranchingObject * createBranch(int way) const;

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
  */
  virtual CbcBranchingObject * preferredNewFeasible() const;
  
  /** \brief Given a valid solution (with reduced costs, etc.),
      return a branching object which would give a new feasible
      point in a bad direction.

    As for preferredNewFeasible(), but the preferred branching object will
    force movement in a direction that degrades the objective.
  */
  virtual CbcBranchingObject * notPreferredNewFeasible() const ;
  
  /** Reset original upper and lower bound values from the solver.
  
    Handy for updating bounds held in this object after bounds held in the
    solver have been tightened.
   */
  virtual void resetBounds();
  
  /// Sequence number
  inline int sequence() const
  {return sequence_;};

  /// Model column number
  inline int modelSequence() const
  {return columnNumber_;};

  /// Original bounds
  inline double originalLowerBound() const
  { return originalLower_;};
  inline double originalUpperBound() const
  { return originalUpper_;};
  /// Breakeven e.g 0.7 -> >= 0.7 go up first
  inline double breakEven() const
  { return breakEven_;};
  /// Set breakeven e.g 0.7 -> >= 0.7 go up first
  inline void setBreakEven(double value)
  { breakEven_=value;};


protected:
  /// data

  /// Sequence
  int sequence_;
  /// Column number in model
  int columnNumber_;
  /// Original lower bound
  double originalLower_;
  /// Original upper bound
  double originalUpper_;
  /// Breakeven i.e. >= this preferred is up
  double breakEven_;
};


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
    branch will be is lb <= x <= floor(x*), the other ceil(x*) <= x <= ub.
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
  
  /** \brief Sets the bounds for the variable according to the current arm
	     of the branch and advances the object state to the next arm.
	     Returns change in guessed objective on next branch
  */
  virtual double branch(bool normalBranch=false);

  /** \brief Print something about branch - only if log level high
  */
  virtual void print(bool normalBranch);

protected:
  /// Lower [0] and upper [1] bounds for the down arm (way_ = -1)
  double down_[2];
  /// Lower [0] and upper [1] bounds for the up arm (way_ = 1)
  double up_[2];
};


/// Define a single integer class but with pseudo costs


class CbcSimpleIntegerPseudoCost : public CbcSimpleInteger {

public:

  // Default Constructor 
  CbcSimpleIntegerPseudoCost ();

  // Useful constructor - passed integer index and model index
  CbcSimpleIntegerPseudoCost (CbcModel * model, int sequence, int iColumn, double breakEven=0.5);
  
  // Useful constructor - passed integer index and model index and pseudo costs
  CbcSimpleIntegerPseudoCost (CbcModel * model, int sequence, int iColumn, 
			      double downPseudoCost, double upPseudoCost);
  
  // Copy constructor 
  CbcSimpleIntegerPseudoCost ( const CbcSimpleIntegerPseudoCost &);
   
  /// Clone
  virtual CbcObject * clone() const;

  // Assignment operator 
  CbcSimpleIntegerPseudoCost & operator=( const CbcSimpleIntegerPseudoCost& rhs);

  // Destructor 
  ~CbcSimpleIntegerPseudoCost ();
  
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const;

  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way) const;

  /// Down pseudo cost
  inline double downPseudoCost() const
  { return downPseudoCost_;};
  /// Set down pseudo cost
  inline void setDownPseudoCost(double value)
  { downPseudoCost_=value;};

  /// Up pseudo cost
  inline double upPseudoCost() const
  { return upPseudoCost_;};
  /// Set up pseudo cost
  inline void setUpPseudoCost(double value)
  { upPseudoCost_=value;};

  /// Return "up" estimate
  virtual double upEstimate() const;
  /// Return "down" estimate (default 1.0e-5)
  virtual double downEstimate() const;
  
  /// method - see below for details
  inline int method() const
  { return method_;};
  /// Set method
  inline void setMethod(int value)
  { method_=value;};

protected:
  /// data

  /// Down pseudo cost
  double downPseudoCost_;
  /// Up pseudo cost
  double upPseudoCost_;
  /** Method - 
      0 - normal - return min (up,down)
      1 - if before any solution return max(up,down)
      2 - if before branched solution return max(up,down)
      3 - always return max(up,down)
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
  
  /** \brief Sets the bounds for the variable according to the current arm
	     of the branch and advances the object state to the next arm.
	     This version also changes guessed objective value
  */
  virtual double branch(bool normalBranch=false);

  /// Change in guessed
  inline double changeInGuessed() const
  { return changeInGuessed_;};
  /// Set change in guessed
  inline void setChangeInGuessed(double value)
  { changeInGuessed_=value;};
protected:
  /// Change in guessed objective value for next branch
  double changeInGuessed_;
};


/** Branching object for unordered cliques

    Intended for cliques which are long enough to make it worthwhile
    but <= 64 members.  There will also be ones for long cliques. 

    Variable_ is the clique id number (redundant, as the object also holds a
    pointer to the clique.
 */
class CbcCliqueBranchingObject : public CbcBranchingObject {

public:

  // Default Constructor 
  CbcCliqueBranchingObject ();

  // Useful constructor
  CbcCliqueBranchingObject (CbcModel * model,  const CbcClique * clique,
			    int way,
			    int numberOnDownSide, const int * down,
			    int numberOnUpSide, const int * up);
  
  // Copy constructor 
  CbcCliqueBranchingObject ( const CbcCliqueBranchingObject &);
   
  // Assignment operator 
  CbcCliqueBranchingObject & operator=( const CbcCliqueBranchingObject& rhs);

  /// Clone
  virtual CbcBranchingObject * clone() const;

  // Destructor 
  virtual ~CbcCliqueBranchingObject ();
  
  /// Does next branch and updates state
  virtual double branch(bool normalBranch=false);

  /** \brief Print something about branch - only if log level high
  */
  virtual void print(bool normalBranch);
private:
  /// data
  const CbcClique * clique_;
  /// downMask - bit set to fix to weak bounds, not set to leave unfixed
  unsigned int downMask_[2];
  /// upMask - bit set to fix to weak bounds, not set to leave unfixed
  unsigned int upMask_[2];
};


/** Unordered Clique Branching Object class.
    These are for cliques which are > 64 members
    Variable is number of clique.
 */
class CbcLongCliqueBranchingObject : public CbcBranchingObject {

public:

  // Default Constructor 
  CbcLongCliqueBranchingObject ();

  // Useful constructor
  CbcLongCliqueBranchingObject (CbcModel * model,  const CbcClique * clique,
				 int way,
			    int numberOnDownSide, const int * down,
			    int numberOnUpSide, const int * up);
  
  // Copy constructor 
  CbcLongCliqueBranchingObject ( const CbcLongCliqueBranchingObject &);
   
  // Assignment operator 
  CbcLongCliqueBranchingObject & operator=( const CbcLongCliqueBranchingObject& rhs);

  /// Clone
  virtual CbcBranchingObject * clone() const;

  // Destructor 
  virtual ~CbcLongCliqueBranchingObject ();
  
  /// Does next branch and updates state
  virtual double branch(bool normalBranch=false);

  /** \brief Print something about branch - only if log level high
  */
  virtual void print(bool normalBranch);
private:
  /// data
  const CbcClique * clique_;
  /// downMask - bit set to fix to weak bounds, not set to leave unfixed
  unsigned int * downMask_;
  /// upMask - bit set to fix to weak bounds, not set to leave unfixed
  unsigned int * upMask_;
};

/** Branching object for Special ordered sets

    Variable_ is the set id number (redundant, as the object also holds a
    pointer to the set.
 */
class CbcSOSBranchingObject : public CbcBranchingObject {

public:

  // Default Constructor 
  CbcSOSBranchingObject ();

  // Useful constructor
  CbcSOSBranchingObject (CbcModel * model,  const CbcSOS * clique,
			    int way,
			 double separator);
  
  // Copy constructor 
  CbcSOSBranchingObject ( const CbcSOSBranchingObject &);
   
  // Assignment operator 
  CbcSOSBranchingObject & operator=( const CbcSOSBranchingObject& rhs);

  /// Clone
  virtual CbcBranchingObject * clone() const;

  // Destructor 
  virtual ~CbcSOSBranchingObject ();
  
  /// Does next branch and updates state
  virtual double branch(bool normalBranch=false);

  /** \brief Print something about branch - only if log level high
  */
  virtual void print(bool normalBranch);
private:
  /// data
  const CbcSOS * set_;
  /// separator
  double separator_;
};

/** Branching decision default class

  This class implements a simple default algorithm
  (betterBranch()) for choosing a branching variable. 
*/

class CbcBranchDefaultDecision : public CbcBranchDecision {
public:
  // Default Constructor 
  CbcBranchDefaultDecision ();

  // Copy constructor 
  CbcBranchDefaultDecision ( const CbcBranchDefaultDecision &);

  virtual ~CbcBranchDefaultDecision();

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
    \c CbcBranchDefaultDecision object. The nonzero return value is +1 if the
    up branch is preferred, -1 if the down branch is preferred.

    As the names imply, the assumption is that the values supplied for
    \p numInfUp and \p numInfDn will be the number of infeasibilities reported
    by the branching object, and \p changeUp and \p changeDn will be the
    estimated change in objective. Other measures can be used if desired.

    Because an \c CbcBranchDefaultDecision object remembers the current best
    branching candidate (#bestObject_) as well as the values used in the
    comparison, the parameter \p bestSoFar is redundant, hence unused.
 */
  virtual int betterBranch(CbcBranchingObject * thisOne,
			    CbcBranchingObject * bestSoFar,
			    double changeUp, int numInfUp,
			    double changeDn, int numInfDn);

private:
  
  /// Illegal Assignment operator 
  CbcBranchDefaultDecision & operator=(const CbcBranchDefaultDecision& rhs);

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

/** Define a follow on class.
    The idea of this is that in air-crew scheduling problems crew may fly in on flight A
    and out on flight B or on some other flight.  A useful branch is one which on one side
    fixes all which go out on flight B to 0, while the other branch fixes all those that do NOT
    go out on flight B to 0.

    This branching rule should be in addition to normal rules and have a high priority.
*/



class CbcFollowOn : public CbcObject {

public:

  // Default Constructor 
  CbcFollowOn ();

  /** Useful constructor
  */
  CbcFollowOn (CbcModel * model);
  
  // Copy constructor 
  CbcFollowOn ( const CbcFollowOn &);
   
  /// Clone
  virtual CbcObject * clone() const;

  // Assignment operator 
  CbcFollowOn & operator=( const CbcFollowOn& rhs);

  // Destructor 
  ~CbcFollowOn ();
  
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const;

  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion();
  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way) const;
  /// As some computation is needed in more than one place - returns row
  int gutsOfFollowOn(int & otherRow) const;

protected:
  /// data
  /// Matrix
  CoinPackedMatrix matrix_;
  /// Matrix by row
  CoinPackedMatrix matrixByRow_; 
  /// Possible rhs (if 0 then not possible)
  int * rhs_;
};
/** General Branching Object class.
    Each way fixes some variables to lower bound
 */
class CbcFixingBranchingObject : public CbcBranchingObject {

public:

  // Default Constructor 
  CbcFixingBranchingObject ();

  // Useful constructor
  CbcFixingBranchingObject (CbcModel * model, 
			    int way,
			    int numberOnDownSide, const int * down,
			    int numberOnUpSide, const int * up);
  
  // Copy constructor 
  CbcFixingBranchingObject ( const CbcFixingBranchingObject &);
   
  // Assignment operator 
  CbcFixingBranchingObject & operator=( const CbcFixingBranchingObject& rhs);

  /// Clone
  virtual CbcBranchingObject * clone() const;

  // Destructor 
  virtual ~CbcFixingBranchingObject ();
  
  /// Does next branch and updates state
  virtual double branch(bool normalBranch=false);

  /** \brief Print something about branch - only if log level high
  */
  virtual void print(bool normalBranch);
private:
  /// data
  /// Number on down list
  int numberDown_;
  /// Number on up list
  int numberUp_;
  /// downList - variables to fix to lb on down branch
  int * downList_;
  /// upList - variables to fix to lb on up branch
  int * upList_;
};
#endif
