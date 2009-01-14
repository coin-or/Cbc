// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcBranchActual_H
#define CbcBranchActual_H

#include "CbcBranchBase.hpp"
#include "CoinPackedMatrix.hpp"
class CbcIntegerBranchingObject;
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
  virtual ~CbcClique ();
  
  using CbcObject::infeasibility ;
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const;

  using CbcObject::feasibleRegion ;
  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion();

  using CbcObject::createBranch ;
  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way) ;
  /// Number of members
  inline int numberMembers() const
  {return numberMembers_;}

  /// Number of Non SOS members i.e. fixing to zero is strong
  inline int numberNonSOSMembers() const
  {return numberNonSOSMembers_;}

  /// Members (indices in range 0 ... numberIntegers_-1)
  inline const int * members() const
  {return members_;}

  /** Type of each member i.e. which way is strong 0=non SOS, 1 =SOS,
      index is 0 ... numberMembers_-1 */
  inline char type(int index) const
  {if (type_) return type_[index]; else return 1;}

  /// Clique type - 0 <=, 1 == 
  inline int cliqueType() const
  {return cliqueType_;}
  /// Redoes data when sequence numbers change
  virtual void redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns);

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

  /** Useful constructor - which are indices
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
  virtual ~CbcSOS ();
  
  using CbcObject::infeasibility ;
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const;
  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiBranchingInformation * info, 
			       int & preferredWay) const;

  using CbcObject::feasibleRegion ;
  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion();

  using CbcObject::createBranch ;
  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way) ;



  /** Pass in information on branch just done and create CbcObjectUpdateData instance.
      If object does not need data then backward pointer will be NULL.
      Assumes can get information from solver */
  virtual CbcObjectUpdateData createUpdateInformation(const OsiSolverInterface * solver, 
							const CbcNode * node,
							const CbcBranchingObject * branchingObject);
  /// Update object by CbcObjectUpdateData
  virtual void updateInformation(const CbcObjectUpdateData & data) ;
  using CbcObject::solverBranch ;
  /** Create an OsiSolverBranch object

      This returns NULL if branch not represented by bound changes
  */
  virtual OsiSolverBranch * solverBranch() const;
  /// Redoes data when sequence numbers change
  virtual void redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns);
  
  /// Construct an OsiSOS object
  OsiSOS * osiObject(const OsiSolverInterface * solver) const;
  /// Number of members
  inline int numberMembers() const
  {return numberMembers_;}

  /// Members (indices in range 0 ... numberColumns-1)
  inline const int * members() const
  {return members_;}

  /// SOS type
  inline int sosType() const
  {return sosType_;}
  /// Down number times
  inline int numberTimesDown() const
  { return numberTimesDown_;}
  /// Up number times
  inline int numberTimesUp() const
  { return numberTimesUp_;}

  /** Array of weights */
  inline const double * weights() const
  { return weights_;}

  /// Set number of members
  inline void setNumberMembers(int n)
  {numberMembers_ = n;}

  /// Members (indices in range 0 ... numberColumns-1)
  inline int * mutableMembers() const
  {return members_;}

  /** Array of weights */
  inline double * mutableWeights() const
  { return weights_;}

  /** \brief Return true if object can take part in normal heuristics
  */
  virtual bool canDoHeuristics() const 
  {return (sosType_==1&&integerValued_);}
  /// Set whether set is integer valued or not
  inline void setIntegerValued(bool yesNo)
  { integerValued_=yesNo;}
private:
  /// data

  /// Members (indices in range 0 ... numberColumns-1)
  int * members_;
  /// Weights
  double * weights_;
  /// Current pseudo-shadow price estimate down
  mutable double shadowEstimateDown_;
  /// Current pseudo-shadow price estimate up
  mutable double shadowEstimateUp_;
  /// Down pseudo ratio
  double downDynamicPseudoRatio_;
  /// Up pseudo ratio
  double upDynamicPseudoRatio_;
  /// Number of times we have gone down
  int numberTimesDown_;
  /// Number of times we have gone up
  int numberTimesUp_;
  /// Number of members
  int numberMembers_;
  /// SOS type
  int sosType_;
  /// Whether integer valued
  bool integerValued_;
};

/// Define a single integer class


class CbcSimpleInteger : public CbcObject {

public:

  // Default Constructor 
  CbcSimpleInteger ();

  // Useful constructor - passed model and index
  CbcSimpleInteger (CbcModel * model,  int iColumn, double breakEven=0.5);
  
  // Useful constructor - passed model and Osi object
  CbcSimpleInteger (CbcModel * model,  const OsiSimpleInteger * object);
  
  // Copy constructor 
  CbcSimpleInteger ( const CbcSimpleInteger &);
   
  /// Clone
  virtual CbcObject * clone() const;

  // Assignment operator 
  CbcSimpleInteger & operator=( const CbcSimpleInteger& rhs);

  // Destructor 
  virtual ~CbcSimpleInteger ();
  /// Construct an OsiSimpleInteger object
  OsiSimpleInteger * osiObject() const;
  using CbcObject::infeasibility ;
  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiSolverInterface * solver, 
			       const OsiBranchingInformation * info, int & preferredWay) const;

  using CbcObject::feasibleRegion ;
  /** Set bounds to fix the variable at the current (integer) value.

    Given an integer value, set the lower and upper bounds to fix the
    variable. Returns amount it had to move variable.
  */
  virtual double feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const;

  using CbcObject::createBranch ;
  /** Create a branching object and indicate which way to branch first.
      
      The branching object has to know how to create branches (fix
      variables, etc.)
  */
  virtual CbcBranchingObject * createBranch(OsiSolverInterface * solver,
					    const OsiBranchingInformation * info, int way) ;
  /// Fills in a created branching object
  void fillCreateBranch(CbcIntegerBranchingObject * branching, const OsiBranchingInformation * info, int way) ;

  using CbcObject::solverBranch ;
  /** Create an OsiSolverBranch object

      This returns NULL if branch not represented by bound changes
  */
  virtual OsiSolverBranch * solverBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info) const;
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const;

  /** Set bounds to fix the variable at the current (integer) value.

    Given an integer value, set the lower and upper bounds to fix the
    variable. The algorithm takes a bit of care in order to compensate for
    minor numerical inaccuracy.
  */
  virtual void feasibleRegion();

  /** Creates a branching object

    The preferred direction is set by \p way, -1 for down, +1 for up.
  */
  virtual CbcBranchingObject * createBranch(int way) ;
  /** Column number if single column object -1 otherwise,
      so returns >= 0
      Used by heuristics
  */
  virtual int columnNumber() const;
  /// Set column number
  inline void setColumnNumber(int value)
  { columnNumber_ = value;}

  /** Reset variable bounds to their original values.
  
    Bounds may be tightened, so it may be good to be able to set this info in object.
   */
  virtual void resetBounds(const OsiSolverInterface * solver) ;
  /**  Change column numbers after preprocessing
   */
  virtual void resetSequenceEtc(int numberColumns, const int * originalColumns) ;
  /// Original bounds
  inline double originalLowerBound() const
  { return originalLower_;}
  inline void setOriginalLowerBound(double value)
  { originalLower_=value;}
  inline double originalUpperBound() const
  { return originalUpper_;}
  inline void setOriginalUpperBound(double value)
  { originalUpper_=value;}
  /// Breakeven e.g 0.7 -> >= 0.7 go up first
  inline double breakEven() const
  { return breakEven_;}
  /// Set breakeven e.g 0.7 -> >= 0.7 go up first
  inline void setBreakEven(double value)
  { breakEven_=value;}


protected:
  /// data

  /// Original lower bound
  double originalLower_;
  /// Original upper bound
  double originalUpper_;
  /// Breakeven i.e. >= this preferred is up
  double breakEven_;
  /// Column number in model
  int columnNumber_;
  /// If -1 down always chosen first, +1 up always, 0 normal
  int preferredWay_;
};
/** Define an n-way class for variables.
    Only valid value is one at UB others at LB
    Normally 0-1
*/


class CbcNWay : public CbcObject {

public:

  // Default Constructor 
  CbcNWay ();

  /** Useful constructor (which are matrix indices)
  */
  CbcNWay (CbcModel * model, int numberMembers,
	     const int * which, int identifier);
  
  // Copy constructor 
  CbcNWay ( const CbcNWay &);
   
  /// Clone
  virtual CbcObject * clone() const;

  /// Assignment operator 
  CbcNWay & operator=( const CbcNWay& rhs);

  /// Destructor 
  virtual ~CbcNWay ();

  /// Set up a consequence for a single member
  void setConsequence(int iColumn, const CbcConsequence & consequence);
  
  /// Applies a consequence for a single member
  void applyConsequence(int iSequence, int state) const;
  
  using CbcObject::infeasibility ;
  /// Infeasibility - large is 0.5 (and 0.5 will give this)
  virtual double infeasibility(int & preferredWay) const;

  using CbcObject::feasibleRegion ;
  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion();

  using CbcObject::createBranch ;
  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way) ;

  /// Number of members
  inline int numberMembers() const
  {return numberMembers_;}

  /// Members (indices in range 0 ... numberColumns-1)
  inline const int * members() const
  {return members_;}
  /// Redoes data when sequence numbers change
  virtual void redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns);

protected:
  /// data
  /// Number of members
  int numberMembers_;

  /// Members (indices in range 0 ... numberColumns-1)
  int * members_;
  /// Consequences (normally NULL)
  CbcConsequence ** consequence_;
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
  inline const double * downBounds() const
  { return down_;}
  /// Lower and upper bounds for up branch
  inline const double * upBounds() const
  { return up_;}
  /// Set lower and upper bounds for down branch
  inline void setDownBounds(const double bounds[2])
  { memcpy(down_,bounds,2*sizeof(double));}
  /// Set lower and upper bounds for up branch
  inline void setUpBounds(const double bounds[2])
  { memcpy(up_,bounds,2*sizeof(double));}
#ifdef FUNNY_BRANCHING
  /** Which variable (top bit if upper bound changing,
      next bit if on down branch */
  inline const int * variables() const
  { return variables_;}
  // New bound
  inline const double * newBounds() const
  { return newBounds_;}
  /// Number of bound changes
  inline int numberExtraChangedBounds() const
  { return numberExtraChangedBounds_;}
  /// Just apply extra bounds to one variable - COIN_DBL_MAX ignore
  int applyExtraBounds(int iColumn, double lower, double upper, int way) ;
  /// Deactivate bounds for branching
  void deactivate();
  /// Are active bounds for branching
  inline bool active() const
  { return (down_[1]!=-COIN_DBL_MAX);}
#endif

  /** Return the type (an integer identifier) of \c this */
  virtual int type() const { return 100; }

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


/// Define a single integer class but with pseudo costs


class CbcSimpleIntegerPseudoCost : public CbcSimpleInteger {

public:

  // Default Constructor 
  CbcSimpleIntegerPseudoCost ();

  // Useful constructor - passed model index
  CbcSimpleIntegerPseudoCost (CbcModel * model, int iColumn, double breakEven=0.5);
  
  // Useful constructor - passed and model index and pseudo costs
  CbcSimpleIntegerPseudoCost (CbcModel * model, int iColumn, 
			      double downPseudoCost, double upPseudoCost);
  // Useful constructor - passed and model index and pseudo costs
  CbcSimpleIntegerPseudoCost (CbcModel * model, int dummy,int iColumn, 
			      double downPseudoCost, double upPseudoCost);
  
  // Copy constructor 
  CbcSimpleIntegerPseudoCost ( const CbcSimpleIntegerPseudoCost &);
   
  /// Clone
  virtual CbcObject * clone() const;

  // Assignment operator 
  CbcSimpleIntegerPseudoCost & operator=( const CbcSimpleIntegerPseudoCost& rhs);

  // Destructor 
  virtual ~CbcSimpleIntegerPseudoCost ();
  
  using CbcObject::infeasibility ;
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const;

  using CbcObject::createBranch ;
  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way) ;

  /// Down pseudo cost
  inline double downPseudoCost() const
  { return downPseudoCost_;}
  /// Set down pseudo cost
  inline void setDownPseudoCost(double value)
  { downPseudoCost_=value;}

  /// Up pseudo cost
  inline double upPseudoCost() const
  { return upPseudoCost_;}
  /// Set up pseudo cost
  inline void setUpPseudoCost(double value)
  { upPseudoCost_=value;}

  /// Up down separator
  inline double upDownSeparator() const
  { return upDownSeparator_;}
  /// Set up down separator
  inline void setUpDownSeparator(double value)
  { upDownSeparator_=value;}

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

protected:
  /// data

  /// Down pseudo cost
  double downPseudoCost_;
  /// Up pseudo cost
  double upPseudoCost_;
  /** Up/down separator
      If >0.0 then do first branch up if value-floor(value)
      >= this value
  */
  double upDownSeparator_;
  /** Method - 
      0 - normal - return min (up,down)
      1 - if before any solution return CoinMax(up,down)
      2 - if before branched solution return CoinMax(up,down)
      3 - always return CoinMax(up,down)
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
  inline double changeInGuessed() const
  { return changeInGuessed_;}
  /// Set change in guessed
  inline void setChangeInGuessed(double value)
  { changeInGuessed_=value;}

  /** Return the type (an integer identifier) of \c this */
  virtual int type() const { return 101; }

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
  
  using CbcBranchingObject::branch ;
  /// Does next branch and updates state
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
  virtual int type() const { return 102; }

  /** Compare the original object of \c this with the original object of \c
      brObj. Assumes that there is an ordering of the original objects.
      This method should be invoked only if \c this and brObj are of the same
      type. 
      Return negative/0/positive depending on whether \c this is
      smaller/same/larger than the argument.
  */
  virtual int compareOriginalObject(const CbcBranchingObject* brObj) const;

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
  
  using CbcBranchingObject::branch ;
  /// Does next branch and updates state
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
  virtual int type() const { return 103; }

  /** Compare the original object of \c this with the original object of \c
      brObj. Assumes that there is an ordering of the original objects.
      This method should be invoked only if \c this and brObj are of the same
      type. 
      Return negative/0/positive depending on whether \c this is
      smaller/same/larger than the argument.
  */
  virtual int compareOriginalObject(const CbcBranchingObject* brObj) const;

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
  
  using CbcBranchingObject::branch ;
  /// Does next branch and updates state
  virtual double branch();
  /** Update bounds in solver as in 'branch' and update given bounds.
      branchState is -1 for 'down' +1 for 'up' */
  virtual void fix(OsiSolverInterface * solver,
		   double * lower, double * upper,
		   int branchState) const ;

  /** Reset every information so that the branching object appears to point to
      the previous child. This method does not need to modify anything in any
      solver. */
  virtual void previousBranch() {
    CbcBranchingObject::previousBranch();
    computeNonzeroRange();
  }

  using CbcBranchingObject::print ;
  /** \brief Print something about branch - only if log level high
  */
  virtual void print();

  /** Return the type (an integer identifier) of \c this */
  virtual int type() const { return 104; }

  /** Compare the original object of \c this with the original object of \c
      brObj. Assumes that there is an ordering of the original objects.
      This method should be invoked only if \c this and brObj are of the same
      type. 
      Return negative/0/positive depending on whether \c this is
      smaller/same/larger than the argument.
  */
  virtual int compareOriginalObject(const CbcBranchingObject* brObj) const;

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

  /** Fill out the \c firstNonzero_ and \c lastNonzero_ data members */
  void computeNonzeroRange();

private:
  /// data
  const CbcSOS * set_;
  /// separator
  double separator_;
  /** The following two members describe the range in the members_ of the
      original object that whose upper bound is not fixed to 0. This is not
      necessary for Cbc to function correctly, this is there for heuristics so
      that separate branching decisions on the same object can be pooled into
      one branching object. */
  int firstNonzero_;
  int lastNonzero_;
};

/** N way branching Object class.
    Variable is number of set.
 */
class CbcNWayBranchingObject : public CbcBranchingObject {

public:

  // Default Constructor 
  CbcNWayBranchingObject ();

  /** Useful constructor - order had matrix indices
      way_ -1 corresponds to setting first, +1 to second, +3 etc.
      this is so -1 and +1 have similarity to normal
  */
  CbcNWayBranchingObject (CbcModel * model,  const CbcNWay * nway,
                          int numberBranches, const int * order);
  
  // Copy constructor 
  CbcNWayBranchingObject ( const CbcNWayBranchingObject &);
   
  // Assignment operator 
  CbcNWayBranchingObject & operator=( const CbcNWayBranchingObject& rhs);

  /// Clone
  virtual CbcBranchingObject * clone() const;

  // Destructor 
  virtual ~CbcNWayBranchingObject ();
  
  using CbcBranchingObject::branch ;
  /// Does next branch and updates state
  virtual double branch();

#if 0
  // FIXME: what do we need to do here?
  /** Reset every information so that the branching object appears to point to
      the previous child. This method does not need to modify anything in any
      solver. */
  virtual void previousBranch();
#endif

  using CbcBranchingObject::print ;
  /** \brief Print something about branch - only if log level high
  */
  virtual void print();
  /** The number of branch arms created for this branching object
  */
  virtual int numberBranches() const
  {return numberInSet_;}
  /// Is this a two way object (-1 down, +1 up)
  virtual bool twoWay() const
  { return false;}

  /** Return the type (an integer identifier) of \c this */
  virtual int type() const { return 105; }

  /** Compare the original object of \c this with the original object of \c
      brObj. Assumes that there is an ordering of the original objects.
      This method should be invoked only if \c this and brObj are of the same
      type. 
      Return negative/0/positive depending on whether \c this is
      smaller/same/larger than the argument.
  */
  virtual int compareOriginalObject(const CbcBranchingObject* brObj) const;

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

private:
  /// order of branching - points back to CbcNWay
  int * order_;
  /// Points back to object
  const CbcNWay * object_;
  /// Number in set
  int numberInSet_;
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
  /** Sets or gets best criterion so far */
  virtual void setBestCriterion(double value);
  virtual double getBestCriterion() const;

  /** \brief Compare N branching objects. Return index of best
      and sets way of branching in chosen object.
    
    This routine is used only after strong branching.
  */

  virtual int
  bestBranch (CbcBranchingObject ** objects, int numberObjects, int numberUnsatisfied,
	      double * changeUp, int * numberInfeasibilitiesUp,
	      double * changeDown, int * numberInfeasibilitiesDown,
	      double objectiveValue) ;
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

  /// Pointer to best branching object
  CbcBranchingObject * bestObject_;

  /// Number of infeasibilities for down
  int bestNumberDown_;

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
  
  using CbcObject::infeasibility ;
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const;

  using CbcObject::feasibleRegion ;
  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion();

  using CbcObject::createBranch ;
  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way) ;
  /// As some computation is needed in more than one place - returns row
  virtual int gutsOfFollowOn(int & otherRow, int & preferredWay) const;

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
  
  using CbcBranchingObject::branch ;
  /// Does next branch and updates state
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
  virtual int type() const { return 106; }

  /** Compare the original object of \c this with the original object of \c
      brObj. Assumes that there is an ordering of the original objects.
      This method should be invoked only if \c this and brObj are of the same
      type. 
      Return negative/0/positive depending on whether \c this is
      smaller/same/larger than the argument.
  */
  virtual int compareOriginalObject(const CbcBranchingObject* brObj) const;

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
/** Class for consequent bounds.
    When a variable is branched on it normally interacts with other variables by
    means of equations.  There are cases where we want to step outside LP and do something
    more directly e.g. fix bounds.  This class is for that.

    A state of -9999 means at LB, +9999 means at UB,
    others mean if fixed to that value.

 */

class CbcFixVariable : public CbcConsequence {

public:

  // Default Constructor 
  CbcFixVariable ();

  // One useful Constructor 
  CbcFixVariable (int numberStates,const int * states, const int * numberNewLower, const int ** newLowerValue,
                  const int ** lowerColumn,
                  const int * numberNewUpper, const int ** newUpperValue,
                  const int ** upperColumn);

  // Copy constructor 
  CbcFixVariable ( const CbcFixVariable & rhs);
   
  // Assignment operator 
  CbcFixVariable & operator=( const CbcFixVariable & rhs);

  /// Clone
  virtual CbcConsequence * clone() const;

  /// Destructor 
  virtual ~CbcFixVariable ();

  /** Apply to an LP solver.  Action depends on state
   */
  virtual void applyToSolver(OsiSolverInterface * solver, int state) const;
  
protected:
  /// Number of states
  int numberStates_;
  /// Values of integers for various states
  int * states_;
  /// Start of information for each state (setting new lower)
  int * startLower_;
  /// Start of information for each state (setting new upper)
  int * startUpper_;
  /// For each variable new bounds
  double * newBound_;
  /// Variable
  int * variable_;
};
/** Dummy branching object

  This object specifies a one-way dummy branch.
  This is so one can carry on branching even when it looks feasible
*/

class CbcDummyBranchingObject : public CbcBranchingObject {

public:

  /// Default constructor 
  CbcDummyBranchingObject (CbcModel * model=NULL);

  /// Copy constructor 
  CbcDummyBranchingObject ( const CbcDummyBranchingObject &);
   
  /// Assignment operator 
  CbcDummyBranchingObject & operator= (const CbcDummyBranchingObject& rhs);

  /// Clone
  virtual CbcBranchingObject * clone() const;

  /// Destructor 
  virtual ~CbcDummyBranchingObject ();
  
  using CbcBranchingObject::branch ;
  /** \brief Dummy branch
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
  virtual int type() const { return 107; }

  /** Compare the original object of \c this with the original object of \c
      brObj. Assumes that there is an ordering of the original objects.
      This method should be invoked only if \c this and brObj are of the same
      type. 
      Return negative/0/positive depending on whether \c this is
      smaller/same/larger than the argument.
  */
  virtual int compareOriginalObject(const CbcBranchingObject* brObj) const;

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

};

/** Define a catch all class.
    This will create a list of subproblems
*/


class CbcGeneral : public CbcObject {

public:

  // Default Constructor 
  CbcGeneral ();

  /** Useful constructor
      Just needs to point to model.
  */
  CbcGeneral (CbcModel * model);
  
  // Copy constructor 
  CbcGeneral ( const CbcGeneral &);
   
  /// Clone
  virtual CbcObject * clone() const=0;

  // Assignment operator 
  CbcGeneral & operator=( const CbcGeneral& rhs);

  // Destructor 
  ~CbcGeneral ();
  
  using CbcObject::infeasibility ;
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const=0;

  using CbcObject::feasibleRegion ;
  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion()=0;

  using CbcObject::createBranch ;
  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way)=0 ;

  /// Redoes data when sequence numbers change
  virtual void redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns)=0;

protected:
  /// data
};
#ifdef COIN_HAS_CLP

/** Define a catch all class.
    This will create a list of subproblems using partial evaluation
*/
#include "ClpSimplex.hpp"
#include "ClpNode.hpp"

class CbcGeneralDepth : public CbcGeneral {

public:

  // Default Constructor 
  CbcGeneralDepth ();

  /** Useful constructor
      Just needs to point to model.
      Initial version does evaluation to depth N
      This is stored in CbcModel but may be
      better here
  */
  CbcGeneralDepth (CbcModel * model, int maximumDepth);
  
  // Copy constructor 
  CbcGeneralDepth ( const CbcGeneralDepth &);
   
  /// Clone
  virtual CbcObject * clone() const;

  // Assignment operator 
  CbcGeneralDepth & operator=( const CbcGeneralDepth& rhs);

  // Destructor 
  ~CbcGeneralDepth ();
  
  using CbcObject::infeasibility ;
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const;

  using CbcObject::feasibleRegion ;
  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion();

  using CbcObject::createBranch ;
  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way) ;
  /// Return maximum number of nodes
  int maximumNodes() const;
  /// Get maximum depth
  inline int maximumDepth() const
  {return maximumDepth_;}
  /// Set maximum depth
  inline void setMaximumDepth(int value)
  {maximumDepth_ = value;}
  /// Get which solution
  inline int whichSolution() const
  {return whichSolution_;}
  /// Get ClpNode info
  inline ClpNode * nodeInfo(int which)
  { return nodeInfo_->nodeInfo_[which];}

  /// Redoes data when sequence numbers change
  virtual void redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns);

protected:
  /// data
  /// Maximum depth
  int maximumDepth_;
  /// Which node has solution (or -1)
  mutable int whichSolution_;
  /// Number of valid nodes (including whichSolution_)
  mutable int numberNodes_;
  /// For solving nodes
  mutable ClpNodeStuff * nodeInfo_;
};

/** Defines a general subproblem
    Basis will be made more compact later
*/
class CoinWarmStartDiff;
class CbcSubProblem {

public:

  /// Default constructor 
  CbcSubProblem ();

  /// Constructor from model
  CbcSubProblem (const OsiSolverInterface * solver,
		 const double * lowerBefore,
		 const double * upperBefore,
		 const unsigned char * status,
		 int depth);

  /// Copy constructor 
  CbcSubProblem ( const CbcSubProblem &);
   
  /// Assignment operator 
  CbcSubProblem & operator= (const CbcSubProblem& rhs);

  /// Destructor 
  virtual ~CbcSubProblem ();

  /// Apply subproblem (1=bounds, 2=basis, 3=both)
  void apply(OsiSolverInterface * model, int what=3) const;

public:
  /// Value of objective
  double objectiveValue_;
  /// Sum of infeasibilities
  double sumInfeasibilities_;
  /** Which variable (top bit if upper bound changing)
      next bit if changing on down branch only */
  int * variables_;
  /// New bound
  double * newBounds_;
  /// Status 
  mutable CoinWarmStartBasis * status_;
  /// Depth
  int depth_;
  /// Number of Extra bound changes
  int numberChangedBounds_;
  /// Number of infeasibilities
  int numberInfeasibilities_;
};

/** Branching object for general objects

 */
class CbcNode;
class CbcGeneralBranchingObject : public CbcBranchingObject {

public:

  // Default Constructor 
  CbcGeneralBranchingObject ();

  // Useful constructor
  CbcGeneralBranchingObject (CbcModel * model);
  
  // Copy constructor 
  CbcGeneralBranchingObject ( const CbcGeneralBranchingObject &);
   
  // Assignment operator 
  CbcGeneralBranchingObject & operator=( const CbcGeneralBranchingObject& rhs);

  /// Clone
  virtual CbcBranchingObject * clone() const;

  // Destructor 
  virtual ~CbcGeneralBranchingObject ();
  
  using CbcBranchingObject::branch ;
  /// Does next branch and updates state
  virtual double branch();
  /** Double checks in case node can change its mind!
      Can change objective etc */
  virtual void checkIsCutoff(double cutoff);

  using CbcBranchingObject::print ;
  /** \brief Print something about branch - only if log level high
  */
  virtual void print();
  /// Fill in current objective etc
  void state(double & objectiveValue,double & sumInfeasibilities,
	     int & numberUnsatisfied,int which) const;
  /// Set CbcNode
  inline void setNode(CbcNode * node)
  { node_ = node;}
  /** Return the type (an integer identifier) of \c this */
  virtual int type() const { return 108; }

  /** Compare the original object of \c this with the original object of \c
      brObj. Assumes that there is an ordering of the original objects.
      This method should be invoked only if \c this and brObj are of the same
      type. 
      Return negative/0/positive depending on whether \c this is
      smaller/same/larger than the argument.
  */
  virtual int compareOriginalObject(const CbcBranchingObject* brObj) const;

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
  /// Number of subproblems
  inline int numberSubProblems() const
  { return numberSubProblems_;}
  /// Decrement number left and return number
  inline int decrementNumberLeft()
  { numberSubLeft_--; return numberSubLeft_;}
  /// Which node we want to use
  inline int whichNode() const
  { return whichNode_;}
  /// Set which node we want to use
  inline void setWhichNode(int value)
  { whichNode_ = value;}
  // Sub problem
  const CbcSubProblem * subProblem(int which) const
  { return subProblems_+which;}

public:
  /// data
  // Sub problems
  CbcSubProblem * subProblems_;
  /// Node
  CbcNode * node_;
  /// Number of subproblems
  int numberSubProblems_;
  /// Number of subproblems left
  int numberSubLeft_;
  /// Which node we want to use (-1 for default)
  int whichNode_;
  /// Number of rows
  int numberRows_;
};

/** Branching object for general objects - just one

 */
class CbcOneGeneralBranchingObject : public CbcBranchingObject {

public:

  // Default Constructor 
  CbcOneGeneralBranchingObject ();

  // Useful constructor
  CbcOneGeneralBranchingObject (CbcModel * model,
				 CbcGeneralBranchingObject * object,
				 int whichOne);
  
  // Copy constructor 
  CbcOneGeneralBranchingObject ( const CbcOneGeneralBranchingObject &);
   
  // Assignment operator 
  CbcOneGeneralBranchingObject & operator=( const CbcOneGeneralBranchingObject& rhs);

  /// Clone
  virtual CbcBranchingObject * clone() const;

  // Destructor 
  virtual ~CbcOneGeneralBranchingObject ();
  
  using CbcBranchingObject::branch ;
  /// Does next branch and updates state
  virtual double branch();
  /** Double checks in case node can change its mind!
      Can change objective etc */
  virtual void checkIsCutoff(double cutoff);

  using CbcBranchingObject::print ;
  /** \brief Print something about branch - only if log level high
  */
  virtual void print();
  /** Return the type (an integer identifier) of \c this */
  virtual int type() const { return 110; }

  /** Compare the original object of \c this with the original object of \c
      brObj. Assumes that there is an ordering of the original objects.
      This method should be invoked only if \c this and brObj are of the same
      type. 
      Return negative/0/positive depending on whether \c this is
      smaller/same/larger than the argument.
  */
  virtual int compareOriginalObject(const CbcBranchingObject* brObj) const;

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

public:
  /// data
  /// Object
  CbcGeneralBranchingObject * object_;
  /// Which one
  int whichOne_;
};
#endif
#endif
