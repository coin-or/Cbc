// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcBranchBase_H
#define CbcBranchBase_H

#include <string>
#include <vector>
#include "OsiBranchingObject.hpp"
class OsiSolverInterface;
class OsiSolverBranch;

class CbcModel;
class CbcNode;
class CbcNodeInfo;
class CbcBranchingObject;
class OsiChooseVariable;
class CbcObjectUpdateData;

//#############################################################################

enum CbcRangeCompare {
  CbcRangeSame,
  CbcRangeDisjoint,
  CbcRangeSubset,
  CbcRangeSuperset,
  CbcRangeOverlap
};

//#############################################################################

/** Abstract base class for `objects'.
    It now just has stuff that OsiObject does not have

  The branching model used in Cbc is based on the idea of an <i>object</i>.
  In the abstract, an object is something that has a feasible region, can be
  evaluated for infeasibility, can be branched on (<i>i.e.</i>, there's some
  constructive action to be taken to move toward feasibility), and allows
  comparison of the effect of branching.

  This class (CbcObject) is the base class for an object. To round out the
  branching model, the class CbcBranchingObject describes how to perform a
  branch, and the class CbcBranchDecision describes how to compare two
  CbcBranchingObjects.

  To create a new type of object you need to provide three methods:
  #infeasibility(), #feasibleRegion(), and #createBranch(), described below.

  This base class is primarily virtual to allow for any form of structure.
  Any form of discontinuity is allowed.

  \todo The notion that all branches are binary (two arms) is wired into the
	implementation of CbcObject, CbcBranchingObject, and
	CbcBranchDecision. Changing this will require a moderate amount of
	recoding.
 */
// This can be used if object wants to skip strong branching
  typedef struct {
    CbcBranchingObject * possibleBranch; // what a branch would do
    double upMovement; // cost going up (and initial away from feasible)
    double downMovement; // cost going down
    int numIntInfeasUp ; // without odd ones
    int numObjInfeasUp ; // just odd ones
    bool finishedUp; // true if solver finished
    int numItersUp ; // number of iterations in solver
    int numIntInfeasDown ; // without odd ones
    int numObjInfeasDown ; // just odd ones
    bool finishedDown; // true if solver finished
    int numItersDown; // number of iterations in solver
    int objectNumber; // Which object it is
    int fix; // 0 if no fix, 1 if we can fix up, -1 if we can fix down
  } CbcStrongInfo;

class CbcObject : public OsiObject {

public:

  // Default Constructor 
  CbcObject ();

  // Useful constructor
  CbcObject (CbcModel * model);
  
  // Copy constructor 
  CbcObject ( const CbcObject &);
   
  // Assignment operator 
  CbcObject & operator=( const CbcObject& rhs);

  /// Clone
  virtual CbcObject * clone() const=0;

  /// Destructor 
  virtual ~CbcObject ();

  /** Infeasibility of the object

      This is some measure of the infeasibility of the object. It should be
      scaled to be in the range [0.0, 0.5], with 0.0 indicating the object
      is satisfied.

      The preferred branching direction is returned in preferredWay,

      This is used to prepare for strong branching but should also think of
      case when no strong branching
      
      The object may also compute an estimate of cost of going "up" or "down".
      This will probably be based on pseudo-cost ideas
  */
  virtual double infeasibility(int &preferredWay) const = 0;
  /// Dummy one for compatibility
  virtual double infeasibility(const OsiBranchingInformation * info,
			       int &preferredWay) const;

  /** For the variable(s) referenced by the object,
      look at the current solution and set bounds to match the solution.
  */
  virtual void feasibleRegion() = 0;
  /// Dummy one for compatibility
  virtual double feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const;

  /** Create a branching object and indicate which way to branch first.

      The branching object has to know how to create branches (fix
      variables, etc.)
  */
  virtual CbcBranchingObject * createBranch(int way) = 0;
  
  /** Infeasibility of the object
      
    This is some measure of the infeasibility of the object. 0.0 
    indicates that the object is satisfied.
  
    The preferred branching direction is returned in way,
  
    This is used to prepare for strong branching but should also think of
    case when no strong branching
  
    The object may also compute an estimate of cost of going "up" or "down".
    This will probably be based on pseudo-cost ideas

    This should also set mutable infeasibility_ and whichWay_
    This is for instant re-use for speed
  */
  virtual double infeasibility(const OsiSolverInterface * solver,int &preferredWay) const;
  
  /** For the variable(s) referenced by the object,
      look at the current solution and set bounds to match the solution.
      Returns measure of how much it had to move solution to make feasible
  */
  virtual double feasibleRegion(OsiSolverInterface * solver) const ;
  
  /** Create a branching object and indicate which way to branch first.
      
      The branching object has to know how to create branches (fix
      variables, etc.)
  */
  virtual OsiBranchingObject * createBranch(OsiSolverInterface * solver, int way) const;
  /** Create a branching object and indicate which way to branch first.
      
      The branching object has to know how to create branches (fix
      variables, etc.)
  */
  virtual OsiBranchingObject * createBranch(OsiSolverInterface * solver,const OsiBranchingInformation * info, int way) const;
  /** Create an OsiSolverBranch object

      This returns NULL if branch not represented by bound changes
  */
  virtual OsiSolverBranch * solverBranch() const;
  
  /** \brief Given a valid solution (with reduced costs, etc.),
      return a branching object which would give a new feasible
      point in a good direction.

      If the method cannot generate a feasible point (because there aren't
      any, or because it isn't bright enough to find one), it should
      return null.
  */
  virtual CbcBranchingObject * preferredNewFeasible() const 
  { return NULL;}
  
  /** \brief Given a valid solution (with reduced costs, etc.),
      return a branching object which would give a new feasible
      point in a bad direction.

      If the method cannot generate a feasible point (because there aren't
      any, or because it isn't bright enough to find one), it should
      return null.
  */
  virtual CbcBranchingObject * notPreferredNewFeasible() const 
  { return NULL;}

  /** Reset variable bounds to their original values.
  
    Bounds may be tightened, so it may be good to be able to set this info in object.
   */
  virtual void resetBounds(const OsiSolverInterface * solver) {}
  
  /** Returns floor and ceiling i.e. closest valid points
  */
  virtual void floorCeiling(double & floorValue, double & ceilingValue, double value,
			    double tolerance) const;

  /** Pass in information on branch just done and create CbcObjectUpdateData instance.
      If object does not need data then backward pointer will be NULL.
      Assumes can get information from solver */
  virtual CbcObjectUpdateData createUpdateInformation(const OsiSolverInterface * solver, 
							const CbcNode * node,
							const CbcBranchingObject * branchingObject);

  /// Update object by CbcObjectUpdateData
  virtual void updateInformation(const CbcObjectUpdateData & data) {}

  /// Identifier (normally column number in matrix)
  inline int id() const
  { return id_;}
  
   /// update model
  inline void setModel(CbcModel * model)
  { model_ = model;}
  
  /// Return model
  inline CbcModel * model() const
  {return  model_;}

  /// If -1 down always chosen first, +1 up always, 0 normal
  inline int preferredWay() const
  { return preferredWay_;}
  /// Set -1 down always chosen first, +1 up always, 0 normal
  inline void setPreferredWay(int value)
  { preferredWay_=value;}
  /// Redoes data when sequence numbers change
  virtual void redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns) {}
  
protected:
  /// data

  /// Model
  CbcModel * model_;
  /// Identifier (normally column number in matrix)
  int id_;
  /// If -1 down always chosen first, +1 up always, 0 normal
  int preferredWay_;

};

/** \brief Abstract branching object base class
    Now just difference with OsiBranchingObject

  In the abstract, an CbcBranchingObject contains instructions for how to
  branch. We want an abstract class so that we can describe how to branch on
  simple objects (<i>e.g.</i>, integers) and more exotic objects
  (<i>e.g.</i>, cliques or hyperplanes).

  The #branch() method is the crucial routine: it is expected to be able to
  step through a set of branch arms, executing the actions required to create
  each subproblem in turn. The base class is primarily virtual to allow for
  a wide range of problem modifications.

  See CbcObject for an overview of the three classes (CbcObject,
  CbcBranchingObject, and CbcBranchDecision) which make up cbc's branching
  model.
*/

class CbcBranchingObject : public OsiBranchingObject {

public:

  /// Default Constructor 
  CbcBranchingObject ();

  /// Constructor 
  CbcBranchingObject (CbcModel * model, int variable, int way , double value);
  
  /// Copy constructor 
  CbcBranchingObject ( const CbcBranchingObject &);
   
  /// Assignment operator 
  CbcBranchingObject & operator=( const CbcBranchingObject& rhs);

  /// Clone
  virtual CbcBranchingObject * clone() const=0;

  /// Destructor 
  virtual ~CbcBranchingObject ();

  /** Some branchingObjects may claim to be able to skip
      strong branching.  If so they ahve to fill in CbcStrongInfo.
      The object mention in incoming CbcStrongInfo must match.
      Returns nonzero if skip is wanted */
  virtual int fillStrongInfo( CbcStrongInfo & info) {return 0;}
  /// Reset number of branches left to original
  inline void resetNumberBranchesLeft()
  { branchIndex_=0;}

  /** \brief Execute the actions required to branch, as specified by the
	     current state of the branching object, and advance the object's
	     state.  Mainly for diagnostics, whether it is true branch or
	     strong branching is also passed.
	     Returns change in guessed objective on next branch
  */
  virtual double branch()=0;
  /** \brief Execute the actions required to branch, as specified by the
	     current state of the branching object, and advance the object's
	     state.  Mainly for diagnostics, whether it is true branch or
	     strong branching is also passed.
	     Returns change in guessed objective on next branch
  */
  virtual double branch(OsiSolverInterface * solver)
  { return branch();}

  /** Reset every information so that the branching object appears to point to
      the previous child. This method does not need to modify anything in any
      solver. */
  virtual void previousBranch() {
    assert(branchIndex_ > 0);
    branchIndex_--;
    way_ = -way_;
  }

  using OsiBranchingObject::print ;
  /** \brief Print something about branch - only if log level high
  */
  virtual void print() const {}

  /** \brief Index identifying the associated CbcObject within its class.
  
    The name is misleading, and typically the index will <i>not</i> refer
    directly to a variable.
    Rather, it identifies an CbcObject within the class of similar
    CbcObjects
    
    <i>E.g.</i>, for an CbcSimpleInteger, variable() is the index of the
    integer variable in the set of integer variables (<i>not</i> the index of
    the variable in the set of all variables).
  */
  inline int variable() const
  {return variable_;}

  /** Get the state of the branching object
  
    Returns a code indicating the active arm of the branching object.
    The precise meaning is defined in the derived class.

    \sa #way_
  */
  inline int way() const
  {return way_;}

  /** Set the state of the branching object.

    See #way()
  */
  inline void way(int way)
  {way_=way;}

   /// update model
  inline void setModel(CbcModel * model)
  { model_ = model;}
  /// Return model
  inline CbcModel * model() const
  {return  model_;}

  /// Return pointer back to object which created
  inline CbcObject * object() const
  {return  originalCbcObject_;}
  /// Set pointer back to object which created
  inline void setOriginalObject(CbcObject * object)
  {originalCbcObject_=object;}

  // Methods used in heuristics
  
  /** Return the type (an integer identifier) of \c this */
  virtual int type() const = 0;

  /** Compare the original object of \c this with the original object of \c
      brObj. Assumes that there is an ordering of the original objects.
      This method should be invoked only if \c this and brObj are of the same
      type. 
      Return negative/0/positive depending on whether \c this is
      smaller/same/larger than the argument.
  */
  virtual int compareOriginalObject(const CbcBranchingObject* brObj) const
  {
    const CbcBranchingObject* br=dynamic_cast<const CbcBranchingObject*>(brObj);
    return variable() - br->variable();
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
  (const CbcBranchingObject* brObj, const bool replaceIfOverlap = false) = 0;

protected:

  /// The model that owns this branching object
  CbcModel * model_;
  /// Pointer back to object which created
  CbcObject * originalCbcObject_;

  /// Branching variable (0 is first integer)
  int variable_;
  // was - Way to branch - -1 down (first), 1 up, -2 down (second), 2 up (second)
  /** The state of the branching object.

    Specifies the active arm of the branching object. Coded as -1 to take
    the `down' arm, +1 for the `up' arm. `Down' and `up' are defined based on
    the natural meaning (floor and ceiling, respectively) for a simple integer.
    The precise meaning is defined in the derived class.
  */
  int way_;

};


/** Abstract branching decision base class

  In the abstract, an CbcBranchDecision object is expected to be able to
  compare two possible branching choices.

  The #betterBranch() method is the crucial routine. It is expected to be able
  to compare two \link CbcBranchingObject CbcBranchingObjects \endlink.

  See CbcObject for an overview of the three classes (CbcObject,
  CbcBranchingObject, and CbcBranchDecision) which make up cbc's branching
  model.
*/

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
  virtual int whichMethod() {return 2;}

  /** Saves a clone of current branching object.  Can be used to update
      information on object causing branch - after branch */
  virtual void saveBranchingObject(OsiBranchingObject * object) {}
  /** Pass in information on branch just done.
      assumes object can get information from solver */
  virtual void updateInformation(OsiSolverInterface * solver, 
                                 const CbcNode * node) {}
  /** Sets or gets best criterion so far */
  virtual void setBestCriterion(double value) {}
  virtual double getBestCriterion() const {return 0.0;}
  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) {}
  /// Model
  inline CbcModel * cbcModel() const
  { return model_;}
  /* If chooseMethod_ id non-null then the rest is fairly pointless
     as choosemethod_ will be doing all work
  */
  OsiChooseVariable * chooseMethod() const
  { return chooseMethod_;}
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
/** Abstract base class for consequent bounds.
    When a variable is branched on it normally interacts with other variables by
    means of equations.  There are cases where we want to step outside LP and do something
    more directly e.g. fix bounds.  This class is for that.

    At present it need not be virtual as only instance is CbcFixVariable, but ...

 */

class CbcConsequence {

public:

  // Default Constructor 
  CbcConsequence ();

  // Copy constructor 
  CbcConsequence ( const CbcConsequence & rhs);
   
  // Assignment operator 
  CbcConsequence & operator=( const CbcConsequence & rhs);

  /// Clone
  virtual CbcConsequence * clone() const=0;

  /// Destructor 
  virtual ~CbcConsequence ();

  /** Apply to an LP solver.  Action depends on state
   */
  virtual void applyToSolver(OsiSolverInterface * solver, int state) const=0;
  
protected:
};
/*  This stores data so an object can be updated
 */
class CbcObjectUpdateData {

public:

  /// Default Constructor 
  CbcObjectUpdateData ();

  /// Useful constructor
  CbcObjectUpdateData (CbcObject * object,
		       int way,
		       double change,
		       int status,
		       int intDecrease_,
		       double branchingValue);
  
  /// Copy constructor 
  CbcObjectUpdateData ( const CbcObjectUpdateData &);
   
  /// Assignment operator 
  CbcObjectUpdateData & operator=( const CbcObjectUpdateData& rhs);

  /// Destructor 
  virtual ~CbcObjectUpdateData ();

  
public:
  /// data

  /// Object
  CbcObject * object_;
  /// Branch as defined by instance of CbcObject
  int way_;
  /// Object number
  int objectNumber_;
  /// Change in objective
  double change_;
  /// Status 0 Optimal, 1 infeasible, 2 unknown
  int status_;
  /// Decrease in number unsatisfied
  int intDecrease_;
  /// Branching value
  double branchingValue_;
  /// Objective value before branching
  double originalObjective_;
  /// Current cutoff
  double cutoff_;

};

//##############################################################################

/** Compare two ranges. The two bounds arrays are both of size two and
    describe closed intervals. Return the appropriate CbcRangeCompare value
    (first argument being the sub/superset if that's the case). In case of
    overlap (and if \c replaceIfOverlap is true) replace the content of thisBd
    with the intersection of the ranges.
*/
static inline CbcRangeCompare
CbcCompareRanges(double* thisBd, const double* otherBd,
		 const bool replaceIfOverlap)
{
  const double lbDiff = thisBd[0] - otherBd[0];
  if (lbDiff < 0) { // lb of this < lb of other
    if (thisBd[1] >= otherBd[1]) { // ub of this >= ub of other
      return CbcRangeSuperset;
    } else if (thisBd[1] < otherBd[0]) {
      return CbcRangeDisjoint;
    } else {
      // overlap
      if (replaceIfOverlap) {
	thisBd[0] = otherBd[0];
      }
      return CbcRangeOverlap;
    }
  } else if (lbDiff > 0) { // lb of this > lb of other
    if (thisBd[1] <= otherBd[1]) { // ub of this <= ub of other
      return CbcRangeSubset;
    } else if (thisBd[0] > otherBd[1]) {
      return CbcRangeDisjoint;
    } else {
      // overlap
      if (replaceIfOverlap) {
	thisBd[1] = otherBd[1];
      }
      return CbcRangeOverlap;
    }
  } else { // lb of this == lb of other
    if (thisBd[1] == otherBd[1]) {
      return CbcRangeSame;
    }
    return thisBd[1] < otherBd[1] ? CbcRangeSubset : CbcRangeSuperset;
  }

  return CbcRangeSame; // fake return

}

//#############################################################################

#endif
