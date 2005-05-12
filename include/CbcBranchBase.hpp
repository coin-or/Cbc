// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcBranchBase_H
#define CbcBranchBase_H

#include <string>
#include <vector>

class OsiSolverInterface;

class CbcModel;
class CbcNode;
class CbcNodeInfo;
class CbcBranchingObject;

//#############################################################################

/** Abstract base class for `objects'.

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

class CbcObject {

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

  /** For the variable(s) referenced by the object,
      look at the current solution and set bounds to match the solution.
  */
  virtual void feasibleRegion() = 0;

  /** Create a branching object and indicate which way to branch first.

      The branching object has to know how to create branches (fix
      variables, etc.)
  */
  virtual CbcBranchingObject * createBranch(int way) = 0;
  
  /** \brief Given a valid solution (with reduced costs, etc.),
      return a branching object which would give a new feasible
      point in a good direction.

      If the method cannot generate a feasible point (because there aren't
      any, or because it isn't bright enough to find one), it should
      return null.
  */
  virtual CbcBranchingObject * preferredNewFeasible() const 
  { return NULL;};
  
  /** \brief Given a valid solution (with reduced costs, etc.),
      return a branching object which would give a new feasible
      point in a bad direction.

      If the method cannot generate a feasible point (because there aren't
      any, or because it isn't bright enough to find one), it should
      return null.
  */
  virtual CbcBranchingObject * notPreferredNewFeasible() const 
  { return NULL;};

  /** Reset variable bounds to their original values.
  
    Bounds may be tightened, so it may be good to be able to reset them to
    their original values.
   */
  virtual void resetBounds() {};
  
  /** \brief Return true if branch created by object should fix variables
  */
  virtual bool boundBranch() const 
  {return true;};
  /** Returns floor and ceiling i.e. closest valid points
  */
  virtual void floorCeiling(double & floorValue, double & ceilingValue, double value,
			    double tolerance) const;

  /// Identifier (normally column number in matrix)
  inline int id() const
  { return id_;};
  /// Return Priority
  inline int priority() const
  { return priority_;};
  /// Set priority
  inline void setPriority(int priority)
  { priority_ = priority;};
  
  /// Column number if single column object -1 otherwise
  virtual int columnNumber() const;
  
   /// update model
  inline void setModel(CbcModel * model)
  { model_ = model;};
  
  /// Return model
  inline CbcModel * model() const
  {return  model_;};

  /// Return "up" estimate (default 1.0e-5)
  virtual double upEstimate() const;
  /// Return "down" estimate (default 1.0e-5)
  virtual double downEstimate() const;
  
protected:
  /// data

  /// Model
  CbcModel * model_;
  /// Identifier (normally column number in matrix)
  int id_;
  /// Priority
  int priority_;

};

/** \brief Abstract branching object base class

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

class CbcBranchingObject {

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
  virtual int fillStrongInfo( CbcStrongInfo & info) {return 0;};
  /** The number of branch arms created for this branching object

    \todo The hardwired `2' has to be changed before cbc can do branches with
	  more than two arms.
  */
  virtual int numberBranches() const
  {return 2;};

  /// The number of branch arms left to be evaluated
  virtual int numberBranchesLeft() const
  {return numberBranchesLeft_;};

  /** \brief Execute the actions required to branch, as specified by the
	     current state of the branching object, and advance the object's
	     state.  Mainly for diagnostics, whether it is true branch or
	     strong branching is also passed.
	     Returns change in guessed objective on next branch
  */
  virtual double branch(bool normalBranch=false)=0;

  /** \brief Print something about branch - only if log level high
  */
  virtual void print(bool normalBranch) {};

  /** \brief Return true if branch should fix variables
  */
  virtual bool boundBranch() const 
  {return true;};

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
  {return variable_;};

  /** Get the state of the branching object
  
    Returns a code indicating the active arm of the branching object.
    The precise meaning is defined in the derived class.

    \sa #way_
  */
  inline int way() const
  {return way_;};

  /** Set the state of the branching object.

    See #way()
  */
  inline void way(int way)
  {way_=way;};

  /// Current value
  inline double value() const
  {return value_;};
  
  /// Return model
  inline CbcModel * model() const
  {return  model_;};

protected:

  /// The model that owns this branching object
  CbcModel * model_;

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

  /// Current value
  double value_;

  /** Number of arms remaining to be evaluated

    \todo Compare with CbcNodeInfo::numberBranchesLeft_, and check for
	  redundancy.
  */
  int numberBranchesLeft_;

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
  virtual int whichMethod() {return 2;};

  /** Saves a clone of current branching object.  Can be used to update
      information on object causing branch - after branch */
  virtual void saveBranchingObject(CbcBranchingObject * object) {};
  /** Pass in information on branch just done.
      assumes object can get information from solver */
  virtual void updateInformation(OsiSolverInterface * solver) {};

protected:
  
  // Clone of branching object
  CbcBranchingObject * object_;
private:
  /// Assignment is illegal
  CbcBranchDecision & operator=(const CbcBranchDecision& rhs);
  
};

#endif
