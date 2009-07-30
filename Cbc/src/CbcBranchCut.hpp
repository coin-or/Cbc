/* $Id$ */
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcBranchCut_H
#define CbcBranchCut_H

#include "CbcBranchBase.hpp"
#include "OsiRowCut.hpp"
#include "CoinPackedMatrix.hpp"

/** Define a cut branching class.
    At present empty - all stuff in descendants
*/

class CbcBranchCut : public CbcObject {

public:

  // Default Constructor 
  CbcBranchCut ();

  /** In to maintain normal methods
  */ 
  CbcBranchCut (CbcModel * model);
  // Copy constructor 
  CbcBranchCut ( const CbcBranchCut &);
   
  /// Clone
  virtual CbcObject * clone() const;

  // Assignment operator 
  CbcBranchCut & operator=( const CbcBranchCut& rhs);

  // Destructor 
  ~CbcBranchCut ();
  
  /// Infeasibility 
  virtual double infeasibility(const OsiBranchingInformation * info,
			       int &preferredWay) const;

  using CbcObject::feasibleRegion ;
  /** Set bounds to contain the current solution.

    More precisely, for the variable associated with this object, take the
    value given in the current solution, force it within the current bounds
    if required, then set the bounds to fix the variable at the integer
    nearest the solution value.

    At present this will do nothing
  */
  virtual void feasibleRegion();

  /** \brief Return true if branch created by object should fix variables
  */
  virtual bool boundBranch() const ;

  /// Creates a branching object
  virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver,const OsiBranchingInformation * info, int way) ;

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

    At present this does nothing
  */
  virtual CbcBranchingObject * preferredNewFeasible() const;
  
  /** \brief Given a valid solution (with reduced costs, etc.),
      return a branching object which would give a new feasible
      point in a bad direction.

    As for preferredNewFeasible(), but the preferred branching object will
    force movement in a direction that degrades the objective.

    At present this does nothing
  */
  virtual CbcBranchingObject * notPreferredNewFeasible() const ;
  
  using CbcObject::resetBounds ;
  /** Reset original upper and lower bound values from the solver.
  
    Handy for updating bounds held in this object after bounds held in the
    solver have been tightened.
   */
  virtual void resetBounds();
  

protected:
  /// data

};


/** Cut branching object

  This object can specify a two-way branch in terms of two cuts
*/

class CbcCutBranchingObject : public CbcBranchingObject {

public:

  /// Default constructor 
  CbcCutBranchingObject ();

  /** Create a cut branching object

      Cut down will applied on way=-1, up on way==1
      Assumed down will be first so way_ set to -1
  */
  CbcCutBranchingObject (CbcModel * model, OsiRowCut & down, OsiRowCut &up, bool canFix);
  
  /// Copy constructor 
  CbcCutBranchingObject ( const CbcCutBranchingObject &);
   
  /// Assignment operator 
  CbcCutBranchingObject & operator= (const CbcCutBranchingObject& rhs);

  /// Clone
  virtual CbcBranchingObject * clone() const;

  /// Destructor 
  virtual ~CbcCutBranchingObject ();
  
  using CbcBranchingObject::branch ;
  /** \brief Sets the bounds for variables or adds a cut depending on the
             current arm of the branch and advances the object state to the next arm.
	     Returns change in guessed objective on next branch
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

  /** \brief Return true if branch should fix variables
  */
  virtual bool boundBranch() const;

  /** Return the type (an integer identifier) of \c this */
  virtual int type() const { return 200; }

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

protected:
  /// Cut for the down arm (way_ = -1)
  OsiRowCut down_;
  /// Cut for the up arm (way_ = 1)
  OsiRowCut up_;
  /// True if one way can fix variables
  bool canFix_;
};


/** Define a branch class that branches so that one way variables are fixed
    while the other way cuts off that solution.
    a) On reduced cost
    b) When enough ==1 or <=1 rows have been satisfied (not fixed - satisfied)
*/


class CbcBranchToFixLots : public CbcBranchCut {

public:

  // Default Constructor 
  CbcBranchToFixLots ();

  /** Useful constructor - passed reduced cost tolerance and fraction we would like fixed.
      Also depth level to do at.
      Also passed number of 1 rows which when clean triggers fix
      Always does if all 1 rows cleaned up and number>0 or if fraction columns reached
      Also whether to create branch if can't reach fraction.
  */ 
  CbcBranchToFixLots (CbcModel * model, double djTolerance,
		      double fractionFixed, int depth,
		      int numberClean=0,
		      const char * mark=NULL,
		      bool alwaysCreate=false);
  
  // Copy constructor 
  CbcBranchToFixLots ( const CbcBranchToFixLots &);
   
  /// Clone
  virtual CbcObject * clone() const;

  // Assignment operator 
  CbcBranchToFixLots & operator=( const CbcBranchToFixLots& rhs);

  // Destructor 
  ~CbcBranchToFixLots ();

  /** Does a lot of the work,
      Returns 0 if no good, 1 if dj, 2 if clean, 3 if both
  */
  int shallWe() const;

  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiBranchingInformation * info,
			       int &preferredWay) const;
  /** \brief Return true if object can take part in normal heuristics
  */
  virtual bool canDoHeuristics() const 
  {return true;}

  /// Creates a branching object
  virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver,const OsiBranchingInformation * info, int way) ;
  /// Redoes data when sequence numbers change
  virtual void redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns);


protected:
  /// data

  /// Reduced cost tolerance i.e. dj has to be >= this before fixed
  double djTolerance_;
  /// We only need to make sure this fraction fixed
  double fractionFixed_;
  /// Never fix ones marked here
  char * mark_;
  /// Matrix by row
  CoinPackedMatrix matrixByRow_; 
  /// Do if depth multiple of this
  int depth_;
  /// number of ==1 rows which need to be clean
  int numberClean_;
  /// If true then always create branch
  bool alwaysCreate_;
};

/** Define a branch class that branches so that it is only satsified if all
    members have different values
    So cut is x <= y-1 or x >= y+1
*/


class CbcBranchAllDifferent : public CbcBranchCut {

public:

  // Default Constructor 
  CbcBranchAllDifferent ();

  /** Useful constructor - passed set of integer variables which must all be different
  */ 
  CbcBranchAllDifferent (CbcModel * model, int number,const int * which);
  
  // Copy constructor 
  CbcBranchAllDifferent ( const CbcBranchAllDifferent &);
   
  /// Clone
  virtual CbcObject * clone() const;

  // Assignment operator 
  CbcBranchAllDifferent & operator=( const CbcBranchAllDifferent& rhs);

  // Destructor 
  ~CbcBranchAllDifferent ();

  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiBranchingInformation * info,
			       int &preferredWay) const;

  /// Creates a branching object
  virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver,const OsiBranchingInformation * info, int way) ;


protected:
  /// data

  /// Number of entries
  int numberInSet_;
  /// Which variables
  int * which_;
};
#endif
