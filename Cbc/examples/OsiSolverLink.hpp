// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef OsiSolverLink_H
#define OsiSolverLink_H

#include "OsiClpSolverInterface.hpp"
#include "CoinModel.hpp"
class OsiModel;
class CoinPackedMatrix;
class OsiLinkedBound;
class OsiObject;
//#############################################################################

/**
   
This is to allow the user to replace initialSolve and resolve
This version changes coefficients
*/

class OsiSolverLink : public OsiClpSolverInterface {
  
public:
  //---------------------------------------------------------------------------
  /**@name Solve methods */
  //@{
  /// Solve initial LP relaxation 
  virtual void initialSolve();
  
  /// Resolve an LP relaxation after problem modification
  virtual void resolve();

  /**
     Problem specific 
     Returns -1 if node fathomed and no solution
              0 if did nothing
	      1 if node fathomed and solution
     allFixed is true if all LinkedBound variables are fixed
  */
  virtual int fathom(bool allFixed) {return 0;};
  //@}
  
  
  /**@name Constructors and destructors */
  //@{
  /// Default Constructor
  OsiSolverLink ();
  
  /** This creates from a coinModel object 

      if errors.then number of sets is -1
      
      This creates linked ordered sets information.  It assumes -

      for product terms syntax is yy*f(zz)
      also just f(zz) is allowed 
      and even a constant

      modelObject not const as may be changed as part of process.
  */
  OsiSolverLink(  CoinModel & modelObject);
  // real constructor for gdb
  void gdb(  CoinModel & modelObject);
  /// Clone
  virtual OsiSolverInterface * clone(bool copyData=true) const;
  
  /// Copy constructor 
  OsiSolverLink (const OsiSolverLink &);
  
  /// Assignment operator 
  OsiSolverLink & operator=(const OsiSolverLink& rhs);
  
  /// Destructor 
  virtual ~OsiSolverLink ();
  
  //@}
  
  
  /**@name Sets and Gets */
  //@{
  /// Add a bound modifier 
  void addBoundModifier(bool upperBoundAffected, bool useUpperBound, int whichVariable, int whichVariableAffected, 
			double multiplier=1.0);
  /// Update coefficients
  void updateCoefficients(ClpSimplex * solver, CoinPackedMatrix * matrix);
  /// Set special options
  inline void setSpecialOptions2(int value)
  { specialOptions2_=value;};
  /// Get special options
  inline int specialOptions2() const
  { return specialOptions2_;};
  /** Clean copy of matrix
      So we can add rows
  */
  CoinPackedMatrix * cleanMatrix() const
  { return matrix_;};
  //@}
  
  //---------------------------------------------------------------------------
  
protected:
  
  
  /**@name functions */
  //@{
  /// Do real work of initialize
  void initialize(ClpSimplex * & solver, OsiObject ** & object) const;
  /// Do real work of delete
  void gutsOfDestructor(bool justNullify=false);
  /// Do real work of copy
  void gutsOfCopy(const OsiSolverLink & rhs) ;
  //@}
  
  /**@name Private member data */
  //@{
  /** Clean copy of matrix
      Marked coefficients will be multiplied by L or U
  */
  CoinPackedMatrix * matrix_;
  /// Model in CoinModel format
  CoinModel coinModel_;
  /// Number of variables in tightening phase
  int numberVariables_;
  /// Information
  OsiLinkedBound * info_;
  /**
     0 bit (1) - don't do mini B&B
  */
  int specialOptions2_;
  //@}
};
/**
   List of bounds which depend on other bounds
*/

class OsiLinkedBound {
  
public:
  //---------------------------------------------------------------------------
  /**@name Action methods */
  //@{
  /// Update other bounds
  void updateBounds(ClpSimplex * solver);
  //@}
  
  
  /**@name Constructors and destructors */
  //@{
  /// Default Constructor
  OsiLinkedBound ();
  /// Useful Constructor
  OsiLinkedBound(OsiSolverInterface * model, int variable,
		 int numberAffected, const int * positionL, 
		 const int * positionU, const double * multiplier);
  
  /// Copy constructor 
  OsiLinkedBound (const OsiLinkedBound &);
  
  /// Assignment operator 
  OsiLinkedBound & operator=(const OsiLinkedBound& rhs);
  
  /// Destructor 
  ~OsiLinkedBound ();
  
  //@}
  
  /**@name Sets and Gets */
  //@{
  /// Get variable
  inline int variable() const
  { return variable_;};
  /// Add a bound modifier
  void addBoundModifier(bool upperBoundAffected, bool useUpperBound, int whichVariable, 
			double multiplier=1.0);
  //@}
  
private:
  typedef struct {
    /*
      0 - LB of variable affected
      1 - UB of variable affected
      2 - element in position (affected) affected
    */
    unsigned int affect:2;
    unsigned int ubUsed:1; // nonzero if UB of this variable is used
    /* 
       0 - use x*multiplier
       1 - use multiplier/x
       2 - if UB use min of current upper and x*multiplier, if LB use max of current lower and x*multiplier
    */
    unsigned int type:4; // type of computation
    unsigned int affected:25; // variable or element affected
    float multiplier; // to use in computation 
  } boundElementAction;
  
  /**@name Private member data */
  //@{
  /// Pointer back to model 
  OsiSolverInterface * model_;
  /// Variable
  int variable_;
  /// Number of variables/elements affected
  int numberAffected_;
  /// Maximum number of variables/elements affected
  int maximumAffected_;
  /// Actions
  boundElementAction * affected_;
  //@}
};
#include "CbcHeuristic.hpp"
/** heuristic - just picks up any good solution
 */

class CbcHeuristicDynamic3 : public CbcHeuristic {
public:
  
  // Default Constructor 
  CbcHeuristicDynamic3 ();
  
  /* Constructor with model
   */
  CbcHeuristicDynamic3 (CbcModel & model);
  
  // Copy constructor 
  CbcHeuristicDynamic3 ( const CbcHeuristicDynamic3 &);
  
  // Destructor 
  ~CbcHeuristicDynamic3 ();
  
  /// Clone
  virtual CbcHeuristic * clone() const;
  
  /// update model
  virtual void setModel(CbcModel * model);
  
  /** returns 0 if no solution, 1 if valid solution.
      Sets solution values if good, sets objective value (only if good)
      We leave all variables which are at one at this node of the
      tree to that value and will
      initially set all others to zero.  We then sort all variables in order of their cost
      divided by the number of entries in rows which are not yet covered.  We randomize that
      value a bit so that ties will be broken in different ways on different runs of the heuristic.
      We then choose the best one and set it to one and repeat the exercise.  
      
  */
  virtual int solution(double & objectiveValue,
		       double * newSolution);
  /// Resets stuff if model changes
  virtual void resetModel(CbcModel * model);
  
protected:
private:
  /// Illegal Assignment operator 
  CbcHeuristicDynamic3 & operator=(const CbcHeuristicDynamic3& rhs);
};

#endif
