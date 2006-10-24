// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef CbcSolverLink_H
#define CbcSolverLink_H

#include "OsiClpSolverInterface.hpp"
#include "CoinModel.hpp"
class CbcModel;
class CoinPackedMatrix;
class CbcLinkedBound;
class CbcObject;
//#############################################################################

/**
   
This is to allow the user to replace initialSolve and resolve
This version changes coefficients
*/

class CbcSolverLink : public OsiClpSolverInterface {
  
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
  CbcSolverLink ();
  
  /** This creates from a coinModel object 

      if errors.then number of sets is -1
      
      This creates linked ordered sets information.  It assumes -

      for product terms syntax is yy*f(zz)
      also just f(zz) is allowed 
      and even a constant

      modelObject not const as may be changed as part of process.
  */
  CbcSolverLink(  CoinModel & modelObject);
  // real constructor for gdb
  void gdb(  CoinModel & modelObject);
  /// Clone
  virtual OsiSolverInterface * clone(bool copyData=true) const;
  
  /// Copy constructor 
  CbcSolverLink (const CbcSolverLink &);
  
  /// Assignment operator 
  CbcSolverLink & operator=(const CbcSolverLink& rhs);
  
  /// Destructor 
  virtual ~CbcSolverLink ();
  
  //@}
  
  
  /**@name Sets and Gets */
  //@{
  /// Setup 
  void initialize(CbcModel * model, const CoinPackedMatrix * matrix,
		  int numberVariables, const int * which,
		  const int *starts, const int * positionL, 
		  const int * positionU);
  /** This creates a solver instance and information and passes to CbcModel
      Granularity can be given for each set
      Negative values mean S1 type grid (using abs(value).
  */
  void initialize ( CbcModel * model, bool top=true);
  /// revert to previous solver
  void revert();
  /// Add a bound modifier 
  void addBoundModifier(bool upperBoundAffected, bool useUpperBound, int whichVariable, int whichVariableAffected, 
			double multiplier=1.0);
  /// Objective value of best solution found internally
  inline double bestObjectiveValue() const
  { return bestObjectiveValue_;};
  /// Best solution found internally
  const double * bestSolution() const
  { return bestSolution_;};
  /// Set special options
  inline void setSpecialOptions2(int value)
  { specialOptions2_=value;};
  /// Get special options
  inline int specialOptions2() const
  { return specialOptions2_;};
  /// Number of Linked Ordered Sets
  inline int numberSets() const
  { return numberSets_;};
  /// Sets granularity
  void setGranularity(double top, double bottom, 
		      const double * topG = NULL, const double * bottomG = NULL);
  //@}
  
  //---------------------------------------------------------------------------
  
protected:
  
  
  /**@name functions */
  //@{
  /// Do real work of initialize
  void initialize(ClpSimplex * & solver, CbcObject ** & object) const;
  /// Do real work of delete
  void gutsOfDestructor(bool justNullify=false);
  /// Do real work of copy
  void gutsOfCopy(const CbcSolverLink & rhs) ;
  //@}
  
  /**@name Private member data */
  //@{
  /// Pointer back to model
  CbcModel * model_;
  /** Clean copy of matrix
      Marked coefficients will be multiplied by L or U
  */
  CoinPackedMatrix * matrix_;
  /// Number of variables in tightening phase
  int numberVariables_;
  /// Information
  CbcLinkedBound * info_;
  /// Objective value of best solution found internally
  double bestObjectiveValue_;
  /// Best solution found internally
  double * bestSolution_;
  /**
     0 bit (1) - don't do mini B&B
  */
  int specialOptions2_;
  /// Number of Sets
  int numberSets_;
  /// Number of nonlinear terms
  int numberNonlinear_;
  /// Start of each linked ordered set
  int * startSet_;
  /// Which row is nonlinear definition row
  int * definitionRow_;
  /// Which row is linear definition row (-1 for just f(z), -2 for numeric)
  int * definitionRowL_;
  /// Nonlinear column
  int * functionColumn_;
  /// Linear column
  int * definitionColumnL_;
  /// Which row is usage row
  int * usageRow_;
  /// Convexity row
  int * convexityRow_;
  /// Copy of truncated model
  ClpSimplex * simplexModel_;
  /// Copy of original model
  ClpSimplex * savedModel_;
  /// Which variables of original were nonlinear
  int * which_;
  /// Copy of CoinModel
  CoinModel coinModel_;
  /// Granularity at top of tree
  double * granularityTop_;
  /// Granularity at bottom of tree (when all integers gone)
  double * granularityBottom_;
  //@}
};
/**
   List of bounds which depend on other bounds
*/

class CbcLinkedBound {
  
public:
  //---------------------------------------------------------------------------
  /**@name Action methods */
  //@{
  /// Update other bounds
  void updateBounds(ClpSimplex * solver);
  /// Update coefficients
  void updateCoefficients(ClpSimplex * solver, CoinPackedMatrix * matrix);
  /// Fix variables in LOS
  void fixLOS(ClpSimplex * solver);
  /// Update this variables bounds from reference row
  void updateBoundsForThis(ClpSimplex * solver);
  //@}
  
  
  /**@name Constructors and destructors */
  //@{
  /// Default Constructor
  CbcLinkedBound ();
  /// Useful Constructor
  CbcLinkedBound(CbcModel * model, int variable,
		 int numberAffected, const int * positionL, 
		 const int * positionU, const double * multiplier);
  
  /// Copy constructor 
  CbcLinkedBound (const CbcLinkedBound &);
  
  /// Assignment operator 
  CbcLinkedBound & operator=(const CbcLinkedBound& rhs);
  
  /// Destructor 
  ~CbcLinkedBound ();
  
  //@}
  
  /**@name Sets and Gets */
  //@{
  /// Get variable
  inline int variable() const
  { return variable_;};
  /// Add a bound modifier
  void addBoundModifier(bool upperBoundAffected, bool useUpperBound, int whichVariable, 
			double multiplier=1.0);
  /// Add a bound modifier
  void addCoefficientModifier( bool useUpperBound, int position, 
			double multiplier=1.0);
  /// Get useObject
  inline int useObject() const
  { return useObject_;};
  /// Set useObject
  inline void setUseObject(int value)
  { useObject_ = value;};
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
  CbcModel * model_;
  /// Variable
  int variable_;
  /// Number of variables/elements affected
  int numberAffected_;
  /// Maximum number of variables/elements affected
  int maximumAffected_;
  /// Actions
  boundElementAction * affected_;
  /// If >=0 object number which can be used to determine bounds on variable
  int useObject_;
  /// Number of sets for variable fixing
  int numberSets_;
  /// Which objects whose reference row entry can be used
  int * whichSets_;
  //@}
};
#include "CbcHeuristic.hpp"
/** heuristic - just picks up any good solution
 */

class CbcHeuristicDynamic2 : public CbcHeuristic {
public:
  
  // Default Constructor 
  CbcHeuristicDynamic2 ();
  
  /* Constructor with model
   */
  CbcHeuristicDynamic2 (CbcModel & model);
  
  // Copy constructor 
  CbcHeuristicDynamic2 ( const CbcHeuristicDynamic2 &);
  
  // Destructor 
  ~CbcHeuristicDynamic2 ();
  
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
  CbcHeuristicDynamic2 & operator=(const CbcHeuristicDynamic2& rhs);
};

#endif
