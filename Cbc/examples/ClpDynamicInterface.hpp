// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef ClpDynamicInterface_H
#define ClpDynamicInterface_H

#include "OsiClpSolverInterface.hpp"
#include "CoinWarmStartBasisDynamic.hpp"

//#############################################################################

/**

    This is to allow the user to replace initialSolve and resolve
*/
class CbcModel;
class ClpDynamicInterface : public OsiClpSolverInterface {

public:
  //---------------------------------------------------------------------------
  /**@name Solve methods */
  //@{
    /// Solve initial LP relaxation 
    virtual void initialSolve();

    /// Resolve an LP relaxation after problem modification
    virtual void resolve();

  //@}


  /**@name Constructors and destructors */
  //@{
    /// Default Constructor
    ClpDynamicInterface ();
    
    /// Clone
    virtual OsiSolverInterface * clone(bool CopyData=true) const;
    
    /// Copy constructor 
    ClpDynamicInterface (const ClpDynamicInterface &);
    
    /// Assignment operator 
    ClpDynamicInterface & operator=(const ClpDynamicInterface& rhs);
    
    /// Destructor 
    virtual ~ClpDynamicInterface ();

  //@}


  /**@name Sets and Gets */
  //@{
  /** Setup */
  void initialize(int * rowBlock,CbcModel * model);
  /// Model
  inline void setCbcModel(CbcModel * model)
  { cbcModel_=model;};
  /** Warm start
      
  NOTE  artificials are treated as +1 elements so for <= rhs
  artificial will be at lower bound if constraint is tight
  
  This means that Clpsimplex flips artificials as it works
  in terms of row activities
  */
  CoinWarmStartBasisDynamic getBasis(ClpSimplex * model, const int * whichColumns,
                                     int numberCommon) const;
  /** Sets up working basis as a copy of input
      
  NOTE  artificials are treated as +1 elements so for <= rhs
  artificial will be at lower bound if constraint is tight
  
  This means that Clpsimplex flips artificials as it works
  in terms of row activities
  */
  void setBasis( const CoinWarmStartBasisDynamic & basis, ClpSimplex * model);

  /// Objective value of best solution found internally
  inline double bestObjectiveValue() const
  { return bestObjectiveValue_;};
  /// Best solution found internally
  const double * bestSolution() const
  { return bestSolution_;};
  //@}

  //---------------------------------------------------------------------------

private:
  /// Adds proposals to model (only if not there)
  void addProposals(ClpSimplex * model, int number, const int * which,
                    bool addEvenIfFixed);
  /** Creates modelPtr_ solution.
      Returns 1 if not feasible, -1 if integer solution
      0 otherwise */
  int setSolution(ClpSimplex * model);
private:
  
  /**@name Private member data */
  //@{
  // Will have decomposition here
  /// Number of blocks
  int numberBlocks_;
  /// Number of artificials to make problem feasible
  int numberArtificials_;
  /// Static part of model plus convexity rows
  ClpSimplex * staticModel_;
  /// Columns for each block
  ClpSimplex * subProblem_;
  /// Master part of each block
  CoinPackedMatrix * masterRow_;
  /** Proposals
      Rows 0-masterRows-1 are elements
      Next row is convexity row (so gives block)
      Next (==masterRows+convexity) is cost
      Next are weights for columns
  */
  CoinPackedMatrix proposals_;
  /// Number of proposals added to current problem
  int numberProposalsAdded_;
  /// Maximum number of proposals added to current problem (size of list)
  int maximumProposalsAdded_;
  /// List of proposals added
  int * proposalsAdded_;
  /// Backward pointers to original for each sub-problem
  int ** backward_;
  /** If nonzero then says all possible columns generated -
      and gives number */
  int * allGenerated_;
  /// Pointer back to model
  CbcModel * cbcModel_;
  /// Objective value of best solution found internally
  double bestObjectiveValue_;
  /// Best solution found internally
  double * bestSolution_;
  /** Warmstart information to be used in resolves. */
  CoinWarmStartBasisDynamic dynamicBasis_;
  //@}
};
#include "CbcHeuristic.hpp"
/** heuristic - just picks up any good solution
 */

class CbcHeuristicDynamic : public CbcHeuristic {
public:

  // Default Constructor 
  CbcHeuristicDynamic ();

  /* Constructor with model
  */
  CbcHeuristicDynamic (CbcModel & model);
  
  // Copy constructor 
  CbcHeuristicDynamic ( const CbcHeuristicDynamic &);
   
  // Destructor 
  ~CbcHeuristicDynamic ();
  
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
  CbcHeuristicDynamic & operator=(const CbcHeuristicDynamic& rhs);
};

#endif
