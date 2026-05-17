// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglGomory_H
#define CglGomory_H

#include <string>

#include "CglCutGenerator.hpp"

class CoinWarmStartBasis;
/** Gomory Cut Generator Class */
class CGLLIB_EXPORT CglGomory : public CglCutGenerator {
   friend CGLLIB_EXPORT void CglGomoryUnitTest(const OsiSolverInterface * siP,
				  const std::string mpdDir );
 
public:
    
  
  /**@name Generate Cuts */
  //@{
  /** Generate Mixed Integer Gomory cuts for the model of the 
      solver interface, si.

      Insert the generated cuts into OsiCut, cs.

      There is a limit option, which will only generate cuts with
      less than this number of entries.

      We can also only look at 0-1 variables a certain distance
      from integer.
  */
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info = CglTreeInfo());
  /** Generates cuts given matrix and solution etc,
      returns number of cuts generated */
  int generateCuts( const OsiRowCutDebugger * debugger, 
		    OsiCuts & cs,
		    const CoinPackedMatrix & columnCopy,
		    const CoinPackedMatrix & rowCopy,
		    const double * colsol,
		    const double * colLower, const double * colUpper,
		    const double * rowLower, const double * rowUpper,
		    const char * intVar ,
		    const CoinWarmStartBasis* warm,
                    const CglTreeInfo info = CglTreeInfo());
  /** Generates cuts given matrix and solution etc,
      returns number of cuts generated (no row copy passed in) */
  int generateCuts( const OsiRowCutDebugger * debugger, 
		    OsiCuts & cs,
		    const CoinPackedMatrix & columnCopy,
		    const double * colsol,
		    const double * colLower, const double * colUpper,
		    const double * rowLower, const double * rowUpper,
		    const char * intVar ,
		    const CoinWarmStartBasis* warm,
                    const CglTreeInfo info = CglTreeInfo());

  /// Return true if needs optimal basis to do cuts (will return true)
  virtual bool needsOptimalBasis() const { return true; }
  //@}

  /**@name Change way Gomory works */
  //@{
  /// Pass in a copy of original solver (clone it)
  void passInOriginalSolver(OsiSolverInterface * solver);
  /// Set type - 0 normal, 1 add original matrix one, 2 replace
  inline void setGomoryType(int type)
  { gomoryType_=type;}
  /// Return type
  inline int gomoryType() const
  { return gomoryType_;}
  //@}

  /**@name Change limit on how many variables in cut (default 50) */
  //@{
  /// Set
  void setLimit(int limit);
  /// Get
  int getLimit() const;
  /// Set at root (if <normal then use normal)
  void setLimitAtRoot(int limit);
  /// Get at root
  int getLimitAtRoot() const;
  /// Return maximum length of cut in tree
  virtual int maximumLengthOfCutInTree() const;
  //@}

  /**@name Change criterion on which variables to look at.  All ones
   more than "away" away from integrality will be investigated 
  (default 0.05) */
  //@{
  /// Set away
  void setAway(double value);
  /// Get away
  double getAway() const;
  /// Set away at root
  void setAwayAtRoot(double value);
  /// Get away at root
  double getAwayAtRoot() const;
  //@}

  /**@name Change criterion on which the cut id relaxed if the code
           thinks the factorization has inaccuracies.  The relaxation to
	   RHS is smallest of -
	   1) 1.0e-4
	   2) conditionNumberMultiplier * condition number of factorization
	   3) largestFactorMultiplier * largest (dual*element) forming tableau
	      row
  */
  //@{
  /// Set ConditionNumberMultiplier
  void setConditionNumberMultiplier(double value);
  /// Get ConditionNumberMultiplier
  double getConditionNumberMultiplier() const;
  /// Set LargestFactorMultiplier
  void setLargestFactorMultiplier(double value);
  /// Get LargestFactorMultiplier
  double getLargestFactorMultiplier() const;
  //@}

  /**@name change factorization */
  //@{
   /// Set/unset alternative factorization
   inline void useAlternativeFactorization(bool yes=true)
   { alternateFactorization_= (yes) ? 1 : 0;} 
   /// Get whether alternative factorization being used
   inline bool alternativeFactorization() const
   { return (alternateFactorization_!=0);} 
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglGomory ();
 
  /// Copy constructor 
  CglGomory (
    const CglGomory &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglGomory &
    operator=(
    const CglGomory& rhs);
  
  /// Destructor 
  virtual
    ~CglGomory ();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  /// This can be used to refresh any inforamtion
  virtual void refreshSolver(OsiSolverInterface * solver);
  //@}
      
private:
  
 // Private member methods

  // Private member data

  /**@name Private member data */
  //@{
  /// Only investigate if more than this away from integrality
  double away_;
  /// Only investigate if more than this away from integrality (at root)
  double awayAtRoot_;
  /// Multiplier for conditionNumber cut relaxation
  double conditionNumberMultiplier_;
  /// Multiplier for largest factor cut relaxation
  double largestFactorMultiplier_;
  /// Limit - only generate if fewer than this in cut
  int limit_;
  /// Limit - only generate if fewer than this in cut (at root)
  int limitAtRoot_;
  /// Dynamic limit in tree
  int dynamicLimitInTree_;
  /// Number of times stalled
  int numberTimesStalled_;
  /// nonzero to use alternative factorization
  int alternateFactorization_;
  /// Type - 0 normal, 1 add original matrix one, 2 replace
  int gomoryType_; // note could add in cutoff as constraint
  //@}
};

//#############################################################################
/** A function that tests the methods in the CglGomory class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
CGLLIB_EXPORT
void CglGomoryUnitTest(const OsiSolverInterface * siP,
			const std::string mpdDir );
  
#endif
