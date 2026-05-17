// Last edit: 4/20/07
//
// Name:     CglRedSplit.hpp
// Author:   Francois Margot
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
//           email: fmargot@andrew.cmu.edu
// Date:     2/6/05
//
//-----------------------------------------------------------------------------
// Copyright (C) 2005, Francois Margot and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglRedSplit_H
#define CglRedSplit_H

#include "CglCutGenerator.hpp"
#include "CglRedSplitParam.hpp"

/** Gomory Reduce-and-Split Cut Generator Class; See method generateCuts().
    Based on the paper by K. Anderson, G. Cornuejols, Yanjun Li, 
    "Reduce-and-Split Cuts: Improving the Performance of Mixed Integer 
    Gomory Cuts", Management Science 51 (2005). */

class CGLLIB_EXPORT CglRedSplit : public CglCutGenerator {

  friend CGLLIB_EXPORT void CglRedSplitUnitTest(const OsiSolverInterface * siP,
				  const std::string mpdDir);
public:
  /**@name generateCuts */
  //@{
  /** Generate Reduce-and-Split Mixed Integer Gomory cuts 
      for the model of the solver interface si.

      Insert the generated cuts into OsiCuts cs.

      Warning: This generator currently works only with the Lp solvers Clp or 
      Cplex9.0 or higher. It requires access to the optimal tableau and 
      optimal basis inverse and makes assumptions on the way slack variables 
      are added by the solver. The Osi implementations for Clp and Cplex 
      verify these assumptions.

      When calling the generator, the solver interface si 
      must contain an optimized
      problem and information related to the optimal basis must be available 
      through the OsiSolverInterface methods (si->optimalBasisIsAvailable()
      must return 'true'). It is also essential that the integrality of
      structural variable i can be obtained using si->isInteger(i).

      Reduce-and-Split cuts are variants of Gomory cuts: Starting from
      the current optimal tableau, linear combinations of the rows of 
      the current optimal simplex tableau are used for generating Gomory
      cuts. The choice of the linear combinations is driven by the objective 
      of reducing the coefficients of the non basic continuous variables
      in the resulting row.
      Note that this generator might not be able to generate cuts for some 
      solutions violating integrality constraints. 

  */
  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());

  /// Return true if needs optimal basis to do cuts (will return true)
  virtual bool needsOptimalBasis() const;
  //@}
  
  
  /**@name Public Methods */
  //@{

  // Set the parameters to the values of the given CglRedSplitParam object.
  void setParam(const CglRedSplitParam &source); 
  // Return the CglRedSplitParam object of the generator. 
  inline CglRedSplitParam getParam() const {return param;}

  // Compute entries of low_is_lub and up_is_lub.
  void compute_is_lub();

  // Compute entries of is_integer.
  void compute_is_integer();

  /// Set given_optsol to the given optimal solution given_sol.
  /// If given_optsol is set using this method, 
  /// the code will stop as soon as
  /// a generated cut is violated by the given solution; exclusively 
  /// for debugging purposes.
  void set_given_optsol(const double *given_sol, const int card_sol);

  /// Print some of the data members  
  void print() const;

  /// Print the current simplex tableau  
  void printOptTab(OsiSolverInterface *solver) const;
  
  //@}

 /**@name Public Methods (soon to be obsolete)*/
  //@{
  //************************************************************
  // TO BE REMOVED
  /** Set limit, the maximum number of non zero coefficients in generated cut;
      Default: 50 */
  void setLimit(int limit);
  /** Get value of limit */
  int getLimit() const;

  /** Set away, the minimum distance from being integer used for selecting 
      rows for cut generation;  all rows whose pivot variable should be 
      integer but is more than away from integrality will be selected; 
      Default: 0.05 */
  void setAway(double value);
  /// Get value of away
  double getAway() const;
 /** Set the value of LUB, value considered large for the absolute value of
      a lower or upper bound on a variable;
      Default: 1000 */
  void setLUB(double value);
  /** Get the value of LUB */
  double getLUB() const;

  /** Set the value of EPS, epsilon for double computations;
      Default: 1e-7 */
  void setEPS(double value);
  /** Get the value of EPS */
  double getEPS() const;

  /** Set the value of EPS_COEFF, epsilon for values of coefficients;
      Default: 1e-8 */
  void setEPS_COEFF(double value);
  /** Get the value of EPS_COEFF */
  double getEPS_COEFF() const;

  /** Set the value of EPS_COEFF_LUB, epsilon for values of coefficients for 
      variables with absolute value of lower or upper bound larger than LUB;
      Default: 1e-13 */
  void setEPS_COEFF_LUB(double value);
  /** Get the value of EPS_COEFF_LUB */
  double getEPS_COEFF_LUB() const;

  /** Set the value of EPS_RELAX, value used for relaxing the right hand side
      of each generated cut;
      Default: 1e-8 */
  void setEPS_RELAX(double value);
  /** Get the value of EPS_RELAX */
  double getEPS_RELAX() const;

  /** Set the value of normIsZero, the threshold for considering a norm to be 
      0; Default: 1e-5 */
  void setNormIsZero(double value);
  /** Get the value of normIsZero */
  double getNormIsZero() const;

  /** Set the value of minReduc, threshold for relative norm improvement for
   performing  a reduction; Default: 0.05 */
  void setMinReduc(double value);
  /// Get the value of minReduc
  double getMinReduc() const;

  /** Set the maximum allowed value for (mTab * mTab * std::max(mTab, nTab)) where 
      mTab is the number of rows used in the combinations and nTab is the 
      number of continuous non basic variables. The work of the generator is 
      proportional to (mTab * mTab * std::max(mTab, nTab)). Reducing the value of 
      maxTab makes the generator faster, but weaker. Default: 1e7. */
  void setMaxTab(double value);
  /// Get the value of maxTab
  double getMaxTab() const;
  // END TO BE REMOVED
  //************************************************************

  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglRedSplit();

  /// Constructor with specified parameters 
  CglRedSplit(const CglRedSplitParam &RS_param);
 
  /// Copy constructor 
  CglRedSplit (const CglRedSplit &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglRedSplit &
    operator=(
    const CglRedSplit& rhs);
  
  /// Destructor 
  virtual
    ~CglRedSplit ();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  //@}
    
private:
  
  // Private member methods

/**@name Private member methods */

  //@{

  // Method generating the cuts after all CglRedSplit members are properly set.
  void generateCuts(OsiCuts & cs);

  /// Compute the fractional part of value, allowing for small error.
  inline double rs_above_integer(double value); 

  /// Perform row r1 of pi := row r1 of pi - step * row r2 of pi.
  void update_pi_mat(int r1, int r2, int step);

  /// Perform row r1 of tab := row r1 of tab - step * row r2 of tab.
  void update_redTab(int r1, int r2, int step);

  /// Find optimal integer step for changing row r1 by adding to it a 
  /// multiple of another row r2.
  void find_step(int r1, int r2, int *step, 
		 double *reduc, double *norm);

  /// Test if an ordered pair of rows yields a reduction. Perform the
  /// reduction if it is acceptable.
  int test_pair(int r1, int r2, double *norm);

  /// Reduce rows of contNonBasicTab.
  void reduce_contNonBasicTab();

  /// Generate a row of the current LP tableau.
  void generate_row(int index_row, double *row);

  /// Generate a mixed integer Chvatal-Gomory cut, when all non basic 
  /// variables are non negative and at their lower bound.
  int generate_cgcut(double *row, double *rhs);

  /// Generate a mixed integer Chvatal-Gomory cut, when all non basic 
  /// variables are non negative and at their lower bound (different formula)
  int generate_cgcut_2(int basic_ind, double *row, double *rhs);

  /// Use multiples of the initial inequalities to cancel out the coefficients
  /// of the slack variables.
  void eliminate_slacks(double *row, 
			const double *elements, 
			const CoinBigIndex *start,
			const int *indices,
			const int *rowLength,
			const double *rhs, double *rowrhs);

  /// Change the sign of the coefficients of the continuous non basic
  /// variables at their upper bound.
  void flip(double *row);

  /// Change the sign of the coefficients of the continuous non basic
  /// variables at their upper bound and do the translations restoring
  /// the original bounds. Modify the right hand side
  /// accordingly.
  void unflip(double *row, double *rowrhs, double *slack_val);

  /// Return the scale factor for the row. 
  /// Compute max_coeff: maximum absolute value of the coefficients.
  /// Compute min_coeff: minimum absolute value of the coefficients
  /// larger than EPS_COEFF.
  /// Return -1 if max_coeff < EPS_COEFF or if max_coeff/min_coeff > MAXDYN
  /// or MAXDYN_LUB (depending if the row has a non zero coeff. for a variable
  /// with large lower/upper bound) */.
  double row_scale_factor(double *row);

  /// Generate the packed cut from the row representation.
  int generate_packed_row(const double *xlp, double *row,
			  int *rowind, double *rowelem, 
			  int *card_row, double & rhs);

  /// Check that the generated cuts do not cut a given optimal solution.
  void check_optsol(const int calling_place,
		    const double *xlp, const double *slack_val,
		    const int do_flip);

  /// Check that the generated cuts do not cut a given optimal solution.
  void check_optsol(const int calling_place,
		    const double *xlp, const double *slack_val,
		    const double *ck_row, const double ck_rhs, 
		    const int cut_number, const int do_flip);

  // Check that two vectors are different.
  bool rs_are_different_vectors(const int *vect1, 
				const int *vect2,
				const int dim);

  // Check that two vectors are different.
  bool rs_are_different_vectors(const double *vect1, 
				const double *vect2,
				const int dim);

  // Check that two matrices are different.
  bool rs_are_different_matrices(const CoinPackedMatrix *mat1, 
				 const CoinPackedMatrix *mat2,
				 const int nmaj,
				 const int nmin);
  //@}

  
  // Private member data

/**@name Private member data */

  //@{

  /// Object with CglRedSplitParam members. 
  CglRedSplitParam param;

  /// Number of rows ( = number of slack variables) in the current LP.
  int nrow; 

  /// Number of structural variables in the current LP.
  int ncol;

  /// Lower bounds for structural variables
  const double *colLower;

  /// Upper bounds for structural variables
  const double *colUpper;
  
  /// Lower bounds for constraints
  const double *rowLower;

  /// Upper bounds for constraints
  const double *rowUpper;

  /// Righ hand side for constraints (upper bound for ranged constraints).
  const double *rowRhs;

  /// Number of integer basic structural variables that are fractional in the
  /// current lp solution (at least param.away_ from being integer).  
  int card_intBasicVar_frac;

  /// Number of integer non basic structural variables in the
  /// current lp solution.  
  int card_intNonBasicVar; 

  /// Number of continuous non basic variables (structural or slack) in the
  /// current lp solution.  
  int card_contNonBasicVar;

  /// Number of non basic variables (structural or slack) at their
  /// upper bound in the current lp solution.
  int card_nonBasicAtUpper; 

  /// Number of non basic variables (structural or slack) at their
  /// lower bound in the current lp solution.
  int card_nonBasicAtLower;

  /// Characteristic vector for integer basic structural variables
  /// with non integer value in the current lp solution.
  int *cv_intBasicVar_frac;  

  /// List of integer structural basic variables 
  /// (in order of pivot in selected rows for cut generation).
  int *intBasicVar_frac;

  /// List of integer structural non basic variables.
  int *intNonBasicVar; 

  /// List of continuous non basic variables (structural or slack). 
  // slacks are considered continuous (no harm if this is not the case).
  int *contNonBasicVar;

  /// List of non basic variables (structural or slack) at their 
  /// upper bound. 
  int *nonBasicAtUpper;

  /// List of non basic variables (structural or slack) at their lower
  /// bound.
  int *nonBasicAtLower;

  /// Number of rows in the reduced tableau (= card_intBasicVar_frac).
  int mTab;

  /// Number of columns in the reduced tableau (= card_contNonBasicVar)
  int nTab;

  /// Tableau of multipliers used to alter the rows used in generation.
  /// Dimensions: mTab by mTab. Initially, pi_mat is the identity matrix.
  int **pi_mat;

  /// Current tableau for continuous non basic variables (structural or slack).
  /// Only rows used for generation.
  /// Dimensions: mTab by nTab.
  double **contNonBasicTab;

  /// Current tableau for integer non basic structural variables.
  /// Only rows used for generation.
  // Dimensions: mTab by card_intNonBasicVar.
  double **intNonBasicTab;

  /// Right hand side of the tableau.
  /// Only rows used for generation.
  double *rhsTab ;

  /// Given optimal solution that should not be cut; only for debug. 
  const double *given_optsol;

  /// Number of entries in given_optsol.
  int card_given_optsol;

  /// Characteristic vectors of structural integer variables or continuous
  /// variables currently fixed to integer values. 
  int *is_integer;

  /// Characteristic vector of the structural variables whose lower bound 
  /// in absolute value is larger than LUB. 
  int *low_is_lub;

  /// Characteristic vector of the structural variables whose upper bound 
  /// in absolute value is larger than LUB. 
  int *up_is_lub;

  /// Pointer on solver. Reset by each call to generateCuts().
  OsiSolverInterface *solver;

  /// Pointer on point to separate. Reset by each call to generateCuts().
  const double *xlp;

  /// Pointer on row activity. Reset by each call to generateCuts().
  const double *rowActivity;

  /// Pointer on column type. Reset by each call to generateCuts().
  const char *colType;

  /// Pointer on matrix of coefficient ordered by rows. 
  /// Reset by each call to generateCuts().
  const CoinPackedMatrix *byRow;

  //@}
};

//#############################################################################
/** A function that tests the methods in the CglRedSplit class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
CGLLIB_EXPORT
void CglRedSplitUnitTest(const OsiSolverInterface * siP,
			 const std::string mpdDir );


#endif
