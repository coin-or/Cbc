// Last edit: 04/03/10
//
// Name:     CglRedSplit2.hpp
// Author:   Giacomo Nannicini
//           Singapore University of Technology and Design
//           Singapore
//           email: nannicini@sutd.edu.sg
//           based on CglRedSplit by Francois Margot
// Date:     03/09/09
//-----------------------------------------------------------------------------
// Copyright (C) 2010, Giacomo Nannicini and others.  All Rights Reserved.

#ifndef CglRedSplit2_H
#define CglRedSplit2_H

#include "CglCutGenerator.hpp"
#include "CglRedSplit2Param.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"

/** Reduce-and-Split Cut Generator Class; See method generateCuts().
    Based on the papers "Practical strategies for generating rank-1
    split cuts in mixed-integer linear programming" by G. Cornuejols
    and G. Nannicini, published on Mathematical Programming
    Computation, and "Combining Lift-and-Project and Reduce-and-Split"
    by E. Balas, G. Cornuejols, T. Kis and G. Nannicini, published on
    INFORMS Journal on Computing. Part of this code is based on
    CglRedSplit by F. Margot. */

class CGLLIB_EXPORT CglRedSplit2 : public CglCutGenerator {

  friend CGLLIB_EXPORT void CglRedSplit2UnitTest(const OsiSolverInterface * siP,
				  const std::string mpdDir);
public:
  /**@name generateCuts */
  //@{
  /** Generate Reduce-and-Split Mixed Integer Gomory cuts 
      for the model of the solver interface si.

      Insert the generated cuts into OsiCuts cs.

      This generator currently works only with the Lp solvers Clp or
      Cplex9.0 or higher. It requires access to the optimal tableau
      and optimal basis inverse and makes assumptions on the way slack
      variables are added by the solver. The Osi implementations for
      Clp and Cplex verify these assumptions.

      When calling the generator, the solver interface si must contain
      an optimized problem and information related to the optimal
      basis must be available through the OsiSolverInterface methods
      (si->optimalBasisIsAvailable() must return 'true'). It is also
      essential that the integrality of structural variable i can be
      obtained using si->isInteger(i).

      Reduce-and-Split cuts are a class of split cuts. We compute
      linear combinations of the rows of the simplex tableau, trying
      to reduce some of the coefficients on the nonbasic continuous
      columns.  We have a large number of heuristics to choose which
      coefficients should be reduced, and by using which rows. The
      paper explains everything in detail.

      Note that this generator can potentially generate a huge number
      of cuts, depending on how it is parametered. Default parameters
      should be good for most situations; if you want to go heavy on
      split cuts, use more row selection strategies or a different
      number of rows in the linear combinations. Again, look at the
      paper for details. If you want to generate a small number of
      cuts, default parameters are not the best choice.

      A combination of Reduce-and-Split with Lift & Project is
      described in the paper "Combining Lift-and-Project and
      Reduce-and-Split". The Reduce-and-Split code for the
      implementation used in that paper is included here.

      This generator does not generate the same cuts as CglRedSplit,
      therefore both generators can be used in conjunction.

  */

  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());

  /// Return true if needs optimal basis to do cuts (will return true)
  virtual bool needsOptimalBasis() const;

  // Generate the row multipliers computed by Reduce-and-Split from the
  // given OsiSolverInterface. The multipliers are written in lambda;
  // lambda should be of size nrow*maxNumMultipliers. We generate at most
  // maxNumMultipliers m-vectors of row multipliers, and return the number
  // of m-vectors that were generated.
  // If the caller wants to know which variables are basic in each row 
  // (same order as lambda), basicVariables should be non-NULL (size nrow).
  // This method can also generate the cuts corresponding to the multipliers
  // returned; it suffices to pass non-NULL OsiCuts.
  // This method is not needed by the typical user; however, it is useful
  // in the context of generating Lift & Project cuts.
  int generateMultipliers(const OsiSolverInterface& si, int* lambda,
			  int maxNumMultipliers, int* basicVariables = NULL,
			  OsiCuts* cs = NULL);

  // Try to improve a Lift & Project cut, by employing the
  // Reduce-and-Split procedure. We start from a row of a L&P tableau,
  // and generate a cut trying to reduce the coefficients on the
  // nonbasic variables.  Note that this L&P tableau will in general
  // have nonbasic variables which are nonzero in the point that we
  // want to cut off, so we should be careful.  Arguments:
  // OsiSolverInterface which contains the simplex tableau, initial
  // row from which the cut is derived, row rhs, row number of the
  // source row (if it is in the simplex tableau; otherwise, a
  // negative number; needed to avoid using duplicate rows), point
  // that we want to cut off (note: this is NOT a basic solution for
  // the OsiSolverInterace!), list of variables which are basic in
  // xbar but are nonbasic in the OsiSolverInterface. The computed cut
  // is written in OsiRowCut* cs. Finally, if a starting disjunction
  // is provided in the vector lambda (of size ncols, i.e. a
  // disjunction on the structural variables), the disjunction is
  // modified according to the cut which is produced.
  int tiltLandPcut(const OsiSolverInterface* si, double* row, 
		   double rowRhs, int rownumber, const double* xbar, 
		   const int* newnonbasics, OsiRowCut* cs, int* lambda = NULL);

  //@}
  
  
  /**@name Public Methods */
  //@{

  // Set the parameters to the values of the given CglRedSplit2Param object.
  void setParam(const CglRedSplit2Param &source); 
  // Return the CglRedSplit2Param object of the generator. 
  inline CglRedSplit2Param& getParam() {return param;}

  /// Print some of the data members; used for debugging
  void print() const;

  /// Print the current simplex tableau  
  void printOptTab(OsiSolverInterface *solver) const;
  
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglRedSplit2();

  /// Constructor with specified parameters 
  CglRedSplit2(const CglRedSplit2Param &RS_param);
 
  /// Copy constructor 
  CglRedSplit2(const CglRedSplit2 &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglRedSplit2 & operator=(const CglRedSplit2& rhs);
  
  /// Destructor 
  virtual ~CglRedSplit2 ();

  //@}

private:
  
  // Private member methods

/**@name Private member methods */

  //@{

  // Method generating the cuts after all CglRedSplit2 members are 
  // properly set. This does the actual work. Returns the number of
  // generated cuts (or multipliers).
  // Will generate cuts if cs != NULL, and will generate multipliers
  // if lambda != NULL. 
  int generateCuts(OsiCuts* cs, int maxNumCuts, int* lambda = NULL);

  /// Compute the fractional part of value, allowing for small error.
  inline double rs_above_integer(const double value) const; 

  /// Fill workNonBasicTab, depending on the column selection strategy.
  /// Accepts a list of variables indices that should be ignored; by
  /// default, this list is empty (it is only used by Lift & Project).
  /// The list ignore_list contains -1 as the last element.
  /// Note that the implementation of the ignore_list is not very efficient
  /// if the list is long, so it should be used only if its short.
  void fill_workNonBasicTab(CglRedSplit2Param::ColumnSelectionStrategy 
			    strategy, const int* ignore_list = NULL);

  /// Fill workNonBasicTab, alternate version for Lift & Project: also
  /// reduces columns which are now nonbasic but are basic in xbar.
  /// This function should be called only when CglRedSplit2 is used in
  /// conjunction with CglLandP to generate L&P+RS cuts.
  void fill_workNonBasicTab(const int* newnonbasics, const double* xbar,
			    CglRedSplit2Param::ColumnScalingStrategy scaling);

  /// Reduce rows of workNonBasicTab, i.e. compute integral linear
  /// combinations of the rows in order to reduce row coefficients on
  /// workNonBasicTab
  void reduce_workNonBasicTab(int numRows, 
			      CglRedSplit2Param::RowSelectionStrategy 
			      rowSelectionStrategy,
			      int maxIterations);

  /// Generate a linear combination of the rows of the current LP
  /// tableau, using the row multipliers stored in the matrix pi_mat
  /// on the row of index index_row
  void generate_row(int index_row, double *row);

  /// Generate a mixed integer Gomory cut, when all non basic 
  /// variables are non negative and at their lower bound.
  int generate_cgcut(double *row, double *rhs);

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
  void unflip(double *row, double *rowrhs);

  /// Returns 1 if the row has acceptable max/min coeff ratio.
  /// Compute max_coeff: maximum absolute value of the coefficients.
  /// Compute min_coeff: minimum absolute value of the coefficients
  /// larger than EPS_COEFF.
  /// Return 0 if max_coeff/min_coeff > MAXDYN.
  int check_dynamism(double *row);

  /// Generate the packed cut from the row representation.
  int generate_packed_row(const double *xlp, double *row,
			  int *rowind, double *rowelem, 
			  int *card_row, double & rhs);

  // Compute entries of is_integer.
  void compute_is_integer();

  // Check that two vectors are different.
  bool rs_are_different_vectors(const int *vect1, 
				const int *vect2,
				const int dim);

  // allocate matrix of integers
  void rs_allocmatINT(int ***v, int m, int n);
  // deallocate matrix of integers
  void rs_deallocmatINT(int ***v, int m);
  // allocate matrix of doubles
  void rs_allocmatDBL(double ***v, int m, int n);
  // deallocate matrix of doubles
  void rs_deallocmatDBL(double ***v, int m);
  // print a vector of integers
  void rs_printvecINT(const char *vecstr, const int *x, int n) const;
  // print a vector of doubles
  void rs_printvecDBL(const char *vecstr, const double *x, int n) const;
  // print a matrix of integers
  void rs_printmatINT(const char *vecstr, const int * const *x, int m, int n) const;
  // print a matrix of doubles
  void rs_printmatDBL(const char *vecstr, const double * const *x, int m, int n) const;
  // dot product
  double rs_dotProd(const double *u, const double *v, int dim) const;
  double rs_dotProd(const int *u, const double *v, int dim) const;
  // From Numerical Recipes in C: LU decomposition
  int ludcmp(double **a, int n, int *indx, double *d, double* vv) const;
  // from Numerical Recipes in C: backward substitution
  void lubksb(double **a, int n, int *indx, double *b) const;

  // Check if the linear combination given by listOfRows with given multipliers
  // improves the norm of row #rowindex; note: multipliers are rounded!
  // Returns the difference with respect to the old norm (if negative there is
  // an improvement, if positive norm increases)
  double compute_norm_change(double oldnorm, const int* listOfRows,
			     int numElemList, const double* multipliers) const;

  // Compute the list of rows that should be used to reduce row #rowIndex
  int get_list_rows_reduction(int rowIndex, int numRowsReduction, 
			      int* list, const double* norm, 
			      CglRedSplit2Param::RowSelectionStrategy 
			      rowSelectionStrategy) const;

  // Sorts the rows by increasing number of nonzeroes with respect to a given
  // row (rowIndex), on the nonbasic variables (whichTab == 0 means only
  // integer, whichTab == 1 means only workTab, whichTab == 2 means both).
  // The array for sorting must be allocated (and deleted) by caller.
  // Corresponds to BRS1 in the paper.
  int sort_rows_by_nonzeroes(struct sortElement* array, int rowIndex, 
			     int maxRows, int whichTab) const;

  // Greedy variant of the previous function; slower but typically
  // more effective. Corresponds to BRS2 in the paper.
  int sort_rows_by_nonzeroes_greedy(struct sortElement* array, int rowIndex, 
				    int maxRows, int whichTab) const;

  // Sorts the rows by decreasing absolute value of the cosine of the
  // angle with respect to a given row (rowIndex), on the nonbasic
  // variables (whichTab == 0 means only integer, whichTab == 1 means
  // only workTab, whichTab == 2 means both).  The array for sorting
  // must be allocated (and deleted) by caller. Very effective
  // strategy in practice. Corresponds to BRS3 in the paper.
  int sort_rows_by_cosine(struct sortElement* array, int rowIndex, 
			  int maxRows, int whichTab) const;
#ifdef CHECK_SPLIT_TIME // incredibly expensive!
  // Did we hit the time limit?
  inline bool checkTime() const{
    if ((CoinCpuTime() - startTime) < param.getTimeLimit()){
      return true;
    }
    return false;
  }
#else
  inline bool checkTime() const{
    return true;
  }
#endif

  //@}

  
  // Private member data

  /**@name Private member data */
  
  //@{

  /// Object with CglRedSplit2Param members. 
  CglRedSplit2Param param;

  /// Number of rows ( = number of slack variables) in the current LP.
  int nrow; 

  /// Number of structural variables in the current LP.
  int ncol;

  /// Number of rows which have been reduced
  int numRedRows;

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

  /// Reduced costs for columns
  const double *reducedCost;

  /// Row price
  const double *rowPrice;

  /// Objective coefficients
  const double* objective;

  /// Number of integer basic structural variables 
  int card_intBasicVar;

  /// Number of integer basic structural variables that are fractional in the
  /// current lp solution (at least param.away_ from being integer).  
  int card_intBasicVar_frac;

  /// Number of integer non basic structural variables in the
  /// current lp solution.  
  int card_intNonBasicVar; 

  /// Number of continuous non basic variables (structural or slack) in the
  /// current lp solution.  
  int card_contNonBasicVar;

  /// Number of continuous non basic variables (structural or slack) in the
  /// current working set for coefficient reduction
  int card_workNonBasicVar;

  /// Number of non basic variables (structural or slack) at their
  /// upper bound in the current lp solution.
  int card_nonBasicAtUpper; 

  /// Number of non basic variables (structural or slack) at their
  /// lower bound in the current lp solution.
  int card_nonBasicAtLower;

  /// Characteristic vector for integer basic structural variables
  int *cv_intBasicVar;  

  /// Characteristic vector for integer basic structural variables
  /// with non integer value in the current lp solution.
  int *cv_intBasicVar_frac;  

  /// Characteristic vector for rows of the tableau selected for reduction
  /// with non integer value in the current lp solution
  int *cv_fracRowsTab;  

  /// List of integer structural basic variables 
  /// (in order of pivot in selected rows for cut generation).
  int *intBasicVar;

  /// List of integer structural basic variables with fractional value
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

  /// Number of rows in the reduced tableau (= card_intBasicVar).
  int mTab;

  /// Number of columns in the reduced tableau (= card_contNonBasicVar)
  int nTab;

  /// Tableau of multipliers used to alter the rows used in generation.
  /// Dimensions: mTab by mTab. Initially, pi_mat is the identity matrix.
  int **pi_mat;

  /// Simplex tableau for continuous non basic variables (structural or slack).
  /// Only rows used for generation.
  /// Dimensions: mTab by card_contNonBasicVar.
  double **contNonBasicTab;

  /// Current tableau for continuous non basic variables (structural or slack).
  /// Only columns used for coefficient reduction.
  /// Dimensions: mTab by card_workNonBasicVar.
  double **workNonBasicTab;

  /// Simplex tableau for integer non basic structural variables.
  /// Only rows used for generation.
  // Dimensions: mTab by card_intNonBasicVar.
  double **intNonBasicTab;

  /// Right hand side of the tableau.
  /// Only rows used for generation.
  double *rhsTab;

  /// Norm of rows in workNonBasicTab; needed for faster computations
  double *norm;

  /// Characteristic vectors of structural integer variables or continuous
  /// variables currently fixed to integer values. 
  int *is_integer;

  /// Pointer on solver. Reset by each call to generateCuts().
  OsiSolverInterface *solver;

  /// Pointer on point to separate. Reset by each call to generateCuts().
  const double *xlp;

  /// Pointer on row activity. Reset by each call to generateCuts().
  const double *rowActivity;

  /// Pointer on matrix of coefficient ordered by rows. 
  /// Reset by each call to generateCuts().
  const CoinPackedMatrix *byRow;

  /// Time at which cut computations began.
  /// Reset by each call to generateCuts().
  double startTime;

  //@}
};

//#############################################################################
/** A function that tests some of the methods in the CglRedSplit2
    class. The only reason for it not to be a member method is that
    this way it doesn't have to be compiled into the library. And
    that's a gain, because the library should be compiled with
    optimization on, but this method should be compiled with
    debugging. */
CGLLIB_EXPORT
void CglRedSplit2UnitTest(const OsiSolverInterface * siP,
			 const std::string mpdDir );


#endif
