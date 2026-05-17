// Last edit: 02/05/2013
//
// Name:     CglGMI.hpp
// Author:   Giacomo Nannicini
//           Singapore University of Technology and Design, Singapore
//           email: nannicini@sutd.edu.sg
// Date:     11/17/09
//-----------------------------------------------------------------------------
// Copyright (C) 2009, Giacomo Nannicini.  All Rights Reserved.

#ifndef CglGMI_H
#define CglGMI_H

#include "CglCutGenerator.hpp"
#include "CglGMIParam.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinFactorization.hpp"

/* Enable tracking of rejection of cutting planes. If this is disabled,
   the cut generator is slightly faster. If defined, it enables proper use
   of setTrackRejection and related functions. */
//#define TRACK_REJECT

/* Debug output */
//#define GMI_TRACE

/* Debug output: print optimal tableau */
//#define GMI_TRACETAB

/* Print reason for cut rejection, whenever a cut is discarded */
//#define GMI_TRACE_CLEAN

/** Gomory cut generator with several cleaning procedures, used to test
 *  the numerical safety of the resulting cuts 
 */

class CGLLIB_EXPORT CglGMI : public CglCutGenerator {

  friend CGLLIB_EXPORT void CglGMIUnitTest(const OsiSolverInterface * siP,
			     const std::string mpdDir);
public:

  /** Public enum: all possible reasons for cut rejection */
  enum RejectionType{
    failureFractionality,
    failureDynamism,
    failureViolation,
    failureSupport,
    failureScale
  };

  /**@name generateCuts */
  //@{
  /** Generate Gomory Mixed-Integer cuts for the model of the solver
      interface si.

      Insert the generated cuts into OsiCuts cs.

      Warning: This generator currently works only with the Lp solvers Clp or 
      Cplex9.0 or higher. It requires access to the optimal tableau and 
      optimal basis inverse and makes assumptions on the way slack variables 
      are added by the solver. The Osi implementations for Clp and Cplex 
      verify these assumptions.

      When calling the generator, the solver interface si must contain
      an optimized problem and information related to the optimal
      basis must be available through the OsiSolverInterface methods
      (si->optimalBasisIsAvailable() must return 'true'). It is also
      essential that the integrality of structural variable i can be
      obtained using si->isInteger(i).

  */
  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());

  /// Return true if needs optimal basis to do cuts (will return true)
  virtual bool needsOptimalBasis() const { return true; }
  //@}

  /**@name Common Methods */
  //@{
  // Function for checking equality with user tolerance
  inline bool areEqual(double x, double y, 
		       double epsAbs = 1e-12, 
		       double epsRel = 1e-12) {
    return (fabs((x) - (y)) <= 
	    std::max(epsAbs, epsRel * std::max(fabs(x), fabs(y))));
  }

  // Function for checking is a number is zero
  inline bool isZero(double x, double epsZero = 1e-20) {
    return (fabs(x) <= epsZero);
  }


  // Function for checking if a number is integer
  inline bool isIntegerValue(double x, 
			     double intEpsAbs = 1e-9,
			     double intEpsRel = 1e-15) {
    return (fabs((x) - floor((x)+0.5)) <= 
	    std::max(intEpsAbs, intEpsRel * fabs(x)));
  }

  
  //@}
  
  
  /**@name Public Methods */
  //@{

  // Set the parameters to the values of the given CglGMIParam object.
  void setParam(const CglGMIParam &source); 
  // Return the CglGMIParam object of the generator. 
  inline CglGMIParam getParam() const {return param;}
  inline CglGMIParam & getParam() {return param;}

  // Compute entries of is_integer.
  void computeIsInteger();

  /// Print the current simplex tableau  
  void printOptTab(OsiSolverInterface *solver) const;

  /// Set/get tracking of the rejection of cutting planes.
  /// Note that all rejection related functions will not do anything 
  /// unless the generator is compiled with the define GMI_TRACK_REJECTION
  void setTrackRejection(bool value);
  bool getTrackRejection();

  /// Get number of cuts rejected for given reason; see above
  int getNumberRejectedCuts(RejectionType reason);

  /// Reset counters for cut rejection tracking; see above
  void resetRejectionCounters();

  /// Get total number of generated cuts since last resetRejectionCounters()
  int getNumberGeneratedCuts();
  
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglGMI();

  /// Constructor with specified parameters 
  CglGMI(const CglGMIParam &param);
 
  /// Copy constructor 
  CglGMI(const CglGMI &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglGMI & operator=(const CglGMI& rhs);
  
  /// Destructor 
  virtual ~CglGMI();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);

  //@}
    
private:
  
  // Private member methods

/**@name Private member methods */

  //@{

  // Method generating the cuts after all CglGMI members are properly set.
  void generateCuts(OsiCuts & cs);

  /// Compute the fractional part of value, allowing for small error.
  inline double aboveInteger(double value) const; 

  /// Compute the fractionalities involved in the cut, and the cut rhs.
  /// Returns true if cut is accepted, false if discarded
  inline bool computeCutFractionality(double varRhs, double& cutRhs);

  /// Compute the cut coefficient on a given variable
  inline double computeCutCoefficient(double rowElem, int index);

  /// Use multiples of the initial inequalities to cancel out the coefficient
  /// on a slack variables. 
  inline void eliminateSlack(double cutElem, int cutIndex, double* cut,
			      double& cutRhs, const double *elements, 
			      const CoinBigIndex *rowStart, const int *indices, 
			      const int *rowLength, const double *rhs);

  /// Change the sign of the coefficients of the non basic
  /// variables at their upper bound.
  inline void flip(double& rowElem, int rowIndex);

  /// Change the sign of the coefficients of the non basic
  /// variables at their upper bound and do the translations restoring
  /// the original bounds. Modify the right hand side
  /// accordingly. Two functions: one for original variables, one for slacks.
  inline void unflipOrig(double& rowElem, int rowIndex, double& rowRhs);
  inline void unflipSlack(double& rowElem, int rowIndex, double& rowRhs,
			   const double* slack_val);

  /// Pack a row of ncol elements
  inline void packRow(double* row, double* rowElem, int* rowIndex,
		       int& rowNz);

  /// Clean the cutting plane; the cleaning procedure does several things
  /// like removing small coefficients, scaling, and checks several
  /// acceptance criteria. If this returns false, the cut should be discarded.
  /// There are several cleaning procedures available, that can be selected
  /// via the parameter param.setCLEANING_PROCEDURE(int value)
  bool cleanCut(double* cutElem, int* cutIndex, int& cutNz,
		 double& cutRhs, const double* xbar);

  /// Cut cleaning procedures: return true if successfull, false if
  /// cut should be discarded by the caller of if problems encountered

  /// Check the violation
  bool checkViolation(const double* cutElem, const int* cutIndex,
		       int cutNz, double cutrhs, const double* xbar);

  /// Check the dynamism
  bool checkDynamism(const double* cutElem, const int* cutIndex,
		      int cutNz);

  /// Check the support
  bool checkSupport(int cutNz);

  /// Remove small coefficients and adjust the rhs accordingly
  bool removeSmallCoefficients(double* cutElem, int* cutIndex, 
				 int& cutNz, double& cutRhs);

  /// Adjust the rhs by relaxing by a small amount (relative or absolute)
  void relaxRhs(double& rhs);

  /// Scale the cutting plane in different ways;
  /// scaling_type possible values:
  /// 0 : scale to obtain integral cut
  /// 1 : scale based on norm, to obtain cut norm equal to ncol
  /// 2 : scale to obtain largest coefficient equal to 1
  bool scaleCut(double* cutElem, int* cutIndex, int cutNz,
		 double& cutRhs, int scalingType);

  /// Scale the cutting plane in order to generate integral coefficients
  bool scaleCutIntegral(double* cutElem, int* cutIndex, int cutNz,
			  double& cutRhs);

  /// Compute the nearest rational number; used by scale_row_integral
  bool nearestRational(double val, double maxdelta, int64_t maxdnom,
			int64_t& numerator, int64_t& denominator);

  /// Compute the greatest common divisor
  int64_t computeGcd(int64_t a, int64_t b);

  /// print a vector of integers
  void printvecINT(const char *vecstr, const int *x, int n) const;
  /// print a vector of doubles: dense form
  void printvecDBL(const char *vecstr, const double *x, int n) const;
  /// print a vector of doubles: sparse form
  void printvecDBL(const char *vecstr, const double *elem, const int * index, 
		   int nz) const;

  /// Recompute the simplex tableau for want of a better accuracy.
  /// Requires an empty CoinFactorization object to do the computations,
  /// and two empty (already allocated) arrays which will contain 
  /// the basis indices on exit. Returns 0 if successfull.
  int factorize(CoinFactorization & factorization,
		int* colBasisIndex, int* rowBasisIndex);


  //@}

  
  // Private member data

/**@name Private member data */

  //@{

  /// Object with CglGMIParam members. 
  CglGMIParam param;

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

  /// Characteristic vectors of structural integer variables or continuous
  /// variables currently fixed to integer values. 
  bool *isInteger;

  /// Current basis status: columns
  int *cstat;

  /// Current basis status: rows
  int *rstat;

  /// Pointer on solver. Reset by each call to generateCuts().
  OsiSolverInterface *solver;

  /// Pointer on point to separate. Reset by each call to generateCuts().
  const double *xlp;

  /// Pointer on row activity. Reset by each call to generateCuts().
  const double *rowActivity;

  /// Pointer on matrix of coefficient ordered by rows. 
  /// Reset by each call to generateCuts().
  const CoinPackedMatrix *byRow;

  /// Pointer on matrix of coefficient ordered by columns. 
  /// Reset by each call to generateCuts().
  const CoinPackedMatrix *byCol;

  /// Fractionality of the cut and related quantities.
  double f0;
  double f0compl;
  double ratiof0compl;

#if defined(TRACK_REJECT) || defined (TRACK_REJECT_SIMPLE)
  /// Should we track the reason of each cut rejection?
  bool trackRejection;
  /// Number of failures by type
  int fracFail;
  int dynFail;
  int violFail;
  int suppFail;
  int smallCoeffFail;
  int scaleFail;
  /// Total number of generated cuts
  int numGeneratedCuts;
#endif

  //@}
};

//#############################################################################
/** A function that tests the methods in the CglGMI class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
CGLLIB_EXPORT
void CglGMIUnitTest(const OsiSolverInterface * siP,
			 const std::string mpdDir );


#endif
