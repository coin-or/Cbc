// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglTwomir_H
#define CglTwomir_H
#include <string>

#include "CglCutGenerator.hpp"
#include "CoinFactorization.hpp"

typedef struct
{

  int nz;             /* current length of arrays index[] and coeff[] */
  int max_nz;         /* max length of arrays index[] and coeff[] */
  double *coeff;      /* coefficient of each variable in the constraint */
  int *index;         /* index of the variable (value in 0 ... nrow+ncol) */
  double rhs;         /* rhs of the constraint */
  char sense;         /* ?? is it necessary */

} DGG_constraint_t;

typedef struct{
  int n;
  DGG_constraint_t **c;
  int *ctype;
  double *alpha;
} DGG_list_t;

/******************** BASIS INFORMATION ADTs **********************************/
typedef struct{
  int q_min;
  int q_max;
  int t_min;
  int t_max;
  int a_max;
  int max_elements;
} cutParams;

#define TWOMIR_LESS_MALLOC
typedef struct
{
  double gomory_threshold; /* factional variable must be this away from int */
  int ncol,        /* number of columns in LP */
    nrow,        /* number of constaints in LP */
    ninteger;    /* number of integer variables in LP */

  int nbasic_col,  /* number of basic columns in the LP */
    nbasic_row;  /* number of basic rows in the LP */

  /* the following arrays are all of size (ncol+nrow) */
  int *info;       /* description of each variable (see below) */
  double *lb;      /* specifies the lower bound (if any) of each variable */
  double *ub;      /* specifies the upper bound (if any) of each variable */
  double *x;       /* current solution */
  double *rc;      /* current reduced cost */
  double *opt_x;
#ifdef TWOMIR_LESS_MALLOC
  /* the following arrays are to avoid many mallocs
   - they are used in DGG_transformConstraint and similar. */
  double *xtemp_;     /* values for cuts */
  double *rctemp_;
  double *tabrowcoeff_;
  int *skalatemp_;
  int *tabrowindex_;
  DGG_constraint_t *constrainttemp0_;
  DGG_constraint_t *constrainttemp1_;
  char *pitemp_;
  int *rowIsBasictemp_;
  int *colIsBasictemp_;
  CoinIndexedVector * vector0_;
  CoinIndexedVector * vector1_;
  double * spareArray0_;
  void * spare1;
  void * spare2;
#endif

  cutParams cparams;
} DGG_data_t;

/* the following macros allow us to decode the info of the DGG_data
   type. The encoding is as follows,
   bit 1 : if the variable is basic or not (non-basic).
   bit 2 : if the variable is integer or or not (rational).
   bit 3 : if the variable is structural or not (artifical). 
   bit 4 : if the variable is non-basic and at its upper bound 
   (else if non-basic at lower bound). */

#define DGG_isBasic(data,idx) ((data->info[idx])&1)
#define DGG_isInteger(data,idx) ((data->info[idx] >> 1)&1)
#define DGG_isStructural(data,idx) ((data->info[idx] >> 2)&1)
#define DGG_isEqualityConstraint(data,idx) ((data->info[idx] >> 3)&1)
#define DGG_isNonBasicAtUB(data,idx) ((data->info[idx] >> 4)&1)
#define DGG_isNonBasicAtLB(data,idx) ((data->info[idx] >> 5)&1)
#define DGG_isConstraintBoundedAbove(data,idx) ((data->info[idx] >> 6)&1)
#define DGG_isConstraintBoundedBelow(data,idx) ((data->info[idx] >> 7)&1)

#define DGG_setIsBasic(data,idx) ((data->info[idx]) |= 1)
#define DGG_setIsInteger(data,idx) ((data->info[idx]) |= (1<<1))
#define DGG_setIsStructural(data,idx) ((data->info[idx]) |= (1<<2))
#define DGG_setEqualityConstraint(data,idx) ((data->info[idx]) |= (1<<3))
#define DGG_setIsNonBasicAtUB(data,idx) ((data->info[idx]) |= (1<<4))
#define DGG_setIsNonBasicAtLB(data,idx) ((data->info[idx]) |= (1<<5))
#define DGG_setIsConstraintBoundedAbove(data,idx) ((data->info[idx]) |= (1<<6))
#define DGG_setIsConstraintBoundedBelow(data,idx) ((data->info[idx]) |= (1<<7))

class CoinWarmStartBasis;
/** Twostep MIR Cut Generator Class */
class CGLLIB_EXPORT CglTwomir : public CglCutGenerator {

  friend CGLLIB_EXPORT void CglTwomirUnitTest(const OsiSolverInterface * siP,
					  const std::string mpdDir );


public:

  /// Problem name
  std::string probname_;
    
  /**@name Generate Cuts */
  //@{
  /** Generate Two step MIR cuts either from the tableau rows or from the
      formulation rows
  */
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs, 
			     const CglTreeInfo info = CglTreeInfo());
  /// Return true if needs optimal basis to do cuts (will return true)
  virtual bool needsOptimalBasis() const;

  /**@name Change criterion on which scalings to use (default = 1,1,1,1) */
  //@{
  /// Set
  void setMirScale (int tmin, int tmax) {t_min_ = tmin; t_max_ = tmax;}
  void setTwomirScale (int qmin, int qmax) {q_min_ = qmin; q_max_ = qmax;}
  void setAMax (int a) {a_max_ = a;}
  void setMaxElements (int n) {max_elements_ = n;}
  void setMaxElementsRoot (int n) {max_elements_root_ = n;}
  void setCutTypes (bool mir, bool twomir, bool tab, bool form)
  { do_mir_ = mir; do_2mir_ = twomir; do_tab_ = tab; do_form_ = form;}
  void setFormulationRows (int n) {form_nrows_ = n;}

  /// Get
  int getTmin() const {return t_min_;}
  int getTmax() const {return t_max_;}
  int getQmin() const {return q_min_;}
  int getQmax() const {return q_max_;}
  int getAmax() const {return a_max_;}
  int getMaxElements() const {return max_elements_;}
  int getMaxElementsRoot() const {return max_elements_root_;}
  int getIfMir() const { return do_mir_;}
  int getIfTwomir() const { return do_2mir_;}
  int getIfTableau() const { return do_tab_;}
  int getIfFormulation() const { return do_form_;}
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
  /// Return maximum length of cut in tree
  virtual int maximumLengthOfCutInTree() const
  { return max_elements_;}
  //@}

  /**@name Change way TwoMir works */
  //@{
  /// Set type - 0 normal, 1 add original matrix one, 2 replace
  inline void setTwomirType(int type)
  { twomirType_=type;}
  /// Return type
  inline int twomirType() const
  { return twomirType_;}
  //@}
  /// Pass in a copy of original solver (clone it)
  void passInOriginalSolver(OsiSolverInterface * solver);

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglTwomir ();

  /// Copy constructor 
  CglTwomir (const CglTwomir &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglTwomir & operator=(const CglTwomir& rhs);
  
  /// Destructor 
  virtual  ~CglTwomir ();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  /// This can be used to refresh any inforamtion
  virtual void refreshSolver(OsiSolverInterface * solver);
  //@}
      
private:
  // Private member data
  /**@name Private member data */
  //@{
  /// Threadsafe random number generator
  CoinThreadRandom randomNumberGenerator_;
  /// Only investigate if more than this away from integrality
  double away_;
  /// Only investigate if more than this away from integrality (at root)
  double awayAtRoot_;
  /// Type - 0 normal, 1 add original matrix one, 2 replace
  int twomirType_;
  bool do_mir_;
  bool do_2mir_;
  bool do_tab_;
  bool do_form_;

  int t_min_;  /// t_min - first value of t to use for tMIR inequalities
  int t_max_;  /// t_max - last value of t to use for tMIR inequalities
  int q_min_;  /// q_min - first value of t to use for 2-Step tMIR inequalities
  int q_max_;  /// q_max - last value of t to use for 2-Step tMIR inequalities
  int a_max_;  /// a_max - maximum value of bhat/alpha
  int max_elements_; /// Maximum number of elements in cut
  int max_elements_root_; /// Maximum number of elements in cut at root
  int form_nrows_; //number of rows on which formulation cuts will be generated
  //@}
};

//#############################################################################

/*
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include <cassert>
#include <iostream.h>
*/

/******************** DEBUG DEFINITIONS ***************************************/

#define DGG_DEBUG_DGG 1
#define DGG_TRACE_ERRORS 0
#define DGG_DISPLAY   0
#define DGG_AUTO_CHECK_CUT_OFF_OPTIMAL 1

/******************** CONFIGURATION DEFAULTS **********************************/

#define DGG_DEFAULT_METHOD 2
#define DGG_DEFAULT_TMIN 1
#define DGG_DEFAULT_TMAX 1
#define DGG_DEFAULT_TAUMIN 2
#define DGG_DEFAULT_TAUMAX 6
#define DGG_DEFAULT_MAX_CUTS 500
#define DGG_DEFAULT_IMPROVEMENT_THRESH 0.001
#define DGG_DEFAULT_NBELOW_THRESH INT_MAX 
#define DGG_DEFAULT_NROOT_ROUNDS 2
#define DGG_DEFAULT_NEGATIVE_SCALED_TWOSTEPS 0
#define DGG_DEFAULT_ALPHA_RULE 0
#define DGG_DEFAULT_CUT_INC 250
#define DGG_DEFAULT_CUT_FORM 0
#define DGG_DEFAULT_NICEFY 0
#define DGG_DEFAULT_ONLY_DELAYED 0
#define DGG_DEFAULT_DELAYED_FREQ 9999999 
#define DGG_DEFAULT_LPROWS_FREQ 9999999
#define DGG_DEFAULT_WHICH_FORMULATION_CUTS 2

/******************** SOLVER CONFIGURATION DEFINITIONS ************************/

#define DGG_OSI 0
#define DGG_CPX 1
#define DGG_QSO 2

/* determines the solver to be used */
#define DGG_SOLVER DGG_OSI

/* adds checking routines to make sure solver works as expected */
#define DGG_DEBUG_SOLVER 0

/* turn off screen output from solver */
#define DGG_SOLVER_SCREEN_FLAG 0

/******************** CUT DEFINITIONS *****************************************/

/* internal names for cut types */
#define DGG_TMIR_CUT 1
#define DGG_2STEP_CUT 2

/* internal names for alpha-selection rules */
#define DGG_ALPHA_MIN_SUM 0
#define DGG_ALPHA_RANDOM_01 1
#define DGG_ALPHA_RANDOM_COEFF 2
#define DGG_ALPHA_ALL 3
#define DGG_ALPHA_MAX_STEEP 5

/******************** PRECISION & NUMERICAL ISSUES DEFINITIONS ****************/

/* how steep a cut must be before adding it to the lp */
#define DGG_MIN_STEEPNESS 1.0e-4
#define DGG_MAX_L2NORM 1.0e7

/* 0 = min steepness, 1 = max norm */
#define DGG_NORM_CRITERIA 1

/* used to define how fractional a basic-integer variable must be
   before choosing to use it to generate a TMIR cut on.
   OSI's default is 1.0e-7 */
#define DGG_GOMORY_THRESH 0.005

#define DGG_RHS_THRESH 0.005

/* used for comparing variables to their upper bounds.
   OSI's default is 1.0e-7.
   We set it to 1.0e6 because e-7 seems too sensitive. 
   In fact, with e-7 the problem dsbmip.mps complains. */
#define DGG_BOUND_THRESH 1.0e-6

/* used for comparing the lhs (activity) value of a tableau row
   with the rhs. This is only used for debugging purposes. */
#define DGG_EQUALITY_THRESH 1.0e-5

/* used for comparing a variable's lower bound to 0.0
   and determining if we need to shift the variable    */
#define DGG_SHIFT_THRESH 1.0e-6

/* used for determing how far from an integer is still an integer.
   This value is used for comparing coefficients to integers.
   OSI's default is 1.0e-10.                               */
#define DGG_INTEGRALITY_THRESH 1.0e-10

/* the min value that a coeff can have in the tableau row
   before being set to zero. */
#define CBC_CHECK_CUT
#ifndef CBC_CHECK_CUT
#define DGG_MIN_TABLEAU_COEFFICIENT 1.0e-8
#else
#define DGG_MIN_TABLEAU_COEFFICIENT 1.0e-12
#endif

/* smallest value rho is allowed to have for a simple 2-step MIR
   (ie: not an extended two-step MIR) */
#define DGG_MIN_RHO 1.0e-7
#define DGG_MIN_ALPHA 1.0e-7

/* when a slack is null: used to check if a cut is satisfied or not. */
#define DGG_NULL_SLACK 1.0e-5

/* nicefy constants */
#define DGG_NICEFY_MIN_ABSVALUE 1.0e-13
#define DGG_NICEFY_MIN_FIX 1.0e-7
#define DGG_NICEFY_MAX_PADDING 1.0e-6
#define DGG_NICEFY_MAX_RATIO 1.0e9


/******************** ERROR-CATCHING MACROS ***********************************/
#if DGG_TRACE_ERRORS > 0

#define __DGG_PRINT_LOC__(F) fprintf(((F==0)?stdout:F), " in %s (%s:%d)\n", __func__, __FILE__, __LINE__)

#define DGG_THROW(A,REST...) {\
 fprintf(stdout, ##REST); \
 __DGG_PRINT_LOC__(stdout); \
 return (A);}

#define DGG_IF_EXIT(A,B,REST...) {\
 if(A) {\
 fprintf(stdout, ##REST); \
 __DGG_PRINT_LOC__(stdout); \
 exit(B);}}

#define DGG_CHECKRVAL(A,B) {\
 if(A) {\
   __DGG_PRINT_LOC__(stdout); \
   return B; } }

#define DGG_CHECKRVAL1(A,B) {\
 if(A) {\
   __DGG_PRINT_LOC__(stdout); \
   rval = B; goto CLEANUP; } }

#define DGG_WARNING(A, REST...) {\
  if(A) {\
	  fprintf(stdout, ##REST); \
		__DGG_PRINT_LOC__(stdout); \
		}}

#define DGG_TEST(A,B,REST...) {\
 if(A) DGG_THROW(B,##REST) }

#define DGG_TEST2(A,B,C,REST)   {DGG_TEST(A,B,C,REST) }
#define DGG_TEST3(A,B,C,D,REST) {DGG_TEST(A,B,C,D,REST) }

#else

#define DGG_IF_EXIT(A,B,REST) {if(A) {fprintf(stdout, REST);exit(B);}}

#define DGG_THROW(A,B) return(A)

#define DGG_CHECKRVAL(A,B) {  if(A) return(B); }
#define DGG_CHECKRVAL1(A,B){ if(A) { rval = B; goto CLEANUP; } }

#define DGG_TEST(A,B,REST) { if(A) return(B);}
#define DGG_TEST2(A,B,REST,C) { DGG_TEST(A,B,REST) }
#define DGG_TEST3(A,B,REST,C,D) { DGG_TEST(A,B,REST) }

#endif

/******************** SIMPLE MACROS AND FUNCTIONS *****************************/

#define DGG_MIN(a,b) ( (a<b)?a:b )
#define DGG_MAX(a,b) ( (a>b)?a:b )
#define KREM(vht,alpha,tau)  (DGG_MIN( ceil(vht / alpha), tau ) - 1)
#define LMIN(vht, d, bht) (DGG_MIN( floor(d*bht/bht), d))
#define ABOV(v) (v - floor(v))
#define QINT(vht,bht,tau) ( (int)floor( (vht*(tau-1))/bht ) )
#define V2I(bht,tau,i) ( ((i+1)*bht / tau) )

int DGG_is_even(double vht, double bht, int tau, int q);
inline double frac_part(double value)
{
  return value-floor(value);
}
int DGG_is_a_multiple_of_b(double a, double b);

/* free function for DGG_data_t. Frees internal arrays and data structure */
int DGG_freeData( DGG_data_t *data );

/******************** CONSTRAINT ADTs *****************************************/
DGG_constraint_t* DGG_newConstraint(int max_arrays);
void DGG_freeConstraint(DGG_constraint_t *c);
DGG_constraint_t *DGG_copyConstraint(DGG_constraint_t *c);
void DGG_scaleConstraint(DGG_constraint_t *c, int t);
#ifdef TWOMIR_LESS_MALLOC
DGG_constraint_t* DGG_newConstraint(DGG_data_t *data, int which); 
DGG_constraint_t *DGG_copyConstraint(DGG_constraint_t* c,
				     DGG_data_t *data, int which);
#endif
/******************** CONFIGURATION *******************************************/
void DGG_list_init (DGG_list_t *l);
int DGG_list_addcut (DGG_list_t *l, DGG_constraint_t *cut, int ctype, double alpha);
void DGG_list_delcut (DGG_list_t *l, int i);
void DGG_list_free(DGG_list_t *l);

/******************* SOLVER SPECIFIC METHODS **********************************/
DGG_data_t *DGG_getData(const void *solver_ptr);

/******************* CONSTRAINT MANIPULATION **********************************/

/* DGG_transformConstraint: manipulates a constraint in the following way: 

packs everything in output

1 - variables at their upper bounds are substituted for their 
complements. This is done by adjusting the coefficients and 
the right hand side (simple substitution). 

2 - variables with non-zero lower bounds are shifted.            */

int DGG_transformConstraint( DGG_data_t *data,
                             double **x_out, 
			     double **rc_out,
                             char **isint_out,
                             DGG_constraint_t *constraint );

/* DGG_unTransformConstraint : 

1 - Undoes step (1) of DGG_transformConstraint 
2 - Undoes step (2) of DGG_transformConstraint                  */
 
int DGG_unTransformConstraint( DGG_data_t *data, 
                               DGG_constraint_t *constraint );

/* substitutes each slack variable by the structural variables which 
   define it. This function, hence, changes the constraint 'cut'.    */

int DGG_substituteSlacks( const void *solver_ptr, 
                          DGG_data_t *data, 
                          DGG_constraint_t *cut );

int DGG_nicefyConstraint( const void *solver_ptr, 
                          DGG_data_t *data,
			  DGG_constraint_t *cut);

/******************* CUT GENERATION *******************************************/
int DGG_getFormulaConstraint( int row_idx,  
                              const void *solver_ptr,   
			      DGG_data_t *data, 
                              DGG_constraint_t* row );

int DGG_getTableauConstraint( int index, 
                              const void *solver_ptr, 
                              DGG_data_t *data, 
                              DGG_constraint_t* tabrow,
                              const int * colIsBasic,
                              const int * rowIsBasic,
                              CoinFactorization & factorization,
                              int mode );

DGG_constraint_t* DGG_getSlackExpression(const void *solver_ptr, DGG_data_t* data, int row_index);

  int DGG_generateTabRowCuts( DGG_list_t *list,
			      DGG_data_t *data,
			      const void *solver_ptr );

  int DGG_generateFormulationCuts( DGG_list_t *list,
				   DGG_data_t *data,
				   const void *solver_ptr,
				   int nrows,
				   CoinThreadRandom & generator);


  int DGG_generateFormulationCutsFromBase( DGG_constraint_t *base,
					   double slack,
					   DGG_list_t *list,
					   DGG_data_t *data,
					   const void *solver_ptr,
					   CoinThreadRandom & generator);

  int DGG_generateCutsFromBase( DGG_constraint_t *base,
				DGG_list_t *list,
				DGG_data_t *data,
				const void *solver_ptr );

int DGG_buildMir( char *isint,
                  DGG_constraint_t *base,
                  DGG_constraint_t **cut_out );

int DGG_build2step( double alpha,
                    char *isint,
                    DGG_constraint_t *base,
                    DGG_constraint_t **cut_out );

  int DGG_addMirToList   ( DGG_constraint_t *base,
			   char *isint,
			   double *x,
			   DGG_list_t *list,
			   DGG_data_t *data,
			   DGG_constraint_t *orig_base );

  int DGG_add2stepToList ( DGG_constraint_t *base,
			   char *isint,
			   double *x,
			   double *rc,
			   DGG_list_t *list,
			   DGG_data_t *data,
			   DGG_constraint_t *orig_base );

/******************* CUT INFORMATION ******************************************/

double DGG_cutLHS(DGG_constraint_t *c, double *x);
int DGG_isCutDesirable(DGG_constraint_t *c, DGG_data_t *d);

/******************* TEST / DEBUGGING ROUTINES ********************************/

int DGG_isConstraintViolated(DGG_data_t *d, DGG_constraint_t *c);

int DGG_isBaseTrivial(DGG_data_t *d, DGG_constraint_t* c);
int DGG_is2stepValid(double alpha, double bht);

int DGG_cutsOffPoint(double *x, DGG_constraint_t *cut);

//#############################################################################
/** A function that tests the methods in the CglTwomir class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
CGLLIB_EXPORT
void CglTwomirUnitTest(const OsiSolverInterface * siP,
		       const std::string mpdDir);


#endif


