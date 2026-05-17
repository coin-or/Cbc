// LAST EDIT: 
//-----------------------------------------------------------------------------
// name: Mixed Integer Rounding Cut Generator
// authors: Joao Goncalves (jog7@lehigh.edu) 
//          Laszlo Ladanyi (ladanyi@us.ibm.com) 
// date: August 11, 2004 
//-----------------------------------------------------------------------------
// Copyright (C) 2004, International Business Machines Corporation and others. 
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef CglMixedIntegerRounding_H
#define CglMixedIntegerRounding_H

#include <iostream>
#include <fstream>
//#include <vector>

#include "CoinError.hpp"

#include "CglCutGenerator.hpp"

//=============================================================================

#ifndef CGL_DEBUG
#define CGL_DEBUG 0
#endif

//=============================================================================

// Class to store variable upper bounds (VUB)
class CGLLIB_EXPORT CglMixIntRoundVUB
{
  // Variable upper bounds have the form x_j <= a y_j, where x_j is
  // a continuous variable and y_j is an integer variable

protected:
  int    var_;            // The index of y_j
  double val_;            // The value of a 

public:
  // Default constructor
  CglMixIntRoundVUB() : var_(-1), val_(-1) {}

  // Copy constructor
  CglMixIntRoundVUB(const CglMixIntRoundVUB& source) { 
    var_ = source.var_; 
    val_ = source.val_; 
  } 

  // Assignment operator
  CglMixIntRoundVUB& operator=(const CglMixIntRoundVUB& rhs) { 
    if (this != &rhs) { 
      var_ = rhs.var_; 
      val_ = rhs.val_; 
    }
    return *this; 
  }

  // Destructor
  ~CglMixIntRoundVUB() {}

  // Query and set functions
  int    getVar() const          { return var_; }
  double getVal() const          { return val_; }
  void   setVar(const int v)     { var_ = v; }
  void   setVal(const double v)  { val_ = v; }
};

//=============================================================================

// Class to store variable lower bounds (VLB).
// It is the same as the class to store variable upper bounds
typedef CglMixIntRoundVUB CglMixIntRoundVLB;

//=============================================================================

/** Mixed Integer Rounding Cut Generator Class */

// Reference: 
//    Hugues Marchand and Laurence A. Wolsey
//    Aggregation and Mixed Integer Rounding to Solve MIPs
//    Operations Research, 49(3), May-June 2001.
//    Also published as CORE Dicusion Paper 9839, June 1998.

class CGLLIB_EXPORT CglMixedIntegerRounding : public CglCutGenerator {

  friend CGLLIB_EXPORT void CglMixedIntegerRoundingUnitTest(const OsiSolverInterface * siP,
					      const std::string mpdDir);


private:
  //---------------------------------------------------------------------------
  // Enumeration constants that describe the various types of rows
  enum RowType {
    // The row type of this row is NOT defined yet.
    ROW_UNDEFINED,
    /** After the row is flipped to 'L', the row has exactly two variables: 
	one is negative binary and the other is a continous, 
	and the RHS is zero.*/
    ROW_VARUB,
    /** After the row is flipped to 'L', the row has exactly two variables: 
	one is positive binary and the other is a continous, 
	and the RHS is zero.*/
    ROW_VARLB,
    /** The row sense is 'E', the row has exactly two variables: 
	one is binary and the other is a continous, and the RHS is zero.*/ 
    ROW_VAREQ,
    // The row contains continuous and integer variables;
    // the total number of variables is at least 2
    ROW_MIX,
    // The row contains only continuous variables
    ROW_CONT,
    // The row contains only integer variables
    ROW_INT,
    // The row contains other types of rows
    ROW_OTHER
  };


public:

  /**@name Generate Cuts */
  //@{
  /** Generate Mixed Integer Rounding cuts for the model data 
      contained in si. The generated cuts are inserted 
      in the collection of cuts cs. 
  */
  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());
  //@}

  //---------------------------------------------------------------------------
  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglMixedIntegerRounding ();

  /// Alternate Constructor 
  CglMixedIntegerRounding (const int maxaggr,
			   const bool multiply,
			   const int criterion,
			   const int preproc = -1);

  /// Copy constructor 
  CglMixedIntegerRounding (
    const CglMixedIntegerRounding &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglMixedIntegerRounding &
    operator=(
    const CglMixedIntegerRounding& rhs);
  
  /// Destructor 
  virtual
    ~CglMixedIntegerRounding ();
  /// This can be used to refresh any inforamtion
  virtual void refreshSolver(OsiSolverInterface * solver);
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  //@}

  //---------------------------------------------------------------------------
  /**@name Set and get methods */
  //@{
  /// Set MAXAGGR_
  inline void setMAXAGGR_ (int maxaggr) {
    if (maxaggr > 0) {
      MAXAGGR_ = maxaggr;
    }
    else {
      throw CoinError("Unallowable value. maxaggr must be > 0",
                      "gutsOfConstruct","CglMixedIntegerRounding");
    }
  }

  /// Get MAXAGGR_
  inline int getMAXAGGR_ () const { return MAXAGGR_; }

  /// Set MULTIPLY_
  inline void setMULTIPLY_ (bool multiply) { MULTIPLY_ = multiply; }

  /// Get MULTIPLY_
  inline bool getMULTIPLY_ () const { return MULTIPLY_; }

  /// Set CRITERION_
  inline void setCRITERION_ (int criterion) {
    if ((criterion >= 1) && (criterion <= 3)) {
      CRITERION_ = criterion;
    }
    else {
      throw CoinError("Unallowable value. criterion must be 1, 2 or 3",
                      "gutsOfConstruct","CglMixedIntegerRounding");
    }
  }

  /// Get CRITERION_
  inline int getCRITERION_ () const { return CRITERION_; }


  /// Set doPreproc
  void setDoPreproc(int value);
  /// Get doPreproc
  bool getDoPreproc() const;

  //@}

private:
  //--------------------------------------------------------------------------
  // Private member methods

  // Construct
  void gutsOfConstruct (const int maxaggr,
			const bool multiply,
			const int criterion,
			const int preproc);

  // Delete
  void gutsOfDelete();

  // Copy
  void gutsOfCopy (const CglMixedIntegerRounding& rhs);

  // Do preprocessing.
  // It determines the type of each row. It also identifies the variable
  // upper bounds and variable lower bounds.
  // It may change sense and RHS for ranged rows
  void mixIntRoundPreprocess(const OsiSolverInterface& si);

  // Determine the type of a given row.
  RowType determineRowType(const OsiSolverInterface& si,
			   const int rowLen, const int* ind, 
			   const double* coef, const char sense, 
			   const double rhs) const;

  // Generate MIR cuts
  void generateMirCuts( const OsiSolverInterface& si,
			const double* xlp,
			const double* colUpperBound,
			const double* colLowerBound,
			const CoinPackedMatrix& matrixByRow,
			const double* LHS,
			const double* coefByRow,
			const int* colInds,
			const CoinBigIndex* rowStarts,
			const int* rowLengths,
			//const CoinPackedMatrix& matrixByCol,
			const double* coefByCol,
			const int* rowInds,
			const CoinBigIndex* colStarts,
			const int* colLengths,
			OsiCuts& cs ) const;

  // Copy row selected to CoinPackedVector
  void copyRowSelected( const int iAggregate,
			const int rowSelected,
			std::set<int>& setRowsAggregated,
			int* listRowsAggregated,
			double* xlpExtra,
			const char sen,
			const double rhs,
			const double lhs,
			const CoinPackedMatrix& matrixByRow,
			CoinPackedVector& rowToAggregate,
			double& rhsToAggregate) const;

  // Select a row to aggregate
  bool selectRowToAggregate( const OsiSolverInterface& si,
			     const CoinPackedVector& rowAggregated,
			     const double* colUpperBound,
			     const double* colLowerBound,
			     const std::set<int>& setRowsAggregated,
			     const double* xlp, const double* coefByCol,
			     const int* rowInds, const CoinBigIndex* colStarts,
			     const int* colLengths,
			     int& rowSelected,
			     int& colSelected ) const;

  // Aggregation heuristic. 
  // Combines one or more rows of the original matrix 
  void aggregateRow( const int colSelected,
		     CoinPackedVector& rowToAggregate, double rhs,
		     CoinPackedVector& rowAggregated, 
		     double& rhsAggregated ) const;

  // Choose the bound substitution based on the criteria defined by the user
  inline bool isLowerSubst(const double inf, 
			   const double aj,
			   const double xlp, 
			   const double LB, 
			   const double UB) const;
    
  // Bound substitution heuristic
  bool boundSubstitution( const OsiSolverInterface& si,
			  const CoinPackedVector& rowAggregated,
			  const double* xlp,
			  const double* xlpExtra,
			  const double* colUpperBound,
			  const double* colLowerBound,
			  CoinPackedVector& mixedKnapsack,
			  double& rhsMixedKnapsack, double& sStar,
			  CoinPackedVector& contVariablesInS ) const;

  // c-MIR separation heuristic
  bool cMirSeparation ( const OsiSolverInterface& si,
			const CoinPackedMatrix& matrixByRow,
			const CoinPackedVector& rowAggregated,
			const int* listRowsAggregated,
			const char* sense, const double* RHS,
			//const double* coefByRow,
			//const int* colInds, const int* rowStarts,
			//const int* rowLengths,
			const double* xlp, const double sStar,
			const double* colUpperBound,
			const double* colLowerBound,
			const CoinPackedVector& mixedKnapsack,
			const double& rhsMixedKnapsack,
			const CoinPackedVector& contVariablesInS,
			OsiRowCut& flowCut ) const;

  // function to create one c-MIR inequality
  void cMirInequality( const int numInt, 
		       const double delta,
		       const double numeratorBeta,
		       const int *knapsackIndices,
		       const double* knapsackElements,
		       const double* xlp, 
		       const double sStar,	       
		       const double* colUpperBound,
		       const std::set<int>& setC,
		       CoinPackedVector& cMIR,
		       double& rhscMIR,
		       double& sCoef,
		       double& violation) const;

  // function to compute G
  inline double functionG( const double d, const double f ) const;

  // function to print statistics (used only in debug mode)
  void printStats(
			    std::ofstream & fout,
			    const bool hasCut,
			    const OsiSolverInterface& si,
			    const CoinPackedVector& rowAggregated,
			    const double& rhsAggregated, const double* xlp,
			    const double* xlpExtra,
			    const int* listRowsAggregated,
			    const int* listColsSelected,
			    const int level,
			    const double* colUpperBound,
			    const double* colLowerBound ) const;


private:
  //---------------------------------------------------------------------------
  // Private member data
  
  // Maximum number of rows to aggregate
  int MAXAGGR_;
  // Flag that indicates if an aggregated row is also multiplied by -1
  bool MULTIPLY_;
  // The criterion to use in the bound substitution
  int CRITERION_;
  // Tolerance used for numerical purposes
  double EPSILON_;
  /// There is no variable upper bound or variable lower bound defined
  int UNDEFINED_;
  // If violation of a cut is greater that this number, the cut is accepted
  double TOLERANCE_;
  /** Controls the preprocessing of the matrix to identify rows suitable for
      cut generation.<UL>
      <LI> -1: preprocess according to solver settings;
      <LI> 0: Do preprocessing only if it has not yet been done;
      <LI> 1: Do preprocessing.
      </UL>
      Default value: -1 **/
  int doPreproc_;
  // The number of rows of the problem.
  int numRows_;
  // The number columns of the problem.
  int numCols_;
  // Indicates whether preprocessing has been done.
  bool doneInitPre_;
  // The array of CglMixIntRoundVUBs.
  CglMixIntRoundVUB* vubs_;
  // The array of CglMixIntRoundVLBs.
  CglMixIntRoundVLB* vlbs_;
  // Array with the row types of the rows in the model.
  RowType* rowTypes_;
  // The indices of the rows of the initial matrix
  int* indRows_;
  // The number of rows of type ROW_MIX
  int numRowMix_;
  // The indices of the rows of type ROW_MIX
  int* indRowMix_;
  // The number of rows of type ROW_CONT
  int numRowCont_;
  // The indices of the rows of type ROW_CONT
  int* indRowCont_;
  // The number of rows of type ROW_INT
  int numRowInt_;
  // The indices of the rows of type ROW_INT
  int* indRowInt_;
  // The number of rows of type ROW_CONT that have at least one variable
  // with variable upper or lower bound
  int numRowContVB_;
  // The indices of the rows of type ROW_CONT that have at least one variable
  // with variable upper or lower bound
  int* indRowContVB_;
  // Sense of rows (modified if ranges)
  char * sense_;
  // RHS of rows (modified if ranges)
  double * RHS_;
  
};

//#############################################################################
// A function that tests the methods in the CglMixedIntegerRounding class. The
// only reason for it not to be a member method is that this way it doesn't
// have to be compiled into the library. And that's a gain, because the
// library should be compiled with optimization on, but this method should be
// compiled with debugging.
CGLLIB_EXPORT
void CglMixedIntegerRoundingUnitTest(const OsiSolverInterface * siP,
				     const std::string mpdDir);
  
#endif
