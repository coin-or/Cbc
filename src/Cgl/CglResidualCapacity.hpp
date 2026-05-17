// LAST EDIT: 
//-----------------------------------------------------------------------------
// Implementation of Residual Capacity Inequalities
// Francisco Barahona (barahon@us.ibm.com)
//         
// date: May 18, 2006
//-----------------------------------------------------------------------------
// Copyright (C) 2004, International Business Machines Corporation and others. 
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef CglResidualCapacity_H
#define CglResidualCapacity_H

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




//=============================================================================

/** Residual Capacity Inequalities Cut Generator Class

 References: 
    T Magnanti, P Mirchandani, R Vachani,
    "The convex hull of two core capacitated network design problems,"
    Math Programming 60 (1993), 233-250.
    
    A Atamturk, D Rajan,
    "On splittable and unsplittable flow capacitated network design 
    arc-set polyhedra," Math Programming 92 (2002), 315-333. **/

class CGLLIB_EXPORT CglResidualCapacity : public CglCutGenerator {
    
    friend CGLLIB_EXPORT void CglResidualCapacityUnitTest(const OsiSolverInterface * siP,
					    const std::string mpdDir );
    
    
private:
    //---------------------------------------------------------------------------
    // Enumeration constants that describe the various types of rows
    enum RowType {
	/**   row of the type a_1 c_1 +  + a_k c_k - d z_1 -  - d z_p <= b,
	      where c_i are continuous variables and z_j are integer variables
	*/
	ROW_L,
	/**  row of the type -a_1 c_1 -  - a_k c_k + d z_1 +  + d z_p >= b,
	     where c_i are continuous variables and z_j are integer variables
	*/
	ROW_G,
	/** equation that can be treated as ROW_L and ROW_G
	 */
	ROW_BOTH,
	/** Other types of rows
	 */
	ROW_OTHER
    };
    
    
public:
    /**@name Get and Set Parameters */
    //@{
    /// Set Epsilon
    void setEpsilon(double value);
    /// Get Epsilon
    double getEpsilon() const;
    /// Set Tolerance
    void setTolerance(double value);
    /// Get Tolerance
    double getTolerance() const;
    /// Set doPreproc
    void setDoPreproc(int value);
    /// Get doPreproc
    bool getDoPreproc() const;
    //@}

    /**@name Generate Cuts */
    //@{
    /** Generate Residual Capacity cuts for the model data 
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
    CglResidualCapacity ();
    
    /// Alternate Constructor 
    CglResidualCapacity ( const double tolerance );
    
    /// Copy constructor 
    CglResidualCapacity (
			 const CglResidualCapacity &);
    
    /// Clone
    virtual CglCutGenerator * clone() const;
    
    /// Assignment operator 
    CglResidualCapacity &
    operator=(
	      const CglResidualCapacity& rhs);
    
    /// Destructor 
    virtual
    ~CglResidualCapacity ();
    /// This is to refresh preprocessing
    virtual void refreshPrep();
    //@}
    
    
    
private:
    //--------------------------------------------------------------------------
    // Private member methods
    
    // Construct
    void gutsOfConstruct ( const double tolerance);
    
    // Delete
    void gutsOfDelete();
    
    // Copy
    void gutsOfCopy (const CglResidualCapacity& rhs);
    
    // Do preprocessing.
    // It determines the type of each row. 
    // It may change sense and RHS for ranged rows
    void resCapPreprocess(const OsiSolverInterface& si);
    
    // Determine the type of a given row.
    RowType determineRowType(const OsiSolverInterface& si,
			     const int rowLen, const int* ind, 
			     const double* coef, const char sense, 
			     const double rhs,
			     const double* colLowerBound,
			     const double* colUpperBound) const;
    // helps the function above
    bool treatAsLessThan(const OsiSolverInterface& si,
			 const int rowLen, const int* ind, 
			 const double* coef,
			 const double rhs,
			 const double* colLowerBound,
			 const double* colUpperBound) const;
    
    // Generate Residual Capacity cuts
    void generateResCapCuts( const OsiSolverInterface& si,
			     const double* xlp,
			     const double* colUpperBound,
			     const double* colLowerBound,
			     const CoinPackedMatrix& matrixByRow,
			     const double* LHS,
			     const double* coefByRow,
			     const int* colInds,
			     const CoinBigIndex* rowStarts,
			     const int* rowLengths,
			     OsiCuts& cs ) const;
    

    // Residual Capacity separation 
    bool resCapSeparation(const OsiSolverInterface& si,
			  const int rowLen, const int* ind, 
			  const double* coef,
			  const double rhs,
			  const double *xlp,  
			  const double* colUpperBound,
			  const double* colLowerBound,
			  OsiRowCut& resCapCut) const;
    
      

private:
    //---------------------------------------------------------------------------
    // Private member data
    /** Tolerance used for numerical purposes, default value: 1.e-6 **/
    double EPSILON_;
    /** If violation of a cut is greater that this number, 
	the cut is accepted, default value: 1.e-4 **/
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
    // Array with the row types of the rows in the model.
    RowType* rowTypes_;
    // The indices of the rows of the initial matrix
    int* indRows_;
    // Sense of rows (modified if ranges)
    char * sense_;
    // RHS of rows (modified if ranges)
    double * RHS_;
    // The number of rows of type ROW_L
    int numRowL_;
    // The indices of the rows of type ROW_L
    int* indRowL_;
    // The number of rows of type ROW_G
    int numRowG_;
    // The indices of the rows of type ROW_G
    int* indRowG_;
};

//#############################################################################
/** A function that tests the methods in the CglResidualCapacity class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
CGLLIB_EXPORT
void CglResidualCapacityUnitTest(const OsiSolverInterface * siP,
				 const std::string mpdDir);

  
#endif
