//-----------------------------------------------------------------------------
// name:     Cgl Lifted Simple Generalized Flow Cover Cut Generator
// author:   Yan Xu                email: yan.xu@sas.com
//           Jeff Linderoth        email: jtl3@lehigh.edu
//           Martin Savelsberg     email: martin.savelsbergh@isye.gatech.edu
// date:     05/01/2003
// comments: please scan this file for '???' and read the comments
//-----------------------------------------------------------------------------
// Copyright (C) 2003, Yan Xu, Jeff Linderoth, Martin Savelsberg and others. 
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef CglFlowCover_H
#define CglFlowCover_H

#include <iostream>

#include "CoinError.hpp"

#include "CglCutGenerator.hpp"

//=============================================================================

//=============================================================================

/** This enumerative constant describes the various col types.*/
enum CglFlowColType {
    /** The column(variable) is a negative binary variable.*/
    CGLFLOW_COL_BINNEG  = -2,
    /** The column is a negative continous variable.*/
    CGLFLOW_COL_CONTNEG,
    /** The column is a positive continous variable.*/
    CGLFLOW_COL_CONTPOS =  1,
    /** The column is a positive binary variable.*/
    CGLFLOW_COL_BINPOS
};

enum CglFlowColStatus{
};

/** This enumerative constant describes the various stati of vars in 
    a cut or not.*/
enum CglFlowColCut{
    /** The column is NOT in cover.*/
    CGLFLOW_COL_OUTCUT = 0,
    /** The column is in cover now. */
    CGLFLOW_COL_INCUT,
    /** The column is decided to be in cover. */
    CGLFLOW_COL_INCUTDONE,
    /** The column is in L-. */
    CGLFLOW_COL_INLMIN,
    /** The column is decided to be in L-. */
    CGLFLOW_COL_INLMINDONE,
    /** The column is in L--.*/
    CGLFLOW_COL_INLMINMIN,
    /** This enumerative constant describes the various stati of vars in 
                   determining the cover.*/
    /** The column is a prime candidate. */
    CGLFLOW_COL_PRIME,
    /** The column is a secondary candidate. */
    CGLFLOW_COL_SECONDARY
};

/** This enumerative constant describes the various row types.*/
enum CglFlowRowType {
    /** The row type of this row is NOT defined yet.*/
    CGLFLOW_ROW_UNDEFINED,
    /** After the row is flipped to 'L', the row has exactly two variables: 
	one is negative binary and the other is continous, and the RHS 
	is zero.*/
    CGLFLOW_ROW_VARUB,
    /** After the row is flipped to 'L', the row has exactlytwo variables: 
	one is positive binary and the other is continous, and the RHS 
	is zero.*/
    CGLFLOW_ROW_VARLB,
    /** The row sense is 'E', the row has exactly two variables: 
	one is binary and the other is a continous, and the RHS is zero.*/ 
    CGLFLOW_ROW_VAREQ,
    /** Rows can not be classfied into other types and the row sense 
	is NOT 'E'.*/
    CGLFLOW_ROW_MIXUB,
    /** Rows can not be classfied into other types and the row sense is 'E'.*/
    CGLFLOW_ROW_MIXEQ,
    /** All variables are NOT binary and the row sense is NOT 'E'. */
    CGLFLOW_ROW_NOBINUB,
    /** All variables are NOT binary and the row sense is 'E'. */
    CGLFLOW_ROW_NOBINEQ,
    /** The row has one binary and 2 or more other types of variables and 
	the row sense is NOT 'E'. */
    CGLFLOW_ROW_SUMVARUB,
    /** The row has one binary and 2 or more other types of variables and 
	the row sense is 'E'. */
    CGLFLOW_ROW_SUMVAREQ,
    /** All variables are binary. */
    CGLFLOW_ROW_UNINTERSTED
};

//=============================================================================

/** Variable upper bound class. */
class CGLLIB_EXPORT CglFlowVUB
{
protected:
    int    varInd_;            /** The index of the associated 0-1 variable.*/
    double upper_;             /** The Value of the associated upper bound.*/ 

public:
    CglFlowVUB() : varInd_(-1), upper_(-1) {}
    
    CglFlowVUB(const CglFlowVUB& source) { 
	varInd_= source.varInd_; 
	upper_ = source.upper_; 
    } 
    
    CglFlowVUB& operator=(const CglFlowVUB& rhs) { 
	if (this == &rhs) 
	    return *this;
	varInd_= rhs.varInd_; 
	upper_ = rhs.upper_; 
	return *this; 
  }
    
    /**@name Query and set functions for associated 0-1 variable index 
       and value.
    */ 
    //@{  
    inline int    getVar() const          { return varInd_; }
    inline double getVal() const          { return upper_; }
    inline void   setVar(const int v)     { varInd_ = v; }
    inline void   setVal(const double v)  { upper_ = v; }
    //@}
};

//=============================================================================

/** Variable lower bound class, which is the same as vub. */
typedef CglFlowVUB CglFlowVLB;

/** Overloaded operator<< for printing VUB and VLB.*/
std::ostream& operator<<( std::ostream& os, const CglFlowVUB &v );

//=============================================================================

/** 
 *  Lifed Simple Generalized Flow Cover Cut Generator Class. 
 */
class CGLLIB_EXPORT CglFlowCover : public CglCutGenerator {
    friend CGLLIB_EXPORT void CglFlowCoverUnitTest(const OsiSolverInterface * siP,
				     const std::string mpdDir );
    
public:
    
    /** 
     *  Do the following tasks:
     *  <ul>
     *  <li> classify row types 
     *  <li> indentify vubs and vlbs
     *  </ul>
     *  This function is called by 
     *  <CODE>generateCuts(const OsiSolverInterface & si, OsiCuts & cs)</CODE>.
   */
    void flowPreprocess(const OsiSolverInterface& si);

    /**@name Generate Cuts */
    //@{
    /** Generate Lifed Simple Generalized flow cover cuts for the model data 
	contained in si. The generated cuts are inserted into and returned 
	in the collection of cuts cs. 
    */
    virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			      const CglTreeInfo info = CglTreeInfo());
    //@}

    /**@name Functions to query and set maximum number of cuts can be 
       generated. */
    //@{
    inline int getMaxNumCuts() const { return maxNumCuts_; }
    inline void setMaxNumCuts(int mc) { maxNumCuts_ = mc; }
    //@}
  
    /**@name Functions to query and set the number of cuts have been
       generated. */
    //@{
    inline int getNumFlowCuts() { return numFlowCuts_; }
    inline void setNumFlowCuts(int fc) { numFlowCuts_ = fc; }
    inline void incNumFlowCuts(int fc = 1) { numFlowCuts_ += fc; } 
    //@}

    //-------------------------------------------------------------------------
    /**@name Constructors and destructors */
    //@{
    /// Default constructor 
    CglFlowCover ();

    /// Copy constructor 
    CglFlowCover (
	const CglFlowCover &);

    /// Clone
    virtual CglCutGenerator * clone() const;

    /// Assignment operator 
    CglFlowCover &
    operator=(
	const CglFlowCover& rhs);
    
    /// Destructor 
    virtual
    ~CglFlowCover ();
    /// Create C++ lines to get to current state
    virtual std::string generateCpp( FILE * fp);
    //@}

private:
    //-------------------------------------------------------------------------
    // Private member functions

    /** Based a given row, a LP solution and other model data, this function
	tries to generate a violated lifted simple generalized flow cover. 
    */
    bool generateOneFlowCut( const OsiSolverInterface & si, 
			     const int rowLen,
			     int* ind,
			     double* coef,
			     char sense,
			     double rhs,
			     OsiRowCut& flowCut,
			     double& violation );


    /** Transform a row from ">=" to "<=", and vice versa. */
    void flipRow(int rowLen, double* coef, double& rhs) const;

    /** Transform a row from ">=" to "<=", and vice versa. Have 'sense'. */
    void flipRow(int rowLen, double* coef, char& sen, double& rhs) const;

    /** Determine the type of a given row. */
    CglFlowRowType determineOneRowType(const OsiSolverInterface& si,
				       int rowLen, int* ind, 
				       double* coef, char sen, 
				       double rhs) const;
    /** Lift functions */
    void liftMinus(double &movement, /* Output */ 
		   int t,
		   int r,
		   double z,
		   double dPrimePrime, 
		   double lambda,
		    double ml,
		   double *M,
		   double *rho) const;

    bool liftPlus(double &alpha, 
		 double &beta,
		 int r,
		 double m_j, 
		 double lambda,
		 double y_j,
		 double x_j,
		 double dPrimePrime,
		 double *M) const;
    

    //-------------------------------------------------------------------------
    //**@name Query and set the row type of a givne row. */
    //@{
    inline const CglFlowRowType* getRowTypes() const 
	{ return rowTypes_; }
    inline CglFlowRowType getRowType(const int i) const 
	{ return rowTypes_[i]; }
    /** Set rowtypes, take over the ownership. */
    inline void setRowTypes(CglFlowRowType* rt) 
	{ rowTypes_ = rt; rt = 0; }  
    inline void setRowTypes(const CglFlowRowType rt, const int i) {
	if (rowTypes_ != 0) 
	    rowTypes_[i] = rt;
	else {
	    std::cout << "ERROR: Should allocate memory for rowType_ before "
		      << "using it " << std::endl;
	    throw CoinError("Forgot to allocate memory for rowType_", 
			    "setRowType", "CglFlowCover");
	}
    }
    //@}
    
    //-------------------------------------------------------------------------
    //**@name Query and set vubs. */
    //@{
    inline const CglFlowVUB* getVubs() const          { return vubs_; }
    inline const CglFlowVUB& getVubs(int i) const     { return vubs_[i]; }
    /** Set CglFlowVUBs,take over the ownership. */
    inline void setVubs(CglFlowVUB* vubs) { vubs_ = vubs; vubs = 0; }
    inline void setVubs(const CglFlowVUB& vub, int i) { 
	if (vubs_ != 0) 
	    vubs_[i] = vub;
	else {
	    std::cout << "ERROR: Should allocate memory for vubs_ before "
		      << "using it " << std::endl;
	    throw CoinError("Forgot to allocate memory for vubs_", "setVubs",
			    "CglFlowCover");
	}
    }
    inline void printVubs(std::ostream& os) const {
	for (int i = 0; i < numCols_; ++i) {
	    os << "ix: " << i << ", " << vubs_[i];
	}
    }
    //@}

    //-------------------------------------------------------------------------
    //**@name Query and set vlbs. */
    //@{
    inline const CglFlowVLB* getVlbs() const          { return vlbs_; }
    inline const CglFlowVLB& getVlbs(int i) const     { return vlbs_[i]; }
    /** Set CglFlowVLBs,take over the ownership. */
    inline void setVlbs(CglFlowVLB* vlbs)          { vlbs_ = vlbs; vlbs = 0; }
    inline void setVlbs(const CglFlowVLB& vlb, int i) { 
	if (vlbs_ != 0) 
	    vlbs_[i] = vlb;
	else {
	    std::cout << "ERROR: Should allocate memory for vlbs_ before "
		      << "using it " << std::endl;
	    throw CoinError("Forgot to allocate memory for vlbs_", "setVlbs",
			    "CglFlowCover");
	}
    }
    //@}

private:
    //------------------------------------------------------------------------
    // Private member data
    
    /** The maximum number of flow cuts to be generated. Default is 1000. */
    int maxNumCuts_;
    /** Tolerance used for numerical purpose. */
    double EPSILON_;
    /** The variable upper bound of a flow is not indentified yet.*/
    int UNDEFINED_;
    /** Very large number. */
    double INFTY_;
    /** If violation of a cut is greater that this number, the cut is useful.*/
    double TOLERANCE_;
    /** First time preprocessing */
    bool firstProcess_;
    /** The number rows of the problem.*/
    int numRows_;
    /** The number columns of the problem.*/
    int numCols_;
    /** The number flow cuts found.*/
    int numFlowCuts_;
    /** Indicate whether initial flow preprecessing has been done. */
    bool doneInitPre_;
    /** The array of CglFlowVUBs. */
    CglFlowVUB* vubs_;
    /** The array of CglFlowVLBs. */
    CglFlowVLB* vlbs_;
    /** CglFlowRowType of the rows in model. */
    CglFlowRowType* rowTypes_;
};

//#############################################################################
/** A function that tests the methods in the CglFlowCover class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
CGLLIB_EXPORT
void CglFlowCoverUnitTest(const OsiSolverInterface * siP,
			  const std::string mpdDir );

#endif
