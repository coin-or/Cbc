// LAST EDIT: 
//-----------------------------------------------------------------------------
// Implementation of Residual Capacity Inequalities
// Francisco Barahona (barahon@us.ibm.com)
//  
// date: May 18 2006
//-----------------------------------------------------------------------------
// Copyright (C) 2004, International Business Machines Corporation and others. 
// All Rights Reserved.
// This code is published under the Eclipse Public License.

//#include <cmath>
//#include <cstdlib>
#include <cassert>

#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"

#include "CglResidualCapacity.hpp"
//#define CGL_DEBUG 1
//-----------------------------------------------------------------------------
// Generate Mixed Integer Rounding inequality
//------------------------------------------------------------------- 
void
CglResidualCapacity::generateCuts(const OsiSolverInterface& si,
				      OsiCuts& cs,
				  const CglTreeInfo /*info*/)
{

  // If the LP or integer presolve is used, then need to redo preprocessing
  // everytime this function is called. Otherwise, just do once.
  bool preInit = false;
  bool preReso = false;
  si.getHintParam(OsiDoPresolveInInitial, preInit);
  si.getHintParam(OsiDoPresolveInResolve, preReso);
  si.getColType(true);
  if (preInit == false &&  preReso == false &&
      doPreproc_ == -1 ) { // Do once
    if (doneInitPre_ == false) {   
      resCapPreprocess(si);
      doneInitPre_ = true;
    }
  }
  else  
      if ( doPreproc_ == 1 ){ // Do everytime       
	  resCapPreprocess(si);
	  doneInitPre_ = true;
      } else
	if (doneInitPre_ == false) {   
	  resCapPreprocess(si);
	  doneInitPre_ = true;
	}  


  const double* xlp        = si.getColSolution();  // LP solution
  const double* colUpperBound = si.getColUpper();  // vector of upper bounds
  const double* colLowerBound = si.getColLower();  // vector of lower bounds

  // get matrix by row
  const CoinPackedMatrix & tempMatrixByRow = *si.getMatrixByRow();
  CoinPackedMatrix matrixByRow;
  matrixByRow.submatrixOf(tempMatrixByRow, numRows_, indRows_);

  const double* LHS        = si.getRowActivity();
  const double* coefByRow  = matrixByRow.getElements();
  const int* colInds       = matrixByRow.getIndices();
  const CoinBigIndex* rowStarts     = matrixByRow.getVectorStarts();
  const int* rowLengths    = matrixByRow.getVectorLengths();


  generateResCapCuts(si, xlp, colUpperBound, colLowerBound,
		     matrixByRow, LHS, coefByRow,
		     colInds, rowStarts, rowLengths, 
		     cs);
}

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglResidualCapacity::CglResidualCapacity ()
    :
    CglCutGenerator()
{ 
    gutsOfConstruct(1.0e-6);
}


//-------------------------------------------------------------------
// Alternate Constructor 
//-------------------------------------------------------------------
CglResidualCapacity::CglResidualCapacity (const double epsilon)
  :
  CglCutGenerator()
{ 
  gutsOfConstruct(epsilon);
}


//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglResidualCapacity::CglResidualCapacity ( 
				 const CglResidualCapacity & rhs)
  :
  CglCutGenerator(rhs)
{ 
  gutsOfCopy(rhs);
}


//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglResidualCapacity::clone() const
{
  return new CglResidualCapacity(*this);
}

//------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglResidualCapacity &
CglResidualCapacity::operator=(const CglResidualCapacity& rhs)
{
  if (this != &rhs) {
    gutsOfDelete();
    CglCutGenerator::operator=(rhs);
    gutsOfCopy(rhs);
  }
  return *this;
}


//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------  
CglResidualCapacity::~CglResidualCapacity ()
{
  gutsOfDelete();
}

//-------------------------------------------------------------------
// Construct
//-------------------------------------------------------------------  
void
CglResidualCapacity::gutsOfConstruct (const double epsilon)
{
  
    EPSILON_ = epsilon;
    TOLERANCE_ = 1.0e-4;
    doPreproc_ = -1;
    numRows_ = 0;
    numCols_ = 0;
    doneInitPre_ = false;
    rowTypes_ = 0;
    indRows_ = 0;
    sense_=NULL;
    RHS_=NULL; 
    numRowL_ = 0;
    indRowL_ = 0;
    numRowG_ = 0;
    indRowG_ = 0;
}

//-------------------------------------------------------------------
// Delete
//-------------------------------------------------------------------  
void
CglResidualCapacity::gutsOfDelete ()
{
  if (rowTypes_ != 0) { delete [] rowTypes_; rowTypes_ = 0; } 
  if (indRows_ != 0) { delete [] indRows_; indRows_ = 0; }
  if (indRowL_ != 0) { delete [] indRowL_; indRowL_ = 0; }
  if (indRowG_ != 0) { delete [] indRowG_; indRowG_ = 0; }
  if (sense_ !=NULL) { delete [] sense_; sense_=NULL;}
  if (RHS_ !=NULL) { delete [] RHS_; RHS_=NULL;}
}

//-------------------------------------------------------------------
// Copy
//-------------------------------------------------------------------  
void
CglResidualCapacity::gutsOfCopy (const CglResidualCapacity& rhs)
{
  EPSILON_ = rhs.EPSILON_;
  TOLERANCE_ = rhs.TOLERANCE_;
  doPreproc_ = rhs.doPreproc_;
  numRows_ = rhs.numRows_;
  numCols_ = rhs.numCols_;
  doneInitPre_ = rhs.doneInitPre_;
  numRowL_ = rhs.numRowL_;
  numRowG_ = rhs.numRowG_;


  if (numRows_ > 0) {
    rowTypes_ = new RowType [numRows_];
    CoinDisjointCopyN(rhs.rowTypes_, numRows_, rowTypes_);
    indRows_ = new int [numRows_];
    CoinDisjointCopyN(rhs.indRows_, numRows_, indRows_);
    sense_ = CoinCopyOfArray(rhs.sense_,numRows_);
    RHS_ = CoinCopyOfArray(rhs.RHS_,numRows_);
  }
  else {
    rowTypes_ = 0;
    indRows_ = 0;
    sense_=NULL;
    RHS_=NULL;
  }

  if (numRowL_ > 0) {
    indRowL_ = new int [numRowL_];
    CoinDisjointCopyN(rhs.indRowL_, numRowL_, indRowL_);
  }
  else {
    indRowL_ = 0;
  }

  if (numRowG_ > 0) {
    indRowG_ = new int [numRowG_];
    CoinDisjointCopyN(rhs.indRowG_, numRowG_, indRowG_);
  }
  else {
    indRowG_ = 0;
  }

}

//-------------------------------------------------------------------
// Do preprocessing
// It determines the type of each row.
//-------------------------------------------------------------------  
void 
CglResidualCapacity::
resCapPreprocess(const OsiSolverInterface& si)
{
    // get matrix stored by row
    const CoinPackedMatrix & matrixByRow = *si.getMatrixByRow();
    numRows_ = si.getNumRows();
    numCols_ = si.getNumCols();
    const double* coefByRow  = matrixByRow.getElements();
    const int* colInds       = matrixByRow.getIndices();
    const CoinBigIndex* rowStarts     = matrixByRow.getVectorStarts();
    const int* rowLengths    = matrixByRow.getVectorLengths();
    const double * colLowerBound = si.getColLower();
    const double * colUpperBound = si.getColUpper();
    // Get copies of sense and RHS so we can modify if ranges
    if (sense_) {
	delete [] sense_;
	delete [] RHS_;
    }
    sense_ = CoinCopyOfArray(si.getRowSense(),numRows_);
    RHS_  = CoinCopyOfArray(si.getRightHandSide(),numRows_);
    
    if (rowTypes_ != 0) {
	delete [] rowTypes_; rowTypes_ = 0;
    }
    rowTypes_ = new RowType [numRows_];     // Destructor will free memory
    
    // Summarize the row type infomation.
    //int numOTHER   = 0;
    int numL       = 0;
    int numG       = 0;
    int numB       = 0;
    
    int iRow;
    const double* rowActivity        = si.getRowActivity();
    const double* rowLower        = si.getRowLower();
    const double* rowUpper        = si.getRowUpper();
    for (iRow = 0; iRow < numRows_; ++iRow) {
	// If range then choose which to use
	if (sense_[iRow]=='R') {
	    if (rowActivity[iRow]-rowLower[iRow]<
		rowUpper[iRow]-rowActivity[iRow]) {
		// treat as G row
		RHS_[iRow]=rowLower[iRow];
		sense_[iRow]='G';
	    } else {
		// treat as L row
		RHS_[iRow]=rowUpper[iRow];
		sense_[iRow]='L';
	    }
	}
	// get the type of a row
	const RowType rowType = 
	    determineRowType(si, rowLengths[iRow], colInds+rowStarts[iRow],
			     coefByRow+rowStarts[iRow], sense_[iRow], RHS_[iRow],
			     colLowerBound, colUpperBound);
	// store the type of the current row
	rowTypes_[iRow] = rowType;
	
	// Summarize information about row types
	switch(rowType) {
	case  ROW_OTHER:
	    //++numOTHER;
	    break;
	case  ROW_L:
	    ++numL; 
	    break;
	case  ROW_G:
	    ++numG; 
	    break;
	case ROW_BOTH:
	    ++numB;
	    break;
	default:
	    throw CoinError("Unknown row type", "ResidualCapacityPreprocess",
			    "CglResidualCapacity");
	}
    }
    
    // allocate memory for vector of indices of all rows
    if (indRows_ != 0) { delete [] indRows_; indRows_ = 0; }
    if (numRows_ > 0)
	indRows_ = new int [numRows_];     // Destructor will free memory
    // allocate memory for vector of indices of rows of type ROW_L and ROW_BOTH
    numRowL_ = numL + numB;
    if (indRowL_ != 0) { delete [] indRowL_; indRowL_ = 0; }
    if (numRowL_ > 0)
	indRowL_ = new int [numRowL_];     // Destructor will free memory
    // allocate memory for vector of indices of rows of type ROW_G and ROW_BOTH
    numRowG_ = numG + numB;
    if (indRowG_ != 0) { delete [] indRowG_; indRowG_ = 0; }
    if (numRowG_ > 0)
	indRowG_ = new int [numRowG_];     // Destructor will free memory
    
    
#if CGL_DEBUG
    std::cout << "The num of rows = "  << numRows_        << std::endl;
    std::cout << "Summary of Row Type" << std::endl;
    std::cout << "numL          = " << numL        << std::endl;
    std::cout << "numG          = " << numG        << std::endl;
#endif
    
    
    int countL = 0;
    int countG = 0;
    for ( iRow = 0; iRow < numRows_; ++iRow) {
	
	RowType rowType = rowTypes_[iRow];
	
	// fill the vector indRows_ with the indices of all rows
	indRows_[iRow] = iRow;
	
	// fill the vector indRowL_ with the indices of the rows of type ROW_L and ROW_BOTH
	if (rowType == ROW_L || rowType == ROW_BOTH) {
	    indRowL_[countL] = iRow;
	    countL++;
	}
	// fill the vector indRowG_ with the indices of rows of type ROW_G and ROW_BOTH
	if (rowType == ROW_G || rowType == ROW_BOTH) {
	    indRowG_[countG] = iRow;
	    countG++;
	}
    }
    
}

//-------------------------------------------------------------------
// Determine the type of a given row 
//-------------------------------------------------------------------
CglResidualCapacity::RowType
CglResidualCapacity::determineRowType(const OsiSolverInterface& si,
				      const int rowLen, const int* ind, 
				      const double* coef, const char sense, 
				      const double rhs,
				      const double* colLowerBound,
				      const double* colUpperBound) const
{
    if (rowLen == 0) 
	return ROW_OTHER;
    RowType rowType = ROW_OTHER;
    double *negCoef;
    bool flagL, flagG, flag1, flag2;
    switch (sense) {
    case 'L':
	flagL=treatAsLessThan(si, rowLen, ind, coef, rhs, colLowerBound,
				  colUpperBound);
	if ( flagL ) rowType=ROW_L;
	break;
    case 'G':
	negCoef = new double[rowLen];
	for ( int i=0; i < rowLen; ++i )
	    negCoef[i]=-coef[i];
	flagG=treatAsLessThan(si, rowLen, ind, negCoef, -rhs, colLowerBound,
				   colUpperBound);
	if ( flagG ) rowType=ROW_G;
	delete [] negCoef;
	break;
    case 'E':
	flag1=treatAsLessThan(si, rowLen, ind, coef, rhs, colLowerBound,
				   colUpperBound);
	
	negCoef = new double[rowLen];
	for ( int i=0; i < rowLen; ++i )
	    negCoef[i]=-coef[i];
	flag2=treatAsLessThan(si, rowLen, ind, negCoef, -rhs, colLowerBound,
				   colUpperBound);
	delete [] negCoef;
	if ( flag1 && !flag2 ) rowType=ROW_L;
	if ( !flag1 && flag2 ) rowType=ROW_G;
	if ( flag1 && flag2  ) rowType=ROW_BOTH;
	break;
    default:
      // presumably range throw CoinError("Unknown sense", "determineRowType",
      //		"CglResidualCapacity");
      break;
    }
    return rowType;
}
//--------------------------------------------
// determine if an ineq of type <= is a good candidate
//--------------------------------------------
	
bool
CglResidualCapacity::treatAsLessThan(const OsiSolverInterface& si,
				     const int rowLen, const int* ind, 
				     const double* coef,
				     const double /*rhs*/,
				     const double* colLowerBound,
				     const double* colUpperBound) const
{
    bool intFound=false;
    bool contFound=false;
    bool goodIneq=true;
    double intCoef=-1;
    const char * intVar = si.getColType();
    // look for a_1 c_1 +   + a_k c_k  - d z_1 -   - d z_p <= b
    // where c_i continuous, z_j integer
    for ( int i = 0; i < rowLen; ++i ) {
	if ( coef[i]  > EPSILON_ || !intVar[ind[i]] ) {
	    if ( colLowerBound[ind[i]] < -EPSILON_ || colUpperBound[ind[i]] > 1.e10 ){
		// cont var with too big bounds
		goodIneq=false;
		break;
	    } else
		contFound=true;
	} else
	    if ( !intFound &&  coef[i] < -EPSILON_  && intVar[ind[i]] ){
		intFound=true;
		intCoef=coef[i];
		continue;
	    } else
		if ( intFound && coef[i] < -EPSILON_  && intVar[ind[i]] &&
		     fabs( coef[i] - intCoef ) > EPSILON_ ){
		    goodIneq=false;
		    break;
		}
    }
    if ( contFound && intFound && goodIneq ) return true;
    else return false;
}

//-------------------------------------------------------------------
// Generate Residual capacity cuts
//-------------------------------------------------------------------
void
CglResidualCapacity::generateResCapCuts( 
				     const OsiSolverInterface& si,
				     const double* xlp,
				     const double* colUpperBound,
				     const double* colLowerBound,
				     const CoinPackedMatrix& /*matrixByRow*/,
				     const double* /*LHS*/,
				     const double* coefByRow,
				     const int* colInds,
				     const CoinBigIndex* rowStarts,
				     const int* rowLengths,
				     OsiCuts& cs ) const
{
    
#if CGL_DEBUG
    // OPEN FILE
    std::ofstream fout("stats.dat");
#endif
    
    for (int iRow = 0; iRow < numRowL_; ++iRow) {
	int rowToUse=indRowL_[iRow];
	OsiRowCut resCapCut;
	// Find a most violated residual capacity ineq
	bool hasCut = resCapSeparation(si, rowLengths[rowToUse],
				       colInds+rowStarts[rowToUse],
				       coefByRow+rowStarts[rowToUse],
				       RHS_[rowToUse],
				       xlp, colUpperBound, colLowerBound, 
				       resCapCut);
	
	// if a cut was found, insert it into cs
	if (hasCut)  {
#if CGL_DEBUG
	    std::cout << "Res. cap. cut generated " << std::endl;
#endif
	    cs.insertIfNotDuplicateAndClean(resCapCut,71);
	}
    }
    
    for (int iRow = 0; iRow < numRowG_; ++iRow) {
	int rowToUse=indRowG_[iRow];
	OsiRowCut resCapCut;
	const int rowLen=rowLengths[rowToUse];
	double *negCoef= new double[rowLen];
	const CoinBigIndex rStart=rowStarts[rowToUse];
	for ( int i=0; i < rowLen; ++i )
	    negCoef[i]=-coefByRow[rStart+i];
	// Find a most violated residual capacity ineq
	bool hasCut = resCapSeparation(si, rowLengths[rowToUse],
				       colInds+rowStarts[rowToUse],
				       negCoef,
				       -RHS_[rowToUse],
				       xlp, colUpperBound, colLowerBound, 
				       resCapCut);
	delete [] negCoef;
	// if a cut was found, insert it into cs
	if (hasCut)  {
#if CGL_DEBUG
	    std::cout << "Res. cap. cut generated " << std::endl;
#endif
	    cs.insertIfNotDuplicateAndClean(resCapCut,72);
	}
    }

#if CGL_DEBUG
    // CLOSE FILE
    fout.close();
#endif
    
    return;
}

//-------------------------------------------------------------------
// separation algorithm
//-------------------------------------------------------------------
bool
CglResidualCapacity::resCapSeparation(const OsiSolverInterface& si,
				      const int rowLen, const int* ind, 
				      const double* coef,
				      const double rhs,
				      const double *xlp,  
				      const double* colUpperBound,
				      const double* /*colLowerBound*/,
				      OsiRowCut& resCapCut) const
{ 
    // process original row to create row in canonical form
    std::vector<int> positionIntVar;
    double ybar=0.0;
    double *xbar;
    double intCoef=-1;
    double *newRowCoef;
    int *positionContVar;
    double newRowRHS;
    int contCount=0;
    const char * intVar = si.getColType();
    for ( int i = 0; i < rowLen; ++i ) {
      if ( coef[i] < -EPSILON_ && intVar[ind[i]] ){
	intCoef=-coef[i];
	ybar+=xlp[ind[i]];
	positionIntVar.push_back(i);
      }
      else 
	++contCount;
    }
    xbar = new double [contCount];
    newRowCoef = new double [contCount];
    positionContVar = new int [contCount];
    contCount=0;
    newRowRHS=rhs;
    for ( int i = 0; i < rowLen; ++i ) {
      if ( coef[i] > EPSILON_ || !intVar[ind[i]] ){
	newRowCoef[contCount]=coef[i]*colUpperBound[ind[i]];
	xbar[contCount]=xlp[ind[i]]/colUpperBound[ind[i]];
	if ( newRowCoef[contCount] < -EPSILON_ ){ // complement
	  newRowCoef[contCount] = -newRowCoef[contCount];
	  xbar[contCount] = 1.0 - xbar[contCount];
	  newRowRHS+= newRowCoef[contCount];
	}
	positionContVar[contCount++]=i;
      }
    } 
    // now separate
    std::vector<int> setSbar;
    const double lambda = ybar - floor(ybar);
    double sumCoef=0.0;
    for ( int i = 0; i < contCount; ++i )
	if ( xbar[i] > lambda ){
	    setSbar.push_back(i);
	    sumCoef+=newRowCoef[i];
	}
    const int sSize = static_cast<int>(setSbar.size());
    bool generated;
    if ( sSize == 0 ) generated=false; // no cut
    else {
	// generate cut
	const double mu= ceil( (sumCoef - newRowRHS)/intCoef );
	double r = sumCoef - newRowRHS - intCoef * floor( (sumCoef - newRowRHS)/intCoef );
	const int numInt = static_cast<int>(positionIntVar.size());
	const int cutLen = sSize + numInt;
	int* cutInd = new int [cutLen];
	double* cutCoef = new double [cutLen];
	double violation=0.0;
	double complCoef=0.0;
	// load continuous variables
	for ( int i = 0; i < sSize; ++i ){
	    const int newRowPosition=setSbar[i];
	    const int originalRowPosition=positionContVar[newRowPosition];
	    cutInd[i]=ind[originalRowPosition];
	    cutCoef[i]=coef[originalRowPosition];
	    if ( cutCoef[i] < -EPSILON_ ) 
		complCoef+= cutCoef[i]*colUpperBound[ind[originalRowPosition]];
	    violation+=cutCoef[i]*xlp[ind[originalRowPosition]];
	}
	// load integer variables
	for ( int i = 0; i < numInt; ++i ){
	    const int originalRowPosition=positionIntVar[i];
	    cutInd[i+sSize]=ind[originalRowPosition];
	    cutCoef[i+sSize]= - r;
	    violation+=cutCoef[i+sSize]*xlp[ind[originalRowPosition]];
	}
	double cutRHS=(sumCoef - r * mu) + complCoef;
	violation-=cutRHS;
	if ( violation > TOLERANCE_ ){
	    resCapCut.setRow(cutLen, cutInd, cutCoef);
	    resCapCut.setLb(-1.0 * si.getInfinity());
	    resCapCut.setUb(cutRHS);
	    resCapCut.setEffectiveness(violation);
	    generated=true;
#if 0
	    std::cout << "coef ";
	    for(int i=0; i<cutLen; ++i)
		std::cout << cutCoef[i] << " ";
	    std::cout << std::endl << " bounds ";
	    for(int i=0; i<cutLen; ++i)
		std::cout << colUpperBound[cutInd[i]] << " ";
	    std::cout << std::endl << " rhs " << cutRHS << std::endl;
#endif	    
	}
	else
	    generated=false;
	delete [] cutCoef;
	delete [] cutInd;
    }
    // free memory	
    delete [] positionContVar;
    delete [] newRowCoef; 
    delete [] xbar;
    return generated;
}



// This can be used to refresh preprocessing
void 
CglResidualCapacity::refreshPrep()
{
  doneInitPre_ = false;
}

//
void CglResidualCapacity::setEpsilon(double value)
{
    EPSILON_ = value;  
}
double CglResidualCapacity::getEpsilon() const
{
    return EPSILON_;
}
//
void CglResidualCapacity::setTolerance(double value)
{
    TOLERANCE_ = value;
}
double CglResidualCapacity::getTolerance() const
{
    return TOLERANCE_;
}
//
void CglResidualCapacity::setDoPreproc(int value)
{
    if ( value != -1 && value != 0 && value != 1 )
	throw CoinError("setDoPrepoc", "invalid value",
			"CglResidualCapacity");
    else
	doPreproc_ = value;  
}
bool CglResidualCapacity::getDoPreproc() const
{
    return (doPreproc_ != 0);
}
