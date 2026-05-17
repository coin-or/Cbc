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

//#include <cmath>
//#include <cstdlib>
#include <cassert>

#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"

#include "CglMixedIntegerRounding.hpp"
//#define CGL_DEBUG 1
//-----------------------------------------------------------------------------
// Generate Mixed Integer Rounding inequality
//------------------------------------------------------------------- 
void
CglMixedIntegerRounding::generateCuts(const OsiSolverInterface& si,
				      OsiCuts& cs,
				      const CglTreeInfo )
{

  // If the LP or integer presolve is used, then need to redo preprocessing
  // everytime this function is called. Otherwise, just do once.
  bool preInit = false;
  bool preReso = false;
  si.getHintParam(OsiDoPresolveInInitial, preInit);
  si.getHintParam(OsiDoPresolveInResolve, preReso);
  if (preInit == false &&  preReso == false && doPreproc_ == -1 ) { // Do once
    if (doneInitPre_ == false) {   
      mixIntRoundPreprocess(si);
      doneInitPre_ = true;
    }
  }
  else {
    if(doPreproc_ == 1){ // Do everytime       
      mixIntRoundPreprocess(si);
      doneInitPre_ = true;
    } 
    else {
      if (doneInitPre_ == false) {   
	mixIntRoundPreprocess(si);
	doneInitPre_ = true;
      }  
    }
  }

  const double* xlp        = si.getColSolution();  // LP solution
  const double* colUpperBound = si.getColUpper();  // vector of upper bounds
  const double* colLowerBound = si.getColLower();  // vector of lower bounds

  // get matrix by row
  const CoinPackedMatrix & tempMatrixByRow = *si.getMatrixByRow();
  CoinPackedMatrix matrixByRow;
  matrixByRow.submatrixOf(tempMatrixByRow, numRows_, indRows_);
  CoinPackedMatrix matrixByCol = matrixByRow;
  matrixByCol.reverseOrdering();
  //const CoinPackedMatrix & matrixByRow = *si.getMatrixByRow();
  const double* LHS        = si.getRowActivity();
  const double* coefByRow  = matrixByRow.getElements();
  const int* colInds       = matrixByRow.getIndices();
  const CoinBigIndex* rowStarts     = matrixByRow.getVectorStarts();
  const int* rowLengths    = matrixByRow.getVectorLengths();

  // get matrix by column
  //const CoinPackedMatrix & matrixByCol = *si.getMatrixByCol();
  const double* coefByCol  = matrixByCol.getElements();
  const int* rowInds       = matrixByCol.getIndices();
  const CoinBigIndex* colStarts     = matrixByCol.getVectorStarts();
  const int* colLengths    = matrixByCol.getVectorLengths();


  generateMirCuts(si, xlp, colUpperBound, colLowerBound,
		  matrixByRow, LHS, coefByRow,
		  colInds, rowStarts, rowLengths, //matrixByCol,
		  coefByCol, rowInds, colStarts, colLengths,
		  cs);
}

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglMixedIntegerRounding::CglMixedIntegerRounding ()
  :
  CglCutGenerator()
{ 
  gutsOfConstruct(1, true, 1, -1);
}


//-------------------------------------------------------------------
// Alternate Constructor 
//-------------------------------------------------------------------
CglMixedIntegerRounding::CglMixedIntegerRounding (const int maxaggr,
						  const bool multiply,
						  const int criterion,
						  const int preproc)
  :
  CglCutGenerator()
{ 
  gutsOfConstruct(maxaggr, multiply, criterion, preproc);
}


//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglMixedIntegerRounding::CglMixedIntegerRounding ( 
				 const CglMixedIntegerRounding & rhs)
  :
  CglCutGenerator(rhs)
{ 
  gutsOfCopy(rhs);
}


//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglMixedIntegerRounding::clone() const
{
  return new CglMixedIntegerRounding(*this);
}

//------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglMixedIntegerRounding &
CglMixedIntegerRounding::operator=(const CglMixedIntegerRounding& rhs)
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
CglMixedIntegerRounding::~CglMixedIntegerRounding ()
{
  gutsOfDelete();
}

//-------------------------------------------------------------------
// Construct
//-------------------------------------------------------------------  
void
CglMixedIntegerRounding::gutsOfConstruct (const int maxaggr,
					  const bool multiply,
					  const int criterion,
					  const int preproc)
{
  if (maxaggr > 0) {
    MAXAGGR_ = maxaggr;
  }
  else {
    throw CoinError("Unallowable value. maxaggr must be > 0",
                      "gutsOfConstruct","CglMixedIntegerRounding");
  }
  MULTIPLY_ = multiply;
  if ((criterion >= 1) && (criterion <= 3)) {
    CRITERION_ = criterion;
  }
  else {
    throw CoinError("Unallowable value. criterion must be 1, 2 or 3",
                      "gutsOfConstruct","CglMixedIntegerRounding");
  }
  if ((preproc >= -1) && (preproc <= 2)) {
    doPreproc_ = preproc;
  }
  else {
    throw CoinError("Unallowable value. preproc must be -1, 0 or 1",
                      "gutsOfConstruct","CglMixedIntegerRounding");
  }
  EPSILON_ = 1.0e-6;
  UNDEFINED_ = -1;
  TOLERANCE_ = 1.0e-4;
  numRows_ = 0;
  numCols_ = 0;
  doneInitPre_ = false;
  vubs_ = 0;
  vlbs_ = 0;
  rowTypes_ = 0;
  indRows_ = 0;
  numRowMix_ = 0;
  indRowMix_ = 0;
  numRowCont_ = 0;
  indRowCont_ = 0;
  numRowInt_ = 0;
  indRowInt_ = 0;
  numRowContVB_ = 0;
  indRowContVB_ = 0;
  sense_=NULL;
  RHS_=NULL;
}

//-------------------------------------------------------------------
// Delete
//-------------------------------------------------------------------  
void
CglMixedIntegerRounding::gutsOfDelete ()
{
  if (vubs_ != 0) { delete [] vubs_; vubs_ = 0; }
  if (vlbs_ != 0) { delete [] vlbs_; vlbs_ = 0; }
  if (rowTypes_ != 0) { delete [] rowTypes_; rowTypes_ = 0; } 
  if (indRows_ != 0) { delete [] indRows_; indRows_ = 0; }
  if (indRowMix_ != 0) { delete [] indRowMix_; indRowMix_ = 0; }
  if (indRowCont_ != 0) { delete [] indRowCont_; indRowCont_ = 0; }
  if (indRowInt_ != 0) { delete [] indRowInt_; indRowInt_ = 0; }
  if (indRowContVB_ != 0) { delete [] indRowContVB_; indRowContVB_ = 0; }
  if (sense_ !=NULL) { delete [] sense_; sense_=NULL;}
  if (RHS_ !=NULL) { delete [] RHS_; RHS_=NULL;}
}

//-------------------------------------------------------------------
// Copy
//-------------------------------------------------------------------  
void
CglMixedIntegerRounding::gutsOfCopy (const CglMixedIntegerRounding& rhs)
{
  MAXAGGR_ = rhs.MAXAGGR_;
  MULTIPLY_ = rhs.MULTIPLY_;
  CRITERION_ = rhs.CRITERION_;
  EPSILON_ = rhs.EPSILON_;
  UNDEFINED_ = rhs.UNDEFINED_;
  TOLERANCE_ = rhs.TOLERANCE_;
  doPreproc_ = rhs.doPreproc_;
  numRows_ = rhs.numRows_;
  numCols_ = rhs.numCols_;
  doneInitPre_ = rhs.doneInitPre_;
  numRowMix_ = rhs.numRowMix_;
  numRowCont_ = rhs.numRowCont_;
  numRowInt_ = rhs.numRowInt_;
  numRowContVB_ = rhs.numRowContVB_;

  if (numCols_ > 0) {
    vubs_ = new CglMixIntRoundVUB [numCols_];
    vlbs_ = new CglMixIntRoundVLB [numCols_];
    CoinDisjointCopyN(rhs.vubs_, numCols_, vubs_);
    CoinDisjointCopyN(rhs.vlbs_, numCols_, vlbs_);
  }
  else {
    vubs_ = 0;
    vlbs_ = 0;
  }

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

  if (numRowMix_ > 0) {
    indRowMix_ = new int [numRowMix_];
    CoinDisjointCopyN(rhs.indRowMix_, numRowMix_, indRowMix_);
  }
  else {
    indRowMix_ = 0;
  }

  if (numRowCont_ > 0) {
    indRowCont_ = new int [numRowCont_];
    CoinDisjointCopyN(rhs.indRowCont_, numRowCont_, indRowCont_);
    indRowContVB_ = new int [numRowCont_];
    CoinDisjointCopyN(rhs.indRowContVB_, numRowCont_, indRowContVB_);
  }
  else {
    indRowCont_ = 0;
    indRowContVB_ = 0;
  }

  if (numRowInt_ > 0) {
    indRowInt_ = new int [numRowInt_];
    CoinDisjointCopyN(rhs.indRowInt_, numRowInt_, indRowInt_);
  }
  else {
    indRowInt_ = 0;
  }

}

//-------------------------------------------------------------------
// Do preprocessing
// It determines the type of each row. It also identifies the variable
// upper bounds and variable lower bounds.
//-------------------------------------------------------------------  
void 
CglMixedIntegerRounding::
mixIntRoundPreprocess(const OsiSolverInterface& si)
{
  // get matrix stored by row
  const CoinPackedMatrix & matrixByRow = *si.getMatrixByRow();
  numRows_ = si.getNumRows();
  numCols_ = si.getNumCols();
  const double* coefByRow  = matrixByRow.getElements();
  const int* colInds       = matrixByRow.getIndices();
  const CoinBigIndex* rowStarts     = matrixByRow.getVectorStarts();
  const int* rowLengths    = matrixByRow.getVectorLengths();
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

  // Summarize the row type information.
#if CGL_DEBUG
  int numUNDEFINED   = 0;
  int numVARUB       = 0;
  int numVARLB       = 0;
  int numVAREQ       = 0;
#endif
  int numMIX         = 0;
  int numCONT        = 0;
  int numINT         = 0;
#if CGL_DEBUG
  int numOTHER       = 0;
#endif

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
		       coefByRow+rowStarts[iRow], sense_[iRow], RHS_[iRow]);
    // store the type of the current row
    rowTypes_[iRow] = rowType;

    // Summarize information about row types
    switch(rowType) {
    case  ROW_UNDEFINED:
#if CGL_DEBUG
      ++numUNDEFINED; 
#endif
      break;
    case  ROW_VARUB:
#if CGL_DEBUG
      ++numVARUB;
#endif
      break;
    case  ROW_VARLB:
#if CGL_DEBUG
      ++numVARLB;
#endif
      break;
    case  ROW_VAREQ:
#if CGL_DEBUG
      ++numVAREQ;
#endif
      break;
    case  ROW_MIX:
      ++numMIX;
      break;
    case  ROW_CONT:
      ++numCONT;
      break;
    case  ROW_INT:
      ++numINT; 
      break;
    case  ROW_OTHER:
#if CGL_DEBUG
      ++numOTHER;
#endif
      break;
    default:
      throw CoinError("Unknown row type", "MixIntRoundPreprocess",
		      "CglMixedIntegerRounding");
    }
  }

  // allocate memory for vector of indices of all rows
  if (indRows_ != 0) { delete [] indRows_; indRows_ = 0; }
  if (numRows_ > 0)
    indRows_ = new int [numRows_];     // Destructor will free memory
  // allocate memory for vector of indices of rows of type ROW_MIX
  numRowMix_ = numMIX;
  if (indRowMix_ != 0) { delete [] indRowMix_; indRowMix_ = 0; }
  if (numRowMix_ > 0)
    indRowMix_ = new int [numRowMix_];     // Destructor will free memory
  // allocate memory for vector of indices of rows of type ROW_CONT
  numRowCont_ = numCONT;
  if (indRowCont_ != 0) { delete [] indRowCont_; indRowCont_ = 0; }
  if (numRowCont_ > 0)
    indRowCont_ = new int [numRowCont_];     // Destructor will free memory
  // allocate memory for vector of indices of rows of type ROW_INT
  numRowInt_ = numINT;
  if (indRowInt_ != 0) { delete [] indRowInt_; indRowInt_ = 0; }
  if (numRowInt_ > 0)
    indRowInt_ = new int [numRowInt_];     // Destructor will free memory

#if CGL_DEBUG
  std::cout << "The num of rows = "  << numRows_        << std::endl;
  std::cout << "Summary of Row Type" << std::endl;
  std::cout << "numUNDEFINED     = " << numUNDEFINED   << std::endl;
  std::cout << "numVARUB         = " << numVARUB       << std::endl;
  std::cout << "numVARLB         = " << numVARLB       << std::endl;
  std::cout << "numVAREQ         = " << numVAREQ       << std::endl;
  std::cout << "numMIX           = " << numMIX         << std::endl;
  std::cout << "numCONT          = " << numCONT        << std::endl;
  std::cout << "numINT           = " << numINT         << std::endl;
  std::cout << "numOTHER         = " << numOTHER       << std::endl;
#endif

  //---------------------------------------------------------------------------
  // Setup  vubs_ and vlbs_
  if (vubs_ != 0) { delete [] vubs_; vubs_ = 0; }
  vubs_ = new CglMixIntRoundVUB [numCols_]; // Destructor will free
  if (vlbs_ != 0) { delete [] vlbs_; vlbs_ = 0; }
  vlbs_ = new CglMixIntRoundVLB [numCols_]; // Destructor will free

  // Initialization. Altough this has been done in constructor, it is needed
  // for the case where the mixIntRoundPreprocess is called more than once
  for (int iCol = 0; iCol < numCols_; ++iCol) {
    vubs_[iCol].setVar(UNDEFINED_);
    vlbs_[iCol].setVar(UNDEFINED_);
  }
  const char * intVar = si.getColType();
  int countM = 0;
  int countC = 0;
  int countI = 0;
  for ( iRow = 0; iRow < numRows_; ++iRow) {

    RowType rowType = rowTypes_[iRow];

    // fill the vector indRows_ with the indices of all rows
    indRows_[iRow] = iRow;

    // fill the vector indRowMix_ with the indices of the rows of type ROW_MIX
    if (rowType == ROW_MIX) {
      indRowMix_[countM] = iRow;
      countM++;
    }
    // fill the vector indRowCont_ with the indices of rows of type ROW_CONT
    else if (rowType == ROW_CONT) {
      indRowCont_[countC] = iRow;
      countC++;
    }
    // fill the vector indRowInt_ with the indices of the rows of type ROW_INT
    else if (rowType == ROW_INT) {
      indRowInt_[countI] = iRow;
      countI++;
    }
    // create vectors with variable lower and upper bounds
    else if ( (rowType == ROW_VARUB) || 
	      (rowType == ROW_VARLB) || 
	      (rowType == ROW_VAREQ) )  { 
      
      CoinBigIndex startPos = rowStarts[iRow];
      CoinBigIndex stopPos  = startPos + rowLengths[iRow];
      int    xInd = 0,  yInd = 0;   // x is continuous, y is integer
      double xCoef = 0.0, yCoef = 0.0;

      for (CoinBigIndex i = startPos; i < stopPos; ++i) {
	if ( fabs(coefByRow[i]) > EPSILON_ ) {
	  if( intVar[colInds[i]] ) {
	    yInd  = colInds[i];
	    yCoef = coefByRow[i];
	  }
	  else {
	    xInd  = colInds[i];
	    xCoef = coefByRow[i];
	  }
	}
      }

      switch (rowType) {
      case ROW_VARUB:       // Inequality: x <= ? * y
	vubs_[xInd].setVar(yInd);
	vubs_[xInd].setVal(-yCoef / xCoef);
	break;
      case ROW_VARLB:       // Inequality: x >= ? * y
	vlbs_[xInd].setVar(yInd);
	vlbs_[xInd].setVal(-yCoef / xCoef);
	break;
      case ROW_VAREQ:       // Inequality: x >= AND <= ? * y
	vubs_[xInd].setVar(yInd);
	vubs_[xInd].setVal(-yCoef / xCoef);
	vlbs_[xInd].setVar(yInd);
	vlbs_[xInd].setVal(-yCoef / xCoef);
	break;
      default:
        // I am getting compiler bug which gets here - I am disabling - JJF
	//throw CoinError("Unknown row type: impossible", 
        //	"MixIntRoundPreprocess",
        //	"CglMixedIntegerRounding");
        break;
      }
    }
  }

  // allocate memory for vector of indices of rows of type ROW_CONT
  // that have at least one variable with variable upper or lower bound
  if (indRowContVB_ != 0) { delete [] indRowContVB_; indRowContVB_ = 0; }
  if (numRowCont_ > 0)
    indRowContVB_ = new int [numRowCont_];     // Destructor will free memory
  // create vector with rows of type ROW_CONT that have at least
  // one variable with variable upper or lower bound
  countC = 0;
  for (int i = 0; i < numRowCont_; ++i) {
    int indRow = indRowCont_[i];
    CoinBigIndex jStart = rowStarts[indRow];
    CoinBigIndex jStop = jStart + rowLengths[indRow];
    for (CoinBigIndex j = jStart; j < jStop; ++j) {
      int indCol = colInds[j];
      CglMixIntRoundVLB VLB = vlbs_[indCol];
      CglMixIntRoundVUB VUB = vubs_[indCol];
      if (( VLB.getVar() != UNDEFINED_ ) || ( VUB.getVar() != UNDEFINED_ ) ){
	indRowContVB_[countC] = indRow;
	countC++;
	break;
      }
    }
  }
  numRowContVB_ = countC;

}

//-------------------------------------------------------------------
// Determine the type of a given row 
//-------------------------------------------------------------------
CglMixedIntegerRounding::RowType
CglMixedIntegerRounding::determineRowType(const OsiSolverInterface& si,
				  const int rowLen, const int* ind, 
				  const double* coef, const char sense, 
				  const double rhs) const
{
  if (rowLen == 0) 
    return ROW_UNDEFINED;

  if (sense == 'N' || rhs == si.getInfinity() || rhs == -si.getInfinity())
    return ROW_OTHER;

  RowType rowType = ROW_UNDEFINED;

  int  numPosInt = 0;      // num of positive integer variables
  int  numNegInt = 0;      // num of negative integer variables
  int  numInt    = 0;      // num of integer variables
  int  numPosCon = 0;      // num of positive continuous variables
  int  numNegCon = 0;      // num of negative continuous variables
  int  numCon    = 0;      // num of continuous variables

  const char * intVar = si.getColType();
  // Summarize the variable types of the given row.
  for ( int i = 0; i < rowLen; ++i ) {
    if ( coef[i] < -EPSILON_ ) {
      if( intVar[ind[i]] )
	++numNegInt;
      else
	++numNegCon;
    }
    else if ( coef[i] > EPSILON_ ) {
      if( intVar[ind[i]] )
	++numPosInt;
      else
	++numPosCon;
    }
  }
  numInt = numNegInt + numPosInt;
  numCon = numNegCon + numPosCon;

#if CGL_DEBUG
  std::cout << "numNegInt = " << numNegInt << std::endl;
  std::cout << "numPosInt = " << numPosInt << std::endl;
  std::cout << "numInt = " << numInt << std::endl;
  std::cout << "numNegCon = " << numNegCon << std::endl;
  std::cout << "numPosCon = " << numPosCon << std::endl;
  std::cout << "numCon = " << numCon << std::endl;
  std::cout << "rowLen = " << rowLen << std::endl;
#endif


  //-------------------------------------------------------------------------
  // Classify row type based on the types of variables.
    
  if ((numInt > 0) && (numCon > 0)) {
    if ((numInt == 1) && (numCon == 1) && (fabs(rhs) <= EPSILON_)) {
      // It's a variable bound constraint
      switch (sense) {
      case 'L':
	rowType = numPosCon == 1 ? ROW_VARUB : ROW_VARLB;
	break;
      case 'G':
	rowType = numPosCon == 1 ? ROW_VARLB : ROW_VARUB;
	break;
      case 'E':
        rowType = ROW_VAREQ;
	break;
      default:
	break;
      }
    }
    else {
      // It's a constraint with continuous and integer variables;
      // The total number of variables is at least 2
      rowType = ROW_MIX;
    }
  }
  else if (numInt == 0) {
    // It's a constraint with only continuous variables
    rowType = ROW_CONT;
  }
  else if ((numCon == 0) && ((sense == 'L') || (sense == 'G'))) {
    // It's a <= or >= constraint with only integer variables 
    rowType = ROW_INT;
  }
  else
    // It's a constraint that does not fit the above categories
    rowType = ROW_OTHER;


  return rowType;
}

//-------------------------------------------------------------------
// Generate MIR cuts
//-------------------------------------------------------------------
void
CglMixedIntegerRounding::generateMirCuts( 
			    const OsiSolverInterface& si,
			    const double* xlp,
			    const double* colUpperBound,
			    const double* colLowerBound,
			    const CoinPackedMatrix& matrixByRow,
			    const double* LHS,
			    const double* /*coefByRow*/,
			    const int* /*colInds*/,
			    const CoinBigIndex* /*rowStarts*/,
			    const int* /*rowLengths*/,
			    //const CoinPackedMatrix& matrixByCol,
			    const double* coefByCol,
			    const int* rowInds,
			    const CoinBigIndex* colStarts,
			    const int* colLengths,
			    OsiCuts& cs ) const
{

#if CGL_DEBUG
  // Open debug data file; incorporate solver name so we get separate files
  // when running unit test.
  std::string dbgFname ;
  si.getStrParam(OsiSolverName,dbgFname) ;
  dbgFname = "mir_"+dbgFname+"_stats.dat" ;
  std::ofstream fout(dbgFname.c_str()) ;
#endif

  // Define upper limit for the loop where the cMIRs are constructed
  int upperLimit;
  if (MULTIPLY_)
    upperLimit = 2;
  else
    upperLimit = 1;
  
  // create a vector with the columns that were used in the aggregation
  int* listColsSelected = new int[MAXAGGR_];
  // create a vector with the rows that were aggregated
  int* listRowsAggregated = new int[MAXAGGR_];
  // create a vector with the LP solutions of the slack variables
  double* xlpExtra = new double[MAXAGGR_];

  // loop until maximum number of aggregated rows is reached or a 
  // violated cut is found
  int numRowMixAndRowContVB = numRowMix_ + numRowContVB_;
  int numRowMixAndRowContVBAndRowInt = numRowMixAndRowContVB + numRowInt_;
  for (int iRow = 0; iRow < numRowMixAndRowContVBAndRowInt; ++iRow) {

    int rowSelected;  // row selected to be aggregated next
    int colSelected;  // column selected for pivot in aggregation
    CoinPackedVector rowAggregated;
    double rhsAggregated;
    // create a set with the indices of rows selected
    std::set<int> setRowsAggregated;

    // loop until the maximum number of aggregated rows is reached
    for (int iAggregate = 0; iAggregate < MAXAGGR_; ++iAggregate) {

      if (iAggregate == 0) {

	// select row
	if (iRow < numRowMix_) {
	  rowSelected = indRowMix_[iRow];
	}
	else if (iRow < numRowMixAndRowContVB) {
	  rowSelected = indRowContVB_[iRow - numRowMix_];
	}
	else {
	  rowSelected = indRowInt_[iRow - numRowMixAndRowContVB];
	}

	copyRowSelected(iAggregate, rowSelected, setRowsAggregated,
			listRowsAggregated, xlpExtra, sense_[rowSelected], 
			RHS_[rowSelected], LHS[rowSelected], 
			matrixByRow, rowAggregated, rhsAggregated);

      } 
      else {

	// search for a row to aggregate
	bool foundRowToAggregate = selectRowToAggregate(
				        si, rowAggregated,
					colUpperBound, colLowerBound, 
					setRowsAggregated, xlp, 
					coefByCol, rowInds, colStarts,
					colLengths, 
					rowSelected, colSelected);

	// if finds row to aggregate, compute aggregated row
	if (foundRowToAggregate) {

	  CoinPackedVector rowToAggregate;
	  double rhsToAggregate;

	  listColsSelected[iAggregate] = colSelected;

	  copyRowSelected(iAggregate, rowSelected, setRowsAggregated,
			  listRowsAggregated, xlpExtra, sense_[rowSelected], 
			  RHS_[rowSelected], LHS[rowSelected], 
			  matrixByRow, rowToAggregate, rhsToAggregate);

	  // call aggregate row heuristic
	  aggregateRow(colSelected, rowToAggregate, rhsToAggregate, 
		       rowAggregated, rhsAggregated);

	}
	else
	  break;
      }


      // construct cMIR with current rowAggregated
      // and, if upperLimit=2 construct also a cMIR with 
      // the current rowAggregated multiplied by -1
      for (int i = 0; i < upperLimit; ++i) {
      
	// create vector for mixed knapsack constraint
	CoinPackedVector rowToUse = rowAggregated;
	double rhsMixedKnapsack = rhsAggregated;
	if (i == 1) {
	  rowToUse *= (-1.0);
	  rhsMixedKnapsack *= (-1.0);
	}	  
	CoinPackedVector mixedKnapsack;
	double sStar = 0.0;

	// create vector for the continuous variables in s
	CoinPackedVector contVariablesInS;

	// call bound substitution heuristic
	bool foundMixedKnapsack = boundSubstitution(
					si, rowToUse, 
					xlp, xlpExtra, 
					colUpperBound, colLowerBound,
					mixedKnapsack, rhsMixedKnapsack, 
					sStar, contVariablesInS);
        // may want some limit?
        if (mixedKnapsack.getNumElements()>25000) {
#if CGL_DEBUG	  
	  std::cout << "mixed knapsack has " 
                    <<mixedKnapsack.getNumElements()<<" elements - rhs is "
                    <<rhsMixedKnapsack
                    << std::endl;
#endif
	  continue;
	}
          
	// if it did not find a mixed knapsack it is because there is at
	// least one integer variable with lower bound different than zero
	// or there are no integer or continuous variables.
	// In this case, we continue without trying to generate a c-MIR
	if (!foundMixedKnapsack) {
#if CGL_DEBUG	  
	  std::cout << "couldn't create mixed knapsack" << std::endl;
#endif
	  continue;
	}

	OsiRowCut cMirCut;

	// Find a c-MIR cut with the current mixed knapsack constraint
	bool hasCut = cMirSeparation(si, matrixByRow, rowToUse,
				     listRowsAggregated, sense_, RHS_,
				     //coefByRow, colInds, rowStarts, rowLengths,
				     xlp, sStar, colUpperBound, colLowerBound, 
				     mixedKnapsack,
				     rhsMixedKnapsack, contVariablesInS,
				     cMirCut);

#if CGL_DEBUG
	// PRINT STATISTICS
	printStats(fout, hasCut, si, rowAggregated, rhsAggregated, xlp,
		   xlpExtra, listRowsAggregated, listColsSelected, 
		   iAggregate+1, colUpperBound, colLowerBound );
#endif

	// if a cut was found, insert it into cs
	if (hasCut)  {
#if CGL_DEBUG
	  std::cout << "MIR cut generated " << std::endl;
#endif
	  cs.insertIfNotDuplicate(cMirCut);
	}

      }
	
    }

  }

  // free memory
  delete [] listColsSelected; listColsSelected = 0;
  delete [] listRowsAggregated; listRowsAggregated = 0;
  delete [] xlpExtra; xlpExtra = 0;
  
#if CGL_DEBUG
  // CLOSE FILE
  fout.close();
#endif

  return;

}

//-------------------------------------------------------------------
// Copy row selected to CoinPackedVector
//-------------------------------------------------------------------
void
CglMixedIntegerRounding::copyRowSelected(
			    const int iAggregate,
			    const int rowSelected,
			    std::set<int>& setRowsAggregated,
			    int* listRowsAggregated,
			    double* xlpExtra,
			    const char sen,
			    const double rhs,
			    const double lhs,
			    const CoinPackedMatrix& matrixByRow,
			    CoinPackedVector& rowToAggregate,
			    double& rhsToAggregate) const
{

  // copy the row selected to a vector of type CoinPackedVector
  const CoinShallowPackedVector reqdBySunCC = matrixByRow.getVector(rowSelected);
  rowToAggregate = reqdBySunCC ;
  rhsToAggregate = rhs;

  // update list of indices of rows selected
  setRowsAggregated.insert(rowSelected);
  listRowsAggregated[iAggregate] = rowSelected;

  // Add a slack variable if needed and compute its current value
  if (sen == 'L') {
    rowToAggregate.insert(numCols_ + iAggregate, 1);
    xlpExtra[iAggregate] = rhs - lhs;
  }
  else if (sen == 'G') {
    rowToAggregate.insert(numCols_ + iAggregate, -1);
    xlpExtra[iAggregate] = lhs - rhs;
  }

}

//-------------------------------------------------------------------
// Construct the set P* and select a row to aggregate
//-------------------------------------------------------------------
bool
CglMixedIntegerRounding::selectRowToAggregate( 
			    const OsiSolverInterface& si,
			    const CoinPackedVector& rowAggregated,
			    const double* colUpperBound,
			    const double* colLowerBound,
			    const std::set<int>& setRowsAggregated,
			    const double* xlp, const double* coefByCol,
			    const int* rowInds, const CoinBigIndex* colStarts,
			    const int* colLengths,
			    int& rowSelected,
			    int& colSelected ) const
{

  bool foundRowToAggregate = false;

  double deltaMax = 0.0;  // maximum delta
  const int numColsAggregated = rowAggregated.getNumElements();
  const int *rowAggregatedIndices = rowAggregated.getIndices();
  const double *rowAggregatedElements = rowAggregated.getElements();  

  for (int j = 0; j < numColsAggregated; ++j) {

    // store the index and coefficient of column j
    int indCol = rowAggregatedIndices[j];
    if (indCol >= numCols_) continue;
    double coefCol = rowAggregatedElements[j];

    // Consider only continuous variables
    if ( (!si.isContinuous(indCol)) || (fabs(coefCol) < EPSILON_)) continue;

    // Compute current lower bound
    CglMixIntRoundVLB VLB = vlbs_[indCol];
    double LB = ( VLB.getVar() != UNDEFINED_ ) ? 
                      VLB.getVal() * xlp[VLB.getVar()] : colLowerBound[indCol];
    
    // Compute current upper bound
    CglMixIntRoundVUB VUB = vubs_[indCol];
    double UB = ( VUB.getVar() != UNDEFINED_ ) ? 
                      VUB.getVal() * xlp[VUB.getVar()] : colUpperBound[indCol];

    // Compute distances from current solution to upper and lower bounds
    double delta = std::min(xlp[indCol] - LB, UB - xlp[indCol]);

    // In case this variable is acceptable look for possible rows
    if (delta > deltaMax) {

      CoinBigIndex iStart = colStarts[indCol];
      CoinBigIndex iStop  = iStart + colLengths[indCol];
      //      int count = 0;

      //      std::vector<int> rowPossible;

      // find a row to use in aggregation
      for (CoinBigIndex i = iStart; i < iStop; ++i) {
	int rowInd = rowInds[i];
	if (setRowsAggregated.find(rowInd) == setRowsAggregated.end()) {
	  // if the row was not already selected, select it
	  RowType rType = rowTypes_[rowInd];
	  if ( ((rType == ROW_MIX) || (rType == ROW_CONT)) 
	       && (fabs(coefByCol[i]) > EPSILON_) ) {
	    //	    rowPossible.push_back(rowInd);
	    rowSelected = rowInd;
	    deltaMax = delta;
	    colSelected = indCol;
	    foundRowToAggregate = true;
	    //count++;
	    break;
	  }
	}
      }

      //      if (count > 0)
      //	rowSelected = rowPossible[rand() % count];
      //      std::cout << count << std::endl;
    }
	
  }

  return foundRowToAggregate;

}
      
//-------------------------------------------------------------------
// Aggregate the selected row with the current aggregated row
//-------------------------------------------------------------------
void
CglMixedIntegerRounding::aggregateRow( 
			    const int colSelected,
			    CoinPackedVector& rowToAggregate, double rhs,
			    CoinPackedVector& rowAggregated, 
			    double& rhsAggregated ) const
{

  // quantity to multiply by the coefficients of the row to aggregate
  double multiCoef = rowAggregated[colSelected] / rowToAggregate[colSelected];

  rowToAggregate *= multiCoef; 
  rhs *= multiCoef;

  rowAggregated = rowAggregated - rowToAggregate;
  rhsAggregated -= rhs;

}

//-------------------------------------------------------------------
// Choose the bound substitution based on the criteria defined by the user
//-------------------------------------------------------------------
inline bool
CglMixedIntegerRounding::isLowerSubst(const double inf, 
				      const double aj,
				      const double xlp, 
				      const double LB, 
				      const double UB) const
{
  if (CRITERION_ == 1) {
    // criterion 1 (the same as criterion (a) in the paper)
    return xlp - LB < UB - xlp;
  }
  else {
    if (UB == inf || xlp == LB) 
      return true;
    if (LB == -inf || xlp == UB)
      return false;
    if (CRITERION_ == 2) 
      // criterion 2 (the same as criterion (b) in the paper)
      return aj < 0;
    else
      // criterion 3 (the same as criterion (c) in the paper)
      return aj > 0;
  }
}



//-------------------------------------------------------------------
// Bound substitution heuristic
//-------------------------------------------------------------------
bool
CglMixedIntegerRounding::boundSubstitution( 
			    const OsiSolverInterface& si,
			    const CoinPackedVector& rowAggregated,
			    const double* xlp,
			    const double* xlpExtra,
			    const double* colUpperBound,
			    const double* colLowerBound,
			    CoinPackedVector& mixedKnapsack,
			    double& rhsMixedKnapsack, double& sStar,
			    CoinPackedVector& contVariablesInS ) const
{

  bool generated = false;
  const int numColsAggregated = rowAggregated.getNumElements();
  const int *rowAggregatedIndices = rowAggregated.getIndices();
  const double *rowAggregatedElements = rowAggregated.getElements();  

  // go through all the variables and if it is continuous and delta is 
  // negative, store variable in the vector contVariablesInS.
  // If it is integer, store variable in the vector mixedKnapsack
  int numCont = 0;
  int j;
  for ( j = 0; j < numColsAggregated; ++j) {

    // get index and coefficient of column j in the aggregated row
    const int indCol = rowAggregatedIndices[j];
    const double coefCol = rowAggregatedElements[j];

    // if the lower bound is equal to the upper bound, remove variable
    if ( (indCol < numCols_) &&
	 (colLowerBound[indCol] == colUpperBound[indCol]) ) {
      rhsMixedKnapsack -= coefCol * colLowerBound[indCol];
      continue;
    }

    if (fabs(coefCol) < EPSILON_) continue;
    // set the coefficients of the integer variables
    if ( (indCol < numCols_)  && (!si.isContinuous(indCol)) ) {
      // Copy the integer variable to the vector mixedKnapsack
      if (mixedKnapsack.isExistingIndex(indCol)) {
	const int index = mixedKnapsack.findIndex(indCol);
	mixedKnapsack.setElement(index, mixedKnapsack[indCol] + coefCol);
      }
      else
	mixedKnapsack.insert(indCol, coefCol);
      continue;
    }

    // Select the continuous variables and copy the ones in s to 
    // the vector contVariablesInS
    if (indCol < numCols_) {  // variable is model variable

      // Compute lower bound for variable indCol
      const CglMixIntRoundVLB VLB = vlbs_[indCol];
      const double LB = ( VLB.getVar() != UNDEFINED_ ) ? 
	        VLB.getVal() * xlp[VLB.getVar()] : colLowerBound[indCol];
    
      // Compute upper bound for variable indCol
      const CglMixIntRoundVUB VUB = vubs_[indCol];
      const double UB = ( VUB.getVar() != UNDEFINED_ ) ? 
	        VUB.getVal() * xlp[VUB.getVar()] : colUpperBound[indCol];

      // if both bounds are infinite, then we cannot form a mixed knapsack
      if ( (LB == -1.0 * si.getInfinity()) &&
	   (UB == si.getInfinity()) ) {
#if CGL_DEBUG
	std::cout << "continuous var with infinite bounds. " <<
                     "Cannot form mixed Knapsack = " << std::endl;
#endif
	return generated;
      }

      // Select the bound substitution
      if (isLowerSubst(si.getInfinity(), rowAggregatedElements[j],
			 xlp[indCol], LB, UB)) {
	if (VLB.getVar() != UNDEFINED_ ) {
	  const int indVLB = VLB.getVar();
	  if (mixedKnapsack.isExistingIndex(indVLB)) {
	    const int index = mixedKnapsack.findIndex(indVLB);
	    mixedKnapsack.setElement(index, mixedKnapsack[indVLB] + 
				     coefCol * VLB.getVal());
	  }
	  else
	    mixedKnapsack.insert(indVLB, coefCol * VLB.getVal());
	}
	else {
	  rhsMixedKnapsack -= coefCol * LB;
	}
	// Update sStar
	if (coefCol < -EPSILON_) {
	  contVariablesInS.insert(indCol, coefCol);
	  sStar -= coefCol * (xlp[indCol] - LB);
	  numCont++;
	}
      }
      else {
	if (VUB.getVar() != UNDEFINED_ ) {
	  const int indVUB = VUB.getVar();
	  if (mixedKnapsack.isExistingIndex(indVUB)) {
	    const int index = mixedKnapsack.findIndex(indVUB);
	    mixedKnapsack.setElement(index, mixedKnapsack[indVUB] + 
				     coefCol * VUB.getVal());
	  }
	  else
	    mixedKnapsack.insert(indVUB, coefCol * VUB.getVal());
	}
	else {
	  rhsMixedKnapsack -= coefCol * UB;
	}
	// Update sStar
	if (coefCol > EPSILON_) {
	  contVariablesInS.insert(indCol, - coefCol);
	  sStar += coefCol * (UB - xlp[indCol]);
	  numCont++;
	}
      }
    }
    else {  // variable is slack variable
      // in this case the LB = 0 and the UB = infinity
      // Update sStar
      const double tLB = xlpExtra[indCol - numCols_];
      if (coefCol < -EPSILON_) {
	contVariablesInS.insert(indCol, coefCol);
	sStar -= coefCol * tLB;
	numCont++;
      }
    }

  }

  // if there are no continuous variables to form s, then we stop
#if CGL_DEBUG
  std::cout << "# of continuous var in mixedKnapsack = " << numCont <<
    std::endl;
#endif
  if (numCont == 0) return generated;

  // check that the integer variables have lower bound equal to zero
  const int numInt = mixedKnapsack.getNumElements();
  // if there are not integer variables in mixedKnapsack, then we stop
  // CAUTION: all the coefficients could be zero
#if CGL_DEBUG
  std::cout << "# of integer var in mixedKnapsack = " << numInt <<
    std::endl;
#endif
  if (numInt == 0) return generated;
  const int *knapsackIndices = mixedKnapsack.getIndices();
  const double *knapsackElements = mixedKnapsack.getElements();  

  for ( j = 0; j < numInt; ++j) {
    // if the coefficient is zero, disregard
    if (fabs(knapsackElements[j]) < EPSILON_) continue;
    // if the lower bound is not zero, then we stop
    if (fabs(colLowerBound[knapsackIndices[j]]) > EPSILON_) return generated;
  }
  // if the lower bounds of all integer variables are zero, proceed
  generated = true;
  return generated;

}

//-------------------------------------------------------------------
// c-MIR separation heuristic
//-------------------------------------------------------------------
bool
CglMixedIntegerRounding::cMirSeparation( 
			    const OsiSolverInterface& si,
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
			    OsiRowCut& cMirCut) const
{

  bool generated = false;
  double numeratorBeta = rhsMixedKnapsack;
  CoinPackedVector cMIR = mixedKnapsack;
  double rhscMIR;
  double maxViolation = 0.0;
  double bestDelta = 0.0;
  CoinPackedVector bestCut;
  double rhsBestCut = 0.0;
  double sCoefBestCut = 0.0;
  const int numInt = mixedKnapsack.getNumElements();  
  const int *knapsackIndices = mixedKnapsack.getIndices();
  const double *knapsackElements = mixedKnapsack.getElements();  
  const int *contVarInSIndices = contVariablesInS.getIndices();
  const double *contVarInSElements = contVariablesInS.getElements();

  // Construct set C, T will be the rest.
  // Also, for T we construct a CoinPackedVector named complT which
  // contains the vars in T that are strictly between their bounds
  std::set<int> setC;
  CoinPackedVector complT;
  int j;
  for ( j = 0; j < numInt; ++j) {
    const int indCol = knapsackIndices[j];
    // if the upper bound is infinity, then indCol is in T and cannot
    // be in complT
    if (colUpperBound[indCol] != si.getInfinity()) {
      if (xlp[indCol] >= colUpperBound[indCol] / 2.0) {
	setC.insert(j);
	numeratorBeta -= knapsackElements[j] * colUpperBound[indCol];
      } else {
	if ( (xlp[indCol] <= EPSILON_) || 
	     (xlp[indCol] >= colUpperBound[indCol] - EPSILON_))
	  continue;
	complT.insert(j, fabs(xlp[indCol] - colUpperBound[indCol]/2));
      }
    }
  }

  // Sort the indices in complT by nondecreasing values 
  // (which are  $|y^*_j-u_j/2|$)
  if (complT.getNumElements() > 0) {
    complT.sortIncrElement();
  }

  // Construct c-MIR inequalities and take the one with the largest violation
  for ( j = 0; j < numInt; ++j) {
    int indCol = knapsackIndices[j];
    if ( (xlp[indCol] <= EPSILON_) || 
	 (xlp[indCol] >= colUpperBound[indCol] - EPSILON_))
      continue;
    double delta = knapsackElements[j];
    // delta has to be positive
    if (delta <= EPSILON_) continue;

    double violation = 0.0;
    double sCoef = 0.0;

    // form a cMIR inequality
    cMirInequality(numInt, delta, numeratorBeta, knapsackIndices, 
		   knapsackElements, xlp, sStar, colUpperBound, setC, cMIR,
		   rhscMIR, sCoef, violation);

    // store cut if it is the best found so far
    if (violation > maxViolation + EPSILON_) {
      bestCut = cMIR;
      rhsBestCut = rhscMIR;
      sCoefBestCut = sCoef;
      maxViolation = violation;
      bestDelta = delta;
    }
  }

  // if no violated inequality has been found, exit now
  if (maxViolation == 0.0) return generated;

  // improve the best violated inequality.
  // try to divide delta by 2, 4 or 8 and see if increases the violation
  double deltaBase = bestDelta;
  for (int multFactor = 2; multFactor <= 8; multFactor *= 2) {
    double delta = deltaBase / multFactor;
    double violation = 0.0;
    double sCoef = 0.0;

    // form a cMIR inequality
    cMirInequality(numInt, delta, numeratorBeta, knapsackIndices, 
		   knapsackElements, xlp, sStar, colUpperBound, setC, cMIR,
		   rhscMIR, sCoef, violation);

    // store cut if it is the best found so far
    if (violation > maxViolation + EPSILON_) {
      bestCut = cMIR;
      rhsBestCut = rhscMIR;
      sCoefBestCut = sCoef;
      maxViolation = violation;
      bestDelta = delta;
    }
  }

  // improve cMIR for the best delta
  // complT contains indices into mixedKnapsack for the variables
  // which may be complemented and they are already appropriately
  // sorted.
  const int complTSize = complT.getNumElements();
  if (complTSize > 0) {
    const int *complTIndices = complT.getIndices();
    for (int j = 0; j < complTSize; ++j) {
      // move variable in set complT from set T to set C
      int jIndex = complTIndices[j];
      int indCol = knapsackIndices[jIndex];
      // do nothing if upper bound is infinity
      if (colUpperBound[indCol] >= si.getInfinity()) continue;
      setC.insert(jIndex);
      double violation = 0.0;
      double sCoef = 0.0;
      double localNumeratorBeta = numeratorBeta -
	mixedKnapsack[indCol] * colUpperBound[indCol];

      // form a cMIR inequality
      cMirInequality(numInt, bestDelta, localNumeratorBeta, knapsackIndices, 
		     knapsackElements, xlp, sStar, colUpperBound, setC, cMIR,
		     rhscMIR, sCoef, violation);

      // store cut if it is the best found so far; otherwise, move the variable
      // that was added to set C back to set T
      if (violation > maxViolation + EPSILON_) {
	bestCut = cMIR;
	rhsBestCut = rhscMIR;
	sCoefBestCut = sCoef;
	maxViolation = violation;
	numeratorBeta = localNumeratorBeta;
      }
      else
	setC.erase(jIndex);
    }
  }

  // write the best cut found with the model variables
  int numCont = contVariablesInS.getNumElements();
  for ( j = 0; j < numCont; ++j) {
    int indCol = contVarInSIndices[j];
    double coefCol = contVarInSElements[j];
      
    if (indCol < numCols_) {  // variable is model variable

      // Compute lower bound for variable indCol
      CglMixIntRoundVLB VLB = vlbs_[indCol];
      double LB = ( VLB.getVar() != UNDEFINED_ ) ? 
	VLB.getVal() * xlp[VLB.getVar()] : colLowerBound[indCol];
    
      // Compute upper bound for variable indCol
      CglMixIntRoundVUB VUB = vubs_[indCol];
      double UB = ( VUB.getVar() != UNDEFINED_ ) ? 
	VUB.getVal() * xlp[VUB.getVar()] : colUpperBound[indCol];

      // Select the bound substitution
      if (isLowerSubst(si.getInfinity(), rowAggregated[indCol],
			 xlp[indCol], LB, UB)) { 
	if (VLB.getVar() != UNDEFINED_ ) {
	  int indVLB = VLB.getVar();
	  if (bestCut.isExistingIndex(indVLB)){
	    int index = bestCut.findIndex(indVLB);
	    bestCut.setElement(index, bestCut[indVLB] - 
			       sCoefBestCut * coefCol * VLB.getVal());
	  }
	  else
	    bestCut.insert(indVLB, - sCoefBestCut * coefCol * VLB.getVal());
	  bestCut.insert(indCol, sCoefBestCut * coefCol);
	}
	else {
	  rhsBestCut += sCoefBestCut * coefCol * colLowerBound[indCol];
	  bestCut.insert(indCol, sCoefBestCut * coefCol);
	}
      }
      else {
	if (VUB.getVar() != UNDEFINED_ ) {
	  int indVUB = VUB.getVar();
	  if (bestCut.isExistingIndex(indVUB)){
	    int index = bestCut.findIndex(indVUB);
	    bestCut.setElement(index, bestCut[indVUB] + 
			       sCoefBestCut * coefCol * VUB.getVal());
	  }
	  else
	    bestCut.insert(indVUB, sCoefBestCut * coefCol * VUB.getVal());
	  bestCut.insert(indCol, - sCoefBestCut * coefCol);
	}
	else {
	  rhsBestCut -= sCoefBestCut * coefCol * colUpperBound[indCol];
	  bestCut.insert(indCol, - sCoefBestCut * coefCol);
	}
      }
    }
    else {  // variable is slack variable
      // in this case the LB = 0 and the UB = infinity
      // copy the row selected to a vector of type CoinPackedVector
      const int iRow = listRowsAggregated[indCol - numCols_];
      const CoinShallowPackedVector reqdBySunCC = matrixByRow.getVector(iRow);
      CoinPackedVector row = reqdBySunCC ;
      double rhs     = RHS[iRow];

      if (sense[iRow] == 'L') {
	// if it is a <= inequality, the coefficient of the slack is 1
	row *= (- sCoefBestCut * coefCol);
	rhs *= (- sCoefBestCut * coefCol);
      }
      else {
        assert (sense[iRow]=='G');
	// if it is a <= inequality, the coefficient of the slack is -1
	row *= (sCoefBestCut * coefCol);
	rhs *= (sCoefBestCut * coefCol);
      }

      rhsBestCut += rhs;
      bestCut = bestCut + row;
    }
  }

  // Check the violation of the cut after it is written with the original
  // variables.
  int cutLen = bestCut.getNumElements();
  int* cutInd = bestCut.getIndices();
  double* cutCoef = bestCut.getElements();
  double cutRHS = rhsBestCut;
  double violation = 0.0;
  double normCut = 0.0;
  // Also weaken by small coefficients
  int n=0;
  for ( j = 0; j < cutLen; ++j) {
    double value = cutCoef[j];
    int column = cutInd[j];
    if (fabs(value)>1.0e-12) {
      violation += cutCoef[j] * xlp[column];
      normCut += cutCoef[j] * cutCoef[j];
      cutCoef[n]=value;
      cutInd[n++]=column;
    } else if (value) {
      // Weaken
      if (value>0.0) {
        // Allow for at lower bound
        cutRHS -= value*colLowerBound[column];
      } else {
        // Allow for at upper bound
        cutRHS -= value*colUpperBound[column];
      }
    }
  }
  cutLen=n;
  violation -= cutRHS;
  violation /= sqrt(normCut);

  if ( violation > TOLERANCE_ ) {
    cMirCut.setRow(cutLen, cutInd, cutCoef);
    cMirCut.setLb(-1.0 * si.getInfinity());
    cMirCut.setUb(cutRHS);
    cMirCut.setEffectiveness(violation);
#ifdef CGL_DEBUG
    {
      for (int k=0; k<cutLen; k++){
	assert(cutInd[k]>=0);
	assert(cutCoef[k]);
        assert (fabs(cutCoef[k])>1.0e-12);
      }
    }
#endif
    generated = true;
  }

  return generated;

}

//-------------------------------------------------------------------
// construct a c-MIR inequality
//-------------------------------------------------------------------
void
CglMixedIntegerRounding::cMirInequality(
				  const int numInt, 
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
				  double& violation) const
{

      // form a cMIR inequality
      double beta = numeratorBeta / delta;
      double f = beta - floor(beta);
      rhscMIR = floor(beta);
      double normCut = 0.0;
      // coefficients of variables in set T
      for (int i = 0; i < numInt; ++i) {
	const int iIndex = knapsackIndices[i];
	double G = 0.0;
	if (setC.find(i) == setC.end()) {
	  // i is not in setC, i.e., it is in T
	  G = functionG(knapsackElements[i] / delta, f);
	  violation += (G * xlp[iIndex]);
	  normCut += G * G;
	  cMIR.setElement(i, G);
	} else {
	  G = functionG( - knapsackElements[i] / delta, f);
	  violation -= (G * xlp[iIndex]);
	  normCut += G * G;
	  rhscMIR -= G * colUpperBound[iIndex];
	  cMIR.setElement(i, -G);	  
	}
      }
      sCoef = 1.0 / (delta * (1.0 - f));
      violation -= (rhscMIR + sCoef * sStar);
      normCut += sCoef * sCoef;
      violation /= sqrt(normCut);

}


//-------------------------------------------------------------------
// function G for computing coefficients in cMIR inequality
//-------------------------------------------------------------------
inline double
CglMixedIntegerRounding::functionG( const double d, const double f ) const
{
  double delta = d - floor(d) - f;
  if (delta > EPSILON_)
    return floor(d) + delta / (1 - f);
  else
    return floor(d);
}

//-------------------------------------------------------------------
// Printing statistics
//-------------------------------------------------------------------
void
CglMixedIntegerRounding::printStats(
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
			    const double* colLowerBound ) const
{


  const int numColsAggregated = rowAggregated.getNumElements();
  const int *rowAggregatedIndices = rowAggregated.getIndices();
  const double *rowAggregatedElements = rowAggregated.getElements();  

  fout << "Rows ";
  for (int i = 0; i < level; ++i) {
    fout << listRowsAggregated[i] << " ";
  }
  fout << std::endl;

  int numColsBack = 0;

  // go through all the variables 
  for (int j = 0; j < numColsAggregated; ++j) {

    // get index and coefficient of column j in the aggregated row
    int indCol = rowAggregatedIndices[j];
    double coefCol = rowAggregatedElements[j];

    // check if a column used in aggregation is back into the aggregated row
    for (int i = 0; i < level-1; ++i) {
      if ( (listColsSelected[i] == indCol) && (coefCol != 0) ) {
	numColsBack++;
	break;
      }
    }



    if (fabs(coefCol) < EPSILON_) {
      // print variable number and coefficient
      fout << indCol << " " << 0.0 << std::endl;
      continue;
    }
    else {
      // print variable number and coefficient
      fout << indCol << " " << coefCol << " ";
    }

    // integer variables
    if ( (indCol < numCols_)  && (!si.isContinuous(indCol)) ) {

      // print 
      fout << "I " << xlp[indCol] << " " << colLowerBound[indCol] <<
	" " << colUpperBound[indCol] << std::endl;

      continue;
    }

    // continuous variables 
    if (indCol < numCols_) {  // variable is model variable

      // print
      fout << "C " << xlp[indCol] << " " << colLowerBound[indCol] <<
	" " << colUpperBound[indCol] << " ";

      // variable lower bound?
      CglMixIntRoundVLB VLB = vlbs_[indCol];
      if (VLB.getVar() != UNDEFINED_) {
	fout << VLB.getVal() << " " << xlp[VLB.getVar()] << " " <<
	  colLowerBound[VLB.getVar()] << " " <<
	  colUpperBound[VLB.getVar()] << " ";
      }
      else {
	fout << "-1 -1 -1 -1 ";
      }

      // variable upper bound?
      CglMixIntRoundVUB VUB = vubs_[indCol];
      if (VUB.getVar() != UNDEFINED_) {
	fout << VUB.getVal() << " " << xlp[VUB.getVar()] << " " <<
	  colLowerBound[VUB.getVar()] << " " <<
	  colUpperBound[VUB.getVar()] << " ";
      }
      else {
	fout << "-1 -1 -1 -1 ";
      }

    }	  
    else {  // variable is slack variable
      // in this case the LB = 0 and the UB = infinity
      // print
      fout << "C " << xlpExtra[indCol-numCols_] << " " << 0.0 <<
	" " << si.getInfinity() << " ";
    }

    fout << std::endl;

  }

  fout << "rhs " << rhsAggregated << std::endl;

  fout << "numColsBack " << numColsBack << std::endl;

  if (hasCut) {
    fout << "CUT: YES" << std::endl;
  }
  else {
    fout << "CUT: NO" << std::endl;
  }

}
// This can be used to refresh any inforamtion
void 
CglMixedIntegerRounding::refreshSolver(OsiSolverInterface * solver)
{
  doneInitPre_ = false;
  // Get integer information
  solver->getColType(true);
}
// Create C++ lines to get to current state
std::string
CglMixedIntegerRounding::generateCpp( FILE * fp) 
{
  CglMixedIntegerRounding other;
  fprintf(fp,"0#include \"CglMixedIntegerRounding.hpp\"\n");
  fprintf(fp,"3  CglMixedIntegerRounding mixedIntegerRounding;\n");
  if (MAXAGGR_!=other.MAXAGGR_)
    fprintf(fp,"3  mixedIntegerRounding.setMAXAGGR_(%d);\n",MAXAGGR_);
  else
    fprintf(fp,"4  mixedIntegerRounding.setMAXAGGR_(%d);\n",MAXAGGR_);
  if (MULTIPLY_!=other.MULTIPLY_)
    fprintf(fp,"3  mixedIntegerRounding.setMULTIPLY_(%d);\n",MULTIPLY_);
  else
    fprintf(fp,"4  mixedIntegerRounding.setMULTIPLY_(%d);\n",MULTIPLY_);
  if (CRITERION_!=other.CRITERION_)
  fprintf(fp,"3  mixedIntegerRounding.setCRITERION_(%d);\n",CRITERION_);
  if (doPreproc_!=other.doPreproc_)
    fprintf(fp,"3  mixedIntegerRounding.setDoPreproc(%d);\n", doPreproc_);
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  mixedIntegerRounding.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  mixedIntegerRounding.setAggressiveness(%d);\n",getAggressiveness());
  return "mixedIntegerRounding";
}
void CglMixedIntegerRounding::setDoPreproc(int value)
{
  if (value != -1 && value != 0 && value != 1) {
    throw CoinError("setDoPrepoc", "invalid value",
		    "CglMixedIntegerRounding2");
  }
  else {
    doPreproc_ = value;
  }  
}

bool CglMixedIntegerRounding::getDoPreproc() const
{
  return (doPreproc_!=0);
}
