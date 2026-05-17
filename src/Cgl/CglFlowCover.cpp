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

#include <cstdlib>
#include <cmath>

#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include "CoinSort.hpp"

#include "CglFlowCover.hpp"

// added #define to get rid of warnings (so uncomment if =true)
//#define CGLFLOW_DEBUG2
//static bool CGLFLOW_DEBUG=false;
#define CGLFLOW_DEBUG 0
#include <iomanip>
//-------------------------------------------------------------------
// Overloaded operator<< for printing VUB and VLB.
//-------------------------------------------------------------------  
std::ostream& operator<<( std::ostream& os, const CglFlowVUB &v ) 
{ 
  os << " VAR = " << v.getVar() << "\t VAL = " << v.getVal() << std::endl; 
  return os; 
}

//-------------------------------------------------------------------
// Determine row types. Find the VUBS and VLBS. 
//-------------------------------------------------------------------  
void 
CglFlowCover::flowPreprocess(const OsiSolverInterface& si)
{
  CoinPackedMatrix matrixByRow(*si.getMatrixByRow());

  int numRows = si.getNumRows();
  int numCols = si.getNumCols();
  
  const char* sense        = si.getRowSense();
  const double* RHS        = si.getRightHandSide();

  const double* coefByRow  = matrixByRow.getElements();
  const int* colInds       = matrixByRow.getIndices();
  const CoinBigIndex* rowStarts     = matrixByRow.getVectorStarts();
  const int* rowLengths    = matrixByRow.getVectorLengths();
  int iRow      = -1; 
  int iCol      = -1;

  numCols_ = numCols;     // Record col and row numbers for copy constructor
  numRows_ = numRows;

  if (rowTypes_ != 0) {
    delete [] rowTypes_; rowTypes_ = 0;
  }
  rowTypes_ = new CglFlowRowType [numRows];// Destructor will free memory
  // Get integer types
  const char * columnType = si.getColType (true);
    
  // Summarize the row type infomation.
  int numUNDEFINED   = 0;
  int numVARUB       = 0;
  int numVARLB       = 0;
  int numVAREQ       = 0;
  int numMIXUB       = 0;
  int numMIXEQ       = 0;
  int numNOBINUB     = 0;
  int numNOBINEQ     = 0;
  int numSUMVARUB    = 0;
  int numSUMVAREQ    = 0;
  int numUNINTERSTED = 0;

  int* ind     = new int [numCols];
  double* coef = new double [numCols];
  for (iRow = 0; iRow < numRows; ++iRow) {
    int rowLen   = rowLengths[iRow];
    char sen     = sense[iRow];
    double rhs   = RHS[iRow];

    CoinDisjointCopyN(colInds + rowStarts[iRow], rowLen, ind);
    CoinDisjointCopyN(coefByRow + rowStarts[iRow], rowLen, coef);
 
    CglFlowRowType rowType = determineOneRowType(si, rowLen, ind, coef, 
						 sen, rhs);

    rowTypes_[iRow] = rowType;

    switch(rowType) {
    case  CGLFLOW_ROW_UNDEFINED:
      ++numUNDEFINED; 
      break;
    case  CGLFLOW_ROW_VARUB:
      ++numVARUB; 
      break;
    case  CGLFLOW_ROW_VARLB:
      ++numVARLB; 
      break;
    case  CGLFLOW_ROW_VAREQ:
      ++numVAREQ; 
      break;
    case  CGLFLOW_ROW_MIXUB:
      ++numMIXUB; 
      break;
    case  CGLFLOW_ROW_MIXEQ:
      ++numMIXEQ; 
      break;
    case  CGLFLOW_ROW_NOBINUB:
      ++numNOBINUB; 
      break;
    case  CGLFLOW_ROW_NOBINEQ:
      ++numNOBINEQ; 
      break;
    case  CGLFLOW_ROW_SUMVARUB:
      ++numSUMVARUB; 
      break;
    case  CGLFLOW_ROW_SUMVAREQ:
      ++numSUMVAREQ; 
      break;
    case  CGLFLOW_ROW_UNINTERSTED:
      ++numUNINTERSTED;
      break;
    default:
      throw CoinError("Unknown row type", "flowPreprocess",
		      "CglFlowCover");
    }
    
  }
  delete [] ind;  ind  = NULL;
  delete [] coef; coef = NULL;

  if(CGLFLOW_DEBUG) {
    std::cout << "The num of rows = "  << numRows        << std::endl;
    std::cout << "Summary of Row Type" << std::endl;
    std::cout << "numUNDEFINED     = " << numUNDEFINED   << std::endl;
    std::cout << "numVARUB         = " << numVARUB       << std::endl;
    std::cout << "numVARLB         = " << numVARLB       << std::endl;
    std::cout << "numVAREQ         = " << numVAREQ       << std::endl;
    std::cout << "numMIXUB         = " << numMIXUB       << std::endl;
    std::cout << "numMIXEQ         = " << numMIXEQ       << std::endl;
    std::cout << "numNOBINUB       = " << numNOBINUB     << std::endl;
    std::cout << "numNOBINEQ       = " << numNOBINEQ     << std::endl;
    std::cout << "numSUMVARUB      = " << numSUMVARUB    << std::endl;
    std::cout << "numSUMVAREQ      = " << numSUMVAREQ    << std::endl;
    std::cout << "numUNINTERSTED   = " << numUNINTERSTED << std::endl;
  }

  //---------------------------------------------------------------------------
  // Setup  vubs_ and vlbs_
  if (vubs_ != 0) { delete [] vubs_; vubs_ = 0; }
  vubs_ = new CglFlowVUB [numCols];      // Destructor will free memory
  if (vlbs_ != 0) { delete [] vlbs_; vlbs_ = 0; }
  vlbs_ = new CglFlowVLB [numCols];      // Destructor will free memory

  for (iCol = 0; iCol < numCols; ++iCol) {   // Initilized in constructor
    vubs_[iCol].setVar(UNDEFINED_);     // but, need redo since may call
    vlbs_[iCol].setVar(UNDEFINED_);     // preprocess(...) more than once
  }
  
  for (iRow = 0; iRow < numRows; ++iRow) {
	
    CglFlowRowType rowType2 = rowTypes_[iRow];
    
    if ( (rowType2 == CGLFLOW_ROW_VARUB) || 
	 (rowType2 == CGLFLOW_ROW_VARLB) || 
	 (rowType2 == CGLFLOW_ROW_VAREQ) )  { 
      
      CoinBigIndex startPos = rowStarts[iRow];
      int index0   = colInds[startPos];
      int index1   = colInds[startPos + 1];
      double coef0 = coefByRow[startPos];
      double coef1 = coefByRow[startPos + 1];
	    
      int    xInd,  yInd;   // x is binary
      double xCoef, yCoef;

      if ( columnType[index0]==1 ) {
	xInd  = index0;   yInd  = index1;
	xCoef = coef0;    yCoef = coef1;
      }
      else {
	xInd  = index1;   yInd  = index0;
	xCoef = coef1;    yCoef = coef0;
      }

      switch (rowType2) {
      case CGLFLOW_ROW_VARUB:       // Inequality: y <= ? * x
	vubs_[yInd].setVar(xInd);
	vubs_[yInd].setVal(-xCoef / yCoef);
	break;
      case CGLFLOW_ROW_VARLB:       // Inequality: y >= ? * x
	vlbs_[yInd].setVar(xInd);
	vlbs_[yInd].setVal(-xCoef / yCoef);
	break;
      case CGLFLOW_ROW_VAREQ:       // Inequality: y >= AND <= ? * x
	vubs_[yInd].setVar(xInd);
	vubs_[yInd].setVal(-xCoef / yCoef);
	vlbs_[yInd].setVar(xInd);
	vlbs_[yInd].setVal(-xCoef / yCoef);
	break;
      default:
	throw CoinError("Unknown row type: impossible", 
			"flowPreprocess", "CglFlowCover");
      }
    }
  }

  if(CGLFLOW_DEBUG) {
    printVubs(std::cout);
  }
}


//-----------------------------------------------------------------------------
// Generate LSGFC cuts
//------------------------------------------------------------------- 
void CglFlowCover::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
				const CglTreeInfo info)
{
  if (getMaxNumCuts() <= 0) return;
    
  if (getNumFlowCuts() >= getMaxNumCuts()) return;

#if 0
  bool preInit = false;
  bool preReso = false;
  si.getHintParam(OsiDoPresolveInInitial, preInit);
  si.getHintParam(OsiDoPresolveInResolve, preReso);

  if (preInit == false &&  preReso == false) { // Do once
    if (doneInitPre_ == false) {   
      flowPreprocess(si);
      doneInitPre_ = true;
    }
  }
  else
#endif
    int numberRowCutsBefore = cs.sizeRowCuts();
    
  flowPreprocess(si);

  CoinPackedMatrix matrixByRow(*si.getMatrixByRow());
  const char* sense = si.getRowSense();
  const double* rhs = si.getRightHandSide();
  const double * colLower = si.getColLower();
  const double * colUpper = si.getColUpper();

  const double* elementByRow = matrixByRow.getElements();
  const int* colInd = matrixByRow.getIndices();
  const CoinBigIndex* rowStart = matrixByRow.getVectorStarts();
  const int* rowLength = matrixByRow.getVectorLengths();
  int numberColumns = si.getNumCols();
  int* ind = new int [numberColumns];
  double* coef = new double [numberColumns];
  int iRow;
  CoinBigIndex iCol;

  CglFlowRowType rType;

  for (iRow = 0; iRow < numRows_; ++iRow) {
    rType = getRowType(iRow);
    if( ( rType != CGLFLOW_ROW_MIXUB ) &&
	( rType != CGLFLOW_ROW_MIXEQ ) &&
	( rType != CGLFLOW_ROW_NOBINUB ) &&
	( rType != CGLFLOW_ROW_NOBINEQ ) &&
	( rType != CGLFLOW_ROW_SUMVARUB ) &&
	( rType != CGLFLOW_ROW_SUMVAREQ ) )
      continue;  

    const CoinBigIndex sta = rowStart[iRow];     // Start position of iRow
    int rowLen = rowLength[iRow]; // iRow length / non-zero elements


    CoinBigIndex lastPos = sta + rowLen;
    double thisRhs = rhs[iRow];
    rowLen=0;
    for (iCol = sta; iCol < lastPos; ++iCol) {
      int jCol=colInd[iCol];
      if (colLower[jCol]<colUpper[jCol]) {
	ind[rowLen]  = jCol;
	coef[rowLen++] = elementByRow[iCol];
      } else {
	thisRhs -= colLower[jCol]*elementByRow[iCol];
      }
    }

    OsiRowCut flowCut1, flowCut2, flowCut3;
    double violation = 0.0;
    bool hasCut = false;

    if (sense[iRow] == 'E') {
      hasCut = generateOneFlowCut(si, rowLen, ind, coef, 'L', 
				  thisRhs, flowCut1, violation);
      if (hasCut)  {                         // If find a cut
	cs.insertIfNotDuplicateAndClean(flowCut1,41);
	incNumFlowCuts();
	if (getNumFlowCuts() >= getMaxNumCuts())
	  break;
      }
      hasCut = false;
      hasCut = generateOneFlowCut(si, rowLen, ind, coef, 'G', 
				  thisRhs, flowCut2, violation);
      if (hasCut)  {
	cs.insertIfNotDuplicateAndClean(flowCut2,42);
	incNumFlowCuts();
	if (getNumFlowCuts() >= getMaxNumCuts())
	  break;
      }
    }
    if (sense[iRow] == 'L' || sense[iRow] == 'G') {
      hasCut = generateOneFlowCut(si, rowLen, ind, coef, sense[iRow], 
				  thisRhs, flowCut3, violation);
      if (hasCut)  {
	cs.insertIfNotDuplicateAndClean(flowCut3,43);
	incNumFlowCuts();
	if (getNumFlowCuts() >= getMaxNumCuts())
	  break;
      }
    }
  }


#ifdef CGLFLOW_DEBUG2
  if(CGLFLOW_DEBUG) {
    std::cout << "\nnumFlowCuts = "<< getNumFlowCuts()  << std::endl;
    std::cout << "CGLFLOW_COL_BINNEG = "<< CGLFLOW_COL_BINNEG  << std::endl;
  }
#endif
  if (!info.inTree&&((info.options&4)==4||((info.options&8)&&!info.pass))) {
    int numberRowCutsAfter = cs.sizeRowCuts();
    for (int i=numberRowCutsBefore;i<numberRowCutsAfter;i++)
      cs.rowCutPtr(i)->setGloballyValid();
  }

  delete [] ind;
  delete [] coef;
}

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglFlowCover::CglFlowCover()
  :
  CglCutGenerator(),
  maxNumCuts_(2000000),
  EPSILON_(1.0e-6),
  UNDEFINED_(-1),
  INFTY_(1.0e30),
  TOLERANCE_(0.05),
  firstProcess_(true),
  numRows_(0),
  numCols_(0),
  numFlowCuts_(0),
  doneInitPre_(false),
  vubs_(0),
  vlbs_(0),
  rowTypes_(0)
{ 
  // DO NOTHING
}


//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglFlowCover::CglFlowCover (const CglFlowCover & source)
  :
  CglCutGenerator(source), 
  maxNumCuts_(source.maxNumCuts_),
  EPSILON_(source.EPSILON_),
  UNDEFINED_(source.UNDEFINED_),
  INFTY_(source.INFTY_),
  TOLERANCE_(source.TOLERANCE_),
  firstProcess_(true),
  numRows_(source.numRows_),
  numCols_(source.numCols_),
  doneInitPre_(source.doneInitPre_)
{ 
  setNumFlowCuts(source.numFlowCuts_);
  if (numCols_ > 0) {
    vubs_ = new CglFlowVUB [numCols_];
    vlbs_ = new CglFlowVLB [numCols_];
    CoinDisjointCopyN(source.vubs_, numCols_, vubs_);
    CoinDisjointCopyN(source.vlbs_, numCols_, vlbs_);
  }
  else {
    vubs_ = 0;
    vlbs_ = 0;
  }
  if (numRows_ > 0) {
    rowTypes_ = new CglFlowRowType [numRows_];
    CoinDisjointCopyN(source.rowTypes_, numRows_, rowTypes_);
  }
  else {
    rowTypes_ = 0;
  }
}


//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglFlowCover::clone() const
{
  return new CglFlowCover(*this);
}

//------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglFlowCover &
CglFlowCover::operator=(const CglFlowCover& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    maxNumCuts_ = rhs.maxNumCuts_;
    EPSILON_ = rhs.EPSILON_;
    UNDEFINED_ = rhs.UNDEFINED_;
    INFTY_ = rhs.INFTY_;
    TOLERANCE_ = rhs.TOLERANCE_;
    numRows_ = rhs.numRows_;
    numCols_ = rhs.numCols_;
    //    numFlowCuts_ = rhs.numFlowCuts_;
    setNumFlowCuts(rhs.numFlowCuts_);
    doneInitPre_ = rhs.doneInitPre_;
    if (numCols_ > 0) {
      vubs_ = new CglFlowVUB [numCols_];
      vlbs_ = new CglFlowVLB [numCols_];
      CoinDisjointCopyN(rhs.vubs_, numCols_, vubs_);
      CoinDisjointCopyN(rhs.vlbs_, numCols_, vlbs_);
    }
    if (numRows_ > 0) {
      rowTypes_ = new CglFlowRowType [numRows_];
      CoinDisjointCopyN(rhs.rowTypes_, numRows_, rowTypes_);
    }
  }
  return *this;
}


//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------  
CglFlowCover::~CglFlowCover ()
{
  if (vubs_ != 0) { delete [] vubs_; vubs_ = 0; }
  if (vlbs_ != 0) { delete [] vlbs_; vlbs_ = 0; }
  if (rowTypes_ != 0) { delete [] rowTypes_; rowTypes_ = 0; } 
}


//-------------------------------------------------------------------
//  Given the model data, a row of the model, and a LP solution, 
//  this function tries to generate a violated lifted simple generalized 
//  flow cover.
//-------------------------------------------------------------------  
bool 
CglFlowCover::generateOneFlowCut( const OsiSolverInterface & si, 
				  const int rowLen,
				  int* ind,
				  double* coef,
				  char sense,
				  double rhs,
				  OsiRowCut& flowCut,
				  double& violation )
{
  bool generated       = false;
  const double* xlp    = si.getColSolution();
  const int numCols    = si.getNumCols();
    
  double* up           = new double [rowLen];
  double* x            = new double [rowLen];
  double* y            = new double [rowLen];
  CglFlowColType* sign = new CglFlowColType [rowLen];
    
  int i, j;  
  double value, LB, UB;
    
  CglFlowVLB VLB;
  CglFlowVUB VUB;
  //CGLFLOW_DEBUG=false;
  bool doLift=true;
  // Get integer types
  const char * columnType = si.getColType ();
  for (i = 0; i < rowLen; ++i) {
    if ( xlp[ind[i]] - floor(xlp[ind[i]]) > EPSILON_ && ceil(xlp[ind[i]]) - xlp[ind[i]] > EPSILON_ )
      break;
  }
  if (i != rowLen) {
    for (int j = 0; j < rowLen; ++j) {
      if ( fabs(coef[j])<=EPSILON_) {
	doLift = false;
	break;
      }
    }
  } else {
    doLift = false;
  }

  if (!doLift)  {
    delete [] sign;
    delete [] up; 
    delete [] x;   
    delete [] y; 
    return generated;
  }

  //-------------------------------------------------------------------------

  if (sense == 'G') flipRow(rowLen, coef, rhs); // flips everything,
  // but the sense
					  

  if(CGLFLOW_DEBUG) {
    std::cout << "***************************" << std::endl;
    std::cout << "Generate Flow cover -- initial constraint, converted to L sense..." << std::endl;
    std::cout << "Rhs = " << rhs << std::endl;
    std::cout << "coef [var_index]" << " -- " <<  "xlp[var_index]" << '\t' << "vub_coef[vub_index] vub_lp_value OR var_index_col_ub" << std::endl;
   
    for(int iD = 0; iD < rowLen; ++iD) {
      VUB = getVubs(ind[iD]);

      std::cout << std::setw(5) << coef[iD] << "["  << std::setw(5)  << ind[iD] << "] -- " 
		<< std::setw(20) << xlp[ind[iD]] << '\t';
      if (VUB.getVar() != UNDEFINED_) {  
	std::cout << std::setw(20) << VUB.getVal() << "[" << std::setw(5) << VUB.getVar() << "]" 
		  << std::setw(20) << xlp[VUB.getVar()] << std::endl; 
      }
      else
	std::cout << std::setw(20) << si.getColUpper()[ind[iD]] << "       " << std::setw(20) << 1.0 << std::endl;
	
    }
  }

  //-------------------------------------------------------------------------
  // Generate conservation inequality and capacity equalities from 
  // the given row.
  
  for (i = 0; i < rowLen; ++i) {
	
    VLB = getVlbs(ind[i]);
    LB = ( VLB.getVar() != UNDEFINED_ ) ? 
      VLB.getVal() : si.getColLower()[ind[i]];

    VUB = getVubs(ind[i]);
    UB = ( VUB.getVar() != UNDEFINED_ ) ? 
      VUB.getVal() : si.getColUpper()[ind[i]];

    if (LB < -EPSILON_) {   // Only consider rows whose variables are all
      delete [] sign;       // non-negative (LB>= 0). 
      delete [] up; 
      delete [] x;   
      delete [] y;  
      return generated;     
    }

    if ( columnType[ind[i]]==1 ) {   // Binary variable
      value = coef[i];
      if (value > 0.0)
	sign[i] = CGLFLOW_COL_BINPOS;
      else {
	sign[i] = CGLFLOW_COL_BINNEG;
	value = -value;
      }    
      up[i] = value;
      x[i] =  xlp[ind[i]];
      y[i] = value * x[i];
    }
    else {   
      value = coef[i];
      if (value > 0.0)
	sign[i] = CGLFLOW_COL_CONTPOS;
      else {
	sign[i] = CGLFLOW_COL_CONTNEG;
	value = -value;
      }
      up[i] = value* UB;
      x[i] = (VUB.getVar() != UNDEFINED_) ? xlp[VUB.getVar()] : 1.0;
      y[i] = value * xlp[ind[i]];
    }
  }

  //-------------------------------------------------------------------------
  // Find a initial cover (C+, C-) in (N+, N-)
  double  knapRHS   = rhs;
  double  tempSum   = 0.0;
  double  tempMin   = INFTY_;
  CglFlowColCut *    candidate = new CglFlowColCut [rowLen];
  CglFlowColCut *    label     = new CglFlowColCut [rowLen];
  double* ratio     = new double [rowLen];
  int t = -1;
  for (i = 0; i < rowLen; ++i) {
    candidate[i] = label[i] = CGLFLOW_COL_OUTCUT;
    ratio[i] = INFTY_;
	
    switch(sign[i]) {
    case CGLFLOW_COL_CONTPOS:
    case CGLFLOW_COL_BINPOS:
      if( y[i] > EPSILON_ ) {
	ratio[i] = (1.0 - x[i]) / up[i];
	if( y[i] > up[i] * x[i] - EPSILON_ ) {       // Violated
	  candidate[i] = CGLFLOW_COL_PRIME;
	  tempSum += up[i];
	}
	else {
	  candidate[i] = CGLFLOW_COL_SECONDARY;
	}
      }
      break;
    case CGLFLOW_COL_CONTNEG:
    case CGLFLOW_COL_BINNEG:
      if( up[i] > ( (1.0 - EPSILON_) * INFTY_ ) ) {  // UB is infty
	label[i] = CGLFLOW_COL_INCUT;
      }
      else {
	knapRHS += up[i];
	if( y[i] < up[i] ) {
	  candidate[i] = CGLFLOW_COL_PRIME;
	  assert (up[i]);
	  ratio[i] = x[i] / up[i];
	  tempSum += up[i];
	}
      }
      break;
    }    
  }
    
  double diff, tempD, lambda;
  int xID = -1;
  if (knapRHS >1.0e10) {
    if(CGLFLOW_DEBUG) {
      std::cout << "knapsack RHS too large. RETURN." << std::endl; 
    }
    delete [] sign;                              
    delete [] up; 
    delete [] x;   
    delete [] y;  
    delete [] candidate;
    delete [] label;
    delete [] ratio;
    return generated;
  }

  while (tempSum < knapRHS + EPSILON_) { // Not a cover yet
    diff = INFTY_;
    for (i = 0; i < rowLen; ++i) {
      if (candidate[i] == CGLFLOW_COL_SECONDARY) {
	tempD = up[i] * x[i] - y[i];
	if (tempD < diff - EPSILON_) {
	  diff = tempD;
	  xID = i;
	}
      }
    }
    
    if( diff > (1.0 - EPSILON_) * INFTY_  ) {   // NO cover exits.
      delete [] sign;                              
      delete [] up; 
      delete [] x;   
      delete [] y;  
      delete [] candidate;
      delete [] label;
      delete [] ratio;
      return generated;
    }
    else {
      tempSum += up[xID];
      candidate[xID] = CGLFLOW_COL_PRIME;
    }
  }

  // Solve the knapsack problem to get an initial cover
  tempSum = 0.0;
  for (i = 0; i < rowLen; ++i) {
    if (candidate[i] == CGLFLOW_COL_PRIME && ratio[i] < EPSILON_) {
      //Zero ratio
      label[i] = CGLFLOW_COL_INCUT;
      tempSum += up[i];
    }
  }
  
  while (tempSum < knapRHS + EPSILON_) {
    tempMin = INFTY_;
    xID=-1;
    for (i = 0; i < rowLen; i++) {   // Search the col with  minimum ratio
      if (candidate[i] == CGLFLOW_COL_PRIME && label[i] == 0 && 
	  ratio[i] < tempMin) {
	tempMin = ratio[i];  
	xID = i; 
      }
    }
    if (xID>=0) {
      label[xID] = CGLFLOW_COL_INCUT;
      tempSum += up[xID];
    } else {
      if(CGLFLOW_DEBUG) {
	std::cout << "knapsack RHS too large B. RETURN." << std::endl; 
      }
      delete [] sign;                              
      delete [] up; 
      delete [] x;   
      delete [] y;  
      delete [] candidate;
      delete [] label;
      delete [] ratio;
      return generated;
    }
  }
  
  // Reduce to a minimal cover
  for (i = 0; i < rowLen; ++i) {
    if (label[i] == CGLFLOW_COL_INCUT && ratio[i] > EPSILON_) {
      if (tempSum - up[i] > knapRHS + EPSILON_) {
	label[i] = CGLFLOW_COL_OUTCUT;
	tempSum -= up[i];
      }
    }
  }
  for (i = 0; i < rowLen; ++i) {
    if (label[i] == CGLFLOW_COL_INCUT && ratio[i] < EPSILON_) {
      if (tempSum - up[i] > knapRHS + EPSILON_) {
	label[i] = CGLFLOW_COL_OUTCUT;
	tempSum -= up[i];
      }
    }
  }
    
  // Due to the way to handle N-
  for(i = 0; i < rowLen; ++i) {
    if( sign[i] < 0 ) 
      label[i] = label[i]==CGLFLOW_COL_OUTCUT?CGLFLOW_COL_INCUT:CGLFLOW_COL_OUTCUT;
  }

  // No cover, no cut. 
  bool emptyCover = true; 
  for (i = 0; i < rowLen; ++i) {
    if (label[i] == CGLFLOW_COL_INCUT) {
      emptyCover = false; 
      break;
    }
  }
  if (emptyCover) {	
    if(CGLFLOW_DEBUG) {
      std::cout << "No cover. RETURN." << std::endl; 
    }
    delete [] sign;                              
    delete [] up; 
    delete [] x;   
    delete [] y;  
    delete [] candidate;
    delete [] label;
    delete [] ratio;
    return generated;  
  }

  lambda = tempSum - knapRHS;

  if(CGLFLOW_DEBUG) {
    //double sum_mj_Cplus = 0.0;
    //double sum_mj_Cminus= 0.0;
    // double checkLambda; // variable not used anywhere (LL)
    // print out the knapsack variables
    std::cout << "Knapsack Cover: C+" << std::endl;
    for (i = 0; i < rowLen; ++i) { 
      if ( label[i] == CGLFLOW_COL_INCUT && sign[i] > 0 ) {
	std::cout << ind[i] << '\t' << up[i] << std::endl;
	//sum_mj_Cplus += up[i];
      }
    } 
    std::cout << "Knapsack Cover: C-" << std::endl;
    for (i = 0; i < rowLen; ++i) { 
      if ( label[i] == CGLFLOW_COL_INCUT && sign[i] < 0 ) {
	std::cout << ind[i] << '\t' << up[i] << std::endl;
	//sum_mj_Cminus += up[i];
      }
    }

    // rlh: verified "lambda" is lambda in the paper.
    // lambda = (sum coefficients in C+) - (sum of VUB
    // coefficients in C-) - rhs-orig-constraint
    std::cout << "lambda = " << lambda << std::endl;
  }

  //-------------------------------------------------------------------------
  // Generate a violated SGFC

  //int numCMinus = 0;
  int numPlusPlus = 0;
  double* rho     = new double [rowLen];
  double* xCoef   = new double [rowLen]; 
  double* yCoef   = new double [rowLen];
  double cutRHS   = rhs;
  double temp     = 0.0;
  double sum      = 0.0;
  double minPlsM  = INFTY_;
  double minNegM  = INFTY_;

  for(i = 0; i < rowLen; ++i) {
    rho[i]   = 0.0;
    xCoef[i] = 0.0;
    yCoef[i] = 0.0;
  }
    
  // Project out variables in C-
  // d^' = d + sum_{i in C^-} m_i. Now cutRHS = d^'
  for (i = 0; i < rowLen; ++i) { 
    if ( label[i] == CGLFLOW_COL_INCUT && sign[i] < 0 ) {
      cutRHS += up[i];
      //++numCMinus;
    }
  }

  // (1) Compute the coefficients of the simple generalized flow cover
  // (2) Compute minPlsM, minNegM and sum
  //
  // sum = sum_{i in C+\C++} m_i + sum_{i in L--} m_i = m. Page 15.
  // minPlsM = min_{i in C++} m_i
  // minNegM = min_{i in L-} m_i

  temp = cutRHS;
    
  for (i = 0; i < rowLen; ++i) {
    if (label[i] == CGLFLOW_COL_INCUT  && sign[i] > 0) { // C+
      yCoef[i] = 1.0;
      if ( up[i] > lambda + EPSILON_ ) { // C++
	++numPlusPlus;
	xCoef[i] = lambda - up[i];
	cutRHS += xCoef[i];
	if( up[i] < minPlsM ) {
	  minPlsM = up[i];
	}
      }
      else {  // C+\C++
	xCoef[i] = 0.0;  // rlh: is this necesarry? (xCoef initialized to zero)
	sum += up[i];
      } 
    }
	
    if (label[i] != CGLFLOW_COL_INCUT && sign[i] < 0) { // N-\C-
      temp += up[i];
      if ( up[i] > lambda) {      // L-
	if(CGLFLOW_DEBUG) {
	  std::cout << "Variable " << ind[i] << " is in L-" << std::endl;
	}
	yCoef[i] = 0.0;
	xCoef[i] = -lambda;
	label[i] = CGLFLOW_COL_INLMIN;
	if ( up[i] < minNegM ) { 
	  minNegM = up[i];
	}
      }
      else  {        // L--
	if(CGLFLOW_DEBUG) {
	  std::cout << "Variable " << ind[i] << " is in L-- " << std::endl;
	}
	yCoef[i] = -1.0;
	xCoef[i] = 0.0; // rlh: is this necesarry? (xCoef initialized to zero)
	label[i] = CGLFLOW_COL_INLMINMIN;
	sum += up[i];
      }
    }
  }
   
  // Sort the upper bounds (m_i) of variables in C++ and L-.

  int     ix;
  int     index  = 0;
  double* mt     = new double [rowLen];
  double* M      = new double [rowLen + 1];
  // order to look at variables
  int * order = new int [rowLen];
  int nLook=0;
  for (int i = 0; i < rowLen; ++i) {
    if ( (label[i] == CGLFLOW_COL_INCUT && sign[i] > 0) || 
	 label[i] == CGLFLOW_COL_INLMIN ) {     //  C+ || L- 
      // possible
      M[nLook]=-up[i];
      order[nLook++]=i;
    }
  }
  CoinSort_2(M,M+nLook,order);
  int kLook=0;
  
  while (kLook<nLook) {
    ix = UNDEFINED_;
    i = order[kLook];
    kLook++;
    if ( (label[i] == CGLFLOW_COL_INCUT && sign[i] > 0) || 
	 label[i] == CGLFLOW_COL_INLMIN ) {     //  C+ || L- 
      if ( up[i] > lambda ) {       // C++ || L-(up[i] > lambda)
	ix = i;
      }
    }
      
    if( ix == UNDEFINED_ )  break;
      
    mt[index++] = up[ix];  // Record m_i in C++ and L-(not all) in descending order.
	
    if( label[ix] == CGLFLOW_COL_INLMIN )  
      label[ix] = CGLFLOW_COL_INLMINDONE;
    else
      label[ix] = CGLFLOW_COL_INCUTDONE;
  }
  //printf("mins %g %g\n",minNegM,minPlsM);
  if( index == 0 || numPlusPlus == 0) {
    // No column in C++ and L-(not all). RETURN.
    if(CGLFLOW_DEBUG) {
      std::cout << "index = 0. RETURN." << std::endl; 
    }
    delete [] sign;
    delete [] up; 
    delete [] x;   
    delete [] y;  
    delete [] candidate;
    delete [] label;
    delete [] ratio;
    delete [] rho;
    delete [] xCoef;
    delete [] yCoef;
    delete [] mt; 
    delete [] M; 
    delete [] order;
    return generated;
  }

  for ( i = 0; i < rowLen; i++ ) {
    switch( label[i] ) {
    case  CGLFLOW_COL_INCUTDONE:
      label[i] = CGLFLOW_COL_INCUT;
      break;
    case  CGLFLOW_COL_INLMIN:
    case  CGLFLOW_COL_INLMINDONE:
    case  CGLFLOW_COL_INLMINMIN:
      label[i] = CGLFLOW_COL_OUTCUT;
      break;
    case CGLFLOW_COL_INCUT:
    case CGLFLOW_COL_OUTCUT:
    case CGLFLOW_COL_PRIME:
    case CGLFLOW_COL_SECONDARY:
      break;
    }
  }
    
  /* Get t */
  t = 0;
  for ( i = 0; i < index; ++i ) {
    if ( mt[i] < minPlsM ) {
      t = i;
      break;
    } 
  }

  if (i == index) {
    t = index;
  }
    
  /* Compute M_i */
  M[0] = 0.0;
  for ( i = 1; i <= index; ++i ) {
    M[i] = M[(i-1)] + mt[(i-1)];
    if(CGLFLOW_DEBUG) {
      std::cout << "t = " << t << std::endl; 
      std::cout << "mt[" << std::setw(5) << (i-1) << "]=" << std::setw(2) << ", M[" << std::setw(5) << i << "]=" << std::setw(20) << M[i] << std::endl;
    }
  }
  // Exit if very big M
  if (M[index]>1.0e30) { // rlh: should test for huge col UB earler 
    // no sense doing all this work in that case.
    if(CGLFLOW_DEBUG) {
      std::cout << "M[index]>1.0e30. RETURN." << std::endl; 
      delete [] sign;
      delete [] up; 
      delete [] x;   
      delete [] y;  
      delete [] candidate;
      delete [] label;
      delete [] ratio;
      delete [] rho;
      delete [] xCoef;
      delete [] yCoef;
      delete [] mt; 
      delete [] M; 
      delete [] order;
      return generated;
    }
  }

  /* Get ml */
  double ml = std::min(sum, lambda);
  if(CGLFLOW_DEBUG) {
    // sum = sum_{i in C+\C++} m_i + sum_{i in L--} m_i = m. Page 15.
    std::cout << "ml = std::min(m, lambda) = std::min(" << sum << ", " << lambda << ") =" << ml << std::endl; 
  }
  /* rho_i = max[0, m_i - (minPlsM - lamda) - ml */
  if (t < index ) { /* rho exits only for t <= index-1 */
    value = (minPlsM - lambda) + ml;
    for (i = t; i < index; ++i) {
      rho[i] =  std::max(0.0, mt[i] - value);
      if(CGLFLOW_DEBUG) {
	std::cout << "rho[" << std::setw(5) << i << "]=" << std::setw(20) << rho[i] << std::endl;
      }
    }
  }
  // Calculate the violation
  violation = -cutRHS;
  for ( i = 0; i < rowLen; ++i ) {
#ifdef CGLFLOW_DEBUG2
    if(CGLFLOW_DEBUG) {
      std::cout << "i = " << i << " ind = " << ind[i] << " sign = " 
		<< sign[i] 
		<< " coef = " << coef[i] << " x = " << x[i] << " xCoef = " 
		<< xCoef[i] << " y = " << y[i] << " yCoef = " << yCoef[i] 
		<< " up = " << up[i] << " label = " << label[i] << std::endl;
    }
#endif
    violation += y[i] * yCoef[i] + x[i] * xCoef[i];
  }

  if(CGLFLOW_DEBUG) {
    std::cout << "violation = " << violation << std::endl;
  }
  //  double violationBeforeLift=violation; // variable not used anywhere (LL)
  if(doLift && fabs(violation) > TOLERANCE_ ) {  // LIFTING
    double estY, estX;
    double movement = 0.0;
    double dPrimePrime = temp + cutRHS; 
    bool lifted = false;
    for( i = 0; i < rowLen; ++i ) {
      if ( (label[i] != CGLFLOW_COL_INCUT) && (sign[i] > 0) ) {/* N+\C+*/
	lifted = liftPlus(estY, estX,
			  index, up[i],
			  lambda,
			  y[i], x[i], 
			  dPrimePrime, M);
	    
	xCoef[i] = -estX;
	yCoef[i] = estY;
	if(CGLFLOW_DEBUG) {
	  if (lifted) {
	    printf("Success: Lifted col %i (up_i=%f,yCoef[i]=%f,xCoef[i]=%f) in N+\\C+\n", 
		   ind[i], up[i], yCoef[i], xCoef[i]);
	  }
	  else {
	    printf("Failed to Lift col %i (m_i=%f) in N+\\C+\n", 
		   ind[i], up[i]);
	  }       
	}
      }
      if (label[i] == CGLFLOW_COL_INCUT && sign[i] < 0) { 
	/* C- */
	liftMinus(movement, t,
		  index, up[i], 
		  dPrimePrime, 
		  lambda, ml,
		  M, rho);
                
	if(movement > EPSILON_) {
	  if(CGLFLOW_DEBUG) {
	    printf("Success: Lifted col %i in C-, movement=%f\n", 
		   ind[i], movement);
	  }
	  lifted = true;
	  xCoef[i] = -movement;
	  cutRHS -= movement;
	}
	else {
	  if(CGLFLOW_DEBUG) {
	    printf("Failed to Lift col %i in C-, g=%f\n",
		   ind[i], movement);
	  }
	}
      }
    }
  }
  //-------------------------------------------------------------------

    
  // Calculate the violation
  violation = -cutRHS;
  for ( i = 0; i < rowLen; ++i ) {
#ifdef CGLFLOW_DEBUG2
    if(CGLFLOW_DEBUG) {
      std::cout << "i = " << i << " ind = " << ind[i] << " sign = " 
		<< sign[i] 
		<< " coef = " << coef[i] << " x = " << x[i] << " xCoef = " 
		<< xCoef[i] << " y = " << y[i] << " yCoef = " << yCoef[i] 
		<< " up = " << up[i] << " label = " << label[i] << std::endl;
    }
#endif
    violation += y[i] * yCoef[i] + x[i] * xCoef[i];
  }

  if(CGLFLOW_DEBUG) {
    std::cout << "violation = " << violation << std::endl;
  }
    
  int     cutLen     = 0;
  int*    cutInd     = 0;
  double* cutCoef    = 0;

  // If violated, transform the inequality back to original system
  if ( violation > TOLERANCE_ ) {
    cutLen = 0;
    cutInd  = new int [3*numCols];
    cutCoef = new double [3*numCols];
      
	  assert (cutLen<numCols);
    for ( i = 0; i < rowLen; ++i )  {
      VUB = getVubs(ind[i]);
      
      if ( ( sign[i] == CGLFLOW_COL_CONTPOS ) || 
	   ( sign[i] == CGLFLOW_COL_CONTNEG ) ) {

	if ( fabs( yCoef[i] ) > EPSILON_ ) {
		    
	  if ( sign[i] == CGLFLOW_COL_CONTPOS ) 
	    cutCoef[cutLen] = coef[i] * yCoef[i];
	  else 
	    cutCoef[cutLen] = -coef[i] * yCoef[i];
	  cutInd[cutLen++] = ind[i];
	}

	if ( fabs( xCoef[i] ) > EPSILON_ ) {
	  if ( VUB.getVar() != UNDEFINED_ ) {
	    cutCoef[cutLen] = xCoef[i];
	    cutInd[cutLen++] = VUB.getVar();
	  }
	  else
	    cutRHS -= xCoef[i];
	}
      }
            
      if ( ( sign[i] == CGLFLOW_COL_BINPOS ) || 
	   ( sign[i] == CGLFLOW_COL_BINNEG ) ) {
	if (fabs(yCoef[i]) > EPSILON_ || fabs(xCoef[i]) > EPSILON_) {
	  if (sign[i] == CGLFLOW_COL_BINPOS) 
	    cutCoef[cutLen] = coef[i] * yCoef[i] + xCoef[i];
	  else 
	    cutCoef[cutLen] = -coef[i] * yCoef[i] + xCoef[i];
	  cutInd[cutLen++] = ind[i];
	}
      }
    }
#if 1
    assert (cutLen);
    CoinShortSort_2(cutInd,cutInd+cutLen,cutCoef);
    j=0;
    int lastInd=cutInd[0];
    double lastCoef=cutCoef[0];
    for ( i = 1; i < cutLen+1; ++i ) {
      if (i==cutLen||cutInd[i]>lastInd) {
	if ( fabs(lastCoef) >= EPSILON_ ) {
	  cutCoef[j]=lastCoef;
	  cutInd[j++]=lastInd;
	}
	lastCoef = cutCoef[i];
	if (i<cutLen)
	  lastInd=cutInd[i];
      } else {
	lastCoef += cutCoef[i];
      }
    }
#else
    for ( i = 0; i < cutLen; ++i ) {
      for ( j = 0; j < i; j++ ) {
	if ( cutInd[j] == cutInd[i] ) { /* Duplicate*/
	  cutCoef[j] += cutCoef[i];
	  cutInd[i] = -1;
	}
      }
    }

    for ( j = 0, i = 0; i < cutLen; ++i ) {
      if ( ( cutInd[i] == -1 ) || ( fabs( cutCoef[i]) < EPSILON_ ) ){
	/* Small coeff*/
      }
      else {
	cutCoef[j] = cutCoef[i];
	cutInd[j] = cutInd[i];
	j++;
      }
    }
#endif
    cutLen = j;
    // Skip if no elements ? - bug somewhere
    if (cutLen == 0) {
        if (xCoef)
            delete[] xCoef;
        if (cutCoef)
            delete[] cutCoef;
        if (yCoef)
            delete[] yCoef;
        if (M)
            delete[] M;
        if (label)
            delete[] label;
        if (sign)
            delete[] sign;
        if (up)
            delete[] up;
        if (candidate)
            delete[] candidate;
        if (y)
            delete[] y;
        if (x)
            delete[] x;
        if (mt)
            delete[] mt;
        if (rho)
            delete[] rho;
        if (order)
            delete[] order;
        return false;
    }
        
    // Recheck the violation.
    double saveViolation = violation;
    violation = 0.0;
    for (i = 0; i < cutLen; ++i) 
      violation += cutCoef[i] * xlp[cutInd[i]];
    
    violation -= cutRHS;
    assert (fabs(violation-saveViolation)<1.0e-2);

    if ( violation > TOLERANCE_ ) {
      flowCut.setRow(cutLen, cutInd, cutCoef);
      flowCut.setLb(-1.0 * si.getInfinity());
      flowCut.setUb(cutRHS);
      flowCut.setEffectiveness(violation);
      generated = true;

      if(CGLFLOW_DEBUG) {
	std::cout << "generateOneFlowCover(): Found a cut" << std::endl;
      }
    }
    else {
      if(CGLFLOW_DEBUG) {
	std::cout << "generateOneFlowCover(): Lost a cut" << std::endl;
      }
    }
  }

  //-------------------------------------------------------------------------
  delete [] sign;
  delete [] up; 
  delete [] x;   
  delete [] y;  
  delete [] candidate;
  delete [] label;
  delete [] ratio;
  delete [] rho;
  delete [] xCoef;
  delete [] yCoef;
  delete [] mt; 
  delete [] M; 
  delete [] order;
  delete [] cutInd;
  delete [] cutCoef;
    
  return generated;
}


//-------------------------------------------------------------------
// Flip a row from ">=" to "<=", and vice versa. 
//-------------------------------------------------------------------
void 
CglFlowCover::flipRow(int rowLen, double* coef, double& rhs) const
{
  for(int i = 0; i < rowLen; ++i) coef[i] = -coef[i]; 
  rhs = -rhs;  
}

//-------------------------------------------------------------------
// Flip a row from ">=" to "<=", and vice versa. Have 'sense'.
//-------------------------------------------------------------------
void 
CglFlowCover::flipRow(int rowLen, double* coef, char& sen,
		      double& rhs) const
{
  for(int i = 0; i < rowLen; ++i) coef[i] = -coef[i]; 
  sen = (sen == 'G') ?  'L' : 'G';
  rhs = -rhs;  
}

//-------------------------------------------------------------------
// Determine the type of a given row 
//-------------------------------------------------------------------
CglFlowRowType
CglFlowCover::determineOneRowType(const OsiSolverInterface& si,
				  int rowLen, int* ind, 
				  double* coef, char sense, 
				  double rhs) const
{
  if (rowLen == 0) 
    return CGLFLOW_ROW_UNDEFINED;
  if (sense == 'R')
    return CGLFLOW_ROW_UNINTERSTED; // Could be fixed
    
  CglFlowRowType rowType = CGLFLOW_ROW_UNDEFINED;
  // Get integer types
  const char * columnType = si.getColType ();
    
  int  numPosBin = 0;      // num of positive binary variables
  int  numNegBin = 0;      // num of negative binary variables
  int  numBin    = 0;      // num of binary variables
  int  numPosCol = 0;      // num of positive variables
  int  numNegCol = 0;      // num of negative variables
  int  i;
  bool flipped = false;

  // Range row will only consider as 'L'
  if (sense == 'G') {        // Transform to " <= "
    flipRow(rowLen, coef, sense, rhs);                
    flipped = true;
  }
    
  // Summarize the variable types of the given row.
  for ( i = 0; i < rowLen; ++i ) {
    if ( coef[i] < -EPSILON_ ) {
      ++numNegCol;
      if( columnType[ind[i]]==1 )
	++numNegBin;
    }
    else {
      ++numPosCol;
      if( columnType[ind[i]]==1 )
	++numPosBin;    
    }
  }
  numBin = numNegBin + numPosBin;

  if(CGLFLOW_DEBUG) {
    std::cout << "numNegBin = " << numNegBin << std::endl;
    std::cout << "numPosBin = " << numPosBin << std::endl;
    std::cout << "numBin = " << numBin << std::endl;
    std::cout << "rowLen = " << rowLen << std::endl;
  }
    
    
  //------------------------------------------------------------------------
  // Classify row type based on the types of variables.
    
  // All variables are binary. NOT interested in this type of row right now
  if (numBin == rowLen) 
    rowType = CGLFLOW_ROW_UNINTERSTED;

  // All variables are NOT binary
  if (rowType == CGLFLOW_ROW_UNDEFINED && numBin == 0) {
    if (sense == 'L')
      rowType = CGLFLOW_ROW_NOBINUB;
    else 
      rowType = CGLFLOW_ROW_NOBINEQ;
  }

  // There are binary and other types of variables   
  if (rowType == CGLFLOW_ROW_UNDEFINED) {  
    if ((rhs < -EPSILON_) || (rhs > EPSILON_) || (numBin != 1)) {
      if (sense == 'L')
	rowType = CGLFLOW_ROW_MIXUB;
      else 
	rowType = CGLFLOW_ROW_MIXEQ;
    }
    else {                               // EXACTLY one binary
      if (rowLen == 2) {               // One binary and one other type
	if (sense == 'L') {
	  if (numNegCol == 1 && numNegBin == 1)
	    rowType = CGLFLOW_ROW_VARUB;
	  if (numPosCol == 1 && numPosBin == 1)
	    rowType = CGLFLOW_ROW_VARLB;
	}
	else
	  rowType = CGLFLOW_ROW_VAREQ;
      }
      else {               // One binary and 2 or more other types
	if (numNegCol==1 && numNegBin==1) {// Binary has neg coef and 
	  if (sense == 'L')  // other are positive
	    rowType = CGLFLOW_ROW_SUMVARUB;
	  else
	    rowType = CGLFLOW_ROW_SUMVAREQ;
	}
      }
    }
  }
  
  // Still undefined
  if (rowType == CGLFLOW_ROW_UNDEFINED) {
    if (sense == 'L') 
      rowType = CGLFLOW_ROW_MIXUB;
    else
      rowType = CGLFLOW_ROW_MIXEQ;
  }
  if (flipped == true) {
    flipRow(rowLen, coef, sense, rhs);                
  }

  return rowType;
}

/*===========================================================================*/

void
CglFlowCover::liftMinus(double &movement, /* Output */ 
			int t,
			int r,
			double z,
			double dPrimePrime, 
			double lambda,
			double ml,
			double *M,
			double *rho) const
{
  int i;
  movement = 0.0;
    
  if (z > dPrimePrime) {
    movement = z - M[r] + r * lambda;
  }
  else {
    for (i = 0; i < t; ++i) {
      if ( (z >= M[i]) && (z <= M[(i+1)] - lambda) ) {
	movement = i * lambda;
	return;
      }
    }
        
    for (i = 1; i < t; ++i) {
      if ( (z >= M[i] - lambda) && (z <= M[i]) ) {
	movement = z - M[i] + i * lambda;
	return;
      }
    }
        
    for (i = t; i < r; ++i) {
      if ( (z >= M[i] - lambda) && (z <= M[i] - lambda + ml + rho[i]) ) {
	movement = z - M[i] + i * lambda;
	return;
      }
    }
        
    for (i = t; i < r; ++i) {
      if ( (z >= M[i]-lambda+ml+rho[i]) && (z <= M[(i+1)]-lambda) ) {
	movement = i * lambda;
	return;
      }
    }
    
    if ((z >= M[r] - lambda) && z <= dPrimePrime) {
      movement = z - M[r] + r * lambda;
    }
  }
    
}

/*===========================================================================*/

bool
CglFlowCover::liftPlus(double &alpha, 
		       double &beta,
		       int r,
		       double m_j, 
		       double lambda,
		       double y_j,
		       double x_j,
		       double dPrimePrime,
		       double *M) const
{
  int i;
  bool status = false;  /* Default: fail to lift */
  double value;
  alpha = 0.0;
  beta = 0.0;
    
  if (m_j > M[r] - lambda + EPSILON_) {
    if (m_j < dPrimePrime - EPSILON_) {
      if ((m_j > (M[r] - lambda)) && (m_j <= M[r])){ /* FIXME: Test */
	value = y_j - x_j * (M[r] - r * lambda);
                
	/* FIXME: Is this "if" useful */
	if (value > 0.0) {
	  status = true;
	  alpha = 1.0;
	  beta = M[r] - r * lambda;
	  if(CGLFLOW_DEBUG) {
	    printf("liftPlus:1: value=%f, alpah=%f, beta=%f\n",
		   value, alpha,beta);
	  }
	}
	else {
	  if(CGLFLOW_DEBUG) {
	    printf("liftPlus:1: value=%f, become worst\n",value);
	  }
	}
      }
    }
    else {    
      if(CGLFLOW_DEBUG) {
	printf("liftPlus:1: too big number\n");
      }
    }
  }
  else {
    for (i = 1; i <= r; ++i) {
      if ((m_j > (M[i] - lambda)) && (m_j <= M[i])){ /* FIXME: Test */
                
	value = y_j - x_j * (M[i] - i * lambda);

	/* FIXME: Is this "if" useful */
	if (value > 0.0) {
	  status = true;
	  alpha = 1.0;
	  beta = M[i] - i * lambda;
	  if(CGLFLOW_DEBUG) {
	    printf("liftPlus:2: value=%f, alpah=%f, beta=%f\n",
		   value, alpha, beta);
	  }
	}
	else {
	  if(CGLFLOW_DEBUG) {
	    printf("liftPlus:2: value=%f, become worst\n",value);
	  }
	}
	return status;
      }
    }
  }
    
  return status;
}
// Create C++ lines to get to current state
std::string
CglFlowCover::generateCpp( FILE * fp) 
{
  CglFlowCover other;
  fprintf(fp,"0#include \"CglFlowCover.hpp\"\n");
  fprintf(fp,"3  CglFlowCover flowCover;\n");
  if (maxNumCuts_!=other.maxNumCuts_)
    fprintf(fp,"3  flowCover.setMaxNumCuts(%d);\n",maxNumCuts_);
  else
    fprintf(fp,"4  flowCover.setMaxNumCuts(%d);\n",maxNumCuts_);
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  flowCover.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  flowCover.setAggressiveness(%d);\n",getAggressiveness());
  return "flowCover";
}

