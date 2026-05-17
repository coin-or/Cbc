// Last edit: 02/05/2013
//
// Name:     CglGMI.cpp
// Author:   Giacomo Nannicini
//           Singapore University of Technology and Design, Singapore
//           email: nannicini@sutd.edu.sg
// Date:     11/17/09
//-----------------------------------------------------------------------------
// Copyright (C) 2009, Giacomo Nannicini.  All Rights Reserved.

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>
#include <climits>

#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CoinFactorization.hpp"
#include "CglGMI.hpp"
#include "CoinFinite.hpp"
#include "CoinRational.hpp"

//-------------------------------------------------------------------
// Generate GMI cuts
//------------------------------------------------------------------- 

/***************************************************************************/
CglGMI::CglGMI() : 
  CglCutGenerator(),
  param(),
  nrow(0),
  ncol(0),
  colLower(NULL),
  colUpper(NULL),
  rowLower(NULL),
  rowUpper(NULL),
  rowRhs(NULL),
  isInteger(NULL),
  cstat(NULL),
  rstat(NULL),
  solver(NULL),
  xlp(NULL),
  rowActivity(NULL),
  byRow(NULL),
  byCol(NULL),
  f0(0.0),
  f0compl(0.0),
  ratiof0compl(0.0)
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
  ,
  trackRejection(false),
  fracFail(0),
  dynFail(0),
  violFail(0),
  suppFail(0),
  scaleFail(0),
  numGeneratedCuts(0)
#endif
{

}

/***************************************************************************/
CglGMI::CglGMI(const CglGMIParam &parameters) : 
  CglCutGenerator(),
  param(parameters),
  nrow(0),
  ncol(0),
  colLower(NULL),
  colUpper(NULL),
  rowLower(NULL),
  rowUpper(NULL),
  rowRhs(NULL),
  isInteger(NULL),
  cstat(NULL),
  rstat(NULL),
  solver(NULL),
  xlp(NULL),
  rowActivity(NULL),
  byRow(NULL),
  byCol(NULL),
  f0(0.0),
  f0compl(0.0),
  ratiof0compl(0.0)
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
  ,
  trackRejection(false),
  fracFail(0),
  dynFail(0),
  violFail(0),
  suppFail(0),
  scaleFail(0),
  numGeneratedCuts(0)
#endif
{

}

/***************************************************************************/
CglGMI::CglGMI(const CglGMI& rhs) : 
  CglCutGenerator(rhs),
  param(rhs.param),
  nrow(rhs.nrow),
  ncol(rhs.ncol),
  colLower(rhs.colLower),
  colUpper(rhs.colUpper),
  rowLower(rhs.rowLower),
  rowUpper(rhs.rowUpper),
  rowRhs(rhs.rowRhs),
  isInteger(rhs.isInteger),
  cstat(rhs.cstat),
  rstat(rhs.rstat),
  solver(rhs.solver),
  xlp(rhs.xlp),
  rowActivity(rhs.rowActivity),
  byRow(rhs.byRow),
  byCol(rhs.byCol),
  f0(rhs.f0),
  f0compl(rhs.f0compl),
  ratiof0compl(rhs.ratiof0compl)
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
  ,
  trackRejection(rhs.trackRejection),
  fracFail(rhs.fracFail),
  dynFail(rhs.dynFail),
  violFail(rhs.violFail),
  suppFail(rhs.suppFail),
  scaleFail(rhs.scaleFail),
  numGeneratedCuts(rhs.numGeneratedCuts)
#endif
{

}

/***************************************************************************/
CglGMI & CglGMI::operator=(const CglGMI& rhs) {
  if(this != &rhs){
    CglCutGenerator::operator=(rhs);
    param = rhs.param;
    nrow = rhs.nrow;
    ncol = rhs.ncol;
    colLower = rhs.colLower;
    colUpper = rhs.colUpper;
    rowLower = rhs.rowLower;
    rowUpper = rhs.rowUpper;
    rowRhs = rhs.rowRhs;
    isInteger = rhs.isInteger;
    cstat = rhs.cstat;
    rstat = rhs.rstat;
    solver = rhs.solver;
    xlp = rhs.xlp;
    rowActivity = rhs.rowActivity;
    byRow = rhs.byRow;
    byCol = rhs.byCol;
    f0 = rhs.f0;
    f0compl = rhs.f0compl;
    ratiof0compl = rhs.ratiof0compl;
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
    trackRejection = rhs.trackRejection;
    fracFail = rhs.fracFail;
    dynFail = rhs.dynFail;
    violFail = rhs.violFail;
    suppFail = rhs.suppFail;
    scaleFail = rhs.scaleFail;
    numGeneratedCuts = rhs.numGeneratedCuts;
#endif						
  }
  return *this;
}

/***************************************************************************/
CglGMI::~CglGMI() {

}

/*********************************************************************/
CglCutGenerator *
CglGMI::clone() const
{
  return new CglGMI(*this);
}

/***************************************************************************/

// Returns (value - floor)
inline double CglGMI::aboveInteger(double value) const {
  return (value - floor(value));
} /* aboveInteger */

/**********************************************************/
void CglGMI::printvecINT(const char *vecstr, const int *x, int n) const {
  int num, fromto, upto;

  num = (n/10) + 1;
  printf("%s :\n", vecstr);
  for (int j = 0; j < num; ++j) {
    fromto = 10*j;
    upto = 10 * (j+1);
    if(n <= upto) upto = n;
    for (int i = fromto; i < upto; ++i)
      printf(" %4d", x[i]);
    printf("\n");
  }
  printf("\n");
} /* printvecINT */

/**********************************************************/
void CglGMI::printvecDBL(const char *vecstr, const double *x, int n) const
{
  int num, fromto, upto;

  num = (n/10) + 1;
  printf("%s :\n", vecstr);
  for (int j = 0; j < num; ++j) {
    fromto = 10*j;
    upto = 10 * (j+1);
    if(n <= upto) upto = n;
    for (int i = fromto; i < upto; ++i)
      printf(" %7.3f", x[i]);
    printf("\n");
  }
  printf("\n");
} /* printvecDBL */

/**********************************************************/
void CglGMI::printvecDBL(const char *vecstr, const double *elem, 
			 const int * index, int nz) const
{
  printf("%s\n", vecstr);
  int written = 0;
  for (int j = 0; j < nz; ++j) {
    written += printf("%d:%.3f ", index[j], elem[j]);
    if (written > 70) {
      printf("\n");
      written = 0;
    }
  }
  if (written > 0) {
    printf("\n");
  }

} /* printvecDBL */

/************************************************************************/
inline bool CglGMI::computeCutFractionality(double varRhs, 
					    double& cutRhs) {
  f0 = aboveInteger(varRhs);
  f0compl = 1 - f0;
  if (f0 < param.getAway() || f0compl < param.getAway())
    return false;
  ratiof0compl = f0/f0compl;
  cutRhs = -f0;
  return true;
} /* computeCutFractionality */

/************************************************************************/
inline double CglGMI::computeCutCoefficient(double rowElem, int index) {
  
  // See Wolsey "Integer Programming" (1998), p. 130, fourth line from top
  // after correcting typo (Proposition 8.8), flipping all signs to get <=.

  if (index < ncol && isInteger[index]) {
    double f = aboveInteger(rowElem);
    if(f > f0) {
      return (-((1-f) * ratiof0compl));
    }
    else {
      return (-f);
    }
  }
  else{
    if(rowElem < 0) {
      return (rowElem*ratiof0compl);
    }
    else {
      return (-rowElem);
    }
  }
  
} /* computeCutCoefficient */

/************************************************************************/
inline void CglGMI::eliminateSlack(double cutElem, int index, double* cut,
				   double& cutRhs, const double *elements, 
				   const CoinBigIndex *rowStart, const int *indices, 
				   const int *rowLength, const double *rhs) {

  // now i is where coefficients on slack variables begin;
  // eliminate the slacks
  int rowpos = index - ncol;
  if(fabs(cutElem) > param.getEPS_ELIM()) {
    if (areEqual(rowLower[rowpos], rowUpper[rowpos], 
		 param.getEPS(), param.getEPS())) {
      // "almost" fixed slack, we'll just skip it
      return;
    }

    CoinBigIndex upto = rowStart[rowpos] + rowLength[rowpos];
    for (CoinBigIndex j = rowStart[rowpos]; j < upto; ++j) {
      cut[indices[j]] -= cutElem * elements[j];      
    }
    cutRhs -= cutElem * rhs[rowpos];
  }

} /* eliminateSlack */

/************************************************************************/
inline void CglGMI::flip(double& rowElem, int index) {  
  if ((index < ncol && cstat[index] == 2) ||
      (index >= ncol && rstat[index-ncol] == 2)) {
      rowElem = -rowElem;
  }
} /* flip */

/************************************************************************/
inline void CglGMI::unflipOrig(double& rowElem, int index, double& rowRhs) {
  if (cstat[index] == 2) {
    // structural variable at upper bound
    rowElem = -rowElem;
    rowRhs += rowElem*colUpper[index];
  }
  else if (cstat[index] == 3) {
    // structural variable at lower bound
    rowRhs += rowElem*colLower[index];
  }
} /* unflipOrig */

/************************************************************************/
inline void CglGMI::unflipSlack(double& rowElem, int index, double& rowRhs, 
				const double* slackVal) {
  if (rstat[index-ncol] == 2) {
    // artificial variable at upper bound
    rowElem = -rowElem;
    rowRhs += rowElem*slackVal[index-ncol];
  }
  else if (rstat[index-ncol] == 3) {
    // artificial variable at lower bound
    rowRhs += rowElem*slackVal[index-ncol];
  }

} /* unflipSlack */

/************************************************************************/
inline void CglGMI::packRow(double* row, double* rowElem, int* rowIndex,
			     int& rowNz) {
  rowNz = 0;
  for (int i = 0; i < ncol; ++i) {
    if (!isZero(fabs(row[i]))) {
      rowElem[rowNz] = row[i];
      rowIndex[rowNz] = i;
      rowNz++;
    }
  }
}

/************************************************************************/
bool CglGMI::cleanCut(double* cutElem, int* cutIndex, int& cutNz,
		       double& cutRhs, const double* xbar) {
  CglGMIParam::CleaningProcedure cleanProc = param.getCLEAN_PROC();
  if (cleanProc == CglGMIParam::CP_CGLLANDP1) {
    if (!checkViolation(cutElem, cutIndex, cutNz, cutRhs, xbar)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad violation\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	violFail++;
      }
#endif
      return false;
    }
    relaxRhs(cutRhs);
    removeSmallCoefficients(cutElem, cutIndex, cutNz, cutRhs);
    if (!checkSupport(cutNz)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: too large support\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	suppFail++;
      }
#endif
      return false;
    }
    if (!checkDynamism(cutElem, cutIndex, cutNz)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad dynamism\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	dynFail++;
      }
#endif
      return false;
    }
    if (!checkViolation(cutElem, cutIndex, cutNz, cutRhs, xbar)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad violation (final check)\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	violFail++;
      }
#endif
      return false;
    }
  } /* end of cleaning procedure CP_CGLLANDP1 */
  else if (cleanProc == CglGMIParam::CP_CGLLANDP2) {
    if (!checkViolation(cutElem, cutIndex, cutNz, cutRhs, xbar)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad violation\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	violFail++;
      }
#endif
      return false;
    }
    relaxRhs(cutRhs);
    if (!checkDynamism(cutElem, cutIndex, cutNz)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad dynamism\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	dynFail++;
      }
#endif
      return false;
    }
    if (!scaleCut(cutElem, cutIndex, cutNz, cutRhs, 1) &&
	param.getENFORCE_SCALING()) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad scaling\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	scaleFail++;
      }
#endif
      return false;
    }
    removeSmallCoefficients(cutElem, cutIndex, cutNz, cutRhs);
    if (!checkSupport(cutNz)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: too large support\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	suppFail++;
      }
#endif
      return false;
    }
    if (!checkViolation(cutElem, cutIndex, cutNz, cutRhs, xbar)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad violation (final check)\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	violFail++;
      }
#endif
      return false;
    }
  } /* end of cleaning procedure CP_CGLLANDP2 */
  else if (cleanProc == CglGMIParam::CP_CGLREDSPLIT) {
    if (!scaleCut(cutElem, cutIndex, cutNz, cutRhs, 3) &&
	param.getENFORCE_SCALING()) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad scaling\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	scaleFail++;
      }
#endif
      return false;
    }
    removeSmallCoefficients(cutElem, cutIndex, cutNz, cutRhs);
    if (!checkDynamism(cutElem, cutIndex, cutNz)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad dynamism\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	dynFail++;
      }
#endif
      return false;
    }
    if (!checkSupport(cutNz)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: too large support\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	suppFail++;
      }
#endif
      return false;
    }
    if (!checkViolation(cutElem, cutIndex, cutNz, cutRhs, xbar)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad violation (final check)\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	violFail++;
      }
#endif
      return false;
    }
    relaxRhs(cutRhs);
  } /* end of cleaning procedure CP_CGLREDSPLIT */
  else if (cleanProc == CglGMIParam::CP_INTEGRAL_CUTS) {
    removeSmallCoefficients(cutElem, cutIndex, cutNz, cutRhs);
    if (!checkSupport(cutNz)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: too large support\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	suppFail++;
      }
#endif
      return false;
    }
    if (!checkDynamism(cutElem, cutIndex, cutNz)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad dynamism\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	dynFail++;
      }
#endif
      return false;
    }
    if (!scaleCut(cutElem, cutIndex, cutNz, cutRhs, 0) &&
	param.getENFORCE_SCALING()) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad scaling\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	scaleFail++;
      }
#endif
      return false;
    }
    if (!checkViolation(cutElem, cutIndex, cutNz, cutRhs, xbar)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad violation (final check)\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	violFail++;
      }
#endif
      return false;
    }
  } /* end of cleaning procedure CP_INTEGRAL_CUTS */
  else if (cleanProc == CglGMIParam::CP_CGLLANDP1_INT) {
    if (!checkViolation(cutElem, cutIndex, cutNz, cutRhs, xbar)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad violation\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	violFail++;
      }
#endif
      return false;
    }
    removeSmallCoefficients(cutElem, cutIndex, cutNz, cutRhs);
    if (!checkSupport(cutNz)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: too large support\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	suppFail++;
      }
#endif
      return false;
    }
    if (!checkDynamism(cutElem, cutIndex, cutNz)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad dynamism\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	dynFail++;
      }
#endif
      return false;
    }
    // scale cut so that it becomes integral, if possible
    if (!scaleCut(cutElem, cutIndex, cutNz, cutRhs, 0)) {
      if (param.getENFORCE_SCALING()){
#if defined GMI_TRACE_CLEAN
	printf("CglGMI::cleanCut(): cut discarded: bad scaling\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
	if (trackRejection) {
	  scaleFail++;
	}
#endif
	return false;
      }
      else {
	// If cannot scale to integral and not enforcing, relax rhs
	// (as per CglLandP cleaning procedure)
	relaxRhs(cutRhs);
      }
    }
    if (!checkViolation(cutElem, cutIndex, cutNz, cutRhs, xbar)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad violation (final check)\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	violFail++;
      }
#endif
      return false;
    }
  } /* end of cleaning procedure CP_CGLLANDP1_INT */
  else if (cleanProc == CglGMIParam::CP_CGLLANDP1_SCALEMAX ||
	   cleanProc == CglGMIParam::CP_CGLLANDP1_SCALERHS) {
    if (!checkViolation(cutElem, cutIndex, cutNz, cutRhs, xbar)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad violation\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	violFail++;
      }
#endif
      return false;
    }
    if (// Try to scale cut, but do not discard if cannot scale
	((cleanProc == CglGMIParam::CP_CGLLANDP1_SCALEMAX &&
	  !scaleCut(cutElem, cutIndex, cutNz, cutRhs, 1)) ||
	 (cleanProc == CglGMIParam::CP_CGLLANDP1_SCALERHS &&
	  !scaleCut(cutElem, cutIndex, cutNz, cutRhs, 2))) &&
	param.getENFORCE_SCALING()) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad scaling\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	scaleFail++;
      }
#endif
      return false;
    }
    relaxRhs(cutRhs);
    removeSmallCoefficients(cutElem, cutIndex, cutNz, cutRhs);
    if (!checkSupport(cutNz)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: too large support\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	suppFail++;
      }
#endif
      return false;
    }
    if (!checkDynamism(cutElem, cutIndex, cutNz)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad dynamism\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	dynFail++;
      }
#endif
      return false;
    }
    if (!checkViolation(cutElem, cutIndex, cutNz, cutRhs, xbar)) {
#if defined GMI_TRACE_CLEAN
      printf("CglGMI::cleanCut(): cut discarded: bad violation (final check)\n");
#endif
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {
	violFail++;
      }
#endif
      return false;
    }
  } /* end of cleaning procedures CP_CGLLANDP1_SCALEMAX and CG_CGLLANDP1_SCALERHS */
  return true;
}

/************************************************************************/
bool CglGMI::checkViolation(const double* cutElem, const int* cutIndex,
			     int cutNz, double cutrhs, const double* xbar) {
  double lhs = 0.0;
  for (int i = 0; i < cutNz; ++i) {
    lhs += cutElem[i]*xbar[cutIndex[i]];
  }
  double violation = lhs - cutrhs;
  if (fabs(cutrhs) > 1) {
    violation /= fabs(cutrhs);
  }
  if (violation >= param.getMINVIOL()) {
    return true;
  }
  else{
#if defined GMI_TRACE_CLEAN
    printf("Cut lhs %g, rhs %g, violation %g; cut discarded\n", lhs, cutrhs, violation);
#endif
    return false;
  }
} /* checkViolation */

/************************************************************************/
bool CglGMI::checkDynamism(const double* cutElem, const int* cutIndex,
			   int cutNz) {
  double min = param.getINFINIT();
  double max = 0.0;
  double val = 0.0;
  for (int i = 0; i < cutNz; ++i) {
    if (!isZero(cutElem[i])) {
      val = fabs(cutElem[i]);
      min = std::min(min, val);
      max = std::max(max, val);
    }
  }
  if (max > min*param.getMAXDYN()) {
#if defined GMI_TRACE_CLEAN
    printf("Max elem %g, min elem %g, dyn %g; cut discarded\n", max, min, max/min);
#endif
    return false;
  }
  else{
    return true;
  }
  
} /* checkDynamism */

/************************************************************************/
bool CglGMI::checkSupport(int cutNz) {
  if (cutNz > param.getMAX_SUPPORT_ABS() + param.getMAX_SUPPORT_REL()*ncol) {
#if defined GMI_TRACE_CLEAN
    printf("Support %d; cut discarded\n", cutNz);
#endif
    return false;
  }
  else{
    return true;
  }
}

/************************************************************************/
bool CglGMI::removeSmallCoefficients(double* cutElem, int* cutIndex, 
				     int& cutNz, double& cutRhs) {
  double value, absval;
  int currPos = 0;
  int col;
  for (int i = 0; i < cutNz; ++i) {
    col = cutIndex[i];
    value = cutElem[i];
    absval = fabs(value);
    if (!isZero(absval) && absval <= param.getEPS_COEFF()) {
      // small coefficient: remove and adjust rhs if possible
      if ((value > 0.0) && (colLower[col] > -param.getINFINIT())) {
        cutRhs -= value * colLower[col];
      } 
      else if ((value < 0.0) && (colUpper[col] < param.getINFINIT())) {
        cutRhs -= value * colUpper[col];      
      }
    }
    else if (absval > param.getEPS_COEFF()) {
      if (currPos < i) {
	cutElem[currPos] = cutElem[i];
	cutIndex[currPos] = cutIndex[i];
      }
      currPos++;
    }
  }
  cutNz = currPos;
  return true;
}

/************************************************************************/
void CglGMI::relaxRhs(double& rhs) {
  if(param.getEPS_RELAX_REL() > 0.0) {
    rhs += fabs(rhs) * param.getEPS_RELAX_REL() + param.getEPS_RELAX_ABS();
  }
  else{
    rhs += param.getEPS_RELAX_ABS();
  }
}

/************************************************************************/
bool CglGMI::scaleCut(double* cutElem, int* cutIndex, int cutNz,
		       double& cutRhs, int scalingType) {
  /// scalingType possible values:
  /// 0 : scale to obtain integral cut
  /// 1 : scale to obtain largest coefficient equal to 1
  /// 2 : scale to obtain rhs equal to 1
  /// 3 : scale based on norm, to obtain cut norm equal to ncol
  /// Returns true if scaling is successful.
  if (scalingType == 0) {
    return scaleCutIntegral(cutElem, cutIndex, cutNz, cutRhs);
  }
  else if (scalingType == 1) {
    double max = fabs(cutRhs);
    for (int i = 0; i < cutNz; ++i) {
      if (!isZero(cutElem[i])) {
	max = std::max(max, fabs(cutElem[i]));
      }
    }
    if (max < param.getEPS() || max > param.getMAXDYN()) {
#if defined GMI_TRACE_CLEAN
      printf("Scale %g; %g %g cut discarded\n", max, param.getEPS(), 1/param.getMAXDYN());
#endif
      return false;
    }
    else{
      for (int i = 0; i < cutNz; ++i) {
	cutElem[i] /= max;
      }
      cutRhs /= max;
      return true;
    }
  }
  else if (scalingType == 2) {
    double max = fabs(cutRhs);
    if (max < param.getEPS() || max > param.getMAXDYN()) {
#if defined GMI_TRACE_CLEAN
      printf("Scale %g; %g %g cut discarded\n", max, param.getEPS(), 1/param.getMAXDYN());
#endif
      return false;
    }
    else{
      for (int i = 0; i < cutNz; ++i) {
	cutElem[i] /= max;
      }
      cutRhs /= max;
      return true;
    }
  }
  else if (scalingType == 3) {
    int support = 0;
    double norm = 0.0;
    for (int i = 0; i < cutNz; ++i) {
      if (!isZero(fabs(cutElem[i]))) {
	support++;
	norm += cutElem[i]*cutElem[i];
      }
    }
    double scale = sqrt(norm / support);
    if ((scale < 0.02) || (scale > 100)) {
#if defined GMI_TRACE_CLEAN
      printf("Scale %g; cut discarded\n", scale);
#endif
      return false;
    }
    else{
      for (int i = 0; i < cutNz; ++i) {
	cutElem[i] /= scale;
      }
      cutRhs /= scale;
      return true;
    }
  }
  return false;
} /* scaleCut */

/************************************************************************/
bool CglGMI::scaleCutIntegral(double* cutElem, int* cutIndex, int cutNz,
			      double& cutRhs) {
  int64_t gcd, lcm;
  double maxdelta = param.getEPS(); 
  double maxscale = 1000; 
  int64_t maxdnom = 1000;
  int64_t numerator = 0, denominator = 0;
  // Initialize gcd and lcm
  CoinRational r = CoinRational(cutRhs, maxdelta, maxdnom);
  if (r.getNumerator() != 0){
     gcd = llabs(r.getNumerator());
     lcm = r.getDenominator();
  }
  else{
#if defined GMI_TRACE_CLEAN
      printf("Cannot compute rational number, scaling procedure aborted\n");
#endif
    return false;
  }
  for (int i = 0; i < cutNz; ++i) {
    if (solver->isContinuous(cutIndex[i]) && !param.getINTEGRAL_SCALE_CONT()) {
      continue;
    }
    CoinRational r = CoinRational(cutElem[i], maxdelta, maxdnom);
    if (r.getNumerator() != 0){
       gcd = computeGcd(gcd, r.getNumerator());
       lcm *= r.getDenominator()/(computeGcd(lcm,r.getDenominator()));
    }
    else{
#if defined GMI_TRACE_CLEAN
      printf("Cannot compute rational number, scaling procedure aborted\n");
#endif
      return false;
    } 
  }
  double scale = ((double)lcm)/((double)gcd);
  if (fabs(scale) > maxscale) {
#if defined GMI_TRACE_CLEAN
      printf("Scaling factor too large, scaling procedure aborted\n");
#endif
      return false;
  }
  // Looks like we have a good scaling factor; scale and return;
  for (int i = 0; i < cutNz; ++i) {
    cutElem[i] *= scale;
  }
  cutRhs *= scale;
  return true;
} /* scaleCutIntegral */

/************************************************************************/
int64_t CglGMI::computeGcd(int64_t a, int64_t b) {
  // This is the standard Euclidean algorithm for gcd
  int64_t remainder = 1;
  // Make sure a<=b (will always remain so)
  if (a > b) {
    // Swap a and b
    int64_t temp = a;
    a = b;
    b = temp;
  }
  // If zero then gcd is nonzero
  if (!a) {
    if (b) {
      return b;
    } 
    else {
      printf("### WARNING: CglGMI::computeGcd() given two zeroes!\n");
      exit(1);
    }
  }
  while (remainder) {
    remainder = b % a;
    b = a;
    a = remainder;
  }
  return b;
} /* computeGcd */


/************************************************************************/
void CglGMI::generateCuts(const OsiSolverInterface &si, OsiCuts & cs,
			  const CglTreeInfo )
{
  solver = const_cast<OsiSolverInterface *>(&si);
  if (solver == NULL) {
    printf("### WARNING: CglGMI::generateCuts(): no solver available.\n");
    return;    
  }  

  if (!solver->optimalBasisIsAvailable()) {
    printf("### WARNING: CglGMI::generateCuts(): no optimal basis available.\n");
    return;
  }

#if defined OSI_TABLEAU
  if (!solver->canDoSimplexInterface()) {
    printf("### WARNING: CglGMI::generateCuts(): solver does not provide simplex tableau.\n");
    printf("### WARNING: CglGMI::generateCuts(): recompile without OSI_TABLEAU.\n");
    return;
  }
#endif
  

  // Get basic problem information from solver
  ncol = solver->getNumCols(); 
  nrow = solver->getNumRows(); 
  colLower = solver->getColLower();
  colUpper = solver->getColUpper();
  rowLower = solver->getRowLower();
  rowUpper = solver->getRowUpper();
  rowRhs = solver->getRightHandSide();

  xlp = solver->getColSolution();
  rowActivity = solver->getRowActivity();
  byRow = solver->getMatrixByRow();
  byCol = solver->getMatrixByCol();
  
  generateCuts(cs);

} /* generateCuts */

/************************************************************************/
void CglGMI::generateCuts(OsiCuts &cs)
{
  isInteger = new bool[ncol]; 
  
  computeIsInteger();

  cstat = new int[ncol];
  rstat = new int[nrow];


  solver->getBasisStatus(cstat, rstat);   // 0: free  1: basic  
                                          // 2: upper 3: lower


#if defined GMI_TRACETAB
  printvecINT("cstat", cstat, ncol);
  printvecINT("rstat", rstat, nrow);
#endif

  // list of basic integer fractional variables
  int *listFracBasic = new int[nrow];
  int numFracBasic = 0;
  for (int i = 0; i < ncol; ++i) {
    // j is the variable which is basic in row i
    if ((cstat[i] == 1) && (isInteger[i])) {
      if (std::min(aboveInteger(xlp[i]),
		  1-aboveInteger(xlp[i])) > param.getAway()) {
	listFracBasic[numFracBasic] = i;
	numFracBasic++;
      }
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      else if (trackRejection) {	
	// Say that we tried to generate a cut, but it was discarded
	// because of small fractionality
	if (!isIntegerValue(xlp[i])) {
	  fracFail++;
	  numGeneratedCuts++;
	}
      }
#endif
    }
  }

#if defined GMI_TRACE
  printf("CglGMI::generateCuts() : %d fractional rows\n", numFracBasic);
#endif
  
  if (numFracBasic == 0) {
    delete[] listFracBasic;
    delete[] cstat;
    delete[] rstat;
    delete[] isInteger;
    return;
  }

  // there are rows with basic integer fractional variables, so we can
  // generate cuts

  // Basis index for columns and rows; each element is -1 if corresponding
  // variable is nonbasic, and contains the basis index if basic.
  // The basis index is the row in which the variable is basic.
  int* colBasisIndex = new int[ncol];
  int* rowBasisIndex = new int[nrow];

#if defined OSI_TABLEAU
  memset(colBasisIndex, -1, ncol*sizeof(int));
  memset(rowBasisIndex, -1, nrow*sizeof(int));
  solver->enableFactorization();
  int* basicVars = new int[nrow];
  solver->getBasics(basicVars);
  for (int i = 0; i < nrow; ++i) {
    if (basicVars[i] < ncol) {
      colBasisIndex[basicVars[i]] = i;
    }
    else {
      rowBasisIndex[basicVars[i] - ncol] = i;
    }
  }
#else
  CoinFactorization factorization;
  if (factorize(factorization, colBasisIndex, rowBasisIndex)) {
    printf("### WARNING: CglGMI::generateCuts(): error during factorization!\n");
    return;
  }
#endif


  // cut in sparse form
  double* cutElem = new double[ncol];
  int* cutIndex = new int[ncol];
  int cutNz = 0;
  double cutRhs;

  // cut in dense form
  double* cut = new double[ncol];

  double *slackVal = new double[nrow];

  for (int i = 0; i < nrow; ++i) {
    slackVal[i] = rowRhs[i] - rowActivity[i];
  }

#if defined OSI_TABLEAU
  // Column part and row part of a row of the simplex tableau
  double* tableauColPart = new double[ncol];
  double* tableauRowPart = new double[nrow];
#else
  // Need some more data for simplex tableau computation
  const int * row = byCol->getIndices();
  const CoinBigIndex * columnStart = byCol->getVectorStarts();
  const int * columnLength = byCol->getVectorLengths(); 
  const double * columnElements = byCol->getElements();

  // Create work arrays for factorization
  // two vectors for updating: the first one is needed to do the computations
  // but we do not use it, the second one contains a row of the basis inverse
  CoinIndexedVector work;
  CoinIndexedVector array;
  // Make sure they large enough
  work.reserve(nrow);
  array.reserve(nrow);
  int * arrayRows = array.getIndices();
  double * arrayElements = array.denseVector();
  // End of code to create work arrays
  double one = 1.0;
#endif

  // Matrix elements by row for slack substitution
  const double *elements = byRow->getElements();
  const CoinBigIndex *rowStart = byRow->getVectorStarts();
  const int *indices = byRow->getIndices();
  const int *rowLength = byRow->getVectorLengths(); 

  // Indices of basic and slack variables, and cut elements
  int iBasic, slackIndex;
  double cutCoeff;
  double rowElem;
  // Now generate the cuts: obtain a row of the simplex tableau
  // where an integer variable is basic and fractional, and compute the cut
  for (int i = 0; i < numFracBasic; ++i) {
    if (!computeCutFractionality(xlp[listFracBasic[i]], cutRhs)) {
      // cut is discarded because of the small fractionalities involved
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
      if (trackRejection) {	
	// Say that we tried to generate a cut, but it was discarded
	// because of small fractionality
	fracFail++;
	numGeneratedCuts++;
      }
#endif
      continue;
    }

    // the variable listFracBasic[i] is basic in row iBasic
    iBasic = colBasisIndex[listFracBasic[i]];

#if defined GMI_TRACE
    printf("Row %d with var %d basic, f0 = %f\n", i, listFracBasic[i], f0);
#endif

#if defined OSI_TABLEAU
    solver->getBInvARow(iBasic, tableauColPart, tableauRowPart);
#else
    array.clear();
    array.setVector(1, &iBasic, &one);

    factorization.updateColumnTranspose (&work, &array);

    int numberInArray=array.getNumElements();
#endif

    // reset the cut
    memset(cut, 0, ncol*sizeof(double));

    // columns
    for (int j = 0; j < ncol; ++j) {
      if ((colBasisIndex[j] >= 0) || 
	  (areEqual(colLower[j], colUpper[j], 
		    param.getEPS(), param.getEPS()))) {
	// Basic or fixed variable -- skip
	continue;
      }
#ifdef OSI_TABLEAU
      rowElem = tableauColPart[j];
#else
      rowElem = 0.0;
      // add in row of tableau
      for (CoinBigIndex h = columnStart[j]; h < columnStart[j]+columnLength[j]; ++h) {
	rowElem += columnElements[h]*arrayElements[row[h]];
      }
#endif
      if (!isZero(fabs(rowElem))) {
	// compute cut coefficient
	flip(rowElem, j);
	cutCoeff = computeCutCoefficient(rowElem, j);
	if (isZero(cutCoeff)) {
	  continue;
	}
	unflipOrig(cutCoeff, j, cutRhs);
	cut[j] = cutCoeff;
#if defined GMI_TRACE
	printf("var %d, row %f, cut %f\n", j, rowElem, cutCoeff);
#endif
      }
    }

    // now do slacks part
#if defined OSI_TABLEAU
    for (int j = 0 ; j < nrow; ++j) {
      // index of the row corresponding to the slack variable
      slackIndex = j;
      if (rowBasisIndex[j] >= 0) {
	// Basic variable -- skip it
	continue;
      }
      rowElem = tableauRowPart[j];
#else
    for (int j = 0 ; j < numberInArray ; ++j) {
      // index of the row corresponding to the slack variable
      slackIndex = arrayRows[j];
      rowElem = arrayElements[slackIndex];
#endif
      if (!isZero(fabs(rowElem))) {
	slackIndex += ncol;
	// compute cut coefficient
	flip(rowElem, slackIndex);
	cutCoeff = computeCutCoefficient(rowElem, slackIndex);
	if (isZero(fabs(cutCoeff))) {
	  continue;
	}
	unflipSlack(cutCoeff, slackIndex, cutRhs, slackVal);
	eliminateSlack(cutCoeff, slackIndex, cut, cutRhs,
		       elements, rowStart, indices, rowLength, rowRhs);
#if defined GMI_TRACE
	printf("var %d, row %f, cut %f\n", slackIndex, rowElem, cutCoeff);
#endif
      }
    }

    packRow(cut, cutElem, cutIndex, cutNz);
    if (cutNz == 0)
      continue;

#if defined GMI_TRACE
    printvecDBL("final cut:", cutElem, cutIndex, cutNz);
    printf("cutRhs: %f\n", cutRhs);
#endif
    
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
    if (trackRejection) {
      numGeneratedCuts++;
    }
#endif
    if (cleanCut(cutElem, cutIndex, cutNz, cutRhs, xlp) && cutNz > 0) {
      // relax if integer values can be large
      for (int i=0;i<cutNz;i++) {
	int j = cutIndex[i];
	if (isInteger[j]) {
	  double difference = colUpper[j]-colLower[j];
	  if (difference > 10.0) 
	    cutRhs += 1.0e-8*difference;
	}
      }
      OsiRowCut rc;
      rc.setRow(cutNz, cutIndex, cutElem);
      rc.setLb(-param.getINFINIT());
      rc.setUb(cutRhs);
      if (!param.getCHECK_DUPLICATES()) {
	cs.insertIfNotDuplicateAndClean(rc,21);
      }
      else{
	cs.insertIfNotDuplicateAndClean(rc, 22, CoinAbsFltEq(param.getEPS_COEFF()));
      }
    }

  }

#if defined GMI_TRACE
  printf("CglGMI::generateCuts() : number of cuts : %d\n", cs.sizeRowCuts());
#endif

#if defined OSI_TABLEAU
  solver->disableFactorization();
  delete[] basicVars;
  delete[] tableauColPart;
  delete[] tableauRowPart;
#endif

  delete[] colBasisIndex;
  delete[] rowBasisIndex;
  delete[] cut;
  delete[] slackVal;
  delete[] cutElem;
  delete[] cutIndex;
  delete[] listFracBasic;
  delete[] cstat;
  delete[] rstat;
  delete[] isInteger;

} /* generateCuts */

/***********************************************************************/
void CglGMI::setParam(const CglGMIParam &source) {
  param = source;
} /* setParam */

/***********************************************************************/
void CglGMI::computeIsInteger() {
  for (int i = 0; i < ncol; ++i) {
    if(solver->isInteger(i)) {
      isInteger[i] = true;
    }
    else {
      if((areEqual(colLower[i], colUpper[i], 
		   param.getEPS(), param.getEPS()))
	 && (isIntegerValue(colUpper[i]))) {	
	// continuous variable fixed to an integer value
	isInteger[i] = true;
      }
      else {
	isInteger[i] = false;
      }
    }
  }    
} /* computeIsInteger */

/***********************************************************************/
void CglGMI::printOptTab(OsiSolverInterface *lclSolver) const
{
  int *cstat = new int[ncol];
  int *rstat = new int[nrow];

  lclSolver->enableFactorization();
  lclSolver->getBasisStatus(cstat, rstat);   // 0: free  1: basic  
  // 2: upper 3: lower

  int *basisIndex = new int[nrow]; // basisIndex[i] = 
                                    //        index of pivot var in row i
                                    //        (slack if number >= ncol) 
  lclSolver->getBasics(basisIndex);

  double *z = new double[ncol];  // workspace to get row of the tableau
  double *slack = new double[nrow];  // workspace to get row of the tableau
  double *slackVal = new double[nrow];

  for (int i = 0; i < nrow; i++) {
    slackVal[i] = rowRhs[i] - rowActivity[i];
  }

  const double *rc = lclSolver->getReducedCost();
  const double *dual = lclSolver->getRowPrice();
  const double *solution = lclSolver->getColSolution();

  printvecINT("cstat", cstat, ncol);
  printvecINT("rstat", rstat, nrow);
  printvecINT("basisIndex", basisIndex, nrow);

  printvecDBL("solution", solution, ncol);
  printvecDBL("slackVal", slackVal, nrow);
  printvecDBL("reduced_costs", rc, ncol);
  printvecDBL("dual solution", dual, nrow);

  printf("Optimal Tableau:\n");

  for (int i = 0; i < nrow; i++) {
    lclSolver->getBInvARow(i, z, slack);
    for (int ii = 0; ii < ncol; ++ii) {
      printf("%5.2f ", z[ii]);
    }
    printf(" | ");
    for (int ii = 0; ii < nrow; ++ii) {
      printf("%5.2f ", slack[ii]);
    }
    printf(" | ");
    if(basisIndex[i] < ncol) {
      printf("%5.2f ", solution[basisIndex[i]]);
    }
    else {
      printf("%5.2f ", slackVal[basisIndex[i]-ncol]);
    }
    printf("\n");
  }
  for (int ii = 0; ii < 7*(ncol+nrow+1); ++ii) {
    printf("-");
  }
  printf("\n");

  for (int ii = 0; ii < ncol; ++ii) {
    printf("%5.2f ", rc[ii]);    
  }
  printf(" | ");
  for (int ii = 0; ii < nrow; ++ii) {
    printf("%5.2f ", -dual[ii]);
  }
  printf(" | ");
  printf("%5.2f\n", -lclSolver->getObjValue());
  lclSolver->disableFactorization();

  delete[] cstat;
  delete[] rstat;
  delete[] basisIndex;
  delete[] slack;
  delete[] z;
  delete[] slackVal;
} /* printOptTab */


/*********************************************************************/
// Create C++ lines to get to current state
std::string
CglGMI::generateCpp(FILE * fp) 
{
  CglGMI other;
  fprintf(fp,"0#include \"CglGMI.hpp\"\n");
  fprintf(fp,"3  CglGMI GMI;\n");
  if (param.getMAX_SUPPORT()!=other.param.getMAX_SUPPORT())
    fprintf(fp,"3  GMI.setLimit(%d);\n",param.getMAX_SUPPORT());
  else
    fprintf(fp,"4  GMI.setLimit(%d);\n",param.getMAX_SUPPORT());
  if (param.getAway()!=other.param.getAway())
    fprintf(fp,"3  GMI.setAway(%g);\n",param.getAway());
  else
    fprintf(fp,"4  GMI.setAway(%g);\n",param.getAway());
  if (param.getEPS()!=other.param.getEPS())
    fprintf(fp,"3  GMI.setEPS(%g);\n",param.getEPS());
  else
    fprintf(fp,"4  GMI.setEPS(%g);\n",param.getEPS());
  if (param.getEPS_COEFF()!=other.param.getEPS_COEFF())
    fprintf(fp,"3  GMI.setEPS_COEFF(%g);\n",param.getEPS_COEFF());
  else
    fprintf(fp,"4  GMI.set.EPS_COEFF(%g);\n",param.getEPS_COEFF());
  if (param.getEPS_RELAX_ABS()!=other.param.getEPS_RELAX_ABS())
    fprintf(fp,"3  GMI.set.EPS_RELAX(%g);\n",param.getEPS_RELAX_ABS());
  else
    fprintf(fp,"4  GMI.set.EPS_RELAX(%g);\n",param.getEPS_RELAX_ABS());
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  GMI.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  GMI.setAggressiveness(%d);\n",getAggressiveness());
  return "GMI";
}

/*********************************************************************/
int
CglGMI::factorize(CoinFactorization & factorization,
		  int* colBasisIndex, int* rowBasisIndex) {
  // Start of code to create a factorization from warm start ====
  // Taken (with small modifications) from CglGomory
  int status=-100;
  for (int i = 0; i < nrow; ++i) {
    if (rstat[i] == 1) {
      rowBasisIndex[i]=1;
    } else {
      rowBasisIndex[i]=-1;
    }
  }
  for (int i = 0; i < ncol; ++i) {
    if (cstat[i] == 1) {
      colBasisIndex[i]=1;
    } else {
      colBasisIndex[i]=-1;
    }
  }
  // returns 0 if okay, -1 singular, -2 too many in basis, -99 memory */
  while (status<-98) {
    status=factorization.factorize(*byCol, rowBasisIndex, colBasisIndex);
    if (status==-99) factorization.areaFactor(factorization.areaFactor()*2.0);
  } 
  if (status) {
    return -1;
  }
#if defined GMI_TRACE
  double condition = 0.0;
  const CoinFactorizationDouble * pivotRegion = factorization.pivotRegion();
  for (int i = 0; i < nrow; ++i) {
    condition += log(fabs(pivotRegion[i]));
  }  
  printf("CglGMI::factorize(): condition number recomputed as sum of log: %g\n", (condition));
#endif
  return 0;
}

/*********************************************************************/
void CglGMI::setTrackRejection(bool value) {
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
  trackRejection = value;
  if (trackRejection) {
    // reset data members
    resetRejectionCounters();
  }
#endif
}

/*********************************************************************/
bool CglGMI::getTrackRejection() {
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
  return trackRejection;
#else
  return false;
#endif
}

/*********************************************************************/
void CglGMI::resetRejectionCounters() {
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
  fracFail = 0;
  dynFail = 0;
  violFail = 0;
  suppFail = 0;
  scaleFail = 0;
  numGeneratedCuts = 0;
#endif
}

/*********************************************************************/
int CglGMI::getNumberRejectedCuts(RejectionType reason) {
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
  switch (reason) {
  case failureFractionality:
    return fracFail;
  case failureDynamism:
    return dynFail;
  case failureViolation:
    return violFail;
  case failureSupport:
    return suppFail;
  case failureScale:
    return scaleFail;
  }
  return 0;
#else
  return 0;
#endif
}

/*********************************************************************/
int CglGMI::getNumberGeneratedCuts() {
#if defined TRACK_REJECT || defined TRACK_REJECT_SIMPLE
  return numGeneratedCuts;
#else
  return 0;
#endif
}
