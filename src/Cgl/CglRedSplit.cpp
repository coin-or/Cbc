// Last edit: 4/20/07
//
// Name:     CglRedSplit.cpp
// Author:   Francois Margot                                                  
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     2/6/05
//
//---------------------------------------------------------------------------
// Copyright (C) 2005, Francois Margot and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>

//#define RS_TRACE
//#define RS_TRACEALL
//#define RS_TRACETAB

#include "OsiSolverInterface.hpp"

#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CoinFactorization.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CglRedSplit.hpp"
#include "CoinFinite.hpp"

//-------------------------------------------------------------------
// Generate Gomory Reduce-and-Split cuts
//------------------------------------------------------------------- 

/***************************************************************************/
// Returns (value - floor) but allowing for small errors
inline double CglRedSplit::rs_above_integer(double value) 
{
  double value2=floor(value);
  double value3=floor(value+0.5);
  if (fabs(value3-value)< param.getEPS() * (fabs(value3)+1.0))
    return 0.0;
  return value-value2;
} /* rs_above_integer */

/**********************************************************/
void rs_allocmatINT(int ***v, const int m, const int n)
{
  int i;

  *v = reinterpret_cast<int **> (calloc (m, sizeof(int *)));
  if (*v == NULL) {
    printf("###ERROR: INTEGER matrix allocation failed\n");
    exit(1);
  }

  for(i=0; i<m; i++) {
    (*v)[i] = reinterpret_cast<int *> (calloc (n, sizeof(int)));
    if ((*v)[i] == NULL) {
      printf("###ERROR: INTEGER matrix allocation failed\n");
      exit(1);
    }
  }
} /* rs_allocmatINT */

/**********************************************************/
void rs_deallocmatINT(int ***v, const int m, const int /*n*/)
{
  int i;

  for(i=0; i<m; i++) {
    free(reinterpret_cast<void *> ((*v)[i]));
  }
  free(reinterpret_cast<void *> (*v));
} /* rs_deallocmatINT */

/**********************************************************/
void rs_allocmatDBL(double ***v, const int m, const int n)
{
  int i;

  *v = reinterpret_cast<double **> (calloc (m, sizeof(double *)));
  if (*v == NULL) {
    printf("###ERROR: DOUBLE matrix allocation failed\n");
    exit(1);
  }

  for(i=0; i<m; i++) {
    (*v)[i] = reinterpret_cast<double *> (calloc (n, sizeof(double)));
    if ((*v)[i] == NULL) {
      printf("###ERROR: DOUBLE matrix allocation failed\n");
      exit(1);
    }
  }
} /* rs_allocmatDBL */

/**********************************************************/
void rs_deallocmatDBL(double ***v, const int m, const int /*n*/)
{
  int i;

  for(i=0; i<m; i++) {
    free(reinterpret_cast<void *> ((*v)[i]));
  }
  free(reinterpret_cast<void *> (*v));
} /* rs_deallocmatDBL */

/**********************************************************/
void rs_printvecINT(const char *vecstr, const int *x, const int n)
{
  int num, fromto, upto, j, i;

  num = (n/10) + 1;
  printf("%s :\n", vecstr);
  for(j=0; j<num; j++){
    fromto = 10*j;
    upto = 10 * (j+1);
    if(n <= upto) upto = n;
    for(i=fromto; i<upto; i++)
      printf(" %4d", x[i]);
    printf("\n");
  }
  printf("\n");
} /* rs_printvecINT */

/**********************************************************/
void rs_printvecDBL(char const *vecstr, 
		    const double *x, const int n)
{
  int num, fromto, upto, j, i;

  num = (n/10) + 1;
  printf("%s :\n", vecstr);
  for(j=0; j<num; j++){
    fromto = 10*j;
    upto = 10 * (j+1);
    if(n <= upto) upto = n;
    for(i=fromto; i<upto; i++)
      printf(" %7.3f", x[i]);
    printf("\n");
  }
  printf("\n");
} /* rs_printvecDBL */

/**********************************************************/
void rs_printmatINT(char const *vecstr, const int **x, 
		    const int m, const int n)
{
  int i, j;

  printf("%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf(" %4d", x[i][j]);
    }
      printf("\n");
  }
  printf("\n");
} /* rs_printmatINT */

/**********************************************************/
void rs_printmatINT(char const *vecstr, int **x, const int m, const int n)
{
  int i, j;

  printf("%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf(" %4d", x[i][j]);
    }
      printf("\n");
  }
  printf("\n");
} /* rs_printmatINT */

/**********************************************************/
void rs_printmatDBL(char const *vecstr, double **x, const int m, const int n)
{
  int i, j;

  printf("%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf(" %7.3f", x[i][j]);
    }
      printf("\n");
  }
  printf("\n");
} /* rs_printmatDBL */

/***************************************************************************/
double rs_dotProd(const double *u, const double *v, const int dim) {

  int i;
  double result = 0;
  for(i=0; i<dim; i++) {
    result += u[i] * v[i];
  }
  return(result);
} /* rs_dotProd */

/***************************************************************************/
double rs_dotProd(const int *u, const double *v, const int dim) {

  int i;
  double result = 0;
  for(i=0; i<dim; i++) {
    result += u[i] * v[i];
  }
  return(result);
} /* rs_dotProd */

/***************************************************************************/
double rs_genalea (int *x0)
{
  int m = 2147483647;
  int a = 16807 ;
  int b = 127773 ;
  int c = 2836 ;
  int x1, k;

  k = static_cast<int> ((*x0)/b) ;
  x1 = a*(*x0 - k*b) - k*c ;
  if(x1 < 0) x1 = x1 + m;
  *x0 = x1;

  return(static_cast<double>(x1)/static_cast<double>(m));
} /* rs_genalea */

/***********************************************************************/
void CglRedSplit::update_pi_mat(int r1, int r2, int step) {

  int j;
  for(j=0; j<mTab; j++) {
    pi_mat[r1][j] = pi_mat[r1][j] - step * pi_mat[r2][j];
  }
} /* update_pi_mat */

/***********************************************************************/
void CglRedSplit::update_redTab(int r1, int r2, int step) {

  int j;
  for(j=0; j<nTab; j++) {
    contNonBasicTab[r1][j] = 
                  contNonBasicTab[r1][j] - step * contNonBasicTab[r2][j];
  }
} /* update_redTab */

/***********************************************************************/
void CglRedSplit::find_step(int r1, int r2, int *step, 
			    double *reduc, double *norm) {

   double btb_val = rs_dotProd(contNonBasicTab[r1], contNonBasicTab[r2], nTab);
   double opt_step = btb_val/norm[r2];

   double f_step= floor(opt_step);
   double c_step = f_step + 1;

   double val_f = norm[r1] + f_step * f_step * norm[r2] - 2 * btb_val * f_step;
   double val_c = norm[r1] + c_step * c_step * norm[r2] - 2 * btb_val * c_step;

   if(val_f <= val_c ) {
     (*step) = static_cast<int>(f_step);
     (*reduc) = norm[r1] - val_f;
   }
   else {
     (*step) = static_cast<int>(c_step);
     (*reduc) = norm[r1] - val_c;
   }
} /* find_step */

/***************************************************************************/
int CglRedSplit::test_pair(int r1, int r2, double *norm) {

   int step;
   double reduc;
   
   find_step(r1, r2, &step, &reduc, norm); 

   if(reduc/norm[r1] >= param.getMinReduc()) {
     update_pi_mat(r1, r2, step);
     update_redTab(r1, r2, step);
     norm[r1] = rs_dotProd(contNonBasicTab[r1], contNonBasicTab[r1], nTab);

#ifdef RS_TRACEALL
     printf("Use %d and %d for reduction (step: %d)\n", r1, r2, step);
#endif

     return(1);
   }

   return(0);
} /* test_pair */

/***************************************************************************/
void CglRedSplit::reduce_contNonBasicTab() {

  int i, j;

  double *norm = new double[mTab];
  for(i=0; i<mTab; i++) {
    norm[i] = rs_dotProd(contNonBasicTab[i], contNonBasicTab[i], nTab);
  }

#ifdef RS_TRACE
  double sum_norms = 0;
  for(i=0; i<mTab; i++) {
    sum_norms += norm[i];
  }

  printf("CglRedSplit::reduce_contNonBasicTab():Initial sum of  norms: %f\n", 
	 sum_norms);
#endif

  int iter = 0, done = 0;
  int *changed = new int[mTab]; // changed[i]: last iter where row i updated
  int **checked; // checked[i][j]: last iter where pair (i, j) checked
  
  rs_allocmatINT(&checked, mTab, mTab);
  for(i=0; i<mTab; i++) {
    changed[i] = 0;
    for(j=0; j<mTab; j++) {
      checked[i][j] = -1;
    }
    checked[i][i] = 0;
  }

  while(!done) {
    done = 1;

#ifdef RS_TRACEALL
    rs_printmatDBL("contNonBasicTab", contNonBasicTab, mTab, nTab);
    rs_printmatINT("checked", checked, mTab, mTab);
    rs_printvecINT("changed", changed, mTab);
    rs_printvecDBL("norm", norm, mTab);
    rs_printmatINT("pi_mat", pi_mat, mTab, mTab);
#endif

    for(i=0; i<mTab; i++) {
      if(norm[i] > param.getNormIsZero()) {
	for(j=i+1; j<mTab; j++) {
	  if(norm[j] > param.getNormIsZero()) {
	    if((checked[i][j] < changed[i]) || (checked[i][j] < changed[j])) {
	      if(test_pair(i, j, norm)) {
		changed[i] = iter+1;
		done = 0;
	      }
	      checked[i][j] = iter;

	      if((checked[j][i] < changed[i]) || 
		 (checked[j][i] < changed[j])) {
		if(test_pair(j, i, norm)) {
		  changed[j] = iter+1;
		  done = 0;
		}
		checked[j][i] = iter;
	      }
	    }
	  }
	}
      }
    }
    iter++;
  }

#ifdef RS_TRACEALL
  rs_printmatDBL("contNonBasicTab", contNonBasicTab, mTab, nTab);
  rs_printmatINT("checked", checked, mTab, mTab);
  rs_printvecINT("changed", changed, mTab);
  rs_printvecDBL("norm", norm, mTab);
  rs_printmatINT("pi_mat", pi_mat, mTab, mTab);
#endif

#ifdef RS_TRACE
  sum_norms = 0;
  for(i=0; i<mTab; i++) {
    sum_norms += norm[i];
  }

  printf("CglRedSplit::reduce_contNonBasicTab():Final sum of norms: %f\n", sum_norms);
#endif

  delete[] norm;
  delete[] changed;
  rs_deallocmatINT(&checked, mTab, mTab);

} /* reduce_contNonBasicTab */

/************************************************************************/
void CglRedSplit::generate_row(int index_row, double *row) {

  int i;
  for(i=0; i<ncol+nrow; i++) {
    row[i] = 0;
  }
  if(!param.getUSE_CG2()) { 
       // coeff will become zero in generate_cgcut_2() anyway
    for(i=0; i<card_intBasicVar_frac; i++) {
      row[intBasicVar_frac[i]] += pi_mat[index_row][i];
    }
  }
  for(i=0; i<card_intNonBasicVar; i++) {
    int locind = intNonBasicVar[i];
    row[locind] = 0;
    int j;
    for(j=0; j<mTab; j++) {
      row[locind] += pi_mat[index_row][j] * intNonBasicTab[j][i];
    }
  }
  for(i=0; i<card_contNonBasicVar; i++) {
    row[contNonBasicVar[i]] = contNonBasicTab[index_row][i];
  }
} /* generate_row */

/************************************************************************/
int CglRedSplit::generate_cgcut(double *row, double *rhs) {
  
  double f0 = rs_above_integer(*rhs);
  double f0compl = 1 - f0;

#ifdef RS_TRACEALL
  rs_printvecDBL("CglRedSplit::generate_cgcut(): starting row", 
		 row, ncol+nrow);
  printf("rhs: %f\n", *rhs);
#endif

  // See Wolsey "Integer Programming" (1998), p. 130, second line of proof of 
  // Proposition 8.8

  if((f0 < param.getAway()) || (f0compl < param.getAway())) {
    return(0);
  }

  int i;
  for(i=0; i<card_intNonBasicVar; i++) {
    int locind = intNonBasicVar[i];
    double f = rs_above_integer(row[locind]);
    row[locind] -= f;
    if(f > f0) {
      row[locind] += (f-f0)/f0compl;
    }
  }

  for(i=0; i<card_contNonBasicVar; i++) {
    if(row[contNonBasicVar[i]] < 0) {
      row[contNonBasicVar[i]] /= f0compl;
    }
    else {
      row[contNonBasicVar[i]] = 0;
    }
  }
  (*rhs) -= f0;
  
#ifdef RS_TRACEALL
  rs_printvecDBL("CglRedSplit::generate_cgcut: row", row, ncol+nrow);
  printf("rhs: %f\n", *rhs);
#endif

  return(1);
} /* generate_cgcut */

/************************************************************************/
int CglRedSplit::generate_cgcut_2(int/* basic_ind*/, double *row, double *rhs) {
  
#ifdef RS_TRACEALL
  rs_printvecDBL("CglRedSplit::generate_cgcut_2(): starting row", 
		 row, ncol+nrow);
  printf("rhs: %f\n", *rhs);
#endif

  // Note that basic integer variables with fractional value
  // other than basic_ind might have a non zero integer coefficient in the 
  // "real" row. However, these coefficients become zero anyway in the
  // cut. It is thus assumed that all basic integer variables have a
  // coefficient of zero in the given row (enforced in generate_row()).
  //
  // Other integer variables and basic continuous variables
  // have a coefficient of zero, as the corresponding row is not selected
  // for combination.

  double f0 = rs_above_integer(*rhs);
  double f0compl = 1 - f0;

  // See Wolsey "Integer Programming" (1998), p. 130, fourth line from top
  // after correcting typo (Proposition 8.8), flipping all signs to get <=.

  if((f0 < param.getAway()) || (f0compl < param.getAway())) {
    return(0);
  }

  double ratf0f0compl = f0/f0compl;
  int i;

  /** Not needed since coeff are already zero 
  for(i=0; i<card_intBasicVar_frac; i++) { 
         // extra integer coefficients compared to formula.
         // They all become zero.
    int locind = intBasicVar_frac[i];
    row[locind] = 0;
  }
  **/

  for(i=0; i<card_intNonBasicVar; i++) {
    int locind = intNonBasicVar[i];
    double f = rs_above_integer(row[locind]);
    double fcompl = 1-f;

    if(f > f0) {
      row[locind] = -ratf0f0compl * fcompl;
    }
    else {
      row[locind] = -f;
    }
  }

  for(i=0; i<card_contNonBasicVar; i++) {
    int locind = contNonBasicVar[i];
    if(row[locind] < 0) {
      row[locind] *= ratf0f0compl;
    }
    else {
      row[locind] = -row[locind];
    }
  }
  (*rhs) = -f0;
  
#ifdef RS_TRACEALL
  rs_printvecDBL("CglRedSplit::generate_cgcut_2(): row", row, ncol+nrow);
  printf("rhs: %f\n", *rhs);
#endif

  return(1);
} /* generate_cgcut_2 */

/************************************************************************/
void CglRedSplit::eliminate_slacks(double *row, 
				   const double *elements, 
				   const CoinBigIndex *rowStart,
				   const int *indices,
				   const int *rowLength,
				   const double *rhs, double *tabrowrhs) {

  for(int i=0; i<nrow; i++) {
    if(fabs(row[ncol+i]) > param.getEPS_ELIM()) {

      if(rowLower[i] > rowUpper[i] - param.getEPS()) {
	row[ncol+i] = 0;
	continue;
      }

      CoinBigIndex upto = rowStart[i] + rowLength[i];
      for(CoinBigIndex j=rowStart[i]; j<upto; j++) {
	row[indices[j]] -= row[ncol+i] * elements[j];
      }
      *tabrowrhs -= row[ncol+i] * rhs[i];
    }
  }

#ifdef RS_TRACEALL
  rs_printvecDBL("CglRedSplit::eliminate_slacks: row", row, ncol+nrow);
  printf("rhs: %f\n", *tabrowrhs);
#endif

} /* eliminate_slacks */

/************************************************************************/
void CglRedSplit::flip(double *row) {
  
  int i;
  for(i=0; i<card_nonBasicAtUpper; i++) {
    row[nonBasicAtUpper[i]] = -row[nonBasicAtUpper[i]];
  }
} /* flip */

/************************************************************************/
void CglRedSplit::unflip(double *row, double *tabrowrhs, double *slack_val) {
  
  int i;
  for(i=0; i<card_nonBasicAtLower; i++) {
    int locind = nonBasicAtLower[i];
    if(locind < ncol) {
      *tabrowrhs += row[locind] * colLower[locind];
    }
    else {
      *tabrowrhs += row[locind] * slack_val[locind-ncol];
    }
  }
  for(i=0; i<card_nonBasicAtUpper; i++) {
    int locind = nonBasicAtUpper[i];
    row[locind] = -row[locind];
    if(locind < ncol) {
      *tabrowrhs += row[locind] * colUpper[locind];
    }
    else {
      *tabrowrhs += row[locind] * slack_val[locind-ncol];
    }
  }

#ifdef RS_TRACEALL
  rs_printvecDBL("After unflip: row", row, ncol+nrow);
  printf("rhs: %f\n", *tabrowrhs);
#endif

} /* unflip */

/************************************************************************/
double CglRedSplit::row_scale_factor(double *row) {

  int i, has_lub = 0, nelem = 0;
  double val,  norm = 0, max_val = 0, min_val = param.getINFINIT();

  for(i=0; i<ncol; i++) {
    val = fabs(row[i]);
    max_val = std::max(max_val, val);
    norm += val * val;

    if(low_is_lub[i] + up_is_lub[i]) {
      if(val > param.getEPS_COEFF_LUB()) {
	min_val = std::min(min_val, val);
	has_lub = 1;
	nelem++;
      }
    }
    else {
      if(val > param.getEPS_COEFF()) {
	min_val = std::min(min_val, val);
	nelem++;
     }
    }
  }

  double retval = 1;

  if(norm > 100 * nelem) {
    retval = 10 * sqrt(norm / nelem);
  }
  if(norm < 0.5 * nelem) {
    retval = 0.5 * sqrt(norm / nelem);
  }

  if((retval < 0.02) || (retval > 50)) {
    return(-1);
  }

  if(has_lub) {
    if((max_val > param.getEPS_COEFF_LUB()) && 
       (max_val < param.getMAXDYN_LUB() * min_val) && 
       (max_val >= min_val)) {
      return(retval);
    }
  }
  else {
    if((max_val > param.getEPS_COEFF()) && 
       (max_val < param.getMAXDYN() * min_val) && 
       (max_val >= min_val)) {
      return(retval);
    }
  }

#ifdef RS_TRACE
  printf("CglRedSplit::scale_row(): max_val: %6.6f   min_val: %6.6f\n",
	 max_val, min_val);
#endif

  return(-1);
} /* row_scale_factor */

/************************************************************************/
int CglRedSplit::generate_packed_row(const double *lclXlp,
				     double *row,
				     int *rowind, double *rowelem, 
				     int *card_row, double & rhs) {
  int i;
  double scale_f = row_scale_factor(row);
  double value;

  if(scale_f < 0) {

#ifdef RS_TRACE
    printf("CglRedSplit::generate_packed_row(): Cut discarded (bad numerical behavior)\n");
#endif	

    return(0);
  }

  *card_row = 0;
  rhs /= scale_f;

  for(i=0; i<ncol; i++) {
    value = row[i]/scale_f;
    if(fabs(value) > param.getEPS_COEFF()) {
      rowind[*card_row] = i;
      rowelem[*card_row] = value;
      (*card_row)++;
      if(*card_row > param.getMAX_SUPPORT()) {

#ifdef RS_TRACE
	printf("CglRedSplit::generate_packed_row(): Cut discarded (too many non zero)\n");
#endif	
	return(0);
      }
    } else {
      if((value > 0.0) && (!low_is_lub[i])) {
        rhs -= value * colLower[i];
      } 
      else if((value < 0.0) && (!up_is_lub[i])) {
        rhs -= value * colUpper[i];
      } 
      else if(fabs(value) > param.getEPS_COEFF_LUB()) {
        // take anyway
        rowind[*card_row] = i;
        rowelem[*card_row] = value;
        (*card_row)++;
        if(*card_row > param.getMAX_SUPPORT()) {
#ifdef RS_TRACE
          printf("CglRedSplit::generate_packed_row(): Cut discarded since too many non zero coefficients\n");
#endif	
          return(0);
        }
      }
    }
  }
  value = 0;
  for(i=0; i<(*card_row); i++) {
    value += lclXlp[rowind[i]] * rowelem[i];
  }

  if(value > rhs) {
    value -= rhs;
    if(value < param.getMINVIOL()*(*card_row)) {

#ifdef RS_TRACE
      printf("CglRedSplit::generate_packed_row(): Cut discarded: violation: %12.10f\n", value);
#endif	
      
      return(0);
    }
  }
  
  return(1);
} /* generate_packed_row */

// TO BE REMOVED
/***********************************************************************/
void CglRedSplit::setLimit(int limit)
{
  if (limit>0)
    param.setMAX_SUPPORT(limit);
} /* setLimit */

/***********************************************************************/
int CglRedSplit::getLimit() const
{
  return param.getMAX_SUPPORT();
} /* getLimit */

/***********************************************************************/
void CglRedSplit::setAway(double value)
{
  if (value>0.0 && value<=0.5)
    param.setAway(value);
}

/***********************************************************************/
double CglRedSplit::getAway() const
{
  return param.getAway();
}

/***********************************************************************/
void CglRedSplit::setMaxTab(double value)
{
  if (value > 10) {
    param.setMaxTab(value);
  }
  else {
    printf("### WARNING: CglRedSplit::setMaxTab(): value: %f ignored\n", 
	   value);
  }
}

/***********************************************************************/
double CglRedSplit::getMaxTab() const
{
  return param.getMaxTab();
}

/***********************************************************************/
void CglRedSplit::setLUB(double value)
{
  if (value > 0.0) {
    param.setLUB(value);
  }
  else {
    printf("### WARNING: CglRedSplit::setLUB(): value: %f ignored\n", value);
  }
} /* setLUB */

/***********************************************************************/
double CglRedSplit::getLUB() const
{
  return param.getLUB();
} /* getLUB */

/***********************************************************************/
void CglRedSplit::setEPS(double value)
{
  if (value>0.0 && value<=0.1) {
    param.setEPS(value);
  }
  else {
    printf("### WARNING: CglRedSplit::setEPS(): value: %f ignored\n", value);
  }
} /* setEPS */

/***********************************************************************/
double CglRedSplit::getEPS() const
{
  return param.getEPS();
} /* getEPS */

/***********************************************************************/
void CglRedSplit::setEPS_COEFF(double value)
{
  if (value > 0.0 && value <= 0.1) {
    param.setEPS_COEFF(value);
  }
  else {
    printf("### WARNING: CglRedSplit::setEPS_COEFF(): value: %f ignored\n", 
	   value);
  }
} /* setEPS_COEFF */

/***********************************************************************/
double CglRedSplit::getEPS_COEFF() const
{
  return param.getEPS_COEFF();
} /* getEPS_COEFF */

/***********************************************************************/
void CglRedSplit::setEPS_COEFF_LUB(double value)
{
  if (value > 0.0 && value <= 0.1) {
    param.setEPS_COEFF_LUB(value);
  }
  else {
    printf("### WARNING: CglRedSplit::setEPS_COEFF_LUB(): value: %f ignored\n", 
	   value);
  }
} /* setEPS_COEFF_LUB */

/***********************************************************************/
double CglRedSplit::getEPS_COEFF_LUB() const
{
  return param.getEPS_COEFF_LUB();
} /* getEPS_COEFF_LUB */

/***********************************************************************/
void CglRedSplit::setEPS_RELAX(double value)
{
  if (value > 0.0 && value <= 0.1) {
    param.setEPS_RELAX_ABS(value);
  }
  else {
    printf("### WARNING: CglRedSplit::setEPS_RELAX(): value: %f ignored\n", 
	   value);
  }
} /* setEPS_RELAX */

/***********************************************************************/
double CglRedSplit::getEPS_RELAX() const
{
  return param.getEPS_RELAX_ABS();
} /* getEPS_RELAX */

/***********************************************************************/
void CglRedSplit::setNormIsZero(double value)
{
  if (value>0.0 && value<=1) {
    param.setNormIsZero(value);
  }
  else {
    printf("### WARNING: CglRedSplit::setNormIsZero(): value: %f ignored\n",
	   value);
  }
} /* setNormIsZero */

/***********************************************************************/
double CglRedSplit::getNormIsZero() const
{
  return param.getNormIsZero();
} /* getNormIsZero */

/***********************************************************************/
void CglRedSplit::setMinReduc(double value)
{
  if (value>0.0 && value<=1) {
    param.setMinReduc(value);
  }
  else {
    printf("### WARNING: CglRedSplit::MinReduc(): value: %f ignored\n",
	   value);
  }
} /* setMinReduc */

/***********************************************************************/
double CglRedSplit::getMinReduc() const
{
  return param.getMinReduc();
} /* getMinReduc */

/************************************************************************/
void CglRedSplit::set_given_optsol(const double *given_sol, int card_sol) {
  given_optsol = given_sol;
  card_given_optsol = card_sol;
} /* set_given_optsol */

/************************************************************************/
void CglRedSplit::check_optsol(const int calling_place,
			       const double *lclXlp, const double *slack_val,
			       const int do_flip) {

  if(card_given_optsol != ncol) {
    printf("### ERROR: CglRedSplit(): card_given_optsol: %d  ncol: %d\n", 
	   card_given_optsol, ncol);
    exit(1);
  }

  int i;
  double *ck_slack = new double[nrow];

#ifdef RS_TRACEALL
  print();
#endif

  byRow->timesMinor(given_optsol, ck_slack);
  int irow;
  for(irow=0; irow<nrow; irow++) {
    ck_slack[irow] = rowRhs[irow] - ck_slack[irow];  
                                        // slack values for optimal solution
  }
  
  double *ck_row = new double[ncol+nrow];
  
  for(irow=0; irow<mTab; irow++) {
    for(i=0; i<ncol+nrow; i++) {
      ck_row[i] = 0;
    }
    for(i=0; i<card_intBasicVar_frac; i++) {
      ck_row[intBasicVar_frac[i]] = pi_mat[irow][i];
    }
    for(i=0; i<card_intNonBasicVar; i++) {
      ck_row[intNonBasicVar[i]] = 0;
      int j;
      for(j=0; j<mTab; j++) {
	ck_row[intNonBasicVar[i]] += pi_mat[irow][j] * intNonBasicTab[j][i];
      }
    }
    for(i=0; i<card_contNonBasicVar; i++) {
      ck_row[contNonBasicVar[i]] = contNonBasicTab[irow][i];
    }

    double adjust_rhs = 0;
    if(do_flip) {
      for(i=0; i<card_nonBasicAtLower; i++) {
	int locind = nonBasicAtLower[i];
	if(locind < ncol) {
	  adjust_rhs += ck_row[locind] * colLower[locind];
	}
	else {
	  adjust_rhs += ck_row[locind] * slack_val[locind-ncol];
	}
      }
      for(i=0; i<card_nonBasicAtUpper; i++) {
	int locind = nonBasicAtUpper[i];
	ck_row[locind] = -ck_row[locind];
	if(locind < ncol) {
	  adjust_rhs += ck_row[locind] * colUpper[locind];
	}
	else {
	  adjust_rhs += ck_row[locind] * slack_val[locind-ncol];
	}
      }
    }

    double ck_lhs = rs_dotProd(ck_row, given_optsol, ncol);
    ck_lhs += rs_dotProd(&(ck_row[ncol]), ck_slack, nrow);
    
    double ck_rhs = adjust_rhs + rs_dotProd(ck_row, lclXlp, ncol);
    ck_rhs += rs_dotProd(&(ck_row[ncol]), slack_val, nrow);
    
#ifdef RS_TRACEALL
    rs_printvecDBL("ck_row", ck_row, ncol);
    rs_printvecDBL("given_optsol", given_optsol, ncol);
    rs_printvecDBL("ck_row(slacks)", &(ck_row[ncol]), nrow);
    rs_printvecDBL("ck_slack", ck_slack, nrow);
    printf("ck_rhs: %12.8f\n", ck_rhs); 
#endif

    if((ck_lhs < ck_rhs - param.getEPS()) || (ck_lhs > ck_rhs + param.getEPS())) {
      printf("### ERROR: CglRedSplit::check_optsol(): Cut %d cuts given_optsol\n", 
	     irow);
      rs_printvecDBL("ck_row", ck_row, ncol+nrow);
      printf("lhs: %f  rhs: %f    calling_place: %d\n", 
	     ck_lhs, ck_rhs, calling_place);
      exit(1);
    }
  }
  delete[] ck_slack;
  delete[] ck_row;
} /* check_optsol */

/************************************************************************/
void CglRedSplit::check_optsol(const int calling_place,
			       const double * /*lclXlp*/, const double *slack_val,
			       const double *ck_row, const double ck_rhs,
			       const int cut_number, const int do_flip) {

  if(card_given_optsol != ncol) {
    printf("### ERROR: CglRedSplit(): card_given_optsol: %d  ncol: %d\n", 
	   card_given_optsol, ncol);
    exit(1);
  }

  double *cpy_row = new double[ncol+nrow];
  double *ck_slack = new double[nrow];

#ifdef RS_TRACEALL
  print();
#endif

  int i, irow;
  for(i=0; i<ncol+nrow; i++) {
    cpy_row[i] = ck_row[i];
  }

  byRow->timesMinor(given_optsol, ck_slack);
  for(irow=0; irow<nrow; irow++) {
    ck_slack[irow] = rowRhs[irow] - ck_slack[irow];  
                                       // slack values for optimal solution
  }
  
  double adjust_rhs = 0;
  if(do_flip) {
    for(i=0; i<card_nonBasicAtLower; i++) {
      int locind = nonBasicAtLower[i];
      if(locind < ncol) {
	adjust_rhs += cpy_row[locind] * colLower[locind];
      }
      else {
	adjust_rhs += cpy_row[locind] * slack_val[locind-ncol];
      }
    }
    for(i=0; i<card_nonBasicAtUpper; i++) {
      int locind = nonBasicAtUpper[i];
      cpy_row[locind] = -cpy_row[locind];
      if(locind < ncol) {
	adjust_rhs += cpy_row[locind] * colUpper[locind];
      }
      else {
	adjust_rhs += cpy_row[locind] * slack_val[locind-ncol];
      }
    }
  }


  double ck_lhs = rs_dotProd(cpy_row, given_optsol, ncol);
  ck_lhs += rs_dotProd(&(cpy_row[ncol]), ck_slack, nrow);
    
  if(ck_lhs > ck_rhs + adjust_rhs + param.getEPS()) {
    printf("### ERROR: CglRedSplit::check_optsol(): Cut %d cuts given_optsol\n", 
	   cut_number);
    rs_printvecDBL("cpy_row", cpy_row, ncol+nrow);
    printf("lhs: %f  rhs: %f    calling_place: %d\n", 
	   ck_lhs, ck_rhs + adjust_rhs, calling_place);
    exit(1);
  }
  delete[] cpy_row;
  delete[] ck_slack;
} /* check_optsol */

/************************************************************************/
bool CglRedSplit::rs_are_different_vectors(const int *vect1, 
					   const int *vect2,
					   const int dim) {
  int i;
  for(i=0; i<dim; i++) {
    if(vect1[i] != vect2[i]) {
      printf("### ERROR: rs_are_different_vectors(): vect1[%d]: %d vect2[%d]: %d\n", i, vect1[i], i, vect2[i]);
      return(0);
    }    
  }
  return(1);
} /* rs_are_different_vectors */

/************************************************************************/
bool CglRedSplit::rs_are_different_vectors(const double *vect1, 
					   const double *vect2,
					   const int dim) {
  int i;
  for(i=0; i<dim; i++) {
    if(fabs(vect1[i] - vect2[i]) > 1e-6) {
      printf("### ERROR: rs_are_different_vectors(): vect1[%d]: %12.8f vect2[%d]: %12.8f\n", i, vect1[i], i, vect2[i]);
      return(0);
    }    
  }
  return(1);
} /* rs_are_different_vectors */

/************************************************************************/
bool CglRedSplit::rs_are_different_matrices(const CoinPackedMatrix *mat1, 
					    const CoinPackedMatrix *mat2,
					    const int nmaj,
					    const int /*nmin*/) {
  
  const CoinBigIndex *matStart1 = mat1->getVectorStarts();
  const double *matElements1 = mat1->getElements();
  const int *matIndices1 = mat1->getIndices();
  const int *matRowLength1 = mat1->getVectorLengths(); 

  const CoinBigIndex *matStart2 = mat2->getVectorStarts();
  const double *matElements2 = mat2->getElements();
  const int *matIndices2 = mat2->getIndices();
  const int *matRowLength2 = mat2->getVectorLengths(); 

  int i;
  CoinBigIndex j;

  for(i=0; i<nmaj; i++) {
    if(matStart1[i] != matStart2[i]) {
      printf("### ERROR: rs_are_different_matrices(): matStart1[%d]: %d matStart2[%d]: %d\n", i, matStart1[i], i, matStart2[i]);
      return(1);
    }
    if(matRowLength1[i] != matRowLength2[i]) {
      printf("### ERROR: rs_are_different_matrices(): matRowLength1[%d]: %d matRowLength2[%d]: %d\n", i, matRowLength1[i], i, matRowLength2[i]);
      return(1);
    }

    for(j=matStart1[i]; j<matStart1[i]+matRowLength1[i]; j++) {
      if(matIndices1[j] != matIndices2[j]) {
	printf("### ERROR: rs_are_different_matrices(): matIndices1[%d]: %d matIndices2[%d]: %d\n", j, matIndices1[j], j, matIndices2[j]);
	return(1);
      }
      if(fabs(matElements1[j] - matElements2[j]) > 1e-6) {
	printf("### ERROR: rs_are_different_matrices(): matElements1[%d]: %12.8f matElements2[%d]: %12.8f\n", j, matElements1[j], j, matElements2[j]);
	return(1);
      }
    }
  }
  return(0);
} /* rs_are_different_matrices */

/************************************************************************/
void CglRedSplit::generateCuts(const OsiSolverInterface &si, OsiCuts & cs,
			       const CglTreeInfo )
{
  solver = const_cast<OsiSolverInterface *>(&si);
  if(solver == NULL) {
    printf("### WARNING: CglRedSplit::generateCuts(): no solver available.\n");
    return;    
  }  

  if(!solver->optimalBasisIsAvailable()) {
    printf("### WARNING: CglRedSplit::generateCuts(): no optimal basis available.\n");
    return;
  }

  // Reset some members of CglRedSplit
  card_intBasicVar_frac = 0;
  card_intNonBasicVar = 0;
  card_contNonBasicVar = 0;
  card_nonBasicAtUpper = 0;
  card_nonBasicAtLower = 0;

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
  colType = NULL;
  byRow = solver->getMatrixByRow();

  solver->enableFactorization();
  generateCuts(cs);
  solver->disableFactorization();
} /* generateCuts */

/************************************************************************/
void CglRedSplit::generateCuts(OsiCuts &cs)
{
  int i;
  low_is_lub = new int[ncol]; 
  up_is_lub = new int[ncol]; 
  is_integer = new int[ncol]; 
  
  compute_is_lub();
  compute_is_integer();

  int *cstat = new int[ncol];
  int *rstat = new int[nrow];
  solver->getBasisStatus(cstat, rstat);   // 0: free  1: basic  
                                          // 2: upper 3: lower

  int *basis_index = new int[nrow]; // basis_index[i] = 
                                    //        index of pivot var in row i
                                    //        (slack if number >= ncol) 

#ifdef RS_TRACETAB
  rs_printvecINT("cstat", cstat, ncol);
  rs_printvecINT("rstat", rstat, nrow);
#endif

  solver->getBasics(basis_index);

  cv_intBasicVar_frac = new int[ncol];  
  intBasicVar_frac = new int[ncol];                                 
  intNonBasicVar = new int[ncol];       
  contNonBasicVar = new int[ncol+nrow]; 
  nonBasicAtUpper = new int[ncol+nrow]; 
  nonBasicAtLower = new int[ncol+nrow]; 
  double dist_int;
  /*
    Even if there are no free variables, clp may mark them as free
    if bounds are very large - which might be generated by a badly
    scaled problem.
   */
  bool goodModel = true;
  for(i=0; i<ncol; i++) {
    cv_intBasicVar_frac[i] = 0;

    switch(cstat[i]) {
    case 1: // basic variable
      
      dist_int = rs_above_integer(xlp[i]);
      if(is_integer[i] && 
	 (dist_int > param.getAway()) && (dist_int < 1 - param.getAway())) {
	cv_intBasicVar_frac[i] = 1;
	card_intBasicVar_frac++;

	// intBasicVar_frac computed below, 
	// as order must be according to selected rows

      }
      break;

    case 2: // Non basic at upper bound: must be flipped and shifted
            // so that it becomes non negative with lower bound 0 
            // It is assumed that bounds for integer variables have been
            // tightend so that non basic integer structural variables 
            // have integer values

      nonBasicAtUpper[card_nonBasicAtUpper] = i;
      card_nonBasicAtUpper++;

      if(is_integer[i]) {
	intNonBasicVar[card_intNonBasicVar] = i;
	card_intNonBasicVar++;
      }
      else {
	contNonBasicVar[card_contNonBasicVar] = i;
	card_contNonBasicVar++;
      }
      break;

    case 3 : // non basic at lower bound: must be shifted so that it becomes
             // non negative with lower bound 0
            // It is assumed that bounds for integer variables have been
            // tightend so that they are integer

      nonBasicAtLower[card_nonBasicAtLower] = i;
      card_nonBasicAtLower++;

      if(is_integer[i]) {
	intNonBasicVar[card_intNonBasicVar] = i;
	card_intNonBasicVar++;
      }
      else {
	contNonBasicVar[card_contNonBasicVar] = i;
	card_contNonBasicVar++;
      }
      break;

    default: // free variable ? Don't know how to handle 
#ifndef NDEBUG
      printf("### ERROR: CglRedSplit::generateCuts(): cstat[%d]: %d\n",
	     i, cstat[i]);
#endif
      goodModel = false;
      break;
    } 
  }

  for(i=0; i<nrow; i++) {
    switch(rstat[i]) {
    case 1: // basic slack
      break;

    case 2: // non basic slack at upper; flipped and shifted
      nonBasicAtUpper[card_nonBasicAtUpper] = ncol+i;
      card_nonBasicAtUpper++;

      contNonBasicVar[card_contNonBasicVar] = ncol+i;
      card_contNonBasicVar++;
      break;

    case 3: // non basic slack at lower; shifted
      nonBasicAtLower[card_nonBasicAtLower] = ncol+i;
      card_nonBasicAtLower++;

      contNonBasicVar[card_contNonBasicVar] = ncol+i;
      card_contNonBasicVar++;
      break;

    default:
#ifndef NDEBUG
      printf("### ERROR: CglRedSplit::generateCuts(): rstat[%d]: %d\n",
	     i, rstat[i]);
#endif
      goodModel = false;
      break;
    }
  }

  if((card_contNonBasicVar == 0) || (card_intBasicVar_frac == 0)
     || !goodModel) {
    delete[] cstat;
    delete[] rstat;
    delete[] basis_index;

    delete[] cv_intBasicVar_frac;  
    delete[] intBasicVar_frac;
    delete[] intNonBasicVar;
    delete[] contNonBasicVar;
    delete[] nonBasicAtUpper;
    delete[] nonBasicAtLower;
    delete[] low_is_lub;
    delete[] up_is_lub;
    delete[] is_integer;

    return; // no cuts can be generated
  }

  /* Loop is mTab * mTab * std::max(mTab, nTab) so may be very expensive. 
     Reduce mTab if the above value is larger than maxTab_ */

  int new_mTab = card_intBasicVar_frac;
  double nc = static_cast<double> (card_contNonBasicVar);
  double nc3 = nc * nc * nc;

  if(nc3 > param.getMaxTab()) {
    new_mTab = static_cast<int> (sqrt(param.getMaxTab()/card_contNonBasicVar));
  }
  else {
#if defined(_MSC_VER)
    new_mTab = static_cast<int> (pow(param.getMaxTab(), 1./3.));
#else
    new_mTab = static_cast<int> (cbrt(param.getMaxTab()));
#endif
  }

  if(new_mTab == 0) {
    delete[] cstat;
    delete[] rstat;
    delete[] basis_index;
    
    delete[] cv_intBasicVar_frac;  
    delete[] intBasicVar_frac;
    delete[] intNonBasicVar;
    delete[] contNonBasicVar;
    delete[] nonBasicAtUpper;
    delete[] nonBasicAtLower;
    delete[] low_is_lub;
    delete[] up_is_lub;
    delete[] is_integer;
    return; // no cuts can be generated
  }

  int start = 0;  // first row for selecting intBasicVar_frac, if too many.
                  // Removing rows used for generation whose pivot is some 
                  // var in intBasicVar_frac is valid since the corresponding 
                  // column in the optimal tableau is a column of the 
                  // identity matrix 

  if(new_mTab < card_intBasicVar_frac) {
    // Remove some of the rows used for generation.
    // Poor randomness; could do better if needed

    int seed = card_intBasicVar_frac;
    start = static_cast<int> (nrow * rs_genalea(&seed));
    
#ifdef RS_TRACE
    printf("CglRedSplit::generateCuts(): mTab: %d  new_mTab: %d\n", 
	   card_intBasicVar_frac, new_mTab);
#endif

    card_intBasicVar_frac = new_mTab;
  }
  
double *slack_val = new double[nrow];

  for(i=0; i<nrow; i++) {
    slack_val[i] = rowRhs[i] - rowActivity[i];
  }

#ifdef RS_DEBUG
  double *solver_rhs = solver->getRightHandSide();
  if(rs_are_different_vectors(rowRhs, solver_rhs, nrow)) {
    printf("### ERROR: CglRedSplit::generateCuts(): rowRhs[%d]: %12.f  solver_rhs[%d]: %12.f\n", i, rowRhs[i], i, solver_rhs[i]);
    exit(1);
  }
#endif

  double *z = new double[ncol];  // workspace to get row of the tableau
  double *slack = new double[nrow];  // workspace to get row of the tableau

#ifdef RS_TRACETAB
  printOptTab(solver);
#endif

  mTab = card_intBasicVar_frac;
  nTab = card_contNonBasicVar;

  rhsTab = new double[mTab];
  int card_rowTab = 0;

  rs_allocmatDBL(&contNonBasicTab, mTab, nTab);
  rs_allocmatDBL(&intNonBasicTab, mTab, card_intNonBasicVar);

  card_intBasicVar_frac = 0; // recompute in pivot order

  for(i=0; i<nrow; i++) {

    if(mTab <= card_intBasicVar_frac) { 
      break;
    }

    int ind_row = start + i;
    if(ind_row > nrow) {
      ind_row -= nrow;
    }

    if(basis_index[ind_row] >= ncol) {
      continue;
    } 

    if(cv_intBasicVar_frac[basis_index[ind_row]] == 1) { 
                                                  // row used in generation
      intBasicVar_frac[card_intBasicVar_frac] = basis_index[ind_row];
      card_intBasicVar_frac++;
      rhsTab[card_rowTab] = xlp[basis_index[ind_row]];
      solver->getBInvARow(ind_row, z, slack);
      int ii;
      for(ii=0; ii<card_contNonBasicVar; ii++) {
	int locind = contNonBasicVar[ii];
	if(locind < ncol) {
	  contNonBasicTab[card_rowTab][ii] = z[locind];
	}
	else {
	  contNonBasicTab[card_rowTab][ii] = slack[locind - ncol];
	}
      }

      for(ii=0; ii<card_intNonBasicVar; ii++) {
	int locind = intNonBasicVar[ii];
	if(locind < ncol) {
	  intNonBasicTab[card_rowTab][ii] = z[locind];
	}
	else {
	  printf("### ERROR: CglRedSplit::generateCuts(): integer slack unexpected\n");
	  exit(1);
	}
      }

      card_rowTab++;
    }
  }

  rs_allocmatINT(&pi_mat, mTab, mTab);
  for(i=0; i<mTab; i++) {
    int ii;
    for(ii=0; ii<mTab; ii++) {
      pi_mat[i][ii] = 0;
    }
    pi_mat[i][i] = 1;
  }

#ifdef RS_TRACE
  printf("intBasicVar_frac:\n");
  for(i=0; i<card_intBasicVar_frac; i++) {
    printf("%d ", intBasicVar_frac[i]);
  }
  printf("\n");
  printf("intNonBasicVar:\n");
  for(i=0; i<card_intNonBasicVar; i++) {
    printf("%d ", intNonBasicVar[i]);
  }
  printf("\n");
  printf("contNonBasicVar:\n");
  for(i=0; i<card_contNonBasicVar; i++) {
    printf("%d ", contNonBasicVar[i]);
  }
  printf("\n");
#endif

  if(given_optsol) {
    check_optsol(1, xlp, slack_val, 0);
  }

  reduce_contNonBasicTab();

  if(given_optsol) {
    check_optsol(2, xlp, slack_val, 0);
  }

  int card_row;
  double *row = new double[ncol+nrow];
  int *rowind = new int[ncol];
  double *rowelem = new double[ncol];

  const double *elements = byRow->getElements();
  const CoinBigIndex *rowStart = byRow->getVectorStarts();
  const int *indices = byRow->getIndices();
  const int *rowLength = byRow->getVectorLengths(); 

  for(i=0; i<mTab; i++) {
    generate_row(i, row);
    flip(row);

    // RHS of equalities after flipping/translating non basic variables
    // is given by the current LP solution (re-ordered according to 
    // basis_index), i.e. what is now in rhsTab. 
    // RHS of row i of (pi_mat * contNonBasicTab) is then simply 
    // pi_mat[i] * rhsTab

    double tabrowrhs = rs_dotProd(pi_mat[i], rhsTab, mTab); 
    int got_one = 0;
    
    if(param.getUSE_CG2()) {
      got_one = generate_cgcut_2(intBasicVar_frac[i], row, &tabrowrhs);
    }
    else {
      got_one = generate_cgcut(row, &tabrowrhs);
    }
    
    if(got_one) {
      unflip(row, &tabrowrhs, slack_val);

      if(given_optsol) {
	check_optsol(3, xlp, slack_val, row, tabrowrhs, i, 0);
      }

      eliminate_slacks(row, elements, rowStart, indices, 
		       rowLength, rowRhs, &tabrowrhs);

      if(given_optsol) {
	check_optsol(4, xlp, slack_val, row, tabrowrhs, i, 0);
      }

      if(generate_packed_row(xlp, row, rowind, rowelem, &card_row, 
			     tabrowrhs)&&card_row<300) {
      	OsiRowCut rc;
	rc.setRow(card_row, rowind, rowelem);
	rc.setLb(-param.getINFINIT());
	double adjust = param.getEPS_RELAX_ABS();
	if(param.getEPS_RELAX_REL() > 0.0) {
	  adjust += fabs(tabrowrhs) * param.getEPS_RELAX_REL();
	}
	rc.setUb(tabrowrhs + adjust);   
                                  // relax the constraint slightly
	cs.insertIfNotDuplicate(rc, CoinAbsFltEq(param.getEPS_COEFF()));
      }
    }
  }

  delete[] cstat;
  delete[] rstat;
  delete[] basis_index;
  delete[] slack;
  delete[] z;
  delete[] slack_val;
  delete[] row;
  delete[] rowind;
  delete[] rowelem;

  delete[] cv_intBasicVar_frac;  
  delete[] intBasicVar_frac;
  delete[] intNonBasicVar;
  delete[] contNonBasicVar;
  delete[] nonBasicAtUpper;
  delete[] nonBasicAtLower;
  delete[] low_is_lub;
  delete[] up_is_lub;
  delete[] is_integer;
  rs_deallocmatDBL(&contNonBasicTab, mTab, nTab);
  rs_deallocmatDBL(&intNonBasicTab, mTab, card_intNonBasicVar);
  rs_deallocmatINT(&pi_mat, mTab, mTab);
  delete[] rhsTab;

  //return(cs.sizeRowCuts());
} /* generateCuts */

/***********************************************************************/
void CglRedSplit::setParam(const CglRedSplitParam &source) {
  param = source;
} /* setParam */

/***********************************************************************/
void CglRedSplit::compute_is_lub() {

  int i;
  for(i=0; i<ncol; i++) {
    low_is_lub[i] = 0;
    up_is_lub[i] = 0;
    if(fabs(colUpper[i]) > param.getLUB()) {
      up_is_lub[i] = 1;
    }
    if(fabs(colLower[i]) > param.getLUB()) {
      low_is_lub[i] = 1;
    }
  }
} /* compute_is_lub */

/***********************************************************************/
void CglRedSplit::compute_is_integer() {

  int i;
  if(colType != NULL) {
    for(i=0; i<ncol; i++) {
      if(colType[i] != 'C') {
	is_integer[i] = 1;
      }
      else {
	if((colUpper[i] - colLower[i] < param.getEPS()) && 
	   (rs_above_integer(colUpper[i]) < param.getEPS())) {
	  
	  // continuous variable fixed to an integer value
	  is_integer[i] = 1;
	}
	else {
	  is_integer[i] = 0;
	}
      }
    }
  }
  else {
    const char * intVar = solver->getColType();
    for(i=0; i<ncol; i++) {
      if(intVar[i]) {
	is_integer[i] = 1;
      }
      else {
	if((colUpper[i] - colLower[i] < param.getEPS()) && 
	   (rs_above_integer(colUpper[i]) < param.getEPS())) {
	  
	  // continuous variable fixed to an integer value
	  is_integer[i] = 1;
	}
	else {
	  is_integer[i] = 0;
	}
      }
    }    
  }
} /* compute_is_integer */

/***********************************************************************/
void CglRedSplit::print() const
{
  rs_printvecINT("intBasicVar_frac", intBasicVar_frac, card_intBasicVar_frac);
  rs_printmatINT("pi_mat", pi_mat, card_intBasicVar_frac, 
		 card_intBasicVar_frac);
  rs_printvecINT("intNonBasicVar", intNonBasicVar, card_intNonBasicVar);
  rs_printmatDBL("intNonBasicTab", intNonBasicTab, card_intBasicVar_frac, 
		 card_intNonBasicVar);
  rs_printvecINT("contNonBasicVar", contNonBasicVar, card_contNonBasicVar);
  rs_printmatDBL("contNonBasicTab", contNonBasicTab, card_intBasicVar_frac, 
		 card_contNonBasicVar);
  rs_printvecINT("nonBasicAtLower", nonBasicAtLower, card_nonBasicAtLower);
  rs_printvecINT("nonBasicAtUpper", nonBasicAtUpper, card_nonBasicAtUpper);

} /* print */

/***********************************************************************/
void CglRedSplit::printOptTab(OsiSolverInterface *lclSolver) const
{
  int i;
  int *cstat = new int[ncol];
  int *rstat = new int[nrow];

  lclSolver->getBasisStatus(cstat, rstat);   // 0: free  1: basic  
                                          // 2: upper 3: lower

  int *basis_index = new int[nrow]; // basis_index[i] = 
                                    //        index of pivot var in row i
                                    //        (slack if number >= ncol) 
  lclSolver->getBasics(basis_index);

  double *z = new double[ncol];  // workspace to get row of the tableau
  double *slack = new double[nrow];  // workspace to get row of the tableau
  double *slack_val = new double[nrow];

  for(i=0; i<nrow; i++) {
    slack_val[i] = rowRhs[i] - rowActivity[i];
  }

  const double *rc = lclSolver->getReducedCost();
  const double *dual = lclSolver->getRowPrice();
  const double *solution = lclSolver->getColSolution();

  rs_printvecINT("cstat", cstat, ncol);
  rs_printvecINT("rstat", rstat, nrow);
  rs_printvecINT("basis_index", basis_index, nrow);

  rs_printvecDBL("solution", solution, ncol);
  rs_printvecDBL("slack_val", slack_val, nrow);
  rs_printvecDBL("reduced_costs", rc, ncol);
  rs_printvecDBL("dual solution", dual, nrow);

  printf("Optimal Tableau:\n");

  for(i=0; i<nrow; i++) {
    lclSolver->getBInvARow(i, z, slack);
    int ii;
    for(ii=0; ii<ncol; ii++) {
      printf("%5.2f ", z[ii]);
    }
    printf(" | ");
    for(ii=0; ii<nrow; ii++) {
      printf("%5.2f ", slack[ii]);
    }
    printf(" | ");
    if(basis_index[i] < ncol) {
      printf("%5.2f ", solution[basis_index[i]]);
    }
    else {
      printf("%5.2f ", slack_val[basis_index[i]-ncol]);
    }
    printf("\n");
  }
  int ii;
  for(ii=0; ii<7*(ncol+nrow+1); ii++) {
    printf("-");
  }
  printf("\n");

  for(ii=0; ii<ncol; ii++) {
    printf("%5.2f ", rc[ii]);    
  }
  printf(" | ");
  for(ii=0; ii<nrow; ii++) {
    printf("%5.2f ", -dual[ii]);
  }
  printf(" | ");
  printf("%5.2f\n", -lclSolver->getObjValue());

  delete[] cstat;
  delete[] rstat;
  delete[] basis_index;
  delete[] slack;
  delete[] z;
  delete[] slack_val;
} /* printOptTab */

/*********************************************************************/
CglRedSplit::CglRedSplit() :
CglCutGenerator(),
nrow(0),
ncol(0),
card_intBasicVar_frac(0),
card_intNonBasicVar(0),
card_contNonBasicVar(0),
card_nonBasicAtUpper(0),
card_nonBasicAtLower(0),
cv_intBasicVar_frac(0),
intBasicVar_frac(0),
intNonBasicVar(0), 
contNonBasicVar(0),
nonBasicAtUpper(0),
nonBasicAtLower(0),
mTab(0),
nTab(0),
pi_mat(0),
contNonBasicTab(0),
intNonBasicTab(0),
rhsTab(0),
given_optsol(0),
card_given_optsol(0)
{}

/*********************************************************************/
CglRedSplit::CglRedSplit(const CglRedSplitParam &RS_param) :
CglCutGenerator(),
nrow(0),
ncol(0),
card_intBasicVar_frac(0),
card_intNonBasicVar(0),
card_contNonBasicVar(0),
card_nonBasicAtUpper(0),
card_nonBasicAtLower(0),
cv_intBasicVar_frac(0),
intBasicVar_frac(0),
intNonBasicVar(0), 
contNonBasicVar(0),
nonBasicAtUpper(0),
nonBasicAtLower(0),
mTab(0),
nTab(0),
pi_mat(0),
contNonBasicTab(0),
intNonBasicTab(0),
rhsTab(0),
given_optsol(0),
card_given_optsol(0)
{
  param = RS_param;
}
 
/*********************************************************************/
CglRedSplit::CglRedSplit (const CglRedSplit & source) :
  CglCutGenerator(source),
  param(source.param),
  nrow(0),
  ncol(0),
  card_intBasicVar_frac(0),
  card_intNonBasicVar(0),
  card_contNonBasicVar(0),
  card_nonBasicAtUpper(0),
  card_nonBasicAtLower(0),
  cv_intBasicVar_frac(NULL),
  intBasicVar_frac(NULL),
  intNonBasicVar(NULL), 
  contNonBasicVar(NULL),
  nonBasicAtUpper(NULL),
  nonBasicAtLower(NULL),
  mTab(0),
  nTab(0),
  pi_mat(NULL),
  contNonBasicTab(NULL),
  intNonBasicTab(NULL),
  rhsTab(NULL),
  given_optsol(source.given_optsol),
  card_given_optsol(source.card_given_optsol)
{}

/*********************************************************************/
CglCutGenerator *
CglRedSplit::clone() const
{
  return new CglRedSplit(*this);
}

/*********************************************************************/
CglRedSplit::~CglRedSplit ()
{}

/*********************************************************************/
CglRedSplit &
CglRedSplit::operator=(const CglRedSplit &source)
{  

  if (this != &source) {
    CglCutGenerator::operator=(source);
    param = source.param;
    given_optsol = source.given_optsol;
    card_given_optsol = source.card_given_optsol;
  }
  return *this;
}
/*********************************************************************/
// Returns true if needs optimal basis to do cuts
bool 
CglRedSplit::needsOptimalBasis() const
{
  return true;
}

/*********************************************************************/
// Create C++ lines to get to current state
std::string
CglRedSplit::generateCpp(FILE * fp) 
{
  CglRedSplit other;
  fprintf(fp,"0#include \"CglRedSplit.hpp\"\n");
  fprintf(fp,"3  CglRedSplit redSplit;\n");
  if (param.getMAX_SUPPORT()!=other.param.getMAX_SUPPORT())
    fprintf(fp,"3  redSplit.setLimit(%d);\n",param.getMAX_SUPPORT());
  else
    fprintf(fp,"4  redSplit.setLimit(%d);\n",param.getMAX_SUPPORT());
  if (param.getAway()!=other.param.getAway())
    fprintf(fp,"3  redSplit.setAway(%g);\n",param.getAway());
  else
    fprintf(fp,"4  redSplit.setAway(%g);\n",param.getAway());
  if (param.getLUB()!=other.param.getLUB())
    fprintf(fp,"3  redSplit.setLUB(%g);\n",param.getLUB());
  else
    fprintf(fp,"4  redSplit.setLUB(%g);\n",param.getLUB());
  if (param.getEPS()!=other.param.getEPS())
    fprintf(fp,"3  redSplit.set.EPS(%g);\n",param.getEPS());
  else
    fprintf(fp,"4  redSplit.setEPS(%g);\n",param.getEPS());
  if (param.getEPS_COEFF()!=other.param.getEPS_COEFF())
    fprintf(fp,"3  redSplit.setEPS_COEFF(%g);\n",param.getEPS_COEFF());
  else
    fprintf(fp,"4  redSplit.set.EPS_COEFF(%g);\n",param.getEPS_COEFF());
  if (param.getEPS_COEFF_LUB()!=other.param.getEPS_COEFF_LUB())
    fprintf(fp,"3  redSplit.set.EPS_COEFF_LUB(%g);\n",param.getEPS_COEFF_LUB());
  else
    fprintf(fp,"4  redSplit.set.EPS_COEFF_LUB(%g);\n",param.getEPS_COEFF_LUB());
  if (param.getEPS_RELAX_ABS()!=other.param.getEPS_RELAX_ABS())
    fprintf(fp,"3  redSplit.set.EPS_RELAX(%g);\n",param.getEPS_RELAX_ABS());
  else
    fprintf(fp,"4  redSplit.set.EPS_RELAX(%g);\n",param.getEPS_RELAX_ABS());
  if (param.getNormIsZero()!=other.param.getNormIsZero())
    fprintf(fp,"3  redSplit.setNormIsZero(%g);\n",param.getNormIsZero());
  else
    fprintf(fp,"4  redSplit.setNormIsZero(%g);\n",param.getNormIsZero());
  if (param.getMinReduc()!=other.param.getMinReduc())
    fprintf(fp,"3  redSplit.setMinReduc(%g);\n",param.getMinReduc());
  else
    fprintf(fp,"4  redSplit.setMinReduc(%g);\n",param.getMinReduc());
  if (param.getMaxTab()!=other.param.getMaxTab())
    fprintf(fp,"3  redSplit.setMaxTab(%g);\n",param.getMaxTab());
  else
    fprintf(fp,"4  redSplit.setMaxTab(%g);\n",param.getMaxTab());
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  redSplit.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  redSplit.setAggressiveness(%d);\n",getAggressiveness());
  return "redSplit";
}
