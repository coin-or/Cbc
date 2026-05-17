// Last edit: 04/03/10
//
// Name:     CglRedSplit2.cpp
// Author:   Giacomo Nannicini
//           Singapore University of Technology and Design
//           Singapore
//           email: nannicini@sutd.edu.sg
//           based on CglRedSplit by Francois Margot
// Date:     03/09/09
//-----------------------------------------------------------------------------
// Copyright (C) 2010, Giacomo Nannicini and others.  All Rights Reserved.

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <climits>
#include <cfloat>
#include <cassert>
#include <iostream>

/* Debug output */
//#define RS2_TRACE

/* Print optimal tableau and basis */
//#define RS2_TRACETAB

/* Use LAPACK instead of own code for solving linear systems */
//#define RS2_USE_LAPACK
// Sparse for intNonBasic 0,1,2
#define RS_FAST_INT 2
// Sparse for contNonBasic 0,1,2
#define RS_FAST_CONT 2
// Sparse for workNonBasic 0,1,2
#define RS_FAST_WORK 2

#include "CoinPragma.hpp"
#include "OsiSolverInterface.hpp"

#include "CglRedSplit2.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CoinFactorization.hpp"
#include "CoinFinite.hpp"

#define TINY 1e-20

#define rs2round(x) (floor((x)+0.5))

//-------------------------------------------------------------------
// Generate Reduce-and-Split cuts
//------------------------------------------------------------------- 

#ifdef RS2_USE_LAPACK
extern "C" void dgesv_( const int * , const int * , double * , const int * , int * , double * , const int * , int * );
#endif

/***************************************************************************/


// Utility functions and definitions for sorting
struct sortElement{
  int index;
  double cost;
};

// Return -1 if firstE has lower cost or same cost but lower index
int rs2_compareElements(const void* firstE, const void* secondE){
  const struct sortElement* a = static_cast<const struct sortElement*>(firstE);
  const struct sortElement* b = static_cast<const struct sortElement*>(secondE);
  if (a->cost < b->cost)
    return -1;
  else if (a->cost > b->cost)
    return 1;
  else if (a->index < b->index)
    return -1;
  else if (a->index > b->index)
    return 1;

  return 0;
}

/***************************************************************************/
// Returns (value - floor) but allowing for small errors;
// taken from the CglGomory cut generator
inline double CglRedSplit2::rs_above_integer(const double value) const
{
  double value2=floor(value);
  double value3=rs2round(value);
  if (fabs(value3-value)< param.getEPS() * (fabs(value3)+1.0))
    return 0.0;
  return value-value2;
} /* rs_above_integer */

/**********************************************************/
void CglRedSplit2::rs_allocmatINT(int ***v, int m, int n)
{
  *v = reinterpret_cast<int **> (calloc (m, sizeof(int *)));
  if (*v == NULL) {
    printf("###ERROR: INTEGER matrix allocation failed\n");
    exit(1);
  }

  for (int i = 0; i < m; ++i) {
    (*v)[i] = reinterpret_cast<int *> (calloc (n, sizeof(int)));
    if ((*v)[i] == NULL) {
      printf("###ERROR: INTEGER matrix allocation failed\n");
      exit(1);
    }
  }
} /* rs_allocmatINT */

/**********************************************************/
void CglRedSplit2::rs_deallocmatINT(int ***v, int m)
{
  for (int i = 0; i < m; ++i) {
    free(reinterpret_cast<void *> ((*v)[i]));
  }
  free(reinterpret_cast<void *> (*v));
} /* rs_deallocmatINT */

/**********************************************************/
void CglRedSplit2::rs_allocmatDBL(double ***v, int m, int n)
{
  *v = reinterpret_cast<double **> (calloc (m, sizeof(double *)));
  if (*v == NULL) {
    printf("###ERROR: DOUBLE matrix allocation failed\n");
    exit(1);
  }

  for (int i = 0; i < m; ++i){
    (*v)[i] = reinterpret_cast<double *> (calloc (n, sizeof(double)));
    if ((*v)[i] == NULL) {
      printf("###ERROR: DOUBLE matrix allocation failed\n");
      exit(1);
    }
  }
} /* rs_allocmatDBL */

/**********************************************************/
void CglRedSplit2::rs_deallocmatDBL(double ***v, int m)
{
  for (int i = 0; i < m; ++i){
    free(reinterpret_cast<void *> ((*v)[i]));
  }
  free(reinterpret_cast<void *> (*v));
} /* rs_deallocmatDBL */

/**********************************************************/
void CglRedSplit2::rs_printvecINT(const char *vecstr, const int *x, int n) const
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
void CglRedSplit2::rs_printvecDBL(const char *vecstr, 
				  const double *x, int n) const
{
  int num, fromto, upto, j, i;

  num = (n/10) + 1;
  printf("%s :\n", vecstr);
  for(j=0; j<num; j++){
    fromto = 10*j;
    upto = 10 * (j+1);
    if(n <= upto) upto = n;
    for(i=fromto; i<upto; i++)
      printf(" %7.5f", x[i]);
    printf("\n");
  }
  printf("\n");
} /* rs_printvecDBL */

/**********************************************************/
void CglRedSplit2::rs_printmatINT(const char *vecstr, const int * const *x, 
				  int m, int n) const
{
  printf("%s :\n", vecstr);
  for (int i = 0; i < m; ++i){
    for (int j = 0; j < n; ++j){
      printf(" %4d", x[i][j]);
    }
    printf("\n");
  }
  printf("\n");
} /* rs_printmatINT */


/**********************************************************/
void CglRedSplit2::rs_printmatDBL(const char *vecstr, const double * const *x, int m, int n) const
{
  printf("%s :\n", vecstr);
  for (int i = 0; i < m; ++i){
    for (int j = 0; j < n; ++j){
      printf(" %7.3f", x[i][j]);
    }
    printf("\n");
  }
  printf("\n");
} /* rs_printmatDBL */

/***************************************************************************/
double CglRedSplit2::rs_dotProd(const double *u, const double *v, int dim) const
{
  double result = 0;
  for (int i = 0; i < dim; ++i){
    result += u[i] * v[i];
  }
  return(result);
} /* rs_dotProd */

/***************************************************************************/
static double rs_sparseDotProd(const double *u, const double *v, 
			       const int * uInd, const int * vInd)
{
  double result = 0;
  if (uInd[0]<vInd[0]) {
    int n=uInd[0];
    uInd++;
    for (int i = 0; i < n; ++i) {
      int j=uInd[i];
      result += u[j] * v[j];
    }
  } else {
    int n=vInd[0];
    vInd++;
    for (int i = 0; i < n; ++i) {
      int j=vInd[i];
      result += u[j] * v[j];
    }
  }
  return(result);
} /* rs_sparseDotProd */

/***************************************************************************/
double CglRedSplit2::rs_dotProd(const int *u, const double *v, int dim) const {
  double result = 0;
  for (int i=0; i < dim; ++i){
    result += u[i] * v[i];
  }
  return(result);
} /* rs_dotProd */

/***************************************************************************/
// From Numerical Recipes in C: LU decomposition
int CglRedSplit2::ludcmp(double **a, int n, int *indx, 
			 double *d, double* vv) const {
  int ret = true;
  int i, imax = 0, j, k;
  double big, dum, sum, temp;

  *d=1.0;
  for (i = 1; i <= n; ++i){
    big = 0.0;
    for(j = 1; j <= n; j++) {
      temp = fabs(a[i-1][j-1]);
      if (temp > big) {
	big = temp;
      }
    }
    if (big == 0.0) {
      return false;
    }
    vv[i-1]=1.0/big;
  }
  for (j = 1; j <= n; ++j){
    for (i=1;i<j;i++) {
      sum = a[i-1][j-1];
      for(k = 1; k < i; k++) {
	sum -= a[i-1][k-1] * a[k-1][j-1];
      }
      a[i-1][j-1]=sum;
    }
    big=0.0;
    for (i = j; i <= n; i++) {
      sum = a[i-1][j-1];
      for (k = 1; k < j; k++) {
	sum -= a[i-1][k-1]*a[k-1][j-1];
      }
      a[i-1][j-1] = sum;
      dum = vv[i-1]*fabs(sum);
      if (dum >= big) {
	big = dum;
	imax = i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax-1][k-1];
	a[imax-1][k-1]=a[j-1][k-1];
	a[j-1][k-1]=dum;
      }
      *d = -(*d);
      vv[imax-1]=vv[j-1];
    }
    indx[j-1]=imax;
    if (a[j-1][j-1] == 0.0) a[j-1][j-1]=TINY;
    if (j != n) {
      dum=1.0/(a[j-1][j-1]);
      for (i=j+1;i<=n;i++) a[i-1][j-1] *= dum;
    }
  }  
  return ret;
}

/***************************************************************************/
// from Numerical Recipes in C: backward substitution
void CglRedSplit2::lubksb(double **a, int n, int *indx, double *b) const {
  int i,ii=0,ip,j;
  double sum;

  for (i = 1; i <= n; ++i){
    ip=indx[i-1];
    sum=b[ip-1];
    b[ip-1]=b[i-1];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i-1][j-1]*b[j-1];
    else if (sum) ii=i;
    b[i-1]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i-1];
    for (j=i+1;j<=n;j++) sum -= a[i-1][j-1]*b[j-1];
    b[i-1]=sum/a[i-1][i-1];
  }
}

/***************************************************************************/
double CglRedSplit2::compute_norm_change(double oldnorm, const int* list, 
					 int numElemList,
					 const double* multipliers) const {
  double newnorm = 0;
  double accumulator;
  for (int j = 0; j < nTab; ++j){
    accumulator = 0;
    for (int i = 0; i < numElemList; ++i){
      accumulator += multipliers[i]*workNonBasicTab[list[i]][j];
    }
    newnorm += accumulator*accumulator;
  }
  return (newnorm-oldnorm);
  
}
/***************************************************************************/
int CglRedSplit2::sort_rows_by_nonzeroes(struct sortElement* array, 
					 int rowIndex, int maxRows,
					 int whichTab) const {
  // Note: this function only takes into account rows which share at least
  // a nonzero column with the row that we want to reduce (rowIndex).
  int counter = 0;
  int numZeroCost = 0;
  for (int i = 0; i < mTab; ++i){
    if (!checkTime()){
      break;
    }
    if (numZeroCost == maxRows){
      // We found enough rows with cost equal to zero, i.e. "perfect" rows
      // No need to continue!
      counter = numZeroCost;
      break;
    }
    if ((i != rowIndex) && (norm[i] > param.getNormIsZero())){
      // sort rows by number of nonzeros on nonbasic columns. 
      // check if they have at least one nonzero in the same place,
      // otherwise skip
      bool match = false;
      for (int j = 0; j < nTab; ++j){
	if ((fabs(workNonBasicTab[rowIndex][j]) > param.getEPS_COEFF()) &&
	    (fabs(workNonBasicTab[i][j]) > param.getEPS_COEFF())){
	  match = true;
	  break;
	}
      }
      if (!match)
	continue;
      array[counter].index = i;
      array[counter].cost = 0;
      // whichTab = 0 means only intNonBasicTab, 1 means only
      // workNonBasicTab, 2 means both
      if (whichTab == 0 || whichTab == 2){
	for (int j = 0; j < card_intNonBasicVar; ++j){
	  if ((fabs(intNonBasicTab[rowIndex][j]) <= param.getEPS_COEFF()) &&
	      (fabs(intNonBasicTab[i][j]) > param.getEPS_COEFF())){
	    array[counter].cost += 1;
	  }
	}
      }
      if (whichTab == 1 || whichTab == 2){
	for (int j = 0; j < nTab; ++j){
	  if ((fabs(workNonBasicTab[rowIndex][j]) <= param.getEPS_COEFF()) &&
	      (fabs(workNonBasicTab[i][j]) > param.getEPS_COEFF())){
	    array[counter].cost += 1;
	  }
	}
      }
      if (array[counter].cost == 0){
	array[counter] = array[numZeroCost];
	array[numZeroCost].index = i;
	array[numZeroCost].cost = 0;
	numZeroCost++;
      }
      counter++;
    }
  }
  if (counter > maxRows){
    qsort(array, counter, sizeof(struct sortElement), rs2_compareElements);
  }
  return counter;
}
/***************************************************************************/
int CglRedSplit2::sort_rows_by_nonzeroes_greedy(struct sortElement* array,
						int rowIndex, int maxRows,
						int whichTab) const{
  int counter = 0;
  counter = sort_rows_by_nonzeroes(array, rowIndex, maxRows, whichTab);
  if (counter > maxRows){
    int i, j;
    // vector of positions of zero elements
    int* z_int = NULL;
    int* z_cont = NULL;
    // whichTab = 0 means only intNonBasicTab, 1 means only
    // workNonBasicTab, 2 means both
    if (whichTab == 0 || whichTab == 2)
      z_int = new int[card_intNonBasicVar];
    if (whichTab == 1 || whichTab == 2)
      z_cont = new int[nTab];
    // number of elements in the vectors above
    int numz_int = 0;
    int numz_cont = 0;
    // compute initial vector of nonzeroes
    if (whichTab == 0 || whichTab == 2){
      for (j = 0; j < card_intNonBasicVar; ++j){
	if (fabs(intNonBasicTab[rowIndex][j]) <= param.getEPS_COEFF()){
	  z_int[numz_int++] = j;
	}
      }
    }
    if (whichTab == 1 || whichTab == 2){
      for (j = 0; j < nTab; ++j){
	if (fabs(workNonBasicTab[rowIndex][j]) <= param.getEPS_COEFF()){
	  z_cont[numz_cont++] = j;
	}
      }
    }
    int numSorted = counter;
    counter = 1;
    // store the minimum number of nz introduced by one row
    int minNz;
    // and the position of that row
    int rowMinNz;
    // maximum number of nonzero that can be introduced
    double maxNz;
    int currentNz;
    while (counter < maxRows && counter < numSorted){
      if (!checkTime()){
	break;
      }
      // initialize minNz to a large value
      minNz = numz_int + numz_cont;;
      maxNz = array[counter].cost + array[counter-1].cost;
      rowMinNz = counter;
      for (i = counter; i < numSorted && array[i].cost < maxNz; ++i){
	int ii = array[i].index;
	// pick the row which introduces the smallest
	// number of nonzeros on the nonbasic integer columns
	currentNz = 0;
	for (j = 0; j < numz_int; ++j){
	  if ((fabs(intNonBasicTab[ii][z_int[j]]) > param.getEPS_COEFF())){
	    currentNz++;
	  }
	}
	for (j = 0; j < numz_cont; ++j){
	  if ((fabs(workNonBasicTab[ii][z_cont[j]]) > param.getEPS_COEFF())){
	    currentNz++;
	  }
	}
	array[i].cost = currentNz;
	if (currentNz < minNz){
	  rowMinNz = i;
	  minNz = currentNz;
	}
	if (currentNz == 0){
	  // perfect matching, no need to look further
	  break;
	}
      }
      // update introduced nonzeroes by the selected row;
      // first, put the row in the correct position in the array
      int swap1 = array[rowMinNz].index;
      double swap2 = array[rowMinNz].cost;
      array[rowMinNz] = array[counter];
      array[counter].index = swap1;
      array[counter].cost = swap2;
      // now count how many nonzeroes were introduces, and update
      // the vectors z_int and z_cont accordingly
      for (j = 0; j < numz_int; ++j){
	if (fabs(intNonBasicTab[swap1][z_int[j]]) > param.getEPS_COEFF()){
	  z_int[j] = z_int[numz_int-1];
	  numz_int--;
	}
      }
      for (j = 0; j < numz_cont; ++j){
	if (fabs(workNonBasicTab[swap1][z_cont[j]]) > param.getEPS_COEFF()){
	  z_cont[j] = z_cont[numz_cont-1];
	  numz_cont--;
	}
      }
      counter++;
    }
    if (z_int)
      delete[] z_int;
    if (z_cont)
      delete[] z_cont;
  }
  return counter;
}

/***************************************************************************/
int CglRedSplit2::sort_rows_by_cosine(struct sortElement* array, 
				      int rowIndex, int maxRows,
				      int whichTab) const {
  // whichTab == 0 means only integer nonbasic variables, whichTab ==
  // 1 means only workTab, whichTab == 2 means both
  int counter = 0;
  double initnorm = 0.0;
  if (whichTab == 0 || whichTab == 2){
#if RS_FAST_INT == 0
    initnorm += rs_dotProd(intNonBasicTab[rowIndex], 
			   intNonBasicTab[rowIndex],
			   card_intNonBasicVar);
#elif RS_FAST_INT == 1
    initnorm += rs_dotProd(intNonBasicTab[rowIndex], 
			   intNonBasicTab[rowIndex],
			   card_intNonBasicVar);
    double value = rs_sparseDotProd(intNonBasicTab[rowIndex],
				       intNonBasicTab[rowIndex], 
				       pi_mat[rowIndex]+mTab,
				       pi_mat[rowIndex]+mTab);
    assert (value==initnorm);
#else
    initnorm = rs_sparseDotProd(intNonBasicTab[rowIndex],
				       intNonBasicTab[rowIndex], 
				       pi_mat[rowIndex]+mTab,
				       pi_mat[rowIndex]+mTab);
#endif
  }
  if (whichTab == 1 || whichTab == 2){
    initnorm += norm[rowIndex];
  }
#if RS_FAST_WORK 
  int workOffset = mTab+card_intNonBasicVar+card_contNonBasicVar+2;
#endif
  for (int i = 0; i < mTab; ++i){
    if ((i != rowIndex) && (norm[i] > param.getNormIsZero())){
      if (!checkTime()){
	break;
      }
      array[counter].index = i;
      array[counter].cost = 0;

      if (whichTab == 0 || whichTab == 2){
#if RS_FAST_INT == 0
	array[counter].cost -= fabs(rs_dotProd(intNonBasicTab[rowIndex],
					       intNonBasicTab[i], 
					       card_intNonBasicVar));
#elif RS_FAST == 1
	array[counter].cost -= fabs(rs_dotProd(intNonBasicTab[rowIndex],
					       intNonBasicTab[i], 
					       card_intNonBasicVar));
	double value = -fabs(rs_sparseDotProd(intNonBasicTab[rowIndex],
				       intNonBasicTab[i], 
				       pi_mat[rowIndex]+mTab,
				       pi_mat[i]+mTab));
	assert (value==array[counter].cost);
#else
	array[counter].cost = -fabs(rs_sparseDotProd(intNonBasicTab[rowIndex],
						     intNonBasicTab[i], 
						     pi_mat[rowIndex]+mTab,
						     pi_mat[i]+mTab));
#endif
      }
      if (whichTab == 1 || whichTab == 2){
#if RS_FAST_WORK == 0
	  array[counter].cost -= fabs(rs_dotProd(workNonBasicTab[rowIndex],
						 workNonBasicTab[i], nTab));
#elif RS_FAST_WORK == 1
	  double value = array[counter].cost;
	  array[counter].cost -= fabs(rs_dotProd(workNonBasicTab[rowIndex],
						 workNonBasicTab[i], nTab));
	  value -= fabs(rs_sparseDotProd(workNonBasicTab[rowIndex],
				       workNonBasicTab[i], 
				       pi_mat[rowIndex]+workOffset,
				       pi_mat[i]+workOffset));
	  assert (value==array[counter].cost);
#else
	  array[counter].cost -= fabs(rs_sparseDotProd(workNonBasicTab[rowIndex],
				       workNonBasicTab[i], 
				       pi_mat[rowIndex]+workOffset,
				       pi_mat[i]+workOffset));
#endif
      }	

      // Now divide by the two norms
      double denom = 0.0;
      if (whichTab == 0 || whichTab == 2){
#if RS_FAST_INT == 0
	denom += initnorm*rs_dotProd(intNonBasicTab[i], intNonBasicTab[i],
				     card_intNonBasicVar);
#elif RS_FAST == 1
	denom += initnorm*rs_dotProd(intNonBasicTab[i], intNonBasicTab[i],
				     card_intNonBasicVar);
	double value = initnorm*rs_sparseDotProd(intNonBasicTab[i],
				       intNonBasicTab[i], 
				       pi_mat[i]+mTab,
				       pi_mat[i]+mTab);
	assert (value==denom);
#else
	denom = initnorm*rs_sparseDotProd(intNonBasicTab[i],
					  intNonBasicTab[i], 
					  pi_mat[i]+mTab,
					  pi_mat[i]+mTab);
#endif
      }
      if (whichTab == 1 || whichTab == 2){
	denom += initnorm*norm[i];
      }
      array[counter].cost /= sqrt(denom);
      if (array[counter].cost != 0.0){
	counter++;
      }
    }
  }
  if (counter >= maxRows){
    qsort(array, counter, sizeof(struct sortElement), rs2_compareElements);
  }
  return counter;
}

/***************************************************************************/
int CglRedSplit2::get_list_rows_reduction(int rowIndex, int maxRowsReduction, 
					  int* list, const double* norm,
					  CglRedSplit2Param::RowSelectionStrategy selectionStrategy) const{
#ifdef RS2_TRACE
  printf("Obtaining list of rows with column strategy %d\n", selectionStrategy);
#endif
  struct sortElement* array = new struct sortElement[mTab];
  int counter = 0;
  int i, j;
  // Look in CglRedSplit2Param for a description of each strategy
  if (selectionStrategy == CglRedSplit2Param::RS1){
    counter = sort_rows_by_nonzeroes(array, rowIndex, maxRowsReduction-1, 0);
  } 
  else if (selectionStrategy == CglRedSplit2Param::RS2){
    counter = sort_rows_by_nonzeroes(array, rowIndex, maxRowsReduction-1, 1);
  } 
  else if (selectionStrategy == CglRedSplit2Param::RS3){
    counter = sort_rows_by_nonzeroes(array, rowIndex, maxRowsReduction-1, 2);
  } 
  else if (selectionStrategy == CglRedSplit2Param::RS4){
    // compute the *true* number of nonzeroes that we introduce
    // on the integer nonbasics
    counter = sort_rows_by_nonzeroes_greedy(array, rowIndex, maxRowsReduction-1, 0);
  }
  else if (selectionStrategy == CglRedSplit2Param::RS5){
    // compute the *true* number of nonzeroes that we introduce
    // on the working set (continuous nonbasics)
    counter = sort_rows_by_nonzeroes_greedy(array, rowIndex, maxRowsReduction-1, 1);
  }
  else if (selectionStrategy == CglRedSplit2Param::RS6){
    // compute the *true* number of nonzeroes that we introduce
    // on all the nonbasics: intNonBasicVar and workNonBasicVar
    counter = sort_rows_by_nonzeroes_greedy(array, rowIndex, maxRowsReduction-1, 2);
  }
  else if (selectionStrategy == CglRedSplit2Param::RS7){
    // sort rows by cosine of the angle in the space of the nonbasic integer
    // and nonbasic continuous (working set) columns, wrt row rowIndex
    counter = sort_rows_by_cosine(array, rowIndex, maxRowsReduction-1, 2);
  }
  else if (selectionStrategy == CglRedSplit2Param::RS8){
    // sort rows by cosine of the angle in the space of the
    // nonbasic continuous (working set) columns, wrt row rowIndex
    counter = sort_rows_by_cosine(array, rowIndex, maxRowsReduction-1, 1);
  }
  list[0] = rowIndex;
  j = 1;
  for (i = 0; i < counter && j < maxRowsReduction; ++i){
    list[j] = array[i].index;
    j++;
  }
  delete[] array;
  return j;
}
/***************************************************************************/
void CglRedSplit2::fill_workNonBasicTab(CglRedSplit2Param::ColumnSelectionStrategy strategy, const int* ignore_list){
  // Sort continuous non basic variables by reduced cost and choose the
  // columns we will work with.
  // We can choose among several strategies; these are described in 
  // CglRedSplit2Param.hpp.
  int i, j;
  if (strategy == CglRedSplit2Param::CS_ALLCONT){
    // select all continuous nonbasic columns
    for (i = 0; i < mTab; ++i){
      memcpy(workNonBasicTab[i], contNonBasicTab[i], 
	     card_contNonBasicVar*sizeof(double));
    }
    nTab = card_contNonBasicVar;
  }
  else{
    // Begin by sorting all continuous nonbasic columns by increasing
    // reduced cost (except those in the ingore_list)
    struct sortElement* array = new struct sortElement[card_contNonBasicVar];
    int pos = 0, iter = 0;
    for (i = 0; i < card_contNonBasicVar; ++i){
      if (ignore_list != NULL){
	iter = 0;
	while (ignore_list[iter] >= 0 && 
	       ignore_list[iter] != contNonBasicVar[i]){
	  iter++;
	}
	if (ignore_list[iter] == contNonBasicVar[i]){
	  continue;
	}
      }
      array[pos].index = i;
      // We take absolute value of reduced costs because sign could be
      // positive or negative depending on minimization/maximization,
      // and which bound the variable is at (upper/lower).
      if (contNonBasicVar[i] < ncol)
	array[pos].cost = fabs(reducedCost[contNonBasicVar[i]]);
      else
	array[pos].cost = fabs(rowPrice[contNonBasicVar[i]-ncol]);
      pos++;
    }
    qsort(array,pos,sizeof(struct sortElement),rs2_compareElements);
    nTab = 0;
    int card_contNonBasicVar = pos;
    // Now choose the variables. These strategies are described in
    // CglRedSplit2Param.hpp.  
    // CS1, CS2 and CS3 are the 3 partitions of C-3P: each one
    // contains one third of the variables, sorted by reduced costs
    if (strategy == CglRedSplit2Param::CS1){
      for (i = 0; i < card_contNonBasicVar/3; ++i){
	for (j = 0; j < mTab; ++j){
	  workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	}
	nTab++;
      }
    }
    else if (strategy == CglRedSplit2Param::CS2){
      for (i = card_contNonBasicVar/3; i < 2*card_contNonBasicVar/3; ++i){
	if (!checkTime()){
	  break;
	}
	for (j = 0; j < mTab; ++j){
	  workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	}
	nTab++;
      }
    }
    else if (strategy == CglRedSplit2Param::CS3){
      for (i = 2*card_contNonBasicVar/3; i < card_contNonBasicVar; ++i){
	if (!checkTime()){
	  break;
	}
	for (j = 0; j < mTab; ++j){
	  workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	}
	nTab++;
      }
    }
    // CS4, CS5, CS6, CS7 and CS8 are the 5 partitions of C-5P; each
    // one contains 1/5 of the variables, sorted by increasing reduced
    // cost: first set has smallest reduced costs, and so on.
    else if (strategy == CglRedSplit2Param::CS4){
      for (i = 0; i < card_contNonBasicVar/5; ++i){
	if (!checkTime()){
	  break;
	}
	for (j = 0; j < mTab; ++j){
	  workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	}
	nTab++;	
      }
    }
    else if (strategy == CglRedSplit2Param::CS5){
      for (i = card_contNonBasicVar/5; i < 2*card_contNonBasicVar/5; ++i){
	if (!checkTime()){
	  break;
	}
	for (j = 0; j < mTab; ++j){
	  workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	}
	nTab++;	
      }
    }
    else if (strategy == CglRedSplit2Param::CS6){
      for (i = 2*card_contNonBasicVar/5; i < 3*card_contNonBasicVar/5; ++i){
	if (!checkTime()){
	  break;
	}
	for (j = 0; j < mTab; ++j){
	  workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	}
	nTab++;	
      }
    }
    else if (strategy == CglRedSplit2Param::CS7){
      for (i = 3*card_contNonBasicVar/5; i < 4*card_contNonBasicVar/5; ++i){
	if (!checkTime()){
	  break;
	}
	for (j = 0; j < mTab; ++j){
	  workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	}
	nTab++;	
      }
    }
    else if (strategy == CglRedSplit2Param::CS8){
      for (i = 4*card_contNonBasicVar/5; i < card_contNonBasicVar; ++i){
	if (!checkTime()){
	  break;
	}
	for (j = 0; j < mTab; ++j){
	  workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	}
	nTab++;	
      }
    }
    // CS9 and CS10 are the 2 partitions of I-2P-2/3; this is a
    // partition with interleaving pattern XX--X- on the 2/3 of the
    // variables with smallest reduced cost
    else if (strategy == CglRedSplit2Param::CS9){
      int pat;
      for (i = 0; i < 2*card_contNonBasicVar/3; ++i){
	if (!checkTime()){
	  break;
	}
	pat = i%6;
	if (pat == 0 || pat == 1 || pat == 4){
	  for (j = 0; j < mTab; ++j){
	    workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	  }
	  nTab++;
	}
      }
    }
    else if (strategy == CglRedSplit2Param::CS10){
      int pat;
      for (i = 0; i < 2*card_contNonBasicVar/3; ++i){
	if (!checkTime()){
	  break;
	}
	pat = i%6;
	if (pat == 2 || pat == 3 || pat == 5){
	  for (j = 0; j < mTab; ++j){
	    workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	  }
	  nTab++;	
	}
      }
    }
    // CS11 and CS12 are the 2 partitions of I-2P-4/5; partition with
    // interleaving pattern X---XXX- on 4/5 of the vars with smallest
    // reduced cost
    else if (strategy == CglRedSplit2Param::CS11){
      int pat;
      for (i = 0; i < 4*card_contNonBasicVar/5; ++i){
	if (!checkTime()){
	  break;
	}
	pat = i%8;
	if (pat == 0 || pat == 4 || pat == 5 || pat == 6){
	  for (j = 0; j < mTab; ++j){
	    workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	  }
	  nTab++;
	}
      }
    }
    else if (strategy == CglRedSplit2Param::CS12){
      int pat;
      for (i = 0; i < 4*card_contNonBasicVar/5; ++i){
	if (!checkTime()){
	  break;
	}
	pat = i%8;
	if (pat == 1 || pat == 2 || pat == 3 || pat == 7){
	  for (j = 0; j < mTab; ++j){
	    workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	  }
	  nTab++;
	}
      }
    }
    // CS13, and CS14 are the 2 partitions of I-2P-1/2; interleaving
    // partition with pattern X--X on 1/2 of the variables
    else if (strategy == CglRedSplit2Param::CS13){
      for (i = 0; i < card_contNonBasicVar/2; ++i){
	if (!checkTime()){
	  break;
	}
	if (i%4 == 0 || i%4 == 3){
	  for (j = 0; j < mTab; ++j){
	    workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	  }
	  nTab++;	
	}
      }
    }
    else if (strategy == CglRedSplit2Param::CS14){
      for (i = 0; i < card_contNonBasicVar/2; ++i){
	if (!checkTime()){
	  break;
	}
	if (i%4 == 1 || i%4 == 2){
	  for (j = 0; j < mTab; ++j){
	    workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	  }
	  nTab++;
	}
      }
    }
    // CS15, CS16 and CS17 are the 3 partitions of I-3P; interleaving
    // partition with pattern X-/ on all the variables
    else if (strategy == CglRedSplit2Param::CS15){
      for (i = 0; i < card_contNonBasicVar; ++i){
	if (!checkTime()){
	  break;
	}
	if (i%3 == 0){
	  for (j = 0; j < mTab; ++j){
	    workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	  }
	  nTab++;
	}
      }
    }
    else if (strategy == CglRedSplit2Param::CS16){
      for (i = 0; i < card_contNonBasicVar; ++i){
	if (!checkTime()){
	  break;
	}
	if (i%3 == 1){
	  for (j = 0; j < mTab; ++j){
	    workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	  }
	  nTab++;
	}
      }
    }
    else if (strategy == CglRedSplit2Param::CS17){
      for (i = 0; i < card_contNonBasicVar; ++i){
	if (!checkTime()){
	  break;
	}
	if (i%3 == 2){
	  for (j = 0; j < mTab; ++j){
	    workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	  }
	  nTab++;
	}
      }
    }
    // CS18, CS19, CS20, CS21 are the 4 partitions of I-4P;
    // interleaving partition with pattern X-X/\\-/ on all the
    // variables
    else if (strategy == CglRedSplit2Param::CS18){
      int pat;
      for (i = 0; i < card_contNonBasicVar; ++i){
	if (!checkTime()){
	  break;
	}
	pat = i%8;
	if (pat == 0 || pat == 2){
	  for (j = 0; j < mTab; ++j){
	    workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	  }
	  nTab++;
	}
      }
    }
    else if (strategy == CglRedSplit2Param::CS19){
      int pat;
      for (i = 0; i < card_contNonBasicVar; ++i){
	if (!checkTime()){
	  break;
	}
	pat = i%8;
	if (pat == 1 || pat == 6){
	  for (j = 0; j < mTab; ++j){
	    workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	  }
	  nTab++;
	}
      }
    }
    else if (strategy == CglRedSplit2Param::CS20){
      int pat;
      for (i = 0; i < card_contNonBasicVar; ++i){
	if (!checkTime()){
	  break;
	}
	pat = i%8;
	if (pat == 3 || pat == 7){
	  for (j = 0; j < mTab; ++j){
	    workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	  }
	  nTab++;
	}
      }
    }
    else if (strategy == CglRedSplit2Param::CS21){
      int pat;
      for (i = 0; i < card_contNonBasicVar; ++i){
	if (!checkTime()){
	  break;
	}
	pat = i%8;
	if (pat == 4 || pat == 5){
	  for (j = 0; j < mTab; ++j){
	    workNonBasicTab[j][nTab] = contNonBasicTab[j][array[i].index];
	  }
	  nTab++;
	}	
      }
    }
    delete[] array;
  }
  //  create sparse work stuff
  for(i=0; i<mTab; i++) {
#if RS_FAST_WORK == 0
    norm[i] = rs_dotProd(workNonBasicTab[i], workNonBasicTab[i], nTab);
#elif RS_FAST_WORK == 1
    int * info = pi_mat[i];
    int * workSparse=info+mTab+2+card_intNonBasicVar+card_contNonBasicVar;
    int n = 0;
    const double * values = workNonBasicTab[i];
    double checkNorm=0.0;
    for (int j=0;j<nTab;j++) {
      if (values[j]) {
	workSparse[++n]=j;
	checkNorm += values[j]*values[j];
      }
    }
    workSparse[0]=n; 
    norm[i] = rs_dotProd(workNonBasicTab[i], workNonBasicTab[i], nTab);
    assert (checkNorm==norm[i]);
#else
    int * info = pi_mat[i];
    int * workSparse=info+mTab+2+card_intNonBasicVar+card_contNonBasicVar;
    int n = 0;
    const double * values = workNonBasicTab[i];
    double checkNorm=0.0;
    for (int j=0;j<nTab;j++) {
      if (values[j]) {
	workSparse[++n]=j;
	checkNorm += values[j]*values[j];
      }
    }
    workSparse[0]=n; 
    norm[i] = checkNorm;
#endif
  }
}
/***************************************************************************/
void CglRedSplit2::reduce_workNonBasicTab(int numRowsReduction, 
					  CglRedSplit2Param::RowSelectionStrategy rowSelectionStrategy,
					  int maxIterations) {
  // Use at most this number of rows
  int maxRowsReduction = std::min(numRowsReduction, mTab);
  if (maxRowsReduction == 1){
    return;
  }
  int i, j, k, h;

  // Allocate space to store the matrix for the linear system
#ifdef RS2_USE_LAPACK
  double* A = new double[maxRowsReduction*maxRowsReduction];
#else
  double** A;
  rs_allocmatDBL(&A, maxRowsReduction, maxRowsReduction);
#endif
  // Right hand side
  double* b = new double[maxRowsReduction];
  // Data for LU decomposition
  int* indexlu = new int[maxRowsReduction];
  double tmpnumlu = 0.0;
  double* tmpveclu = new double[maxRowsReduction];
  // List of rows involved in the linear combination
  int* list = new int[maxRowsReduction];
  
  // Number of rows actually used
  int numUsedRows;
  
  bool resolveWithNormalization = false;

  for (k = 0; k < mTab && k < maxIterations; ++k) {
    if (!checkTime()){
      break;
    }
    if (norm[k] > param.getNormIsZero()) {
      // Obtain the list of rows that should be used to reduce row k
      numUsedRows = get_list_rows_reduction(k, maxRowsReduction, list, norm,
					    rowSelectionStrategy);
#ifdef RS2_TRACE
      rs_printvecINT("rows used for reduction", list, numUsedRows);
#endif
      if (numUsedRows <= 1){
	// This means that the list only contains the current row;
	// Thus, no rows were selected for reduction. Skip.
	continue;
      }      
      // Note: the list must have size maxRowsReduction.
      // Now prepare the linear system according to the paper
      for (i = 0; i < numUsedRows; ++i){
	for (j = 0; j < numUsedRows; ++j){
#ifdef RS2_USE_LAPACK
	  A[i*numUsedRows+j] = 0;
#else
	  A[i][j] = 0;
#endif
	  if (list[i] != k && list[j] != k){
	    for (h = 0; h < nTab; ++h){
#ifdef RS2_USE_LAPACK
	      A[i*numUsedRows+j]+=workNonBasicTab[list[i]][h]*workNonBasicTab[list[j]][h];
#else
	      A[i][j]+=workNonBasicTab[list[i]][h]*workNonBasicTab[list[j]][h];
#endif
	    }
	    if (resolveWithNormalization && i == j){
	      // Penalize the norm of lambda, i.e. the solution
#ifdef RS2_USE_LAPACK
	      A[i*numUsedRows+j] += norm[k]*param.getNormalization();
#else
	      A[i][j] += norm[k]*param.getNormalization();
#endif
	    }
	  }
	}
	if (list[i] == k){
	  b[i] = 1;
#ifdef RS2_USE_LAPACK
	  A[i*numUsedRows+i] = 1;
#else
	  A[i][i] = 1;
#endif
	}
	else{
	  b[i] = 0;
	  for (h = 0; h < nTab; ++h){
	    b[i] -= workNonBasicTab[list[i]][h]*workNonBasicTab[k][h];
	  }
	}
      }
      // Linear system has been written, now solve it
#ifdef RS2_USE_LAPACK
      int info = 0;
      int h = 1;
      dgesv_(&numUsedRows, &h, A, &numUsedRows, indexlu, 
	     b, &numUsedRows, &info);
	
      if (info)
	continue;
#else
      // LU decomposition
      if (!ludcmp(A, numUsedRows, indexlu, &tmpnumlu, tmpveclu)){
	// numerical error: exit
	continue;
      }
      // Backward substitution
      lubksb(A, numUsedRows, indexlu, b);
#endif
      // Check the 1-norm of the solution
      double sumnorm = 0.0;
      for (i = 0; i < numUsedRows; ++i){
	b[i] = rs2round(b[i]);
	sumnorm += fabs(b[i]);
	if (sumnorm > param.getMaxSumMultipliers()){
	  break;
	}
      }
#ifdef RS2_TRACE
      rs_printvecDBL("multipliers", b, numUsedRows);
      printf("sumnorm: %f\n", sumnorm);
#endif
      if (sumnorm == 1){
	// The optimal solution has lambda_k = 1 and all the rest is zero;
	// so we cannot reduce any coefficient, let's skip this
	continue;
      }
      else if (!resolveWithNormalization && 
	       sumnorm > param.getMaxSumMultipliers()){
	// The norm of lambda is too large, resolve with normalization
	// (note that we do not want to do this more than once)
	resolveWithNormalization = true;
	k--;
	continue;
      }
      resolveWithNormalization = false;
      if (sumnorm > param.getMaxSumMultipliers()){
	// If we got this far, even resolving did not help, so we skip
	continue;
      }
      double deltaNorm = compute_norm_change(norm[k], list, 
					     numUsedRows, b);
      if (deltaNorm <= -norm[k]*param.getMinNormReduction()){
	for (i = 0; i < numUsedRows; ++i){
	  pi_mat[k][list[i]] = (int)(b[i]);
	}
	numRedRows++;
      }
    } /* if (norm[k] > param.getNormIsZero()) */
  } /*for (k = 0; k < mTab; ++k) */
  delete[] b;
  delete[] list;
  delete[] indexlu;
  delete[] tmpveclu;
#ifdef RS2_USE_LAPACK
  delete[] A;
#else
  rs_deallocmatDBL(&A, maxRowsReduction);
#endif

#ifdef RS2_TRACE
  double sum_norms = 0;
  for(i=0; i<mTab; i++) {
    sum_norms += norm[i];
  }
  printf("CglRedSplit2::reduce_contNonBasicTab():Final sum of norms: %f\n", sum_norms);
  printf("CglRedSplit2::reduce_contNonBasicTab(): %d rows reduced\n",numRedRows);
#endif

} /* reduce_workNonBasicTab */

/************************************************************************/
void CglRedSplit2::generate_row(int index_row, double *row) {

  memset(row, 0, (ncol+nrow)*sizeof(double));
  // we only deal with nonbasic variables - the cut coefficient will be zero
  // on the basic ones anyway
#if RS_FAST_INT == 0 && RS_FAST_CONT == 0 && RS_FAST_WORK == 0
  int i;
  int locind, j;
  for(i=0; i<card_intNonBasicVar; i++) {
    locind = intNonBasicVar[i];
    row[locind] = 0;
    for(j=0; j<mTab; j++) {
      row[locind] += pi_mat[index_row][j] * intNonBasicTab[j][i];
    }
  }
  for(i=0; i<card_contNonBasicVar; i++) {
    locind = contNonBasicVar[i];
    row[locind] = 0;
    for(j=0; j<mTab; j++) {
      row[locind] += pi_mat[index_row][j] * contNonBasicTab[j][i];
    }
  }
#elif RS_FAST_INT < 2 || RS_FAST_CONT < 2
  double * rowTemp = new double[ncol+nrow];
  memset(rowTemp, 0, (ncol+nrow)*sizeof(double));
  for(int i=0; i<card_intNonBasicVar; i++) {
    int locind = intNonBasicVar[i];
    rowTemp[locind] = 0;
    for(int j=0; j<mTab; j++) {
      rowTemp[locind] += pi_mat[index_row][j] * intNonBasicTab[j][i];
    }
  }
  for(int i=0; i<card_contNonBasicVar; i++) {
    int locind = contNonBasicVar[i];
    rowTemp[locind] = 0;
    for(int j=0; j<mTab; j++) {
      rowTemp[locind] += pi_mat[index_row][j] * contNonBasicTab[j][i];
    }
  }
  const int * pi = pi_mat[index_row];
#if RS_FAST_CONT 
  int contOffset = mTab+card_intNonBasicVar+1;
#endif
  for(int j=0; j<mTab; j++) {
    double value = pi[j];
    if (value) {
      const double * tableau = intNonBasicTab[j];
#if RS_FAST_INT < 1
      for(int i=0; i<card_intNonBasicVar; i++) {
	int locind = intNonBasicVar[i];
	row[locind] += value * tableau[i];
      }
#else
      const int * which = pi_mat[j]+mTab;
      int n = which[0];
      which++;
      for (int k=0;k<n;k++) {
	int i=which[k];
	int locind = intNonBasicVar[i];
	row[locind] += value * tableau[i];
      }
#endif
      tableau = contNonBasicTab[j];
#if RS_FAST_CONT < 1
      for(int i=0; i<card_contNonBasicVar; i++) {
	int locind = contNonBasicVar[i];
	row[locind] += value * tableau[i];
      }
#else
      const int * which2 = pi_mat[j]+contOffset;
      int n2 = which2[0];
      which2++;
      for (int k=0;k<n2;k++) {
	int i=which2[k];
	int locind = contNonBasicVar[i];
	row[locind] += value * tableau[i];
      }
#endif
    }
  }
  for (int i=0;i<nrow+ncol;i++) {
    assert (fabs(row[i]-rowTemp[i])<1.0e-11+1.0e-8*std::max(fabs(row[i]),fabs(rowTemp[i])));
  }
  delete [] rowTemp;
#else
  const int * pi = pi_mat[index_row];
  int contOffset = mTab+card_intNonBasicVar+1;
  for(int j=0; j<mTab; j++) {
    double value = pi[j];
    if (value) {
      const double * tableau = intNonBasicTab[j];
      const int * which = pi_mat[j]+mTab;
      int n = which[0];
      which++;
      for (int k=0;k<n;k++) {
	int i=which[k];
	int locind = intNonBasicVar[i];
	row[locind] += value * tableau[i];
      }
      tableau = contNonBasicTab[j];
      which = pi_mat[j]+contOffset;
      n = which[0];
      which++;
      for (int k=0;k<n;k++) {
	int i=which[k];
	int locind = contNonBasicVar[i];
	row[locind] += value * tableau[i];
      }
    }
  }
#endif
} /* generate_row */

/************************************************************************/
int CglRedSplit2::generate_cgcut(double *row, double *rhs) {
  
#ifdef RS2_TRACE
  rs_printvecDBL("CglRedSplit2::generate_cgcut(): starting row", 
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

  int i, locind;
  double f;

  for(i=0; i<card_intNonBasicVar; i++) {
    locind = intNonBasicVar[i];
    f = rs_above_integer(row[locind]);

    if(f > f0) {
      row[locind] = - ((1-f) * f0);
    }
    else {
      row[locind] = -(f*f0compl);
    }
  }

  for(i=0; i<card_contNonBasicVar; i++) {
    locind = contNonBasicVar[i];
    if(row[locind] < 0) {
      row[locind] *= f0;
    }
    else {
      row[locind] = -(row[locind]*f0compl);
    }
  }
  (*rhs) = -f0*f0compl;
  
#ifdef RS2_TRACE
  rs_printvecDBL("CglRedSplit2::generate_cgcut(): row", row, ncol+nrow);
  printf("rhs: %f\n", *rhs);
#endif

  return(1);
} /* generate_cgcut */

/************************************************************************/
void CglRedSplit2::eliminate_slacks(double *row, 
				   const double *elements, 
				   const CoinBigIndex *rowStart,
				   const int *indices,
				   const int *rowLength,
				   const double *rhs, double *tabrowrhs) {

  for(int i=0; i<nrow; i++) {
    if(fabs(row[ncol+i]) > param.getEPS_ELIM()) {

      CoinBigIndex upto = rowStart[i] + rowLength[i];
      for(CoinBigIndex j=rowStart[i]; j<upto; j++) {
	row[indices[j]] -= row[ncol+i] * elements[j];      
      }
      *tabrowrhs -= row[ncol+i] * rhs[i];
    }
  }

} /* eliminate_slacks */

/************************************************************************/
void CglRedSplit2::flip(double *row) {
  
  int i;
  for(i=0; i<card_nonBasicAtUpper; i++) {
    row[nonBasicAtUpper[i]] = -row[nonBasicAtUpper[i]];
  }
} /* flip */

/************************************************************************/
void CglRedSplit2::unflip(double *row, double *tabrowrhs) {
  
  int i;
  for(i=0; i<card_nonBasicAtLower; i++) {
    int locind = nonBasicAtLower[i];
    if(locind < ncol) {
      *tabrowrhs += row[locind] * colLower[locind];
    }
  }
  for(i=0; i<card_nonBasicAtUpper; i++) {
    int locind = nonBasicAtUpper[i];
    row[locind] = -row[locind];
    if(locind < ncol) {
      *tabrowrhs += row[locind] * colUpper[locind];
    }
  }

} /* unflip */

/************************************************************************/
int CglRedSplit2::check_dynamism(double *row) {

  int i;
  double val, max_val = 0, min_val = param.getINFINIT();

  for(i=0; i<ncol; i++) {
    val = fabs(row[i]);
    max_val = std::max(max_val, val);
    if(val > param.getEPS_COEFF()) {
      min_val = std::min(min_val, val);
    }
  }

  if((max_val < param.getMAXDYN() * min_val) && 
     (max_val >= min_val)) {
    return 1;
  }

#ifdef RS2_TRACE
  printf("CglRedSplit2::check_dynamism(): max_val: %6.6f   min_val: %6.6f   dynamism: %g\n",
	 max_val, min_val, max_val/min_val);
#endif

  return 0;
} /* row_scale_factor */

/************************************************************************/
int CglRedSplit2::generate_packed_row(const double *lclXlp,
				      double *row,
				      int *rowind, double *rowelem, 
				      int *card_row, double & rhs) {
  int i;
  double value;
  int max_support = param.getMAX_SUPP_ABS() + 
    (static_cast<int>(ncol*param.getMAX_SUPP_REL()));

  if(!check_dynamism(row)) {
#ifdef RS2_TRACE
    printf("CglRedSplit2::generate_packed_row(): Cut discarded (bad dynamism)\n");
#endif	
    return(0);
  }

  *card_row = 0;
  
  for(i=0; i<ncol; i++) {
    value = row[i];
    if(fabs(value) > param.getEPS_COEFF()) {
      rowind[*card_row] = i;
      rowelem[*card_row] = value;
      (*card_row)++;
      if(*card_row > max_support) {
#ifdef RS2_TRACE
	printf("CglRedSplit2::generate_packed_row(): Cut discarded (too many non zero)\n");
#endif	
	return(0);
      }
    } else {
      if (value > 0.0) {
        rhs -= value * std::max(colLower[i],-1.0e20);
      } 
      else {
        rhs -= value * std::min(colUpper[i],1.0e20);      
      } 
    }
  }
  value = 0;
  for(i=0; i<(*card_row); i++) {
    value += lclXlp[rowind[i]] * rowelem[i];
  }

  if(value > rhs) {
    value = value-rhs;
    if(value < param.getMINVIOL()*(*card_row)) {

#ifdef RS2_TRACE
      printf("CglRedSplit2::generate_packed_row(): Cut discarded: violation: %12.10f\n", value);
#endif	
      
      return(0);
    }
  }
  
  return(1);
} /* generate_packed_row */

/************************************************************************/
bool CglRedSplit2::rs_are_different_vectors(const int *vect1, 
					    const int *vect2,
					    const int dim) {
  int i;
  for(i=0; i<dim; i++) {
    if(vect1[i] != vect2[i]) {
      return(1);
    }    
  }
  return(0);
} /* rs_are_different_vectors */

/************************************************************************/
void CglRedSplit2::generateCuts(const OsiSolverInterface &si, OsiCuts & cs,
				const CglTreeInfo info)
{
  solver = const_cast<OsiSolverInterface *>(&si);
  if(solver == NULL) {
    printf("### WARNING: CglRedSplit2::generateCuts(): no solver available.\n");
    return;    
  }  

  if(!solver->optimalBasisIsAvailable()) {
    printf("### WARNING: CglRedSplit2::generateCuts(): no optimal basis available.\n");
    return;
  }
  solver->getColType(true);
  // Reset some members of CglRedSplit2
  card_intBasicVar = 0;
  card_intBasicVar_frac = 0;
  card_intNonBasicVar = 0;
  card_contNonBasicVar = 0;
  card_nonBasicAtUpper = 0;
  card_nonBasicAtLower = 0;
  numRedRows = 0;
  startTime = CoinCpuTime();

  // Get basic problem information from solver
  ncol = solver->getNumCols(); 
  nrow = solver->getNumRows(); 
  colLower = solver->getColLower();
  colUpper = solver->getColUpper();
  rowLower = solver->getRowLower();
  rowUpper = solver->getRowUpper();
  rowRhs = solver->getRightHandSide();
  reducedCost = solver->getReducedCost();
  rowPrice = solver->getRowPrice();
  objective = solver->getObjCoefficients();

  xlp = solver->getColSolution();
  rowActivity = solver->getRowActivity();
  byRow = solver->getMatrixByRow();

  solver->enableFactorization();
  generateCuts(&cs, param.getMaxNumCuts());
  solver->disableFactorization();
} /* generateCuts */

/************************************************************************/
int CglRedSplit2::generateCuts(OsiCuts* cs, int maxNumCuts, int* lambda)
{
  int i;
  is_integer = new int[ncol]; 
  
  compute_is_integer();

  int *cstat = new int[ncol];
  int *rstat = new int[nrow];


  solver->getBasisStatus(cstat, rstat);   // 0: free  1: basic  
                                          // 2: upper 3: lower

  int *basis_index = new int[nrow]; // basis_index[i] = 
                                    //        index of pivot var in row i
                                    //        (slack if number >= ncol) 

#ifdef RS2_TRACETAB
  rs_printvecINT("cstat", cstat, ncol);
  rs_printvecINT("rstat", rstat, nrow);
#endif
  solver->getBasics(basis_index);

  cv_intBasicVar = new int[ncol];  
  cv_intBasicVar_frac = new int[ncol];  
  intBasicVar = new int[ncol];                                 
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
    cv_intBasicVar[i] = 0;
    cv_intBasicVar_frac[i] = 0;

    switch(cstat[i]) {
    case 1: // basic variable
      
      dist_int = rs_above_integer(xlp[i]);
      if(is_integer[i]){
	if ((dist_int > param.getAway()) && (dist_int < 1 - param.getAway())) {
	  cv_intBasicVar_frac[i] = 1;
	  card_intBasicVar_frac++;

	  // intBasicVar_frac computed below, 
	  // as order must be according to selected rows
	}
	card_intBasicVar++;
	cv_intBasicVar[i] = 1;

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
      printf("### ERROR: CglRedSplit2::generateCuts(): cstat[%d]: %d\n",
	     i, cstat[i]);
#endif
      goodModel = false;
      break;
    } 
  }

  // Use this instead of rowRhs to allow for ranges
  double *effective_rhs = new double[nrow];

  for(i=0; i<nrow; i++) {
    // effective rhs
    effective_rhs[i] = rowRhs[i];
    switch(rstat[i]) {
    case 1: // basic slack
      break;

    case 2: // non basic slack at upper; flipped and shifted
      effective_rhs[i] = rowLower[i];
      nonBasicAtUpper[card_nonBasicAtUpper] = ncol+i;
      card_nonBasicAtUpper++;

      contNonBasicVar[card_contNonBasicVar] = ncol+i;
      card_contNonBasicVar++;
      break;

    case 3: // non basic slack at lower; shifted
      effective_rhs[i] = rowUpper[i];
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
    assert (fabs(effective_rhs[i])<1.0e100);
  }

#ifdef RS2_TRACE
  printf("CglRedSplit2()::card_intBasicVar_frac %d %d\n", 
	 card_intBasicVar_frac, card_contNonBasicVar);
#endif
  if((card_contNonBasicVar == 0) || (card_intBasicVar_frac == 0)
     || !goodModel) {
    delete[] cstat;
    delete[] rstat;
    delete[] basis_index;

    delete[] cv_intBasicVar;  
    delete[] cv_intBasicVar_frac;  
    delete[] intBasicVar;
    delete[] intNonBasicVar;
    delete[] contNonBasicVar;
    delete[] nonBasicAtUpper;
    delete[] nonBasicAtLower;
    delete[] is_integer;
    delete[] effective_rhs;

    return 0; // no cuts can be generated
  }

  double *z = new double[ncol];  // workspace to get row of the tableau
  double *slack = new double[nrow];  // workspace to get row of the tableau

#ifdef RS2_TRACETAB
  printOptTab(solver);
#endif

  mTab = card_intBasicVar;
  nTab = card_contNonBasicVar;

  rhsTab = new double[mTab];
  cv_fracRowsTab = new int[mTab];
  memset(cv_fracRowsTab, 0, mTab*sizeof(int));
  int card_rowTab = 0;

  rs_allocmatDBL(&contNonBasicTab, mTab, card_contNonBasicVar);
  rs_allocmatDBL(&workNonBasicTab, mTab, card_contNonBasicVar);
  rs_allocmatDBL(&intNonBasicTab, mTab, card_intNonBasicVar);
  norm = new double[mTab];

  intBasicVar_frac = new int[ncol];                                 

  // position of each integer basic variable in the simplex tableau
  int* origRow = new int[nrow];

  card_intBasicVar = 0; // recompute in pivot order
  card_intBasicVar_frac = 0;

  for(i=0; i<nrow; i++) {

    int ind_row = i;

    if(basis_index[ind_row] >= ncol || 
       cv_intBasicVar[basis_index[ind_row]] != 1) {
      // If basic variable is not an integer variable, skip
      continue;
    } 

    origRow[card_intBasicVar] = ind_row;

    // row used in generation
    intBasicVar[card_intBasicVar] = basis_index[ind_row];
    if (cv_intBasicVar_frac[basis_index[ind_row]] == 1){
      // row is fractional
      intBasicVar_frac[card_intBasicVar_frac] = basis_index[ind_row];
      cv_fracRowsTab[card_intBasicVar] = 1;
    }
    // obtain row of simplex tableau
    solver->getBInvARow(ind_row, z, slack);

    rhsTab[card_rowTab] = xlp[basis_index[ind_row]];

    int nonzeroes = 0;
    int ii;

    for(ii=0; ii<card_contNonBasicVar; ii++) {
      int locind = contNonBasicVar[ii];
      if(locind < ncol) {
	contNonBasicTab[card_rowTab][ii] = z[locind];
      }
      else {
	contNonBasicTab[card_rowTab][ii] = slack[locind - ncol];
      }
      if(fabs(contNonBasicTab[card_rowTab][ii]) > TINY)
	nonzeroes++;
    }

    for(ii=0; ii<card_intNonBasicVar; ii++) {
      int locind = intNonBasicVar[ii];
      if(locind < ncol) {
	intNonBasicTab[card_rowTab][ii] = z[locind];
      }
      else {
	printf("### ERROR: CglRedSplit2::generateCuts(): integer slack unexpected\n");
	exit(1);
      }
      if(fabs(intNonBasicTab[card_rowTab][ii]) > TINY)
	nonzeroes++;
    }

#if RS_FAST_INT > 0 || RS_FAST_CONT > 0 || RS_FAST_WORK > 0
    if(nonzeroes < param.getMaxNonzeroesTab()) {
#endif
      // The number of nonzeroes is small enough, we keep the row
      card_rowTab++;
      card_intBasicVar++;
      if (cv_intBasicVar_frac[basis_index[ind_row]] == 1)
	card_intBasicVar_frac++;
#if RS_FAST_INT > 0 || RS_FAST_CONT > 0 || RS_FAST_WORK > 0
    }
#endif
  }
#if RS_FAST_INT == 0 && RS_FAST_CONT == 0 && RS_FAST_WORK == 0
  rs_allocmatINT(&pi_mat, mTab, mTab);
#else
  // Allow for sparse info (and for work)
  rs_allocmatINT(&pi_mat, mTab, mTab+
		 card_intNonBasicVar+2*card_contNonBasicVar+3);
  for (int i=0;i<mTab;i++) {
    int * info = pi_mat[i];
    int * intSparse=info+mTab;
    int * contSparse=intSparse+1+card_intNonBasicVar;
    int n;
    const double * values;
    n=0;
    values=intNonBasicTab[i];
    for (int j=0;j<card_intNonBasicVar;j++) {
      if (fabs(values[j]) > TINY)
	intSparse[++n]=j;
    }
    intSparse[0]=n; 
    n=0;
    values=contNonBasicTab[i];
    for (int j=0;j<card_contNonBasicVar;j++) {
      if (fabs(values[j]))
	contSparse[++n]=j;
    }
    contSparse[0]=n; 
  }
#endif

#ifdef RS2_TRACE
  printf("intBasicVar:\n");
  for(i=0; i<card_intBasicVar; i++) {
    printf("%d ", intBasicVar[i]);
  }
  printf("\n");
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

  int card_row;
  double *row = new double[ncol+nrow];
  int *rowind = new int[ncol];
  double *rowelem = new double[ncol];

  const double *elements = byRow->getElements();
  const CoinBigIndex *rowStart = byRow->getVectorStarts();
  const int *indices = byRow->getIndices();
  const int *rowLength = byRow->getVectorLengths(); 

  CglRedSplit2Param::ColumnSelectionStrategy columnSelection;
  CglRedSplit2Param::RowSelectionStrategy rowSelection;
  int numRows;

  int maxNumComputedCuts = param.getMaxNumComputedCuts();

  int* bufflambda = lambda;
  OsiCuts* buffcs = cs;
  // quality of cuts
  struct sortElement* quality = NULL; 
  if (maxNumCuts < maxNumComputedCuts){
    // user wants the generator to generate several cuts and select only best
    if (lambda != NULL){
      bufflambda = new int[maxNumComputedCuts*nrow];
    }
    if (cs != NULL){
      buffcs = new OsiCuts();
    }
    quality = new struct sortElement[maxNumComputedCuts];
  }

  int initNumCuts = 0;
  if (buffcs != NULL){
    initNumCuts = buffcs->sizeRowCuts();
  }

  // reset vector of lambdas
  if (bufflambda != NULL){
    if (maxNumCuts < maxNumComputedCuts){
      memset(bufflambda, 0, maxNumComputedCuts*nrow*sizeof(int));
    }
    else{
      memset(bufflambda, 0, maxNumCuts*nrow*sizeof(int));
    }
  }


  int numCuts = 0;

  std::vector<CglRedSplit2Param::ColumnSelectionStrategy> listColSel = 
    param.getColumnSelectionStrategy();
  std::vector<CglRedSplit2Param::RowSelectionStrategy> listRowSel = 
    param.getRowSelectionStrategy();
  std::vector<int> listNumRows = param.getNumRowsReduction();
#if 0
  double work= listColSel.size()*listNumRows.size()*listRowSel.size()*nDiag;
  work *= mTab*(card_intNonBasicVar+card_contNonBasicVar);
  printf("listColSel %d listNumRows %d listRowSel %d mTab %d intNon %d contnon %d  - %d generate_rows and dense ops %g\n",
	 listColSel.size(),listNumRows.size(),listRowSel.size(),mTab,
	 card_intNonBasicVar,card_contNonBasicVar,
	 listColSel.size()*listNumRows.size()*listRowSel.size()*nDiag,work);
#endif
  for (unsigned int coliter = 0; coliter < listColSel.size(); ++coliter){
    if (!checkTime() || numCuts >= maxNumComputedCuts){
      break;
    }
    columnSelection = listColSel[coliter];
#ifdef RS2_TRACE
    printf("Filling workNonBasicTab with column strategy %d\n", columnSelection);
#endif
    // select the columns
    fill_workNonBasicTab(columnSelection);
    for (unsigned int nriter = 0; nriter < listNumRows.size(); ++nriter){
      if (!checkTime() || numCuts >= maxNumComputedCuts){
	break;
      }
      numRows = listNumRows[nriter];
#ifdef RS2_TRACE
      printf("Applying reduction algorithm with %d rows\n", numRows);
#endif
      for (unsigned int rowiter = 0; rowiter < listRowSel.size(); ++rowiter){
	if (!checkTime() || numCuts >= maxNumComputedCuts){
	  break;
	}
	rowSelection = listRowSel[rowiter];
#ifdef RS2_TRACE
	printf("Reducing norms with row strategy %d\n", rowSelection);
#endif
	// new iteration: reinitialize pi_mat
	for(i=0; i<mTab; i++) {
	  memset(pi_mat[i], 0, mTab*sizeof(int));
	  if (param.getSkipGomory() == false){
	    pi_mat[i][i] = 1;
	  }
	}
	reduce_workNonBasicTab(numRows, rowSelection, mTab);
	// now generate cuts
	for(i=0; i<mTab && numCuts < maxNumComputedCuts; i++) {
	  if (pi_mat[i][i] != 0){
	    card_row = 0;
	    generate_row(i, row);
	    flip(row);

	    // RHS of equalities after flipping/translating non basic variables
	    // is given by the current LP solution (re-ordered according to 
	    // basis_index), i.e. what is now in rhsTab. 
	    // RHS of row i of (pi_mat * contNonBasicTab) is then simply 
	    // pi_mat[i] * rhsTab

	    double tabrowrhs = rs_dotProd(pi_mat[i], rhsTab, mTab); 
	    int got_one = 0;
    
	    got_one = generate_cgcut(row, &tabrowrhs);

	    if (!got_one && param.getSkipGomory() == false){
	      // Cut improvement failed and we want to generate Gomory
	      // cuts from single rows; so, revert back to one row and try
	      // the corresponding Gomory cut. If we still fail, we give up
	      memset(pi_mat[i], 0, mTab*sizeof(int));
	      pi_mat[i][i] = 1;
	      card_row = 0;
	      generate_row(i, row);
	      flip(row);
	      tabrowrhs = rhsTab[i];
	      got_one = generate_cgcut(row, &tabrowrhs);
	    }


	    if (got_one){
	      bool skip = false;
	      if (bufflambda != NULL){
		for (int k = 0; k < mTab; ++k){
		  bufflambda[numCuts*nrow + origRow[k]] = pi_mat[i][k];
		}
		for (int k = 0; k < numCuts; ++k){
		  if (!rs_are_different_vectors(bufflambda + k*nrow, 
						bufflambda + numCuts*nrow, 
						nrow)){
		    // Cut was previously generated; skip
		    skip = true;
		  }
		}
	      }
	      if (!skip){
		// new multipliers found, check if the generated cut is good
		unflip(row, &tabrowrhs);

		eliminate_slacks(row, elements, rowStart, indices, 
				 rowLength, effective_rhs, &tabrowrhs);

		if(generate_packed_row(xlp, row, rowind, rowelem, &card_row, 
				       tabrowrhs)) {
		  // the cut is good, multipliers should be saved.
		  // record quality of cut if needed
		  if (quality){
		    quality[numCuts].index = numCuts;
		    quality[numCuts].cost = tabrowrhs;
		    double rnorm = 0.0;
		    for (int ii = 0; ii < card_row; ++ii){
		      rnorm += rowelem[ii]*rowelem[ii];
		      quality[numCuts].cost -= rowelem[ii]*xlp[rowind[ii]];
		    }
		    quality[numCuts].cost /= rnorm;
		  }
		  // now, add cut to the collection if needed
		  if (cs != NULL){
#if CBC_CHECK_CUT_LENGTH
		    if (card_row<CBC_CHECK_CUT_LENGTH*nrow) {
#endif
		      OsiRowCut rc;
		      rc.setRow(card_row, rowind, rowelem);
		      rc.setLb(-param.getINFINIT());
		      double adjust = param.getEPS_RELAX_ABS();
		      if(param.getEPS_RELAX_REL() > 0.0) {
			adjust += fabs(tabrowrhs) * param.getEPS_RELAX_REL();
		      }
		      rc.setUb(tabrowrhs + adjust);   
		      // relax the constraint slightly
		      buffcs->insertIfNotDuplicateAndClean(rc, 51, CoinAbsFltEq(param.getEPS()));
		      numCuts = buffcs->sizeRowCuts() - initNumCuts;
#if CBC_CHECK_CUT_LENGTH
		    }
#endif
		  }
		  else{
		    numCuts++;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  if (maxNumCuts < maxNumComputedCuts){
    // select best cuts and copy them back
    if (numCuts > maxNumCuts){
      qsort(quality, numCuts, sizeof(struct sortElement), rs2_compareElements);
    }
    // also delete temp data
    if (buffcs){
      for (int i = 0; i < numCuts && i < maxNumCuts; ++i){
	cs->insertIfNotDuplicateAndClean(buffcs->rowCut(quality[i].index),52,
				 CoinAbsFltEq(param.getEPS_COEFF()));
      }
      delete buffcs;
    }
    if (bufflambda){
      for (int i = 0; i < numCuts && i < maxNumCuts; ++i){
	memcpy(lambda + (i*nrow), bufflambda + (quality[i].index*nrow), 
	       sizeof(int)*nrow);
      }
      delete[] bufflambda;
    }
    if (numCuts > maxNumCuts){
      numCuts = maxNumCuts;
    }
  }

  delete[] cstat;
  delete[] rstat;
  delete[] basis_index;
  delete[] slack;
  delete[] z;
  delete[] effective_rhs;
  delete[] row;
  delete[] rowind;
  delete[] rowelem;

  delete[] origRow;
  delete[] cv_intBasicVar;  
  delete[] cv_intBasicVar_frac;  
  delete[] cv_fracRowsTab;  
  delete[] intBasicVar;
  delete[] intBasicVar_frac;
  delete[] intNonBasicVar;
  delete[] contNonBasicVar;
  delete[] nonBasicAtUpper;
  delete[] nonBasicAtLower;
  delete[] is_integer;
  rs_deallocmatDBL(&contNonBasicTab, mTab);
  rs_deallocmatDBL(&workNonBasicTab, mTab);
  rs_deallocmatDBL(&intNonBasicTab, mTab);
  rs_deallocmatINT(&pi_mat, mTab);
  delete[] rhsTab;
  delete[] norm;

  return numCuts;
} /* generateCuts */

/***********************************************************************/
void CglRedSplit2::setParam(const CglRedSplit2Param &source) {
  param = source;
} /* setParam */

/***********************************************************************/
void CglRedSplit2::compute_is_integer() {

  int i;
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
} /* compute_is_integer */

/***********************************************************************/
void CglRedSplit2::print() const
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
void CglRedSplit2::printOptTab(OsiSolverInterface *lclSolver) const
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
CglRedSplit2::CglRedSplit2() :
  CglCutGenerator(),
  nrow(0),
  ncol(0),
  card_intBasicVar(0),
  card_intBasicVar_frac(0),
  card_intNonBasicVar(0),
  card_contNonBasicVar(0),
  card_nonBasicAtUpper(0),
  card_nonBasicAtLower(0),
  cv_intBasicVar(0),
  cv_intBasicVar_frac(0),
  cv_fracRowsTab(0),
  intBasicVar(0),
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
  rhsTab(0)
{
}

/*********************************************************************/
CglRedSplit2::CglRedSplit2(const CglRedSplit2Param &RS_param) :
  CglCutGenerator(),
  nrow(0),
  ncol(0),
  card_intBasicVar(0),
  card_intBasicVar_frac(0),
  card_intNonBasicVar(0),
  card_contNonBasicVar(0),
  card_nonBasicAtUpper(0),
  card_nonBasicAtLower(0),
  cv_intBasicVar(0),
  cv_intBasicVar_frac(0),
  cv_fracRowsTab(0),
  intBasicVar(0),
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
  rhsTab(0)
{
  param = RS_param;
}
 
/*********************************************************************/
CglRedSplit2::CglRedSplit2 (const CglRedSplit2 & source) :
  CglCutGenerator(source),
  param(source.param),
  nrow(0),
  ncol(0),
  card_intBasicVar(0),
  card_intBasicVar_frac(0),
  card_intNonBasicVar(0),
  card_contNonBasicVar(0),
  card_nonBasicAtUpper(0),
  card_nonBasicAtLower(0),
  cv_intBasicVar(NULL),
  cv_intBasicVar_frac(NULL),
  cv_fracRowsTab(NULL),
  intBasicVar(NULL),
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
  rhsTab(NULL)
{
}

/*********************************************************************/
CglCutGenerator *
CglRedSplit2::clone() const
{
  return new CglRedSplit2(*this);
}

/*********************************************************************/
CglRedSplit2::~CglRedSplit2 ()
{}

/*********************************************************************/
CglRedSplit2 &
CglRedSplit2::operator=(const CglRedSplit2 &source)
{  

  if (this != &source) {
    CglCutGenerator::operator=(source);
    param = source.param;
  }
  return *this;
}
/*********************************************************************/
// Returns true if needs optimal basis to do cuts
bool 
CglRedSplit2::needsOptimalBasis() const
{
  return true;
}

/************************************************************************/

int CglRedSplit2::generateMultipliers(const OsiSolverInterface &si, 
				      int* lambda,
				      int maxNumMultipliers,
				      int* basicVariables,
				      OsiCuts* cs)
{
  solver = const_cast<OsiSolverInterface *>(&si);
  if(solver == NULL) {
    printf("### WARNING: CglRedSplit2::generateCuts(): no solver available.\n");
    return 0;    
  }  

  if(!solver->optimalBasisIsAvailable()) {
    printf("### WARNING: CglRedSplit2::generateCuts(): no optimal basis available.\n");
    return 0;
  }

  // Reset some members of CglRedSplit2
  card_intBasicVar = 0;
  card_intBasicVar_frac = 0;
  card_intNonBasicVar = 0;
  card_contNonBasicVar = 0;
  card_nonBasicAtUpper = 0;
  card_nonBasicAtLower = 0;
  numRedRows = 0;
  startTime = CoinCpuTime();

  // Get basic problem information from solver
  ncol = solver->getNumCols(); 
  nrow = solver->getNumRows(); 
  colLower = solver->getColLower();
  colUpper = solver->getColUpper();
  rowLower = solver->getRowLower();
  rowUpper = solver->getRowUpper();
  rowRhs = solver->getRightHandSide();
  reducedCost = solver->getReducedCost();
  rowPrice = solver->getRowPrice();
  objective = solver->getObjCoefficients();

  xlp = solver->getColSolution();
  rowActivity = solver->getRowActivity();
  byRow = solver->getMatrixByRow();

  solver->enableFactorization();
  if (basicVariables != NULL){
    solver->getBasics(basicVariables);
  }
  int numMultipliers = generateCuts(cs, maxNumMultipliers, lambda);
  solver->disableFactorization();
  return numMultipliers;
} /* generateMultipliers */

/***********************************************************************/

int CglRedSplit2::tiltLandPcut(const OsiSolverInterface* si, 
			       double* landprow, double landprhs, 
			       int rownumber, const double* xbar, 
			       const int* newnonbasics, OsiRowCut* cs,
			       int* lambda){
  solver = const_cast<OsiSolverInterface*>(si);
  if(solver == NULL) {
    printf("### WARNING: CglRedSplit2::tiltLandPcut(): no solver available.\n");
    return 0;   
  }    
 
  // Reset some members of CglRedSplit2
  card_intBasicVar = 0;
  card_intBasicVar_frac = 0;
  card_intNonBasicVar = 0;
  card_contNonBasicVar = 0;
  card_nonBasicAtUpper = 0;
  card_nonBasicAtLower = 0;
  numRedRows = 0;
  startTime = CoinCpuTime();

  // Get basic problem information from solver
  ncol = solver->getNumCols(); 
  nrow = solver->getNumRows(); 
  colLower = solver->getColLower();
  colUpper = solver->getColUpper();
  rowLower = solver->getRowLower();
  rowUpper = solver->getRowUpper();
  rowRhs = solver->getRightHandSide();
  reducedCost = solver->getReducedCost();
  rowPrice = solver->getRowPrice();
  objective = solver->getObjCoefficients();

  xlp = solver->getColSolution();
  rowActivity = solver->getRowActivity();
  byRow = solver->getMatrixByRow();

  int i;
  is_integer = new int[ncol]; 
  
  compute_is_integer();

  int *cstat = new int[ncol];
  int *rstat = new int[nrow];


  solver->getBasisStatus(cstat, rstat);   // 0: free  1: basic  
                                          // 2: upper 3: lower

  int *basis_index = new int[nrow]; // basis_index[i] = 
                                    //        index of pivot var in row i
                                    //        (slack if number >= ncol) 

#ifdef RS2_TRACETAB
  rs_printvecINT("cstat", cstat, ncol);
  rs_printvecINT("rstat", rstat, nrow);
#endif
  solver->getBasics(basis_index);

  cv_intBasicVar = new int[ncol];  
  cv_intBasicVar_frac = new int[ncol];  
  intBasicVar = new int[ncol];                                 
  intNonBasicVar = new int[ncol];       
  contNonBasicVar = new int[ncol+nrow]; 
  nonBasicAtUpper = new int[ncol+nrow]; 
  nonBasicAtLower = new int[ncol+nrow]; 
  double dist_int;

  for(i=0; i<ncol; i++) {
    cv_intBasicVar[i] = 0;
    cv_intBasicVar_frac[i] = 0;

    switch(cstat[i]) {
    case 1: // basic variable
      
      dist_int = rs_above_integer(xlp[i]);
      if(is_integer[i]){
	if ((dist_int > param.getAway()) && (dist_int < 1 - param.getAway())) {
	  cv_intBasicVar_frac[i] = 1;
	  card_intBasicVar_frac++;

	  // intBasicVar_frac computed below, 
	  // as order must be according to selected rows
	}
	card_intBasicVar++;
	cv_intBasicVar[i] = 1;

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
      printf("### ERROR: CglRedSplit2::generateCuts(): cstat[%d]: %d\n",
	     i, cstat[i]);
      exit(1);
      break;
    } 
  }

  // Use this instead of rowRhs to allow for ranges
  double *effective_rhs = new double[nrow];

  for(i=0; i<nrow; i++) {
    // effective rhs
    effective_rhs[i] = rowRhs[i];
    switch(rstat[i]) {
    case 1: // basic slack
      break;

    case 2: // non basic slack at upper; flipped and shifted
      effective_rhs[i] = rowLower[i];
      nonBasicAtUpper[card_nonBasicAtUpper] = ncol+i;
      card_nonBasicAtUpper++;

      contNonBasicVar[card_contNonBasicVar] = ncol+i;
      card_contNonBasicVar++;
      break;

    case 3: // non basic slack at lower; shifted
      effective_rhs[i] = rowUpper[i];
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
      //exit(1);
      break;
    }
    assert (fabs(effective_rhs[i])<1.0e100);
  }

#ifdef RS2_TRACE
  printf("CglRedSplit2()::card_intBasicVar_frac %d %d\n", 
	 card_intBasicVar_frac, card_contNonBasicVar);
#endif
  if((card_contNonBasicVar == 0) || (card_intBasicVar == 0)) {
    delete[] cstat;
    delete[] rstat;
    delete[] basis_index;

    delete[] cv_intBasicVar;  
    delete[] cv_intBasicVar_frac;  
    delete[] intBasicVar;
    delete[] intNonBasicVar;
    delete[] contNonBasicVar;
    delete[] nonBasicAtUpper;
    delete[] nonBasicAtLower;
    delete[] is_integer;
    delete[] effective_rhs;

    printf("No vars to generate cut\n");
    return 0; // no cuts can be generated
  }

  double *z = new double[ncol];  // workspace to get row of the tableau
  double *slack = new double[nrow];  // workspace to get row of the tableau

#ifdef RS2_TRACETAB
  printOptTab(solver);
#endif

  if (rownumber >= 0){
    // this means that the source row is in the given simplex tableau
    mTab = card_intBasicVar;
  }
  else{
    // source row is an extra row, so we need extra space
    mTab = card_intBasicVar + 1;
  }
  nTab = card_contNonBasicVar;

  rhsTab = new double[mTab];
  cv_fracRowsTab = new int[mTab];
  memset(cv_fracRowsTab, 0, mTab*sizeof(int));
  int card_rowTab = 0;

  int extraColumns = 0;
  while (newnonbasics[extraColumns] >= 0){
    extraColumns++;
  }

  rs_allocmatDBL(&contNonBasicTab, mTab, card_contNonBasicVar);
  rs_allocmatDBL(&workNonBasicTab, mTab, card_contNonBasicVar + extraColumns);
  rs_allocmatDBL(&intNonBasicTab, mTab, card_intNonBasicVar);
  norm = new double[mTab];

  intBasicVar_frac = new int[ncol];                                 

  card_intBasicVar = 0; // recompute in pivot order
  card_intBasicVar_frac = 0;

  // write L&P row as first row
  rhsTab[card_rowTab] = landprhs;

  // flip row: L&P flips the nonbasic at upper bound, but we don't want this
  for (int ii=0; ii<card_nonBasicAtUpper; ii++){
    landprow[nonBasicAtUpper[ii]] = -landprow[nonBasicAtUpper[ii]];
  }

  for (int ii=0; ii<card_contNonBasicVar; ii++) {
    int locind = contNonBasicVar[ii];
    contNonBasicTab[card_rowTab][ii] = landprow[locind];
  }

  for (int ii=0; ii<card_intNonBasicVar; ii++) {
    int locind = intNonBasicVar[ii];
    intNonBasicTab[card_rowTab][ii] = landprow[locind];
  }

  // flip back to restore L&P row as it was given
  for (int ii=0; ii<card_nonBasicAtUpper; ii++){
    landprow[nonBasicAtUpper[ii]] = -landprow[nonBasicAtUpper[ii]];
  }

  card_rowTab++;

  for(i=0; i<nrow; i++) {

    int ind_row = i;

    if(basis_index[ind_row] >= ncol || ind_row == rownumber) {
      continue;
    } 

    if(cv_intBasicVar[basis_index[ind_row]] == 1) { 
      // row used in generation
      intBasicVar[card_intBasicVar] = basis_index[ind_row];
      if (cv_intBasicVar_frac[basis_index[ind_row]] == 1){
	// row is fractional
	intBasicVar_frac[card_intBasicVar_frac] = basis_index[ind_row];
	card_intBasicVar_frac++;
	cv_fracRowsTab[card_intBasicVar] = 1;
      }
      card_intBasicVar++;
      int ii;

      rhsTab[card_rowTab] = xlp[basis_index[ind_row]];
      solver->getBInvARow(ind_row, z, slack);

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
	  printf("### ERROR: CglRedSplit2::generateCuts(): integer slack unexpected\n");
	  exit(1);
	}
      }

      card_rowTab++;
    }
  }

  int pi_mat_rows = 1;
  rs_allocmatINT(&pi_mat, pi_mat_rows, mTab);

#ifdef RS2_TRACE
  printf("intBasicVar:\n");
  for(i=0; i<card_intBasicVar; i++) {
    printf("%d ", intBasicVar[i]);
  }
  printf("\n");
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

  int card_row;
  double *row = new double[ncol+nrow];
  int *rowind = new int[ncol];
  double *rowelem = new double[ncol];

  const double *elements = byRow->getElements();
  const CoinBigIndex *rowStart = byRow->getVectorStarts();
  const int *indices = byRow->getIndices();
  const int *rowLength = byRow->getVectorLengths(); 

  CglRedSplit2Param::ColumnSelectionStrategy columnSelection;
  CglRedSplit2Param::RowSelectionStrategy rowSelection;
  int numRows;

  int generatedCuts = 0;

  std::vector<CglRedSplit2Param::ColumnSelectionStrategy> listColSel 
    = param.getColumnSelectionStrategyLAP();
  std::vector<CglRedSplit2Param::RowSelectionStrategy> listRowSel 
    = param.getRowSelectionStrategyLAP();
  std::vector<int> listNumRows = param.getNumRowsReductionLAP();

  for (unsigned int coliter = 0; coliter < listColSel.size(); ++coliter){
    if (!checkTime()){
      break;
    }
    columnSelection = listColSel[coliter];
#ifdef RS2_TRACE
    printf("Filling workNonBasicTab with column strategy %d\n", columnSelection);
#endif
    nTab = 0;
    // select the columns
    if (columnSelection != CglRedSplit2Param::CS_LAP_NONBASICS){
      fill_workNonBasicTab(columnSelection, newnonbasics);
    }
    fill_workNonBasicTab(newnonbasics, xbar, 
			 param.getColumnScalingStrategyLAP());
    for (unsigned int nriter = 0; nriter < listNumRows.size(); ++nriter){
      if (!checkTime()){
	break;
      }
      numRows = listNumRows[nriter];
#ifdef RS2_TRACE
      printf("Applying reduction algorithm with %d rows\n", numRows);
#endif
      for (unsigned int rowiter = 0; rowiter < listRowSel.size(); ++rowiter){
	if (!checkTime()){
	  break;
	}
	rowSelection = listRowSel[rowiter];
#ifdef RS2_TRACE
	printf("Reducing norms with row strategy %d\n", rowSelection);
#endif
	// new iteration: reinitialize pi_mat
	for(i=0; i<pi_mat_rows; i++) {
	  memset(pi_mat[i], 0, mTab*sizeof(int));
	}
	// reduce the coefficients
	reduce_workNonBasicTab(numRows, rowSelection, 1);
	// now generate cuts
	for(i=0; i<pi_mat_rows; i++) {
	  if (pi_mat[i][i] == 0){
	    // means we did not generate a cut from this row
	    continue;
	  }
	  card_row = 0;
 	  generate_row(i, row);
 	  flip(row);

	  // compute RHS
	  double tabrowrhs = rs_dotProd(pi_mat[i], rhsTab, mTab); 
	  int got_one = 0;
    
 	  got_one = generate_cgcut(row, &tabrowrhs);

	  if (got_one){
	    // new multipliers found, check if the generated cut is good
	    unflip(row, &tabrowrhs);

	    eliminate_slacks(row, elements, rowStart, indices, 
			     rowLength, effective_rhs, &tabrowrhs);
	    
	    if(generate_packed_row(xbar, row, rowind, rowelem, &card_row, 
				   tabrowrhs)) {
	      // the cut is good, multipliers should be saved;
	      // first of all, add cut to the collection if needed
	      cs->setRow(card_row, rowind, rowelem);
	      cs->setLb(-param.getINFINIT());
	      double adjust = param.getEPS_RELAX_ABS();
	      if(param.getEPS_RELAX_REL() > 0.0) {
		adjust += fabs(tabrowrhs) * param.getEPS_RELAX_REL();
	      }
	      cs->setUb(tabrowrhs + adjust);   
	      if (lambda){
		// Modify the initial disjunction by adding the new 
		// coefficients
		for (int k = 1; k < mTab; ++k){
		  lambda[intBasicVar[k-1]] += pi_mat[i][k];
		}
	      }
	      generatedCuts++;
	    }
	  }
	}
      }
    }
  }

  delete[] cstat;
  delete[] rstat;
  delete[] basis_index;
  delete[] slack;
  delete[] z;
  delete[] effective_rhs;
  delete[] row;
  delete[] rowind;
  delete[] rowelem;

  delete[] cv_intBasicVar_frac;  
  delete[] cv_fracRowsTab;  
  delete[] intBasicVar;
  delete[] intBasicVar_frac;
  delete[] intNonBasicVar;
  delete[] contNonBasicVar;
  delete[] nonBasicAtUpper;
  delete[] nonBasicAtLower;
  delete[] is_integer;
  rs_deallocmatDBL(&contNonBasicTab, mTab);
  rs_deallocmatDBL(&workNonBasicTab, mTab);
  rs_deallocmatDBL(&intNonBasicTab, mTab);
  rs_deallocmatINT(&pi_mat, pi_mat_rows);
  delete[] rhsTab;
  delete[] norm;

  return generatedCuts;
  
} /* tiltLandPcut */

/***************************************************************************/
void CglRedSplit2::fill_workNonBasicTab(const int* newnonbasics,
					const double* xbar,
					CglRedSplit2Param::ColumnScalingStrategy scaling){
  
  // See CglRedSplit2Param.hpp for a description of the strategies
  int i = 0;
  int currvar, colpos;
#ifdef RS2_TRACE
  while (newnonbasics[i] >= 0){
    printf("Newnonbasic[%d]: %d\n", i, newnonbasics[i]);
    i++;
  }
  i = 0;
#endif
  const char * intVar = solver->getColType();
  // In workNonBasicTab, we write the new nonbasic columns given by the
  // Lift & Project algorithm, scaled by the value of the fractional solution.
  while (newnonbasics[i] >= 0){
    currvar = newnonbasics[i];
    if (currvar < ncol && intVar[currvar]){
      for (colpos = 0; colpos < card_intNonBasicVar; ++colpos){
	if (intNonBasicVar[colpos] == currvar){
#ifdef RS2_TRACE
	  printf("Found int nonbasic at pos %d: %d %d\n", colpos, intNonBasicVar[colpos], currvar);
#endif
	  break;
	}
      }
      double factor = 1.0;
      if (scaling == CglRedSplit2Param::SC_LINEAR){
	if (fabs(xbar[currvar]) > factor){
	  factor = fabs(xbar[currvar]);
	}
      }
      else if (scaling == CglRedSplit2Param::SC_LINEAR_BOUNDED){
	if (fabs(xbar[currvar]) > factor){
	  factor = fabs(xbar[currvar]);
	}
	if (factor < param.getColumnScalingBoundLAP()){
	  factor = param.getColumnScalingBoundLAP();
	}
      }
      else if (scaling == CglRedSplit2Param::SC_LOG_BOUNDED){
	if (log(fabs(xbar[currvar])) > factor){
	  factor = log(fabs(xbar[currvar]));
	}
	if (factor < param.getColumnScalingBoundLAP()){
	  factor = param.getColumnScalingBoundLAP();
	}
      }
      else if (scaling == CglRedSplit2Param::SC_UNIFORM){
	factor = param.getColumnScalingBoundLAP();
      }
      else if (scaling == CglRedSplit2Param::SC_UNIFORM_NZ){
	if (fabs(xbar[currvar]) > param.getEPS()){
	  factor = param.getColumnScalingBoundLAP();
	}
      }
      for (int j = 0; j < mTab; ++j){
	workNonBasicTab[j][nTab] = intNonBasicTab[j][colpos]*factor;
      }
      nTab++;
    }
    else{
      for (colpos = 0; colpos < card_contNonBasicVar; ++colpos){
	if (contNonBasicVar[colpos] == currvar){
#ifdef RS2_TRACE
	  printf("Found cont nonbasic at pos %d: %d %d\n", colpos, contNonBasicVar[colpos], currvar);
#endif
	  break;
	}
      }
      double factor = 1.0;
      if (scaling == CglRedSplit2Param::SC_LINEAR){
	if (fabs(xbar[currvar]) > factor){
	  factor = fabs(xbar[currvar]);
	}
      }
      else if (scaling == CglRedSplit2Param::SC_LINEAR_BOUNDED){
	if (fabs(xbar[currvar]) > factor){
	  factor = fabs(xbar[currvar]);
	}
	if (factor < param.getColumnScalingBoundLAP()){
	  factor = param.getColumnScalingBoundLAP();
	}
      }
      else if (scaling == CglRedSplit2Param::SC_LOG_BOUNDED){
	if (log(fabs(xbar[currvar])) > factor){
	  factor = log(fabs(xbar[currvar]));
	}
	if (factor < param.getColumnScalingBoundLAP()){
	  factor = param.getColumnScalingBoundLAP();
	}
      }
      else if (scaling == CglRedSplit2Param::SC_UNIFORM){
	factor = param.getColumnScalingBoundLAP();
      }
      else if (scaling == CglRedSplit2Param::SC_UNIFORM_NZ){
	if (fabs(xbar[currvar]) > param.getEPS()){
	  factor = param.getColumnScalingBoundLAP();
	}
      }
      for (int j = 0; j < mTab; ++j){
#ifdef RS2_TRACE
	printf("Col %d Row %d: %f, xbar %f\n", colpos, j, contNonBasicTab[j][colpos], xbar[currvar]);
#endif	
	workNonBasicTab[j][nTab] = contNonBasicTab[j][colpos]*factor;
      }
      nTab++;
    }
    i++;
  }
#ifdef RS2_TRACE
  printf("Printing workNonBasicTab:\n");
#endif
#if RS_FAST_WORK 
  int workOffset = mTab+card_intNonBasicVar+card_contNonBasicVar+2;
#endif
  for (i = 0; i < mTab; ++i) {
#ifdef RS2_TRACE
    for (int j = 0; j < nTab; ++j){
      printf("%.6f ", workNonBasicTab[i][j]);
    }
    printf("\n");
#endif
#if RS_FAST_WORK == 0
    norm[i] = rs_dotProd(workNonBasicTab[i], workNonBasicTab[i], nTab);
#elif RS_FAST_WORK == 1
    norm[i] = rs_dotProd(workNonBasicTab[i], workNonBasicTab[i], nTab);
    double value = rs_sparseDotProd(workNonBasicTab[i],
				       workNonBasicTab[i], 
				       pi_mat[i]+workOffset,
				       pi_mat[i]+workOffset);
    assert (value==norm[i]);
#else
    norm[i] = rs_sparseDotProd(workNonBasicTab[i],
				       workNonBasicTab[i], 
				       pi_mat[i]+workOffset,
				       pi_mat[i]+workOffset);
#endif
  }
} /* fill_workNonBasicTab */

/***********************************************************************/

