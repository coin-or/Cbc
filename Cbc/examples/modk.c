// $Id: crew.cpp 2469 2019-01-06 23:17:46Z unxusr $
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).


/*! \file modk.c
  \brief Example of cut generator using callbacks in the CBC C API

  Tries to generate mod-k cuts, i.e. for each constraint
  rounds it down after multiplying it by {1 ... k-1} / k

*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <float.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <Cbc_C_Interface.h>

#define MIN(a,b) ( (a<b) ? (a) : (b) )

static int dblcmp( const void *pd1, const void *pd2 );
static char dbl_is_integer( const double value );
static int intcmp( const void *pd1, const void *pd2 );
static int mcd( int n, const int values[] );
static int dbl_as_int( double value );
static double try_cut( void *osiSolver, int row, const double mult, int num, int den, 
    int *nzCut, int *idx, double *coef, double *rhs, int *icoef );

void cutcallback( void *osiSolver, void *osiCuts, void *appData )
{
  printf("entered callback\n");
  int m = Osi_getNumRows( osiSolver );
  int n = Osi_getNumCols( osiSolver );

  // sorted coefs and rhs values
  int *scoef = (int *) malloc( sizeof(int)*(n+1) );

  // different integer positive coeffs (values for k)
  int *dcoef = (int *) malloc( sizeof(int)*(n+1) );

  // integer coefficients for computing the cut
  int *icoef = (int*) malloc( sizeof(int)*n );

  // double coefficients for storing the cut
  double *coef = (double*) malloc(sizeof(double)*n);

  // variables in this cut 
  int *cidx = (int*) malloc(sizeof(int)*n);

  for ( int i=0 ; i<m ; ++i )
  {
    // checking different coefficients
    int nz = Osi_getRowNz( osiSolver, i );

    double rhs = Osi_getRowRHS( osiSolver, i );
    if (!dbl_is_integer(rhs))
      continue;

    int irhs = dbl_as_int( rhs );

    char allInt = 1;
    const double *coeffs = Osi_getRowCoeffs( osiSolver, i );
    for ( int j=0 ; (j<nz) ; ++j )
    {
      if (!dbl_is_integer(coeffs[j]))
      {
        allInt = 0;
        break;
      }
      scoef[j] = abs(dbl_as_int(coeffs[j]));
    }

    if (!allInt)
      continue;

    scoef[nz++] = abs(irhs);

    qsort( scoef, nz, sizeof(int), intcmp );


    int prev = INT_MAX;
    int ndiff = 0;

    // checking for different values for k 
    for ( int j=0 ; (j<nz) ; ++j )
    {
      if (prev!=scoef[j])
      {
        prev = scoef[j];
        dcoef[ndiff++] = scoef[j];
      } // different coefs
    } // all columns

    if (ndiff<=1)
      continue;

    int nRem = 0;

    // removing those which are divisors of larger values
    for ( int j=0 ; (j<ndiff-1) ; ++j )
    {
      for ( int jl=j+1 ; (jl<ndiff) ; ++jl )
      {
        if (!dcoef[j])
          continue;
        if ( dcoef[jl] % dcoef[j] == 0 )
        {
          dcoef[j] = INT_MAX;
          ++nRem;
          break;
        }
      } // larger
    } // all non-zeros

    qsort( dcoef, ndiff, sizeof(int), intcmp );
    ndiff -= nRem;

    for ( int ik=0 ; ik<ndiff ; ++ ik )
    {
      int k = dcoef[ik];
      for ( int num=1 ; (num<k) ; ++num )
      {
        int nzCut = 0;
        double rhsCut = 0.0;

        int nMult, startMult = 0;
        double mults[] = { 1.0, -1.0 };

        char sense = toupper(Osi_getRowSense( osiSolver, i ));

        switch (sense)
        {
          case 'E':
            {
              nMult = 2;
              break;
            }
          case 'G':
            {
              nMult = 1;
              startMult = 1;
              break;
            }
          case 'L':
            {
              nMult = 1;
              startMult = 0;
              break;
            }
          default:
            {
              fprintf( stderr, "unknown sense!");
              abort();
            }
        }

        for ( int iMult=startMult ; (iMult<nMult) ; ++iMult )
        {
          const double mult = mults[iMult];
          const double viol = try_cut( osiSolver, i, mult, num, k, &nzCut,
              cidx, coef, &rhsCut, icoef );
          if ( viol<0.001 )
            continue;

          OsiCuts_addRowCut( osiCuts, nzCut, cidx, coef, 'L', rhsCut );
        } // constraint in the form <= 
      } // every numerator < k
    } // every denominator k
  } // every row

  free( scoef );
  free( dcoef );
  free( cidx );
  free( coef );
  free( icoef );
}

static int dblcmp( const void *pd1, const void *pd2 )
{
  const double *d1 = (const double *) pd1;
  const double *d2 = (const double *) pd2;

  if ((*d1)<(*d2))
    return -1;
  if ((*d1)>(*d2))
    return 1;

  return 0;
}

static char dbl_is_integer( const double value )
{
  double rv = floor( value+0.5 );

  return fabs(rv-value) < 1e-6;
}

static int intcmp( const void *pi1, const void *pi2 )
{
  const int *i1 = (const int *) pi1;
  const int *i2 = (const int *) pi2;

  return ( (*i1) - (*i2) );
}

static int mcd( int n, const int values[] )
{
  //for ( int j=0 ; (j<n) ; ++j )
  //
  return 0;
}

static int dbl_as_int( double value )
{
  return (int) floor(value+0.5);
}

static double try_cut( void *osiSolver, int row, const double mult, int num, int den, 
    int *nzCut, int *idx, double *coef, double *rhs, int *icoef )
{
  int nz = Osi_getRowNz( osiSolver, row );
  const int *rIdx = Osi_getRowIndices( osiSolver, row );
  const double *rCoef = Osi_getRowCoeffs( osiSolver, row );
  int imult = dbl_as_int(mult);
  int irhs = dbl_as_int( imult*Osi_getRowRHS( osiSolver, row ) );
  for ( int j=0 ; j<nz ; ++j )
    icoef[j] = dbl_as_int( imult*rCoef[j] );

  const double *x = Osi_getColSolution( osiSolver );

  double sumLHS = 0.0;

  (*nzCut) = 0;
  for ( int j=0 ; j<nz ; ++j )
  {
    int col = rIdx[j];
    assert( col >= 0 && col<Osi_getNumCols(osiSolver) );

    int cc = (icoef[j]*num)/den;
    if (!cc)
      continue;

    idx[(*nzCut)] = col;
    coef[(*nzCut)] = cc;

    sumLHS += ((double)cc) * x[col];

    (*nzCut)++;
  }

  *rhs = (irhs*num)/den;

  return (*rhs) - sumLHS;
}

int main( int argc, char **argv )
{
  if (argc<2)
  {
    fprintf( stderr, "enter mps instance name.\n");
    exit(1);
  }

  Cbc_Model *mip = Cbc_newModel();

  Cbc_readMps( mip, argv[1] );

  Cbc_setParameter( mip, "cuts", "off" );
  Cbc_setParameter( mip, "gomory", "ifmove" );
  Cbc_setParameter( mip, "preprocess", "off" );

  Cbc_addCutCallback( mip,  cutcallback, "modk", NULL );


  Cbc_solve( mip );

  Cbc_deleteModel( mip );
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
