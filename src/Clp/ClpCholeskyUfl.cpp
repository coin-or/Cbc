// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "ClpConfig.h"

#ifndef CLP_HAS_AMD
#error "Need AMD to compile ClpCholeskyUfl."
#else
#include "amd.h"
#endif

#include "CoinPragma.hpp"
#include "ClpCholeskyUfl.hpp"
#include "ClpMessage.hpp"
#include "ClpInterior.hpp"
#include "CoinHelperFunctions.hpp"
#include "ClpHelperFunctions.hpp"
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpCholeskyUfl::ClpCholeskyUfl(int denseThreshold)
  : ClpCholeskyBase(denseThreshold)
{
  type_ = 14;
  L_ = NULL;
  c_ = NULL;

}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpCholeskyUfl::ClpCholeskyUfl(const ClpCholeskyUfl &rhs)
  : ClpCholeskyBase(rhs)
{
  abort();
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpCholeskyUfl::~ClpCholeskyUfl()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpCholeskyUfl &
ClpCholeskyUfl::operator=(const ClpCholeskyUfl &rhs)
{
  if (this != &rhs) {
    ClpCholeskyBase::operator=(rhs);
    abort();
  }
  return *this;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpCholeskyBase *ClpCholeskyUfl::clone() const
{
  return new ClpCholeskyUfl(*this);
}

/* Orders rows and saves pointer to matrix.and model */
int ClpCholeskyUfl::order(ClpInterior *model)
{
  int iRow;
  model_ = model;
  if (preOrder(false, true, doKKT_))
    return -1;
  permuteInverse_ = new int[numberRows_];
  permute_ = new int[numberRows_];
  double Control[AMD_CONTROL];
  double Info[AMD_INFO];

  amd_defaults(Control);
  //amd_control(Control);

  int returnCode = amd_order(numberRows_, choleskyStart_, choleskyRow_,
    permute_, Control, Info);
  delete[] choleskyRow_;
  choleskyRow_ = NULL;
  delete[] choleskyStart_;
  choleskyStart_ = NULL;
  //amd_info(Info);

  if (returnCode != AMD_OK) {
    std::cout << "AMD ordering failed" << std::endl;
    return 1;
  }
  for (iRow = 0; iRow < numberRows_; iRow++) {
    permuteInverse_[permute_[iRow]] = iRow;
  }
  return 0;
}

/* Does Symbolic factorization given permutation.
   This is called immediately after order.  If user provides this then
   user must provide factorize and solve.  Otherwise the default factorization is used
   returns non-zero if not enough memory */
int ClpCholeskyUfl::symbolic()
{
  return ClpCholeskyBase::symbolic();
}

/* Factorize - filling in rowsDropped and returning number dropped */
int ClpCholeskyUfl::factorize(const double *diagonal, int *rowsDropped)
{
  return ClpCholeskyBase::factorize(diagonal, rowsDropped);
}

void ClpCholeskyUfl::solve(double *region)
{
  ClpCholeskyBase::solve(region);
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
