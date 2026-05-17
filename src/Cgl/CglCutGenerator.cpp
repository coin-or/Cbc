// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdlib>
#include <cassert>
//#include <cfloat>
//#include <iostream>

#include "CoinPragma.hpp"
#include "CglCutGenerator.hpp"
#include "CoinHelperFunctions.hpp"

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CglCutGenerator::CglCutGenerator()
  : originalSolver_(NULL)
  , aggressive_(0)
  , canDoGlobalCuts_(false)
  , maxSeconds_(0.0)
{
  // nothing to do here
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CglCutGenerator::CglCutGenerator(
  const CglCutGenerator &source)
  : aggressive_(source.aggressive_)
  , canDoGlobalCuts_(source.canDoGlobalCuts_)
  , maxSeconds_(source.maxSeconds_)
{
   if (source.originalSolver_)
     originalSolver_ = source.originalSolver_->clone();
   else
     originalSolver_ = NULL;
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CglCutGenerator::~CglCutGenerator()
{
  delete originalSolver_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CglCutGenerator &
CglCutGenerator::operator=(
  const CglCutGenerator &rhs)
{
  if (this != &rhs) {
    aggressive_ = rhs.aggressive_;
    canDoGlobalCuts_ = rhs.canDoGlobalCuts_;
    maxSeconds_ = rhs.maxSeconds_;
    delete originalSolver_;
    if (rhs.originalSolver_)
      originalSolver_ = rhs.originalSolver_->clone();
    else
      originalSolver_ = NULL;
  }
  return *this;
}
bool CglCutGenerator::mayGenerateRowCutsInTree() const
{
  return true;
}
// Return true if needs optimal basis to do cuts
bool CglCutGenerator::needsOptimalBasis() const
{
  return false;
}
bool CglCutGenerator::needsOriginalModel() const
{
  return false;
}

#ifdef NDEBUG
#undef NDEBUG
#endif

#if 0
//--------------------------------------------------------------------------
// test EKKsolution methods.
//--------------------------------------------------------------------------
void
CglCutGenerator::unitTest()
{
}
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
