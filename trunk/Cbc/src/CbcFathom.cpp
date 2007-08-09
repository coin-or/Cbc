// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>

#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcFathom.hpp"

// Default Constructor
CbcFathom::CbcFathom() 
  :model_(NULL),
   possible_(false)
{
}

// Constructor from model
CbcFathom::CbcFathom(CbcModel & model)
:
  model_(&model),
  possible_(false)
{
}
// Resets stuff if model changes
void 
CbcFathom::resetModel(CbcModel * model)
{
  model_=model;
}

// Destructor 
CbcFathom::~CbcFathom ()
{
}

// update model
void CbcFathom::setModel(CbcModel * model)
{
  model_ = model;
}
#ifdef COIN_HAS_CLP

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CbcOsiSolver::CbcOsiSolver ()
  : OsiClpSolverInterface()
{
  cbcModel_ = NULL;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface * 
CbcOsiSolver::clone(bool copyData) const
{
  assert (copyData);
  return new CbcOsiSolver(*this);
}


//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CbcOsiSolver::CbcOsiSolver (
                  const CbcOsiSolver & rhs)
  : OsiClpSolverInterface(rhs)
{
  cbcModel_ = rhs.cbcModel_;
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CbcOsiSolver::~CbcOsiSolver ()
{
}

//-------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CbcOsiSolver &
CbcOsiSolver::operator=(const CbcOsiSolver& rhs)
{
  if (this != &rhs) { 
    OsiClpSolverInterface::operator=(rhs);
    cbcModel_ = rhs.cbcModel_;
  }
  return *this;
}
#endif
