// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
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
