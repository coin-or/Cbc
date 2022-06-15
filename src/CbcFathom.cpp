// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#ifdef CBC_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcFathom.hpp"

// Default Constructor
CbcFathom::CbcFathom()
  : model_(NULL)
  , possible_(false)
{
}

// Constructor from model
CbcFathom::CbcFathom(CbcModel &model)
  : model_(&model)
  , possible_(false)
{
}
// Resets stuff if model changes
void CbcFathom::resetModel(CbcModel *model)
{
  model_ = model;
}

// Destructor
CbcFathom::~CbcFathom()
{
}

// update model
void CbcFathom::setModel(CbcModel *model)
{
  model_ = model;
}
#ifdef CBC_HAS_CLP

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CbcOsiSolver::CbcOsiSolver()
  : OsiClpSolverInterface()
{
  cbcModel_ = NULL;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface *
CbcOsiSolver::clone(bool /*copyData*/) const
{
  //assert (copyData);
  return new CbcOsiSolver(*this);
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CbcOsiSolver::CbcOsiSolver(
  const CbcOsiSolver &rhs)
  : OsiSolverInterface()
  , // Should not be needed but get warning
  OsiClpSolverInterface(rhs)
{
  cbcModel_ = rhs.cbcModel_;
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CbcOsiSolver::~CbcOsiSolver()
{
}

//-------------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CbcOsiSolver &
CbcOsiSolver::operator=(const CbcOsiSolver &rhs)
{
  if (this != &rhs) {
    OsiClpSolverInterface::operator=(rhs);
    cbcModel_ = rhs.cbcModel_;
  }
  return *this;
}
#endif
#ifdef CBC_HAS_CLP
#include "ClpNode.hpp"
#include "CbcSymmetry.hpp"
/* Do something when called from Clp fathomMany
   needed as Clp does not know about Cbc */
int fromFathomMany(int typeCall, infoForCbc * stuff)
{
  // pick up model
  CbcModel * model = reinterpret_cast<CbcModel *>(stuff->model);
  ClpSimplex * simplex = stuff->simplex;
  OsiClpSolverInterface * clpSolver =
    dynamic_cast<OsiClpSolverInterface *>(model->continuousSolver());
  ClpSimplex * fullSimplex = clpSolver->getModelPtr();
  int returnCode=0;
  if (typeCall==1) {
    int * whichColumn = stuff->info->whichColumn_;
    //printf("Column %d becomes column %d\n",iColumn,whichColumn[iColumn]);
    CbcSymmetry * rootSymmetry = model->rootSymmetryInfo();
    // expand current bounds
    int numberFullColumns = fullSimplex->numberColumns();
    double * fullLower = stuff->info->usefulCbc_->usefulData;
    double * fullUpper = fullLower + numberFullColumns;
    int * backward = reinterpret_cast<int *>(fullLower+4*numberFullColumns);
    memcpy(fullLower,fullSimplex->columnLower(),numberFullColumns*sizeof(double));
    memcpy(fullUpper,fullSimplex->columnUpper(),numberFullColumns*sizeof(double));
    int numberColumns = simplex->numberColumns();
    double * lower = simplex->columnLower();
    double * upper = simplex->columnUpper();
    for (int i=0;i<numberColumns;i++) {
      int jColumn = whichColumn[i];
      fullLower[jColumn]=lower[i];
      fullUpper[jColumn]=upper[i];
    }
    int nFixed = rootSymmetry->changeBounds2(fullLower,fullUpper,
					    clpSolver);
    if (nFixed) {
      int * fixedToZero = rootSymmetry->fixedToZero();
      int nFixedAlready=0;
      for (int i=0;i<nFixed;i++) {
	int jColumn = fixedToZero[i];
	int iColumn = backward[jColumn];
	if (iColumn>=0) {
	  upper[iColumn] = 0.0;
	} else {
	  nFixedAlready++;
	  //printf("column %d already fixed - bounds %g, %g\n",
	  //	 jColumn,fullLower[jColumn],fullUpper[jColumn]);
	  //assert (!upper[iColumn]);
	}
      }
      returnCode = nFixed-nFixedAlready;
    }
  } else if (typeCall==2) {
    int iColumn = stuff->sequence;
    int * whichColumn = stuff->info->whichColumn_;
    //printf("Column %d becomes column %d\n",iColumn,whichColumn[iColumn]);
    CbcSymmetry * rootSymmetry = model->rootSymmetryInfo();
    // expand current bounds
    int numberFullColumns = fullSimplex->numberColumns();
    double * fullLower = CoinCopyOfArray(fullSimplex->columnLower(),numberFullColumns);
    double * fullUpper = CoinCopyOfArray(fullSimplex->columnUpper(),numberFullColumns);
    int * backward = new int[numberFullColumns];
    for (int i=0;i<numberFullColumns;i++) { 
      backward[i]=-1;
    }
    int numberColumns = simplex->numberColumns();
    double * lower = simplex->columnLower();
    double * upper = simplex->columnUpper();
    for (int i=0;i<numberColumns;i++) {
      int jColumn = whichColumn[i];
      backward[jColumn]=i;
      fullLower[jColumn]=lower[i];
      fullUpper[jColumn]=upper[i];
    }
    int nFixed = rootSymmetry->fixSome(iColumn,fullLower,fullUpper);
    if (nFixed) {
      int * fixedToZero = rootSymmetry->fixedToZero();
      int nFixedAlready=0;
      for (int i=0;i<nFixed;i++) {
	int jColumn = fixedToZero[i];
	int iColumn = backward[jColumn];
	if (iColumn>=0) {
	  upper[iColumn] = 0.0;
	} else {
	  nFixedAlready++;
	  //printf("column %d already fixed - bounds %g, %g\n",
	  //	 jColumn,fullLower[jColumn],fullUpper[jColumn]);
	  //assert (!upper[iColumn]);
	}
      }
      returnCode = nFixed-nFixedAlready;
    }
    delete [] fullLower;
    delete [] fullUpper;
    delete [] backward;
  }
  return returnCode;
}
#endif
