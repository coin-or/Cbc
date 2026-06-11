// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*
  Debug compile symbols for CoinPresolve.

  PRESOLVE_CONSISTENCY, PRESOLVE_DEBUG, and PRESOLVE_SUMMARY control
  consistency checking and debugging in the continuous presolve. See the
  comments in CoinPresolvePsdebug.hpp. DO NOT just define the symbols here in
  this file. Unless these symbols are consistent across all presolve code,
  you'll get something between garbage and a core dump.
*/

#include <stdio.h>

#include <cassert>
#include <iostream>

#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

#include "CoinPackedMatrix.hpp"
#include "CoinWarmStartBasis.hpp"
#include "OsiSolverInterface.hpp"

#include "OsiPresolve.hpp"
#include "CoinPresolveMatrix.hpp"

#if PRESOLVE_CONSISTENCY > 0 || PRESOLVE_DEBUG > 0 || PRESOLVE_SUMMARY > 0
#include "CoinPresolvePsdebug.hpp"
#include "CoinPresolveMonitor.hpp"
#endif
#include "CoinPresolveEmpty.hpp"
#include "CoinPresolveFixed.hpp"
#include "CoinPresolveSingleton.hpp"
#include "CoinPresolveDoubleton.hpp"
#include "CoinPresolveTripleton.hpp"
#include "CoinPresolveZeros.hpp"
#include "CoinPresolveSubst.hpp"
#include "CoinPresolveForcing.hpp"
#include "CoinPresolveDual.hpp"
#include "CoinPresolveTighten.hpp"
#include "CoinPresolveUseless.hpp"
#include "CoinPresolveDupcol.hpp"
#include "CoinPresolveImpliedFree.hpp"
#include "CoinPresolveIsolated.hpp"
#include "CoinPresolveParity.hpp"
#include "CoinMessage.hpp"

OsiPresolve::OsiPresolve()
  : originalModel_(NULL)
  , presolvedModel_(NULL)
  , nonLinearValue_(0.0)
  , originalColumn_(NULL)
  , originalRow_(NULL)
  , paction_(0)
  , ncols_(0)
  , nrows_(0)
  , nelems_(0)
  , presolveActions_(0)
  , numberPasses_(5)
{
}

OsiPresolve::~OsiPresolve()
{
  gutsOfDestroy();
}
// Gets rid of presolve actions (e.g.when infeasible)
void OsiPresolve::gutsOfDestroy()
{
  const CoinPresolveAction *paction = paction_;
  while (paction) {
    const CoinPresolveAction *next = paction->next;
    delete paction;
    paction = next;
  }
  delete[] originalColumn_;
  delete[] originalRow_;
  paction_ = NULL;
  originalColumn_ = NULL;
  originalRow_ = NULL;
}
#if DEBUG_PREPROCESS > 1
/*
  This code is intended to allow a known solution to be checked
  against presolve progress. debugSolution is set in CbcSolver
*/
double *debugSolution = NULL;
int debugNumberColumns = -1;
#endif
/* This version of presolve returns a pointer to a new presolved
   model.  NULL if infeasible

   doStatus controls activities required to transform an existing
   solution to match the presolved problem. I'd (lh) argue that this should
   default to false, but to maintain previous behaviour it defaults to true.
   Really, this is only useful if you've already optimised before applying
   presolve and also want to work with the solution after presolve.  I think
   that this is the less common case. The more common situation is to apply
   presolve before optimising.
*/
static thread_local int lastPresolveRows = 0;
static thread_local int lastPresolveCols = 0;

OsiSolverInterface *
OsiPresolve::presolvedModel(OsiSolverInterface &si,
  double feasibilityTolerance,
  bool keepIntegers,
  int numberPasses,
  const char *prohibited,
  bool doStatus,
  const char *rowProhibited,
  const double * scLower)
{
  int newRows = si.getNumRows();
  int newCols = si.getNumCols();
  if (newCols > lastPresolveCols) {
    CoinClearPresolveStats();
  }


  ncols_ = si.getNumCols();
  nrows_ = si.getNumRows();
  nelems_ = si.getNumElements();
  numberPasses_ = numberPasses;

  double maxmin = si.getObjSense();
  originalModel_ = &si;
  delete[] originalColumn_;
  originalColumn_ = new int[ncols_];
  delete[] originalRow_;
  originalRow_ = new int[nrows_];
  int i;
  for (i = 0; i < ncols_; i++)
    originalColumn_[i] = i;
  for (i = 0; i < nrows_; i++)
    originalRow_[i] = i;

  // result is 0 - okay, 1 infeasible, -1 go round again
  int result = -1;

  // User may have deleted - its their responsibility
  presolvedModel_ = NULL;
  // Messages
  CoinMessages msgs = CoinMessage(si.messages().language());
  // Only go round 100 times even if integer preprocessing
  int totalPasses = 100;
  while (result == -1) {

    // make new copy
    delete presolvedModel_;
    presolvedModel_ = si.clone();
    totalPasses--;

    // drop integer information if wanted
    if (!keepIntegers) {
      int i;
      for (i = 0; i < ncols_; i++)
        presolvedModel_->setContinuous(i);
    }

    CoinPresolveMatrix* probptr = construct_CoinPresolveMatrix(ncols_,
      maxmin,
      presolvedModel_,
      nrows_, nelems_, doStatus, nonLinearValue_, prohibited,
      rowProhibited);
    CoinPresolveMatrix& prob(*probptr);
    // make sure row solution correct
    if (doStatus) {
      double *colels = prob.colels_;
      int *hrow = prob.hrow_;
      CoinBigIndex *mcstrt = prob.mcstrt_;
      int *hincol = prob.hincol_;
      int ncols = prob.ncols_;

      double *csol = prob.sol_;
      double *acts = prob.acts_;
      int nrows = prob.nrows_;

      int colx;

      memset(acts, 0, nrows * sizeof(double));

      for (colx = 0; colx < ncols; ++colx) {
        double solutionValue = csol[colx];
        for (CoinBigIndex i = mcstrt[colx]; i < mcstrt[colx] + hincol[colx]; ++i) {
          int row = hrow[i];
          double coeff = colels[i];
          acts[row] += solutionValue * coeff;
        }
      }
    }

    // move across feasibility tolerance
    prob.feasibilityTolerance_ = feasibilityTolerance;

    /*
  Do presolve. Allow for the possibility that presolve might be ineffective
  (i.e., we're feasible but no postsolve actions are queued.
*/
    paction_ = presolve(&prob);
    result = 0;
    // Get rid of useful arrays
    prob.deleteStuff();
    /*
  This we don't need to do unless presolve actually reduced the system.
*/
    if (prob.status_ == 0 && paction_) {
      // Looks feasible but double check to see if anything slipped through
      int n = prob.ncols_;
      double *lo = prob.clo_;
      double *up = prob.cup_;
      int i;

      for (i = 0; i < n; i++) {
        if (up[i] < lo[i]) {
          if (up[i] < lo[i] - 1.0e-8) {
            // infeasible
            prob.status_ = 1;
          } else {
            up[i] = lo[i];
          }
        }
      }

      n = prob.nrows_;
      lo = prob.rlo_;
      up = prob.rup_;

      for (i = 0; i < n; i++) {
        if (up[i] < lo[i]) {
          if (up[i] < lo[i] - 1.0e-8) {
            // infeasible
            prob.status_ = 1;
          } else {
            up[i] = lo[i];
          }
        }
      }
    }
    /*
  If we're feasible, load the presolved system into the solver. Presumably we
  could skip model update and copying of status and solution if presolve took
  no action.
*/
    if (prob.status_ == 0) {

      update_model_CoinPresolveMatrix(prob,presolvedModel_, nrows_, ncols_, nelems_);

#if PRESOLVE_CONSISTENCY > 0
      if (doStatus) {
        int basicCnt = 0;
        int basicColumns = 0;
        int i;
        CoinPresolveMatrix::Status status;
        for (i = 0; i < prob.ncols_; i++) {
          status = prob.getColumnStatus(i);
          if (status == CoinPrePostsolveMatrix::basic)
            basicColumns++;
        }
        basicCnt = basicColumns;
        for (i = 0; i < prob.nrows_; i++) {
          status = prob.getRowStatus(i);
          if (status == CoinPrePostsolveMatrix::basic)
            basicCnt++;
        }

#if PRESOLVE_DEBUG > 0
        presolve_check_nbasic(&prob);
#endif
        if (basicCnt > prob.nrows_) {
          // Take out slacks
          double *acts = prob.acts_;
          double *rlo = prob.rlo_;
          double *rup = prob.rup_;
          double infinity = si.getInfinity();
          for (i = 0; i < prob.nrows_; i++) {
            status = prob.getRowStatus(i);
            if (status == CoinPrePostsolveMatrix::basic) {
              basicCnt--;
              double down = acts[i] - rlo[i];
              double up = rup[i] - acts[i];
              if (std::min(up, down) < infinity) {
                if (down <= up)
                  prob.setRowStatus(i, CoinPrePostsolveMatrix::atLowerBound);
                else
                  prob.setRowStatus(i, CoinPrePostsolveMatrix::atUpperBound);
              } else {
                prob.setRowStatus(i, CoinPrePostsolveMatrix::isFree);
              }
            }
            if (basicCnt == prob.nrows_)
              break;
          }
        }
      }
#endif

      /*
  Install the status and primal solution, if we've been carrying them along.

  The code that copies status is efficient but brittle. The current definitions
  for CoinWarmStartBasis::Status and CoinPrePostsolveMatrix::Status are in
  one-to-one correspondence. This code will fail if that ever changes.
*/
      if (doStatus) {
        presolvedModel_->setColSolution(prob.sol_);
        CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(presolvedModel_->getEmptyWarmStart());
        basis->resize(prob.nrows_, prob.ncols_);
        int i;
        for (i = 0; i < prob.ncols_; i++) {
          CoinWarmStartBasis::Status status = static_cast< CoinWarmStartBasis::Status >(prob.getColumnStatus(i));
          basis->setStructStatus(i, status);
        }
        for (i = 0; i < prob.nrows_; i++) {
          CoinWarmStartBasis::Status status = static_cast< CoinWarmStartBasis::Status >(prob.getRowStatus(i));
          basis->setArtifStatus(i, status);
        }
        presolvedModel_->setWarmStart(basis);
        delete basis;
        delete[] prob.sol_;
        delete[] prob.acts_;
        delete[] prob.colstat_;
        prob.sol_ = NULL;
        prob.acts_ = NULL;
        prob.colstat_ = NULL;
      }
      /*
  Copy original column and row information from the CoinPresolveMatrix object
  so it'll be available for postsolve.
*/
      int ncolsNow = presolvedModel_->getNumCols();
      memcpy(originalColumn_, prob.originalColumn_, ncolsNow * sizeof(int));
      delete[] prob.originalColumn_;
      prob.originalColumn_ = NULL;
      int nrowsNow = presolvedModel_->getNumRows();
      memcpy(originalRow_, prob.originalRow_, nrowsNow * sizeof(int));
      delete[] prob.originalRow_;
      prob.originalRow_ = NULL;

      // now clean up integer variables.  This can modify original
      {
        int numberChanges = 0;
        const double *lower0 = originalModel_->getColLower();
        const double *upper0 = originalModel_->getColUpper();
        const double *lower = presolvedModel_->getColLower();
        const double *upper = presolvedModel_->getColUpper();
        for (i = 0; i < ncolsNow; i++) {
          if (presolvedModel_->isInteger(i)) {
	    int iOriginal = originalColumn_[i];
	    double lowerValue0 = lower0[iOriginal];
	    double upperValue0 = upper0[iOriginal];
	    double lowerValue = ceil(lower[i] - 1.0e-5);
	    double upperValue = floor(upper[i] + 1.0e-5);
	    presolvedModel_->setColBounds(i, lowerValue, upperValue);
	    // need to be careful if dupcols
	    if (lowerValue > upperValue) {
	      numberChanges++;
	      CoinMessageHandler *hdlr = presolvedModel_->messageHandler();
	      hdlr->message(COIN_PRESOLVE_COLINFEAS, msgs)
		<< iOriginal << lowerValue << upperValue << CoinMessageEol;
	      result = 1;
	    } else if ((prob.presolveOptions_ & 0x80000000) == 0) {
	      if (lowerValue > lowerValue0 + 1.0e-8) {
		originalModel_->setColLower(iOriginal, lowerValue);
		numberChanges++;
	      }
	      if (upperValue < upperValue0 - 1.0e-8) {
		originalModel_->setColUpper(iOriginal, upperValue);
		numberChanges++;
	      }
	    }
	  } else if (scLower) {
	    int iOriginal = originalColumn_[i];
	    if (scLower[iOriginal]!= -COIN_DBL_MAX) {
	      double lowerSC = scLower[iOriginal];
	      double lowerValue = lower[i];
	      double upperValue = upper[i];
	      double lowerValue0 = lower0[iOriginal];
	      double upperValue0 = upper0[iOriginal];
	      //printf("SC %d %g <= %g orig %d %g <= %g sclo %g\n",
	      //     i,lowerValue,upperValue,iOriginal,lowerValue0,
	      //     upperValue0,lowerSC);
	      if (upperValue<lowerSC-1.0e-5) {
		lowerValue = 0.0;
		upperValue = 0.0;
	      } else if (lowerValue >1.0e-5) {
		lowerValue = std::max(lowerValue,lowerSC);
	      }
	      presolvedModel_->setColBounds(i, lowerValue, upperValue);
	    }
	  }
        }
        if (numberChanges) {
          CoinMessageHandler *hdlr = presolvedModel_->messageHandler();
          hdlr->message(COIN_PRESOLVE_INTEGERMODS, msgs)
            << numberChanges << CoinMessageEol;
          // we can't go round again in integer if dupcols
          if (!result && totalPasses > 0 && (prob.presolveOptions_ & 0x80000000) == 0) {
            result = -1; // round again
            const CoinPresolveAction *paction = paction_;
            while (paction) {
              const CoinPresolveAction *next = paction->next;
              delete paction;
              paction = next;
            }
            paction_ = NULL;
          }
        }
      }
    } else if (prob.status_ != 0) {
      // infeasible or unbounded
      result = 1;
    }
    delete probptr;
  }
  if (!result) {
    int nrowsAfter = presolvedModel_->getNumRows();
    int ncolsAfter = presolvedModel_->getNumCols();
    CoinBigIndex nelsAfter = presolvedModel_->getNumElements();
    CoinMessageHandler *hdlr = presolvedModel_->messageHandler();
    hdlr->message(COIN_PRESOLVE_STATS, msgs)
      << nrowsAfter << -(nrows_ - nrowsAfter)
      << ncolsAfter << -(ncols_ - ncolsAfter)
      << nelsAfter << -(nelems_ - nelsAfter) << CoinMessageEol;
#if DEBUG_PREPROCESS > 1
    if (debugSolution) {
      for (int i=0;i<ncolsAfter;i++) {
	int iColumn = originalColumn_[i];
	debugSolution[i] = debugSolution[iColumn];
      }
      debugNumberColumns = ncolsAfter;
    }
#endif
  } else {
    gutsOfDestroy();
    delete presolvedModel_;
    presolvedModel_ = NULL;
  }

  if (presolvedModel_) {
    lastPresolveRows = presolvedModel_->getNumRows();
    lastPresolveCols = presolvedModel_->getNumCols();
  } else {
    lastPresolveRows = 0;
    lastPresolveCols = 0;
  }

  return presolvedModel_;
}
OsiSolverInterface *
OsiPresolve::miniPresolvedModel(OsiSolverInterface &si,
  double feasibilityTolerance,
  bool keepIntegers,
  int numberPasses,
  const char *prohibited,
  bool doStatus,
  const char *rowProhibited)
{
  ncols_ = si.getNumCols();
  nrows_ = si.getNumRows();
  nelems_ = si.getNumElements();
  numberPasses_ = numberPasses;

  double maxmin = si.getObjSense();
  originalModel_ = &si;
  delete[] originalColumn_;
  originalColumn_ = new int[ncols_];
  delete[] originalRow_;
  originalRow_ = new int[nrows_];
  int i;
  for (i = 0; i < ncols_; i++)
    originalColumn_[i] = i;
  for (i = 0; i < nrows_; i++)
    originalRow_[i] = i;

  // result is 0 - okay, 1 infeasible, -1 go round again
  int result = -1;

  // User may have deleted - its their responsibility
  presolvedModel_ = NULL;
  // Messages
  CoinMessages msgs = CoinMessage(si.messages().language());

  // make new copy
  delete presolvedModel_;
  presolvedModel_ = si.clone();


  CoinPresolveMatrix* probptr = construct_CoinPresolveMatrix(ncols_,
							     maxmin,
							     presolvedModel_,
							     nrows_, nelems_, doStatus, nonLinearValue_, prohibited,
							     rowProhibited);
  CoinPresolveMatrix& prob(*probptr);
  // make sure row solution correct
  if (doStatus) {
    double *colels = prob.colels_;
    int *hrow = prob.hrow_;
    CoinBigIndex *mcstrt = prob.mcstrt_;
    int *hincol = prob.hincol_;
    int ncols = prob.ncols_;

    double *csol = prob.sol_;
    double *acts = prob.acts_;
    int nrows = prob.nrows_;

    int colx;

    memset(acts, 0, nrows * sizeof(double));

    for (colx = 0; colx < ncols; ++colx) {
      double solutionValue = csol[colx];
      for (CoinBigIndex i = mcstrt[colx]; i < mcstrt[colx] + hincol[colx]; ++i) {
	int row = hrow[i];
	double coeff = colels[i];
	acts[row] += solutionValue * coeff;
      }
    }
  }

  // move across feasibility tolerance
  prob.feasibilityTolerance_ = feasibilityTolerance;

  /*
    Do presolve. Allow for the possibility that presolve might be ineffective
    (i.e., we're feasible but no postsolve actions are queued.
  */
  paction_ = miniPresolve(&prob, &si);
  result = 0;
  // Get rid of useful arrays
  prob.deleteStuff();
  delete probptr;
  return presolvedModel_;
}
const CoinPresolveAction *OsiPresolve::miniPresolve(CoinPresolveMatrix *prob,
						    OsiSolverInterface * solver)
{
  const CoinPresolveAction * paction = NULL;

  prob->status_ = 0; // say feasible

  /*
  Fix variables before we get into the main transform loop.
*/
  paction = make_fixed(prob, paction);
  paction = testRedundant(prob, paction);


    /*
  Set [rows,cols]ToDo to process all rows & cols unless there are
  specific prohibitions.
*/
    prob->initColsToDo();
    prob->initRowsToDo();
    const CoinPresolveAction *const paction0 = paction;
    // this can also make E rows so do one bit here
    paction = remove_dual_action::presolve(prob, paction);
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    int nTotalFix = 0;
    while (paction) {
      const CoinPresolveAction *const pactionNext = paction->next;
      const make_fixed_action * fix = dynamic_cast<const make_fixed_action *>(paction);
      if (fix) {
	int nActions = fix->nActions();
	nTotalFix += nActions;
	bool fixToLower = fix->fixToLower();
	if (fixToLower) {
	  for (int i=0;i<nActions;i++) {
	    int iColumn = fix->fixVariable(i);
	    solver->setColUpper(iColumn,lower[iColumn]);
	  }
	} else {
	  for (int i=0;i<nActions;i++) {
	    int iColumn = fix->fixVariable(i);
	    solver->setColLower(iColumn,upper[iColumn]);
	  }
	}
      }
      delete paction;
      paction = pactionNext;
    }
    if (nTotalFix) {
      CoinMessageHandler *hdlr = prob->messageHandler();
      CoinMessages msgs = CoinMessage(prob->messages().language());
      char line[80];
      sprintf(line,"Mini-presolve fixed %d variables",nTotalFix);
      hdlr->message(COIN_GENERAL_INFO, msgs)
        << line << CoinMessageEol;
    }
    return (NULL);
}

// Return pointer to presolved model
OsiSolverInterface *
OsiPresolve::model() const
{
  return presolvedModel_;
}
// Return pointer to original model
OsiSolverInterface *
OsiPresolve::originalModel() const
{
  return originalModel_;
}

void OsiPresolve::postsolve(bool updateStatus)
{
  // Messages
  CoinMessages msgs = CoinMessage(presolvedModel_->messages().language());
  CoinMessageHandler *hdlr = presolvedModel_->messageHandler();
  if (!presolvedModel_->isProvenOptimal()) {
    hdlr->message(COIN_PRESOLVE_NONOPTIMAL, msgs) << CoinMessageEol;
  }

  // this is the size of the original problem
  const int ncols0 = ncols_;
  const int nrows0 = nrows_;
  const CoinBigIndex nelems0 = nelems_;

  // reality check
  assert(ncols0 == originalModel_->getNumCols());
  assert(nrows0 == originalModel_->getNumRows());

  // this is the reduced problem
  int ncols = presolvedModel_->getNumCols();
  int nrows = presolvedModel_->getNumRows();

  double *acts = new double[nrows0];
  double *sol = new double[ncols0];
  CoinZeroN(acts, nrows0);
  CoinZeroN(sol, ncols0);

  unsigned char *rowstat = NULL;
  unsigned char *colstat = NULL;
  CoinWarmStartBasis *presolvedBasis = dynamic_cast< CoinWarmStartBasis * >(presolvedModel_->getWarmStart());
  if (!presolvedBasis)
    updateStatus = false;
  if (updateStatus) {
    colstat = new unsigned char[ncols0 + nrows0];
#ifdef ZEROFAULT
    memset(colstat, 0, ((ncols0 + nrows0) * sizeof(char)));
#endif
    rowstat = colstat + ncols0;
    for (int i = 0; i < ncols; i++) {
      colstat[i] = presolvedBasis->getStructStatus(i);
    }
    for (int i = 0; i < nrows; i++) {
      rowstat[i] = presolvedBasis->getArtifStatus(i);
    }
  }
  delete presolvedBasis;

#if PRESOLVE_CONSISTENCY > 0
  if (updateStatus) {
    int basicCnt = 0;
    for (int i = 0; i < ncols; i++) {
      if (colstat[i] == CoinWarmStartBasis::basic)
        basicCnt++;
    }
    for (int i = 0; i < nrows; i++) {
      if (rowstat[i] == CoinWarmStartBasis::basic)
        basicCnt++;
    }

    assert(basicCnt == nrows);
  }
#endif

  /*
  Postsolve back to the original problem.  The CoinPostsolveMatrix object
  assumes ownership of sol, acts, colstat, and rowstat.
*/
  CoinPostsolveMatrix* probptr =
      construct_CoinPostsolveMatrix(presolvedModel_,ncols0,nrows0,nelems0,
                                    presolvedModel_->getObjSense(),
                                    sol,acts,colstat,rowstat) ;
  CoinPostsolveMatrix& prob(*probptr);
  probptr->originalRowLower_ = originalModel_->getRowLower();
  probptr->originalRowUpper_ = originalModel_->getRowUpper();
  postsolve(prob) ;

#if PRESOLVE_CONSISTENCY > 0
  if (updateStatus) {
    int basicCnt = 0;
    for (int i = 0; i < ncols0; i++) {
      if (prob.getColumnStatus(i) == CoinPrePostsolveMatrix::basic)
        basicCnt++;
    }
    for (int i = 0; i < nrows0; i++) {
      if (prob.getRowStatus(i) == CoinPrePostsolveMatrix::basic)
        basicCnt++;
    }

    assert(basicCnt == nrows0);
  }
#endif

  originalModel_->setColSolution(sol);
  if (updateStatus) {
    CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(presolvedModel_->getEmptyWarmStart());
    basis->setSize(ncols0, nrows0);
    const double *lower = originalModel_->getColLower();
    const double *upper = originalModel_->getColUpper();
    const double *solution = originalModel_->getColSolution();
    for (int i = 0; i < ncols0; i++) {
      CoinWarmStartBasis::Status status = static_cast< CoinWarmStartBasis::Status >(prob.getColumnStatus(i));
      // Fix obvious mistakes
      if (status != CoinWarmStartBasis::basic && status != CoinWarmStartBasis::isFree) {
        if (solution[i] < lower[i] + 1.0e-8)
          status = CoinWarmStartBasis::atLowerBound;
        else if (solution[i] > upper[i] - 1.0e-8)
          status = CoinWarmStartBasis::atUpperBound;
      }
      assert(status != CoinWarmStartBasis::atLowerBound || originalModel_->getColLower()[i] > -originalModel_->getInfinity());
      assert(status != CoinWarmStartBasis::atUpperBound || originalModel_->getColUpper()[i] < originalModel_->getInfinity());
      basis->setStructStatus(i, status);
    }

#if PRESOLVE_DEBUG > 0
    /*
    Do a thorough check of row and column solutions. There should be no
    inconsistencies at this point.
  */
    std::cout
      << "Checking solution before transferring basis." << std::endl;
    presolve_check_sol(&prob, 2, 2, 2);
    int errs = 0;
#endif
    for (int i = 0; i < nrows0; i++) {
      CoinWarmStartBasis::Status status = static_cast< CoinWarmStartBasis::Status >(prob.getRowStatus(i));
      basis->setArtifStatus(i, status);
    }
    originalModel_->setWarmStart(basis);
    delete basis;
  }
  delete probptr;
}

// return pointer to original columns
const int *
OsiPresolve::originalColumns() const
{
  return originalColumn_;
}
// return pointer to original rows
const int *
OsiPresolve::originalRows() const
{
  return originalRow_;
}
// Set pointer to original model
void OsiPresolve::setOriginalModel(OsiSolverInterface *model)
{
  originalModel_ = model;
}


#if PRESOLVE_DEBUG > 0
// Anonymous namespace for debug routines
namespace {

/*
  A control routine for debug checks --- keeps down the clutter in doPresolve.
  Each time it's called, it prints a list of transforms applied since the last
  call, then does checks.
*/

void check_and_tell(const CoinPresolveMatrix *const prob,
  const CoinPresolveAction *first,
  const CoinPresolveAction *&mark)

{
  const CoinPresolveAction *current;

  if (first != mark) {
    printf("PRESOLVE: applied");
    for (current = first;
         current != mark && current != 0;
         current = current->next) {
      printf(" %s", current->name());
    }
    printf("\n");

    presolve_check_sol(prob);
    presolve_check_nbasic(prob);
    mark = first;
  }

  return;
}
int counter = 1000000;

bool break2(CoinPresolveMatrix *prob)
{
  if (counter > 0)
    printf("break2: counter %d\n", counter);
  counter--;
#if DEBUG_PREPROCESS > 1
  if (debugSolution && prob->ncols_ == debugNumberColumns) {
    for (int i = 0; i < prob->ncols_; i++) {
      double value = debugSolution[i];
      if (value < prob->clo_[i]) {
        printf("%d inf %g %g %g\n", i, prob->clo_[i], value, prob->cup_[i]);
      } else if (value > prob->cup_[i]) {
        printf("%d inf %g %g %g\n", i, prob->clo_[i], value, prob->cup_[i]);
      }
    }
  }
#endif
  if (!counter) {
    printf("skipping next and all\n");
  }
  return (counter <= 0);
}

} // end anonymous namespace for debug routines
#endif

#if PRESOLVE_DEBUG > 0
#define possibleBreak \
  if (break2(prob))   \
  break
#define possibleSkip if (!break2(prob))
#else
#define possibleBreak
#define possibleSkip
#endif

// This is the presolve loop.
// It is a separate virtual function so that it can be easily
// customized by subclassing CoinPresolve.


#include "CoinPresolveMatrix.hpp"
#define TRACK_PRESOLVE(NAME, STMT) \
  do { \
    int rB=0, cB=0, nzB=0; \
    for (int _i=0; _i<prob->nrows_; _i++) if (prob->hinrow_[_i]>0) { rB++; nzB+=prob->hinrow_[_i]; } \
    for (int _j=0; _j<prob->ncols_; _j++) if (prob->hincol_[_j]>0) cB++; \
    STMT; \
    if (prob->status_==0) { \
      int rA=0, cA=0, nzA=0; \
      for (int _i=0; _i<prob->nrows_; _i++) if (prob->hinrow_[_i]>0) { rA++; nzA+=prob->hinrow_[_i]; } \
      for (int _j=0; _j<prob->ncols_; _j++) if (prob->hincol_[_j]>0) cA++; \
      if (rB!=rA || cB!=cA || nzB!=nzA) { \
        CoinAddPresolveStats(NAME, rB-rA, cB-cA, nzB-nzA); \
      } \
    } \
  } while(0)

const CoinPresolveAction *OsiPresolve::presolve(CoinPresolveMatrix *prob)
{
  paction_ = 0;

  prob->status_ = 0; // say feasible

#if PRESOLVE_DEBUG > 0
  const CoinPresolveAction *pactiond = 0;
  presolve_check_sol(prob, 2, 1, 1);

  // CoinPresolveMonitor *monitor = new CoinPresolveMonitor(prob,true,22) ;
  CoinPresolveMonitor *monitor = 0;
#endif
  /*
  Transfer costs off of singleton variables, and also between integer
  variables when advantageous.

  transferCosts is defined in CoinPresolveFixed.cpp
*/
  if ((presolveActions_ & 0x04) != 0) {
    transferCosts(prob);
#if PRESOLVE_DEBUG > 0
    if (monitor)
      monitor->checkAndTell(prob);
#endif
  }
  /*
  Fix variables before we get into the main transform loop.
*/
  TRACK_PRESOLVE("make_fixed", paction_ = make_fixed(prob, paction_););
  TRACK_PRESOLVE("testRedundant", paction_ = testRedundant(prob, paction_););

  // GF(2) parity reduction — run early, before the main loop
  if (prob->anyInteger() && !prob->status_) {
    TRACK_PRESOLVE("parity_action", paction_ = parity_action::presolve(prob, paction_););
    if (prob->status_)
      return (paction_);
  }

#if PRESOLVE_DEBUG > 0
  check_and_tell(prob, paction_, pactiond);
  if (monitor)
    monitor->checkAndTell(prob);
#endif

  // if integers then switch off dual stuff
  // later just do individually
  bool doDualStuff = true;
  if ((presolveActions_ & 0x01) == 0) {
    int ncol = presolvedModel_->getNumCols();
    for (int i = 0; i < ncol; i++)
      if (presolvedModel_->isInteger(i))
        doDualStuff = false;
  }

#if PRESOLVE_CONSISTENCY > 0
  presolve_links_ok(prob);
#endif

  /*
  If we're feasible, set up for the main presolve transform loop.
*/
  if (!prob->status_) {
#if 1
    // normal operation --- all transforms enabled
    bool slackSingleton = true;
    bool slackd = true;
    bool doubleton = true;
    bool tripleton = true;
    bool forcing = true;
    bool ifree = true;
    bool zerocost = true;
    bool dupcol = true;
    bool duprow = true;
    bool dual = doDualStuff;
#else
    // compile time selection of transforms.
    bool slackSingleton = true;
    bool slackd = false;
    bool doubleton = true;
    bool tripleton = true;
    bool forcing = true;
    bool ifree = false;
    bool zerocost = false;
    bool dupcol = false;
    bool duprow = false;
    bool dual = false;
#endif
    /*
  Process OsiPresolve options. Set corresponding CoinPresolve options and
  control variables here.
*/
    // Switch off some stuff if would annoy set partitioning etc
    if ((presolveActions_ & 0x02) != 0) {
      doubleton = false;
      tripleton = false;
      ifree = false;
    }
    // stop x+y+z=1
    if ((presolveActions_ & 0x08) != 0)
      prob->setPresolveOptions(prob->presolveOptions() | 0x04);
    // switch on stuff which can't be unrolled easily
    if ((presolveActions_ & 0x10) != 0)
      prob->setPresolveOptions(prob->presolveOptions() | 0x10);
    // switch on testing stuff for dual
    if ((presolveActions_ & 0x4000) != 0)
      prob->setPresolveOptions(prob->presolveOptions() | 0x200000);
    // switch on gub stuff (unimplemented as of 110605 -- lh --)
    if ((presolveActions_ & 0x20) != 0)
      prob->setPresolveOptions(prob->presolveOptions() | 0x20);
    // allow duplicate column processing for integer columns
    if ((presolveActions_ & 0x01) != 0)
      prob->setPresolveOptions(prob->presolveOptions() | 0x01);
    // allow extra duplicate column processing
    if ((presolveActions_ & 0x40) != 0)
      prob->setPresolveOptions(prob->presolveOptions() | 0x80000);
    // allow extra doubleton column processing
    if ((presolveActions_ & 0x80) != 0)
      prob->setPresolveOptions(prob->presolveOptions() | 0x100000);
    /*
  Set [rows,cols]ToDo to process all rows & cols unless there are
  specific prohibitions.
*/
    prob->initColsToDo();
    prob->initRowsToDo();
    /*
  Try to remove duplicate rows and columns.
*/
    if (dupcol) {
      possibleSkip;
      TRACK_PRESOLVE("dupcol_action", paction_ = dupcol_action::presolve(prob, paction_););
      possibleSkip;
#ifdef CBC_PREPROCESS_EXPERIMENT
      //paction_ = twoxtwo_action::presolve(prob, paction_);
#endif
#if PRESOLVE_DEBUG > 0
      if (monitor)
        monitor->checkAndTell(prob);
#endif
    }
    if (duprow) {
      possibleSkip;
      TRACK_PRESOLVE("duprow_action", paction_ = duprow_action::presolve(prob, paction_););
#if PRESOLVE_DEBUG > 0
      if (monitor)
        monitor->checkAndTell(prob);
#endif
    }
    /*
  The main loop starts with a minor loop that does inexpensive presolve
  transforms until convergence. At each iteration of this loop,
  next[Rows,Cols]ToDo is copied over to [rows,cols]ToDo.

  Then there's a block to set [rows,cols]ToDo to examine all rows & cols,
  followed by executions of expensive transforms. Then we come back around for
  another iteration of the main loop. [rows,cols]ToDo is not reset as we come
  back around, so we dive into the inexpensive loop set up to process all.

  lastDropped is a count of total number of rows dropped by presolve. Used as
  an additional criterion to end the main presolve loop.
*/
    int lastDropped = 0;
    prob->pass_ = 0;
    for (int iLoop = 0; iLoop < numberPasses_; iLoop++) {

#if PRESOLVE_SUMMARY > 0
      std::cout << "Starting major pass " << (iLoop + 1) << std::endl;
#endif

      const CoinPresolveAction *const paction0 = paction_;
// #define IMPLIED 3
#ifdef IMPLIED
      int fill_level = 3;
#define IMPLIED2 1
#if IMPLIED != 3
#if IMPLIED > 0 && IMPLIED < 11
      fill_level = IMPLIED;
      printf("** fill_level == %d !\n", fill_level);
#endif
#if IMPLIED > 11 && IMPLIED < 21
      fill_level = -(IMPLIED - 10);
      printf("** fill_level == %d !\n", fill_level);
#endif
#endif
#else
      // look for substitutions with little fill
      int fill_level = 2;
      if ((presolveActions_&0x300) != 0) {
	fill_level += (presolveActions_&0x300)>>8;
      }
#endif
      int whichPass = 0;
      /*
  Apply inexpensive transforms until convergence or infeasible/unbounded.
*/
      while (true) {
        whichPass++;
        prob->pass_++;
        const CoinPresolveAction *const paction1 = paction_;

        if (slackd) {
          bool notFinished = true;
          while (notFinished) {
            possibleBreak;
            TRACK_PRESOLVE("slack_doubleton_action", paction_ = slack_doubleton_action::presolve(prob, paction_, notFinished););
          }
#if PRESOLVE_DEBUG > 0
          check_and_tell(prob, paction_, pactiond);
          if (monitor)
            monitor->checkAndTell(prob);
#endif
          if (prob->status_)
            break;
        }

        if (zerocost) {
          possibleBreak;
          TRACK_PRESOLVE("do_tighten_action", paction_ = do_tighten_action::presolve(prob, paction_););
#if PRESOLVE_DEBUG > 0
          check_and_tell(prob, paction_, pactiond);
          if (monitor)
            monitor->checkAndTell(prob);
#endif
          if (prob->status_)
            break;
        }

        if (dual && whichPass == 1) {
          possibleBreak;
          // this can also make E rows so do one bit here
          TRACK_PRESOLVE("remove_dual_action", paction_ = remove_dual_action::presolve(prob, paction_););
#if PRESOLVE_DEBUG > 0
          check_and_tell(prob, paction_, pactiond);
          if (monitor)
            monitor->checkAndTell(prob);
#endif
          if (prob->status_)
            break;
        }

        if (doubleton) {
          possibleBreak;
          TRACK_PRESOLVE("doubleton_action", paction_ = doubleton_action::presolve(prob, paction_););
#if PRESOLVE_DEBUG > 0
          check_and_tell(prob, paction_, pactiond);
          if (monitor)
            monitor->checkAndTell(prob);
#endif
          if (prob->status_)
            break;
        }

        if (tripleton) {
          possibleBreak;
          TRACK_PRESOLVE("tripleton_action", paction_ = tripleton_action::presolve(prob, paction_););
#if PRESOLVE_DEBUG > 0
          check_and_tell(prob, paction_, pactiond);
          if (monitor)
            monitor->checkAndTell(prob);
#endif
          if (prob->status_)
            break;
        }

        if (forcing) {
          possibleBreak;
          TRACK_PRESOLVE("forcing_constraint_action", paction_ = forcing_constraint_action::presolve(prob, paction_););
#if PRESOLVE_DEBUG > 0
          check_and_tell(prob, paction_, pactiond);
          if (monitor)
            monitor->checkAndTell(prob);
#endif
          if (prob->status_)
            break;
        }

        if (ifree && (whichPass % 5) == 1) {
          possibleBreak;
          TRACK_PRESOLVE("implied_free_action", paction_ = implied_free_action::presolve(prob, paction_, fill_level););
#if PRESOLVE_DEBUG > 0
          check_and_tell(prob, paction_, pactiond);
          if (monitor)
            monitor->checkAndTell(prob);
#endif
          if (prob->status_)
            break;
        }

#if PRESOLVE_CONSISTENCY > 0
        presolve_links_ok(prob);
        presolve_no_zeros(prob);
        presolve_consistent(prob);
#endif
        /*
  Set up for next pass.

  Original comment adds: later do faster if many changes i.e. memset and memcpy
*/
        prob->stepRowsToDo();

#if PRESOLVE_DEBUG > 0
        int rowCheck = -1;
        bool rowFound = false;
        for (int i = 0; i < prob->numberRowsToDo_; i++) {
          int index = prob->rowsToDo_[i];
          if (index == rowCheck) {
            std::cout
              << "  row " << index << " on list after pass " << whichPass
              << std::endl;
            rowFound = true;
          }
        }
        if (!rowFound && rowCheck >= 0)
          prob->rowsToDo_[prob->numberRowsToDo_++] = rowCheck;
#endif

        prob->stepColsToDo();

#if PRESOLVE_DEBUG > 0
        int colCheck = -1;
        bool colFound = false;
        for (int i = 0; i < prob->numberNextColsToDo_; i++) {
          int index = prob->colsToDo_[i];
          if (index == colCheck) {
            std::cout
              << "  col " << index << " on list after pass " << whichPass
              << std::endl;
            colFound = true;
          }
        }
        if (!colFound && colCheck >= 0)
          prob->colsToDo_[prob->numberColsToDo_++] = colCheck;
#endif
        /*
  Break if nothing happened (no postsolve actions queued).

  The check for fill_level > 0 is a hack to allow repeating the loop with some
  modified fill level (playing with negative values).

  fill_level = 0 (as set in other places) will clearly be a problem.
  -- lh, 110605 --
*/
        if (paction_ == paction1 && fill_level > 0)
          break;
      }
      /*
  End of inexpensive transform loop.
  Reset [rows,cols]ToDo to process all rows and columns unless there are
  specfic prohibitions.
*/
      prob->initRowsToDo();
      prob->initColsToDo();
/*
  Try expensive presolve transforms.

  Original comment adds: this caused world.mps to run into numerical
  			 difficulties
*/
#if PRESOLVE_SUMMARY > 0
      std::cout << "Starting expensive." << std::endl;
#endif
      /*
  Try and fix variables at upper or lower bound by calculating bounds on the
  dual variables and propagating them to the reduced costs. Every other
  iteration, see if this has created free variables.
*/
      if (dual) {
        for (int itry = 0; itry < 5; itry++) {
          const CoinPresolveAction *const paction2 = paction_;
          possibleBreak;
          TRACK_PRESOLVE("remove_dual_action", paction_ = remove_dual_action::presolve(prob, paction_););
#if PRESOLVE_DEBUG > 0
          check_and_tell(prob, paction_, pactiond);
          if (monitor)
            monitor->checkAndTell(prob);
#endif
          if (prob->status_)
            break;
          if (ifree) {
#ifdef IMPLIED
#if IMPLIED2 == 0
            int fill_level = 0; // switches off substitution
#elif IMPLIED2 != 99
            int fill_level = IMPLIED2;
#endif
#endif
            if ((itry & 1) == 0) {
              possibleBreak;
              TRACK_PRESOLVE("implied_free_action", paction_ = implied_free_action::presolve(prob, paction_, fill_level););
            }
#if PRESOLVE_DEBUG > 0
            check_and_tell(prob, paction_, pactiond);
            if (monitor)
              monitor->checkAndTell(prob);
#endif
            if (prob->status_)
              break;
          }
          if (paction_ == paction2)
            break;
        }
      } else if (ifree) {
/*
  Just check for free variables.
*/
#ifdef IMPLIED
#if IMPLIED2 == 0
        int fill_level = 0; // switches off substitution
#elif IMPLIED2 != 99
        int fill_level = IMPLIED2;
#endif
#endif
        possibleBreak;
        TRACK_PRESOLVE("implied_free_action", paction_ = implied_free_action::presolve(prob, paction_, fill_level););
#if PRESOLVE_DEBUG > 0
        check_and_tell(prob, paction_, pactiond);
        if (monitor)
          monitor->checkAndTell(prob);
#endif
        if (prob->status_)
          break;
      }
      /*
  Check if other transformations have produced duplicate rows or columns.
*/
      if (dupcol) {
        possibleBreak;
        TRACK_PRESOLVE("dupcol_action", paction_ = dupcol_action::presolve(prob, paction_););
#if PRESOLVE_DEBUG > 0
        check_and_tell(prob, paction_, pactiond);
        if (monitor)
          monitor->checkAndTell(prob);
#endif
        if (prob->status_)
          break;
      }
      if (duprow) {
        possibleBreak;
        TRACK_PRESOLVE("duprow_action", paction_ = duprow_action::presolve(prob, paction_););
#if PRESOLVE_DEBUG > 0
        check_and_tell(prob, paction_, pactiond);
        if (monitor)
          monitor->checkAndTell(prob);
#endif
        if (prob->status_)
          break;
      }

      if (prob->anyInteger()) {
        possibleBreak;
        TRACK_PRESOLVE("parity_action", paction_ = parity_action::presolve(prob, paction_););
#if PRESOLVE_DEBUG > 0
        check_and_tell(prob, paction_, pactiond);
        if (monitor)
          monitor->checkAndTell(prob);
#endif
        if (prob->status_)
          break;
      }
      // Will trigger abort due to unimplemented postsolve  -- lh, 110605 --
      if ((presolveActions_ & 0x20) != 0) {
        possibleBreak;
        TRACK_PRESOLVE("gubrow_action", paction_ = gubrow_action::presolve(prob, paction_););
      }
      /*
  Count the number of empty rows and see if we've made progress in this pass.
*/
      bool stopLoop = false;
      {
        const int *const hinrow = prob->hinrow_;
        int numberDropped = 0;
        for (int i = 0; i < nrows_; i++)
          if (!hinrow[i])
            numberDropped++;
#if PRESOLVE_DEBUG > 0
        std::cout
          << "  " << (numberDropped - lastDropped)
          << " rows dropped in pass " << iLoop << "." << std::endl;
#endif
        if (numberDropped == lastDropped)
          stopLoop = true;
        else
          lastDropped = numberDropped;
      }
      /*
  Check for singleton variables that can act like a logical, allowing a
  row to be transformed from an equality to an inequality.

  The third parameter allows for costs for the existing logicals. This
  is apparently used by clp; consult the clp presolve before implementing
  it here.  -- lh, 110605 --

  Original comment: Do this here as not very loopy
*/
      if (slackSingleton) {
        possibleBreak;
        TRACK_PRESOLVE("slack_singleton_action", paction_ = slack_singleton_action::presolve(prob, paction_, NULL););
#if PRESOLVE_DEBUG > 0
        check_and_tell(prob, paction_, pactiond);
        if (monitor)
          monitor->checkAndTell(prob);
#endif
      }
#if PRESOLVE_DEBUG > 0
      presolve_check_sol(prob, 1);
#endif

      if (paction_ == paction0 || stopLoop)
        break;

    } // End of major pass loop
  }
  if (!prob->status_) {
    TRACK_PRESOLVE("duprow3_action", paction_ = duprow3_action::presolve(prob, paction_););
  }
  /*
  Final cleanup: drop zero coefficients from the matrix, then drop empty rows
  and columns.
*/
  if (!prob->status_) {
    paction_ = drop_zero_coefficients(prob, paction_);
#if PRESOLVE_DEBUG > 0
    check_and_tell(prob, paction_, pactiond);
    if (monitor)
      monitor->checkAndTell(prob);
#endif

    TRACK_PRESOLVE("drop_empty_cols_action", paction_ = drop_empty_cols_action::presolve(prob, paction_););
#if PRESOLVE_DEBUG > 0
    check_and_tell(prob, paction_, pactiond);
#endif

    TRACK_PRESOLVE("drop_empty_rows_action", paction_ = drop_empty_rows_action::presolve(prob, paction_););
#if PRESOLVE_DEBUG > 0
    check_and_tell(prob, paction_, pactiond);
#endif
  }
  /*
  Not feasible? Say something and clean up.
*/
  CoinMessageHandler *hdlr = prob->messageHandler();
  CoinMessages msgs = CoinMessage(prob->messages().language());
  if (prob->status_) {
    if (prob->status_ == 1)
      hdlr->message(COIN_PRESOLVE_INFEAS, msgs)
        << prob->feasibilityTolerance_ << CoinMessageEol;
    else if (prob->status_ == 2)
      hdlr->message(COIN_PRESOLVE_UNBOUND, msgs) << CoinMessageEol;
    else
      hdlr->message(COIN_PRESOLVE_INFEASUNBOUND, msgs) << CoinMessageEol;
    gutsOfDestroy();
  }
  return (paction_);
}

/*
  We could have implemented this by having each postsolve routine directly
  call the next one, but this makes it easier to add debugging checks.
*/
void OsiPresolve::postsolve(CoinPostsolveMatrix &prob)
{
  const CoinPresolveAction *paction = paction_;

#if PRESOLVE_DEBUG > 0
  std::cout << "Begin POSTSOLVING." << std::endl;
  if (prob.colstat_) {
    presolve_check_nbasic(&prob);
    presolve_check_sol(&prob, 2, 2, 2);
  }
  presolve_check_duals(&prob);
#endif

  while (paction) {
#if PRESOLVE_DEBUG > 0
    std::cout << "POSTSOLVING " << paction->name() << std::endl;
#endif

    paction->postsolve(&prob);

#if PRESOLVE_DEBUG > 0
    if (prob.colstat_) {
      presolve_check_nbasic(&prob);
      presolve_check_sol(&prob, 2, 2, 2);
    }
#endif
    paction = paction->next;
#if PRESOLVE_DEBUG > 0
    presolve_check_duals(&prob);
#endif
  }
#if PRESOLVE_DEBUG > 0
  std::cout << "End POSTSOLVING" << std::endl;
#endif

#if PRESOLVE_DEBUG > 0
  for (int j = 0; j < prob.ncols_; j++) {
    if (!prob.cdone_[j]) {
      printf("!cdone[%d]\n", j);
      abort();
    }
  }
  for (int i = 0; i < prob.nrows_; i++) {
    if (!prob.rdone_[i]) {
      printf("!rdone[%d]\n", i);
      abort();
    }
  }
  for (int j = 0; j < prob.ncols_; j++) {
    if (prob.sol_[j] < -1e10 || prob.sol_[j] > 1e10)
      printf("!!!%d %g\n", j, prob.sol_[j]);
  }
#endif

  /*
    Put back duals. Flip sign for maximisation problems.
  */
  double maxmin = originalModel_->getObjSense();
  if (maxmin < 0.0) {
    double *pi = prob.rowduals_;
    for (int i = 0; i < nrows_; i++)
      pi[i] = -pi[i];
  }
  originalModel_->setRowPrice(prob.rowduals_);
}

static inline double getTolerance(const OsiSolverInterface *si, OsiDblParam key)
{
  double tol;
  if (!si->getDblParam(key, tol)) {
    CoinPresolveAction::throwCoinError("getDblParam failed",
      "CoinPrePostsolveMatrix::CoinPrePostsolveMatrix");
  }
  return (tol);
}

// Assumptions:
// 1. nrows>=m.getNumRows()
// 2. ncols>=m.getNumCols()
//
// In presolve, these values are equal.
// In postsolve, they may be inequal, since the reduced problem
// may be smaller, but we need room for the large problem.
// ncols may be larger than si.getNumCols() in postsolve,
// this at that point si will be the reduced problem,
// but we need to reserve enough space for the original problem.
template <class CoinTsolveMatrix>
CoinTsolveMatrix* construct_CoinPrePostsolveMatrix(const OsiSolverInterface * si,
					     int ncols_in,
					     int nrows_in,
					     CoinBigIndex nelems_in,
                    CoinTsolveMatrix* dummy)
{
  CoinTsolveMatrix* cpm = new CoinTsolveMatrix(ncols_in,nrows_in,si->getNumElements());

  cpm->ncols_ = si->getNumCols();
  cpm->nelems_  = si->getNumElements();

  cpm->mcstrt_ = new CoinBigIndex[ncols_in+1];
  cpm->hincol_ = new int[ncols_in+1];

  cpm->cost_ = new double[ncols_in];
  cpm->clo_ = new double[ncols_in];
  cpm->cup_ = new double[ncols_in];
  cpm->rlo_ = new double[nrows_in];
  cpm->rup_ = new double[nrows_in];
  cpm->originalColumn_ = new int[ncols_in];
  cpm->originalRow_ = new int[nrows_in];

  cpm->ztolzb_ = getTolerance(si, OsiPrimalTolerance);
  cpm->ztoldj_ = getTolerance(si, OsiDualTolerance);

  cpm->maxmin_ = si->getObjSense();

  cpm->bulk0_ = static_cast<CoinBigIndex>(cpm->bulkRatio_*nelems_in+ncols_in) ;
  cpm->hrow_ = new int [cpm->bulk0_+ncols_in] ;
  cpm->colels_ = new double[cpm->bulk0_+ncols_in] ;

  si->getDblParam(OsiObjOffset,cpm->originalOffset_);
  int ncols = si->getNumCols();
  int nrows = si->getNumRows();

  cpm->setMessageHandler(si->messageHandler()) ;

  CoinDisjointCopyN(si->getColLower(), ncols, cpm->clo_);
  CoinDisjointCopyN(si->getColUpper(), ncols, cpm->cup_);
  CoinDisjointCopyN(si->getObjCoefficients(), ncols, cpm->cost_);
  CoinDisjointCopyN(si->getRowLower(), nrows,  cpm->rlo_);
  CoinDisjointCopyN(si->getRowUpper(), nrows,  cpm->rup_);
  int i;
  // initialize and clean up bounds
  double infinity = si->getInfinity();
  if (infinity!=COIN_DBL_MAX) {
    for (i=0;i<ncols;i++) {
      if (cpm->clo_[i]==-infinity)
        cpm->clo_[i]=-COIN_DBL_MAX;
      if (cpm->cup_[i]==infinity)
        cpm->cup_[i]=COIN_DBL_MAX;
    }
    for (i=0;i<nrows;i++) {
      if (cpm->rlo_[i]==-infinity)
        cpm->rlo_[i]=-COIN_DBL_MAX;
      if (cpm->rup_[i]==infinity)
        cpm->rup_[i]=COIN_DBL_MAX;
    }
  }
  for (i=0;i<ncols_in;i++)
    cpm->originalColumn_[i]=i;
  for (i=0;i<nrows_in;i++)
    cpm->originalRow_[i]=i;
  cpm->sol_=NULL;
  cpm->rowduals_=NULL;
  cpm->acts_=NULL;

  cpm->rcosts_=NULL;
  cpm->colstat_=NULL;
  cpm->rowstat_=NULL;

  return cpm;
}

// I am not familiar enough with CoinPackedMatrix to be confident
// that I will implement a row-ordered version of toColumnOrderedGapFree
// properly.
static bool isGapFree(const CoinPackedMatrix &matrix)
{
  const CoinBigIndex *start = matrix.getVectorStarts();
  const int *length = matrix.getVectorLengths();
  int i;
  for (i = matrix.getSizeVectorLengths() - 1; i >= 0; --i) {
    if (start[i + 1] - start[i] != length[i])
      break;
  }
  return (!(i >= 0));
}

OSILIB_EXPORT
CoinPresolveMatrix* construct_CoinPresolveMatrix(int ncols0_in,
				       double maxmin,
				       // end prepost members
				       OsiSolverInterface *si,
				       // rowrep
				       int nrows_in,
				       CoinBigIndex nelems_in,
				       bool doStatus,
				       double nonLinearValue,
               const char *prohibited,
				       const char *rowProhibited)
{
  CoinPresolveMatrix* cpm = construct_CoinPrePostsolveMatrix(si,ncols0_in,nrows_in,nelems_in,(CoinPresolveMatrix*)NULL);
  cpm->clink_ = new presolvehlink[ncols0_in+1];
  cpm->rlink_ = new presolvehlink[nrows_in+1];

  // temporary init
  cpm->mrstrt_ = new CoinBigIndex[nrows_in+1];
  cpm->hinrow_ = new int[nrows_in+1];
  cpm->integerType_ = new unsigned char[ncols0_in];
  cpm->colsToDo_ = new int [ncols0_in];
  cpm->nextColsToDo_ = new int[ncols0_in];
  cpm->rowsToDo_ = new int [nrows_in];
  cpm->nextRowsToDo_ = new int[nrows_in];

  cpm->rowels_ = new double [cpm->bulk0_] ;
  cpm->hcol_ = new int [cpm->bulk0_] ;
  cpm->nrows_ = si->getNumRows() ;
  const CoinBigIndex bufsize = static_cast<CoinBigIndex>(cpm->bulk0_) ;

  // Set up change bits
  cpm->rowChanged_ = new unsigned char[cpm->nrows_];
  memset(cpm->rowChanged_,0,cpm->nrows_);
  cpm->colChanged_ = new unsigned char[cpm->ncols_];
  memset(cpm->colChanged_,0,cpm->ncols_);
  const CoinPackedMatrix * m1 = si->getMatrixByCol();

  // The coefficient matrix is a big hunk of stuff.
  // Do the copy here to try to avoid running out of memory.

  const CoinBigIndex *start = m1->getVectorStarts();
  const int *length = m1->getVectorLengths();
  const int *row = m1->getIndices();
  const double *element = m1->getElements();
  int icol;
  CoinBigIndex nel=0;
  cpm->mcstrt_[0]=0;
  for (icol=0;icol<cpm->ncols_;icol++) {
    CoinBigIndex j;
    for (j=start[icol];j<start[icol]+length[icol];j++) {
      if (fabs(element[j])>ZTOLDP) {
        cpm->hrow_[nel]=row[j];
        cpm->colels_[nel++]=element[j];
      }
    }
    cpm->hincol_[icol]=static_cast<int>(nel-cpm->mcstrt_[icol]);
    cpm->mcstrt_[icol+1]=nel;
  }

  // same thing for row rep
  CoinPackedMatrix *m = new CoinPackedMatrix();
  m->reverseOrderedCopyOf(*si->getMatrixByCol());
  // do by hand because of zeros m->removeGaps();
  CoinDisjointCopyN(m->getVectorStarts(),  cpm->nrows_,  cpm->mrstrt_);
  cpm->mrstrt_[cpm->nrows_] = cpm->nelems_;
  CoinDisjointCopyN(m->getVectorLengths(), cpm->nrows_,  cpm->hinrow_);
  CoinDisjointCopyN(m->getIndices(),       cpm->nelems_, cpm->hcol_);
  CoinDisjointCopyN(m->getElements(),      cpm->nelems_, cpm->rowels_);
  start = m->getVectorStarts();
  length = m->getVectorLengths();
  const int *column = m->getIndices();
  element = m->getElements();
  // out zeros
  int irow;
  nel=0;
  cpm->mrstrt_[0]=0;
  for (irow=0;irow<cpm->nrows_;irow++) {
    CoinBigIndex j;
    for (j=start[irow];j<start[irow]+length[irow];j++) {
      if (fabs(element[j])>ZTOLDP) {
        cpm->hcol_[nel]=column[j];
        cpm->rowels_[nel++]=element[j];
      }
    }
    cpm->hinrow_[irow]=static_cast<int>(nel-cpm->mrstrt_[irow]);
    cpm->mrstrt_[irow+1]=nel;
  }
  cpm->nelems_=nel;

  delete m;
  {
    int i;
    int numberIntegers=0;
    for (i=0;i<cpm->ncols_;i++) {
      if (si->isInteger(i)) {
        cpm->integerType_[i] = 1;
        numberIntegers++;
      } else {
        cpm->integerType_[i] = 0;
      }
    }
    cpm->anyInteger_ = (numberIntegers!=0);
  }

  // Set up prohibited bits if needed
  if (nonLinearValue) {
    cpm->anyProhibited_ = true;
    for (icol=0;icol<cpm->ncols_;icol++) {
      CoinBigIndex j;
      bool nonLinearColumn = false;
      if (cpm->cost_[icol]==nonLinearValue)
      nonLinearColumn=true;
      for (j=cpm->mcstrt_[icol];j<cpm->mcstrt_[icol+1];j++) {
	if (cpm->colels_[j]==nonLinearValue) {
	  nonLinearColumn=true;
	  cpm->setRowProhibited(cpm->hrow_[j]);
	}
      }
      if (nonLinearColumn)
	cpm->setColProhibited(icol);
    }
  } else if (prohibited) {
    cpm->anyProhibited_ = true;
    for (icol=0;icol<cpm->ncols_;icol++) {
      if (prohibited[icol])
        cpm->setColProhibited(icol);
      if (prohibited[icol]==2)
        cpm->setColLeaveTotallyAlone(icol);
    }
  } else {
    cpm->anyProhibited_ = false;
  }
  // Any rows special?
  if (rowProhibited) {
    cpm->anyProhibited_ = true;
    for (int irow=0;irow<cpm->nrows_;irow++) {
      if (rowProhibited[irow])
        cpm->setRowProhibited(irow);
    }
  }
  // Go to minimization
  if (maxmin<0.0) {
    for (int i=0;i<cpm->ncols_;i++)
      cpm->cost_[i]=-cpm->cost_[i];
    cpm->maxmin_=1.0;
  }
  if (doStatus) {
    // allow for status and solution
    cpm->sol_ = new double[cpm->ncols_];
    const double *presol ;
    presol = si->getColSolution() ;
    memcpy(cpm->sol_,presol,cpm->ncols_*sizeof(double));;
    cpm->acts_ = new double [cpm->nrows_];
    memcpy(cpm->acts_,si->getRowActivity(),cpm->nrows_*sizeof(double));
    CoinWarmStartBasis * basis  =
    dynamic_cast<CoinWarmStartBasis*>(si->getWarmStart());
    cpm->colstat_ = new unsigned char [cpm->nrows_+cpm->ncols_];
    cpm->rowstat_ = cpm->colstat_+cpm->ncols_;
    // If basis is NULL then put in all slack basis
    if (basis&&basis->getNumStructural()==cpm->ncols_) {
      int i;
      for (i=0;i<cpm->ncols_;i++) {
        cpm->colstat_[i] = basis->getStructStatus(i);
      }
      for (i=0;i<cpm->nrows_;i++) {
        cpm->rowstat_[i] = basis->getArtifStatus(i);
      }
    } else {
      int i;
      // no basis
      for (i=0;i<cpm->ncols_;i++) {
        cpm->colstat_[i] = 3;
      }
      for (i=0;i<cpm->nrows_;i++) {
        cpm->rowstat_[i] = 1;
      }
    }
    delete basis;
  }


/*
  For building against CoinUtils 2.6, this #if 1 need to be changed into an
  #if 0
*/
  presolve_make_memlists(/*mcstrt_,*/ cpm->hincol_, cpm->clink_, cpm->ncols_);
  presolve_make_memlists(/*mrstrt_,*/ cpm->hinrow_, cpm->rlink_, cpm->nrows_);

  // this allows last col/row to expand up to bufsize-1 (22);
  // this must come after the calls to presolve_prefix
  cpm->mcstrt_[cpm->ncols_] = bufsize-1;
  cpm->mrstrt_[cpm->nrows_] = bufsize-1;
  // Allocate useful arrays
  cpm->initializeStuff();

#if PRESOLVE_CONSISTENCY > 0
  presolve_consistent(cpm);
#endif

  return cpm;
}

OSILIB_EXPORT
void update_model_CoinPresolveMatrix(CoinPresolveMatrix& cpm,
              OsiSolverInterface *si,
				      int nrows0, int ncols0,
				      CoinBigIndex nelems0)
{
  int nels = 0;
  int i;
  if (si->getObjSense() < 0.0) {
    for (int i=0;i<cpm.ncols_;i++)
      cpm.cost_[i]=-cpm.cost_[i];
    cpm.dobias_=-cpm.dobias_;
    cpm.maxmin_ = -1.0;
  }
  for ( i=0; i<cpm.ncols_; i++)
    nels += cpm.hincol_[i];
  CoinPackedMatrix m(true,cpm.nrows_,cpm.ncols_,nels,
                     cpm.colels_,cpm.hrow_,cpm.mcstrt_,cpm.hincol_);
  si->loadProblem(m, cpm.clo_, cpm.cup_, cpm.cost_, cpm.rlo_, cpm.rup_);

  for ( i=0; i<cpm.ncols_; i++) {
    if (cpm.integerType_[i])
      si->setInteger(i);
    else
      si->setContinuous(i);
  }
  si->setDblParam(OsiObjOffset,cpm.originalOffset_-cpm.dobias_);

#if PRESOLVE_SUMMARY > 0
  std::cout
    << "New ncol/nrow/nels: "
    << cpm.ncols_ << "(-" << ncols0-cpm.ncols_ << ") "
    << cpm.nrows_ << "(-" << nrows0-cpm.nrows_ << ") "
    << si->getNumElements() << "(-" << nelems0-si->getNumElements() << ") "
    << std::endl ;
# endif
}

////////////////  POSTSOLVE
OSILIB_EXPORT
CoinPostsolveMatrix* construct_CoinPostsolveMatrix(OsiSolverInterface*  si,
				       int ncols0_in,
				       int nrows0_in,
				       CoinBigIndex nelems0,

				       double maxmin,
				       // end prepost members

				       double *sol_in,
				       double *acts_in,

				       unsigned char *colstat_in,
				       unsigned char *rowstat_in)
{
  CoinPostsolveMatrix* cpm =
      construct_CoinPrePostsolveMatrix(si,ncols0_in,nrows0_in,nelems0,
                                       (CoinPostsolveMatrix*)NULL);
/*
  Used only to mark processed columns and rows so that debugging routines know
  what to check.
*/
# if PRESOLVE_DEBUG > 0 || PRESOLVE_CONSISTENCY > 0
  cpm->cdone_ = new char[ncols0_in];
  cpm->rdone_ = new char[nrows0_in];
# endif

  /*
  The CoinPrePostsolveMatrix constructor will set bulk0_ to bulkRatio_*nelems0.
  By default, bulkRatio_ is 2. This is certainly larger than absolutely
  necessary, but good for efficiency (minimises the need to compress the bulk
  store). The main storage arrays for the threaded column-major representation
  (hrow_, colels_, link_) should be allocated to this size.
*/
  cpm->maxlink_ = cpm->bulk0_ ;
  cpm->link_ = new CoinBigIndex[cpm->maxlink_] ;

  cpm->nrows_ = si->getNumRows() ;
  cpm->ncols_ = si->getNumCols() ;

  cpm->sol_=sol_in;
  cpm->rowduals_=NULL;
  cpm->acts_=acts_in;

  cpm->rcosts_=NULL;
  cpm->colstat_=colstat_in;
  cpm->rowstat_=rowstat_in;

  // this is the *reduced* model, which is probably smaller
  int ncols1 = cpm->ncols_ ;
  int nrows1 = cpm->nrows_ ;

  const CoinPackedMatrix *m = si->getMatrixByCol();
  const CoinBigIndex nelemsr = m->getNumElements();

  if (isGapFree(*m)) {
    CoinDisjointCopyN(m->getVectorStarts(), ncols1, cpm->mcstrt_);
    CoinZeroN(cpm->mcstrt_+ncols1,cpm->ncols0_-ncols1);
    cpm->mcstrt_[cpm->ncols_] = nelems0;	// points to end of bulk store
    CoinDisjointCopyN(m->getVectorLengths(),ncols1,  cpm->hincol_);
    CoinDisjointCopyN(m->getIndices(),      nelemsr, cpm->hrow_);
    CoinDisjointCopyN(m->getElements(),     nelemsr, cpm->colels_);
  }
  else
  {
    CoinPackedMatrix* mm = new CoinPackedMatrix(*m);
    if( mm->hasGaps()) mm->removeGaps();
    assert(nelemsr == mm->getNumElements());
    CoinDisjointCopyN(mm->getVectorStarts(), ncols1, cpm->mcstrt_);
    CoinZeroN(cpm->mcstrt_+ncols1,cpm->ncols0_-ncols1);
    cpm->mcstrt_[cpm->ncols_] = nelems0;  // points to end of bulk store
    CoinDisjointCopyN(mm->getVectorLengths(),ncols1,  cpm->hincol_);
    CoinDisjointCopyN(mm->getIndices(),      nelemsr, cpm->hrow_);
    CoinDisjointCopyN(mm->getElements(),     nelemsr, cpm->colels_);
  }

#if PRESOLVE_DEBUG > 0 || PRESOLVE_CONSISTENCY > 0
  memset(cpm->cdone_, -1, cpm->ncols0_);
  memset(cpm->rdone_, -1, cpm->nrows0_);
#endif

  cpm->rowduals_ = new double[cpm->nrows0_];
  CoinDisjointCopyN(si->getRowPrice(), nrows1, cpm->rowduals_);
  cpm->rcosts_ = new double[cpm->ncols0_];
  CoinDisjointCopyN(si->getReducedCost(), ncols1, cpm->rcosts_);

#if PRESOLVE_DEBUG > 0
  // check accuracy of reduced costs (rcosts_ is recalculated reduced costs)
  si->getMatrixByCol()->transposeTimes(cpm->rowduals_,cpm->rcosts_) ;
  const double *obj = si->getObjCoefficients() ;
  const double *dj = si->getReducedCost() ;
  {
    int i;
    for (i=0;i<ncols1;i++) {
      double newDj = obj[i]-cpm->rcosts_[i];
      cpm->rcosts_[i]=newDj;
      assert (fabs(newDj-dj[i])<1.0e-1);
    }
  }
  // check reduced costs are 0 for basic variables
  {
    int i;
    for (i = 0; i < ncols1; i++)
      if (cpm->columnIsBasic(i)) assert (fabs(cpm->rcosts_[i])<1.0e-5);
    for (i=0;i<nrows1;i++)
      if (cpm->rowIsBasic(i)) assert (fabs(cpm->rowduals_[i])<1.0e-5);
  }
#endif
  /*
    CoinPresolve may, once, have handled both minimisation and maximisation,
    but hard-wired minimisation has crept in.
  */
  if (maxmin<0.0) {
    for (int i = 0 ; i < nrows1 ; i++)
      cpm->rowduals_[i] = -cpm->rowduals_[i] ;
    for (int j = 0 ; j < ncols1 ; j++) {
      cpm->rcosts_[j] = -cpm->rcosts_[j] ;
    }
  }

  /*
    CoinPresolve requires both column solution and row activity for correct
    operation.
  */
  CoinDisjointCopyN(si->getColSolution(), ncols1, cpm->sol_);
  CoinDisjointCopyN(si->getRowActivity(), nrows1, cpm->acts_) ;
  si->setDblParam(OsiObjOffset,cpm->originalOffset_);

  for (int j = 0; j < ncols1; j++) {
    CoinBigIndex kcs = cpm->mcstrt_[j];
    CoinBigIndex kce = kcs + cpm->hincol_[j];
    for (CoinBigIndex k=kcs; k<kce; ++k) {
      cpm->link_[k] = k+1;
    }
    if (kce>0)
      cpm->link_[kce-1] = NO_LINK ;
  }
  if (cpm->maxlink_>0) {
    CoinBigIndex ml = cpm->maxlink_;
    for (CoinBigIndex k=nelemsr; k<ml; ++k)
      cpm->link_[k] = k+1;
    cpm->link_[ml-1] = NO_LINK;
  }
  cpm->free_list_ = nelemsr;

#if PRESOLVE_DEBUG > 0 || PRESOLVE_CONSISTENCY > 0
  /*
    These are used to track the action of postsolve transforms during debugging.
  */
  CoinFillN(cpm->cdone_,ncols1,PRESENT_IN_REDUCED) ;
  CoinZeroN(cpm->cdone_+ncols1,ncols0_in-ncols1) ;
  CoinFillN(cpm->rdone_,nrows1,PRESENT_IN_REDUCED) ;
  CoinZeroN(cpm->rdone_+nrows1,nrows0_in-nrows1) ;
# endif

  return cpm;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
