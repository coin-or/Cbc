// Copyright (C) 2005-2009, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     07/21/05
//
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//---------------------------------------------------------------------------
#include "CglLandP.hpp"
#include "CglLandPSimplex.hpp"
#include "OsiRowCutDebugger.hpp"

#define INT_INFEAS(value) fabs(value - floor(value+0.5))

#include "CglConfig.h"

#ifdef CGL_HAS_OSICLP
#include "OsiClpSolverInterface.hpp"
#endif

#define CLONE_SI //Solver is cloned between two cuts

#include "CoinTime.hpp"
#include "CglGomory.hpp"
#include "CoinFactorization.hpp"
#include <fstream>
namespace LAP
{
//Setup output messages
LapMessages::LapMessages( )
        :CoinMessages(LAP_MESSAGES_DUMMY_END)
{
    strcpy(source_,"Lap");
    addMessage(BEGIN_ROUND,CoinOneMessage( 1, 2,"Starting %s round %d variable considered for separation."));
    addMessage(END_ROUND,CoinOneMessage(2, 2,"End ouf %s round %d cut generated in %g seconds."));
    addMessage(DURING_SEP,CoinOneMessage(3,1,"After %g seconds, separated %d cuts."));
    addMessage(CUT_REJECTED, CoinOneMessage(4,1,"Cut rejected for %s."));
    addMessage(CUT_FAILED,CoinOneMessage(5,1,"Generation failed."));
    addMessage(CUT_GAP, CoinOneMessage(7,1,"CUTGAP after %i pass objective is %g"));
    addMessage(LAP_CUT_FAILED_DO_MIG, CoinOneMessage(3006,1,"Failed to generate a cut generate a Gomory cut instead"));
}
}
using namespace LAP;
CglLandP::Parameters::Parameters():
        CglParam(),
        pivotLimit(20),
        pivotLimitInTree(10),
        maxCutPerRound(5000),
        failedPivotLimit(1),
        degeneratePivotLimit(0),
        extraCutsLimit(5),
	maximumCandidates(1000000),
	maximumCutLength(10000),
        pivotTol(1e-4),
        away(5e-4),
        timeLimit(COIN_DBL_MAX),
        singleCutTimeLimit(COIN_DBL_MAX),
        rhsWeight(1.),
        useTableauRow(true),
        modularize(false),
        strengthen(true),
        countMistakenRc(false),
        sepSpace(Fractional),
        perturb(true),
        normalization(Unweighted),
        rhsWeightType(Fixed),
        lhs_norm(L1),
        generateExtraCuts(none),
        pivotSelection(mostNegativeRc)
{
    EPS = 1e-08;
}

CglLandP::Parameters::Parameters(const Parameters &other):
        CglParam(other),
        pivotLimit(other.pivotLimit),
        pivotLimitInTree(other.pivotLimitInTree),
        maxCutPerRound(other.maxCutPerRound),
        failedPivotLimit(other.failedPivotLimit),
        degeneratePivotLimit(other.degeneratePivotLimit),
        extraCutsLimit(other.extraCutsLimit),
	maximumCandidates(other.maximumCandidates),
	maximumCutLength(other.maximumCutLength),
        pivotTol(other.pivotTol),
        away(other.away),
        timeLimit(other.timeLimit),
        singleCutTimeLimit(other.singleCutTimeLimit),
        rhsWeight(other.rhsWeight),
        useTableauRow(other.useTableauRow),
        modularize(other.modularize),
        strengthen(other.strengthen),
        countMistakenRc(other.countMistakenRc),
        sepSpace(other.sepSpace),
        perturb(other.perturb),
        normalization(other.normalization),
        rhsWeightType(other.rhsWeightType),
        lhs_norm(other.lhs_norm),
        generateExtraCuts(other.generateExtraCuts),
        pivotSelection(other.pivotSelection)
{}

CglLandP::Parameters & CglLandP::Parameters::operator=(const Parameters &other)
{
    if (this != &other)
    {
        CglParam::operator=(other);
        pivotLimit = other.pivotLimit;
        pivotLimitInTree = other.pivotLimitInTree;
        maxCutPerRound = other.maxCutPerRound;
        failedPivotLimit = other.failedPivotLimit;
        degeneratePivotLimit = other.failedPivotLimit;
        extraCutsLimit = other.extraCutsLimit;
	maximumCandidates = other.maximumCandidates;
	maximumCutLength = other.maximumCutLength;
        pivotTol = other.pivotTol;
        away = other.away;
        timeLimit = other.timeLimit;
        singleCutTimeLimit = other.singleCutTimeLimit;
        rhsWeight = other.rhsWeight;
        useTableauRow = other.useTableauRow;
        modularize = other.modularize;
        strengthen = other.strengthen;
        countMistakenRc = other.countMistakenRc;
        sepSpace = other.sepSpace;
        perturb = other.perturb;
        normalization = other.normalization;
        rhsWeightType = other.rhsWeightType;
        lhs_norm = other.lhs_norm;
        generateExtraCuts = other.generateExtraCuts;
        pivotSelection = other.pivotSelection;
    }
    return *this;
}

CglLandP::CachedData::CachedData(int nBasics, int nNonBasics):
        basics_(NULL), nonBasics_(NULL), nBasics_(nBasics),
        nNonBasics_(nNonBasics), basis_(NULL), colsol_(NULL),
        slacks_(NULL), integers_(NULL), solver_(NULL)
{
    if (nBasics_>0)
    {
        basics_ = new int[nBasics_];
        integers_ = new bool [nNonBasics_ + nBasics_];
    }
    if (nNonBasics_>0)
        nonBasics_ = new int[nNonBasics_];
    if (nBasics_ + nNonBasics_ > 0)
    {
        colsol_ = new double[nBasics_ + nNonBasics_];
        slacks_ = &colsol_[nNonBasics_];
    }
}

CglLandP::CachedData::CachedData(const CachedData &source):
  basics_(NULL), nonBasics_(NULL), nBasics_(source.nBasics_),
        nNonBasics_(source.nNonBasics_), basis_(NULL),
        colsol_(NULL), slacks_(NULL), integers_(NULL), solver_(NULL)
{
    if (nBasics_>0)
    {
        basics_ = new int[nBasics_];
        CoinCopyN(source.basics_, nBasics_, basics_);
        integers_ = new bool [nNonBasics_ + nBasics_];
        CoinCopyN(source.integers_, nBasics_ + nNonBasics_, integers_);
    }
    if (nNonBasics_>0)
    {
        nonBasics_ = new int[nNonBasics_];
        CoinCopyN(source.nonBasics_, nBasics_, nonBasics_);
    }
    if (nBasics_ + nNonBasics_ > 0)
    {
        colsol_ = new double[nBasics_ + nNonBasics_];
        slacks_ = &colsol_[nNonBasics_];
        CoinCopyN(source.colsol_, nBasics_ + nNonBasics_, colsol_);
    }
    if (source.basis_!=NULL)
        basis_ = new CoinWarmStartBasis(*source.basis_);
    if (source.solver_!=NULL)
      solver_ = source.solver_->clone();
}

CglLandP::CachedData& CglLandP::CachedData::operator=(const CachedData &source)
{
    if (this != &source)
    {
        nBasics_ = source.nBasics_;
        nNonBasics_ = source.nNonBasics_;
        delete [] basics_;
        basics_ = NULL;
        delete [] nonBasics_;
        nonBasics_ = NULL;
        delete [] basis_;
        basis_ = NULL;
        delete [] colsol_;
        colsol_ = NULL;
        delete [] slacks_;
        slacks_ = NULL;
        delete [] integers_;
        integers_ = NULL;
        if (nBasics_>0)
        {
            basics_ = new int[nBasics_];
            CoinCopyN(source.basics_, nBasics_, basics_);
            integers_ = new bool [nBasics_ + nNonBasics_];
            CoinCopyN(source.integers_, nBasics_ + nNonBasics_, integers_);
        }
        if (nNonBasics_>0)
        {
            nonBasics_ = new int[nNonBasics_];
            CoinCopyN(source.nonBasics_, nBasics_, nonBasics_);
        }
        if (nBasics_ + nNonBasics_ > 0)
        {
            colsol_ = new double[nBasics_ + nNonBasics_];
            slacks_ = &colsol_[nNonBasics_];
            CoinCopyN(source.colsol_, nBasics_ + nNonBasics_, colsol_);
        }
        if (source.basis_!=NULL)
            basis_ = new CoinWarmStartBasis(*source.basis_);
        delete solver_;
	if (source.solver_)
	  solver_ = source.solver_->clone();
    }
    return *this;
}

void
CglLandP::CachedData::getData(const OsiSolverInterface &si)
{
    int nBasics = si.getNumRows();
    int nNonBasics = si.getNumCols();
    if (basis_ != NULL)
        delete basis_;
    basis_ = dynamic_cast<CoinWarmStartBasis *> (si.getWarmStart());
    if (!basis_)
        throw NoBasisError();

    if (nBasics_ > 0 || nBasics != nBasics_)
    {
        delete [] basics_;
        basics_ = NULL;
    }
    if (basics_ == NULL)
    {
        basics_ = new int[nBasics];
        nBasics_ = nBasics;
    }

    if (nNonBasics_ > 0 || nNonBasics != nNonBasics_)
    {
        delete [] nonBasics_;
        nonBasics_ = NULL;
    }
    if (nonBasics_ == NULL)
    {
        nonBasics_ = new int[nNonBasics];
        nNonBasics_ = nNonBasics;
    }
    int n = nBasics + nNonBasics;
    if ( nBasics_ + nNonBasics_ > 0 || nBasics_ + nNonBasics_ != n)
    {
        delete [] colsol_;
        delete [] integers_;
        integers_ = NULL;
        colsol_ = NULL;
        slacks_ = NULL;
    }
    if (colsol_ == NULL)
    {
        colsol_ = new double[n];
        slacks_ = &colsol_[nNonBasics];
    }

    if (integers_ == NULL)
    {
        integers_ = new bool[n];
    }

    
    const double * rowLower = si.getRowLower();
    const double * rowUpper = si.getRowUpper();
    //determine which slacks are integer
    const CoinPackedMatrix * m = si.getMatrixByCol();
    const double * elems = m->getElements();
    const int * inds = m->getIndices();
    const CoinBigIndex * starts = m->getVectorStarts();
    const int * lengths = m->getVectorLengths();
    //    int numElems = m->getNumElements();
    int numCols = m->getNumCols();
    assert(numCols == nNonBasics_);
    //   int numRows = m->getNumRows();
    CoinFillN(integers_ ,n, true);
    for (int i = 0 ;  i < numCols ; i++)
    {
        if (si.isContinuous(i))
            integers_[i] = false;
    }
    bool * integerSlacks = integers_ + numCols;
    for (int i = 0 ; i < nBasics ; i++)
    {
        if (rowLower[i] > -1e50 && INT_INFEAS(rowLower[i]) > 1e-15)
            integerSlacks[i] = false;
        if (rowUpper[i] < 1e50 && INT_INFEAS(rowUpper[i]) > 1e-15)
            integerSlacks[i] = false;
    }
    for (int i = 0 ;  i < numCols ; i++)
    {
        CoinBigIndex end = starts[i] + lengths[i];
        if (integers_[i])
        {
            for (CoinBigIndex k=starts[i] ; k < end; k++)
            {
                if (integerSlacks[inds[k]] && INT_INFEAS(elems[k])>1e-15 )
                    integerSlacks[inds[k]] = false;
            }
        }
        else
        {
            for (CoinBigIndex k=starts[i] ; k < end; k++)
            {
                if (integerSlacks[inds[k]])
                    integerSlacks[inds[k]] = false;
            }
        }
    }

    CoinCopyN(si.getColSolution(), si.getNumCols(), colsol_);
    CoinCopyN(si.getRowActivity(), si.getNumRows(), slacks_);
    for (int i = 0 ; i < si.getNumRows() ; i++)
    {
        slacks_[i]*=-1;
        if (rowLower[i]>-1e50)
        {
            slacks_[i] += rowLower[i];
        }
        else
        {
            slacks_[i] += rowUpper[i];
        }
    }
    //Now fill the arrays;
    nNonBasics = 0;
    nBasics = 0;



    //For having the index variables correctly ordered we need to access to OsiSimplexInterface
    {
        OsiSolverInterface * ncSi = (const_cast<OsiSolverInterface *>(&si));
        ncSi->enableSimplexInterface(0);
        ncSi->getBasics(basics_);
	// Save enabled solver
	solver_ = si.clone();
#ifdef CGL_HAS_OSICLP
	OsiClpSolverInterface * clpSi = getClpSolver(solver_);
	const OsiClpSolverInterface * clpSiRhs = getConstClpSolver(&si);
	if (CBC_SKIP_CLP_TEST||clpSi)
	  clpSi->getModelPtr()->copyEnabledStuff(clpSiRhs->getModelPtr());;
#endif
        ncSi->disableSimplexInterface();
    }

    int numStructural = basis_->getNumStructural();
    for (int i = 0 ; i < numStructural ; i++)
    {
        if (basis_->getStructStatus(i)== CoinWarmStartBasis::basic)
        {
	    nBasics++;
            //Basically do nothing
        }
        else
        {
            nonBasics_[nNonBasics++] = i;
        }
    }

    int numArtificial = basis_->getNumArtificial();
    int numStruct = basis_->getNumStructural();
    for (int i = 0 ; i < numArtificial ; i++)
    {
        if (basis_->getArtifStatus(i)== CoinWarmStartBasis::basic)
        {
            //Just check number of basics
            nBasics++;
        }
        else
        {
            nonBasics_[nNonBasics++] = i + numStruct;
        }
    }
}
void
CglLandP::CachedData::clean(){
    if (basics_!=NULL)
        delete [] basics_;
    basics_ = NULL;
    if (nonBasics_!=NULL)
        delete [] nonBasics_;
    nonBasics_ = NULL;
    if (colsol_ != NULL)
        delete [] colsol_;
    colsol_ = NULL;
    delete basis_;
    basis_ = NULL;
    if (integers_)
        delete [] integers_;
    integers_ = NULL;

   nBasics_ = 0;
   nNonBasics_ = 0;
   delete solver_;
   solver_ = NULL;
}
CglLandP::CachedData::~CachedData()
{
    if (basics_!=NULL)
        delete [] basics_;
    if (nonBasics_!=NULL)
        delete [] nonBasics_;
    if (colsol_ != NULL)
        delete [] colsol_;
    delete basis_;
    if (integers_)
        delete [] integers_;
    delete solver_;
}

CglLandP::CglLandP(const CglLandP::Parameters &params,
                   const LAP::Validator &validator):
        params_(params), cached_(), validator_(validator), numrows_(-1),
        numcols_(-1),originalColLower_(NULL), originalColUpper_(NULL),
        canLift_(false),
        extraCuts_()
{
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(0);
    messages_ = LapMessages();
}


CglLandP::~CglLandP()
{
    delete handler_;
    if (originalColLower_ != NULL)
        delete [] originalColLower_;
    if (originalColUpper_ != NULL)
        delete [] originalColUpper_;
}

CglLandP::CglLandP(const CglLandP & source):
        CglCutGenerator(source),
        params_(source.params_), cached_(source.cached_),
        validator_(source.validator_), numrows_(source.numrows_),numcols_(source.numcols_),
        originalColLower_(NULL), originalColUpper_(NULL),
        canLift_(source.canLift_),
        extraCuts_(source.extraCuts_)
{
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(source.handler_->logLevel());
    messages_ = LapMessages();
    if (numcols_ != -1)
    {
        assert(numcols_ > 0);
        assert(originalColLower_!=NULL);
        assert(originalColUpper_!=NULL);
        originalColLower_ = new double[numcols_];
        originalColUpper_ = new double[numcols_];
        CoinCopyN(source.originalColLower_,numcols_,originalColLower_);
        CoinCopyN(source.originalColUpper_,numcols_,originalColUpper_);
    }
}

/** Assignment operator */
CglLandP& CglLandP::operator=(const CglLandP &rhs)
{
    if (this != &rhs)
    {
        params_ = rhs.params_;
        cached_ = rhs.cached_;
        validator_ = rhs.validator_;
        extraCuts_ = rhs.extraCuts_;
    }
    return *this;
}


CglCutGenerator *
CglLandP::clone() const
{
    return new CglLandP(*this);
}

extern double restaurationTime;

struct cutsCos
{
    int i;
    int j;
    double angle;
    cutsCos(int  i_, int j_ , double angle_):i(i_), j(j_), angle(angle_)
    {
    }
    bool operator<(const cutsCos&other)const
    {
        return angle > other.angle;
    }
};


void
CglLandP::scanExtraCuts(OsiCuts& cs, const double * colsol) const
{
    //int numAdded = 0;
    for (int i = extraCuts_.sizeRowCuts() - 1; i > -1 ; i--)
    {
        double violation = extraCuts_.rowCut(i).violated(colsol);
        if (violation > 0.)
        {
            cs.insert(extraCuts_.rowCut(i));
            //numAdded++;
            //      std::cout<<"A cut computed in a previous iteration is violated by "<<violation<<"."<<std::endl;
            //extraCuts_.eraseRowCut(i);
        }
    }
    //  std::cout<<"Added "<<numAdded<<" previously generated cuts."<<std::endl;
}

void
CglLandP::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
                       const CglTreeInfo info )
{
    int numberRanges = 0;
    if (numrows_<0)
    {
      numrows_ = si.getNumRows();
      // but switch off? if ranges
      const double * rowLower = si.getRowLower();
      const double * rowUpper = si.getRowUpper();
      for (int i=0;i<numrows_;i++) {
	if (rowLower[i]<rowUpper[i]) {
	  if (rowLower[i]> -1.0e50 && rowUpper[i] < 1.0e50) {
	    numberRanges++;
	  }
	}
      }
      if (numberRanges && false) {
	params_.maximumCutLength =-1;
	return;
      }
    } else if (params_.maximumCutLength < 0) {
      return;
    }
// scanExtraCuts(cs, si.getColSolution());
    Parameters params = params_;
    params.rhsWeight = numrows_ + 2;

    handler_->message(CUT_GAP, messages_)<<info.pass<<si.getObjValue() <<CoinMessageEol;

    if (info.inTree)   //put lower pivot limit
    {
        params.pivotLimit = std::min(params.pivotLimit, params.pivotLimitInTree);
        params.countMistakenRc = true;
    }
    if (params.timeLimit < 0)
    {
        params.pivotLimit = 0;
    }

    assert(si.basisIsAvailable());


#ifdef APPEND_ROW
    OsiSolverInterface * t_si = si.clone();
    if (params.modularize)
    {
      if (numberRanges) {
	// modify
      }
      int new_idx = si.getNumCols();
      int v_idx[1] = {new_idx};
      double v_val[1] = {-1};
      CoinPackedVector v(1, v_idx, v_val, false);
      t_si->addCol(CoinPackedVector(), 0, 1, 0);
      t_si->setInteger(new_idx);
      t_si->addRow(v,0, 0);
      t_si->resolve();
    }
#else
    const OsiSolverInterface * t_si = &si;
    int numberCuts = cs.sizeRowCuts(); 
    double * upper = NULL;
    CoinBigIndex * starts = NULL;
    int * row = NULL;
    int nThrownAway=0;
    if (numberRanges) {
      // modify by adding slacks
      OsiSolverInterface *tt_si = si.clone();
      const double * rowLower = si.getRowLower();
      const double * rowUpper = si.getRowUpper();
      const double * rowActivity = si.getRowActivity();
      upper = new double[4*numberRanges];
      double * lower = upper+numberRanges;
      double * obj = lower+numberRanges;
      double * element = obj+numberRanges;
      memset(lower,0,2*numberRanges*sizeof(double));
      starts = new CoinBigIndex[numberRanges+1];
      row = new int[numberRanges];
      // First add columns
      numberRanges = 0;
      starts[0] = 0;
      for (int i=0;i<numrows_;i++) {
	if (rowLower[i]<rowUpper[i]) {
	  if (rowLower[i]> -1.0e50 && rowUpper[i] < 1.0e50) {
	    double value = rowActivity[i];
	    if (value-rowLower[i] < rowUpper[i]-value) {
	      // keep rowLower
	      //tt_si->setRowUpper(i,rowLower[i]);
	      element[numberRanges] = -1.0;
	    } else {
	      // keep rowUpper
	      //tt_si->setRowLower(i,rowUpper[i]);
	      element[numberRanges] = 1.0;
	    }
	    row[numberRanges] = i;
	    upper[numberRanges++] = rowUpper[i]-rowLower[i];
	    starts[numberRanges] = numberRanges;
	  }
	}
      }
      tt_si->addCols(numberRanges,starts,row,element,lower,upper,obj);
      // Now basis and values
      int nCol = si.getNumCols();
      double * solution = CoinCopyOfArray(tt_si->getColSolution(),
					  nCol+numberRanges);
      CoinWarmStartBasis * basis
	= dynamic_cast<CoinWarmStartBasis *>(tt_si->getWarmStart());
      numberRanges = 0;
      for (int i=0;i<numrows_;i++) {
	if (rowLower[i]<rowUpper[i]) {
	  if (rowLower[i]> -1.0e50 && rowUpper[i] < 1.0e50) {
	    double value = rowActivity[i];
	    if (value-rowLower[i] < rowUpper[i]-value) {
	      // keep rowLower
	      tt_si->setRowUpper(i,rowLower[i]);
	      value -= rowLower[i];
	    } else {
	      // keep rowUpper
	      tt_si->setRowLower(i,rowUpper[i]);
	      value = rowUpper[i]-value;
	    }
	    solution[nCol+numberRanges] = value;
	    if (basis->getArtifStatus(i) ==
		CoinWarmStartBasis::basic) {
	      // set basic
	      basis->setStructStatus(nCol+numberRanges,
				     CoinWarmStartBasis::basic);
	      basis->setArtifStatus(i,
				    CoinWarmStartBasis::atLowerBound);
	    }
	    numberRanges++;
	  }
	}
      }
      // update solution and basis
      tt_si->setColSolution(solution);
      tt_si->setWarmStart(basis);
      delete basis;
      delete [] solution;
      tt_si->resolve();
      t_si = tt_si;
    }
#endif

    cached_.getData(*t_si);
    CglLandPSimplex landpSi(*t_si, cached_, params, validator_);
    if (params.generateExtraCuts == CglLandP::AllViolatedMigs)
    {
        landpSi.genThisBasisMigs(cached_, params);
    }
    landpSi.setLogLevel(handler_->logLevel());
    int nCut = 0;

    std::vector<int> indices;
    getSortedFractionalIndices(indices,cached_, params);
    if (indices.size()>params.maximumCandidates)
      indices.resize(params.maximumCandidates);

#ifndef NDEBUG
    int numrows = si.getNumRows();
#endif

#ifdef DO_STAT
    //Get informations on current optimum
    {
        OsiSolverInterface * gapTester = si.clone();
        gapTester->resolve();

        roundsStats_.analyseOptimalBasis(gapTester,info.pass, numrows_);
        delete gapTester;
    }
#endif

    params_.timeLimit += CoinCpuTime();
    CoinRelFltEq eq(1e-04);

    for (unsigned int i = 0; i < indices.size() && nCut < params.maxCutPerRound &&
            nCut < cached_.nBasics_ ; i++)
    {

        //Check for time limit
        int iRow = indices[i];
        assert(iRow < numrows);
        OsiRowCut cut;
        int code=1;
        OsiSolverInterface * ncSi = NULL;

        if (params.pivotLimit != 0)
        {
            ncSi = t_si->clone();
            landpSi.setSi(ncSi);
            ncSi->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
            ncSi->messageHandler()->setLogLevel(0);
        }

        int generated = 0;
        if (params.pivotLimit == 0)
        {
            generated = landpSi.generateMig(iRow, cut, params);
        }
        else
        {
            generated = landpSi.optimize(iRow, cut, cached_, params);
            if (params.generateExtraCuts == CglLandP::AllViolatedMigs)
            {
                landpSi.genThisBasisMigs(cached_, params);
            }
            landpSi.resetSolver(cached_.basis_);
        }
        code = 0;
        if (generated)
            code = validator_(cut, cached_.colsol_, si, params, originalColLower_, originalColUpper_);
        if (!generated || code)
        {
            if (params.pivotLimit !=0)
            {
                handler_->message(LAP_CUT_FAILED_DO_MIG, messages_)<<validator_.failureString(code)<<CoinMessageEol;
                landpSi.freeSi();
                OsiSolverInterface * ncSi = t_si->clone();
                landpSi.setSi(ncSi);
                params.pivotLimit = 0;
                if (landpSi.optimize(iRow, cut, cached_, params))
                {
                    code = validator_(cut, cached_.colsol_, si, params, originalColLower_, originalColUpper_);
                }
                params.pivotLimit = params_.pivotLimit;
            }
        }

        if (params.pivotLimit != 0)
        {
            landpSi.freeSi();
        }
        if (code)
        {
            handler_->message(CUT_REJECTED, messages_)<<
            validator_.failureString(code)<<CoinMessageEol;
        }
        else
        {
  	    if (numberRanges) {
	      int numberColumns = si.getNumCols();
	      int n = cut.row().getNumElements();
	      const int *column = cut.row().getIndices();
	      const double *element = cut.row().getElements();
	      bool oddSlack = false;
	      for (int i=0;i<n;i++) {
		if (column[i]>=numberColumns) {
		  oddSlack=true;
		}
	      }
	      if (oddSlack) {
		nThrownAway++;
#if 0
		// for now throw away
		printf("odd cut\n %d entries ncol=%d, nrange %d %g<=%g\n",
		       n,numberColumns,numberRanges,cut.lb(),cut.ub());
		for (int i=0;i<n;i++) {
		  printf("(%g*x%d) ",column[i],element[i]);
		  if ((i%5)==4)
		    printf("\n");
		}
		if ((n%5))
		  printf("\n");
#endif
		cut = OsiRowCut();
	      }
	    }
            if (canLift_)
            {
                cut.setGloballyValid(true);
            }
#ifdef CHECK_KNOWN_SOLUTION
	    const OsiRowCutDebugger *debugger = si.getRowCutDebugger();
	    if (debugger) {
	      if (debugger->invalidCut(cut)) {
		printf("BAD cut\n");
		exit(0);
	      }
	    }
	    //CoinAssert (!debugger->invalidCut(*cut));
#endif
            cs.insertIfNotDuplicate(cut, eq);
            //cs.insertIfNotDuplicate(cut);
            {
                //std::cout<<"Violation "<<cut.violated(cached_.colsol_)<<std::endl;
                nCut++;
            }
        }
    }

    Cuts& extra = landpSi.extraCuts();
    for (int i = 0 ; i < cached_.nNonBasics_; i++)
    {
        OsiRowCut * cut = extra.rowCut(i);
        if (cut == NULL) continue;
        int code = validator_(*cut, cached_.colsol_, si, params,
                              originalColLower_, originalColUpper_);
        if (code)
        {
            handler_->message(LAP_CUT_FAILED_DO_MIG, messages_)
            <<validator_.failureString(code)<<CoinMessageEol;
        }
        else
        {
            cs.insertIfNotDuplicate(*cut, eq);
            {
                nCut++;
            }
        }
        delete cut;
    }

    landpSi.outPivInfo(nCut);
    params_.timeLimit -= CoinCpuTime();

    cached_.clean();
#ifdef APPEND_ROW
    assert(t_si != &si);
    delete t_si;
#else
    if (t_si != &si) {
      delete [] upper;
      delete [] starts;
      delete [] row;
      delete t_si;
      assert (numberRanges);
      //printf("Ranges (%d) - switching off CglLandP after this %d cuts, %d thrown\n",
      //     numberRanges,
      //     cs.sizeRowCuts()-numberCuts, 
      //     nThrownAway);
      params_.maximumCutLength = -params_.maximumCutLength;
    }
#endif
    return;
}


template < class S, class T, class U >
class StableCompare
{
public:
    inline bool operator()(const CoinTriple<S,T,U>& t1,
                           const CoinTriple<S,T,U>& t2) const
    {
        return (t1.third < t2.third) ||
               ((t1.third == t2.third) && (t1.second < t2.second));
    }

};

template <class T1,class T2>
struct StableExternalComp
{
    const std::vector<T1> &vec_1_;
    const std::vector<T2> &vec_2_;
    StableExternalComp(const std::vector<T1> &vec_1,
                       const std::vector<T2> &vec_2):
            vec_1_(vec_1),
            vec_2_(vec_2)
    {
    }
    CoinRelFltEq eq;
    bool operator()(int i, int j)
    {
        bool result = (vec_1_[i] < vec_1_[j]) ||
                      ( ((vec_1_[i]== vec_1_[j]))
                        && (vec_2_[i] < vec_2_[j]));
        return result;
    }

};
void
CglLandP::getSortedFractionalIndices(std::vector<int> &frac_indices,
                                     const CachedData &data,
                                     const CglLandP::Parameters & params) const
{
    std::vector<int> colIndices;
    std::vector<double> values;
    std::vector<int> indices;
    for (int i = 0 ; i < data.nBasics_ ; i++)
    {
        const int& iCol = data.basics_[i];
        if (iCol >= data.nNonBasics_ ||
                !data.integers_[iCol] ||
                INT_INFEAS(data.colsol_[iCol]) <= params.away)
            continue;
        const double value = INT_INFEAS(data.colsol_[iCol]);

        frac_indices.push_back(i);
        indices.push_back(static_cast<int>(values.size()));
        values.push_back(- value);
        colIndices.push_back(iCol);
    }
    std::sort(indices.begin(), indices.end(),StableExternalComp<double, int>(values,colIndices));
    colIndices = frac_indices;
    for (unsigned int i = 0; i < indices.size() ; i++)
    {
        frac_indices[i] = colIndices[indices[i]];
    }

}


