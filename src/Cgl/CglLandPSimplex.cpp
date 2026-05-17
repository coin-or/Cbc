// Copyright (C) 2005-2009, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     21/07/05
//
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//---------------------------------------------------------------------------
#include "CglLandPSimplex.hpp"
#include "CoinTime.hpp"
#ifdef CGL_HAS_OSICLP
#include "OsiClpSolverInterface.hpp"
#endif
#include "CoinIndexedVector.hpp"
#include <cassert>
#include <iterator>

#include <list>
#include <algorithm>

#define REMOVE_LOG 0

#define RED_COST_CHECK 1e-6


//#define OUT_CGLP_PIVOTS
#ifdef OUT_CGLP_PIVOTS
#include "CglLandPOutput.hpp"
#endif
#ifdef DEBUG_LAP

/* The function is not used anywhere (LL)
static void MyAssertFunc(bool c, const std::string &s, const std::string&  file, unsigned int line){
    if (c != true){
        fprintf(stderr, "Failed MyAssertion: %s in %s line %i.\n", s.c_str(), file.c_str(), line);
        throw -1;
    }
}
*/

static void DblGtAssertFunc(const double& a, const std::string &a_s, const double&b, const std::string& b_s,
                     const std::string&  file, unsigned int line)
{
    if (a<b)
    {
        fprintf(stderr, "Failed comparison: %s = %f < %s =%f in  %s line %i.\n", a_s.c_str(),
                a, b_s.c_str(), b, file.c_str(), line);
        throw -1;
    }
}

static void DblEqAssertFunc(const double& a, const std::string &a_s, const double&b, const std::string& b_s,
                     const std::string&  file, unsigned int line)
{
    CoinRelFltEq eq(1e-7);
    if (!eq(a,b))
    {
        fprintf(stderr, "Failed comparison: %s = %f != %s =%f in  %s line %i.\n", a_s.c_str(),
                a, b_s.c_str(), b, file.c_str(), line);
        throw -1;
    }
}

static void VecModEqAssertFunc(const CoinIndexedVector& a, const std::string a_s,
                        const CoinIndexedVector& b, const std::string b_s,
                        const std::string file, unsigned int line)
{
    CoinRelFltEq eq(1e-7);
    assert(a.capacity()==b.capacity());
    int n = a.capacity();
    const double * a_v = a.denseVector();
    const double * b_v = b.denseVector();
    bool failed=false;
    for (int i = 0 ; i < n ; i++)
    {
        double a = a_v[i] - floor(a_v[i]);
        double b = b_v[i] - floor(b_v[i]);
        if (!eq(a,b) && !eq(a_v[i],b_v[i]))
        {
            fprintf(stderr, "Failed comparison: %s[%i] = %.15f != %s[%i] =%.15f in  %s line %i.\n", a_s.c_str(),i,
                    a_v[i], b_s.c_str(), i, b_v[i], file.c_str(), line);
            failed = true;
        }
    }
    if (failed) throw -1;
}

static void VecEqAssertFunc(const CoinIndexedVector& a, const std::string a_s,
                     const CoinIndexedVector& b, const std::string b_s,
                     const std::string file, unsigned int line)
{
    CoinRelFltEq eq(1e-7);
    assert(a.capacity()==b.capacity());
    int n = a.capacity();
    const double * a_v = a.denseVector();
    const double * b_v = b.denseVector();
    bool failed=false;
    for (int i = 0 ; i < n ; i++)
    {
        if (!eq(a_v[i],b_v[i]))
        {
            fprintf(stderr, "Failed comparison: %s[%i] = %f != %s[%i] =%f in  %s line %i.\n", a_s.c_str(),i,
                    a_v[i], b_s.c_str(), i, b_v[i], file.c_str(), line);
            failed = true;
        }
    }
    if (failed) throw -1;
}


#define MAKE_STRING(exp) std::string(#exp)
#define MyAssert(exp)  MyAssertFunc(exp, MAKE_STRING(exp), __FILE__, __LINE__);
#define DblEqAssert(a,b)  DblEqAssertFunc(a,MAKE_STRING(a),b,MAKE_STRING(b), __FILE__, __LINE__);
#define DblGtAssert(a,b)  DblGtAssertFunc(a,MAKE_STRING(a),b,MAKE_STRING(b), __FILE__, __LINE__);

#define VecEqAssert(a,b) VecEqAssertFunc(a, MAKE_STRING(a), b, MAKE_STRING(b), __FILE__, __LINE__);
#define VecModEqAssert(a,b) VecModEqAssertFunc(a, MAKE_STRING(a), b, MAKE_STRING(b), __FILE__, __LINE__);

void checkVecFunc(const CoinIndexedVector &v)
{
    CoinIndexedVector v2 = v;
    double * x = v2.denseVector();
    const int * idx = v2.getIndices();
    int n = v2.getNumElements();
    for (int i = 0 ; i < n ; i++)
    {
        x[idx[i]] = 0.;
    }
    n = v2.capacity();
    for (int i = 0 ; i < n ; i++)
    {
        DblEqAssert(x[i],0.);
    }
}

#define checkVec(a) checkVecFunc(a)
#else
#define MyAssert(exp)
#define DblEqAssert(a,b)
#define DblGtAssert(a,b)
#define VecEqAssert(a,b)
#define VecModEqAssert(a,b)

#define checkVec(a)
#endif

#include <algorithm>
//#define TEST_M3
namespace LAP
{

void CglLandPSimplex::printTableau(std::ostream & os)
{
    int width = 9;
    os<<"Tableau at current basis"<<std::endl;
    os<<"    ";
    //Head with non basics indices
    for (int i = 0 ; i < ncols_orig_ ; i++)
    {
        os.width(width);
        os.setf(std::ios_base::right, std::ios_base::adjustfield);
        std::cout<<nonBasics_[i]<<" ";
    }

    os.width(width);
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    std::cout<<'b';

    os<<std::endl;

    //print row by row
    for (int i = 0 ; i < nrows_ ; i++)
    {
        //int ind = basics_[i];
        row_i_.num = i;
        pullTableauRow(row_i_);
        row_i_.print(os, width, nonBasics_, ncols_orig_);

    }

}

bool
CglLandPSimplex::checkBasis()
{
    //Check that basics_ is correct
    int * basic2 = new int [nrows_];
    si_->getBasics(basic2);
    for (int i = 0; i < nrows_ ; i++)
        assert(basics_[i]==basic2[i]);
    delete [] basic2;
    return true;
}


CglLandPSimplex::CglLandPSimplex(const OsiSolverInterface &si,
                                 const CglLandP::CachedData &cached,
                                 const CglLandP::Parameters &params,
                                 Validator& validator):
#ifdef CGL_HAS_OSICLP
        clp_(NULL),
#endif
        row_k_(this),
        original_row_k_(this),
        row_i_(this),
#ifndef NDEBUG
        new_row_(this),
#endif
        gammas_(false),
        rowFlags_(NULL),
        col_in_subspace(),
        colCandidateToLeave_(NULL),
        basics_(NULL), nonBasics_(NULL),
        M1_(), M2_(), M3_(),
        sigma_(0), basis_(NULL), colsolToCut_(NULL),
        colsol_(NULL),
        ncols_orig_(0),nrows_orig_(0),
        inDegenerateSequence_(false),
        chosenReducedCostVal_(1e100),
        original_index_(),
        si_(NULL),
        validator_(validator),
        numPivots_(0),
        numSourceRowEntered_(0),
        numIncreased_(0)
{
    ncols_orig_ = si.getNumCols();
    nrows_orig_ = si.getNumRows();
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(2);
    messages_ = LandPMessages();

    si_ = const_cast<OsiSolverInterface *>(&si);

#ifdef CGL_HAS_OSICLP
    OsiClpSolverInterface * clpSi = getClpSolver(si_);
    if (CBC_SKIP_CLP_TEST||clpSi)
    {
        clp_ = clpSi;
    }
#endif

    int rowsize = ncols_orig_ + nrows_orig_ + 1;
    row_k_.reserve(rowsize);
#ifndef NDEBUG
    new_row_.reserve(rowsize);
#endif

    lo_bounds_.resize(ncols_orig_ + nrows_orig_);
    up_bounds_.resize(ncols_orig_ + nrows_orig_);

    CoinCopyN(si.getColLower(),ncols_orig_, &lo_bounds_[0]);
    CoinCopyN(si.getColUpper(),ncols_orig_,&up_bounds_[0]);
    const double * rowUpper = si.getRowUpper();
    const double * rowLower = si.getRowLower();
    double infty = si.getInfinity();
    int i=ncols_orig_;
    for (int iRow = 0; iRow < nrows_orig_ ; iRow++, i++)
    {
        if (rowUpper[iRow] < infty)
            lo_bounds_[i]=0.;
        else lo_bounds_[i]= - infty;
        if (rowLower[iRow] <= - infty)
            up_bounds_[i] = infty;
        else if (rowUpper[iRow] < infty)
        {
            lo_bounds_[i] = rowLower[iRow] - rowUpper[iRow];
            up_bounds_[i] = 0;
        }
        else
            up_bounds_[i] = 0.;
    }
    cuts_.resize(ncols_orig_);
    if (params.pivotLimit != 0)
    {
        own_ = true;
        rWk1_.resize(nrows_orig_);
        rWk2_.resize(nrows_orig_);
        rWk3_.resize(nrows_orig_);
        rWk4_.resize(nrows_orig_);
        rIntWork_.resize(nrows_orig_);

        row_i_.reserve(rowsize);
        rowFlags_ = new bool[nrows_orig_];
        col_in_subspace.resize(ncols_orig_ + nrows_orig_);
        colCandidateToLeave_ = new bool[ncols_orig_];
        basics_ = new int[nrows_orig_];
        nonBasics_ = new int[ncols_orig_];

        colsolToCut_ = new double[ncols_orig_ + nrows_orig_];
        colsol_ = new double[ncols_orig_ + nrows_orig_];
        original_index_.resize(ncols_orig_ + nrows_orig_);
        CoinIotaN(&original_index_[0],ncols_orig_ + nrows_orig_, 0);



    }
    else
    {
        nrows_ = nrows_orig_;
        ncols_ = ncols_orig_;
        original_index_.resize(ncols_orig_ + nrows_orig_);
        CoinIotaN(&original_index_[0],ncols_orig_ + nrows_orig_, 0);
        own_ = false;
        si_->enableSimplexInterface(0);
        basis_ = new CoinWarmStartBasis(*cached.basis_);
    }
    cacheUpdate(cached,params.sepSpace != CglLandP::Full);
    if (params.normalization)
    {
        computeWeights(params.lhs_norm, params.normalization, params.rhsWeightType);
    }
    else rhs_weight_ = 1;
}

CglLandPSimplex::~CglLandPSimplex()
{
    delete handler_;
    handler_ = NULL;
    delete basis_;
    basis_ = NULL;
    if (own_)
    {
        delete [] rowFlags_;
        rowFlags_ = NULL;
        delete [] colCandidateToLeave_;
        colCandidateToLeave_ = NULL;
        delete [] basics_;
        basics_ = NULL;
        delete [] nonBasics_;
        nonBasics_ = NULL;
        delete [] colsolToCut_;
        colsolToCut_ = NULL;
        delete [] colsol_;
        colsol_ = NULL;
    }
    else
    {
        si_->disableSimplexInterface();
    }
}

/** Compute normalization weights.*/
void
CglLandPSimplex::computeWeights(CglLandP::LHSnorm norm, CglLandP::Normalization type,
                                CglLandP::RhsWeightType rhs)
{
    norm_weights_.clear();
    norm_weights_.resize(ncols_orig_ ,1.);
    norm_weights_.resize(ncols_orig_ + nrows_orig_, 0.);

    double * rows_weights = &norm_weights_[ncols_orig_];
#ifndef INTEL_COMPILER
    std::vector<int> nnz(nrows_orig_,0);
#else
    std::vector<int> nnz(nrows_orig_);
#endif
    const CoinPackedMatrix * m = si_->getMatrixByCol();
    const double * val = m->getElements();
    const int * ind = m->getIndices();
    const int * length = m->getVectorLengths();
    const CoinBigIndex * start = m->getVectorStarts();

    rhs_weight_ = 1;

    if (type== CglLandP::WeightRHS)
    {
        if (rhs == CglLandP::Fixed)
            rhs_weight_ = (ncols_orig_ + 1);//params.rhsWeight;
        else if (rhs == CglLandP::Dynamic)
        {
            throw -1;
        }
    }

    if (norm == CglLandP::Infinity)
    {
        for (int i = 0 ; i < ncols_orig_ ; i++)
        {
            CoinBigIndex begin = start[i];
            CoinBigIndex end = begin + length[i];
            for (CoinBigIndex k = begin ; k < end ; k++)
            {
                rows_weights[ind[k]] = std::max(fabs(val[k]), rows_weights[ind[k]]);
                rhs_weight_ += fabs(val[k]);
                nnz[ind[k]] ++;
            }
        }
    }
    else if (norm == CglLandP::L1 ||
             norm == CglLandP::Average)
    {
        for (int i = 0 ; i < ncols_orig_ ; i++)
        {
            CoinBigIndex begin = start[i];
            CoinBigIndex end = begin + length[i];
            for (CoinBigIndex k = begin ; k < end ; k++)
            {
                rows_weights[ind[k]] += fabs(val[k]);
                nnz[ind[k]] ++;
            }
        }
        if (norm == CglLandP::Average)
        {
            for (int i = 0 ; i < nrows_orig_ ; i++)
            {
                rows_weights[i] = static_cast<double>(nnz[i]);
            }
        }
        if (type== CglLandP::WeightBoth)
        {
           rhs_weight_ += (ncols_orig_ + 1);
           std::cout<<"rhs_weight : "<<rhs_weight_<<std::endl;
        }
    }
    else if (norm == CglLandP::L2)
    {
        for (int i = 0 ; i < ncols_orig_ ; i++)
        {
            CoinBigIndex begin = start[i];
            CoinBigIndex end = begin + length[i];
            for (CoinBigIndex k = begin ; k < end ; k++)
            {
                rows_weights[ind[k]] += (val[k])*(val[k]);
                nnz[ind[k]] ++;
                rhs_weight_ += fabs(val[k]);
            }
        }
        for (int i = 0 ; i < nrows_orig_ ; i++)
        {
            rows_weights[i] = sqrt(rows_weights[i]);;
        }
        if (type== CglLandP::WeightBoth)
        {
          rhs_weight_ = (ncols_orig_ + 1);
        }
    }
    else if (norm == CglLandP::SupportSize)
    {
        for (int i = 0 ; i < ncols_orig_ ; i++)
        {
            CoinBigIndex begin = start[i];
            CoinBigIndex end = begin + length[i];
            for (CoinBigIndex k = begin ; k < end ; k++)
            {
                nnz[ind[k]] ++;
            }
        }

        for (int i = 0 ; i < nrows_orig_ ; i++)
        {
            rows_weights[i] = 1./ static_cast<double> (nnz[i]);
        }

       if (type== CglLandP::WeightBoth)
        {
          rhs_weight_ = (ncols_orig_ + 1);
        }

    }
    else if (norm ==CglLandP::Uniform)
    {
        for (int i = 0 ; i < nrows_orig_ ; i++)
        {
            rows_weights[i] = static_cast<double> (1);
        }
       if (type== CglLandP::WeightBoth)
        {
          rhs_weight_ = (ncols_orig_ + 1);
        }

    }

}
void
CglLandPSimplex::cacheUpdate(const CglLandP::CachedData &cached, bool reducedSpace)
{
    integers_ = cached.integers_;
    if (own_)
    {
        CoinCopyN(cached.basics_, nrows_orig_, basics_);
        CoinCopyN(cached.nonBasics_, ncols_orig_, nonBasics_);
        CoinCopyN(cached.colsol_, nrows_orig_ + ncols_orig_, colsol_);
        for (int i = 0 ; i < ncols_orig_ ; i++)
        {
            colsol_[nonBasics_[i]] = 0;
        }
        CoinCopyN(cached.colsol_, nrows_orig_ + ncols_orig_, colsolToCut_);
        //Zero all non basics in colsol setup the reduced space
        col_in_subspace.resize(0);
        col_in_subspace.resize(ncols_orig_+nrows_orig_,true);
        for (int i = 0 ; i < ncols_orig_ ; i++)
        {
            setColsolToCut(nonBasics_[i],0.);
            colsol_[nonBasics_[i]] = 0;
        }
        /** Mark the variables at zero in solution to cut so that we know that their contribution to reduced cost has to be computed*/
        if (reducedSpace)
        {
            for (int ii = 0; ii < ncols_orig_ ; ii++)
            {
                if (getColsolToCut(ii) - up_bounds_[ii] > 1e-08 || getColsolToCut(ii) - lo_bounds_[ii] < 1e-08)
                {
                    col_in_subspace[ii]=false;
                }
            }
        }
    }
    else
    {
        basics_ = cached.basics_;
        nonBasics_ = cached.nonBasics_;
    }
}

bool CglLandPSimplex::resetSolver(const CoinWarmStartBasis * /*basis*/)
{
    si_->disableSimplexInterface();
    return 0;
}

int
CglLandPSimplex::generateExtraCuts(const CglLandP::CachedData & cached,
                                   const CglLandP::Parameters& params)
{
    int ret_val = 0;
    for (int i = 0 ; i < nrows_ && cuts_.numberCuts() < params.extraCutsLimit; i++)
    {
        if (basics_[i] < ncols_)
            ret_val += generateExtraCut(i, cached, params);
    }
    return ret_val;
}

int
CglLandPSimplex::generateExtraCut(int i, const CglLandP::CachedData & cached,
                                  const CglLandP::Parameters& params)
{
    const int & iCol = basics_[i];
    if (!isInteger(iCol) || int_val(colsol_[iCol], params.away) ||
            !int_val(getColsolToCut(iCol), params.away) ||
            colsol_[iCol] < getLoBound(iCol) || colsol_[iCol] > getUpBound(iCol) ||
            (cuts_.rowCut(iCol) != NULL) )
    {
        return false;
    }

#ifdef DBG_OUT
    printf("var: %i, basic in row %i. integer=%i, colsol_=%f, colsolToCut_=%f.\n",iCol, i,
           isInteger(iCol), colsol_[iCol],
           getColsolToCut(iCol));

    printf("generating extra cut....\n");
#endif

    OsiRowCut * cut = new OsiRowCut;
    generateMig(i, *cut, params);
    assert(fabs(row_k_.rhs - colsol_[iCol]) < 1e-10);

    int code = validator_(*cut, cached.colsol_, *si_, params, &lo_bounds_[0], &up_bounds_[0]);
    if (code)
    {
        delete cut;
        return false;
    }
    else
    {
        cuts_.insert(iCol, cut);
        return true;
    }
}


void
CglLandPSimplex::genThisBasisMigs(const CglLandP::CachedData &cached,
                                  const CglLandP::Parameters & params)
{
    for (int i = 0 ; i < cached.nBasics_ ; i++)
    {
        const int iCol = basics_[i];
        if (iCol >= ncols_ ||
                !cached.integers_[iCol] ||
                int_val(colsol_[iCol], params.away))
            continue;
        OsiRowCut * cut = new OsiRowCut;
        generateMig(i, *cut, params);
        int code = validator_(*cut, cached.colsol_, *si_, params, &lo_bounds_[0], &up_bounds_[0]);
        if (code)
        {
            delete cut;
            continue;
        }
        cut->setEffectiveness(cut->violated(cached.colsol_));
        if (cuts_.rowCut(iCol) == NULL || cut->effectiveness() > cuts_.rowCut(iCol)->effectiveness())
        {
            cuts_.insert(iCol,cut);
        }
        else
            delete cut;
    }
}


bool
CglLandPSimplex::generateMig(int row, OsiRowCut & cut,
                             const CglLandP::Parameters & params)
{
    row_k_.num = row;
    pullTableauRow(row_k_);
    row_k_.rhs = row_k_.rhs - floor(row_k_.rhs);
    if (params.strengthen || params.modularize)
        createMIG(row_k_, cut);
    else
        createIntersectionCut(row_k_, cut);

    return 1;//At this point nothing failed, always generate a cut
}

bool
CglLandPSimplex::optimize
(int row, OsiRowCut & cut,const CglLandP::CachedData &cached,const CglLandP::Parameters & params)
{
    bool optimal = false;
    int nRowFailed = 0;

    double timeLimit = std::min(params.timeLimit, params.singleCutTimeLimit);
    timeLimit += CoinCpuTime();
    // double timeBegin = CoinCpuTime();
    int maximumCutLength = params.maximumCutLength;
    /** Copy the cached information */
    nrows_ = nrows_orig_;
    ncols_ = ncols_orig_;
    CoinCopyN(cached.basics_, nrows_, basics_);
    CoinCopyN(cached.nonBasics_, ncols_, nonBasics_);
    CoinCopyN(cached.colsol_, nrows_+ ncols_, colsol_);
    CoinCopyN(cached.colsol_, nrows_+ ncols_, colsolToCut_);

    delete basis_;
    basis_ = new CoinWarmStartBasis(*cached.basis_);
#define CACHED_SOLVER
#ifndef CACHED_SOLVER
    si_->enableSimplexInterface(0);
#else
    delete si_;
    si_ = cached.solver_->clone();
#ifdef CGL_HAS_OSICLP
    OsiClpSolverInterface * clpSi = getClpSolver(si_);
    OsiClpSolverInterface * clpSiRhs = getClpSolver(cached.solver_);
    if (CBC_SKIP_CLP_TEST||clpSi)
    {
        clp_ = clpSi;
	clpSi->getModelPtr()->copyEnabledStuff(clpSiRhs->getModelPtr());;
    }
#endif
#endif
#ifdef APPEND_ROW
    if (params.modularize)
    {
        append_row(row, params.modularize);
        row_k_.modularize(integers_);
    }
#endif
    for (int i = 0 ; i < ncols_; i++)
    {
        setColsolToCut(nonBasics_[i], 0.);
        colsol_[nonBasics_[i]] = 0;
    }

#ifdef APPEND_ROW
    if (params.modularize)
        row_k_.num = nrows_ - 1;
    else
#endif
        row_k_.num = row;

    pullTableauRow(row_k_);
    // give up if too many elements
    if (row_k_.getNumElements()>maximumCutLength)
      return false;
    row_k_.rhs = row_k_.rhs - floor(row_k_.rhs);

    if (params.modularize)
        row_k_.modularize(integers_);

    updateM1_M2_M3(row_k_, 0., params.perturb);
    sigma_ = computeCglpObjective(row_k_);

    handler_->message(Separating,messages_)<<basics_[row]<<sigma_<<CoinMessageEol<<CoinMessageEol;
    handler_->message(LogHead, messages_)<<CoinMessageEol<<CoinMessageEol;

    //Save the variable basic in this row
    //  int var_k = cached.basics_[k_];

    //Put a flag on each row to say if we want to continue trying to use it
    CoinFillN(rowFlags_,nrows_,true);

    int numberConsecutiveDegenerate = 0;
    bool allowDegeneratePivot = numberConsecutiveDegenerate < params.degeneratePivotLimit;
    bool beObstinate = 0;
    int numPivots = 0;
    int saveNumSourceEntered = numSourceRowEntered_;
    int saveNumIncreased = numIncreased_;
    int numCycle = 0;
    int numFailedPivots = 0;
    bool hasFlagedRow = false;
    int maxTryRow = 5;
    while (  !optimal && numPivots < params.pivotLimit)
    {
        if (timeLimit - CoinCpuTime() < 0.) break;

        updateM1_M2_M3(row_k_, 0., params.perturb);
        sigma_ = computeCglpObjective(row_k_);
        int direction = 0;
        int gammaSign = 0;
        int leaving = -1;
        int incoming = -1;
        double bestSigma;
        if (params.pivotSelection != CglLandP::initialReducedCosts || numPivots == 0)
        {
            leaving = fastFindCutImprovingPivotRow(direction, gammaSign, params.pivotTol,
                                                   params.pivotSelection == CglLandP::initialReducedCosts);
#if 0
            plotCGLPobj(direction, params.pivotTol, params.pivotTol, true, true, false);
            exit(1);
#endif
            if (leaving >= 0)
            {
                if (params.pivotSelection == CglLandP::mostNegativeRc ||
                        (params.pivotSelection == CglLandP::initialReducedCosts && numPivots == 0))
                {

                    if (params.pivotSelection == CglLandP::initialReducedCosts)
                        rowFlags_[leaving] = false;
                    incoming = fastFindBestPivotColumn(direction, gammaSign,
                                                       params.pivotTol, params.away,
                                                       (params.sepSpace==CglLandP::Fractional),
                                                       allowDegeneratePivot,
                                                       bestSigma, false//params.modularize
                                                      );
                    while (incoming < 0 && !optimal &&
                            nRowFailed < maxTryRow)   // if no improving was found rescan the tables of reduced cost to find a good one
                    {
                        if (incoming == -1 || params.countMistakenRc) nRowFailed ++;
                        rowFlags_[leaving] = false;
                        hasFlagedRow = true;
                        leaving = rescanReducedCosts(direction, gammaSign, params.pivotTol);
                        if (leaving >= 0)
                        {
                            incoming = fastFindBestPivotColumn(direction, gammaSign,
                                                               params.pivotTol,
                                                               params.away,
                                                               (params.sepSpace==CglLandP::Fractional),
                                                               allowDegeneratePivot,
                                                               bestSigma, false//params.modularize
                                                              );
                        }
                        else optimal = true;
                    }
                }
                else if (params.pivotSelection == CglLandP::bestPivot)
                {
                    incoming = findBestPivot(leaving, direction, params);
                }
            }
        }
        else if (params.pivotSelection == CglLandP::initialReducedCosts)
        {
            assert(numPivots > 0);
            while (incoming < 0 && !optimal)   // if no improving was found rescan the tables of reduced cost to find a good one
            {
                if (!hasFlagedRow)
                    hasFlagedRow = true;
                leaving = rescanReducedCosts(direction, gammaSign, params.pivotTol);
                rowFlags_[leaving] = false;
                if (leaving >= 0)
                {
                    incoming = fastFindBestPivotColumn(direction, gammaSign,
                                                       params.pivotTol,
                                                       params.away,
                                                       (params.sepSpace==CglLandP::Fractional),
                                                       allowDegeneratePivot,
                                                       bestSigma, false//params.modularize
                                                      );
                }
                else optimal = true;
            }
        }
        if (leaving >= 0)
        {
            if ( incoming >= 0 && !optimal)
            {
                if (inDegenerateSequence_)   //flag leaving row
                {
                    numberConsecutiveDegenerate++;
                    allowDegeneratePivot = numberConsecutiveDegenerate < params.degeneratePivotLimit;
                    rowFlags_[leaving] = false;
                }
                else
                {
                    beObstinate = 0;
                    numberConsecutiveDegenerate = 0;
                    allowDegeneratePivot = numberConsecutiveDegenerate < params.degeneratePivotLimit;
                }
                double gamma = - row_k_[nonBasics_[incoming]] / row_i_[nonBasics_[incoming]];
#ifdef CGL_HAS_OSICLP
                if (numPivots && ( numPivots % 40 == 0 ) && clp_)
                {
                    clp_->getModelPtr()->factorize();
                }
#endif

#ifndef OLD_COMPUTATION
                bool recompute_source_row = (numPivots && (numPivots % 10 == 0 ||
                                            fabs(gamma) < 1e-05));
#endif

                std::pair<int, int> cur_pivot(nonBasics_[incoming],basics_[leaving]);

                bool pivoted = changeBasis(incoming,leaving,direction,
#ifndef OLD_COMPUTATION
                                           recompute_source_row, 
#endif
                                           false);//params.modularize);

                if (leaving == row)
                {
                    numSourceRowEntered_++;
                }

                if (params.generateExtraCuts == CglLandP::WhenEnteringBasis &&
                        basics_[leaving] < ncols_ && cuts_.numberCuts() < params.extraCutsLimit)
                    generateExtraCut(leaving, cached, params);
                if (pivoted)
                {
                    numPivots++;

                    double lastSigma = sigma_;
                    if (params.modularize)
                    {
		      row_k_.modularize(integers_);
                    }
                    sigma_ = computeCglpObjective(row_k_);

                    if (sigma_ - lastSigma > -1e-4*(lastSigma))
                    {
                      if(sigma_ > 0) return 0;
#if 0
		      if (sigma_ > 0 || sigma_ - lastSigma > 1e1*(-lastSigma))
			return 0;
#endif
                    }


                    handler_->message(PivotLog,messages_)<<numPivots<<sigma_<<
                    nonBasics_[incoming]<<basics_[leaving]<<direction<<gamma<<inDegenerateSequence_<<CoinMessageEol<<CoinMessageEol;
                }
                else   //pivot failed
                {
                    numFailedPivots++;
                    //check wether sigma has changed if it has exit cut generation if it has not continue
                    double lastSigma = sigma_;
                    sigma_ = computeCglpObjective(row_k_);
                    if ( sigma_-lastSigma>1e-8)
                    {
                        handler_->message(PivotFailedSigmaIncreased,messages_)<<CoinMessageEol<<CoinMessageEol;
                        //break;
                        return 0;
                    }
                    handler_->message(PivotFailedSigmaUnchanged,messages_)<<CoinMessageEol<<CoinMessageEol;
                    numFailedPivots = params.failedPivotLimit + 1;
                    return 0;
                    if (numFailedPivots > params.failedPivotLimit)
                        break;

                }
            }
            else   //attained max number of leaving vars tries with no improvement
            {
                handler_->message(WarnGiveUpRow,messages_)<<nRowFailed<<CoinMessageEol<<CoinMessageEol;
                break;
            }
        }
        else
        {
            if (hasFlagedRow && beObstinate)
            {
                //Reset row flags
                CoinFillN(rowFlags_,nrows_,true);
                hasFlagedRow = false;
                if (inDegenerateSequence_)
                {
                    allowDegeneratePivot = false;
                    beObstinate = false;
                }
            }
            else
            {
                //could perturb but Ionut skipped that will see later
                optimal = true;
                handler_->message(FinishedOptimal, messages_)<<sigma_<<numPivots<<CoinMessageEol<<CoinMessageEol;
            }
        }
    }

    if (!optimal && numPivots >= params.pivotLimit)
    {
        std::string limit="pivots";
        handler_->message(HitLimit, messages_)<<limit<<numPivots<<CoinMessageEol<<CoinMessageEol;
    }
    if (!optimal && numFailedPivots >= params.failedPivotLimit)
    {
        std::string limit="failed pivots";
        handler_->message(HitLimit, messages_)<<limit<<numPivots<<CoinMessageEol<<CoinMessageEol;
    }
    // give up if too many elements
    if (row_k_.getNumElements()>maximumCutLength)
      return false;
    //Create the cut

    //pullTableauRow(row_k_);
    //row_k_.rhs = row_k_.rhs - floor(row_k_.rhs);

    //  double normalization = 100*normCoef(row_k_);
    {
        if (params.strengthen || params.modularize)
            createMIG(row_k_, cut);
        else
            createIntersectionCut(row_k_, cut);
    }

    if (params.generateExtraCuts == CglLandP::AtOptimalBasis)
    {
        generateExtraCuts(cached, params);
    }
    handler_->message(CutStat, messages_)<<row<<numPivots
    <<numSourceRowEntered_ - saveNumSourceEntered
    <<numIncreased_- saveNumIncreased
    <<numCycle<<CoinMessageEol;
    return 1;//At this point nothing failed, always generate a cut
}


bool
CglLandPSimplex::changeBasis(int incoming, int leaving, int leavingStatus,
#ifndef OLD_COMPUTATION
                             bool recompute_source_row,
#endif
                             bool modularize)
{
    double infty = si_->getInfinity();
    int clpLeavingStatus = leavingStatus;

#ifdef CGL_HAS_OSICLP
    if (clp_)
    {
        if (basics_[leaving] >= ncols_)
            clpLeavingStatus = - leavingStatus;
    }
#endif

    int code = 0;

    code = si_->pivot(nonBasics_[incoming],basics_[leaving], clpLeavingStatus);
    if (code)
    {
#ifdef OLD_COMPUTATION
        if (!modularize)
        {
            pullTableauRow(row_k_);
            row_k_.rhs = row_k_.rhs - floor(row_k_.rhs);
        }
	else{
	  int & indexLeaving = basics_[leaving];
	  if (leavingStatus==1)
	  {
	    setColsolToCut(indexLeaving, getUpBound(indexLeaving) - getColsolToCut(indexLeaving));
	  }
	else
	  {
	    setColsolToCut(indexLeaving, getColsolToCut(indexLeaving) + getLoBound(indexLeaving));
	  }
        }
#endif
        return 0;
    }
    numPivots_ ++;
    //swap bounds
    int & indexLeaving = basics_[leaving];
#ifdef OLD_COMPUTATION
    if (!modularize)
      {
	if (leavingStatus==1)
	  {
	    setColsolToCut(indexLeaving, getUpBound(indexLeaving) - getColsolToCut(indexLeaving));
	  }
	else
	  {
	    setColsolToCut(indexLeaving, getColsolToCut(indexLeaving) - getLoBound(indexLeaving));
	  }
      }
#endif
    
    if (indexLeaving < ncols_)
      {
	basis_->setStructStatus(indexLeaving, leavingStatus==1 ? CoinWarmStartBasis::atUpperBound : CoinWarmStartBasis::atLowerBound);
      }
    else
      {
	int iRow = basics_[leaving] - ncols_;
	basis_->setArtifStatus(iRow,  leavingStatus==1 ? CoinWarmStartBasis::atUpperBound : CoinWarmStartBasis::atLowerBound);
	//    assert(leavingStatus==-1 || (rowLower_[iRow]>-1e50 && rowUpper_[iRow] < 1e50));
      }
    
    if (nonBasics_[incoming] < ncols_)
      {
	int & indexIncoming = nonBasics_[incoming];
	CoinWarmStartBasis::Status status = basis_->getStructStatus(indexIncoming);
	if (status==CoinWarmStartBasis::atUpperBound)
	  setColsolToCut(indexIncoming, getUpBound(indexIncoming) - getColsolToCut(indexIncoming));
	else
	  setColsolToCut(indexIncoming, getColsolToCut(indexIncoming) + getLoBound(indexIncoming));
	basis_->setStructStatus(indexIncoming, CoinWarmStartBasis::basic);
      }
    else
      {
	int iRow = nonBasics_[incoming] - ncols_;
	int & indexIncoming = nonBasics_[incoming];
	
	if (basis_->getArtifStatus(iRow)==CoinWarmStartBasis::atUpperBound)
	  setColsolToCut(indexIncoming, getUpBound(indexIncoming) - getColsolToCut(indexIncoming));
	else
	  setColsolToCut(indexIncoming, getColsolToCut(indexIncoming) + getLoBound(indexIncoming));
	
	basis_->setArtifStatus(iRow,  CoinWarmStartBasis::basic);
      }
    
    int swap = basics_[leaving];
    basics_[leaving] = nonBasics_[incoming];
    nonBasics_[incoming] = swap;
    //update solution of leaving variable
    colsol_[nonBasics_[incoming]] = 0;
    
    //update solution for basics
    const double * lpSol = si_->getColSolution();
    const double * rowAct = si_->getRowActivity();
    const double * rowLower = si_->getRowLower();
    const double * rowUpper = si_->getRowUpper();

    for (int i = 0 ; i < nrows_ ; i++)
    {
        int& iCol = basics_[i];
        if (iCol<ncols_)
            colsol_[iCol] = lpSol[iCol];
        else   // Osi does not give direct acces to the value of slacks
        {
            int iRow = iCol - ncols_;
            colsol_[iCol] = - rowAct[iRow];
            if (rowLower[iRow]> -infty)
            {
                colsol_[iCol] += rowLower[iRow];
            }
            else
            {
                colsol_[iCol] += rowUpper[iRow];
            }
        }
    }

    // basics_ may unfortunately change reload
    int k = basics_[row_k_.num];
    si_->getBasics(basics_);
    if (basics_[row_k_.num] != k)
    {
        for (int ii = 0 ; ii < nrows_ ; ii++)
        {
            if (basics_[ii]==k)
            {
                row_k_.num= ii;
                break;
            }
        }
    }

#ifndef OLD_COMPUTATION
    if (!modularize && recompute_source_row)
    {
#else
    if (!modularize)
    {
#endif
        pullTableauRow(row_k_);
        row_k_.rhs =  row_k_.rhs - floor(row_k_.rhs);
    }
#if 0
}
#endif
else //Update row k by hand
{
    double gamma = - row_k_[basics_[leaving]] / row_i_[basics_[leaving]];
    row_k_[basics_[leaving]] = 0;
    row_k_.quickAdd(nonBasics_[incoming], gamma);
    if(1 || fabs(gamma) > 1e-9){
      int nnz = row_i_.getNumElements();
      const int * indices = row_i_.getIndices();
      for (int i = 0 ; i < nnz; i++)
	{
	  if(row_k_.getNumElements() > row_k_.capacity() - 2) row_k_.scan();
	  if (indices[i] != nonBasics_[incoming] && indices[i] != basics_[leaving])
	    row_k_.quickAdd(indices[i], gamma * row_i_[indices[i]]);
	}
      row_k_.rhs += gamma * row_i_.rhs;
    }
    row_k_.scan();
    row_k_.clean(1e-10);
    checkVec(row_k_);
#if 0
    TabRow test_row(this);
    int rowsize = ncols_orig_ + nrows_orig_ + 1;
    test_row.reserve(rowsize);
    test_row.num = row_k_.num;
    pullTableauRow(test_row);
    test_row.rhs =  test_row.rhs - floor(test_row.rhs);
    VecModEqAssert(row_k_, test_row);
#endif
}
    
    return true;
}

/** Find a row which can be used to perform an improving pivot return index of the cut or -1 if none exists
 * (i.e., find the leaving variable).*/
int
CglLandPSimplex::findCutImprovingPivotRow( int &direction, int &gammaSign, double tolerance)
{
    bool bestRed = 0;
    tolerance = -10*tolerance;
    int bestRow = -1;
    int bestDirection = 0;
    int bestGamma = 0;
    double infty = si_->getInfinity();
    for (row_i_.num = 0 ; row_i_.num < nrows_; row_i_.num++)
    {
        //if ( (row_k_.modularized_ || row_i_.num != row_k_.num)//obviously not necessary to combine row k with itself (unless modularized)
        if ( (row_i_.num != row_k_.num)//obviously not necessary to combine row k with itself (unless modularized)
                && rowFlags_[row_i_.num] //row has not been flaged
                //   && fabs(getUpBound(basics_[row_i_.num]) - getLoBound(basics_[row_i_.num]))>1e-09 //variable is not fixed
           )
        {
            pullTableauRow(row_i_);
            double tau = computeRedCostConstantsInRow();

            if (getLoBound(basics_[row_i_.num]) > -infty)
                // variable can leave at its lower bound
                //Compute reduced cost with basics_[i] value decreasing
            {
                direction = -1;

                gammaSign = -1;
                double redCost = computeCglpRedCost(direction, gammaSign, tau);
                if (redCost<tolerance)
                {
                    if (bestRed)
                    {
                        tolerance = redCost;
                        bestRow = row_i_.num;
                        bestDirection = direction;
                        bestGamma = gammaSign;
                    }
                    else return row_i_.num;
                }
                gammaSign = 1;
                redCost = computeCglpRedCost(direction, gammaSign, tau);
                if (redCost<tolerance)
                {
                    if (bestRed)
                    {
                        tolerance = redCost;
                        bestRow = row_i_.num;
                        bestDirection = direction;
                        bestGamma = gammaSign;
                    }
                    else return row_i_.num;
                }
            }
            if ( getUpBound(basics_[row_i_.num])<infty) // variable can leave at its upper bound
                //Compute reduced cost with basics_[i] value decreasing
            {
                direction = 1;
                //         adjustTableauRow(i_, row_i_, rhs_i_, direction);
                gammaSign = -1;
                double redCost = computeCglpRedCost(direction, gammaSign, tau);
                if (redCost<tolerance)
                {
                    if (bestRed)
                    {
                        tolerance = redCost;
                        bestRow = row_i_.num;
                        bestDirection = direction;
                        bestGamma = gammaSign;
                    }
                    else return row_i_.num;
                }
                gammaSign = 1;
                redCost = computeCglpRedCost(direction, gammaSign, tau);
                if (redCost<tolerance)
                {
                    if (bestRed)
                    {
                        tolerance = redCost;
                        bestRow = row_i_.num;
                        bestDirection = direction;
                        bestGamma = gammaSign;
                    }
                    else return row_i_.num;
                }
            }
            rowFlags_[row_i_.num]=false;
        }
    }
    direction = bestDirection;
    gammaSign = bestGamma;
    row_i_.num=bestRow;
    if (row_i_.num >=0 && row_i_.num != bestRow)
    {
        row_i_.num=bestRow;
        pullTableauRow(row_i_);
    }
    assert (bestRow<0||direction!=0); 
    return bestRow;
}


/** Find a row which can be used to perform an improving pivot the fast way
 * (i.e., find the leaving variable).
 \return index of the cut or -1 if none exists. */
int
CglLandPSimplex::fastFindCutImprovingPivotRow( int &direction, int &gammaSign,
        double tolerance, bool flagPositiveRows)
{
    bool modularize = false;
    double sigma = sigma_ /rhs_weight_;
    //Fill vector to compute contribution to reduced cost of variables in M1 and M2 (nz non-basic vars in row_k_.row).
    // 1. Put the values
    // 2. Post multiply by basis inverse
    double * rWk1bis_ =NULL;
    CoinFillN(&rWk1_[0],nrows_,static_cast<double> (0));
    if (modularize)
        CoinFillN(rWk1bis_, nrows_, static_cast<double> (0));
    int capacity = 0;
    const CoinPackedMatrix* mat = si_->getMatrixByCol();

    const CoinBigIndex* starts = mat->getVectorStarts();
    const int * lengths = mat->getVectorLengths();
    const int * indices = mat->getIndices();
    const double * elements = mat->getElements();

    for (unsigned int i = 0 ;  i < M1_.size() ; i++)
    {
        const int& ii = M1_[i];
        if (ii < ncols_)
        {
            const CoinBigIndex& begin = starts[ii];
            const CoinBigIndex end = begin + lengths[ii];
            bool swap = false;
            if (basis_->getStructStatus(ii)==CoinWarmStartBasis::atUpperBound) swap = true;
            for (CoinBigIndex k = begin ; k < end ;  k++)
            {
                if (swap)
                {
                    rWk1_[indices[k]] += normedCoef(elements[k] * sigma, ii);
                    if (modularize)
                    {
                        rWk1bis_[indices[k]] += elements[k] * (getColsolToCut(ii) - normedCoef(sigma, ii));
                    }

                }
                else
                {
                    rWk1_[indices[k]] -= normedCoef(elements[k] * sigma, ii);
                    if (modularize)
                    {
                        rWk1bis_[indices[k]] -= elements[k] * (getColsolToCut(ii) - normedCoef(sigma, ii));
                    }
                }
            }
        }
        else
        {
            bool swap = false;
            if (basis_->getArtifStatus(ii - ncols_orig_)==CoinWarmStartBasis::atUpperBound) swap = true;
            if (swap)
            {
                rWk1_[ii - ncols_] += normedCoef(sigma, ii);
                if (modularize)
                {
                    rWk1bis_[ii - ncols_] += (getColsolToCut(ii) - normedCoef(sigma, ii));
                }
            }
            else
            {
                rWk1_[ii - ncols_] -= normedCoef(sigma, ii);
                if (modularize)
                {
                    rWk1bis_[ii - ncols_] -= (getColsolToCut(ii) - normedCoef(sigma, ii));
                }
            }
        }
    }
    for (unsigned int i = 0 ;  i < M2_.size(); i++)
    {
        const int& ii = M2_[i];
        if (ii<ncols_)
        {
            const CoinBigIndex& begin = starts[ii];
            const CoinBigIndex end = begin + lengths[ii];
            bool swap = false;
            if (basis_->getStructStatus(ii)==CoinWarmStartBasis::atUpperBound) swap = true;
            for (CoinBigIndex k = begin ; k < end ;  k++)
            {
                if (swap)
                {
                    rWk1_[indices[k]] += elements[k] * (getColsolToCut(ii) - normedCoef(sigma, ii));
                    if (modularize)
                        rWk1bis_[indices[k]] += elements[k] * normedCoef(sigma, ii);
                }
                else
                {
                    rWk1_[indices[k]] -= elements[k] * (getColsolToCut(ii) - normedCoef(sigma, ii));
                    if (modularize)
                        rWk1bis_[indices[k]] -= elements[k] * normedCoef(sigma, ii);
                }
            }
        }
        else
        {
            bool swap = false;
            if (basis_->getArtifStatus(M2_[i] - ncols_orig_)==CoinWarmStartBasis::atUpperBound) swap = true;
            if (swap)
            {
                rWk1_[ii - ncols_] += (getColsolToCut(ii) - normedCoef(sigma, ii));
                if (modularize)
                    rWk1bis_[ii - ncols_] += normedCoef(sigma, ii);
            }
            else
            {
                rWk1_[ii - ncols_] -= (getColsolToCut(ii) - normedCoef(sigma, ii));
                if (modularize)
                    rWk1bis_[ii - ncols_] -= normedCoef(sigma, ii);
            }
        }
    }

    for (int i = 0 ; i < nrows_ ; i++)
    {
        if (rWk1_[i])
            rIntWork_[capacity++] = i;
    }
    CoinIndexedVector indexed;
    indexed.borrowVector(nrows_, capacity, &rIntWork_[0], &rWk1_[0]);

#ifdef CGL_HAS_OSICLP
    if (clp_)
        clp_->getBInvACol(&indexed);
    else
#endif
        throw CoinError("Function not implemented in this OsiSolverInterface",
                        "getBInvACol","CglLandpSimplex");
    indexed.returnVector();
    if (modularize)
    {
        capacity = 0;
        for (int i = 0 ; i < nrows_ ; i++)
        {
            if (rWk1bis_[i])
                rIntWork_[capacity++] = i;
        }

        indexed.borrowVector(nrows_, capacity, &rIntWork_[0], rWk1bis_);
#ifdef CGL_HAS_OSICLP
        if (clp_)
            clp_->getBInvACol(&indexed);
        else
#endif
            indexed.returnVector();
    }
    //Now compute the contribution of the variables in M3_
    //Need to get the column of the tableau in rW3_ for each of these and
    //add up with correctly in storage for multiplier for negative gamma (named rW3_) and
    //for positive gamma (which is named rW4_)
    if (!M3_.empty())
    {
        CoinFillN(&rWk3_[0],nrows_,0.);
        CoinFillN(&rWk4_[0],nrows_,0.);
        if (modularize)
        {
            double * rWk3bis_ = NULL;
            double * rWk4bis_ = NULL;
            CoinFillN(rWk3bis_,nrows_,0.);
            CoinFillN(rWk4bis_,nrows_,0.);
        }
    }
    for (unsigned int i = 0 ; i < M3_.size() ; i++)
    {
        const int & ii = M3_[i];
        si_->getBInvACol(ii, &rWk2_[0]);
        bool swap = false;
        if (ii < ncols_orig_ && basis_->getStructStatus(ii)==CoinWarmStartBasis::atUpperBound) swap = true;
        if (ii >= ncols_orig_ && basis_->getArtifStatus(ii - ncols_orig_)==CoinWarmStartBasis::atUpperBound) swap = true;

        for (int j = 0 ; j < nrows_ ; j++)
        {
            if (swap)
                rWk2_[j] = - rWk2_[j];
            if (rWk2_[j] > 0.)
            {
                //is in M1 for multiplier with negative gamma
                rWk3_[j] -= normedCoef(sigma*rWk2_[j], ii);
                //is in M2 for multiplier with positive gamma
                rWk4_[j] -= (getColsolToCut(ii) - normedCoef(sigma, ii))*rWk2_[j];
            }
            else if (rWk2_[j] < 0.)
            {
                //is in M2 for multiplier with negative gamma
                rWk3_[j] -=(getColsolToCut(ii) - normedCoef(sigma, ii))*rWk2_[j];

                //is in M1 for multiplier with positive gamma
                rWk4_[j] -= normedCoef(sigma, ii)*rWk2_[j];
            }
        }

    }
    //Now, we just need to add up everything correctly for each of the reduced
    //cost. Compute the Tau in rWk2_ which is not used anymore then compute the reduced cost in rWk1_ for u^l_j
    // rwk2_ in u^u_j rWk3_ for v^u_j rWk4_ for v^u_j
    // Let's rename not to get too much confused
    double * ul_i = &rWk1_[0];
    double * uu_i = &rWk2_[0];
    double * vl_i = &rWk3_[0];
    double * vu_i = &rWk4_[0];
    int bestRow = -1;
    int bestDirection = 0;
    int bestGammaSign = 0;
    // double infty = si_->getInfinity();


    nNegativeRcRows_ = 0;//counter
    int nZeroRc = 0;
    int nPositiveRc = 0;

    double fzero = getColsolToCut(basics_[row_k_.num]) - floor(getColsolToCut(basics_[row_k_.num]));
    //    fzero = row_k_.rhs;
    //for (int i = 0 ; i < ncols_orig_ ; i++) {
    //  fzero -= getColsolToCut(nonBasics_[i]) * row_k_[nonBasics_[i]];
    //}

    double bestReducedCost = -tolerance;
    for (int i = 0 ; i < nrows_ ; i++)
    {
      //if ((!row_k_.modularized_ && i == row_k_.num)//obviously not necessary to combine row k with itself
        if ((i == row_k_.num)//obviously not necessary to combine row k with itself
                //   && fabs(getUpBound(basics_[row_i_.num]) - getLoBound(basics_[row_i_.num]))>1e-09 //variable is not fixed
                || col_in_subspace[basics_[i]] == false
           )
        {
            ul_i[i]=uu_i[i]=vl_i[i]=vu_i[i]=10.;
            rowFlags_[i] = false;
            continue;
        }

        double tau1 = rWk1_[i];
        double tau2 = rWk1_[i];
        double tau3 = rWk1_[i];
        double tau4 = rWk1_[i];
        if (!M3_.empty())
        {
            tau1 += rWk3_[i];
            tau2 += rWk4_[i];
            tau3 += rWk4_[i];
            tau4 += rWk3_[i];
        }
        if (modularize)
        {
            tau1 = rWk1_[i] + rWk3_[i];
            tau2 = rWk1bis_[i] + rWk3_[i];
            tau3 = - rWk1_[i] - rWk4_[i];
            tau4 = - rWk1bis_[i] - rWk3_[i];
        }


        double redCost;
        bool hasNegativeRc = false;

        double loBound = getLoBound(basics_[i]);

        if (loBound > -1e50)
        {
            redCost = - normedCoef(sigma,basics_[i]) + (tau1)
                      + (1 - fzero) * ( colsol_[basics_[i]]
                                        - loBound);


            if (redCost < -tolerance)
            {
                ul_i[i] = redCost;
                hasNegativeRc = true;
            }
            else
            {
                if (fabs(redCost) < tolerance) nZeroRc++;
                else nPositiveRc ++;
                ul_i[i] = 10.;
            }
            if (redCost < bestReducedCost
                    && rowFlags_[i] )   //row has not been flaged
            {
                bestDirection = -1;
                bestGammaSign = -1;
                bestReducedCost = redCost;
                bestRow = i;
            }


            redCost = -normedCoef(sigma,basics_[i]) - (tau2)
                      - (1 - fzero) * ( colsol_[basics_[i]]
                                        - loBound)
                      - loBound + getColsolToCut(basics_[i]);

            if (redCost < -tolerance)
            {
                vl_i[i] = redCost;
                hasNegativeRc = true;
            }
            else
            {
                if (fabs(redCost) < tolerance) nZeroRc++;
                else nPositiveRc ++;
                vl_i[i]=10.;
            }

            if (redCost < bestReducedCost
                    && rowFlags_[i])   //row has not been flaged
            {
                bestDirection = -1;
                bestGammaSign = 1;
                bestReducedCost = redCost;
                bestRow = i;
            }



        }
        else
        {
            ul_i[i] = 10.;
            vl_i[i] = 10.;
        }
        double upBound = getUpBound(basics_[i]);
        if (getUpBound(basics_[i]) < 1e50)
        {
            redCost = - normedCoef(sigma,basics_[i]) - (tau3)
                      + (1 - fzero) * ( - colsol_[basics_[i]]
                                        + upBound);

            if (redCost < -tolerance)
            {
                uu_i[i] = redCost;
                hasNegativeRc = true;
            }
            else
            {
                if (fabs(redCost) < tolerance) nZeroRc++;
                else nPositiveRc ++;
                uu_i[i] = 10.;
            }

            if (redCost < bestReducedCost
                    && rowFlags_[i])   //row has not been flaged
            {
                bestDirection = 1;
                bestGammaSign = -1;
                bestReducedCost = redCost;
                bestRow = i;
            }

            redCost = -normedCoef(sigma,basics_[i]) + (tau4)
                      - (1 - fzero) * ( - colsol_[basics_[i]]
                                        + upBound)
                      + upBound - getColsolToCut(basics_[i]);

            if (redCost < -tolerance)
            {
                vu_i[i] = redCost;
                hasNegativeRc = true;
            }


            else
            {
                if (fabs(redCost) < tolerance) nZeroRc++;
                else nPositiveRc ++;
                vu_i[i] = 10.;
            }

            if (redCost < bestReducedCost
                    && rowFlags_[i])   //row has not been flaged
            {
                bestDirection = 1;
                bestGammaSign = 1;
                bestReducedCost = redCost;
                bestRow = i;
            }
        }
        else
        {
            uu_i[i] = 10.;
            vu_i[i] = 10.;
        }
        if (hasNegativeRc) nNegativeRcRows_ ++;
        else if (flagPositiveRows) rowFlags_[i] = false;
    }
    handler_->message(NumberNegRc, messages_)<<nNegativeRcRows_<<CoinMessageEol;
    handler_->message(NumberZeroRc, messages_)<<nZeroRc<<CoinMessageEol;
    handler_->message(NumberPosRc, messages_)<<nPositiveRc<<CoinMessageEol;
    //  throw -1;
    direction = bestDirection;
    gammaSign = bestGammaSign;
    //  row_i_.num=bestRow;
    if (bestRow != -1)
    {
        chosenReducedCostVal_ = bestReducedCost;
        row_i_.num=bestRow;
        pullTableauRow(row_i_);
        handler_->message(FoundImprovingRow, messages_)<<
        bestRow<<basics_[bestRow]<<direction<<gammaSign<<bestReducedCost
        <<CoinMessageEol;
    }
    assert (bestRow<0||direction!=0); 
    return bestRow;
}


/** Find a row which can be used to perform an improving pivot tables are already filled.
 \return index of the cut or -1 if none exists. */
int
CglLandPSimplex::rescanReducedCosts( int &direction, int &gammaSign, double tolerance)
{
    // The reduced cost are already here in rWk1_ is u^l_j
    // rwk2_ is u^u_j rWk3_ is v^u_j rWk4_ is v^u_j
    // Let's rename not to get too much confused
    double * ul_i = &rWk1_[0];
    double * uu_i = &rWk2_[0];
    double * vl_i = &rWk3_[0];
    double * vu_i = &rWk4_[0];
    int bestRow = -1;
    int bestDirection = 0;
    int bestGammaSign = 0;
    // double infty = si_->getInfinity();
    double bestReducedCost = -tolerance;
    for (int i = 0 ; i < nrows_ ; i++)
    {
        if (i == row_k_.num//obviously not necessary to combine row k with itself
                || !rowFlags_[i] //row has not been flaged
                //   && fabs(getUpBound(basics_[row_i_.num]) - getLoBound(basics_[row_i_.num]))>1e-09 //variable is not fixed
           )
            continue;

        if (ul_i[i] < bestReducedCost
                && rowFlags_[i])   //row has not been flaged
        {
            bestDirection = -1;
            bestGammaSign = -1;
            bestReducedCost = ul_i[i];
            bestRow = i;
        }

        if (vl_i[i] < bestReducedCost
                && rowFlags_[i])   //row has not been flaged
        {
            bestDirection = -1;
            bestGammaSign = 1;
            bestReducedCost = vl_i[i];
            bestRow = i;
        }


        if (uu_i[i] < bestReducedCost
                && rowFlags_[i])   //row has not been flaged
        {
            bestDirection = 1;
            bestGammaSign = -1;
            bestReducedCost = uu_i[i];
            bestRow = i;
        }
        if (vu_i[i] < bestReducedCost
                && rowFlags_[i])   //row has not been flaged
        {
            bestDirection = 1;
            bestGammaSign = 1;
            bestReducedCost = vu_i[i];
            bestRow = i;
        }

    }
    direction = bestDirection;
    gammaSign = bestGammaSign;
    //  row_i_.num=bestRow;
    if (bestRow != -1)
    {
        chosenReducedCostVal_ = bestReducedCost;
        row_i_.num=bestRow;
        pullTableauRow(row_i_);
        handler_->message(FoundImprovingRow, messages_)<<
        bestRow<<basics_[bestRow]<<direction<<gammaSign<<bestReducedCost
        <<CoinMessageEol;
    }
    assert (bestRow<0||direction!=0); 
    return bestRow;
}


void
CglLandPSimplex::compute_p_q_r_s(double gamma, int gammaSign, double &p, double & q, double & r , double &s)
{
    for (int i = 0 ; i < ncols_ ; i++)
    {
        const int &ii = nonBasics_[i];//True index
        if (colCandidateToLeave_[i]==false) continue;
        const double& val = getColsolToCut(ii); //value in solution to cut
        const double& row_k = row_k_[ii]; // coefficient in row k
        const double& row_i = row_i_[ii]; // coefficient in row i
        double coeff = row_k + gammaSign * gamma*row_i;
        if (coeff>0.)
        {
            if (gammaSign > 0)
            {
                p += row_k * val;
            }
            else
            {
                //        if(fabs(getColsolToCut(nonBasics_[i])) > 0)
                {
                    p += row_k * val;
                    q += row_i * val;
                }
            }
            r += normedCoef(row_k, ii) ;
            s += normedCoef(row_i, ii);
        }
        else if (coeff< 0.)
        {
            if (gammaSign > 0)
            {
                q -= row_i * val;
            }
            r -= normedCoef(row_k,ii);
            s -= normedCoef(row_i,ii);
        }
        else
        {
            if (gammaSign > 0 && row_i < 0)
            {
                q -= row_i * val;
            }
            else if (gammaSign < 0 && row_i < 0)
            {
                q += row_i * val;
            }
            s += normedCoef(gammaSign*fabs(row_i),ii);
        }
    }
}

//  double f_plus(double gamma, int gammaSign){
//  double num = row_k_.
//  for(int i = 0 ; i < ncols_ ; i++){
//  }
//}

/** Find the column which leads to the best cut (i.e., find incoming variable).*/
int
CglLandPSimplex::fastFindBestPivotColumn(int direction, int gammaSign,
        double pivotTol, double rhsTol,
        bool reducedSpace, bool allowDegenerate,
        double & bestSigma, bool modularize)
{
    gammas_.clear();
    pivotTol = 1e-05;
    adjustTableauRow(basics_[row_i_.num], row_i_, direction);

    double fzero = getColsolToCut(basics_[row_k_.num]) - floor(getColsolToCut(basics_[row_k_.num]));

 

    double p = 0;
    double q = 0;
    if(!modularize){//Take a shortcut
      p = -row_k_.rhs * (1 - fzero);
      q = row_i_.rhs * fzero;
      
      if (gammaSign < 0){
	q -= row_i_.rhs;
      }
    }
    double r = 1.;
    double s = normedCoef( static_cast<double> (gammaSign), basics_[row_i_.num]);

    //bool haveSmallGammaPivot = false;
    double gammaTolerance = 0;
    if (allowDegenerate)
        gammaTolerance = 0;
    //fill the array with the gammas of correct sign
    for (int i = 0 ; i < ncols_ ; i++)
    {
        const int &ii = nonBasics_[i];//True index
        const double& val = getColsolToCut(ii); //value in solution to cut
        const double& row_k = row_k_[ii]; // coefficient in row k
        const double& row_i = row_i_[ii]; // coefficient in row i
        if(modularize){
	  p-=row_k_.rhs*row_k*val;
	  q-=row_i_.rhs*row_k*val;
	}

        if (reducedSpace && colCandidateToLeave_[i]==false)
        {
            assert(col_in_subspace[ii]==false);
            continue;
        }
        double gamma = 1;
        if (fabs(row_i) > gammaTolerance && fabs(row_k) > gammaTolerance)
        {
            gamma = - row_k/row_i;
            if (gamma * gammaSign > gammaTolerance)
            {
                gammas_.insert(i,gamma*gammaSign);
            }
        }
        gamma = fabs(gamma); //  we already know the sign of gamma, its absolute value is more usefull
        if (row_k>gammaTolerance)
        {
            if (gammaSign > 0)
            {
                p += row_k * val;
            }
            else
            {
                {
                    p += row_k * val;
                    q += row_i * val;
                }
            }
            r += normedCoef(row_k, ii) ;
            s += normedCoef(row_i, ii);
        }
        else if (row_k< gammaTolerance)
        {
            if (gammaSign > 0)
            {
                q -= row_i * val;
            }
            r -= normedCoef(row_k,ii);
            s -= normedCoef(row_i,ii);
        }
        else
        {
            //haveSmallGammaPivot |= true;
            if (gammaSign > 0 && row_i < 0)
            {
                q -= row_i * val;
            }
            else if (gammaSign < 0 && row_i < 0)
            {
                q += row_i * val;
            }
            s += normedCoef(gammaSign*fabs(row_i),ii);
        }
    }


    if(modularize){
      p -= row_k_.rhs * (1 - row_k_.rhs);
      q += row_i_.rhs * row_k_.rhs;
      if (gammaSign < 0){
	q -= row_i_.rhs;
      }
    }

    int n = gammas_.getNumElements();
    if (n==0)
    {
        resetOriginalTableauRow(basics_[row_i_.num], row_i_, direction);
        return -2;
    }
    gammas_.sortIncrElement();
    const int* inds = gammas_.getIndices();
    const double * elements = gammas_.getElements();
    int bestColumn = -1;
    double newSigma = 1e100;
    DblEqAssert(sigma_, rhs_weight_*p/r);
    bestSigma = sigma_ = rhs_weight_*p/r;
    int lastValid = -1;
#ifndef NDEBUG
    bool rc_positive=false;
    if (M3_.size())
        DblEqAssert( gammaSign*(q * r - p * s)/r, chosenReducedCostVal_);
#endif
    if ( gammaSign*(q * r - p * s) >= 0)
    {
        // after recomputing reduced cost (using exact row) it is found to be >=0
        resetOriginalTableauRow(basics_[row_i_.num], row_i_, direction);
        return -2;
#ifndef NDEBUG
        rc_positive = true;
#endif
    }
    for (int i = 0 ; i < n ; i++)
    {
        double newRhs = row_k_.rhs + gammaSign * elements[i] * row_i_.rhs;
        if (newRhs < rhsTol || newRhs > 1 - rhsTol)
        {
            //	if(i == 0)
            break;
        }
        newSigma = (p + gammaSign * elements[i] * q)*rhs_weight_/(r + gammaSign*elements[i] * s);
#ifdef DEBUG_LAP
        double alt = computeCglpObjective(gammaSign*elements[i], false);
        DblEqAssert(newSigma, alt);
#endif
        if (newSigma > bestSigma - 1e-08*bestSigma)
        {
#ifndef NDEBUG

            if (!rc_positive)
            {
            }


            if (0 && elements[i] <= 1e-05)
            {
            }
#endif
            break;
        }
        else if (newSigma <= bestSigma)  // && colCandidateToLeave_[inds[i]])
        {
            bestColumn = inds[i];
            bestSigma = newSigma;
            lastValid = i;
        }


#ifndef NDEBUG
        if (rc_positive)
        {
            break;
        }
#endif
        int col = nonBasics_[inds[i]];
        if (row_i_[col] *gammaSign > 0)
        {
            p += row_k_[col] * getColsolToCut(col);
            q += row_i_[col] * getColsolToCut(col);
            r += normedCoef(row_k_[col]*2,col);
            s += normedCoef(row_i_[col]*2,col);
        }
        else
        {
            p -= row_k_[col] * getColsolToCut(col);
            q -= row_i_[col] * getColsolToCut(col);
            r -= normedCoef(row_k_[col]*2,col);
            s -= normedCoef(row_i_[col]*2,col);
        }
        if (gammaSign*(q * r - p * s) >= 0)   /* function is starting to increase stop here*/
        {
#ifndef NDEBUG
            if (0 && elements[i] <= 1e-5)
            {
            }
            rc_positive = true;
#endif
            break;
        }

    }
//Get the results did we find a valid pivot? is it degenerate?
    if (bestColumn == -1)   // Apparently no pivot is within the tolerances
    {
        resetOriginalTableauRow(basics_[row_i_.num], row_i_, direction);
        handler_->message(WarnFailedPivotTol, messages_)<<CoinMessageEol<<CoinMessageEol;
        return -1;
    }

    if (fabs(row_i_[nonBasics_[bestColumn]]) < 1e-05)   // Apparently no pivot is within the tolerances
    {
        resetOriginalTableauRow(basics_[row_i_.num], row_i_, direction);
        handler_->message(WarnFailedPivotTol, messages_)<<CoinMessageEol<<CoinMessageEol;
        return -2;
    }

    //std::cout<<"Minimum of f attained at breakpoint "<<lastValid<<std::endl;
    assert(bestSigma <= sigma_);
#ifdef DEBUG_LAP
    bestSigma_ = bestSigma;
    double otherSigma= computeCglpObjective(gammaSign*elements[lastValid], false);
    DblEqAssert(bestSigma, otherSigma);
    DblEqAssert(bestSigma, computeCglpObjective(gammaSign*elements[lastValid], false, new_row_));
    DblEqAssert(bestSigma, computeCglpObjective(new_row_));
    assert(row_k_.modularized_ || row_i_.num != row_k_.num);
#endif

#ifdef OLD_COMPUTATION
    if (!modularize)
        resetOriginalTableauRow(basics_[row_i_.num], row_i_, direction);
#endif
    if (bestSigma < sigma_ - 1e-07)   //everything has gone ok
    {
        handler_->message(FoundBestImprovingCol, messages_)<<nonBasics_[bestColumn]<<gammaSign * elements[lastValid]<<bestSigma<<CoinMessageEol<<CoinMessageEol;
        inDegenerateSequence_ = false;
        assert (bestColumn<0||direction!=0); 
        return bestColumn;
    }
    else if (allowDegenerate)   //Pivot is degenerate and we allow
    {
        inDegenerateSequence_ = true;
        assert (bestColumn<0||direction!=0); 
        return bestColumn;
    }
    else   //we don't accept a degenerate pivot
    {
        handler_->message(WarnFailedBestImprovingCol, messages_)<<chosenReducedCostVal_<<sigma_<<bestSigma<<CoinMessageEol<<CoinMessageEol;
        return -1;
    }
}


int
CglLandPSimplex::findBestPivotColumn(int direction,
                                     double pivotTol, bool reducedSpace, bool allowDegenerate, bool modularize)
{
    TabRow newRow(this);
    newRow.reserve(ncols_ + nrows_);
    int varOut=-1;

    adjustTableauRow(basics_[row_i_.num], row_i_, direction);

    double m = si_->getInfinity();

    int j = 0;
    double gamma = 0.;

    for (; j< ncols_orig_ ; j++)
    {
        if (reducedSpace &&
                !colCandidateToLeave_[j]
           )
            continue;
        if (fabs(row_i_[nonBasics_[j]])< pivotTol)
        {
            continue;
        }
        gamma = - row_k_[nonBasics_[j]]/row_i_[nonBasics_[j]];


        newRow[basics_[row_k_.num]] = 1.;
        newRow.rhs = row_k_.rhs + gamma * row_i_.rhs;
        if (newRow.rhs > 1e-5  && newRow.rhs < 1 - 1e-5 )
        {
            double m_j = computeCglpObjective(gamma, modularize, newRow);
            if (m_j < m)
            {
                varOut = j;
                m = m_j;
            }
        }
    }
    resetOriginalTableauRow(basics_[row_i_.num], row_i_, direction);

    if (m < sigma_ )
    {
        handler_->message(FoundBestImprovingCol, messages_)<<nonBasics_[varOut]<<gamma<<m<<CoinMessageEol<<CoinMessageEol;
        inDegenerateSequence_ = false;
         assert (varOut<0||direction!=0);
        return varOut;
    }
    else if (allowDegenerate && m<=sigma_)
    {
        inDegenerateSequence_ = true;
    }
    else
    {
        return -1;
    }
    return -1;
}

struct reducedCost
{
    /** To avoid computing two times the same row direction will have strange
     values, direction is -1 or 1 if for only one of the two direction
     rc is <0 and is -2 or 2 if the two direction have one <0 rc with the sign
     indicating which one of the two directions is the best.<br>
     Note that by theory only one reduced cost (for u_i, or v_i)
     maybe negative for each direction.
     */
    int direction;
    /** gammSign is the sign of gamma (corresponding to taking rc for u_i or v_i)
     for the best of the two rc for this row.*/
    int gammaSign;
    /** gammaSign2 is the sign of gamma for the worst of the two rc for this row.*/
    int gammaSign2;
    /** if both reduced costs are <0 value is the smallest of the two.*/
    double value;
    /** greatest of the two reduced costs */
    double value2;
    /** index of the row.*/
    int row;
    bool operator<(const reducedCost & other) const
    {
        return (value>other.value);
    }
};

/** Find incoming and leaving variables which lead to the most violated
 adjacent normalized lift-and-project cut.
 \remark At this point reduced costs should be already computed.
 \return incoming variable variable,
 \param leaving variable
 \param direction leaving direction
 \param numTryRows number rows tried
 \param pivotTol pivot tolerance
 \param reducedSpace separaration space (reduced or full)
 \param allowNonStrictlyImproving wether or not to allow non stricly improving pivots.
 */
int CglLandPSimplex::findBestPivot(int &leaving, int & direction,
                                   const CglLandP::Parameters & params)
{
    // 1. Sort <0 reduced costs in increasing order
    // 2. for numTryRows reduced costs call findBestPivotColumn
    // if better record

    // The reduced cost are already here in rWk1_ is u^l_j
    // rwk2_ is u^u_j rWk3_ is v^u_j rWk4_ is v^u_j
    // Let's rename not to get too much confused
    double * ul_i = &rWk1_[0];
    double * uu_i = &rWk2_[0];
    double * vl_i = &rWk3_[0];
    double * vu_i = &rWk4_[0];

    reducedCost * rc = new reducedCost[nNegativeRcRows_];
    int k = 0;
    rc[k].direction = 0;//initialize first rc
    //int k2 = 0;
    for (int i = 0 ; i < nrows_ ; i++)
    {
        if (ul_i[i] < -params.pivotTol)
            //       && rowFlags_[i]) //row has not been flaged
        {
            rc[k].direction = -1;
            rc[k].gammaSign = -1;
            rc[k].value = ul_i[i];
            rc[k].row = i;
            //k2++;
        }

        if (vl_i[i] < -params.pivotTol)
            //&& rowFlags_[i]) //row has not been flaged
        {
            rc[k].direction = -1;
            rc[k].gammaSign = 1;
            rc[k].value = vl_i[i];
            rc[k].row = i;
            //k2++;
        }


        if (uu_i[i] < -params.pivotTol)
            //       && rowFlags_[i]) //row has not been flaged
        {
            if (rc[k].direction == 0)
            {
                rc[k].direction = 1;
                rc[k].gammaSign = -1;
                rc[k].value = uu_i[i];
                rc[k].row = i;
            }
            else
            {
                if (uu_i[i] <rc[k].value)   //this one is better
                {
                    rc[k].direction = 2;
                    rc[k].gammaSign2 = rc[k].gammaSign;
                    rc[k].gammaSign = -1;
                    rc[k].value2 = rc[k].value;
                    rc[k].value = uu_i[i];
                }
                else
                {
                    rc[k].direction = -2;
                    rc[k].gammaSign2 = -1;
                    rc[k].value2 = uu_i[i];
                }
            }
            //k2++;
        }
        if (vu_i[i] < -params.pivotTol)
            //&& rowFlags_[i]) //row has not been flaged
        {
            if (rc[k].direction==0)
            {
                rc[k].direction = 1;
                rc[k].gammaSign = 1;
                rc[k].value = vu_i[i];
                rc[k].row = i;
            }
            else
            {
                if (vu_i[i] < rc[k].value)   //this one is better
                {
                    rc[k].direction = 2;
                    rc[k].gammaSign2 = rc[k].gammaSign;
                    rc[k].gammaSign = 1;
                    rc[k].value2 = rc[k].value;
                    rc[k].value = vu_i[i];
                }
                else   //the other one is better
                {
                    rc[k].direction = -2;
                    rc[k].gammaSign2 = 1;
                    rc[k].value2 = vu_i[i];
                }
            }
            //k2++;
        }
        if (rc[k].direction!=0) //We have added a row with < 0 rc during
            //last iteration
        {
            k++;
            if (k<nNegativeRcRows_)
                rc[k].direction = 0;
            else
                break;
        }

    }

    assert(k==nNegativeRcRows_);

    //now make a heap
    std::make_heap(rc, rc + k);
    //  assert(rc[0].value==chosenReducedCostVal_);
    int bestLeaving = -1;
    int bestIncoming = -1;
    int bestDirection = 0;

    double bestSigma = COIN_DBL_MAX;
    double bestRc = COIN_DBL_MAX;
    //now scan the heap
#ifndef NDEBUG
    int best_l = 0;
#endif
    int notImproved = 0;
    for (int l = 0; l < k && l < 10  ; l++, notImproved++)
    {
        if (!rowFlags_[rc[l].row]) continue;//this row has been marked to be skipped
        //     if(bestLeaving != -1 && rc[l].value > -1e-02) break;
        if (rc[l].value > -1e-02) break;
        row_i_.num=rc[l].row;
        pullTableauRow(row_i_);//Get the tableau row

        //compute f+ or f- for the best negative rc corresponding to this row
        chosenReducedCostVal_ = rc[l].value;
        double sigma;
        int incoming =
            fastFindBestPivotColumn
            (rc[l].direction, rc[l].gammaSign,
             params.pivotTol, params.away,
             (params.sepSpace==CglLandP::Fractional),
             0,
             sigma, params.modularize);
        if (incoming!=-1 && bestSigma > sigma)
        {
            //          std::cout<<"I found a better pivot "<<sigma - sigma_<< " for indice number "<<l<<std::endl;
#ifndef NDEBUG
            best_l = l;
#endif
            bestSigma = sigma;
            bestIncoming = incoming;
            bestLeaving = rc[l].row;
            bestDirection = rc[l].direction > 0 ? 1 : -1;
            bestRc = rc[l].value;
            notImproved = 0;
        }

        //Now evenutally compute f+ or f- for the other negative rc (if if exists)
        if (rc[l].direction == 2 || rc[l].direction == -2)
        {
            rc[l].direction/= -2;//Reverse the direction
            chosenReducedCostVal_ = rc[l].value2;//need to set this for debug double
            //checks
            incoming = fastFindBestPivotColumn
                       (rc[l].direction, rc[l].gammaSign2,
                        params.pivotTol, params.away,
                        (params.sepSpace==CglLandP::Fractional),
                        0,
                        sigma, params.modularize);
            if (incoming!=-1 && bestSigma > sigma)
            {
                // std::cout<<"I found a better pivot "<<sigma - sigma_<<std::endl;
#ifndef NDEBUG
                best_l = l;
#endif
                bestSigma = sigma;
                bestIncoming = incoming;
                bestLeaving = rc[l].row;
                bestDirection = rc[l].direction;
                bestRc = rc[l].value2;
                notImproved = 0;
            }
        }
    }

    row_i_.num = leaving = bestLeaving;
    chosenReducedCostVal_ = bestRc;
    assert(best_l <= nNegativeRcRows_);
    if (bestLeaving!=-1)
    {
        //      std::cout<<"Best pivot pivot "<<best_l<<std::endl;
        pullTableauRow(row_i_);
    }
    direction = bestDirection;
    delete [] rc;
    assert (bestIncoming<0||direction!=0);
    return bestIncoming;

}

double
CglLandPSimplex::computeCglpObjective(const TabRow &row, bool modularize) const
{
    double numerator = -row.rhs * (1 - row.rhs);
    double denominator = 1;

    const int & n = row.getNumElements();
    const int * ind = row.getIndices();
    const double * val = row.denseVector();
    for (int j = 0 ; j < n ; j++)
    {
        const int& jj = ind[j];
        if (col_in_subspace[jj]==false) continue;
        double coeff = val[jj];
        if (modularize && isInteger(jj))
            coeff = modularizedCoef(coeff, row.rhs);
        denominator += normedCoef(fabs(coeff), jj);
        numerator += (coeff > 0 ?
                      coeff *(1- row.rhs):
                      - coeff * row.rhs)*getColsolToCut(jj);
    }
    return numerator*rhs_weight_/denominator;
}

double
CglLandPSimplex::computeCglpObjective(double gamma, bool strengthen)
{
    double rhs = row_k_.rhs + gamma * row_i_.rhs;
    double numerator = - rhs * (1 - rhs);
    double denominator = 1;

    double coeff = gamma;//newRowCoefficient(basics_[row_i_.num], gamma);
    if (strengthen && isInteger(basics_[row_i_.num]))
        coeff = modularizedCoef(coeff, rhs);
    denominator += normedCoef(fabs(coeff), basics_[row_i_.num]);
    numerator += (coeff > 0 ?
                  coeff *(1- rhs):
                  - coeff * rhs)*
                 getColsolToCut(basics_[row_i_.num]);
    for (int j = 0 ; j < ncols_ ; j++)
    {
        if (col_in_subspace[nonBasics_[j]]==false) {
          continue;
        }
        coeff = newRowCoefficient(nonBasics_[j], gamma);
        if (strengthen && nonBasics_[j] < ncols_orig_ &&  isInteger(j))
            coeff = modularizedCoef(coeff, rhs);
        denominator += normedCoef(fabs(coeff), nonBasics_[j]);
        numerator += (coeff > 0 ?
                      coeff *(1- rhs):
                      - coeff * rhs)*
                     getColsolToCut(nonBasics_[j]);
    }
    return numerator*rhs_weight_/denominator;
}


double
CglLandPSimplex::computeCglpObjective(double gamma, bool strengthen, TabRow & newRow)
{
    newRow.clear();
    newRow.rhs = row_k_.rhs + gamma * row_i_.rhs;
    double numerator = -newRow.rhs * (1 - newRow.rhs);
    double denominator = 1;

    int * indices = newRow.getIndices();
    int k = 0;
    {
        if (col_in_subspace[basics_[row_i_.num]]==false)
          {
            DblEqAssert(0.,1.);
          }
        double & val = newRow[ basics_[row_i_.num]] = gamma;//newRowCoefficient(basics_[row_i_.num], gamma);
        indices[k++] = basics_[row_i_.num];
        if (strengthen && row_i_.num < ncols_orig_ && isInteger(row_i_.num))
            newRow[ basics_[row_i_.num]] = modularizedCoef(newRow[ basics_[row_i_.num]], newRow.rhs);
        denominator += normedCoef(fabs(val), basics_[row_i_.num]);
        numerator += (val > 0 ?
                      val *(1- newRow.rhs):
                      - val * newRow.rhs)*
                     getColsolToCut(basics_[row_i_.num]);
    }
    for (int j = 0 ; j < ncols_ ; j++)
    {
        double & val = newRow[nonBasics_[j]] = newRowCoefficient(nonBasics_[j], gamma);
        indices[k++] = nonBasics_[j];
        if (strengthen && nonBasics_[j] < ncols_orig_ &&  isInteger(j))
            newRow[ nonBasics_[j]] = modularizedCoef(val, newRow.rhs);
        if (col_in_subspace[nonBasics_[j]]==false) continue;
        denominator += normedCoef(fabs(val), nonBasics_[j]);
        numerator += (val > 0 ?
                      val *(1- newRow.rhs):
                      - val * newRow.rhs)*
                     getColsolToCut(nonBasics_[j]);
    }
    newRow.setNumElements(k);
    //assert (fabs(numerator/denominator - computeCglpObjective(newRow))<1e-04);
    return numerator*rhs_weight_/denominator;
}



/** Compute the reduced cost of Cglp */
double
CglLandPSimplex::computeCglpRedCost(int direction, int gammaSign, double tau)
{
    double toBound;
    toBound = direction == -1 ? getLoBound(basics_[row_i_.num])  : getUpBound(basics_[row_i_.num]);

    double value =0;
    int sign = gammaSign * direction;
    double tau1 = 0;
    double tau2 = 0;
    for (unsigned int i = 0 ; i < M3_.size() ; i++)
    {
        tau1 += fabs(row_i_ [M3_[i]]);
        if (sign == 1 && row_i_[M3_[i]] < 0)
        {
            tau2 += row_i_[M3_[i]] * getColsolToCut(M3_[i]);
        }
        else if (sign == -1 && row_i_[M3_[i]] > 0)
        {
            tau2 += row_i_[M3_[i]] * getColsolToCut(M3_[i]);
        }
    }
    double Tau = - sign * (tau + tau2) - tau1 * sigma_;
    value = - sigma_ + Tau
            + (1 - getColsolToCut(basics_[row_k_.num])) * sign * (row_i_.rhs
                    -  toBound)
            + (gammaSign == 1)*direction*(toBound - getColsolToCut(basics_[row_i_.num]));

    return value;
}


/** Compute the value of sigma and thau (which are constants for a row i as defined in Mike Perregaard thesis */
double
CglLandPSimplex::computeRedCostConstantsInRow()
{
    double tau1 = 0; //the part which will be multiplied by sigma
    double tau2 = 0;//the rest

    for (unsigned int i = 0 ; i < M1_.size() ; i++)
    {
        tau1 += row_i_[M1_[i]];
    }
    for (unsigned int i = 0 ; i < M2_.size() ; i++)
    {
        tau1 -= row_i_[M2_[i]];
        tau2 += row_i_[M2_[i]] * getColsolToCut(M2_[i]);
    }
    return sigma_ * tau1 + tau2;

}

void
CglLandPSimplex::updateM1_M2_M3(TabRow & row, double tolerance, bool perturb)
{
    M1_.clear();
    M2_.clear();
    M3_.clear();
    tolerance = 0;
    for (int i = 0; i<ncols_ ; i++)
    {
        const int &ii = nonBasics_[i];
        const double &f = row[ii];
        if (f< -tolerance)
        {
            if (col_in_subspace[ii])
            {
                M1_.push_back(ii);
                colCandidateToLeave_[i]=true;
            }
            else
            {
                colCandidateToLeave_[i]=false;
            }
        }
        else if (f>tolerance)
        {
            if (col_in_subspace[ii])
            {
                M2_.push_back(ii);
                colCandidateToLeave_[i]=true;
            }
            else
            {
                colCandidateToLeave_[i]=false;
            }
        }
        else
        {
            if (col_in_subspace[ii])
            {
                if (perturb)   //assign to M1 or M2 at random
                {
                    int sign = CoinDrand48() > 0.5 ? 1 : -1;
                    if (sign == -1)   //put into M1
                    {
                        M1_.push_back(ii);
                        colCandidateToLeave_[i]=true;
                    }
                    else   //put into M2
                    {
                        M2_.push_back(ii);
                        colCandidateToLeave_[i]=true;
                    }
                }
                else
                {
                    M3_.push_back(ii);
                    colCandidateToLeave_[i] = true;
                }
            }
            else
            {
                colCandidateToLeave_[i] = false;
            }
        }
    }
    //std::cout<<"M3 has "<<M3_.size()<<" variables."<<std::endl;
}

/** Create the intersection cut of row k*/
void
CglLandPSimplex::createIntersectionCut(TabRow & row, OsiRowCut &cut) const
{
    const double * colLower = si_->getColLower();
    const double * rowLower = si_->getRowLower();
    const double * colUpper = si_->getColUpper();
    const double * rowUpper = si_->getRowUpper();
    // double f_0 = row.rhs;
    //put the row back into original form
    for (int j = 0; j < ncols_ ; j++)
    {
        if ((nonBasics_[j] < ncols_))
        {
            CoinWarmStartBasis::Status status = getStatus(nonBasics_[j]);

            if (status==CoinWarmStartBasis::atLowerBound)
            {
                //        row.rhs += getLoBound(nonBasics_[j]) * row[nonBasics_[j]];
            }
            else if (status==CoinWarmStartBasis::atUpperBound)
            {
                row[nonBasics_[j]] = - row[nonBasics_[j]];
                //        row.rhs += getUpBound(nonBasics_[j]) * row[nonBasics_[j]];
            }
            else
            {
                throw;
            }
        }
    }



    //  return ;

    cut.setUb(COIN_DBL_MAX);
    double * vec = new double[ncols_orig_+ nrows_orig_ ];
    CoinFillN(vec, ncols_orig_ + nrows_orig_, 0.);
    double infty = si_->getInfinity();
    double cutRhs = row.rhs;
    cutRhs = cutRhs * (1 - cutRhs);
    for (int j = 0; j < ncols_ ; j++)
    {
        if (fabs(row[nonBasics_[j]])>1e-10)
        {
            double value = intersectionCutCoef(row[nonBasics_[j]], row.rhs);

            if (nonBasics_[j]<ncols_)
            {
                CoinWarmStartBasis::Status status = //CoinWarmStartBasis::basic;
                    basis_->getStructStatus(nonBasics_[j]);
                if (status==CoinWarmStartBasis::atUpperBound)
                {
                    value = - intersectionCutCoef(- row[nonBasics_[j]], row.rhs) ;
                    cutRhs += value * colUpper[nonBasics_[j]];
                }
                else
                    cutRhs += value * colLower[nonBasics_[j]];
                vec[original_index_[nonBasics_[j]]] += value;
            }
            else if (nonBasics_[j]>=ncols_)
            {
                int iRow = nonBasics_[j] - ncols_;

                if (rowLower[iRow] > -infty)
                {
                    value = -value;
                    cutRhs -= value*rowLower[iRow];
                    assert(basis_->getArtifStatus(iRow)==CoinWarmStartBasis::atUpperBound ||
                           (fabs(rowLower[iRow] - rowUpper[iRow]) < 1e-08));
                }
                else
                {
                    cutRhs -= value*rowUpper[iRow];
                    assert(basis_->getArtifStatus(iRow)==CoinWarmStartBasis::atLowerBound);
                }
                vec[nonBasics_[j]] = value;
                assert(fabs(cutRhs)<1e100);
            }
        }
    }

    const CoinPackedMatrix * mat = si_->getMatrixByCol();
    const CoinBigIndex * starts = mat->getVectorStarts();
    const int * lengths = mat->getVectorLengths();
    const double * values = mat->getElements();
    const int * indices = mat->getIndices();
    for (int j = 0 ; j < ncols_ ; j++)
    {
        const CoinBigIndex& start = starts[j];
        CoinBigIndex end = start + lengths[j];
        for (CoinBigIndex k = start ; k < end ; k++)
        {
            vec[original_index_[j]] -= vec[original_index_[ncols_ + indices[k]]] * values[k];
        }

    }

    //Pack vec into the cut
    int * inds = new int [ncols_orig_];
    int nelem = 0;
    for (int i = 0 ; i < ncols_orig_ ; i++)
    {
        if (fabs(vec[i]) > COIN_INDEXED_TINY_ELEMENT)
        {
            vec[nelem] = vec[i];
            inds[nelem++] = i;
        }
    }

    cut.setLb(cutRhs);
    cut.setRow(nelem, inds, vec, false);
    delete [] vec;

}

/** Compute the normalization factor of the cut.*/
double
CglLandPSimplex::normalizationFactor(const TabRow & row) const
{
    double numerator = rhs_weight_;
    double denominator = 1.;
    for (int j = 0 ; j < ncols_ ; j++)
    {
        denominator += fabs(normedCoef(row[nonBasics_[j]], nonBasics_[j]));
    }
    return numerator/denominator;
}


/** Create MIG cut from row k*/
void
CglLandPSimplex::createMIG( TabRow &row, OsiRowCut &cut) const
{
    const double * colLower = si_->getColLower();
    const double * rowLower = si_->getRowLower();
    const double * colUpper = si_->getColUpper();
    const double * rowUpper = si_->getRowUpper();

    if (1)
    {
        double f_0 = row.rhs - floor(row.rhs);
        //put the row back into original form
        for (int j = 0; j < ncols_ ; j++)
        {
            if (nonBasics_[j] < ncols_)
            {
                const CoinWarmStartBasis::Status status = basis_->getStructStatus(nonBasics_[j]);

                if (status==CoinWarmStartBasis::atLowerBound)
                {
                    //        row.rhs += getLoBound(nonBasics_[j]) * row[nonBasics_[j]];
                }
                else if (status==CoinWarmStartBasis::atUpperBound)
                {
                    row[nonBasics_[j]] = - row[nonBasics_[j]];
                    //        row.rhs += getUpBound(nonBasics_[j]) * row[nonBasics_[j]];
                }
                else
                {
                    throw;
                }
            }
        }
        //double scaleFactor = normalizationFactor(row);
        row.rhs = f_0;

        cut.setUb(COIN_DBL_MAX);
        double * vec = new double[ncols_orig_ + nrows_orig_];
        CoinFillN(vec, ncols_orig_ + nrows_orig_, 0.);
        //f_0 = row.rhs - floor(row.rhs);
        double infty = si_->getInfinity();
        double cutRhs = row.rhs - floor(row.rhs);
        cutRhs = cutRhs * (1 - cutRhs);
        assert(fabs(cutRhs)<1e100);
        for (int j = 0; j < ncols_ ; j++)
        {
            if (fabs(row[nonBasics_[j]]))
            {
                {
                    if (nonBasics_[j]<ncols_orig_)
                    {
                        const CoinWarmStartBasis::Status status = basis_->getStructStatus(nonBasics_[j]);
                        double value;
                        if (status==CoinWarmStartBasis::atUpperBound)
                        {
                            value = - strengthenedIntersectionCutCoef(nonBasics_[j], - row[nonBasics_[j]], row.rhs) ;
                            cutRhs += value * colUpper[nonBasics_[j]];
                        }
                        else if (status==CoinWarmStartBasis::atLowerBound)
                        {
                            value = strengthenedIntersectionCutCoef(nonBasics_[j], row[nonBasics_[j]], row.rhs);
                            cutRhs += value * colLower[nonBasics_[j]];
                        }
                        else
                        {
                            std::cerr<<"Invalid basis"<<std::endl;
                            throw -1;
                        }

                        assert(fabs(cutRhs)<1e100);
                        vec[original_index_[nonBasics_[j]]] = value;
                    }
                    else
                    {
                        int iRow = nonBasics_[j] - ncols_;
                        double value = strengthenedIntersectionCutCoef(nonBasics_[j], row[nonBasics_[j]], row.rhs);
			//if (rowLower[iRow]<rowUpper[iRow]) 
			  //printf("row %d status %d lower %g sol %g upper %g\n",
			  //	 iRow,basis_->getArtifStatus(iRow),
			  //	 rowLower[iRow],si_->getRowActivity()[iRow],
			  //	 rowUpper[iRow]);
                        if (rowUpper[iRow] < infty)
                        {
			  if (basis_->getArtifStatus(iRow)==CoinWarmStartBasis::atUpperBound)
                            cutRhs -= value*rowLower[iRow];
			  else
                            cutRhs -= value*rowUpper[iRow];
                        }
                        else
                        {
                            value = -value;
                            cutRhs -= value*rowLower[iRow];
                            //assert(basis_->getArtifStatus(iRow)==CoinWarmStartBasis::atUpperBound ||
			    //     (rowUpper[iRow] < infty));
                        }
                        vec[original_index_[nonBasics_[j]]] = value;
                        assert(fabs(cutRhs)<1e100);
                    }
                }
            }
        }

        //Eliminate slacks
        eliminate_slacks(vec);

        //Pack vec into the cut
        int * inds = new int [ncols_orig_];
        int nelem = 0;
        for (int i = 0 ; i < ncols_orig_ ; i++)
        {
            if (fabs(vec[i]) > COIN_INDEXED_TINY_ELEMENT)
            {
                vec[nelem] = vec[i];
                inds[nelem++] = i;
            }
        }

        cut.setLb(cutRhs);
        cut.setRow(nelem, inds, vec, false);
        //std::cout<<"Scale factor: "<<scaleFactor<<" rhs weight "<<rhs_weight_<<std::endl;
        //scaleCut(cut, scaleFactor);
        delete [] vec;
        delete [] inds;

    }
}

void
CglLandPSimplex::eliminate_slacks(double * vec) const
{
    const CoinPackedMatrix * mat = si_->getMatrixByCol();
    const CoinBigIndex * starts = mat->getVectorStarts();
    const int * lengths = mat->getVectorLengths();
    const double * values = mat->getElements();
    const int * indices = mat->getIndices();
    const double * vecSlacks = vec + ncols_orig_;
    for (int j = 0 ; j < ncols_ ; j++)
    {
        const CoinBigIndex& start = starts[j];
        CoinBigIndex end = start + lengths[j];
        double & val = vec[original_index_[j]];
        for (CoinBigIndex k = start ; k < end ; k++)
        {
            val -= vecSlacks[indices[k]] * values[k];
        }
    }
}
void
CglLandPSimplex::scaleCut(OsiRowCut & cut, double factor) const
{
    DblGtAssert(factor, 0.);
    cut *= factor;
    cut.setLb(cut.lb()*factor);
}

#ifdef APPEND_ROW
void
CglLandPSimplex::append_row(int row_num, bool modularize)
{
    int old_idx = basics_[row_num];
    int r_idx = nrows_ - 1;

    if (basis_->getArtifStatus(nrows_ - 1) == CoinWarmStartBasis::basic)
    {
        int i = 0;
        for (; i < ncols_ ; i++)
        {
            if (nonBasics_[i] == ncols_ + nrows_ -1)
                break;
        }
        if (i < ncols_)
        {
            nonBasics_[i] = ncols_ + nrows_ - 1;
            assert(basics_[nrows_ - 1] == ncols_ + nrows_ - 1);
            basics_[nrows_ - 1] = ncols_ - 1;
        }
    }
#if 0
    if (basis_->getStructStatus(ncols_ - 1) != CoinWarmStartBasis::basic)
    {
        basis_->setStructStatus(ncols_ - 1, CoinWarmStartBasis::basic);
        assert(basis_->getArtifStatus(nrows_ - 1) == CoinWarmStartBasis::basic);
        basis_->setArtifStatus(nrows_ - 1, CoinWarmStartBasis::atLowerBound);
        si_->pivot(ncols_ - 1, ncols_ + nrows_ - 1, -1);
    }
#endif
    TabRow row(this);
    int rowsize = ncols_orig_ + nrows_orig_ + 1;
    row.reserve(rowsize);
    row.num = row_num;
    pullTableauRow(row);
    row.rhs -= floor(row.rhs);
    if (basis_->getStructStatus(ncols_ - 1) != CoinWarmStartBasis::basic)
    {
        basis_->setStructStatus(ncols_ - 1, CoinWarmStartBasis::basic);
        assert(basis_->getArtifStatus(nrows_ - 1) == CoinWarmStartBasis::basic);
        basis_->setArtifStatus(nrows_ - 1, CoinWarmStartBasis::atLowerBound);
        si_->pivot(ncols_ - 1, ncols_ + nrows_ - 1, -1);
    }

    int n = row.getNumElements();
    const int* ind = row.getIndices();
#ifdef CGL_HAS_OSICLP
    CoinPackedMatrix * m = clp_->getMutableMatrixByCol();
#endif

#if 1
    const double * rowLower = si_->getRowLower();
    const double * rowUpper = si_->getRowUpper();
    double infty = si_->getInfinity();
    double rhs = floor (colsol_[old_idx]);//floor(row.rhs);
    std::vector<double> vec(rowsize,0.);
    for (int i = 0 ; i < n ; i++)
    {
        const int &ni = ind[i];
        if (ni == old_idx) continue;
        if (integers_[ni])
        {
            double f = modularizedCoef(row[ni],row.rhs);
            f = floor(row[ni] - f + 0.5);
            row[ni] -= f;
            if (ni < ncols_)
            {
                if (basis_->getStructStatus(ni) == CoinWarmStartBasis::atUpperBound)
                {
                    f *= -1;
                    rhs += f * getUpBound(ni);
                }
                if (basis_->getStructStatus(ni) == CoinWarmStartBasis::atLowerBound)
                {
                    //f *= -1;
                    rhs -= f * getLoBound(ni);
                }
                vec[ni]=f;
            }
            else
            {
                int iRow = ind[i] - ncols_;
                if (basis_->getArtifStatus(iRow)==CoinWarmStartBasis::atLowerBound)
                {
                    rhs -= f*rowUpper[iRow];
                }
                else if (basis_->getArtifStatus(iRow)==CoinWarmStartBasis::atUpperBound)
                {
                    f *= -1;
                    rhs -= f*rowLower[iRow];
                    assert(rowLower[iRow] < infty);
                }
                vec[ni] = f;
            }
        }
    }
    vec[old_idx] = 1;

    eliminate_slacks(&vec[0]);

    for (int i = 0 ; i < rowsize ; i++)
    {
        if (vec[i])
            m->modifyCoefficient(r_idx, i, vec[i], true);
    }
    si_->setRowBounds(r_idx, rhs, rhs);
#endif

#if 1
    si_->messageHandler()->setLogLevel(0);
    si_->resolve();
    assert(si_->getIterationCount() == 0);
#endif
    /* Update the cached info.*/
#ifndef NDEBUG
    int * basics = new int[nrows_];
#else
    int * basics = basics_;
#endif
    si_->getBasics(basics);
#ifndef NDEBUG
    for (int i = 0 ; i < r_idx ; i++)
    {
        assert(basics[i] == basics_[i]);
    }
    assert(basics[r_idx] == ncols_ - 1);
    delete [] basics_;
    basics_ = basics;
#endif
    assert(basis_->getStructStatus(old_idx) == CoinWarmStartBasis::basic);
    assert(basis_->getStructStatus(ncols_ - 1) == CoinWarmStartBasis::basic);

    n = 0;
    for (int i = 0 ; i < ncols_ ; i++)
    {
        if (basis_->getStructStatus(i) != CoinWarmStartBasis::basic)
        {
            nonBasics_[n++] = i;
        }
    }
    for (int i = 0 ; i < nrows_ ; i++)
    {
        if (basis_->getArtifStatus(i) != CoinWarmStartBasis::basic)
        {
            nonBasics_[n++] = i + ncols_;
        }
    }
    assert (n == ncols_);
    colsol_[ncols_ - 1] = si_->getColSolution()[ncols_-1];
    colsolToCut_[ncols_ - 1] = si_->getColSolution()[ncols_-1];
}

void
CglLandPSimplex::check_mod_row(TabRow & row)
{
    int rowsize = ncols_orig_ + nrows_orig_;
    for (int i = 0 ; i < rowsize ; i++)
    {
        assert(! integers_[i] || ((row[i] <= row.rhs+1e-08) && (row[i] >= row.rhs - 1- 1e-08)));
    }
}

void
CglLandPSimplex::update_row(TabRow &row)
{
    int r_idx = nrows_ - 1;
    int c_idx = ncols_ - 1;
    assert(basics_[row.num] == ncols_ - 1);
    assert(r_idx == row.num);
    int rowsize = nrows_ + ncols_;
    int n = row.getNumElements();
    const int* ind = row.getIndices();
#ifdef CGL_HAS_OSICLP
    CoinPackedMatrix * m = clp_->getMutableMatrixByCol();
#endif
    //double * m_el = m->getMutableElements();
    //const int * m_id = m->getIndices();
    //const int * m_le = m->getVectorLengths();
    //const CoinBigIndex * m_st = m->getVectorStarts();
    const double * rowLower = si_->getRowLower();
    const double * rowUpper = si_->getRowUpper();
    double infty = si_->getInfinity();
    double rhs = floor(row.rhs);
    std::vector<double> vec(rowsize,0.);
    for (int i = 0 ; i < n ; i++)
    {
        const int &ni = ind[i];
        if (ni == c_idx) continue;
        if (integers_[ni])
        {
            double f = modularizedCoef(row[ni],row.rhs);
            f = floor(row[ni] - f + 0.5);
            //row[ni] -= f;
            if (ni < ncols_)
            {
                assert(basis_->getStructStatus(ni) != CoinWarmStartBasis::basic);
                if (basis_->getStructStatus(ni) == CoinWarmStartBasis::atUpperBound)
                {
                    f *= -1;
                    rhs += f * getUpBound(ni);
                }
                if (basis_->getStructStatus(ni) == CoinWarmStartBasis::atLowerBound)
                {
                    rhs -= f * getLoBound(ni);
                }
                vec[ni]=f;
            }
            else
            {
                //continue;
                int iRow = ind[i] - ncols_;
                if (iRow == r_idx) continue;
                assert(basis_->getArtifStatus(iRow) != CoinWarmStartBasis::basic);
                if (basis_->getArtifStatus(iRow)==CoinWarmStartBasis::atLowerBound)
                {
                    rhs -= f*rowUpper[iRow];
                }
                else if (basis_->getArtifStatus(iRow)==CoinWarmStartBasis::atUpperBound)
                {
                    f *= -1;
                    rhs -= f*rowLower[iRow];
                    assert(rowUpper[iRow] < infty);
                }
                vec[ni] = f;
            }
        }
    }
    std::copy(vec.begin(),vec.end(), std::ostream_iterator<double>(std::cout, "\t"));
    std::cout<<std::endl;

    //vec[c_idx] = 1;
    eliminate_slacks(&vec[0]);
    std::copy(vec.begin(),vec.end(), std::ostream_iterator<double>(std::cout, "\t"));
    std::cout<<std::endl;
//    assert(vec[c_idx] == 0);
    for (int i = 0 ; i < ncols_ ; i++)
    {
        if (vec[i])
        {
            double f = m->getCoefficient(r_idx,i);
            m->modifyCoefficient(r_idx, i, f + vec[i], true);
        }
    }
    si_->setRowBounds(r_idx, rowLower[r_idx] + rhs, rowUpper[r_idx] + rhs);

#ifdef CGL_HAS_OSICLP
    clp_->getModelPtr()->factorize();
#endif

    pullTableauRow(row);
    assert(row.rhs >= 0. && row.rhs <= 1.);
    bool stop = false;
    for (int i = 0 ; i < rowsize ; i++)
    {
        if (integers_[i] && (row[i] > row.rhs || row[i] < row.rhs - 1))
        {
            stop = true;
        }
        //assert(! integers_[i] || (row[i] <= row.rhs && row[i] > row.rhs - 1));
    }
    if (stop) exit(10);

}
#endif

/** Get the row i of the tableau */
void
CglLandPSimplex::pullTableauRow(TabRow &row) const
{
    const double * rowLower = si_->getRowLower();
    const double * rowUpper = si_->getRowUpper();

    row.clear();
    row.modularized_ = false;
    double infty = si_->getInfinity();
    /* Get the row */
#ifdef CGL_HAS_OSICLP
    if (clp_)
    {
        CoinIndexedVector array2;
        array2.borrowVector(nrows_, 0, row.getIndices() + ncols_, row.denseVector() + ncols_);
        clp_->getBInvARow(row.num, &row, &array2);
        {
            int n = array2.getNumElements();
            int * indices1 = row.getIndices() + row.getNumElements();
            int * indices2 = array2.getIndices();
            for ( int i = 0 ; i < n ; i++)
            {
                *indices1 = indices2[i] + ncols_;
                indices1++;
            }
            row.setNumElements(n + row.getNumElements());
            array2.returnVector();
        }
    }
    else
#endif
    {
        si_->getBInvARow(row.num,row.denseVector(),row.denseVector() + ncols_);
    }
    //Clear basic element (it is a trouble for most of the computations)
    row[basics_[row.num]]=0.;
    //  row.row[basics_[row.num]]=1;
    /* get the rhs */
    {
        int iCol = basics_[row.num];
        if (iCol<ncols_)
            row.rhs = si_->getColSolution()[iCol];
        else   // Osi does not give direct acces to the value of slacks
        {
            iCol -= ncols_;
            row.rhs = - si_->getRowActivity()[iCol];
            if (rowLower[iCol]> -infty)
            {
                row.rhs += rowLower[iCol];
            }
            else
            {
                row.rhs+= rowUpper[iCol];
            }
        }
    }
    //Now adjust the row of the tableau to reflect non-basic variables activity
    for (int j = 0; j < ncols_ ; j++)
    {
        if (nonBasics_[j]<ncols_)
        {
            if (basis_->getStructStatus(nonBasics_[j])==CoinWarmStartBasis::atLowerBound)
            {
            }
            else if (basis_->getStructStatus(nonBasics_[j])==CoinWarmStartBasis::atUpperBound)
            {
                row[nonBasics_[j]] = -row[nonBasics_[j]];
            }
            else
            {
                std::cout<<(basis_->getStructStatus(nonBasics_[j])==CoinWarmStartBasis::isFree)<<std::endl;
                throw CoinError("Invalid basis","CglLandPSimplex","pullTableauRow");
            }
        }
        else
        {
            int iRow = nonBasics_[j] - ncols_;

            if (basis_->getArtifStatus(iRow)==CoinWarmStartBasis::atUpperBound)
            {
                row[nonBasics_[j]] = -row[nonBasics_[j]];
            }
        }
    }
    //  row.clean(1e-30);
}

/** Adjust the row of the tableau to reflect leaving variable direction */
void
CglLandPSimplex::adjustTableauRow(int var, TabRow & row, int direction)
{
    double bound = 0;
    assert(direction != 0);
    if (direction > 0)
    {
        for (int j = 0 ; j < ncols_orig_ ; j++)
            row[nonBasics_[j]] = - row[nonBasics_[j]];

        row.rhs = -row.rhs;
        bound = getUpBound(var);
        setColsolToCut(var, bound - getColsolToCut(var));
        row.rhs += bound;
    }
    else if (direction < 0)
    {
        bound = getLoBound(var);
        setColsolToCut(var, getColsolToCut(var) - bound);
        row.rhs -= bound;
    }
    //  assert(fabs(row.rhs)<1e100);
}


/** reset the tableau row after a call to adjustTableauRow */
void
CglLandPSimplex::resetOriginalTableauRow(int var, TabRow & row, int direction)
{
    if (direction > 0)
    {
        adjustTableauRow(var, row, direction);
    }
    else
    {
        double bound = getLoBound(var);
        row.rhs += bound;
        setColsolToCut(var, getColsolToCut(var) + bound);
    }
}

template <class X>
struct SortingOfArray
{
    X * array;
    SortingOfArray(X* a): array(a) {}

    bool operator()(const int i, const int j)
    {
        return array[i] < array[j];
    }
};

void
CglLandPSimplex::removeRows(int nDelete, const int * rowsIdx)
{
    std::vector<int> sortedIdx;
    for (int i = 0 ; i < nDelete ; i++)
        sortedIdx.push_back(rowsIdx[i]);
    std::sort(sortedIdx.end(), sortedIdx.end());
    si_->deleteRows(nDelete, rowsIdx);
    int k = 1;
    int l = sortedIdx[0];
    for (int i = sortedIdx[0] + 1 ; k < nDelete ; i++)
    {
        if (sortedIdx[k] == i)
        {
            k++;
        }
        else
        {
            original_index_[l] = original_index_[i];
            l++;
        }
    }
    delete basis_;
    basis_ = dynamic_cast<CoinWarmStartBasis *> (si_->getWarmStart());
    assert(basis_);

    /* Update rowFlags_ */
    std::vector<int> order(nrows_);
    for (unsigned int i = 0 ; i < order.size() ; i++)
    {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), SortingOfArray<int>(basics_));
    k = 0;
    l = 0;
    for (int i = 0 ; k < nDelete ; i++)
    {
        if (basics_[order[i]] == sortedIdx[k])
        {
            basics_[order[i]] = -1;
            k++;
        }
        else
        {
            order[l] = order[i];
            l++;
        }
    }
    k = 0;
    for (int i = 0 ; i < nrows_ ; i++)
    {
        if (basics_[i] == -1)
        {
            k++;
        }
        else
        {
            basics_[l] = basics_[i];
            rowFlags_[l] = rowFlags_[i];
            rWk1_[l] = rWk1_[i];
            rWk2_[l] = rWk2_[i];
            rWk4_[l] = rWk3_[i];
            rWk4_[l] = rWk4_[i];
            if (row_k_.num == i) row_k_.num = l;

            l++;
        }
    }

    nrows_ = nrows_ - nDelete;
    original_index_.resize(nrows_);

    int numStructural = basis_->getNumStructural();
    assert(ncols_ = numStructural);
    int nNonBasics = 0;
    for (int i = 0 ; i < numStructural ; i++)
    {
        if (basis_->getStructStatus(i) != CoinWarmStartBasis::basic)
        {
            nonBasics_[nNonBasics++] = i;
        }
    }

    int numArtificial = basis_->getNumArtificial();
    assert(nrows_ = numArtificial);
    for (int i = 0 ; i < numArtificial ; i++)
    {
        if (basis_->getArtifStatus(i)!= CoinWarmStartBasis::basic)
        {
            nonBasics_[nNonBasics++] = i + numStructural;
        }
    }
    assert (nNonBasics == ncols_);
}

  void
  CglLandPSimplex::printEverything(){
    row_k_.print(std::cout, 2, nonBasics_, ncols_);
    printf("nonBasics_: ");
    for(int i = 0 ; i < ncols_ ; i++){
      printf("%5i ",nonBasics_[i]);
    }
    printf("\n");

 printf("basics_: ");
    for(int i = 0 ; i < nrows_ ; i++){
      printf("%5i ",basics_[i]);
    }
    printf("\n");

    printf("source row:");
    for(int i = 0 ; i < ncols_ + nrows_ ; i++){
      printf("%10.9g ", row_k_[i]);
    }
    printf("%10.9g", row_k_.rhs);
    printf("\n");

    printf(" source indices: ");
    for(int i = 0 ; i < row_k_.getNumElements() ; i++){
      printf("%5i %20.20g ", row_k_.getIndices()[i], row_k_[row_k_.getIndices()[i]]);
    }
    printf("\n");

    printf("colsolToCut: ");
    for(int i = 0 ; i < ncols_ + nrows_ ; i++){
      printf("%10.6g ", colsolToCut_[i]);
    }
    printf("\n");

    printf("colsol: ");
    for(int i = 0 ; i < ncols_ + nrows_ ; i++){
      printf("%10.6g ", colsol_[i]);
    }
    printf("\n");
  }
}/* Ends LAP namespace.*/

