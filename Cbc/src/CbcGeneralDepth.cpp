// Edwin 11/10/2009-- carved out of CbcBranchActual
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CoinTypes.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcGeneralDepth.hpp" 
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"
#include "ClpNode.hpp"
#include "CbcBranchDynamic.hpp"
// Default Constructor
CbcGeneralDepth::CbcGeneralDepth ()
        : CbcGeneral(),
        maximumDepth_(0),
        maximumNodes_(0),
        whichSolution_(-1),
        numberNodes_(0),
        nodeInfo_(NULL)
{
}

// Useful constructor (which are integer indices)
CbcGeneralDepth::CbcGeneralDepth (CbcModel * model, int maximumDepth)
        : CbcGeneral(model),
        maximumDepth_(maximumDepth),
        maximumNodes_(0),
        whichSolution_(-1),
        numberNodes_(0),
        nodeInfo_(NULL)
{
    assert(maximumDepth_ < 1000000);
    if (maximumDepth_ > 0)
        maximumNodes_ = (1 << maximumDepth_) + 1 + maximumDepth_;
    else if (maximumDepth_ < 0)
        maximumNodes_ = 1 + 1 - maximumDepth_;
    else
        maximumNodes_ = 0;
#define MAX_NODES 100
    maximumNodes_ = CoinMin(maximumNodes_, 1 + maximumDepth_ + MAX_NODES);
    if (maximumNodes_) {
        nodeInfo_ = new ClpNodeStuff();
        nodeInfo_->maximumNodes_ = maximumNodes_;
        ClpNodeStuff * info = nodeInfo_;
        // for reduced costs and duals
        info->solverOptions_ |= 7;
        if (maximumDepth_ > 0) {
            info->nDepth_ = maximumDepth_;
        } else {
            info->nDepth_ = - maximumDepth_;
            info->solverOptions_ |= 32;
        }
        ClpNode ** nodeInfo = new ClpNode * [maximumNodes_];
        for (int i = 0; i < maximumNodes_; i++)
            nodeInfo[i] = NULL;
        info->nodeInfo_ = nodeInfo;
    } else {
        nodeInfo_ = NULL;
    }
}

// Copy constructor
CbcGeneralDepth::CbcGeneralDepth ( const CbcGeneralDepth & rhs)
        : CbcGeneral(rhs)
{
    maximumDepth_ = rhs.maximumDepth_;
    maximumNodes_ = rhs.maximumNodes_;
    whichSolution_ = -1;
    numberNodes_ = 0;
    if (maximumNodes_) {
        assert (rhs.nodeInfo_);
        nodeInfo_ = new ClpNodeStuff(*rhs.nodeInfo_);
        nodeInfo_->maximumNodes_ = maximumNodes_;
        ClpNodeStuff * info = nodeInfo_;
        if (maximumDepth_ > 0) {
            info->nDepth_ = maximumDepth_;
        } else {
            info->nDepth_ = - maximumDepth_;
            info->solverOptions_ |= 32;
        }
        if (!info->nodeInfo_) {
            ClpNode ** nodeInfo = new ClpNode * [maximumNodes_];
            for (int i = 0; i < maximumNodes_; i++)
                nodeInfo[i] = NULL;
            info->nodeInfo_ = nodeInfo;
        }
    } else {
        nodeInfo_ = NULL;
    }
}

// Clone
CbcObject *
CbcGeneralDepth::clone() const
{
    return new CbcGeneralDepth(*this);
}

// Assignment operator
CbcGeneralDepth &
CbcGeneralDepth::operator=( const CbcGeneralDepth & rhs)
{
    if (this != &rhs) {
        CbcGeneral::operator=(rhs);
        delete nodeInfo_;
        maximumDepth_ = rhs.maximumDepth_;
        maximumNodes_ = rhs.maximumNodes_;
        whichSolution_ = -1;
        numberNodes_ = 0;
        if (maximumDepth_) {
            assert (rhs.nodeInfo_);
            nodeInfo_ = new ClpNodeStuff(*rhs.nodeInfo_);
            nodeInfo_->maximumNodes_ = maximumNodes_;
        } else {
            nodeInfo_ = NULL;
        }
    }
    return *this;
}

// Destructor
CbcGeneralDepth::~CbcGeneralDepth ()
{
    delete nodeInfo_;
}
// Infeasibility - large is 0.5
double
CbcGeneralDepth::infeasibility(const OsiBranchingInformation * /*info*/,
                               int &/*preferredWay*/) const
{
    whichSolution_ = -1;
    // should use genuine OsiBranchingInformation usefulInfo = model_->usefulInformation();
    // for now assume only called when correct
    //if (usefulInfo.depth_>=4&&!model_->parentModel()
    //     &&(usefulInfo.depth_%2)==0) {
    if (true) {
        OsiSolverInterface * solver = model_->solver();
        OsiClpSolverInterface * clpSolver
        = dynamic_cast<OsiClpSolverInterface *> (solver);
        if (clpSolver) {
            ClpNodeStuff * info = nodeInfo_;
            info->integerTolerance_ = model_->getIntegerTolerance();
            info->integerIncrement_ = model_->getCutoffIncrement();
            info->numberBeforeTrust_ = model_->numberBeforeTrust();
            info->stateOfSearch_ = model_->stateOfSearch();
            // Compute "small" change in branch
            int nBranches = model_->getIntParam(CbcModel::CbcNumberBranches);
            if (nBranches) {
                double average = model_->getDblParam(CbcModel::CbcSumChange) / static_cast<double>(nBranches);
                info->smallChange_ =
                    CoinMax(average * 1.0e-5, model_->getDblParam(CbcModel::CbcSmallestChange));
                info->smallChange_ = CoinMax(info->smallChange_, 1.0e-8);
            } else {
                info->smallChange_ = 1.0e-8;
            }
            int numberIntegers = model_->numberIntegers();
            double * down = new double[numberIntegers];
            double * up = new double[numberIntegers];
            int * priority = new int[numberIntegers];
            int * numberDown = new int[numberIntegers];
            int * numberUp = new int[numberIntegers];
            int * numberDownInfeasible = new int[numberIntegers];
            int * numberUpInfeasible = new int[numberIntegers];
            model_->fillPseudoCosts(down, up, priority, numberDown, numberUp,
                                    numberDownInfeasible, numberUpInfeasible);
            info->fillPseudoCosts(down, up, priority, numberDown, numberUp,
                                  numberDownInfeasible,
                                  numberUpInfeasible, numberIntegers);
            info->presolveType_ = 1;
            delete [] down;
            delete [] up;
            delete [] numberDown;
            delete [] numberUp;
            delete [] numberDownInfeasible;
            delete [] numberUpInfeasible;
            bool takeHint;
            OsiHintStrength strength;
            solver->getHintParam(OsiDoReducePrint, takeHint, strength);
            ClpSimplex * simplex = clpSolver->getModelPtr();
            int saveLevel = simplex->logLevel();
            if (strength != OsiHintIgnore && takeHint && saveLevel == 1)
                simplex->setLogLevel(0);
            clpSolver->setBasis();
            whichSolution_ = simplex->fathomMany(info);
            //printf("FAT %d nodes, %d iterations\n",
            //info->numberNodesExplored_,info->numberIterations_);
            //printf("CbcBranch %d rows, %d columns\n",clpSolver->getNumRows(),
            //     clpSolver->getNumCols());
            model_->incrementExtra(info->numberNodesExplored_,
                                   info->numberIterations_);
            // update pseudo costs
            double smallest = 1.0e50;
            double largest = -1.0;
            OsiObject ** objects = model_->objects();
#ifndef NDEBUG
            const int * integerVariable = model_->integerVariable();
#endif
            for (int i = 0; i < numberIntegers; i++) {
#ifndef NDEBUG
                CbcSimpleIntegerDynamicPseudoCost * obj =
                    dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(objects[i]) ;
                assert (obj && obj->columnNumber() == integerVariable[i]);
#else
                CbcSimpleIntegerDynamicPseudoCost * obj =
                    static_cast <CbcSimpleIntegerDynamicPseudoCost *>(objects[i]) ;
#endif
                if (info->numberUp_[i] > 0) {
                    if (info->downPseudo_[i] > largest)
                        largest = info->downPseudo_[i];
                    if (info->downPseudo_[i] < smallest)
                        smallest = info->downPseudo_[i];
                    if (info->upPseudo_[i] > largest)
                        largest = info->upPseudo_[i];
                    if (info->upPseudo_[i] < smallest)
                        smallest = info->upPseudo_[i];
                    obj->updateAfterMini(info->numberDown_[i],
                                         info->numberDownInfeasible_[i],
                                         info->downPseudo_[i],
                                         info->numberUp_[i],
                                         info->numberUpInfeasible_[i],
                                         info->upPseudo_[i]);
                }
            }
            //printf("range of costs %g to %g\n",smallest,largest);
            simplex->setLogLevel(saveLevel);
            numberNodes_ = info->nNodes_;
            int numberDo = numberNodes_;
            if (whichSolution_ >= 0)
                numberDo--;
            if (numberDo > 0) {
                return 0.5;
            } else {
                // no solution
                return COIN_DBL_MAX; // say infeasible
            }
        } else {
            return -1.0;
        }
    } else {
        return -1.0;
    }
}

// This looks at solution and sets bounds to contain solution
void
CbcGeneralDepth::feasibleRegion()
{
    // Other stuff should have done this
}
// Redoes data when sequence numbers change
void
CbcGeneralDepth::redoSequenceEtc(CbcModel * /*model*/,
                                 int /*numberColumns*/,
                                 const int * /*originalColumns*/)
{
}

//#define CHECK_PATH
#ifdef CHECK_PATH
extern const double * debuggerSolution_Z;
extern int numberColumns_Z;
extern int gotGoodNode_Z;
#endif
CbcBranchingObject *
CbcGeneralDepth::createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * /*info*/, int /*way*/)
{
    int numberDo = numberNodes_;
    if (whichSolution_ >= 0)
        numberDo--;
    assert (numberDo > 0);
    // create object
    CbcGeneralBranchingObject * branch = new CbcGeneralBranchingObject(model_);
    // skip solution
    branch->numberSubProblems_ = numberDo;
    // If parentBranch_ back in then will have to be 2*
    branch->numberSubLeft_ = numberDo;
    branch->setNumberBranches(numberDo);
    CbcSubProblem * sub = new CbcSubProblem[numberDo];
    int iProb = 0;
    branch->subProblems_ = sub;
    branch->numberRows_ = model_->solver()->getNumRows();
    int iNode;
    //OsiSolverInterface * solver = model_->solver();
    OsiClpSolverInterface * clpSolver
    = dynamic_cast<OsiClpSolverInterface *> (solver);
    assert (clpSolver);
    ClpSimplex * simplex = clpSolver->getModelPtr();
    int numberColumns = simplex->numberColumns();
    double * lowerBefore = CoinCopyOfArray(simplex->getColLower(),
                                           numberColumns);
    double * upperBefore = CoinCopyOfArray(simplex->getColUpper(),
                                           numberColumns);
    ClpNodeStuff * info = nodeInfo_;
    double * weight = new double[numberNodes_];
    int * whichNode = new int [numberNodes_];
    // Sort
    for (iNode = 0; iNode < numberNodes_; iNode++) {
        if (iNode != whichSolution_) {
            double objectiveValue = info->nodeInfo_[iNode]->objectiveValue();
            double sumInfeasibilities = info->nodeInfo_[iNode]->sumInfeasibilities();
            int numberInfeasibilities = info->nodeInfo_[iNode]->numberInfeasibilities();
            double thisWeight = 0.0;
#if 1
            // just closest
            thisWeight = 1.0e9 * numberInfeasibilities;
            thisWeight += sumInfeasibilities;
            thisWeight += 1.0e-7 * objectiveValue;
            // Try estimate
            thisWeight = info->nodeInfo_[iNode]->estimatedSolution();
#else
            thisWeight = 1.0e-3 * numberInfeasibilities;
            thisWeight += 1.0e-5 * sumInfeasibilities;
            thisWeight += objectiveValue;
#endif
            whichNode[iProb] = iNode;
            weight[iProb++] = thisWeight;
        }
    }
    assert (iProb == numberDo);
    CoinSort_2(weight, weight + numberDo, whichNode);
    for (iProb = 0; iProb < numberDo; iProb++) {
        iNode = whichNode[iProb];
        ClpNode * node = info->nodeInfo_[iNode];
        // move bounds
        node->applyNode(simplex, 3);
        // create subproblem
        sub[iProb] = CbcSubProblem(clpSolver, lowerBefore, upperBefore,
                                   node->statusArray(), node->depth());
        sub[iProb].objectiveValue_ = node->objectiveValue();
        sub[iProb].sumInfeasibilities_ = node->sumInfeasibilities();
        sub[iProb].numberInfeasibilities_ = node->numberInfeasibilities();
#ifdef CHECK_PATH
        if (simplex->numberColumns() == numberColumns_Z) {
            bool onOptimal = true;
            const double * columnLower = simplex->columnLower();
            const double * columnUpper = simplex->columnUpper();
            for (int i = 0; i < numberColumns_Z; i++) {
                if (iNode == gotGoodNode_Z)
                    printf("good %d %d %g %g\n", iNode, i, columnLower[i], columnUpper[i]);
                if (columnUpper[i] < debuggerSolution_Z[i] || columnLower[i] > debuggerSolution_Z[i] && simplex->isInteger(i)) {
                    onOptimal = false;
                    break;
                }
            }
            if (onOptimal) {
                printf("adding to node %x as %d - objs\n", this, iProb);
                for (int j = 0; j <= iProb; j++)
                    printf("%d %g\n", j, sub[j].objectiveValue_);
            }
        }
#endif
    }
    delete [] weight;
    delete [] whichNode;
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    // restore bounds
    for ( int j = 0; j < numberColumns; j++) {
        if (lowerBefore[j] != lower[j])
            solver->setColLower(j, lowerBefore[j]);
        if (upperBefore[j] != upper[j])
            solver->setColUpper(j, upperBefore[j]);
    }
    delete [] upperBefore;
    delete [] lowerBefore;
    return branch;
}
#endif