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
#include "CbcGeneralBranchingObject.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"


#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"


// Default Constructor
CbcGeneralBranchingObject::CbcGeneralBranchingObject()
        : CbcBranchingObject(),
        subProblems_(NULL),
        node_(NULL),
        numberSubProblems_(0),
        numberSubLeft_(0),
        whichNode_(-1),
        numberRows_(0)
{
    //  printf("CbcGeneral %x default constructor\n",this);
}

// Useful constructor
CbcGeneralBranchingObject::CbcGeneralBranchingObject (CbcModel * model)
        : CbcBranchingObject(model, -1, -1, 0.5),
        subProblems_(NULL),
        node_(NULL),
        numberSubProblems_(0),
        numberSubLeft_(0),
        whichNode_(-1),
        numberRows_(0)
{
    //printf("CbcGeneral %x useful constructor\n",this);
}

// Copy constructor
CbcGeneralBranchingObject::CbcGeneralBranchingObject ( const CbcGeneralBranchingObject & rhs)
        : CbcBranchingObject(rhs),
        subProblems_(NULL),
        node_(rhs.node_),
        numberSubProblems_(rhs.numberSubProblems_),
        numberSubLeft_(rhs.numberSubLeft_),
        whichNode_(rhs.whichNode_),
        numberRows_(rhs.numberRows_)
{
    abort();
    if (numberSubProblems_) {
        subProblems_ = new CbcSubProblem[numberSubProblems_];
        for (int i = 0; i < numberSubProblems_; i++)
            subProblems_[i] = rhs.subProblems_[i];
    }
}

// Assignment operator
CbcGeneralBranchingObject &
CbcGeneralBranchingObject::operator=( const CbcGeneralBranchingObject & rhs)
{
    if (this != &rhs) {
        abort();
        CbcBranchingObject::operator=(rhs);
        delete [] subProblems_;
        numberSubProblems_ = rhs.numberSubProblems_;
        numberSubLeft_ = rhs.numberSubLeft_;
        whichNode_ = rhs.whichNode_;
        numberRows_ = rhs.numberRows_;
        if (numberSubProblems_) {
            subProblems_ = new CbcSubProblem[numberSubProblems_];
            for (int i = 0; i < numberSubProblems_; i++)
                subProblems_[i] = rhs.subProblems_[i];
        } else {
            subProblems_ = NULL;
        }
        node_ = rhs.node_;
    }
    return *this;
}
CbcBranchingObject *
CbcGeneralBranchingObject::clone() const
{
    return (new CbcGeneralBranchingObject(*this));
}


// Destructor
CbcGeneralBranchingObject::~CbcGeneralBranchingObject ()
{
    //printf("CbcGeneral %x destructor\n",this);
    delete [] subProblems_;
}
bool doingDoneBranch = false;
double
CbcGeneralBranchingObject::branch()
{
    double cutoff = model_->getCutoff();
    //printf("GenB %x whichNode %d numberLeft %d which %d\n",
    // this,whichNode_,numberBranchesLeft(),branchIndex());
    if (whichNode_ < 0) {
        assert (node_);
        bool applied = false;
        while (numberBranchesLeft()) {
            int which = branchIndex();
            decrementNumberBranchesLeft();
            CbcSubProblem * thisProb = subProblems_ + which;
            if (thisProb->objectiveValue_ < cutoff) {
                //printf("branch %x (sub %x) which now %d\n",this,
                //     subProblems_,which);
                OsiSolverInterface * solver = model_->solver();
                thisProb->apply(solver);
                OsiClpSolverInterface * clpSolver
                = dynamic_cast<OsiClpSolverInterface *> (solver);
                assert (clpSolver);
                // Move status to basis
                clpSolver->setWarmStart(NULL);
                //ClpSimplex * simplex = clpSolver->getModelPtr();
                node_->setObjectiveValue(thisProb->objectiveValue_);
                node_->setSumInfeasibilities(thisProb->sumInfeasibilities_);
                node_->setNumberUnsatisfied(thisProb->numberInfeasibilities_);
                applied = true;
                doingDoneBranch = true;
                break;
            } else if (numberBranchesLeft()) {
                node_->nodeInfo()->branchedOn() ;
            }
        }
        if (!applied) {
            // no good one
            node_->setObjectiveValue(cutoff + 1.0e20);
            node_->setSumInfeasibilities(1.0);
            node_->setNumberUnsatisfied(1);
            assert (whichNode_ < 0);
        }
    } else {
        decrementNumberBranchesLeft();
        CbcSubProblem * thisProb = subProblems_ + whichNode_;
        assert (thisProb->objectiveValue_ < cutoff);
        OsiSolverInterface * solver = model_->solver();
        thisProb->apply(solver);
        //OsiClpSolverInterface * clpSolver
        //= dynamic_cast<OsiClpSolverInterface *> (solver);
        //assert (clpSolver);
        // Move status to basis
        //clpSolver->setWarmStart(NULL);
    }
    return 0.0;
}
/* Double checks in case node can change its mind!
   Can change objective etc */
void
CbcGeneralBranchingObject::checkIsCutoff(double cutoff)
{
    assert (node_);
    int first = branchIndex();
    int last = first + numberBranchesLeft();
    for (int which = first; which < last; which++) {
        CbcSubProblem * thisProb = subProblems_ + which;
        if (thisProb->objectiveValue_ < cutoff) {
            node_->setObjectiveValue(thisProb->objectiveValue_);
            node_->setSumInfeasibilities(thisProb->sumInfeasibilities_);
            node_->setNumberUnsatisfied(thisProb->numberInfeasibilities_);
            break;
        }
    }
}
// Print what would happen
void
CbcGeneralBranchingObject::print()
{
    //printf("CbcGeneralObject has %d subproblems\n",numberSubProblems_);
}
// Fill in current objective etc
void
CbcGeneralBranchingObject::state(double & objectiveValue,
                                 double & sumInfeasibilities,
                                 int & numberUnsatisfied, int which) const
{
    assert (which >= 0 && which < numberSubProblems_);
    const CbcSubProblem * thisProb = subProblems_ + which;
    objectiveValue = thisProb->objectiveValue_;
    sumInfeasibilities = thisProb->sumInfeasibilities_;
    numberUnsatisfied = thisProb->numberInfeasibilities_;
}
/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcGeneralBranchingObject::compareOriginalObject
(const CbcBranchingObject* /*brObj*/) const
{
    throw("must implement");
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcGeneralBranchingObject::compareBranchingObject
(const CbcBranchingObject* /*brObj*/, const bool /*replaceIfOverlap*/)
{
    throw("must implement");
}
#endif