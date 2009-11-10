// Edwin 11/9/2009-- carved out of CbcBranchActual
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
#include "CbcSimpleInteger.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"


//##############################################################################

/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
CbcSimpleInteger::CbcSimpleInteger ()
        : CbcObject(),
        originalLower_(0.0),
        originalUpper_(1.0),
        breakEven_(0.5),
        columnNumber_(-1),
        preferredWay_(0)
{
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
CbcSimpleInteger::CbcSimpleInteger ( CbcModel * model, int iColumn, double breakEven)
        : CbcObject(model)
{
    columnNumber_ = iColumn ;
    originalLower_ = model->solver()->getColLower()[columnNumber_] ;
    originalUpper_ = model->solver()->getColUpper()[columnNumber_] ;
    breakEven_ = breakEven;
    assert (breakEven_ > 0.0 && breakEven_ < 1.0);
    preferredWay_ = 0;
}


// Copy constructor
CbcSimpleInteger::CbcSimpleInteger ( const CbcSimpleInteger & rhs)
        : CbcObject(rhs)

{
    columnNumber_ = rhs.columnNumber_;
    originalLower_ = rhs.originalLower_;
    originalUpper_ = rhs.originalUpper_;
    breakEven_ = rhs.breakEven_;
    preferredWay_ = rhs.preferredWay_;
}

// Clone
CbcObject *
CbcSimpleInteger::clone() const
{
    return new CbcSimpleInteger(*this);
}

// Assignment operator
CbcSimpleInteger &
CbcSimpleInteger::operator=( const CbcSimpleInteger & rhs)
{
    if (this != &rhs) {
        CbcObject::operator=(rhs);
        columnNumber_ = rhs.columnNumber_;
        originalLower_ = rhs.originalLower_;
        originalUpper_ = rhs.originalUpper_;
        breakEven_ = rhs.breakEven_;
        preferredWay_ = rhs.preferredWay_;
    }
    return *this;
}

// Destructor
CbcSimpleInteger::~CbcSimpleInteger ()
{
}
// Construct an OsiSimpleInteger object
OsiSimpleInteger *
CbcSimpleInteger::osiObject() const
{
    OsiSimpleInteger * obj = new OsiSimpleInteger(columnNumber_,
            originalLower_, originalUpper_);
    obj->setPriority(priority());
    return obj;
}
double
CbcSimpleInteger::infeasibility(const OsiBranchingInformation * info,
                                int &preferredWay) const
{
    double value = info->solution_[columnNumber_];
    value = CoinMax(value, info->lower_[columnNumber_]);
    value = CoinMin(value, info->upper_[columnNumber_]);
    double nearest = floor(value + (1.0 - breakEven_));
    assert (breakEven_ > 0.0 && breakEven_ < 1.0);
    if (nearest > value)
        preferredWay = 1;
    else
        preferredWay = -1;
    if (preferredWay_)
        preferredWay = preferredWay_;
    double weight = fabs(value - nearest);
    // normalize so weight is 0.5 at break even
    if (nearest < value)
        weight = (0.5 / breakEven_) * weight;
    else
        weight = (0.5 / (1.0 - breakEven_)) * weight;
    if (fabs(value - nearest) <= info->integerTolerance_)
        return 0.0;
    else
        return weight;
}
double
CbcSimpleInteger::feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const
{
    double value = info->solution_[columnNumber_];
#ifdef COIN_DEVELOP
    if (fabs(value - floor(value + 0.5)) > 1.0e-5)
        printf("value for %d away from integer %g\n", columnNumber_, value);
#endif
    double newValue = CoinMax(value, info->lower_[columnNumber_]);
    newValue = CoinMin(newValue, info->upper_[columnNumber_]);
    newValue = floor(newValue + 0.5);
    solver->setColLower(columnNumber_, newValue);
    solver->setColUpper(columnNumber_, newValue);
    return fabs(value - newValue);
}

/* Create an OsiSolverBranch object

This returns NULL if branch not represented by bound changes
*/
OsiSolverBranch *
CbcSimpleInteger::solverBranch(OsiSolverInterface * /*solver*/,
                               const OsiBranchingInformation * info) const
{
    double value = info->solution_[columnNumber_];
    value = CoinMax(value, info->lower_[columnNumber_]);
    value = CoinMin(value, info->upper_[columnNumber_]);
    assert (info->upper_[columnNumber_] > info->lower_[columnNumber_]);
#ifndef NDEBUG
    double nearest = floor(value + 0.5);
    assert (fabs(value - nearest) > info->integerTolerance_);
#endif
    OsiSolverBranch * branch = new OsiSolverBranch();
    branch->addBranch(columnNumber_, value);
    return branch;
}
// Creates a branching object
CbcBranchingObject *
CbcSimpleInteger::createCbcBranch(OsiSolverInterface * /*solver*/,
                                  const OsiBranchingInformation * info, int way)
{
    CbcIntegerBranchingObject * branch = new CbcIntegerBranchingObject(model_, 0, -1, 0.5);
    fillCreateBranch(branch, info, way);
    return branch;
}
// Fills in a created branching object
void
CbcSimpleInteger::fillCreateBranch(CbcIntegerBranchingObject * branch, const OsiBranchingInformation * info, int way)
{
    branch->setOriginalObject(this);
    double value = info->solution_[columnNumber_];
    value = CoinMax(value, info->lower_[columnNumber_]);
    value = CoinMin(value, info->upper_[columnNumber_]);
    assert (info->upper_[columnNumber_] > info->lower_[columnNumber_]);
    if (!info->hotstartSolution_ && priority_ != -999) {
#ifndef NDEBUG
        double nearest = floor(value + 0.5);
        assert (fabs(value - nearest) > info->integerTolerance_);
#endif
    } else if (info->hotstartSolution_) {
        double targetValue = info->hotstartSolution_[columnNumber_];
        if (way > 0)
            value = targetValue - 0.1;
        else
            value = targetValue + 0.1;
    } else {
        if (value <= info->lower_[columnNumber_])
            value += 0.1;
        else if (value >= info->upper_[columnNumber_])
            value -= 0.1;
    }
    assert (value >= info->lower_[columnNumber_] &&
            value <= info->upper_[columnNumber_]);
    branch->fillPart(columnNumber_, way, value);
}
/* Column number if single column object -1 otherwise,
   so returns >= 0
   Used by heuristics
*/
int
CbcSimpleInteger::columnNumber() const
{
    return columnNumber_;
}
/* Reset variable bounds to their original values.

    Bounds may be tightened, so it may be good to be able to set this info in object.
*/
void
CbcSimpleInteger::resetBounds(const OsiSolverInterface * solver)
{
    originalLower_ = solver->getColLower()[columnNumber_] ;
    originalUpper_ = solver->getColUpper()[columnNumber_] ;
}

/*  Change column numbers after preprocessing
 */
void
CbcSimpleInteger::resetSequenceEtc(int /*numberColumns*/,
                                   const int * originalColumns)
{
    //assert (numberColumns>0);
    int iColumn;
#if 0
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (columnNumber_ == originalColumns[iColumn])
            break;
    }
    assert (iColumn < numberColumns);
#else
    iColumn = originalColumns[columnNumber_];
    assert (iColumn >= 0);
#endif
    columnNumber_ = iColumn;
}
// This looks at solution and sets bounds to contain solution
/** More precisely: it first forces the variable within the existing
    bounds, and then tightens the bounds to fix the variable at the
    nearest integer value.
*/
void
CbcSimpleInteger::feasibleRegion()
{
    abort();
}
