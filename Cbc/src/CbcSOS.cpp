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
#include "CbcSOS.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

//##############################################################################

// Default Constructor
CbcSOS::CbcSOS ()
        : CbcObject(),
        members_(NULL),
        weights_(NULL),
        shadowEstimateDown_(1.0),
        shadowEstimateUp_(1.0),
        downDynamicPseudoRatio_(0.0),
        upDynamicPseudoRatio_(0.0),
        numberTimesDown_(0),
        numberTimesUp_(0),
        numberMembers_(0),
        sosType_(-1),
        integerValued_(false)
{
}

// Useful constructor (which are indices)
CbcSOS::CbcSOS (CbcModel * model,  int numberMembers,
                const int * which, const double * weights, int identifier, int type)
        : CbcObject(model),
        shadowEstimateDown_(1.0),
        shadowEstimateUp_(1.0),
        downDynamicPseudoRatio_(0.0),
        upDynamicPseudoRatio_(0.0),
        numberTimesDown_(0),
        numberTimesUp_(0),
        numberMembers_(numberMembers),
        sosType_(type)
{
    id_ = identifier;
    integerValued_ = type == 1;
    if (integerValued_) {
        // check all members integer
        OsiSolverInterface * solver = model->solver();
        if (solver) {
            for (int i = 0; i < numberMembers_; i++) {
                if (!solver->isInteger(which[i]))
                    integerValued_ = false;
            }
        } else {
            // can't tell
            integerValued_ = false;
        }
    }
    if (numberMembers_) {
        members_ = new int[numberMembers_];
        weights_ = new double[numberMembers_];
        memcpy(members_, which, numberMembers_*sizeof(int));
        if (weights) {
            memcpy(weights_, weights, numberMembers_*sizeof(double));
        } else {
            for (int i = 0; i < numberMembers_; i++)
                weights_[i] = i;
        }
        // sort so weights increasing
        CoinSort_2(weights_, weights_ + numberMembers_, members_);
        double last = -COIN_DBL_MAX;
        int i;
        for (i = 0; i < numberMembers_; i++) {
            double possible = CoinMax(last + 1.0e-10, weights_[i]);
            weights_[i] = possible;
            last = possible;
        }
    } else {
        members_ = NULL;
        weights_ = NULL;
    }
    assert (sosType_ > 0 && sosType_ < 3);
}

// Copy constructor
CbcSOS::CbcSOS ( const CbcSOS & rhs)
        : CbcObject(rhs)
{
    shadowEstimateDown_ = rhs.shadowEstimateDown_;
    shadowEstimateUp_ = rhs.shadowEstimateUp_;
    downDynamicPseudoRatio_ = rhs.downDynamicPseudoRatio_;
    upDynamicPseudoRatio_ = rhs.upDynamicPseudoRatio_;
    numberTimesDown_ = rhs.numberTimesDown_;
    numberTimesUp_ = rhs.numberTimesUp_;
    numberMembers_ = rhs.numberMembers_;
    sosType_ = rhs.sosType_;
    integerValued_ = rhs.integerValued_;
    if (numberMembers_) {
        members_ = new int[numberMembers_];
        weights_ = new double[numberMembers_];
        memcpy(members_, rhs.members_, numberMembers_*sizeof(int));
        memcpy(weights_, rhs.weights_, numberMembers_*sizeof(double));
    } else {
        members_ = NULL;
        weights_ = NULL;
    }
}

// Clone
CbcObject *
CbcSOS::clone() const
{
    return new CbcSOS(*this);
}

// Assignment operator
CbcSOS &
CbcSOS::operator=( const CbcSOS & rhs)
{
    if (this != &rhs) {
        CbcObject::operator=(rhs);
        delete [] members_;
        delete [] weights_;
        shadowEstimateDown_ = rhs.shadowEstimateDown_;
        shadowEstimateUp_ = rhs.shadowEstimateUp_;
        downDynamicPseudoRatio_ = rhs.downDynamicPseudoRatio_;
        upDynamicPseudoRatio_ = rhs.upDynamicPseudoRatio_;
        numberTimesDown_ = rhs.numberTimesDown_;
        numberTimesUp_ = rhs.numberTimesUp_;
        numberMembers_ = rhs.numberMembers_;
        sosType_ = rhs.sosType_;
        integerValued_ = rhs.integerValued_;
        if (numberMembers_) {
            members_ = new int[numberMembers_];
            weights_ = new double[numberMembers_];
            memcpy(members_, rhs.members_, numberMembers_*sizeof(int));
            memcpy(weights_, rhs.weights_, numberMembers_*sizeof(double));
        } else {
            members_ = NULL;
            weights_ = NULL;
        }
    }
    return *this;
}

// Destructor
CbcSOS::~CbcSOS ()
{
    delete [] members_;
    delete [] weights_;
}
double
CbcSOS::infeasibility(const OsiBranchingInformation * info,
                      int &preferredWay) const
{
    int j;
    int firstNonZero = -1;
    int lastNonZero = -1;
    OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    //const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    //double largestValue=0.0;
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double weight = 0.0;
    double sum = 0.0;

    // check bounds etc
    double lastWeight = -1.0e100;
    for (j = 0; j < numberMembers_; j++) {
        int iColumn = members_[j];
        if (lastWeight >= weights_[j] - 1.0e-7)
            throw CoinError("Weights too close together in SOS", "infeasibility", "CbcSOS");
        double value = CoinMax(0.0, solution[iColumn]);
        sum += value;
        if (value > integerTolerance && upper[iColumn]) {
            // Possibly due to scaling a fixed variable might slip through
            if (value > upper[iColumn]) {
                value = upper[iColumn];
                // Could change to #ifdef CBC_DEBUG
#ifndef NDEBUG
                if (model_->messageHandler()->logLevel() > 2)
                    printf("** Variable %d (%d) has value %g and upper bound of %g\n",
                           iColumn, j, value, upper[iColumn]);
#endif
            }
            weight += weights_[j] * value;
            if (firstNonZero < 0)
                firstNonZero = j;
            lastNonZero = j;
        }
    }
    preferredWay = 1;
    if (lastNonZero - firstNonZero >= sosType_) {
        // find where to branch
        assert (sum > 0.0);
        weight /= sum;
        if (info->defaultDual_ >= 0.0 && info->usefulRegion_ && info->columnStart_) {
            assert (sosType_ == 1);
            int iWhere;
            for (iWhere = firstNonZero; iWhere < lastNonZero - 1; iWhere++) {
                if (weight < weights_[iWhere+1]) {
                    break;
                }
            }
            int jColumnDown = members_[iWhere];
            int jColumnUp = members_[iWhere+1];
            int n = 0;
            CoinBigIndex j;
            double objMove = info->objective_[jColumnDown];
            for (j = info->columnStart_[jColumnDown];
                    j < info->columnStart_[jColumnDown] + info->columnLength_[jColumnDown]; j++) {
                double value = info->elementByColumn_[j];
                int iRow = info->row_[j];
                info->indexRegion_[n++] = iRow;
                info->usefulRegion_[iRow] = value;
            }
            for (iWhere = firstNonZero; iWhere < lastNonZero; iWhere++) {
                int jColumn = members_[iWhere];
                double solValue = info->solution_[jColumn];
                if (!solValue)
                    continue;
                objMove -= info->objective_[jColumn] * solValue;
                for (j = info->columnStart_[jColumn];
                        j < info->columnStart_[jColumn] + info->columnLength_[jColumn]; j++) {
                    double value = -info->elementByColumn_[j] * solValue;
                    int iRow = info->row_[j];
                    double oldValue = info->usefulRegion_[iRow];
                    if (!oldValue) {
                        info->indexRegion_[n++] = iRow;
                    } else {
                        value += oldValue;
                        if (!value)
                            value = 1.0e-100;
                    }
                    info->usefulRegion_[iRow] = value;
                }
            }
            const double * pi = info->pi_;
            const double * activity = info->rowActivity_;
            const double * lower = info->rowLower_;
            const double * upper = info->rowUpper_;
            double tolerance = info->primalTolerance_;
            double direction = info->direction_;
            shadowEstimateDown_ = objMove * direction;
            bool infeasible = false;
            for (int k = 0; k < n; k++) {
                int iRow = info->indexRegion_[k];
                double movement = info->usefulRegion_[iRow];
                // not this time info->usefulRegion_[iRow]=0.0;
                if (lower[iRow] < -1.0e20)
                    assert (pi[iRow] <= 1.0e-3);
                if (upper[iRow] > 1.0e20)
                    assert (pi[iRow] >= -1.0e-3);
                double valueP = pi[iRow] * direction;
                // if move makes infeasible then make at least default
                double newValue = activity[iRow] + movement;
                if (newValue > upper[iRow] + tolerance || newValue < lower[iRow] - tolerance) {
                    shadowEstimateDown_ += fabs(movement) * CoinMax(fabs(valueP), info->defaultDual_);
                    infeasible = true;
                }
            }
            if (shadowEstimateDown_ < info->integerTolerance_) {
                if (!infeasible) {
                    shadowEstimateDown_ = 1.0e-10;
#ifdef COIN_DEVELOP
                    printf("zero pseudoShadowPrice\n");
#endif
                } else
                    shadowEstimateDown_ = info->integerTolerance_;
            }
            // And other way
            // take off
            objMove -= info->objective_[jColumnDown];
            for (j = info->columnStart_[jColumnDown];
                    j < info->columnStart_[jColumnDown] + info->columnLength_[jColumnDown]; j++) {
                double value = -info->elementByColumn_[j];
                int iRow = info->row_[j];
                double oldValue = info->usefulRegion_[iRow];
                if (!oldValue) {
                    info->indexRegion_[n++] = iRow;
                } else {
                    value += oldValue;
                    if (!value)
                        value = 1.0e-100;
                }
                info->usefulRegion_[iRow] = value;
            }
            // add on
            objMove += info->objective_[jColumnUp];
            for (j = info->columnStart_[jColumnUp];
                    j < info->columnStart_[jColumnUp] + info->columnLength_[jColumnUp]; j++) {
                double value = info->elementByColumn_[j];
                int iRow = info->row_[j];
                double oldValue = info->usefulRegion_[iRow];
                if (!oldValue) {
                    info->indexRegion_[n++] = iRow;
                } else {
                    value += oldValue;
                    if (!value)
                        value = 1.0e-100;
                }
                info->usefulRegion_[iRow] = value;
            }
            shadowEstimateUp_ = objMove * direction;
            infeasible = false;
            for (int k = 0; k < n; k++) {
                int iRow = info->indexRegion_[k];
                double movement = info->usefulRegion_[iRow];
                info->usefulRegion_[iRow] = 0.0;
                if (lower[iRow] < -1.0e20)
                    assert (pi[iRow] <= 1.0e-3);
                if (upper[iRow] > 1.0e20)
                    assert (pi[iRow] >= -1.0e-3);
                double valueP = pi[iRow] * direction;
                // if move makes infeasible then make at least default
                double newValue = activity[iRow] + movement;
                if (newValue > upper[iRow] + tolerance || newValue < lower[iRow] - tolerance) {
                    shadowEstimateUp_ += fabs(movement) * CoinMax(fabs(valueP), info->defaultDual_);
                    infeasible = true;
                }
            }
            if (shadowEstimateUp_ < info->integerTolerance_) {
                if (!infeasible) {
                    shadowEstimateUp_ = 1.0e-10;
#ifdef COIN_DEVELOP
                    printf("zero pseudoShadowPrice\n");
#endif
                } else
                    shadowEstimateUp_ = info->integerTolerance_;
            }
            // adjust
            double downCost = shadowEstimateDown_;
            double upCost = shadowEstimateUp_;
            if (numberTimesDown_)
                downCost *= downDynamicPseudoRatio_ /
                            static_cast<double> (numberTimesDown_);
            if (numberTimesUp_)
                upCost *= upDynamicPseudoRatio_ /
                          static_cast<double> (numberTimesUp_);
#define WEIGHT_AFTER 0.7
#define WEIGHT_BEFORE 0.1
            int stateOfSearch = model_->stateOfSearch() % 10;
            double returnValue = 0.0;
            double minValue = CoinMin(downCost, upCost);
            double maxValue = CoinMax(downCost, upCost);
            if (stateOfSearch <= 2) {
                // no branching solution
                returnValue = WEIGHT_BEFORE * minValue + (1.0 - WEIGHT_BEFORE) * maxValue;
            } else {
                returnValue = WEIGHT_AFTER * minValue + (1.0 - WEIGHT_AFTER) * maxValue;
            }
#ifdef PRINT_SHADOW
            printf("%d id - down %d %g up %d %g shadow %g, %g returned %g\n",
                   id_, numberTimesDown_, downDynamicPseudoRatio_,
                   numberTimesUp_, upDynamicPseudoRatio_, shadowEstimateDown_,
                   shadowEstimateUp_, returnValue);
#endif
            return returnValue;
        } else {
            double value = lastNonZero - firstNonZero + 1;
            value *= 0.5 / static_cast<double> (numberMembers_);
            return value;
        }
    } else {
        return 0.0; // satisfied
    }
}

// This looks at solution and sets bounds to contain solution
void
CbcSOS::feasibleRegion()
{
    int j;
    int firstNonZero = -1;
    int lastNonZero = -1;
    OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    //const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double weight = 0.0;
    double sum = 0.0;

    for (j = 0; j < numberMembers_; j++) {
        int iColumn = members_[j];
        double value = CoinMax(0.0, solution[iColumn]);
        sum += value;
        if (value > integerTolerance && upper[iColumn]) {
            weight += weights_[j] * value;
            if (firstNonZero < 0)
                firstNonZero = j;
            lastNonZero = j;
        }
    }
    assert (lastNonZero - firstNonZero < sosType_) ;
    for (j = 0; j < firstNonZero; j++) {
        int iColumn = members_[j];
        solver->setColUpper(iColumn, 0.0);
    }
    for (j = lastNonZero + 1; j < numberMembers_; j++) {
        int iColumn = members_[j];
        solver->setColUpper(iColumn, 0.0);
    }
}
// Redoes data when sequence numbers change
void
CbcSOS::redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns)
{
    model_ = model;
    int n2 = 0;
    for (int j = 0; j < numberMembers_; j++) {
        int iColumn = members_[j];
        int i;
        for (i = 0; i < numberColumns; i++) {
            if (originalColumns[i] == iColumn)
                break;
        }
        if (i < numberColumns) {
            members_[n2] = i;
            weights_[n2++] = weights_[j];
        }
    }
    if (n2 < numberMembers_) {
        //printf("** SOS number of members reduced from %d to %d!\n",numberMembers_,n2);
        numberMembers_ = n2;
    }
}
CbcBranchingObject *
CbcSOS::createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * /*info*/, int way)
{
    int j;
    const double * solution = model_->testSolution();
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
    //OsiSolverInterface * solver = model_->solver();
    const double * upper = solver->getColUpper();
    int firstNonFixed = -1;
    int lastNonFixed = -1;
    int firstNonZero = -1;
    int lastNonZero = -1;
    double weight = 0.0;
    double sum = 0.0;
    for (j = 0; j < numberMembers_; j++) {
        int iColumn = members_[j];
        if (upper[iColumn]) {
            double value = CoinMax(0.0, solution[iColumn]);
            sum += value;
            if (firstNonFixed < 0)
                firstNonFixed = j;
            lastNonFixed = j;
            if (value > integerTolerance) {
                weight += weights_[j] * value;
                if (firstNonZero < 0)
                    firstNonZero = j;
                lastNonZero = j;
            }
        }
    }
    assert (lastNonZero - firstNonZero >= sosType_) ;
    // find where to branch
    assert (sum > 0.0);
    weight /= sum;
    int iWhere;
    double separator = 0.0;
    for (iWhere = firstNonZero; iWhere < lastNonZero; iWhere++)
        if (weight < weights_[iWhere+1])
            break;
    if (sosType_ == 1) {
        // SOS 1
        separator = 0.5 * (weights_[iWhere] + weights_[iWhere+1]);
    } else {
        // SOS 2
        if (iWhere == firstNonFixed)
            iWhere++;;
        if (iWhere == lastNonFixed - 1)
            iWhere = lastNonFixed - 2;
        separator = weights_[iWhere+1];
    }
    // create object
    CbcBranchingObject * branch;
    branch = new CbcSOSBranchingObject(model_, this, way, separator);
    branch->setOriginalObject(this);
    return branch;
}
/* Pass in information on branch just done and create CbcObjectUpdateData instance.
   If object does not need data then backward pointer will be NULL.
   Assumes can get information from solver */
CbcObjectUpdateData
CbcSOS::createUpdateInformation(const OsiSolverInterface * solver,
                                const CbcNode * node,
                                const CbcBranchingObject * branchingObject)
{
    double originalValue = node->objectiveValue();
    int originalUnsatisfied = node->numberUnsatisfied();
    double objectiveValue = solver->getObjValue() * solver->getObjSense();
    int unsatisfied = 0;
    int i;
    //might be base model - doesn't matter
    int numberIntegers = model_->numberIntegers();;
    const double * solution = solver->getColSolution();
    //const double * lower = solver->getColLower();
    //const double * upper = solver->getColUpper();
    double change = CoinMax(0.0, objectiveValue - originalValue);
    int iStatus;
    if (solver->isProvenOptimal())
        iStatus = 0; // optimal
    else if (solver->isIterationLimitReached()
             && !solver->isDualObjectiveLimitReached())
        iStatus = 2; // unknown
    else
        iStatus = 1; // infeasible

    bool feasible = iStatus != 1;
    if (feasible) {
        double integerTolerance =
            model_->getDblParam(CbcModel::CbcIntegerTolerance);
        const int * integerVariable = model_->integerVariable();
        for (i = 0; i < numberIntegers; i++) {
            int j = integerVariable[i];
            double value = solution[j];
            double nearest = floor(value + 0.5);
            if (fabs(value - nearest) > integerTolerance)
                unsatisfied++;
        }
    }
    int way = branchingObject->way();
    way = - way; // because after branch so moved on
    double value = branchingObject->value();
    CbcObjectUpdateData newData (this, way,
                                 change, iStatus,
                                 originalUnsatisfied - unsatisfied, value);
    newData.originalObjective_ = originalValue;
    // Solvers know about direction
    double direction = solver->getObjSense();
    solver->getDblParam(OsiDualObjectiveLimit, newData.cutoff_);
    newData.cutoff_ *= direction;
    return newData;
}
// Update object by CbcObjectUpdateData
void
CbcSOS::updateInformation(const CbcObjectUpdateData & data)
{
    bool feasible = data.status_ != 1;
    int way = data.way_;
    //double value = data.branchingValue_;
    double originalValue = data.originalObjective_;
    double change = data.change_;
    if (way < 0) {
        // down
        if (!feasible) {
            double distanceToCutoff = 0.0;
            //double objectiveValue = model_->getCurrentMinimizationObjValue();
            distanceToCutoff =  model_->getCutoff()  - originalValue;
            if (distanceToCutoff < 1.0e20)
                change = distanceToCutoff * 2.0;
            else
                change = (downDynamicPseudoRatio_ * shadowEstimateDown_ + 1.0e-3) * 10.0;
        }
        change = CoinMax(1.0e-12 * (1.0 + fabs(originalValue)), change);
#ifdef PRINT_SHADOW
        if (numberTimesDown_)
            printf("Updating id %d - down change %g (true %g) - ndown %d estimated change %g - raw shadow estimate %g\n",
                   id_, change, data.change_, numberTimesDown_, shadowEstimateDown_*
                   (downDynamicPseudoRatio_ / ((double) numberTimesDown_)),
                   shadowEstimateDown_);
        else
            printf("Updating id %d - down change %g (true %g) - shadow estimate %g\n",
                   id_, change, data.change_, shadowEstimateDown_);
#endif
        numberTimesDown_++;
        downDynamicPseudoRatio_ += change / shadowEstimateDown_;
    } else {
        // up
        if (!feasible) {
            double distanceToCutoff = 0.0;
            //double objectiveValue = model_->getCurrentMinimizationObjValue();
            distanceToCutoff =  model_->getCutoff()  - originalValue;
            if (distanceToCutoff < 1.0e20)
                change = distanceToCutoff * 2.0;
            else
                change = (upDynamicPseudoRatio_ * shadowEstimateUp_ + 1.0e-3) * 10.0;
        }
        change = CoinMax(1.0e-12 * (1.0 + fabs(originalValue)), change);
#ifdef PRINT_SHADOW
        if (numberTimesUp_)
            printf("Updating id %d - up change %g (true %g) - nup %d estimated change %g - raw shadow estimate %g\n",
                   id_, change, data.change_, numberTimesUp_, shadowEstimateUp_*
                   (upDynamicPseudoRatio_ / ((double) numberTimesUp_)),
                   shadowEstimateUp_);
        else
            printf("Updating id %d - up change %g (true %g) - shadow estimate %g\n",
                   id_, change, data.change_, shadowEstimateUp_);
#endif
        numberTimesUp_++;
        upDynamicPseudoRatio_ += change / shadowEstimateUp_;
    }
}

/* Create an OsiSolverBranch object

This returns NULL if branch not represented by bound changes
*/
OsiSolverBranch *
CbcSOS::solverBranch() const
{
    int j;
    const double * solution = model_->testSolution();
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
    OsiSolverInterface * solver = model_->solver();
    const double * upper = solver->getColUpper();
    int firstNonFixed = -1;
    int lastNonFixed = -1;
    int firstNonZero = -1;
    int lastNonZero = -1;
    double weight = 0.0;
    double sum = 0.0;
    double * fix = new double[numberMembers_];
    int * which = new int[numberMembers_];
    for (j = 0; j < numberMembers_; j++) {
        int iColumn = members_[j];
        // fix all on one side or other (even if fixed)
        fix[j] = 0.0;
        which[j] = iColumn;
        if (upper[iColumn]) {
            double value = CoinMax(0.0, solution[iColumn]);
            sum += value;
            if (firstNonFixed < 0)
                firstNonFixed = j;
            lastNonFixed = j;
            if (value > integerTolerance) {
                weight += weights_[j] * value;
                if (firstNonZero < 0)
                    firstNonZero = j;
                lastNonZero = j;
            }
        }
    }
    assert (lastNonZero - firstNonZero >= sosType_) ;
    // find where to branch
    assert (sum > 0.0);
    weight /= sum;
    // down branch fixes ones above weight to 0
    int iWhere;
    int iDownStart = 0;
    int iUpEnd = 0;
    for (iWhere = firstNonZero; iWhere < lastNonZero; iWhere++)
        if (weight < weights_[iWhere+1])
            break;
    if (sosType_ == 1) {
        // SOS 1
        iUpEnd = iWhere + 1;
        iDownStart = iUpEnd;
    } else {
        // SOS 2
        if (iWhere == firstNonFixed)
            iWhere++;;
        if (iWhere == lastNonFixed - 1)
            iWhere = lastNonFixed - 2;
        iUpEnd = iWhere + 1;
        iDownStart = iUpEnd + 1;
    }
    //
    OsiSolverBranch * branch = new OsiSolverBranch();
    branch->addBranch(-1, 0, NULL, NULL, numberMembers_ - iDownStart, which + iDownStart, fix);
    branch->addBranch(1, 0, NULL, NULL, iUpEnd, which, fix);
    delete [] fix;
    delete [] which;
    return branch;
}
// Construct an OsiSOS object
OsiSOS *
CbcSOS::osiObject(const OsiSolverInterface * solver) const
{
    OsiSOS * obj = new OsiSOS(solver, numberMembers_, members_, weights_, sosType_);
    obj->setPriority(priority());
    return obj;
}
