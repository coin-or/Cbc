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
#include "CbcClique.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"
//##############################################################################

// Default Constructor
CbcClique::CbcClique ()
        : CbcObject(),
        numberMembers_(0),
        numberNonSOSMembers_(0),
        members_(NULL),
        type_(NULL),
        cliqueType_(-1),
        slack_(-1)
{
}

// Useful constructor (which are integer indices)
CbcClique::CbcClique (CbcModel * model, int cliqueType, int numberMembers,
                      const int * which, const char * type, int identifier, int slack)
        : CbcObject(model)
{
    id_ = identifier;
    numberMembers_ = numberMembers;
    if (numberMembers_) {
        members_ = new int[numberMembers_];
        memcpy(members_, which, numberMembers_*sizeof(int));
        type_ = new char[numberMembers_];
        if (type) {
            memcpy(type_, type, numberMembers_*sizeof(char));
        } else {
            for (int i = 0; i < numberMembers_; i++)
                type_[i] = 1;
        }
    } else {
        members_ = NULL;
        type_ = NULL;
    }
    // Find out how many non sos
    int i;
    numberNonSOSMembers_ = 0;
    for (i = 0; i < numberMembers_; i++)
        if (!type_[i])
            numberNonSOSMembers_++;
    cliqueType_ = cliqueType;
    slack_ = slack;
}

// Copy constructor
CbcClique::CbcClique ( const CbcClique & rhs)
        : CbcObject(rhs)
{
    numberMembers_ = rhs.numberMembers_;
    numberNonSOSMembers_ = rhs.numberNonSOSMembers_;
    if (numberMembers_) {
        members_ = new int[numberMembers_];
        memcpy(members_, rhs.members_, numberMembers_*sizeof(int));
        type_ = new char[numberMembers_];
        memcpy(type_, rhs.type_, numberMembers_*sizeof(char));
    } else {
        members_ = NULL;
        type_ = NULL;
    }
    cliqueType_ = rhs.cliqueType_;
    slack_ = rhs.slack_;
}

// Clone
CbcObject *
CbcClique::clone() const
{
    return new CbcClique(*this);
}

// Assignment operator
CbcClique &
CbcClique::operator=( const CbcClique & rhs)
{
    if (this != &rhs) {
        CbcObject::operator=(rhs);
        delete [] members_;
        delete [] type_;
        numberMembers_ = rhs.numberMembers_;
        numberNonSOSMembers_ = rhs.numberNonSOSMembers_;
        if (numberMembers_) {
            members_ = new int[numberMembers_];
            memcpy(members_, rhs.members_, numberMembers_*sizeof(int));
            type_ = new char[numberMembers_];
            memcpy(type_, rhs.type_, numberMembers_*sizeof(char));
        } else {
            members_ = NULL;
            type_ = NULL;
        }
        cliqueType_ = rhs.cliqueType_;
        slack_ = rhs.slack_;
    }
    return *this;
}

// Destructor
CbcClique::~CbcClique ()
{
    delete [] members_;
    delete [] type_;
}
double
CbcClique::infeasibility(const OsiBranchingInformation * /*info*/,
                         int &preferredWay) const
{
    int numberUnsatis = 0, numberFree = 0;
    int j;
    const int * integer = model_->integerVariable();
    OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    double largestValue = 0.0;
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double * sort = new double[numberMembers_];

    double slackValue = 0.0;
    for (j = 0; j < numberMembers_; j++) {
        int sequence = members_[j];
        int iColumn = integer[sequence];
        double value = solution[iColumn];
        value = CoinMax(value, lower[iColumn]);
        value = CoinMin(value, upper[iColumn]);
        double nearest = floor(value + 0.5);
        double distance = fabs(value - nearest);
        if (distance > integerTolerance) {
            if (!type_[j])
                value = 1.0 - value; // non SOS
            // if slack then choose that
            if (j == slack_ && value > 0.05)
                slackValue = value;
            largestValue = CoinMax(value, largestValue);
            sort[numberUnsatis++] = -value;
        } else if (upper[iColumn] > lower[iColumn]) {
            numberFree++;
        }
    }
    preferredWay = 1;
    double otherWay = 0.0;
    if (numberUnsatis) {
        // sort
        std::sort(sort, sort + numberUnsatis);
        for (j = 0; j < numberUnsatis; j++) {
            if ((j&1) != 0)
                otherWay += -sort[j];
        }
        // Need to think more
        double value = 0.2 * numberUnsatis + 0.01 * (numberMembers_ - numberFree);
        if (fabs(largestValue - 0.5) < 0.1) {
            // close to half
            value += 0.1;
        }
        if (slackValue) {
            // branching on slack
            value += slackValue;
        }
        // scale other way
        otherWay *= value / (1.0 - otherWay);
        delete [] sort;
        return value;
    } else {
        delete [] sort;
        return 0.0; // satisfied
    }
}

// This looks at solution and sets bounds to contain solution
void
CbcClique::feasibleRegion()
{
    int j;
    const int * integer = model_->integerVariable();
    OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
#ifndef NDEBUG
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
#endif
    for (j = 0; j < numberMembers_; j++) {
        int sequence = members_[j];
        int iColumn = integer[sequence];
        double value = solution[iColumn];
        value = CoinMax(value, lower[iColumn]);
        value = CoinMin(value, upper[iColumn]);
        double nearest = floor(value + 0.5);
#ifndef NDEBUG
        double distance = fabs(value - nearest);
        assert(distance <= integerTolerance);
#endif
        solver->setColLower(iColumn, nearest);
        solver->setColUpper(iColumn, nearest);
    }
}
// Redoes data when sequence numbers change
void
CbcClique::redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns)
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
            type_[n2++] = type_[j];
        }
    }
    if (n2 < numberMembers_) {
        //printf("** SOS number of members reduced from %d to %d!\n",numberMembers_,n2);
        numberMembers_ = n2;
    }
    // Find out how many non sos
    int i;
    numberNonSOSMembers_ = 0;
    for (i = 0; i < numberMembers_; i++)
        if (!type_[i])
            numberNonSOSMembers_++;
}
CbcBranchingObject *
CbcClique::createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * /*info*/, int way)
{
    int numberUnsatis = 0;
    int j;
    int nUp = 0;
    int nDown = 0;
    int numberFree = numberMembers_;
    const int * integer = model_->integerVariable();
    //OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    int * upList = new int[numberMembers_];
    int * downList = new int[numberMembers_];
    double * sort = new double[numberMembers_];
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);

    double slackValue = 0.0;
    for (j = 0; j < numberMembers_; j++) {
        int sequence = members_[j];
        int iColumn = integer[sequence];
        double value = solution[iColumn];
        value = CoinMax(value, lower[iColumn]);
        value = CoinMin(value, upper[iColumn]);
        double nearest = floor(value + 0.5);
        double distance = fabs(value - nearest);
        if (distance > integerTolerance) {
            if (!type_[j])
                value = 1.0 - value; // non SOS
            // if slack then choose that
            if (j == slack_ && value > 0.05)
                slackValue = value;
            value = -value; // for sort
            upList[numberUnsatis] = j;
            sort[numberUnsatis++] = value;
        } else if (upper[iColumn] > lower[iColumn]) {
            upList[--numberFree] = j;
        }
    }
    assert (numberUnsatis);
    if (!slackValue) {
        // sort
        CoinSort_2(sort, sort + numberUnsatis, upList);
        // put first in up etc
        int kWay = 1;
        for (j = 0; j < numberUnsatis; j++) {
            if (kWay > 0)
                upList[nUp++] = upList[j];
            else
                downList[nDown++] = upList[j];
            kWay = -kWay;
        }
        for (j = numberFree; j < numberMembers_; j++) {
            if (kWay > 0)
                upList[nUp++] = upList[j];
            else
                downList[nDown++] = upList[j];
            kWay = -kWay;
        }
    } else {
        // put slack to 0 in first way
        nUp = 1;
        upList[0] = slack_;
        for (j = 0; j < numberUnsatis; j++) {
            downList[nDown++] = upList[j];
        }
        for (j = numberFree; j < numberMembers_; j++) {
            downList[nDown++] = upList[j];
        }
    }
    // create object
    CbcBranchingObject * branch;
    if (numberMembers_ <= 64)
        branch = new CbcCliqueBranchingObject(model_, this, way,
                                              nDown, downList, nUp, upList);
    else
        branch = new CbcLongCliqueBranchingObject(model_, this, way,
                nDown, downList, nUp, upList);
    delete [] upList;
    delete [] downList;
    delete [] sort;
    return branch;
}
