// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

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
#include "CbcNWay.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

//##############################################################################

// Default Constructor
CbcNWay::CbcNWay ()
        : CbcObject(),
        numberMembers_(0),
        members_(NULL),
        consequence_(NULL)
{
}

// Useful constructor (which are integer indices)
CbcNWay::CbcNWay (CbcModel * model, int numberMembers,
                  const int * which, int identifier)
        : CbcObject(model)
{
    id_ = identifier;
    numberMembers_ = numberMembers;
    if (numberMembers_) {
        members_ = new int[numberMembers_];
        memcpy(members_, which, numberMembers_*sizeof(int));
    } else {
        members_ = NULL;
    }
    consequence_ = NULL;
}

// Copy constructor
CbcNWay::CbcNWay ( const CbcNWay & rhs)
        : CbcObject(rhs)
{
    numberMembers_ = rhs.numberMembers_;
    consequence_ = NULL;
    if (numberMembers_) {
        members_ = new int[numberMembers_];
        memcpy(members_, rhs.members_, numberMembers_*sizeof(int));
        if (rhs.consequence_) {
            consequence_ = new CbcConsequence * [numberMembers_];
            for (int i = 0; i < numberMembers_; i++) {
                if (rhs.consequence_[i])
                    consequence_[i] = rhs.consequence_[i]->clone();
                else
                    consequence_[i] = NULL;
            }
        }
    } else {
        members_ = NULL;
    }
}

// Clone
CbcObject *
CbcNWay::clone() const
{
    return new CbcNWay(*this);
}

// Assignment operator
CbcNWay &
CbcNWay::operator=( const CbcNWay & rhs)
{
    if (this != &rhs) {
        CbcObject::operator=(rhs);
        delete [] members_;
        numberMembers_ = rhs.numberMembers_;
        if (consequence_) {
            for (int i = 0; i < numberMembers_; i++)
                delete consequence_[i];
            delete [] consequence_;
            consequence_ = NULL;
        }
        if (numberMembers_) {
            members_ = new int[numberMembers_];
            memcpy(members_, rhs.members_, numberMembers_*sizeof(int));
        } else {
            members_ = NULL;
        }
        if (rhs.consequence_) {
            consequence_ = new CbcConsequence * [numberMembers_];
            for (int i = 0; i < numberMembers_; i++) {
                if (rhs.consequence_[i])
                    consequence_[i] = rhs.consequence_[i]->clone();
                else
                    consequence_[i] = NULL;
            }
        }
    }
    return *this;
}

// Destructor
CbcNWay::~CbcNWay ()
{
    delete [] members_;
    if (consequence_) {
        for (int i = 0; i < numberMembers_; i++)
            delete consequence_[i];
        delete [] consequence_;
    }
}
// Set up a consequence for a single member
void
CbcNWay::setConsequence(int iColumn, const CbcConsequence & consequence)
{
    if (!consequence_) {
        consequence_ = new CbcConsequence * [numberMembers_];
        for (int i = 0; i < numberMembers_; i++)
            consequence_[i] = NULL;
    }
    for (int i = 0; i < numberMembers_; i++) {
        if (members_[i] == iColumn) {
            consequence_[i] = consequence.clone();
            break;
        }
    }
}

// Applies a consequence for a single member
void
CbcNWay::applyConsequence(int iSequence, int state) const
{
    assert (state == -9999 || state == 9999);
    if (consequence_) {
        CbcConsequence * consequence = consequence_[iSequence];
        if (consequence)
            consequence->applyToSolver(model_->solver(), state);
    }
}
double
CbcNWay::infeasibility(const OsiBranchingInformation * /*info*/,
                       int &preferredWay) const
{
    int numberUnsatis = 0;
    int j;
    OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    double largestValue = 0.0;

    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);

    for (j = 0; j < numberMembers_; j++) {
        int iColumn = members_[j];
        double value = solution[iColumn];
        value = CoinMax(value, lower[iColumn]);
        value = CoinMin(value, upper[iColumn]);
        double distance = CoinMin(value - lower[iColumn], upper[iColumn] - value);
        if (distance > integerTolerance) {
            numberUnsatis++;
            largestValue = CoinMax(distance, largestValue);
        }
    }
    preferredWay = 1;
    if (numberUnsatis) {
        return largestValue;
    } else {
        return 0.0; // satisfied
    }
}

// This looks at solution and sets bounds to contain solution
void
CbcNWay::feasibleRegion()
{
    int j;
    OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    double integerTolerance =
        model_->getDblParam(CbcModel::CbcIntegerTolerance);
    for (j = 0; j < numberMembers_; j++) {
        int iColumn = members_[j];
        double value = solution[iColumn];
        value = CoinMax(value, lower[iColumn]);
        value = CoinMin(value, upper[iColumn]);
        if (value >= upper[iColumn] - integerTolerance) {
            solver->setColLower(iColumn, upper[iColumn]);
        } else {
            assert (value <= lower[iColumn] + integerTolerance);
            solver->setColUpper(iColumn, lower[iColumn]);
        }
    }
}
// Redoes data when sequence numbers change
void
CbcNWay::redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns)
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
            consequence_[n2++] = consequence_[j];
        } else {
            delete consequence_[j];
        }
    }
    if (n2 < numberMembers_) {
        printf("** NWay number of members reduced from %d to %d!\n", numberMembers_, n2);
        numberMembers_ = n2;
    }
}
CbcBranchingObject *
CbcNWay::createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * /*info*/, int /*way*/)
{
    int numberFree = 0;
    int j;

    //OsiSolverInterface * solver = model_->solver();
    const double * solution = model_->testSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    int * list = new int[numberMembers_];
    double * sort = new double[numberMembers_];

    for (j = 0; j < numberMembers_; j++) {
        int iColumn = members_[j];
        double value = solution[iColumn];
        value = CoinMax(value, lower[iColumn]);
        value = CoinMin(value, upper[iColumn]);
        if (upper[iColumn] > lower[iColumn]) {
            double distance = upper[iColumn] - value;
            list[numberFree] = j;
            sort[numberFree++] = distance;
        }
    }
    assert (numberFree);
    // sort
    CoinSort_2(sort, sort + numberFree, list);
    // create object
    CbcBranchingObject * branch;
    branch = new CbcNWayBranchingObject(model_, this, numberFree, list);
    branch->setOriginalObject(this);
    delete [] list;
    delete [] sort;
    return branch;
}

// Default Constructor
CbcNWayBranchingObject::CbcNWayBranchingObject()
        : CbcBranchingObject()
{
    order_ = NULL;
    object_ = NULL;
    numberInSet_ = 0;
    way_ = 0;
}

// Useful constructor
CbcNWayBranchingObject::CbcNWayBranchingObject (CbcModel * model,
        const CbcNWay * nway,
        int number, const int * order)
        : CbcBranchingObject(model, nway->id(), -1, 0.5)
{
    numberBranches_ = number;
    order_ = new int [number];
    object_ = nway;
    numberInSet_ = number;
    memcpy(order_, order, number*sizeof(int));
}

// Copy constructor
CbcNWayBranchingObject::CbcNWayBranchingObject ( const CbcNWayBranchingObject & rhs) : CbcBranchingObject(rhs)
{
    numberInSet_ = rhs.numberInSet_;
    object_ = rhs.object_;
    if (numberInSet_) {
        order_ = new int [numberInSet_];
        memcpy(order_, rhs.order_, numberInSet_*sizeof(int));
    } else {
        order_ = NULL;
    }
}

// Assignment operator
CbcNWayBranchingObject &
CbcNWayBranchingObject::operator=( const CbcNWayBranchingObject & rhs)
{
    if (this != &rhs) {
        CbcBranchingObject::operator=(rhs);
        object_ = rhs.object_;
        delete [] order_;
        numberInSet_ = rhs.numberInSet_;
        if (numberInSet_) {
            order_ = new int [numberInSet_];
            memcpy(order_, rhs.order_, numberInSet_*sizeof(int));
        } else {
            order_ = NULL;
        }
    }
    return *this;
}
CbcBranchingObject *
CbcNWayBranchingObject::clone() const
{
    return (new CbcNWayBranchingObject(*this));
}


// Destructor
CbcNWayBranchingObject::~CbcNWayBranchingObject ()
{
    delete [] order_;
}
double
CbcNWayBranchingObject::branch()
{
    int which = branchIndex_;
    branchIndex_++;
    assert (numberBranchesLeft() >= 0);
    if (which == 0) {
        // first branch so way_ may mean something
        assert (way_ == -1 || way_ == 1);
        if (way_ == -1)
            which++;
    } else if (which == 1) {
        // second branch so way_ may mean something
        assert (way_ == -1 || way_ == 1);
        if (way_ == -1)
            which--;
        // switch way off
        way_ = 0;
    }
    const double * lower = model_->solver()->getColLower();
    const double * upper = model_->solver()->getColUpper();
    const int * members = object_->members();
    for (int j = 0; j < numberInSet_; j++) {
        int iSequence = order_[j];
        int iColumn = members[iSequence];
        if (j != which) {
            model_->solver()->setColUpper(iColumn, lower[iColumn]);
            //model_->solver()->setColLower(iColumn,lower[iColumn]);
            assert (lower[iColumn] > -1.0e20);
            // apply any consequences
            object_->applyConsequence(iSequence, -9999);
        } else {
            model_->solver()->setColLower(iColumn, upper[iColumn]);
            //model_->solver()->setColUpper(iColumn,upper[iColumn]);
#ifdef FULL_PRINT
            printf("Up Fix %d to %g\n", iColumn, upper[iColumn]);
#endif
            assert (upper[iColumn] < 1.0e20);
            // apply any consequences
            object_->applyConsequence(iSequence, 9999);
        }
    }
    return 0.0;
}
void
CbcNWayBranchingObject::print()
{
    printf("NWay - Up Fix ");
    const int * members = object_->members();
    for (int j = 0; j < way_; j++) {
        int iColumn = members[order_[j]];
        printf("%d ", iColumn);
    }
    printf("\n");
}

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcNWayBranchingObject::compareOriginalObject
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
CbcNWayBranchingObject::compareBranchingObject
(const CbcBranchingObject* /*brObj*/, const bool /*replaceIfOverlap*/)
{
    throw("must implement");
}

