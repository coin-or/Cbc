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
/*
  Unfortunately, that comment is untrue. And there are other issues. This
  routine is clearly an unfinished work.
*/
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
    /*
      Calculate integer infeasibility and fill an array. Pick off the infeasibility
      of the slack and the max infeasibility while we're at it. You can see here
      the conversion of `non-SOS' (strong value of 0, negative coefficient) to
      `SOS' (strong value of 1, positive coefficient). Also count the number of
      variables that have integral values but are not fixed.
    */
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
    /*
      preferredWay will not change. The calculation of otherWay is an expensive
      noop --- value is ultimately unused. Same for the sort of sort. It looks like
      there was some notion of branching by splitting the set using even and odd
      indices (as opposed to first and second half).
    */
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
        /*
          Here we have the actual infeasibility calculation. Most previous work is
          discarded, and we calculate a value using various counts, adjusted by the
          max value and slack value. This is not scaled to [0, .5].
        */

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

// Default Constructor
CbcCliqueBranchingObject::CbcCliqueBranchingObject()
        : CbcBranchingObject()
{
    clique_ = NULL;
    downMask_[0] = 0;
    downMask_[1] = 0;
    upMask_[0] = 0;
    upMask_[1] = 0;
}

// Useful constructor
CbcCliqueBranchingObject::CbcCliqueBranchingObject (CbcModel * model,
        const CbcClique * clique,
        int way ,
        int numberOnDownSide, const int * down,
        int numberOnUpSide, const int * up)
        : CbcBranchingObject(model, clique->id(), way, 0.5)
{
    clique_ = clique;
    downMask_[0] = 0;
    downMask_[1] = 0;
    upMask_[0] = 0;
    upMask_[1] = 0;
    int i;
    for (i = 0; i < numberOnDownSide; i++) {
        int sequence = down[i];
        int iWord = sequence >> 5;
        int iBit = sequence - 32 * iWord;
        unsigned int k = 1 << iBit;
        downMask_[iWord] |= k;
    }
    for (i = 0; i < numberOnUpSide; i++) {
        int sequence = up[i];
        int iWord = sequence >> 5;
        int iBit = sequence - 32 * iWord;
        unsigned int k = 1 << iBit;
        upMask_[iWord] |= k;
    }
}

// Copy constructor
CbcCliqueBranchingObject::CbcCliqueBranchingObject ( const CbcCliqueBranchingObject & rhs) : CbcBranchingObject(rhs)
{
    clique_ = rhs.clique_;
    downMask_[0] = rhs.downMask_[0];
    downMask_[1] = rhs.downMask_[1];
    upMask_[0] = rhs.upMask_[0];
    upMask_[1] = rhs.upMask_[1];
}

// Assignment operator
CbcCliqueBranchingObject &
CbcCliqueBranchingObject::operator=( const CbcCliqueBranchingObject & rhs)
{
    if (this != &rhs) {
        CbcBranchingObject::operator=(rhs);
        clique_ = rhs.clique_;
        downMask_[0] = rhs.downMask_[0];
        downMask_[1] = rhs.downMask_[1];
        upMask_[0] = rhs.upMask_[0];
        upMask_[1] = rhs.upMask_[1];
    }
    return *this;
}
CbcBranchingObject *
CbcCliqueBranchingObject::clone() const
{
    return (new CbcCliqueBranchingObject(*this));
}


// Destructor
CbcCliqueBranchingObject::~CbcCliqueBranchingObject ()
{
}
double
CbcCliqueBranchingObject::branch()
{
    decrementNumberBranchesLeft();
    int iWord;
    int numberMembers = clique_->numberMembers();
    const int * which = clique_->members();
    const int * integerVariables = model_->integerVariable();
    int numberWords = (numberMembers + 31) >> 5;
    // *** for way - up means fix all those in down section
    if (way_ < 0) {
#ifdef FULL_PRINT
        printf("Down Fix ");
#endif
        for (iWord = 0; iWord < numberWords; iWord++) {
            int i;
            for (i = 0; i < 32; i++) {
                unsigned int k = 1 << i;
                if ((upMask_[iWord]&k) != 0) {
                    int iColumn = which[i+32*iWord];
#ifdef FULL_PRINT
                    printf("%d ", i + 32*iWord);
#endif
                    // fix weak way
                    if (clique_->type(i + 32*iWord))
                        model_->solver()->setColUpper(integerVariables[iColumn], 0.0);
                    else
                        model_->solver()->setColLower(integerVariables[iColumn], 1.0);
                }
            }
        }
        way_ = 1;	  // Swap direction
    } else {
#ifdef FULL_PRINT
        printf("Up Fix ");
#endif
        for (iWord = 0; iWord < numberWords; iWord++) {
            int i;
            for (i = 0; i < 32; i++) {
                unsigned int k = 1 << i;
                if ((downMask_[iWord]&k) != 0) {
                    int iColumn = which[i+32*iWord];
#ifdef FULL_PRINT
                    printf("%d ", i + 32*iWord);
#endif
                    // fix weak way
                    if (clique_->type(i + 32*iWord))
                        model_->solver()->setColUpper(integerVariables[iColumn], 0.0);
                    else
                        model_->solver()->setColLower(integerVariables[iColumn], 1.0);
                }
            }
        }
        way_ = -1;	  // Swap direction
    }
#ifdef FULL_PRINT
    printf("\n");
#endif
    return 0.0;
}
// Print what would happen
void
CbcCliqueBranchingObject::print()
{
    int iWord;
    int numberMembers = clique_->numberMembers();
    const int * which = clique_->members();
    const int * integerVariables = model_->integerVariable();
    int numberWords = (numberMembers + 31) >> 5;
    // *** for way - up means fix all those in down section
    if (way_ < 0) {
        printf("Clique - Down Fix ");
        for (iWord = 0; iWord < numberWords; iWord++) {
            int i;
            for (i = 0; i < 32; i++) {
                unsigned int k = 1 << i;
                if ((upMask_[iWord]&k) != 0) {
                    int iColumn = which[i+32*iWord];
                    printf("%d ", integerVariables[iColumn]);
                }
            }
        }
    } else {
        printf("Clique - Up Fix ");
        for (iWord = 0; iWord < numberWords; iWord++) {
            int i;
            for (i = 0; i < 32; i++) {
                unsigned int k = 1 << i;
                if ((downMask_[iWord]&k) != 0) {
                    int iColumn = which[i+32*iWord];
                    printf("%d ", integerVariables[iColumn]);
                }
            }
        }
    }
    printf("\n");
}

static inline int
CbcCompareCliques(const CbcClique* cl0, const CbcClique* cl1)
{
    if (cl0->cliqueType() < cl1->cliqueType()) {
        return -1;
    }
    if (cl0->cliqueType() > cl1->cliqueType()) {
        return 1;
    }
    if (cl0->numberMembers() != cl1->numberMembers()) {
        return cl0->numberMembers() - cl1->numberMembers();
    }
    if (cl0->numberNonSOSMembers() != cl1->numberNonSOSMembers()) {
        return cl0->numberNonSOSMembers() - cl1->numberNonSOSMembers();
    }
    return memcmp(cl0->members(), cl1->members(),
                  cl0->numberMembers() * sizeof(int));
}

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcCliqueBranchingObject::compareOriginalObject
(const CbcBranchingObject* brObj) const
{
    const CbcCliqueBranchingObject* br =
        dynamic_cast<const CbcCliqueBranchingObject*>(brObj);
    assert(br);
    return CbcCompareCliques(clique_, br->clique_);
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
CbcCliqueBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool /*replaceIfOverlap*/)
{
    const CbcCliqueBranchingObject* br =
        dynamic_cast<const CbcCliqueBranchingObject*>(brObj);
    assert(br);
    unsigned int* thisMask = way_ < 0 ? upMask_ : downMask_;
    const unsigned int* otherMask = br->way_ < 0 ? br->upMask_ : br->downMask_;
    const CoinUInt64 cl0 =
        (static_cast<CoinUInt64>(thisMask[0]) << 32) | thisMask[1];
    const CoinUInt64 cl1 =
        (static_cast<CoinUInt64>(otherMask[0]) << 32) | otherMask[1];
    if (cl0 == cl1) {
        return CbcRangeSame;
    }
    const CoinUInt64 cl_intersection = (cl0 & cl1);
    if (cl_intersection == cl0) {
        return CbcRangeSuperset;
    }
    if (cl_intersection == cl1) {
        return CbcRangeSubset;
    }
    const CoinUInt64 cl_xor = (cl0 ^ cl1);
    if (cl_intersection == 0 && cl_xor == 0) {
        return CbcRangeDisjoint;
    }
    const CoinUInt64 cl_union = (cl0 | cl1);
    thisMask[0] = static_cast<unsigned int>(cl_union >> 32);
    thisMask[1] = static_cast<unsigned int>(cl_union & 0xffffffff);
    return CbcRangeOverlap;
}

//##############################################################################

// Default Constructor
CbcLongCliqueBranchingObject::CbcLongCliqueBranchingObject()
        : CbcBranchingObject()
{
    clique_ = NULL;
    downMask_ = NULL;
    upMask_ = NULL;
}

// Useful constructor
CbcLongCliqueBranchingObject::CbcLongCliqueBranchingObject (CbcModel * model,
        const CbcClique * clique,
        int way ,
        int numberOnDownSide, const int * down,
        int numberOnUpSide, const int * up)
        : CbcBranchingObject(model, clique->id(), way, 0.5)
{
    clique_ = clique;
    int numberMembers = clique_->numberMembers();
    int numberWords = (numberMembers + 31) >> 5;
    downMask_ = new unsigned int [numberWords];
    upMask_ = new unsigned int [numberWords];
    memset(downMask_, 0, numberWords*sizeof(unsigned int));
    memset(upMask_, 0, numberWords*sizeof(unsigned int));
    int i;
    for (i = 0; i < numberOnDownSide; i++) {
        int sequence = down[i];
        int iWord = sequence >> 5;
        int iBit = sequence - 32 * iWord;
        unsigned int k = 1 << iBit;
        downMask_[iWord] |= k;
    }
    for (i = 0; i < numberOnUpSide; i++) {
        int sequence = up[i];
        int iWord = sequence >> 5;
        int iBit = sequence - 32 * iWord;
        unsigned int k = 1 << iBit;
        upMask_[iWord] |= k;
    }
}

// Copy constructor
CbcLongCliqueBranchingObject::CbcLongCliqueBranchingObject ( const CbcLongCliqueBranchingObject & rhs) : CbcBranchingObject(rhs)
{
    clique_ = rhs.clique_;
    if (rhs.downMask_) {
        int numberMembers = clique_->numberMembers();
        int numberWords = (numberMembers + 31) >> 5;
        downMask_ = new unsigned int [numberWords];
        memcpy(downMask_, rhs.downMask_, numberWords*sizeof(unsigned int));
        upMask_ = new unsigned int [numberWords];
        memcpy(upMask_, rhs.upMask_, numberWords*sizeof(unsigned int));
    } else {
        downMask_ = NULL;
        upMask_ = NULL;
    }
}

// Assignment operator
CbcLongCliqueBranchingObject &
CbcLongCliqueBranchingObject::operator=( const CbcLongCliqueBranchingObject & rhs)
{
    if (this != &rhs) {
        CbcBranchingObject::operator=(rhs);
        clique_ = rhs.clique_;
        delete [] downMask_;
        delete [] upMask_;
        if (rhs.downMask_) {
            int numberMembers = clique_->numberMembers();
            int numberWords = (numberMembers + 31) >> 5;
            downMask_ = new unsigned int [numberWords];
            memcpy(downMask_, rhs.downMask_, numberWords*sizeof(unsigned int));
            upMask_ = new unsigned int [numberWords];
            memcpy(upMask_, rhs.upMask_, numberWords*sizeof(unsigned int));
        } else {
            downMask_ = NULL;
            upMask_ = NULL;
        }
    }
    return *this;
}
CbcBranchingObject *
CbcLongCliqueBranchingObject::clone() const
{
    return (new CbcLongCliqueBranchingObject(*this));
}


// Destructor
CbcLongCliqueBranchingObject::~CbcLongCliqueBranchingObject ()
{
    delete [] downMask_;
    delete [] upMask_;
}
double
CbcLongCliqueBranchingObject::branch()
{
    decrementNumberBranchesLeft();
    int iWord;
    int numberMembers = clique_->numberMembers();
    const int * which = clique_->members();
    const int * integerVariables = model_->integerVariable();
    int numberWords = (numberMembers + 31) >> 5;
    // *** for way - up means fix all those in down section
    if (way_ < 0) {
#ifdef FULL_PRINT
        printf("Down Fix ");
#endif
        for (iWord = 0; iWord < numberWords; iWord++) {
            int i;
            for (i = 0; i < 32; i++) {
                unsigned int k = 1 << i;
                if ((upMask_[iWord]&k) != 0) {
                    int iColumn = which[i+32*iWord];
#ifdef FULL_PRINT
                    printf("%d ", i + 32*iWord);
#endif
                    // fix weak way
                    if (clique_->type(i + 32*iWord))
                        model_->solver()->setColUpper(integerVariables[iColumn], 0.0);
                    else
                        model_->solver()->setColLower(integerVariables[iColumn], 1.0);
                }
            }
        }
        way_ = 1;	  // Swap direction
    } else {
#ifdef FULL_PRINT
        printf("Up Fix ");
#endif
        for (iWord = 0; iWord < numberWords; iWord++) {
            int i;
            for (i = 0; i < 32; i++) {
                unsigned int k = 1 << i;
                if ((downMask_[iWord]&k) != 0) {
                    int iColumn = which[i+32*iWord];
#ifdef FULL_PRINT
                    printf("%d ", i + 32*iWord);
#endif
                    // fix weak way
                    if (clique_->type(i + 32*iWord))
                        model_->solver()->setColUpper(integerVariables[iColumn], 0.0);
                    else
                        model_->solver()->setColLower(integerVariables[iColumn], 1.0);
                }
            }
        }
        way_ = -1;	  // Swap direction
    }
#ifdef FULL_PRINT
    printf("\n");
#endif
    return 0.0;
}
void
CbcLongCliqueBranchingObject::print()
{
    int iWord;
    int numberMembers = clique_->numberMembers();
    const int * which = clique_->members();
    const int * integerVariables = model_->integerVariable();
    int numberWords = (numberMembers + 31) >> 5;
    // *** for way - up means fix all those in down section
    if (way_ < 0) {
        printf("Clique - Down Fix ");
        for (iWord = 0; iWord < numberWords; iWord++) {
            int i;
            for (i = 0; i < 32; i++) {
                unsigned int k = 1 << i;
                if ((upMask_[iWord]&k) != 0) {
                    int iColumn = which[i+32*iWord];
                    printf("%d ", integerVariables[iColumn]);
                }
            }
        }
    } else {
        printf("Clique - Up Fix ");
        for (iWord = 0; iWord < numberWords; iWord++) {
            int i;
            for (i = 0; i < 32; i++) {
                unsigned int k = 1 << i;
                if ((downMask_[iWord]&k) != 0) {
                    int iColumn = which[i+32*iWord];
                    printf("%d ", integerVariables[iColumn]);
                }
            }
        }
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
CbcLongCliqueBranchingObject::compareOriginalObject
(const CbcBranchingObject* brObj) const
{
    const CbcLongCliqueBranchingObject* br =
        dynamic_cast<const CbcLongCliqueBranchingObject*>(brObj);
    assert(br);
    return CbcCompareCliques(clique_, br->clique_);
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
CbcLongCliqueBranchingObject::compareBranchingObject
(const CbcBranchingObject* brObj, const bool /*replaceIfOverlap*/)
{
    const CbcLongCliqueBranchingObject* br =
        dynamic_cast<const CbcLongCliqueBranchingObject*>(brObj);
    assert(br);
    const int numberMembers = clique_->numberMembers();
    const int numberWords = (numberMembers + 31) >> 5;
    unsigned int* thisMask = way_ < 0 ? upMask_ : downMask_;
    const unsigned int* otherMask = br->way_ < 0 ? br->upMask_ : br->downMask_;

    if (memcmp(thisMask, otherMask, numberWords * sizeof(unsigned int)) == 0) {
        return CbcRangeSame;
    }
    bool canBeSuperset = true;
    bool canBeSubset = true;
    int i;
    for (i = numberWords - 1; i >= 0 && (canBeSuperset || canBeSubset); --i) {
        const unsigned int both = (thisMask[i] & otherMask[i]);
        canBeSuperset &= (both == thisMask[i]);
        canBeSubset &= (both == otherMask[i]);
    }
    if (canBeSuperset) {
        return CbcRangeSuperset;
    }
    if (canBeSubset) {
        return CbcRangeSubset;
    }

    for (i = numberWords - 1; i >= 0; --i) {
        if ((thisMask[i] ^ otherMask[i]) != 0) {
            break;
        }
    }
    if (i == -1) { // complement
        return CbcRangeDisjoint;
    }
    // must be overlap
    for (i = numberWords - 1; i >= 0; --i) {
        thisMask[i] |= otherMask[i];
    }
    return CbcRangeOverlap;
}

