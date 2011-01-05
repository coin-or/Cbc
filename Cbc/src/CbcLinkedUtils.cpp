// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/* $Id$ */

/*! \file CbcAugmentClpSimplex.cpp
    \brief Hooks to Ampl (for CbcLinked)

    This code is a condensation of ClpAmplStuff.cpp, renamed to better
    reflect its current place in cbc.

  The code here had ties to NEW_STYLE_SOLVER code. During the 091209 Watson
  meeting, NEW_STYLE_SOLVER code was eliminated. The code here was condensed
  from ClpAmplStuff.cpp. The hook into CbcLinked is loadNonLinear. Once you
  bring that in, all the rest follows. Still, we're down about 400 lines of
  code. In the process, it appears that ClpAmplObjective.cpp was never needed
  here; the code was hooked into ClpAmplStuff.cpp.  --lh, 091209 --
*/

#include "ClpConfig.h"
#include "CbcConfig.h"
#ifdef COIN_HAS_ASL
#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpFactorization.hpp"
#include "ClpSimplex.hpp"
#include "ClpAmplObjective.hpp"
#include "ClpConstraintAmpl.hpp"
#include "ClpMessage.hpp"
#include "CoinUtilsConfig.h"
#include "CoinHelperFunctions.hpp"
#include "CoinWarmStartBasis.hpp"
#include "OsiSolverInterface.hpp"
#include "Cbc_ampl.h"
#include "CoinTime.hpp"
#include "CglStored.hpp"
#include "CoinModel.hpp"
#include "CbcLinked.hpp"

extern "C" {
    //# include "getstub.h"
# include "asl_pfgh.h"
}

// stolen from IPopt with changes
typedef struct {
    double obj_sign_;
    ASL_pfgh * asl_;
    double * non_const_x_;
    int * column_; // for jacobian
    int * rowStart_;
    double * gradient_;
    double * constraintValues_;
    int nz_h_full_; // number of nonzeros in hessian
    int nerror_;
    bool objval_called_with_current_x_;
    bool conval_called_with_current_x_;
    bool jacval_called_with_current_x_;
} CbcAmplInfo;

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpAmplObjective::ClpAmplObjective ()
        : ClpObjective()
{
    type_ = 12;
    objective_ = NULL;
    amplObjective_ = NULL;
    gradient_ = NULL;
    offset_ = 0.0;
}

bool get_constraints_linearity(void * amplInfo, int  n,
                               int * const_types)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    ASL_pfgh* asl = info->asl_;
    //check that n is good
    assert(n == n_con);
    // check that there are no network constraints
    assert(nlnc == 0 && lnc == 0);
    //the first nlc constraints are non linear the rest is linear
    int i;
    for (i = 0; i < nlc; i++) {
        const_types[i] = 1;
    }
    // the rest is linear
    for (i = nlc; i < n_con; i++)
        const_types[i] = 0;
    return true;
}
static bool internal_objval(CbcAmplInfo * info , double & obj_val)
{
    ASL_pfgh* asl = info->asl_;
    info->objval_called_with_current_x_ = false; // in case the call below fails

    if (n_obj == 0) {
        obj_val = 0;
        info->objval_called_with_current_x_ = true;
        return true;
    }  else {
        double  retval = objval(0, info->non_const_x_, (fint*)info->nerror_);
        if (!info->nerror_) {
            obj_val = info->obj_sign_ * retval;
            info->objval_called_with_current_x_ = true;
            return true;
        } else {
            abort();
        }
    }

    return false;
}

static bool internal_conval(CbcAmplInfo * info , double * g)
{
    ASL_pfgh* asl = info->asl_;
    info->conval_called_with_current_x_ = false; // in case the call below fails
    assert (g);

    conval(info->non_const_x_, g, (fint*)info->nerror_);

    if (!info->nerror_) {
        info->conval_called_with_current_x_ = true;
        return true;
    } else {
        abort();
    }
    return false;
}

static bool apply_new_x(CbcAmplInfo * info  , bool new_x, int  n, const double * x)
{
    ASL_pfgh* asl = info->asl_;

    if (new_x) {
        // update the flags so these methods are called
        // before evaluating the hessian
        info->conval_called_with_current_x_ = false;
        info->objval_called_with_current_x_ = false;
        info->jacval_called_with_current_x_ = false;

        //copy the data to the non_const_x_
        if (!info->non_const_x_) {
            info->non_const_x_ = new double [n];
        }

        for (int  i = 0; i < n; i++) {
            info->non_const_x_[i] = x[i];
        }

        // tell ampl that we have a new x
        xknowne(info->non_const_x_, (fint*)info->nerror_);
        return info->nerror_ ? false : true;
    }

    return true;
}

static bool eval_f(void * amplInfo, int  n, const double * x, bool new_x, double & obj_value)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    if (!apply_new_x(info, new_x, n, x)) {
        return false;
    }

    return internal_objval(info, obj_value);
}

static bool eval_grad_f(void * amplInfo, int  n, const double * x, bool new_x, double * grad_f)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    ASL_pfgh* asl = info->asl_;
    if (!apply_new_x(info, new_x, n, x)) {
        return false;
    }
    int i;

    if (n_obj == 0) {
        for (i = 0; i < n; i++) {
            grad_f[i] = 0.;
        }
    } else {
        objgrd(0, info->non_const_x_, grad_f, (fint*)info->nerror_);
        if (info->nerror_) {
            return false;
        }

        if (info->obj_sign_ == -1) {
            for (i = 0; i < n; i++) {
                grad_f[i] = -grad_f[i];
            }
        }
    }
    return true;
}

static bool eval_g(void * amplInfo, int  n, const double * x, bool new_x, double * g)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
#ifndef NDEBUG
    ASL_pfgh* asl = info->asl_;
#endif
    // warning: n_var is a macro that assumes we have a variable called asl
    assert(n == n_var);

    if (!apply_new_x(info, new_x, n, x)) {
        return false;
    }

    return internal_conval(info, g);
}

static bool eval_jac_g(void * amplInfo, int  n, const double * x, bool new_x,
                       double * values)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    ASL_pfgh* asl = info->asl_;
    assert(n == n_var);

    assert (values);
    if (!apply_new_x(info, new_x, n, x)) {
        return false;
    }

    jacval(info->non_const_x_, values, (fint*)info->nerror_);
    if (!info->nerror_) {
        return true;
    } else {
        abort();
    }
    return false;
}
//-------------------------------------------------------------------
// Useful Constructor
//-------------------------------------------------------------------
ClpAmplObjective::ClpAmplObjective (void * amplInfo)
        : ClpObjective()
{
    type_ = 12;
    activated_ = 1;
    gradient_ = NULL;
    objective_ = NULL;
    offset_ = 0.0;
    amplObjective_ = amplInfo;
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpAmplObjective::ClpAmplObjective (const ClpAmplObjective & rhs)
        : ClpObjective(rhs)
{
    amplObjective_ = rhs.amplObjective_;
    offset_ = rhs.offset_;
    type_ = rhs.type_;
    if (!amplObjective_) {
        objective_ = NULL;
        gradient_ = NULL;
    } else {
        CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
        ASL_pfgh* asl = info->asl_;

        int numberColumns = n_var;;
        if (rhs.objective_) {
            objective_ = new double [numberColumns];
            memcpy(objective_, rhs.objective_, numberColumns*sizeof(double));
        } else {
            objective_ = NULL;
        }
        if (rhs.gradient_) {
            gradient_ = new double [numberColumns];
            memcpy(gradient_, rhs.gradient_, numberColumns*sizeof(double));
        } else {
            gradient_ = NULL;
        }
    }
}


//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpAmplObjective::~ClpAmplObjective ()
{
    delete [] objective_;
    delete [] gradient_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpAmplObjective &
ClpAmplObjective::operator=(const ClpAmplObjective & rhs)
{
    if (this != &rhs) {
        delete [] objective_;
        delete [] gradient_;
        amplObjective_ = rhs.amplObjective_;
        offset_ = rhs.offset_;
        type_ = rhs.type_;
        if (!amplObjective_) {
            objective_ = NULL;
            gradient_ = NULL;
        } else {
            CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
            ASL_pfgh* asl = info->asl_;

            int numberColumns = n_var;;
            if (rhs.objective_) {
                objective_ = new double [numberColumns];
                memcpy(objective_, rhs.objective_, numberColumns*sizeof(double));
            } else {
                objective_ = NULL;
            }
            if (rhs.gradient_) {
                gradient_ = new double [numberColumns];
                memcpy(gradient_, rhs.gradient_, numberColumns*sizeof(double));
            } else {
                gradient_ = NULL;
            }
        }
    }
    return *this;
}

// Returns gradient
double *
ClpAmplObjective::gradient(const ClpSimplex * model,
                           const double * solution, double & offset, bool refresh,
                           int includeLinear)
{
    if (model)
        assert (model->optimizationDirection() == 1.0);
    bool scaling = false;
    if (model && (model->rowScale() ||
                  model->objectiveScale() != 1.0 || model->optimizationDirection() != 1.0))
        scaling = true;
    const double * cost = NULL;
    if (model)
        cost = model->costRegion();
    if (!cost) {
        // not in solve
        cost = objective_;
        scaling = false;
    }
    assert (!scaling);
    if (!amplObjective_ || !solution || !activated_) {
        offset = offset_;
        return objective_;
    } else {
        if (refresh || !gradient_) {
            CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
            ASL_pfgh* asl = info->asl_;
            int numberColumns = n_var;;

            if (!gradient_)
                gradient_ = new double[numberColumns];
            assert (solution);
            eval_grad_f(amplObjective_, numberColumns, solution, true, gradient_);
            // Is this best way?
            double objValue = 0.0;
            eval_f(amplObjective_, numberColumns, solution, false, objValue);
            double objValue2 = 0.0;
            for (int i = 0; i < numberColumns; i++)
                objValue2 += gradient_[i] * solution[i];
            offset_ = objValue2 - objValue; // or other way???
            if (model && model->optimizationDirection() != 1.0) {
                offset *= model->optimizationDirection();
                for (int i = 0; i < numberColumns; i++)
                    gradient_[i] *= -1.0;
            }
        }
        offset = offset_;
        return gradient_;
    }
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpObjective * ClpAmplObjective::clone() const
{
    return new ClpAmplObjective(*this);
}
// Resize objective
void
ClpAmplObjective::resize(int newNumberColumns)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
    ASL_pfgh* asl = info->asl_;
    int numberColumns = n_var;;
    if (numberColumns != newNumberColumns) {
        abort();
    }

}
// Delete columns in  objective
void
ClpAmplObjective::deleteSome(int numberToDelete, const int * which)
{
    if (numberToDelete)
        abort();
}
/* Returns reduced gradient.Returns an offset (to be added to current one).
 */
double
ClpAmplObjective::reducedGradient(ClpSimplex * model, double * region,
                                  bool useFeasibleCosts)
{
    int numberRows = model->numberRows();
    int numberColumns = model->numberColumns();

    //work space
    CoinIndexedVector  * workSpace = model->rowArray(0);

    CoinIndexedVector arrayVector;
    arrayVector.reserve(numberRows + 1);

    int iRow;
#ifdef CLP_DEBUG
    workSpace->checkClear();
#endif
    double * array = arrayVector.denseVector();
    int * index = arrayVector.getIndices();
    int number = 0;
    const double * costNow = gradient(model, model->solutionRegion(), offset_,
                                      true, useFeasibleCosts ? 2 : 1);
    double * cost = model->costRegion();
    const int * pivotVariable = model->pivotVariable();
    for (iRow = 0; iRow < numberRows; iRow++) {
        int iPivot = pivotVariable[iRow];
        double value;
        if (iPivot < numberColumns)
            value = costNow[iPivot];
        else if (!useFeasibleCosts)
            value = cost[iPivot];
        else
            value = 0.0;
        if (value) {
            array[iRow] = value;
            index[number++] = iRow;
        }
    }
    arrayVector.setNumElements(number);

    // Btran basic costs
    model->factorization()->updateColumnTranspose(workSpace, &arrayVector);
    double * work = workSpace->denseVector();
    ClpFillN(work, numberRows, 0.0);
    // now look at dual solution
    double * rowReducedCost = region + numberColumns;
    double * dual = rowReducedCost;
    const double * rowCost = cost + numberColumns;
    for (iRow = 0; iRow < numberRows; iRow++) {
        dual[iRow] = array[iRow];
    }
    double * dj = region;
    ClpDisjointCopyN(costNow, numberColumns, dj);

    model->transposeTimes(-1.0, dual, dj);
    for (iRow = 0; iRow < numberRows; iRow++) {
        // slack
        double value = dual[iRow];
        value += rowCost[iRow];
        rowReducedCost[iRow] = value;
    }
    return offset_;
}
/* Returns step length which gives minimum of objective for
   solution + theta * change vector up to maximum theta.

   arrays are numberColumns+numberRows
*/
double
ClpAmplObjective::stepLength(ClpSimplex * model,
                             const double * solution,
                             const double * change,
                             double maximumTheta,
                             double & currentObj,
                             double & predictedObj,
                             double & thetaObj)
{
    // Assume convex
    CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
    ASL_pfgh* asl = info->asl_;

    int numberColumns = n_var;;
    double * tempSolution = new double [numberColumns];
    double * tempGradient = new double [numberColumns];
    // current
    eval_f(amplObjective_, numberColumns, solution, true, currentObj);
    double objA = currentObj;
    double thetaA = 0.0;
    // at maximum
    int i;
    for (i = 0; i < numberColumns; i++)
        tempSolution[i] = solution[i] + maximumTheta * change[i];
    eval_f(amplObjective_, numberColumns, tempSolution, true, thetaObj);
    double objC = thetaObj;
    double thetaC = maximumTheta;
    double objB = 0.5 * (objA + objC);
    double thetaB = 0.5 * maximumTheta;
    double gradientNorm = 1.0e6;
    while (gradientNorm > 1.0e-6 && thetaC - thetaA > 1.0e-8) {
        for (i = 0; i < numberColumns; i++)
            tempSolution[i] = solution[i] + thetaB * change[i];
        eval_grad_f(amplObjective_, numberColumns, tempSolution, true, tempGradient);
        eval_f(amplObjective_, numberColumns, tempSolution, false, objB);
        double changeObj = 0.0;
        gradientNorm = 0.0;
        for (i = 0; i < numberColumns; i++) {
            changeObj += tempGradient[i] * change[i];
            gradientNorm += tempGradient[i] * tempGradient[i];
        }
        gradientNorm = fabs(changeObj) / sqrt(gradientNorm);
        // Should try and get quadratic convergence by interpolation
        if (changeObj < 0.0) {
            // increasing is good
            thetaA = thetaB;
        } else {
            // decreasing is good
            thetaC = thetaB;
        }
        thetaB = 0.5 * (thetaA + thetaC);
    }
    delete [] tempSolution;
    delete [] tempGradient;
    predictedObj = objB;
    return thetaB;
}
// Return objective value (without any ClpModel offset) (model may be NULL)
double
ClpAmplObjective::objectiveValue(const ClpSimplex * model, const double * solution) const
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
    ASL_pfgh* asl = info->asl_;

    int numberColumns = n_var;;
    // current
    double currentObj = 0.0;
    eval_f(amplObjective_, numberColumns, solution, true, currentObj);
    return currentObj;
}
// Scale objective
void
ClpAmplObjective::reallyScale(const double * columnScale)
{
    abort();
}
/* Given a zeroed array sets nonlinear columns to 1.
   Returns number of nonlinear columns
*/
int
ClpAmplObjective::markNonlinear(char * which)
{
    int iColumn;
    CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
    ASL_pfgh* asl = info->asl_;
    int nonLinear = CoinMax(nlvc, nlvo);
    for (iColumn = 0; iColumn < nonLinear; iColumn++) {
        which[iColumn] = 1;
    }
    int numberNonLinearColumns = 0;
    int numberColumns = n_var;;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (which[iColumn])
            numberNonLinearColumns++;
    }
    return numberNonLinearColumns;
}
// Say we have new primal solution - so may need to recompute
void
ClpAmplObjective::newXValues()
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
    info->conval_called_with_current_x_ = false;
    info->objval_called_with_current_x_ = false;
    info->jacval_called_with_current_x_ = false;
}

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################
//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpConstraintAmpl::ClpConstraintAmpl ()
        : ClpConstraint()
{
    type_ = 3;
    column_ = NULL;
    coefficient_ = NULL;
    numberCoefficients_ = 0;
    amplInfo_ = NULL;
}

//-------------------------------------------------------------------
// Useful Constructor
//-------------------------------------------------------------------
ClpConstraintAmpl::ClpConstraintAmpl (int row, void * amplInfo)
        : ClpConstraint()
{
    type_ = 3;
    rowNumber_ = row;
    amplInfo_ = amplInfo;
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo_;
#ifndef NDEBUG
    ASL_pfgh* asl = info->asl_;
#endif
    // warning: nlc is a macro that assumes we have a variable called asl
    assert (rowNumber_ < nlc);
    numberCoefficients_ = info->rowStart_[rowNumber_+1] - info->rowStart_[rowNumber_];
    column_ = CoinCopyOfArray(info->column_ + info->rowStart_[rowNumber_], numberCoefficients_);
    coefficient_ = new double [numberCoefficients_];;
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpConstraintAmpl::ClpConstraintAmpl (const ClpConstraintAmpl & rhs)
        : ClpConstraint(rhs)
{
    numberCoefficients_ = rhs.numberCoefficients_;
    column_ = CoinCopyOfArray(rhs.column_, numberCoefficients_);
    coefficient_ = CoinCopyOfArray(rhs.coefficient_, numberCoefficients_);
}


//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpConstraintAmpl::~ClpConstraintAmpl ()
{
    delete [] column_;
    delete [] coefficient_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpConstraintAmpl &
ClpConstraintAmpl::operator=(const ClpConstraintAmpl & rhs)
{
    if (this != &rhs) {
        delete [] column_;
        delete [] coefficient_;
        numberCoefficients_ = rhs.numberCoefficients_;
        column_ = CoinCopyOfArray(rhs.column_, numberCoefficients_);
        coefficient_ = CoinCopyOfArray(rhs.coefficient_, numberCoefficients_);
    }
    return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpConstraint * ClpConstraintAmpl::clone() const
{
    return new ClpConstraintAmpl(*this);
}

// Returns gradient
int
ClpConstraintAmpl::gradient(const ClpSimplex * model,
                            const double * solution,
                            double * gradient,
                            double & functionValue,
                            double & offset,
                            bool useScaling,
                            bool refresh) const
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo_;
    ASL_pfgh* asl = info->asl_;
    int numberColumns = n_var;;
    // If not done then do all
    if (!info->jacval_called_with_current_x_) {
        bool getStuff = eval_g(amplInfo_, numberColumns, solution, true, info->constraintValues_);
        assert (getStuff);
        getStuff = eval_jac_g(amplInfo_, numberColumns, solution, false, info->gradient_);
        assert (getStuff);
        info->jacval_called_with_current_x_ = getStuff;
    }
    if (refresh || !lastGradient_) {
        functionValue_ = info->constraintValues_[rowNumber_];
        offset_ = functionValue_; // sign??
        if (!lastGradient_)
            lastGradient_ = new double[numberColumns];
        CoinZeroN(lastGradient_, numberColumns);
        assert (!(model && model->rowScale() && useScaling));
        int i;
        int start = info->rowStart_[rowNumber_];
        assert (numberCoefficients_ == info->rowStart_[rowNumber_+1] - start);
        for (i = 0; i < numberCoefficients_; i++) {
            int iColumn = column_[i];
            double valueS = solution[iColumn];
            double valueG = info->gradient_[start+i];
            lastGradient_[iColumn] = valueG;
            offset_ -= valueS * valueG;
        }
    }
    functionValue = functionValue_;
    offset = offset_;
    memcpy(gradient, lastGradient_, numberColumns*sizeof(double));
    return 0;
}
// Resize constraint
void
ClpConstraintAmpl::resize(int newNumberColumns)
{
    abort();
}
// Delete columns in  constraint
void
ClpConstraintAmpl::deleteSome(int numberToDelete, const int * which)
{
    if (numberToDelete) {
        abort();
    }
}
// Scale constraint
void
ClpConstraintAmpl::reallyScale(const double * columnScale)
{
    abort();
}
/* Given a zeroed array sets nonlinear columns to 1.
   Returns number of nonlinear columns
*/
int
ClpConstraintAmpl::markNonlinear(char * which) const
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo_;
    ASL_pfgh* asl = info->asl_;
    int iColumn;
    int numberNon = 0;
    int nonLinear = CoinMax(nlvc, nlvo);
    for (iColumn = 0; iColumn < numberCoefficients_; iColumn++) {
        int jColumn = column_[iColumn];
        if (jColumn < nonLinear) {
            which[jColumn] = 1;
            numberNon++;
        }
    }
    return numberNon;
}
/* Given a zeroed array sets possible nonzero coefficients to 1.
   Returns number of nonzeros
*/
int
ClpConstraintAmpl::markNonzero(char * which) const
{
    int iColumn;
    for (iColumn = 0; iColumn < numberCoefficients_; iColumn++) {
        which[column_[iColumn]] = 1;
    }
    return numberCoefficients_;
}
// Number of coefficients
int
ClpConstraintAmpl::numberCoefficients() const
{
    return numberCoefficients_;
}
// Say we have new primal solution - so may need to recompute
void
ClpConstraintAmpl::newXValues()
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo_;
    info->conval_called_with_current_x_ = false;
    info->objval_called_with_current_x_ = false;
    info->jacval_called_with_current_x_ = false;
}

/* Load nonlinear part of problem from AMPL info
   Returns 0 if linear
   1 if quadratic objective
   2 if quadratic constraints
   3 if nonlinear objective
   4 if nonlinear constraints
   -1 on failure
*/
int
ClpSimplex::loadNonLinear(void * amplInfo, int & numberConstraints,
                          ClpConstraint ** & constraints)
{
    numberConstraints = 0;
    constraints = NULL;
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    ASL_pfgh* asl = info->asl_;
    // For moment don't say quadratic
    int type = 0;
    if (nlo + nlc) {
        // nonlinear
        if (!nlc) {
            type = 3;
            delete objective_;
            objective_ = new ClpAmplObjective(amplInfo);
        } else {
            type = 4;
            numberConstraints = nlc;
            constraints = new ClpConstraint * [numberConstraints];
            if (nlo) {
                delete objective_;
                objective_ = new ClpAmplObjective(amplInfo);
            }
            for (int i = 0; i < numberConstraints; i++) {
                constraints[i] = new ClpConstraintAmpl(i, amplInfo);
            }
        }
    }
    return type;
}
#else
#include "ClpSimplex.hpp"
#include "ClpConstraint.hpp"
int
ClpSimplex::loadNonLinear(void * , int & ,
                          ClpConstraint ** & )
{
    abort();
    return 0;
}
#endif

