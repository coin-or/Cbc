/* $Id: ClpAmplObjective.hpp 1173 2009-06-04 09:44:10Z forrest $ */
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef ClpAmplObjective_H
#define ClpAmplObjective_H

#include "ClpObjective.hpp"
#include "CoinPackedMatrix.hpp"

//#############################################################################

/** Ampl Objective Class

*/

class ClpAmplObjective : public ClpObjective {

public:

    ///@name Stuff
    //@{

    /** Returns gradient.  If Ampl then solution may be NULL,
        also returns an offset (to be added to current one)
        If refresh is false then uses last solution
        Uses model for scaling
        includeLinear 0 - no, 1 as is, 2 as feasible
    */
    virtual double * gradient(const ClpSimplex * model,
                              const double * solution, double & offset, bool refresh,
                              int includeLinear = 2);
    /// Resize objective
    /** Returns reduced gradient.Returns an offset (to be added to current one).
    */
    virtual double reducedGradient(ClpSimplex * model, double * region,
                                   bool useFeasibleCosts);
    /** Returns step length which gives minimum of objective for
        solution + theta * change vector up to maximum theta.

        arrays are numberColumns+numberRows
        Also sets current objective, predicted and at maximumTheta
    */
    virtual double stepLength(ClpSimplex * model,
                              const double * solution,
                              const double * change,
                              double maximumTheta,
                              double & currentObj,
                              double & predictedObj,
                              double & thetaObj);
    /// Return objective value (without any ClpModel offset) (model may be NULL)
    virtual double objectiveValue(const ClpSimplex * model, const double * solution) const ;
    virtual void resize(int newNumberColumns) ;
    /// Delete columns in  objective
    virtual void deleteSome(int numberToDelete, const int * which) ;
    /// Scale objective
    virtual void reallyScale(const double * columnScale) ;
    /** Given a zeroed array sets nonlinear columns to 1.
        Returns number of nonlinear columns
     */
    virtual int markNonlinear(char * which);

    /// Say we have new primal solution - so may need to recompute
    virtual void newXValues() ;
    //@}


    ///@name Constructors and destructors
    //@{
    /// Default Constructor
    ClpAmplObjective();

    /// Constructor from ampl info
    ClpAmplObjective(void * amplInfo);

    /** Copy constructor .
    */
    ClpAmplObjective(const ClpAmplObjective & rhs);

    /// Assignment operator
    ClpAmplObjective & operator=(const ClpAmplObjective& rhs);

    /// Destructor
    virtual ~ClpAmplObjective ();

    /// Clone
    virtual ClpObjective * clone() const;

    //@}
    ///@name Gets and sets
    //@{
    /// Linear objective
    double * linearObjective() const;
    //@}

    //---------------------------------------------------------------------------

private:
    ///@name Private member data
    /// Saved offset
    double offset_;
    /// Ampl info
    void * amplObjective_;
    /// Objective
    double * objective_;
    /// Gradient
    double * gradient_;
    //@}
};

#endif

