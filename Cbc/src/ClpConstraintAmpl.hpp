/* $Id: ClpConstraintAmpl.hpp 1173 2009-06-04 09:44:10Z forrest $ */
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef ClpConstraintAmpl_H
#define ClpConstraintAmpl_H

#include "ClpConstraint.hpp"

//#############################################################################

/** Ampl Constraint Class

*/

class ClpConstraintAmpl : public ClpConstraint {

public:

    ///@name Stuff
    //@{


    /** Fills gradient.  If Ampl then solution may be NULL,
        also returns true value of function and offset so we can use x not deltaX in constraint
        If refresh is false then uses last solution
        Uses model for scaling
        Returns non-zero if gradient udefined at current solution
    */
    virtual int gradient(const ClpSimplex * model,
                         const double * solution,
                         double * gradient,
                         double & functionValue ,
                         double & offset,
                         bool useScaling = false,
                         bool refresh = true) const ;
    /// Resize constraint
    virtual void resize(int newNumberColumns) ;
    /// Delete columns in  constraint
    virtual void deleteSome(int numberToDelete, const int * which) ;
    /// Scale constraint
    virtual void reallyScale(const double * columnScale) ;
    /** Given a zeroed array sets nonampl columns to 1.
        Returns number of nonampl columns
     */
    virtual int markNonlinear(char * which) const ;
    /** Given a zeroed array sets possible nonzero coefficients to 1.
        Returns number of nonzeros
     */
    virtual int markNonzero(char * which) const;
    /// Say we have new primal solution - so may need to recompute
    virtual void newXValues() ;
    //@}


    ///@name Constructors and destructors
    //@{
    /// Default Constructor
    ClpConstraintAmpl();

    /// Constructor from ampl
    ClpConstraintAmpl(int row, void * amplInfo);

    /** Copy constructor .
    */
    ClpConstraintAmpl(const ClpConstraintAmpl & rhs);

    /// Assignment operator
    ClpConstraintAmpl & operator=(const ClpConstraintAmpl& rhs);

    /// Destructor
    virtual ~ClpConstraintAmpl ();

    /// Clone
    virtual ClpConstraint * clone() const;
    //@}
    ///@name Gets and sets
    //@{
    /// Number of coefficients
    virtual int numberCoefficients() const;
    /// Columns
    inline const int * column() const {
        return column_;
    }
    /// Coefficients
    inline const double * coefficient() const {
        return coefficient_;
    }
    //@}

    //---------------------------------------------------------------------------

private:
    ///@name Private member data
    /// Ampl info
    void * amplInfo_;
    /// Column
    int * column_;
    /// Coefficients
    double * coefficient_;
    /// Number of coefficients in gradient
    int numberCoefficients_;
    //@}
};

#endif

