// Copyright (C) 2005-2009, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami
//           LIF
//           CNRS, Aix-Marseille Universites
// Date:     02/23/08
//
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//---------------------------------------------------------------------------

#ifndef CglLandPTabRow_H
#define CglLandPTabRow_H

#include "CoinIndexedVector.hpp"
#include <iostream>

namespace LAP
{
class CglLandPSimplex;
struct TabRow: public CoinIndexedVector
{
    /** Row number.*/
    int num;
    /** Row right-hand-side.*/
    double rhs;
    /** Row of what?*/
    const CglLandPSimplex * si_;

    /** Flag to say if row is modularized.*/
    bool modularized_;

    TabRow():
            CoinIndexedVector(), si_(NULL), modularized_(false) {}

    TabRow(const CglLandPSimplex *si):
            CoinIndexedVector(), num(-1), rhs(0), si_(si), modularized_(false) {}

    TabRow(const TabRow & source):CoinIndexedVector(source),
            num(source.num), rhs(source.rhs), si_(source.si_)
    {
    }

    TabRow& operator=(const TabRow & r)
    {
        if (this != &r)
        {
            CoinIndexedVector::operator=(r);
            num = r.num;
            rhs = r.rhs;
            si_ = r.si_;
        }
        return *this;
    }

    bool operator==(const TabRow &r) const;
    ~TabRow()
    {
    }

    void modularize(const bool * integerVar);

    void print(std::ostream & os, int width = 9, const int * nonBasics = NULL,
               int m = 0);
    inline
    const double& operator[](const int &index) const
    {
        return denseVector()[index];
    }

    inline
    double& operator[](const int &index)
    {
        return denseVector()[index];
    }
};
}/* Ends LAP Namespace.*/

#endif

