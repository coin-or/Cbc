// Copyright (C) 2005-2009, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami
//           LIF
//           CNRS, Aix-Marseille Universites
// Date:     02/23/08
//
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//---------------------------------------------------------------------------

#include "CglLandPTabRow.hpp"
#include "CglLandPSimplex.hpp"
namespace LAP
{
void
TabRow::print(std::ostream & os, int width, const int * nonBasics,
              int m)
{
    os.width(3);
    os.precision(4);
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    os<<"idx: ";
    const double * dense = denseVector();
    for (int j = 0 ; j < m ; j++)
    {
        os.width(width);
        os.setf(std::ios_base::right, std::ios_base::adjustfield);
        os<<nonBasics[j]<<" ";
    }

    os<<std::endl;
    os.width(3);
    os.precision(4);
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    os<< num <<": ";
    for (int j = 0 ; j < m ; j++)
    {
        os.width(width);
        os.precision(3);
        //      os.setf(std::ios_base::fixed, std::ios_base::floatfield);
        os.setf(std::ios_base::right, std::ios_base::adjustfield);
        os<<dense[nonBasics[j]]<<" ";
    }

    os.width(width);
    os.precision(4);
    //    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    os<<rhs;

    os<<std::endl;

}

bool
TabRow::operator==(const TabRow &r) const
{
    return CoinIndexedVector::operator==(r);
}


/** Modularize row.*/
void
TabRow::modularize(const bool * integerVar)
{
    const int& n = getNumElements();
    const int* ind = getIndices();
    double * el = denseVector();
    for (int i = 0 ; i < n ; i++)
    {
        const int &ni = ind[i];
        if (integerVar[ni])
        {
            el[ni] = modularizedCoef(el[ni], rhs);
        }
    }
    modularized_ = true;
}
}/* Ends namespace LAP.*/
