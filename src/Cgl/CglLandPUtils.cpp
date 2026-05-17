// Copyright (C) 2005-2009, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami
//           LIF
//           CNRS, Aix-Marseille Universites
// Date:     02/23/08
//
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//---------------------------------------------------------------------------

#include "CglLandPUtils.hpp"
#include "OsiRowCut.hpp"
#include "OsiCuts.hpp"
namespace LAP
{
double
normCoef(TabRow &row, int ncols, const int * nonBasics)
{
    double res = 1;
    for (int i = 0 ; i < ncols ; i++)
        res += fabs(row[nonBasics[i]]);
    return res/(1-row.rhs);
}


/** scale the cut passed as argument using provided normalization factor*/
void scale(OsiRowCut &cut, double norma)
{
    assert(norma >0.);
    CoinPackedVector row;
    row.reserve(cut.row().getNumElements());
    for (int i = 0 ; i < cut.row().getNumElements() ; i++)
    {
        row.insert(cut.row().getIndices()[i], cut.row().getElements()[i]/norma);
    }
    cut.setLb(cut.lb()/norma);
    cut.setRow(row);
}

/** scale the cut passed as argument*/
void scale(OsiRowCut &cut)
{
    double rhs = fabs(cut.lb());
    CoinPackedVector row;
    row.reserve(cut.row().getNumElements());
    for (int i = 0 ; i < cut.row().getNumElements() ; i++)
    {
        row.insert(cut.row().getIndices()[i], cut.row().getElements()[i]/rhs);
    }
    cut.setLb(cut.lb()/rhs);
    cut.setRow(row);
}


/** Modularize row.*/
void modularizeRow(TabRow & row, const bool * integerVar)
{
    const int& n = row.getNumElements();
    const int* ind = row.getIndices();
    for (int i = 0 ; i < n ; i++)
    {
        const int &ni = ind[i];
        if (integerVar[ni])
            row[ni] = modularizedCoef(row[ni],row.rhs);
    }
}


int
Cuts::insertAll(OsiCuts & cs, CoinRelFltEq& eq)
{
    int r_val = 0;
    for (unsigned int i = 0 ; i < cuts_.size() ; i++)
    {
        if (cuts_[i] != NULL)
        {
            cs.insertIfNotDuplicate(*cuts_[i], eq);
            delete cuts_[i];
            cuts_[i] = NULL;
            r_val++;
        }
    }
    return r_val;
}

/** insert a cut for variable i and count number of cuts.*/
void
Cuts::insert(int i, OsiRowCut * cut)
{
    if (cuts_[i] == NULL) numberCuts_++;
    else
    {
        printf("Replacing cut with violation %g with one from optimal basis with violation %g.\n",
               cuts_[i]->effectiveness(), cut->effectiveness());
        delete cuts_[i];
    }
    cuts_[i] = cut;
}

}

