// Name:     CglParam.cpp
// Author:   Francois Margot
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     11/24/06
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
// Copyright (C) 2006, Francois Margot and others.  All Rights Reserved.
//---------------------------------------------------------------------------

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>

#include "CoinPragma.hpp"
#include "CglParam.hpp"

/***********************************************************************/
void CglParam::setINFINIT(const double inf)
{
  if (inf > 0)
    INFINIT = inf;
} /* setINFINIT */

/***********************************************************************/
void CglParam::setEPS(const double eps)
{
  if (eps >= 0)
    EPS = eps;
} /* setEPS */

/***********************************************************************/
void CglParam::setEPS_COEFF(const double eps_c)
{
  if (eps_c >= 0)
    EPS_COEFF = eps_c;
} /* setEPS_COEFF */

/***********************************************************************/
void CglParam::setMAX_SUPPORT(const int max_s)
{
  if (max_s > 0)
    MAX_SUPPORT = max_s;
} /* setMAX_SUPPORT */

/***********************************************************************/
CglParam::CglParam(const double inf, const double eps, const double eps_c,
  const int max_s)
  : INFINIT(inf)
  , EPS(eps)
  , EPS_COEFF(eps_c)
  , MAX_SUPPORT(max_s)
{
}

/***********************************************************************/
CglParam::CglParam(const CglParam &source)
  : INFINIT(source.INFINIT)
  , EPS(source.EPS)
  , EPS_COEFF(source.EPS_COEFF)
  , MAX_SUPPORT(source.MAX_SUPPORT)
{
}

/***********************************************************************/
CglParam *CglParam::clone() const
{
  return new CglParam(*this);
}

/***********************************************************************/
CglParam &CglParam::operator=(const CglParam &rhs)
{
  if (this != &rhs) {
    INFINIT = rhs.INFINIT;
    EPS = rhs.EPS;
    EPS_COEFF = rhs.EPS_COEFF;
    MAX_SUPPORT = rhs.MAX_SUPPORT;
  }
  return *this;
}

/***********************************************************************/
CglParam::~CglParam()
{
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
