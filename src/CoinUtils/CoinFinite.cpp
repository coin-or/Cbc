// Copyright (C) 2011, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinFinite.hpp"
#include "CoinUtilsConfig.h"

#include <cfloat>
#include <cmath>

#ifdef HAVE_CIEEEFP
#include <cieeefp>
#else
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif
#endif

bool CoinFinite(double val)
{
#ifdef COINUTILS_C_FINITE
  return COINUTILS_C_FINITE(val) != 0;
#else
  return val != COIN_DBL_MAX && val != -COIN_DBL_MAX;
#endif
}

bool CoinIsnan(double val)
{
#ifdef COINUTILS_C_ISNAN
  return COINUTILS_C_ISNAN(val) != 0;
#else
  return false;
#endif
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
