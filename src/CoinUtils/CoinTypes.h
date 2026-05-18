/* Copyright (C) 2004, International Business Machines
 * Corporation and others.  All Rights Reserved.
 * This code is licensed under the terms of the Eclipse Public License (EPL).
 */

#ifndef _CoinTypes_h
#define _CoinTypes_h

#include "CoinUtilsConfig.h"

/* On some systems, we require cstdint/stdint.h to have the 64bit integer type defined. */
#if defined(COINUTILS_HAS_CSTDINT) && defined(__cplusplus)
#include <cstdint>
#elif defined(COINUTILS_HAS_STDINT_H)
#include <stdint.h>
#endif

#define CoinInt64 COINUTILS_INT64_T
#define CoinUInt64 COINUTILS_UINT64_T
#define CoinIntPtr COINUTILS_INTPTR_T

/* ============================================================================= */

typedef COINUTILS_BIGINDEX_T CoinBigIndex;
/* CoinByteArray should be an integer that is long enough for any size of array in bytes */
typedef CoinIntPtr CoinByteArray;

/* ============================================================================= */
#ifndef COIN_BIG_DOUBLE
#define COIN_BIG_DOUBLE 0
#endif

/* See if we want the ability to have long double work arrays */
#if COIN_BIG_DOUBLE == 2
#undef COIN_BIG_DOUBLE
#define COIN_BIG_DOUBLE 0
#define COIN_LONG_WORK 1
typedef long double CoinWorkDouble;
#elif COIN_BIG_DOUBLE == 3
#undef COIN_BIG_DOUBLE
#define COIN_BIG_DOUBLE 1
#define COIN_LONG_WORK 1
typedef long double CoinWorkDouble;
#else
#define COIN_LONG_WORK 0
typedef double CoinWorkDouble;
#endif
/** For factorizations inheriting from CoinDenseFactorization -
    leave partial conversions but back to double. */
typedef double CoinFactorizationDouble2;
#if COIN_BIG_DOUBLE == 0
typedef double CoinFactorizationDouble;
#elif COIN_BIG_DOUBLE == 1
typedef long double CoinFactorizationDouble;
#else
typedef double CoinFactorizationDouble;
#endif

#endif
