// Authors: Matthew Saltzman and Ted Ralphs
// Copyright 2015, Matthew Saltzman and Ted Ralphs
// Licensed under the Eclipse Public License 

#ifndef CoinRational_H
#define CoinRational_H

#include <cmath>

#include "CoinUtilsConfig.h"

#if defined(COINUTILS_HAS_CSTDINT)
#include <cstdint>
#elif defined(COINUTILS_HAS_STDINT_H)
#include <stdint.h>
#endif

//Small class for rational numbers
class COINUTILSLIB_EXPORT CoinRational
{

public:
  int64_t getDenominator() { return denominator_; }
  int64_t getNumerator() { return numerator_; }

  CoinRational()
    : numerator_(0)
    , denominator_(1) {};

  CoinRational(int64_t n, int64_t d)
    : numerator_(n)
    , denominator_(d) {};

  CoinRational(double val, double maxdelta, int64_t maxdnom)
  {
    if (!nearestRational_(val, maxdelta, maxdnom)) {
      numerator_ = 0;
      denominator_ = 1;
    }
  };

private:
  int64_t numerator_;
  int64_t denominator_;

  bool nearestRational_(double val, double maxdelta, int64_t maxdnom);
};

//#############################################################################
/** A function that tests the methods in the class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void CoinRationalUnitTest();

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
