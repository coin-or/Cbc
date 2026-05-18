/* -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; -*-
 *
 * This file is part of the COIN-OR CoinUtils package.
 *
 * @file   CoinTerm.hpp
 * @brief  Lightweight structure holding a column index and coefficient.
 *
 */

#ifndef COIN_TERM_HPP
#define COIN_TERM_HPP
#include "CoinUtilsConfig.h"

/**
 * @brief Simple POD used to store a matrix column index and its coefficient.
 */
struct COINUTILSLIB_EXPORT CoinTerm {
  int index;
  double value;

  CoinTerm()
    : index(0)
    , value(0.0)
  {
  }

  CoinTerm(int idx, double val)
    : index(idx)
    , value(val)
  {
  }

  struct COINUTILSLIB_EXPORT AscendingValueThenIndexComparator {
    bool operator()(const CoinTerm &left, const CoinTerm &right) const {
      if (left.value == right.value) {
        return left.index < right.index;
      }
      return left.value < right.value;
    }
  };

  static constexpr AscendingValueThenIndexComparator AscendingValueThenIndex = AscendingValueThenIndexComparator();
};

#endif // COIN_TERM_HPP
