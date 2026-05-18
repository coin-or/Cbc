// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CoinDistance_H
#define CoinDistance_H

#include <iterator>

//-------------------------------------------------------------------
//
// Attempt to provide an std::distance function
// that will work on multiple platforms
//
//-------------------------------------------------------------------

/** CoinDistance

This is the Coin implementation of the std::function that is 
designed to work on multiple platforms.
*/
template < class ForwardIterator, class Distance >
void coinDistance(ForwardIterator first, ForwardIterator last,
  Distance &n)
{
  n = std::distance(first, last);
}

template < class ForwardIterator >
size_t coinDistance(ForwardIterator first, ForwardIterator last)
{
  return std::distance(first, last);
}

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
