/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcBranchBase_H
#define CbcBranchBase_H

#include <string>
#include <vector>
#include "OsiBranchingObject.hpp"

enum CbcRangeCompare {
  CbcRangeSame,
  CbcRangeDisjoint,
  CbcRangeSubset,
  CbcRangeSuperset,
  CbcRangeOverlap
};

#include "CbcObject.hpp"
#include "CbcBranchingObject.hpp"
#include "CbcBranchDecision.hpp"
#include "CbcConsequence.hpp"
#include "CbcObjectUpdateData.hpp"

//##############################################################################

/** Compare two ranges. The two bounds arrays are both of size two and
    describe closed intervals. Return the appropriate CbcRangeCompare value
    (first argument being the sub/superset if that's the case). In case of
    overlap (and if \c replaceIfOverlap is true) replace the content of thisBd
    with the intersection of the ranges.
*/
static inline CbcRangeCompare
CbcCompareRanges(double *thisBd, const double *otherBd,
  const bool replaceIfOverlap)
{
  const double lbDiff = thisBd[0] - otherBd[0];
  if (lbDiff < 0) { // lb of this < lb of other
    if (thisBd[1] >= otherBd[1]) { // ub of this >= ub of other
      return CbcRangeSuperset;
    } else if (thisBd[1] < otherBd[0]) {
      return CbcRangeDisjoint;
    } else {
      // overlap
      if (replaceIfOverlap) {
        thisBd[0] = otherBd[0];
      }
      return CbcRangeOverlap;
    }
  } else if (lbDiff > 0) { // lb of this > lb of other
    if (thisBd[1] <= otherBd[1]) { // ub of this <= ub of other
      return CbcRangeSubset;
    } else if (thisBd[0] > otherBd[1]) {
      return CbcRangeDisjoint;
    } else {
      // overlap
      if (replaceIfOverlap) {
        thisBd[1] = otherBd[1];
      }
      return CbcRangeOverlap;
    }
  } else { // lb of this == lb of other
    if (thisBd[1] == otherBd[1]) {
      return CbcRangeSame;
    }
    return thisBd[1] < otherBd[1] ? CbcRangeSubset : CbcRangeSuperset;
  }

  return CbcRangeSame; // fake return
}

//#############################################################################

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
