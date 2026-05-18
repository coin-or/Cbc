// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
/*
   Authors

   Jeremy Omer, Mehdi Towhidi

   Last update: april 10, 2015

 */

#ifndef ClpPEPrimalColumnSteepest_H
#define ClpPEPrimalColumnSteepest_H

#include "ClpPrimalColumnSteepest.hpp"
#include "ClpFactorization.hpp"
#include "ClpPESimplex.hpp"
#include <bitset>

//#############################################################################
class CoinIndexedVector;

/** Primal Column Pivot Steepest Edge Algorithm Class

See Forrest-Goldfarb paper for algorithm

*/

class CLPLIB_EXPORT ClpPEPrimalColumnSteepest : public ClpPrimalColumnSteepest {
public:
  ///@name Constructors and destructors
  //@{
  /** Default Constructor
         0 is exact devex, 1 full steepest, 2 is partial exact devex
         3 switches between 0 and 2 depending on factorization
         4 starts as partial dantzig/devex but then may switch between 0 and 2.
         By partial exact devex is meant that the weights are updated as normal
         but only part of the nonbasic variables are scanned.
         This can be faster on very easy problems.
     */
  ClpPEPrimalColumnSteepest(double psi = 0.5, int mode = 3);

  /// Copy constructor
  ClpPEPrimalColumnSteepest(const ClpPEPrimalColumnSteepest &rhs);

  /// Assignment operator
  ClpPEPrimalColumnSteepest &operator=(const ClpPEPrimalColumnSteepest &rhs);

  /// Destructor
  virtual ~ClpPEPrimalColumnSteepest();

  /// Clone
  virtual ClpPrimalColumnPivot *clone(bool copyData = true) const;

public:
  ///@name Algorithmic methods
  //@{

  /** Returns pivot column, -1 if none.
         The Packed CoinIndexedVector updates has cost updates - for normal LP
         that is just +-weight where a feasibility changed.  It also has
         reduced cost from last iteration in pivot row
         Parts of operation split out into separate functions for
         profiling and speed
     */
  virtual int pivotColumn(CoinIndexedVector *updates,
    CoinIndexedVector *spareRow1,
    CoinIndexedVector *spareRow2,
    CoinIndexedVector *spareColumn1,
    CoinIndexedVector *spareColumn2);

  //@}
  /** Save weights - this may initialize weights as well
	 This is as parent but may initialize ClpPESimplex
     */
  virtual void saveWeights(ClpSimplex *model, int mode);
  /// Updates weights - as ordinary but checks for zero moves
  virtual void updateWeights(CoinIndexedVector *input);
  //---------------------------------------------------------------------------
  // Psi
  inline double psi() const
  {
    return psi_;
  }

private:
  /* this PESimplex object is used to identify the compatible variables */
  ClpPESimplex *modelPE_;

  /* psi is the factor used in the bi-dimensional pricing, it is < 1 and
       1/psi grows with the priority given to compatible variables */
  double psi_;

  /* useful counters for the update of the set of compatible variables */
  int iCurrent_;
  int iInterval_;

  /* record if previous iterations concluded that compatibles should not be checked */
  int coDegenCompatibles_;
  int coConsecutiveCompatibles_;
  bool updateCompatibles_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
