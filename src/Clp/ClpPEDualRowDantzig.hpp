// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
/*
   Authors

   Jeremy Omer

   Last update: april 10, 2015

 */

#ifndef ClpPEDualRowDantzig_H
#define ClpPEDualRowDantzig_H

#include "ClpDualRowPivot.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpSimplex.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpFactorization.hpp"
#include "ClpNonLinearCost.hpp"
#include "ClpSimplexDual.hpp"
#include "ClpPackedMatrix.hpp"
#include "ClpPESimplex.hpp"

class CLPLIB_EXPORT ClpPEDualRowDantzig : public ClpDualRowDantzig {

public:
  /// Default Constructor
  ClpPEDualRowDantzig(double psi = 0.5);

  /// Copy constructor
  ClpPEDualRowDantzig(const ClpPEDualRowDantzig &);

  /// Assignment operator
  ClpPEDualRowDantzig &operator=(const ClpPEDualRowDantzig &rhs);

  /// Destructor
  virtual ~ClpPEDualRowDantzig();

  /// Clone
  virtual ClpDualRowPivot *clone(bool copyData = true) const;

public:
  ///@name Algorithmic methods
  //@{

  /// Returns pivot row, -1 if none
  virtual int pivotRow();

  /// Update the compatible variables and
  /// call the base class method to update weights
  virtual double updateWeights(CoinIndexedVector *input,
    CoinIndexedVector *spare,
    CoinIndexedVector *spare2,
    CoinIndexedVector *updatedColumn);

  /** Save weights - this may initialize weights as well
	 This is as parent but may initialize ClpPESimplex
     */
  virtual void saveWeights(ClpSimplex *model, int mode);
  //@}

  //---------------------------------------------------------------------------

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
  bool updateCompatibles_;
  int coDegenCompatibles_, coConsecutiveCompatibles_;
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
