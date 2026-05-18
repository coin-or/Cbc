// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
/*
   Authors

   Jeremy Omer

   Last update: april 10, 2015

 */

#ifndef ClpPEDualRowSteepest_H
#define ClpPEDualRowSteepest_H

#include "ClpDualRowSteepest.hpp"
#include "ClpPESimplex.hpp"
class CoinIndexedVector;

//#############################################################################

/** Dual Row Pivot Steepest Edge Algorithm Class

See Forrest-Goldfarb paper for algorithm

*/

class CLPLIB_EXPORT ClpPEDualRowSteepest : public ClpDualRowSteepest {

public:
  /** Default Constructor
         mode: 0 is uninitialized, 1 full, 2 is partial uninitialized,
         3 starts as 2 but may switch to 1.
         By partial is meant that the weights are updated as normal
         but only part of the infeasible basic variables are scanned.
         This can be faster on very easy problems.
     */
  ClpPEDualRowSteepest(double psi = 0.5, int mode = 3);

  /// Copy constructor
  ClpPEDualRowSteepest(const ClpPEDualRowSteepest &);

  /// Assignment operator
  ClpPEDualRowSteepest &operator=(const ClpPEDualRowSteepest &rhs);

  /// Destructor
  virtual ~ClpPEDualRowSteepest();

  /// Clone
  virtual ClpDualRowPivot *clone(bool copyData = true) const;

public:
  ///@name Algorithmic methods
  //@{

  /// Returns pivot row, -1 if none
  virtual int pivotRow();

  /** Save weights - this may initialize weights as well
	 This is as parent but may initialize ClpPESimplex
     */
  virtual void saveWeights(ClpSimplex *model, int mode);
  /** Updates primal solution (and maybe list of candidates)
         Uses input vector which it deletes
         Computes change in objective function
	 As ordinary steepest but checks for zero moves
     */
  virtual void updatePrimalSolution(CoinIndexedVector *input,
    double theta,
    double &changeInObjective);
  //@}

  // Psi
  inline double psi() const
  {
    return psi_;
  }

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
