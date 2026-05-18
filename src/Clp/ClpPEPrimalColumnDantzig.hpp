/*
   Authors

   Jeremy Omer, Mehdi Towhidi

   Last update: april 10, 2015

 */

#ifndef ClpPEPrimalColumnDantzig_H
#define ClpPEPrimalColumnDantzig_H

#include "ClpPrimalColumnDantzig.hpp"
#include "ClpPrimalColumnPivot.hpp"
#include "ClpSimplex.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpFactorization.hpp"
#include "ClpNonLinearCost.hpp"
#include "ClpSimplexPrimal.hpp"
#include "ClpPackedMatrix.hpp"
#include "ClpPESimplex.hpp"

class CLPLIB_EXPORT ClpPEPrimalColumnDantzig : public ClpPrimalColumnDantzig {

public:
  /** constructors */
  ClpPEPrimalColumnDantzig(double psi = 0.5);
  ClpPEPrimalColumnDantzig(const ClpPEPrimalColumnDantzig &); //copy constructor

  /** destructor */
  ~ClpPEPrimalColumnDantzig();

  /** assignment operator */
  ClpPEPrimalColumnDantzig &operator=(const ClpPEPrimalColumnDantzig &rhs);

  /** clone */
  ClpPrimalColumnPivot *clone(bool copyData = true) const;

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
  int coDegenCompatibles_;
  int coConsecutiveCompatibles_;
  bool updateCompatibles_;
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
