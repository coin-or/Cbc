// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <stdio.h>
#include <math.h>

#include "CoinPresolveMatrix.hpp"
#include "CoinPresolveFixed.hpp"
#include "CoinPresolveDual.hpp"
#include "CoinMessage.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinFloatEqual.hpp"

/*
  Define PRESOLVE_DEBUG and PRESOLVE_CONSISTENCY as compile flags! If not
  uniformly on across all uses of presolve code you'll get something between
  garbage and a core dump. See comments in CoinPresolvePsdebug.hpp
*/
#if PRESOLVE_DEBUG > 0 || PRESOLVE_CONSISTENCY > 0
#include "CoinPresolvePsdebug.hpp"
#endif

/*
  Set to > 0 to enable bound propagation on dual variables? This has been
  disabled for a good while. See notes in code.

  #define PRESOLVE_TIGHTEN_DUALS 1
*/
  #define PRESOLVE_TIGHTEN_DUALS 1
/*
  Guards incorrect code that attempts to adjust a column solution. See
  comments with the code. Could possibly be fixed, which is why it hasn't
  been chopped out.

  #define REMOVE_DUAL_ACTION_REPAIR_SOLN 1
*/

/*
  In this transform we're looking to prove bounds on the duals y<i> and
  reduced costs cbar<j>. We can use this in two ways.
  
  First, if we can prove a bound on the reduced cost cbar<j> with strict
  inequality,
    * cbar<j> > 0 ==> x<j> NBLB at optimality
    * cbar<j> < 0 ==> x<j> NBUB at optimality
  If the required bound is not finite, the problem is unbounded.  Andersen &
  Andersen call this a dominated column.

  Second, suppose we can show that cbar<i> = -y<i> is strictly nonzero
  for the logical s<i> associated with some inequality i. (There must be
  exactly one finite row bound, so that the logical has exactly one finite
  bound which is 0).  If cbar<i> demands the logical be nonbasic at bound,
  it's zero, and we can convert the inequality to an equality.

  Not based on duals and reduced costs, but in the same vein, if we can
  identify an architectural x<j> that'll accomplish the same goal (bring
  one constraint tight when moved in a favourable direction), we can make
  that constraint an equality. The conditions are different, however:
    * Moving x<j> in the favourable direction tightens exactly one
      constraint.
    * The bound (l<j> or u<j>) in the favourable direction is infinite.
    * If x<j> is entangled in any other constraints, moving x<j> in the
      favourable direction loosens those constraints.
  A bit of thought and linear algebra is required to convince yourself that in
  the circumstances, cbar<j> = c<j>, hence we need only look at c<j> to
  determine the favourable direction.

  
  To get from bounds on the duals to bounds on the reduced costs, note that
    cbar<j> = c<j> - (SUM{P}a<ij>y<i> + SUM{M}a<ij>y<i>)
  for P = {a<ij> > 0} and M = {a<ij> < 0}. Then
    cbarmin<j> = c<j> - (SUM{P}a<ij>ymax<i> + SUM{M}a<ij>ymin<i>)
    cbarmax<j> = c<j> - (SUM{P}a<ij>ymin<i> + SUM{M}a<ij>ymax<i>)

  As A&A note, the reverse implication also holds:
    * if l<j> = -infty, cbar<j> <= 0 at optimality
    * if u<j> =  infty, cbar<j> >= 0 at optimality

  We can use this to run bound propagation on the duals in an attempt to
  tighten bounds and force more reduced costs to strict inequality. Suppose
  u<j> = infty. Then cbar<j> >= 0. It must be possible to achieve this, so
  cbarmax<j> >= 0:
    0 <= c<j> - (SUM{P}a<ij>ymin<i> + SUM{M}a<ij>ymax<i>)
  Solve for y<t>, a<tj> > 0:
    y<t> <= 1/a<tj>[c<j> - (SUM{P\t}a<ij>ymin<i> + SUM{M}a<ij>ymax<i>)]
  If a<tj> < 0,
    y<t> >= 1/a<tj>[c<j> - (SUM{P}a<ij>ymin<i> + SUM{M\t}a<ij>ymax<i>)]
  For l<j> = -infty, cbar<j> <= 0, hence cbarmin <= 0: 
    0 >= c<j> - (SUM{P}a<ij>ymax<i> + SUM{M}a<ij>ymin<i>) 
  Solve for y<t>, a<tj> > 0:
    y<t> >= 1/a<tj>[c<j> - (SUM{P\t}a<ij>ymax<i> + SUM{M}a<ij>ymin<i>)]
  If a<tj> < 0,
    y<t> <= 1/a<tj>[c<j> - (SUM{P}a<ij>ymax<i> + SUM{M\t}a<ij>ymin<i>)]

  We can get initial bounds on ymin<i> and ymax<i> from column singletons
  x<j>, where
    cbar<j> = c<j> - a<ij>y<i>
  If u<j> = infty, then at optimality
    0 <= cbar<j> = c<j> - a<ij>y<i>
  For a<ij> > 0 we have
    y<i> <= c<j>/a<ij>
  and for a<ij> < 0 we have
    y<i> >= c<j>/a<ij>
  We can do a similar calculation for l<j> = -infty, so that
    0 >= cbar<j> = c<j> - a<ij>y<i>
  For a<ij> > 0 we have
    y<i> >= c<j>/a<ij>
  and for a<ij> < 0 we have
    y<i> <= c<j>/a<ij>
  A logical is a column singleton with c<j> = 0.0 and |a<ij>| = 1.0. One or
  both bounds can be finite, depending on constraint representation.

  This code is hardwired for minimisation.

  Extensive rework of the second part of the routine Fall 2011 to add a
  postsolve transform to fix bug #67, incorrect status for logicals.
  -- lh, 111208 --
*/
/*
  Original comments from 040916:

  This routine looks to be something of a work in progress.  Down in the
  bound propagation loop, why do we only work with variables with u_j =
  infty? The corresponding section of code for l_j = -infty is ifdef'd away.
  l<j> = -infty is uncommon; perhaps it's not worth the effort?

  Why exclude the code protected by PRESOLVE_TIGHTEN_DUALS? Why are we
  using ekkinf instead of PRESOLVE_INF?
*/
const CoinPresolveAction
  *
  remove_dual_action::presolve(CoinPresolveMatrix *prob,
    const CoinPresolveAction *next)
{
#if PRESOLVE_DEBUG > 0
  std::cout
    << "Entering remove_dual_action::presolve, " << prob->nrows_
    << "x" << prob->ncols_ << "." << std::endl;
#endif
#if PRESOLVE_DEBUG > 0 || PRESOLVE_CONSISTENCY > 0
  presolve_check_sol(prob, 2, 1, 1);
  presolve_check_nbasic(prob);
#endif
#if PRESOLVE_DEBUG > 0 || COIN_PRESOLVE_TUNING > 0
  int startEmptyRows = 0;
  int startEmptyColumns = 0;
  startEmptyRows = prob->countEmptyRows();
  startEmptyColumns = prob->countEmptyCols();
#if COIN_PRESOLVE_TUNING > 0
  double startTime = 0.0;
  if (prob->tuning_)
    startTime = CoinCpuTime();
#endif
#endif

  // column-major representation
  const int ncols = prob->ncols_;
  const CoinBigIndex *const mcstrt = prob->mcstrt_;
  const int *const hincol = prob->hincol_;
  const int *const hrow = prob->hrow_;
  const double *const colels = prob->colels_;
  const double *const cost = prob->cost_;

  // column type, bounds, solution, and status
  const unsigned char *const integerType = prob->integerType_;
  const double *const clo = prob->clo_;
  const double *const cup = prob->cup_;
  double *const csol = prob->sol_;
  unsigned char *&colstat = prob->colstat_;

  // row-major representation
  const int nrows = prob->nrows_;
  const CoinBigIndex *const mrstrt = prob->mrstrt_;
  const int *const hinrow = prob->hinrow_;
  const int *const hcol = prob->hcol_;
#if 1 //REMOVE_DUAL_ACTION_REPAIR_SOLN > 0
  const double *const rowels = prob->rowels_;
#endif

  // row bounds
  double *const rlo = prob->rlo_;
  double *const rup = prob->rup_;

  // tolerances
  const double ekkinf2 = PRESOLVE_SMALL_INF;
  const double ekkinf = ekkinf2 * 1.0e8;
  const double ztolcbarj = prob->ztoldj_;
  const CoinRelFltEq relEq(prob->ztolzb_);

  /*
  Grab one of the preallocated scratch arrays to hold min and max values of
  row duals.
*/
  double *ymin = prob->usefulRowDouble_;
  double *ymax = ymin + nrows;
  /* use arrays to say which bound active
     1 - lower
     2 - upper
     4 - fixed
     8 - don't touch
  */
  char *active = reinterpret_cast< char * >(prob->usefulColumnInt_+ncols);
  memset(active, 0, ncols);
  //int nOneBound = 0;
  int numberLook = prob->numberColsToDo_;
  int *look = prob->colsToDo_;
  for (int iLook = 0; iLook < numberLook; iLook++) {
    int j = look[iLook];
    char type;
    if (cup[j] >= ekkinf) {
      if (clo[j] <= -ekkinf)
        type = 0;
      else
        type = 1;
    } else {
      if (clo[j] <= -ekkinf)
        type = 2;
      else if (clo[j] < cup[j])
        type = 3;
      else
        type = 4;
    }
    //if (type == 1 || type == 2)
    //  nOneBound++;
    active[j] = type;
  }
  //int nFreed = 0;
#define USE_ACTIVE 1
  //#define PRESOLVE_DEBUG 2
  for (int i = 0; i < nrows; i++) {
    const CoinBigIndex krs = mrstrt[i];
    const CoinBigIndex kre = krs + hinrow[i];
    int positive = 0;
    int negative = 0;
    CoinBigIndex onlyPositive = -1;
    CoinBigIndex onlyNegative = -1;
    bool allSingletons = (prob->presolveOptions_&0x200000)!=0;
    for (CoinBigIndex k = krs; k < kre; k++) {
      const double coeff = rowels[k];
      const int icol = hcol[k];
      char type = active[icol];
      //const double lb = clo[icol] ;
      //const double ub = cup[icol] ;
      if (hincol[icol]!=1)
	allSingletons = false;
      if (type == 0 || (type & 8) != 0) {
        // free or already used
        positive = 2;
        negative = 2;
        break;
      } else if (type == 1) {
        // lower bound active
        if (coeff > 0.0) {
          positive++;
        } else {
          negative++;
        }
      } else if (type == 2) {
        // upper bound active
        if (coeff < 0.0) {
          positive++;
        } else {
          negative++;
        }
      } else if (type == 3) {
        // both bounds active
        if (coeff > 0.0) {
          onlyPositive = k;
          positive++;
        } else {
          onlyNegative = k;
          negative++;
        }
      } else {
        // fixed
      }
    }
    bool doneSomething = false;
    if (onlyPositive >= 0 && positive == 1) {
      const double coeff = rowels[onlyPositive];
      const int icol = hcol[onlyPositive];
      const double lb = clo[icol];
      const double ub = cup[icol];
      // Get possible without icol
      double maxUp = 0.0;
      double maxDown = 0.0;
#if USE_ACTIVE > 2
      printf("%d pos coeff %g bounds %g,%g rhs %g,%g - %d negative",
        icol, coeff, lb, ub, rlo[i], rup[i], negative);
      printf(", row is %g <= ", rlo[i]);
#endif
      for (CoinBigIndex k = krs; k < kre; k++) {
        const int icol2 = hcol[k];
#if USE_ACTIVE > 2
        printf("(%d,%g - bounds %g,%g) ", icol2, rowels[k],
          clo[icol2], cup[icol2]);
#endif
        if (icol2 == icol)
          continue;
        const double coeff2 = rowels[k];
        const double lb2 = clo[icol2];
        const double ub2 = cup[icol2];
        if (coeff2 > 0.0) {
          if (ub2 > 1.0e30)
            maxUp = COIN_DBL_MAX;
          else
            maxUp += coeff2 * ub2;
          if (lb2 < -1.0e30)
            maxDown = -COIN_DBL_MAX;
          else
            maxDown += coeff2 * lb2;
        } else {
          if (lb2 < -1.0e30)
            maxUp = COIN_DBL_MAX;
          else
            maxUp += coeff2 * lb2;
          if (ub2 > 1.0e30)
            maxDown = -COIN_DBL_MAX;
          else
            maxDown += coeff2 * ub2;
        }
      }
      double impliedLower = (rlo[i] - maxUp) / coeff;
      double impliedUpper = (rup[i] - maxDown) / coeff;
#if USE_ACTIVE > 2
      printf("<= %g - implied l,u %g,%g\n", rup[i], impliedLower, impliedUpper);
#endif
      if (impliedLower > lb - 1.0e-7) {
#if USE_ACTIVE == 2
        printf("%d pos coeff %g bounds %g,%g maxDown %g rhs %g,%g - %d negative\n",
          icol, coeff, lb, ub, maxDown, rlo[i], rup[i], negative);
#endif
#if USE_ACTIVE > 1
        printf("can take off lb as implied lower %g > %g\n", impliedLower, lb);
#endif
        doneSomething = true;
        active[icol] = 2 + 8;
        //nFreed++;
#define TRY_UPPER 3
#if TRY_UPPER > 1
      } else if (impliedUpper < ub + 1.0e-7) {
#if TRY_UPPER > 2
        doneSomething = true;
        active[icol] = 1 + 8;
        //nFreed++;
#endif
#endif
      } else {
      }
    } else if (allSingletons && onlyPositive >= 0) {
      //printf("Possible pos\n");
      // clean up and simplify later
      // for now - only if all 1.0 and integer stuff ok - or not integer
      double dlo = rlo[i];
      bool possible = true;
      int nIntegers=0;
      for (CoinBigIndex k = krs; k < kre; k++) {
	const int icol = hcol[k];
	double value = rowels[k];
	if (integerType[icol]) {
	  nIntegers++;
	  if (value!=1.0)
	    possible = false;
	}
	if (clo[icol])
	  dlo -= clo[icol]*value;
      }
      if (nIntegers&&nIntegers<hinrow[i])
	possible = false;
      if (possible) {
	while (dlo>1.0e-9) {
	  // get cheapest
	  double cheapest = 1.0e50;
	  int iCheap = -1;
	  for (CoinBigIndex k = krs; k < kre; k++) {
	    const int icol = hcol[k];
	    if (cup[icol]>clo[icol]) {
	      double value = rowels[k];
	      double thisCost = cost[icol];
	      if (thisCost<value*cheapest) {
		cheapest = thisCost/value;
		iCheap = k;
	      }
	    }
	  }
	  if (iCheap<0)
	    break;
	  double value = rowels[iCheap];
	  iCheap = hcol[iCheap]; 
	  double diff2 = dlo/value;
	  double diff = std::min(cup[iCheap]-clo[iCheap],diff2);
	  // fix
	  prob->clo_[iCheap] += diff;
	  prob->cup_[iCheap] = prob->clo_[iCheap];
	  dlo -= diff*value;
	}
	// fix rest to lower bound
	for (CoinBigIndex k = krs; k < kre; k++) {
	  const int icol = hcol[k];
	  if (cup[icol]>clo[icol]) {
	    prob->cup_[icol] = clo[icol];
	    csol[icol] = clo[icol];
	    colstat[icol] = CoinPrePostsolveMatrix::atLowerBound;
	  }
	}
      }
    } 
    if (onlyNegative >= 0 && negative == 1 && !doneSomething) {
      const double coeff = rowels[onlyNegative];
      const int icol = hcol[onlyNegative];
      const double lb = clo[icol];
      const double ub = cup[icol];
      // Get possible without icol
      double maxUp = 0.0;
      double maxDown = 0.0;
#if USE_ACTIVE > 2
      printf("%d neg coeff %g bounds %g,%g rhs %g,%g - %d positive",
        icol, coeff, lb, ub, rlo[i], rup[i], positive);
      printf(", row is %g <= ", rlo[i]);
#endif
      for (CoinBigIndex k = krs; k < kre; k++) {
        const int icol2 = hcol[k];
#if USE_ACTIVE > 2
        printf("(%d,%g - bounds %g,%g) ", icol2, rowels[k],
          clo[icol2], cup[icol2]);
#endif
        if (icol2 == icol)
          continue;
        const double coeff2 = rowels[k];
        const double lb2 = clo[icol2];
        const double ub2 = cup[icol2];
        if (coeff2 > 0.0) {
          if (ub2 > 1.0e30)
            maxUp = COIN_DBL_MAX;
          else
            maxUp += coeff2 * ub2;
          if (lb2 < -1.0e30)
            maxDown = -COIN_DBL_MAX;
          else
            maxDown += coeff2 * lb2;
        } else {
          if (lb2 < -1.0e30)
            maxUp = COIN_DBL_MAX;
          else
            maxUp += coeff2 * lb2;
          if (ub2 > 1.0e30)
            maxDown = -COIN_DBL_MAX;
          else
            maxDown += coeff2 * ub2;
        }
      }
      double impliedLower = (rup[i] - maxDown) / coeff;
      double impliedUpper = (rlo[i] - maxUp) / coeff;
#if USE_ACTIVE > 2
      printf("<= %g - implied l,u %g,%g\n", rup[i], impliedLower, impliedUpper);
#endif
      int nOne = 0;
      for (CoinBigIndex k = krs; k < kre; k++) {
        const double coeff = rowels[k];
        const int icol = hcol[k];
        if (hincol[icol] == 1 && coeff > 0.0)
          nOne++;
      }
      if (nOne >= hinrow[i] - 1) {
        for (CoinBigIndex k = krs; k < kre; k++) {
          const int icol = hcol[k];
          const double coeff = rowels[k];
          if (hincol[icol] == 1 && coeff > 0.0) {
#if USE_ACTIVE > 2
            printf("(%d has one entry) ", icol);
#endif
            if (active[icol] == 3 || clo[icol] < -1.0e30)
              nOne = 0; // no good
          }
        }
#if USE_ACTIVE > 2
        printf("\n");
#endif
      }
      if (impliedLower > lb - 1.0e-7 /*&& nOne*/) {
#if USE_ACTIVE > 1
        printf("second type can take off lb of %g on column %d as implied lower %g - effective upper %g, coeff %g on row %d\n", lb, icol, impliedLower, maxUp, coeff, i);
#endif
        active[icol] = 2 + 8;
        //nFreed++;
#if TRY_UPPER > 1
      } else if (impliedUpper < ub + 1.0e-7) {
#if TRY_UPPER > 2
        doneSomething = true;
        active[icol] = 1 + 8;
        //nFreed++;
#endif
#endif
      }
    } else if (allSingletons  && !doneSomething && 
	       onlyNegative >= 0 && negative >= 1) {
      //printf("Possible neg\n");
    }
  }
  //printf("%d had only one bound - %d added\n",nOneBound,nFreed);
  /*
  Initialise row dual min/max. The defaults are +/- infty, but if we know
  that the logical has no upper (lower) bound, it can only be nonbasic at
  the other bound. Given minimisation, we can conclude:
    * <= constraint ==> [0,infty] ==> NBLB ==> cbar<i> > 0 ==> y<i> < 0
    * >= constraint ==> [-infty,0] ==> NBUB ==> cbar<i> < 0 ==> y<i> > 0
  Range constraints are not helpful here, the dual can be either sign. There's
  no calculation because we assume the objective coefficient c<i> = 0  and the
  coefficient a<ii> = 1.0, hence cbar<i> = -y<i>.
*/
  for (int i = 0; i < nrows; i++) {
    const bool no_lb = (rup[i] >= ekkinf);
    const bool no_ub = (rlo[i] <= -ekkinf);

    ymin[i] = ((no_lb && !no_ub) ? 0.0 : (-PRESOLVE_INF));
    ymax[i] = ((no_ub && !no_lb) ? 0.0 : PRESOLVE_INF);
  }
  /*
  We can do a similar calculation with singleton columns where the variable
  has only one finite bound, but we have to work a bit harder, as the cost
  coefficient c<j> is not necessarily 0.0 and a<ij> is not necessarily 1.0.

  cbar<j> = c<j> - y<i>a<ij>, hence y<i> = c<j>/a<ij> at cbar<j> = 0. The
  question is whether this is an upper or lower bound on y<i>.
    * x<j> NBLB ==> cbar<j> >= 0
      a<ij> > 0 ==> increasing y<i> decreases cbar<j> ==> upper bound 
      a<ij> < 0 ==> increasing y<i> increases cbar<j> ==> lower bound 
    * x<j> NBUB ==> cbar<j> <= 0
      a<ij> > 0 ==> increasing y<i> decreases cbar<j> ==> lower bound 
      a<ij> < 0 ==> increasing y<i> increases cbar<j> ==> upper bound 
  The condition below (simple test for equality to choose the bound) looks
  a bit odd, but a bit of boolean algebra should convince you it's correct.
  We have a bound; the only question is whether it's upper or lower.

  Skip integer variables; it's far too likely that we'll tighten infinite
  bounds elsewhere in presolve.

  NOTE: If bound propagation is applied to continuous variables, the same
  	hazard will apply.  -- lh, 110611 --
*/
  for (int j = 0; j < ncols; j++) {
    if (integerType[j])
      continue;
    if (hincol[j] != 1)
      continue;
#ifndef USE_ACTIVE
    const bool no_ub = (cup[j] >= ekkinf);
    const bool no_lb = (clo[j] <= -ekkinf);
#else
    const bool no_ub = (active[j] & 6) == 0;
    const bool no_lb = (active[j] & 5) == 0;
#endif
    if (no_ub != no_lb) {
      const int &i = hrow[mcstrt[j]];
      double aij = colels[mcstrt[j]];
      PRESOLVEASSERT(fabs(aij) > ZTOLDP);
      const double yzero = cost[j] / aij;
      if ((aij > 0.0) == no_ub) {
        if (ymax[i] > yzero)
          ymax[i] = yzero;
      } else {
        if (ymin[i] < yzero)
          ymin[i] = yzero;
      }
    }
  }
  int nfixup_cols = 0;
  int nfixdown_cols = ncols;
  // Grab another work array, sized to ncols_
  int *fix_cols = prob->usefulColumnInt_;
#if PRESOLVE_TIGHTEN_DUALS > 0
  double *cbarmin = new double[ncols];
  double *cbarmax = new double[ncols];
#endif
  /*
  Now we have (admittedly weak) bounds on the dual for each row. We can use
  these to calculate upper and lower bounds on cbar<j>. Open loops to
  take multiple passes over all columns.
*/
  int nPass = 0;
  while (nPass++ < 100) {
    int tightened = 0;
    for (int j = 0; j < ncols; j++) {
      if (hincol[j] <= 0)
        continue;
      /*
  Calculate min cbar<j> and max cbar<j> for the column by calculating the
  contribution to c<j>-ya<j> using ymin<i> and ymax<i> from above.
*/
      const CoinBigIndex &kcs = mcstrt[j];
      const CoinBigIndex kce = kcs + hincol[j];
      // Number of infinite rows
      int nflagu = 0;
      int nflagl = 0;
      // Number of ordinary rows
#if PROCESS_INFINITE_LB
      int nordu = 0;
#endif
      int nordl = 0;
      double cbarjmin = cost[j];
      double cbarjmax = cbarjmin;
      for (CoinBigIndex k = kcs; k < kce; k++) {
        const int &i = hrow[k];
        const double &aij = colels[k];
        const double mindelta = aij * ymin[i];
        const double maxdelta = aij * ymax[i];

        if (aij > 0.0) {
          if (ymin[i] >= -ekkinf2) {
            cbarjmax -= mindelta;
#if PROCESS_INFINITE_LB
            nordu++;
#endif
          } else {
            nflagu++;
          }
          if (ymax[i] <= ekkinf2) {
            cbarjmin -= maxdelta;
            nordl++;
          } else {
            nflagl++;
          }
        } else {
          if (ymax[i] <= ekkinf2) {
            cbarjmax -= maxdelta;
#if PROCESS_INFINITE_LB
            nordu++;
#endif
          } else {
            nflagu++;
          }
          if (ymin[i] >= -ekkinf2) {
            cbarjmin -= mindelta;
            nordl++;
          } else {
            nflagl++;
          }
        }
      }
      /*
  See if we can tighten bounds on a dual y<t>. See the comments at the head of
  the file for the linear algebra.

  At this point, I don't understand the restrictions on propagation. Neither
  are necessary. The net effect is to severely restrict the circumstances
  where we'll propagate. Requiring nflagu == 1 excludes the case where all
  duals have finite bounds (unlikely?). And why cbarjmax < -ztolcbarj?

  In any event, we're looking for the row that's contributing the infinity. If
  a<tj> > 0, the contribution is due to ymin<t> and we tighten ymax<t>;
  a<tj> < 0, ymax<t> and ymin<t>, respectively.

  The requirement cbarjmax < (ymax[i]*aij-ztolcbarj) ensures we don't propagate
  tiny changes. If we make a change, it will affect cbarjmin. Make the
  adjustment immediately.

  Continuous variables only.
*/
      if (!integerType[j]) {
#ifndef USE_ACTIVE
        const bool no_cub = (cup[j] >= ekkinf);
        //const bool no_clb = (clo[j] <= -ekkinf) ;
#else
        const bool no_cub = (active[j] & 6) == 0;
        //const bool no_clb = (active[j]&5)==0 ;
#endif
        if (no_cub) {
          if (nflagu == 1 && cbarjmax < -ztolcbarj) {
            for (CoinBigIndex k = kcs; k < kce; k++) {
              const int i = hrow[k];
              const double aij = colels[k];
              if (aij > 0.0 && ymin[i] < -ekkinf2) {
                if (cbarjmax < (ymax[i] * aij - ztolcbarj)) {
                  const double newValue = cbarjmax / aij;
                  if (ymax[i] > ekkinf2 && newValue <= ekkinf2) {
                    nflagl--;
                    cbarjmin -= aij * newValue;
                  } else if (ymax[i] <= ekkinf2) {
                    cbarjmin -= aij * (newValue - ymax[i]);
                  }
#if PRESOLVE_DEBUG > 1
                  std::cout
                    << "NDUAL(infu/inf): u(" << j << ") = " << cup[j]
                    << ": max y(" << i << ") was " << ymax[i]
                    << " now " << newValue << "." << std::endl;
#endif
                  ymax[i] = newValue;
                  tightened++;
                }
              } else if (aij < 0.0 && ymax[i] > ekkinf2) {
                if (cbarjmax < (ymin[i] * aij - ztolcbarj)) {
                  const double newValue = cbarjmax / aij;
                  if (ymin[i] < -ekkinf2 && newValue >= -ekkinf2) {
                    nflagl--;
                    cbarjmin -= aij * newValue;
                  } else if (ymin[i] >= -ekkinf2) {
                    cbarjmin -= aij * (newValue - ymin[i]);
                  }
#if PRESOLVE_DEBUG > 1
                  std::cout
                    << "NDUAL(infu/inf): u(" << j << ") = " << cup[j]
                    << ": min y(" << i << ") was " << ymin[i]
                    << " now " << newValue << "." << std::endl;
#endif
                  ymin[i] = newValue;
                  tightened++;
                  // Huh? asymmetric
                  // cbarjmin = 0.0 ;
                }
              }
            }
          } else if (nflagl == 0 && nordl == 1 && cbarjmin < -ztolcbarj) {
            /*
  This is a column singleton. Why are we doing this? It's not like changes
  to other y will affect this.
*/
            /* In fact a bug as ymax/min changed - but as I don't
	       wish to totally understand coding will let it abort
	       if not a singleton */
            if (kce != kcs + 1)
              abort();
          }
        } // end u<j> = infty
#if PROCESS_INFINITE_LB
        /*
  Unclear why this section is commented out, except for the possibility that
  bounds of -infty < x < something are rare and it likely suffered fromm the
  same errors. Consider the likelihood that this whole block needs to be
  edited to sway min/max, & similar.
*/
        if (no_clb) {
          // cbarj can not be positive
          if (cbarjmin > ztolcbarj && nflagl == 1) {
            // We can make bound finite one way
            for (CoinBigIndex k = kcs; k < kce; k++) {
              const int i = hrow[k];
              const double coeff = colels[k];

              if (coeff < 0.0 && ymin[i] < -ekkinf2) {
                // ymax[i] has upper bound
                if (cbarjmin > ymax[i] * coeff + ztolcbarj) {
                  const double newValue = cbarjmin / coeff;
                  // re-compute hi
                  if (ymax[i] > ekkinf2 && newValue <= ekkinf2) {
                    nflagu--;
                    cbarjmax -= coeff * newValue;
                  } else if (ymax[i] <= ekkinf2) {
                    cbarjmax -= coeff * (newValue - ymax[i]);
                  }
                  ymax[i] = newValue;
                  tightened++;
#if PRESOLVE_DEBUG > 1
                  printf("Col %d, row %d max pi now %g\n", j, i, ymax[i]);
#endif
                }
              } else if (coeff > 0.0 && ymax[i] > ekkinf2) {
                // ymin[i] has lower bound
                if (cbarjmin > ymin[i] * coeff + ztolcbarj) {
                  const double newValue = cbarjmin / coeff;
                  // re-compute lo
                  if (ymin[i] < -ekkinf2 && newValue >= -ekkinf2) {
                    nflagu--;
                    cbarjmax -= coeff * newValue;
                  } else if (ymin[i] >= -ekkinf2) {
                    cbarjmax -= coeff * (newValue - ymin[i]);
                  }
                  ymin[i] = newValue;
                  tightened++;
#if PRESOLVE_DEBUG > 1
                  printf("Col %d, row %d min pi now %g\n", j, i, ymin[i]);
#endif
                }
              }
            }
          } else if (nflagu == 0 && nordu == 1 && cbarjmax > ztolcbarj) {
            // We may be able to tighten
            for (CoinBigIndex k = kcs; k < kce; k++) {
              const int i = hrow[k];
              const double coeff = colels[k];

              if (coeff < 0.0) {
                ymax[i] += cbarjmax / coeff;
                cbarjmax = 0.0;
                tightened++;
#if PRESOLVE_DEBUG > 1
                printf("Col %d, row %d max pi now %g\n", j, i, ymax[i]);
#endif
              } else if (coeff > 0.0) {
                ymin[i] += cbarjmax / coeff;
                cbarjmax = 0.0;
                tightened++;
#if PRESOLVE_DEBUG > 1
                printf("Col %d, row %d min pi now %g\n", j, i, ymin[i]);
#endif
              }
            }
          }
        } // end l<j> < -infty
#endif // PROCESS_INFINITE_LB
      }
/*
  That's the end of propagation of bounds for dual variables.
*/
#if PRESOLVE_TIGHTEN_DUALS > 0
      cbarmin[j] = (nflagl ? (-PRESOLVE_INF) : cbarjmin);
      cbarmax[j] = (nflagu ? PRESOLVE_INF : cbarjmax);
#endif
      /*
  If cbarmin<j> > 0 (strict inequality) then x<j> NBLB at optimality. If l<j>
  is -infinity, notify the user and set the status to unbounded.
*/
      if (cbarjmin > ztolcbarj && nflagl == 0 && !prob->colProhibited2(j)) {
        if (clo[j] <= -ekkinf) {
          CoinMessageHandler *msghdlr = prob->messageHandler();
          msghdlr->message(COIN_PRESOLVE_COLUMNBOUNDB, prob->messages())
            << j << CoinMessageEol;
          prob->status_ |= 2;
          break;
        } else {
          fix_cols[--nfixdown_cols] = j;
#if PRESOLVE_DEBUG > 1
          std::cout << "NDUAL(fix l): fix x(" << j << ")";
          if (csol)
            std::cout << " = " << csol[j];
          std::cout
            << " at l(" << j << ") = " << clo[j] << "; cbar("
            << j << ") > " << cbarjmin << "." << std::endl;
#endif
          if (csol) {
#if REMOVE_DUAL_ACTION_REPAIR_SOLN > 0
            /*
  Original comment: User may have given us feasible solution - move if simple

  Except it's not simple. The net result is that we end up with an excess of
  basic variables. Mark x<j> NBLB and let the client recalculate the solution
  after establishing the basis.
*/
            if (csol[j] - clo[j] > 1.0e-7 && hincol[j] == 1) {
              double value_j = colels[mcstrt[j]];
              double distance_j = csol[j] - clo[j];
              int row = hrow[mcstrt[j]];
              // See if another column can take value
              for (CoinBigIndex kk = mrstrt[row]; kk < mrstrt[row] + hinrow[row]; kk++) {
                const int k = hcol[kk];
                if (colstat[k] == CoinPrePostsolveMatrix::superBasic)
                  continue;

                if (hincol[k] == 1 && k != j) {
                  const double value_k = rowels[kk];
                  double movement;
                  if (value_k * value_j > 0.0) {
                    // k needs to increase
                    double distance_k = cup[k] - csol[k];
                    movement = std::min((distance_j * value_j) / value_k, distance_k);
                  } else {
                    // k needs to decrease
                    double distance_k = clo[k] - csol[k];
                    movement = std::max((distance_j * value_j) / value_k, distance_k);
                  }
                  if (relEq(movement, 0))
                    continue;

                  csol[k] += movement;
                  if (relEq(csol[k], clo[k])) {
                    colstat[k] = CoinPrePostsolveMatrix::atLowerBound;
                  } else if (relEq(csol[k], cup[k])) {
                    colstat[k] = CoinPrePostsolveMatrix::atUpperBound;
                  } else if (colstat[k] != CoinPrePostsolveMatrix::isFree) {
                    colstat[k] = CoinPrePostsolveMatrix::basic;
                  }
                  printf("NDUAL: x<%d> moved %g to %g; ",
                    k, movement, csol[k]);
                  printf("lb = %g, ub = %g, status now %s.\n",
                    clo[k], cup[k], columnStatusString(k));
                  distance_j -= (movement * value_k) / value_j;
                  csol[j] -= (movement * value_k) / value_j;
                  if (distance_j < 1.0e-7)
                    break;
                }
              }
            }
#endif // repair solution.

            csol[j] = clo[j];
            colstat[j] = CoinPrePostsolveMatrix::atLowerBound;
          }
        }
      } else if (cbarjmax < -ztolcbarj && nflagu == 0 && !prob->colProhibited2(j)) {
        /*
  If cbarmax<j> < 0 (strict inequality) then x<j> NBUB at optimality. If u<j>
  is infinity, notify the user and set the status to unbounded.
*/
        if (cup[j] >= ekkinf) {
          CoinMessageHandler *msghdlr = prob->messageHandler();
          msghdlr->message(COIN_PRESOLVE_COLUMNBOUNDA, prob->messages())
            << j << CoinMessageEol;
          prob->status_ |= 2;
          break;
        } else {
          fix_cols[nfixup_cols++] = j;
#if PRESOLVE_DEBUG > 1
          std::cout << "NDUAL(fix u): fix x(" << j << ")";
          if (csol)
            std::cout << " = " << csol[j];
          std::cout
            << " at u(" << j << ") = " << cup[j] << "; cbar("
            << j << ") < " << cbarjmax << "." << std::endl;
#endif
          if (csol) {
            csol[j] = cup[j];
            colstat[j] = CoinPrePostsolveMatrix::atUpperBound;
          }
        }
      } // end cbar<j> < 0
    }
    /*
  That's the end of this walk through the columns.
*/
    // I don't know why I stopped doing this.
#if PRESOLVE_TIGHTEN_DUALS > 0
    const double *rowels = prob->rowels_;
    const int *hcol = prob->hcol_;
    const CoinBigIndex *mrstrt = prob->mrstrt_;
    int *hinrow = prob->hinrow_;
    // tighten row dual bounds, as described on p. 229
    for (int i = 0; i < nrows; i++) {
      const bool no_ub = (rup[i] >= ekkinf);
      const bool no_lb = (rlo[i] <= -ekkinf);

      if ((no_ub ^ no_lb) == true) {
        const CoinBigIndex krs = mrstrt[i];
        const CoinBigIndex kre = krs + hinrow[i];
        double rmax = ymax[i];
        double rmin = ymin[i];

        // all row columns are non-empty
        for (CoinBigIndex k = krs; k < kre; k++) {
          const double coeff = rowels[k];
          const int icol = hcol[k];
          const double cbarmax0 = cbarmax[icol];
          const double cbarmin0 = cbarmin[icol];

          if (no_ub) {
            // cbarj must not be negative
            if (coeff > ZTOLDP2 && cbarmax0 < PRESOLVE_INF && cup[icol] >= ekkinf) {
              const double bnd = cbarmax0 / coeff;
              if (rmax > bnd) {
#if PRESOLVE_DEBUG > 1
                printf("MAX TIGHT[%d,%d]: %g --> %g\n", i, hrow[k], ymax[i], bnd);
#endif
                ymax[i] = rmax = bnd;
                tightened++;
                ;
              }
            } else if (coeff < -ZTOLDP2 && cbarmax0 < PRESOLVE_INF && cup[icol] >= ekkinf) {
              const double bnd = cbarmax0 / coeff;
              if (rmin < bnd) {
#if PRESOLVE_DEBUG > 1
                printf("MIN TIGHT[%d,%d]: %g --> %g\n", i, hrow[k], ymin[i], bnd);
#endif
                ymin[i] = rmin = bnd;
                tightened++;
                ;
              }
            }
          } else { // no_lb
            // cbarj must not be positive
            if (coeff > ZTOLDP2 && cbarmin0 > -PRESOLVE_INF && clo[icol] <= -ekkinf) {
              const double bnd = cbarmin0 / coeff;
              if (rmin < bnd) {
#if PRESOLVE_DEBUG > 1
                printf("MIN1 TIGHT[%d,%d]: %g --> %g\n", i, hrow[k], ymin[i], bnd);
#endif
                ymin[i] = rmin = bnd;
                tightened++;
                ;
              }
            } else if (coeff < -ZTOLDP2 && cbarmin0 > -PRESOLVE_INF && clo[icol] <= -ekkinf) {
              const double bnd = cbarmin0 / coeff;
              if (rmax > bnd) {
#if PRESOLVE_DEBUG > 1
                printf("MAX TIGHT1[%d,%d]: %g --> %g\n", i, hrow[k], ymax[i], bnd);
#endif
                ymax[i] = rmax = bnd;
                tightened++;
                ;
              }
            }
          }
        }
      }
    }
#endif // PRESOLVE_TIGHTEN_DUALS
/*
  Is it productive to continue with another pass? Essentially, we need lots of
  tightening but no fixing. If we fixed any variables, break and process them.
*/
#if PRESOLVE_DEBUG > 1
    std::cout
      << "NDUAL: pass " << nPass << ": fixed " << (ncols - nfixdown_cols)
      << " down, " << nfixup_cols << " up, tightened " << tightened
      << " duals." << std::endl;
#endif
    if (tightened < 100 || nfixdown_cols < ncols || nfixup_cols)
      break;
  }
  assert(nfixup_cols <= nfixdown_cols);
  /*
  Process columns fixed at upper bound.
*/
  if (nfixup_cols) {
#if PRESOLVE_DEBUG > 1
    std::cout << "NDUAL(upper):";
    for (int k = 0; k < nfixup_cols; k++)
      std::cout << " " << fix_cols[k];
    std::cout << "." << std::endl;
#endif
    next = make_fixed_action::presolve(prob, fix_cols, nfixup_cols, false, next);
  }
  /*
  Process columns fixed at lower bound.
*/
  if (nfixdown_cols < ncols) {
    int *fixdown_cols = fix_cols + nfixdown_cols;
    nfixdown_cols = ncols - nfixdown_cols;
#if PRESOLVE_DEBUG > 1
    std::cout << "NDUAL(lower):";
    for (int k = 0; k < nfixdown_cols; k++)
      std::cout << " " << fixdown_cols[k];
    std::cout << "." << std::endl;
#endif
    next = make_fixed_action::presolve(prob, fixdown_cols,
      nfixdown_cols, true, next);
  }
  /*
  Now look for variables that, when moved in the favourable direction
  according to reduced cost, will naturally tighten an inequality to an
  equality. We can convert that inequality to an equality. See the comments
  at the head of the routine.

  Start with logicals. Open a loop and look for suitable rows. At the
  end of this loop, rows marked with +/-1 will be forced to equality by
  their logical. Rows marked with +/-2 are inequalities that (perhaps)
  can be forced to equality using architecturals. Rows marked with 0 are
  not suitable (range or nonbinding).

  Note that usefulRowInt_ is 3*nrows_; we'll use the second partition below.
*/
  int *canFix = prob->usefulRowInt_;
  for (int i = 0; i < nrows; i++) {
    const bool no_rlb = (rlo[i] <= -ekkinf);
    const bool no_rub = (rup[i] >= ekkinf);
    canFix[i] = 0;
    if (no_rub && !no_rlb) {
      if (ymin[i] > 0.0)
        canFix[i] = -1;
      else
        canFix[i] = -2;
    } else if (no_rlb && !no_rub) {
      if (ymax[i] < 0.0)
        canFix[i] = 1;
      else
        canFix[i] = 2;
    }
#if PRESOLVE_DEBUG > 1
    if (abs(canFix[i]) == 1) {
      std::cout
        << "NDUAL(eq): candidate row (" << i << ") (" << rlo[i]
        << "," << rup[i] << "), logical must be "
        << ((canFix[i] == -1) ? "NBUB" : "NBLB") << "." << std::endl;
    }
#endif
  }
  /*
  Can we do a similar trick with architectural variables? Here, we're looking
  for x<j> such that
  (1) Exactly one of l<j> or u<j> is infinite.
  (2) Moving x<j> towards the infinite bound can tighten exactly one
      constraint i to equality.  If x<j> is entangled with other constraints,
      moving x<j> towards the infinite bound will loosen those constraints.
  (3) Moving x<j> towards the infinite bound is a good idea according to the
      cost c<j> (note we don't have to consider reduced cost here).
  If we can find a suitable x<j>, constraint i can become an equality.

  This is yet another instance of bound propagation, but we're looking for
  a very specific pattern: A variable that can be increased arbitrarily
  in all rows it's entangled with, except for one, which bounds it. And
  we're going to push the variable so as to make that row an equality.
  But note what we're *not* doing: No actual comparison of finite bound
  values to the amount necessary to force an equality. So no worries about
  accuracy, the bane of bound propagation.

  Open a loop to scan columns. bindingUp and bindingDown indicate the result
  of the analysis; -1 says `possible', -2 is ruled out. Test first for
  condition (1). Column singletons are presumably handled elsewhere. Integer
  variables need not apply. If both bounds are finite, no need to look
  further.
*/
  for (int j = 0; j < ncols; j++) {
    if (hincol[j] <= 1)
      continue;
    if (integerType[j])
      continue;
    int bindingUp = -1;
    int bindingDown = -1;
    if (cup[j] < ekkinf)
      bindingUp = -2;
    if (clo[j] > -ekkinf)
      bindingDown = -2;
    if (bindingUp == -2 && bindingDown == -2)
      continue;
    /*
  Open a loop to walk the column and check for condition (2).

  The test for |canFix[i]| != 2 is a non-interference check. We don't want to
  mess with constraints where we've already decided to use the logical to
  force equality. Nor do we want to deal with range or nonbinding constraints.
*/
    const CoinBigIndex &kcs = mcstrt[j];
    const CoinBigIndex kce = kcs + hincol[j];
    for (CoinBigIndex k = kcs; k < kce; k++) {
      const int &i = hrow[k];
      if (abs(canFix[i]) != 2) {
        bindingUp = -2;
        bindingDown = -2;
        break;
      }
      double aij = colels[k];
      /*
  For a<ij> > 0 in a <= constraint (canFix = 2), the up direction is
  binding. For a >= constraint, it'll be the down direction. If the relevant
  binding code is still -1, set it to the index of the row. Similarly for
  a<ij> < 0.

  If this is the second or subsequent binding constraint in that direction,
  set binding[Up,Down] to -2 (we don't want to get into the business of
  calculating which constraint is actually binding).
*/
      if (aij > 0.0) {
        if (canFix[i] == 2) {
          if (bindingUp == -1)
            bindingUp = i;
          else
            bindingUp = -2;
        } else {
          if (bindingDown == -1)
            bindingDown = i;
          else
            bindingDown = -2;
        }
      } else {
        if (canFix[i] == 2) {
          if (bindingDown == -1)
            bindingDown = i;
          else
            bindingDown = -2;
        } else {
          if (bindingUp == -1)
            bindingUp = i;
          else
            bindingUp = -2;
        }
      }
    }
    if (bindingUp == -2 && bindingDown == -2)
      continue;
    /*
  If bindingUp > -2, then either no constraint provided a bound (-1) or
  there's a single constraint (0 <= i < m) that bounds x<j>.  If we have
  just one binding constraint, check that the reduced cost is favourable
  (c<j> <= 0 for x<j> NBUB at optimum for minimisation). If so, declare
  that we will force the row to equality (canFix[i] = +/-1). Note that we
  don't adjust the primal solution value for x<j>.

  If no constraint provided a bound, we might be headed for unboundedness,
  but leave that for some other code to determine.
*/
    double cj = cost[j];
    if (bindingUp > -2 && cj <= 0.0) {
      if (bindingUp >= 0) {
        canFix[bindingUp] /= 2;
#if PRESOLVE_DEBUG > 1
        std::cout
          << "NDUAL(eq): candidate row (" << bindingUp << ") ("
          << rlo[bindingUp] << "," << rup[bindingUp] << "), "
          << " increasing x(" << j << "), cbar = " << cj << "."
          << std::endl;
      } else {
        std::cout
          << "NDUAL(eq): no binding upper bound for x(" << j
          << "), cbar = " << cj << "." << std::endl;
#endif
      }
    } else if (bindingDown > -2 && cj >= 0.0) {
      if (bindingDown >= 0) {
        canFix[bindingDown] /= 2;
#if PRESOLVE_DEBUG > 1
        std::cout
          << "NDUAL(eq): candidate row (" << bindingDown << ") ("
          << rlo[bindingDown] << "," << rup[bindingDown] << "), "
          << " decreasing x(" << j << "), cbar = " << cj << "."
          << std::endl;
      } else {
        std::cout
          << "NDUAL(eq): no binding lower bound for x(" << j
          << "), cbar = " << cj << "." << std::endl;
#endif
      }
    }
  }
  /*
  We have candidate rows. We've avoided scanning full rows until now,
  but there's one remaining hazard: if the row contains unfixed integer
  variables then we don't want to just pin the row to a fixed rhs; that
  might prevent us from achieving integrality. Scan canFix, count and
  record suitable rows (use the second partition of usefulRowInt_).
*/
  // get min max stuff - we can re-use ymin/ymax
  int *infCount = prob->usefulRowInt_ + 2 * nrows;
  {
    // copied from CglProbing
    int i, j;
    CoinBigIndex k, krs, kre;
    int iflagu, iflagl;
    double dmaxup, dmaxdown;

    for (i = 0; i < nrows; ++i) {
      if (rlo[i] == rup[i]) {
        infCount[i] = 10 | (10 << 16);
      } else if (rlo[i] > -1.0e20 || rup[i] < 1.0e20) {
        iflagu = 0;
        iflagl = 0;
        dmaxup = 0.0;
        dmaxdown = 0.0;
        krs = mrstrt[i];
        kre = mrstrt[i] + hinrow[i];

        /* ------------------------------------------------------------*/
        /* Compute L(i) and U(i) */
        /* ------------------------------------------------------------*/
        for (k = krs; k < kre; ++k) {
          double value = rowels[k];
          j = hcol[k];
          if (value > 0.0) {
            if (cup[j] < 1.0e12)
              dmaxup += cup[j] * value;
            else
              ++iflagu;
            if (clo[j] > -1.0e12)
              dmaxdown += clo[j] * value;
            else
              ++iflagl;
          } else if (value < 0.0) {
            if (cup[j] < 1.0e12)
              dmaxdown += cup[j] * value;
            else
              ++iflagl;
            if (clo[j] > -1.0e12)
              dmaxup += clo[j] * value;
            else
              ++iflagu;
          }
        }
        iflagl = std::min(iflagl, 2);
        iflagu = std::min(iflagu, 2);
        infCount[i] = iflagl | (iflagu << 16);
        ymax[i] = dmaxup;
        ymin[i] = dmaxdown;
      } else {
        infCount[i] = 10 | (10 << 16);
      }
    }
  }
  for (int j = 0; j < ncols; j++) {
    if (hincol[j] == 1 && cost[j] > -1.0e10) {
      // can we make equality row
      double coeff = colels[mcstrt[j]];
      int irow = hrow[mcstrt[j]];
      int iflagl = infCount[irow] & 63;
      int iflagu = infCount[irow] >> 16;
      if (iflagu > 1 && iflagl > 1)
        continue;
      double dmaxdown = ymin[irow];
      double dmaxup = ymax[irow];
#if USE_ACTIVE > 1
      printf("singleton col %d has %g on row %d - rhs %g,%g infs %d,%d maxes %g,%g\n",
        j, coeff, irow, rlo[irow], rup[irow],
        iflagl, iflagu,
        dmaxdown, dmaxup);
#endif
      if (coeff > 0.0) {
        if (clo[j] > -1.0e12) {
          dmaxdown -= coeff * clo[j];
          if (iflagl)
            dmaxdown = -1.0e60;
        } else if (iflagl > 1) {
          dmaxdown = -1.0e60;
        }
        if (cup[j] < 1.0e12) {
          dmaxup -= coeff * cup[j];
          if (iflagu)
            dmaxup = 1.0e60;
        } else if (iflagu > 1) {
          dmaxup = 1.0e60;
        }
      } else {
        if (cup[j] < 1.0e12) {
          dmaxdown -= coeff * cup[j];
          if (iflagl)
            dmaxdown = -1.0e60;
        } else if (iflagl > 1) {
          dmaxdown = -1.0e60;
        }
        if (clo[j] > -1.0e12) {
          dmaxup -= coeff * clo[j];
          if (iflagu)
            dmaxup = 1.0e60;
        } else if (iflagu > 1) {
          dmaxup = 1.0e60;
        }
      }
      double rhs;
      if (cost[j] > 0.0) {
        // we want as small as possible
        if (coeff > 0) {
          if (rlo[irow] > -1.0e12) {
	    rhs = rlo[irow];
            if (clo[j] * coeff + dmaxdown < rhs + 1.0e-7
              && cup[j] * coeff + dmaxup > rhs - 1.0e-7) {
              //printf("making rlo equality\n");
              canFix[irow] = -11;
              infCount[irow] = 10 | (10 << 16);
            }
	  }
	} else {
	  // negative coefficient
	  if (rup[irow] < 1.0e12) {
	    rhs = rup[irow];
	    if (clo[j] * coeff + dmaxdown > rhs - 1.0e-7
		&& cup[j] * coeff + dmaxup < rhs + 1.0e-7) {
	      //printf("making rup equality\n");
	      canFix[irow] = 1;
	      infCount[irow] = 10 | (10 << 16);
	    }
	  }
	}
      } else {
        // we want as large as possible
        if (coeff > 0) {
          if (rup[irow] < 1.0e12) {
            rhs = rup[irow];
            if (cup[j] * coeff + dmaxdown > rhs - 1.0e-7 && clo[j] * coeff + dmaxup < rhs + 1.0e-7) {
              //printf("making rup equality\n");
              canFix[irow] = 1;
              infCount[irow] = 10 | (10 << 16);
            }
            rhs /= coeff;
          } else {
            if (cup[j] * coeff + dmaxdown > 1.0e30) {
              //printf("unbounded?\n");
            }
          }
        } else {
          if (rlo[irow] > -1.0e12) {
            rhs = rlo[irow];
            if (cup[j] * coeff + dmaxdown < rhs + 1.0e-7
              && clo[j] * coeff + dmaxup > rhs - 1.0e-7) {
#if USE_ACTIVE > 1
              printf("singleton col %d has %g on row %d - rhs %g,%g infs %d,%d maxes %g,%g\n",
                j, coeff, irow, rlo[irow], rup[irow],
                iflagl, iflagu,
                dmaxdown, dmaxup);
              printf("making rlo equality\n");
#endif
              canFix[irow] = -11;
              infCount[irow] = 10 | (10 << 16);
            }
            rhs /= coeff;
          } else {
            if (cup[j] * coeff + dmaxdown > 1.0e30) {
              //printf("unbounded?\n");
            }
          }
        }
      }
#if USE_ACTIVE > 1
      printf("adjusted maxes %g,%g - rhs %g\n", dmaxdown, dmaxup, rhs);
#endif
    }
  }
#if PRESOLVE_DEBUG > 0
  int makeEqCandCnt = 0;
  for (int i = 0; i < nrows; i++) {
    if (abs(canFix[i]) == 1)
      makeEqCandCnt++;
  }
#endif
  int makeEqCnt = nrows;
  for (int i = 0; i < nrows; i++) {
    if (abs(canFix[i]) == 1 && hinrow[i] > 1) {
      const CoinBigIndex krs = mrstrt[i];
      const CoinBigIndex kre = krs + hinrow[i];
      int nBinary = 0;
      int nOther = 0;
      int nSameCoeff = 0;
      double rhs = canFix[i] == 1 ? rup[i] : rlo[i];
      double firstCoeff = 0.0;
      int j = -1;
      for (CoinBigIndex k = krs; k < kre; k++) {
        j = hcol[k];
        double coeff = rowels[k];
        if (cup[j] > clo[j]) {
          if (!firstCoeff) {
            firstCoeff = coeff;
            nSameCoeff = 1;
          } else {
            if (coeff == firstCoeff)
              nSameCoeff++;
            else if (coeff != -firstCoeff)
              nOther = 1;
          }
          if (prob->colProhibited2(j)) {
            nBinary = 1;
            nOther = 1;
            break;
          } else if (integerType[j]) {
            if (cup[j] - clo[j] == 1.0) {
              nBinary++;
            } else {
              nBinary = 1;
              nOther = 1;
              break;
            }
          } else {
            nOther++;
          }
        } else {
          rhs -= coeff * clo[j];
        }
      }
      bool canFixThis = true;
      if (nBinary) {
#if PRESOLVE_DEBUG > 1
	printf("poss eq row %d is %g <= ", i,rlo[i]);
	for (CoinBigIndex k = krs; k < kre; k++) {
	  const double coeff = rowels[k];
	  const int icol = hcol[k];
	  printf("(%d,%g - bounds %g,%g%s nel %d) ", icol, coeff, clo[icol], cup[icol],integerType[icol] ? "(I)" : "",hincol[icol]);
	}
	printf("<= %g\n", rup[i]);
#endif
        if (nOther) {
          canFixThis = false;
#if PRESOLVE_DEBUG > 1
          std::cout
            << "NDUAL(eq): cannot convert row " << i << " to equality; "
            << "unfixed integer variable x(" << j << ")." << std::endl;
#endif
        } else {
          // may be able to if all easy integer
          if (!rhs) {
            // for now just two of opposite sign
            if (nSameCoeff == 1 && nBinary == 2) {
              // ok
#if PRESOLVE_DEBUG > 1
              printf("int eq row is %g <= ", rlo[i]);
              for (CoinBigIndex k = krs; k < kre; k++) {
                const double coeff = rowels[k];
                const int icol = hcol[k];
                printf("(%d,%g - bounds %g,%g) ", icol, coeff, clo[icol], cup[icol]);
              }
              printf("<= %g\n", rup[i]);
#endif
            } else {
              canFixThis = false;
#if PRESOLVE_DEBUG > 1
              std::cout
                << "NDUAL(eq2): cannot convert row " << i << " to equality; "
                << "unfixed integer variable x(" << j << ")." << std::endl;
#endif
            }
          } else if (rhs == 1.0 && canFix[i] == 1 && nSameCoeff == nBinary && firstCoeff == 1.0) {
            // ok
#if PRESOLVE_DEBUG > 1
            printf("int2 eq row is %g <= ", rlo[i]);
            for (CoinBigIndex k = krs; k < kre; k++) {
              const double coeff = rowels[k];
              const int icol = hcol[k];
              printf("(%d,%g - bounds %g,%g) ", icol, coeff, clo[icol], cup[icol]);
            }
            printf("<= %g\n", rup[i]);
#endif
          } else {
            canFixThis = false;
#if PRESOLVE_DEBUG > 1
            std::cout
              << "NDUAL(eq3): cannot convert row " << i << " to equality; "
              << "unfixed integer variable x(" << j << ")." << std::endl;
#endif
          }
        }
      }
      if (canFixThis) {
#if PRESOLVE_DEBUG > 1
	if (nBinary)
	  printf(" -- good row %d --\n",i);
#endif
        canFix[makeEqCnt++] = i;
      } else {
	// try a different way
	int nInteger = 0;
	int nOther = 0;
	int nSameCoeff = 0;
	int nSingletons = 0;
	// Singletons - all continuous 1 - all binary 2 - other (general or mixed) 3
	int nMixed = 0;
	double rhs = canFix[i] == 1 ? rup[i] : rlo[i];
	double firstCoeff = 0.0;
	int j = -1;
#if PRESOLVE_DEBUG > 1
	printf("poss eq row %d is %g <= ", i,rlo[i]);
	for (CoinBigIndex k = krs; k < kre; k++) {
	  const double coeff = rowels[k];
	  const int icol = hcol[k];
	  printf("(%d,%g - bounds %g,%g%s nel %d) ", icol, coeff, clo[icol], cup[icol],integerType[icol] ? "(I)" : "",hincol[icol]);
	}
	printf("<= %g\n", rup[i]);
#endif
	for (CoinBigIndex k = krs; k < kre; k++) {
	  j = hcol[k];
	  if (hincol[j]==1) {
	    nSingletons++;
	  }
	  double coeff = fabs(rowels[k]);
	  if (cup[j] > clo[j]) {
	    if (!firstCoeff) {
	      firstCoeff = coeff;
	      nSameCoeff = 1;
	    } else {
	      if (coeff == firstCoeff)
		nSameCoeff++;
	      else if (coeff != -firstCoeff)
		nOther = 1;
	    }
	    if (prob->colProhibited2(j)) {
	      nInteger = 1;
	      nOther = 1;
	      nMixed = 3;
	      break;
	    } else if (integerType[j]) {
	      nInteger++;
	      if (clo[j]==0.0&&cup[j]==1.0) {
		if (!nMixed)
		  nMixed = 2;
		else if (nMixed!=2)
		  nMixed = 3; // be on safe side
	      } else {
		nMixed = 3;
	      }
	    } else {
	      nOther++;
	      if (!nMixed)
		nMixed = 1;
	      else if (nMixed!=1)
		nMixed = 3; // be on safe side
	    }
	  } else {
	    rhs -= coeff * clo[j];
	  }
	}
	if (nSingletons>1) {
	  //printf("%d singletons nOther %d\n",nSingletons,nOther);
	  if (nMixed<3)
	    nOther=0;
	}
	if (!nOther) {
	  if (nSameCoeff == nInteger) {
	    rhs /= firstCoeff;
	    if (fabs(rhs-floor(rhs+0.5))<1.0e-9) {
	    // ok
#if PRESOLVE_DEBUG > 1
              printf("int eq row is %g <= ", rlo[i]);
              for (CoinBigIndex k = krs; k < kre; k++) {
                const double coeff = rowels[k];
                const int icol = hcol[k];
                printf("(%d,%g - bounds %g,%g%s) ", icol, coeff, clo[icol], cup[icol],integerType[icol] ? "(I)" : "");
              }
              printf("<= %g\n", rup[i]);
#endif
	      canFix[makeEqCnt++] = i;
	    }
	  }
	}
      }
    }
  }
  makeEqCnt -= nrows;
#if PRESOLVE_DEBUG > 0
  if ((makeEqCandCnt - makeEqCnt) > 0) {
    std::cout
      << "NDUAL(eq): rejected " << (makeEqCandCnt - makeEqCnt)
      << " rows due to unfixed integer variables." << std::endl;
  }
#endif
  /*
  If we've identified inequalities to convert, do the conversion, record
  the information needed to restore bounds in postsolve, and finally create
  the postsolve object.
*/
  if (makeEqCnt > 0) {
    action *bndRecords = new action[makeEqCnt];
    for (int k = 0; k < makeEqCnt; k++) {
      const int i = canFix[k + nrows];
#if PRESOLVE_DEBUG > 1
      std::cout << "NDUAL(eq): forcing row " << i << " to equality;";
      if (canFix[i] == -1)
        std::cout << " dropping b = " << rup[i] << " to " << rlo[i];
      else
        std::cout << " raising blow = " << rlo[i] << " to " << rup[i];
      std::cout << "." << std::endl;
#endif
      action &bndRec = bndRecords[(k)];
      bndRec.rlo_ = rlo[i];
      bndRec.rup_ = rup[i];
      bndRec.ndx_ = i;
      if (canFix[i] == 1) {
        rlo[i] = rup[i];
        prob->addRow(i);
      } else if (canFix[i] == -1) {
        rup[i] = rlo[i];
        prob->addRow(i);
      }
    }
    next = new remove_dual_action(makeEqCnt, bndRecords, next);
  }

#if PRESOLVE_TIGHTEN_DUALS > 0
  delete[] cbarmin;
  delete[] cbarmax;
#endif
#undef PRESOLVE_DEBUG
#if COIN_PRESOLVE_TUNING > 0
  double thisTime = 0.0;
  if (prob->tuning_)
    thisTime = CoinCpuTime();
#endif
#if PRESOLVE_DEBUG > 0 || PRESOLVE_CONSISTENCY > 0
  presolve_check_sol(prob, 2, 1, 1);
  presolve_check_nbasic(prob);
#endif
#if PRESOLVE_DEBUG > 0 || COIN_PRESOLVE_TUNING > 0
  int droppedRows = prob->countEmptyRows() - startEmptyRows;
  int droppedColumns = prob->countEmptyCols() - startEmptyColumns;
  std::cout
    << "Leaving remove_dual_action::presolve, dropped " << droppedRows
    << " rows, " << droppedColumns << " columns, forced "
    << makeEqCnt << " equalities";
#if COIN_PRESOLVE_TUNING > 0
  if (prob->tuning_)
    std::cout << " in " << (thisTime - prob->startTime_) << "s";
#endif
  std::cout << "." << std::endl;
#endif

  return (next);
}

/*
  Postsolve: replace the original row bounds.

  The catch here is that each constraint was an equality in the presolved
  problem, with a logical s<i> that had l<i> = u<i> = 0. We're about to
  convert the equality back to an inequality. One row bound will go to
  infinity, as will one of the bounds of the logical. We may need to patch the
  basis. The logical for a <= constraint cannot be NBUB, and the logical for a
  >= constraint cannot be NBLB.
*/
void remove_dual_action::postsolve(CoinPostsolveMatrix *prob) const
{
  const action *const &bndRecords = actions_;
  const int &numRecs = nactions_;

  double *&rlo = prob->rlo_;
  double *&rup = prob->rup_;
  unsigned char *&rowstat = prob->rowstat_;

#if PRESOLVE_CONSISTENCY > 0 || PRESOLVE_DEBUG > 0
#if PRESOLVE_DEBUG > 0
  std::cout
    << "Entering remove_dual_action::postsolve, " << numRecs
    << " bounds to restore." << std::endl;
#endif
  presolve_check_threads(prob);
  presolve_check_sol(prob, 2, 2, 2);
  presolve_check_nbasic(prob);
#endif

  /*
  For each record, restore the row bounds. If we have status arrays, check
  the status of the logical and adjust if necessary.

  In spite of the fact that the status array is an unsigned char array,
  we still need to use getRowStatus to make sure we're only looking at the
  bottom three bits. Why is this an issue? Because the status array isn't
  necessarily cleared to zeros, and setRowStatus carefully changes only
  the bottom three bits!
*/
  for (int k = 0; k < numRecs; k++) {
    const action &bndRec = bndRecords[k];
    const int &i = bndRec.ndx_;
    const double &rloi = bndRec.rlo_;
    const double &rupi = bndRec.rup_;

#if PRESOLVE_DEBUG > 1
    std::cout << "NDUAL(eq): row(" << i << ")";
    if (rlo[i] != rloi)
      std::cout << " LB " << rlo[i] << " -> " << rloi;
    if (rup[i] != rupi)
      std::cout << " UB " << rup[i] << " -> " << rupi;
#endif

    rlo[i] = rloi;
    rup[i] = rupi;
    if (rowstat) {
      unsigned char stati = prob->getRowStatus(i);
      if (stati == CoinPresolveMatrix::atLowerBound) {
        if (rloi <= -PRESOLVE_INF) {
          rowstat[i] = CoinPresolveMatrix::atUpperBound;
#if PRESOLVE_DEBUG > 1
          std::cout
            << ", status forced to "
            << statusName(static_cast< CoinPresolveMatrix::Status >(rowstat[i]));
#endif
        }
      } else if (stati == CoinPresolveMatrix::atUpperBound) {
        if (rupi >= PRESOLVE_INF) {
          rowstat[i] = CoinPresolveMatrix::atLowerBound;
#if PRESOLVE_DEBUG > 1
          std::cout
            << ", status forced to "
            << statusName(static_cast< CoinPresolveMatrix::Status >(rowstat[i]));
#endif
        }
      }
#if PRESOLVE_DEBUG > 2
      else if (stati == CoinPresolveMatrix::basic) {
        std::cout << ", status is basic.";
      } else if (stati == CoinPresolveMatrix::isFree) {
        std::cout << ", status is free?!";
      } else {
        unsigned int tmp = static_cast< unsigned int >(stati);
        std::cout << ", status is invalid (" << tmp << ")!";
      }
#endif
    }
#if PRESOLVE_DEBUG > 1
    std::cout << "." << std::endl;
#endif
  }

#if PRESOLVE_CONSISTENCY > 0 || PRESOLVE_DEBUG > 0
  presolve_check_threads(prob);
  presolve_check_sol(prob, 2, 2, 2);
  presolve_check_nbasic(prob);
#if PRESOLVE_DEBUG > 0
  std::cout << "Leaving remove_dual_action::postsolve." << std::endl;
#endif
#endif

  return;
}

/*
  Destructor
*/
remove_dual_action::~remove_dual_action()
{
  deleteAction(actions_, action *);
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
