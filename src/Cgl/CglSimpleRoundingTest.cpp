// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>

#include "CoinPragma.hpp"
#include "CglSimpleRounding.hpp" 
#include <stdio.h>

//--------------------------------------------------------------------------
// test the simple rounding cut generators methods.
void
CglSimpleRoundingUnitTest(
  const OsiSolverInterface * baseSiP,
  const std::string mpsDir )
{

  // Test default constructor
  {
    CglSimpleRounding cg;
  }

  // Test copy & assignment
  {
    CglSimpleRounding rhs;
    {
      CglSimpleRounding cg;
      CglSimpleRounding cgC(cg);
      rhs=cg;
    }
  }

  // Test gcd and gcdn
  {
    CglSimpleRounding cg;
    int v = cg.gcd(122,356);
    assert(v==2);
    v=cg.gcd(356,122);
    assert(v==2);
    v=cg.gcd(54,67);
    assert(v==1);
    v=cg.gcd(67,54);
    assert(v==1);
    v=cg.gcd(485,485);
    assert(v==485);
    v=cg.gcd(17*13,17*23);
    assert( v==17);
    v=cg.gcd(17*13*5,17*23);
    assert( v==17);
    v=cg.gcd(17*13*23,17*23);
    assert(v==17*23);

    int a[4] = {12, 20, 32, 400};
    v= cg.gcdv(4,a);
    assert(v== 4);
    int b[4] = {782, 4692, 51, 2754};
    v= cg.gcdv(4,b);
    assert(v== 17);
    int c[4] = {50, 40, 30, 10};
    v= cg.gcdv(4,c);
    assert(v== 10);
  }


  // Test generate cuts method on exmip1.5.mps
  {
    CglSimpleRounding cg;
    
    OsiSolverInterface * siP = baseSiP->clone();
    std::string fn = mpsDir+"exmip1.5.mps";
    siP->readMps(fn.c_str(),"");
    OsiCuts cuts;
    cg.generateCuts(*siP,cuts);

    // there should be 3 cuts
    int nRowCuts = cuts.sizeRowCuts();
    assert(nRowCuts==3);

    // get the last "sr"=simple rounding cut that was derived
    OsiRowCut srRowCut2 = cuts.rowCut(2); 
    CoinPackedVector srRowCutPV2 = srRowCut2.row();

    // this is what the last cut should look like: i.e. the "solution"
    const int solSize = 2;
    int solCols[solSize]={2,3};
    double solCoefs[solSize]={5.0, 4.0};
    OsiRowCut solRowCut;
    solRowCut.setRow(solSize,solCols,solCoefs);
    solRowCut.setLb(-COIN_DBL_MAX);
    solRowCut.setUb(2.0);

    // Test for equality between the derived cut and the solution cut

    // Note: testing two OsiRowCuts are equal invokes testing two
    // CoinPackedVectors are equal which invokes testing two doubles
    // are equal.  Usually not a good idea to test that two doubles are equal, 
    // but in this cut the "doubles" represent integer values. Also allow that
    // different solvers have different orderings in packed vectors, which may
    // not match the ordering defined for solRowCut.

    assert(srRowCut2.OsiCut::operator==(solRowCut)) ;
    assert(srRowCut2.row().isEquivalent(solRowCut.row())) ;
    assert(srRowCut2.lb() == solRowCut.lb()) ;
    assert(srRowCut2.ub() == solRowCut.ub()) ;

    delete siP;
  }

  // Test generate cuts method on p0033
  {
    CglSimpleRounding cg;
    
    OsiSolverInterface * siP = baseSiP->clone();
    std::string fn = mpsDir+"p0033";
    siP->readMps(fn.c_str(),"mps");
    OsiCuts cuts;
    cg.generateCuts(*siP,cuts);

    // p0033 is the optimal solution to p0033
    int objIndices[14] = { 
       0,  6,  7,  9, 13, 17, 18,
      22, 24, 25, 26, 27, 28, 29 };
    CoinPackedVector p0033(14,objIndices,1.0);

    // test that none of the generated cuts
    // chops off the optimal solution
    int nRowCuts = cuts.sizeRowCuts();
    OsiRowCut rcut;
    CoinPackedVector rpv;
    int i;
    for (i=0; i<nRowCuts; i++){
      rcut = cuts.rowCut(i);
      rpv = rcut.row();
      double p0033Sum = (rpv*p0033).sum();
      double rcutub = rcut.ub();
      assert (p0033Sum <= rcutub);
    }

    // test that the cuts improve the 
    // lp objective function value
    siP->initialSolve();
    double lpRelaxBefore=siP->getObjValue();
    OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);
    siP->resolve();
    double lpRelaxAfter=siP->getObjValue(); 
#ifdef CGL_DEBUG
    printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
    printf("Final LP min=%f\n\n",lpRelaxAfter);
#endif
    assert( lpRelaxBefore < lpRelaxAfter );

    delete siP;

  }


}

