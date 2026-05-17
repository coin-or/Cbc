// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdio>

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>
#include "CoinPragma.hpp"
#include "CglMixedIntegerRounding2.hpp"


void
CglMixedIntegerRounding2UnitTest(const OsiSolverInterface *baseSiP,
			    const std::string mpsDir)
{
  // Test default constructor
  {
    CglMixedIntegerRounding2 aGenerator;
  }
  
  // Test copy & assignment
  {
    CglMixedIntegerRounding2 rhs;
    {
      CglMixedIntegerRounding2 bGenerator;
      CglMixedIntegerRounding2 cGenerator(bGenerator);
      rhs=bGenerator;
    }
  }

  // Test get/set methods
  {
    CglMixedIntegerRounding2 getset;

    int gagg = 10 * getset.getMAXAGGR_();
    getset.setMAXAGGR_(gagg);
    int gagg2 = getset.getMAXAGGR_();
    assert(gagg == gagg2);

    bool gmult = !getset.getMULTIPLY_();
    getset.setMULTIPLY_(gmult);
    bool gmult2 = getset.getMULTIPLY_();
    assert(gmult == gmult2);

    int gcrit = getset.getCRITERION_();
    gcrit = (gcrit) % 3 + 1;
    getset.setCRITERION_(gcrit);
    int gcrit2 = getset.getCRITERION_();
    assert(gcrit == gcrit2);

    int gpre = getset.getDoPreproc();
    gpre = (gpre + 1) % 3 - 1;
    getset.setDoPreproc(gpre);
    int gpre2 = getset.getDoPreproc();
    assert(gpre == gpre2);
  }

  // Test generateCuts
  {
    CglMixedIntegerRounding2 gct;
    OsiSolverInterface  *siP = baseSiP->clone();
    std::string fn = mpsDir+"capPlan1";
    std::string fn2 = mpsDir+"capPlan1.mps";
    FILE *in_f = fopen(fn2.c_str(), "r");
    if(in_f == NULL) {
      std::cout<<"Can not open file "<<fn2<<std::endl<<"Skip test of CglMixedIntegerRounding2::generateCuts()"<<std::endl;
    }
    else {
      fclose(in_f);
      siP->readMps(fn.c_str(),"mps");
 
      siP->initialSolve();
      double lpRelax = siP->getObjValue();
      
      OsiCuts cs;
      gct.setDoPreproc(1); // Needed for DyLP
      gct.generateCuts(*siP, cs);
      int nRowCuts = cs.sizeRowCuts();
      std::cout<<"There are "<<nRowCuts<<" MIR2 cuts"<<std::endl;
      assert(cs.sizeRowCuts() > 0);
      OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cs);
      
      siP->resolve();
      
      double lpRelaxAfter= siP->getObjValue(); 
      printf("Initial LP value: %f\n", lpRelax);
      printf("LP value with cuts: %f\n", lpRelaxAfter);
      assert( lpRelax < lpRelaxAfter );
      assert(lpRelaxAfter < 964);
    }
    delete siP;
  }

}

