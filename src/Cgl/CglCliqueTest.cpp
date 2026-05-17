// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdlib>
#include <cstdio>

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>
#include "CoinPragma.hpp"
#include "CglClique.hpp"


void
CglCliqueUnitTest(const OsiSolverInterface *baseSiP,
			    const std::string mpsDir)
{
  // Test default constructor
  {
    CglClique aGenerator;
  }
  
  // Test copy & assignment
  {
    CglClique rhs;
    {
      CglClique bGenerator;
      CglClique cGenerator(bGenerator);
      //rhs=bGenerator;
    }
  }

  // Test get/set methods
  {
    CglClique getset;
    // None to test
  }

  // Test generateCuts
  {
    CglClique gct;
    OsiSolverInterface  *siP = baseSiP->clone();
    std::string fn = mpsDir+"l152lav";
    std::string fn2 = mpsDir+"l152lav.mps";
    FILE *in_f = fopen(fn2.c_str(), "r");
    if(in_f == NULL) {
      std::cout<<"Can not open file "<<fn2<<std::endl<<"Skip test of CglClique::generateCuts()"<<std::endl;
    }
    else {
      fclose(in_f);
      siP->readMps(fn.c_str(),"mps");
 
      siP->initialSolve();
      double lpRelax = siP->getObjValue();
      
      OsiCuts cs;
      gct.generateCuts(*siP, cs);
      int nRowCuts = cs.sizeRowCuts();
      std::cout<<"There are "<<nRowCuts<<" Clique cuts"<<std::endl;
      assert(cs.sizeRowCuts() > 0);
      OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cs);
      
      siP->resolve();
      
      double lpRelaxAfter= siP->getObjValue(); 
      std::cout<<"Initial LP value: "<<lpRelax<<std::endl;
      std::cout<<"LP value with cuts: "<<lpRelaxAfter<<std::endl;
      assert( lpRelax < lpRelaxAfter );
      assert(lpRelaxAfter < 4722.1);
    }
    delete siP;
  }

}

