// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdio>

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>
#include "CoinPragma.hpp"
#include "CglResidualCapacity.hpp"


void
CglResidualCapacityUnitTest(const OsiSolverInterface *baseSiP,
			    const std::string mpsDir)
{
  // Test default constructor
  {
    CglResidualCapacity aGenerator;
  }
  
  // Test copy & assignment
  {
    CglResidualCapacity rhs;
    {
      CglResidualCapacity bGenerator;
      CglResidualCapacity cGenerator(bGenerator);
      rhs=bGenerator;
    }
  }

  // Test get/set methods
  {
    CglResidualCapacity getset;
    
    double geps = 10 * getset.getEpsilon();
    getset.setEpsilon(geps);
    double geps2 = getset.getEpsilon();
    assert(geps == geps2);

    double gtol = 10 * getset.getTolerance();
    getset.setTolerance(gtol);
    double gtol2 = getset.getTolerance();
    assert(gtol == gtol2);

    int gpre = getset.getDoPreproc();
    gpre = (gpre + 1) % 3 - 1;
    getset.setDoPreproc(gpre);
    int gpre2 = getset.getDoPreproc();
    assert(gpre == gpre2);
  }

  // Test generateCuts
  {
    CglResidualCapacity gct;
    OsiSolverInterface  *siP = baseSiP->clone();
    std::string fn = mpsDir+"capPlan1";
    std::string fn2 = mpsDir+"capPlan1.mps";
    FILE *in_f = fopen(fn2.c_str(), "r");
    if(in_f == NULL) {
      std::cout<<"Can not open file "<<fn2<<std::endl<<"Skip test of CglResidualCapacity::generateCuts()"<<std::endl;
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
      std::cout<<"There are "<<nRowCuts<<" Residual Capacity cuts"<<std::endl;
      assert(cs.sizeRowCuts() > 0);
      OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cs);
      
      siP->resolve();
      
      double lpRelaxAfter= siP->getObjValue(); 
      std::cout<<"Initial LP value: "<<lpRelax<<std::endl;
      std::cout<<"LP value with cuts: "<<lpRelaxAfter<<std::endl;
      assert( lpRelax < lpRelaxAfter );
      assert(lpRelaxAfter < 964);
    }
    delete siP;
  }

}

