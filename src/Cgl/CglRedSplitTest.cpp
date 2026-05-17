//
// Name:     CglRedSplit.cpp
// Author:   Francois Margot                                                  
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     2/6/05
//
//---------------------------------------------------------------------------
// Copyright (C) 2005, Francois Margot and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdlib>
#include <cstdio>

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>
#include "CoinPragma.hpp"
#include "CglRedSplit.hpp"


void
CglRedSplitUnitTest(const OsiSolverInterface *baseSiP,
		    const std::string mpsDir)
{
  // Test default constructor
  {
    CglRedSplit aGenerator;
  }
  
  // Test copy & assignment
  {
    CglRedSplit rhs;
    {
      CglRedSplit bGenerator;
      CglRedSplit cGenerator(bGenerator);
      rhs=bGenerator;
    }
  }

  // Test get/set methods
  {
    CglRedSplit getset;
    CglRedSplitParam gsparam = getset.getParam();
    
    double geps = 10 * gsparam.getEPS();
    gsparam.setEPS(geps);
    double geps2 = gsparam.getEPS();
    assert(geps == geps2);

    double gepse = 10 * gsparam.getEPS_ELIM();
    gsparam.setEPS_ELIM(gepse);
    double gepse2 = gsparam.getEPS_ELIM();
    assert(gepse == gepse2);

    double gmv = 10 * gsparam.getMINVIOL();
    gsparam.setMINVIOL(gmv);
    double gmv2 = gsparam.getMINVIOL();
    assert(gmv == gmv2);

    int gucg = gsparam.getUSE_CG2();
    gucg = 1 - gucg;
    gsparam.setUSE_CG2(gucg);
    int gucg2 = gsparam.getUSE_CG2();
    assert(gucg == gucg2);
  }

  // Test generateCuts
  {
    CglRedSplit gct;
    OsiSolverInterface  *siP = baseSiP->clone();
    std::string fn = mpsDir+"p0033";
    std::string fn2 = mpsDir+"p0033.mps";
    FILE *in_f = fopen(fn2.c_str(), "r");
    if(in_f == NULL) {
      std::cout<<"Can not open file "<<fn2<<std::endl<<"Skip test of CglRedSplit::generateCuts()"<<std::endl;
    }
    else {
      fclose(in_f);
      siP->readMps(fn.c_str(),"mps");
 
      siP->initialSolve();
      double lpRelax = siP->getObjValue();
      
      OsiCuts cs;
      gct.getParam().setMAX_SUPPORT(34);
      gct.getParam().setUSE_CG2(1);
      //      gct.getParam().setUSE_CG2(1);
      gct.generateCuts(*siP, cs);
      int nRowCuts = cs.sizeRowCuts();
      std::cout<<"There are "<<nRowCuts<<" Reduce-and-Split cuts"<<std::endl;
      assert(cs.sizeRowCuts() > 0);
      OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cs);
      
      siP->resolve();
      
      double lpRelaxAfter= siP->getObjValue(); 
      std::cout<<"Initial LP value: "<<lpRelax<<std::endl;
      std::cout<<"LP value with cuts: "<<lpRelaxAfter<<std::endl;
      assert( lpRelax < lpRelaxAfter );
      assert(lpRelaxAfter < 3089.1);
    }
    delete siP;
  }

}

