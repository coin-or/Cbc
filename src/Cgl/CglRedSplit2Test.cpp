#include <cstdlib>
#include <cstdio>

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>
#include "CoinPragma.hpp"
#include "CglRedSplit2.hpp"


void
CglRedSplit2UnitTest(const OsiSolverInterface *baseSiP,
		     const std::string mpsDir)
{
  // Test default constructor
  {
    CglRedSplit2 aGenerator;
  }
  
  // Test copy & assignment
  {
    CglRedSplit2 rhs;
    {
      CglRedSplit2 bGenerator;
      CglRedSplit2 cGenerator(bGenerator);
      rhs=bGenerator;
    }
  }

  // Test get/set methods
  {
    CglRedSplit2 getset;
    CglRedSplit2Param gsparam = getset.getParam();
    
    double geps = 1.1 * gsparam.getEPS();
    gsparam.setEPS(geps);
    double geps2 = gsparam.getEPS();
    assert(geps == geps2);

    double gepse = 1.1 * gsparam.getEPS_ELIM();
    gsparam.setEPS_ELIM(gepse);
    double gepse2 = gsparam.getEPS_ELIM();
    assert(gepse == gepse2);

    double gmv = 1.1 * gsparam.getMINVIOL();
    gsparam.setMINVIOL(gmv);
    double gmv2 = gsparam.getMINVIOL();
    assert(gmv == gmv2);

  }

  // Test generateCuts
  {
    CglRedSplit2 gct;
    OsiSolverInterface  *siP = baseSiP->clone();
    std::string fn = mpsDir+"p0033";
    std::string fn2 = mpsDir+"p0033.mps";
    FILE *in_f = fopen(fn2.c_str(), "r");
    if(in_f == NULL) {
      std::cout<<"Can not open file "<<fn2<<std::endl<<"Skip test of CglRedSplit2::generateCuts()"<<std::endl;
    }
    else {
      fclose(in_f);
      siP->readMps(fn.c_str(),"mps");
 
      siP->initialSolve();
      double lpRelax = siP->getObjValue();
      
      OsiCuts cs;
      gct.getParam().setMAX_SUPPORT(34);
      gct.generateCuts(*siP, cs);
      int nRowCuts = cs.sizeRowCuts();
      std::cout<<"There are "<<nRowCuts<<" Reduce-and-Split2 cuts"<<std::endl;
      assert(cs.sizeRowCuts() > 0);
      OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cs);
      
      siP->resolve();
      
      double lpRelaxAfter= siP->getObjValue(); 
      std::cout<<"Initial LP value: "<<lpRelax<<std::endl;
      std::cout<<"LP value with cuts: "<<lpRelaxAfter<<std::endl;
      assert(lpRelax < lpRelaxAfter);
      assert(lpRelaxAfter < 3089.1);
    }
    delete siP;
  }

}

