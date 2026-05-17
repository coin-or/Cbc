// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdio>

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>

#include "CoinPragma.hpp"
#include "CoinPackedMatrix.hpp"
#include "CglOddHole.hpp"


//--------------------------------------------------------------------------
// test EKKsolution methods.
void
CglOddHoleUnitTest(
  const OsiSolverInterface * baseSiP,
  const std::string mpsDir )
{
  CoinRelFltEq eq(0.000001);

  // Test default constructor
  {
    CglOddHole aGenerator;
  }
  
  // Test copy & assignment
  {
    CglOddHole rhs;
    {
      CglOddHole bGenerator;
      CglOddHole cGenerator(bGenerator);
      rhs=bGenerator;
    }
  }


  // test on simple case
  {  
    const int nRows=3;
    const int nCols=3;
    const int nEls=6;
    const double elem[]={1.0,1.0,1.0,1.0,1.0,1.0};
    const int row[]={0,1,0,2,1,2};
    const CoinBigIndex start[]={0,2,4};
    const int len[]={2,2,2};
    CoinPackedMatrix matrix(true,nRows,nCols,nEls,elem,row,start,len);
    const double sol[]={0.5,0.5,0.5};
    const double dj[]={0,0,0};
    const int which[]={1,1,1};
    const int fixed[]={0,0,0};
    OsiCuts cs;
    CglOddHole test1;
    CglTreeInfo info;
    info.randomNumberGenerator=NULL;
    test1.generateCuts(NULL,matrix,sol,dj,cs,which,fixed,info,true);
    CoinPackedVector check;
    int index[] = {0,1,2};
    double el[] = {1,1,1};
    check.setVector(3,index,el);
    //assert (cs.sizeRowCuts()==2);
    assert (cs.sizeRowCuts()==1);
    // sort Elements in increasing order
    CoinPackedVector rpv=cs.rowCut(0).row();
    rpv.sortIncrIndex();
    assert (check==rpv);
  }
  
  // Testcase /u/rlh/osl2/mps/scOneInt.mps
  // Model has 3 continous, 2 binary, and 1 general
  // integer variable.
  {
    OsiSolverInterface  * siP = baseSiP->clone();
    
    std::string fn = mpsDir+"scOneInt";
    siP->readMps(fn.c_str(),"mps");
#if 0
    CglOddHole cg;
    int nCols=siP->getNumCols();
    
    // Test the siP methods for detecting
    // variable type
    int numCont=0, numBinary=0, numIntNonBinary=0, numInt=0;
    for (int thisCol=0; thisCol<nCols; thisCol++) {
      if ( siP->isContinuous(thisCol) ) numCont++;
      if ( siP->isBinary(thisCol) ) numBinary++;
      if ( siP->isIntegerNonBinary(thisCol) ) numIntNonBinary++;
      if ( siP->isInteger(thisCol) ) numInt++;
    }
    assert(numCont==3);
    assert(numBinary==2);
    assert(numIntNonBinary==1);
    assert(numInt==3);
    
    
    // Test initializeCutGenerator
    siP->initialSolve();
    assert(xstar !=NULL);
    for (i=0; i<nCols; i++){
      assert(complement[i]==0);
    }
    int nRows=siP->getNumRows();
    for (i=0; i<nRows; i++){
    int vectorsize = siP->getMatrixByRow()->getVectorSize(i);
    assert(vectorsize==2);
    }
    
    kccg.cleanUpCutGenerator(complement,xstar);
#endif  
    delete siP;
  }

}

