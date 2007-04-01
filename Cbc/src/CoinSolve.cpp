// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
   
#include "CbcConfig.h"
#include "CoinPragma.hpp"
#include "CbcModel.hpp"
#include "OsiClpSolverInterface.hpp"

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <iostream>
//#############################################################################
// To use USERCBC or USERCLP uncomment the following define and add in your fake main program here
//#define USER_HAS_FAKE_MAIN
//  Start any fake main program
#ifdef USER_HAS_FAKE_MAIN
#endif
//  End any fake main program
//#############################################################################
void CbcClpUnitTest (const CbcModel & saveModel);
int main (int argc, const char *argv[])
{
  OsiClpSolverInterface solver1;
  CbcModel * babSolver = NULL;
  int returnCode = CbcMain (argc, argv,solver1,&babSolver);
  if (returnCode!=777) {
    return returnCode;
  } else {
    // do miplib
    CbcClpUnitTest(*babSolver);
    return 0;
  }
}
/*
  Version 1.00.00 November 16 2005.
  This is to stop me (JJF) messing about too much.
  Tuning changes should be noted here.
  The testing next version may be activated by CBC_NEXT_VERSION
  This applies to OsiClp, Clp etc
  Version 1.00.01 November 24 2005
  Added several classes for advanced users.  This can't affect code (if you don't use it)
  Made some tiny changes (for N way branching) which should not change anything.
  CbcNWay object class - for N way branching this also allows use of CbcConsequence class.
  CbcBranchAllDifferent object class - for branching on general integer variables
  to stop them having same value so branches are x >= y+1 and x <= y-1.
  Added two new Cgl classes - CglAllDifferent which does column fixing (too slowly)
  and CglStored which just has a list of cuts which can be activated.
  Modified preprocess option to SOS
  Version 1.00.02 December 9 2005
  Added use of CbcStrategy to do clean preprocessing
  Added use of referenceSolver for cleaner repetition of Cbc
  Version 1.01.00 February 2 2006
  Added first try at Ampl interface
  Made dummy program so real main can be called from other programs
*/
