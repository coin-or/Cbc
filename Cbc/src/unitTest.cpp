/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Test individual classes or groups of classes

#include "CbcConfig.h"
#include "CoinPragma.hpp"

#include <cstdlib>
#include <cassert>
#include <vector>
#include <iostream>
#include <cstdio>

#include "CoinHelperFunctions.hpp"
#include "CoinFileIO.hpp"

#ifdef COIN_HAS_CBC
#include "OsiCbcSolverInterface.hpp"
#endif
#ifdef COIN_HAS_OSL
#include "OsiOslSolverInterface.hpp"
#endif
#ifdef COIN_HAS_SPX
#include "OsiSpxSolverInterface.hpp"
#endif
#ifdef COIN_HAS_DYLP
#include "OsiDylpSolverInterface.hpp"
#endif
#ifdef COIN_HAS_GLPK
#include "OsiGlpkSolverInterface.hpp"
#endif
#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif
#ifdef NDEBUG
#undef NDEBUG
#endif

#include "CoinTime.hpp"
// Function Prototypes. Function definitions is in this file.
void testingMessage(const char *const msg);

#ifdef COIN_HAS_CBC
void CbcMiplibTest(const std::vector< OsiCbcSolverInterface * > &vecEmptySiP,
  const std::string &mpsDir)
{
  int i;
  unsigned int m;
  // See if files exist
  FILE *fp;
  bool doTest = false;
  std::string test1 = mpsDir + "p0033";
  fp = fopen(test1.c_str(), "r");
  if (fp) {
    doTest = true;
    fclose(fp);
  }
  if (!doTest && CoinFileInput::haveGzipSupport()) {
    test1 += ".gz";
    fp = fopen(test1.c_str(), "r");
    if (fp) {
      doTest = true;
      fclose(fp);
    }
  }
  if (!doTest)
    return;
  /*
      Vectors to hold test problem names and characteristics. The objective value
      after optimization (objValue) must agree to the specified tolerance
      (objValueTol).
    */
  std::vector< std::string > mpsName;
  std::vector< int > nRows;
  std::vector< int > nCols;
  std::vector< double > objValueC;
  std::vector< double > objValue;
  std::vector< int > strategy;
  /*
      And a macro to make the vector creation marginally readable.
    */
#define PUSH_MPS(zz_mpsName_zz,                              \
  zz_nRows_zz, zz_nCols_zz, zz_objValue_zz, zz_objValueC_zz, \
  zz_strategy_zz)                                            \
  mpsName.push_back(zz_mpsName_zz);                          \
  nRows.push_back(zz_nRows_zz);                              \
  nCols.push_back(zz_nCols_zz);                              \
  objValueC.push_back(zz_objValueC_zz);                      \
  strategy.push_back(zz_strategy_zz);                        \
  objValue.push_back(zz_objValue_zz);

  /*
      Load up the problem vector. Note that the row counts here include the
      objective function.

      Set HOWMANY to 0 for no test, 1 for some, 2 for many, 3 for all.
    */
#define HOWMANY 1
#if HOWMANY
#if HOWMANY > 1
  PUSH_MPS("10teams", 230, 2025, 924, 917, 7)
#endif
  PUSH_MPS("air03", 124, 10757, 340160, 338864.25, 7)
#if HOWMANY == 3
  PUSH_MPS("air04", 823, 8904, 56137, 55535.436, 8)
  PUSH_MPS("air05", 426, 7195, 26374, 25877.609, 8)
#endif
  //    PUSH_MPS("arki001",1048,1388,7580813.0459,7579599.80787,7)
  PUSH_MPS("bell3a", 123, 133, 878430.32, 862578.64, 7)
#if HOWMANY > 1
  PUSH_MPS("bell5", 91, 104, 8966406.49, 8608417.95, 7)
#endif
  PUSH_MPS("blend2", 274, 353, 7.598985, 6.9156751140, 7)
#if HOWMANY > 1
  PUSH_MPS("cap6000", 2176, 6000, -2451377, -2451537.325, 7)
#endif
  //    PUSH_MPS("dano3mip",3202,13873,728.1111,576.23162474,7)
  //    PUSH_MPS("danoint",664,521,65.67,62.637280418,7)
  PUSH_MPS("dcmulti", 290, 548, 188182, 183975.5397, 7)
  PUSH_MPS("dsbmip", 1182, 1886, -305.19817501, -305.19817501, 7)
  PUSH_MPS("egout", 98, 141, 568.101, 149.589, 7)
  PUSH_MPS("enigma", 21, 100, 0.0, 0.0, 7)
#if HOWMANY == 3
  PUSH_MPS("fast0507", 507, 63009, 174, 172.14556668, 7)
#endif
  PUSH_MPS("fiber", 363, 1298, 405935.18000, 156082.51759, 7)
#if HOWMANY > 1
  PUSH_MPS("fixnet6", 478, 878, 3983, 1200.88, 7)
#endif
  PUSH_MPS("flugpl", 18, 18, 1201500, 1167185.7, 7)
  PUSH_MPS("gen", 780, 870, 112313, 112130.0, 7)
#if HOWMANY > 1
  PUSH_MPS("gesa2", 1392, 1224, 25779856.372, 25476489.678, 7)
  PUSH_MPS("gesa2_o", 1248, 1224, 25779856.372, 25476489.678, 7)
#endif
  PUSH_MPS("gesa3", 1368, 1152, 27991042.648, 27833632.451, 7)
  PUSH_MPS("gesa3_o", 1224, 1152, 27991042.648, 27833632.451, 7)
  PUSH_MPS("gt2", 29, 188, 21166.000, 13460.233074, 7)
#if HOWMANY == 3
  PUSH_MPS("harp2", 112, 2993, -73899798.00, -74353341.502, 7)
#endif
  PUSH_MPS("khb05250", 101, 1350, 106940226, 95919464.0, 7)
#if HOWMANY > 1
  PUSH_MPS("l152lav", 97, 1989, 4722, 4656.36, 7)
#endif
  PUSH_MPS("lseu", 28, 89, 1120, 834.68, 7)
  PUSH_MPS("misc03", 96, 160, 3360, 1910., 7)
  PUSH_MPS("misc06", 820, 1808, 12850.8607, 12841.6, 7)
#if HOWMANY > 1
  PUSH_MPS("misc07", 212, 260, 2810, 1415.0, 7)
  PUSH_MPS("mitre", 2054, 10724, 115155, 114740.5184, 7)
#endif
  PUSH_MPS("mod008", 6, 319, 307, 290.9, 7)
  PUSH_MPS("mod010", 146, 2655, 6548, 6532.08, 7)
#if HOWMANY == 3
  PUSH_MPS("mod011", 4480, 10958, -54558535, -62121982.55, 7)
  PUSH_MPS("modglob", 291, 422, 20740508, 20430947., 7)
  PUSH_MPS("noswot", 182, 128, -43, -43.0, 7)
#endif
#if HOWMANY > 1
  PUSH_MPS("nw04", 36, 87482, 16862, 16310.66667, 7)
#endif
  PUSH_MPS("p0033", 16, 33, 3089, 2520.57, 7)
  PUSH_MPS("p0201", 133, 201, 7615, 6875.0, 7)
  PUSH_MPS("p0282", 241, 282, 258411, 176867.50, 7)
  PUSH_MPS("p0548", 176, 548, 8691, 315.29, 7)
  PUSH_MPS("p2756", 755, 2756, 3124, 2688.75, 7)
#if HOWMANY == 3
  PUSH_MPS("pk1", 45, 86, 11.0, 0.0, 7)
#endif
#if HOWMANY > 1
  PUSH_MPS("pp08a", 136, 240, 7350.0, 2748.3452381, 7)
  PUSH_MPS("pp08aCUTS", 246, 240, 7350.0, 5480.6061563, 7)
#endif
#if HOWMANY == 3
  PUSH_MPS("qiu", 1192, 840, -132.873137, -931.638857, 7)
#endif
  PUSH_MPS("qnet1", 503, 1541, 16029.692681, 14274.102667, 7)
  PUSH_MPS("qnet1_o", 456, 1541, 16029.692681, 12095.571667, 7)
  PUSH_MPS("rentacar", 6803, 9557, 30356761, 28806137.644, 7)
  PUSH_MPS("rgn", 24, 180, 82.1999, 48.7999, 7)
#if HOWMANY == 3
  PUSH_MPS("rout", 291, 556, 1077.56, 981.86428571, 7)
  PUSH_MPS("set1ch", 492, 712, 54537.75, 32007.73, 7)
#endif
  //    PUSH_MPS("seymour",4944,1372,423,403.84647413,7)
  PUSH_MPS("stein27", 118, 27, 18, 13.0, 7)
#if HOWMANY > 1
  PUSH_MPS("stein45", 331, 45, 30, 22.0, 7)
#endif
  PUSH_MPS("vpm1", 234, 378, 20, 15.4167, 7)
  PUSH_MPS("vpm2", 234, 378, 13.75, 9.8892645972, 7)
#endif
#undef PUSH_MPS

  /*
      Create a vector of solver interfaces that we can use to run the test
      problems. The strategy is to create a fresh clone of the `empty' solvers
      from vecEmptySiP for each problem, then proceed in stages: read the MPS
      file, solve the problem, check the solution. If there are multiple
      solvers in vecSiP, the results of each solver are compared with its
      neighbors in the vector.
    */
  int numberSolvers = vecEmptySiP.size();
  std::vector< OsiSolverInterface * > vecSiP(numberSolvers);

  // Create vector to store a name for each solver interface
  // and a count on the number of problems the solver interface solved.
  std::vector< std::string > siName;
  std::vector< int > numProbSolved;
  std::vector< double > timeTaken;
  for (i = 0; i < numberSolvers; i++) {
    std::string name;
    vecEmptySiP[i]->getStrParam(OsiSolverName, name);
    siName.push_back(name);
    numProbSolved.push_back(0);
    timeTaken.push_back(0.0);
  }

  /*
      Open the main loops. Outer loop steps through MPS problems, inner loop
      steps through solvers.
    */
  for (m = 0; m < mpsName.size(); m++) {
    std::cerr << "  processing mps file: " << mpsName[m]
              << " (" << m + 1 << " out of " << mpsName.size() << ")"
              << std::endl;
    for (i = vecSiP.size() - 1; i >= 0; --i) {
      vecSiP[i] = vecEmptySiP[i]->clone();
      /*
              Stage 1: Read the MPS file into the solver interface.

              As a basic check, make sure the size of the constraint matrix is correct.
            */

      std::string fn = mpsDir + mpsName[m];
      vecSiP[i]->readMps(fn.c_str(), "");

      vecSiP[i]->setObjSense(1.0);

      int nr = vecSiP[i]->getNumRows();
      int nc = vecSiP[i]->getNumCols();
      assert(nr == nRows[m]);
      assert(nc == nCols[m]);
      /*
              Stage 2: Call the solver to get a solution for the LP relaxation.
            */
      double startTime = CoinCpuTime();
      OsiCbcSolverInterface *integerSolver = dynamic_cast< OsiCbcSolverInterface * >(vecSiP[i]);
      assert(integerSolver);
      integerSolver->initialSolve();
      /*
              Stage 3: Call the solver to perform branch and cut.

              We call each solver, then check the return code and objective.
              Limits are 50000 nodes and one hour of time.
            */

      integerSolver->setMaximumNodes(50000);
      integerSolver->setMaximumSeconds(60 * 60);
      integerSolver->getModelPtr()->messageHandler()->setLogLevel(1);
      integerSolver->branchAndBound();

      double timeOfSolution = CoinCpuTime() - startTime;
      if (!integerSolver->status()) {
        double soln = integerSolver->getObjValue();
        CoinRelFltEq eq(1.0e-3);
        if (eq(soln, objValue[m])) {
          std::cerr
            << siName[i] << " "
            << soln << " = " << objValue[m] << " ; okay";
          numProbSolved[i]++;
        } else {
          std::cerr << siName[i] << " " << soln << " != " << objValue[m] << "; error=";
          std::cerr << fabs(objValue[m] - soln);
        }
      } else {
        std::cerr << "error; too many nodes";
      }
      std::cerr << " - took " << timeOfSolution << " seconds." << std::endl;
      timeTaken[i] += timeOfSolution;
      delete integerSolver;
    }
  }

  const int siName_size = siName.size();
  for (i = 0; i < siName_size; i++) {
    std::cerr
      << siName[i]
      << " solved "
      << numProbSolved[i]
      << " out of "
      << objValue.size()
      << " and took "
      << timeTaken[i]
      << " seconds."
      << std::endl;
  }
}
#endif // COIN_HAS_CBC

//----------------------------------------------------------------
// unitTest  [-miplibDir=V2]
//
// where:
//   -miplibDir: directory containing miplib files
//       Default value V2="./examples/miplib3"
//
// All parameters are optional.
//----------------------------------------------------------------

int mainTest(int argc, const char *argv[])
{
  int i;

  // define valid parameter keywords
  std::set< std::string > definedKeyWords;
  definedKeyWords.insert("-miplibDir");

  // Create a map of parameter keys and associated data
  std::map< std::string, std::string > parms;
  for (i = 1; i < argc; i++) {
    std::string parm(argv[i]);
    std::string key, value;
    unsigned int eqPos = parm.find('=');

    // Does parm contain and '='
    if (eqPos == std::string::npos) {
      //Parm does not contain '='
      key = parm;
    } else {
      key = parm.substr(0, eqPos);
      value = parm.substr(eqPos + 1);
    }

    // Is specifed key valid?
    if (definedKeyWords.find(key) == definedKeyWords.end()) {
      // invalid key word.
      // Write help text
      std::cerr << "Undefined parameter \"" << key << "\".\n";
      std::cerr << "Correct usage: \n";
      std::cerr << "  unitTest [-miplibDir=V2] \n";
      std::cerr << "  where:\n";
      std::cerr << "    -miplibDir: directory containing miplib files\n";
      std::cerr << "        Default value V2=\"./Data/miplib3\"\n";
      return 1;
    }
    parms[key] = value;
  }

  const char dirsep = CoinFindDirSeparator();

  // Set directory containing miplib data files.
  std::string miplibDir;
  if (parms.find("-miplibDir") != parms.end())
    miplibDir = parms["-miplibDir"] + dirsep;
  else
    miplibDir = dirsep == '/' ? "./Data/miplib3/" : ".\\Data\\miplib3\\";
#ifdef COIN_HAS_CBC

  {
    // Create vector of solver interfaces
    std::vector< OsiCbcSolverInterface * > vecSi;
    CbcStrategyDefault strategy(0);
#if COIN_HAS_OSL
    OsiSolverInterface *oslSi = new OsiOslSolverInterface;
    vecSi.push_back(new OsiCbcSolverInterface(oslSi, &strategy));
#endif
#if COIN_HAS_SPX
    OsiSolverInterface *spxSi = new OsiSpxSolverInterface;
    vecSi.push_back(new OsiCbcSolverInterface(spxSi, &strategy));
#endif
#if COIN_HAS_CLP
    OsiSolverInterface *clpSi = new OsiClpSolverInterface;
    /* Quiet, already! */
    clpSi->setHintParam(OsiDoReducePrint, true, OsiHintDo);
    vecSi.push_back(new OsiCbcSolverInterface(clpSi, &strategy));
#endif
#if COIN_HAS_DYLP
    OsiSolverInterface *dylpSi = new OsiDylpSolverInterface;
    vecSi.push_back(new OsiCbcSolverInterface(dylpSi, &strategy));
#endif
#if COIN_HAS_GLPK
    OsiSolverInterface *glpkSi = new OsiGlpkSolverInterface;
    vecSi.push_back(new OsiCbcSolverInterface(glpkSi, &strategy));
#endif

    testingMessage("Testing some miplib stuff\n");
    CbcMiplibTest(vecSi, miplibDir);

    unsigned int i;
    for (i = 0; i < vecSi.size(); i++)
      delete vecSi[i];
  }
#else // COIN_HAS_CBC
  std::cerr
    << "cbc has been built without OsiCbc support. To enable the -miplib\n"
    << "option, you must enable libOsiCbc in Makefile.location, then\n"
    << "execute the command `make clean cbc' to rebuild the cbc program."
    << std::endl;
#endif // COIN_HAS_CBC
  testingMessage("All tests completed successfully\n");
  return 0;
}

// Display message on stdout and stderr
void testingMessage(const char *const msg)
{
  std::cerr << msg << std::endl;
  // std::cout << msg << std::endl ;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
