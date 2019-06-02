/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdio>
#include <string>
#include <iostream>
#include <iomanip>

#include "CoinTime.hpp"
#include "CoinFileIO.hpp"
#include "CbcModel.hpp"
#include "CbcHeuristic.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcBranchCut.hpp"
#include "CglProbing.hpp"
#include "OsiClpSolverInterface.hpp"
#include "ClpFactorization.hpp"
#include "OsiRowCutDebugger.hpp"
//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif

//#############################################################################

// Display message on stdout and stderr
void testingMessage(const char *const msg)
{
  std::cout << msg;
  //cout <<endl <<"*****************************************"
  //     <<endl <<msg <<endl;
}

//#############################################################################

static inline bool CbcTestFile(const std::string name)
{
  FILE *fp = fopen(name.c_str(), "r");
  if (fp) {
    fclose(fp);
    return true;
  }
  return false;
}

//#############################################################################

bool CbcTestMpsFile(std::string &fname)
{
  if (CbcTestFile(fname)) {
    return true;
  }
  if (CbcTestFile(fname + ".mps")) {
    fname += ".mps";
    return true;
  }
  if (CbcTestFile(fname + ".MPS")) {
    fname += ".MPS";
    return true;
  }
  if (CoinFileInput::haveGzipSupport()) {
    if (CbcTestFile(fname + ".gz")) {
      return true;
    }
    if (CbcTestFile(fname + ".mps.gz")) {
      fname += ".mps";
      return true;
    }
    if (CbcTestFile(fname + ".MPS.gz")) {
      fname += ".MPS";
      return true;
    }
    if (CbcTestFile(fname + ".MPS.GZ")) {
      fname += ".MPS";
      return true;
    }
  }
  if (CoinFileInput::haveBzip2Support()) {
    if (CbcTestFile(fname + ".bz2")) {
      return true;
    }
    if (CbcTestFile(fname + ".mps.bz2")) {
      fname += ".mps";
      return true;
    }
    if (CbcTestFile(fname + ".MPS.bz2")) {
      fname += ".MPS";
      return true;
    }
    if (CbcTestFile(fname + ".MPS.BZ2")) {
      fname += ".MPS";
      return true;
    }
  }
  return false;
}
//#############################################################################
/*
  jjf: testSwitch -2 unitTest, -1 normal (==2)

  MiplibTest might be more appropriate.

  TestSwitch and stuff[6] together control how much of miplib is executed:
    For testSwitch set to:
      -2: solve p0033 and p0201 only (the unit test)
      -1: solve miplib sets #0 and #1
       0: solve nothing
       k: execute sets j:k, where j is determined by the value of stuff[6]
  The last parameter of PUSH_MPS specifies the test set membership.

  For -miplib, -extra2 sets testSwitch, -extra3 sets stuff[6]. The command
    cbc -extra2 -2 -miplib
  will execute the unit test on the miplib directory.

  dirMiplib should end in the directory separator character for the platform.

  If you want to activate the row cut debugger for a given problem, change the
  last parameter of the PUSH_MPS macro for the problem to true.

  Returns 0 if all goes well, -1 if the Miplib directory is missing, otherwise
     100*(number with bad objective)+(number that exceeded node limit)
*/
int CbcClpUnitTest(const CbcModel &saveModel, const std::string &dirMiplib,
  int testSwitch, const double *stuff)
{
  // Stop Windows popup
  WindowsErrorPopupBlocker();
  unsigned int m;

  // Do an existence check.
  std::string test1 = dirMiplib + "p0033";
  bool doTest = CbcTestMpsFile(test1);
  if (!doTest) {
    std::cout
      << "Not doing miplib run as can't find mps files." << std::endl
      << "Perhaps you're trying to read gzipped (.gz) files without libz?"
      << std::endl;
    return (0);
  }
  int dfltPrecision = static_cast< int >(std::cout.precision());
  /*
  Set the range of problems to be tested. testSwitch = -2 is special and is
  picked up below.
*/
  int loSet = 0;
  int hiSet = 0;
  if (testSwitch == -1) {
    loSet = 0;
    hiSet = 1;
  } else if (testSwitch >= 0) {
    loSet = static_cast< int >(stuff[6]);
    hiSet = testSwitch;
    std::cout
      << "Solving miplib problems in sets " << loSet
      << ":" << hiSet << "." << std::endl;
  }
  /*
  Vectors to hold test problem names and characteristics.
*/
  std::vector< std::string > mpsName;
  std::vector< int > nRows;
  std::vector< int > nCols;
  std::vector< double > objValueC;
  std::vector< double > objValue;
  std::vector< int > testSet;
  std::vector< int > rowCutDebugger;
/*
  A macro to make the vector creation marginally readable. Parameters are
  name, rows, columns, integer objective, continuous objective, set ID,
  row cut debugger

  To enable the row cut debugger for a given problem, change the last
  parameter to true. Don't forget to turn it off before committing changes!
*/
#define PUSH_MPS(zz_mpsName_zz,                              \
  zz_nRows_zz, zz_nCols_zz, zz_objValue_zz, zz_objValueC_zz, \
  zz_testSet_zz, zz_rcDbg_zz)                                \
  mpsName.push_back(zz_mpsName_zz);                          \
  nRows.push_back(zz_nRows_zz);                              \
  nCols.push_back(zz_nCols_zz);                              \
  objValueC.push_back(zz_objValueC_zz);                      \
  testSet.push_back(zz_testSet_zz);                          \
  objValue.push_back(zz_objValue_zz);                        \
  rowCutDebugger.push_back(zz_rcDbg_zz);
  /*
  Push the miplib problems. Except for -2 (unitTest), push all, even if we're
  not going to do all of them.
*/
  if (testSwitch == -2) {
    PUSH_MPS("p0033", 16, 33, 3089, 2520.57, 0, false);
    PUSH_MPS("p0201", 133, 201, 7615, 6875.0, 0, false);
    // PUSH_MPS("flugpl", 18, 18, 1201500, 1167185.7, 0, false);
  } else {
/*
  Load up the problem vector. Note that the row counts here include the
  objective function.
*/
#if 1
    PUSH_MPS("10teams", 230, 2025, 924, 917, 1, false);
    PUSH_MPS("air03", 124, 10757, 340160, 338864.25, 0, false);
    PUSH_MPS("air04", 823, 8904, 56137, 55535.436, 2, false);
    PUSH_MPS("air05", 426, 7195, 26374, 25877.609, 2, false);
    PUSH_MPS("arki001", 1048, 1388, 7580813.0459, 7579599.80787, 7, false);
    PUSH_MPS("bell3a", 123, 133, 878430.32, 862578.64, 0, false);
    PUSH_MPS("bell5", 91, 104, 8966406.49, 8608417.95, 1, false);
    PUSH_MPS("blend2", 274, 353, 7.598985, 6.9156751140, 0, false);
    PUSH_MPS("cap6000", 2176, 6000, -2451377, -2451537.325, 1, false);
    PUSH_MPS("dano3mip", 3202, 13873, 728.1111, 576.23162474, 7, false);
    PUSH_MPS("danoint", 664, 521, 65.67, 62.637280418, 6, false);
    PUSH_MPS("dcmulti", 290, 548, 188182, 183975.5397, 0, false);
    PUSH_MPS("dsbmip", 1182, 1886, -305.19817501, -305.19817501, 0, false);
    PUSH_MPS("egout", 98, 141, 568.101, 149.589, 0, false);
    PUSH_MPS("enigma", 21, 100, 0.0, 0.0, 0, false);
    PUSH_MPS("fast0507", 507, 63009, 174, 172.14556668, 5, false);
    PUSH_MPS("fiber", 363, 1298, 405935.18000, 156082.51759, 0, false);
    PUSH_MPS("fixnet6", 478, 878, 3983, 1200.88, 1, false);
    PUSH_MPS("flugpl", 18, 18, 1201500, 1167185.7, 0, false);
    PUSH_MPS("gen", 780, 870, 112313, 112130.0, 0, false);
    PUSH_MPS("gesa2", 1392, 1224, 25779856.372, 25476489.678, 1, false);
    PUSH_MPS("gesa2_o", 1248, 1224, 25779856.372, 25476489.678, 1, false);
    PUSH_MPS("gesa3", 1368, 1152, 27991042.648, 27833632.451, 0, false);
    PUSH_MPS("gesa3_o", 1224, 1152, 27991042.648, 27833632.451, 0, false);
    PUSH_MPS("gt2", 29, 188, 21166.000, 13460.233074, 0, false);
    PUSH_MPS("harp2", 112, 2993, -73899798.00, -74353341.502, 6, false);
    PUSH_MPS("khb05250", 101, 1350, 106940226, 95919464.0, 0, false);
    PUSH_MPS("l152lav", 97, 1989, 4722, 4656.36, 1, false);
    PUSH_MPS("lseu", 28, 89, 1120, 834.68, 0, false);
    PUSH_MPS("mas74", 13, 151, 11801.18573, 10482.79528, 3, false);
    PUSH_MPS("mas76", 12, 151, 40005.05414, 38893.9036, 2, false);
    PUSH_MPS("misc03", 96, 160, 3360, 1910., 0, false);
    PUSH_MPS("misc06", 820, 1808, 12850.8607, 12841.6, 0, false);
    PUSH_MPS("misc07", 212, 260, 2810, 1415.0, 1, false);
    PUSH_MPS("mitre", 2054, 10724, 115155, 114740.5184, 1, false);
    PUSH_MPS("mkc", 3411, 5325, -553.75, -611.85, 7, false); // suboptimal
    PUSH_MPS("mod008", 6, 319, 307, 290.9, 0, false);
    PUSH_MPS("mod010", 146, 2655, 6548, 6532.08, 0, false);
    PUSH_MPS("mod011", 4480, 10958, -54558535, -62121982.55, 2, false);
    PUSH_MPS("modglob", 291, 422, 20740508, 20430947., 2, false);
    PUSH_MPS("noswot", 182, 128, -43, -43.0, 6, false);
    PUSH_MPS("nw04", 36, 87482, 16862, 16310.66667, 1, false);
    PUSH_MPS("p0033", 16, 33, 3089, 2520.57, 0, false);
    PUSH_MPS("p0201", 133, 201, 7615, 6875.0, 0, false);
    PUSH_MPS("p0282", 241, 282, 258411, 176867.50, 0, false);
    PUSH_MPS("p0548", 176, 548, 8691, 315.29, 0, false);
    PUSH_MPS("p2756", 755, 2756, 3124, 2688.75, 0, false);
    PUSH_MPS("pk1", 45, 86, 11.0, 0.0, 2, false);
    PUSH_MPS("pp08a", 136, 240, 7350.0, 2748.3452381, 1, false);
    PUSH_MPS("pp08aCUTS", 246, 240, 7350.0, 5480.6061563, 1, false);
    PUSH_MPS("qiu", 1192, 840, -132.873137, -931.638857, 3, false);
    PUSH_MPS("qnet1", 503, 1541, 16029.692681, 14274.102667, 0, false);
    PUSH_MPS("qnet1_o", 456, 1541, 16029.692681, 12095.571667, 0, false);
    PUSH_MPS("rentacar", 6803, 9557, 30356761, 28806137.644, 0, false);
    PUSH_MPS("rgn", 24, 180, 82.1999, 48.7999, 0, false);
    PUSH_MPS("rout", 291, 556, 1077.56, 981.86428571, 3, false);
    PUSH_MPS("set1ch", 492, 712, 54537.75, 32007.73, 5, false);
    PUSH_MPS("seymour", 4944, 1372, 423, 403.84647413, 7, false);
    PUSH_MPS("seymour_1", 4944, 1372, 410.76370, 403.84647413, 5, false);
    PUSH_MPS("stein27", 118, 27, 18, 13.0, 0, false);
    PUSH_MPS("stein45", 331, 45, 30, 22.0, 1, false);
    PUSH_MPS("swath", 884, 6805, 497.603, 334.4968581, 7, false);
    PUSH_MPS("vpm1", 234, 378, 20, 15.4167, 0, false);
    PUSH_MPS("vpm2", 234, 378, 13.75, 9.8892645972, 0, false);
#endif
  }
#undef PUSH_MPS

  /*
  Normally the problems are executed in order. Define RANDOM_ORDER below to
  randomize.

  #define RANDOM_ORDER
*/
  int which[100];
  int nLoop = static_cast< int >(mpsName.size());
  assert(nLoop <= 100);
  for (int i = 0; i < nLoop; i++)
    which[i] = i;

#ifdef RANDOM_ORDER
  unsigned int iTime = static_cast< unsigned int >(CoinGetTimeOfDay() - 1.256e9);
  std::cout << "Time (seed) " << iTime << "." << std::endl;
  double sort[100];
  CoinDrand48(true, iTime);
  for (int i = 0; i < nLoop; i++)
    sort[i] = CoinDrand48();
  CoinSort_2(sort, sort + nLoop, which);
#endif

  int problemCnt = 0;
  for (m = 0; m < mpsName.size(); m++) {
    int setID = testSet[m];
    if (loSet <= setID && setID <= hiSet)
      problemCnt++;
  }

  int numberFailures = 0;
  int numberAttempts = 0;
  int numProbSolved = 0;
  double timeTaken = 0.0;

//#define CLP_FACTORIZATION_INSTRUMENT
#ifdef CLP_FACTORIZATION_INSTRUMENT
  double timeTakenFac = 0.0;
#endif

  /*
  Open the main loop to step through the MPS problems.
*/
  for (unsigned int mw = 0; mw < mpsName.size(); mw++) {
    m = which[mw];
    int setID = testSet[m];
    // Skip if problem is not in specified problem set(s)
    if (!(loSet <= setID && setID <= hiSet))
      continue;

    numberAttempts++;
    std::cout << "  processing mps file: " << mpsName[m]
              << " (" << numberAttempts << " out of "
              << problemCnt << ")" << std::endl;
    /*
  Stage 1: Read the MPS and make sure the size of the constraint matrix
	   is correct.
*/
    CbcModel *model = new CbcModel(saveModel);

    std::string fn = dirMiplib + mpsName[m];
    if (!CbcTestMpsFile(fn)) {
      std::cout << "ERROR: Cannot find MPS file " << fn << "." << std::endl;
      continue;
    }
    model->solver()->readMps(fn.c_str(), "");
    assert(model->getNumRows() == nRows[m]);
    assert(model->getNumCols() == nCols[m]);

    // Careful! We're initialising for the benefit of other code.
    CoinDrand48(true, 1234567);
    //printf("RAND1 %g %g\n",CoinDrand48(true,1234567),model->randomNumberGenerator()->randomDouble());
    //printf("RAND1 %g\n",CoinDrand48(true,1234567));

    // Higher limits for the serious problems.
    int testMaximumNodes = 200000;
    if (hiSet > 1)
      testMaximumNodes = 20000000;
    if (model->getMaximumNodes() > testMaximumNodes) {
      model->setMaximumNodes(testMaximumNodes);
    }
    /*
  Stage 2: Call solver to solve the problem.
*/

#ifdef CLP_FACTORIZATION_INSTRUMENT
    extern double factorization_instrument(int type);
    double facTime1 = factorization_instrument(0);
    std::cout
      << "Factorization - initial solve " << facTime1 << " seconds."
      << std::endl;
    timeTakenFac += facTime1;
#endif

    double startTime = CoinCpuTime() + CoinCpuTimeJustChildren();

    // Setup specific to clp
    OsiClpSolverInterface *siClp = dynamic_cast< OsiClpSolverInterface * >(model->solver());
    ClpSimplex *modelC = NULL;
    if (siClp) {
      modelC = siClp->getModelPtr();
      ClpMatrixBase *matrix = modelC->clpMatrix();
      ClpPackedMatrix *clpMatrix = dynamic_cast< ClpPackedMatrix * >(matrix);
      if (stuff && stuff[9] && clpMatrix) {
        // vector matrix!
        clpMatrix->makeSpecialColumnCopy();
      }

#ifdef JJF_ZERO
      if (clpMatrix) {
        int numberRows = clpMatrix->getNumRows();
        int numberColumns = clpMatrix->getNumCols();
        double *elements = clpMatrix->getMutableElements();
        const int *row = clpMatrix->getIndices();
        const CoinBigIndex *columnStart = clpMatrix->getVectorStarts();
        const int *columnLength = clpMatrix->getVectorLengths();
        double *smallest = new double[numberRows];
        double *largest = new double[numberRows];
        char *flag = new char[numberRows];
        CoinZeroN(flag, numberRows);
        for (int i = 0; i < numberRows; i++) {
          smallest[i] = COIN_DBL_MAX;
          largest[i] = 0.0;
        }
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          bool isInteger = modelC->isInteger(iColumn);
          CoinBigIndex j;
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            double value = fabs(elements[j]);
            if (!isInteger)
              flag[iRow] = 1;
            smallest[iRow] = CoinMin(smallest[iRow], value);
            largest[iRow] = CoinMax(largest[iRow], value);
          }
        }
        double *rowLower = modelC->rowLower();
        double *rowUpper = modelC->rowUpper();
        bool changed = false;
        for (int i = 0; i < numberRows; i++) {
          if (flag[i] && smallest[i] > 10.0 && false) {
            smallest[i] = 1.0 / smallest[i];
            if (rowLower[i] > -1.0e20)
              rowLower[i] *= smallest[i];
            if (rowUpper[i] < 1.0e20)
              rowUpper[i] *= smallest[i];
            changed = true;
          } else {
            smallest[i] = 0.0;
          }
        }
        if (changed) {
          printf("SCALED\n");
          for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
            CoinBigIndex j;
            for (j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              if (smallest[iRow])
                elements[j] *= smallest[iRow];
            }
          }
        }
        delete[] smallest;
        delete[] largest;
        delete[] flag;
      }
#endif // JJF_ZERO

      model->checkModel();
      modelC->tightenPrimalBounds(0.0, 0, true);
      model->initialSolve();
      if (modelC->dualBound() == 1.0e10) {
        // user did not set - so modify
        // get largest scaled away from bound
        ClpSimplex temp = *modelC;
        temp.dual(0, 7);
        double largestScaled = 1.0e-12;
        double largest = 1.0e-12;
        int numberRows = temp.numberRows();
        const double *rowPrimal = temp.primalRowSolution();
        const double *rowLower = temp.rowLower();
        const double *rowUpper = temp.rowUpper();
        const double *rowScale = temp.rowScale();
        int iRow;
        for (iRow = 0; iRow < numberRows; iRow++) {
          double value = rowPrimal[iRow];
          double above = value - rowLower[iRow];
          double below = rowUpper[iRow] - value;
          if (above < 1.0e12) {
            largest = CoinMax(largest, above);
          }
          if (below < 1.0e12) {
            largest = CoinMax(largest, below);
          }
          if (rowScale) {
            double multiplier = rowScale[iRow];
            above *= multiplier;
            below *= multiplier;
          }
          if (above < 1.0e12) {
            largestScaled = CoinMax(largestScaled, above);
          }
          if (below < 1.0e12) {
            largestScaled = CoinMax(largestScaled, below);
          }
        }

        int numberColumns = temp.numberColumns();
        const double *columnPrimal = temp.primalColumnSolution();
        const double *columnLower = temp.columnLower();
        const double *columnUpper = temp.columnUpper();
        const double *columnScale = temp.columnScale();
        int iColumn;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          double value = columnPrimal[iColumn];
          double above = value - columnLower[iColumn];
          double below = columnUpper[iColumn] - value;
          if (above < 1.0e12) {
            largest = CoinMax(largest, above);
          }
          if (below < 1.0e12) {
            largest = CoinMax(largest, below);
          }
          if (columnScale) {
            double multiplier = 1.0 / columnScale[iColumn];
            above *= multiplier;
            below *= multiplier;
          }
          if (above < 1.0e12) {
            largestScaled = CoinMax(largestScaled, above);
          }
          if (below < 1.0e12) {
            largestScaled = CoinMax(largestScaled, below);
          }
        }
        std::cout << "Largest (scaled) away from bound " << largestScaled
                  << " unscaled " << largest << std::endl;
#ifdef JJF_ZERO
        modelC->setDualBound(CoinMax(1.0001e8,
          CoinMin(1000.0 * largestScaled, 1.00001e10)));
#else
        modelC->setDualBound(CoinMax(1.0001e9,
          CoinMin(1000.0 * largestScaled, 1.0001e10)));
#endif
      }
    } // end clp-specific setup
    /*
  Cut passes: For small models (n < 500) always do 100 passes, if possible
  (-100). For larger models, use minimum drop to stop (100, 20).
*/
    model->setMinimumDrop(CoinMin(5.0e-2,
      fabs(model->getMinimizationObjValue()) * 1.0e-3 + 1.0e-4));
    if (CoinAbs(model->getMaximumCutPassesAtRoot()) <= 100) {
      if (model->getNumCols() < 500) {
        model->setMaximumCutPassesAtRoot(-100);
      } else if (model->getNumCols() < 5000) {
        model->setMaximumCutPassesAtRoot(100);
      } else {
        model->setMaximumCutPassesAtRoot(20);
      }
    }
    // If defaults then increase trust for small models
    if (model->numberStrong() == 5 && model->numberBeforeTrust() == 10) {
      int numberColumns = model->getNumCols();
      if (numberColumns <= 50) {
        model->setNumberBeforeTrust(1000);
      } else if (numberColumns <= 100) {
        model->setNumberBeforeTrust(100);
      } else if (numberColumns <= 300) {
        model->setNumberBeforeTrust(50);
      }
    }
    //if (model->getNumCols()>=500) {
    // switch off Clp stuff
    //model->setFastNodeDepth(-1);
    //}
    /*
  Activate the row cut debugger, if requested.
*/
    if (rowCutDebugger[m] == true) {
      std::string probName;
      model->solver()->getStrParam(OsiProbName, probName);
      model->solver()->activateRowCutDebugger(probName.c_str());
      if (model->solver()->getRowCutDebugger())
        std::cout << "Row cut debugger activated for ";
      else
        std::cout << "Failed to activate row cut debugger for ";
      std::cout << mpsName[m] << "." << std::endl;
    }
    setCutAndHeuristicOptions(*model);
    /*
  More clp-specific setup.
*/
    if (siClp) {
#ifdef CLP_MULTIPLE_FACTORIZATIONS
      if (!modelC->factorization()->isDenseOrSmall()) {
        int denseCode = stuff ? static_cast< int >(stuff[4]) : -1;
        int smallCode = stuff ? static_cast< int >(stuff[10]) : -1;
        if (stuff && stuff[8] >= 1) {
          if (denseCode < 0)
            denseCode = 40;
          if (smallCode < 0)
            smallCode = 40;
        }
        if (denseCode > 0)
          modelC->factorization()->setGoDenseThreshold(denseCode);
        if (smallCode > 0)
          modelC->factorization()->setGoSmallThreshold(smallCode);
        if (denseCode >= modelC->numberRows()) {
          //printf("problem going dense\n");
          //modelC->factorization()->goDenseOrSmall(modelC->numberRows());
        }
      }
#endif
      if (stuff && stuff[8] >= 1) {
        printf("Fast node size Columns %d rows %d - depth %d\n",
          modelC->numberColumns(), modelC->numberRows(),
          model->fastNodeDepth());
        if (modelC->numberColumns() + modelC->numberRows() <= 10000 && model->fastNodeDepth() == -1)
          model->setFastNodeDepth(-10 /*-9*/);
      }
    }
#ifdef CONFLICT_CUTS
    {
      model->setCutoffAsConstraint(true); // very slow on bell5 ??
      int moreOptions = model->moreSpecialOptions();
      model->setMoreSpecialOptions(moreOptions | 4194304);
    }
#endif
    /*
  Finally, the actual call to solve the MIP with branch-and-cut.
*/
    model->branchAndBound();

#ifdef CLP_FACTORIZATION_INSTRUMENT
    double facTime = factorization_instrument(0);
    std::cout << "Factorization " << facTime << " seconds." << std::endl,
      timeTakenFac += facTime;
#endif

    /*
  Stage 3: Do the statistics and check the answer.
*/
    double timeOfSolution = CoinCpuTime() + CoinCpuTimeJustChildren() - startTime;
    std::cout
      << "Cuts at root node changed objective from "
      << model->getContinuousObjective() << " to "
      << model->rootObjectiveAfterCuts() << std::endl;
    int numberGenerators = model->numberCutGenerators();
    for (int iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
      CbcCutGenerator *generator = model->cutGenerator(iGenerator);
#ifdef CLIQUE_ANALYSIS
#ifndef CLP_INVESTIGATE
      CglImplication *implication = dynamic_cast< CglImplication * >(generator->generator());
      if (implication)
        continue;
#endif
#endif
      std::cout
        << generator->cutGeneratorName() << " was tried "
        << generator->numberTimesEntered() << " times and created "
        << generator->numberCutsInTotal() << " cuts of which "
        << generator->numberCutsActive()
        << " were active after adding rounds of cuts";
      if (generator->timing())
        std::cout << " (" << generator->timeInCutGenerator() << " seconds)";
      std::cout << "." << std::endl;
    }
    std::cout
      << model->getNumberHeuristicSolutions()
      << " solutions found by heuristics." << std::endl;
    int numberHeuristics = model->numberHeuristics();
    for (int iHeuristic = 0; iHeuristic < numberHeuristics; iHeuristic++) {
      CbcHeuristic *heuristic = model->heuristic(iHeuristic);
      if (heuristic->numRuns()) {
        std::cout
          << heuristic->heuristicName() << " was tried "
          << heuristic->numRuns() << " times out of "
          << heuristic->numCouldRun() << " and created "
          << heuristic->numberSolutionsFound() << " solutions." << std::endl;
      }
    }
    /*
  Check for the correct answer.
*/
    if (!model->status()) {

      double objActual = model->getObjValue();
      double objExpect = objValue[m];
      double tolerance = CoinMin(fabs(objActual), fabs(objExpect));
      tolerance = CoinMax(1.0e-5, 1.0e-5 * tolerance);
      //CoinRelFltEq eq(1.0e-3) ;

      std::cout
        << "cbc_clp (" << mpsName[m] << ") "
        << std::setprecision(10) << objActual;
      if (fabs(objActual - objExpect) < tolerance) {
        std::cout << std::setprecision(dfltPrecision) << "; okay";
        numProbSolved++;
      } else {
        std::cout
          << " != " << objExpect << std::setprecision(dfltPrecision)
          << "; error = " << fabs(objExpect - objActual);
        numberFailures++;
        //#ifdef COIN_DEVELOP
        //abort();
        //#endif
      }
    } else {
      std::cout
        << "cbc_clp (" << mpsName[m] << ") status not optimal; "
        << "assuming too many nodes";
    }
    timeTaken += timeOfSolution;
    std::cout
      << " -- (" << model->getNodeCount() << " n / "
      << model->getIterationCount() << " i / "
      << timeOfSolution << " s) (subtotal " << timeTaken << " seconds)"
      << std::endl;
    delete model;
  }
  /*
  End main loop on MPS problems. Print a summary and calculate the return
  value.
*/
  int returnCode = 0;
  std::cout
    << "cbc_clp solved " << numProbSolved << " out of " << numberAttempts;
  int numberOnNodes = numberAttempts - numProbSolved - numberFailures;
  if (numberFailures || numberOnNodes) {
    if (numberOnNodes) {
      std::cout << " (" << numberOnNodes << " stopped on nodes)";
      returnCode = numberOnNodes;
    }
    if (numberFailures) {
      std::cout << " (" << numberFailures << " gave bad answer!)";
      returnCode += 100 * numberFailures;
    }
  }
  std::cout
    << " and took " << timeTaken << " seconds." << std::endl;

  if (testSwitch == -2) {
    if (numberFailures || numberOnNodes) {
      std::cout << "****** Unit Test failed." << std::endl;
      std::cerr << "****** Unit Test failed." << std::endl;
    } else {
      std::cerr << "****** Unit Test succeeded." << std::endl;
    }
  }
#ifdef CLP_FACTORIZATION_INSTRUMENT
  std::cout
    << "Total factorization time " << timeTakenFac << "seconds." << std::endl;
#endif
  return (returnCode);
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
