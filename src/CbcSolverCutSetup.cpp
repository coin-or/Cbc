// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file CbcSolverCutSetup.cpp
    \brief Routine for installing cut generators on a CbcModel.
*/

#include "CbcConfig.h"

#include "CoinPragma.hpp"

#include "CbcCutGenerator.hpp"
#include "CbcModel.hpp"
#include "CbcParam.hpp"
#include "CbcParameters.hpp"
#include "CbcSolverCutSetup.hpp"

#include "CglBKClique.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglGMI.hpp"
#include "CglGomory.hpp"
#include "CglKnapsackCover.hpp"
#include "CglLandP.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglOddWheel.hpp"
#include "CglProbing.hpp"
#include "CglRedSplit.hpp"
#include "CglRedSplit2.hpp"
#include "CglResidualCapacity.hpp"
#include "CglStored.hpp"
#include "CglTwomir.hpp"
#include "CglZeroHalf.hpp"

// Register all cut generators on babModel based on parameter settings,
// then apply per-generator tuning (switches, accuracy, timing, cutDepth).
void installCutGenerators(
  CbcModel &babModel,
  CbcParameters &parameters,
  int complicatedInteger,
  bool dominatedCuts,
  bool miplib,
  const std::string &cgraphMode,
  int oldCliqueMode,
  int maxCallsBK,
  int bkClqExtMethod,
  CoinBronKerbosch::PivotingStrategy bkPivotingStrategy,
  int oddWExtMethod,
  int mixedRoundStrategy)
{
  int switches[30] = {};
  int accuracyFlag[30] = {};
  char doAtEnd[30] = {};
  int lagrangeanFlag = (parameters[CbcParam::MOREMOREMIPOPTIONS]->intVal() & (7 * 33554432)) >> 9;
#define ALL_LAGRANGEAN 1
  int numberGenerators = 0;
  std::map<int, int> translate;
  translate[CbcParameters::CGOff] = -100;
  translate[CbcParameters::CGOn] = -1;
  translate[CbcParameters::CGRoot] = -99;
  translate[CbcParameters::CGIfMove] = -98;
  translate[CbcParameters::CGForceOn] = 1;
  translate[CbcParameters::CGOnGlobal] = -1098;
  translate[CbcParameters::CGForceOnGlobal] = -999;
  translate[CbcParameters::CGForceOnBut] = 1;
  translate[CbcParameters::CGForceOnStrong] = 1;
  translate[CbcParameters::CGForceOnButStrong] = 1;
  translate[CbcParameters::CGStrongRoot] = -1;
  std::map<int, int> laTranslate;
  laTranslate[CbcParameters::CGEndOnlyRoot] = 1;
  laTranslate[CbcParameters::CGEndCleanRoot] = 2;
  laTranslate[CbcParameters::CGEndBothRoot] = 3;
  laTranslate[CbcParameters::CGEndOnly] = 4;
  laTranslate[CbcParameters::CGEndClean] = 5;
  laTranslate[CbcParameters::CGEndBoth] = 6;
  laTranslate[CbcParameters::CGOnlyAsWell] = 7;
  laTranslate[CbcParameters::CGOnlyAsWellRoot] = 13;
  laTranslate[CbcParameters::CGCleanAsWell] = 8;
  laTranslate[CbcParameters::CGCleanAsWellRoot] = 14;
  laTranslate[CbcParameters::CGBothAsWell] = 9;
  laTranslate[CbcParameters::CGBothAsWellRoot] = 15;
  laTranslate[CbcParameters::CGOnlyInstead] = 10;
  laTranslate[CbcParameters::CGCleanInstead] = 11;
  laTranslate[CbcParameters::CGBothInstead] = 12;
  int maximumSlowPasses = parameters[CbcParam::MAXSLOWCUTS]->intVal();

  // --- Probing ---
  int probingMode = parameters[CbcParam::PROBINGCUTS]->modeVal();
  if (probingMode) {
    CglProbing probingGen;
    probingGen.setUsingObjective(1);
    probingGen.setMaxPass(1);
    probingGen.setMaxPassRoot(1);
    probingGen.setMaxProbe(10);
    probingGen.setMaxProbeRoot(50);
    probingGen.setMaxLook(10);
    probingGen.setMaxLookRoot(50);
    probingGen.setMaxLookRoot(10);
    probingGen.setMaxElements(200);
    probingGen.setMaxElementsRoot(300);
    probingGen.setRowCuts(3);
    int numberColumns = babModel.solver()->getNumCols();
    if (probingMode > CbcParameters::CGForceOnBut) {
      probingGen.setMaxElements(numberColumns);
      probingGen.setMaxElementsRoot(numberColumns);
    }
    probingGen.setMaxProbeRoot(std::min(2000, numberColumns));
    probingGen.setMaxProbeRoot(123);
    probingGen.setMaxProbe(123);
    probingGen.setMaxLookRoot(20);
    if (probingMode == CbcParameters::CGForceOnBut || probingMode == CbcParameters::CGForceOnButStrong)
      probingGen.setRowCuts(-3);
    if (probingMode == CbcParameters::CGForceOnStrong || probingMode == CbcParameters::CGForceOnButStrong) {
      probingGen.setMaxProbeRoot(numberColumns);
      probingGen.setMaxProbe(numberColumns);
      probingGen.setMaxLook(50);
      probingGen.setMaxLookRoot(50);
    }
    if (probingMode == CbcParameters::CGStrongRoot) {
      probingGen.setMaxPassRoot(2);
      probingGen.setMaxProbeRoot(numberColumns);
      probingGen.setMaxLookRoot(numberColumns);
    }
    int iMode = translate[probingMode];
    if (probingMode == CbcParameters::CGOnGlobal)
      iMode = 1;
    babModel.addCutGenerator(&probingGen, iMode, "Probing");
    accuracyFlag[numberGenerators] = 5;
    switches[numberGenerators++] = 0;
  }

  // --- Gomory ---
  int gomoryMode = parameters[CbcParam::GOMORYCUTS]->modeVal();
  if (gomoryMode && (complicatedInteger != 1 || (gomoryMode == 1 || gomoryMode >= 4))) {
    CglGomory gomoryGen;
    gomoryGen.setLimitAtRoot(1000);
    gomoryGen.setLimit(50);
    // MORE_CUTS defaults (applied in the "Default strategy stuff" block)
#define MORE_CUTS
#ifdef MORE_CUTS
    gomoryGen.setAwayAtRoot(0.005);
#else
    gomoryGen.setAwayAtRoot(0.01);
#endif
    // Strategy 0 overrides awayAtRoot
    if (parameters[CbcParam::STRATEGY]->modeVal() == 0)
      gomoryGen.setAwayAtRoot(0.05);
    int numberColumns = babModel.getNumCols();
    if (gomoryMode == CbcParameters::CGForceOnBut) {
      gomoryMode = CbcParameters::CGForceOn;
      gomoryGen.setLimitAtRoot(numberColumns);
      gomoryGen.setLimit(numberColumns);
    } else if (gomoryMode == CbcParameters::CGForceOnStrong) {
      gomoryMode = CbcParameters::CGIfMove;
      gomoryGen.setLimitAtRoot(numberColumns);
      gomoryGen.setLimit(200);
    } else if (gomoryMode == CbcParameters::CGForceOnButStrong) {
      gomoryMode = CbcParameters::CGIfMove;
      gomoryGen.setLimitAtRoot(500);
      gomoryGen.setLimit(200);
    } else if (numberColumns > 5000) {
#ifdef MORE_CUTS2
      gomoryGen.setLimitAtRoot(numberColumns);
      gomoryGen.setLimit(200);
#else
      gomoryGen.setLimitAtRoot(2000);
#endif
    } else {
#ifdef MORE_CUTS2
      gomoryGen.setLimitAtRoot(numberColumns);
      gomoryGen.setLimit(200);
#endif
    }
    int cutLength = parameters[CbcParam::CUTLENGTH]->intVal();
    if (cutLength != -1) {
      gomoryGen.setLimitAtRoot(cutLength);
      if (cutLength < 10000000) {
        gomoryGen.setLimit(cutLength);
      } else {
        gomoryGen.setLimit(cutLength % 10000000);
      }
    }
    int laGomory = parameters[CbcParam::LAGOMORYCUTS]->modeVal();
    int gType = translate[gomoryMode];
    if (!laGomory) {
      babModel.addCutGenerator(&gomoryGen, translate[gomoryMode], "Gomory");
      accuracyFlag[numberGenerators] = 3;
      switches[numberGenerators++] = lagrangeanFlag;
    } else {
      laGomory = laTranslate[laGomory] - 1;
      int type = (laGomory % 3) + 1;
      int when = laGomory / 3;
      char atEnd = (when < 2) ? 1 : 0;
      int gomoryTypeMajor = 10;
      if (when != 3) {
        babModel.addCutGenerator(&gomoryGen, gType, "Gomory");
        accuracyFlag[numberGenerators] = 3;
        switches[numberGenerators++] = 0;
        if (when == 2) {
          gomoryTypeMajor = 20;
        } else if (when == 4) {
          gomoryTypeMajor = 20;
          when = 0;
        }
      } else {
        when--;
        gomoryTypeMajor = 20;
      }
      if (!when)
        gType = -99;
      gomoryGen.passInOriginalSolver(babModel.solver());
      if ((type & 1) != 0) {
        gomoryGen.setGomoryType(gomoryTypeMajor + 1);
        babModel.addCutGenerator(&gomoryGen, gType, "GomoryL1");
        accuracyFlag[numberGenerators] = 3;
        doAtEnd[numberGenerators] = atEnd;
        if (atEnd) {
          babModel.cutGenerator(numberGenerators)
            ->setMaximumTries(99999999);
          babModel.cutGenerator(numberGenerators)->setHowOften(1);
        }
        switches[numberGenerators++] = 16384;
      }
      if ((type & 2) != 0) {
        gomoryGen.setGomoryType(gomoryTypeMajor + 2);
        babModel.addCutGenerator(&gomoryGen, gType, "GomoryL2");
        accuracyFlag[numberGenerators] = 3;
        doAtEnd[numberGenerators] = atEnd;
        if (atEnd) {
          babModel.cutGenerator(numberGenerators)
            ->setMaximumTries(99999999);
          babModel.cutGenerator(numberGenerators)->setHowOften(1);
        }
        switches[numberGenerators++] = 32768;
      }
    }
  }

#ifdef CLIQUE_ANALYSIS
  if (miplib) {
    CglStored storedAmpl;
    if (!storedAmpl.sizeRowCuts()) {
      printf("looking at probing\n");
      babModel.addCutGenerator(&storedAmpl, 1, "Stored");
      accuracyFlag[numberGenerators] = 0;
      switches[numberGenerators++] = 0;
    }
  }
#endif

  // --- Knapsack ---
  int knapsackMode = parameters[CbcParam::KNAPSACKCUTS]->modeVal();
  if (knapsackMode) {
    CglKnapsackCover knapsackGen;
    babModel.addCutGenerator(&knapsackGen, translate[knapsackMode], "Knapsack");
    accuracyFlag[numberGenerators] = 1;
    switches[numberGenerators++] = -2;
  }

  // --- RedSplit ---
  int redsplitMode = parameters[CbcParam::REDSPLITCUTS]->modeVal();
  if (redsplitMode && !complicatedInteger) {
    CglRedSplit redsplitGen;
    babModel.addCutGenerator(&redsplitGen, translate[redsplitMode], "Reduce-and-split");
    accuracyFlag[numberGenerators] = 5;
    if (redsplitMode != CbcParameters::CGOn) {
      babModel.cutGenerator(numberGenerators)
        ->setMaximumTries(maximumSlowPasses);
      babModel.cutGenerator(numberGenerators)->setHowOften(10);
    }
    switches[numberGenerators++] = 1 | (ALL_LAGRANGEAN * lagrangeanFlag);
  }

  // --- RedSplit2 ---
  int redsplit2Mode = parameters[CbcParam::REDSPLIT2CUTS]->modeVal();
  if (redsplit2Mode && !complicatedInteger) {
    CglRedSplit2 redsplit2Gen;
    int maxLength = 256;
    if (redsplit2Mode > CbcParameters::CGRoot) {
      redsplit2Mode -= 2;
      maxLength = COIN_INT_MAX;
    }
    CglRedSplit2Param &rs2params = redsplit2Gen.getParam();
    rs2params.setMaxNonzeroesTab(maxLength);
    babModel.addCutGenerator(&redsplit2Gen, translate[redsplit2Mode], "Reduce-and-split(2)");
    accuracyFlag[numberGenerators] = 5;
    if (redsplit2Mode != CbcParameters::CGOn) {
      babModel.cutGenerator(numberGenerators)
        ->setHowOften(maximumSlowPasses);
      babModel.cutGenerator(numberGenerators)
        ->setMaximumTries(maximumSlowPasses);
      babModel.cutGenerator(numberGenerators)->setHowOften(5);
    }
    switches[numberGenerators++] = 1 | (ALL_LAGRANGEAN * lagrangeanFlag);
  }

  // --- GMI ---
  int GMIMode = parameters[CbcParam::GMICUTS]->modeVal();
  if (GMIMode && !complicatedInteger) {
    CglGMI GMIGen;
    if (GMIMode > CbcParameters::CGOnGlobal) {
      GMIMode -= 5;
      CglGMIParam &gmiParams = GMIGen.getParam();
      gmiParams.setMaxSupportRel(1.0);
    }
    babModel.addCutGenerator(&GMIGen, translate[GMIMode], "Gomory(2)");
    if (GMIMode == CbcParameters::CGOnGlobal) {
      GMIMode = CbcParameters::CGRoot;
      doAtEnd[numberGenerators] = 1;
      babModel.cutGenerator(numberGenerators)
        ->setMaximumTries(99999999);
      babModel.cutGenerator(numberGenerators)->setHowOften(1);
    }
    accuracyFlag[numberGenerators] = 5;
    switches[numberGenerators++] = 0 | (ALL_LAGRANGEAN * lagrangeanFlag);
  }

  // --- Clique ---
  int cliqueMode = parameters[CbcParam::CLIQUECUTS]->modeVal();
  int oddWheelMode = parameters[CbcParam::ODDWHEELCUTS]->modeVal();
  if (cliqueMode) {
    CglBKClique bkCliqueGen;
    bkCliqueGen.setMaxCallsBK(maxCallsBK);
    bkCliqueGen.setExtendingMethod(bkClqExtMethod);
    bkCliqueGen.setPivotingStrategy(bkPivotingStrategy);
    babModel.addCutGenerator(&bkCliqueGen, translate[cliqueMode], "Clique");
    accuracyFlag[numberGenerators] = 0;
    switches[numberGenerators++] = 0;
  } else if (cgraphMode == "clq") {
    CglClique clique;
    clique.setStarCliqueReport(false);
    clique.setRowCliqueReport(false);
    clique.setMinViolation(0.05);
    oddWheelMode = 0;
    parameters[CbcParam::ODDWHEELCUTS]->setModeVal(CbcParameters::CGOff);
    babModel.addCutGenerator(&clique, translate[oldCliqueMode], "Clique");
    accuracyFlag[numberGenerators] = 0;
    switches[numberGenerators++] = 0;
  }

  // --- OddWheel ---
  if (oddWheelMode) {
    CglOddWheel oddWheelGen;
    oddWheelGen.setExtendingMethod(oddWExtMethod);
    babModel.addCutGenerator(&oddWheelGen, translate[oddWheelMode], "OddWheel");
    accuracyFlag[numberGenerators] = 0;
    switches[numberGenerators++] = 0;
  }

  // --- MIR ---
  int mixedMode = parameters[CbcParam::MIRCUTS]->modeVal();
  if (mixedMode) {
    CglMixedIntegerRounding2 mixedGen(1, true, 1);
    mixedGen.setDoPreproc(1);
    if (mixedRoundStrategy != 1)
      mixedGen.setMAXAGGR_(mixedRoundStrategy);
    babModel.addCutGenerator(&mixedGen, translate[mixedMode], "MixedIntegerRounding2");
    accuracyFlag[numberGenerators] = 2;
    switches[numberGenerators++] = 0 | (ALL_LAGRANGEAN * lagrangeanFlag);
  }

  // --- FlowCover ---
  int flowMode = parameters[CbcParam::FLOWCUTS]->modeVal();
  if (flowMode) {
    CglFlowCover flowGen;
    babModel.addCutGenerator(&flowGen, translate[flowMode], "FlowCover");
    accuracyFlag[numberGenerators] = 2;
    switches[numberGenerators++] = 0 | (ALL_LAGRANGEAN * lagrangeanFlag);
  }

  // --- TwoMir ---
  int twomirMode = parameters[CbcParam::TWOMIRCUTS]->modeVal();
  if (twomirMode && (complicatedInteger != 1 || (twomirMode == CbcParameters::CGOn || twomirMode >= CbcParameters::CGForceOn))) {
    CglTwomir twomirGen;
    twomirGen.setMaxElements(250);
    // MORE_CUTS defaults
#ifdef MORE_CUTS
    twomirGen.setAwayAtRoot(0.005);
    twomirGen.setAway(0.01);
#else
    twomirGen.setAwayAtRoot(0.01);
    twomirGen.setAway(0.01);
#endif
    int numberColumns = babModel.getNumCols();
    if (twomirMode == CbcParameters::CGForceOnBut) {
      twomirMode = CbcParameters::CGForceOn;
      twomirGen.setMaxElements(numberColumns);
    } else if (numberColumns > 5000 && twomirMode == CbcParameters::CGForceOn) {
      twomirGen.setMaxElements(2000);
    }
    int laTwomir = parameters[CbcParam::LATWOMIRCUTS]->modeVal();
    int twomirType = translate[twomirMode];
    if (!laTwomir) {
      babModel.addCutGenerator(&twomirGen, translate[twomirMode], "TwoMirCuts");
      accuracyFlag[numberGenerators] = 4;
      switches[numberGenerators++] = 1 | lagrangeanFlag;
    } else {
      laTwomir = laTranslate[laTwomir] - 1;
      int type = (laTwomir % 3) + 1;
      int when = laTwomir / 3;
      char atEnd = (when < 2) ? 1 : 0;
      int twomirTypeMajor = 10;
      if (when < 3) {
        babModel.addCutGenerator(&twomirGen, translate[twomirMode], "TwoMirCuts");
        accuracyFlag[numberGenerators] = 4;
        switches[numberGenerators++] = 1;
        if (when == 2)
          twomirTypeMajor = 10;
      } else {
        when--;
        twomirTypeMajor = 20;
      }
      if (!when)
        twomirType = -99;
      twomirGen.passInOriginalSolver(babModel.solver());
      if ((type & 1) != 0) {
        twomirGen.setTwomirType(twomirTypeMajor + 1);
        babModel.addCutGenerator(&twomirGen, twomirType, "TwoMirCutsL1");
        accuracyFlag[numberGenerators] = 4;
        doAtEnd[numberGenerators] = atEnd;
        switches[numberGenerators++] = (atEnd ? 0 : 1) | 16384;
      }
      if ((type & 2) != 0) {
        twomirGen.setTwomirType(twomirTypeMajor + 2);
        babModel.addCutGenerator(&twomirGen, twomirType, "TwoMirCutsL2");
        accuracyFlag[numberGenerators] = 4;
        doAtEnd[numberGenerators] = atEnd;
        switches[numberGenerators++] = (atEnd ? 0 : 1) | 32768;
      }
    }
  }

  // --- LandP ---
#ifndef DEBUG_MALLOC
  int landpMode = parameters[CbcParam::LANDPCUTS]->modeVal();
  if (landpMode) {
    CglLandP landpGen;
    landpGen.parameter().maximumCutLength = 2000;
    landpGen.validator().setMinViolation(1.0e-4);
    if (landpMode == CbcParameters::CGOnGlobal) {
      landpGen.parameter().maximumCutLength = 2000000;
      landpMode = CbcParameters::CGIfMove;
    }
    babModel.addCutGenerator(&landpGen, translate[landpMode], "LiftAndProject");
    accuracyFlag[numberGenerators] = 5;
    if (landpMode != CbcParameters::CGOn) {
      babModel.cutGenerator(numberGenerators)
        ->setMaximumTries(maximumSlowPasses);
      babModel.cutGenerator(numberGenerators)->setHowOften(10);
    }
    switches[numberGenerators++] = 1 | (ALL_LAGRANGEAN * lagrangeanFlag);
  }
#endif

  // --- ResidualCapacity ---
  int residualCapacityMode = parameters[CbcParam::RESIDCAPCUTS]->modeVal();
  if (residualCapacityMode) {
    CglResidualCapacity residualCapacityGen;
    residualCapacityGen.setDoPreproc(1);
    babModel.addCutGenerator(&residualCapacityGen,
      translate[residualCapacityMode], "ResidualCapacity");
    accuracyFlag[numberGenerators] = 5;
    switches[numberGenerators++] = 1 | (ALL_LAGRANGEAN * lagrangeanFlag);
  }

  // --- ZeroHalf ---
  int zerohalfMode = parameters[CbcParam::ZEROHALFCUTS]->modeVal();
  if (zerohalfMode) {
    CglZeroHalf zerohalfGen;
    zerohalfGen.setSepGraphSparseThreshold(parameters[CbcParam::ZEROHALFSPARSETHRESH]->intVal());
    zerohalfGen.setRowMaxPairCount(parameters[CbcParam::ZEROHALFROWMAXPAIRCOUNT]->intVal());
    zerohalfGen.setRowMaxFractionalCount(parameters[CbcParam::ZEROHALFROWMAXFRACTIONALCOUNT]->intVal());
    if (zerohalfMode > CbcParameters::CGForceOn)
      zerohalfGen.setFlags(1);
    babModel.addCutGenerator(&zerohalfGen, translate[zerohalfMode], "ZeroHalf");
    accuracyFlag[numberGenerators] = 5;
    CglZeroHalf *storedZeroHalf =
      dynamic_cast<CglZeroHalf *>(babModel.cutGenerator(numberGenerators)->generator());
    if (storedZeroHalf)
      storedZeroHalf->setSepGraphSparseThreshold(
        parameters[CbcParam::ZEROHALFSPARSETHRESH]->intVal());
    if (storedZeroHalf)
      storedZeroHalf->setRowMaxPairCount(
        parameters[CbcParam::ZEROHALFROWMAXPAIRCOUNT]->intVal());
    if (storedZeroHalf)
      storedZeroHalf->setRowMaxFractionalCount(
        parameters[CbcParam::ZEROHALFROWMAXFRACTIONALCOUNT]->intVal());
    babModel.cutGenerator(numberGenerators)->setNeedsRefresh(true);
    switches[numberGenerators++] = 2;
  }

  if (dominatedCuts)
    babModel.setSpecialOptions(babModel.specialOptions() | 64);

  // Per-generator tuning
  numberGenerators = babModel.numberCutGenerators();
  int cutDepth = parameters[CbcParam::CUTDEPTH]->intVal();
  for (int iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CbcCutGenerator *generator = babModel.cutGenerator(iGenerator);
    int howOften = generator->howOften();
    int iSwitch = switches[iGenerator];
    int iSwitch2, iSwitch3;
    if (iSwitch >= 0) {
      iSwitch2 = iSwitch & 127;
      iSwitch3 = iSwitch & ~16383;
      generator->setSwitches(generator->switches() | iSwitch3);
    } else {
      iSwitch2 = iSwitch;
    }
    if (howOften == -98 || howOften == -99 || generator->maximumTries() > 0)
      generator->setSwitchOffIfLessThan(iSwitch2);
    generator->setInaccuracy(accuracyFlag[iGenerator]);
    if (doAtEnd[iGenerator]) {
      generator->setWhetherCallAtEnd(true);
    }
    generator->setTiming(true);
    if (cutDepth >= 0)
      generator->setWhatDepth(cutDepth);
  }
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
