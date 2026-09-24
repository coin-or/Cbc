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
#include "CglImpliedClique.hpp"
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

#include <cstdlib>

namespace {
int envInt(const char *name, int def)
{
  const char *v = getenv(name);
  return v ? atoi(v) : def;
}
double envDouble(const char *name, double def)
{
  const char *v = getenv(name);
  return v ? atof(v) : def;
}
} // namespace

// Register all cut generators on babModel based on parameter settings,
// then apply per-generator tuning (switches, accuracy, timing, cutDepth).
void installCutGenerators(
  CbcModel &babModel,
  CbcParameters &parameters,
  int complicatedInteger,
  bool dominatedCuts,
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
  std::map< int, int > translate;
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
  std::map< int, int > laTranslate;
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
  // See CbcSolver.cpp's REDSPLIT2CUTS/GMICUTS/LANDPCUTS default ("root"): a
  // 2026-09 sanity-suite sweep found these only pay for themselves in
  // aggregate within a "sweet spot" instance-size window. Below it, the
  // adaptive root cut-generator skip's own size floor
  // (ADAPTIVE_SKIP_MIN_COLS, default 500 -- see CbcModel::serialCuts())
  // never engages, so the generator just burns time every pass with no
  // throttling; above it, even a single call can be expensive enough to
  // blow the node/time budget before enough misses accumulate for backoff
  // to kick in. This only narrows the "root" default itself -- a user who
  // wants a generator regardless of size can still force it on explicitly
  // (e.g. -redsplit2Cuts=on/ifmove).
  const int numberColumnsForGate = babModel.getNumCols();
  // Window bounds are env-var overridable (no rebuild) purely to speed up
  // experimentation/tuning sweeps; the defaults below (500/50000) are the
  // ones actually shipped.
  int gateMinCols = 500;
  int gateMaxCols = 50000;
  if (const char *s = getenv("CBC_CUT_ROOT_GATE_MIN_COLS"))
    gateMinCols = atoi(s);
  if (const char *s = getenv("CBC_CUT_ROOT_GATE_MAX_COLS"))
    gateMaxCols = atoi(s);
  // Row-count companion to the column gate above. Found 2026-09 while
  // investigating an OOM crash: CglRedSplit2::generateCuts() allocates a
  // "bufflambda" work array sized maxNumComputedCuts*nrow ints whenever
  // maxNumCuts<maxNumComputedCuts (see CglRedSplit2.cpp ~L2011,2029) --
  // nrow there is solver->getNumRows(), i.e. the *actual* preprocessed row
  // count, completely unrelated to the column-count gate above. A handful
  // of real 2017+spp instances have millions of rows despite a modest
  // column count (e.g. neos-3402454-bohle: 2496 cols, comfortably inside
  // the [500,50000) column window, but 2,881,228 rows) -- for those, even
  // this file's own reduced maxNumComputedCuts default (2000, see below)
  // allocates ~23GB in that one array alone, enough to OOM-kill the whole
  // process (observed directly: a real crash at total-vm ~73GB,
  // anon-rss ~68GB, see kernel oom-killer log). GMI/LandP's own
  // per-cut work arrays are only O(nrow) (no maxNumComputedCuts-like
  // multiplier), so they don't share this specific failure mode and are
  // deliberately not gated on rows here -- only RedSplit2 needs it.
  int gateMaxRowsRedsplit2 = 200000;
  if (const char *s = getenv("CBC_CUT_ROOT_GATE_MAX_ROWS_REDSPLIT2"))
    gateMaxRowsRedsplit2 = atoi(s);
  auto gateRootDefault = [numberColumnsForGate, gateMinCols, gateMaxCols](int mode) {
    if (mode == CbcParameters::CGRoot
      && (numberColumnsForGate < gateMinCols || numberColumnsForGate >= gateMaxCols))
      return static_cast< int >(CbcParameters::CGOff);
    return mode;
  };
  const int numberRowsForGate = babModel.getNumRows();
  auto gateRootDefaultRedsplit2 = [numberRowsForGate, gateMaxRowsRedsplit2](int mode) {
    if (numberRowsForGate >= gateMaxRowsRedsplit2)
      return static_cast< int >(CbcParameters::CGOff);
    return mode;
  };
  // Per-round/per-lifetime cost caps for RedSplit2/GMI/LandP -- see the
  // comments at each generator's setup below for the exact semantics
  // (RedSplit2's timeLimit is per-call, LandP's is a cumulative
  // whole-solve budget). Env-var overridable, no rebuild, purely to let a
  // sweep A/B these against their un-throttled CglXxxParam defaults
  // without needing two binaries -- same pattern as
  // CBC_CUT_ROOT_GATE_MIN_COLS/MAX_COLS above and CBC_CUTPOOL_FILTER_*.
  // Defaults match the "old-unthrottled" config from the cutgen-throttle-sweep
  // experiment (348 root-fixture instances, 128 nodes, 2h cap): it edged out
  // the tighter-throttle variants on both dual gap (30.66% vs 30.76%/31.20%)
  // and primal gap (12.85% vs 13.11%/13.58%), with no meaningful time cost
  // at this node budget -- see session notes for 2026-09-24. The independent
  // RedSplit2 memory-safety fixes (row-gate, dynamic buffer cap on
  // maxNumComputedCuts*nrow, and the mTab*card_contNonBasicVar tableau-size
  // guard in CglRedSplit2.cpp) remain in effect regardless of these values.
  const double redsplit2TimeLimit = envDouble("CBC_REDSPLIT2_TIME_LIMIT", 60.0);
  const int redsplit2MaxNumCuts = envInt("CBC_REDSPLIT2_MAX_NUM_CUTS", 10000);
  const int redsplit2MaxNumComputedCuts = envInt("CBC_REDSPLIT2_MAX_NUM_COMPUTED_CUTS", 10000);
  // Defense-in-depth for the bufflambda blowup described above: even if a
  // user forces RedSplit2 on explicitly (-redsplit2Cuts=on/ifmove, which
  // bypasses gateRootDefaultRedsplit2 -- that gate only narrows the
  // "root" default) or the row count sits just under
  // gateMaxRowsRedsplit2, cap maxNumComputedCuts so the
  // maxNumComputedCuts*nrow allocation can never exceed a fixed memory
  // budget, shrinking gracefully as nrow grows instead of either being
  // unbounded or an all-or-nothing gate. 50,000,000 ints (~200MB) is
  // generous next to Gomory's near-instant footprint while making even a
  // million-plus-row instance's worst case a bounded, known quantity.
  const long redsplit2MaxBufferInts = envInt("CBC_REDSPLIT2_MAX_BUFFER_INTS", 50000000);
  const int gmiHowOften = envInt("CBC_GMI_HOW_OFTEN", 1);
  const double landpTimeLimit = envDouble("CBC_LANDP_TIME_LIMIT", 1e30);
  const double landpSingleCutTimeLimit = envDouble("CBC_LANDP_SINGLE_CUT_TIME_LIMIT", 1e30);
  const int landpMaxCutPerRound = envInt("CBC_LANDP_MAX_CUT_PER_ROUND", 5000);

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
    probingGen.setRowCuts(-3);
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
  int redsplit2Mode = gateRootDefaultRedsplit2(
    gateRootDefault(parameters[CbcParam::REDSPLIT2CUTS]->modeVal()));
  if (redsplit2Mode && !complicatedInteger) {
    CglRedSplit2 redsplit2Gen;
    int maxLength = 256;
    if (redsplit2Mode > CbcParameters::CGRoot) {
      redsplit2Mode -= 2;
      maxLength = COIN_INT_MAX;
    }
    CglRedSplit2Param &rs2params = redsplit2Gen.getParam();
    rs2params.setMaxNonzeroesTab(maxLength);
    // CglRedSplit2Param's own defaults (timeLimit=60s, maxNumCuts=
    // maxNumComputedCuts=10000) are unbounded relative to what Gomory
    // allows per round (setLimitAtRoot(1000)/setLimit(50)) -- a single
    // slow round could burn a large fraction of a node/time budget before
    // boundStallAware()'s backoff ever gets a chance to react. Bring this
    // generator's per-round cost onto the same order of magnitude as
    // Gomory's, still generous enough for its typically much smaller cut
    // yield per pass.
    //
    // On top of that: CglRedSplit2::generateCuts() allocates a
    // "bufflambda" work array sized maxNumComputedCuts*nrow ints whenever
    // maxNumCuts<maxNumComputedCuts (see CglRedSplit2.cpp), where nrow is
    // this model's actual row count -- unrelated to maxNumComputedCuts
    // itself. gateRootDefaultRedsplit2 above already turns this generator
    // off outright for its "root" default once rows get extreme, but that
    // only applies to the default mode; also shrink maxNumComputedCuts/
    // maxNumCuts here so the product with this model's real row count
    // never exceeds redsplit2MaxBufferInts, regardless of how the
    // generator was enabled (including -redsplit2Cuts=on/ifmove, or
    // "root" instances just under the row gate's threshold).
    const int numberRows = babModel.getNumRows();
    int cappedMaxNumComputedCuts = redsplit2MaxNumComputedCuts;
    if (numberRows > 0) {
      const long byMemory = redsplit2MaxBufferInts / numberRows;
      cappedMaxNumComputedCuts = static_cast< int >(
        std::max< long >(1, std::min< long >(redsplit2MaxNumComputedCuts, byMemory)));
    }
    const int cappedMaxNumCuts = std::min(redsplit2MaxNumCuts, cappedMaxNumComputedCuts);
    rs2params.setTimeLimit(redsplit2TimeLimit);
    rs2params.setMaxNumComputedCuts(cappedMaxNumComputedCuts);
    rs2params.setMaxNumCuts(cappedMaxNumCuts);
    babModel.addCutGenerator(&redsplit2Gen, translate[redsplit2Mode], "Reduce-and-split(2)");
    accuracyFlag[numberGenerators] = 5;
    if (redsplit2Mode != CbcParameters::CGOn) {
      babModel.cutGenerator(numberGenerators)
        ->setHowOften(maximumSlowPasses);
      babModel.cutGenerator(numberGenerators)
        ->setMaximumTries(maximumSlowPasses);
      babModel.cutGenerator(numberGenerators)->setHowOften(5);
    }
    // Expensive generator: only counts as "productive" for the adaptive
    // root skip (CbcModel::serialCuts()) when it actually moves the bound,
    // not just when it emits some (possibly near-useless) cuts.
    babModel.cutGenerator(numberGenerators)->setBoundStallAware(true);
    switches[numberGenerators++] = 1 | (ALL_LAGRANGEAN * lagrangeanFlag);
  }

  // --- GMI ---
  int GMIMode = gateRootDefault(parameters[CbcParam::GMICUTS]->modeVal());
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
    // Unlike RedSplit2/LandP just above, GMI's normal ("root") path was
    // never given the same setMaximumTries()/setHowOften() throttle --
    // an oversight, not a deliberate exemption, since GMI is no cheaper
    // per call. Bring it in line: skip most passes past the root the same
    // way RedSplit2/LandP already do.
    if (GMIMode != CbcParameters::CGOn && GMIMode != CbcParameters::CGOnGlobal) {
      babModel.cutGenerator(numberGenerators)
        ->setMaximumTries(maximumSlowPasses);
      babModel.cutGenerator(numberGenerators)->setHowOften(gmiHowOften);
    }
    // See RedSplit2 comment above: only count as productive when the bound
    // actually moves, not merely when some cut is produced.
    babModel.cutGenerator(numberGenerators)->setBoundStallAware(true);
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
    // Also expensive at root on dense/large conflict graphs; only count as
    // productive for the adaptive root skip when it actually moves the
    // bound (see RedSplit2/GMI/LandP comment above).
    babModel.cutGenerator(numberGenerators)->setBoundStallAware(true);
    switches[numberGenerators++] = 0;
  }

  // --- ImpliedClique ---
  int impliedCliqueMode = parameters[CbcParam::IMPLIEDCLIQUECUTS]->modeVal();
  if (impliedCliqueMode) {
    CglImpliedClique impliedCliqueGen;
    babModel.addCutGenerator(&impliedCliqueGen, translate[impliedCliqueMode], "ImpliedClique");
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
  int landpMode = gateRootDefault(parameters[CbcParam::LANDPCUTS]->modeVal());
  if (landpMode) {
    CglLandP landpGen;
    // Lowered from 2000: a 2026-09 fixture-replay sweep found capping cut
    // length here cuts LandP's worst-case root time by ~50x for only
    // ~1-5% of its bound value (see BENCHMARKING-CUT-GENERATORS.md).
    landpGen.parameter().maximumCutLength = 200;
    // CglLandP::Parameters' own defaults leave timeLimit/singleCutTimeLimit
    // at COIN_DBL_MAX (i.e. genuinely unbounded) and maxCutPerRound at
    // 5000. Unlike CglRedSplit2's timeLimit (reset every generateCuts()
    // call -- a true per-round budget), CglLandP::timeLimit_ is a
    // *cumulative* CPU-time allowance consumed across the generator's
    // entire lifetime in this solve (CglLandP.cpp adds/subtracts
    // CoinCpuTime() around each call, decrementing a persistent running
    // total; once it goes negative, pivotLimit is forced to 0 and the
    // generator effectively self-disables for the rest of the solve --
    // see CglLandP.cpp ~L1047/1221/904). So this value should be sized as
    // a total time budget for the whole B&B, not a per-call one; 30s is
    // generous next to Gomory's near-instant per-round cost while still
    // capping the worst case where LandP's pivot search fails to converge
    // repeatedly. singleCutTimeLimit additionally caps any single cut
    // attempt's own pivot search (used via
    // std::min(timeLimit, singleCutTimeLimit) in CglLandPSimplex.cpp), so
    // one degenerate candidate can't consume the whole remaining budget.
    landpGen.parameter().timeLimit = landpTimeLimit;
    landpGen.parameter().singleCutTimeLimit = landpSingleCutTimeLimit;
    landpGen.parameter().maxCutPerRound = landpMaxCutPerRound;
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
    // See RedSplit2 comment above: only count as productive when the bound
    // actually moves, not merely when some cut is produced.
    babModel.cutGenerator(numberGenerators)->setBoundStallAware(true);
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
    CglZeroHalf *storedZeroHalf = dynamic_cast< CglZeroHalf * >(babModel.cutGenerator(numberGenerators)->generator());
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
    // Also expensive on dense graphs at root; only count as productive for
    // the adaptive root skip when it actually moves the bound (see
    // RedSplit2/GMI/LandP comment above).
    babModel.cutGenerator(numberGenerators)->setBoundStallAware(true);
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
