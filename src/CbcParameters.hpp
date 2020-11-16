/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/

#ifndef CbcSettings_H
#define CbcSettings_H

/* \file CbcSettings.hpp
   \brief Declarations for parameters of Cbc.
*/

#include "CoinMessageHandler.hpp"
#include "CoinParam.hpp"

#include "ClpParam.hpp"

#include "CglClique.hpp"
#include "CglCutGenerator.hpp"
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
#include "CglTwomir.hpp"
#include "CglZeroHalf.hpp"

#include "CbcModel.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicDINS.hpp"
#include "CbcHeuristicDW.hpp"
#include "CbcHeuristicDive.hpp"
#include "CbcHeuristicDiveCoefficient.hpp"
#include "CbcHeuristicDiveFractional.hpp"
#include "CbcHeuristicDiveGuided.hpp"
#include "CbcHeuristicDiveLineSearch.hpp"
#include "CbcHeuristicDivePseudoCost.hpp"
#include "CbcHeuristicDiveVectorLength.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcHeuristicPivotAndFix.hpp"
#include "CbcHeuristicRENS.hpp"
#include "CbcHeuristicRINS.hpp"
#include "CbcHeuristicRandRound.hpp"
#include "CbcHeuristicVND.hpp"
#include "CbcParam.hpp"
#include "CbcTreeLocal.hpp"

#include "CbcGenMessages.hpp"

enum CbcParamCode {
   CBC_FIRSTPARAM = 0,
   
   // Unused paramters that we may delete
   CBC_FIRSTUNUSEDPARAM,
   BRANCHSTRATEGY,
   CLEARCUTS,
   SOLVER,
   CBC_LASTUNUSEDPARAM,
   
   // Help and Information Parameters
   CBC_FIRSTHELPPARAM,
   GENERALQUERY,
   FULLGENERALQUERY,
   HELP,
   PRINTVERSION,
   CBC_LASTHELPPARAM,
   
   // Action Parameters
   CBC_FIRSTACTIONPARAM,
   BAB,
   DUMMY,
   ENVIRONMENT,
   EXIT,
   EXPORT,
   GMPL_SOLUTION,
   IMPORT,
   MIPLIB,
   OUTDUPROWS,
   RESTORE,
   SAVE,
   SHOWUNIMP,
   SOLVECONTINUOUS,
   STATISTICS,
   STDIN,
   STRENGTHEN,
   UNITTEST,
   CBC_LASTACTIONPARAM,
   
   // String (Directory) Parameters
   CBC_FIRSTSTRPARAM,
   CSVSTATISTICS,
   DEBUG,
   DIRECTORY,
   DIRSAMPLE,
   DIRNETLIB,
   DIRMIPLIB,
   MIPSTART,
   NEXTBESTSOLUTION,
   PRINTMASK,
   PRIORITYIN,
   SOLUTION,
   SAVESOL,
   CBC_LASTSTRPARAM,
   
   // Cut Parameters
   CBC_FIRSTCUTPARAM,
   CUTSTRATEGY,
   CLIQUECUTS,
   FLOWCUTS,
   GMICUTS,
   GOMORYCUTS,
   KNAPSACKCUTS,
   LAGOMORYCUTS,
   LANDPCUTS,
   LATWOMIRCUTS,
   MIRCUTS,
   ODDHOLECUTS, // Not used
   ODDWHEELCUTS,
   PROBINGCUTS,
   REDSPLITCUTS,
   REDSPLIT2CUTS,
   RESIDCAPCUTS,
   TWOMIRCUTS,
   ZEROHALFCUTS,
   CBC_LASTCUTPARAM,
   
   // Heuristic Parameters
   CBC_FIRSTHEURPARAM,
   COMBINE,
   CROSSOVER,
   DINS,
   DIVINGC,
   DIVINGF,
   DIVINGG,
   DIVINGL,
   DIVINGP,
   DIVINGS,
   DIVINGV,
   DW,
   FPUMP,
   GREEDY,
   HEURISTICSTRATEGY,
   LOCALTREE,
   NAIVE,
   PIVOTANDFIX,
   PIVOTANDCOMPLEMENT,
   PROXIMITY,
   RANDROUND,
   RENS,
   RINS,
   ROUNDING,
   VND,
   CBC_LASTHEURPARAM,
   
   // On/Off Parameters
   CBC_FIRSTBOOLPARAM,
   CPX,
   DOHEURISTIC,
   ERRORSALLOWED,
   EXTRAVARIABLES,
   MESSAGES,
   PREPROCNAMES,
   SOS,
   USESOLUTION,
   CBC_LASTBOOLPARAM,
   
   // Keyword Parameters
   CBC_FIRSTKWDPARAM,
   COMMANDPRINTLEVEL,
   CLQSTRENGTHENING,
   BRANCHPRIORITY,
   CUTOFFCONSTRAINT,
   INTPRINT,
   NODESTRATEGY,
   ORBITAL,
   PREPROCESS,
   SOSPRIORITIZE,
   STRATEGY,
   TIMEMODE,
   USECGRAPH,
   CBC_LASTKWDPARAM,
   
   // Integer Parameters
   CBC_FIRSTINTPARAM,
   BKPIVOTINGSTRATEGY,
   BKMAXCALLS,
   BKCLQEXTMETHOD,
   CPP,
   CUTDEPTH,
   CUTLENGTH,
   CUTPASSINTREE,
   DEPTHMINIBAB,
   DIVEOPT,
   DIVEOPTSOLVES,
   EXPERIMENT,
   EXTRA1,
   EXTRA2,
   EXTRA3,
   EXTRA4,
   FPUMPITS,
   FPUMPTUNE,
   FPUMPTUNE2,
   HEUROPTIONS,
   LOGLEVEL,
   LPLOGLEVEL,
   MAXHOTITS,
   MAXSAVEDSOLS,
   MAXSLOWCUTS,
   MOREMOREMIPOPTIONS,
   MULTIPLEROOTS,
   ODDWEXTMETHOD,
   OUTPUTFORMAT,
   PRINTOPTIONS,
   PROCESSTUNE,
   RANDOMSEED,
   STRONGSTRATEGY,
   TESTOSI,
   THREADS,
   USERCBC,
   VERBOSE,
   VUBTRY,
   CBC_LASTINTPARAM,
   
   // Double Parameters
   CBC_FIRSTDBLPARAM,
   ARTIFICIALCOST,
   DEXTRA3,
   DEXTRA4,
   DEXTRA5,
   DJFIX,
   FAKECUTOFF,
   FAKEINCREMENT,
   SMALLBAB,
   TIGHTENFACTOR,
   CBC_LASTDBLPARAM,
   
   CBC_FIRSTMODELPARAM,
   ALLOWABLEGAP,
   COSTSTRATEGY,
   //CUTDEPTH,
   CUTOFF,
   CUTPASS,
   DIRECTION,
   GAPRATIO,
   INCREMENT,
   INFEASIBILITYWEIGHT,
   INTEGERTOLERANCE,
   //LOGLEVEL,
   MAXIMIZE,
   MAXNODES,
   MAXNODESNOTIMPROVING,   // Added
   MAXSECONDSNOTIMPROVING, // Added
   MAXSOLS,                // Added
   MINIMIZE,
   MIPOPTIONS,
   MOREMIPOPTIONS,
   NUMBERANALYZE,
   NUMBERBEFORE,
   REVERSE,
   STRONGBRANCHING,
   TIMELIMIT,
   CBC_LASTMODELPARAM,

   CBC_LASTPARAM
};

/*
  It turns out that doxygen is not good with anonymous structures. Hence the
  `struct name_struct' style used for structured fields in CbcSettings.
*/

/*
 */

/* \brief Cbc algorithm settings class

  This class defines and stores the parameters used to control the operation 
  of Cbc.
*/

class CbcSettings {

public:
  /*! \name Constructors and destructors */
  //@{

  /*! \brief Default constructor */
  CbcSettings();

  /*! \brief Destructor */
  ~CbcSettings();

  //@}

  /*! \name Enumeration types used for CBC keyword parameters */
  //@{

  /*! \brief Codes to control integer preprocessing

    - IPPOff: Integer preprocessing is off.
    - IPPOn:  Integer preprocessing is on.
    - IPPSave: IPPOn, plus preprocessed system will be saved to presolved.mps.
    - IPPEqual: IPPOn, plus `<=' cliques are converted to `=' cliques.
    - IPPSOS: IPPOn, plus will create SOS sets (see below).
    - IPPTrySOS: IPPOn, plus will create SOS sets (see below).
    - IPPEqualAll: IPPOn, plus turns all valid inequalities into equalities
    with integer slacks.
    - IPPStrategy: look to CbcStrategy object for instructions.

    IPPSOS will create SOS sets if all binary variables (except perhaps one)
    can be covered by SOS sets with no overlap between sets. IPPTrySOS will
    allow any number of binary variables to be uncovered.
  */

  enum IPPMode {
    IPPOff = 0,
    IPPOn,
    IPPSave,
    IPPEqual,
    IPPSOS,
    IPPTrySOS,
    IPPEqualAll,
    IPPStrategy,
    IPPEndMarker
  };

  /*! \brief What parameters to print

    - displayAll:
    - displayLowHigh:
    - displayHigh:
    - displayEndMarker

   */

  enum CommandDisplayMode {
     displayAll = 0,
     displayLowHigh,
     displayHigh,
     displayEndMarker };

  /*! \brief When to do clique strengthening

    - ClqOff:
    - ClqAfter:
    - ClqBefore:
    - ClqEndMarker:

  */

  enum ClqStrMode { ClqStrOff = 0, ClqStrAfter, ClqStrBefore, ClqStrEndMarker };

  /*! \brief What node strategy to use

    - NSHybrid:
    - NSFewest:
    - NSDepth:
    - NSUpFewest:
    - NSUpDepth:
    - NSDownDepth

   */

  enum NodeStrategy {
    NSHybrid = 0,
    NSFewest,
    NSDepth,
    NSUpFewest,
    NSUpDepth,
    NSDownFewest,
    NSDownDepth
  };

  /*! \brief What orbital branching strategy to use

    - OBOff:
    - OBOn:
    - OBSlowish:
    - OBStrong:
    - OBForce:
    - OBSimple:
    - OBMorePrinting:
    - OBEndMarker
   */

  enum OrbitalStrategy {
    OBOff = 0,
    OBOn,
    OBSlowish,
    OBStrong,
    OBForce,
    OBSimple,
    OBMorePrinting,
    OBEndMarker
  };

  /*! \brief What SOS prioritization strategy to use

    - SOSOff:
    - SOSHigh:
    - SOSLow:
    - SOSOrderHigh:
    - SOSOrderLow:
    - SOSEndMarker
   */

  enum SOSStrategy {
    SOSOff = 0,
    SOSHigh,
    SOSLow,
    SOSOrderHigh,
    SOSOrderLow,
    SOSEndMarker
  };

  /*! \brief What overall strategy to use

    - StrategyDefaulty = 0,
    - StrategyEasy,
    - StrategyAggressive,
    - StrategyEndMarker

  */

  enum StrategyMode {
    StrategyDefault = 0,
    StrategyEasy,
    StrategyAggressive,
    StrategyEndMarker
  };

  /*! \brief What clock type to use

    - Clock Cpu: Use CPU time
    - ClockElapsed: Use elapsed time

  */

  enum ClockType { ClockCpu = 0, ClockElapsed, ClockEndMarker };

  /*! \brief What clock type to use

    - CGraphOff:
    - CGraphOn:
    - CGraphClique:
    - CGraphEndMarker

  */

  enum CGraphMode { CGraphOff = 0, CGraphOn, CGraphClique, CGraphEndMarker };

  /*
      In order for initialisation to work properly, the order of declaration of
      the enum constants here must match the order of keyword declaration for
      the various cut and heuristic control parameters in
      CbcGenParamUtils::addCbcGenParams
    */
  /*! \brief Codes to control the use of cut generators and heuristics

      - CGOff: the cut generator will not be installed
      - CGOn:  the cut generator will be installed; exactly how often it's
           activated depends on the settings at installation
      - CGRoot: the cut generator will be installed with settings that restrict
           it to activation at the root node only.
      - CGIfMove: the cut generator will be installed with settings that allow
           it to remain active only so long as it's generating cuts that
           tighten the relaxation.
      - CGForceOn: the cut generator will be installed with settings that force
           it to be called at every node
      - CGForceOnBut: the cut generator will be installed with settings that
     force it to be called at every node, but more active at root (probing only)
      - CGEndOnly:
      - CGEndClean:
      - CGEndBoth:
      - CGLlong:
      - CGLongRoot:
      - CGLongIfMove:
      - CGForceLongOn:
      - CGLongEndOnly:
      - CGOnlyAsWell:
      - CGOnlyAsWellRoot:
      - CGCleanAsWell:
      - CGCleanAsWellRoot:
      - CGBothAsWell:
      - CGBothAsWellRoot:
      - CGOnlyInstead:
      - CGCleanInstead:
      - CGBothInstead:
      - CGOnGlobal:
      - CGForceAndGlobal:
      - CGShorter:
      - CGEndMarker: a convenience to mark the end of the codes.

    */

  enum CGMode {
    CGOff = 0,
    CGOn,
    CGRoot,
    CGIfMove,
    CGForceOn,
    CGForceOnBut,
    CGEndOnly,
    CGEndOnlyRoot,
    CGEndClean,
    CGEndCleanRoot,
    CGEndBoth,
    CGLong,
    CGLongOn,
    CGLongRoot,
    CGLongIfMove,
    CGForceLongOn,
    CGLongEndOnly,
    CGOnlyAsWell,
    CGOnlyAsWellRoot,
    CGCleanAsWell,
    CGCleanAsWellRoot,
    CGBothAsWell,
    CGBothAsWellRoot,
    CGCleanBothAsWellRoot,
    CGOnlyInstead,
    CGCleanInstead,
    CGBothInstead,
    CGOnGlobal,
    CGForceAndGlobal,
    CGLonger,
    CGShorter,
    CGEndMarker
  };

  /*! \brief Codes to specify whether to use a cutoff constraint
   */

  enum CutoffMode {
    COOff = 0,
    COOn,
    COVariable,
    COForceVariable,
    COConflict,
    COEndMarker
  };

  /*! \brief Codes to specify the assignment of branching priorities

      - HeurOff:
      - HeurOn:
      - HeurBoth,
      - HeurBefore,
      - HuerOften,
      - HeurTen,
      - HeurOneHundred,
      - HeurTwoHundred,
      - HeurThreeHundred,
      - HeurOneThousand,
      - HeurTenThousand,
      - HeurDj,
      - HeurDjBefore,
      - HeurUseSolution,
      - HeurInTree,
    */

  enum HeurMode {
    HeurOff = 0,
    HeurOn,
    HeurRoot,
    HeurBoth,
    HeurBefore,
    HeurOften,
    HeurTen,
    HeurOneHundred,
    HeurTwoHundred,
    HeurThreeHundred,
    HeurOneThousand,
    HeurTenThousand,
    HeurDj,
    HeurDjBefore,
    HeurUseSolution,
    HeurInTree,
    HeurEndMarker
  };

  /*! \brief Codes to specify one or off for binary parameters

     - ParamOff: Capability is switched off
     - ParamOn: Capability is switched on
   */

  enum OnOffMode { ParamOff = 0, ParamOn, ParamEndMarker };

  /*! \brief Codes to specify the assignment of branching priorities

      - BPOff: no priorities are passed to cbc
      - BPCost: a priority vector is constructed based on objective coefficients
      - BPOrder: a priority vector is constructed based on column order
      - BPExt: the user has provided a priority vector
    */

  enum BPMode { BPOff = 0, BPCost, BPOrder, BPExt, BPEndMarker };

  /*! \brief Codes tos pecify mode for printing integers

    - PMNormal = 0,
    - PMInteger,
    - PMSpecial,
    - PMRows,
    - PMAll,
    - PMEndMarker

  */

  enum IntPrintMode {
    PMNormal = 0,
    PMInteger,
    PMSpecial,
    PMRows,
    PMAll,
    PMEndMarker
  };

  /*! \brief Major status codes for branch-and-cut

      - BACInvalid: status not yet set
      - BACNotRun: branch-and-cut has not yet run for the current problem
      - BACFinish: branch-and-cut has finished normally
      - BACStop: branch-and-cut has stopped on a limit
      - BACAbandon: branch-and-cut abandoned the problem
      - BACUser: branch-and-cut stopped on user signal

      Consult minorStatus_ for details.

      These codes are (mostly) set to match the codes used by CbcModel.
     Additions to CbcModel codes should be reflected here and in translateMajor.
    */

  enum BACMajorStatus {
    BACInvalid = -1,
    BACFinish = 0,
    BACStop = 1,
    BACAbandon = 2,
    BACNotRun,
    BACUser = 5,
    BacEndMarker
  };

  /*! \brief Minor status codes

      - BACmInvalid		status not yet set
      - BACmFinish		search exhausted the tree; optimal solution
     found
      - BACmInfeas		problem is infeasible
      - BACmUbnd		problem is unbounded
      - BACmGap		stopped on integrality gap
      - BACmNodeLimit	stopped on node limit
      - BACmTimeLimit	stopped on time limit
      - BACmSolnLimit	stopped on number of solutions limit
      - BACmUser		stopped due to user event
      - BACmOther		nothing else is appropriate

      It's not possible to make these codes agree with CbcModel. The meaning
     varies according to context: if the BACWhere code specifies a relaxation,
     then the minor status reflects the underlying OSI solver. Otherwise, it
     reflects the integer problem.
    */

  enum BACMinorStatus {
    BACmInvalid = -1,
    BACmFinish = 0,
    BACmInfeas,
    BACmUbnd,
    BACmGap,
    BACmNodeLimit,
    BACmTimeLimit,
    BACmSolnLimit,
    BACmUser,
    BACmOther,
    BACmEndMarker
  };

  /*! \brief Codes to specify where branch-and-cut stopped

      - BACwNotStarted	stopped before we ever got going
      - BACwBareRoot	stopped after initial solve of root relaxation
      - BACwIPP		stopped after integer preprocessing
      - BACwIPPRelax	stopped after initial solve of preprocessed problem
      - BACwBAC		stopped at some point in branch-and-cut
    */

  enum BACWhere {
    BACwInvalid = -1,
    BACwNotStarted = 0,
    BACwBareRoot,
    BACwIPP,
    BACwIPPRelax,
    BACwBAC,
    BACwEndMarker
  };

  enum OptimizationDirection {
    OptDirMaximize = 0,
    OptDirMinimize,
    OptDirZero,
    OptDireEndMarker
  };

  //@}

  /*! \name Functions for Setting Up Parameters
      \brief Functions that populate the parameters objects. 
  */

  //@{

  /*! set up the solver parameter vector */
  void addCbcParams();
  void addCbcSolverStrParams();
  void addCbcSolverHelpParams();
  void addCbcSolverActionParams();
  void addCbcSolverKwdParams();
  void addCbcSolverDblParams();
  void addCbcSolverIntParams();
  void addCbcSolverBoolParams();
  void addCbcSolverCutParams();
  void addCbcSolverHeurParams();
  void addCbcModelParams();

  //@}

  /*! \name Access and Control Functions for Cut Generators
        \brief Control functions, plus lazy creation functions for cut
     generators

        Cbc avoids creating objects for cut generators and heuristics
        unless they're actually used. For cut generators, a prototype is created
        and reused. For heuristics, the default is to create a new object with
     each call, because the model may have changed. The object is returned
     through the reference parameter. The return value of the function is the
     current mode.

        Cut generator and heuristic objects created by these calls will be
     deleted with the destruction of the CbcSettings object.
    */
  //@{

  /*! \brief Get cut depth setting

      The name is a bit of a misnomer. Essentially, this overrides the
      `every so many nodes' control with `execute when (depth in tree)
      mod (cut depth) == 0'.
    */
  inline int getCutDepth() { return cutDepth_; }

  /*! \brief Get cut length setting */
  inline int getCutLength() { return cutLength_; }

  /*! \brief Get cut length setting */
  inline int getCutPassInTree() { return cutPassInTree_; }

  /*! \brief Set cut depth setting.

      See comments for getCutDepth().
  */
  inline void setCutDepth(int cutDepth) { cutDepth_ = cutDepth; }

  /*! \brief Set cut length setting. */
  inline void setCutLength(int cutLength) { cutLength_ = cutLength; }

  /*! \brief Set cut pass in tree setting. */
  inline void setCutPassInTree(int cutPassInTree) {
    cutPassInTree_ = cutPassInTree;
  }

  /*! \brief Obtain a prototype for a clique cut generator. */
  CbcSettings::CGMode getClique(CglCutGenerator *&gen);

  /*! \brief Set mode for use of clique cut generator. */
  inline void setCliqueMode(CbcSettings::CGMode mode) {
    clique_.mode_ = mode;
  }

  /*! \brief Get mode for use of clique cut generator. */
  inline CbcSettings::CGMode getCliqueMode() { return (clique_.mode_); }

  /*! \brief Obtain a prototype for a flow cover cut generator. */
  CbcSettings::CGMode getFlow(CglCutGenerator *&gen);

  /*! \brief Set mode for use of flow cover cut generator. */
  inline void setFlowMode(CbcSettings::CGMode mode) { flow_.mode_ = mode; }

  /*! \brief Get mode for use of flow cover cut generator. */
  inline CbcSettings::CGMode getFlowMode() { return (flow_.mode_); }

  /*! \brief Obtain a prototype for a GMI cut generator. */
  CbcSettings::CGMode getGMI(CglCutGenerator *&gen);

  /*! \brief Set mode for use of GMI cut generator. */
  inline void setGMIMode(CbcSettings::CGMode mode) { gmi_.mode_ = mode; }

  /*! \brief Get mode for use of GMI cut generator. */
  inline CbcSettings::CGMode getGMIMode() { return (gmi_.mode_); }

  /*! \brief Obtain a prototype for a Gomory cut generator. */
  CbcSettings::CGMode getGomory(CglCutGenerator *&gen);

  /*! \brief Set mode for use of Gomory cut generator. */
  inline void setGomoryMode(CbcSettings::CGMode mode) {
    gomory_.mode_ = mode;
  }

  /*! \brief Get mode for use of Gomory cut generator. */
  inline CbcSettings::CGMode getGomoryMode() { return (gomory_.mode_); }

  /*! \brief Obtain a prototype for a knapsack cover cut generator. */
  CbcSettings::CGMode getKnapsack(CglCutGenerator *&gen);

  /*! \brief Set mode for use of knapsack cut generator. */
  inline void setKnapsackMode(CbcSettings::CGMode mode) {
    knapsack_.mode_ = mode;
  }

  /*! \brief Get mode for use of knapsack cut generator. */
  inline CbcSettings::CGMode getKnapsackMode() { return (knapsack_.mode_); }

  /*! \brief Obtain a prototype for a LaGomory cut generator. */
  CbcSettings::CGMode getLaGomory(CglCutGenerator *&gen);

  /*! \brief Set mode for use of LaGomory cut generator. */
  inline void setLaGomoryMode(CbcSettings::CGMode mode) {
    laGomory_.mode_ = mode;
  }

  /*! \brief Get mode for use of LaGomory cut generator. */
  inline CbcSettings::CGMode getLaGomoryMode() { return (laGomory_.mode_); }

  /*  \brief Obtain a prototype for a lift-and-project cut generator. */
  CbcSettings::CGMode getLandP(CglCutGenerator *&gen);

  /* \brief Set mode for use of lift-and-project cut generator. */
  inline void setLandPMode(CbcSettings::CGMode mode) { landP_.mode_ = mode; }

  /*! \brief Get mode for use of lift-and-project cut generator. */
  inline CbcSettings::CGMode getLandPMode() { return (landP_.mode_); }

  /*! \brief Obtain a prototype for a LaTwomir cut generator. */
  CbcSettings::CGMode getLaTwomir(CglCutGenerator *&gen);

  /*! \brief Set mode for use of LaTwomir cut generator. */
  inline void setLaTwomirMode(CbcSettings::CGMode mode) {
    laTwomir_.mode_ = mode;
  }

  /*! \brief Get mode for use of LaTwomir cut generator. */
  inline CbcSettings::CGMode getLaTwomirMode() { return (laTwomir_.mode_); }

  /*! \brief Obtain a prototype for a mixed integer rounding (MIR)
          cut generator.
   */
  CbcSettings::CGMode getMir(CglCutGenerator *&gen);

  /*! \brief Set mode for use of MIR cut generator. */
  inline void setMirMode(CbcSettings::CGMode mode) { mir_.mode_ = mode; }

  /*! \brief Get mode for use of MIR cut generator. */
  inline CbcSettings::CGMode getMirMode() { return (mir_.mode_); }

  /*! \brief Obtain a prototype for an odd wheel cut generator. */
  CbcSettings::CGMode getOddWheel(CglCutGenerator *&gen);

  /*! \brief Set mode for use of odd wheel cut generator. */
  inline void setOddWheelMode(CbcSettings::CGMode mode) {
    oddWheel_.mode_ = mode;
  }

  /*! \brief Get mode for use of odd wheel cut generator. */
  inline CbcSettings::CGMode getOddWheelMode() { return (oddWheel_.mode_); }

  /*! \brief Obtain a prototype for a probing cut generator. */
  CbcSettings::CGMode getProbing(CglCutGenerator *&gen);

  /*! \brief Set mode for use of probing cut generator. */
  inline void setProbingMode(CbcSettings::CGMode mode) {
    probing_.mode_ = mode;
  }

  /*! \brief Get mode for use of probing cut generator. */
  inline CbcSettings::CGMode getProbingMode() { return (probing_.mode_); }

  /*! \brief Obtain a prototype for a reduce and split cut generator. */
  CbcSettings::CGMode getRedSplit(CglCutGenerator *&gen);

  /*! \brief Set mode for use of reduce and split cut generator. */
  inline void setRedSplitMode(CbcSettings::CGMode mode) {
    redSplit_.mode_ = mode;
  }

  /*! \brief Get mode for use of reduce and split cut generator. */
  inline CbcSettings::CGMode getRedSplitMode() { return (redSplit_.mode_); }

  /*! \brief Obtain a prototype for a reduce and split 2 cut generator. */
  CbcSettings::CGMode getRedSplit2(CglCutGenerator *&gen);

  /*! \brief Set mode for use of reduce and split 2 cut generator. */
  inline void setRedSplit2Mode(CbcSettings::CGMode mode) {
    redSplit2_.mode_ = mode;
  }

  /*! \brief Get mode for use of reduce and split 2 cut generator. */
  inline CbcSettings::CGMode getRedSplit2Mode() {
    return (redSplit2_.mode_);
  }

  /*! \brief Obtain a prototype for a residual capacity cut generator. */
  CbcSettings::CGMode getResidCap(CglCutGenerator *&gen);

  /*! \brief Set mode for use of residual capacity cut generator. */
  inline void setResidCapMode(CbcSettings::CGMode mode) {
    residCap_.mode_ = mode;
  }

  /*! \brief Get mode for use of residual capacity cut generator. */
  inline CbcSettings::CGMode getResidCapMode() { return (residCap_.mode_); }

  /*! \brief Obtain a prototype for a 2-MIR cut generator. */
  CbcSettings::CGMode getTwomir(CglCutGenerator *&gen);

  /*! \brief Set mode for use of 2-MIR cut generator. */
  inline void setTwomirMode(CbcSettings::CGMode mode) {
    twomir_.mode_ = mode;
  }

  /*! \brief Get mode for use of Twomir cut generator. */
  inline CbcSettings::CGMode getTwomirMode() { return (twomir_.mode_); }

  /*! \brief Obtain a prototype for a zero-half cut generator. */
  CbcSettings::CGMode getZeroHalf(CglCutGenerator *&gen);

  /*! \brief Set mode for use of zero-half cut generator. */
  inline void setZeroHalfMode(CbcSettings::CGMode mode) {
    zeroHalf_.mode_ = mode;
  }

  /*! \brief Get mode for use of reduce and split cut generator. */
  inline CbcSettings::CGMode getZeroHalfMode() { return (zeroHalf_.mode_); }

  //@}

  /*! \name Access and Control Functions for Heuristics
        \brief Control functions, plus lazy creation functions for cut
     generators

        CBC avoids creating objects for heuristics
        unless they're actually used. For cut generators, a prototype is created
        and reused. For heuristics, the default is to create a new object with
     each call, because the model may have changed. The object is returned
     through the reference parameter. The return value of the function is the
     current mode.

        Heuristic objects created by these calls will be deleted
        with the destruction of the CbcSettings object.
  */
  //@{

  /*! \brief Set mode for overall control of heuristics. */
  inline CbcSettings::OnOffMode getDoHeuristicMode() {
    return (doHeuristicMode_);
  }

  /*! \brief Get mode for overall control of heuristics. */
  inline void
  setDoHeuristicMode(CbcSettings::OnOffMode mode) {
    doHeuristicMode_ = mode;
  }

  /*! \brief Obtain a local search/combine heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getCombine(CbcHeuristic *&gen,
                                      bool alwaysCreate = true);

  /*! \brief Set mode for use of local search/combine heuristic. */
  inline void setCombineMode(CbcSettings::HeurMode mode) {
    combine_.mode_ = mode;
  }

  /*! \brief Get mode for use of local search/combine heuristic. */
  CbcSettings::HeurMode getCombineMode() { return (combine_.mode_); }

  /*! \brief Obtain a crossover heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getCrossover(CbcHeuristic *&gen, 
                                        bool alwaysCreate = true);

  /*! \brief Set mode for use of crossover heuristic. */
  inline void setCrossoverMode(CbcSettings::HeurMode mode) {
    crossover_.mode_ = mode;
  }

  /*! \brief Get mode for use of crossover heuristic. */
  CbcSettings::HeurMode getCrossoverMode() { return (crossover_.mode_); }

  /*! \brief Obtain a DINS heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getDins(CbcHeuristic *&gen, 
                                   bool alwaysCreate = true);

  /*! \brief Set mode for use of DINS heuristic. */
  inline void setDinsMode(CbcSettings::HeurMode mode) { dins_.mode_ = mode; }

  /*! \brief Get mode for use of DINS heuristic. */
  CbcSettings::HeurMode getDinsMode() { return (dins_.mode_); }

  /*! \brief Obtain a Diving Coefficient heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getDiveCofficient(CbcHeuristic *&gen,
                                             bool alwaysCreate = true);

  /*! \brief Set mode for use of Diving Coefficient heuristic. */
  inline void setDiveCoefficientMode(CbcSettings::HeurMode mode) {
    divingc_.mode_ = mode;
  }

  /*! \brief Get mode for use of Diving Coefficient heuristic. */
  CbcSettings::HeurMode getDiveCoefficientMode() { return (divingc_.mode_); }

  /*! \brief Obtain a Diving Fractional heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getDiveFractional(CbcHeuristic *&gen,
                                             bool alwaysCreate = true);

  /*! \brief Set mode for use of Diving Coefficient heuristic. */
  inline void setDiveFractionalMode(CbcSettings::HeurMode mode) {
    divingf_.mode_ = mode;
  }

  /*! \brief Get mode for use of Diving Fractional heuristic. */
  CbcSettings::HeurMode getDiveFractionalMode() { return (divingf_.mode_); }

  /*! \brief Obtain a Diving Guided heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getDiveGuided(CbcHeuristic *&gen, 
                                         bool alwaysCreate = true);

  /*! \brief Set mode for use of Diving Guided heuristic. */
  inline void setDiveGuidedMode(CbcSettings::HeurMode mode) {
    divingg_.mode_ = mode;
  }

  /*! \brief Get mode for use of Diving Guided heuristic. */
  CbcSettings::HeurMode getDiveGuidedMode() { return (divingg_.mode_); }

  /*! \brief Obtain a Diving Line Search heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getDiveLineSearch(CbcHeuristic *&gen,
                                             bool alwaysCreate = true);

  /*! \brief Set mode for use of Diving Line Search heuristic. */
  inline void setDiveLineSearchMode(CbcSettings::HeurMode mode) {
    divingl_.mode_ = mode;
  }

  /*! \brief Get mode for use of Diving LineSearch heuristic. */
  CbcSettings::HeurMode getDiveLineSearchMode() { return (divingl_.mode_); }

  /*! \brief Obtain a Diving Pseudocost heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getDivePseudocost(CbcHeuristic *&gen,
                                             bool alwaysCreate = true);

  /*! \brief Set mode for use of Diving Pseudocost heuristic. */
  inline void setDivePseudocostMode(CbcSettings::HeurMode mode) {
    divingp_.mode_ = mode;
  }

  /*! \brief Get mode for use of Diving Pseudocost heuristic. */
  CbcSettings::HeurMode getDivePseudocostMode() { return (divingp_.mode_); }

  /*! \brief Set mode for use of rand. */
  inline void setDiveRandomMode(CbcSettings::HeurMode mode) {
    randomDivingMode_ = mode;
  }

  /*! \brief Get mode for use of Diving Pseudocost heuristic. */
  CbcSettings::HeurMode getDiveRandomMode() { return (randomDivingMode_); }

  /*! \brief Obtain a Diving Vector Length heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getDiveVectorLength(CbcHeuristic *&gen,
                                               bool alwaysCreate = true);

  /*! \brief Set mode for use of Diving Vector Length heuristic. */
  inline void setDiveVectorLengthMode(CbcSettings::HeurMode mode) {
    divingv_.mode_ = mode;
  }

  /*! \brief Get mode for use of Diving Vector Length heuristic. */
  CbcSettings::HeurMode getDiveVectorLengthMode() {
    return (divingv_.mode_);
  }

  /*! \brief Obtain a DW heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getDW(CbcHeuristic *&gen, 
                                 bool alwaysCreate = true);

  /*! \brief Set mode for use of DW heuristic. */
  inline void setDWMode(CbcSettings::HeurMode mode) { dw_.mode_ = mode; }

  /*! \brief Get mode for use of DW heuristic. */
  CbcSettings::HeurMode getDWMode() { return (dw_.mode_); }

  /*! \brief Obtain a feasibility pump heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
  */
  CbcSettings::HeurMode getFeasPump(CbcHeuristic *&gen, 
                                       bool alwaysCreate = true);

  /*! \brief Set mode for use of feasibility pump heuristic. */
  inline void setFeasPumpMode(CbcSettings::HeurMode mode) {
    fpump_.mode_ = mode;
  }

  /*! \brief Get mode for use of feasibility pump heuristic. */
  CbcSettings::HeurMode getFeasPumpMode() { return (fpump_.mode_); }

  /*! \brief Set iterations for use of feasibility pump heuristic. */
  inline void setFeasPumpIters(int iters) { fpump_.iters_ = iters; }

  /*! \brief Get iterations for use of feasibility pump heuristic. */
  inline int getFeasPumpIters() { return (fpump_.iters_); }

  /*! \brief Set tune mode for use of feasibility pump heuristic. */
  inline void setFeasPumpTune(int tune) { fpump_.tune_ = tune; }

  /*! \brief Get tune mode for use of feasibility pump heuristic. */
  inline int getFeasPumpTune() { return (fpump_.tune_); }

  /*! \brief Set second tune mode for use of feasibility pump heuristic. */
  inline void setFeasPumpTune2(int tune2) { fpump_.tune2_ = tune2; }

  /*! \brief Get second tune mode for use of feasibility pump heuristic. */
  inline int getFeasPumpTune2(int tune2) { return (fpump_.tune2_); }

  /*! \brief Set fake cutoff for use of feasibility pump heuristic. */
  inline void setFeasPumpFakeCutoff(double cutoff) { fpump_.cutoff_ = cutoff; }

  /*! \brief Get fake cutoff for use of feasibility pump heuristic. */
  inline double getFeasPumpFakeCutoff() { return (fpump_.cutoff_); }

  /*! \brief Set fake increment for use of feasibility pump heuristic. */
  inline void setFeasPumpFakeIncrement(double increment) {
    fpump_.increment_ = increment;
  }

  /*! \brief Get fake increment for use of feasibility pump heuristic. */
  inline double getFeasPumpFakeIncrement() { return (fpump_.increment_); }

  /*! \brief Obtain a greedy cover heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getGreedyCover(CbcHeuristic *&gen, 
                                          bool alwaysCreate = true);

  /*! \brief Set mode for use of greedy cover heuristic. */
  inline void setGreedyCoverMode(CbcSettings::HeurMode mode) {
    greedyCover_.mode_ = mode;
  }

  /*! \brief Get mode for use of greedy cover heuristic. */
  CbcSettings::HeurMode getGreedyCoverMode() { return (greedyCover_.mode_); }

  /*! \brief Obtain a greedy equality heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getGreedyEquality(CbcHeuristic *&gen,
                                             bool alwaysCreate = true);

  /*! \brief Set mode for use of greedy equality heuristic. */
  inline void setGreedyEqualityMode(CbcSettings::HeurMode mode) {
    greedyEquality_.mode_ = mode;
  }

  /*! \brief Get mode for use of greedy equality heuristic. */
  CbcSettings::HeurMode getGreedyEqualityMode() {
    return (greedyEquality_.mode_);
  }

  /*! \brief Obtain a Naive heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
  */
  CbcSettings::HeurMode getNaiveHeur(CbcHeuristic *&gen,
                                        bool alwaysCreate = true);

  /*! \brief Set mode for use of Naive heuristic. */
  inline void setNaiveHeurMode(CbcSettings::HeurMode mode) {
    naive_.mode_ = mode;
  }

  /*! \brief Get mode for use of Naive heuristic. */
  CbcSettings::HeurMode getNaiveHeurMode() { return (divingc_.mode_); }

  /*! \brief Obtain a Pivot And Fix heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getPivotAndFix(CbcHeuristic *&gen, 
                                          bool alwaysCreate = true);

  /*! \brief Set mode for use of Pivot and Fix heuristic. */

  inline void setPivotAndFixMode(CbcSettings::HeurMode mode) {
    pivotAndFix_.mode_ = mode;
  }

  /*! \brief Get mode for use of pivot and fix heuristic. */
  CbcSettings::HeurMode getPivotAndFixMode() { return (pivotAndFix_.mode_); }

#if 0
  /*! \brief Obtain a Pivot and Complement heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */
  CbcSettings::HeurMode getPivotAndComplement(CbcHeuristic *&gen,
                                                 bool alwaysCreate = true);

  /*! \brief Set mode for use of Pivot and Complement heuristic. */
  inline void setPivotAndComplementMode(CbcSettings::HeurMode mode)
  {
    pivotandcomplement_.mode_ = mode;
  }

  /*! \brief Get mode for use of pivot and complement heuristic. */
  CbcSettings::HeurMode getPivotAndComplementMode()
  {
     return(pivotAndComplement_.mode_);
  }
#endif

  /*! \brief Obtain a Proximity heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getProximity(CbcHeuristic *&gen, 
                                        bool alwaysCreate = true);

  /*! \brief Set mode for use of Proximity heuristic. */
  inline void setProximityMode(CbcSettings::HeurMode mode) {
    proximity_.mode_ = mode;
  }

  /*! \brief Get mode for use of proximity heuristic. */
  CbcSettings::HeurMode getProximityMode() { return (proximity_.mode_); }

  /*! \brief Obtain a Randomized Rounding heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getRandRound(CbcHeuristic *&gen,
                                        bool alwaysCreate = true);

  /*! \brief Set mode for use of Randomized Rounding heuristic. */
  inline void setRandRoundMode(CbcSettings::HeurMode mode) {
    randRound_.mode_ = mode;
  }

  /*! \brief Get mode for use of randomized rounding heuristic. */
  CbcSettings::HeurMode getRandRoundMode() { return (randRound_.mode_); }

  /*! \brief Obtain a RENS heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getRens(CbcHeuristic *&gen,
                                   bool alwaysCreate = true);

  /*! \brief Set mode for use of RENS heuristic. */
  inline void setRensMode(CbcSettings::HeurMode mode) { rens_.mode_ = mode; }

  /*! \brief Get mode for use of RENS heuristic. */
  CbcSettings::HeurMode getRensMode() { return (rens_.mode_); }

  /*! \brief Obtain a RINS heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getRins(CbcHeuristic *&gen, 
                                   bool alwaysCreate = true);

  /*! \brief Set mode for use of RINS heuristic. */
  inline void setRinsMode(CbcSettings::HeurMode mode) { rins_.mode_ = mode; }

  /*! \brief Get mode for use of RINS heuristic. */
  CbcSettings::HeurMode getRinsMode() { return (rins_.mode_); }

  /*! \brief Obtain a simple rounding heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getRounding(CbcHeuristic *&gen,
                                       bool alwaysCreate = true);

  /*! \brief Set mode for use of simple rounding heuristic. */
  inline void setRoundingMode(CbcSettings::HeurMode mode) {
    rounding_.mode_ = mode;
  }

  /*! \brief Get mode for use of rounding heuristic. */
  CbcSettings::HeurMode getRoundingMode() { return (rounding_.mode_); }

  /*! \brief Obtain a variable neighborhood heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getVnd(CbcHeuristic *&gen,
                                  bool alwaysCreate = true);

  /*! \brief Set mode for use of variable neighborhood heuristic. */
  inline void setVndMode(CbcSettings::HeurMode mode) { vnd_.mode_ = mode; }

  /*! \brief Get mode for use of variable neighborhood heuristic. */
  CbcSettings::HeurMode getVndMode() { return (vnd_.mode_); }

  /*! \brief Obtain a local search tree object

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing
     object if one exists.
    */
  CbcSettings::HeurMode getLocalTree(CbcTreeLocal *&localTree,
                                        bool alwaysCreate = true);

  /*! \brief Set mode for use of local tree. */
  inline void setLocalTreeMode(CbcSettings::HeurMode mode) {
    localTree_.mode_ = mode;
  }

  /*! \brief Get mode for use of Diving Coefficient heuristic. */
  CbcSettings::HeurMode getLocalTreeMode() { return (localTree_.mode_); }

  //@}

  /*! \name Miscellaneous Integer Parameters

  */

  //@{

  /*! \brief Get BkPivotStrategy setting */
  inline int getBkPivotStrategy() { return bkPivotStrategy_; }

  /*! \brief Set BkPivotStrategy setting */
  inline void setBkPivotStrategy(int bkPivotStrategy) {
    bkPivotStrategy_ = bkPivotStrategy;
  }

  /*! \brief Get BkMaxCalls setting */
  inline int getBkMaxCalls() { return bkMaxCalls_; }

  /*! \brief Set BkMaxCalls setting */
  inline void setBkMaxCalls(int bkMaxCalls) { bkMaxCalls_ = bkMaxCalls; }

  /*! \brief Get BkClqExtMethod setting */
  inline int getBkClqExtMethod() { return bkClqExtMethod_; }

  /*! \brief Set BkClqExtMethod setting */
  inline void setBkClqExtMethod(int bkClqExtMethod) {
    bkClqExtMethod_ = bkClqExtMethod;
  }

  /*! \brief Get CppMode setting */
  inline int getCppMode() { return cppMode_; }

  /*! \brief Set CppMode setting */
  inline void setCppMode(int cppMode) { cppMode_ = cppMode; }

  /*! \brief Get DepthMiniBaB setting */
  inline int getDepthMiniBaB() { return depthMiniBaB_; }

  /*! \brief Set DepthMiniBaB setting */
  inline void setDepthMiniBaB(int depthMiniBaB) {
    depthMiniBaB_ = depthMiniBaB;
  }

  /*! \brief Get DiveOpt setting */
  inline int getDiveOpt() { return diveOpt_; }

  /*! \brief Set DiveOpt setting */
  inline void setDiveOpt(int diveOpt) { diveOpt_ = diveOpt; }

  /*! \brief Get DiveOptSolves setting */
  inline int getDiveOptSolves() { return diveOptSolves_; }

  /*! \brief Set DiveOptSolves setting */
  inline void setDiveOptSolves(int diveOptSolves) {
    diveOptSolves_ = diveOptSolves;
  }

  /*! \brief Get ExperimentMode setting */
  inline int getExperimentMode() { return experiment_; }

  /*! \brief Set ExperimentMode setting */
  inline void setExperimentMode(int experimentMode) {
    experiment_ = experimentMode;
  }

  /*! \brief Get ExtraIntParam1 setting */
  inline int getExtraIntParam1() { return extraInt1_; }

  /*! \brief Set ExtraIntParam1 setting */
  inline void setExtraIntParam1(int extraInt1) { extraInt1_ = extraInt1; }

  /*! \brief Get ExtraIntParam2 setting */
  inline int getExtraIntParam2() { return extraInt2_; }

  /*! \brief Set ExtraIntParam2 setting */
  inline void setExtraIntParam2(int extraInt2) { extraInt2_ = extraInt2; }

  /*! \brief Get ExtraIntParam3 setting */
  inline int getExtraIntParam3() { return extraInt3_; }

  /*! \brief Set ExtraIntParam3 setting */
  inline void setExtraIntParam3(int extraInt3) { extraInt3_ = extraInt3; }

  /*! \brief Get ExtraIntParam4 setting */
  inline int getExtraIntParam4() { return extraInt4_; }

  /*! \brief Set ExtraIntParam4 setting */
  inline void setExtraIntParam4(int extraInt4) { extraInt4_ = extraInt4; }

  /*! \brief Get HeurOptions setting */
  inline int getHeurOptions() { return heurOptions_; }

  /*! \brief Set HeurOptions setting */
  inline void setHeurOptions(int heurOptions) { heurOptions_ = heurOptions; }

  /*! \brief Get MaxSavedSols setting */
  inline int getMaxSavedSols() { return maxSavedSols_; }

  /*! \brief Set MaxSavedSols setting */
  inline void setMaxSavedSols(int maxSavedSols) {
    maxSavedSols_ = maxSavedSols;
  }

  /*! \brief Get MaxSlowCuts setting */
  inline int getMaxSlowCuts() { return maxSlowCuts_; }

  /*! \brief Set MaxSlowCuts setting */
  inline void setMaxSlowCuts(int maxSlowCuts) { maxSlowCuts_ = maxSlowCuts; }

  /*! \brief Get MoreMoreOptions setting */
  inline int getMoreMoreOptions() { return moreMoreOptions_; }

  /*! \brief Set MoreMoreOptions setting */
  inline void setMoreMoreOptions(int moreMoreOptions) {
    moreMoreOptions_ = moreMoreOptions;
  }

  /*! \brief Get MultipleRoots setting */
  inline int getMultipleRoots() { return multipleRoots_; }

  /*! \brief Set MultipleRoots setting */
  inline void setMultipleRoots(int multipleRoots) {
    multipleRoots_ = multipleRoots;
  }

  /*! \brief Get OddWextMethod setting */
  inline int getOddWextMethod() { return oddWextMethod_; }

  /*! \brief Set OddWextMethod setting */
  inline void setOddWextMethod(int oddWextMethod) {
    oddWextMethod_ = oddWextMethod;
  }

  /*! \brief Get OutputFormat setting */
  inline int getOutputFormat() { return outputFormat_; }

  /*! \brief Set OutputFormat setting */
  inline void setOutputFormat(int outputFormat) {
    outputFormat_ = outputFormat;
  }

  /*! \brief Get PrintOptions setting */
  inline int getPrintOptions() { return printOpt_; }

  /*! \brief Set PrintOptions setting */
  inline void setPrintOptions(int printOpt) { printOpt_ = printOpt; }

  /*! \brief Enable printing */
  inline void enablePrinting() { noPrinting_ = false; }
   
  /*! \brief Disable printing */
  inline void disablePrinting() { noPrinting_ = true; }

  /*! \brief Get printing status*/
  inline bool noPrinting() { return noPrinting_; }

  /*! \brief Enable welcome printing */
  inline void enableWelcomePrinting() { printWelcome_ = true; }
   
  /*! \brief Disable welcome printing */
  inline void disableWelcomePrinting() { printWelcome_ = false; }

  /*! \brief Get welcome printing status */
  inline bool printWelcome() { return printWelcome_; }

  /*! \brief Enable signal handler */
  inline void enableSignalHandler() { useSignalHandler_ = true; }
   
  /*! \brief Disable signal handler */
  inline void disableSignalHandler() { useSignalHandler_ = false; }

  /*! \brief Get signal handler status*/
  inline bool useSignalHandler() { return useSignalHandler_; }

  /*! \brief Get ProcessTune setting */
  inline int getProcessTune() { return processTune_; }

  /*! \brief Set ProcessTune setting */
  inline void setProcessTune(int processTune) { processTune_ = processTune; }

  /*! \brief Get RandomSeed setting */
  inline int getRandomSeed() { return randomSeed_; }

  /*! \brief Set RandomSeed setting */
  inline void setRandomSeed(int randomSeed) { randomSeed_ = randomSeed; }

  /*! \brief Get StrongStrategy setting */
  inline int getStrongStrategy() { return strongStrategy_; }

  /*! \brief Set StrongStrategy setting */
  inline void setStrongStrategy(int strongStrategy) {
    strongStrategy_ = strongStrategy;
  }

  /*! \brief Get TestOsi setting */
  inline int getTestOsi() { return testOsi_; }

  /*! \brief Set TestOsi setting */
  inline void setTestOsi(int testOsi) { testOsi_ = testOsi; }

  /*! \brief Get Threads setting */
  inline int getThreads() { return threads_; }

  /*! \brief Set Threads setting */
  inline void setThreads(int threads) { threads_ = threads; }

  /*! \brief Get UserCbc setting */
  inline int getUserCbc() { return userCbc_; }

  /*! \brief Set UserCbc setting */
  inline void setUserCbc(int userCbc) { userCbc_ = userCbc; }

  /*! \brief Get Verbose setting */
  inline int getVerbose() { return verbose_; }

  /*! \brief Set Verbose setting */
  inline void setVerbose(int verbose) { verbose_ = verbose; }

  /*! \brief Get VubTry setting */
  inline int getVubTry() { return vubTry_; }

  /*! \brief Set VubTry setting */
  inline void setVubTry(int vubTry) { vubTry_ = vubTry; }

  //@}

  /*! \name Miscellaneous Double Parameters

  */

  //@{

  /*! \brief Get extra double 3 setting */
  inline double getExtraDbl3() { return extraDbl3_; }

  /*! \brief Set extra double 3 setting */
  inline void setExtraDbl3(double extraDbl3) { extraDbl3_ = extraDbl3; }

  /*! \brief Get extra double 4 setting */
  inline double getExtraDbl4() { return extraDbl4_; }

  /*! \brief Set extra double 4 setting */
  inline void setExtraDbl4(double extraDbl4) { extraDbl4_ = extraDbl4; }

  /*! \brief Get extra double 3 setting */
  inline double getExtraDbl5() { return extraDbl5_; }

  /*! \brief Set extra double 5 setting */
  inline void setExtraDbl5(double extraDbl5) { extraDbl5_ = extraDbl5; }

  /*! \brief Get small branch and bound setting */
  inline double getSmallBaB() { return smallBaB_; }

  /*! \brief Set small branch and bound setting */
  inline void setSmallBab(double smallBaB) { smallBaB_ = smallBaB; }

  /*! \brief Get tighten factor */
  inline double getTightenFactor() { return tightenFactor_; }

  /*! \brief Set tighten factor

  */

  inline void setTightenFactor(double tightenFactor) {
    tightenFactor_ = tightenFactor;
  }

  //@}

  /*! \name Miscellaneous Keyword Parameters

  */

  //@{

  /*! \brief Get command mode */
  inline CbcSettings::CommandDisplayMode getCommandDisplayMode() {
     return (commandDisplayMode_);
  }

  /*! \brief Set command mode */
  inline void setCommandDisplayMode(int mode) {
     commandDisplayMode_ = static_cast<CbcSettings::CommandDisplayMode>(mode);
  }

  /*! \brief Get mode for clique strengthening */
  inline CbcSettings::ClqStrMode getClqStrMode() { return (clqStrMode_); }

  /*! \brief Set mode for use of integer preprocessing */
  inline void setClqStrMode(CbcSettings::ClqStrMode mode) { clqStrMode_ = mode;
  }

  /*! \brief Get mode for branching priorities */
  inline CbcSettings::BPMode getBranchPriority() { return (branchPriority_);
  }

  /*! \brief Set mode for branching priorities */
  inline void setBranchPriority(CbcSettings::BPMode mode) { branchPriority_ = mode;
  }

  /*! \brief Get mode for use of cutoff constraint */
  inline CbcSettings::CutoffMode getCutoffMode() { return (cutoffMode_); }

  /*! \brief Set mode for use of cutoff constraint */
  inline void setCutoffMode(CbcSettings::CutoffMode mode) { cutoffMode_ = mode;
  }

  /*! \brief Get mode for printing integers */
  inline CbcSettings::IntPrintMode getIntPrintMode() { return (intPrintMode_);
  }

  /*! \brief Set  mode for printing integers */
  inline void setIntPrintMode(CbcSettings::IntPrintMode mode) { intPrintMode_ = mode;
  }

  /*! \brief Get print mode */
  inline int getPrintMode() { return (printMode_); }

  /*! \brief Set print mode */
  inline void setIntPrintMode(int mode) { printMode_ = mode; }

  /*! \brief Get print mask */
   inline std::string getPrintMask() { return (printMask_);
  }

  /*! \brief Set print mask */
   inline void setPrintMask(std::string mask) { printMask_ = mask; }

  /*! \brief Get node search strategy */
  inline CbcSettings::NodeStrategy getNodeStrategy() { return (nodeStrategy_);
  }

  /*! \brief Set node search strategy */
  inline void setNodeStrategy(CbcSettings::NodeStrategy mode) { nodeStrategy_ = mode;
  }

  /*! \brief Get strategy for orbital branching */
  inline CbcSettings::OrbitalStrategy getOrbitalStrategy() { return (orbitalStrategy_);
  }

  /*! \brief Set  strategy for orbital branching */
  inline void setOrbitalStrategy(CbcSettings::OrbitalStrategy mode) {
     orbitalStrategy_ = mode;
  }

  /*! \brief Get mode for use of integer preprocessing */
  inline CbcSettings::IPPMode getIPPMode() { return (preProcess_); }

  /*! \brief Set mode for use of integer preprocessing */
  inline void setIPPMode(CbcSettings::IPPMode mode) { preProcess_ = mode; }

  /*! \brief Get priority mode for SOS */
  inline CbcSettings::SOSStrategy getSOSStrategy() { return (sosStrategy_); }

  /*! \brief Set mode state for use of integer preprocessing */
  inline void setSOSStrategy(CbcSettings::SOSStrategy mode) {
    sosStrategy_ = mode;
  }

  /*! \brief Get overall strategy mode */
  inline CbcSettings::StrategyMode getStrategyMode() { return (strategyMode_); }

  /*! \brief Set overall strategy mode */
  inline void setStrategyMode(CbcSettings::StrategyMode mode) { strategyMode_ = mode;
  }

  /*! \brief Get clock type */
  inline CbcSettings::ClockType getClockType() { return (clockType_); }

  /*! \brief Set clock type */
  inline void setClockType(CbcSettings::ClockType type) { clockType_ = type;
  }

  /*! \brief Get mode for CGraph */
  inline CbcSettings::CGraphMode getCGraphMode() { return (cgraphMode_); }

  /*! \brief Set mode for CGraph */
  inline void setCGraphMode(CbcSettings::CGraphMode mode) { cgraphMode_ = mode;
  }

  /*! \brief Get threshold for artificial costs */
  inline double getArtVarThreshold() { return (artVar_.threshold_); }

  /*! \brief Get mode for artificial costs */
  inline double getArtVarMode() { return (artVar_.mode_); }

  /*! \brief Set mode for artificial costs */
  inline void setArtVarMode(CbcSettings::OnOffMode mode, double threshold) {
    artVar_.threshold_ = threshold;
    artVar_.mode_ = mode;
  }

  /*! \brief Get threshold for reduced costs fixing */
  inline double getDjFixThreshold() { return (djFix_.threshold_); }

  /*! \brief Get mode for reduced cost fixing */
  inline double getDjFixMode() { return (djFix_.mode_); }

  /*! \brief Set mode for reduced cost fixing */
  inline void setDjFixMode(CbcSettings::OnOffMode mode, double threshold) {
    djFix_.threshold_ = threshold;
    djFix_.mode_ = mode;
  }

  /*! \brief Get mode for branching priorities */
  inline CbcSettings::BPMode getBranchingPriorityMode() {
    return (priorityMode_);
  }

  /*! \brief Set mode for reduced cost fixing */
  inline void setBranchingPriorityMode(CbcSettings::BPMode mode) {
    priorityMode_ = mode;
  }

  //@}

  /*! \name Miscellaneous Bool Parameters

  */

  //@{

  /*! \brief Get CPX mode */
  inline CbcSettings::OnOffMode getCpxMode() { return (CPXMode_); }

  /*! \brief Set CPX mode */
  inline void setCPXMode(CbcSettings::OnOffMode mode) { CPXMode_ = mode; }

  /*! \brief Get import errors mode */
  inline CbcSettings::OnOffMode getImportErrorsMode() {
    return (importErrorsMode_);
  }

  /*! \brief Set import errors mode */
  inline void setImportErrorsMode(CbcSettings::OnOffMode mode) {
    importErrorsMode_ = mode;
  }

  /*! \brief Get message prefix mode */
  inline CbcSettings::OnOffMode getMessagePrefixMode() {
    return (messagePrefixMode_);
  }

  /*! \brief Set message prefix mode */
  inline void setMessagePrefixMode(CbcSettings::OnOffMode mode) {
    messagePrefixMode_ = mode;
  }

  /*! \brief Get preprocess names mode */
  inline CbcSettings::OnOffMode getPreProcNamesMode() {
    return (preProcNamesMode_);
  }

  /*! \brief Set preprocess names mode */
  inline void setPreProcNamesMode(CbcSettings::OnOffMode mode) {
    preProcNamesMode_ = mode;
  }

  /*! \brief Get SOS mode */
  inline CbcSettings::OnOffMode getSOSMode() { return (SOSMode_); }

  /*! \brief Set SOS mode */
  inline void setSOSMode(CbcSettings::OnOffMode mode) { SOSMode_ = mode; }

  /*! \brief Get use solution mode */
  inline CbcSettings::OnOffMode getUseSolutionMode() {
    return (useSolutionMode_);
  }

  /*! \brief Set use solution mode */
  inline void setUseSolutionMode(CbcSettings::OnOffMode mode) {
    useSolutionMode_ = mode;
  }

  //@}

  /*! \name Status Functions
        \brief Convenience routines for status codes.
    */
  //@{

  /*! \brief Set the result of branch-and-cut search */
  inline void setBaBStatus(CbcSettings::BACMajorStatus majorStatus,
                           CbcSettings::BACMinorStatus minorStatus,
                           CbcSettings::BACWhere where, bool haveAnswer,
                           OsiSolverInterface *answerSolver) {
    bab_.majorStatus_ = majorStatus;
    bab_.minorStatus_ = minorStatus;
    bab_.where_ = where;
    bab_.haveAnswer_ = haveAnswer;
    bab_.answerSolver_ = answerSolver;
  }

  /*! \brief Set the result of branch-and-cut search

      This version will extract the necessary information from the CbcModel
      object and set appropriate status based on the value passed for where.
    */
  void setBaBStatus(CbcSettings::BACWhere where,
                    bool haveAnswer = false,
                    OsiSolverInterface *answerSolver = 0);

  /*! \brief Get status of branch-and-cut search */
  inline bool haveAnswer() { return  bab_.haveAnswer_; }
  inline OsiSolverInterface *answerSolver() { return bab_.answerSolver_; }

  /*! \brief Translate CbcModel major status to #BACMajorStatus

      See the #BACMajorStatus enum for details.
    */
  CbcSettings::BACMajorStatus translateMajor(int status);

  /*!\brief Translate CbcModel minor status to #BACMinorStatus

      See the #BACMinorStatus enum for details.
    */
  CbcSettings::BACMinorStatus translateMinor(int status);

  /*!\brief Translate OsiSolverInterface status to #BACMinorStatus

      See the #BACMinorStatus enum for details. Optimal, infeasible, and
      unbounded get their own codes; everything else maps to BACmOther.
    */
  CbcSettings::BACMinorStatus translateMinor(const OsiSolverInterface *osi);

  /*! \brief Print the status block */

  void printBaBStatus();

  //@}

  /*! \name Messages and statistics */
  //@{

  /*! \brief Print a message
      Uses the current message handler and messages.
    */
  CoinMessageHandler &message(CbcGenMsgCode inID);

  /*! \brief Supply a new message handler.

      Replaces the current message handler. The current handler is destroyed
      if ourMsgHandler_ is true, and the call will set ourMsgHandler_ = true.
    */
  void passInMessageHandler(CoinMessageHandler *handler);

  /*! \brief Return a pointer to the message handler */
  inline CoinMessageHandler *messageHandler() const { return msgHandler_; }

  /*! \brief Set up messages in the specified language.

      Building a set of messages in a given language implies rebuilding the
      whole set of messages, for reasons explained in the body of the code.
      Hence there's no separate setLanguage routine. Use this routine for the
      initial setup of messages and any subsequent change in language. Note
      that the constructor gives you a message handler by default, but \e not
      messages. You need to call setMessages explicitly.

      The default value specified here for lang effectively sets the default
      language.
    */
  void setMessages(CoinMessages::Language lang = CoinMessages::us_en);

  /*! \brief Set log level */
  inline void setLogLevel(int lvl) {
    logLvl_ = lvl;
    if (msgHandler_)
      msgHandler_->setLogLevel(lvl);
  }

  /*! \brief Get log level */
  inline int getLogLevel() const { return (logLvl_); }

  /*! \brief Set LP log level */
  inline void setLpLogLevel(int lvl) { lpLogLvl_ = lvl; }

  /*! \brief Get LP log level */
  inline int getLpLogLevel() const { return (lpLogLvl_); }

  /*! \brief Set LP log level */
  inline void setModel(CbcModel *model) { model_ = model; }

  /*! \brief Get LP log level */
  inline CbcModel *getModel() const { return (model_); }

  /*! \brief Get total time */
  inline double getTotalTime() const { return totalTime_; }

  /*! \brief Set total time */
  inline void setTotalTime(double t) { totalTime_ = t; }

  /* \brief Get Clp solver parameter vector */
  inline std::vector < ClpParam > &getClpParameters() { return clpParameters_; }
   
  /* \brief Get Cbc solver parameter vector */
  inline CoinParamVec &getParameters() { return parameters_; }

  /* \brief Get version */
  inline std::string getVersion() { return version_; }

  /* \brief Get default directory */
  inline std::string getDefaultDirectory() { return dfltDirectory_; }

  /* \brief Set default directory */
   inline void setDefaultDirectory(std::string dir) { dfltDirectory_ = dir; }

  /* \brief Get model status */
  inline bool goodModel() { return goodModel_; }

  /* \brief Set model status */
   inline void setGoodModel(bool gm) { goodModel_ = gm; }

  /* \brief Get debug file */
  inline std::string getDebugFile() { return debugFile_; }

  /* \brief Set debug file */
  inline void setDebugFile(std::string f) { debugFile_ = f; }

  /* \brief Get debug creation file */
  inline std::string getDebugCreate() { return debugCreate_; }

  /* \brief Set debug creation file */
  inline void setDebugCreate(std::string d) { debugCreate_ = d; }

  /* \brief Get last MPS file */
  inline std::string getLastMpsIn() { return lastMpsIn_; }

  /* \brief Set last MPS file */
  inline void setLastMpsIn(std::string f) { lastMpsIn_ = f; }

  /* \brief Get last solution file */
  inline std::string getLastSolnOut() { return lastSolnOut_; }

  /* \brief Set last solution file */
  inline void setLastSolnOut(std::string f) { lastSolnOut_ = f; }

  /* \brief Toggle set by user for given parameter */
  inline void toggleSetByUser(CbcParamCode code) { setByUser_[code] = true; }

  /* \brief Toggle set by user for given parameter */
  inline bool setByUser(CbcParamCode code) { return setByUser_[code]; }

  /* \brief Set debug solution */
  void setDebugSol(int numCols, double* values) {
     if (debugSol_.values_){
        delete[] debugSol_.values_;
     }
     debugSol_.numCols_ = numCols;
     debugSol_.values_ = values;
  }

  /* \brief Get debug solution */
  inline int getDebugSolNumCols() const { return debugSol_.numCols_; }
  inline double *getDebugSolValues() const { return debugSol_.values_; }
   
  /* \brief getShadowPriceMode */
  inline int getShadowPriceMode() { return chooseStrong_.shadowPriceMode_; }

  //@}

private:

  friend class CbcParam; 
   
  /*! \name Parameter parsing and input/output. */
  //@{
  /*! \brief cbc-generic version */
  std::string version_;

  /*! \brief Default directory prefix */
  std::string dfltDirectory_;

  /*! \brief Last MPS input file */
  std::string lastMpsIn_;

  /*! \brief Allow/disallow errors when importing a model */
  bool allowImportErrors_;

  /*! \brief The Cbc parameter vector (parameters stored by their index) */
  CoinParamVec parameters_;

  /*! \brief The Clp parameter vector */
  std::vector < ClpParam > clpParameters_;

  /*! \brief Last solution output file */
  std::string lastSolnOut_;

  /*! \brief Solution printing mode
      Controls the amount of information printed when printing a solution.
      Coding is set by the keyword declarations for the printingOptions
      command.
    */
  int printMode_;

  /*! \brief When greater than 0, integer presolve gives more information and
           branch-and-cut provides statistics.
    */
  int printOpt_;

  /*! \brief Print mask

      Used to specify row/column names to be printed. Not implemented as of
      060920.
    */
  std::string printMask_;

  /*! \brief Disable printing altogether */
  bool noPrinting_;

  /*! \brief Disable printing of welcome message */
  bool printWelcome_;

  /*! \brief Whether to use ginal handler */
  bool useSignalHandler_;
   
  /*! \brief Verbosity level for help messages.

      Interpretation is bitwise:
      - (0): short help
      - (1): long help
      - (2): unused (for compatibility with cbc; indicates AMPL)
      - (3): show parameters with display = false.
    */
  int verbose_;

  /*! \brief Number of parameters processed */
  int paramsProcessed_;

  /*! \brief Record of parameters changed by user command */
  std::vector<bool> setByUser_;

  /*! \brief False if the user has made nontrivial modifications to the
           default control settings.

      Initially true. Specifying DJFIX, TIGHTENFACTOR, or any cut or heuristic
      parameter will set this to false.
    */
  bool defaultSettings_;

  /*! \brief Control debug file creation

      At the conclusion of branch-and-cut, dump the full solution in a binary
      format to debug.file in the current directory.  When set to
      "createAfterPre", the solution is dumped before integer presolve
      transforms are removed.  When set to "create", the solution is dumped
      after integer presolve transforms are backed out.
    */
  std::string debugCreate_;

  /*! \brief Last debug input file

      The file is expected to be in a binary format understood by
      activateRowCutDebugger.
    */
  std::string debugFile_;

  /*! \brief Array of primal variable values for debugging

      Used to provide a known optimal solution to activateRowCutDebugger().
    */
  struct debugSolInfo_struct {
    int numCols_;
    double *values_;
  } debugSol_;
  //@}

  /* \name Timing */
  //@{

  /*! \brief Total elapsed time for this run. */
  double totalTime_;

  //@}

  /*! \name Models of various flavours */
  //@{

  /*! \brief The reference CbcModel object.

      This is the CbcModel created when cbc-generic boots up. It holds the
      default solver with the current constraint system. CbcCbcParam parameters
      are applied here, and CbcOsiParam parameters are applied to the solver.
      Major modifications for branch-and-cut (integer preprocessing,
      installation of heuristics and cut generators) are performed on a clone.
      The solution is transferred back into this object.
    */
  CbcModel *model_;

  /*! \brief The current default LP solver

      This is a pointer to a reference copy. If you want the solver associated
      with #model_, ask for it directly.
    */
  OsiSolverInterface *dfltSolver_;

  /*! \brief True if we have a valid model loaded, false otherwise. */
  bool goodModel_;

  /*! \brief State of branch-and-cut

      Major and minor status codes, and a solver holding the answer, assuming
      we have a valid answer. See the documentation with the BACMajorStatus,
      BACMinorStatus, and BACWhere enums for the meaning of the codes.
    */
  struct babState_struct {
    CbcSettings::BACMajorStatus majorStatus_;
    CbcSettings::BACMinorStatus minorStatus_;
    CbcSettings::BACWhere where_;
    bool haveAnswer_;
    OsiSolverInterface *answerSolver_;
  } bab_;

  //@}

  /*! \name Various algorithm control variables and settings */
  //@{

  /*! \brief Control use of reduced cost fixing prior to B&C

      This heuristic fixes variables whose reduced cost for the root
      relaxtion exceeds the specified threshold. This is purely a heuristic,
      performed before there's any incumbent solution. It may well fix variables
      at the wrong bound!
    */
  struct djFix_struct {
    CbcSettings::OnOffMode mode_;
    double threshold_;
  } djFix_;

  /*! \brief Treat some variables as artificial in feasibility pump

    */
  struct artVar_struct {
    CbcSettings::OnOffMode mode_;
    double threshold_;
  } artVar_;

  /*! \brief Control the assignment of branching priorities to integer
           variables.
    */
  CbcSettings::BPMode priorityMode_;

  //@}

  /*! \name Branching Method Control
        \brief Usage control and prototypes for branching methods.

        Looking to the future, this covers only OsiChoose methods.
    */
  //@{

  /*! \brief Control variables for a strong branching method.

      Consult OsiChooseVariable and CbcModel for details. An artifact of the
      changeover from CbcObjects to OsiObjects is that the number of uses before
      pseudo costs are trusted (numBeforeTrust_) and the number of variables
      evaluated with strong branching (numStrong_) are parameters of CbcModel.
    */
  struct chooseStrong_struct {
    int numBeforeTrust_;
    int numStrong_;
    int shadowPriceMode_;
  } chooseStrong_;
  //@}

  /*! \name Cut Generator Control
        \brief Usage control and prototypes for cut generators.
    */
  //@{

  /*! \brief Control integer preprocessing. */
  CbcSettings::IPPMode preProcess_;

  /*! \brief Control cut generator activity

      Generators that are active in the tree will be activated when
      (depth) mod (cutDepth) == 0.
    */
  int cutDepth_;
  int cutLength_;
  int cutPassInTree_;

  /*! \brief Control variable and prototype for clique cut generator */
  struct clique_struct {
    CbcSettings::CGMode mode_;
    CglClique *proto_;
    bool starCliqueReport_;
    bool rowCliqueReport_;
    double minViolation_;
  } clique_;

  /*! \brief Control variable and prototype for flow cover cut generator */
  struct flow_struct {
    CbcSettings::CGMode mode_;
    CglFlowCover *proto_;
  } flow_;

  /*! \brief Control variable and prototype for Gomory cut generator */
  struct gmi_struct {
    CbcSettings::CGMode mode_;
    CglGomory *proto_;
    int limit_;
    int limitAtRoot_;
  } gmi_;

  /*! \brief Control variable and prototype for Gomory cut generator */
  struct gomory_struct {
    CbcSettings::CGMode mode_;
    CglGomory *proto_;
    int limit_;
    int limitAtRoot_;
  } gomory_;

  /*! \brief Control variable and prototype for knapsack cover cut generator */
  struct knapsack_struct {
    CbcSettings::CGMode mode_;
    CglKnapsackCover *proto_;
  } knapsack_;

  /*   \brief Control variable and prototype for LaGomory cut generator */
  struct lagomory_struct {
    CbcSettings::CGMode mode_;
    CglGomory *proto_;
  } laGomory_;

  /*   \brief Control variable and prototype for lift-and-project cut generator */
  struct landp_struct {
    CbcSettings::CGMode mode_;
    CglLandP *proto_;
  } landP_;

  /*   \brief Control variable and prototype for lift-and-project cut generator */
  struct laTwomir_struct {
    CbcSettings::CGMode mode_;
    CglTwomir *proto_;
  } laTwomir_;

  /*! \brief Control variable and prototype for MIR cut generator */
  struct mir_struct {
    CbcSettings::CGMode mode_;
    CglMixedIntegerRounding2 *proto_;
  } mir_;

  /*! \brief Control variable and prototype for odd hole cut generator */
  struct oddWheel_struct {
    CbcSettings::CGMode mode_;
    CglOddWheel *proto_;
  } oddWheel_;

  /*! \brief Control variable and prototype for probing cut generator */
  struct probing_struct {
    CbcSettings::CGMode mode_;
    CglProbing *proto_;
    bool usingObjective_;
    int maxPass_;
    int maxPassRoot_;
    int maxProbe_;
    int maxProbeRoot_;
    int maxLook_;
    int maxLookRoot_;
    int maxElements_;
    int rowCuts_;
  } probing_;

  /*! \brief Control variable and prototype for reduce-and-split cut generator */
  struct redSplit_struct {
    CbcSettings::CGMode mode_;
    CglRedSplit *proto_;
  } redSplit_;

  /*! \brief Control variable and prototype for reduce-and-split 2 cut generator */
  struct redSplit2_struct {
    CbcSettings::CGMode mode_;
    CglRedSplit2 *proto_;
  } redSplit2_;

  /*! \brief Control variable and prototype for residual capacity cut generator */
  struct residCap_struct {
    CbcSettings::CGMode mode_;
    CglResidualCapacity *proto_;
  } residCap_;

  /*! \brief Control variable and prototype for Two-MIR cut generator */
  struct twomir_struct {
    CbcSettings::CGMode mode_;
    CglTwomir *proto_;
    int maxElements_;
  } twomir_;

  /*! \brief Control variable and prototype for residual capacity cut generator */
  struct zeroHalf_struct {
    CbcSettings::CGMode mode_;
    CglZeroHalf *proto_;
  } zeroHalf_;

  //@}

  /*! \name Heuristic Control
        \brief Usage control and prototypes for heuristics.
  */
  //@{

  /*! \brief Overall control variable for heuristics */
  CbcSettings::OnOffMode doHeuristicMode_;

  /*! \brief Control variable and prototype for combine heuristic */
  struct combine_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicLocal *proto_;
    int trySwap_;
  } combine_;

  /*! \brief Control variable and prototype for crossover heuristic */
  struct crossover_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicCrossover *proto_;
  } crossover_;

  /*! \brief Control variable and prototype for heuristic */
  struct dins_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicDINS *proto_;
  } dins_;

  /*! \brief Control variable and prototype for Dive Coefficient heuristic */
  struct diveCoefficient_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicDiveCoefficient *proto_;
  } divingc_;

  /*! \brief Control variable and prototype for Dive Fractional heuristic */
  struct diveFractional_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicDiveFractional *proto_;
  } divingf_;

  /*! \brief Control variable and prototype for Dive Guided heuristic */
  struct diveGuided_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicDiveGuided *proto_;
  } divingg_;

  /*! \brief Control variable and prototype for Dive Line Sarch heuristic */
  struct diveLineSearch_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicDiveLineSearch *proto_;
  } divingl_;

  /*! \brief Control variable and prototype for Dive Pseudocost heuristic */
  struct _struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicDivePseudoCost *proto_;
  } divingp_;

  /*! \brief Control variable and prototype for Dive Vector Lengthheuristic */
  struct diveVectorLength_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicDiveVectorLength *proto_;
  } divingv_;

  /*! \brief Control variable and prototype for DW heuristic */
  struct DW_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicDW *proto_;
  } dw_;

  /*! \brief Control variable and prototype for feasibility pump heuristic */
  struct fpump_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicFPump *proto_;
    int iters_;
    int initialTune_; 
    int tune_;
    int tune2_;
    double cutoff_;
    double increment_;
  } fpump_;

  /*! \brief Control variable and prototype for greedy cover heuristic */
  struct greedyCover_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicGreedyCover *proto_;
  } greedyCover_;

  /*! \brief Control variable and prototype for greedy equality heuristic */
  struct greedyEquality_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicGreedyEquality *proto_;
  } greedyEquality_;

  /*! \brief Control variable and prototype for Naive heuristic */
  struct naive_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicNaive *proto_;
  } naive_;

  /*! \brief Control variable and prototype for Pivot and Fix heuristic */
  struct pivotAndFix_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicPivotAndFix *proto_;
  } pivotAndFix_;

#if 0
  /*! \brief Control variable and prototype for Pivot and Complement heuristic */
  struct pivotAndComp_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicPivotAndComplement *proto_;
  } pivotAndComplement_;
#endif

  /*! \brief Control variable and prototype for Proximity heuristic */
  struct proximity_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicProximity *proto_;
  } proximity_;

  /*! \brief Control variable and prototype for Randomized Rounding heuristic */
  struct randRound_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristic *proto_;
  } randRound_;

  /*! \brief Control variable and prototype for RENS heuristic */
  struct Rens_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicRENS *proto_;
  } rens_;

  /*! \brief Control variable and prototype for RINS heuristic */
  struct Rins_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicRINS *proto_;
  } rins_;

  /*! \brief Control variable and prototype for simple rounding heuristic */
  struct rounding_struct {
    CbcSettings::HeurMode mode_;
    CbcRounding *proto_;
  } rounding_;

  /*! \brief Control variable and prototype for Variable Neighborhood
    heuristic
  */
  struct Vnd_struct {
    CbcSettings::HeurMode mode_;
    CbcHeuristicVND *proto_;
  } vnd_;

  CbcSettings::HeurMode randomDivingMode_;

  /*! \brief Control variables for local tree

      This is a bit different --- getTreeLocal() takes a CbcModel as a parameter
      and installs a local tree object. But we can keep the parameters here and
      hide the details. Consult CbcTreeLocal.hpp for details.
    */
  struct localTree_struct {
    CbcSettings::HeurMode mode_;
    CbcTreeLocal *proto_;
    double *soln_;
    int range_;
    int typeCuts_;
    int maxDiverge_;
    int timeLimit_;
    int nodeLimit_;
    bool refine_;
  } localTree_;

  //@}

  /*! \name Miscellaneous double paramters

  */
  //@{

  double extraDbl3_;
  double extraDbl4_;
  double extraDbl5_;
  double smallBaB_;
  double tightenFactor_;

  //@}

  /*! \name Miscellaneous integer paramters

  */
  //@{

  int bkPivotStrategy_;
  int bkMaxCalls_;
  int bkClqExtMethod_;
  int cppMode_;
  int depthMiniBaB_;
  int diveOpt_;
  int diveOptSolves_;
  int experiment_;
  int extraInt1_;
  int extraInt2_;
  int extraInt3_;
  int extraInt4_;
  int heurOptions_;
  int maxSavedSols_;
  int maxSlowCuts_;
  int moreMoreOptions_;
  int multipleRoots_;
  int oddWextMethod_;
  int outputFormat_;
  int processTune_;
  int randomSeed_;
  int strongStrategy_;
  int testOsi_;
  int threads_;
  int userCbc_;
  int vubTry_;

  //@}

  /*! \name Miscellaneous keyword paramters

  */
  //@{

  CbcSettings::CommandDisplayMode commandDisplayMode_;
  CbcSettings::ClqStrMode clqStrMode_;
  CbcSettings::BPMode branchPriority_;
  CbcSettings::CutoffMode cutoffMode_;
  CbcSettings::IntPrintMode intPrintMode_;
  CbcSettings::NodeStrategy nodeStrategy_;
  CbcSettings::OrbitalStrategy orbitalStrategy_;
  CbcSettings::SOSStrategy sosStrategy_;
  CbcSettings::StrategyMode strategyMode_;
  CbcSettings::ClockType clockType_;
  CbcSettings::CGraphMode cgraphMode_;

  //@}

  /*! \name Miscellaneous bool paramters

  */
  //@{

  CbcSettings::OnOffMode CPXMode_;
  CbcSettings::OnOffMode importErrorsMode_;
  CbcSettings::OnOffMode messagePrefixMode_;
  CbcSettings::OnOffMode preProcNamesMode_;
  CbcSettings::OnOffMode SOSMode_;
  CbcSettings::OnOffMode useSolutionMode_;

  //@}

  /*! \name Messages and statistics (private)
        \brief Data and objects related to messages and statistics that should
     be protected from direct manipulation.
    */
  //@{

  /*! \brief Message handler. */
  CoinMessageHandler *msgHandler_;

  /*! \brief Ownership of message handler.

      If true, the control block owns the message handler and it will be
     destroyed with the control block. If false, the client is responsible for
     the message handler.
    */
  bool ourMsgHandler_;

  /*! \brief The current language */
  CoinMessages::Language cur_lang_;

  /*! \brief The current set of messages. */
  CoinMessages *msgs_;

  /*! \brief The current log level */
  int logLvl_;

  /*! \brief The current LP log level */
  int lpLogLvl_;

  //@}
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
