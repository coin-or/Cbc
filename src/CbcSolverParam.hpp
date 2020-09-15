/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/
/*
  This file is part of cbc-generic.
*/

#ifndef CbcSolverParam_H
#define CbcSolverParam_H

#include "CoinParam.hpp"

#include "OsiSolverInterface.hpp"

/* \file CbcSolverParam.hpp
   \brief Declarations for parameters that control the cbc-generic main
          program.
*/

/*
 */

class CbcSolverSettings;

/*! \class CbcSolverParam
    \brief Class for cbc-generic control parameters

    Adds parameter type codes and push/pull functions to the generic parameter
    object.
*/

class CbcSolverParam : public CoinParam {

public:
  /*! \name Subtypes */
  //@{

  enum CbcSolverParamCode {
    CBCSOLVER_FIRSTPARAM = 0,

    // Unused paramters that we may delete
    CBCSOLVER_FIRSTUNUSEDPARAM,
    BRANCHSTRATEGY,
    CLEARCUTS,
    SOLVER,
    CBCSOLVER_LASTUNUSEDPARAM,

    // Help and Information Parameters
    CBCSOLVER_FIRSTHELPPARAM,
    GENERALQUERY,
    FULLGENERALQUERY,
    HELP,
    PRINTVERSION,
    CBCSOLVER_LASTHELPPARAM,

    // Action Parameters
    CBCSOLVER_FIRSTACTIONPARAM,
    BAB,
    DUMMY,
    ENVIRONMENT,
    EXIT,
    EXPORT,
    GMPL_SOLUTION,
    IMPORT,
    MIPLIB,
    OUTDUPROWS,
    SHOWUNIMP,
    SOLVECONTINUOUS,
    STATISTICS,
    STDIN,
    STRENGTHEN,
    UNITTEST,
    CBCSOLVER_LASTACTIONPARAM,

    // String (Directory) Parameters
    CBCSOLVER_FIRSTSTRPARAM,
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
    CBCSOLVER_LASTSTRPARAM,

    // Cut Parameters
    CBCSOLVER_FIRSTCUTPARAM,
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
    CBCSOLVER_LASTCUTPARAM,

    // Heuristic Parameters
    CBCSOLVER_FIRSTHEURPARAM,
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
    CBCSOLVER_LASTHEURPARAM,

    // On/Off Parameters
    CBCSOLVER_FIRSTBOOLPARAM,
    CPX,
    DOHEURISTIC,
    ERRORSALLOWED,
    EXTRAVARIABLES,
    MESSAGES,
    PREPROCNAMES,
    SOS,
    USESOLUTION,
    CBCSOLVER_LASTBOOLPARAM,

    // Keyword Parameters
    CBCSOLVER_FIRSTKWDPARAM,
    ALLCOMMANDS,
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
    CBCSOLVER_LASTKWDPARAM,

    // Integer Parameters
    CBCSOLVER_FIRSTINTPARAM,
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
    CBCSOLVER_LASTINTPARAM,

    // Double Parameters
    CBCSOLVER_FIRSTDBLPARAM,
    ARTIFICIALCOST,
    DEXTRA3,
    DEXTRA4,
    DEXTRA5,
    DJFIX,
    FAKECUTOFF,
    FAKEINCREMENT,
    SMALLBAB,
    TIGHTENFACTOR,
    CBCSOLVER_LASTDBLPARAM,

    CBCSOLVER_LASTPARAM
  };

  //@}

  /*! \name Enumeration types used for CBC keyword parameters */
  //@{

  /*
      In order for initialisation to work properly, the order of declaration of
      the enum constants here must match the order of keyword declaration for
      the PREPROCESS parameter in CbcGenParamUtils::addCbcGenParams
  */
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

    - CMOff:
    - CMMore:
    - CMAll:
    - CMEndMarker

   */

  enum CommandMode { CMOff = 0, CMMore, CMAll, CMEndMarker };

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

  //@}

  /*! \name Constructors and Destructors

      Be careful how you specify parameters for the constructors! There's great
      potential for confusion.
    */
  //@{
  /*! \brief Default constructor */

  CbcSolverParam();

  /*! \brief Constructor for a parameter with a double value

      The default value is 0.0. Be careful to clearly indicate that \p lower and
      \p upper are real (double) values to distinguish this constructor from the
      constructor for an integer parameter.
    */
  CbcSolverParam(CbcSolverParamCode code, std::string name, std::string help,
                 double lower = -COIN_DBL_MAX, double upper = COIN_DBL_MAX,
                 double defaultValue = 0.0, std::string longHelp = "",
                 CoinDisplayPriority display = displayPriorityHigh);

  /*! \brief Constructor for a parameter with an integer value

      The default value is 0.
    */
  CbcSolverParam(CbcSolverParamCode code, std::string name, std::string help,
                 int lower = -COIN_INT_MAX, int upper = COIN_INT_MAX,
                 int defaultValue = 0, std::string longHelp = "",
                 CoinDisplayPriority display = displayPriorityHigh);

  /*! \brief Constructor for a parameter with keyword values

      The string supplied as \p firstValue becomes the first keyword.
      Additional keywords can be added using appendKwd(). Keywords are numbered
      from zero. It's necessary to specify both the first keyword (\p
      firstValue) and the default keyword index (\p dflt) in order to
      distinguish this constructor from the string and action parameter
      constructors.
    */
  CbcSolverParam(CbcSolverParamCode code, std::string name, std::string help,
                 std::string defaultKwd, int defaultMode,
                 std::string longHelp = "",
                 CoinDisplayPriority display = displayPriorityHigh);

  /*! \brief Constructor for a string parameter

      The default string value must be specified explicitly to distinguish
      a string constructor from an action parameter constructor.
    */

  CbcSolverParam(CbcSolverParamCode code, std::string name, std::string help,
                 std::string defaultValue, std::string longHelp = "",
                 CoinDisplayPriority display = displayPriorityHigh);

  /*! \brief Constructor for an action parameter */

  // No defaults to resolve ambiguity
  CbcSolverParam(CbcSolverParamCode code, std::string name, std::string help,
                 std::string longHelp, CoinDisplayPriority display);

  /*! \brief Copy constructor */

  CbcSolverParam(const CbcSolverParam &orig);

  /*! \brief Clone */

  CbcSolverParam *clone();

  /*! \brief Assignment */

  CbcSolverParam &operator=(const CbcSolverParam &rhs);

  /*! \brief  Destructor */

  ~CbcSolverParam();

  //@}

  /*! \name Methods to query and manipulate a parameter object */
  //@{

  /*! \brief Get the parameter code  */

  inline CbcSolverParamCode paramCode() const { return (paramCode_); }

  /*! \brief Set the parameter code */

  inline void setParamCode(CbcSolverParamCode code) { paramCode_ = code; }

  /*! \brief Get the underlying cbc-generic control object */

  inline CbcSolverSettings *obj() const { return (obj_); }

  /*! \brief Set the underlying cbc-generic control object */

  inline void setObj(CbcSolverSettings *obj) { obj_ = obj; }

  //@}

private:
  /*! \name Data */
  //@{

  /// Parameter code
  CbcSolverParamCode paramCode_;

  /// cbc-generic control object
  CbcSolverSettings *obj_;

  //@}
};

/*
  Declare the utility functions.
*/

namespace CbcSolverParamUtils {

void addCbcSolverParams(CoinParamVec &paramVec, CbcSolverSettings *cbcSettings);
void addCbcSolverStrParams(CoinParamVec &paramVec,
                           CbcSolverSettings *cbcSettings);
void addCbcSolverHelpParams(CoinParamVec &paramVec,
                            CbcSolverSettings *cbcSettings);
void addCbcSolverActionParams(CoinParamVec &paramVec,
                              CbcSolverSettings *cbcSettings);
void addCbcSolverKwdParams(CoinParamVec &paramVec,
                           CbcSolverSettings *cbcSettings);
void addCbcSolverDblParams(CoinParamVec &paramVec,
                           CbcSolverSettings *cbcSettings);
void addCbcSolverIntParams(CoinParamVec &paramVec,
                           CbcSolverSettings *cbcSettings);
void addCbcSolverBoolParams(CoinParamVec &paramVec,
                            CbcSolverSettings *cbcSettings);
void addCbcSolverCutParams(CoinParamVec &paramVec,
                           CbcSolverSettings *cbcSettings);
void addCbcSolverHeurParams(CoinParamVec &paramVec,
                            CbcSolverSettings *cbcSettings);

void loadGenParamObj(const CoinParamVec paramVec, int first, int last,
                     CbcSolverSettings *cbcSetting);

void saveSolution(const OsiSolverInterface *osi, std::string fileName);
bool readSolution(std::string fileName, int &numRows, int &numCols,
                  double &objVal, double **rowActivity, double **dualVars,
                  double **primalVars, double **reducedCosts);

int doBaCParam(CoinParam &param);
int doDebugParam(CoinParam &param);
int doExitParam(CoinParam &param);
int doHelpParam(CoinParam &param);
int doImportParam(CoinParam &param);
int doPrintMaskParam(CoinParam &param);
int doNothingParam(CoinParam &param);
int doSolutionParam(CoinParam &param);
int doUnimplementedParam(CoinParam &param);
int doVersionParam(CoinParam &param);

int pushCbcSolverDblParam(CoinParam &param);
int pushCbcSolverIntParam(CoinParam &param);
int pushCbcSolverKwdParam(CoinParam &param);
int pushCbcSolverStrParam(CoinParam &param);
int pushCbcSolverBoolParam(CoinParam &param);

int pushCbcSolverHeurParam(CoinParam &param);
int pushCbcSolverCutParam(CoinParam &param);
} // namespace CbcSolverParamUtils

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
