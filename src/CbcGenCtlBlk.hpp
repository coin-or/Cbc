/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/
/*
  This file is part of cbc-generic.
*/

#ifndef CbcSettings_H
#define CbcSettings_H

/* \file CbcSettings.hpp
   \brief Declarations for parameters of the cbc-generic main program.
*/

#include "CoinParam.hpp"
#include "CoinMessageHandler.hpp"

#include "CglCutGenerator.hpp"
#include "CglProbing.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglGomory.hpp"
#include "CglKnapsackCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglOddHole.hpp"
#include "CglRedSplit.hpp"
#include "CglTwomir.hpp"

#include "CbcModel.hpp"

#include "CbcHeuristic.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcTreeLocal.hpp"

#include "CbcGenMessages.hpp"

/*
  It turns out that doxygen is not good with anonymous structures. Hence the
  `struct name_struct' style used for structured fields in CbcSettings.
*/

/*
*/

#define CBC_GENERIC_VERSION "00.01.00"

class CbcSettings;
namespace CbcGenParamUtils {
void addCbcGenParams(int &numParams, CoinParamVec &paramVec,
  CbcSettings *settings);
}

/* \brief CBC algorithm control class

  This class defines values and methods used to control the operation of the
  CBC main program.
*/

class CbcSettings {

  friend void CbcGenParamUtils::addCbcGenParams(int &numParams,
    CoinParamVec &paramVec, CbcSettings *settings);

public:
  /*! \name Enumeration types used for CBC control variables */
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

  typedef enum {
     IPPOff = 0,
     IPPon,
     IPPsave,
     IPPequal,
     IPPsos,
     IPPtrySOS,
     IPPequalAll,
     IPPstrategy,
     IPPEdnMarker
  } IPPMode;

  /*! \brief What parameters to print

    - CMOff:
    - CMMore:
    - CMAll:
    - CMEndMarker

   */

   tyepdef enum{
      CMOff = 0,
      CMMore,
      CMAll,
      CMEndMarker
   } CommandMode;

  /*! \brief When to do clique strengthening

    - ClqOff:
    - ClqAfter:
    - ClqBefore:
    - ClqEndMarker:

  */
   
  typedef enum{
     ClqStrOff = 0,
     ClqStrAfter,
     ClqStrBefore,
     ClqStrEndMarker
  } ClqStrMode;
     
  /*! \brief What node strategy to use

    - NSHybrid:
    - NSFewest:
    - NSDepth:
    - NSUpFewest:
    - NSUpDepth:
    - NSDownDepth

   */

  typedef enum{
     NSHybrid = 0,
     NSFewest,
     NSDepth,
     NSUpFewest,
     NSUpDepth,
     NSDownDepth
  } NodeStrategy;

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
   
  typedef enum{
     OBOff = 0,
     OBOn,
     OBSlowish,
     OBStrong,
     OBForce,
     OBSimple,
     OBMorePrinting,
     OBEndMarker
  } OrbitalStrategy;

  /*! \brief What SOS prioritization strategy to use

    - SOSOff:
    - SOSHigh:
    - SOSLow:
    - SOSOrderHigh:
    - SOSOrderLow:
    - SOSEndMarker
   */

   typedef enum{
      SOSOff = 0,
      SOSHigh,
      SOSLow,
      SOSOrderHigh,
      SOSOrderLow,
      SOSEndMarker
   } SOSStrategy;

  /*! \brief What overall strategy to use

    - StrategyDefaulty = 0,
    - StrategyEasy,
    - StrategyAggressive,
    - StrategyEndMarker

  */

  typedef enum{
     StrategyDefault = 0,
     StrategyEasy,
     StrategyAggressive,
     StrategyEndMarker
  } StrategyMode;

  /*! \brief What clock type to use

    - Clock Cpu: Use CPU time
    - ClockElapsed: Use elapsed time

  */

   typedef enum{
      ClockCpu = 0,
      ClockElapsed,
      ClockEndMarker
   } ClockType;

  /*! \brief What clock type to use

    - CGraphOff:
    - CGraphOn:
    - CGraphClique:
    - CGraphEndMarker

  */

  typedef enum{
     CGraphOff = 0,
     CGraphOn,
     CGraphClique,
     CGraphEndMarker
  } CGraphMode;
   
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
      - CGForceBut: the cut generator will be installed with settings that force
    	   it to be called at every node, but more active at root (probing
    	   only)
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

  typedef enum {
     CGOff = 0,
     CGOn,
     CGRoot,
     CGIfMove,
     CGForceOn,
     CGForceBut,
     CGEndOnly,
     CGEndClean,
     CGEndBoth,
     CGLlong,
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
     CGOnlyInstead,
     CGCleanInstead,
     CGBothInstead,
     CGOnGlobal,
     CGForceAndGlobal,
     CGShorter,
     CGEndMarker
  } CGMode;

  /*! \brief Codes to specify whether to use a cutoff constraint
   */

   typedef enum{
      COOff = 0,
      COOn,
      COVariable,
      COForceVariable,
      COConflict,
      COEndMarker
   } CutoffMode;

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

   typedef enum {
      HeurOff = 0,
      HeurOn,
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
   } HeurMode;

   /*! \brief Codes to specify one or off for binary parameters

      - ParamOff: Capability is switched off 
      - ParamOn: Capability is switched on 
    */

   typedef enum {
      ParamOff = 0,
      ParamOn,
      ParamEndMarker
   } OnOffMode;
   
  /*! \brief Codes to specify the assignment of branching priorities

      - BPOff: no priorities are passed to cbc
      - BPCost: a priority vector is constructed based on objective coefficients
      - BPOrder: a priority vector is constructed based on column order
      - BPExt: the user has provided a priority vector
    */

  typedef enum {
     BPOff = 0,
     BPCost,
     BPOrder,
     BPExt,
     BPEndMarker
  } BPMode;

  /*! \brief Codes tos pecify mode for printing integers

    - PMNormal = 0,
    - PMInteger,
    - PMSpecial,
    - PMRows,
    - PMAll,
    - PMEndMarker

  */

  typedef enum{
     PMNormal = 0,
     PMInteger,
     PMSpecial,
     PMRows,
     PMAll,
     PMEndMarker
  } IntPrintMode;
   
  /*! \brief Major status codes for branch-and-cut

      - BACInvalid: status not yet set
      - BACNotRun: branch-and-cut has not yet run for the current problem
      - BACFinish: branch-and-cut has finished normally
      - BACStop: branch-and-cut has stopped on a limit
      - BACAbandon: branch-and-cut abandoned the problem
      - BACUser: branch-and-cut stopped on user signal

      Consult minorStatus_ for details.

      These codes are (mostly) set to match the codes used by CbcModel. Additions
      to CbcModel codes should be reflected here and in translateMajor.
    */

  typedef enum {
     BACInvalid = -1,
     BACFinish = 0,
     BACStop = 1,
     BACAbandon = 2,
     BACNotRun,
     BACUser = 5,
     BacEndMarker
  } BACMajorStatus;

  /*! \brief Minor status codes

      - BACmInvalid		status not yet set
      - BACmFinish		search exhausted the tree; optimal solution found
      - BACmInfeas		problem is infeasible
      - BACmUbnd		problem is unbounded
      - BACmGap		stopped on integrality gap
      - BACmNodeLimit	stopped on node limit
      - BACmTimeLimit	stopped on time limit
      - BACmSolnLimit	stopped on number of solutions limit
      - BACmUser		stopped due to user event
      - BACmOther		nothing else is appropriate

      It's not possible to make these codes agree with CbcModel. The meaning varies
      according to context: if the BACWhere code specifies a relaxation, then the
      minor status reflects the underlying OSI solver. Otherwise, it reflects the
      integer problem.
    */

  typedef enum {
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
  } BACMinorStatus;

  /*! \brief Codes to specify where branch-and-cut stopped

      - BACwNotStarted	stopped before we ever got going
      - BACwBareRoot	stopped after initial solve of root relaxation
      - BACwIPP		stopped after integer preprocessing
      - BACwIPPRelax	stopped after initial solve of preprocessed problem
      - BACwBAC		stopped at some point in branch-and-cut
    */

  typedef enum {
     BACwInvalid = -1,
     BACwNotStarted = 0,
     BACwBareRoot,
     BACwIPP,
     BACwIPPRelax,
     BACwBAC,
     BACwEndMarker
  } BACWhere;

  //@}

  /*! \name Constructors and destructors */
  //@{

  /*! \brief Default constructor */

  CbcSettings();

  /*! \brief Destructor */

  ~CbcSettings();
  //@}

  /*! \name Access and Control Functions for Cut Generators
        \brief Control functions, plus lazy creation functions for cut generators

        CBC avoids creating objects for cut generators and heuristics
        unless they're actually used. For cut generators, a prototype is created
        and reused. For heuristics, the default is to create a new object with each
        call, because the model may have changed. The object is returned through
        the reference parameter. The return value of the function is the current
        mode.

        Cut generator and heuristic objects created by these calls will be deleted
        with the destruction of the CbcSettings object.
    */
  //@{

  /*! \brief Get cut depth setting

      The name is a bit of a misnomer. Essentially, this overrides the
      `every so many nodes' control with `execute when (depth in tree)
      mod (cut depth) == 0'.
    */

  inline int getCutDepth()
  {
    return cutDepth_;
  }

  /*! \brief Get cut length setting

  */

  inline int getCutLength()
  {
    return cutLength_;
  }

  /*! \brief Get cut length setting

  */

  inline int getCutPassInTree()
  {
    return cutPassInTree_;
  }

  /*! \brief Set cut depth setting.

      See comments for getCutDepth().
  */

  inline void setCutDepth(int cutDepth)
  {
    cutDepth_ = cutDepth;
  }

  /*! \brief Set cut length setting.

    */

  inline void setCutLength(int cutLength)
  {
    cutLength_ = cutLength;
  }

  /*! \brief Set cut pass in tree setting.

    */

  inline void setCutPassInTree(int cutPassInTree)
  {
    cutPassInTree_ = cutPassInTree;
  }

  /*! \brief Obtain a prototype for a probing cut generator. */

  CGControl getProbing(CglCutGenerator *&gen);

  /*! \brief Set mode for use of probing cut generator. */

  inline void setProbingMode(CGControl mode)
  {
    probing_.mode_ = mode;
  }

  /*! \brief Obtain a prototype for a clique cut generator. */

  CGControl getClique(CglCutGenerator *&gen);

  /*! \brief Set mode for use of clique cut generator. */

  inline void setCliqueMode(CGControl mode)
  {
    clique_.mode_ = mode;
  }

  /*! \brief Obtain a prototype for a flow cover cut generator. */

  CGControl getFlow(CglCutGenerator *&gen);

  /*! \brief Set mode for use of flow cover cut generator. */

  inline void setFlowMode(CGControl mode)
  {
    flow_.mode_ = mode;
  }

  /*! \brief Obtain a prototype for a Gomory cut generator. */

  CGControl getGomory(CglCutGenerator *&gen);

  /*! \brief Set mode for use of Gomory cut generator. */

  inline void setGomoryMode(CGControl mode)
  {
    gomory_.mode_ = mode;
  }

  /*! \brief Obtain a prototype for a knapsack cover cut generator. */

  CGControl getKnapsack(CglCutGenerator *&gen);

  /*! \brief Set mode for use of knapsack cut generator. */

  inline void setKnapsackMode(CGControl mode)
  {
    knapsack_.mode_ = mode;
  }

   /*  \brief Obtain a prototype for a lift-and-project cut generator. */

  CGControl getLandP(CglCutGenerator *&gen) ;

   /* \brief Set mode for use of lift-and-project cut generator. */

   inline void setLandPMode(CGControl mode)
   {
      landp_.mode_ = mode ;
   }

  /*! \brief Obtain a prototype for a mixed integer rounding (MIR)
           cut generator.
    */

  CGControl getMir(CglCutGenerator *&gen);

  /*! \brief Set mode for use of MIR cut generator. */

  inline void setMirMode(CGControl mode)
  {
    mir_.mode_ = mode;
  }

  /*! \brief Obtain a prototype for a reduce and split cut generator. */

  CGControl getRedSplit(CglCutGenerator *&gen);

  /*! \brief Set mode for use of reduce and split cut generator. */

  inline void setRedSplitMode(CGControl mode)
  {
    redSplit_.mode_ = mode;
  }

  /*! \brief Obtain a prototype for a 2-MIR cut generator. */

  CGControl getTwomir(CglCutGenerator *&gen);

  /*! \brief Set mode for use of 2-MIR cut generator. */

  inline void setTwomirMode(CGControl mode)
  {
    twomir_.mode_ = mode;
  }
   
  //@}

  /*! \name Access and Control Functions for Cut Generators
        \brief Control functions, plus lazy creation functions for cut generators

        CBC avoids creating objects for cut generators and heuristics
        unless they're actually used. For cut generators, a prototype is created
        and reused. For heuristics, the default is to create a new object with each
        call, because the model may have changed. The object is returned through
        the reference parameter. The return value of the function is the current
        mode.

        Cut generator and heuristic objects created by these calls will be deleted
        with the destruction of the CbcSettings object.
  */
  //@{
  
  /*! \brief Obtain a local search/combine heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getCombine(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of local search/combine heuristic. */

  inline void setCombineMode(CGControl mode)
  {
    combine_.mode_ = mode;
  }

  /*! \brief Obtain a crossover heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getCrossover(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of crossover heuristic. */

  inline void setCrossoverMode(CGControl mode)
  {
    crossover_.mode_ = mode;
  }

  /*! \brief Obtain a DINS heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getDins(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of DINS heuristic. */

  inline void setDinsMode(CGControl mode)
  {
    dins_.mode_ = mode;
  }

  /*! \brief Obtain a Diving Coefficient heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getDiveCofficient(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of Diving Coefficient heuristic. */

  inline void setDiveCoefficientMode(CGControl mode)
  {
    divingc_.mode_ = mode;
  }

  /*! \brief Obtain a Diving Fractional heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getDiveFractional(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of Diving Coefficient heuristic. */

  inline void setDiveFractionalMode(CGControl mode)
  {
    divingf_.mode_ = mode;
  }

  /*! \brief Obtain a Diving Guided heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getDiveGuided(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of Diving Guided heuristic. */

  inline void setDiveGuidedMode(CGControl mode)
  {
    divingg_.mode_ = mode;
  }

  /*! \brief Obtain a Diving Line Search heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getDiveLineSearch(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of Diving Line Search heuristic. */

  inline void setDiveLineSearchMode(CGControl mode)
  {
    divingl_.mode_ = mode;
  }

  /*! \brief Obtain a Diving Pseudocost heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getDivePseudocost(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of Diving Pseudocost heuristic. */

  inline void setDivePseudocostMode(CGControl mode)
  {
    divingp_.mode_ = mode;
  }

  /*! \brief Obtain a Diving Vector Length heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getDiveVectorLength(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of Diving Vector Length heuristic. */

  inline void setDiveVectorLengthMode(CGControl mode)
  {
    divingv_.mode_ = mode;
  }

  /*! \brief Obtain a DW heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getDW(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of DW heuristic. */

  inline void setDWMode(CGControl mode)
  {
    dw_.mode_ = mode;
  }

  /*! \brief Obtain a feasibility pump heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
  */

  CGControl getFPump(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of feasibility pump heuristic. */

  inline void setFPumpMode(CGControl mode)
  {
    fpump_.mode_ = mode;
  }

  /*! \brief Set iterations for use of feasibility pump heuristic. */

  inline void setFPumpIters(int iters)
  {
    fpump_.iters_ = iters;
  }

  /*! \brief Set tune mode for use of feasibility pump heuristic. */

  inline void setFPumpTune(int tune)
  {
    fpump_.tune_ = tune;
  }

  /*! \brief Set second tune mode for use of feasibility pump heuristic. */

  inline void setFPumpTune2(int tune2)
  {
    fpump_.tune2_ = tune2;
  }

  /*! \brief Set fake cutoff for use of feasibility pump heuristic. */

  inline void setFPumpFakeCutoff(double cutoff)
  {
    fpump_.cutoff_ = cutoff;
  }
  /*! \brief Set fake increment for use of feasibility pump heuristic. */

  inline void setFPumpFakeIncrement(double increment)
  {
    fpump_.increment_ = increment;
  }

  /*! \brief Obtain a greedy cover heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getGreedyCover(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of greedy cover heuristic. */

  inline void setGreedyCoverMode(CGControl mode)
  {
    greedyCover_.mode_ = mode;
  }

  /*! \brief Obtain a greedy equality heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getGreedyEquality(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of greedy equality heuristic. */

  inline void setGreedyEqualityMode(CGControl mode)
  {
    greedyEquality_.mode_ = mode;
  }

  /*! \brief Obtain a Naive heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
  */

  CGControl getNaive(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of Naive heuristic. */

  inline void setNaiveMode(CGControl mode)
  {
    naive_.mode_ = mode;
  }

  /*! \brief Obtain a Pivot And Fix heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getPivotAndFix(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of Pivot and Fix heuristic. */

  inline void setPivotAndFixMode(CGControl mode)
  {
    pivotandfix_.mode_ = mode;
  }

  /*! \brief Obtain a Pivot and Complement heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getPivotAndComplement(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of Pivot and Complement heuristic. */

  inline void setPivotAndComplementMode(CGControl mode)
  {
    pivotandcomplement_.mode_ = mode;
  }

  /*! \brief Obtain a Proximity heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getProximity(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of Proximity heuristic. */

  inline void setProximityMode(CGControl mode)
  {
    proximity_.mode_ = mode;
  }

  /*! \brief Obtain a Randomized Rounding heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getRandRound(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of Randomized Rounding heuristic. */

  inline void setRandRoundMode(CGControl mode)
  {
    randround_.mode_ = mode;
  }

  /*! \brief Obtain a RENS heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getRens(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of RENS heuristic. */

  inline void setRensMode(CGControl mode)
  {
    rens_.mode_ = mode;
  }

  /*! \brief Obtain a RINS heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getRins(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of RENS heuristic. */

  inline void setRensMode(CGControl mode)
  {
    rens_.mode_ = mode;
  }

  /*! \brief Obtain a simple rounding heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getRounding(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of simple rounding heuristic. */

  inline void setRoundingMode(CGControl mode)
  {
    rounding_.mode_ = mode;
  }

  /*! \brief Obtain a Variable Neighborhood heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getVnd(CbcHeuristic *&gen, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of Variable Neighborhood heuristic. */

  inline void setVndMode(CGControl mode)
  {
    vnd_.mode_ = mode;
  }

  /*! \brief Obtain a local search tree object

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

  CGControl getTreeLocal(CbcTreeLocal *&localTree, CbcModel *model,
    bool alwaysCreate = true);

  /*! \brief Set mode for use of local tree. */

  inline void setTreeLocalMode(CGControl mode)
  {
    localTree_.mode_ = mode;
  }

  //@}

  /*! \name Miscellaneous Integer Parameters

  */

  //@{

  /*! \brief Get BkPivotStrategy setting

    */

  inline int getBkPivotStrategy()
  {
    return bkPivotStrategy_;
  }

  /*! \brief Set BkPivotStrategy setting

    */

  inline void setBkPivotStrategy(int bkPivotStrategy)
  {
    bkPivotStrategy_ = bkPivotStrategy;
  }

  /*! \brief Get BkMaxCalls setting

    */

  inline int getBkMaxCalls()
  {
    return bkMaxCalls_;
  }

  /*! \brief Set BkMaxCalls setting

    */

  inline void setBkMaxCalls(int bkMaxCalls)
  {
    bkMaxCalls_ = bkMaxCalls;
  }

  /*! \brief Get BkClqExtMethod setting

    */

  inline int getBkClqExtMethod()
  {
    return bkClqExtMethod_;
  }

  /*! \brief Set BkClqExtMethod setting

    */

  inline void setBkClqExtMethod(int bkClqExtMethod)
  {
    bkClqExtMethod_ = bkClqExtMethod;
  }

  /*! \brief Get CppMode setting

    */

  inline int getCppMode()
  {
    return cppMode_;
  }

  /*! \brief Set CppMode setting

    */

  inline void setCppMode(int cppMode)
  {
    cppMode_ = cppMode;
  }

  /*! \brief Get DepthMiniBaB setting

    */

  inline int getDepthMiniBaB()
  {
    return depthMiniBaB_;
  }

  /*! \brief Set DepthMiniBaB setting

    */

  inline void setDepthMiniBaB(int depthMiniBaB)
  {
    depthMiniBaB_ = depthMiniBaB;
  }

  /*! \brief Get DiveOpt setting

    */

  inline int getDiveOpt()
  {
    return diveOpt_;
  }

  /*! \brief Set DiveOpt setting

    */

  inline void setdiveOpt(int diveOpt)
  {
    diveOpt_ = diveOpt;
  }

  /*! \brief Get DiveOptSolves setting

    */

  inline int getDiveOptSolves()
  {
    return diveOptSolves_;
  }

  /*! \brief Set DiveOptSolves setting

    */

  inline void setDiveOptSolves(int diveOptSolves)
  {
    diveOptSolves_ = diveOptSolves;
  }

  /*! \brief Get ExperimentMode setting

    */

  inline int getExperimentMode()
  {
    return experimentMode_;
  }

  /*! \brief Set ExperimentMode setting

    */

  inline void setExperimentMode(int experimentMode)
  {
    experimentMode_ = experimentMode;
  }

  /*! \brief Get ExtraIntParam1 setting

    */

  inline int getExtraIntParam1()
  {
    return extraInt1_;
  }

  /*! \brief Set ExtraIntParam1 setting

    */

  inline void setExtraIntParam1(int extraInt1)
  {
    extraInt1_ = extraInt1;
  }

  /*! \brief Get ExtraIntParam2 setting

    */

  inline int getExtraIntParam2()
  {
    return extraInt2_;
  }

  /*! \brief Set ExtraIntParam2 setting

    */

  inline void setExtraIntParam2(int extraInt2)
  {
    extraInt2_ = extraInt2;
  }

  /*! \brief Get ExtraIntParam3 setting

    */

  inline int getExtraIntParam3()
  {
    return extraInt3_;
  }

  /*! \brief Set ExtraIntParam3 setting

    */

  inline void setExtraIntParam3(int extraInt3)
  {
    extraInt3_ = extraInt3;
  }

  /*! \brief Get ExtraIntParam4 setting

    */

  inline int getExtraIntParam4()
  {
    return extraInt4_;
  }

  /*! \brief Set ExtraIntParam4 setting

    */

  inline void setExtraIntParam4(int extraInt4)
  {
    extraInt4_ = extraInt4;
  }

  /*! \brief Get HeurOptions setting

    */

  inline int getHeurOptions()
  {
    return heurOptions_;
  }

  /*! \brief Set HeurOptions setting

    */

  inline void setHeurOptions(int heurOptions)
  {
    heurOptions_ = heurOptions;
  }

  /*! \brief Get MaxSavedSols setting

    */

  inline int getMaxSavedSols()
  {
    return maxSavedSols_;
  }

  /*! \brief Set MaxSavedSols setting

    */

  inline void setMaxSavedSols(int maxSavedSols)
  {
    maxSavedSols_ = maxSavedSols;
  }

  /*! \brief Get MaxSlowCuts setting

    */

  inline int getMaxSlowCuts()
  {
    return maxSlowCuts_;
  }

  /*! \brief Set MaxSlowCuts setting

    */

  inline void setMaxSlowCuts(int maxSlowCuts)
  {
    maxSlowCuts_ = maxSlowCuts;
  }

  /*! \brief Get MoreMoreOptions setting

    */

  inline int getMoreMoreOptions()
  {
    return moreMoreOptions_;
  }

  /*! \brief Set MoreMoreOptions setting

    */

  inline void setMoreMoreOptions(int moreMoreOptions)
  {
    moreMoreOptions_ = moreMoreOptions;
  }

  /*! \brief Get MultipleRoots setting

    */

  inline int getMultipleRoots()
  {
    return multipleRoots_;
  }

  /*! \brief Set MultipleRoots setting

    */

  inline void setMultipleRoots(int multipleRoots)
  {
    multipleRoots_ = multipleRoots;
  }

  /*! \brief Get OddWextMethod setting

    */

  inline int getOddWextMethod()
  {
    return oddWextMethod_;
  }

  /*! \brief Set OddWextMethod setting

    */

  inline void setOddWextMethod(int oddWextMethod)
  {
    oddWextMethod_ = oddWextMethod;
  }

  /*! \brief Get OutputFormat setting

    */

  inline int getOutputFormat()
  {
    return outputFormat_;
  }

  /*! \brief Set OutputFormat setting

    */

  inline void setOutputFormat(int outputFormat)
  {
    outputFormat_ = outputFormat;
  }

  /*! \brief Get PrintOptions setting

    */

  inline int getPrintOptions()
  {
    return printOpt_;
  }

  /*! \brief Set PrintOptions setting

    */

  inline void setPrintOptions(int printOpt)
  {
    printOpt_ = printOpt;
  }

  /*! \brief Get ProcessTune setting

    */

  inline int getProcessTune()
  {
    return processTune_;
  }

  /*! \brief Set ProcessTune setting

    */

  inline void setProcessTune(int processTune)
  {
    processTune_ = processTune;
  }

  /*! \brief Get RandomSeed setting

    */

  inline int getRandomSeed()
  {
    return randomSeed_;
  }

  /*! \brief Set RandomSeed setting

    */

  inline void setRandomSeed(int randomSeed)
  {
    randomSeed_ = randomSeed;
  }

  /*! \brief Get StrongStrategy setting

    */

  inline int getStrongStrategy()
  {
    return strongStrategy_;
  }

  /*! \brief Set StrongStrategy setting

    */

  inline void setStrongStrategy(int strongStrategy)
  {
    strongStrategy_ = strongStrategy;
  }
   
  /*! \brief Get TestOsi setting

    */

  inline int getTestOsi()
  {
    return testOsi_;
  }

  /*! \brief Set TestOsi setting

    */

  inline void setTestOsi(int testOsi)
  {
    testOsi_ = testOsi;
  }
   
  /*! \brief Get Threads setting

    */

  inline int getThreads()
  {
    return threads_;
  }

  /*! \brief Set Threads setting

    */

  inline void setThreads(int threads)
  {
    threads_ = threads;
  }
   
  /*! \brief Get UserCbc setting

    */

  inline int getUserCbc()
  {
    return userCbc_;
  }

  /*! \brief Set UserCbc setting

  */

  inline void setUserCbc(int userCbc)
  {
    userCbc_ = userCbc;
  }

  /*! \brief Get Verbose setting

    */

  inline int getVerbose()
  {
    return verbose_;
  }

  /*! \brief Set Verbose setting

    */

  inline void setVerbose(int verbose)
  {
    verbose_ = verbose;
  }
   
  /*! \brief Get VubTry setting

    */

  inline int getVubTry()
  {
    return vubTry_;
  }

  /*! \brief Set VubTry setting

    */

  inline void setVubTry(int vubTry)
  {
    vubTry_ = vubTry;
  }

  //@}
   
  /*! \name Miscellaneous Double Parameters

  */

  //@{

  /*! \brief Get extra double 3 setting

    */

  inline double getExtraDbl3()
  {
    return extraDbl3_;
  }

  /*! \brief Set extra double 3 setting

    */

  inline void setExtraDbl3(double extraDbl3)
  {
    extraDbl3_ = extraDbl3;
  }
   
  /*! \brief Get extra double 3 setting

    */

  inline double getExtraDbl3()
  {
    return extraDbl3_;
  }

  /*! \brief Set extra double 4 setting

    */

  inline void setExtraDbl4(double extraDbl4)
  {
    extraDbl4_ = extraDbl4;
  }
   
  /*! \brief Get extra double 3 setting

    */

  inline double getExtraDbl5()
  {
    return extraDbl5_;
  }

  /*! \brief Set extra double 5 setting

    */

  inline void setExtraDbl5(double extraDbl5)
  {
    extraDbl5_ = extraDbl5;
  }
   
  /*! \brief Get small branch and bound setting

    */

  inline double getSmallBaB()
  {
    return smallBaB_;
  }

  /*! \brief Set small branch and bound setting

    */

  inline void setSmallBab(double smallBab)
  {
    smallBab_ = smallBab;
  }
   
  /*! \brief Get tighten factor

    */

  inline double getTightenFactor()
  {
    return tightenFactor_;
  }

  /*! \brief Set tighten factor

  */

  inline void set(double tightenFactor)
  {
    tightenFactor_ = tightenFactor;
  }
   
  //@}
   
  /*! \name Miscellaneous Keyword Parameters

  */

  //@{

  /*! \brief Get command mode */

  inline CommandMode getCommandMode()
  {
    return (commandMode_);
  }

  /*! \brief Set command mode */

  inline void setCommandMode(CommandMode mode)
  {
    commandMode_ = mode;
  }

  /*! \brief Get mode for clique strengthening */

  inline ClqStrMode getClqStrMode()
  {
    return (clqStrMode_);
  }

  /*! \brief Set mode for use of integer preprocessing */

  inline void setClqStrMode(ClqStrMode mode)
  {
    clqStrMode_ = mode;
  }

  /*! \brief Get mode for branching priorities */

  inline BPMode getBranchPriority()
  {
    return (branchPriority_);
  }

  /*! \brief Set mode for branching priorities */

  inline void setBranchPriority(BPMode mode)
  {
    branchPriority_ = mode;
  }

  /*! \brief Get mode for use of cutoff constraint */

  inline CutoffMode getCutoffMode()
  {
    return (cutoffMode_);
  }

  /*! \brief Set mode for use of cutoff constraint */

  inline void setCutoff(CutoffMode mode)
  {
    cutoffMode_ = mode;
  }

  /*! \brief Get mode for printing integers */

  inline IntPrintMode getIntPrintMode()
  {
    return (intPrintMode_);
  }

  /*! \brief Set  mode for printing integers */

  inline void setIntPrintMode(IntPrintMode mode)
  {
    intPrintMode_ = mode;
  }

  /*! \brief Get node search strategy */

  inline NodeStrategy getNodeStrategy()
  {
    return (nodeStrategy_);
  }

  /*! \brief Set node search strategy */

  inline void setNodeStrategy(NodeStrategy mode)
  {
    nodeStrategy_ = mode;
  }

  /*! \brief Get strategy for orbital branching */

   inline OrbitalStrategy getOrbitalStrategy()
  {
    return (orbitalStrategy_);
  }

  /*! \brief Set  strategy for orbital branching */

  inline void setOrbitalStrategy(OrbitalStrategy mode)
  {
    orbitalStrategy_ = mode;
  }
   
  /*! \brief Get mode for use of integer preprocessing */

  inline IPPMode getIPPMode()
  {
    return (preProcess_);
  }

  /*! \brief Set mode for use of integer preprocessing */

  inline void setIPPMode(IPPMode mode)
  {
    preProcess_ = mode;
  }
   
  /*! \brief Get priority mode for SOS */

  inline SOSMode getSOSMode()
  {
    return (sosMode_);
  }

  /*! \brief Set mode state for use of integer preprocessing */

  inline void setSOSMode(SOSMode mode)
  {
    SOSMode_ = mode;
  }
   
  /*! \brief Get overall strategy mode */

  inline StrategyMode getStrategyMode()
  {
    return (strategyMode_);
  }

  /*! \brief Set overall strategy mode */

  inline void setStrategyMode(StrategyMode mode)
  {
    strategyMode_ = mode;
  }
   
  /*! \brief Get clock type */

  inline ClockType getClockType()
  {
    return (clockType_);
  }

  /*! \brief Set clock type */

  inline void setClockType(ClockType type)
  {
    clockType_ = type;
  }
   
  /*! \brief Get mode for CGraph */

  inline CGraphMode getCGraphMode()
  {
    return (cgraphMode_);
  }

  /*! \brief Set mode for CGraph */

  inline void setCGraphMode(CGraphMode mode)
  {
    CGraphMode_ = mode;
  }
  //@}


  /*! \name Status Functions
        \brief Convenience routines for status codes.
    */
  //@{

  /*! \brief Set the result of branch-and-cut search */

  inline void setBaBStatus(BACMajor majorStatus, BACMinor minorStatus,
    BACWhere where, bool haveAnswer,
    OsiSolverInterface *answerSolver)
  {
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
  void setBaBStatus(const CbcModel *model, BACWhere where,
    bool haveAnswer = false,
    OsiSolverInterface *answerSolver = 0);

  /*! \brief Translate CbcModel major status to #BACMajor

      See the #BACMajor enum for details.
    */
  BACMajor translateMajor(int status);

  /*!\brief Translate CbcModel minor status to #BACMinor

      See the #BACMinor enum for details.
    */
  BACMinor translateMinor(int status);

  /*!\brief Translate OsiSolverInterface status to #BACMinor

      See the #BACMinor enum for details. Optimal, infeasible, and unbounded
      get their own codes; everything else maps to BACmOther.
    */
  BACMinor translateMinor(const OsiSolverInterface *osi);

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
  inline CoinMessageHandler *messageHandler() const
  {
    return msgHandler_;
  }

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
  inline void setLogLevel(int lvl)
  {
    logLvl_ = lvl;
    if (msgHandler_)
      msgHandler_->setLogLevel(lvl);
  }

  /*! \brief Get log level */
  inline int logLevel() const
  {
    return (logLvl_);
  }

  /*! \brief When greater than 0, integer presolve gives more information and
           branch-and-cut provides statistics.
    */
  int printOpt_;

  //@}

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

  /*! \brief Last solution output file */

  std::string lastSolnOut_;

  /*! \brief Solution printing mode

      Controls the amount of information printed when printing a solution.
      Coding is set by the keyword declarations for the printingOptions
      command.
    */
  int printMode_;

  /*! \brief Print mask

      Used to specify row/column names to be printed. Not implemented as of
      060920.
    */
  std::string printMask_;

  /*! \brief The parameter vector */

  CoinParamVec *paramVec_;

  /*! \brief Start and end of cbc-generic parameters in parameter vector */

  struct genParamsInfo_struct {
    int first_;
    int last_;
  } genParams_;

  /*! \brief Start and end of CbcModel parameters in parameter vector */

  struct cbcParamsInfo_struct {
    int first_;
    int last_;
  } cbcParams_;

  /*! \brief Start and end of OsiSolverInterface  parameters in parameter
           vector
    */

  struct osiParamsInfo_struct {
    int first_;
    int last_;
  } osiParams_;

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

  std::vector< bool > setByUser_;

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
      we have a valid answer. See the documentation with the BACMajor,
      BACMinor, and BACWhere enums for the meaning of the codes.
    */

  struct babState_struct {
    BACMajor majorStatus_;
    BACMinor minorStatus_;
    BACWhere where_;
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
    bool on_;
    double threshold_;
  } djFix_;

  /*! \brief Treat some variables as artificial in feasibility pump

    */

  struct artVar_struct {
    bool on_;
    double threshold_;
  } artVar_;

  /*! \brief Control the assignment of branching priorities to integer
           variables.
    */
  BPControl priorityMode_;

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

private:
  /*! \name Cut Generator Control
        \brief Usage control and prototypes for cut generators.
    */
  //@{

  /*! \brief Control integer preprocessing. */

  IPPMode preProcess_;

  /*! \brief Control cut generator activity

      Generators that are active in the tree will be activated when
      (depth) mod (cutDepth) == 0.
    */

  int cutDepth_;
  int cutLength_;
  int cutPassInTree_;

  /*! \brief Control variable and prototype for clique cut generator */
  struct clique_struct {
    CGControl mode_;
    CglClique *proto_;
    bool starCliqueReport_;
    bool rowCliqueReport_;
    double minViolation_;
  } clique_;

  /*! \brief Control variable and prototype for flow cover cut generator */
  struct flow_struct {
    CGControl mode_;
    CglFlowCover *proto_;
  } flow_;

  /*! \brief Control variable and prototype for Gomory cut generator */
  struct gmi_struct {
    CGControl mode_;
    CglGomory *proto_;
    int limit_;
    int limitAtRoot_;
  } gmi_;

  /*! \brief Control variable and prototype for Gomory cut generator */
  struct gomory_struct {
    CGControl mode_;
    CglGomory *proto_;
    int limit_;
    int limitAtRoot_;
  } gomory_;

  /*! \brief Control variable and prototype for knapsack cover cut generator */
  struct knapsack_struct {
    CGControl mode_;
    CglKnapsackCover *proto_;
  } knapsack_;

  /*   \brief Control variable and prototype for LaGomory cut
    	     generator
    */
  struct lagomory_struct{
     CGControl mode_ ;
     CglGomory *proto_ ;
  } lagomory_ ;
   
  /*   \brief Control variable and prototype for lift-and-project cut
    	     generator
    */
  struct landp_struct{
     CGControl mode_ ;
     CglLandP *proto_ ;
  } landp_ ;
   
  /*   \brief Control variable and prototype for lift-and-project cut
    	     generator
    */
  struct laTwoMir_struct{
     CGControl mode_ ;
     CglTwoMir *proto_ ;
  } latwomir_ ;
   
  /*! \brief Control variable and prototype for MIR cut generator */
  struct mir_struct {
    CGControl mode_;
    CglMixedIntegerRounding2 *proto_;
  } mir_;

  /*! \brief Control variable and prototype for odd hole cut generator */
  struct oddWheel_struct {
    CGControl mode_;
    CglOddWheel *proto_;
  } oddWheel_;

   /*! \brief Control variable and prototype for probing cut generator */
  struct probing_struct {
    CGControl mode_;
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

  /*! \brief Control variable and prototype for reduce-and-split
           cut generator
    */
  struct redSplit_struct {
    CGControl mode_;
    CglRedSplit *proto_;
  } redSplit_;

  /*! \brief Control variable and prototype for reduce-and-split 2
           cut generator
    */
  struct redSplit2_struct {
    CGControl mode_;
    CglRedSplit2 *proto_;
  } redSplit2_;

  /*! \brief Control variable and prototype for residual capacity
           cut generator
    */
  struct residCap_struct {
    CGControl mode_;
    CglResidualCapacity *proto_;
  } residcap_;

  /*! \brief Control variable and prototype for Two-MIR cut generator */
  struct twomir_struct {
    CGControl mode_;
    CglTwomir *proto_;
    int maxElements_;
  } twomir_;

  /*! \brief Control variable and prototype for residual capacity
           cut generator
    */
  struct zeroHalf_struct {
    CGControl mode_;
    CglZeroHalf *proto_;
  } zerohalf_;

  //@}
   
  /*! \name Heuristic Control
        \brief Usage control and prototypes for heuristics.
  */
  //@{

  /*! \brief Control variable and prototype for combine heuristic */
  struct combine_struct {
    HeurControl mode_;
    CbcHeuristicLocal *proto_;
    int trySwap_;
  } combine_;

  /*! \brief Control variable and prototype for crossover heuristic */
  struct crossover_struct {
    HeurControl mode_;
    CbcHeuristicCrossover *proto_;
  } crossover_;

  /*! \brief Control variable and prototype for heuristic */
  struct dins_struct {
    HeurControl mode_;
    CbcHeuristicDINS *proto_;
  } dins_;

  /*! \brief Control variable and prototype for Dive Coefficient heuristic */
  struct diveCoefficient_struct {
    HeurControl mode_;
    CbcHeuristicDiveCoefficient *proto_;
  } divec_;

  /*! \brief Control variable and prototype for Dive Fractional heuristic */
  struct diveFractional_struct {
    HeurControl mode_;
    CbcHeuristicDiveFractional *proto_;
  } divef_;

  /*! \brief Control variable and prototype for Dive Guided heuristic */
  struct diveGuided_struct {
    HeurControl mode_;
    CbcHeuristicDiveGuided *proto_;
  } diveg_;

  /*! \brief Control variable and prototype for Dive Line Sarch heuristic */
  struct diveLineSearch_struct {
    HeurControl mode_;
    CbcHeuristicDiveLineSearch *proto_;
  } divel_;

  /*! \brief Control variable and prototype for Dive Pseudocost heuristic */
  struct _struct {
    HeurControl mode_;
    CbcHeuristicDivePseudocost *proto_;
  } divep_;

  /*! \brief Control variable and prototype for Dive Vector Lengthheuristic */
  struct diveVectorLength_struct {
    HeurControl mode_;
    CbcHeuristicDiveVectorLength *proto_;
  } divev_;

  /*! \brief Control variable and prototype for DW heuristic */
  struct DW_struct {
    HeurControl mode_;
    CbcHeuristicDW *proto_;
  } dw_;

  /*! \brief Control variable and prototype for feasibility pump heuristic */
  struct fpump_struct {
    HeurControl mode_;
    CbcHeuristicFPump *proto_;
    int iters_;
    int tune_;
    int tune2_;
    double fakeCutoff_;
    doublefakeIncrement_
  } fpump_;

  /*! \brief Control variable and prototype for greedy cover heuristic */
  struct greedyCover_struct {
    HeurControl mode_;
    CbcHeuristicGreedyCover *proto_;
  } greedyCover_;

  /*! \brief Control variable and prototype for greedy equality heuristic */
  struct greedyEquality_struct {
    HeurControl mode_;
    CbcHeuristicGreedyEquality *proto_;
  } greedyEquality_;

  /*! \brief Control variable and prototype for Naive heuristic */
  struct naive_struct {
    HeurControl mode_;
    CbcHeuristicNaive *proto_;
  } naive_;

  /*! \brief Control variable and prototype for Pivot and Fix heuristic */
  struct pivotAndFix_struct {
    HeurControl mode_;
    CbcHeuristicPivotAndFix *proto_;
  } pivotf_;

  /*! \brief Control variable and prototype for Pivot and Complement heuristic */
  struct pivotAndComp_struct {
    HeurControl mode_;
    CbcHeuristicPivtoAndComplement *proto_;
  } pivotc_;

  /*! \brief Control variable and prototype for Proximity heuristic */
  struct proximity_struct {
    HeurControl mode_;
    CbcHeuristicProximity *proto_;
  } proximity_;

  /*! \brief Control variable and prototype for Randomized Rounding heuristic */
  struct randRound_struct {
    HeurControl mode_;
    CbcHeuristic *proto_;
  } randround_;

  /*! \brief Control variable and prototype for RENS heuristic */
  struct Rens_struct {
    HeurControl mode_;
    CbcHeuristicRENS *proto_;
  } rens_;

  /*! \brief Control variable and prototype for RINS heuristic */
  struct Rins_struct {
    HeurControl mode_;
    CbcHeuristicRINS *proto_;
  } rins_;

  /*! \brief Control variable and prototype for simple rounding heuristic */
  struct rounding_struct {
    HeurControl mode_;
    CbcRounding *proto_;
  } rounding_;

  /*! \brief Control variable and prototype for Variable Neighborhood heuristic */
  struct Vnd_struct {
    HeurControl mode_;
    CbcHeuristicVND *proto_;
  } vnd_;

  /*! \brief Control variables for local tree

      This is a bit different --- getTreeLocal() takes a CbcModel as a parameter
      and installs a local tree object. But we can keep the parameters here and
      hide the details. Consult CbcTreeLocal.hpp for details.
    */
  struct localTree_struct {
    HeurControl mode_;
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
  double smallBab_;
  double tightenFactor_;

  //@}

  /*! \name Miscellaneous integer paramters

  */
  //@{

  int bkPivotStrategy_;
  int bkMaxCalls_;
  int bkClqExtMethod_;
  int cppMode_;
  int depthMinBab_;
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
  int StrongStrategy_;
  int testOsi_;
  int threads_;
  int userCbc_;
  int verbose_;
  int vubTry_;

  //@}
   
  /*! \name Miscellaneous keyword paramters

  */
  //@{

  CommandMode commandmode_;
  ClqStrMode clqStrMode_;
  BPMode branchPriority_;
  CutoffMode cutoffMode_;
  IntPrintMode intPrintMode_;
  NodeStrategy nodeStrategy_;
  OrbitalStrategy orbitalStrategy_;
  SOSMode sosMode_;
  StrategyMode strategyMode_;
  ClockType clockType_;
  CGraphMode cgraphMode_;
   
  //@}
   
  /*! \name Messages and statistics (private)
        \brief Data and objects related to messages and statistics that should be
    	   protected from direct manipulation.
    */
  //@{

  /*! \brief Message handler. */
  CoinMessageHandler *msgHandler_;

  /*! \brief Ownership of message handler.

      If true, the control block owns the message handler and it will be destroyed
      with the control block. If false, the client is responsible for the message
      handler.
    */
  bool ourMsgHandler_;

  /*! \brief The current language */
  CoinMessages::Language cur_lang_;

  /*! \brief The current set of messages. */
  CoinMessages *msgs_;

  /*! \brief The current log level */
  int logLvl_;

  //@}
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
