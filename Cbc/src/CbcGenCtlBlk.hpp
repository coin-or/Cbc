/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file is part of cbc-generic.
*/

#ifndef CbcGenCtlBlk_H
#define CbcGenCtlBlk_H

/* \file CbcGenCtlBlk.hpp
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
  `struct nameCtl_struct' style used for structured fields in CbcGenCtlBlk.
*/

/*
  $Id: CbcGenCtlBlk.hpp 1173 2009-06-04 09:44:10Z forrest $
*/

#define CBC_GENERIC_VERSION "00.01.00"

class CbcGenCtlBlk ;
namespace CbcGenParamUtils {
void addCbcGenParams(int &numParams, CoinParamVec &paramVec,
                     CbcGenCtlBlk *ctlBlk) ;
}

/* \brief cbc-generic algorithm control class

  This class defines values and methods used to control the operation of the
  cbc-generic main program.
*/

class CbcGenCtlBlk {

    friend void CbcGenParamUtils::addCbcGenParams(int &numParams,
            CoinParamVec &paramVec, CbcGenCtlBlk *ctlBlk) ;

public:

    /*! \name Enumeration types used for cbc-generic control variables */
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

    typedef enum { IPPOff = 0, IPPOn, IPPSave, IPPEqual,
                   IPPSOS, IPPTrySOS, IPPEqualAll, IPPStrategy
                 } IPPControl ;

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
      - CGMarker: a convenience to mark the end of the codes.

      The same codes are used for heuristics.
    */

    typedef enum { CGOff, CGOn, CGRoot, CGIfMove,
                   CGForceOn, CGForceBut, CGMarker
                 } CGControl ;

    /*! \brief Codes to specify the assignment of branching priorities

      - BPOff: no priorities are passed to cbc
      - BPCost: a priority vector is constructed based on objective coefficients
      - BPOrder: a priority vector is constructed based on column order
      - BPExt: the user has provided a priority vector
    */

    typedef enum { BPOff, BPCost, BPOrder, BPExt } BPControl ;

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

    typedef enum { BACInvalid = -1, BACFinish = 0,
                   BACStop = 1, BACAbandon = 2, BACNotRun, BACUser = 5
                 } BACMajor ;

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

    typedef enum { BACmInvalid = -1, BACmFinish = 0, BACmInfeas, BACmUbnd,
                   BACmGap, BACmNodeLimit, BACmTimeLimit, BACmSolnLimit,
                   BACmUser, BACmOther
                 } BACMinor ;

    /*! \brief Codes to specify where branch-and-cut stopped

      - BACwNotStarted	stopped before we ever got going
      - BACwBareRoot	stopped after initial solve of root relaxation
      - BACwIPP		stopped after integer preprocessing
      - BACwIPPRelax	stopped after initial solve of preprocessed problem
      - BACwBAC		stopped at some point in branch-and-cut
    */

    typedef enum { BACwInvalid = -1, BACwNotStarted = 0, BACwBareRoot,
                   BACwIPP, BACwIPPRelax, BACwBAC
                 } BACWhere ;

//@}

    /*! \name Constructors and destructors */
//@{

    /*! \brief Default constructor */

    CbcGenCtlBlk() ;

    /*! \brief Destructor */

    ~CbcGenCtlBlk() ;
//@}

    /*! \name Access and Control Functions for Cut Generators and Heuristics
        \brief Control functions, plus lazy creation functions for cut generators
    	   and heuristics

        cbc-generic avoids creating objects for cut generators and heuristics
        unless they're actually used. For cut generators, a prototype is created
        and reused. For heuristics, the default is to create a new object with each
        call, because the model may have changed. The object is returned through
        the reference parameter. The return value of the function is the current
        action state.

        Cut generator and heuristic objects created by these calls will be deleted
        with the destruction of the CbcGenCtlBlk object.
    */
//@{

    /*! \brief Get cut depth setting

      The name is a bit of a misnomer. Essentially, this overrides the
      `every so many nodes' control with `execute when (depth in tree)
      mod (cut depth) == 0'.
    */

    inline int getCutDepth() {
        return cutDepth_ ;
    }

    /*! \brief Set cut depth setting.

      See comments for getCutDepth().
    */

    inline void setCutDepth(int cutDepth) {
        cutDepth_ = cutDepth ;
    }

    /*1 \brief Get action state for use of integer preprocessing */

    inline IPPControl getIPPAction() {
        return (preProcess_) ;
    }

    /*! \brief Set action state for use of integer preprocessing */

    inline void setIPPAction(IPPControl action) {
        preProcess_ = action ;
    }

    /*! \brief Obtain a prototype for a probing cut generator. */

    CGControl getProbing(CglCutGenerator *&gen) ;

    /*! \brief Set action state for use of probing cut generator. */

    inline void setProbingAction(CGControl action) {
        probing_.action_ = action ;
    }

    /*! \brief Obtain a prototype for a clique cut generator. */

    CGControl getClique(CglCutGenerator *&gen) ;

    /*! \brief Set action state for use of clique cut generator. */

    inline void setCliqueAction(CGControl action) {
        clique_.action_ = action ;
    }

    /*! \brief Obtain a prototype for a flow cover cut generator. */

    CGControl getFlow(CglCutGenerator *&gen) ;

    /*! \brief Set action state for use of flow cover cut generator. */

    inline void setFlowAction(CGControl action) {
        flow_.action_ = action ;
    }

    /*! \brief Obtain a prototype for a Gomory cut generator. */

    CGControl getGomory(CglCutGenerator *&gen) ;

    /*! \brief Set action state for use of Gomory cut generator. */

    inline void setGomoryAction(CGControl action) {
        gomory_.action_ = action ;
    }

    /*! \brief Obtain a prototype for a knapsack cover cut generator. */

    CGControl getKnapsack(CglCutGenerator *&gen) ;

    /*! \brief Set action state for use of knapsack cut generator. */

    inline void setKnapsackAction(CGControl action) {
        knapsack_.action_ = action ;
    }

    /*  \brief Obtain a prototype for a lift-and-project cut generator.

      CGControl getLandP(CglCutGenerator *&gen) ;

       \brief Set action state for use of lift-and-project cut generator.

      inline void setLandPAction(CGControl action)
      { landp_.action_ = action ; }
    */

    /*! \brief Obtain a prototype for a mixed integer rounding (MIR)
           cut generator.
    */

    CGControl getMir(CglCutGenerator *&gen) ;

    /*! \brief Set action state for use of MIR cut generator. */

    inline void setMirAction(CGControl action) {
        mir_.action_ = action ;
    }

    /*! \brief Obtain a prototype for a reduce and split cut generator. */

    CGControl getRedSplit(CglCutGenerator *&gen) ;

    /*! \brief Set action state for use of reduce and split cut generator. */

    inline void setRedSplitAction(CGControl action) {
        redSplit_.action_ = action ;
    }

    /*! \brief Obtain a prototype for a 2-MIR cut generator. */

    CGControl getTwomir(CglCutGenerator *&gen) ;

    /*! \brief Set action state for use of 2-MIR cut generator. */

    inline void setTwomirAction(CGControl action) {
        twomir_.action_ = action ;
    }


    /*! \brief Obtain a feasibility pump heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

    CGControl getFPump(CbcHeuristic *&gen, CbcModel *model,
                       bool alwaysCreate = true) ;

    /*! \brief Set action state for use of feasibility pump heuristic. */

    inline void setFPumpAction(CGControl action) {
        fpump_.action_ = action ;
    }

    /*! \brief Obtain a local search/combine heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

    CGControl getCombine(CbcHeuristic *&gen, CbcModel *model,
                         bool alwaysCreate = true) ;

    /*! \brief Set action state for use of local search/combine heuristic. */

    inline void setCombineAction(CGControl action) {
        combine_.action_ = action ;
    }

    /*! \brief Obtain a greedy cover heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

    CGControl getGreedyCover(CbcHeuristic *&gen, CbcModel *model,
                             bool alwaysCreate = true) ;

    /*! \brief Set action state for use of greedy cover heuristic. */

    inline void setGreedyCoverAction(CGControl action) {
        greedyCover_.action_ = action ;
    }

    /*! \brief Obtain a greedy equality heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

    CGControl getGreedyEquality(CbcHeuristic *&gen, CbcModel *model,
                                bool alwaysCreate = true) ;

    /*! \brief Set action state for use of greedy equality heuristic. */

    inline void setGreedyEqualityAction(CGControl action) {
        greedyEquality_.action_ = action ;
    }

    /*! \brief Obtain a simple rounding heuristic.

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

    CGControl getRounding(CbcHeuristic *&gen, CbcModel *model,
                          bool alwaysCreate = true) ;

    /*! \brief Set action state for use of simple rounding heuristic. */

    inline void setRoundingAction(CGControl action) {
        rounding_.action_ = action ;
    }

    /*! \brief Obtain a local search tree object

      By default, any existing object is deleted and a new object is created and
      loaded with \c model. Set alwaysCreate = false to return an existing object
      if one exists.
    */

    CGControl getTreeLocal(CbcTreeLocal *&localTree, CbcModel *model,
                           bool alwaysCreate = true) ;

    /*! \brief Set action state for use of local tree. */

    inline void setTreeLocalAction(CGControl action) {
        localTree_.action_ = action ;
    }

//@}

    /*! \name Status Functions
        \brief Convenience routines for status codes.
    */
//@{

    /*! \brief Set the result of branch-and-cut search */

    inline void setBaBStatus(BACMajor majorStatus, BACMinor minorStatus,
                             BACWhere where, bool haveAnswer,
                             OsiSolverInterface *answerSolver) {
        bab_.majorStatus_ = majorStatus ;
        bab_.minorStatus_ = minorStatus ;
        bab_.where_ = where ;
        bab_.haveAnswer_ = haveAnswer ;
        bab_.answerSolver_ = answerSolver ;
    }

    /*! \brief Set the result of branch-and-cut search

      This version will extract the necessary information from the CbcModel
      object and set appropriate status based on the value passed for where.
    */
    void setBaBStatus(const CbcModel *model, BACWhere where,
                      bool haveAnswer = false,
                      OsiSolverInterface *answerSolver = 0) ;

    /*! \brief Translate CbcModel major status to #BACMajor

      See the #BACMajor enum for details.
    */
    BACMajor translateMajor(int status) ;

    /*!\brief Translate CbcModel minor status to #BACMinor

      See the #BACMinor enum for details.
    */
    BACMinor translateMinor(int status) ;

    /*!\brief Translate OsiSolverInterface status to #BACMinor

      See the #BACMinor enum for details. Optimal, infeasible, and unbounded
      get their own codes; everything else maps to BACmOther.
    */
    BACMinor translateMinor(const OsiSolverInterface *osi) ;

    /*! \brief Print the status block */

    void printBaBStatus() ;

//@}

    /*! \name Messages and statistics */
//@{

    /*! \brief Print a message

      Uses the current message handler and messages.
    */
    CoinMessageHandler &message(CbcGenMsgCode inID) ;

    /*! \brief Supply a new message handler.

      Replaces the current message handler. The current handler is destroyed
      if ourMsgHandler_ is true, and the call will set ourMsgHandler_ = true.
    */
    void passInMessageHandler(CoinMessageHandler *handler) ;

    /*! \brief Return a pointer to the message handler */
    inline CoinMessageHandler *messageHandler() const {
        return msgHandler_ ;
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
    void setMessages(CoinMessages::Language lang = CoinMessages::us_en) ;

    /*! \brief Set log level */
    inline void setLogLevel(int lvl) {
        logLvl_ = lvl ;
        if (msgHandler_) msgHandler_->setLogLevel(lvl) ;
    }

    /*! \brief Get log level */
    inline int logLevel() const {
        return (logLvl_) ;
    }

    /*! \brief When greater than 0, integer presolve gives more information and
           branch-and-cut provides statistics.
    */
    int printOpt_ ;

//@}

    /*! \name Parameter parsing and input/output. */
//@{
    /*! \brief cbc-generic version */

    std::string version_ ;

    /*! \brief Default directory prefix */

    std::string dfltDirectory_ ;

    /*! \brief Last MPS input file */

    std::string lastMpsIn_ ;

    /*! \brief Allow/disallow errors when importing a model */
    bool allowImportErrors_ ;

    /*! \brief Last solution output file */

    std::string lastSolnOut_ ;

    /*! \brief Solution printing mode

      Controls the amount of information printed when printing a solution.
      Coding is set by the keyword declarations for the printingOptions
      command.
    */
    int printMode_ ;

    /*! \brief Print mask

      Used to specify row/column names to be printed. Not implemented as of
      060920.
    */
    std::string printMask_ ;

    /*! \brief The parameter vector */

    CoinParamVec *paramVec_ ;

    /*! \brief Start and end of cbc-generic parameters in parameter vector */

    struct genParamsInfo_struct {
        int first_ ;
        int last_ ;
    } genParams_ ;

    /*! \brief Start and end of CbcModel parameters in parameter vector */

    struct cbcParamsInfo_struct {
        int first_ ;
        int last_ ;
    } cbcParams_ ;

    /*! \brief Start and end of OsiSolverInterface  parameters in parameter
           vector
    */

    struct osiParamsInfo_struct {
        int first_ ;
        int last_ ;
    } osiParams_ ;

    /*! \brief Verbosity level for help messages.

      Interpretation is bitwise:
      - (0): short help
      - (1): long help
      - (2): unused (for compatibility with cbc; indicates AMPL)
      - (3): show parameters with display = false.
    */

    int verbose_ ;

    /*! \brief Number of parameters processed */

    int paramsProcessed_ ;

    /*! \brief Record of parameters changed by user command */

    std::vector<bool> setByUser_ ;

    /*! \brief False if the user has made nontrivial modifications to the
           default control settings.

      Initially true. Specifying DJFIX, TIGHTENFACTOR, or any cut or heuristic
      parameter will set this to false.
    */
    bool defaultSettings_ ;

    /*! \brief Control debug file creation

      At the conclusion of branch-and-cut, dump the full solution in a binary
      format to debug.file in the current directory.  When set to
      "createAfterPre", the solution is dumped before integer presolve
      transforms are removed.  When set to "create", the solution is dumped
      after integer presolve transforms are backed out.
    */
    std::string debugCreate_ ;

    /*! \brief Last debug input file

      The file is expected to be in a binary format understood by
      activateRowCutDebugger.
    */

    std::string debugFile_ ;

    /*! \brief Array of primal variable values for debugging

      Used to provide a known optimal solution to activateRowCutDebugger().
    */

    struct debugSolInfo_struct {
        int numCols_ ;
        double *values_ ;
    } debugSol_ ;
//@}

    /* \name Timing */
//@{

    /*! \brief Total elapsed time for this run. */

    double totalTime_ ;

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

    CbcModel *model_ ;

    /*! \brief The current default LP solver

      This is a pointer to a reference copy. If you want the solver associated
      with #model_, ask for it directly.
    */

    OsiSolverInterface *dfltSolver_ ;

    /*! \brief True if we have a valid model loaded, false otherwise. */

    bool goodModel_ ;

    /*! \brief State of branch-and-cut

      Major and minor status codes, and a solver holding the answer, assuming
      we have a valid answer. See the documentation with the BACMajor,
      BACMinor, and BACWhere enums for the meaning of the codes.
    */

    struct babState_struct {
        BACMajor majorStatus_ ;
        BACMinor minorStatus_ ;
        BACWhere where_ ;
        bool haveAnswer_ ;
        OsiSolverInterface *answerSolver_ ;
    } bab_ ;

//@}

    /*! \name Various algorithm control variables and settings */
//@{

    /*! \brief Control use of reduced cost fixing prior to B&C

      This heuristic fixes variables whose reduced cost for the root
      relaxtion exceeds the specified threshold. This is purely a heuristic,
      performed before there's any incumbent solution. It may well fix variables
      at the wrong bound!
    */

    struct djFixCtl_struct {
        bool action_ ;
        double threshold_ ;
    } djFix_ ;

    /*! \brief Control the assignment of branching priorities to integer
           variables.
    */
    BPControl priorityAction_ ;

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
    struct chooseStrongCtl_struct {
        int numBeforeTrust_ ;
        int numStrong_ ;
        int shadowPriceMode_ ;
    } chooseStrong_ ;
//@}

private:

    /*! \name Cut Generator and Heuristic Control
        \brief Usage control and prototypes for cut generators and heuristics.
    */
//@{

    /*! \brief Control integer preprocessing. */

    IPPControl preProcess_ ;

    /*! \brief Control cut generator activity

      Generators that are active in the tree will be activated when
      (depth) mod (cutDepth) == 0.
    */

    int cutDepth_ ;

    /*! \brief Control variable and prototype for probing cut generator */
    struct probingCtl_struct {
        CGControl action_ ;
        CglProbing *proto_ ;
        bool usingObjective_ ;
        int maxPass_ ;
        int maxPassRoot_ ;
        int maxProbe_ ;
        int maxProbeRoot_ ;
        int maxLook_ ;
        int maxLookRoot_ ;
        int maxElements_ ;
        int rowCuts_ ;
    } probing_ ;

    /*! \brief Control variable and prototype for clique cut generator */
    struct cliqueCtl_struct {
        CGControl action_ ;
        CglClique *proto_ ;
        bool starCliqueReport_ ;
        bool rowCliqueReport_ ;
        double minViolation_ ;
    } clique_ ;

    /*! \brief Control variable and prototype for flow cover cut generator */
    struct flowCtl_struct {
        CGControl action_ ;
        CglFlowCover *proto_ ;
    } flow_ ;

    /*! \brief Control variable and prototype for Gomory cut generator */
    struct gomoryCtl_struct {
        CGControl action_ ;
        CglGomory *proto_ ;
        int limit_ ;
        int limitAtRoot_ ;
    } gomory_ ;

    /*   \brief Control variable and prototype for lift-and-project cut
    	     generator
       struct landpCtl_struct
       { CGControl action_ ;
         CglLandP *proto_ ; } landp_ ;
    */

    /*! \brief Control variable and prototype for knapsack cover cut generator */
    struct knapsackCtl_struct {
        CGControl action_ ;
        CglKnapsackCover *proto_ ;
    } knapsack_ ;

    /*! \brief Control variable and prototype for MIR cut generator */
    struct mirCtl_struct {
        CGControl action_ ;
        CglMixedIntegerRounding2 *proto_ ;
    } mir_ ;

    /*! \brief Control variable and prototype for odd hole cut generator */
    struct oddHoleCtl_struct {
        CGControl action_ ;
        CglOddHole *proto_ ;
    } oddHole_ ;

    /*! \brief Control variable and prototype for reduce-and-split
           cut generator
    */
    struct redSplitCtl_struct {
        CGControl action_ ;
        CglRedSplit *proto_ ;
    } redSplit_ ;

    /*! \brief Control variable and prototype for Two-MIR cut generator */
    struct twomirCtl_struct {
        CGControl action_ ;
        CglTwomir *proto_ ;
        int maxElements_ ;
    } twomir_ ;

    /*! \brief Control variable and prototype for feasibility pump heuristic */
    struct fpumpCtl_struct {
        CGControl action_ ;
        CbcHeuristicFPump *proto_ ;
        int iters_ ;
    } fpump_ ;

    /*! \brief Control variable and prototype for combine heuristic */
    struct combineCtl_struct {
        CGControl action_ ;
        CbcHeuristicLocal *proto_ ;
        int trySwap_ ;
    } combine_ ;

    /*! \brief Control variable and prototype for greedy cover heuristic */
    struct greedyCoverCtl_struct {
        CGControl action_ ;
        CbcHeuristicGreedyCover *proto_ ;
    } greedyCover_ ;

    /*! \brief Control variable and prototype for greedy equality heuristic */
    struct greedyEqualityCtl_struct {
        CGControl action_ ;
        CbcHeuristicGreedyEquality *proto_ ;
    } greedyEquality_ ;

    /*! \brief Control variable and prototype for simple rounding heuristic */
    struct roundingCtl_struct {
        CGControl action_ ;
        CbcRounding *proto_ ;
    } rounding_ ;


    /*! \brief Control variables for local tree

      This is a bit different --- getTreeLocal() takes a CbcModel as a parameter
      and installs a local tree object. But we can keep the parameters here and
      hide the details. Consult CbcTreeLocal.hpp for details.
    */
    struct localTreeCtl_struct {
        CGControl action_ ;
        CbcTreeLocal *proto_ ;
        double *soln_ ;
        int range_ ;
        int typeCuts_ ;
        int maxDiverge_ ;
        int timeLimit_ ;
        int nodeLimit_ ;
        bool refine_ ;
    } localTree_ ;

//@}

    /*! \name Messages and statistics (private)
        \brief Data and objects related to messages and statistics that should be
    	   protected from direct manipulation.
    */
//@{

    /*! \brief Message handler. */
    CoinMessageHandler *msgHandler_ ;

    /*! \brief Ownership of message handler.

      If true, the control block owns the message handler and it will be destroyed
      with the control block. If false, the client is responsible for the message
      handler.
    */
    bool ourMsgHandler_ ;

    /*! \brief The current language */
    CoinMessages::Language cur_lang_ ;

    /*! \brief The current set of messages. */
    CoinMessages *msgs_ ;

    /*! \brief The current log level */
    int logLvl_ ;

//@}

} ;


#endif

