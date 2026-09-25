// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcHeuristic_H
#define CbcHeuristic_H

#include <string>
#include <vector>
#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"
#include "CoinHelperFunctions.hpp"
#include "OsiBranchingObject.hpp"
#include "CbcConfig.h"
#include "CbcModel.hpp"

class OsiSolverInterface;

class CbcModel;
#include "OsiClpSolverInterface.hpp"
//#############################################################################

class CbcHeuristicNodeList;
class CbcBranchingObject;

/** A class describing the branching decisions that were made to get
    to the node where a heuristic was invoked from */

class CBCLIB_EXPORT CbcHeuristicNode {
private:
  void gutsOfConstructor(CbcModel &model);
  CbcHeuristicNode();
  CbcHeuristicNode &operator=(const CbcHeuristicNode &);

private:
  /// The number of branching decisions made
  int numObjects_;
  /** The indices of the branching objects. Note: an index may be
        listed multiple times. E.g., a general integer variable that has
        been branched on multiple times. */
  CbcBranchingObject **brObj_;

public:
  CbcHeuristicNode(CbcModel &model);

  CbcHeuristicNode(const CbcHeuristicNode &rhs);
  ~CbcHeuristicNode();
  double distance(const CbcHeuristicNode *node) const;
  double minDistance(const CbcHeuristicNodeList &nodeList) const;
  bool minDistanceIsSmall(const CbcHeuristicNodeList &nodeList,
    const double threshold) const;
  double avgDistance(const CbcHeuristicNodeList &nodeList) const;
};

class CBCLIB_EXPORT CbcHeuristicNodeList {
private:
  void gutsOfDelete();
  void gutsOfCopy(const CbcHeuristicNodeList &rhs);

private:
  std::vector< CbcHeuristicNode * > nodes_;

public:
  CbcHeuristicNodeList() {}
  CbcHeuristicNodeList(const CbcHeuristicNodeList &rhs);
  CbcHeuristicNodeList &operator=(const CbcHeuristicNodeList &rhs);
  ~CbcHeuristicNodeList();

  void append(CbcHeuristicNode *&node);
  void append(const CbcHeuristicNodeList &nodes);
  inline const CbcHeuristicNode *node(int i) const
  {
    return nodes_[i];
  }
  inline int size() const
  {
    return static_cast< int >(nodes_.size());
  }
};

//#############################################################################

/** Classification of heuristics by solution requirements.
    - CONSTRUCTIVE: does not need an existing feasible solution (e.g. FPump, Diving, Rounding)
    - IMPROVEMENT: needs at least one feasible solution (e.g. RINS, VND, DINS, Proximity)
    - IMPROVEMENT_2: needs at least two feasible solutions (e.g. Combine, Crossover)
*/
enum class HeuristicCategory {
  CONSTRUCTIVE = 0,
  IMPROVEMENT = 1,
  IMPROVEMENT_2 = 2
};

/** Node-count/depth-based scheduling policy for periodic re-invocation of a
    heuristic deep in the tree. Several heuristics (notably
    CbcHeuristicRINS/CbcHeuristicVND) historically reimplemented their own
    ad-hoc node-count math for this (magic constants like "wait at least 12
    nodes", "force a run at node 50 or 100", "double the interval if 2x it
    has passed with no trigger", decay the interval by a fixed fraction
    every 10 tries if success rate is below 30% -- see
    CbcHeuristicRINS::solution()/CbcHeuristicVND::solution() pre-refactor).
    This enum gives that policy an explicit, named, reusable form:
      - Legacy: unchanged, heuristic-specific historical behavior (default;
        shouldRunBySchedule() below is a no-op observer in this mode --
        each heuristic's own pre-existing gating logic remains solely
        responsible).
      - EveryKDepth: run only when the current tree depth is a multiple of
        K (i.e. at depths 0, K, 2K, 3K, ...).
      - EveryKNodes: run once every K nodes explored overall, independent
        of depth.
      - EveryKNodesNoImprove: run every K nodes since the last incumbent
        improvement, AND immediately again right after any improvement
        (which resets the K-node window). Suited to improvement
        heuristics like RINS/VND: a fresh incumbent is itself a natural
        trigger (there is a new target to compare the LP relaxation
        against), while a long stall (no improvement for K nodes) should
        prompt another attempt.
    See RINS-FIXTURES.md for the benchmarking behind these modes. */
enum class HeuristicScheduleMode {
  Legacy = 0,
  EveryKDepth = 1,
  EveryKNodes = 2,
  EveryKNodesNoImprove = 3
};

/** Heuristic base class */

class CBCLIB_EXPORT CbcHeuristic {
private:
  void gutsOfDelete() {}
  void gutsOfCopy(const CbcHeuristic &rhs);

public:
  // Default Constructor
  CbcHeuristic();

  // Constructor with model - assumed before cuts
  CbcHeuristic(CbcModel &model);

  // Copy constructor
  CbcHeuristic(const CbcHeuristic &);

  virtual ~CbcHeuristic();

  /// Clone
  virtual CbcHeuristic *clone() const = 0;

  /// Classification: constructive vs improvement heuristic
  virtual HeuristicCategory category() const { return HeuristicCategory::CONSTRUCTIVE; }

  /** Orchestration order used by CbcRootHeuristicSchedule (and any other
      scheduler that wants it) to decide, among heuristics eligible to run
      at the same moment, which ones to try first. Lower runs first; ties
      keep insertion order (stable sort). Purely a hint used by whichever
      scheduler groups multiple heuristics together -- it has no effect on
      when a heuristic runs on its own (see rootHooks()/shouldHeurRun()). */
  inline void setOrder(int value) { order_ = value; }
  inline int order() const { return order_; }

  /** Minimum number of incumbent solutions CBC must already have before
      this heuristic is allowed to start. 0 (the default for constructive
      heuristics, e.g. FeasibilityPump/Diving/FeasibilityJump) means "no
      incumbent needed yet". Improvement heuristics that require at least
      one/two solutions to build a neighbourhood (RINS/VND: 1, Combine2: 2)
      should set this explicitly; category() supplies a sensible default
      (see effectiveMinSolutionsToStart()) when it is left unset. */
  inline void setMinSolutionsToStart(int value) { minSolutionsToStart_ = value; }
  inline int minSolutionsToStart() const { return minSolutionsToStart_; }

  /** Maximum number of incumbent solutions CBC may already have for this
      heuristic to still be worth starting. -1 (default) means unlimited.
      Set to 0 for heuristics that should only be tried while no feasible
      solution has been found at all (e.g. an expensive constructive
      heuristic that stops being useful once cheaper ones already
      succeeded). */
  inline void setMaxSolutionsToStart(int value) { maxSolutionsToStart_ = value; }
  inline int maxSolutionsToStart() const { return maxSolutionsToStart_; }

  /** category()-derived default for minSolutionsToStart() when the user has
      not explicitly overridden it (minSolutionsToStart_ == -1 sentinel). */
  inline int effectiveMinSolutionsToStart() const
  {
    if (minSolutionsToStart_ >= 0)
      return minSolutionsToStart_;
    switch (category()) {
    case HeuristicCategory::IMPROVEMENT:
      return 1;
    case HeuristicCategory::IMPROVEMENT_2:
      return 2;
    case HeuristicCategory::CONSTRUCTIVE:
    default:
      return 0;
    }
  }

  /** Whether, given that CBC currently has numberSolutions incumbent(s),
      this heuristic is allowed to start at all (independent of *where* --
      see runsAtRootHook()). Centralises the
      minSolutionsToStart()/maxSolutionsToStart() gate so schedulers don't
      each reimplement it. */
  inline bool solutionCountAllowsStart(int numberSolutions) const
  {
    if (numberSolutions < effectiveMinSolutionsToStart())
      return false;
    if (maxSolutionsToStart_ >= 0 && numberSolutions > maxSolutionsToStart_)
      return false;
    return true;
  }

  /** Root-node "hook" moments (see CbcRootHeuristicSchedule) at which this
      heuristic is allowed to run, encoded as a string of single-letter
      codes:
        l - first LP solution (before pre-processing)
        L - pre-processed LP solution (default; matches historical
            behaviour of "root heuristics before cuts")
        c - an intermediate round of cut generation
        C - after cut generation is finished at the root
        P - after FeasibilityPump, seeded with its least-fractional
            attempt (currently implemented directly inside
            CbcHeuristicFPump's failure-recovery fallback -- see
            CbcHeuristicFPump::setFeasibilityJumpFallback() -- rather than
            through this generic hook, since it needs FPump's internal
            per-pass tracking; the code is reserved here for future use by
            a scheduler-level equivalent)
      Default is "LcC", i.e. runs before cuts, during intermediate cut
      rounds, and after root cuts, which is what every root heuristic has
      always done (the 'c' gate at CbcModel's intermediate-cut-round
      heuristics call was only added alongside this attribute, so the
      default had to include it to avoid changing behaviour). */
  inline void setRootHooks(const char *hooks) { rootHooks_ = hooks ? hooks : ""; }
  inline void setRootHooks(const std::string &hooks) { rootHooks_ = hooks; }
  inline const std::string &rootHooks() const { return rootHooks_; }
  inline bool runsAtRootHook(char hook) const
  {
    return rootHooks_.find(hook) != std::string::npos;
  }

  /// Assignment operator
  CbcHeuristic &operator=(const CbcHeuristic &rhs);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel *model);

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel *model) = 0;

  /** returns 0 if no solution, 1 if valid solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value
        This is called after cuts have been added - so can not add cuts
    */
  virtual int solution(double &objectiveValue,
    double *newSolution)
    = 0;

  /** returns 0 if no solution, 1 if valid solution, -1 if just
        returning an estimate of best possible solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value (only if nonzero code)
        This is called at same time as cut generators - so can add cuts
        Default is do nothing
    */
  virtual int solution2(double & /*objectiveValue*/,
    double * /*newSolution*/,
    OsiCuts & /*cs*/)
  {
    return 0;
  }

  /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
  virtual void validate() {}

  /** Sets "when" flag - 0 off, 1 at root, 2 other than root, 3 always.
        If 10 added then don't worry if validate says there are funny objects
        as user knows it will be fine
    */
  inline void setWhen(int value)
  {
    when_ = value;
  }
  /// Gets "when" flag - 0 off, 1 at root, 2 other than root, 3 always
  inline int when() const
  {
    return when_;
  }

  /// Sets number of nodes in subtree (default 200)
  inline void setNumberNodes(int value)
  {
    numberNodes_ = value;
  }
  /// Gets number of nodes in a subtree (default 200)
  inline int numberNodes() const
  {
    return numberNodes_;
  }
  /** Switches (does not apply equally to all heuristics)
        1 bit - stop once allowable gap on objective reached
        2 bit - always do given number of passes
        4 bit - weaken cutoff by 5% every 50 passes?
        8 bit - if has cutoff and suminf bobbling for 20 passes then
                first try halving distance to best possible then
                try keep halving distance to known cutoff
        16 bit - needs new solution to run
        1024 bit - stop all heuristics on max time
    */
  inline void setSwitches(int value)
  {
    switches_ = value;
  }
  /** Switches (does not apply equally to all heuristics)
        1 bit - stop once allowable gap on objective reached
        2 bit - always do given number of passes
        4 bit - weaken cutoff by 5% every 50 passes?
        8 bit - if has cutoff and suminf bobbling for 20 passes then
                first try halving distance to best possible then
                try keep halving distance to known cutoff
        16 bit - needs new solution to run
        1024 bit - stop all heuristics on max time
	65536 bit and above used for temporary communication
    */
  inline int switches() const
  {
    return switches_;
  }
  /// Whether to exit at once on gap
  bool exitNow(double bestObjective) const;
  /// Sets feasibility pump options (-1 is off)
  inline void setFeasibilityPumpOptions(int value)
  {
    feasibilityPumpOptions_ = value;
  }
  /// Gets feasibility pump options (-1 is off)
  inline int feasibilityPumpOptions() const
  {
    return feasibilityPumpOptions_;
  }
  /// Just set model - do not do anything else
  inline void setModelOnly(CbcModel *model)
  {
    model_ = model;
  }

  /// Sets fraction of new(rows+columns)/old(rows+columns) before doing small branch and bound (default 1.0)
  inline void setFractionSmall(double value)
  {
    fractionSmall_ = value;
  }
  /// Gets fraction of new(rows+columns)/old(rows+columns) before doing small branch and bound (default 1.0)
  inline double fractionSmall() const
  {
    return fractionSmall_;
  }
  /// Get how many solutions the heuristic thought it got
  inline int numberSolutionsFound() const
  {
    return numberSolutionsFound_;
  }
  /// Increment how many solutions the heuristic thought it got
  inline void incrementNumberSolutionsFound()
  {
    numberSolutionsFound_++;
  }

  /** Do mini branch and bound - return
        0 not finished - no solution
        1 not finished - solution
        2 finished - no solution
        3 finished - solution
        (could add global cut if finished)
        -1 returned on size
        -2 time or user event
    */
  int smallBranchAndBound(OsiSolverInterface *solver, int numberNodes,
    double *newSolution, double &newSolutionValue,
    double cutoff, std::string name) const;
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *) {}
  /// Create C++ lines to get to current state - does work for base class
  void generateCpp(FILE *fp, const char *heuristic);
  /// Returns true if can deal with "odd" problems e.g. sos type 2
  virtual bool canDealWithOdd() const
  {
    return false;
  }
  /// return name of heuristic
  inline const char *heuristicName() const
  {
    return heuristicName_.c_str();
  }
  /// set name of heuristic
  inline void setHeuristicName(const char *name)
  {
    heuristicName_ = name;
  }
  /// Set random number generator seed
  void setSeed(int value);
  /// Get random number generator seed
  int getSeed() const;
  /// Sets decay factor (for howOften) on failure
  inline void setDecayFactor(double value)
  {
    decayFactor_ = value;
  }
  /// Set input solution
  void setInputSolution(const double *solution, double objValue);
  /* Runs if bit set
        0 - before cuts at root node (or from doHeuristics)
        1 - during cuts at root
        2 - after root node cuts
        3 - after cuts at other nodes
        4 - during cuts at other nodes
            8 added if previous heuristic in loop found solution
     */
  inline void setWhereFrom(int value)
  {
    whereFrom_ = value;
  }
  inline int whereFrom() const
  {
    return whereFrom_;
  }
  /** Upto this depth we call the tree shallow and the heuristic can be called
        multiple times. That is, the test whether the current node is far from
        the others where the jeuristic was invoked will not be done, only the
        frequency will be tested. After that depth the heuristic will can be
        invoked only once per node, right before branching. That's when it'll be
        tested whether the heur should run at all. */
  inline void setShallowDepth(int value)
  {
    shallowDepth_ = value;
  }
  /** How often to invoke the heuristics in the shallow part of the tree */
  inline void setHowOftenShallow(int value)
  {
    howOftenShallow_ = value;
  }
  /** How "far" should this node be from every other where the heuristic was
        run in order to allow the heuristic to run in this node, too. Currently
        this is tested, but we may switch to avgDistanceToRun_ in the future. */
  inline void setMinDistanceToRun(int value)
  {
    minDistanceToRun_ = value;
  }
  /** Select a HeuristicScheduleMode and its K parameter. Only takes
      effect for heuristics that consult shouldRunBySchedule() from their
      own gating logic (currently CbcHeuristicRINS and CbcHeuristicVND);
      other heuristics ignore it. Default is Legacy (no behavior change). */
  inline void setScheduleMode(HeuristicScheduleMode mode, int k = 1)
  {
    scheduleMode_ = mode;
    scheduleK_ = k > 0 ? k : 1;
  }
  inline HeuristicScheduleMode scheduleMode() const
  {
    return scheduleMode_;
  }
  inline int scheduleK() const
  {
    return scheduleK_;
  }
  /** Returns true if, per scheduleMode()/scheduleK(), this is a node at
      which the heuristic should attempt to run. Has side effects
      (advances internal node/solution-count bookkeeping), so a heuristic
      should call it at most once per candidate invocation. Always
      returns true when scheduleMode() is Legacy -- in that mode it is a
      pure no-op observer and gates nothing; the heuristic's own
      pre-existing logic remains solely responsible. */
  bool shouldRunBySchedule();

  /** Check whether the heuristic should run at all
        0 - before cuts at root node (or from doHeuristics)
        1 - during cuts at root
        2 - after root node cuts
        3 - after cuts at other nodes
        4 - during cuts at other nodes
            8 added if previous heuristic in loop found solution
    */
  virtual bool shouldHeurRun(int whereFrom);
  /** Check whether the heuristic should run this time */
  bool shouldHeurRun_randomChoice();
  void debugNodes();
  void printDistanceToNodes();
  /// how many times the heuristic has actually run
  inline int numRuns() const
  {
    return numRuns_;
  }

  /// How many times the heuristic could run
  inline int numCouldRun() const
  {
    return numCouldRun_;
  }
  /// Return objective function value with sign corrected
  inline double trueObjValue(double value) const
  {
    return (model_->moreSpecialOptions2()&67108864)==0 ? value : -value;
  }
  /// Is it integer for heuristics?
  inline bool isHeuristicInteger(const OsiSolverInterface *solver, int iColumn) const
  {
    const OsiClpSolverInterface *clpSolver = getConstClpSolver(solver);
#ifndef CBC_OTHER_SOLVER
    return clpSolver->isHeuristicInteger(iColumn);
#else
    if (clpSolver)
      return clpSolver->isHeuristicInteger(iColumn);
    else
      return solver->isInteger(iColumn);
#endif
  }
  /*! \brief Clone, but ...

      If type is
	- 0 clone the solver for the model,
	- 1 clone the continuous solver for the model
        - Add 2 to say without integer variables which are at low priority
        - Add 4 to say quite likely infeasible so give up easily (clp only).
    */
  OsiSolverInterface *cloneBut(int type);

protected:
  /// Model
  CbcModel *model_;
  /// When flag - 0 off, 1 at root, 2 other than root, 3 always
  int when_;
  /// Number of nodes in any sub tree
  int numberNodes_;
  /** Feasibility pump options , -1 is off
	>=0 for feasibility pump itself
        -2 quick proximity search
        -3 longer proximity search
    */
  int feasibilityPumpOptions_;
  /// Fraction of new(rows+columns)/old(rows+columns) before doing small branch and bound
  mutable double fractionSmall_;
  /// Thread specific random number generator
  CoinThreadRandom randomNumberGenerator_;
  /// Name for printing
  std::string heuristicName_;

  /// See setOrder()
  int order_;
  /// See setMinSolutionsToStart(); -1 sentinel means "use category() default"
  int minSolutionsToStart_;
  /// See setMaxSolutionsToStart(); -1 means unlimited
  int maxSolutionsToStart_;
  /// See setRootHooks()
  std::string rootHooks_;

  /// How often to do (code can change)
  mutable int howOften_;
  /// How much to increase how often
  double decayFactor_;
  /** Switches (does not apply equally to all heuristics)
        1 bit - stop once allowable gap on objective reached
        2 bit - always do given number of passes
        4 bit - weaken cutoff by 5% every 50 passes?
        8 bit - if has cutoff and suminf bobbling for 20 passes then
                first try halving distance to best possible then
                try keep halving distance to known cutoff
        16 bit - needs new solution to run
        1024 bit - stop all heuristics on max time
    */
  mutable int switches_;
  /* Runs if bit set
        0 - before cuts at root node (or from doHeuristics)
        1 - during cuts at root
        2 - after root node cuts
        3 - after cuts at other nodes
        4 - during cuts at other nodes
            8 added if previous heuristic in loop found solution
     */
  int whereFrom_;
  /** Upto this depth we call the tree shallow and the heuristic can be called
        multiple times. That is, the test whether the current node is far from
        the others where the jeuristic was invoked will not be done, only the
        frequency will be tested. After that depth the heuristic will can be
        invoked only once per node, right before branching. That's when it'll be
        tested whether the heur should run at all. */
  int shallowDepth_;
  /** How often to invoke the heuristics in the shallow part of the tree */
  int howOftenShallow_;
  /** How many invocations happened within the same node when in a shallow
        part of the tree. */
  int numInvocationsInShallow_;
  /** How many invocations happened when in the deep part of the tree. For
        every node we count only one invocation. */
  int numInvocationsInDeep_;
  /** After how many deep invocations was the heuristic run last time */
  int lastRunDeep_;
  /// how many times the heuristic has actually run
  int numRuns_;
  /** How "far" should this node be from every other where the heuristic was
        run in order to allow the heuristic to run in this node, too. Currently
        this is tested, but we may switch to avgDistanceToRun_ in the future. */
  int minDistanceToRun_;

  /// See setScheduleMode()
  HeuristicScheduleMode scheduleMode_;
  /// See setScheduleMode()
  int scheduleK_;
  /// Node count at which this heuristic last ran under the schedule
  int scheduleLastRunNode_;
  /// model_->getSolutionCount() observed at the last shouldRunBySchedule()
  /// call, used by EveryKNodesNoImprove to detect a fresh incumbent
  int scheduleLastSolutionCount_;

  /// The description of the nodes where this heuristic has been applied
  CbcHeuristicNodeList runNodes_;

  /// How many times the heuristic could run
  int numCouldRun_;

  /// How many solutions the heuristic thought it got
  int numberSolutionsFound_;

  /// How many nodes the heuristic did this go
  mutable int numberNodesDone_;

  // Input solution - so can be used as seed
  double *inputSolution_;

#ifdef JJF_ZERO
  /// Lower bounds of last node where the heuristic found a solution
  double *lowerBoundLastNode_;
  /// Upper bounds of last node where the heuristic found a solution
  double *upperBoundLastNode_;
#endif
};
/** Rounding class
 */

class CBCLIB_EXPORT CbcRounding : public CbcHeuristic {
public:
  // Default Constructor
  CbcRounding();

  // Constructor with model - assumed before cuts
  CbcRounding(CbcModel &model);

  // Copy constructor
  CbcRounding(const CbcRounding &);

  // Destructor
  ~CbcRounding();

  /// Assignment operator
  CbcRounding &operator=(const CbcRounding &rhs);

  /// Clone
  virtual CbcHeuristic *clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp);

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel *model);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel *model);

  using CbcHeuristic::solution;
  /** returns 0 if no solution, 1 if valid solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value (only if good)
        This is called after cuts have been added - so can not add cuts
    */
  virtual int solution(double &objectiveValue,
    double *newSolution);
  /** returns 0 if no solution, 1 if valid solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value (only if good)
        This is called after cuts have been added - so can not add cuts
        Use solutionValue rather than solvers one
    */
  virtual int solution(double &objectiveValue,
    double *newSolution,
    double solutionValue);
  /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
  virtual void validate();

  /// Set seed
  void setSeed(int value)
  {
    seed_ = value;
  }
  /** Check whether the heuristic should run at all
        0 - before cuts at root node (or from doHeuristics)
        1 - during cuts at root
        2 - after root node cuts
        3 - after cuts at other nodes
        4 - during cuts at other nodes
            8 added if previous heuristic in loop found solution
    */
  virtual bool shouldHeurRun(int whereFrom);

protected:
  // Data

  // Original matrix by column
  CoinPackedMatrix matrix_;

  // Original matrix by
  CoinPackedMatrix matrixByRow_;

  // Down locks
  unsigned short *down_;

  // Up locks
  unsigned short *up_;

  // Equality locks
  unsigned short *equal_;

  // Seed for random stuff
  int seed_;
};

/** Partial solution class
    If user knows a partial solution this tries to get an integer solution
    it uses hotstart information
 */

class CBCLIB_EXPORT CbcHeuristicPartial : public CbcHeuristic {
public:
  // Default Constructor
  CbcHeuristicPartial();

  /** Constructor with model - assumed before cuts
        Fixes all variables with priority <= given
        and does given number of nodes
    */
  CbcHeuristicPartial(CbcModel &model, int fixPriority = 10000, int numberNodes = 200);

  // Copy constructor
  CbcHeuristicPartial(const CbcHeuristicPartial &);

  // Destructor
  ~CbcHeuristicPartial();

  /// Assignment operator
  CbcHeuristicPartial &operator=(const CbcHeuristicPartial &rhs);

  /// Clone
  virtual CbcHeuristic *clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp);

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel *model);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel *model);

  using CbcHeuristic::solution;
  /** returns 0 if no solution, 1 if valid solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value (only if good)
        This is called after cuts have been added - so can not add cuts
    */
  virtual int solution(double &objectiveValue,
    double *newSolution);
  /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
  virtual void validate();

  /// Set priority level
  void setFixPriority(int value)
  {
    fixPriority_ = value;
  }

  /** Check whether the heuristic should run at all */
  virtual bool shouldHeurRun(int whereFrom);

protected:
  // Data

  // All variables with abs priority <= this will be fixed
  int fixPriority_;
};

/** heuristic - just picks up any good solution
    found by solver - see OsiBabSolver
 */

class CBCLIB_EXPORT CbcSerendipity : public CbcHeuristic {
public:
  // Default Constructor
  CbcSerendipity();

  /* Constructor with model
    */
  CbcSerendipity(CbcModel &model);

  // Copy constructor
  CbcSerendipity(const CbcSerendipity &);

  // Destructor
  ~CbcSerendipity();

  /// Assignment operator
  CbcSerendipity &operator=(const CbcSerendipity &rhs);

  /// Clone
  virtual CbcHeuristic *clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp);

  /// update model
  virtual void setModel(CbcModel *model);

  using CbcHeuristic::solution;
  /** returns 0 if no solution, 1 if valid solution.
        Sets solution values if good, sets objective value (only if good)
        We leave all variables which are at one at this node of the
        tree to that value and will
        initially set all others to zero.  We then sort all variables in order of their cost
        divided by the number of entries in rows which are not yet covered.  We randomize that
        value a bit so that ties will be broken in different ways on different runs of the heuristic.
        We then choose the best one and set it to one and repeat the exercise.

    */
  virtual int solution(double &objectiveValue,
    double *newSolution);
  /// Resets stuff if model changes
  virtual void resetModel(CbcModel *model);

protected:
};

/** Just One class - this chooses one at random
 */

class CBCLIB_EXPORT CbcHeuristicJustOne : public CbcHeuristic {
public:
  // Default Constructor
  CbcHeuristicJustOne();

  // Constructor with model - assumed before cuts
  CbcHeuristicJustOne(CbcModel &model);

  // Copy constructor
  CbcHeuristicJustOne(const CbcHeuristicJustOne &);

  // Destructor
  ~CbcHeuristicJustOne();

  /// Clone
  virtual CbcHeuristicJustOne *clone() const;

  /// Assignment operator
  CbcHeuristicJustOne &operator=(const CbcHeuristicJustOne &rhs);

  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp);

  /** returns 0 if no solution, 1 if valid solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value (only if good)
        This is called after cuts have been added - so can not add cuts
        This does Fractional Diving
    */
  virtual int solution(double &objectiveValue,
    double *newSolution);
  /// Resets stuff if model changes
  virtual void resetModel(CbcModel *model);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel *model);
  /// Selects the next variable to branch on
  /** Returns true if all the fractional variables can be trivially
        rounded. Returns false, if there is at least one fractional variable
        that is not trivially roundable. In this case, the bestColumn
        returned will not be trivially roundable.
        This is dummy as never called
    */
  virtual bool selectVariableToBranch(OsiSolverInterface * /*solver*/,
    const double * /*newSolution*/,
    int & /*bestColumn*/,
    int & /*bestRound*/)
  {
    return true;
  }
  /// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
  virtual void validate();
  /// Adds an heuristic with probability
  void addHeuristic(const CbcHeuristic *heuristic, double probability);
  /// Normalize probabilities
  void normalizeProbabilities();

protected:
  // Data

  // Probability of running a heuristic
  double *probabilities_;

  // Heuristics
  CbcHeuristic **heuristic_;

  // Number of heuristics
  int numberHeuristics_;
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
