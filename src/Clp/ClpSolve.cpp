// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// This file has higher level solve functions

#include "CoinPragma.hpp"
#include "ClpConfig.h"

#if defined(CLP_HAS_AMD) || defined(CLP_HAS_CHOLMOD)
#define UFL_BARRIER
#endif

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpHelperFunctions.hpp"
#include "CoinSort.hpp"
#include "ClpFactorization.hpp"
#include "ClpSimplex.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpSimplexDual.hpp"
#ifndef SLIM_CLP
#include "ClpQuadraticObjective.hpp"
#include "ClpInterior.hpp"
#include "ClpCholeskyDense.hpp"
#include "ClpCholeskyBase.hpp"
#include "ClpPlusMinusOneMatrix.hpp"
#include "ClpNetworkMatrix.hpp"
#endif
#include "ClpEventHandler.hpp"
#include "ClpLinearObjective.hpp"
#include "ClpSolve.hpp"
#include "ClpPackedMatrix.hpp"
#include "ClpMessage.hpp"
#include "CoinTime.hpp"
#include "CoinStructuredModel.hpp"
double zz_slack_value = 0.0;
#ifdef CLP_USEFUL_PRINTOUT
double debugDouble[10];
int debugInt[24];
#endif

#include "ClpPresolve.hpp"
#ifndef SLIM_CLP
#include "Idiot.hpp"
#if defined(UFL_BARRIER) && (defined(CLP_HAS_AMD) || defined(CLP_HAS_CHOLMOD))
#include "ClpCholeskyUfl.hpp"
#endif
#ifdef PARDISO_BARRIER
#include "ClpCholeskyPardiso.hpp"
#endif
#ifdef COIN_HAS_VOL
#include "VolVolume.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinMpsIO.hpp"
//#############################################################################

class lpHook : public VOL_user_hooks {
private:
  lpHook(const lpHook &);
  lpHook &operator=(const lpHook &);

private:
  /// Pointer to dense vector of structural variable upper bounds
  double *colupper_;
  /// Pointer to dense vector of structural variable lower bounds
  double *collower_;
  /// Pointer to dense vector of objective coefficients
  double *objcoeffs_;
  /// Pointer to dense vector of right hand sides
  double *rhs_;
  /// Pointer to dense vector of senses
  char *sense_;

  /// The problem matrix in a row ordered form
  CoinPackedMatrix rowMatrix_;
  /// The problem matrix in a column ordered form
  CoinPackedMatrix colMatrix_;

public:
  lpHook(double *clb, double *cub, double *obj,
    double *rhs, char *sense, const CoinPackedMatrix &mat);
  virtual ~lpHook();

public:
  // for all hooks: return value of -1 means that volume should quit
  /** compute reduced costs
         @param u (IN) the dual variables
         @param rc (OUT) the reduced cost with respect to the dual values
     */
  virtual int compute_rc(const VOL_dvector &u, VOL_dvector &rc);

  /** Solve the subproblem for the subgradient step.
         @param dual (IN) the dual variables
         @param rc (IN) the reduced cost with respect to the dual values
         @param lcost (OUT) the lagrangean cost with respect to the dual values
         @param x (OUT) the primal result of solving the subproblem
         @param v (OUT) b-Ax for the relaxed constraints
         @param pcost (OUT) the primal objective value of <code>x</code>
     */
  virtual int solve_subproblem(const VOL_dvector &dual, const VOL_dvector &rc,
    double &lcost, VOL_dvector &x, VOL_dvector &v,
    double &pcost);
  /** Starting from the primal vector x, run a heuristic to produce
         an integer solution
         @param x (IN) the primal vector
         @param heur_val (OUT) the value of the integer solution (return
         <code>DBL_MAX</code> here if no feas sol was found
     */
  virtual int heuristics(const VOL_problem &p,
    const VOL_dvector &x, double &heur_val)
  {
    return 0;
  }
};

//#############################################################################

lpHook::lpHook(double *clb, double *cub, double *obj,
  double *rhs, char *sense,
  const CoinPackedMatrix &mat)
{
  colupper_ = cub;
  collower_ = clb;
  objcoeffs_ = obj;
  rhs_ = rhs;
  sense_ = sense;
  assert(mat.isColOrdered());
  colMatrix_.copyOf(mat);
  rowMatrix_.reverseOrderedCopyOf(mat);
}

//-----------------------------------------------------------------------------

lpHook::~lpHook()
{
}

//#############################################################################

int lpHook::compute_rc(const VOL_dvector &u, VOL_dvector &rc)
{
  rowMatrix_.transposeTimes(u.v, rc.v);
  const int psize = rowMatrix_.getNumCols();

  for (int i = 0; i < psize; ++i)
    rc[i] = objcoeffs_[i] - rc[i];
  return 0;
}

//-----------------------------------------------------------------------------

int lpHook::solve_subproblem(const VOL_dvector &dual, const VOL_dvector &rc,
  double &lcost, VOL_dvector &x, VOL_dvector &v,
  double &pcost)
{
  int i;
  const int psize = x.size();
  const int dsize = v.size();

  // compute the lagrangean solution corresponding to the reduced costs
  for (i = 0; i < psize; ++i)
    x[i] = (rc[i] >= 0.0) ? collower_[i] : colupper_[i];

  // compute the lagrangean value (rhs*dual + primal*rc)
  lcost = 0;
  for (i = 0; i < dsize; ++i)
    lcost += rhs_[i] * dual[i];
  for (i = 0; i < psize; ++i)
    lcost += x[i] * rc[i];

  // compute the rhs - lhs
  colMatrix_.times(x.v, v.v);
  for (i = 0; i < dsize; ++i)
    v[i] = rhs_[i] - v[i];

  // compute the lagrangean primal objective
  pcost = 0;
  for (i = 0; i < psize; ++i)
    pcost += x[i] * objcoeffs_[i];

  return 0;
}

//#############################################################################
/** A quick inlined function to convert from lb/ub style constraint
    definition to sense/rhs/range style */
inline void
convertBoundToSense(const double lower, const double upper,
  char &sense, double &right,
  double &range)
{
  range = 0.0;
  if (lower > -1.0e20) {
    if (upper < 1.0e20) {
      right = upper;
      if (upper == lower) {
        sense = 'E';
      } else {
        sense = 'R';
        range = upper - lower;
      }
    } else {
      sense = 'G';
      right = lower;
    }
  } else {
    if (upper < 1.0e20) {
      sense = 'L';
      right = upper;
    } else {
      sense = 'N';
      right = 0.0;
    }
  }
}

static int
solveWithVolume(ClpSimplex *model, int numberPasses, int doIdiot)
{
  VOL_problem volprob;
  volprob.parm.gap_rel_precision = 0.00001;
  volprob.parm.maxsgriters = 3000;
  if (numberPasses > 3000) {
    volprob.parm.maxsgriters = numberPasses;
    volprob.parm.primal_abs_precision = 0.0;
    volprob.parm.minimum_rel_ascent = 0.00001;
  } else if (doIdiot > 0) {
    volprob.parm.maxsgriters = doIdiot;
  }
  if (model->logLevel() < 2)
    volprob.parm.printflag = 0;
  else
    volprob.parm.printflag = 3;
  const CoinPackedMatrix *mat = model->matrix();
  int psize = model->numberColumns();
  int dsize = model->numberRows();
  char *sense = new char[dsize];
  double *rhs = new double[dsize];

  // Set the lb/ub on the duals
  volprob.dsize = dsize;
  volprob.psize = psize;
  volprob.dual_lb.allocate(dsize);
  volprob.dual_ub.allocate(dsize);
  int i;
  const double *rowLower = model->rowLower();
  const double *rowUpper = model->rowUpper();
  for (i = 0; i < dsize; ++i) {
    double range;
    convertBoundToSense(rowLower[i], rowUpper[i],
      sense[i], rhs[i], range);
    switch (sense[i]) {
    case 'E':
      volprob.dual_lb[i] = -1.0e31;
      volprob.dual_ub[i] = 1.0e31;
      break;
    case 'L':
      volprob.dual_lb[i] = -1.0e31;
      volprob.dual_ub[i] = 0.0;
      break;
    case 'G':
      volprob.dual_lb[i] = 0.0;
      volprob.dual_ub[i] = 1.0e31;
      break;
    default:
      printf("Volume Algorithm can't work if there is a non ELG row\n");
      return 1;
    }
  }
  // Check bounds
  double *saveLower = model->columnLower();
  double *saveUpper = model->columnUpper();
  bool good = true;
  for (i = 0; i < psize; i++) {
    if (saveLower[i] < -1.0e20 || saveUpper[i] > 1.0e20) {
      good = false;
      break;
    }
  }
  if (!good) {
    saveLower = CoinCopyOfArray(model->columnLower(), psize);
    saveUpper = CoinCopyOfArray(model->columnUpper(), psize);
    for (i = 0; i < psize; i++) {
      if (saveLower[i] < -1.0e20)
        saveLower[i] = -1.0e20;
      if (saveUpper[i] > 1.0e20)
        saveUpper[i] = 1.0e20;
    }
  }
  lpHook myHook(saveLower, saveUpper,
    model->objective(),
    rhs, sense, *mat);

  volprob.solve(myHook, false /* no warmstart */);

  if (saveLower != model->columnLower()) {
    delete[] saveLower;
    delete[] saveUpper;
  }
  //------------- extract the solution ---------------------------

  //printf("Best lagrangean value: %f\n", volprob.value);

  double avg = 0;
  for (i = 0; i < dsize; ++i) {
    switch (sense[i]) {
    case 'E':
      avg += std::abs(volprob.viol[i]);
      break;
    case 'L':
      if (volprob.viol[i] < 0)
        avg += (-volprob.viol[i]);
      break;
    case 'G':
      if (volprob.viol[i] > 0)
        avg += volprob.viol[i];
      break;
    }
  }

  //printf("Average primal constraint violation: %f\n", avg/dsize);

  // volprob.dsol contains the dual solution (dual feasible)
  // volprob.psol contains the primal solution
  //              (NOT necessarily primal feasible)
  CoinMemcpyN(volprob.dsol.v, dsize, model->dualRowSolution());
  CoinMemcpyN(volprob.psol.v, psize, model->primalColumnSolution());
  return 0;
}
#endif
static ClpInterior *currentModel2 = NULL;
#endif
//#############################################################################
// Allow for interrupts
// But is this threadsafe ? (so switched off by option)

#include "CoinSignal.hpp"
static ClpSimplex *currentModel = NULL;

extern "C" {
static void
#if defined(_MSC_VER)
  __cdecl
#endif // _MSC_VER
  signal_handler(int /*whichSignal*/)
{
  if (currentModel != NULL)
    currentModel->setMaximumIterations(0); // stop at next iterations
#ifndef SLIM_CLP
  if (currentModel2 != NULL)
    currentModel2->setMaximumBarrierIterations(0); // stop at next iterations
#endif
  return;
}
}
#if ABC_INSTRUMENT > 1
int abcPricing[20];
int abcPricingDense[20];
static int trueNumberRows;
static int numberTypes;
#define MAX_TYPES 25
#define MAX_COUNT 20
#define MAX_FRACTION 101
static char *types[MAX_TYPES];
static double counts[MAX_TYPES][MAX_COUNT];
static double countsFraction[MAX_TYPES][MAX_FRACTION];
static double *currentCounts;
static double *currentCountsFraction;
static int currentType;
static double workMultiplier[MAX_TYPES];
static double work[MAX_TYPES];
static double currentWork;
static double otherWork[MAX_TYPES];
static int timesCalled[MAX_TYPES];
static int timesStarted[MAX_TYPES];
static int fractionDivider;
void instrument_initialize(int numberRows)
{
  trueNumberRows = numberRows;
  numberTypes = 0;
  memset(counts, 0, sizeof(counts));
  currentCounts = NULL;
  memset(countsFraction, 0, sizeof(countsFraction));
  currentCountsFraction = NULL;
  memset(workMultiplier, 0, sizeof(workMultiplier));
  memset(work, 0, sizeof(work));
  memset(otherWork, 0, sizeof(otherWork));
  memset(timesCalled, 0, sizeof(timesCalled));
  memset(timesStarted, 0, sizeof(timesStarted));
  currentType = -1;
  fractionDivider = (numberRows + MAX_FRACTION - 2) / (MAX_FRACTION - 1);
}
void instrument_start(const char *type, int numberRowsEtc)
{
  if (currentType >= 0)
    instrument_end();
  currentType = -1;
  currentWork = 0.0;
  for (int i = 0; i < numberTypes; i++) {
    if (!strcmp(types[i], type)) {
      currentType = i;
      break;
    }
  }
  if (currentType == -1) {
    assert(numberTypes < MAX_TYPES);
    currentType = numberTypes;
    types[numberTypes++] = strdup(type);
  }
  currentCounts = &counts[currentType][0];
  currentCountsFraction = &countsFraction[currentType][0];
  timesStarted[currentType]++;
  assert(trueNumberRows);
  workMultiplier[currentType] += static_cast< double >(numberRowsEtc) / static_cast< double >(trueNumberRows);
}
void instrument_add(int count)
{
  assert(currentType >= 0);
  currentWork += count;
  timesCalled[currentType]++;
  if (count < MAX_COUNT - 1)
    currentCounts[count]++;
  else
    currentCounts[MAX_COUNT - 1]++;
  assert(count / fractionDivider >= 0 && count / fractionDivider < MAX_FRACTION);
  currentCountsFraction[count / fractionDivider]++;
}
void instrument_do(const char *type, double count)
{
  int iType = -1;
  for (int i = 0; i < numberTypes; i++) {
    if (!strcmp(types[i], type)) {
      iType = i;
      break;
    }
  }
  if (iType == -1) {
    assert(numberTypes < MAX_TYPES);
    iType = numberTypes;
    types[numberTypes++] = strdup(type);
  }
  timesStarted[iType]++;
  otherWork[iType] += count;
}
void instrument_end()
{
  work[currentType] += currentWork;
  currentType = -1;
}
void instrument_end_and_adjust(double factor)
{
  work[currentType] += currentWork * factor;
  currentType = -1;
}
void instrument_print()
{
  for (int iType = 0; iType < numberTypes; iType++) {
    currentCounts = &counts[iType][0];
    currentCountsFraction = &countsFraction[iType][0];
    if (!otherWork[iType]) {
      printf("%s started %d times, used %d times, work %g (average length %.1f) multiplier %g\n",
        types[iType], timesStarted[iType], timesCalled[iType],
        work[iType], work[iType] / (timesCalled[iType] + 1.0e-100), workMultiplier[iType] / (timesStarted[iType] + 1.0e-100));
      int n = 0;
      for (int i = 0; i < MAX_COUNT - 1; i++) {
        if (currentCounts[i]) {
          if (n == 5) {
            n = 0;
            printf("\n");
          }
          n++;
          printf("(%d els,%.0f times) ", i, currentCounts[i]);
        }
      }
      if (currentCounts[MAX_COUNT - 1]) {
        if (n == 5) {
          n = 0;
          printf("\n");
        }
        n++;
        printf("(>=%d els,%.0f times) ", MAX_COUNT - 1, currentCounts[MAX_COUNT - 1]);
      }
      printf("\n");
      int largestFraction;
      int nBig = 0;
      for (largestFraction = MAX_FRACTION - 1; largestFraction >= 10; largestFraction--) {
        double count = currentCountsFraction[largestFraction];
        if (count && largestFraction > 10)
          nBig++;
        if (nBig > 4)
          break;
      }
      int chunk = (largestFraction + 5) / 10;
      int lo = 0;
      for (int iChunk = 0; iChunk < largestFraction; iChunk += chunk) {
        int hi = std::min(lo + chunk * fractionDivider, trueNumberRows);
        double sum = 0.0;
        for (int i = iChunk; i < std::min(iChunk + chunk, MAX_FRACTION); i++)
          sum += currentCountsFraction[i];
        if (sum)
          printf("(%d-%d %.0f) ", lo, hi, sum);
        lo = hi;
      }
      for (int i = lo / fractionDivider; i < MAX_FRACTION; i++) {
        if (currentCountsFraction[i])
          printf("(%d %.0f) ", i * fractionDivider, currentCountsFraction[i]);
      }
      printf("\n");
    } else {
      printf("%s started %d times, used %d times, work %g multiplier %g other work %g\n",
        types[iType], timesStarted[iType], timesCalled[iType],
        work[iType], workMultiplier[iType], otherWork[iType]);
    }
    free(types[iType]);
  }
}
#endif
//#if ABC_PARALLEL == 2
/** General solve algorithm which can do presolve
    special options (bits)
    1 - do not perturb
    2 - do not scale
    4 - use crash (default allslack in dual, idiot in primal)
    8 - all slack basis in primal
    16 - switch off interrupt handling
    32 - do not try and make plus minus one matrix
    64 - do not use sprint even if problem looks good
 */
int ClpSimplex::initialSolve(ClpSolve &options)
{
  presolveTime_ = 0.0;
  presolveRows_ = -1;
  presolveCols_ = -1;
  ClpSolve::SolveType method = options.getSolveType();
  //ClpSolve::SolveType originalMethod=method;
  ClpSolve::PresolveType presolve = options.getPresolveType();
  int saveMaxIterations = maximumIterations();
  int finalStatus = -1;
  int numberIterations = 0;
  // Select timing basis: wall time when a wall-clock limit is active (the
  // default in CBC), CPU time otherwise.  All interval and summary messages
  // produced by this function use the same basis so they stay consistent with
  // whatever limit is being enforced by hitMaximumIterations().
  const bool clpUseWallTime = (dblParam_[ClpMaxWallSeconds] >= 0.0);
  auto clpGetTime = [clpUseWallTime]() -> double {
    return clpUseWallTime ? CoinGetTimeOfDay() : CoinCpuTime();
  };
  double time1 = clpGetTime();
  double timeX = time1;
  double time2 = 0.0;
  ClpMatrixBase *saveMatrix = NULL;
  ClpObjective *savedObjective = NULL;
  int idiotOptions = 0;
  if (options.getSpecialOption(6))
    idiotOptions = options.getExtraInfo(6) * 32768;
#ifdef CLP_USEFUL_PRINTOUT
  debugInt[0] = numberRows();
  debugInt[1] = numberColumns();
  debugInt[2] = matrix()->getNumElements();
#endif
  if (!objective_ || !matrix_) {
    // totally empty
    handler_->message(CLP_EMPTY_PROBLEM, messages_)
      << 0
      << 0
      << 0
      << CoinMessageEol;
    return -1;
  } else if (!numberRows_ || !numberColumns_ || !getNumElements()) {
    presolve = ClpSolve::presolveOff;
  }
  if (objective_->type() >= 2 && optimizationDirection_ == 0) {
    // pretend linear
    savedObjective = objective_;
    // make up objective
    double *obj = new double[numberColumns_];
    for (int i = 0; i < numberColumns_; i++) {
      double l = fabs(columnLower_[i]);
      double u = fabs(columnUpper_[i]);
      obj[i] = 0.0;
      if (std::min(l, u) < 1.0e20) {
        if (l < u)
          obj[i] = 1.0 + randomNumberGenerator_.randomDouble() * 1.0e-2;
        else
          obj[i] = -1.0 - randomNumberGenerator_.randomDouble() * 1.0e-2;
      }
    }
    objective_ = new ClpLinearObjective(obj, numberColumns_);
    delete[] obj;
  }
  ClpSimplex *model2 = this;
  bool interrupt = (options.getSpecialOption(2) == 0);
  CoinSighandler_t saveSignal = static_cast< CoinSighandler_t >(0);
  if (interrupt) {
    currentModel = model2;
    // register signal handler
    saveSignal = signal(SIGINT, signal_handler);
  }
  // If no status array - set up basis
  if (!status_)
    allSlackBasis();
  ClpPresolve *pinfo = new ClpPresolve();
  pinfo->setSubstitution(options.substitution());
  int presolveOptions = options.presolveActions();
  bool presolveToFile = (presolveOptions & 0x40000000) != 0;
  presolveOptions &= ~0x40000000;
  if ((presolveOptions & 0xffffff) != 0)
    pinfo->setPresolveActions(presolveOptions);
  // switch off singletons to slacks
  //pinfo->setDoSingletonColumn(false); // done by bits
  int printOptions = options.getSpecialOption(5);
  if ((printOptions & 1) != 0)
    pinfo->statistics();
  double timePresolve = 0.0;
  double timeIdiot = 0.0;
  double timeCore = 0.0;
  eventHandler()->event(ClpEventHandler::presolveStart);
  int savePerturbation = perturbation_;
  int saveScaling = scalingFlag_;
#ifndef SLIM_CLP
#ifndef NO_RTTI
  if (dynamic_cast< ClpNetworkMatrix * >(matrix_)) {
    // network - switch off stuff
    presolve = ClpSolve::presolveOff;
  }
#else
  if (matrix_->type() == 11) {
    // network - switch off stuff
    presolve = ClpSolve::presolveOff;
  }
#endif
#endif
#ifndef CLPSOLVE_ACTIONS
#define CLPSOLVE_ACTIONS 2
#endif
#if CLPSOLVE_ACTIONS
  bool wasAutomatic = (method == ClpSolve::automatic);
#endif
  if (presolve != ClpSolve::presolveOff) {
    bool costedSlacks = false;
#if CLP_INHERIT_MODE > 1
    int numberPasses = 20;
#else
    int numberPasses = 5;
#endif
    if (presolve == ClpSolve::presolveNumber) {
      numberPasses = options.getPresolvePasses();
      presolve = ClpSolve::presolveOn;
    } else if (presolve == ClpSolve::presolveNumberCost) {
      numberPasses = options.getPresolvePasses();
      presolve = ClpSolve::presolveOn;
      costedSlacks = true;
      // switch on singletons to slacks
      pinfo->setDoSingletonColumn(true);
      // gub stuff for testing
      //pinfo->setDoGubrow(true);
    }
#ifndef CLP_NO_STD
    if (presolveToFile) {
      // PreSolve to file - not fully tested
      printf("Presolving to file - presolve.save\n");
      pinfo->presolvedModelToFile(*this, "presolve.save", dblParam_[ClpPresolveTolerance],
        false, numberPasses);
      model2 = this;
    } else {
#endif
      model2 = pinfo->presolvedModel(*this, dblParam_[ClpPresolveTolerance],
        false, numberPasses, true, costedSlacks);
#ifndef CLP_NO_STD
    }
#endif
    time2 = clpGetTime();
    timePresolve = time2 - timeX;
    presolveTime_ = timePresolve;
    if (model2 && model2 != this) {
      presolveRows_ = model2->numberRows();
      presolveCols_ = model2->numberColumns();
      // Store on model2 so event handler can report presolve reduction.
      // presolveRows_/Cols_ on model2 stores the presolved dimensions;
      // the event handler compares with origRows/origCols from state.
      model2->presolveTime_ = timePresolve;
      model2->presolveRows_ = model2->numberRows();
      model2->presolveCols_ = model2->numberColumns();
    }
    handler_->message(CLP_INTERVAL_TIMING, messages_)
      << "Presolve" << timePresolve << time2 - time1
      << CoinMessageEol;
    timeX = time2;
    if (!model2) {
      handler_->message(CLP_INFEASIBLE, messages_)
        << CoinMessageEol;
      model2 = this;
      eventHandler()->event(ClpEventHandler::presolveInfeasible);
      problemStatus_ = pinfo->presolveStatus();
      secondaryStatus_ = 11;
      if (options.infeasibleReturn() || (moreSpecialOptions_ & 1) != 0) {
        delete pinfo;
        return -1;
      }
      presolve = ClpSolve::presolveOff;
      // Barrier no good at infeasible problems
      if (method == ClpSolve::useBarrier ||
	  method == ClpSolve::useBarrierNoCross) {
	method = ClpSolve::usePrimal;
	handler_->message(CLP_GENERAL, messages_)
	  << "Looks infeasible - using Simplex rather than barrier"
	  << CoinMessageEol;
      }
    } else {
      //ClpModel::stopPermanentArrays();
      //setSpecialOptions(specialOptions()&~65536);
      // try setting tolerances up
#if CLPSOLVE_ACTIONS
      bool changeTolerances = wasAutomatic;
#if CLPSOLVE_ACTIONS > 1
      changeTolerances = true;
#endif
      if (changeTolerances && model2 != this) {
#define CLP_NEW_TOLERANCE 1.0e-6
        if (model2->primalTolerance() == 1.0e-7 && model2->dualTolerance() == 1.0e-7) {
          model2->setPrimalTolerance(CLP_NEW_TOLERANCE);
          model2->setDualTolerance(CLP_NEW_TOLERANCE);
        }
      }
#endif
      model2->eventHandler()->setSimplex(model2);
      int rcode = model2->eventHandler()->event(ClpEventHandler::presolveSize);
      // see if too big or small
      if (rcode == 2) {
        delete model2;
        delete pinfo;
        return -2;
      } else if (rcode == 3) {
        delete model2;
        delete pinfo;
        return -3;
      }
    }
    model2->setMoreSpecialOptions(model2->moreSpecialOptions() & (~1024));
    model2->eventHandler()->setSimplex(model2);
    // We may be better off using original (but if dual leave because of bounds)
    if (presolve != ClpSolve::presolveOff && numberRows_ < 1.01 * model2->numberRows_ && numberColumns_ < 1.01 * model2->numberColumns_
      && model2 != this) {
      if (method != ClpSolve::useDual || (numberRows_ == model2->numberRows_ && numberColumns_ == model2->numberColumns_)) {
        delete model2;
        model2 = this;
        presolve = ClpSolve::presolveOff;
      }
    }
  }
#ifdef CLP_USEFUL_PRINTOUT
  debugInt[3] = model2->numberRows();
  debugInt[4] = model2->numberColumns();
  debugInt[5] = model2->matrix()->getNumElements();
  // analyze
  {
    double time1 = CoinCpuTime();
    int numberColumns = model2->numberColumns();
    const double *columnLower = model2->columnLower();
    const double *columnUpper = model2->columnUpper();
    int numberRows = model2->numberRows();
    const double *rowLower = model2->rowLower();
    const double *rowUpper = model2->rowUpper();
    const double *objective = model2->objective();
    CoinPackedMatrix *matrix = model2->matrix();
    CoinBigIndex numberElements = matrix->getNumElements();
    const int *columnLength = matrix->getVectorLengths();
    //const CoinBigIndex * columnStart = matrix->getVectorStarts();
    const double *elementByColumn = matrix->getElements();
    const int *row = matrix->getIndices();
    int *rowCount = new int[numberRows];
    memset(rowCount, 0, numberRows * sizeof(int));
    int n = std::max(2 * numberRows, numberElements);
    n = std::max(2 * numberColumns, n);
    double *check = new double[n];
    memcpy(check, elementByColumn, numberElements * sizeof(double));
    for (int i = 0; i < numberElements; i++) {
      check[i] = fabs(check[i]);
      rowCount[row[i]]++;
    }
    int largestIndex = 0;
    for (int i = 0; i < numberColumns; i++) {
      largestIndex = std::max(largestIndex, columnLength[i]);
    }
    debugInt[12] = largestIndex;
    largestIndex = 0;
    for (int i = 0; i < numberRows; i++) {
      largestIndex = std::max(largestIndex, rowCount[i]);
    }
    n = numberElements;
    delete[] rowCount;
    debugInt[11] = largestIndex;
    std::sort(check, check + n);
    debugDouble[4] = check[0];
    debugDouble[5] = check[n - 1];
    int nAtOne = 0;
    int nDifferent = 0;
    double last = -COIN_DBL_MAX;
    for (int i = 0; i < n; i++) {
      if (fabs(last - check[i]) > 1.0e-12) {
        nDifferent++;
        last = check[i];
      }
      if (check[i] == 1.0)
        nAtOne++;
    }
    debugInt[10] = nDifferent;
    debugInt[15] = nAtOne;
    int nInf = 0;
    int nZero = 0;
    n = 0;
    for (int i = 0; i < numberRows; i++) {
      double value = fabs(rowLower[i]);
      if (!value)
        nZero++;
      else if (value != COIN_DBL_MAX)
        check[n++] = value;
      else
        nInf++;
    }
    for (int i = 0; i < numberRows; i++) {
      double value = fabs(rowUpper[i]);
      if (!value)
        nZero++;
      else if (value != COIN_DBL_MAX)
        check[n++] = value;
      else
        nInf++;
    }
    debugInt[16] = nInf;
    debugInt[20] = nZero;
    if (n) {
      std::sort(check, check + n);
      debugDouble[0] = check[0];
      debugDouble[1] = check[n - 1];
    } else {
      debugDouble[0] = 0.0;
      debugDouble[1] = 0.0;
    }
    nAtOne = 0;
    nDifferent = 0;
    last = -COIN_DBL_MAX;
    for (int i = 0; i < n; i++) {
      if (fabs(last - check[i]) > 1.0e-12) {
        nDifferent++;
        last = check[i];
      }
      if (check[i] == 1.0)
        nAtOne++;
    }
    debugInt[8] = nDifferent;
    debugInt[13] = nAtOne;
    nZero = 0;
    n = 0;
    for (int i = 0; i < numberColumns; i++) {
      double value = fabs(objective[i]);
      if (value)
        check[n++] = value;
      else
        nZero++;
    }
    debugInt[21] = nZero;
    if (n) {
      std::sort(check, check + n);
      debugDouble[2] = check[0];
      debugDouble[3] = check[n - 1];
    } else {
      debugDouble[2] = 0.0;
      debugDouble[3] = 0.0;
    }
    nAtOne = 0;
    nDifferent = 0;
    last = -COIN_DBL_MAX;
    for (int i = 0; i < n; i++) {
      if (fabs(last - check[i]) > 1.0e-12) {
        nDifferent++;
        last = check[i];
      }
      if (check[i] == 1.0)
        nAtOne++;
    }
    debugInt[9] = nDifferent;
    debugInt[14] = nAtOne;
    nInf = 0;
    nZero = 0;
    n = 0;
    for (int i = 0; i < numberColumns; i++) {
      double value = fabs(columnLower[i]);
      if (!value)
        nZero++;
      else if (value != COIN_DBL_MAX)
        check[n++] = value;
      else
        nInf++;
    }
    for (int i = 0; i < numberColumns; i++) {
      double value = fabs(columnUpper[i]);
      if (!value)
        nZero++;
      else if (value != COIN_DBL_MAX)
        check[n++] = value;
      else
        nInf++;
    }
    debugInt[17] = nInf;
    double smallestColBound;
    double largestColBound;
    if (n) {
      std::sort(check, check + n);
      smallestColBound = check[0];
      largestColBound = check[n - 1];
    } else {
      smallestColBound = 0.0;
      largestColBound = 0.0;
    }
    nAtOne = 0;
    nDifferent = 0;
    last = -COIN_DBL_MAX;
    for (int i = 0; i < n; i++) {
      if (fabs(last - check[i]) > 1.0e-12) {
        nDifferent++;
        last = check[i];
      }
      if (check[i] == 1.0)
        nAtOne++;
    }
    //debugInt[8]=nDifferent;
    //debugInt[13]=nAtOne;
    printf("BENCHMARK_STATS rhs %d different - %g -> %g (%d at one, %d infinite, %d zero)\n",
      debugInt[8], debugDouble[0], debugDouble[1], debugInt[13], debugInt[16], debugInt[20]);
    printf("BENCHMARK_STATS col %d different - %g -> %g (%d at one, %d infinite, %d zero)\n",
      nDifferent, smallestColBound, largestColBound, nAtOne, nInf, nZero);
    printf("BENCHMARK_STATS els %d different - %g -> %g (%d at one) - longest r,c %d,%d\n",
      debugInt[10], debugDouble[4], debugDouble[5], debugInt[15],
      debugInt[11], debugInt[12]);
    printf("BENCHMARK_STATS obj %d different - %g -> %g (%d at one, %d zero) - time %g\n",
      debugInt[9], debugDouble[2], debugDouble[3], debugInt[14], debugInt[21],
      CoinCpuTime() - time1);
    delete[] check;
  }
#endif
  if (interrupt)
    currentModel = model2;
  int saveMoreOptions = moreSpecialOptions_;
  // For below >0 overrides
  // 0 means no, -1 means maybe
  int doIdiot = 0;
  int doCrash = 0;
  int doSprint = 0;
  int doSlp = 0;
  int primalStartup = 1;
  model2->eventHandler()->event(ClpEventHandler::presolveBeforeSolve);
#if CLP_POOL_MATRIX
  if (vectorMode() >= 10) {
    ClpPoolMatrix *poolMatrix = new ClpPoolMatrix(*model2->matrix());
    char output[80];
    int numberDifferent = poolMatrix->getNumDifferentElements();
    if (numberDifferent > 0) {
      sprintf(output, "Pool matrix has %d different values",
        numberDifferent);
      model2->replaceMatrix(poolMatrix, true);
    } else {
      delete poolMatrix;
      sprintf(output, "Pool matrix has more than %d different values - no good",
        -numberDifferent);
    }
    handler_->message(CLP_GENERAL, messages_) << output
                                              << CoinMessageEol;
  }
#endif
  int tryItSave = 0;
#if CLPSOLVE_ACTIONS
  if (method == ClpSolve::automatic)
    model2->moreSpecialOptions_ |= 8192; // stop switch over
#endif
  // switch to primal from automatic if just one cost entry
  if (method == ClpSolve::automatic && model2->numberColumns() > 5000
#ifndef CLPSOLVE_ACTIONS
    && (specialOptions_ & 1024) != 0
#endif
  ) {
    // look at original model for objective
    int numberColumns = model2->numberColumns();
    int numberRows = model2->numberRows();
    int numberColumnsOrig = this->numberColumns();
    const double *obj = this->objective();
    int nNon = 0;
    bool allOnes = true;
    for (int i = 0; i < numberColumnsOrig; i++) {
      if (obj[i]) {
        nNon++;
        if (fabs(obj[i]) != 1.0)
          allOnes = false;
      }
    }
    if (nNon <= 1 || allOnes || (options.getExtraInfo(1) > 0 && options.getSpecialOption(1) == 2)) {
#ifdef COIN_DEVELOP
      printf("Forcing primal\n");
#endif
      method = ClpSolve::usePrimal;
#ifndef CLPSOLVE_ACTIONS
      tryItSave = (numberRows > 200 && numberColumns > 2000 && (numberColumns > 2 * numberRows || (specialOptions_ & 1024) != 0)) ? 3 : 0;
#else
      if (numberRows > 200 && numberColumns > 2000 && numberColumns > 2 * numberRows) {
        tryItSave = 3;
      } else {
        // If rhs also rubbish then maybe
        int numberRowsOrig = this->numberRows();
        const double *rowLower = this->rowLower();
        const double *rowUpper = this->rowUpper();
        double last = COIN_DBL_MAX;
        int nDifferent = 0;
        for (int i = 0; i < numberRowsOrig; i++) {
          double value = fabs(rowLower[i]);
          if (value && value != COIN_DBL_MAX) {
            if (value != last) {
              nDifferent++;
              last = value;
            }
          }
          value = fabs(rowUpper[i]);
          if (value && value != COIN_DBL_MAX) {
            if (value != last) {
              nDifferent++;
              last = value;
            }
          }
        }
        if (nDifferent < 2)
          tryItSave = 1;
      }
#endif
    }
  }
  if (method != ClpSolve::useDual && method != ClpSolve::useBarrier
    && method != ClpSolve::tryBenders && method != ClpSolve::tryDantzigWolfe
    && method != ClpSolve::useBarrierNoCross) {
    switch (options.getSpecialOption(1)) {
    case 0:
      doIdiot = -1;
      doCrash = -1;
      doSprint = -1;
      break;
    case 1:
      doIdiot = 0;
      doCrash = 1;
      if (options.getExtraInfo(1) > 0)
        doCrash = options.getExtraInfo(1);
      doSprint = 0;
      break;
    case 2:
      doIdiot = 1;
      if (options.getExtraInfo(1) > 0)
        doIdiot = options.getExtraInfo(1);
      doCrash = 0;
      doSprint = 0;
      break;
    case 3:
      doIdiot = 0;
      doCrash = 0;
      doSprint = 1;
      break;
    case 4:
      doIdiot = 0;
      doCrash = 0;
      doSprint = 0;
      break;
    case 5:
      doIdiot = 0;
      doCrash = -1;
      doSprint = -1;
      break;
    case 6:
      doIdiot = -1;
      doCrash = -1;
      doSprint = 0;
      break;
    case 7:
      doIdiot = -1;
      doCrash = 0;
      doSprint = -1;
      break;
    case 8:
      doIdiot = -1;
      doCrash = 0;
      doSprint = 0;
      break;
    case 9:
      doIdiot = 0;
      doCrash = 0;
      doSprint = -1;
      break;
    case 10:
      doIdiot = 0;
      doCrash = 0;
      doSprint = 0;
      if (options.getExtraInfo(1))
        doSlp = options.getExtraInfo(1);
      break;
    case 11:
      doIdiot = 0;
      doCrash = 0;
      doSprint = 0;
      primalStartup = 0;
      break;
    default:
      abort();
    }
  } else if (method != ClpSolve::tryBenders && method != ClpSolve::tryDantzigWolfe) {
    // Dual
    switch (options.getSpecialOption(0)) {
    case 0:
      doIdiot = 0;
      doCrash = 0;
      doSprint = 0;
      break;
    case 1:
      doIdiot = 0;
      doCrash = 1;
      if (options.getExtraInfo(0) > 0)
        doCrash = options.getExtraInfo(0);
      doSprint = 0;
      break;
    case 2:
      doIdiot = -1;
      if (options.getExtraInfo(0) > 0)
        doIdiot = options.getExtraInfo(0);
      doCrash = 0;
      doSprint = 0;
      break;
    default:
      abort();
    }
  } else {
    // decomposition
  }
#ifndef NO_RTTI
  ClpQuadraticObjective *quadraticObj = (dynamic_cast< ClpQuadraticObjective * >(objectiveAsObject()));
#else
  ClpQuadraticObjective *quadraticObj = NULL;
  if (objective_->type() == 2)
    quadraticObj = (static_cast< ClpQuadraticObjective * >(objective_));
#endif
  // If quadratic then primal or barrier or slp
  if (quadraticObj) {
    doSprint = 0;
    //doIdiot = 0;
    // off
    if (method == ClpSolve::useBarrier)
      method = ClpSolve::useBarrierNoCross;
    else if (method != ClpSolve::useBarrierNoCross)
      method = ClpSolve::usePrimal;
  }
#ifdef COIN_HAS_VOL
  // Save number of idiot
  int saveDoIdiot = doIdiot;
#endif
  // Just do this number of passes in Sprint
  int maxSprintPass = 100;
  // See if worth trying +- one matrix
  bool plusMinus = false;
  CoinBigIndex numberElements = model2->getNumElements();
#ifndef SLIM_CLP
#ifndef NO_RTTI
  if (dynamic_cast< ClpNetworkMatrix * >(matrix_)) {
    // network - switch off stuff
    doIdiot = 0;
    if (doSprint < 0)
      doSprint = 0;
  }
#else
  if (matrix_->type() == 11) {
    // network - switch off stuff
    doIdiot = 0;
    //doSprint=0;
  }
#endif
#endif
  int numberColumns = model2->numberColumns();
  int numberRows = model2->numberRows();
  // If not all slack basis - switch off all except sprint
  int numberRowsBasic = 0;
  int iRow;
  for (iRow = 0; iRow < numberRows; iRow++)
    if (model2->getRowStatus(iRow) == basic)
      numberRowsBasic++;
  if (numberRowsBasic < numberRows && objective_->type() < 2) {
    doIdiot = 0;
    doCrash = 0;
    //doSprint=0;
  }
  if (options.getSpecialOption(3) == 0) {
    if (numberElements > 100000)
      plusMinus = true;
    if (numberElements > 10000 && (doIdiot || doSprint))
      plusMinus = true;
  } else if ((specialOptions_ & 1024) != 0) {
    plusMinus = true;
  }
#ifndef SLIM_CLP
  // Statistics (+1,-1, other) - used to decide on strategy if not +-1
  CoinBigIndex statistics[3] = { -1, 0, 0 };
  if (plusMinus) {
    saveMatrix = model2->clpMatrix();
#ifndef NO_RTTI
    ClpPackedMatrix *clpMatrix = dynamic_cast< ClpPackedMatrix * >(saveMatrix);
#else
    ClpPackedMatrix *clpMatrix = NULL;
    if (saveMatrix->type() == 1)
      clpMatrix = static_cast< ClpPackedMatrix * >(saveMatrix);
#endif
    if (clpMatrix) {
      ClpPlusMinusOneMatrix *newMatrix = new ClpPlusMinusOneMatrix(*(clpMatrix->matrix()));
      if (newMatrix->getIndices()) {
        // CHECKME This test of specialOptions and the one above
        // don't seem compatible.
        if ((specialOptions_ & 1024) == 0) {
          model2->replaceMatrix(newMatrix);
        } else {
          // in integer (or abc) - just use for sprint/idiot
          saveMatrix = NULL;
          delete newMatrix;
        }
      } else {
        handler_->message(CLP_MATRIX_CHANGE, messages_)
          << "+- 1"
          << CoinMessageEol;
        CoinMemcpyN(newMatrix->startPositive(), 3, statistics);
        saveMatrix = NULL;
        plusMinus = false;
        delete newMatrix;
      }
    } else {
      saveMatrix = NULL;
      plusMinus = false;
    }
  }
#endif
  if (this->factorizationFrequency() == 200) {
    // User did not touch preset
    model2->defaultFactorizationFrequency();
  } else if (model2 != this) {
    // make sure model2 has correct value
    model2->setFactorizationFrequency(this->factorizationFrequency());
  }
  if (method == ClpSolve::automatic) {
    if (doSprint == 0 && doIdiot == 0) {
      // off
      method = ClpSolve::useDual;
    } else {
      // only do primal if sprint or idiot
      if (doSprint > 0) {
        method = ClpSolve::usePrimalorSprint;
      } else if (doIdiot > 0) {
        method = ClpSolve::usePrimal;
      } else {
        if (numberElements < 500000) {
          // Small problem
          if (numberRows * 10 > numberColumns || numberColumns < 6000
            || (numberRows * 20 > numberColumns && !plusMinus))
            doSprint = 0; // switch off sprint
        } else {
          // larger problem
          if (numberRows * 8 > numberColumns)
            doSprint = 0; // switch off sprint
        }
        // switch off idiot or sprint if any free variable
        // switch off sprint if very few with costs
        // or great variation in cost
        int iColumn;
        const double *columnLower = model2->columnLower();
        const double *columnUpper = model2->columnUpper();
        const double *objective = model2->objective();
        int nObj = 0;
        int nFree = 0;
        double smallestObj = COIN_DBL_MAX;
        double largestObj = 0.0;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (columnLower[iColumn] < -1.0e10 && columnUpper[iColumn] > 1.0e10) {
            nFree++;
          } else if (objective[iColumn]) {
            nObj++;
            smallestObj = std::min(smallestObj, objective[iColumn]);
            largestObj = std::max(largestObj, objective[iColumn]);
          }
        }
        if (nObj * 10 < numberColumns || smallestObj * 10.0 < largestObj)
          doSprint = 0;
        if (nFree)
          doIdiot = 0;
        if (nFree * 10 > numberRows)
          doSprint = 0;
        int nPasses = 0;
        // look at rhs
        int iRow;
        double largest = 0.0;
        double smallest = 1.0e30;
        double largestGap = 0.0;
        //int numberNotE = 0;
        bool notInteger = false;
        for (iRow = 0; iRow < numberRows; iRow++) {
          double value1 = model2->rowLower_[iRow];
          if (value1 && value1 > -1.0e31) {
            largest = std::max(largest, fabs(value1));
            smallest = std::min(smallest, fabs(value1));
            if (fabs(value1 - floor(value1 + 0.5)) > 1.0e-8) {
              notInteger = true;
              break;
            }
          }
          double value2 = model2->rowUpper_[iRow];
          if (value2 && value2 < 1.0e31) {
            largest = std::max(largest, fabs(value2));
            smallest = std::min(smallest, fabs(value2));
            if (fabs(value2 - floor(value2 + 0.5)) > 1.0e-8) {
              notInteger = true;
              break;
            }
          }
          // CHECKME This next bit can't be right...
          if (value2 > value1) {
            //numberNotE++;
            //if (value2 > 1.0e31 || value1 < -1.0e31)
            //   largestGap = COIN_DBL_MAX;
            //else
            //   largestGap = value2 - value1;
          }
        }
        int tryIt = (numberRows > 200 && numberColumns > 2000 && (numberColumns > 2 * numberRows || (method != ClpSolve::useDual && (specialOptions_ & 1024) != 0))) ? 3 : 0;
        tryItSave = tryIt;
        if (numberRows < 1000 && numberColumns < 3000)
          tryIt = 0;
        if (notInteger)
          tryIt = 0;
        if (largest / smallest > 10 || (largest / smallest > 2.0 && largest > 50))
          tryIt = 0;
        if (tryIt) {
          if (largest / smallest > 2.0) {
            nPasses = 10 + numberColumns / 100000;
            nPasses = std::min(nPasses, 50);
            nPasses = std::max(nPasses, 15);
            if (numberRows > 20000 && nPasses > 5) {
              // Might as well go for it
              nPasses = std::max(nPasses, 71);
            } else if (numberRows > 2000 && nPasses > 5) {
              nPasses = std::max(nPasses, 50);
            } else if (numberElements < 3 * numberColumns) {
              nPasses = std::min(nPasses, 10); // probably not worh it
            }
          } else if (largest / smallest > 1.01 || numberElements <= 3 * numberColumns) {
            nPasses = 10 + numberColumns / 1000;
            nPasses = std::min(nPasses, 100);
            nPasses = std::max(nPasses, 30);
            if (numberRows > 25000) {
              // Might as well go for it
              nPasses = std::max(nPasses, 71);
            }
            if (!largestGap)
              nPasses *= 2;
          } else {
            nPasses = 10 + numberColumns / 1000;
            nPasses = std::max(nPasses, 100);
            if (!largestGap)
              nPasses *= 2;
            nPasses = std::min(nPasses, 200);
          }
        }
        //printf("%d rows %d cols plus %c tryIt %c largest %g smallest %g largestGap %g npasses %d sprint %c\n",
        //     numberRows,numberColumns,plusMinus ? 'Y' : 'N',
        //     tryIt ? 'Y' :'N',largest,smallest,largestGap,nPasses,doSprint ? 'Y' :'N');
        //exit(0);
        if (!tryIt || nPasses <= 5)
          doIdiot = 0;
#if CLPSOLVE_ACTIONS
        if (doIdiot && doSprint < 0 && wasAutomatic && 20 * model2->numberRows() > model2->numberColumns())
          doSprint = 0; // switch off sprint
#endif
        if (doSprint) {
          method = ClpSolve::usePrimalorSprint;
        } else if (doIdiot) {
          method = ClpSolve::usePrimal;
        } else {
          method = ClpSolve::useDual;
        }
      }
    }
  }
  if (method == ClpSolve::tryBenders) {
    // Now build model
    int lengthNames = model2->lengthNames();
    model2->setLengthNames(0);
    CoinModel *build = model2->createCoinModel();
    model2->setLengthNames(lengthNames);
    CoinStructuredModel benders;
    build->convertMatrix();
    int numberBlocks = options.independentOption(0);
    benders.setMessageHandler(handler_);
    numberBlocks = benders.decompose(*build, 2, numberBlocks, NULL);
    delete build;
    //exit(0);
    if (numberBlocks) {
      options.setIndependentOption(1, 1); // don't do final clean up
      model2->solveBenders(&benders, options);
      //move solution
      method = ClpSolve::notImplemented;
      time2 = clpGetTime();
      timeCore = time2 - timeX;
      handler_->message(CLP_INTERVAL_TIMING, messages_)
        << "Crossover" << timeCore << time2 - time1
        << CoinMessageEol;
      timeX = time2;
    } else {
      printf("No structure\n");
      method = ClpSolve::useDual;
    }
  } else if (method == ClpSolve::tryDantzigWolfe) {
    abort();
  }
  if (method == ClpSolve::usePrimalorSprint) {
    if (doSprint < 0) {
      if (numberElements < 500000) {
        // Small problem
        if (numberRows * 10 > numberColumns || numberColumns < 6000
          || (numberRows * 20 > numberColumns && !plusMinus))
          method = ClpSolve::usePrimal; // switch off sprint
      } else {
        // larger problem
        if (numberRows * 8 > numberColumns)
          method = ClpSolve::usePrimal; // switch off sprint
        // but make lightweight
        if (numberRows * 10 > numberColumns || numberColumns < 6000
          || (numberRows * 20 > numberColumns && !plusMinus))
          maxSprintPass = 10;
      }
    } else if (doSprint == 0) {
      method = ClpSolve::usePrimal; // switch off sprint
    }
  }
  if (method == ClpSolve::useDual) {
#ifdef CLP_USEFUL_PRINTOUT
    debugInt[6] = 1;
#endif
    double *saveLower = NULL;
    double *saveUpper = NULL;
    if (presolve == ClpSolve::presolveOn) {
      int numberInfeasibilities = model2->tightenPrimalBounds(0.0, 0);
      if (numberInfeasibilities) {
        handler_->message(CLP_INFEASIBLE, messages_)
          << CoinMessageEol;
        delete model2;
        model2 = this;
        presolve = ClpSolve::presolveOff;
      }
    } else if (numberRows_ + numberColumns_ > 5000) {
      // do anyway
      saveLower = new double[numberRows_ + numberColumns_];
      CoinMemcpyN(model2->columnLower(), numberColumns_, saveLower);
      CoinMemcpyN(model2->rowLower(), numberRows_, saveLower + numberColumns_);
      saveUpper = new double[numberRows_ + numberColumns_];
      CoinMemcpyN(model2->columnUpper(), numberColumns_, saveUpper);
      CoinMemcpyN(model2->rowUpper(), numberRows_, saveUpper + numberColumns_);
      int numberInfeasibilities = model2->tightenPrimalBounds();
      if (numberInfeasibilities) {
        handler_->message(CLP_INFEASIBLE, messages_)
          << CoinMessageEol;
        CoinMemcpyN(saveLower, numberColumns_, model2->columnLower());
        CoinMemcpyN(saveLower + numberColumns_, numberRows_, model2->rowLower());
        delete[] saveLower;
        saveLower = NULL;
        CoinMemcpyN(saveUpper, numberColumns_, model2->columnUpper());
        CoinMemcpyN(saveUpper + numberColumns_, numberRows_, model2->rowUpper());
        delete[] saveUpper;
        saveUpper = NULL;
        // return if wanted
        if (options.infeasibleReturn() || (moreSpecialOptions_ & 1) != 0)
          return -1;
      }
    }
#ifndef COIN_HAS_VOL
    // switch off idiot and volume for now
    doIdiot = 0;
#endif
    // pick up number passes
    int nPasses = 0;
    int numberNotE = 0;
#ifndef SLIM_CLP
    if ((doIdiot < 0 && plusMinus) || doIdiot > 0) {
      // See if candidate for idiot
      nPasses = 0;
      Idiot info(*model2);
      info.setMinIntervalStatusUpdate(model2->getMinIntervalProgressUpdate());
      info.setStrategy(idiotOptions | info.getStrategy());
      // Get average number of elements per column
      double ratio = static_cast< double >(numberElements) / static_cast< double >(numberColumns);
      // look at rhs
      int iRow;
      double largest = 0.0;
      double smallest = 1.0e30;
      for (iRow = 0; iRow < numberRows; iRow++) {
        double value1 = model2->rowLower_[iRow];
        if (value1 && value1 > -1.0e31) {
          largest = std::max(largest, fabs(value1));
          smallest = std::min(smallest, fabs(value1));
        }
        double value2 = model2->rowUpper_[iRow];
        if (value2 && value2 < 1.0e31) {
          largest = std::max(largest, fabs(value2));
          smallest = std::min(smallest, fabs(value2));
        }
        if (value2 > value1) {
          numberNotE++;
        }
      }
      if (doIdiot < 0) {
        if (numberRows > 200 && numberColumns > 5000 && ratio >= 3.0 && largest / smallest < 1.1 && !numberNotE) {
          nPasses = 71;
        }
      }
      if (doIdiot > 0) {
        nPasses = std::max(nPasses, doIdiot);
        if (nPasses > 70) {
          info.setStartingWeight(1.0e3);
          info.setDropEnoughFeasibility(0.01);
        }
      }
      if (nPasses > 20) {
#ifdef COIN_HAS_VOL
        int returnCode = solveWithVolume(model2, nPasses, saveDoIdiot);
        if (!returnCode) {
          time2 = clpGetTime();
          timeIdiot = time2 - timeX;
          handler_->message(CLP_INTERVAL_TIMING, messages_)
            << "Idiot Crash" << timeIdiot << time2 - time1
            << CoinMessageEol;
          timeX = time2;
        } else {
          nPasses = 0;
        }
#else
        nPasses = 0;
#endif
      } else {
        nPasses = 0;
      }
    }
#endif
    if (doCrash) {
        switch (doCrash) {
          // standard
        case 1:
          model2->crash(1000, 1);
          break;
          // As in paper by Solow and Halim (approx)
        case 2:
        case 3:
          model2->crash(model2->dualBound(), 0);
          break;
          // Just put free in basis
        case 4:
          model2->crash(0.0, 3);
          break;
          // Move zero cost variables to minimize primal infeasibility
        case 5:
          model2->crash(0.0, 4);
          break;
          // Put singletons in basis to minimize primal infeasibility
        case 6:
          model2->crash(0.0, 5);
          break;
        }
    }
    if (!nPasses) {
      int saveOptions = model2->specialOptions();
      if (model2->numberRows() > 100)
        model2->setSpecialOptions(saveOptions | 64); // go as far as possible
      //int numberRows = model2->numberRows();
      //int numberColumns = model2->numberColumns();
      if (dynamic_cast< ClpPackedMatrix * >(matrix_)) {
        // See if original wanted vector
        ClpPackedMatrix *clpMatrixO = dynamic_cast< ClpPackedMatrix * >(matrix_);
        ClpMatrixBase *matrix = model2->clpMatrix();
        if (dynamic_cast< ClpPackedMatrix * >(matrix) && clpMatrixO->wantsSpecialColumnCopy()) {
          ClpPackedMatrix *clpMatrix = dynamic_cast< ClpPackedMatrix * >(matrix);
          clpMatrix->makeSpecialColumnCopy();
          //model2->setSpecialOptions(model2->specialOptions()|256); // to say no row copy for comparisons
          model2->dual(0);
          clpMatrix->releaseSpecialColumnCopy();
        } else {
          model2->dual(0);
        }
      } else {
        model2->dual(0);
      }
    } else if (!numberNotE && 0) {
      // E so we can do in another way
      double *pi = model2->dualRowSolution();
      int i;
      int numberColumns = model2->numberColumns();
      int numberRows = model2->numberRows();
      double *saveObj = new double[numberColumns];
      CoinMemcpyN(model2->objective(), numberColumns, saveObj);
      CoinMemcpyN(model2->objective(),
        numberColumns, model2->dualColumnSolution());
      model2->clpMatrix()->transposeTimes(-1.0, pi, model2->dualColumnSolution());
      CoinMemcpyN(model2->dualColumnSolution(),
        numberColumns, model2->objective());
      const double *rowsol = model2->primalRowSolution();
      double offset = 0.0;
      for (i = 0; i < numberRows; i++) {
        offset += pi[i] * rowsol[i];
      }
      double value2;
      model2->getDblParam(ClpObjOffset, value2);
      //printf("Offset %g %g\n",offset,value2);
      model2->setDblParam(ClpObjOffset, value2 - offset);
      model2->setPerturbation(51);
      //model2->setRowObjective(pi);
      // zero out pi
      //memset(pi,0,numberRows*sizeof(double));
      // Could put some in basis - only partially tested
      model2->allSlackBasis();
      //model2->factorization()->maximumPivots(200);
      //model2->setLogLevel(63);
      // solve
      model2->dual(0);
      model2->setDblParam(ClpObjOffset, value2);
      CoinMemcpyN(saveObj, numberColumns, model2->objective());
      // zero out pi
      //memset(pi,0,numberRows*sizeof(double));
      //model2->setRowObjective(pi);
      delete[] saveObj;
      //model2->dual(0);
      model2->setPerturbation(50);
      model2->primal();
    } else {
      // solve
      model2->setPerturbation(100);
      model2->dual(2);
      model2->setPerturbation(50);
      model2->dual(0);
    }
    if (saveLower) {
      CoinMemcpyN(saveLower, numberColumns_, model2->columnLower());
      CoinMemcpyN(saveLower + numberColumns_, numberRows_, model2->rowLower());
      delete[] saveLower;
      saveLower = NULL;
      CoinMemcpyN(saveUpper, numberColumns_, model2->columnUpper());
      CoinMemcpyN(saveUpper + numberColumns_, numberRows_, model2->rowUpper());
      delete[] saveUpper;
      saveUpper = NULL;
    }
    time2 = clpGetTime();
    timeCore = time2 - timeX;
    handler_->message(CLP_INTERVAL_TIMING, messages_)
      << "Dual" << timeCore << time2 - time1
      << CoinMessageEol;
    timeX = time2;
  } else if (method == ClpSolve::usePrimal) {
#ifdef CLP_USEFUL_PRINTOUT
    debugInt[6] = 2;
#endif
#ifndef SLIM_CLP
    if (doIdiot) {
      int nPasses = 0;
      Idiot info(*model2);
      info.setMinIntervalStatusUpdate(model2->getMinIntervalProgressUpdate());
      info.setStrategy(idiotOptions | info.getStrategy());
      // Get average number of elements per column
      double ratio = static_cast< double >(numberElements) / static_cast< double >(numberColumns);
      // look at rhs
      int iRow;
      double largest = 0.0;
      double smallest = 1.0e30;
      double largestGap = 0.0;
      int numberNotE = 0;
      for (iRow = 0; iRow < numberRows; iRow++) {
        double value1 = model2->rowLower_[iRow];
        if (value1 && value1 > -1.0e31) {
          largest = std::max(largest, fabs(value1));
          smallest = std::min(smallest, fabs(value1));
        }
        double value2 = model2->rowUpper_[iRow];
        if (value2 && value2 < 1.0e31) {
          largest = std::max(largest, fabs(value2));
          smallest = std::min(smallest, fabs(value2));
        }
        if (value2 > value1) {
          numberNotE++;
          if (value2 > 1.0e31 || value1 < -1.0e31)
            largestGap = COIN_DBL_MAX;
          else
            largestGap = value2 - value1;
        }
      }
      bool increaseSprint = plusMinus;
      if ((specialOptions_ & 1024) != 0)
        increaseSprint = false;
      if (!plusMinus) {
        // If 90% +- 1 then go for sprint
        if (statistics[0] >= 0 && 10 * statistics[2] < statistics[0] + statistics[1])
          increaseSprint = true;
      }
      int tryIt = tryItSave;
      if (numberRows < 1000 && numberColumns < 3000)
        tryIt = 0;
      if (tryIt) {
        if (increaseSprint) {
          info.setStartingWeight(1.0e3);
          info.setReduceIterations(6);
          // also be more lenient on infeasibilities
          info.setDropEnoughFeasibility(0.5 * info.getDropEnoughFeasibility());
          info.setDropEnoughWeighted(-2.0);
          if (largest / smallest > 2.0) {
            nPasses = 10 + numberColumns / 100000;
            nPasses = std::min(nPasses, 50);
            nPasses = std::max(nPasses, 15);
            if (numberRows > 20000 && nPasses > 5) {
              // Might as well go for it
              nPasses = std::max(nPasses, 71);
            } else if (numberRows > 2000 && nPasses > 5) {
              nPasses = std::max(nPasses, 50);
            } else if (numberElements < 3 * numberColumns) {
              nPasses = std::min(nPasses, 10); // probably not worh it
              if (doIdiot < 0)
                info.setLightweight(1); // say lightweight idiot
            } else {
              if (doIdiot < 0)
                info.setLightweight(1); // say lightweight idiot
            }
          } else if (largest / smallest > 1.01 || numberElements <= 3 * numberColumns) {
            nPasses = 10 + numberColumns / 1000;
            nPasses = std::min(nPasses, 100);
            nPasses = std::max(nPasses, 30);
            if (numberRows > 25000) {
              // Might as well go for it
              nPasses = std::max(nPasses, 71);
            }
            if (!largestGap)
              nPasses *= 2;
          } else {
            nPasses = 10 + numberColumns / 1000;
            nPasses = std::min(nPasses, 200);
            nPasses = std::max(nPasses, 100);
            info.setStartingWeight(1.0e-1);
            info.setReduceIterations(6);
            if (!largestGap && nPasses <= 50)
              nPasses *= 2;
            //info.setFeasibilityTolerance(1.0e-7);
          }
          // If few passes - don't bother
          if (nPasses <= 5 && !plusMinus)
            nPasses = 0;
        } else {
          if (doIdiot < 0)
            info.setLightweight(1); // say lightweight idiot
          if (largest / smallest > 1.01 || numberNotE || statistics[2] > statistics[0] + statistics[1]) {
            if (numberRows > 25000 || numberColumns > 5 * numberRows) {
              nPasses = 50;
            } else if (numberColumns > 4 * numberRows) {
              nPasses = 20;
            } else {
              nPasses = 5;
            }
          } else {
            if (numberRows > 25000 || numberColumns > 5 * numberRows) {
              nPasses = 50;
              info.setLightweight(0); // say not lightweight idiot
            } else if (numberColumns > 4 * numberRows) {
              nPasses = 20;
            } else {
              nPasses = 15;
            }
          }
          if (ratio < 3.0) {
            nPasses = static_cast< int >(ratio * static_cast< double >(nPasses) / 4.0); // probably not worth it
          } else {
            nPasses = std::max(nPasses, 5);
          }
          if (numberRows > 25000 && nPasses > 5) {
            // Might as well go for it
            nPasses = std::max(nPasses, 71);
          } else if (increaseSprint) {
            nPasses *= 2;
            nPasses = std::min(nPasses, 71);
          } else if (nPasses == 5 && ratio > 5.0) {
            nPasses = static_cast< int >(static_cast< double >(nPasses) * (ratio / 5.0)); // increase if lots of elements per column
          }
          if (nPasses <= 5 && !plusMinus)
            nPasses = 0;
          //info.setStartingWeight(1.0e-1);
        }
        if (tryIt == 1) {
          idiotOptions |= 262144;
          info.setStrategy(idiotOptions | info.getStrategy());
          //model2->setSpecialOptions(model2->specialOptions()
          //			|8388608);
        }
      }
      if (doIdiot > 0) {
        // pick up number passes
        nPasses = options.getExtraInfo(1) % 1000000;
#ifdef COIN_HAS_VOL
        int returnCode = solveWithVolume(model2, nPasses, saveDoIdiot);
        nPasses = 0;
        if (!returnCode) {
          time2 = clpGetTime();
          timeIdiot = time2 - timeX;
          handler_->message(CLP_INTERVAL_TIMING, messages_)
            << "Idiot Crash" << timeIdiot << time2 - time1
            << CoinMessageEol;
          timeX = time2;
        }
#endif
#ifdef CLP_USEFUL_PRINTOUT
        debugInt[6] = 4;
        debugInt[7] = nPasses;
#endif
        if (nPasses > 70) {
          info.setStartingWeight(1.0e3);
          info.setReduceIterations(6);
          //if (nPasses > 200)
          //info.setFeasibilityTolerance(1.0e-9);
          //if (nPasses > 1900)
          //info.setWeightFactor(0.93);
          if (nPasses > 900) {
            double reductions = nPasses / 6.0;
            if (nPasses < 5000) {
              reductions /= 12.0;
            } else {
              reductions /= 13.0;
              info.setStartingWeight(1.0e4);
            }
            double ratio = 1.0 / std::pow(10.0, (1.0 / reductions));
            printf("%d passes reduction factor %g\n", nPasses, ratio);
            info.setWeightFactor(ratio);
          } else if (nPasses > 500) {
            info.setWeightFactor(0.7);
          } else if (nPasses > 200) {
            info.setWeightFactor(0.5);
          }
          if (maximumIterations() < nPasses) {
            printf("Presuming maximumIterations is just for Idiot\n");
            nPasses = maximumIterations();
            setMaximumIterations(COIN_INT_MAX);
            model2->setMaximumIterations(COIN_INT_MAX);
          }
          if (nPasses >= 10000 && nPasses < 100000) {
            int k = nPasses % 100;
            nPasses /= 200;
            info.setReduceIterations(3);
            if (k)
              info.setStartingWeight(1.0e2);
          }
          // also be more lenient on infeasibilities
          info.setDropEnoughFeasibility(0.5 * info.getDropEnoughFeasibility());
          info.setDropEnoughWeighted(-2.0);
        } else if (nPasses >= 50) {
          info.setStartingWeight(1.0e3);
          //info.setReduceIterations(6);
        }
        // For experimenting
        if (nPasses < 70 && (nPasses % 10) > 0 && (nPasses % 10) < 4) {
          info.setStartingWeight(1.0e3);
          info.setLightweight(nPasses % 10); // special testing
#ifdef COIN_DEVELOP
          printf("warning - odd lightweight %d\n", nPasses % 10);
          //info.setReduceIterations(6);
#endif
        }
      }
      if (options.getExtraInfo(1) > 1000000)
        nPasses += 1000000;
      if (nPasses) {
        doCrash = 0;
        // Allow for crossover
        //#define LACI_TRY
#ifndef LACI_TRY
        //if (doIdiot>0)
        info.setStrategy(512 | info.getStrategy());
#endif
        // Allow for scaling
        info.setStrategy(32 | info.getStrategy());
        int saveScalingFlag = model2->scalingFlag();
        bool linearObjective = objective_->type() < 2;
        if (!linearObjective)
          model2->scaling(0);
#define CLP_DUAL_IDIOT
#ifdef CLP_DUAL_IDIOT
        bool doubleIdiot = false;
        if (nPasses == 99 && linearObjective)
          doubleIdiot = true;
        if (doubleIdiot) {
          ClpSimplex *dualModel2 = static_cast< ClpSimplexOther * >(model2)->dualOfModel(1.0, 1.0);
          if (dualModel2) {
            //printf("Dual of model has %d rows and %d columns\n",
            //  dualModel2->numberRows(), dualModel2->numberColumns());
            dualModel2->setOptimizationDirection(1.0);
            Idiot infoDual(info);
	    info.setMinIntervalStatusUpdate(dualModel2->getMinIntervalProgressUpdate());
            infoDual.setModel(dualModel2);
            info.crash(nPasses, model2->messageHandler(),
              model2->messagesPointer(), false);
            infoDual.crash(nPasses, model2->messageHandler(),
              model2->messagesPointer(), false);
            // two copies of solutions
            ClpSimplex temp(*model2);
            // move duals and just copy primal
            memcpy(temp.dualRowSolution(), dualModel2->primalColumnSolution(),
              numberRows * sizeof(double));
            memcpy(temp.primalColumnSolution(), model2->primalColumnSolution(),
              numberColumns * sizeof(double));
            delete dualModel2;
            int numberRows = model2->numberRows();
            int numberColumns = model2->numberColumns();
            ClpSimplex *tempModel[2];
            tempModel[0] = model2;
            tempModel[1] = &temp;
            const double *primalColumn[2];
            const double *dualRow[2];
            double *dualColumn[2];
            double *primalRow[2];
            for (int i = 0; i < 2; i++) {
              primalColumn[i] = tempModel[i]->primalColumnSolution();
              dualRow[i] = tempModel[i]->dualRowSolution();
              dualColumn[i] = tempModel[i]->dualColumnSolution();
              primalRow[i] = tempModel[i]->primalRowSolution();
              memcpy(dualColumn[i], model2->objective(),
                numberColumns * sizeof(double));
              memset(primalRow[i], 0, numberRows * sizeof(double));
              tempModel[i]->clpMatrix()->transposeTimes(-1.0,
                dualRow[i],
                dualColumn[i]);
              tempModel[i]->clpMatrix()->times(1.0,
                primalColumn[i],
                primalRow[i]);
              tempModel[i]->checkSolutionInternal();
              //printf("model %d - dual inf %g primal inf %g\n",
              //  i, tempModel[i]->sumDualInfeasibilities(),
              //  tempModel[i]->sumPrimalInfeasibilities());
            }
            //printf("What now\n");
          } else {
            doubleIdiot = false;
          }
        }
        if (!doubleIdiot)
          info.crash(nPasses, model2->messageHandler(), model2->messagesPointer(),
            (objective_->type() < 2));
#else
        info.crash(nPasses, model2->messageHandler(), model2->messagesPointer(), (objective_->type() < 2));
#endif
        model2->scaling(saveScalingFlag);
        time2 = clpGetTime();
        timeIdiot = time2 - timeX;
        handler_->message(CLP_INTERVAL_TIMING, messages_)
          << "Idiot Crash" << timeIdiot << time2 - time1
          << CoinMessageEol;
        timeX = time2;
        if (nPasses > 100000 && nPasses < 100500) {
          // make sure no status left
          model2->createStatus();
          // solve
          if (model2->factorizationFrequency() == 200) {
            // User did not touch preset
            model2->defaultFactorizationFrequency();
          }
          //int numberRows = model2->numberRows();
          int numberColumns = model2->numberColumns();
          // save duals
          //double * saveDuals = CoinCopyOfArray(model2->dualRowSolution(),numberRows);
          // for moment this only works on nug etc (i.e. all ==)
          // needs row objective
          double *saveObj = CoinCopyOfArray(model2->objective(), numberColumns);
          double *pi = model2->dualRowSolution();
          model2->clpMatrix()->transposeTimes(-1.0, pi, model2->objective());
          // just primal values pass
          double saveScale = model2->objectiveScale();
          model2->setObjectiveScale(1.0e-3);
          model2->primal(2);
          model2->writeMps("xx.mps");
          double *solution = model2->primalColumnSolution();
          double *upper = model2->columnUpper();
          for (int i = 0; i < numberColumns; i++) {
            if (solution[i] < 100.0)
              upper[i] = 1000.0;
          }
          model2->setProblemStatus(-1);
          model2->setObjectiveScale(saveScale);
          memcpy(model2->objective(), saveObj, numberColumns * sizeof(double));
          //delete [] saveDuals;
          delete[] saveObj;
          model2->dual(2);
        } // end dubious idiot
      }
    }
#endif
    // ?
    if (doCrash) {
      switch (doCrash) {
        // standard
      case 1:
        model2->crash(1000, 1);
        break;
        // As in paper by Solow and Halim (approx)
      case 2:
        model2->crash(model2->dualBound(), 0);
        break;
        // My take on it
      case 3:
        model2->crash(model2->dualBound(), -1);
        break;
        // Just put free in basis
      case 4:
        model2->crash(0.0, 3);
        break;
	// Move zero cost variables to minimize primal infeasibility
      case 5:
	model2->crash(0.0, 4);
	break;
	// Put singletons in basis to minimize primal infeasibility
      case 6:
	model2->crash(0.0, 5);
	break;
      }
    }
#ifndef SLIM_CLP
    if (doSlp && objective_->type() == 2) {
      model2->nonlinearSLP(doSlp, 1.0e-5);
    }
#endif
#ifndef LACI_TRY
    if (model2->status()&&
	(options.getSpecialOption(1) != 2 || options.getExtraInfo(1) < 1000000)) {
      if (dynamic_cast< ClpPackedMatrix * >(matrix_)) {
        // See if original wanted vector
        ClpPackedMatrix *clpMatrixO = dynamic_cast< ClpPackedMatrix * >(matrix_);
        ClpMatrixBase *matrix = model2->clpMatrix();
        if (dynamic_cast< ClpPackedMatrix * >(matrix) && clpMatrixO->wantsSpecialColumnCopy()) {
          ClpPackedMatrix *clpMatrix = dynamic_cast< ClpPackedMatrix * >(matrix);
          clpMatrix->makeSpecialColumnCopy();
          //model2->setSpecialOptions(model2->specialOptions()|256); // to say no row copy for comparisons
          model2->primal(primalStartup);
          clpMatrix->releaseSpecialColumnCopy();
        } else {
	  model2->primal(primalStartup);
        }
      } else {
        model2->primal(primalStartup);
      }
    }
#endif
    time2 = clpGetTime();
    timeCore = time2 - timeX;
    handler_->message(CLP_INTERVAL_TIMING, messages_)
      << "Primal" << timeCore << time2 - time1
      << CoinMessageEol;
    timeX = time2;
  } else if (method == ClpSolve::usePrimalorSprint) {
    // Sprint
    /*
            This driver implements what I called Sprint when I introduced the idea
            many years ago.  Cplex calls it "sifting" which I think is just as silly.
            When I thought of this trivial idea
            it reminded me of an LP code of the 60's called sprint which after
            every factorization took a subset of the matrix into memory (all
            64K words!) and then iterated very fast on that subset.  On the
            problems of those days it did not work very well, but it worked very
            well on aircrew scheduling problems where there were very large numbers
            of columns all with the same flavor.
          */

    /* The idea works best if you can get feasible easily.  To make it
             more general we can add in costed slacks */

    int originalNumberColumns = model2->numberColumns();
    int numberRows = model2->numberRows();
    ClpSimplex *originalModel2 = model2;

    // We will need arrays to choose variables.  These are too big but ..
    double *weight = new double[numberRows + originalNumberColumns];
    int *sort = new int[numberRows + originalNumberColumns];
    int numberSort = 0;
    // We are going to add slacks to get feasible.
    // initial list will just be artificials
    int iColumn;
    const double *columnLower = model2->columnLower();
    const double *columnUpper = model2->columnUpper();
    double *columnSolution = model2->primalColumnSolution();

    // See if we have costed slacks
    int *negSlack = new int[numberRows + numberColumns];
    int *posSlack = new int[numberRows + numberColumns];
    int *nextNegSlack = negSlack + numberRows;
    int *nextPosSlack = posSlack + numberRows;
    int iRow;
    for (iRow = 0; iRow < numberRows + numberColumns; iRow++) {
      negSlack[iRow] = -1;
      posSlack[iRow] = -1;
    }
    const double *element = model2->matrix()->getElements();
    const int *row = model2->matrix()->getIndices();
    const CoinBigIndex *columnStart = model2->matrix()->getVectorStarts();
    const int *columnLength = model2->matrix()->getVectorLengths();
    //bool allSlack = (numberRowsBasic==numberRows);
    for (iColumn = 0; iColumn < originalNumberColumns; iColumn++) {
      if (!columnSolution[iColumn] || fabs(columnSolution[iColumn]) > 1.0e20) {
        double value = 0.0;
        if (columnLower[iColumn] > 0.0)
          value = columnLower[iColumn];
        else if (columnUpper[iColumn] < 0.0)
          value = columnUpper[iColumn];
        columnSolution[iColumn] = value;
      }
      if (columnLength[iColumn] == 1) {
        int jRow = row[columnStart[iColumn]];
        if (!columnLower[iColumn]) {
          if (element[columnStart[iColumn]] > 0.0) {
            if (posSlack[jRow] < 0) {
              posSlack[jRow] = iColumn;
            } else {
              int jColumn = posSlack[jRow];
              while (nextPosSlack[jColumn] >= 0)
                jColumn = nextPosSlack[jColumn];
              nextPosSlack[jColumn] = iColumn;
            }
          } else if (element[columnStart[iColumn]] < 0.0) {
            if (negSlack[jRow] < 0) {
              negSlack[jRow] = iColumn;
            } else {
              int jColumn = negSlack[jRow];
              while (nextNegSlack[jColumn] >= 0)
                jColumn = nextNegSlack[jColumn];
              nextNegSlack[jColumn] = iColumn;
            }
          }
        } else if (!columnUpper[iColumn] && false) { // out for testing
          if (element[columnStart[iColumn]] < 0.0 && posSlack[jRow] < 0)
            posSlack[jRow] = iColumn;
          else if (element[columnStart[iColumn]] > 0.0 && negSlack[jRow] < 0)
            negSlack[jRow] = iColumn;
        }
      }
    }
    // now see what that does to row solution
    const double *objective = model2->objective();
    double *rowSolution = model2->primalRowSolution();
    CoinZeroN(rowSolution, numberRows);
    model2->clpMatrix()->times(1.0, columnSolution, rowSolution);
    // See if we can adjust using costed slacks
    double penalty = std::max(1.0e5, std::min(infeasibilityCost_ * 0.01, 1.0e10)) * optimizationDirection_;
    const double *lower = model2->rowLower();
    const double *upper = model2->rowUpper();
    for (iRow = 0; iRow < numberRows; iRow++) {
      if (lower[iRow] > rowSolution[iRow] + 1.0e-8) {
        int jColumn = posSlack[iRow];
        if (jColumn >= 0) {
          // sort if more than one
          int nPos = 1;
          sort[0] = jColumn;
          weight[0] = objective[jColumn];
          while (nextPosSlack[jColumn] >= 0) {
            jColumn = nextPosSlack[jColumn];
            sort[nPos] = jColumn;
            weight[nPos++] = objective[jColumn];
          }
          if (nPos > 1) {
            CoinSort_2(weight, weight + nPos, sort);
            for (int i = 0; i < nPos; i++) {
              double difference = lower[iRow] - rowSolution[iRow];
              jColumn = sort[i];
              double elementValue = element[columnStart[jColumn]];
              assert(elementValue > 0.0);
              double value = columnSolution[jColumn];
              double movement = std::min(difference / elementValue, columnUpper[jColumn] - value);
              columnSolution[jColumn] += movement;
              rowSolution[iRow] += movement * elementValue;
            }
            continue;
          }
          if (jColumn < 0 || columnSolution[jColumn])
            continue;
          double difference = lower[iRow] - rowSolution[iRow];
          double elementValue = element[columnStart[jColumn]];
          if (elementValue > 0.0) {
            double movement = std::min(difference / elementValue, columnUpper[jColumn]);
            columnSolution[jColumn] = movement;
            rowSolution[iRow] += movement * elementValue;
          } else {
            double movement = std::max(difference / elementValue, columnLower[jColumn]);
            columnSolution[jColumn] = movement;
            rowSolution[iRow] += movement * elementValue;
          }
        }
      } else if (upper[iRow] < rowSolution[iRow] - 1.0e-8) {
        int jColumn = negSlack[iRow];
        if (jColumn >= 0) {
          // sort if more than one
          int nNeg = 1;
          sort[0] = jColumn;
          weight[0] = objective[jColumn];
          while (nextNegSlack[jColumn] >= 0) {
            jColumn = nextNegSlack[jColumn];
            sort[nNeg] = jColumn;
            weight[nNeg++] = objective[jColumn];
          }
          if (nNeg > 1) {
            CoinSort_2(weight, weight + nNeg, sort);
            for (int i = 0; i < nNeg; i++) {
              double difference = rowSolution[iRow] - upper[iRow];
              jColumn = sort[i];
              double elementValue = element[columnStart[jColumn]];
              assert(elementValue < 0.0);
              double value = columnSolution[jColumn];
              double movement = std::min(difference / -elementValue, columnUpper[jColumn] - value);
              columnSolution[jColumn] += movement;
              rowSolution[iRow] += movement * elementValue;
            }
            continue;
          }
          if (jColumn < 0 || columnSolution[jColumn])
            continue;
          double difference = upper[iRow] - rowSolution[iRow];
          double elementValue = element[columnStart[jColumn]];
          if (elementValue < 0.0) {
            double movement = std::min(difference / elementValue, columnUpper[jColumn]);
            columnSolution[jColumn] = movement;
            rowSolution[iRow] += movement * elementValue;
          } else {
            double movement = std::max(difference / elementValue, columnLower[jColumn]);
            columnSolution[jColumn] = movement;
            rowSolution[iRow] += movement * elementValue;
          }
        }
      }
    }
    delete[] negSlack;
    delete[] posSlack;
    int nRow = numberRows;
    bool network = false;
    if (dynamic_cast< ClpNetworkMatrix * >(matrix_)) {
      network = true;
      nRow *= 2;
    }
    CoinBigIndex *addStarts = new CoinBigIndex[nRow + 1];
    int *addRow = new int[nRow];
    double *addElement = new double[nRow];
    addStarts[0] = 0;
    int numberArtificials = 0;
    int numberAdd = 0;
    double *addCost = new double[numberRows];
    for (iRow = 0; iRow < numberRows; iRow++) {
      if (lower[iRow] > rowSolution[iRow] + 1.0e-8) {
        addRow[numberAdd] = iRow;
        addElement[numberAdd++] = 1.0;
        if (network) {
          addRow[numberAdd] = numberRows;
          addElement[numberAdd++] = -1.0;
        }
        addCost[numberArtificials] = penalty;
        numberArtificials++;
        addStarts[numberArtificials] = numberAdd;
      } else if (upper[iRow] < rowSolution[iRow] - 1.0e-8) {
        addRow[numberAdd] = iRow;
        addElement[numberAdd++] = -1.0;
        if (network) {
          addRow[numberAdd] = numberRows;
          addElement[numberAdd++] = 1.0;
        }
        addCost[numberArtificials] = penalty;
        numberArtificials++;
        addStarts[numberArtificials] = numberAdd;
      }
    }
    if (numberArtificials) {
      // need copy so as not to disturb original
      model2 = new ClpSimplex(*model2);
      if (network) {
        // network - add a null row
        model2->addRow(0, NULL, NULL, -COIN_DBL_MAX, COIN_DBL_MAX);
        numberRows++;
      }
      model2->addColumns(numberArtificials, NULL, NULL, addCost,
        addStarts, addRow, addElement);
    }
    // redo
    element = model2->matrix()->getElements();
    row = model2->matrix()->getIndices();
    columnStart = model2->matrix()->getVectorStarts();
    columnLength = model2->matrix()->getVectorLengths();
    delete[] addStarts;
    delete[] addRow;
    delete[] addElement;
    delete[] addCost;
    // look at rhs to see if to perturb
    double largest = 0.0;
    double smallest = 1.0e30;
    for (iRow = 0; iRow < numberRows; iRow++) {
      double value;
      value = fabs(model2->rowLower_[iRow]);
      if (value && value < 1.0e30) {
        largest = std::max(largest, value);
        smallest = std::min(smallest, value);
      }
      value = fabs(model2->rowUpper_[iRow]);
      if (value && value < 1.0e30) {
        largest = std::max(largest, value);
        smallest = std::min(smallest, value);
      }
    }
    double *saveLower = NULL;
    double *saveUpper = NULL;
    if (largest < -2.01 * smallest) { // switch off for now
      // perturb - so switch off standard
      model2->setPerturbation(100);
      saveLower = new double[numberRows];
      CoinMemcpyN(model2->rowLower_, numberRows, saveLower);
      saveUpper = new double[numberRows];
      CoinMemcpyN(model2->rowUpper_, numberRows, saveUpper);
      double *lower = model2->rowLower();
      double *upper = model2->rowUpper();
      for (iRow = 0; iRow < numberRows; iRow++) {
        double lowerValue = lower[iRow], upperValue = upper[iRow];
        double value = randomNumberGenerator_.randomDouble();
        if (upperValue > lowerValue + primalTolerance_) {
          if (lowerValue > -1.0e20 && lowerValue)
            lowerValue -= value * 1.0e-4 * fabs(lowerValue);
          if (upperValue < 1.0e20 && upperValue)
            upperValue += value * 1.0e-4 * fabs(upperValue);
        } else if (upperValue > 0.0) {
          upperValue -= value * 1.0e-4 * fabs(lowerValue);
          lowerValue -= value * 1.0e-4 * fabs(lowerValue);
        } else if (upperValue < 0.0) {
          upperValue += value * 1.0e-4 * fabs(lowerValue);
          lowerValue += value * 1.0e-4 * fabs(lowerValue);
        } else {
        }
        lower[iRow] = lowerValue;
        upper[iRow] = upperValue;
      }
    }
    int i;
    // Just do this number of passes in Sprint
    if (doSprint > 0)
      maxSprintPass = options.getExtraInfo(1);
    // but if big use to get ratio
    double ratio = 3;
#ifdef CLP_USEFUL_PRINTOUT
    debugInt[6] = 3;
    debugInt[7] = maxSprintPass;
#endif
    if (maxSprintPass > 1000) {
      ratio = static_cast< double >(maxSprintPass) * 0.0001;
      ratio = std::max(ratio, 1.1);
      maxSprintPass = maxSprintPass % 1000;
#ifdef COIN_DEVELOP
      printf("%d passes wanted with ratio of %g\n", maxSprintPass, ratio);
#endif
    }
    // Just take this number of columns in small problem
    int smallNumberColumns = static_cast< int >(std::min(ratio * numberRows, static_cast< double >(numberColumns)));
    smallNumberColumns = std::max(smallNumberColumns, 3000);
    smallNumberColumns = std::min(smallNumberColumns, numberColumns);
    int saveSmallNumber = smallNumberColumns;
    bool emergencyMode = false;
    //int smallNumberColumns = std::min(12*numberRows/10,numberColumns);
    //smallNumberColumns = std::max(smallNumberColumns,3000);
    //smallNumberColumns = std::max(smallNumberColumns,numberRows+1000);
    // redo as may have changed
    columnLower = model2->columnLower();
    columnUpper = model2->columnUpper();
    columnSolution = model2->primalColumnSolution();
    // Set up initial list
    numberSort = 0;
    if (numberArtificials) {
      numberSort = numberArtificials;
      for (i = 0; i < numberSort; i++)
        sort[i] = i + originalNumberColumns;
    }
    // put in free
    // maybe a solution there already
    for (iColumn = 0; iColumn < originalNumberColumns; iColumn++) {
      if (model2->getColumnStatus(iColumn) == basic || (columnLower[iColumn] < -1.0e30 && columnUpper[iColumn] > 1.0e30))
        sort[numberSort++] = iColumn;
    }
    for (iColumn = 0; iColumn < originalNumberColumns; iColumn++) {
      if (model2->getColumnStatus(iColumn) != basic) {
        if (columnSolution[iColumn] > columnLower[iColumn] && columnSolution[iColumn] < columnUpper[iColumn] && columnSolution[iColumn])
          sort[numberSort++] = iColumn;
      }
    }
    numberSort = std::min(numberSort, smallNumberColumns);

    int numberColumns = model2->numberColumns();
    double *fullSolution = model2->primalColumnSolution();

    int iPass;
    double lastObjective[] = { 1.0e31, 1.0e31 };
    // It will be safe to allow dense
    model2->setInitialDenseFactorization(true);

    // We will be using all rows
    int *whichRows = new int[numberRows];
    for (iRow = 0; iRow < numberRows; iRow++)
      whichRows[iRow] = iRow;
    double originalOffset;
    model2->getDblParam(ClpObjOffset, originalOffset);
    int totalIterations = 0;
    double lastSumArtificials = COIN_DBL_MAX;
    int originalMaxSprintPass = maxSprintPass;
    maxSprintPass = std::max(0,maxSprintPass); // so we do that many if infeasible
    for (iPass = 0; iPass < maxSprintPass; iPass++) {
      //printf("Bug until submodel new version\n");
      //CoinSort_2(sort,sort+numberSort,weight);
      // Create small problem
      ClpSimplex small(model2, numberRows, whichRows, numberSort, sort);
      small.setPerturbation(model2->perturbation());
      small.setInfeasibilityCost(model2->infeasibilityCost());
      if (model2->factorizationFrequency() == 200) {
        // User did not touch preset
        small.defaultFactorizationFrequency();
      }
      // now see what variables left out do to row solution
      double *rowSolution = model2->primalRowSolution();
      double *sumFixed = new double[numberRows];
      CoinZeroN(sumFixed, numberRows);
      int iRow, iColumn;
      // zero out ones in small problem
      for (iColumn = 0; iColumn < numberSort; iColumn++) {
        int kColumn = sort[iColumn];
        fullSolution[kColumn] = 0.0;
      }
      // Get objective offset
      const double *objective = model2->objective();
      double offset = 0.0;
      for (iColumn = 0; iColumn < originalNumberColumns; iColumn++)
        offset += fullSolution[iColumn] * objective[iColumn];
      small.setDblParam(ClpObjOffset, originalOffset - offset);
      int smallMore = small.moreSpecialOptions();
      smallMore &= ~1048576; // make sure can't stop early
      small.setMoreSpecialOptions(smallMore);
      model2->clpMatrix()->times(1.0, fullSolution, sumFixed);

      double *lower = small.rowLower();
      double *upper = small.rowUpper();
      for (iRow = 0; iRow < numberRows; iRow++) {
        if (lower[iRow] > -1.0e50)
          lower[iRow] -= sumFixed[iRow];
        if (upper[iRow] < 1.0e50)
          upper[iRow] -= sumFixed[iRow];
        rowSolution[iRow] -= sumFixed[iRow];
      }
      delete[] sumFixed;
      // Solve
      if (interrupt)
        currentModel = &small;
      small.defaultFactorizationFrequency();
      if (emergencyMode) {
        // not much happening so big model
        int options = small.moreSpecialOptions();
        small.setMoreSpecialOptions(options | 1048576);
        small.setMaximumIterations(1000400);
        small.setPerturbation(100);
      }
      if (dynamic_cast< ClpPackedMatrix * >(matrix_)) {
        // See if original wanted vector
        ClpPackedMatrix *clpMatrixO = dynamic_cast< ClpPackedMatrix * >(matrix_);
        ClpMatrixBase *matrix = small.clpMatrix();
        if (dynamic_cast< ClpPackedMatrix * >(matrix) && clpMatrixO->wantsSpecialColumnCopy()) {
          ClpPackedMatrix *clpMatrix = dynamic_cast< ClpPackedMatrix * >(matrix);
          clpMatrix->makeSpecialColumnCopy();
	  small.setMaximumIterations(std::max(small.numberRows(),1000));
          small.primal(1);
          clpMatrix->releaseSpecialColumnCopy();
        } else {
#if 1
          if (iPass || !numberArtificials)
            small.primal(1);
          else
            small.dual(0);
          if (emergencyMode) {
            if (small.problemStatus() == 3)
              small.setProblemStatus(0);
            smallNumberColumns = saveSmallNumber;
            emergencyMode = false;
            double *temp = new double[numberRows];
            memset(temp, 0, numberRows * sizeof(double));
            double *solution = small.primalColumnSolution();
            small.matrix()->times(solution, temp);
            double sumInf = 0.0;
            double *lower = small.rowLower();
            double *upper = small.rowUpper();
            for (int iRow = 0; iRow < numberRows; iRow++) {
              if (temp[iRow] > upper[iRow])
                sumInf += temp[iRow] - upper[iRow];
              else if (temp[iRow] < lower[iRow])
                sumInf += lower[iRow] - temp[iRow];
            }
            printf("row inf %g\n", sumInf);
            sumInf = 0.0;
            lower = small.columnLower();
            upper = small.columnUpper();
            for (int iColumn = 0; iColumn < small.numberColumns(); iColumn++) {
              if (solution[iColumn] > upper[iColumn])
                sumInf += solution[iColumn] - upper[iColumn];
              else if (solution[iColumn] < lower[iColumn])
                sumInf += lower[iColumn] - solution[iColumn];
            }
            printf("column inf %g\n", sumInf);
            delete[] temp;
          }
          if (small.problemStatus()) {
            int numberIterations = small.numberIterations();
            small.dual(0);
            small.setNumberIterations(small.numberIterations() + numberIterations);
          }
#else
          int numberColumns = small.numberColumns();
          int numberRows = small.numberRows();
          // Use dual region
          double *rhs = small.dualRowSolution();
          int *whichRow = new int[3 * numberRows];
          int *whichColumn = new int[2 * numberColumns];
          int nBound;
          ClpSimplex *small2 = ((ClpSimplexOther *)(&small))->crunch(rhs, whichRow, whichColumn, nBound, false, false);
          if (small2) {
            small.primal(1);
            if (small2->problemStatus() == 0) {
              small.setProblemStatus(0);
              ((ClpSimplexOther *)(&small))->afterCrunch(*small2, whichRow, whichColumn, nBound);
            } else {
              small.primal(1);
              if (small2->problemStatus())
                small.primal(1);
            }
            delete small2;
          } else {
            small.primal(1);
          }
          delete[] whichRow;
          delete[] whichColumn;
#endif
        }
      } else {
        small.primal(1);
      }
      int smallIterations = small.numberIterations();
      totalIterations += smallIterations;
      if ((2 * smallIterations < std::min(numberRows, 1000)||small.status()==3) && iPass) {
        int oldNumber = smallNumberColumns;
        if (smallIterations < 100)
          smallNumberColumns *= 1.2;
        else
          smallNumberColumns *= 1.1;
        smallNumberColumns = std::min(smallNumberColumns, numberColumns);
        if (smallIterations < 200) {
          // try kicking it
          smallNumberColumns *= 1.5;
        }
        char line[100];
        sprintf(line, "sample size increased from %d to %d",
          oldNumber, smallNumberColumns);
        handler_->message(CLP_GENERAL, messages_)
          << line
          << CoinMessageEol;
      }
      // move solution back
      const double *solution = small.primalColumnSolution();
      for (iColumn = 0; iColumn < numberSort; iColumn++) {
        int kColumn = sort[iColumn];
        model2->setColumnStatus(kColumn, small.getColumnStatus(iColumn));
        fullSolution[kColumn] = solution[iColumn];
      }
      for (iRow = 0; iRow < numberRows; iRow++)
        model2->setRowStatus(iRow, small.getRowStatus(iRow));
      CoinMemcpyN(small.primalRowSolution(),
        numberRows, model2->primalRowSolution());
      double sumArtificials = 0.0;
      for (i = 0; i < numberArtificials; i++)
        sumArtificials += fullSolution[i + originalNumberColumns];
      if (sumArtificials && iPass > 5 && sumArtificials >= lastSumArtificials) {
        // increase costs
        double *cost = model2->objective() + originalNumberColumns;
        double newCost = std::min(1.0e10, cost[0] * 1.5);
        for (i = 0; i < numberArtificials; i++)
          cost[i] = newCost;
      }
      lastSumArtificials = sumArtificials;
      // get reduced cost for large problem
      double *djs = model2->dualColumnSolution();
      CoinMemcpyN(model2->objective(), numberColumns, djs);
      model2->clpMatrix()->transposeTimes(-1.0, small.dualRowSolution(), djs);
      int numberNegative = 0;
      double sumNegative = 0.0;
      // now massage weight so all basic in plus good djs
      // first count and do basic
      numberSort = 0;
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        double dj = djs[iColumn] * optimizationDirection_;
        double value = fullSolution[iColumn];
        if (model2->getColumnStatus(iColumn) == ClpSimplex::basic) {
          sort[numberSort++] = iColumn;
        } else if (dj < -dualTolerance_ && value < columnUpper[iColumn]) {
          numberNegative++;
          sumNegative -= dj;
        } else if (dj > dualTolerance_ && value > columnLower[iColumn]) {
          numberNegative++;
          sumNegative += dj;
        }
      }
      handler_->message(CLP_SPRINT, messages_)
        << iPass + 1 << small.numberIterations() << small.objectiveValue() << sumNegative
        << numberNegative << sumArtificials << small.numberColumns()
        << CoinMessageEol;
      if (sumArtificials < 1.0e-8 && originalMaxSprintPass >= 0) {
        maxSprintPass = iPass + originalMaxSprintPass;
        originalMaxSprintPass = -1;
      }
      if ((small.objectiveValue() * optimizationDirection_ > lastObjective[1] - 1.0e-7 && iPass > 15 && sumArtificials < 1.0e-8 && maxSprintPass < 200) || (!small.numberIterations() && iPass) || iPass == maxSprintPass - 1) {

        break; // finished
      } else {
        lastObjective[1] = lastObjective[0];
        lastObjective[0] = small.objectiveValue() * optimizationDirection_;
        double tolerance;
        double averageNegDj = sumNegative / static_cast< double >(numberNegative + 1);
        if (numberNegative + numberSort > smallNumberColumns && false)
          tolerance = -dualTolerance_;
        else
          tolerance = 10.0 * averageNegDj;
        if (emergencyMode)
          tolerance = 1.0e100;
        int saveN = numberSort;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          double dj = djs[iColumn] * optimizationDirection_;
	  ClpSimplex::Status colStatus = model2->getColumnStatus(iColumn);
	  if (colStatus == ClpSimplex::isFree || colStatus == ClpSimplex::superBasic) {
	    dj = - fabs(dj);
	  } else if (colStatus == ClpSimplex::atLowerBound) {
	    dj = dj;
	  } else if (colStatus == ClpSimplex::atUpperBound) {
	    dj = -dj;
	  } else {
	    continue;
	  }
	  if (dj < tolerance) {
	    weight[numberSort] = dj;
	    sort[numberSort++] = iColumn;
	  }
	}
        // sort
        CoinSort_2(weight + saveN, weight + numberSort, sort + saveN);
        //if (numberSort < smallNumberColumns)
	//printf("using %d columns not %d\n", numberSort, smallNumberColumns);
        numberSort = std::min(smallNumberColumns, numberSort);
        // try singletons
        char *markX = new char[numberColumns];
        memset(markX, 0, numberColumns);
        for (int i = 0; i < numberSort; i++)
          markX[sort[i]] = 1;
        int n = numberSort;
        for (int i = 0; i < numberColumns; i++) {
          if (columnLength[i] == 1 && !markX[i])
            sort[numberSort++] = i;
        }
        //if (n < numberSort)
	//printf("%d slacks added\n", numberSort - n);
        delete[] markX;
      }
    }
    if (interrupt)
      currentModel = model2;
    for (i = 0; i < numberArtificials; i++)
      sort[i] = i + originalNumberColumns;
    model2->deleteColumns(numberArtificials, sort);
    if (network) {
      int iRow = numberRows - 1;
      model2->deleteRows(1, &iRow);
    }
    delete[] weight;
    delete[] sort;
    delete[] whichRows;
    if (saveLower) {
      // unperturb and clean
      for (iRow = 0; iRow < numberRows; iRow++) {
        model2->rowLower_[iRow] = saveLower[iRow];
        model2->rowUpper_[iRow] = saveUpper[iRow];
      }
      delete[] saveLower;
      delete[] saveUpper;
    }
    model2->primal(1);
    model2->setPerturbation(savePerturbation);
    if (model2 != originalModel2) {
      originalModel2->moveInfo(*model2);
      delete model2;
      model2 = originalModel2;
    }
    time2 = clpGetTime();
    timeCore = time2 - timeX;
    handler_->message(CLP_INTERVAL_TIMING, messages_)
      << "Sprint" << timeCore << time2 - time1
      << CoinMessageEol;
    timeX = time2;
    model2->setNumberIterations(model2->numberIterations() + totalIterations);
  } else if (method == ClpSolve::useBarrier || method == ClpSolve::useBarrierNoCross) {
    if (presolve == ClpSolve::presolveOn) {
      int numberInfeasibilities = model2->tightenPrimalBounds(0.0, 0);
      if (numberInfeasibilities) {
        handler_->message(CLP_INFEASIBLE, messages_)
          << CoinMessageEol;
        delete model2;
        model2 = this;
        presolve = ClpSolve::presolveOff;
      }
    }
#ifndef SLIM_CLP
    //printf("***** experimental pretty crude barrier\n");
    //#define SAVEIT 2
#ifndef SAVEIT
#define BORROW
#endif
#ifdef BORROW
    ClpInterior barrier;
    barrier.borrowModel(*model2);
#else
    ClpInterior barrier(*model2);
#endif
    barrier.eventHandler()->setSimplex(NULL);
    if (interrupt)
      currentModel2 = &barrier;
    if (barrier.numberRows() + barrier.numberColumns() > 10000)
      barrier.setMaximumBarrierIterations(1000);
    int barrierOptions = options.getSpecialOption(4);
    int aggressiveGamma = 0;
    bool presolveInCrossover = false;
    bool scale = false;
    bool doKKT = false;
    bool forceFixing = false;
    int speed = 0;
    if (barrierOptions & 16) {
      barrierOptions &= ~16;
      doKKT = true;
    }
    if (barrierOptions & (32 + 64 + 128)) {
      aggressiveGamma = (barrierOptions & (32 + 64 + 128)) >> 5;
      barrierOptions &= ~(32 + 64 + 128);
    }
    if (barrierOptions & 256) {
      barrierOptions &= ~256;
      presolveInCrossover = true;
    }
    if (barrierOptions & 512) {
      barrierOptions &= ~512;
      forceFixing = true;
    }
    if (barrierOptions & 1024) {
      barrierOptions &= ~1024;
      barrier.setProjectionTolerance(1.0e-9);
    }
    if (barrierOptions & (2048 | 4096)) {
      speed = (barrierOptions & (2048 | 4096)) >> 11;
      barrierOptions &= ~(2048 | 4096);
    }
    if (barrierOptions & 8) {
      barrierOptions &= ~8;
      scale = true;
    }
    // If quadratic force KKT
    if (quadraticObj) {
      doKKT = true;
    }
    switch (barrierOptions) {
    case 0:
    default:
      if (!doKKT) {
        ClpCholeskyBase *cholesky = new ClpCholeskyBase(options.getExtraInfo(1));
        cholesky->setIntegerParameter(0, speed);
        barrier.setCholesky(cholesky);
      } else {
        ClpCholeskyBase *cholesky = new ClpCholeskyBase();
        cholesky->setKKT(true);
        barrier.setCholesky(cholesky);
      }
      break;
    case 1:
      if (!doKKT) {
        ClpCholeskyDense *cholesky = new ClpCholeskyDense();
        barrier.setCholesky(cholesky);
      } else {
        ClpCholeskyDense *cholesky = new ClpCholeskyDense();
        cholesky->setKKT(true);
        barrier.setCholesky(cholesky);
      }
      break;
#if defined(UFL_BARRIER) && (defined(CLP_HAS_AMD) || defined(CLP_HAS_CHOLMOD))
    case 4:
      if (!doKKT) {
        ClpCholeskyUfl *cholesky = new ClpCholeskyUfl(options.getExtraInfo(1));
        barrier.setCholesky(cholesky);
      } else {
        ClpCholeskyUfl *cholesky = new ClpCholeskyUfl();
        cholesky->setKKT(true);
        barrier.setCholesky(cholesky);
      }
      break;
#endif
#ifdef PARDISO_BARRIER
    case 7: {
      ClpCholeskyPardiso *cholesky = new ClpCholeskyPardiso();
      barrier.setCholesky(cholesky);
      assert(!doKKT);
    } break;
#endif
    }
    int numberRows = model2->numberRows();
    int numberColumns = model2->numberColumns();
    int saveMaxIts = model2->maximumIterations();
    if (saveMaxIts < 1000) {
      barrier.setMaximumBarrierIterations(saveMaxIts);
      model2->setMaximumIterations(10000000);
    }
#ifndef SAVEIT
    //barrier.setDiagonalPerturbation(1.0e-25);
    if (aggressiveGamma) {
      switch (aggressiveGamma) {
      case 1:
        barrier.setGamma(1.0e-5);
        barrier.setDelta(1.0e-5);
        break;
      case 2:
        barrier.setGamma(1.0e-7);
        break;
      case 3:
        barrier.setDelta(1.0e-5);
        break;
      case 4:
        barrier.setGamma(1.0e-3);
        barrier.setDelta(1.0e-3);
        break;
      case 5:
        barrier.setGamma(1.0e-3);
        break;
      case 6:
        barrier.setDelta(1.0e-3);
        break;
      }
    }
    if (scale)
      barrier.scaling(1);
    else
      barrier.scaling(0);
    barrier.primalDual();
#elif SAVEIT == 1
    barrier.primalDual();
#else
    model2->restoreModel("xx.save");
    // move solutions
    CoinMemcpyN(model2->primalRowSolution(),
      numberRows, barrier.primalRowSolution());
    CoinMemcpyN(model2->dualRowSolution(),
      numberRows, barrier.dualRowSolution());
    CoinMemcpyN(model2->primalColumnSolution(),
      numberColumns, barrier.primalColumnSolution());
    CoinMemcpyN(model2->dualColumnSolution(),
      numberColumns, barrier.dualColumnSolution());
#endif
    time2 = clpGetTime();
    timeCore = time2 - timeX;
    handler_->message(CLP_INTERVAL_TIMING, messages_)
      << "Barrier" << timeCore << time2 - time1
      << CoinMessageEol;
    timeX = time2;
    int maxIts = barrier.maximumBarrierIterations();
    int barrierStatus = barrier.status();
    double gap = barrier.complementarityGap();
    // get which variables are fixed
    double *saveLower = NULL;
    double *saveUpper = NULL;
    ClpPresolve pinfo2;
    ClpSimplex *saveModel2 = NULL;
    bool extraPresolve = false;
    int numberFixed = barrier.numberFixed();
    if (numberFixed) {
      int numberRows = barrier.numberRows();
      int numberColumns = barrier.numberColumns();
      int numberTotal = numberRows + numberColumns;
      saveLower = new double[numberTotal];
      saveUpper = new double[numberTotal];
      CoinMemcpyN(barrier.columnLower(), numberColumns, saveLower);
      CoinMemcpyN(barrier.rowLower(), numberRows, saveLower + numberColumns);
      CoinMemcpyN(barrier.columnUpper(), numberColumns, saveUpper);
      CoinMemcpyN(barrier.rowUpper(), numberRows, saveUpper + numberColumns);
    }
    if (((numberFixed * 20 > barrier.numberRows() && numberFixed > 5000) || forceFixing) && presolveInCrossover) {
      // may as well do presolve
      if (!forceFixing) {
        barrier.fixFixed();
      } else {
        // Fix
        int n = barrier.numberColumns();
        double *lower = barrier.columnLower();
        double *upper = barrier.columnUpper();
        double *solution = barrier.primalColumnSolution();
#ifdef CLP_INVESTIGATE
        int nFix = 0;
#endif
        for (int i = 0; i < n; i++) {
          if (barrier.fixedOrFree(i) && lower[i] < upper[i]) {
            double value = solution[i];
            if (value < lower[i] + 1.0e-6 && value - lower[i] < upper[i] - value) {
              solution[i] = lower[i];
              upper[i] = lower[i];
#ifdef CLP_INVESTIGATE
              nFix++;
#endif
            } else if (value > upper[i] - 1.0e-6 && value - lower[i] > upper[i] - value) {
              solution[i] = upper[i];
              lower[i] = upper[i];
#ifdef CLP_INVESTIGATE
              nFix++;
#endif
            }
          }
        }
#ifdef CLP_INVESTIGATE
        printf("%d columns fixed\n", nFix);
#endif
        int nr = barrier.numberRows();
        lower = barrier.rowLower();
        upper = barrier.rowUpper();
        solution = barrier.primalRowSolution();
#ifdef CLP_INVESTIGATE
        nFix = 0;
#endif
        for (int i = 0; i < nr; i++) {
          if (barrier.fixedOrFree(i + n) && lower[i] < upper[i]) {
            double value = solution[i];
            if (value < lower[i] + 1.0e-6 && value - lower[i] < upper[i] - value) {
              solution[i] = lower[i];
              upper[i] = lower[i];
#ifdef CLP_INVESTIGATE
              nFix++;
#endif
            } else if (value > upper[i] - 1.0e-6 && value - lower[i] > upper[i] - value) {
              solution[i] = upper[i];
              lower[i] = upper[i];
#ifdef CLP_INVESTIGATE
              nFix++;
#endif
            }
          }
        }
#ifdef CLP_INVESTIGATE
        printf("%d row slacks fixed\n", nFix);
#endif
      }
      saveModel2 = model2;
      extraPresolve = true;
    } else if (numberFixed) {
      // Set fixed to bounds (may have restored earlier solution)
      if (!forceFixing) {
        barrier.fixFixed(false);
      } else {
        // Fix
        int n = barrier.numberColumns();
        double *lower = barrier.columnLower();
        double *upper = barrier.columnUpper();
        double *solution = barrier.primalColumnSolution();
#ifdef CLP_INVESTIGATE
        int nFix = 0;
#endif
        for (int i = 0; i < n; i++) {
          if (barrier.fixedOrFree(i) && lower[i] < upper[i]) {
            double value = solution[i];
            if (value < lower[i] + 1.0e-8 && value - lower[i] < upper[i] - value) {
              solution[i] = lower[i];
              upper[i] = lower[i];
#ifdef CLP_INVESTIGATE
              nFix++;
#endif
            } else if (value > upper[i] - 1.0e-8 && value - lower[i] > upper[i] - value) {
              solution[i] = upper[i];
              lower[i] = upper[i];
#ifdef CLP_INVESTIGATE
              nFix++;
#endif
            } else {
              //printf("fixcol %d %g <= %g <= %g\n",
              //     i,lower[i],solution[i],upper[i]);
            }
          }
        }
#ifdef CLP_INVESTIGATE
        printf("%d columns fixed\n", nFix);
#endif
        int nr = barrier.numberRows();
        lower = barrier.rowLower();
        upper = barrier.rowUpper();
        solution = barrier.primalRowSolution();
#ifdef CLP_INVESTIGATE
        nFix = 0;
#endif
        for (int i = 0; i < nr; i++) {
          if (barrier.fixedOrFree(i + n) && lower[i] < upper[i]) {
            double value = solution[i];
            if (value < lower[i] + 1.0e-5 && value - lower[i] < upper[i] - value) {
              solution[i] = lower[i];
              upper[i] = lower[i];
#ifdef CLP_INVESTIGATE
              nFix++;
#endif
            } else if (value > upper[i] - 1.0e-5 && value - lower[i] > upper[i] - value) {
              solution[i] = upper[i];
              lower[i] = upper[i];
#ifdef CLP_INVESTIGATE
              nFix++;
#endif
            } else {
              //printf("fixrow %d %g <= %g <= %g\n",
              //     i,lower[i],solution[i],upper[i]);
            }
          }
        }
#ifdef CLP_INVESTIGATE
        printf("%d row slacks fixed\n", nFix);
#endif
      }
    }
#ifdef BORROW
    int saveNumberIterations = barrier.numberIterations();
    barrier.returnModel(*model2);
    double *rowPrimal = new double[numberRows];
    double *columnPrimal = new double[numberColumns];
    double *rowDual = new double[numberRows];
    double *columnDual = new double[numberColumns];
    // move solutions other way
    CoinMemcpyN(model2->primalRowSolution(),
      numberRows, rowPrimal);
    CoinMemcpyN(model2->dualRowSolution(),
      numberRows, rowDual);
    CoinMemcpyN(model2->primalColumnSolution(),
      numberColumns, columnPrimal);
    CoinMemcpyN(model2->dualColumnSolution(),
      numberColumns, columnDual);
#else
    double *rowPrimal = barrier.primalRowSolution();
    double *columnPrimal = barrier.primalColumnSolution();
    double *rowDual = barrier.dualRowSolution();
    double *columnDual = barrier.dualColumnSolution();
    // move solutions
    CoinMemcpyN(rowPrimal,
      numberRows, model2->primalRowSolution());
    CoinMemcpyN(rowDual,
      numberRows, model2->dualRowSolution());
    CoinMemcpyN(columnPrimal,
      numberColumns, model2->primalColumnSolution());
    CoinMemcpyN(columnDual,
      numberColumns, model2->dualColumnSolution());
#endif
    if (saveModel2) {
      // do presolve
      model2 = pinfo2.presolvedModel(*model2, dblParam_[ClpPresolveTolerance],
        false, 5, true);
      if (!model2) {
        model2 = saveModel2;
        saveModel2 = NULL;
        int numberRows = model2->numberRows();
        int numberColumns = model2->numberColumns();
        CoinMemcpyN(saveLower, numberColumns, model2->columnLower());
        CoinMemcpyN(saveLower + numberColumns, numberRows, model2->rowLower());
        delete[] saveLower;
        CoinMemcpyN(saveUpper, numberColumns, model2->columnUpper());
        CoinMemcpyN(saveUpper + numberColumns, numberRows, model2->rowUpper());
        delete[] saveUpper;
        saveLower = NULL;
        saveUpper = NULL;
      }
    }
    if (method == ClpSolve::useBarrier || barrierStatus < 0) {
      if (maxIts && barrierStatus < 4 && !quadraticObj) {
        //printf("***** crossover - needs more thought on difficult models\n");
#if SAVEIT == 1
        model2->ClpSimplex::saveModel("xx.save");
#endif
        // make sure no status left
        model2->createStatus();
        // solve
        if (!forceFixing)
          model2->setPerturbation(100);
        if (model2->factorizationFrequency() == 200) {
          // User did not touch preset
          model2->defaultFactorizationFrequency();
        }
  // throw some into basis
        if (!forceFixing) {
          int numberRows = model2->numberRows();
          int numberColumns = model2->numberColumns();
          double *dsort = new double[numberColumns];
          int *sort = new int[numberColumns];
          int n = 0;
          const double *columnLower = model2->columnLower();
          const double *columnUpper = model2->columnUpper();
          double *primalSolution = model2->primalColumnSolution();
          const double *dualSolution = model2->dualColumnSolution();
          double tolerance = 10.0 * primalTolerance_;
          int i;
          for (i = 0; i < numberRows; i++)
            model2->setRowStatus(i, superBasic);
          for (i = 0; i < numberColumns; i++) {
            double distance = std::min(columnUpper[i] - primalSolution[i],
              primalSolution[i] - columnLower[i]);
            if (distance > tolerance) {
              if (fabs(dualSolution[i]) < 1.0e-5)
                distance *= 100.0;
              dsort[n] = -distance;
              sort[n++] = i;
              model2->setStatus(i, superBasic);
            } else if (distance > primalTolerance_) {
              model2->setStatus(i, superBasic);
            } else if (primalSolution[i] <= columnLower[i] + primalTolerance_) {
              model2->setStatus(i, atLowerBound);
              primalSolution[i] = columnLower[i];
            } else {
              model2->setStatus(i, atUpperBound);
              primalSolution[i] = columnUpper[i];
            }
          }
          CoinSort_2(dsort, dsort + n, sort);
          n = std::min(numberRows, n);
          for (i = 0; i < n; i++) {
            int iColumn = sort[i];
            model2->setStatus(iColumn, basic);
          }
          delete[] sort;
          delete[] dsort;
          // model2->allSlackBasis();
          if (gap < 1.0e-3 * static_cast< double >(numberRows + numberColumns)) {
            if (saveUpper) {
              int numberRows = model2->numberRows();
              int numberColumns = model2->numberColumns();
              CoinMemcpyN(saveLower, numberColumns, model2->columnLower());
              CoinMemcpyN(saveLower + numberColumns, numberRows, model2->rowLower());
              CoinMemcpyN(saveUpper, numberColumns, model2->columnUpper());
              CoinMemcpyN(saveUpper + numberColumns, numberRows, model2->rowUpper());
              delete[] saveLower;
              delete[] saveUpper;
              saveLower = NULL;
              saveUpper = NULL;
            }
            //int numberRows = model2->numberRows();
            //int numberColumns = model2->numberColumns();
            // just primal values pass
            double saveScale = model2->objectiveScale();
            model2->setObjectiveScale(1.0e-3);
            model2->primal(2);
            model2->setObjectiveScale(saveScale);
            // save primal solution and copy back dual
            CoinMemcpyN(model2->primalRowSolution(),
              numberRows, rowPrimal);
            CoinMemcpyN(rowDual,
              numberRows, model2->dualRowSolution());
            CoinMemcpyN(model2->primalColumnSolution(),
              numberColumns, columnPrimal);
            CoinMemcpyN(columnDual,
              numberColumns, model2->dualColumnSolution());
            //model2->primal(1);
            // clean up reduced costs and flag variables
            {
              double *dj = model2->dualColumnSolution();
              double *cost = model2->objective();
              double *saveCost = new double[numberColumns];
              CoinMemcpyN(cost, numberColumns, saveCost);
              double *saveLower = new double[numberColumns];
              double *lower = model2->columnLower();
              CoinMemcpyN(lower, numberColumns, saveLower);
              double *saveUpper = new double[numberColumns];
              double *upper = model2->columnUpper();
              CoinMemcpyN(upper, numberColumns, saveUpper);
              int i;
              double tolerance = 10.0 * dualTolerance_;
              for (i = 0; i < numberColumns; i++) {
                if (model2->getStatus(i) == basic) {
                  dj[i] = 0.0;
                } else if (model2->getStatus(i) == atLowerBound) {
                  if (optimizationDirection_ * dj[i] < tolerance) {
                    if (optimizationDirection_ * dj[i] < 0.0) {
                      //if (dj[i]<-1.0e-3)
                      //printf("bad dj at lb %d %g\n",i,dj[i]);
                      cost[i] -= dj[i];
                      dj[i] = 0.0;
                    }
                  } else {
                    upper[i] = lower[i];
                  }
                } else if (model2->getStatus(i) == atUpperBound) {
                  if (optimizationDirection_ * dj[i] > tolerance) {
                    if (optimizationDirection_ * dj[i] > 0.0) {
                      //if (dj[i]>1.0e-3)
                      //printf("bad dj at ub %d %g\n",i,dj[i]);
                      cost[i] -= dj[i];
                      dj[i] = 0.0;
                    }
                  } else {
                    lower[i] = upper[i];
                  }
                }
              }
              // just dual values pass
              //model2->setLogLevel(63);
              //model2->setFactorizationFrequency(1);
              if (!gap)
                model2->dual(2);
              CoinMemcpyN(saveCost, numberColumns, cost);
              delete[] saveCost;
              CoinMemcpyN(saveLower, numberColumns, lower);
              delete[] saveLower;
              CoinMemcpyN(saveUpper, numberColumns, upper);
              delete[] saveUpper;
            }
          }
          // and finish
          // move solutions
          CoinMemcpyN(rowPrimal,
            numberRows, model2->primalRowSolution());
          CoinMemcpyN(columnPrimal,
            numberColumns, model2->primalColumnSolution());
        }
        double saveScale = model2->objectiveScale();
        model2->setObjectiveScale(1.0e-3);
        model2->primal(2);
        model2->setObjectiveScale(saveScale);
        model2->primal(1);
      } else if (barrierStatus == 4) {
        // memory problems
        model2->setPerturbation(savePerturbation);
        model2->createStatus();
        model2->dual();
      } else if (maxIts && quadraticObj) {
        // make sure no status left
        model2->createStatus();
        // solve
        model2->setPerturbation(100);
        model2->reducedGradient(1);
      }
    }

    //model2->setMaximumIterations(saveMaxIts);
#ifdef BORROW
    model2->setNumberIterations(model2->numberIterations() + saveNumberIterations);
    delete[] rowPrimal;
    delete[] columnPrimal;
    delete[] rowDual;
    delete[] columnDual;
#endif
    if (extraPresolve) {
      pinfo2.postsolve(true);
      delete model2;
      model2 = saveModel2;
    }
    if (saveUpper) {
      if (!forceFixing) {
        int numberRows = model2->numberRows();
        int numberColumns = model2->numberColumns();
        CoinMemcpyN(saveLower, numberColumns, model2->columnLower());
        CoinMemcpyN(saveLower + numberColumns, numberRows, model2->rowLower());
        CoinMemcpyN(saveUpper, numberColumns, model2->columnUpper());
        CoinMemcpyN(saveUpper + numberColumns, numberRows, model2->rowUpper());
      }
      delete[] saveLower;
      delete[] saveUpper;
      saveLower = NULL;
      saveUpper = NULL;
      if (method != ClpSolve::useBarrierNoCross)
        model2->primal(1);
    }
    model2->setPerturbation(savePerturbation);
    time2 = clpGetTime();
    timeCore = time2 - timeX;
    handler_->message(CLP_INTERVAL_TIMING, messages_)
      << "Crossover" << timeCore << time2 - time1
      << CoinMessageEol;
    timeX = time2;
#else
    abort();
#endif
  } else if (method == ClpSolve::notImplemented) {
    printf("done decomposition\n");
  } else {
    assert(method != ClpSolve::automatic); // later
    time2 = 0.0;
  }
  if (saveMatrix) {
    if (model2 == this) {
      // delete and replace
      delete model2->clpMatrix();
      model2->replaceMatrix(saveMatrix);
    } else {
      delete saveMatrix;
    }
  }
  numberIterations = model2->numberIterations();
  finalStatus = model2->status();
  int finalSecondaryStatus = model2->secondaryStatus();
  if (presolve == ClpSolve::presolveOn) {
    int saveLevel = logLevel();
    if ((specialOptions_ & 1024) == 0)
      setLogLevel(std::min(1, saveLevel));
    else
      setLogLevel(std::min(0, saveLevel));
    pinfo->postsolve(true);
    numberIterations_ = 0;
    delete pinfo;
    pinfo = NULL;
    factorization_->areaFactor(model2->factorization()->adjustedAreaFactor());
    time2 = clpGetTime();
    timePresolve += time2 - timeX;
    handler_->message(CLP_INTERVAL_TIMING, messages_)
      << "Postsolve" << time2 - timeX << time2 - time1
      << CoinMessageEol;
    timeX = time2;
    if (!presolveToFile) {
      delete model2;
    }
    if (interrupt)
      currentModel = this;
    // checkSolution(); already done by postSolve
    setLogLevel(saveLevel);
    int oldStatus = problemStatus_;
    setProblemStatus(finalStatus);
    setSecondaryStatus(finalSecondaryStatus);
    /* Code modified so rcode -1 as normal, >0 clean up, 0 say optimal */
    int rcode = eventHandler()->event(ClpEventHandler::presolveAfterFirstSolve);
    //#define TREAT_AS_OPTIMAL_TOLERANCE 1.0e-4
#ifdef TREAT_AS_OPTIMAL_TOLERANCE
    if (rcode == -1 && sumPrimalInfeasibilities_ < TREAT_AS_OPTIMAL_TOLERANCE && sumDualInfeasibilities_ < TREAT_AS_OPTIMAL_TOLERANCE)
      rcode = 0;
#endif
    if (finalStatus != 3 && rcode < 0 && (finalStatus || oldStatus == -1)) {
      double sumPrimal = sumPrimalInfeasibilities_;
      double sumDual = sumDualInfeasibilities_;
      if (sumDual > 1.0e-6 && sumPrimal > 1.0e-6)
        moreSpecialOptions_ &= ~2; // be safe and do final solve
      // ignore some parts of solution
      if (finalStatus == 1) {
        // infeasible
        sumDual = 0.0;
      } else if (finalStatus == 2) {
        sumPrimal = 0.0;
      }
      int savePerturbation = perturbation();
      if (savePerturbation == 50)
        setPerturbation(51); // small
      if ((finalStatus>=0 && finalStatus <= 2) ||
	  (moreSpecialOptions_ & 2) == 0 ||
	  fabs(sumDual) + fabs(sumPrimal) < 1.0e-3) {
        if (finalStatus == 2) {
          if (sumDual > 1.0e-4) {
            // unbounded - get feasible first
            double save = optimizationDirection_;
            optimizationDirection_ = 0.0;
            primal(1);
            optimizationDirection_ = save;
          }
          primal(1);
        } else if (finalStatus == 1) {
          dual();
        } else {
          if ((moreSpecialOptions_ & 65536) == 0) {
            if (numberRows_ < 10000 || true)
              setPerturbation(100); // probably better to perturb after n its
            else if (savePerturbation < 100)
              setPerturbation(51); // probably better to perturb after n its
          }
          // use method thought suitable
          int numberSuperBasic = 0;
          for (int i = 0; i < numberColumns_; i++) {
            if (getColumnStatus(i) == superBasic)
              numberSuperBasic++;
          }
	  // double check if looks odd
          if (sumDual > 1000.0 * sumPrimal || numberSuperBasic) {
            primal(1);
	    if (!finalStatus&&problemStatus_) {
	      handler_->message(CLP_GENERAL, messages_)
		<< "Primal after postsolve not optimal! - trying dual"
		<< CoinMessageEol;
	      dual();
	    }
          } else if (sumPrimal > 1000.0 * sumDual) {
            dual();
	    if (!finalStatus&&problemStatus_) {
	      handler_->message(CLP_GENERAL, messages_)
		<< "Dual after postsolve not optimal! - trying primal"
		<< CoinMessageEol;
	      primal(1);
	    }
          } else {
            if (method != ClpSolve::useDual) {
              primal(1);
	      if (!finalStatus&&problemStatus_) {
		handler_->message(CLP_GENERAL, messages_)
		  << "Primal after postsolve not optimal! - trying dual"
		  << CoinMessageEol;
		dual();
	      }
            } else {
              dual();
	      if (!finalStatus&&problemStatus_) {
		handler_->message(CLP_GENERAL, messages_)
		  << "Dual after postsolve not optimal! - trying primal"
		  << CoinMessageEol;
		primal(1);
	      }
	    }
          }
        }
      } else {
        // just set status
        problemStatus_ = finalStatus;
      }
      setPerturbation(savePerturbation);
      numberIterations += numberIterations_;
      numberIterations_ = numberIterations;
      finalStatus = status();
      time2 = clpGetTime();
      handler_->message(CLP_INTERVAL_TIMING, messages_)
        << "Cleanup" << time2 - timeX << time2 - time1
        << CoinMessageEol;
      timeX = time2;
    } else if (rcode > 0) { // was >= 0
      primal(1);
    } else {
      secondaryStatus_ = finalSecondaryStatus;
    }
  } else if (model2 != this) {
    // not presolved - but different model used (sprint probably)
    CoinMemcpyN(model2->primalRowSolution(),
      numberRows_, this->primalRowSolution());
    CoinMemcpyN(model2->dualRowSolution(),
      numberRows_, this->dualRowSolution());
    CoinMemcpyN(model2->primalColumnSolution(),
      numberColumns_, this->primalColumnSolution());
    CoinMemcpyN(model2->dualColumnSolution(),
      numberColumns_, this->dualColumnSolution());
    CoinMemcpyN(model2->statusArray(),
      numberColumns_ + numberRows_, this->statusArray());
    objectiveValue_ = model2->objectiveValue_;
    numberIterations_ = model2->numberIterations_;
    problemStatus_ = model2->problemStatus_;
    secondaryStatus_ = model2->secondaryStatus_;
    delete model2;
  }
  if (method != ClpSolve::useBarrierNoCross && method != ClpSolve::useBarrier)
    setMaximumIterations(saveMaxIterations);
  std::string statusMessage[] = { "Unknown", "Optimal", "PrimalInfeasible", "DualInfeasible", "Stopped",
    "Errors", "User stopped" };
  assert(finalStatus >= -1 && finalStatus <= 5);
  numberIterations_ = numberIterations;
  handler_->message(CLP_TIMING, messages_)
    << statusMessage[finalStatus + 1] << objectiveValue() << numberIterations << time2 - time1;
  handler_->printing(presolve == ClpSolve::presolveOn)
    << timePresolve;
  handler_->printing(timeIdiot != 0.0)
    << timeIdiot;
  handler_->message() << CoinMessageEol;
  if (finalStatus==0) {
    if (secondaryStatus_>1 && secondaryStatus_<5) {
      std::string scaledMessage[] = {
	    "Unscaled problem has primal infeasibilities",
	    "Unscaled problem has dual infeasibilities",
	    "Unscaled problem has primal and dual infeasibilities"};
      handler_->message(CLP_GENERAL,messages_)
	<< scaledMessage[secondaryStatus_-2];
      handler_->message() << CoinMessageEol;
      if ((moreSpecialOptions_&134217728)!=0) {
	// solve without scaling
	scalingFlag_ = 0;
	delete[] rowScale_;
	delete[] columnScale_;
	rowScale_ = NULL;
	columnScale_ = NULL;
	inverseRowScale_ = NULL;
	inverseColumnScale_ = NULL;
	primal(1);
      }
    }
  }
  if (interrupt)
    signal(SIGINT, saveSignal);
  perturbation_ = savePerturbation;
  scalingFlag_ = saveScaling;
  // If faking objective - put back correct one
  if (savedObjective) {
    delete objective_;
    objective_ = savedObjective;
  }
  if (options.getSpecialOption(1) == 2 && options.getExtraInfo(1) > 1000000
      && options.getExtraInfo(1) < 2000000) {
    ClpObjective *savedObjective = objective_;
    // make up zero objective
    double *obj = new double[numberColumns_];
    for (int i = 0; i < numberColumns_; i++)
      obj[i] = 0.0;
    objective_ = new ClpLinearObjective(obj, numberColumns_);
    delete[] obj;
    primal(1);
    delete objective_;
    objective_ = savedObjective;
    finalStatus = status();
  }
  eventHandler()->event(ClpEventHandler::presolveEnd);
  delete pinfo;
  moreSpecialOptions_ = saveMoreOptions;
#ifdef CLP_USEFUL_PRINTOUT
  debugInt[23] = numberIterations_;
#endif
  return finalStatus;
}
// General solve
int ClpSimplex::initialSolve()
{
  // Default so use dual
  ClpSolve options;
  return initialSolve(options);
}
// General dual solve
int ClpSimplex::initialDualSolve()
{
  ClpSolve options;
  // Use dual
  options.setSolveType(ClpSolve::useDual);
  return initialSolve(options);
}
// General primal solve
int ClpSimplex::initialPrimalSolve()
{
  ClpSolve options;
  // Use primal
  options.setSolveType(ClpSolve::usePrimal);
  return initialSolve(options);
}
// barrier solve, not to be followed by crossover
int ClpSimplex::initialBarrierNoCrossSolve()
{
  ClpSolve options;
  // Use primal
  options.setSolveType(ClpSolve::useBarrierNoCross);
  int returnCode = initialSolve(options);
  // clean for simplex and put slacks in basis
  for (int i=0;i<numberRows_;i++)
    status_[i+numberColumns_] = ClpSimplex::basic;
  for (int i=0;i<numberColumns_;i++) {
    if (columnLower_[i]==columnUpper_[i]) {
      columnActivity_[i] = columnLower_[i];
      status_[i] = ClpSimplex::isFixed;
    } else if (columnActivity_[i]<=columnLower_[i]) {
      columnActivity_[i] = columnLower_[i];
      status_[i] = ClpSimplex::atLowerBound;
    } else if (columnActivity_[i]>=columnUpper_[i]) {
      columnActivity_[i] = columnUpper_[i];
      status_[i] = ClpSimplex::atUpperBound;
    } else {
      status_[i] = ClpSimplex::superBasic;
    }
  }
  return returnCode;
}

// General barrier solve
int ClpSimplex::initialBarrierSolve()
{
  ClpSolve options;
  // Use primal
  options.setSolveType(ClpSolve::useBarrier);
  return initialSolve(options);
}

// Default constructor
ClpSolve::ClpSolve()
{
  method_ = automatic;
  presolveType_ = presolveOn;
  numberPasses_ = 5;
  int i;
  for (i = 0; i < 7; i++)
    options_[i] = 0;
  // say no +-1 matrix
  options_[3] = 1;
  for (i = 0; i < 7; i++)
    extraInfo_[i] = -1;
  independentOptions_[0] = 0;
  // But switch off slacks
  independentOptions_[1] = 512;
  // Substitute up to 3
  independentOptions_[2] = 3;
}
// Constructor when you really know what you are doing
ClpSolve::ClpSolve(SolveType method, PresolveType presolveType,
  int numberPasses, int options[6],
  int extraInfo[6], int independentOptions[3])
{
  method_ = method;
  presolveType_ = presolveType;
  numberPasses_ = numberPasses;
  int i;
  for (i = 0; i < 6; i++)
    options_[i] = options[i];
  options_[6] = 0;
  for (i = 0; i < 6; i++)
    extraInfo_[i] = extraInfo[i];
  extraInfo_[6] = 0;
  for (i = 0; i < 3; i++)
    independentOptions_[i] = independentOptions[i];
}

// Copy constructor.
ClpSolve::ClpSolve(const ClpSolve &rhs)
{
  method_ = rhs.method_;
  presolveType_ = rhs.presolveType_;
  numberPasses_ = rhs.numberPasses_;
  int i;
  for (i = 0; i < 7; i++)
    options_[i] = rhs.options_[i];
  for (i = 0; i < 7; i++)
    extraInfo_[i] = rhs.extraInfo_[i];
  for (i = 0; i < 3; i++)
    independentOptions_[i] = rhs.independentOptions_[i];
}
// Assignment operator. This copies the data
ClpSolve &
ClpSolve::operator=(const ClpSolve &rhs)
{
  if (this != &rhs) {
    method_ = rhs.method_;
    presolveType_ = rhs.presolveType_;
    numberPasses_ = rhs.numberPasses_;
    int i;
    for (i = 0; i < 7; i++)
      options_[i] = rhs.options_[i];
    for (i = 0; i < 7; i++)
      extraInfo_[i] = rhs.extraInfo_[i];
    for (i = 0; i < 3; i++)
      independentOptions_[i] = rhs.independentOptions_[i];
  }
  return *this;
}
// Destructor
ClpSolve::~ClpSolve()
{
}
// See header file for details
void ClpSolve::setSpecialOption(int which, int value, int extraInfo)
{
  options_[which] = value;
  extraInfo_[which] = extraInfo;
}
int ClpSolve::getSpecialOption(int which) const
{
  return options_[which];
}

// Solve types
void ClpSolve::setSolveType(SolveType method, int /*extraInfo*/)
{
  method_ = method;
}

ClpSolve::SolveType
ClpSolve::getSolveType()
{
  return method_;
}

// Presolve types
void ClpSolve::setPresolveType(PresolveType amount, int extraInfo)
{
  presolveType_ = amount;
  numberPasses_ = extraInfo;
}
ClpSolve::PresolveType
ClpSolve::getPresolveType()
{
  return presolveType_;
}
// Extra info for idiot (or sprint)
int ClpSolve::getExtraInfo(int which) const
{
  return extraInfo_[which];
}
int ClpSolve::getPresolvePasses() const
{
  return numberPasses_;
}
/* Say to return at once if infeasible,
   default is to solve */
void ClpSolve::setInfeasibleReturn(bool trueFalse)
{
  independentOptions_[0] = trueFalse ? 1 : 0;
}
#include <string>
// Generates code for above constructor
void ClpSolve::generateCpp(FILE *fp)
{
  std::string solveType[] = {
    "ClpSolve::useDual",
    "ClpSolve::usePrimal",
    "ClpSolve::usePrimalorSprint",
    "ClpSolve::useBarrier",
    "ClpSolve::useBarrierNoCross",
    "ClpSolve::automatic",
    "ClpSolve::notImplemented"
  };
  std::string presolveType[] = {
    "ClpSolve::presolveOn",
    "ClpSolve::presolveOff",
    "ClpSolve::presolveNumber",
    "ClpSolve::presolveNumberCost"
  };
  fprintf(fp, "3  ClpSolve::SolveType method = %s;\n", solveType[method_].c_str());
  fprintf(fp, "3  ClpSolve::PresolveType presolveType = %s;\n",
    presolveType[presolveType_].c_str());
  fprintf(fp, "3  int numberPasses = %d;\n", numberPasses_);
  fprintf(fp, "3  int options[] = {%d,%d,%d,%d,%d,%d};\n",
    options_[0], options_[1], options_[2],
    options_[3], options_[4], options_[5]);
  fprintf(fp, "3  int extraInfo[] = {%d,%d,%d,%d,%d,%d};\n",
    extraInfo_[0], extraInfo_[1], extraInfo_[2],
    extraInfo_[3], extraInfo_[4], extraInfo_[5]);
  fprintf(fp, "3  int independentOptions[] = {%d,%d,%d};\n",
    independentOptions_[0], independentOptions_[1], independentOptions_[2]);
  fprintf(fp, "3  ClpSolve clpSolve(method,presolveType,numberPasses,\n");
  fprintf(fp, "3                    options,extraInfo,independentOptions);\n");
}
//#############################################################################
#include "ClpNonLinearCost.hpp"

ClpSimplexProgress::ClpSimplexProgress()
{
  int i;
  for (i = 0; i < CLP_PROGRESS; i++) {
    objective_[i] = COIN_DBL_MAX * 1.0e-50;
    infeasibility_[i] = -1.0; // set to an impossible value
    realInfeasibility_[i] = COIN_DBL_MAX * 1.0e-50;
    numberInfeasibilities_[i] = -1;
    iterationNumber_[i] = -1;
  }
#ifdef CLP_PROGRESS_WEIGHT
  for (i = 0; i < CLP_PROGRESS_WEIGHT; i++) {
    objectiveWeight_[i] = COIN_DBL_MAX * 1.0e-50;
    infeasibilityWeight_[i] = -1.0; // set to an impossible value
    realInfeasibilityWeight_[i] = COIN_DBL_MAX * 1.0e-50;
    numberInfeasibilitiesWeight_[i] = -1;
    iterationNumberWeight_[i] = -1;
  }
  drop_ = 0.0;
  best_ = 0.0;
#endif
  initialWeight_ = 0.0;
  for (i = 0; i < CLP_CYCLE; i++) {
    //obj_[i]=COIN_DBL_MAX*1.0e-50;
    in_[i] = -1;
    out_[i] = -1;
    way_[i] = 0;
  }
  numberTimes_ = 0;
  numberBadTimes_ = 0;
  numberReallyBadTimes_ = 0;
  numberTimesFlagged_ = 0;
  model_ = NULL;
  oddState_ = 0;
  checkScalingAfter_ = 0;
}

//-----------------------------------------------------------------------------

ClpSimplexProgress::~ClpSimplexProgress()
{
}
// Copy constructor.
ClpSimplexProgress::ClpSimplexProgress(const ClpSimplexProgress &rhs)
{
  int i;
  for (i = 0; i < CLP_PROGRESS; i++) {
    objective_[i] = rhs.objective_[i];
    infeasibility_[i] = rhs.infeasibility_[i];
    realInfeasibility_[i] = rhs.realInfeasibility_[i];
    numberInfeasibilities_[i] = rhs.numberInfeasibilities_[i];
    iterationNumber_[i] = rhs.iterationNumber_[i];
  }
#ifdef CLP_PROGRESS_WEIGHT
  for (i = 0; i < CLP_PROGRESS_WEIGHT; i++) {
    objectiveWeight_[i] = rhs.objectiveWeight_[i];
    infeasibilityWeight_[i] = rhs.infeasibilityWeight_[i];
    realInfeasibilityWeight_[i] = rhs.realInfeasibilityWeight_[i];
    numberInfeasibilitiesWeight_[i] = rhs.numberInfeasibilitiesWeight_[i];
    iterationNumberWeight_[i] = rhs.iterationNumberWeight_[i];
  }
  drop_ = rhs.drop_;
  best_ = rhs.best_;
#endif
  initialWeight_ = rhs.initialWeight_;
  for (i = 0; i < CLP_CYCLE; i++) {
    //obj_[i]=rhs.obj_[i];
    in_[i] = rhs.in_[i];
    out_[i] = rhs.out_[i];
    way_[i] = rhs.way_[i];
  }
  numberTimes_ = rhs.numberTimes_;
  numberBadTimes_ = rhs.numberBadTimes_;
  numberReallyBadTimes_ = rhs.numberReallyBadTimes_;
  numberTimesFlagged_ = rhs.numberTimesFlagged_;
  model_ = rhs.model_;
  oddState_ = rhs.oddState_;
  checkScalingAfter_ = rhs.checkScalingAfter_;
}
// Copy constructor.from model
ClpSimplexProgress::ClpSimplexProgress(ClpSimplex *model)
{
  model_ = model;
  reset();
  initialWeight_ = 0.0;
}
// Fill from model
void ClpSimplexProgress::fillFromModel(ClpSimplex *model)
{
  model_ = model;
  reset();
  initialWeight_ = 0.0;
}
// Assignment operator. This copies the data
ClpSimplexProgress &
ClpSimplexProgress::operator=(const ClpSimplexProgress &rhs)
{
  if (this != &rhs) {
    int i;
    for (i = 0; i < CLP_PROGRESS; i++) {
      objective_[i] = rhs.objective_[i];
      infeasibility_[i] = rhs.infeasibility_[i];
      realInfeasibility_[i] = rhs.realInfeasibility_[i];
      numberInfeasibilities_[i] = rhs.numberInfeasibilities_[i];
      iterationNumber_[i] = rhs.iterationNumber_[i];
    }
#ifdef CLP_PROGRESS_WEIGHT
    for (i = 0; i < CLP_PROGRESS_WEIGHT; i++) {
      objectiveWeight_[i] = rhs.objectiveWeight_[i];
      infeasibilityWeight_[i] = rhs.infeasibilityWeight_[i];
      realInfeasibilityWeight_[i] = rhs.realInfeasibilityWeight_[i];
      numberInfeasibilitiesWeight_[i] = rhs.numberInfeasibilitiesWeight_[i];
      iterationNumberWeight_[i] = rhs.iterationNumberWeight_[i];
    }
    drop_ = rhs.drop_;
    best_ = rhs.best_;
#endif
    initialWeight_ = rhs.initialWeight_;
    for (i = 0; i < CLP_CYCLE; i++) {
      //obj_[i]=rhs.obj_[i];
      in_[i] = rhs.in_[i];
      out_[i] = rhs.out_[i];
      way_[i] = rhs.way_[i];
    }
    numberTimes_ = rhs.numberTimes_;
    numberBadTimes_ = rhs.numberBadTimes_;
    numberReallyBadTimes_ = rhs.numberReallyBadTimes_;
    numberTimesFlagged_ = rhs.numberTimesFlagged_;
    model_ = rhs.model_;
    oddState_ = rhs.oddState_;
    checkScalingAfter_ = rhs.checkScalingAfter_;
  }
  return *this;
}
// Seems to be something odd about exact comparison of doubles on linux
static bool equalDouble(double value1, double value2)
{

  union {
    double d;
    int i[2];
  } v1, v2;
  v1.d = value1;
  v2.d = value2;
  if (sizeof(int) * 2 == sizeof(double))
    return (v1.i[0] == v2.i[0] && v1.i[1] == v2.i[1]);
  else
    return (v1.i[0] == v2.i[0]);
}
int ClpSimplexProgress::looping()
{
  if (!model_)
    return -1;
  double objective;
  if (model_->algorithm() < 0) {
    objective = model_->rawObjectiveValue();
    objective -= model_->bestPossibleImprovement();
  } else {
    objective = model_->nonLinearCost()->feasibleReportCost();
  }
  double infeasibility;
  double realInfeasibility = 0.0;
  int numberInfeasibilities;
  int iterationNumber = model_->numberIterations();
  //numberTimesFlagged_ = 0;
  if (model_->algorithm() < 0) {
    // dual
    infeasibility = model_->sumPrimalInfeasibilities();
    numberInfeasibilities = model_->numberPrimalInfeasibilities();
  } else {
    //primal
    infeasibility = model_->sumDualInfeasibilities();
    realInfeasibility = model_->nonLinearCost()->sumInfeasibilities();
    numberInfeasibilities = model_->numberDualInfeasibilities();
    if (iterationNumber>3*model_->numberRows()+3*model_->numberColumns()) {
      // should I put out a message
      return 1;
    }
  }
  int i;
  int numberMatched = 0;
  int matched = 0;
  int nsame = 0;
  for (i = 0; i < CLP_PROGRESS; i++) {
    bool matchedOnObjective = equalDouble(objective, objective_[i]);
    bool matchedOnInfeasibility = equalDouble(infeasibility, infeasibility_[i]);
    bool matchedOnInfeasibilities = (numberInfeasibilities == numberInfeasibilities_[i]);

    if (matchedOnObjective && matchedOnInfeasibility && matchedOnInfeasibilities) {
      matched |= (1 << i);
      // Check not same iteration
      if (iterationNumber != iterationNumber_[i]) {
        numberMatched++;
        // here mainly to get over compiler bug?
        if (model_->messageHandler()->logLevel() > 10)
          printf("%d %d %d %d %d loop check\n", i, numberMatched,
            matchedOnObjective, matchedOnInfeasibility,
            matchedOnInfeasibilities);
      } else {
        // stuck but code should notice
        nsame++;
      }
    }
    if (i) {
      objective_[i - 1] = objective_[i];
      infeasibility_[i - 1] = infeasibility_[i];
      realInfeasibility_[i - 1] = realInfeasibility_[i];
      numberInfeasibilities_[i - 1] = numberInfeasibilities_[i];
      iterationNumber_[i - 1] = iterationNumber_[i];
    }
  }
  objective_[CLP_PROGRESS - 1] = objective;
  infeasibility_[CLP_PROGRESS - 1] = infeasibility;
  realInfeasibility_[CLP_PROGRESS - 1] = realInfeasibility;
  numberInfeasibilities_[CLP_PROGRESS - 1] = numberInfeasibilities;
  iterationNumber_[CLP_PROGRESS - 1] = iterationNumber;
  if (nsame == CLP_PROGRESS)
    numberMatched = CLP_PROGRESS; // really stuck
  if (model_->progressFlag())
    numberMatched = 0;
  numberTimes_++;
  if (numberTimes_ < 10)
    numberMatched = 0;
  // skip if just last time as may be checking something
  if (matched == (1 << (CLP_PROGRESS - 1)))
    numberMatched = 0;
  if (model_->numberIterations()>20*model_->numberRows()
      +5*model_->numberColumns()+100 && (model_->specialOptions()&0x03000000)!=0) {
    // pretty bad
    // make factorize more often
    if (model_->numberIterations()<25*model_->numberRows()
	+8*model_->numberColumns()+300 && numberReallyBadTimes_<100) {
      model_->forceFactorization(std::min(model_->forceFactorization(),5));
      numberReallyBadTimes_++;
    } else {
      // give up
      numberMatched = 1000;
      numberBadTimes_ = 100;
    }
  }
  if (numberMatched && model_->clpMatrix()->type() < 15) {
    model_->messageHandler()->message(CLP_POSSIBLELOOP, model_->messages())
      << numberMatched
      << matched
      << numberTimes_
      << CoinMessageEol;
    numberBadTimes_++;
    if (numberBadTimes_ < 10) {
      // make factorize every iteration
      model_->forceFactorization(1);
      if (numberBadTimes_ < 2) {
        startCheck(); // clear other loop check
        if (model_->algorithm() < 0) {
          // dual - change tolerance
          model_->setCurrentDualTolerance(model_->currentDualTolerance() * 1.05);
          // if infeasible increase dual bound
          if (model_->dualBound() < 1.0e17) {
            model_->setDualBound(model_->dualBound() * 1.1);
            static_cast< ClpSimplexDual * >(model_)->resetFakeBounds(0);
          }
        } else {
          // primal - change tolerance
          if (numberBadTimes_ > 3)
            model_->setCurrentPrimalTolerance(model_->currentPrimalTolerance() * 1.05);
          // if infeasible increase infeasibility cost
          if (model_->nonLinearCost()->numberInfeasibilities() && model_->infeasibilityCost() < 1.0e17) {
            model_->setInfeasibilityCost(model_->infeasibilityCost() * 1.1);
          }
        }
      } else {
        // flag
        int iSequence;
        if (model_->algorithm() < 0) {
          // dual
          if (model_->dualBound() > 1.0e14)
            model_->setDualBound(1.0e14);
          iSequence = in_[CLP_CYCLE - 1];
        } else {
          // primal
          //if (model_->infeasibilityCost() > 1.0e14)
          //   model_->setInfeasibilityCost(1.0e14);
          iSequence = out_[CLP_CYCLE - 1];
        }
        if (iSequence >= 0) {
          char x = model_->isColumn(iSequence) ? 'C' : 'R';
          if (model_->messageHandler()->logLevel() >= 63)
            model_->messageHandler()->message(CLP_SIMPLEX_FLAG, model_->messages())
              << x << model_->sequenceWithin(iSequence)
              << CoinMessageEol;
          // if Gub then needs to be sequenceIn_
          int save = model_->sequenceIn();
          model_->setSequenceIn(iSequence);
          model_->setFlagged(iSequence);
          model_->setSequenceIn(save);
          //printf("flagging %d from loop\n",iSequence);
          startCheck();
        } else {
          // Give up
          if (model_->messageHandler()->logLevel() >= 63)
            printf("***** All flagged?\n");
          return 4;
        }
        // reset
        numberBadTimes_ = 2;
      }
      return -2;
    } else {
      // look at solution and maybe declare victory
      if (infeasibility < 1.0e-4) {
        return 0;
      } else {
        model_->messageHandler()->message(CLP_LOOP, model_->messages())
          << CoinMessageEol;
#ifndef NDEBUG
        printf("debug loop ClpSimplex A\n");
        abort();
#endif
        return 3;
      }
    }
  }
  return -1;
}
// Resets as much as possible
void ClpSimplexProgress::reset()
{
  int i;
  for (i = 0; i < CLP_PROGRESS; i++) {
    if (model_->algorithm() >= 0)
      objective_[i] = COIN_DBL_MAX * 1.0e-50;
    else
      objective_[i] = -COIN_DBL_MAX * 1.0e-50;
    infeasibility_[i] = -1.0; // set to an impossible value
    realInfeasibility_[i] = COIN_DBL_MAX * 1.0e-50;
    numberInfeasibilities_[i] = -1;
    iterationNumber_[i] = -1;
  }
#ifdef CLP_PROGRESS_WEIGHT
  for (i = 0; i < CLP_PROGRESS_WEIGHT; i++) {
    objectiveWeight_[i] = COIN_DBL_MAX * 1.0e-50;
    infeasibilityWeight_[i] = -1.0; // set to an impossible value
    realInfeasibilityWeight_[i] = COIN_DBL_MAX * 1.0e-50;
    numberInfeasibilitiesWeight_[i] = -1;
    iterationNumberWeight_[i] = -1;
  }
  drop_ = 0.0;
  best_ = 0.0;
#endif
  for (i = 0; i < CLP_CYCLE; i++) {
    //obj_[i]=COIN_DBL_MAX*1.0e-50;
    in_[i] = -1;
    out_[i] = -1;
    way_[i] = 0;
  }
  numberTimes_ = 0;
  numberBadTimes_ = 0;
  numberReallyBadTimes_ = 0;
  numberTimesFlagged_ = 0;
  oddState_ = 0;
  checkScalingAfter_ = 0;
}
#if CLP_CHECK_SCALING
// Checks if all going well - may rescale
int
ClpSimplexProgress::checkScalingEtc()
{
  if (model_->numberIterations()<checkScalingAfter_||
      (model_->specialOptions() & 0x03000000) != 0)
    return 0;
  if (model_->numberIterations()) {
    checkScalingAfter_ = model_->checkScaling();
    if (checkScalingAfter_!=COIN_INT_MAX) {
      for (int i = 0; i < CLP_PROGRESS; i++) {
	if (model_->algorithm() >= 0)
	  objective_[i] = COIN_DBL_MAX * 1.0e-50;
	else
	  objective_[i] = -COIN_DBL_MAX * 1.0e-50;
      }
    }
    return 1;
  } else {
    return 0;
  }
}
#endif
// Returns previous objective (if -1) - current if (0)
double
ClpSimplexProgress::lastObjective(int back) const
{
  return objective_[CLP_PROGRESS - 1 - back];
}
// Returns previous infeasibility (if -1) - current if (0)
double
ClpSimplexProgress::lastInfeasibility(int back) const
{
  return realInfeasibility_[CLP_PROGRESS - 1 - back];
}
// Sets real primal infeasibility
void ClpSimplexProgress::setInfeasibility(double value)
{
  for (int i = 1; i < CLP_PROGRESS; i++)
    realInfeasibility_[i - 1] = realInfeasibility_[i];
  realInfeasibility_[CLP_PROGRESS - 1] = value;
}
// Returns number of primal infeasibilities (if +1) - current if (0)
int ClpSimplexProgress::numberInfeasibilities(int back) const
{
  return numberInfeasibilities_[CLP_PROGRESS - 1 - back];
}
// Modify objective e.g. if dual infeasible in dual
void ClpSimplexProgress::modifyObjective(double value)
{
  objective_[CLP_PROGRESS - 1] = value;
}
// Returns previous iteration number (if +1) - current if (0)
int ClpSimplexProgress::lastIterationNumber(int back) const
{
  return iterationNumber_[CLP_PROGRESS - 1 - back];
}
// clears iteration numbers (to switch off panic)
void ClpSimplexProgress::clearIterationNumbers()
{
  for (int i = 0; i < CLP_PROGRESS; i++)
    iterationNumber_[i] = -1;
}
// Start check at beginning of whileIterating
void ClpSimplexProgress::startCheck()
{
  int i;
  for (i = 0; i < CLP_CYCLE; i++) {
    //obj_[i]=COIN_DBL_MAX*1.0e-50;
    in_[i] = -1;
    out_[i] = -1;
    way_[i] = 0;
  }
}
// Returns cycle length in whileIterating
int ClpSimplexProgress::cycle(int in, int out, int wayIn, int wayOut)
{
  int i;
  int matched = 0;
  // first see if in matches any out
  for (i = 1; i < CLP_CYCLE; i++) {
    if (in == out_[i]) {
      // even if flip then suspicious
      matched = -1;
      break;
    }
  }
  if (matched && in_[0] >= 0) {
    // possible cycle - only check [0] against all
    matched = 0;
    int nMatched = 0;
    char way0 = way_[0];
    int in0 = in_[0];
    int out0 = out_[0];
    //double obj0 = obj_[i];
    for (int k = 1; k < CLP_CYCLE - 4; k++) {
      if (in0 == in_[k] && out0 == out_[k] && way0 == way_[k]) {
        nMatched++;
        // See if repeats
        int end = CLP_CYCLE - k;
        int j;
        for (j = 1; j < end; j++) {
          if (in_[j + k] != in_[j] || out_[j + k] != out_[j] || way_[j + k] != way_[j])
            break;
        }
        if (j == end) {
          matched = k;
          break;
        }
      }
    }
    // If three times then that is too much even if not regular
    if (matched <= 0 && nMatched > 1)
      matched = 100;
  }
  for (i = 0; i < CLP_CYCLE - 1; i++) {
    //obj_[i]=obj_[i+1];
    in_[i] = in_[i + 1];
    out_[i] = out_[i + 1];
    way_[i] = way_[i + 1];
  }
  int way = 1 - wayIn + 4 * (1 - wayOut);
  //obj_[i]=model_->objectiveValue();
  in_[CLP_CYCLE - 1] = in;
  out_[CLP_CYCLE - 1] = out;
  way_[CLP_CYCLE - 1] = static_cast< char >(way);
  return matched;
}
#include "CoinStructuredModel.hpp"
// Solve using structure of model and maybe in parallel
int ClpSimplex::solve(CoinStructuredModel *model)
{
  // analyze structure
  int numberRowBlocks = model->numberRowBlocks();
  int numberColumnBlocks = model->numberColumnBlocks();
  int numberElementBlocks = model->numberElementBlocks();
  if (numberElementBlocks == 1) {
    loadProblem(*model, false);
    return dual();
  }
  // For now just get top level structure
  CoinModelBlockInfo *blockInfo = new CoinModelBlockInfo[numberElementBlocks];
  for (int i = 0; i < numberElementBlocks; i++) {
    CoinStructuredModel *subModel = dynamic_cast< CoinStructuredModel * >(model->block(i));
    CoinModel *thisBlock;
    if (subModel) {
      thisBlock = subModel->coinModelBlock(blockInfo[i]);
      model->setCoinModel(thisBlock, i);
    } else {
      thisBlock = dynamic_cast< CoinModel * >(model->block(i));
      assert(thisBlock);
      // just fill in info
      CoinModelBlockInfo info = CoinModelBlockInfo();
      int whatsSet = thisBlock->whatIsSet();
      info.matrix = static_cast< char >(((whatsSet & 1) != 0) ? 1 : 0);
      info.rhs = static_cast< char >(((whatsSet & 2) != 0) ? 1 : 0);
      info.rowName = static_cast< char >(((whatsSet & 4) != 0) ? 1 : 0);
      info.integer = static_cast< char >(((whatsSet & 32) != 0) ? 1 : 0);
      info.bounds = static_cast< char >(((whatsSet & 8) != 0) ? 1 : 0);
      info.columnName = static_cast< char >(((whatsSet & 16) != 0) ? 1 : 0);
      // Which block
      int iRowBlock = model->rowBlock(thisBlock->getRowBlock());
      info.rowBlock = iRowBlock;
      int iColumnBlock = model->columnBlock(thisBlock->getColumnBlock());
      info.columnBlock = iColumnBlock;
      blockInfo[i] = info;
    }
  }
  int *rowCounts = new int[numberRowBlocks];
  CoinZeroN(rowCounts, numberRowBlocks);
  int *columnCounts = new int[numberColumnBlocks + 1];
  CoinZeroN(columnCounts, numberColumnBlocks);
  int decomposeType = 0;
  for (int i = 0; i < numberElementBlocks; i++) {
    int iRowBlock = blockInfo[i].rowBlock;
    int iColumnBlock = blockInfo[i].columnBlock;
    rowCounts[iRowBlock]++;
    columnCounts[iColumnBlock]++;
  }
  if (numberRowBlocks == numberColumnBlocks || numberRowBlocks == numberColumnBlocks + 1) {
    // could be Dantzig-Wolfe
    int numberG1 = 0;
    for (int i = 0; i < numberRowBlocks; i++) {
      if (rowCounts[i] > 1)
        numberG1++;
    }
    bool masterColumns = (numberColumnBlocks == numberRowBlocks);
    if ((masterColumns && numberElementBlocks == 2 * numberRowBlocks - 1)
      || (!masterColumns && numberElementBlocks == 2 * numberRowBlocks)) {
      if (numberG1 < 2)
        decomposeType = 1;
    }
  }
  if (!decomposeType && (numberRowBlocks == numberColumnBlocks || numberRowBlocks == numberColumnBlocks - 1)) {
    // could be Benders
    int numberG1 = 0;
    for (int i = 0; i < numberColumnBlocks; i++) {
      if (columnCounts[i] > 1)
        numberG1++;
    }
    bool masterRows = (numberColumnBlocks == numberRowBlocks);
    if ((masterRows && numberElementBlocks == 2 * numberColumnBlocks - 1)
      || (!masterRows && numberElementBlocks == 2 * numberColumnBlocks)) {
      if (numberG1 < 2)
        decomposeType = 2;
    }
  }
  delete[] rowCounts;
  delete[] columnCounts;
  delete[] blockInfo;
  // decide what to do
  ClpSolve options;
  options.setIndependentOption(2, 100);
  switch (decomposeType) {
    // No good
  case 0:
    loadProblem(*model, false);
    return dual();
    // DW
  case 1:
    return solveDW(model, options);
    // Benders
  case 2:
    return solveBenders(model, options);
  }
  return 0; // to stop compiler warning
}
/* This loads a model from a CoinStructuredModel object - returns number of errors.
   If originalOrder then keep to order stored in blocks,
   otherwise first column/rows correspond to first block - etc.
   If keepSolution true and size is same as current then
   keeps current status and solution
*/
int ClpSimplex::loadProblem(CoinStructuredModel &coinModel,
  bool originalOrder,
  bool keepSolution)
{
  unsigned char *status = NULL;
  double *psol = NULL;
  double *dsol = NULL;
  int numberRows = coinModel.numberRows();
  int numberColumns = coinModel.numberColumns();
  int numberRowBlocks = coinModel.numberRowBlocks();
  int numberColumnBlocks = coinModel.numberColumnBlocks();
  int numberElementBlocks = coinModel.numberElementBlocks();
  if (status_ && numberRows_ && numberRows_ == numberRows && numberColumns_ == numberColumns && keepSolution) {
    status = new unsigned char[numberRows_ + numberColumns_];
    CoinMemcpyN(status_, numberRows_ + numberColumns_, status);
    psol = new double[numberRows_ + numberColumns_];
    CoinMemcpyN(columnActivity_, numberColumns_, psol);
    CoinMemcpyN(rowActivity_, numberRows_, psol + numberColumns_);
    dsol = new double[numberRows_ + numberColumns_];
    CoinMemcpyN(reducedCost_, numberColumns_, dsol);
    CoinMemcpyN(dual_, numberRows_, dsol + numberColumns_);
  }
  int returnCode = 0;
  double *rowLower = new double[numberRows];
  double *rowUpper = new double[numberRows];
  double *columnLower = new double[numberColumns];
  double *columnUpper = new double[numberColumns];
  double *objective = new double[numberColumns];
  int *integerType = new int[numberColumns];
  CoinBigIndex numberElements = 0;
  // Bases for blocks
  int *rowBase = new int[numberRowBlocks];
  CoinFillN(rowBase, numberRowBlocks, -1);
  // And row to put it
  int *whichRow = new int[numberRows + numberRowBlocks];
  int *columnBase = new int[numberColumnBlocks];
  CoinFillN(columnBase, numberColumnBlocks, -1);
  // And column to put it
  int *whichColumn = new int[numberColumns + numberColumnBlocks];
  for (int iBlock = 0; iBlock < numberElementBlocks; iBlock++) {
    CoinModel *block = coinModel.coinBlock(iBlock);
    numberElements += block->numberElements();
    //and set up elements etc
    double *associated = block->associatedArray();
    // If strings then do copies
    if (block->stringsExist())
      returnCode += block->createArrays(rowLower, rowUpper, columnLower, columnUpper,
        objective, integerType, associated);
    const CoinModelBlockInfo &info = coinModel.blockType(iBlock);
    int iRowBlock = info.rowBlock;
    int iColumnBlock = info.columnBlock;
    if (rowBase[iRowBlock] < 0) {
      rowBase[iRowBlock] = block->numberRows();
      // Save block number
      whichRow[numberRows + iRowBlock] = iBlock;
    } else {
      assert(rowBase[iRowBlock] == block->numberRows());
    }
    if (columnBase[iColumnBlock] < 0) {
      columnBase[iColumnBlock] = block->numberColumns();
      // Save block number
      whichColumn[numberColumns + iColumnBlock] = iBlock;
    } else {
      assert(columnBase[iColumnBlock] == block->numberColumns());
    }
  }
  // Fill arrays with defaults
  CoinFillN(rowLower, numberRows, -COIN_DBL_MAX);
  CoinFillN(rowUpper, numberRows, COIN_DBL_MAX);
  CoinFillN(columnLower, numberColumns, 0.0);
  CoinFillN(columnUpper, numberColumns, COIN_DBL_MAX);
  CoinFillN(objective, numberColumns, 0.0);
  CoinFillN(integerType, numberColumns, 0);
  int n = 0;
  for (int iBlock = 0; iBlock < numberRowBlocks; iBlock++) {
    int k = rowBase[iBlock];
    rowBase[iBlock] = n;
    assert(k >= 0);
    // block number
    int jBlock = whichRow[numberRows + iBlock];
    if (originalOrder) {
      memcpy(whichRow + n, coinModel.coinBlock(jBlock)->originalRows(), k * sizeof(int));
    } else {
      CoinIotaN(whichRow + n, k, n);
    }
    n += k;
  }
  assert(n == numberRows);
  n = 0;
  for (int iBlock = 0; iBlock < numberColumnBlocks; iBlock++) {
    int k = columnBase[iBlock];
    columnBase[iBlock] = n;
    assert(k >= 0);
    if (k) {
      // block number
      int jBlock = whichColumn[numberColumns + iBlock];
      if (originalOrder) {
        memcpy(whichColumn + n, coinModel.coinBlock(jBlock)->originalColumns(),
          k * sizeof(int));
      } else {
        CoinIotaN(whichColumn + n, k, n);
      }
      n += k;
    }
  }
  assert(n == numberColumns);
  bool gotIntegers = false;
  for (int iBlock = 0; iBlock < numberElementBlocks; iBlock++) {
    CoinModel *block = coinModel.coinBlock(iBlock);
    const CoinModelBlockInfo &info = coinModel.blockType(iBlock);
    int iRowBlock = info.rowBlock;
    int iRowBase = rowBase[iRowBlock];
    int iColumnBlock = info.columnBlock;
    int iColumnBase = columnBase[iColumnBlock];
    if (info.rhs) {
      int nRows = block->numberRows();
      const double *lower = block->rowLowerArray();
      const double *upper = block->rowUpperArray();
      for (int i = 0; i < nRows; i++) {
        int put = whichRow[i + iRowBase];
        rowLower[put] = lower[i];
        rowUpper[put] = upper[i];
      }
    }
    if (info.bounds) {
      int nColumns = block->numberColumns();
      const double *lower = block->columnLowerArray();
      const double *upper = block->columnUpperArray();
      const double *obj = block->objectiveArray();
      for (int i = 0; i < nColumns; i++) {
        int put = whichColumn[i + iColumnBase];
        columnLower[put] = lower[i];
        columnUpper[put] = upper[i];
        objective[put] = obj[i];
      }
    }
    if (info.integer) {
      gotIntegers = true;
      int nColumns = block->numberColumns();
      const int *type = block->integerTypeArray();
      for (int i = 0; i < nColumns; i++) {
        int put = whichColumn[i + iColumnBase];
        integerType[put] = type[i];
      }
    }
  }
  gutsOfLoadModel(numberRows, numberColumns,
    columnLower, columnUpper, objective, rowLower, rowUpper, NULL);
  delete[] rowLower;
  delete[] rowUpper;
  delete[] columnLower;
  delete[] columnUpper;
  delete[] objective;
  // Do integers if wanted
  if (gotIntegers) {
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (integerType[iColumn])
        setInteger(iColumn);
    }
  }
  delete[] integerType;
  setObjectiveOffset(coinModel.objectiveOffset());
  // Space for elements
  int *row = new int[numberElements];
  int *column = new int[numberElements];
  double *element = new double[numberElements];
  numberElements = 0;
  for (int iBlock = 0; iBlock < numberElementBlocks; iBlock++) {
    CoinModel *block = coinModel.coinBlock(iBlock);
    const CoinModelBlockInfo &info = coinModel.blockType(iBlock);
    int iRowBlock = info.rowBlock;
    int iRowBase = rowBase[iRowBlock];
    int iColumnBlock = info.columnBlock;
    int iColumnBase = columnBase[iColumnBlock];
    if (info.rowName) {
      int numberItems = block->rowNames()->numberItems();
      assert(block->numberRows() >= numberItems);
      if (numberItems) {
        const char *const *rowNames = block->rowNames()->names();
        for (int i = 0; i < numberItems; i++) {
          int put = whichRow[i + iRowBase];
          std::string name = rowNames[i];
          setRowName(put, name);
        }
      }
    }
    if (info.columnName) {
      int numberItems = block->columnNames()->numberItems();
      assert(block->numberColumns() >= numberItems);
      if (numberItems) {
        const char *const *columnNames = block->columnNames()->names();
        for (int i = 0; i < numberItems; i++) {
          int put = whichColumn[i + iColumnBase];
          std::string name = columnNames[i];
          setColumnName(put, name);
        }
      }
    }
    if (info.matrix) {
      CoinPackedMatrix matrix2;
      const CoinPackedMatrix *matrix = block->packedMatrix();
      if (!matrix) {
        double *associated = block->associatedArray();
        block->createPackedMatrix(matrix2, associated);
        matrix = &matrix2;
      }
      // get matrix data pointers
      const int *row2 = matrix->getIndices();
      const CoinBigIndex *columnStart = matrix->getVectorStarts();
      const double *elementByColumn = matrix->getElements();
      const int *columnLength = matrix->getVectorLengths();
      int n = matrix->getNumCols();
      assert(matrix->isColOrdered());
      for (int iColumn = 0; iColumn < n; iColumn++) {
        CoinBigIndex j;
        int jColumn = whichColumn[iColumn + iColumnBase];
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          row[numberElements] = whichRow[row2[j] + iRowBase];
          column[numberElements] = jColumn;
          element[numberElements++] = elementByColumn[j];
        }
      }
    }
  }
  delete[] whichRow;
  delete[] whichColumn;
  delete[] rowBase;
  delete[] columnBase;
  CoinPackedMatrix *matrix = new CoinPackedMatrix(true, row, column, element, numberElements);
  matrix_ = new ClpPackedMatrix(matrix);
  matrix_->setDimensions(numberRows, numberColumns);
  delete[] row;
  delete[] column;
  delete[] element;
  createStatus();
  if (status) {
    // copy back
    CoinMemcpyN(status, numberRows_ + numberColumns_, status_);
    CoinMemcpyN(psol, numberColumns_, columnActivity_);
    CoinMemcpyN(psol + numberColumns_, numberRows_, rowActivity_);
    CoinMemcpyN(dsol, numberColumns_, reducedCost_);
    CoinMemcpyN(dsol + numberColumns_, numberRows_, dual_);
    delete[] status;
    delete[] psol;
    delete[] dsol;
  }
  optimizationDirection_ = coinModel.optimizationDirection();
  return returnCode;
}
/*  If input negative scales objective so maximum <= -value
    and returns scale factor used.  If positive unscales and also
    redoes dual stuff
*/
double
ClpSimplex::scaleObjective(double value)
{
  double *obj = objective();
  double largest = 0.0;
  if (value < 0.0) {
    value = -value;
    for (int i = 0; i < numberColumns_; i++) {
      largest = std::max(largest, fabs(obj[i]));
    }
    if (largest > value) {
      double scaleFactor = value / largest;
      for (int i = 0; i < numberColumns_; i++) {
        obj[i] *= scaleFactor;
        reducedCost_[i] *= scaleFactor;
      }
      for (int i = 0; i < numberRows_; i++) {
        dual_[i] *= scaleFactor;
      }
      largest /= value;
    } else {
      // no need
      largest = 1.0;
    }
  } else {
    // at end
    if (value != 1.0) {
      for (int i = 0; i < numberColumns_; i++) {
        obj[i] *= value;
        reducedCost_[i] *= value;
      }
      for (int i = 0; i < numberRows_; i++) {
        dual_[i] *= value;
      }
      computeObjectiveValue();
    }
  }
  return largest;
}
// Solve using Dantzig-Wolfe decomposition and maybe in parallel
int ClpSimplex::solveDW(CoinStructuredModel *model, ClpSolve &options)
{
  double time1 = CoinCpuTime();
  int numberColumns = model->numberColumns();
  int numberRowBlocks = model->numberRowBlocks();
  int numberColumnBlocks = model->numberColumnBlocks();
  int numberElementBlocks = model->numberElementBlocks();
  // We already have top level structure
  CoinModelBlockInfo *blockInfo = new CoinModelBlockInfo[numberElementBlocks];
  for (int i = 0; i < numberElementBlocks; i++) {
    CoinModel *thisBlock = model->coinBlock(i);
    assert(thisBlock);
    // just fill in info
    CoinModelBlockInfo info = CoinModelBlockInfo();
    int whatsSet = thisBlock->whatIsSet();
    info.matrix = static_cast< char >(((whatsSet & 1) != 0) ? 1 : 0);
    info.rhs = static_cast< char >(((whatsSet & 2) != 0) ? 1 : 0);
    info.rowName = static_cast< char >(((whatsSet & 4) != 0) ? 1 : 0);
    info.integer = static_cast< char >(((whatsSet & 32) != 0) ? 1 : 0);
    info.bounds = static_cast< char >(((whatsSet & 8) != 0) ? 1 : 0);
    info.columnName = static_cast< char >(((whatsSet & 16) != 0) ? 1 : 0);
    // Which block
    int iRowBlock = model->rowBlock(thisBlock->getRowBlock());
    info.rowBlock = iRowBlock;
    int iColumnBlock = model->columnBlock(thisBlock->getColumnBlock());
    info.columnBlock = iColumnBlock;
    blockInfo[i] = info;
  }
  // make up problems
  int numberBlocks = numberRowBlocks - 1;
  // Find master rows and columns
  int *rowCounts = new int[numberRowBlocks];
  CoinZeroN(rowCounts, numberRowBlocks);
  int *columnCounts = new int[numberColumnBlocks + 1];
  CoinZeroN(columnCounts, numberColumnBlocks);
  int iBlock;
  for (iBlock = 0; iBlock < numberElementBlocks; iBlock++) {
    int iRowBlock = blockInfo[iBlock].rowBlock;
    rowCounts[iRowBlock]++;
    int iColumnBlock = blockInfo[iBlock].columnBlock;
    columnCounts[iColumnBlock]++;
  }
  int *whichBlock = new int[numberElementBlocks];
  int masterRowBlock = -1;
  for (iBlock = 0; iBlock < numberRowBlocks; iBlock++) {
    if (rowCounts[iBlock] > 1) {
      if (masterRowBlock == -1) {
        masterRowBlock = iBlock;
      } else {
        // Can't decode
        masterRowBlock = -2;
        break;
      }
    }
  }
  int masterColumnBlock = -1;
  int kBlock = 0;
  for (iBlock = 0; iBlock < numberColumnBlocks; iBlock++) {
    int count = columnCounts[iBlock];
    columnCounts[iBlock] = kBlock;
    kBlock += count;
  }
  for (iBlock = 0; iBlock < numberElementBlocks; iBlock++) {
    int iColumnBlock = blockInfo[iBlock].columnBlock;
    whichBlock[columnCounts[iColumnBlock]] = iBlock;
    columnCounts[iColumnBlock]++;
  }
  for (iBlock = numberColumnBlocks - 1; iBlock >= 0; iBlock--)
    columnCounts[iBlock + 1] = columnCounts[iBlock];
  columnCounts[0] = 0;
  for (iBlock = 0; iBlock < numberColumnBlocks; iBlock++) {
    int count = columnCounts[iBlock + 1] - columnCounts[iBlock];
    if (count == 1) {
      int kBlock = whichBlock[columnCounts[iBlock]];
      int iRowBlock = blockInfo[kBlock].rowBlock;
      if (iRowBlock == masterRowBlock) {
        if (masterColumnBlock == -1) {
          masterColumnBlock = iBlock;
        } else {
          // Can't decode
          masterColumnBlock = -2;
          break;
        }
      }
    }
  }
  if (masterRowBlock < 0 || masterColumnBlock == -2) {
    // What now
    abort();
  }
  delete[] rowCounts;
  // create all data
  const CoinPackedMatrix **top = new const CoinPackedMatrix *[numberColumnBlocks];
  ClpSimplex *sub = new ClpSimplex[numberBlocks];
  ClpSimplex master;
  // Set offset
  master.setObjectiveOffset(model->objectiveOffset());
  bool reducePrint = logLevel() == 7;
  if (reducePrint) {
    // special
    setLogLevel(1);
    master.setLogLevel(1);
  }
  kBlock = 0;
  int masterBlock = -1;
  for (iBlock = 0; iBlock < numberColumnBlocks; iBlock++) {
    top[kBlock] = NULL;
    int start = columnCounts[iBlock];
    int end = columnCounts[iBlock + 1];
    assert(end - start <= 2);
    for (int j = start; j < end; j++) {
      int jBlock = whichBlock[j];
      int iRowBlock = blockInfo[jBlock].rowBlock;
      int iColumnBlock = blockInfo[jBlock].columnBlock;
      assert(iColumnBlock == iBlock);
      if (iColumnBlock != masterColumnBlock && iRowBlock == masterRowBlock) {
        // top matrix
        top[kBlock] = model->coinBlock(jBlock)->packedMatrix();
      } else {
        const CoinPackedMatrix *matrix
          = model->coinBlock(jBlock)->packedMatrix();
        // Get pointers to arrays
        const double *rowLower;
        const double *rowUpper;
        const double *columnLower;
        const double *columnUpper;
        const double *objective;
        model->block(iRowBlock, iColumnBlock, rowLower, rowUpper,
          columnLower, columnUpper, objective);
        if (iColumnBlock != masterColumnBlock) {
          // diagonal block
          sub[kBlock].loadProblem(*matrix, columnLower, columnUpper,
            objective, rowLower, rowUpper);
          if (true) {
            double *lower = sub[kBlock].columnLower();
            double *upper = sub[kBlock].columnUpper();
            int n = sub[kBlock].numberColumns();
            for (int i = 0; i < n; i++) {
              lower[i] = std::max(-1.0e8, lower[i]);
              upper[i] = std::min(1.0e8, upper[i]);
            }
          }
          if (optimizationDirection_ < 0.0) {
            double *obj = sub[kBlock].objective();
            int n = sub[kBlock].numberColumns();
            for (int i = 0; i < n; i++)
              obj[i] = -obj[i];
          }
          if (this->factorizationFrequency() == 200) {
            // User did not touch preset
            sub[kBlock].defaultFactorizationFrequency();
          } else {
            // make sure model has correct value
            sub[kBlock].setFactorizationFrequency(this->factorizationFrequency());
          }
          sub[kBlock].setPerturbation(50);
          // Set columnCounts to be diagonal block index for cleanup
          columnCounts[kBlock] = jBlock;
        } else {
          // master
          masterBlock = jBlock;
          master.loadProblem(*matrix, columnLower, columnUpper,
            objective, rowLower, rowUpper);
          if (optimizationDirection_ < 0.0) {
            double *obj = master.objective();
            int n = master.numberColumns();
            for (int i = 0; i < n; i++)
              obj[i] = -obj[i];
          }
        }
      }
    }
    if (iBlock != masterColumnBlock)
      kBlock++;
  }
  delete[] whichBlock;
  delete[] blockInfo;
  // For now master must have been defined (does not have to have columns)
  assert(master.numberRows());
  assert(masterBlock >= 0);
  int numberMasterRows = master.numberRows();
  // Overkill in terms of space
  int spaceNeeded = std::max(numberBlocks * (numberMasterRows + 1),
    2 * numberMasterRows);
  int *rowAdd = new int[spaceNeeded];
  double *elementAdd = new double[spaceNeeded];
  spaceNeeded = numberBlocks;
  CoinBigIndex *columnAdd = new CoinBigIndex[spaceNeeded + 1];
  double *objective = new double[spaceNeeded];
  // Add in costed slacks
  int firstArtificial = master.numberColumns();
  int lastArtificial = firstArtificial;
  if (true) {
    const double *lower = master.rowLower();
    const double *upper = master.rowUpper();
    int kCol = 0;
    for (int iRow = 0; iRow < numberMasterRows; iRow++) {
      if (lower[iRow] > -1.0e10) {
        rowAdd[kCol] = iRow;
        elementAdd[kCol++] = 1.0;
      }
      if (upper[iRow] < 1.0e10) {
        rowAdd[kCol] = iRow;
        elementAdd[kCol++] = -1.0;
      }
    }
    if (kCol > spaceNeeded) {
      spaceNeeded = kCol;
      delete[] columnAdd;
      delete[] objective;
      columnAdd = new CoinBigIndex[spaceNeeded + 1];
      objective = new double[spaceNeeded];
    }
    for (int i = 0; i < kCol; i++) {
      columnAdd[i] = i;
      objective[i] = 1.0e13;
    }
    columnAdd[kCol] = kCol;
    master.addColumns(kCol, NULL, NULL, objective,
      columnAdd, rowAdd, elementAdd);
    lastArtificial = master.numberColumns();
  }
  int maxPass = options.independentOption(2);
  if (maxPass < 2)
    maxPass = 100;
  int iPass;
  double lastObjective = 1.0e31;
  // Create convexity rows for proposals
  int numberMasterColumns = master.numberColumns();
  master.resize(numberMasterRows + numberBlocks, numberMasterColumns);
  if (this->factorizationFrequency() == 200) {
    // User did not touch preset
    master.defaultFactorizationFrequency();
  } else {
    // make sure model has correct value
    master.setFactorizationFrequency(this->factorizationFrequency());
  }
  master.setPerturbation(50);
  // Arrays to say which block and when created
  int maximumColumns = 2 * numberMasterRows + 10 * numberBlocks;
  whichBlock = new int[maximumColumns];
  int *when = new int[maximumColumns];
  int numberColumnsGenerated = numberBlocks;
  // fill in rhs and add in artificials
  {
    double *rowLower = master.rowLower();
    double *rowUpper = master.rowUpper();
    int iBlock;
    columnAdd[0] = 0;
    for (iBlock = 0; iBlock < numberBlocks; iBlock++) {
      int iRow = iBlock + numberMasterRows;
      ;
      rowLower[iRow] = 1.0;
      rowUpper[iRow] = 1.0;
      rowAdd[iBlock] = iRow;
      elementAdd[iBlock] = 1.0;
      objective[iBlock] = 1.0e13;
      columnAdd[iBlock + 1] = iBlock + 1;
      when[iBlock] = -1;
      whichBlock[iBlock] = iBlock;
    }
    master.addColumns(numberBlocks, NULL, NULL, objective,
      columnAdd, rowAdd, elementAdd);
  }
  char generalPrint[200];
  // and resize matrix to double check clp will be happy
  //master.matrix()->setDimensions(numberMasterRows+numberBlocks,
  //			 numberMasterColumns+numberBlocks);
  sprintf(generalPrint, "Time to decompose %.2f seconds", CoinCpuTime() - time1);
  handler_->message(CLP_GENERAL, messages_)
    << generalPrint
    << CoinMessageEol;
  for (iPass = 0; iPass < maxPass; iPass++) {
    sprintf(generalPrint, "Start of pass %d", iPass);
    handler_->message(CLP_GENERAL, messages_)
      << generalPrint
      << CoinMessageEol;
    // Solve master - may be infeasible
    //master.scaling(0);
    if (0) {
      master.writeMps("yy.mps");
    }
    // Correct artificials
    double sumArtificials = 0.0;
    if (iPass) {
      double *upper = master.columnUpper();
      double *solution = master.primalColumnSolution();
      double *obj = master.objective();
      sumArtificials = 0.0;
      for (int i = firstArtificial; i < lastArtificial; i++) {
        sumArtificials += solution[i];
        //assert (solution[i]>-1.0e-2);
        if (solution[i] < 1.0e-6) {
          obj[i] = 1.0e7;
          upper[i] = 1.0e-1;
          solution[i] = 0.0;
          master.setColumnStatus(i, isFixed);
        } else {
          upper[i] = solution[i] + 1.0e-5 * (1.0 + solution[i]);
        }
      }
      sprintf(generalPrint, "Sum of artificials before solve is %g", sumArtificials);
      handler_->message(CLP_GENERAL, messages_)
        << generalPrint
        << CoinMessageEol;
    }
    // scale objective to be reasonable
    double scaleFactor = master.scaleObjective(-1.0e9);
    {
      double *dual = master.dualRowSolution();
      int n = master.numberRows();
      memset(dual, 0, n * sizeof(double));
      double *solution = master.primalColumnSolution();
      master.clpMatrix()->times(1.0, solution, dual);
      double sum = 0.0;
      double *lower = master.rowLower();
      double *upper = master.rowUpper();
      for (int iRow = 0; iRow < n; iRow++) {
        double value = dual[iRow];
        if (value > upper[iRow])
          sum += value - upper[iRow];
        else if (value < lower[iRow])
          sum -= value - lower[iRow];
      }
      printf("** suminf %g\n", sum);
      lower = master.columnLower();
      upper = master.columnUpper();
      n = master.numberColumns();
      for (int iColumn = 0; iColumn < n; iColumn++) {
        double value = solution[iColumn];
        if (value > upper[iColumn] + 1.0e-5)
          sum += value - upper[iColumn];
        else if (value < lower[iColumn] - 1.0e-5)
          sum -= value - lower[iColumn];
      }
      printf("** suminf %g\n", sum);
    }
    master.primal();
    //master.primal(1);
    // Correct artificials
    sumArtificials = 0.0;
    {
      double *solution = master.primalColumnSolution();
      for (int i = firstArtificial; i < lastArtificial; i++) {
        sumArtificials += solution[i];
      }
      printf("** Sum of artificials after solve is %g\n", sumArtificials);
    }
    master.scaleObjective(scaleFactor);
    int problemStatus = master.status(); // do here as can change (delcols)
    if (problemStatus == 2 && master.numberColumns()) {
      master.primal(1);
      //master.primal(1);
      if (problemStatus == 2) {
        int numberColumns = master.numberColumns();
        double *lower = master.columnLower();
        double *upper = master.columnUpper();
        for (int i = 0; i < numberColumns; i++) {
          lower[i] = std::max(lower[i], -1.0e10);
          upper[i] = std::min(upper[i], 1.0e10);
        }
        master.primal(1);
        //master.primal(1);
        assert(problemStatus != 2);
      }
    }
    if (master.numberIterations() == 0 && iPass)
      break; // finished
    if (master.objectiveValue() > lastObjective - 1.0e-7 && iPass > 555)
      break; // finished
    lastObjective = master.objectiveValue();
    // mark basic ones and delete if necessary
    int iColumn;
    numberColumnsGenerated = master.numberColumns() - numberMasterColumns;
    for (iColumn = 0; iColumn < numberColumnsGenerated; iColumn++) {
      if (master.getStatus(iColumn + numberMasterColumns) == ClpSimplex::basic)
        when[iColumn] = iPass;
    }
    if (numberColumnsGenerated + numberBlocks > maximumColumns) {
      // delete
      int numberKeep = 0;
      int numberDelete = 0;
      int *whichDelete = new int[numberColumnsGenerated];
      for (iColumn = 0; iColumn < numberColumnsGenerated; iColumn++) {
        if (when[iColumn] > iPass - 7) {
          // keep
          when[numberKeep] = when[iColumn];
          whichBlock[numberKeep++] = whichBlock[iColumn];
        } else {
          // delete
          whichDelete[numberDelete++] = iColumn + numberMasterColumns;
        }
      }
      numberColumnsGenerated -= numberDelete;
      master.deleteColumns(numberDelete, whichDelete);
      delete[] whichDelete;
    }
    const double *dual = NULL;
    bool deleteDual = false;
    if (problemStatus == 0) {
      dual = master.dualRowSolution();
    } else if (problemStatus == 1) {
      // could do composite objective
      dual = master.infeasibilityRay();
      deleteDual = true;
      printf("** The sum of infeasibilities is %g\n",
        master.sumPrimalInfeasibilities());
    } else if (!master.numberColumns()) {
      assert(!iPass);
      dual = master.dualRowSolution();
      memset(master.dualRowSolution(),
        0, (numberMasterRows + numberBlocks) * sizeof(double));
    } else {
      master.writeMps("unbounded.mps");
      abort();
    }
    // Scale back on first time
    if (!iPass) {
      double *dual2 = master.dualRowSolution();
      for (int i = 0; i < numberMasterRows + numberBlocks; i++) {
        dual2[i] *= 1.0e-7;
      }
      dual = master.dualRowSolution();
    }
    // Create objective for sub problems and solve
    columnAdd[0] = 0;
    int numberProposals = 0;
    double **saveObj2 = new double *[numberBlocks];
    for (iBlock = 0; iBlock < numberBlocks; iBlock++) {
      int numberColumns2 = sub[iBlock].numberColumns();
      saveObj2[iBlock] = new double[numberColumns2];
    }
    for (iBlock = 0; iBlock < numberBlocks; iBlock++) {
      int numberColumns2 = sub[iBlock].numberColumns();
      double *saveObj = saveObj2[iBlock];
      double *objective2 = sub[iBlock].objective();
      memcpy(saveObj, objective2, numberColumns2 * sizeof(double));
      // new objective
      top[iBlock]->transposeTimes(dual, objective2);
      int i;
      if (problemStatus == 0) {
        for (i = 0; i < numberColumns2; i++)
          objective2[i] = saveObj[i] - objective2[i];
      } else {
        for (i = 0; i < numberColumns2; i++)
          objective2[i] = -objective2[i];
      }
      // scale objective to be reasonable
      //double scaleFactor =
      //   sub[iBlock].scaleObjective((sumArtificials > 1.0e-5) ? -1.0e-4 : -1.0e9);

      if (reducePrint)
        sub[iBlock].setLogLevel(0);
    }
      for (iBlock = 0; iBlock < numberBlocks; iBlock++) {
        if (iPass) {
          sub[iBlock].primal();
        } else {
          sub[iBlock].dual();
        }
      }
    for (iBlock = 0; iBlock < numberBlocks; iBlock++) {
      int numberColumns2 = sub[iBlock].numberColumns();
      double *saveObj = saveObj2[iBlock];
      double *objective2 = sub[iBlock].objective();
      int i;
      sub[iBlock].scaleObjective(scaleFactor);
      if (!sub[iBlock].isProvenOptimal() && !sub[iBlock].isProvenDualInfeasible()) {
        memset(objective2, 0, numberColumns2 * sizeof(double));
        sub[iBlock].primal();
        if (problemStatus == 0) {
          for (int i = 0; i < numberColumns2; i++)
            objective2[i] = saveObj[i] - objective2[i];
        } else {
          for (i = 0; i < numberColumns2; i++)
            objective2[i] = -objective2[i];
        }
        double scaleFactor = sub[iBlock].scaleObjective(-1.0e9);
        sub[iBlock].primal(1);
        sub[iBlock].scaleObjective(scaleFactor);
      }
      memcpy(objective2, saveObj, numberColumns2 * sizeof(double));
      // get proposal
      if (sub[iBlock].numberIterations() || !iPass) {
        double objValue = 0.0;
        CoinBigIndex start = columnAdd[numberProposals];
        // proposal
        if (sub[iBlock].isProvenOptimal()) {
          const double *solution = sub[iBlock].primalColumnSolution();
          top[iBlock]->times(solution, elementAdd + start);
          for (i = 0; i < numberColumns2; i++)
            objValue += solution[i] * saveObj[i];
          // See if good dj and pack down
          CoinBigIndex number = start;
          double dj = objValue;
          if (problemStatus)
            dj = 0.0;
          double smallest = 1.0e100;
          double largest = 0.0;
          for (i = 0; i < numberMasterRows; i++) {
            double value = elementAdd[start + i];
            if (fabs(value) > 1.0e-15) {
              dj -= dual[i] * value;
              smallest = std::min(smallest, fabs(value));
              largest = std::max(largest, fabs(value));
              rowAdd[number] = i;
              elementAdd[number++] = value;
            }
          }
          // and convexity
          dj -= dual[numberMasterRows + iBlock];
          rowAdd[number] = numberMasterRows + iBlock;
          elementAdd[number++] = 1.0;
          // if elements large then scale?
          //if (largest>1.0e8||smallest<1.0e-8)
          sprintf(generalPrint, "For subproblem %d smallest - %g, largest %g - dj %g",
            iBlock, smallest, largest, dj);
          handler_->message(CLP_GENERAL2, messages_)
            << generalPrint
            << CoinMessageEol;
          if (dj < -1.0e-6 || !iPass) {
            // take
            objective[numberProposals] = objValue;
            columnAdd[++numberProposals] = number;
            when[numberColumnsGenerated] = iPass;
            whichBlock[numberColumnsGenerated++] = iBlock;
          }
        } else if (sub[iBlock].isProvenDualInfeasible()) {
          // use ray
          const double *solution = sub[iBlock].unboundedRay();
          top[iBlock]->times(solution, elementAdd + start);
          for (i = 0; i < numberColumns2; i++)
            objValue += solution[i] * saveObj[i];
          // See if good dj and pack down
          CoinBigIndex number = start;
          double dj = objValue;
          double smallest = 1.0e100;
          double largest = 0.0;
          for (i = 0; i < numberMasterRows; i++) {
            double value = elementAdd[start + i];
            if (fabs(value) > 1.0e-15) {
              dj -= dual[i] * value;
              smallest = std::min(smallest, fabs(value));
              largest = std::max(largest, fabs(value));
              rowAdd[number] = i;
              elementAdd[number++] = value;
            }
          }
          // if elements large or small then scale?
          //if (largest>1.0e8||smallest<1.0e-8)
          sprintf(generalPrint, "For subproblem ray %d smallest - %g, largest %g - dj %g",
            iBlock, smallest, largest, dj);
          handler_->message(CLP_GENERAL2, messages_)
            << generalPrint
            << CoinMessageEol;
          if (dj < -1.0e-6) {
            // take
            objective[numberProposals] = objValue;
            columnAdd[++numberProposals] = number;
            when[numberColumnsGenerated] = iPass;
            whichBlock[numberColumnsGenerated++] = iBlock;
          }
        } else {
          abort();
        }
      }
    }
    for (iBlock = 0; iBlock < numberBlocks; iBlock++)
      delete[] saveObj2[iBlock];
    delete[] saveObj2;
    if (deleteDual)
      delete[] dual;
    if (numberProposals)
      master.addColumns(numberProposals, NULL, NULL, objective,
        columnAdd, rowAdd, elementAdd);
  }
  sprintf(generalPrint, "Time at end of D-W %.2f seconds", CoinCpuTime() - time1);
  handler_->message(CLP_GENERAL, messages_)
    << generalPrint
    << CoinMessageEol;
  //master.scaling(0);
  //master.primal(1);
  loadProblem(*model);
  // now put back a good solution
  double *lower = new double[numberMasterRows + numberBlocks];
  double *upper = new double[numberMasterRows + numberBlocks];
  numberColumnsGenerated += numberMasterColumns;
  double *sol = new double[numberColumnsGenerated];
  const double *solution = master.primalColumnSolution();
  const double *masterLower = master.rowLower();
  const double *masterUpper = master.rowUpper();
  double *fullSolution = primalColumnSolution();
  const double *fullLower = columnLower();
  const double *fullUpper = columnUpper();
  const double *rowSolution = master.primalRowSolution();
  double *fullRowSolution = primalRowSolution();
  const int *rowBack = model->coinBlock(masterBlock)->originalRows();
  int numberRows2 = model->coinBlock(masterBlock)->numberRows();
  const int *columnBack = model->coinBlock(masterBlock)->originalColumns();
  int numberColumns2 = model->coinBlock(masterBlock)->numberColumns();
  for (int iRow = 0; iRow < numberRows2; iRow++) {
    int kRow = rowBack[iRow];
    setRowStatus(kRow, master.getRowStatus(iRow));
    fullRowSolution[kRow] = rowSolution[iRow];
  }
  for (int iColumn = 0; iColumn < numberColumns2; iColumn++) {
    int kColumn = columnBack[iColumn];
    setStatus(kColumn, master.getStatus(iColumn));
    fullSolution[kColumn] = solution[iColumn];
  }
  for (iBlock = 0; iBlock < numberBlocks; iBlock++) {
    // move basis
    int kBlock = columnCounts[iBlock];
    const int *rowBack = model->coinBlock(kBlock)->originalRows();
    int numberRows2 = model->coinBlock(kBlock)->numberRows();
    const int *columnBack = model->coinBlock(kBlock)->originalColumns();
    int numberColumns2 = model->coinBlock(kBlock)->numberColumns();
    for (int iRow = 0; iRow < numberRows2; iRow++) {
      int kRow = rowBack[iRow];
      setRowStatus(kRow, sub[iBlock].getRowStatus(iRow));
    }
    for (int iColumn = 0; iColumn < numberColumns2; iColumn++) {
      int kColumn = columnBack[iColumn];
      setStatus(kColumn, sub[iBlock].getStatus(iColumn));
    }
    // convert top bit to by rows
    CoinPackedMatrix topMatrix = *top[iBlock];
    topMatrix.reverseOrdering();
    // zero solution
    memset(sol, 0, numberColumnsGenerated * sizeof(double));

    for (int i = numberMasterColumns; i < numberColumnsGenerated; i++) {
      if (whichBlock[i - numberMasterColumns] == iBlock)
        sol[i] = solution[i];
    }
    memset(lower, 0, (numberMasterRows + numberBlocks) * sizeof(double));
    master.clpMatrix()->times(1.0, sol, lower);
    for (int iRow = 0; iRow < numberMasterRows; iRow++) {
      double value = lower[iRow];
      if (masterUpper[iRow] < 1.0e20)
        upper[iRow] = value;
      else
        upper[iRow] = COIN_DBL_MAX;
      if (masterLower[iRow] > -1.0e20)
        lower[iRow] = value;
      else
        lower[iRow] = -COIN_DBL_MAX;
    }
    sub[iBlock].addRows(numberMasterRows, lower, upper,
      topMatrix.getVectorStarts(),
      topMatrix.getVectorLengths(),
      topMatrix.getIndices(),
      topMatrix.getElements());
    sub[iBlock].primal(1);
    const double *subSolution = sub[iBlock].primalColumnSolution();
    const double *subRowSolution = sub[iBlock].primalRowSolution();
    // move solution
    for (int iRow = 0; iRow < numberRows2; iRow++) {
      int kRow = rowBack[iRow];
      fullRowSolution[kRow] = subRowSolution[iRow];
    }
    for (int iColumn = 0; iColumn < numberColumns2; iColumn++) {
      int kColumn = columnBack[iColumn];
      fullSolution[kColumn] = subSolution[iColumn];
    }
  }
  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (fullSolution[iColumn] < fullUpper[iColumn] - 1.0e-8 && fullSolution[iColumn] > fullLower[iColumn] + 1.0e-8) {
      if (getStatus(iColumn) != ClpSimplex::basic) {
        if (columnLower_[iColumn] > -1.0e30 || columnUpper_[iColumn] < 1.0e30)
          setStatus(iColumn, ClpSimplex::superBasic);
        else
          setStatus(iColumn, ClpSimplex::isFree);
      }
    } else if (fullSolution[iColumn] >= fullUpper[iColumn] - 1.0e-8) {
      // may help to make rest non basic
      if (getStatus(iColumn) != ClpSimplex::basic)
        setStatus(iColumn, ClpSimplex::atUpperBound);
    } else if (fullSolution[iColumn] <= fullLower[iColumn] + 1.0e-8) {
      // may help to make rest non basic
      if (getStatus(iColumn) != ClpSimplex::basic)
        setStatus(iColumn, ClpSimplex::atLowerBound);
    }
  }
  //int numberRows=model->numberRows();
  //for (int iRow=0;iRow<numberRows;iRow++)
  //setRowStatus(iRow,ClpSimplex::superBasic);
  sprintf(generalPrint, "Time before cleanup of full model %.2f seconds", CoinCpuTime() - time1);
  handler_->message(CLP_GENERAL, messages_)
    << generalPrint
    << CoinMessageEol;
  primal(1);
  sprintf(generalPrint, "Total time %.2f seconds", CoinCpuTime() - time1);
  handler_->message(CLP_GENERAL, messages_)
    << generalPrint
    << CoinMessageEol;
  delete[] columnCounts;
  delete[] sol;
  delete[] lower;
  delete[] upper;
  delete[] whichBlock;
  delete[] when;
  delete[] columnAdd;
  delete[] rowAdd;
  delete[] elementAdd;
  delete[] objective;
  delete[] top;
  delete[] sub;
  return 0;
}
static ClpSimplex *deBound(ClpSimplex *oldModel)
{
  ClpSimplex *model = new ClpSimplex(*oldModel);
  int numberRows = model->numberRows();
  CoinPackedMatrix *matrix = model->matrix();
  const int *row = matrix->getIndices();
  const int *columnLength = matrix->getVectorLengths();
  const CoinBigIndex *columnStart = matrix->getVectorStarts();
  double *elementByColumn = matrix->getMutableElements();
  int numberColumns = model->numberColumns();
  double *rowLower = model->rowLower();
  double *rowUpper = model->rowUpper();
  double *columnLower = model->columnLower();
  double *columnUpper = model->columnUpper();
  double *objective = model->objective();
  double *change = new double[std::max(numberRows, numberColumns) + numberColumns];
  CoinBigIndex *rowStart = new CoinBigIndex[2 * numberColumns + 1];
  memset(change, 0, numberRows * sizeof(double));
  // first swap ones with infinite lower bounds
  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (columnLower[iColumn] == -COIN_DBL_MAX && columnUpper[iColumn] != COIN_DBL_MAX) {
      for (CoinBigIndex j = columnStart[iColumn]; j < columnStart[iColumn]
             + columnLength[iColumn];
           j++)
        elementByColumn[j] *= -1.0;
      objective[iColumn] *= -1.0;
      columnLower[iColumn] = -columnUpper[iColumn];
      columnUpper[iColumn] = COIN_DBL_MAX;
    }
  }
  // Out nonzero LB's
  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (columnLower[iColumn]) {
      double value = columnLower[iColumn];
      for (CoinBigIndex j = columnStart[iColumn]; j < columnStart[iColumn]
             + columnLength[iColumn];
           j++) {
        int iRow = row[j];
        change[iRow] -= value * elementByColumn[j];
      }
    }
  }
  for (int iRow = 0; iRow < numberRows; iRow++) {
    double value = change[iRow];
    if (rowLower[iRow] > -COIN_DBL_MAX)
      rowLower[iRow] -= value;
    if (rowUpper[iRow] < COIN_DBL_MAX)
      rowUpper[iRow] -= value;
  }
  int nExtra = 0;
  int *columnNew = reinterpret_cast< int * >(rowStart + numberColumns + 1);
  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (columnUpper[iColumn] < COIN_DBL_MAX && columnUpper[iColumn]) {
      columnNew[nExtra] = iColumn;
      change[nExtra++] = columnUpper[iColumn];
      columnUpper[iColumn] = COIN_DBL_MAX;
    }
  }
  double *elementNew = change + numberColumns;
  for (int i = 0; i < nExtra; i++) {
    rowStart[i] = i;
    elementNew[i] = 1.0;
  }
  rowStart[nExtra] = nExtra;
  model->addRows(nExtra, NULL, change,
    rowStart, columnNew, elementNew);
  delete[] rowStart;
  delete[] change;
  return model;
}
// Solve using Benders decomposition and maybe in parallel
int ClpSimplex::solveBenders(CoinStructuredModel *model, ClpSolve &options)
{
  double time1 = CoinCpuTime();
  //ClpSimplex * xxxx = deBound(this);
  //xxxx->writeMps("nobounds.mps");
  //delete xxxx;
  //int numberColumns = model->numberColumns();
  int numberRowBlocks = model->numberRowBlocks();
  int numberColumnBlocks = model->numberColumnBlocks();
  int numberElementBlocks = model->numberElementBlocks();
  char generalPrint[200];
  // We already have top level structure
  CoinModelBlockInfo *blockInfo = new CoinModelBlockInfo[numberElementBlocks];
  for (int i = 0; i < numberElementBlocks; i++) {
    CoinModel *thisBlock = model->coinBlock(i);
    assert(thisBlock);
    // just fill in info
    CoinModelBlockInfo info = CoinModelBlockInfo();
    int whatsSet = thisBlock->whatIsSet();
    info.matrix = static_cast< char >(((whatsSet & 1) != 0) ? 1 : 0);
    info.rhs = static_cast< char >(((whatsSet & 2) != 0) ? 1 : 0);
    info.rowName = static_cast< char >(((whatsSet & 4) != 0) ? 1 : 0);
    info.integer = static_cast< char >(((whatsSet & 32) != 0) ? 1 : 0);
    info.bounds = static_cast< char >(((whatsSet & 8) != 0) ? 1 : 0);
    info.columnName = static_cast< char >(((whatsSet & 16) != 0) ? 1 : 0);
    // Which block
    int iRowBlock = model->rowBlock(thisBlock->getRowBlock());
    info.rowBlock = iRowBlock;
    int iColumnBlock = model->columnBlock(thisBlock->getColumnBlock());
    info.columnBlock = iColumnBlock;
    blockInfo[i] = info;
  }
  // make up problems
  int numberBlocks = numberColumnBlocks - 1;
  // Find master columns and rows
  int *columnCounts = new int[numberColumnBlocks];
  CoinZeroN(columnCounts, numberColumnBlocks);
  int *rowCounts = new int[numberRowBlocks + 1];
  CoinZeroN(rowCounts, numberRowBlocks);
  int iBlock;
  for (iBlock = 0; iBlock < numberElementBlocks; iBlock++) {
    int iColumnBlock = blockInfo[iBlock].columnBlock;
    columnCounts[iColumnBlock]++;
    int iRowBlock = blockInfo[iBlock].rowBlock;
    rowCounts[iRowBlock]++;
  }
  int *whichBlock = new int[numberElementBlocks];
  int masterColumnBlock = -1;
  for (iBlock = 0; iBlock < numberColumnBlocks; iBlock++) {
    if (columnCounts[iBlock] > 1) {
      if (masterColumnBlock == -1) {
        masterColumnBlock = iBlock;
      } else {
        // Can't decode
        masterColumnBlock = -2;
        break;
      }
    }
  }
  int masterRowBlock = -1;
  int kBlock = 0;
  for (iBlock = 0; iBlock < numberRowBlocks; iBlock++) {
    int count = rowCounts[iBlock];
    rowCounts[iBlock] = kBlock;
    kBlock += count;
  }
  for (iBlock = 0; iBlock < numberElementBlocks; iBlock++) {
    int iRowBlock = blockInfo[iBlock].rowBlock;
    whichBlock[rowCounts[iRowBlock]] = iBlock;
    rowCounts[iRowBlock]++;
  }
  for (iBlock = numberRowBlocks - 1; iBlock >= 0; iBlock--)
    rowCounts[iBlock + 1] = rowCounts[iBlock];
  rowCounts[0] = 0;
  for (iBlock = 0; iBlock < numberRowBlocks; iBlock++) {
    int count = rowCounts[iBlock + 1] - rowCounts[iBlock];
    if (count == 1) {
      int kBlock = whichBlock[rowCounts[iBlock]];
      int iColumnBlock = blockInfo[kBlock].columnBlock;
      if (iColumnBlock == masterColumnBlock) {
        if (masterRowBlock == -1) {
          masterRowBlock = iBlock;
        } else {
          // Can't decode
          masterRowBlock = -2;
          break;
        }
      }
    }
  }
  if (masterColumnBlock < 0 || masterRowBlock == -2) {
    // What now
    abort();
  }
  delete[] columnCounts;
  // create all data
  const CoinPackedMatrix **first = new const CoinPackedMatrix *[numberRowBlocks];
  ClpSimplex *sub = new ClpSimplex[numberBlocks];
  ClpSimplex masterModel;
  // Set offset
  masterModel.setObjectiveOffset(model->objectiveOffset());
  kBlock = 0;
#define ADD_ARTIFICIALS
#ifdef ADD_ARTIFICIALS
  int *originalSubColumns = new int[numberBlocks];
  for (iBlock = 0; iBlock < numberBlocks; iBlock++)
    originalSubColumns[iBlock] = 9999999;
#endif
  int masterBlock = -1;
  for (iBlock = 0; iBlock < numberRowBlocks; iBlock++) {
    first[kBlock] = NULL;
    int start = rowCounts[iBlock];
    int end = rowCounts[iBlock + 1];
    assert(end - start <= 2);
    for (int j = start; j < end; j++) {
      int jBlock = whichBlock[j];
      int iColumnBlock = blockInfo[jBlock].columnBlock;
      int iRowBlock = blockInfo[jBlock].rowBlock;
      assert(iRowBlock == iBlock);
      if (iRowBlock != masterRowBlock && iColumnBlock == masterColumnBlock) {
        // first matrix
        first[kBlock] = model->coinBlock(jBlock)->packedMatrix();
      } else {
        const CoinPackedMatrix *matrix
          = model->coinBlock(jBlock)->packedMatrix();
        // Get pointers to arrays
        const double *columnLower;
        const double *columnUpper;
        const double *rowLower;
        const double *rowUpper;
        const double *objective;
        model->block(iRowBlock, iColumnBlock, rowLower, rowUpper,
          columnLower, columnUpper, objective);
        if (iRowBlock != masterRowBlock) {
          // diagonal block
          sub[kBlock].loadProblem(*matrix, columnLower, columnUpper,
            objective, rowLower, rowUpper);
          if (optimizationDirection_ < 0.0) {
            double *obj = sub[kBlock].objective();
            int n = sub[kBlock].numberColumns();
            for (int i = 0; i < n; i++)
              obj[i] = -obj[i];
          }
          if (this->factorizationFrequency() == 200) {
            // User did not touch preset
            sub[kBlock].defaultFactorizationFrequency();
          } else {
            // make sure model has correct value
            sub[kBlock].setFactorizationFrequency(this->factorizationFrequency());
          }
          sub[kBlock].setPerturbation(50);
#ifdef ADD_ARTIFICIALS
          originalSubColumns[kBlock] = sub[kBlock].numberColumns();
          if (intParam_[0] < 1000000) {
            printf("** Adding artificials\n");
            int nRow = sub[kBlock].numberRows();
            CoinBigIndex *addStarts = new CoinBigIndex[2 * nRow + 1];
            int *addRow = new int[2 * nRow];
            double *addElement = new double[2 * nRow];
            addStarts[0] = 0;
            int numberArtificials = 0;
            double penalty = 1.0e5;
            double *addCost = new double[2 * nRow];
            const double *lower = sub[kBlock].rowLower();
            const double *upper = sub[kBlock].rowUpper();
            for (int iRow = 0; iRow < nRow; iRow++) {
              if (lower[iRow] > -1.0e20) {
                addRow[numberArtificials] = iRow;
                addElement[numberArtificials] = 1.0;
                addCost[numberArtificials] = penalty;
                numberArtificials++;
                addStarts[numberArtificials] = numberArtificials;
              }
              if (upper[iRow] < 1.0e20) {
                addRow[numberArtificials] = iRow;
                addElement[numberArtificials] = -1.0;
                addCost[numberArtificials] = penalty;
                numberArtificials++;
                addStarts[numberArtificials] = numberArtificials;
              }
            }
            if (numberArtificials) {
              sub[kBlock].addColumns(numberArtificials, NULL, NULL, addCost,
                addStarts, addRow, addElement);
            }
            delete[] addStarts;
            delete[] addRow;
            delete[] addElement;
            delete[] addCost;
          }
#endif
          // Set rowCounts to be diagonal block index for cleanup
          rowCounts[kBlock] = jBlock;
        } else {
          // master
          masterBlock = jBlock;
          masterModel.loadProblem(*matrix, columnLower, columnUpper,
            objective, rowLower, rowUpper);
          if (optimizationDirection_ < 0.0) {
            double *obj = masterModel.objective();
            int n = masterModel.numberColumns();
            for (int i = 0; i < n; i++)
              obj[i] = -obj[i];
          }
        }
      }
    }
    if (iBlock != masterRowBlock)
      kBlock++;
  }
  delete[] whichBlock;
  delete[] blockInfo;
  assert(masterBlock >= 0);
  int numberMasterColumns = masterModel.numberColumns();
  masterModel.setStrParam(ClpProbName, "Master");
  // Overkill in terms of space
  int spaceNeeded = std::max(numberBlocks * (numberMasterColumns + 1),
    2 * numberMasterColumns);
  CoinBigIndex *columnAdd = new CoinBigIndex[spaceNeeded];
  int *indexColumnAdd = reinterpret_cast< int * >(columnAdd);
  double *elementAdd = new double[spaceNeeded];
  spaceNeeded = numberBlocks;
  CoinBigIndex *rowAdd = new CoinBigIndex[2 * spaceNeeded + 1]; // temp for block info
  int *blockPrint = reinterpret_cast< int * >(rowAdd + spaceNeeded + 1);
  double *objective = new double[spaceNeeded];
  int logLevel = handler_->logLevel();
  //#define TEST_MODEL
#ifdef TEST_MODEL
  double goodValue = COIN_DBL_MAX;
  ClpSimplex goodModel;
  if (logLevel > 3) {
    // temp - create copy with master at front
    const int *columnBack = model->coinBlock(masterBlock)->originalColumns();
    int numberColumns2 = model->coinBlock(masterBlock)->numberColumns();
    const int *rowBack = model->coinBlock(masterBlock)->originalRows();
    int numberRows2 = model->coinBlock(masterBlock)->numberRows();
    int *whichColumn = new int[numberColumns_ + numberBlocks];
    int *whichRow = new int[numberRows_];
    int nColumn = 0;
    for (int iColumn = 0; iColumn < numberColumns2; iColumn++) {
      int kColumn = columnBack[iColumn];
      whichColumn[nColumn++] = kColumn;
    }
    int nRow = 0;
    for (int iRow = 0; iRow < numberRows2; iRow++) {
      int kRow = rowBack[iRow];
      whichRow[nRow++] = kRow;
    }
    for (iBlock = 0; iBlock < numberBlocks; iBlock++) {
      int kBlock = rowCounts[iBlock];
      const int *columnBack = model->coinBlock(kBlock)->originalColumns();
      int numberColumns2 = model->coinBlock(kBlock)->numberColumns();
      const int *rowBack = model->coinBlock(kBlock)->originalRows();
      int numberRows2 = model->coinBlock(kBlock)->numberRows();
      for (int iColumn = 0; iColumn < numberColumns2; iColumn++) {
        int kColumn = columnBack[iColumn];
        whichColumn[nColumn++] = kColumn;
      }
      for (int iRow = 0; iRow < numberRows2; iRow++) {
        int kRow = rowBack[iRow];
        whichRow[nRow++] = kRow;
      }
    }
    ClpSimplex temp(this, nRow, whichRow, nColumn, whichColumn);
    temp.writeMps("ordered.mps");
    for (int i = numberMasterColumns; i < numberColumns_; i++)
      whichColumn[i + numberBlocks] = i;
    for (int i = 0; i < numberMasterColumns; i++)
      whichColumn[i] = i;
    for (int i = 0; i < numberBlocks; i++)
      whichColumn[i + numberMasterColumns] = i + numberColumns_;
    double *lower = new double[numberBlocks];
    // Add extra variables
    columnAdd[0] = 0;
    for (int iBlock = 0; iBlock < numberBlocks; iBlock++) {
      objective[iBlock] = 0.0;
      lower[iBlock] = -COIN_DBL_MAX;
      columnAdd[iBlock + 1] = 0;
    }
    temp.addColumns(numberBlocks, lower, NULL, objective,
      columnAdd, NULL, NULL);
    delete[] lower;
    ClpSimplex temp2(&temp, nRow, whichRow, nColumn + numberBlocks, whichColumn);
    goodModel = temp2;
    goodModel.dual();
    goodValue = goodModel.objectiveValue();
    delete[] whichColumn;
    delete[] whichRow;
  }
#endif
  int maxPass = options.independentOption(2);
  if (maxPass < 2)
    maxPass = 100;
  int iPass;
  double lastObjective = -1.0e31;
  // Create columns for proposals
  int numberMasterRows = masterModel.numberRows();
  //masterModel.resize(numberMasterColumns + numberBlocks, numberMasterRows);
  if (this->factorizationFrequency() == 200) {
    // User did not touch preset
    masterModel.defaultFactorizationFrequency();
  } else {
    // make sure model has correct value
    masterModel.setFactorizationFrequency(this->factorizationFrequency());
  }
  masterModel.setPerturbation(50);
  // temp bounds
  if (0) {
    printf("temp bounds\n");
    double *lower = masterModel.columnLower();
    double *upper = masterModel.columnUpper();
    for (int i = 0; i < numberMasterColumns; i++) {
      lower[i] = std::max(lower[i], -1.0e8);
      upper[i] = std::min(upper[i], 1.0e8);
    }
  }
  //printf("take out debound\n");
  //master=*deBound(&master);
  // Arrays to say which block and when created
  int maximumRows = 2 * numberMasterColumns + 10 * numberBlocks;
  whichBlock = new int[maximumRows];
  int *when = new int[maximumRows];
  // state of each problem (0 first or infeas, 1 feasible, add 2 if extra variables freed)
  int *problemState = new int[numberBlocks];
  memset(problemState, 0, numberBlocks * sizeof(int));
  int numberRowsGenerated = numberBlocks;
  // space for rhs modifications
  double **modification = new double *[numberBlocks];
  // Add extra variables
  columnAdd[0] = 0;
  for (int iBlock = 0; iBlock < numberBlocks; iBlock++) {
    int numberRows2 = sub[iBlock].numberRows();
    int numberColumns2 = sub[iBlock].numberColumns();
    int numberTotal2 = numberRows2 + numberColumns2;
    modification[iBlock] = new double[2 * numberTotal2 + numberRows2];
    double *save = modification[iBlock];
    memcpy(save, sub[iBlock].rowLower(), numberRows2 * sizeof(double));
    save += numberRows2;
    memcpy(save, sub[iBlock].columnLower(), numberColumns2 * sizeof(double));
    save += numberColumns2;
    memcpy(save, sub[iBlock].rowUpper(), numberRows2 * sizeof(double));
    save += numberRows2;
    memcpy(save, sub[iBlock].columnUpper(), numberColumns2 * sizeof(double));
    objective[iBlock] = 1.0;
    columnAdd[iBlock + 1] = 0;
    when[iBlock] = -1;
    whichBlock[iBlock] = iBlock;
    if (logLevel < 2) {
      // modify printing
      sub[iBlock].messagesPointer()->setDetailMessage(2, 6);
      sub[iBlock].messagesPointer()->setDetailMessage(2, 14);
      sub[iBlock].messagesPointer()->setDetailMessage(2, 0);
    }
    // scaling
    sub[iBlock].scaling(scalingFlag_);
    //sub[iBlock].scaling(0);
  }
  masterModel.addColumns(numberBlocks, NULL, NULL, objective,
    columnAdd, NULL, NULL);
  sprintf(generalPrint, "Time to decompose %.2f seconds", CoinCpuTime() - time1);
  handler_->message(CLP_GENERAL, messages_)
    << generalPrint
    << CoinMessageEol;
  int ixxxxxx = 0;
  if (intParam_[0] >= 100000 && intParam_[0] < 100999) {
    ixxxxxx = intParam_[0] - 100000;
    printf("ixxxxxx %d\n", ixxxxxx);
  }
#ifdef ADD_ARTIFICIALS
  double **saveObjective = new double *[numberBlocks];
  for (int i = 0; i < numberBlocks; i++)
    saveObjective[i] = CoinCopyOfArray(sub[i].objective(),
      originalSubColumns[i]);
  int iFudge = 1;
  int lowFudge = 0;
  int highFudge = 0;
#endif
#define UNBOUNDED
  if (ixxxxxx > 0) {
    for (iBlock = 0; iBlock < std::min(numberBlocks, ixxxxxx); iBlock++) {
      ClpSimplex *temp = deBound(sub + iBlock);
      sub[iBlock] = *temp;
      delete temp;
      delete[] modification[iBlock];
      int numberRows2 = sub[iBlock].numberRows();
      int numberColumns2 = sub[iBlock].numberColumns();
      int numberTotal2 = numberRows2 + numberColumns2;
      modification[iBlock] = new double[2 * numberTotal2 + numberRows2];
      double *save = modification[iBlock];
      memcpy(save, sub[iBlock].rowLower(), numberRows2 * sizeof(double));
      save += numberRows2;
      memcpy(save, sub[iBlock].columnLower(), numberColumns2 * sizeof(double));
      save += numberColumns2;
      memcpy(save, sub[iBlock].rowUpper(), numberRows2 * sizeof(double));
      save += numberRows2;
      memcpy(save, sub[iBlock].columnUpper(), numberColumns2 * sizeof(double));
    }
  }
  double treatSubAsFeasible = 1.0e-6;
  int numberSubInfeasible = 0;
  bool canSkipSubSolve = false;
  int numberProposals = 999;
  for (iPass = 0; iPass < maxPass; iPass++) {
    sprintf(generalPrint, "Start of pass %d", iPass);
    handler_->message(CLP_GENERAL, messages_)
      << generalPrint
      << CoinMessageEol;
    // Solve master - may be unbounded
    //masterModel.scaling(0);
    // get obj for debug
    //double objSum = masterModel.objectiveValue();
    //for (int i = 0; i < numberBlocks; i++)
    //  objSum += sub[i].objectiveValue();
    //printf("objsum %g\n",objSum);
    if (0) {
      masterModel.writeMps("yy.mps");
      masterModel.writeBasis("yy.bas", true, 2);
    }
    {
      // free up extra variables
      double *lower = masterModel.columnLower();
      int numberFreed = 0;
      for (int i = 0; i < numberBlocks; i++) {
        if (problemState[i] == 1) {
          // ? need trust region ?
          lower[i + numberMasterColumns] = -COIN_DBL_MAX;
          //if (problemState[i]!=2) {
          numberFreed++;
          problemState[i] = 3;
          //}
        }
      }
      if (numberFreed)
        lastObjective = -1.0e31;
    }
#ifdef TRY_NO_SCALING
    masterModel.scaling(0);
#endif
    masterModel.dual();
    if ((maxPass == 5000 && scalingFlag_) || (maxPass == 4000 && !scalingFlag_)) {
      int n = masterModel.numberIterations();
      masterModel.scaling(0);
      masterModel.primal();
      masterModel.setNumberIterations(n + masterModel.numberIterations());
      masterModel.scaling(scalingFlag_);
    }
    int masterStatus = masterModel.status(); // do here as can change
    sprintf(generalPrint, "Pass %d objective %g change %g",
      iPass, masterModel.objectiveValue(),
      masterModel.objectiveValue() - lastObjective);
    handler_->message(CLP_GENERAL, messages_)
      << generalPrint
      << CoinMessageEol;
#ifndef UNBOUNDED
    if (masterStatus == 2) {
      // unbounded
      masterModel.writeMps("unbounded.mps");
      // get primal feasible
      masterModel.primal();
      const double *fullLower = columnLower();
      const double *fullUpper = columnUpper();
      double *lower = masterModel.columnLower();
      double *upper = masterModel.columnUpper();
      double *solution = masterModel.primalColumnSolution();
      const int *columnBack = model->coinBlock(masterBlock)->originalColumns();
      if (!numberTimesUnbounded) {
        for (int iColumn = 0; iColumn < numberMasterColumns; iColumn++) {
          int kColumn = columnBack[iColumn];
          double value = solution[iColumn];
          double lowerValue = std::max(fullLower[kColumn],
            std::min(value, fullUpper[kColumn]) - trust);
          lower[iColumn] = lowerValue;
          double upperValue = std::min(fullUpper[kColumn],
            std::max(value, fullLower[kColumn]) + trust);
          upper[iColumn] = upperValue;
        }
#ifdef TEST_MODEL
        if (logLevel > 3) {
          //use ones from solved
          const double *solutionGood = goodModel.primalColumnSolution();
          for (int iColumn = 0; iColumn < numberMasterColumns; iColumn++) {
            double value = solutionGood[iColumn];
            lower[iColumn] = std::min(value, -trust);
            upper[iColumn] = std::max(value, trust);
          }
        }
#endif
      } else {
        abort(); // probably can happen
      }
      numberTimesUnbounded++;
      masterModel.dual();
      masterStatus = masterModel.status();
      assert(!masterStatus);
      masterModel.setNumberIterations(1); // so will continue
    }
    if (numberTimesUnbounded > 1 && trust > 0.0) {
      assert(!masterStatus);
      const double *fullLower = columnLower();
      const double *fullUpper = columnUpper();
      double *lower = masterModel.columnLower();
      double *upper = masterModel.columnUpper();
      //double * solution = masterModel.primalColumnSolution();
      const int *columnBack = model->coinBlock(masterBlock)->originalColumns();
      int nTrusted = 0;
      for (int iColumn = 0; iColumn < numberMasterColumns; iColumn++) {
        int kColumn = columnBack[iColumn];
        if (lower[iColumn] > fullLower[kColumn]) {
          if (masterModel.getColumnStatus(iColumn) != atLowerBound) {
            lower[iColumn] = fullLower[kColumn];
          } else {
            nTrusted++;
          }
        }
        if (upper[iColumn] < fullUpper[kColumn]) {
          if (masterModel.getColumnStatus(iColumn) != atUpperBound) {
            upper[iColumn] = fullUpper[kColumn];
          } else {
            nTrusted++;
          }
        }
      }
      if (nTrusted) {
        sprintf(generalPrint, "%d at artificial bound", nTrusted);
      } else {
        sprintf(generalPrint, "All at natural bounds");
        trust = 0.0;
      }
      handler_->message(CLP_GENERAL2, messages_)
        << generalPrint
        << CoinMessageEol;
    }
#endif
    if (!masterStatus) {
      if (masterModel.numberIterations() == 0 && iPass) {
        if ((!numberSubInfeasible && !numberProposals) || treatSubAsFeasible > 1.0e-2 || iPass > 5555)
          break; // finished
        if (!numberProposals && numberSubInfeasible) {
          treatSubAsFeasible *= 2.0;
          printf("Doubling sub primal tolerance to %g\n", treatSubAsFeasible);
        } else {
          treatSubAsFeasible *= 1.2;
          printf("Increasing sub primal tolerance to %g\n", treatSubAsFeasible);
        }
        canSkipSubSolve = false;
      } else if (!numberSubInfeasible) {
        if (treatSubAsFeasible > 1.0e-6) {
          treatSubAsFeasible = std::max(0.9 * treatSubAsFeasible, 1.0e-6);
          printf("Reducing sub primal tolerance to %g\n", treatSubAsFeasible);
        }
      }
      if (masterModel.objectiveValue() < lastObjective + 1.0e-7 && iPass > 5555)
        break; // finished
      lastObjective = masterModel.objectiveValue();
    }
    // mark non-basic rows and delete if necessary
    int iRow;
    numberRowsGenerated = masterModel.numberRows() - numberMasterRows;
    for (iRow = 0; iRow < numberRowsGenerated; iRow++) {
      if (masterModel.getStatus(iRow + numberMasterRows) != ClpSimplex::basic)
        when[iRow] = iPass;
    }
    if (numberRowsGenerated > maximumRows - numberBlocks) {
      // delete
      int numberKeep = 0;
      int numberDelete = 0;
      int *whichDelete = new int[numberRowsGenerated];
      for (iRow = 0; iRow < numberRowsGenerated; iRow++) {
        if (masterModel.getRowStatus(iRow + numberMasterRows) != basic) {
          // keep
          when[numberKeep] = when[iRow];
          whichBlock[numberKeep++] = whichBlock[iRow];
        } else {
          // delete
          whichDelete[numberDelete++] = iRow + numberMasterRows;
        }
      }
      if (numberRowsGenerated - numberDelete > maximumRows - numberBlocks) {
        for (iRow = 0; iRow < numberRowsGenerated; iRow++) {
          if (when[iRow] > iPass - 7) {
            // keep
            when[numberKeep] = when[iRow];
            whichBlock[numberKeep++] = whichBlock[iRow];
          } else {
            // delete
            whichDelete[numberDelete++] = iRow + numberMasterRows;
          }
        }
      }
      numberRowsGenerated -= numberDelete;
      masterModel.deleteRows(numberDelete, whichDelete);
      delete[] whichDelete;
    }
    double *primal = NULL;
    bool deletePrimal = false;
    if (masterStatus == 0) {
      primal = masterModel.primalColumnSolution();
    } else if (masterStatus == 2 && masterModel.numberRows()) {
      // scale back ray (1.0e20?)
      primal = masterModel.ray();
      //deletePrimal = true;
      sprintf(generalPrint, "The sum of dual infeasibilities is %g",
        masterModel.sumDualInfeasibilities());
      handler_->message(CLP_GENERAL, messages_)
        << generalPrint
        << CoinMessageEol;
    } else if (!masterModel.numberRows()) {
      assert(!iPass);
      primal = masterModel.primalColumnSolution();
      memset(masterModel.primalColumnSolution(),
        0, numberMasterColumns * sizeof(double));
    } else {
      printf("Master infeasible - sum %g\n",
        masterModel.sumPrimalInfeasibilities());
      masterModel.setProblemStatus(0);
      primal = masterModel.primalColumnSolution();
      masterModel.writeMps("inf.mps");
      //abort();
    }
#ifndef UNBOUNDED
    if (masterStatus == 2) {
      // adjust variables with no elements
      const int *columnLength = masterModel.matrix()->getVectorLengths();
      const double *lower = masterModel.columnLower();
      const double *upper = masterModel.columnUpper();
      const double *obj = masterModel.objective();
      for (int i = 0; i < numberMasterColumns; i++) {
        double value = primal[i];
        if (!columnLength[i]) {
          if (obj[i] < 0.0)
            value += 1.0e10;
          else if (obj[i] > 0.0)
            value -= 1.0e10;
        }
        // make sure feasible
        primal[i] = std::max(-1.0e10, std::min(1.0e10, value));
        primal[i] = std::max(lower[i], std::min(upper[i], primal[i]));
      }
    }
#endif
    // Create rhs for sub problems and solve
    for (iBlock = 0; iBlock < numberBlocks; iBlock++) {
#ifdef ADD_ARTIFICIALS
      {
        double *columnUpper2 = sub[iBlock].columnUpper();
        double *obj = sub[iBlock].objective();
        int start = originalSubColumns[iBlock];
        int numberColumns2 = sub[iBlock].numberColumns();
        if (iFudge >= lowFudge && iFudge <= abs(highFudge)) {
          for (int i = start; i < numberColumns2; i++)
            columnUpper2[i] = COIN_DBL_MAX;
          memset(obj, 0, originalSubColumns[iBlock] * sizeof(double));
        } else {
          for (int i = start; i < numberColumns2; i++)
            columnUpper2[i] = 0.0;
          memcpy(obj, saveObjective[iBlock], originalSubColumns[iBlock] * sizeof(double));
          if (highFudge < 0) {
            sub[iBlock].allSlackBasis(true);
          }
        }
      }
#endif
      int numberRows2 = sub[iBlock].numberRows();
      int numberColumns2 = sub[iBlock].numberColumns();
      double *saveLower = modification[iBlock];
      double *lower2 = sub[iBlock].rowLower();
      memcpy(lower2, saveLower, numberRows2 * sizeof(double));
      double *saveColumnLower = saveLower + numberRows2;
      double *columnLower2 = sub[iBlock].columnLower();
      memcpy(columnLower2, saveColumnLower, numberColumns2 * sizeof(double));
      double *saveUpper = saveColumnLower + numberColumns2;
      double *upper2 = sub[iBlock].rowUpper();
      memcpy(upper2, saveUpper, numberRows2 * sizeof(double));
      double *saveColumnUpper = saveUpper + numberRows2;
      double *columnUpper2 = sub[iBlock].columnUpper();
      memcpy(columnUpper2, saveColumnUpper, numberColumns2 * sizeof(double));
      double *lastMod = saveColumnUpper + numberColumns2;
#ifdef UNBOUNDED
      if (masterStatus == 2) {
        for (int i = 0; i < numberRows2; i++) {
          if (lower2[i] > -COIN_DBL_MAX)
            lower2[i] = 0.0;
          if (upper2[i] < COIN_DBL_MAX)
            upper2[i] = 0.0;
        }
        for (int i = 0; i < numberColumns2; i++) {
          if (columnLower2[i] > -COIN_DBL_MAX)
            columnLower2[i] = 0.0;
          if (columnUpper2[i] < COIN_DBL_MAX)
            columnUpper2[i] = 0.0;
        }
      }
#endif
      // new rhs
      double *rhs = sub[iBlock].dualRowSolution();
      CoinZeroN(rhs, numberRows2);
      first[iBlock]->times(primal, rhs);
      for (int i = 0; i < numberRows2; i++) {
        double value = rhs[i];
        if (lower2[i] > -1.0e30)
          lower2[i] -= value;
        if (upper2[i] < 1.0e30)
          upper2[i] -= value;
      }
      bool canSkip = false;
      if (canSkipSubSolve) {
        canSkip = true;
        const double *rowSolution = sub[iBlock].primalRowSolution();
        for (int i = 0; i < numberRows2; i++) {
          double value = lastMod[i] - rhs[i];
          if (fabs(value) > primalTolerance_) {
            // see if we can adjust
            double rowValue = rowSolution[i];
            if (rowValue < lower2[i] - primalTolerance_
              || rowValue > upper2[i] + primalTolerance_) {
              canSkip = false;
              break;
            } else if (sub[iBlock].getRowStatus(i) != basic) {
              canSkip = false;
              break;
            }
          }
        }
      }
      if (!canSkip) {
        memcpy(lastMod, rhs, numberRows2 * sizeof(double));
      } else {
        // mark
        sub[iBlock].setSecondaryStatus(99);
      }
    }
    canSkipSubSolve = true;
    if (!iPass) {
      // do first and then copy status?
#ifdef TRY_NO_SCALING
      sub[0].scaling(0);
#endif
      sub[0].dual();
      if ((maxPass == 5000 && scalingFlag_) || (maxPass == 4000 && !scalingFlag_)) {
        int n = sub[0].numberIterations();
        sub[0].scaling(0);
        sub[0].primal();
        sub[0].setNumberIterations(n + sub[0].numberIterations());
        sub[0].scaling(scalingFlag_);
      }
      int numberIterations = sub[0].numberIterations();
      if (sub[0].problemStatus()) {
        sub[0].primal();
        numberIterations += sub[0].numberIterations();
      }
      sprintf(generalPrint, "First block - initial solve - %d iterations, objective %g",
        numberIterations, sub[0].objectiveValue());
      handler_->message(CLP_GENERAL, messages_)
        << generalPrint
        << CoinMessageEol;
      // copy status if same size
      int numberRows2 = sub[0].numberRows();
      int numberColumns2 = sub[0].numberColumns();
      for (int iBlock = 1; iBlock < numberBlocks; iBlock++) {
        int numberRows2a = sub[iBlock].numberRows();
        int numberColumns2a = sub[iBlock].numberColumns();
        if (numberRows2 == numberRows2a && numberColumns2 == numberColumns2a) {
          memcpy(sub[iBlock].primalColumnSolution(),
            sub[0].primalColumnSolution(),
            numberColumns2 * sizeof(double));
          memcpy(sub[iBlock].statusArray(), sub[0].statusArray(),
            numberRows2 + numberColumns2);
        }
      }
      // mark
      sub[0].setSecondaryStatus(99);
    }
      numberSubInfeasible = 0;
      for (iBlock = 0; iBlock < numberBlocks; iBlock++) {
#ifdef TRY_NO_SCALING
        sub[iBlock].scaling(0);
#endif
        if (sub[iBlock].secondaryStatus() != 99 || false) {
          //int ix=sub[iBlock].secondaryStatus();
          int lastStatus = sub[iBlock].problemStatus();
          // was do dual unless unbounded
          double saveTolerance = sub[iBlock].primalTolerance();
          if (lastStatus == 0 || !iPass) {
            //if (lastStatus<2||!iPass) {
            //sub[iBlock].dual();
            sub[iBlock].primal();
            if (!sub[iBlock].isProvenOptimal() && sub[iBlock].sumPrimalInfeasibilities() < treatSubAsFeasible) {
              printf("Block %d was feasible now has small infeasibility %g\n", iBlock,
                sub[iBlock].sumPrimalInfeasibilities());
              sub[iBlock].setPrimalTolerance(std::min(treatSubAsFeasible, 1.0e-4));
              sub[iBlock].setCurrentPrimalTolerance(std::min(treatSubAsFeasible, 1.0e-4));
              sub[iBlock].primal();
              sub[iBlock].setProblemStatus(0);
              problemState[iBlock] |= 4; // force actions
            }
            if ((maxPass == 5000 && scalingFlag_) || (maxPass == 4000 && !scalingFlag_)) {
              int n = sub[iBlock].numberIterations();
              sub[iBlock].scaling(0);
              sub[iBlock].primal();
              sub[iBlock].setNumberIterations(n + sub[iBlock].numberIterations());
              sub[iBlock].scaling(scalingFlag_);
            }
          } else if (lastStatus == 1) {
            // zero out objective
            double saveScale = sub[iBlock].infeasibilityCost();
            ClpObjective *saveObjective = sub[iBlock].objectiveAsObject();
            int numberColumns = sub[iBlock].numberColumns();
            ClpLinearObjective fake(NULL, numberColumns);
            sub[iBlock].setObjectivePointer(&fake);
            int saveOptions = sub[iBlock].specialOptions();
            sub[iBlock].setSpecialOptions(saveOptions | 8192);
            sub[iBlock].primal();
            if ((maxPass == 5000 && scalingFlag_) || (maxPass == 4000 && !scalingFlag_)) {
              int n = sub[iBlock].numberIterations();
              sub[iBlock].scaling(0);
              sub[iBlock].primal();
              sub[iBlock].setNumberIterations(n + sub[iBlock].numberIterations());
              sub[iBlock].scaling(scalingFlag_);
            }
            sub[iBlock].setObjectivePointer(saveObjective);
            sub[iBlock].setInfeasibilityCost(saveScale);
            if (!sub[iBlock].isProvenOptimal() && sub[iBlock].sumPrimalInfeasibilities() < treatSubAsFeasible) {
              printf("Block %d was infeasible now has small infeasibility %g\n", iBlock,
                sub[iBlock].sumPrimalInfeasibilities());
              sub[iBlock].setProblemStatus(0);
              sub[iBlock].setPrimalTolerance(std::min(treatSubAsFeasible, 1.0e-4));
              sub[iBlock].setCurrentPrimalTolerance(std::min(treatSubAsFeasible, 1.0e-4));
            }
            if (sub[iBlock].isProvenOptimal()) {
              sub[iBlock].primal();
              if ((maxPass == 5000 && scalingFlag_) || (maxPass == 4000 && !scalingFlag_)) {
                int n = sub[iBlock].numberIterations();
                sub[iBlock].scaling(0);
                sub[iBlock].primal();
                sub[iBlock].setNumberIterations(n + sub[iBlock].numberIterations());
                sub[iBlock].scaling(scalingFlag_);
              }
              if (!sub[iBlock].isProvenOptimal()) {
                printf("Block %d infeasible on second go has small infeasibility %g\n", iBlock,
                  sub[iBlock].sumPrimalInfeasibilities());
                sub[iBlock].setProblemStatus(0);
              }
              problemState[iBlock] |= 4; // force actions
            } else {
              printf("Block %d still infeasible - sum %g - %d iterations\n", iBlock,
                sub[iBlock].sumPrimalInfeasibilities(),
                sub[iBlock].numberIterations());
              numberSubInfeasible++;
              if (!sub[iBlock].ray()) {
                printf("Block %d has no ray!\n", iBlock);
                sub[iBlock].primal();
                assert(sub[iBlock].ray()); // otherwise declare optimal
              }
            }
            sub[iBlock].setSpecialOptions(saveOptions);
          } else {
            sub[iBlock].primal();
          }
          sub[iBlock].setPrimalTolerance(saveTolerance);
          sub[iBlock].setCurrentPrimalTolerance(saveTolerance);
          if (!sub[iBlock].isProvenOptimal() && !sub[iBlock].isProvenPrimalInfeasible()) {
            printf("!!!Block %d has bad status %d\n", iBlock, sub[iBlock].problemStatus());
            sub[iBlock].primal(); // last go
          }
          //#define WRITE_ALL
#ifdef WRITE_ALL
          char name[20];
          sprintf(name, "pass_%d_block_%d.mps", iPass, iBlock);
          sub[iBlock].writeMps(name);
          sprintf(name, "pass_%d_block_%d.bas", iPass, iBlock);
          sub[iBlock].writeBasis(name, true);
          if (sub[iBlock].problemStatus() == 1) {
            sub[iBlock].readBasis(name);
            sub[iBlock].primal();
          }
#endif
          //assert (!sub[iBlock].numberIterations()||ix!=99);
        }
      }
    for (iBlock = 0; iBlock < numberBlocks; iBlock++) {
      if (!iPass)
        problemState[iBlock] |= 4; // force actions
      // if state changed then fake number of iterations
      if ((problemState[iBlock] & 1) == 0) {
        // was infeasible
        if (sub[iBlock].isProvenOptimal()) {
          // say feasible and changed
          problemState[iBlock] |= 1 + 4;
        }
      } else {
        // was feasible
        if (sub[iBlock].isProvenPrimalInfeasible()) {
          // say infeasible and changed
          problemState[iBlock] &= ~1;
          problemState[iBlock] |= 4;
        }
      }
      if (sub[iBlock].secondaryStatus() != 99) {
        if (logLevel > 1) {
          sprintf(generalPrint, "Block %d - %d iterations, objective %g",
            iBlock, sub[iBlock].numberIterations(),
            sub[iBlock].objectiveValue());
          handler_->message(CLP_GENERAL2, messages_)
            << generalPrint
            << CoinMessageEol;
        }
        if (sub[iBlock].problemStatus() && sub[iBlock].algorithm() < 0 && false) {
          int numberRows2 = sub[iBlock].numberRows();
          double *ray = sub[iBlock].infeasibilityRay();
          double *saveRay = CoinCopyOfArray(ray, numberRows2);
          double *obj = sub[iBlock].objective();
          int numberColumns2 = sub[iBlock].numberColumns();
          double *saveObj = CoinCopyOfArray(obj, numberColumns2);
          memset(obj, 0, numberColumns2 * sizeof(double));
          sub[iBlock].allSlackBasis(true);
          sub[iBlock].primal();
          memcpy(obj, saveObj, numberColumns2 * sizeof(double));
          delete[] saveObj;
          ray = sub[iBlock].infeasibilityRay();
          for (int i = 0; i < numberRows2; i++) {
            if (fabs(ray[i] - saveRay[i]) > 1.0e-4 + 1.0e20)
              printf("** diffray block %d row %d first %g second %g\n",
                iBlock, i, saveRay[i], ray[i]);
          }
          delete[] saveRay;
        }
      }
    }
    if (!iPass)
      sub[0].setSecondaryStatus(0);
    if (logLevel > 2) {
      for (iBlock = 0; iBlock < numberBlocks; iBlock++) {
        printf("block %d obj %g thetaC %g\n", iBlock, sub[iBlock].objectiveValue(), masterModel.primalColumnSolution()[numberMasterColumns + iBlock]);
      }
    }
    rowAdd[0] = 0;
    numberProposals = 0;
    for (iBlock = 0; iBlock < numberBlocks; iBlock++) {
      int numberRows2 = sub[iBlock].numberRows();
      int numberColumns2 = sub[iBlock].numberColumns();
      double *saveLower = modification[iBlock];
      //double *lower2 = sub[iBlock].rowLower();
      double *saveUpper = saveLower + numberRows2 + numberColumns2;
      //double *upper2 = sub[iBlock].rowUpper();
      int typeRun = sub[iBlock].secondaryStatus();
      sub[iBlock].setSecondaryStatus(0);
      if (typeRun != 99) {
        if (0) {
          //double objValue = 0.0;
          const double *solution = sub[iBlock].dualRowSolution();
          for (int i = 0; i < numberRows2; i++) {
            if (solution[i] < -dualTolerance_) {
              // at upper
              assert(saveUpper[i] < 1.0e30);
              //objValue += solution[i] * upper2[i];
            } else if (solution[i] > dualTolerance_) {
              // at lower
              assert(saveLower[i] > -1.0e30);
              //objValue += solution[i] * lower2[i];
            }
          }
          //printf("obj %g\n",objValue);
        }
        // temp
        if (sub[iBlock].isProvenPrimalInfeasible() && !sub[iBlock].numberIterations())
          problemState[iBlock] |= 4;
        // get proposal
        if (sub[iBlock].numberIterations() || (problemState[iBlock] & 4) != 0) {
          double objValue = 0.0;
          int start = static_cast< int >(rowAdd[numberProposals]);
          // proposal
          if (sub[iBlock].isProvenOptimal()) {
            double *solution = sub[iBlock].dualRowSolution();
            first[iBlock]->transposeTimes(solution, elementAdd + start);
            for (int i = 0; i < numberRows2; i++) {
              if (sub[iBlock].getRowStatus(i) == basic)
                solution[i] = 0.0;
              if (saveUpper[i] > saveLower[i]) {
                if (sub[iBlock].getRowStatus(i) == atUpperBound)
                  objValue += solution[i] * saveUpper[i];
                else
                  objValue += solution[i] * saveLower[i];
              } else {
                // fixed
                objValue += solution[i] * saveLower[i];
              }
            }
            const double *dj = sub[iBlock].dualColumnSolution();
            const double *columnLower = sub[iBlock].columnLower();
            const double *columnUpper = sub[iBlock].columnUpper();
            double objValue2 = 0.0;
            int numberColumns2 = sub[iBlock].numberColumns();
            for (int i = 0; i < numberColumns2; i++) {
              if (logLevel > 2) {
                if (sub[iBlock].getColumnStatus(i) != basic && fabs(sub[iBlock].primalColumnSolution()[i]) > 1.0e-5)
                  printf("zz %d has value %g\n", i, sub[iBlock].primalColumnSolution()[i]);
              }
              if (sub[iBlock].getColumnStatus(i) == isFixed) {
                objValue2 += columnLower[i] * dj[i];
              } else if (sub[iBlock].getColumnStatus(i) == atLowerBound) {
                objValue2 += columnLower[i] * dj[i];
              } else if (sub[iBlock].getColumnStatus(i) == atUpperBound) {
                objValue2 += columnUpper[i] * dj[i];
              }
            }
#if 1
            double objValue3 = 0.0;
            const double *cost = sub[iBlock].objective();
            for (int i = 0; i < numberColumns2; i++) {
              double value = dj[i] - cost[i];
              if (sub[iBlock].getColumnStatus(i) == isFixed) {
                objValue3 += columnLower[i] * value;
              } else if (sub[iBlock].getColumnStatus(i) == atLowerBound) {
                objValue3 += columnLower[i] * value;
              } else if (sub[iBlock].getColumnStatus(i) == atUpperBound) {
                objValue3 += columnUpper[i] * value;
              }
            }
#endif
            // recompute
            if (logLevel > 3) {
              printf("objValue %g from lp %g, obj2 %g, obj3 %g\n",
                objValue, sub[iBlock].objectiveValue(),
                objValue2, objValue3);
            }
            objValue += objValue2;
            //objValue=sub[iBlock].objectiveValue();
            // See if cuts off and pack down
            int number = start;
            double infeas = -objValue;
            double smallest = 1.0e100;
            double largest = 0.0;
            for (int i = 0; i < numberMasterColumns; i++) {
              double value = elementAdd[start + i];
              if (fabs(value) > 1.0e-12) {
                infeas += primal[i] * value;
                smallest = std::min(smallest, fabs(value));
                largest = std::max(largest, fabs(value));
                indexColumnAdd[number] = i;
                elementAdd[number++] = -value;
              }
            }

            infeas += primal[numberMasterColumns + iBlock];
            indexColumnAdd[number] = numberMasterColumns + iBlock;
            elementAdd[number++] = -1.0;
            // if elements large then scale?
            if (largest > 1.0e8 || smallest < 1.0e-8) {
              sprintf(generalPrint, "For subproblem %d smallest - %g, largest %g - infeas %g",
                iBlock, smallest, largest, infeas);
              handler_->message(CLP_GENERAL2, messages_)
                << generalPrint
                << CoinMessageEol;
              if (smallest < 1.0e-12 * largest) {
                sprintf(generalPrint, "Removing small elements");
                handler_->message(CLP_GENERAL2, messages_)
                  << generalPrint
                  << CoinMessageEol;
                double target = 1.0e-12 * largest;
                smallest = largest;
                int number2 = number - 1;
                number = start;
                for (int i = start; i < number2; i++) {
                  double value = elementAdd[i];
                  if (fabs(value) > target) {
                    smallest = std::min(smallest, fabs(value));
                    indexColumnAdd[number] = indexColumnAdd[i];
                    elementAdd[number++] = value;
                  }
                }
                indexColumnAdd[number] = numberMasterColumns + iBlock;
                elementAdd[number++] = -1.0;
              }
            }
            // if smallest >1.0 then scale
            if ((smallest > 1.0e6 || fabs(objValue) > 0.01 * largeValue_)
              && number > start + 1) {
              double scale = 1.0 / smallest;
              if (fabs(scale * objValue) > 0.01 * largeValue_) {
                printf("** scale before obj scale %g\n", scale);
                scale = (0.01 * largeValue_) / fabs(objValue);
              }
              printf("** scale %g infeas %g\n", scale, infeas);
              objValue *= scale;
              infeas = -objValue;
              for (int i = start; i < number - 1; i++) {
                double value = elementAdd[i] * scale;
                elementAdd[i] = value;
                int iColumn = indexColumnAdd[i];
                infeas -= primal[iColumn] * value;
              }
              elementAdd[number - 1] *= scale;
              infeas += primal[numberMasterColumns + iBlock] * scale;
              printf("** new infeas %g - scales to %g\n", infeas, infeas / scale);
            }
            if (infeas < -1.0e-6 || (problemState[iBlock] & 4) != 0) {
              // take
              // double check infeasibility
              if (logLevel > 3)
                printf("objValue %g objectiveValue() %g\n",
                  objValue, sub[iBlock].objectiveValue());
              double sum = 0.0;
              for (int i = start; i < number; i++) {
                int iColumn = indexColumnAdd[i];
                sum += primal[iColumn] * elementAdd[i];
              }
              if (logLevel > 3)
                printf("Sum %g rhs %g\n", sum, -objValue);
              if (logLevel > 1)
                printf("Cut for block %d has %d elements\n", iBlock, number - 1 - start);
              blockPrint[numberProposals] = iBlock;
              objective[numberProposals] = -objValue;
              rowAdd[++numberProposals] = number;
              when[numberRowsGenerated] = iPass;
              whichBlock[numberRowsGenerated++] = iBlock;
            }
          } else if (sub[iBlock].isProvenPrimalInfeasible()) {
            // use ray
            double *solution = sub[iBlock].infeasibilityRay();
            if (0) {
              double trueOffset = 0.0;
              int numberRows = sub[iBlock].numberRows();
              int numberColumns = sub[iBlock].numberColumns();
              double *farkas = new double[std::max(2 * numberColumns + numberRows, numberMasterColumns)];
              double *bound = farkas + numberColumns;
              double *effectiveRhs = bound + numberColumns;
              // get ray as user would
              double *ray = solution; //sub[iBlock].infeasibilityRay();
              // get farkas row
              memset(farkas, 0, (2 * numberColumns + numberRows) * sizeof(double));
              // Looks to me as if ray should be flipped according to mosek
              sub[iBlock].clpMatrix()->transposeTimes(-1.0, ray, farkas);
              // now farkas has A_T_y
              // Put nonzero bounds in bound
              const double *columnLower = sub[iBlock].columnLower();
              const double *columnUpper = sub[iBlock].columnUpper();
              int numberBad = 0;
              // For sum in mosek
              double ySum = 0.0;
              for (int i = 0; i < numberColumns; i++) {
                double value = farkas[i];
                double boundValue = 0.0;
                if (sub[iBlock].getStatus(i) == ClpSimplex::basic) {
                  // treat as zero if small
                  if (fabs(value) < 1.0e-8) {
                    value = 0.0;
                    farkas[i] = 0.0;
                  }
                  if (value) {
                    //printf("basic %d direction %d farkas %g\n",
                    //	   i,sub[iBlock].directionOut(),value);
                    if (value < 0.0)
                      boundValue = columnLower[i];
                    else
                      boundValue = columnUpper[i];
                  }
                } else if (fabs(value) > 1.0e-8) {
                  if (value < 0.0)
                    boundValue = columnLower[i];
                  else
                    boundValue = columnUpper[i];
                  if (fabs(boundValue) > 1.0e12 && fabs(value) < 1.0e-8) {
                    boundValue = 0.0;
                    value = 0.0;
                  }
                } else {
                  value = 0.0;
                }
                if (fabs(boundValue) > 1.0e20) {
                  numberBad++;
                  boundValue = 0.0;
                  value = 0.0;
                  farkas[i] = 0.0;
                }
                // mosek way
                // A_T_y + s_x_l -s_x_u == 0
                // So if value >0 s_x_l->0 s_x_u->value
                // otherwise s_x_l->-value, s_x_u->0
                double s_x_l = 0.0;
                double s_x_u = 0.0;
                if (value > 0)
                  s_x_u = value;
                else
                  s_x_l = -value;
                ySum += columnLower[i] * s_x_l;
                ySum -= columnUpper[i] * s_x_u;
                bound[i] = boundValue;
              }
              const double *rowLower = sub[iBlock].rowLower();
              const double *rowUpper = sub[iBlock].rowUpper();
              //int pivotRow = sub[iBlock].spareIntArray_[3];
              //bool badPivot=pivotRow<0;
              for (int i = 0; i < numberRows; i++) {
                double value = ray[i];
                double rhsValue = 0.0;
                if (sub[iBlock].getRowStatus(i) == ClpSimplex::basic) {
                  // treat as zero if small
                  if (fabs(value) < 1.0e-7) {
                    value = 0.0;
                  }
                  if (value) {
                    //printf("row basic %d direction %d ray %g\n",
                    //	   i,sub[iBlock].directionOut(),value);
                    if (value < 0.0)
                      rhsValue = rowLower[i];
                    else
                      rhsValue = rowUpper[i];
                  }
                } else if (fabs(value) > 1.0e-10) {
                  if (value < 0.0)
                    rhsValue = rowLower[i];
                  else
                    rhsValue = rowUpper[i];
                } else {
                  value = 0.0;
                }
                if (fabs(rhsValue) > 1.0e20) {
                  numberBad++;
                  value = 0.0;
                }
                ray[i] = value;
                if (!value)
                  rhsValue = 0.0;
                // for mosek flip value back
                double yvalue = -value;
                // -y + s_c_l - s_c_u==0
                double s_c_l = 0.0;
                double s_c_u = 0.0;
                if (yvalue > 0)
                  s_c_l = yvalue;
                else
                  s_c_u = -yvalue;
                ySum += rowLower[i] * s_c_l;
                ySum -= rowUpper[i] * s_c_u;
                effectiveRhs[i] = rhsValue;
                if (fabs(effectiveRhs[i]) > 1.0e10)
                  printf("Large rhs row %d %g\n",
                    i, effectiveRhs[i]);
              }
              sub[iBlock].clpMatrix()->times(-1.0, bound, effectiveRhs);
              double bSum = 0.0;
              for (int i = 0; i < numberRows; i++) {
                bSum += effectiveRhs[i] * ray[i];
                if (fabs(effectiveRhs[i]) > 1.0e10)
                  printf("Large rhs row %d %g after\n",
                    i, effectiveRhs[i]);
              }
              if (logLevel > 1)
                printf("Block %d Mosek user manual wants %g to be positive so bSum should be negative %g\n",
                  iBlock, ySum, bSum);
              if (numberBad || bSum > 1.0e-6) {
                printf("Bad infeasibility ray %g  - %d bad\n",
                  bSum, numberBad);
              } else {
                //printf("Good ray - infeasibility %g\n",
                //     -bSum);
              }
              /*
		      wanted cut is
		      plus or minus! (temp2 * x - temp2 *x_bar) <= bSum
		      first[iBlock]->transposeTimes(ray, temp2);
		     */
              memset(farkas, 0, numberColumns * sizeof(double));
              first[iBlock]->transposeTimes(ray, farkas);
              double offset = 0.0;
              const double *masterSolution = masterModel.primalColumnSolution();
              for (int i = 0; i < numberMasterColumns; i++) {
                double value = farkas[i];
                if (fabs(value) > 1.0e-9) {
                  offset += value * masterSolution[i];
                  if (logLevel > 2)
                    printf("(%d,%g) ", i, value);
                } else {
                  farkas[i] = 0.0;
                }
              }
              trueOffset = bSum + offset;
              if (sub[iBlock].algorithm() > 0)
                trueOffset *= 1.0e-5;
              if (logLevel > 2)
                printf(" - offset %g - ? rhs of %g\n", offset, trueOffset);
              //delete [] ray;
              delete[] farkas;
            }
            // if primal then scale
            if (sub[iBlock].algorithm() > 0) {
              for (int i = 0; i < numberRows2; i++)
                solution[i] = -1.0e-5 * solution[i];
            } else {
              for (int i = 0; i < numberRows2; i++)
                solution[i] = -solution[i];
            }
            first[iBlock]->transposeTimes(solution, elementAdd + start);
            for (int i = 0; i < numberRows2; i++)
              solution[i] = -solution[i];
            for (int i = 0; i < numberRows2; i++) {
              if (sub[iBlock].getRowStatus(i) == basic && fabs(solution[i]) < 1.0e-7)
                solution[i] = 0.0;
              if (solution[i] > dualTolerance_) {
                // at upper
                if (saveUpper[i] > 1.0e20)
                  solution[i] = 0.0;
                objValue += solution[i] * saveUpper[i];
              } else if (solution[i] < -dualTolerance_) {
                // at lower
                if (saveLower[i] < -1.0e20)
                  solution[i] = 0.0;
                objValue += solution[i] * saveLower[i];
              } else {
                solution[i] = 0.0;
              }
            }
            //objValue=-objValue;
            {
              int numberColumns2 = sub[iBlock].numberColumns();
              double *temp = new double[numberColumns2];
              memset(temp, 0, numberColumns2 * sizeof(double));
              sub[iBlock].clpMatrix()->transposeTimes(-1.0, solution, temp);
              double loX = 0.0;
              double upX = 0.0;
              const double *lower = sub[iBlock].columnLower();
              const double *upper = sub[iBlock].columnUpper();
              const double *primal = sub[iBlock].primalColumnSolution();
              for (int i = 0; i < numberColumns2; i++) {
                double value = temp[i];
                if (sub[iBlock].getColumnStatus(i) == basic && fabs(value) < 1.0e-7)
                  value = 0.0;
                if (logLevel > 2) {
                  //if (sub[iBlock].getColumnStatus(i)!=basic&&
                  //  primal[i]>1.0e-5)
                  //printf("zz_inf %d has value %g\n",i,primal[i]);
                }
                if (sub[iBlock].getStatus(i) == atLowerBound || sub[iBlock].getStatus(i) == isFixed) {
                  loX += lower[i] * value;
                } else if (sub[iBlock].getStatus(i) == atUpperBound) {
                  upX += upper[i] * value;
                } else if (sub[iBlock].getStatus(i) == basic) {
                  double value2 = primal[i] * value;
                  if (logLevel > 2) {
                    if (fabs(value2) > 1.0e-3)
                      printf("Basic %d arrayval %g primal %g bounds %g %g\n",
                        i, value, primal[i],
                        lower[i], upper[i]);
                  }
                  if (value < 0.0) {
                    assert(primal[i] < lower[i]);
                    value2 = value * lower[i];
                  } else if (value > 0.0) {
                    assert(primal[i] > upper[i]);
                    if (primal[i] - upper[i] < 1.0e-3) {
                      if (logLevel > 2)
                        printf("small diff %g\n", primal[i] - upper[i]);
                      //handler_->message(CLP_GENERAL2, messages_)
                      //<< generalPrint
                      //<< CoinMessageEol;
                      //value=0.0;
                      //elementAdd[start+i]=0.0;
                    }
                    value2 = value * upper[i];
                  }
                  loX += value2;
                }
              }
              objValue += loX + upX;
              if (logLevel > 2)
                printf("Inf Offsets %g %g - new Objvalue %g\n", loX, upX, objValue);
#define OBJ_OFFSET 0
#if OBJ_OFFSET == 1
              objValue -= loX + upX;
              objValue -= loX + upX;
#elif OBJ_OFFSET == 2
              objValue -= loX + upX;
              objValue = -objValue;
              objValue += loX + upX;
#elif OBJ_OFFSET == 3
              objValue -= loX + upX;
              objValue = -objValue;
              objValue -= loX + upX;
#endif
              if (iBlock == -3) {
                ClpSimplex *temp = deBound(sub + iBlock);
                //temp->allSlackBasis();
                temp->primal();
                // use ray
                double *solution = temp->infeasibilityRay();
                int numberRows2 = temp->numberRows();
                // bug somewhere - if primal then flip
                if (temp->algorithm_ > 0) {
                  for (int i = 0; i < numberRows2; i++)
                    solution[i] = -1.0e-5 * solution[i];
                }
                double objValue7 = 0.0;
                const double *lower = temp->rowLower();
                const double *upper = temp->rowUpper();
                for (int i = 0; i < numberRows2; i++) {
                  if (solution[i] < -dualTolerance_) {
                    // at upper
                    assert(upper[i] < 1.0e30);
                    if (i < sub[iBlock].numberRows())
                      objValue7 += solution[i] * saveUpper[i];
                    else
                      objValue7 += solution[i] * upper[i];
                  } else if (solution[i] > dualTolerance_) {
                    // at lower
                    assert(lower[i] > -1.0e30);
                    if (i < sub[iBlock].numberRows())
                      objValue7 += solution[i] * saveLower[i];
                    else
                      objValue7 += solution[i] * lower[i];
                  }
                }
                sprintf(generalPrint, "new objValue %g - old %g", objValue7, objValue);
                handler_->message(CLP_GENERAL2, messages_)
                  << generalPrint
                  << CoinMessageEol;
                //objValue=-objValue7;
                //loX=1.0e-2*objValue;
                //first[iBlock]->transposeTimes(solution, elementAdd + start);
                double *temp2 = new double[numberMasterColumns];
                memset(temp2, 0, numberMasterColumns * sizeof(double));
                first[iBlock]->transposeTimes(solution, temp2);
                double *temp3 = elementAdd + start;
                for (int i = 0; i < numberMasterColumns; i++) {
                  if (fabs(temp2[i] - temp3[i]) > 1.0e-4)
                    printf("** %d bound el %g nobound %g\n", i, temp3[i], temp2[i]);
                }
                memcpy(temp3, temp2, numberMasterColumns * sizeof(double));
                delete[] temp2;
                delete temp;
              }
              // relax slightly
              objValue += 1.0e-9 * fabs(loX);
              delete[] temp;
            }
            delete[] solution;
            // See if good infeas and pack down (signs on infeas,value changed)
            int number = start;
            //printf("Not changing objValue from %g to %g\n",objValue,trueOffset);
            //printf("Changing objValue from %g to %g\n",objValue,trueOffset);
            //objValue=trueOffset;
            double infeas = objValue;
            double smallest = 1.0e100;
            double largest = 0.0;
            for (int i = 0; i < numberMasterColumns; i++) {
              double value = -elementAdd[start + i];
              if (fabs(value) > 1.0e-12) {
                infeas -= primal[i] * value;
                smallest = std::min(smallest, fabs(value));
                largest = std::max(largest, fabs(value));
                indexColumnAdd[number] = i;
                elementAdd[number++] = value;
              }
            }
            // if elements large or small then scale?
            if (largest > 1.0e8 || smallest < 1.0e-8) {
              sprintf(generalPrint, "For subproblem ray %d smallest - %g, largest %g - infeas %g",
                iBlock, smallest, largest, infeas);
              handler_->message(CLP_GENERAL2, messages_)
                << generalPrint
                << CoinMessageEol;
            }
            if (smallest < 1.0e-12 * largest) {
              sprintf(generalPrint, "Removing small elements");
              handler_->message(CLP_GENERAL2, messages_)
                << generalPrint
                << CoinMessageEol;
              smallest = 1.0e-12 * largest;
              int number2 = number - 1;
              number = start;
              for (int i = start; i < number2; i++) {
                double value = elementAdd[i];
                if (fabs(value) > smallest) {
                  indexColumnAdd[number] = indexColumnAdd[i];
                  elementAdd[number++] = value;
                }
              }
              indexColumnAdd[number] = numberMasterColumns + iBlock;
              elementAdd[number++] = -1.0;
            }
            if (infeas < -1.0e-6 || false) {
              double sum = 0.0;
              for (int i = start; i < number; i++) {
                int iColumn = indexColumnAdd[i];
                sum += primal[iColumn] * elementAdd[i];
              }
              if (logLevel > 2)
                printf("Sum %g rhs %g\n", sum, objValue);
              if (logLevel > 1)
                printf("Cut for block %d has %d elements (infeasibility)\n", iBlock, number - start);
              blockPrint[numberProposals] = iBlock;
              // take
              objective[numberProposals] = objValue;
              rowAdd[++numberProposals] = number;
              when[numberRowsGenerated] = iPass;
              whichBlock[numberRowsGenerated++] = iBlock;
            }
          } else {
            abort();
          }
        }
      } else {
        //printf("Can skip\n");
      }
      problemState[iBlock] &= ~4;
    }
    if (deletePrimal)
      delete[] primal;
    if (numberProposals) {
      sprintf(generalPrint, "%d cuts added with %d elements",
        numberProposals, rowAdd[numberProposals]);
      handler_->message(CLP_GENERAL, messages_)
        << generalPrint
        << CoinMessageEol;
      if (logLevel > 2) {
        for (int i = 0; i < numberProposals; i++) {
          printf("Cut %d block %d thetac %d ", i, blockPrint[i],
            blockPrint[i] + numberMasterColumns);
          int k = 0;
          for (CoinBigIndex j = rowAdd[i]; j < rowAdd[i + 1]; j++) {
            if (k == 12) {
              printf("\n");
              k = 0;
            }
            k++;
            printf("(%d,%g) ", indexColumnAdd[j], elementAdd[j]);
          }
          printf(" <= %g\n", objective[i]);
          //if (k)
          //printf("\n");
        }
      }
#ifdef TEST_MODEL
      if (logLevel > 3) {
        const double *solution = goodModel.primalColumnSolution();
        const double *solution2 = masterModel.primalColumnSolution();
        for (int i = 0; i < numberProposals; i++) {
          double sum = 0.0;
          double sum2 = 0.0;
          for (int j = rowAdd[i]; j < rowAdd[i + 1]; j++) {
            int iColumn = indexColumnAdd[j];
            double value = elementAdd[j];
            sum += value * solution[iColumn];
            sum2 += value * solution2[iColumn];
          }
          if (sum2 < objective[i] - 1.0e-4) {
            sprintf(generalPrint, "Rhs for cut %d (from block %d) does not cutoff sum2 %g sum %g rhs %g)",
              i, blockPrint[i], sum2, sum, objective[i]);
            handler_->message(CLP_GENERAL2, messages_)
              << generalPrint
              << CoinMessageEol;
          }
#define FIXUP_RHS 0
          if (sum > objective[i] + 1.0e-4) {
            sprintf(generalPrint, "Rhs for cut %d (from block %d) is %g too low (rhs is %g)",
              i, blockPrint[i], sum - objective[i], objective[i]);
            handler_->message(CLP_GENERAL2, messages_)
              << generalPrint
              << CoinMessageEol;
#if FIXUP_RHS == 1 || FIXUP_RHS == 3
            objective[i] = sum;
#endif
          } else if (sum < objective[i] - 1.0e-4) {
            sprintf(generalPrint, "Rhs for cut %d (from block %d) is %g ineffective (rhs is %g)",
              i, blockPrint[i], objective[i] - sum, objective[i]);
            handler_->message(CLP_GENERAL2, messages_)
              << generalPrint
              << CoinMessageEol;
#if FIXUP_RHS == 2 || FIXUP_RHS == 3
            objective[i] = sum;
#endif
          }
        }
      }
      if (logLevel > 3) {
        goodModel.addRows(numberProposals, NULL, objective,
          rowAdd, indexColumnAdd, elementAdd);
        goodModel.dual();
        if (goodModel.problemStatus() == 1 || goodModel.objectiveValue() > goodValue + 1.0e-5 * fabs(goodValue)) {
          int numberRows = goodModel.numberRows();
          int numberStart = numberRows - numberProposals;
          double *upper = goodModel.rowUpper();
          for (int iRow = numberStart; iRow < numberRows; iRow++) {
            upper[iRow] = COIN_DBL_MAX;
          }
          for (int iRow = numberStart; iRow < numberRows; iRow++) {
            upper[iRow] = objective[iRow - numberStart];
            goodModel.allSlackBasis(true);
            goodModel.dual();
            if (goodModel.problemStatus() == 1) {
              sprintf(generalPrint, "Cut %d makes infeasible - upper=%g",
                iRow - numberStart, upper[iRow]);
              handler_->message(CLP_GENERAL, messages_)
                << generalPrint
                << CoinMessageEol;
            } else if (goodModel.objectiveValue() > goodValue + 1.0e-5 * fabs(goodValue)) {
              sprintf(generalPrint, "Cut %d makes too expensive - upper=%g",
                iRow - numberStart, upper[iRow]);
              handler_->message(CLP_GENERAL, messages_)
                << generalPrint
                << CoinMessageEol;
              int iBlock = blockPrint[iRow - numberStart];
              ClpSimplex *temp = deBound(sub + iBlock);
              temp->allSlackBasis();
              temp->primal();
              // use ray
              double *solution = temp->infeasibilityRay();
              int numberRows2 = temp->numberRows();
              // bug somewhere - if primal then flip
              if (true) {
                for (int i = 0; i < numberRows2; i++)
                  solution[i] = -1.0e-5 * solution[i];
              }
              double objValue = 0.0;
              const double *lower = temp->rowLower();
              const double *upper = temp->rowUpper();
              for (int i = 0; i < numberRows2; i++) {
                if (solution[i] < -dualTolerance_) {
                  // at upper
                  assert(upper[i] < 1.0e30);
                  objValue += solution[i] * upper[i];
                } else if (solution[i] > dualTolerance_) {
                  // at lower
                  assert(lower[i] > -1.0e30);
                  objValue += solution[i] * lower[i];
                }
              }
              //printf("new objValue %g\n",objValue);
              int numberColumns2 = sub[iBlock].numberColumns();
              double *temp2 = new double[numberColumns2];
              memset(temp2, 0, numberColumns2 * sizeof(double));
              sub[iBlock].clpMatrix()->transposeTimes(1.0, sub[iBlock].infeasibilityRay(), temp2);
              double loX = 0.0;
              double upX = 0.0;
              const double *lower2 = sub[iBlock].columnLower();
              const double *upper2 = sub[iBlock].columnUpper();
              for (int i = 0; i < numberColumns2; i++) {
                if (sub[iBlock].getColumnStatus(i) != basic && fabs(sub[iBlock].primalColumnSolution()[i]) > 1.0e-5)
                  printf("zz_inf %d has value %g\n", i, sub[iBlock].primalColumnSolution()[i]);
                if (sub[iBlock].getStatus(i) == atLowerBound) {
                  loX += lower2[i] * temp2[i];
                } else if (sub[iBlock].getStatus(i) == atUpperBound) {
                  upX += upper2[i] * temp2[i];
                }
              }
              printf("Offsets %g %g\n", loX, upX);
              memset(temp2, 0, numberColumns2 * sizeof(double));
              temp->clpMatrix()->transposeTimes(1.0, solution, temp2);
              loX = 0.0;
              upX = 0.0;
              lower2 = temp->columnLower();
              upper2 = temp->columnUpper();
              for (int i = 0; i < numberColumns2; i++) {
                if (temp->getColumnStatus(i) != basic && fabs(temp->primalColumnSolution()[i]) > 1.0e-5)
                  printf("zz_inf %d has value %g\n", i, temp->primalColumnSolution()[i]);
                if (temp->getStatus(i) == atLowerBound) {
                  loX += lower2[i] * temp2[i];
                } else if (temp->getStatus(i) == atUpperBound) {
                  upX += upper2[i] * temp2[i];
                }
              }
              printf("Offsets %g %g\n", loX, upX);
              delete[] temp2;
              delete temp;
            }
            upper[iRow] = COIN_DBL_MAX;
          }
        }
        double objValue = goodModel.objectiveValue();
        const double *obj = goodModel.objective();
        const double *solution = goodModel.primalColumnSolution();
        double obj1 = 0.0;
        for (int i = 0; i < numberMasterColumns; i++)
          obj1 += obj[i] * solution[i];
        double obj2 = 0.0;
        for (int i = numberMasterColumns; i < numberMasterColumns + numberBlocks; i++)
          obj2 += obj[i] * solution[i];
        double obj3 = 0.0;
        for (int i = numberMasterColumns + numberBlocks; i < goodModel.numberColumns(); i++)
          obj3 += obj[i] * solution[i];
        //assert (fabs(goodValue-objValue)<1.0e-3+1.0e-7*fabs(goodValue));
        printf("XXXX good %g this %g difference %g - objs %g, %g, %g\n",
          goodValue, objValue, goodValue - objValue, obj1, obj2, obj3);
      }
#endif
      masterModel.addRows(numberProposals, NULL, objective,
        rowAdd, indexColumnAdd, elementAdd);
    }
  }
  sprintf(generalPrint, "Time at end of Benders %.2f seconds", CoinCpuTime() - time1);
  handler_->message(CLP_GENERAL, messages_)
    << generalPrint
    << CoinMessageEol;
  delete[] problemState;
  for (int iBlock = 0; iBlock < numberBlocks; iBlock++) {
    delete[] modification[iBlock];
  }
  delete[] modification;
#ifdef ADD_ARTIFICIALS
  delete[] originalSubColumns;
  for (int i = 0; i < numberBlocks; i++)
    delete[] saveObjective[i];
  delete[] saveObjective;
#endif
  //masterModel.scaling(0);
  //masterModel.primal(1);
  if (!options.independentOption(1))
    loadProblem(*model);
  // now put back a good solution
  const double *columnSolution = masterModel.primalColumnSolution();
  double *fullColumnSolution = primalColumnSolution();
  const int *columnBack = model->coinBlock(masterBlock)->originalColumns();
  int numberColumns2 = model->coinBlock(masterBlock)->numberColumns();
  const int *rowBack = model->coinBlock(masterBlock)->originalRows();
  int numberRows2 = model->coinBlock(masterBlock)->numberRows();
#ifndef NDEBUG
  for (int iColumn = 0; iColumn < numberColumns_; iColumn++)
    fullColumnSolution[iColumn] = COIN_DBL_MAX;
#endif
  for (int iColumn = 0; iColumn < numberColumns2; iColumn++) {
    int kColumn = columnBack[iColumn];
    setColumnStatus(kColumn, masterModel.getColumnStatus(iColumn));
    fullColumnSolution[kColumn] = columnSolution[iColumn];
  }
  for (int iRow = 0; iRow < numberRows2; iRow++) {
    int kRow = rowBack[iRow];
    setRowStatus(kRow, masterModel.getRowStatus(iRow));
    //fullSolution[kRow]=solution[iRow];
  }
  for (iBlock = 0; iBlock < numberBlocks; iBlock++) {
    // move basis
    int kBlock = rowCounts[iBlock];
    const int *columnBack = model->coinBlock(kBlock)->originalColumns();
    int numberColumns2 = model->coinBlock(kBlock)->numberColumns();
    const int *rowBack = model->coinBlock(kBlock)->originalRows();
    int numberRows2 = model->coinBlock(kBlock)->numberRows();
    const double *subColumnSolution = sub[iBlock].primalColumnSolution();
    for (int iColumn = 0; iColumn < numberColumns2; iColumn++) {
      int kColumn = columnBack[iColumn];
      setColumnStatus(kColumn, sub[iBlock].getColumnStatus(iColumn));
#ifndef NDEBUG
      assert(fullColumnSolution[kColumn] == COIN_DBL_MAX);
#endif
      fullColumnSolution[kColumn] = subColumnSolution[iColumn];
    }
    for (int iRow = 0; iRow < numberRows2; iRow++) {
      int kRow = rowBack[iRow];
      setRowStatus(kRow, sub[iBlock].getRowStatus(iRow));
      //setStatus(kRow, atLowerBound);
    }
  }
#ifndef NDEBUG
  for (int iColumn = 0; iColumn < numberColumns_; iColumn++)
    assert(fullColumnSolution[iColumn] != COIN_DBL_MAX);
#endif
  double *fullSolution = primalRowSolution();
  CoinZeroN(fullSolution, numberRows_);
  times(1.0, fullColumnSolution, fullSolution);
  int numberRowBasic = 0;
#ifndef NDEBUG
  int numberInfeasibilities = 0;
  double sumInfeasibilities = 0.0;
#endif
  for (int iRow = 0; iRow < numberRows_; iRow++) {
    if (getRowStatus(iRow) == ClpSimplex::basic)
      numberRowBasic++;
#ifndef NDEBUG
    if (fullSolution[iRow] < rowLower_[iRow] - primalTolerance_) {
      numberInfeasibilities++;
      sumInfeasibilities -= fullSolution[iRow] - rowLower_[iRow];
      if (getRowStatus(iRow) != basic)
        setRowStatus(iRow, superBasic);
    } else if (fullSolution[iRow] > rowUpper_[iRow] + primalTolerance_) {
      numberInfeasibilities++;
      sumInfeasibilities += fullSolution[iRow] - rowUpper_[iRow];
      if (getRowStatus(iRow) != basic)
        setRowStatus(iRow, superBasic);
    }
#endif
  }
  int numberColumnBasic = 0;
  for (int iColumn = 0; iColumn < numberColumns_; iColumn++)
    if (getColumnStatus(iColumn) == ClpSimplex::basic)
      numberColumnBasic++;
  sprintf(generalPrint, "%d row basic %d col basic (total %d) - wanted %d",
    numberRowBasic, numberColumnBasic,
    numberRowBasic + numberColumnBasic, numberRows_);
  handler_->message(CLP_GENERAL2, messages_)
    << generalPrint
    << CoinMessageEol;
#ifndef NDEBUG
  sprintf(generalPrint, "%d infeasibilities summing to %g",
    numberInfeasibilities, sumInfeasibilities);
  handler_->message(CLP_GENERAL2, messages_)
    << generalPrint
    << CoinMessageEol;
#endif
  //for (int i=0;i<numberRows_;i++) setRowStatus(i,basic);
  sprintf(generalPrint, "Time before cleanup of full model %.2f seconds", CoinCpuTime() - time1);
  handler_->message(CLP_GENERAL, messages_)
    << generalPrint
    << CoinMessageEol;
  this->primal(1);
  sprintf(generalPrint, "Total time %.2f seconds", CoinCpuTime() - time1);
  handler_->message(CLP_GENERAL, messages_)
    << generalPrint
    << CoinMessageEol;
  delete[] rowCounts;
  //delete [] sol;
  //delete [] lower;
  //delete [] upper;
  delete[] whichBlock;
  delete[] when;
  delete[] rowAdd;
  delete[] columnAdd;
  delete[] elementAdd;
  delete[] objective;
  delete[] first;
  delete[] sub;
  return 0;
}
