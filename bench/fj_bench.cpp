/* fj_bench.cpp - standalone Feasibility Jump benchmark
 *
 * Loads an MPS instance, runs the full root-LP preprocessing pipeline
 * (bound propagation, clique merging, LP recommender) via CbcSolver, then
 * runs the FJ heuristic directly (bypassing CbcModel) with a custom callback
 * that captures per-solution timing and effort metrics.
 *
 * Usage:
 *   fj_bench [options] <instance.mps[.gz]>
 *
 * Options:
 *   --effort N          Iteration budget (default: 10000000)
 *   --maxsol N          Stop after N solutions (default: 5)
 *   --seed N            PRNG seed (default: 0)
 *   --timelimit S       Wall-clock limit in seconds (default: 300)
 *   --relax-continuous  Relax continuous variable constraints (default: off)
 *   --loglevel N        0=silent 1=FJ progress 2=verbose (default: 1)
 *   --header            Print TSV header line then exit
 *
 * Output: one TSV line per run to stdout:
 *   instance  status  obj  lp_time  fj_time  time_first  effort_first
 *   total_effort  solutions  effort_budget  timelimit  seed
 *
 * Designed for parallel invocation with GNU parallel, one instance per job.
 */

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "CbcSolver.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wsign-compare"
#endif
#include "feasibilityjump.hh"
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

static double now_secs()
{
  using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();
}

static std::string basename_noext(const std::string &path)
{
  size_t slash = path.rfind('/');
  std::string name = (slash == std::string::npos) ? path : path.substr(slash + 1);
  // strip .gz
  if (name.size() > 3 && name.substr(name.size() - 3) == ".gz")
    name = name.substr(0, name.size() - 3);
  // strip .mps
  if (name.size() > 4 && name.substr(name.size() - 4) == ".mps")
    name = name.substr(0, name.size() - 4);
  return name;
}

// ---------------------------------------------------------------------------
// result struct
// ---------------------------------------------------------------------------

struct FJResult {
  bool found = false;
  double bestObj = DBL_MAX;
  double timeToFirstSol = 0.0;
  double totalTime = 0.0;
  int64_t effortAtFirstSol = 0;
  int64_t totalEffort = 0;
  int solutionsFound = 0;
};

// ---------------------------------------------------------------------------
// run FJ directly from a preprocessed, LP-solved Osi model
// ---------------------------------------------------------------------------

static FJResult runFJ(
  OsiSolverInterface &osi,
  int64_t maxEffort,
  int maxSol,
  int seed,
  double timelimit,
  bool relaxContinuous,
  int loglevel)
{
  const int numCols = osi.getNumCols();
  const int numRows = osi.getNumRows();

  FeasibilityJumpSolver fj(seed, loglevel >= 2 ? 1 : 0, 1.0);

  // add variables
  const double *colLower = osi.getColLower();
  const double *colUpper = osi.getColUpper();
  const double *objCoeff = osi.getObjCoefficients();
  for (int j = 0; j < numCols; ++j) {
    VarType vt = osi.isInteger(j) ? VarType::Integer : VarType::Continuous;
    fj.addVar(vt, colLower[j], colUpper[j], objCoeff[j]);
  }

  // add constraints
  const CoinPackedMatrix *mat = osi.getMatrixByRow();
  const double *elements = mat->getElements();
  const int *colIdxs = mat->getIndices();
  const CoinBigIndex *rowStarts = mat->getVectorStarts();
  const int *rowLengths = mat->getVectorLengths();
  const char *rowSense = osi.getRowSense();
  const double *rhs = osi.getRightHandSide();
  const double *rowRange = osi.getRowRange();
  int relax = relaxContinuous ? 1 : 0;

  for (int i = 0; i < numRows; ++i) {
    int len = rowLengths[i];
    std::vector<int> idxBuf(colIdxs + rowStarts[i], colIdxs + rowStarts[i] + len);
    std::vector<double> elBuf(elements + rowStarts[i], elements + rowStarts[i] + len);
    char sense = rowSense[i];
    if (sense == 'E') {
      fj.addConstraint(RowType::Equal, rhs[i], len, idxBuf.data(), elBuf.data(), relax);
    } else if (sense == 'L') {
      fj.addConstraint(RowType::Lte, rhs[i], len, idxBuf.data(), elBuf.data(), relax);
    } else if (sense == 'G') {
      fj.addConstraint(RowType::Gte, rhs[i], len, idxBuf.data(), elBuf.data(), relax);
    } else if (sense == 'R') {
      fj.addConstraint(RowType::Lte, rhs[i], len, idxBuf.data(), elBuf.data(), relax);
      fj.addConstraint(RowType::Gte, rhs[i] - rowRange[i], len, idxBuf.data(), elBuf.data(), relax);
    }
    // 'N' free row -- skip
  }

  // initial point: LP solution rounded for integers
  const double *lpSol = osi.getColSolution();
  std::vector<double> init(numCols);
  for (int j = 0; j < numCols; ++j) {
    double v = osi.isInteger(j) ? std::floor(lpSol[j] + 0.5) : lpSol[j];
    init[j] = std::max(colLower[j], std::min(colUpper[j], v));
  }

  FJResult result;
  double sense = osi.getObjSense(); // +1 min, -1 max
  double startTime = now_secs();
  double deadline = timelimit > 0.0 ? startTime + timelimit : DBL_MAX;

  auto callback = [&](FJStatus status) -> CallbackControlFlow {
    double t = now_secs();

    if (t >= deadline)
      return CallbackControlFlow::Terminate;

    if (status.totalEffort >= maxEffort)
      return CallbackControlFlow::Terminate;

    result.totalEffort = status.totalEffort;

    if (status.solution != nullptr) {
      double obj = status.solutionObjectiveValue * sense;
      if (!result.found || obj < result.bestObj) {
        result.bestObj = obj;
      }
      if (!result.found) {
        result.timeToFirstSol = t - startTime;
        result.effortAtFirstSol = status.totalEffort;
        result.found = true;
      }
      result.solutionsFound++;

      if (loglevel >= 1) {
        fprintf(stdout, "  sol %d  obj=%.6g  effort=%lldM  t=%.2fs\n",
          result.solutionsFound, obj,
          (long long)(status.totalEffort / 1000000),
          t - startTime);
        fflush(stdout);
      }

      if (result.solutionsFound >= maxSol)
        return CallbackControlFlow::Terminate;
    }

    return CallbackControlFlow::Continue;
  };

  fj.solve(init.data(), callback);

  result.totalTime = now_secs() - startTime;
  return result;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

static void print_header()
{
  printf("instance\tstatus\tobj\tlp_time\tfj_time\t"
         "time_first\teffort_first\ttotal_effort\t"
         "solutions\teffort_budget\ttimelimit\tseed\n");
}

int main(int argc, char *argv[])
{
  // defaults
  int64_t maxEffort = 10000000;
  int maxSol = 5;
  int seed = 0;
  double timelimit = 300.0;
  bool relaxContinuous = false;
  int loglevel = 1;
  const char *mpsFile = nullptr;

  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--header") == 0) {
      print_header();
      return 0;
    } else if (strcmp(argv[i], "--effort") == 0 && i + 1 < argc) {
      maxEffort = (int64_t)atoll(argv[++i]);
    } else if (strcmp(argv[i], "--maxsol") == 0 && i + 1 < argc) {
      maxSol = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
      seed = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--timelimit") == 0 && i + 1 < argc) {
      timelimit = atof(argv[++i]);
    } else if (strcmp(argv[i], "--relax-continuous") == 0) {
      relaxContinuous = true;
    } else if (strcmp(argv[i], "--loglevel") == 0 && i + 1 < argc) {
      loglevel = atoi(argv[++i]);
    } else if (argv[i][0] != '-') {
      mpsFile = argv[i];
    } else {
      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
    }
  }

  if (!mpsFile) {
    fprintf(stderr,
      "Usage: fj_bench [options] <instance.mps[.gz]>\n"
      "  --effort N          iteration budget (default: 10000000)\n"
      "  --maxsol N          stop after N solutions (default: 5)\n"
      "  --seed N            PRNG seed (default: 0)\n"
      "  --timelimit S       wall-clock limit in seconds (default: 300)\n"
      "  --relax-continuous  relax continuous constraints\n"
      "  --loglevel N        0=silent 1=solutions 2=verbose (default: 1)\n"
      "  --header            print TSV header and exit\n");
    return 1;
  }

  std::string instName = basename_noext(mpsFile);

  // Load MPS file with a silenced OsiClp (controlled output)
  OsiClpSolverInterface osi;
  osi.messageHandler()->setLogLevel(0);
  // readMps(filename, keepNames, allowErrors) handles .mps.gz directly
  if (osi.readMps(mpsFile, false, false) != 0) {
    fprintf(stderr, "ERROR: could not read %s\n", mpsFile);
    return 1;
  }

  // check for integer variables before committing to solve
  const int numCols = osi.getNumCols();
  bool hasIntegers = false;
  for (int j = 0; j < numCols; ++j)
    if (osi.isInteger(j)) { hasIntegers = true; break; }
  if (!hasIntegers) {
    fprintf(stderr, "ERROR: no integer variables in %s\n", mpsFile);
    return 1;
  }

  // Run full root-LP pipeline via CbcSolver:
  //   CbcSolver(osi)  -- copies model, no file-import output
  //   initialize()    -- sets up parameter defaults
  //   solveRootLp()   -- bound propagation + clique merging + LP solve
  CbcSolver cbcSolver(osi);
  cbcSolver.initialize();

  if (loglevel == 0) {
    // parameters_.noPrinting() is what applyLpMethod checks to silence Clp
    cbcSolver.parameters().disablePrinting();
    // message handler covers bound propagation, clique merging, printGeneralMessage
    cbcSolver.model()->messageHandler()->setLogLevel(0);
  }

  double lpStart = now_secs();
  int lpRc = cbcSolver.solveRootLp(); // 0=optimal, 1=stopped, 2=infeasible, 3=unbounded
  double lpTime = now_secs() - lpStart;

  if (lpRc == 2 || lpRc == 3) {
    fprintf(stderr, "WARNING: root LP result %d for %s\n", lpRc, mpsFile);
    printf("%s\tINFEASIBLE\tNA\t%.4f\tNA\tNA\tNA\tNA\t0\t%lld\t%.0f\t%d\n",
      instName.c_str(), lpTime,
      (long long)maxEffort, timelimit, seed);
    return 2;
  }

  // LP solution lives in the (possibly preprocessed) model's solver.
  // Bounds may have been tightened; integer flags are preserved.
  OsiSolverInterface *preparedOsi = cbcSolver.model()->solver();

  if (loglevel >= 1) {
    fprintf(stdout, "Instance: %s  cols=%d rows=%d  lp_time=%.2fs\n",
      instName.c_str(), numCols, osi.getNumRows(), lpTime);
    fflush(stdout);
  }

  FJResult r = runFJ(*preparedOsi, maxEffort, maxSol, seed, timelimit,
                     relaxContinuous, loglevel);

  // TSV output line
  printf("%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%lld\t%lld\t%d\t%lld\t%.0f\t%d\n",
    instName.c_str(),
    r.found ? "FOUND" : "NOTFOUND",
    r.found ? [&]{ static char buf[32]; snprintf(buf, sizeof(buf), "%.6g", r.bestObj); return buf; }() : "NA",
    lpTime,
    r.totalTime,
    r.found ? r.timeToFirstSol : -1.0,
    r.found ? (long long)r.effortAtFirstSol : -1LL,
    (long long)r.totalEffort,
    r.solutionsFound,
    (long long)maxEffort,
    timelimit,
    seed);

  return r.found ? 0 : 2;
}
