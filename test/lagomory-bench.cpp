/*
 * lagomory-bench -- replay one Lagrangean Gomory (lagomory) separation call from
 * a fixture and emit one CSV row describing what it cost and what it produced.
 *
 * WHY THIS IS NOT gomory-bench --orig-solver. gomory-bench has a --orig-solver
 * flag, and it does call passInOriginalSolver, but it passes **the fixture solver
 * itself**:
 *
 *     if (origSolver) gomory.passInOriginalSolver(&f.si);
 *
 * so `originalSolver_` and `si` are the same LP and `numberOriginalRows ==
 * si.getNumRows()`. Every Lagrangean quantity is then trivial -- there are no cut
 * rows to dualize, the perturbation is empty, and the resolve lands back on the
 * same vertex. Measured over the 330 gomoryFixtures: `--gomory-type=11` and `=12`
 * yield 0 cuts in ~4 us (the whenToDo==1 gate is simply false), and `=21`/`=22`
 * yield statistics **byte-identical to plain `=0`**. That is not lagomory being
 * fast, it is lagomory not running.
 *
 * The cause is upstream of the bench: the Cbc-side fixture dump fires at
 * `currentPassNumber_ == 1`, before any generator has contributed a row, and
 * every gomoryFixture's `.meta` confirms it -- `rows` equals
 * `infoFormulationRows` in all 330. A lagomory fixture has to be captured where
 * the LP *already holds cut rows*, which is what CglLagomoryFixtureDump.hpp does
 * from inside CglGomory::generateCuts.
 *
 * THE TWO SOLVERS, AND WHY BOTH COME FROM DISK. A lagomory call is defined by a
 * pair: the augmented LP `si` (formulation + the cut rows to be dualized away,
 * plus the duals that weight them) and the original formulation
 * `originalSolver_` (the matrix the cuts will actually be rows of). Only the
 * first is a normal fixture. The second cannot be reconstructed by deleting
 * si's cut rows, because the Lagrangean pass copies si's *column* bounds into it
 * on every call and never its *row* bounds, so the two drift -- `rowBoundDrift`
 * in the `.meta` reports how far, and a nonzero value is proof the shortcut
 * would have replayed the wrong LP. So `.orig.mps.gz` is stored separately and
 * loaded separately here.
 *
 * WHAT IS TIMED. The CSV's sepTime is the generateCuts call and nothing else:
 * both solvers are loaded, warm-started and validated first, and the generator
 * is constructed outside the timed region. --repeat takes the minimum, because
 * the thing being measured is the code and the noise is only ever additive.
 * Timing is meaningful ONLY serially -- see the note in
 * BENCHMARKING-CUT-GENERATORS.md about never timing under a parallel driver, and
 * note that the stage accumulators behind --profile are file-static in
 * CglGomory.cpp and would be summed across threads.
 *
 * WHAT THE FIDELITY COLUMNS ARE FOR. piMaxDev and solMaxDev compare the replayed
 * LP's duals and primal solution against the `.pi`/`.sol` captured at dump time.
 * They are *checks*, never injected: the replay reaches its own duals by reading
 * the `.bas` and resolving, exactly as CBC would, and if those duals disagree
 * with the captured ones then the perturbed objective is different and the run
 * is separating a different problem. A fixture with a large piMaxDev is not noisy,
 * it is invalid -- filter on it rather than averaging it in.
 *
 * Usage:
 *   lagomory-bench --header
 *   lagomory-bench <fixture-stem> [options]
 *
 * Exit codes: 0 ok, 1 load/replay failure, 2 usage.
 */

#include "CglGomory.hpp"
#include "CglTreeInfo.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCuts.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <sys/stat.h>
#include <sys/time.h>
#include <vector>

static double wallClock()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + 1.0e-6 * (double)tv.tv_usec;
}

/*
 * Stage attribution inside the generator, off unless Cgl was built
 * -DCGL_GOMORY_PROFILE. Declared by hand rather than pulled from a Cgl header on
 * purpose: this is diagnostic scaffolding, not API, and shipping the prototypes
 * would turn a normal Cgl build into a link error for everyone. The cost is that
 * these two signatures must be kept in step with the definitions near the top of
 * CglGomory.cpp.
 */
#ifdef CGL_GOMORY_PROFILE
void cglGomoryProfileReset();
void cglGomoryProfilePrint(const char *tag);
#define LAG_PROF_RESET() cglGomoryProfileReset()
#define LAG_PROF_PRINT(tag) cglGomoryProfilePrint(tag)
#else
#define LAG_PROF_RESET()
#define LAG_PROF_PRINT(tag)
#endif

static bool fileExists(const std::string &path)
{
  struct stat st;
  return stat(path.c_str(), &st) == 0;
}

/// Reduce any of the fixture's file names to the shared stem, so a caller can
/// pass whichever one tab-completion produced. The `.orig.*` forms are listed
/// before the bare ones so that "x.orig.mps.gz" does not stem to "x.orig".
static std::string fixtureStem(const char *arg)
{
  std::string s(arg);
  static const char *suffixes[] = { ".orig.mps.gz", ".orig.mps", ".orig.ctype",
    ".mps.gz", ".mps", ".bas", ".sol", ".pi", ".ctype", ".meta", ".bas.status" };
  for (size_t i = 0; i < sizeof(suffixes) / sizeof(suffixes[0]); ++i) {
    const std::string suf(suffixes[i]);
    if (s.size() > suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0)
      return s.substr(0, s.size() - suf.size());
  }
  return s;
}

static std::string baseName(const std::string &path)
{
  const size_t slash = path.rfind('/');
  return slash == std::string::npos ? path : path.substr(slash + 1);
}

/// One numeric key out of the `.meta`; `dflt` when the file or key is absent, so
/// a fixture written before a key existed still loads.
static double metaNum(const std::string &path, const char *key, double dflt)
{
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp)
    return dflt;
  char line[512];
  double value = dflt;
  while (fgets(line, sizeof(line), fp)) {
    char k[256];
    double v = 0.0;
    if (sscanf(line, "%255s %lf", k, &v) == 2 && strcmp(k, key) == 0) {
      value = v;
      break;
    }
  }
  fclose(fp);
  return value;
}

static long metaInt(const std::string &path, const char *key, long dflt)
{
  return (long)metaNum(path, key, (double)dflt);
}

/**
 * Remove the pad row, if this fixture has one.
 *
 * A capture whose matrix held empty columns carries one extra final row, without
 * which writeMps drops those columns and shifts every index after them. The row
 * is redundant but not cosmetic: it is an extra basis position, so leaving it in
 * changes the factorization the cuts come from, and it leaves the `.bas` one
 * artificial short of the matrix.
 *
 * Keyed on the meta field names because this file loads two matrices from one
 * `.meta`: `paddedColumns`/`rows` describe si, `origPaddedColumns`/`origRows`
 * describe the original formulation, and they pad independently.
 *
 * Both conditions are required before deleting anything -- meta must say padding
 * happened *and* the loaded row count must be exactly one more than the captured
 * one. Either alone could delete a real row from a fixture whose meta is stale.
 */
static bool dropPadRow(OsiSolverInterface &si, const std::string &stem,
  const char *padKey, const char *rowsKey, const char *label, bool quiet)
{
  const std::string meta = stem + ".meta";
  if (metaInt(meta, padKey, 0) <= 0)
    return false;

  const long capturedRows = metaInt(meta, rowsKey, -1);
  if (capturedRows < 0 || si.getNumRows() != (int)capturedRows + 1) {
    if (!quiet)
      fprintf(stderr, "WARNING: %s (%s): meta says padded but rows=%d against "
                      "captured %ld; leaving the matrix alone\n",
        baseName(stem).c_str(), label, si.getNumRows(), capturedRows);
    return false;
  }
  const int last = si.getNumRows() - 1;
  si.deleteRows(1, &last);
  return true;
}

/**
 * Restore integrality from a `.ctype` sidecar.
 *
 * MPS conveys integrality only through the bound type, and a column with
 * lb == ub takes writeMps's " FX " branch, which has no integer form -- so every
 * integer column CBC had fixed by bound tightening reads back **continuous**.
 *
 * For Gomory a lost marker changes the experiment in three separate places:
 * `intVar[]` gates candidate selection; it selects the mixed-integer coefficient
 * formula rather than the continuous one, and a column taking the continuous
 * branch increments `numberNonInteger`, which moves the rhs relaxation ladder --
 * so the right-hand side changes even for cuts that still get generated; and it
 * gates the integer-slack test that decides whether a slack may contribute an
 * integral coefficient. Skipping this does not measure a smaller problem, it
 * measures different cuts.
 *
 * Applied to BOTH solvers here. The original formulation is where the cuts'
 * coefficients are actually computed, so a marker missing there is if anything
 * worse than one missing in si.
 *
 * A sidecar whose column count disagrees with the model is refused rather than
 * partly applied: marking arbitrary columns integer is worse than the loss it
 * repairs.
 */
static int restoreColTypes(OsiSolverInterface &si, const std::string &path,
  const std::string &label, bool quiet)
{
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp) {
    if (!quiet)
      fprintf(stderr, "WARNING: %s: no %s; integer columns fixed at capture will "
                      "read back continuous\n",
        label.c_str(), baseName(path).c_str());
    return -1;
  }
  int sidecarCols = -1;
  if (fscanf(fp, "cols %d\n", &sidecarCols) != 1 || sidecarCols != si.getNumCols()) {
    fprintf(stderr, "ERROR: %s: %s is for %d columns, model has %d; ignoring it\n",
      label.c_str(), baseName(path).c_str(), sidecarCols, si.getNumCols());
    fclose(fp);
    return -1;
  }
  int idx = 0, type = 0, restored = 0;
  while (fscanf(fp, "%d %d\n", &idx, &type) == 2) {
    if (idx < 0 || idx >= si.getNumCols()) {
      fprintf(stderr, "ERROR: %s: %s names column %d, out of range\n",
        label.c_str(), baseName(path).c_str(), idx);
      fclose(fp);
      return -1;
    }
    if (si.isContinuous(idx)) {
      si.setInteger(idx);
      ++restored;
    }
  }
  fclose(fp);
  // getColType() caches, and CglGomory reads it through si.getColType() to build
  // intVar[], so it must be recomputed or the generator sees the pre-restore view.
  si.getColType(true);
  return restored;
}

/**
 * Largest absolute deviation between the replayed row duals and the `.pi` written
 * at capture. Returns -1.0 when there is no `.pi` to compare against, which is
 * reported as such rather than silently as a pass.
 *
 * Only the CUT rows are compared, and that is the point rather than a shortcut:
 * pi over rows [formulationRows, rows) is precisely the vector that perturbs the
 * objective, so it is the only part whose disagreement changes the separation.
 * Duals on the formulation rows are never read by the Lagrangean pass.
 */
static double piDeviation(const OsiSolverInterface &si, const std::string &stem,
  int formulationRows, int &compared)
{
  compared = 0;
  const std::string path = stem + ".pi";
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp)
    return -1.0;

  int m = -1, formulation = -1;
  if (fscanf(fp, "pi %d formulation %d\n", &m, &formulation) != 2) {
    fclose(fp);
    return -1.0;
  }
  const double *pi = si.getRowPrice();
  if (!pi || m != si.getNumRows()) {
    fclose(fp);
    return -1.0;
  }
  double worst = 0.0;
  int idx = 0;
  double v = 0.0;
  while (fscanf(fp, "%d %lf\n", &idx, &v) == 2) {
    if (idx < formulationRows || idx >= m)
      continue;
    worst = std::max(worst, fabs(v - pi[idx]));
    ++compared;
  }
  fclose(fp);
  return worst;
}

/// Largest absolute deviation between the replayed primal solution and the
/// `.sol` written at capture. The `.sol` lists only nonzeros, so a captured zero
/// that came back nonzero would be missed -- which is why the pi comparison is
/// the primary check and this one is corroboration.
static double solDeviation(const OsiSolverInterface &si, const std::string &stem,
  int &compared)
{
  compared = 0;
  const std::string path = stem + ".sol";
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp)
    return -1.0;
  const double *x = si.getColSolution();
  if (!x) {
    fclose(fp);
    return -1.0;
  }
  char line[1024];
  double worst = 0.0;
  const int n = si.getNumCols();
  while (fgets(line, sizeof(line), fp)) {
    if (line[0] == '=')
      continue;
    int idx = 0;
    char name[256];
    double v = 0.0;
    if (sscanf(line, "%d %255s %lf", &idx, name, &v) != 3)
      continue;
    if (idx < 0 || idx >= n)
      continue;
    worst = std::max(worst, fabs(v - x[idx]));
    ++compared;
  }
  fclose(fp);
  return worst;
}

/*
 * Emit every cut, coefficient by coefficient, in a form where textual equality
 * is bit equality.
 *
 * The CSV fields are a fingerprint, not a proof: totalViol and avgCutLen are sums
 * over all cuts, so two different cut sets can agree on every one of them. When
 * the claim is "this optimization changes nothing", a fingerprint match is the
 * weaker statement. %a prints the exact IEEE double, so a one-ulp difference in a
 * single coefficient of a single cut shows up as a diff rather than rounding away
 * at the 15th digit. Order is as generated, deliberately unsorted: the order cuts
 * land in OsiCuts is itself part of what must not change.
 */
static void dumpCuts(const OsiCuts &cs, const std::string &name, int round)
{
  for (int c = 0; c < cs.sizeRowCuts(); ++c) {
    const OsiRowCut &rc = cs.rowCut(c);
    const CoinPackedVector &row = rc.row();
    const int len = row.getNumElements();
    const int *idx = row.getIndices();
    const double *el = row.getElements();
    printf("[cut] %s r%d row %d n %d lb %a ub %a gv %d", name.c_str(), round, c,
      len, rc.lb(), rc.ub(), rc.globallyValid() ? 1 : 0);
    for (int k = 0; k < len; ++k)
      printf(" %d:%a", idx[k], el[k]);
    printf("\n");
  }
  for (int c = 0; c < cs.sizeColCuts(); ++c) {
    const OsiColCut &cc = cs.colCut(c);
    const CoinPackedVector &lbs = cc.lbs();
    const CoinPackedVector &ubs = cc.ubs();
    printf("[cut] %s r%d col %d nlb %d nub %d", name.c_str(), round, c,
      lbs.getNumElements(), ubs.getNumElements());
    for (int k = 0; k < lbs.getNumElements(); ++k)
      printf(" L%d:%a", lbs.getIndices()[k], lbs.getElements()[k]);
    for (int k = 0; k < ubs.getNumElements(); ++k)
      printf(" U%d:%a", ubs.getIndices()[k], ubs.getElements()[k]);
    printf("\n");
  }
}

struct CutStats {
  int rowCuts = 0;
  int colCuts = 0;
  int globallyValid = 0;
  double totalViol = 0.0;
  double maxViol = 0.0;
  long totalLen = 0;
  int maxLen = 0;
  int minLen = 0;
};

/// Violation is measured at si's own solution, which is also the point the
/// generator's own final filter uses -- so a cut counted here with zero violation
/// would have been erased by that filter and cannot appear.
static CutStats scoreCuts(const OsiCuts &cs, const double *x)
{
  CutStats s;
  s.rowCuts = cs.sizeRowCuts();
  s.colCuts = cs.sizeColCuts();
  for (int c = 0; c < s.rowCuts; ++c) {
    const OsiRowCut &rc = cs.rowCut(c);
    const CoinPackedVector &row = rc.row();
    const int len = row.getNumElements();
    const int *idx = row.getIndices();
    const double *el = row.getElements();
    double sum = 0.0;
    for (int k = 0; k < len; ++k)
      sum += el[k] * x[idx[k]];
    double viol = 0.0;
    if (sum > rc.ub())
      viol = sum - rc.ub();
    else if (sum < rc.lb())
      viol = rc.lb() - sum;
    s.totalViol += viol;
    s.maxViol = std::max(s.maxViol, viol);
    s.totalLen += len;
    s.maxLen = std::max(s.maxLen, len);
    s.minLen = (c == 0) ? len : std::min(s.minLen, len);
    if (rc.globallyValid())
      ++s.globallyValid;
  }
  return s;
}

static void usage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s --header\n"
    "       %s <fixture-stem> [options]\n"
    "\n"
    "Replays one Lagrangean Gomory separation call from a fixture written by\n"
    "CglLagomoryFixtureDump.hpp and prints one CSV row.\n"
    "\n"
    "The fixture is a pair of LPs: <stem>.mps.gz is the augmented LP whose\n"
    "trailing rows are the cuts to dualize, <stem>.orig.mps.gz is the original\n"
    "formulation handed to passInOriginalSolver. Both are required, along with\n"
    "<stem>.bas -- a Gomory cut is a row of the tableau at a specific basis, so\n"
    "there is nothing to replay without it.\n"
    "\n"
    "Options:\n"
    "  --gomory-type=N     composite type (default: the fixture's, else 21).\n"
    "                      N%%10: 1 dualize every cut row, 2 also copy the\n"
    "                      all-integral ones back in as real rows.\n"
    "                      N/10:  1 only when cut rows exist, 2 always.\n"
    "  --plain             ignore the original solver entirely and run plain\n"
    "                      Gomory on the augmented LP. This is the baseline the\n"
    "                      Lagrangean variants have to beat, and the only\n"
    "                      like-for-like one: same si, same basis, same knobs.\n"
    "  --rounds=N          separation rounds (default 1). Each round applies the\n"
    "                      previous round's cuts and resolves, as CBC does.\n"
    "  --repeat=N          time the call N times and report the minimum (default 1)\n"
    "  --limit=N           setLimit (max nonzeros in a cut in the tree)\n"
    "  --limit-at-root=N   setLimitAtRoot\n"
    "  --away=X            setAway (min fractionality in the tree)\n"
    "  --away-at-root=X    setAwayAtRoot\n"
    "  --alt-factorization=N  useAlternativeFactorization(N!=0)\n"
    "  --condition-mult=X  setConditionNumberMultiplier\n"
    "  --largest-mult=X    setLargestFactorMultiplier\n"
    "  --pass=N            CglTreeInfo::pass (default: the fixture's)\n"
    "  --options=N         CglTreeInfo::options (default: the fixture's)\n"
    "  --in-tree=0|1       CglTreeInfo::inTree (default: the fixture's)\n"
    "  --meta-knobs        take every generator knob from the .meta (default:\n"
    "                      library defaults, so a sweep is comparable across\n"
    "                      fixtures captured under different CBC settings)\n"
    "  --dump-cuts         print every cut with %%a coefficients. This, not the\n"
    "                      CSV, is the exactness gate for an optimization.\n"
    "  --profile           print per-stage times (needs Cgl -DCGL_GOMORY_PROFILE)\n"
    "  --quiet             suppress load warnings\n"
    "  --header            print the CSV header and exit\n",
    prog, prog);
}

static const char *CSV_HEADER
  = "name,rows,cols,elements,formulationRows,cutRows,cutRowNz,dualizedRows,"
    "integralCutRows,piNorm,fractionalInts,rowBoundDrift,origRows,origCols,"
    "origElements,mode,gomoryType,rounds,repeat,pass,options,inTree,limit,"
    "limitAtRoot,away,awayAtRoot,warmStartIters,piMaxDev,piCompared,solMaxDev,"
    "solCompared,objStart,objEnd,objImprove,rowCuts,colCuts,totalViol,maxViol,"
    "avgCutLen,maxCutLen,minCutLen,globallyValid,loadTime,sepTime";

int main(int argc, char *argv[])
{
  if (argc >= 2 && strcmp(argv[1], "--header") == 0) {
    printf("%s\n", CSV_HEADER);
    return 0;
  }
  if (argc < 2) {
    usage(argv[0]);
    return 2;
  }

  std::string stem;
  int gomoryType = -1;
  int rounds = 1, repeat = 1;
  int limit = -1, limitAtRoot = -1, altFactorization = -1;
  double away = -1.0, awayAtRoot = -1.0;
  double conditionMult = -1.0, largestMult = -1.0;
  int passOverride = -1, optionsOverride = -1, inTreeOverride = -1;
  bool plain = false, metaKnobs = false, dumpCutsFlag = false;
  bool profile = false, quiet = false;

  for (int i = 1; i < argc; ++i) {
    const std::string a(argv[i]);
    const char *eq = strchr(argv[i], '=');
    const char *val = eq ? eq + 1 : "";
    if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    } else if (a.rfind("--gomory-type=", 0) == 0) {
      gomoryType = atoi(val);
    } else if (a == "--plain") {
      plain = true;
    } else if (a.rfind("--rounds=", 0) == 0) {
      rounds = std::max(1, atoi(val));
    } else if (a.rfind("--repeat=", 0) == 0) {
      repeat = std::max(1, atoi(val));
    } else if (a.rfind("--limit=", 0) == 0) {
      limit = atoi(val);
    } else if (a.rfind("--limit-at-root=", 0) == 0) {
      limitAtRoot = atoi(val);
    } else if (a.rfind("--away=", 0) == 0) {
      away = atof(val);
    } else if (a.rfind("--away-at-root=", 0) == 0) {
      awayAtRoot = atof(val);
    } else if (a.rfind("--alt-factorization=", 0) == 0) {
      altFactorization = atoi(val);
    } else if (a.rfind("--condition-mult=", 0) == 0) {
      conditionMult = atof(val);
    } else if (a.rfind("--largest-mult=", 0) == 0) {
      largestMult = atof(val);
    } else if (a.rfind("--pass=", 0) == 0) {
      passOverride = atoi(val);
    } else if (a.rfind("--options=", 0) == 0) {
      optionsOverride = atoi(val);
    } else if (a.rfind("--in-tree=", 0) == 0) {
      inTreeOverride = atoi(val);
    } else if (a == "--meta-knobs") {
      metaKnobs = true;
    } else if (a == "--dump-cuts") {
      dumpCutsFlag = true;
    } else if (a == "--profile") {
      profile = true;
    } else if (a == "--quiet") {
      quiet = true;
    } else if (!a.empty() && a[0] == '-') {
      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      usage(argv[0]);
      return 2;
    } else if (stem.empty()) {
      stem = fixtureStem(argv[i]);
    } else {
      fprintf(stderr, "Only one fixture stem, got a second: %s\n", argv[i]);
      return 2;
    }
  }
  if (stem.empty()) {
    usage(argv[0]);
    return 2;
  }

  const std::string meta = stem + ".meta";
  const std::string name = baseName(stem);
  if (!fileExists(meta) && !quiet)
    fprintf(stderr, "WARNING: %s: no .meta; every default falls back to the "
                    "library's, not the capture's\n",
      name.c_str());

  const long metaFormulationRows = metaInt(meta, "formulationRows", -1);
  const double t0load = wallClock();

  // ---- the augmented LP -----------------------------------------------------
  const std::string mps
    = fileExists(stem + ".mps.gz") ? stem + ".mps.gz" : stem + ".mps";
  if (!fileExists(mps)) {
    fprintf(stderr, "ERROR: no problem file for stem %s\n", stem.c_str());
    return 1;
  }
  OsiClpSolverInterface si;
  ClpSimplex *lp = si.getModelPtr();
  lp->setLogLevel(0);
  if (si.readMps(mps.c_str())) {
    fprintf(stderr, "ERROR: failed to read %s\n", mps.c_str());
    return 1;
  }
  // Before the basis: the basis was written from the unpadded model, so with the
  // pad row still in place it is one artificial short and would not line up.
  dropPadRow(si, stem, "paddedColumns", "rows", "si", quiet);
  restoreColTypes(si, stem + ".ctype", name + " (si)", quiet);

  const std::string bas = stem + ".bas";
  if (!fileExists(bas)) {
    fprintf(stderr, "ERROR: %s: no .bas. A Gomory cut is a row of the tableau at "
                    "a specific basis; there is nothing to replay without it\n",
      name.c_str());
    return 1;
  }
  if (lp->readBasis(bas.c_str()) < 0) {
    fprintf(stderr, "ERROR: %s: failed to read basis %s\n", name.c_str(), bas.c_str());
    return 1;
  }
  // readBasis writes ClpSimplex::status_, but OsiClp caches a separate
  // CoinWarmStartBasis basis_ and resolve() overwrites the model from it.
  // setWarmStart(NULL) refreshes that cache from the model, which is what makes
  // the file's basis survive the solve -- and what makes si.getWarmStart(), where
  // CglGomory reads the basis, return the captured one rather than the cache.
  si.setWarmStart(NULL);
  lp->setPerturbation(50);
  si.setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
  si.setHintParam(OsiDoDualInResolve, true, OsiHintDo);
  si.resolve();
  const int warmStartIters = si.getIterationCount();
  if (!si.isProvenOptimal()) {
    fprintf(stderr, "ERROR: %s: augmented LP not optimal after warm start\n",
      name.c_str());
    return 1;
  }
  if (warmStartIters > 0 && !quiet)
    fprintf(stderr, "WARNING: %s: warm start took %d iterations; the captured "
                    "basis did not survive, so these are cuts from a DIFFERENT "
                    "basis and not comparable across runs\n",
      name.c_str(), warmStartIters);

  // ---- the original formulation --------------------------------------------
  // Loaded even in --plain mode, so the two modes report identical structural
  // columns and differ only in what they ran.
  OsiClpSolverInterface orig;
  bool haveOrig = false;
  {
    const std::string omps = fileExists(stem + ".orig.mps.gz")
      ? stem + ".orig.mps.gz"
      : stem + ".orig.mps";
    if (!fileExists(omps)) {
      if (!plain) {
        fprintf(stderr, "ERROR: %s: no .orig.mps.gz. Without the original\n"
                        "formulation there is no Lagrangean call to replay -- and\n"
                        "passing the augmented LP instead is exactly the mistake\n"
                        "that makes gomory-bench --orig-solver measure plain\n"
                        "Gomory. Use --plain if that baseline is what you want.\n",
          name.c_str());
        return 1;
      }
    } else {
      orig.getModelPtr()->setLogLevel(0);
      if (orig.readMps(omps.c_str())) {
        fprintf(stderr, "ERROR: failed to read %s\n", omps.c_str());
        return 1;
      }
      dropPadRow(orig, stem, "origPaddedColumns", "origRows", "orig", quiet);
      restoreColTypes(orig, stem + ".orig.ctype", name + " (orig)", quiet);
      // Solve once, cold. Not for its answer -- the Lagrangean pass overwrites
      // the objective, the column bounds, the primal solution and the basis
      // before it re-solves -- but because it allocates ClpSimplex's solution
      // arrays. memcpy(simplex->primalColumnSolution(), ...) in the generator has
      // nowhere to write on a model that has never been solved. It is also the
      // faithful state: CBC calls passInOriginalSolver with a solver it has
      // already solved.
      orig.initialSolve();
      haveOrig = true;
    }
  }

  if (haveOrig && metaFormulationRows >= 0
      && orig.getNumRows() != (int)metaFormulationRows && !quiet)
    fprintf(stderr, "WARNING: %s: orig has %d rows but meta says formulationRows "
                    "%ld; the dualized row range will not be the captured one\n",
      name.c_str(), orig.getNumRows(), metaFormulationRows);

  const double loadTime = wallClock() - t0load;

  // ---- fidelity -------------------------------------------------------------
  const int formulationRows = (metaFormulationRows >= 0)
    ? (int)metaFormulationRows
    : (haveOrig ? orig.getNumRows() : si.getNumRows());
  int piCompared = 0, solCompared = 0;
  const double piMaxDev = piDeviation(si, stem, formulationRows, piCompared);
  const double solMaxDev = solDeviation(si, stem, solCompared);

  // ---- the call -------------------------------------------------------------
  if (gomoryType < 0)
    gomoryType = (int)metaInt(meta, "gomoryType", 21);
  const int pass = (passOverride >= 0) ? passOverride : (int)metaInt(meta, "infoPass", 0);
  const int options
    = (optionsOverride >= 0) ? optionsOverride : (int)metaInt(meta, "infoOptions", 0);
  const int inTree
    = (inTreeOverride >= 0) ? inTreeOverride : (int)metaInt(meta, "infoInTree", 0);

  if (metaKnobs) {
    if (limit < 0)
      limit = (int)metaInt(meta, "limit", -1);
    if (limitAtRoot < 0)
      limitAtRoot = (int)metaInt(meta, "limitAtRoot", -1);
    /* dynamicLimitInTree is recorded in the .meta but deliberately not replayed:
       CglGomory recomputes it from numberColumns whenever limit_ is 0, and that
       is the only case in which it is ever read, so an injected value is always
       overwritten before use. It is a derived internal, not a knob. */
    if (altFactorization < 0)
      altFactorization = (int)metaInt(meta, "alternateFactorization", -1);
    if (away < 0.0)
      away = metaNum(meta, "away", -1.0);
    if (awayAtRoot < 0.0)
      awayAtRoot = metaNum(meta, "awayAtRoot", -1.0);
    if (conditionMult < 0.0)
      conditionMult = metaNum(meta, "conditionNumberMultiplier", -1.0);
    if (largestMult < 0.0)
      largestMult = metaNum(meta, "largestFactorMultiplier", -1.0);
  }

  // A working copy, so multi-round runs can apply cuts without disturbing the
  // solver the CSV's objStart and the violation scores refer to.
  OsiClpSolverInterface *work = dynamic_cast<OsiClpSolverInterface *>(si.clone());
  const double objStart = si.getObjValue();
  double objEnd = objStart;
  double sepTime = 0.0;
  CutStats stats;
  int effLimit = 0, effLimitAtRoot = 0;
  double effAway = 0.0, effAwayAtRoot = 0.0;

  for (int round = 0; round < rounds; ++round) {
    OsiCuts best;
    double bestTime = 0.0;

    for (int rep = 0; rep < repeat; ++rep) {
      // A fresh generator every repetition. CglGomory carries numberTimesStalled_
      // across calls and it appears in the whenToDo gate, so reusing one would
      // make repetition 3 a different experiment from repetition 1.
      CglGomory gomory;
      if (!plain && haveOrig) {
        // Order matters: passInOriginalSolver forces gomoryType_ to 1 when it is
        // still 0, so setting the type first would have it silently overwritten
        // for the default-constructed generator.
        gomory.passInOriginalSolver(&orig);
        gomory.setGomoryType(gomoryType);
      }
      if (limit >= 0)
        gomory.setLimit(limit);
      if (limitAtRoot >= 0)
        gomory.setLimitAtRoot(limitAtRoot);
      if (altFactorization >= 0)
        gomory.useAlternativeFactorization(altFactorization != 0);
      if (away >= 0.0)
        gomory.setAway(away);
      if (awayAtRoot >= 0.0)
        gomory.setAwayAtRoot(awayAtRoot);
      if (conditionMult >= 0.0)
        gomory.setConditionNumberMultiplier(conditionMult);
      if (largestMult >= 0.0)
        gomory.setLargestFactorMultiplier(largestMult);
      effLimit = gomory.getLimit();
      effLimitAtRoot = gomory.getLimitAtRoot();
      effAway = gomory.getAway();
      effAwayAtRoot = gomory.getAwayAtRoot();

      CglTreeInfo info;
      info.level = 0;
      info.pass = (round == 0) ? pass : pass + round;
      info.formulation_rows = formulationRows;
      info.inTree = (inTree != 0);
      info.options = options;

      OsiCuts cs;
      if (profile)
        LAG_PROF_RESET();
      const double t0 = wallClock();
      gomory.generateCuts(*work, cs, info);
      const double dt = wallClock() - t0;
      if (rep == 0 || dt < bestTime) {
        bestTime = dt;
        best = cs;
      }
      if (profile) {
        char tag[256];
        snprintf(tag, sizeof(tag), "%s round %d rep %d", name.c_str(), round, rep);
        LAG_PROF_PRINT(tag);
      }
    }

    sepTime += bestTime;
    const CutStats rs = scoreCuts(best, work->getColSolution());
    stats.rowCuts += rs.rowCuts;
    stats.colCuts += rs.colCuts;
    stats.globallyValid += rs.globallyValid;
    stats.totalViol += rs.totalViol;
    stats.maxViol = std::max(stats.maxViol, rs.maxViol);
    stats.totalLen += rs.totalLen;
    stats.maxLen = std::max(stats.maxLen, rs.maxLen);
    stats.minLen = (round == 0) ? rs.minLen
                                : (rs.rowCuts ? std::min(stats.minLen, rs.minLen)
                                              : stats.minLen);
    if (dumpCutsFlag)
      dumpCuts(best, name, round);

    // Apply and resolve, so objImprove is the bound movement the cuts actually
    // buy -- the metric that decides whether a change is an improvement. Cut
    // count is not: a change that produces more, weaker cuts looks better by
    // count and worse by bound.
    if (best.sizeRowCuts() || best.sizeColCuts()) {
      work->applyCuts(best);
      work->resolve();
      if (work->isProvenOptimal())
        objEnd = work->getObjValue();
      else
        break; // an infeasible or unbounded relaxation ends the sequence
    } else {
      break; // nothing to apply, so later rounds would repeat this one exactly
    }
  }

  const double objImprove = objEnd - objStart;
  const double avgCutLen = stats.rowCuts ? (double)stats.totalLen / stats.rowCuts : 0.0;

  printf("%s,%d,%d,%d,%d,%ld,%ld,%ld,%ld,%.15g,%ld,%ld,%d,%d,%d,%s,%d,%d,%d,%d,%d,"
         "%d,%d,%d,%.10g,%.10g,%d,%.10g,%d,%.10g,%d,%.15g,%.15g,%.15g,%d,%d,%.15g,"
         "%.15g,%.4f,%d,%d,%d,%.6f,%.6f\n",
    name.c_str(), si.getNumRows(), si.getNumCols(), si.getNumElements(),
    formulationRows, metaInt(meta, "cutRows", -1), metaInt(meta, "cutRowNz", -1),
    metaInt(meta, "dualizedRows", -1), metaInt(meta, "integralCutRows", -1),
    metaNum(meta, "piNorm", -1.0), metaInt(meta, "fractionalInts", -1),
    metaInt(meta, "rowBoundDrift", -1),
    haveOrig ? orig.getNumRows() : -1, haveOrig ? orig.getNumCols() : -1,
    haveOrig ? orig.getNumElements() : -1,
    (plain || !haveOrig) ? "plain" : "lagomory",
    (plain || !haveOrig) ? 0 : gomoryType, rounds, repeat, pass, options, inTree,
    effLimit, effLimitAtRoot, effAway, effAwayAtRoot, warmStartIters,
    piMaxDev, piCompared, solMaxDev, solCompared,
    objStart, objEnd, objImprove, stats.rowCuts, stats.colCuts, stats.totalViol,
    stats.maxViol, avgCutLen, stats.maxLen, stats.minLen, stats.globallyValid,
    loadTime, sepTime);

  delete work;
  return 0;
}
