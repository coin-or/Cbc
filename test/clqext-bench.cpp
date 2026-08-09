/**
 * Stand-alone benchmark for clique extension / merging.
 *
 * @file clqext-bench.cpp
 * @brief replay a CglCliqueStrengthening call captured from CBC, without CBC
 *
 * The companion to bkclique-bench: same fixture format, the other consumer of
 * Bron-Kerbosch inside CBC. Loads a fixture written by CbcClqFixtureDump.hpp at
 * a "clqstr-*" tag and runs CglCliqueStrengthening on it.
 *
 * Extension does not produce violated cuts -- it rewrites set-packing rows in
 * place and drops the ones the extended rows dominate -- so there is no violation
 * sum here. That leaves the bound as the only real measure of quality, and it is
 * the headline:
 *
 *   - **objImprove**, the LP bound movement on re-solving the strengthened
 *     formulation, is the definitive measure. `objImproveRel` scales it by the
 *     starting bound so instances are comparable, and `boundMoved` says whether
 *     it clears the noise floor at all -- see boundMoved() for why a raw
 *     objImprove must never be read against zero.
 *   - rows extended and rows dominated (the generator's own counters) say how
 *     much work happened, not how much good it did. High counts at
 *     `boundMoved` 0 mean the rewrite was busywork.
 *   - the change in matrix size, since extension adds nonzeros and merging
 *     removes rows -- a strengthening that adds many nonzeros for few dominated
 *     rows is a worse trade even at equal counts.
 *   - time, the remaining concern, and only meaningful once the bound holds.
 *
 * Extension methods 4 and 5 need reduced costs, and
 * CglCliqueStrengthening::getReducedCost() returns NULL unless the model is
 * proven optimal -- in which case those methods silently degrade to 2.
 *
 * Cbc sidesteps that by asking for different methods in its two passes: mode
 * "before" runs on the original problem with no LP solved and requests method 2,
 * mode "after" runs on the preprocessed problem with an optimal LP and requests
 * method 4. So defaulting this tool to 4 would misrepresent the "before" pass
 * twice over -- wrong method, and silently degraded at that. Instead both the
 * requested method and the LP status are read from the fixture's `.meta`, and the
 * default is whatever the captured site asked for. --ext-method overrides it to
 * compare methods on a fixture, and --force-lp solves an LP the capture did not
 * have so 4 and 5 get their reduced costs; effMethod reports what actually ran.
 *
 * A no-LP fixture still needs a "before" bound, or the bound columns say nothing.
 * That one is taken from a *clone*, so the solver being measured is left exactly
 * as captured -- reduced-cost-free. `lpOptimal` in the output is 0 whenever
 * objBefore came from that clone rather than from the fixture itself.
 *
 * Usage:
 *   clqext-bench <stem> [options]           (stem = dir/name.tag)
 */

#include "CglCliqueStrengthening.hpp"
#include "ClpSimplex.hpp"
#include "CoinStaticConflictGraph.hpp"
#include "OsiClpSolverInterface.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/stat.h>

static double wallClock()
{
  return CoinGetTimeOfDay();
}

static bool fileExists(const std::string &path)
{
  struct stat st;
  return stat(path.c_str(), &st) == 0;
}

/// Reduce any of the fixture's file names to the shared stem.
static std::string fixtureStem(const char *arg)
{
  std::string s(arg);
  static const char *suffixes[]
    = { ".mps.gz", ".mps", ".cgraph", ".bas", ".sol", ".bas.status" };
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

/**
 * Read one integer key out of the fixture's `.meta` file.
 * Returns `dflt` when the file or the key is absent, so a fixture written before
 * `.meta` existed still runs.
 */
static long metaInt(const std::string &path, const char *key, long dflt)
{
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp)
    return dflt;

  char line[512];
  long value = dflt;
  while (fgets(line, sizeof(line), fp)) {
    char k[256];
    double v = 0.0;
    if (sscanf(line, "%255s %lf", k, &v) == 2 && strcmp(k, key) == 0) {
      value = (long)v;
      break;
    }
  }
  fclose(fp);
  return value;
}

/**
 * Remove the pad row, if this fixture has one.
 *
 * A capture whose matrix held empty columns carries one extra final row, without
 * which writeMps would drop those columns and shift every index after them --
 * see cbcClqFixtureWriteMps. The row is redundant, so it would not change the LP,
 * but it is a row the captured model did not have, and row selection inside
 * CglCliqueStrengthening walks every row -- so leaving it in would change what
 * gets considered. Only clqstr-before captures are ever padded, since
 * preprocessing removes empty columns itself.
 *
 * Both conditions are required before deleting anything: `.meta` must say padding
 * happened, *and* the loaded row count must be exactly one more than the captured
 * one. Either alone could delete a real row from a fixture whose meta is stale.
 */
static bool dropPadRow(OsiSolverInterface &si, const std::string &stem, bool quiet)
{
  const std::string meta = stem + ".meta";
  if (metaInt(meta, "paddedColumns", 0) <= 0)
    return false;

  const long capturedRows = metaInt(meta, "rows", -1);
  if (capturedRows < 0 || si.getNumRows() != (int)capturedRows + 1) {
    if (!quiet)
      fprintf(stderr, "WARNING: %s: meta says padded but rows=%d against captured %ld; "
                      "leaving the matrix alone\n",
        baseName(stem).c_str(), si.getNumRows(), capturedRows);
    return false;
  }

  const int last = si.getNumRows() - 1;
  si.deleteRows(1, &last);
  return true;
}

static void usage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s <fixture-stem> [options]\n"
    "\n"
    "<fixture-stem> is dir/name.tag; .mps.gz/.cgraph/.bas are appended.\n"
    "\n"
    "Options:\n"
    "  --ext-method=N      0=none 1=random 2=max degree 3=max modified degree\n"
    "                      4=reduced cost 5=reduced cost+mdegree. Default is the\n"
    "                      method the captured site asked for (2 for clqstr-before,\n"
    "                      4 for clqstr-after), from .meta; 4 if unrecorded\n"
    "  --rebuild-cgraph    rebuild the graph from the matrix instead of loading\n"
    "                      the captured one (not faithful; for comparison only)\n"
    "  --max-seconds=F     strengthening wall-clock limit (0 = none, default)\n"
    "  --force-lp          solve the LP even if the fixture was captured without\n"
    "                      one, so methods 4/5 get the reduced costs they want\n"
    "                      (diverges from what Cbc did; reports effMethod)\n"
    "  --no-resolve        skip the after-LP; report objAfter as objBefore\n"
    "  --header            print the CSV header line and exit\n"
    "  --csv-header        print the CSV header before the data line\n"
    "  --quiet             suppress warnings\n",
    prog);
}

/**
 * Did the bound actually move, or is this floating-point noise?
 *
 * objImprove is a difference of two LP objectives, so it inherits their absolute
 * error, which scales with their magnitude -- an LP around 1e6 carries roughly
 * 1e-10 of slop. Reading such a difference against zero therefore reports an
 * improvement on every fixture: measured, brazil3's clqstr-before pass gives
 * 2.7e-12 over 29317 iterations, which is nothing at all, while mzzv11 gives
 * 171.49, which is real. A relative test separates the two, with an absolute
 * floor so that an objective near zero does not make the relative test
 * hypersensitive.
 *
 * Kept identical to bkclique-bench's, so a bound gain means the same thing in
 * both tools.
 */
static bool boundMoved(double objBefore, double objImprove)
{
  const double scale = fabs(objBefore) > 1.0 ? fabs(objBefore) : 1.0;
  return objImprove > 1.0e-9 * scale && objImprove > 1.0e-9;
}

static const char *CSV_HEADER
  = "name,extMethod,effMethod,lpOptimal,extended,dominated,strTime,rowsBefore,"
    "rowsAfter,nzBefore,nzAfter,objBefore,objAfter,objImprove,objImproveRel,"
    "boundMoved,resolveTime,resolveIters,warmStartTime,cgraphTime,cgNodes,"
    "cgDirectConf,cgCliques,cgDensity";

int main(int argc, char *argv[])
{
  if (argc < 2) {
    usage(argv[0]);
    return 1;
  }
  if (strcmp(argv[1], "--header") == 0) {
    printf("%s\n", CSV_HEADER);
    return 0;
  }
  if (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
    usage(argv[0]);
    return 0;
  }

  bool rebuildCgraph = false, csvHeader = false, quiet = false, resolveAfter = true;
  bool forceLp = false;
  // -1 until either --ext-method or the fixture's .meta supplies one, so an
  // explicit request can be told apart from the default.
  long extMethodArg = -1;
  double maxSeconds = 0.0;
  const char *stemArg = NULL;

  for (int i = 1; i < argc; ++i) {
    const char *a = argv[i];
    if (strcmp(a, "--rebuild-cgraph") == 0) {
      rebuildCgraph = true;
    } else if (strcmp(a, "--csv-header") == 0) {
      csvHeader = true;
    } else if (strcmp(a, "--quiet") == 0) {
      quiet = true;
    } else if (strcmp(a, "--no-resolve") == 0) {
      resolveAfter = false;
    } else if (strcmp(a, "--force-lp") == 0) {
      forceLp = true;
    } else if (strncmp(a, "--ext-method=", 13) == 0) {
      extMethodArg = atol(a + 13);
    } else if (strncmp(a, "--max-seconds=", 14) == 0) {
      maxSeconds = atof(a + 14);
    } else if (a[0] == '-') {
      fprintf(stderr, "ERROR: unknown option %s\n", a);
      usage(argv[0]);
      return 1;
    } else {
      stemArg = a;
    }
  }

  if (!stemArg) {
    fprintf(stderr, "ERROR: no fixture stem given\n");
    usage(argv[0]);
    return 1;
  }

  const std::string stem = fixtureStem(stemArg);
  const std::string mps
    = fileExists(stem + ".mps.gz") ? stem + ".mps.gz" : stem + ".mps";
  const std::string bas = stem + ".bas";
  const std::string cgr = stem + ".cgraph";

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

  // Two things the fixture knows and the files alone do not: whether the captured
  // solver had an optimal LP, and which extension method that site asked for.
  // Extension on the original problem ("before") runs before any LP is solved, so
  // solving one here would hand methods 4/5 reduced costs Cbc did not have.
  // Before anything reads a row count or a basis: the pad row is not part of the
  // captured model.
  dropPadRow(si, stem, quiet);

  const std::string meta = stem + ".meta";
  const int lpOptimalInFixture = (int)metaInt(meta, "lpOptimal", -1);
  const long metaExtMethod = metaInt(meta, "extMethodRequested", -1);
  // Explicit flag first, then the fixture's own request, then Cbc's "after"
  // default for fixtures written before .meta recorded the method.
  const size_t extMethod = (size_t)(extMethodArg >= 0 ? extMethodArg
                                                      : (metaExtMethod >= 0 ? metaExtMethod : 4));
  const bool wantLp
    = forceLp || lpOptimalInFixture == 1 || (lpOptimalInFixture < 0 && fileExists(bas));

  double warmStartTime = 0.0;
  if (wantLp) {
    if (fileExists(bas)) {
      if (lp->readBasis(bas.c_str()) < 0)
        fprintf(stderr, "WARNING: failed to read basis %s; solving cold\n", bas.c_str());
    } else if (!quiet) {
      fprintf(stderr, "WARNING: no basis %s; solving cold\n", bas.c_str());
    }

    lp->setPerturbation(50);
    si.setHintParam(OsiDoDualInInitial, true, OsiHintDo);

    const double t0 = wallClock();
    si.initialSolve();
    warmStartTime = wallClock() - t0;

    if (!si.isProvenOptimal()) {
      // An error rather than a warning: the fixture says an optimal LP was
      // available, so failing to reach it means methods 4/5 would silently
      // become 2 and the run would measure something other than what was asked.
      fprintf(stderr, "ERROR: LP not optimal after warm start (%s); reduced costs "
                      "unavailable, ext-method %lu would downgrade to 2\n",
        stem.c_str(), (unsigned long)extMethod);
      return 1;
    }
  } else if (!quiet) {
    fprintf(stderr, "NOTE: %s captured without an optimal LP; not solving one "
                    "(pass --force-lp to override)\n",
      baseName(stem).c_str());
  }

  // What CglCliqueStrengthening will actually use, after its own downgrade.
  const size_t effMethod
    = (!si.isProvenOptimal() && (extMethod == 4 || extMethod == 5)) ? 2 : extMethod;

  const double t1 = wallClock();
  if (rebuildCgraph) {
    si.checkCGraph(NULL);
  } else {
    if (!fileExists(cgr)) {
      fprintf(stderr, "ERROR: no conflict graph %s (use --rebuild-cgraph to build one)\n",
        cgr.c_str());
      return 1;
    }
    CoinStaticConflictGraph *cg = CoinStaticConflictGraph::load(cgr.c_str());
    if (!cg) {
      fprintf(stderr, "ERROR: failed to load conflict graph %s\n", cgr.c_str());
      return 1;
    }
    if (cg->size() != (size_t)si.getNumCols() * 2) {
      fprintf(stderr, "ERROR: graph/model mismatch for %s: graph %lu nodes, model %d cols\n",
        stem.c_str(), (unsigned long)cg->size(), si.getNumCols());
      delete cg;
      return 1;
    }
    si.setCGraph(cg);
  }
  const double cgraphTime = wallClock() - t1;

  const CoinConflictGraph *cg = si.getCGraph();
  if (!cg) {
    fprintf(stderr, "ERROR: no conflict graph available for %s\n", stem.c_str());
    return 1;
  }

  // Captured before strengthening: it deletes dominated rows and adds extended
  // ones, so none of these survive the call.
  const int rowsBefore = si.getNumRows();
  const int nzBefore = si.getNumElements();
  const size_t cgNodes = cg->size();
  const size_t cgDirectConf = cg->nTotalDirectConflicts();
  const size_t cgCliques = cg->nCliques();
  const double cgDensity = cg->density();

  // The before-bound. When the fixture carried an LP, si already holds it. When
  // it did not, solve a *clone*: si itself must stay reduced-cost-free, or
  // methods 4/5 would silently stop downgrading and the run would measure a
  // different algorithm than Cbc ran. objBefore is then only a yardstick for
  // objAfter, which is why lpOptimal is reported alongside it.
  double objBefore = si.getObjValue();
  if (!si.isProvenOptimal() && resolveAfter) {
    OsiSolverInterface *probe = si.clone();
    probe->messageHandler()->setLogLevel(0);
    probe->initialSolve();
    if (probe->isProvenOptimal())
      objBefore = probe->getObjValue();
    else if (!quiet)
      fprintf(stderr, "WARNING: %s: LP not optimal before strengthening\n", stem.c_str());
    delete probe;
  }

  CglCliqueStrengthening clqStr(&si);
  clqStr.setMaximumSeconds(maxSeconds);

  const double t2 = wallClock();
  clqStr.strengthenCliques(extMethod);
  const double strTime = wallClock() - t2;

  const int extended = clqStr.constraintsExtended();
  const int dominated = clqStr.constraintsDominated();

  double objAfter = objBefore, resolveTime = 0.0;
  int resolveIters = 0;
  // Worth solving even when the fixture had no LP: without a bound before and
  // after there is nothing to say about whether the rewrite tightened anything.
  // Strengthening is already done by this point, so warm-starting si here cannot
  // feed reduced costs back into it.
  if (resolveAfter && (extended > 0 || dominated > 0)) {
    const double t3 = wallClock();
    if (wantLp)
      si.resolve();
    else
      si.initialSolve();
    resolveTime = wallClock() - t3;
    resolveIters = si.getIterationCount();
    if (si.isProvenOptimal())
      objAfter = si.getObjValue();
    else if (!quiet)
      fprintf(stderr, "WARNING: %s: LP not optimal after strengthening\n", stem.c_str());
  }

  const double objImprove = si.getObjSenseInCbc() * (objAfter - objBefore);
  // Relative to the starting bound's magnitude, which is what makes gains
  // comparable across instances whose objectives differ by orders of magnitude.
  const double objImproveRel
    = objImprove / (fabs(objBefore) > 1.0 ? fabs(objBefore) : 1.0);

  if (csvHeader)
    printf("%s\n", CSV_HEADER);

  printf("%s,%lu,%lu,%d,%d,%d,%.6f,%d,%d,%d,%d,%.15g,%.15g,%.6g,%.6g,%d,%.6f,%d,"
         "%.6f,%.6f,%lu,%lu,%lu,%.10f\n",
    baseName(stem).c_str(), (unsigned long)extMethod, (unsigned long)effMethod,
    (int)wantLp, extended, dominated,
    strTime, rowsBefore, si.getNumRows(), nzBefore, si.getNumElements(),
    objBefore, objAfter, objImprove, objImproveRel,
    (int)boundMoved(objBefore, objImprove), resolveTime, resolveIters,
    warmStartTime, cgraphTime,
    (unsigned long)cgNodes, (unsigned long)cgDirectConf,
    (unsigned long)cgCliques, cgDensity);

  return 0;
}
