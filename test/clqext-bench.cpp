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
 * **Under method 4 the counts must reproduce exactly, and a mismatch means the
 * basis was lost.** An earlier version of this comment claimed the opposite --
 * that these LPs are degenerate, so re-solving legitimately lands on a different
 * optimal vertex with different reduced costs, and that the 5 fixtures of 660
 * which disagreed with Cbc's log were therefore expected noise. The mechanism was
 * real but the conclusion was wrong: the basis was not being restored *at all*
 * (an MPS/`.bas` name mismatch, then OsiClp's cached `basis_` overwriting the
 * model -- see the warm-start block below), so every replay started from an
 * all-slack basis and the "seven ways" evidence was seven cold solves. With both
 * halves fixed all five match to the row: momentum1 73/118, mzzv11 34/67,
 * neos-3216931-puriri 26/43, neos8 146/188, trdHardCA4500-8 41/109.
 *
 * What survives is the *sensitivity*, which is why warmStartIters is reported and
 * warned about: reduced costs really do drive CoinCliqueExtender's candidate
 * ordering, and replaying momentum1 with --no-lp (no reduced costs, so method 4
 * degrades to 2) gives 72/135 instead of 73/118. So treat a nonzero
 * warmStartIters as a defect to chase, not as tolerance to budget for. One
 * caveat: 0 is the expectation, not a guarantee -- momentum1 needs 11 pivots and
 * still reproduces Cbc exactly, repeatably. Nonzero means "verify before
 * trusting", not "wrong".
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
 * Read one numeric key out of the fixture's `.meta` file.
 * Returns `dflt` when the file or the key is absent, so a fixture written before
 * `.meta` existed -- or before it carried that key -- still runs.
 */
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

/**
 * Restore integrality from the `.ctype` sidecar.
 *
 * MPS cannot express "fixed and integer": CoinMpsIO::writeMps conveys
 * integrality only through the bound type (BV/UI/LI/MI+UI) and a column with
 * lb == ub takes the " FX " branch, which has no integer form -- so every integer
 * column CBC had fixed by bound tightening reads back continuous. That silently
 * changes the experiment rather than merely annotating it, because integrality is
 * a gate on both BK consumers: CglCliqueStrengthening::detectCliqueRows rejects
 * any row containing a non-binary column, so lost markers shrink the clique-row
 * set. Measured on physiciansched3-3: 39994 markers lost, clique rows 188486 ->
 * 83611, replayed extensions 4731 -> 2972 against what CBC did.
 *
 * Returns the number of columns re-marked, or -1 when no usable sidecar was
 * found. A sidecar whose column count disagrees with the loaded model is refused
 * rather than partly applied: a wrong-model sidecar would mark arbitrary columns
 * integer, which is worse than the loss it is meant to repair.
 */
static int restoreColTypes(OsiSolverInterface &si, const std::string &stem, bool quiet)
{
  const std::string path = stem + ".ctype";
  FILE *fp = fopen(path.c_str(), "r");
  if (!fp) {
    if (!quiet)
      fprintf(stderr, "WARNING: %s: no .ctype sidecar; integer columns that were "
                      "fixed at capture will read back continuous\n",
        baseName(stem).c_str());
    return -1;
  }

  int sidecarCols = -1;
  if (fscanf(fp, "cols %d\n", &sidecarCols) != 1 || sidecarCols != si.getNumCols()) {
    fprintf(stderr, "ERROR: %s: .ctype is for %d columns, model has %d; ignoring it\n",
      baseName(stem).c_str(), sidecarCols, si.getNumCols());
    fclose(fp);
    return -1;
  }

  int idx = 0, type = 0, restored = 0, alreadyInteger = 0;
  while (fscanf(fp, "%d %d\n", &idx, &type) == 2) {
    if (idx < 0 || idx >= si.getNumCols()) {
      fprintf(stderr, "ERROR: %s: .ctype names column %d, out of range\n",
        baseName(stem).c_str(), idx);
      fclose(fp);
      return -1;
    }
    if (si.isContinuous(idx)) {
      si.setInteger(idx);
      ++restored;
    } else {
      ++alreadyInteger;
    }
  }
  fclose(fp);

  // Bounds already carried the marker for these; only the fixed ones needed help.
  (void)alreadyInteger;
  // getColType() caches, and derives Binary vs GeneralInteger from the bounds, so
  // it has to be recomputed after this or consumers would see the pre-restore view.
  si.getColType(true);
  return restored;
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
    "  --max-seconds=F     strengthening wall-clock limit; 0 = none. Default is the\n"
    "                      limit the captured site imposed, from .meta (Cbc gives\n"
    "                      strengthening its remaining budget), or 0 if unrecorded\n"
    "  --force-lp          solve the LP even if the fixture was captured without\n"
    "                      one, so methods 4/5 get the reduced costs they want\n"
    "                      (diverges from what Cbc did; reports effMethod)\n"
    "  --no-lp             never solve the warm-start LP, even when the fixture\n"
    "                      has one. Methods 4/5 then downgrade to 2, so this\n"
    "                      measures a different algorithm -- but it is the way to\n"
    "                      time strengthening on a fixture whose LP does not\n"
    "                      converge (z26, cdc7-4-3-2, scpn2)\n"
    "  --lp-seconds=F      cap the warm-start LP; 0 = none (default 300). Those\n"
    "                      same fixtures do not reach optimality at any cap, cold\n"
    "                      or warm, so the cap turns a hang into a clear error\n"
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
    "boundMoved,resolveTime,resolveIters,warmStartTime,warmStartIters,cgraphTime,cgNodes,"
    "cgDirectConf,cgCliques,cgDensity,restoredInt,maxSeconds";

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
  bool forceLp = false, noLp = false;
  // -1 until either --ext-method or the fixture's .meta supplies one, so an
  // explicit request can be told apart from the default.
  long extMethodArg = -1;
  // Negative until either --max-seconds or the fixture's .meta supplies one; 0 is
  // a meaningful value here (no limit) and so cannot double as "unset".
  double maxSecondsArg = -1.0;
  // Cap on the warm-start LP. Not unlimited by default: see the solve site.
  double lpSeconds = 300.0;
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
    } else if (strcmp(a, "--no-lp") == 0) {
      noLp = true;
    } else if (strncmp(a, "--ext-method=", 13) == 0) {
      extMethodArg = atol(a + 13);
    } else if (strncmp(a, "--max-seconds=", 14) == 0) {
      maxSecondsArg = atof(a + 14);
    } else if (strncmp(a, "--lp-seconds=", 13) == 0) {
      lpSeconds = atof(a + 13);
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

  // Before anything reads a column type -- which detectCliqueRows does, through
  // the CglCliqueStrengthening constructor.
  const int restoredCols = restoreColTypes(si, stem, quiet);

  const std::string meta = stem + ".meta";
  const int lpOptimalInFixture = (int)metaInt(meta, "lpOptimal", -1);
  const long metaExtMethod = metaInt(meta, "extMethodRequested", -1);
  // Explicit flag first, then the fixture's own request, then Cbc's "after"
  // default for fixtures written before .meta recorded the method.
  const size_t extMethod = (size_t)(extMethodArg >= 0 ? extMethodArg
                                                      : (metaExtMethod >= 0 ? metaExtMethod : 4));
  // Cbc gives strengthening whatever is left of the model's time budget, so an
  // unbounded replay is not the same experiment: measured, cdc7-4-3-2, scpm1,
  // scpn2 and z26 were all cut short inside Cbc and logged 0 extended / 0
  // dominated, while an unbounded replay of the same fixture ran past 900s and
  // found extensions. Honour the captured limit by default, so the replay matches
  // the capture; --max-seconds overrides it, including back to 0 for no limit.
  const double maxSeconds
    = maxSecondsArg >= 0.0 ? maxSecondsArg : metaNum(meta, "maxSeconds", 0.0);
  const bool wantLp = !noLp
    && (forceLp || lpOptimalInFixture == 1 || (lpOptimalInFixture < 0 && fileExists(bas)));

  double warmStartTime = 0.0;
  int warmStartIters = 0;
  if (wantLp) {
    // Two steps, and skipping either leaves the stored basis inert while still
    // looking like success. `readBasis` writes ClpSimplex's own `status_`, but
    // OsiClp caches a separate `CoinWarmStartBasis basis_` and installs it over
    // the model in `resolve()` (OsiClpSolverInterface.cpp:1199), while
    // `initialSolve()` presolves from scratch and discards the basis outright.
    // `setWarmStart(NULL)` refreshes the cache from the model (`basis_ =
    // getBasis(modelPtr_)`, applying the slack flip), and `resolve()` then pivots
    // from it. With both in place an already-optimal fixture costs 0 iterations,
    // which is the cheap check that the basis survived -- anything else is a
    // different vertex, hence different reduced costs and a different candidate
    // ordering inside CoinCliqueExtender.
    bool haveBasis = false;
    if (fileExists(bas)) {
      if (lp->readBasis(bas.c_str()) < 0) {
        fprintf(stderr, "WARNING: failed to read basis %s; solving cold\n", bas.c_str());
      } else {
        haveBasis = true;
        si.setWarmStart(NULL);
      }
    } else if (!quiet) {
      fprintf(stderr, "WARNING: no basis %s; solving cold\n", bas.c_str());
    }

    lp->setPerturbation(50);
    // Bounded, because on a few fixtures this LP does not converge in any
    // practical time and an unbounded solve makes the whole tool look hung on
    // them. Measured on z26 (38223x17752): the warm start runs 62948 iterations
    // in 102s without reaching optimality, and a *cold* solve of the same model
    // does no better (152s, also not optimal) -- so it is the LP that is
    // expensive, not the stored basis that is bad. cdc7-4-3-2 behaves the same
    // way (65s warm, 75s cold, neither optimal). Strengthening on both takes
    // 0.00s and reproduces Cbc's 0 extended / 0 dominated exactly, so the LP is
    // the only slow part and capping it costs nothing but the reduced costs.
    if (lpSeconds > 0.0)
      lp->setMaximumSeconds(lpSeconds);

    const double t0 = wallClock();
    if (haveBasis) {
      si.setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
      si.setHintParam(OsiDoDualInResolve, true, OsiHintDo);
      si.resolve();
    } else {
      si.setHintParam(OsiDoDualInInitial, true, OsiHintDo);
      si.initialSolve();
    }
    warmStartTime = wallClock() - t0;
    warmStartIters = si.getIterationCount();

    if (haveBasis && si.isProvenOptimal() && warmStartIters > 0 && !quiet) {
      fprintf(stderr, "WARNING: %s: warm start took %d iterations; the captured "
                      "basis did not survive, so this is a different vertex\n",
        baseName(stem).c_str(), warmStartIters);
    }

    if (!si.isProvenOptimal()) {
      // An error rather than a warning: the fixture says an optimal LP was
      // available, so failing to reach it means methods 4/5 would silently
      // become 2 and the run would measure something other than what was asked.
      fprintf(stderr, "ERROR: LP not optimal after warm start (%s) in %.1fs; "
                      "reduced costs unavailable, ext-method %lu would downgrade "
                      "to 2. Raise --lp-seconds, or use --no-lp to measure "
                      "strengthening without them.\n",
        stem.c_str(), warmStartTime, (unsigned long)extMethod);
      return 1;
    }
  } else if (!quiet) {
    if (noLp)
      fprintf(stderr, "NOTE: %s: --no-lp, so ext-method %lu runs without reduced "
                      "costs (4 and 5 downgrade to 2)\n",
        baseName(stem).c_str(), (unsigned long)extMethod);
    else
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
         "%.6f,%d,%.6f,%lu,%lu,%lu,%.10f,%d,%.6f\n",
    baseName(stem).c_str(), (unsigned long)extMethod, (unsigned long)effMethod,
    (int)wantLp, extended, dominated,
    strTime, rowsBefore, si.getNumRows(), nzBefore, si.getNumElements(),
    objBefore, objAfter, objImprove, objImproveRel,
    (int)boundMoved(objBefore, objImprove), resolveTime, resolveIters,
    warmStartTime, warmStartIters, cgraphTime,
    (unsigned long)cgNodes, (unsigned long)cgDirectConf,
    (unsigned long)cgCliques, cgDensity, restoredCols, maxSeconds);

  return 0;
}
