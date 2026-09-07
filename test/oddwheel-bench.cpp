/**
 * Stand-alone benchmark for odd-wheel (odd-cycle) separation.
 *
 * @file oddwheel-bench.cpp
 * @brief replay a CglOddWheel call captured from CBC, without CBC
 *
 * Loads a fixture written by CbcClqFixtureDump.hpp -- preprocessed problem,
 * optimal LP basis, serialized conflict graph -- warm-starts the LP, and runs N
 * rounds of CglOddWheel. The fixtures were captured for clique separation, and
 * they serve here unchanged because odd-wheel separation consumes exactly the
 * same three things: the conflict graph, the LP solution and the reduced costs.
 * Reaching that state inside CBC costs a full solve, which is what this exists
 * to avoid paying on every iteration of an optimization loop.
 *
 * Figures of merit, in order of authority -- the same order as bkclique-bench,
 * and for the same reasons:
 *
 *   - **bound improvement on reoptimizing** (`objImprove`). The only figure that
 *     says the cuts tightened the relaxation, which is the point of finding them.
 *   - total **violation** (`totalViol`), the proxy available before a re-solve.
 *   - separation **time**, since separation runs many times per solve.
 *
 * The cut count is reported and is the weakest of the four: more cuts at equal
 * bound movement is strictly worse.
 *
 * ## Why this drives CglOddWheel and not CoinOddWheelSeparator
 *
 * The separator could be constructed directly, and that would be a mistake. A
 * standalone probe measures whatever structure the probe itself builds, which is
 * not necessarily the one the real call site builds -- for clique separation a
 * probe of this shape reported 0.95s against ~19ms actually spent, because it
 * built the full conflict graph where the real path uses a small induced
 * subgraph. So the generator is driven through its public entry point, and the
 * separator's internals are read back through CglOddWheel::stats(), which
 * carries CoinOddWheelSeparator::Stats from the last call.
 *
 * ## The one flag that matters most
 *
 * `--ext-method=0` disables the wheel-center lifting entirely
 * (`searchOddWheels()` gates the lifting loop on `extMethod_ > 0`), so
 * `--ext-method=2` minus `--ext-method=0` isolates the cost of lifting exactly.
 * The analogous control is what located the real bottleneck in clique
 * separation in a single sweep -- there it was extension, not the search.
 *
 * ## Reproducibility
 *
 * `--max-seconds` defaults to 0, i.e. no limit, deliberately. CbcCutGenerator
 * hands the generator whatever wall clock remains, and under such a budget the
 * separator aborts partway through graph preparation on a large graph -- so a
 * nonzero limit makes the result depend on machine load and the run stops being
 * reproducible. With the limit off, the counters are a function of the fixture
 * alone.
 *
 * There is no `--min-viol`: the threshold is a compile-time constant in
 * CoinOddWheelSeparator.cpp (ODDWHEEL_SEP_DEF_MIN_VIOL) with no setter, and
 * inventing one just to have a flag would change the code under measurement.
 *
 * The graph is loaded rather than rebuilt by default; a rebuild does not
 * reproduce the graph CBC separated on (measured: trdtaunimep 190900 direct
 * conflicts rebuilt against 191028 used). --rebuild-cgraph selects the old
 * behaviour so the two can be compared directly.
 *
 * Usage:
 *   oddwheel-bench <stem> [options]          (stem = dir/name.tag)
 *   oddwheel-bench --self-test <stem>
 *
 * <stem> is the fixture prefix: <stem>.mps.gz, <stem>.cgraph, <stem>.bas.
 * A full path to any one of those files also works; the suffix is stripped.
 */

#include "CglOddWheel.hpp"
#include "CoinOddWheelSeparator.hpp"
#include "ClpSimplex.hpp"
#include "CoinStaticConflictGraph.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCuts.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/stat.h>
#include <vector>

static int g_failures = 0;

static void check(bool ok, const char *msg)
{
  if (!ok) {
    printf("FAIL: %s\n", msg);
    ++g_failures;
  }
}

/// Wall-clock, not CPU: separation cost as a user would feel it.
static double wallClock()
{
  return CoinGetTimeOfDay();
}

static bool fileExists(const std::string &path)
{
  struct stat st;
  return stat(path.c_str(), &st) == 0;
}

/**
 * Reduce any of the fixture's file names to the shared stem, so a caller can
 * pass whichever one tab-completion produced.
 */
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

/// Last path component, for the CSV name column.
static std::string baseName(const std::string &path)
{
  const size_t slash = path.rfind('/');
  return slash == std::string::npos ? path : path.substr(slash + 1);
}

/**
 * Read one integer key out of the fixture's `.meta`.
 * Returns `dflt` when the file or the key is absent, so a fixture written before
 * `.meta` carried that key still loads.
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
 * see cbcClqFixtureWriteMps. The row is redundant, so leaving it in place would
 * not change the LP, but it would change the row count the separator sees and
 * so the cut and row-count columns; drop it and the model is the captured one.
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
 * column CBC had fixed by bound tightening reads back continuous.
 *
 * CglOddWheel itself reads only the conflict graph, so integrality does not gate
 * it the way it gates the clique-row detector. It is restored anyway, because
 * refreshSolver() recomputes column types and because a fixture that loads
 * differently here than in bkclique-bench would not be the same experiment.
 *
 * Returns the number of columns re-marked, or -1 when no usable sidecar was
 * found. A sidecar whose column count disagrees with the loaded model is refused
 * rather than partly applied.
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

  int idx = 0, type = 0, restored = 0;
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
    }
  }
  fclose(fp);

  // getColType() caches, and derives Binary vs GeneralInteger from the bounds, so
  // it has to be recomputed after this or consumers would see the pre-restore view.
  si.getColType(true);
  return restored;
}

/// Everything the fixture determines, plus how long loading it took.
struct Fixture {
  OsiClpSolverInterface si;
  double warmStartTime = 0.0;
  /// Pivots the warm start needed. Zero is the expected value and the check that
  /// the captured basis actually took; see loadFixture below.
  int warmStartIters = 0;
  double cgraphTime = 0.0;
  bool paddedRowDropped = false;
  int restoredColTypes = -1;
  bool ok = false;
};

/**
 * Load the fixture and warm-start to the captured optimum.
 *
 * Two things are needed to actually *land* on the captured vertex, and getting
 * either wrong looks like success while silently changing the experiment -- the
 * odd cycles found are those violated by the current fractional point, so a
 * different optimal vertex means a different set of cuts.
 *
 * First, `readBasis` writes into ClpSimplex's own `status_` array, but OsiClp
 * caches a separate `CoinWarmStartBasis basis_` and both entry points overwrite
 * the model from it -- `resolve()` at OsiClpSolverInterface.cpp:1199, and
 * `initialSolve()` by presolving from scratch. `setWarmStart(NULL)` refreshes
 * `basis_` from the model (`basis_ = getBasis(modelPtr_)`, applying the slack
 * flip), which is what makes the file's basis survive into the solve.
 *
 * Second, `resolve()` rather than `initialSolve()`: presolve discards the basis.
 * With both in place an already-optimal fixture costs 0 iterations, which is the
 * cheap self-check that the warm start worked. `setPerturbation(50)` matches the
 * generator; perturbation left on moves the vertex even from a correct basis
 * (decomp2: 846 columns moved versus 803).
 */
static bool loadFixture(Fixture &f, const std::string &stem, bool rebuildCgraph,
  bool quiet)
{
  const std::string mps = fileExists(stem + ".mps.gz") ? stem + ".mps.gz" : stem + ".mps";
  const std::string bas = stem + ".bas";
  const std::string cgr = stem + ".cgraph";

  if (!fileExists(mps)) {
    fprintf(stderr, "ERROR: no problem file for stem %s\n", stem.c_str());
    return false;
  }

  ClpSimplex *lp = f.si.getModelPtr();
  lp->setLogLevel(0);
  if (f.si.readMps(mps.c_str())) {
    fprintf(stderr, "ERROR: failed to read %s\n", mps.c_str());
    return false;
  }

  // Before the basis: the basis was written from the captured model, so it has
  // one artificial fewer than the padded matrix and would not line up.
  f.paddedRowDropped = dropPadRow(f.si, stem, quiet);

  f.restoredColTypes = restoreColTypes(f.si, stem, quiet);

  bool haveBasis = false;
  if (fileExists(bas)) {
    if (lp->readBasis(bas.c_str()) < 0) {
      fprintf(stderr, "WARNING: failed to read basis %s; solving cold\n", bas.c_str());
    } else {
      haveBasis = true;
      // Push the model's freshly-read status into OsiClp's cached basis_, or the
      // solve below installs the stale cache over it.
      f.si.setWarmStart(NULL);
    }
  } else if (!quiet) {
    fprintf(stderr, "WARNING: no basis %s; solving cold\n", bas.c_str());
  }

  lp->setPerturbation(50);

  const double t0 = wallClock();
  if (haveBasis) {
    f.si.setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
    f.si.setHintParam(OsiDoDualInResolve, true, OsiHintDo);
    f.si.resolve();
  } else {
    f.si.setHintParam(OsiDoDualInInitial, true, OsiHintDo);
    f.si.initialSolve();
  }
  f.warmStartTime = wallClock() - t0;
  f.warmStartIters = f.si.getIterationCount();

  if (!f.si.isProvenOptimal()) {
    fprintf(stderr, "ERROR: LP not optimal after warm start (%s)\n", stem.c_str());
    return false;
  }
  // A correct warm start from an optimal basis costs no pivots. Anything else means
  // the fixture landed on a different vertex, so say so rather than quietly
  // measuring a different LP.
  if (haveBasis && f.warmStartIters > 0 && !quiet) {
    fprintf(stderr, "WARNING: %s: warm start took %d iterations; the captured basis "
                    "did not survive, so this is a different vertex\n",
      baseName(stem).c_str(), f.warmStartIters);
  }

  const double t1 = wallClock();
  if (rebuildCgraph) {
    f.si.checkCGraph(NULL);
  } else {
    if (!fileExists(cgr)) {
      fprintf(stderr, "ERROR: no conflict graph %s (use --rebuild-cgraph to build one)\n",
        cgr.c_str());
      return false;
    }
    CoinStaticConflictGraph *cg = CoinStaticConflictGraph::load(cgr.c_str());
    if (!cg) {
      fprintf(stderr, "ERROR: failed to load conflict graph %s\n", cgr.c_str());
      return false;
    }
    if (cg->size() != (size_t)f.si.getNumCols() * 2) {
      // CglOddWheel calls exit() on this, so catch it here with a message that
      // says which fixture is inconsistent.
      fprintf(stderr, "ERROR: graph/model mismatch for %s: graph %lu nodes, model %d cols\n",
        stem.c_str(), (unsigned long)cg->size(), f.si.getNumCols());
      delete cg;
      return false;
    }
    f.si.setCGraph(cg);
  }
  f.cgraphTime = wallClock() - t1;

  f.ok = true;
  return true;
}

#define OWF_CSV_HEADER                                                        \
  "name,nAct,nFrac,xMean,xMax,actArcs,asymArcs,actDegLt2,core2,bipNodes,"      \
  "triAny,triDomGlobal,triDom,certShort,certCore,skipAll,actMin,featTime"

/**
 * Cheap structural features of the *active* subgraph, to answer "could we have
 * known in advance that this call would find nothing?".
 *
 * Motivation: 0-cut fixtures burn the overwhelming majority of separation time,
 * and the counters say the rejecter is ohShort -- the shortest odd closed walk
 * through the node is a triangle -- not a violation test and not a failed
 * search. So the useful feature is not "how fractional is this node" but
 * "does this node's conflict neighbourhood already close a cheap triangle".
 *
 * TWO OF THESE ARE SOUND SKIP CERTIFICATES, not heuristics. A node they reject
 * provably cannot yield a cut, so gating on them cannot lose a valid cut:
 *
 *  - actDegLt2: fewer than 2 neighbours *among active nodes*. The separator's
 *    own gate tests degree in the whole graph (CoinOddWheelSeparator.cpp:175),
 *    which is weaker: a node can pass it and still lie on no cycle of the
 *    subgraph actually searched.
 *
 *  - bipNodes: the node sits in a 2-colourable component of the active
 *    subgraph. prepareGraph builds the bipartite double cover -- arcs
 *    (i1, n+i2) at :317 and their mirrors (n+i1, i2) at :460 -- so a path
 *    v' -> v'' is exactly an odd-length closed walk through v. Two-colour the
 *    component treating every arc as an undirected constraint; if that
 *    succeeds, every arc joins unlike colours, so every closed walk has even
 *    length and no such path exists. find() then returns nothing, path()
 *    yields oddSize 0, and the call is charged to oddHolesShort -- which is
 *    why "short" cannot be read as "found a triangle". Sound whether or not
 *    the adjacency is symmetric, since 2-colouring uses each arc as a
 *    constraint in both directions.
 *
 *  - triDom: arcs are pushed as (i1 -> icaCount_+i2, icaActivity_[i2]) at :317
 *    and :324, i.e. an arc carries the *destination* node's activity, so the
 *    weight of a closed odd walk is the sum of the activities of the nodes it
 *    enters -- 3 terms for a triangle, >= 5 for anything longer. Activities are
 *    1001 - 1000x in [1, 1001], strictly positive, so a longer walk cannot be
 *    cheaper per node. An odd closed walk v -> u1 -> ... -> u_{L-1} -> v with
 *    L >= 5 therefore weighs at least
 *
 *        acti(v) + 2*minNbrActi(v) + 2*minActi
 *
 *    because positions 1 and L-1 are both neighbours of v (two entries in the
 *    multiset even if they are the same node) and positions 2..L-2 are at least
 *    two further entries. If the lightest triangle through v beats that bound,
 *    the minimum-weight odd walk through v IS a triangle, find() must return
 *    one, and the separator discards it as too short. Skipping v is then free.
 *
 *    triDomGlobal is the same idea with the weaker bound "sum of the five
 *    smallest activities anywhere", kept because it shows how much the per-node
 *    neighbour term matters -- it fires on nothing at all.
 *
 * The triangle probe is deliberately incomplete (top-K lightest neighbours):
 * missing a triangle costs a skip we could have taken, never a cut. The gate is
 * one-sided by construction, which is what makes it safe to be approximate.
 *
 * core2 and the fractionality columns are the *hypotheses* being tested against
 * these, not proposals -- they are here to be compared, and may well be inert.
 */
static void printNodeFeatures(const CoinConflictGraph *cg, const double *xCols,
  int numCols, const std::string &name, bool header)
{
  const double t0 = wallClock();
  const size_t n = cg->size();

  // The doubled graph: node j is "x_j = 1", node j+numCols is its complement.
  std::vector< double > x(n, 0.0);
  for (int j = 0; j < numCols; ++j) {
    x[(size_t)j] = xCols[j];
    if ((size_t)j + (size_t)numCols < n)
      x[(size_t)j + (size_t)numCols] = 1.0 - xCols[j];
  }

  // Exactly fillActiveColumns' rule -- anything else would describe a
  // different call than the one being measured.
  std::vector< size_t > act;
  for (size_t j = 0; j < n; ++j) {
    if (cg->degree(j) < 2)
      continue;
    if (x[j] + 1e-6 <= 0.001)
      continue;
    act.push_back(j);
  }
  const size_t nAct = act.size();

  if (header)
    printf("%s\n", OWF_CSV_HEADER);
  // searchOddWheels() returns immediately on icaCount_ <= 4, so the separator
  // makes no shortest-path call at all and there is nothing here to describe.
  // Without this the six fixtures with nAct in {1,2,4} report a certificate
  // against zero calls and read as soundness failures.
  if (nAct <= 4) {
    printf("%s,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,%.3f\n", name.c_str(),
      wallClock() - t0);
    return;
  }

  std::vector< size_t > pos(n, (size_t)-1);
  for (size_t i = 0; i < nAct; ++i)
    pos[act[i]] = i;

  std::vector< double > acti(nAct);
  double xMax = 0.0, xSum = 0.0;
  size_t nFrac = 0;
  for (size_t i = 0; i < nAct; ++i) {
    const double v = x[act[i]];
    acti[i] = 1001.0 - 1000.0 * v;
    xSum += v;
    if (v > xMax)
      xMax = v;
    if (v >= 0.001 && v <= 0.999)
      ++nFrac;
  }

  // Lower bound on the weight of any odd closed walk of length >= 5.
  std::vector< double > byWeight(acti);
  std::sort(byWeight.begin(), byWeight.end());
  double w5min = 0.0;
  for (size_t k = 0; k < 5 && k < byWeight.size(); ++k)
    w5min += byWeight[k];

  // Adjacency restricted to active nodes.
  std::vector< size_t > temp(n);
  std::vector< char > iv(n, 0);
  std::vector< std::vector< size_t > > adj(nAct);
  size_t arcs = 0;
  for (size_t i = 0; i < nAct; ++i) {
    const std::pair< size_t, const size_t * > conf
      = cg->conflictingNodes(act[i], temp.data(), iv.data());
    for (size_t k = 0; k < conf.first; ++k) {
      const size_t p = pos[conf.second[k]];
      if (p == (size_t)-1 || p == i)
        continue;
      adj[i].push_back(p);
      ++arcs;
    }
  }

  // Sort each list so membership -- "is this arc really in the graph the search
  // walks?" -- is a binary search rather than an appeal to symmetry.
  for (size_t i = 0; i < nAct; ++i)
    std::sort(adj[i].begin(), adj[i].end());

  // Is the active adjacency symmetric? Every gate below would be simpler if it
  // were, and the O(n^2) branch of prepareGraph implicitly assumes it, so count
  // the arcs whose reverse is missing rather than assuming either way.
  size_t asymArcs = 0;
  for (size_t i = 0; i < nAct; ++i)
    for (size_t k = 0; k < adj[i].size(); ++k) {
      const std::vector< size_t > &back = adj[adj[i][k]];
      if (!std::binary_search(back.begin(), back.end(), i))
        ++asymArcs;
    }

  // Symmetrised adjacency. Both endpoints of a closed walk's first and last step
  // are neighbours of v in *this* relation whichever way the arcs point, so
  // using it keeps the degree, 2-core and lb5 bounds sound on a directed graph.
  std::vector< std::vector< size_t > > und(adj);
  for (size_t i = 0; i < nAct; ++i)
    for (size_t k = 0; k < adj[i].size(); ++k)
      und[adj[i][k]].push_back(i);
  for (size_t i = 0; i < nAct; ++i) {
    std::sort(und[i].begin(), und[i].end());
    und[i].erase(std::unique(und[i].begin(), und[i].end()), und[i].end());
  }

  size_t degLt2 = 0;
  for (size_t i = 0; i < nAct; ++i)
    if (und[i].size() < 2)
      ++degLt2;

  // 2-core by peeling the symmetrised graph. A node outside it lies on no simple
  // cycle, hence in no odd *hole*, so findOddHolesWithNode can never turn it into
  // a cut -- but note what it does NOT imply: the shortest path v' -> v'' may
  // still exist, because the doubled graph lets a walk traverse one edge twice.
  // v of degree 1 with neighbour a, a on a triangle a-b-c, gives the legal path
  // v' -> a'' -> b' -> c'' -> a' -> v'' (five distinct doubled nodes) of odd
  // length 5. So such a call lands in oddHolesRepeatedNode, not oddHolesShort --
  // which is exactly what made the first version of the skipAll <= ohShort check
  // fail on 31 fixtures. certCore is tallied separately for that reason.
  std::vector< size_t > deg(nAct);
  std::vector< char > dead(nAct, 0);
  std::vector< size_t > stack;
  for (size_t i = 0; i < nAct; ++i) {
    deg[i] = und[i].size();
    if (deg[i] < 2) {
      dead[i] = 1;
      stack.push_back(i);
    }
  }
  while (!stack.empty()) {
    const size_t v = stack.back();
    stack.pop_back();
    for (size_t k = 0; k < und[v].size(); ++k) {
      const size_t u = und[v][k];
      if (dead[u])
        continue;
      if (deg[u])
        --deg[u];
      if (deg[u] < 2) {
        dead[u] = 1;
        stack.push_back(u);
      }
    }
  }
  size_t core2 = 0;
  for (size_t i = 0; i < nAct; ++i)
    if (!dead[i])
      ++core2;

  // Two-colour each component of the symmetrised active subgraph. A node in a
  // 2-colourable component lies on no odd closed walk, so the search cannot
  // return anything for it.

  std::vector< signed char > colour(nAct, -1);
  std::vector< char > bip(nAct, 0);
  std::vector< size_t > comp, queue;
  for (size_t s = 0; s < nAct; ++s) {
    if (colour[s] >= 0)
      continue;
    comp.clear();
    queue.clear();
    colour[s] = 0;
    queue.push_back(s);
    bool twoColourable = true;
    for (size_t head = 0; head < queue.size(); ++head) {
      const size_t v = queue[head];
      comp.push_back(v);
      for (size_t k = 0; k < und[v].size(); ++k) {
        const size_t u = und[v][k];
        if (colour[u] < 0) {
          colour[u] = colour[v] ^ 1;
          queue.push_back(u);
        } else if (colour[u] == colour[v]) {
          twoColourable = false;
        }
      }
    }
    if (twoColourable)
      for (size_t k = 0; k < comp.size(); ++k)
        bip[comp[k]] = 1;
  }
  size_t bipNodes = 0;
  for (size_t i = 0; i < nAct; ++i)
    if (bip[i])
      ++bipNodes;

  // Triangle dominance. The closure test walks the arc lists prepareGraph
  // actually builds -- v -> a, a -> b, b -> v -- rather than asking
  // cg->conflicting(a, b), which is a claim about the conflict graph and not
  // about the directed graph the search walks. (Measured afterwards: asymArcs is
  // 0 on all 336 fixtures, so on this fixture set the two agree -- but the arc
  // form costs nothing and does not rest on that.)
  //
  // The bound is a bound on odd closed *walks*, not cycles, so it stays sound
  // under the repeated-edge shape described above certCore.
  const double actMin = byWeight[0];
  const size_t K = 16;
  size_t triAny = 0, triDom = 0, triDomGlobal = 0;
  size_t certShort = 0, certCore = 0;
  std::vector< std::pair< double, size_t > > cand;
  for (size_t i = 0; i < nAct; ++i) {
    if (bip[i]) {
      ++certShort;
      continue;
    }
    // Lower bound on any odd closed walk of length >= 5 through this node. Its
    // first and last steps both touch symmetrised neighbours of i (two entries
    // in the multiset even when they are the same node), and at least two
    // further entries lie somewhere in the active set.
    double minNbr = acti[und[i][0]];
    for (size_t k = 1; k < und[i].size(); ++k)
      if (acti[und[i][k]] < minNbr)
        minNbr = acti[und[i][k]];
    const double lb5 = acti[i] + 2.0 * minNbr + 2.0 * actMin;

    cand.clear();
    cand.reserve(adj[i].size());
    for (size_t k = 0; k < adj[i].size(); ++k)
      cand.push_back(std::make_pair(acti[adj[i][k]], adj[i][k]));
    std::sort(cand.begin(), cand.end());
    if (cand.size() > K)
      cand.resize(K);

    double best = 0.0;
    bool found = false;
    for (size_t a = 0; a < cand.size() && !(found && best < lb5); ++a) {
      const size_t va = cand[a].second;
      for (size_t k = 0; k < adj[va].size(); ++k) {
        const size_t vb = adj[va][k];
        if (vb == i)
          continue;
        if (!std::binary_search(adj[vb].begin(), adj[vb].end(), i))
          continue; // the closing arc vb -> i is absent: not a walk back to i
        const double w = acti[i] + cand[a].first + acti[vb];
        if (!found || w < best) {
          best = w;
          found = true;
        }
      }
    }
    if (found) {
      ++triAny;
      if (best < lb5)
        ++triDom;
      if (best < w5min)
        ++triDomGlobal;
    }
    if (found && best < lb5)
      ++certShort;
    else if (dead[i])
      ++certCore;
  }
  const size_t skipAll = certShort + certCore;

  printf(
    "%s,%lu,%lu,%.6f,%.6f,%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,%.3f,%.3f\n",
    name.c_str(), (unsigned long)nAct, (unsigned long)nFrac,
    xSum / (double)nAct, xMax, (unsigned long)arcs, (unsigned long)asymArcs,
    (unsigned long)degLt2, (unsigned long)core2, (unsigned long)bipNodes,
    (unsigned long)triAny, (unsigned long)triDomGlobal, (unsigned long)triDom,
    (unsigned long)certShort, (unsigned long)certCore, (unsigned long)skipAll,
    actMin, wallClock() - t0);
}

/* One row per active node for the --node-outcomes dump. The label comes from the
 * separator itself, not from anything reconstructed here. */
#define OWO_CSV_HEADER                                                        \
  "name,node,col,isComp,x,degOut,degIn,degSym,"                               \
  "xNbrMin,xNbrMax,xNbrMean,nNbrMid,xTop4Nbr,viol5,"                          \
  "compSize,bip,core2,triAny,triBestW,lb5,triDom,"                            \
  "outcome"

/**
 * Dump one labelled row per active node, so a classifier can be fitted on the
 * ~640k individual shortest-path calls instead of on 336 fixture aggregates.
 *
 * Why per node and not per fixture: 84.1% of all calls land in oddHolesShort and
 * 0.5% are kept, so a fixture-level counter tells you almost nothing about which
 * *nodes* were the waste. The label here is written by
 * CoinOddWheelSeparator::findOddHolesWithNode() at each of its five exits (see
 * setRecordNodeOutcomes), so it is the real outcome and not an inference from
 * aggregate counters.
 *
 * The gate is forced OFF and maxSeconds to 0: with the gate on, every skipped
 * node reports OUTCOME_NOT_CALLED and the labels go missing exactly where the
 * gate fires, i.e. on the rows the study is about.
 *
 * Two honest caveats about anything fitted on this.
 *
 *  - Three of the features here (bip, core2, triDom) are the *certificates*
 *    already implemented in buildFutilityGate(). They are emitted so a tree can
 *    be measured against them rather than credited for rediscovering them; any
 *    claimed improvement must come from the other columns.
 *  - A tree is a heuristic, not a certificate. Skipping on a learned rule can
 *    lose a cut, which the three certificates provably cannot. So the only
 *    defensible use of a tree here is either (a) to find a *new* rule that can
 *    then be proved, or (b) to bound how much residual waste is left for a
 *    certificate to reach at all. Do not wire a fitted threshold into the gate.
 *
 * viol5 is the arithmetic the user's fractionality question turns on: an odd hole
 * of size k is violated only when sum(x) > (k-1)/2, so a 5-hole through this node
 * needs x + (its four heaviest neighbours) > 2.0 + MIN_VIOL. It is an optimistic
 * bound -- the four heaviest neighbours need not form a hole -- so viol5 <= 0 is
 * a genuine "no 5-hole through here can be violated" statement, while viol5 > 0
 * says nothing. It is *not* sound as a skip on its own, because a longer hole can
 * be violated where a 5-hole cannot.
 */
static void printNodeOutcomes(const CoinConflictGraph *cg, const double *xCols,
  const double *rcCols, int numCols, size_t extMethod, const std::string &name,
  bool header)
{
  const size_t n = cg->size();

  std::vector< double > x(n, 0.0), rc(n, 0.0);
  for (int j = 0; j < numCols; ++j) {
    x[(size_t)j] = xCols[j];
    rc[(size_t)j] = rcCols[j];
    if ((size_t)j + (size_t)numCols < n) {
      x[(size_t)j + (size_t)numCols] = 1.0 - xCols[j];
      rc[(size_t)j + (size_t)numCols] = -rcCols[j];
    }
  }

  // Run the separator first, so the labels exist before any feature work. Same
  // rule as fillActiveColumns is replicated below to index them; the two agree
  // because both walk j ascending over the doubled graph with the same filter.
  CoinOddWheelSeparator sep(cg, x.data(), rc.data(), extMethod);
  sep.setUseFutilityGate(false);
  sep.setRecordNodeOutcomes(true);
  sep.searchOddWheels();
  const std::vector< unsigned char > &lab = sep.nodeOutcomes();

  std::vector< size_t > act;
  for (size_t j = 0; j < n; ++j) {
    if (cg->degree(j) < 2)
      continue;
    if (x[j] + 1e-6 <= 0.001)
      continue;
    act.push_back(j);
  }
  const size_t nAct = act.size();

  if (header)
    printf("%s\n", OWO_CSV_HEADER);

  // searchOddWheels() returns immediately on icaCount_ <= 4, so no call was made
  // from any of these nodes and there is no outcome to report.
  if (nAct <= 4)
    return;
  if (lab.size() != nAct) {
    fprintf(stderr, "ERROR: %s: %lu labels for %lu active nodes -- the active "
                    "rule here and fillActiveColumns' have diverged\n",
      name.c_str(), (unsigned long)lab.size(), (unsigned long)nAct);
    return;
  }

  std::vector< size_t > pos(n, (size_t)-1);
  for (size_t i = 0; i < nAct; ++i)
    pos[act[i]] = i;

  std::vector< double > acti(nAct), xa(nAct);
  double actMin = 1e100;
  for (size_t i = 0; i < nAct; ++i) {
    xa[i] = x[act[i]];
    acti[i] = 1001.0 - 1000.0 * xa[i];
    if (acti[i] < actMin)
      actMin = acti[i];
  }

  std::vector< size_t > temp(n);
  std::vector< char > iv(n, 0);
  std::vector< std::vector< size_t > > adj(nAct);
  for (size_t i = 0; i < nAct; ++i) {
    const std::pair< size_t, const size_t * > conf
      = cg->conflictingNodes(act[i], temp.data(), iv.data());
    for (size_t k = 0; k < conf.first; ++k) {
      const size_t p = pos[conf.second[k]];
      if (p == (size_t)-1 || p == i)
        continue;
      adj[i].push_back(p);
    }
    std::sort(adj[i].begin(), adj[i].end());
  }

  std::vector< size_t > degIn(nAct, 0);
  for (size_t i = 0; i < nAct; ++i)
    for (size_t k = 0; k < adj[i].size(); ++k)
      ++degIn[adj[i][k]];

  std::vector< std::vector< size_t > > und(adj);
  for (size_t i = 0; i < nAct; ++i)
    for (size_t k = 0; k < adj[i].size(); ++k)
      und[adj[i][k]].push_back(i);
  for (size_t i = 0; i < nAct; ++i) {
    std::sort(und[i].begin(), und[i].end());
    und[i].erase(std::unique(und[i].begin(), und[i].end()), und[i].end());
  }

  // 2-core of the symmetrised graph.
  std::vector< size_t > deg(nAct);
  std::vector< char > dead(nAct, 0);
  std::vector< size_t > stack;
  for (size_t i = 0; i < nAct; ++i) {
    deg[i] = und[i].size();
    if (deg[i] < 2) {
      dead[i] = 1;
      stack.push_back(i);
    }
  }
  while (!stack.empty()) {
    const size_t v = stack.back();
    stack.pop_back();
    for (size_t k = 0; k < und[v].size(); ++k) {
      const size_t u = und[v][k];
      if (dead[u])
        continue;
      if (deg[u])
        --deg[u];
      if (deg[u] < 2) {
        dead[u] = 1;
        stack.push_back(u);
      }
    }
  }

  // Components and 2-colourability of the symmetrised graph. compSize is carried
  // per node because "how big is the piece I am in" is a plausible feature that
  // the fixture-level sweep could not express at all.
  std::vector< signed char > colour(nAct, -1);
  std::vector< char > bip(nAct, 0);
  std::vector< size_t > csize(nAct, 0);
  std::vector< size_t > comp, queue;
  for (size_t s = 0; s < nAct; ++s) {
    if (colour[s] >= 0)
      continue;
    comp.clear();
    queue.clear();
    colour[s] = 0;
    queue.push_back(s);
    bool twoColourable = true;
    for (size_t head = 0; head < queue.size(); ++head) {
      const size_t v = queue[head];
      comp.push_back(v);
      for (size_t k = 0; k < und[v].size(); ++k) {
        const size_t u = und[v][k];
        if (colour[u] < 0) {
          colour[u] = colour[v] ^ 1;
          queue.push_back(u);
        } else if (colour[u] == colour[v]) {
          twoColourable = false;
        }
      }
    }
    for (size_t k = 0; k < comp.size(); ++k) {
      csize[comp[k]] = comp.size();
      if (twoColourable)
        bip[comp[k]] = 1;
    }
  }

  const size_t K = 16;
  std::vector< std::pair< double, size_t > > cand;
  std::vector< double > nbrX;
  for (size_t i = 0; i < nAct; ++i) {
    double xMin = 1e100, xMax = -1e100, xSum = 0.0;
    size_t nMid = 0;
    nbrX.clear();
    for (size_t k = 0; k < und[i].size(); ++k) {
      const double v = xa[und[i][k]];
      nbrX.push_back(v);
      xSum += v;
      if (v < xMin)
        xMin = v;
      if (v > xMax)
        xMax = v;
      if (v >= 0.2 && v <= 0.8)
        ++nMid;
    }
    const size_t d = und[i].size();
    if (!d) {
      xMin = xMax = 0.0;
    }

    // Sum of the four heaviest neighbours: the optimistic ceiling on sum(x) over
    // a 5-hole through i. See viol5 in the comment above.
    std::sort(nbrX.begin(), nbrX.end(), std::greater< double >());
    double top4 = 0.0;
    for (size_t k = 0; k < 4 && k < nbrX.size(); ++k)
      top4 += nbrX[k];
    const double viol5 = xa[i] + top4 - 2.0;

    // Same bound and the same probe as buildFutilityGate's certificate 2, with
    // in- and out-minima kept separate for the reason recorded there.
    double minOut = 1e100, minIn = 1e100;
    for (size_t k = 0; k < adj[i].size(); ++k)
      if (acti[adj[i][k]] < minOut)
        minOut = acti[adj[i][k]];
    for (size_t k = 0; k < und[i].size(); ++k) {
      const size_t u = und[i][k];
      if (std::binary_search(adj[u].begin(), adj[u].end(), i)
        && acti[u] < minIn)
        minIn = acti[u];
    }
    const double lb5 = (minOut < 1e99 && minIn < 1e99)
      ? acti[i] + minOut + minIn + 2.0 * actMin
      : 0.0;

    cand.clear();
    for (size_t k = 0; k < adj[i].size(); ++k)
      cand.push_back(std::make_pair(acti[adj[i][k]], adj[i][k]));
    std::sort(cand.begin(), cand.end());
    if (cand.size() > K)
      cand.resize(K);
    double best = 0.0;
    bool found = false;
    for (size_t a = 0; a < cand.size() && !(found && best < lb5); ++a) {
      const size_t va = cand[a].second;
      for (size_t k = 0; k < adj[va].size(); ++k) {
        const size_t vb = adj[va][k];
        if (vb == i || vb == va)
          continue;
        if (!std::binary_search(adj[vb].begin(), adj[vb].end(), i))
          continue;
        const double w = acti[i] + cand[a].first + acti[vb];
        if (!found || w < best) {
          best = w;
          found = true;
        }
      }
    }

    printf("%s,%lu,%lu,%d,%.6f,%lu,%lu,%lu,"
           "%.6f,%.6f,%.6f,%lu,%.6f,%.6f,"
           "%lu,%d,%d,%d,%.3f,%.3f,%d,"
           "%d\n",
      name.c_str(), (unsigned long)i, (unsigned long)act[i],
      act[i] >= (size_t)numCols ? 1 : 0, xa[i],
      (unsigned long)adj[i].size(), (unsigned long)degIn[i], (unsigned long)d,
      xMin, xMax, d ? xSum / (double)d : 0.0, (unsigned long)nMid, top4, viol5,
      (unsigned long)csize[i], (int)bip[i], dead[i] ? 0 : 1, found ? 1 : 0,
      found ? best : 0.0, lb5, (found && best < lb5) ? 1 : 0,
      (int)lab[i]);
  }
}

/**
 * Round-trip the graph through save()/load() and prove the reconstruction is
 * identical, field by field. "Close" is not good enough here: the whole point of
 * serializing is that a fixture reproduces what the separator saw exactly,
 * including the deliberately *approximate* degrees, which a rebuild would
 * recompute differently -- and degree gates which nodes are active at all
 * (fillActiveColumns skips degree < 2) as well as which hole node is picked as
 * the wheel-center seed.
 */
static int selfTest(const std::string &stem)
{
  const std::string cgr = stem + ".cgraph";
  if (!fileExists(cgr)) {
    fprintf(stderr, "ERROR: no conflict graph %s\n", cgr.c_str());
    return 1;
  }

  CoinStaticConflictGraph *a = CoinStaticConflictGraph::load(cgr.c_str());
  if (!a) {
    printf("FAIL: could not load %s\n", cgr.c_str());
    return 1;
  }

  const std::string t1 = "/tmp/oddwheel-selftest-1.cgraph";
  const std::string t2 = "/tmp/oddwheel-selftest-2.cgraph";
  check(a->save(t1.c_str()) == 0, "save of loaded graph");

  CoinStaticConflictGraph *b = CoinStaticConflictGraph::load(t1.c_str());
  if (!b) {
    printf("FAIL: could not re-load saved graph\n");
    delete a;
    return 1;
  }
  check(b->save(t2.c_str()) == 0, "re-save of re-loaded graph");

  // Byte identity of save(load(save(x))) vs save(x): catches any field that is
  // written but silently not read back, which a getter-by-getter comparison
  // would miss if the getter does not expose it.
  {
    FILE *f1 = fopen(t1.c_str(), "rb");
    FILE *f2 = fopen(t2.c_str(), "rb");
    bool same = f1 && f2;
    long n = 0;
    while (same) {
      char b1[65536], b2[65536];
      const size_t r1 = fread(b1, 1, sizeof(b1), f1);
      const size_t r2 = fread(b2, 1, sizeof(b2), f2);
      if (r1 != r2 || memcmp(b1, b2, r1) != 0) {
        same = false;
        break;
      }
      n += (long)r1;
      if (r1 == 0)
        break;
    }
    if (f1)
      fclose(f1);
    if (f2)
      fclose(f2);
    check(same, "save/load/save is byte-identical");
    if (same)
      printf("  round-trip stable over %ld bytes\n", n);
  }

  check(a->size() == b->size(), "size preserved");
  check(a->density() == b->density(), "density preserved");
  check(a->minDegree() == b->minDegree(), "minDegree preserved");
  check(a->maxDegree() == b->maxDegree(), "maxDegree preserved");
  check(a->nCliques() == b->nCliques(), "nCliques preserved");
  check(a->nTotalDirectConflicts() == b->nTotalDirectConflicts(),
    "total direct conflicts preserved");
  check(a->nTotalCliqueElements() == b->nTotalCliqueElements(),
    "total clique elements preserved");

  if (a->size() == b->size()) {
    size_t badDeg = 0, badMDeg = 0, badConf = 0, badNodeClq = 0;
    for (size_t i = 0; i < a->size(); ++i) {
      if (a->degree(i) != b->degree(i))
        ++badDeg;
      if (a->modifiedDegree(i) != b->modifiedDegree(i))
        ++badMDeg;
      if (a->nNodeCliques(i) != b->nNodeCliques(i)) {
        ++badNodeClq;
        continue;
      }
      for (size_t k = 0; k < a->nNodeCliques(i); ++k)
        if (a->nodeCliques(i)[k] != b->nodeCliques(i)[k])
          ++badNodeClq;
      if (a->nDirectConflicts(i) != b->nDirectConflicts(i)) {
        ++badConf;
        continue;
      }
      for (size_t k = 0; k < a->nDirectConflicts(i); ++k)
        if (a->directConflicts(i)[k] != b->directConflicts(i)[k])
          ++badConf;
    }
    check(badDeg == 0, "per-node degree preserved");
    check(badMDeg == 0, "per-node modified degree preserved");
    check(badConf == 0, "per-node direct conflicts preserved");
    check(badNodeClq == 0, "per-node clique index rebuilt correctly");
  }

  if (a->nCliques() == b->nCliques()) {
    size_t badClq = 0;
    for (size_t c = 0; c < a->nCliques(); ++c) {
      if (a->cliqueSize(c) != b->cliqueSize(c)) {
        ++badClq;
        continue;
      }
      for (size_t k = 0; k < a->cliqueSize(c); ++k)
        if (a->cliqueElements(c)[k] != b->cliqueElements(c)[k])
          ++badClq;
    }
    check(badClq == 0, "clique elements preserved");
  }

  printf("%s: size=%lu directConf=%lu cliques=%lu clqElems=%lu density=%.10f\n",
    baseName(stem).c_str(), (unsigned long)a->size(),
    (unsigned long)a->nTotalDirectConflicts(), (unsigned long)a->nCliques(),
    (unsigned long)a->nTotalCliqueElements(), a->density());

  delete a;
  delete b;

  printf(g_failures ? "self-test: %d FAILURE(S)\n" : "self-test: all checks passed\n",
    g_failures);
  return g_failures ? 1 : 0;
}

static void usage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s <fixture-stem> [options]\n"
    "       %s --self-test <fixture-stem>\n"
    "\n"
    "<fixture-stem> is dir/name.tag; .mps.gz/.cgraph/.bas are appended.\n"
    "Passing any one of those files also works.\n"
    "\n"
    "Options:\n"
    "  --rounds=N          separation rounds, LP re-solved between them (default 4)\n"
    "  --ext-method=N      wheel-center lifting: 0 = none, 1 = one variable,\n"
    "                      2 = a clique (default 2, CBC's). 0 isolates the cost\n"
    "                      of lifting, since searchOddWheels() skips the whole\n"
    "                      lifting loop when it is 0.\n"
    "  --max-seconds=F     separator wall-clock budget (default 0 = none; a\n"
    "                      nonzero value makes the run load-dependent)\n"
    "  --rebuild-cgraph    rebuild the graph from the matrix instead of loading\n"
    "                      the captured one (not faithful; for comparison only)\n"
    "  --stage-times       also print a human-readable stage breakdown to stderr\n"
    "  --verify-prepare    build the auxiliary graph's arcs both ways and check\n"
    "                      they agree; prepMismatch must come out 0. Diagnostic\n"
    "                      only -- it roughly doubles graph preparation time\n"
    "  --check-validity    certify every odd wheel against the conflict graph\n"
    "                      before it becomes a row cut. certBad* must all be 0.\n"
    "                      Needs no reference solution, so unlike the row-cut\n"
    "                      debugger it works on every fixture\n"
    "  --no-gate           disable the separator's futility gate (default on).\n"
    "                      Every field except spFindCalls, ohShort, ohRepeated,\n"
    "                      the gate* counters and the times must be identical\n"
    "                      with and without it -- that is the no-cut-lost check.\n"
    "  --node-features     one row per fixture describing the active subgraph\n"
    "                      (bipartite / triangle-dominant / 2-core shares).\n"
    "                      Emits a different CSV schema and exits\n"
    "  --node-outcomes     one row per *active node*: features plus the outcome\n"
    "                      the separator actually reached for it (1 short,\n"
    "                      2 repeated, 3 not violated, 4 duplicate, 5 kept).\n"
    "                      Forces the gate off, or the rows the study is about\n"
    "                      are the ones with no label. Different CSV schema\n"
    "  --header            print the CSV header line and exit\n"
    "  --csv-header        print the CSV header before the data line\n"
    "  --quiet             suppress warnings\n",
    prog, prog);
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
 */
static bool boundMoved(double objStart, double objImprove)
{
  const double scale = fabs(objStart) > 1.0 ? fabs(objStart) : 1.0;
  return objImprove > 1.0e-9 * scale && objImprove > 1.0e-9;
}

static const char *CSV_HEADER
  = "name,extMethod,rounds,rowsAdded,totalCuts,totalViol,maxViol,avgCutLen,"
    "sepTime,warmStartTime,warmStartIters,cgraphTime,resolveTime,resolveIters,"
    "objStart,objEnd,objImprove,objImproveRel,boundMoved,"
    "cgNodes,cgDirectConf,cgCliques,cgDensity,restoredInt,"
    "icaCount,arcs,spFindCalls,oddHoles,ohShort,ohRepeated,ohNotViol,ohDuplicate,"
    "wheelCenters,wcElements,avgWcSize,cutsBeforePool,cutsDupIdx,cutsZeroCoefs,"
    "cutsEmpty,cutsAfterPool,"
    "timeLimitHit,tSetup,tSeparator,tActive,tPrepArcs,tPrepRev,tPrepSpf,tSearch,"
    "tWheelCenter,tCutPool,cutsPerRound,violPerRound,objImprovePerRound,"
    "prepMethod,prepWalkCost,prepUnsorted,prepVerifyArcs,prepMismatch,"
    "prepWalkOnly,prepPairOnly,"
    "certChecked,certBadCycle,certBadCenterAdj,certBadCenterClq,certBadAlpha,"
    "certBadTranslate,certComplCycle,certComplCenter,certComplPair,"
    "certCenterOnComplCycle,certComplAtLeastK,"
    // searchWheelCenter() filter attribution; appended so earlier column
    // positions are unchanged and older sweep CSVs stay comparable.
    "wcCalls,wcPool,wcRejInCycle,wcRejDegree,wcRejAdjacency,wcRejCost,"
    "wcCandidates,wcCliqueDropped,"
    // buildFutilityGate() attribution. gateSkipped is shortest-path calls proved
    // unable to yield a cut and therefore not made, so it moves spFindCalls,
    // ohShort, ohRepeated and tSearch and nothing else -- which is exactly what
    // --no-gate is for checking.
    //
    // gateBipartite/gateNoCycle/gateBlockOnly are one decision split three ways:
    // the skip is made by the block certificate ("some biconnected block
    // containing this node is non-bipartite" decides "lies on an odd cycle"
    // exactly), and the count goes to whichever weaker certificate would also
    // have caught the node. So gateBlockOnly is the column that prices what the
    // block certificate adds over the bipartite-component and 2-core tests it
    // subsumes; if it is 0 everywhere, the exact certificate bought nothing.
    "gateSkipped,gateBipartite,gateNoCycle,gateBlockOnly,gateTriangle,tGate";

/**
 * Sums of CglOddWheel::stats() over the rounds.
 *
 * Kept separate from the per-call struct rather than accumulated into it: the
 * generator is rebuilt every round (see the loop) so its own stats are always
 * one call's worth, which is what makes a single round's numbers readable in
 * isolation when --rounds=1.
 */
struct Totals {
  size_t icaCount = 0, arcs = 0, spFindCalls = 0;
  // prepareGraph() picks its method per call, so a run where the rounds disagree
  // is worth seeing rather than averaging away: 3 means both methods were used.
  size_t prepMethod = 0, prepWalkCost = 0, prepVerifyArcs = 0, prepMismatch = 0;
  size_t prepUnsorted = 0, prepWalkOnly = 0, prepPairOnly = 0;
  // setCheckValidity(): certChecked is coverage, the five certBad* are the
  // verdict and must be 0, the three certCompl* measure how much of the cut
  // set actually depends on the complemented half of the doubled graph.
  size_t certChecked = 0, certBadCycle = 0, certBadCenterAdj = 0;
  size_t certBadCenterClq = 0, certBadAlpha = 0, certBadTranslate = 0;
  size_t certComplCycle = 0, certComplCenter = 0, certComplPair = 0;
  size_t certCenterOnComplCycle = 0, certComplAtLeastK = 0;
  size_t oddHoles = 0, ohShort = 0, ohRepeated = 0, ohNotViol = 0, ohDuplicate = 0;
  size_t wheelCenters = 0, wcElements = 0;
  size_t wcCalls = 0, wcPool = 0, wcRejInCycle = 0, wcRejDegree = 0;
  size_t wcRejAdjacency = 0, wcRejCost = 0, wcCandidates = 0, wcCliqueDropped = 0;
  size_t gateSkipped = 0, gateBipartite = 0, gateTriangle = 0, gateNoCycle = 0;
  size_t gateBlockOnly = 0;
  size_t cutsBeforePool = 0, cutsDupIdx = 0, cutsAfterPool = 0;
  size_t cutsZeroCoefs = 0, cutsEmpty = 0;
  bool timeLimitHit = false;
  double tSetup = 0.0, tSeparator = 0.0, tActive = 0.0, tPrepArcs = 0.0;
  double tPrepRev = 0.0, tPrepSpf = 0.0, tSearch = 0.0, tWheelCenter = 0.0;
  double tCutPool = 0.0, tGate = 0.0;

  void add(const CglOddWheel::Stats &s)
  {
    // icaCount is a size, not a count of events, so the max over the rounds is
    // the meaningful summary: rounds after the first separate a *larger* LP but
    // usually a less fractional one, and summing would report a node count no
    // single call ever saw.
    if (s.sep.activeColumns > icaCount)
      icaCount = s.sep.activeColumns;
    if (s.sep.arcs > arcs)
      arcs = s.sep.arcs;
    prepMethod |= s.sep.prepareMethod;
    if (s.sep.prepareWalkCost > prepWalkCost)
      prepWalkCost = s.sep.prepareWalkCost;
    prepVerifyArcs += s.sep.prepareVerifyArcs;
    prepMismatch += s.sep.prepareMismatches;
    prepUnsorted += s.sep.prepareUnsorted;
    prepWalkOnly += s.sep.prepareWalkOnly;
    prepPairOnly += s.sep.preparePairOnly;

    certChecked += s.certChecked;
    certBadCycle += s.certBadCycle;
    certBadCenterAdj += s.certBadCenterAdj;
    certBadCenterClq += s.certBadCenterClq;
    certBadAlpha += s.certBadAlpha;
    certBadTranslate += s.certBadTranslate;
    certComplCycle += s.certComplCycle;
    certComplCenter += s.certComplCenter;
    certComplPair += s.certComplPair;
    certCenterOnComplCycle += s.certCenterOnComplCycle;
    certComplAtLeastK += s.certComplAtLeastK;

    spFindCalls += s.sep.spFindCalls;
    oddHoles += s.sep.oddHolesFound;
    ohShort += s.sep.oddHolesShort;
    ohRepeated += s.sep.oddHolesRepeatedNode;
    ohNotViol += s.sep.oddHolesNotViolated;
    ohDuplicate += s.sep.oddHolesDuplicate;
    wheelCenters += s.sep.wheelCenters;
    wcElements += s.sep.wheelCenterElements;
    wcCalls += s.sep.wcCalls;
    wcPool += s.sep.wcPool;
    wcRejInCycle += s.sep.wcRejInCycle;
    wcRejDegree += s.sep.wcRejDegree;
    wcRejAdjacency += s.sep.wcRejAdjacency;
    wcRejCost += s.sep.wcRejCost;
    wcCandidates += s.sep.wcCandidates;
    wcCliqueDropped += s.sep.wcCliqueDropped;
    gateSkipped += s.sep.gateSkipped;
    gateBipartite += s.sep.gateBipartite;
    gateTriangle += s.sep.gateTriangle;
    gateNoCycle += s.sep.gateNoCycle;
    gateBlockOnly += s.sep.gateBlockOnly;
    timeLimitHit = timeLimitHit || s.sep.timeLimitReached;

    cutsBeforePool += s.cutsBeforePool;
    cutsDupIdx += s.cutsDuplicatedIdx;
    cutsZeroCoefs += s.cutsZeroCoefs;
    cutsEmpty += s.cutsEmpty;
    cutsAfterPool += s.cutsAfterPool;

    tSetup += s.tSetup;
    tSeparator += s.tSeparator;
    tActive += s.sep.tActiveColumns;
    tPrepArcs += s.sep.tPrepareArcs;
    tPrepRev += s.sep.tPrepareReverse;
    tPrepSpf += s.sep.tPrepareShortestPath;
    tSearch += s.sep.tSearch;
    tWheelCenter += s.sep.tWheelCenter;
    tCutPool += s.tCutPool;
    tGate += s.sep.tGate;
  }
};

/**
 * The stage breakdown, on stderr so it never contaminates the CSV line.
 *
 * tSearch is reported as one figure and deliberately not split per
 * CoinShortestPath::find() call: two clock reads times icaCount_ calls is
 * hundreds of thousands of reads on a large fixture, which would perturb the
 * very loop being measured. The Dijkstra share comes out by difference instead.
 */
static void printStageTimes(const Totals &t, double totalSepTime)
{
  const double stages = t.tSetup + t.tActive + t.tPrepArcs + t.tPrepRev + t.tPrepSpf
    + t.tGate + t.tSearch + t.tWheelCenter + t.tCutPool;
  const double pct = totalSepTime > 0.0 ? 100.0 / totalSepTime : 0.0;

  fprintf(stderr, "\n  stage                       seconds     %% of sepTime\n");
  fprintf(stderr, "  doubled x/rc setup       %10.4f   %8.2f\n", t.tSetup, t.tSetup * pct);
  fprintf(stderr, "  fillActiveColumns        %10.4f   %8.2f\n", t.tActive, t.tActive * pct);
  fprintf(stderr, "  prepareGraph: conflicts  %10.4f   %8.2f\n", t.tPrepArcs, t.tPrepArcs * pct);
  fprintf(stderr, "  prepareGraph: reverse    %10.4f   %8.2f\n", t.tPrepRev, t.tPrepRev * pct);
  fprintf(stderr, "  prepareGraph: sp ctor    %10.4f   %8.2f\n", t.tPrepSpf, t.tPrepSpf * pct);
  fprintf(stderr, "  odd-hole search          %10.4f   %8.2f\n", t.tSearch, t.tSearch * pct);
  fprintf(stderr, "  wheel-center lifting     %10.4f   %8.2f\n", t.tWheelCenter,
    t.tWheelCenter * pct);
  fprintf(stderr, "  futility gate            %10.4f   %8.2f\n", t.tGate, t.tGate * pct);
  fprintf(stderr, "  cut pool + insertion     %10.4f   %8.2f\n", t.tCutPool, t.tCutPool * pct);
  fprintf(stderr, "  ------------------------------------------------\n");
  fprintf(stderr, "  accounted                %10.4f   %8.2f\n", stages, stages * pct);
  fprintf(stderr, "  sepTime (wall)           %10.4f   %8.2f\n", totalSepTime, 100.0);
  fprintf(stderr, "\n  shortest-path calls %lu, odd holes kept %lu of %lu paths"
                  " (%lu short, %lu repeated node, %lu not violated, %lu duplicate)\n",
    (unsigned long)t.spFindCalls, (unsigned long)t.oddHoles,
    (unsigned long)t.spFindCalls, (unsigned long)t.ohShort,
    (unsigned long)t.ohRepeated, (unsigned long)t.ohNotViol,
    (unsigned long)t.ohDuplicate);
  fprintf(stderr, "  wheel centers %lu with %lu elements; cuts %lu pooled,"
                  " %lu merged (%lu coefs cancelled, %lu emptied), %lu survived\n\n",
    (unsigned long)t.wheelCenters, (unsigned long)t.wcElements,
    (unsigned long)t.cutsBeforePool, (unsigned long)t.cutsDupIdx,
    (unsigned long)t.cutsZeroCoefs, (unsigned long)t.cutsEmpty,
    (unsigned long)t.cutsAfterPool);

  // Every pooled wheel-centre candidate is dropped by exactly one of the four
  // filters or survives, so this identity must hold. It is what makes the
  // per-filter percentages a partition rather than four unrelated tallies --
  // without it, a candidate rejected by two filters at once (or a counter
  // placed on the wrong side of a `continue`) would silently skew the split.
  const size_t wcAccounted = t.wcRejInCycle + t.wcRejDegree + t.wcRejAdjacency
                           + t.wcRejCost + t.wcCandidates;
  if (wcAccounted != t.wcPool)
    fprintf(stderr, "  ** wheel-centre attribution does not partition the pool:"
                    " %lu accounted vs %lu pooled (delta %ld)\n",
      (unsigned long)wcAccounted, (unsigned long)t.wcPool,
      (long)wcAccounted - (long)t.wcPool);
  fprintf(stderr, "  wheel-centre pool %lu = %lu in cycle + %lu low degree"
                  " + %lu not adjacent to all of C + %lu below cost gate"
                  " + %lu candidates (%lu dropped by the clique test)\n\n",
    (unsigned long)t.wcPool, (unsigned long)t.wcRejInCycle,
    (unsigned long)t.wcRejDegree, (unsigned long)t.wcRejAdjacency,
    (unsigned long)t.wcRejCost, (unsigned long)t.wcCandidates,
    (unsigned long)t.wcCliqueDropped);
}

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

  bool doSelfTest = false;
  bool doNodeFeatures = false;
  bool rebuildCgraph = false;
  bool csvHeader = false;
  bool quiet = false;
  bool stageTimes = false;
  bool verifyPrepare = false;
  bool checkValidity = false;
  bool useGate = true;
  bool doNodeOutcomes = false;
  int maxRounds = 4;
  size_t extMethod = 2;
  double maxSeconds = 0.0;
  const char *stemArg = NULL;

  for (int i = 1; i < argc; ++i) {
    const char *a = argv[i];
    if (strcmp(a, "--self-test") == 0) {
      doSelfTest = true;
    } else if (strcmp(a, "--node-features") == 0) {
      doNodeFeatures = true;
    } else if (strcmp(a, "--node-outcomes") == 0) {
      doNodeOutcomes = true;
    } else if (strcmp(a, "--rebuild-cgraph") == 0) {
      rebuildCgraph = true;
    } else if (strcmp(a, "--csv-header") == 0) {
      csvHeader = true;
    } else if (strcmp(a, "--quiet") == 0) {
      quiet = true;
    } else if (strcmp(a, "--stage-times") == 0) {
      stageTimes = true;
    } else if (strcmp(a, "--verify-prepare") == 0) {
      verifyPrepare = true;
    } else if (strcmp(a, "--check-validity") == 0) {
      checkValidity = true;
    } else if (strcmp(a, "--no-gate") == 0) {
      useGate = false;
    } else if (strncmp(a, "--rounds=", 9) == 0) {
      maxRounds = atoi(a + 9);
    } else if (strncmp(a, "--ext-method=", 13) == 0) {
      extMethod = (size_t)atol(a + 13);
      if (extMethod > 2) {
        fprintf(stderr, "ERROR: --ext-method must be 0, 1 or 2\n");
        return 1;
      }
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

  if (doSelfTest)
    return selfTest(stem);

  Fixture f;
  if (!loadFixture(f, stem, rebuildCgraph, quiet))
    return 1;

  const CoinConflictGraph *cg = f.si.getCGraph();
  if (!cg) {
    fprintf(stderr, "ERROR: no conflict graph available for %s\n", stem.c_str());
    return 1;
  }

  if (doNodeFeatures) {
    printNodeFeatures(cg, f.si.getColSolution(), f.si.getNumCols(),
      baseName(stem), csvHeader);
    return 0;
  }

  if (doNodeOutcomes) {
    printNodeOutcomes(cg, f.si.getColSolution(), f.si.getReducedCost(),
      f.si.getNumCols(), extMethod, baseName(stem), csvHeader);
    return 0;
  }

  CglTreeInfo info;
  info.level = 0;
  info.pass = 0;
  info.formulation_rows = f.si.getNumRows();
  info.inTree = false;
  info.options = 0;

  const double objStart = f.si.getObjValue();
  const int nRows0 = f.si.getNumRows();

  int totalCuts = 0;
  size_t totalCutLen = 0;
  double totalViol = 0.0, maxViol = 0.0;
  double totalSepTime = 0.0, totalResolveTime = 0.0;
  int totalResolveIters = 0;
  Totals tot;
  std::string cutsPerRound, violPerRound, objImprovePerRound;
  // The bound this round starts from, so each round's own contribution can be
  // reported: a separator whose whole gain arrives in round 1 behaves differently
  // under Cbc's repeated calls than one that keeps paying off.
  double objRoundStart = objStart;

  int round = 0;
  for (; round < maxRounds; ++round) {
    // A fresh generator per round. CglOddWheel keeps no ratchet of its own, but
    // it does keep the scratch arrays and the stats of the last call, and a fresh
    // object makes each round independent of the previous one by construction
    // rather than by inspection.
    CglOddWheel oddWheel(extMethod);
    if (maxSeconds > 0.0)
      oddWheel.setMaxSeconds(maxSeconds);
    if (verifyPrepare)
      oddWheel.setVerifyPrepare(true);
    if (!useGate)
      oddWheel.setUseGate(false);
    if (checkValidity)
      oddWheel.setCheckValidity(true);

    // The solution the cuts are generated against, kept for violation scoring:
    // getColSolution() moves under applyCuts/resolve.
    const std::vector< double > xRound(f.si.getColSolution(),
      f.si.getColSolution() + f.si.getNumCols());

    OsiCuts cs;
    const double t0 = wallClock();
    oddWheel.generateCuts(f.si, cs, info);
    totalSepTime += wallClock() - t0;
    tot.add(oddWheel.stats());

    const int nCuts = cs.sizeRowCuts();
    if (nCuts == 0)
      break;

    double roundViol = 0.0;
    for (int c = 0; c < nCuts; ++c) {
      const OsiRowCut &rc = cs.rowCut(c);
      const double v = rc.violated(xRound.data());
      roundViol += v;
      if (v > maxViol)
        maxViol = v;
      totalCutLen += (size_t)rc.row().getNumElements();
    }
    totalViol += roundViol;
    totalCuts += nCuts;

    char buf[64];
    snprintf(buf, sizeof(buf), "%s%d", round ? "+" : "", nCuts);
    cutsPerRound += buf;
    snprintf(buf, sizeof(buf), "%s%.6g", round ? "+" : "", roundViol);
    violPerRound += buf;

    f.si.applyCuts(cs);

    const double t1 = wallClock();
    f.si.resolve();
    totalResolveTime += wallClock() - t1;
    totalResolveIters += f.si.getIterationCount();

    // Per-round bound gain. A round whose LP does not reach optimality
    // contributes 0 rather than a bound taken from an unsolved LP.
    const double objRoundEnd
      = f.si.isProvenOptimal() ? f.si.getObjValue() : objRoundStart;
    snprintf(buf, sizeof(buf), "%s%.6g", round ? "+" : "",
      f.si.getObjSenseInCbc() * (objRoundEnd - objRoundStart));
    objImprovePerRound += buf;
    objRoundStart = objRoundEnd;

    info.pass = round + 1;
  }

  if (cutsPerRound.empty()) {
    cutsPerRound = "0";
    violPerRound = "0";
    objImprovePerRound = "0";
  }

  const double objEnd = f.si.isProvenOptimal() ? f.si.getObjValue() : objStart;
  const double objImprove = f.si.getObjSenseInCbc() * (objEnd - objStart);
  // Relative to the starting bound's magnitude, which is what makes gains
  // comparable across instances whose objectives differ by orders of magnitude.
  const double objImproveRel
    = objImprove / (fabs(objStart) > 1.0 ? fabs(objStart) : 1.0);

  if (csvHeader)
    printf("%s\n", CSV_HEADER);

  printf("%s,%lu,%d,%d,%d,%.10g,%.10g,%.3f,%.6f,%.6f,%d,%.6f,%.6f,%d,"
         "%.15g,%.15g,%.6g,%.6g,%d,"
         "%lu,%lu,%lu,%.10f,%d,"
         "%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,"
         "%lu,%lu,%.3f,%lu,%lu,%lu,%lu,%lu,"
         "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
         "%.6f,%.6f,%s,%s,%s,"
         "%lu,%lu,%lu,%lu,%lu,%lu,%lu,"
         "%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,"
         "%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,"
         "%lu,%lu,%lu,%lu,%lu,%.6f\n",
    baseName(stem).c_str(), (unsigned long)extMethod, round,
    f.si.getNumRows() - nRows0, totalCuts, totalViol, maxViol,
    totalCuts ? (double)totalCutLen / totalCuts : 0.0,
    totalSepTime, f.warmStartTime, f.warmStartIters, f.cgraphTime,
    totalResolveTime, totalResolveIters, objStart, objEnd, objImprove,
    objImproveRel, (int)boundMoved(objStart, objImprove),
    (unsigned long)cg->size(), (unsigned long)cg->nTotalDirectConflicts(),
    (unsigned long)cg->nCliques(), cg->density(), f.restoredColTypes,
    (unsigned long)tot.icaCount, (unsigned long)tot.arcs,
    (unsigned long)tot.spFindCalls, (unsigned long)tot.oddHoles,
    (unsigned long)tot.ohShort, (unsigned long)tot.ohRepeated,
    (unsigned long)tot.ohNotViol, (unsigned long)tot.ohDuplicate,
    (unsigned long)tot.wheelCenters, (unsigned long)tot.wcElements,
    tot.wheelCenters ? (double)tot.wcElements / tot.wheelCenters : 0.0,
    (unsigned long)tot.cutsBeforePool, (unsigned long)tot.cutsDupIdx,
    (unsigned long)tot.cutsZeroCoefs, (unsigned long)tot.cutsEmpty,
    (unsigned long)tot.cutsAfterPool, (int)tot.timeLimitHit,
    tot.tSetup, tot.tSeparator, tot.tActive, tot.tPrepArcs, tot.tPrepRev,
    tot.tPrepSpf, tot.tSearch, tot.tWheelCenter, tot.tCutPool,
    cutsPerRound.c_str(), violPerRound.c_str(), objImprovePerRound.c_str(),
    (unsigned long)tot.prepMethod, (unsigned long)tot.prepWalkCost,
    (unsigned long)tot.prepUnsorted, (unsigned long)tot.prepVerifyArcs,
    (unsigned long)tot.prepMismatch, (unsigned long)tot.prepWalkOnly,
    (unsigned long)tot.prepPairOnly,
    (unsigned long)tot.certChecked, (unsigned long)tot.certBadCycle,
    (unsigned long)tot.certBadCenterAdj, (unsigned long)tot.certBadCenterClq,
    (unsigned long)tot.certBadAlpha, (unsigned long)tot.certBadTranslate,
    (unsigned long)tot.certComplCycle, (unsigned long)tot.certComplCenter,
    (unsigned long)tot.certComplPair,
    (unsigned long)tot.certCenterOnComplCycle,
    (unsigned long)tot.certComplAtLeastK,
    (unsigned long)tot.wcCalls, (unsigned long)tot.wcPool,
    (unsigned long)tot.wcRejInCycle, (unsigned long)tot.wcRejDegree,
    (unsigned long)tot.wcRejAdjacency, (unsigned long)tot.wcRejCost,
    (unsigned long)tot.wcCandidates, (unsigned long)tot.wcCliqueDropped,
    (unsigned long)tot.gateSkipped, (unsigned long)tot.gateBipartite,
    (unsigned long)tot.gateNoCycle, (unsigned long)tot.gateBlockOnly,
    (unsigned long)tot.gateTriangle, tot.tGate);

  if (stageTimes)
    printStageTimes(tot, totalSepTime);

  return 0;
}
