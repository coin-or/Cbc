/**
 * impclique-explore -- exploratory separator for "implication clique" cuts.
 *
 * Motivating pattern: a row
 *
 *     x_1 + x_2 + ... + x_k <= M * y        (all binary)
 *
 * modelling "x_1 OR x_2 OR ... -> y" is usually left with M == k after plain
 * coefficient tightening (CbcCoefficientStrengthening). That is still weaker
 * than disaggregating it into k rows x_i <= y, and weaker still than grouping
 * mutually-conflicting x_i's together: if x_i and x_j can never both be 1,
 *
 *     x_i + x_j <= y
 *
 * is valid and dominates both x_i <= y and x_j <= y whenever x_i*+x_j* > y*.
 *
 * The key algebraic fact this tool exploits: x_i <= y  <=>  x_i and (1-y)
 * conflict (cannot both be 1), i.e. x_i's literal node conflicts with y's
 * *complement* node in the conflict graph. And a whole clique C of mutually
 * conflicting x_i's, all of which also conflict with (1-y), gives exactly
 * the valid cut  sum_{i in C} x_i <= y  -- because {C, complement(y)} is then
 * itself a clique of literals, i.e. the standard set-packing inequality
 * sum_{i in C} x_i + (1 - y) <= 1 rewritten.
 *
 * So the separator is: for every column y, root the search at node
 * complement(y) in the conflict graph, restrict candidates to y's neighbours
 * that are themselves *primal* column literals with positive LP value, and
 * greedily grow a clique among them (weighted by LP value). This is exactly
 * CglBKClique's clique machinery but *rooted* at a forced hub node, with a
 * variable (not constant-1) right-hand side -- it is not currently produced
 * by CglBKClique itself, whose candidate pool is the global set of fractional
 * primal/complement literals (a y sitting at an integer bound with a nonzero
 * complement is not visited unless it happens to be fractional too).
 *
 * This tool is deliberately NOT a CglCutGenerator yet. It follows the
 * BENCHMARKING-CUT-GENERATORS.md methodology: measure whether the idea has
 * any signal (violation, bound improvement) before writing the production
 * generator. Loads a raw .mps.gz, solves the LP relaxation, builds the
 * conflict graph via checkCGraph(), separates, adds the found cuts, and
 * reports the bound improvement on reoptimizing -- the only trustworthy
 * metric per that document.
 *
 * Usage:
 *   impclique-explore <file.mps.gz> [--csv] [--singleton-only] [--quiet]
 *
 * --singleton-only disables clique growth (k=1 disaggregation only), so
 * running the same instance with/without it isolates what the clique
 * extension itself buys, exactly like CglBKClique's --ext-method=0 control.
 */

#include "CoinConflictGraph.hpp"
#include "CoinStaticConflictGraph.hpp"
#include "CoinTime.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

const double EPS = 1e-6;

struct FoundCut {
  int hubCol;            // y
  std::vector< int > members; // x_i's (>=1)
  double lhs;            // sum x_i* at separation time
  double rhs;            // y* at separation time
  double violation() const { return lhs - rhs; }
};

/**
 * Greedily grow a clique among `candidates` (already known to all conflict
 * with hubComp), ordered by descending LP weight, keeping only members that
 * conflict with everything already accepted. O(|candidates|^2) conflict
 * checks -- fine for exploration; a production generator would reuse
 * CoinBronKerbosch for an exact/pivoted search instead.
 */
std::vector< int > greedyClique(const CoinConflictGraph *cg,
  const std::vector< std::pair< double, int > > &sortedCandidates)
{
  std::vector< int > clique;
  clique.reserve(sortedCandidates.size());
  for (const auto &wc : sortedCandidates) {
    const int cand = wc.second;
    bool ok = true;
    for (int m : clique) {
      if (!cg->conflicting(static_cast< size_t >(cand), static_cast< size_t >(m))) {
        ok = false;
        break;
      }
    }
    if (ok)
      clique.push_back(cand);
  }
  return clique;
}

std::vector< FoundCut > separate(const OsiSolverInterface &si, bool singletonOnly)
{
  std::vector< FoundCut > cuts;
  const CoinConflictGraph *cg = si.getCGraph();
  if (!cg)
    return cuts;

  const int numCols = si.getNumCols();
  const double *sol = si.getColSolution();
  const char *colType = si.getColType(true);

  std::vector< size_t > temp(cg->size());
  std::vector< char > iv(cg->size(), 0);

  for (int y = 0; y < numCols; y++) {
    if (colType[y] == 0) // not binary/integer-01
      continue;
    const double yVal = sol[y];
    if (yVal > 1.0 - EPS) // y already at 1: sum x_i <= 1 is never violated by this route
      continue;

    const size_t yComp = static_cast< size_t >(y) + numCols;
    if (yComp >= cg->size())
      continue;

    auto conf = cg->conflictingNodes(yComp, temp.data(), iv.data());

    std::vector< std::pair< double, int > > candidates; // (weight, col) desc by weight
    for (size_t k = 0; k < conf.first; k++) {
      const size_t node = conf.second[k];
      if (node >= static_cast< size_t >(numCols))
        continue; // skip complement literals; keep primal x_i's only
      const int col = static_cast< int >(node);
      if (col == y)
        continue;
      if (colType[col] == 0)
        continue;
      const double xv = sol[col];
      if (xv > EPS)
        candidates.push_back(std::make_pair(xv, col));
    }
    if (candidates.empty())
      continue;

    std::sort(candidates.begin(), candidates.end(),
      [](const std::pair< double, int > &a, const std::pair< double, int > &b) {
        return a.first > b.first;
      });

    std::vector< int > clique;
    if (singletonOnly) {
      clique.push_back(candidates.front().second);
    } else {
      clique = greedyClique(cg, candidates);
    }

    double lhs = 0.0;
    for (int c : clique)
      lhs += sol[c];

    if (lhs - yVal > EPS) {
      FoundCut fc;
      fc.hubCol = y;
      fc.members = clique;
      fc.lhs = lhs;
      fc.rhs = yVal;
      cuts.push_back(fc);
    }
  }
  return cuts;
}

OsiRowCut toRowCut(const FoundCut &fc, int numCols)
{
  // sum_{i in members} x_i - y <= 0
  std::vector< int > idx(fc.members);
  std::vector< double > val(idx.size(), 1.0);
  idx.push_back(fc.hubCol);
  val.push_back(-1.0);
  OsiRowCut rc;
  rc.setRow(static_cast< int >(idx.size()), idx.data(), val.data());
  rc.setLb(-COIN_DBL_MAX);
  rc.setUb(0.0);
  return rc;
}

} // namespace

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "usage: %s <file.mps.gz> [--csv] [--singleton-only] [--quiet]\n", argv[0]);
    return 2;
  }
  std::string file = argv[1];
  bool csv = false, singletonOnly = false, quiet = false;
  std::string dumpAugmented;
  for (int i = 2; i < argc; i++) {
    if (!strcmp(argv[i], "--csv"))
      csv = true;
    else if (!strcmp(argv[i], "--singleton-only"))
      singletonOnly = true;
    else if (!strcmp(argv[i], "--quiet"))
      quiet = true;
    else if (!strncmp(argv[i], "--dump-augmented=", 17))
      dumpAugmented = argv[i] + 17;
  }

  OsiClpSolverInterface si;
  si.getModelPtr()->setLogLevel(0);
  si.messageHandler()->setLogLevel(0);
  if (si.readMps(file.c_str()) < 0) {
    fprintf(stderr, "failed to read %s\n", file.c_str());
    return 2;
  }

  const double t0 = CoinGetTimeOfDay();
  si.initialSolve();
  if (!si.isProvenOptimal()) {
    fprintf(stderr, "%s: root LP not optimal, skipping\n", file.c_str());
    return 0;
  }
  const double objStart = si.getObjValue();

  si.checkCGraph(NULL);
  const double cgraphTime = CoinGetTimeOfDay() - t0;
  const CoinConflictGraph *cg = si.getCGraph();

  const double tSep0 = CoinGetTimeOfDay();
  std::vector< FoundCut > cuts = separate(si, singletonOnly);
  const double sepTime = CoinGetTimeOfDay() - tSep0;

  double totalViol = 0.0, maxViol = 0.0;
  size_t sumLen = 0, nMultiMember = 0;
  for (const auto &c : cuts) {
    totalViol += c.violation();
    maxViol = std::max(maxViol, c.violation());
    sumLen += c.members.size();
    if (c.members.size() > 1)
      nMultiMember++;
  }
  const double avgLen = cuts.empty() ? 0.0 : static_cast< double >(sumLen) / cuts.size();

  double objEnd = objStart;
  double resolveTime = 0.0;
  if (!cuts.empty()) {
    OsiCuts osc;
    for (const auto &c : cuts)
      osc.insert(toRowCut(c, si.getNumCols()));
    si.applyCuts(osc);
    const double tR0 = CoinGetTimeOfDay();
    si.resolve();
    resolveTime = CoinGetTimeOfDay() - tR0;
    if (si.isProvenOptimal())
      objEnd = si.getObjValue();
  }

  if (!dumpAugmented.empty())
    si.writeMps(dumpAugmented.c_str(), "");

  const double objImprove = objEnd - objStart;
  const double scale = std::max(1.0, std::fabs(objStart));
  const double objImproveRel = objImprove / scale;
  const bool boundMoved = objImproveRel > 1e-9;

  const size_t slash = file.rfind('/');
  const std::string name = slash == std::string::npos ? file : file.substr(slash + 1);

  if (csv) {
    if (!quiet)
      printf("name,cgraphSize,cgraphNodes,cgraphTime,totalCuts,multiMemberCuts,avgCutLen,"
             "totalViol,maxViol,sepTime,resolveTime,objStart,objEnd,objImprove,objImproveRel,boundMoved\n");
    printf("%s,%zu,%zu,%.4f,%zu,%zu,%.3f,%.6f,%.6f,%.4f,%.4f,%.6f,%.6f,%.6g,%.6g,%d\n",
      name.c_str(), cg ? cg->size() : 0, cg ? cg->nTotalDirectConflicts() : 0, cgraphTime,
      cuts.size(), nMultiMember, avgLen, totalViol, maxViol, sepTime, resolveTime,
      objStart, objEnd, objImprove, objImproveRel, boundMoved ? 1 : 0);
  } else {
    printf("%-40s cgraph=%zu nodes cuts=%-4zu (multi=%-3zu avgLen=%.2f) "
           "totalViol=%.4f objStart=%.6f objEnd=%.6f objImprove=%.6g%s\n",
      name.c_str(), cg ? cg->size() : 0, cuts.size(), nMultiMember, avgLen,
      totalViol, objStart, objEnd, objImprove, boundMoved ? " (MOVED)" : "");
  }

  return 0;
}
