// Copyright (C) 2026, COIN-OR Foundation and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// spp-grasp-bench: runs the set partitioning/packing GRASP heuristic
// (CbcSppGrasp, see CbcHeuristicSppGrasp.hpp) standalone on one instance,
// for tuning its criteria offline without a full Cbc run.
//
//   spp-grasp-bench <file.mps[.gz]> [options]
//     --row=R      row rule (0 min options, 1 max regret, 2 options/regret,
//                  3 log(1+options)/regret, 4 min best, 5 max best, 6 random)
//     --col=C      column rule (0 cost per row, 1 reduced cost, 2 profit)
//     --beta=B     cost-per-row exponent
//     --calpha=A   column RCL relative tolerance     --crcl=K  column RCL size
//     --ralpha=A   row RCL relative tolerance
//     --lp=W       LP value bias weight              --conf=W  conflict weight
//     --ls=0|1|2   local search (2: greedy repair)
//     --iters=N    GRASP iterations (first one is deterministic)
//     --time=T     time limit in seconds
//     --lns=N      LNS (destroy/rebuild) iterations after GRASP
//     --drows=N    slack rows around which LNS destroys
//     --dcols=N    best options per such row whose conflicts are removed
//     --dfrac=F    probability of removing any other column
//     --accept=F   LNS accepts solutions within relative F of the best
//     --stall=N    LNS iterations without a new best before doubling --drows
//                  (up to --maxdrows, then restart from the best)
//     --pool=K     seed columns drawn among the K*dcols best of the row
//     --rrows=F    probability that a seed row is any partitioning row
//     --amode=M    LNS acceptance: 0 threshold (--accept), 1 LAHC,
//                  2 LAHC with history entries only replaced by better values
//     --lahc=L     LAHC history length
//     --noimp=N    stop the LNS after N iterations without a new best (0 never)
//     --maxwork=W  stop the LNS after W units of work (CbcSppGrasp::work(),
//                  0 never); tuning sweeps use --noimp=0 --maxwork=0 and
//                  replay the --trace offline
//     --lcalpha=A --lralpha=A --lls=L --lrow=R  column/row RCL tolerance,
//                  local search and row rule of the LNS repair (default:
//                  those of the construction)
//     --focus=0|1  LNS local search only from rows changed by destroy/repair
//     --scan=K     greedy repair examines at most K columns per freed row
//     --trace=FILE write instance, lnsIterations:work:time, then
//                  iteration:value:time:work of every new best to FILE
//   Cbc's settings (CbcHeuristicSppGrasp): --iters=1 --lns=100000000
//     (other defaults are Cbc's already, except --time, whose Cbc default
//     sppGraspTime is 60)
//     --seed=S
//     --lpcache=DIR  cache the root LP solution there (solved with Clp)
//     --nolp       do not solve/use the LP
//
// Prints one tab-separated line:
//   instance trivial first best iterations time time_to_best lp_obj slacks

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "CbcHeuristicSppGrasp.hpp"
#include "CoinTime.hpp"
#include "OsiClpSolverInterface.hpp"

static bool readCache(const std::string &file, int n, std::vector< double > &x,
  std::vector< double > &dj, double &obj)
{
  FILE *fp = fopen(file.c_str(), "rb");
  if (!fp)
    return false;
  int n2 = 0;
  bool ok = fread(&n2, sizeof(int), 1, fp) == 1 && n2 == n;
  x.resize(n);
  dj.resize(n);
  ok = ok && fread(&obj, sizeof(double), 1, fp) == 1;
  ok = ok && fread(x.data(), sizeof(double), n, fp) == static_cast< size_t >(n);
  ok = ok && fread(dj.data(), sizeof(double), n, fp) == static_cast< size_t >(n);
  fclose(fp);
  return ok;
}

static void writeCache(const std::string &file, int n, const double *x,
  const double *dj, double obj)
{
  FILE *fp = fopen(file.c_str(), "wb");
  if (!fp)
    return;
  fwrite(&n, sizeof(int), 1, fp);
  fwrite(&obj, sizeof(double), 1, fp);
  fwrite(x, sizeof(double), n, fp);
  fwrite(dj, sizeof(double), n, fp);
  fclose(fp);
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "usage: %s <file.mps[.gz]> [options] (see source header)\n", argv[0]);
    return 2;
  }
  std::string file = argv[1];
  CbcSppGrasp::Config config;
  int iterations = 100;
  double timeLimit = 1.0e10; // iterations column: grasp*100000 + lns
  unsigned int seed = 1;
  std::string lpCache;
  bool useLp = true;
  CbcSppGrasp::RunParams defaults;
  int lnsIterations = 0, destroyRows = defaults.destroyRows, seedColumns = defaults.seedColumns;
  double destroyFraction = defaults.destroyFraction, acceptRelative = defaults.acceptRelative;
  int stallIterations = defaults.stallIterations, maxDestroyRows = defaults.maxDestroyRows;
  int seedPool = defaults.seedPool;
  double randomRowFraction = defaults.randomRowFraction;
  CbcSppGrasp::RunParams params = defaults;
  std::string traceFile;
  for (int a = 2; a < argc; a++) {
    std::string arg = argv[a];
    std::string::size_type eq = arg.find('=');
    std::string key = arg.substr(0, eq);
    const char *val = eq == std::string::npos ? "" : argv[a] + eq + 1;
    if (key == "--row")
      config.rowRule = atoi(val);
    else if (key == "--col")
      config.colRule = atoi(val);
    else if (key == "--beta")
      config.beta = atof(val);
    else if (key == "--calpha")
      config.colAlpha = atof(val);
    else if (key == "--crcl")
      config.colRclMax = atoi(val);
    else if (key == "--ralpha")
      config.rowAlpha = atof(val);
    else if (key == "--lp")
      config.lpWeight = atof(val);
    else if (key == "--conf")
      config.conflictWeight = atof(val);
    else if (key == "--ls")
      config.localSearch = atoi(val);
    else if (key == "--iters")
      iterations = atoi(val);
    else if (key == "--time")
      timeLimit = atof(val);
    else if (key == "--seed")
      seed = static_cast< unsigned int >(atoi(val));
    else if (key == "--lpcache")
      lpCache = val;
    else if (key == "--lns")
      lnsIterations = atoi(val);
    else if (key == "--drows")
      destroyRows = atoi(val);
    else if (key == "--dcols")
      seedColumns = atoi(val);
    else if (key == "--accept")
      acceptRelative = atof(val);
    else if (key == "--dfrac")
      destroyFraction = atof(val);
    else if (key == "--stall")
      stallIterations = atoi(val);
    else if (key == "--maxdrows")
      maxDestroyRows = atoi(val);
    else if (key == "--pool")
      seedPool = atoi(val);
    else if (key == "--rrows")
      randomRowFraction = atof(val);
    else if (key == "--amode")
      params.acceptMode = atoi(val);
    else if (key == "--lahc")
      params.lahcLength = atoi(val);
    else if (key == "--noimp")
      params.maxNoImprove = atoi(val);
    else if (key == "--maxwork")
      params.maxWork = atof(val);
    else if (key == "--lcalpha")
      params.lnsColAlpha = atof(val);
    else if (key == "--lralpha")
      params.lnsRowAlpha = atof(val);
    else if (key == "--lls")
      params.lnsLocalSearch = atoi(val);
    else if (key == "--lrow")
      params.lnsRowRule = atoi(val);
    else if (key == "--focus")
      params.lnsFocus = atoi(val) != 0;
    else if (key == "--scan")
      config.repairScan = atoi(val);
    else if (key == "--trace")
      traceFile = val;
    else if (key == "--nolp")
      useLp = false;
    else {
      fprintf(stderr, "unknown option %s\n", arg.c_str());
      return 2;
    }
  }
  OsiClpSolverInterface solver;
  solver.messageHandler()->setLogLevel(0);
  solver.getModelPtr()->messageHandler()->setLogLevel(0);
  if (solver.readMps(file.c_str(), "") != 0) {
    fprintf(stderr, "cannot read %s\n", file.c_str());
    return 2;
  }
  std::string name = file.substr(file.find_last_of('/') + 1);
  name = name.substr(0, name.find('.'));
  const int n = solver.getNumCols();
  std::vector< double > x, dj;
  double lpObj = 0.0;
  if (useLp) {
    std::string cacheFile = lpCache.empty() ? "" : lpCache + "/" + name + ".lpsol";
    if (cacheFile.empty() || !readCache(cacheFile, n, x, dj, lpObj)) {
      solver.initialSolve();
      if (!solver.isProvenOptimal()) {
        fprintf(stderr, "LP not optimal\n");
        return 2;
      }
      lpObj = solver.getObjValue();
      x.assign(solver.getColSolution(), solver.getColSolution() + n);
      dj.assign(solver.getReducedCost(), solver.getReducedCost() + n);
      if (!cacheFile.empty())
        writeCache(cacheFile, n, x.data(), dj.data(), lpObj);
    }
  }
  double start = CoinGetTimeOfDay();
  CbcSppGrasp grasp;
  if (!grasp.load(solver)) {
    printf("%s\tNA\t%s\n", name.c_str(), grasp.reason().c_str());
    return 1;
  }
  if (useLp)
    grasp.setLpSolution(x.data(), dj.data());
  grasp.setSeed(seed);
  std::vector< char > best;
  params.graspIterations = iterations;
  params.lnsIterations = lnsIterations;
  params.timeLimit = timeLimit;
  params.destroyRows = destroyRows;
  params.seedColumns = seedColumns;
  params.destroyFraction = destroyFraction;
  params.acceptRelative = acceptRelative;
  params.stallIterations = stallIterations;
  params.maxDestroyRows = maxDestroyRows;
  params.seedPool = seedPool;
  params.randomRowFraction = randomRowFraction;
  CbcSppGrasp::RunStats stats;
  double bestValue = grasp.run(config, params, best, &stats);
  double trivial = stats.trivial, first = stats.first, timeToBest = stats.timeToBest;
  int it = stats.graspIterations * 100000 + stats.lnsIterations;
  if (!traceFile.empty()) {
    // one line: instance, then iteration:value:time of every new best
    FILE *fp = fopen(traceFile.c_str(), "w");
    if (fp) {
      fprintf(fp, "%s\t%d:%.0f:%.4f", name.c_str(), stats.lnsIterations, stats.work, stats.time);
      for (size_t k = 0; k < stats.improvements.size(); k++)
        fprintf(fp, "\t%d:%.6f:%.4f:%.0f", stats.improvements[k].iteration,
          stats.improvements[k].value, stats.improvements[k].time, stats.improvements[k].work);
      fprintf(fp, "\n");
      fclose(fp);
    }
  }
  double check = grasp.evaluate(best);
  if (check == COIN_DBL_MAX || fabs(check - bestValue) > 1.0e-6 * (1.0 + fabs(bestValue))) {
    fprintf(stderr, "%s: invalid best solution (%g vs %g)\n", name.c_str(), check, bestValue);
    return 1;
  }
  printf("%s\t%.6f\t%.6f\t%.6f\t%d\t%.4f\t%.4f\t%.6f\t%d\n", name.c_str(), trivial, first, bestValue,
    it, CoinGetTimeOfDay() - start, timeToBest, lpObj, grasp.numberSlacksUsed(best));
  return 0;
}
