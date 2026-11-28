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
//     --dgrow=P    each LNS iteration doubles --drows with probability P,
//                  repeatedly (up to --maxdrows); 0 = fixed size
//     --related=P  probability that a further seed row is one freed by the
//                  previous seeds (connected neighbourhood)
//     --relatedafter=0|1  enable related rows only after the first
//                  artificial-free solution (default 0)
//     --btprob=P   probability of supplementing greedy repair with backtracking
//     --btwidth=N --btdepth=N --btnodes=N --btwork=W  search limits
//     --btstats=FILE  aggregate repair diagnostics (TSV, separate from trace)
//     --btonlyart=0|1  skip search when greedy is already artificial-free
//     --btrealfirst=0|1  branch on real columns before artificial columns
//     --gaint=K    seed rows by tournament of K on the gain (0/1 uniform)
//     --pool=K     seed columns drawn among the K*dcols best of the row
//     --rrows=F    probability that a seed row is any partitioning row
//     --amode=M    LNS acceptance: 0 threshold, 1 LAHC, 2 LAHC-better,
//                  3 simulated annealing
//     --lahc=L     LAHC history length
//     --satemp=T   SA initial temperature (fraction of objective scale)
//     --sacool=W   SA temperature half-life in work units
//     --sareheat=N reheat after N consecutive rejected trials (0 disables)
//     --saschedule=0|1 exponential (default) or linear-to-zero cooling
//     --sacalibrate=N threshold-warmup trials for empirical temperature (0 off)
//     --satarget=P --saquantile=Q target acceptance and positive-delta quantile
//     --acceptstats=FILE acceptance counts (separate from repair diagnostics)
//     --noimp=N    stop the LNS after N iterations without a new best (0 never)
//     --maxwork=W  stop the LNS after W units of work (CbcSppGrasp::work(),
//                  0 never); tuning sweeps use --noimp=0 --maxwork=0 and
//                  replay the --trace offline
//     --totalwork=W hard aggregate work cap including artificial-penalty rounds
//     --realstallwork=W  stop after W work without a better real-feasible
//                  incumbent (0 disables; no real incumbent means no stop)
//     --lcalpha=A --lralpha=A --lls=L --lrow=R  column/row RCL tolerance,
//                  local search and row rule of the LNS repair (default:
//                  those of the construction)
//                  --lls=3: 1 before first real feasible incumbent, then 2;
//                  --lls=4: 0 then 2; --lls=5: 2 then 1
//     --focus=0|1  LNS local search only from rows changed by destroy/repair
//     --scan=K     greedy repair examines at most K columns per freed row
//     --artificial=F  artificial slack factor (CbcSppGrasp::setArtificialSlacks),
//                  0 (default) none
//     --artgrow=G --artrounds=R  penalty growth and number of escalations
//                  (RunParams::artificialGrowth/artificialRounds)
//     --artwork=W  escalate after W work per round without a real solution,
//                  independently of the stop rules (artificialRoundWork);
//                  needed to replay traces of open-ended runs
//     --solfile=FILE  write the best solution (if free of artificials) in
//                  Cbc's .sol format, for cbc_validate_sol / mipstart
//     --trace=FILE write instance, lnsIterations:work:time, then
//                  iteration:value:time:work of every new best to FILE
//     --progress=FILE write best/progress/penalty-epoch events separately
//     --features=FILE write cheap structural features and optional cached LP data
//   Cbc's settings (CbcHeuristicSppGrasp): --iters=1 --lns=100000000
//     (other defaults are Cbc's already, except --time, whose Cbc default
//     sppGraspTime is 60)
//     --seed=S
//     --lpcache=DIR  cache the root LP solution there (solved with Clp)
//     --nolp       do not solve/use the LP
//
// A root fixture (<name>.root.mps.gz with its .ctype/.bas/.meta sidecars,
// see ROOT-FIXTURES.md) is loaded as the preprocessed model: integrality
// from .ctype, LP from the saved basis; the instance is named <name>.root.
// Objective values include the model's objective offset.
//
// Prints one tab-separated line:
//   instance trivial first best iterations time time_to_best lp_obj slacks
//   artificial_columns artificial_used (best has artificials if > 0)
//   restricting_rows restricted_columns (CbcSppGrasp::numberRestrictingRows)
//   partitioning_rows slack_rows columns
// --info prints only those structure counts after loading (all with
// --artificial; see numberArtificial):
//   instance INFO partitioning_rows slack_rows artificial_columns
//   restricting_rows restricted_columns columns rows

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "ClpSimplex.hpp"

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

static long metaValue(const std::string &file, const char *key, long missing)
{
  FILE *fp = fopen(file.c_str(), "r");
  if (!fp)
    return missing;
  char k[256];
  double v;
  long value = missing;
  while (fscanf(fp, "%255s %lf", k, &v) == 2)
    if (!strcmp(k, key))
      value = static_cast< long >(v);
  fclose(fp);
  return value;
}

/// root fixture: drop the pad row, restore integrality, install the basis
static bool loadFixture(OsiClpSolverInterface &solver, const std::string &stem)
{
  const std::string meta = stem + ".meta";
  if (metaValue(meta, "paddedColumns", 0) > 0
    && solver.getNumRows() == metaValue(meta, "rows", -1) + 1) {
    int last = solver.getNumRows() - 1;
    solver.deleteRows(1, &last);
  }
  FILE *fp = fopen((stem + ".ctype").c_str(), "r");
  if (!fp) {
    fprintf(stderr, "no %s.ctype\n", stem.c_str());
    return false;
  }
  int cols = -1, idx, type;
  if (fscanf(fp, "cols %d\n", &cols) != 1 || cols != solver.getNumCols()) {
    fclose(fp);
    fprintf(stderr, "bad %s.ctype\n", stem.c_str());
    return false;
  }
  while (fscanf(fp, "%d %d\n", &idx, &type) == 2)
    if (type && idx >= 0 && idx < cols)
      solver.setInteger(idx);
  fclose(fp);
  if (solver.getModelPtr()->readBasis((stem + ".bas").c_str()) >= 0)
    solver.setWarmStart(NULL);
  return true;
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
  double destroyGrowProbability = defaults.destroyGrowProbability;
  double relatedRowProbability = defaults.relatedRowProbability;
  int gainTournament = defaults.gainTournament;
  int seedPool = defaults.seedPool;
  double randomRowFraction = defaults.randomRowFraction;
  CbcSppGrasp::RunParams params = defaults;
  std::string traceFile, progressFile, featuresFile, solFile, btStatsFile, acceptStatsFile;
  double artificial = 0.0;
  bool infoOnly = false;
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
    else if (key == "--btonlyart" || key == "--btrealfirst") {
      if (strcmp(val, "0") != 0 && strcmp(val, "1") != 0) {
        fprintf(stderr, "%s expects 0 or 1\n", key.c_str());
        return 2;
      }
      if (key == "--btonlyart")
        params.backtrackArtificialOnly = strcmp(val, "1") == 0;
      else
        params.backtrackRealFirst = strcmp(val, "1") == 0;
    }
    else if (key == "--btprob")
      params.backtrackProbability = atof(val);
    else if (key == "--btwidth")
      params.backtrackWidth = atoi(val);
    else if (key == "--btdepth")
      params.backtrackDepth = atoi(val);
    else if (key == "--btnodes")
      params.backtrackNodes = atoi(val);
    else if (key == "--btwork")
      params.backtrackWork = atof(val);
    else if (key == "--related")
      relatedRowProbability = atof(val);
    else if (key == "--relatedafter") {
      if (strcmp(val, "0") != 0 && strcmp(val, "1") != 0) {
        fprintf(stderr, "--relatedafter expects 0 or 1\n");
        return 2;
      }
      params.relatedAfterFeasible = strcmp(val, "1") == 0;
    }
    else if (key == "--gaint")
      gainTournament = atoi(val);
    else if (key == "--dgrow")
      destroyGrowProbability = atof(val);
    else if (key == "--pool")
      seedPool = atoi(val);
    else if (key == "--rrows")
      randomRowFraction = atof(val);
    else if (key == "--amode")
      params.acceptMode = atoi(val);
    else if (key == "--lahc")
      params.lahcLength = atoi(val);
    else if (key == "--satemp")
      params.saInitialTemperature = atof(val);
    else if (key == "--sacool")
      params.saCoolingWork = atof(val);
    else if (key == "--sareheat")
      params.saReheatAfter = atoi(val);
    else if (key == "--saschedule")
      params.saSchedule = atoi(val);
    else if (key == "--sacalibrate")
      params.saCalibrationIterations = atoi(val);
    else if (key == "--satarget")
      params.saCalibrationTarget = atof(val);
    else if (key == "--saquantile")
      params.saCalibrationQuantile = atof(val);
    else if (key == "--noimp")
      params.maxNoImprove = atoi(val);
    else if (key == "--maxwork")
      params.maxWork = atof(val);
    else if (key == "--totalwork") {
      char *end = NULL;
      params.totalWorkLimit = strtod(val, &end);
      if (end == val || *end) {
        fprintf(stderr, "--totalwork requires a numeric aggregate work budget\n");
        return 2;
      }
    }
    else if (key == "--realstallwork")
      params.maxRealStallWork = atof(val);
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
    else if (key == "--artificial")
      artificial = atof(val);
    else if (key == "--artgrow")
      params.artificialGrowth = atof(val);
    else if (key == "--artrounds")
      params.artificialRounds = atoi(val);
    else if (key == "--artwork")
      params.artificialRoundWork = atof(val);
    else if (key == "--solfile")
      solFile = val;
    else if (key == "--trace")
      traceFile = val;
    else if (key == "--progress")
      progressFile = val;
    else if (key == "--features")
      featuresFile = val;
    else if (key == "--btstats")
      btStatsFile = val;
    else if (key == "--acceptstats")
      acceptStatsFile = val;
    else if (key == "--info")
      infoOnly = true;
    else if (key == "--nolp")
      useLp = false;
    else {
      fprintf(stderr, "unknown option %s\n", arg.c_str());
      return 2;
    }
  }
  if (!(params.backtrackProbability >= 0.0 && params.backtrackProbability <= 1.0)
    || params.backtrackWidth < 1 || params.backtrackDepth < 0
    || params.backtrackNodes < 1 || !(params.backtrackWork > 0.0)) {
    fprintf(stderr, "invalid backtracking probability or search limits\n");
    return 2;
  }
  if (!std::isfinite(params.maxRealStallWork) || params.maxRealStallWork < 0.0) {
    fprintf(stderr, "invalid real-stall work limit\n");
    return 2;
  }
  if (!std::isfinite(params.totalWorkLimit) || params.totalWorkLimit < 0.0) {
    fprintf(stderr, "invalid aggregate work limit\n");
    return 2;
  }
  if (params.acceptMode < CbcSppGrasp::AcceptThreshold
    || params.acceptMode > CbcSppGrasp::AcceptSimulatedAnnealing
    || !std::isfinite(params.saInitialTemperature) || params.saInitialTemperature < 0.0
    || !std::isfinite(params.saCoolingWork) || params.saCoolingWork <= 0.0
    || params.saReheatAfter < 0 || (params.saSchedule != 0 && params.saSchedule != 1)
    || params.saCalibrationIterations < 0
    || !std::isfinite(params.saCalibrationTarget)
    || params.saCalibrationTarget <= 0.0 || params.saCalibrationTarget >= 1.0
    || !std::isfinite(params.saCalibrationQuantile)
    || params.saCalibrationQuantile < 0.0 || params.saCalibrationQuantile > 1.0) {
    fprintf(stderr, "invalid acceptance mode or simulated-annealing settings\n");
    return 2;
  }
  OsiClpSolverInterface solver;
  solver.messageHandler()->setLogLevel(0);
  solver.getModelPtr()->messageHandler()->setLogLevel(0);
  if (solver.readMps(file.c_str(), "") != 0) {
    fprintf(stderr, "cannot read %s\n", file.c_str());
    return 2;
  }
  // name: file name without .mps[.gz] (fixtures keep their .root)
  std::string name = file.substr(file.find_last_of('/') + 1);
  for (const char *suffix : { ".gz", ".mps" })
    if (name.size() > strlen(suffix) && name.compare(name.size() - strlen(suffix), std::string::npos, suffix) == 0)
      name.erase(name.size() - strlen(suffix));
  const std::string::size_type rootPos = file.rfind(".root.mps");
  const bool fixture = rootPos != std::string::npos;
  if (fixture && !loadFixture(solver, file.substr(0, rootPos + 5)))
    return 2;
  double offset = 0.0;
  solver.getDblParam(OsiObjOffset, offset);
  const int n = solver.getNumCols();
  if (infoOnly) {
    CbcSppGrasp grasp;
    grasp.setArtificialSlacks(artificial);
    double loadStart = CoinGetTimeOfDay();
    if (!grasp.load(solver)) {
      printf("%s\tNA\t%s\n", name.c_str(), grasp.reason().c_str());
      return 1;
    }
    double loadSeconds = CoinGetTimeOfDay() - loadStart;
    printf("%s\tINFO\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", name.c_str(), grasp.numberPartitioningRows(),
      grasp.numberSlackRows() - grasp.numberArtificial(), grasp.numberArtificial(),
      grasp.numberRestrictingRows(), grasp.numberRestrictedColumns(), n, solver.getNumRows());
    if (!featuresFile.empty()) {
      CbcSppGrasp::FeatureSummary features;
      double featureStart = CoinGetTimeOfDay();
      grasp.featureSummary(features);
      double featureSeconds = CoinGetTimeOfDay() - featureStart;
      FILE *fp = fopen(featuresFile.c_str(), "w");
      if (!fp) {
        perror(featuresFile.c_str());
        return 2;
      }
      fprintf(fp, "instance\tpartition_rows\tnative_singleton_rows\tnative_singleton_missing\trow_options_min\trow_options_q25\trow_options_median\trow_options_q75\trow_options_mean\tcolumn_degree_q25\tcolumn_degree_median\tcolumn_degree_q75\tcolumn_degree_mean\tsingleton_premium_q25\tsingleton_premium_median\tsingleton_premium_q75\trestricting_rows\trestricted_columns\tcomplemented_columns\tforced_columns\trestricted_row_fraction\trestricted_column_fraction\tcomplemented_fraction\tforced_fraction\tlp_available\tlp_fractionality_mean\tnative_singleton_lp_mass\tload_seconds\tfeature_seconds\n");
      fprintf(fp, "%s\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.9g\t%.9g\n",
        name.c_str(), features.partitionRows, features.nativeSingletonRows,
        features.nativeSingletonMissing, features.rowOptionsMin, features.rowOptionsQ25,
        features.rowOptionsMedian, features.rowOptionsQ75, features.rowOptionsMean,
        features.columnDegreeQ25, features.columnDegreeMedian, features.columnDegreeQ75,
        features.columnDegreeMean, features.singletonPremiumQ25,
        features.singletonPremiumMedian, features.singletonPremiumQ75,
        features.restrictedRows, features.restrictedColumns,
        features.complementedColumns, features.forcedColumns,
        features.restrictedRowFraction, features.restrictedColumnFraction,
        features.complementedFraction, features.forcedFraction, features.lpAvailable,
        features.lpFractionalityMean, features.nativeSingletonLpMass, loadSeconds,
        featureSeconds);
      fclose(fp);
    }
    return 0;
  }
  std::vector< double > x, dj;
  double lpObj = 0.0;
  if (useLp) {
    std::string cacheFile = lpCache.empty() ? "" : lpCache + "/" + name + ".lpsol";
    if (cacheFile.empty() || !readCache(cacheFile, n, x, dj, lpObj)) {
      if (fixture) {
        solver.setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
        solver.setHintParam(OsiDoDualInResolve, true, OsiHintDo);
        solver.resolve();
      } else {
        solver.initialSolve();
      }
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
  grasp.setArtificialSlacks(artificial);
  double loadStart = CoinGetTimeOfDay();
  if (!grasp.load(solver)) {
    printf("%s\tNA\t%s\n", name.c_str(), grasp.reason().c_str());
    return 1;
  }
  double loadSeconds = CoinGetTimeOfDay() - loadStart;
  if (useLp)
    grasp.setLpSolution(x.data(), dj.data());
  if (!featuresFile.empty()) {
    CbcSppGrasp::FeatureSummary features;
    double featureStart = CoinGetTimeOfDay();
    grasp.featureSummary(features);
    double featureSeconds = CoinGetTimeOfDay() - featureStart;
    FILE *fp = fopen(featuresFile.c_str(), "w");
    if (!fp) {
      perror(featuresFile.c_str());
      return 2;
    }
    fprintf(fp, "instance\tpartition_rows\tnative_singleton_rows\tnative_singleton_missing\trow_options_min\trow_options_q25\trow_options_median\trow_options_q75\trow_options_mean\tcolumn_degree_q25\tcolumn_degree_median\tcolumn_degree_q75\tcolumn_degree_mean\tsingleton_premium_q25\tsingleton_premium_median\tsingleton_premium_q75\trestricting_rows\trestricted_columns\tcomplemented_columns\tforced_columns\trestricted_row_fraction\trestricted_column_fraction\tcomplemented_fraction\tforced_fraction\tlp_available\tlp_fractionality_mean\tnative_singleton_lp_mass\tload_seconds\tfeature_seconds\n");
    fprintf(fp, "%s\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.9g\t%.9g\n",
      name.c_str(), features.partitionRows, features.nativeSingletonRows,
      features.nativeSingletonMissing, features.rowOptionsMin, features.rowOptionsQ25,
      features.rowOptionsMedian, features.rowOptionsQ75, features.rowOptionsMean,
      features.columnDegreeQ25, features.columnDegreeMedian, features.columnDegreeQ75,
      features.columnDegreeMean, features.singletonPremiumQ25,
      features.singletonPremiumMedian, features.singletonPremiumQ75,
      features.restrictedRows, features.restrictedColumns,
      features.complementedColumns, features.forcedColumns,
      features.restrictedRowFraction, features.restrictedColumnFraction,
      features.complementedFraction, features.forcedFraction, features.lpAvailable,
      features.lpFractionalityMean, features.nativeSingletonLpMass, loadSeconds,
      featureSeconds);
    fclose(fp);
  }
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
  params.destroyGrowProbability = destroyGrowProbability;
  params.relatedRowProbability = relatedRowProbability;
  params.gainTournament = gainTournament;
  params.seedPool = seedPool;
  params.randomRowFraction = randomRowFraction;
  CbcSppGrasp::RunStats stats;
  double bestValue = grasp.run(config, params, best, &stats);
  if (!acceptStatsFile.empty()) {
    FILE *fp = fopen(acceptStatsFile.c_str(), "w");
    if (!fp) {
      perror(acceptStatsFile.c_str());
      return 2;
    }
    const CbcSppGrasp::AcceptanceStats &a = stats.acceptance;
    fprintf(fp, "instance\ttrials\taccepted\taccepted_worse\trejected\tinvalid\treheats\tcalibration_trials\tcalibration_positive\tinitial_temperature\n");
    fprintf(fp, "%s\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%d\t%d\t%.12g\n", name.c_str(),
      a.trials, a.accepted, a.acceptedWorse, a.rejected, a.invalid, a.reheats,
      a.calibrationTrials, a.calibrationPositive, a.calibratedTemperature);
    if (fclose(fp) != 0) {
      perror(acceptStatsFile.c_str());
      return 2;
    }
  }
  if (!btStatsFile.empty()) {
    FILE *fp = fopen(btStatsFile.c_str(), "w");
    if (!fp) {
      perror(btStatsFile.c_str());
      return 2;
    }
    const CbcSppGrasp::RepairStats &r = stats.repair;
    fprintf(fp, "instance\tcalls\tnodes\tleaves\treal_leaves\tdead_ends\twork_exits\tnode_exits\texhausted\tinvalid_partial\timproved\tfeasibility_gains\tgreedy_artificial\tno_leaf\topen_rows\tmax_open_rows\tmax_depth\textra_work\tskipped_feasible\n");
    fprintf(fp, "%s\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%.0f\t%lld\n",
      name.c_str(), r.calls, r.nodes, r.leaves, r.realLeaves, r.deadEnds,
      r.workExits, r.nodeExits, r.exhausted, r.invalidPartial, r.improved,
      r.feasibilityGains, r.greedyArtificial, r.noLeaf, r.openRows,
      r.maxOpenRows, r.maxDepth, r.work, r.skippedFeasible);
    if (fclose(fp) != 0) {
      perror(btStatsFile.c_str());
      return 2;
    }
  }
  double trivial = stats.trivial, first = stats.first, timeToBest = stats.timeToBest;
  int it = stats.graspIterations * 100000 + stats.lnsIterations;
  if (!traceFile.empty()) {
    // one line: instance, then iteration:value:time of every new best
    FILE *fp = fopen(traceFile.c_str(), "w");
    if (fp) {
      fprintf(fp, "%s\t%d:%.0f:%.4f", name.c_str(), stats.lnsIterations, stats.work, stats.time);
      for (size_t k = 0; k < stats.improvements.size(); k++)
        fprintf(fp, "\t%d:%.6f:%.4f:%.0f", stats.improvements[k].iteration,
          stats.improvements[k].value - offset, stats.improvements[k].time, stats.improvements[k].work);
      fprintf(fp, "\n");
      fclose(fp);
    }
  }
  if (!progressFile.empty()) {
    FILE *fp = fopen(progressFile.c_str(), "w");
    if (!fp) {
      perror(progressFile.c_str());
      return 2;
    }
    fprintf(fp, "iteration\tkind\tpenalty_round\tartificial_used\tvalue\ttime\twork\n");
    for (size_t k = 0; k < stats.progress.size(); k++) {
      const CbcSppGrasp::ProgressEvent &event = stats.progress[k];
      fprintf(fp, "%d\t%d\t%d\t%d\t%.12g\t%.6f\t%.0f\n", event.iteration,
        event.kind, event.penaltyRound, event.artificialUsed, event.value - offset,
        event.time, event.work);
    }
    fclose(fp);
  }
  if (bestValue != COIN_DBL_MAX && useLp) {
    double cx = 0.0;
    for (int j = 0; j < n; j++)
      cx += x[j] * solver.getObjCoefficients()[j];
    if (fabs(cx - offset - lpObj) > 1.0e-5 * (1.0 + fabs(lpObj)))
      fprintf(stderr, "%s: objective offset mismatch (c x %g, offset %g, LP %g)\n",
        name.c_str(), cx, offset, lpObj);
  }
  if (bestValue != COIN_DBL_MAX && !grasp.numberArtificialUsed(best)) {
    // the solution mapped back must satisfy every row of the model
    std::vector< double > xs(n), activity(solver.getNumRows());
    grasp.toOriginal(best, xs.data());
    solver.getMatrixByCol()->times(xs.data(), activity.data());
    double cx = 0.0;
    for (int j = 0; j < n; j++) {
      cx += xs[j] * solver.getObjCoefficients()[j];
      if (xs[j] < solver.getColLower()[j] - 1.0e-9 || xs[j] > solver.getColUpper()[j] + 1.0e-9)
        fprintf(stderr, "%s: column %d at %g outside bounds\n", name.c_str(), j, xs[j]);
    }
    int bad = 0;
    for (int i = 0; i < solver.getNumRows(); i++)
      if (activity[i] < solver.getRowLower()[i] - 1.0e-7 || activity[i] > solver.getRowUpper()[i] + 1.0e-7)
        bad++;
    if (bad || fabs(cx * solver.getObjSense() - bestValue) > 1.0e-6 * (1.0 + fabs(bestValue)))
      fprintf(stderr, "%s: mapped solution: %d rows violated, c x %.6f vs %.6f\n", name.c_str(),
        bad, cx, bestValue);
    FILE *fp = solFile.empty() ? NULL : fopen(solFile.c_str(), "w");
    if (fp) {
      fprintf(fp, "Stopped on iterations - objective value %.12g\n", cx - offset);
      for (int j = 0; j < n; j++)
        fprintf(fp, "%d %s %.0f %.12g\n", j, solver.getColName(j).c_str(), xs[j],
          solver.getObjCoefficients()[j]);
      fclose(fp);
    }
  }
  double check = grasp.evaluate(best);
  if (check == COIN_DBL_MAX || fabs(check - bestValue) > 1.0e-6 * (1.0 + fabs(bestValue))) {
    fprintf(stderr, "%s: invalid best solution (%g vs %g)\n", name.c_str(), check, bestValue);
    return 1;
  }
  // the value of a solution is c x - offset
  auto shift = [&](double v) { return v == COIN_DBL_MAX ? v : v - offset; };
  printf("%s\t%.6f\t%.6f\t%.6f\t%d\t%.4f\t%.4f\t%.6f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", name.c_str(), shift(trivial),
    shift(first), shift(bestValue), it, CoinGetTimeOfDay() - start, timeToBest, lpObj,
    grasp.numberSlacksUsed(best), grasp.numberArtificial(), grasp.numberArtificialUsed(best),
    grasp.numberRestrictingRows(), grasp.numberRestrictedColumns(), grasp.numberPartitioningRows(),
    grasp.numberSlackRows() - grasp.numberArtificial(), n);
  return 0;
}
