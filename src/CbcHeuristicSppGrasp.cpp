// Copyright (C) 2026, COIN-OR Foundation and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <sstream>

#include "CbcHeuristicSppGrasp.hpp"
#include "CbcSppGraspAcceptance.hpp"
#include "CbcModel.hpp"
#include "CbcOutput.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinTime.hpp"
#include "OsiSolverInterface.hpp"

namespace {
const double kNoOption = 1.0e20;

bool applyCoverageTuningValue(const std::string &key, double value,
  CbcSppGrasp::Config &config, CbcSppGrasp::RunParams &params, double &artificial)
{
  if (!std::isfinite(value))
    return false;
#define TUNE_REAL(KEY, FIELD, LOW, HIGH) \
  if (key == KEY) {                      \
    if (value < LOW || value > HIGH)     \
      return false;                      \
    FIELD = value;                       \
    return true;                         \
  }
#define TUNE_INT(KEY, FIELD, LOW, HIGH)                       \
  if (key == KEY) {                                           \
    if (value < LOW || value > HIGH || value != floor(value)) \
      return false;                                           \
    FIELD = static_cast< int >(value);                        \
    return true;                                              \
  }
  TUNE_INT("row", config.rowRule, 0, 6)
  TUNE_INT("col", config.colRule, 0, 2)
  TUNE_REAL("beta", config.beta, 0, 10)
  TUNE_REAL("calpha", config.colAlpha, 0, 10)
  TUNE_INT("crcl", config.colRclMax, 1, COIN_INT_MAX)
  TUNE_REAL("ralpha", config.rowAlpha, 0, 10)
  TUNE_REAL("lp", config.lpWeight, 0, 1.0e6)
  TUNE_REAL("conf", config.conflictWeight, 0, 1.0e6)
  TUNE_INT("ls", config.localSearch, 0, 2)
  TUNE_INT("scan", config.repairScan, 1, COIN_INT_MAX)
  TUNE_INT("iters", params.graspIterations, 1, 10000)
  TUNE_INT("lns", params.lnsIterations, 0, COIN_INT_MAX)
  TUNE_INT("drows", params.destroyRows, 1, COIN_INT_MAX)
  TUNE_INT("dcols", params.seedColumns, 1, COIN_INT_MAX)
  TUNE_REAL("dfrac", params.destroyFraction, 0, 1)
  TUNE_REAL("accept", params.acceptRelative, 0, 1)
  TUNE_INT("stall", params.stallIterations, 0, COIN_INT_MAX)
  TUNE_INT("maxdrows", params.maxDestroyRows, 1, COIN_INT_MAX)
  TUNE_REAL("dgrow", params.destroyGrowProbability, 0, 1)
  TUNE_INT("pool", params.seedPool, 1, 10000)
  TUNE_REAL("rrows", params.randomRowFraction, 0, 1)
  TUNE_REAL("related", params.relatedRowProbability, 0, 1)
  TUNE_INT("relatedafter", params.relatedAfterFeasible, 0, 1)
  TUNE_INT("gaint", params.gainTournament, 0, 10000)
  TUNE_REAL("btprob", params.backtrackProbability, 0, 1)
  TUNE_INT("btonlyart", params.backtrackArtificialOnly, 0, 1)
  TUNE_INT("btrealfirst", params.backtrackRealFirst, 0, 1)
  TUNE_INT("btwidth", params.backtrackWidth, 1, 10000)
  TUNE_INT("btdepth", params.backtrackDepth, 0, 10000)
  TUNE_INT("btnodes", params.backtrackNodes, 1, COIN_INT_MAX)
  TUNE_REAL("btwork", params.backtrackWork, 1, 1.0e12)
  TUNE_INT("amode", params.acceptMode, 0, 3)
  TUNE_INT("lahc", params.lahcLength, 1, 1000000)
  TUNE_REAL("satemp", params.saInitialTemperature, 0, 1.0e6)
  TUNE_REAL("sacool", params.saCoolingWork, 1, 1.0e15)
  TUNE_INT("sareheat", params.saReheatAfter, 0, COIN_INT_MAX)
  TUNE_INT("saschedule", params.saSchedule, 0, 1)
  TUNE_INT("sacalibrate", params.saCalibrationIterations, 0, 1000000)
  TUNE_REAL("satarget", params.saCalibrationTarget, 1.0e-9, 1 - 1.0e-9)
  TUNE_REAL("saquantile", params.saCalibrationQuantile, 0, 1)
  TUNE_INT("noimp", params.maxNoImprove, 0, COIN_INT_MAX)
  TUNE_REAL("maxwork", params.maxWork, 0, 1.0e15)
  TUNE_REAL("totalwork", params.totalWorkLimit, 0, 1.0e15)
  TUNE_REAL("realstallwork", params.maxRealStallWork, 0, 1.0e15)
  TUNE_REAL("lcalpha", params.lnsColAlpha, -1, 10)
  TUNE_REAL("lralpha", params.lnsRowAlpha, -1, 10)
  TUNE_INT("lls", params.lnsLocalSearch, -1, 5)
  TUNE_INT("lrow", params.lnsRowRule, -1, 6)
  TUNE_INT("focus", params.lnsFocus, 0, 1)
  TUNE_REAL("artificial", artificial, 0, 1.0e10)
  TUNE_REAL("artgrow", params.artificialGrowth, 1, 1.0e6)
  TUNE_INT("artrounds", params.artificialRounds, 0, 10000)
  TUNE_REAL("artwork", params.artificialRoundWork, 0, 1.0e15)
#undef TUNE_REAL
#undef TUNE_INT
  return false;
}
}

CbcSppGrasp::CbcSppGrasp()
  : applicable_(false)
  , numberColumns_(0)
  , numberOriginalColumns_(0)
  , artificialFactor_(0.0)
  , numberRestrictingRows_(0)
  , numberRestrictedColumns_(0)
  , numberComplemented_(0)
  , numberRows_(0)
  , numberEqual_(0)
  , numberSlackRows_(0)
  , objOffset_(0.0)
  , costTolerance_(1.0e-9)
  , scoreRule_(-1)
  , scoreBeta_(0.0)
  , scoreLpWeight_(0.0)
  , scoreConflictWeight_(0.0)
  , scoresValid_(false)
  , seed_(1)
  , markStamp_(0)
  , work_(0.0)
{
}

bool CbcSppGrasp::load(const OsiSolverInterface &solver, int numberRows)
{
  applicable_ = false;
  scoresValid_ = false;
  lpX_.clear();
  lpDj_.clear();
  int n = solver.getNumCols();
  int m = solver.getNumRows();
  if (numberRows >= 0 && numberRows < m)
    m = numberRows;
  numberColumns_ = n;
  numberOriginalColumns_ = n;
  numberRows_ = m;
  numberEqual_ = 0;
  numberSlackRows_ = 0;
  if (!n || !m) {
    reason_ = "empty problem";
    return false;
  }
  const double *objective = solver.getObjCoefficients();
  const double direction = solver.getObjSense();
  const double inf = solver.getInfinity();
  const double tol = 1.0e-9;
  const CoinPackedMatrix *matrix = solver.getMatrixByCol();
  const double *element = matrix->getElements();
  const int *row = matrix->getIndices();
  const CoinBigIndex *start = matrix->getVectorStarts();
  const int *length = matrix->getVectorLengths();
  // binary columns with only -1 coefficients are complemented (y' = 1 - y):
  // preprocessing turns a singleton slack s of sum x + s = 1 into y = 1 - s,
  // leaving sum x - y = 0 (or <= 0)
  complemented_.assign(n, 0);
  numberComplemented_ = 0;
  std::vector< double > lowerT(solver.getColLower(), solver.getColLower() + n);
  std::vector< double > upperT(solver.getColUpper(), solver.getColUpper() + n);
  std::vector< double > rowLowerT(solver.getRowLower(), solver.getRowLower() + m);
  std::vector< double > rowUpperT(solver.getRowUpper(), solver.getRowUpper() + m);
  for (int j = 0; j < n; j++) {
    if (!solver.isInteger(j) || lowerT[j] < -tol || upperT[j] > 1.0 + tol)
      continue;
    int number = 0;
    bool allMinusOne = true;
    for (CoinBigIndex k = start[j]; k < start[j] + length[j]; k++) {
      if (row[k] >= m || fabs(element[k]) < 1.0e-12)
        continue;
      number++;
      if (fabs(element[k] + 1.0) > tol)
        allMinusOne = false;
    }
    if (!number || !allMinusOne)
      continue;
    complemented_[j] = 1;
    numberComplemented_++;
    double l = lowerT[j];
    lowerT[j] = 1.0 - upperT[j];
    upperT[j] = 1.0 - l;
    for (CoinBigIndex k = start[j]; k < start[j] + length[j]; k++) {
      int i = row[k];
      if (i >= m || fabs(element[k]) < 1.0e-12)
        continue;
      // -y = -1 + y'
      if (rowLowerT[i] > -inf)
        rowLowerT[i] += 1.0;
      if (rowUpperT[i] < inf)
        rowUpperT[i] += 1.0;
    }
  }
  const double *lower = lowerT.data();
  const double *upper = upperT.data();
  const double *rowLower = rowLowerT.data();
  const double *rowUpper = rowUpperT.data();
  // coefficient in the complemented space
  auto coef = [&](int j, double e) { return complemented_[j] ? -e : e; };

  isEqual_.assign(m, 0);
  std::vector< char > freeRow(m, 0);
  const CoinPackedMatrix *byRow = solver.getMatrixByRow();
  numberRestrictingRows_ = 0;
  numberRestrictedColumns_ = 0;
  for (int i = 0; i < m; i++) {
    if (rowUpper[i] >= inf && rowLower[i] <= -inf) {
      freeRow[i] = 1;
      continue;
    }
    // coefficients on columns not fixed at zero; activity of those at one
    bool ones = true;
    double fixedActivity = 0.0;
    for (CoinBigIndex k = byRow->getVectorFirst(i); k < byRow->getVectorLast(i); k++) {
      int j = byRow->getIndices()[k];
      double e = coef(j, byRow->getElements()[k]);
      if (upper[j] >= 0.5 && fabs(e) >= 1.0e-12 && fabs(e - 1.0) > tol)
        ones = false;
      if (lower[j] > 0.5)
        fixedActivity += e;
    }
    if (ones && fabs(rowUpper[i] - 1.0) <= tol && rowLower[i] <= 1.0 + tol) {
      if (rowLower[i] > 1.0e-6) {
        isEqual_[i] = 1;
        numberEqual_++;
      }
      continue;
    }
    if (rowLower[i] - fixedActivity <= tol && rowUpper[i] - fixedActivity >= -tol) {
      // satisfied by zero (e.g. x - y <= 0 left by preprocessing): its
      // columns that could violate it are not used
      freeRow[i] = 2;
      numberRestrictingRows_++;
      continue;
    }
    reason_ = ones ? "row with right-hand side other than 1" : "coefficient other than 1";
    return false;
  }
  if (!numberEqual_) {
    reason_ = "no partitioning rows";
    return false;
  }
  usable_.assign(n, 0);
  forced_.assign(n, 0);
  cost_.resize(n);
  double offset = 0.0;
  for (int j = 0; j < n; j++) {
    cost_[j] = direction * objective[j];
    if (complemented_[j]) {
      offset += cost_[j];
      cost_[j] = -cost_[j];
    }
    if (upper[j] < 0.5) {
      if (lower[j] < -tol) {
        reason_ = "column with negative lower bound";
        return false;
      }
      continue;
    }
    if (!solver.isInteger(j) || lower[j] < -tol || upper[j] > 1.0 + tol) {
      reason_ = "non-binary column";
      return false;
    }
    usable_[j] = 1;
    if (lower[j] > 0.5)
      forced_[j] = 1;
  }
  for (int j = 0; j < n; j++) {
    for (CoinBigIndex k = start[j]; k < start[j] + length[j]; k++) {
      int i = row[k];
      if (i >= m || freeRow[i] != 2 || fabs(element[k]) < 1.0e-12)
        continue;
      bool harmful = coef(j, element[k]) > 0.0 ? rowUpper[i] < inf : rowLower[i] > -inf;
      if (!harmful || forced_[j])
        continue;
      if (usable_[j]) {
        usable_[j] = 0;
        numberRestrictedColumns_++;
      }
    }
  }
  colStart_.assign(n + 1, 0);
  colRows_.clear();
  colRows_.reserve(matrix->getNumElements());
  std::vector< int > rowCount(m, 0);
  for (int j = 0; j < n; j++) {
    colStart_[j] = static_cast< int >(colRows_.size());
    for (CoinBigIndex k = start[j]; k < start[j] + length[j]; k++) {
      int i = row[k];
      if (i >= m || freeRow[i] || fabs(element[k]) < 1.0e-12)
        continue;
      if (fabs(coef(j, element[k]) - 1.0) > tol) {
        if (usable_[j]) {
          reason_ = "coefficient other than 1";
          return false;
        }
        continue;
      }
      if (!usable_[j])
        continue;
      colRows_.push_back(i);
      rowCount[i]++;
    }
  }
  colStart_[n] = static_cast< int >(colRows_.size());
  if (artificialFactor_ > 0.0) {
    // artificial singleton columns on partitioning rows without a singleton
    std::vector< char > hasSingleton(m, 0);
    std::vector< double > maxCost(m, -COIN_DBL_MAX);
    for (int j = 0; j < n; j++) {
      if (!usable_[j] || colStart_[j + 1] == colStart_[j])
        continue;
      if (colStart_[j + 1] - colStart_[j] == 1)
        hasSingleton[colRows_[colStart_[j]]] = 1;
      for (int k = colStart_[j]; k < colStart_[j + 1]; k++)
        maxCost[colRows_[k]] = std::max(maxCost[colRows_[k]], cost_[j]);
    }
    double maxAbs = 0.0;
    for (int j = 0; j < n; j++)
      if (usable_[j])
        maxAbs = std::max(maxAbs, fabs(cost_[j]));
    if (maxAbs == 0.0)
      maxAbs = 1.0;
    artificialBase_.clear();
    artificialUnit_.clear();
    for (int i = 0; i < m; i++) {
      if (!isEqual_[i] || hasSingleton[i] || !rowCount[i])
        continue;
      double c = maxCost[i];
      artificialBase_.push_back(c);
      artificialUnit_.push_back(std::max(fabs(c), maxAbs));
      cost_.push_back(c + (artificialFactor_ - 1.0) * artificialUnit_.back());
      usable_.push_back(1);
      forced_.push_back(0);
      colRows_.push_back(i);
      colStart_.push_back(static_cast< int >(colRows_.size()));
      rowCount[i]++;
    }
    numberColumns_ = static_cast< int >(cost_.size());
  }
  // from here on, n counts the artificial columns too
  n = numberColumns_;
  rowStart_.assign(m + 1, 0);
  for (int i = 0; i < m; i++)
    rowStart_[i + 1] = rowStart_[i] + rowCount[i];
  rowCols_.resize(colRows_.size());
  std::vector< int > fill(rowStart_.begin(), rowStart_.end() - 1);
  for (int j = 0; j < n; j++)
    for (int k = colStart_[j]; k < colStart_[j + 1]; k++)
      rowCols_[fill[colRows_[k]]++] = j;
  numberEqualInCol_.assign(n, 0);
  degree_.assign(n, 0);
  slackCol_.assign(m, -1);
  for (int j = 0; j < n; j++) {
    if (!usable_[j])
      continue;
    for (int k = colStart_[j]; k < colStart_[j + 1]; k++) {
      int i = colRows_[k];
      numberEqualInCol_[j] += isEqual_[i];
      degree_[j] += rowCount[i] - 1;
    }
    if (colStart_[j + 1] - colStart_[j] == 1) {
      int i = colRows_[colStart_[j]];
      if (isEqual_[i] && (slackCol_[i] < 0 || cost_[j] < cost_[slackCol_[i]]))
        slackCol_[i] = j;
    }
  }
  for (int i = 0; i < m; i++) {
    if (isEqual_[i] && slackCol_[i] >= 0)
      numberSlackRows_++;
    if (isEqual_[i] && !rowCount[i]) {
      reason_ = "partitioning row without columns";
      return false;
    }
  }
  sortRepairOrder();
  objOffset_ = offset;
  alive_.assign(n, 0);
  count_.assign(m, 0);
  pointer_.assign(2 * m, 0);
  owner_.assign(m, -1);
  mark_.assign(std::max(n, m), 0);
  markStamp_ = 0;
  applicable_ = true;
  reason_.clear();
  return true;
}

void CbcSppGrasp::sortRepairOrder()
{
  // move deltas sum a few costs: tolerate their rounding error, which
  // matters once escalated artificial penalties are huge
  double maxCost = 0.0;
  for (int j = 0; j < numberColumns_; j++)
    maxCost = std::max(maxCost, fabs(cost_[j]));
  costTolerance_ = std::max(1.0e-9, 1.0e-14 * maxCost);
  // columns of each partitioning row by cost per partitioning row
  repairOrder_ = rowCols_;
  minCostPerRow_.assign(numberRows_, 0.0);
  auto perRow = [&](int c) { return cost_[c] / numberEqualInCol_[c]; };
  for (int i = 0; i < numberRows_; i++) {
    if (!isEqual_[i])
      continue;
    std::sort(repairOrder_.begin() + rowStart_[i], repairOrder_.begin() + rowStart_[i + 1],
      [&](int a, int b) {
        double pa = perRow(a), pb = perRow(b);
        return pa != pb ? pa < pb : a < b;
      });
    minCostPerRow_[i] = perRow(repairOrder_[rowStart_[i]]);
  }
}

void CbcSppGrasp::scaleArtificialPenalty(double scale)
{
  if (numberColumns_ == numberOriginalColumns_)
    return;
  for (int k = 0; k < numberColumns_ - numberOriginalColumns_; k++)
    cost_[numberOriginalColumns_ + k] = artificialBase_[k]
      + scale * (artificialFactor_ - 1.0) * artificialUnit_[k];
  sortRepairOrder();
  scoresValid_ = false;
}

void CbcSppGrasp::toOriginal(const std::vector< char > &solution, double *x) const
{
  for (int j = 0; j < numberOriginalColumns_; j++)
    x[j] = (solution[j] != 0) != (complemented_[j] != 0) ? 1.0 : 0.0;
}

void CbcSppGrasp::fromOriginal(const double *x, std::vector< char > &solution) const
{
  solution.assign(numberColumns_, 0);
  for (int j = 0; j < numberOriginalColumns_; j++)
    solution[j] = (x[j] > 0.5) != (complemented_[j] != 0) ? 1 : 0;
}

int CbcSppGrasp::numberArtificialUsed(const std::vector< char > &solution) const
{
  int number = 0;
  for (int j = numberOriginalColumns_; j < numberColumns_; j++)
    number += solution[j] ? 1 : 0;
  return number;
}

CbcSppGrasp::FeatureSummary::FeatureSummary()
  : partitionRows(0)
  , nativeSingletonRows(0)
  , nativeSingletonMissing(0)
  , restrictedRows(0)
  , restrictedColumns(0)
  , complementedColumns(0)
  , forcedColumns(0)
  , lpAvailable(0)
  , rowOptionsMin(0.0)
  , rowOptionsQ25(0.0)
  , rowOptionsMedian(0.0)
  , rowOptionsQ75(0.0)
  , rowOptionsMean(0.0)
  , columnDegreeQ25(0.0)
  , columnDegreeMedian(0.0)
  , columnDegreeQ75(0.0)
  , columnDegreeMean(0.0)
  , singletonPremiumQ25(0.0)
  , singletonPremiumMedian(0.0)
  , singletonPremiumQ75(0.0)
  , restrictedRowFraction(0.0)
  , restrictedColumnFraction(0.0)
  , complementedFraction(0.0)
  , forcedFraction(0.0)
  , lpFractionalityMean(0.0)
  , nativeSingletonLpMass(0.0)
{
}

void CbcSppGrasp::featureSummary(FeatureSummary &features) const
{
  features = FeatureSummary();
  features.partitionRows = numberEqual_;
  features.restrictedRows = numberRestrictingRows_;
  features.restrictedColumns = numberRestrictedColumns_;
  features.complementedColumns = numberComplemented_;
  if (!numberOriginalColumns_ || !numberRows_)
    return;
  std::vector< double > rowOptions;
  std::vector< double > columnDegree;
  std::vector< double > singletonPremium;
  rowOptions.reserve(numberEqual_);
  columnDegree.reserve(numberOriginalColumns_);
  singletonPremium.reserve(numberEqual_);
  double optionSum = 0.0;
  double degreeSum = 0.0;
  for (int i = 0; i < numberRows_; i++) {
    if (!isEqual_[i])
      continue;
    int count = 0;
    int bestSingleton = slackCol_[i] >= 0 && slackCol_[i] < numberOriginalColumns_
      ? slackCol_[i]
      : -1;
    if (bestSingleton >= 0)
      features.nativeSingletonRows++;
    else
      features.nativeSingletonMissing++;
    int bestMulti = -1;
    for (int k = rowStart_[i]; k < rowStart_[i + 1]; k++) {
      int j = rowCols_[k];
      if (j >= numberOriginalColumns_)
        continue;
      count++;
      if (colStart_[j + 1] - colStart_[j] > 1
        && (bestMulti < 0 || cost_[j] < cost_[bestMulti]))
        bestMulti = j;
    }
    rowOptions.push_back(count);
    optionSum += count;
    if (bestSingleton >= 0 && bestMulti >= 0) {
      const double scale = std::max(1.0,
        std::max(fabs(cost_[bestSingleton]), fabs(cost_[bestMulti])));
      singletonPremium.push_back(cost_[bestMulti] / scale
        - cost_[bestSingleton] / scale);
    }
  }
  for (int j = 0; j < numberOriginalColumns_; j++) {
    if (forced_[j])
      features.forcedColumns++;
    if (!usable_[j])
      continue;
    const int degree = colStart_[j + 1] - colStart_[j];
    columnDegree.push_back(degree);
    degreeSum += degree;
  }
  auto quantile = [](std::vector< double > values, double p) {
    if (values.empty())
      return 0.0;
    const double rank = p * (values.size() - 1);
    const size_t low = static_cast< size_t >(rank);
    const size_t high = std::min(low + 1, values.size() - 1);
    std::nth_element(values.begin(), values.begin() + low, values.end());
    const double lower = values[low];
    if (high == low)
      return lower;
    std::nth_element(values.begin() + low + 1, values.begin() + high, values.end());
    return lower + (rank - low) * (values[high] - lower);
  };
  if (!rowOptions.empty()) {
    features.rowOptionsMin = *std::min_element(rowOptions.begin(), rowOptions.end());
    features.rowOptionsQ25 = quantile(rowOptions, 0.25);
    features.rowOptionsMedian = quantile(rowOptions, 0.5);
    features.rowOptionsQ75 = quantile(rowOptions, 0.75);
    features.rowOptionsMean = optionSum / rowOptions.size();
  }
  if (!columnDegree.empty()) {
    features.columnDegreeQ25 = quantile(columnDegree, 0.25);
    features.columnDegreeMedian = quantile(columnDegree, 0.5);
    features.columnDegreeQ75 = quantile(columnDegree, 0.75);
    features.columnDegreeMean = degreeSum / columnDegree.size();
  }
  if (!singletonPremium.empty()) {
    features.singletonPremiumQ25 = quantile(singletonPremium, 0.25);
    features.singletonPremiumMedian = quantile(singletonPremium, 0.5);
    features.singletonPremiumQ75 = quantile(singletonPremium, 0.75);
  }
  features.restrictedRowFraction = static_cast< double >(numberRestrictingRows_) / numberRows_;
  features.restrictedColumnFraction = static_cast< double >(numberRestrictedColumns_) / numberOriginalColumns_;
  features.complementedFraction = static_cast< double >(numberComplemented_) / numberOriginalColumns_;
  features.forcedFraction = static_cast< double >(features.forcedColumns) / numberOriginalColumns_;
  if (lpX_.size() >= static_cast< size_t >(numberOriginalColumns_)) {
    double fractionality = 0.0;
    int singletonCount = 0;
    for (int j = 0; j < numberOriginalColumns_; j++) {
      const double x = std::min(1.0, std::max(0.0, lpX_[j]));
      fractionality += x * (1.0 - x);
    }
    for (int i = 0; i < numberRows_; i++)
      if (isEqual_[i] && slackCol_[i] >= 0 && slackCol_[i] < numberOriginalColumns_) {
        features.nativeSingletonLpMass += lpX_[slackCol_[i]];
        singletonCount++;
      }
    features.lpAvailable = 1;
    features.lpFractionalityMean = fractionality / numberOriginalColumns_;
    if (singletonCount)
      features.nativeSingletonLpMass /= singletonCount;
  }
}

bool CbcSppGrasp::meetsNativeCoverage(double minimum, std::string *description) const
{
  if (!applicable_) {
    if (description)
      *description = "native coverage unavailable: skip";
    return false;
  }
  int covered = 0;
  int forced = 0;
  for (int i = 0; i < numberRows_; i++)
    if (isEqual_[i] && slackCol_[i] >= 0 && slackCol_[i] < numberOriginalColumns_)
      covered++;
  for (int j = 0; j < numberOriginalColumns_; j++)
    forced += forced_[j] != 0;
  const double coverage = numberEqual_ > 0
    ? static_cast< double >(covered) / numberEqual_
    : 0.0;
  const bool accepted = minimum <= 0.0
    || (applicable_ && numberEqual_ > 0 && !forced && coverage + 1.0e-12 >= minimum);
  if (description) {
    char line[256];
    sprintf(line, "native coverage %d/%d (%.8g), forced %d, minimum %.8g: %s",
      covered, numberEqual_, coverage, forced, minimum, accepted ? "run" : "skip");
    *description = line;
  }
  return accepted;
}

void CbcSppGrasp::setLpSolution(const double *x, const double *reducedCost)
{
  scoresValid_ = false;
  // artificial columns: at 0 in the LP, reduced cost = their cost
  if (x) {
    lpX_.assign(x, x + numberOriginalColumns_);
    for (int j = 0; j < numberOriginalColumns_; j++)
      if (complemented_[j])
        lpX_[j] = 1.0 - lpX_[j];
    lpX_.resize(numberColumns_, 0.0);
  } else {
    lpX_.clear();
  }
  if (reducedCost) {
    lpDj_.assign(reducedCost, reducedCost + numberOriginalColumns_);
    for (int j = 0; j < numberOriginalColumns_; j++)
      if (complemented_[j])
        lpDj_[j] = -lpDj_[j];
    for (int j = numberOriginalColumns_; j < numberColumns_; j++)
      lpDj_.push_back(cost_[j]);
  } else {
    lpDj_.clear();
  }
}

void CbcSppGrasp::setLpGuide(const double *x)
{
  scoresValid_ = false;
  if (!x) {
    lpX_.clear();
    return;
  }
  lpX_.assign(x, x + numberOriginalColumns_);
  for (int j = 0; j < numberOriginalColumns_; j++)
    if (complemented_[j])
      lpX_[j] = 1.0 - lpX_[j];
  lpX_.resize(numberColumns_, 0.0);
}

int CbcSppGrasp::randomInt(int n)
{
  // xorshift32
  seed_ ^= seed_ << 13;
  seed_ ^= seed_ >> 17;
  seed_ ^= seed_ << 5;
  return static_cast< int >(seed_ % static_cast< unsigned int >(n));
}

void CbcSppGrasp::computeScores(const Config &config)
{
  int colRule = config.colRule;
  if (colRule == ColReducedCost && lpDj_.empty())
    colRule = ColCostPerRow;
  double lpWeight = lpX_.empty() ? 0.0 : config.lpWeight;
  if (scoresValid_ && scoreRule_ == colRule && scoreBeta_ == config.beta
    && scoreLpWeight_ == lpWeight && scoreConflictWeight_ == config.conflictWeight)
    return;
  scoreRule_ = colRule;
  scoreBeta_ = config.beta;
  scoreLpWeight_ = lpWeight;
  scoreConflictWeight_ = config.conflictWeight;
  scoresValid_ = true;
  const int n = numberColumns_;
  score_.assign(n, kNoOption);
  std::vector< double > perRow(n, kNoOption);
  double sumAbs = 0.0, sumDegree = 0.0;
  int number = 0;
  for (int j = 0; j < n; j++) {
    int k = numberEqualInCol_[j];
    if (!usable_[j] || !k)
      continue;
    perRow[j] = cost_[j] / (config.beta == 1.0 ? k : pow(static_cast< double >(k), config.beta));
    if (colRule == ColReducedCost) {
      score_[j] = lpDj_[j];
    } else if (colRule == ColProfit) {
      double saved = 0.0;
      for (int p = colStart_[j]; p < colStart_[j + 1]; p++) {
        int i = colRows_[p];
        if (isEqual_[i] && slackCol_[i] >= 0)
          saved += cost_[slackCol_[i]];
      }
      score_[j] = (cost_[j] - saved) / (config.beta == 1.0 ? k : pow(static_cast< double >(k), config.beta));
    } else {
      score_[j] = perRow[j];
    }
    sumAbs += fabs(score_[j]);
    sumDegree += degree_[j];
    number++;
  }
  double scale = number ? sumAbs / number : 1.0;
  if (scale < 1.0e-9)
    scale = 1.0;
  double meanDegree = (number && sumDegree > 0.0) ? sumDegree / number : 1.0;
  if (lpWeight || config.conflictWeight) {
    for (int j = 0; j < n; j++) {
      if (score_[j] == kNoOption)
        continue;
      if (lpWeight)
        score_[j] -= lpWeight * lpX_[j] * scale;
      if (config.conflictWeight)
        score_[j] += config.conflictWeight * (degree_[j] / meanDegree) * scale;
    }
  }
  sortedRowCols_ = rowCols_;
  for (int i = 0; i < numberRows_; i++) {
    if (!isEqual_[i])
      continue;
    std::sort(sortedRowCols_.begin() + rowStart_[i], sortedRowCols_.begin() + rowStart_[i + 1],
      [&](int a, int b) {
        if (score_[a] != score_[b])
          return score_[a] < score_[b];
        if (perRow[a] != perRow[b])
          return perRow[a] < perRow[b];
        return a < b;
      });
  }
}

double CbcSppGrasp::evaluate(const std::vector< char > &solution,
  bool accountWork) const
{
  std::vector< int > cover(numberRows_, 0);
  double value = objOffset_;
  if (accountWork)
    work_ += numberColumns_ + numberRows_;
  for (int j = 0; j < numberColumns_; j++) {
    if (!solution[j]) {
      if (forced_[j])
        return COIN_DBL_MAX;
      continue;
    }
    if (!usable_[j])
      return COIN_DBL_MAX;
    value += cost_[j];
    if (accountWork)
      work_ += colStart_[j + 1] - colStart_[j];
    for (int k = colStart_[j]; k < colStart_[j + 1]; k++)
      cover[colRows_[k]]++;
  }
  for (int i = 0; i < numberRows_; i++) {
    if (cover[i] > 1 || (isEqual_[i] && cover[i] != 1))
      return COIN_DBL_MAX;
  }
  return value;
}

int CbcSppGrasp::numberSlacksUsed(const std::vector< char > &solution) const
{
  int number = 0;
  for (int i = 0; i < numberRows_; i++)
    if (slackCol_[i] >= 0 && solution[slackCol_[i]])
      number++;
  return number;
}

double CbcSppGrasp::trivialSolution(std::vector< char > &solution)
{
  solution.assign(numberColumns_, 0);
  std::fill(owner_.begin(), owner_.end(), -1);
  for (int j = 0; j < numberColumns_; j++) {
    if (!forced_[j])
      continue;
    solution[j] = 1;
    for (int k = colStart_[j]; k < colStart_[j + 1]; k++)
      owner_[colRows_[k]] = j;
  }
  for (int i = 0; i < numberRows_; i++) {
    if (isEqual_[i] && owner_[i] < 0) {
      if (slackCol_[i] < 0)
        return COIN_DBL_MAX;
      solution[slackCol_[i]] = 1;
      owner_[i] = slackCol_[i];
    }
  }
  return evaluate(solution);
}

double CbcSppGrasp::construct(const Config &config, std::vector< char > &solution,
  const std::vector< char > *partial, const std::vector< char > *previous)
{
  if (!applicable_)
    return COIN_DBL_MAX;
  computeScores(config);
  const int n = numberColumns_;
  const int m = numberRows_;
  // lazy: a column is alive if usable and none of its rows is covered,
  // checked on demand; otherwise selections kill conflicting columns
  // eagerly, which also maintains count_ (needed by the options rules)
  const bool lazy = config.rowRule != RowMinOptions && config.rowRule != RowOptionsOverRegret
    && config.rowRule != RowLogOptionsOverRegret;
  solution.assign(n, 0);
  // basic operations, added to work_ on exit
  struct WorkCount {
    double &total;
    double w;
    ~WorkCount() { total += w; }
  } work = { work_, static_cast< double >(n + m) };
  if (!lazy) {
    for (int j = 0; j < n; j++)
      alive_[j] = usable_[j];
  }
  for (int i = 0; i < m; i++) {
    owner_[i] = -1;
    pointer_[i] = rowStart_[i];
    pointer_[m + i] = rowStart_[i];
    if (!lazy)
      count_[i] = rowStart_[i + 1] - rowStart_[i];
  }
  auto isAlive = [&](int c) {
    if (!lazy)
      return alive_[c] != 0;
    if (!usable_[c])
      return false;
    work.w++;
    for (int k = colStart_[c]; k < colStart_[c + 1]; k++) {
      if (owner_[colRows_[k]] >= 0)
        return false;
      work.w++;
    }
    return true;
  };
  auto kill = [&](int c) {
    alive_[c] = 0;
    work.w += colStart_[c + 1] - colStart_[c];
    for (int k = colStart_[c]; k < colStart_[c + 1]; k++)
      count_[colRows_[k]]--;
  };
  auto select = [&](int c) {
#ifndef NDEBUG
    const double saved = work.w; // keep work identical in debug builds
    assert(isAlive(c));
    work.w = saved;
#endif
    solution[c] = 1;
    work.w += colStart_[c + 1] - colStart_[c];
    for (int k = colStart_[c]; k < colStart_[c + 1]; k++) {
      int i = colRows_[k];
      owner_[i] = c;
      if (lazy)
        continue;
      work.w += rowStart_[i + 1] - rowStart_[i];
      for (int l = rowStart_[i]; l < rowStart_[i + 1]; l++) {
        int c2 = rowCols_[l];
        if (alive_[c2])
          kill(c2);
      }
    }
  };
  for (int j = 0; j < n; j++) {
    if (forced_[j]) {
      if (!isAlive(j))
        return COIN_DBL_MAX;
      select(j);
    }
  }
  if (partial) {
    for (int j = 0; j < n; j++) {
      if ((*partial)[j] && !solution[j]) {
        if (!isAlive(j))
          return COIN_DBL_MAX;
        select(j);
      }
    }
  }
  std::vector< int > open;
  open.reserve(numberEqual_);
  for (int i = 0; i < m; i++)
    if (isEqual_[i] && owner_[i] < 0)
      open.push_back(i);
  std::vector< double > key(open.size());
  // first and second alive option of row i (by score), -1 if none
  auto firstAlive = [&](int i) {
    int &p = pointer_[i];
    while (p < rowStart_[i + 1] && !isAlive(sortedRowCols_[p]))
      p++;
    return p < rowStart_[i + 1] ? sortedRowCols_[p] : -1;
  };
  auto secondAlive = [&](int i) {
    int &p = pointer_[m + i];
    if (p <= pointer_[i])
      p = pointer_[i] + 1;
    while (p < rowStart_[i + 1] && !isAlive(sortedRowCols_[p]))
      p++;
    return p < rowStart_[i + 1] ? sortedRowCols_[p] : -1;
  };
  while (!open.empty()) {
    // drop covered rows, settle rows whose best option is a singleton
    bool settled = true;
    while (settled) {
      settled = false;
      work.w += open.size();
      size_t nOpen = 0;
      for (size_t t = 0; t < open.size(); t++) {
        int i = open[t];
        if (owner_[i] >= 0)
          continue;
        int c = firstAlive(i);
        if (c < 0)
          return COIN_DBL_MAX; // dead end
        if (colStart_[c + 1] - colStart_[c] == 1) {
          select(c);
          settled = true;
          continue;
        }
        open[nOpen++] = i;
      }
      open.resize(nOpen);
    }
    if (open.empty())
      break;
    double bestKey = COIN_DBL_MAX;
    key.resize(open.size());
    work.w += 2 * open.size();
    for (size_t t = 0; t < open.size(); t++) {
      int i = open[t];
      int c1 = firstAlive(i);
      int c2 = secondAlive(i);
      double s1 = score_[c1];
      double regret = (c2 >= 0) ? score_[c2] - s1 : kNoOption;
      double options = count_[i];
      double k;
      if (c2 < 0) {
        k = -COIN_DBL_MAX; // only one option left: take it now
      } else {
        switch (config.rowRule) {
        case RowMinOptions:
          k = options;
          break;
        case RowMaxRegret:
          k = -regret;
          break;
        case RowOptionsOverRegret:
          k = options / (regret + 1.0e-6);
          break;
        case RowLogOptionsOverRegret:
          k = log(1.0 + options) / (regret + 1.0e-6);
          break;
        case RowMinBest:
          k = s1;
          break;
        case RowMaxBest:
          k = -s1;
          break;
        case RowRandom:
        default:
          k = 0.0;
          break;
        }
      }
      key[t] = k;
      bestKey = std::min(bestKey, k);
    }
    double threshold = bestKey + config.rowAlpha * fabs(bestKey) + 1.0e-12;
    int chosen = -1, numberCandidates = 0;
    for (size_t t = 0; t < open.size(); t++) {
      if (key[t] <= threshold) {
        numberCandidates++;
        if (numberCandidates == 1 || !randomInt(numberCandidates))
          chosen = static_cast< int >(t);
      }
    }
    int i = open[chosen];
    int c1 = firstAlive(i);
    int col = c1;
    if (config.colAlpha > 0.0 && config.colRclMax > 1) {
      double s1 = score_[c1];
      double limit = s1 + config.colAlpha * std::max(fabs(s1), 1.0e-6) + 1.0e-12;
      numberCandidates = 0;
      for (int p = pointer_[i]; p < rowStart_[i + 1] && numberCandidates < config.colRclMax; p++) {
        int c = sortedRowCols_[p];
        if (!isAlive(c))
          continue;
        if (score_[c] > limit)
          break;
        numberCandidates++;
        if (numberCandidates == 1 || !randomInt(numberCandidates))
          col = c;
      }
    }
    select(col);
    open[chosen] = open.back();
    open.pop_back();
  }
  // packing rows only: add any profitable column still possible
  for (int j = 0; j < n; j++)
    if (cost_[j] < 0.0 && !solution[j] && isAlive(j))
      select(j);
  if (config.localSearch)
    return localSearch(solution, config.localSearch, partial, config.repairScan, previous);
  return evaluate(solution);
}

double CbcSppGrasp::backtrackingRepair(const Config &config, const RunParams &params,
  std::vector< char > &solution, const std::vector< char > &partial,
  const std::vector< char > *previous, double bestValue, RepairStats &stats)
{
  int bestArtificial = bestValue == COIN_DBL_MAX ? COIN_INT_MAX : numberArtificialUsed(solution);
  const double greedyValue = bestValue;
  const int greedyArtificial = bestArtificial;
  const double startWork = work_;
  stats.calls++;
  stats.greedyArtificial += greedyArtificial > 0;
  int nodes = 0;
  long long leaves = 0;
  auto finish = [&]() {
    stats.nodes += nodes;
    stats.work += work_ - startWork;
    stats.noLeaf += leaves == 0;
    stats.improved += bestArtificial < greedyArtificial
      || (bestArtificial == greedyArtificial && bestValue < greedyValue);
    stats.feasibilityGains += greedyArtificial > 0 && bestArtificial == 0;
    return bestValue;
  };
  std::vector< int > owners(numberRows_, -1);
  std::vector< char > trial(numberColumns_, 0);
  auto fits = [&](int j) {
    work_++;
    for (int k = colStart_[j]; k < colStart_[j + 1]; k++) {
      work_++;
      if (owners[colRows_[k]] >= 0)
        return false;
    }
    return true;
  };
  auto set = [&](int j, bool selected) {
    trial[j] = selected;
    for (int k = colStart_[j]; k < colStart_[j + 1]; k++) {
      work_++;
      owners[colRows_[k]] = selected ? j : -1;
    }
  };
  for (int j = 0; j < numberColumns_; j++) {
    work_++;
    if (forced_[j] || partial[j]) {
      if (!usable_[j] || !fits(j)) {
        stats.invalidPartial++;
        return finish();
      }
      set(j, true);
    }
  }
  std::vector< int > openRows;
  for (int i = 0; i < numberRows_; i++) {
    work_++;
    if (isEqual_[i] && owners[i] < 0)
      openRows.push_back(i);
  }
  stats.openRows += openRows.size();
  stats.maxOpenRows = std::max(stats.maxOpenRows, static_cast< long long >(openRows.size()));
  struct Frame {
    std::vector< int > options;
    size_t next;
    int selected;
  };
  std::vector< Frame > stack;
  bool descend = true;
  while (nodes < params.backtrackNodes && work_ - startWork < params.backtrackWork) {
    if (descend) {
      int openRow = -1;
      std::vector< int > options;
      std::vector< int > candidates;
      const int width = static_cast< int >(stack.size()) < params.backtrackDepth ? params.backtrackWidth : 1;
      // Scarce-option rows first; candidates retain the construction score order.
      for (int i : openRows) {
        work_++;
        if (owners[i] >= 0)
          continue;
        candidates.clear();
        const int passes = params.backtrackRealFirst ? 2 : 1;
        for (int pass = 0; pass < passes && static_cast< int >(candidates.size()) < width + 1; pass++) {
          for (int p = rowStart_[i]; p < rowStart_[i + 1]; p++) {
            if (work_ - startWork >= params.backtrackWork) {
              stats.workExits++;
              return finish();
            }
            int j = sortedRowCols_[p];
            if (params.backtrackRealFirst) {
              work_++;
              if ((j >= numberOriginalColumns_) != (pass == 1))
                continue;
            }
            if (usable_[j] && fits(j))
              candidates.push_back(j);
            if (static_cast< int >(candidates.size()) >= width + 1)
              break;
          }
        }
        if (openRow < 0 || candidates.size() < options.size()) {
          openRow = i;
          options.swap(candidates);
        }
        if (options.empty() || (params.backtrackRealFirst && options.size() == 1))
          break;
      }
      if (openRow < 0) {
        // Complete profitable packing-only columns without changing the DFS state.
        std::vector< int > added;
        for (int j = 0; j < numberColumns_; j++) {
          work_++;
          if (work_ - startWork >= params.backtrackWork) {
            stats.workExits++;
            return finish();
          }
          if (usable_[j] && cost_[j] < 0.0 && !trial[j] && fits(j)) {
            set(j, true);
            added.push_back(j);
          }
        }
        std::vector< char > completed = trial;
        for (int j : added)
          set(j, false);
        work_ += numberColumns_ + numberRows_ + colRows_.size();
        double value = evaluate(completed);
        int artificial = numberArtificialUsed(completed);
        leaves++;
        stats.leaves++;
        stats.realLeaves += value != COIN_DBL_MAX && artificial == 0;
        if (value != COIN_DBL_MAX && (artificial < bestArtificial || (artificial == bestArtificial && value < bestValue))) {
          if (config.localSearch)
            value = localSearch(completed, config.localSearch, &partial, config.repairScan, previous);
          artificial = numberArtificialUsed(completed);
          if (artificial < bestArtificial || (artificial == bestArtificial && value < bestValue)) {
            bestArtificial = artificial;
            bestValue = value;
            solution.swap(completed);
          }
        }
      } else if (!options.empty()) {
        if (static_cast< int >(options.size()) > width)
          options.resize(width);
        stack.push_back(Frame { options, 0, -1 });
        stats.maxDepth = std::max(stats.maxDepth, static_cast< long long >(stack.size()));
      } else {
        stats.deadEnds++;
      }
      descend = false;
    }
    while (!stack.empty()) {
      Frame &frame = stack.back();
      if (frame.selected >= 0) {
        set(frame.selected, false);
        frame.selected = -1;
      }
      if (frame.next < frame.options.size()) {
        frame.selected = frame.options[frame.next++];
        set(frame.selected, true);
        nodes++;
        descend = true;
        break;
      }
      stack.pop_back();
    }
    if (stack.empty()) {
      stats.exhausted++;
      return finish();
    }
  }
  if (work_ - startWork >= params.backtrackWork)
    stats.workExits++;
  else
    stats.nodeExits++;
  return finish();
}

double CbcSppGrasp::localSearch(std::vector< char > &solution, int level,
  const std::vector< char > *reference, int repairScan,
  const std::vector< char > *previous)
{
  const int n = numberColumns_;
  const int m = numberRows_;
  // basic operations, added to work_ at the end
  double w = n + m;
  std::fill(owner_.begin(), owner_.end(), -1);
  for (int j = 0; j < n; j++) {
    if (!solution[j])
      continue;
    w += colStart_[j + 1] - colStart_[j];
    for (int k = colStart_[j]; k < colStart_[j + 1]; k++)
      owner_[colRows_[k]] = j;
  }
  // worklist of candidate columns to insert; a column is re-examined only
  // when the owner of one of its rows changed
  std::vector< int > queue;
  queue.reserve(n);
  std::vector< char > queued(n, 0);
  auto push = [&](int j) {
    if (!queued[j] && !solution[j] && usable_[j] && !forced_[j]) {
      queued[j] = 1;
      queue.push_back(j);
    }
  };
  if (previous) {
    std::vector< int > previousOwner(m, -1);
    w += n + m;
    for (int j = 0; j < n; j++) {
      if (!(*previous)[j])
        continue;
      w += colStart_[j + 1] - colStart_[j];
      for (int k = colStart_[j]; k < colStart_[j + 1]; k++)
        previousOwner[colRows_[k]] = j;
    }
    for (int i = 0; i < m; i++) {
      if (owner_[i] == previousOwner[i])
        continue;
      w += rowStart_[i + 1] - rowStart_[i];
      for (int l = rowStart_[i]; l < rowStart_[i + 1]; l++)
        push(rowCols_[l]);
    }
  } else if (reference) {
    // only rows not covered by a column of reference changed
    for (int i = 0; i < m; i++) {
      int o = owner_[i];
      if (o >= 0 && (*reference)[o])
        continue;
      w += rowStart_[i + 1] - rowStart_[i];
      for (int l = rowStart_[i]; l < rowStart_[i + 1]; l++)
        push(rowCols_[l]);
    }
  } else {
    w += n;
    for (int j = 0; j < n; j++)
      if (colStart_[j + 1] > colStart_[j])
        push(j);
  }
  w += queue.size();
  for (size_t t = queue.size(); t > 1; t--)
    std::swap(queue[t - 1], queue[randomInt(static_cast< int >(t))]);
  // tentative changes: (row, previous owner)
  std::vector< std::pair< int, int > > changed;
  std::vector< int > ejected, freed, added;
  const int tentative = n; // owner value of rows taken by the move
  auto setOwner = [&](int i, int o) {
    changed.push_back(std::make_pair(i, owner_[i]));
    owner_[i] = o;
  };
  {
    size_t head = 0;
    while (head < queue.size()) {
      int j = queue[head++];
      queued[j] = 0;
      if (head > 1000000 && head * 2 > queue.size()) {
        queue.erase(queue.begin(), queue.begin() + head);
        head = 0;
      }
      if (solution[j])
        continue;
      if (++markStamp_ == 0) {
        std::fill(mark_.begin(), mark_.end(), 0);
        markStamp_ = 1;
      }
      double delta = cost_[j];
      ejected.clear();
      bool feasible = true;
      w += 1 + colStart_[j + 1] - colStart_[j];
      for (int k = colStart_[j]; k < colStart_[j + 1]; k++) {
        int o = owner_[colRows_[k]];
        if (o >= 0 && mark_[o] != markStamp_) {
          if (forced_[o]) {
            feasible = false;
            break;
          }
          mark_[o] = markStamp_;
          ejected.push_back(o);
          delta -= cost_[o];
        }
      }
      if (!feasible || delta >= -costTolerance_)
        continue;
      // rows of ejected columns left open by j; bound repair cost from below
      changed.clear();
      for (int k = colStart_[j]; k < colStart_[j + 1]; k++)
        setOwner(colRows_[k], tentative);
      freed.clear();
      double bound = delta;
      for (size_t e = 0; e < ejected.size(); e++) {
        int o = ejected[e];
        w += colStart_[o + 1] - colStart_[o];
        for (int k = colStart_[o]; k < colStart_[o + 1]; k++) {
          int i = colRows_[k];
          if (owner_[i] == o) {
            setOwner(i, -1);
            if (isEqual_[i]) {
              freed.push_back(i);
              if (level < 2) {
                if (slackCol_[i] < 0)
                  feasible = false;
                else
                  bound += cost_[slackCol_[i]];
              } else {
                bound += minCostPerRow_[i];
              }
            }
          }
        }
      }
      added.clear();
      if (feasible && bound < -costTolerance_) {
        if (level < 2) {
          delta = bound;
          for (size_t f = 0; f < freed.size(); f++) {
            added.push_back(slackCol_[freed[f]]);
            setOwner(freed[f], tentative);
          }
        } else {
          // greedy repair: cheapest (per row) column whose rows are all open;
          // abandoned once delta plus a lower bound on covering the rows
          // still open cannot be negative (the move would be rejected)
          double remaining = bound - delta;
          for (size_t f = 0; f < freed.size() && feasible; f++) {
            int i = freed[f];
            if (owner_[i] != -1)
              continue;
            int chosen = -1, scanned = 0;
            for (int p = rowStart_[i]; p < rowStart_[i + 1] && scanned < repairScan; p++) {
              int c = repairOrder_[p];
              scanned++;
              bool free = true;
              w++;
              for (int k = colStart_[c]; k < colStart_[c + 1]; k++) {
                if (owner_[colRows_[k]] != -1) {
                  free = false;
                  break;
                }
                w++;
              }
              if (free) {
                chosen = c;
                break;
              }
            }
            if (chosen < 0)
              chosen = slackCol_[i];
            if (chosen < 0) {
              feasible = false;
              break;
            }
            added.push_back(chosen);
            delta += cost_[chosen];
            w += colStart_[chosen + 1] - colStart_[chosen];
            for (int k = colStart_[chosen]; k < colStart_[chosen + 1]; k++) {
              int r = colRows_[k];
              if (isEqual_[r] && owner_[r] == -1)
                remaining -= minCostPerRow_[r];
              setOwner(r, tentative);
            }
            if (delta + remaining > std::max(1.0e-7, costTolerance_))
              feasible = false;
          }
        }
      } else {
        feasible = false;
      }
      w += changed.size();
      if (!feasible || delta >= -costTolerance_) {
        for (size_t c = changed.size(); c > 0; c--)
          owner_[changed[c - 1].first] = changed[c - 1].second;
        continue;
      }
      // commit
      for (size_t e = 0; e < ejected.size(); e++)
        solution[ejected[e]] = 0;
      added.push_back(j);
      for (size_t a = 0; a < added.size(); a++) {
        int c = added[a];
        solution[c] = 1;
        for (int k = colStart_[c]; k < colStart_[c + 1]; k++)
          owner_[colRows_[k]] = c;
      }
      for (size_t c = 0; c < changed.size(); c++) {
        int i = changed[c].first;
        w += rowStart_[i + 1] - rowStart_[i];
        for (int l = rowStart_[i]; l < rowStart_[i + 1]; l++)
          push(rowCols_[l]);
      }
    }
  }
  work_ += w;
  return evaluate(solution);
}

double CbcSppGrasp::run(const Config &config, const RunParams &params,
  std::vector< char > &best, RunStats *stats, const std::vector< char > *initial)
{
  const double start = CoinGetTimeOfDay();
  const double startWork = work_;
  const int n = numberColumns_;
  RunStats local;
  RunStats &st = stats ? *stats : local;
  st = RunStats();
  // with artificial columns the search keeps the best solution overall
  // (best) but reports only solutions without artificials (realBest)
  const bool artificial = numberColumns_ > numberOriginalColumns_;
  std::vector< char > realBest;
  double realBestValue = COIN_DBL_MAX;
  double realBestWork = -1.0;
  double bestValue = trivialSolution(best);
  st.trivial = bestValue;
  int penaltyRound = 0;
  std::vector< char > current, trial;
  double currentValue = COIN_DBL_MAX;
  auto elapsed = [&]() { return CoinGetTimeOfDay() - start; };
  auto addProgress = [&](int kind, double value, int artificialUsed) {
    ProgressEvent event = { st.lnsIterations, kind, penaltyRound, artificialUsed,
      value, elapsed(), work_ - startWork };
    st.progress.push_back(event);
  };
  if (bestValue != COIN_DBL_MAX)
    addProgress(1, bestValue, numberArtificialUsed(best));
  auto record = [&](double value, std::vector< char > &sol) {
    if (value < bestValue - 1.0e-7) {
      bestValue = value;
      best = sol;
      const int artificialUsed = artificial ? numberArtificialUsed(sol) : 0;
      addProgress(artificial ? 1 : 2, value, artificialUsed);
      if (!artificial) {
        st.timeToBest = elapsed();
        realBestWork = work_ - startWork;
        Improvement imp = { st.lnsIterations, value, st.timeToBest, work_ - startWork };
        st.improvements.push_back(imp);
      }
    }
    if (artificial && value < realBestValue - 1.0e-7 && !numberArtificialUsed(sol)) {
      realBestValue = value;
      realBest = sol;
      st.timeToBest = elapsed();
      realBestWork = work_ - startWork;
      addProgress(2, value, 0);
      Improvement imp = { st.lnsIterations, value, st.timeToBest, work_ - startWork };
      st.improvements.push_back(imp);
    }
  };
  if (artificial && bestValue != COIN_DBL_MAX && !numberArtificialUsed(best)) {
    realBestValue = bestValue;
    realBest = best;
    realBestWork = work_ - startWork;
  } else if (!artificial && bestValue != COIN_DBL_MAX) {
    realBestWork = work_ - startWork;
  }
  double penaltyScale = 1.0;
  auto finish = [&]() {
    st.time = elapsed();
    st.work = work_ - startWork;
    if (artificial) {
      st.artificialUsed = bestValue == COIN_DBL_MAX ? -1 : numberArtificialUsed(best);
      if (penaltyScale != 1.0) {
        scaleArtificialPenalty(1.0);
        if (bestValue != COIN_DBL_MAX)
          bestValue = evaluate(best);
      }
      if (realBestValue != COIN_DBL_MAX) {
        best.swap(realBest);
        return realBestValue;
      }
    }
    return bestValue;
  };
  // GRASP: first construction deterministic
  for (int it = 0; it < params.graspIterations; it++) {
    if (it && elapsed() > params.timeLimit)
      break;
    Config c = config;
    if (!it) {
      c.colAlpha = 0.0;
      c.rowAlpha = 0.0;
    }
    double value = construct(c, trial);
    st.graspIterations++;
    if (!it)
      st.first = value;
    if (value < currentValue) {
      currentValue = value;
      current = trial;
    }
    record(value, trial);
  }
  if (initial) {
    std::vector< char > start = *initial;
    double value = evaluate(start);
    if (value < currentValue) {
      currentValue = value;
      current = start;
    }
    record(value, start);
  }
  if (currentValue == COIN_DBL_MAX) {
    current = best;
    currentValue = bestValue;
  }
  if (currentValue == COIN_DBL_MAX)
    return finish();
  // LNS: destroy around partitioning rows left to slacks, rebuild, improve
  std::vector< int > owner(numberRows_), seeds, equalRows, pool, freed, chosen;
  const int configuredRelatedThreshold = static_cast< int >(params.relatedRowProbability * 1000000.0);
  for (int i = 0; i < numberRows_; i++)
    if (isEqual_[i])
      equalRows.push_back(i);
  std::vector< char > partial;
  Config c = config;
  if (params.lnsColAlpha >= 0.0)
    c.colAlpha = params.lnsColAlpha;
  if (params.lnsRowAlpha >= 0.0)
    c.rowAlpha = params.lnsRowAlpha;
  if (params.lnsLocalSearch >= 0 && params.lnsLocalSearch <= 2)
    c.localSearch = params.lnsLocalSearch;
  if (params.lnsRowRule >= 0)
    c.rowRule = params.lnsRowRule;
  int destroyRows = params.destroyRows;
  int sinceBest = 0, noImprove = 0;
  const int randomRowThreshold = static_cast< int >(params.randomRowFraction * 1000000.0);
  CbcSppGraspDetail::LnsAcceptance acceptance(params.acceptMode,
    params.acceptRelative, params.lahcLength, params.saInitialTemperature,
    params.saCoolingWork, params.saReheatAfter, currentValue, work_ - startWork,
    params.saSchedule, params.saCalibrationIterations,
    params.saCalibrationTarget, params.saCalibrationQuantile);
  double roundStartWork = startWork;
  int round = 0;
  // artificialRoundWork > 0: escalate after that much work per round
  // without a real solution; maxWork then limits the total work
  const bool roundBudget = params.artificialRoundWork > 0.0;
  for (int it = 0; it < params.lnsIterations; it++) {
    const bool totalWorkStop = params.totalWorkLimit > 0.0
      && work_ - startWork >= params.totalWorkLimit;
    const bool stop = (params.maxNoImprove > 0 && noImprove >= params.maxNoImprove)
      || (params.maxWork > 0.0 && work_ - (roundBudget ? startWork : roundStartWork) >= params.maxWork)
      || totalWorkStop;
    if (totalWorkStop)
      break;
    if (stop && (roundBudget || !artificial || realBestValue != COIN_DBL_MAX || round >= params.artificialRounds))
      break;
    if (stop || (roundBudget && artificial && realBestValue == COIN_DBL_MAX && round < params.artificialRounds && work_ - roundStartWork >= params.artificialRoundWork)) {
      // best still needs artificials: make them dearer and go on from it
      round++;
      penaltyRound = round;
      penaltyScale *= params.artificialGrowth;
      scaleArtificialPenalty(penaltyScale);
      bestValue = evaluate(best);
      current = best;
      currentValue = bestValue;
      addProgress(3, currentValue, numberArtificialUsed(best));
      acceptance.reset(currentValue, work_ - startWork);
      roundStartWork = work_;
      noImprove = 0;
      sinceBest = 0;
      destroyRows = params.destroyRows;
    }
    if (elapsed() > params.timeLimit)
      break;
    const bool foundFeasible = artificial ? realBestValue != COIN_DBL_MAX : bestValue != COIN_DBL_MAX;
    if (params.lnsLocalSearch >= 3 && params.lnsLocalSearch <= 5) {
      const int before[] = { 1, 0, 2 };
      const int after[] = { 2, 2, 1 };
      c.localSearch = foundFeasible ? after[params.lnsLocalSearch - 3] : before[params.lnsLocalSearch - 3];
    }
    const int relatedThreshold = params.relatedAfterFeasible && !foundFeasible ? 0 : configuredRelatedThreshold;
    work_ += 2 * n + numberRows_;
    std::fill(owner.begin(), owner.end(), -1);
    seeds.clear();
    for (int j = 0; j < n; j++) {
      if (!current[j])
        continue;
      for (int k = colStart_[j]; k < colStart_[j + 1]; k++)
        owner[colRows_[k]] = j;
    }
    for (int i = 0; i < numberRows_; i++)
      if (isEqual_[i] && owner[i] >= 0 && owner[i] == slackCol_[i])
        seeds.push_back(i);
    if (seeds.empty())
      seeds = equalRows;
    partial = current;
    int size = destroyRows;
    if (params.destroyGrowProbability > 0.0) {
      const int maxDestroyRows = std::max(destroyRows, params.maxDestroyRows);
      const int threshold = static_cast< int >(params.destroyGrowProbability * 1000000.0);
      while (size < maxDestroyRows && randomInt(1000000) < threshold)
        size = std::min(2 * size, maxDestroyRows);
    }
    int numberSeeds = std::min(size, static_cast< int >(seeds.size()));
    freed.clear();
    chosen.clear();
    auto gain = [&](int i) {
      int o = owner[i];
      return o < 0 ? COIN_DBL_MAX : score_[o] - score_[sortedRowCols_[rowStart_[i]]];
    };
    for (int s = 0; s < numberSeeds; s++) {
      int r = -1;
      if (relatedThreshold && !freed.empty() && randomInt(1000000) < relatedThreshold) {
        // a row freed by the previous seeds: keep the neighbourhood connected
        for (int a = 0; a < 4 && r < 0; a++) {
          int t = freed[randomInt(static_cast< int >(freed.size()))];
          if (std::find(chosen.begin(), chosen.end(), t) == chosen.end())
            r = t;
        }
      }
      if (r >= 0) {
        // related row found
      } else if (randomRowThreshold && !equalRows.empty() && randomInt(1000000) < randomRowThreshold) {
        r = equalRows[randomInt(static_cast< int >(equalRows.size()))];
      } else {
        int left = static_cast< int >(seeds.size()) - s;
        int t = s + randomInt(left);
        if (params.gainTournament > 1) {
          double bestGain = gain(seeds[t]);
          work_ += params.gainTournament;
          for (int q = 1; q < params.gainTournament; q++) {
            int t2 = s + randomInt(left);
            double g = gain(seeds[t2]);
            if (g > bestGain) {
              bestGain = g;
              t = t2;
            }
          }
        }
        std::swap(seeds[s], seeds[t]);
        r = seeds[s];
      }
      chosen.push_back(r);
      int last = std::min(rowStart_[r] + params.seedColumns * std::max(1, params.seedPool), rowStart_[r + 1]);
      pool.assign(sortedRowCols_.begin() + rowStart_[r], sortedRowCols_.begin() + last);
      int numberPicked = std::min(params.seedColumns, static_cast< int >(pool.size()));
      for (int q = 0; q < numberPicked; q++) {
        if (static_cast< int >(pool.size()) > numberPicked) {
          int t = q + randomInt(static_cast< int >(pool.size()) - q);
          std::swap(pool[q], pool[t]);
        }
        int j = pool[q];
        work_ += colStart_[j + 1] - colStart_[j];
        for (int k = colStart_[j]; k < colStart_[j + 1]; k++) {
          int o = owner[colRows_[k]];
          if (o >= 0 && !forced_[o] && partial[o]) {
            partial[o] = 0;
            if (relatedThreshold) {
              work_ += colStart_[o + 1] - colStart_[o];
              for (int k2 = colStart_[o]; k2 < colStart_[o + 1]; k2++)
                if (isEqual_[colRows_[k2]])
                  freed.push_back(colRows_[k2]);
            }
          }
        }
      }
    }
    if (params.destroyFraction > 0.0) {
      work_ += n;
      int threshold = static_cast< int >(params.destroyFraction * 1000000.0);
      for (int j = 0; j < n; j++)
        if (partial[j] && !forced_[j] && randomInt(1000000) < threshold)
          partial[j] = 0;
    }
    double previousBest = bestValue;
    auto drawBacktrack = [&]() {
      return params.backtrackProbability > 0.0
        && (params.backtrackProbability >= 1.0 || randomInt(1000000) < params.backtrackProbability * 1000000.0);
    };
    bool backtrack = !params.backtrackArtificialOnly && drawBacktrack();
    double value = construct(c, trial, &partial, params.lnsFocus ? &current : NULL);
    if (params.backtrackArtificialOnly && params.backtrackProbability > 0.0) {
      if (value != COIN_DBL_MAX && !numberArtificialUsed(trial))
        st.repair.skippedFeasible++;
      else
        backtrack = drawBacktrack();
    }
    if (backtrack)
      value = backtrackingRepair(c, params, trial, partial, params.lnsFocus ? &current : NULL, value, st.repair);
    st.lnsIterations++;
    record(value, trial);
    if (bestValue < previousBest) {
      sinceBest = 0;
      noImprove = 0;
      destroyRows = params.destroyRows;
    } else {
      noImprove++;
      if (params.stallIterations > 0 && ++sinceBest >= params.stallIterations) {
        sinceBest = 0;
        int maxDestroyRows = std::max(params.destroyRows, params.maxDestroyRows);
        if (destroyRows >= maxDestroyRows) {
          // largest neighbourhood also stalled: restart from the best
          destroyRows = params.destroyRows;
          current = best;
          currentValue = bestValue;
        } else {
          destroyRows = std::min(2 * destroyRows, maxDestroyRows);
        }
      }
    }
    const double previousCurrentValue = currentValue;
    bool accepted = acceptance.accepts(value, currentValue, bestValue, it,
      work_ - startWork, [this](int n) { return randomInt(n); });
    st.acceptance.trials++;
    if (value == COIN_DBL_MAX || !std::isfinite(value))
      st.acceptance.invalid++;
    if (accepted) {
      st.acceptance.accepted++;
      if (value > currentValue + 1.0e-7)
        st.acceptance.acceptedWorse++;
    } else {
      st.acceptance.rejected++;
    }
    if (accepted) {
      currentValue = value;
      current.swap(trial);
    }
    if (acceptance.update(accepted, currentValue, it, work_ - startWork,
          value, previousCurrentValue))
      st.acceptance.reheats++;
    const double roundWork = work_ - (roundBudget ? startWork : roundStartWork);
    const bool workLimitReached = params.maxWork > 0.0 && roundWork >= params.maxWork;
    const bool noImproveLimitReached = params.maxNoImprove > 0 && noImprove >= params.maxNoImprove;
    if (params.maxRealStallWork > 0.0 && realBestWork >= 0.0
      && work_ - startWork - realBestWork >= params.maxRealStallWork
      && !workLimitReached && !noImproveLimitReached
      && elapsed() <= params.timeLimit && it + 1 < params.lnsIterations) {
      st.realStallStop = true;
      addProgress(4, artificial ? realBestValue : bestValue, 0);
      break;
    }
  }
  st.acceptance.calibrationTrials = acceptance.calibrationSeen();
  st.acceptance.calibrationPositive = acceptance.calibrationPositive();
  st.acceptance.calibratedTemperature = acceptance.initialTemperature();
  return finish();
}

// ---------------------------------------------------------------- heuristic

CbcHeuristicSppGrasp::CbcHeuristicSppGrasp()
  : CbcHeuristic()
  , pumpSourceFractionality_(0.0)
  , pumpSourceFractionalCount_(0)
  , rootLpFractionality_(0.0)
  , rootLpFractionalCount_(0)
  , needSlacks_(true)
  , artificialFactor_(0.0)
  , minimumNativeCoverage_(0.0)
  , numberCalls_(0)
{
  setHeuristicName("SppGrasp");
}

CbcHeuristicSppGrasp::CbcHeuristicSppGrasp(CbcModel &model)
  : CbcHeuristic(model)
  , pumpSourceFractionality_(0.0)
  , pumpSourceFractionalCount_(0)
  , rootLpFractionality_(0.0)
  , rootLpFractionalCount_(0)
  , needSlacks_(true)
  , artificialFactor_(0.0)
  , minimumNativeCoverage_(0.0)
  , numberCalls_(0)
{
  setHeuristicName("SppGrasp");
}

CbcHeuristicSppGrasp::CbcHeuristicSppGrasp(const CbcHeuristicSppGrasp &rhs)
  : CbcHeuristic(rhs)
  , config_(rhs.config_)
  , params_(rhs.params_)
  , pumpSourceFractionality_(0.0)
  , pumpSourceFractionalCount_(0)
  , rootLpFractionality_(0.0)
  , rootLpFractionalCount_(0)
  , needSlacks_(rhs.needSlacks_)
  , artificialFactor_(rhs.artificialFactor_)
  , minimumNativeCoverage_(rhs.minimumNativeCoverage_)
  , coverageTuningProfile_(rhs.coverageTuningProfile_)
  , coverageTuningStage_(rhs.coverageTuningStage_)
  , coverageTuningValues_(rhs.coverageTuningValues_)
  , numberCalls_(rhs.numberCalls_)
{
}

CbcHeuristicSppGrasp::~CbcHeuristicSppGrasp() {}

CbcHeuristicSppGrasp &CbcHeuristicSppGrasp::operator=(const CbcHeuristicSppGrasp &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    config_ = rhs.config_;
    params_ = rhs.params_;
    pumpFractionalGuide_.clear();
    rootLpSnapshot_.clear();
    rootReducedCostsSnapshot_.clear();
    pumpSourceFractionality_ = 0.0;
    pumpSourceFractionalCount_ = 0;
    rootLpFractionality_ = 0.0;
    rootLpFractionalCount_ = 0;
    needSlacks_ = rhs.needSlacks_;
    artificialFactor_ = rhs.artificialFactor_;
    minimumNativeCoverage_ = rhs.minimumNativeCoverage_;
    coverageTuningProfile_ = rhs.coverageTuningProfile_;
    coverageTuningStage_ = rhs.coverageTuningStage_;
    coverageTuningValues_ = rhs.coverageTuningValues_;
    numberCalls_ = rhs.numberCalls_;
  }
  return *this;
}

void CbcHeuristicSppGrasp::setPumpFractionalGuide(const double *values,
  const double *rootLp, const double *rootReducedCosts, int n,
  double sourceFractionality, int sourceFractionalCount,
  double rootFractionality, int rootFractionalCount)
{
  pumpFractionalGuide_.clear();
  rootLpSnapshot_.clear();
  rootReducedCostsSnapshot_.clear();
  if (!values || n <= 0)
    return;
  pumpFractionalGuide_.assign(values, values + n);
  if (rootLp && rootReducedCosts) {
    rootLpSnapshot_.assign(rootLp, rootLp + n);
    rootReducedCostsSnapshot_.assign(rootReducedCosts, rootReducedCosts + n);
  }
  pumpSourceFractionality_ = sourceFractionality;
  pumpSourceFractionalCount_ = sourceFractionalCount;
  rootLpFractionality_ = rootFractionality;
  rootLpFractionalCount_ = rootFractionalCount;
}

bool CbcHeuristicSppGrasp::setCoverageTuningProfile(const std::string &path,
  const std::vector< std::string > &explicitOverrides, std::string *reason,
  const std::string &stage)
{
  coverageTuningValues_.clear();
  coverageTuningProfile_.clear();
  coverageTuningStage_ = stage.empty() ? "active" : stage;
  if (path.empty())
    return true;
  std::ifstream input(path.c_str());
  if (!input) {
    if (reason)
      *reason = "cannot read SppGrasp coverage tuning profile: " + path;
    return false;
  }
  CbcSppGrasp::Config config;
  CbcSppGrasp::RunParams params;
  double artificial = 1.1;
  std::vector< std::pair< std::string, double > > values;
  std::string line;
  int lineNumber = 0;
  while (std::getline(input, line)) {
    ++lineNumber;
    line = line.substr(0, line.find('#'));
    std::istringstream fields(line);
    std::string key, extra;
    double value;
    if (!(fields >> key))
      continue;
    if (!(fields >> value) || fields >> extra
      || !applyCoverageTuningValue(key, value, config, params, artificial)) {
      if (reason) {
        std::ostringstream message;
        message << "invalid SppGrasp coverage tuning profile " << path
                << " at line " << lineNumber;
        *reason = message.str();
      }
      return false;
    }
    if (std::find(explicitOverrides.begin(), explicitOverrides.end(), key)
      == explicitOverrides.end())
      values.push_back(std::make_pair(key, value));
  }
  if (input.bad() || (input.fail() && !input.eof())) {
    if (reason)
      *reason = "cannot read SppGrasp coverage tuning profile: " + path;
    return false;
  }
  coverageTuningProfile_ = path;
  coverageTuningStage_ = stage;
  coverageTuningValues_.swap(values);
  return true;
}
CbcHeuristic *CbcHeuristicSppGrasp::clone() const
{
  return new CbcHeuristicSppGrasp(*this);
}

void CbcHeuristicSppGrasp::resetModel(CbcModel *model)
{
  setModel(model);
}

void CbcHeuristicSppGrasp::setModel(CbcModel *model)
{
  model_ = model;
}

void CbcHeuristicSppGrasp::generateCpp(FILE *fp)
{
  CbcHeuristicSppGrasp other;
  fprintf(fp, "0#include \"CbcHeuristicSppGrasp.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicSppGrasp heuristicSppGrasp(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "heuristicSppGrasp");
  fprintf(fp, "3  cbcModel->addHeuristic(&heuristicSppGrasp);\n");
}

int CbcHeuristicSppGrasp::applicability(const OsiSolverInterface &solver,
  double artificialFactor, std::string *reason, double minimumNativeCoverage)
{
  if (minimumNativeCoverage > 0.0) {
    CbcSppGrasp grasp;
    grasp.setArtificialSlacks(artificialFactor);
    if (!grasp.load(solver) || grasp.tooRestricted()) {
      if (reason)
        *reason = grasp.reason().empty() ? "mapped model is too restricted" : grasp.reason();
      return 0;
    }
    if (!grasp.meetsNativeCoverage(minimumNativeCoverage, reason))
      return 0;
    if (!grasp.trivialSolutionExists()) {
      if (reason)
        *reason += "; incomplete coverage and artificial slacks disabled: skip";
      return 0;
    }
    return grasp.numberArtificial() ? 2 : 1;
  }
  std::string why;
  if (isApplicable(solver, true, &why))
    return 1;
  if (reason)
    *reason = why;
  if (artificialFactor <= 0.0)
    return 0;
  CbcSppGrasp grasp;
  grasp.setArtificialSlacks(artificialFactor);
  if (!grasp.load(solver) || grasp.tooRestricted())
    return 0;
  return 2;
}

bool CbcHeuristicSppGrasp::isApplicable(const OsiSolverInterface &solver, bool needSlacks,
  std::string *reason)
{
  CbcSppGrasp grasp;
  if (!grasp.load(solver)) {
    if (reason)
      *reason = grasp.reason();
    return false;
  }
  if (grasp.tooRestricted()) {
    if (reason) {
      char line[200];
      sprintf(line, "%d of %d columns kept at zero by %d restricting rows",
        grasp.numberRestrictedColumns(), grasp.numberOriginalColumns(), grasp.numberRestrictingRows());
      *reason = line;
    }
    return false;
  }
  if (needSlacks && !grasp.trivialSolutionExists()) {
    if (reason) {
      char line[200];
      sprintf(line, "only %d of %d partitioning rows have a slack column",
        grasp.numberSlackRows(), grasp.numberPartitioningRows());
      *reason = line;
    }
    return false;
  }
  return true;
}

int CbcHeuristicSppGrasp::solution(double &objectiveValue, double *newSolution)
{
  return solution(objectiveValue, newSolution, false);
}

int CbcHeuristicSppGrasp::solution(double &objectiveValue, double *newSolution, bool standalone)
{
  std::vector< double > pumpGuide;
  pumpGuide.swap(pumpFractionalGuide_);
  std::vector< double > rootLpSnapshot;
  rootLpSnapshot.swap(rootLpSnapshot_);
  std::vector< double > rootReducedCostsSnapshot;
  rootReducedCostsSnapshot.swap(rootReducedCostsSnapshot_);
  const double pumpSourceFractionality = pumpSourceFractionality_;
  const int pumpSourceFractionalCount = pumpSourceFractionalCount_;
  const double rootLpFractionality = rootLpFractionality_;
  const int rootLpFractionalCount = rootLpFractionalCount_;
  if (!when() || !model_)
    return 0;
  OsiSolverInterface *solver = model_->solver();
  if (!solver)
    return 0;
  numberCalls_++;
  CbcSppGrasp grasp;
  grasp.setArtificialSlacks(artificialFactor_);
  int numberRows = model_->numberRowsAtContinuous();
  if (numberRows <= 0)
    numberRows = solver->getNumRows();
  if (!grasp.load(*solver, minimumNativeCoverage_ > 0.0 ? -1 : numberRows) || grasp.tooRestricted()) {
    if (model_->sppGraspStagedWorkEnabled()
      && model_->messageHandler()->logLevel() > 1) {
      model_->messageHandler()->message(CBC_GENERAL, model_->messages())
        << "SppGrasp staged pass skipped: "
        << (grasp.reason().empty() ? "mapped model is too restricted"
                                   : grasp.reason())
        << CoinMessageEol;
    }
    return 0;
  }
  if (minimumNativeCoverage_ > 0.0) {
    std::string coverageDescription;
    const bool accepted = grasp.meetsNativeCoverage(minimumNativeCoverage_, &coverageDescription);
    if (model_->messageHandler()->logLevel() > 1)
      model_->messageHandler()->message(CBC_GENERAL, model_->messages())
        << "SppGrasp active-call " + coverageDescription << CoinMessageEol;
    if (!accepted)
      return 0;
  }
  CbcSppGrasp::Config config = config_;
  CbcSppGrasp::RunParams params = params_;
  if (!coverageTuningProfile_.empty()) {
    const double profileStart = CoinGetTimeOfDay();
    // The qualification check uses all current constraints, including cuts.
    // If they are not an eligible SPP model, preserve the legacy call rather
    // than disabling it or applying a certificate from another stage.
    CbcSppGrasp active;
    CbcSppGrasp *qualifier = &grasp;
    if (numberRows < solver->getNumRows()) {
      active.setArtificialSlacks(artificialFactor_);
      if (active.load(*solver) && !active.tooRestricted())
        qualifier = &active;
      else
        qualifier = NULL;
    }
    std::string description;
    bool qualifies = qualifier && qualifier->meetsNativeCoverage(0.5, &description);
    if (!qualifier)
      description = active.reason().empty() ? "all-active model too restricted" : active.reason();
    if (qualifies) {
      double artificial = artificialFactor_;
      for (size_t i = 0; i < coverageTuningValues_.size(); ++i)
        applyCoverageTuningValue(coverageTuningValues_[i].first,
          coverageTuningValues_[i].second, config, params, artificial);
      if (artificial != artificialFactor_) {
        CbcSppGrasp tuned;
        tuned.setArtificialSlacks(artificial);
        if (!tuned.load(*solver) || tuned.tooRestricted()
          || (needSlacks_ && !tuned.trivialSolutionExists())) {
          qualifies = false;
          config = config_;
          params = params_;
          description += "; incompatible artificial setting, retaining legacy model";
        } else
          grasp = tuned;
      } else if (qualifier == &active)
        grasp = active;
    }
    if (params_.totalWorkLimit > 0.0)
      params.totalWorkLimit = params.totalWorkLimit > 0.0
        ? std::min(params.totalWorkLimit, params_.totalWorkLimit)
        : params_.totalWorkLimit;
    if (model_->messageHandler()->logLevel() > 1) {
      std::ostringstream message;
      message << "SppGrasp coverage tuning " << (qualifies ? "applied" : "legacy fallback")
              << ": stage " << coverageTuningStage_
              << "; " << (description.empty() ? "active model not eligible" : description)
              << "; profile " << coverageTuningProfile_
              << "; drows " << params.destroyRows << "; dgrow " << params.destroyGrowProbability
              << "; seedColumns " << params.seedColumns << "; acceptance " << params.acceptMode
              << "; lnsLocalSearch " << params.lnsLocalSearch
              << "; maxWork " << params.maxWork << "; noImprove " << params.maxNoImprove
              << "; totalWork " << params.totalWorkLimit
              << "; decisionSeconds " << CoinGetTimeOfDay() - profileStart;
      model_->messageHandler()->message(CBC_GENERAL, model_->messages())
        << message.str() << CoinMessageEol;
    }
  }
  const double *bestSolution = model_->bestSolution();
  const bool hasIncumbent = bestSolution != NULL;
  std::vector< char > incumbent;
  bool mappedIncumbentValid = false;
  if (bestSolution) {
    grasp.fromOriginal(bestSolution, incumbent);
    const double mappedObjective = grasp.evaluate(incumbent, false);
    mappedIncumbentValid = std::isfinite(mappedObjective)
      && mappedObjective != COIN_DBL_MAX
      && grasp.numberArtificialUsed(incumbent) == 0;
    if (!mappedIncumbentValid) {
      if (model_->sppGraspStagedWorkEnabled()
        && model_->messageHandler()->logLevel() > 1) {
        char line[256];
        sprintf(line,
          "SppGrasp staged incumbent map rejected (value %.12g, artificial "
          "rows %d, solver cols %d, SPP cols %d)",
          mappedObjective, grasp.numberArtificialUsed(incumbent),
          solver->getNumCols(), grasp.numberOriginalColumns());
        model_->messageHandler()->message(CBC_GENERAL, model_->messages())
          << line << CoinMessageEol;
      }
      incumbent.clear();
      if (model_->messageHandler()->logLevel() > 1) {
        model_->messageHandler()->message(CBC_GENERAL, model_->messages())
          << "SppGrasp: existing incumbent is not feasible in the mapped SPP model"
          << CoinMessageEol;
      }
    }
  }
  if (needSlacks_ && !grasp.trivialSolutionExists()
    && !(model_->sppGraspStagedWorkEnabled() && mappedIncumbentValid)) {
    if (model_->sppGraspStagedWorkEnabled()
      && model_->messageHandler()->logLevel() > 1) {
      model_->messageHandler()->message(CBC_GENERAL, model_->messages())
        << "SppGrasp staged pass skipped: no trivial slack solution"
        << CoinMessageEol;
    }
    return 0;
  }
  const bool noValidSppIncumbent = !hasIncumbent
    || (model_->sppGraspStagedWorkEnabled() && !mappedIncumbentValid);
  if (noValidSppIncumbent && solver->isProvenOptimal()
    && rootLpSnapshot.size() == static_cast< size_t >(grasp.numberOriginalColumns())
    && rootReducedCostsSnapshot.size() == rootLpSnapshot.size())
    grasp.setLpSolution(rootLpSnapshot.data(), rootReducedCostsSnapshot.data());
  else if (solver->isProvenOptimal())
    grasp.setLpSolution(solver->getColSolution(), solver->getReducedCost());
  bool usedPumpGuide = false;
  if (noValidSppIncumbent && solver->isProvenOptimal()
    && static_cast< int >(pumpGuide.size()) == grasp.numberOriginalColumns()) {
    grasp.setLpGuide(pumpGuide.data());
    usedPumpGuide = true;
    if (model_->messageHandler()->logLevel() > 1) {
      char line[256];
      sprintf(line, "SppGrasp: FPump guide selected (source fractionality %.6g, %d fractional; root LP %.6g, %d fractional)",
        pumpSourceFractionality, pumpSourceFractionalCount,
        rootLpFractionality, rootLpFractionalCount);
      model_->messageHandler()->message(CBC_GENERAL, model_->messages())
        << line << CoinMessageEol;
    }
  }
  grasp.setSeed(static_cast< unsigned int >(1234567 + 7919 * numberCalls_ + model_->getNodeCount()));
  const int n = grasp.numberOriginalColumns();
  if (model_->getMaximumSeconds() < 1.0e10)
    params.timeLimit = std::min(params.timeLimit,
      std::max(0.0, model_->getMaximumSeconds() - model_->getCurrentSeconds()));
  if (!coverageTuningProfile_.empty() && params.timeLimit <= 0.0) {
    if (model_->messageHandler()->logLevel() > 1)
      model_->messageHandler()->message(CBC_GENERAL, model_->messages())
        << "SppGrasp coverage tuning skipped: stage " + coverageTuningStage_ + "; model time limit exhausted"
        << CoinMessageEol;
    return 0;
  }
  // improve on the incumbent if there is one
  std::vector< char > best;
  CbcSppGrasp::RunStats stats;
  const bool inlineSummary = standalone && model_->messageHandler()->logLevel() == 1;
  FILE *fp = model_->messageHandler()->filePointer();
  const bool utf8 = CbcOutput::useUtf8();
  if (inlineSummary) {
    fprintf(fp, "\n%s SPPGrasp ... ", utf8 ? "▶" : "[>]");
    fflush(fp);
  }
  double bestValue = grasp.run(config, params, best, &stats,
    mappedIncumbentValid ? &incumbent : NULL);
  if (model_->messageHandler()->logLevel() > 1) {
    std::ostringstream message;
    message << (coverageTuningProfile_.empty() ? "SppGrasp run result: stage " : "SppGrasp coverage tuning result: stage ")
            << (coverageTuningStage_.empty() ? "active" : coverageTuningStage_)
            << "; work " << stats.work
            << "; iterations " << stats.lnsIterations << "; seconds " << stats.time
            << "; real " << (bestValue != COIN_DBL_MAX && !grasp.numberArtificialUsed(best) ? "yes" : "no");
    model_->messageHandler()->message(CBC_GENERAL, model_->messages())
      << message.str() << CoinMessageEol;
  }
  model_->addSppGraspStagedWork(stats.work);
  if (bestValue == COIN_DBL_MAX || grasp.numberArtificialUsed(best)) {
    if (model_->messageHandler()->logLevel() > 1) {
      char line[200];
      sprintf(line, "SppGrasp: no feasible solution (%d artificial rows) after %d LNS iterations in %.2f seconds",
        bestValue == COIN_DBL_MAX ? -1 : grasp.numberArtificialUsed(best), stats.lnsIterations, stats.time);
      model_->messageHandler()->message(CBC_GENERAL, model_->messages())
        << line << CoinMessageEol;
    }
    if (inlineSummary) {
      fprintf(fp, "%s No feasible solution in %.2f seconds.\n",
        utf8 ? "✘" : "[x]", stats.time);
      fflush(fp);
    }
    return 0;
  }
  const double direction = solver->getObjSenseInCbc();
  const double *objective = solver->getObjCoefficients();
  double offset = 0.0;
  solver->getDblParam(OsiObjOffset, offset);
  std::vector< double > x(n);
  grasp.toOriginal(best, x.data());
  double value = -offset;
  for (int j = 0; j < n; j++)
    value += x[j] * direction * objective[j];
  if (model_->messageHandler()->logLevel() > 1) {
    char line[200];
    sprintf(line, "SppGrasp: first %g, best %g (%d slack rows) after %d LNS iterations in %.2f seconds",
      stats.first - offset, value, grasp.numberSlacksUsed(best), stats.lnsIterations, stats.time);
    model_->messageHandler()->message(CBC_GENERAL, model_->messages())
      << line << CoinMessageEol;
  }
  if (inlineSummary) {
    fprintf(fp, "%s %.15g %s %.15g in %.2f seconds.\n",
      utf8 ? "✔" : "[+]",
      stats.first - offset, utf8 ? "→" : "->", value, stats.time);
    fflush(fp);
  }
  if (stats.realStallStop && model_->messageHandler()->logLevel() > 1) {
    char line[200];
    sprintf(line, "SppGrasp: progress gate stopped after %.0f work",
      stats.work);
    model_->messageHandler()->message(CBC_GENERAL, model_->messages())
      << line << CoinMessageEol;
  }
  if (usedPumpGuide && model_->messageHandler()->logLevel() > 1) {
    model_->messageHandler()->message(CBC_GENERAL, model_->messages())
      << "SppGrasp: FPump guide used as a soft LP-value bias (not an incumbent)"
      << CoinMessageEol;
  }
  if (value >= objectiveValue - 1.0e-7)
    return 0;
  std::copy(x.begin(), x.end(), newSolution);
  objectiveValue = value;
  return 1;
}
