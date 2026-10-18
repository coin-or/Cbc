// Copyright (C) 2026, COIN-OR Foundation and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <numeric>

#include "CbcHeuristicSppGrasp.hpp"
#include "CbcModel.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinTime.hpp"
#include "OsiSolverInterface.hpp"

namespace {
const double kNoOption = 1.0e20;
}

CbcSppGrasp::CbcSppGrasp()
  : applicable_(false)
  , numberColumns_(0)
  , numberRows_(0)
  , numberEqual_(0)
  , numberSlackRows_(0)
  , objOffset_(0.0)
  , scoreRule_(-1)
  , scoreBeta_(0.0)
  , scoreLpWeight_(0.0)
  , scoreConflictWeight_(0.0)
  , scoresValid_(false)
  , seed_(1)
  , markStamp_(0)
{
}

bool CbcSppGrasp::load(const OsiSolverInterface &solver, int numberRows)
{
  applicable_ = false;
  scoresValid_ = false;
  lpX_.clear();
  lpDj_.clear();
  const int n = solver.getNumCols();
  int m = solver.getNumRows();
  if (numberRows >= 0 && numberRows < m)
    m = numberRows;
  numberColumns_ = n;
  numberRows_ = m;
  numberEqual_ = 0;
  numberSlackRows_ = 0;
  if (!n || !m) {
    reason_ = "empty problem";
    return false;
  }
  const double *lower = solver.getColLower();
  const double *upper = solver.getColUpper();
  const double *rowLower = solver.getRowLower();
  const double *rowUpper = solver.getRowUpper();
  const double *objective = solver.getObjCoefficients();
  const double direction = solver.getObjSense();
  const double inf = solver.getInfinity();
  const double tol = 1.0e-9;

  isEqual_.assign(m, 0);
  std::vector< char > freeRow(m, 0);
  for (int i = 0; i < m; i++) {
    if (rowUpper[i] >= inf && rowLower[i] <= -inf) {
      freeRow[i] = 1;
      continue;
    }
    if (fabs(rowUpper[i]) <= tol && rowLower[i] <= tol) {
      // left by preprocessing: every column in it must be zero
      freeRow[i] = 2;
      continue;
    }
    if (fabs(rowUpper[i] - 1.0) > tol) {
      reason_ = "row with right-hand side other than 1";
      return false;
    }
    if (rowLower[i] > 1.0 + tol) {
      reason_ = "row with lower bound above 1";
      return false;
    }
    if (rowLower[i] > 1.0e-6) {
      isEqual_[i] = 1;
      numberEqual_++;
    }
  }
  if (!numberEqual_) {
    reason_ = "no partitioning rows";
    return false;
  }
  usable_.assign(n, 0);
  forced_.assign(n, 0);
  cost_.resize(n);
  for (int j = 0; j < n; j++) {
    cost_[j] = direction * objective[j];
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
  const CoinPackedMatrix *matrix = solver.getMatrixByCol();
  const double *element = matrix->getElements();
  const int *row = matrix->getIndices();
  const CoinBigIndex *start = matrix->getVectorStarts();
  const int *length = matrix->getVectorLengths();
  for (int j = 0; j < n; j++) {
    for (CoinBigIndex k = start[j]; k < start[j] + length[j]; k++) {
      int i = row[k];
      if (i >= m || freeRow[i] != 2 || fabs(element[k]) < 1.0e-12)
        continue;
      if (element[k] < 0.0 || forced_[j]) {
        reason_ = "row with right-hand side 0";
        return false;
      }
      usable_[j] = 0;
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
      if (fabs(element[k] - 1.0) > tol) {
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
  // columns of each partitioning row by cost per partitioning row
  repairOrder_ = rowCols_;
  minCostPerRow_.assign(m, 0.0);
  for (int i = 0; i < m; i++) {
    if (!isEqual_[i])
      continue;
    auto perRow = [&](int c) { return cost_[c] / numberEqualInCol_[c]; };
    std::sort(repairOrder_.begin() + rowStart_[i], repairOrder_.begin() + rowStart_[i + 1],
      [&](int a, int b) {
        double pa = perRow(a), pb = perRow(b);
        return pa != pb ? pa < pb : a < b;
      });
    minCostPerRow_[i] = perRow(repairOrder_[rowStart_[i]]);
  }
  objOffset_ = 0.0;
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

void CbcSppGrasp::setLpSolution(const double *x, const double *reducedCost)
{
  scoresValid_ = false;
  if (x)
    lpX_.assign(x, x + numberColumns_);
  else
    lpX_.clear();
  if (reducedCost)
    lpDj_.assign(reducedCost, reducedCost + numberColumns_);
  else
    lpDj_.clear();
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

double CbcSppGrasp::evaluate(const std::vector< char > &solution) const
{
  std::vector< int > cover(numberRows_, 0);
  double value = objOffset_;
  for (int j = 0; j < numberColumns_; j++) {
    if (!solution[j]) {
      if (forced_[j])
        return COIN_DBL_MAX;
      continue;
    }
    if (!usable_[j])
      return COIN_DBL_MAX;
    value += cost_[j];
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
  const std::vector< char > *partial)
{
  if (!applicable_)
    return COIN_DBL_MAX;
  computeScores(config);
  const int n = numberColumns_;
  const int m = numberRows_;
  solution.assign(n, 0);
  for (int j = 0; j < n; j++)
    alive_[j] = usable_[j];
  for (int i = 0; i < m; i++) {
    owner_[i] = -1;
    pointer_[i] = rowStart_[i];
    pointer_[m + i] = rowStart_[i];
    count_[i] = rowStart_[i + 1] - rowStart_[i];
  }
  auto kill = [&](int c) {
    alive_[c] = 0;
    for (int k = colStart_[c]; k < colStart_[c + 1]; k++)
      count_[colRows_[k]]--;
  };
  auto select = [&](int c) {
    assert(alive_[c]);
    solution[c] = 1;
    for (int k = colStart_[c]; k < colStart_[c + 1]; k++) {
      int i = colRows_[k];
      owner_[i] = c;
      for (int l = rowStart_[i]; l < rowStart_[i + 1]; l++) {
        int c2 = rowCols_[l];
        if (alive_[c2])
          kill(c2);
      }
    }
  };
  for (int j = 0; j < n; j++) {
    if (forced_[j]) {
      if (!alive_[j])
        return COIN_DBL_MAX;
      select(j);
    }
  }
  if (partial) {
    for (int j = 0; j < n; j++) {
      if ((*partial)[j] && !solution[j]) {
        if (!alive_[j])
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
    while (p < rowStart_[i + 1] && !alive_[sortedRowCols_[p]])
      p++;
    return p < rowStart_[i + 1] ? sortedRowCols_[p] : -1;
  };
  auto secondAlive = [&](int i) {
    int &p = pointer_[m + i];
    if (p <= pointer_[i])
      p = pointer_[i] + 1;
    while (p < rowStart_[i + 1] && !alive_[sortedRowCols_[p]])
      p++;
    return p < rowStart_[i + 1] ? sortedRowCols_[p] : -1;
  };
  while (!open.empty()) {
    // drop covered rows, settle rows whose best option is a singleton
    bool settled = true;
    while (settled) {
      settled = false;
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
        if (!alive_[c])
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
    if (alive_[j] && cost_[j] < 0.0)
      select(j);
  if (config.localSearch)
    return localSearch(solution, config.localSearch, partial);
  return evaluate(solution);
}

double CbcSppGrasp::localSearch(std::vector< char > &solution, int level,
  const std::vector< char > *reference)
{
  const int n = numberColumns_;
  const int m = numberRows_;
  std::fill(owner_.begin(), owner_.end(), -1);
  for (int j = 0; j < n; j++) {
    if (!solution[j])
      continue;
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
  if (reference) {
    // only rows not covered by a column of reference changed
    for (int i = 0; i < m; i++) {
      int o = owner_[i];
      if (o >= 0 && (*reference)[o])
        continue;
      for (int l = rowStart_[i]; l < rowStart_[i + 1]; l++)
        push(rowCols_[l]);
    }
  } else {
    for (int j = 0; j < n; j++)
      if (colStart_[j + 1] > colStart_[j])
        push(j);
  }
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
      if (!feasible || delta >= -1.0e-9)
        continue;
      // rows of ejected columns left open by j; bound repair cost from below
      changed.clear();
      for (int k = colStart_[j]; k < colStart_[j + 1]; k++)
        setOwner(colRows_[k], tentative);
      freed.clear();
      double bound = delta;
      for (size_t e = 0; e < ejected.size(); e++) {
        int o = ejected[e];
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
      if (feasible && bound < -1.0e-9) {
        if (level < 2) {
          delta = bound;
          for (size_t f = 0; f < freed.size(); f++) {
            added.push_back(slackCol_[freed[f]]);
            setOwner(freed[f], tentative);
          }
        } else {
          // greedy repair: cheapest (per row) column whose rows are all open
          for (size_t f = 0; f < freed.size() && feasible; f++) {
            int i = freed[f];
            if (owner_[i] != -1)
              continue;
            int chosen = -1, scanned = 0;
            for (int p = rowStart_[i]; p < rowStart_[i + 1] && scanned < 64; p++) {
              int c = repairOrder_[p];
              scanned++;
              bool free = true;
              for (int k = colStart_[c]; k < colStart_[c + 1]; k++) {
                if (owner_[colRows_[k]] != -1) {
                  free = false;
                  break;
                }
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
            for (int k = colStart_[chosen]; k < colStart_[chosen + 1]; k++)
              setOwner(colRows_[k], tentative);
          }
        }
      } else {
        feasible = false;
      }
      if (!feasible || delta >= -1.0e-9) {
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
        for (int l = rowStart_[i]; l < rowStart_[i + 1]; l++)
          push(rowCols_[l]);
      }
    }
  }
  return evaluate(solution);
}

double CbcSppGrasp::run(const Config &config, const RunParams &params,
  std::vector< char > &best, RunStats *stats, const std::vector< char > *initial)
{
  const double start = CoinGetTimeOfDay();
  const int n = numberColumns_;
  RunStats local;
  RunStats &st = stats ? *stats : local;
  st = RunStats();
  double bestValue = trivialSolution(best);
  st.trivial = bestValue;
  std::vector< char > current, trial;
  double currentValue = COIN_DBL_MAX;
  auto elapsed = [&]() { return CoinGetTimeOfDay() - start; };
  auto record = [&](double value, std::vector< char > &sol) {
    if (value < bestValue - 1.0e-7) {
      bestValue = value;
      best = sol;
      st.timeToBest = elapsed();
    }
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
  if (currentValue == COIN_DBL_MAX) {
    st.time = elapsed();
    return bestValue;
  }
  // LNS: destroy around partitioning rows left to slacks, rebuild, improve
  std::vector< int > owner(numberRows_), seeds, equalRows, pool;
  for (int i = 0; i < numberRows_; i++)
    if (isEqual_[i])
      equalRows.push_back(i);
  std::vector< char > partial;
  Config c = config;
  int destroyRows = params.destroyRows;
  int sinceBest = 0;
  const int randomRowThreshold = static_cast< int >(params.randomRowFraction * 1000000.0);
  for (int it = 0; it < params.lnsIterations; it++) {
    if (elapsed() > params.timeLimit)
      break;
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
    int numberSeeds = std::min(destroyRows, static_cast< int >(seeds.size()));
    for (int s = 0; s < numberSeeds; s++) {
      int r;
      if (randomRowThreshold && !equalRows.empty() && randomInt(1000000) < randomRowThreshold) {
        r = equalRows[randomInt(static_cast< int >(equalRows.size()))];
      } else {
        int t = s + randomInt(static_cast< int >(seeds.size()) - s);
        std::swap(seeds[s], seeds[t]);
        r = seeds[s];
      }
      int last = std::min(rowStart_[r] + params.seedColumns * std::max(1, params.seedPool), rowStart_[r + 1]);
      pool.assign(sortedRowCols_.begin() + rowStart_[r], sortedRowCols_.begin() + last);
      int numberPicked = std::min(params.seedColumns, static_cast< int >(pool.size()));
      for (int q = 0; q < numberPicked; q++) {
        if (static_cast< int >(pool.size()) > numberPicked) {
          int t = q + randomInt(static_cast< int >(pool.size()) - q);
          std::swap(pool[q], pool[t]);
        }
        int j = pool[q];
        for (int k = colStart_[j]; k < colStart_[j + 1]; k++) {
          int o = owner[colRows_[k]];
          if (o >= 0 && !forced_[o])
            partial[o] = 0;
        }
      }
    }
    if (params.destroyFraction > 0.0) {
      int threshold = static_cast< int >(params.destroyFraction * 1000000.0);
      for (int j = 0; j < n; j++)
        if (partial[j] && !forced_[j] && randomInt(1000000) < threshold)
          partial[j] = 0;
    }
    double previousBest = bestValue;
    double value = construct(c, trial, &partial);
    st.lnsIterations++;
    record(value, trial);
    if (bestValue < previousBest) {
      sinceBest = 0;
      destroyRows = params.destroyRows;
    } else if (params.stallIterations > 0 && ++sinceBest >= params.stallIterations) {
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
    if (value <= currentValue + 1.0e-7
      || (value != COIN_DBL_MAX && value <= bestValue + params.acceptRelative * fabs(bestValue))) {
      currentValue = value;
      current.swap(trial);
    }
  }
  st.time = elapsed();
  return bestValue;
}

// ---------------------------------------------------------------- heuristic

CbcHeuristicSppGrasp::CbcHeuristicSppGrasp()
  : CbcHeuristic()
  , needSlacks_(true)
  , numberCalls_(0)
{
  setHeuristicName("SppGrasp");
}

CbcHeuristicSppGrasp::CbcHeuristicSppGrasp(CbcModel &model)
  : CbcHeuristic(model)
  , needSlacks_(true)
  , numberCalls_(0)
{
  setHeuristicName("SppGrasp");
}

CbcHeuristicSppGrasp::CbcHeuristicSppGrasp(const CbcHeuristicSppGrasp &rhs)
  : CbcHeuristic(rhs)
  , config_(rhs.config_)
  , params_(rhs.params_)
  , needSlacks_(rhs.needSlacks_)
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
    needSlacks_ = rhs.needSlacks_;
    numberCalls_ = rhs.numberCalls_;
  }
  return *this;
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

bool CbcHeuristicSppGrasp::isApplicable(const OsiSolverInterface &solver, bool needSlacks,
  std::string *reason)
{
  CbcSppGrasp grasp;
  if (!grasp.load(solver)) {
    if (reason)
      *reason = grasp.reason();
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
  if (!when() || !model_)
    return 0;
  OsiSolverInterface *solver = model_->solver();
  if (!solver)
    return 0;
  numberCalls_++;
  CbcSppGrasp grasp;
  int numberRows = model_->numberRowsAtContinuous();
  if (numberRows <= 0)
    numberRows = solver->getNumRows();
  if (!grasp.load(*solver, numberRows))
    return 0;
  if (needSlacks_ && !grasp.trivialSolutionExists())
    return 0;
  if (solver->isProvenOptimal())
    grasp.setLpSolution(solver->getColSolution(), solver->getReducedCost());
  grasp.setSeed(static_cast< unsigned int >(1234567 + 7919 * numberCalls_ + model_->getNodeCount()));
  const int n = grasp.numberColumns();
  CbcSppGrasp::RunParams params = params_;
  if (model_->getMaximumSeconds() < 1.0e10)
    params.timeLimit = std::min(params.timeLimit,
      std::max(0.0, model_->getMaximumSeconds() - model_->getCurrentSeconds()));
  // improve on the incumbent if there is one
  std::vector< char > incumbent;
  const double *bestSolution = model_->bestSolution();
  if (bestSolution && model_->getSolutionCount()) {
    incumbent.resize(n);
    for (int j = 0; j < n; j++)
      incumbent[j] = bestSolution[j] > 0.5 ? 1 : 0;
  }
  std::vector< char > best;
  CbcSppGrasp::RunStats stats;
  double bestValue = grasp.run(config_, params, best, &stats,
    incumbent.empty() ? NULL : &incumbent);
  if (bestValue == COIN_DBL_MAX)
    return 0;
  const double direction = solver->getObjSenseInCbc();
  const double *objective = solver->getObjCoefficients();
  double offset = 0.0;
  solver->getDblParam(OsiObjOffset, offset);
  double value = -offset;
  for (int j = 0; j < n; j++)
    value += best[j] ? direction * objective[j] : 0.0;
  if (model_->messageHandler()->logLevel() > 0) {
    char line[200];
    sprintf(line, "SppGrasp: first %g, best %g (%d slack rows) after %d LNS iterations in %.2f seconds",
      stats.first - offset, value, grasp.numberSlacksUsed(best), stats.lnsIterations, stats.time);
    model_->messageHandler()->message(CBC_GENERAL, model_->messages())
      << line << CoinMessageEol;
  }
  if (value >= objectiveValue - 1.0e-7)
    return 0;
  for (int j = 0; j < n; j++)
    newSolution[j] = best[j] ? 1.0 : 0.0;
  objectiveValue = value;
  return 1;
}
