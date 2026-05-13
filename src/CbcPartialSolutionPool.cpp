// Copyright (C) 2026, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcPartialSolutionPool.hpp"
#include <algorithm>
#include <cassert>

CbcPartialSolutionPool::CbcPartialSolutionPool(int maxSize)
  : maxSize_(maxSize)
{
  pool_.reserve(maxSize);
}

bool CbcPartialSolutionPool::add(std::vector<CbcFixing> &&fixings,
  CoinWarmStart *basis, double score)
{
  if (static_cast<int>(pool_.size()) < maxSize_) {
    // Room available — always accept
    CbcPartialSolution ps;
    ps.fixings = std::move(fixings);
    ps.basis.reset(basis);
    ps.score = score;
    pool_.push_back(std::move(ps));
    return true;
  }

  // Pool full — find worst entry
  int worstIdx = 0;
  for (int i = 1; i < static_cast<int>(pool_.size()); i++) {
    if (pool_[i].score < pool_[worstIdx].score)
      worstIdx = i;
  }

  if (score <= pool_[worstIdx].score) {
    // Not better than worst — reject
    delete basis;
    return false;
  }

  // Replace worst
  CbcPartialSolution ps;
  ps.fixings = std::move(fixings);
  ps.basis.reset(basis);
  ps.score = score;
  pool_[worstIdx] = std::move(ps);
  return true;
}

const CbcPartialSolution &CbcPartialSolutionPool::select(int seed) const
{
  assert(!pool_.empty());
  if (pool_.size() == 1)
    return pool_[0];

  // Sort indices by score descending, pick from top half
  std::vector<int> indices(pool_.size());
  for (int i = 0; i < static_cast<int>(pool_.size()); i++)
    indices[i] = i;
  std::sort(indices.begin(), indices.end(),
    [&](int a, int b) { return pool_[a].score > pool_[b].score; });

  int topHalf = std::max(1, static_cast<int>(pool_.size()) / 2);
  int pick = ((unsigned)seed) % topHalf;
  return pool_[indices[pick]];
}

double CbcPartialSolutionPool::worstScore() const
{
  if (pool_.empty())
    return 0.0;
  double worst = pool_[0].score;
  for (int i = 1; i < static_cast<int>(pool_.size()); i++)
    worst = std::min(worst, pool_[i].score);
  return worst;
}
