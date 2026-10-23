// Copyright (C) 2026, COIN-OR Foundation and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include "CbcSppGraspAcceptance.hpp"

using CbcSppGraspDetail::LnsAcceptance;

static bool legacyAccept(int mode, double value, double current, double best,
  double threshold, double history)
{
  if (mode == LnsAcceptance::Threshold) {
    return value <= current + 1.0e-7
      || (value != std::numeric_limits< double >::max()
        && value <= best + threshold * fabs(best));
  }
  return value != std::numeric_limits< double >::max()
    && (value <= history + 1.0e-7 || value <= current + 1.0e-7);
}

static void checkLegacyPolicy(int mode)
{
  const double noValue = std::numeric_limits< double >::max();
  const double values[] = { 8.0, 10.5, 12.0, 9.0, 11.0, noValue, -2.0, 0.0 };
  const int length = mode == LnsAcceptance::Threshold ? 0 : 3;
  std::vector< double > history(length ? length : 1, 10.0);
  double current = 10.0;
  double best = 10.0;
  LnsAcceptance policy(mode, 0.1, length, 0.05, 100.0, 0, current, 0.0);
  for (int it = 0; it < static_cast< int >(sizeof(values) / sizeof(values[0])); it++) {
    const double value = values[it];
    double h = history[it % history.size()];
    bool expected = legacyAccept(mode, value, current, best, 0.1, h);
    bool actual = policy.accepts(value, current, best, it, it * 10.0,
      [](int) { assert(false); return 0; });
    assert(actual == expected);
    if (actual) {
      current = value;
      best = std::min(best, value);
    }
    if (mode != LnsAcceptance::Threshold) {
      if (mode == LnsAcceptance::Lahc || current < h)
        history[it % history.size()] = current;
    }
    policy.update(actual, current, it, it * 10.0);
  }
}

static void checkInvalidTrials()
{
  const double invalid[] = {
    std::numeric_limits< double >::max(),
    std::numeric_limits< double >::quiet_NaN(),
    std::numeric_limits< double >::infinity()
  };
  for (int mode = LnsAcceptance::Threshold; mode <= LnsAcceptance::SimulatedAnnealing; mode++) {
    LnsAcceptance policy(mode, 0.1, 4, 1.0, 100.0, 0, 4.0, 0.0);
    int draws = 0;
    for (double value : invalid) {
      assert(!policy.accepts(value, 4.0, 4.0, 0, 0.0,
        [&draws](int) { draws++; return 0; }));
    }
    assert(draws == 0);
  }
}

static void checkNegativeAndZeroObjectives()
{
  LnsAcceptance negative(LnsAcceptance::SimulatedAnnealing, 0.1, 1,
    0.2, 1000.0, 0, -100.0, 0.0);
  int draws = 0;
  assert(negative.accepts(-99.0, -100.0, -100.0, 0, 0.0,
    [&draws](int) { draws++; return 0; }));
  assert(draws == 1);

  LnsAcceptance zero(LnsAcceptance::SimulatedAnnealing, 0.1, 1,
    2.0, 1000.0, 0, 0.0, 0.0);
  assert(zero.accepts(1.0, 0.0, 0.0, 0, 0.0,
    [](int) { return 0; }));

  LnsAcceptance cold(LnsAcceptance::SimulatedAnnealing, 0.0, 1,
    0.0, 1000.0, 0, 0.0, 0.0);
  assert(!cold.accepts(1.0, 0.0, 0.0, 0, 0.0,
    [](int) { return 0; }));
}

static void checkLateAcceptanceReset()
{
  LnsAcceptance policy(LnsAcceptance::LahcBetter, 0.0, 2,
    0.0, 1.0, 0, 10.0, 0.0);
  assert(!policy.accepts(11.0, 10.0, 10.0, 0, 0.0,
    [](int) { return 0; }));
  policy.update(false, 10.0, 0, 0.0);
  policy.reset(20.0, 100.0);
  assert(policy.accepts(19.0, 30.0, 30.0, 1, 100.0,
    [](int) { return 0; }));
}

static void checkCoolingAndReheating()
{
  LnsAcceptance policy(LnsAcceptance::SimulatedAnnealing, 0.1, 1,
    2.0, 10.0, 1, 0.0, 0.0);
  const auto drawHalf = [](int n) { return n / 2; };
  assert(policy.accepts(20.0, 0.0, 0.0, 0, 0.0, drawHalf));
  assert(!policy.accepts(20.0, 0.0, 0.0, 1, 10.0, drawHalf));
  assert(policy.update(false, 0.0, 1, 10.0));
  assert(policy.accepts(20.0, 0.0, 0.0, 2, 10.0, drawHalf));
  policy.reset(-5.0, 20.0);
  assert(policy.accepts(-4.0, -5.0, -5.0, 3, 20.0, drawHalf));
}

static void checkLinearCooling()
{
  LnsAcceptance policy(LnsAcceptance::SimulatedAnnealing, 0.0, 1,
    2.0, 10.0, 0, 0.0, 0.0, 1);
  assert(policy.accepts(1.0, 0.0, 0.0, 0, 0.0,
    [](int) { return 0; }));
  assert(policy.accepts(1.0, 0.0, 0.0, 1, 10.0,
    [](int) { return 0; }));
  assert(!policy.accepts(1.0, 0.0, 0.0, 2, 20.0,
    [](int) { return 0; }));
}

static void checkEmpiricalCalibration()
{
  LnsAcceptance policy(LnsAcceptance::SimulatedAnnealing, 0.0, 1,
    99.0, 1000.0, 0, -100.0, 0.0, 0, 2, 0.5, 0.5);
  assert(!policy.accepts(-90.0, -100.0, -100.0, 0, 0.0,
    [](int) { return 0; }));
  policy.update(false, -100.0, 0, 1.0, -90.0, -100.0);
  assert(!policy.accepts(-70.0, -100.0, -100.0, 1, 1.0,
    [](int) { return 0; }));
  policy.update(false, -100.0, 1, 2.0, -70.0, -100.0);
  assert(policy.calibrationSeen() == 2);
  assert(policy.calibrationPositive() == 2);
  assert(fabs(policy.initialTemperature() - 0.2 / 0.6931471805599453) < 1.0e-12);
  assert(policy.accepts(-90.0, -100.0, -100.0, 2, 2.0,
    [](int) { return 0; }));

  LnsAcceptance noWorsening(LnsAcceptance::SimulatedAnnealing, 0.1, 1,
    99.0, 100.0, 0, 5.0, 0.0, 0, 3, 0.5, 0.5);
  noWorsening.update(true, 4.0, 0, 1.0, 4.0, 5.0);
  noWorsening.update(true, 3.0, 1, 2.0, 3.0, 4.0);
  noWorsening.update(true, 2.0, 2, 3.0, 2.0, 3.0);
  assert(noWorsening.calibrationSeen() == 3);
  assert(noWorsening.calibrationPositive() == 0);
  assert(noWorsening.initialTemperature() == 0.0);
  noWorsening.reset(2.0, 10.0);
  assert(noWorsening.calibrationSeen() == 0);
  assert(noWorsening.calibrationPositive() == 0);
}

static void checkDeterministicTrajectory()
{
  auto run = []() {
    LnsAcceptance policy(LnsAcceptance::SimulatedAnnealing, 0.0, 1,
      0.2, 40.0, 4, -10.0, 0.0);
    unsigned int seed = 12345;
    std::vector< int > decisions;
    double current = -10.0;
    for (int it = 0; it < 20; it++) {
      const double trial = it % 3 == 0 ? -9.0 : -10.5 + (it % 5);
      bool accepted = policy.accepts(trial, current, -12.0, it, it * 2.0,
        [&seed](int n) {
          seed ^= seed << 13;
          seed ^= seed >> 17;
          seed ^= seed << 5;
          return static_cast< int >(seed % static_cast< unsigned int >(n));
        });
      decisions.push_back(accepted);
      if (accepted)
        current = trial;
      policy.update(accepted, current, it, it * 2.0);
    }
    return decisions;
  };
  assert(run() == run());
}

int main()
{
  checkLegacyPolicy(LnsAcceptance::Threshold);
  checkLegacyPolicy(LnsAcceptance::Lahc);
  checkLegacyPolicy(LnsAcceptance::LahcBetter);
  checkInvalidTrials();
  checkNegativeAndZeroObjectives();
  checkLateAcceptanceReset();
  checkCoolingAndReheating();
  checkLinearCooling();
  checkEmpiricalCalibration();
  checkDeterministicTrajectory();
  puts("PASS: legacy acceptance equivalence, invalid trials, SA scaling/cooling/reheat/reset and seeded trajectory");
  return 0;
}
