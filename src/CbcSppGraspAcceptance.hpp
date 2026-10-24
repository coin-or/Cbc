// Copyright (C) 2026, COIN-OR Foundation and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcSppGraspAcceptance_H
#define CbcSppGraspAcceptance_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace CbcSppGraspDetail {

class LnsAcceptance {
public:
  enum Mode {
    Threshold = 0,
    Lahc = 1,
    LahcBetter = 2,
    SimulatedAnnealing = 3
  };

  LnsAcceptance(int mode, double acceptRelative, int lahcLength,
    double saInitialTemperature, double saCoolingWork, int saReheatAfter,
    double initialValue, double work, int saSchedule = 0,
    int calibrationIterations = 0, double calibrationTarget = 0.5,
    double calibrationQuantile = 0.5)
    : mode_(mode)
    , acceptRelative_(acceptRelative)
    , saInitialTemperature_(std::isfinite(saInitialTemperature)
          ? std::max(0.0, saInitialTemperature)
          : 0.0)
    , saTemperature_(saInitialTemperature_)
    , saCoolingWork_(std::isfinite(saCoolingWork) && saCoolingWork > 0.0
          ? saCoolingWork
          : 1.0)
    , saReheatAfter_(std::max(0, saReheatAfter))
    , saSchedule_(saSchedule == 1 ? 1 : 0)
    , calibrationIterations_(std::max(0, calibrationIterations))
    , calibrationTarget_(std::isfinite(calibrationTarget)
          ? std::min(0.99, std::max(0.01, calibrationTarget))
          : 0.5)
    , calibrationQuantile_(std::isfinite(calibrationQuantile)
          ? std::min(1.0, std::max(0.0, calibrationQuantile))
          : 0.5)
    , calibrationSeen_(0)
    , calibrationPositive_(0)
    , saPhaseStartWork_(work)
    , saRejected_(0)
  {
    if (mode_ == Lahc || mode_ == LahcBetter)
      history_.assign(std::max(1, lahcLength), initialValue);
  }

  template < typename RandomInt >
  bool accepts(double value, double currentValue, double bestValue,
    int iteration, double work, RandomInt randomInt) const
  {
    if (value == std::numeric_limits< double >::max() || !std::isfinite(value))
      return false;
    if (mode_ == Threshold) {
      return value <= currentValue + 1.0e-7
        || value <= bestValue + acceptRelative_ * fabs(bestValue);
    }
    if (mode_ == Lahc || mode_ == LahcBetter) {
      const double historyValue = history_[iteration % history_.size()];
      return value <= historyValue + 1.0e-7 || value <= currentValue + 1.0e-7;
    }
    if (mode_ != SimulatedAnnealing)
      return false;
    if (calibrationIterations_ > calibrationSeen_)
      return thresholdAccepts(value, currentValue, bestValue);
    if (value <= currentValue + 1.0e-7)
      return true;

    const double normalizedDelta = normalizeDelta(value, currentValue);
    const double workSinceReheat = std::max(0.0, work - saPhaseStartWork_);
    double temperature = saTemperature_;
    if (saSchedule_ == 1)
      temperature *= std::max(0.0, 1.0 - 0.5 * workSinceReheat / saCoolingWork_);
    else
      temperature *= exp(-0.6931471805599453 * workSinceReheat / saCoolingWork_);
    if (!(temperature > 0.0) || !std::isfinite(temperature))
      return false;
    const double probability = exp(-normalizedDelta / temperature);
    return probability >= 1.0
      || (probability > 0.0 && randomInt(1000000) < probability * 1000000.0);
  }

  bool update(bool accepted, double currentValue, int iteration, double work,
    double trialValue = std::numeric_limits< double >::quiet_NaN(),
    double previousCurrentValue = std::numeric_limits< double >::quiet_NaN())
  {
    const bool calibrationTrial = mode_ == SimulatedAnnealing
      && calibrationSeen_ < calibrationIterations_;
    if (calibrationTrial) {
      if (std::isfinite(trialValue) && std::isfinite(previousCurrentValue)
        && trialValue != std::numeric_limits< double >::max()
        && previousCurrentValue != std::numeric_limits< double >::max()
        && trialValue > previousCurrentValue + 1.0e-7) {
        calibrationDeltas_.push_back(normalizeDelta(trialValue, previousCurrentValue));
        calibrationPositive_++;
      }
      calibrationSeen_++;
      if (calibrationSeen_ == calibrationIterations_)
        calibrateTemperature(work);
    }
    if (mode_ == Lahc || mode_ == LahcBetter) {
      double &historyValue = history_[iteration % history_.size()];
      if (mode_ == Lahc || currentValue < historyValue)
        historyValue = currentValue;
    } else if (mode_ == SimulatedAnnealing && !calibrationTrial) {
      if (accepted) {
        saRejected_ = 0;
      } else if (saReheatAfter_ > 0 && ++saRejected_ >= saReheatAfter_) {
        saPhaseStartWork_ = work;
        saRejected_ = 0;
        return true;
      }
    }
    return false;
  }

  void reset(double currentValue, double work)
  {
    std::fill(history_.begin(), history_.end(), currentValue);
    saPhaseStartWork_ = work;
    saRejected_ = 0;
    if (calibrationIterations_ > 0) {
      calibrationSeen_ = 0;
      calibrationPositive_ = 0;
      calibrationDeltas_.clear();
      saTemperature_ = saInitialTemperature_;
    }
  }

  int calibrationSeen() const { return calibrationSeen_; }
  int calibrationPositive() const { return calibrationPositive_; }
  double initialTemperature() const { return saTemperature_; }

private:
  static double normalizeDelta(double value, double currentValue)
  {
    const double scale = std::max(1.0, std::max(fabs(currentValue), fabs(value)));
    return value / scale - currentValue / scale;
  }

  bool thresholdAccepts(double value, double currentValue, double bestValue) const
  {
    return value <= currentValue + 1.0e-7
      || value <= bestValue + acceptRelative_ * fabs(bestValue);
  }

  void calibrateTemperature(double work)
  {
    saPhaseStartWork_ = work;
    if (calibrationDeltas_.empty()) {
      saTemperature_ = 0.0;
      return;
    }
    std::sort(calibrationDeltas_.begin(), calibrationDeltas_.end());
    const double rank = calibrationQuantile_ * (calibrationDeltas_.size() - 1);
    const size_t low = static_cast< size_t >(rank);
    const size_t high = std::min(low + 1, calibrationDeltas_.size() - 1);
    const double delta = calibrationDeltas_[low]
      + (rank - low) * (calibrationDeltas_[high] - calibrationDeltas_[low]);
    saTemperature_ = delta > 0.0
      ? -delta / std::log(calibrationTarget_)
      : 0.0;
  }

  int mode_;
  double acceptRelative_;
  double saInitialTemperature_;
  double saTemperature_;
  double saCoolingWork_;
  int saReheatAfter_;
  int saSchedule_;
  int calibrationIterations_;
  double calibrationTarget_;
  double calibrationQuantile_;
  int calibrationSeen_;
  int calibrationPositive_;
  std::vector< double > calibrationDeltas_;
  double saPhaseStartWork_;
  int saRejected_;
  std::vector< double > history_;
};

} // namespace CbcSppGraspDetail

#endif
