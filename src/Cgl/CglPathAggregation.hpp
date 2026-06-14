// Copyright (C) 2026 MIPster contributors
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglPathAggregation_H
#define CglPathAggregation_H

#include <cstdio>

#include "CglCutGenerator.hpp"

class CoinPackedMatrix;
class OsiRowCut;

void CglPathAggregationUnitTest();

class CGLLIB_EXPORT CglPathAggregation : public CglCutGenerator {
  friend void CglPathAggregationUnitTest();

public:
  CglPathAggregation();
  CglPathAggregation(const CglPathAggregation& rhs);
  CglPathAggregation& operator=(const CglPathAggregation& rhs);
  virtual ~CglPathAggregation();

  virtual CglCutGenerator* clone() const;

  virtual void generateCuts(const OsiSolverInterface& si, OsiCuts& cs,
    const CglTreeInfo info = CglTreeInfo());

  inline int getMaxPathLength() const { return maxPathLength_; }
  inline void setMaxPathLength(int value) { maxPathLength_ = value > 0 ? value : 1; }

  inline int getMaxCuts() const { return maxCuts_; }
  inline void setMaxCuts(int value) { maxCuts_ = value > 0 ? value : 1; }

  inline double getFeasibilityTolerance() const { return feasibilityTolerance_; }
  inline void setFeasibilityTolerance(double value) { feasibilityTolerance_ = value > 0.0 ? value : 1.0e-7; }

  inline double getMinViolation() const { return minViolation_; }
  inline void setMinViolation(double value) { minViolation_ = value > 0.0 ? value : 1.0e-6; }

  virtual std::string generateCpp(FILE* fp);

  // Exposed for unit testing
  bool generateMirCutFromRow(const OsiSolverInterface& si,
    const std::vector<int>& active, const std::vector<double>& dense,
    double rhs, const double* xlp, OsiRowCut& cut) const;

private:
  void gutsOfCopy(const CglPathAggregation& rhs);

  int maxPathLength_;
  int maxCuts_;
  double feasibilityTolerance_;
  double minViolation_;
};

#endif
