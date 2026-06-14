// Copyright (C) 2026 MIPster contributors
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// This cut generator is based on the path aggregation algorithm in HiGHS
// (HighsPathSeparator.cpp). It aggregates binding LP rows along paths of
// continuous variables and calls MIR / cover cut generation on each
// aggregated row.

#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdio>
#include <vector>

#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPragma.hpp"
#include "OsiRowCut.hpp"
#include "OsiSolverInterface.hpp"

#include "CglPathAggregation.hpp"

namespace {

inline bool finiteBound(double value, double infinity)
{
  return value > -0.5 * infinity && value < 0.5 * infinity;
}

inline double cleanFloor(double value, double tolerance)
{
  double f = std::floor(value);
  if (value - f <= tolerance)
    return f;
  if (f + 1.0 - value <= tolerance)
    return f + 1.0;
  return f;
}

// Row types after checking binding at the LP solution
enum RowType { kUnusable = 0,
  kLeq,
  kGeq,
  kEq };

} // namespace

CglPathAggregation::CglPathAggregation()
  : CglCutGenerator()
  , maxPathLength_(6)
  , maxCuts_(200)
  , feasibilityTolerance_(1.0e-7)
  , minViolation_(1.0e-6)
{
}

CglPathAggregation::CglPathAggregation(const CglPathAggregation &rhs)
  : CglCutGenerator(rhs)
{
  gutsOfCopy(rhs);
}

CglPathAggregation &
CglPathAggregation::operator=(const CglPathAggregation &rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    gutsOfCopy(rhs);
  }
  return *this;
}

CglPathAggregation::~CglPathAggregation()
{
}

CglCutGenerator *
CglPathAggregation::clone() const
{
  return new CglPathAggregation(*this);
}

void CglPathAggregation::gutsOfCopy(const CglPathAggregation &rhs)
{
  maxPathLength_ = rhs.maxPathLength_;
  maxCuts_ = rhs.maxCuts_;
  feasibilityTolerance_ = rhs.feasibilityTolerance_;
  minViolation_ = rhs.minViolation_;
}

// ============================================================
// Main entry point
// ============================================================

void CglPathAggregation::generateCuts(const OsiSolverInterface &si,
  OsiCuts &cs, const CglTreeInfo info)
{
  (void)info;
  const int nRows = si.getNumRows();
  const int nCols = si.getNumCols();
  if (!nRows || !nCols)
    return;

  const double *xlp = si.getColSolution();
  if (!xlp)
    return;

  const double *rowAct = si.getRowActivity();
  const double *rowLb = si.getRowLower();
  const double *rowUb = si.getRowUpper();
  const double *colLb = si.getColLower();
  const double *colUb = si.getColUpper();
  const double infinity = si.getInfinity();
  const double feastol = feasibilityTolerance_;

  const CoinPackedMatrix &matByRow = *si.getMatrixByRow();
  const CoinPackedMatrix &matByCol = *si.getMatrixByCol();

  // ----------------------------------------------------------
  // 1. Classify rows: kLeq, kGeq, kEq, kUnusable
  //    A row is usable if it is binding (slack <= feastol)
  // ----------------------------------------------------------
  std::vector< RowType > rowType(nRows, kUnusable);
  for (int r = 0; r < nRows; ++r) {
    const bool hasLb = finiteBound(rowLb[r], infinity);
    const bool hasUb = finiteBound(rowUb[r], infinity);
    if (!hasLb && !hasUb)
      continue;
    if (hasLb && hasUb && fabs(rowUb[r] - rowLb[r]) <= feastol) {
      // Equality
      if (!rowAct || fabs(rowAct[r] - rowUb[r]) <= 10.0 * feastol)
        rowType[r] = kEq;
      continue;
    }
    if (hasUb && rowAct && rowUb[r] - rowAct[r] <= 10.0 * feastol) {
      rowType[r] = kLeq;
    } else if (hasLb && rowAct && rowAct[r] - rowLb[r] <= 10.0 * feastol) {
      rowType[r] = kGeq;
    } else if (!rowAct) {
      // No activity available — treat as potentially usable
      if (hasUb)
        rowType[r] = kLeq;
      else if (hasLb)
        rowType[r] = kGeq;
    }
  }

  // ----------------------------------------------------------
  // 2. Count continuous columns per row
  // ----------------------------------------------------------
  std::vector< int > numCont(nRows, 0);
  for (int c = 0; c < nCols; ++c) {
    if (si.isInteger(c))
      continue;
    // Only count if strictly between bounds (non-zero bound distance)
    const bool hasCL = finiteBound(colLb[c], infinity);
    const bool hasCU = finiteBound(colUb[c], infinity);
    double bd = infinity;
    if (hasCL)
      bd = std::min(bd, xlp[c] - colLb[c]);
    if (hasCU)
      bd = std::min(bd, colUb[c] - xlp[c]);
    if (bd <= feastol)
      continue; // at its bound — skip
    CoinShallowPackedVector cv = matByCol.getVector(c);
    for (int k = 0; k < cv.getNumElements(); ++k)
      ++numCont[cv.getIndices()[k]];
  }

  // ----------------------------------------------------------
  // 3. Identify free-substitution columns:
  //    equality rows with exactly one fractional continuous variable.
  //    Mark those rows kUnusable (consumed) and store substitution row/coeff.
  // ----------------------------------------------------------
  // colSubstRow[c] = row index to use for substituting col c (-1 if none)
  // colSubstCoeff[c] = the coefficient of col c in that row (with sign factor)
  std::vector< int > colSubstRow(nCols, -1);
  std::vector< double > colSubstCoeff(nCols, 0.0);

  for (int r = 0; r < nRows; ++r) {
    if (rowType[r] != kEq)
      continue;
    if (numCont[r] != 1)
      continue;

    CoinShallowPackedVector rv = matByRow.getVector(r);
    const int *idx = rv.getIndices();
    const double *val = rv.getElements();
    int col = -1;
    double coeff = 0.0;
    for (int k = 0; k < rv.getNumElements(); ++k) {
      if (si.isInteger(idx[k]))
        continue;
      // Check fractional (bound distance > 0)
      const bool hCL = finiteBound(colLb[idx[k]], infinity);
      const bool hCU = finiteBound(colUb[idx[k]], infinity);
      double bd = infinity;
      if (hCL)
        bd = std::min(bd, xlp[idx[k]] - colLb[idx[k]]);
      if (hCU)
        bd = std::min(bd, colUb[idx[k]] - xlp[idx[k]]);
      if (bd <= feastol)
        continue;
      col = idx[k];
      coeff = val[k];
      break;
    }
    if (col < 0)
      continue;
    if (colSubstRow[col] != -1)
      continue; // already claimed

    colSubstRow[col] = r;
    colSubstCoeff[col] = coeff;
    rowType[r] = kUnusable; // consumed — cannot be used as seed
  }

  // ----------------------------------------------------------
  // 4. Build inArc / outArc lists for each fractional continuous column.
  //
  //    Following HiGHS convention:
  //      LEQ row, positive coeff → outArc
  //      LEQ row, negative coeff → inArc
  //      GEQ row, positive coeff → inArc
  //      GEQ row, negative coeff → outArc
  //      EQ  row                 → both inArc and outArc
  //    (Only for columns that are fractional and NOT handled by free substitution)
  // ----------------------------------------------------------

  // inArcs[c] = list of (row, coeff) pairs where col c is an in-arc
  // outArcs[c] = list of (row, coeff) pairs where col c is an out-arc
  struct ArcEntry {
    int row;
    double coeff;
  };
  std::vector< std::vector< ArcEntry > > inArcs(nCols);
  std::vector< std::vector< ArcEntry > > outArcs(nCols);

  for (int c = 0; c < nCols; ++c) {
    if (si.isInteger(c))
      continue;
    if (colSubstRow[c] != -1)
      continue; // handled by free substitution
    // Check fractional
    const bool hasCL = finiteBound(colLb[c], infinity);
    const bool hasCU = finiteBound(colUb[c], infinity);
    double bd = infinity;
    if (hasCL)
      bd = std::min(bd, xlp[c] - colLb[c]);
    if (hasCU)
      bd = std::min(bd, colUb[c] - xlp[c]);
    if (bd <= feastol)
      continue;

    CoinShallowPackedVector cv = matByCol.getVector(c);
    const int *ridx = cv.getIndices();
    const double *rval = cv.getElements();
    for (int k = 0; k < cv.getNumElements(); ++k) {
      int r = ridx[k];
      double v = rval[k];
      switch (rowType[r]) {
      case kUnusable:
        break;
      case kLeq:
        if (v < 0.0)
          inArcs[c].push_back({ r, v });
        else
          outArcs[c].push_back({ r, v });
        break;
      case kGeq:
        if (v > 0.0)
          inArcs[c].push_back({ r, v });
        else
          outArcs[c].push_back({ r, v });
        break;
      case kEq:
        inArcs[c].push_back({ r, v });
        outArcs[c].push_back({ r, v });
        break;
      }
    }
  }

  // ----------------------------------------------------------
  // 5. Dense accumulator for aggregation
  // ----------------------------------------------------------
  std::vector< double > dense(nCols, 0.0);
  std::vector< int > activeList;
  activeList.reserve(256);
  std::vector< char > inActiveList(nCols, 0);
  std::vector< char > usedRows(nRows, 0);

  int generated = 0;

  // Helper: add a row to the dense accumulator with a given scale
  auto addRowScaled = [&](int r, double scale) {
    double sign = 1.0;
    if (rowType[r] == kGeq)
      sign = -1.0; // standardise to LEQ
    CoinShallowPackedVector rv = matByRow.getVector(r);
    const int *idx = rv.getIndices();
    const double *val = rv.getElements();
    double rhs_delta = (rowType[r] == kGeq) ? (-rowLb[r]) : rowUb[r];
    // Apply scale * sign
    double sc = scale * sign;
    for (int k = 0; k < rv.getNumElements(); ++k) {
      int c = idx[k];
      double newval = dense[c] + sc * val[k];
      if (fabs(newval) <= feastol)
        newval = 0.0;
      dense[c] = newval;
      if (!inActiveList[c]) {
        activeList.push_back(c);
        inActiveList[c] = 1;
      }
    }
    // return the rhs contribution
    return sc * rhs_delta;
  };

  // tryMirCut lambda removed (now generateMirCutFromRow is called directly)

  // Helper: select best continuous column to eliminate from current aggregation.
  // Returns column index (-1 = none found).
  // Returns sign: +1 means use inArc, -1 means use outArc.
  // Prefers largest bound distance.
  // Following HiGHS: positive coeff in current dense → look for inArc row
  //                  negative coeff in current dense → look for outArc row
  auto selectBestCol = [&](int &bestSign) -> int {
    int bestCol = -1;
    double bestBD = -1.0;
    int bSign = 0;
    for (int c : activeList) {
      double a = dense[c];
      if (fabs(a) <= feastol)
        continue;
      if (si.isInteger(c))
        continue;
      if (colSubstRow[c] != -1)
        continue; // handled as free substitution
      const bool hCL = finiteBound(colLb[c], infinity);
      const bool hCU = finiteBound(colUb[c], infinity);
      double bd = infinity;
      if (hCL)
        bd = std::min(bd, xlp[c] - colLb[c]);
      if (hCU)
        bd = std::min(bd, colUb[c] - xlp[c]);
      if (bd <= feastol)
        continue;

      // positive coeff → want inArc (to subtract it out)
      // negative coeff → want outArc
      const auto &arcs = (a > 0.0) ? inArcs[c] : outArcs[c];
      bool hasArc = false;
      for (const auto &ae : arcs) {
        if (!usedRows[ae.row]) {
          hasArc = true;
          break;
        }
      }
      if (!hasArc)
        continue;

      if (bd > bestBD) {
        bestBD = bd;
        bestCol = c;
        bSign = (a > 0.0) ? 1 : -1;
      }
    }
    bestSign = bSign;
    return bestCol;
  };

  // Helper: find best arc row for a given column+sign
  // Returns row index (-1 if none), sets weight
  auto findArcRow = [&](int col, int colSign, double &weight) -> int {
    double colCoeff = dense[col];
    const auto &arcs = (colSign > 0) ? inArcs[col] : outArcs[col];
    int bestRow = -1;
    int bestNewNZ = INT_MAX;
    double bestPivot = 0.0;

    for (const auto &ae : arcs) {
      int r = ae.row;
      if (usedRows[r])
        continue;
      double pivotCoeff = ae.coeff;
      // For GEQ rows, the sign is already flipped in the standardised form,
      // so negate pivot
      if (rowType[r] == kGeq)
        pivotCoeff = -pivotCoeff;
      if (fabs(pivotCoeff) <= feastol)
        continue;

      double w = -colCoeff / pivotCoeff;
      if (fabs(w) < feastol || fabs(w) > 1.0 / feastol)
        continue;

      // Count new non-zeros this row would add
      CoinShallowPackedVector rv = matByRow.getVector(r);
      int newNZ = 0;
      for (int k = 0; k < rv.getNumElements(); ++k) {
        int c = rv.getIndices()[k];
        if (!inActiveList[c]) {
          double contrib = w * (rowType[r] == kGeq ? -rv.getElements()[k] : rv.getElements()[k]);
          if (fabs(contrib) > feastol)
            ++newNZ;
        }
      }

      bool better = (bestRow < 0) || (newNZ < bestNewNZ);
      if (!better && newNZ == bestNewNZ)
        better = fabs(pivotCoeff) > bestPivot + feastol;
      if (!better && newNZ == bestNewNZ && fabs(fabs(pivotCoeff) - bestPivot) <= feastol)
        better = r < bestRow;

      if (better) {
        bestRow = r;
        bestNewNZ = newNZ;
        bestPivot = fabs(pivotCoeff);
        weight = w;
      }
    }
    return bestRow;
  };

  // ----------------------------------------------------------
  // 6. Main loop: iterate over seed rows
  // ----------------------------------------------------------
  for (int seedRow = 0; seedRow < nRows && generated < maxCuts_; ++seedRow) {
    if (rowType[seedRow] == kUnusable)
      continue;

    std::vector< double > seedScales;
    if (rowType[seedRow] == kLeq || rowType[seedRow] == kGeq) {
      seedScales.push_back(1.0);
    } else { // kEq
      seedScales.push_back(1.0);
      seedScales.push_back(-1.0);
    }

    for (double seedScale : seedScales) {
      // Reset accumulators
      for (int c : activeList) {
        dense[c] = 0.0;
        inActiveList[c] = 0;
      }
      activeList.clear();
      std::fill(usedRows.begin(), usedRows.end(), 0);

      // Unpack seed row
      double aggrRhs = addRowScaled(seedRow, seedScale);
      usedRows[seedRow] = 1;

      // Now extend the path step by step
      for (int depth = 0; depth < maxPathLength_ && generated < maxCuts_; ++depth) {

        // First: apply all free substitutions (equality rows with single continuous)
        bool appliedFreeSubst = false;
        for (int ci = 0; ci < (int)activeList.size(); ++ci) {
          int c = activeList[ci];
          if (fabs(dense[c]) <= feastol)
            continue;
          if (si.isInteger(c))
            continue;
          if (colSubstRow[c] < 0)
            continue;
          int sr = colSubstRow[c];
          if (usedRows[sr])
            continue;
          double pivot = colSubstCoeff[c];
          double w = -dense[c] / pivot;
          if (fabs(w) < feastol || fabs(w) > 1.0 / feastol)
            continue;
          aggrRhs += addRowScaled(sr, w);
          usedRows[sr] = 1;
          appliedFreeSubst = true;
        }

        // Try a cut at this aggregation level
        OsiRowCut cut;
        bool cutMade = generateMirCutFromRow(si, activeList, dense, aggrRhs, xlp, cut);
        if (cutMade) {
          cs.insertIfNotDuplicate(cut);
          ++generated;
          break;
        }

        // Check whether any fractional continuous variable remains that has arc rows
        int colSign = 0;
        int col = selectBestCol(colSign);
        if (col < 0)
          break;

        double weight = 0.0;
        int arcRow = findArcRow(col, colSign, weight);
        if (arcRow < 0)
          break;

        aggrRhs += addRowScaled(arcRow, weight);
        usedRows[arcRow] = 1;
      }
    }
  }
}

// ============================================================
std::string
CglPathAggregation::generateCpp(FILE *fp)
{
  CglPathAggregation other;
  fprintf(fp, "0#include \"CglPathAggregation.hpp\"\n");
  fprintf(fp, "3  CglPathAggregation pathAggregation;\n");
  if (maxPathLength_ != other.maxPathLength_)
    fprintf(fp, "3  pathAggregation.setMaxPathLength(%d);\n", maxPathLength_);
  if (maxCuts_ != other.maxCuts_)
    fprintf(fp, "3  pathAggregation.setMaxCuts(%d);\n", maxCuts_);
  if (feasibilityTolerance_ != other.feasibilityTolerance_)
    fprintf(fp, "3  pathAggregation.setFeasibilityTolerance(%.17g);\n", feasibilityTolerance_);
  if (minViolation_ != other.minViolation_)
    fprintf(fp, "3  pathAggregation.setMinViolation(%.17g);\n", minViolation_);
  return "pathAggregation";
}

// ============================================================
// generateMirCutFromRow — exposed for unit tests
// ============================================================
bool CglPathAggregation::generateMirCutFromRow(const OsiSolverInterface &si,
  const std::vector< int > &active, const std::vector< double > &dense,
  double rhs, const double *xlp, OsiRowCut &cut) const
{
  const int nCols = si.getNumCols();
  const double infinity = si.getInfinity();
  const double *colLb = si.getColLower();
  const double *colUb = si.getColUpper();
  const double feastol = feasibilityTolerance_;

  // 1. Gather all variables with non-zero coefficient in dense
  struct TransformedVar {
    int col;
    double a_orig;
    double a_trans;
    bool is_integer;
    bool complemented;
    double y_val;
    double bnd_val; // the bound value used
  };
  std::vector< TransformedVar > tVars;
  double yRhs = rhs;

  for (int c : active) {
    if (c < 0 || c >= nCols)
      return false;
    double a = dense[c];
    if (fabs(a) <= feastol)
      continue;

    const bool is_int = si.isInteger(c);
    const bool hCL = finiteBound(colLb[c], infinity);
    const bool hCU = finiteBound(colUb[c], infinity);

    if (!hCL && !hCU)
      return false;

    bool comp = false;
    if (hCU && (!hCL || colUb[c] - xlp[c] < xlp[c] - colLb[c]))
      comp = true;

    double a_trans = a;
    double y_val = 0.0;
    double bnd_val = 0.0;
    if (comp) {
      a_trans = -a;
      y_val = colUb[c] - xlp[c];
      bnd_val = colUb[c];
      yRhs -= a * colUb[c];
    } else {
      a_trans = a;
      y_val = xlp[c] - colLb[c];
      bnd_val = colLb[c];
      yRhs -= a * colLb[c];
    }

    tVars.push_back({ c, a, a_trans, is_int, comp, y_val, bnd_val });
  }

  if (tVars.empty())
    return false;

  // 2. Collect delta candidates for scaling
  std::vector< double > deltas;
  for (const auto &tv : tVars) {
    if (tv.is_integer) {
      double d = fabs(tv.a_trans);
      if (d > 1e-4) {
        deltas.push_back(d);
      }
    }
  }
  deltas.push_back(1.0);
  int num_candidates = deltas.size();
  for (int i = 0; i < num_candidates; ++i) {
    double d = deltas[i];
    if (d > 1e-4) {
      deltas.push_back(d * 2.0);
      deltas.push_back(d * 4.0);
      deltas.push_back(d * 8.0);
      deltas.push_back(d / 2.0);
      deltas.push_back(d / 4.0);
      deltas.push_back(d / 8.0);
    }
  }

  std::sort(deltas.begin(), deltas.end());
  std::vector< double > uniqueDeltas;
  for (double d : deltas) {
    if (d <= 1e-4)
      continue;
    if (uniqueDeltas.empty() || d - uniqueDeltas.back() > 10.0 * feastol) {
      uniqueDeltas.push_back(d);
    }
  }

  double bestEfficacy = minViolation_;
  double bestDelta = -1.0;
  double bestBetaFloor = 0.0;
  std::vector< double > bestAlphas;

  // 3. Try each delta candidate
  for (double delta : uniqueDeltas) {
    double scale = 1.0 / delta;
    double scalRhs = yRhs * scale;
    double betaFloor = cleanFloor(scalRhs, feastol);
    double f0 = scalRhs - betaFloor;
    if (f0 <= 0.005 || f0 >= 0.995)
      continue;

    double oneOverOneMinusF0 = 1.0 / (1.0 - f0);
    double violation = -betaFloor * delta;
    double normSq = 0.0;
    std::vector< double > alphas;
    alphas.reserve(tVars.size());

    for (const auto &tv : tVars) {
      double alpha = 0.0;
      if (tv.is_integer) {
        double scalaj = tv.a_trans * scale;
        double downaj = cleanFloor(scalaj, feastol);
        double fj = scalaj - downaj;
        double aj = downaj;
        if (fj > f0) {
          aj += (fj - f0) * oneOverOneMinusF0;
        }
        alpha = aj * delta;
      } else {
        if (tv.a_trans > 0.0) {
          alpha = 0.0;
        } else {
          alpha = tv.a_trans * oneOverOneMinusF0;
        }
      }
      alphas.push_back(alpha);
      violation += alpha * tv.y_val;
      normSq += alpha * alpha;
    }

    if (normSq <= 0.0)
      continue;

    double efficacy = violation / std::sqrt(normSq);
    if (efficacy > bestEfficacy) {
      bestEfficacy = efficacy;
      bestDelta = delta;
      bestBetaFloor = betaFloor;
      bestAlphas = alphas;
    }
  }

  if (bestDelta <= 0.0)
    return false;

  // 4. Construct cut
  std::vector< int > cutIdx;
  std::vector< double > cutCoef;
  double cutRhs = bestBetaFloor * bestDelta;

  for (size_t i = 0; i < tVars.size(); ++i) {
    double alpha = bestAlphas[i];
    if (fabs(alpha) <= feastol)
      continue;

    const auto &tv = tVars[i];
    if (tv.complemented) {
      cutIdx.push_back(tv.col);
      cutCoef.push_back(-alpha);
      cutRhs -= alpha * tv.bnd_val;
    } else {
      cutIdx.push_back(tv.col);
      cutCoef.push_back(alpha);
      cutRhs += alpha * tv.bnd_val;
    }
  }

  if (cutIdx.empty())
    return false;

  double activity = 0.0;
  double norm = 0.0;
  for (size_t i = 0; i < cutIdx.size(); ++i) {
    activity += cutCoef[i] * xlp[cutIdx[i]];
    norm += cutCoef[i] * cutCoef[i];
  }
  if (norm <= 0.0)
    return false;
  double finalViol = (activity - cutRhs) / std::sqrt(norm);
  if (finalViol <= minViolation_)
    return false;

  cut.setRow((int)cutIdx.size(), cutIdx.data(), cutCoef.data());
  cut.setLb(-infinity);
  cut.setUb(cutRhs);
  cut.setEffectiveness(finalViol);
  return true;
}
