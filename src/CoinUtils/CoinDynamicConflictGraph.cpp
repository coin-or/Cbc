/**
 *
 * This file is part of the COIN-OR CoinUtils project
 *
 * CoinConflictGraph implementation which supports modifications.
 * For a static conflict graph implemenation with faster queries
 * check CoinStaticConflictGraph.
 *
 * @file CoinDynamicConflictGraph.cpp
 * @brief CoinConflictGraph implementation which supports modifications.
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include <algorithm>
#include <climits>
#include <limits>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <algorithm>
#include "CoinDynamicConflictGraph.hpp"
#include "CoinStaticConflictGraph.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinCliqueList.hpp"
#include "CoinColumnType.hpp"
#include "CoinKnapsackRow.hpp"

#ifdef CGRAPH_STATS
#include <chrono>
#include <cstring>
#include "CoinRowType.hpp"
#endif

#define EPS 1e-6

#define CG_INI_SPACE_NODE_CONFLICTS 32

#define CG_INI_SPACE_ADJACENCY_VECTOR 32

// default sizes for clique lists (initial size, growth/reserve)
#define CG_LARGE_CLIQUE_INIT 4096
#define CG_LARGE_CLIQUE_GROW 32768
#define CG_SMALL_CLIQUE_INIT 4096
#define CG_SMALL_CLIQUE_GROW 3276

#ifdef CGRAPH_DEEP_DIVE
// enter here column/row to deep dive
#define CGRAPH_DEEP_DIVE_COLUMN_NAME "C1583"
#define CGRAPH_DEEP_DIVE_ROW_INDEX 6849
#endif // CGRAPH_DEEP_DIVE
// Tests involving nan seem to fail with some optimizations and cpu's
// This should be safe
inline bool isGoodNumber(double rhs)
{
  return (rhs<=std::numeric_limits<double>::max()&&
	  rhs>=std::numeric_limits<double>::min());
}
/**
 * @brief Locates the earliest column index that must belong to a clique.
 *
 * Runs a binary search over the sorted row copy to find the first position whose
 * sum with the subsequent element violates the tightened RHS.
 */
static size_t clique_start(const CoinTerm *columns, size_t nz, double rhs);

/**
 * @brief Finds the first column that conflicts with a fixed reference column.
 *
 * Searches within [colStart, colEnd] for the left-most index that, together with
 * columns[pos], exceeds the tightened RHS and therefore belongs to the same clique.
 */
static size_t binary_search(const CoinTerm *columns, size_t pos, double rhs, size_t colStart, size_t colEnd);

CoinDynamicConflictGraph::CoinDynamicConflictGraph(size_t _size)
  : CoinConflictGraph(_size)
  , conflicts(new CoinAdjacencyVector(_size, CG_INI_SPACE_NODE_CONFLICTS))
  , degree_(std::vector< size_t >(_size))
  , modifiedDegree_(std::vector< size_t >(_size))
  , largeClqs(new CoinCliqueList(CG_LARGE_CLIQUE_INIT, CG_LARGE_CLIQUE_GROW))
{
}

/* Build a dynamic conflict graph for a MIP model in several phases:
 *   1. Scan every constraint row once, tightening the RHS with fixed or
 *      continuous/general integer variables, complementing negative
 *      coefficients, and immediately storing explicit cliques or small
 *      cliques when they are detected.
 *   2. For rows that are not full cliques, keep their relevant columns and
 *      senses so they can be revisited later through clique detection on the
 *      processed (sorted) row copy.
 *   3. After all rows are inspected, add trivial variable/complement conflicts,
 *      analyze the stored partial rows to detect additional cliques, and then
 *      merge small cliques into the adjacency structure before finalizing
 *      node degrees.
 **/
CoinDynamicConflictGraph::CoinDynamicConflictGraph(
  const int numCols,
  const char *colType,
  const double *colLB,
  const double *colUB,
  const CoinPackedMatrix *matrixByRow,
  const char *sense,
  const double *rowRHS,
  const double *rowRange,
  const double primalTolerance,
  const double infinity,
  const std::vector< std::string > &colNames,
  const std::vector< std::string > &rowNames)
  : conflicts(new CoinAdjacencyVector(numCols * 2, CG_INI_SPACE_ADJACENCY_VECTOR))
  , largeClqs(new CoinCliqueList(CG_LARGE_CLIQUE_INIT, CG_LARGE_CLIQUE_GROW))
  , degree_(std::vector< size_t >(numCols * 2))
  , modifiedDegree_(std::vector< size_t >(numCols * 2))

{
  iniCoinConflictGraph(numCols * 2);

  // Buffers that keep "interesting" rows (and their tightened RHS) for the
  // delayed clique detection pass performed after the first matrix scan.
  this->tRowElements = std::vector< std::vector< CoinTerm > >();
  this->tRowRHS = std::vector< double >();

  // Mutable copies of column bounds — updated as fixings are discovered
  // so that subsequent rows benefit from tighter bounds.
  std::vector<double> mutableLB(colLB, colLB + numCols);
  std::vector<double> mutableUB(colUB, colUB + numCols);
  double *mColLB = mutableLB.data();
  double *mColUB = mutableUB.data();

  // Each equality/ranged row can spawn up to two temporary rows, so reserve 2x.
  this->tRowElements.reserve(matrixByRow->getNumRows() * 2);
  this->tRowRHS.reserve(matrixByRow->getNumRows() * 2);

  // constraint matrix input data
  const int *idxs = matrixByRow->getIndices();
  const double *coefs = matrixByRow->getElements();
  const CoinBigIndex *start = matrixByRow->getVectorStarts();
  const int *length = matrixByRow->getVectorLengths();

  // Reserve the staging area that accumulates the “small” cliques uncovered
  // while scanning each row; they will later be expanded into explicit
  // pairwise conflicts once all large cliques are indexed.
  smallCliques = new CoinCliqueList(CG_SMALL_CLIQUE_INIT, CG_SMALL_CLIQUE_GROW);

  // Workspace to copy indices of non-explicit cliques before sending them to smallCliques.
  size_t *tmpClq = new size_t[size_];

  CoinKnapsackRow knapsackRow(numCols,
    colType,
    mColLB,
    mColUB,
    primalTolerance,
    infinity);

  // arrays to hold multiplier/rhs adjustments
  // for equality/ranged rows
  double multiplierAdjustments[2];
  double rhsAdjustments[2];

  // inspecting all rows, compute initially tighthened rhs and
  // two largest and smallest values to check if constraint is
  // worth deeper inspection

  const size_t nRows = static_cast<size_t>(matrixByRow->getNumRows());

#ifdef CGRAPH_STATS
  memset(rowTypeStats_, 0, sizeof(rowTypeStats_));
#endif
  for (size_t idxRow = 0; idxRow < nRows; idxRow++) {
    const char rowSense = sense[idxRow];
    const CoinBigIndex rowStart = start[idxRow];
    const size_t rowLength = static_cast< size_t >(length[idxRow]);
    const int *indexesRow = idxs + rowStart;
    const double *coefficientsRow = coefs + rowStart;

#ifdef CGRAPH_STATS
    CoinRowType rtype = classifyRow(
      static_cast<int>(rowLength), indexesRow, coefficientsRow,
      rowSense, rowRHS[idxRow], colType, mColLB, mColUB);
    auto statsT0 = std::chrono::high_resolution_clock::now();
    bool rowFoundConflict = false;
    bool rowFoundFixing = false;
    size_t rowConflictsFound = 0;
    // Compute profile key: nz bucket, sense, abs(rhs) bucket
    int nzBkt = 0;
    { size_t rl = rowLength;
      if (rl>=256) nzBkt=7; else if (rl>=128) nzBkt=6; else if (rl>=64) nzBkt=5;
      else if (rl>=32) nzBkt=4; else if (rl>=16) nzBkt=3; else if (rl>=8) nzBkt=2;
      else if (rl>=4) nzBkt=1; }
    double absRhs = fabs(rowRHS[idxRow]);
    int rhsBkt = 0;
    if (absRhs > 100) rhsBkt=5; else if (absRhs > 20) rhsBkt=4;
    else if (absRhs > 5) rhsBkt=3; else if (absRhs > 1) rhsBkt=2;
    else if (absRhs > 0.5) rhsBkt=1;
#endif

#ifdef CGRAPH_DEEP_DIVE
    if (idxRow == CGRAPH_DEEP_DIVE_ROW_INDEX) {
      debugRowDetails(idxRow,
        sense,
        rowRHS,
        rowRange,
        idxs,
        coefs,
        start,
        length,
        colNames,
        colType,
        mColLB,
        mColUB);
    }
#endif // CGRAPH_DEEP_DIVE

    // equality / ranged constraints can generate two knapsack rows
    const int numRowIterations = CoinKnapsackRow::rowIterations(rowSense, rowRHS[idxRow],
      rowRange[idxRow], multiplierAdjustments, rhsAdjustments);

    for (int rowIter = 0; rowIter < numRowIterations; ++rowIter) {
      const double aMult = multiplierAdjustments[rowIter];
      const double aRHS = rhsAdjustments[rowIter];

      // basic processing
      knapsackRow.processRow(
        indexesRow, coefficientsRow, rowLength, rowSense, aMult, aRHS);

      processKnapsackRowFixings(knapsackRow, numCols, idxRow, colNames, rowNames);

      // Propagate fixings to mutable bounds so subsequent rows see them
      if (knapsackRow.nFixedVariables() > 0) {
        const int *fixedVars = knapsackRow.fixedVariables();
        for (size_t fi = 0; fi < knapsackRow.nFixedVariables(); ++fi) {
          size_t cidx = static_cast<size_t>(fixedVars[fi]);
          size_t origCol = cidx % numCols;
          if (cidx < (size_t)numCols) {
            mColLB[origCol] = 0.0; mColUB[origCol] = 0.0;
          } else {
            mColLB[origCol] = 1.0; mColUB[origCol] = 1.0;
          }
        }
      }

#ifdef CGRAPH_STATS
      if (knapsackRow.nFixedVariables() > 0)
        rowFoundFixing = true;
#endif

      const size_t nz = knapsackRow.nzs();

      // before going deeper
      if (knapsackRow.isUnbounded() || (nz == 0)) {
        continue;
      }

      double rhs = knapsackRow.rhs();
      const double *twoLargest = knapsackRow.twoLargest();
      const double *twoSmallest = knapsackRow.twoSmallest();

      if ((nz <= 1) || (twoLargest[0] + twoLargest[1] <= rhs + primalTolerance) || !isGoodNumber(rhs)) {
        continue; // no conflicts to search here
      }

      // explicit clique found - all variables conflict with each other
      if (knapsackRow.isExplicitClique()) {
        knapsackRow.copyColumnIndices(tmpClq);

        if (nz >= CoinConflictGraph::minClqRow_) {
          if (largeClqs->nCliques() < CoinConflictGraph::maxCliques_)
            processClique(tmpClq, nz);
        } else {
          smallCliques->addClique(nz, tmpClq);
        }
#ifdef CGRAPH_STATS
        rowFoundConflict = true;
        rowConflictsFound += nz * (nz - 1) / 2; // edges in clique
#endif
      } else {
        // partial clique - need to sort columns by coefficient
        if (twoLargest[0] != twoSmallest[0])
          knapsackRow.sortColumns();

        addTmpRow(nz, knapsackRow.columns(), rhs);
#ifdef CGRAPH_STATS
        rowFoundConflict = true;
        rowConflictsFound += nz; // approximate: at least nz vars involved
#endif
      } // not explicit clique
    } // row iterations (changed multiplier/rhs)

#ifdef CGRAPH_STATS
    {
      auto statsT1 = std::chrono::high_resolution_clock::now();
      double rowTime = std::chrono::duration<double>(statsT1 - statsT0).count();
      RowTypeStats &s = rowTypeStats_[rtype];
      s.nRows++;
      s.totalTime += rowTime;
      if (rowFoundConflict) s.rowsWithConflicts++;
      if (rowFoundFixing) s.rowsWithFixings++;
      // Row profile: find or create entry for this (nzBkt, sense, rhsBkt)
      RowProfileStats *rp = nullptr;
      for (auto &p : rowProfileStats_) {
        if (p.nzBucket == nzBkt && p.sense == rowSense && p.rhsBucket == rhsBkt) {
          rp = &p; break;
        }
      }
      if (!rp) {
        rowProfileStats_.push_back({nzBkt, rowSense, rhsBkt});
        rp = &rowProfileStats_.back();
      }
      rp->nRows++;
      rp->totalTime += rowTime;
      if (rowFoundConflict) rp->rowsWithConflicts++;
      if (rowFoundFixing) rp->rowsWithFixings++;
      rp->totalConflictsFound += rowConflictsFound;
    }
#endif

  } // all rows

  addVariableComplementConflicts(numCols, colType);

  // detecting cliques in less-structured constraints
  for (size_t idxTR = 0; (idxTR < tRowElements.size()); ++idxTR) {
    double rhs = tRowRHS[idxTR];
    if (isGoodNumber(rhs)) 
      cliqueDetection(tRowElements[idxTR], tRowElements[idxTR].size(), rhs + primalTolerance);
  }

  // at this point large cliques will already be include
  this->largeClqs->computeNodeOccurrences(size_);

  // processing small cliques
  if (smallCliques->nCliques()) {
    smallCliques->computeNodeOccurrences(size_);

    std::vector< char > iv(size_);

    for (size_t k = 0; (k < smallCliques->nDifferentNodes()); ++k) {
      size_t idxNode = smallCliques->differentNodes()[k];
      const size_t nNodeCliques = smallCliques->nNodeOccurrences(idxNode);
      const size_t *nodeCliques = smallCliques->nodeOccurrences(idxNode);
      processSmallCliquesNode(idxNode, nodeCliques, nNodeCliques, smallCliques, iv.data());
    }

  } // small cliques
  delete smallCliques;
  smallCliques = NULL;

  delete[] tmpClq;
  tmpClq = NULL;

  conflicts->flush();
  recomputeDegree();
}

void CoinDynamicConflictGraph::addClique(size_t size, const size_t elements[])
{
  this->largeClqs->addClique(size, elements);
}

void CoinDynamicConflictGraph::addVariableComplementConflicts(int numCols, const char *colType)
{
  for (size_t i = 0; i < static_cast<size_t>(numCols); ++i) {
    if (colType[i] == CoinColumnType::Binary) {
      conflicts->addNeighbor(i, i + numCols);
      conflicts->addNeighbor(numCols + i, i);
    }
  }
}

void CoinDynamicConflictGraph::addNodeConflicts(const size_t node, const size_t nodeConflicts[], const size_t nConflicts)
{
  for (size_t i = 0; i < nConflicts; ++i) {
    if (nodeConflicts[i] != node && (!conflictInCliques(node, nodeConflicts[i])))
      conflicts->addNeighbor(node, nodeConflicts[i], true); // also checks for repeated entries
  }
}

void CoinDynamicConflictGraph::addCliqueAsNormalConflicts(const size_t idxs[], const size_t len)
{
  for (size_t i = 0; (i < len); ++i)
    this->conflicts->addNeighborsBuffer(idxs[i], len, idxs);
}

void CoinDynamicConflictGraph::processClique(const size_t idxs[], const size_t size)
{
  if (size >= CoinConflictGraph::minClqRow_) {
    if (largeClqs->nCliques() < CoinConflictGraph::maxCliques_)
      addClique(size, idxs);
  } else {
    addCliqueAsNormalConflicts(idxs, size);
  }
}

/**
 * @brief Finds the first column index that forms an infeasible pair with a fixed column.
 *
 * The routine scans the sorted range $[\text{colStart}, \text{colEnd}]$ to locate the
 * left-most position whose coefficient, when added to @p columns[pos], exceeds the
 * constraint right-hand side. This allows clique detection to grow a clique starting
 * from variable @p pos.
 *
 * @param columns   Sorted (nondecreasing) array of columns belonging to the row.
 * @param pos       Index of the reference column that is being paired.
 * @param rhs       Tightened right-hand side of the knapsack row.
 * @param colStart  First index (inclusive) to be considered during the search.
 * @param colEnd    Last index (inclusive) to be considered during the search.
 *
 * @return The index of the first column in the search interval that violates the RHS
 *         when combined with @p columns[pos]. If all pairs remain feasible, returns
 *         @p colEnd + 1.
 */
size_t binary_search(const CoinTerm *columns, size_t pos, double rhs, size_t colStart, size_t colEnd)
{
  size_t left = colStart, right = colEnd;
  const double prevLHS = columns[pos].value;

#ifdef DEBUGCG
  assert(pos >= 0);
#endif

  while (left <= right) {
    const size_t mid = (left + right) / 2;
    const double lhs = prevLHS + columns[mid].value;

#ifdef DEBUGCG
    assert(mid >= 0);
    assert(mid <= colEnd);
    assert(pos <= colEnd);
#endif

    if (lhs <= rhs) {
      left = mid + 1;
    } else {
      if (mid > 0) {
        right = mid - 1;
      } else {
        return 0;
      }
    }
  }

  return right + 1;
}

/**
 * @brief Locates the earliest position where a clique must begin in a sorted row copy.
 *
 * Given a sorted row (nondecreasing coefficients), the function identifies the
 * smallest index @f$i@f$ such that the sum of @p columns[i] and the immediately
 * following term exceeds the row right-hand side. All indices from that position
 * onward must belong to at least one clique implied by the constraint.
 *
 * @param columns Sorted (nondecreasing) array of row elements.
 * @param nz      Number of non-zero elements in @p columns.
 * @param rhs     Tightened right-hand side of the knapsack row.
 *
 * @return The index where the first clique starts. A return value of zero indicates
 *         that the clique begins at the first element.
 */
size_t clique_start(const CoinTerm *columns, size_t nz, double rhs)
{
#ifdef DEBUGCG
  assert(nz > 1);
#endif

  size_t left = 0, right = nz - 2;

  while (left <= right) {
    const size_t mid = (left + right) / 2;
#ifdef DEBUGCG
    assert(mid >= 0);
    assert(mid + 1 < nz);
#endif
    const double lhs = columns[mid].value + columns[mid + 1].value;

    if (lhs <= rhs) {
      left = mid + 1;
    } else {
      if (mid > 0) {
        right = mid - 1;
      } else {
        return 0;
      }
    }
  }

  return right + 1;
}

void CoinDynamicConflictGraph::cliqueDetection(const std::vector< CoinTerm > &columns, size_t nz, const double rhs)
{
#ifdef DEBUGCG
  assert(nz > 1);
#endif

  if (columns[nz - 1].value + columns[nz - 2].value <= rhs) {
    return; // there is no clique in this constraint
  }

  size_t cliqueSize = 0;
  size_t *idxs = new size_t[nz];
  const size_t cliqueStart = clique_start(columns.data(), nz, rhs);

#ifdef DEBUGCG
  assert(cliqueStart >= 0);
  assert(cliqueStart <= nz - 2);
#endif

  for (size_t j = cliqueStart; j < nz; j++) {
    idxs[cliqueSize++] = static_cast<size_t>(columns[j].index);
  }

  // process the first clique found
  processClique(idxs, cliqueSize);

  // now we have to check the variables that are outside of the clique found.
  for (size_t j = cliqueStart; j-- > 0;) {
    const size_t idx = static_cast<size_t>(columns[j].index);

    if (columns[j].value + columns[nz - 1].value <= rhs) {
      break;
    }

    size_t position = binary_search(columns.data(), j, rhs, cliqueStart, nz - 1);

#ifdef DEBUGCG
    assert(position >= cliqueStart && position <= nz - 1);
#endif

    // new clique detected
    cliqueSize = 0;
    idxs[cliqueSize++] = idx;

    for (size_t k = position; k < nz; k++) {
      idxs[cliqueSize++] = static_cast<size_t>(columns[k].index);
    }

    processClique(idxs, cliqueSize);
  }

  delete[] idxs;
}

const std::vector< std::pair< size_t, std::pair< double, double > > > &CoinDynamicConflictGraph::updatedBounds()
{
  return newBounds_;
}

void CoinDynamicConflictGraph::addFixingBound(size_t columnIdx,
  size_t numCols,
  size_t idxRow,
  const std::vector< std::string > &colNames,
  const std::vector< std::string > &rowNames)
{
  size_t originalIdx = columnIdx;


  {
    // storing which row caused the implication, avoiding duplicates
    // and reporting if conflicting implications were found
    int varIdx = static_cast<int>(columnIdx % numCols);
    auto &implicationsForVar = columnRowImplications_[varIdx];
    if (implicationsForVar.find(static_cast<int>(columnIdx)) != implicationsForVar.end()) {
      // same implication already discovered
      return;
    }

    std::string rowName = (idxRow < rowNames.size()) ? rowNames[idxRow] : std::to_string(idxRow);
    implicationsForVar[static_cast<int>(columnIdx)] = std::make_pair(rowName, static_cast<int>(idxRow));

    int otherIdx = (columnIdx < numCols) ? varIdx + numCols : varIdx;

    if (implicationsForVar.find(otherIdx) != implicationsForVar.end()) {
      const auto &otherDetail = implicationsForVar[otherIdx];
      const bool currentFixesToZero = columnIdx < numCols;

      BinaryBoundInfeasibility infeasibility;
      infeasibility.variableIndex = static_cast<size_t>(varIdx);
      infeasibility.variableName = (varIdx < colNames.size()) ? colNames[varIdx] : std::to_string(varIdx);

      if (currentFixesToZero) {
        infeasibility.fixedToZero = {rowName, static_cast<int>(idxRow)};
        infeasibility.fixedToOne = {otherDetail.first, otherDetail.second};
      } else {
        infeasibility.fixedToOne = {rowName, static_cast<int>(idxRow)};
        infeasibility.fixedToZero = {otherDetail.first, otherDetail.second};
      }

      registerBoundImplicationInfeasibility(infeasibility);
    }
  }

  double lb = 0.0;
  double ub = 0.0;
  if (columnIdx >= (size_t)numCols) {
    originalIdx = columnIdx - numCols;
    lb = 1.0;
    ub = 1.0;
  }
#ifdef CGRAPH_DEEP_DIVE
  if (colNames[originalIdx] == CGRAPH_DEEP_DIVE_COLUMN_NAME) {
    printf("CGraph: column %s was fixed to %g based on analysis of row %zu\n", colNames[originalIdx].c_str(), lb, idxRow);
  }
#endif // CGRAPH_DEEP_DIVE
  newBounds_.push_back(std::make_pair(originalIdx, std::make_pair(lb, ub)));
}

void CoinDynamicConflictGraph::processKnapsackRowFixings(
  const CoinKnapsackRow &knapsackRowHelper,
  size_t numCols,
  size_t idxRow,
  const std::vector< std::string > &colNames,
  const std::vector< std::string > &rowNames)
{
  const size_t nFixings = knapsackRowHelper.nFixedVariables();
  const int *fixedVars = knapsackRowHelper.fixedVariables();

  for (size_t fixIdx = 0; fixIdx < nFixings; ++fixIdx) {
    addFixingBound(static_cast< size_t >(fixedVars[fixIdx]), numCols, idxRow, colNames, rowNames);
  }
}

const size_t *CoinDynamicConflictGraph::cliqueElements(size_t idxClique) const
{
  return largeClqs->cliqueElements(idxClique);
}

size_t CoinDynamicConflictGraph::getClqSize(size_t idxClique) const
{
  return largeClqs->cliqueSize(idxClique);
}

CoinDynamicConflictGraph::~CoinDynamicConflictGraph()
{
  delete conflicts;
  delete largeClqs;
  if (smallCliques)
    delete smallCliques;
}

#ifdef CGRAPH_DEEP_DIVE
void CoinDynamicConflictGraph::debugRowDetails(
  size_t idxRow,
  const char *sense,
  const double *rowRHS,
  const double *rowRange,
  const int *idxs,
  const double *coefs,
  const CoinBigIndex *start,
  const int *length,
  const std::vector< std::string > &colNames,
  const char *colType,
  const double *colLB,
  const double *colUB) const
{
  printf("Debugging row %zu\n", idxRow);
  const CoinBigIndex rowStart = start[idxRow];
  const CoinBigIndex rowEnd = rowStart + length[idxRow];
  std::ostringstream rowStream;
  rowStream << "Row " << idxRow << ": ";
  bool firstTerm = true;
  std::vector< size_t > rowVariableIdxs;
  rowVariableIdxs.reserve(length[idxRow]);

  for (CoinBigIndex j = rowStart; j < rowEnd; ++j) {
    const size_t idxCol = static_cast< size_t >(idxs[j]);
    const double coef = coefs[j];
    const double absCoef = fabs(coef);

    if (!firstTerm) {
      rowStream << (coef >= 0.0 ? "+ " : "- ");
    } else if (coef < 0.0) {
      rowStream << "- ";
    }

    if (fabs(absCoef - 1.0) > EPS) {
      rowStream << absCoef << " ";
    }

    rowStream << colNames[idxCol];
    firstTerm = false;

    if (std::find(rowVariableIdxs.begin(), rowVariableIdxs.end(), idxCol) == rowVariableIdxs.end()) {
      rowVariableIdxs.push_back(idxCol);
    }
  }

  rowStream << " ";
  const char rowSense = sense[idxRow];
  switch (rowSense) {
  case 'L':
    rowStream << "<= " << rowRHS[idxRow];
    break;
  case 'G':
    rowStream << ">= " << rowRHS[idxRow];
    break;
  case 'E':
    rowStream << "= " << rowRHS[idxRow];
    break;
  case 'R': {
    const double upper = rowRHS[idxRow];
    const double range = (rowRange) ? rowRange[idxRow] : 0.0;
    rowStream << "in [" << (upper - range) << ", " << upper << "]";
    break;
  }
  default:
    rowStream << "(sense " << rowSense << ") " << rowRHS[idxRow];
    break;
  }

  printf("%s\n", rowStream.str().c_str());
  printf("Variables in row %zu:\n", idxRow);
  for (size_t idxCol : rowVariableIdxs) {
    const char *typeName = CoinColumnType::nameFromChar(colType[idxCol]);
    printf("  - %s (type=%s, lb=%g, ub=%g)\n",
      colNames[idxCol].c_str(),
      typeName,
      colLB[idxCol],
      colUB[idxCol]);
  }
}
#endif // CGRAPH_DEEP_DIVE

void CoinDynamicConflictGraph::printInfo() const
{

  size_t minClq = INT_MAX;
  size_t maxClq = 0;

  for (size_t i = 0; (i < this->nCliques()); ++i) {
    minClq = std::min(minClq, cliqueSize(i));
    maxClq = std::max(maxClq, cliqueSize(i));
  }

  double totalDegree = 0.0;
  size_t minD = INT_MAX;
  size_t maxD = 0;
  for (size_t i = 0; (i < size()); ++i) {
    totalDegree += conflicts->rowSize(i);
    minD = std::min(minD, conflicts->rowSize(i));
    maxD = std::max(maxD, conflicts->rowSize(i));
  }
  double avd = totalDegree / ((double)size_);

  printf("Conflict graph info:\n");
  printf("\tnodes: %zu\n", this->size());
  printf("\tdensity: %.4f\n", this->density());
  printf("\tdegree min/max/av: %zu, %zu, %g\n", minD, maxD, avd);
  printf("\tncliques: %zu (min: %zu, max: %zu)\n", nCliques(), minClq, maxClq);
}

void CoinDynamicConflictGraph::addTmpRow(size_t nz, const CoinTerm els[], double rhs)
{
  tRowRHS.push_back(rhs);
  tRowElements.emplace_back(std::vector< CoinTerm >(els, els + nz));
}

size_t CoinDynamicConflictGraph::nDirectConflicts(size_t idxNode) const
{
  return this->conflicts->rowSize(idxNode);
}

const size_t *CoinDynamicConflictGraph::directConflicts(size_t idxNode) const
{
  return this->conflicts->getRow(idxNode);
}

std::vector<size_t> CoinDynamicConflictGraph::moveDirectConflicts(size_t idxNode)
{
  return this->conflicts->moveRow(idxNode);
}

std::vector<size_t> CoinDynamicConflictGraph::moveClique(size_t idxClique)
{
  return this->largeClqs->moveClique(idxClique);
}

size_t CoinDynamicConflictGraph::nCliques() const
{
  return this->largeClqs->nCliques();
}

size_t CoinDynamicConflictGraph::cliqueSize(size_t idxClique) const
{
  return this->largeClqs->cliqueSize(idxClique);
}

const size_t *CoinDynamicConflictGraph::nodeCliques(size_t idxNode) const
{
  return this->largeClqs->nodeOccurrences(idxNode);
}

size_t CoinDynamicConflictGraph::nNodeCliques(size_t idxNode) const
{
  return this->largeClqs->nNodeOccurrences(idxNode);
}

void CoinDynamicConflictGraph::setDegree(size_t idxNode, size_t deg)
{
  this->degree_[idxNode] = deg;
}

void CoinDynamicConflictGraph::setModifiedDegree(size_t idxNode, size_t mdegree)
{
  this->modifiedDegree_[idxNode] = mdegree;
}

void CoinDynamicConflictGraph::processSmallCliquesNode(
  size_t node,
  const size_t scn[],
  const size_t nscn,
  const CoinCliqueList *smallCliques,
  char *iv)
{
  size_t newConf = 0;
  size_t prevConf = conflicts->rowSize(node);
  const size_t *oldConfs = conflicts->getRow(node);
  for (size_t j = 0; (j < prevConf); ++j)
    iv[oldConfs[j]] = true;
  iv[node] = true;

  // marking conflicts in large cliques
  for (size_t j = 0; (j < nNodeCliques(node)); ++j) {
    size_t idxClq = nodeCliques(node)[j];
    for (size_t k = 0; (k < cliqueSize(idxClq)); ++k) {
      iv[cliqueElements(idxClq)[k]] = true;
    }
  }

  for (size_t k = 0; (k < nscn); ++k) {
    size_t idxClq = scn[k];
    for (size_t j = 0; (j < smallCliques->cliqueSize(idxClq)); ++j) {
      size_t clqEl = smallCliques->cliqueElements(idxClq)[j];
      if (!iv[clqEl]) {
        ++newConf;
        iv[clqEl] = true;
        conflicts->fastAddNeighbor(node, clqEl);
      }
    }
  }

  // clearing iv
  iv[node] = false;
  for (size_t i = 0; (i < nDirectConflicts(node)); ++i)
    iv[directConflicts(node)[i]] = false;

  for (size_t j = 0; (j < nNodeCliques(node)); ++j) {
    size_t idxClq = nodeCliques(node)[j];
    for (size_t k = 0; (k < cliqueSize(idxClq)); ++k) {
      iv[cliqueElements(idxClq)[k]] = false;
    }
  }

  if (newConf >= 1)
    conflicts->sort(node);
}

size_t CoinDynamicConflictGraph::degree(const size_t node) const
{
  return degree_[node];
}

size_t CoinDynamicConflictGraph::modifiedDegree(const size_t node) const
{
  return modifiedDegree_[node];
}

size_t CoinDynamicConflictGraph::nTotalCliqueElements() const
{
  return this->largeClqs->totalElements();
}

size_t CoinDynamicConflictGraph::nTotalDirectConflicts() const
{
  return this->conflicts->totalElements();
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
