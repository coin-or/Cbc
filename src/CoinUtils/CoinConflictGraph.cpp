/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Abstract class for a Conflict Graph, see CoinStaticConflictGraph and
 * CoinDynamicConflictGraph for concrete implementations.
 *
 * @file CoinConflictGraph.cpp
 * @brief Abstract class for conflict graph
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include <cassert>
#include <cstring>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <climits>
#include <sstream>
#include <cmath>
#include <limits>
#include "CoinConflictGraph.hpp"
#include "CoinAdjacencyVector.hpp"
#include "CoinTime.hpp"

// Minimum row length required for the row to be stored as an explicit clique instead of
// being expanded into pairwise conflicts (tunable via setMinCliqueRow).
size_t CoinConflictGraph::minClqRow_ = 256;
size_t CoinConflictGraph::maxCliques_ = 1200;

CoinConflictGraph::CoinConflictGraph(size_t _size) {
    iniCoinConflictGraph(_size);
}

CoinConflictGraph::CoinConflictGraph(const CoinConflictGraph *other) {
    iniCoinConflictGraph(other);
}

double CoinConflictGraph::density() const {
    return density_;
}

size_t CoinConflictGraph::size() const {
    return this->size_;
}

CoinConflictGraph::~CoinConflictGraph() {
}

size_t CoinConflictGraph::minDegree() const {
    return minDegree_;
}

size_t CoinConflictGraph::maxDegree() const {
    return maxDegree_;
}

bool CoinConflictGraph::conflicting(size_t n1, size_t n2) const {
    if (n1 == n2) {
        return false;
    }

    size_t ndc1 = nDirectConflicts(n1);
    size_t ndc2 = nDirectConflicts(n2);
    const size_t *dc;
    size_t ndc, nodeToSearch;
    // checking direct conflicts
    if (ndc1 < ndc2) {
        ndc = ndc1;
        dc = directConflicts(n1);
        nodeToSearch = n2;
    } else {
        ndc = ndc2;
        dc = directConflicts(n2);
        nodeToSearch = n1;
    }

    if (std::binary_search(dc, dc + ndc, nodeToSearch))
        return true;

    return conflictInCliques(n1, n2);
}

void CoinConflictGraph::recomputeDegree() {
    double start = CoinCpuTime();
    this->nConflicts_ = 0;
    minDegree_ = std::numeric_limits<size_t>::max();
    maxDegree_ = std::numeric_limits<size_t>::min();

    const size_t maxDg = size_ - 1;

    for (size_t i = 0; (i < size_); ++i) {
        size_t dg = nDirectConflicts(i);

        // Approximate degree contribution from cliques:
        // instead of deduplicating element-by-element (expensive),
        // add (cliqueSize - 1) per clique. This may overcount when
        // nodes appear in multiple cliques or overlap with direct
        // conflicts, but degree is only used for prioritization.
        const size_t nnc = this->nNodeCliques(i);
        const size_t *nc = this->nodeCliques(i);
        for (size_t k = 0; (k < nnc); ++k) {
            dg += this->cliqueSize(nc[k]) - 1;
        }

        // cap at maximum possible degree
        if (dg > maxDg)
            dg = maxDg;

        setDegree(i, dg);
        setModifiedDegree(i, dg);
        minDegree_ = std::min(minDegree_, dg);
        maxDegree_ = std::max(maxDegree_, dg);
        nConflicts_ += dg;
    }

    density_ = (double) nConflicts_ / maxConflicts_;
    double secs = CoinCpuTime() - start;
    if (secs>1.0)
      printf("recompute degree took %.3f seconds!\n", secs);
}

void CoinConflictGraph::computeModifiedDegree() {
    if (!updateMDegree) {
        return;
    }

    std::vector<char> iv = std::vector<char>(size_);

    for (size_t i = 0; i < size_; i++) {
        const size_t ndc = nDirectConflicts(i);
        const size_t *dc = directConflicts(i);
        size_t mdegree = degree(i);

        iv[i] = 1;
        for (size_t k = 0; k < ndc; k++) {
            mdegree += degree(dc[k]);
            iv[dc[k]] = 1;
        }

        const size_t nnc = nNodeCliques(i);
        const size_t *nc = nodeCliques(i);
        for (size_t k = 0; k < nnc; k++) {
            const size_t *clqEls = this->cliqueElements(nc[k]);
            for (size_t l = 0; l < cliqueSize(nc[k]); l++) {
                if (!iv[clqEls[l]]) {
                    mdegree += degree(clqEls[l]);
                    iv[clqEls[l]] = 1;
                }
            }
        }

        setModifiedDegree(i, mdegree);

        //clearing iv
        iv[i] = 0;
        for (size_t k = 0; k < ndc; k++) {
            iv[dc[k]] = 0;
        }
        for (size_t k = 0; k < nnc; k++) {
            const size_t *clqEls = cliqueElements(nc[k]);
            for (size_t l = 0; l < cliqueSize(nc[k]); l++) {
                iv[clqEls[l]] = 0;
            }
        }
    }

    updateMDegree = false;
}

std::pair<size_t, const size_t *> CoinConflictGraph::conflictingNodes(size_t node, size_t *temp, char *iv) const {
    if (nNodeCliques(node)) {
        const size_t ndc = nDirectConflicts(node);
        const size_t *dc = directConflicts(node);

        // adding direct conflicts and after conflicts from cliques
        iv[node] = 1;
        for (size_t k = 0; k < ndc; k++) {
            temp[k] = dc[k];
            iv[dc[k]] = 1;
        }

        size_t nConf = ndc;

        // traversing node cliques
        for (size_t ic = 0; (ic < nNodeCliques(node)); ++ic) {
            size_t idxClq = nodeCliques(node)[ic];
            // elements of clique
            for (size_t j = 0; (j < cliqueSize(idxClq)); ++j) {
                const size_t neigh = cliqueElements(idxClq)[j];
                if (!iv[neigh]) {
                    temp[nConf++] = neigh;
                    iv[neigh] = 1;
                }
            }
        }

#ifdef DEBUGCG
        assert(nConf <= degree(node));
#endif

        // clearing iv
        iv[node] = 0;
        for (size_t i = 0; (i < nConf); ++i)
            iv[temp[i]] = 0;

        std::sort(temp, temp + nConf);
        return std::pair<size_t, const size_t *>(nConf, temp);
    } else {
#ifdef DEBUGCG
        assert(nDirectConflicts(node) <= degree(node));
#endif
        // easy, node does not appears on explicit cliques
        return std::pair<size_t, const size_t *>(nDirectConflicts(node), directConflicts(node));
    }
}

bool CoinConflictGraph::conflictInCliques(size_t n1, size_t n2) const {
    size_t nc1 = nNodeCliques(n1);
    size_t nc2 = nNodeCliques(n2);
    size_t nnc;
    size_t nodeToSearch;
    const size_t *clqList;
    if (nc1 < nc2) {
        nnc = nc1;
        clqList = nodeCliques(n1);
        nodeToSearch = n2;
    } else {
        nnc = nc2;
        clqList = nodeCliques(n2);
        nodeToSearch = n1;
    }

    for (size_t i = 0; i < nnc; ++i) {
        size_t idxClq = clqList[i];
        size_t cs = cliqueSize(idxClq);
        const size_t *ce = cliqueElements(idxClq);
        if (std::binary_search(ce, ce + cs, nodeToSearch))
            return true;
    }

    return false;
}

void CoinConflictGraph::iniCoinConflictGraph(size_t _size) {
    size_ = _size;
    nConflicts_ = 0;
    maxConflicts_ = ((double) size_) + ((double) size_) * ((double) size_);
    density_ = 0.0;
    minDegree_ = UINT_MAX;
    maxDegree_ = 0;
    updateMDegree = true;
}


void CoinConflictGraph::iniCoinConflictGraph(const CoinConflictGraph *other) {
    size_ = other->size_;
    nConflicts_ = other->nConflicts_;
    maxConflicts_ = other->maxConflicts_;
    density_ = other->density_;
    minDegree_ = other->minDegree_;
    maxDegree_ = other->maxDegree_;
    updateMDegree = other->updateMDegree;
    infeasibleImplications_ = infeasibleImplications_;
}

void CoinConflictGraph::setMinCliqueRow(size_t minClqRow) {
	CoinConflictGraph::minClqRow_ = minClqRow;
}

size_t CoinConflictGraph::getMinCliqueRow() {
	return CoinConflictGraph::minClqRow_;
}

void CoinConflictGraph::setMaxCliques(size_t maxClq) {
  CoinConflictGraph::maxCliques_ = maxClq;
}

size_t CoinConflictGraph::getMaxCliques() {
  return CoinConflictGraph::maxCliques_;
}

void CoinConflictGraph::printSummary() const {
    size_t numVertices = 0;
    size_t numEdges = 0;
    size_t minDegree = std::numeric_limits<size_t>::max(), maxDegree = std::numeric_limits<size_t>::min();
    const size_t numCols = size_ / 2;
    assert(size_ % 2 == 0);

    double avgDegree = 0.0, density = 0.0;
    double confsActiveVars = 0.0, confsCompVars = 0.0, mixedConfs = 0.0, trivialConflicts = 0.0;

    for (size_t i = 0; i < size_; i++) {
        const size_t dg = degree(i);

        minDegree = std::min(minDegree, dg);
        maxDegree = std::max(maxDegree, dg);
        numEdges += dg;
        numVertices++;
    }

    if (numEdges) {
        std::vector<size_t> neighs = std::vector<size_t>(size_);
        std::vector<char> iv = std::vector<char>(size_);

        avgDegree = ((double) numEdges) / ((double) numVertices);
        density = (2.0 * ((double) numEdges)) / (((double) numVertices) * (((double) numVertices) - 1.0));

        for (size_t i = 0; i < size_; i++) {
            const std::pair<size_t, const size_t*> rescg = conflictingNodes(i, neighs.data(), iv.data());

            for (size_t j = 0; j < rescg.first; j++) {
                const size_t vertexNeighbor = rescg.second[j];

                if (vertexNeighbor == i + numCols || vertexNeighbor + numCols == i) {
                    trivialConflicts += 1.0;
                    continue;
                }

                if (i < numCols && vertexNeighbor < numCols) {
                    confsActiveVars += 1.0;
                } else if (i >= numCols && vertexNeighbor >= numCols) {
                    confsCompVars += 1.0;
                } else {
                    mixedConfs += 1.0;
                }
            }
        }

        confsActiveVars = (confsActiveVars / ((double)numEdges));
        confsCompVars = (confsCompVars / ((double)numEdges));
        mixedConfs = (mixedConfs / ((double)numEdges));
        trivialConflicts = (trivialConflicts / ((double)numEdges));
    }

    printf("%ld;%ld;%lf;%ld;%ld;%lf;%lf;%lf;%lf;%lf;", numVertices, numEdges, density, minDegree, maxDegree,
            avgDegree, confsActiveVars, confsCompVars, mixedConfs, trivialConflicts);
}

const std::vector<CoinConflictGraph::BinaryBoundInfeasibility> &CoinConflictGraph::infeasibleImplications() const
{
    return infeasibleImplications_;
}

void CoinConflictGraph::registerBoundImplicationInfeasibility(const BinaryBoundInfeasibility &info)
{
    infeasibleImplications_.push_back(info);
}

#ifdef CGRAPH_DEEP_DIVE
#include <unordered_map>
#include "CoinColumnType.hpp"
void CoinConflictGraph::validateConflictGraphUsingFeasibleSolution(
    size_t numCols,
    const char *colTypes,
    const double *colLower,
    const double *colUpper,
    const std::vector<std::string> &colNames,
    const std::vector< std::pair< std::string, double > > &mipStart) const {
#define MAX_WARNINGS_COLUMNS 10
    static int nWarningsColumns = 0;
    if (!numCols || !colTypes || !colLower || !colUpper) {
        std::cerr << "CoinConflictGraph::validateConflictGraph: incomplete model data." << std::endl;
        return;
    }

    if (size_ != 2 * numCols) {
        std::cerr << "CoinConflictGraph::validateConflictGraph: conflict graph size (" << size_
                  << ") is incompatible with solver columns (" << numCols << ")." << std::endl;
        return;
    }

    // Validate that only binary variables should have conflicts recorded.
    for (size_t col = 0; col < numCols; ++col) {
        const bool isBinary = colTypes[col] == CoinColumnType::Binary;
        const size_t degActive = degree(col);
        const size_t degComplement = degree(col + numCols);
        if (!isBinary && (degActive > 0 || degComplement > 0)) {
            std::string varName;
            if (col < colNames.size() && !colNames[col].empty()) {
                varName = colNames[col];
            } else {
                std::ostringstream os;
                os << "col_" << col;
                varName = os.str();
            }
            if (nWarningsColumns  < MAX_WARNINGS_COLUMNS) {
                std::cerr << "Warning: variable '" << varName
                        << "' (index=" << col << ", type=" << colTypes[col]
                        << ", bounds=[" << colLower[col] << ',' << colUpper[col]
                        << "]) has conflicts but is not binary." << std::endl;
                nWarningsColumns++;
            }
        }
    }

    // Build map from column names to indices for quick lookup.
    std::unordered_map<std::string, size_t> nameToIndex;
    nameToIndex.reserve(numCols * 2);
    for (size_t col = 0; col < numCols; ++col) {
        if (col < colNames.size() && !colNames[col].empty()) {
            nameToIndex[colNames[col]] = col;
        }
    }

    struct StartInfo {
        size_t index;
        double value;
        std::string name;
    };

    std::vector<StartInfo> binaryAssignments;
    binaryAssignments.reserve(mipStart.size());
    std::unordered_map<size_t, size_t> indexToPos;

    nWarningsColumns = 0;

    const double tol = 1e-6;
    for (const auto &assignment : mipStart) {
        auto it = nameToIndex.find(assignment.first);
        if (it == nameToIndex.end()) {
            continue;
        }

        const size_t idx = it->second;
        if (colTypes[idx] != CoinColumnType::Binary) {
            continue;
        }

        const double rawValue = assignment.second;
        double normalizedValue = rawValue;
        if (std::fabs(rawValue) <= tol) {
            normalizedValue = 0.0;
        } else if (std::fabs(rawValue - 1.0) <= tol) {
            normalizedValue = 1.0;
        }

        if (!(normalizedValue == 0.0 || normalizedValue == 1.0)) {
            if (nWarningsColumns < MAX_WARNINGS_COLUMNS) {
                std::cerr << "Warning: mipStart assigns value " << rawValue << " to binary variable '"
                        << assignment.first << "' (index=" << idx << ", type=" << colTypes[idx]
                        << ", bounds=[" << colLower[idx] << ',' << colUpper[idx] << "]). Expected 0/1.\n";
                nWarningsColumns++;
            }
            continue;
        }

        const auto insRes = indexToPos.emplace(idx, binaryAssignments.size());
        if (insRes.second) {
            binaryAssignments.push_back({ idx, normalizedValue, assignment.first });
        } else {
            binaryAssignments[insRes.first->second].value = normalizedValue;
        }
    }

    nWarningsColumns = 0;

    // Verify that mipStart assignments do not contradict the conflict graph.
    const size_t nBinaryAssignments = binaryAssignments.size();
    for (size_t i = 0; i < nBinaryAssignments; ++i) {
        const StartInfo &vi = binaryAssignments[i];
        const size_t nodeI = (vi.value == 1.0) ? vi.index : vi.index + numCols;
        for (size_t j = i + 1; j < nBinaryAssignments; ++j) {
            const StartInfo &vj = binaryAssignments[j];
            const size_t nodeJ = (vj.value == 1.0) ? vj.index : vj.index + numCols;
            if (conflicting(nodeI, nodeJ)) {

              if (nWarningsColumns  < MAX_WARNINGS_COLUMNS) {
                std::cerr << "Warning: mipStart assigns incompatible binary variables: '"
                          << vi.name << "' (idx=" << vi.index << ", type=" << colTypes[vi.index]
                          << ", bounds=[" << colLower[vi.index] << ',' << colUpper[vi.index]
                          << "], value=" << vi.value << ") and '" << vj.name << "' (idx="
                          << vj.index << ", type=" << colTypes[vj.index] << ", bounds=["
                          << colLower[vj.index] << ',' << colUpper[vj.index] << "], value="
                          << vj.value << ").\n";
                nWarningsColumns++;
              }
            }
        }
    }
}
#endif //CGRAPH_DEEP_DIVE


/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
