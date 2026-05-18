/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class for separating violated odd-cycles. It contains
 * a lifting module that tries to transform the odd-cycles
 * into odd-wheels.
 *
 * @file CoinOddWheelSeparator.cpp
 * @brief Odd-cycle cut separator
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <limits>
#include "CoinOddWheelSeparator.hpp"
#include "CoinConflictGraph.hpp"
#include "CoinShortestPath.hpp"
#include "CoinTime.hpp"

#define ODDWHEEL_SEP_DEF_MIN_FRAC               0.001
#define ODDWHEEL_SEP_DEF_EPS                    1e-6
#define ODDWHEEL_SEP_DEF_MAX_RC                 100.0
#define ODDWHEEL_SEP_DEF_MIN_VIOL               0.02
#define ODDWHEEL_SEP_DEF_MAX_WHEEL_CENTERS      ((size_t)256)

struct CompareCost {
    explicit CompareCost(const double *costs) { this->costs_ = costs; }

    bool operator () (const size_t &i, const size_t &j) {
        if(fabs(costs_[i] - costs_[j]) >= ODDWHEEL_SEP_DEF_EPS) {
            return costs_[i] + ODDWHEEL_SEP_DEF_EPS <= costs_[j];
        }

        return i < j;
    }

    const double *costs_;
};

CoinOddWheelSeparator::CoinOddWheelSeparator(const CoinConflictGraph *cgraph, const double *x, const double *rc, size_t extMethod) {
    const size_t cgSize = cgraph->size();

    cgraph_ = cgraph;
    x_ = x;
    rc_ = rc;
    icaCount_ = 0;
    icaIdx_ = std::vector<size_t>(cgSize);
    icaActivity_ = std::vector<double>(cgSize);
    fillActiveColumns();
    extMethod_ = extMethod;
    maxSeconds_ = 0.0;
    spf_ = NULL;

    if (icaCount_ > 4) {
        spArcStart_ = std::vector<size_t>((icaCount_ * 2) + 1);
        spArcCap_ = icaCount_ * 2;
        spArcTo_ = std::vector<size_t>(spArcCap_);
        spArcDist_ = std::vector<double>(spArcCap_);

        tmp_ = std::vector<size_t>(cgSize + 1);

        costs_ = std::vector<double>(cgSize);
        for (size_t i = 0; i < cgSize; i++) {
            if (x_[i] >= ODDWHEEL_SEP_DEF_EPS) {
                costs_[i] = (x_[i] * 1000.0);
            } else if (rc_[i] <= ODDWHEEL_SEP_DEF_MAX_RC) {
                costs_[i] = (1000000.0 + rc_[i]);
            } else {
                costs_[i] = std::numeric_limits<double>::max();
            }
        }

        iv_ = std::vector<char>(cgSize);
        iv2_ = std::vector<char>(cgSize);

        ohIdxs_ = std::vector<std::vector<size_t> >();
        ohIdxs_.reserve(icaCount_);
        wcIdxs_ = std::vector<std::vector<size_t> >(icaCount_);
    }
}

CoinOddWheelSeparator::~CoinOddWheelSeparator() {
    if (spf_) {
        delete spf_;
    }
}

void CoinOddWheelSeparator::searchOddWheels() {
    if (icaCount_ <= 4) {
        return;
    }

    const double startTime = (maxSeconds_ > 0.0) ? CoinGetTimeOfDay() : 0.0;

    if (!prepareGraph(startTime))
        return;

    // Check time every 128 nodes so we don't check too frequently on small
    // instances or too rarely on large ones.
    for (size_t i = 0; i < icaCount_; i++) {
        if (maxSeconds_ > 0.0 && (i & 127) == 0 && i > 0) {
            if (CoinGetTimeOfDay() - startTime >= maxSeconds_)
                break;
        }
        findOddHolesWithNode(i);
    }

    if (extMethod_ > 0) {
        //try to insert a wheel center
        for (size_t i = 0; i < ohIdxs_.size(); i++) {
            searchWheelCenter(i);
        }
    }
}

void CoinOddWheelSeparator::fillActiveColumns() {
    const size_t cgSize = cgraph_->size();

    icaCount_ = 0;

    if (cgSize <= 4) {
        return;
    }

    for (size_t j = 0; j < cgSize; j++) {
        if(cgraph_->degree(j) < 2) {
            continue;
        }

        if(x_[j] + ODDWHEEL_SEP_DEF_EPS <= ODDWHEEL_SEP_DEF_MIN_FRAC) {
            continue;
        }

#ifdef DEBUGCG
        assert(x_[j] >= -0.001 && x_[j] <= 1.001);
#endif

        icaIdx_[icaCount_] = j;
        icaActivity_[icaCount_] = 1001.0 - (1000.0 * x_[j]);
        icaCount_++;
    }

#ifdef DEBUGCG
    assert(icaCount_ <= cgraph_->size());
#endif
}

bool CoinOddWheelSeparator::prepareGraph(double startTime) {
    size_t idxArc = 0;
    const size_t nodes = icaCount_ * 2;

    //Conflicts: (x', y'')
    for (size_t i1 = 0; i1 < icaCount_; i1++) {
        // Check time every 64 outer iterations (each does icaCount_ inner steps).
        if (maxSeconds_ > 0.0 && (i1 & 63) == 0 && i1 > 0) {
            if (CoinGetTimeOfDay() - startTime >= maxSeconds_)
                return false;
        }
        spArcStart_[i1] = idxArc;
        const size_t idx1 = icaIdx_[i1];

        for (size_t i2 = 0; i2 < icaCount_; i2++) {
            const size_t idx2 = icaIdx_[i2];

            if (cgraph_->conflicting(idx1, idx2)) {
                if(idxArc + 1 > spArcCap_) {
                    spArcCap_ *= 2;
                    spArcTo_.resize(spArcCap_);
                    spArcDist_.resize(spArcCap_);
                }
                spArcTo_[idxArc] = icaCount_ + i2;
                spArcDist_[idxArc] = icaActivity_[i2];
                idxArc++;
            } // conflict found
        } // i2
    } // i1

    //Conflicts: (x'', y')
    for (size_t i1 = 0; i1 < icaCount_; i1++) {
        spArcStart_[icaCount_ + i1] = idxArc;

        for (size_t i2 = spArcStart_[i1]; i2 < spArcStart_[i1 + 1]; i2++) {
#ifdef DEBUGCG
            assert(spArcTo_[i2] >= icaCount_);
#endif
            const size_t arcTo = spArcTo_[i2] - icaCount_;
            const size_t arcDist = spArcDist_[i2];

            if(idxArc + 1 > spArcCap_) {
                spArcCap_ *= 2;
                spArcTo_.resize(spArcCap_);
                spArcDist_.resize(spArcCap_);
            }

            spArcTo_[idxArc] = arcTo;
            spArcDist_[idxArc] = arcDist;
            idxArc++;
        }
    }

    spArcStart_[icaCount_ * 2] = idxArc;
    spf_ = new CoinShortestPath(nodes, idxArc, spArcStart_.data(), spArcTo_.data(), spArcDist_.data());
    return true;
}

void CoinOddWheelSeparator::findOddHolesWithNode(size_t node) {
    const size_t dest = icaCount_ + node;

    spf_->find(node, dest);
    size_t oddSize = spf_->path(dest, tmp_.data());

#ifdef DEBUGCG
    assert(oddSize > 0);
#endif

    // first and last indexes are equal
    oddSize--;

    if (oddSize < 5) {
        return;
    }

    // translating indexes and checking for repeated entries
    for (size_t i = 0; i < oddSize; i++) {
        const size_t pos = tmp_[i] % icaCount_;
        tmp_[i] = icaIdx_[pos];

#ifdef DEBUGCG
        assert(pos < cgraph_->size());
        assert(icaIdx_[pos] < cgraph_->size());
#endif

        if (iv_[tmp_[i]]) { //repeated entry
            for (size_t j = 0; j <= i; j++) {
                iv_[tmp_[j]] = 0;
            }
            return;
        }

        iv_[tmp_[i]] = 1;
    }
    // clearing iv
    for (size_t i = 0; i < oddSize; i++) {
        iv_[tmp_[i]] = 0;
    }

    /* checking if it is violated */
    double lhs = 0.0;
    for (size_t i = 0; i < oddSize; i++) {
#ifdef DEBUGCG
        assert(tmp_[i] < cgraph_->size());
#endif
        lhs += x_[tmp_[i]];
    }
    const double rhs = floor(oddSize / 2.0);
    const double viol = lhs - rhs;
    if (viol + ODDWHEEL_SEP_DEF_EPS <= ODDWHEEL_SEP_DEF_MIN_VIOL) {
        return;
    }

    addOddHole(oddSize, tmp_);
}

bool CoinOddWheelSeparator::addOddHole(size_t nz, const std::vector<size_t> &idxs) {
    // checking for repeated entries
    if (alreadyInserted(nz, idxs)) {
        return false;
    }

    ohIdxs_.push_back(idxs);

    return true;
}

bool CoinOddWheelSeparator::alreadyInserted(size_t nz, const std::vector<size_t> &idxs) {
    bool repeated = false;

    for (size_t i = 0; i < nz; i++) {
        iv_[idxs[i]] = 1;
    }

    for (size_t idxOH = 0; idxOH < ohIdxs_.size(); idxOH++) {
        // checking size
        if (nz != ohIdxs_[idxOH].size()) {
            continue;
        }

        // checking indexes
        bool isEqual = true;
        const size_t *ohIdx = ohIdxs_[idxOH].data();
        for (size_t j = 0; j < nz; j++) {
            if (!iv_[ohIdx[j]]) {
                isEqual = false;
                break;
            }
        }
        if (isEqual) {
            repeated = true;
            break;
        }
    }

    // clearing iv
    for (size_t i = 0; i < nz; i++) {
        iv_[idxs[i]] = 0;
    }

    return repeated;
}

void CoinOddWheelSeparator::searchWheelCenter(size_t idxOH) {
#ifdef DEBUGCG
    assert(idxOH < ohIdxs_.size());
#endif

    const size_t *ohIdxs = ohIdxs_[idxOH].data();
    const size_t ohSize = ohIdxs_[idxOH].size();

#ifdef DEBUGCG
    assert(ohSize <= cgraph_->size());
#endif

    /* picking node with the smallest degree */
    size_t nodeSD = ohIdxs[0], minDegree = cgraph_->degree(ohIdxs[0]);
    iv_[ohIdxs[0]] = 1;
    for (size_t i = 1; i < ohSize; i++) {
        const size_t dg = cgraph_->degree(ohIdxs[i]);
        if (dg < minDegree) {
            minDegree = dg;
            nodeSD = ohIdxs[i];
        }

        iv_[ohIdxs[i]] = 1;
    }

    // generating candidates
    const std::pair<size_t, const size_t*> rescg = cgraph_->conflictingNodes(nodeSD, tmp_.data(), iv2_.data());
    size_t numCandidates = 0;
    for (size_t i = 0; i < rescg.first; i++) {
        const size_t node = rescg.second[i];

        //already inserted
        if (iv_[node]) {
            continue;
        }

        if (cgraph_->degree(node) < ohSize) {
            continue;
        }

        bool insert = true;
        for (size_t j = 0; j < ohSize; j++) {
            if (!cgraph_->conflicting(node, ohIdxs[j])) {
                insert = false;
                break;
            }
        }
        if (!insert) {
            continue;
        }

        //new candidate
        if (x_[node] >= ODDWHEEL_SEP_DEF_EPS || rc_[node] <= ODDWHEEL_SEP_DEF_MAX_RC) {
            tmp_[numCandidates++] = node;
        }
    }

    if (numCandidates != 0) {
    	size_t sizeWC = 0;

    	if (extMethod_ == 1) { //wheel center with only one variable
    		size_t bestCandidate = tmp_[0];
            double bestCost = costs_[tmp_[0]];
            for (size_t i = 1; i < numCandidates; i++) {
                if (costs_[tmp_[i]] + ODDWHEEL_SEP_DEF_EPS <= bestCost) {
                    bestCandidate = tmp_[i];
                    bestCost = costs_[tmp_[i]];
                }
            }
            tmp_[sizeWC++] = bestCandidate;
    	} else { //wheel center formed by a clique
    		assert(extMethod_ == 2);
    		const size_t n = numCandidates;
	        numCandidates = std::min(numCandidates, ODDWHEEL_SEP_DEF_MAX_WHEEL_CENTERS);
	        std::partial_sort(tmp_.begin(), tmp_.begin() + numCandidates, tmp_.begin() + n, CompareCost(costs_.data()));

	        for (size_t i = 0; i < numCandidates; i++) {
	            /* need to have conflict with all nodes in clique */
	            const size_t selected = tmp_[i];
	            bool insert = true;
	            for (size_t j = 0; j < sizeWC; j++) {
	                if (!cgraph_->conflicting(tmp_[j], selected)) {
	                    insert = false;
	                    break;
	                }
	            }
	            if (insert) {
	                tmp_[sizeWC++] = selected;
	            }
	        }
    	}

        wcIdxs_[idxOH] = std::vector<size_t>(tmp_.begin(), tmp_.begin() + sizeWC);
    }

    // clearing iv
    for (size_t i = 0; i < ohSize; i++) {
        iv_[ohIdxs[i]] = 0;
    }
}

const size_t* CoinOddWheelSeparator::oddHole(size_t idxOH) const {
#ifdef DEBUGCG
    assert(idxOH < ohIdxs_.size());
#endif

    return ohIdxs_[idxOH].data();
}

size_t CoinOddWheelSeparator::oddHoleSize(size_t idxOH) const {
#ifdef DEBUGCG
    assert(idxOH < ohIdxs_.size());
#endif

    return ohIdxs_[idxOH].size();
}

double CoinOddWheelSeparator::oddWheelRHS(size_t idxOH) const {
#ifdef DEBUGCG
    assert(idxOH < ohIdxs_.size());
#endif

    return floor(static_cast<double>(ohIdxs_[idxOH].size()) / 2.0);
}

const size_t* CoinOddWheelSeparator::wheelCenter(const size_t idxOH) const {
#ifdef DEBUGCG
    assert(idxOH < ohIdxs_.size());
#endif

    return wcIdxs_[idxOH].data();
}

size_t CoinOddWheelSeparator::wheelCenterSize(const size_t idxOH) const {
#ifdef DEBUGCG
    assert(idxOH < ohIdxs_.size());
#endif
    return wcIdxs_[idxOH].size();
}

