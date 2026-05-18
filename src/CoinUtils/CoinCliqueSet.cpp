/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * A class to store a set of cliques. It is an
 * extension of class CoinCliqueList.
 *
 * @file CoinCliqueSet.cpp
 * @brief Set of cliques
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include "CoinCliqueSet.hpp"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

const size_t CoinCliqueSet::sequence_[] = {37, 31, 29, 17, 13, 11, 7, 1};
const size_t CoinCliqueSet::nSequence_ = 8;
const size_t CoinCliqueSet::nBuckets_ = 8192;

#define INI_SPACE_BUCKETS 32

CoinCliqueSet::CoinCliqueSet(size_t _iniClqCap, size_t _iniClqElCap)
  : CoinCliqueList(_iniClqCap, _iniClqElCap)
  , hash_(std::vector<std::vector<size_t> >(nBuckets_)) {
    for (size_t i = 1; i < nBuckets_; i++) {
        hash_[i].reserve(INI_SPACE_BUCKETS);
    }
}

CoinCliqueSet::~CoinCliqueSet() {}

bool CoinCliqueSet::insertIfNotDuplicate(size_t size, const size_t *els) {
    std::vector<size_t> tmpClq = std::vector<size_t>(els, els + size);
    std::sort(tmpClq.begin(), tmpClq.end());

    size_t idxBucket = vectorHashCode(tmpClq);

#ifdef DEBUGCG
    assert(idxBucket >= 0);
    assert(idxBucket < nBuckets_);
#endif

    if (alreadyInserted(tmpClq, idxBucket)) {
        return false;
    }

    const size_t currCap = hash_[idxBucket].capacity();
    const size_t currSize = hash_[idxBucket].size();

    /*updating capacity*/
    if (currSize + 1 > currCap) {
        const size_t newCap = currCap * 2;
        hash_[idxBucket].reserve(newCap);
    }

    /*updating hash*/
    hash_[idxBucket].push_back(nCliques());

    /* inserting into CliqueList */
    addClique(size, tmpClq.data());

    return true;
}

size_t CoinCliqueSet::vectorHashCode(const std::vector<size_t> &els) {

    size_t size = els.size();
#ifdef DEBUGCG
    assert(size > 0);
#endif
    size_t code = (size * sequence_[0]);
    code += (els[0] * sequence_[1]);

    for (size_t i = 1; i < size; i++) {
        code += (sequence_[i % nSequence_] * els[i]);
    }

    code = (code % nBuckets_);

#ifdef DEBUGCG
    assert(code >= 0);
    assert(code < nBuckets_);
#endif

    return code;
}

bool CoinCliqueSet::alreadyInserted(const std::vector<size_t> &els, size_t hashCode) {
#ifdef DEBUGCG
    assert(hashCode >= 0);
    assert(hashCode < nBuckets_);
#endif

    size_t size = els.size();

    for(size_t i = 0; i < hash_[hashCode].size(); i++) {
        const size_t idxClq = hash_[hashCode][i];

#ifdef DEBUGCG
        assert(idxClq >= 0);
        assert(idxClq < nCliques());
#endif

        if (size != cliqueSize(idxClq)) {
            continue;
        }

        const size_t *clqEls = cliqueElements(idxClq);

        bool isEqual = true;
        for (size_t j = 0; j < size; j++) {
            if (els[j] != clqEls[j]) {
                isEqual = false;
                break;
            }
        }
        if (isEqual) {
            return true;
        }
    }

    return false;
}

