/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class for separating violated odd-cycles. It contains
 * a lifting module that tries to transform the odd-cycles
 * into odd-wheels.
 *
 * @file CglOddWheel.cpp
 * @brief Odd-wheel cut separator
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include <cstdio>
#include <cassert>

#include "CglOddWheel.hpp"
#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "CoinTime.hpp"
#include "CoinCutPool.hpp"
#include "CoinOddWheelSeparator.hpp"

#define ODDHWC_EPS 1e-6

size_t CglOddWheel::sepCuts = 0;
double CglOddWheel::sepTime = 0.0;

static void *xmalloc( const size_t size );

CglOddWheel::CglOddWheel(size_t extMethod) : cap_(0), extMethod_(extMethod) {
    idxs_ = NULL;
    idxMap_ = NULL;
    coefs_ = NULL;
    x_ = NULL;
    rc_ = NULL;
}

CglOddWheel::CglOddWheel(const CglOddWheel& rhs) {
    this->cap_ = rhs.cap_;
    this->extMethod_ = rhs.extMethod_;

    if (this->cap_ > 0) {
        this->idxs_ = (int*)xmalloc(sizeof(int) * this->cap_);
        this->idxMap_ = (int*)xmalloc(sizeof(int) * this->cap_);
        this->coefs_ = (double*)xmalloc(sizeof(double) * this->cap_);
        this->x_ = (double*)xmalloc(sizeof(double) * this->cap_ * 2);
        this->rc_ = (double*)xmalloc(sizeof(double) * this->cap_ * 2);
    } else {
        this->idxs_ = NULL;
        this->idxMap_ = NULL;
        this->coefs_ = NULL;
        this->x_ = NULL;
        this->rc_ = NULL;
    }
}

CglOddWheel::~CglOddWheel() {
    if (this->idxs_) {
        free(this->idxs_);
    }
    if (this->idxMap_) {
        free(this->idxMap_);
    }
    if (this->coefs_) {
        free(this->coefs_);
    }
    if (this->x_) {
        free(this->x_);
    }
    if (this->rc_) {
        free(this->rc_);
    }
}

void CglOddWheel::refreshSolver(OsiSolverInterface *solver) {
    solver->checkCGraph();
  // Get integer information
    solver->getColType(true);
}

CglCutGenerator * CglOddWheel::clone() const {
    return new CglOddWheel(*this);
}

void CglOddWheel::generateCuts( const OsiSolverInterface & si, OsiCuts & cs, const CglTreeInfo info ) {
    if (si.getNumCols() == 0 || si.getNumRows() == 0) {
        return;
    }

    double startSep = CoinCpuTime();
    const size_t numCols = si.getNumCols();
    const CoinConflictGraph *cgraph = si.getCGraph();

	if(numCols != cgraph->size() / 2) {
        fprintf(stderr, "Invalid conflict graph! Number of columns %ld ... in graph %ld\n", numCols, cgraph->size() / 2);
        exit(EXIT_FAILURE);
    }

    checkMemory(numCols);

    const double *colSol = si.getColSolution();
    const double *rCost = si.getReducedCost();
    for(size_t i = 0; i < numCols; i++) {
        x_[i] = colSol[i];
        rc_[i] = rCost[i];
        x_[i + numCols] = 1.0 - x_[i];
        rc_[i + numCols] = -rc_[i];
    }

    CoinOddWheelSeparator oddH(cgraph, x_, rc_, extMethod_);
    if (maxSeconds_ > 0.0)
        oddH.setMaxSeconds(maxSeconds_);
    CoinCutPool cutPool(x_, numCols);

    oddH.searchOddWheels();

    /* adding odd holes */
    for(size_t j = 0; j < oddH.numOddWheels(); j++) {
        const size_t *oddEl = oddH.oddHole(j);
        const size_t oddSize = oddH.oddHoleSize(j);
        double rhs = oddH.oddWheelRHS(j);

        if(oddSize < 5) {
            fprintf(stderr, "Invalid size of cut: %lu\n", oddSize);
            exit(EXIT_FAILURE);
        }

        int realSize = 0;
        bool duplicated = false;
        std::fill(idxMap_, idxMap_ + numCols, -1);

        for(size_t k = 0; k < oddSize; k++) {
            if(oddEl[k] < numCols) {
                if(idxMap_[oddEl[k]] == -1) {
                    idxMap_[oddEl[k]] = realSize;
                    idxs_[realSize] = oddEl[k];
                    coefs_[realSize] = 1.0;
                    realSize++;
                } else {
                    duplicated = true;
                    break;
                }
            }
            else {
                if(idxMap_[oddEl[k]-numCols] == -1) {
                    idxMap_[oddEl[k]-numCols] = realSize;
                    idxs_[realSize] = ((int)(oddEl[k] - numCols));
                    coefs_[realSize] = -1.0;
                    rhs -= 1.0;
                    realSize++;
                }
                else {
                    duplicated = true;
                    break;
                }
            }
        }

        if (duplicated) {
            continue;
        }

        const size_t centerSize = oddH.wheelCenterSize(j);
        const size_t *centerIdx = oddH.wheelCenter(j);
        if (centerSize && fabs(rhs) >= ODDHWC_EPS) {
            const double oldRhs = rhs;
            for (size_t k = 0; k < centerSize; k++) {
                if (centerIdx[k] < numCols) {
                    if (idxMap_[centerIdx[k]] == -1) {
                        idxMap_[centerIdx[k]] = realSize;
                        idxs_[realSize] = centerIdx[k];
                        coefs_[realSize] = oldRhs;
                        realSize++;
                    } else {
                        duplicated = true;
                        break;
                    }
                } else {
                    if (idxMap_[centerIdx[k] - numCols] == -1) {
                        idxMap_[centerIdx[k] - numCols] = realSize;
                        idxs_[realSize] = ((int) (centerIdx[k] - numCols));
                        coefs_[realSize] = -1.0 * oldRhs;
                        rhs = rhs - oldRhs;
                        realSize++;
                    } else {
                        duplicated = true;
                        break;
                    }
                }
            }

            if (duplicated) {
                continue;
            }
        }

        cutPool.add(idxs_, coefs_, realSize, rhs);
    }

    cutPool.removeNullCuts();

    const size_t numberRowCutsBefore = cs.sizeRowCuts();
    for(size_t i = 0; i < cutPool.numCuts(); i++) {
        osrc_.setRow(cutPool.cutSize(i) , cutPool.cutIdxs(i), cutPool.cutCoefs(i));
        osrc_.setUb(cutPool.cutRHS(i));
        cs.insertIfNotDuplicate(osrc_);
    }

    int numberRowCutsAfter = cs.sizeRowCuts();
    CglOddWheel::sepCuts += numberRowCutsAfter - numberRowCutsBefore;

    if(!info.inTree && ((info.options & 4) == 4 || ((info.options & 8) && !info.pass))) {
        numberRowCutsAfter = cs.sizeRowCuts();
        for(int i = numberRowCutsBefore; i < numberRowCutsAfter; i++) {
            cs.rowCutPtr(i)->setGloballyValid();
        }
    }

	CglOddWheel::sepTime += (CoinCpuTime() - startSep);
}

void CglOddWheel::checkMemory(const size_t newNumCols) {
    if (cap_ < newNumCols) {
        if (cap_ > 0) {
#ifdef DEBUGCG
            assert(idxs_);
            assert(idxMap_);
            assert(coefs_);
            assert(x_);
            assert(rc_);
#endif
            free(idxs_);
            free(idxMap_);
            free(coefs_);
            free(x_);
            free(rc_);
        }

        idxs_ = (int*)xmalloc(sizeof(int) * newNumCols);
        idxMap_ = (int*)xmalloc(sizeof(int) * newNumCols);
        coefs_ = (double*)xmalloc(sizeof(double) * newNumCols);
        x_ = (double*)xmalloc(sizeof(double) * newNumCols * 2);
        rc_ = (double*)xmalloc(sizeof(double) * newNumCols * 2);
        cap_ = newNumCols;
    }
}

static void *xmalloc( const size_t size ) {
    void *result = malloc( size );
    if (!result) {
        fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
        abort();
    }

    return result;
}

void CglOddWheel::setExtendingMethod(size_t extMethod) {
    if(extMethod > 2) {
        fprintf(stderr, "Invalid value for parameter extMethod (%lu).\n", extMethod);
        abort();
    }

    this->extMethod_ = extMethod;
}
