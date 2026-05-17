/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class that implements a conflict-based preprocessing.
 * It tries to extend set packing constraints considering
 * the conflict graph and using a greedy strategy.
 *
 * @file CglCliqueStrengthening.cpp
 * @brief Conflict-based preprocessing
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include <cfloat>
#include <cassert>
#include <algorithm>
#include "CglCliqueStrengthening.hpp"
#include "CoinStaticConflictGraph.hpp"
#include "CoinCliqueExtender.hpp"
#include "CglMessage.hpp"
#include "CoinTime.hpp"

#define CLQ_STR_EPS 1e-6
#define MAX_SIZE_CLIQUE_TO_BE_EXTENDED 256

static void *xmalloc( const size_t size );
static void *xcalloc( const size_t elements, const size_t size );

CliqueRows::CliqueRows(size_t linesToReserve, size_t nzsToReserve) {
  nRows_ = 0;
  starts_ = (size_t*)xmalloc(sizeof(size_t) * (linesToReserve + 1));
  starts_[0] = 0;
  rowIdx_ = (size_t*)xmalloc(sizeof(size_t) * linesToReserve);
  rowStatus_ = (CliqueRowStatus*)xmalloc(sizeof(CliqueRowStatus) * linesToReserve);
  elements_ = (size_t*)xmalloc(sizeof(size_t) * nzsToReserve);
}

CliqueRows::~CliqueRows() {
    free(starts_);
    free(elements_);
    free(rowIdx_);
    free(rowStatus_);
}

void CliqueRows::addRow(size_t nz, const size_t els[], size_t rowIdx, CliqueRowStatus status) {
  std::copy(els, els + nz, elements_ + starts_[nRows_]);
  rowIdx_[nRows_] = rowIdx;
  rowStatus_[nRows_] = status;
  nRows_++;
  starts_[nRows_] = starts_[nRows_ - 1] + nz;
}

const size_t* CliqueRows::row(size_t idxRow) const {
#ifdef DEBUGCG
  assert(idxRow < nRows_);
#endif
  return elements_ + starts_[idxRow];
}

size_t CliqueRows::origIdxRow(size_t idxRow) const {
#ifdef DEBUGCG
  assert(idxRow < nRows_);
#endif
  return rowIdx_[idxRow];
}

size_t CliqueRows::nz(size_t idxRow) const {
#ifdef DEBUGCG
  assert(idxRow < nRows_);
#endif
  return starts_[idxRow + 1] - starts_[idxRow];
}

CliqueRowStatus CliqueRows::status(size_t idxRow) const {
#ifdef DEBUGCG
  assert(idxRow < nRows_);
#endif
  return rowStatus_[idxRow];
}

void CliqueRows::setStatus(size_t idxRow, CliqueRowStatus status) const {
#ifdef DEBUGCG
  assert(idxRow < nRows_);
#endif
  rowStatus_[idxRow] = status;
}

size_t CliqueRows::rows() const {
  return nRows_;
}

CglCliqueStrengthening::CglCliqueStrengthening(OsiSolverInterface *model, CoinMessageHandler *dhandler) :
nExtended_(0), nDominated_(0), maxSeconds_(0.0), startTime_(0.0), handler_(NULL), defaultHandler_(true) {

  if (dhandler)
    this->passInMessageHandler(dhandler);
  else
    handler_ = new CoinMessageHandler();

  messages_ = CglMessage();
  model_ = model;
  model_->checkCGraph(dhandler);
  cgraph_ = model_->getCGraph();
  cliqueRows_ = NULL;
  posInClqRows_ = NULL;
  nColClqs_ = NULL;
  colClqs_ = NULL;
  dirtyRows_ = NULL;

  if (model_->getNumElements() > 0) {
    cliqueRows_ = new CliqueRows(model_->getNumRows(), model_->getNumElements());
    posInClqRows_ = (size_t*)xmalloc(sizeof(size_t) * model_->getNumRows());
    dirtyRows_ = (size_t*)xmalloc(sizeof(size_t) * model_->getNumRows());

    detectCliqueRows();
    fillCliquesByColumn();
  }
}

CglCliqueStrengthening::~CglCliqueStrengthening() {
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }

  if (cliqueRows_) {
    delete cliqueRows_;
    free(nColClqs_);
    free(colClqs_[0]);
    free(colClqs_);
    free(posInClqRows_);
    free(dirtyRows_);
  }
}

void CglCliqueStrengthening::detectCliqueRows() {
  const int numRows = model_->getNumRows();
  const int numCols = model_->getNumCols();
  const CoinPackedMatrix *cpmRow = model_->getMatrixByRow();
  const CoinBigIndex *starts = cpmRow->getVectorStarts();
  const double *Arhs = model_->getRightHandSide();
  const char *Asense = model_->getRowSense();
  const double *colLB = model_->getColLower();
  const double *colUB = model_->getColUpper();
  const char *colType = model_->getColType();
  const int *lengths = cpmRow->getVectorLengths();
  size_t *tmpRow = (size_t*)xmalloc(sizeof(size_t) * numCols);

  for (size_t i = 0; i < numRows; i++) {
    const size_t nz = lengths[i];
    const char sense = Asense[i];

    posInClqRows_[i] = numRows; //just initializing

    if (nz <= 1 || nz > MAX_SIZE_CLIQUE_TO_BE_EXTENDED) {
      continue;
    }

    if (sense != 'L' && sense != 'G') {
      continue;
    }

    const int *idxs = cpmRow->getIndices() + starts[i];
    const double *coefs = cpmRow->getElements() + starts[i];
    const double mult = (sense == 'G') ? -1.0 : 1.0;
    double rhs = mult * Arhs[i];

    bool testRow = true;
    double minCoef1 = std::numeric_limits< double >::max();
    double minCoef2 = std::numeric_limits< double >::max();
    for (size_t j = 0; j < nz; j++) {
      tmpRow[j] = idxs[j];

      double coefCol = coefs[j] * mult;
      const bool isBinary = (colType[tmpRow[j]] != 0) && (colLB[tmpRow[j]] == 1.0 || colLB[tmpRow[j]] == 0.0)
                            && (colUB[tmpRow[j]] == 0.0 || colUB[tmpRow[j]] == 1.0);

      if (!isBinary) {
        testRow = false;
        break;
      }

      if (coefCol <=- CLQ_STR_EPS) {
        tmpRow[j] += numCols;
        coefCol = -coefCol;
        rhs = rhs + coefCol;
      }

      if (coefCol + CLQ_STR_EPS <= minCoef1) {
        minCoef2 = minCoef1;
        minCoef1 = coefCol;
      } else if (coefCol + CLQ_STR_EPS <= minCoef2) {
        minCoef2 = coefCol;
      }
    }

    if (!testRow) {
      continue;
    }

    if (minCoef1 + minCoef2 >= rhs + CLQ_STR_EPS) {
      posInClqRows_[i] = cliqueRows_->rows();
      cliqueRows_->addRow(nz, tmpRow, i, NotDominated);
    }
  }

  free(tmpRow);
}

void CglCliqueStrengthening::fillCliquesByColumn() {
  size_t nElements = 0;
  const int numCols = model_->getNumCols();
  nColClqs_ = (size_t *) xcalloc(numCols * 2, sizeof(size_t));

  for (size_t i = 0; i < cliqueRows_->rows(); i++) {
    const size_t *clqEl = cliqueRows_->row(i);
    const size_t clqSize = cliqueRows_->nz(i);
#ifdef DEBUGCG
    assert(clqSize >= 2);
#endif

    for (size_t j = 0; j < clqSize; j++) {
      nColClqs_[clqEl[j]]++;
      nElements++;
    }
  }

  colClqs_ = (size_t **) xmalloc(sizeof(size_t *) * numCols * 2);
  colClqs_[0] = (size_t *) xmalloc(sizeof(size_t) * nElements);

  for (size_t i = 1; i < numCols * 2; i++) {
    colClqs_[i] = colClqs_[i - 1] + nColClqs_[i - 1];
    nColClqs_[i - 1] = 0;
  }

  nColClqs_[(2 * numCols) - 1] = 0;

  for (size_t i = 0; i < cliqueRows_->rows(); i++) {
    const size_t *clqEl = cliqueRows_->row(i);

    for (size_t j = 0; j < cliqueRows_->nz(i); j++) {
      const size_t col = clqEl[j];
      colClqs_[col][nColClqs_[col]++] = i;
    }
  }
}

void CglCliqueStrengthening::strengthenCliques(size_t extMethod) {
  nExtended_ = nDominated_ = 0;
  startTime_ = CoinWallclockTime();

  if (model_->getNumCols() == 0 || model_->getNumRows() == 0 || cliqueRows_->rows() == 0) {
    if (handler_->logLevel())
        handler_->message(CGL_PROCESS_CLQSTR, messages_) << nExtended_ << nDominated_ << CoinMessageEol;
    return;
  }

  CoinCliqueSet *newCliques = new CoinCliqueSet(4096, 32768);
  cliqueExtension(extMethod, newCliques);

  if (newCliques->nCliques() > 0) {
    // removing dominated cliques
    removeDominatedRows();

    // adding stronger cliques
    addStrongerCliques(newCliques);
  }

  if (handler_->logLevel())
      handler_->message(CGL_PROCESS_CLQSTR, messages_) << nExtended_ << nDominated_ << CoinMessageEol;
  delete newCliques;
}

void CglCliqueStrengthening::strengthenCliques(size_t n, const size_t rows[], size_t extMethod) {
  nExtended_ = nDominated_ = 0;
  startTime_ = CoinWallclockTime();

  if (model_->getNumCols() == 0 || model_->getNumRows() == 0 || cliqueRows_->rows() == 0) {
    if (handler_->logLevel())
        handler_->message(CGL_PROCESS_CLQSTR, messages_) << nExtended_ << nDominated_ << CoinMessageEol;
    return;
  }

  CoinCliqueSet *newCliques = new CoinCliqueSet(4096, 32768);
  cliqueExtension(extMethod, newCliques, n, rows);

  if (newCliques->nCliques() > 0) {
    // removing dominated cliques
    removeDominatedRows();

    // adding stronger cliques
    addStrongerCliques(newCliques);
  }

  if (handler_->logLevel())
    handler_->message(CGL_PROCESS_CLQSTR, messages_) << nExtended_ << nDominated_ << CoinMessageEol;
  delete newCliques;
}

void CglCliqueStrengthening::cliqueExtension(size_t extMethod, CoinCliqueSet *newCliques) {
  const int numCols = model_->getNumCols();
  bool *ivRow = (bool*)xcalloc(model_->getNumRows(), sizeof(bool));
  bool *ivCol = (bool*)xcalloc(numCols * 2, sizeof(bool));
  char name[256];

  // filling reduced costs
  double *rc = getReducedCost();

  // if reduced costs are not available, change the
  // extension method
  if (rc == NULL && (extMethod == 4 || extMethod == 5)) {
    extMethod = 2;
  }

  CoinCliqueExtender clqe(cgraph_, extMethod, rc);
  clqe.setMaxCandidates(512);

  for (size_t i = 0; i < cliqueRows_->rows(); i++) {
    if (maxSeconds_ > 0.0 && (i & 63) == 0) {
      if (CoinWallclockTime() - startTime_ >= maxSeconds_)
        break;
    }
    const size_t clqOrigRowIdx = cliqueRows_->origIdxRow(i);
    const size_t *clqIdx = cliqueRows_->row(i);
    const size_t clqSize = cliqueRows_->nz(i);

#ifdef DEBUGCG
    bool goodClique = CoinCliqueList::validateClique(cgraph_, clqIdx, clqSize);
    if (!goodClique)
      continue;
#endif

    if (cliqueRows_->status(i) == Dominated) {
        continue;
    }

    bool extended = clqe.extendClique(clqIdx, clqSize);

    if (extended) {
      cliqueRows_->setStatus(i, Dominated);
      nExtended_++;

      const size_t lastClq = clqe.nCliques() - 1;
      const bool inserted = newCliques->insertIfNotDuplicate(clqe.getCliqueSize(lastClq), clqe.getClique(lastClq));

      if (inserted) {
        checkDominance(clqe.getClique(lastClq), clqe.getCliqueSize(lastClq), ivRow, ivCol);
        sprintf(name, "%s_ext", model_->getRowName(clqOrigRowIdx).c_str());
        rowClqNames_.push_back(name);
      }
    }
  }

  // freeing memory
  if (rc) {
    free(rc);
  }
  free(ivRow);
  free(ivCol);
}

void CglCliqueStrengthening::cliqueExtension(size_t extMethod, CoinCliqueSet *newCliques, size_t n, const size_t rows[]) {
  const int numCols = model_->getNumCols();
  bool *ivRow = (bool*)xcalloc(model_->getNumRows(), sizeof(bool));
  bool *ivCol = (bool*)xcalloc(numCols * 2, sizeof(bool));
  char name[256];

  // filling reduced costs
  double *rc = getReducedCost();

  // if reduced costs are not available, change the
  // extension method
  if (rc == NULL) {
    extMethod = 2;
  }

  CoinCliqueExtender clqe(cgraph_, extMethod, rc);
  clqe.setMaxCandidates(512);

  for (size_t i = 0; i < n; i++) {
    if (maxSeconds_ > 0.0 && (i & 63) == 0) {
      if (CoinWallclockTime() - startTime_ >= maxSeconds_)
        break;
    }
    const size_t rowIdx = posInClqRows_[rows[i]];
    const size_t clqOrigRowIdx = cliqueRows_->origIdxRow(rowIdx);
    const size_t *clqIdx = cliqueRows_->row(rowIdx);
    const size_t clqSize = cliqueRows_->nz(rowIdx);

#ifdef DEBUGCG
    CoinCliqueList::validateClique(cgraph_, clqIdx, clqSize);
#endif

    if (cliqueRows_->status(rowIdx) == Dominated) {
        continue;
    }

    bool extended = clqe.extendClique(clqIdx, clqSize);

    if (extended) {
      cliqueRows_->setStatus(rowIdx, Dominated);
      nExtended_++;

      const size_t lastClq = clqe.nCliques() - 1;
      const bool inserted = newCliques->insertIfNotDuplicate(clqe.getCliqueSize(lastClq), clqe.getClique(lastClq));

      if (inserted) {
        checkDominance(clqe.getClique(lastClq), clqe.getCliqueSize(lastClq), ivRow, ivCol);
        sprintf(name, "%s_ext", model_->getRowName(clqOrigRowIdx).c_str());
        rowClqNames_.push_back(name);
      }
    }
  }

  // freeing memory
  if (rc) {
    free(rc);
  }
  free(ivRow);
  free(ivCol);
}

double* CglCliqueStrengthening::getReducedCost() {
  double *rc = NULL;

  if (model_->isProvenOptimal()) {
    const int numCols = model_->getNumCols();
    const double *redCost = model_->getReducedCost();
    rc = (double*)xmalloc(sizeof(double) * numCols * 2);

      for (size_t i = 0; i < numCols; i++) {
        rc[i] = redCost[i];
        rc[i + numCols] = -rc[i];
      }
  }

  return rc;
}

void CglCliqueStrengthening::checkDominance(const size_t *extClqEl, size_t extClqSize, bool *ivRow, bool *ivCol) {
#ifdef DEBUGCG
  for (size_t i = 0; i < cliqueRows_->rows(); i++) {
    assert(!ivRow[i]);
  }
  for (size_t i = 0; i < model_->getNumCols() * 2; i++) {
    assert(!ivCol[i]);
  }
#endif

  for (size_t i = 0; i < extClqSize; i++) {
    ivCol[extClqEl[i]] = true;
  }

  size_t nDirty = 0;
  for (size_t i = 0; i < extClqSize; i++) {
    size_t col = extClqEl[i];

    // checking cliques where column col appear
    for (size_t j = 0; j < nColClqs_[col]; j++) {
      const size_t clqRowIdx = colClqs_[col][j];

      // skipping already dominated or already tested rows
      if (cliqueRows_->status(clqRowIdx) == Dominated || ivRow[clqRowIdx]) {
        continue;
      }

      ivRow[clqRowIdx] = true;
      dirtyRows_[nDirty++] = clqRowIdx;

      const size_t *clqEl = cliqueRows_->row(clqRowIdx);
      const size_t clqNZ = cliqueRows_->nz(clqRowIdx);
      bool dominates = true;

      for (size_t k = 0; k < clqNZ; k++) {
        if (!ivCol[clqEl[k]]) {
          dominates = false;
          break;
        }
      }

      if (dominates) {
        cliqueRows_->setStatus(clqRowIdx, Dominated);
      }
    }
  }

  //clearing ivCol
  for (size_t i = 0; i < extClqSize; i++) {
    ivCol[extClqEl[i]] = false;
  }

  //clearing ivRow — only entries that were touched
  for (size_t i = 0; i < nDirty; i++) {
    ivRow[dirtyRows_[i]] = false;
  }
}

void CglCliqueStrengthening::removeDominatedRows() {
  int *toRemove = (int*)xmalloc(sizeof(int) * model_->getNumRows());

  nDominated_ = 0;

  for (size_t i = 0; i < cliqueRows_->rows(); i++) {
    if (cliqueRows_->status(i) == Dominated) {
      toRemove[nDominated_++] = (int)cliqueRows_->origIdxRow(i);
    }
  }

  if (nDominated_ > 0) {
    model_->deleteRows(nDominated_, toRemove);
  }

  free(toRemove);
}

void CglCliqueStrengthening::addStrongerCliques(const CoinCliqueSet *newCliques) {
  char name[256];
  const int numCols = model_->getNumCols();
  const size_t nCliques = newCliques->nCliques();
  size_t numVars = 0;

  int *nrIdx = (int*)xmalloc(sizeof(int) * newCliques->totalElements());
  int *idxMap = (int*)xmalloc(sizeof(int) * numCols);//controls duplicated indexes (var and complement)
  double *nrCoef = (double*)xmalloc(sizeof(double) * newCliques->totalElements());
  CoinBigIndex *nrStart = (CoinBigIndex*)xmalloc(sizeof(CoinBigIndex) * (nCliques + 1)); nrStart[0] = 0;
  double *nrLB = (double*)xmalloc(sizeof(double) * nCliques);
  double *nrUB = (double*)xmalloc(sizeof(double) * nCliques);
  bool feasible = true;
  for (size_t ic = 0; ic < nCliques; ic++) {
    const size_t extClqSize = newCliques->cliqueSize(ic);
    const size_t *extClqEl = newCliques->cliqueElements(ic);
    double rhs = 1.0;
    size_t duplicated = 0;

    std::fill(idxMap, idxMap + numCols, -1);

    for (size_t i = 0; i < extClqSize; i++) {
      if (extClqEl[i] < numCols) {
        if(idxMap[extClqEl[i]] == -1) {
          idxMap[extClqEl[i]] = numVars;
          nrIdx[numVars] = (int)extClqEl[i];
          nrCoef[numVars] = 1.0;
          numVars++;
        } else {
          nrCoef[idxMap[extClqEl[i]]] += 1.0;
          assert(nrCoef[idxMap[extClqEl[i]]] == 0.0);
          duplicated++;
        }
      } else {
        rhs -= 1.0;
        if(idxMap[extClqEl[i]-numCols] == -1) {
          idxMap[extClqEl[i]-numCols] = numVars;
          nrIdx[numVars] = ((int)extClqEl[i] - numCols);
          nrCoef[numVars] = -1.0;
          numVars++;
        } else {
          nrCoef[idxMap[extClqEl[i]-numCols]] -= 1.0;
          assert(nrCoef[idxMap[extClqEl[i]-numCols]] == 0.0);
          duplicated++;
        }
      }
    }

    if (duplicated>1) {
      // infeasible - I think - but better than assert anyway
      feasible = false;
      break;
    }
    if(duplicated == 1) {
      CoinBigIndex last = nrStart[ic];
      rhs = 0.0;
      for(CoinBigIndex k = nrStart[ic]; k < numVars; k++) {
        assert(nrCoef[k] == -1.0 || nrCoef[k] == 0.0 || nrCoef[k] == 1.0);
        if(nrCoef[k] == -1.0 || nrCoef[k] == 1.0) {
          nrIdx[last] = nrIdx[k];
          nrCoef[last] = nrCoef[k];
          last++;
          if (nrCoef[k] == -1.0) {
            rhs -= 1.0;
          }
        }
      }
      numVars = last;

      sprintf(name, "%s_dup", rowClqNames_[ic].c_str());
      rowClqNames_[ic] = name;
    }

    nrLB[ic] = -DBL_MAX;
    nrUB[ic] = rhs;
    nrStart[ic + 1] = numVars;
  }
  if (feasible) {
    // adding rows
    const int lastOrigIdx = model_->getNumRows();
    model_->addRows(nCliques, nrStart, nrIdx, nrCoef, nrLB, nrUB);

    // setting names
    for (int i = 0; i < (int)nCliques; i++) {
      model_->setRowName(lastOrigIdx + i, rowClqNames_[i]);
    }
  } else {
    // if infeasible add a bad row
    model_->addRow(0,NULL,NULL,1.0,0.0);
  }

  free(nrIdx);
  free(nrCoef);
  free(nrLB);
  free(nrUB);
  free(nrStart);
  free(idxMap);
}

static void *xmalloc( const size_t size ) {
  void *result = malloc( size );
  if (!result) {
    fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
    abort();
  }

  return result;
}

static void *xcalloc( const size_t elements, const size_t size ) {
  void *result = calloc( elements, size );
  if (!result) {
    fprintf(stderr, "No more memory available. Trying to callocate %zu bytes.", size * elements);
    abort();
  }

  return result;
}

// Pass in Message handler (not deleted at end)
void CglCliqueStrengthening::passInMessageHandler(CoinMessageHandler *handler) {
  if (defaultHandler_)
    delete handler_;
  defaultHandler_ = false;
  handler_ = handler;
}

// Set language
void CglCliqueStrengthening::newLanguage(CoinMessages::Language language) {
  messages_ = CglMessage(language);
}
