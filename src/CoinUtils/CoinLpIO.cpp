// Last edit: 11/5/08
//
// Name:     CoinLpIO.cpp; Support for Lp files
// Author:   Francois Margot
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     12/28/03
//-----------------------------------------------------------------------------
// Copyright (C) 2003, Francois Margot, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinUtilsConfig.h"

#include <cmath>
#include <cfloat>
#include <cctype>
#include <cassert>
#include <string>
#include <cstdarg>

#include "CoinError.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinLpIO.hpp"
#include "CoinMpsIO.hpp"
#include "CoinFinite.hpp"
#include "CoinSort.hpp"

#define LPIO_MODIFY_MESSAGES 1
#if LPIO_MODIFY_MESSAGES == 0
#define CoinLpIOError(a,b,c,d,e) throw CoinError(a,b,c,d,e)
#else
#define CoinLpIOError(a,b,c,d,e) throwError(a,b,c,d,e)
#endif

//#define LPIO_DEBUG

/************************************************************************/

CoinLpIO::CoinLpIO()
  : problemName_(CoinStrdup(""))
  , defaultHandler_(true)
  , numberRows_(0)
  , numberColumns_(0)
  , numberElements_(0)
  , matrixByColumn_(NULL)
  , matrixByRow_(NULL)
  , quadraticObjective_(NULL)
  , quadraticList_()
  , rowlower_(NULL)
  , rowupper_(NULL)
  , collower_(NULL)
  , colupper_(NULL)
  , rhs_(NULL)
  , rowrange_(NULL)
  , rowsense_(NULL)
  , num_objectives_(0)
  , integerType_(NULL)
  , set_(NULL)
  , numberSets_(0)
  , fileName_(NULL)
  , infinity_(COIN_DBL_MAX)
  , epsilon_(1e-5)
  , numberAcross_(10)
  , decimals_(9)
  , wasMaximization_(0)
  , lineNumber_(0)
  , input_(NULL)
{
  for (int j = 0; j < MAX_OBJECTIVES; j++) {
    objective_[j] = NULL;
    objName_[j] = NULL;
    objectiveOffset_[j] = 0;
  }
  card_previous_names_[0] = 0;
  card_previous_names_[1] = 0;
  previous_names_[0] = NULL;
  previous_names_[1] = NULL;

  maxHash_[0] = 0;
  numberHash_[0] = 0;
  hash_[0] = NULL;
  names_[0] = NULL;
  maxHash_[1] = 0;
  numberHash_[1] = 0;
  hash_[1] = NULL;
  names_[1] = NULL;
  handler_ = new CoinMessageHandler();
  messages_ = CoinMessage();
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CoinLpIO::CoinLpIO(const CoinLpIO &rhs)
  : problemName_(CoinStrdup(""))
  , defaultHandler_(true)
  , numberRows_(0)
  , numberColumns_(0)
  , numberElements_(0)
  , matrixByColumn_(NULL)
  , matrixByRow_(NULL)
  , quadraticObjective_(NULL)
  , quadraticList_()
  , rowlower_(NULL)
  , rowupper_(NULL)
  , collower_(NULL)
  , colupper_(NULL)
  , rhs_(NULL)
  , rowrange_(NULL)
  , rowsense_(NULL)
  , integerType_(NULL)
  , set_(NULL)
  , numberSets_(0)
  , fileName_(CoinStrdup(""))
  , infinity_(COIN_DBL_MAX)
  , epsilon_(1e-5)
  , numberAcross_(10)
  , input_(NULL)
{
  num_objectives_ = rhs.num_objectives_;
  for (int j = 0; j < MAX_OBJECTIVES; j++) {
    objective_[j] = NULL;
    if (j < num_objectives_) {
      objName_[j] = CoinStrdup(rhs.objName_[j]);
    } else {
      objName_[j] = NULL;
    }
    objectiveOffset_[j] = 0;
  }
  card_previous_names_[0] = 0;
  card_previous_names_[1] = 0;
  previous_names_[0] = NULL;
  previous_names_[1] = NULL;
  maxHash_[0] = 0;
  numberHash_[0] = 0;
  hash_[0] = NULL;
  names_[0] = NULL;
  maxHash_[1] = 0;
  numberHash_[1] = 0;
  hash_[1] = NULL;
  names_[1] = NULL;

  if (rhs.rowlower_ != NULL || rhs.collower_ != NULL) {
    gutsOfCopy(rhs);
  }

  defaultHandler_ = rhs.defaultHandler_;

  if (defaultHandler_) {
    handler_ = new CoinMessageHandler(*rhs.handler_);
  } else {
    handler_ = rhs.handler_;
  }

  messages_ = CoinMessage();
}

void CoinLpIO::gutsOfCopy(const CoinLpIO &rhs)
{
  defaultHandler_ = rhs.defaultHandler_;

  if (rhs.matrixByRow_) {
    matrixByRow_ = new CoinPackedMatrix(*(rhs.matrixByRow_));
  }

  numberElements_ = rhs.numberElements_;
  numberRows_ = rhs.numberRows_;
  numberColumns_ = rhs.numberColumns_;
  decimals_ = rhs.decimals_;
  wasMaximization_ = rhs.wasMaximization_;
  lineNumber_ = rhs.lineNumber_;

  if (rhs.rowlower_) {
    rowlower_ = reinterpret_cast< double * >(malloc(numberRows_ * sizeof(double)));
    rowupper_ = reinterpret_cast< double * >(malloc(numberRows_ * sizeof(double)));
    memcpy(rowlower_, rhs.rowlower_, numberRows_ * sizeof(double));
    memcpy(rowupper_, rhs.rowupper_, numberRows_ * sizeof(double));
    rowrange_ = reinterpret_cast< double * >(malloc(numberRows_ * sizeof(double)));
    rowsense_ = reinterpret_cast< char * >(malloc(numberRows_ * sizeof(char)));
    rhs_ = reinterpret_cast< double * >(malloc(numberRows_ * sizeof(double)));
    memcpy(rowrange_, rhs.getRowRange(), numberRows_ * sizeof(double));
    memcpy(rowsense_, rhs.getRowSense(), numberRows_ * sizeof(char));
    memcpy(rhs_, rhs.getRightHandSide(), numberRows_ * sizeof(double));
  }

  if (rhs.collower_) {
    collower_ = reinterpret_cast< double * >(malloc(numberColumns_ * sizeof(double)));
    colupper_ = reinterpret_cast< double * >(malloc(numberColumns_ * sizeof(double)));
    memcpy(collower_, rhs.collower_, numberColumns_ * sizeof(double));
    memcpy(colupper_, rhs.colupper_, numberColumns_ * sizeof(double));
    for (int j = 0; j < num_objectives_; j++) {
      objective_[j] = reinterpret_cast< double * >(malloc(numberColumns_ * sizeof(double)));
      memcpy(objective_[j], rhs.objective_[j], numberColumns_ * sizeof(double));
    }
  }

  if (rhs.integerType_) {
    integerType_ = reinterpret_cast< char * >(malloc(numberColumns_ * sizeof(char)));
    memcpy(integerType_, rhs.integerType_, numberColumns_ * sizeof(char));
  }
  numberSets_ = rhs.numberSets_;
  if (numberSets_) {
    set_ = new CoinSet *[numberSets_];
    for (int j = 0; j < numberSets_; j++)
      set_[j] = new CoinSet(*rhs.set_[j]);
  }

  if (rhs.quadraticObjective_)
    quadraticObjective_ = new CoinPackedMatrix(*rhs.quadraticObjective_);

  free(fileName_);
  free(problemName_);
  fileName_ = CoinStrdup(rhs.fileName_);
  problemName_ = CoinStrdup(rhs.problemName_);
  numberHash_[0] = rhs.numberHash_[0];
  numberHash_[1] = rhs.numberHash_[1];
  maxHash_[0] = rhs.maxHash_[0];
  maxHash_[1] = rhs.maxHash_[1];
  infinity_ = rhs.infinity_;
  numberAcross_ = rhs.numberAcross_;
  for (int j = 0; j < num_objectives_; j++) {
    objectiveOffset_[j] = rhs.objectiveOffset_[j];
  }
  int section;

  for (section = 0; section < 2; section++) {
    if (numberHash_[section]) {
      char **names2 = rhs.names_[section];
      names_[section] = reinterpret_cast< char ** >(malloc(maxHash_[section] * sizeof(char *)));
      char **names = names_[section];
      int i;

      for (i = 0; i < numberHash_[section]; i++) {
        names[i] = CoinStrdup(names2[i]);
      }

      hash_[section] = new CoinHashLink[maxHash_[section]];
      std::memcpy(hash_[section], rhs.hash_[section], maxHash_[section] * sizeof(CoinHashLink));
    }
  }
}

CoinLpIO &
CoinLpIO::operator=(const CoinLpIO &rhs)
{
  if (this != &rhs) {
    gutsOfDestructor();

    if (rhs.rowlower_ != NULL || rhs.collower_ != NULL) {
      gutsOfCopy(rhs);
    }

    defaultHandler_ = rhs.defaultHandler_;

    if (defaultHandler_) {
      handler_ = new CoinMessageHandler(*rhs.handler_);
    } else {
      handler_ = rhs.handler_;
    }

    messages_ = CoinMessage();
  }

  return *this;
}

void CoinLpIO::gutsOfDestructor()
{
  freeAll();

  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
}

/************************************************************************/
CoinLpIO::~CoinLpIO()
{
  stopHash(0);
  stopHash(1);
  freeAll();
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
}

/************************************************************************/
void CoinLpIO::freePreviousNames(const int section)
{

  int j;

  if (previous_names_[section] != NULL) {
    for (j = 0; j < card_previous_names_[section]; j++) {
      free(previous_names_[section][j]);
    }
    free(previous_names_[section]);
  }
  previous_names_[section] = NULL;
  card_previous_names_[section] = 0;
} /* freePreviousNames */

/************************************************************************/
void CoinLpIO::freeAll()
{

  delete matrixByColumn_;
  matrixByColumn_ = NULL;
  delete matrixByRow_;
  matrixByRow_ = NULL;
  delete quadraticObjective_;
  quadraticObjective_ = NULL;
  free(rowupper_);
  rowupper_ = NULL;
  free(rowlower_);
  rowlower_ = NULL;
  free(colupper_);
  colupper_ = NULL;
  free(collower_);
  collower_ = NULL;
  free(rhs_);
  rhs_ = NULL;
  free(rowrange_);
  rowrange_ = NULL;
  free(rowsense_);
  rowsense_ = NULL;
  for (int j = 0; j < num_objectives_; j++) {
    free(objective_[j]);
    objective_[j] = NULL;
  }
  free(integerType_);
  integerType_ = NULL;
  for (int j = 0; j < numberSets_; j++)
    delete set_[j];
  delete[] set_;
  set_ = NULL;
  numberSets_ = 0;
  free(problemName_);
  problemName_ = NULL;
  free(fileName_);
  fileName_ = NULL;

  freePreviousNames(0);
  freePreviousNames(1);
  delete input_;
  input_ = NULL;
}

/*************************************************************************/
const char *CoinLpIO::getProblemName() const
{
  return problemName_;
}

void CoinLpIO::setProblemName(const char *name)
{
  free(problemName_);
  problemName_ = CoinStrdup(name);
}

/*************************************************************************/
int CoinLpIO::getNumCols() const
{
  return numberColumns_;
}

/*************************************************************************/
int CoinLpIO::getNumRows() const
{
  return numberRows_;
}

/*************************************************************************/
CoinBigIndex CoinLpIO::getNumElements() const
{
  return numberElements_;
}

/*************************************************************************/
const double *CoinLpIO::getColLower() const
{
  return collower_;
}

/*************************************************************************/
const double *CoinLpIO::getColUpper() const
{
  return colupper_;
}

/*************************************************************************/
const double *CoinLpIO::getRowLower() const
{
  return rowlower_;
}

/*************************************************************************/
const double *CoinLpIO::getRowUpper() const
{
  return rowupper_;
}

/*************************************************************************/
/** A quick inlined function to convert from lb/ub style constraint
    definition to sense/rhs/range style */
inline void
CoinLpIO::convertBoundToSense(const double lower, const double upper,
  char &sense, double &right,
  double &range) const
{
  range = 0.0;
  if (lower > -infinity_) {
    if (upper < infinity_) {
      right = upper;
      if (upper == lower) {
        sense = 'E';
      } else {
        sense = 'R';
        range = upper - lower;
      }
    } else {
      sense = 'G';
      right = lower;
    }
  } else {
    if (upper < infinity_) {
      sense = 'L';
      right = upper;
    } else {
      sense = 'N';
      right = 0.0;
    }
  }
}

/*************************************************************************/
const char *CoinLpIO::getRowSense() const
{
  if (rowsense_ == NULL) {
    int nr = numberRows_;
    rowsense_ = reinterpret_cast< char * >(malloc(nr * sizeof(char)));

    double dum1, dum2;
    int i;
    for (i = 0; i < nr; i++) {
      convertBoundToSense(rowlower_[i], rowupper_[i], rowsense_[i], dum1, dum2);
    }
  }
  return rowsense_;
}

/*************************************************************************/
const double *CoinLpIO::getRightHandSide() const
{
  if (rhs_ == NULL) {
    int nr = numberRows_;
    rhs_ = reinterpret_cast< double * >(malloc(nr * sizeof(double)));

    char dum1;
    double dum2;
    int i;
    for (i = 0; i < nr; i++) {
      convertBoundToSense(rowlower_[i], rowupper_[i], dum1, rhs_[i], dum2);
    }
  }
  return rhs_;
}

/*************************************************************************/
const double *CoinLpIO::getRowRange() const
{
  if (rowrange_ == NULL) {
    int nr = numberRows_;
    rowrange_ = reinterpret_cast< double * >(malloc(nr * sizeof(double)));
    std::fill(rowrange_, rowrange_ + nr, 0.0);

    char dum1;
    double dum2;
    int i;
    for (i = 0; i < nr; i++) {
      convertBoundToSense(rowlower_[i], rowupper_[i], dum1, dum2, rowrange_[i]);
    }
  }
  return rowrange_;
}

/*************************************************************************/
const int CoinLpIO::getNumObjectives() const
{
  return num_objectives_;
}

/*************************************************************************/
const double *CoinLpIO::getObjCoefficients() const
{
  return objective_[0];
}

/*************************************************************************/
const double *CoinLpIO::getObjCoefficients(int j) const
{
  return objective_[j];
}

/*************************************************************************/
const CoinPackedMatrix *CoinLpIO::getMatrixByRow() const
{
  return matrixByRow_;
}

/*************************************************************************/
const CoinPackedMatrix *CoinLpIO::getMatrixByCol() const
{
  if (matrixByColumn_ == NULL && matrixByRow_) {
    matrixByColumn_ = new CoinPackedMatrix(*matrixByRow_);
    matrixByColumn_->reverseOrdering();
  }
  return matrixByColumn_;
}

/*************************************************************************/
const char *CoinLpIO::getObjName() const
{
  return objName_[0];
}

/*************************************************************************/
const char *CoinLpIO::getObjName(int j) const
{
  return objName_[j];
}

/*************************************************************************/
void CoinLpIO::checkRowNames()
{

  int i, nrow = getNumRows();

  if (numberHash_[0] != nrow + 1) {
    setDefaultRowNames();
    handler_->message(COIN_GENERAL_WARNING, messages_) << "### CoinLpIO::checkRowNames(): non distinct or missing row names or objective function name.\nNow using default row names."
                                                       << CoinMessageEol;
  }

  char const *const *rowNames = getRowNames();
  const char *rSense = getRowSense();
  char rName[256];

  // Check that row names and objective function name are all distinct,
  /// even after adding "_low" to ranged constraint names

  for (i = 0; i < nrow; i++) {
    if (rSense[i] == 'R') {
      sprintf(rName, "%s_low", rowNames[i]);
      if (findHash(rName, 0) != -1) {
        setDefaultRowNames();
        char printBuffer[512];
        sprintf(printBuffer, "### CoinLpIO::checkRowNames(): ranged constraint %d has a name %s identical to another constraint name or objective function name.\nUse getPreviousNames() to get the old row names.\nNow using default row names.", i, rName);
        handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
                                                           << CoinMessageEol;
        break;
      }
    }
  }
} /* checkRowNames */

/*************************************************************************/
void CoinLpIO::checkColNames()
{

  int ncol = getNumCols();

  if (numberHash_[1] != ncol) {
    setDefaultColNames();
    handler_->message(COIN_GENERAL_WARNING, messages_) << "### CoinLpIO::checkColNames(): non distinct or missing column names.\nNow using default column names."
                                                       << CoinMessageEol;
  }
} /* checkColNames */

/*************************************************************************/
void CoinLpIO::getPreviousRowNames(char const *const *prev,
  int *card_prev) const
{
  *card_prev = card_previous_names_[0];
  prev = previous_names_[0];
}

/*************************************************************************/
void CoinLpIO::getPreviousColNames(char const *const *prev,
  int *card_prev) const
{
  *card_prev = card_previous_names_[1];
  prev = previous_names_[1];
}

/*************************************************************************/
char const *const *CoinLpIO::getRowNames() const
{
  return names_[0];
}

/*************************************************************************/
char const *const *CoinLpIO::getColNames() const
{
  return names_[1];
}

/*************************************************************************/
const char *CoinLpIO::rowName(int index) const
{

  if ((names_[0] != NULL) && (index >= 0) && (index < numberRows_ + 1)) {
    return names_[0][index];
  } else {
    return NULL;
  }
}

/*************************************************************************/
const char *CoinLpIO::columnName(int index) const
{

  if ((names_[1] != NULL) && (index >= 0) && (index < numberColumns_)) {
    return names_[1][index];
  } else {
    return NULL;
  }
}

/*************************************************************************/
int CoinLpIO::rowIndex(const char *name) const
{

  if (!hash_[0]) {
    return -1;
  }
  return findHash(name, 0);
}

/*************************************************************************/
int CoinLpIO::columnIndex(const char *name) const
{

  if (!hash_[1]) {
    return -1;
  }
  return findHash(name, 1);
}

/************************************************************************/
double CoinLpIO::getInfinity() const
{
  return infinity_;
}

/************************************************************************/
void CoinLpIO::setInfinity(const double value)
{
  if (value >= 1.0e20) {
    infinity_ = value;
  } else {
    char str[8192];
    sprintf(str, "### ERROR: value: %f\n", value);
    throw CoinError(str, "setInfinity", "CoinLpIO", __FILE__, __LINE__);
  }
}

/************************************************************************/
double CoinLpIO::getEpsilon() const
{
  return epsilon_;
}

/************************************************************************/
void CoinLpIO::setEpsilon(const double value)
{
  if (value < 0.1) {
    epsilon_ = value;
  } else {
    char str[8192];
    sprintf(str, "### ERROR: value: %f\n", value);
    throw CoinError(str, "setEpsilon", "CoinLpIO", __FILE__, __LINE__);
  }
}

/************************************************************************/
int CoinLpIO::getNumberAcross() const
{
  return numberAcross_;
}

/************************************************************************/
void CoinLpIO::setNumberAcross(const int value)
{
  if (value > 0) {
    numberAcross_ = value;
  } else {
    char str[8192];
    sprintf(str, "### ERROR: value: %d\n", value);
    throw CoinError(str, "setNumberAcross", "CoinLpIO", __FILE__, __LINE__);
  }
}

/************************************************************************/
int CoinLpIO::getDecimals() const
{
  return decimals_;
}

/************************************************************************/
void CoinLpIO::setDecimals(const int value)
{
  if (value > 0) {
    decimals_ = value;
  } else {
    char str[8192];
    sprintf(str, "### ERROR: value: %d\n", value);
    throw CoinError(str, "setDecimals", "CoinLpIO", __FILE__, __LINE__);
  }
}

/************************************************************************/
double CoinLpIO::objectiveOffset() const
{
  return objectiveOffset_[0];
}

/************************************************************************/
double CoinLpIO::objectiveOffset(int j) const
{
  return objectiveOffset_[j];
}

/************************************************************************/
bool CoinLpIO::isInteger(int columnNumber) const
{
  const char *intType = integerType_;
  if (intType == NULL)
    return false;
  assert(columnNumber >= 0 && columnNumber < numberColumns_);
  if (intType[columnNumber] != 0)
    return true;
  return false;
}

/************************************************************************/
const char *CoinLpIO::integerColumns() const
{
  return integerType_;
}

/************************************************************************/
void CoinLpIO::setLpDataWithoutRowAndColNames(
  const CoinPackedMatrix &m,
  const double *collb, const double *colub,
  const double *obj_coeff,
  const char *is_integer,
  const double *rowlb, const double *rowub)
{

  setLpDataWithoutRowAndColNames(m, collb, colub, &obj_coeff, 1, is_integer,
    rowlb, rowub);
}

/************************************************************************/
void CoinLpIO::setLpDataWithoutRowAndColNames(
  const CoinPackedMatrix &m,
  const double *collb, const double *colub,
  const double *obj_coeff[MAX_OBJECTIVES],
  int num_objectives,
  const char *is_integer,
  const double *rowlb, const double *rowub)
{

  freeAll();
  problemName_ = CoinStrdup("");

  if (m.isColOrdered()) {
    matrixByRow_ = new CoinPackedMatrix();
    matrixByRow_->reverseOrderedCopyOf(m);
  } else {
    matrixByRow_ = new CoinPackedMatrix(m);
  }
  numberColumns_ = matrixByRow_->getNumCols();
  numberRows_ = matrixByRow_->getNumRows();

  rowlower_ = reinterpret_cast< double * >(malloc(numberRows_ * sizeof(double)));
  rowupper_ = reinterpret_cast< double * >(malloc(numberRows_ * sizeof(double)));
  collower_ = reinterpret_cast< double * >(malloc(numberColumns_ * sizeof(double)));
  colupper_ = reinterpret_cast< double * >(malloc(numberColumns_ * sizeof(double)));
  std::copy(rowlb, rowlb + numberRows_, rowlower_);
  std::copy(rowub, rowub + numberRows_, rowupper_);
  std::copy(collb, collb + numberColumns_, collower_);
  std::copy(colub, colub + numberColumns_, colupper_);
  num_objectives_ = num_objectives;
  for (int j = 0; j < num_objectives; j++) {
    objective_[j] = reinterpret_cast< double * >(malloc(numberColumns_ * sizeof(double)));
    std::copy(obj_coeff[j], obj_coeff[j] + numberColumns_, objective_[j]);
  }

  if (is_integer) {
    integerType_ = reinterpret_cast< char * >(malloc(numberColumns_ * sizeof(char)));
    std::copy(is_integer, is_integer + numberColumns_, integerType_);
  } else {
    integerType_ = 0;
  }

  if ((numberHash_[0] > 0) && (numberHash_[0] != numberRows_ + 1)) {
    stopHash(0);
  }
  if ((numberHash_[1] > 0) && (numberHash_[1] != numberColumns_)) {
    stopHash(1);
  }
} /* SetLpDataWithoutRowAndColNames */

/*************************************************************************/
void CoinLpIO::setDefaultRowNames()
{

  int i, nrow = getNumRows();
  char **defaultRowNames = reinterpret_cast< char ** >(malloc((nrow + 1) * sizeof(char *)));
  char buff[1024];

  for (i = 0; i < nrow; i++) {
    sprintf(buff, "cons%d", i);
    defaultRowNames[i] = CoinStrdup(buff);
  }
  sprintf(buff, "obj");
  defaultRowNames[nrow] = CoinStrdup(buff);

  stopHash(0);
  startHash(defaultRowNames, nrow + 1, 0);
  objName_[0] = CoinStrdup("obj");

  for (i = 0; i < nrow + 1; i++) {
    free(defaultRowNames[i]);
  }
  free(defaultRowNames);

} /* setDefaultRowNames */

/*************************************************************************/
void CoinLpIO::setDefaultColNames()
{

  int j, ncol = getNumCols();
  char **defaultColNames = reinterpret_cast< char ** >(malloc(ncol * sizeof(char *)));
  char buff[256];

  for (j = 0; j < ncol; j++) {
    sprintf(buff, "x%d", j);
    defaultColNames[j] = CoinStrdup(buff);
  }
  stopHash(1);
  startHash(defaultColNames, ncol, 1);

  for (j = 0; j < ncol; j++) {
    free(defaultColNames[j]);
  }
  free(defaultColNames);

} /* setDefaultColNames */

/*************************************************************************/
void CoinLpIO::setLpDataRowAndColNames(char const *const *const rownames,
  char const *const *const colnames)
{

  int nrow = getNumRows();
  int ncol = getNumCols();

  if (rownames != NULL) {
    if (are_invalid_names(rownames, nrow + 1, true)) {
      setDefaultRowNames();
      handler_->message(COIN_GENERAL_WARNING, messages_) << "### CoinLpIO::setLpDataRowAndColNames(): Invalid row names\nUse getPreviousNames() to get the old row names.\nNow using default row names."
                                                         << CoinMessageEol;
    } else {
      stopHash(0);
      startHash(rownames, nrow + 1, 0);
      objName_[0] = CoinStrdup(rownames[nrow]);
      checkRowNames();
    }
  } else {
    if (objName_[0] == NULL) {
      objName_[0] = CoinStrdup("obj");
    }
  }

  if (colnames != NULL) {
    if (are_invalid_names(colnames, ncol, false)) {
      setDefaultColNames();
      handler_->message(COIN_GENERAL_WARNING, messages_) << "### CoinLpIO::setLpDataRowAndColNames(): Invalid column names\nNow using default row names."
                                                         << CoinMessageEol;
    } else {
      stopHash(1);
      startHash(colnames, ncol, 1);
      checkColNames();
    }
  }
} /* setLpDataColAndRowNames */
// Load in SOS stuff
void CoinLpIO::loadSOS(int numberSets, const CoinSet *sets)
{
  if (numberSets_) {
    for (int i = 0; i < numberSets_; i++)
      delete set_[i];
    delete[] set_;
    set_ = NULL;
    numberSets_ = 0;
  }
  if (numberSets) {
    numberSets_ = numberSets;
    set_ = new CoinSet *[numberSets_];
    for (int i = 0; i < numberSets_; i++)
      set_[i] = new CoinSet(sets[i]);
  }
}

// Load in SOS stuff
void CoinLpIO::loadSOS(int numberSets, const CoinSet **sets)
{
  if (numberSets_) {
    for (int i = 0; i < numberSets_; i++)
      delete set_[i];
    delete[] set_;
    set_ = NULL;
    numberSets_ = 0;
  }
  if (numberSets) {
    numberSets_ = numberSets;
    set_ = new CoinSet *[numberSets_];
    for (int i = 0; i < numberSets_; i++)
      set_[i] = new CoinSet(*sets[i]);
  }
}

/************************************************************************/
void CoinLpIO::out_coeff(FILE *fp, const double v, const int print_1) const
{

  double lp_eps = getEpsilon();

  if (!print_1) {
    if (fabs(v - 1) < lp_eps) {
      return;
    }
    if (fabs(v + 1) < lp_eps) {
      fprintf(fp, " -");
      return;
    }
  }

  double frac = v - floor(v);

  if (frac < lp_eps) {
    fprintf(fp, " %.0f", floor(v));
  } else {
    if (frac > 1 - lp_eps) {
      fprintf(fp, " %.0f", floor(v + 0.5));
    } else {
      int decimals = getDecimals();
      char form[15];
      sprintf(form, " %%.%df", decimals);
      fprintf(fp, form, v);
    }
  }
} /* out_coeff */

/************************************************************************/
int CoinLpIO::writeLp(const char *filename, const double epsilon,
  const int numberAcross, const int decimals,
  const bool useRowNames)
{

  FILE *fp = NULL;
  fp = fopen(filename, "w");
  if (!fp) {
    char str[8192];
    sprintf(str, "### ERROR: unable to open file %s\n", filename);
    throw CoinError(str, "writeLP", "CoinLpIO", __FILE__, __LINE__);
  }
  int nerr = writeLp(fp, epsilon, numberAcross, decimals, useRowNames);
  fclose(fp);
  return (nerr);
}

/************************************************************************/
int CoinLpIO::writeLp(FILE *fp, const double epsilon,
  const int numberAcross, const int decimals,
  const bool useRowNames)
{

  setEpsilon(epsilon);
  setNumberAcross(numberAcross);
  setDecimals(decimals);
  return writeLp(fp, useRowNames);
}

/************************************************************************/
int CoinLpIO::writeLp(const char *filename, const bool useRowNames)
{
  FILE *fp = NULL;
  fp = fopen(filename, "w");
  if (!fp) {
    char str[8192];
    sprintf(str, "### ERROR: unable to open file %s\n", filename);
    throw CoinError(str, "writeLP", "CoinLpIO", __FILE__, __LINE__);
  }
  int nerr = writeLp(fp, useRowNames);
  fclose(fp);
  return (nerr);
}

/************************************************************************/
int CoinLpIO::writeLp(FILE *fp, const bool useRowNames)
{
  double lp_eps = getEpsilon();
  double lp_inf = getInfinity();
  int numberAcross = getNumberAcross();

  int i, cnt_print, loc_row_names = 0, loc_col_names = 0;
  CoinBigIndex j;
  char **prowNames = NULL, **pcolNames = NULL;

  const int *indices = matrixByRow_->getIndices();
  const double *elements = matrixByRow_->getElements();
  int ncol = getNumCols();
  int nrow = getNumRows();
  const double *collow = getColLower();
  const double *colup = getColUpper();
  const double *rowlow = getRowLower();
  const double *rowup = getRowUpper();
  const char *integerType = integerColumns();
  char const *const *rowNames = getRowNames();
  char const *const *colNames = getColNames();

  char buff[256];

  if (rowNames == NULL) {
    loc_row_names = 1;
    prowNames = reinterpret_cast< char ** >(malloc((nrow + 1) * sizeof(char *)));

    for (j = 0; j < nrow; j++) {
      sprintf(buff, "cons%d", j);
      prowNames[j] = CoinStrdup(buff);
    }
    prowNames[nrow] = CoinStrdup("obj");
    rowNames = prowNames;
  }

  if (colNames == NULL) {
    loc_col_names = 1;
    pcolNames = reinterpret_cast< char ** >(malloc(ncol * sizeof(char *)));

    for (j = 0; j < ncol; j++) {
      sprintf(buff, "x%d", j);
      pcolNames[j] = CoinStrdup(buff);
    }
    colNames = pcolNames;
  }

#ifdef LPIO_DEBUG
  printf("CoinLpIO::writeLp(): numberRows: %d numberColumns: %d\n",
    nrow, ncol);
#endif

  fprintf(fp, "\\Problem name: %s\n\n", getProblemName());
  fprintf(fp, "Minimize\n");

  for (int k = 0; k < num_objectives_; k++) {
    if (useRowNames) {
      fprintf(fp, "%s:", objName_[k]);
    }

    cnt_print = 0;
    for (j = 0; j < ncol; j++) {
      if ((cnt_print > 0) && (objective_[k][j] > lp_eps)) {
        fprintf(fp, " +");
      }
      if (fabs(objective_[k][j]) > lp_eps) {
        out_coeff(fp, objective_[k][j], 0);
        fprintf(fp, " %s", colNames[j]);
        cnt_print++;
        if (cnt_print % numberAcross == 0) {
          fprintf(fp, "\n");
        }
      }
    }

    if ((cnt_print > 0) && (objectiveOffset_[k] > lp_eps)) {
      fprintf(fp, " +");
    }
    if (fabs(objectiveOffset_[k]) > lp_eps) {
      out_coeff(fp, objectiveOffset_[k], 1);
      cnt_print++;
    }
    if ((cnt_print == 0) || (cnt_print % numberAcross != 0)) {
      fprintf(fp, "\n");
    }
  }
  // Quadratic objective
  if (quadraticObjective_) {
    const CoinBigIndex * start = quadraticObjective_->getVectorStarts();
    const int * column = quadraticObjective_->getIndices();
    const double * element = quadraticObjective_->getElements();
    fprintf (fp," + [ ");
    cnt_print++;
    if (cnt_print % numberAcross == 0) 
      fprintf(fp, "\n");
    bool first = true;
    double lp_eps = getEpsilon();
    char coeff[24];
    for (int iCol = 0;iCol < ncol;iCol++) {
      for (CoinBigIndex j=start[iCol]; j < start[iCol+1];j++) {
	int iRow = column[j];
	double coefficient = element[j];
	if (fabs(coefficient)<lp_eps)
	  continue;
	if (!first) {
	  if (coefficient > 0.0)
	    fprintf(fp,"+");
	} else {
	  first = false;
	}
	if (iCol != iRow)
	  coefficient *= 2.0;
	void CoinConvertDouble(int section, int formatType, double value, char outputValue[24]);
	CoinConvertDouble(0,0,coefficient,coeff);
	int n = strlen(coeff)-1;
	while (coeff[n]==' ')
	  n--;
	coeff[n+1] = '\0';
	if (iCol == iRow) {
	  fprintf(fp,"%s %s^2 ",coeff,colNames[iCol]);
	} else {
	  fprintf(fp,"%s %s * %s ",coeff,colNames[iCol],
		  colNames[iRow]);
	}
	cnt_print++;
	if (cnt_print % numberAcross == 0) 
	  fprintf(fp, "\n");
      }
    }
    fprintf (fp,"] /2\n");
  }

  fprintf(fp, "Subject To\n");

  //int cnt_out_rows = 0;

  for (i = 0; i < nrow; i++) {
    cnt_print = 0;

    if (useRowNames) {
      fprintf(fp, "%s: ", rowNames[i]);
    }
    //cnt_out_rows++;

    for (j = matrixByRow_->getVectorFirst(i);
         j < matrixByRow_->getVectorLast(i); j++) {
      if ((cnt_print > 0) && (elements[j] > lp_eps)) {
        fprintf(fp, " +");
      }
      if (fabs(elements[j]) > lp_eps) {
        out_coeff(fp, elements[j], 0);
        fprintf(fp, " %s", colNames[indices[j]]);
        cnt_print++;
        if (cnt_print % numberAcross == 0) {
          fprintf(fp, "\n");
        }
      }
    }

    if (rowup[i] - rowlow[i] < lp_eps) {
      fprintf(fp, " =");
      out_coeff(fp, rowlow[i], 1);
      fprintf(fp, "\n");
    } else {
      if (rowup[i] < lp_inf) {
        fprintf(fp, " <=");
        out_coeff(fp, rowup[i], 1);
        fprintf(fp, "\n");

        if (rowlower_[i] > -lp_inf) {

          cnt_print = 0;

          if (useRowNames) {
            fprintf(fp, "%s_low:", rowNames[i]);
          }
          //cnt_out_rows++;

          for (j = matrixByRow_->getVectorFirst(i);
               j < matrixByRow_->getVectorLast(i); j++) {
            if ((cnt_print > 0) && (elements[j] > lp_eps)) {
              fprintf(fp, " +");
            }
            if (fabs(elements[j]) > lp_eps) {
              out_coeff(fp, elements[j], 0);
              fprintf(fp, " %s", colNames[indices[j]]);
              cnt_print++;
              if (cnt_print % numberAcross == 0) {
                fprintf(fp, "\n");
              }
            }
          }
          fprintf(fp, " >=");
          out_coeff(fp, rowlow[i], 1);
          fprintf(fp, "\n");
        }

      } else {
        fprintf(fp, " >=");
        out_coeff(fp, rowlow[i], 1);
        fprintf(fp, "\n");
      }
    }
  }

#ifdef LPIO_DEBUG
  printf("CoinLpIO::writeLp(): Done with constraints\n");
#endif

  fprintf(fp, "Bounds\n");

  for (j = 0; j < ncol; j++) {
    if ((collow[j] > -lp_inf) && (colup[j] < lp_inf)) {
      out_coeff(fp, collow[j], 1);
      fprintf(fp, " <= %s <=", colNames[j]);
      out_coeff(fp, colup[j], 1);
      fprintf(fp, "\n");
    }
    if ((collow[j] == -lp_inf) && (colup[j] < lp_inf)) {
      fprintf(fp, "%s <=", colNames[j]);
      out_coeff(fp, colup[j], 1);
      fprintf(fp, "\n");
    }
    if ((collow[j] > -lp_inf) && (colup[j] == lp_inf)) {
      if (fabs(collow[j]) > lp_eps) {
        out_coeff(fp, collow[j], 1);
        fprintf(fp, " <= %s\n", colNames[j]);
      }
    }
    if (collow[j] == -lp_inf) {
      fprintf(fp, " %s Free\n", colNames[j]);
    }
  }

#ifdef LPIO_DEBUG
  printf("CoinLpIO::writeLp(): Done with bounds\n");
#endif

  bool semis = false;
  if (integerType != NULL) {
    int first_int = 1;
    cnt_print = 0;
    for (j = 0; j < ncol; j++) {
      if (integerType[j] == 1 || integerType[j] == 4) {

        if (first_int) {
          fprintf(fp, "Integers\n");
          first_int = 0;
        }

        fprintf(fp, "%s ", colNames[j]);
        cnt_print++;
        if (cnt_print % numberAcross == 0) {
          fprintf(fp, "\n");
        }
      }
      if (integerType[j] > 1)
        semis = true;
    }

    if (cnt_print % numberAcross != 0) {
      fprintf(fp, "\n");
    }
    if (semis) {
      int first_int = 1;
      cnt_print = 0;
      for (j = 0; j < ncol; j++) {
        if (integerType[j] > 2) {

          if (first_int) {
            fprintf(fp, "Semis\n");
            first_int = 0;
          }

          fprintf(fp, "%s ", colNames[j]);
          cnt_print++;
          if (cnt_print % numberAcross == 0) {
            fprintf(fp, "\n");
          }
        }
      }

      if (cnt_print % numberAcross != 0) {
        fprintf(fp, "\n");
      }
    }
  }

#ifdef LPIO_DEBUG
  printf("CoinLpIO::writeLp(): Done with integers\n");
#endif

  if (set_ != NULL) {
    fprintf(fp, "SOS\n");
    double lp_eps = getEpsilon();
    int decimals = getDecimals();
    char form[15];
    sprintf(form, "%%.%df", decimals);
    for (int iSet = 0; iSet < numberSets_; iSet++) {
      cnt_print = 0;
      const CoinSet *set = set_[iSet];
      // no space as readLp gets marginally confused
      fprintf(fp, "set%d:S%c::", iSet, '0' + set->setType());
      const int *which = set->which();
      const double *weights = set->weights();
      int numberEntries = set->numberEntries();
      for (j = 0; j < numberEntries; j++) {
        int iColumn = which[j];
        fprintf(fp, " %s:", colNames[iColumn]);
        // modified out_coeff (no leading space)
        double v = weights[j];
        double frac = v - floor(v);

        if (frac < lp_eps) {
          fprintf(fp, "%.0f", floor(v));
        } else {
          if (frac > 1 - lp_eps) {
            fprintf(fp, "%.0f", floor(v + 0.5));
          } else {
            fprintf(fp, form, v);
          }
        }
        cnt_print++;
        if (cnt_print % numberAcross == 0) {
          fprintf(fp, "\n");
        }
      }

      if (cnt_print % numberAcross != 0) {
        fprintf(fp, "\n");
      }
    }
  }

#ifdef LPIO_DEBUG
  printf("CoinLpIO::writeLp(): Done with SOS\n");
#endif

  fprintf(fp, "End\n");

  if (loc_row_names) {
    for (j = 0; j < nrow + 1; j++) {
      free(prowNames[j]);
    }
    free(prowNames);
  }

  if (loc_col_names) {
    for (j = 0; j < ncol; j++) {
      free(pcolNames[j]);
    }
    free(pcolNames);
  }
  return 0;
} /* writeLp */

/*************************************************************************/
int CoinLpIO::find_obj()
{

  char buff[1024];

  sprintf(buff, "aa");
  size_t lbuff = strlen(buff);

  while (((lbuff != 8) || (CoinStrNCaseCmp(buff, "minimize", 8) != 0)) && ((lbuff != 3) || (CoinStrNCaseCmp(buff, "min", 3) != 0)) && ((lbuff != 8) || (CoinStrNCaseCmp(buff, "maximize", 8) != 0)) && ((lbuff != 3) || (CoinStrNCaseCmp(buff, "max", 3) != 0))) {

    int x = fscanfLpIO(buff);
    lbuff = strlen(buff);

    if (x <= 0) {
      char str[8192];
      sprintf(str, "### ERROR: Unable to locate objective function\n");
      throw CoinError(str, "find_obj", "CoinLpIO", __FILE__, __LINE__);
    }
  }

  if (((lbuff == 8) && (CoinStrNCaseCmp(buff, "minimize", 8) == 0)) || ((lbuff == 3) && (CoinStrNCaseCmp(buff, "min", 3) == 0))) {
    return (1);
  }
  wasMaximization_ = 1; // say maximize (for qp)
  return (-1);
} /* find_obj */

/*************************************************************************/
int CoinLpIO::is_subject_to(const char *buff) const
{

  size_t lbuff = strlen(buff);

  if (((lbuff == 4) && (CoinStrNCaseCmp(buff, "s.t.", 4) == 0)) || ((lbuff == 3) && (CoinStrNCaseCmp(buff, "st.", 3) == 0)) || ((lbuff == 2) && (CoinStrNCaseCmp(buff, "st", 2) == 0))) {
    return (1);
  }
  if ((lbuff == 7) && (CoinStrNCaseCmp(buff, "subject", 7) == 0)) {
    return (2);
  }
  return (0);
} /* is_subject_to */

/*************************************************************************/
int CoinLpIO::first_is_number(const char *buff) const
{

  size_t pos;
  char str_num[] = "1234567890";

  pos = strcspn(buff, str_num);
  if (pos == 0) {
    return (1);
  }
  return (0);
} /* first_is_number */

/*************************************************************************/
int CoinLpIO::is_sense(const char *buff) const
{

  size_t pos;
  char str_sense[] = "<>=";

  pos = strcspn(buff, str_sense);
  if (pos == 0) {
    if (strcmp(buff, "<=") == 0) {
      return (0);
    }
    if (strcmp(buff, "=") == 0) {
      return (1);
    }
    if (strcmp(buff, ">=") == 0) {
      return (2);
    }
    char generalPrint[1000];
    sprintf(generalPrint,"### ERROR: CoinLpIO: is_sense(): string: %s", buff);
    warnError(generalPrint,lineNumber_);
  }
  return (-1);
} /* is_sense */

/*************************************************************************/
int CoinLpIO::is_free(const char *buff) const
{

  size_t lbuff = strlen(buff);

  if ((lbuff == 4) && (CoinStrNCaseCmp(buff, "free", 4) == 0)) {
    return (1);
  }
  return (0);
} /* is_free */

/*************************************************************************/
int CoinLpIO::is_inf(const char *buff) const
{
  size_t lbuff = strlen(buff);
#if 1
  switch (lbuff) {
  case 3:
    return (CoinStrNCaseCmp(buff, "inf", 3) == 0) ? 1 : 0;
  case 4:
    return (CoinStrNCaseCmp(buff, "-inf", 4) == 0) ? -1 : 0;
    break;
  case 8:
    return (CoinStrNCaseCmp(buff, "infinity", 8) == 0) ? 1 : 0;
    break;
  case 9:
    return (CoinStrNCaseCmp(buff, "-infinity", 9) == 0) ? -1 : 0;
    break;
  default:
    return 0;
  }
#else

  if ((lbuff == 3) && (CoinStrNCaseCmp(buff, "inf", 3) == 0)) {
    return (1);
  }
  return (0);
#endif
} /* is_inf */

/*************************************************************************/
int CoinLpIO::is_comment(const char *buff) const
{

  if ((buff[0] == '/') || (buff[0] == '\\')) {
    return (1);
  }
  return (0);
} /* is_comment */

/*************************************************************************/
void CoinLpIO::skip_comment(char *buff) const
{

  while (strcspn(buff, "\n") == strlen(buff)) { // end of line not read yet
    // keep going until in correct buffer
    while (bufferLength_ < 0) {
      if (fscanfLpIO(buff) == 0)
        throw("bad fgets");
    }
    // and throw away
    bufferPosition_ = bufferLength_;
    break;
  }
} /* skip_comment */

/*************************************************************************/
int CoinLpIO::is_invalid_name(const char *name,
  const bool ranged) const
{

  size_t pos, lname, valid_lname = 100;
  char str_valid[] = "1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ\"!#$%&(),.;?@_'`{}~";

  if (ranged) {
    valid_lname -= 4; // will add "_low" when writing the Lp file
  }

  if (name == NULL) {
    lname = 0;
  } else {
    lname = strlen(name);
  }
  if (lname < 1) {
    handler_->message(COIN_GENERAL_WARNING, messages_) << "### CoinLpIO::is_invalid_name(): Name is empty"
                                                       << CoinMessageEol;
    return (5);
  }
  if (lname > valid_lname) {
    char printBuffer[512];
    sprintf(printBuffer, "### CoinLpIO::is_invalid_name(): Name %s is too long",
      name);
    handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
                                                       << CoinMessageEol;
    return (1);
  }
  if (first_is_number(name)) {
    char printBuffer[512];
    sprintf(printBuffer, "### CoinLpIO::is_invalid_name(): Name %s should not start with a number", name);
    handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
                                                       << CoinMessageEol;
    return (2);
  }
  pos = strspn(name, str_valid);
  if (pos != lname) {
    char printBuffer[512];
    sprintf(printBuffer, "### CoinLpIO::is_invalid_name(): Name %s contains illegal character '%c'", name, name[pos]);
    handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
                                                       << CoinMessageEol;
    return (3);
  }

  if ((is_keyword(name)) || (is_free(name) || (is_inf(name)))) {
    return (4);
  }

  return (0);
} /* is_invalid_name */

/*************************************************************************/
int CoinLpIO::are_invalid_names(char const *const *const vnames,
  const int card_vnames,
  const bool check_ranged) const
{

  int i, invalid = 0, flag, nrows = getNumRows();
  bool is_ranged = 0;
  const char *rSense = getRowSense();

  if ((check_ranged) && (card_vnames != nrows + 1)) {
    char str[8192];
    sprintf(str, "### ERROR: card_vnames: %d   number of rows: %d\n",
      card_vnames, getNumRows());
    CoinLpIOError(str, "are_invalid_names", "CoinLpIO", __FILE__, __LINE__);
  }
  for (i = 0; i < card_vnames; i++) {

    if ((check_ranged) && (i < nrows) && (rSense[i] == 'R')) {
      is_ranged = true;
    } else {
      is_ranged = false;
    }
    flag = is_invalid_name(vnames[i], is_ranged);
    if (flag) {
      invalid = flag;
    }
  }
  return (invalid);
} /* are_invalid_names */

/*************************************************************************/
int CoinLpIO::read_monom_obj(double *coeff, char **name, int *cnt,
  char **obj_name, int *num_objectives, int *obj_starts)
{

  double mult;
  char buff[1024] = "aa", loc_name[1024], *start;
  int read_st = 0;
  int savePos = bufferPosition_;
  
  int x = fscanfLpIO(buff);

  if (x <= 0) {
    char str[8192];
    sprintf(str, "### ERROR: Unable to read objective function\n");
    CoinLpIOError(str, "read_monom_obj", "CoinLpIO", __FILE__, __LINE__);
  }

  if (buff[strlen(buff) - 1] == ':') {
    buff[strlen(buff) - 1] = '\0';

#ifdef LPIO_DEBUG
    printf("CoinLpIO: read_monom_obj(): obj_name: %s\n", buff);
#endif

    if (*num_objectives == MAX_OBJECTIVES) {
      char str[8192];
      sprintf(str, "### ERROR: Too many objective functions.\n");
      sprintf(str, "### ERROR: Change MAX_OBJECTIVES to larger number.\n");
      throw CoinError(str, "read_monom_obj", "CoinLpIO", __FILE__, __LINE__);
    }
    obj_name[*num_objectives] = CoinStrdup(buff);
    obj_starts[(*num_objectives)++] = *cnt;
    return (0);
  }

  if (*num_objectives == 0) {
    obj_starts[(*num_objectives)++] = *cnt;
  }

  read_st = is_subject_to(buff);

#ifdef LPIO_DEBUG
  printf("read_monom_obj: first buff: (%s)\n", buff);
#endif

  if (read_st > 0) {
    return (read_st);
  }

  start = buff;
  mult = 1;
  if (buff[0] == '+') {
    mult = 1;
    if (strlen(buff) == 1) {
      fscanfLpIO(buff);
      start = buff;
    } else {
      start = &(buff[1]);
    }
  }

  if (buff[0] == '-') {
    mult = -1;
    if (strlen(buff) == 1) {
      fscanfLpIO(buff);
      start = buff;
    } else {
      start = &(buff[1]);
    }
  }

  // If first - check not offset
  int offset = 0;
  if ((*cnt)==0 && first_is_number(start)) {
    // see if number to start
    char *temp = inputBuffer_+bufferPosition_;
    while (*temp==' ')
      temp++;
    if (*temp=='+'||*temp=='-') {
      offset = atof(start);
      setObjectiveOffset(-mult*offset);
      start = temp;
      offset = (*temp=='+') ? 1 : -1;
      if (!first_is_number(start)) {
	start--;
	start =temp;
	*start = '1';
      }
    }
  }
  if (first_is_number(start)) {
    coeff[*cnt] = atof(start);
    if ((*cnt)==0&&offset) {
      coeff[*cnt] *= offset;
      bufferPosition_++;
      while(inputBuffer_[bufferPosition_]==' ')
	bufferPosition_++;
    }
    sprintf(loc_name, "aa");
    fscanfLpIO(loc_name);
  } else {
    coeff[*cnt] = 1;
    strcpy(loc_name, start);
    // check if quadratic objective
    if (start[0]=='[' && strlen(start)==1) {
      // Save mult
      double overallMult = mult;
      if (wasMaximization_ > 0)
	overallMult = - mult;
      while (true) {
	double coefficient = 1.0;
	std::string x="x";
	std::string y="y";
	fscanfLpIO(buff);
	if (buff[0]==']') 
	  break; // finished
	start = buff;
	mult = 1;
	if (buff[0] == '+') {
	  mult = 1;
	  if (strlen(buff) == 1) {
	    fscanfLpIO(buff);
	    start = buff;
	  } else {
	    start = &(buff[1]);
	  }
	}

	if (buff[0] == '-') {
	  mult = -1;
	  if (strlen(buff) == 1) {
	    fscanfLpIO(buff);
	    start = buff;
	  } else {
	    start = &(buff[1]);
	  }
	}

	if (first_is_number(start)) {
	  coefficient = atof(start);
	  sprintf(buff, "aa");
	  fscanfLpIO(buff);
	  start = buff;
	}
	// see if ^ or * in string
	int gotxx = 0;
	if (strchr(start,'^')) {
	  gotxx = 1;
	  char * zap = strchr(start,'^');
	  *zap = '\0';
	  x = start;
	  start = zap+1;
	} else if (strchr(start,'*')) {
	  gotxx = 2;
	  char * zap = strchr(start,'*');
	  *zap = '\0';
	  x = start;
	  start = zap+1;
	} else {
	  x = start;
	}
	if (!gotxx) {
	  // next must be ^ or *
	  fscanfLpIO(buff);
	  int length = strlen(buff);
	  if (buff[0]=='^') {
	    if (length==1) {
	      fscanfLpIO(buff);
	      assert (buff[0]=='2' && strlen(buff)==1);
	    } else {
	      assert (buff[1]=='2' && length==2);
	    }
	    y = x;
	  } else if (buff[0]=='*') {
	    sprintf(loc_name, "bb");
	    fscanfLpIO(loc_name);
	    y = loc_name;
	  } else {
	    char str[100];
	    sprintf(str, "### ERROR:expected ^ or * in quadratic objective\n");
	    throw CoinError(str, "read_monom_obj", "CoinLpIO", __FILE__, __LINE__);
	  }
	} else if (gotxx==1) {
	  // ^
	  int length = strlen(start);
	  if (length==0) {
	    fscanfLpIO(buff);
	    assert (buff[0]=='2' && strlen(buff)==1);
	  } else if (length==1) {
	    assert (start[0]=='2');
	  }
	  y = x;
	} else {
	  // *
	  int length = strlen(start);
	  if (length==0) {
	    fscanfLpIO(buff);
	    y = buff;
	  } else {
	    y = start;
	  }
	}
	CoinLpQuadratic quadObj(coefficient*mult*overallMult,x,y);
	quadraticList_.push_back(quadObj);
      }
      fscanfLpIO(buff); // get Subject
      return 2; // expect To (of Subject To)
    }
  }

  read_st = is_subject_to(loc_name);

#ifdef LPIO_DEBUG
  printf("read_monom_obj: second buff: (%s)\n", buff);
#endif
  /* First one may just be objective offset
     - but I don't think it worth fixing as I
     think more natural at end */
  if (read_st > 0) {
    // offset is stored with sign other way round
    setObjectiveOffset(-mult * coeff[*cnt]+objectiveOffset());

#ifdef LPIO_DEBUG
    printf("read_monom_obj: objectiveOffset: %f\n", objectiveOffset_);
#endif

    return (read_st);
  }

  coeff[*cnt] *= mult;
  name[*cnt] = CoinStrdup(loc_name);

#ifdef LPIO_DEBUG
  printf("read_monom_obj: (%f)  (%s)\n", coeff[*cnt], name[*cnt]);
#endif

  (*cnt)++;

  return (read_st);
} /* read_monom_obj */

/*************************************************************************/
int CoinLpIO::read_monom_row(char *start_str,
  double *coeff, char **name,
  int cnt_coeff) const
{

  double mult;
  char buff[1024], loc_name[1024], *start;
  int read_sense = -1;

  sprintf(buff, "%s", start_str);
  read_sense = is_sense(buff);
  if (read_sense > -1) {
    return (read_sense);
  }

  start = buff;
  mult = 1;
  if (buff[0] == '+') {
    mult = 1;
    if (strlen(buff) == 1) {
      fscanfLpIO(buff);
      start = buff;
    } else {
      start = &(buff[1]);
    }
  }

  if (buff[0] == '-') {
    mult = -1;
    if (strlen(buff) == 1) {
      fscanfLpIO(buff);
      start = buff;
    } else {
      start = &(buff[1]);
    }
  }

  if (first_is_number(start)) {
    coeff[cnt_coeff] = atof(start);
    fscanfLpIO(loc_name);
  } else {
    coeff[cnt_coeff] = 1;
    strcpy(loc_name, start);
  }

  coeff[cnt_coeff] *= mult;
#ifdef KILL_ZERO_READLP
  if (fabs(coeff[cnt_coeff]) > epsilon_)
    name[cnt_coeff] = CoinStrdup(loc_name);
  else
    read_sense = -2; // effectively zero
#else
  name[cnt_coeff] = CoinStrdup(loc_name);
#endif

#ifdef LPIO_DEBUG
  printf("CoinLpIO: read_monom_row: (%f)  (%s)\n",
    coeff[cnt_coeff], name[cnt_coeff]);
#endif
  return (read_sense);
} /* read_monom_row */

/*************************************************************************/
void CoinLpIO::realloc_coeff(double **coeff, char ***colNames,
  int *maxcoeff) const
{

  *maxcoeff *= 5;

  *colNames = reinterpret_cast< char ** >(realloc((*colNames), (*maxcoeff + 1) * sizeof(char *)));
  *coeff = reinterpret_cast< double * >(realloc((*coeff), (*maxcoeff + 1) * sizeof(double)));

} /* realloc_coeff */

/*************************************************************************/
void CoinLpIO::realloc_row(char ***rowNames, CoinBigIndex **start, double **rhs,
  double **rowlow, double **rowup, int *maxrow) const
{

  *maxrow *= 5;
  *rowNames = reinterpret_cast< char ** >(realloc((*rowNames), (*maxrow + 1) * sizeof(char *)));
  *start = reinterpret_cast< CoinBigIndex * >(realloc((*start), (*maxrow + 1) * sizeof(CoinBigIndex)));
  *rhs = reinterpret_cast< double * >(realloc((*rhs), (*maxrow + 1) * sizeof(double)));
  *rowlow = reinterpret_cast< double * >(realloc((*rowlow), (*maxrow + 1) * sizeof(double)));
  *rowup = reinterpret_cast< double * >(realloc((*rowup), (*maxrow + 1) * sizeof(double)));

} /* realloc_row */

/*************************************************************************/
void CoinLpIO::realloc_col(double **collow, double **colup, char **is_int,
  int *maxcol) const
{
  *maxcol += 100 + (*maxcol)/10;
  *collow = reinterpret_cast< double * >(realloc((*collow), (*maxcol + 1) * sizeof(double)));
  *colup = reinterpret_cast< double * >(realloc((*colup), (*maxcol + 1) * sizeof(double)));
  *is_int = reinterpret_cast< char * >(realloc((*is_int), (*maxcol + 1) * sizeof(char)));
  // clean values
  double lp_inf = getInfinity();
  for (int i = (*maxcol) - 100; i < *maxcol; i++) {
    (*collow)[i] = 0;
    (*colup)[i] = lp_inf;
    (*is_int)[i] = 0;
  }

} /* realloc_col */

/*************************************************************************/
void CoinLpIO::read_row(char *buff,
  double **pcoeff, char ***pcolNames,
  int *cnt_coeff,
  int *maxcoeff,
  double *rhs, double *rowlow, double *rowup,
  int *cnt_row, double inf) const
{

  int read_sense = -1;
  char start_str[1024];

  sprintf(start_str, "%s", buff);

  while (read_sense < 0) {

    if ((*cnt_coeff) == (*maxcoeff)) {
      realloc_coeff(pcoeff, pcolNames, maxcoeff);
    }
    read_sense = read_monom_row(start_str,
      *pcoeff, *pcolNames, *cnt_coeff);
#ifdef KILL_ZERO_READLP
    if (read_sense != -2) // see if zero
#endif
      (*cnt_coeff)++;

    int x = fscanfLpIO(start_str);

    if (x <= 0) {
      char str[8192];
      sprintf(str, "### ERROR: Unable to read row monomial\n");
      CoinLpIOError(str, "read_monom_row", "CoinLpIO", __FILE__, __LINE__);
    }
  }
  (*cnt_coeff)--;

  rhs[*cnt_row] = atof(start_str);

  switch (read_sense) {
  case 0:
    rowlow[*cnt_row] = -inf;
    rowup[*cnt_row] = rhs[*cnt_row];
    break;
  case 1:
    rowlow[*cnt_row] = rhs[*cnt_row];
    rowup[*cnt_row] = rhs[*cnt_row];
    break;
  case 2:
    rowlow[*cnt_row] = rhs[*cnt_row];
    rowup[*cnt_row] = inf;
    break;
  default:
    break;
  }
  (*cnt_row)++;

} /* read_row */

/*************************************************************************/
int CoinLpIO::is_keyword(const char *buff) const
{

  size_t lbuff = strlen(buff);

  if (((lbuff == 5) && (CoinStrNCaseCmp(buff, "bound", 5) == 0)) || ((lbuff == 6) && (CoinStrNCaseCmp(buff, "bounds", 6) == 0))) {
    return (1);
  }

  if (((lbuff == 7) && (CoinStrNCaseCmp(buff, "integer", 7) == 0)) || ((lbuff == 8) && (CoinStrNCaseCmp(buff, "integers", 8) == 0))) {
    return (2);
  }

  if (((lbuff == 7) && (CoinStrNCaseCmp(buff, "general", 7) == 0)) || ((lbuff == 8) && (CoinStrNCaseCmp(buff, "generals", 8) == 0))) {
    return (2);
  }

  if (((lbuff == 6) && (CoinStrNCaseCmp(buff, "binary", 6) == 0)) || ((lbuff == 8) && (CoinStrNCaseCmp(buff, "binaries", 8) == 0))) {
    return (3);
  }

  if (((lbuff == 15) && (CoinStrNCaseCmp(buff, "semi-continuous", 15) == 0)) || ((lbuff == 4) && (CoinStrNCaseCmp(buff, "semi", 4) == 0)) || ((lbuff == 5) && (CoinStrNCaseCmp(buff, "semis", 5) == 0))) {
    return (4);
  }

  if ((lbuff == 3) && (CoinStrNCaseCmp(buff, "sos", 3) == 0)) {
    return (5);
  }

  if ((lbuff == 3) && (CoinStrNCaseCmp(buff, "end", 3) == 0)) {
    return (6);
  }

  return (0);

} /* is_keyword */

/*************************************************************************/
void CoinLpIO::readLp(const char *filename, const double epsilon)
{
  setEpsilon(epsilon);
  readLp(filename);
}

/*************************************************************************/
void CoinLpIO::readLp(const char *filename)
{
  delete input_;
  input_ = NULL;
  bool readable = false;
  // see if just .lp
  int length = strlen(filename);
  if ((length > 3 && !strncmp(filename + length - 3, ".lp", 3))) {
#if defined(_MSC_VER) || defined(__MINGW32__) || defined(__CYGWIN32__)
    // leave problemName_ as null
#else
    // find last /
    length--;
    while (length>=0) {
      if (filename[length]=='/')
	break;
      length--;
    }
    problemName_ = CoinStrdup(filename+length+1);
#endif
    FILE *fp = fopen(filename, "r");
    if (fp) {
      readable = true;
      input_ = new CoinPlainFileInput(fp);
    }
  } else if (strstr(filename, ".lp")) {
    std::string fname(filename);
    readable = fileCoinReadable(fname);
    if (readable)
      input_ = CoinFileInput::create(fname);
  } else if (!strcmp(filename,"-")) {
    input_ = new CoinPlainFileInput(stdin);
    readable = true;
  }
  if (!readable) {
    char str[8192];
    sprintf(str, "### ERROR: Unable to open file %s for reading\n", filename);
    throw CoinError(str, "readLp", "CoinLpIO", __FILE__, __LINE__);
  }
  readLp();
}

/*************************************************************************/
void CoinLpIO::readLp(FILE *fp, const double epsilon)
{
  setEpsilon(epsilon);
  readLp(fp);
}

/*************************************************************************/
void CoinLpIO::readLp(FILE *fp)
{
  delete input_;
  input_ = new CoinPlainFileInput(fp);
  readLp();
}
/*************************************************************************/
void CoinLpIO::readLp()
{

  int maxrow = 1000;
  int maxcoeff = 40000;
  double lp_eps = getEpsilon();
  double lp_inf = getInfinity();

  char buff[1024];
  bufferPosition_ = 0;
  bufferLength_ = 0;
  fakeBufferLength_ = 0;
  eofFound_ = false;
  lineNumber_ = 0;

  int objsense, cnt_coeff = 0, cnt_row = 0, cnt_obj = 0;
  int num_objectives = 0;
  char *objName[MAX_OBJECTIVES] = { NULL, NULL };
  int obj_starts[MAX_OBJECTIVES + 1];
  char **colNames = reinterpret_cast< char ** >(malloc((maxcoeff + 1) * sizeof(char *)));
  double *coeff = reinterpret_cast< double * >(malloc((maxcoeff + 1) * sizeof(double)));
  char **rowNames = reinterpret_cast< char ** >(malloc((maxrow + MAX_OBJECTIVES) * sizeof(char *)));
  CoinBigIndex *start = reinterpret_cast< CoinBigIndex * >(malloc((maxrow + MAX_OBJECTIVES) * sizeof(CoinBigIndex)));
  double *rhs = reinterpret_cast< double * >(malloc((maxrow + 1) * sizeof(double)));
  double *rowlow = reinterpret_cast< double * >(malloc((maxrow + 1) * sizeof(double)));
  double *rowup = reinterpret_cast< double * >(malloc((maxrow + 1) * sizeof(double)));

  int i;

  objsense = find_obj();

  int read_st = 0;
  while (!read_st) {
    read_st = read_monom_obj(coeff, colNames, &cnt_obj, objName, &num_objectives, obj_starts);

    if (cnt_obj == maxcoeff) {
      realloc_coeff(&coeff, &colNames, &maxcoeff);
    }
  }

  obj_starts[num_objectives] = cnt_obj;
  start[0] = cnt_obj;
  cnt_coeff = cnt_obj;

  if (read_st == 2) {
    int x = fscanfLpIO(buff);
    if (x <= 0)
      throw("bad fscanf");
    size_t lbuff = strlen(buff);

    if ((lbuff != 2) || (CoinStrNCaseCmp(buff, "to", 2) != 0)) {
      char str[8192];
      sprintf(str, "### ERROR: Can not locate keyword 'Subject To'\n");
      CoinLpIOError(str, "readLp", "CoinLpIO", __FILE__, __LINE__);
    }
  }

  fscanfLpIO(buff);

  while (!is_keyword(buff)) {
    if (buff[strlen(buff) - 1] == ':') {
      buff[strlen(buff) - 1] = '\0';

#ifdef LPIO_DEBUG
      printf("CoinLpIO::readLp(): rowName[%d]: %s\n", cnt_row, buff);
#endif

      rowNames[cnt_row] = CoinStrdup(buff);
      fscanfLpIO(buff);
    } else {
      char rname[15];
      sprintf(rname, "cons%d", cnt_row);
      rowNames[cnt_row] = CoinStrdup(rname);
    }
    read_row(buff,
      &coeff, &colNames, &cnt_coeff, &maxcoeff, rhs, rowlow, rowup,
      &cnt_row, lp_inf);
    fscanfLpIO(buff);
    start[cnt_row] = cnt_coeff;

    if (cnt_row == maxrow) {
      realloc_row(&rowNames, &start, &rhs, &rowlow, &rowup, &maxrow);
    }
  }

  numberRows_ = cnt_row;

  stopHash(1);
  startHash(colNames, cnt_coeff, 1);

  COINColumnIndex icol;
  int read_sense1, read_sense2;
  double bnd1 = 0, bnd2 = 0;

  int maxcol = numberHash_[1] + 100;

  double *collow = reinterpret_cast< double * >(malloc((maxcol + 1) * sizeof(double)));
  double *colup = reinterpret_cast< double * >(malloc((maxcol + 1) * sizeof(double)));
  char *is_int = reinterpret_cast< char * >(malloc((maxcol + 1) * sizeof(char)));
  int has_int = 0;

  for (i = 0; i < maxcol; i++) {
    collow[i] = 0;
    colup[i] = lp_inf;
    is_int[i] = 0;
  }

  int done = 0;
  // limit number of warnings
  int warnings = 0;
  while (!done) {
    switch (is_keyword(buff)) {

    case 1: /* Bounds section */
      fscanfLpIO(buff);

      while (is_keyword(buff) == 0) {

        read_sense1 = -1;
        read_sense2 = -1;
        int mult = 1;
        char *start_str = buff;

        if (buff[0] == '-' || buff[0] == '+') {
          mult = (buff[0] == '-') ? -1 : +1;
          if (strlen(buff) == 1) {
            fscanfLpIO(buff);
            start_str = buff;
          } else {
            start_str = &(buff[1]);
          }
        }

        int scan_sense = 0;
	// look for infinity
	int isInf = is_inf(start_str);
        if (!isInf && first_is_number(start_str)) {
          bnd1 = mult * atof(start_str);
          scan_sense = 1;
        } else if (isInf) {
	  bnd1 = isInf * mult * lp_inf;
	  scan_sense = 1;
        }
        if (scan_sense) {
          fscanfLpIO(buff);
          read_sense1 = is_sense(buff);
          if (read_sense1 < 0) {
            char str[8192];
            sprintf(str, "### ERROR: Bounds; expect a sense, get: %s\n", buff);
            CoinLpIOError(str, "readLp", "CoinLpIO", __FILE__, __LINE__);
          }
          fscanfLpIO(buff);
        }

        icol = findHash(buff, 1);
        if (icol < 0) {
	  warnings++;
	  if (warnings<50) {
	    char printBuffer[512];
	    sprintf(printBuffer, "### CoinLpIO::readLp(): Variable %s does not appear in objective function or constraints", buff);
	    handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
							       << CoinMessageEol;
	  }
          insertHash(buff, 1);
          icol = findHash(buff, 1);
          if (icol == maxcol) {
            realloc_col(&collow, &colup, &is_int, &maxcol);
          }
        }

        fscanfLpIO(buff);
        if (is_free(buff)) {
          collow[icol] = -lp_inf;
          fscanfLpIO(buff);
        } else {
          read_sense2 = is_sense(buff);
          if (read_sense2 > -1) {
            fscanfLpIO(buff);
            mult = 1;
            start_str = buff;

            if (buff[0] == '-' || buff[0] == '+') {
              mult = (buff[0] == '-') ? -1 : +1;
              if (strlen(buff) == 1) {
                fscanfLpIO(buff);
                start_str = buff;
              } else {
                start_str = &(buff[1]);
              }
            }
            if (first_is_number(start_str)) {
              bnd2 = mult * atof(start_str);
              fscanfLpIO(buff);
            } else {
              if (is_inf(start_str)) {
                bnd2 = mult * lp_inf;
                fscanfLpIO(buff);
              } else {
                char str[8192];
                sprintf(str, "### ERROR: Bounds; expect a number, get: %s\n",
                  buff);
                CoinLpIOError(str, "readLp", "CoinLpIO", __FILE__, __LINE__);
              }
            }
          }

          if ((read_sense1 > -1) && (read_sense2 > -1)) {
            if (read_sense1 != read_sense2) {
              char str[8192];
              sprintf(str, "### ERROR: Bounds; variable: %s read_sense1: %d  read_sense2: %d\n",
                buff, read_sense1, read_sense2);
              CoinLpIOError(str, "readLp", "CoinLpIO", __FILE__, __LINE__);
            } else {
              if (read_sense1 == 1) {
                if (fabs(bnd1 - bnd2) > lp_eps) {
                  char str[8192];
                  sprintf(str, "### ERROR: Bounds; variable: %s read_sense1: %d  read_sense2: %d  bnd1: %f  bnd2: %f\n",
                    buff, read_sense1, read_sense2, bnd1, bnd2);
                  CoinLpIOError(str, "readLp", "CoinLpIO", __FILE__, __LINE__);
                }
                collow[icol] = bnd1;
                colup[icol] = bnd1;
              }
              if (read_sense1 == 0) {
                collow[icol] = bnd1;
                colup[icol] = bnd2;
              }
              if (read_sense1 == 2) {
                colup[icol] = bnd1;
                collow[icol] = bnd2;
              }
            }
          } else {
            if (read_sense1 > -1) {
              switch (read_sense1) {
              case 0:
                collow[icol] = bnd1;
                break;
              case 1:
                collow[icol] = bnd1;
                colup[icol] = bnd1;
                break;
              case 2:
                colup[icol] = bnd1;
                break;
              }
            }
            if (read_sense2 > -1) {
              switch (read_sense2) {
              case 0:
                colup[icol] = bnd2;
                break;
              case 1:
                collow[icol] = bnd2;
                colup[icol] = bnd2;
                break;
              case 2:
                collow[icol] = bnd2;
                break;
              }
            }
          }
        }
      }
      break;

    case 2: /* Integers/Generals section */

      fscanfLpIO(buff);

      while (is_keyword(buff) == 0) {

        icol = findHash(buff, 1);

#ifdef LPIO_DEBUG
        printf("CoinLpIO::readLp(): Integer: colname: (%s)  icol: %d\n",
          buff, icol);
#endif

        if (icol < 0) {
          char printBuffer[512];
          sprintf(printBuffer, "### CoinLpIO::readLp(): Integer variable %s does not appear in objective function or constraints", buff);
          handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
                                                             << CoinMessageEol;
          insertHash(buff, 1);
          icol = findHash(buff, 1);
          if (icol == maxcol) {
            realloc_col(&collow, &colup, &is_int, &maxcol);
          }

#ifdef LPIO_DEBUG
          printf("CoinLpIO::readLp(): Integer: colname: (%s)  icol: %d\n",
            buff, icol);
#endif
        }
        if (is_int[icol] == 3)
          is_int[icol] = 4;
        else
          is_int[icol] = 1;
        has_int = 1;
        fscanfLpIO(buff);
      };
      break;

    case 3: /* Binaries section */

      fscanfLpIO(buff);

      while (is_keyword(buff) == 0) {

        icol = findHash(buff, 1);

#ifdef LPIO_DEBUG
        printf("CoinLpIO::readLp(): binary: colname: (%s)  icol: %d\n",
          buff, icol);
#endif

        if (icol < 0) {
          char printBuffer[512];
          sprintf(printBuffer, "### CoinLpIO::readLp(): Binary variable %s does not appear in objective function or constraints", buff);
          handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
                                                             << CoinMessageEol;
          insertHash(buff, 1);
          icol = findHash(buff, 1);
          if (icol == maxcol) {
            realloc_col(&collow, &colup, &is_int, &maxcol);
          }
#ifdef LPIO_DEBUG
          printf("CoinLpIO::readLp(): binary: colname: (%s)  icol: %d\n",
            buff, icol);
#endif
        }

        is_int[icol] = 1;
        has_int = 1;
        if (collow[icol] < 0) {
          collow[icol] = 0;
        }
        if (colup[icol] > 1) {
          colup[icol] = 1;
        }
        fscanfLpIO(buff);
      }
      break;
    case 4: /* Semis section */

      fscanfLpIO(buff);

      while (is_keyword(buff) == 0) {

        icol = findHash(buff, 1);

#ifdef LPIO_DEBUG
        printf("CoinLpIO::readLp(): Semi: colname: (%s)  icol: %d\n",
          buff, icol);
#endif

        if (icol < 0) {
          char printBuffer[512];
          sprintf(printBuffer, "### CoinLpIO::readLp(): Semi-continuous variable %s does not appear in objective function or constraints", buff);
          handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
                                                             << CoinMessageEol;
          insertHash(buff, 1);
          icol = findHash(buff, 1);
          if (icol == maxcol) {
            realloc_col(&collow, &colup, &is_int, &maxcol);
          }

#ifdef LPIO_DEBUG
          printf("CoinLpIO::readLp(): Semi-continuous: colname: (%s)  icol: %d\n",
            buff, icol);
#endif
        }
        if (is_int[icol] == 1)
          is_int[icol] = 4;
        else
          is_int[icol] = 3;
        has_int = 1;
        fscanfLpIO(buff);
      };
      break;

    case 5: /* sos section */
    {
      numberSets_ = 0;
      int maxSets = 10;
      CoinSet **set = new CoinSet *[maxSets];
      int maxEntries = 100;
      double *weights = new double[maxEntries];
      int *which = new int[maxEntries];
      char printBuffer[512];
      int numberBad = 0;
      bool maybeSetName;
      fscanfLpIO(buff);
      while (true) {
        int numberEntries = 0;
        int setType = -1;
        int goodLine = 2;
        bool endLine = false;
        bool gotStart = false;
        maybeSetName = false;
        while (!endLine) {
          if (is_keyword(buff) == 0) {
            // see if ::
            char *next = strstr(buff, "::");
            if (!gotStart) {
              if (!next) {
                // OK first time - may be name of set
                if (goodLine == 2) {
                  int length = strlen(buff);
                  if (buff[length - 1] == ':') {
                    goodLine = 1;
                    // merge to get rid of space
                    char temp[200];
                    strcpy(temp, buff);
                    fscanfLpIO(buff); // try again
                    strcat(temp, buff);
                    strcpy(buff, temp);
                    if (maybeSetName) {
                      // get rid of error
                      numberBad--;
                      maybeSetName = false;
                    }
                    continue;
                  } else {
                    goodLine = 0;
                  }
                } else {
                  goodLine = 0;
                }
              } else {
                // check nothing or set name
                if (strchr(buff, ':') < next || next == buff + 2) {
                  // be lazy and assume set name
                  // get type
                  next -= 2;
                  if (next >= buff && (!strncmp(next, "S1::", 4) || !strncmp(next, "S2::", 4))) {
                    setType = next[1] - '0';
                    gotStart = true;
                  } else {
                    // error
                    goodLine = 0;
                  }
                } else {
                  goodLine = 0;
                }
              }
            } else if (next) {
              // end of set
              endLine = true;
            }
          } else {
            endLine = true;
          }
          if (!goodLine) {
            endLine = true;
          }
          while (!endLine) {
            fscanfLpIO(buff);
            if (is_keyword(buff) == 0 && !strstr(buff, "::")) {
              // expect pair
              char *start_str = buff;
              char *next = strchr(start_str, ':');
              if (!next) {
                endLine = true;
                goodLine = 0;
              } else {
                *next = '\0';
                int iColumn = columnIndex(start_str);
                *next = ':';
                if (iColumn >= 0) {
                  int length = strlen(next + 1);
                  if (!length) {
                    // assume user put in spaces
                    fscanfLpIO(buff);
                    if (is_keyword(buff) != 0 || strstr(buff, "::")) {
                      goodLine = 0;
                    } else {
                      next = buff - 1;
                    }
                  }
                  double value = atof(next + 1);
                  if (numberEntries == maxEntries) {
                    maxEntries = 2 * maxEntries;
                    double *tempD = new double[maxEntries];
                    int *tempI = new int[maxEntries];
                    memcpy(tempD, weights, numberEntries * sizeof(double));
                    memcpy(tempI, which, numberEntries * sizeof(int));
                    delete[] weights;
                    weights = tempD;
                    delete[] which;
                    which = tempI;
                  }
                  weights[numberEntries] = value;
                  which[numberEntries++] = iColumn;
                } else {
                  // no match - assume start of next set
                  if (!numberBad) {
                    sprintf(printBuffer, "### CoinLpIO::readLp(): Variable %s not found or no weight",
                      buff);
                    handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
                                                                       << CoinMessageEol;
                    sprintf(printBuffer,
                      "Assuming next set name - consider no set names or use setnn:S1:: (no spaces)");
                    handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
                                                                       << CoinMessageEol;
                    maybeSetName = true;
                  }
                  numberBad++;
                  endLine = true;
                }
              }
            } else {
              endLine = true;
            }
          }
          if (!goodLine) {
            // print bad line
            setType = 3;
            sprintf(printBuffer, "### CoinLpIO::readLp(): bad SOS item %s", buff);
            handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
                                                               << CoinMessageEol;
            numberBad++;
            endLine = true;
            break;
          }
        }
        if (setType == 1 || setType == 2) {
          if (!numberEntries) {
            // empty set - error
            sprintf(printBuffer, "### CoinLpIO::readLp(): set %d is empty", numberSets_);
            handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
                                                               << CoinMessageEol;
          } else {
            CoinSort_2(weights, weights + numberEntries, which);
            double last = weights[0];
            for (int i = 1; i < numberEntries; i++) {
              if (fabs(last - weights[i]) < 1.0e-12) {
                setType = 3;
                break;
              }
            }
            if (setType != 3) {
              if (numberSets_ == maxSets) {
                maxSets *= 2;
                CoinSet **temp = new CoinSet *[maxSets];
                memcpy(temp, set, numberSets_ * sizeof(CoinSet *));
                delete[] set;
                set = temp;
              }
              CoinSosSet *newSet = new CoinSosSet(numberEntries, which, weights, setType);
              set[numberSets_++] = newSet;
            } else {
              sprintf(printBuffer, "### CoinLpIO::readLp(): set %d has duplicate weights",
                numberSets_);
              handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
                                                                 << CoinMessageEol;
            }
          }
        }
        if (is_keyword(buff) || (numberBad && !maybeSetName))
          break; // end
      }
      delete[] weights;
      delete[] which;
      if (numberSets_) {
        set_ = new CoinSet *[numberSets_];
        memcpy(set_, set, numberSets_ * sizeof(CoinSet *));
        delete[] set;
      }
    } break;

    case 6:
      done = 1;
      break;

    default:
      char str[8192];
      sprintf(str, "### ERROR: Lost while reading: (%s)\n", buff);
      throw CoinError(str, "readLp", "CoinLpIO", __FILE__, __LINE__);
      break;
    }
  }
  if (warnings>50) {
    char printBuffer[512];
    sprintf(printBuffer, "### CoinLpIO::readLp(): %d Variables did not appear in objective function or constraints", warnings);
    handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
						       << CoinMessageEol;
  }

#ifdef LPIO_DEBUG
  printf("CoinLpIO::readLp(): Done with reading the Lp file\n");
#endif

  int *ind = reinterpret_cast< int * >(malloc((maxcoeff + 1) * sizeof(int)));

  for (i = 0; i < cnt_coeff; i++) {
    ind[i] = findHash(colNames[i], 1);

#ifdef LPIO_DEBUG
    printf("CoinLpIO::readLp(): name[%d]: (%s)   ind: %d\n",
      i, colNames[i], ind[i]);
#endif

    if (ind[i] < 0) {
      char str[8192];
      sprintf(str, "### ERROR: Hash table: %s not found\n", colNames[i]);
      throw CoinError(str, "readLp", "CoinLpIO", __FILE__, __LINE__);
    }
  }

  numberColumns_ = numberHash_[1];
  numberElements_ = cnt_coeff - start[0];

  double *obj[MAX_OBJECTIVES];
  // Check for duplicates - first in objectives
  int * whichColumn = new int [numberColumns_];
  char * inRow = new char[numberColumns_];
  memset(inRow,0,numberColumns_);
  int numberDuplicates = 0;

  for (int j = 0; j < num_objectives; j++) {
    obj[j] = reinterpret_cast< double * >(malloc(numberColumns_ * sizeof(double)));
    memset(obj[j], 0, numberColumns_ * sizeof(double));

    for (i = obj_starts[j]; i < obj_starts[j + 1]; i++) {
      icol = findHash(colNames[i], 1);
      if (icol < 0) {
        char str[8192];
        sprintf(str, "### ERROR: Hash table: %s (obj) not found\n", colNames[i]);
        throw CoinError(str, "readLp", "CoinLpIO", __FILE__, __LINE__);
      }
      if (!inRow[icol])
	inRow[icol]=1;
      else
	numberDuplicates++;
      obj[j][icol] += objsense * coeff[i];
    }
    memset(inRow,0,numberColumns_);
  }
  if (numberDuplicates) {
    char str[8192];
    sprintf(str, "### ERROR: %d duplicates in objective\n",
	    numberDuplicates);
    handler_->message(COIN_GENERAL_INFO, messages_) << str
                                                    << CoinMessageEol;
  }

  if (objsense == -1) {
    handler_->message(COIN_GENERAL_INFO, messages_) << " CoinLpIO::readLp(): Maximization problem reformulated as minimization"
                                                    << CoinMessageEol;
    for (int j = 0; j < num_objectives_; j++) {
      objectiveOffset_[j] = -objectiveOffset_[j];
    }
    wasMaximization_ = -1;
  }

  
  for (i = 0; i < cnt_row + 1; i++) {
    start[i] -= cnt_obj;
  }
  // Check for duplicates
  for (int iRow = 0;iRow<numberRows_;iRow++) {
    CoinBigIndex startRow = start[iRow]+cnt_obj;
    CoinBigIndex endRow = start[iRow+1]+cnt_obj;
    for (CoinBigIndex j=startRow;j<endRow;j++) {
      int iColumn = ind[j];
      if (!inRow[iColumn])
	inRow[iColumn]=1;
      else
	numberDuplicates ++;
    }
    for (CoinBigIndex j=startRow;j<endRow;j++) {
      int iColumn = ind[j];
      inRow[iColumn] = 0;
    }
  }
  delete [] inRow;
  delete [] whichColumn;
  if (numberDuplicates) {
    char str[8192];
    sprintf(str, "### ERROR: %d duplicates in objective and matrix\n",
	    numberDuplicates);
    handler_->message(COIN_GENERAL_INFO, messages_) << str
                                                    << CoinMessageEol;
    throw CoinError(str, "readLp", "CoinLpIO", __FILE__, __LINE__);
  }
  CoinPackedMatrix *matrix = new CoinPackedMatrix(false,
    numberColumns_, numberRows_, numberElements_,
    &(coeff[cnt_obj]), &(ind[cnt_obj]), start, NULL);

#ifdef LPIO_DEBUG
  matrix->dumpMatrix();
#endif
  // save sets
  CoinSet **saveSet = set_;
  int saveNumberSets = numberSets_;
  set_ = NULL;
  numberSets_ = 0;
  if (strlen(problemName_))
    handler_->message(COIN_MPS_STATS, messages_) << problemName_
						 << numberRows_
						 << numberColumns_
						 << numberElements_
						 << CoinMessageEol;
  setLpDataWithoutRowAndColNames(*matrix, collow, colup,
    const_cast< const double ** >(obj),
    num_objectives, has_int ? is_int : 0, rowlow, rowup);

  set_ = saveSet;
  numberSets_ = saveNumberSets;

  for (int j = 0; j < num_objectives; j++) {
    if (objName[j] == NULL) {
      rowNames[cnt_row + j] = CoinStrdup("obj");
    } else {
      rowNames[cnt_row + j] = CoinStrdup(objName[j]);
    }
  }
  // Hash tables for column names are already set up
  setLpDataRowAndColNames(rowNames, NULL);

  if (are_invalid_names(names_[1], numberHash_[1], false)) {
    setDefaultColNames();
    handler_->message(COIN_GENERAL_WARNING, messages_) << "### CoinLpIO::readLp(): Invalid column names\nNow using default column names."
                                                       << CoinMessageEol;
  }

  for (i = 0; i < cnt_coeff; i++) {
    free(colNames[i]);
  }
  free(colNames);

  for (i = 0; i < cnt_row + 1; i++) {
    free(rowNames[i]);
  }
  free(rowNames);

  for (int j = 0; j < num_objectives; j++) {
    free(objName[j]);
  }

#ifdef LPIO_DEBUG
  writeLp("readlp.xxx");
  printf("CoinLpIO::readLp(): read Lp file written in file readlp.xxx\n");
#endif

  free(coeff);
  free(start);
  free(ind);
  free(colup);
  free(collow);
  free(rhs);
  free(rowlow);
  free(rowup);
  free(is_int);
  for (int j = 0; j < num_objectives; j++) {
    free(obj[j]);
  }
  delete matrix;

} /* readLp */

/*************************************************************************/
void CoinLpIO::print() const
{

  printf("problemName_: %s\n", problemName_);
  printf("numberRows_: %d\n", numberRows_);
  printf("numberColumns_: %d\n", numberColumns_);

  printf("matrixByRows_:\n");
  matrixByRow_->dumpMatrix();

  int i;
  printf("rowlower_:\n");
  for (i = 0; i < numberRows_; i++) {
    printf("%.5f ", rowlower_[i]);
  }
  printf("\n");

  printf("rowupper_:\n");
  for (i = 0; i < numberRows_; i++) {
    printf("%.5f ", rowupper_[i]);
  }
  printf("\n");

  printf("collower_:\n");
  for (i = 0; i < numberColumns_; i++) {
    printf("%.5f ", collower_[i]);
  }
  printf("\n");

  printf("colupper_:\n");
  for (i = 0; i < numberColumns_; i++) {
    printf("%.5f ", colupper_[i]);
  }
  printf("\n");

  for (int j = 0; j < num_objectives_; j++) {
    printf("objective_[%i]:\n", j);
    for (i = 0; i < numberColumns_; i++) {
      printf("%.5f ", objective_[j][i]);
    }
  }
  printf("\n");

  if (integerType_ != NULL) {
    printf("integerType_:\n");
    for (i = 0; i < numberColumns_; i++) {
      printf("%c ", integerType_[i]);
    }
  } else {
    printf("integerType_: NULL\n");
  }

  printf("\n");
  if (fileName_ != NULL) {
    printf("fileName_: %s\n", fileName_);
  }
  printf("infinity_: %.5f\n", infinity_);
} /* print */

/*************************************************************************/
// Hash functions slightly modified from CoinMpsIO.cpp

namespace {
const int mmult[] = {
  262139, 259459, 256889, 254291, 251701, 249133, 246709, 244247,
  241667, 239179, 236609, 233983, 231289, 228859, 226357, 223829,
  221281, 218849, 216319, 213721, 211093, 208673, 206263, 203773,
  201233, 198637, 196159, 193603, 191161, 188701, 186149, 183761,
  181303, 178873, 176389, 173897, 171469, 169049, 166471, 163871,
  161387, 158941, 156437, 153949, 151531, 149159, 146749, 144299,
  141709, 139369, 136889, 134591, 132169, 129641, 127343, 124853,
  122477, 120163, 117757, 115361, 112979, 110567, 108179, 105727,
  103387, 101021, 98639, 96179, 93911, 91583, 89317, 86939, 84521,
  82183, 79939, 77587, 75307, 72959, 70793, 68447, 66103
};
int compute_hash(const char *name, int maxsiz, int length)
{

  int n = 0;
  int j;

  const int nEntriesMMult = sizeof(mmult) / sizeof(int);

  for (j = 0; j < length; ++j) {
    int iname = name[j];

    n += mmult[j % nEntriesMMult] * iname;
  }
  return (abs(n) % maxsiz); /* integer abs */
}
} // end file-local namespace

/************************************************************************/
//  startHash.  Creates hash list for names
//  setup names_[section] with names in the same order as in the parameter,
//  but removes duplicates

void CoinLpIO::startHash(char const *const *const names,
  const COINColumnIndex number, int section)
{
  maxHash_[section] = 4 * number;
  int maxhash = maxHash_[section];
  COINColumnIndex i, ipos, iput;

  names_[section] = reinterpret_cast< char ** >(malloc(maxhash * sizeof(char *)));
  hash_[section] = new CoinHashLink[maxhash];

  CoinHashLink *hashThis = hash_[section];
  char **hashNames = names_[section];

  for (i = 0; i < maxhash; i++) {
    hashThis[i].index = -1;
    hashThis[i].next = -1;
  }

  /*
   * Initialize the hash table.  Only the index of the first name that
   * hashes to a value is entered in the table; subsequent names that
   * collide with it are not entered.
   */

  for (i = 0; i < number; i++) {
    const char *thisName = names[i];
    int length = CoinStrlenAsInt(thisName);

    ipos = compute_hash(thisName, maxhash, length);
    if (hashThis[ipos].index == -1) {
      hashThis[ipos].index = i; // will be changed below
    }
  }

  /*
   * Now take care of the names that collided in the preceding loop,
   * by finding some other entry in the table for them.
   * Since there are as many entries in the table as there are names,
   * there must be room for them.
   * Also setting up hashNames.
   */

  int cnt_distinct = 0;

  iput = -1;
  for (i = 0; i < number; i++) {
    const char *thisName = names[i];
    int length = CoinStrlenAsInt(thisName);

    ipos = compute_hash(thisName, maxhash, length);

    while (1) {
      COINColumnIndex j1 = hashThis[ipos].index;

      if (j1 == i) {

        // first occurence of thisName in the parameter "names"

        hashThis[ipos].index = cnt_distinct;
        hashNames[cnt_distinct] = CoinStrdup(thisName);
        cnt_distinct++;
        break;
      } else {

#ifdef LPIO_DEBUG
        if (j1 > i) {
          char str[8192];
          sprintf(str, "### ERROR: Hash table: j1: %d  i: %d\n", j1, i);
          throw CoinError(str, "startHash", "CoinLpIO", __FILE__, __LINE__);
        }
#endif

        if (strcmp(thisName, hashNames[j1]) == 0) {

          // thisName already entered

          break;
        } else {
          // Collision; check if thisName already entered

          COINColumnIndex k = hashThis[ipos].next;

          if (k == -1) {

            // thisName not found; enter it

            while (1) {
              ++iput;
              if (iput > maxhash) {
                char str[8192];
                sprintf(str, "### ERROR: Hash table: too many names\n");
                throw CoinError(str, "startHash", "CoinLpIO", __FILE__, __LINE__);
                break;
              }
              if (hashThis[iput].index == -1) {
                break;
              }
            }
            hashThis[ipos].next = iput;
            hashThis[iput].index = cnt_distinct;
            hashNames[cnt_distinct] = CoinStrdup(thisName);
            cnt_distinct++;
            break;
          } else {
            ipos = k;

            // continue the check with names in collision
          }
        }
      }
    }
  }

  numberHash_[section] = cnt_distinct;

} /* startHash */

/**************************************************************************/
//  stopHash.  Deletes hash storage
void CoinLpIO::stopHash(int section)
{
  freePreviousNames(section);
  previous_names_[section] = names_[section];
  card_previous_names_[section] = numberHash_[section];

  delete[] hash_[section];
  hash_[section] = NULL;

  maxHash_[section] = 0;
  numberHash_[section] = 0;

  if (section == 0) {
    for (int j = 0; j < num_objectives_; j++) {
      if (objName_[j] != NULL) {
        free(objName_[j]);
        objName_[j] = NULL;
      }
    }
  }
} /* stopHash */

/**********************************************************************/
//  findHash.  -1 not found
COINColumnIndex
CoinLpIO::findHash(const char *name, int section) const
{
  COINColumnIndex found = -1;

  char **names = names_[section];
  CoinHashLink *hashThis = hash_[section];
  COINColumnIndex maxhash = maxHash_[section];
  COINColumnIndex ipos;

  /* default if we don't find anything */
  if (!maxhash)
    return -1;

  int length = CoinStrlenAsInt(name);

  ipos = compute_hash(name, maxhash, length);
  while (1) {
    COINColumnIndex j1 = hashThis[ipos].index;

    if (j1 >= 0) {
      char *thisName2 = names[j1];

      if (strcmp(name, thisName2) != 0) {
        COINColumnIndex k = hashThis[ipos].next;

        if (k != -1)
          ipos = k;
        else
          break;
      } else {
        found = j1;
        break;
      }
    } else {
      found = -1;
      break;
    }
  }
  return found;
} /* findHash */

/*********************************************************************/
void CoinLpIO::insertHash(const char *thisName, int section)
{

  int number = numberHash_[section];
  int maxhash = maxHash_[section];

  CoinHashLink *hashThis = hash_[section];
  char **hashNames = names_[section];

  int iput = -1;
  int length = CoinStrlenAsInt(thisName);

  int ipos = compute_hash(thisName, maxhash, length);

  while (1) {
    COINColumnIndex j1 = hashThis[ipos].index;

    if (j1 == -1) {
      hashThis[ipos].index = number;
      break;
    } else {
      char *thisName2 = hashNames[j1];

      if (strcmp(thisName, thisName2) != 0) {
        COINColumnIndex k = hashThis[ipos].next;

        if (k == -1) {
          while (1) {
            ++iput;
            if (iput == maxhash) {
              char str[8192];
              sprintf(str, "### ERROR: Hash table: too many names\n");
              throw CoinError(str, "insertHash", "CoinLpIO", __FILE__, __LINE__);
              break;
            }
            if (hashThis[iput].index == -1) {
              break;
            }
          }
          hashThis[ipos].next = iput;
          hashThis[iput].index = number;
          break;
        } else {
          ipos = k;
          /* nothing worked - try it again */
        }
      }
    }
  }

  hashNames[number] = CoinStrdup(thisName);
  (numberHash_[section])++;
}
// Pass in Message handler (not deleted at end)
void CoinLpIO::passInMessageHandler(CoinMessageHandler *handler)
{
  if (defaultHandler_)
    delete handler_;
  defaultHandler_ = false;
  handler_ = handler;
}
// Set language
void CoinLpIO::newLanguage(CoinMessages::Language language)
{
  messages_ = CoinMessage(language);
}

// Get next line into inputBuffer_ (returns number in)
int CoinLpIO::newCardLpIO() const
{
  while (bufferPosition_ == bufferLength_) {
    // check if really new line
    if (!fakeBufferLength_)
      lineNumber_++;
    // new line
    bufferPosition_ = 0;
    bufferLength_ = 0;
    // maximum length of name 100? so this can be 900?
#define BUFFER_LENGTH 900
    char * buffer = inputBuffer_;
    int getLength = 1024;
    if (fakeBufferLength_) {
      // still stuff in fake buffer
      strcpy(inputBuffer_,fakeBuffer_);
      buffer = inputBuffer_+fakeBufferLength_;
      getLength -= fakeBufferLength_;
      fakeBufferLength_ = 0;
    }
    char * ok = input_->gets(buffer, getLength);
    if (!ok)
      return 0;
    int length = strlen(inputBuffer_);
    if (length == 1023) {
      // find first space after 900
      char * space = strchr(inputBuffer_+900,' ');
      assert (space);
      // copy rest
      strcpy(fakeBuffer_,space+1);
      space[1] = '\n';
      space[2] = '\0';
      fakeBufferLength_ = strlen(fakeBuffer_);
      assert (fakeBufferLength_<=130);
    } else {
      // complete line in - make sure \n at end
      int usefulChar;
      for (usefulChar=length-1;usefulChar>=0;usefulChar--) {
	char thisChar = inputBuffer_[usefulChar];
	if (thisChar>=' ') {
	  break;
	}
      }
      // add \n\0
      inputBuffer_[usefulChar+1]='\n';
      inputBuffer_[usefulChar+2]='\0';
      length = usefulChar+3;
    }
    memcpy(originalBuffer_,inputBuffer_,length+1);
    originalBuffer_[length+1] = '\0'; 
    // go to single blanks and remove all blanks before :: or :
    char *colons = strstr(inputBuffer_, "::");
    int nn = 0;
    if (colons) {
      nn = colons - inputBuffer_;
      for (int i = 0; i < nn; i++) {
        if (inputBuffer_[i] != ' ')
          inputBuffer_[bufferLength_++] = inputBuffer_[i];
      }
    }
    while (nn < 1024) {
      if (inputBuffer_[nn] == ':') {
        // take out blank before
        if (inputBuffer_[bufferLength_ - 1] == ' ')
          bufferLength_--;
      }
      if (inputBuffer_[nn] == '\t')
        inputBuffer_[nn] = ' ';
      if (inputBuffer_[nn] == '\0' || inputBuffer_[nn] == '\n' || inputBuffer_[nn] == '\r') {
        if (inputBuffer_[nn] == '\n' || inputBuffer_[nn] == '\r')
          inputBuffer_[nn] = '\n';
        break;
      }
      if (inputBuffer_[nn] != ' ' || inputBuffer_[nn + 1] != ' ')
        inputBuffer_[bufferLength_++] = inputBuffer_[nn];
      nn++;
    }
    //inputBuffer_[bufferLength_]='\n';
    inputBuffer_[bufferLength_] = '\0';
    if (inputBuffer_[0] == ' ')
      bufferPosition_++;
    if (fakeBufferLength_)
      bufferLength_ = -bufferLength_;
  }
  return abs(bufferLength_);
}

// Get next string (returns number in)
int CoinLpIO::fscanfLpIO(char *buff) const
{
  assert(input_);
  if (bufferPosition_ == bufferLength_) {
    int returnCode = newCardLpIO();
    if (!returnCode) {
      if (eofFound_)
        return 0;
      eofFound_ = true;
      warnError("scan_next(): End inserted");
      strcpy(buff, "End");
    }
  }
  char *space = strchr(inputBuffer_ + bufferPosition_, ' ');
  int n = 0;
  int start = 0;
  if (space)
    n = space - (inputBuffer_ + bufferPosition_);
  if (n == 0) {
    if (bufferLength_ >= 0) {
      n = bufferLength_ - bufferPosition_;
    } else {
      // partial line - get more
      start = std::max(abs(bufferLength_) - bufferPosition_, 0);
      memcpy(buff, inputBuffer_ + bufferPosition_, start);
      bufferPosition_ = bufferLength_;
      int returnCode = newCardLpIO();
      if (!returnCode)
        return 0;
      if (inputBuffer_[0] != ' ') {
        space = strchr(inputBuffer_, ' ');
        assert(space || bufferLength_ > 0);
        if (space)
          n = space - (inputBuffer_ + bufferPosition_);
        else
          n = bufferLength_ - bufferPosition_;
      } else {
        n = 0;
      }
    }
  }
  memcpy(buff + start, inputBuffer_ + bufferPosition_, n);
  bufferPosition_ += n;
  if (inputBuffer_[bufferPosition_] == ' ')
    bufferPosition_++;
  buff[start + n] = '\0';
  while (is_comment(buff)) {
    skip_comment(buff);
    int x = fscanfLpIO(buff);
    if (x <= 0) {
      warnError("scan_next(): field expected");
      throw("bad fscanf");
    }
  }
  return n + start;
}
// Throw an error after printing message
void
CoinLpIO::throwError(const char * str,const char * methodName,
		     const char * className, const char * fileName,
		     int line) const
{
  char generalPrint[1200];
  // could tailor message using LPIO_MODIFY_MESSAGES
  sprintf(generalPrint,"Line %d %s",lineNumber_,originalBuffer_);
  handler_->message(COIN_GENERAL_WARNING, messages_) <<
    generalPrint << CoinMessageEol;
  throw CoinError(str,methodName,className,fileName,line);
}
void
CoinLpIO::warnError(const char * printBuffer, int line) const
{
#if LPIO_MODIFY_MESSAGES == 0
  handler_->message(COIN_GENERAL_WARNING, messages_) << printBuffer
						     << CoinMessageEol;
#else
  char generalPrint[1200];
  const char * functions[] = {"is_invalid_name():"};
  const char * copyFrom = printBuffer;
  for (int i=0;i<sizeof(functions)/sizeof(char *);i++) {
    const char * find = strstr(printBuffer,functions[i]);
    if (find) {
      copyFrom = find + strlen(functions[i]);
      break;
    }
  }
  sprintf(generalPrint,"#CoinLpIO %s",copyFrom);
  handler_->message(COIN_GENERAL_WARNING, messages_) <<
    generalPrint << CoinMessageEol;
  if (line>0) {
    sprintf(generalPrint,"Line %d %s",line,originalBuffer_);
    handler_->message(COIN_GENERAL_WARNING, messages_) <<
      generalPrint << CoinMessageEol;
  }
#endif
}
// Get pointer to quadratic objective (or NULL)
CoinPackedMatrix *
CoinLpIO::getQuadraticObjective()
{
  if (quadraticObjective_) {
    return quadraticObjective_;
  } else if (!quadraticList_.size()) {
    return NULL;
  } else {
    // create objective
    int numberElements = quadraticList_.size();
    int * row = new int[2*numberElements];
    int * col = new int[numberElements];
    int * counts = new int[numberColumns_];
    memset(counts,0,numberColumns_*sizeof(int));
    CoinBigIndex * starts = new CoinBigIndex[numberColumns_+1];
    double * element = new double[2*numberElements];
    for (int i=0;i<numberElements;i++) {
      CoinLpQuadratic triple = quadraticList_[i];
      double coefficient = triple.first;
      std::string x = triple.second;
      std::string y = triple.third;
      int ix = columnIndex(x.c_str());
      if (ix<0) {
	char str[8192];
	sprintf(str, "### ERROR: unknown column in quadratic objective %s\n",
		x.c_str());
	throw CoinError(str, "insertHash", "CoinLpIO", __FILE__, __LINE__);
      }
      counts[ix]++;
      int iy = columnIndex(y.c_str());
      if (iy<0) {
	char str[8192];
	sprintf(str, "### ERROR: unknown column in quadratic objective %s\n",
		y.c_str());
	throw CoinError(str, "insertHash", "CoinLpIO", __FILE__, __LINE__);
      }
      col[i] = ix;
      row[i] = iy;
      if (ix!=iy)
	coefficient*=0.5;
      element[i] = coefficient;
    }
    // sort
    starts[0] = 0;
    int n = 0;
    for (int i=0;i<numberColumns_;i++) {
      n += counts[i];
      counts[i] = 0;
      starts[i+1]=n;
    }
    double * element2 = element+numberElements;
    int * row2 = row+numberElements;
    for (int i=0;i<numberElements;i++) {
      int iCol = col[i];
      CoinBigIndex iPut = starts[iCol] + counts[iCol];
      counts[iCol]++;
      row2[iPut] = row[i];
      element2[iPut] = element[i];
    }
    quadraticObjective_ =
      new CoinPackedMatrix(true,numberColumns_,numberColumns_,numberElements,
			   element2,row2,starts,NULL);
    delete [] row;
    delete [] col;
    delete [] element;
    delete [] starts;
    delete [] counts;
    return quadraticObjective_;
  }
}
void CoinLpIO::setQuadraticObjective(CoinPackedMatrix * matrix)
{
  quadraticObjective_ = new CoinPackedMatrix(*matrix);
}
