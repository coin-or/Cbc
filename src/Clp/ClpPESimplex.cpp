// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
/*
   Authors

   Jeremy Omer

   Last update: april 10, 2015

 */

#include <iostream>
#include "ClpPESimplex.hpp"
#include "ClpFactorization.hpp"
#include "ClpMessage.hpp"

/** SHARED METHODS FOR USEFUL ALGEBRAIC OPERATIONS */

/** inner product between a coin vector and a pointer */
double PEdot(CoinIndexedVector &v1, const double *v2)
{
  double sum = 0;
  int size = v1.getNumElements();
  int *indices = v1.getIndices();

  for (int i = 0; i < size; i++)
    sum += v1[indices[i]] * v2[indices[i]];
  return sum;
}

/** inner product between two coin vectors
    call the function with the sparser vector first for efficiency */
double PEdot(CoinIndexedVector &v1, CoinIndexedVector &v2)
{
  double sum = 0;
  int size = v1.getNumElements();
  int *indices = v1.getIndices();

  for (int i = 0; i < size; i++)
    sum += v1[indices[i]] * v2[indices[i]];
  return sum;
}

/** compute the product y + x^T*[A I] for the indices "which" of [A I] and store it into y */
void PEtransposeTimesSubsetAll(ClpSimplex *model, int number, const int *which,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y,
  const double *COIN_RESTRICT rowScale,
  const double *COIN_RESTRICT columnScale)
{
  // get the packed matrix
  CoinPackedMatrix *clpMatrix = model->matrix();

  // get the matrix data pointers
  const int *row = clpMatrix->getIndices();
  const CoinBigIndex *columnStart = clpMatrix->getVectorStarts();
  const int *columnLength = clpMatrix->getVectorLengths();
  const double *elementByColumn = clpMatrix->getElements();

  // there is scaling iff rowScale is not null
  if (rowScale) {
    for (int jColumn = 0; jColumn < number; jColumn++) {
      int iColumn = which[jColumn];
      CoinBigIndex j;
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex next = start + columnLength[iColumn];
      double value = 0.0;

      // if looking at a slack variable, there is no scaling (rowScale=1/columnScale)
      if (iColumn > model->getNumCols()) {
        int jRow = iColumn - model->getNumCols();
        y[iColumn] = x[jRow] * (-1);
      }
      // don't forget the scaling for the decision variables
      else {
        for (j = start; j < next; j++) {
          int jRow = row[j];
          value += x[jRow] * elementByColumn[j] * rowScale[jRow];
        }
        y[iColumn] += value * columnScale[iColumn];
      }
    }
  } else {
    for (int jColumn = 0; jColumn < number; jColumn++) {
      int iColumn = which[jColumn];
      CoinBigIndex j;
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex next = start + columnLength[iColumn];
      double value = 0.0;
      if (iColumn > model->getNumCols()) {
        int jRow = iColumn - model->getNumCols();
        value = x[jRow] * (-1);
      } else {
        for (j = start; j < next; j++) {
          int jRow = row[j];
          value += x[jRow] * elementByColumn[j];
        }
      }
      y[iColumn] += value;
    }
  }
}

/** BASE CLASS FOR THE IMPROVED SIMPLEX
*/

/** Constructor */
ClpPESimplex::ClpPESimplex(ClpSimplex *model)
  : coPrimalDegenerates_(0)
  , primalDegenerates_(NULL)
  , isPrimalDegenerate_(NULL)
  , coDualDegenerates_(0)
  , dualDegenerates_(NULL)
  , isDualDegenerate_(NULL)
  , coCompatibleCols_(0)
  , isCompatibleCol_(NULL)
  , coCompatibleRows_(0)
  , isCompatibleRow_(NULL)
  , model_(model)
  , epsDegeneracy_(1.0e-07)
  , epsCompatibility_(1.0e-07)
  , tempRandom_(NULL)
  , coPrimalDegeneratesAvg_(0)
  , coDualDegeneratesAvg_(0)
  , coCompatibleColsAvg_(0)
  , coCompatibleRowsAvg_(0)
  , coUpdateDegenerates_(0)
  //, coIdentifyCompatibles_(0)
  , coDegeneratePivots_(0)
  , coCompatiblePivots_(0)
  , coDegenerateCompatiblePivots_(0)
  , coDegeneratePivotsConsecutive_(0)
  , coPriorityPivots_(0)
  , doStatistics_(0)
  , lastObjectiveValue_(COIN_DBL_MAX)
  , isLastPivotCompatible_(false)
  , timeCompatibility_(0.0)
  , timeMultRandom_(0.0)
  , timeLinearSystem_(0.0)
  , timeTmp_(0.0)
{

  // the improved simplex object should only be implemented after loading the model
  assert(model_->numberColumns() > 0);

  // size of the original model
  numberColumns_ = model_->numberColumns();
  numberRows_ = model_->numberRows();

  // allocate memory to the arrays
  primalDegenerates_ = reinterpret_cast< int * >(malloc(numberRows_ * sizeof(int)));
  isPrimalDegenerate_ = reinterpret_cast< bool * >(malloc((numberRows_ + numberColumns_) * sizeof(bool)));

  dualDegenerates_ = reinterpret_cast< int * >(malloc(numberColumns_ * sizeof(int)));
  isDualDegenerate_ = reinterpret_cast< bool * >(malloc((numberRows_ + numberColumns_) * sizeof(bool)));

  compatibilityCol_ = reinterpret_cast< double * >(malloc((numberRows_ + numberColumns_) * sizeof(double)));
  isCompatibleCol_ = reinterpret_cast< bool * >(malloc((numberRows_ + numberColumns_) * sizeof(bool)));
  std::fill(isCompatibleCol_, isCompatibleCol_ + numberRows_ + numberColumns_, false);

  compatibilityRow_ = reinterpret_cast< double * >(malloc(numberRows_ * sizeof(double)));
  isCompatibleRow_ = reinterpret_cast< bool * >(malloc(numberRows_ * sizeof(bool)));
  std::fill(isCompatibleRow_, isCompatibleRow_ + numberRows_, false);

  tempRandom_ = reinterpret_cast< double * >(malloc(std::max(numberColumns_, numberRows_) * sizeof(double)));
  // fill
  CoinThreadRandom generator = *model_->randomNumberGenerator();
  for (int i = 0; i < std::max(numberColumns_, numberRows_); i++) {
    double random;
    do
      random = static_cast< int >(generator.randomDouble() * 1.0e6) - 5.0e5;
    while (random == 0.0);
    //random = static_cast<int>(model_->randomNumberGenerator()->randomDouble()*1.0e7)-5.0e6;
    tempRandom_[i] = random;
  }
  if (model_->logLevel() > 2)
    doStatistics_ = model_->logLevel();
}

/** Destructor */
ClpPESimplex::~ClpPESimplex()
{

  // free memory
  if (primalDegenerates_)
    free(primalDegenerates_);
  if (isPrimalDegenerate_)
    free(isPrimalDegenerate_);
  if (dualDegenerates_)
    free(dualDegenerates_);
  if (isDualDegenerate_)
    free(isDualDegenerate_);
  if (isCompatibleCol_)
    free(isCompatibleCol_);
  if (compatibilityCol_)
    free(compatibilityCol_);
  if (isCompatibleRow_)
    free(isCompatibleRow_);
  if (compatibilityRow_)
    free(compatibilityRow_);
  if (tempRandom_)
    free(tempRandom_);

  if (doStatistics_ && model_ && model_->numberIterations()) {
    char generalPrint[200];
    sprintf(generalPrint, "Degenerate pivots   : %d, compatibility time %.2f",
      coDegeneratePivots(),
      timeCompatibility());
    model_->messageHandler()->message(CLP_GENERAL,
      *model_->messagesPointer())
      << generalPrint << CoinMessageEol;

#if 1
    int numberPivots = model_->numberIterations();
    if (coDualDegeneratesAvg()) {
      sprintf(generalPrint, "coDegenAvg/rows %g coCompatAvg/rows %g",
        coDualDegeneratesAvg() / numberRows_,
        coCompatibleRowsAvg() / numberRows_);
      model_->messageHandler()->message(CLP_GENERAL,
        *model_->messagesPointer())
        << generalPrint << CoinMessageEol;
    } else if (coPrimalDegeneratesAvg()) {
      sprintf(generalPrint, "coDegenAvg/columns %g coCompatAvg/columns %g",
        coPrimalDegeneratesAvg() / numberColumns_,
        coCompatibleColsAvg() / numberColumns_);
      model_->messageHandler()->message(CLP_GENERAL,
        *model_->messagesPointer())
        << generalPrint << CoinMessageEol;
    }
    if (numberPivots - coCompatiblePivots()) {
      sprintf(generalPrint, "(coDegeneratePivots()-coDegenerateCompatiblePivots())/( (numberPivots-coCompatiblePivots()) %g", (static_cast< double >(coDegeneratePivots() - coDegenerateCompatiblePivots())) / (static_cast< double >(numberPivots - coCompatiblePivots())));
      model_->messageHandler()->message(CLP_GENERAL,
        *model_->messagesPointer())
        << generalPrint << CoinMessageEol;
    }
    if (coCompatiblePivots()) {
      sprintf(generalPrint, "coDegenerateCompatiblePivots()/coCompatiblePivots() %g", static_cast< double >(coDegenerateCompatiblePivots()) / static_cast< double >(coCompatiblePivots()));
      model_->messageHandler()->message(CLP_GENERAL,
        *model_->messagesPointer())
        << generalPrint << CoinMessageEol;
    }
    sprintf(generalPrint, "coDegeneratePivots()/ numberPivots %g", static_cast< double >(coDegeneratePivots()) / static_cast< double >(numberPivots));
    model_->messageHandler()->message(CLP_GENERAL,
      *model_->messagesPointer())
      << generalPrint << CoinMessageEol;
    sprintf(generalPrint, "coCompatiblePivots() %d coPriorityPivots() %d",
      coCompatiblePivots(), coPriorityPivots());
    model_->messageHandler()->message(CLP_GENERAL,
      *model_->messagesPointer())
      << generalPrint << CoinMessageEol;
    //sprintf(generalPrint, numberPivots);
    // sprintf(generalPrint, timeMultRandom());
    //sprintf(generalPrint, timeLinearSystem());
#endif
  }
}

/** Updates the set of variables that are not at their bounds */
void ClpPESimplex::updatePrimalDegenerates()
{

  coPrimalDegenerates_ = 0;
  epsDegeneracy_ = 1.0e-04; //std::min(1.0e-2, 1.0e-7*fabs(model_->objectiveValue()));
  std::fill(isPrimalDegenerate_, isPrimalDegenerate_ + numberRows_ + numberColumns_, false);
  int *pivotVariable = model_->pivotVariable();

  // Only go over the basic variables, since non basic varaibles are
  // set to one of their bounds (or zero if the variable is free)
  // An epsDegeneracy_ tolerance is used to detect variables at a bound
  for (int i = 0; i < numberRows_; i++) {
    int iVar = pivotVariable[i]; // index of basic variable at row i

    // std::cout << "solution[" << iVar << "] = " << model_->solution(iVar) << " ; lb = " << model_->lower(iVar) << " ; ub = "<< model_->upper(iVar) << "\n" ;

    double dVal = model_->solution(iVar);
    double dUb = model_->upper(iVar);
    double dLb = model_->lower(iVar);
    // if (fabs(dVal) <= epsDegeneracy_) {
    if ((dLb > -COIN_DBL_MAX && fabs(dVal - dLb) <= std::max(1.0, fabs(dLb)) * epsDegeneracy_)
      || (dUb < COIN_DBL_MAX && fabs(dVal - dUb) <= std::max(1.0, fabs(dUb)) * epsDegeneracy_)) {
      primalDegenerates_[coPrimalDegenerates_++] = i;
      isPrimalDegenerate_[iVar] = true;
    }
  }
  coUpdateDegenerates_++;
}

/** Updates the set of dual degenerate variables */
void ClpPESimplex::updateDualDegenerates()
{

  coDualDegenerates_ = 0;
  std::fill(isDualDegenerate_, isDualDegenerate_ + numberRows_ + numberColumns_, false);

  // The dual degenerate variables are the nonbasic variables with a zero reduced costs
  // An epsDegeneracy_ tolerance is used to detect zero reduced costs
  epsDegeneracy_ = 1.0e-04; //std::min(1.0e-2, 1.0e-7*fabs(model_->objectiveValue()));
  double maxDegen = 0.0;
  for (int i = 0; i < numberColumns_ + numberRows_; i++) {

    if (model_->getStatus(i) != ClpSimplex::basic && fabs(model_->reducedCost(i)) <= epsDegeneracy_) {
      dualDegenerates_[coDualDegenerates_++] = i;
      isDualDegenerate_[i] = true;
      maxDegen = std::max(maxDegen, fabs(model_->reducedCost(i)));
    }
  }
  coUpdateDegenerates_++;
}

/** Identify the primal compatible columns in a subset of columns
  The input argument is a temporary array that is needed for the Clp's BTRAN
  if which is NULL, every variable is tested */
void ClpPESimplex::identifyCompatibleCols(int number, const int *which,
  CoinIndexedVector *spareRow2,
  CoinIndexedVector *wPrimal)
{

  // the update of variables within their bounds must have been called at least once
  assert(primalDegenerates_);

  // initialize every variable to incompatible
  coCompatibleCols_ = 0;
  std::fill(isCompatibleCol_, isCompatibleCol_ + numberRows_ + numberColumns_, false);
  std::fill(compatibilityCol_, compatibilityCol_ + numberRows_ + numberColumns_, -1.0);

  // treat the two obvious cases:
  // - no primal degenerate variable => every column is compatible
  // - only primal degenerate variables => no compatible column
  if (coPrimalDegenerates_ == 0) {
    if (!which) {
      std::fill(isCompatibleCol_, isCompatibleCol_ + numberRows_ + numberColumns_, true);
      coCompatibleCols_ = numberRows_ + numberColumns_;
    } else {
      for (int j = 0; j < number; j++) {
        isCompatibleCol_[which[j]] = true;
      }
      coCompatibleCols_ = number;
    }
    return;
  } else if (coPrimalDegenerates_ == numberRows_) {
    return;
  }

  // fill the elements of w corresponding to the degenerate variables
  // with random values and leave the other to zero
  // no longer using wPrimal_
  //wPrimal_->clear();
  wPrimal->checkClear();
  assert(coPrimalDegenerates_ <= std::max(numberColumns_, numberRows_));
  for (int i = 0; i < coPrimalDegenerates_; i++) {
#if 0
    double random;
      do
        random = static_cast<double> ((rand() %(static_cast<int> (1.0e6)+1))) -5.0e5;
      while (random == 0.0);
    wPrimal->quickInsert(primalDegenerates_[i], 0.5+random);
#else
    wPrimal->quickInsert(primalDegenerates_[i], tempRandom_[i]);
#endif
  }

  // compute wTran * Binv and store it into wTran
#ifdef PE_TEST
  model_->factorization()->doStatistics(false);
#endif
  model_->factorization()->updateColumnTranspose(spareRow2, wPrimal);
#ifdef PE_TEST
  model_->factorization()->doStatistics(true);
#endif

  // compute wTran * AN, where AN is the matrix of nonbasic variables
  // the zero elements (with a tolerance epsCompatibility_) are the compatibles
  coCompatibleCols_ = 0;
  number = which ? number : numberColumns_ + numberRows_;
  int jColumn;
  assert(!wPrimal->packedMode());
  double *values = wPrimal->denseVector();
  const double *rowScale = model_->rowScale();
  CoinPackedMatrix *clpMatrix = model_->matrix();
  const int *row = clpMatrix->getIndices();
  const CoinBigIndex *columnStart = clpMatrix->getVectorStarts();
  const int *columnLength = clpMatrix->getVectorLengths();
  const double *elementByColumn = clpMatrix->getElements();
  for (int j = 0; j < number; j++) {
    if (which)
      jColumn = which[j];
    else
      jColumn = j;

    if (model_->getStatus(jColumn) == ClpSimplex::basic
      /* && !isPrimalDegenerate_[jColumn]*/) {
      isCompatibleCol_[jColumn] = false;
      /*coCompatibleCols_++;*/
    } else {
      double dotProduct = 0.0;
      if (jColumn < numberColumns_) {
        if (!rowScale) {
          for (CoinBigIndex i = columnStart[jColumn]; i < columnStart[jColumn] + columnLength[jColumn]; i++) {
            int iRow = row[i];
            dotProduct += values[iRow] * elementByColumn[i];
          }
        } else {
          // apply scaling
          double scale = model_->columnScale()[jColumn];
          for (CoinBigIndex i = columnStart[jColumn]; i < columnStart[jColumn] + columnLength[jColumn]; i++) {
            int iRow = row[i];
            dotProduct += values[iRow] * elementByColumn[i] * rowScale[iRow];
          }
          dotProduct *= scale;
        }
      } else {
        // slack
        dotProduct = values[jColumn - numberColumns_];
      }
      dotProduct = fabs(dotProduct);
#if 0
      model_->unpack(tempColumn_, jColumn);
      // perform the inner product <tempColumn_,wPrimal_>
      double dotProduct2 = fabs( PEdot(*tempColumn_, *wPrimal) );
      assert (fabs(dotProduct-dotProduct2)<1.0e-6);
#endif
      compatibilityCol_[jColumn] = dotProduct;

      if (dotProduct < epsCompatibility_) {
        isCompatibleCol_[jColumn] = true;
        coCompatibleCols_++;
      }
    }
  }
  wPrimal->clear();
}

/** Identify the dual compatible rows */
void ClpPESimplex::identifyCompatibleRows(CoinIndexedVector *spare,
  CoinIndexedVector *wDual)
{

  // the update of dual degenerate variables must have been called at least once
  assert(dualDegenerates_);

  // only one case is obvious here:
  // - no dual degenerate variable => every row is compatible
  if (coDualDegenerates_ == 0) {
    std::fill(isCompatibleRow_, isCompatibleRow_ + numberRows_, false);
    coCompatibleRows_ = numberRows_;
    return;
  }
  assert(coDualDegenerates_ <= std::max(numberColumns_, numberRows_));
#if 0
  // fill the elements of tempRandom with as many random elements as dual degenerate variables
  for (int j = 0; j < coDualDegenerates_; j++) {
    //do
    tempRandom_[j] = static_cast<double> ((rand() %(static_cast<int> (1.0e7)+1))) -5.0e6;
    //while (tempRandom_[j] == 0.0);
  }
#endif

  // compute the product A_D * tempRandom and store it into wDual
  // No longer using wDual_
  wDual->checkClear();
  //wDual_->clear();
  double timeTmp = 0.0;
  if (doStatistics_)
    timeTmp = CoinCpuTime();
  double *values = wDual->denseVector();
  const double *rowScale = model_->rowScale();
  CoinPackedMatrix *clpMatrix = model_->matrix();
  const int *row = clpMatrix->getIndices();
  const CoinBigIndex *columnStart = clpMatrix->getVectorStarts();
  const int *columnLength = clpMatrix->getVectorLengths();
  const double *elementByColumn = clpMatrix->getElements();
  for (int j = 0; j < coDualDegenerates_; j++) {
    // I try and save time to avoid calling unpack
    int sequence = dualDegenerates_[j];
    if (sequence >= numberColumns_) {
      //slack
      values[sequence - numberColumns_] -= tempRandom_[j];
    } else {
      // column
      CoinBigIndex i;
      if (!rowScale) {
        for (i = columnStart[sequence]; i < columnStart[sequence] + columnLength[sequence]; i++) {
          int iRow = row[i];
          values[iRow] += tempRandom_[j] * elementByColumn[i];
        }
      } else {
        // apply scaling
        double scale = model_->columnScale()[sequence];
        for (i = columnStart[sequence]; i < columnStart[sequence] + columnLength[sequence]; i++) {
          int iRow = row[i];
          values[iRow] += tempRandom_[j] * elementByColumn[i] * scale * rowScale[iRow];
        }
      }
    }
    //       matrix_->unpack(this, rowArray, sequenceIn_);
    //  }
    // model_->unpack(tempColumn_, dualDegenerates_[j]);

    // int nz = tempColumn_->getNumElements();
    // int *ind = tempColumn_->getIndices();
    // double *val = tempColumn_->denseVector();

    // for (int k = 0; k < nz ; k++) {
    //   values[ind[k]] = values[ind[k]] + tempRandom_[j]*val[ind[k]];
    // }
  }
  int n = 0;
  int *indices = wDual->getIndices();
  for (int i = 0; i < numberRows_; i++) {
    if (values[i])
      indices[n++] = i;
  }
  wDual->setNumElements(n);
  wDual->setPackedMode(false);

  // compute Binv*wDual and store it into wDual
#ifdef PE_TEST
  if (doStatistics_) {
    timeMultRandom_ += CoinCpuTime() - timeTmp;
    timeTmp = CoinCpuTime();
  }
  model_->factorization()->doStatistics(false);
#endif
  model_->factorization()->updateColumn(spare, wDual);
#ifdef PE_TEST
  model_->factorization()->doStatistics(true);
  if (doStatistics_) {
    timeLinearSystem_ += CoinCpuTime() - timeTmp;
  }
#endif

  // the zero elements of wDual (with a tolerance epsCompatibility_) are the
  // dual-compatible variables
  assert(!wDual->packedMode());
  int size = wDual->getNumElements();
  std::fill(isCompatibleRow_, isCompatibleRow_ + numberRows_, true);
  coCompatibleRows_ = numberRows_;
  for (int i = 0; i < size; i++) {
    int iRow = indices[i];
    double value = values[iRow];
    if (fabs(value) >= epsCompatibility_ * 100.0) {
      isCompatibleRow_[iRow] = false;
      coCompatibleRows_--;
    }
  }
  wDual->clear();
}

/* Update the average number of primal degenerate columns
	The input is the number of pivots since last update.

*/
void ClpPESimplex::updatePrimalDegeneratesAvg(int coPivots)
{
  int totalPivots = model_->numberIterations() + 1;
  double fracPivots = static_cast< double >(coPivots) / totalPivots;
  coPrimalDegeneratesAvg_ = floor((1.0 - fracPivots) * (coPrimalDegeneratesAvg_ + fracPivots * coPrimalDegenerates_));
}

/* Update the average number of dual degenerate columns
	The input is the number of pivots since last update.

*/
void ClpPESimplex::updateDualDegeneratesAvg(int coPivots)
{
  int totalPivots = model_->numberIterations() + 1;
  double fracPivots = static_cast< double >(coPivots) / totalPivots;
  coDualDegeneratesAvg_ = floor((1.0 - fracPivots) * coDualDegeneratesAvg_ + fracPivots * coDualDegenerates_);
}

/* Update the average number of compatible columns.
	The input is the number of pivots since last update.
*/
void ClpPESimplex::updateCompatibleColsAvg(int coPivots)
{
  int totalPivots = model_->numberIterations() + 1;
  double fracPivots = static_cast< double >(coPivots) / totalPivots;
  coCompatibleColsAvg_ = floor((1.0 - fracPivots) * coCompatibleColsAvg_ + fracPivots * (coCompatibleCols_ /*-(numberRows_-coPrimalDegenerates_)*/));
}

/* Update the average number of compatible rows.
	The input is the number of pivots since last update.
*/
void ClpPESimplex::updateCompatibleRowsAvg(int coPivots)
{
  int totalPivots = model_->numberIterations() + 1;
  double fracPivots = static_cast< double >(coPivots) / totalPivots;
  coCompatibleRowsAvg_ = floor((1.0 - fracPivots) * coCompatibleRowsAvg_ + fracPivots * coCompatibleRows_);
}

/** DEBUG METHODS */

#if PE_DEBUG >= 1
/** Print the set of variables within their bounds */
void ClpPESimplex::printPrimalDegenerates()
{
  assert(primalDegenerates_);

  if (!coPrimalDegenerates_) {
    std::cout << "There is no primal degenerate variable!" << std::endl;
    return;
  }

  int *pivotVariable = model_->pivotVariable();

  std::cout << "List of primal degenerate variables:";
  for (int i = 0; i < coPrimalDegenerates_; i++) {
    std::cout << " " << pivotVariable[primalDegenerates_[i]] << " ;";
  }
  std::cout << std::endl;
}

/** Print the set of primal compatible variables */
void ClpPESimplex::printCompatibleCols()
{
  assert(isCompatibleCol_);

  if (!coCompatibleCols_) {
    std::cout << "There is no compatible column!" << std::endl;
    return;
  }

  std::cout << "List of compatible columns:";
  for (int i = 0; i < numberColumns_ + numberRows_; i++) {
    if (isCompatibleCol_[i])
      std::cout << " " << i << " ;";
  }
  std::cout << std::endl;
}

/** Check that a nonbasic variable is indeed contemptible */
bool ClpPESimplex::checkCompatibilityCol(int sequence,
  CoinIndexedVector *spareRow2)
{

  bool isCompatible = true;
  tempColumn_ = new CoinIndexedVector();
  tempColumn_->reserve(numberRows_);

  model_->unpack(tempColumn_, sequence);
#ifdef PE_TEST
  model_->factorization()->doStatistics(false);
#endif
  model_->factorization()->updateColumn(spareRow2, tempColumn_);
#ifdef PE_TEST
  model_->factorization()->doStatistics(true);
#endif

  for (int j = 0; j < coPrimalDegenerates_; j++) {
    if (fabs((*tempColumn_)[primalDegenerates_[j]]) > epsDegeneracy_) {
      return false;
    }
  }
  delete tempColumn_;
  return isCompatible;
}
#endif
/** Check that a basic row is indeed compatible */
bool ClpPESimplex::checkCompatibilityRow(int pivotRow)
{

  bool isCompatible = true;
  double direction = 1;
  model_->rowArray(0)->createPacked(1, &pivotRow, &direction);
#ifdef PE_TEST
  model_->factorization()->doStatistics(false);
#endif
  model_->factorization()->updateColumnTranspose(model_->rowArray(1), model_->rowArray(0));
#ifdef PE_TEST
  model_->factorization()->doStatistics(true);
#endif
  model_->clpMatrix()->transposeTimes(model_, -1.0, model_->rowArray(0), model_->rowArray(1), model_->columnArray(0));

  CoinIndexedVector *columnArray = model_->columnArray(0);
  CoinIndexedVector *rowArray = model_->rowArray(0);
  int nzCol = columnArray->getNumElements();
  int *indCol = columnArray->getIndices();
  double *valCol = columnArray->denseVector();
  int nzRow = rowArray->getNumElements();
  int *indRow = rowArray->getIndices();
  double *valRow = rowArray->denseVector();

  if (columnArray->packedMode()) {
    for (int j = 0; j < nzCol; j++) {
      int iCol = indCol[j];
      if (isDualDegenerate_[iCol] && fabs(valCol[j]) > epsDegeneracy_) {
        std::cout << "Dual degenerate column: " << valCol[j] << std::endl;
      }
    }
  } else {
    for (int j = 0; j < nzCol; j++) {
      int iCol = indCol[j];
      if (isDualDegenerate_[iCol] && fabs(valCol[iCol]) > epsDegeneracy_) {
        std::cout << "Dual degenerate column: " << valCol[iCol] << std::endl;
      }
    }
  }
  if (rowArray->packedMode()) {
    for (int j = 0; j < nzRow; j++) {
      int iRow = indRow[j];
      if (isDualDegenerate_[iRow + numberColumns_] && fabs(valRow[j]) > epsDegeneracy_) {
        std::cout << "Dual degenerate row: " << valRow[j] << std::endl;
      }
    }
  } else {
    for (int j = 0; j < nzRow; j++) {
      int iRow = indRow[j];
      if (isDualDegenerate_[iRow + numberColumns_] && fabs(valRow[iRow]) > epsDegeneracy_) {
        std::cout << "Dual degenerate row: " << valRow[iRow] << std::endl;
      }
    }
  }

  return isCompatible;
}
// checks size
bool ClpPESimplex::checkSize()
{
  return (numberRows_ == model_->numberRows() && numberColumns_ == model_->numberColumns());
}
/* Update the dual compatible rows */
void ClpPESimplex::updateCompatibleRows(int iColumn)
{
  if (iColumn < numberColumns_) {
    // get the packed matrix
    CoinPackedMatrix *clpMatrix = model_->matrix();

    // get the matrix data pointers
    const int *row = clpMatrix->getIndices();
    const CoinBigIndex *columnStart = clpMatrix->getVectorStarts();
    const int *columnLength = clpMatrix->getVectorLengths();
    CoinBigIndex start = columnStart[iColumn];
    CoinBigIndex next = start + columnLength[iColumn];
    for (CoinBigIndex j = start; j < next; j++) {
      int jRow = row[j];
      if (isCompatibleRow_[jRow]) {
        isCompatibleRow_[jRow] = false;
        coCompatibleRows_--;
      }
    }
  } else {
    int jRow = iColumn - numberColumns_;
    if (isCompatibleRow_[jRow]) {
      isCompatibleRow_[jRow] = false;
      coCompatibleRows_--;
    }
  }
}

/** Print the total recorded time */
void ClpPESimplex::printTimer(std::ostream &out)
{
  out << "Cpu in compatibility: " << timeCompatibility_ << " s" << std::endl;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
