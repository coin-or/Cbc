#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <OsiSolverInterface.hpp>
#include "CbcMessage.hpp"
#include "CbcHeuristic.hpp"
#include <CbcModel.hpp>
#include "CbcMipStartIO.hpp"
#include "CbcSOS.hpp"
#include "CoinTime.hpp"

using namespace std;

static
bool isNumericStr(const char *str)
{
  const size_t l = strlen(str);

  for (size_t i = 0; i < l; ++i)
    if (!(isdigit(str[i]) || (str[i] == '.') || (str[i] == '-') || (str[i] == '+') || (str[i] == 'e')))
      return false;

  return true;
}

int readMIPStart(CbcModel *model, const char *fileName,
  vector< pair< string, double > > &colValues,
  double & /*solObj*/)
{
#define STR_SIZE 256
  FILE *f = fopen(fileName, "r");
  if (!f)
    return 1;
  char line[STR_SIZE];

  int nLine = 0;
  char printLine[STR_SIZE];
  while (fgets(line, STR_SIZE, f)) {
    ++nLine;
    char col[4][STR_SIZE];
    int nread = sscanf(line, "%s %s %s %s", col[0], col[1], col[2], col[3]);
    if (!nread)
      continue;
    /* line with variable value */
    if (strlen(col[0]) && isdigit(col[0][0]) && (nread >= 3)) {
      if (!isNumericStr(col[0])) {
        sprintf(printLine, "Reading: %s, line %d - first column in mipstart file should be numeric, ignoring.", fileName, nLine);
        model->messageHandler()->message(CBC_GENERAL, model->messages()) << printLine << CoinMessageEol;
        continue;
      }
      if (!isNumericStr(col[2])) {
        sprintf(printLine, "Reading: %s, line %d - Third column in mipstart file should be numeric, ignoring.", fileName, nLine);
        model->messageHandler()->message(CBC_GENERAL, model->messages()) << printLine << CoinMessageEol;
        continue;
      }

      char *name = col[1];
      double value = atof(col[2]);

      colValues.push_back(pair< string, double >(string(name), value));
    }
  }

  if (colValues.size()) {
    sprintf(printLine, "MIPStart values read for %d variables.", static_cast< int >(colValues.size()));
    model->messageHandler()->message(CBC_GENERAL, model->messages()) << printLine << CoinMessageEol;
    if ((int)colValues.size() < model->getNumCols()) {
      int numberColumns = model->getNumCols();
      OsiSolverInterface *solver = model->solver();
      vector< pair< string, double > > fullValues;
      /* for fast search of column names */
      map< string, int > colIdx;
      for (int i = 0; i < numberColumns; i++) {
        fullValues.push_back(pair< string, double >(solver->getColName(i), 0.0));
        colIdx[solver->getColName(i)] = i;
      }
      for (int i = 0; (i < static_cast< int >(colValues.size())); ++i) {
        map< string, int >::const_iterator mIt = colIdx.find(colValues[i].first);
        if (mIt != colIdx.end()) {
          const int idx = mIt->second;
          double v = colValues[i].second;
          fullValues[idx].second = v;
        }
      }
      colValues = fullValues;
    }
  } else {
    sprintf(printLine, "No mipstart solution read from %s", fileName);
    model->messageHandler()->message(CBC_GENERAL, model->messages()) << printLine << CoinMessageEol;
    fclose(f);
    return 1;
  }

  fclose(f);
  return 0;
}

int computeCompleteSolution(CbcModel *model,
  const vector< string > colNames,
  const std::vector< std::pair< std::string, double > > &colValues,
  double *sol, double &obj)
{
  if (!model->getNumCols())
    return 0;

  int status = 0;
  double compObj = COIN_DBL_MAX;
  bool foundIntegerSol = false;
  OsiSolverInterface *lp = model->solver()->clone();
  map< string, int > colIdx;
  assert((static_cast< int >(colNames.size())) == lp->getNumCols());
  /* for fast search of column names */
  for (int i = 0; (i < static_cast< int >(colNames.size())); ++i)
    colIdx[colNames[i]] = i;

  char printLine[STR_SIZE+100];
  int fixed = 0;
  int notFound = 0;
  char colNotFound[256] = "";
  int nContinuousFixed = 0;
  double *realObj = new double[lp->getNumCols()];
  memcpy(realObj, lp->getObjCoefficients(), sizeof(double)*lp->getNumCols());

#ifndef JUST_FIX_INTEGER
#define JUST_FIX_INTEGER 0
#endif

#if JUST_FIX_INTEGER > 1
  // all not mentioned are at zero
  for (int i = 0; (i < lp->getNumCols()); ++i) {
    if (lp->isInteger(i))
      lp->setColBounds(i, 0.0, 0.0);
  }
#endif
  const double *lower = lp->getColLower();
  const double *upper = lp->getColUpper();
  int nBadValues = 0;
  for (int i = 0; (i < static_cast< int >(colValues.size())); ++i) {
    map< string, int >::const_iterator mIt = colIdx.find(colValues[i].first);
    if (mIt == colIdx.end()) {
      if (!notFound)
        strcpy(colNotFound, colValues[i].first.c_str());
      notFound++;
    } else {
      const int idx = mIt->second;
      double v = colValues[i].second;
      if (v>upper[idx]) {
	nBadValues++;
	v = upper[idx];
      } else if (v<lower[idx]) {
	nBadValues++;
	v = lower[idx];
      }
#if JUST_FIX_INTEGER
      if (!lp->isInteger(idx))
        continue;
#endif
      if (fabs(v) < 1e-8)
        v = 0.0;
      if (lp->isInteger(idx)) // just to avoid small
        v = floor(v + 0.5); // fractional garbage
      else
        nContinuousFixed++;

      lp->setColBounds(idx, v, v);
      ++fixed;
    }
  }
  if (nBadValues) {
    char line[100];
    sprintf(line,"Warning: modifying %d solution values outside bounds",
	    nBadValues);
    model->messageHandler()->message(CBC_GENERAL, model->messages())
      << line << CoinMessageEol;
  }

  if (!fixed) {
    model->messageHandler()->message(CBC_GENERAL, model->messages())
      << "Warning: MIPstart solution is not valid, column names do not match, ignoring it."
      << CoinMessageEol;
    goto TERMINATE;
  }

  if (notFound >= ((static_cast< double >(colNames.size())) * 0.5)) {
    sprintf(printLine, "Warning: %d column names were not found (e.g. %s) while filling solution.", notFound, colNotFound);
    model->messageHandler()->message(CBC_GENERAL, model->messages())
      << printLine << CoinMessageEol;
  }
#if JUST_FIX_INTEGER
  lp->setHintParam(OsiDoPresolveInInitial, true, OsiHintDo);
#endif
  lp->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
  lp->initialSolve();

  if ((lp->isProvenPrimalInfeasible()) || (lp->isProvenDualInfeasible())) {
    if (nContinuousFixed) {
      model->messageHandler()->message(CBC_GENERAL, model->messages())
        << "Trying just fixing integer variables (and fixingish SOS)." << CoinMessageEol;
      int numberColumns = lp->getNumCols();
      const double *oldLower = model->solver()->getColLower();
      const double *oldUpper = model->solver()->getColUpper();
      double *savedSol = CoinCopyOfArray(lp->getColLower(), numberColumns);
      for (int i = 0; i < numberColumns; ++i) {
        if (!lp->isInteger(i)) {
          lp->setColLower(i, oldLower[i]);
          lp->setColUpper(i, oldUpper[i]);
        }
      }
      // but look at SOS
      int numberObjects = model->numberObjects();
      for (int i = 0; i < numberObjects; i++) {
        const CbcSOS *object = dynamic_cast< const CbcSOS * >(model->object(i));
        if (object) {
          int n = object->numberMembers();
          const int *members = object->members();
          int sosType = object->sosType();
          if (sosType == 1) {
            // non zero can take any value - others zero
            int iColumn = -1;
            for (int j = 0; j < n; j++) {
              int jColumn = members[j];
              if (savedSol[jColumn])
                iColumn = jColumn;
            }
            for (int j = 0; j < n; j++) {
              int jColumn = members[j];
              if (jColumn != iColumn) {
                lp->setColLower(jColumn, 0.0);
                lp->setColUpper(jColumn, 0.0);
              }
            }
          } else if (sosType == 2) {
            // SOS 2 - make a guess if just one nonzero
            int jA = -1;
            int jB = -1;
            for (int j = 0; j < n; j++) {
              int jColumn = members[j];
              if (savedSol[jColumn]) {
                if (jA == -1)
                  jA = j;
                jB = j;
              }
            }
            if (jB > jA + 1) {
              jB = jA + 1;
            } else if (jA == jB) {
              if (jA == n - 1)
                jA--;
              else
                jB++;
            }
            for (int j = 0; j < n; j++) {
              if (j != jA && j != jB) {
                int jColumn = members[j];
                lp->setColLower(jColumn, 0.0);
                lp->setColUpper(jColumn, 0.0);
              }
            }
          }
        }
      }
      delete[] savedSol;
      lp->initialSolve();
    } else {
      model->messageHandler()->message(CBC_GENERAL, model->messages())
        << "Fixing only non-zero variables." << CoinMessageEol;
      /* unfix all variables which are zero */
      int notZeroAnymore = 0;
      for (int i = 0; (i < lp->getNumCols()); ++i)
        if (((fabs(lp->getColLower()[i])) <= 1e-8) && (fabs(lp->getColLower()[i] - lp->getColUpper()[i]) <= 1e-8)) {
          const double *oldLower = model->solver()->getColLower();
          const double *oldUpper = model->solver()->getColUpper();
          lp->setColLower(i, oldLower[i]);
          lp->setColUpper(i, oldUpper[i]);
          notZeroAnymore++;
        }
      if (notZeroAnymore)
        lp->initialSolve();
    }
  }

  if (!lp->isProvenOptimal()) {
    model->messageHandler()->message(CBC_GENERAL, model->messages())
      << "Warning: mipstart values could not be used to build a solution." << CoinMessageEol;
    status = 1;
    goto TERMINATE;
  }

  /* some additional effort is needed to provide an integer solution */
  if (lp->getFractionalIndices().size() > 0) {
    sprintf(printLine, "MIPStart solution provided values for %d of %d integer variables, %d variables are still fractional.", fixed, lp->getNumIntegers(), static_cast< int >(lp->getFractionalIndices().size()));
    model->messageHandler()->message(CBC_GENERAL, model->messages())
      << printLine << CoinMessageEol;
    double start = CoinCpuTime();
#if 1
    CbcSerendipity heuristic(*model);
    heuristic.setFractionSmall(2.0);
    heuristic.setFeasibilityPumpOptions(1008013);
    int returnCode = heuristic.smallBranchAndBound(lp,
      1000, sol,
      compObj,
      model->getCutoff(),
      "ReduceInMIPStart");
    if ((returnCode & 1) != 0) {
      sprintf(printLine, "Mini branch and bound defined values for remaining variables in %.2f seconds.",
        CoinCpuTime() - start);
      model->messageHandler()->message(CBC_GENERAL, model->messages())
        << printLine << CoinMessageEol;
      foundIntegerSol = true;
      obj = compObj;
    }
#else
    CbcModel babModel(*lp);
    lp->writeLp("lessFix");
    babModel.setLogLevel(2);
    babModel.setMaximumNodes(1000);
    babModel.setMaximumSeconds(60);
    babModel.branchAndBound();
    if (babModel.bestSolution()) {
      sprintf(printLine, "Mini branch and bound defined values for remaining variables in %.2f seconds.",
        CoinCpuTime() - start);
      model->messageHandler()->message(CBC_GENERAL, model->messages())
        << printLine << CoinMessageEol;
      copy(babModel.bestSolution(), babModel.bestSolution() + babModel.getNumCols(), sol);
      foundIntegerSol = true;
      obj = compObj = babModel.getObjValue();
    }
#endif
    else {
      model->messageHandler()->message(CBC_GENERAL, model->messages())
        << "Warning: mipstart values could not be used to build a solution." << CoinMessageEol;
      status = 1;
      goto TERMINATE;
    }
  } else {
    foundIntegerSol = true;
    
    obj = 0.0;
    for ( int i=0 ; (i<lp->getNumCols()) ; ++i )
        obj += realObj[i]*lp->getColSolution()[i];
    compObj = obj;
    copy(lp->getColSolution(), lp->getColSolution() + lp->getNumCols(), sol);
  }

  if (foundIntegerSol) {
    sprintf(printLine, "MIPStart provided solution with cost %g", compObj);
    model->messageHandler()->message(CBC_GENERAL, model->messages())
      << printLine << CoinMessageEol;
#if 0
      {
	int numberColumns=lp->getNumCols();
	double largestInfeasibility = 0.0;
	double primalTolerance ;
	double offset;
	lp->getDblParam(OsiObjOffset, offset);
	lp->getDblParam(OsiPrimalTolerance, primalTolerance) ;
	const double *objective = lp->getObjCoefficients() ;
	const double * rowLower = lp->getRowLower() ;
	const double * rowUpper = lp->getRowUpper() ;
	const double * columnLower = lp->getColLower() ;
	const double * columnUpper = lp->getColUpper() ;
	int numberRows = lp->getNumRows() ;
	double *rowActivity = new double[numberRows] ;
	memset(rowActivity, 0, numberRows*sizeof(double)) ;
	double *rowSum = new double[numberRows] ;
	memset(rowSum, 0, numberRows*sizeof(double)) ;
	const double * element = lp->getMatrixByCol()->getElements();
	const int * row = lp->getMatrixByCol()->getIndices();
	const CoinBigIndex * columnStart = lp->getMatrixByCol()->getVectorStarts();
	const int * columnLength = lp->getMatrixByCol()->getVectorLengths();
	const CoinPackedMatrix * rowCopy = lp->getMatrixByRow();
	const int * column = rowCopy->getIndices();
	const int * rowLength = rowCopy->getVectorLengths();
	const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
	const double * elementByRow = rowCopy->getElements();
	double objValue=-offset;
	for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	  double value = sol[iColumn];
	  if (lp->isInteger(iColumn))
	    assert (fabs(value-floor(value+0.5))<1.0e-6);
	  objValue += value*objective[iColumn];
	  if (value>columnUpper[iColumn]) {
	    if (value-columnUpper[iColumn]>1.0e-8)
	      printf("column %d has value %.12g above %.12g\n",iColumn,value,columnUpper[iColumn]);
	    value=columnUpper[iColumn];
	  } else if (value<columnLower[iColumn]) {
	    if (value-columnLower[iColumn]<-1.0e-8)
	      printf("column %d has value %.12g below %.12g\n",iColumn,value,columnLower[iColumn]);
	    value=columnLower[iColumn];
	  }
	  if (value) {
	    CoinBigIndex start = columnStart[iColumn];
	    CoinBigIndex end = start + columnLength[iColumn];
	    for (CoinBigIndex j = start; j < end; j++) {
	      int iRow = row[j];
	      if (fabs(value)<1.0e-6&&fabs(value*element[j])>1.0e-5)
		printf("Column %d row %d value %.8g element %g %s\n",
		       iColumn,iRow,value,element[j],lp->isInteger(iColumn) ? "integer" : "");
	      rowActivity[iRow] += value * element[j];
	      rowSum[iRow] += fabs(value * element[j]);
	    }
	  }
	}
	for (int i = 0 ; i < numberRows ; i++) {
#if 0 //def CLP_INVESTIGATE
	  double inf;
	  inf = rowLower[i] - rowActivity[i];
	  if (inf > primalTolerance)
	    printf("Row %d inf %g sum %g %g <= %g <= %g\n",
		   i, inf, rowSum[i], rowLower[i], rowActivity[i], rowUpper[i]);
	  inf = rowActivity[i] - rowUpper[i];
	  if (inf > primalTolerance)
	    printf("Row %d inf %g sum %g %g <= %g <= %g\n",
		   i, inf, rowSum[i], rowLower[i], rowActivity[i], rowUpper[i]);
#endif
	  double infeasibility = CoinMax(rowActivity[i]-rowUpper[i],
					 rowLower[i]-rowActivity[i]);
	  // but allow for errors
	  double factor = CoinMax(1.0,rowSum[i]*1.0e-3);
	  if (infeasibility>largestInfeasibility*factor) {
	    largestInfeasibility = infeasibility/factor;
	    printf("Cinf of %g on row %d sum %g scaled %g\n",
		   infeasibility,i,rowSum[i],largestInfeasibility);
	    if (infeasibility>1.0e10) {
	      for (CoinBigIndex j=rowStart[i];
		   j<rowStart[i]+rowLength[i];j++) {
		printf("col %d element %g\n",
		       column[j],elementByRow[j]);
	      }
	    }
	  }
	}
	delete [] rowActivity ;
	delete [] rowSum;
	if (largestInfeasibility > 10.0*primalTolerance)
	  printf("Clargest infeasibility is %g - obj %g\n", largestInfeasibility,objValue);
	else
	  printf("Cfeasible (%g) - obj %g\n", largestInfeasibility,objValue);
      }
#endif
    for (int i = 0; (i < lp->getNumCols()); ++i) {
#if 0
         if (sol[i]<1e-8)
            sol[i] = 0.0;
         else
            if (lp->isInteger(i))
               sol[i] = floor( sol[i]+0.5 );
#else
      if (lp->isInteger(i)) {
        //if (fabs(sol[i] - floor( sol[i]+0.5 ))>1.0e-8)
        //printf("bad sol for %d - %.12g\n",i,sol[i]);
        sol[i] = floor(sol[i] + 0.5);
      }
#endif
    }
#if 0
      {
	int numberColumns=lp->getNumCols();
	double largestInfeasibility = 0.0;
	double primalTolerance ;
	double offset;
	lp->getDblParam(OsiObjOffset, offset);
	lp->getDblParam(OsiPrimalTolerance, primalTolerance) ;
	const double *objective = lp->getObjCoefficients() ;
	const double * rowLower = lp->getRowLower() ;
	const double * rowUpper = lp->getRowUpper() ;
	const double * columnLower = lp->getColLower() ;
	const double * columnUpper = lp->getColUpper() ;
	int numberRows = lp->getNumRows() ;
	double *rowActivity = new double[numberRows] ;
	memset(rowActivity, 0, numberRows*sizeof(double)) ;
	double *rowSum = new double[numberRows] ;
	memset(rowSum, 0, numberRows*sizeof(double)) ;
	const double * element = lp->getMatrixByCol()->getElements();
	const int * row = lp->getMatrixByCol()->getIndices();
	const CoinBigIndex * columnStart = lp->getMatrixByCol()->getVectorStarts();
	const int * columnLength = lp->getMatrixByCol()->getVectorLengths();
	const CoinPackedMatrix * rowCopy = lp->getMatrixByRow();
	const int * column = rowCopy->getIndices();
	const int * rowLength = rowCopy->getVectorLengths();
	const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
	const double * elementByRow = rowCopy->getElements();
	double objValue=-offset;
	for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	  double value = sol[iColumn];
	  if (lp->isInteger(iColumn))
	    assert (fabs(value-floor(value+0.5))<1.0e-6);
	  objValue += value*objective[iColumn];
	  if (value>columnUpper[iColumn]) {
	    if (value-columnUpper[iColumn]>1.0e-8)
	      printf("column %d has value %.12g above %.12g\n",iColumn,value,columnUpper[iColumn]);
	    value=columnUpper[iColumn];
	  } else if (value<columnLower[iColumn]) {
	    if (value-columnLower[iColumn]<-1.0e-8)
	      printf("column %d has value %.12g below %.12g\n",iColumn,value,columnLower[iColumn]);
	    value=columnLower[iColumn];
	  }
	  if (value) {
	    CoinBigIndex start = columnStart[iColumn];
	    CoinBigIndex end = start + columnLength[iColumn];
	    for (CoinBigIndex j = start; j < end; j++) {
	      int iRow = row[j];
	      rowActivity[iRow] += value * element[j];
	      rowSum[iRow] += fabs(value * element[j]);
	    }
	  }
	}
	for (int i = 0 ; i < numberRows ; i++) {
#if 0 //def CLP_INVESTIGATE
	  double inf;
	  inf = rowLower[i] - rowActivity[i];
	  if (inf > primalTolerance)
	    printf("Row %d inf %g sum %g %g <= %g <= %g\n",
		   i, inf, rowSum[i], rowLower[i], rowActivity[i], rowUpper[i]);
	  inf = rowActivity[i] - rowUpper[i];
	  if (inf > primalTolerance)
	    printf("Row %d inf %g sum %g %g <= %g <= %g\n",
		   i, inf, rowSum[i], rowLower[i], rowActivity[i], rowUpper[i]);
#endif
	  double infeasibility = CoinMax(rowActivity[i]-rowUpper[i],
					 rowLower[i]-rowActivity[i]);
	  // but allow for errors
	  double factor = CoinMax(1.0,rowSum[i]*1.0e-3);
	  if (infeasibility>largestInfeasibility*factor) {
	    largestInfeasibility = infeasibility/factor;
	    printf("Dinf of %g on row %d sum %g scaled %g\n",
		   infeasibility,i,rowSum[i],largestInfeasibility);
	    if (infeasibility>1.0e10) {
	      for (CoinBigIndex j=rowStart[i];
		   j<rowStart[i]+rowLength[i];j++) {
		printf("col %d element %g\n",
		       column[j],elementByRow[j]);
	      }
	    }
	  }
	}
	delete [] rowActivity ;
	delete [] rowSum;
	if (largestInfeasibility > 10.0*primalTolerance)
	  printf("Dlargest infeasibility is %g - obj %g\n", largestInfeasibility,objValue);
	else
	  printf("Dfeasible (%g) - obj %g\n", largestInfeasibility,objValue);
      }
#endif
#if JUST_FIX_INTEGER
    const double *oldLower = model->solver()->getColLower();
    const double *oldUpper = model->solver()->getColUpper();
    const double *dj = lp->getReducedCost();
    int nNaturalLB = 0;
    int nMaybeLB = 0;
    int nForcedLB = 0;
    int nNaturalUB = 0;
    int nMaybeUB = 0;
    int nForcedUB = 0;
    int nOther = 0;
    for (int i = 0; i < lp->getNumCols(); ++i) {
      if (lp->isInteger(i)) {
        if (sol[i] == oldLower[i]) {
          if (dj[i] > 1.0e-5)
            nNaturalLB++;
          else if (dj[i] < -1.0e-5)
            nForcedLB++;
          else
            nMaybeLB++;
        } else if (sol[i] == oldUpper[i]) {
          if (dj[i] < -1.0e-5)
            nNaturalUB++;
          else if (dj[i] > 1.0e-5)
            nForcedUB++;
          else
            nMaybeUB++;
        } else {
          nOther++;
        }
      }
    }
    printf("%d other, LB %d natural, %d neutral, %d forced, UB %d natural, %d neutral, %d forced\n",
      nOther, nNaturalLB, nMaybeLB, nForcedLB,
      nNaturalUB, nMaybeUB, nForcedUB = 0);
#endif
  }

TERMINATE:
  delete[] realObj;
  delete lp;
  return status;
}
#undef STR_SIZE

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
