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

bool isNumericStr(const char *str)
{
  const size_t l = strlen(str);

  for (size_t i = 0; i < l; ++i)
    if (!(isdigit(str[i]) || (str[i] == '.') || (str[i] == '-') || (str[i] == '+') || (str[i] == 'e')))
      return false;

  return true;
}

int CbcMipStartIO::read(OsiSolverInterface *solver, const char *fileName,
  std::vector< std::pair< std::string, double > > &colValues,
  double &solObj, CoinMessageHandler *messHandler, CoinMessages *pcoinmsgs)
{
  CoinMessages &messages = *pcoinmsgs;
#define STR_SIZE 256
  char printLine[STR_SIZE] = "";
  FILE *f = fopen(fileName, "r");
  if (!f) {
    sprintf(printLine, "Unable to open file %s.", fileName);
    messHandler->message(CBC_GENERAL, messages) << printLine << CoinMessageEol;
    return 1;
  }
  char line[STR_SIZE] = "";

  int nLine = 0;
  // check if psv format!
  int lengthFilename = strlen(fileName);
  // separator
  char separator = ' ';
  if (strstr(fileName,".psv") == fileName+lengthFilename-4)
    separator = '|';
  else if (strstr(fileName,".csv") == fileName+lengthFilename-4) 
    separator = ',';
  if (separator==' ') {
    // ordinary
    while (fgets(line, STR_SIZE, f)) {
      ++nLine;
      char col[4][STR_SIZE] = {"", "", "", ""};
      int nread = sscanf(line, "%s %s %s %s", col[0], col[1], col[2], col[3]);
      if (!nread)
	continue;
      /* line with variable value */
      if (strlen(col[0]) && isdigit(col[0][0]) && (nread >= 3)) {
	if (!isNumericStr(col[0])) {
	  sprintf(printLine, "Reading: %s, line %d - first column in mipstart file should be numeric, ignoring.", fileName, nLine);
	  messHandler->message(CBC_GENERAL, messages) << printLine << CoinMessageEol;
	  continue;
	}
	if (!isNumericStr(col[2])) {
	  sprintf(printLine, "Reading: %s, line %d - Third column in mipstart file should be numeric, ignoring.", fileName, nLine);
	  messHandler->message(CBC_GENERAL, messages) << printLine << CoinMessageEol;
	  continue;
	}
	
	char *name = col[1];
	double value = atof(col[2]);
	
	colValues.push_back(pair< string, double >(string(name), value));
      }
    }
  } else {
    // csv or psv
    int nBad1 = 0;
    int nBad2 = 0;
    while (fgets(line, STR_SIZE, f)) {
      ++nLine;
      // clean line
      // out \n \r and blanks
      int n = strlen(line);
      int nNew = 0;
      for (int i=0;i<n;i++) {
	char charX = line[i];
	if (charX==' ') {
	  continue;
	} else if (charX=='\n'||charX=='\r') {
	  line[nNew]='\0';
	  break;
	} else {
	  line[nNew++] = charX;
	}
      }
      char * pipeorcomma = strchr(line,separator);
      if (!pipeorcomma) {
	if (!nBad1) {
	  if (nLine>1) {
	    sprintf(printLine, "Reading: %s, line %d (%s) - mipstart file should contain |.", fileName, nLine,line);
	    messHandler->message(CBC_GENERAL, messages) << printLine << CoinMessageEol;
	  } else {
	    // may be OK
	    nBad1--;
	  }
	}
	nBad1++;
	continue;
      }
      *pipeorcomma = '\0';
      if (!isNumericStr(pipeorcomma+1)) {
	if (!nBad2) {
	  if (nLine>1) {
	    sprintf(printLine, "Reading: %s, line %d (%s) - Second column in mipstart file should be numeric.", fileName, nLine,line);
	    messHandler->message(CBC_GENERAL, messages) << printLine << CoinMessageEol;
	  } else {
	    // may be OK
	    nBad2--;
	  }
	}
	nBad2++;
	continue;
      }
	
      double value = atof(pipeorcomma+1);
      
      colValues.push_back(pair< string, double >(string(line), value));
    }
    if (nBad1||nBad2) {
      sprintf(printLine, "Reading: %s, %d errors.", fileName, nBad1+nBad2);
      messHandler->message(CBC_GENERAL, messages) << printLine << CoinMessageEol;
      return 1;
    }
  }

  const int numCols = solver->getNumCols();
  if (colValues.size()) {
    sprintf(printLine, "MIPStart values read for %d variables.", static_cast< int >(colValues.size()));
    messHandler->message(CBC_GENERAL, messages) << printLine << CoinMessageEol;
    vector< pair< string, double > > fullValues;
    /* for fast search of column names */
    map< string, int > colIdx;
    for (int i = 0; i < numCols; i++) {
      fullValues.push_back(pair< string, double >(solver->getColName(i), 0.0));
      colIdx[solver->getColName(i)] = i;
    }
    const double *lower = solver->getColLower();
    const double *upper = solver->getColUpper();
    int nBadValues = 0;
    for (int i = 0; (i < static_cast< int >(colValues.size())); ++i) {
      map< string, int >::const_iterator mIt = colIdx.find(colValues[i].first);
      if (mIt != colIdx.end()) {
        const int idx = mIt->second;
        double v = colValues[i].second;
	if (v>upper[idx]) {
	  nBadValues++;
	  v = upper[idx];
	} else if (v<lower[idx]) {
	  nBadValues++;
	  v = lower[idx];
	}
        fullValues[idx].second = v;
      }
    }
    if (nBadValues) {
      sprintf(printLine,"Warning: modifying %d solution values outside bounds",
	      nBadValues);
      messHandler->message(CBC_GENERAL, messages) << printLine << CoinMessageEol;
    }
    colValues = fullValues;
  } else {
    sprintf(printLine, "File %s does not contains a solution.", fileName);
    messHandler->message(CBC_GENERAL, messages) << printLine << CoinMessageEol;
    fclose(f);
    return 1;
  }

  fclose(f);
  return 0;
}

int CbcMipStartIO::computeCompleteSolution(CbcModel *model, OsiSolverInterface *solver,
  const std::vector< std::string > colNames,
  const std::vector< std::pair< std::string, double > > &colValues,
  double *sol, double &obj, int extraActions, CoinMessageHandler *messHandler, CoinMessages *pmessages)
{
  if (!solver->getNumCols())
    return 0;
  
  CoinMessages &messages = *pmessages;

  int status = 0;
  double compObj = COIN_DBL_MAX;
  bool foundIntegerSol = false;
  OsiSolverInterface *lp = solver->clone();

  map< string, int > colIdx;
  assert((static_cast< int >(colNames.size())) == lp->getNumCols());
  /* for fast search of column names */
  for (int i = 0; (i < static_cast< int >(colNames.size())); ++i)
    colIdx[colNames[i]] = i;

  char printLine[STR_SIZE];
  int fixed = 0;
  int notFound = 0;
  char colNotFound[256] = "";
  int nContinuousFixed = 0;
  double *realObj = new double[lp->getNumCols()];
  memcpy(realObj, lp->getObjCoefficients(), sizeof(double)*lp->getNumCols());

  // assuming that variables not fixed are more likely to have zero as value,
  // inserting as default objective function 1
  if (0) { // to get more accurate answers
    vector< double > obj(lp->getNumCols(), lp->getObjSense());
    lp->setObjective(&obj[0]);
  }

#ifndef JUST_FIX_INTEGER
#define JUST_FIX_INTEGER 2
#endif

#if JUST_FIX_INTEGER > 1
  // all not mentioned are at zero
  for (int i = 0; (i < lp->getNumCols()); ++i) {
    if (lp->isInteger(i))
      lp->setColBounds(i, 0.0, 0.0);
  }
#endif
  if (extraActions) {
    const double * objective = lp->getObjCoefficients();
    const double * lower = lp->getColLower();
    const double * upper = lp->getColUpper();
    for (int i = 0; (i < lp->getNumCols()); ++i) {
      if (lp->isInteger(i)) {
        double objValue = objective[i];
        double lowerValue = lower[i];
        double upperValue = upper[i];
        switch (extraActions) {
        case 1:
          lp->setColBounds(i, lowerValue, lowerValue);
          break;
        case 2:
          lp->setColBounds(i, upperValue, upperValue);
          break;
        case 3:
          lp->setColBounds(i, lowerValue, lowerValue);
          if (objValue<0.0)
            lp->setColBounds(i, upperValue, upperValue);
          break;
        case 4:
          lp->setColBounds(i, upperValue, upperValue);
          if (objValue>0.0)
            lp->setColBounds(i, lowerValue, lowerValue);
          break;
        case 5:
          lp->setColBounds(i, lowerValue, lowerValue);
          if (objValue>0.0)
            lp->setColBounds(i, upperValue, upperValue);
          break;
        case 6:
          lp->setColBounds(i, upperValue, upperValue);
          if (objValue<0.0)
            lp->setColBounds(i, lowerValue, lowerValue);
          break;
        }
      }
    }
  }
  for (int i = 0; (i < static_cast< int >(colValues.size())); ++i) {
    map< string, int >::const_iterator mIt = colIdx.find(colValues[i].first);
    if (mIt == colIdx.end()) {
      if (!notFound)
        strcpy(colNotFound, colValues[i].first.c_str());
      notFound++;
    } else {
      const int idx = mIt->second;
      double v = colValues[i].second;
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

  if (extraActions)
    fixed = lp->getNumIntegers();
  if (!fixed) {
    messHandler->message(CBC_GENERAL, messages)
      << "Warning: MIPstart solution is not valid, column names do not match, ignoring it."
      << CoinMessageEol;
    goto TERMINATE;
  }

  if (notFound >= ((static_cast< double >(colNames.size())) * 0.5)) {
    sprintf(printLine, "Warning: %d column names were not found (e.g. %s) while filling solution.", notFound, colNotFound);
    messHandler->message(CBC_GENERAL, messages)
      << printLine << CoinMessageEol;
  }
#if JUST_FIX_INTEGER
  lp->setHintParam(OsiDoPresolveInInitial, true, OsiHintDo);
#endif

  //lp->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
  lp->initialSolve();

  if ((lp->isProvenPrimalInfeasible()) || (lp->isProvenDualInfeasible())) {
    if (nContinuousFixed) {
      messHandler->message(CBC_GENERAL, messages)
        << "Trying just fixing integer variables (and fixingish SOS)." << CoinMessageEol;
      int numberColumns = lp->getNumCols();
      const double *oldLower = solver->getColLower();
      const double *oldUpper = solver->getColUpper();
      double *savedSol = CoinCopyOfArray(lp->getColLower(), numberColumns);
      for (int i = 0; i < numberColumns; ++i) {
        if (!lp->isInteger(i)) {
          lp->setColLower(i, oldLower[i]);
          lp->setColUpper(i, oldUpper[i]);
        }
      }
      // but look at SOS
      if (model) {
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
      }

      delete[] savedSol;
      lp->initialSolve();
    } else {
      messHandler->message(CBC_GENERAL, messages)
        << "Fixing only non-zero variables." << CoinMessageEol;
      /* unfix all variables which are zero */
      int notZeroAnymore = 0;
      for (int i = 0; (i < lp->getNumCols()); ++i)
        if (((fabs(lp->getColLower()[i])) <= 1e-8) && (fabs(lp->getColLower()[i] - lp->getColUpper()[i]) <= 1e-8)) {
          const double *oldLower = solver->getColLower();
          const double *oldUpper = solver->getColUpper();
          lp->setColLower(i, oldLower[i]);
          lp->setColUpper(i, oldUpper[i]);
          notZeroAnymore++;
        }
      if (notZeroAnymore)
        lp->initialSolve();
    }
  }

  if (!lp->isProvenOptimal()) {
    messHandler->message(CBC_GENERAL, messages)
      << "Warning: mipstart values could not be used to build a solution." << CoinMessageEol;
    status = 1;
    goto TERMINATE;
  }
  /* some additional effort is needed to provide an integer solution */
  if (lp->getFractionalIndices().size() > 0) {
    sprintf(printLine, "MIPStart solution provided values for %d of %d integer variables, %d variables are still fractional.", fixed, lp->getNumIntegers(), static_cast< int >(lp->getFractionalIndices().size()));
    messHandler->message(CBC_GENERAL, messages)
      << printLine << CoinMessageEol;
    if (lp->getFractionalIndices().size()<5) {
      for (int i=0;i<lp->getFractionalIndices().size();i++) {
	int iColumn = lp->getFractionalIndices()[i];
	sprintf(printLine, "Variable %d %s has value %g",iColumn,
		colNames[iColumn].c_str(),lp->getColSolution()[iColumn]);
	messHandler->message(CBC_GENERAL, messages)
	  << printLine << CoinMessageEol;
      }
    }
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
      messHandler->message(CBC_GENERAL, messages)
        << printLine << CoinMessageEol;
      foundIntegerSol = true;
      lp->getDblParam(OsiObjOffset, obj);
      obj = -obj;
      for ( int i=0 ; (i<lp->getNumCols()) ; ++i )
          obj += realObj[i]*sol[i];
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
      messHandler->message(CBC_GENERAL, messages)
        << printLine << CoinMessageEol;
      copy(babModel.bestSolution(), babModel.bestSolution() + babModel.getNumCols(), sol);
      foundIntegerSol = true;
      obj = compObj = babModel.getObjValue();
    }
#endif
    else {
      messHandler->message(CBC_GENERAL, messages)
        << "Warning: mipstart values could not be used to build a solution." << CoinMessageEol;
      status = 1;
      goto TERMINATE;
    }
  } else {
    foundIntegerSol = true;
    
    lp->getDblParam(OsiObjOffset, obj);
    obj = -obj;
    for ( int i=0 ; (i<lp->getNumCols()) ; ++i )
      obj += realObj[i]*lp->getColSolution()[i];
    compObj = obj;
    copy(lp->getColSolution(), lp->getColSolution() + lp->getNumCols(), sol);
  }

  if (foundIntegerSol) {
    sprintf(printLine, "MIPStart provided solution with cost %g", compObj);
    messHandler->message(CBC_GENERAL, messages)
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
    const double *oldLower = solver->getColLower();
    const double *oldUpper = solver->getColUpper();
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
    //printf("%d other, LB %d natural, %d neutral, %d forced, UB %d natural, %d neutral, %d forced\n",
    //nOther, nNaturalLB, nMaybeLB, nForcedLB,
    //nNaturalUB, nMaybeUB, nForcedUB = 0);
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
