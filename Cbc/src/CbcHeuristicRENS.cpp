// $Id$
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// edwin 12/5/09 carved out of CbcHeuristicRINS

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcHeuristicRENS.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinSort.hpp"
#include "CbcBranchActual.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"

// Default Constructor
CbcHeuristicRENS::CbcHeuristicRENS()
        : CbcHeuristic()
{
    numberTries_ = 0;
    rensType_ = 0;
    whereFrom_ = 256 + 1;
}

// Constructor with model - assumed before cuts

CbcHeuristicRENS::CbcHeuristicRENS(CbcModel & model)
        : CbcHeuristic(model)
{
    numberTries_ = 0;
    rensType_ = 0;
    whereFrom_ = 256 + 1;
}

// Destructor
CbcHeuristicRENS::~CbcHeuristicRENS ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicRENS::clone() const
{
    return new CbcHeuristicRENS(*this);
}

// Assignment operator
CbcHeuristicRENS &
CbcHeuristicRENS::operator=( const CbcHeuristicRENS & rhs)
{
    if (this != &rhs) {
        CbcHeuristic::operator=(rhs);
        numberTries_ = rhs.numberTries_;
	rensType_ = rhs.rensType_;
    }
    return *this;
}

// Copy constructor
CbcHeuristicRENS::CbcHeuristicRENS(const CbcHeuristicRENS & rhs)
        :
        CbcHeuristic(rhs),
        numberTries_(rhs.numberTries_),
	rensType_(rhs.rensType_)
{
}
// Resets stuff if model changes
void
CbcHeuristicRENS::resetModel(CbcModel * )
{
}
int
CbcHeuristicRENS::solution(double & solutionValue,
                           double * betterSolution)
{
    int returnCode = 0;
    const double * bestSolution = model_->bestSolution();
    if ((numberTries_&&(rensType_&16)==0) || numberTries_>1 || (when() < 2 && bestSolution))
        return 0;
    numberTries_++;
    double saveFractionSmall=fractionSmall_;
    OsiSolverInterface * solver = model_->solver();

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();

    OsiSolverInterface * newSolver = cloneBut(3); // was model_->continuousSolver()->clone();
    const double * currentSolution = newSolver->getColSolution();
    int type = rensType_&15;
    if (type<12)
      newSolver->resolve();
    double direction = newSolver->getObjSense();
    double cutoff=model_->getCutoff();
    newSolver->setDblParam(OsiDualObjectiveLimit, 1.0e100);
    //cutoff *= direction;
    double gap = cutoff - newSolver->getObjValue() * direction ;
    double tolerance;
    newSolver->getDblParam(OsiDualTolerance, tolerance) ;
    if ((gap > 0.0 || !newSolver->isProvenOptimal())&&type<12) {
      gap += 100.0 * tolerance;
      int nFix = newSolver->reducedCostFix(gap);
      if (nFix) {
	char line [200];
	sprintf(line, "Reduced cost fixing fixed %d variables", nFix);
	model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
	  << line
	  << CoinMessageEol;
      }
    } else if (type<12) {
      return 0; // finished?
    }
    int numberColumns = solver->getNumCols();
    double * dj = CoinCopyOfArray(solver->getReducedCost(),numberColumns);
    double djTolerance = (type!=1) ? -1.0e30 : 1.0e-4;
    const double * colLower = newSolver->getColLower();
    const double * colUpper = newSolver->getColUpper();
    double * contribution = NULL;
    int numberFixed = 0;
    if (type==3) {
      double total=0.0;
      int n=0;
      CoinWarmStartBasis * basis =
	dynamic_cast<CoinWarmStartBasis *>(solver->getWarmStart()) ;
      if (basis&&basis->getNumArtificial()) {
	for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	  if (colUpper[iColumn]>colLower[iColumn]&&
	      basis->getStructStatus(iColumn) !=
	      CoinWarmStartBasis::basic) {
	    n++;
	    total += fabs(dj[iColumn]);
	  }
	}
	if (n)
	  djTolerance = (0.01*total)/static_cast<double>(n);
	delete basis;
      }
    } else if (type>=5&&type<=12) {
      /* 5 fix sets at one
	 6 fix on dj but leave unfixed SOS slacks
	 7 fix sets at one but use pi
	 8 fix all at zero but leave unfixed SOS slacks
	 9 as 8 but only fix all at zero if just one in set nonzero
	 10 fix all "stable" ones
	 11 fix all "stable" ones - approach 2
	 12 layered approach
      */
      // SOS type fixing
      bool fixSets = (type==5)||(type==7)||(type==10)||(type==11);
      CoinWarmStartBasis * basis =
	dynamic_cast<CoinWarmStartBasis *>(solver->getWarmStart()) ;
      if (basis&&basis->getNumArtificial()) {
	//const double * rowLower = solver->getRowLower();
	const double * rowUpper = solver->getRowUpper();
	
	int numberRows = solver->getNumRows();
	// Column copy
	const CoinPackedMatrix * matrix = solver->getMatrixByCol();
	const double * element = matrix->getElements();
	const int * row = matrix->getIndices();
	const CoinBigIndex * columnStart = matrix->getVectorStarts();
	const int * columnLength = matrix->getVectorLengths();
	double * bestDj = new double [numberRows];
	for (int i=0;i<numberRows;i++) {
	  if (rowUpper[i]==1.0)
	    bestDj[i]=1.0e20;
	  else
	    bestDj[i]=1.0e30;
	}
	for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	  if (colUpper[iColumn]>colLower[iColumn]) {
	    CoinBigIndex j;
	    if (currentSolution[iColumn]>1.0e-6&&
		currentSolution[iColumn]<0.999999) {
	      for (j = columnStart[iColumn];
		   j < columnStart[iColumn] + columnLength[iColumn]; j++) {
		int iRow = row[j];
		if (bestDj[iRow]<1.0e30) {
		  if (element[j] != 1.0)
		    bestDj[iRow]=1.0e30;
		  else
		    bestDj[iRow]=1.0e25;
		}
	      }
	    } else if ( basis->getStructStatus(iColumn) !=
	      CoinWarmStartBasis::basic) {
	      for (j = columnStart[iColumn];
		   j < columnStart[iColumn] + columnLength[iColumn]; j++) {
		int iRow = row[j];
		if (bestDj[iRow]<1.0e25) {
		  if (element[j] != 1.0)
		    bestDj[iRow]=1.0e30;
		  else
		    bestDj[iRow]=CoinMin(fabs(dj[iColumn]),bestDj[iRow]);
		}
	      }
	    }
	  }
	}
	// Just leave one slack in each set
	{
	  const double * objective = newSolver->getObjCoefficients();
	  int * best = new int [numberRows];
	  double * cheapest = new double[numberRows];
	  for (int i=0;i<numberRows;i++) {
	    best[i]=-1;
	    cheapest[i]=COIN_DBL_MAX;
	  }
	  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	    if (colUpper[iColumn]>colLower[iColumn]) {
	      if (columnLength[iColumn]==1) {
		CoinBigIndex j = columnStart[iColumn];
		int iRow = row[j];
		if (bestDj[iRow]<1.0e30) {
		  double obj = direction*objective[iColumn];
		  if (obj<cheapest[iRow]) {
		    cheapest[iRow]=obj;
		    best[iRow]=iColumn;
		  }
		}
	      }
	    }
	  }
	  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	    if (colUpper[iColumn]>colLower[iColumn]) {
	      if (columnLength[iColumn]==1) {
		CoinBigIndex j = columnStart[iColumn];
		int iRow = row[j];
		if (bestDj[iRow]<1.0e30) {
		  if (best[iRow]!=-1&&iColumn!=best[iRow]) {
		    newSolver->setColUpper(iColumn,0.0);
		  }
		}
	      }
	    }
	  }
	  delete [] best;
	  delete [] cheapest;
	}
	int nSOS=0;
	double * sort = new double [numberRows];
	const double * pi = newSolver->getRowPrice();
	if (type==12) {
	  contribution = new double [numberRows];
	  for (int i=0;i<numberRows;i++) {
	    if (bestDj[i]<1.0e30) 
	      contribution[i]=0.0;
	    else
	      contribution[i]=-1.0;
	  }
	}
	for (int i=0;i<numberRows;i++) {
	  if (bestDj[i]<1.0e30) {
	    if (type==5)
	      sort[nSOS++]=bestDj[i];
	    else if (type==7)
	      sort[nSOS++]=-fabs(pi[i]);
	    else
	      sort[nSOS++]=fabs(pi[i]);
	  }
	}
	if (10*nSOS>8*numberRows) {
	  if (type<10) {
	    std::sort(sort,sort+nSOS);
	    int last = static_cast<int>(nSOS*0.9*fractionSmall_);
	    double tolerance = sort[last];
	    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	      if (colUpper[iColumn]>colLower[iColumn]) {
		CoinBigIndex j;
		if (currentSolution[iColumn]<=1.0e-6||
		    currentSolution[iColumn]>=0.999999) {
		  if (fixSets) {
		    for (j = columnStart[iColumn];
			 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
		      int iRow = row[j];
		      double useDj;
		      if (type==5) 
			useDj = bestDj[iRow];
		      else if (type==7)
			useDj= -fabs(pi[iRow]);
		      else
			useDj= fabs(pi[iRow]);
		      if (bestDj[iRow]<1.0e30&&useDj>=tolerance) {
			numberFixed++;
			if (currentSolution[iColumn]<=1.0e-6)
			  newSolver->setColUpper(iColumn,0.0);
			else if (currentSolution[iColumn]>=0.999999) 
			  newSolver->setColLower(iColumn,1.0);
		      }
		    }
		  } else if (columnLength[iColumn]==1) {
		    // leave more slacks
		    int iRow = row[columnStart[iColumn]];
		    if (bestDj[iRow]<1.0e30) {
		      // fake dj
		      dj[iColumn] *= 0.000001;
		    }
		  } else if (type==8||type==9) {
		    if (currentSolution[iColumn]<=1.0e-6) {
		      if (type==8) {
			dj[iColumn] *= 1.0e6;
		      } else {
			bool fix=false;
			for (j = columnStart[iColumn];
			     j < columnStart[iColumn] + columnLength[iColumn]; j++) {
			  int iRow = row[j];
			  if (bestDj[iRow]<1.0e25) {
			    fix=true;
			    break;
			  }
			}
			if (fix) {
			  dj[iColumn] *= 1.0e6;
			}
		      }
		    } else {
		      dj[iColumn] *= 0.000001;
		    }
		  }
		}
	      }
	    }
	    if (fixSets)
	      djTolerance = 1.0e30;
	  } else if (type==10) {
	    double * saveUpper = new double [numberRows];
	    memset(saveUpper,0,numberRows*sizeof(double));
	    char * mark = new char [numberColumns];
	    char * nonzero = new char [numberColumns];
	    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	      if (colUpper[iColumn]>colLower[iColumn]) {
		CoinBigIndex j;
		for (j = columnStart[iColumn];
		     j < columnStart[iColumn] + columnLength[iColumn]; j++) {
		  int iRow = row[j];
		  saveUpper[iRow] += element[j];
		}
	      }
	    }
	    double sum=0.0;
	    double sumRhs=0.0;
	    const double * rowUpper = newSolver->getRowUpper();
	    for (int i=0;i<numberRows;i++) {
	      if (bestDj[i]>=1.0e30) {
		sum += saveUpper[i];
		sumRhs += rowUpper[i];
	      }
	    }
	    double averagePerSet = sum/static_cast<double>(numberRows);
	    // allow this extra
	    double factor = averagePerSet*fractionSmall_*numberRows;
	    factor = 1.0+factor/sumRhs;
	    fractionSmall_ = 0.5;
	    memcpy(saveUpper,rowUpper,numberRows*sizeof(double));
	    // loosen up
	    for (int i=0;i<numberRows;i++) {
	      if (bestDj[i]>=1.0e30) {
		newSolver->setRowUpper(i,factor*saveUpper[i]);
	      }
	    }
	    newSolver->resolve();
	    const double * solution = newSolver->getColSolution();
	    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	      mark[iColumn]=0;
	      nonzero[iColumn]=0;
	      if (colUpper[iColumn]>colLower[iColumn]&&
		  solution[iColumn]>0.9999)
		mark[iColumn]=1;
	      else if (solution[iColumn]>0.00001)
		nonzero[iColumn]=1;
	    }
	    // slightly small
	    for (int i=0;i<numberRows;i++) {
	      if (bestDj[i]>=1.0e30) {
		newSolver->setRowUpper(i,saveUpper[i]*0.9999);
	      }
	    }
	    newSolver->resolve();
	    int nCheck=2;
	    if (newSolver->isProvenOptimal()) {
	      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
		if (colUpper[iColumn]>colLower[iColumn]&&
		    solution[iColumn]>0.9999)
		  mark[iColumn]++;
		else if (solution[iColumn]>0.00001)
		  nonzero[iColumn]=1;
	      }
	    } else {
	      nCheck=1;
	    }
	    // correct values
	    for (int i=0;i<numberRows;i++) {
	      if (bestDj[i]>=1.0e30) {
		newSolver->setRowUpper(i,saveUpper[i]);
	      }
	    }
	    newSolver->resolve();
	    int nFixed=0;
	    int nFixedToZero=0;
	    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	      if (colUpper[iColumn]>colLower[iColumn]) {
		if (solution[iColumn]>0.9999&&mark[iColumn]==nCheck) {
		  newSolver->setColLower(iColumn,1.0);
		  nFixed++;
		} else if (!mark[iColumn]&&!nonzero[iColumn]&&
			   columnLength[iColumn]>1&&solution[iColumn]<0.00001) {
		  newSolver->setColUpper(iColumn,0.0);
		  nFixedToZero++;
		}
	      }
	    }
	    char line[100];
	    sprintf(line,"Heuristic %s fixed %d to one (and %d to zero)",
		    heuristicName(),
		    nFixed,nFixedToZero);
	    model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
	      << line
	      << CoinMessageEol;
	    delete [] mark;
	    delete []nonzero;
	    delete [] saveUpper;
	    numberFixed=numberColumns;
	    djTolerance = 1.0e30;
	  } else if (type==11) {
	    double * saveUpper = CoinCopyOfArray(newSolver->getRowUpper(),numberRows);
	    char * mark = new char [numberColumns];
	    char * nonzero = new char [numberColumns];
	    // save basis and solution
	    CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis*>(newSolver->getWarmStart()) ;
	    assert(basis != NULL);
	    double * saveSolution = 
	      CoinCopyOfArray(newSolver->getColSolution(), 
			      numberColumns);
	    double factors[] = {1.1,1.05,1.01,0.98};
	    int nPass = (sizeof(factors)/sizeof(double))-1;
	    double factor=factors[0];
	    double proportion = fractionSmall_;
	    fractionSmall_ = 0.5;
	    // loosen up
	    for (int i=0;i<numberRows;i++) {
	      if (bestDj[i]>=1.0e30) {
		newSolver->setRowUpper(i,factor*saveUpper[i]);
	      }
	    }
            bool takeHint;
            OsiHintStrength strength;
            newSolver->getHintParam(OsiDoDualInResolve, takeHint, strength);
            newSolver->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
            newSolver->resolve();
            newSolver->setHintParam(OsiDoDualInResolve, true, OsiHintDo);
	    const double * solution = newSolver->getColSolution();
	    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	      mark[iColumn]=0;
	      nonzero[iColumn]=0;
	      if (colUpper[iColumn]>colLower[iColumn]&&
		  solution[iColumn]>0.9999)
		mark[iColumn]=1;
	      else if (solution[iColumn]>0.00001)
		nonzero[iColumn]=1;
	    }
	    int nCheck=2;
	    for (int iPass=0;iPass<nPass;iPass++) {
	      // smaller
	      factor = factors[iPass+1];
	      for (int i=0;i<numberRows;i++) {
		if (bestDj[i]>=1.0e30) {
		  newSolver->setRowUpper(i,saveUpper[i]*factor);
		}
	      }
	      newSolver->resolve();
	      if (newSolver->isProvenOptimal()) {
		nCheck++;
		for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
		  if (colUpper[iColumn]>colLower[iColumn]&&
		      solution[iColumn]>0.9999)
		    mark[iColumn]++;
		  else if (solution[iColumn]>0.00001)
		    nonzero[iColumn]++;
		}
	      }
	    }
	    // correct values
	    for (int i=0;i<numberRows;i++) {
	      if (bestDj[i]>=1.0e30) {
		newSolver->setRowUpper(i,saveUpper[i]);
	      }
	    }
	    newSolver->setColSolution(saveSolution);
	    delete [] saveSolution;
	    newSolver->setWarmStart(basis);
	    delete basis ;
            newSolver->setHintParam(OsiDoDualInResolve, takeHint, strength);
	    newSolver->resolve();
	    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	      if (colUpper[iColumn]>colLower[iColumn]&&
		  solution[iColumn]>0.9999)
		mark[iColumn]++;
	      else if (solution[iColumn]>0.00001)
		nonzero[iColumn]++;
	    }
	    int nFixed=0;
	    int numberSetsToFix = static_cast<int>(nSOS*(1.0-proportion));
	    int * mixed = new int[numberRows];
	    memset(mixed,0,numberRows*sizeof(int));
	    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	      if (colUpper[iColumn]>colLower[iColumn]) {
		int iSOS=-1;
		for (CoinBigIndex j = columnStart[iColumn];
		     j < columnStart[iColumn] + columnLength[iColumn]; j++) {
		  int iRow = row[j];
		  if (bestDj[iRow]<1.0e25) {
		    iSOS=iRow;
		    break;
		  }
		}
		if (iSOS>=0) {
		  int numberTimesAtOne = mark[iColumn];
		  int numberTimesNonZero = nonzero[iColumn]+
		    numberTimesAtOne;
		  if (numberTimesAtOne<nCheck&&
		      numberTimesNonZero) {
		    mixed[iSOS]+=
		      CoinMin(numberTimesNonZero,
			      nCheck-numberTimesNonZero);
		  } 
		}
	      }
	    }
	    int nFix=0;
	    for (int i=0;i<numberRows;i++) {
	      if (bestDj[i]<1.0e25) {
		sort[nFix] = -bestDj[i]+1.0e8*mixed[i];
		mixed[nFix++]=i;
	      }
	    }
	    CoinSort_2(sort,sort+nFix,mixed);
	    nFix = CoinMin(nFix,numberSetsToFix);
	    memset(sort,0,sizeof(double)*numberRows);
	    for (int i=0;i<nFix;i++)
	      sort[mixed[i]]=1.0;
	    delete [] mixed;
	    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	      if (colUpper[iColumn]>colLower[iColumn]) {
		if (solution[iColumn]>0.9999) {
		  int iSOS=-1;
		  for (CoinBigIndex j = columnStart[iColumn];
		       j < columnStart[iColumn] + columnLength[iColumn]; j++) {
		    int iRow = row[j];
		    if (bestDj[iRow]<1.0e25) {
		      iSOS=iRow;
		      break;
		    }
		  }
		  if (iSOS>=0&&sort[iSOS]) {
		    newSolver->setColLower(iColumn,1.0);
		    nFixed++;
		  }
		}
	      }
	    }
	    char line[100];
	    sprintf(line,"Heuristic %s fixed %d to one (%d sets)",
		    heuristicName(),
		    nFixed,nSOS);
	    model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
	      << line
	      << CoinMessageEol;
	    delete [] mark;
	    delete [] nonzero;
	    delete [] saveUpper;
	    numberFixed=numberColumns;
	    djTolerance = 1.0e30;
	  }
	}
	delete basis;
	delete [] sort;
	delete [] bestDj;
	if (10*nSOS<=8*numberRows) {
	  // give up
	  delete [] contribution;
	  delete newSolver;
	  return 0;
	}
      }
    }
    // Do dj to get right number
    if (type==4||type==6||(type>7&&type<10)) {
      double * sort = new double [numberColumns];
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	sort[iColumn]=1.0e30;
	if (colUpper[iColumn]>colLower[iColumn]) {
	  sort[iColumn] = fabs(dj[iColumn]);
	}
      }
      std::sort(sort,sort+numberColumns);
      int last = static_cast<int>(numberColumns*fractionSmall_);
      djTolerance = CoinMax(sort[last],1.0e-5);
      delete [] sort;
    } else if (type==12) {
      // Do layered in a different way
      int numberRows = solver->getNumRows();
      // Column copy
      const CoinPackedMatrix * matrix = newSolver->getMatrixByCol();
      const double * element = matrix->getElements();
      const int * row = matrix->getIndices();
      const CoinBigIndex * columnStart = matrix->getVectorStarts();
      const int * columnLength = matrix->getVectorLengths();
      int * whichRow = new int[numberRows];
      int * whichSet = new int [numberColumns];
      int nSOS=0;
      for (int i=0;i<numberRows;i++) {
	whichRow[i]=0;
	if (!contribution[i])
	  nSOS++;
      }
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	whichSet[iColumn]=-2;
	if (colUpper[iColumn]>colLower[iColumn]) {
	  CoinBigIndex j;
	  double sum=0.0;
	  int iSOS=-1;
	  int n=0;
	  for (j = columnStart[iColumn];
	       j < columnStart[iColumn] + columnLength[iColumn]; j++) {
	    int iRow = row[j];
	    if (contribution[iRow]>=0.0) {
	      iSOS=iRow;
	      n++;
	    } else {
	      sum += fabs(element[j]);
	    }
	  }
	  if (n>1)
	    COIN_DETAIL_PRINT(printf("Too many SOS entries (%d) for column %d\n",
				     n,iColumn));
	  if (sum) {
	    assert (iSOS>=0);
	    contribution[iSOS] += sum;
	    whichRow[iSOS]++;
	    whichSet[iColumn]=iSOS;
	  } else {
	    whichSet[iColumn]=iSOS+numberRows;
	  }
	}
      }
      int * chunk = new int [numberRows];
      for (int i=0;i<numberRows;i++) {
	chunk[i]=-1;
	if (whichRow[i]) {
	  contribution[i]= - contribution[i]/static_cast<double>(whichRow[i]);
	} else {
	  contribution[i] = COIN_DBL_MAX;
	}
	whichRow[i]=i;
      }
      newSolver->setDblParam(OsiDualObjectiveLimit, 1.0e100);
      double * saveLower = CoinCopyOfArray(colLower,numberColumns);
      double * saveUpper = CoinCopyOfArray(colUpper,numberColumns);
      CoinSort_2(contribution,contribution+numberRows,whichRow);
      // Set do nothing solution
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	if(whichSet[iColumn]>=numberRows)
	  newSolver->setColLower(iColumn,1.0);
      }
      newSolver->resolve();
      int nChunk = (nSOS+9)/10;
      int nPass=0;
      int inChunk=0;
      for (int i=0;i<nSOS;i++) {
	chunk[whichRow[i]]=nPass;
	inChunk++;
	if (inChunk==nChunk) {
	  inChunk=0;
	  // last two together
	  if (i+nChunk<nSOS)
	    nPass++;
	}
      }
      // adjust
      nPass++;
      for (int iPass=0;iPass<nPass;iPass++) {
	// fix last chunk and unfix this chunk
	for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	  int iSOS = whichSet[iColumn];
	  if (iSOS>=0) {
	    if (iSOS>=numberRows)
	      iSOS-=numberRows;
	    if (chunk[iSOS]==iPass-1&&betterSolution[iColumn]>0.9999) {
	      newSolver->setColLower(iColumn,1.0);
	    } else if (chunk[iSOS]==iPass) {
	      newSolver->setColLower(iColumn,saveLower[iColumn]);
	      newSolver->setColUpper(iColumn,saveUpper[iColumn]);
	    }
	  }
	}
	// solve
        returnCode = smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue,
                                         model_->getCutoff(), "CbcHeuristicRENS");
        if (returnCode < 0) {
            returnCode = 0; // returned on size
	    break;
	} else if ((returnCode&1)==0) {
	  // no good
	  break;
	}
      }
      if ((returnCode&2) != 0) {
	// could add cut
	returnCode &= ~2;
      }
      delete [] chunk;
      delete [] saveLower;
      delete [] saveUpper;
      delete [] whichRow;
      delete [] whichSet;
      delete [] contribution;
      delete newSolver;
      return returnCode;
    }
    
    double primalTolerance;
    solver->getDblParam(OsiPrimalTolerance, primalTolerance);

    int i;
    int numberTightened = 0;
    int numberAtBound = 0;
    int numberContinuous = numberColumns - numberIntegers;

    for (i = 0; i < numberIntegers; i++) {
        int iColumn = integerVariable[i];
        double value = currentSolution[iColumn];
        double lower = colLower[iColumn];
        double upper = colUpper[iColumn];
        value = CoinMax(value, lower);
        value = CoinMin(value, upper);
	double djValue=dj[iColumn]*direction;
#define RENS_FIX_ONLY_LOWER
#ifndef RENS_FIX_ONLY_LOWER
        if (fabs(value - floor(value + 0.5)) < 1.0e-8) {
            value = floor(value + 0.5);
            if (value == lower || value == upper)
                numberAtBound++;
            newSolver->setColLower(iColumn, value);
            newSolver->setColUpper(iColumn, value);
            numberFixed++;
        } else if (colUpper[iColumn] - colLower[iColumn] >= 2.0) {
            numberTightened++;
            newSolver->setColLower(iColumn, floor(value));
            newSolver->setColUpper(iColumn, ceil(value));
        }
#else
        if (fabs(value - floor(value + 0.5)) < 1.0e-8 &&
                floor(value + 0.5) == lower &&
	    djValue > djTolerance ) {
	  value = floor(value + 0.5);
	  numberAtBound++;
	  newSolver->setColLower(iColumn, value);
	  newSolver->setColUpper(iColumn, value);
	  numberFixed++;
        } else if (fabs(value - floor(value + 0.5)) < 1.0e-8 &&
                floor(value + 0.5) == upper &&
		   -djValue > djTolerance && (djTolerance > 0.0||type==2)) {
	  value = floor(value + 0.5);
	  numberAtBound++;
	  newSolver->setColLower(iColumn, value);
	  newSolver->setColUpper(iColumn, value);
	  numberFixed++;
        } else if (colUpper[iColumn] - colLower[iColumn] >= 2.0 &&
		   djTolerance <0.0) {
            numberTightened++;
            if (fabs(value - floor(value + 0.5)) < 1.0e-8) {
                value = floor(value + 0.5);
                if (value < upper) {
                    newSolver->setColLower(iColumn, CoinMax(value - 1.0, lower));
                    newSolver->setColUpper(iColumn, CoinMin(value + 1.0, upper));
                } else {
                    newSolver->setColLower(iColumn, upper - 1.0);
                }
            } else {
                newSolver->setColLower(iColumn, floor(value));
                newSolver->setColUpper(iColumn, ceil(value));
            }
        }
#endif
    }
    delete [] dj;
    if (numberFixed > numberIntegers / 5) {
        if ( numberFixed < numberColumns / 5) {
#define RENS_FIX_CONTINUOUS
#ifdef RENS_FIX_CONTINUOUS
            const double * colLower = newSolver->getColLower();
            //const double * colUpper = newSolver->getColUpper();
            int nAtLb = 0;
            double sumDj = 0.0;
            const double * dj = newSolver->getReducedCost();
            for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                if (!newSolver->isInteger(iColumn)) {
                    double value = currentSolution[iColumn];
                    if (value < colLower[iColumn] + 1.0e-8) {
                        double djValue = dj[iColumn] * direction;
                        nAtLb++;
                        sumDj += djValue;
                    }
                }
            }
            if (nAtLb) {
                // fix some continuous
                double * sort = new double[nAtLb];
                int * which = new int [nAtLb];
                double threshold = CoinMax((0.01 * sumDj) / static_cast<double>(nAtLb), 1.0e-6);
                int nFix2 = 0;
                for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (!newSolver->isInteger(iColumn)) {
                        double value = currentSolution[iColumn];
                        if (value < colLower[iColumn] + 1.0e-8) {
                            double djValue = dj[iColumn] * direction;
                            if (djValue > threshold) {
                                sort[nFix2] = -djValue;
                                which[nFix2++] = iColumn;
                            }
                        }
                    }
                }
                CoinSort_2(sort, sort + nFix2, which);
                nFix2 = CoinMin(nFix2, (numberColumns - numberFixed) / 2);
                for (int i = 0; i < nFix2; i++) {
                    int iColumn = which[i];
                    newSolver->setColUpper(iColumn, colLower[iColumn]);
                }
                delete [] sort;
                delete [] which;
#ifdef CLP_INVESTIGATE2
                printf("%d integers fixed (%d tightened) (%d at bound), and %d continuous fixed at lb\n",
                       numberFixed, numberTightened, numberAtBound, nFix2);
#endif
            }
#endif
        }
#ifdef COIN_DEVELOP
        printf("%d integers fixed and %d tightened\n", numberFixed, numberTightened);
#endif
        returnCode = smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue,
                                         model_->getCutoff(), "CbcHeuristicRENS");
        if (returnCode < 0 || returnCode == 0) {
#ifdef RENS_FIX_CONTINUOUS
            if (numberContinuous > numberIntegers && numberFixed >= numberColumns / 5) {
                const double * colLower = newSolver->getColLower();
                //const double * colUpper = newSolver->getColUpper();
                int nAtLb = 0;
                double sumDj = 0.0;
                const double * dj = newSolver->getReducedCost();
                double direction = newSolver->getObjSense();
                for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (!newSolver->isInteger(iColumn)) {
                        double value = currentSolution[iColumn];
                        if (value < colLower[iColumn] + 1.0e-8) {
                            double djValue = dj[iColumn] * direction;
                            nAtLb++;
                            sumDj += djValue;
                        }
                    }
                }
                if (nAtLb) {
                    // fix some continuous
                    double * sort = new double[nAtLb];
                    int * which = new int [nAtLb];
                    double threshold = CoinMax((0.01 * sumDj) / static_cast<double>(nAtLb), 1.0e-6);
                    int nFix2 = 0;
                    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                        if (!newSolver->isInteger(iColumn)) {
                            double value = currentSolution[iColumn];
                            if (value < colLower[iColumn] + 1.0e-8) {
                                double djValue = dj[iColumn] * direction;
                                if (djValue > threshold) {
                                    sort[nFix2] = -djValue;
                                    which[nFix2++] = iColumn;
                                }
                            }
                        }
                    }
                    CoinSort_2(sort, sort + nFix2, which);
                    nFix2 = CoinMin(nFix2, (numberColumns - numberFixed) / 2);
                    for (int i = 0; i < nFix2; i++) {
                        int iColumn = which[i];
                        newSolver->setColUpper(iColumn, colLower[iColumn]);
                    }
                    delete [] sort;
                    delete [] which;
#ifdef CLP_INVESTIGATE2
                    printf("%d integers fixed (%d tightened) (%d at bound), and %d continuous fixed at lb\n",
                           numberFixed, numberTightened, numberAtBound, nFix2);
#endif
                }
                returnCode = smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue,
                                                 model_->getCutoff(), "CbcHeuristicRENS");
	    }
#endif
	    if (returnCode < 0 || returnCode == 0) {
		// Do passes fixing up those >0.9 and
		// down those < 0.05
#define RENS_PASS 3
	      //#define KEEP_GOING
#ifdef KEEP_GOING
 	        double * saveLower = CoinCopyOfArray(colLower,numberColumns);
	        double * saveUpper = CoinCopyOfArray(colUpper,numberColumns);
		bool badPass=false;
		int nSolved=0;
#endif
		for (int iPass=0;iPass<RENS_PASS;iPass++) {
		  int nFixed=0;
		  int nFixedAlready=0;
		  int nFixedContinuous=0;
		  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
		    if (colUpper[iColumn]>colLower[iColumn]) {
		      if (newSolver->isInteger(iColumn)) {
			double value = currentSolution[iColumn];
			double fixTo = floor(value+0.1);
			if (fixTo>value || value-fixTo < 0.05) {
			  // above 0.9 or below 0.05
			  nFixed++;
			  newSolver->setColLower(iColumn, fixTo);
			  newSolver->setColUpper(iColumn, fixTo);
			}
		      }
		    } else if (newSolver->isInteger(iColumn)) {
		      nFixedAlready++;
		    } else {
		      nFixedContinuous++;
		    }
		  }
#ifdef CLP_INVESTIGATE2
		  printf("%d more integers fixed (total %d) plus %d continuous\n",
			 nFixed,nFixed+nFixedAlready,nFixedContinuous);
#endif
#ifdef KEEP_GOING
		  if (nFixed) {
		    newSolver->resolve();
		    if (!newSolver->isProvenOptimal()) {
		      badPass=true;
		      break;
		    } else {
		      nSolved++;
		      memcpy(saveLower,colLower,numberColumns*sizeof(double));
		      memcpy(saveUpper,colUpper,numberColumns*sizeof(double));
		    }
		  } else {
		    break;
		  }
#else
		  if (nFixed) {
		    newSolver->resolve();
		    if (!newSolver->isProvenOptimal()) {
		      returnCode=0;
		      break;
		    }
		    returnCode = smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue,
						     model_->getCutoff(), "CbcHeuristicRENS");
		  } else {
		    returnCode=0;
		  }
		  if (returnCode>=0)
		    break;
		}
		if (returnCode < 0) 
                returnCode = 0; // returned on size
#endif
	    }
#ifdef KEEP_GOING
	    if (badPass) {
	      newSolver->setColLower(saveLower);
	      newSolver->setColUpper(saveUpper);
	      newSolver->resolve();
	    }
	    delete [] saveLower;
	    delete [] saveUpper;
	    if (nSolved)
	      returnCode = 
		smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue,
				    model_->getCutoff(), "CbcHeuristicRENS");
	    else
	      returnCode=0;
	}
#endif
        }
        //printf("return code %d",returnCode);
        if ((returnCode&2) != 0) {
            // could add cut
            returnCode &= ~2;
#ifdef COIN_DEVELOP
            if (!numberTightened && numberFixed == numberAtBound)
                printf("could add cut with %d elements\n", numberFixed);
#endif
        } else {
            //printf("\n");
        }
    }
    //delete [] whichRow;
    //delete [] contribution;
    delete newSolver;
    fractionSmall_ = saveFractionSmall;
    return returnCode;
}
// update model
void CbcHeuristicRENS::setModel(CbcModel * model)
{
    model_ = model;
}

