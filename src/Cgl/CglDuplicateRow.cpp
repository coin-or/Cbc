// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <climits>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>
//#define PRINT_DEBUG
//#define CGL_DEBUG
#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinFinite.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CglDuplicateRow.hpp"
#include "CglStored.hpp"
//-------------------------------------------------------------------
// Generate duplicate row column cuts
//------------------------------------------------------------------- 
void CglDuplicateRow::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			      const CglTreeInfo info)
{
#ifdef CGL_DEBUG
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
  if (debugger&&debugger->onOptimalPath(si)) {
    printf("On optimal path\n");
  }
#endif
  // Don't do in tree ?
  if (info.inTree) {
    // but do any stored cuts
    if (storedCuts_)
      storedCuts_->generateCuts(si,cs,info);
    return;
  }
  if ((mode_&3)!=0) {
    // bug generateCuts12(si,cs,info);
  } else if ((mode_&4)!=0) {
    generateCuts4(si,cs,info);
  } else {
    assert ((mode_&8)!=0);
    generateCuts8(si,cs,info);
  }
}
void CglDuplicateRow::generateCuts12(const OsiSolverInterface & si, OsiCuts & cs,
			      const CglTreeInfo info)
{
  int numberColumns = matrix_.getNumCols();
  CoinPackedVector ubs;
  
  // Column copy
  const double * element = matrix_.getElements();
  const int * row = matrix_.getIndices();
  const CoinBigIndex * columnStart = matrix_.getVectorStarts();
  const int * columnLength = matrix_.getVectorLengths();
  // Row copy
  const double * elementByRow = matrixByRow_.getElements();
  const int * column = matrixByRow_.getIndices();
  const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
  const int * rowLength = matrixByRow_.getVectorLengths();
  const double * columnLower = si.getColLower();
  const double * columnUpper = si.getColUpper();
  int nFree=0;
  int nOut=0;
  int nFixed=0;
  int i;
  int numberRows=matrix_.getNumRows();
  const double * rowLower = si.getRowLower();
  const double * rowUpper = si.getRowUpper();
  int * effectiveRhs = CoinCopyOfArray(rhs_,numberRows);
  int * effectiveLower = CoinCopyOfArray(lower_,numberRows);
  double * effectiveRhs2 = new double [numberRows];
  const char * intVar = si.getColType();
  /* For L or G rows - compute effective lower we have to reach */
  // mark bad rows - also used for domination
  for (i=0;i<numberRows;i++) {
    int duplicate=-1;
    int LorG=0;
    double rhs=0.0;
    if (rowLower[i]<-1.0e20) {
      LorG=1;
      rhs=rowUpper[i];
    } else if (rowUpper[i]>1.0e20) {
      LorG=2;
      rhs=rowLower[i];
    }
    CoinBigIndex j;
    for (j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
      int iColumn = column[j];
      double value=elementByRow[j];
      if (LorG==1) {
	// need lowest contribution
	if (value>0.0) {
	  double bound=columnLower[iColumn];
	  if (bound>-1.0e20) 
	    rhs -= bound*value;
	  else
	    LorG=-1;
	} else {
	  double bound=columnUpper[iColumn];
	  if (bound<1.0e20) 
	    rhs -= bound*value;
	  else
	    LorG=-1;
	}
      } else if (LorG==2) {
	// need highest contribution
	if (value<0.0) {
	  double bound=columnLower[iColumn];
	  if (bound>-1.0e20) 
	    rhs -= bound*value;
	  else
	    LorG=-2;
	} else {
	  double bound=columnUpper[iColumn];
	  if (bound<1.0e20) 
	    rhs -= bound*value;
	  else
	    LorG=-2;
	}
      }
      if (duplicate!=-3) {
	if (value!=1.0) {
	  duplicate=-3;
	  rhs_[i]=-1000000;
	  //break;
	} else if (!intVar[iColumn]) {
	  duplicate=-5;
	}
      }
    }
    duplicate_[i]=duplicate;
    if (!LorG) {
      effectiveRhs2[i]=0.0;
    } else if (LorG<0) {
      // weak
      effectiveRhs2[i]=-COIN_DBL_MAX;
    } else if (LorG==1) {
      effectiveRhs2[i]=rhs;
    } else {
      effectiveRhs2[i]=rhs;
    }
  }
  double * colUpper2 = CoinCopyOfArray(columnUpper,numberColumns);
  if (!info.pass&&(mode_&2)!=0) {
    // First look at duplicate or dominated columns
    double * random = new double[numberRows];
    double * sort = new double[numberColumns+1];
    if (info.randomNumberGenerator) {
      const CoinThreadRandom * randomGenerator = info.randomNumberGenerator;
      for (i=0;i<numberRows;i++) {
	if (rowLower[i]<-1.0e20||rowUpper[i]>1.0e20)
	  random[i]=0.0;
	else
	  random[i] = randomGenerator->randomDouble();
      }
    } else {
      for (i=0;i<numberRows;i++) {
	if (rowLower[i]<-1.0e20||rowUpper[i]>1.0e20)
	  random[i]=0.0;
	else
	  random[i] = CoinDrand48();
      }
    }
    int * which = new int[numberColumns];
    int nPossible=0;
    const char * intVar = si.getColType();
    for ( i=0;i<numberColumns;i++) {
      if (intVar[i]==1) {
	double value = 0.0;
	for (CoinBigIndex jj=columnStart[i];jj<columnStart[i]+columnLength[i];jj++) {
	  int iRow = row[jj];
	  value += element[jj]*random[iRow];
	}
	sort[nPossible]=value;
	which[nPossible++]=i;
      }
    }
    sort[nPossible]=COIN_DBL_MAX;
    CoinSort_2(sort,sort+nPossible,which);
    int last=maximumDominated_-1;
    double lastValue=-1.0;
    const double *objective = si.getObjCoefficients() ;
    double direction = si.getObjSense();
    // arrays for checking
    double * elementEqualJ = new double [2*numberRows];
    CoinZeroN(elementEqualJ,numberRows); // for memory checkers
    double * elementGeJ = elementEqualJ + numberRows;
    CoinZeroN(elementGeJ,numberRows);
    int * rowEqualJ = new int[2*numberRows];
    CoinZeroN(rowEqualJ,numberRows); // for memory checkers
    int * rowGeJ = rowEqualJ + numberRows;
    char * mark = new char[numberRows];
    CoinZeroN(mark,numberRows);
#if 1
    for (i=0;i<nPossible+1;i++) {
      if (sort[i]>lastValue) {
	if (i-last<=maximumDominated_&&i>last+1) {
	  // look to see if dominated
	  for (int j=last;j<i;j++) {
	    int jColumn = which[j];
	    // skip if already fixed
	    if (!colUpper2[jColumn]||columnLower[jColumn])
	      continue;
	    int nGeJ=0;
	    int nEqualJ=0;
	    CoinBigIndex jj;
	    int nJ=columnLength[jColumn];
	    for (jj=columnStart[jColumn];jj<columnStart[jColumn]+columnLength[jColumn];jj++) {
	      int iRow = row[jj];
	      if (random[iRow]) {
		elementEqualJ[nEqualJ]=element[jj];
		rowEqualJ[nEqualJ++]=iRow;
	      } else {
		// swap sign so all rows look like G
		elementGeJ[iRow]=(rowUpper[iRow]>1.0e20) ? element[jj] : -element[jj];
		rowGeJ[nGeJ++]=iRow;
	      }
	    }
	    double objValueJ = objective[jColumn]*direction;
	    for (int k=j+1;k<i;k++) {
	      int kColumn = which[k];
	      // skip if already fixed
	      if (!colUpper2[kColumn]||columnLower[kColumn])
		continue;
	      int nK=columnLength[kColumn];
	      double objValueK = objective[kColumn]*direction;
	      if ((nJ-nK)*(objValueK-objValueJ)>0.0)
		continue;
	      int nEqualK=0;
	      // -2 no good, -1 J dominates K, 0 unknown or equal, 1 K dominates J
	      int dominate=0;
	      // mark
	      CoinBigIndex kk;
	      for (kk=0;kk<nGeJ;kk++)
		mark[rowGeJ[kk]]=1;
	      for (kk=columnStart[kColumn];kk<columnStart[kColumn]+columnLength[kColumn];kk++) {
		int iRow = row[kk];
		if (random[iRow]) {
		  if (iRow!=rowEqualJ[nEqualK]||
		      element[kk]!=elementEqualJ[nEqualK]) {
		    dominate=-2;
		    break;
		  } else {
		    nEqualK++;
		  }
		} else {
		  // swap sign so all rows look like G
		  double valueK = (rowUpper[iRow]>1.0e20) ? element[kk] : -element[kk];
		  double valueJ = elementGeJ[iRow];
		  mark[iRow]=0;
		  if (valueJ==valueK) {
		    // equal
		  } else if (valueJ>valueK) {
		    // J would dominate K
		    if (dominate==1) {
		      // no good
		      dominate=-2;
		      break;
		    } else {
		      dominate=-1;
		    }
		  } else {
		    // K would dominate J
		    if (dominate==-1) {
		      // no good
		      dominate=-2;
		      break;
		    } else {
		      dominate=1;
		    }
		  }
		}
	      }
	      kk=0;
	      if (dominate!=-2) {
		// unmark and check
		for (;kk<nGeJ;kk++) {
		  int iRow = rowGeJ[kk];
		  if (mark[iRow]) {
		    double valueK = 0.0;
		    double valueJ = elementGeJ[iRow];
		    if (valueJ>valueK) {
		      // J would dominate K
		      if (dominate==1) {
			// no good
			dominate=-2;
			break;
		      } else {
			dominate=-1;
		      }
		    } else {
		      // K would dominate J
		      if (dominate==-1) {
			// no good
			dominate=-2;
			break;
		      } else {
			dominate=1;
		      }
		    }
		  }
		  mark[iRow]=0;
		}
	      } 
	      // just unmark rest
	      for (;kk<nGeJ;kk++)
		mark[rowGeJ[kk]]=0;
	      if (nEqualK==nEqualJ&&dominate!=-2) {
		if (objValueJ==objValueK) {
		  if (dominate<=0) {
		    // say J dominates
		    assert (colUpper2[kColumn]);
		    dominate=-1;
		  } else {
		    // say K dominates
		    assert (colUpper2[jColumn]);
		    dominate=1;
		  }
		} else if (objValueJ<objValueK&&dominate<=0) {
		  // say J dominates
		  assert (colUpper2[kColumn]);
		  dominate=-1;
		} else if (objValueJ>objValueK&&dominate==1) {
		  // say K dominates
		  assert (colUpper2[jColumn]);
		  dominate=1;
		} else {
		  dominate=0;
		}
		if (dominate) {
		  // see if both can be 1
		  bool canFix=false;
		  for (int jj=0;jj<nEqualJ;jj++) {
		    double value = 2.0*elementEqualJ[jj];
		    int iRow = rowEqualJ[jj];
		    if (duplicate_[iRow]==-1&&rowUpper[iRow]<1.999999) {
		      canFix=true;
		    } else {
		      double minSum=0.0;
		      double maxSum=0.0;
		      for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
			int iColumn = column[j];
			if (iColumn!=jColumn&&iColumn!=kColumn) {
			  double elValue = elementByRow[j];
			  double lo = columnLower[iColumn];
			  double up = colUpper2[iColumn];
			  if (elValue>0.0) {
			    minSum += lo*elValue;
			    maxSum += up*elValue;
			  } else {
			    maxSum += lo*elValue;
			    minSum += up*elValue;
			  }
			}
		      }
		      if (minSum+value>rowUpper[iRow]+1.0e-5)
			canFix=true;
		      else if (maxSum+value<rowLower[iRow]-1.0e-5)
			canFix=true;
		    }
		    if (canFix)
		      break;
		  }
		  if (!canFix) {
		    for (kk=columnStart[kColumn];kk<columnStart[kColumn]+columnLength[kColumn];kk++) {
		      int iRow = row[kk];
		      if (!random[iRow]) {
			if (rowUpper[iRow]<1.0e20) {
			  // just <= row
			  double valueK = element[kk] - elementGeJ[iRow];
			  if (valueK>effectiveRhs2[iRow]+1.0e-4) {
			    canFix=true;
			    break;
			  }
			} else {
			  // >= row
			  double valueK = element[kk] + elementGeJ[iRow];
			  if (valueK<effectiveRhs2[iRow]-1.0e-4) {
			    canFix=true;
			    break;
			  }
			}
		      }
		    }
		  }
		  if (canFix) {
		    int iColumn = (dominate>0) ? jColumn : kColumn;
		    nFixed++;
		    assert (!columnLower[iColumn]);
		    colUpper2[iColumn]=0.0;
		    ubs.insert(iColumn,0.0);
		    if (iColumn==jColumn)
		      break; // no need to carry on on jColumn
		  } else {
		    int iDominated = (dominate>0) ? jColumn : kColumn;
		    int iDominating = (dominate<0) ? jColumn : kColumn;
		    double els[]={1.0,-1.0};
		    int inds[2];
		    inds[0]=iDominating;
		    inds[1]=iDominated;
		    if (!storedCuts_)
		      storedCuts_ = new CglStored();
		    storedCuts_->addCut(0.0,COIN_DBL_MAX,2,inds,els);
		  }
		}
	      }
	    }
	    for (jj=0;jj<nGeJ;jj++) {
	      int iRow = rowGeJ[jj];
	      elementGeJ[iRow]=0.0;
	    }
	  }
	}
	last=i;
	lastValue = sort[i];
      }
    }
#endif
    delete [] mark;
    delete [] elementEqualJ;
    delete [] rowEqualJ;
    delete [] random;
    delete [] sort;
    delete [] which;
#ifdef COIN_DEVELOP
    int numberCuts = storedCuts_ ? storedCuts_->sizeRowCuts() : 0;
    if (nFixed||numberCuts) 
      printf("** %d fixed and %d cuts from domination\n",nFixed,numberCuts);
#endif
  }
  delete [] effectiveRhs2;
  bool infeasible=false;
  // if we were just doing columns - mark all as bad
  if ((mode_&1)==0) {
    for (i=0;i<numberRows;i++) {
      duplicate_[i]=-3;
      rhs_[i]=-1000000;
      effectiveLower[i]=-1000000;
    }
  }
  for ( i=0;i<numberColumns;i++) {
    if (columnLower[i]) {
      double value = columnLower[i];
      for (CoinBigIndex jj=columnStart[i];jj<columnStart[i]+columnLength[i];jj++) {
        int iRow = row[jj];
        nOut += static_cast<int> (element[jj]*value);
        effectiveRhs[iRow] -= static_cast<int> (element[jj]*value);
        effectiveLower[iRow] -= static_cast<int> (element[jj]*value);
      }
    }
  }
  for ( i=0;i<numberColumns;i++) {
    if (columnLower[i]!=colUpper2[i]) {
      bool fixed=false;
      for (CoinBigIndex jj=columnStart[i];jj<columnStart[i]+columnLength[i];jj++) {
        int iRow = row[jj];
        if (rhs_[iRow]>=0&&element[jj]>effectiveRhs[iRow]) 
          fixed=true;
      }
      if (fixed) {
        nFixed++;
        colUpper2[i]=columnLower[i];
        ubs.insert(i,columnLower[i]);
      } else {
        nFree++;
      }
    }
  }
  // See if anything odd
  char * check = new char[numberColumns];
  memset(check,0,numberColumns);
  int * which2 = new int[numberColumns];
  for (i=0;i<numberRows;i++) {
    if (duplicate_[i]==-5) {
      if ((rowLower[i]<=0.0||rowLower[i]==rowUpper[i])&&
	  rowUpper[i]==floor(rowUpper[i])) {
	effectiveRhs[i]= static_cast<int> (rowUpper[i]);
	effectiveLower[i] = static_cast<int> (std::max(0.0,rowLower[i]));
	bool goodRow=true;
	for (CoinBigIndex j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	  int iColumn = column[j];
	  double value=columnLower[iColumn];
	  if (value) {
	    if (value==floor(value)) {
	      effectiveRhs[i] -= static_cast<int> (value);
	      effectiveLower[i] -= static_cast<int> (value);
	    } else {
	      goodRow=false;
	    }
	  }
	}
	if (goodRow)
	  duplicate_[i] = -1; // can have continuous variables now
	else
	  duplicate_[i] = -3;
      } else {
	duplicate_[i] = -3;
      }
    }
    if (duplicate_[i]==-1) {
      if (effectiveRhs[i]>0) {
	// leave
      } else if (effectiveRhs[i]==0) {
        duplicate_[i]=-2;
      } else {
        duplicate_[i]=-3;
	// leave unless >=1 row
	if (effectiveLower[i]==1&&rhs_[i]<0.0)
	  duplicate_[i]=-4;
      }
    } else {
      effectiveRhs[i]=-1000;
    }
  }
  // Look at <= rows
  double maxLook=100*numberRows;
  double nLook=maxLook;
  for (i=0;i<numberRows;i++) {
    // initially just one
    if (effectiveRhs[i]==1&&duplicate_[i]==-1) {
      int nn=0;
      CoinBigIndex j;
      int k;
      nLook -= numberRows-i;
      if (nLook<0)
	break;
      for (j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	int iColumn = column[j];
	if (columnLower[iColumn]!=colUpper2[iColumn]) {
#ifndef NDEBUG
          assert (elementByRow[j]==1.0);
#endif
          check[iColumn]=1;
          which2[nn++]=iColumn;
        }
      }
      for ( k=i+1;k<numberRows;k++) {
        if (effectiveRhs[k]==1&&duplicate_[k]==-1) {
          int nn2=0;
          int nnsame=0;
          for ( j=rowStart[k];j<rowStart[k]+rowLength[k];j++) {
            int iColumn = column[j];
            if (columnLower[iColumn]!=colUpper2[iColumn]) {
#ifndef NDEBUG
              assert (elementByRow[j]==1.0);
#endif
              nn2++;
              if (check[iColumn]) 
                nnsame++;
            }
          }
	  //if (nnsame)
	  //printf("rows %d and %d, %d same - %d %d\n",
	  //   i,k,nnsame,nn,nn2);
	  bool checked=false;
          if (nnsame==nn2) {
            if (nn2<nn&&effectiveLower[k]==rhs_[k]&&rhs_[i]==rhs_[k]) {
              if (logLevel_)
                printf("row %d strict subset of row %d, fix some in row %d\n",
                       k,i,i);
              // treat i as duplicate
              duplicate_[i]=k;
	      nLook=maxLook; // reset
              // zero out check so we can see what is extra
              for ( j=rowStart[k];j<rowStart[k]+rowLength[k];j++) {
                int iColumn = column[j];
                check[iColumn]=0; 
              }
              // now redo and fix
              nn=0;
              for (j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
                int iColumn = column[j];
                if (columnLower[iColumn]!=colUpper2[iColumn]) {
                  if (check[iColumn]) {
                    // fix
                    colUpper2[iColumn]=columnLower[iColumn];
                    nFixed++;
                    ubs.insert(iColumn,columnLower[iColumn]);
                    check[iColumn]=0;
                  } else {
                    check[iColumn]=1;
                    which2[nn++]=iColumn;
                  }
                }
              }
	      checked=true;
            } else if (nn2==nn&&effectiveLower[i]==rhs_[i]&&effectiveLower[k]==rhs_[k]) {
              if (logLevel_)
                printf("row %d identical to row %d\n",
                       k,i);
              duplicate_[k]=i;
	      nLook=maxLook; // reset
	      checked=true;
            } else if (nn2>=nn&&effectiveLower[i]==rhs_[i]&&effectiveLower[k]==rhs_[k]) {
              abort();
            }
          } else if (nnsame==nn&&nn2>nn&&effectiveLower[i]==rhs_[i]&&rhs_[i]<=rhs_[k]) {
            if (logLevel_)
              printf("row %d strict superset of row %d, fix some in row %d\n",
                     k,i,k);
            // treat k as duplicate
            duplicate_[k]=i;
	    nLook=maxLook; // reset
            // set check for k
            for ( j=rowStart[k];j<rowStart[k]+rowLength[k];j++) {
              int iColumn = column[j];
              if (columnLower[iColumn]!=colUpper2[iColumn]) 
                check[iColumn]=1; 
            }
            // zero out check so we can see what is extra
            for ( j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
              int iColumn = column[j];
              check[iColumn]=0; 
            }
            //  fix
            for (j=rowStart[k];j<rowStart[k]+rowLength[k];j++) {
              int iColumn = column[j];
              if (check[iColumn]) {
                // fix
                colUpper2[iColumn]=columnLower[iColumn];
                nFixed++;
                ubs.insert(iColumn,columnLower[iColumn]);
                check[iColumn]=0;
              }
            }
            // redo
            nn=0;
            for (j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
              int iColumn = column[j];
              if (columnLower[iColumn]!=colUpper2[iColumn]) {
                check[iColumn]=1;
                which2[nn++]=iColumn;
              }
            }
	    checked=true;
          }
	  if (!checked) {
            // may be redundant
            if (nnsame==nn2) {
              // k redundant ?
              if (nn2<nn&&effectiveLower[k]<=0&&rhs_[i]<=rhs_[k]) {
                if (logLevel_)
                  printf("row %d slack subset of row %d, drop row %d\n",
                         k,i,k);
                // treat k as duplicate
                duplicate_[k]=i;
		nLook=maxLook; // reset
              }
            } else if (nnsame==nn) {
              // i redundant ?
              if (nn2>nn&&effectiveLower[i]<=0&&rhs_[k]<=rhs_[i]) {
                if (logLevel_)
                  printf("row %d slack subset of row %d, drop row %d\n",
                         i,k,i);
                // treat i as duplicate
                duplicate_[i]=k;
		nLook=maxLook; // reset
              }
            }
          }
        }
      }
      for (k=0;k<nn;k++) 
        check[which2[k]]=0;
      
    }
  }
  // Look at >=1 rows
  for (i=0;i<numberRows;i++) {
    if (duplicate_[i]==-4) {
      int nn=0;
      CoinBigIndex j;
      int k;
      for (j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	int iColumn = column[j];
	if (columnLower[iColumn]!=colUpper2[iColumn]) {
#ifndef NDEBUG
          assert (elementByRow[j]==1.0);
#endif
          check[iColumn]=1;
          which2[nn++]=iColumn;
        }
      }
      for ( k=i+1;k<numberRows;k++) {
        if (duplicate_[k]==-4) {
          int nn2=0;
          int nnsame=0;
          for ( j=rowStart[k];j<rowStart[k]+rowLength[k];j++) {
            int iColumn = column[j];
            if (columnLower[iColumn]!=colUpper2[iColumn]) {
#ifndef NDEBUG
              assert (elementByRow[j]==1.0);
#endif
              nn2++;
              if (check[iColumn]) 
                nnsame++;
            }
          }
	  // may be redundant
	  if (nnsame==nn||nnsame==nn2) {
	    if (nn2>nn) {
	      // k redundant
	      if (logLevel_) 
		printf("row %d slack superset of row %d, drop row %d\n",
		       k,i,k);
	      // treat k as duplicate
	      duplicate_[k]=i;
	    } else if (nn2<nn) {
	      // i redundant ?
	      if (logLevel_)
		printf("row %d slack superset of row %d, drop row %d\n",
		     i,k,i);
	      // treat i as duplicate
	      duplicate_[i]=k;
	    } else {
	      if (logLevel_) 
		printf("row %d same as row %d, drop row %d\n",
		       k,i,k);
	      // treat k as duplicate
	      duplicate_[k]=i;
	    }
          }
        }
      }
      for (k=0;k<nn;k++) 
        check[which2[k]]=0;
      
    }
  }
  if ((mode_&1)!=0&&true) {
    // look at doubletons
    const double * rowLower = si.getRowLower();
    const double * rowUpper = si.getRowUpper();
    int i;
    int nPossible=0;
    const char * intVar = si.getColType();
    for (i=0;i<numberRows;i++) {
      if (rowLength[i]==2&&(duplicate_[i]<0&&duplicate_[i]!=-2)) {
	bool possible=true;
	CoinBigIndex j;
	for (j=rowStart[i];j<rowStart[i]+2;j++) {
	  int iColumn = column[j];
	  if (fabs(elementByRow[j])!=1.0||intVar[iColumn]!=1) {
	    possible=false;
	    break;
	  }
	}
	if (possible) {
	  CoinBigIndex j = rowStart[i];
	  int column0 = column[j];
	  double element0 = elementByRow[j];
	  int column1 = column[j+1];
	  double element1 = elementByRow[j+1];
	  if (element0==1.0&&element1==1.0&&rowLower[i]==1.0&&
	      rowUpper[i]>1.0e30) {
	    if (logLevel_) {
	      printf("Cover row %d %g <= ",i,rowLower[i]);
	      printf("(%d,%g) (%d,%g) ",column0,element0,column1,element1);
	      printf(" <= %g\n",rowUpper[i]);
	    }
	    effectiveRhs[nPossible++]=i;
	  } else {
	    // not at present
	    //printf("NON Cover row %d %g <= ",i,rowLower[i]);
	    //printf("(%d,%g) (%d,%g) ",column0,element0,column1,element1);
	    //printf(" <= %g\n",rowUpper[i]);
	  }
	}
      }
    }
    if (nPossible) {
      int * check2 = new int [numberColumns];
      CoinFillN(check2,numberColumns,-1);
      for (int iPossible=0;iPossible<nPossible;iPossible++) {
#ifndef NDEBUG
	for (i=0;i<numberColumns;i++)
	  assert (check2[i]==-1);
#endif
	i=effectiveRhs[iPossible];
	CoinBigIndex j = rowStart[i];
	int column0 = column[j];
	int column1 = column[j+1];
	int k;
	int nMarked=0;
	for (int kPossible=iPossible+1;kPossible<nPossible;kPossible++) {
	  k=effectiveRhs[kPossible];
	  CoinBigIndex j = rowStart[k];
	  int columnB0 = column[j];
	  int columnB1 = column[j+1];
	  if (column0==columnB1||column1==columnB1) {
	    columnB1=columnB0;
	    columnB0=column[j+1];
	  }
	  bool good = false;
	  if (column0==columnB0) {
	    if (column1==columnB1) {
	      // probably should have been picked up
	      // safest to ignore
	    } else {
	      good=true;
	    }
	  } else if (column1==columnB0) {
	    if (column0==columnB1) {
	      // probably should have been picked up
	      // safest to ignore
	    } else {
	      good=true;
	    }
	  }
	  if (good) {
	    if (check2[columnB1]<0) {
	      check2[columnB1]=k;
	      which2[nMarked++]=columnB1;
	    } else {
	      // found
#ifndef COIN_DEVELOP
	      if (logLevel_>1)
#endif 
		printf("***Make %d %d %d >=2 and take out rows %d %d %d\n",
		       columnB1,column0,column1,
		       i,k,check2[columnB1]);
	      OsiRowCut rc;
	      rc.setLb(2.0);
	      rc.setUb(COIN_DBL_MAX);   
	      int index[3];
	      double element[3]={1.0,1.0,1.0};
	      index[0]=column0;
	      index[1]=column1;
	      index[2]=columnB1;
	      rc.setRow(3,index,element,false);
	      cs.insertIfNotDuplicate(rc);
	      // drop rows
	      duplicate_[i]=-2;
	      duplicate_[k]=-2;
	      duplicate_[check2[columnB1]]=-2;
	    }
	  }
	}
	for (k=0;k<nMarked;k++) {
	  int iColumn = which2[k];
	  check2[iColumn]=-1;
	}
      }
      delete [] check2;
    }
  }
  delete [] check;
  delete [] which2;
  delete [] colUpper2;
  int nRow=0;
  sizeDynamic_=1;
  for (i=0;i<numberRows;i++) {
    if (duplicate_[i]!=-3) {
      if (duplicate_[i]==-1) {
        nRow++;
        int k=effectiveRhs[i];
        while (k) {
          if (sizeDynamic_<1000000000)
            sizeDynamic_ = sizeDynamic_<<1;
          k = k >>1;
        }
      }
    } else {
      duplicate_[i]=-1;
    }
  }
  delete [] effectiveRhs;
  delete [] effectiveLower;

  if (logLevel_)
    printf("%d free (but %d fixed this time), %d out of rhs, DP size %d, %d rows\n",
           nFree,nFixed,nOut,sizeDynamic_,nRow);
  if (nFixed) {
    OsiColCut cc;
    cc.setUbs(ubs);
    cc.setEffectiveness(100.0);
    cs.insert(cc);
  }
  if (infeasible) {
    // generate infeasible cut and return
    OsiRowCut rc;
    rc.setLb(COIN_DBL_MAX);
    rc.setUb(0.0);   
    cs.insertIfNotDuplicate(rc);
  }
}
void CglDuplicateRow::generateCuts4(const OsiSolverInterface & si, OsiCuts & cs,
			      const CglTreeInfo )
{
  int numberColumns = matrix_.getNumCols();
  
  // Column copy
  const double * element = matrix_.getElements();
  const int * row = matrix_.getIndices();
  const CoinBigIndex * columnStart = matrix_.getVectorStarts();
  const int * columnLength = matrix_.getVectorLengths();
  const double * columnLower = si.getColLower();
  const double * columnUpper = si.getColUpper();
  int nFixed=0;
  int numberRows=matrix_.getNumRows();
  const double * rowLower = si.getRowLower();
  const double * rowUpper = si.getRowUpper();
  bool infeasible=false;
  // try more complicated domination
  int * originalColumns = new int [numberColumns];
  int * originalRows = new int [2*numberRows];
  int * rowCount = originalRows+numberRows;
  memset(rowCount,0,numberRows*sizeof(int));
  memset(originalRows,0,numberRows*sizeof(int));
  unsigned char * rowFlag = new unsigned char[numberRows];
  unsigned char * columnFlag = new unsigned char[2*numberColumns];
  double * newBound = new double[numberColumns];
  double * trueLower = new double[numberColumns];
  double * effectiveRhs = new double [2*numberRows];
  double *rhs2 = effectiveRhs+numberRows;
  // first take out fixed stuff
  memset(rhs2,0,numberRows*sizeof(double));
  memset(rowFlag,0,numberRows);
  memset(columnFlag,0,numberColumns);
  const char * intVar = si.getColType();
  int nCol2=0;
  for (int i=0;i<numberColumns;i++) {
    if (columnLower[i]<-1.0e20&&columnUpper[i]>-1.0e20) {
      for (CoinBigIndex jj=columnStart[i];jj<columnStart[i]+columnLength[i];jj++) {
	int iRow = row[jj];
	rowFlag[iRow] |= 8; // say no good
      }
    } else {
      double lo=columnLower[i];
      double up=columnUpper[i];
      double value=lo;
      if (lo<-1.0e20) {
	columnFlag[nCol2]=1; // say flipped
	lo=-up;
	up=-value;
	value=-lo;
	abort(); // double check
      }
      if (up>lo) {
	int add=0;
	if (intVar[i]) {
	  columnFlag[nCol2]|=2;
	  if (up>lo+1.5) {
	    add=1; // only allow one general
	    columnFlag[nCol2]|=4;
	  }
	} else {
	  add=1;
	}
	newBound[nCol2]=up-lo;
	trueLower[nCol2]=lo;
	originalColumns[nCol2++]=i;
	for (CoinBigIndex jj=columnStart[i];jj<columnStart[i]+columnLength[i];jj++) {
	  int iRow = row[jj];
	  rowCount[iRow]++;
	  originalRows[iRow] += add;
	}
      }
      for (CoinBigIndex jj=columnStart[i];jj<columnStart[i]+columnLength[i];jj++) {
	int iRow = row[jj];
	rhs2[iRow] += element[jj]*value;
      }
    }
  }
  int nRow2=0;
  for (int i=0;i<numberRows;i++) {
    int nCont=originalRows[i];
    int nInt=rowCount[i]-nCont;
    unsigned char flag=rowFlag[i];
    if (nCont>1||!nInt)
      flag=8; // don't look at for now
    if (rowLower[i]==rowUpper[i]) 
      flag |= 1;
    else if (rowLower[i]>-1.0e20&&rowUpper[i]<1.0e20)
      flag |=8;
    else if (rowUpper[i]>1.0e20)
      flag |= 2;
    if ((flag&8)==0) {
      rowCount[nRow2]=rowCount[i];
      originalRows[nRow2++]=i;
    }
  }
  CoinSort_2(rowCount,rowCount+nRow2,originalRows);
  for (int i=0;i<nRow2;i++) {
    int k=originalRows[i];
    unsigned char flag=0;
    if (rowLower[k]==rowUpper[k]) 
      flag |= 1;
    else if (rowUpper[k]>1.0e20)
      flag |= 2;
    rowFlag[i]=flag;
  }
  if (nRow2&&nCol2) {
    CoinPackedMatrix small(matrix_,nRow2,originalRows,
			   nCol2,originalColumns);
    // Column copy
    small.removeGaps();
    double * element = small.getMutableElements();
    const int * row = small.getIndices();
    const CoinBigIndex * columnStart = small.getVectorStarts();
    //const int * columnLength = small.getVectorLengths();
    for (int i=0;i<nCol2;i++) {
      for (CoinBigIndex jj=columnStart[i];jj<columnStart[i+1];jj++) {
	int iRow=row[jj];
	if ((rowFlag[iRow]&2)!=0)
	  element[jj] = -element[jj];
      }
      if ((columnFlag[i]&1)!=0) {
	for (CoinBigIndex jj=columnStart[i];jj<columnStart[i+1];jj++) {
	  element[jj] = -element[jj];
	}
      }
    }
    // same thing for row matrix
    CoinPackedMatrix smallRow;
    smallRow.setExtraGap(0.0);
    smallRow.setExtraMajor(0.0);
    smallRow.reverseOrderedCopyOf(small);
    // Row copy
    double * elementByRow = smallRow.getMutableElements();
    int * column = smallRow.getMutableIndices();
    const CoinBigIndex * rowStart = smallRow.getVectorStarts();
    //const int * rowLength = smallRow.getVectorLengths();
    for (int i=0;i<nRow2;i++) {
      CoinBigIndex start=rowStart[i];
      CoinBigIndex end=rowStart[i+1];
      CoinSort_2(column+start,column+end,elementByRow+start);
      int k=originalRows[i];
      double rhs;
      if ((rowFlag[i]&2)==0) 
	rhs = rowUpper[k]-rhs2[k];
      else
	rhs = - (rowLower[k]-rhs2[k]);
      effectiveRhs[i]=rhs;
    }
    int nRowLook=0;
    int nRowStart=-1;
#define MAX_IN_BASE 3
#define MAX_IN_COMP 3
    for (nRowLook=0;nRowLook<nRow2;nRowLook++) {
      CoinBigIndex start1 = rowStart[nRowLook];
      int n=static_cast<int>(rowStart[nRowLook+1]-start1);
      if (n>=2&&nRowStart<0)
	nRowStart=nRowLook;
      if (n>MAX_IN_BASE) 
	break;
    }
    // cut back nRow2
    int nnRow2=nRowLook;
    for (nnRow2=0;nnRow2<nRow2;nnRow2++) {
      CoinBigIndex start1 = rowStart[nnRow2];
      int n=static_cast<int>(rowStart[nnRow2+1]-start1);
      if (n>MAX_IN_COMP) 
	break;
    }
    nRow2=nnRow2;
    nRowStart=std::max(0,nRowStart);
    unsigned char * mark = columnFlag+nCol2;
    memset(mark,0,nCol2);
    /* at most 3 0-1 integers -
       if all 0-1 then see if same allowed
       if one other then get bounds
    */
    double loC0[4];
    double upC0[4];
    int allowed0[8];
    double loC1[4];
    double upC1[4];
    int allowed1[8];
    for (int i=nRowStart;i<nRowLook;i++) {
      CoinBigIndex start0 = rowStart[i];
      int n=static_cast<int>(rowStart[i+1]-start0);
      const int * column0 = column+start0;
      const double * element0 = elementByRow+start0;
      int nInt=0;
      for (int j=0;j<n;j++) {
	int iColumn = column0[j];
	if ((columnFlag[iColumn]&(2+4))==2)
	  nInt++;
	mark[iColumn] =1;
      }
      for (int k=i+1;k<nRow2;k++) {
	if (duplicate_[k]==-2)
	  continue;
	if (duplicate_[i]==-2)
	  break;
	CoinBigIndex start1 = rowStart[k];
	int n1=static_cast<int>(rowStart[k+1]-start1);
	const int * column1 = column+start1;
	const double * element1 = elementByRow+start1;
	int nMatch=0;
	for (int j=0;j<n1;j++) {
	  if (mark[column1[j]])
	    nMatch++;
	}
	if (nMatch==n) {
#define CGL_INVESTIGATE
	  if (n==n1) {
	    // same - look at all 0-1 integers
	    if (nInt==n) {
	      // crude - should go stack based
	      if (nInt==2) {
		double upRhs = effectiveRhs[i];
		double loRhs = ((rowFlag[i]&1)!=0) ? effectiveRhs[i] : -1.0e30;
		double tolerance = std::max(1.0e-5,fabs(upRhs)*1.0e-10);
		for (int j0=0;j0<2;j0++) {
		  for (int j1=0;j1<2;j1++) {
		    double value = element0[0]*j0+element0[1]*j1;
		    int put = j0+2*j1;
		    if (value<upRhs+tolerance&&value>loRhs-tolerance)
		      allowed0[put]=1;
		    else
		      allowed0[put]=0;
		  }
		}
		upRhs = effectiveRhs[k];
		loRhs = ((rowFlag[k]&1)!=0) ? effectiveRhs[k] : -1.0e30;
		tolerance = std::max(1.0e-5,fabs(upRhs)*1.0e-10);
		for (int j0=0;j0<2;j0++) {
		  for (int j1=0;j1<2;j1++) {
		    double value = element1[0]*j0+element1[1]*j1;
		    int put = j0+2*j1;
		    if (value<upRhs+tolerance&&value>loRhs-tolerance)
		      allowed1[put]=1;
		    else
		      allowed1[put]=0;
		  }
		}
		/* interesting cases are when -
		   one forces fixing (but probably found in probing)
		   and of two is zero
		   and of two forces fixing
		   two same - this is probably only one
		*/
		int intersect[4];
		bool same=true;
		bool tighter0=true;
		bool tighter1=true;
		bool feasible=false;
		bool redundant=true;
		for (int j=0;j<4;j++) {
		  intersect[j]=allowed0[j]&allowed1[j];
		  if (intersect[j]<allowed0[j])
		    tighter0=false;
		  if (intersect[j]<allowed1[j])
		    tighter1=false;
		  if (allowed0[j]!=allowed1[j])
		    same=false;
		  if (!intersect[j])
		    redundant=false;
		  if (intersect[j])
		    feasible=true;
		}
		int fixed[2]={0,0};
		if (feasible) {
		  int count=2;
		  for (int jj=0;jj<count;jj++) {
		    int multiplier=1<<jj;
		    int increment=1<<(count-1-jj);
		    bool zeroOk=false;
		    bool oneOk=false;
		    for (int j=0;j<(1<<(count-1));j++) {
		      if (intersect[0*multiplier+increment*j])
			zeroOk=true;
		      if (intersect[1*multiplier+increment*j])
			oneOk=true;
		    }
		    if (!zeroOk) {
		      fixed[jj]=1;
		      assert (oneOk);
		    } else if (!oneOk) {
		      fixed[jj]=-1;
		    }
		  }
		}
		if (same||redundant||!feasible||fixed[0]||fixed[1]||
		    tighter0||tighter1) {
#ifdef CGL_INVESTIGATE
		  /* start debug print */
		  {
		    printf("Base %d (orig %d) ",i,originalRows[i]);
		    for (int j=0;j<n;j++) {
		      double value = element0[j];
		      if (j) {
			if(value>0.0)
			  printf(" +");
			else
			  printf(" ");
		      }
		      int iColumn=column0[j];
		      if ((columnFlag[iColumn]&2)==0)
			printf("%g*X%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		      else if ((columnFlag[iColumn]&(2+4))==2)
			printf("%g*B%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		      else
			printf("%g*I%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		    } 
		    if ((rowFlag[i]&1)!=0)
		      printf(" == ");
		    else
		      printf(" <= ");
		    printf("%g\n",effectiveRhs[i]);
		    printf("Comp %d (orig %d) ",k,originalRows[k]);
		    for (int j=0;j<n1;j++) {
		      double value = element1[j];
		      if (j) {
			if(value>0.0)
			  printf(" +");
			else
			  printf(" ");
		      }
		      int iColumn=column1[j];
		      if ((columnFlag[iColumn]&2)==0)
			printf("%g*X%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		      else if ((columnFlag[iColumn]&(2+4))==2)
			printf("%g*B%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		      else
			printf("%g*I%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		    } 
		    if ((rowFlag[k]&1)!=0)
		      printf(" == ");
		    else
		      printf(" <= ");
		    printf("%g\n",effectiveRhs[k]);
		  }
		  /* end debug print */
#endif
#ifdef CGL_INVESTIGATE
		  printf("**same %c redundant %c feasible %c tight %c,%c fixed %d,%d\n",
			 same ? 'Y' : 'N',
			 redundant ? 'Y' : 'N',
			 feasible ? 'Y' : 'N',
			 tighter0 ? 'Y' : 'N',
			 tighter1 ? 'Y' : 'N',
			 fixed[0],fixed[1]);
#endif
		  if (!feasible) {
#ifdef CGL_INVESTIGATE
		    printf("QQ infeasible\n");
#endif
		    infeasible=true;
		  } else if (fixed[0]||fixed[1]) {
#ifdef CGL_INVESTIGATE
		    printf("QQ fixed\n");
#endif
		    for (int k=0;k<2;k++) {
		      if (fixed[k]) {
			int iColumn=column0[k];
			int kColumn=originalColumns[iColumn];
#ifdef CGL_INVESTIGATE
			printf("true bounds %g %g\n",columnLower[kColumn],
			       columnUpper[kColumn]);
#endif
			double lo;
			if (fixed[k]>0) {
			  lo=1.0;
			} else {
			  lo=0.0;
			}
			if ((columnFlag[iColumn]&1)==0) {
			  columnFlag[iColumn] |= 16;
			  trueLower[iColumn] += lo;
			  newBound[iColumn] = 0.0;
			  for (CoinBigIndex jj=columnStart[iColumn];jj<columnStart[iColumn+1];
			       jj++) {
			    int iRow=row[jj];
			    effectiveRhs[iRow] -= lo*element[jj];
			  }
			} else {
			  abort();
			}
		      }
		    }
		  } else if (!same&&(tighter0||tighter1)) {
		    assert (!tighter0||!tighter1);
		    if (tighter0) {
#ifdef CGL_INVESTIGATE
		      printf("QQ discard oneT k\n");
#endif
		      duplicate_[k]=-2;
		    } else {
#ifdef CGL_INVESTIGATE
		      printf("QQ discard oneT i\n");
#endif
		      duplicate_[i]=-2;
		    }
		  } else if (redundant) {
#ifdef CGL_INVESTIGATE
		    printf("QQ discard both\n");
#endif
		    duplicate_[i]=-2;
		    duplicate_[k]=-2;
		  } else {
		    assert (same);
		    if (fabs(effectiveRhs[i]-effectiveRhs[k])<1.0e-7&&
			element0[1]==element1[1]&&
			element0[0]==element1[0]) {
#ifdef CGL_INVESTIGATE
		      printf("QQ discard identical k I2\n");
#endif
		      duplicate_[k]=-2;
		    } else {
#ifdef CGL_INVESTIGATE
		      printf("QQ Don't know what to do nintI==2 I\n");
#endif
		    }
		  }
		}
	      } else {
		assert (nInt==3);
		double upRhs = effectiveRhs[i];
		double loRhs = ((rowFlag[i]&1)!=0) ? effectiveRhs[i] : -1.0e30;
		double tolerance = std::max(1.0e-5,fabs(upRhs)*1.0e-10);
		for (int j0=0;j0<2;j0++) {
		  for (int j1=0;j1<2;j1++) {
		    for (int j2=0;j2<2;j2++) {
		      double value = element0[0]*j0+element0[1]*j1
			+element0[2]*j2;
		      int put = j0+2*j1+4*j2;
		      if (value<upRhs+tolerance&&value>loRhs-tolerance)
			allowed0[put]=1;
		      else
			allowed0[put]=0;
		    }
		  }
		}
		upRhs = effectiveRhs[k];
		loRhs = ((rowFlag[k]&1)!=0) ? effectiveRhs[k] : -1.0e30;
		tolerance = std::max(1.0e-5,fabs(upRhs)*1.0e-10);
		for (int j0=0;j0<2;j0++) {
		  for (int j1=0;j1<2;j1++) {
		    for (int j2=0;j2<2;j2++) {
		      double value = element1[0]*j0+element1[1]*j1
			+element1[2]*j2;
		      int put = j0+2*j1+4*j2;
		      if (value<upRhs+tolerance&&value>loRhs-tolerance)
			allowed1[put]=1;
		      else
			allowed1[put]=0;
		    }
		  }
		}
		/* interesting cases are when -
		   one forces fixing (but probably found in probing)
		   and of two is zero
		   and of two forces fixing
		   two same - this is probably only one
		*/
		int intersect[8];
		bool same=true;
		bool tighter0=true;
		bool tighter1=true;
		bool feasible=false;
		bool redundant=true;
		for (int j=0;j<8;j++) {
		  intersect[j]=allowed0[j]&allowed1[j];
		  if (intersect[j]<allowed0[j])
		    tighter0=false;
		  if (intersect[j]<allowed1[j])
		    tighter1=false;
		  if (allowed0[j]!=allowed1[j])
		    same=false;
		  if (!intersect[j])
		    redundant=false;
		  if (intersect[j])
		    feasible=true;
		}
		int fixed[3]={0,0,0};
		if (feasible) {
		  bool zeroOk[3]={false,false,false};
		  bool oneOk[3]={false,false,false};
		  for (int j0=0;j0<2;j0++) {
		    for (int j1=0;j1<2;j1++) {
		      for (int j2=0;j2<2;j2++) {
			int get = j0+2*j1+4*j2;
			if (intersect[get]) {
			  if (j0)
			    oneOk[0]=true;
			  else
			    zeroOk[0]=true;
			  if (j1)
			    oneOk[1]=true;
			  else
			    zeroOk[1]=true;
			  if (j2)
			    oneOk[2]=true;
			  else
			    zeroOk[2]=true;
			}
		      }
		    }
		  }
		  for (int jj=0;jj<3;jj++) {
		    if (!zeroOk[jj]) {
		      fixed[jj]=1;
		      assert (oneOk[jj]);
		    } else if (!oneOk[jj]) {
		      fixed[jj]=-1;
		    }
		  }
		}
		if (same||redundant||!feasible||fixed[0]||fixed[1]||fixed[2]||
		    tighter0||tighter1) {
#ifdef CGL_INVESTIGATE
		  /* start debug print */
		  {
		    printf("Base %d (orig %d) ",i,originalRows[i]);
		    for (int j=0;j<n;j++) {
		      double value = element0[j];
		      if (j) {
			if(value>0.0)
			  printf(" +");
			else
			  printf(" ");
		      }
		      int iColumn=column0[j];
		      if ((columnFlag[iColumn]&2)==0)
			printf("%g*X%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		      else if ((columnFlag[iColumn]&(2+4))==2)
			printf("%g*B%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		      else
			printf("%g*I%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		    } 
		    if ((rowFlag[i]&1)!=0)
		      printf(" == ");
		    else
		      printf(" <= ");
		    printf("%g\n",effectiveRhs[i]);
		    printf("Comp %d (orig %d) ",k,originalRows[k]);
		    for (int j=0;j<n1;j++) {
		      double value = element1[j];
		      if (j) {
			if(value>0.0)
			  printf(" +");
			else
			  printf(" ");
		      }
		      int iColumn=column1[j];
		      if ((columnFlag[iColumn]&2)==0)
			printf("%g*X%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		      else if ((columnFlag[iColumn]&(2+4))==2)
			printf("%g*B%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		      else
			printf("%g*I%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		    } 
		    if ((rowFlag[k]&1)!=0)
		      printf(" == ");
		    else
		      printf(" <= ");
		    printf("%g\n",effectiveRhs[k]);
		  }
		  /* end debug print */
#endif
#ifdef CGL_INVESTIGATE
		  printf("**same3 %c redundant %c feasible %c tight %c,%c fixed %d,%d,%d\n",
			 same ? 'Y' : 'N',
			 redundant ? 'Y' : 'N',
			 feasible ? 'Y' : 'N',
			 tighter0 ? 'Y' : 'N',
			 tighter1 ? 'Y' : 'N',
			 fixed[0],fixed[1],fixed[2]);
#endif
		  if (!feasible) {
#ifdef CGL_INVESTIGATE
		    printf("QQ infeasible\n");
#endif
		    infeasible=true;
		  } else if (fixed[0]||fixed[1]||fixed[2]) {
#ifdef CGL_INVESTIGATE
		    printf("QQ fixed\n");
#endif
		    for (int k=0;k<3;k++) {
		      if (fixed[k]) {
			int iColumn=column0[k];
			int kColumn=originalColumns[iColumn];
#ifdef CGL_INVESTIGATE
			printf("true bounds %g %g\n",columnLower[kColumn],
			       columnUpper[kColumn]);
#endif
			double lo;
			if (fixed[k]>0) {
			  lo=1.0;
			} else {
			  lo=0.0;
			}
			if ((columnFlag[iColumn]&1)==0) {
			  columnFlag[iColumn] |= 16;
			  trueLower[iColumn] += lo;
			  newBound[iColumn] = 0.0;
			  for (CoinBigIndex jj=columnStart[iColumn];jj<columnStart[iColumn+1];
			       jj++) {
			    int iRow=row[jj];
			    effectiveRhs[iRow] -= lo*element[jj];
			  }
			} else {
			  abort();
			}
		      }
		    }
		  } else if (!same&&(tighter0||tighter1)) {
		    assert (!tighter0||!tighter1);
		    if (tighter0) {
#ifdef CGL_INVESTIGATE
		      printf("QQ discard oneT k\n");
#endif
		      duplicate_[k]=-2;
		    } else {
#ifdef CGL_INVESTIGATE
		      printf("QQ discard oneT i\n");
#endif
		      duplicate_[i]=-2;
		    }
		  } else if (redundant) {
#ifdef CGL_INVESTIGATE
		    printf("QQ discard both\n");
#endif
		    duplicate_[i]=-2;
		    duplicate_[k]=-2;
		  } else {
		    assert (same);
		    if (fabs(effectiveRhs[i]-effectiveRhs[k])<1.0e-7&&
			element0[2]==element1[2]&&
			element0[1]==element1[1]&&
			element0[0]==element1[0]) {
#ifdef CGL_INVESTIGATE
		      printf("QQ discard identical k I\n");
#endif
		      duplicate_[k]=-2;
		    } else {
#ifdef CGL_INVESTIGATE
		      printf("QQ Don't know what to do nintI==3 I\n");
#endif
		    }
		  }
		}
	      }
	    } else {
	      // one other - put last
	      double el0[3];
	      double el1[3];
	      int col[3];
	      double bound=0.0;
	      int kk=0;
	      el0[1]=0.0;
	      el1[1]=0.0;
	      for (int j=0;j<n;j++) {
		int iColumn=column0[j];
		if ((columnFlag[iColumn]&(2+4))==2) {
		  el0[kk]=element0[j];
		  col[kk++]=iColumn;
		} else {
		  el0[2]=element0[j];
		  col[2]=iColumn;
		  bound=std::min(newBound[iColumn],1.0e30);
		}
	      }
	      kk=0;
	      for (int j=0;j<n;j++) {
		int iColumn=column1[j];
		if ((columnFlag[iColumn]&(2+4))==2) {
		  el1[kk++]=element1[j];
		} else {
		  el1[2]=element1[j];
		}
	      }
	      double gap0=bound*el0[2];
	      double gap1=bound*el1[2];
	      for (kk=0;kk<4;kk++) {
		loC0[kk]=0.0;
		loC1[kk]=0.0;
		upC0[kk]=bound;
		upC1[kk]=bound;
	      }
	      // crude - should go stack based
	      double upRhs = effectiveRhs[i];
	      double loRhs = ((rowFlag[i]&1)!=0) ? effectiveRhs[i] : -1.0e30;
	      double tolerance = std::max(1.0e-5,fabs(upRhs)*1.0e-10);
	      for (int j0=0;j0<2;j0++) {
		for (int j1=0;j1<2;j1++) {
		  double value = el0[0]*j0+el0[1]*j1;
		  int put = j0+2*j1;
		  double valueLo,valueHi;
		  if (gap0>0.0) {
		    valueLo=value;
		    valueHi=value+gap0;
		  } else {
		    valueLo=value+gap0;
		    valueHi=value;
		  }
		  if (valueLo<upRhs+tolerance&&valueHi>loRhs-tolerance)
		    allowed0[put]=1;
		  else
		    allowed0[put]=0;
		  if (valueLo<loRhs-tolerance) {
		    if (gap0>0.0) {
		      loC0[put]=(loRhs-valueLo)/el0[2];
		    } else {
		      upC0[put]=bound-((valueLo-loRhs)/el0[2]);
		    }
		  }
		  if (valueHi>upRhs+tolerance) {
		    if (gap0>0.0) {
		      upC0[put]=bound-((valueHi-upRhs)/el0[2]);
		    } else {
		      loC0[put]=(upRhs-valueHi)/el0[2];
		    }
		  }
		}
	      }
	      upRhs = effectiveRhs[k];
	      loRhs = ((rowFlag[k]&1)!=0) ? effectiveRhs[k] : -1.0e30;
	      tolerance = std::max(1.0e-5,fabs(upRhs)*1.0e-10);
	      for (int j0=0;j0<2;j0++) {
		for (int j1=0;j1<2;j1++) {
		  double value = el1[0]*j0+el1[1]*j1;
		  int put = j0+2*j1;
		  double valueLo,valueHi;
		  if (gap1>0.0) {
		    valueLo=value;
		    valueHi=value+gap1;
		  } else {
		    valueLo=value+gap1;
		    valueHi=value;
		  }
		  if (valueLo<upRhs+tolerance&&valueHi>loRhs-tolerance)
		    allowed1[put]=1;
		  else
		    allowed1[put]=0;
		  if (valueLo<loRhs-tolerance) {
		    if (gap1>0.0) {
		      loC1[put]=(loRhs-valueLo)/el1[2];
		    } else {
		      upC1[put]=bound-((valueLo-loRhs)/el1[2]);
		    }
		  }
		  if (valueHi>upRhs+tolerance) {
		    if (gap1>0.0) {
		      upC1[put]=bound-((valueHi-upRhs)/el1[2]);
		    } else {
		      loC1[put]=(upRhs-valueHi)/el1[2];
		    }
		  }
		}
	      }
	      /* interesting cases are when -
		 one forces fixing (but probably found in probing)
		 and of two is zero
		 and of two forces fixing
		 two same - this is probably only one
	      */
	      int intersect[4];
	      bool same=true;
	      bool tighter0=true;
	      bool tighter1=true;
	      bool feasible=false;
	      bool redundant=true;
	      int count=nInt;
	      for (int j=0;j<(1<<count);j++) {
		intersect[j]=allowed0[j]&allowed1[j];
		if (intersect[j]<allowed0[j])
		  tighter0=false;
		if (intersect[j]<allowed1[j])
		  tighter1=false;
		if (allowed0[j]!=allowed1[j])
		  same=false;
		if (!intersect[j])
		  redundant=false;
		if (intersect[j])
		  feasible=true;
	      }
	      int fixed[2]={0,0};
	      if (feasible) {
		for (int jj=0;jj<count;jj++) {
		  int multiplier=1<<jj;
		  int increment=1<<(count-1-jj);
		  bool zeroOk=false;
		  bool oneOk=false;
		  for (int j=0;j<(1<<(count-1));j++) {
		    if (intersect[0*multiplier+increment*j])
		      zeroOk=true;
		    if (intersect[1*multiplier+increment*j])
		      oneOk=true;
		  }
		  if (!zeroOk) {
		    fixed[jj]=1;
		    assert (oneOk);
		  } else if (!oneOk) {
		    fixed[jj]=-1;
		  }
		}
	      }
	      double newLo=bound;
	      double newUp=0.0;
	      if (nInt==1) {
		for (int jj=0;jj<2;jj++) {
#ifdef CGL_INVESTIGATE
		  printf("int at %d -> lo0 %g lo1 %g up0 %g up1 %g",
			 jj,loC0[jj],loC1[jj],upC0[jj],upC1[jj]);
#endif
		  if (intersect[jj]) {
#ifdef CGL_INVESTIGATE
		    printf("\n");
#endif
		    newLo=std::min(newLo,std::max(loC0[jj],loC1[jj]));
		    newUp=std::max(newUp,std::min(upC0[jj],upC1[jj]));
		  } else {
#ifdef CGL_INVESTIGATE
		    printf(" INF\n");
#endif
		    loC0[jj]=0.0;
		    loC1[jj]=0.0;
		    upC0[jj]=bound;
		    upC1[jj]=bound;
		  }
		}
	      } else {
		for (int jj=0;jj<2;jj++) { 
		  for (int jj1=0;jj1<2;jj1++) {
		    int k=jj+2*jj1;
#ifdef CGL_INVESTIGATE
		    printf("first int at %d, second at %d -> lo0 %g lo1 %g up0 %g up1 %g",
			   jj,jj1,loC0[k],loC1[k],upC0[k],upC1[k]);
#endif
		    if (intersect[k]) {
#ifdef CGL_INVESTIGATE
		      printf("\n");
#endif
		      newLo=std::min(newLo,std::max(loC0[k],loC1[k]));
		      newUp=std::max(newUp,std::min(upC0[k],upC1[k]));
		    } else {
#ifdef CGL_INVESTIGATE
		      printf(" INF\n");
#endif
		      loC0[k]=0.0;
		      loC1[k]=0.0;
		      upC0[k]=bound;
		      upC1[k]=bound;
		    }
		  }
		}
	      }
	      if (newLo>0.0||newUp<bound) {
#ifdef CGL_INVESTIGATE
		printf("Can tighten bounds to %g,%g\n",newLo,newUp);
#endif
		int iColumn=col[2];
		int kColumn=originalColumns[iColumn];
#ifdef CGL_INVESTIGATE
		printf("true bounds %g %g\n",columnLower[kColumn],
		       columnUpper[kColumn]);
#endif
		if ((columnFlag[iColumn]&1)==0) {
		  columnFlag[iColumn] |= 16;
		  trueLower[iColumn] += newLo;
		  newBound[iColumn] = newUp-newLo;
		  for (CoinBigIndex jj=columnStart[iColumn];jj<columnStart[iColumn+1];
		       jj++) {
		    int iRow=row[jj];
		    effectiveRhs[iRow] -= newLo*element[jj];
		  }
		} else {
		  abort();
		}
	      }
	      for (int jj=0;jj<(1<<nInt);jj++) {
		loC0[jj]=std::max(loC0[jj],newLo);
		loC1[jj]=std::max(loC1[jj],newLo);
		upC0[jj]=std::min(upC0[jj],newUp);
		upC1[jj]=std::min(upC1[jj],newUp);
	      }
	      for (int jj=0;jj<(1<<nInt);jj++) {
		if (fabs(loC0[jj]-loC1[jj])>1.0e-8)
		  same=false;
		if (fabs(upC0[jj]-upC1[jj])>1.0e-8)
		  same=false;
		if (loC0[jj]<loC1[jj]-1.0e-12||
		    upC0[jj]>upC1[jj]+1.0e-12) {
		  tighter0=false;
		  redundant=false;
		}
		if (loC1[jj]<loC0[jj]-1.0e-12||
		    upC1[jj]>upC0[jj]+1.0e-12) {
		  tighter1=false;
		  redundant=false;
		}
		if (fabs(newLo-loC0[jj])>1.0e-12)
		  redundant=false;
		if (fabs(newLo-loC1[jj])>1.0e-12)
		  redundant=false;
		if (fabs(newUp-upC0[jj])>1.0e-12)
		  redundant=false;
		if (fabs(newUp-upC1[jj])>1.0e-12)
		  redundant=false;
	      }
	      if (same||redundant||!feasible||fixed[0]||fixed[1]||
		  tighter0||tighter1) {
#ifdef CGL_INVESTIGATE
		/* start debug print */
		{
		  printf("Base %d (orig %d) ",i,originalRows[i]);
		  for (int j=0;j<n;j++) {
		    double value = element0[j];
		    if (j) {
		      if(value>0.0)
			printf(" +");
		      else
			printf(" ");
		    }
		    int iColumn=column0[j];
		    if ((columnFlag[iColumn]&2)==0)
		      printf("%g*X%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		    else if ((columnFlag[iColumn]&(2+4))==2)
		      printf("%g*B%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		    else
		      printf("%g*I%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		  } 
		  if ((rowFlag[i]&1)!=0)
		    printf(" == ");
		  else
		    printf(" <= ");
		  printf("%g\n",effectiveRhs[i]);
		  printf("Comp %d (orig %d) ",k,originalRows[k]);
		  for (int j=0;j<n1;j++) {
		    double value = element1[j];
		    if (j) {
		      if(value>0.0)
			printf(" +");
		      else
			printf(" ");
		    }
		    int iColumn=column1[j];
		    if ((columnFlag[iColumn]&2)==0)
		      printf("%g*X%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		    else if ((columnFlag[iColumn]&(2+4))==2)
		      printf("%g*B%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		    else
		      printf("%g*I%d(%d) (<=%g)",value,iColumn,originalColumns[iColumn],newBound[iColumn]);
		  } 
		  if ((rowFlag[k]&1)!=0)
		    printf(" == ");
		  else
		    printf(" <= ");
		  printf("%g\n",effectiveRhs[k]);
		}
		/* end debug print */
#endif
#ifdef CGL_INVESTIGATE
		printf("**same %c redundant %c feasible %c tight %c,%c fixed %d,%d\n",
		       same ? 'Y' : 'N',
		       redundant ? 'Y' : 'N',
		       feasible ? 'Y' : 'N',
		       tighter0 ? 'Y' : 'N',
		       tighter1 ? 'Y' : 'N',
		       fixed[0],fixed[1]);
#endif
		//if (nInt>1)
		//continue;
		if (!feasible) {
#ifdef CGL_INVESTIGATE
		  printf("QQ infeasible\n");
#endif
		  infeasible=true;
		} else if (fixed[0]||fixed[1]) {
#ifdef CGL_INVESTIGATE
		  printf("QQ fixed\n");
#endif
		  for (int k=0;k<2;k++) {
		    if (fixed[k]) {
		      int iColumn=col[k];
		      int kColumn=originalColumns[iColumn];
#ifdef CGL_INVESTIGATE
		      printf("true bounds %g %g\n",columnLower[kColumn],
			     columnUpper[kColumn]);
#endif
		      double lo;
		      if (fixed[k]>0) {
			lo=1.0;
		      } else {
			lo=0.0;
		      }
		      if ((columnFlag[iColumn]&1)==0) {
			columnFlag[iColumn] |= 16;
			trueLower[iColumn] += lo;
			newBound[iColumn] = 0.0;
			for (CoinBigIndex jj=columnStart[iColumn];jj<columnStart[iColumn+1];
			     jj++) {
			  int iRow=row[jj];
			  effectiveRhs[iRow] -= lo*element[jj];
			}
		      } else {
			abort();
		      }
		    }
		  }
		} else if (!same&&(tighter0||tighter1)) {
		  assert (!tighter0||!tighter1);
		  if (tighter0) {
#ifdef CGL_INVESTIGATE
		    printf("QQ discard oneT k\n");
#endif
		    duplicate_[k]=-2;
		  } else {
#ifdef CGL_INVESTIGATE
		    printf("QQ discard oneT i\n");
#endif
		    duplicate_[i]=-2;
		  }
		} else if (redundant) {
#ifdef CGL_INVESTIGATE
		  printf("QQ discard both\n");
#endif
		  duplicate_[i]=-2;
		  duplicate_[k]=-2;
		} else {
		  assert (same);
		  if (fabs(effectiveRhs[i]-effectiveRhs[k])<1.0e-7&&
		      el0[2]==el1[2]&&
		      el0[0]==el1[0]) {
		    if (nInt==1||el0[1]==el1[1]) {
#ifdef CGL_INVESTIGATE
		      printf("QQ discard identical k\n");
#endif
		      duplicate_[k]=-2;
		    }
		  }
		  if (duplicate_[k]!=-2) {
		    if (nInt==1) {
		      // one may be stronger
		      if (el0[2]>0.0&&el1[2]>0.0&&
			  el0[0]>0.0&&el1[0]>0.0&&
			  (rowFlag[k]&1)==0&&
			  (rowFlag[i]&1)==0) {
			// bounds same at 0 and 1
			double up0 = (effectiveRhs[i]-el0[0])/el0[2];
			double up1 = (effectiveRhs[k]-el1[0])/el1[2];
			if (up0<up1) {
#ifdef CGL_INVESTIGATE
			  printf("QQ discard oneS k\n");
#endif
			  duplicate_[k]=-2;
			} else {
#ifdef CGL_INVESTIGATE
			  printf("QQ discard oneS i\n");
#endif
			  duplicate_[i]=-2;
			}
		      } else {
			if ((rowFlag[k]&1)==0&&(rowFlag[i]&1)==0) {
			  int which=0;
			  for (int i0=0;i<=10;i++) {
			    double value0=0.05*i0;
			    double rhs0 =effectiveRhs[i]-el0[0]*value0;
			    double lo0=0.0;
			    double up0=1.0e30;
			    double bound0=rhs0/el0[2];
			    if (el0[2]>0.0) 
			      up0=bound0;
			    else
			      lo0=std::max(0.0,bound0);
			    double rhs1 =effectiveRhs[k]-el1[0]*value0;
			    double lo1=0.0;
			    double up1=1.0e30;
			    double bound1=rhs1/el1[2];
			    if (el1[2]>0.0) 
			      up1=bound1;
			    else
			      lo1=std::max(0.0,bound1);
			    if (fabs(lo0-lo1)>1.0e-8||
				fabs(up0-up1)>1.0e-8*(1.0+fabs(up1))) {
			      if (lo0>lo1+1.0e-8) {
				if (up0<up1-1.0e-8) {
				  // 0 tighter
				  if (which==1) {
				    which=-2;
				    break;
				  } else {
				    which=-1;
				  }
				} else if (up0>up1+1.0e-8) {
				  which=-2;
				  break;
				}
			      } else if (lo0<lo1-1.0e-8) {
				if (up1<up0-1.0e-8) {
				  // 1 tighter
				  if (which==-1) {
				    which=-2;
				    break;
				  } else {
				    which=1;
				  }
				} else if (up1>up0+1.0e-8) {
				  which=-2;
				  break;
				}
			      } else {
				if (up1<up0-1.0e-8) {
				  // 1 tighter
				  if (which==-1) {
				    which=-2;
				    break;
				  } else {
				    which=1;
				  }
				} else if (up1>up0+1.0e-8) {
				  // 0 tighter
				  if (which==1) {
				    which=-2;
				    break;
				  } else {
				    which=-1;
				  }
				}
			      }
			    }
			  }
			  if (which==0) {
			    duplicate_[k]=-2;
#ifdef CGL_INVESTIGATE
			  printf("QQ discard one same same k\n");
#endif
			  } else if (which==-1) {
			    duplicate_[k]=-2;
#ifdef CGL_INVESTIGATE
			  printf("QQ discard one (i tighter) k\n");
#endif
			  } else if (which==1) {
			    duplicate_[i]=-2;
#ifdef CGL_INVESTIGATE
			  printf("QQ discard one (k tighter) i\n");
#endif
			  } else {
			    printf("QQ Can't decide what to do nint==1\n");
			  }
			} else {
			  printf("QQ Don't know what to do nint==1\n");
			}
		      }
		    } else {
		      if ((rowFlag[k]&1)==0&&(rowFlag[i]&1)==0) {
			int which=0;
			for (int i0=0;i0<=10;i0++) {
			  double value0=0.05*i0;
			  for (int i1=0;i1<=10;i1++) {
			    double value1=0.05*i1;
			    double rhs0 =effectiveRhs[i]-el0[0]*value0
			      -el0[1]*value1;
			    double lo0=0.0;
			    double up0=1.0e30;
			    double bound0=rhs0/el0[2];
			    if (el0[2]>0.0) 
			      up0=bound0;
			    else
			      lo0=std::max(0.0,bound0);
			    double rhs1 =effectiveRhs[k]-el1[0]*value0
			      -el1[1]*value1;
			    double lo1=0.0;
			    double up1=1.0e30;
			    double bound1=rhs1/el1[2];
			    if (el1[2]>0.0) 
			      up1=bound1;
			    else
			      lo1=std::max(0.0,bound1);
			    if (fabs(lo0-lo1)>1.0e-8||
				fabs(up0-up1)>1.0e-8*(1.0+fabs(up1))) {
			      if (lo0>lo1+1.0e-8) {
				if (up0<up1-1.0e-8) {
				  // 0 tighter
				  if (which==1) {
				    which=-2;
				    break;
				  } else {
				    which=-1;
				  }
				} else if (up0>up1+1.0e-8) {
				  which=-2;
				  break;
				}
			      } else if (lo0<lo1-1.0e-8) {
				if (up1<up0-1.0e-8) {
				  // 1 tighter
				  if (which==-1) {
				    which=-2;
				    break;
				  } else {
				    which=1;
				  }
				} else if (up1>up0+1.0e-8) {
				  which=-2;
				  break;
				}
			      } else {
				if (up1<up0-1.0e-8) {
				  // 1 tighter
				  if (which==-1) {
				    which=-2;
				    break;
				  } else {
				    which=1;
				  }
				} else if (up1>up0+1.0e-8) {
				  // 0 tighter
				  if (which==1) {
				    which=-2;
				    break;
				  } else {
				    which=-1;
				  }
				}
			      }
			    }
			  }
			  if (which==-2)
			    break;
			}
			if (which==0) {
			  duplicate_[k]=-2;
#ifdef CGL_INVESTIGATE
			  printf("QQ discard one same same k\n");
#endif
			} else if (which==-1) {
			  duplicate_[k]=-2;
#ifdef CGL_INVESTIGATE
			  printf("QQ discard one (i tighter) k\n");
#endif
			} else if (which==1) {
			  duplicate_[i]=-2;
#ifdef CGL_INVESTIGATE
			  printf("QQ discard one (k tighter) i\n");
#endif
			} else {
			  printf("QQ Can't decide what to do nint==2\n");
			}
		      } else {
#ifdef CGL_INVESTIGATE
		      printf("QQ Don't know what to do nint==2\n");
#endif
		      }
		    }
		  }
		}
	      }
	    }	      
	  }
	}
      }
      for (int j=0;j<n;j++) 
	mark[column0[j]] =0;
    }
  }
  if (0) {
    // Column copy
    const double * element = matrix_.getElements();
    const int * row = matrix_.getIndices();
    const CoinBigIndex * columnStart = matrix_.getVectorStarts();
    //const int * columnLength = matrix_.getVectorLengths();
    // Row copy
    const double * elementByRow = matrixByRow_.getElements();
    const int * column = matrixByRow_.getIndices();
    const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
    //const int * rowLength = matrixByRow_.getVectorLengths();
    for (int i=1231;i<1235;i++) {
      for (CoinBigIndex jj=columnStart[i];jj<columnStart[i+1];
	   jj++) {
	int iRow=row[jj];
	printf("row %d el %g\n",iRow,element[jj]);
      }
    }
    for (int i=283;i<284;i++) {
      for (CoinBigIndex jj=rowStart[i];jj<rowStart[i+1];
	   jj++) {
	printf("col %d el %g\n",column[jj],elementByRow[jj]);
      }
    }
    printf("OK?\n");
  }
  // See if any fixed
  CoinPackedVector lbs;
  CoinPackedVector ubs;
  //OsiSolverInterface * xx = si.clone();
  for (int i=0;i<nCol2;i++) {
    if ((columnFlag[i]&16)!=0) {
      assert ((columnFlag[i]&1)==0);
      int kColumn=originalColumns[i];
      double lower=trueLower[i];
      double upper=lower+newBound[i];
      if (fabs(lower-floor(lower+0.5))<1.0e-8) 
	lower=floor(lower+0.5);
      if (fabs(upper-floor(upper+0.5))<1.0e-8) 
	upper=floor(upper+0.5);
      if ((columnFlag[i]&2)==0) {
	// continuous
	if (lower>columnLower[kColumn]+1.0e-7) {
	  nFixed++;
	  lbs.insert(kColumn,lower);
	  //xx->setColLower(kColumn,lower);
	}
	if (upper<columnUpper[kColumn]-1.0e-7) {
	  nFixed++;
	  ubs.insert(kColumn,upper);
	  //xx->setColUpper(kColumn,upper);
	}
      } else {
	// integer
	if (lower>columnLower[kColumn]+1.0e-7) {
	  nFixed++;
	  lower = ceil(lower);
	  //xx->setColLower(kColumn,lower);
	  lbs.insert(kColumn,lower);
	}
	if (upper<columnUpper[kColumn]-1.0e-7) {
	  upper=floor(upper);
	  nFixed++;
	  //xx->setColUpper(kColumn,upper);
	  ubs.insert(kColumn,upper);
	}
      }
      if (lower>upper+1.0e-7)
	infeasible=true;
      //printf("Bounds for %d are %g and %g, were %g %g\n",
      //     kColumn,
      //     xx->getColLower()[kColumn],
      //     xx->getColUpper()[kColumn],
      //     si.getColLower()[kColumn],
      //     si.getColUpper()[kColumn]);
    }
  }
  //if (!infeasible) {
  //si.writeMps("si");
  //xx->writeMps("xx");
  //xx->resolve();
  //assert (xx->isProvenOptimal());
  //printf("obj value %g %g\n",si.getObjValue(),xx->getObjValue());
  //}
  //delete xx;
  // Move duplicate flags
  int * temp = reinterpret_cast<int *>(effectiveRhs);
  for (int i=0;i<numberRows;i++)
    temp[i]=-1;
  for (int i=0;i<nRow2;i++)
    temp[originalRows[i]]=duplicate_[i];
  memcpy(duplicate_,temp,numberRows*sizeof(int));
  delete [] effectiveRhs;
  delete [] newBound;
  delete [] trueLower;
  delete [] rowFlag;
  delete [] columnFlag;
  delete [] originalColumns;
  delete [] originalRows;
  if (nFixed) {
#ifdef CGL_INVESTIGATE
    printf("QQ - %d bounds changed\n",nFixed);
#endif
    OsiColCut cc;
    cc.setLbs(lbs);
    cc.setUbs(ubs);
    cc.setEffectiveness(100.0);
    cs.insert(cc);
  }
  if (infeasible) {
    // generate infeasible cut and return
    printf("QQ**** infeasible cut\n");
    OsiRowCut rc;
    rc.setLb(COIN_DBL_MAX);
    rc.setUb(0.0);   
    cs.insertIfNotDuplicate(rc);
  }
}
#if 0
class CglOneRow {

public:
  /// data
  /// Row number
  int row_;
  /// Start
  const int * start_;
  /// End
  const int * end_;

public:

    // Default Constructor
  inline CglOneRow () : row_(-1),start_(0), end_(0) {} 

  // Useful constructor
  inline CglOneRow (int iRow,const int * start, const int * end) : row_(iRow),start_(start), end_(end) {}
  // Destructor
  inline ~CglOneRow () {}

};
class CglCompare {
public:
  /// Compare function
  inline bool operator()(const CglOneRow & row1,
			 const CglOneRow & row2) const
  { 
    const int * where1 = row1.start_;
    const int * where2 = row2.start_;
    while (where1 != row1.end_ && where2 != row2.end_) {
      int iColumn1 = *where1;
      int iColumn2 = *where2;
      if (iColumn1<iColumn2) {
	return true;
      } else if (iColumn1>iColumn2) {
	return false;
      } else {
	where1++;
	where2++;
      }
    }
    if (where1==row1.end_)
      return false;
    else
      return true;
  }
};
static int * lexSort(int numberCliques,
		    int * cliqueStart, int * entry)
{
  CglOneRow * rows = new CglOneRow [numberCliques];
  for (int i=0;i<numberCliques;i++) {
    rows[i]=CglOneRow(i,entry+cliqueStart[i],entry+cliqueStart[i+1]);
  }
  std::sort(rows,rows+numberCliques,CglCompare());
  int * sorted = new int [numberCliques];
  for (int i=0;i<numberCliques;i++) {
    sorted[i]=rows[i].row_;
  }
  delete [] rows;
  return sorted;
}
static int outDupsEtc2(int numberIntegers, int numberCliques, int * statusClique,
		      int * cliqueStart, char * cliqueType, int * entry, 
		      int printit)
{
  int * sorted = lexSort(numberCliques,cliqueStart,entry);
  delete [] sorted;
  return 0;
}
#endif
static int outDupsEtc(int numberIntegers, int numberCliques, int * statusClique,
		      int * cliqueStart, char * cliqueType, int * entry,
		      int * fixed,
		      int printit)
#if 0
{
  //outDupsEtc2(numberIntegers,numberCliques,statusClique,
  //      cliqueStart,cliqueType,entry,printit);
  int * whichP = new int [numberIntegers];
  int iClique;
  assert (sizeof(int)==4);
  assert (sizeof(int)==4);
  // sort
  for (iClique=0;iClique<numberCliques;iClique++) {
    int j = cliqueStart[iClique];
    int n = cliqueStart[iClique+1]-j;
    for (int i=0;i<n;i++) 
      whichP[i]=entry[i+j];
    CoinSort_2(whichP,whichP+n,entry+j);
  }
  // lexicographic sort
  int * which = new int [numberCliques];
  int * position = new int [numberCliques];
  int * sort = new int [numberCliques];
  for (iClique=0;iClique<numberCliques;iClique++) {
    which[iClique]=iClique;
    sort[iClique]=entry[cliqueStart[iClique]];
    statusClique[iClique]=sort[iClique];
    position[iClique]=0;
  }
  CoinSort_2(sort,sort+numberCliques,which);
  int lastDone=-1;
  int nDup=0;
  int nSave=0;
  while (lastDone<numberCliques-1) {
    int jClique=lastDone+1;
    int jFirst = jClique;
    int iFirst = which[jFirst];
    int iValue = statusClique[iFirst];
    int iPos = position[iFirst];
    jClique++;
    for (;jClique<numberCliques;jClique++) {
      int kClique = which[jClique];
      int jValue = statusClique[kClique];
      if (jValue>iValue||position[kClique]<iPos)
	break;
    }
    if (jClique==jFirst+1) {
      // done that bit
      lastDone++;
    } else {
      // use next bit to sort and then repeat
      int jLast=jClique;
      for (jClique=jFirst;jClique<jLast;jClique++) {
	int kClique = which[jClique];
	int iValue = statusClique[kClique];
	// put at end if finished
	if (iValue<numberIntegers) {
	  int kPos=position[kClique]+1;
	  position[kClique]=kPos;
	  kPos += cliqueStart[kClique];
	  if (kPos==cliqueStart[kClique+1]) {
	    iValue = numberIntegers;
	  } else {
	    iValue = entry[kPos];
	  }
	  statusClique[kClique]=iValue;
	}
	sort[jClique]=iValue;
      }
      CoinSort_2(sort+jFirst,sort+jLast,which+jFirst);
      // if duplicate mark and move on
      int iLowest=numberCliques;
      char type='S';
      for (jClique=jFirst;jClique<jLast;jClique++) {
	int kClique = which [jClique];
	int iValue = statusClique[kClique];
	if (iValue<numberIntegers) 
	  break;
	if (cliqueType[kClique]=='E') {
	  iLowest = std::min(iLowest,kClique);
	  type='E';
	} else if (type=='S') {
	  iLowest = std::min(iLowest,kClique);
	}
      }
      if (jClique>jFirst) {
	// mark all apart from lowest number as duplicate and move on
	// use cliqueType
	lastDone =jClique-1;
	for (jClique=jFirst;jClique<=lastDone;jClique++) {
	  int kClique = which [jClique];
	  if (kClique!=iLowest) {
	    statusClique[kClique]=-2;
	    nDup++;
	    nSave += cliqueStart[kClique+1]-cliqueStart[kClique];
	  }
	}
      }
    }
  }
#if 1
  for (int jClique=0;jClique<numberCliques;jClique++) {
    int iClique=which[jClique];
    printf("clique %d %d ",jClique,iClique);
    for (int j=cliqueStart[iClique];j<cliqueStart[iClique+1];j++) {
      int iColumn = entry[j];
      printf("%d ",iColumn);
    }
    printf("\n");
  }
#endif
  if (printit)
    printf("%d duplicates\n",nDup); 
  // For column version
  int numberElements=cliqueStart[numberCliques];
  int * start = new int [numberIntegers];
  int * end = new int [numberIntegers];
  int * clique = new int [numberElements];
  int * marked = new int [numberCliques];
  int * count = new int [numberCliques];
  int * fixed = new int [numberIntegers];
  memset(count,0,numberCliques*sizeof(int));
  memset(end,0,numberIntegers*sizeof(int));
  memset(fixed,0,numberIntegers*sizeof(int));
  nSave=0;
  for (int jClique=0;jClique<numberCliques;jClique++) {
    if (statusClique[jClique]!=-2) {
      for (int j=cliqueStart[jClique];j<cliqueStart[jClique+1];j++) {
	int iColumn = entry[j];
	end[iColumn]++;
      }
    } else {
      nSave += cliqueStart[jClique+1]-cliqueStart[jClique];
    }
  }
  numberElements=0;
  for (int i=0;i<numberIntegers;i++) {
    start[i]=numberElements;
    int n=end[i];
    end[i]=numberElements;
    numberElements += n;
  }
  for (int jClique=0;jClique<numberCliques;jClique++) {
    int iClique=which[jClique];
    if (statusClique[iClique]!=-2) {
      for (int j=cliqueStart[iClique];j<cliqueStart[iClique+1];j++) {
	int iColumn = entry[j];
	int put=end[iColumn];
	end[iColumn]=put+1;
	clique[put]=iClique;
      }
    }
  }

  // Now see if any subset
  int nOut=0;
  int nFixed=0;
  for (int jClique=0;jClique<numberCliques;jClique++) {
    int kClique = which[jClique];
    if (statusClique[kClique]==-2) 
      continue;
    // do first
    int nMarked=0;
    int iEnd=cliqueStart[kClique+1];
    int iEl = cliqueStart[kClique];
    int iColumn=entry[iEl++];
    while (fixed[iColumn]) {
      iColumn=-1;
      if(iEl<iEnd) {
	iColumn=entry[iEl++];
      } else {
	break;
      }
    }
    if (iColumn<0) {
      // now empty
      printf("now empty clique %d !\n",kClique);
      statusClique[kClique]=-2;
      continue;
    }
    int i;
    for ( i=start[iColumn];i<end[iColumn];i++) {
      int iClique=clique[i];
      // faster to shuffle up -2's?
      if (statusClique[iClique]!=-2) {
	if (iClique!=kClique) {
	  count[iClique]=1;
	  marked[nMarked++]=iClique;
	} else {
	  break;
	}
      }
    }
    assert (i<end[iColumn]);
    if (nMarked) {
      int n=nMarked;
      int sizeClique=1;
      for (;iEl<iEnd;iEl++) {
	int iColumn=entry[iEl];
	for ( i=start[iColumn];i<end[iColumn];i++) {
	  int iClique=clique[i];
	  // faster to shuffle up -2's?
	  if (statusClique[iClique]!=-2) {
	    if (iClique!=kClique) {
	      if (count[iClique]) {
		if (count[iClique]==sizeClique) {
		  count[iClique]++;
		} else {
		  count[iClique]=0;
		  n--;
		}
	      }
	    } else {
	      break;
	    }
	  }
	}
	sizeClique++;
	assert (i<end[iColumn]);
	if (!n)
	  break;
      }
      if (n) {
	// still some left
	// But need to look at type
	// when might be able to fix variables
	bool subset=false;
	if (cliqueType[kClique]=='E') {
	  for (i=0;i<nMarked;i++) {
	    int iClique=marked[i];
	    if (count[iClique]==sizeClique) {
	      subset=true;
	      // can fix all in iClique not in kClique
	      int iEl=cliqueStart[kClique];
	      int jEl=cliqueStart[iClique];
	      int jEnd=cliqueStart[iClique+1];
	      int iColumn=entry[iEl];
	      for (int j=jEl;j<jEnd;j++) {
		int jColumn=entry[j];
		if (jColumn==iColumn) {
		  iEl++;
		  if (iEl<iEnd)
		    iColumn=entry[iEl];
		  else
		    iColumn=numberIntegers;
		} else if (!fixed[jColumn]) {
		  // fix
		  nFixed++;
		  fixed[jColumn]=1;
		}
	      }
	    }
	  }
	} else {
	  // print first for now
	  for (i=0;i<nMarked;i++) {
	    int iClique=marked[i];
	    if (count[iClique]==sizeClique) {
	      subset=true;
#if 0
	      if (printit>10
		  printf("clique %d is subset of %d\n",kClique,iClique);
	      printf("Kclique %d ",kClique);
	      for (int j=cliqueStart[kClique];j<cliqueStart[kClique+1];j++) {
		int kColumn = entry[j];
		printf("%d ",kColumn);
	      }
	      printf("\n");
	      printf("Iclique %d ",iClique);
	      for (int j=cliqueStart[iClique];j<cliqueStart[iClique+1];j++) {
		int iColumn = entry[j];
		printf("%d ",iColumn);
	      }
	      printf("\n");
#endif
	      break;
	    }
	  }
	}
	if (subset) {
	  nOut++;
	  statusClique[kClique]=-2;
	}
      }
      for (i=0;i<nMarked;i++)
	count[marked[i]]=0;
    }
  }
  if (nOut) {
    if(printit) 
      printf("Can get rid of %d cliques\n",nOut);
  }
  if (nFixed) {
    printf("Can fix to zero ");
    // numbers are subset of column variables
    for (int i=0;i<numberIntegers;i++) {
      if (fixed[i])
	printf("%d ",i);
    }
    printf("\n");
    abort();
  }
  delete [] sort;
  delete [] which;
  delete [] position;
  delete [] whichP;
  delete [] start;
  delete [] end;
  delete [] clique;
  delete [] marked;
  delete [] count;
  delete [] fixed;
  return nOut;
}
#else
{
  //outDupsEtc2(numberIntegers,numberCliques,statusClique,
  //      cliqueStart,cliqueType,entry,printit);
  int * whichP = new int [numberIntegers];
  int iClique;
  assert (sizeof(int)==4);
  assert (sizeof(int)==4);
  // sort
  for (iClique=0;iClique<numberCliques;iClique++) {
    int j = cliqueStart[iClique];
    int n = cliqueStart[iClique+1]-j;
    for (int i=0;i<n;i++) 
      whichP[i]=entry[i+j];
    CoinSort_2(whichP,whichP+n,entry+j);
  }
  // lexicographic sort
  int * which = new int [numberCliques];
  int * position = new int [numberCliques];
  int * sort = new int [numberCliques];
  for (iClique=0;iClique<numberCliques;iClique++) {
    which[iClique]=iClique;
    sort[iClique]=entry[cliqueStart[iClique]];
    statusClique[iClique]=sort[iClique];
    position[iClique]=0;
  }
  CoinSort_2(sort,sort+numberCliques,which);
  int lastDone=-1;
  int nDup=0;
  //int nSave=0;
  while (lastDone<numberCliques-1) {
    int jClique=lastDone+1;
    int jFirst = jClique;
    int iFirst = which[jFirst];
    int iValue = statusClique[iFirst];
    int iPos = position[iFirst];
    jClique++;
    for (;jClique<numberCliques;jClique++) {
      int kClique = which[jClique];
      int jValue = statusClique[kClique];
      if (jValue>iValue||position[kClique]<iPos)
	break;
    }
    if (jClique==jFirst+1) {
      // done that bit
      lastDone++;
    } else {
      // use next bit to sort and then repeat
      int jLast=jClique;
      for (jClique=jFirst;jClique<jLast;jClique++) {
	int kClique = which[jClique];
	int iValue = statusClique[kClique];
	// put at end if finished
	if (iValue<numberIntegers) {
	  int kPos=position[kClique]+1;
	  position[kClique]=kPos;
	  kPos += cliqueStart[kClique];
	  if (kPos==cliqueStart[kClique+1]) {
	    iValue = numberIntegers;
	  } else {
	    iValue = entry[kPos];
	  }
	  statusClique[kClique]=iValue;
	}
	sort[jClique]=iValue;
      }
      CoinSort_2(sort+jFirst,sort+jLast,which+jFirst);
      // if duplicate mark and move on
      int iLowest=numberCliques;
      char type='S';
      for (jClique=jFirst;jClique<jLast;jClique++) {
	int kClique = which [jClique];
	int iValue = statusClique[kClique];
	if (iValue<numberIntegers) 
	  break;
	if (cliqueType[kClique]=='E') {
	  iLowest = std::min(iLowest,kClique);
	  type='E';
	} else if (type=='S') {
	  iLowest = std::min(iLowest,kClique);
	}
      }
      if (jClique>jFirst) {
	// mark all apart from lowest number as duplicate and move on
	// use cliqueType
	lastDone =jClique-1;
	for (jClique=jFirst;jClique<=lastDone;jClique++) {
	  int kClique = which [jClique];
	  if (kClique!=iLowest) {
	    statusClique[kClique]=-2;
	    nDup++;
	    //nSave += cliqueStart[kClique+1]-cliqueStart[kClique];
	  }
	}
      }
    }
  }
#if 0
  for (int jClique=0;jClique<numberCliques;jClique++) {
    int iClique=which[jClique];
    printf("clique %d %d ",jClique,iClique);
    for (int j=cliqueStart[iClique];j<cliqueStart[iClique+1];j++) {
      int iColumn = entry[j];
      printf("%d ",iColumn);
    }
    printf("\n");
  }
#endif
  if (printit)
    printf("%d duplicates\n",nDup);
  int nOutMax=2000000;
  // mark cliques used to remove other cliques
  int * used = statusClique+numberCliques;
  // Now see if any subset
  int nOut=0;
  for (int jClique=0;jClique<numberCliques;jClique++) {
    used[jClique]=numberCliques;
    if (statusClique[jClique]!=-2) {
      position[jClique]=cliqueStart[jClique];
      statusClique[jClique]=entry[cliqueStart[jClique]];
    }
  }
  //nSave=0;
  int startLooking=0;
  for (int jClique=0;jClique<numberCliques;jClique++) {
    int kClique = which[jClique];
    if (statusClique[kClique]==-2) {
      nOut++;
      //nSave += cliqueStart[kClique+1]-cliqueStart[kClique];
      if (jClique==startLooking)
	startLooking++;
      continue;
    }
    int kValue =statusClique[kClique];
    bool ppp=false;
    for (int iiClique=startLooking;iiClique<jClique;iiClique++) {
      int iClique = which[iiClique];
      int iValue = statusClique[iClique];
      if (iValue==-2||iValue==numberIntegers) {
	if (iiClique==startLooking)
	  startLooking++;
	continue;
      } else {
	if (kValue>entry[cliqueStart[iClique+1]-1]) {
	  statusClique[iClique]=numberIntegers;
	  continue;
	}
      }
      if (iValue<kValue) {
	while (iValue<kValue) {
	  int iPos=position[iClique]+1;
	  position[iClique]=iPos;
	  if (iPos==cliqueStart[iClique+1]) {
	    iValue = numberIntegers;
	  } else {
	    iValue = entry[iPos];
	  }
	  statusClique[iClique]=iValue;
	}
      } 
      if (iValue>kValue) 
	continue; // not a candidate
      // See if subset (remember duplicates have gone)
      if (cliqueStart[iClique+1]-position[iClique]>
	  cliqueStart[kClique+1]-cliqueStart[kClique]) {
	// could be subset ?
	int offset = cliqueStart[iClique]-position[kClique];
	int j;
	bool subset=true;
	// what about different fixes bool odd=false;
	for (j=cliqueStart[kClique];j<cliqueStart[kClique+1];j++) {
	  int kColumn = entry[j];
	  int iColumn = entry[j+offset];
	  while (iColumn<kColumn) {
	    offset++;
	    if (j+offset<cliqueStart[iClique+1]) {
	      iColumn = entry[j+offset];
	    } else {
	      iColumn=numberIntegers;
	    }
	  }
	  if (iColumn!=kColumn) {
	    subset=false;
	    break;
	  }
	}
	if (subset&&nOut<=nOutMax) {
	  int kSave=statusClique[kClique];
	  statusClique[kClique]=-2;
	  if (printit>1)
	    printf("clique %d is subset of %d\n",kClique,iClique);
	  if (!ppp&&false) {
	    ppp=true;
	    printf("Kclique %d ",kClique);
	    for (int j=cliqueStart[kClique];j<cliqueStart[kClique+1];j++) {
	      int kColumn = entry[j];
	      printf("%d ",kColumn);
	    }
	    printf("\n");
	  }
	  if (false) {
	    printf("Iclique %d ",iClique);
	    for (int j=cliqueStart[iClique];j<cliqueStart[iClique+1];j++) {
	      int iColumn = entry[j];
	      printf("%d ",iColumn);
	    }
	    printf("\n");
	  }
	  nOut++;
	  used[iClique]=std::min(used[iClique],kClique);;
	  used[kClique]=std::min(used[kClique],iClique);;
	  // But need to look at type
	  // when might be able to fix variables
	  if (cliqueType[kClique]=='E') {
	    statusClique[kClique]=kSave;
	    //statusClique[iClique]=-2;
	    nOut--;
	    // can fix all in iClique not in kClique
	    //printf("ZZ clique %d E, %d S\n",kClique,iClique);
	    int offset = cliqueStart[iClique]-position[kClique];
	    int j;
	    for (j=cliqueStart[kClique];j<cliqueStart[kClique+1];j++) {
	      int kColumn = entry[j];
	      int iColumn = entry[j+offset];
	      while (iColumn<kColumn) {
		if (!fixed[iColumn]) {
		  //printf("ZZ fixing %d to zero\n",iColumn);
		  fixed[iColumn]=-1;
		} else {
		  assert (fixed[iColumn]==-1);
		}
		offset++;
		if (j+offset<cliqueStart[iClique+1]) {
		  iColumn = entry[j+offset];
		} else {
		  iColumn=numberIntegers;
		}
	      }
	    }
#if 0
	  } else if (cliqueType[iClique]=='E') {
	    printf("ZZ clique %d S, %d E\n",kClique,iClique);
	    statusClique[kClique]=-1;
	    nOut--;
#endif
	  }
	  break;
	}
      }
    }
  }
  if (nOut) {
    if(printit) 
      printf("Can get rid of %d cliques\n",nOut);
  }
  for (int i=0;i<numberCliques;i++) {
    if (statusClique[i]!=-2) {
      statusClique[i]=-1;
    }
  }
  delete [] sort;
  delete [] which;
  delete [] position;
  delete [] whichP;
  return nOut;
}
#endif
void CglDuplicateRow::generateCuts8(const OsiSolverInterface & si, OsiCuts & cs,
				    const CglTreeInfo )
{
  bool printit=false;
  bool feasible=true;
  int numberCliques=0;
  int numberEntries=0;
  int * cliqueStart = NULL;
  int * entry = NULL;
  char * cliqueType=NULL;
  int numberRows=si.getNumRows(); 
  const CoinPackedMatrix * rowCopy = si.getMatrixByRow();
  assert(numberRows&&si.getNumCols());
  int iRow;
  const int * column = rowCopy->getIndices();
  const double * elementByRow = rowCopy->getElements();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * lower = si.getColLower();
  const double * upper = si.getColUpper();
  const double * rowLower = si.getRowLower();
  const double * rowUpper = si.getRowUpper();
  // Find 0-1 variables
  int numberIntegers=0;
  int numberColumns=si.getNumCols();
  int * backward = new int [numberColumns];
  for (int i=0;i<numberColumns;i++) {
    if (lower[i]==0.0&&upper[i]==1.0) 
      backward[i]=numberIntegers++;
    else
      backward[i]=-1;
  }
  int * whichP = new int [numberIntegers];
  for (int iPass=0;iPass<2;iPass++) {
    if (iPass) {
      cliqueStart = new int [numberCliques+1];
      cliqueStart[0]=0;
      entry = new int [numberEntries];
      cliqueType = new char [numberCliques];
      numberCliques=0;
      numberEntries=0;
    }
    for (iRow=0;iRow<numberRows;iRow++) {
      duplicate_[iRow]=-1;
      int numberP1=0;
      //int numberTotal=0;
      CoinBigIndex j;
      double upperValue=rowUpper[iRow];
      double lowerValue=rowLower[iRow];
      bool good=true;
      for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
	int iColumn = column[j];
	double value = elementByRow[j];
	if (upper[iColumn]-lower[iColumn]<1.0e-8) {
	  // fixed
	  upperValue -= lower[iColumn]*value;
	  lowerValue -= lower[iColumn]*value;
	  continue;
	} else if (backward[iColumn]<0) {
	  good = false;
	  break;
	} else {
	  iColumn = backward[iColumn];
	  //numberTotal++;
	}
	if (value!=1.0) {
	  good=false;
	} else {
	  assert (numberP1<numberIntegers);
	  whichP[numberP1++]=iColumn;;
	}
      }
      int iUpper = upperValue > INT_MAX ? INT_MAX : static_cast<int> (floor(upperValue+1.0e-5));
      int iLower = lowerValue < INT_MIN ? INT_MIN : static_cast<int> (ceil(lowerValue-1.0e-5));
      int state=0;
      if (upperValue<1.0e6) {
	if (iUpper==1)
	  state=1;
	else if (iUpper==0)
	  state=2;
	else if (iUpper<0)
	  state=3;
	if (fabs(static_cast<double> (iUpper)-upperValue)>1.0e-9)
	  state =-1;
      }
      if (!state&&lowerValue>-1.0e6) {
	if (-iLower==1-numberP1)
	  state=-1;
	else if (-iLower==-numberP1)
	  state=-2;
	else if (-iLower<-numberP1)
	  state=-3;
	if (fabs(static_cast<double> (iLower)-lowerValue)>1.0e-9)
	  state =-1;
      }
      if (numberP1<2)
	state=-1;
      if (good&&state>0) {
	if (abs(state)==3) {
	  // infeasible
	  printf("FFF Infeasible\n");;
	  feasible=false;
	  break;
	} else if (abs(state)==2) {
	  // we can fix all
	  //numberFixed += numberP1+numberM1;
	  printf("FFF can fix %d\n",numberP1);
	} else {
	  for (j=0;j<numberP1;j++) {
	    int iColumn = whichP[j];
	    if (iPass) {
	      entry[numberEntries]=iColumn;
	    }
	    numberEntries++;
	  }
	  if (iPass) {
	    if (iLower!=iUpper) {
	      // slack
	      cliqueType[numberCliques]='S';
	    } else {
	      cliqueType[numberCliques]='E';
	    }
	    cliqueStart[numberCliques+1]=numberEntries;
	    duplicate_[iRow]=numberCliques;
	  }
	  numberCliques++;
	}
      }
    }
  }
  delete[] whichP;
  int * dups = new int [2*numberCliques];
  int * fixed = new int[std::max(numberIntegers,numberCliques)];
  memset(fixed,0,numberIntegers*sizeof(int));
  outDupsEtc(numberIntegers, numberCliques, dups,
	     cliqueStart, cliqueType, entry, fixed, printit ? 2 : 0);
  delete[] cliqueStart;
  delete[] entry;
  delete[] cliqueType;
  int nFixed=0;
  CoinPackedVector ubs;
  for (int i=0;i<numberColumns;i++) {
    int i01=backward[i];
    if (i01>=0&&fixed[i01]==-1) {
      ubs.insert(i,0.0);
      nFixed++;
    }
  }
  for (iRow=0;iRow<numberRows;iRow++) {
    int iClique = duplicate_[iRow];
    if (iClique>=0) 
      fixed[iClique]=iRow;
  }
  int * dup2 = new int [2*numberRows];
  int * used2 = dup2+numberRows;
  int * used = dups+numberCliques;
  for (iRow=0;iRow<numberRows;iRow++) {
    dup2[iRow]=-3; // say not clique
    used2[iRow]=-1; // say not used
  }
  for (iRow=0;iRow<numberRows;iRow++) {
    int iClique = duplicate_[iRow];
    if (iClique>=0) {
      dup2[iRow] = dups[iClique];
      int which = used[iClique];
      if (which>=0&&which<numberCliques)
	which = fixed[which];
      used2[iRow]=which;
    }
  }
  delete [] duplicate_;
  duplicate_=dup2;
  delete [] fixed;
  delete [] dups;
  delete [] backward;
  if (nFixed) {
    OsiColCut cc;
    cc.setUbs(ubs);
    cc.setEffectiveness(100.0);
    cs.insert(cc);
  }
  if (!feasible) {
    // generate infeasible cut and return
    printf("QQ**** infeasible cut\n");
    OsiRowCut rc;
    rc.setLb(COIN_DBL_MAX);
    rc.setUb(0.0);   
    cs.insertIfNotDuplicate(rc);
  }
}
//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglDuplicateRow::CglDuplicateRow ()
:
  CglCutGenerator(),
  rhs_(NULL),
  duplicate_(NULL),
  lower_(NULL),
  storedCuts_(NULL),
  maximumDominated_(1000),
  maximumRhs_(1),
  sizeDynamic_(COIN_INT_MAX),
  mode_(3),
  logLevel_(0)
{
}
// Useful constructor
CglDuplicateRow::CglDuplicateRow(OsiSolverInterface * solver)
  : CglCutGenerator(),
    rhs_(NULL),
    duplicate_(NULL),
    lower_(NULL),
    storedCuts_(NULL),
    maximumDominated_(1000),
    maximumRhs_(1),
    sizeDynamic_(COIN_INT_MAX),
    mode_(3),
    logLevel_(0)
{
  refreshSolver(solver);
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglDuplicateRow::CglDuplicateRow (  const CglDuplicateRow & rhs)
                                                              :
  CglCutGenerator(rhs),
  matrix_(rhs.matrix_),
  matrixByRow_(rhs.matrixByRow_),
  storedCuts_(NULL),
  maximumDominated_(rhs.maximumDominated_),
  maximumRhs_(rhs.maximumRhs_),
  sizeDynamic_(rhs.sizeDynamic_),
  mode_(rhs.mode_),
  logLevel_(rhs.logLevel_)
{  
  int numberRows=matrix_.getNumRows();
  rhs_ = CoinCopyOfArray(rhs.rhs_,numberRows);
  duplicate_ = CoinCopyOfArray(rhs.duplicate_,numberRows);
  lower_ = CoinCopyOfArray(rhs.lower_,numberRows);
  if (rhs.storedCuts_)
    storedCuts_ = new CglStored(*rhs.storedCuts_);
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglDuplicateRow::clone() const
{
  return new CglDuplicateRow(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglDuplicateRow::~CglDuplicateRow ()
{
  // free memory
  delete [] rhs_;
  delete [] duplicate_;
  delete [] lower_;
  delete storedCuts_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglDuplicateRow &
CglDuplicateRow::operator=(
                                         const CglDuplicateRow& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    delete [] rhs_;
    delete [] duplicate_;
    delete [] lower_;
    delete storedCuts_;
    storedCuts_ = NULL;
    matrix_=rhs.matrix_;
    matrixByRow_=rhs.matrixByRow_;
    maximumDominated_ = rhs.maximumDominated_;
    maximumRhs_=rhs.maximumRhs_;
    sizeDynamic_ = rhs.sizeDynamic_;
    mode_ = rhs.mode_;
    logLevel_ = rhs.logLevel_;
    int numberRows=matrix_.getNumRows();
    rhs_ = CoinCopyOfArray(rhs.rhs_,numberRows);
    duplicate_ = CoinCopyOfArray(rhs.duplicate_,numberRows);
    lower_ = CoinCopyOfArray(rhs.lower_,numberRows);
  if (rhs.storedCuts_)
    storedCuts_ = new CglStored(*rhs.storedCuts_);
  }
  return *this;
}

// This can be used to refresh any information
void 
CglDuplicateRow::refreshSolver(OsiSolverInterface * solver)
{
  delete [] rhs_;
  delete [] duplicate_;
  delete [] lower_;
  matrix_ = *solver->getMatrixByCol();
  matrix_.removeGaps();
  matrix_.orderMatrix();
  matrixByRow_ = *solver->getMatrixByRow();
  int numberRows=matrix_.getNumRows();
  rhs_ = new int[numberRows];
  duplicate_ = new int[numberRows];
  lower_ = new int[numberRows];
  const double * columnLower = solver->getColLower();
  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();
  // Row copy
  const double * elementByRow = matrixByRow_.getElements();
  const int * column = matrixByRow_.getIndices();
  const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
  const int * rowLength = matrixByRow_.getVectorLengths();
  int iRow;
  const char * intVar = solver->getColType(true);
  //int numberGood=0;
  int markBad = -(solver->getNumCols()+1);
  for (iRow=0;iRow<numberRows;iRow++) {
    rhs_[iRow]=markBad;
    lower_[iRow]=markBad;
    duplicate_[iRow]=-1;
    if (rowUpper[iRow]<100) {
      int iRhs= static_cast<int> (floor(rowUpper[iRow]));
      // check elements
      bool good=true;
      for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
        int iColumn = column[j];
        if (!intVar[iColumn])
	  good=false;
        double value = elementByRow[j];
        if (floor(value)!=value||value<1.0) {
          good=false;
        }
      }
      if (good) {
        lower_[iRow] = static_cast<int> (std::max(0.0,ceil(rowLower[iRow])));
        if (iRhs>=lower_[iRow]) {
          rhs_[iRow]=iRhs;
          //numberGood++;
        } else {
          // infeasible ?
          lower_[iRow]=markBad;
          rhs_[iRow]=markBad;
        }
      } else {
        lower_[iRow]=markBad;
        rhs_[iRow]=markBad;
      }
    } else if (rowUpper[iRow]>1.0e30&&rowLower[iRow]==1.0) {
      // may be OK to look for dominated in >=1 rows
      // check elements
      bool good=true;
      for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
        int iColumn = column[j];
        if (!intVar[iColumn])
	  good=false;
        double value = elementByRow[j];
        if (floor(value)!=value||value<1.0) {
          good=false;
        }
	if (columnLower[iColumn]!=0.0)
	  good=false;
      }
      if (good) {
        lower_[iRow] = 1;
      }
    }
  }
}
  /** Fix variables and find duplicate/dominated rows for the model of the 
      solver interface, si.

      This is a very simple minded idea but I (JJF) am using it in a project so thought
      I might as well add it.  It should really be called before first solve and I may
      modify CBC to allow for that.

      This is designed for problems with few rows and many integer variables where the rhs
      are <= or == and all coefficients and rhs are small integers.

      If effective rhs is K then we can fix all variables with coefficients > K to their lower bounds
      (effective rhs just means original with variables with nonzero lower bounds subtracted out).

      If one row is a subset of another and the effective rhs are same we can fix some variables
      and then the two rows are identical.

      This version does deletions and fixings and may return stored cuts for
      dominated columns 
  */
CglStored * 
CglDuplicateRow::outDuplicates( OsiSolverInterface * solver)
{
  
  CglTreeInfo info;
  info.level = 0;
  info.pass = 0;
  int numberRows = solver->getNumRows();
  info.formulation_rows = numberRows;
  info.inTree = false;
  info.strengthenRow= NULL;
  info.pass = 0;
  OsiCuts cs;
  generateCuts(*solver,cs,info);
  // Get rid of duplicate rows
  int * which = new int[numberRows]; 
  int numberDrop=0;
  for (int iRow=0;iRow<numberRows;iRow++) {
    if (duplicate_[iRow]==-2||duplicate_[iRow]>=0) 
      which[numberDrop++]=iRow;
  }
  if (numberDrop) {
    solver->deleteRows(numberDrop,which);
  }
  delete [] which;
  // see if we have any column cuts
  int numberColumnCuts = cs.sizeColCuts() ;
  const double * columnLower = solver->getColLower();
  const double * columnUpper = solver->getColUpper();
  for (int k = 0;k<numberColumnCuts;k++) {
    OsiColCut * thisCut = cs.colCutPtr(k) ;
    const CoinPackedVector & lbs = thisCut->lbs() ;
    const CoinPackedVector & ubs = thisCut->ubs() ;
    int j ;
    int n ;
    const int * which ;
    const double * values ;
    n = lbs.getNumElements() ;
    which = lbs.getIndices() ;
    values = lbs.getElements() ;
    for (j = 0;j<n;j++) {
      int iColumn = which[j] ;
      if (values[j]>columnLower[iColumn]) 
        solver->setColLower(iColumn,values[j]) ;
    }
    n = ubs.getNumElements() ;
    which = ubs.getIndices() ;
    values = ubs.getElements() ;
    for (j = 0;j<n;j++) {
      int iColumn = which[j] ;
      if (values[j]<columnUpper[iColumn]) 
        solver->setColUpper(iColumn,values[j]) ;
    }
  }
  return storedCuts_;
}
// Create C++ lines to get to current state
std::string
CglDuplicateRow::generateCpp( FILE * fp) 
{
  CglDuplicateRow other;
  fprintf(fp,"0#include \"CglDuplicateRow.hpp\"\n");
  fprintf(fp,"3  CglDuplicateRow duplicateRow;\n");
  if (logLevel_!=other.logLevel_)
    fprintf(fp,"3  duplicateRow.setLogLevel(%d);\n",logLevel_);
  else
    fprintf(fp,"4  duplicateRow.setLogLevel(%d);\n",logLevel_);
  if (maximumRhs_!=other.maximumRhs_)
    fprintf(fp,"3  duplicateRow.setMaximumRhs(%d);\n",maximumRhs_);
  else
    fprintf(fp,"4  duplicateRow.setMaximumRhs(%d);\n",maximumRhs_);
  if (maximumDominated_!=other.maximumDominated_)
    fprintf(fp,"3  duplicateRow.setMaximumDominated(%d);\n",maximumDominated_);
  else
    fprintf(fp,"4  duplicateRow.setMaximumDominated(%d);\n",maximumDominated_);
  if (mode_!=other.mode_)
    fprintf(fp,"3  duplicateRow.setMode(%d);\n",mode_);
  else
    fprintf(fp,"4  duplicateRow.setMode(%d);\n",mode_);
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  duplicateRow.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  duplicateRow.setAggressiveness(%d);\n",getAggressiveness());
  return "duplicateRow";
}
