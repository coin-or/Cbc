// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpFactorization.hpp"
#include "ClpSimplex.hpp"
#include "ClpAmplObjective.hpp"
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpAmplObjective::ClpAmplObjective () 
: ClpObjective()
{
  type_=12;
  objective_=NULL;
  amplObjective_=NULL;
  gradient_ = NULL;
}

//-------------------------------------------------------------------
// Useful Constructor 
//-------------------------------------------------------------------
ClpAmplObjective::ClpAmplObjective (void * amplInfo)
  : ClpObjective()
{
  type_=12;
  loadAmplObjective(amplInfo);
  numberColumns_ = numberColumns;
  if (numberExtendedColumns>=0)
    numberExtendedColumns_= max(numberColumns_,numberExtendedColumns);
  else
    numberExtendedColumns_= numberColumns_;
  if (objective) {
    objective_ = new double [numberExtendedColumns_];
    memcpy(objective_,objective,numberColumns_*sizeof(double));
    memset(objective_+numberColumns_,0,(numberExtendedColumns_-numberColumns_)*sizeof(double));
  } else {
    objective_ = new double [numberExtendedColumns_];
    memset(objective_,0,numberExtendedColumns_*sizeof(double));
  }
  if (start) 
    amplObjective_ = new CoinPackedMatrix(true,numberColumns,numberColumns,
					     start[numberColumns],element,column,start,NULL);
  else
  amplObjective_=NULL;
  gradient_ = NULL;
  activated_=1;
  fullMatrix_=false;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpAmplObjective::ClpAmplObjective (const ClpAmplObjective & rhs) 
: ClpObjective(rhs)
{  
  numberColumns_=rhs.numberColumns_;
  numberExtendedColumns_=rhs.numberExtendedColumns_;
  fullMatrix_=rhs.fullMatrix_;
  if (rhs.objective_) {
    objective_ = new double [numberExtendedColumns_];
    memcpy(objective_,rhs.objective_,numberExtendedColumns_*sizeof(double));
  } else {
    objective_=NULL;
  }
  if (rhs.gradient_) {
    gradient_ = new double [numberExtendedColumns_];
    memcpy(gradient_,rhs.gradient_,numberExtendedColumns_*sizeof(double));
  } else {
    gradient_=NULL;
  }
  if (rhs.amplObjective_) {
    // see what type of matrix wanted
    if (type==0) {
      // just copy
      amplObjective_ = new CoinPackedMatrix(*rhs.amplObjective_);
    } else if (type==1) {
      // expand to full symmetric
      fullMatrix_=true;
      const int * columnAmpl1 = rhs.amplObjective_->getIndices();
      const CoinBigIndex * columnAmplStart1 = rhs.amplObjective_->getVectorStarts();
      const int * columnAmplLength1 = rhs.amplObjective_->getVectorLengths();
      const double * amplElement1 = rhs.amplObjective_->getElements();
      CoinBigIndex * columnAmplStart2 = new CoinBigIndex [numberExtendedColumns_+1];
      int * columnAmplLength2 = new int [numberExtendedColumns_];
      int iColumn;
      int numberColumns = rhs.amplObjective_->getNumCols();
      int numberBelow=0;
      int numberAbove=0;
      int numberDiagonal=0;
      CoinZeroN(columnAmplLength2,numberExtendedColumns_);
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	for (CoinBigIndex j=columnAmplStart1[iColumn];
	     j<columnAmplStart1[iColumn]+columnAmplLength1[iColumn];j++) {
	  int jColumn = columnAmpl1[j];
	  if (jColumn>iColumn) {
	    numberBelow++;
	    columnAmplLength2[jColumn]++;
	    columnAmplLength2[iColumn]++;
	  } else if (jColumn==iColumn) {
	    numberDiagonal++;
	    columnAmplLength2[iColumn]++;
	  } else {
	    numberAbove++;
	  }
	}
      }
      if (numberAbove>0) {
	if (numberAbove==numberBelow) {
	  // already done
	  amplObjective_ = new CoinPackedMatrix(*rhs.amplObjective_);
	  delete [] columnAmplStart2;
	  delete [] columnAmplLength2;
	} else {
	  printf("number above = %d, number below = %d, error\n",
		 numberAbove,numberBelow);
	  abort();
	}
      } else {
	int numberElements=numberDiagonal+2*numberBelow;
	int * columnAmpl2 = new int [numberElements];
	double * amplElement2 = new double [numberElements];
	columnAmplStart2[0]=0;
	numberElements=0;
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  int n=columnAmplLength2[iColumn];
	  columnAmplLength2[iColumn]=0;
	  numberElements += n;
	  columnAmplStart2[iColumn+1]=numberElements;
	}
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  for (CoinBigIndex j=columnAmplStart1[iColumn];
	       j<columnAmplStart1[iColumn]+columnAmplLength1[iColumn];j++) {
	    int jColumn = columnAmpl1[j];
	    if (jColumn>iColumn) {
	      // put in two places
	      CoinBigIndex put=columnAmplLength2[jColumn]+columnAmplStart2[jColumn];
	      columnAmplLength2[jColumn]++;
	      amplElement2[put]=amplElement1[j];
	      columnAmpl2[put]=iColumn;
	      put=columnAmplLength2[iColumn]+columnAmplStart2[iColumn];
	      columnAmplLength2[iColumn]++;
	      amplElement2[put]=amplElement1[j];
	      columnAmpl2[put]=jColumn;
	    } else if (jColumn==iColumn) {
	      CoinBigIndex put=columnAmplLength2[iColumn]+columnAmplStart2[iColumn];
	      columnAmplLength2[iColumn]++;
	      amplElement2[put]=amplElement1[j];
	      columnAmpl2[put]=iColumn;
	    } else {
	      abort();
	    }
	  }
	}
	// Now create
	amplObjective_ = 
	  new CoinPackedMatrix (true,
				rhs.numberExtendedColumns_,
				rhs.numberExtendedColumns_,
				numberElements,
				amplElement2,
				columnAmpl2,
				columnAmplStart2,
				columnAmplLength2,0.0,0.0);
	delete [] columnAmplStart2;
	delete [] columnAmplLength2;
	delete [] columnAmpl2;
	delete [] amplElement2;
      }
    } else {
      fullMatrix_=false;
      abort(); // code when needed
    }
	    
  } else {
    amplObjective_=NULL;
  }
}
  

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpAmplObjective::~ClpAmplObjective ()
{
  delete [] objective_;
  delete [] gradient_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpAmplObjective &
ClpAmplObjective::operator=(const ClpAmplObjective& rhs)
{
  if (this != &rhs) {
    amplObjective_ = NULL;
    ClpObjective::operator=(rhs);
    numberColumns_=rhs.numberColumns_;
    numberExtendedColumns_=rhs.numberExtendedColumns_;
    if (rhs.objective_) {
      objective_ = new double [numberExtendedColumns_];
      memcpy(objective_,rhs.objective_,numberExtendedColumns_*sizeof(double));
    } else {
      objective_=NULL;
    }
    if (rhs.gradient_) {
      gradient_ = new double [numberExtendedColumns_];
      memcpy(gradient_,rhs.gradient_,numberExtendedColumns_*sizeof(double));
    } else {
      gradient_=NULL;
    }
    if (rhs.amplObjective_) {
      amplObjective_ = new CoinPackedMatrix(*rhs.amplObjective_);
    } else {
      amplObjective_=NULL;
    }
  }
  return *this;
}

// Returns gradient
double *  
ClpAmplObjective::gradient(const ClpSimplex * model,
				const double * solution, double & offset,bool refresh,
				int includeLinear)
{
  offset=0.0;
  bool scaling=false;
  if (model&&(model->rowScale()||
	      model->objectiveScale()!=1.0||model->optimizationDirection()!=1.0))
    scaling=true;
  const double * cost = NULL;
  if (model)
    cost = model->costRegion();
  if (!cost) {
    // not in solve
    cost = objective_;
    scaling=false;
  }
  if (!scaling) {
    if (!amplObjective_||!solution||!activated_) {
      return objective_;
    } else {
      if (refresh||!gradient_) {
	if (!gradient_) 
	  gradient_ = new double[numberExtendedColumns_];
	const int * columnAmpl = amplObjective_->getIndices();
	const CoinBigIndex * columnAmplStart = amplObjective_->getVectorStarts();
	const int * columnAmplLength = amplObjective_->getVectorLengths();
	const double * amplElement = amplObjective_->getElements();
	offset=0.0;
	// use current linear cost region 
	if (includeLinear==1)
	  memcpy(gradient_,cost,numberExtendedColumns_*sizeof(double));
	else if (includeLinear==2)
	  memcpy(gradient_,objective_,numberExtendedColumns_*sizeof(double));
	else
	  memset(gradient_,0,numberExtendedColumns_*sizeof(double));
	if (activated_) {
	  if (!fullMatrix_) {
	    int iColumn;
	    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	      double valueI = solution[iColumn];
	      CoinBigIndex j;
	      for (j=columnAmplStart[iColumn];
		   j<columnAmplStart[iColumn]+columnAmplLength[iColumn];j++) {
		int jColumn = columnAmpl[j];
		double valueJ = solution[jColumn];
		double elementValue = amplElement[j];
		if (iColumn!=jColumn) {
		  offset += valueI*valueJ*elementValue;
		  //if (fabs(valueI*valueJ*elementValue)>1.0e-12)
                  //printf("%d %d %g %g %g -> %g\n",
                  //       iColumn,jColumn,valueI,valueJ,elementValue,
                  //       valueI*valueJ*elementValue);
		  double gradientI = valueJ*elementValue;
		  double gradientJ = valueI*elementValue;
		  gradient_[iColumn] += gradientI;
		  gradient_[jColumn] += gradientJ;
		} else {
		  offset += 0.5*valueI*valueI*elementValue;
		  //if (fabs(valueI*valueI*elementValue)>1.0e-12)
                  //printf("XX %d %g %g -> %g\n",
                  //       iColumn,valueI,elementValue,
                  //       0.5*valueI*valueI*elementValue);
		  double gradientI = valueI*elementValue;
		  gradient_[iColumn] += gradientI;
		}
	      }
	    }
	  } else {
	    // full matrix
	    int iColumn;
	    offset *= 2.0;
	    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	      CoinBigIndex j;
	      double value=0.0;
	      double current = gradient_[iColumn];
	      for (j=columnAmplStart[iColumn];
		   j<columnAmplStart[iColumn]+columnAmplLength[iColumn];j++) {
		int jColumn = columnAmpl[j];
		double valueJ = solution[jColumn]*amplElement[j];
		value += valueJ;
	      }
	      offset += value*solution[iColumn];
	      gradient_[iColumn] = current+value;
	    }
	    offset *= 0.5;
	  }
	}
      }
      if (model)
        offset *= model->optimizationDirection()*model->objectiveScale();
      return gradient_;
    }
  } else {
    // do scaling
    assert (solution);
    // for now only if half
    assert (!fullMatrix_);
    if (refresh||!gradient_) {
      if (!gradient_) 
	gradient_ = new double[numberExtendedColumns_];
      double direction = model->optimizationDirection()*model->objectiveScale();
      // direction is actually scale out not scale in
      //if (direction)
      //direction = 1.0/direction;
      const int * columnAmpl = amplObjective_->getIndices();
      const CoinBigIndex * columnAmplStart = amplObjective_->getVectorStarts();
      const int * columnAmplLength = amplObjective_->getVectorLengths();
      const double * amplElement = amplObjective_->getElements();
      int iColumn;
      const double * columnScale = model->columnScale();
      // use current linear cost region (already scaled)
      if (includeLinear==1) {
	memcpy(gradient_,model->costRegion(),numberExtendedColumns_*sizeof(double));
      }	else if (includeLinear==2) {
	memset(gradient_+numberColumns_,0,(numberExtendedColumns_-numberColumns_)*sizeof(double));
	if (!columnScale) {
	  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	    gradient_[iColumn]= objective_[iColumn]*direction;
	  }
	} else {
	  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	    gradient_[iColumn]= objective_[iColumn]*direction*columnScale[iColumn];
	  }
	}
      } else {
	memset(gradient_,0,numberExtendedColumns_*sizeof(double));
      }
      if (!columnScale) {
	if (activated_) {
	  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	    double valueI = solution[iColumn];
	    CoinBigIndex j;
	    for (j=columnAmplStart[iColumn];
		 j<columnAmplStart[iColumn]+columnAmplLength[iColumn];j++) {
	      int jColumn = columnAmpl[j];
	      double valueJ = solution[jColumn];
	      double elementValue = amplElement[j];
	      elementValue *= direction;
	      if (iColumn!=jColumn) {
		offset += valueI*valueJ*elementValue;
		double gradientI = valueJ*elementValue;
		double gradientJ = valueI*elementValue;
		gradient_[iColumn] += gradientI;
		gradient_[jColumn] += gradientJ;
	      } else {
		offset += 0.5*valueI*valueI*elementValue;
		double gradientI = valueI*elementValue;
		gradient_[iColumn] += gradientI;
	      }
	    }
	  }
	}
      } else {
	// scaling
	if (activated_) {
	  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	    double valueI = solution[iColumn];
	    double scaleI = columnScale[iColumn]*direction;
	    CoinBigIndex j;
	    for (j=columnAmplStart[iColumn];
		 j<columnAmplStart[iColumn]+columnAmplLength[iColumn];j++) {
	      int jColumn = columnAmpl[j];
	      double valueJ = solution[jColumn];
	      double elementValue = amplElement[j];
	      double scaleJ = columnScale[jColumn];
	      elementValue *= scaleI*scaleJ;
	      if (iColumn!=jColumn) {
		offset += valueI*valueJ*elementValue;
		double gradientI = valueJ*elementValue;
		double gradientJ = valueI*elementValue;
		gradient_[iColumn] += gradientI;
		gradient_[jColumn] += gradientJ;
	      } else {
		offset += 0.5*valueI*valueI*elementValue;
		double gradientI = valueI*elementValue;
		gradient_[iColumn] += gradientI;
	      }
	    }
	  }
	}
      }
    }
    if (model)
      offset *= model->optimizationDirection();
    return gradient_;
  }
}
  
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpObjective * ClpAmplObjective::clone() const
{
  return new ClpAmplObjective(*this);
}
// Resize objective
void 
ClpAmplObjective::resize(int newNumberColumns)
{
  if (numberColumns_!=newNumberColumns) {
    abort();
  } 
  
}
// Delete columns in  objective
void 
ClpAmplObjective::deleteSome(int numberToDelete, const int * which) 
{
  if (numberToDelete)
    abort();
}

// Load up ampl objective 
void 
ClpAmplObjective::loadAmplObjective(void * amplInfo)
{
  fullMatrix_=false;
  delete amplObjective_;
  amplObjective_ = new CoinPackedMatrix(true,numberColumns,numberColumns,
					     start[numberColumns],element,column,start,NULL);
  numberColumns_=numberColumns;
  if (numberExtended>numberExtendedColumns_) {
    if (objective_) {
      // make correct size
      double * newArray = new double[numberExtended];
      memcpy(newArray,objective_,numberColumns_*sizeof(double));
      delete [] objective_;
      objective_ = newArray;
      memset(objective_+numberColumns_,0,(numberExtended-numberColumns_)*sizeof(double));
    }
    if (gradient_) {
      // make correct size
      double * newArray = new double[numberExtended];
      memcpy(newArray,gradient_,numberColumns_*sizeof(double));
      delete [] gradient_;
      gradient_ = newArray;
      memset(gradient_+numberColumns_,0,(numberExtended-numberColumns_)*sizeof(double));
    }
    numberExtendedColumns_ = numberExtended;
  } else {
    numberExtendedColumns_ = numberColumns_;
  }
}
/* Returns reduced gradient.Returns an offset (to be added to current one).
 */
double 
ClpAmplObjective::reducedGradient(ClpSimplex * model, double * region,
				       bool useFeasibleCosts)
{
  int numberRows = model->numberRows();
  int numberColumns=model->numberColumns();
  
  //work space
  CoinIndexedVector  * workSpace = model->rowArray(0);
  
  CoinIndexedVector arrayVector;
  arrayVector.reserve(numberRows+1);
  
  int iRow;
#ifdef CLP_DEBUG
  workSpace->checkClear();
#endif
  double * array = arrayVector.denseVector();
  int * index = arrayVector.getIndices();
  int number=0;
  const double * costNow = gradient(model,model->solutionRegion(),offset_,
				    true,useFeasibleCosts ? 2 : 1);
  double * cost = model->costRegion();
  const int * pivotVariable = model->pivotVariable();
  for (iRow=0;iRow<numberRows;iRow++) {
    int iPivot=pivotVariable[iRow];
    double value;
    if (iPivot<numberColumns)
      value = costNow[iPivot];
    else if (!useFeasibleCosts) 
      value = cost[iPivot];
    else 
      value=0.0;
    if (value) {
      array[iRow]=value;
      index[number++]=iRow;
    }
  }
  arrayVector.setNumElements(number);

  // Btran basic costs
  model->factorization()->updateColumnTranspose(workSpace,&arrayVector);
  double * work = workSpace->denseVector();
  ClpFillN(work,numberRows,0.0);
  // now look at dual solution
  double * rowReducedCost = region+numberColumns;
  double * dual = rowReducedCost;
  const double * rowCost = cost+numberColumns;
  for (iRow=0;iRow<numberRows;iRow++) {
    dual[iRow]=array[iRow];
  }
  double * dj = region;
  ClpDisjointCopyN(costNow,numberColumns,dj);
  
  model->transposeTimes(-1.0,dual,dj);
  for (iRow=0;iRow<numberRows;iRow++) {
    // slack
    double value = dual[iRow];
    value += rowCost[iRow];
    rowReducedCost[iRow]=value;
  }
  return offset_;
}
/* Returns step length which gives minimum of objective for
   solution + theta * change vector up to maximum theta.
   
   arrays are numberColumns+numberRows
*/
double 
ClpAmplObjective::stepLength(ClpSimplex * model,
				  const double * solution,
				  const double * change,
				  double maximumTheta,
				  double & currentObj,
				  double & predictedObj,
				  double & thetaObj)
{
  const double * cost = model->costRegion();
  bool inSolve=true;
  if (!cost) {
    // not in solve
    cost = objective_;
    inSolve=false;
  }
  double delta=0.0;
  double linearCost =0.0;
  int numberRows = model->numberRows();
  int numberColumns = model->numberColumns();
  int numberTotal = numberColumns;
  if (inSolve)
    numberTotal += numberRows;
  currentObj=0.0;
  thetaObj=0.0;
  for (int iColumn=0;iColumn<numberTotal;iColumn++) {
    delta += cost[iColumn]*change[iColumn];
    linearCost += cost[iColumn]*solution[iColumn];
  }
  if (!activated_||!amplObjective_) {
    currentObj=linearCost;
    thetaObj =currentObj + delta*maximumTheta;
    if (delta<0.0) {
      return maximumTheta;
    } else {
      printf("odd linear direction %g\n",delta);
      return 0.0;
    }
  }
  assert (model);
  bool scaling=false;
  if ((model->rowScale()||
       model->objectiveScale()!=1.0||model->optimizationDirection()!=1.0)&&inSolve)
    scaling=true;
  const int * columnAmpl = amplObjective_->getIndices();
  const CoinBigIndex * columnAmplStart = amplObjective_->getVectorStarts();
  const int * columnAmplLength = amplObjective_->getVectorLengths();
  const double * amplElement = amplObjective_->getElements();
  double a=0.0;
  double b=delta;
  double c=0.0;
  if (!scaling) {
    if (!fullMatrix_) {
      int iColumn;
      for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	double valueI = solution[iColumn];
	double changeI = change[iColumn];
	CoinBigIndex j;
	for (j=columnAmplStart[iColumn];
	     j<columnAmplStart[iColumn]+columnAmplLength[iColumn];j++) {
	  int jColumn = columnAmpl[j];
	  double valueJ = solution[jColumn];
	  double changeJ = change[jColumn];
	  double elementValue = amplElement[j];
	  if (iColumn!=jColumn) {
	    a += changeI*changeJ*elementValue;
	    b += (changeI*valueJ+changeJ*valueI)*elementValue;
	    c += valueI*valueJ*elementValue;
	  } else {
	    a += 0.5*changeI*changeI*elementValue;
	    b += changeI*valueI*elementValue;
	    c += 0.5*valueI*valueI*elementValue;
	  }
	}
      }
    } else {
      // full matrix stored
      int iColumn;
      for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	double valueI = solution[iColumn];
	double changeI = change[iColumn];
	CoinBigIndex j;
	for (j=columnAmplStart[iColumn];
	     j<columnAmplStart[iColumn]+columnAmplLength[iColumn];j++) {
	  int jColumn = columnAmpl[j];
	  double valueJ = solution[jColumn];
	  double changeJ = change[jColumn];
	  double elementValue = amplElement[j];
	  valueJ *= elementValue;
	  a += changeI*changeJ*elementValue;
	  b += changeI*valueJ;
	  c += valueI*valueJ;
	}
      }
      a *= 0.5;
      c *= 0.5;
    }
  } else {
    // scaling
    // for now only if half
    assert (!fullMatrix_);
    const double * columnScale = model->columnScale();
    double direction = model->optimizationDirection()*model->objectiveScale();
    // direction is actually scale out not scale in
    if (direction)
      direction = 1.0/direction;
    if (!columnScale) {
      for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
	double valueI = solution[iColumn];
	double changeI = change[iColumn];
	CoinBigIndex j;
	for (j=columnAmplStart[iColumn];
	     j<columnAmplStart[iColumn]+columnAmplLength[iColumn];j++) {
	  int jColumn = columnAmpl[j];
	  double valueJ = solution[jColumn];
	  double changeJ = change[jColumn];
	  double elementValue = amplElement[j];
	  elementValue *= direction;
	  if (iColumn!=jColumn) {
	    a += changeI*changeJ*elementValue;
	    b += (changeI*valueJ+changeJ*valueI)*elementValue;
	    c += valueI*valueJ*elementValue;
	  } else {
	    a += 0.5*changeI*changeI*elementValue;
	    b += changeI*valueI*elementValue;
	    c += 0.5*valueI*valueI*elementValue;
	  }
	}
      }
    } else {
      // scaling
      for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
	double valueI = solution[iColumn];
	double changeI = change[iColumn];
	double scaleI = columnScale[iColumn]*direction;
	CoinBigIndex j;
	for (j=columnAmplStart[iColumn];
	     j<columnAmplStart[iColumn]+columnAmplLength[iColumn];j++) {
	  int jColumn = columnAmpl[j];
	  double valueJ = solution[jColumn];
	  double changeJ = change[jColumn];
	  double elementValue = amplElement[j];
	  elementValue *= scaleI*columnScale[jColumn];
	  if (iColumn!=jColumn) {
	    a += changeI*changeJ*elementValue;
	    b += (changeI*valueJ+changeJ*valueI)*elementValue;
	    c += valueI*valueJ*elementValue;
	  } else {
	    a += 0.5*changeI*changeI*elementValue;
	    b += changeI*valueI*elementValue;
	    c += 0.5*valueI*valueI*elementValue;
	  }
	}
      }
    }
  }
  double theta;
  //printf("Current cost %g\n",c+linearCost);
  currentObj = c+linearCost;
  thetaObj = currentObj + a*maximumTheta*maximumTheta+b*maximumTheta;
  // minimize a*x*x + b*x + c
  if (a<=0.0) {
    theta = maximumTheta;
  } else {
    theta = -0.5*b/a;
  }
  predictedObj = currentObj + a*theta*theta+b*theta;
  if (b>0.0) {
    if (model->messageHandler()->logLevel()&32)
      printf("a %g b %g c %g => %g\n",a,b,c,theta); 
    b=0.0;
  }
  return CoinMin(theta,maximumTheta);
}
// Scale objective 
void 
ClpAmplObjective::reallyScale(const double * columnScale) 
{
  abort();
}
/* Given a zeroed array sets nonlinear columns to 1.
   Returns number of nonlinear columns
*/
int 
ClpAmplObjective::markNonlinear(char * which)
{
  int iColumn;
  const int * columnAmpl = amplObjective_->getIndices();
  const CoinBigIndex * columnAmplStart = amplObjective_->getVectorStarts();
  const int * columnAmplLength = amplObjective_->getVectorLengths();
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    CoinBigIndex j;
    for (j=columnAmplStart[iColumn];
	 j<columnAmplStart[iColumn]+columnAmplLength[iColumn];j++) {
      int jColumn = columnAmpl[j];
      which[jColumn]=1;
      which[iColumn]=1;
    }
  }
  int numberNonLinearColumns = 0;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    if(which[iColumn])
      numberNonLinearColumns++;
  }
  return numberNonLinearColumns;
}
