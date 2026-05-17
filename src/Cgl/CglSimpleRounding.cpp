// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cfloat> 
#include <cassert>

#include "CoinPragma.hpp"
#include "CglSimpleRounding.hpp" 
#include "CoinPackedVector.hpp"
#include "CoinSort.hpp"
#include "CoinPackedMatrix.hpp"

//-------------------------------------------------------------
void
CglSimpleRounding::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
				const CglTreeInfo /*info*/)
{
  int nRows=si.getNumRows(); // number of rows in the coefficient matrix
  int nCols=si.getNumCols(); // number of columns in the coefficient matrix
  int rowIndex;             // index into the constraint matrix stored in row
                            // order 
  CoinPackedVector irow;     // "integer row": working space to hold the integer
                            // <= inequality derived from the rowIndex-th
                            // constraint 
  double b=0;             // working space for the rhs of integer <= inequality
  bool * negative= new bool[nCols]; // negative[i]= true if coefficient of the 
                                    // ith variable is negative and false
                                    // otherwise 
  int k;                  // dummy iterator variable 
  for ( k=0; k<nCols; k++ ) negative[k] = false;
  
  const CoinPackedMatrix * rowCopy = 
    si.getMatrixByRow(); // row copy: matrix stored in row order

  /////////////////////////////////////////////////////////////////////////////
  // Main loop:                                                              //
  // For every row in the matrix,                                            //
  //     if we can derive a valid <= inequality in integer variables, then   //
  //     try to construct a simple rounding cut from the integer inequality. //
  //     Add the resulting cut to the set of cuts.                           //
  /////////////////////////////////////////////////////////////////////////////

  for (rowIndex=0; rowIndex<nRows; rowIndex++){

    // Only look at tight rows
    // double * pi=ekk_rowduals(model); 
    // if (fabs(pi[row]) < epsilon_){
   //  continue;
    // }

    // Try to derive an <= inequality in integer variables from the row 
    // by netting out the continuous variables.
    // Store the value and the sign of the coefficients separately:
    // irow.getElements() contains the absolute values of the coefficients.
    // negative is a boolean vector indicating the sign of the coeffcients.
    // b is the rhs of the <= integer inequality

    if (!deriveAnIntegerRow( si, 
                             rowIndex, 
                             rowCopy->getVector(rowIndex),
                             irow, b, negative))
    {

      // Reset local data for the next iteration of the rowIndex-loop
      for(k=0; k<irow.getNumElements(); k++) negative[irow.getIndices()[k]]=false;
      irow.setVector(0,NULL,NULL);
      continue;
    } 
 
    // Euclid's greatest common divisor (gcd) algorithm applies to positive
    // INTEGERS. 
    // Determine the power of 10 needed, so that multipylying the integer
    // inequality through by 10**power makes all coefficients essentially
    // integral. 
    int power = power10ToMakeDoubleAnInt(irow.getNumElements(),irow.getElements(),epsilon_*1.0e-4);

    // Now a vector to store the integer-ized values. For instance, 
    // if x[i] is .66 and power is 1000 then xInt[i] will be 660
    int * xInt = NULL;
    if (power >=0) {

      xInt = new int[irow.getNumElements()]; 
      double dxInt; // a double version of xInt for error trapping
      
      
#ifdef CGL_DEBUG      
      printf("The (double) coefficients and their integer-ized counterparts:\n");
#endif

      for (k=0; k<irow.getNumElements(); k++){
	dxInt = irow.getElements()[k]*pow(10.0,power);
	xInt[k]= static_cast<int> (dxInt+0.5); // Need to add the 0.5 
	// so that a dxInt=9.999 will give a xInt=1

#ifdef CGL_DEBUG
	printf("%g     %g   \n",irow.getElements()[k],dxInt);
#endif

      }

    } else {

      // If overflow is detected, one warning message is printed and 
      // the row is skipped.
#ifdef CGL_DEBUG
      printf("SimpleRounding: Warning: Overflow detected \n");
      printf("      on %i of vars in processing row %i. Row skipped.\n",
	     -power, rowIndex);
#endif
      // reset local data for next iteration
      for(k=0; k<irow.getNumElements(); k++) negative[irow.getIndices()[k]]=false;
      irow.setVector(0,NULL,NULL);
      continue;
    }

    // find greatest common divisor of the irow.elements
    int gcd = gcdv(irow.getNumElements(), xInt);

#ifdef CGL_DEBUG
    printf("The gcd of xInt is %i\n",gcd);    
#endif

    // construct new cut by dividing through by gcd and 
    // rounding down rhs and accounting for negatives
    CoinPackedVector cut;
    for (k=0; k<irow.getNumElements(); k++){
        cut.insert(irow.getIndices()[k],xInt[k]/gcd);
    }
    double cutRhs = floor((b*pow(10.0,power))/gcd);

    // un-negate the negated variables in the cut
    {
       const int s = cut.getNumElements();
       const int * indices = cut.getIndices();
       double* elements = cut.getElements();
       for (k=0; k<s; k++){
	 int column=indices[k];
	  if (negative[column]) {
	     elements[k] *= -1;
	  }
       }
    }

    // Create the row cut and add it to the set of cuts
    // It may not be violated
    if (fabs(cutRhs*gcd-b)> epsilon_){ // if the cut and row are different. 
      OsiRowCut rc;
      rc.setRow(cut.getNumElements(),cut.getIndices(),cut.getElements());
      rc.setLb(-COIN_DBL_MAX);
      rc.setUb(cutRhs);   
      cs.insertIfNotDuplicate(rc);

#ifdef CGL_DEBUG
      printf("Row %i had a simple rounding cut:\n",rowIndex);
      printf("Cut size: %i Cut rhs: %g  Index       Element \n",
	     cut.getNumElements(), cutRhs);
      for (k=0; k<cut.getNumElements(); k++){
        printf("%i      %g\n",cut.getIndices()[k], cut.getElements()[k]);
      }
      printf("\n");
#endif
    }

    // Reset local data for the next iteration of the rowIndex-loop
    for(k=0; k<irow.getNumElements(); k++) negative[irow.getIndices()[k]]=false;
    irow.setVector(0,NULL,NULL);
    delete [] xInt;


  }

  delete [] negative;
}


//-------------------------------------------------------------------
// deriveAnIntegerRow:  dervies a <=  inequality
//                  in integer variables of the form ax<=b 
//                  from a row in the model, if possible by
//                  netting out the continuous variables
//-------------------------------------------------------------------
bool
CglSimpleRounding::deriveAnIntegerRow(
       const OsiSolverInterface & si, 
       int rowIndex,
       const CoinShallowPackedVector & matrixRow,
       CoinPackedVector & irow, 
       double & b,
       bool * negative) const
{
  irow.clear();
  int i;           // dummy iterator variable
  double sign=1.0; // +1 if le row, -1 if ge row  

  // number of columns in the row
  int sizeOfRow=matrixRow.getNumElements();

  // Get the sense of the row constraint
  const char  rowsense = si.getRowSense()[rowIndex];

  // Skip equality rows  
  if  (rowsense=='E' || rowsense=='N') {
    return 0; 
  }
  // le row  
  if (rowsense=='L'){
    b=si.getRightHandSide()[rowIndex];
  }
  // ge row
  // Multiply through by -1 to convert it to a le row 
  if (rowsense=='G'){
    b=-si.getRightHandSide()[rowIndex];
    sign=-1.0;
  }
  
  // Finite, but unequal row bounds  
  // Could derive an simple rounding inequality from either 
  // (or from both!) but for expediency, 
  // use the le relationship as the default for now  
  if  (rowsense=='R') {
    b=si.getRightHandSide()[rowIndex];
  }
  
   // Try to net out the continuous variables from the constraint.
  // Multipy through by sign to convert the inequality to a le inequality  
  // If the coefficient on a continuous variable is positive, replace
  // the continous variable with its lower bound 
  // If the coefficient on a continuous variable is negative, replace
  // the continuous variable with its upper bound.
  // example:
  //                   2.5 <= 3x0-2.8x1+4x2,  0<=x0<=0.2, 0.4<=x1, x2 integer
  //                         -3x0+2.8x1-4x2 <= -2.5
  // -3(0.2)+2.8(0.4)-4x3 <= -3x0+2.8x1-4x2 <= -2.5
  // gives the (weaker, valid) integer inequality
  //                                   -4x2 <= -2.5+3(0.2)-2.8(0.4)
  // sign = -1
  // irow.elements = 4
  // irow.indices = 2
  // negative = (true, true, false)
  // b=-2.5+3(0.2)-2.8(0.4)= -3.02

  const double * colupper = si.getColUpper();
  const double * collower = si.getColLower();
  const char * intVar = si.getColType();
  for (i=0; i<sizeOfRow; i++){
    // if the variable is continuous
    if ( !intVar[matrixRow.getIndices()[i]] ) {
      // and the coefficient is strictly negative
      if((sign*matrixRow.getElements()[i])<-epsilon_){
        // and the continuous variable has a fintite upper bound
        if (colupper[matrixRow.getIndices()[i]] < si.getInfinity()){
          // then replace the variable with its upper bound.
          b=b-(sign*matrixRow.getElements()[i]*colupper[matrixRow.getIndices()[i]]);
        } 
        else 
          return 0;
      }
      // if the coefficient in strictly positive
      else if((sign*matrixRow.getElements()[i])>epsilon_){
        // and the continuous variable has a finite lower bound
        if (collower[matrixRow.getIndices()[i]] > -si.getInfinity()){
          // then replace the variable with its lower bound.
          b=b-(sign*matrixRow.getElements()[i]*collower[matrixRow.getIndices()[i]]);
        }
        else
          return 0;
      }
      // else the coefficient is essentially an explicitly stored zero; do
      // nothing   
    }
    // else: the variable is integer
    else{
      // if the integer variable is fixed, net it out of the integer inequality
      if (colupper[matrixRow.getIndices()[i]]- collower[matrixRow.getIndices()[i]]<
	  epsilon_){
          b=b-(sign*matrixRow.getElements()[i]*colupper[matrixRow.getIndices()[i]]);
      }
      // else the variable is a free integer variable and it becomes
      // part of the integer inequality
      else {
        irow.insert(matrixRow.getIndices()[i],sign*matrixRow.getElements()[i]);
      }
    }
  }
  
  // if there are no free integer variables, then abandon this row;
  if(irow.getNumElements() == 0){
    return 0;
  }
  
  // Store the values and the signs of the coefficients separately.
  // irow.elements stores the absolute values of the coefficients
  // negative indicates the sign.
  // Note: after this point b is essentially the effecitve rhs of a le
  // contraint
  {
     const int s = irow.getNumElements();
     const int * indices = irow.getIndices();
     double * elements = irow.getElements();
     for(i=0; i<s; i++){
	if (elements[i] < -epsilon_) {
	   negative[indices[i]]= true; // store indicator of the sign 
	   elements[i] *= -1;          // store only positive values
	}
    }
  }

  return 1;
}


//-------------------------------------------------------------------
// power10ToMakeDoubleAnInt: 
//   given a vector of positive doubles x_i, i=1, size, and a positive
//   tolerance dataTol, determine the smallest power of 10 needed so that
//   x[i]*10**power is integer for all i.

//   dataTol_ should be correlated to the accuracy of the data,
//   and choosen to be the largest value that's tolerable.
//   
//   (Easily extended to take an input vector of arbitrary sign)
//-------------------------------------------------------------------
//
int
CglSimpleRounding::power10ToMakeDoubleAnInt( 
    int size,             // the length of the input vector x
    const double * x,     // the input vector of postive values  
    double dataTol) const // the (strictly postive) precision of the data

{
  // Assumption: data precision is positive
  assert( dataTol > 0 );


  int i;           // loop iterator 
  int maxPower=0;  // maximum power of 10 used to convert any x[i] to an
                   // integer 
                   // this is the number we are after.
  int power = 0;   // power of 10 used to convert a particular x[i] to an
                   // integer 

#ifdef OLD_MULT
  double intPart;  // the integer part of the number
#endif
  double fracPart; // the fractional part of the number
                   // we keep multiplying by 10 until the fractional part is 0
                   // (well, really just until the factional part is less than
                   // dataTol) 

  // JJF - code seems to fail sometimes as multiplying by 10 - so
  // definition of dataTol changed - see header file

  const double multiplier[16]={1.0,1.0e1,1.0e2,1.0e3,1.0e4,1.0e5,
			       1.0e6,1.0e7,1.0e8,1.0e9,1.0e10,1.0e11,
			       1.0e12,1.0e13,1.0e14,1.0e15};

  // Loop through every element in the array in x
  for (i=0; i<size; i++){
    power = 0;

#ifdef OLD_MULT 
    // look at the fractional part of x[i]
    // FYI: if you want to modify this member function to take an input
    // vector x of arbitary sign, change this line below to 
    // fracPart = modf(fabs(x[i]),&intPart);
    fracPart = modf(x[i],&intPart);

    // if the fractional part is close enough to 0 or 1, we're done with this
    // value
    while(!(fracPart < dataTol || 1-fracPart < dataTol )) {
       // otherwise, multiply by 10 and look at the fractional part of the
       // result. 
       ++power;
       fracPart = fracPart*10.0;
       fracPart = modf(fracPart,&intPart);     
    }
#else
    // use fabs as safer and does no harm
    double value = fabs(x[i]);
    double scaledValue;
    // Do loop - always using original value to stop round off error.
    // If we don't find in 15 goes give up
    for (power=0;power<16;power++) {
      double tolerance = dataTol*multiplier[power];
      scaledValue = value*multiplier[power];
      fracPart = scaledValue-floor(scaledValue);
      if(fracPart < tolerance || 1.0-fracPart < tolerance ) {
	break;
      }
    }
    if (power==16||scaledValue>2147483647) {
#ifdef CGL_DEBUG
      printf("Overflow %g => %g, power %d\n",x[i],scaledValue,power);
#endif
      return -1;
    }
#endif    
#ifdef CGL_DEBUG
    printf("The smallest power of 10 to make %g  integral = %i\n",x[i],power);
#endif

    
    // keep track of the largest power needed so that at the end of the for
    // loop
    // x[i]*10**maxPower will be integral for all i
    if (maxPower < power) maxPower=power;
  }

  return maxPower;
}

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglSimpleRounding::CglSimpleRounding ()
:
CglCutGenerator(),
epsilon_(1.0e-08)
{
  // nothing to do here
}
//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglSimpleRounding::CglSimpleRounding (
                  const CglSimpleRounding & source)
:
CglCutGenerator(source),
epsilon_(source.epsilon_)
{  
  // Nothing to do here
}


//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglSimpleRounding::clone() const
{
  return new CglSimpleRounding(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglSimpleRounding::~CglSimpleRounding ()
{
  // Nothing to do here
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglSimpleRounding &
CglSimpleRounding::operator=(
                   const CglSimpleRounding& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    epsilon_=rhs.epsilon_;
  }
  return *this;
}
// Create C++ lines to get to current state
std::string
CglSimpleRounding::generateCpp( FILE * fp) 
{
  CglSimpleRounding other;
  fprintf(fp,"0#include \"CglSimpleRounding.hpp\"\n");
  fprintf(fp,"3  CglSimpleRounding simpleRounding;\n");
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  simpleRounding.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  simpleRounding.setAggressiveness(%d);\n",getAggressiveness());
  return "simpleRounding";
}
