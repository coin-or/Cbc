// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglSimpleRounding_H
#define CglSimpleRounding_H

#include <string>

#include "CglCutGenerator.hpp"
#include "CoinPackedMatrix.hpp"

/** Simple Rounding Cut Generator Class

 This class generates simple rounding cuts via the following method:
    For each contraint,
      attempt to derive a <= inequality in all integer variables
      by netting out any continuous variables.
      Divide the resulting integer inequality through by 
      the greatest common denomimator (gcd) of the lhs coefficients.
      Round down the rhs.

 Warning: Use with careful attention to data precision.

 (Reference: Nemhauser and Wolsey, Integer and Combinatorial Optimization, 1988, pg 211.)
*/

class CGLLIB_EXPORT CglSimpleRounding : public CglCutGenerator {
   friend CGLLIB_EXPORT void CglSimpleRoundingUnitTest(const OsiSolverInterface * siP,
					 const std::string mpdDir );
 
public:

  /**@name Generate Cuts */
  //@{
  /** Generate simple rounding cuts for the model accessed through the solver interface. 
  Insert generated cuts into the cut set cs.
  */
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info = CglTreeInfo());
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglSimpleRounding ();
 
  /// Copy constructor 
  CglSimpleRounding (
    const CglSimpleRounding &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglSimpleRounding &
    operator=(
    const CglSimpleRounding& rhs);
  
  /// Destructor 
  virtual
    ~CglSimpleRounding ();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  //@}

private:
  
  // Private member methods
   
  /**@name Private methods */
  //@{
  
  /// Derive a <= inequality in integer variables from the rowIndex-th constraint
  bool deriveAnIntegerRow(
                          const OsiSolverInterface & si,
                          int rowIndex,
                          const CoinShallowPackedVector & matrixRow, 
                          CoinPackedVector & irow,
                          double & b,
                          bool * negative) const;
  

  /** Given a vector of doubles, x, with size elements and a positive tolerance,
     dataTol, this method returns the smallest power of 10 needed so that
     x[i]*10**power "is integer" for all i=0,...,size-1.
  
     ** change of definition of dataTol so that it refers to original
     data, not to scaled data as that seems to lead to problems.

     So if xScaled is x[i]*10**power and xInt is rounded(xScaled)
     then fabs(xScaled-xInt) <= dataTol*10**power.  This means that
     dataTol should be smaller - say 1.0e-12 rather tahn 1.0e-8

     Returns -number of times overflowed  if the power is so big that it will
     cause overflow (i.e. integer stored will be bigger than 2**31).
     Test in cut generator.
  */ 
  int power10ToMakeDoubleAnInt( 
       int size,               // the length of the vector x
       const double * x,   
       double dataTol ) const; // the precision of the data, i.e. the positive
                               // epsilon, which is equivalent to zero

  /**@name Greatest common denominators methods */
  //@{
  /// Returns the greatest common denominator of two positive integers, a and b.
  inline  int gcd(int a, int b) const; 
  
  /** Returns the greatest common denominator of a vector of
      positive integers, vi, of length n.
  */
  inline  int gcdv(int n, const int * const vi) const; 
  //@}

  //@}
  
  /**@name Private member data */
  //@{
  /// A value within an epsilon_ neighborhood of 0  is considered to be 0.
  double epsilon_;
  //@}
};


//-------------------------------------------------------------------
// Returns the greatest common denominator of two 
// positive integers, a and b, found using Euclid's algorithm 
//-------------------------------------------------------------------
int 
CglSimpleRounding::gcd(int a, int b) const
{
  if(a > b) {
    // Swap a and b
    int temp = a;
    a = b;
    b = temp;
  }
  int remainder = b % a;
  if (remainder == 0) return a;
  else return gcd(remainder,a);
}

//-------------------------------------------------------------------
// Returns the greatest common denominator of a vector of
// positive integers, vi, of length n.
//-------------------------------------------------------------------
int 
CglSimpleRounding::gcdv(int n, const int* const vi) const
{
  if (n==0)
    abort();

  if (n==1)
    return vi[0];

  int retval=gcd(vi[0], vi[1]);
  for (int i=2; i<n; i++){
     retval=gcd(retval,vi[i]);
  }
  return retval;
}

//#############################################################################
/** A function that tests the methods in the CglSimpleRounding class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
CGLLIB_EXPORT
void CglSimpleRoundingUnitTest(const OsiSolverInterface * siP,
			       const std::string mpdDir );
  
#endif
