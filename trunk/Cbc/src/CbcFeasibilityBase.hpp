// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcFeasibilityBase_H
#define CbcFeasibilityBase_H


//#############################################################################
/*  There are cases where the user wants to control how CBC sees the problems feasibility.
    The user may want to examine the problem and say :
    a) The default looks OK
    b) Pretend this problem is Integer feasible 
    c) Pretend this problem is infeasible even though it looks feasible
    
    This simple class allows user to do that.

*/

class CbcModel;
class CbcFeasibilityBase {
public:
  // Default Constructor 
  CbcFeasibilityBase () {}

  /**
     On input mode:
     0 - called after a solve but before any cuts
     -1 - called after strong branching
     Returns :
     0 - no opinion
     -1 pretend infeasible
     1 pretend integer solution
  */
  virtual int feasible(CbcModel * model, int mode) {return 0;}

  virtual ~CbcFeasibilityBase() {}

  // Copy constructor 
  CbcFeasibilityBase ( const CbcFeasibilityBase & rhs)
  {}
   
  // Assignment operator 
  CbcFeasibilityBase & operator=( const CbcFeasibilityBase& rhs)
  {  
    return *this;
  }

  /// Clone
  virtual CbcFeasibilityBase * clone() const
  { return new CbcFeasibilityBase(*this);}
};
#endif
