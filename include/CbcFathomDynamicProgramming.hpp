// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcFathomDynamicProgramming_H
#define CbcFathomDynamicProgramming_H

#include "CbcFathom.hpp"

//#############################################################################
/** FathomDynamicProgramming base class.

    The idea is that after some branching the problem will be effectively smaller than
    the original problem and maybe there will be a more specialized technique which can completely
    fathom this branch quickly.

    One method is to presolve the problem to give a much smaller new problem and then do branch 
    and cut on that.  Another might be dynamic programming.

 */

class CbcFathomDynamicProgramming : public CbcFathom {
public:
  // Default Constructor 
  CbcFathomDynamicProgramming ();

  // Constructor with model - assumed before cuts
  CbcFathomDynamicProgramming (CbcModel & model);
  // Copy constructor 
  CbcFathomDynamicProgramming(const CbcFathomDynamicProgramming & rhs);

  virtual ~CbcFathomDynamicProgramming();

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel * model);
  
  /// Clone
  virtual CbcFathom * clone() const;

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel * model);

  /** returns 0 if no fathoming attempted, 1 fully fathomed ,
      2 incomplete search, 3 incomplete search but treat as complete.
      If solution then newSolution will not be NULL and
      will be freed by CbcModel.  It is expected that the solution is better
      than best so far but CbcModel will double check.

      If returns 3 then of course there is no guarantee of global optimum
  */
  virtual int fathom(double *& newSolution);

  /// Maximum size allowed
  inline int maximumSize() const
  { return maximumSizeAllowed_;};
  inline void setMaximumSize(int value)
  { maximumSizeAllowed_=value;};
private:
  /// Returns type
  int gutsOfCheckPossible(int allowableSize=0);
  /// Does deleteions
  void gutsOfDelete();

  /** Adds one attempt of one column of type 0,
      returns true if was used in making any changes
  */
  bool addOneColumn0(int id,int numberElements, const int * rows,
                     double cost);
  /** Adds one attempt of one column of type 1,
      returns true if was used in making any changes.
      At present the user has to call it once for each possible value
  */
  bool addOneColumn1(int id,int numberElements, const int * rows,
                     const int * coefficients, double cost);
  /** Adds one attempt of one column of type 1,
      returns true if was used in making any changes.
      At present the user has to call it once for each possible value.
      This version is when there are enough 1 rhs to do faster
  */
  bool addOneColumn1A(int id,int numberElements, const int * rows,
                     const int * coefficients, double cost);

protected:

  /// Size of states (power of 2 unless just one constraint)
  int size_;
  /** Type - 0 coefficients and rhs all 1,
      1 - coefficients > 1 or rhs > 1
  */
  int type_;
  /// Space for states (? could be float)
  double * cost_;
  /// Which state produced this cheapest one
  int * back_;
  /// Some id as to which variable was used
  int * id_;
  /// Some rows may be satisified so we need a lookup
  int * lookup_;
  /// Number of active rows
  int numberActive_;
  /// Maximum size allowed
  int maximumSizeAllowed_;
  /// Start bit for each active row
  int * startBit_;
  /// Number bits for each active row
  int * numberBits_;
  /// Effective rhs
  int * rhs_;
  /// Number of Non 1 rhs
  int numberNonOne_;
private:
  
  /// Illegal Assignment operator 
  CbcFathomDynamicProgramming & operator=(const CbcFathomDynamicProgramming& rhs);
  
};

#endif
