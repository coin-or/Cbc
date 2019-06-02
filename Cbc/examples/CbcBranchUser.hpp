// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcBranchUser_H
#define CbcBranchUser_H

#include "CbcBranchBase.hpp"
#include "CbcBranchActual.hpp"

/** Branching decision user class */

class CbcBranchUserDecision : public CbcBranchDecision {
public:
  // Default Constructor
  CbcBranchUserDecision();

  // Copy constructor
  CbcBranchUserDecision(const CbcBranchUserDecision &);

  virtual ~CbcBranchUserDecision();

  /// Clone
  virtual CbcBranchDecision *clone() const;

  /// Initialize i.e. before start of choosing at a node
  virtual void initialize(CbcModel *model);

  /** Returns nonzero if branching on first object is "better" than on
      second (if second NULL first wins).
      This is only used after strong branching.  The initial selection
      is done by infeasibility() for each CbcObject
      return code +1 for up branch preferred, -1 for down
      
 */
  virtual int betterBranch(CbcBranchingObject *thisOne,
    CbcBranchingObject *bestSoFar,
    double changeUp, int numberInfeasibilitiesUp,
    double changeDown, int numberInfeasibilitiesDown);

  /** \brief Compare N branching objects. Return index of best
      and sets way of branching in chosen object.
    
    This routine is used only after strong branching.
    This is reccommended version as it can be more sophisticated
  */

  virtual int
  bestBranch(CbcBranchingObject **objects, int numberObjects, int numberUnsatisfied,
    double *changeUp, int *numberInfeasibilitiesUp,
    double *changeDown, int *numberInfeasibilitiesDown,
    double objectiveValue);

private:
  /// Illegal Assignment operator
  CbcBranchUserDecision &operator=(const CbcBranchUserDecision &rhs);
};

/// Define a single integer class where branching is forced until fixed

class CbcSimpleIntegerFixed : public CbcSimpleInteger {

public:
  // Default Constructor
  CbcSimpleIntegerFixed();

  // Useful constructor - passed integer index and model index
  CbcSimpleIntegerFixed(CbcModel *model, int iColumn, double breakEven = 0.5);

  // Constructor from simple
  CbcSimpleIntegerFixed(const CbcSimpleInteger &simple);

  // Copy constructor
  CbcSimpleIntegerFixed(const CbcSimpleIntegerFixed &);

  /// Clone
  virtual CbcObject *clone() const;

  // Assignment operator
  CbcSimpleIntegerFixed &operator=(const CbcSimpleIntegerFixed &rhs);

  // Destructor
  ~CbcSimpleIntegerFixed();

  /// Infeasibility - large is 0.5
  virtual double infeasibility(int &preferredWay) const;

  /** Creates a branching object

    The preferred direction is set by \p way, -1 for down, +1 for up.
  */
  //virtual CbcBranchingObject * createBranch(int way) ;
  /** Create a branching object and indicate which way to branch first.
      
      The branching object has to know how to create branches (fix
      variables, etc.)
  */
  virtual CbcBranchingObject *createBranch(OsiSolverInterface *solver,
    const OsiBranchingInformation *info, int way);

protected:
  /// data
};

#endif
