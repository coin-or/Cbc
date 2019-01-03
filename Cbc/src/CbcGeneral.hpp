// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#ifndef CbcGeneral_H
#define CbcGeneral_H

#include "CbcBranchBase.hpp"
/** Define a catch all class.
    This will create a list of subproblems
*/

class CbcGeneral : public CbcObject {

public:
  // Default Constructor
  CbcGeneral();

  /** Useful constructor
        Just needs to point to model.
    */
  CbcGeneral(CbcModel *model);

  // Copy constructor
  CbcGeneral(const CbcGeneral &);

  /// Clone
  virtual CbcObject *clone() const = 0;

  // Assignment operator
  CbcGeneral &operator=(const CbcGeneral &rhs);

  // Destructor
  ~CbcGeneral();

  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiBranchingInformation *info,
    int &preferredWay) const;

  using CbcObject::feasibleRegion;
  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion() = 0;

  /// Creates a branching object
  virtual CbcBranchingObject *createCbcBranch(OsiSolverInterface *solver, const OsiBranchingInformation *info, int way);

  /// Redoes data when sequence numbers change
  virtual void redoSequenceEtc(CbcModel *model, int numberColumns, const int *originalColumns) = 0;

protected:
  /// data
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
