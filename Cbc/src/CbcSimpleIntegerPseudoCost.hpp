// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#ifndef CbcSimpleIntegerPseudoCost_H
#define CbcSimpleIntegerPseudoCost_H

#include "CbcSimpleInteger.hpp"
/// Define a single integer class but with pseudo costs

class CbcSimpleIntegerPseudoCost : public CbcSimpleInteger {

public:
  // Default Constructor
  CbcSimpleIntegerPseudoCost();

  // Useful constructor - passed model index
  CbcSimpleIntegerPseudoCost(CbcModel *model, int iColumn, double breakEven = 0.5);

  // Useful constructor - passed and model index and pseudo costs
  CbcSimpleIntegerPseudoCost(CbcModel *model, int iColumn,
    double downPseudoCost, double upPseudoCost);
  // Useful constructor - passed and model index and pseudo costs
  CbcSimpleIntegerPseudoCost(CbcModel *model, int dummy, int iColumn,
    double downPseudoCost, double upPseudoCost);

  // Copy constructor
  CbcSimpleIntegerPseudoCost(const CbcSimpleIntegerPseudoCost &);

  /// Clone
  virtual CbcObject *clone() const;

  // Assignment operator
  CbcSimpleIntegerPseudoCost &operator=(const CbcSimpleIntegerPseudoCost &rhs);

  // Destructor
  virtual ~CbcSimpleIntegerPseudoCost();

  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiBranchingInformation *info,
    int &preferredWay) const;

  /// Creates a branching object
  virtual CbcBranchingObject *createCbcBranch(OsiSolverInterface *solver, const OsiBranchingInformation *info, int way);

  /// Down pseudo cost
  inline double downPseudoCost() const
  {
    return downPseudoCost_;
  }
  /// Set down pseudo cost
  inline void setDownPseudoCost(double value)
  {
    downPseudoCost_ = value;
  }

  /// Up pseudo cost
  inline double upPseudoCost() const
  {
    return upPseudoCost_;
  }
  /// Set up pseudo cost
  inline void setUpPseudoCost(double value)
  {
    upPseudoCost_ = value;
  }

  /// Up down separator
  inline double upDownSeparator() const
  {
    return upDownSeparator_;
  }
  /// Set up down separator
  inline void setUpDownSeparator(double value)
  {
    upDownSeparator_ = value;
  }

  /// Return "up" estimate
  virtual double upEstimate() const;
  /// Return "down" estimate (default 1.0e-5)
  virtual double downEstimate() const;

  /// method - see below for details
  inline int method() const
  {
    return method_;
  }
  /// Set method
  inline void setMethod(int value)
  {
    method_ = value;
  }

protected:
  /// data

  /// Down pseudo cost
  double downPseudoCost_;
  /// Up pseudo cost
  double upPseudoCost_;
  /** Up/down separator
        If >0.0 then do first branch up if value-floor(value)
        >= this value
    */
  double upDownSeparator_;
  /** Method -
        0 - normal - return min (up,down)
        1 - if before any solution return CoinMax(up,down)
        2 - if before branched solution return CoinMax(up,down)
        3 - always return CoinMax(up,down)
    */
  int method_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
