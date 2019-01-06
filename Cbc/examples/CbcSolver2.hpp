// $Id$
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcSolver2_H
#define CbcSolver2_H

#include "OsiClpSolverInterface.hpp"
class CbcModel;
//#############################################################################

/**

    This is to allow the user to replace initialSolve and resolve.
    This version is to try and speed up long thin problems.

    This particular version assumes unit elements and rhs
    Can be E or G rhs
*/

class CbcSolver2 : public OsiClpSolverInterface {

public:
  //---------------------------------------------------------------------------
  /**@name Solve methods */
  //@{
  /// Solve initial LP relaxation
  virtual void initialSolve();

  /// Resolve an LP relaxation after problem modification
  virtual void resolve();

  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default Constructor
  CbcSolver2();

  /// Clone
  virtual OsiSolverInterface *clone(bool CopyData = true) const;

  /// Copy constructor
  CbcSolver2(const CbcSolver2 &);

  /// Assignment operator
  CbcSolver2 &operator=(const CbcSolver2 &rhs);

  /// Destructor
  virtual ~CbcSolver2();

  //@}

  /**@name Sets and Getss */
  //@{
  /// Setup arrays - ones in keep will always be in
  void initialize(CbcModel *model, const char *keep);
  /// get which ones have been used
  inline const int *when() const
  {
    return node_;
  }
  /// Get memory (i.e. how recent use should be)
  inline int getMemory() const
  {
    return memory_;
  }
  /// Get current count
  inline int getCount() const
  {
    return count_;
  }
  /// Set memory (i.e. how recent use should be)
  inline void setMemory(int value)
  {
    memory_ = value;
  }
  /// Say whether to just count usage
  inline void setAlgorithm(int value)
  {
    algorithm_ = value;
  }
  /// Say whether to just count usage
  inline int getAlgorithm() const
  {
    return algorithm_;
  }
  /// Strategy
  inline void setStrategy(int value)
  {
    strategy_ = value;
  }
  /// Strategy
  inline int getStrategy() const
  {
    return strategy_;
  }
  //@}

  //---------------------------------------------------------------------------

private:
  /**@name Private member data */
  //@{
  /// Node number when variable last in problem
  int *node_;
  /// How many times in problem
  int *howMany_;
  /// Pointer back to model
  CbcModel *model_;
  /// Counter
  int count_;
  /// How recently it must have been used
  int memory_;
  /// If 0 nothing, 1 compress and fix, 2 long thin
  int algorithm_;
  /// If 0 get rid of rows, 1 keep rows (to stay dual feasible)
  int strategy_;
  //@}
};

#endif
