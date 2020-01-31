// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcSolverLongThin_H
#define CbcSolverLongThin_H

#include "OsiClpSolverInterface.hpp"
class CbcModel;
//#############################################################################

/**

    This is to allow the user to replace initialSolve and resolve
*/

class CbcSolverLongThin : public OsiClpSolverInterface {

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
  CbcSolverLongThin();

  /// Clone
  virtual OsiSolverInterface *clone(bool CopyData = true) const;

  /// Copy constructor
  CbcSolverLongThin(const CbcSolverLongThin &);

  /// Assignment operator
  CbcSolverLongThin &operator=(const CbcSolverLongThin &rhs);

  /// Destructor
  virtual ~CbcSolverLongThin();

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
  /// Say whether to believe infeasible
  inline void setBelieveInfeasible(bool yesNo)
  {
    believeInfeasible_ = yesNo;
  }
  /// Say whether to just count usage
  inline void setAlgorithm(int value)
  {
    algorithm_ = value;
  }
  /// Do nested search if this fraction fixed
  inline void setNested(double value)
  {
    nestedSearch_ = value;
  }
  /// Say whether to just count usage
  inline int getAlgorithm() const
  {
    return algorithm_;
  }
  /// Do nested search if this fraction fixed
  inline double getNested() const
  {
    return nestedSearch_;
  }
  //@}

  //---------------------------------------------------------------------------

private:
  /**@name Private member data */
  //@{
  /// Do nested search if this fraction fixed
  double nestedSearch_;
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
  /// If infeasible on subset means infeasible
  bool believeInfeasible_;
  /// If 0 nothing, 1 compress and fix, 2 long thin
  bool algorithm_;
  //@}
};

#endif
