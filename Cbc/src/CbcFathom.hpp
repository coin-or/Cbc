/* $Id$ */
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcFathom_H
#define CbcFathom_H
#include "CbcConfig.h"

/*
  This file contains two classes, CbcFathom and CbcOsiSolver. It's unclear why
  they're in the same file. CbcOsiSolver is a base class for CbcLinked.

  --lh, 071031 --
*/

class CbcModel;

//#############################################################################
/** Fathom base class.

    The idea is that after some branching the problem will be effectively smaller than
    the original problem and maybe there will be a more specialized technique which can completely
    fathom this branch quickly.

    One method is to presolve the problem to give a much smaller new problem and then do branch
    and cut on that.  Another might be dynamic programming.

 */

class CbcFathom {
public:
  // Default Constructor
  CbcFathom();

  // Constructor with model - assumed before cuts
  CbcFathom(CbcModel &model);

  virtual ~CbcFathom();

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel *model);

  /// Clone
  virtual CbcFathom *clone() const = 0;

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel *model) = 0;

  /** returns 0 if no fathoming attempted, 1 fully fathomed,
        2 incomplete search, 3 incomplete search but treat as complete.
        If solution then newSolution will not be NULL and
        will be freed by CbcModel.  It is expected that the solution is better
        than best so far but CbcModel will double check.

        If returns 3 then of course there is no guarantee of global optimum
    */
  virtual int fathom(double *&newSolution) = 0;

  // Is this method possible
  inline bool possible() const
  {
    return possible_;
  }

protected:
  /// Model
  CbcModel *model_;
  /// Possible - if this method of fathoming can be used
  bool possible_;

private:
  /// Illegal Assignment operator
  CbcFathom &operator=(const CbcFathom &rhs);
};

#include "OsiClpSolverInterface.hpp"

//#############################################################################

/**

This is for codes where solver needs to know about CbcModel
  Seems to provide only one value-added feature, a CbcModel object.

*/

class CbcOsiSolver : public OsiClpSolverInterface {

public:
  /**@name Constructors and destructors */
  //@{
  /// Default Constructor
  CbcOsiSolver();

  /// Clone
  virtual OsiSolverInterface *clone(bool copyData = true) const;

  /// Copy constructor
  CbcOsiSolver(const CbcOsiSolver &);

  /// Assignment operator
  CbcOsiSolver &operator=(const CbcOsiSolver &rhs);

  /// Destructor
  virtual ~CbcOsiSolver();

  //@}

  /**@name Sets and Gets */
  //@{
  /// Set Cbc Model
  inline void setCbcModel(CbcModel *model)
  {
    cbcModel_ = model;
  }
  /// Return Cbc Model
  inline CbcModel *cbcModel() const
  {
    return cbcModel_;
  }
  //@}

  //---------------------------------------------------------------------------

protected:
  /**@name Private member data */
  //@{
  /// Pointer back to CbcModel
  CbcModel *cbcModel_;
  //@}
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
