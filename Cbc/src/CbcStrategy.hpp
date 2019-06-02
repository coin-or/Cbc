/* $Id$ */
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcStrategy_H
#define CbcStrategy_H

#include "CbcModel.hpp"
class CglPreProcess;
class CbcNodeInfo;
class CbcNode;
class CoinWarmStartDiff;

//#############################################################################
/** Strategy base class */

class CbcStrategy {
public:
  // Default Constructor
  CbcStrategy();

  virtual ~CbcStrategy();

  /// Clone
  virtual CbcStrategy *clone() const = 0;

  /// Setup cut generators
  virtual void setupCutGenerators(CbcModel &model) = 0;
  /// Setup heuristics
  virtual void setupHeuristics(CbcModel &model) = 0;
  /// Do printing stuff
  virtual void setupPrinting(CbcModel &model, int modelLogLevel) = 0;
  /// Other stuff e.g. strong branching and preprocessing
  virtual void setupOther(CbcModel &model) = 0;
  /// Set model depth (i.e. how nested)
  inline void setNested(int depth)
  {
    depth_ = depth;
  }
  /// Get model depth (i.e. how nested)
  inline int getNested() const
  {
    return depth_;
  }
  /// Say preProcessing done
  inline void setPreProcessState(int state)
  {
    preProcessState_ = state;
  }
  /// See what sort of preprocessing was done
  inline int preProcessState() const
  {
    return preProcessState_;
  }
  /// Pre-processing object
  inline CglPreProcess *process() const
  {
    return process_;
  }
  /// Delete pre-processing object to save memory
  void deletePreProcess();
  /// Return a new Full node information pointer (descendant of CbcFullNodeInfo)
  virtual CbcNodeInfo *fullNodeInfo(CbcModel *model, int numberRowsAtContinuous) const;
  /// Return a new Partial node information pointer (descendant of CbcPartialNodeInfo)
  virtual CbcNodeInfo *partialNodeInfo(CbcModel *model, CbcNodeInfo *parent, CbcNode *owner,
    int numberChangedBounds, const int *variables,
    const double *boundChanges,
    const CoinWarmStartDiff *basisDiff) const;
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *) {}
  /** After a CbcModel::resolve this can return a status
        -1 no effect
        0 treat as optimal
        1 as 0 but do not do any more resolves (i.e. no more cuts)
        2 treat as infeasible
    */
  virtual int status(CbcModel *model, CbcNodeInfo *parent, int whereFrom);

private:
  /// Illegal Assignment operator
  CbcStrategy &operator=(const CbcStrategy &rhs);

protected:
  // Data
  /// Model depth
  int depth_;
  /** PreProcessing state -
        -1 infeasible
        0 off
        1 was done (so need post-processing)
    */
  int preProcessState_;
  /// If preprocessing then this is object
  CglPreProcess *process_;
};

/** Null class
 */

class CbcStrategyNull : public CbcStrategy {
public:
  // Default Constructor
  CbcStrategyNull() {}

  // Copy constructor
  CbcStrategyNull(const CbcStrategyNull &rhs)
    : CbcStrategy(rhs)
  {
  }

  // Destructor
  ~CbcStrategyNull() {}

  /// Clone
  virtual CbcStrategy *clone() const
  {
    return new CbcStrategyNull(*this);
  }

  /// Setup cut generators
  virtual void setupCutGenerators(CbcModel &) {}
  /// Setup heuristics
  virtual void setupHeuristics(CbcModel &) {}
  /// Do printing stuff
  virtual void setupPrinting(CbcModel &, int) {}
  /// Other stuff e.g. strong branching
  virtual void setupOther(CbcModel &) {}

protected:
  // Data
private:
  /// Illegal Assignment operator
  CbcStrategyNull &operator=(const CbcStrategyNull &rhs);
};

/** Default class
 */

class CbcStrategyDefault : public CbcStrategy {
public:
  // Default Constructor
  CbcStrategyDefault(int cutsOnlyAtRoot = 1,
    int numberStrong = 5,
    int numberBeforeTrust = 0,
    int printLevel = 0);

  // Copy constructor
  CbcStrategyDefault(const CbcStrategyDefault &);

  // Destructor
  ~CbcStrategyDefault();

  /// Clone
  virtual CbcStrategy *clone() const;

  /// Setup cut generators
  virtual void setupCutGenerators(CbcModel &model);
  /// Setup heuristics
  virtual void setupHeuristics(CbcModel &model);
  /// Do printing stuff
  virtual void setupPrinting(CbcModel &model, int modelLogLevel);
  /// Other stuff e.g. strong branching
  virtual void setupOther(CbcModel &model);
  /// Set up preProcessing - see below
  inline void setupPreProcessing(int desired = 1, int passes = 10)
  {
    desiredPreProcess_ = desired;
    preProcessPasses_ = passes;
  }
  /// See what sort of preprocessing wanted
  inline int desiredPreProcess() const
  {
    return desiredPreProcess_;
  }
  /// See how many passes wanted
  inline int preProcessPasses() const
  {
    return preProcessPasses_;
  }
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp);

protected:
  // Data

  // Whether to do cuts only at root (-1 -> switch off totally)
  int cutsOnlyAtRoot_;

  // How much strong branching to do
  int numberStrong_;

  // Number branches needed to trust with dynamic pseudo costs
  int numberBeforeTrust_;

  // Print level 0 little, 1 medium
  int printLevel_;

  /** Desired pre-processing
        0 - none
        1 - ordinary
        2 - find sos
        3 - find cliques
        4 - more aggressive sos
        5 - add integer slacks
    */
  int desiredPreProcess_;
  /// Number of pre-processing passes
  int preProcessPasses_;

private:
  /// Illegal Assignment operator
  CbcStrategyDefault &operator=(const CbcStrategyDefault &rhs);
};

/** Default class for sub trees
 */

class CbcStrategyDefaultSubTree : public CbcStrategy {
public:
  // Default Constructor
  CbcStrategyDefaultSubTree(CbcModel *parent = NULL, int cutsOnlyAtRoot = 1,
    int numberStrong = 5,
    int numberBeforeTrust = 0,
    int printLevel = 0);

  // Copy constructor
  CbcStrategyDefaultSubTree(const CbcStrategyDefaultSubTree &);

  // Destructor
  ~CbcStrategyDefaultSubTree();

  /// Clone
  virtual CbcStrategy *clone() const;

  /// Setup cut generators
  virtual void setupCutGenerators(CbcModel &model);
  /// Setup heuristics
  virtual void setupHeuristics(CbcModel &model);
  /// Do printing stuff
  virtual void setupPrinting(CbcModel &model, int modelLogLevel);
  /// Other stuff e.g. strong branching
  virtual void setupOther(CbcModel &model);

protected:
  // Data
  // Parent model
  CbcModel *parentModel_;
  // Whether to do cuts only at root (-1 -> switch off totally)
  int cutsOnlyAtRoot_;

  // How much strong branching to do
  int numberStrong_;

  // Number branches needed to trust with dynamic pseudo costs
  int numberBeforeTrust_;

  // Print level 0 little, 1 medium
  int printLevel_;

private:
  /// Illegal Assignment operator
  CbcStrategyDefaultSubTree &operator=(const CbcStrategyDefaultSubTree &rhs);
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
