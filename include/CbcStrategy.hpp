// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcStrategy_H
#define CbcStrategy_H

#include "CbcModel.hpp"


//#############################################################################
/** Strategy base class */

class CbcStrategy {
public:
  // Default Constructor 
  CbcStrategy ();

  virtual ~CbcStrategy();

  /// Clone
  virtual CbcStrategy * clone() const=0;

  /// Setup cut generators
  virtual void setupCutGenerators(CbcModel & model);
  /// Setup heuristics
  virtual void setupHeuristics(CbcModel & model);
  /// Do printing stuff
  virtual void setupPrinting(CbcModel & model);
  /// Other stuff e.g. strong branching
  virtual void setupOther(CbcModel & model);
  /// Set model depth (i.e. how nested)
  inline void setNested(int depth)
  { depth_=depth;};
  /// Get model depth (i.e. how nested)
  inline int getNested() const
  { return depth_;};
private:
  
  /// Illegal Assignment operator 
  CbcStrategy & operator=(const CbcStrategy& rhs);
protected:
  // Data
  /// Model depth
  int depth_;
};

/** Default class
 */

class CbcStrategyDefault : public CbcStrategy {
public:

  // Default Constructor 
  CbcStrategyDefault (bool cutsOnlyAtRoot=true,
                      int numberStrong=5,
                      int printLevel=0);

  // Copy constructor 
  CbcStrategyDefault ( const CbcStrategyDefault &);
   
  // Destructor 
  ~CbcStrategyDefault ();
  
  /// Clone
  virtual CbcStrategy * clone() const;

  /// Setup cut generators
  virtual void setupCutGenerators(CbcModel & model);
  /// Setup heuristics
  virtual void setupHeuristics(CbcModel & model);
  /// Do printing stuff
  virtual void setupPrinting(CbcModel & model);
  /// Other stuff e.g. strong branching
  virtual void setupOther(CbcModel & model);

protected:
  // Data

  // Whether to do cuts only at root
  bool cutsOnlyAtRoot_;

  // How much strong branching to do 
  int numberStrong_;

  // Print level 0 little, 1 medium
  int printLevel_;

private:
  /// Illegal Assignment operator 
  CbcStrategyDefault & operator=(const CbcStrategyDefault& rhs);
};


#endif
