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
  virtual void setupCutGenerators(CbcModel & model)=0;
  /// Setup heuristics
  virtual void setupHeuristics(CbcModel & model)=0;
  /// Do printing stuff
  virtual void setupPrinting(CbcModel & model,int modelLogLevel)=0;
  /// Other stuff e.g. strong branching
  virtual void setupOther(CbcModel & model)=0;
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

/** Null class
 */

class CbcStrategyNull : public CbcStrategy {
public:

  // Default Constructor 
  CbcStrategyNull () {}

  // Copy constructor 
  CbcStrategyNull ( const CbcStrategyNull & rhs) : CbcStrategy(rhs) {};
   
  // Destructor 
  ~CbcStrategyNull () {};
  
  /// Clone
  virtual CbcStrategy * clone() const { return new CbcStrategyNull(*this);};

  /// Setup cut generators
  virtual void setupCutGenerators(CbcModel & model) {};
  /// Setup heuristics
  virtual void setupHeuristics(CbcModel & model) {};
  /// Do printing stuff
  virtual void setupPrinting(CbcModel & model,int modelLogLevel) {};
  /// Other stuff e.g. strong branching
  virtual void setupOther(CbcModel & model) {};

protected:
  // Data
private:
  /// Illegal Assignment operator 
  CbcStrategyNull & operator=(const CbcStrategyNull& rhs);
};

/** Default class
 */

class CbcStrategyDefault : public CbcStrategy {
public:

  // Default Constructor 
  CbcStrategyDefault (bool cutsOnlyAtRoot=true,
                      int numberStrong=5,
                      int numberBeforeTrust=0,
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
  virtual void setupPrinting(CbcModel & model,int modelLogLevel) ;
  /// Other stuff e.g. strong branching
  virtual void setupOther(CbcModel & model);

protected:
  // Data

  // Whether to do cuts only at root
  bool cutsOnlyAtRoot_;

  // How much strong branching to do 
  int numberStrong_;

  // Number branches needed to trust with dynamic pseudo costs
  int numberBeforeTrust_;

  // Print level 0 little, 1 medium
  int printLevel_;

private:
  /// Illegal Assignment operator 
  CbcStrategyDefault & operator=(const CbcStrategyDefault& rhs);
};


#endif
