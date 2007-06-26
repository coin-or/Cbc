// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcCutGenerator_H
#define CbcCutGenerator_H

#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "CglCutGenerator.hpp"

class CbcModel;
class OsiRowCut;
class OsiRowCutDebugger;

//#############################################################################

/** Interface between Cbc and Cut Generation Library.

  \c CbcCutGenerator is intended to provide an intelligent interface between
  Cbc and the cutting plane algorithms in the CGL. A \c CbcCutGenerator is
  bound to a \c CglCutGenerator and to an \c CbcModel. It contains parameters
  which control when and how the \c generateCuts method of the
  \c CglCutGenerator will be called.

  The builtin decision criteria available to use when deciding whether to
  generate cuts are limited: every <i>X</i> nodes, when a solution is found,
  and when a subproblem is found to be infeasible. The idea is that the class
  will grow more intelligent with time.

  \todo Add a pointer to function member which will allow a client to install
	their own decision algorithm to decide whether or not to call the CGL
	\p generateCuts method. Create a default decision method that looks
	at the builtin criteria.

  \todo It strikes me as not good that generateCuts contains code specific to
	individual CGL algorithms. Another set of pointer to function members,
	so that the client can specify the cut generation method as well as
	pre- and post-generation methods? Taken a bit further, should this
	class contain a bunch of pointer to function members, one for each
	of the places where the cut generator might be referenced?
	Initialization, root node, search tree node, discovery of solution,
	and termination all come to mind. Initialization and termination would
	also be useful for instrumenting cbc.
*/

class CbcCutGenerator  {
  
public:
    
  /** \name Generate Cuts */
  //@{
  /** Generate cuts for the client model.

    Evaluate the state of the client model and decide whether to generate cuts.
    The generated cuts are inserted into and returned in the collection of cuts
    \p cs.

    If \p fullScan is true, the generator is obliged to call the CGL
    \c generateCuts routine.  Otherwise, it is free to make a local decision.
    The current implementation uses \c whenCutGenerator_ to decide.

    The routine returns true if reoptimisation is needed (because the state of
    the solver interface has been modified).

    If node then can find out depth
  */
  bool generateCuts( OsiCuts &cs, bool fullScan,CbcNode * node); 
  //@}

    
  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CbcCutGenerator (); 

  /// Normal constructor
  CbcCutGenerator(CbcModel * model,CglCutGenerator * generator,
		  int howOften=1, const char * name=NULL,
		  bool normal=true, bool atSolution=false, 
		  bool infeasible=false,int howOftenInsub=-100,
		  int whatDepth=-1, int whatDepthInSub=-1,int switchOffIfLessThan=0);
 
  /// Copy constructor 
  CbcCutGenerator (const CbcCutGenerator &);

  /// Assignment operator 
  CbcCutGenerator & operator=(const CbcCutGenerator& rhs);

  /// Destructor 
  ~CbcCutGenerator ();
  //@}

  /**@name Gets and sets */
  //@{
  /** Set the client model.
  
    In addition to setting the client model, refreshModel also calls
    the \c refreshSolver method of the CglCutGenerator object.
  */
  void refreshModel(CbcModel * model);

  /// return name of generator
  inline const char * cutGeneratorName() const
  { return generatorName_;};

  /** Set the cut generation interval

    Set the number of nodes evaluated between calls to the Cgl object's
    \p generateCuts routine.

    If \p value is positive, cuts will always be generated at the specified
    interval.
    If \p value is negative, cuts will initially be generated at the specified
    interval, but Cbc may adjust the value depending on the success of cuts
    produced by this generator.

    A value of -100 disables the generator, while a value of -99 means
    just at root.
  */
  void setHowOften(int value) ;

  /// Get the cut generation interval.
  inline int howOften() const
  { return whenCutGenerator_;};
  /// Get the cut generation interval.in sub tree
  inline int howOftenInSub() const
  { return whenCutGeneratorInSub_;};

  /** Set the cut generation depth

    Set the depth criterion for calls to the Cgl object's
    \p generateCuts routine.  Only active if > 0.

    If whenCutGenerator is positive and this is positive then this overrides.  
    If whenCutGenerator is -1 then this is used as criterion if any cuts
    were generated at root node.
    If whenCutGenerator is anything else this is ignored.
  */
  void setWhatDepth(int value) ;
  /// Set the cut generation depth in sub tree
  void setWhatDepthInSub(int value) ;
  /// Get the cut generation depth criterion.
  inline int whatDepth() const
  { return depthCutGenerator_;};
  /// Get the cut generation depth criterion.in sub tree
  inline int whatDepthInSub() const
  { return depthCutGeneratorInSub_;};

  /// Get whether the cut generator should be called in the normal place
  inline bool normal() const
  { return normal_;};
  /// Set whether the cut generator should be called in the normal place
  inline void setNormal(bool value) 
  { normal_=value;};
  /// Get whether the cut generator should be called when a solution is found
  inline bool atSolution() const
  { return atSolution_;};
  /// Set whether the cut generator should be called when a solution is found
  inline void setAtSolution(bool value) 
  { atSolution_=value;};
  /** Get whether the cut generator should be called when the subproblem is
      found to be infeasible.
  */
  inline bool whenInfeasible() const
  { return whenInfeasible_;};
  /** Set whether the cut generator should be called when the subproblem is
      found to be infeasible.
  */
  inline void setWhenInfeasible(bool value) 
  { whenInfeasible_=value;};
  /// Get whether the cut generator is being timed
  inline bool timing() const
  { return timing_;};
  /// Set whether the cut generator is being timed
  inline void setTiming(bool value) 
  { timing_=value; timeInCutGenerator_=0.0;};
  /// Return time taken in cut generator
  inline double timeInCutGenerator() const
  { return timeInCutGenerator_;};
  /// Get the \c CglCutGenerator corresponding to this \c CbcCutGenerator.
  inline CglCutGenerator * generator() const
  { return generator_;};
  /// Number times cut generator entered
  inline int numberTimesEntered() const
  { return numberTimes_;};
  inline void setNumberTimesEntered(int value)
  { numberTimes_ = value;};
  inline void incrementNumberTimesEntered(int value=1)
  { numberTimes_ += value;};
  /// Total number of cuts added
  inline int numberCutsInTotal() const
  { return numberCuts_;};
  inline void setNumberCutsInTotal(int value)
  { numberCuts_ = value;};
  inline void incrementNumberCutsInTotal(int value=1)
  { numberCuts_ += value;};
  /// Total number of column cuts
  inline int numberColumnCuts() const
  { return numberColumnCuts_;};
  inline void setNumberColumnCuts(int value)
  { numberColumnCuts_ = value;};
  inline void incrementNumberColumnCuts(int value=1)
  { numberColumnCuts_ += value;};
  /// Total number of cuts active after (at end of n cut passes at each node)
  inline int numberCutsActive() const
  { return numberCutsActive_;};
  inline void setNumberCutsActive(int value)
  { numberCutsActive_ = value;};
  inline void incrementNumberCutsActive(int value=1)
  { numberCutsActive_ += value;};
  inline void setSwitchOffIfLessThan(int value) 
  { switchOffIfLessThan_ = value;};
  inline int switchOffIfLessThan() const
  { return switchOffIfLessThan_;};
  /// Say if optimal basis needed
  inline bool needsOptimalBasis() const
  { return generator_->needsOptimalBasis();};
  /// Whether generator MUST be called again if any cuts (i.e. ignore break from loop)
  inline bool mustCallAgain() const
  { return mustCallAgain_;};
  /// Set whether generator MUST be called again if any cuts (i.e. ignore break from loop)
  inline void setMustCallAgain(bool yesNo)
  { mustCallAgain_=yesNo;};
  /// Whether generator switched off for moment
  inline bool switchedOff() const
  { return switchedOff_;};
  /// Set whether generator switched off for moment
  inline void setSwitchedOff(bool yesNo)
  { switchedOff_=yesNo;};
  /// Number of cuts generated at root
  inline int numberCutsAtRoot() const
  { return numberCutsAtRoot_;};
  inline void setNumberCutsAtRoot(int value)
  { numberCutsAtRoot_ = value;};
  /// Number of cuts active at root
  inline int numberActiveCutsAtRoot() const
  { return numberActiveCutsAtRoot_;};
  inline void setNumberActiveCutsAtRoot(int value)
  { numberActiveCutsAtRoot_ = value;};
  //@}
  
private:
  /// The client model
  CbcModel *model_;

  // The CglCutGenerator object
  CglCutGenerator * generator_;

  /** Number of nodes between calls to the CglCutGenerator::generateCuts
     routine.
  */
  int whenCutGenerator_;
  /** Number of nodes between calls to the CglCutGenerator::generateCuts
     routine in sub tree.
  */
  int whenCutGeneratorInSub_;
  /** If first pass at root produces fewer than this cuts then switch off
   */
  int switchOffIfLessThan_;

  /** Depth at which to call the CglCutGenerator::generateCuts
     routine (If >0 then overrides when and is called if depth%depthCutGenerator==0).
  */
  int depthCutGenerator_;

  /** Depth at which to call the CglCutGenerator::generateCuts
     routine (If >0 then overrides when and is called if depth%depthCutGenerator==0).
     In sub tree.
  */
  int depthCutGeneratorInSub_;

  /// Name of generator
  char * generatorName_;

  /// Whether to call the generator in the normal place
  bool normal_;

  /// Whether to call the generator when a new solution is found
  bool atSolution_;

  /// Whether to call generator when a subproblem is found to be infeasible
  bool whenInfeasible_;
  /// Whether generator MUST be called again if any cuts (i.e. ignore break from loop)
  bool mustCallAgain_;
  /// Temporary switch off marker
  bool switchedOff_;
  /// Whether call generator being timed
  bool timing_;
  /// Time in cut generator
  double timeInCutGenerator_;
  
  /// Number times cut generator entered
  int numberTimes_;
  /// Total number of cuts added
  int numberCuts_;
  /// Total number of column cuts added
  int numberColumnCuts_;
  /// Total number of cuts active after (at end of n cut passes at each node)
  int numberCutsActive_;
  /// Number of cuts generated at root
  int numberCutsAtRoot_;
  /// Number of cuts active at root
  int numberActiveCutsAtRoot_;
};
/** Abstract cut modifier base class

    In exotic circumstances - cuts may need to be modified
    a) strengthened - changed
    b) weakened - changed
    c) deleted - set to NULL
    d) unchanged
*/

class CbcCutModifier {
public:
  /// Default Constructor 
  CbcCutModifier ();

  // Copy constructor 
  CbcCutModifier ( const CbcCutModifier &);
   
  /// Destructor
  virtual ~CbcCutModifier();

  /// Assignment 
  CbcCutModifier & operator=(const CbcCutModifier& rhs);
 /// Clone
  virtual CbcCutModifier * clone() const = 0;

  /** Returns
      0 unchanged
      1 strengthened
      2 weakened
      3 deleted
  */
  virtual int modify(const OsiSolverInterface * solver, OsiRowCut & cut) =0;
  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) {};
protected:
  
};

/** Simple cut modifier base class

    In exotic circumstances - cuts may need to be modified
    a) strengthened - changed
    b) weakened - changed
    c) deleted - set to NULL
    d) unchanged

    initially get rid of cuts with variables >= k
    could weaken
*/

class CbcCutSubsetModifier  : public CbcCutModifier {
public:
  /// Default Constructor 
  CbcCutSubsetModifier ();

  /// Useful Constructor 
  CbcCutSubsetModifier (int firstOdd);

  // Copy constructor 
  CbcCutSubsetModifier ( const CbcCutSubsetModifier &);
   
  /// Destructor
  virtual ~CbcCutSubsetModifier();

  /// Assignment 
  CbcCutSubsetModifier & operator=(const CbcCutSubsetModifier& rhs);
 /// Clone
  virtual CbcCutModifier * clone() const ;

  /** Returns
      0 unchanged
      1 strengthened
      2 weakened
      3 deleted
  */
  virtual int modify(const OsiSolverInterface * solver, OsiRowCut & cut) ;
  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) {};
protected:
  /// data
  /// First odd variable
  int firstOdd_;
};

#endif
