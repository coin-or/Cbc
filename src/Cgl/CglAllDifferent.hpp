// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglAllDifferent_H
#define CglAllDifferent_H

#include <string>

#include "CglCutGenerator.hpp"

/** AllDifferent Cut Generator Class 
    This has a number of sets.  All the members in each set are general integer
    variables which have to be different from all others in the set.

    At present this only generates column cuts

    At present it is very primitive compared to proper CSP implementations
 */
class CGLLIB_EXPORT CglAllDifferent : public CglCutGenerator {
 
public:
    
  
  /**@name Generate Cuts */
  //@{
  /** This fixes (or reduces bounds) on sets of all different variables
  */
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info = CglTreeInfo());
  //@}


  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglAllDifferent ();

  /// Useful constructot
  CglAllDifferent(int numberSets, const int * starts, const int * which);
 
  /// Copy constructor 
  CglAllDifferent (
    const CglAllDifferent &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglAllDifferent &
    operator=(
    const CglAllDifferent& rhs);
  
  /// Destructor 
  virtual
    ~CglAllDifferent ();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);

  /// This can be used to refresh any inforamtion
  virtual void refreshSolver(OsiSolverInterface * solver);
  /**
     Returns true if may generate Row cuts in tree (rather than root node).
     Used so know if matrix will change in tree.  Really
     meant so column cut generators can still be active
     without worrying code.
     Default is true
  */
  virtual bool mayGenerateRowCutsInTree() const
  { return false;}
  //@}
  /**@name Sets and Gets */
  //@{
  /// Set log level
  inline void setLogLevel(int value)
  { logLevel_=value;}
  /// Get log level
  inline int getLogLevel() const
  { return logLevel_;}
  /// Set Maximum number of sets to look at at once
  inline void setMaxLook(int value)
  { maxLook_=value;}
  /// Get Maximum number of sets to look at at once
  inline int getMaxLook() const
  { return maxLook_;}
  //@}
      
private:
  
 // Private member methods
  /**@name  */
  //@{
  //@}

  // Private member data

  /**@name Private member data */
  //@{
  /// Number of sets
  int numberSets_;
  /// Total number of variables in all different sets
  int numberDifferent_;
  /// Maximum number of sets to look at at once
  int maxLook_;
  /// Log level - 0 none, 1 - a bit, 2 - more details
  int logLevel_;
  /// Start of each set
  int * start_;
  /// Members (0,1,....) not as in original model
  int * which_;
  /// Original members
  int * originalWhich_;
  //@}
};
#endif
