// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglDuplicateRow_H
#define CglDuplicateRow_H

#include <string>

#include "CglCutGenerator.hpp"
class CglStored;

/** DuplicateRow Cut Generator Class */
class CGLLIB_EXPORT CglDuplicateRow : public CglCutGenerator {
 
public:
    
  
  /**@name Generate Cuts */
  //@{
  /** Fix variables and find duplicate/dominated rows for the model of the 
      solver interface, si.

      This is a very simple minded idea but I (JJF) am using it in a project so thought
      I might as well add it.  It should really be called before first solve and I may
      modify CBC to allow for that.

      This is designed for problems with few rows and many integer variables where the rhs
      are <= or == and all coefficients and rhs are small integers.

      If effective rhs is K then we can fix all variables with coefficients > K to their lower bounds
      (effective rhs just means original with variables with nonzero lower bounds subtracted out).

      If one row is a subset of another and the effective rhs are same we can fix some variables
      and then the two rows are identical.

      The generator marks identical rows so can be taken out in solve
  */
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info = CglTreeInfo());
private:
  /// Does work for modes 1,2
  void generateCuts12( const OsiSolverInterface & si, OsiCuts & cs,
		       const CglTreeInfo info = CglTreeInfo());
  /// Does work for mode 4
  void generateCuts4( const OsiSolverInterface & si, OsiCuts & cs,
		       const CglTreeInfo info = CglTreeInfo());
  /// Does work for mode 8
  void generateCuts8( const OsiSolverInterface & si, OsiCuts & cs,
		       const CglTreeInfo info = CglTreeInfo());
public:
  /** Fix variables and find duplicate/dominated rows for the model of the 
      solver interface, si.

      This is a very simple minded idea but I (JJF) am using it in a project so thought
      I might as well add it.  It should really be called before first solve and I may
      modify CBC to allow for that.

      This is designed for problems with few rows and many integer variables where the rhs
      are <= or == and all coefficients and rhs are small integers.

      If effective rhs is K then we can fix all variables with coefficients > K to their lower bounds
      (effective rhs just means original with variables with nonzero lower bounds subtracted out).

      If one row is a subset of another and the effective rhs are same we can fix some variables
      and then the two rows are identical.

      This version does deletions and fixings and may return stored cuts for
      dominated columns 
  */
  CglStored * outDuplicates( OsiSolverInterface * solver);

  //@}

  /**@name Get information on size of problem */
  //@{
  /// Get duplicate row list, -1 means still in, -2 means out (all fixed), k>= means same as row k 
  inline const int * duplicate() const
  { return duplicate_;}
  /// Size of dynamic program
  inline int sizeDynamic() const
  { return sizeDynamic_;}
  /// Number of rows in original problem
  inline int numberOriginalRows() const
  { return matrix_.getNumRows();}
  //@}

  /**@name Get information on size of problem */
  //@{
  /// logLevel
  inline int logLevel() const
  { return logLevel_;}
  inline void setLogLevel(int value)
  { logLevel_ = value;}
  //@}


  /**@name We only check for duplicates amongst rows with effective rhs <= this */
  //@{
  /// Get
  inline int maximumRhs() const
  { return maximumRhs_;}
  /// Set
  inline void setMaximumRhs(int value)
  { maximumRhs_=value;}
  //@}

  /**@name We only check for dominated amongst groups of columns whose size <= this */
  //@{
  /// Get
  inline int maximumDominated() const
  { return maximumDominated_;}
  /// Set
  inline void setMaximumDominated(int value)
  { maximumDominated_=value;}
  //@}
  /**@name gets and sets */
  //@{
  /// Get mode
  inline int mode() const
  { return mode_;}
  /// Set mode
  inline void setMode(int value)
  { mode_=value;}
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglDuplicateRow ();
 
  /// Useful constructor 
  CglDuplicateRow (OsiSolverInterface * solver);
 
  /// Copy constructor 
  CglDuplicateRow (
    const CglDuplicateRow & rhs);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglDuplicateRow &
    operator=(
    const CglDuplicateRow& rhs);
  
  /// Destructor 
  virtual
    ~CglDuplicateRow ();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);

  /// This can be used to refresh any information
  virtual void refreshSolver(OsiSolverInterface * solver);
  //@}
      
protected:
  

  // Protected member data

  /**@name Protected member data */
  //@{
  /// Matrix
  CoinPackedMatrix matrix_;
  /// Matrix by row
  CoinPackedMatrix matrixByRow_; 
  /// Possible rhs (if 0 then not possible)
  int * rhs_;
  /// Marks duplicate rows
  int * duplicate_;
  /// To allow for <= rows
  int * lower_;
  /// Stored cuts if we found dominance cuts
  CglStored * storedCuts_;
  /// Check dominated columns if less than this number of candidates
  int maximumDominated_;
  /// Check duplicates if effective rhs <= this
  int maximumRhs_;
  /// Size of dynamic program
  int sizeDynamic_;
  /// 1 do rows, 2 do columns, 3 do both
  int mode_;
  /// Controls print out
  int logLevel_;
  //@}
};
#endif
