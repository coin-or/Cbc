// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglProbing_H
#define CglProbing_H

#include <string>

#include "CglCutGenerator.hpp"
/** Only useful type of disaggregation is most normal
    For now just done for 0-1 variables
    Can be used for building cliques
 */
typedef struct {
  // unsigned int zeroOne:1; // nonzero if affected variable is 0-1
  // unsigned int whenAtUB:1; // nonzero if fixing happens when this variable at 1
  // unsigned int affectedToUB:1; // nonzero if affected variable fixed to UB
  // unsigned int affected:29; // If 0-1 then 0-1 sequence, otherwise true
  unsigned int affected;
} disaggregationAction;

/** Probing Cut Generator Class */
class CGLLIB_EXPORT CglProbing : public CglCutGenerator {
  friend CGLLIB_EXPORT void CglProbingUnitTest(const OsiSolverInterface *siP,
    const std::string mpdDir);

public:
  /**@name Generate Cuts */
  //@{
  /** Generate probing/disaggregation cuts for the model of the
      solver interface, si.

      This is a simplification of probing ideas put into OSL about
      ten years ago.  The only known documentation is a copy of a
      talk handout - we think Robin Lougee-Heimer has a copy!

      For selected integer variables (e.g. unsatisfied ones) the effect of
      setting them up or down is investigated.  Setting a variable up
      may in turn set other variables (continuous as well as integer).
      There are various possible results:

      1) It is shown that problem is infeasible (this may also be
         because objective function or reduced costs show worse than
         best solution).  If the other way is feasible we can generate
         a column cut (and continue probing), if not feasible we can
         say problem infeasible.

      2) If both ways are feasible, it can happen that x to 0 implies y to 1
         ** and x to 1 implies y to 1 (again a column cut).  More common
         is that x to 0 implies y to 1 and x to 1 implies y to 0 so we could
         substitute for y which might lead later to more powerful cuts.
         ** This is not done in this code as there is no mechanism for
         returning information.

      3) When x to 1 a constraint went slack by c.  We can tighten the
         constraint ax + .... <= b (where a may be zero) to
         (a+c)x + .... <= b.  If this cut is violated then it is
         generated.

      4) Similarly we can generate implied disaggregation cuts

      Note - differences to cuts in OSL.

      a) OSL had structures intended to make this faster.
      b) The "chaining" in 2) was done
      c) Row cuts modified original constraint rather than adding cut
      b) This code can cope with general integer variables.

      Insert the generated cuts into OsiCut, cs.

      If a "snapshot" of a matrix exists then this will be used.
      Presumably this will give global cuts and will be faster.
      No check is done to see if cuts will be global.

      Otherwise use current matrix.

      Both row cuts and column cuts may be returned

      The mode options are:
      0) Only unsatisfied integer variables will be looked at.
         If no information exists for that variable then
         probing will be done so as a by-product you "may" get a fixing
         or infeasibility.  This will be fast and is only available
         if a snapshot exists (otherwise as 1).
         The bounds in the snapshot are the ones used.
      1) Look at unsatisfied integer variables, using current bounds.
         Probing will be done on all looked at.
      2) Look at all integer variables, using current bounds.
         Probing will be done on all

         ** If generateCutsAndModify is used then new relaxed
         row bounds and tightened column bounds are generated
         Returns number of infeasibilities
  */
  virtual void generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
    const CglTreeInfo info = CglTreeInfo());
  int generateCutsAndModify(const OsiSolverInterface &si, OsiCuts &cs,
    CglTreeInfo *info);
  //@}

  /**@name snapshot etc */
  //@{
  /** Create a copy of matrix which is to be used
      this is to speed up process and to give global cuts
      Can give an array with 1 set to select, 0 to ignore
      column bounds are tightened
      If array given then values of 1 will be set to 0 if redundant.
      Objective may be added as constraint
      Returns 1 if infeasible otherwise 0
  */
  int snapshot(const OsiSolverInterface &si,
    char *possible = NULL,
    bool withObjective = true);
  /// Deletes snapshot
  void deleteSnapshot();
  /** Creates cliques for use by probing.
      Only cliques >= minimumSize and < maximumSize created
      Can also try and extend cliques as a result of probing (root node).
      Returns number of cliques found.
  */
  int createCliques(OsiSolverInterface &si,
    int minimumSize = 2, int maximumSize = 100);
  /// Delete all clique information
  void deleteCliques();
  /** Create a fake model by adding cliques
      if type&4 then delete rest of model first,
      if 1 then add proper cliques, 2 add fake cliques */
  OsiSolverInterface *cliqueModel(const OsiSolverInterface *model,
    int type);
  //@}

  /**@name Get tighter column bounds */
  //@{
  /// Lower
  const double *tightLower() const;
  /// Upper
  const double *tightUpper() const;
  /// Array which says tighten continuous
  const char *tightenBounds() const
  {
    return tightenBounds_;
  }
  //@}

  /**@name Get possible freed up row bounds - only valid after mode==3 */
  //@{
  /// Lower
  const double *relaxedRowLower() const;
  /// Upper
  const double *relaxedRowUpper() const;
  //@}

  /**@name Change mode */
  //@{
  /// Set
  void setMode(int mode);
  /// Get
  int getMode() const;
  //@}

  /**@name Change maxima */
  //@{
  /// Set maximum number of passes per node
  void setMaxPass(int value);
  /// Get maximum number of passes per node
  int getMaxPass() const;
  /// Set log level - 0 none, 1 - a bit, 2 - more details
  void setLogLevel(int value);
  /// Get log level
  int getLogLevel() const;
  /// Set maximum number of unsatisfied variables to look at
  void setMaxProbe(int value);
  /// Get maximum number of unsatisfied variables to look at
  int getMaxProbe() const;
  /// Set maximum number of variables to look at in one probe
  void setMaxLook(int value);
  /// Get maximum number of variables to look at in one probe
  int getMaxLook() const;
  /// Set maximum number of elements in row for it to be considered
  void setMaxElements(int value);
  /// Get maximum number of elements in row for it to be considered
  int getMaxElements() const;
  /// Set maximum number of passes per node  (root node)
  void setMaxPassRoot(int value);
  /// Get maximum number of passes per node (root node)
  int getMaxPassRoot() const;
  /// Set maximum number of unsatisfied variables to look at (root node)
  void setMaxProbeRoot(int value);
  /// Get maximum number of unsatisfied variables to look at (root node)
  int getMaxProbeRoot() const;
  /// Set maximum number of variables to look at in one probe (root node)
  void setMaxLookRoot(int value);
  /// Get maximum number of variables to look at in one probe (root node)
  int getMaxLookRoot() const;
  /// Set maximum number of elements in row for it to be considered (root node)
  void setMaxElementsRoot(int value);
  /// Get maximum number of elements in row for it to be considered (root node)
  int getMaxElementsRoot() const;
  /**
     Returns true if may generate Row cuts in tree (rather than root node).
     Used so know if matrix will change in tree.  Really
     meant so column cut generators can still be active
     without worrying code.
     Default is true
  */
  virtual bool mayGenerateRowCutsInTree() const;
  //@}

  /**@name Get information back from probing */
  //@{
  /// Number looked at this time
  inline int numberThisTime() const
  {
    return numberThisTime_;
  }
  /// Which ones looked at this time
  inline const int *lookedAt() const
  {
    return lookedAt_;
  }
  /// Modified upper bounds (if set)
  inline const double *colUpper() const
  {
    return colUpper_;
  }
  /// Modified lower bounds (if set)
  inline const double *colLower() const
  {
    return colLower_;
  }
  //@}

  /**@name Stop or restart row cuts (otherwise just fixing from probing) */
  //@{
  /// Set
  /// 0 no cuts, 1 just disaggregation type, 2 coefficient ( 3 both)
  void setRowCuts(int type);
  /// Get
  int rowCuts() const;
  //@}
  /// Clique type
  typedef struct {
    unsigned int equality : 1; //  nonzero if clique is ==
  } CliqueType;

  /**@name Information on cliques */
  //@{
  /// Number of cliques
  inline int numberCliques() const
  {
    return numberCliques_;
  }
  /// Clique type
  inline CliqueType *cliqueType() const
  {
    return cliqueType_;
  }
  /// Start of each clique
  inline CoinBigIndex *cliqueStart() const
  {
    return cliqueStart_;
  }
  /// Entries for clique
  inline CliqueEntry *cliqueEntry() const
  {
    return cliqueEntry_;
  }
  /// Changed up/down (if used)
  inline int *endFixStart() const
  {
    return endFixStart_;
  }
  //@}

  /**@name Whether use objective as constraint */
  //@{
  /** Set
      0 don't
      1 do
      -1 don't even think about it
  */
  void setUsingObjective(int yesNo);
  /// Get
  int getUsingObjective() const;
  //@}

  /**@name Wall-clock time limit for a single generateCuts() call */
  //@{
  /// Set maximum wall-clock seconds allowed (0 = no limit)
  void setMaxSeconds(double value) { maxSeconds_ = value; }
  /// Get maximum wall-clock seconds allowed
  double getMaxSeconds() const { return maxSeconds_; }
  //@}

  /**@name Mark which continuous variables are to be tightened */
  //@{
  /// Mark variables to be tightened
  void tightenThese(const OsiSolverInterface &solver, int number, const int *which);
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor
  CglProbing();

  /// Copy constructor
  CglProbing(
    const CglProbing &);

  /// Clone
  virtual CglCutGenerator *clone() const;

  /// Assignment operator
  CglProbing &
  operator=(
    const CglProbing &rhs);

  /// Destructor
  virtual ~CglProbing();

  /// This can be used to refresh any inforamtion
  virtual void refreshSolver(OsiSolverInterface *solver);
  /// Create C++ lines to get to current state
  virtual std::string generateCpp(FILE *fp);
  //@}

private:
  // Private member methods
  /**@name probe */
  //@{
  /// Does probing and adding cuts (without cliques and mode_!=0)
  int probe(const OsiSolverInterface &si,
    const OsiRowCutDebugger *debugger,
    OsiCuts &cs,
    double *colLower, double *colUpper, CoinPackedMatrix *rowCopy,
    CoinPackedMatrix *columnCopy, const CoinBigIndex *rowStartPos,
    const int *realRow, const double *rowLower, const double *rowUpper,
    const char *intVar, double *minR, double *maxR, int *markR,
    CglTreeInfo *info);
  /// Does probing and adding cuts (without cliques and mode_!=0)
  int probeFast(const OsiSolverInterface &si,
    const OsiRowCutDebugger *debugger,
    OsiCuts &cs,
    double *colLower, double *colUpper, CoinPackedMatrix *rowCopy,
    CoinPackedMatrix *columnCopy, const CoinBigIndex *rowStartPos,
    const int *realRow, const double *rowLower, const double *rowUpper,
    const char *intVar, double *minR, double *maxR, int *markR,
    CglTreeInfo *info);
  /// Does probing and adding cuts (with cliques)
  int probeCliques(const OsiSolverInterface &si,
    const OsiRowCutDebugger *debugger,
    OsiCuts &cs,
    double *colLower, double *colUpper, CoinPackedMatrix *rowCopy,
    CoinPackedMatrix *columnCopy, const int *realRow,
    double *rowLower, double *rowUpper,
    char *intVar, double *minR, double *maxR, int *markR,
    CglTreeInfo *info);
  /// Does probing and adding cuts for clique slacks
  int probeSlacks(const OsiSolverInterface &si,
    const OsiRowCutDebugger *debugger,
    OsiCuts &cs,
    double *colLower, double *colUpper, CoinPackedMatrix *rowCopy,
    CoinPackedMatrix *columnCopy,
    double *rowLower, double *rowUpper,
    char *intVar, double *minR, double *maxR, int *markR,
    CglTreeInfo *info);
  /** Does most of work of generateCuts
      Returns number of infeasibilities */
  int gutsOfGenerateCuts(const OsiSolverInterface &si,
    OsiCuts &cs,
    double *rowLower, double *rowUpper,
    double *colLower, double *colUpper,
    CglTreeInfo *info);
  /// Sets up clique information for each row
  void setupRowCliqueInformation(const OsiSolverInterface &si);
  /** This tightens column bounds (and can declare infeasibility)
      It may also declare rows to be redundant */
  int tighten(double *colLower, double *colUpper,
    const int *column, const double *rowElements,
    const CoinBigIndex *rowStart, const CoinBigIndex *rowStartPos,
    const int *rowLength,
    const CoinPackedMatrix *,
    double *rowLower, double *rowUpper,
    int nRows, int nCols, const char *intVar, int maxpass,
    double tolerance);
  /// This just sets minima and maxima on rows
  void tighten2(double *colLower, double *colUpper,
    const int *column, const double *rowElements,
    const CoinBigIndex *rowStart,
    const int *rowLength,
    double *rowLower, double *rowUpper,
    double *minR, double *maxR, int *markR,
    int nRows);
  //@}

  // Private member data

  struct disaggregation_struct_tag;
  friend struct CglProbing::disaggregation_struct_tag;

  /**@name Private member data */
  //@{
  /// Row copy (only if snapshot)
  CoinPackedMatrix *rowCopy_;
  /// Column copy (only if snapshot)
  CoinPackedMatrix *columnCopy_;
  /// Lower bounds on rows
  double *rowLower_;
  /// Upper bounds on rows
  double *rowUpper_;
  /// Lower bounds on columns
  double *colLower_;
  /// Upper bounds on columns
  double *colUpper_;
  /// Pass in from probe to tighten
  double *minR_;
  double *maxR_;
  /// If row looked at
  char *markRow_;
  /// If column looked at and interesting
  char *markColumn_;
  /// Columns to update
  int *lookColumn_;
  /// Rows to look at (alternate passes)
  int *lookRow_[2];
  /// Number of rows in snapshot (or when cliqueRow stuff computed)
  int numberRows_;
  /// Number of columns in problem ( must == current)
  int numberColumns_;
  /// Tolerance to see if infeasible
  double primalTolerance_;
  /** Mode - 0 lazy using snapshot, 1 just unsatisfied, 2 all.
      16 bit set if want to extend cliques at root node
  */
  int mode_;
  /** Row cuts flag
      0 no cuts, 1 just disaggregation type, 2 coefficient ( 3 both), 4 just column cuts
      -n as +n but just fixes variables unless at root
  */
  int rowCuts_;
  /// Maximum number of passes to do in probing
  int maxPass_;
  /// Log level - 0 none, 1 - a bit, 2 - more details
  int logLevel_;
  /// Maximum number of unsatisfied variables to probe
  int maxProbe_;
  /// Maximum number of variables to look at in one probe
  int maxStack_;
  /// Maximum number of elements in row for scan
  int maxElements_;
  /// Maximum number of passes to do in probing at root
  int maxPassRoot_;
  /// Maximum number of unsatisfied variables to probe at root
  int maxProbeRoot_;
  /// Maximum number of variables to look at in one probe at root
  int maxStackRoot_;
  /// Maximum number of elements in row for scan at root
  int maxElementsRoot_;
  /// Whether to include objective as constraint
  int usingObjective_;
  /// Maximum wall-clock seconds per generateCuts() call (0 = no limit)
  double maxSeconds_;
  /// Number of integer variables
  int numberIntegers_;
  /// Number of 0-1 integer variables
  int number01Integers_;
  /// Number looked at this time
  int numberThisTime_;
  /// Total number of times called
  int totalTimesCalled_;
  /// Which ones looked at this time
  int *lookedAt_;
  /// Disaggregation cuts and for building cliques
  typedef struct disaggregation_struct_tag {
    int sequence; // integer variable
    // index will be NULL if no probing done yet
    int length; // length of newValue
    disaggregationAction *index; // columns whose bounds will be changed
  } disaggregation;
  disaggregation *cutVector_;
  /// Cliques
  /// Number of cliques
  int numberCliques_;
  /// Clique type
  CliqueType *cliqueType_;
  /// Start of each clique
  CoinBigIndex *cliqueStart_;
  /// Entries for clique
  CliqueEntry *cliqueEntry_;
  /** Start of oneFixes cliques for a column in matrix or -1 if not
      in any clique */
  int *oneFixStart_;
  /** Start of zeroFixes cliques for a column in matrix or -1 if not
      in any clique */
  int *zeroFixStart_;
  /// End of fixes for a column
  int *endFixStart_;
  /// Clique numbers for one or zero fixes
  int *whichClique_;
  /** For each column with nonzero in row copy this gives a clique "number".
      So first clique mentioned in row is always 0.  If no entries for row
      then no cliques.  If sequence > numberColumns then not in clique.
  */
  CliqueEntry *cliqueRow_;
  /// cliqueRow_ starts for each row
  int *cliqueRowStart_;
  /// If not null and [i] !=0 then also tighten even if continuous
  char *tightenBounds_;
  //@}
};
inline int affectedInDisaggregation(const disaggregationAction &dis)
{
  return dis.affected & 0x1fffffff;
}
inline void setAffectedInDisaggregation(disaggregationAction &dis,
  int affected)
{
  dis.affected = affected | (dis.affected & 0xe0000000);
}
#ifdef NDEBUG
inline bool zeroOneInDisaggregation(const disaggregationAction &)
{
  return true;
}
#else
inline bool zeroOneInDisaggregation(const disaggregationAction &dis)
//{ return (dis.affected&0x80000000)!=0;}
{
  assert((dis.affected & 0x80000000) != 0);
  return true;
}
#endif
inline void setZeroOneInDisaggregation(disaggregationAction &dis, bool zeroOne)
{
  dis.affected = (zeroOne ? 0x80000000 : 0) | (dis.affected & 0x7fffffff);
}
inline bool whenAtUBInDisaggregation(const disaggregationAction &dis)
{
  return (dis.affected & 0x40000000) != 0;
}
inline void setWhenAtUBInDisaggregation(disaggregationAction &dis, bool whenAtUB)
{
  dis.affected = (whenAtUB ? 0x40000000 : 0) | (dis.affected & 0xbfffffff);
}
inline bool affectedToUBInDisaggregation(const disaggregationAction &dis)
{
  return (dis.affected & 0x20000000) != 0;
}
inline void setAffectedToUBInDisaggregation(disaggregationAction &dis, bool affectedToUB)
{
  dis.affected = (affectedToUB ? 0x20000000 : 0) | (dis.affected & 0xdfffffff);
}

// #############################################################################
/** A function that tests the methods in the CglProbing class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
CGLLIB_EXPORT
void CglProbingUnitTest(const OsiSolverInterface *siP,
  const std::string mpdDir);

/// This just uses implication info
class CGLLIB_EXPORT CglImplication : public CglCutGenerator {

public:
  /**@name Generate Cuts */
  //@{
  /** Generate cuts from implication table
  Insert generated cuts into the cut set cs.
  */
  virtual void generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
    const CglTreeInfo info = CglTreeInfo());
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor
  CglImplication();

  /// Constructor with info
  CglImplication(CglTreeProbingInfo *info);

  /// Copy constructor
  CglImplication(
    const CglImplication &);

  /// Clone
  virtual CglCutGenerator *clone() const;

  /// Assignment operator
  CglImplication &
  operator=(
    const CglImplication &rhs);

  /// Destructor
  virtual ~CglImplication();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp(FILE *fp);
  //@}
  /**@name Set implication */
  //@{
  /// Set implication
  inline void setProbingInfo(CglTreeProbingInfo *info)
  {
    probingInfo_ = info;
  }
  //@}

private:
  /**@name Private member data */
  //@{
  /// Pointer to tree probing info
  CglTreeProbingInfo *probingInfo_;
  //@}
};
#endif
