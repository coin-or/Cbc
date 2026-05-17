// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglPreProcess_H
#define CglPreProcess_H

#include <string>
#include <vector>

#include "CoinMessageHandler.hpp"
#include "OsiSolverInterface.hpp"
#include "CglStored.hpp"
#include "OsiPresolve.hpp"
#include "CglCutGenerator.hpp"

//#############################################################################

/** Class for preProcessing and postProcessing.

    While cuts can be added at any time in the tree, some cuts are actually just
    stronger versions of existing constraints.  In this case they can replace those
    constraints rather than being added as new constraints.  This is awkward in the
    tree but reasonable at the root node.

    This is a general process class which uses other cut generators to strengthen
    constraints, establish that constraints are redundant, fix variables and
    find relationships such as x + y == 1.

    Presolve will also be done.

    If row names existed they may be replaced by R0000000..., unless 
    setKeepColumnNames(true) is set.

*/

class CGLLIB_EXPORT CglPreProcess {

public:
  ///@name Main methods
  //@{
  /** preProcess problem - returning new problem.
      If makeEquality true then <= cliques converted to ==.
      Presolve will be done numberPasses times.

      Returns NULL if infeasible

      This version uses default strategy.  For more control copy and edit
      code from this function i.e. call preProcessNonDefault
  */
  OsiSolverInterface *preProcess(OsiSolverInterface &model,
    bool makeEquality = false, int numberPasses = 5);
  /** preProcess problem - returning new problem.
      If makeEquality true then <= cliques converted to ==.
      Presolve will be done numberPasses times.

      Returns NULL if infeasible

      This version assumes user has added cut generators to CglPreProcess object
      before calling it.  As an example use coding in preProcess
      If makeEquality is 1 add slacks to get cliques,
      if 2 add slacks to get sos (but only if looks plausible) and keep sos info
  */
  OsiSolverInterface *preProcessNonDefault(OsiSolverInterface &model,
    int makeEquality = 0, int numberPasses = 5,
    int tuning = 0);
  /** Creates solution in original model
      deleteStuff 0 - don't, 1 do (but not if infeasible), 2 always */
  void postProcess(OsiSolverInterface &model, int deleteStuff = 2);
  /** Tightens primal bounds to make dual and branch and cutfaster.  Unless
      fixed or integral, bounds are slightly looser than they could be.
      Returns negative if problem infeasible, number tightened if feasible
  */
  int tightenPrimalBounds(OsiSolverInterface &model,
			  bool tightenRowBounds=false,
			  double * scBound=NULL);
  /** Fix some of problem - returning new problem.
      Uses reduced costs.
      Optional signed character array
      1 always keep, -1 always discard, 0 use djs

  */
  OsiSolverInterface *someFixed(OsiSolverInterface &model,
    double fractionToKeep = 0.25,
    bool fixContinuousAsWell = false,
    char *keep = NULL) const;
  /** Replace cliques by more maximal cliques
      Returns NULL if rows not reduced by greater than cliquesNeeded*rows

  */
  OsiSolverInterface *cliqueIt(OsiSolverInterface &model,
    double cliquesNeeded = 0.0) const;
  /// If we have a cutoff - fix variables
  int reducedCostFix(OsiSolverInterface &model);
  //@}

  //---------------------------------------------------------------------------

  /**@name Parameter set/get methods

     The set methods return true if the parameter was set to the given value,
     false if the value of the parameter is out of range.

     The get methods return the value of the parameter.

  */
  //@{
  /** Set cutoff bound on the objective function.

    When using strict comparison, the bound is adjusted by a tolerance to
    avoid accidentally cutting off the optimal solution.
  */
  void setCutoff(double value);

  /// Get the cutoff bound on the objective function - always as minimize
  double getCutoff() const;
  /// The original solver associated with this model.
  inline OsiSolverInterface *originalModel() const
  {
    return originalModel_;
  }
  /// Set original model (probably acopy)
  inline void setOriginalModel(OsiSolverInterface *originalModel)
  {
    originalModel_ = originalModel;
  }
  /// Solver after making clique equalities (may == original)
  inline OsiSolverInterface *startModel() const
  {
    return startModel_;
  }
  /// Number of solvers
  inline int numberSolvers() const
  {
    return numberSolvers_;
  }
  /// Copies of solver at various stages after presolve
  inline OsiSolverInterface *modelAtPass(int iPass) const
  {
    if (iPass >= 0 && iPass < numberSolvers_)
      return model_[iPass];
    else
      return NULL;
  }
  /// Copies of solver at various stages after presolve after modifications
  inline OsiSolverInterface *modifiedModel(int iPass) const
  {
    if (iPass >= 0 && iPass < numberSolvers_)
      return modifiedModel_[iPass];
    else
      return NULL;
  }
  /// Matching presolve information
  inline OsiPresolve *presolve(int iPass) const
  {
    if (iPass >= 0 && iPass < numberSolvers_)
      return presolve_[iPass];
    else
      return NULL;
  }
  /// Set matching presolve information
  inline void setPresolve(int iPass, OsiPresolve *presolve)
  {
    if (iPass >= 0 && iPass < numberSolvers_)
      presolve_[iPass] = presolve;
  }
  /** Return a pointer to the original columns (with possible  clique slacks)
      MUST be called before postProcess otherwise you just get 0,1,2.. */
  const int *originalColumns();
  /** Return a pointer to the original rows
      MUST be called before postProcess otherwise you just get 0,1,2.. */
  const int *originalRows();
  /// Number of SOS if found
  inline int numberSOS() const
  {
    return numberSOS_;
  }
  /// Type of each SOS
  inline const int *typeSOS() const
  {
    return typeSOS_;
  }
  /// Start of each SOS
  inline const int *startSOS() const
  {
    return startSOS_;
  }
  /// Columns in SOS
  inline const int *whichSOS() const
  {
    return whichSOS_;
  }
  /// Weights for each SOS column
  inline const double *weightSOS() const
  {
    return weightSOS_;
  }
  /// Pass in prohibited columns
  void passInProhibited(const char *prohibited, int numberColumns);
  /// Updated prohibited columns
  inline const char *prohibited()
  {
    return prohibited_;
  }
  /// Number of iterations PreProcessing
  inline int numberIterationsPre() const
  {
    return numberIterationsPre_;
  }
  /// Number of iterations PostProcessing
  inline int numberIterationsPost() const
  {
    return numberIterationsPost_;
  }
  /** Pass in row types
      -2 do not touch row
      0 normal
      1 cut rows - will be dropped if remain in
      At end of preprocess cut rows will be dropped
      and put into cuts
  */
  void passInRowTypes(const char *rowTypes, int numberRows);
  /** Updated row types - may be NULL
      Carried around and corresponds to existing rows
      -1 added by preprocess e.g. x+y=1
      0 normal
      1 cut rows - can be dropped if wanted
  */
  inline const char *rowTypes()
  {
    return rowType_;
  }
  /// Return cuts from dropped rows
  inline const CglStored &cuts() const
  {
    return cuts_;
  }
  /// Return pointer to cuts from dropped rows
  inline const CglStored *cutsPointer() const
  {
    return &cuts_;
  }
  /// Update prohibited and rowType (and SC stuff)
  void update(const OsiPresolve *pinfo,
	      const OsiSolverInterface *solver,double * scBound);
  /// Get options
  inline int options() const
  {
    return options_;
  }
  /// Set options
  inline void setOptions(int value)
  {
    options_ = value;
  }
  //@}

  ///@name Cut generator methods
  //@{
  /// Get the number of cut generators
  inline int numberCutGenerators() const
  {
    return numberCutGenerators_;
  }
  /// Get the list of cut generators
  inline CglCutGenerator **cutGenerators() const
  {
    return generator_;
  }
  ///Get the specified cut generator
  inline CglCutGenerator *cutGenerator(int i) const
  {
    return generator_[i];
  }
  /** Add one generator - up to user to delete generators.
  */
  void addCutGenerator(CglCutGenerator *generator);
  //@}

  /**@name Setting/Accessing application data */
  //@{
  /** Set application data.

	This is a pointer that the application can store into and
	retrieve.
	This field is available for the application to optionally
	define and use.
    */
  void setApplicationData(void *appData);

  /// Get application data
  void *getApplicationData() const;
  //@}

  //---------------------------------------------------------------------------

  /**@name Message handling */
  //@{
  /// Pass in Message handler (not deleted at end)
  void passInMessageHandler(CoinMessageHandler *handler);
  /// Set language
  void newLanguage(CoinMessages::Language language);
  inline void setLanguage(CoinMessages::Language language)
  {
    newLanguage(language);
  }
  /// Return handler
  inline CoinMessageHandler *messageHandler() const
  {
    return handler_;
  }
  /// Return messages
  inline CoinMessages messages()
  {
    return messages_;
  }
  /// Return pointer to messages
  inline CoinMessages *messagesPointer()
  {
    return &messages_;
  }
  //@}
  //---------------------------------------------------------------------------

  ///@name Constructors and destructors etc
  //@{
  /// Constructor
  CglPreProcess();

  /// Copy constructor .
  CglPreProcess(const CglPreProcess &rhs);

  /// Assignment operator
  CglPreProcess &operator=(const CglPreProcess &rhs);

  /// Destructor
  ~CglPreProcess();

  /// Clears out as much as possible
  void gutsOfDestructor();

  /// Gets rid of models
  void clean();

  /// Set time limit
  void setTimeLimit(const double timeLimit, const bool useElapsedTime);

  /// Keeps original column names
  void setKeepColumnNames(const bool keep);

  //@}
private:
  ///@name private methods
  //@{
  /** Return model with useful modifications.  
      If constraints true then adds any x+y=1 or x-y=0 constraints
      If NULL infeasible
  */
  OsiSolverInterface *modified(OsiSolverInterface *model,
    bool constraints,
    int &numberChanges,
    int iBigPass,
    int numberPasses);
  /// create original columns and rows
  void createOriginalIndices();
  /// Make continuous variables integer
  void makeInteger();
  //@}

  //---------------------------------------------------------------------------

private:
  ///@name Private member data
  //@{

  /// The original solver associated with this model.
  OsiSolverInterface *originalModel_;
  /// Solver after making clique equalities (may == original)
  OsiSolverInterface *startModel_;
  /// Number of solvers at various stages
  int numberSolvers_;
  /// Copies of solver at various stages after presolve
  OsiSolverInterface **model_;
  /// Copies of solver at various stages after presolve after modifications
  OsiSolverInterface **modifiedModel_;
  /// Matching presolve information
  OsiPresolve **presolve_;

  /// Message handler
  CoinMessageHandler *handler_;

  /** Flag to say if handler_ is the default handler.
  
    The default handler is deleted when the model is deleted. Other
    handlers (supplied by the client) will not be deleted.
  */
  bool defaultHandler_;

  /// Cgl messages
  CoinMessages messages_;

  /// Pointer to user-defined data structure
  void *appData_;
  /// Original column numbers
  int *originalColumn_;
  /// Original row numbers
  int *originalRow_;
  /// Number of cut generators
  int numberCutGenerators_;
  /// Cut generators
  CglCutGenerator **generator_;
  /// Number of SOS if found
  int numberSOS_;
  /// Type of each SOS
  int *typeSOS_;
  /// Start of each SOS
  int *startSOS_;
  /// Columns in SOS
  int *whichSOS_;
  /// Weights for each SOS column
  double *weightSOS_;
  /// Number of columns in original prohibition set
  int numberProhibited_;
  /// Number of iterations done in PreProcessing
  int numberIterationsPre_;
  /// Number of iterations done in PostProcessing
  int numberIterationsPost_;
  /** Columns which should not be presolved e.g. SOS
    2 says really really leave alone e.g. sc or lot-sizing */
  char *prohibited_;
  /// Number of rows in original row types
  int numberRowType_;
  /** Options
      1 - original model had integer bounds before tightening
      2 - don't do probing
      4 - don't do duplicate rows
      8 - don't do cliques
      16 - some heavy probing options
      64 - very heavy probing
      128 - appData_ has SC data
      256 - clumsy postprocessing
      512 - try for dominated columns
  */
  int options_;
  /** Row types (may be NULL) 
      Carried around and corresponds to existing rows
      -1 added by preprocess e.g. x+y=1
      0 normal
      1 cut rows - can be dropped if wanted
  */
  char *rowType_;
  /// Cuts from dropped rows
  CglStored cuts_;

  /// use elapsed (wallclock time) or cpu time
  bool useElapsedTime_;

  /// time limit (default COIN_DBL_MAX)
  double timeLimit_;

  /// keep column names
  bool keepColumnNames_;

  /// current elapsed or cpu time
  double getCurrentCPUTime() const;

  //@}
};

/// For Bron-Kerbosch
class CGLLIB_EXPORT CglBK {

public:
  ///@name Main methods
  //@{
  /// For recursive Bron-Kerbosch
  void bronKerbosch();
  /// Creates strengthened smaller model
  OsiSolverInterface *newSolver(const OsiSolverInterface &model);
  //@}

  //---------------------------------------------------------------------------

  /**@name Parameter set/get methods

     The set methods return true if the parameter was set to the given value,
     false if the value of the parameter is out of range.

     The get methods return the value of the parameter.

  */
  //@{
  //@}

  //---------------------------------------------------------------------------

  ///@name Constructors and destructors etc
  //@{
  /// Default constructor
  CglBK();

  /// Useful constructor
  CglBK(const OsiSolverInterface &model, const char *rowType,
    int numberElements);

  /// Copy constructor .
  CglBK(const CglBK &rhs);

  /// Assignment operator
  CglBK &operator=(const CglBK &rhs);

  /// Destructor
  ~CglBK();

  //@}

  //---------------------------------------------------------------------------

private:
  ///@name Private member data
  //@{
  /// Current candidates (created at each level)
  int *candidates_;
  /// Array to mark stuff
  char *mark_;
  /// Starts for graph (numberPossible+1)
  CoinBigIndex *start_;
  /// Other column/node
  int *otherColumn_;
  /// Original row (in parallel with otherColumn_)
  int *originalRow_;
  /// How many times each original row dominated
  int *dominated_;
  /// Clique entries
  CoinPackedMatrix *cliqueMatrix_;
  /// points to row types
  const char *rowType_;
  /// Number of original columns
  int numberColumns_;
  /// Number of original rows
  int numberRows_;
  /// Number possible
  int numberPossible_;
  /// Current number of candidates
  int numberCandidates_;
  /// First not (stored backwards from numberPossible_)
  int firstNot_;
  /// Current number in clique
  int numberIn_;
  /// For acceleration
  int left_;
  int lastColumn_;
  //@}
};
/**
   Only store unique row cuts
*/
// for hashing
typedef struct {
  int index, next;
} CglHashLink;
class OsiRowCut;
class CGLLIB_EXPORT CglUniqueRowCuts {
public:
  CglUniqueRowCuts(int initialMaxSize = 0, int hashMultiplier = 4);
  ~CglUniqueRowCuts();
  CglUniqueRowCuts(const CglUniqueRowCuts &rhs);
  CglUniqueRowCuts &operator=(const CglUniqueRowCuts &rhs);
  inline OsiRowCut *cut(int sequence) const
  {
    return rowCut_[sequence];
  }
  inline int numberCuts() const
  {
    return numberCuts_;
  }
  inline int sizeRowCuts() const
  {
    return numberCuts_;
  }
  inline OsiRowCut *rowCutPtr(int sequence)
  {
    return rowCut_[sequence];
  }
  void eraseRowCut(int sequence);
  // insert cut
  inline void insert(const OsiRowCut &cut)
  {
    insertIfNotDuplicate(cut);
  }
  // Return 0 if added, 1 if not
  int insertIfNotDuplicate(const OsiRowCut &cut);
  // Add in cuts as normal cuts (and delete)
  void addCuts(OsiCuts &cs);

private:
  OsiRowCut **rowCut_;
  /// Hash table
  CglHashLink *hash_;
  int size_;
  int hashMultiplier_;
  int numberCuts_;
  int lastHash_;
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
