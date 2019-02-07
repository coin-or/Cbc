/* $Id$ */
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file CbcSolver.hpp
    \brief Defines CbcSolver, the proposed top-level class for the new-style
    cbc solver.

    This class is currently an orphan. With the removal of all code flagged
    with the NEW_STYLE_SOLVER, this class is never instantiated (and cannot
    be instantiated). It is available to be coopted as a top-level object
    wrapping the current CbcMain0 and CbcMain1, should that appear to be a
    desireable path forward.  -- lh, 091211 --
*/

#ifndef CbcSolver_H
#define CbcSolver_H

#include <string>
#include <vector>
#include "CoinMessageHandler.hpp"
#include "OsiClpSolverInterface.hpp"

#if CBC_OTHER_SOLVER == 1
#include "OsiCpxSolverInterface.hpp"
#endif

#include "CbcModel.hpp"
#include "CbcOrClpParam.hpp"

class CbcUser;
class CbcStopNow;
class CglCutGenerator;

//#############################################################################

/*! \brief This allows the use of the standalone solver in a flexible manner.

    It has an original OsiClpSolverInterface and CbcModel which it can use
    repeatedly, e.g., to get a heuristic solution and then start again.

    So I [jjf] will need a primitive scripting language which can then call
    solve and manipulate solution value and solution arrays.

    Also provides for user callback functions. Currently two ideas in
    gestation, CbcUser and CbcStopNow. The latter seems limited to deciding
    whether or not to stop. The former seems completely general, with a notion
    of importing and exporting, and a `solve', which should be interpreted as
    `do whatever this user function does'.

    Parameter initialisation is at last centralised in fillParameters().
*/

class CbcSolver {

public:
  ///@name Solve method
  //@{
  /** This takes a list of commands, does "stuff" and returns
        returnMode -
        0 model and solver untouched - babModel updated
        1 model updated - just with solution basis etc
        2 model updated i.e. as babModel (babModel NULL) (only use without preprocessing)
    */
  int solve(int argc, const char *argv[], int returnMode);
  /** This takes a list of commands, does "stuff" and returns
        returnMode -
        0 model and solver untouched - babModel updated
        1 model updated - just with solution basis etc
        2 model updated i.e. as babModel (babModel NULL) (only use without preprocessing)
    */
  int solve(const char *input, int returnMode);
  //@}
  ///@name Constructors and destructors etc
  //@{
  /// Default Constructor
  CbcSolver();

  /// Constructor from solver
  CbcSolver(const OsiClpSolverInterface &);

  /// Constructor from model
  CbcSolver(const CbcModel &);

  /** Copy constructor .
     */
  CbcSolver(const CbcSolver &rhs);

  /// Assignment operator
  CbcSolver &operator=(const CbcSolver &rhs);

  /// Destructor
  ~CbcSolver();
  /// Fill with standard parameters
  void fillParameters();
  /*! \brief Set default values in solvers from parameters

      Misleading. The current code actually reads default values from
      the underlying solvers and installs them as default values for a subset of
      parameters in #parameters_.
    */
  void fillValuesInSolver();
  /// Add user function
  void addUserFunction(CbcUser *function);
  /// Set user call back
  void setUserCallBack(CbcStopNow *function);
  /// Add cut generator
  void addCutGenerator(CglCutGenerator *generator);
  //@}
  ///@name miscellaneous methods to line up with old
  //@{
  // analyze model
  int *analyze(OsiClpSolverInterface *solverMod, int &numberChanged, double &increment,
    bool changeInt, CoinMessageHandler *generalMessageHandler);
  /** 1 - add heuristics to model
        2 - do heuristics (and set cutoff and best solution)
        3 - for miplib test so skip some
        (out model later)
    */
  //int doHeuristics(CbcModel * model, int type);
  /** Updates model_ from babModel_ according to returnMode
        returnMode -
        0 model and solver untouched - babModel updated
        1 model updated - just with solution basis etc
        2 model updated i.e. as babModel (babModel NULL) (only use without preprocessing)
    */
  void updateModel(ClpSimplex *model2, int returnMode);
  //@}
  ///@name useful stuff
  //@{
  /// Get int value
  int intValue(CbcOrClpParameterType type) const;
  /// Set int value
  void setIntValue(CbcOrClpParameterType type, int value);
  /// Get double value
  double doubleValue(CbcOrClpParameterType type) const;
  /// Set double value
  void setDoubleValue(CbcOrClpParameterType type, double value);
  /// User function (NULL if no match)
  CbcUser *userFunction(const char *name) const;
  /// Return original Cbc model
  inline CbcModel *model()
  {
    return &model_;
  }
  /// Return updated Cbc model
  inline CbcModel *babModel()
  {
    return babModel_;
  }
  /// Number of userFunctions
  inline int numberUserFunctions() const
  {
    return numberUserFunctions_;
  }
  /// User function array
  inline CbcUser **userFunctionArray() const
  {
    return userFunction_;
  }
  /// Copy of model on initial load (will contain output solutions)
  inline OsiClpSolverInterface *originalSolver() const
  {
    return originalSolver_;
  }
  /// Copy of model on initial load
  inline CoinModel *originalCoinModel() const
  {
    return originalCoinModel_;
  }
  /// Copy of model on initial load (will contain output solutions)
  void setOriginalSolver(OsiClpSolverInterface *originalSolver);
  /// Copy of model on initial load
  void setOriginalCoinModel(CoinModel *originalCoinModel);
  /// Number of cutgenerators
  inline int numberCutGenerators() const
  {
    return numberCutGenerators_;
  }
  /// Cut generator array
  inline CglCutGenerator **cutGeneratorArray() const
  {
    return cutGenerator_;
  }
  /// Start time
  inline double startTime() const
  {
    return startTime_;
  }
  /// Whether to print to std::cout
  inline void setPrinting(bool onOff)
  {
    noPrinting_ = !onOff;
  }
  /// Where to start reading commands
  inline void setReadMode(int value)
  {
    readMode_ = value;
  }
  //@}
private:
  ///@name Private member data
  //@{

  /// Reference model
  CbcModel model_;

  /// Updated model
  CbcModel *babModel_;

  /// User functions
  CbcUser **userFunction_;
  /** Status of user functions
        0 - not used
        1 - needs cbc_load
        2 - available - data in coinModel
        3 - data loaded - can do cbc_save
    */
  int *statusUserFunction_;
  /// Copy of model on initial load (will contain output solutions)
  OsiClpSolverInterface *originalSolver_;
  /// Copy of model on initial load
  CoinModel *originalCoinModel_;
  /// Cut generators
  CglCutGenerator **cutGenerator_;
  /// Number of user functions
  int numberUserFunctions_;
  /// Number of cut generators
  int numberCutGenerators_;
  /// Stop now stuff
  CbcStopNow *callBack_;
  /// Cpu time at instantiation
  double startTime_;
  /// Parameters and values
  std::vector< CbcOrClpParam > parameters_;
  /// Whether to do miplib test
  bool doMiplib_;
  /// Whether to print to std::cout
  bool noPrinting_;
  /// Where to start reading commands
  int readMode_;
  //@}
};
//#############################################################################

/// Structure to hold useful arrays
typedef struct {
  // Priorities
  int *priorities_;
  // SOS priorities
  int *sosPriority_;
  // Direction to branch first
  int *branchDirection_;
  // Input solution
  double *primalSolution_;
  // Down pseudo costs
  double *pseudoDown_;
  // Up pseudo costs
  double *pseudoUp_;
} CbcSolverUsefulData2;

//#############################################################################

/**
   The CbcSolver class was taken out at a 9/12/09 meeting
   This is a feeble replacement.
   At present everything is public
*/
class CbcSolverUsefulData {

public:
  ///@name Constructors and destructors etc
  //@{
  /// Default Constructor
  CbcSolverUsefulData();

  /** Copy constructor .
     */
  CbcSolverUsefulData(const CbcSolverUsefulData &rhs);

  /// Assignment operator
  CbcSolverUsefulData &operator=(const CbcSolverUsefulData &rhs);

  /// Destructor
  ~CbcSolverUsefulData();
  //@}

  ///@name Member data
  //@{
  // For time
  double totalTime_;
  // Parameters
  std::vector<CbcOrClpParam> parameters_;
  // Printing
  bool noPrinting_;
  // Whether to use signal handler
  bool useSignalHandler_;
  // Default pump tuning
  int initialPumpTune_;
  //@}
};
/// And this uses it
// When we want to load up CbcModel with options first
void CbcMain0(CbcModel &babSolver, CbcSolverUsefulData &solverData);
int CbcMain1(int argc, const char *argv[], CbcModel &babSolver, int(CbcModel *currentSolver, int whereFrom), CbcSolverUsefulData &solverData);

//#############################################################################

/*! \brief A class to allow the use of unknown user functionality

    For example, access to a modelling language (CbcAmpl).
*/
class CbcUser {

public:
  ///@name import/export methods
  //@{
  /*! \brief Import - gets full command arguments

      \return
      - -1 - no action
      -  0 - data read in without error
      -  1 - errors
    */
  virtual int importData(CbcSolver * /*model*/, int & /*argc*/, char ** /*argv[]*/)
  {
    return -1;
  }

  /*! \brief Export

      Values for mode:
      - 1 OsiClpSolver
      - 2 CbcModel
      - add 10 if infeasible from odd situation
    */
  virtual void exportSolution(CbcSolver * /*model*/,
    int /*mode*/, const char * /*message*/ = NULL) {}

  /// Export Data (i.e. at very end)
  virtual void exportData(CbcSolver * /*model*/) {}

  /// Get useful stuff
  virtual void fillInformation(CbcSolver * /*model*/,
    CbcSolverUsefulData & /*info*/) {}
  //@}

  ///@name usage methods
  //@{
  /// CoinModel if valid
  inline CoinModel *coinModel() const
  {
    return coinModel_;
  }
  /// Other info - needs expanding
  virtual void *stuff()
  {
    return NULL;
  }
  /// Name
  inline std::string name() const
  {
    return userName_;
  }
  /// Solve (whatever that means)
  virtual void solve(CbcSolver *model, const char *options) = 0;
  /// Returns true if function knows about option
  virtual bool canDo(const char *options) = 0;
  //@}

  ///@name Constructors and destructors etc
  //@{
  /// Default Constructor
  CbcUser();

  /// Copy constructor
  CbcUser(const CbcUser &rhs);

  /// Assignment operator
  CbcUser &operator=(const CbcUser &rhs);

  /// Clone
  virtual CbcUser *clone() const = 0;

  /// Destructor
  virtual ~CbcUser();
  //@}

protected:
  ///@name Private member data
  //@{

  /// CoinModel
  CoinModel *coinModel_;

  /// Name of user function
  std::string userName_;

  //@}
};
//#############################################################################

/*! \brief Support the use of a call back class to decide whether to stop

  Definitely under construction.
*/

class CbcStopNow {

public:
  ///@name Decision methods
  //@{
  /*! \brief Import

      Values for whereFrom:
       - 1 after initial solve by dualsimplex etc
       - 2 after preprocessing
       - 3 just before branchAndBound (so user can override)
       - 4 just after branchAndBound (before postprocessing)
       - 5 after postprocessing
       - 6 after a user called heuristic phase

      \return 0 if good
       nonzero return code to stop
    */
  virtual int callBack(CbcModel * /*currentSolver*/, int /*whereFrom*/)
  {
    return 0;
  }
  //@}

  ///@name Constructors and destructors etc
  //@{
  /// Default Constructor
  CbcStopNow();

  /** Copy constructor .
     */
  CbcStopNow(const CbcStopNow &rhs);

  /// Assignment operator
  CbcStopNow &operator=(const CbcStopNow &rhs);

  /// Clone
  virtual CbcStopNow *clone() const;

  /// Destructor
  virtual ~CbcStopNow();
  //@}

private:
  ///@name Private member data
  //@{
  //@}
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
