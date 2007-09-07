// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcSolver_H
#define CbcSolver_H
#include <string>
#include <vector>
#include "CoinFinite.hpp"
#include "CoinMessageHandler.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcOrClpParam.hpp"
class CbcUser;
class CbcStopNow;

//#############################################################################

/** This allows the use of the standalone solver in a flexible manner
    It has an original OsiClpSolverInterface and CbcModel which
    it can use repeatedly e.g. get a heuristic solution and then start again

    So I will need a primitive scripting language which can then call solve
    and manipulate solution value and solution arrays.
    
*/

class CbcSolver  {
  
public:
  ///@name Solve method 
  //@{
  /** This takes a list of commands, does "stuff" and returns 
      returnMode - 
      0 model and solver untouched - babModel updated
      1 model updated - just with solution basis etc
      2 model updated i.e. as babModel (babModel NULL)
  */
  int solve(int argc, const char * argv[], int returnMode); 
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
  CbcSolver(const CbcSolver & rhs);
  
  /// Assignment operator 
  CbcSolver & operator=(const CbcSolver& rhs);
  
  /// Destructor 
  ~CbcSolver ();
  /// Fill with standard parameters
  void fillParameters();
  /// Set default values in solvers from parameters
  void fillValuesInSolver();
  //@}
  ///@name miscellaneous methods to line up with old
  //@{
  // analyze model
  int * analyze(OsiClpSolverInterface * solverMod, int & numberChanged, double & increment,
		bool changeInt,  CoinMessageHandler * generalMessageHandler);
  /** 1 - add heuristics to model
      2 - do heuristics (and set cutoff and best solution)
      3 - for miplib test so skip some
      (out model later)
  */
  int doHeuristics(CbcModel * model, int type);
  /** Updates model_ from babModel_ according to returnMode
      returnMode - 
      0 model and solver untouched - babModel updated
      1 model updated - just with solution basis etc
      2 model updated i.e. as babModel (babModel NULL)
  */
  void updateModel(ClpSimplex * model2, int returnMode);
  //@}
  ///@name useful stuff
  //@{
  /// Get int value
  int intValue(CbcOrClpParameterType type) const;
  /// Set int value
  void setIntValue(CbcOrClpParameterType type,int value);
  /// Get double value
  double doubleValue(CbcOrClpParameterType type) const;
  /// Set double value
  void setDoubleValue(CbcOrClpParameterType type,double value);
  /// User function (NULL if no match)
  CbcUser * userFunction(const char * name) const;
  inline CbcModel * model() 
  { return &model_;}
  //@}
private:
  ///@name Private member data 
  //@{
  
  /// Reference model
  CbcModel model_;
  
  /// Updated model
  CbcModel * babModel_;
  
  /// User functions
  CbcUser ** userFunction_;
  /// Number of user functions
  int numberUserFunctions_;
  /// Stop now stuff
  CbcStopNow * callBack_;
  /// Cpu time at instantiation
  double startTime_;
  /// Parameters and values
  CbcOrClpParam * parameters_;
  /// Number of parameters
  int numberParameters_ ;
  /// Whether to do miplib test
  bool doMiplib_;
  /// Whether to print to std::cout
  bool noPrinting_;
  
  //@}
};
//#############################################################################

/** This allows the use of an unknown user stuff including modeling languages
 */

class CbcUser  {
  
public:
  ///@name import/export methods 
  //@{
  /// Import - 0 if good
  virtual int importData(std::string fileName) {return 0;}
  /// Export
  virtual void exportSolution() {}
  //@}
  ///@name usage methods 
  //@{
  /// CoinModel if valid
  inline CoinModel * coinModel() const
  { return coinModel_;}
  /// Other info - needs expanding
  virtual void * stuff() {return NULL;}
  /// Name
  inline std::string name() const
  { return userName_;}
  /// Solve (whatever that means)
  void solve(CbcSolver * model, const char * options) {}
  //@}
  ///@name Constructors and destructors etc
  //@{
  /// Default Constructor
  CbcUser(); 
  
  /** Copy constructor .
   */  
  CbcUser(const CbcUser & rhs);
  
  /// Assignment operator 
  CbcUser & operator=(const CbcUser& rhs);

  /// Clone
  CbcUser * clone() const;
  
  /// Destructor 
  virtual ~CbcUser ();
  //@}
  
protected:
  ///@name Private member data 
  //@{
  
  /// CoinModel
  CoinModel * coinModel_;
  
  /// Name of user function
  std::string userName_;

 //@}
};
//#############################################################################

/** This allows the use of a call back class to decide whether to stop
 */

class CbcStopNow  {
  
public:
  ///@name Decision methods
  //@{
  /// Import - 0 if good
  /** Meaning of whereFrom:
     1 after initial solve by dualsimplex etc
     2 after preprocessing
     3 just before branchAndBound (so user can override)
     4 just after branchAndBound (before postprocessing)
     5 after postprocessing
     6 after a user called heuristic phase
     nonzero return code to stop
  */
  virtual int callBack(CbcModel * currentSolver, int whereFrom) {return 0;}
  //@}
  ///@name Constructors and destructors etc
  //@{
  /// Default Constructor
  CbcStopNow(); 
  
  /** Copy constructor .
   */  
  CbcStopNow(const CbcStopNow & rhs);
  
  /// Assignment operator 
  CbcStopNow & operator=(const CbcStopNow& rhs);

  /// Clone
  CbcStopNow * clone() const;
  
  /// Destructor 
  virtual ~CbcStopNow ();
  //@}
  
private:
  ///@name Private member data 
  //@{
 //@}
};
#endif
