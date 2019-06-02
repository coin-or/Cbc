/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcParam_H
#define CbcParam_H

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
class ClpSimplex;
/*! \brief Parameter codes

  Parameter type ranges are allocated as follows
  <ul>
    <li>   1 -- 100	double parameters
    <li> 101 -- 200	integer parameters
    <li> 201 -- 250	string parameters
    <li> 251 -- 300	cuts etc(string but broken out for clarity)
    <li> 301 -- 400	`actions'
  </ul>

  `Actions' do not necessarily invoke an immediate action; it's just that they
  don't fit neatly into the parameters array.

  This coding scheme is in flux.
  CBC_PARAM_STR_NODESTRATEGY, 
  CBC_PARAM_STR_BRANCHSTRATEGY, 
  CBC_PARAM_NOTUSED_ADDCUTSSTRATEGY,
  CLP_PARAM_ACTION_CLEARCUTS, 
  CBC_PARAM_NOTUSED_OSLSTUFF, 
  CBC_PARAM_NOTUSED_CBCSTUFF are not used at present (03.10.24).
*/

enum CbcParameterType

{ CBC_PARAM_GENERALQUERY = -100,
  CBC_PARAM_FULLGENERALQUERY,

  CLP_PARAM_DBL_PRIMALTOLERANCE = 1,
  CLP_PARAM_DBL_DUALTOLERANCE,
  CBC_PARAM_DBL_CUTOFF,
  CLP_PARAM_DBL_TIMELIMIT,
  CLP_PARAM_DBL_DUALBOUND,
  CLP_PARAM_DBL_PRIMALWEIGHT,
  CLP_PARAM_DBL_OBJSCALE,
  CLP_PARAM_DBL_RHSSCALE,

  CBC_PARAM_DBL_INFEASIBILITYWEIGHT = 51,
  CBC_PARAM_DBL_INTEGERTOLERANCE,
  CBC_PARAM_DBL_INCREMENT,
  CBC_PARAM_DBL_ALLOWABLEGAP,

  CBC_PARAM_DBL_DJFIX = 81,
  CBC_PARAM_DBL_GAPRATIO,
  CBC_PARAM_DBL_TIGHTENFACTOR,

  CLP_PARAM_INT_LOGLEVEL = 101,
  CLP_PARAM_INT_SOLVERLOGLEVEL,
  CBC_PARAM_INT_MAXNODES,
  CBC_PARAM_INT_STRONGBRANCHING,
  CLP_PARAM_INT_MAXFACTOR,
  CLP_PARAM_INT_PERTVALUE,
  CLP_PARAM_INT_MAXITERATION,
  CLP_PARAM_INT_PRESOLVEPASS,
  CLP_PARAM_INT_IDIOT,
  CLP_PARAM_INT_SPRINT,
  CLP_PARAM_INT_OUTPUTFORMAT,
  CLP_PARAM_INT_SLPVALUE,
  CLP_PARAM_INT_PRESOLVEOPTIONS,
  CLP_PARAM_INT_PRINTOPTIONS,
  CLP_PARAM_INT_SPECIALOPTIONS,

  CLP_PARAM_STR_DIRECTION = 201,
  CLP_PARAM_STR_DUALPIVOT,
  CLP_PARAM_STR_SCALING,
  CLP_PARAM_STR_ERRORSALLOWED,
  CLP_PARAM_STR_KEEPNAMES,
  CLP_PARAM_STR_SPARSEFACTOR,
  CLP_PARAM_STR_PRIMALPIVOT,
  CLP_PARAM_STR_PRESOLVE,
  CLP_PARAM_STR_CRASH,
  CLP_PARAM_STR_BIASLU,
  CLP_PARAM_STR_PERTURBATION,
  CLP_PARAM_STR_MESSAGES,
  CLP_PARAM_STR_AUTOSCALE,
  CLP_PARAM_STR_CHOLESKY,
  CLP_PARAM_STR_KKT,
  CLP_PARAM_STR_BARRIERSCALE,
  CLP_PARAM_STR_GAMMA,
  CLP_PARAM_STR_CROSSOVER,
  CLP_PARAM_STR_PFI,
  CLP_PARAM_NOTUSED_ALGORITHM,

  CBC_PARAM_STR_NODESTRATEGY = 251,
  CBC_PARAM_STR_BRANCHSTRATEGY,
  CBC_PARAM_NOTUSED_ADDCUTSSTRATEGY,
  CBC_PARAM_STR_GOMORYCUTS,
  CBC_PARAM_STR_PROBINGCUTS,
  CBC_PARAM_STR_KNAPSACKCUTS,
  CBC_PARAM_NOTUSED_ODDHOLECUTS,
  CBC_PARAM_STR_ROUNDING,
  CBC_PARAM_STR_SOLVER,
  CBC_PARAM_STR_CLIQUECUTS,
  CBC_PARAM_STR_COSTSTRATEGY,
  CBC_PARAM_STR_FLOWCUTS,
  CBC_PARAM_STR_MIXEDCUTS,
  CBC_PARAM_STR_TWOMIRCUTS,
  CBC_PARAM_STR_PREPROCESS,

  CLP_PARAM_ACTION_DIRECTORY = 301,
  CLP_PARAM_ACTION_IMPORT,
  CLP_PARAM_ACTION_EXPORT,
  CLP_PARAM_ACTION_RESTORE,
  CLP_PARAM_ACTION_SAVE,
  CLP_PARAM_ACTION_DUALSIMPLEX,
  CLP_PARAM_ACTION_PRIMALSIMPLEX,
  CLP_PARAM_ACTION_MAXIMIZE,
  CLP_PARAM_ACTION_MINIMIZE,
  CLP_PARAM_ACTION_EXIT,
  CLP_PARAM_ACTION_STDIN,
  CLP_PARAM_ACTION_UNITTEST,
  CLP_PARAM_ACTION_NETLIB_DUAL,
  CLP_PARAM_ACTION_NETLIB_PRIMAL,
  CLP_PARAM_ACTION_SOLUTION,
  CLP_PARAM_ACTION_TIGHTEN,
  CLP_PARAM_ACTION_FAKEBOUND,
  CLP_PARAM_ACTION_HELP,
  CLP_PARAM_ACTION_PLUSMINUS,
  CLP_PARAM_ACTION_NETWORK,
  CLP_PARAM_ACTION_ALLSLACK,
  CLP_PARAM_ACTION_REVERSE,
  CLP_PARAM_ACTION_BARRIER,
  CLP_PARAM_ACTION_NETLIB_BARRIER,
  CLP_PARAM_ACTION_REALLY_SCALE,
  CLP_PARAM_ACTION_BASISIN,
  CLP_PARAM_ACTION_BASISOUT,
  CLP_PARAM_ACTION_SOLVECONTINUOUS,
  CBC_PARAM_ACTION_BAB,
  CBC_PARAM_ACTION_MIPLIB,
  CLP_PARAM_ACTION_CLEARCUTS,
  CLP_VERSION_NOTUSED_PRINTVERSION,

  CBC_PARAM_NOTUSED_OSLSTUFF = 401,
  CBC_PARAM_NOTUSED_CBCSTUFF,

  CBC_PARAM_NOTUSED_INVALID = 1000
};

/// Very simple class for setting parameters

class CbcParam {

public:
  /**@name Constructor and destructor */
  //@{
  /// Constructors
  CbcParam();
  CbcParam(std::string name, std::string help,
    double lower, double upper, CbcParameterType type, bool display = true);
  CbcParam(std::string name, std::string help,
    int lower, int upper, CbcParameterType type, bool display = true);
  // Other strings will be added by insert
  CbcParam(std::string name, std::string help, std::string firstValue,
    CbcParameterType type, int defaultIndex = 0, bool display = true);
  // Action
  CbcParam(std::string name, std::string help,
    CbcParameterType type, int indexNumber = -1, bool display = true);
  /// Copy constructor.
  CbcParam(const CbcParam &);
  /// Assignment operator. This copies the data
  CbcParam &operator=(const CbcParam &rhs);
  /// Destructor
  ~CbcParam();
  //@}

  /**@name stuff */
  //@{
  /// Insert string (only valid for keywords)
  void append(std::string keyWord);
  /// Adds one help line
  void addHelp(std::string keyWord);
  /// Returns name
  inline std::string name() const
  {
    return name_;
  };
  /// Returns short help
  inline std::string shortHelp() const
  {
    return shortHelp_;
  };
  /// Sets a double parameter (nonzero code if error)
  int setDoubleParameter(CbcModel &model, double value) const;
  /// Gets a double parameter
  double doubleParameter(CbcModel &model) const;
  /// Sets a int parameter (nonzero code if error)
  int setIntParameter(CbcModel &model, int value) const;
  /// Gets a int parameter
  int intParameter(CbcModel &model) const;
  /// Sets a double parameter (nonzero code if error)
  int setDoubleParameter(ClpSimplex *model, double value) const;
  /// Gets a double parameter
  double doubleParameter(ClpSimplex *model) const;
  /// Sets a int parameter (nonzero code if error)
  int setIntParameter(ClpSimplex *model, int value) const;
  /// Gets a int parameter
  int intParameter(ClpSimplex *model) const;
  /// Sets a double parameter (nonzero code if error)
  int setDoubleParameter(OsiSolverInterface *model, double value) const;
  /// Gets a double parameter
  double doubleParameter(OsiSolverInterface *model) const;
  /// Sets a int parameter (nonzero code if error)
  int setIntParameter(OsiSolverInterface *model, int value) const;
  /// Gets a int parameter
  int intParameter(OsiSolverInterface *model) const;
  /// Checks a double parameter (nonzero code if error)
  int checkDoubleParameter(double value) const;
  /// Returns name which could match
  std::string matchName() const;
  /// Returns parameter option which matches (-1 if none)
  int parameterOption(std::string check) const;
  /// Prints parameter options
  void printOptions() const;
  /// Returns current parameter option
  inline std::string currentOption() const
  {
    return definedKeyWords_[currentKeyWord_];
  }
  /// Sets current parameter option
  inline void setCurrentOption(int value)
  {
    currentKeyWord_ = value;
  }
  /// Sets int value
  inline void setIntValue(int value)
  {
    intValue_ = value;
  }
  inline int intValue() const
  {
    return intValue_;
  }
  /// Sets double value
  inline void setDoubleValue(double value)
  {
    doubleValue_ = value;
  }
  inline double doubleValue() const
  {
    return doubleValue_;
  }
  /// Sets string value
  inline void setStringValue(std::string value)
  {
    stringValue_ = value;
  }
  inline std::string stringValue() const
  {
    return stringValue_;
  }
  /// Returns 1 if matches minimum, 2 if matches less, 0 if not matched
  int matches(std::string input) const;
  /// type
  inline CbcParameterType type() const
  {
    return type_;
  }
  /// whether to display
  inline bool displayThis() const
  {
    return display_;
  }
  /// Set Long help
  inline void setLonghelp(const std::string help)
  {
    longHelp_ = help;
  }
  /// Print Long help
  void printLongHelp() const;
  /// Print action and string
  void printString() const;
  /// type for classification
  inline int indexNumber() const
  {
    return indexNumber_;
  }

private:
  /// gutsOfConstructor
  void gutsOfConstructor();
  //@}
  ////////////////// data //////////////////
private:
  /**@name data
     We might as well throw all type data in - could derive?
    */
  //@{
  // Type see CbcParameterType
  CbcParameterType type_;
  /// If double == okay
  double lowerDoubleValue_;
  double upperDoubleValue_;
  /// If int == okay
  int lowerIntValue_;
  int upperIntValue_;
  // Length of name
  unsigned int lengthName_;
  // Minimum match
  unsigned int lengthMatch_;
  /// set of valid strings
  std::vector< std::string > definedKeyWords_;
  /// Name
  std::string name_;
  /// Short help
  std::string shortHelp_;
  /// Long help
  std::string longHelp_;
  /// Action
  CbcParameterType action_;
  /// Current keyWord (if a keyword parameter)
  int currentKeyWord_;
  /// Display on ?
  bool display_;
  /// Integer parameter - current value
  int intValue_;
  /// Double parameter - current value
  double doubleValue_;
  /// String parameter - current value
  std::string stringValue_;
  /// index number to use for display purposes
  int indexNumber_;
  //@}
};
#endif /* CbcParam_H */

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
