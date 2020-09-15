// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcParam_H
#define CbcParam_H

#include "CbcModel.hpp"
#include "ClpParam.hpp"
#include "OsiSolverInterface.hpp"
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

{ // Begin used in Clp
  // Begin CbcGen
  CBC_PARAM_GENERALQUERY = -100,
  CBC_PARAM_FULLGENERALQUERY,
  CBC_PARAM_STR_ALLCOMMANDS,
  // End CbcGen
  // End used in Clp

  // Begin CbcGenCbc
  CBC_PARAM_DBL_INFEASIBILITYWEIGHT = 51,
  CBC_PARAM_DBL_CUTOFF,
  CBC_PARAM_DBL_INTEGERTOLERANCE,
  CBC_PARAM_DBL_INCREMENT,
  CBC_PARAM_DBL_ALLOWABLEGAP,
  CBC_PARAM_DBL_TIMELIMIT_BAB,
  CBC_PARAM_DBL_MAXSECONDSNIFS,
  CBC_PARAM_DBL_GAPRATIO,
  // End CbcGenCbc

  // Begin CbcGen
  CBC_PARAM_DBL_DJFIX = 81,
  CBC_PARAM_DBL_TIGHTENFACTOR,
  CBC_PARAM_DBL_FAKEINCREMENT,
  CBC_PARAM_DBL_FAKECUTOFF,
  CBC_PARAM_DBL_ARTIFICIALCOST,
  CBC_PARAM_DBL_DEXTRA3,
  CBC_PARAM_DBL_SMALLBAB,
  CBC_PARAM_DBL_DEXTRA4,
  CBC_PARAM_DBL_DEXTRA5,

  // Begin Used in Clp
  CBC_PARAM_INT_SOLVERLOGLEVEL = 101,
  CBC_PARAM_INT_OUTPUTFORMAT,
  CBC_PARAM_INT_PRINTOPTIONS,
  CBC_PARAM_INT_VERBOSE,

  CBC_PARAM_INT_THREADS = 151,
  // End used in Clp
  CBC_PARAM_INT_LPLOGLEVEL,
  CBC_PARAM_INT_CUTDEPTH,
  CBC_PARAM_INT_BKPIVOTINGSTRATEGY,
  CBC_PARAM_INT_BKMAXCALLS,
  CBC_PARAM_INT_BKCLQEXTMETHOD,
  CBC_PARAM_INT_ODDWEXTMETHOD,
  CBC_PARAM_INT_FPUMPITS,
  CBC_PARAM_INT_FPUMPTUNE,
  CBC_PARAM_INT_TESTOSI,
  CBC_PARAM_INT_EXTRA1,
  CBC_PARAM_INT_EXTRA2,
  CBC_PARAM_INT_EXTRA3,
  CBC_PARAM_INT_EXTRA4,
  CBC_PARAM_INT_DEPTHMINIBAB,
  CBC_PARAM_INT_CUTPASSINTREE,
  CBC_PARAM_INT_VUBTRY,
  CBC_PARAM_INT_EXPERIMENT,
  CBC_PARAM_INT_DIVEOPT,
  CBC_PARAM_INT_DIVEOPTSOLVES,
  CBC_PARAM_INT_STRATEGY,
  CBC_PARAM_INT_HEUROPTIONS,
  CBC_PARAM_INT_CUTLENGTH,
  CBC_PARAM_INT_FPUMPTUNE2,
  CBC_PARAM_INT_MAXSAVEDSOLS,
  CBC_PARAM_INT_RANDOMSEED,
  CBC_PARAM_INT_MULTIPLEROOTS,
  CBC_PARAM_INT_STRONG_STRATEGY,
  CBC_PARAM_INT_EXTRA_VARIABLES,
  CBC_PARAM_INT_MAXSLOWCUTS,
  CBC_PARAM_INT_MOREMOREMIPOPTIONS,
  CBC_PARAM_INT_PROCESSTUNE,
  CBC_PARAM_INT_USESOLUTION,
  // End CbcGen

  // Begin CbcGenCbc
  CBC_PARAM_INT_STRONGBRANCHING,
  CBC_PARAM_INT_MAXNODES,
  CBC_PARAM_INT_MAXNODESNOTIMPROVINGFS,
  CBC_PARAM_INT_NUMBERBEFORE,
  CBC_PARAM_INT_NUMBERANALYZE,
  CBC_PARAM_INT_MAXSOLS,
  CBC_PARAM_INT_MIPOPTIONS,
  CBC_PARAM_INT_MOREMIPOPTIONS,
  CBC_PARAM_INT_CUTPASS,
  // End CbcGenCbc

  // Begin CbcGenOsi
  CBC_PARAM_INT_MAXHOTITS,
  // End CbcGenOsi

  // Begin used in Cbc
  CBC_PARAM_STR_DIRECTION = 201,
  CBC_PARAM_STR_ERRORSALLOWED,
  CBC_PARAM_STR_INTPRINT,
  // End used in Cbc

  // Begin CbcGen
  CBC_PARAM_STR_NODESTRATEGY = 301,
  CBC_PARAM_STR_BRANCHSTRATEGY,
  CBC_PARAM_STR_CUTSSTRATEGY,
  CBC_PARAM_STR_HEURISTICSTRATEGY,
  CBC_PARAM_STR_USECGRAPH,
  CBC_PARAM_STR_CLIQUECUTS,
  CBC_PARAM_STR_ODDWHEELCUTS,
  CBC_PARAM_STR_CLQSTRENGTHENING,
  CBC_PARAM_STR_GOMORYCUTS,
  CBC_PARAM_STR_PROBINGCUTS,
  CBC_PARAM_STR_KNAPSACKCUTS,
  CBC_PARAM_STR_REDSPLITCUTS,
  CBC_PARAM_STR_ROUNDING,
  CBC_PARAM_STR_SOLVER,
  CBC_PARAM_STR_BRANCHPRIORITY,
  CBC_PARAM_STR_FLOWCUTS,
  CBC_PARAM_STR_MIXEDCUTS,
  CBC_PARAM_STR_TWOMIRCUTS,
  CBC_PARAM_STR_PREPROCESS,
  CBC_PARAM_STR_FPUMP,
  CBC_PARAM_STR_GREEDY,
  CBC_PARAM_STR_COMBINE,
  CBC_PARAM_STR_PROXIMITY,
  CBC_PARAM_STR_LOCALTREE,
  CBC_PARAM_STR_SOS,
  CBC_PARAM_STR_LANDPCUTS,
  CBC_PARAM_STR_RINS,
  CBC_PARAM_STR_RESIDCUTS,
  CBC_PARAM_STR_RENS,
  CBC_PARAM_STR_DIVINGS,
  CBC_PARAM_STR_DIVINGC,
  CBC_PARAM_STR_DIVINGF,
  CBC_PARAM_STR_DIVINGG,
  CBC_PARAM_STR_DIVINGL,
  CBC_PARAM_STR_DIVINGP,
  CBC_PARAM_STR_DIVINGV,
  CBC_PARAM_STR_DINS,
  CBC_PARAM_STR_PIVOTANDFIX,
  CBC_PARAM_STR_RANDROUND,
  CBC_PARAM_STR_NAIVE,
  CBC_PARAM_STR_ZEROHALFCUTS,
  CBC_PARAM_STR_CPX,
  CBC_PARAM_STR_CROSSOVER2,
  CBC_PARAM_STR_PIVOTANDCOMPLEMENT,
  CBC_PARAM_STR_VND,
  CBC_PARAM_STR_LAGOMORYCUTS,
  CBC_PARAM_STR_LATWOMIRCUTS,
  CBC_PARAM_STR_REDSPLIT2CUTS,
  CBC_PARAM_STR_GMICUTS,
  CBC_PARAM_STR_CUTOFF_CONSTRAINT,
  CBC_PARAM_STR_DW,
  CBC_PARAM_STR_ORBITAL,
  CBC_PARAM_STR_PREPROCNAMES,
  CBC_PARAM_STR_SOSPRIORITIZE,
  CBC_PARAM_STR_TIME_MODE,
  // End CbcGen

  // Begin used in Clp
  CBC_PARAM_ACTION_DIRECTORY = 401,
  CBC_PARAM_ACTION_DIRSAMPLE,
  CBC_PARAM_ACTION_DIRNETLIB,
  CBC_PARAM_ACTION_DIRMIPLIB,
  CBC_PARAM_ACTION_IMPORT,
  CBC_PARAM_ACTION_EXPORT,
  CBC_PARAM_ACTION_RESTORE,
  CBC_PARAM_ACTION_SAVE,
  CBC_PARAM_ACTION_MAXIMIZE,
  CBC_PARAM_ACTION_MINIMIZE,
  CBC_PARAM_ACTION_EXIT,
  CBC_PARAM_ACTION_STDIN,
  CBC_PARAM_ACTION_UNITTEST,
  CBC_PARAM_ACTION_SOLUTION,
  CBC_PARAM_ACTION_SAVESOL,
  CBC_PARAM_ACTION_HELP,
  CBC_PARAM_ACTION_REVERSE,
  CBC_PARAM_ACTION_STATISTICS,
  CBC_PARAM_ACTION_PRINTMASK,
  CBC_PARAM_ACTION_ENVIRONMENT,
  CBC_PARAM_ACTION_GMPL_SOLUTION,
  CBC_PARAM_ACTION_DUMMY,

  CBC_PARAM_ACTION_BAB = 501,
  // End Used in Clp
  CBC_PARAM_ACTION_DEBUG,
  CBC_PARAM_ACTION_MIPLIB,
  CBC_PARAM_ACTION_STRENGTHEN,
  CBC_PARAM_ACTION_PRIORITYIN,
  CBC_PARAM_ACTION_MIPSTART,
  CBC_PARAM_ACTION_USERCBC,
  CBC_PARAM_ACTION_DOHEURISTIC,
  CBC_PARAM_ACTION_NEXTBESTSOLUTION,
  CBC_PARAM_ACTION_CSVSTATISTICS,

  CBC_PARAM_NOTUSED_OSLSTUFF = 601,
  CBC_PARAM_NOTUSED_CBCSTUFF,

  CBC_PARAM_NOTUSED_INVALID = 1000
};

/// Very simple class for setting parameters

class CBCLIB_EXPORT CbcParam {

public:
  /**@name Constructor and destructor */
  //@{
  /// Constructors
  CbcParam();
  CbcParam(std::string name, std::string help, double lower, double upper,
           CbcParameterType type, int display = 2);
  CbcParam(std::string name, std::string help, int lower, int upper,
           CbcParameterType type, int display = 2);
  // Other strings will be added by insert
  CbcParam(std::string name, std::string help, std::string firstValue,
           CbcParameterType type, int whereUsed = 7, int display = 2);
  // Action
  CbcParam(std::string name, std::string help, CbcParameterType type,
           int whereUsed = 7, int display = 2);
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
  /// Returns name
  inline std::string name() const { return name_; };
  /// Returns short help
  inline std::string shortHelp() const { return shortHelp_; };
  /// Returns long help
  inline std::string longHelp() const { return longHelp_; }
  /// Returns set of valid strings
  inline const std::vector<std::string> &definedKeywords() const {
    return definedKeyWords_;
  }
  /// Returns the lower bound for a double-valued parameter
  inline double lowerDoubleValue() const { return lowerDoubleValue_; }
  /// Returns the upper bound for a double-valued parameter
  inline double upperDoubleValue() const { return upperDoubleValue_; }
  /// Returns the lower bound for an int-valued parameter
  inline int lowerIntValue() const { return lowerIntValue_; }
  /// Returns the upper bound for an int-valued parameter
  inline int upperIntValue() const { return upperIntValue_; }
  /// Gets a double parameter
  double doubleParameter(OsiSolverInterface *model) const;
  /// Gets a double parameter
  double doubleParameter(CbcModel &model) const;
  /// Sets a double parameter (nonzero code if error)
  int setDoubleParameter(OsiSolverInterface *model, double value);
  /// Sets a double parameter (nonzero code if error)
  int setDoubleParameter(CbcModel &model, double value);
  /// Sets double parameter and returns printable string and error code
  const char *setDoubleParameterWithMessage(OsiSolverInterface *model,
                                            double value, int &returnCode);
  /// Sets double parameter and returns printable string and error code
  const char *setDoubleParameterWithMessage(CbcModel &model, double value,
                                            int &returnCode);
  inline double doubleValue() const { return doubleValue_; }
  /// Sets double value
  void setDoubleValue(double value);
  /// Sets double value with message
  const char *setDoubleValueWithMessage(double value);
  /// Checks a double parameter (nonzero code if error)
  int checkDoubleParameter(double value) const;
  /// Gets a int parameter
  int intParameter(OsiSolverInterface *model) const;
  /// Gets a int parameter
  int intParameter(CbcModel &model) const;
  /// Sets a int parameter (nonzero code if error)
  int setIntParameter(CbcModel &model, int value);
  /// Sets a int parameter (nonzero code if error)
  int setIntParameter(OsiSolverInterface *model, int value);
  /// Sets int parameter and returns printable string and error code
  const char *setIntParameterWithMessage(OsiSolverInterface *model, int value,
                                         int &returnCode);
  /// Sets int parameter and returns printable string and error code
  const char *setIntParameterWithMessage(CbcModel &model, int value,
                                         int &returnCode);
  /// Gets int value
  inline int intValue() const { return intValue_; }
  /// Sets int value
  void setIntValue(int value);
  /// Sets int value with message
  const char *setIntValueWithMessage(int value);
  /// Gets string value
  inline std::string stringValue() const { return stringValue_; }
  /// Sets string value
  void setStringValue(std::string value);
  /// Returns name which could match
  std::string matchName() const;
  /// Returns length of name for ptinting
  int lengthMatchName() const;
  /// Returns 1 if matches minimum, 2 if matches less, 0 if not matched
  int matches(std::string input) const;
  /// Returns parameter option which matches (-1 if none)
  int parameterOption(std::string check) const;
  /// Prints parameter options
  void printOptions() const;
  /// Returns current parameter option
  inline std::string currentOption() const {
    return definedKeyWords_[currentKeyWord_];
  }
  /// Sets current parameter option
  void setCurrentOption(int value, bool printIt = false);
  /// Sets current parameter option using string
  void setCurrentOption(const std::string value);
  /// Sets current parameter option and returns printable string
  const char *setCurrentOptionWithMessage(int value);
  /// Sets current parameter option using string with message
  const char *setCurrentOptionWithMessage(const std::string value);
  /// Returns current parameter option position
  int currentOptionAsInteger() const;
  /** Returns current parameter option position
         but if fake keyword returns a fake value and sets
         fakeInteger to true value.  If not fake then fakeInteger is
     -COIN_INT_MAX
      */
  int currentOptionAsInteger(int &fakeInteger) const;
  /// type
  inline CbcParameterType type() const { return type_; }
  /// whether to display
  inline bool displayThis() const { return display_; }
  /// Set Long help
  inline void setLonghelp(const std::string help) { longHelp_ = help; }
  /// Print Long help
  void printLongHelp() const;
  /// Print action and string
  void printString() const;
  /// type for classification
  inline int indexNumber() const { return indexNumber_; }
  /** 7 if used everywhere,
         1 - used by clp
         2 - used by cbc
         4 - used by ampl
     */
  inline int whereUsed() const { return whereUsed_; }
  /// Gets value of fake keyword
  inline int fakeKeyWord() const { return fakeKeyWord_; }
  /// Sets value of fake keyword
  inline void setFakeKeyWord(int value, int fakeValue) {
    fakeKeyWord_ = value;
    fakeValue_ = fakeValue;
  }
  /// Sets value of fake keyword to current size of keywords
  void setFakeKeyWord(int fakeValue);

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
  std::vector<std::string> definedKeyWords_;
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
  /** 7 if used everywhere,
         1 - used by clp
         2 - used by cbc
         4 - used by ampl
     */
  int whereUsed_;
  /** If >=0 then integers allowed as a fake keyword
         So minusnnnn would got to -nnnn in currentKeyword_
         and plusnnnn would go to fakeKeyword_+nnnn
     */
  int fakeKeyWord_;
  /// Return this as main value if an integer
  int fakeValue_;
  //@}
  //@}
};

CBCLIB_EXPORT
void CbcReadPrintit(const char *input);
/// Simple read stuff
CLPLIB_EXPORT
std::string CbcReadNextField();
CBCLIB_EXPORT
std::string CbcReadGetCommand(int &whichArgument, int argc, const char *argv[]);
CBCLIB_EXPORT
std::string CbcReadGetString(int &whichArgument, int argc, const char *argv[]);
// valid 0 - okay, 1 bad, 2 not there
CBCLIB_EXPORT
int CbcReadGetIntField(int &whichArgument, int argc, const char *argv[],
                       int *valid);
CBCLIB_EXPORT
double CbcReadGetDoubleField(int &whichArgument, int argc, const char *argv[],
                             int *valid);
CLPLIB_EXPORT
void setCbcPrinting(bool yesNo);
CBCLIB_EXPORT
void setCbcReadCommand(FILE *f);

#define CBCMAXPARAMETERS 250
/*
  Subroutine to establish the cbc parameter array. See the description of
  class CbcParam for details.
*/
CLPLIB_EXPORT
void establishCbcParams(std::vector<CbcParam> &params);
// Given a parameter type - returns its number in list
CLPLIB_EXPORT
int whichCbcParam(const CbcParameterType &name,
                  const std::vector<CbcParam> &parameters);

#endif /* CbcParam_H */

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
