#ifndef CoinParam_H
#define CoinParam_H

/*
  Copyright (C) 2002, International Business Machines
  Corporation and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).
*/

/*! \file CoinParam.hpp
    \brief Declaration of a class for command line parameters.
*/

#include <vector>
#include <deque>
#include <string>
#include <cstdio>
#include <map>

#include "CoinUtilsConfig.h"
#include "CoinFinite.hpp"
#include "CoinHelperFunctions.hpp"

/*! \class CoinParam
    \brief A base class for `keyword value' command line parameters.

  The underlying paradigm is that a parameter specifies an action to be
  performed on a target object. The base class provides two function
  pointers, a `push' function and a `pull' function.  By convention, a push
  function will set some value in the target object or perform some action
  using the target object.  A `pull' function will retrieve some value from
  the target object.  This is only a convention, however; CoinParam and
  associated utilities make no use of these functions and have no hardcoded
  notion of how they should be used.

  The action to be performed, and the target object, will be specific to a
  particular application. It is expected that users will derive
  application-specific parameter classes from this base class. A derived
  class will typically add fields and methods to set/get a code for the
  action to be performed (often, an enum class) and the target object (often,
  a pointer or reference).

  Facilities provided by the base class and associated utility routines
  include:
  <ul>
    <li> Support for common parameter types with numeric, string, or
	 keyword values.
    <li> Support for short and long help messages.
    <li> Pointers to `push' and `pull' functions as described above.
    <li> Command line parsing and keyword matching.
  </ul>
  All utility routines are declared in the #CoinParamUtils namespace.

  The base class recognises five types of parameters: actions (which require
  no value); numeric parameters with integer or real (double) values; keyword
  parameters, where the value is one of a defined set of value-keywords;
  and string parameters (where the value is a string).
  The base class supports the definition of a valid range, a default value,
  and short and long help messages for a parameter.

  As defined by the #CoinParamFunc typedef, push and pull functions
  should take a single parameter, a pointer to a CoinParam. Typically this
  object will actually be a derived class as described above, and the
  implementation function will have access to all capabilities of CoinParam and
  of the derived class.

  When specified as command line parameters, the expected syntax is `-keyword
  value' or `-keyword=value'. You can also use the Gnu double-dash style,
  `--keyword'. Spaces around the `=' will \e not work.

  The keyword (name) for a parameter can be defined with an `!' to mark the
  minimal match point. For example, allow!ableGap will be considered matched
  by the strings `allow', `allowa', `allowab', \e etc. Similarly, the
  value-keyword strings for keyword parameters can be defined with `!' to
  mark the minimal match point.  Matching of keywords and value-keywords is
  \e not case sensitive.
*/

class COINUTILSLIB_EXPORT CoinParam
{
 
public:
  /*! \name Subtypes */
  //@{

  /*! \brief Enumeration for the types of parameters supported by CoinParam

    CoinParam provides support for several types of parameters:
    <ul>
      <li> Action parameters, which require no value.
      <li> Integer and double numeric parameters, with upper and lower bounds.
      <li> String parameters that take an arbitrary string value.
      <li> Keyword parameters that take a defined set of string (value-keyword)
	   values. Value-keywords are associated with integers in the order in
	   which they are added, starting from zero.
    </ul>
  */
  enum CoinParamType {
     paramInvalid = 0,
     paramAct,
     paramInt,
     paramDbl,
     paramStr,
     paramDir,
     paramFile,
     paramKwd
  };

  enum CoinDisplayPriority {
     displayPriorityNone = 0,
     displayPriorityLow,
     displayPriorityHigh
  };

  enum ParamPushMode {
     pushDefault = 0,
     pushOn = 1,
     pushOff = 2
  };

  /*! \brief Type declaration for push and pull functions.

    By convention, a return code of  0 indicates execution without error, >0
    indicates nonfatal error, and <0 indicates fatal error. This is only
    convention, however; the base class makes no use of the push and pull
    functions and has no hardcoded interpretation of the return code.
  */
  typedef int (*CoinParamFunc)(CoinParam &param);

  //@}

  /*! \name Constructors and Destructors

  Be careful how you specify parameters for the constructors! Some compilers
  are entirely too willing to convert almost anything to bool.
*/
  //@{

  /*! \brief Default constructor */

  CoinParam();

  /*! \brief Constructor for a parameter with a double value
  
    The default value is 0.0. Be careful to clearly indicate that \p lower and
    \p upper are real (double) values to distinguish this constructor from the
    constructor for an integer parameter.
  */
  CoinParam(std::string name, std::string help,
            double lower = -COIN_DBL_MAX, double upper = COIN_DBL_MAX,
            std::string longHelp = "",
            CoinDisplayPriority displayPriority = displayPriorityHigh);

  /*! \brief Constructor for a parameter with an integer value
  
    The default value is 0.
  */
  CoinParam(std::string name, std::string help,
            int lower = -COIN_INT_MAX, int upper = COIN_INT_MAX,
            std::string longHelp = "",
            CoinDisplayPriority displayPriority = displayPriorityHigh);

  /*! \brief Constructor for any parameter taking a string (or no) value 
    (all others). Type is not optional to resolve ambiguity.

    The default value is "" for all such parameter types
  */
  CoinParam(std::string name, CoinParamType type, std::string help,
            std::string longHelp = "",
            CoinDisplayPriority displayPriority = displayPriorityHigh);

  /*! \brief Copy constructor */
  CoinParam(const CoinParam &orig);

  /*! \brief Clone */
  virtual CoinParam *clone();

  /*! \brief Assignment */
  CoinParam &operator=(const CoinParam &rhs);

  /*! \brief  Destructor */
  virtual ~CoinParam();

  //@}

  /*! \name Methods to initialize parameters already constructed */
  //@{

  /*! \brief These mirror the constructors with the addition of long help */
  void setup(std::string name, std::string help,
             double lower = -COIN_DBL_MAX, double upper = COIN_DBL_MAX,
             std::string longHelp = "",
             CoinDisplayPriority display = displayPriorityHigh);

  void setup(std::string name, std::string help,
             int lower = -COIN_INT_MAX, int upper = COIN_INT_MAX,
             std::string longHelp = "",
             CoinDisplayPriority display = displayPriorityHigh);

  void setup(std::string name, std::string help, std::string longHelp = "",
             CoinDisplayPriority display = displayPriorityHigh);
  
  //@}
  /*! \name Methods to query and manipulate the value(s) of a parameter */
  //@{

  /*! Set the value of a parameter of any type */
  int setVal(std::string value, std::string *message = NULL,
             ParamPushMode pMode = pushDefault);
  int setVal(double value, std::string *message = NULL,
             ParamPushMode pMode = pushDefault);
  int setVal(int value, std::string *message = NULL,
             ParamPushMode pMode = pushDefault);
  
  /*! Set the default value of a parameter of any type */
  int setDefault(std::string value, std::string *message = NULL);
  int setDefault(double value, std::string *message = NULL);
  int setDefault(int value, std::string *message = NULL);
  
  /*! Restore the default value of a parameter of any type */
  int restoreDefault();

  /*! Get the value of a parameter of any type */ 
  int getVal(std::string &value) const;
  int getVal(double &value) const;
  int getVal(int &value) const;

  /*! Try to read the value of a parameter from the input queue*/ 
  int readValue(std::deque<std::string> &inputQueue, std::string &value,
               std::string *message = NULL);
  int readValue(std::deque<std::string> &inputQueue, int &value,
                std::string *message = NULL);
  int readValue(std::deque<std::string> &inputQueue, double &value,
                std::string *message = NULL);

  /*! \brief Add an additional value-keyword to a keyword parameter */
  void appendKwd(std::string kwd, int index);

  /*! \brief Add an additional value-keyword to a keyword parameter,
    generating sequential mode values */
  void appendKwd(std::string kwd);

  /*! \brief Return the integer associated with the specified value-keyword
  
    Returns -1 if no value-keywords match the specified string.
  */
  int kwdToMode(std::string kwd) const;

  /*! \brief Return the value-keyword that is the current value of the
	     keyword parameter
  */
  std::string kwdVal() const;

  /*! \brief Return the integer mode value of the current keyword parameter
  */
  int modeVal() const;

  /*! \brief Return the string for an integer mode of keyword parameter
  */
  std::string  modeString(int value) const;

  /*! \brief Set the value of the keyword parameter using the integer
	     associated with a value-keyword.
  
    If \p printIt is true, the corresponding value-keyword string will be
    echoed to std::cout.
  */
  int setModeVal(int value, std::string *message = NULL,
                 ParamPushMode pMode = pushDefault);

  /*! \brief Set the default value of the keyword parameter using the integer
	     associated with a value-keyword.
  
    If \p printIt is true, the corresponding value-keyword string will be
    echoed to std::cout.
  */
  int setModeValDefault(int value, std::string *message = NULL);

  /*! \brief Set the value of the keyword parameter using a value-keyword
	     string.
  
    The given string will be tested against the set of value-keywords for
    the parameter using the shortest match rules.
  */
  int setKwdVal(const std::string value, std::string *message = NULL,
                ParamPushMode pMode = pushDefault);

  /*! \brief Set the default value of the keyword parameter using a value-keyword
	     string.
  
    The given string will be tested against the set of value-keywords for
    the parameter using the shortest match rules.
  */
  int setKwdValDefault(const std::string value, std::string *message = NULL);

  /*! \brief Prints the set of value-keywords defined for this keyword
	     parameter
  */
  std::string printKwds() const;

  /// Set of valid value-keywords for a keyword parameter
  const std::map<std::string, int>& definedKwds() const
  {
    return definedKwds_;
  }

  /// Valid value-keywords for a keyword parameter sorted by integer mode
  std::vector<std::string> definedKwdsSorted() const;

  /*! \brief Set the value of a string parameter */
  int setStrVal(std::string value, std::string *message = NULL,
                ParamPushMode pMode = pushDefault);

  /*! \brief Set the default value of a string parameter */
  int setStrValDefault(std::string value, std::string *message = NULL);

  /*! \brief Get the value of a string parameter */
  std::string strVal() const;

  /*! \brief Set the value of a directory parameter */
  int setDirName(std::string value, std::string *message = NULL,
                 ParamPushMode pMode = pushDefault);

  /*! \brief Set the default value of a directory parameter */
  int setDirNameDefault(std::string value, std::string *message = NULL);

  /*! \brief Get the value of a directory parameter */
  std::string dirName() const;

  /*! \brief Set the value of a file name parameter */
  int setFileName(std::string value, std::string *message = NULL,
                  ParamPushMode pMode = pushDefault);

  /*! \brief Set the default value of a file name parameter */
  int setFileNameDefault(std::string value, std::string *message = NULL);

  /*! \brief Get the value of a file name parameter */
  std::string fileName() const;

  /*! \brief Set the value of a double parameter */
  int setDblVal(double value, std::string *message = NULL,
                ParamPushMode pMode = pushDefault);

  /*! \brief Set the default value of a double parameter */
  int setDblValDefault(double value, std::string *message = NULL);

  /*! \brief Get the value of a double parameter */
  double dblVal() const;

  /*! \brief Set the lower value of an double parameter */
  void setLowerDblVal(double value);

  /*! \brief Get the lower value of a double parameter */
  double lowerDblVal() const;

  /*! \brief Set the upper value of a double parameter */
  void setUpperDblVal(double value);

  /*! \brief Get the upper value of a double parameter */
  double upperDblVal() const;

  /*! \brief Set the value of a integer parameter */
  int setIntVal(int value, std::string *message = NULL,
                ParamPushMode pMode = pushDefault);

  /*! \brief Set the value of a integer parameter */
  int setIntValDefault(int value, std::string *message = NULL);

  /*! \brief Get the value of a integer parameter */
  int intVal() const;

  /*! \brief Set the lower value of an integer parameter */
  void setLowerIntVal(int value);

  /*! \brief Get the lower value of a integer parameter */
  int lowerIntVal() const;

  /*! \brief Set the upper value of a integer parameter */
  void setUpperIntVal(int value);

  /*! \brief Get the upper value of a integer parameter */
  int upperIntVal() const;

  /*! \brief Print possible options for a keyword parameter */
  void printOptions(std::string *message = NULL);

  /*! \brief Add a short help string to a parameter */
  inline void setShortHelp(const std::string help) { shortHelp_ = help; }

  /*! \brief Retrieve the short help string */
  inline std::string shortHelp() const { return (shortHelp_); }

  /*! \brief Add a long help message to a parameter
  
    See printLongHelp() for a description of how messages are broken into
    lines.
  */
  inline void setLongHelp(const std::string help) { longHelp_ = help; }

  /*! \brief Retrieve the long help message */
  inline std::string longHelp() const { return (longHelp_); }

  /*! \brief Set the topic (semantic category) for this parameter.

    Topics group parameters by function (e.g., "Cuts", "Heuristics",
    "Stopping", "Output") independently of their data type.  Useful for
    generating documentation grouped by subject and for topic-filtered
    help listings.
  */
  inline void setTopic(const std::string &topic) { topic_ = topic; }

  /*! \brief Retrieve the topic string (empty if not set) */
  inline const std::string &topic() const { return (topic_); }

  /*! \brief  Print long help

    Prints the long help string, plus the valid range and/or keywords if
    appropriate. The routine makes a best effort to break the message into
    lines appropriate for an 80-character line. Explicit line breaks in the
    message will be observed. The short help string will be used if
    long help is not available.
  */
  std::string printLongHelp() const;

  //@}

  /*! \name Methods to query and manipulate a parameter object */
  //@{

  /*! \brief Return the type of the parameter */
  inline CoinParamType type() const { return (type_); }

  /*! \brief Set the type of the parameter */
  inline void setType(CoinParamType type) { type_ = type; }

  /*! \brief Return the parameter keyword (name) string */
  inline std::string name() const { return (name_); }

  /*! \brief Set the parameter keyword (name) string */
  inline void setName(std::string name)
  {
    name_ = name;
    processName();
  }

  /*! \brief Check if the specified string matches the parameter keyword (name)
	     string
  
    Returns 1 if the string matches and meets the minimum match length,
    2 if the string matches but doesn't meet the minimum match length,
    and 0 if the string doesn't match. Matches are \e not case-sensitive.
  */
  int matches(std::string input) const;

  /*! \brief Return the parameter keyword (name) string formatted to show
	     the minimum match length
  
    For example, if the parameter name was defined as allow!ableGap, the
    string returned by matchName would be allow(ableGap).
  */
  std::string matchName() const;

  /*! \brief Return the length of the formatted paramter name for printing. */
  int lengthMatchName() const;

  /*! \brief Set visibility of parameter

    Intended to control whether the parameter is shown when a list of
    parameters is processed. Used by CoinParamUtils::printHelp when printing
    help messages for a list of parameters.
  */
  inline void setDisplayPriority(CoinDisplayPriority display) {
     display_ = display;
  }

  /*! \brief Get visibility of parameter */
  inline CoinDisplayPriority getDisplayPriority() const { return (display_); }

  /*! \brief Get push function */
  inline CoinParamFunc pushFunc() { return (pushFunc_); }

  /*! \brief Set push function */
  inline void setPushFunc(CoinParamFunc func) { pushFunc_ = func; }

  /*! \brief Get pull function */
  inline CoinParamFunc pullFunc() { return (pullFunc_); }

  /*! \brief Set pull function */
  inline void setPullFunc(CoinParamFunc func) { pullFunc_ = func; }

  //@}

protected:
  /*! \name Protected methods */
  //@{

  /*! Process a name for efficient matching */
  void processName();

  //@}

  /*! \name Protected parameter data */
  //@{
  /// Parameter type (see #CoinParamType)
  CoinParamType type_;

  /// Parameter name
  std::string name_;

  /// Length of parameter name
  size_t lengthName_;

  /*! \brief  Minimum length required to declare a match for the parameter
	      name.
  */
  size_t lengthMatch_;

  /// Lower bound on value for a double parameter
  double lowerDblValue_;

  /// Upper bound on value for a double parameter
  double upperDblValue_;

  /// Double parameter - current value
  double dblValue_;

  /// Double parameter - current value
  double dblDefaultValue_;

  /// Lower bound on value for an integer parameter
  int lowerIntValue_;

  /// Upper bound on value for an integer parameter
  int upperIntValue_;

  /// Integer parameter - current value
  int intValue_;

  /// Integer parameter - default value
  int intDefaultValue_;

  /// String parameter - current value
  std::string strValue_;

  /// String parameter - default value
  std::string strDefaultValue_;

  /** Set of valid value-keywords for a keyword parameter
      (can also be used for keyword/bit for integer values */
  std::map< std::string, int > definedKwds_;

  /*! \brief Current integer value (mode) for a keyword parameter 
  */
  int currentMode_;

  /*! \brief Default integer value (mode) for a keyword parameter 
  */
  int defaultMode_;

  /*! \brief Current string value for a keyword parameter 
  */
  std::string currentKwd_;

  /*! \brief Default string value for a keyword parameter 
  */
  std::string defaultKwd_;

  /// Push function
  CoinParamFunc pushFunc_;

  /// Pull function
  CoinParamFunc pullFunc_;

  /// Short help
  std::string shortHelp_;

  /// Long help
  std::string longHelp_;

  /// Display when processing lists of parameters?
  CoinDisplayPriority display_;

  /// Topic (semantic category, e.g. "Cuts", "Heuristics", "Output")
  std::string topic_;
  //@}
};

/*! \relatesalso CoinParam
    \brief A type for a parameter vector.
*/
typedef std::vector< CoinParam * > CoinParamVec;

/*! \relatesalso CoinParam
    \brief A stream output function for a CoinParam object.
*/
COINUTILSLIB_EXPORT
std::ostream &operator<<(std::ostream &s, const CoinParam &param);

/*
  Bring in the utility functions for parameter handling (CbcParamUtils).
*/

/*! \brief Utility functions for processing CoinParam parameters.

  The functions in CoinParamUtils support command line or interactive
  parameter processing and a help facility. Consult the `Related Functions'
  section of the CoinParam class documentation for individual function
  documentation.
*/
namespace CoinParamUtils {

/*! \relatesalso CoinParam
      \brief Scan a parameter vector for parameters whose keyword (name) string
	     matches \p name using minimal match rules.
      
       \p matchNdx is set to the index of the last parameter that meets the
       minimal match criteria (but note there should be at most one matching
       parameter if the parameter vector is properly configured). \p shortCnt
       is set to the number of short matches (should be zero for a properly
       configured parameter vector if a minimal match is found). The return
       value is the number of matches satisfying the minimal match requirement
       (should be 0 or 1 in a properly configured vector).
  */
COINUTILSLIB_EXPORT
int matchParam(const CoinParamVec &paramVec, std::string name,
               int &matchNdx, int &shortCnt);

/*! \relatesalso CoinParam
      \brief Look up the command keyword (name) in the parameter vector.
      	     Print help if requested.

    In the most straightforward use, \p name is a string without `?', and the
    value returned is the index in \p paramVec of the single parameter that
    matched \p name. One or more '?' characters at the end of \p name is a
    query for information. The routine prints short (one '?') or long (more
    than one '?') help messages for a query.  Help is also printed in the case
    where the name is ambiguous (some of the matches did not meet the minimal
    match length requirement).

    Note that multiple matches meeting the minimal match requirement is a
    configuration error. The mimimal match length for the parameters
    involved is too short.

    If provided as parameters, on return
    <ul>
      <li> \p matchCnt will be set to the number of matches meeting the
	   minimal match requirement
      <li> \p shortCnt will be set to the number of matches that did not
	   meet the miminal match requirement
      <li> \p queryCnt will be set to the number of '?' characters at the
	   end of the name
    </ul>

    The return values are:
    <ul>
      <li> >0: index in \p paramVec of the single unique match for \p name
      <li> -1: a query was detected (one or more '?' characters at the end
	       of \p name
      <li> -2: one or more short matches, not a query
      <li> -3: no matches, not a query
      <li> -4: multiple matches meeting the minimal match requirement
	       (configuration error)
    </ul>
  */
COINUTILSLIB_EXPORT
int lookupParam(std::string name, CoinParamVec &paramVec,
                int *matchCnt = 0, int *shortCnt = 0, int *queryCnt = 0);

/*! \relatesalso CoinParam
      \brief Utility to print a long message as filled lines of text

      The routine makes a best effort to break lines without exceeding the
      standard 80 character line length. Explicit newlines in \p msg will
      be obeyed.
  */
COINUTILSLIB_EXPORT
void printIt(const char *msg);

/*! \relatesalso CoinParam
      \brief Utility routine to print help given a short match or explicit
	     request for help.

      The two really are related, in that a query (a string that ends with
      one or more `?' characters) will often result in a short match. The
      routine expects that \p name matches a single parameter, and does not
      look for multiple matches.
      
      If called with \p matchNdx < 0, the routine will look up \p name in \p
      paramVec and print the full name from the parameter. If called with \p
      matchNdx > 0, it just prints the name from the specified parameter.  If
      the name is a query, short (one '?') or long (more than one '?') help
      is printed.

  */
COINUTILSLIB_EXPORT
void shortOrHelpOne(CoinParamVec &paramVec, int matchNdx, std::string name, int numQuery);

/*! \relatesalso CoinParam
      \brief Utility routine to print help given multiple matches.

      If the name is not a query, or asks for short help (\e i.e., contains
      zero or one '?' characters), the list of matching names is printed. If
      the name asks for long help (contains two or more '?' characters),
      short help is printed for each matching name.
  */
COINUTILSLIB_EXPORT
void shortOrHelpMany(CoinParamVec &paramVec, std::string name, int numQuery);

/*! \relatesalso CoinParam
      \brief Print a generic `how to use the command interface' help message.

    The message is hard coded to match the behaviour of the parsing utilities.
  */
COINUTILSLIB_EXPORT
void printGenericHelp();

/*! \relatesalso CoinParam
      \brief Utility routine to print help messages for one or more
	     parameters.
    
    Intended as a utility to implement explicit `help' commands. Help will be
    printed for all parameters in \p paramVec from \p firstParam to \p
    lastParam, inclusive. If \p shortHelp is true, short help messages will
    be printed. If \p longHelp is true, long help messages are printed. \p
    shortHelp overrules \p longHelp. If neither is true, only command
    keywords are printed. \p prefix is printed before each line; it's an
    imperfect attempt at indentation.
  */
COINUTILSLIB_EXPORT
void printHelp(CoinParamVec &paramVec, int firstParam, int lastParam,
               std::string prefix,
               bool shortHelp, bool longHelp, bool hidden);

COINUTILSLIB_EXPORT
std::string printString(const std::string &input, int maxWidth = 65);

COINUTILSLIB_EXPORT
std::string printString(const char *input, int maxWidth = 65);

COINUTILSLIB_EXPORT
void readFromStream(std::deque<std::string> &inputQueue,
                    std::istream &inputStream);

COINUTILSLIB_EXPORT
std::string getNextField(std::deque<std::string> &inputQueue,
                         bool interactiveMode = false,
                         std::string prompt="");

COINUTILSLIB_EXPORT
void processFile(std::string &fileName, std::string dirName,
                 bool *fileExists = NULL);

COINUTILSLIB_EXPORT
void formInputQueue(std::deque<std::string> &inputQueue,
                    std::string commandName,
                    int argc, char **argv);
}
#endif /* CoinParam_H */

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
