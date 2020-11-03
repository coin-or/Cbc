// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include "CbcConfig.h"

#include <string>
#include <sstream>
#include <iostream>
#include <cassert>

#include "CoinParam.hpp"

#include "ClpFactorization.hpp"

#include "CbcParam.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CbcParam::CbcParam()
    : type_(CBC_PARAM_NOTUSED_INVALID), lowerDoubleValue_(0.0),
      upperDoubleValue_(0.0), lowerIntValue_(0), upperIntValue_(0),
      lengthName_(0), lengthMatch_(0), definedKeyWords_(), name_(),
      shortHelp_(), longHelp_(), action_(CBC_PARAM_NOTUSED_INVALID),
      currentKeyWord_(-1), display_(0), intValue_(-1), doubleValue_(-1.0),
      stringValue_(""), whereUsed_(7), fakeKeyWord_(-1), fakeValue_(0) {}
// Other constructors
CbcParam::CbcParam(std::string name, std::string help, double lower,
                   double upper, CbcParameterType type, int display)
    : type_(type), lowerIntValue_(0), upperIntValue_(0), definedKeyWords_(),
      name_(name), shortHelp_(help), longHelp_(), action_(type),
      currentKeyWord_(-1), display_(display), intValue_(-1), doubleValue_(-1.0),
      stringValue_(""), whereUsed_(7), fakeKeyWord_(-1), fakeValue_(0) {
  lowerDoubleValue_ = lower;
  upperDoubleValue_ = upper;
  gutsOfConstructor();
}
CbcParam::CbcParam(std::string name, std::string help, int lower, int upper,
                   CbcParameterType type, int display)
    : type_(type), lowerDoubleValue_(0.0), upperDoubleValue_(0.0),
      definedKeyWords_(), name_(name), shortHelp_(help), longHelp_(),
      action_(type), currentKeyWord_(-1), display_(display), intValue_(-1),
      doubleValue_(-1.0), stringValue_(""), whereUsed_(7), fakeKeyWord_(-1),
      fakeValue_(0) {
  gutsOfConstructor();
  lowerIntValue_ = lower;
  upperIntValue_ = upper;
}
// Other strings will be added by append
CbcParam::CbcParam(std::string name, std::string help, std::string firstValue,
                   CbcParameterType type, int whereUsed, int display)
    : type_(type), lowerDoubleValue_(0.0), upperDoubleValue_(0.0),
      lowerIntValue_(0), upperIntValue_(0), definedKeyWords_(), name_(name),
      shortHelp_(help), longHelp_(), action_(type), currentKeyWord_(0),
      display_(display), intValue_(-1), doubleValue_(-1.0), stringValue_(""),
      whereUsed_(whereUsed), fakeKeyWord_(-1), fakeValue_(0) {
  gutsOfConstructor();
  definedKeyWords_.push_back(firstValue);
}
// Action
CbcParam::CbcParam(std::string name, std::string help, CbcParameterType type,
                   int whereUsed, int display)
    : type_(type), lowerDoubleValue_(0.0), upperDoubleValue_(0.0),
      lowerIntValue_(0), upperIntValue_(0), definedKeyWords_(), name_(name),
      shortHelp_(help), longHelp_(), action_(type), currentKeyWord_(-1),
      display_(display), intValue_(-1), doubleValue_(-1.0), stringValue_(""),
      fakeKeyWord_(-1), fakeValue_(0) {
  whereUsed_ = whereUsed;
  gutsOfConstructor();
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CbcParam::CbcParam(const CbcParam &rhs) {
  type_ = rhs.type_;
  lowerDoubleValue_ = rhs.lowerDoubleValue_;
  upperDoubleValue_ = rhs.upperDoubleValue_;
  lowerIntValue_ = rhs.lowerIntValue_;
  upperIntValue_ = rhs.upperIntValue_;
  lengthName_ = rhs.lengthName_;
  lengthMatch_ = rhs.lengthMatch_;
  definedKeyWords_ = rhs.definedKeyWords_;
  name_ = rhs.name_;
  shortHelp_ = rhs.shortHelp_;
  longHelp_ = rhs.longHelp_;
  action_ = rhs.action_;
  currentKeyWord_ = rhs.currentKeyWord_;
  display_ = rhs.display_;
  intValue_ = rhs.intValue_;
  doubleValue_ = rhs.doubleValue_;
  stringValue_ = rhs.stringValue_;
  whereUsed_ = rhs.whereUsed_;
  fakeKeyWord_ = rhs.fakeKeyWord_;
  fakeValue_ = rhs.fakeValue_;
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CbcParam::~CbcParam() {}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CbcParam &CbcParam::operator=(const CbcParam &rhs) {
  if (this != &rhs) {
    type_ = rhs.type_;
    lowerDoubleValue_ = rhs.lowerDoubleValue_;
    upperDoubleValue_ = rhs.upperDoubleValue_;
    lowerIntValue_ = rhs.lowerIntValue_;
    upperIntValue_ = rhs.upperIntValue_;
    lengthName_ = rhs.lengthName_;
    lengthMatch_ = rhs.lengthMatch_;
    definedKeyWords_ = rhs.definedKeyWords_;
    name_ = rhs.name_;
    shortHelp_ = rhs.shortHelp_;
    longHelp_ = rhs.longHelp_;
    action_ = rhs.action_;
    currentKeyWord_ = rhs.currentKeyWord_;
    display_ = rhs.display_;
    intValue_ = rhs.intValue_;
    doubleValue_ = rhs.doubleValue_;
    stringValue_ = rhs.stringValue_;
    whereUsed_ = rhs.whereUsed_;
    fakeKeyWord_ = rhs.fakeKeyWord_;
    fakeValue_ = rhs.fakeValue_;
  }
  return *this;
}
void CbcParam::gutsOfConstructor() {
  std::string::size_type shriekPos = name_.find('!');
  lengthName_ = static_cast<unsigned int>(name_.length());
  if (shriekPos == std::string::npos) {
    // does not contain '!'
    lengthMatch_ = lengthName_;
  } else {
    lengthMatch_ = static_cast<unsigned int>(shriekPos);
    name_ = name_.substr(0, shriekPos) + name_.substr(shriekPos + 1);
    lengthName_--;
  }
}

// Insert string (only valid for keywords)
void CbcParam::append(std::string keyWord) {
  definedKeyWords_.push_back(keyWord);
}

double CbcParam::doubleParameter(OsiSolverInterface *model) const {
  double value = 0.0;
  bool getDblParamRetValue;
  switch (type_) {
  case CLP_PARAM_DBL_DUALTOLERANCE:
    getDblParamRetValue = model->getDblParam(OsiDualTolerance, value);
    assert(getDblParamRetValue);
    break;
  case CLP_PARAM_DBL_PRIMALTOLERANCE:
    getDblParamRetValue = model->getDblParam(OsiPrimalTolerance, value);
    assert(getDblParamRetValue);
    break;
  default:
    return doubleValue_;
    break;
  }
  return value;
}
double CbcParam::doubleParameter(CbcModel &model) const {
  double value;
  switch (type_) {
  case CBC_PARAM_DBL_INFEASIBILITYWEIGHT:
    value = model.getDblParam(CbcModel::CbcInfeasibilityWeight);
    break;
  case CBC_PARAM_DBL_INTEGERTOLERANCE:
    value = model.getDblParam(CbcModel::CbcIntegerTolerance);
    break;
  case CBC_PARAM_DBL_INCREMENT:
    value = model.getDblParam(CbcModel::CbcCutoffIncrement);
    break;
  case CBC_PARAM_DBL_ALLOWABLEGAP:
    value = model.getDblParam(CbcModel::CbcAllowableGap);
    break;
  case CBC_PARAM_DBL_GAPRATIO:
    value = model.getDblParam(CbcModel::CbcAllowableFractionGap);
    break;
  case CBC_PARAM_DBL_CUTOFF:
    value = model.getCutoff();
    break;
  case CBC_PARAM_DBL_TIMELIMIT_BAB:
    value = model.getDblParam(CbcModel::CbcMaximumSeconds);
    break;
  case CBC_PARAM_DBL_MAXSECONDSNIFS:
    value = model.getDblParam(CbcModel::CbcMaxSecondsNotImproving);
    break;

  case CLP_PARAM_DBL_DUALTOLERANCE:
  case CLP_PARAM_DBL_PRIMALTOLERANCE:
    value = doubleParameter(model.solver());
    break;
  default:
    value = doubleValue_;
    break;
  }
  return value;
}
int CbcParam::setDoubleParameter(OsiSolverInterface *model, double value,
                                 bool doPrinting) {
  int returnCode;
  std::string message =
     setDoubleParameterWithMessage(model, value, returnCode);
  if (doPrinting){
    std::cout << message << std::endl;
  }
  return returnCode;
}
int CbcParam::setDoubleParameter(CbcModel &model, double value,
                                 bool doPrinting) {
  int returnCode = 0;
  std::string message =
     setDoubleParameterWithMessage(model, value, returnCode);
  if (doPrinting){
    std::cout << message << std::endl;
  }
  return returnCode;
}
// Sets double parameter and returns printable string and error code
std::string CbcParam::setDoubleParameterWithMessage(OsiSolverInterface *model,
                                                    double value,
                                                    int &returnCode) {
  std::ostringstream buffer;
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
     buffer << value << " was provided for " << name_;
     buffer << " - valid range is " << lowerDoubleValue_;
     buffer << " to " << upperDoubleValue_ << std::endl;
     returnCode = 1;
  } else {
    double oldValue = doubleValue_;
    buffer << name_ << " was changed from ";
    buffer << oldValue << " to " << value << std::endl;
    returnCode = 0;
    doubleValue_ = value;
    switch (type_) {
    case CLP_PARAM_DBL_DUALTOLERANCE:
      model->getDblParam(OsiDualTolerance, oldValue);
      model->setDblParam(OsiDualTolerance, value);
      break;
    case CLP_PARAM_DBL_PRIMALTOLERANCE:
      model->getDblParam(OsiPrimalTolerance, oldValue);
      model->setDblParam(OsiPrimalTolerance, value);
      break;
    default:
      break;
    }
  }
  return buffer.str();
}
// TODO: This function doesn't set all the possible parameters. Does this
// matter?
// Sets double parameter and returns printable string and error code
std::string CbcParam::setDoubleParameterWithMessage(CbcModel &model,
                                                    double value,
                                                    int &returnCode) {
  std::ostringstream buffer;
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
     buffer << value << " was provided for " << name_;
     buffer << " - valid range is " << lowerDoubleValue_;
     buffer << " to " << upperDoubleValue_ << std::endl;
     returnCode = 1;
  } else {
    double oldValue = doubleValue_;
    buffer << name_ << " was changed from ";
    buffer << oldValue << " to " << value << std::endl;
    returnCode = 0;
    doubleValue_ = value;
    switch (type_) {
    case CBC_PARAM_DBL_INFEASIBILITYWEIGHT:
      oldValue = model.getDblParam(CbcModel::CbcInfeasibilityWeight);
      model.setDblParam(CbcModel::CbcInfeasibilityWeight, value);
      break;
    case CBC_PARAM_DBL_INTEGERTOLERANCE:
      oldValue = model.getDblParam(CbcModel::CbcIntegerTolerance);
      model.setDblParam(CbcModel::CbcIntegerTolerance, value);
      break;
    case CBC_PARAM_DBL_INCREMENT:
      oldValue = model.getDblParam(CbcModel::CbcCutoffIncrement);
      model.setDblParam(CbcModel::CbcCutoffIncrement, value);
    case CBC_PARAM_DBL_ALLOWABLEGAP:
      oldValue = model.getDblParam(CbcModel::CbcAllowableGap);
      model.setDblParam(CbcModel::CbcAllowableGap, value);
      break;
    case CBC_PARAM_DBL_GAPRATIO:
      oldValue = model.getDblParam(CbcModel::CbcAllowableFractionGap);
      model.setDblParam(CbcModel::CbcAllowableFractionGap, value);
      break;
    case CBC_PARAM_DBL_CUTOFF:
      oldValue = model.getCutoff();
      model.setCutoff(value);
      break;
    case CBC_PARAM_DBL_TIMELIMIT_BAB:
      oldValue = model.getDblParam(CbcModel::CbcMaximumSeconds);
      model.setDblParam(CbcModel::CbcMaximumSeconds, value);
      break;
    case CBC_PARAM_DBL_MAXSECONDSNIFS:
      oldValue = model.getDblParam(CbcModel::CbcMaxSecondsNotImproving);
      model.setDblParam(CbcModel::CbcMaxSecondsNotImproving, value);
      break;
    case CLP_PARAM_DBL_DUALTOLERANCE:
    case CLP_PARAM_DBL_PRIMALTOLERANCE:
      setDoubleParameter(model.solver(), value);
      return 0; // to avoid message
    default:
      break;
    }
  }
  return buffer.str();
}
void CbcParam::setDoubleValue(double value) {
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is "
              << lowerDoubleValue_ << " to " << upperDoubleValue_ << std::endl;
  } else {
    doubleValue_ = value;
  }
}
std::string CbcParam::setDoubleValueWithMessage(double value) {
  std::ostringstream buffer;
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
     buffer << value << " was provided for " << name_;
     buffer << " - valid range is " << lowerDoubleValue_;
     buffer << " to " << upperDoubleValue_ << std::endl;
  } else {
    double oldValue = doubleValue_;
    buffer << name_ << " was changed from ";
    buffer << oldValue << " to " << value << std::endl;
    doubleValue_ = value;
  }
  return buffer.str();
}
int CbcParam::checkDoubleParameter(double value) const {
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is "
              << lowerDoubleValue_ << " to " << upperDoubleValue_ << std::endl;
    return 1;
  } else {
    return 0;
  }
}
int CbcParam::intParameter(OsiSolverInterface *model) const {
  int value = 0;
  switch (type_) {
  case CBC_PARAM_INT_LPLOGLEVEL:
    // value=model->logLevel();
    break;
  default:
    abort();
  }
  return value;
}
int CbcParam::intParameter(CbcModel &model) const {
  int value;
  switch (type_) {
  case CBC_PARAM_INT_SOLVERLOGLEVEL:
    value = model.messageHandler()->logLevel();
    break;
  case CBC_PARAM_INT_LPLOGLEVEL:
    value = model.solver()->messageHandler()->logLevel();
    break;
  case CBC_PARAM_INT_MAXNODES:
    value = model.getIntParam(CbcModel::CbcMaxNumNode);
    break;
  case CBC_PARAM_INT_MAXNODESNOTIMPROVINGFS:
    value = model.getIntParam(CbcModel::CbcMaxNodesNotImproving);
    break;
  case CBC_PARAM_INT_MAXSOLS:
    value = model.getIntParam(CbcModel::CbcMaxNumSol);
    break;
  case CBC_PARAM_INT_MAXSAVEDSOLS:
    value = model.maximumSavedSolutions();
    break;
  case CBC_PARAM_INT_STRONGBRANCHING:
    value = model.numberStrong();
    break;
  case CBC_PARAM_INT_NUMBERBEFORE:
    value = model.numberBeforeTrust();
    break;
  case CBC_PARAM_INT_NUMBERANALYZE:
    value = model.numberAnalyzeIterations();
    break;
  case CBC_PARAM_INT_CUTPASSINTREE:
    value = model.getMaximumCutPasses();
    break;
  case CBC_PARAM_INT_CUTPASS:
    value = model.getMaximumCutPassesAtRoot();
    break;
#ifdef COIN_HAS_CBC
#ifdef CBC_THREAD
  case CBC_PARAM_INT_THREADS:
    value = model.getNumberThreads();
    break;
#endif
  case CBC_PARAM_INT_RANDOMSEED:
    value = model.getRandomSeed();
    break;
#endif
  default:
    value = intValue_;
    break;
  }
  return value;
}
int CbcParam::setIntParameter(OsiSolverInterface *model, int value,
                              bool doPrinting) {
  int returnCode;
  std::string message =
     setIntParameterWithMessage(model, value, returnCode);
  if (doPrinting){
    std::cout << message << std::endl;
  }
  return returnCode;
}
int CbcParam::setIntParameter(CbcModel &model, int value,
                              bool doPrinting) {
  int returnCode;
  std::string message =
     setIntParameterWithMessage(model, value, returnCode);
  if (doPrinting){
    std::cout << message << std::endl;
  }
  return returnCode;
}
// Sets int parameter and returns printable string and error code
std::string CbcParam::setIntParameterWithMessage(OsiSolverInterface *model,
                                                 int value, int &returnCode) {
  std::ostringstream buffer;
  if (value < lowerIntValue_ || value > upperIntValue_) {
     buffer << value << " was provided for " << name_;
     buffer << " - valid range is " << lowerIntValue_ << " to ";
     buffer << upperIntValue_ << std::endl;
     returnCode = 1;
  } else {
    int oldValue = intValue_;
    buffer << name_ << " was changed from " << oldValue;
    buffer << " to " << value << std::endl;
    returnCode = 0;
    intValue_ = value;
    switch (type_) {
    case CBC_PARAM_INT_LPLOGLEVEL:
      model->messageHandler()->setLogLevel(value);
      break;
    default:
      break;
    }
  }
  return buffer.str();
}
// Sets int parameter and returns printable string and error code
std::string CbcParam::setIntParameterWithMessage(CbcModel &model, int value,
                                                 int &returnCode) {
  std::ostringstream buffer;
  if (value < lowerIntValue_ || value > upperIntValue_) {
     buffer << value << " was provided for " << name_;
     buffer << " - valid range is " << lowerIntValue_ << " to ";
     buffer << upperIntValue_ << std::endl;
     returnCode = 1;
  } else {
    int oldValue = intValue_;
    buffer << name_ << " was changed from " << oldValue;
    buffer << " to " << value << std::endl;
    returnCode = 0;
    intValue_ = value;
    switch (type_) {
    case CBC_PARAM_INT_SOLVERLOGLEVEL:
      oldValue = model.messageHandler()->logLevel();
      model.messageHandler()->setLogLevel(CoinAbs(value));
      break;
    case CBC_PARAM_INT_LPLOGLEVEL:
      oldValue = model.solver()->messageHandler()->logLevel();
      model.solver()->messageHandler()->setLogLevel(value);
      break;
    case CBC_PARAM_INT_MAXNODES:
      oldValue = model.getIntParam(CbcModel::CbcMaxNumNode);
      model.setIntParam(CbcModel::CbcMaxNumNode, value);
      break;
    case CBC_PARAM_INT_MAXNODESNOTIMPROVINGFS:
      oldValue = model.getIntParam(CbcModel::CbcMaxNodesNotImproving);
      model.setIntParam(CbcModel::CbcMaxNodesNotImproving, value);
      break;

    case CBC_PARAM_INT_MAXSOLS:
      oldValue = model.getIntParam(CbcModel::CbcMaxNumSol);
      model.setIntParam(CbcModel::CbcMaxNumSol, value);
      break;
    case CBC_PARAM_INT_MAXSAVEDSOLS:
      oldValue = model.maximumSavedSolutions();
      model.setMaximumSavedSolutions(value);
      break;
    case CBC_PARAM_INT_STRONGBRANCHING:
      oldValue = model.numberStrong();
      model.setNumberStrong(value);
      break;
    case CBC_PARAM_INT_NUMBERBEFORE:
      oldValue = model.numberBeforeTrust();
      model.setNumberBeforeTrust(value);
      break;
    case CBC_PARAM_INT_NUMBERANALYZE:
      oldValue = model.numberAnalyzeIterations();
      model.setNumberAnalyzeIterations(value);
      break;
    case CBC_PARAM_INT_CUTPASSINTREE:
      oldValue = model.getMaximumCutPasses();
      model.setMaximumCutPasses(value);
      break;
    case CBC_PARAM_INT_CUTPASS:
      oldValue = model.getMaximumCutPassesAtRoot();
      model.setMaximumCutPassesAtRoot(value);
      break;
#ifdef COIN_HAS_CBC
#ifdef CBC_THREAD
    case CBC_PARAM_INT_THREADS:
      oldValue = model.getNumberThreads();
      model.setNumberThreads(value);
      break;
#endif
    case CBC_PARAM_INT_RANDOMSEED:
      oldValue = model.getRandomSeed();
      model.setRandomSeed(value);
      break;
#endif
    default:
      break;
    }
  }
  return buffer.str();
}
void CbcParam::setIntValue(int value) {
  if (value < lowerIntValue_ || value > upperIntValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is "
              << lowerIntValue_ << " to " << upperIntValue_ << std::endl;
  } else {
    intValue_ = value;
  }
}
std::string CbcParam::setIntValueWithMessage(int value) {
   std::ostringstream buffer;
   if (value < lowerIntValue_ || value > upperIntValue_) {
      buffer << value << " was provided for " << name_;
      buffer << " - valid range is " << lowerIntValue_ << " to ";
      buffer << upperIntValue_ << std::endl;
   } else {
      int oldValue = intValue_;
      buffer << name_ << " was changed from " << oldValue;
      buffer << " to " << value << std::endl;
      intValue_ = value;
   }
   return buffer.str();
}
void CbcParam::setStringValue(std::string value) { stringValue_ = value; }

// Returns name which could match
std::string CbcParam::matchName() const {
  if (lengthMatch_ == lengthName_)
    return name_;
  else
    return name_.substr(0, lengthMatch_) + "(" + name_.substr(lengthMatch_) +
           ")";
}

// Returns length of name for printing
int CbcParam::lengthMatchName() const {
  if (lengthName_ == lengthMatch_)
    return lengthName_;
  else
    return lengthName_ + 2;
}

int CbcParam::matches(std::string input) const {
  // look up strings to do more elegantly
  if (input.length() > lengthName_) {
    return 0;
  } else {
    unsigned int i;
    for (i = 0; i < input.length(); i++) {
      if (tolower(name_[i]) != tolower(input[i]))
        break;
    }
    if (i < input.length()) {
      return 0;
    } else if (i >= lengthMatch_) {
      return 1;
    } else {
      // matched but too short
      return 2;
    }
  }
}

// Returns parameter option which matches (-1 if none)
int CbcParam::parameterOption(std::string check) const {
  int numberItems = static_cast<int>(definedKeyWords_.size());
  if (!numberItems) {
    return -1;
  } else {
    int whichItem = 0;
    unsigned int it;
    for (it = 0; it < definedKeyWords_.size(); it++) {
      std::string thisOne = definedKeyWords_[it];
      std::string::size_type shriekPos = thisOne.find('!');
      size_t length1 = thisOne.length();
      size_t length2 = length1;
      if (shriekPos != std::string::npos) {
        // contains '!'
        length2 = shriekPos;
        thisOne = thisOne.substr(0, shriekPos) + thisOne.substr(shriekPos + 1);
        length1 = thisOne.length();
      }
      if (check.length() <= length1 && length2 <= check.length()) {
        unsigned int i;
        for (i = 0; i < check.length(); i++) {
          if (tolower(thisOne[i]) != tolower(check[i]))
            break;
        }
        if (i < check.length()) {
          whichItem++;
        } else if (i >= length2) {
          break;
        }
      } else {
        whichItem++;
      }
    }
    if (whichItem < numberItems) {
      return whichItem;
    } else {
      if (fakeKeyWord_ <= 0)
        return -1;
      // allow plus or minus
      int n;
      if (check.substr(0, 4) == "plus" || check.substr(0, 4) == "PLUS") {
        n = 4;
      } else if (check.substr(0, 5) == "minus" ||
                 check.substr(0, 5) == "MINUS") {
        n = 5;
      } else {
        return -1;
      }
      int value = 0;
      std::string field = check.substr(n);
      if (field != "EOL") {
        const char *start = field.c_str();
        char *endPointer = NULL;
        // check valid
        value = static_cast<int>(strtol(start, &endPointer, 10));
        if (*endPointer != '\0') {
          return -1;
        }
        if (n == 4)
          return value + 1000;
        else
          return -value - 1000;
      } else {
        return -1;
      }
    }
  }
}
// Prints parameter options
void CbcParam::printOptions() const {
  std::cout << "<Possible options for " << name_ << " are:";
  unsigned int it;
  for (it = 0; it < definedKeyWords_.size(); it++) {
    std::string thisOne = definedKeyWords_[it];
    std::string::size_type shriekPos = thisOne.find('!');
    if (shriekPos != std::string::npos) {
      // contains '!'
      thisOne = thisOne.substr(0, shriekPos) + "(" +
                thisOne.substr(shriekPos + 1) + ")";
    }
    std::cout << " " << thisOne;
  }
  assert(currentKeyWord_ >= 0 &&
         currentKeyWord_ < static_cast<int>(definedKeyWords_.size()));
  std::string current = definedKeyWords_[currentKeyWord_];
  std::string::size_type shriekPos = current.find('!');
  if (shriekPos != std::string::npos) {
    // contains '!'
    current = current.substr(0, shriekPos) + "(" +
              current.substr(shriekPos + 1) + ")";
  }
  std::cout << ";\n\tcurrent  " << current << ">" << std::endl;
}

// Sets current parameter option using string
void CbcParam::setCurrentOption(const std::string value) {
  int action = parameterOption(value);
  if (action >= 0)
    currentKeyWord_ = action;
#if FLUSH_PRINT_BUFFER > 2
  if (name_ == "bufferedMode")
    coinFlushBufferFlag = action;
#endif
}
// Sets current parameter option
void CbcParam::setCurrentOption(int value, bool printIt) {
  if (printIt && value != currentKeyWord_)
    std::cout << "Option for " << name_ << " changed from "
              << definedKeyWords_[currentKeyWord_] << " to "
              << definedKeyWords_[value] << std::endl;

#if FLUSH_PRINT_BUFFER > 2
  if (name_ == "bufferedMode")
    coinFlushBufferFlag = value;
#endif
  currentKeyWord_ = value;
}
// Sets current parameter option and returns printable string
std::string CbcParam::setCurrentOptionWithMessage(int value) {
  std::ostringstream buffer;
  if (value != currentKeyWord_) {
    buffer << "Option for " << name_ << " changed from ";
    char current[100];
    char newString[100];
    if (currentKeyWord_ >= 0 && (fakeKeyWord_ <= 0 ||
                                 currentKeyWord_ < fakeKeyWord_)){
       buffer << definedKeyWords_[currentKeyWord_];
    }else if (currentKeyWord_ < 0){
       buffer << "minus" << -currentKeyWord_ - 1000;
    }else{
       buffer << "plus" << currentKeyWord_ - 1000;
    }
    buffer << " to ";
    if (value >= 0 && (fakeKeyWord_ <= 0 || value < fakeKeyWord_)){
       buffer << definedKeyWords_[value];
    }else if (value < 0){
       buffer << "minus" << -value - 1000;
    }else{
       buffer << "plus" << value - 1000;
    }

#if FLUSH_PRINT_BUFFER > 2
    if (name_ == "bufferedMode")
      coinFlushBufferFlag = value;
#endif
    currentKeyWord_ = value;
  }
  return buffer.str();
}
// Sets current parameter option using string with message
std::string CbcParam::setCurrentOptionWithMessage(const std::string value) {
  std::ostringstream buffer;
  int action = parameterOption(value);
  char current[100];
  if (action >= 0) {
#if FLUSH_PRINT_BUFFER > 2
    if (name_ == "bufferedMode")
      coinFlushBufferFlag = action;
#endif
    buffer << "Option for " << name_ << " changed from ";
    if (currentKeyWord_ >= 0 && (fakeKeyWord_ <= 0 ||
                                 currentKeyWord_ < fakeKeyWord_)){
       buffer << definedKeyWords_[currentKeyWord_];
    }else if (currentKeyWord_ < 0){
       buffer << "minus" << -currentKeyWord_ - 1000;
    }else{
       buffer << "plus" << currentKeyWord_ - 1000;
    }
    buffer << " to " << value;
    currentKeyWord_ = action;
  } else {
     buffer << "Option for " << name_ << " given illegal value " << value;
  }
  return buffer.str();
}
/* Returns current parameter option position
   but if fake keyword returns fakeValue_
*/
int CbcParam::currentOptionAsInteger() const {
  int fakeInteger;
  return currentOptionAsInteger(fakeInteger);
}
/* Returns current parameter option position
   but if fake keyword returns fakeValue_ and sets
   fakeInteger to value
*/
int CbcParam::currentOptionAsInteger(int &fakeInteger) const {
  fakeInteger = -COIN_INT_MAX;
  if (fakeKeyWord_ < 0) {
    return currentKeyWord_;
  } else if (currentKeyWord_ >= 0 && currentKeyWord_ < fakeKeyWord_) {
    return currentKeyWord_;
  } else {
    // fake
    if (currentKeyWord_ < 0)
      fakeInteger = currentKeyWord_ + 1000;
    else
      fakeInteger = currentKeyWord_ - 1000;
    return fakeValue_;
  }
}

// TODO: Fix this
// Print Long help
void CbcParam::printLongHelp() const {
  if (type_ >= 1 && type_ < 600) {
    CoinPrintString(longHelp_);
    if (type_ < CLP_PARAM_INT_LOGLEVEL) {
      printf("<Range of values is %g to %g;\n\tcurrent %g>\n",
             lowerDoubleValue_, upperDoubleValue_, doubleValue_);
      assert(upperDoubleValue_ > lowerDoubleValue_);
    } else if (type_ < CLP_PARAM_STR_DIRECTION) {
      printf("<Range of values is %d to %d;\n\tcurrent %d>\n", lowerIntValue_,
             upperIntValue_, intValue_);
      assert(upperIntValue_ > lowerIntValue_);
    } else if (type_ < CLP_PARAM_ACTION_DIRECTORY) {
      printOptions();
    }
  }
}

void CbcParam::printString() const {
  if (name_ == "directory")
    std::cout << "Current working directory is " << stringValue_ << std::endl;
  else if (name_.substr(0, 6) == "printM")
    std::cout << "Current value of printMask is " << stringValue_ << std::endl;
  else
    std::cout << "Current default (if $ as parameter) for " << name_ << " is "
              << stringValue_ << std::endl;
}

// Sets value of fake keyword to current size of keywords
void CbcParam::setFakeKeyWord(int fakeValue) {
  fakeKeyWord_ = static_cast<int>(definedKeyWords_.size());
  assert(fakeKeyWord_ > 0);
  fakeValue_ = fakeValue;
  assert(fakeValue_ >= 0);
}

//###########################################################################
//###########################################################################

/*
  Subroutine to establish the cbc parameter array. See the description of
  class CbcParam for details.
*/
void establishCbcParams(std::vector<CbcParam> &parameters) {
  parameters.clear();
  parameters.push_back(CbcParam("?", "For help", CBC_PARAM_GENERALQUERY, 7, 0));
  parameters.push_back(
      CbcParam("???", "For help", CBC_PARAM_FULLGENERALQUERY, 7, 0));
  parameters.push_back(
      CbcParam("-", "From stdin", CBC_PARAM_ACTION_STDIN, 3, 0));

// some help strings that repeat for many options
#define CUTS_LONGHELP                                                          \
  "Value 'on' enables the cut generator and CBC will try it in the branch "    \
  "and cut tree (see cutDepth on how to fine tune the behavior). "             \
  "Value 'root' lets CBC run the cut generator generate only at the root "     \
  "node. "                                                                     \
  "Value 'ifmove' lets CBC use the cut generator in the tree if it looks as "  \
  "if it is doing some good and moves the objective value. "                   \
  "Value 'forceon' turns on the cut generator and forces CBC to use it at "    \
  "every node."
#define HEURISTICS_LONGHELP                                                    \
  "Value 'on' means to use the heuristic in each node of the tree, i.e. "      \
  "after preprocessing. "                                                      \
  "Value 'before' means use the heuristic only if option doHeuristics is "     \
  "used. "                                                                     \
  "Value 'both' means to use the heuristic if option doHeuristics is used "    \
  "and during solve."

  {
    CbcParam p("allC!ommands", "Whether to print less used commands", "no",
               CBC_PARAM_STR_ALLCOMMANDS);

    p.append("more");
    p.append("all");
    p.setLonghelp(
        "For the sake of your sanity, only the more useful and simple commands \
are printed out on ?.");
    parameters.push_back(p);
  }
  {
    CbcParam p("allow!ableGap", "Stop when gap between best possible and \
best less than this",
               0.0, COIN_DBL_MAX, CBC_PARAM_DBL_ALLOWABLEGAP);
    p.setDoubleValue(1e-10);
    p.setLonghelp(
        "If the gap between best known solution and the best possible solution is less than this \
value, then the search will be terminated.  Also see ratioGap.");
    parameters.push_back(p);
  }
  {
    CbcParam p("artif!icialCost",
               "Costs >= this treated as artificials in feasibility pump", 0.0,
               COIN_DBL_MAX, CBC_PARAM_DBL_ARTIFICIALCOST, 1);
    p.setDoubleValue(0.0);
    p.setLonghelp("A value of 0.0 means off. Otherwise, variables with costs "
                  ">= this are treated as artificial variables and fixed to "
                  "lower bound in feasibility pump.");
    parameters.push_back(p);
  }
  {
    CbcParam p("branch!AndCut", "Do Branch and Cut", CBC_PARAM_ACTION_BAB);
    p.setLonghelp(
        "This does branch and cut.  There are many parameters which can affect the performance.  \
First just try with default settings and look carefully at the log file.  Did cuts help?  Did they take too long?  \
Look at output to see which cuts were effective and then do some tuning.  You will see that the \
options for cuts are off, on, root and ifmove, forceon.  Off is \
obvious. " CUTS_LONGHELP
        " For probing, forceonbut just does fixing probing in tree - not strengthening etc.  \
If pre-processing reduced the size of the \
problem or strengthened many coefficients then it is probably wise to leave it on.  Switch off heuristics \
which did not provide solutions.  The other major area to look at is the search.  Hopefully good solutions \
were obtained fairly early in the search so the important point is to select the best variable to branch on.  \
See whether strong branching did a good job - or did it just take a lot of iterations.  Adjust the strongBranching \
and trustPseudoCosts parameters.  If cuts did a good job, then you may wish to \
have more rounds of cuts - see passC!uts and passT!ree.");
    parameters.push_back(p);
  }
  {
    CbcParam p("combine!Solutions", "Whether to use combine solution heuristic",
               "off", CBC_PARAM_STR_COMBINE);

    p.append("on");
    p.append("both");
    p.append("before");
    p.append("onquick");
    p.append("bothquick");
    p.append("beforequick");
    p.setLonghelp("This heuristic does branch and cut on given problem by just \
using variables which have appeared in one or more solutions. \
It is obviously only tried after two or more solutions have been found. " HEURISTICS_LONGHELP);

    parameters.push_back(p);
  }
  {
    CbcParam p("combine2!Solutions",
               "Whether to use crossover solution heuristic", "off",
               CBC_PARAM_STR_CROSSOVER2);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp("This heuristic does branch and cut on the problem given by \
fixing variables which have the same value in two or more solutions. \
It obviously only tries after two or more solutions. " HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("constraint!fromCutoff", "Whether to use cutoff as constraint",
               "off", CBC_PARAM_STR_CUTOFF_CONSTRAINT);

    p.append("on");
    p.append("variable");
    p.append("forcevariable");
    p.append("conflict");
    p.setLonghelp(
        "For some problems, cut generators and general branching work better "
        "if the problem would be infeasible if the cost is too high. "
        "If this option is enabled, the objective function is added as a "
        "constraint which right hand side is set to the current cutoff value "
        "(objective value of best known solution)");
    parameters.push_back(p);
  }
  {
    CbcParam p("cost!Strategy", "How to use costs for branching priorities",
               "off", CBC_PARAM_STR_BRANCHPRIORITY);

    p.append("pri!orities");
    p.append("column!Order?");
    p.append("01f!irst?");
    p.append("01l!ast?");
    p.append("length!?");
    p.append("singletons");
    p.append("nonzero");
    p.append("general!Force?");
    p.setLonghelp(
        "Value 'priorities' assigns highest priority to variables with largest "
        "absolute cost. This primitive strategy can be surprisingly effective. "
        "Value 'columnorder' assigns the priorities 1, 2, 3, ... with respect "
        "to the column ordering. "
        "Value '01first' ('01last') assignes two sets of priorities such that "
        "binary variables get high (low) priority. "
        "Value 'length' assigns high priority to variables that occur in many "
        "equations. ");

    parameters.push_back(p);
  }
  {
    CbcParam p("cplex!Use", "Whether to use Cplex!", "off", CBC_PARAM_STR_CPX);
    p.append("on");
    p.setLonghelp(
        " If the user has Cplex, but wants to use some of Cbc's heuristics \
then you can!  If this is on, then Cbc will get to the root node and then \
hand over to Cplex.  If heuristics find a solution this can be significantly \
quicker.  You will probably want to switch off Cbc's cuts as Cplex thinks \
they are genuine constraints.  It is also probable that you want to switch \
off preprocessing, although for difficult problems it is worth trying \
both.");
    parameters.push_back(p);
  }
  {
    CbcParam p("csv!Statistics", "Create one line of statistics",
               CBC_PARAM_ACTION_CSVSTATISTICS, 2, 1);
    p.setLonghelp(
        "This appends statistics to given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  Adds header if file empty or does not exist.");
    parameters.push_back(p);
  }
  {
    CbcParam p("cutD!epth", "Depth in tree at which to do cuts", -1,
               COIN_INT_MAX, CBC_PARAM_INT_CUTDEPTH);
    p.setLonghelp(
        "Cut generators may be off, on, on only at the root node, or on if they look useful. \
      Setting this option to a positive value K let CBC call a cutgenerator on a node whenever the depth in the tree is a multiple of K. \
      The default of -1 lets CBC decide.");
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcParam p("cutL!ength", "Length of a cut", -1, COIN_INT_MAX,
               CBC_PARAM_INT_CUTLENGTH);
    p.setLonghelp(
        "At present this only applies to Gomory cuts. -1 (default) leaves as is. \
Any value >0 says that all cuts <= this length can be generated both at \
root node and in tree. 0 says to use some dynamic lengths.  If value >=10,000,000 \
then the length in tree is value%10000000 - so 10000100 means unlimited length \
at root and 100 in tree.");
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcParam p("cuto!ff", "Bound on the objective value for all solutions",
               -COIN_DBL_MAX, COIN_DBL_MAX, CBC_PARAM_DBL_CUTOFF);
    p.setDoubleValue(1.0e50);
    p.setLonghelp(
        "All solutions must have a better objective value (in a minimization sense) than the value of this option.  \
CBC also updates this value whenever it obtains a solution to the value of \
the objective function of the solution minus the cutoff increment.");
    parameters.push_back(p);
  }
  {
    CbcParam p("cuts!OnOff", "Switches all cut generators on or off", "off",
               CBC_PARAM_STR_CUTSSTRATEGY);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.setLonghelp(
        "This can be used to switch on or off all cut generators (apart from "
        "Reduce and Split). "
        "Then one can turn individual ones off or on. " CUTS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("debug!In", "read valid solution from file",
               CBC_PARAM_ACTION_DEBUG, 7, 1);

    p.setLonghelp(
        "This will read a solution file from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.\n\n\
If set to create it will create a file called debug.file  after search.\n\n\
The idea is that if you suspect a bad cut generator \
you can do a good run with debug set to 'create' and then switch on the cuts you suspect and \
re-run with debug set to 'debug.file'  The create case has same effect as saveSolution.");
    parameters.push_back(p);
  }
  {
    CbcParam p("depth!MiniBab", "Depth at which to try mini branch-and-bound",
               -COIN_INT_MAX, COIN_INT_MAX, CBC_PARAM_INT_DEPTHMINIBAB);

    p.setIntValue(-1);
    p.setLonghelp(
        "Rather a complicated parameter but can be useful. -1 means off for large problems but on as if -12 for problems where rows+columns<500, -2 \
means use Cplex if it is linked in.  Otherwise if negative then go into depth first complete search fast branch and bound when depth>= -value-2 \
(so -3 will use this at depth>=1).  This mode is only switched on after 500 nodes.  If you really want to switch it off for small problems then set \
this to -999.  If >=0 the value doesn't matter very much.  The code will do approximately 100 nodes of fast branch and bound every now and then at depth>=5.  \
The actual logic is too twisted to describe here.");
    parameters.push_back(p);
  }
  {
    CbcParam p("dextra3", "Extra double parameter 3", -COIN_DBL_MAX,
               COIN_DBL_MAX, CBC_PARAM_DBL_DEXTRA3, 0);
    p.setDoubleValue(0.0);
    parameters.push_back(p);
  }
  {
    CbcParam p("dextra4", "Extra double parameter 4", -COIN_DBL_MAX,
               COIN_DBL_MAX, CBC_PARAM_DBL_DEXTRA4, 0);
    p.setDoubleValue(0.0);
    parameters.push_back(p);
  }
  {
    CbcParam p("dextra4", "Extra double parameter 5", -COIN_DBL_MAX,
               COIN_DBL_MAX, CBC_PARAM_DBL_DEXTRA5, 0);
    p.setDoubleValue(0.0);
    parameters.push_back(p);
  }
  {
    CbcParam p("Dins", "Whether to try Distance Induced Neighborhood Search",
               "off", CBC_PARAM_STR_DINS);

    p.append("on");
    p.append("both");
    p.append("before");
    p.append("often");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("direction", "Minimize or Maximize", "min!imize",
               CBC_PARAM_STR_DIRECTION);
    p.append("max!imize");
    p.append("zero");
    p.setLonghelp(
        "The default is minimize - use 'direction maximize' for maximization.\n\
You can also use the parameters 'maximize' or 'minimize'.");
    parameters.push_back(p);
  }
  {
    CbcParam p("directory", "Set Default directory for import etc.",
               CBC_PARAM_ACTION_DIRECTORY);
    p.setLonghelp(
        "This sets the directory which import, export, saveModel, restoreModel etc will use.\
  It is initialized to './'");
    parameters.push_back(p);
  }
  {
    CbcParam p("dirSample",
               "Set directory where the COIN-OR sample problems are.",
               CBC_PARAM_ACTION_DIRSAMPLE, 7, 1);

    p.setLonghelp(
        "This sets the directory where the COIN-OR sample problems reside. It is\
 used only when -unitTest is passed to cbc. cbc will pick up the test problems\
 from this directory.\
 It is initialized to '../../Data/Sample'");
    parameters.push_back(p);
  }
  {
    CbcParam p("dirNetlib", "Set directory where the netlib problems are.",
               CBC_PARAM_ACTION_DIRNETLIB, 7, 1);

    p.setLonghelp(
        "This sets the directory where the netlib problems reside. One can get\
 the netlib problems from COIN-OR or from the main netlib site. This\
 parameter is used only when -netlib is passed to cbc. cbc will pick up the\
 netlib problems from this directory. If clp is built without zlib support\
 then the problems must be uncompressed.\
 It is initialized to '../../Data/Netlib'");
    parameters.push_back(p);
  }
  {
    CbcParam p("dirMiplib", "Set directory where the miplib 2003 problems are.",
               CBC_PARAM_ACTION_DIRMIPLIB, 7, 1);

    p.setLonghelp(
        "This sets the directory where the miplib 2003 problems reside. One can\
 get the miplib problems from COIN-OR or from the main miplib site. This\
 parameter is used only when -miplib is passed to cbc. cbc will pick up the\
 miplib problems from this directory. If cbc is built without zlib support\
 then the problems must be uncompressed.\
 It is initialized to '../../Data/miplib3'");
    parameters.push_back(p);
  }
  {
    CbcParam p("diveO!pt", "Diving options", -1, 200000, CBC_PARAM_INT_DIVEOPT,
               1);
    p.setLonghelp("If >2 && <20 then modify diving options - \
	 \n\t3 only at root and if no solution,  \
	 \n\t4 only at root and if this heuristic has not got solution, \
	 \n\t5 decay only if no solution, \
	 \n\t6 if depth <3 or decay, \
	 \n\t7 run up to 2 times if solution found 4 otherwise, \
	 \n\t>10 All only at root (DivingC normal as value-10), \
	 \n\t>20 All with value-20).");
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcParam p("diveS!olves", "Diving solve option", -1, 200000,
               CBC_PARAM_INT_DIVEOPTSOLVES, 1);

    p.setLonghelp(
        "If >0 then do up to this many solves. However, the last digit is ignored \
and used for extra options: \
      1-3 enables fixing of satisfied integer variables (but not at bound), \
      where 1 switches this off for that dive if the dive goes infeasible, \
      and 2 switches it off permanently if the dive goes infeasible.");
    p.setIntValue(100);
    parameters.push_back(p);
  }
  {
    CbcParam p("DivingS!ome", "Whether to try Diving heuristics", "off",
               CBC_PARAM_STR_DIVINGS);

    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(
        "This switches on a random diving heuristic at various times. \
One may prefer to individually turn diving heuristics on or off. " HEURISTICS_LONGHELP);
    // C - Coefficient, F - Fractional, G - Guided, L - LineSearch, P -
    // PseudoCost, V - VectorLength.
    parameters.push_back(p);
  }
  {
    CbcParam p("DivingC!oefficient",
               "Whether to try Coefficient diving heuristic", "off",
               CBC_PARAM_STR_DIVINGC);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("DivingF!ractional",
               "Whether to try Fractional diving heuristic", "off",
               CBC_PARAM_STR_DIVINGF);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("DivingG!uided", "Whether to try Guided diving heuristic", "off",
               CBC_PARAM_STR_DIVINGG);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("DivingL!ineSearch",
               "Whether to try Linesearch diving heuristic", "off",
               CBC_PARAM_STR_DIVINGL);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("DivingP!seudoCost",
               "Whether to try Pseudocost diving heuristic", "off",
               CBC_PARAM_STR_DIVINGP);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("DivingV!ectorLength",
               "Whether to try Vectorlength diving heuristic", "off",
               CBC_PARAM_STR_DIVINGV);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("doH!euristic", "Do heuristics before any preprocessing",
               CBC_PARAM_ACTION_DOHEURISTIC, 3);
    p.setLonghelp(
        "Normally heuristics are done in branch and bound.  It may be useful to do them outside. \
Only those heuristics with 'both' or 'before' set will run.  \
Doing this may also set cutoff, which can help with preprocessing.");
    parameters.push_back(p);
  }
  {
    CbcParam p("dw!Heuristic", "Whether to try Dantzig Wolfe heuristic", "off",
               CBC_PARAM_STR_DW);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(
        "This heuristic is very very compute intensive. It tries to find a "
        "Dantzig Wolfe structure and use that. " HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("end", "Stops clp execution", CBC_PARAM_ACTION_EXIT);
    p.setLonghelp(
        "This stops execution ; end, exit, quit and stop are synonyms");
    parameters.push_back(p);
  }
  {
    CbcParam p("environ!ment", "Read commands from environment",
               CBC_PARAM_ACTION_ENVIRONMENT, 7, 0);
    p.setLonghelp(
        "This starts reading from environment variable CBC_CLP_ENVIRONMENT.");
    parameters.push_back(p);
  }
  {
    CbcParam p("error!sAllowed", "Whether to allow import errors", "off",
               CBC_PARAM_STR_ERRORSALLOWED, 3);

    p.append("on");
    p.setLonghelp(
        "The default is not to use any model which had errors when reading the mps file.\
  Setting this to 'on' will allow all errors from which the code can recover\
 simply by ignoring the error.  There are some errors from which the code can not recover \
e.g. no ENDATA.  This has to be set before import i.e. -errorsAllowed on -import xxxxxx.mps.");
    parameters.push_back(p);
  }
  {
    CbcParam p("exper!iment", "Whether to use testing features", -1, 200000,
               CBC_PARAM_INT_EXPERIMENT, 0);
    p.setLonghelp("Defines how adventurous you want to be in using new ideas. \
0 then no new ideas, 1 fairly sensible, 2 a bit dubious, 3 you are on your own!");

    p.setIntValue(0);
    parameters.push_back(p);
  }
  {
    CbcParam p("expensive!Strong", "Whether to do even more strong branching",
               0, COIN_INT_MAX, CBC_PARAM_INT_STRONG_STRATEGY, 0);

    p.setLonghelp("Strategy for extra strong branching. \
0 is normal strong branching. \
1, 2, 4, and 6 does strong branching on all fractional variables if \
at the root node (1), \
at depth less than modifier (2), \
objective equals best possible (4), or \
at depth less than modifier and objective equals best possible (6). \
11, 12, 14, and 16 are like 1, 2, 4, and 6, respecitively, but do strong branching on all integer (incl. non-fractional) variables. \
Values >= 100 are used to specify a depth limit (value/100), otherwise 5 is used. \
If the values >= 100, then above rules are applied to value%100.");
    p.setIntValue(0);
    parameters.push_back(p);
  }
  {
    CbcParam p("export", "Export model as mps file", CBC_PARAM_ACTION_EXPORT);

    p.setLonghelp(
        "This will write an MPS format file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.mps'.  \
It can be useful to get rid of the original names and go over to using Rnnnnnnn and Cnnnnnnn.  This can be done by setting 'keepnames' off before importing mps file.");
    parameters.push_back(p);
  }
  {
    CbcParam p("extra1", "Extra integer parameter 1", -COIN_INT_MAX,
               COIN_INT_MAX, CBC_PARAM_INT_EXTRA1, 0);
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcParam p("extra2", "Extra integer parameter 2", -COIN_INT_MAX,
               COIN_INT_MAX, CBC_PARAM_INT_EXTRA2, 0);
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcParam p("extra3", "Extra integer parameter 3", -COIN_INT_MAX,
               COIN_INT_MAX, CBC_PARAM_INT_EXTRA3, 0);
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcParam p("extra4", "Extra integer parameter 4", -1, COIN_INT_MAX,
               CBC_PARAM_INT_EXTRA4, 0);

    p.setIntValue(-1);
    p.setLonghelp("This switches on yet more special options!! \
The bottom digit is a strategy when to used shadow price stuff e.g. 3 \
means use until a solution is found.  The next two digits say what sort \
of dual information to use.  After that it goes back to powers of 2 so -\n\
\n\t1000 - switches on experimental hotstart\n\
\n\t2,4,6000 - switches on experimental methods of stopping cuts\n\
\n\t8000 - increase minimum drop gradually\n\
\n\t16000 - switches on alternate gomory criterion");
    parameters.push_back(p);
  }
  {
    CbcParam p("extraV!ariables", "Allow creation of extra integer variables",
               -COIN_INT_MAX, COIN_INT_MAX, CBC_PARAM_INT_EXTRA_VARIABLES, 0);
    p.setIntValue(0);
    p.setLonghelp(
        "Switches on a trivial re-formulation that introduces extra integer "
        "variables to group together variables with same cost.");
    parameters.push_back(p);
  }
  {
    CbcParam p("feas!ibilityPump",
               "Whether to try the Feasibility Pump heuristic", "off",
               CBC_PARAM_STR_FPUMP);

    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp("This heuristic is due to Fischetti, Glover, and Lodi \
and uses a sequence of LPs to try and get an integer feasible solution. \
Some fine tuning is available by options passFeasibilityPump and pumpTune. " HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("fix!OnDj", "Try heuristic based on fixing variables with \
reduced costs greater than this",
               -COIN_DBL_MAX, COIN_DBL_MAX, CBC_PARAM_DBL_DJFIX, 1);
    p.setLonghelp(
        "If this is set integer variables with reduced costs greater than this will be fixed \
before branch and bound - use with extreme caution!");
    parameters.push_back(p);
  }
  {
    CbcParam p("flow!CoverCuts", "Whether to use Flow Cover cuts", "off",
               CBC_PARAM_STR_FLOWCUTS);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.setFakeKeyWord(3);
    p.setLonghelp(CUTS_LONGHELP " Reference: "
                                "https://github.com/coin-or/Cgl/wiki/"
                                "CglFlowCover"); // Can also enter testing
                                                 // values by plusnn (==ifmove)
    parameters.push_back(p);
  }
  {
    CbcParam p("force!Solution",
               "Whether to use given solution as crash for BAB", -1, 20000000,
               CBC_PARAM_INT_USESOLUTION);
    p.setIntValue(-1);
    p.setLonghelp(
        "-1 off.  If 1 then tries to branch to solution given by AMPL or priorities file. \
If 0 then just tries to set as best solution \
If >1 then also does that many nodes on fixed problem.");
    parameters.push_back(p);
  }
  {
    CbcParam p("fraction!forBAB", "Fraction in feasibility pump", 1.0e-5, 1.1,
               CBC_PARAM_DBL_SMALLBAB, 1);
    p.setDoubleValue(0.5);
    p.setLonghelp(
        "After a pass in the feasibility pump, variables which have not moved \
about are fixed and if the preprocessed model is smaller than this fraction of the original problem, \
a few nodes of branch and bound are done on the reduced problem.");
    parameters.push_back(p);
  }
  {
    CbcParam p("GMI!Cuts", "Whether to use alternative Gomory cuts", "off",
               CBC_PARAM_STR_GMICUTS);

    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("endonly");
    p.append("long");
    p.append("longroot");
    p.append("longifmove");
    p.append("forceLongOn");
    p.append("longendonly");
    p.setLonghelp(CUTS_LONGHELP " This version is by Giacomo Nannicini and may "
                                "be more robust than gomoryCuts.");
    parameters.push_back(p);
  }
  {
    CbcParam p("gomory!Cuts", "Whether to use Gomory cuts", "off",
               CBC_PARAM_STR_GOMORYCUTS);

    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.append("forceandglobal");
    p.append("forceLongOn");
    p.append("long");
    p.append("shorter");
    p.setLonghelp(
        "The original cuts - beware of imitations!  Having gone out of favor, \
they are now more fashionable as LP solvers are more robust and they interact well \
with other cuts.  They will almost always give cuts (although in this executable \
they are limited as to number of variables in cut).  However the cuts may be dense \
so it is worth experimenting (Long allows any length). " CUTS_LONGHELP
        " Reference: https://github.com/coin-or/Cgl/wiki/CglGomory");
    parameters.push_back(p);
  }
  {
    CbcParam p("greedy!Heuristic", "Whether to use a greedy heuristic", "off",
               CBC_PARAM_STR_GREEDY);

    p.append("on");
    p.append("both");
    p.append("before");
    // p.append("root");
    p.setLonghelp("This heuristic tries to obtain a feasible solution by just "
                  "fixing a percentage of variables and then try a small "
                  "branch and cut run. " HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("gsolu!tion", "Puts glpk solution to file",
               CBC_PARAM_ACTION_GMPL_SOLUTION);

    p.setLonghelp(
        "Will write a glpk solution file to the given file name.  It will use the default \
directory given by 'directory'.  A name of '$' will use the previous value for the \
name.  This is initialized to 'stdout' (this defaults to ordinary solution if stdout). \
If problem created from gmpl model - will do any reports.");
    parameters.push_back(p);
  }
  {
    CbcParam p("heur!isticsOnOff", "Switches most primal heuristics on or off",
               "off", CBC_PARAM_STR_HEURISTICSTRATEGY);
    p.append("on");
    p.setLonghelp(
        "This option can be used to switch on or off all heuristics that search for feasible solutions,\
      except for the local tree search, as it dramatically alters the search.\
      Then individual heuristics can be turned off or on.");
    parameters.push_back(p);
  }
  {
    CbcParam p("clique!Cuts", "Whether to use Clique cuts", "off",
               CBC_PARAM_STR_CLIQUECUTS);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.setLonghelp(
        "This switches on clique cuts (either at root or in entire tree). \
An improved version of the Bron-Kerbosch algorithm is used to separate cliques.");
    parameters.push_back(p);
  }
  {
    CbcParam p("cgraph",
               "Whether to use the conflict graph-based preprocessing and cut "
               "separation routines.",
               "on", CBC_PARAM_STR_USECGRAPH);
    p.append("off");
    p.append("clq");
    p.setLonghelp(
        "This switches the conflict graph-based preprocessing and cut separation routines \
(CglBKClique, CglOddWheel and CliqueStrengthening) on or off. Values:\
\n\toff: turns these routines off;\
\n\ton: turns these routines on;\
\n\tclq: turns these routines off and enables the cut separator of CglClique.");
    parameters.push_back(p);
  }
  {
    CbcParam p("bkpivot!ing",
               "Pivoting strategy used in Bron-Kerbosch algorithm", 0, 6,
               CBC_PARAM_INT_BKPIVOTINGSTRATEGY);
    p.setIntValue(3);
    parameters.push_back(p);
  }
  {
    CbcParam p(
        "bkmaxcalls",
        "Maximum number of recursive calls made by Bron-Kerbosch algorithm", 1,
        2147483647, CBC_PARAM_INT_BKMAXCALLS);
    p.setIntValue(1000);
    parameters.push_back(p);
  }
  {
    CbcParam p("bkclqext!method",
               "Strategy used to extend violated cliques found by BK Clique "
               "Cut Separation routine",
               0, 5, CBC_PARAM_INT_BKCLQEXTMETHOD);
    p.setLonghelp(
        "Sets the method used in the extension module of BK Clique Cut Separation routine: \
0=no extension; 1=random; 2=degree; 3=modified degree; 4=reduced cost(inversely proportional); 5=reduced cost(inversely proportional) + modified degree");
    p.setIntValue(4);
    parameters.push_back(p);
  }
  {
    CbcParam p("oddwheel!Cuts", "Whether to use odd wheel cuts", "off",
               CBC_PARAM_STR_ODDWHEELCUTS);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.setLonghelp("This switches on odd-wheel inequalities (either at root or "
                  "in entire tree).");
    parameters.push_back(p);
  }
  {
    CbcParam p(
        "oddwext!method",
        "Strategy used to search for wheel centers for the cuts found by \
Odd Wheel Cut Separation routine",
        0, 2, CBC_PARAM_INT_ODDWEXTMETHOD);
    p.setLonghelp(
        "Sets the method used in the extension module of Odd Wheel Cut Separation routine: \
0=no extension; 1=one variable; 2=clique");
    p.setIntValue(2);
    parameters.push_back(p);
  }
  {
    CbcParam p("clqstr!engthen",
               "Whether to perform Clique Strengthening preprocessing routine",
               "after", CBC_PARAM_STR_CLQSTRENGTHENING);
    p.setLonghelp(
        "Sets the method used in the Clique Strengthening Preprocessing routine:\
\n\toff: do not perform clique strengthening;\
\n\tbefore: perform clique strengthening before initialSolve;\
\n\tafter: perform clique strengthening after initialSolve.");
    p.append("off");
    p.append("before");
    parameters.push_back(p);
  }
  {
    CbcParam p("help", "Print out version, non-standard options and some help",
               CBC_PARAM_ACTION_HELP, 3);
    p.setLonghelp(
        "This prints out some help to get user started.  If you have printed this then \
you should be past that stage:-)");
    parameters.push_back(p);
  }
  {
    CbcParam p("hOp!tions", "Heuristic options", -COIN_INT_MAX, COIN_INT_MAX,
               CBC_PARAM_INT_HEUROPTIONS, 1);
    p.setIntValue(0);
    p.setLonghelp(
        "Value 1 stops heuristics immediately if the allowable gap has been reached. \
Other values are for the feasibility pump - \
2 says do exact number of passes given, \
4 only applies if an initial cutoff has been given and says relax after 50 passes, \
while 8 will adapt the cutoff rhs after the first solution if it looks as if the code is stalling.");
    parameters.push_back(p);
  }
  {
    CbcParam p("hot!StartMaxIts", "Maximum iterations on hot start", 0,
               COIN_INT_MAX, CBC_PARAM_INT_MAXHOTITS);
    parameters.push_back(p);
  }
  {
    CbcParam p("import", "Import model from mps file", CBC_PARAM_ACTION_IMPORT,
               3);
    p.setLonghelp(
        "This will read an MPS format file from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  If you have libgz then it can read compressed\
 files 'xxxxxxxx.gz' or 'xxxxxxxx.bz2'.  \
If 'keepnames' is off, then names are dropped -> Rnnnnnnn and Cnnnnnnn.");
    parameters.push_back(p);
  }
  {
    CbcParam p("inc!rement", "A valid solution must be at least this \
much better than last integer solution",
               -COIN_DBL_MAX, COIN_DBL_MAX, CBC_PARAM_DBL_INCREMENT);

    p.setLonghelp(
        "Whenever a solution is found the bound on the objective value for new solutions is set to the\
      objective function of the found solution (in a minimization sense) plus this.  If it is not set then CBC will try and work one out, e.g. if \
all objective coefficients are multiples of 0.01 and only integer variables have entries in \
the objective function, then the increment can be set to 0.01.  Be careful if setting this to a negative value!");

    parameters.push_back(p);
  }
  {
    CbcParam p("inf!easibilityWeight", "Each integer infeasibility is expected \
to cost this much",
               0.0, COIN_DBL_MAX, CBC_PARAM_DBL_INFEASIBILITYWEIGHT, 1);
    p.setLonghelp(
        "A primitive way of deciding which node to explore next.  Satisfying each integer infeasibility is \
expected to cost this much.");
    parameters.push_back(p);
  }
  {
    CbcParam p("integerT!olerance", "For a feasible solution \
no integer variable may be more than this away from an integer value",
               1.0e-20, 0.5, CBC_PARAM_DBL_INTEGERTOLERANCE);
    p.setLonghelp("Beware of setting this smaller than the primal feasibility "
                  "tolerance.");
    parameters.push_back(p);
  }
  {
    CbcParam p("knapsack!Cuts", "Whether to use Knapsack cuts", "off",
               CBC_PARAM_STR_KNAPSACKCUTS);

    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.append("forceandglobal");
    p.setLonghelp(
        CUTS_LONGHELP
        " Reference: https://github.com/coin-or/Cgl/wiki/CglKnapsackCover");
    parameters.push_back(p);
  }
  {
    CbcParam p("lagomory!Cuts", "Whether to use Lagrangean Gomory cuts", "off",
               CBC_PARAM_STR_LAGOMORYCUTS);
    p.append("endonlyroot");
    p.append("endcleanroot");
    p.append("root");
    p.append("endonly");
    p.append("endclean");
    p.append("endboth");
    p.append("onlyaswell");
    p.append("cleanaswell");
    p.append("bothaswell");
    p.append("onlyinstead");
    p.append("cleaninstead");
    p.append("bothinstead");
    p.append("onlyaswellroot");
    p.append("cleanaswellroot");
    p.append("bothaswellroot");
    p.setLonghelp(
        "This is a gross simplification of 'A Relax-and-Cut Framework for Gomory's Mixed-Integer Cuts' \
by Matteo Fischetti & Domenico Salvagnin.  This simplification \
just uses original constraints while modifying objective using other cuts. \
So you don't use messy constraints generated by Gomory etc. \
A variant is to allow non messy cuts e.g. clique cuts. \
So 'only' does this while 'clean' also allows integral valued cuts.  \
'End' is recommended and waits until other cuts have finished before it \
does a few passes. \
The length options for gomory cuts are used.");
    parameters.push_back(p);
  }
  {
    CbcParam p("latwomir!Cuts", "Whether to use Lagrangean TwoMir cuts", "off",
               CBC_PARAM_STR_LATWOMIRCUTS);

    p.append("endonlyroot");
    p.append("endcleanroot");
    p.append("endbothroot");
    p.append("endonly");
    p.append("endclean");
    p.append("endboth");
    p.append("onlyaswell");
    p.append("cleanaswell");
    p.append("bothaswell");
    p.append("onlyinstead");
    p.append("cleaninstead");
    p.append("bothinstead");
    p.setLonghelp("This is a Lagrangean relaxation for TwoMir cuts.  See \
  lagomoryCuts for description of options.");
    parameters.push_back(p);
  }
  {
    CbcParam p("lift!AndProjectCuts", "Whether to use Lift and Project cuts",
               "off", CBC_PARAM_STR_LANDPCUTS);

    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("iflongon");
    p.setLonghelp("These cuts may be expensive to compute. " CUTS_LONGHELP
                  " Reference: https://github.com/coin-or/Cgl/wiki/CglLandP");
    parameters.push_back(p);
  }
  {
    CbcParam p("local!TreeSearch",
               "Whether to use local tree search when a solution is found",
               "off", CBC_PARAM_STR_LOCALTREE);
    p.append("on");
    p.setLonghelp(
        "The heuristic is from Fischetti and Lodi and is not really a heuristic although it can be used as one \
(with limited functionality).  It is not switched on when heuristics are switched on.");
    parameters.push_back(p);
  }
  {
    CbcParam p("max!imize", "Set optimization direction to maximize",
               CBC_PARAM_ACTION_MAXIMIZE, 7);
    p.setLonghelp("The default is minimize - use 'maximize' for maximization.\n\
You can also use the parameters 'direction maximize'.");
    parameters.push_back(p);
  }
  {
    CbcParam p("maxN!odes", "Maximum number of nodes to do", -1, COIN_INT_MAX,
               CBC_PARAM_INT_MAXNODES);
    p.setLonghelp(
        "This is a repeatable way to limit search.  Normally using time is easier \
but then the results may not be repeatable.");
    parameters.push_back(p);
  }
  {
    CbcParam p("maxNI!FS",
               "Maximum number of nodes to be processed without improving the "
               "incumbent solution.",
               -1, COIN_INT_MAX, CBC_PARAM_INT_MAXNODESNOTIMPROVINGFS);
    p.setLonghelp(
        "This criterion specifies that when a feasible solution is available, the search should continue\
only if better feasible solutions were produced in the last nodes.");
    parameters.push_back(p);
  }
  {
    CbcParam p("maxSaved!Solutions", "Maximum number of solutions to save", 0,
               COIN_INT_MAX, CBC_PARAM_INT_MAXSAVEDSOLS);
    p.setLonghelp("Number of solutions to save.");
    parameters.push_back(p);
  }
  {
    CbcParam p("maxSo!lutions", "Maximum number of feasible solutions to get",
               1, COIN_INT_MAX, CBC_PARAM_INT_MAXSOLS);
    p.setLonghelp("You may want to stop after (say) two solutions or an hour.  \
This is checked every node in tree, so it is possible to get more solutions from heuristics.");
    parameters.push_back(p);
  }
  {
    CbcParam p("min!imize", "Set optimization direction to minimize",
               CBC_PARAM_ACTION_MINIMIZE, 7);
    p.setLonghelp("The default is minimize - use 'maximize' for maximization.\n\
This should only be necessary if you have previously set maximization \
You can also use the parameters 'direction minimize'.");
    parameters.push_back(p);
  }
  {
    CbcParam p("mipO!ptions", "Dubious options for mip", 0, COIN_INT_MAX,
               CBC_PARAM_INT_MIPOPTIONS, 0);
    p.setIntValue(1057);
    parameters.push_back(p);
  }
  {
    CbcParam p("more!MipOptions", "More dubious options for mip", -1,
               COIN_INT_MAX, CBC_PARAM_INT_MOREMIPOPTIONS, 0);
    parameters.push_back(p);
  }
  {
    CbcParam p("more2!MipOptions", "More more dubious options for mip", -1,
               COIN_INT_MAX, CBC_PARAM_INT_MOREMOREMIPOPTIONS, 0);
    p.setIntValue(0);
    parameters.push_back(p);
  }
  {
    CbcParam p("mixed!IntegerRoundingCuts",
               "Whether to use Mixed Integer Rounding cuts", "off",
               CBC_PARAM_STR_MIXEDCUTS);

    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.setLonghelp(
        CUTS_LONGHELP
        " Reference: "
        "https://github.com/coin-or/Cgl/wiki/CglMixedIntegerRounding2");
    parameters.push_back(p);
  }
  {
    CbcParam p("miplib", "Do some of miplib test set", CBC_PARAM_ACTION_MIPLIB,
               3, 1);
    parameters.push_back(p);
  }
  {
    CbcParam p("mips!tart", "reads an initial feasible solution from file",
               CBC_PARAM_ACTION_MIPSTART);
    p.setLonghelp("\
The MIPStart allows one to enter an initial integer feasible solution \
to CBC. Values of the main decision variables which are active (have \
non-zero values) in this solution are specified in a text  file. The \
text file format used is the same of the solutions saved by CBC, but \
not all fields are required to be filled. First line may contain the \
solution status and will be ignored, remaining lines contain column \
indexes, names and values as in this example:\n\
\n\
Stopped on iterations - objective value 57597.00000000\n\
      0  x(1,1,2,2)               1 \n\
      1  x(3,1,3,2)               1 \n\
      5  v(5,1)                   2 \n\
      33 x(8,1,5,2)               1 \n\
      ...\n\
\n\
Column indexes are also ignored since pre-processing can change them. \
There is no need to include values for continuous or integer auxiliary \
variables, since they can be computed based on main decision variables. \
Starting CBC with an integer feasible solution can dramatically improve \
its performance: several MIP heuristics (e.g. RINS) rely on having at \
least one feasible solution available and can start immediately if the \
user provides one. Feasibility Pump (FP) is a heuristic which tries to \
overcome the problem of taking too long to find feasible solution (or \
not finding at all), but it not always succeeds. If you provide one \
starting solution you will probably save some time by disabling FP. \
\n\n\
Knowledge specific to your problem can be considered to write an \
external module to quickly produce an initial feasible solution - some \
alternatives are the implementation of simple greedy heuristics or the \
solution (by CBC for example) of a simpler model created just to find \
a feasible solution. \
\n\n\
Silly options added.  If filename ends .low then integers not mentioned \
are set low - also .high, .lowcheap, .highcheap, .lowexpensive, .highexpensive \
where .lowexpensive sets costed ones to make expensive others low. Also if \
filename starts empty. then no file is read at all - just actions done. \
\n\n\
Question and suggestions regarding MIPStart can be directed to\n\
haroldo.santos@gmail.com. ");
    parameters.push_back(p);
  }
  {
    CbcParam p("moreT!une", "Yet more dubious ideas for feasibility pump", 0,
               100000000, CBC_PARAM_INT_FPUMPTUNE2, 0);

    p.setLonghelp("Yet more ideas for Feasibility Pump \n\
\t/100000 == 1 use box constraints and original obj in cleanup\n\
\t/1000 == 1 Pump will run twice if no solution found\n\
\t/1000 == 2 Pump will only run after root cuts if no solution found\n\
\t/1000 >10 as above but even if solution found\n\
\t/100 == 1,3.. exact 1.0 for objective values\n\
\t/100 == 2,3.. allow more iterations per pass\n\
\t n fix if value of variable same for last n iterations.");
    p.setIntValue(0);
    parameters.push_back(p);
  }
  {
    CbcParam p("multiple!RootPasses",
               "Do multiple root passes to collect cuts and solutions", 0,
               COIN_INT_MAX, CBC_PARAM_INT_MULTIPLEROOTS, 0);
    p.setIntValue(0);
    p.setLonghelp(
        "Solve (in parallel, if enabled) the root phase this number of times, \
      each with its own different seed, and collect all solutions and cuts generated. \
      The actual format is aabbcc where aa is the number of extra passes; \
      if bb is non zero, then it is number of threads to use (otherwise uses threads setting); \
      and cc is the number of times to do root phase. \
The solvers do not interact with each other.  However if extra passes are specified \
then cuts are collected and used in later passes - so there is interaction there. \
Some parts of this implementation have their origin in idea of \
Andrea Lodi, Matteo Fischetti, Michele Monaci, Domenico Salvagnin, and Andrea Tramontani.");
    parameters.push_back(p);
  }
  {
    CbcParam p("naive!Heuristics", "Whether to try some stupid heuristic",
               "off", CBC_PARAM_STR_NAIVE, 7, 1);

    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp("This is naive heuristics which, e.g., fix all integers with "
                  "costs to zero!. " HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("nextB!estSolution", "Prints next best saved solution to file",
               CBC_PARAM_ACTION_NEXTBESTSOLUTION);

    p.setLonghelp(
        "To write best solution, just use solution.  This prints next best (if exists) \
 and then deletes it. \
 This will write a primitive solution file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'stdout'.  The amount of output can be varied using printi!ngOptions or printMask.");
    parameters.push_back(p);
  }
  {
    CbcParam p("node!Strategy",
               "What strategy to use to select the next node from the branch "
               "and cut tree",
               "hybrid", CBC_PARAM_STR_NODESTRATEGY);
    p.append("fewest");
    p.append("depth");
    p.append("upfewest");
    p.append("downfewest");
    p.append("updepth");
    p.append("downdepth");
    p.setLonghelp(
        "Normally before a feasible solution is found, CBC will choose a node with fewest infeasibilities. \
  Alternatively, one may choose tree-depth as the criterion. This requires the minimal amount of memory, but may take a long time to find the best solution.\
  Additionally, one may specify whether up or down branches must \
be selected first (the up-down choice will carry on after a first solution has been bound). \
The choice 'hybrid' does breadth first on small depth nodes and then switches to 'fewest'.");
    parameters.push_back(p);
  }
  {
    CbcParam p("numberA!nalyze", "Number of analysis iterations", -COIN_INT_MAX,
               COIN_INT_MAX, CBC_PARAM_INT_NUMBERANALYZE, 0);
    p.setLonghelp(
        "This says how many iterations to spend at root node analyzing problem. \
This is a first try and will hopefully become more sophisticated.");
    parameters.push_back(p);
  }
#ifdef CBC_HAS_NAUTY
  {
    CbcParam p("Orbit!alBranching", "Whether to try orbital branching", "off",
               CBC_PARAM_STR_ORBITAL);
    p.append("slow!ish");
    p.append("strong");
    p.append("force");
    p.append("simple");
    p.append("on");
    p.append("more!printing");
    p.setLonghelp("This switches on Orbital branching. \
Value 'on' just adds orbital, 'strong' tries extra fixing in strong branching.");
    parameters.push_back(p);
  }
#endif
  {
    CbcParam p("PrepN!ames",
               "If column names will be kept in pre-processed model", "off",
               CBC_PARAM_STR_PREPROCNAMES);
    p.append("on");
    p.setLonghelp(
        "Normally the preprocessed model has column names replaced by new names C0000...\
Setting this option to on keeps original names in variables which still exist in the preprocessed problem");
    parameters.push_back(p);
  }

  {
    CbcParam p("output!Format", "Which output format to use", 1, 6,
               CBC_PARAM_INT_OUTPUTFORMAT);
    p.setLonghelp(
        "Normally export will be done using normal representation for numbers and two values\
 per line.  You may want to do just one per line (for grep or suchlike) and you may wish\
 to save with absolute accuracy using a coded version of the IEEE value. A value of 2 is normal.\
 otherwise odd values gives one value per line, even two.  Values 1,2 give normal format, 3,4\
 gives greater precision, while 5,6 give IEEE values.  When used for exporting a basis 1 does not save \
values, 2 saves values, 3 with greater accuracy and 4 in IEEE.");
    parameters.push_back(p);
  }
  {
    CbcParam p(
        "passC!uts",
        "Number of rounds that cut generators are applied in the root node",
        -COIN_INT_MAX, COIN_INT_MAX, CBC_PARAM_INT_CUTPASS);

    p.setIntValue(20);
    p.setLonghelp(
        "The default is to do 100 passes if the problem has less than 500 columns, 100 passes (but \
stop if the drop in the objective function value is small) if the problem has less than 5000 columns, and 20 passes otherwise. \
A negative value -n means that n passes are also applied if the objective does not drop.");
    parameters.push_back(p);
  }
  {
    CbcParam p("passF!easibilityPump",
               "How many passes to do in the Feasibility Pump heuristic", 0,
               10000, CBC_PARAM_INT_FPUMPITS);
    p.setIntValue(20);
    parameters.push_back(p);
  }
  {
    CbcParam p("passT!reeCuts",
               "Number of rounds that cut generators are applied in the tree",
               -COIN_INT_MAX, COIN_INT_MAX, CBC_PARAM_INT_CUTPASSINTREE);
    p.setIntValue(1);
    p.setLonghelp(
        "The default is to do one pass. A negative value -n means that n "
        "passes are also applied if the objective does not drop.");
    parameters.push_back(p);
  }
  {
    CbcParam p("pivotAndC!omplement",
               "Whether to try Pivot and Complement heuristic", "off",
               CBC_PARAM_STR_PIVOTANDCOMPLEMENT);

    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("pivotAndF!ix", "Whether to try Pivot and Fix heuristic", "off",
               CBC_PARAM_STR_PIVOTANDFIX);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("pO!ptions", "Dubious print options", 0, COIN_INT_MAX,
               CBC_PARAM_INT_PRINTOPTIONS, 1);
    p.setIntValue(0);
    p.setLonghelp("If this is > 0 then presolve will give more information and "
                  "branch and cut will give statistics");
    parameters.push_back(p);
  }
  {
    CbcParam p("preprocess", "Whether to use integer preprocessing", "off",
               CBC_PARAM_STR_PREPROCESS);

    p.append("on");
    p.append("save");
    p.append("equal");
    p.append("sos");
    p.append("trysos");
    p.append("equalall");
    p.append("strategy");
    p.append("aggregate");
    p.append("forcesos");
    p.append("stop!aftersaving");
    p.setLonghelp(
        "This tries to reduce size of model in a similar way to presolve and \
it also tries to strengthen the model - this can be very useful and is worth trying. \
 Value 'save' saves the presolved problem to a file presolved.mps.\
 Value 'equal' will turn inequality-cliques into equalities.\
 Value 'sos' lets CBC search for rows with upper bound 1 and where all nonzero coefficients are 1 and creates special ordered sets if the sets are not overlapping and all integer variables (except for at most one) are in the sets.\
 Value 'trysos' is same as 'sos', but allows any number of integer variables outside of sets.\
 Value 'equalall' lets CBC turn all valid inequalities into equalities by adding integer slack variables."); // Value 'strategy' is as on but uses CbcStrategy.
    parameters.push_back(p);
  }
  {
    CbcParam p("printi!ngOptions", "Print options", "normal",
               CBC_PARAM_STR_INTPRINT, 3);
    p.append("integer");
    p.append("special");
    p.append("rows");
    p.append("all");
    p.append("csv");
    p.append("bound!ranging");
    p.append("rhs!ranging");
    p.append("objective!ranging");
    p.append("stats");
    p.append("boundsint");
    p.append("boundsall");
    p.append("fixint");
    p.append("fixall");
    p.setLonghelp(
        "This changes the amount and format of printing a solution:\nnormal - nonzero column variables \n\
integer - nonzero integer column variables\n\
special - in format suitable for OsiRowCutDebugger\n\
rows - nonzero column variables and row activities\n\
all - all column variables and row activities.\n\
\nFor non-integer problems 'integer' and 'special' act like 'normal'.  \
Also see printMask for controlling output.");
    parameters.push_back(p);
  }
  {
    CbcParam p("printM!ask", "Control printing of solution on a  mask",
               CBC_PARAM_ACTION_PRINTMASK, 3);

    p.setLonghelp(
        "If set then only those names which match mask are printed in a solution. \
'?' matches any character and '*' matches any set of characters. \
 The default is '' i.e. unset so all variables are printed. \
This is only active if model has names.");
    parameters.push_back(p);
  }

  {
    CbcParam p("prio!rityIn", "Import priorities etc from file",
               CBC_PARAM_ACTION_PRIORITYIN, 3);
    p.setLonghelp(
        "This will read a file with priorities from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  This can not read from compressed files. \
File is in csv format with allowed headings - name, number, priority, direction, up, down, solution.  Exactly one of\
 name and number must be given.");
    parameters.push_back(p);
  }

  {
    CbcParam p("probing!Cuts", "Whether to use Probing cuts", "off",
               CBC_PARAM_STR_PROBINGCUTS);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.append("forceonglobal");
    p.append("forceOnBut");
    p.append("forceOnStrong");
    p.append("forceOnButStrong");
    p.append("strongRoot");
    p.setLonghelp(
        CUTS_LONGHELP
        " Value 'forceOnBut' turns on probing and forces CBC to do probing at every node, but does only probing, not strengthening etc. \
    Value 'strong' forces CBC to strongly do probing at every node, that is, also when CBC would usually turn it off because it hasn't found something. \
    Value 'forceonbutstrong' is like 'forceonstrong', but does only probing (column fixing) and turns off row strengthening, so the matrix will not change inside the branch and bound. \
    Reference: https://github.com/coin-or/Cgl/wiki/CglProbing");
    parameters.push_back(p);
  }
  {
    CbcParam p("proximity!Search", "Whether to do proximity search heuristic",
               "off", CBC_PARAM_STR_PROXIMITY);

    p.append("on");
    p.append("both");
    p.append("before");
    p.append("10");
    p.append("100");
    p.append("300");
    // but allow numbers after this (returning 1)
    p.setFakeKeyWord(1);
    p.setLonghelp(
        "This heuristic looks for a solution close to the incumbent solution (Fischetti and Monaci, 2012). \
The idea is to define a sub-MIP without additional constraints but with a modified objective function intended to attract the search \
in the proximity of the incumbent. \
The approach works well for 0-1 MIPs whose solution landscape is not too irregular (meaning the there is reasonable probability of \
finding an improved solution by flipping a small number of binary variables), in particular when it is applied to the first heuristic solutions \
found at the root node. " HEURISTICS_LONGHELP); // Can also set different
                                                // maxNode settings by plusnnnn
                                                // (and are 'on'(on==30)).
    parameters.push_back(p);
  }
  {
    CbcParam p("pumpC!utoff", "Fake cutoff for use in feasibility pump",
               -COIN_DBL_MAX, COIN_DBL_MAX, CBC_PARAM_DBL_FAKECUTOFF);
    p.setDoubleValue(0.0);
    p.setLonghelp(
        "A value of 0.0 means off. Otherwise, add a constraint forcing objective below this value\
 in feasibility pump");
    parameters.push_back(p);
  }
  {
    CbcParam p("pumpI!ncrement", "Fake increment for use in feasibility pump",
               -COIN_DBL_MAX, COIN_DBL_MAX, CBC_PARAM_DBL_FAKEINCREMENT, 1);
    p.setDoubleValue(0.0);
    p.setLonghelp(
        "A value of 0.0 means off. Otherwise use as absolute increment to cutoff \
when solution found in feasibility pump");
    parameters.push_back(p);
  }
  {
    CbcParam p("pumpT!une", "Dubious ideas for feasibility pump", 0, 100000000,
               CBC_PARAM_INT_FPUMPTUNE);
    p.setIntValue(1003);
    p.setLonghelp("This fine tunes Feasibility Pump \n\
\t>=10000000 use as objective weight switch\n\
\t>=1000000 use as accumulate switch\n\
\t>=1000 use index+1 as number of large loops\n\
\t==100 use objvalue +0.05*fabs(objvalue) as cutoff OR fakeCutoff if set\n\
\t%100 == 10,20 affects how each solve is done\n\
\t1 == fix ints at bounds, 2 fix all integral ints, 3 and continuous at bounds. \
If accumulate is on then after a major pass, variables which have not moved \
are fixed and a small branch and bound is tried.");
    p.setIntValue(0);
    parameters.push_back(p);
  }
  {
    CbcParam p("randomC!bcSeed", "Random seed for Cbc", -1, COIN_INT_MAX,
               CBC_PARAM_INT_RANDOMSEED);

    p.setLonghelp("Allows initialization of the random seed for pseudo-random "
                  "numbers used in heuristics such as the Feasibility Pump to "
                  "decide whether to round up or down. "
                  "The special value of 0 lets Cbc use the time of the day for "
                  "the initial seed.");
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcParam p("randomi!zedRounding",
               "Whether to try randomized rounding heuristic", "off",
               CBC_PARAM_STR_RANDROUND);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("ratio!Gap", "Stop when gap between best possible and \
best known is less than this fraction of larger of two",
               0.0, COIN_DBL_MAX, CBC_PARAM_DBL_GAPRATIO);
    p.setDoubleValue(1e-4);
    p.setLonghelp(
        "If the gap between the best known solution and the best possible solution is less than this fraction \
of the objective value at the root node then the search will terminate.  See 'allowableGap' for a \
way of using absolute value rather than fraction.");
    parameters.push_back(p);
  }
  {
    CbcParam p("reduce!AndSplitCuts", "Whether to use Reduce-and-Split cuts",
               "off", CBC_PARAM_STR_REDSPLITCUTS);

    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.setLonghelp(
        "These cuts may be expensive to generate. " CUTS_LONGHELP
        " Reference: https://github.com/coin-or/Cgl/wiki/CglRedSplit");
    parameters.push_back(p);
  }
  {
    CbcParam p("reduce2!AndSplitCuts",
               "Whether to use Reduce-and-Split cuts - style 2", "off",
               CBC_PARAM_STR_REDSPLIT2CUTS);
    p.append("on");
    p.append("root");
    p.append("longOn");
    p.append("longRoot");
    p.setLonghelp(
        "This switches on reduce and split  cuts (either at root or in entire tree). \
This version is by Giacomo Nannicini based on Francois Margot's version. \
Standard setting only uses rows in tableau <= 256, long uses all. \
These cuts may be expensive to generate. \
See option cuts for more information on the possible values.");
    parameters.push_back(p);
  }
  {
    CbcParam p("residual!CapacityCuts", "Whether to use Residual Capacity cuts",
               "off", CBC_PARAM_STR_RESIDCUTS);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.setLonghelp(
        CUTS_LONGHELP
        " Reference: https://github.com/coin-or/Cgl/wiki/CglResidualCapacity");

    parameters.push_back(p);
  }
  {
    CbcParam p("restore!Model", "Restore model from binary file",
               CBC_PARAM_ACTION_RESTORE, 7, 1);
    p.setLonghelp(
        "This reads data save by saveModel from the given file.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.prob'.");

    parameters.push_back(p);
  }
  {
    CbcParam p("reverse", "Reverses sign of objective",
               CBC_PARAM_ACTION_REVERSE, 7, 0);
    p.setLonghelp("Useful for testing if maximization works correctly");
    parameters.push_back(p);
  }
  {
    CbcParam p("Rens", "Whether to try Relaxation Enforced Neighborhood Search",
               "off", CBC_PARAM_STR_RENS);
    p.append("on");
    p.append("both");
    p.append("before");
    p.append("200");
    p.append("1000");
    p.append("10000");
    p.append("dj");
    p.append("djbefore");
    p.append("usesolution");
    p.setLonghelp(HEURISTICS_LONGHELP " Value 'on' just does 50 nodes. 200, "
                                      "1000, and 10000 does that many nodes.");
    parameters.push_back(p);
  }
  {
    CbcParam p("Rins", "Whether to try Relaxed Induced Neighborhood Search",
               "off", CBC_PARAM_STR_RINS);
    p.append("on");
    p.append("both");
    p.append("before");
    p.append("often");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("round!ingHeuristic",
               "Whether to use simple (but effective) Rounding heuristic",
               "off", CBC_PARAM_STR_ROUNDING);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("saveM!odel", "Save model to binary file", CBC_PARAM_ACTION_SAVE,
               7, 1);
    p.setLonghelp(
        "This will save the problem to the given file name for future use\
 by restoreModel.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.prob'.");
    parameters.push_back(p);
  }
  {
    CbcParam p("saveS!olution", "saves solution to file",
               CBC_PARAM_ACTION_SAVESOL);

    p.setLonghelp(
        "This will write a binary solution file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'solution.file'.  To read the file use fread(int) twice to pick up number of rows \
and columns, then fread(double) to pick up objective value, then pick up row activities, row duals, column \
activities and reduced costs - see bottom of CbcParam.cpp for code that reads or writes file. \
If name contains '_fix_read_' then does not write but reads and will fix all variables");
    parameters.push_back(p);
  }
  {
    CbcParam p("sec!onds", "maximum seconds", -1.0, COIN_DBL_MAX,
               CBC_PARAM_DBL_TIMELIMIT_BAB);
    // Meaning 0 - start at very beginning
    // 1 start at beginning of preprocessing
    // 2 start at beginning of branch and bound
#ifndef CBC_USE_INITIAL_TIME
#define CBC_USE_INITIAL_TIME 1
#endif
#if CBC_USE_INITIAL_TIME == 0
    p.setLonghelp("After this many seconds in Branch and Bound coin solver "
                  "will act as if maximum nodes had been reached (time in "
                  "initial solve and preprocessing is included).");
#elif CBC_USE_INITIAL_TIME == 1
    p.setLonghelp("After this many seconds in Branch and Bound coin solver "
                  "will act as if maximum nodes had been reached (time in "
                  "initial solve not included).");
#else
    p.setLonghelp("After this many seconds in Branch and Bound coin solver "
                  "will act as if maximum nodes had been reached (time in "
                  "initial solve and preprocessing not included).");
#endif
    parameters.push_back(p);
  }
  {
    CbcParam p("secni!fs",
               "maximum seconds without improving the incumbent solution", -1.0,
               COIN_DBL_MAX, CBC_PARAM_DBL_MAXSECONDSNIFS);
    p.setLonghelp(
        "With this stopping criterion, after a feasible solution is found, the search should continue only if the incumbent solution was updated recently, \
the tolerance is specified here. A discussion on why this criterion can be useful is included here: \
https://yetanothermathprogrammingconsultant.blogspot.com/2019/11/mip-solver-stopping-criteria.html .");
    parameters.push_back(p);
  }
  {
    CbcParam p("sleep", "for debug", CBC_PARAM_ACTION_DUMMY, 7, 0);

    p.setLonghelp("If passed to solver fom ampl, then ampl will wait so that "
                  "you can copy .nl file for debug.");
    parameters.push_back(p);
  }
  {
    CbcParam p("slow!cutpasses",
               "Maximum number of rounds for slower cut generators", -1,
               COIN_INT_MAX, CBC_PARAM_INT_MAXSLOWCUTS);
    p.setLonghelp(
        "Some cut generators are fairly slow - this limits the number of times they are tried.\
      The cut generators identified as 'may be slow' at present are Lift and project cuts and both versions of Reduce and Split cuts.");
    p.setIntValue(10);
    parameters.push_back(p);
  }
  {
    CbcParam p("solu!tion", "Prints solution to file",
               CBC_PARAM_ACTION_SOLUTION);
    p.setLonghelp(
        "This will write a primitive solution file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'stdout'.  The amount of output can be varied using printi!ngOptions or printMask.");
    parameters.push_back(p);
  }
  {
    CbcParam p("solv!e", "Solve problem", CBC_PARAM_ACTION_BAB);
    p.setLonghelp(
        "If there are no integer variables then this just solves LP.  If there are integer variables \
this does branch and cut.");
    parameters.push_back(p);
  }
  {
    CbcParam p("sosO!ptions", "Whether to use SOS from AMPL", "off",
               CBC_PARAM_STR_SOS);
    p.append("on");
    p.setCurrentOption("on");
    p.setLonghelp(
        "Normally if AMPL says there are SOS variables they should be used, but sometime sthey should\
 be turned off - this does so.");
    parameters.push_back(p);
  }
  {
    CbcParam p("slog!Level", "Level of detail in main solver output", -1, 63,
               CBC_PARAM_INT_SOLVERLOGLEVEL);
    p.setLonghelp(
        "If 0 then there should be no output in normal circumstances.  1 is probably the best\
 value for most uses, while 2 and 3 give more information.");
    parameters.push_back(p);
  }
  {
    CbcParam p("lplog!Level", "Level of detail in LP solver output", -1, 63,
               CBC_PARAM_INT_LPLOGLEVEL);
    p.setLonghelp(
        "If 0 then there should be no output in normal circumstances.  1 is probably the best\
 value for most uses, while 2 and 3 give more information.  This parameter is only for controlling the output of the LP solver");
    parameters.push_back(p);
  }
  {
    // Due to James Howey
    CbcParam p("sosP!rioritize", "How to deal with SOS priorities", "off",
               CBC_PARAM_STR_SOSPRIORITIZE);
    p.append("high");
    p.append("low");
    p.append("orderhigh");
    p.append("orderlow");
    p.setLonghelp(
        "This sets priorities for SOS.  Values 'high' and 'low' just set a priority \
    relative to the for integer variables.  Value 'orderhigh' gives first highest priority to the first SOS and integer variables \
    a low priority.  Value 'orderlow' gives integer variables a high priority then SOS in order.");
    parameters.push_back(p);
  }
  {
    CbcParam p("timeM!ode", "Whether to use CPU or elapsed time", "cpu",
               CBC_PARAM_STR_TIME_MODE);
    p.append("elapsed");
    p.setLonghelp(
        "cpu uses CPU time for stopping, while elapsed uses elapsed time. \
(On Windows, elapsed time is always used).");
    parameters.push_back(p);
  }
  {
    CbcParam p("stat!istics", "Print some statistics",
               CBC_PARAM_ACTION_STATISTICS);
    p.setLonghelp("This command prints some statistics for the current model.\
 If log level >1 then more is printed.\
 These are for presolved model if presolve on (and unscaled).");
    parameters.push_back(p);
  }
  {
    CbcParam p("strat!egy", "Switches on groups of features", 0, 2,
               CBC_PARAM_INT_STRATEGY);
    p.setLonghelp("This turns on newer features. \
Use 0 for easy problems, 1 is default, 2 is aggressive. \
1 uses Gomory cuts with a tolerance of 0.01 at the root node, \
does a possible restart after 100 nodes if many variables could be fixed, \
activates a diving and RINS heuristic, and makes the feasibility pump \
more aggressive."); // This does not apply to unit tests (where 'experiment' may
                    // have similar effects
    p.setIntValue(1);
    parameters.push_back(p);
  }
#ifdef CBC_KEEP_DEPRECATED
  {
    CbcParam p("strengthen", "Create strengthened problem",
               CBC_PARAM_ACTION_STRENGTHEN, 3);
    p.setLonghelp(
        "This creates a new problem by applying the root node cuts.  All tight constraints \
will be in resulting problem");
    parameters.push_back(p);
  }
#endif
  {
    CbcParam p("strong!Branching",
               "Number of variables to look at in strong branching", 0,
               COIN_INT_MAX, CBC_PARAM_INT_STRONGBRANCHING);
    p.setIntValue(20);
    p.setLonghelp(
        "In order to decide which variable to branch on, the code will choose up to this number \
of unsatisfied variables to try minimal up and down branches on.  Then the most effective one is chosen. \
If a variable is branched on many times then the previous average up and down costs may be used - \
see also option trustPseudoCosts.");
    parameters.push_back(p);
  }
  {
    CbcParam p("testO!si", "Test OsiObject stuff", -1, COIN_INT_MAX,
               CBC_PARAM_INT_TESTOSI, 0);
    parameters.push_back(p);
  }
#ifdef CBC_THREAD
  {
    CbcParam p("thread!s", "Number of threads to try and use", -100, 100000,
               CBC_PARAM_INT_THREADS, 1);
    p.setIntValue(0);
    p.setLonghelp(
        "To use multiple threads, set threads to number wanted.  It may be better \
to use one or two more than number of cpus available.  If 100+n then n threads and \
search is repeatable (maybe be somewhat slower), \
if 200+n use threads for root cuts, 400+n threads used in sub-trees.");
    parameters.push_back(p);
  }
#endif
  {
    CbcParam p("tighten!Factor", "Tighten bounds using this times largest \
activity at continuous solution",
               1.0e-3, COIN_DBL_MAX, CBC_PARAM_DBL_TIGHTENFACTOR, 0);
    p.setLonghelp("This sleazy trick can help on some problems.");
    parameters.push_back(p);
  }
  {
    CbcParam p("trust!PseudoCosts",
               "Number of branches before we trust pseudocosts", -3,
               COIN_INT_MAX, CBC_PARAM_INT_NUMBERBEFORE);
    p.setLonghelp(
        "Using strong branching computes pseudo-costs.  This parameter determines after how many branches for a variable we just \
trust the pseudo costs and do not do any more strong branching.");
    p.setIntValue(10);
    parameters.push_back(p);
  }
  {
    CbcParam p("tune!PreProcess", "Dubious tuning parameters for preprocessing",
               0, COIN_INT_MAX, CBC_PARAM_INT_PROCESSTUNE, 1);
    p.setLonghelp(
        "Format aabbcccc - \n If aa then this is number of major passes (i.e. with presolve) \n \
If bb and bb>0 then this is number of minor passes (if unset or 0 then 10) \n \
cccc is bit set \n 0 - 1 Heavy probing \n 1 - 2 Make variables integer if possible (if obj value)\n \
2 - 4 As above but even if zero objective value\n \
7 - 128 Try and create cliques\n 8 - 256 If all +1 try hard for dominated rows\n \
9 - 512 Even heavier probing \n \
10 - 1024 Use a larger feasibility tolerance in presolve\n \
11 - 2048 Try probing before creating cliques\n \
12 - 4096 Switch off duplicate column checking for integers \n \
13 - 8192 Allow scaled duplicate column checking \n \n \
     Now aa 99 has special meaning i.e. just one simple presolve.");
    parameters.push_back(p);
  }
  {
    CbcParam p("two!MirCuts",
               "Whether to use Two phase Mixed Integer Rounding cuts", "off",
               CBC_PARAM_STR_TWOMIRCUTS);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.append("forceandglobal");
    p.append("forceLongOn");
    p.setLonghelp(CUTS_LONGHELP
                  " Reference: https://github.com/coin-or/Cgl/wiki/CglTwomir");
    parameters.push_back(p);
  }
  {
    CbcParam p("unitTest", "Do unit test", CBC_PARAM_ACTION_UNITTEST, 3, 1);
    p.setLonghelp("This exercises the unit test for clp");
    parameters.push_back(p);
  }
  {
    CbcParam p("userCbc", "Hand coded Cbc stuff", CBC_PARAM_ACTION_USERCBC, 0,
               0);
    p.setLonghelp(
        "There are times e.g. when using AMPL interface when you may wish to do something unusual.  \
Look for USERCBC in main driver and modify sample code. \
It is possible you can get same effect by using example driver4.cpp.");
    parameters.push_back(p);
  }
  {
    CbcParam p("Vnd!VariableNeighborhoodSearch",
               "Whether to try Variable Neighborhood Search", "off",
               CBC_PARAM_STR_VND);
    p.append("on");
    p.append("both");
    p.append("before");
    p.append("intree");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcParam p("verbose", "Switches on longer help on single ?", 0, 31,
               CBC_PARAM_INT_VERBOSE, 0);
    p.setLonghelp("Set to 1 to get short help with ? list, 2 to get long help, "
                  "3 for both.  (add 4 to just get ampl ones).");
    p.setIntValue(0);
    parameters.push_back(p);
  }
  {
    CbcParam p("vub!heuristic", "Type of VUB heuristic", -2, 20,
               CBC_PARAM_INT_VUBTRY, 0);
    p.setLonghelp("This heuristic tries and fix some integer variables.");
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcParam p("zero!HalfCuts", "Whether to use zero half cuts", "off",
               CBC_PARAM_STR_ZEROHALFCUTS);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.setLonghelp(CUTS_LONGHELP
                  " This implementation was written by Alberto Caprara.");
    parameters.push_back(p);
  }
}

// Given a parameter type - returns its number in list
int whichCbcParam(const CbcParameterType &name,
                  const std::vector<CbcParam> &parameters) {
  for (int i = 0; i < (int)parameters.size(); i++) {
    if (parameters[i].type() == name)
      return i;
  }
  return std::numeric_limits<int>::max(); // should not arrive here
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
