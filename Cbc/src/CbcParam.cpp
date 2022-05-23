/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include <string>
#include <iostream>
#include <cassert>

#include "CbcParam.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CbcParam::CbcParam()
  : type_(CBC_PARAM_NOTUSED_INVALID)
  , lowerDoubleValue_(0.0)
  , upperDoubleValue_(0.0)
  , lowerIntValue_(0)
  , upperIntValue_(0)
  , lengthName_(0)
  , lengthMatch_(0)
  , definedKeyWords_()
  , name_()
  , shortHelp_()
  , longHelp_()
  , action_(CBC_PARAM_NOTUSED_INVALID)
  , currentKeyWord_(-1)
  , display_(false)
  , intValue_(0)
  , doubleValue_(0)
  , indexNumber_(0)
{
}
// Other constructors
CbcParam::CbcParam(std::string name, std::string help,
  double lower, double upper, CbcParameterType type,
  bool display)
  : type_(type)
  , lowerIntValue_(0)
  , upperIntValue_(0)
  , definedKeyWords_()
  , name_(name)
  , shortHelp_(help)
  , longHelp_()
  , action_(type)
  , currentKeyWord_(-1)
  , display_(false)
  , intValue_(0)
  , doubleValue_(0)
  , indexNumber_(0)
{
  lowerDoubleValue_ = lower;
  upperDoubleValue_ = upper;
  gutsOfConstructor();
}
CbcParam::CbcParam(std::string name, std::string help,
  int lower, int upper, CbcParameterType type,
  bool display)
  : type_(type)
  , lowerDoubleValue_(0.0)
  , upperDoubleValue_(0.0)
  , definedKeyWords_()
  , name_(name)
  , shortHelp_(help)
  , longHelp_()
  , action_(type)
  , currentKeyWord_(-1)
  , display_(false)
  , intValue_(0)
  , doubleValue_(0)
  , indexNumber_(0)
{
  gutsOfConstructor();
  lowerIntValue_ = lower;
  upperIntValue_ = upper;
}
// Other strings will be added by append
CbcParam::CbcParam(std::string name, std::string help,
  std::string firstValue,
  CbcParameterType type, int defaultIndex,
  bool display)
  : type_(type)
  , lowerDoubleValue_(0.0)
  , upperDoubleValue_(0.0)
  , lowerIntValue_(0)
  , upperIntValue_(0)
  , definedKeyWords_()
  , name_(name)
  , shortHelp_(help)
  , longHelp_()
  , action_(type)
  , currentKeyWord_(0)
  , display_(false)
  , intValue_(0)
  , doubleValue_(0)
  , indexNumber_(0)
{
  gutsOfConstructor();
  definedKeyWords_.push_back(firstValue);
}
// Action
CbcParam::CbcParam(std::string name, std::string help,
  CbcParameterType type, int indexNumber,
  bool display)
  : type_(type)
  , lowerDoubleValue_(0.0)
  , upperDoubleValue_(0.0)
  , lowerIntValue_(0)
  , upperIntValue_(0)
  , definedKeyWords_()
  , name_(name)
  , shortHelp_(help)
  , longHelp_()
  , action_(type)
  , currentKeyWord_(-1)
  , display_(false)
  , intValue_(0)
  , doubleValue_(0)
  , indexNumber_(0)
{
  gutsOfConstructor();
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CbcParam::CbcParam(const CbcParam &rhs)
{
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
  indexNumber_ = rhs.indexNumber_;
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CbcParam::~CbcParam()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CbcParam &
CbcParam::operator=(const CbcParam &rhs)
{
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
    indexNumber_ = rhs.indexNumber_;
  }
  return *this;
}
void CbcParam::gutsOfConstructor()
{
  std::string::size_type shriekPos = name_.find('!');
  lengthName_ = name_.length();
  if (shriekPos == std::string::npos) {
    //does not contain '!'
    lengthMatch_ = lengthName_;
  } else {
    lengthMatch_ = shriekPos;
    name_ = name_.substr(0, shriekPos) + name_.substr(shriekPos + 1);
    lengthName_--;
  }
}
// Insert string (only valid for keywords)
void CbcParam::append(std::string keyWord)
{
  definedKeyWords_.push_back(keyWord);
}

int CbcParam::matches(std::string input) const
{
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
// Returns name which could match
std::string
CbcParam::matchName() const
{
  if (lengthMatch_ == lengthName_)
    return name_;
  else
    return name_.substr(0, lengthMatch_) + "(" + name_.substr(lengthMatch_) + ")";
}

// Returns parameter option which matches (-1 if none)
int CbcParam::parameterOption(std::string check) const
{
  int numberItems = definedKeyWords_.size();
  if (!numberItems) {
    return -1;
  } else {
    int whichItem = 0;
    unsigned int it;
    for (it = 0; it < definedKeyWords_.size(); it++) {
      std::string thisOne = definedKeyWords_[it];
      std::string::size_type shriekPos = thisOne.find('!');
      unsigned int length1 = thisOne.length();
      unsigned int length2 = length1;
      if (shriekPos != std::string::npos) {
        //contains '!'
        length2 = shriekPos;
        thisOne = thisOne.substr(0, shriekPos) + thisOne.substr(shriekPos + 1);
        length1 = thisOne.length();
      }
      if (check.length() <= length1) {
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
    if (whichItem < numberItems)
      return whichItem;
    else
      return -1;
  }
}
// Prints parameter options
void CbcParam::printOptions() const
{
  std::cout << "Possible options for " << name_ << " are:" << std::endl;
  unsigned int it;
  for (it = 0; it < definedKeyWords_.size(); it++) {
    std::string thisOne = definedKeyWords_[it];
    std::string::size_type shriekPos = thisOne.find('!');
    if (shriekPos != std::string::npos) {
      //contains '!'
      thisOne = thisOne.substr(0, shriekPos) + "(" + thisOne.substr(shriekPos + 1) + ")";
    }
    std::cout << thisOne << std::endl;
  }
}
int CbcParam::setDoubleParameter(OsiSolverInterface *model, double value) const
{
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is " << lowerDoubleValue_ << " to " << upperDoubleValue_ << std::endl;
    return 1;
  } else {
    double oldValue;
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
      oldValue = 0.0; // to avoid compiler message
      abort();
    }
    std::cout << name_ << " was changed from " << oldValue << " to "
              << value << std::endl;
    return 0;
  }
}
int CbcParam::checkDoubleParameter(double value) const
{
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is " << lowerDoubleValue_ << " to " << upperDoubleValue_ << std::endl;
    return 1;
  } else {
    return 0;
  }
}
double
CbcParam::doubleParameter(OsiSolverInterface *model) const
{
  double value;
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
    abort();
  }
  (void)getDblParamRetValue;
  return value;
}
int CbcParam::setIntParameter(OsiSolverInterface *model, int value) const
{
  if (value < lowerIntValue_ || value > upperIntValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is " << lowerIntValue_ << " to " << upperIntValue_ << std::endl;
    return 1;
  } else {
    int oldValue;
    switch (type_) {
    case CLP_PARAM_INT_LOGLEVEL:
      oldValue = model->messageHandler()->logLevel();
      model->messageHandler()->setLogLevel(value);
      break;
    default:
      oldValue = 0; // to avoid compiler message
      abort();
    }
    std::cout << name_ << " was changed from " << oldValue << " to "
              << value << std::endl;
    return 0;
  }
}
int CbcParam::intParameter(OsiSolverInterface *model) const
{
  int value = 0;
  switch (type_) {
  case CLP_PARAM_INT_LOGLEVEL:
    //value=model->logLevel();
    break;
  default:
    abort();
  }
  return value;
}
int CbcParam::setDoubleParameter(CbcModel &model, double value) const
{
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is " << lowerDoubleValue_ << " to " << upperDoubleValue_ << std::endl;
    return 1;
  } else {
    double oldValue;
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
      break;
    case CBC_PARAM_DBL_ALLOWABLEGAP:
      oldValue = model.getDblParam(CbcModel::CbcAllowableGap);
      model.setDblParam(CbcModel::CbcAllowableGap, value);
      break;
    case CLP_PARAM_DBL_TIMELIMIT: {
      oldValue = model.getDblParam(CbcModel::CbcMaximumSeconds);
      model.setDblParam(CbcModel::CbcMaximumSeconds, value);
      break;
    }
    case CLP_PARAM_DBL_DUALTOLERANCE:
    case CLP_PARAM_DBL_PRIMALTOLERANCE:
      setDoubleParameter(model.solver(), value);
      return 0; // to avoid message
    default:
      oldValue = 0.0; // to avoid compiler message
      abort();
    }
    std::cout << name_ << " was changed from " << oldValue << " to "
              << value << std::endl;
    return 0;
  }
}
double
CbcParam::doubleParameter(CbcModel &model) const
{
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
  case CBC_PARAM_DBL_ALLOWABLEGAP:
    value = model.getDblParam(CbcModel::CbcAllowableGap);
    break;
  case CLP_PARAM_DBL_TIMELIMIT: {
    value = model.getDblParam(CbcModel::CbcMaximumSeconds);
    break;
  }
  case CLP_PARAM_DBL_DUALTOLERANCE:
  case CLP_PARAM_DBL_PRIMALTOLERANCE:
    value = doubleParameter(model.solver());
    break;
  default:
    abort();
  }
  return value;
}
int CbcParam::setIntParameter(CbcModel &model, int value) const
{
  if (value < lowerIntValue_ || value > upperIntValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is " << lowerIntValue_ << " to " << upperIntValue_ << std::endl;
    return 1;
  } else {
    int oldValue;
    switch (type_) {
    case CLP_PARAM_INT_LOGLEVEL:
      oldValue = model.messageHandler()->logLevel();
      model.messageHandler()->setLogLevel(value);
      break;
    case CLP_PARAM_INT_SOLVERLOGLEVEL:
      oldValue = model.solver()->messageHandler()->logLevel();
      model.solver()->messageHandler()->setLogLevel(value);
      break;
    case CBC_PARAM_INT_MAXNODES:
      oldValue = model.getIntParam(CbcModel::CbcMaxNumNode);
      model.setIntParam(CbcModel::CbcMaxNumNode, value);
      break;
    case CBC_PARAM_INT_STRONGBRANCHING:
      oldValue = model.numberStrong();
      model.setNumberStrong(value);
      break;
    default:
      oldValue = 0; // to avoid compiler message
      abort();
    }
    std::cout << name_ << " was changed from " << oldValue << " to "
              << value << std::endl;
    return 0;
  }
}
int CbcParam::intParameter(CbcModel &model) const
{
  int value;
  switch (type_) {
  case CLP_PARAM_INT_LOGLEVEL:
    value = model.messageHandler()->logLevel();
    break;
  case CLP_PARAM_INT_SOLVERLOGLEVEL:
    value = model.solver()->messageHandler()->logLevel();
    break;
  case CBC_PARAM_INT_MAXNODES:
    value = model.getIntParam(CbcModel::CbcMaxNumNode);
    break;
  case CBC_PARAM_INT_STRONGBRANCHING:
    value = model.numberStrong();
    break;
  default:
    abort();
  }
  return value;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
