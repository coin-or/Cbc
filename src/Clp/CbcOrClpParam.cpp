// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinPragma.hpp"
#include "CoinTime.hpp"
//#define COIN_HAS_CBC
#define COIN_HAS_CLP
#include "CbcOrClpParam.hpp"
#include "CoinHelperFunctions.hpp"

#include <string>
#include <iostream>
#include <cassert>
#include "CoinFinite.hpp"

#ifdef COIN_HAS_CBC
#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#include "ClpSimplex.hpp"
#endif
#include "CbcModel.hpp"
#endif

#ifdef COIN_HAS_CLP
#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"
#endif

#ifdef COIN_HAS_CBC
// from CoinSolve
static char coin_prompt[] = "Coin:";
#else
static char coin_prompt[] = "Clp:";
#endif

#ifdef CLP_CILK
#ifndef CBC_THREAD
#define CBC_THREAD
#endif
#endif



static bool doPrinting = true;
static std::string afterEquals = "";
static char printArray[250];
#if COIN_INT_MAX == 0
#undef COIN_INT_MAX
#define COIN_INT_MAX 2147483647
#endif
#if FLUSH_PRINT_BUFFER > 2
int coinFlushBufferFlag = 0;
#endif
void setCbcOrClpPrinting(bool yesNo)
{
  doPrinting = yesNo;
}
// Returns next valid field
static int CbcOrClpRead_mode = 1;
static FILE *CbcOrClpReadCommand = stdin;
int getCbcOrClpReadMode() { return CbcOrClpRead_mode; }
void setCbcOrClpReadMode(int mode) { CbcOrClpRead_mode = mode; }
void setCbcOrClpReadCommand(FILE* f) { CbcOrClpReadCommand = f; }
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CbcOrClpParam::CbcOrClpParam()
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
  , display_(0)
  , intValue_(-1)
  , doubleValue_(-1.0)
  , stringValue_("")
  , whereUsed_(7)
  , fakeKeyWord_(-1)
  , fakeValue_(0)
{
}
// Other constructors
CbcOrClpParam::CbcOrClpParam(std::string name, std::string help,
  double lower, double upper, CbcOrClpParameterType type,
  int display)
  : type_(type)
  , lowerIntValue_(0)
  , upperIntValue_(0)
  , definedKeyWords_()
  , name_(name)
  , shortHelp_(help)
  , longHelp_()
  , action_(type)
  , currentKeyWord_(-1)
  , display_(display)
  , intValue_(-1)
  , doubleValue_(-1.0)
  , stringValue_("")
  , whereUsed_(7)
  , fakeKeyWord_(-1)
  , fakeValue_(0)
{
  lowerDoubleValue_ = lower;
  upperDoubleValue_ = upper;
  gutsOfConstructor();
}
CbcOrClpParam::CbcOrClpParam(std::string name, std::string help,
  int lower, int upper, CbcOrClpParameterType type,
  int display)
  : type_(type)
  , lowerDoubleValue_(0.0)
  , upperDoubleValue_(0.0)
  , definedKeyWords_()
  , name_(name)
  , shortHelp_(help)
  , longHelp_()
  , action_(type)
  , currentKeyWord_(-1)
  , display_(display)
  , intValue_(-1)
  , doubleValue_(-1.0)
  , stringValue_("")
  , whereUsed_(7)
  , fakeKeyWord_(-1)
  , fakeValue_(0)
{
  gutsOfConstructor();
  lowerIntValue_ = lower;
  upperIntValue_ = upper;
}
// Other strings will be added by append
CbcOrClpParam::CbcOrClpParam(std::string name, std::string help,
  std::string firstValue,
  CbcOrClpParameterType type, int whereUsed,
  int display)
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
  , display_(display)
  , intValue_(-1)
  , doubleValue_(-1.0)
  , stringValue_("")
  , whereUsed_(whereUsed)
  , fakeKeyWord_(-1)
  , fakeValue_(0)
{
  gutsOfConstructor();
  definedKeyWords_.push_back(firstValue);
}
// Action
CbcOrClpParam::CbcOrClpParam(std::string name, std::string help,
  CbcOrClpParameterType type, int whereUsed,
  int display)
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
  , display_(display)
  , intValue_(-1)
  , doubleValue_(-1.0)
  , stringValue_("")
  , fakeKeyWord_(-1)
  , fakeValue_(0)
{
  whereUsed_ = whereUsed;
  gutsOfConstructor();
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CbcOrClpParam::CbcOrClpParam(const CbcOrClpParam &rhs)
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
  stringValue_ = rhs.stringValue_;
  whereUsed_ = rhs.whereUsed_;
  fakeKeyWord_ = rhs.fakeKeyWord_;
  fakeValue_ = rhs.fakeValue_;
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CbcOrClpParam::~CbcOrClpParam()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CbcOrClpParam &
CbcOrClpParam::operator=(const CbcOrClpParam &rhs)
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
    stringValue_ = rhs.stringValue_;
    whereUsed_ = rhs.whereUsed_;
    fakeKeyWord_ = rhs.fakeKeyWord_;
    fakeValue_ = rhs.fakeValue_;
  }
  return *this;
}
void CbcOrClpParam::gutsOfConstructor()
{
  std::string::size_type shriekPos = name_.find('!');
  lengthName_ = static_cast< unsigned int >(name_.length());
  if (shriekPos == std::string::npos) {
    //does not contain '!'
    lengthMatch_ = lengthName_;
  } else {
    lengthMatch_ = static_cast< unsigned int >(shriekPos);
    name_ = name_.substr(0, shriekPos) + name_.substr(shriekPos + 1);
    lengthName_--;
  }
}
// Sets value of fake keyword to current size of keywords
void CbcOrClpParam::setFakeKeyWord(int fakeValue)
{
  fakeKeyWord_ = static_cast< int >(definedKeyWords_.size());
  assert(fakeKeyWord_ > 0);
  fakeValue_ = fakeValue;
  assert(fakeValue_ >= 0);
}
/* Returns current parameter option position
   but if fake keyword returns fakeValue_
*/
int CbcOrClpParam::currentOptionAsInteger() const
{
  int fakeInteger;
  return currentOptionAsInteger(fakeInteger);
}
/* Returns current parameter option position
   but if fake keyword returns fakeValue_ and sets
   fakeInteger to value
*/
int CbcOrClpParam::currentOptionAsInteger(int &fakeInteger) const
{
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
// Returns length of name for printing
int CbcOrClpParam::lengthMatchName() const
{
  if (lengthName_ == lengthMatch_)
    return lengthName_;
  else
    return lengthName_ + 2;
}
// Insert string (only valid for keywords)
void CbcOrClpParam::append(std::string keyWord)
{
  definedKeyWords_.push_back(keyWord);
}

int CbcOrClpParam::matches(std::string input) const
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
CbcOrClpParam::matchName() const
{
  if (lengthMatch_ == lengthName_)
    return name_;
  else
    return name_.substr(0, lengthMatch_) + "(" + name_.substr(lengthMatch_) + ")";
}

// Returns parameter option which matches (-1 if none)
int CbcOrClpParam::parameterOption(std::string check) const
{
  int numberItems = static_cast< int >(definedKeyWords_.size());
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
        //contains '!'
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
      } else if (check.substr(0, 5) == "minus" || check.substr(0, 5) == "MINUS") {
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
        value = static_cast< int >(strtol(start, &endPointer, 10));
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
void CbcOrClpParam::printOptions() const
{
  std::cout << "<Possible options for " << name_ << " are:";
  unsigned int it;
  for (it = 0; it < definedKeyWords_.size(); it++) {
    std::string thisOne = definedKeyWords_[it];
    std::string::size_type shriekPos = thisOne.find('!');
    if (shriekPos != std::string::npos) {
      //contains '!'
      thisOne = thisOne.substr(0, shriekPos) + "(" + thisOne.substr(shriekPos + 1) + ")";
    }
    std::cout << " " << thisOne;
  }
  assert(currentKeyWord_ >= 0 && currentKeyWord_ < static_cast< int >(definedKeyWords_.size()));
  std::string current = definedKeyWords_[currentKeyWord_];
  std::string::size_type shriekPos = current.find('!');
  if (shriekPos != std::string::npos) {
    //contains '!'
    current = current.substr(0, shriekPos) + "(" + current.substr(shriekPos + 1) + ")";
  }
  std::cout << ";\n\tcurrent  " << current << ">" << std::endl;
}
// Print action and string
void CbcOrClpParam::printString() const
{
  if (name_ == "directory")
    std::cout << "Current working directory is " << stringValue_ << std::endl;
  else if (name_.substr(0, 6) == "printM")
    std::cout << "Current value of printMask is " << stringValue_ << std::endl;
  else
    std::cout << "Current default (if $ as parameter) for " << name_
              << " is " << stringValue_ << std::endl;
}
void CoinReadPrintit(const char *input)
{
  int length = static_cast< int >(strlen(input));
  assert(length <= 1000);
  char temp[1001];
  int i;
  int n = 0;
  for (i = 0; i < length; i++) {
    if (input[i] == '\n') {
      temp[n] = '\0';
      std::cout << temp << std::endl;
      n = 0;
    } else if (n >= 65 && input[i] == ' ') {
      temp[n] = '\0';
      std::cout << temp << std::endl;
      n = 0;
    } else if (n || input[i] != ' ') {
      temp[n++] = input[i];
    }
  }
  if (n) {
    temp[n] = '\0';
    std::cout << temp << std::endl;
  }
}
// Print Long help
void CbcOrClpParam::printLongHelp() const
{
  if (type_ >= 1 && type_ < 600) {
    CoinReadPrintit(longHelp_.c_str());
    if (type_ < CLP_PARAM_INT_SOLVERLOGLEVEL) {
      printf("<Range of values is %g to %g;\n\tcurrent %g>\n", lowerDoubleValue_, upperDoubleValue_, doubleValue_);
      assert(upperDoubleValue_ > lowerDoubleValue_);
    } else if (type_ < CLP_PARAM_STR_DIRECTION) {
      printf("<Range of values is %d to %d; current %d>\n", lowerIntValue_, upperIntValue_, intValue_);
      assert(upperIntValue_ > lowerIntValue_);
      if (stringValue_!="") {
	// print options
	size_t last = stringValue_.find_last_of('#');
	if (stringValue_[last+1]=='+')
	  std::cout << "Also keywords (optionally with + ) - " << std::endl;
	else
	  std::cout << "Also keywords - " << std::endl;
	size_t current = 0;
	while (current<last) {
	  size_t next = stringValue_.find_first_of('#',current);
	  std::string temp = stringValue_.substr(current,next-current);
	  current=next+1;
	  size_t end1 = temp.find_first_of('[');
	  std::string name = temp.substr(0,end1);
	  end1 = temp.find_first_of('[',1);
	  end1 = temp.find_first_of('[',end1+1);
	  std::string comment =
	    temp.substr(end1+1,temp.find_first_of(']',end1)-end1-1);
	  std::cout << name << " - " << comment << std::endl;
	}
      }
    } else if (type_ < CLP_PARAM_ACTION_DIRECTORY) {
      printOptions();
    }
  }
}
#ifdef COIN_HAS_CBC
int CbcOrClpParam::setDoubleParameter(OsiSolverInterface *model, double value)
{
  int returnCode;
  setDoubleParameterWithMessage(model, value, returnCode);
  if (doPrinting && strlen(printArray))
    std::cout << printArray << std::endl;
  return returnCode;
}
// Sets double parameter and returns printable string and error code
const char *
CbcOrClpParam::setDoubleParameterWithMessage(OsiSolverInterface *model, double value, int &returnCode)
{
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
    sprintf(printArray, "%g was provided for %s - valid range is %g to %g",
      value, name_.c_str(), lowerDoubleValue_, upperDoubleValue_);
    std::cout << value << " was provided for " << name_ << " - valid range is " << lowerDoubleValue_ << " to " << upperDoubleValue_ << std::endl;
    returnCode = 1;
  } else {
    double oldValue = doubleValue_;
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
    sprintf(printArray, "%s was changed from %g to %g",
      name_.c_str(), oldValue, value);
    returnCode = 0;
  }
  return printArray;
}
#endif
#ifdef COIN_HAS_CLP
int CbcOrClpParam::setDoubleParameter(ClpSimplex *model, double value)
{
  int returnCode;
  setDoubleParameterWithMessage(model, value, returnCode);
  if (doPrinting && strlen(printArray))
    std::cout << printArray << std::endl;
  return returnCode;
}
// Sets int parameter and returns printable string and error code
const char *
CbcOrClpParam::setDoubleParameterWithMessage(ClpSimplex *model, double value, int &returnCode)
{
  double oldValue = doubleValue_;
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
    sprintf(printArray, "%g was provided for %s - valid range is %g to %g",
      value, name_.c_str(), lowerDoubleValue_, upperDoubleValue_);
    returnCode = 1;
  } else {
    sprintf(printArray, "%s was changed from %g to %g",
      name_.c_str(), oldValue, value);
    returnCode = 0;
    doubleValue_ = value;
    switch (type_) {
    case CLP_PARAM_DBL_DUALTOLERANCE:
      model->setDualTolerance(value);
      break;
    case CLP_PARAM_DBL_PRIMALTOLERANCE:
      model->setPrimalTolerance(value);
      break;
    case CLP_PARAM_DBL_ZEROTOLERANCE:
      model->setSmallElementValue(value);
      break;
    case CLP_PARAM_DBL_DUALBOUND:
      model->setDualBound(value);
      break;
    case CLP_PARAM_DBL_PRIMALWEIGHT:
      model->setInfeasibilityCost(value);
      break;
#ifndef COIN_HAS_CBC
    case CLP_PARAM_DBL_TIMELIMIT:
      model->setMaximumSeconds(value);
      break;
#endif
    case CLP_PARAM_DBL_OBJSCALE:
      model->setObjectiveScale(value);
      break;
    case CLP_PARAM_DBL_RHSSCALE:
      model->setRhsScale(value);
      break;
    case CLP_PARAM_DBL_PRESOLVETOLERANCE:
      model->setDblParam(ClpPresolveTolerance, value);
      break;
    case CLP_PARAM_DBL_PROGRESS:
      model->setMinIntervalProgressUpdate(value);
      break;
    default:
      break;
    }
  }
  return printArray;
}
double
CbcOrClpParam::doubleParameter(ClpSimplex *model) const
{
  double value;
  switch (type_) {
#ifndef COIN_HAS_CBC
  case CLP_PARAM_DBL_DUALTOLERANCE:
    value = model->dualTolerance();
    break;
  case CLP_PARAM_DBL_PRIMALTOLERANCE:
    value = model->primalTolerance();
    break;
#endif
  case CLP_PARAM_DBL_ZEROTOLERANCE:
    value = model->getSmallElementValue();
    break;
  case CLP_PARAM_DBL_DUALBOUND:
    value = model->dualBound();
    break;
  case CLP_PARAM_DBL_PRIMALWEIGHT:
    value = model->infeasibilityCost();
    break;
#ifndef COIN_HAS_CBC
  case CLP_PARAM_DBL_TIMELIMIT:
    value = model->maximumSeconds();
    break;
#endif
  case CLP_PARAM_DBL_OBJSCALE:
    value = model->objectiveScale();
    break;
  case CLP_PARAM_DBL_RHSSCALE:
    value = model->rhsScale();
    break;
  case CLP_PARAM_DBL_PRESOLVETOLERANCE:
    value = model->presolveTolerance();
    break;
  default:
    value = doubleValue_;
    break;
  }
  return value;
}
int CbcOrClpParam::setIntParameter(ClpSimplex *model, int value)
{
  int returnCode;
  setIntParameterWithMessage(model, value, returnCode);
  if (doPrinting && strlen(printArray))
    std::cout << printArray << std::endl;
  return returnCode;
}
// Sets int parameter and returns printable string and error code
const char *
CbcOrClpParam::setIntParameterWithMessage(ClpSimplex *model, int value, int &returnCode)
{
  int oldValue = intValue_;
  if (value < lowerIntValue_ || value > upperIntValue_) {
    sprintf(printArray, "%d was provided for %s - valid range is %d to %d",
      value, name_.c_str(), lowerIntValue_, upperIntValue_);
    returnCode = 1;
  } else {
    intValue_ = value;
    sprintf(printArray, "%s was changed from %d to %d",
      name_.c_str(), oldValue, value);
    returnCode = 0;
    switch (type_) {
    case CLP_PARAM_INT_SOLVERLOGLEVEL:
      model->setLogLevel(value);
      if (value > 2)
        model->factorization()->messageLevel(8);
      else
        model->factorization()->messageLevel(0);
      break;
    case CLP_PARAM_INT_MAXFACTOR:
      model->factorization()->maximumPivots(value);
      break;
    case CLP_PARAM_INT_PERTVALUE:
      model->setPerturbation(value);
      break;
    case CLP_PARAM_INT_MAXITERATION:
      model->setMaximumIterations(value);
      break;
    case CLP_PARAM_INT_SPECIALOPTIONS:
      model->setSpecialOptions(value);
      break;
    case CLP_PARAM_INT_RANDOMSEED: {
      if (value == 0) {
        double time = fabs(CoinGetTimeOfDay());
        while (time >= COIN_INT_MAX)
          time *= 0.5;
        value = static_cast< int >(time);
        sprintf(printArray, "using time of day %s was changed from %d to %d",
          name_.c_str(), oldValue, value);
      }
      model->setRandomSeed(value);
    } break;
    case CLP_PARAM_INT_MORESPECIALOPTIONS:
      model->setMoreSpecialOptions(value);
      break;
    case CLP_PARAM_INT_VECTOR_MODE:
      model->setVectorMode(value);
      break;
#ifndef COIN_HAS_CBC
#ifdef CBC_THREAD
    case CBC_PARAM_INT_THREADS:
      model->setNumberThreads(value);
      break;
#endif
#endif
    default:
      break;
    }
  }
  return printArray;
}
int CbcOrClpParam::intParameter(ClpSimplex *model) const
{
  int value;
  switch (type_) {
#ifndef COIN_HAS_CBC
  case CLP_PARAM_INT_SOLVERLOGLEVEL:
    value = model->logLevel();
    break;
#endif
  case CLP_PARAM_INT_MAXFACTOR:
    value = model->factorization()->maximumPivots();
    break;
    break;
  case CLP_PARAM_INT_PERTVALUE:
    value = model->perturbation();
    break;
  case CLP_PARAM_INT_MAXITERATION:
    value = model->maximumIterations();
    break;
  case CLP_PARAM_INT_SPECIALOPTIONS:
    value = model->specialOptions();
    break;
  case CLP_PARAM_INT_RANDOMSEED:
    value = model->randomNumberGenerator()->getSeed();
    break;
  case CLP_PARAM_INT_MORESPECIALOPTIONS:
    value = model->moreSpecialOptions();
    break;
  case CLP_PARAM_INT_VECTOR_MODE:
    value = model->vectorMode();
    break;
#ifndef COIN_HAS_CBC
#ifdef CBC_THREAD
  case CBC_PARAM_INT_THREADS:
    value = model->numberThreads();
    break;
#endif
#endif
  default:
    value = intValue_;
    break;
  }
  return value;
}
#endif
int CbcOrClpParam::checkDoubleParameter(double value) const
{
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is " << lowerDoubleValue_ << " to " << upperDoubleValue_ << std::endl;
    return 1;
  } else {
    return 0;
  }
}
#ifdef COIN_HAS_CBC
double
CbcOrClpParam::doubleParameter(OsiSolverInterface *model) const
{
  double value = 0.0;
  switch (type_) {
  case CLP_PARAM_DBL_DUALTOLERANCE:
    model->getDblParam(OsiDualTolerance, value);
    break;
  case CLP_PARAM_DBL_PRIMALTOLERANCE:
    model->getDblParam(OsiPrimalTolerance, value);
    break;
  default:
    return doubleValue_;
    break;
  }
  return value;
}
int CbcOrClpParam::setIntParameter(OsiSolverInterface *model, int value)
{
  int returnCode;
  setIntParameterWithMessage(model, value, returnCode);
  if (doPrinting && strlen(printArray))
    std::cout << printArray << std::endl;
  return returnCode;
}
// Sets int parameter and returns printable string and error code
const char *
CbcOrClpParam::setIntParameterWithMessage(OsiSolverInterface *model, int value, int &returnCode)
{
  if (value < lowerIntValue_ || value > upperIntValue_) {
    sprintf(printArray, "%d was provided for %s - valid range is %d to %d",
      value, name_.c_str(), lowerIntValue_, upperIntValue_);
    returnCode = 1;
  } else {
    int oldValue = intValue_;
    intValue_ = oldValue;
    switch (type_) {
    case CLP_PARAM_INT_SOLVERLOGLEVEL:
      model->messageHandler()->setLogLevel(value);
      break;
    default:
      break;
    }
    sprintf(printArray, "%s was changed from %d to %d",
      name_.c_str(), oldValue, value);
    returnCode = 0;
  }
  return printArray;
}
int CbcOrClpParam::intParameter(OsiSolverInterface *model) const
{
  int value = 0;
  switch (type_) {
  case CLP_PARAM_INT_SOLVERLOGLEVEL:
    value = model->messageHandler()->logLevel();
    break;
  default:
    value = intValue_;
    break;
  }
  return value;
}
int CbcOrClpParam::setDoubleParameter(CbcModel &model, double value)
{
  int returnCode = 0;
  setDoubleParameterWithMessage(model, value, returnCode);
  if (doPrinting && strlen(printArray))
    std::cout << printArray << std::endl;
  return returnCode;
}
// Sets double parameter and returns printable string and error code
const char *
CbcOrClpParam::setDoubleParameterWithMessage(CbcModel &model, double value, int &returnCode)
{
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
    sprintf(printArray, "%g was provided for %s - valid range is %g to %g",
      value, name_.c_str(), lowerDoubleValue_, upperDoubleValue_);
    returnCode = 1;
  } else {
    double oldValue = doubleValue_;
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
    sprintf(printArray, "%s was changed from %g to %g",
      name_.c_str(), oldValue, value);
    returnCode = 0;
  }
  return printArray;
}
double
CbcOrClpParam::doubleParameter(CbcModel &model) const
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
int CbcOrClpParam::setIntParameter(CbcModel &model, int value)
{
  int returnCode;
  setIntParameterWithMessage(model, value, returnCode);
  if (doPrinting && strlen(printArray))
    std::cout << printArray << std::endl;
  return returnCode;
}
// Sets int parameter and returns printable string and error code
const char *
CbcOrClpParam::setIntParameterWithMessage(CbcModel &model, int value, int &returnCode)
{
  if (value < lowerIntValue_ || value > upperIntValue_) {
    sprintf(printArray, "%d was provided for %s - valid range is %d to %d",
      value, name_.c_str(), lowerIntValue_, upperIntValue_);
    returnCode = 1;
  } else {
    printArray[0] = '\0';
    if (value == intValue_)
      return printArray;
    int oldValue = intValue_;
    intValue_ = value;
    switch (type_) {
    case CLP_PARAM_INT_LOGLEVEL:
      oldValue = model.messageHandler()->logLevel();
      model.messageHandler()->setLogLevel(std::abs(value));
      break;
    case CLP_PARAM_INT_SOLVERLOGLEVEL:
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
    sprintf(printArray, "%s was changed from %d to %d",
      name_.c_str(), oldValue, value);
    returnCode = 0;
  }
  return printArray;
}
int CbcOrClpParam::intParameter(CbcModel &model) const
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
#endif
// Sets current parameter option using string
void CbcOrClpParam::setCurrentOption(const std::string value)
{
  int action = parameterOption(value);
  if (action >= 0)
    currentKeyWord_ = action;
#if FLUSH_PRINT_BUFFER > 2
  if (name_ == "bufferedMode")
    coinFlushBufferFlag = action;
#endif
}
// Sets current parameter option
void CbcOrClpParam::setCurrentOption(int value, bool printIt)
{
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
const char *
CbcOrClpParam::setCurrentOptionWithMessage(int value)
{
  if (value != currentKeyWord_) {
    char current[100];
    char newString[100];
    if (currentKeyWord_ >= 0 && (fakeKeyWord_ <= 0 || currentKeyWord_ < fakeKeyWord_))
      strcpy(current, definedKeyWords_[currentKeyWord_].c_str());
    else if (currentKeyWord_ < 0)
      sprintf(current, "minus%d", -currentKeyWord_ - 1000);
    else
      sprintf(current, "plus%d", currentKeyWord_ - 1000);
    if (value >= 0 && (fakeKeyWord_ <= 0 || value < fakeKeyWord_))
      strcpy(newString, definedKeyWords_[value].c_str());
    else if (value < 0)
      sprintf(newString, "minus%d", -value - 1000);
    else
      sprintf(newString, "plus%d", value - 1000);
    sprintf(printArray, "Option for %s changed from %s to %s",
      name_.c_str(), current, newString);
#if FLUSH_PRINT_BUFFER > 2
    if (name_ == "bufferedMode")
      coinFlushBufferFlag = value;
#endif
    currentKeyWord_ = value;
  } else {
    printArray[0] = '\0';
  }
  return printArray;
}
// Sets current parameter option using string with message
const char *
CbcOrClpParam::setCurrentOptionWithMessage(const std::string value)
{
  int action = parameterOption(value);
  char current[100];
  printArray[0] = '\0';
  if (action >= 0) {
#if FLUSH_PRINT_BUFFER > 2
    if (name_ == "bufferedMode")
      coinFlushBufferFlag = action;
#endif
    if (action == currentKeyWord_)
      return NULL;
    if (currentKeyWord_ >= 0 && (fakeKeyWord_ <= 0 || currentKeyWord_ < fakeKeyWord_))
      strcpy(current, definedKeyWords_[currentKeyWord_].c_str());
    else if (currentKeyWord_ < 0)
      sprintf(current, "minus%d", -currentKeyWord_ - 1000);
    else
      sprintf(current, "plus%d", currentKeyWord_ - 1000);
    sprintf(printArray, "Option for %s changed from %s to %s",
      name_.c_str(), current, value.c_str());
    currentKeyWord_ = action;
  } else {
    sprintf(printArray, "Option for %s given illegal value %s",
      name_.c_str(), value.c_str());
  }
  return printArray;
}
void CbcOrClpParam::setIntValue(int value)
{
  if (value < lowerIntValue_ || value > upperIntValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is " << lowerIntValue_ << " to " << upperIntValue_ << std::endl;
  } else {
    intValue_ = value;
  }
}
const char *
CbcOrClpParam::setIntValueWithMessage(int value)
{
  printArray[0] = '\0';
  if (value < lowerIntValue_ || value > upperIntValue_) {
    sprintf(printArray, "%d was provided for %s - valid range is %d to %d",
      value, name_.c_str(), lowerIntValue_, upperIntValue_);
  } else {
    if (value == intValue_)
      return NULL;
    sprintf(printArray, "%s was changed from %d to %d",
      name_.c_str(), intValue_, value);
    intValue_ = value;
  }
  return printArray;
}
void CbcOrClpParam::setDoubleValue(double value)
{
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is " << lowerDoubleValue_ << " to " << upperDoubleValue_ << std::endl;
  } else {
    doubleValue_ = value;
  }
}
const char *
CbcOrClpParam::setDoubleValueWithMessage(double value)
{
  printArray[0] = '\0';
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
    sprintf(printArray, "%g was provided for %s - valid range is %g to %g",
      value, name_.c_str(), lowerDoubleValue_, upperDoubleValue_);
  } else {
    if (value == doubleValue_)
      return NULL;
    sprintf(printArray, "%s was changed from %g to %g",
      name_.c_str(), doubleValue_, value);
    doubleValue_ = value;
  }
  return printArray;
}
void CbcOrClpParam::setStringValue(std::string value)
{
  stringValue_ = value;
}
void CbcOrClpParam::appendStringValue(std::string value)
{
  stringValue_ = stringValue_ + value;
}
static char line[1000];
static char *where = NULL;
extern int CbcOrClpRead_mode;
int CbcOrClpEnvironmentIndex = -1;
// Alternative to environment
char *alternativeEnvironment = NULL;
static size_t fillEnv()
{
#if defined(_MSC_VER) || defined(__MSVCRT__)
  return 0;
#else
  // Don't think it will work on Windows
  char *environ;
  if (!alternativeEnvironment)
    environ = getenv("CBC_CLP_ENVIRONMENT");
  else
    environ = alternativeEnvironment;
  size_t length = 0;
  if (environ) {
    length = strlen(environ);
    if (CbcOrClpEnvironmentIndex < static_cast< int >(length)) {
      // find next non blank
      char *whereEnv = environ + CbcOrClpEnvironmentIndex;
      // munch white space
      while (*whereEnv == ' ' || *whereEnv == '\t' || *whereEnv < ' ')
        whereEnv++;
      // copy
      char *put = line;
      while (*whereEnv != '\0') {
        if (*whereEnv == ' ' || *whereEnv == '\t' || *whereEnv < ' ') {
          break;
        }
        *put = *whereEnv;
        put++;
        assert(put - line < 1000);
        whereEnv++;
      }
      CbcOrClpEnvironmentIndex = static_cast< int >(whereEnv - environ);
      *put = '\0';
      length = strlen(line);
    } else {
      length = 0;
    }
  }
  if (!length) {
    CbcOrClpEnvironmentIndex = -1;
    if (alternativeEnvironment) {
      delete[] alternativeEnvironment;
      alternativeEnvironment = NULL;
    }
  }
  return length;
#endif
}
extern FILE *CbcOrClpReadCommand;
// Simple read stuff
std::string
CoinReadNextField()
{
  std::string field;
  if (!where) {
    // need new line
    if (CbcOrClpReadCommand == stdin) {
      fputs(coin_prompt, stdout);
      fflush(stdout);
    }
    where = fgets(line, 1000, CbcOrClpReadCommand);
    if (!where)
      return field; // EOF
    where = line;
    // clean image
    char *lastNonBlank = line - 1;
    while (*where != '\0') {
      if (*where != '\t' && *where < ' ') {
        break;
      } else if (*where != '\t' && *where != ' ') {
        lastNonBlank = where;
      }
      where++;
    }
    where = line;
    *(lastNonBlank + 1) = '\0';
  }
  // munch white space
  while (*where == ' ' || *where == '\t')
    where++;
  char *saveWhere = where;
  while (*where != ' ' && *where != '\t' && *where != '\0')
    where++;
  if (where != saveWhere) {
    char save = *where;
    *where = '\0';
    //convert to string
    field = saveWhere;
    *where = save;
  } else {
    where = NULL;
    field = "EOL";
  }
  return field;
}

std::string
CoinReadGetCommand(int argc, const char *argv[])
{
  std::string field = "EOL";
  // say no =
  afterEquals = "";
  while (field == "EOL") {
    if (CbcOrClpRead_mode > 0) {
      if ((CbcOrClpRead_mode < argc && argv[CbcOrClpRead_mode]) || CbcOrClpEnvironmentIndex >= 0) {
        if (CbcOrClpEnvironmentIndex < 0) {
          field = argv[CbcOrClpRead_mode++];
        } else {
          if (fillEnv()) {
            field = line;
          } else {
            // not there
            continue;
          }
        }
        if (field == "-") {
          std::cout << "Switching to line mode" << std::endl;
          CbcOrClpRead_mode = -1;
          field = CoinReadNextField();
        } else if (field[0] != '-') {
          if (CbcOrClpRead_mode != 2) {
            // now allow std::cout<<"skipping non-command "<<field<<std::endl;
            // field="EOL"; // skip
          } else if (CbcOrClpEnvironmentIndex < 0) {
            // special dispensation - taken as -import name
            CbcOrClpRead_mode--;
            field = "import";
          }
        } else {
          if (field != "--") {
            // take off -
            field = field.substr(1);
          } else {
            // special dispensation - taken as -import --
            CbcOrClpRead_mode--;
            field = "import";
          }
        }
      } else {
        field = "";
      }
    } else {
      field = CoinReadNextField();
    }
  }
  // if = then modify and save
  std::string::size_type found = field.find('=');
  if (found != std::string::npos) {
    afterEquals = field.substr(found + 1);
    field = field.substr(0, found);
  }
  //std::cout<<field<<std::endl;
  return field;
}
std::string
CoinReadGetString(int argc, const char *argv[])
{
  std::string field = "EOL";
  if (afterEquals == "") {
    if (CbcOrClpRead_mode > 0) {
      if (CbcOrClpRead_mode < argc || CbcOrClpEnvironmentIndex >= 0) {
        if (CbcOrClpEnvironmentIndex < 0) {
	  const char * input = argv[CbcOrClpRead_mode];
	  if (strcmp(input,"--")&&strcmp(input,"stdin")&&
	      strcmp(input,"stdin_lp")) {
            field = argv[CbcOrClpRead_mode++];
          } else {
            CbcOrClpRead_mode++;
            // -- means import from stdin
	    // but allow for other than mps files
	    // Clp does things in different way !!
	    if (!strcmp(input,"--"))
	      field = "-";
	    else if (!strcmp(input,"stdin"))
	      field = "-";
	    else if (!strcmp(input,"stdin_lp"))
	      field = "-lp";
	  }
        } else {
          fillEnv();
          field = line;
        }
      }
    } else {
      field = CoinReadNextField();
    }
  } else {
    field = afterEquals;
    afterEquals = "";
  }
  //std::cout<<field<<std::endl;
  return field;
}
static std::string errorField="";
// valid 0 - okay, 1 bad, 2 not there
int CoinReadGetIntField(int argc, const char *argv[], int *valid)
{
  std::string field = "EOL";
  if (afterEquals == "") {
    if (CbcOrClpRead_mode > 0) {
      if (CbcOrClpRead_mode < argc || CbcOrClpEnvironmentIndex >= 0) {
        if (CbcOrClpEnvironmentIndex < 0) {
          // may be negative value so do not check for -
          field = argv[CbcOrClpRead_mode++];
        } else {
          fillEnv();
          field = line;
        }
      }
    } else {
      field = CoinReadNextField();
    }
  } else {
    field = afterEquals;
    afterEquals = "";
  }
  long int value = 0;
  //std::cout<<field<<std::endl;
  if (field != "EOL") {
    const char *start = field.c_str();
    char *endPointer = NULL;
    // check valid
    value = strtol(start, &endPointer, 10);
    if (*endPointer == '\0') {
      *valid = 0;
    } else {
      *valid = 1;
      errorField = field;
      //std::cout << "String of " << field;
    }
  } else {
    *valid = 2;
  }
  return static_cast< int >(value);
}
std::string getCoinErrorField()
{
  return errorField;
}
// Decodes options
int
CbcOrClpParam::optionIntField(std::string field, int *valid)
{
  size_t last = stringValue_.find_last_of('#');
  char allowed = stringValue_[last+1];
  *valid = 0;
  int value = 0;
  while (field.length()) {
    std::string thisPart;
    size_t findSep = field.find_first_of(allowed);
    if (findSep != std::string::npos) {
      thisPart = field.substr(0,findSep);
      field = field.substr(findSep+1);
    } else {
      thisPart = field;
      field = "";
    }
    for (int i=0;i<field.length();i++) {
      char x = thisPart[i];
      x = tolower(static_cast<unsigned char>(x));
      thisPart[i] = x;
    }
    size_t current = 0;
    bool found = false;
    while (current<last) {
      size_t next = stringValue_.find_first_of('#',current);
      std::string temp = stringValue_.substr(current,next-current);
      int n = 0;
      size_t shriek = temp.find_first_of('[');
      for (int i=0;i<temp.length();i++) {
	if (temp[i]!='!') {
	  char x = temp[i];
	  x = tolower(static_cast<unsigned char>(x));
	  temp[n++] = x;
	} else {
	  shriek = n;
	}
      }
      current=next+1;
      bool foundThis = true;
      for (int i=0;i<thisPart.length();i++) {
	if (thisPart[i] != temp[i]) {
	  foundThis = false;
	  break;
	}
      }
      if (thisPart.length()<shriek)
	foundThis = false;
      if (foundThis) {
	found = true;
	size_t end1 = temp.find_first_of('[');
	size_t end2 = temp.find_first_of(']');
	// must be better way
	int thisValue = atoi(temp.substr(end1+1,end2-end1-1).c_str());
	assert (allowed=='+'||allowed=='=');
	if (allowed=='+') {
	  value += thisValue;
	} else {
	  if (value) {
	    std::cout << "Only one = item allowed" << std::endl;
	    found = false;
	  } else {
	    value = thisValue;
	  }
	}
	break;
      }
    }
    if (!found) {
      *valid = 1;
      return -1;
    }
  }
  return value;
}
double
CoinReadGetDoubleField(int argc, const char *argv[], int *valid)
{
  std::string field = "EOL";
  if (afterEquals == "") {
    if (CbcOrClpRead_mode > 0) {
      if (CbcOrClpRead_mode < argc || CbcOrClpEnvironmentIndex >= 0) {
        if (CbcOrClpEnvironmentIndex < 0) {
          // may be negative value so do not check for -
          field = argv[CbcOrClpRead_mode++];
        } else {
          fillEnv();
          field = line;
        }
      }
    } else {
      field = CoinReadNextField();
    }
  } else {
    field = afterEquals;
    afterEquals = "";
  }
  double value = 0.0;
  //std::cout<<field<<std::endl;
  if (field != "EOL") {
    const char *start = field.c_str();
    char *endPointer = NULL;
    // check valid
    value = strtod(start, &endPointer);
    if (*endPointer == '\0') {
      *valid = 0;
    } else {
      *valid = 1;
      std::cout << "String of " << field;
    }
  } else {
    *valid = 2;
  }
  return value;
}
/*
  Subroutine to establish the cbc parameter array. See the description of
  class CbcOrClpParam for details. Pulled from C..Main() for clarity.
*/
void establishParams(std::vector< CbcOrClpParam > &parameters)
{
  parameters.clear();
  parameters.push_back(CbcOrClpParam("?", "For help", CBC_PARAM_GENERALQUERY, 7, 0));
  parameters.push_back(CbcOrClpParam("???", "For help", CBC_PARAM_FULLGENERALQUERY, 7, 0));
  parameters.push_back(CbcOrClpParam("-", "From stdin", CLP_PARAM_ACTION_STDIN, 3, 0));

// some help strings that repeat for many options
#define CUTS_LONGHELP \
  "Value 'on' enables the cut generator and CBC will try it in the branch and cut tree (see cutDepth on how to fine tune the behavior). " \
  "Value 'root' lets CBC run the cut generator generate only at the root node. " \
  "Value 'ifmove' lets CBC use the cut generator in the tree if it looks as if it is doing some good and moves the objective value. " \
  "Value 'forceon' turns on the cut generator and forces CBC to use it at every node."
#define HEURISTICS_LONGHELP \
  "Value 'on' means to use the heuristic in each node of the tree, i.e. after preprocessing. " \
  "Value 'before' means use the heuristic only if option doHeuristics is used. " \
  "Value 'both' means to use the heuristic if option doHeuristics is used and during solve."

  {
    CbcOrClpParam p("allC!ommands", "Whether to print less used commands",
      "no", CLP_PARAM_STR_ALLCOMMANDS);

    p.append("more");
    p.append("all");
    p.setLonghelp(
      "For the sake of your sanity, only the more useful and simple commands \
are printed out on ?.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("allow!ableGap", "Stop when gap between best possible and \
best less than this",
      0.0, COIN_DBL_MAX, CBC_PARAM_DBL_ALLOWABLEGAP);
    p.setDoubleValue(1e-10);
    p.setLonghelp(
      "If the gap between best known solution and the best possible solution is less than this \
value, then the search will be terminated.  Also see ratioGap.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("allS!lack", "Set basis back to all slack and reset solution",
      CLP_PARAM_ACTION_ALLSLACK, 3);
    p.setLonghelp(
      "Mainly useful for tuning purposes.  Normally the first dual or primal will be using an all slack \
basis anyway.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("artif!icialCost", "Costs >= this treated as artificials in feasibility pump",
      0.0, COIN_DBL_MAX, CBC_PARAM_DBL_ARTIFICIALCOST, 1);
    p.setDoubleValue(0.0);
    p.setLonghelp(
      "A value of 0.0 means off. Otherwise, variables with costs >= this are treated as artificial variables and fixed to lower bound in feasibility pump.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("auto!Scale", "Whether to scale objective, rhs and bounds of problem if they look odd",
      "off", CLP_PARAM_STR_AUTOSCALE, 7, 0);
    p.append("on");
    p.setLonghelp(
      "If you think you may get odd objective values or large equality rows etc then\
 it may be worth setting this true.  It is still experimental and you may prefer\
 to use objective!Scale and rhs!Scale.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("barr!ier", "Solve using primal dual predictor corrector algorithm",
      CLP_PARAM_ACTION_BARRIER);
    p.setLonghelp(
      "This command solves the current model using the  primal dual predictor \
corrector algorithm.  You may want to link in an alternative \
ordering and factorization.  It will also solve models \
with quadratic objectives.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("basisI!n", "Import basis from bas file",
      CLP_PARAM_ACTION_BASISIN, 3);
    p.setLonghelp(
      "This will read an MPS format basis file from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  If you have libz then it can read compressed\
 files 'xxxxxxxx.gz' or xxxxxxxx.bz2.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("basisO!ut", "Export basis as bas file",
      CLP_PARAM_ACTION_BASISOUT);

    p.setLonghelp(
      "This will write an MPS format basis file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.bas'.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("biasLU", "Whether factorization biased towards U",
      "UU", CLP_PARAM_STR_BIASLU, 2, 0);

    p.append("UX");
    p.append("LX");
    p.append("LL");
    p.setCurrentOption("LX");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("branch!AndCut", "Do Branch and Cut",
      CBC_PARAM_ACTION_BAB);
    p.setLonghelp(
      "This does branch and cut.  There are many parameters which can affect the performance.  \
First just try with default settings and look carefully at the log file.  Did cuts help?  Did they take too long?  \
Look at output to see which cuts were effective and then do some tuning.  You will see that the \
options for cuts are off, on, root and ifmove, forceon.  Off is \
obvious. " CUTS_LONGHELP " For probing, forceonbut just does fixing probing in tree - not strengthening etc.  \
If pre-processing reduced the size of the \
problem or strengthened many coefficients then it is probably wise to leave it on.  Switch off heuristics \
which did not provide solutions.  The other major area to look at is the search.  Hopefully good solutions \
were obtained fairly early in the search so the important point is to select the best variable to branch on.  \
See whether strong branching did a good job - or did it just take a lot of iterations.  Adjust the strongBranching \
and trustPseudoCosts parameters.  If cuts did a good job, then you may wish to \
have more rounds of cuts - see passC!uts and passT!ree.");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("bscale", "Whether to scale in barrier (and ordering speed)",
      "off", CLP_PARAM_STR_BARRIERSCALE, 7, 0);
    p.append("on");
    p.append("off1");
    p.append("on1");
    p.append("off2");
    p.append("on2");
    parameters.push_back(p);
  }
#if FLUSH_PRINT_BUFFER > 2
  {
    CbcOrClpParam p("buff!eredMode", "Whether to flush print buffer",
      "on", CLP_PARAM_STR_BUFFER_MODE);

    p.append("off");
    p.setLonghelp(
      "Default is on, off switches on unbuffered output");
    p.setIntValue(0);
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("chol!esky", "Which cholesky algorithm",
      "native", CLP_PARAM_STR_CHOLESKY, 7);

    p.append("dense");
    //#ifdef FOREIGN_BARRIER
    p.append("fudge!Long_dummy");
    p.append("wssmp_dummy");
#if defined(CLP_HAS_AMD) || defined(CLP_HAS_CHOLMOD)
    p.append("Uni!versityOfFlorida");
#else
    p.append("Uni!versityOfFlorida_dummy");
#endif
    p.append("Taucs_dummy");
    p.append("Mumps_dummy");
#ifdef PARDISO_BARRIER
    p.append("Pardiso");
#else
    p.append("Pardiso_dummy");
#endif
    p.setLonghelp(
      "For a barrier code to be effective it needs a good Cholesky ordering and factorization.  \
The native ordering and factorization is not state of the art, although acceptable.  \
You may want to link in one from another source.  See Makefile.locations for some \
possibilities.");

    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("combine!Solutions", "Whether to use combine solution heuristic",
      "off", CBC_PARAM_STR_COMBINE);

    p.append("on");
    p.append("both");
    p.append("before");
    p.append("onquick");
    p.append("bothquick");
    p.append("beforequick");
    p.setLonghelp(
      "This heuristic does branch and cut on given problem by just \
using variables which have appeared in one or more solutions. \
It is obviously only tried after two or more solutions have been found. "
      HEURISTICS_LONGHELP);

    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("combine2!Solutions", "Whether to use crossover solution heuristic",
      "off", CBC_PARAM_STR_CROSSOVER2);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(
      "This heuristic does branch and cut on the problem given by \
fixing variables which have the same value in two or more solutions. \
It obviously only tries after two or more solutions. "
      HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("constraint!fromCutoff", "Whether to use cutoff as constraint",
      "off", CBC_PARAM_STR_CUTOFF_CONSTRAINT);

    p.append("on");
    p.append("variable");
    p.append("forcevariable");
    p.append("conflict");
    p.setLonghelp(
      "For some problems, cut generators and general branching work better if the problem would be infeasible if the cost is too high. "
      "If this option is enabled, the objective function is added as a constraint which right hand side is set to the current cutoff value (objective value of best known solution)");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("cost!Strategy", "How to use costs for branching priorities",
      "off", CBC_PARAM_STR_COSTSTRATEGY);

    p.append("pri!orities");
    p.append("column!Order?");
    p.append("01f!irst?");
    p.append("01l!ast?");
    p.append("length!?");
    p.append("singletons");
    p.append("nonzero");
    p.append("general!Force?");
    p.setLonghelp(
      "Value 'priorities' assigns highest priority to variables with largest absolute cost. This primitive strategy can be surprisingly effective. "
      "Value 'columnorder' assigns the priorities 1, 2, 3, ... with respect to the column ordering. "
      "Value '01first' ('01last') assignes two sets of priorities such that binary variables get high (low) priority. "
      "Value 'length' assigns high priority to variables that occur in many equations. ");

    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("cplex!Use", "Whether to use Cplex!",
      "off", CBC_PARAM_STR_CPX);
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
#endif
  {
    CbcOrClpParam p("cpp!Generate", "Generates C++ code",
      -1, 50000, CLP_PARAM_INT_CPP, 1);
    p.setLonghelp(
      "Once you like what the stand-alone solver does then this allows \
you to generate user_driver.cpp which approximates the code.  \
0 gives simplest driver, 1 generates saves and restores, 2 \
generates saves and restores even for variables at default value. \
4 bit in cbc generates size dependent code rather than computed values.  \
This is now deprecated as you can call stand-alone solver - see \
Cbc/examples/driver4.cpp.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("crash", "Whether to create basis for problem",
      "off", CLP_PARAM_STR_CRASH);

    p.append("on");
    p.append("so!low_halim");
    p.append("lots");
    p.append("free");
    p.append("zero");
    p.append("single!ton");
#ifdef CLP_INHERIT_MODE
    p.append("dual");
    p.append("dw");
    p.append("idiot");
#else
    p.append("idiot1");
    p.append("idiot2");
    p.append("idiot3");
    p.append("idiot4");
    p.append("idiot5");
    p.append("idiot6");
    p.append("idiot7");
#endif
    p.setLonghelp(
      "If crash is set to 'on' and there is an all slack basis then Clp will flip or put structural\
     variables into the basis with the aim of getting dual feasible.  On average, dual simplex seems to perform\
     better without it and there are alternative types of 'crash' for primal simplex, e.g. 'idiot' or 'sprint'. \
    A variant due to Solow and Halim which is as 'on' but just flips is also available.");

    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("cross!over", "Whether to get a basic solution with the simplex algorithm after the barrier algorithm finished",
      "on", CLP_PARAM_STR_CROSSOVER);
    p.append("off");
    p.append("maybe");
    p.append("presolve");
    p.setLonghelp(
      "Interior point algorithms do not obtain a basic solution.\
 This option will crossover to a basic solution suitable for ranging or branch and cut.  With the current state \
of the solver for quadratic programs it may be a good idea to switch off crossover for this case (and maybe \
presolve as well) - the option 'maybe' does this.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("csv!Statistics", "Create one line of statistics",
      CLP_PARAM_ACTION_CSVSTATISTICS, 2, 1);
    p.setLonghelp(
      "This appends statistics to given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  Adds header if file empty or does not exist.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("cutD!epth", "Depth in tree at which to do cuts",
      -1, COIN_INT_MAX, CBC_PARAM_INT_CUTDEPTH);
    p.setLonghelp(
      "Cut generators may be off, on, on only at the root node, or on if they look useful. \
      Setting this option to a positive value K let CBC call a cutgenerator on a node whenever the depth in the tree is a multiple of K. \
      The default of -1 lets CBC decide.");
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("cutL!ength", "Length of a cut",
      -1, COIN_INT_MAX, CBC_PARAM_INT_CUTLENGTH);
    p.setLonghelp("At present this only applies to Gomory cuts. -1 (default) leaves as is. \
Any value >0 says that all cuts <= this length can be generated both at \
root node and in tree. 0 says to use some dynamic lengths.  If value >=10,000,000 \
then the length in tree is value%10000000 - so 10000100 means unlimited length \
at root and 100 in tree.");
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("cuto!ff", "Bound on the objective value for all solutions",
      -COIN_DBL_MAX, COIN_DBL_MAX, CBC_PARAM_DBL_CUTOFF);
    p.setDoubleValue(1.0e50);
    p.setLonghelp(
      "All solutions must have a better objective value (in a minimization sense) than the value of this option.  \
CBC also updates this value whenever it obtains a solution to the value of \
the objective function of the solution minus the cutoff increment.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("cuts!OnOff", "Switches all cut generators on or off",
      "off", CBC_PARAM_STR_CUTSSTRATEGY);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.setLonghelp(
      "This can be used to switch on or off all cut generators (apart from Reduce and Split). "
      "Then one can turn individual ones off or on. "
      CUTS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("debug!In", "read valid solution from file",
      CLP_PARAM_ACTION_DEBUG, 7, 1);

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
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("decomp!ose", "Whether to try decomposition",
      -COIN_INT_MAX, COIN_INT_MAX, CLP_PARAM_INT_DECOMPOSE_BLOCKS, 1);
    p.setLonghelp(
      "0 - off, 1 choose blocks >1 use as blocks \
Dantzig Wolfe if primal, Benders if dual \
- uses sprint pass for number of passes");
    p.setIntValue(0);
    parameters.push_back(p);
  }
#if CLP_MULTIPLE_FACTORIZATIONS > 0
  {
    CbcOrClpParam p("dense!Threshold", "Threshold for using dense factorization",
      -1, 10000, CBC_PARAM_INT_DENSE, 1);
    p.setLonghelp(
      "If processed problem <= this use dense factorization");
    p.setIntValue(-1);
    parameters.push_back(p);
  }
#endif
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("depth!MiniBab", "Depth at which to try mini branch-and-bound",
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
    CbcOrClpParam p("dextra3", "Extra double parameter 3",
      -COIN_DBL_MAX, COIN_DBL_MAX, CBC_PARAM_DBL_DEXTRA3, 0);
    p.setDoubleValue(0.0);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("dextra4", "Extra double parameter 4",
      -COIN_DBL_MAX, COIN_DBL_MAX, CBC_PARAM_DBL_DEXTRA4, 0);
    p.setDoubleValue(0.0);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("dextra4", "Extra double parameter 5",
      -COIN_DBL_MAX, COIN_DBL_MAX, CBC_PARAM_DBL_DEXTRA5, 0);
    p.setDoubleValue(0.0);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("Dins", "Whether to try Distance Induced Neighborhood Search",
      "off", CBC_PARAM_STR_DINS);

    p.append("on");
    p.append("both");
    p.append("before");
    p.append("often");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("direction", "Minimize or Maximize",
      "min!imize", CLP_PARAM_STR_DIRECTION);
    p.append("max!imize");
    p.append("zero");
    p.setLonghelp(
      "The default is minimize - use 'direction maximize' for maximization.\n\
You can also use the parameters 'maximize' or 'minimize'.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("directory", "Set Default directory for import etc.",
      CLP_PARAM_ACTION_DIRECTORY);
    p.setLonghelp(
      "This sets the directory which import, export, saveModel, restoreModel etc will use.\
  It is initialized to './'");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("dirSample", "Set directory where the COIN-OR sample problems are.",
      CLP_PARAM_ACTION_DIRSAMPLE, 7, 1);

    p.setLonghelp(
      "This sets the directory where the COIN-OR sample problems reside. It is\
 used only when -unitTest is passed to clp. clp will pick up the test problems\
 from this directory.\
 It is initialized to '../../Data/Sample'");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("dirNetlib", "Set directory where the netlib problems are.",
      CLP_PARAM_ACTION_DIRNETLIB, 7, 1);

    p.setLonghelp(
      "This sets the directory where the netlib problems reside. One can get\
 the netlib problems from COIN-OR or from the main netlib site. This\
 parameter is used only when -netlib is passed to clp. clp will pick up the\
 netlib problems from this directory. If clp is built without zlib support\
 then the problems must be uncompressed.\
 It is initialized to '../../Data/Netlib'");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("dirMiplib", "Set directory where the miplib 2003 problems are.",
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
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("diveO!pt", "Diving options",
      -1, 200000, CBC_PARAM_INT_DIVEOPT, 1);
    p.setLonghelp(
      "If >2 && <20 then modify diving options - \
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
    CbcOrClpParam p("diveS!olves", "Diving solve option",
      -1, 200000, CBC_PARAM_INT_DIVEOPTSOLVES, 1);

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
    CbcOrClpParam p("DivingS!ome", "Whether to try Diving heuristics",
      "off", CBC_PARAM_STR_DIVINGS);

    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(
      "This switches on a random diving heuristic at various times. \
One may prefer to individually turn diving heuristics on or off. "
      HEURISTICS_LONGHELP);
// C - Coefficient, F - Fractional, G - Guided, L - LineSearch, P - PseudoCost, V - VectorLength.
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("DivingC!oefficient", "Whether to try Coefficient diving heuristic",
      "off", CBC_PARAM_STR_DIVINGC);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("DivingF!ractional", "Whether to try Fractional diving heuristic",
      "off", CBC_PARAM_STR_DIVINGF);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("DivingG!uided", "Whether to try Guided diving heuristic",
      "off", CBC_PARAM_STR_DIVINGG);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("DivingL!ineSearch", "Whether to try Linesearch diving heuristic",
      "off", CBC_PARAM_STR_DIVINGL);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("DivingP!seudoCost", "Whether to try Pseudocost diving heuristic",
      "off", CBC_PARAM_STR_DIVINGP);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("DivingV!ectorLength", "Whether to try Vectorlength diving heuristic",
      "off", CBC_PARAM_STR_DIVINGV);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("doH!euristic", "Do heuristics before any preprocessing",
      CBC_PARAM_ACTION_DOHEURISTIC, 3);
    p.setLonghelp(
      "Normally heuristics are done in branch and bound.  It may be useful to do them outside. \
Only those heuristics with 'both' or 'before' set will run.  \
Doing this may also set cutoff, which can help with preprocessing.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("dualB!ound", "Initially algorithm acts as if no \
gap between bounds exceeds this value",
      1.0e-20, 1.0e20, CLP_PARAM_DBL_DUALBOUND);
    p.setLonghelp(
      "The dual algorithm in Clp is a single phase algorithm as opposed to a two phase\
 algorithm where you first get feasible then optimal.  If a problem has both upper and\
 lower bounds then it is trivial to get dual feasible by setting non basic variables\
 to correct bound.  If the gap between the upper and lower bounds of a variable is more\
 than the value of dualBound Clp introduces fake bounds so that it can make the problem\
 dual feasible.  This has the same effect as a composite objective function in the\
 primal algorithm.  Too high a value may mean more iterations, while too low a bound means\
 the code may go all the way and then have to increase the bounds.  OSL had a heuristic to\
 adjust bounds, maybe we need that here.");

    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("dualize", "Solves dual reformulation",
      0, 4, CLP_PARAM_INT_DUALIZE, 1);
    p.setLonghelp(
      "Don't even think about it.");

    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("dualP!ivot", "Dual pivot choice algorithm",
      "auto!matic", CLP_PARAM_STR_DUALPIVOT, 7, 1);
    p.append("dant!zig");
    p.append("partial");
    p.append("steep!est");
    p.append("PEsteep!est");
    p.append("PEdantzig");
    p.setLonghelp(
      "The Dantzig method is simple but its use is deprecated.  Steepest is the method of choice and there\
 are two variants which keep all weights updated but only scan a subset each iteration.\
 Partial switches this on while automatic decides at each iteration based on information\
 about the factorization.\
 The PE variants add the Positive Edge criterion. \
 This selects incoming variables to try to avoid degenerate moves. See also option psi.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("dualS!implex", "Do dual simplex algorithm",
      CLP_PARAM_ACTION_DUALSIMPLEX);
    p.setLonghelp(
      "This command solves the continuous relaxation of the current model using the dual steepest edge algorithm.\
The time and iterations may be affected by settings such as presolve, scaling, crash\
 and also by dual pivot method, fake bound on variables and dual and primal tolerances.");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("dualT!olerance", "For an optimal solution \
no dual infeasibility may exceed this value",
      1.0e-20, COIN_DBL_MAX, CLP_PARAM_DBL_DUALTOLERANCE);
    p.setLonghelp(
      "Normally the default tolerance is fine, but one may want to increase it a\
 bit if the dual simplex algorithm seems to be having a hard time.  One method which can be faster is \
to use a large tolerance e.g. 1.0e-4 and the dual simplex algorithm and then to clean up the problem using the primal simplex algorithm with the \
correct tolerance (remembering to switch off presolve for this final short clean up phase).");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("dw!Heuristic", "Whether to try Dantzig Wolfe heuristic",
      "off", CBC_PARAM_STR_DW);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(
      "This heuristic is very very compute intensive. It tries to find a Dantzig Wolfe structure and use that. "
      HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("either!Simplex", "Do dual or primal simplex algorithm",
      CLP_PARAM_ACTION_EITHERSIMPLEX);
    p.setLonghelp(
      "This command solves the continuous relaxation of the current model using the dual or primal algorithm,\
 based on a dubious analysis of model.");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("end", "Stops clp execution",
      CLP_PARAM_ACTION_EXIT);
    p.setLonghelp(
      "This stops execution ; end, exit, quit and stop are synonyms");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("environ!ment", "Read commands from environment",
      CLP_PARAM_ACTION_ENVIRONMENT, 7, 0);
    p.setLonghelp(
      "This starts reading from environment variable CBC_CLP_ENVIRONMENT.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("error!sAllowed", "Whether to allow import errors",
      "off", CLP_PARAM_STR_ERRORSALLOWED, 3);

    p.append("on");
    p.setLonghelp(
      "The default is not to use any model which had errors when reading the mps file.\
  Setting this to 'on' will allow all errors from which the code can recover\
 simply by ignoring the error.  There are some errors from which the code can not recover \
e.g. no ENDATA.  This has to be set before import i.e. -errorsAllowed on -import xxxxxx.mps.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("exit", "Stops clp execution",
      CLP_PARAM_ACTION_EXIT);
    p.setLonghelp(
      "This stops the execution of Clp, end, exit, quit and stop are synonyms");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("exper!iment", "Whether to use testing features",
      -1, 200000, CBC_PARAM_INT_EXPERIMENT, 0);
    p.setLonghelp(
      "Defines how adventurous you want to be in using new ideas. \
0 then no new ideas, 1 fairly sensible, 2 a bit dubious, 3 you are on your own!");

    p.setIntValue(0);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("expensive!Strong", "Whether to do even more strong branching",
      0, COIN_INT_MAX, CBC_PARAM_INT_STRONG_STRATEGY, 0);

    p.setLonghelp(
      "Strategy for extra strong branching. \
0 is normal strong branching. \
1, 2, 4, and 6 does strong branching on all fractional variables if \
at the root node (1), \
at depth less than modifier (2), \
objective equals best possible (4), or \
at depth less than modifier and objective equals best possible (6). \
11, 12, 14, and 16 are like 1, 2, 4, and 6, respectively, but do strong branching on all integer (incl. non-fractional) variables. \
Values >= 100 are used to specify a depth limit (value/100), otherwise 5 is used. \
If the values >= 100, then above rules are applied to value%100.");
    p.setIntValue(0);
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("export", "Export model as mps file",
      CLP_PARAM_ACTION_EXPORT);

    p.setLonghelp(
      "This will write an MPS format file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.mps'.  \
It can be useful to get rid of the original names and go over to using Rnnnnnnn and Cnnnnnnn.  This can be done by setting 'keepnames' off before importing mps file.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("extra1", "Extra integer parameter 1",
      -COIN_INT_MAX, COIN_INT_MAX, CBC_PARAM_INT_EXTRA1, 0);
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("extra2", "Extra integer parameter 2",
      -COIN_INT_MAX, COIN_INT_MAX, CBC_PARAM_INT_EXTRA2, 0);
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("extra3", "Extra integer parameter 3",
      -COIN_INT_MAX, COIN_INT_MAX, CBC_PARAM_INT_EXTRA3, 0);
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("extra4", "Extra integer parameter 4",
      -1, COIN_INT_MAX, CBC_PARAM_INT_EXTRA4, 0);

    p.setIntValue(-1);
    p.setLonghelp(
      "This switches on yet more special options!! \
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
    CbcOrClpParam p("extraV!ariables", "Allow creation of extra integer variables",
      -COIN_INT_MAX, COIN_INT_MAX, CBC_PARAM_INT_EXTRA_VARIABLES, 0);
    p.setIntValue(0);
    p.setLonghelp(
      "Switches on a trivial re-formulation that introduces extra integer variables to group together variables with same cost.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("fact!orization", "Which factorization to use",
      "normal", CLP_PARAM_STR_FACTORIZATION);
    p.append("dense");
    p.append("simple");
    p.append("osl");
    p.setLonghelp(
      "The default is to use the normal CoinFactorization, but \
other choices are a dense one, OSL's, or one designed for small problems."
    );
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("fakeB!ound", "All bounds <= this value - DEBUG",
      1.0, 1.0e15, CLP_PARAM_ACTION_FAKEBOUND, 0);
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("feas!ibilityPump", "Whether to try the Feasibility Pump heuristic",
      "off", CBC_PARAM_STR_FPUMP);

    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(
      "This heuristic is due to Fischetti, Glover, and Lodi \
and uses a sequence of LPs to try and get an integer feasible solution. \
Some fine tuning is available by options passFeasibilityPump and pumpTune. "
      HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("fix!OnDj", "Try heuristic based on fixing variables with \
reduced costs greater than this",
      -COIN_DBL_MAX, COIN_DBL_MAX, CBC_PARAM_DBL_DJFIX, 1);
    p.setLonghelp(
      "If this is set integer variables with reduced costs greater than this will be fixed \
before branch and bound - use with extreme caution!");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("flow!CoverCuts", "Whether to use Flow Cover cuts",
      "off", CBC_PARAM_STR_FLOWCUTS);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.setFakeKeyWord(3);
    p.setLonghelp(CUTS_LONGHELP
      " Reference: https://github.com/coin-or/Cgl/wiki/CglFlowCover"); // Can also enter testing values by plusnn (==ifmove)
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("force!Solution", "Whether to use given solution as crash for BAB",
      -1, 20000000, CLP_PARAM_INT_USESOLUTION);
    p.setIntValue(-1);
    p.setLonghelp(
      "-1 off.  If 1 then tries to branch to solution given by AMPL or priorities file. \
If 0 then just tries to set as best solution \
If >1 then also does that many nodes on fixed problem.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("fraction!forBAB", "Fraction in feasibility pump",
      1.0e-5, 1.1, CBC_PARAM_DBL_SMALLBAB, 1);
    p.setDoubleValue(0.5);
    p.setLonghelp(
      "After a pass in the feasibility pump, variables which have not moved \
about are fixed and if the preprocessed model is smaller than this fraction of the original problem, \
a few nodes of branch and bound are done on the reduced problem.");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("gamma!(Delta)", "Whether to regularize barrier",
      "off", CLP_PARAM_STR_GAMMA, 7, 1);
    p.append("on");
    p.append("gamma");
    p.append("delta");
    p.append("onstrong");
    p.append("gammastrong");
    p.append("deltastrong");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("GMI!Cuts", "Whether to use alternative Gomory cuts",
      "off", CBC_PARAM_STR_GMICUTS);

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
    p.setLonghelp(CUTS_LONGHELP
      " This version is by Giacomo Nannicini and may be more robust than gomoryCuts.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("gomory!Cuts", "Whether to use Gomory cuts",
      "off", CBC_PARAM_STR_GOMORYCUTS);

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
so it is worth experimenting (Long allows any length). "
    CUTS_LONGHELP
    " Reference: https://github.com/coin-or/Cgl/wiki/CglGomory");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("greedy!Heuristic", "Whether to use a greedy heuristic",
      "off", CBC_PARAM_STR_GREEDY);

    p.append("on");
    p.append("both");
    p.append("before");
    //p.append("root");
    p.setLonghelp(
      "This heuristic tries to obtain a feasible solution by just fixing a percentage of variables and then try a small branch and cut run. "
      HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("gsolu!tion", "Puts glpk solution to file",
      CLP_PARAM_ACTION_GMPL_SOLUTION);

    p.setLonghelp(
      "Will write a glpk solution file to the given file name.  It will use the default \
directory given by 'directory'.  A name of '$' will use the previous value for the \
name.  This is initialized to 'stdout' (this defaults to ordinary solution if stdout). \
If problem created from gmpl model - will do any reports.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("guess", "Guesses at good parameters", CLP_PARAM_ACTION_GUESS, 7);
    p.setLonghelp(
    "This looks at model statistics and does an initial solve \
setting some parameters which may help you to think of possibilities.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("heur!isticsOnOff", "Switches most primal heuristics on or off",
      "off", CBC_PARAM_STR_HEURISTICSTRATEGY);
    p.append("on");
    p.setLonghelp(
      "This option can be used to switch on or off all heuristics that search for feasible solutions,\
      except for the local tree search, as it dramatically alters the search.\
      Then individual heuristics can be turned off or on.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("clique!Cuts", "Whether to use Clique cuts",
      "off", CBC_PARAM_STR_CLIQUECUTS);
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
    CbcOrClpParam p("cgraph",
        "Whether to use the conflict graph-based preprocessing and cut separation routines.",
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
    CbcOrClpParam p("bkpivot!ing", "Pivoting strategy used in Bron-Kerbosch algorithm",
                        0, 6, CBC_PARAM_INT_BKPIVOTINGSTRATEGY);
    p.setIntValue(3);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("bkmaxcalls", "Maximum number of recursive calls made by Bron-Kerbosch algorithm",
                        1, 2147483647, CBC_PARAM_INT_BKMAXCALLS);
    p.setIntValue(1000);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("bkclqext!method",
        "Strategy used to extend violated cliques found by BK Clique Cut Separation routine",
        0, 5, CBC_PARAM_INT_BKCLQEXTMETHOD);
    p.setLonghelp(
        "Sets the method used in the extension module of BK Clique Cut Separation routine: \
0=no extension; 1=random; 2=degree; 3=modified degree; 4=reduced cost(inversely proportional); 5=reduced cost(inversely proportional) + modified degree");
    p.setIntValue(4);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("oddwheel!Cuts", "Whether to use odd wheel cuts",
                      "off", CBC_PARAM_STR_ODDWHEELCUTS);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.setLonghelp("This switches on odd-wheel inequalities (either at root or in entire tree).");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p(
        "oddwext!method", "Strategy used to search for wheel centers for the cuts found by \
Odd Wheel Cut Separation routine", 0, 2, CBC_PARAM_INT_ODDWEXTMETHOD);
    p.setLonghelp(
        "Sets the method used in the extension module of Odd Wheel Cut Separation routine: \
0=no extension; 1=one variable; 2=clique");
    p.setIntValue(2);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("clqstr!engthen", "Whether to perform Clique Strengthening preprocessing routine",
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
    CbcOrClpParam p("help", "Print out version, non-standard options and some help",
      CLP_PARAM_ACTION_HELP, 3);
    p.setLonghelp(
      "This prints out some help to get user started.  If you have printed this then \
you should be past that stage:-)");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("hOp!tions", "Heuristic options",
      -COIN_INT_MAX, COIN_INT_MAX, CBC_PARAM_INT_HOPTIONS, 1);
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
    CbcOrClpParam p("hot!StartMaxIts", "Maximum iterations on hot start",
      0, COIN_INT_MAX, CBC_PARAM_INT_MAXHOTITS);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("idiot!Crash", "Whether to try idiot crash",
      -1, COIN_INT_MAX, CLP_PARAM_INT_IDIOT);

    p.setLonghelp(
      "This is a type of 'crash' which works well on some homogeneous problems.\
 It works best on problems with unit elements and rhs but will do something to any \
 model.  It should only be used before the primal simplex algorithm.  It can be set to -1 when the code \
 decides for itself whether to use it, 0 to switch off, or n > 0 to do n passes.");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("import", "Import model from mps file",
      CLP_PARAM_ACTION_IMPORT, 3);
    p.setLonghelp(
      "This will read an MPS format file from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  If you have libgz then it can read compressed\
 files 'xxxxxxxx.gz' or 'xxxxxxxx.bz2'.  \
If 'keepnames' is off, then names are dropped -> Rnnnnnnn and Cnnnnnnn.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("inc!rement", "A valid solution must be at least this \
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
    CbcOrClpParam p("inf!easibilityWeight", "Each integer infeasibility is expected \
to cost this much",
      0.0, COIN_DBL_MAX, CBC_PARAM_DBL_INFEASIBILITYWEIGHT, 1);
    p.setLonghelp(
      "A primitive way of deciding which node to explore next.  Satisfying each integer infeasibility is \
expected to cost this much.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("initialS!olve", "Solve to continuous",
      CLP_PARAM_ACTION_SOLVECONTINUOUS);

    p.setLonghelp(
      "This just solves the problem to continuous - without adding any cuts");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("integerT!olerance", "For a feasible solution \
no integer variable may be more than this away from an integer value",
      1.0e-20, 0.5, CBC_PARAM_DBL_INTEGERTOLERANCE);
    p.setLonghelp(
      "Beware of setting this smaller than the primal feasibility tolerance.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("keepN!ames", "Whether to keep names from import",
      "on", CLP_PARAM_STR_KEEPNAMES);
    p.append("off");
    p.setLonghelp(
      "It saves space to get rid of names so if you need to you can set this to off.  \
This needs to be set before the import of model - so -keepnames off -import xxxxx.mps.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("KKT", "Whether to use KKT factorization in barrier",
      "off", CLP_PARAM_STR_KKT, 7, 1);
    p.append("on");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("knapsack!Cuts", "Whether to use Knapsack cuts",
      "off", CBC_PARAM_STR_KNAPSACKCUTS);

    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.append("forceandglobal");
    p.setLonghelp(CUTS_LONGHELP
      " Reference: https://github.com/coin-or/Cgl/wiki/CglKnapsackCover");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("lagomory!Cuts", "Whether to use Lagrangean Gomory cuts",
      "off", CBC_PARAM_STR_LAGOMORYCUTS);
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
    CbcOrClpParam p("latwomir!Cuts", "Whether to use Lagrangean TwoMir cuts",
      "off", CBC_PARAM_STR_LATWOMIRCUTS);

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
    p.setLonghelp(
      "This is a Lagrangean relaxation for TwoMir cuts.  See \
  lagomoryCuts for description of options.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("lift!AndProjectCuts", "Whether to use Lift and Project cuts",
      "off", CBC_PARAM_STR_LANDPCUTS);

    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("iflongon");
    p.setLonghelp(
      "These cuts may be expensive to compute. "
      CUTS_LONGHELP
      " Reference: https://github.com/coin-or/Cgl/wiki/CglLandP");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("local!TreeSearch", "Whether to use local tree search when a solution is found",
      "off", CBC_PARAM_STR_LOCALTREE);
    p.append("on");
    p.setLonghelp(
      "The heuristic is from Fischetti and Lodi and is not really a heuristic although it can be used as one \
(with limited functionality).  It is not switched on when heuristics are switched on.");
    parameters.push_back(p);
  }
#endif
  {
#ifndef COIN_HAS_CBC
    CbcOrClpParam p("log!Level", "Level of detail in Solver output",
      -1, COIN_INT_MAX, CLP_PARAM_INT_SOLVERLOGLEVEL);
    parameters.push_back(p);
#else
    CbcOrClpParam p("log!Level", "Level of detail in Coin branch and Cut output",
      -63, 63, CLP_PARAM_INT_LOGLEVEL);
    p.setIntValue(1);
#endif
    p.setLonghelp(
      "If 0 then there should be no output in normal circumstances.  1 is probably the best\
 value for most uses, while 2 and 3 give more information.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("max!imize", "Set optimization direction to maximize",
      CLP_PARAM_ACTION_MAXIMIZE, 7);
    p.setLonghelp(
      "The default is minimize - use 'maximize' for maximization.\n\
You can also use the parameters 'direction maximize'.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("maxF!actor", "Maximum number of iterations between \
refactorizations",
      1, COIN_INT_MAX, CLP_PARAM_INT_MAXFACTOR);
    p.setLonghelp(
      "If this is left at its default value of 200 then CLP will guess a\
 value to use.  CLP may decide to re-factorize earlier for accuracy.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("maxIt!erations", "Maximum number of iterations before \
stopping",
      0, COIN_INT_MAX, CLP_PARAM_INT_MAXITERATION);
    p.setLonghelp(
      "This can be used for testing purposes.  The corresponding library call\n\
      \tsetMaximumIterations(value)\n can be useful.  If the code stops on\
 seconds or by an interrupt this will be treated as stopping on maximum iterations.  This is ignored in branchAndCut - use maxN!odes.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("maxN!odes", "Maximum number of nodes to do",
      -1, COIN_INT_MAX, CBC_PARAM_INT_MAXNODES);
    p.setLonghelp(
      "This is a repeatable way to limit search.  Normally using time is easier \
but then the results may not be repeatable.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("maxNI!FS", "Maximum number of nodes to be processed without improving the incumbent solution.",
      -1, COIN_INT_MAX, CBC_PARAM_INT_MAXNODESNOTIMPROVINGFS);
    p.setLonghelp(
      "This criterion specifies that when a feasible solution is available, the search should continue\
only if better feasible solutions were produced in the last nodes.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("maxSaved!Solutions", "Maximum number of solutions to save",
      0, COIN_INT_MAX, CBC_PARAM_INT_MAXSAVEDSOLS);
    p.setLonghelp(
      "Number of solutions to save.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("maxSo!lutions", "Maximum number of feasible solutions to get",
      1, COIN_INT_MAX, CBC_PARAM_INT_MAXSOLS);
    p.setLonghelp(
      "You may want to stop after (say) two solutions or an hour.  \
This is checked every node in tree, so it is possible to get more solutions from heuristics.");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("min!imize", "Set optimization direction to minimize",
      CLP_PARAM_ACTION_MINIMIZE, 7);
    p.setLonghelp(
      "The default is minimize - use 'maximize' for maximization.\n\
This should only be necessary if you have previously set maximization \
You can also use the parameters 'direction minimize'.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("mipO!ptions", "Dubious options for mip",
      0, COIN_INT_MAX, CBC_PARAM_INT_MIPOPTIONS, 0);
    p.setIntValue(1057);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("more!MipOptions", "More dubious options for mip",
      -1, COIN_INT_MAX, CBC_PARAM_INT_MOREMIPOPTIONS, 0);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("more2!MipOptions", "More more dubious options for mip",
      -1, COIN_INT_MAX, CBC_PARAM_INT_MOREMOREMIPOPTIONS, 0);
    p.appendStringValue("nodezero1[8192][More strong branching at root node]#");
    p.appendStringValue("nodezero2[16384][More strong branching at root node - more]#");
    p.appendStringValue("nodezero3[24576][More strong branching at root node - yet more]#");
    p.appendStringValue("lagrangean1[234881024][lagrangean cuts at end of root cuts]#");
    p.appendStringValue("lagrangean2[268435456][lagrangean cuts at end of root cuts -alt]#");
    p.appendStringValue("lessused[536870912][less used cuts at beginning of root cuts]#");
    p.appendStringValue("#+"); // + allowed
    p.setIntValue(0);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("mixed!IntegerRoundingCuts", "Whether to use Mixed Integer Rounding cuts",
      "off", CBC_PARAM_STR_MIXEDCUTS);

    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.setLonghelp(CUTS_LONGHELP
      " Reference: https://github.com/coin-or/Cgl/wiki/CglMixedIntegerRounding2");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("mess!ages", "Controls if Clpnnnn is printed",
      "off", CLP_PARAM_STR_MESSAGES);

    p.append("on");
    p.setLonghelp("The default behavior is to put out messages such as:\n\
   Clp0005 2261  Objective 109.024 Primal infeas 944413 (758)\n\
but this program turns this off to make it look more friendly.  It can be useful\
 to turn them back on if you want to be able to 'grep' for particular messages or if\
 you intend to override the behavior of a particular message.  This only affects Clp not Cbc.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("miplib", "Do some of miplib test set",
      CBC_PARAM_ACTION_MIPLIB, 3, 1);
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("mips!tart", "reads an initial feasible solution from file",
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
#endif
  {
    CbcOrClpParam p("moreS!pecialOptions", "Yet more dubious options for Simplex - see ClpSimplex.hpp",
		    0, COIN_INT_MAX, CLP_PARAM_INT_MORESPECIALOPTIONS, 0);
    p.appendStringValue("keep!DualOrPrimal[8192][If you ask for dual you will always get dual (and for primal)]#");
    p.appendStringValue("clean!Scaled[134217728][Make sure unscaled problem is feasible if scaled problem is feasible]#");
    p.appendStringValue("#+"); // + allowed
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("moreT!une", "Yet more dubious ideas for feasibility pump",
      0, 100000000, CBC_PARAM_INT_FPUMPTUNE2, 0);

    p.setLonghelp(
      "Yet more ideas for Feasibility Pump \n\
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
    CbcOrClpParam p("multiple!RootPasses", "Do multiple root passes to collect cuts and solutions",
      0, COIN_INT_MAX, CBC_PARAM_INT_MULTIPLEROOTS, 0);
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
    CbcOrClpParam p("naive!Heuristics", "Whether to try some stupid heuristic",
      "off", CBC_PARAM_STR_NAIVE, 7, 1);

    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(
      "This is naive heuristics which, e.g., fix all integers with costs to zero!. "
      HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("netlib", "Solve entire netlib test set",
      CLP_PARAM_ACTION_NETLIB_EITHER, 3, 1);
    p.setLonghelp(
      "This exercises the unit test for clp and then solves the netlib test set using dual or primal.\
The user can set options before e.g. clp -presolve off -netlib");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("netlibB!arrier", "Solve entire netlib test set with barrier",
      CLP_PARAM_ACTION_NETLIB_BARRIER, 3, 1);
    p.setLonghelp(
      "This exercises the unit test for clp and then solves the netlib test set using barrier.\
The user can set options before e.g. clp -kkt on -netlib");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("netlibD!ual", "Solve entire netlib test set (dual)",
      CLP_PARAM_ACTION_NETLIB_DUAL, 3, 1);

    p.setLonghelp(
      "This exercises the unit test for clp and then solves the netlib test set using dual.\
The user can set options before e.g. clp -presolve off -netlib");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("netlibP!rimal", "Solve entire netlib test set (primal)",
      CLP_PARAM_ACTION_NETLIB_PRIMAL, 3, 1);
    p.setLonghelp(
      "This exercises the unit test for clp and then solves the netlib test set using primal.\
The user can set options before e.g. clp -presolve off -netlibp");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("netlibT!une", "Solve entire netlib test set with 'best' algorithm",
      CLP_PARAM_ACTION_NETLIB_TUNE, 3, 1);
    p.setLonghelp(
      "This exercises the unit test for clp and then solves the netlib test set using whatever \
works best.  I know this is cheating but it also stresses the code better by doing a \
mixture of stuff.  The best algorithm was chosen on a Linux ThinkPad using native cholesky \
with University of Florida ordering.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("network", "Tries to make network matrix",
      CLP_PARAM_ACTION_NETWORK, 7, 0);
    p.setLonghelp(
      "Clp will go faster if the matrix can be converted to a network.  The matrix\
 operations may be a bit faster with more efficient storage, but the main advantage\
 comes from using a network factorization.  It will probably not be as fast as a \
specialized network code.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("nextB!estSolution", "Prints next best saved solution to file",
      CLP_PARAM_ACTION_NEXTBESTSOLUTION);

    p.setLonghelp(
      "To write best solution, just use solution.  This prints next best (if exists) \
 and then deletes it. \
 This will write a primitive solution file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'stdout'.  The amount of output can be varied using printi!ngOptions or printMask.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("node!Strategy", "What strategy to use to select the next node from the branch and cut tree",
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
    CbcOrClpParam p("numberA!nalyze", "Number of analysis iterations",
      -COIN_INT_MAX, COIN_INT_MAX, CBC_PARAM_INT_NUMBERANALYZE, 0);
    p.setLonghelp(
      "This says how many iterations to spend at root node analyzing problem. \
This is a first try and will hopefully become more sophisticated.");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("objective!Scale", "Scale factor to apply to objective",
      -COIN_DBL_MAX, COIN_DBL_MAX, CLP_PARAM_DBL_OBJSCALE, 1);
    p.setLonghelp(
      "If the objective function has some very large values, you may wish to scale them\
 internally by this amount.  It can also be set by autoscale.  It is applied after scaling.  You are unlikely to need this.");
    p.setDoubleValue(1.0);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("PrepN!ames", "If column names will be kept in pre-processed model",
      "off", CBC_PARAM_STR_PREPROCNAMES);
    p.append("on");
    p.setLonghelp(
      "Normally the preprocessed model has column names replaced by new names C0000...\
Setting this option to on keeps original names in variables which still exist in the preprocessed problem");
    parameters.push_back(p);
  }

  {
    CbcOrClpParam p("outDup!licates", "takes duplicate rows etc out of integer model",
      CLP_PARAM_ACTION_OUTDUPROWS, 7, 0);
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("output!Format", "Which output format to use",
      1, 6, CLP_PARAM_INT_OUTPUTFORMAT);
    p.setLonghelp(
      "Normally export will be done using normal representation for numbers and two values\
 per line.  You may want to do just one per line (for grep or suchlike) and you may wish\
 to save with absolute accuracy using a coded version of the IEEE value. A value of 2 is normal.\
 otherwise odd values gives one value per line, even two.  Values 1,2 give normal format, 3,4\
 gives greater precision, while 5,6 give IEEE values.  When used for exporting a basis 1 does not save \
values, 2 saves values, 3 with greater accuracy and 4 in IEEE.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("para!metrics", "Import data from file and do parametrics",
      CLP_PARAM_ACTION_PARAMETRICS, 3);
    p.setLonghelp(
      "This will read a file with parametric data from the given file name \
and then do parametrics.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  This can not read from compressed files. \
File is in modified csv format - a line ROWS will be followed by rows data \
while a line COLUMNS will be followed by column data.  The last line \
should be ENDATA. The ROWS line must exist and is in the format \
ROWS, initial theta, final theta, interval theta, n where n is 0 to get \
CLPI0062 message at interval or at each change of theta \
and 1 to get CLPI0063 message at each iteration.  If interval theta is 0.0 \
or >= final theta then no interval reporting.  n may be missed out when it is \
taken as 0.  If there is Row data then \
there is a headings line with allowed headings - name, number, \
lower(rhs change), upper(rhs change), rhs(change).  Either the lower and upper \
fields should be given or the rhs field. \
The optional COLUMNS line is followed by a headings line with allowed \
headings - name, number, objective(change), lower(change), upper(change). \
 Exactly one of name and number must be given for either section and \
missing ones have value 0.0.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("passC!uts", "Number of rounds that cut generators are applied in the root node",
      -COIN_INT_MAX, COIN_INT_MAX, CBC_PARAM_INT_CUTPASS);

    p.setIntValue(20);
    p.setLonghelp(
      "The default is to do 100 passes if the problem has less than 500 columns, 100 passes (but \
stop if the drop in the objective function value is small) if the problem has less than 5000 columns, and 20 passes otherwise. \
A negative value -n means that n passes are also applied if the objective does not drop.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("passF!easibilityPump", "How many passes to do in the Feasibility Pump heuristic",
      0, 10000, CBC_PARAM_INT_FPUMPITS);
    p.setIntValue(20);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("passP!resolve", "How many passes in presolve",
      -200, 100, CLP_PARAM_INT_PRESOLVEPASS, 1);
    p.setLonghelp(
      "Normally Presolve does 10 passes but you may want to do less to make it\
 more lightweight or do more if improvements are still being made.  As Presolve will return\
 if nothing is being taken out, you should not normally need to use this fine tuning.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("passT!reeCuts", "Number of rounds that cut generators are applied in the tree",
      -COIN_INT_MAX, COIN_INT_MAX, CBC_PARAM_INT_CUTPASSINTREE);
    p.setIntValue(1);
    p.setLonghelp("The default is to do one pass. A negative value -n means that n passes are also applied if the objective does not drop.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("pertV!alue", "Method of perturbation",
      -5000, 102, CLP_PARAM_INT_PERTVALUE, 1);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("perturb!ation", "Whether to perturb the problem",
      "on", CLP_PARAM_STR_PERTURBATION);
    p.append("off");
    p.setLonghelp(
      "Perturbation helps to stop cycling, but CLP uses other measures for this.\
  However, large problems and especially ones with unit elements and unit right hand sides or costs\
 benefit from perturbation.  Normally CLP tries to be intelligent, but one can switch this off.");
  // The Clp library has this off by default.  This program has it on by default.
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("PFI", "Whether to use Product Form of Inverse in simplex",
      "off", CLP_PARAM_STR_PFI, 7, 0);
    p.append("on");
    p.setLonghelp(
      "By default clp uses Forrest-Tomlin L-U update.  If you are masochistic you can switch it off.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("pivotAndC!omplement", "Whether to try Pivot and Complement heuristic",
      "off", CBC_PARAM_STR_PIVOTANDCOMPLEMENT);

    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("pivotAndF!ix", "Whether to try Pivot and Fix heuristic",
      "off", CBC_PARAM_STR_PIVOTANDFIX);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("plus!Minus", "Tries to make +- 1 matrix",
      CLP_PARAM_ACTION_PLUSMINUS, 7, 0);
    p.setLonghelp(
      "Clp will go slightly faster if the matrix can be converted so that the elements are\
 not stored and are known to be unit.  The main advantage is memory use.  Clp may automatically\
 see if it can convert the problem so you should not need to use this.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("pO!ptions", "Dubious print options",
      0, COIN_INT_MAX, CLP_PARAM_INT_PRINTOPTIONS, 1);
    p.setIntValue(0);
    p.setLonghelp(
      "If this is > 0 then presolve will give more information and branch and cut will give statistics");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("preO!pt", "Presolve options",
      0, COIN_INT_MAX, CLP_PARAM_INT_PRESOLVEOPTIONS, 0);
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("presolve", "Whether to presolve problem",
      "on", CLP_PARAM_STR_PRESOLVE);
    p.append("off");
    p.append("more");
    p.append("file");
    p.setLonghelp("Presolve analyzes the model to find such things as redundant equations, equations\
 which fix some variables, equations which can be transformed into bounds, etc.  For the\
 initial solve of any problem this is worth doing unless one knows that it will have no effect.  \
Option 'on' will normally do 5 passes, while using 'more' will do 10.  If the problem is very large one can \
let CLP write the original problem to file by using 'file'.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("preprocess", "Whether to use integer preprocessing",
      "off", CBC_PARAM_STR_PREPROCESS);

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
    p.append("equalallstop");
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
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("preT!olerance", "Tolerance to use in presolve",
      1.0e-20, COIN_DBL_MAX, CLP_PARAM_DBL_PRESOLVETOLERANCE);
    p.setLonghelp(
      "One may want to increase this tolerance if presolve says the problem is \
infeasible and one has awkward numbers and is sure that the problem is really feasible.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("primalP!ivot", "Primal pivot choice algorithm",
      "auto!matic", CLP_PARAM_STR_PRIMALPIVOT, 7, 1);

    p.append("exa!ct");
    p.append("dant!zig");
    p.append("part!ial");
    p.append("steep!est");
    p.append("change");
    p.append("sprint");
    p.append("PEsteep!est");
    p.append("PEdantzig");
    p.setLonghelp(
      "The Dantzig method is simple but its use is deprecated.  Exact devex is the method of choice and there\
 are two variants which keep all weights updated but only scan a subset each iteration.\
 Partial switches this on while 'change' initially does 'dantzig' until the factorization\
 becomes denser.  This is still a work in progress.\
 The PE variants add the Positive Edge criterion.\
 This selects incoming variables to try to avoid degenerate moves. \
 See also Towhidi, M., Desrosiers, J., Soumis, F., The positive edge criterion within COIN-OR's CLP;\
 Omer, J., Towhidi, M., Soumis, F., The positive edge pricing rule for the dual simplex.");

    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("primalS!implex", "Do primal simplex algorithm",
      CLP_PARAM_ACTION_PRIMALSIMPLEX);

    p.setLonghelp(
      "This command solves the continuous relaxation of the current model using the primal algorithm.\
  The default is to use exact devex.\
 The time and iterations may be affected by settings such as presolve, scaling, crash\
 and also by column selection  method, infeasibility weight and dual and primal tolerances.");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("primalT!olerance", "For a feasible solution \
no primal infeasibility, i.e., constraint violation, may exceed this value",
      1.0e-20, COIN_DBL_MAX, CLP_PARAM_DBL_PRIMALTOLERANCE);
    p.setLonghelp(
      "Normally the default tolerance is fine, but one may want to increase it a\
 bit if the primal simplex algorithm seems to be having a hard time.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("primalW!eight", "Initially algorithm acts as if it \
costs this much to be infeasible",
      1.0e-20, COIN_DBL_MAX, CLP_PARAM_DBL_PRIMALWEIGHT);
    p.setLonghelp(
      "The primal algorithm in Clp is a single phase algorithm as opposed to a two phase\
 algorithm where you first get feasible then optimal.  So Clp is minimizing this weight times\
 the sum of primal infeasibilities plus the true objective function (in minimization sense).\
  Too high a value may mean more iterations, while too low a value means\
 the algorithm may iterate into the wrong directory for long and then has to increase the weight in order to get feasible."); // OSL had a heuristic to adjust bounds, maybe we need that here.
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("psi", "Two-dimension pricing factor for Positive Edge criterion",
      -1.1, 1.1, CLP_PARAM_DBL_PSI);

    p.setDoubleValue(-0.5);
    p.setLonghelp(
      "The Positive Edge criterion has been added to \
select incoming variables to try and avoid degenerate moves. \
Variables not in the promising set have their infeasibility weight multiplied by psi, \
so 0.01 would mean that if there were any promising variables, then they would always be chosen, \
while 1.0 effectively switches the algorithm off. \
There are two ways of switching this feature on.  One way is to set psi to a positive value and then \
the Positive Edge criterion will be used for both primal and dual simplex.  The other way is to select PEsteepest \
in dualpivot choice (for example), then the absolute value of psi is used. \
Code donated by Jeremy Omer.  See \
Towhidi, M., Desrosiers, J., Soumis, F., The positive edge criterion within COIN-OR's CLP; \
Omer, J., Towhidi, M., Soumis, F., The positive edge pricing rule for the dual simplex.");  // Until this settles down it is only implemented in CLP.
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("printi!ngOptions", "Print options",
      "normal", CLP_PARAM_STR_INTPRINT, 3);
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
    CbcOrClpParam p("printM!ask", "Control printing of solution on a  mask",
      CLP_PARAM_ACTION_PRINTMASK, 3);

    p.setLonghelp(
      "If set then only those names which match mask are printed in a solution. \
'?' matches any character and '*' matches any set of characters. \
 The default is '' i.e. unset so all variables are printed. \
This is only active if model has names.");
    parameters.push_back(p);
  }

#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("prio!rityIn", "Import priorities etc from file",
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
    CbcOrClpParam p("probing!Cuts", "Whether to use Probing cuts",
      "off", CBC_PARAM_STR_PROBINGCUTS);
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
    p.setLonghelp(CUTS_LONGHELP
      " Value 'forceOnBut' turns on probing and forces CBC to do probing at every node, but does only probing, not strengthening etc. \
    Value 'strong' forces CBC to strongly do probing at every node, that is, also when CBC would usually turn it off because it hasn't found something. \
    Value 'forceonbutstrong' is like 'forceonstrong', but does only probing (column fixing) and turns off row strengthening, so the matrix will not change inside the branch and bound. \
    Reference: https://github.com/coin-or/Cgl/wiki/CglProbing");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("progress!(Interval)", "Time interval for printing progress",
		    -COIN_DBL_MAX,COIN_DBL_MAX,
		    CLP_PARAM_DBL_PROGRESS);
    p.setLonghelp(
      "This sets a minimum interval for some printing - elapsed seconds");

    p.setDoubleValue(0.7);
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("proximity!Search", "Whether to do proximity search heuristic",
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
found at the root node. "
      HEURISTICS_LONGHELP); // Can also set different maxNode settings by plusnnnn (and are 'on'(on==30)).
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("pumpC!utoff", "Fake cutoff for use in feasibility pump",
      -COIN_DBL_MAX, COIN_DBL_MAX, CBC_PARAM_DBL_FAKECUTOFF);
    p.setDoubleValue(0.0);
    p.setLonghelp(
      "A value of 0.0 means off. Otherwise, add a constraint forcing objective below this value\
 in feasibility pump");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("pumpI!ncrement", "Fake increment for use in feasibility pump",
      -COIN_DBL_MAX, COIN_DBL_MAX, CBC_PARAM_DBL_FAKEINCREMENT, 1);
    p.setDoubleValue(0.0);
    p.setLonghelp(
      "A value of 0.0 means off. Otherwise use as absolute increment to cutoff \
when solution found in feasibility pump");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("pumpT!une", "Dubious ideas for feasibility pump",
      0, 100000000, CBC_PARAM_INT_FPUMPTUNE);
    p.setIntValue(1003);
    p.setLonghelp(
      "This fine tunes Feasibility Pump \n\
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
#endif
  {
    CbcOrClpParam p("quit", "Stops clp execution",
      CLP_PARAM_ACTION_EXIT);
    p.setLonghelp(
      "This stops the execution of Clp, end, exit, quit and stop are synonyms");

    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("randomC!bcSeed", "Random seed for Cbc",
      -1, COIN_INT_MAX, CBC_PARAM_INT_RANDOMSEED);

    p.setLonghelp(
      "Allows initialization of the random seed for pseudo-random numbers used in heuristics such as the Feasibility Pump to decide whether to round up or down. "
      "The special value of 0 lets Cbc use the time of the day for the initial seed.");
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("randomi!zedRounding", "Whether to try randomized rounding heuristic",
      "off", CBC_PARAM_STR_RANDROUND);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("randomS!eed", "Random seed for Clp",
      0, COIN_INT_MAX, CLP_PARAM_INT_RANDOMSEED);

    p.setLonghelp(
      "Initialization of the random seed for pseudo-random numbers used to break ties in degenerate problems. "
      "This may yield a different continuous optimum and, in the context of Cbc, different cuts and heuristic solutions. "
      "The special value of 0 lets CLP use the time of the day for the initial seed.");
    p.setIntValue(1234567);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("ratio!Gap", "Stop when gap between best possible and \
best known is less than this fraction of larger of two",
      0.0, COIN_DBL_MAX, CBC_PARAM_DBL_GAPRATIO);
    p.setDoubleValue(1e-4);
    p.setLonghelp(
      "If the gap between the best known solution and the best possible solution is less than this fraction \
of the objective value at the root node then the search will terminate.  See 'allowableGap' for a \
way of using absolute value rather than fraction.");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("restoreS!olution", "reads solution from file",
      CLP_PARAM_ACTION_RESTORESOL);

    p.setLonghelp(
      "This will read a binary solution file from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'solution.file'.  This reads in a file from saveSolution");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("readSt!ored", "Import stored cuts from file",
      CLP_PARAM_ACTION_STOREDFILE, 3, 0);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("reallyO!bjectiveScale", "Scale factor to apply to objective in place",
      -COIN_DBL_MAX, COIN_DBL_MAX, CLP_PARAM_DBL_OBJSCALE2, 0);
    p.setLonghelp("You can set this to -1.0 to test maximization or other to stress code");
    p.setDoubleValue(1.0);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("reallyS!cale", "Scales model in place",
      CLP_PARAM_ACTION_REALLY_SCALE, 7, 0);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("reduce!AndSplitCuts", "Whether to use Reduce-and-Split cuts",
      "off", CBC_PARAM_STR_REDSPLITCUTS);

    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.setLonghelp(
      "These cuts may be expensive to generate. "
      CUTS_LONGHELP
      " Reference: https://github.com/coin-or/Cgl/wiki/CglRedSplit");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("reduce2!AndSplitCuts", "Whether to use Reduce-and-Split cuts - style 2",
      "off", CBC_PARAM_STR_REDSPLIT2CUTS);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("longOn");
    p.append("longRoot");
    p.setLonghelp("This switches on reduce and split  cuts (either at root or in entire tree). \
This version is by Giacomo Nannicini based on Francois Margot's version. \
Standard setting only uses rows in tableau <= 256, long uses all. \
These cuts may be expensive to generate. \
See option cuts for more information on the possible values.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("residual!CapacityCuts", "Whether to use Residual Capacity cuts",
      "off", CBC_PARAM_STR_RESIDCUTS);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.setLonghelp(CUTS_LONGHELP
      " Reference: https://github.com/coin-or/Cgl/wiki/CglResidualCapacity");

    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("restore!Model", "Restore model from binary file",
      CLP_PARAM_ACTION_RESTORE, 7, 1);
    p.setLonghelp(
      "This reads data save by saveModel from the given file.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.prob'.");

    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("reverse", "Reverses sign of objective",
      CLP_PARAM_ACTION_REVERSE, 7, 0);
    p.setLonghelp(
      "Useful for testing if maximization works correctly");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("rhs!Scale", "Scale factor to apply to rhs and bounds",
      -COIN_DBL_MAX, COIN_DBL_MAX, CLP_PARAM_DBL_RHSSCALE, 0);
    p.setLonghelp(
      "If the rhs or bounds have some very large meaningful values, you may wish to scale them\
 internally by this amount.  It can also be set by autoscale.  This should not be needed.");
    p.setDoubleValue(1.0);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("Rens", "Whether to try Relaxation Enforced Neighborhood Search",
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
    p.setLonghelp(HEURISTICS_LONGHELP
      " Value 'on' just does 50 nodes. 200, 1000, and 10000 does that many nodes.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("Rins", "Whether to try Relaxed Induced Neighborhood Search",
      "off", CBC_PARAM_STR_RINS);
    p.append("on");
    p.append("both");
    p.append("before");
    p.append("often");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("round!ingHeuristic", "Whether to use simple (but effective) Rounding heuristic",
      "off", CBC_PARAM_STR_ROUNDING);
    p.append("on");
    p.append("both");
    p.append("before");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("saveM!odel", "Save model to binary file",
      CLP_PARAM_ACTION_SAVE, 7, 1);
    p.setLonghelp(
      "This will save the problem to the given file name for future use\
 by restoreModel.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.prob'.");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("saveS!olution", "saves solution to file",
      CLP_PARAM_ACTION_SAVESOL);

    p.setLonghelp(
      "This will write a binary solution file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'solution.file'.  To read the file use fread(int) twice to pick up number of rows \
and columns, then fread(double) to pick up objective value, then pick up row activities, row duals, column \
activities and reduced costs - see bottom of CbcOrClpParam.cpp for code that reads or writes file. \
If name contains '_fix_read_' then does not write but reads and will fix all variables");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("scal!ing", "Whether to scale problem",
      "off", CLP_PARAM_STR_SCALING);
    p.append("equi!librium");
    p.append("geo!metric");
    p.append("auto!matic");
    p.append("dynamic");
    p.append("rows!only");
    p.setLonghelp(
      "Scaling can help in solving problems which might otherwise fail because of lack of\
 accuracy.  It can also reduce the number of iterations.  It is not applied if the range\
 of elements is small.  When the solution is evaluated in the unscaled problem, it is possible that small primal and/or\
 dual infeasibilities occur. "
 "Option 'equilibrium' uses the largest element for scaling. "
 "Option 'geometric' uses the squareroot of the product of largest and smallest element. "
 "Option 'auto' let CLP choose a method that gives the best ratio of the largest element to the smallest one.");
    p.setCurrentOption(3); // say auto
    parameters.push_back(p);
  }
#ifndef COIN_HAS_CBC
  {
    CbcOrClpParam p("sec!onds", "Maximum seconds",
      -1.0, COIN_DBL_MAX, CLP_PARAM_DBL_TIMELIMIT);

    p.setLonghelp("After this many seconds clp will act as if maximum iterations had been reached \
(if value >=0).");
    parameters.push_back(p);
  }
#else
  {
    CbcOrClpParam p("sec!onds", "maximum seconds",
      -1.0, COIN_DBL_MAX, CBC_PARAM_DBL_TIMELIMIT_BAB);
  // Meaning 0 - start at very beginning
  // 1 start at beginning of preprocessing
  // 2 start at beginning of branch and bound
#ifndef CBC_USE_INITIAL_TIME
#define CBC_USE_INITIAL_TIME 1
#endif
#if CBC_USE_INITIAL_TIME==0
    p.setLonghelp(
      "After this many seconds in Branch and Bound coin solver will act as if maximum nodes had been reached (time in initial solve and preprocessing is included).");
#elif CBC_USE_INITIAL_TIME==1
    p.setLonghelp(
      "After this many seconds in Branch and Bound coin solver will act as if maximum nodes had been reached (time in initial solve not included).");
#else
    p.setLonghelp(
      "After this many seconds in Branch and Bound coin solver will act as if maximum nodes had been reached (time in initial solve and preprocessing not included).");
#endif
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("secni!fs", "maximum seconds without improving the incumbent solution",
    -1.0, COIN_DBL_MAX, CBC_PARAM_DBL_MAXSECONDSNIFS);
    p.setLonghelp(
      "With this stopping criterion, after a feasible solution is found, the search should continue only if the incumbent solution was updated recently, \
the tolerance is specified here. A discussion on why this criterion can be useful is included here: \
https://yetanothermathprogrammingconsultant.blogspot.com/2019/11/mip-solver-stopping-criteria.html .");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("sleep", "for debug",
      CLP_PARAM_ACTION_DUMMY, 7, 0);

    p.setLonghelp(
      "If passed to solver fom ampl, then ampl will wait so that you can copy .nl file for debug.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("slow!cutpasses", "Maximum number of rounds for slower cut generators",
      -1, COIN_INT_MAX, CBC_PARAM_INT_MAX_SLOW_CUTS);
    p.setLonghelp(
      "Some cut generators are fairly slow - this limits the number of times they are tried.\
      The cut generators identified as 'may be slow' at present are Lift and project cuts and both versions of Reduce and Split cuts.");
    p.setIntValue(10);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("slp!Value", "Number of slp passes before primal",
      -50000, 50000, CLP_PARAM_INT_SLPVALUE, 1);
    p.setLonghelp(
      "If you are solving a quadratic problem using primal then it may be helpful to do some \
sequential Lps to get a good approximate solution.");
    parameters.push_back(p);
  }
#if CLP_MULTIPLE_FACTORIZATIONS > 0
  {
    CbcOrClpParam p("small!Factorization", "Threshold for using small factorization",
      -1, 10000, CBC_PARAM_INT_SMALLFACT, 1);
    p.setLonghelp(
      "If processed problem <= this use small factorization");
    p.setIntValue(-1);
    parameters.push_back(p);
  }
#endif
#endif
  {
    CbcOrClpParam p("solu!tion", "Prints solution to file",
      CLP_PARAM_ACTION_SOLUTION);
    p.setLonghelp(
      "This will write a primitive solution file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'stdout'.  The amount of output can be varied using printi!ngOptions or printMask.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CLP
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("solv!e", "Solve problem",
      CBC_PARAM_ACTION_BAB);
         p.setLonghelp(
          "If there are no integer variables then this just solves LP.  If there are integer variables \
this does branch and cut." );
         parameters.push_back( p );
  }
  {
    CbcOrClpParam p("sosO!ptions", "Whether to use SOS from AMPL",  "off", CBC_PARAM_STR_SOS);
    p.append("on");
    p.setCurrentOption("on");
         p.setLonghelp(
          "Normally if AMPL says there are SOS variables they should be used, but sometime sthey should\
 be turned off - this does so." );
         parameters.push_back( p );
  }
  {
    CbcOrClpParam p("slog!Level", "Level of detail in (LP) Solver output", -1, 63, CLP_PARAM_INT_SOLVERLOGLEVEL);
    p.setLonghelp(
      "If 0 then there should be no output in normal circumstances.  1 is probably the best\
 value for most uses, while 2 and 3 give more information.  This parameter is only used inside MIP - for Clp use 'log'");
    parameters.push_back(p);
  }
  {
     // Due to James Howey
     CbcOrClpParam p("sosP!rioritize", "How to deal with SOS priorities",
       "off", CBC_PARAM_STR_SOSPRIORITIZE);
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
#else
  // allow solve as synonym for possible dual
  {
    CbcOrClpParam p("solv!e", "Solve problem using dual simplex (probably)",
      CLP_PARAM_ACTION_EITHERSIMPLEX);
    p.setLonghelp(
      "Just so can use solve for clp as well as in cbc");
    parameters.push_back(p);
  }
#endif
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("spars!eFactor", "Whether factorization treated as sparse",
      "on", CLP_PARAM_STR_SPARSEFACTOR, 7, 0);
    p.append("off");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("special!Options", "Dubious options for Simplex - see ClpSimplex.hpp",
      0, COIN_INT_MAX, CLP_PARAM_INT_SPECIALOPTIONS, 0);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("sprint!Crash", "Whether to try sprint crash",
      -1, COIN_INT_MAX, CLP_PARAM_INT_SPRINT);
    p.setLonghelp(
      "For long and thin problems this method may solve a series of small problems\
 created by taking a subset of the columns.  The idea as 'Sprint' was introduced by J. Forrest after\
 an LP code of that name of the 60's which tried the same tactic (not totally successfully).\
  CPLEX calls it 'sifting'.  -1 lets CLP automatically choose the number of passes, 0 is off, n is number of passes");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("stat!istics", "Print some statistics",
      CLP_PARAM_ACTION_STATISTICS);
    p.setLonghelp(
      "This command prints some statistics for the current model.\
 If log level >1 then more is printed.\
 These are for presolved model if presolve on (and unscaled).");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("stop", "Stops clp execution",
      CLP_PARAM_ACTION_EXIT);
    p.setLonghelp(
      "This stops the execution of Clp, end, exit, quit and stop are synonyms");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("strat!egy", "Switches on groups of features",
      0, 2, CBC_PARAM_INT_STRATEGY);
    p.setLonghelp(
      "This turns on newer features. \
Use 0 for easy problems, 1 is default, 2 is aggressive. \
1 uses Gomory cuts with a tolerance of 0.01 at the root node, \
does a possible restart after 100 nodes if many variables could be fixed, \
activates a diving and RINS heuristic, and makes the feasibility pump \
more aggressive."); // This does not apply to unit tests (where 'experiment' may have similar effects
    p.setIntValue(1);
    parameters.push_back(p);
  }
#ifdef CBC_KEEP_DEPRECATED
  {
    CbcOrClpParam p("strengthen", "Create strengthened problem",
      CBC_PARAM_ACTION_STRENGTHEN, 3);
    p.setLonghelp(
      "This creates a new problem by applying the root node cuts.  All tight constraints \
will be in resulting problem");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("strong!Branching", "Number of variables to look at in strong branching",
      0, COIN_INT_MAX, CBC_PARAM_INT_STRONGBRANCHING);
    p.setIntValue(20);
    p.setLonghelp(
      "In order to decide which variable to branch on, the code will choose up to this number \
of unsatisfied variables to try minimal up and down branches on.  Then the most effective one is chosen. \
If a variable is branched on many times then the previous average up and down costs may be used - \
see also option trustPseudoCosts.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("subs!titution", "How long a column to substitute for in presolve",
      0, 10000, CLP_PARAM_INT_SUBSTITUTION, 0);
    p.setLonghelp(
      "Normally Presolve gets rid of 'free' variables when there are no more than 3 \
 coefficients in a row.  If you increase this, the number of rows may decrease but the number of \
 coefficients may increase.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("testO!si", "Test OsiObject stuff",
      -1, COIN_INT_MAX, CBC_PARAM_INT_TESTOSI, 0);
    parameters.push_back(p);
  }
#endif
#ifdef CBC_THREAD
  {
    CbcOrClpParam p("thread!s", "Number of threads to try and use",
      -100, 100000, CBC_PARAM_INT_THREADS, 1);
    p.setIntValue(0);
    p.setLonghelp(
      "To use multiple threads, set threads to number wanted.  It may be better \
to use one or two more than number of cpus available.  If 100+n then n threads and \
search is repeatable (maybe be somewhat slower), \
if 200+n use threads for root cuts, 400+n threads used in sub-trees.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("tighten!Factor", "Tighten bounds using this times largest \
activity at continuous solution",
      1.0e-3, COIN_DBL_MAX, CBC_PARAM_DBL_TIGHTENFACTOR, 0);
    p.setLonghelp(
      "This sleazy trick can help on some problems.");
    parameters.push_back(p);
  }

#endif
#ifdef COIN_HAS_CLP
  {
    CbcOrClpParam p("tightLP", "Poor person's preSolve for now",
      CLP_PARAM_ACTION_TIGHTEN, 7, 0);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("trust!PseudoCosts", "Number of branches before we trust pseudocosts",
      -3, COIN_INT_MAX, CBC_PARAM_INT_NUMBERBEFORE);
    p.setLonghelp(
      "Using strong branching computes pseudo-costs.  This parameter determines after how many branches for a variable we just \
trust the pseudo costs and do not do any more strong branching.");
    p.setIntValue(10);
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("tune!PreProcess", "Dubious tuning parameters for preprocessing",
      0, COIN_INT_MAX, CLP_PARAM_INT_PROCESSTUNE, 1);
    p.appendStringValue("heavy!Probing[7][Do more probing]#");
    p.appendStringValue("heavier!Probing[519][Do yet more probing]#");
    p.appendStringValue("#="); // = allowed (so only one)
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
    CbcOrClpParam p("two!MirCuts", "Whether to use Two phase Mixed Integer Rounding cuts",
      "off", CBC_PARAM_STR_TWOMIRCUTS);
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
#endif
  {
    CbcOrClpParam p("unitTest", "Do unit test",
      CLP_PARAM_ACTION_UNITTEST, 3, 1);
    p.setLonghelp(
      "This exercises the unit test for clp");
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("userClp", "Hand coded Clp stuff",
      CLP_PARAM_ACTION_USERCLP, 0, 0);
    p.setLonghelp(
      "There are times e.g. when using AMPL interface when you may wish to do something unusual.  \
Look for USERCLP in main driver and modify sample code.");
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("userCbc", "Hand coded Cbc stuff",
      CBC_PARAM_ACTION_USERCBC, 0, 0);
    p.setLonghelp(
      "There are times e.g. when using AMPL interface when you may wish to do something unusual.  \
Look for USERCBC in main driver and modify sample code. \
It is possible you can get same effect by using example driver4.cpp.");
    parameters.push_back(p);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("Vnd!VariableNeighborhoodSearch", "Whether to try Variable Neighborhood Search",
      "off", CBC_PARAM_STR_VND);
    p.append("on");
    p.append("both");
    p.append("before");
    p.append("intree");
    p.setLonghelp(HEURISTICS_LONGHELP);
    parameters.push_back(p);
  }
#endif
  {
#ifndef COIN_AVX2
    CbcOrClpParam p("vector", "Whether to use vector? Form of matrix in simplex",
      "off", CLP_PARAM_STR_VECTOR, 7, 0);
    p.append("on");
    p.setLonghelp(
      "If this is on ClpPackedMatrix uses extra column copy in odd format.");
    parameters.push_back(p);
#else
    CbcOrClpParam p("vector", "Try and use vector instructions in simplex",
      "off", CLP_PARAM_STR_VECTOR, 7, 0);
    p.append("on");
    p.append("ones");
    p.setLonghelp(
      "At present only for Intel architectures - but could be extended.  \
Uses avx2 or avx512 instructions. Uses different storage for matrix - can be \
of benefit without instruction set on some problems.  \
I may add pool to switch on a pool matrix");
    parameters.push_back(p);
#endif
  }
  {
    CbcOrClpParam p("verbose", "Switches on longer help on single ?",
      0, 31, CLP_PARAM_INT_VERBOSE, 0);
    p.setLonghelp(
      "Set to 1 to get short help with ? list, 2 to get long help, 3 for both.  (add 4 to just get ampl ones).");
    p.setIntValue(0);
    parameters.push_back(p);
  }
#ifdef COIN_HAS_CBC
  {
    CbcOrClpParam p("vub!heuristic", "Type of VUB heuristic",
      -2, 20, CBC_PARAM_INT_VUBTRY, 0);
    p.setLonghelp(
      "This heuristic tries and fix some integer variables.");
    p.setIntValue(-1);
    parameters.push_back(p);
  }
  {
    CbcOrClpParam p("zero!HalfCuts", "Whether to use zero half cuts",
      "off", CBC_PARAM_STR_ZEROHALFCUTS);
    p.append("on");
    p.append("root");
    p.append("ifmove");
    p.append("forceOn");
    p.append("onglobal");
    p.setLonghelp(CUTS_LONGHELP
      " This implementation was written by Alberto Caprara.");
    parameters.push_back(p);
  }
#endif
  {
    CbcOrClpParam p("zeroT!olerance", "Kill all coefficients \
whose absolute value is less than this value",
      1.0e-100, 1.0e-5, CLP_PARAM_DBL_ZEROTOLERANCE);
    p.setLonghelp(
      "This applies to reading mps files (and also lp files \
if KILL_ZERO_READLP defined)");
    p.setDoubleValue(1.0e-20);
    parameters.push_back(p);
  }
}
// Given a parameter type - returns its number in list
int whichParam(const CbcOrClpParameterType &name,
  const std::vector< CbcOrClpParam > &parameters)
{
  for (int i = 0; i < (int)parameters.size(); i++) {
    if (parameters[i].type() == name)
      return i;
  }
  return std::numeric_limits< int >::max(); // should not arrive here
}
#ifdef COIN_HAS_CLP
/* Restore a solution from file.
   mode 0 normal, 1 swap rows and columns and primal and dual
   if 2 set then also change signs
*/
void restoreSolution(ClpSimplex *lpSolver, std::string fileName, int mode)
{
  FILE *fp = fopen(fileName.c_str(), "rb");
  if (fp) {
    int numberRows = lpSolver->numberRows();
    int numberColumns = lpSolver->numberColumns();
    int numberRowsFile;
    int numberColumnsFile;
    double objectiveValue;
    size_t nRead;
    nRead = fread(&numberRowsFile, sizeof(int), 1, fp);
    if (nRead != 1)
      throw("Error in fread");
    nRead = fread(&numberColumnsFile, sizeof(int), 1, fp);
    if (nRead != 1)
      throw("Error in fread");
    nRead = fread(&objectiveValue, sizeof(double), 1, fp);
    if (nRead != 1)
      throw("Error in fread");
    double *dualRowSolution = lpSolver->dualRowSolution();
    double *primalRowSolution = lpSolver->primalRowSolution();
    double *dualColumnSolution = lpSolver->dualColumnSolution();
    double *primalColumnSolution = lpSolver->primalColumnSolution();
    if (mode) {
      // swap
      int k = numberRows;
      numberRows = numberColumns;
      numberColumns = k;
      double *temp;
      temp = dualRowSolution;
      dualRowSolution = primalColumnSolution;
      primalColumnSolution = temp;
      temp = dualColumnSolution;
      dualColumnSolution = primalRowSolution;
      primalRowSolution = temp;
    }
    if (numberRows > numberRowsFile || numberColumns > numberColumnsFile) {
      std::cout << "Mismatch on rows and/or columns - giving up" << std::endl;
    } else {
      lpSolver->setObjectiveValue(objectiveValue);
      if (numberRows == numberRowsFile && numberColumns == numberColumnsFile) {
        nRead = fread(primalRowSolution, sizeof(double), numberRows, fp);
        if (nRead != static_cast< size_t >(numberRows))
          throw("Error in fread");
        nRead = fread(dualRowSolution, sizeof(double), numberRows, fp);
        if (nRead != static_cast< size_t >(numberRows))
          throw("Error in fread");
        nRead = fread(primalColumnSolution, sizeof(double), numberColumns, fp);
        if (nRead != static_cast< size_t >(numberColumns))
          throw("Error in fread");
        nRead = fread(dualColumnSolution, sizeof(double), numberColumns, fp);
        if (nRead != static_cast< size_t >(numberColumns))
          throw("Error in fread");
      } else {
        std::cout << "Mismatch on rows and/or columns - truncating" << std::endl;
        double *temp = new double[std::max(numberRowsFile, numberColumnsFile)];
        nRead = fread(temp, sizeof(double), numberRowsFile, fp);
        if (nRead != static_cast< size_t >(numberRowsFile))
          throw("Error in fread");
        CoinMemcpyN(temp, numberRows, primalRowSolution);
        nRead = fread(temp, sizeof(double), numberRowsFile, fp);
        if (nRead != static_cast< size_t >(numberRowsFile))
          throw("Error in fread");
        CoinMemcpyN(temp, numberRows, dualRowSolution);
        nRead = fread(temp, sizeof(double), numberColumnsFile, fp);
        if (nRead != static_cast< size_t >(numberColumnsFile))
          throw("Error in fread");
        CoinMemcpyN(temp, numberColumns, primalColumnSolution);
        nRead = fread(temp, sizeof(double), numberColumnsFile, fp);
        if (nRead != static_cast< size_t >(numberColumnsFile))
          throw("Error in fread");
        CoinMemcpyN(temp, numberColumns, dualColumnSolution);
        delete[] temp;
      }
      if (mode == 3) {
        int i;
        for (i = 0; i < numberRows; i++) {
          primalRowSolution[i] = -primalRowSolution[i];
          dualRowSolution[i] = -dualRowSolution[i];
        }
        for (i = 0; i < numberColumns; i++) {
          primalColumnSolution[i] = -primalColumnSolution[i];
          dualColumnSolution[i] = -dualColumnSolution[i];
        }
      }
    }
    fclose(fp);
  } else {
    std::cout << "Unable to open file " << fileName << std::endl;
  }
}
// Dump a solution to file
void saveSolution(const ClpSimplex *lpSolver, std::string fileName)
{
  if (strstr(fileName.c_str(), "_fix_read_")) {
    FILE *fp = fopen(fileName.c_str(), "rb");
    if (fp) {
      ClpSimplex *solver = const_cast< ClpSimplex * >(lpSolver);
      restoreSolution(solver, fileName, 0);
      // fix all
      int logLevel = solver->logLevel();
      int iColumn;
      int numberColumns = solver->numberColumns();
      double *primalColumnSolution = solver->primalColumnSolution();
      double *columnLower = solver->columnLower();
      double *columnUpper = solver->columnUpper();
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = primalColumnSolution[iColumn];
        if (value > columnUpper[iColumn]) {
          if (value > columnUpper[iColumn] + 1.0e-6 && logLevel > 1)
            printf("%d value of %g - bounds %g %g\n",
              iColumn, value, columnLower[iColumn], columnUpper[iColumn]);
          value = columnUpper[iColumn];
        } else if (value < columnLower[iColumn]) {
          if (value < columnLower[iColumn] - 1.0e-6 && logLevel > 1)
            printf("%d value of %g - bounds %g %g\n",
              iColumn, value, columnLower[iColumn], columnUpper[iColumn]);
          value = columnLower[iColumn];
        }
        columnLower[iColumn] = value;
        columnUpper[iColumn] = value;
      }
      return;
    }
  }
  FILE *fp = fopen(fileName.c_str(), "wb");
  if (fp) {
    int numberRows = lpSolver->numberRows();
    int numberColumns = lpSolver->numberColumns();
    double objectiveValue = lpSolver->objectiveValue();
    size_t nWrite;
    nWrite = fwrite(&numberRows, sizeof(int), 1, fp);
    if (nWrite != 1)
      throw("Error in fwrite");
    nWrite = fwrite(&numberColumns, sizeof(int), 1, fp);
    if (nWrite != 1)
      throw("Error in fwrite");
    nWrite = fwrite(&objectiveValue, sizeof(double), 1, fp);
    if (nWrite != 1)
      throw("Error in fwrite");
    double *dualRowSolution = lpSolver->dualRowSolution();
    double *primalRowSolution = lpSolver->primalRowSolution();
    nWrite = fwrite(primalRowSolution, sizeof(double), numberRows, fp);
    if (nWrite != static_cast< size_t >(numberRows))
      throw("Error in fwrite");
    nWrite = fwrite(dualRowSolution, sizeof(double), numberRows, fp);
    if (nWrite != static_cast< size_t >(numberRows))
      throw("Error in fwrite");
    double *dualColumnSolution = lpSolver->dualColumnSolution();
    double *primalColumnSolution = lpSolver->primalColumnSolution();
    nWrite = fwrite(primalColumnSolution, sizeof(double), numberColumns, fp);
    if (nWrite != static_cast< size_t >(numberColumns))
      throw("Error in fwrite");
    nWrite = fwrite(dualColumnSolution, sizeof(double), numberColumns, fp);
    if (nWrite != static_cast< size_t >(numberColumns))
      throw("Error in fwrite");
    fclose(fp);
  } else {
    std::cout << "Unable to open file " << fileName << std::endl;
  }
}
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
