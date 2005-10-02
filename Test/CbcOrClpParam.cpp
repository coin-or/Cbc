// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <string>
#include <iostream>
#include <cassert>

#include "CbcOrClpParam.hpp"
#ifdef COIN_USE_CBC
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#endif
#ifdef COIN_USE_CLP
#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"
#endif
#ifdef COIN_USE_READLINE     
#include <readline/readline.h>
#include <readline/history.h>
#endif

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CbcOrClpParam::CbcOrClpParam () 
  : type_(INVALID),
    lowerDoubleValue_(0.0),
    upperDoubleValue_(0.0),
    lowerIntValue_(0),
    upperIntValue_(0),
    lengthName_(0),
    lengthMatch_(0),
    definedKeyWords_(),
    name_(),
    shortHelp_(),
    longHelp_(),
    action_(INVALID),
    currentKeyWord_(-1),
    display_(false),
    intValue_(-1),
    doubleValue_(-1.0),
    stringValue_(""),
    indexNumber_(INVALID)
{
}
// Other constructors
CbcOrClpParam::CbcOrClpParam (std::string name, std::string help,
	   double lower, double upper, CbcOrClpParameterType type,
		    bool display)
  : type_(type),
    lowerIntValue_(0),
    upperIntValue_(0),
    definedKeyWords_(),
    name_(name),
    shortHelp_(help),
    longHelp_(),
    action_(type),
    currentKeyWord_(-1),
    display_(display),
    intValue_(-1),
    doubleValue_(-1.0),
    stringValue_(""),
    indexNumber_(type)
{
  lowerDoubleValue_ = lower;
  upperDoubleValue_ = upper;
  gutsOfConstructor();
}
CbcOrClpParam::CbcOrClpParam (std::string name, std::string help,
	   int lower, int upper, CbcOrClpParameterType type,
		    bool display)
  : type_(type),
    lowerDoubleValue_(0.0),
    upperDoubleValue_(0.0),
    definedKeyWords_(),
    name_(name),
    shortHelp_(help),
    longHelp_(),
    action_(type),
    currentKeyWord_(-1),
    display_(display),
    intValue_(-1),
    doubleValue_(-1.0),
    stringValue_(""),
    indexNumber_(type)
{
  gutsOfConstructor();
  lowerIntValue_ = lower;
  upperIntValue_ = upper;
}
// Other strings will be added by append
CbcOrClpParam::CbcOrClpParam (std::string name, std::string help, 
		    std::string firstValue,
		    CbcOrClpParameterType type,int defaultIndex,
		    bool display)
  : type_(type),
    lowerDoubleValue_(0.0),
    upperDoubleValue_(0.0),
    lowerIntValue_(0),
    upperIntValue_(0),
    definedKeyWords_(),
    name_(name),
    shortHelp_(help),
    longHelp_(),
    action_(type),
    currentKeyWord_(defaultIndex),
    display_(display),
    intValue_(-1),
    doubleValue_(-1.0),
    stringValue_(""),
    indexNumber_(type)
{
  gutsOfConstructor();
  definedKeyWords_.push_back(firstValue);
}
// Action
CbcOrClpParam::CbcOrClpParam (std::string name, std::string help,
		    CbcOrClpParameterType type,int indexNumber,
		    bool display)
  : type_(type),
    lowerDoubleValue_(0.0),
    upperDoubleValue_(0.0),
    lowerIntValue_(0),
    upperIntValue_(0),
    definedKeyWords_(),
    name_(name),
    shortHelp_(help),
    longHelp_(),
    action_(type),
    currentKeyWord_(-1),
    display_(display),
    intValue_(-1),
    doubleValue_(-1.0),
    stringValue_("")
{
  if (indexNumber<0)
    indexNumber_=type;
  else
    indexNumber_=indexNumber;
  gutsOfConstructor();
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CbcOrClpParam::CbcOrClpParam (const CbcOrClpParam & rhs) 
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
  display_=rhs.display_;
  intValue_=rhs.intValue_;
  doubleValue_=rhs.doubleValue_;
  stringValue_=rhs.stringValue_;
  indexNumber_=rhs.indexNumber_;
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CbcOrClpParam::~CbcOrClpParam ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CbcOrClpParam &
CbcOrClpParam::operator=(const CbcOrClpParam& rhs)
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
    display_=rhs.display_;
    intValue_=rhs.intValue_;
    doubleValue_=rhs.doubleValue_;
    stringValue_=rhs.stringValue_;
    indexNumber_=rhs.indexNumber_;
  }
  return *this;
}
void 
CbcOrClpParam::gutsOfConstructor()
{
  std::string::size_type  shriekPos = name_.find('!');
  lengthName_ = name_.length();
  if ( shriekPos==std::string::npos ) {
    //does not contain '!'
    lengthMatch_= lengthName_;
  } else {
    lengthMatch_=shriekPos;
    name_ = name_.substr(0,shriekPos)+name_.substr(shriekPos+1);
    lengthName_--;
  }
}
// Insert string (only valid for keywords)
void 
CbcOrClpParam::append(std::string keyWord)
{
  definedKeyWords_.push_back(keyWord);
}

int 
CbcOrClpParam::matches (std::string input) const
{
  // look up strings to do more elegantly
  if (input.length()>lengthName_) {
    return 0;
  } else {
    unsigned int i;
    for (i=0;i<input.length();i++) {
      if (tolower(name_[i])!=tolower(input[i])) 
	break;
    }
    if (i<input.length()) {
      return 0;
    } else if (i>=lengthMatch_) {
      return 1;
    } else {
      // matched but too short
      return 2;
    }
  }
}
// Returns name which could match
std::string 
CbcOrClpParam::matchName (  ) const
{ 
  if (lengthMatch_==lengthName_) 
    return name_;
  else
    return name_.substr(0,lengthMatch_)+"("+name_.substr(lengthMatch_)+")";
}

// Returns parameter option which matches (-1 if none)
int 
CbcOrClpParam::parameterOption ( std::string check ) const
{
  int numberItems = definedKeyWords_.size();
  if (!numberItems) {
    return -1;
  } else {
    int whichItem=0;
    unsigned int it;
    for (it=0;it<definedKeyWords_.size();it++) {
      std::string thisOne = definedKeyWords_[it];
      std::string::size_type  shriekPos = thisOne.find('!');
      unsigned int length1 = thisOne.length();
      unsigned int length2 = length1;
      if ( shriekPos!=std::string::npos ) {
	//contains '!'
	length2 = shriekPos;
	thisOne = thisOne.substr(0,shriekPos)+
	  thisOne.substr(shriekPos+1);
	length1 = thisOne.length();
      }
      if (check.length()<=length1&&length2<=check.length()) {
	unsigned int i;
	for (i=0;i<check.length();i++) {
	  if (tolower(thisOne[i])!=tolower(check[i])) 
	    break;
	}
	if (i<check.length()) {
	  whichItem++;
	} else if (i>=length2) {
	  break;
	} 
      } else {
	whichItem++;
      }
    }
    if (whichItem<numberItems)
      return whichItem;
    else
      return -1;
  }
}
// Prints parameter options
void 
CbcOrClpParam::printOptions (  ) const
{
  std::cout<<"Possible options for "<<name_<<" are:"<<std::endl;
  unsigned int it;
  for (it=0;it<definedKeyWords_.size();it++) {
    std::string thisOne = definedKeyWords_[it];
    std::string::size_type  shriekPos = thisOne.find('!');
    if ( shriekPos!=std::string::npos ) {
      //contains '!'
      thisOne = thisOne.substr(0,shriekPos)+
	"("+thisOne.substr(shriekPos+1)+")";
    }
    std::cout<<thisOne<<std::endl;
  }
}
// Print action and string
void 
CbcOrClpParam::printString() const
{
  if (name_=="directory")
    std::cout<<"Current working directory is "<<stringValue_<<std::endl;
  else
    std::cout<<"Current default (if $ as parameter) for "<<name_
	     <<" is "<<stringValue_<<std::endl;
}
void CoinReadPrintit(const char * input)
{
  int length =strlen(input);
  char temp[101];
  int i;
  int n=0;
  for (i=0;i<length;i++) {
    if (input[i]=='\n') {
      temp[n]='\0';
      std::cout<<temp<<std::endl;
      n=0;
    } else if (n>=65&&input[i]==' ') {
      temp[n]='\0';
      std::cout<<temp<<std::endl;
      n=0;
    } else if (n||input[i]!=' ') {
      temp[n++]=input[i];
    }
  }
  if (n) {
    temp[n]='\0';
    std::cout<<temp<<std::endl;
  }
}
// Print Long help
void 
CbcOrClpParam::printLongHelp() const
{
  if (type_>=1&&type_<400) {
    if (type_<LOGLEVEL) {
      printf("Range of values is %g to %g\n",lowerDoubleValue_,upperDoubleValue_);
    } else if (type_<DIRECTION) {
      printf("Range of values is %d to %d\n",lowerIntValue_,upperIntValue_);
    } else if (type_<DIRECTORY) {
      printOptions();
    }
    CoinReadPrintit(longHelp_.c_str());
  }
}
#ifdef COIN_USE_CBC
int
CbcOrClpParam::setDoubleParameter (OsiSolverInterface * model,double value) 
{
  if (value<lowerDoubleValue_||value>upperDoubleValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerDoubleValue_<<" to "<<
      upperDoubleValue_<<std::endl;
    return 1;
  } else {
    double oldValue;
    switch(type_) {
    case DUALTOLERANCE:
      model->getDblParam(OsiDualTolerance,oldValue);
      model->setDblParam(OsiDualTolerance,value);
      break;
    case PRIMALTOLERANCE:
      model->getDblParam(OsiPrimalTolerance,oldValue);
      model->setDblParam(OsiPrimalTolerance,value);
      break;
    default:
      oldValue=0.0; // to avoid compiler message
      abort();
    }
    std::cout<<name_<<" was changed from "<<oldValue<<" to "
	     <<value<<std::endl;
    return 0;
  }
}
#endif
#ifdef COIN_USE_CLP
int
CbcOrClpParam::setDoubleParameter (ClpSimplex * model,double value) 
{
  double oldValue = doubleParameter(model);
  if (value<lowerDoubleValue_||value>upperDoubleValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerDoubleValue_<<" to "<<
      upperDoubleValue_<<std::endl;
    return 1;
  } else {
    std::cout<<name_<<" was changed from "<<oldValue<<" to "
	     <<value<<std::endl;
    switch(type_) {
#ifndef COIN_USE_CBC
    case DUALTOLERANCE:
      model->setDualTolerance(value);
      break;
    case PRIMALTOLERANCE:
      model->setPrimalTolerance(value);
      break;
#endif
    case DUALBOUND:
      model->setDualBound(value);
      break;
    case PRIMALWEIGHT:
      model->setInfeasibilityCost(value);
      break;
    case TIMELIMIT:
      model->setMaximumSeconds(value);
      break;
    case OBJSCALE:
      model->setObjectiveScale(value);
      break;
    case RHSSCALE:
      model->setRhsScale(value);
      break;
    default:
      abort();
    }
    return 0;
  }
}
double 
CbcOrClpParam::doubleParameter (ClpSimplex * model) const
{
  double value;
  switch(type_) {
#ifndef COIN_USE_CBC
  case DUALTOLERANCE:
    value=model->dualTolerance();
    break;
  case PRIMALTOLERANCE:
    value=model->primalTolerance();
    break;
#endif
  case DUALBOUND:
    value=model->dualBound();
    break;
  case PRIMALWEIGHT:
    value=model->infeasibilityCost();
    break;
  case TIMELIMIT:
    value=model->maximumSeconds();
    break;
  case OBJSCALE:
    value=model->objectiveScale();
    break;
  case RHSSCALE:
    value=model->rhsScale();
    break;
  default:
    abort();
  }
  return value;
}
int 
CbcOrClpParam::setIntParameter (ClpSimplex * model,int value) 
{
  int oldValue = intParameter(model);
  if (value<lowerIntValue_||value>upperIntValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerIntValue_<<" to "<<
      upperIntValue_<<std::endl;
    return 1;
  } else {
    std::cout<<name_<<" was changed from "<<oldValue<<" to "
	     <<value<<std::endl;
    switch(type_) {
    case SOLVERLOGLEVEL:
      model->setLogLevel(value);
      if (value>2)
	model->factorization()->messageLevel(8);
      else
	model->factorization()->messageLevel(0);
      break;
    case MAXFACTOR:
      model->factorization()->maximumPivots(value);
      break;
    case PERTVALUE:
      model->setPerturbation(value);
      break;
    case MAXITERATION:
      model->setMaximumIterations(value);
      break;
    case SPECIALOPTIONS:
      model->setSpecialOptions(value);
      break;
    default:
      abort();
    }
    return 0;
  }
}
int 
CbcOrClpParam::intParameter (ClpSimplex * model) const
{
  int value;
  switch(type_) {
#ifndef COIN_USE_CBC
  case SOLVERLOGLEVEL:
    value=model->logLevel();
    break;
#endif
  case MAXFACTOR:
    value=model->factorization()->maximumPivots();
    break;
    break;
  case PERTVALUE:
    value=model->perturbation();
    break;
  case MAXITERATION:
    value=model->maximumIterations();
    break;
  case SPECIALOPTIONS:
    value=model->specialOptions();
    break;
  default:
    value=-1;
    break;
  }
  return value;
}
#endif
int
CbcOrClpParam::checkDoubleParameter (double value) const
{
  if (value<lowerDoubleValue_||value>upperDoubleValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerDoubleValue_<<" to "<<
      upperDoubleValue_<<std::endl;
    return 1;
  } else {
    return 0;
  }
}
#ifdef COIN_USE_CBC
double 
CbcOrClpParam::doubleParameter (OsiSolverInterface * model) const
{
  double value;
  switch(type_) {
  case DUALTOLERANCE:
    assert(model->getDblParam(OsiDualTolerance,value));
    break;
  case PRIMALTOLERANCE:
    assert(model->getDblParam(OsiPrimalTolerance,value));
    break;
  default:
    abort();
  }
  return value;
}
int 
CbcOrClpParam::setIntParameter (OsiSolverInterface * model,int value) 
{
  if (value<lowerIntValue_||value>upperIntValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerIntValue_<<" to "<<
      upperIntValue_<<std::endl;
    return 1;
  } else {
    int oldValue;
    switch(type_) {
    case SOLVERLOGLEVEL:
      model->messageHandler()->setLogLevel(value);
      break;
    default:
      oldValue=0; // to avoid compiler message
      abort();
    }
    std::cout<<name_<<" was changed from "<<oldValue<<" to "
	     <<value<<std::endl;
    return 0;
  }
}
int 
CbcOrClpParam::intParameter (OsiSolverInterface * model) const
{
  int value=0;
  switch(type_) {
  case SOLVERLOGLEVEL:
    value=model->messageHandler()->logLevel();
    break;
  default:
    abort();
  }
  return value;
}
int
CbcOrClpParam::setDoubleParameter (CbcModel &model,double value) 
{
  if (value<lowerDoubleValue_||value>upperDoubleValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerDoubleValue_<<" to "<<
      upperDoubleValue_<<std::endl;
    return 1;
  } else {
    double oldValue;
    setDoubleValue(value);
    switch(type_) {
    case INFEASIBILITYWEIGHT:
      oldValue=model.getDblParam(CbcModel::CbcInfeasibilityWeight);
      model.setDblParam(CbcModel::CbcInfeasibilityWeight,value);
      break;
    case INTEGERTOLERANCE:
      oldValue=model.getDblParam(CbcModel::CbcIntegerTolerance);
      model.setDblParam(CbcModel::CbcIntegerTolerance,value);
      break;
    case INCREMENT:
      oldValue=model.getDblParam(CbcModel::CbcCutoffIncrement);
      model.setDblParam(CbcModel::CbcCutoffIncrement,value);
    case ALLOWABLEGAP:
      oldValue=model.getDblParam(CbcModel::CbcAllowableGap);
      model.setDblParam(CbcModel::CbcAllowableGap,value);
      break;
    case TIMELIMIT:
    { oldValue = model.getDblParam(CbcModel::CbcMaximumSeconds) ;
      model.setDblParam(CbcModel::CbcMaximumSeconds,value) ;
      break ; }
    case DUALTOLERANCE:
    case PRIMALTOLERANCE:
      setDoubleParameter(model.solver(),value);
      return 0; // to avoid message
    default:
      oldValue=0.0; // to avoid compiler message
      break;
    }
    std::cout<<name_<<" was changed from "<<oldValue<<" to "
	     <<value<<std::endl;
    return 0;
  }
}
double 
CbcOrClpParam::doubleParameter (CbcModel &model) const
{
  double value;
  switch(type_) {
  case INFEASIBILITYWEIGHT:
    value=model.getDblParam(CbcModel::CbcInfeasibilityWeight);
    break;
  case INTEGERTOLERANCE:
    value=model.getDblParam(CbcModel::CbcIntegerTolerance);
    break;
  case INCREMENT:
    value=model.getDblParam(CbcModel::CbcCutoffIncrement);
  case ALLOWABLEGAP:
    value=model.getDblParam(CbcModel::CbcAllowableGap);
    break;
  case TIMELIMIT:
  { value = model.getDblParam(CbcModel::CbcMaximumSeconds) ;
    break ; }
  case DUALTOLERANCE:
  case PRIMALTOLERANCE:
    value=doubleParameter(model.solver());
    break;
  default:
    abort();
  }
  return value;
}
int 
CbcOrClpParam::setIntParameter (CbcModel &model,int value) 
{
  if (value<lowerIntValue_||value>upperIntValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerIntValue_<<" to "<<
      upperIntValue_<<std::endl;
    return 1;
  } else {
    setIntValue(value);
    int oldValue;
    switch(type_) {
    case LOGLEVEL:
      oldValue = model.messageHandler()->logLevel();
      model.messageHandler()->setLogLevel(value);
      break;
    case SOLVERLOGLEVEL:
      oldValue = model.solver()->messageHandler()->logLevel();
      model.solver()->messageHandler()->setLogLevel(value);
      break;
    case MAXNODES:
      oldValue=model.getIntParam(CbcModel::CbcMaxNumNode);
      model.setIntParam(CbcModel::CbcMaxNumNode,value);
      break;
    case STRONGBRANCHING:
      oldValue=model.numberStrong();
      model.setNumberStrong(value);
      break;
    case NUMBERBEFORE:
      oldValue=model.numberBeforeTrust();
      model.setNumberBeforeTrust(value);
      break;
    default:
      oldValue=0; // to avoid compiler message
      break;
    }
    std::cout<<name_<<" was changed from "<<oldValue<<" to "
	     <<value<<std::endl;
    return 0;
  }
}
int 
CbcOrClpParam::intParameter (CbcModel &model) const
{
  int value;
  switch(type_) {
  case LOGLEVEL:
    value = model.messageHandler()->logLevel();
      break;
  case SOLVERLOGLEVEL:
    value = model.solver()->messageHandler()->logLevel();
      break;
  case MAXNODES:
    value = model.getIntParam(CbcModel::CbcMaxNumNode);
    break;
  case STRONGBRANCHING:
    value=model.numberStrong();
    break;
  case NUMBERBEFORE:
    value=model.numberBeforeTrust();
    break;
  default:
    abort();
  }
  return value;
}
#endif
void 
CbcOrClpParam::setIntValue ( int value )
{ 
  if (value<lowerIntValue_||value>upperIntValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerIntValue_<<" to "<<
      upperIntValue_<<std::endl;
  } else {
    intValue_=value;
  }
}
void 
CbcOrClpParam::setDoubleValue ( double value )
{ 
  if (value<lowerDoubleValue_||value>upperDoubleValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerDoubleValue_<<" to "<<
      upperDoubleValue_<<std::endl;
  } else {
    doubleValue_=value;
  }
}
void 
CbcOrClpParam::setStringValue ( std::string value )
{ 
  stringValue_=value;
}
static char line[1000];
static char * where=NULL;
extern int CbcOrClpRead_mode;
extern FILE * CbcOrClpReadCommand;
// Simple read stuff
std::string
CoinReadNextField()
{
  std::string field;
  if (!where) {
    // need new line
#ifdef COIN_USE_READLINE     
    if (CbcOrClpReadCommand==stdin) {
      // Get a line from the user. 
      where = readline ("Clp:");
      
      // If the line has any text in it, save it on the history.
      if (where) {
	if ( *where)
	  add_history (where);
	strcpy(line,where);
	free(where);
      }
    } else {
      where = fgets(line,1000,CbcOrClpReadCommand);
    }
#else
    if (CbcOrClpReadCommand==stdin) {
      fprintf(stdout,"Clp:");
      fflush(stdout);
    }
    where = fgets(line,1000,CbcOrClpReadCommand);
#endif
    if (!where)
      return field; // EOF
    where = line;
    // clean image
    char * lastNonBlank = line-1;
    while ( *where != '\0' ) {
      if ( *where != '\t' && *where < ' ' ) {
	break;
      } else if ( *where != '\t' && *where != ' ') {
	lastNonBlank = where;
      }
      where++;
    }
    where=line;
    *(lastNonBlank+1)='\0';
  }
  // munch white space
  while(*where==' '||*where=='\t')
    where++;
  char * saveWhere = where;
  while (*where!=' '&&*where!='\t'&&*where!='\0')
    where++;
  if (where!=saveWhere) {
    char save = *where;
    *where='\0';
    //convert to string
    field=saveWhere;
    *where=save;
  } else {
    where=NULL;
    field="EOL";
  }
  return field;
}

std::string
CoinReadGetCommand(int argc, const char *argv[])
{
  std::string field="EOL";
  while (field=="EOL") {
    if (CbcOrClpRead_mode>0) {
      if (CbcOrClpRead_mode<argc) {
	field = argv[CbcOrClpRead_mode++];
	if (field=="-") {
	  std::cout<<"Switching to line mode"<<std::endl;
	  CbcOrClpRead_mode=-1;
	  field=CoinReadNextField();
	} else if (field[0]!='-') {
	  if (CbcOrClpRead_mode!=2) {
	    std::cout<<"skipping non-command "<<field<<std::endl;
	    field="EOL"; // skip
	  } else {
	    // special dispensation - taken as -import name
	    CbcOrClpRead_mode--;
	    field="import";
	  }
	} else {
	  if (field!="--") {
	    // take off -
	    field = field.substr(1);
	  } else {
	    // special dispensation - taken as -import --
	    CbcOrClpRead_mode--;
	    field="import";
	  }
	}
      } else {
	field="";
      }
    } else {
      field=CoinReadNextField();
    }
  }
  //std::cout<<field<<std::endl;
  return field;
}
std::string
CoinReadGetString(int argc, const char *argv[])
{
  std::string field="EOL";
  if (CbcOrClpRead_mode>0) {
    if (CbcOrClpRead_mode<argc) {
      if (argv[CbcOrClpRead_mode][0]!='-') { 
	field = argv[CbcOrClpRead_mode++];
      } else if (!strcmp(argv[CbcOrClpRead_mode],"--")) {
	field = argv[CbcOrClpRead_mode++];
	// -- means import from stdin
	field = "-";
      }
    }
  } else {
    field=CoinReadNextField();
  }
  //std::cout<<field<<std::endl;
  return field;
}
// valid 0 - okay, 1 bad, 2 not there
int
CoinReadGetIntField(int argc, const char *argv[],int * valid)
{
  std::string field="EOL";
  if (CbcOrClpRead_mode>0) {
    if (CbcOrClpRead_mode<argc) {
      // may be negative value so do not check for -
      field = argv[CbcOrClpRead_mode++];
    }
  } else {
    field=CoinReadNextField();
  }
  int value=0;
  //std::cout<<field<<std::endl;
  if (field!="EOL") {
    // how do I check valid
    value =  atoi(field.c_str());
    *valid=0;
  } else {
    *valid=2;
  }
  return value;
}
double
CoinReadGetDoubleField(int argc, const char *argv[],int * valid)
{
  std::string field="EOL";
  if (CbcOrClpRead_mode>0) {
    if (CbcOrClpRead_mode<argc) {
      // may be negative value so do not check for -
      field = argv[CbcOrClpRead_mode++];
    }
  } else {
    field=CoinReadNextField();
  }
  double value=0.0;
  //std::cout<<field<<std::endl;
  if (field!="EOL") {
    // how do I check valid
    value = atof(field.c_str());
    *valid=0;
  } else {
    *valid=2;
  }
  return value;
}
/*
  Subroutine to establish the cbc parameter array. See the description of
  class CbcOrClpParam for details. Pulled from C..Main() for clarity. 
*/
void 
establishParams (int &numberParameters, CbcOrClpParam *const parameters)
{
  numberParameters=0;
  parameters[numberParameters++]=
    CbcOrClpParam("?","For help",GENERALQUERY,-1,false);
  parameters[numberParameters++]=
    CbcOrClpParam("???","For help",FULLGENERALQUERY,-1,false);
  parameters[numberParameters++]=
    CbcOrClpParam("-","From stdin",
		  STDIN,299,false);
#ifdef COIN_USE_CBC
    parameters[numberParameters++]=
      CbcOrClpParam("allow!ableGap","Stop when gap between best possible and \
best less than this",
	      0.0,1.0e20,ALLOWABLEGAP);
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("allS!lack","Set basis back to all slack and reset solution",
		  ALLSLACK,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Useful for playing around"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("auto!Scale","Whether to scale objective, rhs and bounds of problem if they look odd",
		  "off",AUTOSCALE,0,false);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    (
     "If you think you may get odd objective values or large equality rows etc then\
 it may be worth setting this true.  It is still experimental and you may prefer\
 to use objective!Scale and rhs!Scale"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("barr!ier","Solve using primal dual predictor corrector algorithm",
		  BARRIER);
  parameters[numberParameters-1].setLonghelp
    (
     "This command solves the current model using the  primal dual predictor \
corrector algorithm."
     
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("basisI!n","Import basis from bas file",
		  BASISIN);
  parameters[numberParameters-1].setLonghelp
    (
     "This will read an MPS format basis file from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  If you have libgz then it can read compressed\
 files 'xxxxxxxx.gz'.."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("basisO!ut","Export basis as bas file",
		  BASISOUT);
  parameters[numberParameters-1].setLonghelp
    (
     "This will write an MPS format basis file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.bas'."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("biasLU","Whether factorization biased towards U",
		  "UU",BIASLU,2,false);
  parameters[numberParameters-1].append("UX");
  parameters[numberParameters-1].append("LX");
  parameters[numberParameters-1].append("LL");
#endif
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("branch!AndCut","Do Branch and Cut",
		  BAB);
  parameters[numberParameters-1].setLonghelp
    (
     "This does branch and cut.  Far too many parameters apply to give a full description here.  \
the main thing is to think about which cuts to apply.  .. expand ..."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("bscale","Whether to scale in barrier",
		  "off",BARRIERSCALE,0,false);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters++]=
    CbcOrClpParam("chol!esky","Which cholesky algorithm",
		  "native",CHOLESKY,false);
  parameters[numberParameters-1].append("dense");
  //#ifdef FOREIGN_BARRIER
#ifdef WSSMP_BARRIER
  parameters[numberParameters-1].append("fudge!Long");
  parameters[numberParameters-1].append("wssmp");
#define REAL_BARRIER
#else
  parameters[numberParameters-1].append("fudge!Long_dummy");
  parameters[numberParameters-1].append("wssmp_dummy");
#endif
#ifdef UFL_BARRIER
  parameters[numberParameters-1].append("Uni!versityOfFlorida");
#define REAL_BARRIER
#else
  parameters[numberParameters-1].append("Uni!versityOfFlorida_dummy");
#endif
#ifdef TAUCS_BARRIER
  parameters[numberParameters-1].append("Taucs");
#define REAL_BARRIER
#else
  parameters[numberParameters-1].append("Taucs_dummy");
#endif
  //#endif
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("clique!Cuts","Whether to use Clique cuts",
		  "off",CLIQUECUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("cost!Strategy","How to use costs",
		  "off",COSTSTRATEGY);
  parameters[numberParameters-1].append("pri!orities");
  parameters[numberParameters-1].append("pseudo!costs(not implemented yet)");
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("crash","Whether to create basis for problem",
		  "off",CRASH);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("so!low_halim");
  parameters[numberParameters-1].append("ha!lim_solow(JJF mods)");
  parameters[numberParameters-1].append("4");
  parameters[numberParameters-1].append("5");
  parameters[numberParameters-1].setLonghelp
    (
     "If crash is set on and there is an all slack basis then Clp will flip or put structural\
 variables into basis with the aim of getting dual feasible.  On the whole dual seems to be\
 better without it and there alernative types of 'crash' for primal e.g. 'idiot' or 'sprint'. \
I have also added a variant due to Solow and Halim which is as on but just flip."); 
  parameters[numberParameters++]=
    CbcOrClpParam("cross!over","Whether to get a basic solution after barrier",
		  "on",CROSSOVER);
  parameters[numberParameters-1].append("off");
  parameters[numberParameters-1].append("maybe");
  parameters[numberParameters-1].setLonghelp
    (
     "Interior point algorithms do not obtain a basic solution (and \
the feasibility criterion is a bit suspect (JJF)).  This option will crossover \
to a basic solution suitable for ranging or branch and cut.  With the current state \
of quadratic it may be a good idea to switch off crossover for quadratic (and maybe \
presolve as well) - the option maybe does this."
     );
#endif
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("cutD!epth","Depth in tree at which to do cuts",
		  -1,999999,CUTDEPTH);
  parameters[numberParameters-1].setLonghelp
    (
     "Cut generators may be - off, on only at root, on if they look possible \
and on.  If they are done every node then that is that, but it may be worth doing them \
every so often.  The original method was every so many node but it may be more logical \
to do it whenever depth in tree is a multiple of K.  This option does that and defaults \
to 5."
     );
  parameters[numberParameters-1].setIntValue(5);
#endif 
  parameters[numberParameters++]=
    CbcOrClpParam("direction","Minimize or Maximize",
		  "min!imize",DIRECTION);
  parameters[numberParameters-1].append("max!imize");
  parameters[numberParameters-1].append("zero");
  parameters[numberParameters-1].setLonghelp
    (
     "The default is minimize - use 'direction maximize' for maximization.\n\
You can also use the parameters 'maximize' or 'minimize'."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("directory","Set Default directory for import etc.",
		  DIRECTORY,299);
  parameters[numberParameters-1].setLonghelp
    (
     "This sets the directory which import, export, saveModel and restoreModel will use.\
  It is initialized to './'"
     ); 
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("dualB!ound","Initially algorithm acts as if no \
gap between bounds exceeds this value",
		  1.0e-20,1.0e12,DUALBOUND);
  parameters[numberParameters-1].setLonghelp
    (
     "The dual algorithm in Clp is a single phase algorithm as opposed to a two phase\
 algorithm where you first get feasible then optimal.  If a problem has both upper and\
 lower bounds then it is trivial to get dual feasible by setting non basic variables\
 to correct bound.  If the gap between the upper and lower bounds of a variable is more\
 than the value of dualBound Clp introduces fake bounds so that it can make the problem\
 dual feasible.  This has the same effect as a composite objective function in the\
 primal algorithm.  Too high a value may mean more iterations, while too low a bound means\
 the code may go all the way and then have to increase the bounds.  OSL had a heuristic to\
 adjust bounds, maybe we need that here."
     );
  parameters[numberParameters++]=
    CbcOrClpParam("dualize","Solves dual reformulation",
		  0,1,DUALIZE,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Don't even think about it."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("dualP!ivot","Dual pivot choice algorithm",
		  "auto!matic",DUALPIVOT);
  parameters[numberParameters-1].append("dant!zig");
  parameters[numberParameters-1].append("partial");
  parameters[numberParameters-1].append("steep!est");
  parameters[numberParameters-1].setLonghelp
    (
     "Clp can use any pivot selection algorithm which the user codes as long as it\
 implements the features in the abstract pivot base class.  The Dantzig method is implemented\
 to show a simple method but its use is deprecated.  Steepest is the method of choice and there\
 are two variants which keep all weights updated but only scan a subset each iteration.\
 Partial switches this on while automatic decides at each iteration based on information\
 about the factorization."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("dualS!implex","Do dual simplex algorithm",
		  DUALSIMPLEX);
  parameters[numberParameters-1].setLonghelp
    (
     "This command solves the current model using the dual steepest algorithm.\
The time and iterations may be affected by settings such as presolve, scaling, crash\
 and also by dual pivot method, fake bound on variables and dual and primal tolerances."
     );
#endif 
  parameters[numberParameters++]=
    CbcOrClpParam("dualT!olerance","For an optimal solution \
no dual infeasibility may exceed this value",
		  1.0e-20,1.0e12,DUALTOLERANCE);
  parameters[numberParameters-1].setLonghelp
    (
     "Normally the default tolerance is fine, but you may want to increase it a\
 bit if a dual run seems to be having a hard time"
     ); 
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("either!Simplex","Do dual or primal simplex algorithm",
		  EITHERSIMPLEX);
  parameters[numberParameters-1].setLonghelp
    (
     "This command solves the current model using the dual or primal algorithm,\
 based on a dubious analysis of model."
     );
#endif 
  parameters[numberParameters++]=
    CbcOrClpParam("end","Stops clp execution",
		  EXIT);
  parameters[numberParameters-1].setLonghelp
    (
     "This stops execution ; end, exit, quit and stop are synonyms"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("error!sAllowed","Whether to allow import errors",
		  "off",ERRORSALLOWED);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    (
     "The default is not to use any model which had errors when reading the mps file.\
  Setting this to 'on' will allow all errors from which the code can recover\
 by ignoring the error."
     );
  parameters[numberParameters++]=
    CbcOrClpParam("exit","Stops clp execution",
		  EXIT);
  parameters[numberParameters-1].setLonghelp
    (
     "This stops the execution of Clp, end, exit, quit and stop are synonyms"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("export","Export model as mps file",
		  EXPORT);
  parameters[numberParameters-1].setLonghelp
    (
     "This will write an MPS format file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.mps'."
     ); 
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("fakeB!ound","All bounds <= this value - DEBUG",
		  1.0,1.0e15,FAKEBOUND,false);
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("fix!OnDj","Try heuristic based on fixing variables with \
reduced costs greater than this",
		  -1.0e20,1.0e20,DJFIX,false);
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
    parameters[numberParameters++]=
      CbcOrClpParam("flow!CoverCuts","Whether to use Flow Cover cuts",
		    "off",FLOWCUTS);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("gamma","Whether to regularize barrier",
		  "off",GAMMA,0,false);
  parameters[numberParameters-1].append("on");
#endif
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("gap!Ratio","Stop when gap between best possible and \
best less than this fraction of larger of two",
		  0.0,1.0e20,GAPRATIO);
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("gomory!Cuts","Whether to use Gomory cuts",
		  "off",GOMORYCUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].setLonghelp
    (
     "The original cuts - beware of imitations!  Having gone out of favo(u)r, they are now more \
fashionable as LP solvers are more robust and they interact well with other cuts.  They will almost always \
give cuts (although in this executable they are limited as to number of variables in cut)."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("help","Print out version, non-standard options and some help",
		  HELP);
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("idiot!Crash","Whether to try idiot crash",
		  -1,999999,IDIOT);
  parameters[numberParameters-1].setLonghelp
    (
     "This is a type of 'crash' which works well on some homogeneous problems.\
 It works best on problems with unit elements and rhs but will do something to any model.  It should only be\
 used before primal.  It can be set to -1 when the code decides for itself whether to use it,\
 0 to switch off or n > 0 to do n passes."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("import","Import model from mps file",
		  IMPORT);
  parameters[numberParameters-1].setLonghelp
    (
     "This will read an MPS format file from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  If you have libgz then it can read compressed\
 files 'xxxxxxxx.gz'.."
     );
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("inc!rement","A valid solution must be at least this \
much better than last integer solution",
		  -1.0e20,1.0e20,INCREMENT);
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("inf!easibilityWeight","Each integer infeasibility is expected \
to cost this much",
		  0.0,1.0e20,INFEASIBILITYWEIGHT);
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("initialS!olve","Solve to continuous",
		  SOLVECONTINUOUS);
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("integerT!olerance","For an optimal solution \
no integer variable may be this away from an integer value",
	      1.0e-20,0.5,INTEGERTOLERANCE);
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif 
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("keepN!ames","Whether to keep names from import",
		  "on",KEEPNAMES);
  parameters[numberParameters-1].append("off");
  parameters[numberParameters-1].setLonghelp
    (
     "It saves space to get rid of names so if you need to you can set this to off."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("KKT","Whether to use KKT factorization",
		  "off",KKT,0,false);
  parameters[numberParameters-1].append("on");
#endif
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("knapsack!Cuts","Whether to use Knapsack cuts",
		  "off",KNAPSACKCUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif
#ifndef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("log!Level","Level of detail in Solver output",
		  -1,63,SOLVERLOGLEVEL);
#else
  parameters[numberParameters++]=
    CbcOrClpParam("log!Level","Level of detail in Coin branch and Cut output",
		  -1,63,LOGLEVEL);
  parameters[numberParameters-1].setIntValue(1);
#endif
  parameters[numberParameters-1].setLonghelp
    (
     "If 0 then there should be no output in normal circumstances.  1 is probably the best\
 value for most uses, while 2 and 3 give more information."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("max!imize","Set optimization direction to maximize",
		  MAXIMIZE,299);
  parameters[numberParameters-1].setLonghelp
    (
     "The default is minimize - use 'maximize' for maximization.\n\
You can also use the parameters 'direction maximize'."
     ); 
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("maxF!actor","Maximum number of iterations between \
refactorizations",
		  1,999999,MAXFACTOR);
  parameters[numberParameters-1].setLonghelp
    (
     "If this is at its initial value of 201 then in this executable clp will guess at a\
 value to use.  Otherwise the user can set a value.  The code may decide to re-factorize\
 earlier."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("maxIt!erations","Maximum number of iterations before \
stopping",
		  0,99999999,MAXITERATION);
  parameters[numberParameters-1].setLonghelp
    (
     "This can be used for testing purposes.  The corresponding library call\n\
      \tsetMaximumIterations(value)\n can be useful.  If the code stops on\
 seconds or by an interrupt this will be treated as stopping on maximum iterations"
     ); 
#endif
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("maxN!odes","Maximum number of nodes to do",
		  1,999999,MAXNODES);
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("min!imize","Set optimization direction to minimize",
		  MINIMIZE,299);
  parameters[numberParameters-1].setLonghelp
    (
     "The default is minimize - use 'maximize' for maximization.\n\
This should only be necessary if you have previously set maximization \
You can also use the parameters 'direction minimize'."
     );
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("mixed!IntegerRoundingCuts","Whether to use Mixed Integer Rounding cuts",
		  "off",MIXEDCUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif 
  parameters[numberParameters++]=
    CbcOrClpParam("mess!ages","Controls if Clpnnnn is printed",
		  "off",MESSAGES);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    ("The default for the Clp library is to put out messages such as:\n\
   Clp0005 2261  Objective 109.024 Primal infeas 944413 (758)\n\
but this program turns this off to make it look more friendly.  It can be useful\
 to turn them back on if you want to be able 'grep' for particular messages or if\
 you intend to override the behavior of a particular message."
     );
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("miplib","Do some of miplib test set",
		  MIPLIB);
#endif 
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("netlib","Solve entire netlib test set",
		  NETLIB_EITHER);
  parameters[numberParameters-1].setLonghelp
    (
     "This exercises the unit test for clp and then solves the netlib test set using dual or primal.\
The user can set options before e.g. clp -presolve off -netlib"
     ); 
#ifdef REAL_BARRIER
  parameters[numberParameters++]=
    CbcOrClpParam("netlibB!arrier","Solve entire netlib test set with barrier",
		  NETLIB_BARRIER);
  parameters[numberParameters-1].setLonghelp
    (
     "This exercises the unit test for clp and then solves the netlib test set using barrier.\
The user can set options before e.g. clp -kkt on -netlib"
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("netlibD!ual","Solve entire netlib test set (dual)",
		  NETLIB_DUAL);
  parameters[numberParameters-1].setLonghelp
    (
     "This exercises the unit test for clp and then solves the netlib test set using dual.\
The user can set options before e.g. clp -presolve off -netlib"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("netlibP!rimal","Solve entire netlib test set (primal)",
		  NETLIB_PRIMAL);
  parameters[numberParameters-1].setLonghelp
    (
     "This exercises the unit test for clp and then solves the netlib test set using primal.\
The user can set options before e.g. clp -presolve off -netlibp"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("netlibT!une","Solve entire netlib test set with 'best' algorithm",
		  NETLIB_TUNE);
  parameters[numberParameters-1].setLonghelp
    (
     "This exercises the unit test for clp and then solves the netlib test set using whatever \
works best.  I know this is cheating but it also stresses the code better by doing a \
mixture of stuff.  The best algorithm was chosen on a Linux ThinkPad using native cholesky \
with University of Florida ordering."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("network","Tries to make network matrix",
		  NETWORK,-1,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Clp will go faster if the matrix can be converted to a network.  The matrix\
 operations may be a bit faster with more efficient storage, but the main advantage\
 comes from using a network factorization.  It will probably not be as fast as a \
specialized network code."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("objective!Scale","Scale factor to apply to objective",
		  -1.0e20,1.0e20,OBJSCALE,false);
  parameters[numberParameters-1].setLonghelp
    (
     "If the objective function has some very large values, you may wish to scale them\
 internally by this amount.  It can also be set by autoscale.  It is applied after scaling"
     ); 
  parameters[numberParameters-1].setDoubleValue(1.0);
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("output!Format","Which output format to use",
		  1,6,OUTPUTFORMAT);
  parameters[numberParameters-1].setLonghelp
    (
     "Normally export will be done using normal representation for numbers and two values\
 per line.  You may want to do just one per line (for grep or suchlike) and you may wish\
 to save with absolute accuracy using a coded version of the IEEE value. A value of 2 is normal.\
 otherwise odd values gives one value per line, even two.  Values 1,2 give normal format, 3,4\
 gives greater precision, while 5,6 give IEEE values.  When used for exporting a basis 1 does not save \
values, 2 saves values, 3 with greater accuracy and 4 in IEEE."
     ); 
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("passP!resolve","How many passes in presolve",
		  -200,100,PRESOLVEPASS,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Normally Presolve does 5 passes but you may want to do less to make it\
 more lightweight or do more if improvements are still being made.  As Presolve will return\
 if nothing is being taken out, then you should not need to use this fine tuning."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("pertV!alue","Method of perturbation",
		  -5000,102,PERTVALUE,false);
  parameters[numberParameters++]=
    CbcOrClpParam("perturb!ation","Whether to perturb problem",
		  "on",PERTURBATION);
  parameters[numberParameters-1].append("off");
  parameters[numberParameters-1].setLonghelp
    (
     "Perturbation helps to stop cycling, but Clp uses other measures for this.\
  However large problems and especially ones with unit elements and unit rhs or costs\
 benefit from perturbation.  Normally Clp tries to be intelligent, but you can switch this off.\
  The Clp library has this off by default.  This program has it on."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("PFI","Whether to use Product Form of Inverse in simplex",
		  "off",PFI,0,false);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    (
     "By default clp uses Forrest-Tomlin L-U update.  If you are masochistic you can switch it off."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("plus!Minus","Tries to make +- 1 matrix",
		  PLUSMINUS,-1,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Clp will go slightly faster if the matrix can be converted so that the elements are\
 not stored and are known to be unit.  The main advantage is memory use.  Clp may automatically\
 see if it can convert the problem so you should not need to use this."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("preO!pt","Presolve options",
		  0,INT_MAX,PRESOLVEOPTIONS,false);
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("presolve","Whether to presolve problem",
		  "on",PRESOLVE);
  parameters[numberParameters-1].append("off");
  parameters[numberParameters-1].append("more");
  parameters[numberParameters-1].append("file");
  parameters[numberParameters-1].setLonghelp
    (
     "Presolve analyzes the model to find such things as redundant equations, equations\
 which fix some variables, equations which can be transformed into bounds etc etc.  For the\
 initial solve of any problem this is worth doing unless you know that it will have no effect."
     ); 
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("preprocess","Whether to use integer preprocessing",
                  "on",PREPROCESS);
  parameters[numberParameters-1].append("off");
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("primalP!ivot","Primal pivot choice algorithm",
		  "auto!matic",PRIMALPIVOT);
  parameters[numberParameters-1].append("exa!ct");
  parameters[numberParameters-1].append("dant!zig");
  parameters[numberParameters-1].append("part!ial");
  parameters[numberParameters-1].append("steep!est");
  parameters[numberParameters-1].append("change");
  parameters[numberParameters-1].append("sprint");
  parameters[numberParameters-1].setLonghelp
    (
     "Clp can use any pivot selection algorithm which the user codes as long as it\
 implements the features in the abstract pivot base class.  The Dantzig method is implemented\
 to show a simple method but its use is deprecated.  Exact devex is the method of choice and there\
 are two variants which keep all weights updated but only scan a subset each iteration.\
 Partial switches this on while change initially does dantzig until the factorization\
 becomes denser.  This is still a work in progress."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("primalS!implex","Do primal simplex algorithm",
		  PRIMALSIMPLEX);
  parameters[numberParameters-1].setLonghelp
    (
     "This command solves the current model using the primal algorithm.\
  The default is to use exact devex.\
 The time and iterations may be affected by settings such as presolve, scaling, crash\
 and also by column selection  method, infeasibility weight and dual and primal tolerances."
     );
#endif 
  parameters[numberParameters++]=
    CbcOrClpParam("primalT!olerance","For an optimal solution \
no primal infeasibility may exceed this value",
		  1.0e-20,1.0e12,PRIMALTOLERANCE);
  parameters[numberParameters-1].setLonghelp
    (
     "Normally the default tolerance is fine, but you may want to increase it a\
 bit if a primal run seems to be having a hard time"
     ); 
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("primalW!eight","Initially algorithm acts as if it \
costs this much to be infeasible",
		  1.0e-20,1.0e20,PRIMALWEIGHT);
  parameters[numberParameters-1].setLonghelp
    (
     "The primal algorithm in Clp is a single phase algorithm as opposed to a two phase\
 algorithm where you first get feasible then optimal.  So Clp is minimizing this weight times\
 the sum of primal infeasibilities plus the true objective function (in minimization sense).\
  Too high a value may mean more iterations, while too low a bound means\
 the code may go all the way and then have to increase the weight in order to get feasible.\
  OSL had a heuristic to\
 adjust bounds, maybe we need that here."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("printO!ptions","Print options",
		  0,INT_MAX,PRINTOPTIONS,false);
#endif
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("probing!Cuts","Whether to use Probing cuts",
		  "off",PROBINGCUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("quit","Stops clp execution",
		  EXIT);
  parameters[numberParameters-1].setLonghelp
    (
     "This stops the execution of Clp, end, exit, quit and stop are synonyms"
     ); 
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("reallyS!cale","Scales model in place",
		  REALLY_SCALE,false);
#endif
#ifdef COIN_USE_CBC
    parameters[numberParameters++]=
      CbcOrClpParam("reduce!AndSplitCuts","Whether to use Reduce-and-Split cuts",
	      "off",REDSPLITCUTS);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("restore!Model","Restore model from binary file",
		  RESTORE);
  parameters[numberParameters-1].setLonghelp
    (
     "This reads data save by saveModel from the given file.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.prob'."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("reverse","Reverses sign of objective",
		  REVERSE,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Useful for testing if maximization works correctly"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("rhs!Scale","Scale factor to apply to rhs and bounds",
		  -1.0e20,1.0e20,RHSSCALE,false);
  parameters[numberParameters-1].setLonghelp
    (
     "If the rhs or bounds have some very large meaningful values, you may wish to scale them\
 internally by this amount.  It can also be set by autoscale"
     ); 
  parameters[numberParameters-1].setDoubleValue(1.0);
#endif
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("round!ingHeuristic","Whether to use Rounding heuristic",
		  "off",ROUNDING);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("save!Model","Save model to binary file",
		  SAVE);
  parameters[numberParameters-1].setLonghelp
    (
     "This will save the problem to the given file name for future use\
 by restoreModel.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.prob'."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("scal!ing","Whether to scale problem",
		  "off",SCALING);
  parameters[numberParameters-1].append("equi!librium");
  parameters[numberParameters-1].append("geo!metric");
  parameters[numberParameters-1].append("auto!matic");
  parameters[numberParameters-1].setLonghelp
    (
     "Scaling can help in solving problems which might otherwise fail because of lack of\
 accuracy.  It can also reduce the number of iterations.  It is not applied if the range\
 of elements is small.  When unscaled it is possible that there may be small primal and/or\
 infeasibilities."
     ); 
  parameters[numberParameters-1].setCurrentOption(3); // say auto
#ifndef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("sec!onds","maximum seconds",
		  -1.0,1.0e12,TIMELIMIT);
  parameters[numberParameters-1].setLonghelp
    (
     "After this many seconds clp will act as if maximum iterations had been reached.\
  In this program it is really only useful for testing but the library function\n\
      \tsetMaximumSeconds(value)\n can be useful."
     ); 
#endif
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("slp!Value","Number of slp passes before primal",
		  -1,50000,SLPVALUE,false);
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("solu!tion","Prints solution to file",
		  SOLUTION);
  parameters[numberParameters-1].setLonghelp
    (
     "This will write a primitive solution file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'stdout'."
     ); 
#ifdef COIN_USE_CBC
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("solv!e","Solve problem",
		  BAB);
  parameters[numberParameters-1].setLonghelp
    (
     "If there are no integer variables then this just solves LP.  If there are integer variables \
this does branch and cut.  Far too many parameters apply to give a full description here.  \
the main thing is to think about which cuts to apply.  .. expand ..."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("slog!Level","Level of detail in Solver output",
		  -1,63,SOLVERLOGLEVEL);
  parameters[numberParameters-1].setLonghelp
    (
     "If 0 then there should be no output in normal circumstances.  1 is probably the best\
 value for most uses, while 2 and 3 give more information."
     ); 
#endif
#endif
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("spars!eFactor","Whether factorization treated as sparse",
		  "on",SPARSEFACTOR,0,false);
  parameters[numberParameters-1].append("off");
  parameters[numberParameters++]=
    CbcOrClpParam("special!Options","Dubious options for Simplex - see ClpSimplex.hpp",
		  0,INT_MAX,SPECIALOPTIONS,false);
  parameters[numberParameters++]=
    CbcOrClpParam("sprint!Crash","Whether to try sprint crash",
		  -1,500,SPRINT);
  parameters[numberParameters-1].setLonghelp
    (
     "For long and thin problems this program may solve a series of small problems\
 created by taking a subset of the columns.  I introduced the idea as 'Sprint' after\
 an LP code of that name of the 60's which tried the same tactic (not totally successfully).\
  Cplex calls it 'sifting'.  -1 is automatic choice, 0 is off, n is number of passes"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("stat!istics","Print some statistics",
		  STATISTICS);
  parameters[numberParameters-1].setLonghelp
    (
     "This command prints crude statistics for the current model.\
 If log level >1 then more is printed.\
 These are for presolved model if presolve on (and unscaled)."
     );
#endif
  CbcOrClpParam("stdin","From stdin",
		STDIN,-1,false);
  parameters[numberParameters++]=
    CbcOrClpParam("stop","Stops clp execution",
		  EXIT);
  parameters[numberParameters-1].setLonghelp
    (
     "This stops the execution of Clp, end, exit, quit and stop are synonyms"
     ); 
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("strong!Branching","Number of variables to look at in strong branching",
		  0,999999,STRONGBRANCHING);
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("subs!titution","How long a column to substitute for in presolve",
		  0,10000,SUBSTITUTION,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Normally Presolve gets rid of 'free' variables when there are no more than 3 \
 variables in column.  If you increase this the number of rows may decrease but number of \
 elements may increase."
     ); 
#endif
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("tighten!Factor","Tighten bounds using this times largest \
activity at continuous solution",
		  1.0,1.0e20,TIGHTENFACTOR,false);
  parameters[numberParameters++] =
    CbcOrClpParam("time!Limit","Set a time limit for solving this problem",
		  -1.0,(double)(60*60*24*365*10),TIMELIMIT) ;
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif
#ifdef COIN_USE_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("tightLP","Poor person's preSolve for now",
		  TIGHTEN,-1,false);
#endif
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("trust!PseudoCosts","Number of branches before we trust pseudocosts",
		  0,999999,NUMBERBEFORE);
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif
#ifdef COIN_USE_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("two!MirCuts","Whether to use Two phase Mixed Integer Rounding cuts",
		  "off",TWOMIRCUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].setLonghelp
    (
     "TODO"
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("unitTest","Do unit test",
		  UNITTEST);
  parameters[numberParameters-1].setLonghelp
    (
     "This exercises the unit test for clp"
     ); 
  assert(numberParameters<CBCMAXPARAMETERS);
}
// Given a parameter type - returns its number in list
int whichParam (CbcOrClpParameterType name, 
		int numberParameters, CbcOrClpParam *const parameters)
{
  int i;
  for (i=0;i<numberParameters;i++) {
    if (parameters[i].type()==name)
      break;
  }
  assert (i<numberParameters);
  return i;
}
