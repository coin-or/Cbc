/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include <cassert>
#include <string>
#include <sstream>

#include "CoinUtilsConfig.h"

#include "CoinFileIO.hpp"
#include "CoinFinite.hpp"

#include "CbcModel.hpp"

#include "CbcParam.hpp"
#include "CbcParameters.hpp"

/*
  Constructors and destructors

  There's a generic constructor and one for integer, double, keyword, string,
  and action parameters.
*/

//###########################################################################
//###########################################################################

/*
  Default constructor.
*/
CbcParam::CbcParam()
   : CoinParam(), paramCode_(CbcParamCode(0)), parameters_(0), model_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Constructor for double parameter
*/
CbcParam::CbcParam(int code, std::string name,
                   std::string help, double lower, double upper,
                   double defaultValue, std::string longHelp,
                   CoinDisplayPriority displayPriority)
   : CoinParam(name, help, lower, upper, defaultValue, longHelp,
               displayPriority),
     paramCode_(code), parameters_(0), model_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Constructor for integer parameter
*/
CbcParam::CbcParam(int code, std::string name,
                   std::string help, int lower, int upper,
                   int defaultValue, std::string longHelp,
                   CoinDisplayPriority displayPriority)
   : CoinParam(name, help, lower, upper, defaultValue, longHelp,
                displayPriority),
     paramCode_(code), parameters_(0), model_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Constructor for keyword parameter.
*/
CbcParam::CbcParam(int code, std::string name,
                   std::string help, std::string defaultKwd,
                   int defaultMode, std::string longHelp,
                   CoinDisplayPriority displayPriority)
    : CoinParam(name, help, defaultKwd, defaultMode, longHelp, displayPriority),
      paramCode_(code), parameters_(0), model_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Constructor for string parameter.
*/
CbcParam::CbcParam(int code, std::string name,
                   std::string help, std::string defaultValue,
                   std::string longHelp,
                   CoinDisplayPriority displayPriority)
    : CoinParam(name, help, defaultValue, longHelp, displayPriority),
      paramCode_(code), parameters_(0), model_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Constructor for action parameter.
*/
CbcParam::CbcParam(int code, std::string name,
                   std::string help, std::string longHelp,
                   CoinDisplayPriority displayPriority)
    : CoinParam(name, help, longHelp, displayPriority), paramCode_(code),
      parameters_(0), model_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Copy constructor.
*/
CbcParam::CbcParam(const CbcParam &orig)
   : CoinParam(orig), paramCode_(orig.paramCode_), parameters_(orig.parameters_),
     model_(orig.model_) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Clone
*/

CbcParam *CbcParam::clone() { return (new CbcParam(*this)); }

CbcParam &CbcParam::operator=(const CbcParam &rhs) {
  if (this != &rhs) {
    CoinParam::operator=(rhs);

    paramCode_ = rhs.paramCode_;
    model_ = rhs.model_;
    parameters_ = rhs.parameters_;
  }

  return *this;
}

//###########################################################################
//###########################################################################

/*
  Destructor
*/
CbcParam::~CbcParam() { /* Nothing more to do */
}

/*
  Utility routine to save the current solution to a file. No formatting, and
  not intended to be portable in any way, shape, or form.
*/


/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
