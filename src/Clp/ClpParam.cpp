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

#include "ClpParam.hpp"

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
ClpParam::ClpParam()
    : CoinParam(), paramCode_(ClpParamCode(0)), parameters_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Constructor for double parameter

  The default value is 0.
*/
ClpParam::ClpParam(int code, std::string name,
                   std::string help, double lower, double upper,
                   std::string longHelp, CoinDisplayPriority displayPriority)
   : CoinParam(name, help, lower, upper, longHelp, displayPriority),
     paramCode_(code), parameters_(0) {
   /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Constructor for integer parameter

  The default value is 0.
*/
ClpParam::ClpParam(int code, std::string name,
                   std::string help, int lower, int upper,
                   std::string longHelp, CoinDisplayPriority displayPriority)
    : CoinParam(name, help, lower, upper, longHelp, displayPriority),
      paramCode_(code), parameters_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Constructor for parameters with a string (or no) value (all others).
  Type is not optional to resolve ambiguity.

  The default value is "" for all such parameter types
*/
ClpParam::ClpParam(int code, std::string name,
                   CoinParam::CoinParamType type,
                   std::string help, std::string longHelp,
                   CoinDisplayPriority displayPriority)
   : CoinParam(name, type, help, longHelp, displayPriority),
      paramCode_(code), parameters_(0) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

/*
  Copy constructor.
*/
ClpParam::ClpParam(const ClpParam &orig)
    : CoinParam(orig), paramCode_(orig.paramCode_), parameters_(orig.parameters_) {
  /* Nothing to be done here */
}

//###########################################################################
//###########################################################################

std::string ClpParam::printString() const {
   std::ostringstream buffer;
   if (name_ == "directory") {
      buffer << "Current working directory is " << strValue_ << std::endl;
   } else if (name_.substr(0, 6) == "printM") {
      buffer << "Current value of printMask is " << strValue_ << std::endl;
   } else {
      buffer << "Current default (if $ as parameter) for " << name_ << " is "
             << strValue_ << std::endl;
   }
   return buffer.str();
}
//###########################################################################
//###########################################################################

/*
  Clone
*/

ClpParam *ClpParam::clone() { return (new ClpParam(*this)); }

ClpParam &ClpParam::operator=(const ClpParam &rhs) {
  if (this != &rhs) {
    CoinParam::operator=(rhs);

    paramCode_ = rhs.paramCode_;
    parameters_ = rhs.parameters_;
  }

  return *this;
}

//###########################################################################
//###########################################################################

/*
  Destructor
*/
ClpParam::~ClpParam() { /* Nothing more to do */
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
