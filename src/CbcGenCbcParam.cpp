/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/
/*
  This file is part of cbc-generic.
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include <string>
#include <cassert>

#include "CoinParam.hpp"

#include "CbcModel.hpp"

#include "CbcSolverSettings.hpp"
#include "CbcSolverParam.hpp"
#include "CbcModelParam.hpp"

namespace {


}

/*
  Constructors and destructors

  There's a generic constructor and one for integer, double, keyword, string,
  and action parameters.
*/

/*
  Default constructor.
*/
CbcModelParam::CbcModelParam()
  : CoinParam()
  , paramCode_(CbcModelParamCode(0))
  , obj_(0)
{
  /* Nothing to be done here */
}

/*
  Constructor for double parameter
*/
CbcModelParam::CbcModelParam(CbcModelParamCode code,
  std::string name, std::string help,
  double lower, double upper, double dflt,
  bool display)
  : CoinParam(name, help, lower, upper, dflt, display)
  , paramCode_(code)
  , obj_(0)
{
  /* Nothing to be done here */
}

/*
  Constructor for integer parameter
*/
CbcModelParam::CbcModelParam(CbcModelParamCode code,
  std::string name, std::string help,
  int lower, int upper, int dflt,
  bool display)
  : CoinParam(name, help, lower, upper, dflt, display)
  , paramCode_(code)
  , obj_(0)
{
  /* Nothing to be done here */
}

/*
  Constructor for keyword parameter.
*/
CbcModelParam::CbcModelParam(CbcModelParamCode code,
  std::string name, std::string help,
  std::string firstValue, int dflt,
  bool display)
  : CoinParam(name, help, firstValue, dflt, display)
  , paramCode_(code)
  , obj_(0)
{
  /* Nothing to be done here */
}

/*
  Constructor for string parameter.
*/
CbcModelParam::CbcModelParam(CbcModelParamCode code,
  std::string name, std::string help,
  std::string dflt,
  bool display)
  : CoinParam(name, help, dflt, display)
  , paramCode_(code)
  , obj_(0)
{
  /* Nothing to be done here */
}

/*
  Constructor for action parameter.
*/
CbcModelParam::CbcModelParam(CbcModelParamCode code,
  std::string name, std::string help,
  bool display)
  : CoinParam(name, help, display)
  , paramCode_(code)
  , obj_(0)
{
  /* Nothing to be done here */
}

/*
  Copy constructor.
*/
CbcModelParam::CbcModelParam(const CbcModelParam &orig)
  : CoinParam(orig)
  , paramCode_(orig.paramCode_)
  , obj_(orig.obj_)
{
  /* Nothing to be done here */
}

/*
  Clone
*/

CbcModelParam *CbcModelParam::clone()
{
  return (new CbcModelParam(*this));
}

CbcModelParam &CbcModelParam::operator=(const CbcModelParam &rhs)
{
  if (this != &rhs) {
    CoinParam::operator=(rhs);

    paramCode_ = rhs.paramCode_;
    obj_ = rhs.obj_;
  }

  return *this;
}

/*
  Destructor
*/
CbcModelParam::~CbcModelParam()
{ /* Nothing more to do */
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
