/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file is part of cbc-generic.
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <string>
#include <cassert>

#include "CoinParam.hpp"

#include "CbcModel.hpp"

#include "CbcGenCtlBlk.hpp"
#include "CbcGenParam.hpp"

namespace {

char svnid[] = "$Id: CbcGenParam.cpp 1173 2009-06-04 09:44:10Z forrest $" ;

}

/*
  Constructors and destructors

  There's a generic constructor and one for integer, double, keyword, string,
  and action parameters.
*/

/*
  Default constructor.
*/
CbcGenParam::CbcGenParam ()
        : CoinParam(),
        paramCode_(CbcGenParamCode(0)),
        obj_(0)
{
    /* Nothing to be done here */
}


/*
  Constructor for double parameter
*/
CbcGenParam::CbcGenParam (CbcGenParamCode code,
                          std::string name, std::string help,
                          double lower, double upper, double dflt,
                          bool display)
        : CoinParam(name, help, lower, upper, dflt, display),
        paramCode_(code),
        obj_(0)
{
    /* Nothing to be done here */
}

/*
  Constructor for integer parameter
*/
CbcGenParam::CbcGenParam (CbcGenParamCode code,
                          std::string name, std::string help,
                          int lower, int upper, int dflt,
                          bool display)
        : CoinParam(name, help, lower, upper, dflt, display),
        paramCode_(code),
        obj_(0)
{
    /* Nothing to be done here */
}

/*
  Constructor for keyword parameter.
*/
CbcGenParam::CbcGenParam (CbcGenParamCode code,
                          std::string name, std::string help,
                          std::string firstValue, int dflt,
                          bool display)
        : CoinParam(name, help, firstValue, dflt, display),
        paramCode_(code),
        obj_(0)
{
    /* Nothing to be done here */
}

/*
  Constructor for string parameter.
*/
CbcGenParam::CbcGenParam (CbcGenParamCode code,
                          std::string name, std::string help,
                          std::string dflt,
                          bool display)
        : CoinParam(name, help, dflt, display),
        paramCode_(code),
        obj_(0)
{
    /* Nothing to be done here */
}

/*
  Constructor for action parameter.
*/
CbcGenParam::CbcGenParam (CbcGenParamCode code,
                          std::string name, std::string help,
                          bool display)
        : CoinParam(name, help, display),
        paramCode_(code),
        obj_(0)
{
    /* Nothing to be done here */
}


/*
  Copy constructor.
*/
CbcGenParam::CbcGenParam (const CbcGenParam &orig)
        : CoinParam(orig),
        paramCode_(orig.paramCode_),
        obj_(orig.obj_)
{
    /* Nothing to be done here */
}

/*
  Clone
*/

CbcGenParam *CbcGenParam::clone ()
{
    return (new CbcGenParam(*this)) ;
}

CbcGenParam &CbcGenParam::operator= (const CbcGenParam & rhs)
{
    if (this != &rhs) {
        CoinParam::operator=(rhs) ;

        paramCode_ = rhs.paramCode_ ;
        obj_ = rhs.obj_ ;
    }

    return *this ;
}

/*
  Destructor
*/
CbcGenParam::~CbcGenParam ()
{ /* Nothing more to do */ }

