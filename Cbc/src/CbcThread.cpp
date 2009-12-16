/* $Id: CbcModel.cpp 1261 2009-10-30 12:45:20Z forrest $ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CbcConfig.h"

#include <string>
#include <cassert>
#include <cmath>
#include <cfloat>

#include "CbcEventHandler.hpp"

#include "OsiSolverInterface.hpp"
#include "CbcThread.hpp"
#include "CbcModel.hpp"
#include "CbcFathom.hpp"

#include "CoinTime.hpp"
CbcThread::CbcThread()
{
}
// Constructor with base model
CbcThread::CbcThread (CbcModel & model, int deterministic, CbcModel * baseModel)
{
}
// Destructor
CbcThread::~CbcThread()
{
}
