/* $Id$ */
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Shamelessly adapted from ClpEventHandler.

#include "CoinPragma.hpp"

#include "CbcEventHandler.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------

CbcEventHandler::CbcEventHandler(CbcModel *model)
  : model_(model)
  , dfltAction_(CbcEventHandler::noAction)
  , eaMap_(0)
{ /* nothing more required */
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
/*
  Here we need to clone the event/action map, if it exists
*/
CbcEventHandler::CbcEventHandler(const CbcEventHandler &rhs)
  : model_(rhs.model_)
  , dfltAction_(rhs.dfltAction_)
  , eaMap_(0)
{
  if (rhs.eaMap_ != 0) {
    eaMap_ = new eaMapPair(*rhs.eaMap_);
  }
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CbcEventHandler &
CbcEventHandler::operator=(const CbcEventHandler &rhs)
{
  if (this != &rhs) {
    model_ = rhs.model_;
    dfltAction_ = rhs.dfltAction_;
    if (rhs.eaMap_ != 0) {
      eaMap_ = new eaMapPair(*rhs.eaMap_);
    } else {
      eaMap_ = 0;
    }
  }
  return (*this);
}

//----------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CbcEventHandler *
CbcEventHandler::clone() const
{
  return (new CbcEventHandler(*this));
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
/*
  Take care to free the event/action map.
*/
CbcEventHandler::~CbcEventHandler()
{
  if (eaMap_ != 0)
    delete eaMap_;
}

//-------------------------------------------------------------------
// event() -- return the action for an event.
//-------------------------------------------------------------------

CbcEventHandler::CbcAction CbcEventHandler::event(CbcEvent event)
/*
  If an event/action map exists and contains an entry for the event, return it.
  Otherwise return the default action.
*/
{
  if (eaMap_ != 0) {
    eaMapPair::iterator entry = eaMap_->find(event);
    if (entry != eaMap_->end()) {
      return (entry->second);
    } else {
      return (dfltAction_);
    }
  } else {
    return (dfltAction_);
  }
}

//-------------------------------------------------------------------
// event() -- return the action for an event.
//-------------------------------------------------------------------

CbcEventHandler::CbcAction CbcEventHandler::event(CbcEvent event,
  void * /*data*/)
/*
  If an event/action map exists and contains an entry for the event, return it.
  Otherwise return the default action.
*/
{
  if (eaMap_ != 0) {
    eaMapPair::iterator entry = eaMap_->find(event);
    if (entry != eaMap_->end()) {
      return (entry->second);
    } else {
      return (dfltAction_);
    }
  } else {
    return (dfltAction_);
  }
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
