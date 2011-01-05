/*! \legal
  Copyright (C) 2007
  Lou Hafer, International Business Machines Corporation and others. All
  Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file is part of cbc-generic.
*/

#include "CbcGenMessages.hpp"

#include "CbcGenCtlBlk.hpp"

namespace {

char svnid[] = "$Id: CbcGenMessages.cpp 1173 2009-06-04 09:44:10Z forrest $" ;

}

/*
  Begin file local namespace
*/
namespace {

/*
  Message definitions.

  The precise form isn't important here, so long as the method that loads them
  into a CoinMessages object can come up with values for external ID,
  detail level, and format string. The use of an enum to provide an internal ID
  for each message is mainly useful with the internationalisation feature. It
  makes it easy to slap the same ID on alternate versions of a message.
*/


typedef struct {
    CbcGenMsgCode inID ;
    int exID ;
    int lvl ;
    const char *fmt ;
} MsgDefn ;

static MsgDefn us_en_defns[] = {
    // informational (0 -- 2999)
    { CBCGEN_TEST_MSG, 1, 2, "This is the us_en test message, eh." },
    { CBCGEN_NEW_SOLVER, 2, 2, "Solver is now \"%s\"." },
    // warning (3000 -- 5999)
    // Non-fatal errors (6000 -- 8999)
    // Fatal errors (9000 and up)
    { CBCGEN_CONFUSION, 9001, 1, "Internal confusion, line %d." },
    { CBCGEN_DUMMY_END, 999999, 0, "" }
} ;

/*
  We seem to need a dummy CoinMessages object to prevent the compiler from
  complaining that CoinMessages::Language is unintialised.

  const CoinMessages dummy(0) ;
*/
/*
  The author is Canadian, eh. But we'll go with us_en anyways.
*/
const CoinMessages::Language default_language = CoinMessages::us_en ;

} /* End file local namespace */

/*!
  This function constructs a CoinMessages object filled with a default set of
  messages, overlaid with whatever is available for the specified language.
  It is used to establish the initial set of messages, and is also called
  whenever the language is changed. The latter, because there's no way of
  guaranteeing that the message sets for alternate languages will all
  replace the same messages. This approach guarantees that the set of
  messages is always composed of the default language overlaid with any
  messages for an alternate language.

  The default for lang is us_en, specified up in CbcGenCtlBlk.hpp. If you want
  to change the default language, change the declaration there. That said,
  you'll also have to provide the necessary message definitions and augment the
  case statements below.
*/

void CbcGenCtlBlk::setMessages (CoinMessages::Language lang)

{
    /*
      If messages exist, in the correct language, we have nothing more to do.
    */
    if (msgs_ && cur_lang_ == lang) {
        return ;
    }
    /*
      Otherwise, we need to do a wholesale rebuild. Create a new object of the
      appropriate size.
    */
    CoinMessages *msgs = new CoinMessages(sizeof(us_en_defns) / sizeof(MsgDefn)) ;

    msgs->setLanguage(lang) ;
    strcpy(msgs->source_, "CbcG");
    /*
      Yes, this is gloriously redundant, but it's set up in anticipation of
      future extensions.
    */
    MsgDefn *msgdefn ;
    switch (lang) {
    case CoinMessages::us_en: {
        msgdefn = us_en_defns ;
        break ;
    }
    default: {
        msgdefn = us_en_defns ;
        break ;
    }
    }
    /*
      Open a loop to create and load the messages.
    */
    while (msgdefn->inID != CBCGEN_DUMMY_END) {
        CoinOneMessage msg(msgdefn->exID, msgdefn->lvl, msgdefn->fmt) ;
        msgs->addMessage(msgdefn->inID, msg) ;
        msgdefn++ ;
    }
    /*
      Now, if the local language differs from the default language, load any
      overrides. Again, useless now, but maybe in the future ...
    */
    if (lang != cur_lang_) {
        switch (lang) {
        case CoinMessages::us_en: {
            msgdefn = us_en_defns ;
            break;
        }
        default: {
            msgdefn = us_en_defns ;
            break;
        }
        }

        while (msgdefn->inID != CBCGEN_DUMMY_END) {
            msgs->replaceMessage(msgdefn->inID, msgdefn->fmt) ;
            msgdefn++ ;
        }
    }
    /*
      Each CoinOneMessage has a fixed-length array to hold the message; by default
      this is 400 chars. Convert to `compressed' CoinOneMessage objects where the
      array is only as large as necessary. Any attempt to replace a message, or the
      message text, will automatically trigger a decompress operation before doing
      the replacement, but the messages will *not* be automatically recompressed.
    */
    msgs->toCompact() ;
    msgs_ = msgs ;

    return ;
}


/*
  Replaces the current message handler with the handler supplied as a
  parameter. If ourMsgHandler_ is true, the existing handler is destroyed.
*/

void CbcGenCtlBlk::passInMessageHandler (CoinMessageHandler *newMsgHandler)

{
    if (msgHandler_ && ourMsgHandler_) {
        delete msgHandler_ ;
    }

    msgHandler_ = newMsgHandler ;
    ourMsgHandler_ = false ;

    return ;
}

/*
  Start a message. This routine buries the whole business of locating the
  message handler and messages, getting the log level right, etc. If, by some
  chance, messages are not yet loaded, do so.
*/

CoinMessageHandler &CbcGenCtlBlk::message (CbcGenMsgCode inID)

{
    if (!msgs_) {
        setMessages() ;
    }
    msgHandler_->setLogLevel(logLvl_) ;

    msgHandler_->message(inID, *msgs_) ;

    return (*msgHandler_) ;
}

