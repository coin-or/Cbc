/*
  Copyright (C) 2007
  Lou Hafer, International Business Machines Corporation and others. All
  Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).
*/
/*
  This file is part of cbc-generic.
*/

#ifndef CbcGenMessages_H
#define CbcGenMessages_H


/*! \file

  This file contains the enum that defines symbolic names for for cbc-generic
  messages.
*/

/*
  $Id: CbcGenMessages.hpp 1173 2009-06-04 09:44:10Z forrest $
*/

/*
  There's arguably not enough content here to justify a separate file, but it
  maintains the common pattern for COIN code.
*/

/*! \brief Symbolic names for cbc-generic messages

  These are the `internal IDs' for cbc-generic messages.
*/

typedef enum {
    CBCGEN_TEST_MSG = 1,
    CBCGEN_NEW_SOLVER,
    CBCGEN_CONFUSION,
    CBCGEN_DUMMY_END
} CbcGenMsgCode ;

#endif

