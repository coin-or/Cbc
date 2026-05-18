// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinError.hpp"

bool COINUTILSLIB_EXPORT CoinError::printErrors_ = false;

/** A function to block the popup windows that windows creates when the code
    crashes */
#ifdef _MSC_VER
#include <windows.h>
#include <stdlib.h>
COINUTILSLIB_EXPORT
void WindowsErrorPopupBlocker()
{
  SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
  _set_abort_behavior(0, _WRITE_ABORT_MSG);
}
#else
COINUTILSLIB_EXPORT
void WindowsErrorPopupBlocker() {}
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
