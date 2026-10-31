// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file CbcMakeTest.cpp
    \brief Main routine for running quick tests that need no tools beyond the build environment, similar to 'make test'.
*/

// The cbc 'make test' builds multiple executables and runs them.
// For MSVS builds, the CbcMakeTest project is a single executable that calls each of the unit tests in turn. 
// This file is the main routine for that executable.
// As in cbc make QUICK_TESTS, the following tests are run in order, 
// and the first failure aborts the test suite with a non-zero exit code:
// 
// lazy-constraint-test
// mipstart-test
// coefstr-test
// row-activity-test

#pragma once

#include "..\src\CbcConfig.h"

#include <cstdio>

#include "CoinError.hpp"

// The following #defines rename the global variables and functions in each test file 
// so that they do not conflict with each other when linked into a single executable.

#define g_failures g_failureslazy_constraint_test
#define check checklazy_constraint_test
#define checkClose checkCloselazy_constraint_test
#define checkEqInt checkEqIntlazy_constraint_test
#define main mainlazy_constraint_test
#include "..\test\lazy-constraint-test.cpp"

#define g_failures g_failuresmipstart_test
#define check checkmipstart_test
#define checkClose checkClosemipstart_test
#define checkEqInt checkEqIntmipstart_test
#define main mainmipstart_test
#include "..\test\mipstart-test.cpp"

#define g_failures g_failurescoefstr_test
#define check checkcoefstr_test
#define checkClose checkClosecoefstr_test
#define checkEqInt checkEqIntcoefstr_test
#define main maincoefstr_test
#include "..\test\coefstr-test.cpp"

#define g_failures g_failuresrow_activity_test
#define check checkrow_activity_test
#define checkClose checkCloserow_activity_test
#define checkEqInt checkEqIntrow_activity_testt
#define main mainrow_activity_test
#include "..\test\row-activity-test.cpp"


#undef g_failures
#undef check
#undef checkClose
#undef checkEqInt
#undef main

// Run the 'test' suite in order, aborting on the first failure.
int test(int argc, char **argv)
{
  int returnCode = 0;

  //"Exit code: 0 = all checks passed, 1 = one or more failed.\n");
  printf("Running lazy-contraint-test");
  returnCode = ::mainlazy_constraint_test(argc, argv);
  if (returnCode != 0)
    return returnCode;

  printf("Running mipstart-test");
  returnCode = ::mainmipstart_test(argc, argv);
  if (returnCode != 0)
    return returnCode;

  printf("Running coefstr-test");
  returnCode = ::maincoefstr_test(argc, argv);
  if (returnCode != 0)
    return returnCode;

  printf("Running row-activity-test");
  returnCode = ::mainrow_activity_test(argc, argv);
  if (returnCode != 0)
    return returnCode;

  return returnCode;
}

int main(int argc, char **argv)
{
  // With MSVC, make a crash or failed assertion end the program instead of
  // opening a dialog that waits for a click (and hangs scripts and CI).
  WindowsErrorPopupBlocker();

  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "--help" || arg == "-h") {
      printf(
        "CbcMakeTest - Runs tests for Cbc aligned to 'make test'.\n"
        "\n"
        "Usage: CbcMakeTest [--help]\n"
        "\n"
        "Runs the 'test' suite of tests.\n"
        "The 'sanity-test' suite is not supported.\n"
        "\n"
        "Exit code: 0 = all checks passed, 1 = one or more checks failed.\n");
      return 0;
    }
    fprintf(stderr, "CbcMakeTest: unrecognized argument '%s' "
                    "(try --help)\n",
      arg.c_str());
    return 1;
  }

  return test(argc, argv);
}


