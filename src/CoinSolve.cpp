// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file CbcSolver.cpp
    \brief Main routine for the cbc stand-alone solver.
*/

#include "CbcConfig.h"

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <iostream>

#include "CoinPragma.hpp"
#include "CoinModel.hpp"
#include "CbcModel.hpp"
#include "CbcSolver.hpp"
#include "CbcDebug.hpp"
#include "OsiClpSolverInterface.hpp"

/*
  We have the following compile-time symbols.

  CBC_OTHER_SOLVER	CoinSolve.cpp, CbcSolver.[cpp,hpp], CbcModel.cpp

    A value of 1 says `cplex'. Other values not defined. The intent with
    cplex is to apply all of cbc's smarts at the root, then hand the problem
    over to cplex to finish. Cplex is not used as an alternate lp solver
    under cbc control.

    Usage in CbcModel is a fake; a small bit of code that's now `#if 0'.


  CLP_DEBUG_MALLOC

    This ties in with the functions clp_malloc, clp_free, and clp_memory,
    which are defined in CoinOslFactorization.cpp. (Right where you'd expect
    to find them, eh?).  Looks to be a relatively nice debugging wrapper for
    standard C malloc.  Calls standard C malloc/free directly if
    CLP_DEBUG_MALLOC is not defined.  Worth consideration for breaking out as
    a separate utility.  The hooks for new and delete defined here should be
    incorporated.

    Absolutely not thread safe --- lots of static variables.

    Hmmm ... is it still the case that standard C malloc and C++ new/delete
    do not play well together? 'Cause the hooks here for new and delete will
    not escape from this file.
*/

#if CBC_OTHER_SOLVER == 1
#ifndef CBC_HAS_OSICPX
#error "Configuration did not detect OsiCpx installation."
#else
#include "OsiCpxSolverInterface.hpp"
#endif
#endif

/*
  Hooks for a debugging wrapper for malloc/free. This bit of definition hooks
  C++ new / delete and diverts them into the debugging wrapper.
*/

//#define CLP_DEBUG_MALLOC
#ifdef CLP_DEBUG_MALLOC
void *operator new(size_t size) throw(std::bad_alloc)
{
  void *p = clp_malloc(size);
  return p;
}
void operator delete(void *p) throw()
{
  clp_free(p);
}
#endif // CLP_DEBUG_MALLOC

//#############################################################################
//#############################################################################

// To use USERCBC or USERCLP change 0 to 1 in defines and add in your fake
// main program(s) and any other code
//#define USER_HAS_FAKE_CBC
//#define USER_HAS_FAKE_CLP

#ifdef USER_HAS_FAKE_CBC
#endif
void fakeMain(ClpSimplex &model, OsiSolverInterface & /*osiSolver*/, CbcModel &babSolver)
{
#ifdef USER_HAS_FAKE_CBC
#else
  printf("Dummy user cbc code - model has %d rows and %d columns\n",
    model.getNumRows(), model.getNumCols());
  // Reduce printout
  babSolver.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  // Do complete search
  babSolver.branchAndBound();
#endif
}

// Clp stuff
#ifdef USER_HAS_FAKE_CLP
#endif
void fakeMain2(ClpSimplex & /*model*/,
  OsiClpSolverInterface &osiSolver,
  int /*options*/)
{
#ifdef USER_HAS_FAKE_CLP
#else
  ClpSimplex *lpSolver = osiSolver.getModelPtr();
  printf("Dummy user clp code - model has %d rows and %d columns\n",
    lpSolver->numberRows(), lpSolver->numberColumns());
  osiSolver.initialSolve();
#endif
}
//  End any fake main program

//#############################################################################
//#############################################################################

// void CbcClpUnitTest (const CbcModel & saveModel);

#ifdef CBC_STATISTICS
int osi_crunch = 0;
static int cbc_resolve = 0;
int osi_primal = 0;
int osi_dual = 0;
int osi_hot = 0;
void cbc_resolve_check(const OsiSolverInterface *solver)
{
  cbc_resolve++;
  printf("R %d stats %d %d %d\n",
    cbc_resolve, solver->getNumRows(), solver->getNumCols(),
    solver->getMatrixByCol()->getNumElements());
  if ((cbc_resolve % 1000) == 0)
    printf("RR %d resolve crunch %d primal %d dual %d hot %d\n",
      cbc_resolve, osi_crunch, osi_primal, osi_dual, osi_hot);
}
#endif

//#############################################################################
//#############################################################################

/*
  Somehow with some BLAS we get multithreaded by default
  For 99.99% of problems this is not a good idea.
  The openblas_set_num_threads(1) seems to work even with other blas
 */
#if CLP_USE_OPENBLAS
extern "C" {
void openblas_set_num_threads(int num_threads);
}
#endif

//#############################################################################
//#############################################################################

#ifdef TEST_MESSAGE_HANDLER
void test_message_handler(CbcModel &model){
    MyMessageHandler2 messageHandler(&model);
    std::cout << "Testing derived message handler" << std::endl;
    model.passInMessageHandler(&messageHandler);
    OsiClpSolverInterface *clpSolver =
       dynamic_cast< OsiClpSolverInterface * >(model.solver());
    // Could use different handlers (if different log levels)
    clpSolver->passInMessageHandler(&messageHandler);
    //clpSolver->getModelPtr()->passInMessageHandler(&messageHandler);

    // Set log levels same so can use one message handler
    clpSolver->messageHandler()->setLogLevel(1);
    model.messageHandler()->setLogLevel(1);
}
#endif

//#############################################################################
//#############################################################################

void formInputQueue(std::deque<std::string> &inputQueue,
                    int argc, char **argv)
{
   for (int i = 1; i < argc; i++){
      std::string tmp(argv[i]);
      std::string::size_type found = tmp.find('=');
      if (found != std::string::npos) {
         inputQueue.push_back(tmp.substr(0, found));
         inputQueue.push_back(tmp.substr(found + 1));
      } else {
         inputQueue.push_back(tmp);
      }
   }
}

//#############################################################################
//#############################################################################

int main(int argc, const char *argv[])
{
  int returnCode = 0;
#ifdef CLP_DEBUG_MALLOC
  clp_memory(0);
#endif
#ifndef CBC_OTHER_SOLVER
  OsiClpSolverInterface solver1;
#if CLP_USE_OPENBLAS
  openblas_set_num_threads(CLP_USE_OPENBLAS);
#endif
#elif CBC_OTHER_SOLVER == 1
  OsiCpxSolverInterface solver1;
#endif
  CbcModel model(solver1);
  
  CbcParameters parameters;
#ifndef CBC_NO_INTERRUPT
  parameters.enableSignalHandler();
#endif
  // initialize
  CbcMain0(model, parameters);
  
  // define TEST_MESSAGE_HANDLER at top of file to check works on all messages
#ifdef TEST_MESSAGE_HANDLER
  test_messaage_handler(&model);
#endif
     
  std::deque<std::string> inputQueue;
  
  if (argc > 2 && !strcmp(argv[2], "-AMPL")) {
     ampl_info info;
     returnCode = cbcReadAmpl(&info, argc, const_cast< char ** >(argv), model);
     if (!returnCode) {
        // Put arguments into a queue.
        // This should be moved to constructor of ClpSolver
        formInputQueue(inputQueue, info.numberArguments, info.arguments);
        // We don't need to first two arguments from here on
        inputQueue.pop_front();
        inputQueue.pop_front();
        returnCode = CbcMain1(inputQueue, model, parameters, dummyCallback,
                              &info);
     }
  } else {
     // Put arguments into a queue.
     formInputQueue(inputQueue, argc, const_cast< char ** >(argv));
     returnCode = CbcMain1(inputQueue, model, parameters);
  }     

#ifdef CLP_DEBUG_MALLOC
  clp_memory(1);
#endif

  if (returnCode != 777) {
    return returnCode;
  } else {
    return 0;
  }
}

//#############################################################################
// Some old comments for historical interest
//#############################################################################
/*
  Version 1.00.00 November 16 2005.
  This is to stop me (JJF) messing about too much.
  Tuning changes should be noted here.
  The testing next version may be activated by CBC_NEXT_VERSION
  This applies to OsiClp, Clp etc
  Version 1.00.01 November 24 2005
  Added several classes for advanced users.  This can't affect code (if you don't use it)
  Made some tiny changes (for N way branching) which should not change anything.
  CbcNWay object class - for N way branching this also allows use of CbcConsequence class.
  CbcBranchAllDifferent object class - for branching on general integer variables
  to stop them having same value so branches are x >= y+1 and x <= y-1.
  Added two new Cgl classes - CglAllDifferent which does column fixing (too slowly)
  and CglStored which just has a list of cuts which can be activated.
  Modified preprocess option to SOS
  Version 1.00.02 December 9 2005
  Added use of CbcStrategy to do clean preprocessing
  Added use of referenceSolver for cleaner repetition of Cbc
  Version 1.01.00 February 2 2006
  Added first try at Ampl interface
  Made dummy program so real main can be called from other programs
*/

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
