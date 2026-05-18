// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef ClpSolver_H
#define ClpSolver_H

#include "ClpConfig.h"

#include <string>
#include <vector>

#include "ClpSimplex.hpp"
#include "ClpMessage.hpp"

// for printing
#ifndef CLP_OUTPUT_FORMAT
#define CLP_OUTPUT_FORMAT % 15.8g
#endif
#define CLP_QUOTE(s) CLP_STRING(s)
#define CLP_STRING(s) #s

#ifdef CLP_USEFUL_PRINTOUT
extern double debugDouble[10];
extern int debugInt[24];
#endif
#if defined(CLP_HAS_AMD)
#define FOREIGN_BARRIER
#endif

static bool maskMatches(const int *starts, char **masks, std::string &check);
static ClpSimplex *currentModel = NULL;

extern "C" {
   static void
#if defined(_MSC_VER)
   __cdecl
#endif // _MSC_VER
   signal_handler(int /*whichSignal*/)
   {
      if (currentModel != NULL)
         currentModel->setMaximumIterations(0); // stop at next iterations
      return;
   }
   /*
     Somehow with some BLAS we get multithreaded by default
     For 99.99% of problems this is not a good idea.
     The openblas_set_num_threads(1) seems to work even with other blas
   */
   void openblas_set_num_threads(int num_threads);
}

int mainTest(int argc, const char *argv[], int algorithm,
  ClpSimplex empty, ClpSolve solveOptions, int switchOff, bool doVector);

CLPLIB_EXPORT
void ClpMain0(ClpSimplex &model);
CLPLIB_EXPORT
int ClpMain1(std::deque<std::string> inputQueue, ClpSimplex &model,
             ampl_info *info = NULL);
CLPLIB_EXPORT
int ClpMain1(int argc, const char *argv[], ClpSimplex *model);

void printGeneralMessage(ClpSimplex &model, std::string message, int type = CLP_GENERAL);

void printGeneralWarning(ClpSimplex &model, std::string message, int type = CLP_GENERAL_WARNING);

CLPLIB_EXPORT
int clpReadAmpl(ampl_info *info, int argc, char **argv, ClpSimplex &model);
static void statistics(ClpSimplex *originalModel, ClpSimplex *model);
static void generateCode(const char *fileName, int type);
#ifdef CLP_USER_DRIVEN1
/* Returns true if variable sequenceOut can leave basis when
   model->sequenceIn() enters.
   This function may be entered several times for each sequenceOut.
   The first time realAlpha will be positive if going to lower bound
   and negative if going to upper bound (scaled bounds in lower,upper) - then will be zero.
   currentValue is distance to bound.
   currentTheta is current theta.
   alpha is fabs(pivot element).
   Variable will change theta if currentValue - currentTheta*alpha < 0.0
*/
bool userChoiceValid1(const ClpSimplex *model,
  int sequenceOut,
  double currentValue,
  double currentTheta,
  double alpha,
  double realAlpha)
{
  return true;
}
/* This returns true if chosen in/out pair valid.
   The main thing to check would be variable flipping bounds may be
   OK.  This would be signaled by reasonable theta_ and valueOut_.
   If you return false sequenceIn_ will be flagged as ineligible.
*/
bool userChoiceValid2(const ClpSimplex *model)
{
  return true;
}
/* If a good pivot then you may wish to unflag some variables.
 */
void userChoiceWasGood(ClpSimplex *model)
{
}
#endif

#endif
