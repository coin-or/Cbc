// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinPragma.hpp"
#include "CglMessage.hpp"
#include <cstring>
/// Structure for use by CglMessage.cpp
typedef struct {
  CGL_Message internalNumber;
  int externalNumber; // or continuation
  char detail;
  const char *message;
} Cgl_message;
static Cgl_message us_english[] = {
  { CGL_INFEASIBLE, 0, 1, "Cut generators found to be infeasible! (or unbounded)" },
  { CGL_CLIQUES, 1, 2, "%d cliques of average size %g" },
  { CGL_FIXED, 2, 1, "%d variables fixed" },
  { CGL_PROCESS_STATS, 3, 1, "%d fixed, %d tightened bounds, %d strengthened rows, %d substitutions" },
  { CGL_SLACKS, 8, 1, "%d inequality constraints converted to equality constraints" },
  { CGL_PROCESS_STATS2, 4, 1, "processed model has %d rows, %d columns (%d integer (%d of which binary)) and %d elements" },
  { CGL_PROCESS_SOS1, 5, 1, "%d SOS with %d members" },
  { CGL_PROCESS_SOS2, 6, 1, "%d SOS (%d members out of %d) with %d overlaps - too much overlap or too many others" },
  { CGL_UNBOUNDED, 7, 1, "Continuous relaxation is unbounded!" },
  { CGL_ELEMENTS_CHANGED1, 9, 2, "%d elements changed" },
  { CGL_ELEMENTS_CHANGED2, 10, 3, "element in row %d for column %d changed from %g to %g" },
  { CGL_MADE_INTEGER, 11, 1, "%d variables made integer" },
  { CGL_ADDED_INTEGERS, 12, 1, "Added %d variables (from %d rows) with %d elements" },
  { CGL_POST_INFEASIBLE, 13, 1, "Postprocessed model is infeasible - possible tolerance issue - try without preprocessing" },
  { CGL_POST_CHANGED, 14, 1, "Postprocessing changed objective from %g to %g - possible tolerance issue - try without preprocessing" },
  { CGL_PROCESS_CLQSTR, 15, 2,"Clique Strengthening extended %ld cliques, %ld were dominated" },
  { CGL_WARNING_CLQSTR, 16, 1, "Warning: reduced costs not available in clique strengthening - changing the extension method to 'max degree'" },
  { CGL_GENERAL, 1000, 1, "%s" },
  { CGL_DUMMY_END, 999999, 0, "" }
};
/* Constructor */
CglMessage::CglMessage(Language language)
  : CoinMessages(sizeof(us_english) / sizeof(Cgl_message))
{
  language_ = language;
  strcpy(source_, "Cgl");
  class_ = 3; // Cuts
  Cgl_message *message = us_english;

  while (message->internalNumber != CGL_DUMMY_END) {
    CoinOneMessage oneMessage(message->externalNumber, message->detail,
      message->message);
    addMessage(message->internalNumber, oneMessage);
    message++;
  }
  // Put into compact form
  toCompact();
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
