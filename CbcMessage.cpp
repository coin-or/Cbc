// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CbcMessage.hpp"

typedef struct {
  CBC_Message internalNumber;
  int externalNumber; // or continuation
  char detail;
  const char * message;
} Cbc_message;
static Cbc_message us_english[]=
{
  {CBC_END_GOOD,1,1,"Search completed - best objective %g, took %d iterations and %d nodes (%.2f seconds)"},
  {CBC_MAXNODES,3,1,"Exiting on maximum nodes"},
  {CBC_MAXTIME,20,1,"Exiting on maximum time"},
  {CBC_MAXSOLS,19,1,"Exiting on maximum solutions"},
  {CBC_EVENT,27,1,"Exiting on user event"},
  {CBC_SOLUTION,4,1,"Integer solution of %g found after %d iterations and %d nodes (%.2f seconds)"},
  {CBC_END_SOLUTION,34,2,"Final check on integer solution of %g found after %d iterations and %d nodes (%.2f seconds)"},
  {CBC_SOLUTION2,33,1,"Integer solution of %g found (by alternate solver) after %d iterations and %d nodes (%.2f seconds)"},
  {CBC_END,5,1,"Partial search - best objective %g (best possible %g), took %d iterations and %d nodes (%.2f seconds)"},
  {CBC_INFEAS,6,1,"The LP relaxation is infeasible or too expensive"},
  {CBC_STRONG,7,3,"Strong branching on %d (%d), down %g (%d) up %g (%d) value %g"},
  {CBC_SOLINDIVIDUAL,8,2,"%d has value %g"},
  {CBC_INTEGERINCREMENT,9,1,"Objective coefficients multiple of %g"},
  {CBC_STATUS,10,1,"After %d nodes, %d on tree, %g best solution, best possible %g (%.2f seconds)"},
  {CBC_GAP,11,1,"Exiting as integer gap of %g less than %g or %g%%"},
  {CBC_ROUNDING,12,1,"Integer solution of %g found by heuristic after %d iterations and %d nodes (%.2f seconds)"},
  {CBC_TREE_SOL,24,1,"Integer solution of %g found by subtree after %d iterations and %d nodes (%.2f seconds)"},
  {CBC_ROOT,13,1,"At root node, %d cuts changed objective from %g to %g in %d passes"},
  {CBC_GENERATOR,14,1,"Cut generator %d (%s) - %d row cuts (%d active), %d column cuts %? in %g seconds - new frequency is %d"},
  {CBC_BRANCH,15,2,"Node %d Obj %g Unsat %d depth %d"},
  {CBC_STRONGSOL,16,1,"Integer solution of %g found by strong branching after %d iterations and %d nodes (%.2f seconds)"},
  {CBC_NOINT,3007,1,"No integer variables - nothing to do"},
  {CBC_VUB_PASS,17,1,"%d solved, %d variables fixed, %d tightened"},
  {CBC_VUB_END,18,1,"After tightenVubs, %d variables fixed, %d tightened"},
  {CBC_NOTFEAS1,21,2,"On closer inspection node is infeasible"},
  {CBC_NOTFEAS2,22,2,"On closer inspection objective value of %g above cutoff of %g"},
  {CBC_NOTFEAS3,23,2,"Allowing solution, even though largest row infeasibility is %g"},
  {CBC_CUTOFF_WARNING1,23,1,"Cutoff set to %g - equivalent to best solution of %g"},
  {CBC_ITERATE_STRONG,25,3,"%d cleanup iterations before strong branching"},
  {CBC_PRIORITY,26,1,"Setting priorities for objects %d to %d inclusive (out of %d)"},
  {CBC_WARNING_STRONG,3008,1,"Strong branching is fixing too many variables, too expensively!"},
  {CBC_START_SUB,28,1,"Starting sub-tree for %s - maximum nodes %d"},
  {CBC_END_SUB,29,1,"Ending sub-tree for %s"},
  {CBC_HEURISTIC_SOLUTION,30,1,"solution of %g found by %s after %2.f seconds"},
  {CBC_CUTS_STATS,31,1,"%d added rows had average density of %g"},
  {CBC_STRONG_STATS,32,1,"Strong branching done %d times (%d iterations), fathomed %d nodes and fixed %d variables"},
  {CBC_DUMMY_END,999999,0,""}
};
/* Constructor */
CbcMessage::CbcMessage(Language language) :
  CoinMessages(sizeof(us_english)/sizeof(Cbc_message))
{
  language_=language;
  strcpy(source_,"Cbc");
  class_ = 0; // branch and bound
  Cbc_message * message = us_english;

  while (message->internalNumber!=CBC_DUMMY_END) {
    CoinOneMessage oneMessage(message->externalNumber,message->detail,
		message->message);
    addMessage(message->internalNumber,oneMessage);
    message ++;
  }

  // now override any language ones

  switch (language) {

  default:
    message=NULL;
    break;
  }

  // replace if any found
  if (message) {
    while (message->internalNumber!=CBC_DUMMY_END) {
      replaceMessage(message->internalNumber,message->message);
      message ++;
    }
  }
}
