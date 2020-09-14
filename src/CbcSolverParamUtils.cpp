/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/
/*
  This file is part of cbc-generic.
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include <string>
#include <cassert>

#include "CoinUtilsConfig.h"
#include "CoinFileIO.hpp"

#include "CoinFinite.hpp"
#include "CoinParam.hpp"

#include "CbcModel.hpp"

#include "CbcSolverParam.hpp"
#include "CbcSolverSettings.hpp"

/*! \file CbcSolverParamUtils
    \brief Implementation functions for CbcSolverParam parameters.
*/

namespace CbcGenSolvers {
   void setupSolverParam(CbcSolverParam &solverParam);
}

// some help strings that repeat for many options
#define CUTS_LONGHELP \
  "Value 'on' enables the cut generator and CBC will try it in the branch and "\
  "cut tree (see cutDepth on how to fine tune the behavior). Value 'root' "\
  "lets CBC run the cut generator generate only at the root node. "\
  "Value 'ifmove' lets CBC use the cut generator in the tree if it looks as "\
  "if it is doing some good and moves the objective value. Value 'forceon' "\
  "turns on the cut generator and forces CBC to use it at every node."

  
#define HEURISTICS_LONGHELP \
  "Value 'on' means to use the heuristic in each node of the tree, i.e. "\
  "after preprocessing. Value 'before' means use the heuristic only if "\
  "option doHeuristics is used. Value 'both' means to use the heuristic if "\
  "option doHeuristics is used and during solve."

namespace CbcSolverParamUtils {

/*
  Function to add CBC parameters to the CBC parameter
  vector. Where needed, defaults are drawn from cbcSettings-> This function is a
  friend of CbcSolverSettings.
*/
   
void addCbcSolverParams(CoinParamVec &parameters, CbcSolverSettings *cbcSettings)
{

  addCbcSolverStrParams(parameters, cbcSettings);
  addCbcSolverHelpParams(parameters, cbcSettings);
  addCbcSolverActionParams(parameters, cbcSettings);
  addCbcSolverKwdParams(parameters, cbcSettings);
  addCbcSolverDblParams(parameters, cbcSettings);
  addCbcSolverIntParams(parameters, cbcSettings);
  addCbcSolverBoolParams(parameters, cbcSettings);
  addCbcSolverCutParams(parameters, cbcSettings);
  addCbcSolverHeurParams(parameters, cbcSettings);

  return;
}
   
//###########################################################################
//###########################################################################
   
void addCbcSolverStrParams(CoinParamVec &parameters, CbcSolverSettings *cbcSettings)
{
  CbcSolverParam *param;

  for (int code = CbcSolverParam::CBCSOLVER_FIRSTSTRPARAM + 1;
       code < CbcSolverParam::CBCSOLVER_LASTSTRPARAM; code++){
     param = new CbcSolverParam;
     param->setObj(cbcSettings);
     param->setType(CoinParam::coinParamStr);
     parameters[code] = param;
  }

  parameters[CbcSolverParam::CSVSTATISTICS]->setup("csv!Statistics",
                                                   "Create one line of statistics",
                                                   cbcSettings->dfltDirectory_,
                                                   "This appends statistics to given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to '', i.e. it must be set.  Adds header if file empty or does not exist.",
                                                   CoinParam::displayPriorityLow);
  parameters[CbcSolverParam::CSVSTATISTICS]->setPushFunc(pushCbcSolverStrParam);

  parameters[CbcSolverParam::DEBUG]->setup("debug!In",
                                           "Read/write valid solution from/to file", "", 
                                           "This will read a solution file from the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to '', i.e. it must be set.\n\nIf set to create it will create a file called debug.file after B&C search; if set to createAfterPre it will create the file before undoing preprocessing.\n\nThe idea is that if you suspect a bad cut generator and you did not use preprocessing you can do a good run with debug set to 'create' and then switch on the cuts you suspect and re-run with debug set to 'debug.file'  Similarly if you do use preprocessing, but use createAfterPre.  The create case has the same effect as saveSolution.",
                                           CoinParam::displayPriorityNone);
  parameters[CbcSolverParam::DEBUG]->setPushFunc(doDebugParam);

  parameters[CbcSolverParam::DIRECTORY]->setup("directory",
                                               "Set Default directory for import etc.",
                                               cbcSettings->dfltDirectory_,
                                               "This sets the directory which import, export, saveModel, restoreModel etc. will use. It is initialized to the current directory.");
  parameters[CbcSolverParam::DIRECTORY]->setPushFunc(pushCbcSolverStrParam);

  parameters[CbcSolverParam::DIRSAMPLE]->setup("dirSample",
                                               "Set directory where the COIN-OR sample problems are.",
                                               cbcSettings->dfltDirectory_, "This sets the directory where the COIN-OR sample problems reside. It is used only when -unitTest is passed to clp. clp will pick up the test problems from this directory. It is initialized to '../../Data/Sample'",
                                               CoinParam::displayPriorityLow);
  parameters[CbcSolverParam::DIRSAMPLE]->setPushFunc(pushCbcSolverStrParam);

  parameters[CbcSolverParam::DIRNETLIB]->setup("dirNetlib",
                                               "Set directory where the netlib problems are.",
                                               cbcSettings->dfltDirectory_, 
                                               "This sets the directory where the netlib problems reside. One can get the netlib problems from COIN-OR or from the main netlib site. This parameter is used only when -netlib is passed to cbc. cbc will pick up the netlib problems from this directory. If clp is built without zlib support then the problems must be uncompressed. It is initialized to '../../Data/Netlib'",
                                               CoinParam::displayPriorityLow);
  parameters[CbcSolverParam::DIRNETLIB]->setPushFunc(pushCbcSolverStrParam);

   parameters[CbcSolverParam::DIRMIPLIB]->setup("dirMiplib",
                                                "Set directory where the miplib 2003 problems are.",
                                                cbcSettings->dfltDirectory_,
                                                "This sets the directory where the miplib 2003 problems reside. One can get the miplib problems from COIN-OR or from the main miplib site. This parameter is used only when -miplib is passed to cbc. cbc will pick up the miplib problems from this directory. If cbc is built without zlib support then the problems must be uncompressed. It is initialized to '../../Data/miplib3'",
                                                 CoinParam::displayPriorityLow);
  parameters[CbcSolverParam::DIRMIPLIB]->setPushFunc(pushCbcSolverStrParam);

  parameters[CbcSolverParam::MIPSTART]->setup("mips!tart",
                                              "reads an initial feasible solution from file",
                                              std::string("mipstart.sln"),
                                              "The MIPStart allows one to enter an initial integer feasible solution \
to CBC. Values of the main decision variables which are active (have \
non-zero values) in this solution are specified in a text  file. The \
text file format used is the same of the solutions saved by CBC, but \
not all fields are required to be filled. First line may contain the \
solution status and will be ignored, remaining lines contain column \
indexes, names and values as in this example:\n\
\n\
Stopped on iterations - objective value 57597.00000000\n\
      0  x(1,1,2,2)               1 \n\
      1  x(3,1,3,2)               1 \n\
      5  v(5,1)                   2 \n\
      33 x(8,1,5,2)               1 \n\
      ...\n\
\n\
Column indexes are also ignored since pre-processing can change them. \
There is no need to include values for continuous or integer auxiliary \
variables, since they can be computed based on main decision variables. \
Starting CBC with an integer feasible solution can dramatically improve \
its performance: several MIP heuristics (e.g. RINS) rely on having at \
least one feasible solution available and can start immediately if the \
user provides one. Feasibility Pump (FP) is a heuristic which tries to \
overcome the problem of taking too long to find feasible solution (or \
not finding at all), but it not always succeeds. If you provide one \
starting solution you will probably save some time by disabling FP. \
\n\n\
Knowledge specific to your problem can be considered to write an \
external module to quickly produce an initial feasible solution - some \
alternatives are the implementation of simple greedy heuristics or the \
solution (by CBC for example) of a simpler model created just to find \
a feasible solution. \
\n\n\
Silly options added.  If filename ends .low then integers not mentioned \
are set low - also .high, .lowcheap, .highcheap, .lowexpensive, .highexpensive \
where .lowexpensive sets costed ones to make expensive others low. Also if \
filename starts empty. then no file is read at all - just actions done. \
\n\n\
Question and suggestions regarding MIPStart can be directed to\n\
haroldo.santos@gmail.com. ");
  parameters[CbcSolverParam::MIPSTART]->setPushFunc(pushCbcSolverStrParam);

  parameters[CbcSolverParam::NEXTBESTSOLUTION]->setup("nextB!estSolution",
                                                      "Prints next best saved solution to file",
                                                      "",
                                                      "To write best solution, just use solution.  This prints next best (if exists) and then deletes it. This will write a primitive solution file to the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to 'stdout'.  The amount of output can be varied using printi!ngOptions or printMask.");
  parameters[CbcSolverParam::NEXTBESTSOLUTION]->setPushFunc(doNothingParam);

  parameters[CbcSolverParam::PRINTMASK]->setup("printM!ask",
                                               "Control printing of solution with a regular expression",
                                               "",
                                               "If set then only those names which match mask are printed in a solution. '?' matches any character and '*' matches any set of characters.  The default is '' (unset) so all variables are printed. This is only active if model has names.");
  parameters[CbcSolverParam::PRINTMASK]->setPushFunc(doPrintMaskParam);

  parameters[CbcSolverParam::SAVESOL]->setup("saveS!olution",
                                             "saves solution to file",
                                             std::string("solution.sln"),
                                             "This will write a binary solution file to the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to 'solution.file'.  To read the file use fread(int) twice to pick up number of rows and columns, then fread(double) to pick up objective value, then pick up row activities, row duals, column activities and reduced costs - see bottom of CbcParam.cpp for code that reads or writes file. If name contains '_fix_read_' then does not write but reads and will fix all variables");
  parameters[CbcSolverParam::SAVESOL]->setPushFunc(pushCbcSolverStrParam);

  parameters[CbcSolverParam::SOLUTION]->setup("solu!tion",
                                              "Prints solution to file",
                                              std::string("stdout"),
                                              "This will write a primitive solution file to the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to 'stdout'.  The amount of output can be varied using printi!ngOptions or printMask.");
  parameters[CbcSolverParam::SOLUTION]->setPushFunc(doSolutionParam);

  parameters[CbcSolverParam::PRIORITYIN]->setup("prio!rityIn",
                                                "Import priorities etc from file",
                                                std::string("priorities.txt"),
                                                "This will read a file with priorities from the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to '', i.e. it must be set.  This can not read from compressed files. File is in csv format with allowed headings - name, number, priority, direction, up, down, solution.  Exactly one of name and number must be given.");
  parameters[CbcSolverParam::PRIORITYIN]->setPushFunc(pushCbcSolverStrParam);

}

//###########################################################################
//###########################################################################
   
void addCbcSolverHelpParams(CoinParamVec &parameters, CbcSolverSettings *cbcSettings)
{
   CbcSolverParam *param;
   
  for (int code = CbcSolverParam::CBCSOLVER_FIRSTHELPPARAM + 1;
       code < CbcSolverParam::CBCSOLVER_LASTHELPPARAM; code++){
     param = new CbcSolverParam;
     param->setObj(cbcSettings);
     param->setType(CoinParam::coinParamAct);
     param->setPushFunc(doHelpParam);
     parameters[code] = param;
  }
  parameters[CbcSolverParam::GENERALQUERY]->setup("?",
                                                  "Print a list of commands",
                                                  CoinParam::displayPriorityNone);

  parameters[CbcSolverParam::FULLGENERALQUERY]->setup("???",
                                                      "Print a list with *all* commands, even those hidden with `?'",
                                                      CoinParam::displayPriorityNone);

  //Need display parameter to resolve ambiguity
  parameters[CbcSolverParam::HELP]->setup("help",
                                          "Print out version, non-standard options and some help",
                                          "This prints out some help to get a user started. If you're seeing this message, you should be past that stage.",
                                          CoinParam::displayPriorityHigh);

}
   
//###########################################################################
//###########################################################################
   
void addCbcSolverActionParams(CoinParamVec &parameters, CbcSolverSettings *cbcSettings)
{
  CbcSolverParam *param;
   
  for (int code = CbcSolverParam::CBCSOLVER_FIRSTACTIONPARAM + 1;
       code < CbcSolverParam::CBCSOLVER_LASTACTIONPARAM; code++){
     param = new CbcSolverParam;
     param->setObj(cbcSettings);
     param->setType(CoinParam::coinParamAct);
     parameters[code] = param;
  }
  
  parameters[CbcSolverParam::DUMMY]->setup("sleep",
                                           "for debug",
                                           0, 9999, 0,
                                           "If passed to solver from ampl, then ampl will wait so that you can copy .nl file for debug.",
                                           CoinParam::displayPriorityNone);

  parameters[CbcSolverParam::OUTDUPROWS]->setup("outDup!licates",
                                                "Takes duplicate rows, etc., out of the integer model",
                                                "", CoinParam::displayPriorityNone);
  
  parameters[CbcSolverParam::SHOWUNIMP]->setup("unimp!lemented",
                                               "Report unimplemented commands.", "",
                                               CoinParam::displayPriorityNone);
  parameters[CbcSolverParam::SHOWUNIMP]->setPushFunc(doUnimplementedParam);

  parameters[CbcSolverParam::STRENGTHEN]->setup("strengthen",
                                                "Create strengthened problem",
                                                "This creates a new problem by applying the root node cuts. All tight constraints will be in resulting problem.",
                                                 CoinParam::displayPriorityHigh);

  parameters[CbcSolverParam::BAB]->setup("branch!AndCut",
                                         "Do Branch and Cut",
                                         "This does branch and cut. There are many parameters which can affect the performance.  First just try with default cbcSettings and look carefully at the log file.  Did cuts help?  Did they take too long?  Look at output to see which cuts were effective and then do some tuning.  You will see that the options for cuts are off, on, root and ifmove.  Off is obvious, on means that this cut generator will be tried in the branch and cut tree (you can fine tune using 'depth').  Root means just at the root node while 'ifmove' means that cuts will be used in the tree if they look as if they are doing some good and moving the objective value.  If pre-processing reduced the size of the problem or strengthened many coefficients then it is probably wise to leave it on.  Switch off heuristics which did not provide solutions.  The other major area to look at is the search.  Hopefully good solutions were obtained fairly early in the search so the important point is to select the best variable to branch on.  See whether strong branching did a good job - or did it just take a lot of iterations.  Adjust the strongBranching and trustPseudoCosts parameters.",
                                         CoinParam::displayPriorityHigh);
  parameters[CbcSolverParam::BAB]->setPushFunc(doBaCParam);

  parameters[CbcSolverParam::ENVIRONMENT]->setup("environ!ment",
                                                 "Read commands from environment",
                                                 "This starts reading from environment variable COIN_ENVIRONMENT.",
                                                 CoinParam::displayPriorityNone);
  parameters[CbcSolverParam::ENVIRONMENT]->setPushFunc(doNothingParam);

  parameters[CbcSolverParam::EXIT]->setup("end",
                                          "Stops execution",
                                          "This stops execution; end, exit, quit and stop are synonyms.",
                                          CoinParam::displayPriorityHigh);
  parameters[CbcSolverParam::EXIT]->setPushFunc(doExitParam);

  parameters[CbcSolverParam::EXPORT]->setup("export",
                                            "Export model as mps file",
                                            std::string("default.mps"),
                                            "This will write an MPS format file to the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to 'default.mps'. It can be useful to get rid of the original names and go over to using Rnnnnnnn and Cnnnnnnn.  This can be done by setting 'keepnames' off before importing mps file.",
                                             CoinParam::displayPriorityHigh);
  parameters[CbcSolverParam::EXPORT]->setPushFunc(pushCbcSolverStrParam);

  parameters[CbcSolverParam::IMPORT]->setup("import",
                                            "Import model from file",
                                            cbcSettings->lastMpsIn_,
                                            "This will read an MPS format file from the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to '', i.e., it must be set.  If you have libgz then it can read compressed files 'xxxxxxxx.gz'.",
                                             CoinParam::displayPriorityHigh);
  parameters[CbcSolverParam::IMPORT]->setPushFunc(doImportParam);

  parameters[CbcSolverParam::MIPLIB]->setup("miplib", "Do some of miplib test set", "",
                                             CoinParam::displayPriorityHigh);

  parameters[CbcSolverParam::PRINTVERSION]->setup("version", "Print version", "",
                                                  CoinParam::displayPriorityHigh);
  parameters[CbcSolverParam::PRINTVERSION]->setPushFunc(doVersionParam);

  parameters[CbcSolverParam::SOLVECONTINUOUS]->setup("initialS!olve",
                                                     "Solve to continuous optimum",
                                                     "This just solves the problem to the continuous optimum, without adding any cuts.",
                                                      CoinParam::displayPriorityHigh);

  parameters[CbcSolverParam::STATISTICS]->setup( "stat!istics",
                                                 "Print some statistics",
                                                 "This command prints some statistics for the current model. If log level >1 then more is printed. These are for presolved model if presolve on (and unscaled).",
                                                  CoinParam::displayPriorityHigh);

#if 0
  // Need to figure out what to do here. Same parameter can't have two names...
  parameters[CbcSolverParam::STDIN]->setup( "-",
                                            "Switch to interactive command line mode", ""
                                            CoinParam::displayPriorityNone);
  parameters[CbcSolverParam::STDIN]->setPushFunc(doNothingParam);
#endif
  
  parameters[CbcSolverParam::STDIN]->setup("stdin",
                                           "Switch to interactive command line mode", "",
                                           CoinParam::displayPriorityNone);
  parameters[CbcSolverParam::STDIN]->setPushFunc(doNothingParam);

  parameters[CbcSolverParam::UNITTEST]->setup("unitTest", "Do unit test",
                                              "This exercises the unit test.", "",
                                               CoinParam::displayPriorityHigh);
}

//###########################################################################
//###########################################################################
   
   void addCbcSolverKwdParams(CoinParamVec &parameters, CbcSolverSettings *cbcSettings)
{
  CbcSolverParam *param;
   
  for (int code = CbcSolverParam::CBCSOLVER_FIRSTKWDPARAM + 1;
       code < CbcSolverParam::CBCSOLVER_LASTKWDPARAM; code++){

     param = new CbcSolverParam;
     param->setObj(cbcSettings);
     param->setType(CoinParam::coinParamKwd);
     param->setPushFunc(pushCbcSolverKwdParam);
     parameters[code] = param;
  }

  parameters[CbcSolverParam::ALLCOMMANDS]->setup("allC!ommands",
                             "Whether to print less used commands",
                             "off", CbcSolverParam::CMOff,
                             "For the sake of your sanity, only the more useful and simple commands are printed out on ?.",
                             CoinParam::displayPriorityNone);
  param->appendKwd("more", CbcSolverParam::CMMore);
  param->appendKwd("all", CbcSolverParam::CMAll);

  parameters[CbcSolverParam::CLQSTRENGTHENING]->setup("clqstr!engthen",
                             "Whether and when to perform Clique Strengthening preprocessing routine",
                             "after", CbcSolverParam::ClqStrAfter);
  param->appendKwd("off", CbcSolverParam::ClqStrOff);
  param->appendKwd("before", CbcSolverParam::ClqStrBefore);

  parameters[CbcSolverParam::BRANCHPRIORITY]->setup("cost!Strategy",
                                                    "Whether to use costs or column order as priorities",
                                                    "off", CbcSolverParam::BPOff,
                                                    "This orders the variables in order of their absolute costs - with largest cost ones being branched on first.  This primitive strategy can be surprisingly effective.  The column order option is obviously not on costs but it's easy to implement.");
  param->appendKwd("pri!orities", CbcSolverParam::BPCost);
  param->appendKwd("column!Order", CbcSolverParam::BPOrder);

  parameters[CbcSolverParam::CUTOFFCONSTRAINT]->setup("constraint!fromCutoff",
                                                      "Whether to use cutoff as constraint",
                                                      "off", CbcSolverParam::COOff,
                                                      "For some problems, cut generators and general branching work better if the problem would be infeasible if the cost is too high. "
      "If this option is enabled, the objective function is added as a constraint which right hand side is set to the current cutoff value (objective value of best known solution)");
  param->appendKwd("on", CbcSolverParam::COOn);
  param->appendKwd("variable", CbcSolverParam::COVariable);
  param->appendKwd("forcevariable", CbcSolverParam::COForceVariable);
  param->appendKwd("conflict", CbcSolverParam::COConflict);

  parameters[CbcSolverParam::INTPRINT]->setup("printi!ngOptions",
                                              "Print options",
                                              "normal", CbcSolverParam::PMNormal,
                                              "This changes the amount and format of printing a solution:\nnormal - nonzero column variables \ninteger - nonzero integer column variables\nspecial - in format suitable for OsiRowCutDebugger\nrows - nonzero column variables and row activities\nall - all column variables and row activities.\n\nFor non-integer problems 'integer' and 'special' act like 'normal'.  Also see printMask for controlling output.");
  param->appendKwd("integer", CbcSolverParam::PMInteger);
  param->appendKwd("special", CbcSolverParam::PMSpecial);
  param->appendKwd("rows", CbcSolverParam::PMRows);
  param->appendKwd("all", CbcSolverParam::PMAll);

  parameters[CbcSolverParam::NODESTRATEGY]->setup("node!Strategy",
                                                  "What strategy to use to select the next node from the branch and cut tree",
                                                  "hybrid", CbcSolverParam::NSHybrid,
                                                  "Normally before a feasible solution is found, CBC will choose a node with fewest infeasibilities. Alternatively, one may choose tree-depth as the criterion. This requires the minimal amount of memory, but may take a long time to find the best solution. Additionally, one may specify whether up or down branches must be selected first (the up-down choice will carry on after a first solution has been bound). The choice 'hybrid' does breadth first on small depth nodes and then switches to 'fewest'.");
  param->appendKwd("fewest", CbcSolverParam::NSFewest);
  param->appendKwd("depth", CbcSolverParam::NSDepth);
  param->appendKwd("upfewest", CbcSolverParam::NSUpFewest);
  param->appendKwd("downfewest", CbcSolverParam::NSDownFewest);
  param->appendKwd("updepth", CbcSolverParam::NSUpDepth);
  param->appendKwd("downdepth", CbcSolverParam::NSDownDepth);

  parameters[CbcSolverParam::ORBITAL]->setup("Orbit!alBranching",
                                             "Whether to try orbital branching",
                                             "off", CbcSolverParam::OBOff,
                                             "This switches on Orbital branching. Value 'on' just adds orbital, 'strong' tries extra fixing in strong branching.");
  param->appendKwd("on", CbcSolverParam::OBOn);
  param->appendKwd("slowish", CbcSolverParam::OBSlowish);
  param->appendKwd("strong", CbcSolverParam::OBStrong);
  param->appendKwd("force", CbcSolverParam::OBForce);
  param->appendKwd("simple", CbcSolverParam::OBSimple);
  param->appendKwd("moreprinting", CbcSolverParam::OBMorePrinting);

  parameters[CbcSolverParam::PREPROCESS]->setup("preprocess",
                                                "Whether to use integer preprocessing",
                                                "off", CbcSolverParam::IPPOff,
                                                "This tries to reduce size of the model in a similar way to presolve and it also tries to strengthen the model. This can be very useful and is worth trying.  save option saves on file presolved.mps.  equal will turn <= cliques into ==.  sos will create sos sets if all 0-1 in sets (well one extra is allowed) and no overlaps.  trysos is same but allows any number extra. equalall will turn all valid inequalities into equalities with integer slacks. strategy is as on but uses CbcStrategy.");
  param->appendKwd("on", CbcSolverParam::IPPOn);
  param->appendKwd("save", CbcSolverParam::IPPSave);
  param->appendKwd("equal", CbcSolverParam::IPPEqual);
  param->appendKwd("sos", CbcSolverParam::IPPSOS);
  param->appendKwd("trysos", CbcSolverParam::IPPTrySOS);
  param->appendKwd("equalall", CbcSolverParam::IPPEqualAll);
  param->appendKwd("strategy", CbcSolverParam::IPPStrategy);

  parameters[CbcSolverParam::SOSPRIORITIZE]->setup("sosP!rioritize",
                                                   "How to deal with SOS priorities",
                                                   "off", CbcSolverParam::SOSOff,
                                                   "This sets priorities for SOS.  Values 'high' and 'low' just set a priority relative to the for integer variables.  Value 'orderhigh' gives first highest priority to the first SOS and integer variables a low priority.  Value 'orderlow' gives integer variables a high priority then SOS in order.");
  param->appendKwd("high", CbcSolverParam::SOSHigh);
  param->appendKwd("low", CbcSolverParam::SOSLow);
  param->appendKwd("orderhigh", CbcSolverParam::SOSOrderHigh);
  param->appendKwd("orderlow", CbcSolverParam::SOSOrderLow);

  parameters[CbcSolverParam::STRATEGY]->setup("strat!egy",
                                              "Switches on groups of features",
                                              "default", CbcSolverParam::StrategyDefault,
                                              "This turns on newer features. Use 0 for easy problems, 1 is default, 2 is aggressive. 1 uses Gomory cuts with a tolerance of 0.01 at the root node, does a possible restart after 100 nodes if many variables could be fixed, activates a diving and RINS heuristic, and makes the feasibility pump more aggressive."); // This does not apply to unit tests (where 'experiment' may have similar effects)
  param->appendKwd("easy", CbcSolverParam::StrategyEasy);
  param->appendKwd("aggressive", CbcSolverParam::StrategyAggressive);

  parameters[CbcSolverParam::TIMEMODE]->setup("timeM!ode",
                                              "Whether to use CPU or elapsed time",
                                              "cpu", CbcSolverParam::ClockCpu,
                                              "cpu uses CPU time for stopping, while elapsed uses elapsed time. (On Windows, elapsed time is always used).");
  param->appendKwd("elapsed", CbcSolverParam::ClockElapsed);

  parameters[CbcSolverParam::USECGRAPH]->setup("cgraph",
                                               "Whether to use the conflict graph-based preprocessing and cut separation routines.",
                                               "on", CbcSolverParam::CGraphOn,
                                               "This switches the conflict graph-based preprocessing and cut separation routines (CglBKClique, CglOddWheel and CliqueStrengthening) on or off. Values:\
\n\toff: turns these routines off;\
\n\ton: turns these routines on;\
\n\tclq: turns these routines off and enables the cut separator of CglClique.");
  param->appendKwd("off", CbcSolverParam::CGraphOff);
  param->appendKwd("clq", CbcSolverParam::CGraphClique);

}
   
//###########################################################################
//###########################################################################
   
void addCbcSolverDblParams(CoinParamVec &parameters, CbcSolverSettings *cbcSettings)
{
  CbcSolverParam *param;
   
  for (int code = CbcSolverParam::CBCSOLVER_FIRSTDBLPARAM + 1;
       code < CbcSolverParam::CBCSOLVER_LASTDBLPARAM; code++){
     param = new CbcSolverParam;
     param->setObj(cbcSettings);
     param->setType(CoinParam::coinParamDbl);
     param->setPushFunc(pushCbcSolverDblParam);
     parameters[code] = param;
  }

  parameters[CbcSolverParam::ARTIFICIALCOST]->setup("artif!icialCost",
                                                    "Costs >= this treated as artificials in feasibility pump",
                                                    0.0, COIN_DBL_MAX,
                                                    cbcSettings->getArtVarThreshold(), "",
                                                    CoinParam::displayPriorityLow);

    parameters[CbcSolverParam::DEXTRA3]->setup("dextra3",
                                               "Extra double parameter 3",
                                               -COIN_DBL_MAX, COIN_DBL_MAX, 0.0, "",
                                               CoinParam::displayPriorityNone);

  parameters[CbcSolverParam::DEXTRA4]->setup("dextra4",
                                             "Extra double parameter 4",
                                             -COIN_DBL_MAX, COIN_DBL_MAX, 0.0, "",
                                             CoinParam::displayPriorityNone);

  parameters[CbcSolverParam::DEXTRA5]->setup("dextra5",
                                             "Extra double parameter 5",
                                             -COIN_DBL_MAX, COIN_DBL_MAX, 0.0, "",
                                             CoinParam::displayPriorityNone);
 
  parameters[CbcSolverParam::DJFIX]->setup("fix!OnDj",
                                           "Try heuristic that fixes variables based on reduced costs",
                                           -1.0e20, 1.0e20, cbcSettings->getDjFixThreshold(),
                                           "If set, integer variables with reduced costs greater than the specified value will be fixed before branch and bound - use with extreme caution!");

  parameters[CbcSolverParam::FAKECUTOFF]->setup("pumpC!utoff",
                                                "Fake cutoff for use in feasibility pump",
                                                -COIN_DBL_MAX, COIN_DBL_MAX, 0.0,
                                                "A value of 0.0 means off. Otherwise, add a constraint forcing objective below this value in feasibility pump",
                                                CoinParam::displayPriorityLow);

  parameters[CbcSolverParam::FAKEINCREMENT]->setup("pumpI!ncrement",
                                                   "Fake increment for use in feasibility pump",
                                                   -COIN_DBL_MAX, COIN_DBL_MAX, 0.0,
                                                   "A value of 0.0 means off. Otherwise, add a constraint forcing objective below this value in feasibility pump",
                                                   CoinParam::displayPriorityLow);

  parameters[CbcSolverParam::SMALLBAB]->setup("fraction!forBAB",
                                              "Fraction in feasibility pump",
                                              1.0e-5, 1.1, 0.5,
                                              "After a pass in the feasibility pump, variables which have not moved about are fixed and if the preprocessed model is smaller than this fraction of the original problem, a few nodes of branch and bound are done on the reduced problem.",
                                              CoinParam::displayPriorityLow);

  parameters[CbcSolverParam::TIGHTENFACTOR]->setup("tighten!Factor",
                                                   "Tighten bounds using value times largest activity at continuous solution",
                                                   0.0, COIN_DBL_MAX, 0.0,
                                                   "This sleazy trick can help on some problems.");


}
   
//###########################################################################
//###########################################################################
   
void addCbcSolverIntParams(CoinParamVec &parameters, CbcSolverSettings *cbcSettings)
{
  CbcSolverParam *param;
   
  for (int code = CbcSolverParam::CBCSOLVER_FIRSTINTPARAM + 1;
       code < CbcSolverParam::CBCSOLVER_LASTINTPARAM; code++){
     param = new CbcSolverParam;
     param->setObj(cbcSettings);
     param->setType(CoinParam::coinParamInt);
     param->setPushFunc(pushCbcSolverIntParam);
     parameters[code] = param;
  }
   
  parameters[CbcSolverParam::BKPIVOTINGSTRATEGY]->setup("bkpivot!ing",
                                                        "Pivoting strategy used in Bron-Kerbosch algorithm",
                                                        0, 6, 3);

  parameters[CbcSolverParam::BKMAXCALLS]->setup("bkmaxcalls",
                                                "Maximum number of recursive calls made by Bron-Kerbosch algorithm",
                                                1, COIN_INT_MAX, 1000);

  parameters[CbcSolverParam::BKCLQEXTMETHOD]->setup("bkclqext!method",
                                                    "Strategy used to extend violated cliques found by BK Clique Cut Separation routine",
                                                    0, 5, 4,
                                                    "Sets the method used in the extension module of BK Clique Cut Separation routine: 0=no extension; 1=random; 2=degree; 3=modified degree; 4=reduced cost(inversely proportional); 5=reduced cost(inversely proportional) + modified degree");

  parameters[CbcSolverParam::CPP]->setup("cpp!Generate",
                                         "Generates C++ code",
                                         0, 4, 0, "Once you like what the stand-alone solver does then this allows you to generate user_driver.cpp which approximates the code.  0 gives simplest driver, 1 generates saves and restores, 2 generates saves and restores even for variables at default value. 4 bit in cbc generates size dependent code rather than computed values.");

  parameters[CbcSolverParam::CUTDEPTH]->setup("cutD!epth",
                                              "Depth in tree at which to do cuts",
                                              -1, 999999, cbcSettings->getCutDepth(),
                                              "Cut generators may be off, on only at the root, on if they look useful, and on at some interval.  If they are done every node then that is that, but it may be worth doing them every so often.  The original method was every so many nodes but it is more logical to do it whenever depth in tree is a multiple of K.  This option does that and defaults to -1 (off).");

  parameters[CbcSolverParam::CUTLENGTH]->setup("cutL!ength",
                                               "Length of a cut",
                                               -1, COIN_INT_MAX, -1,
                                               "At present this only applies to Gomory cuts. -1 (default) leaves as is. Any value >0 says that all cuts <= this length can be generated both at root node and in tree. 0 says to use some dynamic lengths.  If value >=10,000,000 then the length in tree is value%10000000 - so 10000100 means unlimited length at root and 100 in tree.");

  parameters[CbcSolverParam::CUTPASSINTREE]->setup("passT!reeCuts",
                                                   "Number of rounds that cut generators are applied in the tree",
                                                   -COIN_INT_MAX, COIN_INT_MAX);

  parameters[CbcSolverParam::DEPTHMINIBAB]->setup("depth!MiniBab",
                                                  "Depth at which to try mini branch-and-bound",
                                                  -COIN_INT_MAX, COIN_INT_MAX, -1,
                                                  "Rather a complicated parameter but can be useful. -1 means off for large problems but on as if -12 for problems where rows+columns<500, -2 means use Cplex if it is linked in.  Otherwise if negative then go into depth first complete search fast branch and bound when depth>= -value-2 (so -3 will use this at depth>=1).  This mode is only switched on after 500 nodes.  If you really want to switch it off for small problems then set this to -999.  If >=0 the value doesn't matter very much.  The code will do approximately 100 nodes of fast branch and bound every now and then at depth>=5. The actual logic is too twisted to describe here.");

   parameters[CbcSolverParam::DIVEOPT]->setup("diveO!pt",
                                              "Diving options",
                                              -1, 20, -1,"If >2 && <20 then modify diving options - \
	 \n\t3 only at root and if no solution,  \
	 \n\t4 only at root and if this heuristic has not got solution, \
	 \n\t5 decay only if no solution, \
	 \n\t6 if depth <3 or decay, \
	 \n\t7 run up to 2 times if solution found 4 otherwise, \
	 \n\t>10 All only at root (DivingC normal as value-10), \
	 \n\t>20 All with value-20).",
                                              CoinParam::displayPriorityLow);

  parameters[CbcSolverParam::DIVEOPTSOLVES]->setup("diveS!olves",
                                                   "Diving solve option",
                                                   -1, 200000, 100, "If >0 then do up to this many solves. However, the last digit is ignored and used for extra options: 1-3 enables fixing of satisfied integer variables (but not at bound), where 1 switches this off for that dive if the dive goes infeasible, and 2 switches it off permanently if the dive goes infeasible.",
                                                   CoinParam::displayPriorityLow);

  parameters[CbcSolverParam::EXPERIMENT]->setup("exper!iment",
                                                "Whether to use testing features",
                                                -1, 200000, 0,
                                                "Defines how adventurous you want to be in using new ideas. 0 then no new ideas, 1 fairly sensible, 2 a bit dubious, 3 you are on your own!",
                                                CoinParam::displayPriorityLow);

  parameters[CbcSolverParam::EXTRA1]->setup("extra1",
                                            "Extra integer parameter 1",
                                            -COIN_INT_MAX, COIN_INT_MAX, -1, "",
                                            CoinParam::displayPriorityLow);

  parameters[CbcSolverParam::EXTRA2]->setup("extra2",
                                            "Extra integer parameter 2",
                                            -COIN_INT_MAX, COIN_INT_MAX, -1, "",
                                            CoinParam::displayPriorityLow);

  parameters[CbcSolverParam::EXTRA3]->setup("extra3",
                                            "Extra integer parameter 3",
                                            -COIN_INT_MAX, COIN_INT_MAX, -1, "",
                                            CoinParam::displayPriorityLow);

  parameters[CbcSolverParam::EXTRA4]->setup("extra4",
                                            "Extra integer parameter 4",
                                            -COIN_INT_MAX, COIN_INT_MAX, -1, "",
                                            CoinParam::displayPriorityLow);

  parameters[CbcSolverParam::FPUMPITS]->setup("passF!easibilityPump",
                                              "How many passes in feasibility pump",
                                              0, 10000, cbcSettings->getFeasPumpIters(),
                                              "This fine tunes the Feasibility Pump heuristic by doing more or fewer passes.");

  parameters[CbcSolverParam::FPUMPTUNE]->setup("pumpT!une",
                                               "Dubious ideas for feasibility pump",
                                               0, 100000000, 0,
                                               "This fine tunes Feasibility Pump \
     \n\t>=10000000 use as objective weight switch \
     \n\t>=1000000 use as accumulate switch \
     \n\t>=1000 use index+1 as number of large loops \
     \n\t==100 use objvalue +0.05*fabs(objvalue) as cutoff OR fakeCutoff if set \
     \n\t%100 == 10,20 affects how each solve is done \
     \n\t1 == fix ints at bounds, 2 fix all integral ints, 3 and continuous at bounds. If accumulate is on then after a major pass, variables which have not moved are fixed and a small branch and bound is tried.");

  parameters[CbcSolverParam::FPUMPTUNE2]->setup("moreT!une",
                                                "Yet more dubious ideas for feasibility pump",
                                                0, 100000000, 0, 
                                                "Yet more ideas for Feasibility Pump \
     \n\t/100000 == 1 use box constraints and original obj in cleanup \
     \n\t/1000 == 1 Pump will run twice if no solution found \
     \n\t/1000 == 2 Pump will only run after root cuts if no solution found \
     \n\t/1000 >10 as above but even if solution found \
     \n\t/100 == 1,3.. exact 1.0 for objective values \
     \n\t/100 == 2,3.. allow more iterations per pass \
     \n\t n fix if value of variable same for last n iterations.",
                                                CoinParam::displayPriorityNone);

  parameters[CbcSolverParam::HEUROPTIONS]->setup("hOp!tions",
                                                 "Heuristic options",
                                                 -COIN_INT_MAX, COIN_INT_MAX, 0,
                                                 "Value 1 stops heuristics immediately if the allowable gap has been reached. Other values are for the feasibility pump - 2 says do exact number of passes given, 4 only applies if an initial cutoff has been given and says relax after 50 passes, while 8 will adapt the cutoff rhs after the first solution if it looks as if the code is stalling.",
                                                 CoinParam::displayPriorityLow);

  parameters[CbcSolverParam::LOGLEVEL]->setup("log!Level",
                                              "Level of detail in CBC output.",
                                              -1, 999999, cbcSettings->getLogLevel(),
                                              "If set to 0 then there should be no output in normal circumstances. A value of 1 is probably the best value for most uses, while 2 and 3 give more information.");

  parameters[CbcSolverParam::LPLOGLEVEL]->setup("log!Level",
                                                "Level of detail in LP solver output.",
                                                -1, 999999, cbcSettings->getLpLogLevel(),
                                                "If set to 0 then there should be no output in normal circumstances. A value of 1 is probably the best value for most uses, while 2 and 3 give more information.");

  parameters[CbcSolverParam::MAXSAVEDSOLS]->setup("maxSaved!Solutions",
                                                  "Maximum number of solutions to save",
                                                  0, COIN_INT_MAX, 1,
                                                  "Number of solutions to save.");

  parameters[CbcSolverParam::MAXSLOWCUTS]->setup("slow!cutpasses",
                                                 "Maximum number of rounds for slower cut generators",
                                                 -1, COIN_INT_MAX, 10,
                                                 "Some cut generators are fairly slow - this limits the number of times they are tried. The cut generators identified as 'may be slow' at present are Lift and project cuts and both versions of Reduce and Split cuts.");

  parameters[CbcSolverParam::MOREMOREMIPOPTIONS]->setup("more2!MipOptions",
                                                        "More more dubious options for mip",
                                                        -1, COIN_INT_MAX, 0, "",
                                                        CoinParam::displayPriorityNone);

  parameters[CbcSolverParam::MULTIPLEROOTS]->setup("multiple!RootPasses",
                                                   "Do multiple root passes to collect cuts and solutions",
                                                   0, COIN_INT_MAX, 0,   
                                                   "Solve (in parallel, if enabled) the root phase this number of times, each with its own different seed, and collect all solutions and cuts generated. The actual format is aabbcc where aa is the number of extra passes; if bb is non zero, then it is number of threads to use (otherwise uses threads setting); and cc is the number of times to do root phase. The solvers do not interact with each other.  However if extra passes are specified then cuts are collected and used in later passes - so there is interaction there. Some parts of this implementation have their origin in idea of Andrea Lodi, Matteo Fischetti, Michele Monaci, Domenico Salvagnin, and Andrea Tramontani.",
                                                   CoinParam::displayPriorityNone);

  parameters[CbcSolverParam::ODDWEXTMETHOD]->setup("oddwext!method",
                                                   "Strategy used to search for wheel centers for the cuts found by Odd Wheel Cut Separation routine",
                                                   0, 2, 2,
                                                   "Sets the method used in the extension module of Odd Wheel Cut Separation routine: 0=no extension; 1=one variable; 2=clique");

  parameters[CbcSolverParam::OUTPUTFORMAT]->setup("output!Format",
                                                  "Which output format to use",
                                                  1, 6, 2,
                                                  "Normally export will be done using normal representation for numbers and two values per line.  You may want to do just one per line (for grep or suchlike) and you may wish to save with absolute accuracy using a coded version of the IEEE value. A value of 2 is normal. Otherwise, odd values give one value per line, even values two.  Values of 1 and 2 give normal format, 3 and 4 give greater precision, 5 and 6 give IEEE values.  When exporting a basis, 1 does not save values, 2 saves values, 3 saves with greater accuracy and 4 saves in IEEE format.");

  parameters[CbcSolverParam::PRINTOPTIONS]->setup("pO!ptions",
                                                  "Dubious print options",
                                                  0, COIN_INT_MAX, 0,
                                                  "If this is greater than 0 then presolve will give more information and branch and cut will give statistics");

    parameters[CbcSolverParam::PROCESSTUNE]->setup("tune!PreProcess",
                                                   "Dubious tuning parameters for preprocessing",
                                                   0, COIN_INT_MAX, 0,
                                                   "Format aabbcccc - \n If aa then this is number of major passes (i.e. with presolve) \n \
If bb and bb>0 then this is number of minor passes (if unset or 0 then 10) \n \
cccc is bit set \n 0 - 1 Heavy probing \n 1 - 2 Make variables integer if possible (if obj value)\n \
2 - 4 As above but even if zero objective value\n \
7 - 128 Try and create cliques\n 8 - 256 If all +1 try hard for dominated rows\n \
9 - 512 Even heavier probing \n \
10 - 1024 Use a larger feasibility tolerance in presolve\n \
11 - 2048 Try probing before creating cliques\n \
12 - 4096 Switch off duplicate column checking for integers \n \
13 - 8192 Allow scaled duplicate column checking \n \n \
     Now aa 99 has special meaning i.e. just one simple presolve.",
                                                   CoinParam::displayPriorityLow);

  parameters[CbcSolverParam::RANDOMSEED]->setup("randomC!bcSeed",
                                                "Random seed for Cbc",
                                                -1, COIN_INT_MAX, -1,
                                                "Allows initialization of the random seed for pseudo-random numbers used in heuristics such as the Feasibility Pump to decide whether to round up or down. The special value of 0 lets Cbc use the time of the day for the initial seed.");

  parameters[CbcSolverParam::STRONGSTRATEGY]->setup("expensive!Strong",
                                                    "Whether to do even more strong branching",
                                                    0, COIN_INT_MAX, 0,
                                                    "Strategy for extra strong branching. 0 is normal strong branching. 1, 2, 4, and 6 does strong branching on all fractional variables if at the root node (1), at depth less than modifier (2), objective equals best possible (4), or at depth less than modifier and objective equals best possible (6). 11, 12, 14, and 16 are like 1, 2, 4, and 6, respecitively, but do strong branching on all integer (incl. non-fractional) variables. Values >= 100 are used to specify a depth limit (value/100), otherwise 5 is used. If the values >= 100, then above rules are applied to value%100.",
                                                    CoinParam::displayPriorityNone);

  parameters[CbcSolverParam::TESTOSI]->setup("testO!si", 
                                             "Test OsiObject stuff",
                                             -1, COIN_INT_MAX, -1, "",
                                             CoinParam::displayPriorityNone);

#ifdef CBC_THREAD
  parameters[CbcSolverParam::THREADS]->setup("thread!s",
                                             "Number of threads to try and use",
                                             -100, 100000, 0,
                                             "To use multiple threads, set threads to number wanted.  It may be better to use one or two more than number of cpus available.  If 100+n then n threads and search is repeatable (maybe be somewhat slower), if 200+n use threads for root cuts, 400+n threads used in sub-trees.",
                                             CoinParam::displayPriorityLow);
#endif
    
  parameters[CbcSolverParam::USERCBC]->setup("userCbc",
                                             "Hand coded Cbc stuff",
                                             0, COIN_INT_MAX, 0,
                                             "There are times (e.g., when using AMPL interface) when you may wish to do something unusual.  Look for USERCBC in main driver and modify sample code.",
                                             CoinParam::displayPriorityNone);

  parameters[CbcSolverParam::VERBOSE]->setup("verbose",
                                             "Switches on longer help on single ?",
                                             0, 15, cbcSettings->verbose_,
                                             "Set to 1 to get short help with ? list, 2 to get long help.",
                                             CoinParam::displayPriorityNone);

  parameters[CbcSolverParam::VUBTRY]->setup("vub!heuristic",
                                            "Type of VUB heuristic",
                                            -2, 20, -1,
                                            "This heuristic tries to fix some integer variables.",
                                            CoinParam::displayPriorityNone);
}

//###########################################################################
//###########################################################################
   
void addCbcSolverBoolParams(CoinParamVec &parameters, CbcSolverSettings *cbcSettings)
{
  CbcSolverParam *param;
   
  // Populate the keyword lists
  for (int code = CbcSolverParam::CBCSOLVER_FIRSTBOOLPARAM + 1;
       code < CbcSolverParam::CBCSOLVER_LASTBOOLPARAM; code++){
     param = new CbcSolverParam;
     param->setObj(cbcSettings);
     param->setType(CoinParam::coinParamKwd);
     param->setPushFunc(pushCbcSolverKwdParam);
     parameters[code] = param;
  }

  parameters[CbcSolverParam::CPX]->setup("cplex!Use",
                                         "Whether to use Cplex!",
                                         "off", 0,
                                         "If the user has Cplex, but wants to use some of Cbc's heuristics then you can!  If this is on, then Cbc will get to the root node and then hand over to Cplex.  If heuristics find a solution this can be significantly quicker.  You will probably want to switch off Cbc's cuts as Cplex thinks they are genuine constraints.  It is also probable that you want to switch off preprocessing, although for difficult problems it is worth trying both.");

  parameters[CbcSolverParam::DOHEURISTIC]->setup("doH!euristic",
                                                 "Do heuristics before any preprocessing",
                                                 "off", cbcSettings->getDiveCoefficientMode(),
                                                 "Normally heuristics are done in branch and bound.  It may be useful to do them outside. Only those heuristics with 'both' or 'before' set will run. Doing this may also set cutoff, which can help with preprocessing.");

  parameters[CbcSolverParam::ERRORSALLOWED]->setup("error!sAllowed",
                                                   "Whether to allow import errors",
                                                   "off", 0, "The default is not to use any model which had errors when reading the mps file.  Setting this to 'on' will allow all errors from which the code can recover simply by ignoring the error.  There are some errors from which the code can not recover, e.g., no ENDATA.  This has to be set before import, i.e., -errorsAllowed on -import xxxxxx.mps.");

  parameters[CbcSolverParam::EXTRAVARIABLES]->setup("extraV!ariables",
                                                    "Allow creation of extra integer variables",
                                                    "off", 0,
                                                    "Switches on a trivial re-formulation that introduces extra integer variables to group together variables with same cost.",
                                                    CoinParam::displayPriorityLow);

  parameters[CbcSolverParam::MESSAGES]->setup("mess!ages",
                                              "Controls whether standardised message prefix is printed",
                                              "off", 0, "By default, messages have a standard prefix, such as:\n   Clp0005 2261  Objective 109.024 Primal infeas 944413 (758)\nbut this program turns this off to make it look more friendly.  It can be useful to turn them back on if you want to be able to 'grep' for particular messages or if you intend to override the behavior of a particular message.");

  parameters[CbcSolverParam::PREPROCNAMES]->setup("PrepN!ames",
                                                  "If column names will be kept in pre-processed model",
                                                  "off", 0, "Normally the preprocessed model has column names replaced by new names C0000... Setting this option to on keeps original names in variables which still exist in the preprocessed problem");

  parameters[CbcSolverParam::SOS]->setup("sos!Options",
                                         "Whether to use SOS from AMPL",
                                         "off", 0,
                                         "Normally if AMPL says there are SOS variables they should be used, but sometimes they should be turned off - this does so.");

  parameters[CbcSolverParam::USESOLUTION]->setup("force!Solution",
                                                 "Whether to use given solution as crash for BAB",
                                                 "off", 0,
                                                 "If on then tries to branch to solution given by AMPL or priorities file.");

}

//###########################################################################
//###########################################################################
   
void addCbcSolverCutParams(CoinParamVec &parameters, CbcSolverSettings *cbcSettings)
{
   CbcSolverParam *param;
   
   for (int code = CbcSolverParam::CBCSOLVER_FIRSTCUTPARAM + 1;
        code < CbcSolverParam::CBCSOLVER_LASTCUTPARAM; code++){
     param = new CbcSolverParam;
     param->setObj(cbcSettings);
     param->setType(CoinParam::coinParamKwd);
     param->setPushFunc(pushCbcSolverKwdParam);
     parameters[code] = param;
   }
     
   parameters[CbcSolverParam::CLIQUECUTS]->setup("clique!Cuts",
                                                 "Whether to use clique cuts",
                                                 "off", CbcSolverParam::CGOff,
                                                 "This switches on clique cuts (either at root or in entire tree). See branchAndCut for information on options.");
   
   parameters[CbcSolverParam::CUTSTRATEGY]->setup("cuts!OnOff",
                                                  "Switches all cuts on or off",
                                                  "off", CbcSolverParam::CGOff,
                                                  "This can be used to switch on or off all cuts (apart from Reduce and Split).  Then you can set individual ones off or on.  See branchAndCut for information on options.");

   parameters[CbcSolverParam::FLOWCUTS]->setup("flow!CoverCuts",
                                               "Whether to use Flow Cover cuts",
                                               "off", cbcSettings->getFlowMode(),
                                               "This switches on flow cover cuts (either at root or in entire tree)." CUTS_LONGHELP);

  parameters[CbcSolverParam::GMICUTS]->setup("GMI!Cuts",
                                             "Whether to use alternative Gomory cuts",
                                             "off", cbcSettings->getGMIMode(),
                                             CUTS_LONGHELP " This version is by Giacomo Nannicini and may be more robust than gomoryCuts.");

  parameters[CbcSolverParam::GOMORYCUTS]->setup("gomory!Cuts",
                                                "Whether to use Gomory cuts",
                                                "off", cbcSettings->getGomoryMode(),
                                                "The original cuts - beware of imitations!  Having gone out of favor, they are now more fashionable as LP solvers are more robust and they interact well with other cuts.  They will almost always give cuts (although in this executable they are limited as to number of variables in cut).  However the cuts may be dense so it is worth experimenting (Long allows any length). " CUTS_LONGHELP " Reference: https://github.com/coin-or/Cgl/wiki/CglGomory");

  parameters[CbcSolverParam::KNAPSACKCUTS]->setup("knapsack!Cuts",
                                                  "Whether to use Knapsack cuts",
                                                  "off", cbcSettings->getKnapsackMode(),
                                                  "This switches on knapsack cuts (either at root or in entire tree). See branchAndCut for information on options.");

  parameters[CbcSolverParam::LANDPCUTS]->setup("lift!AndProjectCuts",
                                               "Whether to use lift-and-project cuts",
                                               "off", cbcSettings->getLandPMode(),
                                               "This switches on lift-and-project cuts (either at root or in entire tree). See branchAndCut for information on options.");

  parameters[CbcSolverParam::LAGOMORYCUTS]->setup("lagomory!Cuts",
                                                  "Whether to use Lagrangean Gomory cuts",
                                                  "off", cbcSettings->getLaGomoryMode(),
                                                  "This is a gross simplification of 'A Relax-and-Cut Framework for Gomory's Mixed-Integer Cuts' by Matteo Fischetti & Domenico Salvagnin.  This simplification just uses original constraints while modifying objective using other cuts. So you don't use messy constraints generated by Gomory etc. A variant is to allow non messy cuts e.g. clique cuts. So 'only' does this while 'clean' also allows integral valued cuts.  'End' is recommended and waits until other cuts have finished before it does a few passes. The length options for gomory cuts are used.");

  parameters[CbcSolverParam::LATWOMIRCUTS]->setup("latwomir!Cuts",
                                                  "Whether to use Lagrangean Twomir cuts",
                                                  "off", cbcSettings->getLaTwomirMode(),
                                                  "This is a Lagrangean relaxation for Twomir cuts.  See lagomoryCuts for description of options.");

  parameters[CbcSolverParam::MIRCUTS]->setup("mixed!IntegerRoundingCuts",
                                             "Whether to use Mixed Integer Rounding cuts",
                                             "off", cbcSettings->getMirMode(),
                                             "This switches on mixed integer rounding cuts (either at root or in entire tree).  See branchAndCut for information on options.");

  parameters[CbcSolverParam::ODDWHEELCUTS]->setup("oddwheel!Cuts",
                                                  "Whether to use odd wheel cuts",
                                                  "off", cbcSettings->getOddWheelMode(),
                                                  "This switches on odd-wheel inequalities (either at root or in entire tree).");

  parameters[CbcSolverParam::PROBINGCUTS]->setup("probing!Cuts",
                                                 "Whether to use Probing cuts", "off",
                                                 cbcSettings->getProbingMode(),
                                                 "This switches on probing cuts (either at root or in entire tree). See branchAndCut for information on options.");

  parameters[CbcSolverParam::REDSPLITCUTS]->setup("reduce!AndSplitCuts",
                                                  "Whether to use Reduce-and-Split cuts",
                                                  "off", cbcSettings->getRedSplitMode(),
                                                  "This switches on reduce and split cuts (either at root or in entire tree). See branchAndCut for information on options.");

  parameters[CbcSolverParam::REDSPLIT2CUTS]->setup("reduce2!AndSplitCuts",
                                                   "Whether to use Reduce-and-Split cuts - style 2",
                                                   "off", cbcSettings->getRedSplit2Mode(),
                                                   "This switches on reduce and split cuts (either at root or in entire tree). See branchAndCut for information on options.");

  parameters[CbcSolverParam::RESIDCAPCUTS]->setup("residual!CapacityCuts",
                                                  "Whether to use Residual Capacity cuts",
                                                  "off", cbcSettings->getResidCapMode(),
                                                  CUTS_LONGHELP " Reference: https://github.com/coin-or/Cgl/wiki/CglResidualCapacity");

  parameters[CbcSolverParam::TWOMIRCUTS]->setup("two!MirCuts",
                                                "Whether to use Two phase Mixed Integer Rounding cuts",
                                                "off", cbcSettings->getTwomirMode(),
                                                "This switches on two phase mixed integer rounding cuts (either at root or in entire tree). See branchAndCut for information on options.");

  parameters[CbcSolverParam::ZEROHALFCUTS]->setup("zero!HalfCuts",
                                                  "Whether to use zero half cuts",
                                                  "off", cbcSettings->getZeroHalfMode(),
                                                  CUTS_LONGHELP " This implementation was written by Alberto Caprara.");
  
  // Populate the keyword lists
  for (int code = CbcSolverParam::CBCSOLVER_FIRSTCUTPARAM + 1;
       code < CbcSolverParam::CBCSOLVER_LASTCUTPARAM; code++){
     // First the common keywords
     switch (code) {
      case CbcSolverParam::CUTSTRATEGY:
      case CbcSolverParam::CLIQUECUTS:
      case CbcSolverParam::FLOWCUTS:
      case CbcSolverParam::GMICUTS:
      case CbcSolverParam::GOMORYCUTS:
      case CbcSolverParam::KNAPSACKCUTS:
      case CbcSolverParam::LANDPCUTS:
      case CbcSolverParam::MIRCUTS:
      case CbcSolverParam::ODDWHEELCUTS:
      case CbcSolverParam::PROBINGCUTS:
      case CbcSolverParam::RESIDCAPCUTS:
      case CbcSolverParam::TWOMIRCUTS:
      case CbcSolverParam::ZEROHALFCUTS: {
         parameters[code]->appendKwd("on", CbcSolverParam::CGOn);
         parameters[code]->appendKwd("root", CbcSolverParam::CGRoot);
         parameters[code]->appendKwd("ifmove", CbcSolverParam::CGIfMove);
         parameters[code]->appendKwd("forceon", CbcSolverParam::CGForceOn);
         break;
      }
      default: {
         break;
     }
     
     // Now, add some additional keywords for different classes
     switch (code) {
      case CbcSolverParam::GOMORYCUTS: {
         parameters[code]->appendKwd("onglobal", CbcSolverParam::CGOnGlobal);
         parameters[code]->appendKwd("forceandglobal",
                                     CbcSolverParam::CGForceAndGlobal);
         parameters[code]->appendKwd("forcelongon", CbcSolverParam::CGForceLongOn);
         parameters[code]->appendKwd("longer", CbcSolverParam::CGLonger);
         parameters[code]->appendKwd("shorter", CbcSolverParam::CGShorter);
         break;
      }
      case CbcSolverParam::GMICUTS: {
         parameters[code]->appendKwd("long", CbcSolverParam::CGLong);
         parameters[code]->appendKwd("longroot", CbcSolverParam::CGLongRoot);
         parameters[code]->appendKwd("longifmove", CbcSolverParam::CGLongIfMove);
         parameters[code]->appendKwd("forcelongon", CbcSolverParam::CGForceLongOn);
         parameters[code]->appendKwd("longendonly", CbcSolverParam::CGLongEndOnly);
         break;
      }
      case CbcSolverParam::LAGOMORYCUTS: {
         parameters[code]->appendKwd("root", CbcSolverParam::CGRoot);
         parameters[code]->appendKwd("onlyaswellroot",
                                     CbcSolverParam::CGOnlyAsWellRoot);
         parameters[code]->appendKwd("cleanaswellroot",
                                     CbcSolverParam::CGCleanAsWellRoot);
         parameters[code]->appendKwd("bothaswellroot",
                                     CbcSolverParam::CGCleanBothAsWellRoot);
         // Here, we intentionally drop through to the next set
      }
      case CbcSolverParam::LATWOMIRCUTS: {
         parameters[code]->appendKwd("endonlyroot", CbcSolverParam::CGEndOnlyRoot);
         parameters[code]->appendKwd("endcleanroot", CbcSolverParam::CGEndCleanRoot);
         parameters[code]->appendKwd("endonly", CbcSolverParam::CGEndOnly);
         parameters[code]->appendKwd("endclean", CbcSolverParam::CGEndClean);
         parameters[code]->appendKwd("endboth", CbcSolverParam::CGEndBoth);
         parameters[code]->appendKwd("onlyaswell", CbcSolverParam::CGOnlyAsWell);
         parameters[code]->appendKwd("cleanaswell", CbcSolverParam::CGCleanAsWell);
         parameters[code]->appendKwd("bothaswell", CbcSolverParam::CGBothAsWell);
         parameters[code]->appendKwd("onlyinstead", CbcSolverParam::CGOnlyInstead);
         parameters[code]->appendKwd("cleaninstead",CbcSolverParam::CGCleanInstead);
         parameters[code]->appendKwd("bothinstead", CbcSolverParam::CGBothInstead);
         break;
      }
      case CbcSolverParam::ODDWHEELCUTS: {
         parameters[code]->appendKwd("onglobal", CbcSolverParam::CGOnGlobal);
         break;
      }
      case CbcSolverParam::PROBINGCUTS: {
         parameters[code]->appendKwd("forceonbut", CbcSolverParam::CGForceOnBut);
         break;
      }
      case CbcSolverParam::REDSPLITCUTS:
      case CbcSolverParam::REDSPLIT2CUTS: {
        parameters[code]->appendKwd("on", CbcSolverParam::CGOn);
        parameters[code]->appendKwd("root", CbcSolverParam::CGRoot);
        parameters[code]->appendKwd("longon", CbcSolverParam::CGLongOn);
        parameters[code]->appendKwd("longroot", CbcSolverParam::CGLongRoot);
        break;
      }
      default:
        break;
     }
     }
  }
}

//###########################################################################
//###########################################################################
   
void addCbcSolverHeurParams(CoinParamVec &parameters, CbcSolverSettings *cbcSettings)
{
  CbcSolverParam *param;
   
  for (int code = CbcSolverParam::CBCSOLVER_FIRSTHEURPARAM + 1;
       code < CbcSolverParam::CBCSOLVER_LASTHEURPARAM; code++){
     param = new CbcSolverParam;
     param->setObj(cbcSettings);
     param->setType(CoinParam::coinParamKwd);
     param->setPushFunc(pushCbcSolverKwdParam);
     parameters[code] = param;
  }

  parameters[CbcSolverParam::COMBINE]->setup("combine!Solutions",
                                             "Whether to use combine solution heuristic",
                                             "off", cbcSettings->getCombineMode(),
                                             "This switches on a heuristic which does branch and cut on the problem given by just using variables which have appeared in one or more solutions. It is obviously only tried after two or more solutions.");

  parameters[CbcSolverParam::CROSSOVER]->setup("combine2!Solutions",
                                                "Whether to use crossover solution heuristic",
                                                "off", cbcSettings->getCrossoverMode());

  parameters[CbcSolverParam::DINS]->setup("Dins",
                                          "Whether to try Distance Induced Neighborhood Search",
                                          "off", cbcSettings->getDinsMode(),
                                          HEURISTICS_LONGHELP);


  parameters[CbcSolverParam::DIVINGC]->setup("DivingC!oefficient",
                                             "Whether to try Coefficient diving heuristic",
                                             "off", cbcSettings->getDiveCoefficientMode(),
                                             HEURISTICS_LONGHELP);
  
  parameters[CbcSolverParam::DIVINGF]->setup("DivingF!ractional",
                                             "Whether to try Fractional diving heuristic",
                                             "off", cbcSettings->getDiveFractionalMode(),
                                             HEURISTICS_LONGHELP);

  parameters[CbcSolverParam::DIVINGG]->setup("DivingG!uided",
                                             "Whether to try Guided diving heuristic",
                                             "off", cbcSettings->getDiveGuidedMode(),
                                             HEURISTICS_LONGHELP);

  parameters[CbcSolverParam::DIVINGL]->setup("DivingL!ineSearch",
                                             "Whether to try Linesearch diving heuristic",
                                             "off", cbcSettings->getDiveLineSearchMode(),
                                             HEURISTICS_LONGHELP);

  parameters[CbcSolverParam::DIVINGP]->setup("DivingP!seudocost",
                                             "Whether to try Pseudocost diving heuristic",
                                             "off", cbcSettings->getDivePseudocostMode(),
                                             HEURISTICS_LONGHELP);
  
  parameters[CbcSolverParam::DIVINGS]->setup("DivingS!ome",
                                             "Whether to try Diving heuristics",
                                             "off", cbcSettings->getDiveRandomMode(),
                                             "This switches on a random diving heuristic at various times. One may prefer to individually turn diving heuristics on or off. ");

  parameters[CbcSolverParam::DIVINGV]->setup("DivingV!ectorLength",
                                             "Whether to try Vectorlength diving heuristic",
                                             "off", cbcSettings->getDiveVectorLengthMode(),
                                             HEURISTICS_LONGHELP);

  parameters[CbcSolverParam::DW]->setup("dw!Heuristic",
                                        "Whether to try Dantzig Wolfe heuristic",
                                        "off", cbcSettings->getDWMode(),
                                        "This heuristic is very very compute intensive. It tries to find a Dantzig Wolfe structure and use that. " HEURISTICS_LONGHELP);

  parameters[CbcSolverParam::FPUMP]->setup("feas!ibilityPump",
                                           "Whether to try Feasibility Pump",
                                           "off", cbcSettings->getFeasPumpMode(),
                                           "This switches on feasibility pump heuristic at root. This is due to Fischetti and Lodi and uses a sequence of LPs to try and get an integer feasible solution.  Some fine tuning is available by passFeasibilityPump.");

  parameters[CbcSolverParam::GREEDY]->setup("greedy!Heuristic",
                                            "Whether to use a greedy heuristic", "off",
                                            cbcSettings->getGreedyCoverMode(),
                                            "Switches on a pair of greedy heuristic which will try and obtain a solution.  It may just fix a percentage of variables and then try a small branch and cut run.");

  parameters[CbcSolverParam::HEURISTICSTRATEGY]->setup("heur!isticsOnOff",
                                                       "Switches most heuristics on or off",
                                                       "off", 0,
                                                       "This can be used to switch on or off all heuristics.  Then you can set individual ones off or on.  CbcTreeLocal is not included as it dramatically alters search.");

  parameters[CbcSolverParam::NAIVE]->setup("naive!Heuristics",
                                           "Whether to try some stupid heuristic",
                                           "off", cbcSettings->getNaiveHeurMode(),
                                           "This is naive heuristics which, e.g., fix all integers with costs to zero!. " HEURISTICS_LONGHELP,
                                           CoinParam::displayPriorityLow
);

  parameters[CbcSolverParam::PIVOTANDFIX]->setup("pivotAndF!ix",
                                                 "Whether to try Pivot and Fix heuristic",
                                                 "off", cbcSettings->getPivotAndFixMode(),
                                                 HEURISTICS_LONGHELP);

#if 0
  parameters[CbcSolverParam::PIVOTANDCOMPLEMENT]->setup("pivotAndC!omplement",
                                                        "Whether to try Pivot and Complement heuristic",
                                                        "off", cbcSettings->getPivotAndComplementMode(),
                                                        HEURISTICS_LONGHELP);
#endif

  parameters[CbcSolverParam::PROXIMITY]->setup("proximity!Search",
                                               "Whether to do proximity search heuristic",
                                               "off", cbcSettings->getProximityMode(),
                                               "This heuristic looks for a solution close to the incumbent solution (Fischetti and Monaci, 2012). The idea is to define a sub-MIP without additional constraints but with a modified objective function intended to attract the search in the proximity of the incumbent. The approach works well for 0-1 MIPs whose solution landscape is not too irregular (meaning the there is reasonable probability of finding an improved solution by flipping a small number of binary variables), in particular when it is applied to the first heuristic solutions found at the root node. "
      HEURISTICS_LONGHELP); // Can also set different maxNode cbcSettings by plusnnnn (and are 'on'(on==30)).

  parameters[CbcSolverParam::RANDROUND]->setup("randomi!zedRounding",
                                               "Whether to try randomized rounding heuristic",
                                               "off", cbcSettings->getRandRoundMode(),
                                               HEURISTICS_LONGHELP);

  parameters[CbcSolverParam::RENS]->setup("Rens",
                                          "Whether to try Relaxation Enforced Neighborhood Search",
                                          "off", cbcSettings->getRensMode(),
                                          HEURISTICS_LONGHELP " Value 'on' just does 50 nodes. 200, 1000, and 10000 does that many nodes.");

  parameters[CbcSolverParam::RINS]->setup("Rins",
                                          "Whether to try Relaxed Induced Neighborhood Search",
                                          "off", cbcSettings->getRinsMode(),
                                          HEURISTICS_LONGHELP);

  parameters[CbcSolverParam::ROUNDING]->setup("round!ingHeuristic",
                                              "Whether to use Rounding heuristic",
                                              "off", cbcSettings->getRoundingMode(),
                                              "This switches on a simple (but effective) rounding heuristic at each node of tree.");

  parameters[CbcSolverParam::LOCALTREE]->setup("local!TreeSearch",
                                               "Whether to use local tree search", "off",
                                               cbcSettings->getLocalTreeMode(), "This switches on a local search algorithm when a solution is found.  This is from Fischetti and Lodi and is not really a heuristic although it can be used as one. When used from this program it has limited functionality.  It is not controlled by heuristicsOnOff.");

  parameters[CbcSolverParam::VND]->setup("Vnd!VariableNeighborhoodSearch",
                                         "Whether to try Variable Neighborhood Search",
                                         "off", cbcSettings->getVndMode(),
                                         HEURISTICS_LONGHELP);
  
  // Populate the keyword lists
  for (int code = CbcSolverParam::CBCSOLVER_FIRSTHEURPARAM + 1;
       code < CbcSolverParam::CBCSOLVER_LASTHEURPARAM; code++){
     // First the common keywords
     switch (code) {
      case CbcSolverParam::HEURISTICSTRATEGY:
      case CbcSolverParam::CROSSOVER: 
      case CbcSolverParam::DINS: 
      case CbcSolverParam::DIVINGC:
      case CbcSolverParam::DIVINGF:
      case CbcSolverParam::DIVINGG:
      case CbcSolverParam::DIVINGL:
      case CbcSolverParam::DIVINGP:
      case CbcSolverParam::DIVINGS:
      case CbcSolverParam::DIVINGV: 
      case CbcSolverParam::DW: 
      case CbcSolverParam::NAIVE: 
      case CbcSolverParam::PIVOTANDFIX: 
      case CbcSolverParam::PIVOTANDCOMPLEMENT:
      case CbcSolverParam::PROXIMITY: 
      case CbcSolverParam::RANDROUND: 
      case CbcSolverParam::RENS: 
      case CbcSolverParam::RINS: 
      case CbcSolverParam::ROUNDING: 
      case CbcSolverParam::VND: {
        parameters[code]->appendKwd("off", CbcSolverParam::CGOff);
        parameters[code]->appendKwd("on", CbcSolverParam::CGOn);
        parameters[code]->appendKwd("both", CbcSolverParam::CGRoot);
        parameters[code]->appendKwd("before", CbcSolverParam::CGIfMove);
        break;
      }
      default: 
         break;
     }
     // Check the unique keywords
     switch (code) {
      case CbcSolverParam::COMBINE: {
        parameters[code]->appendKwd("off", CbcSolverParam::HeurOff);
        parameters[code]->appendKwd("on", CbcSolverParam::HeurOn);
        break;
      }
      case CbcSolverParam::DINS: {
        parameters[code]->appendKwd("often", CbcSolverParam::HeurOften);
        break;
      }
      case CbcSolverParam::FPUMP: {
        parameters[code]->appendKwd("off", CbcSolverParam::HeurOff);
        parameters[code]->appendKwd("on", CbcSolverParam::HeurOn);
        break;
      }
      case CbcSolverParam::GREEDY: {
        parameters[code]->appendKwd("off", CbcSolverParam::HeurOff);
        parameters[code]->appendKwd("on", CbcSolverParam::HeurOn);
        parameters[code]->appendKwd("root", CbcSolverParam::HeurRoot);
        break;
      }
      case CbcSolverParam::LOCALTREE: {
        parameters[code]->appendKwd("off", CbcSolverParam::HeurOff);
        parameters[code]->appendKwd("on", CbcSolverParam::HeurOn);
      }
      case CbcSolverParam::PROXIMITY: {
        parameters[code]->appendKwd("10", CbcSolverParam::HeurTen);
        parameters[code]->appendKwd("100", CbcSolverParam::HeurOneHundred);
        parameters[code]->appendKwd("300", CbcSolverParam::HeurThreeHundred);
        break;
      }
      case CbcSolverParam::RENS: { 
        parameters[code]->appendKwd("200", CbcSolverParam::HeurTwoHundred);
        parameters[code]->appendKwd("1000", CbcSolverParam::HeurOneThousand);
        parameters[code]->appendKwd("10000", CbcSolverParam::HeurTenThousand);
        parameters[code]->appendKwd("dj", CbcSolverParam::HeurDj);
        parameters[code]->appendKwd("djbefore", CbcSolverParam::HeurDjBefore);
        parameters[code]->appendKwd("usesolution", CbcSolverParam::HeurUseSolution);
        break;
      }
      case CbcSolverParam::RINS: {
        parameters[code]->appendKwd("often", CbcSolverParam::HeurOften);
        break;
      }
      case CbcSolverParam::VND: { 
        parameters[code]->appendKwd("intree", CbcSolverParam::HeurInTree);
        break;
      }
      default:
        break;
     }
  }  
}

//###########################################################################
//###########################################################################
   
void loadGenParamObj(const CoinParamVec paramVec, int first, int last,
  CbcSolverSettings *cbcSettings)

{
  int i;
  /*
      Load the cbc-generic object into the parameters
    */
  for (i = first; i <= last; i++) {
    CbcSolverParam *genParam = dynamic_cast< CbcSolverParam * >(paramVec[i]);
    assert(genParam != 0);
    genParam->setObj(cbcSettings);
  }

  return;
}

//###########################################################################
//###########################################################################
   
/* Functions to implement  cbc-generic (CbcSolverParam) parameters */

/*
  Maintainer's utility, scan the parameters and report the ones that are
  unimplemented (i.e., have no pushFunc).
*/

int doUnimplementedParam(CoinParam *param)

{
  assert(param != 0);

  CbcSolverParam *genParam = dynamic_cast< CbcSolverParam * >(param);
  assert(genParam != 0);

  CbcSolverSettings *cbcSettings = genParam->obj();
  assert(cbcSettings != 0);

  assert(cbcSettings->paramVec_ != 0);
  CoinParamVec &paramVec = *cbcSettings->paramVec_;

  int unimpCnt = 0;
  int maxAcross = 5;
  for (unsigned i = 0; i < paramVec.size(); i++) {
    CoinParam *param = paramVec[i];
    if (param->pushFunc() == 0) {
      if (unimpCnt % maxAcross == 0) {
        std::cout << std::endl;
      } else {
        std::cout << " ";
      }
      std::cout << param->name();
      unimpCnt++;
    }
  }
  if (unimpCnt % maxAcross != 1) {
    std::cout << std::endl;
  }
  std::cout << unimpCnt << " unimplemented parameters." << std::endl;

  return (0);
}

//###########################################################################
//###########################################################################
   
/*
  Noop function. Mainly to eliminate commands from the list returned by
  doUnimplmentedParam.
*/

int doNothingParam(CoinParam *param)
{
  return (0);
}

/*
  Function to terminate command parsing by returning -1.
*/

int doExitParam(CoinParam *param)

{
  return (-1);
}

//###########################################################################
//###########################################################################
   
/*
  Function to print the current version.
*/

int doVersionParam(CoinParam *param)

{
  assert(param != 0);
  CbcSolverParam *genParam = dynamic_cast< CbcSolverParam * >(param);
  assert(genParam != 0);
  CbcSolverSettings *cbcSettings = genParam->obj();
  assert(cbcSettings != 0);

  std::cout << "cbc-generic version " << cbcSettings->version_ << std::endl;
  std::cout
    << "cbc-generic is experimental software. If you want a stable MIP "
    << "solver, please" << std::endl
    << "use cbc. If you discover bugs while using cbc-generic "
    << "please specify" << std::endl
    << "cbc-generic in the ticket description or email subject line."
    << std::endl;

  return (0);
}

//###########################################################################
//###########################################################################
   
/*
  Function to handle help (HELP), `?' (GENERALQUERY), and `???'
  (FULLGENERALQUERY).
*/

int doHelpParam(CoinParam *param)

{
  assert(param != 0);
  CbcSolverParam *genParam = dynamic_cast< CbcSolverParam * >(param);
  assert(genParam != 0);
  CbcSolverSettings *cbcSettings = genParam->obj();
  assert(cbcSettings != 0);

  CbcSolverParam::CbcSolverParamCode code = genParam->paramCode();

  int verbose = cbcSettings->verbose_;
  bool shortHelp = ((verbose & 0x01) ? true : false);
  bool longHelp = ((verbose & 0x02) ? true : false);
  bool hidden = ((verbose & 0x08) ? true : false);

  CoinParamVec *paramVec = cbcSettings->paramVec_;
  assert(paramVec != 0);
  /*
      Tune up the initial cbcSettings. FULLGENERALQUERY will print normally hidden
      params, and a request for long help overrules a request for short help.
    */
  if (code == CbcSolverParam::FULLGENERALQUERY) {
    hidden = true;
  }
  if (longHelp) {
    shortHelp = false;
  }

  CoinParamUtils::printGenericHelp();

  std::cout << "\nAvailable commands are:";
  std::string pfx("  ");
  CoinParamUtils::printHelp(*paramVec, 0, paramVec->size() - 1, pfx,
    shortHelp, longHelp, hidden);

  return (0);
}

//###########################################################################
//###########################################################################
   
/*
  Function to push a double-valued parameter.
*/

int pushCbcSolverDblParam(CoinParam *param)

{
  assert(param != 0);
  CbcSolverParam *genParam = dynamic_cast< CbcSolverParam * >(param);
  assert(genParam != 0);
  CbcSolverSettings *cbcSettings = genParam->obj();
  assert(cbcSettings != 0);

  double val = genParam->dblVal();
  CbcSolverParam::CbcSolverParamCode code = genParam->paramCode();

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
    */
  switch (code) {
   case CbcSolverParam::ARTIFICIALCOST: {
      cbcSettings->setArtVarMode(CbcSolverParam::ParamOn, val);
      break;
   }
   case CbcSolverParam::DEXTRA3: {
      cbcSettings->setExtraDbl3(val);
      break;
   }
   case CbcSolverParam::DEXTRA4: {
      cbcSettings->setExtraDbl4(val);
      break;
   }
   case CbcSolverParam::DEXTRA5: {
      cbcSettings->setExtraDbl5(val);
      break;
   }
   case CbcSolverParam::DJFIX: {
      cbcSettings->setDjFixMode(CbcSolverParam::ParamOn, val);
      break;
   }
   case CbcSolverParam::FAKECUTOFF: {
      cbcSettings->setFeasPumpFakeCutoff(val);
      break;
   }
   case CbcSolverParam::FAKEINCREMENT: {
      cbcSettings->setFeasPumpFakeIncrement(val);
      break;
   }
   case CbcSolverParam::SMALLBAB: {
      cbcSettings->setSmallBab(val);
      break;
   }
   case CbcSolverParam::TIGHTENFACTOR: {
      cbcSettings->setTightenFactor(val);
      break;
   }

   default: {
      std::cerr << "pushCbcSolverDbl: no equivalent CbcSolverSettings field for "
              << "parameter code `" << code << "'." << std::endl;
    retval = -1;
    break;
  }
  }

  return (retval);
}

//###########################################################################
//###########################################################################
   
/*
  Function to push an integer-valued parameter.
*/

int pushCbcSolverIntParam(CoinParam *param)

{
  assert(param != 0);
  CbcSolverParam *genParam = dynamic_cast< CbcSolverParam * >(param);
  assert(genParam != 0);
  CbcSolverSettings *cbcSettings = genParam->obj();
  assert(cbcSettings != 0);

  int val = genParam->intVal();
  CbcSolverParam::CbcSolverParamCode code = genParam->paramCode();

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
    */
  switch (code) {
  case CbcSolverParam::BKPIVOTINGSTRATEGY: {
    cbcSettings->setBkPivotStrategy(val);
    break;
  }
  case CbcSolverParam::BKMAXCALLS: {
    cbcSettings->setBkMaxCalls(val);
     break;
  }
  case CbcSolverParam::BKCLQEXTMETHOD: {
    cbcSettings->setBkClqExtMethod(val);
     break;
  }
  case CbcSolverParam::CPP: {
    cbcSettings->setCppMode(val);
     break;
  }
  case CbcSolverParam::CUTDEPTH: {
    cbcSettings->setCutDepth(val);
    break;
  }
  case CbcSolverParam::CUTLENGTH: {
    cbcSettings->setCutLength(val);
     break;
  }
  case CbcSolverParam::CUTPASSINTREE: {
    cbcSettings->setCutPassInTree(val);
     break;
  }
  case CbcSolverParam::DEPTHMINIBAB: {
    cbcSettings->setDepthMiniBaB(val);
     break;
  }
  case CbcSolverParam::DIVEOPT: {
    cbcSettings->setDiveOpt(val);
     break;
  }
  case CbcSolverParam::DIVEOPTSOLVES: {
    cbcSettings->setDiveOptSolves(val);
     break;
  }
  case CbcSolverParam::EXPERIMENT: {
    cbcSettings->setExperimentMode(val);
     break;
  }
  case CbcSolverParam::EXTRA1: {
    cbcSettings->setExtraIntParam1(val);
     break;
  }
  case CbcSolverParam::EXTRA2: {
    cbcSettings->setExtraIntParam2(val);
     break;
  }
  case CbcSolverParam::EXTRA3: {
    cbcSettings->setExtraIntParam3(val);
     break;
  }
  case CbcSolverParam::EXTRA4: {
    cbcSettings->setExtraIntParam4(val);
     break;
  }
  case CbcSolverParam::FPUMPITS: {
    cbcSettings->setFeasPumpIters(val);
     break;
  }
  case CbcSolverParam::FPUMPTUNE: {
    cbcSettings->setFeasPumpTune(val);
     break;
  }
  case CbcSolverParam::FPUMPTUNE2: {
    cbcSettings->setFeasPumpTune2(val);
     break;
  }
  case CbcSolverParam::HEUROPTIONS: {
    cbcSettings->setHeurOptions(val);
     break;
  }
  case CbcSolverParam::LOGLEVEL: {
    cbcSettings->setLogLevel(val);
    break;
  }
  case CbcSolverParam::LPLOGLEVEL: {
    cbcSettings->setLpLogLevel(val);
     break;
  }
  case CbcSolverParam::MAXSAVEDSOLS: {
     cbcSettings->setMaxSavedSols(val);
    break;
  }
  case CbcSolverParam::MAXSLOWCUTS: {
    cbcSettings->setMaxSlowCuts(val);
     break;
  }
  case CbcSolverParam::MOREMOREMIPOPTIONS: {
    cbcSettings->setMoreMoreOptions(val);
     break;
  }
  case CbcSolverParam::MULTIPLEROOTS: {
    cbcSettings->setMultipleRoots(val);
     break;
  }
  case CbcSolverParam::ODDWEXTMETHOD: {
    cbcSettings->setOddWextMethod(val);
     break;
  }
  case CbcSolverParam::OUTPUTFORMAT: {
    cbcSettings->setOutputFormat(val);
     break;
  }
  case CbcSolverParam::PRINTOPTIONS: {
     cbcSettings->setPrintOptions(val);
    break;
  }
  case CbcSolverParam::PROCESSTUNE: {
    cbcSettings->setProcessTune(val);
     break;
  }
  case CbcSolverParam::RANDOMSEED: {
    cbcSettings->setRandomSeed(val);
     break;
  }
  case CbcSolverParam::STRONGSTRATEGY: {
    cbcSettings->setStrongStrategy(val);
     break;
  }
  case CbcSolverParam::TESTOSI: {
    cbcSettings->setTestOsi(val);
     break;
  }
#ifdef CBC_THREADS
  case CbcSolverParam::THREADS: {
    cbcSettings->setThreads(val);
     break;
  }
#endif
  case CbcSolverParam::USERCBC: {
    cbcSettings->setUserCbc(val);
     break;
  }
  case CbcSolverParam::VERBOSE: {
     cbcSettings->setVerbose(val);
    break;
  }
  case CbcSolverParam::VUBTRY: {
    cbcSettings->setVubTry(val);
     break;
  }
  default: {
    std::cerr << "pushCbcSolverInt: no method for storing "
              << "parameter code `" << code << "'." << std::endl;
    retval = -1;
    break;
  }
  }

  return (retval);
}

//###########################################################################
//###########################################################################
   
/*
  Function to push a keyword-valued parameter. This is the catch-all function
  for keyword parameters that don't belong to any other useful grouping.
*/

int pushCbcSolverKwdParam(CoinParam *param)
{
  assert(param != 0);
  CbcSolverParam *genParam = dynamic_cast< CbcSolverParam * >(param);
  assert(genParam != 0);
  CbcSolverSettings *cbcSettings = genParam->obj();
  assert(cbcSettings != 0);

  int mode = genParam->modeVal();
  CbcSolverParam::CbcSolverParamCode code = genParam->paramCode();

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
  */

  switch (code) {
   case CbcSolverParam::ALLCOMMANDS: {
      cbcSettings->setCommandMode(static_cast<CbcSolverParam::CommandMode>(mode));
      break;
   }
   case CbcSolverParam::CLQSTRENGTHENING: {
      cbcSettings->setClqStrMode(static_cast<CbcSolverParam::ClqStrMode>(mode));
      break;
   }
   case CbcSolverParam::BRANCHPRIORITY: {
      cbcSettings->priorityMode_ = CbcSolverParam::BPOff;
      break;
   }
   case CbcSolverParam::CUTOFFCONSTRAINT: {
      cbcSettings->setCutoffMode(static_cast<CbcSolverParam::CutoffMode>(mode));
      break;
   }
   case CbcSolverParam::INTPRINT: {
      cbcSettings->setIntPrintMode(static_cast<CbcSolverParam::IntPrintMode>(mode));
      break;
   }
   case CbcSolverParam::NODESTRATEGY: {
      cbcSettings->setNodeStrategy(static_cast<CbcSolverParam::NodeStrategy>(mode));
      break;
   }
   case CbcSolverParam::ORBITAL: {
      cbcSettings->setOrbitalStrategy(static_cast<CbcSolverParam::OrbitalStrategy>(mode));
      break;
   }
   case CbcSolverParam::PREPROCESS: {
      cbcSettings->setIPPMode(static_cast<CbcSolverParam::IPPMode>(mode));
      break;
   }
   case CbcSolverParam::SOSPRIORITIZE: {
      cbcSettings->setSOSStrategy(static_cast<CbcSolverParam::SOSStrategy>(mode));
      break;
   }
   case CbcSolverParam::STRATEGY: {
      cbcSettings->setStrategyMode(static_cast<CbcSolverParam::StrategyMode>(mode));
      break;
   }
   case CbcSolverParam::TIMEMODE: {
      cbcSettings->setClockType(static_cast<CbcSolverParam::ClockType>(mode));
      break;
   }
   case CbcSolverParam::USECGRAPH: {
      cbcSettings->setCGraphMode(static_cast<CbcSolverParam::CGraphMode>(mode));
      break;
   }
   default: 
      break;
  }

  return (retval);
}

//###########################################################################
//###########################################################################
   
/*
  Function to push a bool-valued parameter. These are really just keyword
  parameters that take values "on" and "off"

*/

int pushCbcSolverBoolParam(CoinParam *param)
{
  assert(param != 0);
  CbcSolverParam *genParam = dynamic_cast< CbcSolverParam * >(param);
  assert(genParam != 0);
  CbcSolverSettings *cbcSettings = genParam->obj();
  assert(cbcSettings != 0);

  // This is ugly, get keyword and set parameter with it instead.
  CbcSolverParam::OnOffMode mode = static_cast<CbcSolverParam::OnOffMode>(genParam->modeVal());
  CbcSolverParam::CbcSolverParamCode code = genParam->paramCode();

  int retval = 0;

  switch (code) {
   case CbcSolverParam::CPX: {
      cbcSettings->setCPXMode(mode);
      break;
   }
   case CbcSolverParam::DOHEURISTIC: {
      cbcSettings->setDoHeuristicMode(mode);
      break;
   }
   case CbcSolverParam::ERRORSALLOWED: {
      cbcSettings->setImportErrorsMode(mode);
      break;
   }
   case CbcSolverParam::MESSAGES: {
      cbcSettings->setMessagePrefixMode(mode);
      break;
   }
   case CbcSolverParam::PREPROCNAMES: {
      cbcSettings->setPreProcNamesMode(mode);
      break;
   }
   case CbcSolverParam::SOS: {
      cbcSettings->setSOSMode(mode);
      break;
   }
   case CbcSolverParam::USESOLUTION: {
      cbcSettings->setUseSolutionMode(mode);
      break;
   }
   default:
     break;
  }
  
  return (retval);
}

//###########################################################################
//###########################################################################
   
/*
  Function to push a keyword-valued parameter related to heuristics. 
*/

int pushCbcSolverHeurParam(CoinParam *param)
{
  assert(param != 0);
  CbcSolverParam *genParam = dynamic_cast< CbcSolverParam * >(param);
  assert(genParam != 0);
  CbcSolverSettings *cbcSettings = genParam->obj();
  assert(cbcSettings != 0);

  // This is ugly, get keyword and set parameter with it instead.
  CbcSolverParam::HeurMode mode = static_cast<CbcSolverParam::HeurMode>(genParam->modeVal());
  CbcSolverParam::CbcSolverParamCode code = genParam->paramCode();

  int retval = 0;

  /*
      We've done the basic checks; go ahead and set the relevant field in the
      control block. We shouldn't need the default case, but some compilers will
      complain if it's missing.
  */
  switch (code) {
   case CbcSolverParam::HEURISTICSTRATEGY: {
      cbcSettings->setGreedyCoverMode(mode);
      cbcSettings->setGreedyEqualityMode(mode);
      cbcSettings->setCrossoverMode(mode);
      cbcSettings->setDinsMode(mode);
      cbcSettings->setDiveCoefficientMode(mode);
      cbcSettings->setDiveFractionalMode(mode);
      cbcSettings->setDiveGuidedMode(mode);
      cbcSettings->setDiveLineSearchMode(mode);
      cbcSettings->setDivePseudocostMode(mode);
      cbcSettings->setDiveVectorLengthMode(mode);
      cbcSettings->setDWMode(mode);
      cbcSettings->setNaiveHeurMode(mode);
      cbcSettings->setPivotAndFixMode(mode);
#if 0
      cbcSettings->setPivotAndComplementMode(mode);
#endif
      cbcSettings->setProximityMode(mode);
      cbcSettings->setRandRoundMode(mode);
      cbcSettings->setRensMode(mode);
      cbcSettings->setRinsMode(mode);
      cbcSettings->setRoundingMode(mode);
      cbcSettings->setVndMode(mode);
      break;
   }
   case CbcSolverParam::COMBINE:  {
      cbcSettings->setCombineMode(mode);
      break;
   }
   case CbcSolverParam::CROSSOVER:  {
      cbcSettings->setCrossoverMode(mode);
      break;
   }
   case CbcSolverParam::DINS:  {
      cbcSettings->setDinsMode(mode);
      break;
   }
   case CbcSolverParam::DIVINGC: {
      cbcSettings->setDiveCoefficientMode(mode);
      break;
   }
   case CbcSolverParam::DIVINGF: {
      cbcSettings->setDiveFractionalMode(mode);
      break;
   }
   case CbcSolverParam::DIVINGG: {
      cbcSettings->setDiveGuidedMode(mode);
      break;
   }
   case CbcSolverParam::DIVINGL: {
      cbcSettings->setDiveLineSearchMode(mode);
      break;
   }
   case CbcSolverParam::DIVINGP: {
      cbcSettings->setDivePseudocostMode(mode);
      break;
   }
   case CbcSolverParam::DIVINGS: {
      cbcSettings->setDiveRandomMode(mode);
      break;
   }
   case CbcSolverParam::DIVINGV:  {
      cbcSettings->setDiveVectorLengthMode(mode);
      break;
   }
   case CbcSolverParam::DW:  {
      cbcSettings->setDWMode(mode);
      break;
   }
   case CbcSolverParam::FPUMP: {
      cbcSettings->setFeasPumpMode(mode);
      break;
   }
   case CbcSolverParam::GREEDY: {
      cbcSettings->setGreedyCoverMode(mode);
      cbcSettings->setGreedyEqualityMode(mode);
      break;
   }
   case CbcSolverParam::NAIVE: {
      cbcSettings->setNaiveHeurMode(mode);
      break;
   }
   case CbcSolverParam::PIVOTANDFIX: {
      cbcSettings->setPivotAndFixMode(mode);
      break;
   }
#if 0
   case CbcSolverParam::PIVOTANDCOMPLEMENT: {
      cbcSettings->setPivotAndComplementMode(mode);
      break;
   }
#endif
   case CbcSolverParam::PROXIMITY: {
      cbcSettings->setProximityMode(mode);
      break;
   }
   case CbcSolverParam::RANDROUND: {
      cbcSettings->setRandRoundMode(mode);
      break;
   }
   case CbcSolverParam::RENS: {
      cbcSettings->setRensMode(mode);
      break;
   }
   case CbcSolverParam::RINS: {
      cbcSettings->setRinsMode(mode);
      break;
   }
   case CbcSolverParam::ROUNDING: {
      cbcSettings->setRoundingMode(mode);
      break;
   }
   case CbcSolverParam::VND: {
      cbcSettings->setVndMode(mode);
      break;
   }
   default:
     break;
  }

  return (retval);
}

//###########################################################################
//###########################################################################
   
/*
  Function to push a string-valued parameter
*/

int pushCbcSolverStrParam(CoinParam *param)

{
  assert(param != 0);
  CbcSolverParam *genParam = dynamic_cast< CbcSolverParam * >(param);
  assert(genParam != 0);
  CbcSolverSettings *cbcSettings = genParam->obj();
  assert(cbcSettings != 0);

  std::string str = genParam->strVal();
  CbcSolverParam::CbcSolverParamCode code = genParam->paramCode();

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
    */
  switch (code) {
  case CbcSolverParam::DIRECTORY: {
    char dirSep = CoinFindDirSeparator();
    if (str[str.length() - 1] != dirSep) {
      str += dirSep;
    }
    cbcSettings->dfltDirectory_ = str;
    break;
  }
  default: {
    std::cerr << "pushCbcSolverStr: no equivalent CbcSolverSettings field for "
              << "parameter code `" << code << "'." << std::endl;
    retval = -1;
    break;
  }
  }

  return (retval);
}

//###########################################################################
//###########################################################################
   
/*
  The various parameters to control cut generators can be
  grouped, as they all use the same set of keywords.
*/
int pushCbcSolverCutParam(CoinParam *param)

{
  assert(param != 0);
  CbcSolverParam *genParam = dynamic_cast< CbcSolverParam * >(param);
  assert(genParam != 0);
  CbcSolverSettings *cbcSettings = genParam->obj();
  assert(cbcSettings != 0);

  CbcSolverParam::CGMode mode = static_cast<CbcSolverParam::CGMode>(genParam->modeVal());
  CbcSolverParam::CbcSolverParamCode code = genParam->paramCode();
  /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we need
      to do is correct to 0 (no error) if we're successful.
    */
  int retval;

  if (CoinParamUtils::isInteractive()) {
    retval = 1;
  } else {
    retval = -1;
  }

  /*
      We've done the basic checks; go ahead and set the relevant field in the
      control block. We shouldn't need the default case, but some compilers will
      complain if it's missing.
  */
  switch (code) {
   case CbcSolverParam::CLIQUECUTS: {
      cbcSettings->setCliqueMode(mode);
      break;
   }
   case CbcSolverParam::FLOWCUTS: {
      cbcSettings->setFlowMode(mode);
      break;
   }
   case CbcSolverParam::GMICUTS: {
      cbcSettings->setGMIMode(mode);
      break;
   }
   case CbcSolverParam::GOMORYCUTS: {
      cbcSettings->setGomoryMode(mode);
      break;
   }
   case CbcSolverParam::KNAPSACKCUTS: {
      cbcSettings->setKnapsackMode(mode);
      break;
   }
   case CbcSolverParam::LAGOMORYCUTS: {
      cbcSettings->setLaGomoryMode(mode);
      break;
   }
   case CbcSolverParam::LANDPCUTS: {
      cbcSettings->setLandPMode(mode);
      break;
   }
   case CbcSolverParam::LATWOMIRCUTS: {
      cbcSettings->setLaTwomirMode(mode);
      break;
   }
   case CbcSolverParam::MIRCUTS: {
      cbcSettings->setMirMode(mode);
      break;
   }
   case CbcSolverParam::ODDWHEELCUTS: {
      cbcSettings->setOddWheelMode(mode);
      break;
   }
   case CbcSolverParam::PROBINGCUTS: {
      cbcSettings->setProbingMode(mode);
      break;
   }
   case CbcSolverParam::REDSPLITCUTS: {
      cbcSettings->setRedSplitMode(mode);
      break;
   }
   case CbcSolverParam::REDSPLIT2CUTS: {
      cbcSettings->setRedSplit2Mode(mode);
      break;
   }
   case CbcSolverParam::RESIDCAPCUTS: {
      cbcSettings->setResidCapMode(mode);
      break;
   }
   case CbcSolverParam::TWOMIRCUTS: {
      cbcSettings->setTwomirMode(mode);
      break;
   }
   case CbcSolverParam::ZEROHALFCUTS: {
      cbcSettings->setZeroHalfMode(mode);
      break;
   }
   case CbcSolverParam::CUTSTRATEGY: {
      cbcSettings->setCliqueMode(mode);
      cbcSettings->setFlowMode(mode);
      cbcSettings->setGMIMode(mode);
      cbcSettings->setGomoryMode(mode);
      cbcSettings->setKnapsackMode(mode);
      cbcSettings->setLaGomoryMode(mode);
      cbcSettings->setLandPMode(mode);
      cbcSettings->setLaTwomirMode(mode);
      cbcSettings->setMirMode(mode);
      cbcSettings->setOddWheelMode(mode);
      cbcSettings->setProbingMode(mode);
      cbcSettings->setRedSplitMode(mode);
      cbcSettings->setRedSplit2Mode(mode);
      cbcSettings->setRedSplit2Mode(mode);
      cbcSettings->setTwomirMode(mode);
      cbcSettings->setZeroHalfMode(mode);
      break;
   }
   default:
     break;
  }

  return (0);
}

//###########################################################################
//###########################################################################
   
/*
  This routine imports a new constraint system into the solver.
*/

int doImportParam(CoinParam *param)

{
  assert(param != 0);
  CbcSolverParam *genParam = dynamic_cast< CbcSolverParam * >(param);
  assert(genParam != 0);
  CbcSolverSettings *cbcSettings = genParam->obj();
  assert(cbcSettings != 0);
  /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we need
      to do is correct to 0 (no error) if we're successful.
    */
  int retval;
  if (CoinParamUtils::isInteractive()) {
    retval = 1;
  } else {
    retval = -1;
  }
  /*
      Figure out where we're going to acquire this new model. As special cases,
      `$' says `use the previous input source' and `-' says `use stdin'.
    */
  std::string field = genParam->strVal();
  std::string fileName;
  if (field == "$") {
    fileName = cbcSettings->lastMpsIn_;
    field = fileName;
  } else if (field == "-") {
    fileName = "stdin";
    field = fileName;
  } else {
    fileName = field;
  }
  /*
      See if we can open a file. fileCoinReadable understands a fair bit about
      platforms and compressed files and will try variations of the file name
      (see the doxygen doc'n for details). The file name returned in field wil
      be the one that actually worked.
    */
  bool canOpen = fileCoinReadable(fileName, cbcSettings->dfltDirectory_);
  if (canOpen == false) {
    std::cout
      << "Unable to open file `" << fileName
      << "', original name '" << genParam->strVal() << "'." << std::endl;
    return (retval);
  }
  /*
      We can find the file. Record the name. This requires a little finesse: what
      we want is the base file name (and extension(s), if present) but not the
      prefix, unless it's an absolute path.
    */
  if (!fileAbsPath(fileName)) {
    std::string::size_type pos = fileName.rfind(field);
    cbcSettings->lastMpsIn_ = fileName.substr(pos);
  } else {
    cbcSettings->lastMpsIn_ = fileName;
  }
  /*
      Try to read the file. Standard OSI doesn't support the Clp extensions for
      keepImportNames and allowImportErrors. It should at least support
      keepImportNames. Status will be zero for a successful read.
    */
  OsiSolverInterface *lpSolver = cbcSettings->model_->solver();
  int status = lpSolver->readMps(fileName.c_str(), "");
  if (status) {
    std::cout
      << "There were " << status << " errors on input." << std::endl;
    return (retval);
  }
  /*
      We have a model! Return success.
    */
  cbcSettings->goodModel_ = true;

  return (0);
}

//###########################################################################
//###########################################################################
   
/*
  This routine imports a debug file into the solver, or arranges for its
  creation. Import works in the standard way, using the file name provided
  with the command.

  As special cases, if the file name is `create' or `createAfterPre', the
  action here sets up to cause a debug file containing the solution to be
  dumped after branch-and-cut is completed.  `createAfterPre' will dump the
  solution before undoing the presolve transforms.  `create' will dump the
  solution after integer presolve is backed out.
*/

int doDebugParam(CoinParam *param)

{
  assert(param != 0);
  CbcSolverParam *genParam = dynamic_cast< CbcSolverParam * >(param);
  assert(genParam != 0);
  CbcSolverSettings *cbcSettings = genParam->obj();
  assert(cbcSettings != 0);
  /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we need
      to do is correct to 0 (no error) if we're successful.
    */
  int retval;
  if (CoinParamUtils::isInteractive()) {
    retval = 1;
  } else {
    retval = -1;
  }
  /*
       If the file name is `create' or `createAfterPre', we're just setting up to
       make a debug file the next time we do branch-and-cut.
    */
  std::string field = genParam->strVal();
  if (field == "create" || field == "createAfterPre") {
    cbcSettings->debugCreate_ = field;
    return (0);
  }
  /*
      Figure out where we're going to acquire the debug file. As special cases,
      `$' says `use the previous input source' and `-' says `use stdin'.
    */
  std::string fileName;
  if (field == "$") {
    fileName = cbcSettings->debugFile_;
    field = fileName;
  } else if (field == "-") {
    fileName = "stdin";
    field = fileName;
  } else {
    fileName = field;
  }
  /*
      See if we can open a file. fileCoinReadable understands a fair bit about
      platforms and compressed files and will try variations of the file name
      (see the doxygen doc'n for details). The file name returned in field wil be
      the one that actually worked. No default prefix --- a debug file is assumed
      to always be in the current directory.
    */
  bool canOpen = fileCoinReadable(fileName, cbcSettings->dfltDirectory_);
  if (canOpen == false) {
    std::cout
      << "Unable to open file `" << fileName
      << "', original name '" << genParam->strVal() << "'." << std::endl;
    return (retval);
  }
  /*
      We can find the file. Record the name. This requires a little finesse: what
      we want is the base file name (and extension(s), if present) but not the
      prefix, unless it's an absolute path.
    */
  if (!fileAbsPath(fileName)) {
    std::string::size_type pos = fileName.rfind(field);
    cbcSettings->lastMpsIn_ = fileName.substr(pos);
  } else {
    cbcSettings->lastMpsIn_ = fileName;
  }
  /*
      Load the primal variable values into the debug solution vector.
    */
  int intUnused, numCols;
  double dblUnused;
  double *primals;

  bool readOK = readSolution(fileName,
    intUnused, numCols, dblUnused, 0, 0, &primals, 0);

  if (readOK) {
    if (cbcSettings->debugSol_.values_) {
      delete[] cbcSettings->debugSol_.values_;
    }
    cbcSettings->debugSol_.numCols_ = numCols;
    cbcSettings->debugSol_.values_ = primals;
    retval = 0;
  } else {
    if (primals) {
      delete[] primals;
    }
  }

  return (retval);
}

//###########################################################################
//###########################################################################
   
/*
  Utility routine to save the current solution to a file. No formatting, and
  not intended to be portable in any way, shape, or form.
*/

void saveSolution(const OsiSolverInterface *osi, std::string fileName)

{
  FILE *fp = fopen(fileName.c_str(), "wb");

  if (fp) {
    int numberRows = osi->getNumRows();
    int numberColumns = osi->getNumCols();
    double objectiveValue = osi->getObjValue();

    fwrite(&numberRows, sizeof(int), 1, fp);
    fwrite(&numberColumns, sizeof(int), 1, fp);
    fwrite(&objectiveValue, sizeof(double), 1, fp);

    const double *primalRowSolution = osi->getRowActivity();
    const double *dualRowSolution = osi->getRowPrice();
    const double *primalColumnSolution = osi->getColSolution();
    const double *dualColumnSolution = osi->getReducedCost();

    fwrite(primalRowSolution, sizeof(double), numberRows, fp);
    fwrite(dualRowSolution, sizeof(double), numberRows, fp);
    fwrite(primalColumnSolution, sizeof(double), numberColumns, fp);
    fwrite(dualColumnSolution, sizeof(double), numberColumns, fp);

    fclose(fp);
  } else {
    std::cout
      << "saveSolution: Unable to open file `"
      << fileName << "'." << std::endl;
  }

  return;
}

//###########################################################################
//###########################################################################
   
/*
  Utility routine to read in a solution dump created by saveSolution. Generally
  we don't need all the info in this file, so the routine accepts a bunch of
  reference/pointer paramaters and fills in any that are non-null. It's the
  client's responsibility to dispose of space allocated for solution vectors.
  The parameters fileName, numRows, numCols, and objVal are mandatory. The rest
  can be null.
*/
bool readSolution(std::string fileName,
  int &numRows, int &numCols, double &objVal,
  double **rowActivity, double **dualVars,
  double **primalVars, double **reducedCosts)

{
  FILE *fp = fopen(fileName.c_str(), "rb");
  bool retval = true;

  numRows = -1;
  numCols = -1;
  objVal = 0;
  *rowActivity = 0;
  *dualVars = 0;
  *primalVars = 0;
  *reducedCosts = 0;

  if (fp) {
    fread(&numRows, sizeof(int), 1, fp);
    fread(&numCols, sizeof(int), 1, fp);
    fread(&objVal, sizeof(double), 1, fp);

    if (rowActivity != NULL) {
      *rowActivity = new double[numRows];
      fread(*rowActivity, sizeof(double), numRows, fp);
    } else {
      fseek(fp, numRows * sizeof(double), SEEK_CUR);
    }

    if (dualVars != NULL) {
      *dualVars = new double[numRows];
      fread(*dualVars, sizeof(double), numRows, fp);
    } else {
      fseek(fp, numRows * sizeof(double), SEEK_CUR);
    }

    if (primalVars != NULL) {
      *primalVars = new double[numCols];
      fread(*primalVars, sizeof(double), numCols, fp);
    } else {
      fseek(fp, numCols * sizeof(double), SEEK_CUR);
    }

    if (reducedCosts != NULL) {
      *reducedCosts = new double[numCols];
      fread(*reducedCosts, sizeof(double), numCols, fp);
    } else {
      fseek(fp, numCols * sizeof(double), SEEK_CUR);
    }

    fclose(fp);
  } else {
    std::cout
      << "readSolution: Unable to open file `"
      << fileName << "'." << std::endl;
    retval = false;
  }

  return (retval);
}

} // end namespace CbcSolverParamUtils

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
