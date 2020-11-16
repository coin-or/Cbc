/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/

#include "CoinPragma.hpp"

#include <cassert>

#include "CbcSettings.hpp"
#include "CbcParamUtils.hpp"

/*
  Constructor for settings class.

  Set up defaults. Note that prototypes for
  cut generators and heuristics will be created on demand; see the access
  functions.

  Once this structure settles down, simple intialisation should move up to
  the standard `:' block. In the meantime, this avoids complaints about
  ordering.
*/

CbcSettings::CbcSettings() : parameters_(CBC_LASTPARAM)
{
  version_ = CBC_VERSION;
  /*
      It's unclear to me that this is a good choice for dfltDirectory. Makes
      sense for commands, but seems unnecessary for data files. Perhaps a null
      string instead?
    */
  char dirsep = CoinFindDirSeparator();
  dfltDirectory_ = (dirsep == '/' ? "./" : ".\\");
  lastMpsIn_ = "";
  allowImportErrors_ = false;
  lastSolnOut_ = "stdout";
  printMode_ = 0;
  printMask_ = "";
  noPrinting_ = false;
  printWelcome_ = true;
  useSignalHandler_ = false;

  verbose_ = 0;
  paramsProcessed_ = 0;
  defaultSettings_ = true;

  debugCreate_ = "";
  debugFile_ = "";
  debugSol_.numCols_ = -1;
  debugSol_.values_ = 0;

  printOpt_ = 0;

  /*
      Assigning us_en to cur_lang_ is entirely bogus, but CoinMessages::Language
      does not provide an `unspecified' code.
    */
  msgHandler_ = new CoinMessageHandler();
  ourMsgHandler_ = true;
  cur_lang_ = CoinMessages::us_en;
  msgs_ = 0;
  logLvl_ = 0;

  totalTime_ = 0.0;

  model_ = 0;
  dfltSolver_ = 0;
  goodModel_ = false;
  bab_.majorStatus_ = CbcSettings::BACInvalid;
  bab_.minorStatus_ = CbcSettings::BACmInvalid;
  bab_.where_ = CbcSettings::BACwInvalid;
  bab_.haveAnswer_ = false;
  bab_.answerSolver_ = 0;

  preProcess_ = CbcSettings::IPPSOS;
  cutDepth_ = -1;

  probing_.mode_ = CbcSettings::CGIfMove;
  probing_.proto_ = 0;
  probing_.usingObjective_ = true;
  probing_.maxPass_ = 3;
  probing_.maxPassRoot_ = 3;
  probing_.maxProbe_ = 10;
  probing_.maxProbeRoot_ = 50;
  probing_.maxLook_ = 10;
  probing_.maxLookRoot_ = 50;
  probing_.maxElements_ = 200;
  probing_.rowCuts_ = 3;

  clique_.mode_ = CbcSettings::CGIfMove;
  clique_.proto_ = 0;
  clique_.starCliqueReport_ = false;
  clique_.rowCliqueReport_ = false;
  clique_.minViolation_ = 0.1;

  flow_.mode_ = CbcSettings::CGIfMove;
  flow_.proto_ = 0;

  gomory_.mode_ = CbcSettings::CGIfMove;
  gomory_.proto_ = 0;
  gomory_.limit_ = 50;
  gomory_.limitAtRoot_ = 512;

  knapsack_.mode_ = CbcSettings::CGIfMove;
  knapsack_.proto_ = 0;

  // landp_mode_ = CbcSettings::CGOff ;
  // landp_.proto_ = 0 ;

  mir_.mode_ = CbcSettings::CGIfMove;
  mir_.proto_ = 0;

#if 0
  oddHole_.mode_ = CbcSettings::CGOff;
  oddHole_.proto_ = 0;
#endif

  redSplit_.mode_ = CbcSettings::CGRoot;
  redSplit_.proto_ = 0;

  twomir_.mode_ = CbcSettings::CGRoot;
  twomir_.proto_ = 0;
  twomir_.maxElements_ = 250;

  fpump_.mode_ = CbcSettings::HeurOn;
  fpump_.proto_ = 0;
  fpump_.initialTune_ = -1;

  combine_.mode_ = CbcSettings::HeurOn;
  combine_.proto_ = 0;
  combine_.trySwap_ = 1;

  greedyCover_.mode_ = CbcSettings::HeurOn;
  greedyCover_.proto_ = 0;
  greedyEquality_.mode_ = CbcSettings::HeurOn;
  greedyEquality_.proto_ = 0;

  localTree_.mode_ = CbcSettings::HeurOff;
  localTree_.proto_ = 0;
  localTree_.soln_ = 0;
  localTree_.range_ = 10;
  localTree_.typeCuts_ = 0;
  localTree_.maxDiverge_ = 0;
  localTree_.timeLimit_ = 10000;
  localTree_.nodeLimit_ = 2000;
  localTree_.refine_ = true;

  rounding_.mode_ = CbcSettings::HeurOn;
  rounding_.proto_ = 0;

  djFix_.mode_ = CbcSettings::ParamOff;
  djFix_.threshold_ = 1.0e100;

  artVar_.mode_ = CbcSettings::ParamOff;
  artVar_.threshold_ = 0.0;

  priorityMode_ = CbcSettings::BPOff;
  /*
      The value for numBeforeTrust is as recommended by Achterberg. Cbc's
      implementation doesn't really have a parameter equivalent to Achterberg's
      dynamic limit on number of strong branching evaluations, so go with a
     fairly large default. As of 06.12.16, the magic number for shadow price
     mode meant `use shadow prices (penalties, I think) if there's no strong
     branching info'.
    */
  chooseStrong_.numBeforeTrust_ = 8;
  chooseStrong_.numStrong_ = 100;
  chooseStrong_.shadowPriceMode_ = 1;

  addCbcParams();
  addCbcModelParams();

  establishClpParams(clpParameters_);
  
  return;
}

/*
  Note that we don't want to delete dfltSolver_ here because it's just a
  copy of the pointer held in the solvers map over in CbcGenSolvers.cpp.
*/
CbcSettings::~CbcSettings() {
   // TODO Do we own pointer here?
   if (model_)
    delete model_;
  if (bab_.answerSolver_)
    delete bab_.answerSolver_;

  if (probing_.proto_)
    delete probing_.proto_;
  if (clique_.proto_)
    delete clique_.proto_;
  if (flow_.proto_)
    delete flow_.proto_;
  if (gomory_.proto_)
    delete gomory_.proto_;
  if (knapsack_.proto_)
    delete knapsack_.proto_;
  if (mir_.proto_)
    delete mir_.proto_;
#if 0
  if (oddHole_.proto_)
    delete oddHole_.proto_;
#endif
  if (redSplit_.proto_)
    delete redSplit_.proto_;
  if (twomir_.proto_)
    delete twomir_.proto_;

  if (fpump_.proto_)
    delete fpump_.proto_;
  if (combine_.proto_)
    delete combine_.proto_;
  if (greedyCover_.proto_)
    delete greedyCover_.proto_;
  if (greedyEquality_.proto_)
    delete greedyEquality_.proto_;
  if (rounding_.proto_)
    delete rounding_.proto_;

  if (msgHandler_ && ourMsgHandler_)
    delete msgHandler_;
  if (msgs_)
    delete msgs_;

  return;
}

//###########################################################################
//###########################################################################

// some help strings that repeat for many options
#define CUTS_LONGHELP                                                          \
  "Value 'on' enables the cut generator and CBC will try it in the branch "    \
  "and cut tree (see cutDepth on how to fine tune the behavior). Value "       \
  "'root' lets CBC run the cut generator generate only at the root node. "     \
  "Value 'ifmove' lets CBC use the cut generator in the tree if it looks as "  \
  "if it is doing some good and moves the objective value. Value 'forceon' "   \
  "turns on the cut generator and forces CBC to use it at every node."

#define HEURISTICS_LONGHELP                                                    \
  "Value 'on' means to use the heuristic in each node of the tree, i.e. "      \
  "after preprocessing. Value 'before' means use the heuristic only if "       \
  "option doHeuristics is used. Value 'both' means to use the heuristic if "   \
  "option doHeuristics is used and during solve."

/*
  Function to add Cbc parameters to the Cbc parameter
  vector. Where needed, defaults are drawn from CbcSettings.
  This function is a friend of CbcSettings.
*/

void CbcSettings::addCbcParams() {

  addCbcSolverStrParams();
  addCbcSolverHelpParams();
  addCbcSolverActionParams();
  addCbcSolverKwdParams();
  addCbcSolverDblParams();
  addCbcSolverIntParams();
  addCbcSolverBoolParams();
  addCbcSolverCutParams();
  addCbcSolverHeurParams();

  for (int code = CBC_FIRSTPARAM + 1;
       code < CBC_LASTPARAM; code++) {
    CoinParam &param = parameters_[code];
    static_cast<CbcParam &>(param).setSettings(this);
    static_cast<CbcParam &>(param).setModel(model_);
    static_cast<CbcParam &>(param).setParamCode(code);
  }

  return;
}

//###########################################################################
//###########################################################################

void CbcSettings::addCbcSolverStrParams() {

  parameters_[CSVSTATISTICS].setup(
      "csv!Statistics", "Create one line of statistics",
      dfltDirectory_,
      "This appends statistics to given file name.  It will use the default "
      "directory given by 'directory'.  A name of '$' will use the previous "
      "value for the name.  This is initialized to '', i.e. it must be set.  "
      "Adds header if file empty or does not exist.",
      CoinParam::displayPriorityLow);
  parameters_[CSVSTATISTICS].setPushFunc(CbcParamUtils::pushCbcSolverStrParam);

  parameters_[DEBUG].setup(
      "debug!In", "Read/write valid solution from/to file", "",
      "This will read a solution file from the given file name.  It will use "
      "the default directory given by 'directory'.  A name of '$' will use the "
      "previous value for the name.  This is initialized to '', i.e. it must "
      "be set.\n\nIf set to create it will create a file called debug.file "
      "after B&C search; if set to createAfterPre it will create the file "
      "before undoing preprocessing.\n\nThe idea is that if you suspect a bad "
      "cut generator and you did not use preprocessing you can do a good run "
      "with debug set to 'create' and then switch on the cuts you suspect and "
      "re-run with debug set to 'debug.file'  Similarly if you do use "
      "preprocessing, but use createAfterPre.  The create case has the same "
      "effect as saveSolution.",
      CoinParam::displayPriorityNone);
  parameters_[DEBUG].setPushFunc(CbcParamUtils::doDebugParam);

  parameters_[DIRECTORY].setup(
      "directory", "Set Default directory for import etc.",
      dfltDirectory_,
      "This sets the directory which import, export, saveModel, restoreModel "
      "etc. will use. It is initialized to the current directory.");
  parameters_[DIRECTORY].setPushFunc(CbcParamUtils::pushCbcSolverStrParam);

  parameters_[DIRSAMPLE].setup(
      "dirSample", "Set directory where the COIN-OR sample problems are.",
      dfltDirectory_,
      "This sets the directory where the COIN-OR sample problems reside. It is "
      "used only when -unitTest is passed to clp. clp will pick up the test "
      "problems from this directory. It is initialized to '../../Data/Sample'",
      CoinParam::displayPriorityLow);
  parameters_[DIRSAMPLE].setPushFunc(CbcParamUtils::pushCbcSolverStrParam);

  parameters_[DIRNETLIB].setup(
      "dirNetlib", "Set directory where the netlib problems are.",
      dfltDirectory_,
      "This sets the directory where the netlib problems reside. One can get "
      "the netlib problems from COIN-OR or from the main netlib site. This "
      "parameter is used only when -netlib is passed to cbc. cbc will pick up "
      "the netlib problems from this directory. If clp is built without zlib "
      "support then the problems must be uncompressed. It is initialized to "
      "'../../Data/Netlib'",
      CoinParam::displayPriorityLow);
  parameters_[DIRNETLIB].setPushFunc(CbcParamUtils::pushCbcSolverStrParam);

  parameters_[DIRMIPLIB].setup(
      "dirMiplib", "Set directory where the miplib 2003 problems are.",
      dfltDirectory_,
      "This sets the directory where the miplib 2003 problems reside. One can "
      "get the miplib problems from COIN-OR or from the main miplib site. This "
      "parameter is used only when -miplib is passed to cbc. cbc will pick up "
      "the miplib problems from this directory. If cbc is built without zlib "
      "support then the problems must be uncompressed. It is initialized to "
      "'../../Data/miplib3'",
      CoinParam::displayPriorityLow);
  parameters_[DIRMIPLIB].setPushFunc(CbcParamUtils::pushCbcSolverStrParam);

  parameters_[MIPSTART].setup(
      "mips!tart", "reads an initial feasible solution from file",
      std::string("mipstart.sln"),
      "The MIPStart allows one to enter an initial integer feasible solution "
      "to CBC. Values of the main decision variables which are active (have "
      "non-zero values) in this solution are specified in a text  file. The "
      "text file format used is the same of the solutions saved by CBC, but "
      "not all fields are required to be filled. First line may contain the "
      "solution status and will be ignored, remaining lines contain column "
      "indexes, names and values as in this example:\n\n Stopped on iterations "
      "- objective value 57597.00000000\n      0  x(1,1,2,2)               1 "
      "\n      1  x(3,1,3,2)               1 \n      5  v(5,1)                 "
      "  2 \n      33 x(8,1,5,2)               1 \n      ...\n\n Column "
      "indexes are also ignored since pre-processing can change them. There is "
      "no need to include values for continuous or integer auxiliary "
      "variables, since they can be computed based on main decision variables. "
      "Starting CBC with an integer feasible solution can dramatically improve "
      "its performance: several MIP heuristics (e.g. RINS) rely on having at "
      "least one feasible solution available and can start immediately if the "
      "user provides one. Feasibility Pump (FP) is a heuristic which tries to "
      "overcome the problem of taking too long to find feasible solution (or "
      "not finding at all), but it not always succeeds. If you provide one "
      "starting solution you will probably save some time by disabling FP. "
      "\n\n Knowledge specific to your problem can be considered to write an "
      "external module to quickly produce an initial feasible solution - some "
      "alternatives are the implementation of simple greedy heuristics or the "
      "solution (by CBC for example) of a simpler model created just to find a "
      "feasible solution. \n\n Silly options added.  If filename ends .low "
      "then integers not mentioned are set low - also .high, .lowcheap, "
      ".highcheap, .lowexpensive, .highexpensive where .lowexpensive sets "
      "costed ones to make expensive others low. Also if filename starts "
      "empty. then no file is read at all - just actions done. \n\n Question "
      "and suggestions regarding MIPStart can be directed to\n "
      "haroldo.santos@gmail.com. ");
  parameters_[MIPSTART].setPushFunc(CbcParamUtils::pushCbcSolverStrParam);

  parameters_[NEXTBESTSOLUTION].setup(
      "nextB!estSolution", "Prints next best saved solution to file", "",
      "To write best solution, just use solution.  This prints next best (if "
      "exists) and then deletes it. This will write a primitive solution file "
      "to the given file name.  It will use the default directory given by "
      "'directory'.  A name of '$' will use the previous value for the name.  "
      "This is initialized to 'stdout'.  The amount of output can be varied "
      "using printi!ngOptions or printMask.");
  parameters_[NEXTBESTSOLUTION].setPushFunc(CbcParamUtils::doNothingParam);

  parameters_[PRINTMASK].setup(
      "printM!ask", "Control printing of solution with a regular expression",
      "",
      "If set then only those names which match mask are printed in a "
      "solution. '?' matches any character and '*' matches any set of "
      "characters.  The default is '' (unset) so all variables are printed. "
      "This is only active if model has names.");
  parameters_[PRINTMASK].setPushFunc(CbcParamUtils::doPrintMaskParam);

  parameters_[SAVESOL].setup(
      "saveS!olution", "saves solution to file", std::string("solution.sln"),
      "This will write a binary solution file to the given file name.  It will "
      "use the default directory given by 'directory'.  A name of '$' will use "
      "the previous value for the name.  This is initialized to "
      "'solution.file'.  To read the file use fread(int) twice to pick up "
      "number of rows and columns, then fread(double) to pick up objective "
      "value, then pick up row activities, row duals, column activities and "
      "reduced costs - see bottom of CbcParam.cpp for code that reads or "
      "writes file. If name contains '_fix_read_' then does not write but "
      "reads and will fix all variables");
  parameters_[SAVESOL].setPushFunc(CbcParamUtils::pushCbcSolverStrParam);

  parameters_[SOLUTION].setup(
      "solu!tion", "Prints solution to file", std::string("stdout"),
      "This will write a primitive solution file to the given file name.  It "
      "will use the default directory given by 'directory'.  A name of '$' "
      "will use the previous value for the name.  This is initialized to "
      "'stdout'.  The amount of output can be varied using printi!ngOptions or "
      "printMask.");
  parameters_[SOLUTION].setPushFunc(CbcParamUtils::doSolutionParam);

  parameters_[PRIORITYIN].setup(
      "prio!rityIn", "Import priorities etc from file",
      std::string("priorities.txt"),
      "This will read a file with priorities from the given file name.  It "
      "will use the default directory given by 'directory'.  A name of '$' "
      "will use the previous value for the name.  This is initialized to '', "
      "i.e. it must be set.  This can not read from compressed files. File is "
      "in csv format with allowed headings - name, number, priority, "
      "direction, up, down, solution.  Exactly one of name and number must be "
      "given.");
  parameters_[PRIORITYIN].setPushFunc(CbcParamUtils::pushCbcSolverStrParam);
}

//###########################################################################
//###########################################################################

void CbcSettings::addCbcSolverHelpParams() {
  for (int code = CBC_FIRSTHELPPARAM + 1;
       code < CBC_LASTHELPPARAM; code++) {
    parameters_[code].setPushFunc(CbcParamUtils::doHelpParam);
  }
  parameters_[GENERALQUERY].setup(
      "?", "Print a list of commands", CoinParam::displayPriorityNone);

  parameters_[FULLGENERALQUERY].setup(
      "???", "Print a list with *all* commands, even those hidden with `?'",
      CoinParam::displayPriorityNone);

  // Need display parameter to resolve ambiguity
  parameters_[HELP].setup(
      "help", "Print out version, non-standard options and some help",
      "This prints out some help to get a user started. If you're seeing this "
      "message, you should be past that stage.",
      CoinParam::displayPriorityHigh);
}

//###########################################################################
//###########################################################################

void CbcSettings::addCbcSolverActionParams() {

  parameters_[DUMMY].setup(
      "sleep", "for debug", 0, 9999, 0,
      "If passed to solver from ampl, then ampl will wait so that you can copy "
      ".nl file for debug.",
      CoinParam::displayPriorityNone);

  parameters_[OUTDUPROWS].setup(
      "outDup!licates", "Takes duplicate rows, etc., out of the integer model",
      "", CoinParam::displayPriorityNone);

  parameters_[SHOWUNIMP].setup(
      "unimp!lemented", "Report unimplemented commands.", "",
      CoinParam::displayPriorityNone);
  parameters_[SHOWUNIMP].setPushFunc(CbcParamUtils::doUnimplementedParam);

  parameters_[STRENGTHEN].setup(
      "strengthen", "Create strengthened problem",
      "This creates a new problem by applying the root node cuts. All tight "
      "constraints will be in resulting problem.",
      CoinParam::displayPriorityHigh);

  parameters_[BAB].setup(
      "branch!AndCut", "Do Branch and Cut",
      "This does branch and cut. There are many parameters which can affect "
      "the performance.  First just try with default cbcSettings and look "
      "carefully at the log file.  Did cuts help?  Did they take too long?  "
      "Look at output to see which cuts were effective and then do some "
      "tuning.  You will see that the options for cuts are off, on, root and "
      "ifmove.  Off is obvious, on means that this cut generator will be tried "
      "in the branch and cut tree (you can fine tune using 'depth').  Root "
      "means just at the root node while 'ifmove' means that cuts will be used "
      "in the tree if they look as if they are doing some good and moving the "
      "objective value.  If pre-processing reduced the size of the problem or "
      "strengthened many coefficients then it is probably wise to leave it on. "
      " Switch off heuristics which did not provide solutions.  The other "
      "major area to look at is the search.  Hopefully good solutions were "
      "obtained fairly early in the search so the important point is to select "
      "the best variable to branch on.  See whether strong branching did a "
      "good job - or did it just take a lot of iterations.  Adjust the "
      "strongBranching and trustPseudoCosts parameters.",
      CoinParam::displayPriorityHigh);
  parameters_[BAB].setPushFunc(CbcParamUtils::doBaCParam);

  parameters_[ENVIRONMENT].setup(
      "environ!ment", "Read commands from environment",
      "This starts reading from environment variable COIN_ENVIRONMENT.",
      CoinParam::displayPriorityNone);
  parameters_[ENVIRONMENT].setPushFunc(CbcParamUtils::doNothingParam);

  parameters_[EXIT].setup(
      "end", "Stops execution",
      "This stops execution; end, exit, quit and stop are synonyms.",
      CoinParam::displayPriorityHigh);
  parameters_[EXIT].setPushFunc(CbcParamUtils::doExitParam);

  parameters_[EXPORT].setup(
      "export", "Export model as mps file", std::string("default.mps"),
      "This will write an MPS format file to the given file name.  It will use "
      "the default directory given by 'directory'.  A name of '$' will use the "
      "previous value for the name.  This is initialized to 'default.mps'. It "
      "can be useful to get rid of the original names and go over to using "
      "Rnnnnnnn and Cnnnnnnn.  This can be done by setting 'keepnames' off "
      "before importing mps file.",
      CoinParam::displayPriorityHigh);
  parameters_[EXPORT].setPushFunc(CbcParamUtils::pushCbcSolverStrParam);

  parameters_[IMPORT].setup(
      "import", "Import model from file", lastMpsIn_,
      "This will read an MPS format file from the given file name.  It will "
      "use the default directory given by 'directory'.  A name of '$' will use "
      "the previous value for the name.  This is initialized to '', i.e., it "
      "must be set.  If you have libgz then it can read compressed files "
      "'xxxxxxxx.gz'.",
      CoinParam::displayPriorityHigh);
  parameters_[IMPORT].setPushFunc(CbcParamUtils::doImportParam);

  parameters_[MIPLIB].setup("miplib", "Do some of miplib test set", "",
                            CoinParam::displayPriorityHigh);

  parameters_[PRINTVERSION].setup(
      "version", "Print version", "", CoinParam::displayPriorityHigh);
  parameters_[PRINTVERSION].setPushFunc(CbcParamUtils::doVersionParam);

  parameters_[SOLVECONTINUOUS].setup(
      "initialS!olve", "Solve to continuous optimum",
      "This just solves the problem to the continuous optimum, without adding "
      "any cuts.",
      CoinParam::displayPriorityHigh);

  parameters_[STATISTICS].setup(
      "stat!istics", "Print some statistics",
      "This command prints some statistics for the current model. If log level "
      ">1 then more is printed. These are for presolved model if presolve on "
      "(and unscaled).",
      CoinParam::displayPriorityHigh);

#if 0
  // Need to figure out what to do here. Same parameter can't have two names...
  parameters_[STDIN].setup( "-", "Switch to interactive command line mode", ""
                            CoinParam::displayPriorityNone);
  parameters_[STDIN].setPushFunc(CbcParamUtils::doNothingParam);
#endif

  parameters_[STDIN].setup(
      "stdin", "Switch to interactive command line mode", "",
      CoinParam::displayPriorityNone);
  parameters_[STDIN].setPushFunc(CbcParamUtils::doNothingParam);

  parameters_[UNITTEST].setup(
      "unitTest", "Do unit test", "This exercises the unit test.", "",
      CoinParam::displayPriorityHigh);
}

//###########################################################################
//###########################################################################

void CbcSettings::addCbcSolverKwdParams() {
  for (int code = CBC_FIRSTKWDPARAM + 1;
       code < CBC_LASTKWDPARAM; code++) {
    parameters_[code].setPushFunc(CbcParamUtils::pushCbcSolverKwdParam);
  }

  parameters_[COMMANDPRINTLEVEL].setup(
      "allC!ommands", "What priority level of commands to print", "high",
      CbcSettings::displayHigh,
      "For the sake of your sanity, only the more useful and simple commands "
      "are printed out on ?.",
      CoinParam::displayPriorityNone);
  parameters_[COMMANDPRINTLEVEL].appendKwd("all", CbcSettings::displayAll);
  parameters_[COMMANDPRINTLEVEL].appendKwd("highlow", CbcSettings::displayLowHigh);

  parameters_[CLQSTRENGTHENING].setup(
      "clqstr!engthen",
      "Whether and when to perform Clique Strengthening preprocessing routine",
      "after", CbcSettings::ClqStrAfter);
  parameters_[CLQSTRENGTHENING].appendKwd("off", CbcSettings::ClqStrOff);
  parameters_[CLQSTRENGTHENING].appendKwd("before", CbcSettings::ClqStrBefore);

  parameters_[BRANCHPRIORITY].setup(
      "cost!Strategy", "Whether to use costs or column order as priorities",
      "off", CbcSettings::BPOff,
      "This orders the variables in order of their absolute costs - with "
      "largest cost ones being branched on first.  This primitive strategy can "
      "be surprisingly effective.  The column order option is obviously not on "
      "costs but it's easy to implement.");
  parameters_[BRANCHPRIORITY].appendKwd("pri!orities", CbcSettings::BPCost);
  parameters_[BRANCHPRIORITY].appendKwd("column!Order", CbcSettings::BPOrder);

  parameters_[CUTOFFCONSTRAINT].setup(
      "constraint!fromCutoff", "Whether to use cutoff as constraint", "off",
      CbcSettings::COOff,
      "For some problems, cut generators and general branching work better if "
      "the problem would be infeasible if the cost is too high. "
      "If this option is enabled, the objective function is added as a "
      "constraint which right hand side is set to the current cutoff value "
      "(objective value of best known solution)");
  parameters_[CUTOFFCONSTRAINT].appendKwd("on", CbcSettings::COOn);
  parameters_[CUTOFFCONSTRAINT].appendKwd("variable", CbcSettings::COVariable);
  parameters_[CUTOFFCONSTRAINT].appendKwd("forcevariable", CbcSettings::COForceVariable);
  parameters_[CUTOFFCONSTRAINT].appendKwd("conflict", CbcSettings::COConflict);

  parameters_[INTPRINT].setup(
      "printi!ngOptions", "Print options", "normal", CbcSettings::PMNormal,
      "This changes the amount and format of printing a solution:\n normal - "
      "nonzero column variables \ninteger - nonzero integer column variables\n "
      "special - in format suitable for OsiRowCutDebugger\n rows - nonzero "
      "column variables and row activities\n all - all column variables and "
      "row activities.\n\n For non-integer problems 'integer' and 'special' "
      "act like 'normal'.  Also see printMask for controlling output.");
  parameters_[INTPRINT].appendKwd("integer", CbcSettings::PMInteger);
  parameters_[INTPRINT].appendKwd("special", CbcSettings::PMSpecial);
  parameters_[INTPRINT].appendKwd("rows", CbcSettings::PMRows);
  parameters_[INTPRINT].appendKwd("all", CbcSettings::PMAll);

  parameters_[NODESTRATEGY].setup(
      "node!Strategy",
      "What strategy to use to select the next node from the branch and cut "
      "tree",
      "hybrid", CbcSettings::NSHybrid,
      "Normally before a feasible solution is found, CBC will choose a node "
      "with fewest infeasibilities. Alternatively, one may choose tree-depth "
      "as the criterion. This requires the minimal amount of memory, but may "
      "take a long time to find the best solution. Additionally, one may "
      "specify whether up or down branches must be selected first (the up-down "
      "choice will carry on after a first solution has been bound). The choice "
      "'hybrid' does breadth first on small depth nodes and then switches to "
      "'fewest'.");
  parameters_[NODESTRATEGY].appendKwd("fewest", CbcSettings::NSFewest);
  parameters_[NODESTRATEGY].appendKwd("depth", CbcSettings::NSDepth);
  parameters_[NODESTRATEGY].appendKwd("upfewest", CbcSettings::NSUpFewest);
  parameters_[NODESTRATEGY].appendKwd("downfewest", CbcSettings::NSDownFewest);
  parameters_[NODESTRATEGY].appendKwd("updepth", CbcSettings::NSUpDepth);
  parameters_[NODESTRATEGY].appendKwd("downdepth", CbcSettings::NSDownDepth);

  parameters_[ORBITAL].setup(
      "Orbit!alBranching", "Whether to try orbital branching", "off",
      CbcSettings::OBOff,
      "This switches on Orbital branching. Value 'on' just adds orbital, "
      "'strong' tries extra fixing in strong branching.");
  parameters_[ORBITAL].appendKwd("on", CbcSettings::OBOn);
  parameters_[ORBITAL].appendKwd("slowish", CbcSettings::OBSlowish);
  parameters_[ORBITAL].appendKwd("strong", CbcSettings::OBStrong);
  parameters_[ORBITAL].appendKwd("force", CbcSettings::OBForce);
  parameters_[ORBITAL].appendKwd("simple", CbcSettings::OBSimple);
  parameters_[ORBITAL].appendKwd("moreprinting", CbcSettings::OBMorePrinting);

  parameters_[PREPROCESS].setup(
      "preprocess", "Whether to use integer preprocessing", "off",
      CbcSettings::IPPOff,
      "This tries to reduce size of the model in a similar way to presolve and "
      "it also tries to strengthen the model. This can be very useful and is "
      "worth trying.  save option saves on file presolved.mps.  equal will "
      "turn <= cliques into ==.  sos will create sos sets if all 0-1 in sets "
      "(well one extra is allowed) and no overlaps.  trysos is same but allows "
      "any number extra. equalall will turn all valid inequalities into "
      "equalities with integer slacks. strategy is as on but uses "
      "CbcStrategy.");
  parameters_[PREPROCESS].appendKwd("on", CbcSettings::IPPOn);
  parameters_[PREPROCESS].appendKwd("save", CbcSettings::IPPSave);
  parameters_[PREPROCESS].appendKwd("equal", CbcSettings::IPPEqual);
  parameters_[PREPROCESS].appendKwd("sos", CbcSettings::IPPSOS);
  parameters_[PREPROCESS].appendKwd("trysos", CbcSettings::IPPTrySOS);
  parameters_[PREPROCESS].appendKwd("equalall", CbcSettings::IPPEqualAll);
  parameters_[PREPROCESS].appendKwd("strategy", CbcSettings::IPPStrategy);

  parameters_[SOSPRIORITIZE].setup(
      "sosP!rioritize", "How to deal with SOS priorities", "off",
      CbcSettings::SOSOff,
      "This sets priorities for SOS.  Values 'high' and 'low' just set a "
      "priority relative to the for integer variables.  Value 'orderhigh' "
      "gives first highest priority to the first SOS and integer variables a "
      "low priority.  Value 'orderlow' gives integer variables a high priority "
      "then SOS in order.");
  parameters_[SOSPRIORITIZE].appendKwd("high", CbcSettings::SOSHigh);
  parameters_[SOSPRIORITIZE].appendKwd("low", CbcSettings::SOSLow);
  parameters_[SOSPRIORITIZE].appendKwd("orderhigh", CbcSettings::SOSOrderHigh);
  parameters_[SOSPRIORITIZE].appendKwd("orderlow", CbcSettings::SOSOrderLow);

  parameters_[STRATEGY].setup(
      "strat!egy", "Switches on groups of features", "default",
      CbcSettings::StrategyDefault,
      "This turns on newer features. Use 0 for easy problems, 1 is default, 2 "
      "is aggressive. 1 uses Gomory cuts with a tolerance of 0.01 at the root "
      "node, does a possible restart after 100 nodes if many variables could "
      "be fixed, activates a diving and RINS heuristic, and makes the "
      "feasibility pump more aggressive."); // This does not apply to unit tests
                                            // (where 'experiment' may have
                                            // similar effects)
  parameters_[STRATEGY].appendKwd("easy", CbcSettings::StrategyEasy);
  parameters_[STRATEGY].appendKwd("aggressive", CbcSettings::StrategyAggressive);

  parameters_[TIMEMODE].setup(
      "timeM!ode", "Whether to use CPU or elapsed time", "cpu",
      CbcSettings::ClockCpu,
      "cpu uses CPU time for stopping, while elapsed uses elapsed time. (On "
      "Windows, elapsed time is always used).");
  parameters_[TIMEMODE].appendKwd("elapsed", CbcSettings::ClockElapsed);

  parameters_[USECGRAPH].setup(
      "cgraph",
      "Whether to use the conflict graph-based preprocessing and cut "
      "separation routines.",
      "on", CbcSettings::CGraphOn,
      "This switches the conflict graph-based preprocessing and cut separation "
      "routines (CglBKClique, CglOddWheel and CliqueStrengthening) on or off. "
      "Values: \n\t off: turns these routines off;\n\t on: turns these "
      "routines on; \n\t clq: turns these routines off and enables the cut "
      "separator of CglClique.");
  parameters_[USECGRAPH].appendKwd("off", CbcSettings::CGraphOff);
  parameters_[USECGRAPH].appendKwd("clq", CbcSettings::CGraphClique);
}

//###########################################################################
//###########################################################################

void CbcSettings::addCbcSolverDblParams() {
  for (int code = CBC_FIRSTDBLPARAM + 1;
       code < CBC_LASTDBLPARAM; code++) {
    parameters_[code].setPushFunc(CbcParamUtils::pushCbcSolverDblParam);
  }

  parameters_[ARTIFICIALCOST].setup(
      "artif!icialCost",
      "Costs >= this treated as artificials in feasibility pump", 0.0,
      COIN_DBL_MAX, getArtVarThreshold(), "",
      CoinParam::displayPriorityLow);

  parameters_[DEXTRA3].setup(
      "dextra3", "Extra double parameter 3", -COIN_DBL_MAX, COIN_DBL_MAX, 0.0,
      "", CoinParam::displayPriorityNone);

  parameters_[DEXTRA4].setup(
      "dextra4", "Extra double parameter 4", -COIN_DBL_MAX, COIN_DBL_MAX, 0.0,
      "", CoinParam::displayPriorityNone);

  parameters_[DEXTRA5].setup(
      "dextra5", "Extra double parameter 5", -COIN_DBL_MAX, COIN_DBL_MAX, 0.0,
      "", CoinParam::displayPriorityNone);

  parameters_[DJFIX].setup(
      "fix!OnDj", "Try heuristic that fixes variables based on reduced costs",
      -1.0e20, 1.0e20, getDjFixThreshold(),
      "If set, integer variables with reduced costs greater than the specified "
      "value will be fixed before branch and bound - use with extreme "
      "caution!");

  parameters_[FAKECUTOFF].setup(
      "pumpC!utoff", "Fake cutoff for use in feasibility pump", -COIN_DBL_MAX,
      COIN_DBL_MAX, 0.0,
      "A value of 0.0 means off. Otherwise, add a constraint forcing objective "
      "below this value in feasibility pump",
      CoinParam::displayPriorityLow);

  parameters_[FAKEINCREMENT].setup(
      "pumpI!ncrement", "Fake increment for use in feasibility pump",
      -COIN_DBL_MAX, COIN_DBL_MAX, 0.0,
      "A value of 0.0 means off. Otherwise, add a constraint forcing objective "
      "below this value in feasibility pump",
      CoinParam::displayPriorityLow);

  parameters_[SMALLBAB].setup(
      "fraction!forBAB", "Fraction in feasibility pump", 1.0e-5, 1.1, 0.5,
      "After a pass in the feasibility pump, variables which have not moved "
      "about are fixed and if the preprocessed model is smaller than this "
      "fraction of the original problem, a few nodes of branch and bound are "
      "done on the reduced problem.",
      CoinParam::displayPriorityLow);

  parameters_[TIGHTENFACTOR].setup(
      "tighten!Factor",
      "Tighten bounds using value times largest activity at continuous "
      "solution",
      0.0, COIN_DBL_MAX, 0.0, "This sleazy trick can help on some problems.");
}

//###########################################################################
//###########################################################################

void CbcSettings::addCbcSolverIntParams() {
  for (int code = CBC_FIRSTINTPARAM + 1;
       code < CBC_LASTINTPARAM; code++) {
    parameters_[code].setPushFunc(CbcParamUtils::pushCbcSolverIntParam);
  }

  parameters_[BKPIVOTINGSTRATEGY].setup(
      "bkpivot!ing", "Pivoting strategy used in Bron-Kerbosch algorithm", 0, 6,
      3);

  parameters_[BKMAXCALLS].setup(
      "bkmaxcalls",
      "Maximum number of recursive calls made by Bron-Kerbosch algorithm", 1,
      COIN_INT_MAX, 1000);

  parameters_[BKCLQEXTMETHOD].setup(
      "bkclqext!method",
      "Strategy used to extend violated cliques found by BK Clique Cut "
      "Separation routine",
      0, 5, 4,
      "Sets the method used in the extension module of BK Clique Cut "
      "Separation routine: 0=no extension; 1=random; 2=degree; 3=modified "
      "degree; 4=reduced cost(inversely proportional); 5=reduced "
      "cost(inversely proportional) + modified degree");

  parameters_[CPP].setup(
      "cpp!Generate", "Generates C++ code", 0, 4, 0,
      "Once you like what the stand-alone solver does then this allows you to "
      "generate user_driver.cpp which approximates the code.  0 gives simplest "
      "driver, 1 generates saves and restores, 2 generates saves and restores "
      "even for variables at default value. 4 bit in cbc generates size "
      "dependent code rather than computed values.");

  parameters_[CUTDEPTH].setup(
      "cutD!epth", "Depth in tree at which to do cuts", -1, 999999,
      getCutDepth(),
      "Cut generators may be off, on only at the root, on if they look useful, "
      "and on at some interval.  If they are done every node then that is "
      "that, but it may be worth doing them every so often.  The original "
      "method was every so many nodes but it is more logical to do it whenever "
      "depth in tree is a multiple of K.  This option does that and defaults "
      "to -1 (off).");

  parameters_[CUTLENGTH].setup(
      "cutL!ength", "Length of a cut", -1, COIN_INT_MAX, -1,
      "At present this only applies to Gomory cuts. -1 (default) leaves as is. "
      "Any value >0 says that all cuts <= this length can be generated both at "
      "root node and in tree. 0 says to use some dynamic lengths.  If value "
      ">=10,000,000 then the length in tree is value%10000000 - so 10000100 "
      "means unlimited length at root and 100 in tree.");

  parameters_[CUTPASSINTREE].setup(
      "passT!reeCuts",
      "Number of rounds that cut generators are applied in the tree",
      -COIN_INT_MAX, COIN_INT_MAX);

  parameters_[DEPTHMINIBAB].setup(
      "depth!MiniBab", "Depth at which to try mini branch-and-bound",
      -COIN_INT_MAX, COIN_INT_MAX, -1,
      "Rather a complicated parameter but can be useful. -1 means off for "
      "large problems but on as if -12 for problems where rows+columns<500, -2 "
      "means use Cplex if it is linked in.  Otherwise if negative then go into "
      "depth first complete search fast branch and bound when depth>= -value-2 "
      "(so -3 will use this at depth>=1).  This mode is only switched on after "
      "500 nodes.  If you really want to switch it off for small problems then "
      "set this to -999.  If >=0 the value doesn't matter very much.  The code "
      "will do approximately 100 nodes of fast branch and bound every now and "
      "then at depth>=5. The actual logic is too twisted to describe here.");

  parameters_[DIVEOPT].setup(
      "diveO!pt", "Diving options", -1, 20, -1,
      "If >2 && <20 then modify diving options -	 \n\t3 only at root "
      "and if no solution,	 \n\t4 only at root and if this heuristic has "
      "not got solution,	 \n\t5 decay only if no solution,	 \n\t6 "
      "if depth <3 or decay,	 \n\t7 run up to 2 times if solution found 4 "
      "otherwise,	 \n\t>10 All only at root (DivingC normal as "
      "value-10),	 \n\t>20 All with value-20).",
      CoinParam::displayPriorityLow);

  parameters_[DIVEOPTSOLVES].setup(
      "diveS!olves", "Diving solve option", -1, 200000, 100,
      "If >0 then do up to this many solves. However, the last digit is "
      "ignored and used for extra options: 1-3 enables fixing of satisfied "
      "integer variables (but not at bound), where 1 switches this off for "
      "that dive if the dive goes infeasible, and 2 switches it off "
      "permanently if the dive goes infeasible.",
      CoinParam::displayPriorityLow);

  parameters_[EXPERIMENT].setup(
      "exper!iment", "Whether to use testing features", -1, 200000, 0,
      "Defines how adventurous you want to be in using new ideas. 0 then no "
      "new ideas, 1 fairly sensible, 2 a bit dubious, 3 you are on your own!",
      CoinParam::displayPriorityLow);

  parameters_[EXTRA1].setup(
      "extra1", "Extra integer parameter 1", -COIN_INT_MAX, COIN_INT_MAX, -1,
      "", CoinParam::displayPriorityLow);

  parameters_[EXTRA2].setup(
      "extra2", "Extra integer parameter 2", -COIN_INT_MAX, COIN_INT_MAX, -1,
      "", CoinParam::displayPriorityLow);

  parameters_[EXTRA3].setup(
      "extra3", "Extra integer parameter 3", -COIN_INT_MAX, COIN_INT_MAX, -1,
      "", CoinParam::displayPriorityLow);

  parameters_[EXTRA4].setup(
      "extra4", "Extra integer parameter 4", -COIN_INT_MAX, COIN_INT_MAX, -1,
      "", CoinParam::displayPriorityLow);

  parameters_[FPUMPITS].setup(
      "passF!easibilityPump", "How many passes in feasibility pump", 0, 10000,
      getFeasPumpIters(),
      "This fine tunes the Feasibility Pump heuristic by doing more or fewer "
      "passes.");

  parameters_[FPUMPTUNE].setup(
      "pumpT!une", "Dubious ideas for feasibility pump", 0, 100000000, 0,
      "This fine tunes Feasibility Pump     \n\t>=10000000 use as objective "
      "weight switch     \n\t>=1000000 use as accumulate switch     \n\t>=1000 "
      "use index+1 as number of large loops     \n\t==100 use objvalue "
      "+0.05*fabs(objvalue) as cutoff OR fakeCutoff if set     \n\t%100 == "
      "10,20 affects how each solve is done     \n\t1 == fix ints at bounds, 2 "
      "fix all integral ints, 3 and continuous at bounds. If accumulate is on "
      "then after a major pass, variables which have not moved are fixed and a "
      "small branch and bound is tried.");

  parameters_[FPUMPTUNE2].setup(
      "moreT!une", "Yet more dubious ideas for feasibility pump", 0, 100000000,
      0,
      "Yet more ideas for Feasibility Pump     \n\t/100000 == 1 use box "
      "constraints and original obj in cleanup     \n\t/1000 == 1 Pump will "
      "run twice if no solution found     \n\t/1000 == 2 Pump will only run "
      "after root cuts if no solution found     \n\t/1000 >10 as above but "
      "even if solution found     \n\t/100 == 1,3.. exact 1.0 for objective "
      "values     \n\t/100 == 2,3.. allow more iterations per pass     \n\t n "
      "fix if value of variable same for last n iterations.",
      CoinParam::displayPriorityNone);

  parameters_[HEUROPTIONS].setup(
      "hOp!tions", "Heuristic options", -COIN_INT_MAX, COIN_INT_MAX, 0,
      "Value 1 stops heuristics immediately if the allowable gap has been "
      "reached. Other values are for the feasibility pump - 2 says do exact "
      "number of passes given, 4 only applies if an initial cutoff has been "
      "given and says relax after 50 passes, while 8 will adapt the cutoff rhs "
      "after the first solution if it looks as if the code is stalling.",
      CoinParam::displayPriorityLow);

  parameters_[MAXHOTITS].setup(
      "hot!StartMaxIts", "Maximum iterations on hot start",
      0, COIN_INT_MAX);

  parameters_[LOGLEVEL].setup(
      "log!Level", "Level of detail in CBC output.", -1, 999999,
      getLogLevel(),
      "If set to 0 then there should be no output in normal circumstances. A "
      "value of 1 is probably the best value for most uses, while 2 and 3 give "
      "more information.");

  parameters_[LPLOGLEVEL].setup(
      "log!Level", "Level of detail in LP solver output.", -1, 999999,
      getLpLogLevel(),
      "If set to 0 then there should be no output in normal circumstances. A "
      "value of 1 is probably the best value for most uses, while 2 and 3 give "
      "more information.");

  parameters_[MAXSAVEDSOLS].setup(
      "maxSaved!Solutions", "Maximum number of solutions to save", 0,
      COIN_INT_MAX, 1, "Number of solutions to save.");

  parameters_[MAXSLOWCUTS].setup(
      "slow!cutpasses", "Maximum number of rounds for slower cut generators",
      -1, COIN_INT_MAX, 10,
      "Some cut generators are fairly slow - this limits the number of times "
      "they are tried. The cut generators identified as 'may be slow' at "
      "present are Lift and project cuts and both versions of Reduce and Split "
      "cuts.");

  parameters_[MOREMOREMIPOPTIONS].setup(
      "more2!MipOptions", "More more dubious options for mip", -1, COIN_INT_MAX,
      0, "", CoinParam::displayPriorityNone);

  parameters_[MULTIPLEROOTS].setup(
      "multiple!RootPasses",
      "Do multiple root passes to collect cuts and solutions", 0, COIN_INT_MAX,
      0,
      "Solve (in parallel, if enabled) the root phase this number of times, "
      "each with its own different seed, and collect all solutions and cuts "
      "generated. The actual format is aabbcc where aa is the number of extra "
      "passes; if bb is non zero, then it is number of threads to use "
      "(otherwise uses threads setting); and cc is the number of times to do "
      "root phase. The solvers do not interact with each other.  However if "
      "extra passes are specified then cuts are collected and used in later "
      "passes - so there is interaction there. Some parts of this "
      "implementation have their origin in idea of Andrea Lodi, Matteo "
      "Fischetti, Michele Monaci, Domenico Salvagnin, and Andrea Tramontani.",
      CoinParam::displayPriorityNone);

  parameters_[ODDWEXTMETHOD].setup(
      "oddwext!method",
      "Strategy used to search for wheel centers for the cuts found by Odd "
      "Wheel Cut Separation routine",
      0, 2, 2,
      "Sets the method used in the extension module of Odd Wheel Cut "
      "Separation routine: 0=no extension; 1=one variable; 2=clique");

  parameters_[OUTPUTFORMAT].setup(
      "output!Format", "Which output format to use", 1, 6, 2,
      "Normally export will be done using normal representation for numbers "
      "and two values per line.  You may want to do just one per line (for "
      "grep or suchlike) and you may wish to save with absolute accuracy using "
      "a coded version of the IEEE value. A value of 2 is normal. Otherwise, "
      "odd values give one value per line, even values two.  Values of 1 and 2 "
      "give normal format, 3 and 4 give greater precision, 5 and 6 give IEEE "
      "values.  When exporting a basis, 1 does not save values, 2 saves "
      "values, 3 saves with greater accuracy and 4 saves in IEEE format.");

  parameters_[PRINTOPTIONS].setup(
      "pO!ptions", "Dubious print options", 0, COIN_INT_MAX, 0,
      "If this is greater than 0 then presolve will give more information and "
      "branch and cut will give statistics");

  parameters_[PROCESSTUNE].setup(
      "tune!PreProcess", "Dubious tuning parameters for preprocessing", 0,
      COIN_INT_MAX, 0,
      "Format aabbcccc - \n If aa then this is number of major passes (i.e. "
      "with presolve) \n If bb and bb>0 then this is number of minor passes "
      "(if unset or 0 then 10) \n cccc is bit set \n 0 - 1 Heavy probing \n 1 "
      "- 2 Make variables integer if possible (if obj value)\n 2 - 4 As above "
      "but even if zero objective value\n 7 - 128 Try and create cliques\n 8 - "
      "256 If all +1 try hard for dominated rows\n 9 - 512 Even heavier "
      "probing \n 10 - 1024 Use a larger feasibility tolerance in presolve\n "
      "11 - 2048 Try probing before creating cliques\n 12 - 4096 Switch off "
      "duplicate column checking for integers \n 13 - 8192 Allow scaled "
      "duplicate column checking \n \n     Now aa 99 has special meaning i.e. "
      "just one simple presolve.",
      CoinParam::displayPriorityLow);

  parameters_[RANDOMSEED].setup(
      "randomC!bcSeed", "Random seed for Cbc", -1, COIN_INT_MAX, -1,
      "Allows initialization of the random seed for pseudo-random numbers used "
      "in heuristics such as the Feasibility Pump to decide whether to round "
      "up or down. The special value of 0 lets Cbc use the time of the day for "
      "the initial seed.");

  parameters_[STRONGSTRATEGY].setup(
      "expensive!Strong", "Whether to do even more strong branching", 0,
      COIN_INT_MAX, 0,
      "Strategy for extra strong branching. 0 is normal strong branching. 1, "
      "2, 4, and 6 does strong branching on all fractional variables if at the "
      "root node (1), at depth less than modifier (2), objective equals best "
      "possible (4), or at depth less than modifier and objective equals best "
      "possible (6). 11, 12, 14, and 16 are like 1, 2, 4, and 6, "
      "respecitively, but do strong branching on all integer (incl. "
      "non-fractional) variables. Values >= 100 are used to specify a depth "
      "limit (value/100), otherwise 5 is used. If the values >= 100, then "
      "above rules are applied to value%100.",
      CoinParam::displayPriorityNone);

  parameters_[TESTOSI].setup("testO!si", "Test OsiObject stuff",
                                            -1, COIN_INT_MAX, -1, "",
                                            CoinParam::displayPriorityNone);

#ifdef CBC_THREAD
  parameters_[THREADS].setup(
      "thread!s", "Number of threads to try and use", -100, 100000, 0,
      "To use multiple threads, set threads to number wanted.  It may be "
      "better to use one or two more than number of cpus available.  If 100+n "
      "then n threads and search is repeatable (maybe be somewhat slower), if "
      "200+n use threads for root cuts, 400+n threads used in sub-trees.",
      CoinParam::displayPriorityLow);
#endif

  parameters_[USERCBC].setup(
      "userCbc", "Hand coded Cbc stuff", 0, COIN_INT_MAX, 0,
      "There are times (e.g., when using AMPL interface) when you may wish to "
      "do something unusual.  Look for USERCBC in main driver and modify "
      "sample code.",
      CoinParam::displayPriorityNone);

  parameters_[VERBOSE].setup(
      "verbose", "Switches on longer help on single ?", 0, 15,
      verbose_,
      "Set to 1 to get short help with ? list, 2 to get long help.",
      CoinParam::displayPriorityNone);

  parameters_[VUBTRY].setup(
      "vub!heuristic", "Type of VUB heuristic", -2, 20, -1,
      "This heuristic tries to fix some integer variables.",
      CoinParam::displayPriorityNone);
}

//###########################################################################
//###########################################################################

void CbcSettings::addCbcSolverBoolParams() {
  for (int code = CBC_FIRSTBOOLPARAM + 1;
       code < CBC_LASTBOOLPARAM; code++) {
    parameters_[code].appendKwd("off", CbcSettings::ParamOff);
    parameters_[code].appendKwd("on", CbcSettings::ParamOn);
    parameters_[code].setPushFunc(CbcParamUtils::pushCbcSolverKwdParam);
  }

  parameters_[CPX].setup(
      "cplex!Use", "Whether to use Cplex!", "off", 0,
      "If the user has Cplex, but wants to use some of Cbc's heuristics then "
      "you can!  If this is on, then Cbc will get to the root node and then "
      "hand over to Cplex.  If heuristics find a solution this can be "
      "significantly quicker.  You will probably want to switch off Cbc's cuts "
      "as Cplex thinks they are genuine constraints.  It is also probable that "
      "you want to switch off preprocessing, although for difficult problems "
      "it is worth trying both.");

  parameters_[DOHEURISTIC].setup(
      "doH!euristic", "Do heuristics before any preprocessing", "off",
      getDiveCoefficientMode(),
      "Normally heuristics are done in branch and bound.  It may be useful to "
      "do them outside. Only those heuristics with 'both' or 'before' set will "
      "run. Doing this may also set cutoff, which can help with "
      "preprocessing.");

  parameters_[ERRORSALLOWED].setup(
      "error!sAllowed", "Whether to allow import errors", "off", 0,
      "The default is not to use any model which had errors when reading the "
      "mps file.  Setting this to 'on' will allow all errors from which the "
      "code can recover simply by ignoring the error.  There are some errors "
      "from which the code can not recover, e.g., no ENDATA.  This has to be "
      "set before import, i.e., -errorsAllowed on -import xxxxxx.mps.");

  parameters_[EXTRAVARIABLES].setup(
      "extraV!ariables", "Allow creation of extra integer variables", "off", 0,
      "Switches on a trivial re-formulation that introduces extra integer "
      "variables to group together variables with same cost.",
      CoinParam::displayPriorityLow);

  parameters_[MESSAGES].setup(
      "mess!ages", "Controls whether standardised message prefix is printed",
      "off", 0,
      "By default, messages have a standard prefix, such as:\n   Clp0005 2261  "
      "Objective 109.024 Primal infeas 944413 (758)\nbut this program turns "
      "this off to make it look more friendly.  It can be useful to turn them "
      "back on if you want to be able to 'grep' for particular messages or if "
      "you intend to override the behavior of a particular message.");

  parameters_[PREPROCNAMES].setup(
      "PrepN!ames", "If column names will be kept in pre-processed model",
      "off", 0,
      "Normally the preprocessed model has column names replaced by new names "
      "C0000... Setting this option to on keeps original names in variables "
      "which still exist in the preprocessed problem");

  parameters_[SOS].setup(
      "sos!Options", "Whether to use SOS from AMPL", "off", 0,
      "Normally if AMPL says there are SOS variables they should be used, but "
      "sometimes they should be turned off - this does so.");

  parameters_[USESOLUTION].setup(
      "force!Solution", "Whether to use given solution as crash for BAB", "off",
      0,
      "If on then tries to branch to solution given by AMPL or priorities "
      "file.");
}

//###########################################################################
//###########################################################################

void CbcSettings::addCbcSolverCutParams() {
  for (int code = CBC_FIRSTCUTPARAM + 1;
       code < CBC_LASTCUTPARAM; code++) {
    parameters_[code].setPushFunc(CbcParamUtils::pushCbcSolverKwdParam);
  }

  parameters_[CLIQUECUTS].setup(
      "clique!Cuts", "Whether to use clique cuts", "off", CbcSettings::CGOff,
      "This switches on clique cuts (either at root or in entire tree). See "
      "branchAndCut for information on options.");

  parameters_[CUTSTRATEGY].setup(
      "cuts!OnOff", "Switches all cuts on or off", "off", CbcSettings::CGOff,
      "This can be used to switch on or off all cuts (apart from Reduce and "
      "Split).  Then you can set individual ones off or on.  See branchAndCut "
      "for information on options.");

  parameters_[FLOWCUTS].setup(
      "flow!CoverCuts", "Whether to use Flow Cover cuts", "off",
      getFlowMode(),
      "This switches on flow cover cuts (either at root or in entire "
      "tree)." CUTS_LONGHELP);

  parameters_[GMICUTS].setup(
      "GMI!Cuts", "Whether to use alternative Gomory cuts", "off",
      getGMIMode(),
      CUTS_LONGHELP " This version is by Giacomo Nannicini and may be more "
                    "robust than gomoryCuts.");

  parameters_[GOMORYCUTS].setup(
      "gomory!Cuts", "Whether to use Gomory cuts", "off",
      getGomoryMode(),
      "The original cuts - beware of imitations!  Having gone out of favor, "
      "they are now more fashionable as LP solvers are more robust and they "
      "interact well with other cuts.  They will almost always give cuts "
      "(although in this executable they are limited as to number of variables "
      "in cut).  However the cuts may be dense so it is worth experimenting "
      "(Long allows any length). " CUTS_LONGHELP
      " Reference: https://github.com/coin-or/Cgl/wiki/CglGomory");

  parameters_[KNAPSACKCUTS].setup(
      "knapsack!Cuts", "Whether to use Knapsack cuts", "off",
      getKnapsackMode(),
      "This switches on knapsack cuts (either at root or in entire tree). See "
      "branchAndCut for information on options.");

  parameters_[LANDPCUTS].setup(
      "lift!AndProjectCuts", "Whether to use lift-and-project cuts", "off",
      getLandPMode(),
      "This switches on lift-and-project cuts (either at root or in entire "
      "tree). See branchAndCut for information on options.");

  parameters_[LAGOMORYCUTS].setup(
      "lagomory!Cuts", "Whether to use Lagrangean Gomory cuts", "off",
      getLaGomoryMode(),
      "This is a gross simplification of 'A Relax-and-Cut Framework for "
      "Gomory's Mixed-Integer Cuts' by Matteo Fischetti & Domenico Salvagnin.  "
      "This simplification just uses original constraints while modifying "
      "objective using other cuts. So you don't use messy constraints "
      "generated by Gomory etc. A variant is to allow non messy cuts e.g. "
      "clique cuts. So 'only' does this while 'clean' also allows integral "
      "valued cuts.  'End' is recommended and waits until other cuts have "
      "finished before it does a few passes. The length options for gomory "
      "cuts are used.");

  parameters_[LATWOMIRCUTS].setup(
      "latwomir!Cuts", "Whether to use Lagrangean Twomir cuts", "off",
      getLaTwomirMode(),
      "This is a Lagrangean relaxation for Twomir cuts.  See lagomoryCuts for "
      "description of options.");

  parameters_[MIRCUTS].setup(
      "mixed!IntegerRoundingCuts", "Whether to use Mixed Integer Rounding cuts",
      "off", getMirMode(),
      "This switches on mixed integer rounding cuts (either at root or in "
      "entire tree).  See branchAndCut for information on options.");

  parameters_[ODDWHEELCUTS].setup(
      "oddwheel!Cuts", "Whether to use odd wheel cuts", "off",
      getOddWheelMode(),
      "This switches on odd-wheel inequalities (either at root or in entire "
      "tree).");

  parameters_[PROBINGCUTS].setup(
      "probing!Cuts", "Whether to use Probing cuts", "off",
      getProbingMode(),
      "This switches on probing cuts (either at root or in entire tree). See "
      "branchAndCut for information on options.");

  parameters_[REDSPLITCUTS].setup(
      "reduce!AndSplitCuts", "Whether to use Reduce-and-Split cuts", "off",
      getRedSplitMode(),
      "This switches on reduce and split cuts (either at root or in entire "
      "tree). See branchAndCut for information on options.");

  parameters_[REDSPLIT2CUTS].setup(
      "reduce2!AndSplitCuts", "Whether to use Reduce-and-Split cuts - style 2",
      "off", getRedSplit2Mode(),
      "This switches on reduce and split cuts (either at root or in entire "
      "tree). See branchAndCut for information on options.");

  parameters_[RESIDCAPCUTS].setup(
      "residual!CapacityCuts", "Whether to use Residual Capacity cuts", "off",
      getResidCapMode(),
      CUTS_LONGHELP
      " Reference: https://github.com/coin-or/Cgl/wiki/CglResidualCapacity");

  parameters_[TWOMIRCUTS].setup(
      "two!MirCuts", "Whether to use Two phase Mixed Integer Rounding cuts",
      "off", getTwomirMode(),
      "This switches on two phase mixed integer rounding cuts (either at root "
      "or in entire tree). See branchAndCut for information on options.");

  parameters_[ZEROHALFCUTS].setup(
      "zero!HalfCuts", "Whether to use zero half cuts", "off",
      getZeroHalfMode(),
      CUTS_LONGHELP " This implementation was written by Alberto Caprara.");

  // Populate the keyword lists
  for (int code = CBC_FIRSTCUTPARAM + 1;
       code < CBC_LASTCUTPARAM; code++) {
    // First the common keywords
    switch (code) {
    case CUTSTRATEGY:
    case CLIQUECUTS:
    case FLOWCUTS:
    case GMICUTS:
    case GOMORYCUTS:
    case KNAPSACKCUTS:
    case LANDPCUTS:
    case MIRCUTS:
    case ODDWHEELCUTS:
    case PROBINGCUTS:
    case RESIDCAPCUTS:
    case TWOMIRCUTS:
    case ZEROHALFCUTS: {
      parameters_[code].appendKwd("on", CbcSettings::CGOn);
      parameters_[code].appendKwd("root", CbcSettings::CGRoot);
      parameters_[code].appendKwd("ifmove", CbcSettings::CGIfMove);
      parameters_[code].appendKwd("forceon", CbcSettings::CGForceOn);
      break;
    }
    default: { break; }

      // Now, add some additional keywords for different classes
      switch (code) {
      case GOMORYCUTS: {
        parameters_[code].appendKwd("onglobal", CbcSettings::CGOnGlobal);
        parameters_[code].appendKwd("forceandglobal", CbcSettings::CGForceAndGlobal);
        parameters_[code].appendKwd("forcelongon", CbcSettings::CGForceLongOn);
        parameters_[code].appendKwd("longer", CbcSettings::CGLonger);
        parameters_[code].appendKwd("shorter", CbcSettings::CGShorter);
        break;
      }
      case GMICUTS: {
        parameters_[code].appendKwd("long", CbcSettings::CGLong);
        parameters_[code].appendKwd("longroot", CbcSettings::CGLongRoot);
        parameters_[code].appendKwd("longifmove", CbcSettings::CGLongIfMove);
        parameters_[code].appendKwd("forcelongon", CbcSettings::CGForceLongOn);
        parameters_[code].appendKwd("longendonly", CbcSettings::CGLongEndOnly);
        break;
      }
      case LAGOMORYCUTS: {
        parameters_[code].appendKwd("root", CbcSettings::CGRoot);
        parameters_[code].appendKwd("onlyaswellroot", CbcSettings::CGOnlyAsWellRoot);
        parameters_[code].appendKwd("cleanaswellroot", CbcSettings::CGCleanAsWellRoot);
        parameters_[code].appendKwd("bothaswellroot", CbcSettings::CGCleanBothAsWellRoot);
        // Here, we intentionally drop through to the next set
      }
      case LATWOMIRCUTS: {
        parameters_[code].appendKwd("endonlyroot", CbcSettings::CGEndOnlyRoot);
        parameters_[code].appendKwd("endcleanroot", CbcSettings::CGEndCleanRoot);
        parameters_[code].appendKwd("endonly", CbcSettings::CGEndOnly);
        parameters_[code].appendKwd("endclean", CbcSettings::CGEndClean);
        parameters_[code].appendKwd("endboth", CbcSettings::CGEndBoth);
        parameters_[code].appendKwd("onlyaswell", CbcSettings::CGOnlyAsWell);
        parameters_[code].appendKwd("cleanaswell", CbcSettings::CGCleanAsWell);
        parameters_[code].appendKwd("bothaswell", CbcSettings::CGBothAsWell);
        parameters_[code].appendKwd("onlyinstead", CbcSettings::CGOnlyInstead);
        parameters_[code].appendKwd("cleaninstead", CbcSettings::CGCleanInstead);
        parameters_[code].appendKwd("bothinstead", CbcSettings::CGBothInstead);
        break;
      }
      case ODDWHEELCUTS: {
        parameters_[code].appendKwd("onglobal", CbcSettings::CGOnGlobal);
        break;
      }
      case PROBINGCUTS: {
        parameters_[code].appendKwd("forceonbut", CbcSettings::CGForceOnBut);
        break;
      }
      case REDSPLITCUTS:
      case REDSPLIT2CUTS: {
        parameters_[code].appendKwd("on", CbcSettings::CGOn);
        parameters_[code].appendKwd("root", CbcSettings::CGRoot);
        parameters_[code].appendKwd("longon", CbcSettings::CGLongOn);
        parameters_[code].appendKwd("longroot", CbcSettings::CGLongRoot);
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

void CbcSettings::addCbcSolverHeurParams() {
  for (int code = CBC_FIRSTHEURPARAM + 1;
       code < CBC_LASTHEURPARAM; code++) {
    parameters_[code].setPushFunc(CbcParamUtils::pushCbcSolverKwdParam);
  }

  parameters_[COMBINE].setup(
      "combine!Solutions", "Whether to use combine solution heuristic", "off",
      getCombineMode(),
      "This switches on a heuristic which does branch and cut on the problem "
      "given by just using variables which have appeared in one or more "
      "solutions. It is obviously only tried after two or more solutions.");

  parameters_[CROSSOVER].setup(
      "combine2!Solutions", "Whether to use crossover solution heuristic",
      "off", getCrossoverMode());

  parameters_[DINS].setup(
      "Dins", "Whether to try Distance Induced Neighborhood Search", "off",
      getDinsMode(), HEURISTICS_LONGHELP);

  parameters_[DIVINGC].setup(
      "DivingC!oefficient", "Whether to try Coefficient diving heuristic",
      "off", getDiveCoefficientMode(), HEURISTICS_LONGHELP);

  parameters_[DIVINGF].setup(
      "DivingF!ractional", "Whether to try Fractional diving heuristic", "off",
      getDiveFractionalMode(), HEURISTICS_LONGHELP);

  parameters_[DIVINGG].setup(
      "DivingG!uided", "Whether to try Guided diving heuristic", "off",
      getDiveGuidedMode(), HEURISTICS_LONGHELP);

  parameters_[DIVINGL].setup(
      "DivingL!ineSearch", "Whether to try Linesearch diving heuristic", "off",
      getDiveLineSearchMode(), HEURISTICS_LONGHELP);

  parameters_[DIVINGP].setup(
      "DivingP!seudocost", "Whether to try Pseudocost diving heuristic", "off",
      getDivePseudocostMode(), HEURISTICS_LONGHELP);

  parameters_[DIVINGS].setup(
      "DivingS!ome", "Whether to try Diving heuristics", "off",
      getDiveRandomMode(),
      "This switches on a random diving heuristic at various times. One may "
      "prefer to individually turn diving heuristics on or off. ");

  parameters_[DIVINGV].setup(
      "DivingV!ectorLength", "Whether to try Vectorlength diving heuristic",
      "off", getDiveVectorLengthMode(), HEURISTICS_LONGHELP);

  parameters_[DW].setup(
      "dw!Heuristic", "Whether to try Dantzig Wolfe heuristic", "off",
      getDWMode(),
      "This heuristic is very very compute intensive. It tries to find a "
      "Dantzig Wolfe structure and use that. " HEURISTICS_LONGHELP);

  parameters_[FPUMP].setup(
      "feas!ibilityPump", "Whether to try Feasibility Pump", "off",
      getFeasPumpMode(),
      "This switches on feasibility pump heuristic at root. This is due to "
      "Fischetti and Lodi and uses a sequence of LPs to try and get an integer "
      "feasible solution.  Some fine tuning is available by "
      "passFeasibilityPump.");

  parameters_[GREEDY].setup(
      "greedy!Heuristic", "Whether to use a greedy heuristic", "off",
      getGreedyCoverMode(),
      "Switches on a pair of greedy heuristic which will try and obtain a "
      "solution.  It may just fix a percentage of variables and then try a "
      "small branch and cut run.");

  parameters_[HEURISTICSTRATEGY].setup(
      "heur!isticsOnOff", "Switches most heuristics on or off", "off", 0,
      "This can be used to switch on or off all heuristics.  Then you can set "
      "individual ones off or on.  CbcTreeLocal is not included as it "
      "dramatically alters search.");

  parameters_[NAIVE].setup(
      "naive!Heuristics", "Whether to try some stupid heuristic", "off",
      getNaiveHeurMode(),
      "This is naive heuristics which, e.g., fix all integers with costs to "
      "zero!. " HEURISTICS_LONGHELP,
      CoinParam::displayPriorityLow);

  parameters_[PIVOTANDFIX].setup(
      "pivotAndF!ix", "Whether to try Pivot and Fix heuristic", "off",
      getPivotAndFixMode(), HEURISTICS_LONGHELP);

#if 0
  parameters_[PIVOTANDCOMPLEMENT].setup("pivotAndC!omplement",
                                                        "Whether to try Pivot and Complement heuristic",
                                                        "off", getPivotAndComplementMode(),
                                                        HEURISTICS_LONGHELP);
#endif

  parameters_[PROXIMITY].setup(
      "proximity!Search", "Whether to do proximity search heuristic", "off",
      getProximityMode(),
      "This heuristic looks for a solution close to the incumbent solution "
      "(Fischetti and Monaci, 2012). The idea is to define a sub-MIP without "
      "additional constraints but with a modified objective function intended "
      "to attract the search in the proximity of the incumbent. The approach "
      "works well for 0-1 MIPs whose solution landscape is not too irregular "
      "(meaning the there is reasonable probability of finding an improved "
      "solution by flipping a small number of binary variables), in particular "
      "when it is applied to the first heuristic solutions found at the root "
      "node. " HEURISTICS_LONGHELP); // Can also set different maxNode
                                     // cbcSettings by plusnnnn (and are
                                     // 'on'(on==30)).

  parameters_[RANDROUND].setup(
      "randomi!zedRounding", "Whether to try randomized rounding heuristic",
      "off", getRandRoundMode(), HEURISTICS_LONGHELP);

  parameters_[RENS].setup(
      "Rens", "Whether to try Relaxation Enforced Neighborhood Search", "off",
      getRensMode(),
      HEURISTICS_LONGHELP " Value 'on' just does 50 nodes. 200, 1000, and "
                          "10000 does that many nodes.");

  parameters_[RINS].setup(
      "Rins", "Whether to try Relaxed Induced Neighborhood Search", "off",
      getRinsMode(), HEURISTICS_LONGHELP);

  parameters_[ROUNDING].setup(
      "round!ingHeuristic", "Whether to use Rounding heuristic", "off",
      getRoundingMode(),
      "This switches on a simple (but effective) rounding heuristic at each "
      "node of tree.");

  parameters_[LOCALTREE].setup(
      "local!TreeSearch", "Whether to use local tree search", "off",
      getLocalTreeMode(),
      "This switches on a local search algorithm when a solution is found.  "
      "This is from Fischetti and Lodi and is not really a heuristic although "
      "it can be used as one. When used from this program it has limited "
      "functionality.  It is not controlled by heuristicsOnOff.");

  parameters_[VND].setup(
      "Vnd!VariableNeighborhoodSearch",
      "Whether to try Variable Neighborhood Search", "off",
       getVndMode(), HEURISTICS_LONGHELP);

  // Populate the keyword lists
  for (int code = CBC_FIRSTHEURPARAM + 1;
       code < CBC_LASTHEURPARAM; code++) {
    // First the common keywords
    switch (code) {
    case HEURISTICSTRATEGY:
    case CROSSOVER:
    case DINS:
    case DIVINGC:
    case DIVINGF:
    case DIVINGG:
    case DIVINGL:
    case DIVINGP:
    case DIVINGS:
    case DIVINGV:
    case DW:
    case NAIVE:
    case PIVOTANDFIX:
    case PIVOTANDCOMPLEMENT:
    case PROXIMITY:
    case RANDROUND:
    case RENS:
    case RINS:
    case ROUNDING:
    case VND: {
      parameters_[code].appendKwd("off", CbcSettings::CGOff);
      parameters_[code].appendKwd("on", CbcSettings::CGOn);
      parameters_[code].appendKwd("both", CbcSettings::CGRoot);
      parameters_[code].appendKwd("before", CbcSettings::CGIfMove);
      break;
    }
    default:
      break;
    }
    // Check the unique keywords
    switch (code) {
    case COMBINE: {
      parameters_[code].appendKwd("off", CbcSettings::HeurOff);
      parameters_[code].appendKwd("on", CbcSettings::HeurOn);
      break;
    }
    case DINS: {
      parameters_[code].appendKwd("often", CbcSettings::HeurOften);
      break;
    }
    case FPUMP: {
      parameters_[code].appendKwd("off", CbcSettings::HeurOff);
      parameters_[code].appendKwd("on", CbcSettings::HeurOn);
      break;
    }
    case GREEDY: {
      parameters_[code].appendKwd("off", CbcSettings::HeurOff);
      parameters_[code].appendKwd("on", CbcSettings::HeurOn);
      parameters_[code].appendKwd("root", CbcSettings::HeurRoot);
      break;
    }
    case LOCALTREE: {
      parameters_[code].appendKwd("off", CbcSettings::HeurOff);
      parameters_[code].appendKwd("on", CbcSettings::HeurOn);
    }
    case PROXIMITY: {
      parameters_[code].appendKwd("10", CbcSettings::HeurTen);
      parameters_[code].appendKwd("100", CbcSettings::HeurOneHundred);
      parameters_[code].appendKwd("300", CbcSettings::HeurThreeHundred);
      break;
    }
    case RENS: {
      parameters_[code].appendKwd("200", CbcSettings::HeurTwoHundred);
      parameters_[code].appendKwd("1000", CbcSettings::HeurOneThousand);
      parameters_[code].appendKwd("10000", CbcSettings::HeurTenThousand);
      parameters_[code].appendKwd("dj", CbcSettings::HeurDj);
      parameters_[code].appendKwd("djbefore", CbcSettings::HeurDjBefore);
      parameters_[code].appendKwd("usesolution", CbcSettings::HeurUseSolution);
      break;
    }
    case RINS: {
      parameters_[code].appendKwd("often", CbcSettings::HeurOften);
      break;
    }
    case VND: {
      parameters_[code].appendKwd("intree", CbcSettings::HeurInTree);
      break;
    }
    default:
      break;
    }
  }
}

//###########################################################################
//###########################################################################

/* Function to set up cbc (CbcModel) parameters.  */

void CbcSettings::addCbcModelParams()
{

  parameters_[ALLOWABLEGAP].setup(
      "allow!ableGap",
      "Stop when gap between best possible and incumbent is less than this",
      0.0, 1.0e20, 0.0,
      "If the gap between best solution and best possible solution is less "
      "than this then the search will be terminated. Also see ratioGap.");

  parameters_[CUTOFF].setup(
      "cuto!ff", "All solutions must be better than this", -1.0e60, 1.0e60,
      1.0e50,
      "All solutions must be better than this value (in a minimization sense). "
      " This is also set by cbc whenever it obtains a solution and is set to "
      "the value of the objective for the solution minus the cutoff "
      "increment.");

  parameters_[DIRECTION].setup(
      "direction", "Minimize or maximize", "min!imize",
      CbcSettings::OptDirMinimize,
      "The default is minimize - use 'direction maximize' for "
      "maximization.\nYou can also use the parameters_ 'maximize' or "
      "'minimize'.");
  parameters_[DIRECTION].appendKwd("max!imize", CbcSettings::OptDirMaximize);
  parameters_[DIRECTION].appendKwd("zero", CbcSettings::OptDirZero);

  parameters_[INCREMENT].setup(
      "inc!rement",
      "A new solution must be at least this much better than the incumbent",
      -1.0e20, 1.0e20, model_->getDblParam(CbcModel::CbcCutoffIncrement),
      "Whenever a solution is found the bound on future solutions is set to "
      "the objective of the solution (in a minimization sense) plus the "
      "specified increment.  If this option is not specified, the code will "
      "try and work out an increment.  E.g., if all objective coefficients are "
      "multiples of 0.01 and only integer variables have entries in objective "
      "then the increment can be set to 0.01.  Be careful if you set this "
      "negative!");

  parameters_[INFEASIBILITYWEIGHT].setup(
      "inf!easibilityWeight",
      "Each integer infeasibility is expected to cost this much", 0.0, 1.0e20,
      model_->getDblParam(CbcModel::CbcInfeasibilityWeight),
      "A primitive way of deciding which node to explore next.  Satisfying "
      "each integer infeasibility is expected to cost this much.");

  parameters_[INTEGERTOLERANCE].setup(
      "integerT!olerance",
      "For an optimal solution, no integer variable may be farther than this "
      "from an integer value",
      1.0e-20, 0.5, model_->getDblParam(CbcModel::CbcIntegerTolerance),
      "When checking a solution for feasibility, if the difference between the "
      "value of a variable and the nearest integer is less than the integer "
      "tolerance, the value is considered to be integral. Beware of setting "
      "this smaller than the primal tolerance.");

  parameters_[LOGLEVEL].setup(
      "bclog!Level", "Level of detail in Coin branch and Cut output", -1, 63,
      model_->messageHandler()->logLevel(),
      "If set to 0 then there should be no output in normal circumstances. A "
      "value of 1 is probably the best value for most uses, while 2 and 3 give "
      "more information.");

  parameters_[MAXIMIZE].setup(
      "max!imize", "Set optimization direction to maximize",
      "The default is minimize - use 'maximize' for maximization.\n A synonym "
      "for 'direction maximize'.",
      CoinParam::displayPriorityHigh);

  parameters_[MAXNODES].setup(
      "maxN!odes", "Maximum number of nodes to evaluate", 1, 2147483647, 1,
      "This is a repeatable way to limit search.  Normally using time is "
      "easier but then the results may not be repeatable.");

  parameters_[MAXNODESNOTIMPROVING].setup(
      "maxNNI!FS",
      "Maximum number of nodes to be processed without improving the incumbent "
      "solution.",
      -1, COIN_INT_MAX, -1,
      "This criterion specifies that when a feasible solution is available, "
      "the search should continue only if better feasible solutions were "
      "produced in the last nodes.");

  parameters_[MAXSECONDSNOTIMPROVING].setup(
      "secni!fs", "maximum seconds without improving the incumbent solution",
      -1.0, COIN_DBL_MAX, -1.0,
      "With this stopping criterion, after a feasible solution is found, the "
      "search should continue only if the incumbent solution was updated "
      "recently, the tolerance is specified here. A discussion on why this "
      "criterion can be useful is included here: "
      "https://yetanothermathprogrammingconsultant.blogspot.com/2019/11/"
      "mip-solver-stopping-criteria.html .");

  parameters_[MAXSOLS].setup(
      "maxSo!lutions", "Maximum number of feasible solutions to get", 1,
      COIN_INT_MAX, COIN_INT_MAX,
      "You may want to stop after (say) two solutions or an hour. This is "
      "checked every node in tree, so it is possible to get more solutions "
      "from heuristics.");

  parameters_[MINIMIZE].setup(
      "min!imize", "Set optimization direction to minimize",
      "The default is minimize - use 'maximize' for maximization.\nThis should "
      "only be necessary if you have previously set maximization. A synonym "
      "for 'direction minimize'.");

  parameters_[MIPOPTIONS].setup(
      "mipO!ptions", "Dubious options for mip", 0, COIN_INT_MAX, 0, "",
      CoinParam::displayPriorityNone);

  parameters_[MOREMIPOPTIONS].setup(
      "more!MipOptions", "More dubious options for mip", -1, COIN_INT_MAX, 0,
      "", CoinParam::displayPriorityNone);

#if 0
  parameters_[NUMBERMINI].setup("miniT!ree",
      "Size of fast mini tree", 0, COIN_INT_MAX, 0,
      "The idea is that I can do a small tree fast. This is a first try and will" 
      "hopefully become more sophisticated.", CoinParam::displayPriorityNone);
#endif

  parameters_[NUMBERANALYZE].setup(
      "numberA!nalyze", "Number of analysis iterations", -COIN_INT_MAX,
      COIN_INT_MAX, 0,
      "This says how many iterations to spend at the root node analyzing the "
      "problem.  This is a first try and will hopefully become more "
      "sophisticated.",
      CoinParam::displayPriorityNone);

  parameters_[REVERSE].setup(
      "reverse", "Reverses sign of objective",
      "Useful for testing if maximization works correctly",
      CoinParam::displayPriorityNone);

  parameters_[CUTPASS].setup(
      "passC!uts", "Number of cut passes at root node", -999999, 999999,
      model_->getMaximumCutPassesAtRoot(),
      "The default is 100 passes if less than 500 columns, 100 passes (but "
      "stop if the drop is small) if less than 5000 columns, 20 otherwise.");

  parameters_[GAPRATIO].setup(
      "ratio!Gap",
      "Stop when the gap between the best possible solution and the incumbent "
      "is less than this fraction of the larger of the two",
      0.0, 1.0e20, model_->getDblParam(CbcModel::CbcAllowableFractionGap),
      "If the gap between the best solution and the best possible solution is "
      "less than this fraction of the objective value at the root node then "
      "the search will terminate.  See 'allowableGap' for a way of using "
      "absolute value rather than fraction.");

  parameters_[TIMELIMIT].setup(
      "sec!onds", "Maximum seconds for branch and cut", -1.0, 1.0e12, -1.0,
      "After this many seconds the program will act as if maximum nodes had "
      "been reached.");

  parameters_[STRONGBRANCHING].setup(
      "strong!Branching", "Number of variables to look at in strong branching",
      0, 999999, model_->numberStrong(),
      "In order to decide which variable to branch on, the code will choose up "
      "to this number of unsatisfied variables and try mini up and down "
      "branches.  The most effective one is chosen. If a variable is branched "
      "on many times then the previous average up and down costs may be used - "
      "see number before trust.");

  parameters_[NUMBERBEFORE].setup(
      "trust!PseudoCosts", "Number of branches before we trust pseudocosts", -1,
      2000000, model_->numberBeforeTrust(),
      "Using strong branching computes pseudo-costs.  After this many times "
      "for a variable we just trust the pseudo costs and do not do any more "
      "strong branching.");

  for (int code = CBC_FIRSTMODELPARAM + 1;
       code < CBC_LASTMODELPARAM; code++) {
     if (parameters_[code].type() == CoinParam::coinParamInt){
        parameters_[code].setPushFunc(CbcParamUtils::pushCbcModelIntParam);
     }else{
        parameters_[code].setPushFunc(CbcParamUtils::pushCbcModelDblParam);
     }
  }
}

//###########################################################################
//###########################################################################

/*
  Access functions for cut generators. These support lazy
  creation --- if mode_ is other than CGOff, an object is created if
  necessary and a pointer is stored in proto_. The pointer is returned as
  a generic CglCutGenerator or CbcHeuristic. The return value of the function
  is the value of mode_.

  Because the model may have changed, the default for heuristics is to delete
  any existing object and create a new one. This can be suppressed if desired.
*/

CbcSettings::CGMode CbcSettings::getClique(CglCutGenerator *&gen) {
  if (clique_.mode_ != CbcSettings::CGOff && clique_.proto_ == 0) {
    clique_.proto_ = new CglClique();
    clique_.proto_->setStarCliqueReport(clique_.starCliqueReport_);
    clique_.proto_->setRowCliqueReport(clique_.rowCliqueReport_);
    clique_.proto_->setMinViolation(clique_.minViolation_);
  }
  gen = dynamic_cast<CglCutGenerator *>(clique_.proto_);

  return (clique_.mode_);
}

CbcSettings::CGMode CbcSettings::getFlow(CglCutGenerator *&gen)

{
  if (flow_.mode_ != CbcSettings::CGOff && flow_.proto_ == 0) {
    flow_.proto_ = new CglFlowCover();
  }
  gen = dynamic_cast<CglCutGenerator *>(flow_.proto_);

  return (flow_.mode_);
}

CbcSettings::CGMode CbcSettings::getGomory(CglCutGenerator *&gen) {
  if (gomory_.mode_ != CbcSettings::CGOff && gomory_.proto_ == 0) {
    gomory_.proto_ = new CglGomory();
    gomory_.proto_->setLimitAtRoot(gomory_.limitAtRoot_);
    gomory_.proto_->setLimit(gomory_.limit_);
  }
  gen = dynamic_cast<CglCutGenerator *>(gomory_.proto_);

  return (gomory_.mode_);
}

CbcSettings::CGMode CbcSettings::getKnapsack(CglCutGenerator *&gen) {
  if (knapsack_.mode_ != CbcSettings::CGOff && knapsack_.proto_ == 0) {
    knapsack_.proto_ = new CglKnapsackCover();
  }
  gen = dynamic_cast<CglCutGenerator *>(knapsack_.proto_);

  return (knapsack_.mode_);
}

CbcSettings::CGMode CbcSettings::getMir(CglCutGenerator *&gen) {
  if (mir_.mode_ != CbcSettings::CGOff && mir_.proto_ == 0) {
    mir_.proto_ = new CglMixedIntegerRounding2();
  }
  gen = dynamic_cast<CglCutGenerator *>(mir_.proto_);

  return (mir_.mode_);
}

CbcSettings::CGMode CbcSettings::getProbing(CglCutGenerator *&gen) {
  if (probing_.mode_ != CbcSettings::CGOff && probing_.proto_ == 0) {
    probing_.proto_ = new CglProbing();
    probing_.proto_->setUsingObjective(probing_.usingObjective_);
    probing_.proto_->setMaxPass(probing_.maxPass_);
    probing_.proto_->setMaxPassRoot(probing_.maxPassRoot_);
    probing_.proto_->setMaxProbe(probing_.maxProbe_);
    probing_.proto_->setMaxProbeRoot(probing_.maxProbeRoot_);
    probing_.proto_->setMaxLook(probing_.maxLook_);
    probing_.proto_->setMaxLookRoot(probing_.maxLookRoot_);
    probing_.proto_->setMaxElements(probing_.maxElements_);
    probing_.proto_->setRowCuts(probing_.rowCuts_);
  }
  gen = dynamic_cast<CglCutGenerator *>(probing_.proto_);

  return (probing_.mode_);
}

CbcSettings::CGMode CbcSettings::getRedSplit(CglCutGenerator *&gen)

{
  if (redSplit_.mode_ != CbcSettings::CGOff && redSplit_.proto_ == 0) {
    redSplit_.proto_ = new CglRedSplit();
  }
  gen = dynamic_cast<CglCutGenerator *>(redSplit_.proto_);

  return (redSplit_.mode_);
}

CbcSettings::CGMode CbcSettings::getTwomir(CglCutGenerator *&gen) {
  if (twomir_.mode_ != CbcSettings::CGOff && twomir_.proto_ == 0) {
    twomir_.proto_ = new CglTwomir();
    twomir_.proto_->setMaxElements(twomir_.maxElements_);
  }
  gen = dynamic_cast<CglCutGenerator *>(twomir_.proto_);

  return (twomir_.mode_);
}

CbcSettings::HeurMode CbcSettings::getFeasPump(CbcHeuristic *&gen,
                                                        bool alwaysCreate)

{
  if (fpump_.mode_ != CbcSettings::HeurOff &&
      (fpump_.proto_ == 0 || alwaysCreate)) {
    if (fpump_.proto_) {
      delete fpump_.proto_;
    }
    fpump_.proto_ = new CbcHeuristicFPump(*model_);
    fpump_.proto_->setMaximumPasses(fpump_.iters_);
  }
  gen = dynamic_cast<CbcHeuristic *>(fpump_.proto_);

  return (fpump_.mode_);
}

/*
  Access functions for heuristics. These support lazy
  creation --- if mode_ is other than HeurOff, an object is created if
  necessary and a pointer is stored in proto_. The pointer is returned as
  a generic CbcHeuristic. The return value of the function
  is the value of mode_.

  Because the model may have changed, the default for heuristics is to delete
  any existing object and create a new one. This can be suppressed if desired.
*/

CbcSettings::HeurMode CbcSettings::getCombine(CbcHeuristic *&gen,
                                                       bool alwaysCreate)

{
  if (combine_.mode_ != CbcSettings::HeurOff &&
      (combine_.proto_ == 0 || alwaysCreate)) {
    if (combine_.proto_) {
      delete combine_.proto_;
    }
    //TODO Should we be passing a pointer here? Otherwise, making a copy
    combine_.proto_ = new CbcHeuristicLocal(*model_);
    combine_.proto_->setSearchType(combine_.trySwap_);
  }
  gen = dynamic_cast<CbcHeuristic *>(combine_.proto_);

  return (combine_.mode_);
}

CbcSettings::HeurMode CbcSettings::getGreedyCover(CbcHeuristic *&gen,
                                                           bool alwaysCreate)

{
  if (greedyCover_.mode_ != CbcSettings::HeurOff &&
      (greedyCover_.proto_ == 0 || alwaysCreate)) {
    if (greedyCover_.proto_) {
      delete greedyCover_.proto_;
    }
    greedyCover_.proto_ = new CbcHeuristicGreedyCover(*model_);
  }
  gen = dynamic_cast<CbcHeuristic *>(greedyCover_.proto_);

  return (greedyCover_.mode_);
}

CbcSettings::HeurMode
CbcSettings::getGreedyEquality(CbcHeuristic *&gen,
                                     bool alwaysCreate)

{
  if (greedyEquality_.mode_ != CbcSettings::HeurOff &&
      (greedyEquality_.proto_ == 0 || alwaysCreate)) {
    if (greedyEquality_.proto_) {
      delete greedyEquality_.proto_;
    }
    greedyEquality_.proto_ = new CbcHeuristicGreedyEquality(*model_);
  }
  gen = dynamic_cast<CbcHeuristic *>(greedyEquality_.proto_);

  return (greedyEquality_.mode_);
}

CbcSettings::HeurMode CbcSettings::getRounding(CbcHeuristic *&gen,
                                                        bool alwaysCreate)

{
  if (rounding_.mode_ != CbcSettings::HeurOff &&
      (rounding_.proto_ == 0 || alwaysCreate)) {
    if (rounding_.proto_) {
      delete rounding_.proto_;
    }
    rounding_.proto_ = new CbcRounding(*model_);
  }
  gen = dynamic_cast<CbcHeuristic *>(rounding_.proto_);

  return (rounding_.mode_);
}

CbcSettings::HeurMode
CbcSettings::getLocalTree(CbcTreeLocal *&localTree, 
                                bool alwaysCreate)

{
  if (localTree_.mode_ != CbcSettings::HeurOff &&
      (localTree_.proto_ == 0 || alwaysCreate)) {
    if (localTree_.proto_) {
      delete localTree_.proto_;
    }
    localTree_.proto_ = new CbcTreeLocal(
        model_, localTree_.soln_, localTree_.range_, localTree_.typeCuts_,
        localTree_.maxDiverge_, localTree_.timeLimit_, localTree_.nodeLimit_,
        localTree_.refine_);
  }
  localTree = localTree_.proto_;

  return (localTree_.mode_);
}

/*
  A bunch of little translation helper routines leading up to a version of
  setBaBStatus that figures it all out given a CbcModel and BACWhere code.
  This translation needs to be centralised to avoid sprinkling magic numbers
  all through the code.

  Be a bit careful with the translation routines --- they aren't sensitive to
  where the search stopped.
*/

CbcSettings::BACMajorStatus CbcSettings::translateMajor(int status)

{
  switch (status) {
  case -1: {
    return (CbcSettings::BACNotRun);
  }
  case 0: {
    return (CbcSettings::BACFinish);
  }
  case 1: {
    return (CbcSettings::BACStop);
  }
  case 2: {
    return (CbcSettings::BACAbandon);
  }
  case 5: {
    return (CbcSettings::BACUser);
  }
  default: { return (CbcSettings::BACInvalid); }
  }
}

CbcSettings::BACMinorStatus CbcSettings::translateMinor(int status)

{
  switch (status) {
  case -1: {
    return (CbcSettings::BACmInvalid);
  }
  case 0: {
    return (CbcSettings::BACmFinish);
  }
  case 1: {
    return (CbcSettings::BACmInfeas);
  }
  case 2: {
    return (CbcSettings::BACmGap);
  }
  case 3: {
    return (CbcSettings::BACmNodeLimit);
  }
  case 4: {
    return (CbcSettings::BACmTimeLimit);
  }
  case 5: {
    return (CbcSettings::BACmUser);
  }
  case 6: {
    return (CbcSettings::BACmSolnLimit);
  }
  case 7: {
    return (CbcSettings::BACmUbnd);
  }
  default: { return (CbcSettings::BACmOther); }
  }
}

/*
  A bit different --- given an OSI, use its interrogation functions to choose
  an appropriate BACMinorStatus code. Not everything matches up, eh?
*/
CbcSettings::BACMinorStatus
CbcSettings::translateMinor(const OsiSolverInterface *osi)

{
  if (osi->isProvenOptimal()) {
    return (CbcSettings::BACmFinish);
  } else if (osi->isProvenPrimalInfeasible()) {
    return (CbcSettings::BACmInfeas);
  } else if (osi->isProvenDualInfeasible()) {
    return (CbcSettings::BACmUbnd);
  } else {
    return (CbcSettings::BACmOther);
  }
}

/*
  A routine to set the bab_ status block given a CbcModel and an indication
  of where we're at in the search. Really, this is just a big mapping from
  CbcModel codes to CbcGeneric codes.
*/

void CbcSettings::setBaBStatus(CbcSettings::BACWhere where,
                                     bool haveAnswer,
                                     OsiSolverInterface *answerSolver)

{
  CbcSettings::BACMajorStatus major;
  CbcSettings::BACMinorStatus minor;

  major = translateMajor(model_->status());

  if (where == CbcSettings::BACwBareRoot ||
      where == CbcSettings::BACwIPPRelax) {
    minor = translateMinor(model_->solver());
  } else {
    minor = translateMinor(model_->secondaryStatus());
  }

  setBaBStatus(major, minor, where, haveAnswer, answerSolver);

  return;
}

/*
  Last, but not least, a routine to print the result.
*/

void CbcSettings::printBaBStatus()

{
  std::cout << "BAC result: stopped ";

  switch (bab_.where_) {
  case CbcSettings::BACwNotStarted: {
    std::cout << "before root relaxation";
    break;
  }
  case CbcSettings::BACwBareRoot: {
    std::cout << "after root relaxation";
    break;
  }
  case CbcSettings::BACwIPP: {
    std::cout << "after integer preprocessing";
    break;
  }
  case CbcSettings::BACwIPPRelax: {
    std::cout << "after solving preprocessed relaxation";
    break;
  }
  case CbcSettings::BACwBAC: {
    std::cout << "after branch-and-cut";
    break;
  }
  default: {
    std::cout << "!!invalid phase code!!";
    break;
  }
  }

  std::cout << std::endl << "    Branch-and-cut ";

  switch (bab_.majorStatus_) {
  case CbcSettings::BACNotRun: {
    std::cout << "never got started";
    break;
  }
  case CbcSettings::BACFinish: {
    std::cout << "finished";
    break;
  }
  case CbcSettings::BACStop: {
    std::cout << "stopped on a limit";
    break;
  }
  case CbcSettings::BACAbandon: {
    std::cout << "was abandoned";
    break;
  }
  case CbcSettings::BACUser: {
    std::cout << "stopped due to a user event";
    break;
  }
  default: {
    std::cout << "!!invalid major status code!!";
    break;
  }
  }

  std::cout << "; minor status is ";

  switch (bab_.minorStatus_) {
  case CbcSettings::BACmFinish: {
    std::cout << "optimal";
    break;
  }
  case CbcSettings::BACmInfeas: {
    std::cout << "infeasible";
    break;
  }
  case CbcSettings::BACmUbnd: {
    std::cout << "unbounded";
    break;
  }
  case CbcSettings::BACmGap: {
    std::cout << "reached specified integrality gap.";
    break;
  }
  case CbcSettings::BACmNodeLimit: {
    std::cout << "reached node limit";
    break;
  }
  case CbcSettings::BACmTimeLimit: {
    std::cout << "reached time limit";
    break;
  }
  case CbcSettings::BACmSolnLimit: {
    std::cout << "reached limit on number of solutions";
    break;
  }
  case CbcSettings::BACmUser: {
    std::cout << "stopped due to a user event";
    break;
  }
  case CbcSettings::BACmOther: {
    std::cout << "other";
    break;
  }
  default: {
    std::cout << "!!invalid minor status code!!";
    break;
  }
  }

  std::cout << "." << std::endl;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
