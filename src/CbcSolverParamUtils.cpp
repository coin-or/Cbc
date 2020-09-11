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

void addCbcSolverParams(int &numberParameters, CoinParamVec &parameters,
  CbcSolverSettings *cbcSettings)

{
  CbcSolverParam *param;
  std::string empty = "";

  param = new CbcSolverParam(CbcSolverParam::GENERALQUERY, "?",
                             "Print a list of commands",
                             CoinParam::displayPriorityNone);
  param->setPushFunc(doHelpParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::FULLGENERALQUERY, "???",
    "Print a list with *all* commands, even those hidden with `?'",
                          CoinParam::displayPriorityNone);
  param->setPushFunc(doHelpParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::ALLCOMMANDS, "allC!ommands",
                             "Whether to print less used commands",
                             "off", 0, CoinParam::displayPriorityNone);
  param->appendKwd("more");
  param->appendKwd("all");
  param->setPushFunc(pushCbcSolverKwdParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "For the sake of your sanity, only the more useful and simple commands are printed out on ?.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::HELP, "help",
    "Print out version, non-standard options and some help");
  param->setPushFunc(doHelpParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This prints out some help to get a user started. If you're seeing this message, you should be past that stage.");
  parameters.push_back(param);

  /*
      Built into CoinParam parsing. No additional actions required. doNothingParam
      simply prevents them from being reported as unimplemented.
    */
  param = new CbcSolverParam(CbcSolverParam::STDIN, "-",
    "Switch to interactive command line mode", CoinParam::displayPriorityNone);
  param->setPushFunc(doNothingParam);
  parameters.push_back(param);
  
  param = new CbcSolverParam(CbcSolverParam::STDIN, "stdin",
    "Switch to interactive command line mode", CoinParam::displayPriorityNone);
  param->setPushFunc(doNothingParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::PRINTVERSION,
    "version", "Print version");
  param->setPushFunc(doVersionParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::ARTIFICIALCOST, "artif!icialCost",
                             "Costs >= this treated as artificials in"
                             "feasibility pump",
                             0.0, COIN_DBL_MAX,
                             cbcSettings->getArtVarThreshold(),
                             CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverDblParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::BAB,
    "branch!AndCut", "Do Branch and Cut");
  param->setPushFunc(doBaCParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This does branch and cut. There are many parameters which can affect the performance.  First just try with default cbcSettings and look carefully at the log file.  Did cuts help?  Did they take too long?  Look at output to see which cuts were effective and then do some tuning.  You will see that the options for cuts are off, on, root and ifmove.  Off is obvious, on means that this cut generator will be tried in the branch and cut tree (you can fine tune using 'depth').  Root means just at the root node while 'ifmove' means that cuts will be used in the tree if they look as if they are doing some good and moving the objective value.  If pre-processing reduced the size of the problem or strengthened many coefficients then it is probably wise to leave it on.  Switch off heuristics which did not provide solutions.  The other major area to look at is the search.  Hopefully good solutions were obtained fairly early in the search so the important point is to select the best variable to branch on.  See whether strong branching did a good job - or did it just take a lot of iterations.  Adjust the strongBranching and trustPseudoCosts parameters.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::BKPIVOTINGSTRATEGY, "bkpivot!ing",
                          "Pivoting strategy used in Bron-Kerbosch algorithm",
                          0, 6, 3);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::BKMAXCALLS, "bkmaxcalls",
          "Maximum number of recursive calls made by Bron-Kerbosch algorithm",
                          1, COIN_INT_MAX, 1000);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::BKCLQEXTMETHOD, "bkclqext!method",
                             "Strategy used to extend violated cliques found by"
                             "BK Clique Cut Separation routine",
                             0, 5, 4);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "Sets the method used in the extension module of BK Clique Cut Separation routine: 0=no extension; 1=random; 2=degree; 3=modified degree; 4=reduced cost(inversely proportional); 5=reduced cost(inversely proportional) + modified degree");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::CLIQUECUTS,
    "clique!Cuts", "Whether to use clique cuts", "off",
                             cbcSettings->getCliqueMode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("ifmove");
  param->appendKwd("forceon");
  param->setObj(cbcSettings);
  param->setPushFunc(pushCbcSolverCutParam);
  param->setLongHelp(
    "This switches on clique cuts (either at root or in entire tree). See branchAndCut for information on options.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::CLQSTRENGTHENING, "clqstr!engthen",
         "Whether and when to perform Clique Strengthening preprocessing routine",
                          "after", 0);
  param->appendKwd("off");
  param->appendKwd("before");
  param->setPushFunc(pushCbcSolverKwdParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::COMBINE, "combine!Solutions",
                             "Whether to use combine solution heuristic",
                             "off", cbcSettings->getCombineMode());
  param->appendKwd("on");
  param->setPushFunc(pushCbcSolverBoolParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This switches on a heuristic which does branch and cut on the problem given by just using variables which have appeared in one or more solutions. It is obviously only tried after two or more solutions.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::BRANCHPRIORITY,
    "cost!Strategy", "Whether to use costs or column order as priorities",
    "off", 0);
  param->appendKwd("pri!orities");
  param->appendKwd("column!Order");
  param->setPushFunc(pushCbcSolverKwdParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This orders the variables in order of their absolute costs - with largest cost ones being branched on first.  This primitive strategy can be surprisingly effective.  The column order option is obviously not on costs but it's easy to implement.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::CPP,
                          "cpp!Generate", "Generates C++ code", 0, 4, 0);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "Once you like what the stand-alone solver does then this allows you to generate user_driver.cpp which approximates the code.  0 gives simplest driver, 1 generates saves and restores, 2 generates saves and restores even for variables at default value. 4 bit in cbc generates size dependent code rather than computed values.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::CPX, "cplex!Use",
                          "Whether to use Cplex!", "off", 0);
  param->appendKwd("on");
  param->setPushFunc(pushCbcSolverBoolParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    " If the user has Cplex, but wants to use some of Cbc's heuristics then you can!  If this is on, then Cbc will get to the root node and then hand over to Cplex.  If heuristics find a solution this can be significantly quicker.  You will probably want to switch off Cbc's cuts as Cplex thinks they are genuine constraints.  It is also probable that you want to switch off preprocessing, although for difficult problems it is worth trying both.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::CROSSOVER, "combine2!Solutions",
                             "Whether to use crossover solution heuristic",
                             "off", cbcSettings->getCrossoverMode());
  param->appendKwd("on");
  param->appendKwd("both");
  param->appendKwd("before");
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::CSVSTATISTICS, "csv!Statistics",
                          "Create one line of statistics",
                          cbcSettings->dfltDirectory_, CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverStrParam);
  param->setObj(cbcSettings);
   param->setLongHelp(
     "This appends statistics to given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to '', i.e. it must be set.  Adds header if file empty or does not exist.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::CUTDEPTH, "cutD!epth",
                             "Depth in tree at which to do cuts",
                             -1, 999999, cbcSettings->getCutDepth());
  param->setObj(cbcSettings);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setLongHelp(
    "Cut generators may be off, on only at the root, on if they look useful, and on at some interval.  If they are done every node then that is that, but it may be worth doing them every so often.  The original method was every so many nodes but it is more logical to do it whenever depth in tree is a multiple of K.  This option does that and defaults to -1 (off).");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::CUTLENGTH, "cutL!ength", "Length of a cut",
                          -1, COIN_INT_MAX, -1);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp("At present this only applies to Gomory cuts. -1 (default) leaves as is. Any value >0 says that all cuts <= this length can be generated both at root node and in tree. 0 says to use some dynamic lengths.  If value >=10,000,000 then the length in tree is value%10000000 - so 10000100 means unlimited length at root and 100 in tree.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::CUTOFFCONSTRAINT, "constraint!fromCutoff",
                          "Whether to use cutoff as constraint", "off", 0);
  param->appendKwd("on");
  param->appendKwd("variable");
  param->appendKwd("forceVariable");
  param->appendKwd("conflict");
  param->setPushFunc(pushCbcSolverKwdParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "For some problems, cut generators and general branching work better if the problem would be infeasible if the cost is too high. "
      "If this option is enabled, the objective function is added as a constraint which right hand side is set to the current cutoff value (objective value of best known solution)");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::CUTPASSINTREE, "passT!reeCuts",
                  "Number of rounds that cut generators are applied in the tree",
                          -COIN_INT_MAX, COIN_INT_MAX);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::CUTSTRATEGY,
    "cuts!OnOff", "Switches all cuts on or off", "off", 0);
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("ifmove");
  param->appendKwd("forceon");
  param->setPushFunc(pushCbcSolverCutParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This can be used to switch on or off all cuts (apart from Reduce and Split).  Then you can set individual ones off or on.  See branchAndCut for information on options.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DEBUG,
    "debug!In", "Read/write valid solution from/to file", "", CoinParam::displayPriorityNone);
  param->setObj(cbcSettings);
  param->setPushFunc(doDebugParam);
  param->setLongHelp(
    "This will read a solution file from the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to '', i.e. it must be set.\n\nIf set to create it will create a file called debug.file after B&C search; if set to createAfterPre it will create the file before undoing preprocessing.\n\nThe idea is that if you suspect a bad cut generator and you did not use preprocessing you can do a good run with debug set to 'create' and then switch on the cuts you suspect and re-run with debug set to 'debug.file'  Similarly if you do use preprocessing, but use createAfterPre.  The create case has the same effect as saveSolution.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DEPTHMINIBAB, "depth!MiniBab",
                          "Depth at which to try mini branch-and-bound",
                          -COIN_INT_MAX, COIN_INT_MAX, -1);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "Rather a complicated parameter but can be useful. -1 means off for large problems but on as if -12 for problems where rows+columns<500, -2 means use Cplex if it is linked in.  Otherwise if negative then go into depth first complete search fast branch and bound when depth>= -value-2 (so -3 will use this at depth>=1).  This mode is only switched on after 500 nodes.  If you really want to switch it off for small problems then set this to -999.  If >=0 the value doesn't matter very much.  The code will do approximately 100 nodes of fast branch and bound every now and then at depth>=5. The actual logic is too twisted to describe here.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DEXTRA3, "dextra3",
                          "Extra double parameter 3",
                          -COIN_DBL_MAX, COIN_DBL_MAX, 0.0,
                          CoinParam::displayPriorityNone);
  param->setPushFunc(pushCbcSolverDblParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DEXTRA4, "dextra4",
                          "Extra double parameter 4",
                          -COIN_DBL_MAX, COIN_DBL_MAX, 0.0,
                          CoinParam::displayPriorityNone);
  param->setPushFunc(pushCbcSolverDblParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DEXTRA5, "dextra5",
                          "Extra double parameter 5",
                          -COIN_DBL_MAX, COIN_DBL_MAX, 0.0,
                          CoinParam::displayPriorityNone);
  param->setPushFunc(pushCbcSolverDblParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DINS, "Dins",
                          "Whether to try Distance Induced Neighborhood Search",
                          "off", cbcSettings->getDinsMode());
  param->appendKwd("on");
  param->appendKwd("both");
  param->appendKwd("before");
  param->appendKwd("often");
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DIRECTORY,
    "directory", "Set Default directory for import etc.",
    cbcSettings->dfltDirectory_);
  param->setPushFunc(pushCbcSolverStrParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This sets the directory which import, export, saveModel, restoreModel etc. will use. It is initialized to the current directory.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DIRSAMPLE, "dirSample",
                          "Set directory where the COIN-OR sample problems are.",
                          cbcSettings->dfltDirectory_, CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverStrParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This sets the directory where the COIN-OR sample problems reside. It is used only when -unitTest is passed to clp. clp will pick up the test problems from this directory. It is initialized to '../../Data/Sample'");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DIRNETLIB, "dirNetlib",
                          "Set directory where the netlib problems are.",
                          cbcSettings->dfltDirectory_, CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverStrParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "This sets the directory where the netlib problems reside. One can get the netlib problems from COIN-OR or from the main netlib site. This parameter is used only when -netlib is passed to cbc. cbc will pick up the netlib problems from this directory. If clp is built without zlib support then the problems must be uncompressed. It is initialized to '../../Data/Netlib'");
  parameters.push_back(param);

   param = new CbcSolverParam(CbcSolverParam::DIRMIPLIB, "dirMiplib",
                          "Set directory where the miplib 2003 problems are.",
                          cbcSettings->dfltDirectory_, CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverStrParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This sets the directory where the miplib 2003 problems reside. One can get the miplib problems from COIN-OR or from the main miplib site. This parameter is used only when -miplib is passed to cbc. cbc will pick up the miplib problems from this directory. If cbc is built without zlib support then the problems must be uncompressed. It is initialized to '../../Data/miplib3'");
  parameters.push_back(param);

   param = new CbcSolverParam(CbcSolverParam::DIVEOPT, "diveO!pt", "Diving options",
                           -1, 20, -1, CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "If >2 && <20 then modify diving options - \
	 \n\t3 only at root and if no solution,  \
	 \n\t4 only at root and if this heuristic has not got solution, \
	 \n\t5 decay only if no solution, \
	 \n\t6 if depth <3 or decay, \
	 \n\t7 run up to 2 times if solution found 4 otherwise, \
	 \n\t>10 All only at root (DivingC normal as value-10), \
	 \n\t>20 All with value-20).");
   parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DIVEOPTSOLVES, "diveS!olves",
                          "Diving solve option",
                          -1, 200000, 100, CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "If >0 then do up to this many solves. However, the last digit is ignored and used for extra options: 1-3 enables fixing of satisfied integer variables (but not at bound), where 1 switches this off for that dive if the dive goes infeasible, and 2 switches it off permanently if the dive goes infeasible.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DIVINGC, "DivingC!oefficient",
                          "Whether to try Coefficient diving heuristic",
                          "off", cbcSettings->getDiveCoefficientMode());
  param->appendKwd("on");   
  param->appendKwd("both");   
  param->appendKwd("before");   
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DIVINGF, "DivingF!ractional",
                          "Whether to try Fractional diving heuristic",
                          "off", cbcSettings->getDiveFractionalMode());
  param->appendKwd("on");   
  param->appendKwd("both");   
  param->appendKwd("before");   
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DIVINGG, "DivingG!uided",
                          "Whether to try Guided diving heuristic",
                          "off", cbcSettings->getDiveGuidedMode());
  param->appendKwd("on");   
  param->appendKwd("both");   
  param->appendKwd("before");   
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DIVINGL, "DivingL!ineSearch",
                          "Whether to try Linesearch diving heuristic",
                          "off", cbcSettings->getDiveLineSearchMode());
  param->appendKwd("on");   
  param->appendKwd("both");   
  param->appendKwd("before");   
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DIVINGP, "DivingP!seudocost",
                          "Whether to try Pseudocost diving heuristic",
                          "off", cbcSettings->getDivePseudocostMode());
  param->appendKwd("on");   
  param->appendKwd("both");   
  param->appendKwd("before");   
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DIVINGS, "DivingS!ome",
                          "Whether to try Diving heuristics",
                          "off", cbcSettings->getDiveRandomMode());
  param->appendKwd("on");   
  param->appendKwd("both");   
  param->appendKwd("before");   
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
                     "This switches on a random diving heuristic at various times. One may prefer to individually turn diving heuristics on or off. ");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DIVINGV, "DivingV!ectorLength",
                          "Whether to try Vectorlength diving heuristic",
                          "off", cbcSettings->getDiveVectorLengthMode());
  param->appendKwd("on");   
  param->appendKwd("both");   
  param->appendKwd("before");   
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DJFIX, "fix!OnDj",
    "Try heuristic that fixes variables based on reduced costs",
                             -1.0e20, 1.0e20, cbcSettings->getDjFixThreshold());
  param->setPushFunc(pushCbcSolverDblParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "If set, integer variables with reduced costs greater than the specified value will be fixed before branch and bound - use with extreme caution!");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DOHEURISTIC, "doH!euristic",
                             "Do heuristics before any preprocessing",
                             "off", cbcSettings->getDiveCoefficientMode());
  param->appendKwd("on");   
  param->setPushFunc(pushCbcSolverBoolParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "Normally heuristics are done in branch and bound.  It may be useful to do them outside. Only those heuristics with 'both' or 'before' set will run. Doing this may also set cutoff, which can help with preprocessing.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DUMMY,
    "sleep", "for debug", 0, 9999, 0, CoinParam::displayPriorityNone);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "If passed to solver from ampl, then ampl will wait so that you can copy .nl file for debug.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::DW, "dw!Heuristic",
                          "Whether to try Dantzig Wolfe heuristic",
                          "off", cbcSettings->getDWMode());
  param->appendKwd("on");   
  param->appendKwd("both");   
  param->appendKwd("before");   
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This heuristic is very very compute intensive. It tries to find a Dantzig Wolfe structure and use that. " HEURISTICS_LONGHELP);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::ENVIRONMENT, "environ!ment",
                          "Read commands from environment",
                          CoinParam::displayPriorityNone);
  param->setPushFunc(doNothingParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "This starts reading from environment variable COIN_ENVIRONMENT.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::ERRORSALLOWED,
    "error!sAllowed", "Whether to allow import errors", "off", 0);
  param->appendKwd("on");
  param->setObj(cbcSettings);
  param->setLongHelp(
    "The default is not to use any model which had errors when reading the mps file.  Setting this to 'on' will allow all errors from which the code can recover simply by ignoring the error.  There are some errors from which the code can not recover, e.g., no ENDATA.  This has to be set before import, i.e., -errorsAllowed on -import xxxxxx.mps.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::EXIT, "end", "Stops execution");
  param->setPushFunc(doExitParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This stops execution; end, exit, quit and stop are synonyms.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::EXPERIMENT, "exper!iment",
                          "Whether to use testing features",
                          -1, 200000, 0, CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "Defines how adventurous you want to be in using new ideas. 0 then no new ideas, 1 fairly sensible, 2 a bit dubious, 3 you are on your own!");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::EXPORT,
    "export", "Export model as mps file",
    std::string("default.mps"));
  param->setPushFunc(pushCbcSolverStrParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This will write an MPS format file to the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to 'default.mps'. It can be useful to get rid of the original names and go over to using Rnnnnnnn and Cnnnnnnn.  This can be done by setting 'keepnames' off before importing mps file.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::EXTRA1, "extra1",
                          "Extra integer parameter 1",
                          -COIN_INT_MAX, COIN_INT_MAX, -1,
                          CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::EXTRA2, "extra2",
                          "Extra integer parameter 2",
                          -COIN_INT_MAX, COIN_INT_MAX, -1,
                          CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::EXTRA3, "extra3",
                          "Extra integer parameter 3",
                          -COIN_INT_MAX, COIN_INT_MAX, -1,
                          CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::EXTRA4, "extra4",
                          "Extra integer parameter 4",
                          -COIN_INT_MAX, COIN_INT_MAX, -1,
                          CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::EXTRA_VARIABLES, "extraV!ariables",
                          "Allow creation of extra integer variables",
                          "off", 0, CoinParam::displayPriorityLow);
  param->appendKwd("on");   
  param->setPushFunc(pushCbcSolverStrParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "Switches on a trivial re-formulation that introduces extra integer variables to group together variables with same cost.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::FAKECUTOFF, "pumpC!utoff",
                          "Fake cutoff for use in feasibility pump",
                          -COIN_DBL_MAX, COIN_DBL_MAX, 0.0,
                          CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverDblParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "A value of 0.0 means off. Otherwise, add a constraint forcing objective below this value in feasibility pump");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::FAKEINCREMENT, "pumpI!ncrement",
                          "Fake increment for use in feasibility pump",
                          -COIN_DBL_MAX, COIN_DBL_MAX, 0.0,
                          CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverDblParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "A value of 0.0 means off. Otherwise, add a constraint forcing objective below this value in feasibility pump");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::FLOWCUTS,
    "flow!CoverCuts", "Whether to use Flow Cover cuts", "off",
                             cbcSettings->getFlowMode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("ifmove");
  param->appendKwd("forceon");
  param->setObj(cbcSettings);
  param->setPushFunc(pushCbcSolverCutParam);
  param->setLongHelp(
    "This switches on flow cover cuts (either at root or in entire tree).  See branchAndCut for information on options.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::FPUMP,
    "feas!ibilityPump", "Whether to try Feasibility Pump", "off",
                             cbcSettings->getFeasPumpMode());
  param->appendKwd("on");
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This switches on feasibility pump heuristic at root. This is due to Fischetti and Lodi and uses a sequence of LPs to try and get an integer feasible solution.  Some fine tuning is available by passFeasibilityPump.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::FPUMPITS,
    "passF!easibilityPump", "How many passes in feasibility pump",
                             0, 10000, cbcSettings->getFeasPumpIters());
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This fine tunes the Feasibility Pump heuristic by doing more or fewer passes.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::FPUMPTUNE, "pumpT!une",
                          "Dubious ideas for feasibility pump",
                          0, 100000000, 0);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This fine tunes Feasibility Pump \
     \n\t>=10000000 use as objective weight switch \
     \n\t>=1000000 use as accumulate switch \
     \n\t>=1000 use index+1 as number of large loops \
     \n\t==100 use objvalue +0.05*fabs(objvalue) as cutoff OR fakeCutoff if set \
     \n\t%100 == 10,20 affects how each solve is done \
     \n\t1 == fix ints at bounds, 2 fix all integral ints, 3 and continuous at bounds. If accumulate is on then after a major pass, variables which have not moved are fixed and a small branch and bound is tried.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::FPUMPTUNE2, "moreT!une",
                          "Yet more dubious ideas for feasibility pump",
                          0, 100000000, 0, CoinParam::displayPriorityNone);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setLongHelp(
    "Yet more ideas for Feasibility Pump \
     \n\t/100000 == 1 use box constraints and original obj in cleanup \
     \n\t/1000 == 1 Pump will run twice if no solution found \
     \n\t/1000 == 2 Pump will only run after root cuts if no solution found \
     \n\t/1000 >10 as above but even if solution found \
     \n\t/100 == 1,3.. exact 1.0 for objective values \
     \n\t/100 == 2,3.. allow more iterations per pass \
     \n\t n fix if value of variable same for last n iterations.");
  param->setObj(cbcSettings);

  param = new CbcSolverParam(CbcSolverParam::GMICUTS, "GMI!Cuts",
                          "Whether to use alternative Gomory cuts",
                          "off", cbcSettings->getGMIMode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("ifmove");
  param->appendKwd("forceon");
  param->appendKwd("endonly");
  param->appendKwd("long");
  param->appendKwd("longroot");
  param->appendKwd("longifmove");
  param->appendKwd("forcelongon");
  param->appendKwd("longendonly");
  param->setPushFunc(pushCbcSolverCutParam);
  param->setObj(cbcSettings);
  param->setLongHelp(CUTS_LONGHELP
    " This version is by Giacomo Nannicini and may be more robust than gomoryCuts.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::GOMORYCUTS, "gomory!Cuts",
    "Whether to use Gomory cuts", "off",
                             cbcSettings->getGomoryMode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("ifMove");
  param->appendKwd("forceon");
  param->appendKwd("onGlobal");
  param->appendKwd("forceandglobal");
  param->appendKwd("forcelongon");
  param->appendKwd("long");
  param->appendKwd("shorter");
  param->setObj(cbcSettings);
  param->setPushFunc(pushCbcSolverCutParam);
  param->setLongHelp(
      "The original cuts - beware of imitations!  Having gone out of favor, \
they are now more fashionable as LP solvers are more robust and they interact well \
with other cuts.  They will almost always give cuts (although in this executable \
they are limited as to number of variables in cut).  However the cuts may be dense \
so it is worth experimenting (Long allows any length). "
    CUTS_LONGHELP
    " Reference: https://github.com/coin-or/Cgl/wiki/CglGomory");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::GREEDY,
    "greedy!Heuristic", "Whether to use a greedy heuristic", "off",
                             cbcSettings->getGreedyCoverMode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "Switches on a pair of greedy heuristic which will try and obtain a solution.  It may just fix a percentage of variables and then try a small branch and cut run.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::HEURISTICSTRATEGY,
    "heur!isticsOnOff", "Switches most heuristics on or off", "off", 0);
  param->appendKwd("on");
  param->setPushFunc(pushCbcSolverBoolParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This can be used to switch on or off all heuristics.  Then you can set individual ones off or on.  CbcTreeLocal is not included as it dramatically alters search.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::HEUROPTIONS, "hOp!tions",
                             "Heuristic options",
                          -COIN_INT_MAX, COIN_INT_MAX, 0,
                          CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "Value 1 stops heuristics immediately if the allowable gap has been reached. Other values are for the feasibility pump - 2 says do exact number of passes given, 4 only applies if an initial cutoff has been given and says relax after 50 passes, while 8 will adapt the cutoff rhs after the first solution if it looks as if the code is stalling.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::IMPORT,
    "import", "Import model from mps file",
    cbcSettings->lastMpsIn_);
  param->setPushFunc(doImportParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This will read an MPS format file from the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to '', i.e., it must be set.  If you have libgz then it can read compressed files 'xxxxxxxx.gz'.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::INTPRINT,
    "printi!ngOptions", "Print options", "normal", 0);
  param->appendKwd("integer");
  param->appendKwd("special");
  param->appendKwd("rows");
  param->appendKwd("all");
  param->setPushFunc(pushCbcSolverKwdParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This changes the amount and format of printing a solution:\nnormal - nonzero column variables \ninteger - nonzero integer column variables\nspecial - in format suitable for OsiRowCutDebugger\nrows - nonzero column variables and row activities\nall - all column variables and row activities.\n\nFor non-integer problems 'integer' and 'special' act like 'normal'.  Also see printMask for controlling output.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::KNAPSACKCUTS,
    "knapsack!Cuts", "Whether to use Knapsack cuts", "off",
                             cbcSettings->getKnapsackMode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("ifmove");
  param->appendKwd("forceon");
  param->setObj(cbcSettings);
  param->setPushFunc(pushCbcSolverCutParam);
  param->setLongHelp(
    "This switches on knapsack cuts (either at root or in entire tree). See branchAndCut for information on options.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::LANDPCUTS, "lift!AndProjectCuts",
                          "Whether to use lift-and-project cuts", "off",
                             cbcSettings->getLandPMode()) ;
  param->appendKwd("on") ;
  param->appendKwd("root") ;
  param->appendKwd("ifmove") ;
  param->appendKwd("forceon") ;
  param->setPushFunc(pushCbcSolverCutParam);
  param->setObj(cbcSettings) ;
  param->setLongHelp(
    "This switches on lift-and-project cuts (either at root or in entire tree). See branchAndCut for information on options."
    	) ;
  parameters.push_back(param) ;

  param = new CbcSolverParam(CbcSolverParam::LAGOMORYCUTS, "lagomory!Cuts",
                          "Whether to use Lagrangean Gomory cuts",
                          "off", cbcSettings->getLaGomoryMode());
  param->appendKwd("endonlyroot");
  param->appendKwd("endcleanroot");
  param->appendKwd("root");
  param->appendKwd("endonly");
  param->appendKwd("endclean");
  param->appendKwd("endboth");
  param->appendKwd("onlyaswell");
  param->appendKwd("cleanaswell");
  param->appendKwd("bothaswell");
  param->appendKwd("onlyinstead");
  param->appendKwd("cleaninstead");
  param->appendKwd("bothinstead");
  param->appendKwd("onlyaswellroot");
  param->appendKwd("cleanaswellroot");
  param->appendKwd("bothaswellroot");
  param->setPushFunc(pushCbcSolverCutParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This is a gross simplification of 'A Relax-and-Cut Framework for Gomory's Mixed-Integer Cuts' by Matteo Fischetti & Domenico Salvagnin.  This simplification just uses original constraints while modifying objective using other cuts. So you don't use messy constraints generated by Gomory etc. A variant is to allow non messy cuts e.g. clique cuts. So 'only' does this while 'clean' also allows integral valued cuts.  'End' is recommended and waits until other cuts have finished before it does a few passes. The length options for gomory cuts are used.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::LATWOMIRCUTS, "latwomir!Cuts",
                          "Whether to use Lagrangean Twomir cuts",
                          "off", cbcSettings->getLaTwomirMode());
  param->appendKwd("endonlyroot");
  param->appendKwd("endcleanroot");
  param->appendKwd("endbothroot");
  param->appendKwd("endonly");
  param->appendKwd("endclean");
  param->appendKwd("endboth");
  param->appendKwd("onlyaswell");
  param->appendKwd("cleanaswell");
  param->appendKwd("bothaswell");
  param->appendKwd("onlyinstead");
  param->appendKwd("cleaninstead");
  param->appendKwd("bothinstead");
  param->setPushFunc(pushCbcSolverCutParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "This is a Lagrangean relaxation for Twomir cuts.  See \
  lagomoryCuts for description of options.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::LOCALTREE,
    "local!TreeSearch", "Whether to use local tree search", "off",
                             cbcSettings->getLocalTreeMode());
  param->appendKwd("on");
  param->setPushFunc(pushCbcSolverStrParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This switches on a local search algorithm when a solution is found.  This is from Fischetti and Lodi and is not really a heuristic although it can be used as one. When used from this program it has limited functionality.  It is not controlled by heuristicsOnOff.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::LOGLEVEL,
    "log!Level", "Level of detail in CBC output.",
                             -1, 999999, cbcSettings->getLogLevel());
  param->setObj(cbcSettings);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setLongHelp(
    "If set to 0 then there should be no output in normal circumstances. A value of 1 is probably the best value for most uses, while 2 and 3 give more information.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::LPLOGLEVEL,
    "log!Level", "Level of detail in LP solver output.",
                             -1, 999999, cbcSettings->getLpLogLevel());
  param->setObj(cbcSettings);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setLongHelp(
    "If set to 0 then there should be no output in normal circumstances. A value of 1 is probably the best value for most uses, while 2 and 3 give more information.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::MAXSAVEDSOLS, "maxSaved!Solutions",
                          "Maximum number of solutions to save",
                          0, COIN_INT_MAX, 1);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "Number of solutions to save.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::MAXSLOWCUTS, "slow!cutpasses",
                          "Maximum number of rounds for slower cut generators",
                          -1, COIN_INT_MAX, 10);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "Some cut generators are fairly slow - this limits the number of times they are tried. The cut generators identified as 'may be slow' at present are Lift and project cuts and both versions of Reduce and Split cuts.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::MESSAGES, "mess!ages",
    "Controls whether standardised message prefix is printed", "off", 0);
  param->appendKwd("on");
  param->setPushFunc(pushCbcSolverBoolParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "By default, messages have a standard prefix, such as:\n   Clp0005 2261  Objective 109.024 Primal infeas 944413 (758)\nbut this program turns this off to make it look more friendly.  It can be useful to turn them back on if you want to be able to 'grep' for particular messages or if you intend to override the behavior of a particular message.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::MIPLIB,
    "miplib", "Do some of miplib test set");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::MIPSTART, "mips!tart",
                          "reads an initial feasible solution from file",
                          std::string("mipstart.sln"));
  param->setPushFunc(pushCbcSolverStrParam);
  param->setObj(cbcSettings);
  param->setLongHelp("\
The MIPStart allows one to enter an initial integer feasible solution \
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
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::MIRCUTS, "mixed!IntegerRoundingCuts",
    "Whether to use Mixed Integer Rounding cuts", "off", 
                             cbcSettings->getMirMode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("ifmove");
  param->appendKwd("forceon");
  param->setObj(cbcSettings);
  param->setPushFunc(pushCbcSolverCutParam);
  param->setLongHelp(
    "This switches on mixed integer rounding cuts (either at root or in entire tree).  See branchAndCut for information on options.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::MOREMOREMIPOPTIONS, "more2!MipOptions",
                          "More more dubious options for mip",
                          -1, COIN_INT_MAX, 0, CoinParam::displayPriorityNone);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::MULTIPLEROOTS, "multiple!RootPasses",
                          "Do multiple root passes to collect cuts and solutions",
                          0, COIN_INT_MAX, 0, CoinParam::displayPriorityNone);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "Solve (in parallel, if enabled) the root phase this number of times, each with its own different seed, and collect all solutions and cuts generated. The actual format is aabbcc where aa is the number of extra passes; if bb is non zero, then it is number of threads to use (otherwise uses threads setting); and cc is the number of times to do root phase. The solvers do not interact with each other.  However if extra passes are specified then cuts are collected and used in later passes - so there is interaction there. Some parts of this implementation have their origin in idea of Andrea Lodi, Matteo Fischetti, Michele Monaci, Domenico Salvagnin, and Andrea Tramontani.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::NAIVE, "naive!Heuristics",
                             "Whether to try some stupid heuristic",
                             "off", cbcSettings->getNaiveHeurMode(),
                             CoinParam::displayPriorityLow);
  param->appendKwd("on");
  param->appendKwd("both");
  param->appendKwd("before");
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "This is naive heuristics which, e.g., fix all integers with costs to zero!. "
      HEURISTICS_LONGHELP);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::NEXTBESTSOLUTION, "nextB!estSolution",
                          "Prints next best saved solution to file");
  param->setPushFunc(doNothingParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "To write best solution, just use solution.  This prints next best (if exists) and then deletes it. This will write a primitive solution file to the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to 'stdout'.  The amount of output can be varied using printi!ngOptions or printMask.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::NODESTRATEGY, "node!Strategy",
        "What strategy to use to select the next node from the branch and cut tree",
                          "hybrid", 0);
  param->appendKwd("fewest");
  param->appendKwd("depth");
  param->appendKwd("upfewest");
  param->appendKwd("downfewest");
  param->appendKwd("updepth");
  param->appendKwd("downdepth");
  param->setPushFunc(pushCbcSolverKwdParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "Normally before a feasible solution is found, CBC will choose a node with fewest infeasibilities. Alternatively, one may choose tree-depth as the criterion. This requires the minimal amount of memory, but may take a long time to find the best solution. Additionally, one may specify whether up or down branches must be selected first (the up-down choice will carry on after a first solution has been bound). The choice 'hybrid' does breadth first on small depth nodes and then switches to 'fewest'.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::ODDWHEELCUTS, "oddwheel!Cuts",
                          "Whether to use odd wheel cuts",
                          "off", cbcSettings->getOddWheelMode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("ifmove");
  param->appendKwd("forceon");
  param->appendKwd("onglobal");
  param->setPushFunc(pushCbcSolverCutParam);
  param->setObj(cbcSettings);
  param->setLongHelp("This switches on odd-wheel inequalities (either at root or in entire tree).");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::ODDWEXTMETHOD, "oddwext!method",
                          "Strategy used to search for wheel centers for the cuts found by Odd Wheel Cut Separation routine",
                          0, 2, 2);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
        "Sets the method used in the extension module of Odd Wheel Cut Separation routine: 0=no extension; 1=one variable; 2=clique");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::ORBITAL, "Orbit!alBranching",
                          "Whether to try orbital branching",
                          "off", 0);
  param->appendKwd("on");
  param->appendKwd("slow!ish");
  param->appendKwd("strong");
  param->appendKwd("force");
  param->appendKwd("simple");
  param->appendKwd("more!printing");
  param->setPushFunc(pushCbcSolverKwdParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "This switches on Orbital branching. Value 'on' just adds orbital, 'strong' tries extra fixing in strong branching.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::OUTDUPROWS,
    "outDup!licates",
    "Takes duplicate rows, etc., out of the integer model", CoinParam::displayPriorityNone);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::OUTPUTFORMAT,
    "output!Format", "Which output format to use", 1, 6, 2);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "Normally export will be done using normal representation for numbers and two values per line.  You may want to do just one per line (for grep or suchlike) and you may wish to save with absolute accuracy using a coded version of the IEEE value. A value of 2 is normal. Otherwise, odd values give one value per line, even values two.  Values of 1 and 2 give normal format, 3 and 4 give greater precision, 5 and 6 give IEEE values.  When exporting a basis, 1 does not save values, 2 saves values, 3 saves with greater accuracy and 4 saves in IEEE format.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::PIVOTANDFIX, "pivotAndF!ix",
                          "Whether to try Pivot and Fix heuristic",
                          "off", cbcSettings->getPivotAndFixMode());
  param->appendKwd("on");
  param->appendKwd("both");
  param->appendKwd("before");
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP);
  parameters.push_back(param);

#if 0
  param = new CbcSolverParam(CbcSolverParam::PIVOTANDCOMPLEMENT, "pivotAndC!omplement",
                          "Whether to try Pivot and Complement heuristic",
                          "off", cbcSettings->getPivotAndComplementMode());
  param->appendKwd("on");
  param->appendKwd("both");
  param->appendKwd("before");
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP);
  parameters.push_back(param);
#endif
  /*
      In order for initialisation to work properly, the order of the keywords here
      must match the order of the enum IPPControl in CbcSolverSettings.hpp.
    */
  param = new CbcSolverParam(CbcSolverParam::PREPROCESS,
    "preprocess", "Whether to use integer preprocessing", "off",
                             cbcSettings->getIPPMode());
  param->appendKwd("on");
  param->appendKwd("save");
  param->appendKwd("equal");
  param->appendKwd("sos");
  param->appendKwd("trysos");
  param->appendKwd("equalall");
  param->appendKwd("strategy");
  param->setPushFunc(pushCbcSolverKwdParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This tries to reduce size of the model in a similar way to presolve and it also tries to strengthen the model. This can be very useful and is worth trying.  save option saves on file presolved.mps.  equal will turn <= cliques into ==.  sos will create sos sets if all 0-1 in sets (well one extra is allowed) and no overlaps.  trysos is same but allows any number extra. equalall will turn all valid inequalities into equalities with integer slacks. strategy is as on but uses CbcStrategy.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::PREPROCNAMES, "PrepN!ames",
                          "If column names will be kept in pre-processed model",
                          "off", 0);
  param->appendKwd("on");
  param->setPushFunc(pushCbcSolverBoolParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "Normally the preprocessed model has column names replaced by new names C0000...\
Setting this option to on keeps original names in variables which still exist in the preprocessed problem");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::PRINTMASK,
    "printM!ask",
    "Control printing of solution with a regular expression", empty);
  param->setPushFunc(doPrintMaskParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "If set then only those names which match mask are printed in a solution. '?' matches any character and '*' matches any set of characters.  The default is '' (unset) so all variables are printed. This is only active if model has names.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::PRINTOPTIONS,
    "pO!ptions", "Dubious print options", 0, COIN_INT_MAX, 0, CoinParam::displayPriorityNone);
  param->setObj(cbcSettings);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setLongHelp(
    "If this is greater than 0 then presolve will give more information and branch and cut will give statistics");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::PRIORITYIN, "prio!rityIn",
                          "Import priorities etc from file",
                          std::string("priorities.txt"));
  param->setPushFunc(pushCbcSolverStrParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This will read a file with priorities from the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to '', i.e. it must be set.  This can not read from compressed files. File is in csv format with allowed headings - name, number, priority, direction, up, down, solution.  Exactly one of name and number must be given.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::PROBINGCUTS,
    "probing!Cuts", "Whether to use Probing cuts", "off",
                             cbcSettings->getProbingMode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("ifmove");
  param->appendKwd("forceon");
  param->appendKwd("forceonbut");
  param->setObj(cbcSettings);
  param->setPushFunc(pushCbcSolverCutParam);
  param->setLongHelp(
    "This switches on probing cuts (either at root or in entire tree). See branchAndCut for information on options.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::PROCESSTUNE, "tune!PreProcess",
                          "Dubious tuning parameters for preprocessing",
                          0, COIN_INT_MAX, 0, CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
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
     Now aa 99 has special meaning i.e. just one simple presolve.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::PROXIMITY, "proximity!Search",
                          "Whether to do proximity search heuristic",
                          "off", cbcSettings->getProximityMode());
  param->appendKwd("on");
  param->appendKwd("both");
  param->appendKwd("before");
  param->appendKwd("10");
  param->appendKwd("100");
  param->appendKwd("300");
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "This heuristic looks for a solution close to the incumbent solution (Fischetti and Monaci, 2012). The idea is to define a sub-MIP without additional constraints but with a modified objective function intended to attract the search in the proximity of the incumbent. The approach works well for 0-1 MIPs whose solution landscape is not too irregular (meaning the there is reasonable probability of finding an improved solution by flipping a small number of binary variables), in particular when it is applied to the first heuristic solutions found at the root node. "
      HEURISTICS_LONGHELP); // Can also set different maxNode cbcSettings by plusnnnn (and are 'on'(on==30)).
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::RANDOMSEED, "randomC!bcSeed",
                          "Random seed for Cbc",
                          -1, COIN_INT_MAX, -1);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "Allows initialization of the random seed for pseudo-random numbers used in heuristics such as the Feasibility Pump to decide whether to round up or down. "
      "The special value of 0 lets Cbc use the time of the day for the initial seed.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::RANDROUND, "randomi!zedRounding",
                          "Whether to try randomized rounding heuristic",
                          "off", cbcSettings->getRandRoundMode());
  param->appendKwd("on");
  param->appendKwd("both");
  param->appendKwd("before");
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::REDSPLITCUTS,
    "reduce!AndSplitCuts",
    "Whether to use Reduce-and-Split cuts", "off",
                             cbcSettings->getRedSplitMode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("longOn");
  param->appendKwd("longroot");
  param->setObj(cbcSettings);
  param->setPushFunc(pushCbcSolverCutParam);
  param->setLongHelp(
    "This switches on reduce and split cuts (either at root or in entire tree). See branchAndCut for information on options.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::REDSPLIT2CUTS,
    "reduce2!AndSplitCuts",
    "Whether to use Reduce-and-Split cuts - style 2", "off",
                             cbcSettings->getRedSplit2Mode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("longOn");
  param->appendKwd("longroot");
  param->setObj(cbcSettings);
  param->setPushFunc(pushCbcSolverCutParam);
  param->setLongHelp(
    "This switches on reduce and split cuts (either at root or in entire tree). See branchAndCut for information on options.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::RENS, "Rens",
                       "Whether to try Relaxation Enforced Neighborhood Search",
                          "off", cbcSettings->getRensMode());
  param->appendKwd("on");
  param->appendKwd("both");
  param->appendKwd("before");
  param->appendKwd("200");
  param->appendKwd("1000");
  param->appendKwd("10000");
  param->appendKwd("dj");
  param->appendKwd("djbefore");
  param->appendKwd("usesolution");
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP
      " Value 'on' just does 50 nodes. 200, 1000, and 10000 does that many nodes.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::RESIDCAPCUTS, "residual!CapacityCuts",
                          "Whether to use Residual Capacity cuts",
                          "off", cbcSettings->getResidCapMode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("ifmove");
  param->appendKwd("forceon");
  param->setPushFunc(pushCbcSolverCutParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP
      " Value 'on' just does 50 nodes. 200, 1000, and 10000 does that many nodes.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::RINS, "Rins",
                    "Whether to try Relaxed Induced Neighborhood Search",
                             "off", cbcSettings->getRinsMode());
  param->appendKwd("on");
  param->appendKwd("both");
  param->appendKwd("before");
  param->appendKwd("often");
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::ROUNDING,
    "round!ingHeuristic", "Whether to use Rounding heuristic", "off",
                             cbcSettings->getRoundingMode());
  param->appendKwd("on");
  param->setPushFunc(pushCbcSolverBoolParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This switches on a simple (but effective) rounding heuristic at each node of tree.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::SAVESOL, "saveS!olution",
                          "saves solution to file", std::string("solution.sln"));
  param->setPushFunc(pushCbcSolverStrParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "This will write a binary solution file to the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to 'solution.file'.  To read the file use fread(int) twice to pick up number of rows and columns, then fread(double) to pick up objective value, then pick up row activities, row duals, column activities and reduced costs - see bottom of CbcParam.cpp for code that reads or writes file. If name contains '_fix_read_' then does not write but reads and will fix all variables");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::SHOWUNIMP,
    "unimp!lemented", "Report unimplemented commands.", CoinParam::displayPriorityNone);
  param->setPushFunc(doUnimplementedParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::SMALLBAB, "fraction!forBAB",
                          "Fraction in feasibility pump",
                          1.0e-5, 1.1, 0.5, CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverDblParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "After a pass in the feasibility pump, variables which have not moved about are fixed and if the preprocessed model is smaller than this fraction of the original problem, a few nodes of branch and bound are done on the reduced problem.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::SOLUTION,
    "solu!tion", "Prints solution to file",
    std::string("stdout"));
  param->setPushFunc(doSolutionParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This will write a primitive solution file to the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to 'stdout'.  The amount of output can be varied using printi!ngOptions or printMask.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::SOLVECONTINUOUS, "initialS!olve",
    "Solve to continuous optimum");
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This just solves the problem to the continuous optimum, without adding any cuts.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::SOS,
    "sos!Options", "Whether to use SOS from AMPL", "off", 0);
  param->appendKwd("on");
  param->setPushFunc(pushCbcSolverBoolParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "Normally if AMPL says there are SOS variables they should be used, but sometimes they should be turned off - this does so.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::SOSPRIORITIZE, "sosP!rioritize",
                          "How to deal with SOS priorities",
                          "off", 0);
  param->appendKwd("high");
  param->appendKwd("low");
  param->appendKwd("orderhigh");
  param->appendKwd("orderlow");
  param->setPushFunc(pushCbcSolverKwdParam);
  param->setObj(cbcSettings);
   param->setLongHelp(
       "This sets priorities for SOS.  Values 'high' and 'low' just set a priority relative to the for integer variables.  Value 'orderhigh' gives first highest priority to the first SOS and integer variables a low priority.  Value 'orderlow' gives integer variables a high priority then SOS in order.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::STATISTICS, "stat!istics",
                          "Print some statistics");
  param->setObj(cbcSettings);
  param->setLongHelp(
      "This command prints some statistics for the current model. If log level >1 then more is printed. These are for presolved model if presolve on (and unscaled).");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::STRATEGY, "strat!egy",
                          "Switches on groups of features",
                          "default", 0);
  param->appendKwd("easy");
  param->appendKwd("aggressive");
  param->setPushFunc(pushCbcSolverKwdParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "This turns on newer features. Use 0 for easy problems, 1 is default, 2 is aggressive. 1 uses Gomory cuts with a tolerance of 0.01 at the root node, does a possible restart after 100 nodes if many variables could be fixed, activates a diving and RINS heuristic, and makes the feasibility pump more aggressive."); // This does not apply to unit tests (where 'experiment' may have similar effects)
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::STRENGTHEN,
    "strengthen", "Create strengthened problem");
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This creates a new problem by applying the root node cuts. All tight constraints will be in resulting problem.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::STRONGSTRATEGY, "expensive!Strong",
                          "Whether to do even more strong branching",
                          0, COIN_INT_MAX, 0, CoinParam::displayPriorityNone);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "Strategy for extra strong branching. 0 is normal strong branching. 1, 2, 4, and 6 does strong branching on all fractional variables if at the root node (1), \
at depth less than modifier (2), objective equals best possible (4), or at depth less than modifier and objective equals best possible (6). 11, 12, 14, and 16 are like 1, 2, 4, and 6, respecitively, but do strong branching on all integer (incl. non-fractional) variables. Values >= 100 are used to specify a depth limit (value/100), otherwise 5 is used. If the values >= 100, then above rules are applied to value%100.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::TESTOSI, "testO!si", 
                          "Test OsiObject stuff",
                          -1, COIN_INT_MAX, -1, CoinParam::displayPriorityNone);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  parameters.push_back(param);

#ifdef CBC_THREAD
  param = new CbcSolverParam(CbcSolverParam::THREADS, "thread!s",
                          "Number of threads to try and use",
                          -100, 100000, 0, CoinParam::displayPriorityLow);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "To use multiple threads, set threads to number wanted.  It may be better to use one or two more than number of cpus available.  If 100+n then n threads and search is repeatable (maybe be somewhat slower), if 200+n use threads for root cuts, 400+n threads used in sub-trees.");
  parameters.push_back(param);
#endif
    
  param = new CbcSolverParam(CbcSolverParam::TIGHTENFACTOR, "tighten!Factor",
    "Tighten bounds using value times largest activity at continuous solution",
    1.0, 1.0e20);
  param->setPushFunc(pushCbcSolverDblParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This sleazy trick can help on some problems.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::TIMEMODE, "timeM!ode",
                          "Whether to use CPU or elapsed time",
                          "cpu", 0);
  param->appendKwd("elapsed");
  param->setPushFunc(pushCbcSolverKwdParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "cpu uses CPU time for stopping, while elapsed uses elapsed time. \
(On Windows, elapsed time is always used).");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::TWOMIRCUTS,
    "two!MirCuts",
    "Whether to use Two phase Mixed Integer Rounding cuts", "off",
                             cbcSettings->getTwomirMode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("ifmove");
  param->appendKwd("forceon");
  param->setObj(cbcSettings);
  param->setPushFunc(pushCbcSolverCutParam);
  param->setLongHelp(
    "This switches on two phase mixed integer rounding cuts (either at root or in entire tree). See branchAndCut for information on options.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::UNITTEST, "unitTest", "Do unit test");
  param->setObj(cbcSettings);
  param->setLongHelp(
    "This exercises the unit test.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::USECGRAPH, "cgraph",
                          "Whether to use the conflict graph-based preprocessing and cut separation routines.", "on", 0);
  param->appendKwd("off");
  param->appendKwd("clq");
  param->setPushFunc(pushCbcSolverKwdParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
        "This switches the conflict graph-based preprocessing and cut separation routines (CglBKClique, CglOddWheel and CliqueStrengthening) on or off. Values:\
\n\toff: turns these routines off;\
\n\ton: turns these routines on;\
\n\tclq: turns these routines off and enables the cut separator of CglClique.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::USERCBC, "userCbc",
                          "Hand coded Cbc stuff",
                          0, COIN_INT_MAX, 0, CoinParam::displayPriorityNone);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "There are times (e.g., when using AMPL interface) when you may wish to do something unusual.  Look for USERCBC in main driver and modify sample code.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::USESOLUTION,
    "force!Solution",
    "Whether to use given solution as crash for BAB", "off", 0);
  param->appendKwd("on");
  param->setPushFunc(pushCbcSolverBoolParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "If on then tries to branch to solution given by AMPL or priorities file.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::VERBOSE,
    "verbose", "Switches on longer help on single ?",
    0, 15, cbcSettings->verbose_, CoinParam::displayPriorityNone);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
    "Set to 1 to get short help with ? list, 2 to get long help.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::VND, "Vnd!VariableNeighborhoodSearch",
                          "Whether to try Variable Neighborhood Search",
                          "off", cbcSettings->getVndMode());
  param->appendKwd("on");
  param->appendKwd("both");
  param->appendKwd("before");
  param->appendKwd("intree");
  param->setPushFunc(pushCbcSolverHeurParam);
  param->setObj(cbcSettings);
  param->setLongHelp(HEURISTICS_LONGHELP);
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::VUBTRY, "vub!heuristic",
                          "Type of VUB heuristic",
                          -2, 20, -1, CoinParam::displayPriorityNone);
  param->setPushFunc(pushCbcSolverIntParam);
  param->setObj(cbcSettings);
  param->setLongHelp(
      "This heuristic tries to fix some integer variables.");
  parameters.push_back(param);

  param = new CbcSolverParam(CbcSolverParam::ZEROHALFCUTS, "zero!HalfCuts",
                          "Whether to use zero half cuts",
                          "off", cbcSettings->getZeroHalfMode());
  param->appendKwd("on");
  param->appendKwd("root");
  param->appendKwd("ifmove");
  param->appendKwd("forceon");
  param->appendKwd("onglobal");
  param->setPushFunc(pushCbcSolverCutParam);
  param->setObj(cbcSettings);
  param->setLongHelp(CUTS_LONGHELP
      " This implementation was written by Alberto Caprara.");
  parameters.push_back(param);

  numberParameters = parameters.size();
  assert(((unsigned)numberParameters) <= parameters.capacity());

#if 0
  param = new CbcSolverParam;
  CbcGenSolvers::setupSolverParam(*param);
  param->setObj(cbcSettings);
  parameters.push_back(param);
#endif
  
  return;
}

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
   case CbcSolverParam::DJFIX: {
      cbcSettings->setDjFixMode(CbcSolverSettings::ParamOn, val);
      break;
   }
   case CbcSolverParam::ARTIFICIALCOST: {
      cbcSettings->setArtVarMode(CbcSolverSettings::ParamOn, val);
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

  std::string str = genParam->kwdVal();
  CbcSolverParam::CbcSolverParamCode code = genParam->paramCode();

  int retval = 0;
  /*
      Figure out what we're doing and set the relevant field.
    */
  switch (code) {
   case CbcSolverParam::ALLCOMMANDS: {
      if (str == "off") {
         cbcSettings->setCommandMode(CbcSolverSettings::CMOff);
      } else if (str == "more") {
         cbcSettings->setCommandMode(CbcSolverSettings::CMMore);
      } else if (str == "all") {
         cbcSettings->setCommandMode(CbcSolverSettings::CMAll);
      } else {
         std::cerr
            << "pushCbcSolverKwdParam(ALLCOMMANDS): unrecognized keyword `"
            << str << "'." << std::endl;
         retval = -1;
      }
      break;
   }
   case CbcSolverParam::CLQSTRENGTHENING: {
      if (str == "off") {
         cbcSettings->setClqStrMode(CbcSolverSettings::ClqStrOff);
      } else if (str == "after") {
         cbcSettings->setClqStrMode(CbcSolverSettings::ClqStrAfter);
      } else if (str == "before") {
         cbcSettings->setClqStrMode(CbcSolverSettings::ClqStrBefore);
      } else {
         std::cerr
            << "pushCbcSolverKwdParam(CLQSTRENGTHENING): unrecognized keyword `"
            << str << "'." << std::endl;
         retval = -1;
      }
      break;
   }
   case CbcSolverParam::BRANCHPRIORITY: {
      if (str == "off") {
         cbcSettings->priorityMode_ = CbcSolverSettings::BPOff;
      } else if (str == "priorities") {
         cbcSettings->priorityMode_ = CbcSolverSettings::BPCost;
      }
      if (str == "columnOrder") {
         cbcSettings->priorityMode_ = CbcSolverSettings::BPOrder;
      } else {
         std::cerr
            << "pushCbcSolverKwdParam(BRANCHPRIORITY): unrecognized keyword `"
            << str << "'." << std::endl;
         retval = -1;
      }
      break;
   }
   case CbcSolverParam::CUTOFFCONSTRAINT: {
      if (str == "off") {
         cbcSettings->setCutoffMode(CbcSolverSettings::COOff);
      } else if (str == "on") {
         cbcSettings->setCutoffMode(CbcSolverSettings::COOn);
      } else if (str == "variable") {
         cbcSettings->setCutoffMode(CbcSolverSettings::COVariable);
      } else if (str == "forceVariable") {
         cbcSettings->setCutoffMode(CbcSolverSettings::COForceVariable);
      } else if (str == "conflict") {
         cbcSettings->setCutoffMode(CbcSolverSettings::COConflict);
      } else {
         std::cerr
            << "pushCbcSolverKwdParam(INTPRINT): unrecognized keyword `"
            << str << "'." << std::endl;
         retval = -1;
      }
      break;
   }
   case CbcSolverParam::INTPRINT: {
      if (str == "normal") {
         cbcSettings->setIntPrintMode(CbcSolverSettings::PMNormal);
      } else if (str == "integer") {
         cbcSettings->setIntPrintMode(CbcSolverSettings::PMInteger);
      } else if (str == "special") {
         cbcSettings->setIntPrintMode(CbcSolverSettings::PMSpecial);
      } else if (str == "rows") {
         cbcSettings->setIntPrintMode(CbcSolverSettings::PMRows);
      } else if (str == "all") {
         cbcSettings->setIntPrintMode(CbcSolverSettings::PMAll);
      } else {
         std::cerr
            << "pushCbcSolverKwdParam(INTPRINT): unrecognized keyword `"
            << str << "'." << std::endl;
         retval = -1;
      }
      break;
   }
   case CbcSolverParam::NODESTRATEGY: {
      if (str == "hybrid") {
         cbcSettings->setNodeStrategy(CbcSolverSettings::NSHybrid);
      } else if (str == "fewest") {
         cbcSettings->setNodeStrategy(CbcSolverSettings::NSFewest);
      } else if (str == "depth") {
         cbcSettings->setNodeStrategy(CbcSolverSettings::NSDepth);
      } else if (str == "upfewest") {
         cbcSettings->setNodeStrategy(CbcSolverSettings::NSUpFewest);
      } else if (str == "updepth") {
         cbcSettings->setNodeStrategy(CbcSolverSettings::NSUpDepth);
      } else if (str == "downdepth") {
         cbcSettings->setNodeStrategy(CbcSolverSettings::NSDownDepth);
      } else {
         std::cerr
            << "pushCbcSolverKwdParam(NODESTRATEGY): unrecognized keyword `"
            << str << "'." << std::endl;
         retval = -1;
      }
      break;
   }
   case CbcSolverParam::ORBITAL: {
      if (str == "off") {
         cbcSettings->setOrbitalStrategy(CbcSolverSettings::OBOff);
      } else if (str == "on") {
         cbcSettings->setOrbitalStrategy(CbcSolverSettings::OBOn);
      } else if (str == "slowish") {
         cbcSettings->setOrbitalStrategy(CbcSolverSettings::OBSlowish);
      } else if (str == "strong") {
         cbcSettings->setOrbitalStrategy(CbcSolverSettings::OBStrong);
      } else if (str == "force") {
         cbcSettings->setOrbitalStrategy(CbcSolverSettings::OBForce);
      } else if (str == "simple") {
         cbcSettings->setOrbitalStrategy(CbcSolverSettings::OBSimple);
      } else if (str == "moreprinting") {
         cbcSettings->setOrbitalStrategy(CbcSolverSettings::OBMorePrinting);
      } else {
         std::cerr
            << "pushCbcSolverKwdParam(NODESTRATEGY): unrecognized keyword `"
            << str << "'." << std::endl;
         retval = -1;
      }
      break;
   }
   case CbcSolverParam::PREPROCESS: {
      if (str == "off") {
         cbcSettings->setIPPMode(CbcSolverSettings::IPPOff);
      } else if (str == "on") {
         cbcSettings->setIPPMode(CbcSolverSettings::IPPOn);
      } else if (str == "save") {
         cbcSettings->setIPPMode(CbcSolverSettings::IPPSave);
      } else if (str == "equal") {
         cbcSettings->setIPPMode(CbcSolverSettings::IPPEqual);
      } else if (str == "sos") {
         cbcSettings->setIPPMode(CbcSolverSettings::IPPSOS);
      } else if (str == "trysos") {
         cbcSettings->setIPPMode(CbcSolverSettings::IPPTrySOS);
      } else if (str == "equalall") {
         cbcSettings->setIPPMode(CbcSolverSettings::IPPEqualAll);
      } else if (str == "strategy") {
         cbcSettings->setIPPMode(CbcSolverSettings::IPPStrategy);
      } else {
         std::cerr
            << "pushCbcSolverKwdParam(PREPROCESS): unrecognized keyword `"
            << str << "'." << std::endl;
         retval = -1;
      }
      break;
   }
   case CbcSolverParam::SOSPRIORITIZE: {
      if (str == "off") {
         cbcSettings->setSOSStrategy(CbcSolverSettings::SOSOff);
      } else if (str == "high") {
         cbcSettings->setSOSStrategy(CbcSolverSettings::SOSHigh);
      } else if (str == "low") {
         cbcSettings->setSOSStrategy(CbcSolverSettings::SOSLow);
      } else if (str == "orderhigh") {
         cbcSettings->setSOSStrategy(CbcSolverSettings::SOSOrderHigh);
      } else if (str == "orderlow") {
         cbcSettings->setSOSStrategy(CbcSolverSettings::SOSOrderLow);
      } else {
         std::cerr
            << "pushCbcSolverKwdParam(SOSPRIORITIZE): unrecognized keyword `"
            << str << "'." << std::endl;
         retval = -1;
      }
      break;
   }
   case CbcSolverParam::STRATEGY: {
      if (str == "default") {
         cbcSettings->setStrategyMode(CbcSolverSettings::StrategyDefault);
      } else if (str == "easy") {
         cbcSettings->setStrategyMode(CbcSolverSettings::StrategyEasy);
      } else if (str == "aggressive") {
         cbcSettings->setStrategyMode(CbcSolverSettings::StrategyAggressive);
      } else {
         std::cerr
            << "pushCbcSolverKwdParam(STRATEGY): unrecognized keyword `"
            << str << "'." << std::endl;
         retval = -1;
      }
      break;
   }
   case CbcSolverParam::TIMEMODE: {
      if (str == "cpu") {
         cbcSettings->setClockType(CbcSolverSettings::ClockCpu);
      } else if (str == "elapsed") {
         cbcSettings->setClockType(CbcSolverSettings::ClockElapsed);
      } else {
         std::cerr
            << "pushCbcSolverKwdParam(TIMEMODE): unrecognized keyword `"
            << str << "'." << std::endl;
         retval = -1;
      }
      break;
   }
   case CbcSolverParam::USECGRAPH: {
      if (str == "off") {
         cbcSettings->setCGraphMode(CbcSolverSettings::CGraphOff);
      } else if (str == "on") {
         cbcSettings->setCGraphMode(CbcSolverSettings::CGraphOn);
      } else if (str == "clique") {
         cbcSettings->setCGraphMode(CbcSolverSettings::CGraphClique);
      } else {
         std::cerr
            << "pushCbcSolverKwdParam(USECGRAPH): unrecognized keyword `"
            << str << "'." << std::endl;
         retval = -1;
      }
      break;
   }
   default: {
      std::cerr
         << "pushCbcSolverKwdParam: unrecognized parameter code `"
         << code << "'." << std::endl;
      retval = -1;
      break;
   }
  }

  return (retval);
}

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

  std::string str = genParam->kwdVal();
  CbcSolverParam::CbcSolverParamCode code = genParam->paramCode();

  int retval = 0;

  CbcSolverSettings::OnOffMode mode = CbcSolverSettings::ParamEndMarker;

  /*
      Figure out what we're doing and set the relevant field.
    */
  switch (code) {
   case CbcSolverParam::CPX:
   case CbcSolverParam::DOHEURISTIC: 
   case CbcSolverParam::MESSAGES:
   case CbcSolverParam::PREPROCNAMES:
   case CbcSolverParam::SOS:
   case CbcSolverParam::USESOLUTION: {
      if (str == "off" ){
         mode = CbcSolverSettings::ParamOff;
      } else if (str == "on") {
         mode = CbcSolverSettings::ParamOn;
      } else {
         break;
      }
   }
   default: {
      std::cerr
         << "pushCbcSolverBoolParam: unrecognized parameter `"
         << code << "'." << std::endl;
      retval = -1;
      break;
   }
  }

  switch (code) {
   case CbcSolverParam::CPX: {
      cbcSettings->setCPXMode(mode);
      break;
   }
   case CbcSolverParam::DOHEURISTIC: {
      cbcSettings->setDoHeuristicMode(mode);
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
  
  if (mode == CbcSolverSettings::HeurEndMarker){
     std::cerr
        << "pushCbcSolverBoolParam: unrecognized keyword `"
        << str << "'." << std::endl;
     return (retval);
  }

  return (retval);
}

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

  std::string str = genParam->kwdVal();
  CbcSolverParam::CbcSolverParamCode code = genParam->paramCode();

  int retval = 0;

  CbcSolverSettings::HeurMode mode = CbcSolverSettings::HeurEndMarker;

  /*
      First translate the keyword into the correct CGControl enum value.
  */

  // Check the common keywords that mos heuristics use
  switch (code) {
   case CbcSolverParam::GREEDY:
   case CbcSolverParam::COMBINE: 
     break;
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
   case CbcSolverParam::VND:
      if (str == "off" ){
         mode = CbcSolverSettings::HeurOff;
      } else if (str == "on") {
         mode = CbcSolverSettings::HeurOn;
      } else if (str == "both") {
         mode = CbcSolverSettings::HeurBoth;
      } else if (str == "before") {
         mode = CbcSolverSettings::HeurBefore;
      } else {
         break;
      }
   default: {
      std::cerr
         << "pushCbcSolverHeurParam: unrecognized heuristic type `"
         << code << "'." << std::endl;
      retval = -1;
      break;
   }
  }

  // Check the common parameters
  switch (code) {
   case CbcSolverParam::COMBINE: 
     if (str == "off" ){
         mode = CbcSolverSettings::HeurOff;
     } else if (str == "on" ){
         mode = CbcSolverSettings::HeurOn;
     } else {
        break;
     }
   case CbcSolverParam::DINS: 
     if (str == "often" ){
         mode = CbcSolverSettings::HeurOften;
     } else {
        break;
     }
   case CbcSolverParam::GREEDY: 
     if (str == "off" ){
         mode = CbcSolverSettings::HeurOff;
     } else if (str == "on" ){
         mode = CbcSolverSettings::HeurOn;
     } else if (str == "root" ){
         mode = CbcSolverSettings::HeurRoot;
     } else {
        break;
     }
   case CbcSolverParam::PROXIMITY: 
     if (str == "10" ){
         mode = CbcSolverSettings::HeurTen;
     } else if (str == "100" ){
         mode = CbcSolverSettings::HeurOneHundred;
     } else if (str == "300" ){
         mode = CbcSolverSettings::HeurThreeHundred;
     } else {
        break;
     }
   case CbcSolverParam::RENS: 
     if (str == "200" ){
         mode = CbcSolverSettings::HeurTwoHundred;
     } else if (str == "1000" ){
         mode = CbcSolverSettings::HeurOneThousand;
     } else if (str == "10000" ){
         mode = CbcSolverSettings::HeurTenThousand;
     } else if (str == "dj" ){
         mode = CbcSolverSettings::HeurDj;
     } else if (str == "djbefore" ){
         mode = CbcSolverSettings::HeurDjBefore;
     } else if (str == "usesolution" ){
         mode = CbcSolverSettings::HeurUseSolution;
     } else {
        break;
     }
   case CbcSolverParam::RINS: 
     if (str == "often" ){
         mode = CbcSolverSettings::HeurOften;
     } else {
        break;
     }
   case CbcSolverParam::VND: 
     if (str == "intree" ){
         mode = CbcSolverSettings::HeurInTree;
     } else {
        break;
     }
   default:
     break;
  }
     
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


  if (mode == CbcSolverSettings::HeurEndMarker){
     std::cerr
        << "pushCbcSolverHeurParam: unrecognized keyword `"
        << str << "'." << std::endl;
     return (retval);
  }

  return (retval);
}

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

  std::string str = genParam->kwdVal();
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

  CbcSolverSettings::CGMode mode = CbcSolverSettings::CGEndMarker;

  /*
      First translate the keyword into the correct CGControl enum value.
  */

  // Check the common parameters that all but two cut classes use
  switch (code) {
   case CbcSolverParam::LAGOMORYCUTS:
   case CbcSolverParam::LATWOMIRCUTS: 
   case CbcSolverParam::REDSPLITCUTS:
   case CbcSolverParam::REDSPLIT2CUTS:
      break;
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
   case CbcSolverParam::ZEROHALFCUTS:
      if (str == "off") {
         mode = CbcSolverSettings::CGOff;
      } else if (str == "on") {
         mode = CbcSolverSettings::CGOn;
      } else if (str == "root") {
         mode = CbcSolverSettings::CGRoot;
      } else if (str == "ifmove") {
         mode = CbcSolverSettings::CGIfMove;
      } else if (str == "forceon") {
         mode = CbcSolverSettings::CGForceOn;
      }
   default: {
     std::cerr
        << "pushCbcSolverCutParam: unrecognized cut type `"
        << str << "'." << std::endl;
     return (retval);
   }
  }
     
  // Now, add some additional parameters for different classes
  switch (code) {
   case CbcSolverParam::GOMORYCUTS:
     if (str == "onglobal") {
        mode = CbcSolverSettings::CGOnGlobal;
     } else if (str == "forceandglobal") {
        mode = CbcSolverSettings::CGForceAndGlobal;
     } else if (str == "forcelongon") {
        mode = CbcSolverSettings::CGForceLongOn;
     } else if (str == "longer") {
        mode = CbcSolverSettings::CGLonger;
     } else if (str == "shorter") {
        mode = CbcSolverSettings::CGShorter;
     } else{
        break;
     }
   case CbcSolverParam::GMICUTS:
     if (str == "long") {
        mode = CbcSolverSettings::CGLong;
     } else if (str == "longroot") {
        mode = CbcSolverSettings::CGLongRoot;
     } else if (str == "longifmove") {
        mode = CbcSolverSettings::CGLongIfMove;
     } else if (str == "forcelongon") {
        mode = CbcSolverSettings::CGForceLongOn;
     } else if (str == "longendonly") {
        mode = CbcSolverSettings::CGLongEndOnly;
     } else{
        break;
     }
   case CbcSolverParam::LAGOMORYCUTS:
     if (str == "off") {
        mode = CbcSolverSettings::CGOff;
     } else if (str == "root") {
        mode = CbcSolverSettings::CGRoot;
     } else if (str == "onlyaswellroot") {
        mode = CbcSolverSettings::CGOnlyAsWellRoot;
     } else if (str == "cleanaswellroot") {
        mode = CbcSolverSettings::CGCleanAsWellRoot;
     } else if (str == "bothaswellroot") {
        mode = CbcSolverSettings::CGCleanBothAsWellRoot;
     }
     // Here, we intentionally drop through to the next set
   case CbcSolverParam::LATWOMIRCUTS:
     if (str == "endonlyroot") {
        mode = CbcSolverSettings::CGEndOnlyRoot;
     } else if (str == "endcleanroot") {
        mode = CbcSolverSettings::CGEndCleanRoot;
     } else if (str == "endonly") {
        mode = CbcSolverSettings::CGEndOnly;
     } else if (str == "endclean") {
        mode = CbcSolverSettings::CGEndClean;
     } else if (str == "endboth") {
        mode = CbcSolverSettings::CGEndBoth;
     } else if (str == "onlyaswell") {
        mode = CbcSolverSettings::CGOnlyAsWell;
     } else if (str == "cleanaswell") {
        mode = CbcSolverSettings::CGCleanAsWell;
     } else if (str == "bothaswell") {
        mode = CbcSolverSettings::CGBothAsWell;
     } else if (str == "onlyinstead") {
        mode = CbcSolverSettings::CGOnlyInstead;
     } else if (str == "cleaninstead") {
        mode = CbcSolverSettings::CGCleanInstead;
     } else if (str == "bothinstead") {
        mode = CbcSolverSettings::CGBothInstead;
     } else {
        break;
     }
   case CbcSolverParam::ODDWHEELCUTS:
     if (str == "onglobal") {
        mode = CbcSolverSettings::CGOnGlobal;
     } else {
        break;
     }
   case CbcSolverParam::PROBINGCUTS:
     if (str == "forceonbut") {
        mode = CbcSolverSettings::CGForceOnBut;
     } else {
        break;
     }
   case CbcSolverParam::REDSPLITCUTS:
   case CbcSolverParam::REDSPLIT2CUTS:
     if (str == "off") {
        mode = CbcSolverSettings::CGOff;
     } else if (str == "on") {
        mode = CbcSolverSettings::CGOn;
     } else if (str == "root") {
        mode = CbcSolverSettings::CGRoot;
     } else if (str == "longon") {
        mode = CbcSolverSettings::CGLongOn;
     } else if (str == "longroot") {
        mode = CbcSolverSettings::CGLongRoot;
     } else {
        break;
     }
   default:
     break;
  }
     
  
  if (mode == CbcSolverSettings::CGEndMarker){
     std::cerr
        << "pushCbcSolverCutParam: unrecognized keyword `"
        << str << "'." << std::endl;
     return (retval);
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
  default: {
    std::cerr
      << "pushCbcSolverCutParam: internal confusion!" << std::endl;
    return (-1);
  }
  }

  return (0);
}

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
