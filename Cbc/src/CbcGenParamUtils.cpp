/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file is part of cbc-generic.
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <string>
#include <cassert>

#include "CoinFileIO.hpp"

#include "CoinFinite.hpp"
#include "CoinParam.hpp"

#include "CbcModel.hpp"

#include "CbcGenParam.hpp"
#include "CbcGenCtlBlk.hpp"

/*! \file CbcGenParamUtils
    \brief Implementation functions for CbcGenParam parameters.
*/

namespace {

char svnid[] = "$Id: CbcGenParamUtils.cpp 1173 2009-06-04 09:44:10Z forrest $" ;

}

namespace CbcGenSolvers {
void setupSolverParam(CbcGenParam &solverParam) ;
}

namespace CbcGenParamUtils {

/*
  Function to add cbc-generic control parameters to the cbc-generic parameter
  vector. Were needed, defaults are drawn from ctlBlk-> This function is a
  friend of CbcGenCtlBlk.
*/

void addCbcGenParams (int &numberParameters, CoinParamVec &parameters,
                      CbcGenCtlBlk *ctlBlk)

{
    CbcGenParam *param ;
    std::string empty = "" ;

    param = new CbcGenParam(CbcGenParam::GENERALQUERY, "?",
                            "Print a list of commands", false) ;
    param->setPushFunc(doHelpParam) ;
    param->setObj(ctlBlk) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::FULLGENERALQUERY, "???",
                            "Print a list with *all* commands, even those hidden with `?'", false) ;
    param->setPushFunc(doHelpParam) ;
    param->setObj(ctlBlk) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::PRINTVERSION,
                            "version", "Print version") ;
    param->setPushFunc(doVersionParam) ;
    param->setObj(ctlBlk) ;
    parameters.push_back(param) ;

    /*
      Built into CoinParam parsing. No additional actions required. doNothingParam
      simply prevents them from being reported as unimplemented.
    */
    param = new CbcGenParam(CbcGenParam::STDIN, "-",
                            "Switch to interactive command line mode", false) ;
    param->setPushFunc(doNothingParam) ;
    parameters.push_back(param) ;
    param = new CbcGenParam(CbcGenParam::STDIN, "stdin",
                            "Switch to interactive command line mode", false) ;
    param->setPushFunc(doNothingParam) ;
    param->setObj(ctlBlk) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::BAB,
                            "branch!AndCut", "Do Branch and Cut") ;
    param->setPushFunc(doBaCParam) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This does branch and cut. There are many parameters which can affect the performance.  First just try with default settings and look carefully at the log file.  Did cuts help?  Did they take too long?  Look at output to see which cuts were effective and then do some tuning.  You will see that the options for cuts are off, on, root and ifmove.  Off is obvious, on means that this cut generator will be tried in the branch and cut tree (you can fine tune using 'depth').  Root means just at the root node while 'ifmove' means that cuts will be used in the tree if they look as if they are doing some good and moving the objective value.  If pre-processing reduced the size of the problem or strengthened many coefficients then it is probably wise to leave it on.  Switch off heuristics which did not provide solutions.  The other major area to look at is the search.  Hopefully good solutions were obtained fairly early in the search so the important point is to select the best variable to branch on.  See whether strong branching did a good job - or did it just take a lot of iterations.  Adjust the strongBranching and trustPseudoCosts parameters."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::CPP,
                            "cpp!Generate", "Generates C++ code", -1, 50000) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "Once you like what the stand-alone solver does then this allows you to generate user_driver.cpp which approximates the code.  0 gives simplest driver, 1 generates saves and restores, 2 generates saves and restores even for variables at default value. 4 bit in cbc generates size dependent code rather than computed values."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::CLIQUECUTS,
                            "clique!Cuts", "Whether to use clique cuts", "off",
                            ctlBlk->clique_.action_) ;
    param->appendKwd("on") ;
    param->appendKwd("root") ;
    param->appendKwd("ifmove") ;
    param->appendKwd("forceOn") ;
    param->setObj(ctlBlk) ;
    param->setPushFunc(pushCbcGenCutParam) ;
    param->setLongHelp(
        "This switches on clique cuts (either at root or in entire tree). See branchAndCut for information on options."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::CUTDEPTH,
                            "cutD!epth", "Depth in tree at which to do cuts",
                            -1, 999999, ctlBlk->cutDepth_) ;
    param->setObj(ctlBlk) ;
    param->setPushFunc(pushCbcGenIntParam) ;
    param->setLongHelp(
        "Cut generators may be off, on only at the root, on if they look useful, and on at some interval.  If they are done every node then that is that, but it may be worth doing them every so often.  The original method was every so many nodes but it is more logical to do it whenever depth in tree is a multiple of K.  This option does that and defaults to -1 (off)."
    ) ;
    parameters.push_back(param) ;


    param = new CbcGenParam(CbcGenParam::CUTSTRATEGY,
                            "cuts!OnOff", "Switches all cuts on or off", "off", 0) ;
    param->appendKwd("on") ;
    param->appendKwd("root") ;
    param->appendKwd("ifmove") ;
    param->appendKwd("forceOn") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This can be used to switch on or off all cuts (apart from Reduce and Split).  Then you can set individual ones off or on.  See branchAndCut for information on options."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::COMBINE,
                            "combine!Solutions",
                            "Whether to use combine solution heuristic", "off",
                            ctlBlk->combine_.action_) ;
    param->appendKwd("on") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This switches on a heuristic which does branch and cut on the problem given by just using variables which have appeared in one or more solutions. It is obviously only tried after two or more solutions."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::COSTSTRATEGY,
                            "cost!Strategy", "Whether to use costs or column order as priorities",
                            "off", 0) ;
    param->appendKwd("pri!orities") ;
    param->appendKwd("column!Order") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This orders the variables in order of their absolute costs - with largest cost ones being branched on first.  This primitive strategy can be surprisingly effective.  The column order option is obviously not on costs but it's easy to implement."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::DEBUG,
                            "debug!In", "Read/write valid solution from/to file", "", false) ;
    param->setObj(ctlBlk) ;
    param->setPushFunc(doDebugParam) ;
    param->setLongHelp(
        "This will read a solution file from the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to '', i.e. it must be set.\n\nIf set to create it will create a file called debug.file after B&C search; if set to createAfterPre it will create the file before undoing preprocessing.\n\nThe idea is that if you suspect a bad cut generator and you did not use preprocessing you can do a good run with debug set to 'create' and then switch on the cuts you suspect and re-run with debug set to 'debug.file'  Similarly if you do use preprocessing, but use createAfterPre.  The create case has the same effect as saveSolution."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::DIRECTORY,
                            "directory", "Set Default directory for import etc.",
                            ctlBlk->dfltDirectory_) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This sets the directory which import, export, saveModel, restoreModel etc. will use. It is initialized to the current directory."
    ) ;
    param->setPushFunc(pushCbcGenStrParam) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::EXIT, "end", "Stops execution") ;
    param->setPushFunc(doExitParam) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This stops execution; end, exit, quit and stop are synonyms."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::ERRORSALLOWED,
                            "error!sAllowed", "Whether to allow import errors", "off", 0) ;
    param->appendKwd("on") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "The default is not to use any model which had errors when reading the mps file.  Setting this to 'on' will allow all errors from which the code can recover simply by ignoring the error.  There are some errors from which the code can not recover, e.g., no ENDATA.  This has to be set before import, i.e., -errorsAllowed on -import xxxxxx.mps."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::EXIT, "exit", "Stops execution") ;
    param->setPushFunc(doExitParam) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This stops execution; end, exit, quit and stop are synonyms."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::EXPORT,
                            "export", "Export model as mps file",
                            std::string("default.mps")) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This will write an MPS format file to the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to 'default.mps'."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::DJFIX, "fix!OnDj",
                            "Try heuristic that fixes variables based on reduced costs",
                            -1.0e20, 1.0e20, ctlBlk->djFix_.threshold_) ;
    param->setPushFunc(pushCbcGenDblParam) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "If set, integer variables with reduced costs greater than the specified value will be fixed before branch and bound - use with extreme caution!"
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::FPUMP,
                            "feas!ibilityPump", "Whether to try Feasibility Pump", "off",
                            ctlBlk->fpump_.action_) ;
    param->appendKwd("on") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This switches on feasibility pump heuristic at root. This is due to Fischetti and Lodi and uses a sequence of LPs to try and get an integer feasible solution.  Some fine tuning is available by passFeasibilityPump."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::FPUMPITS,
                            "passF!easibilityPump", "How many passes in feasibility pump",
                            0, 10000, ctlBlk->fpump_.iters_) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This fine tunes the Feasibility Pump heuristic by doing more or fewer passes."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::FLOWCUTS,
                            "flow!CoverCuts", "Whether to use Flow Cover cuts", "off",
                            ctlBlk->flow_.action_) ;
    param->appendKwd("on") ;
    param->appendKwd("root") ;
    param->appendKwd("ifmove") ;
    param->appendKwd("forceOn") ;
    param->setObj(ctlBlk) ;
    param->setPushFunc(pushCbcGenCutParam) ;
    param->setLongHelp(
        "This switches on flow cover cuts (either at root or in entire tree).  See branchAndCut for information on options."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::GOMORYCUTS, "gomory!Cuts",
                            "Whether to use Gomory cuts", "off",
                            ctlBlk->gomory_.action_) ;
    param->appendKwd("on") ;
    param->appendKwd("root") ;
    param->appendKwd("ifmove") ;
    param->appendKwd("forceOn") ;
    param->setObj(ctlBlk) ;
    param->setPushFunc(pushCbcGenCutParam) ;
    param->setLongHelp(
        "The original cuts - beware of imitations!  Having gone out of favor, they are now more fashionable as LP solvers are more robust and they interact well with other cuts.  They will almost always give cuts (although in this executable they are limited as to number of variables in cut).  However the cuts may be dense so it is worth experimenting. See branchAndCut for information on options."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::GREEDY,
                            "greedy!Heuristic", "Whether to use a greedy heuristic", "off",
                            ctlBlk->greedyCover_.action_) ;
    param->appendKwd("on") ;
    param->appendKwd("root") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "Switches on a pair of greedy heuristic which will try and obtain a solution.  It may just fix a percentage of variables and then try a small branch and cut run."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::HEURISTICSTRATEGY,
                            "heur!isticsOnOff", "Switches most heuristics on or off", "off", 0) ;
    param->appendKwd("on") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This can be used to switch on or off all heuristics.  Then you can set individual ones off or on.  CbcTreeLocal is not included as it dramatically alters search."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::KNAPSACKCUTS,
                            "knapsack!Cuts", "Whether to use Knapsack cuts", "off",
                            ctlBlk->knapsack_.action_) ;
    param->appendKwd("on") ;
    param->appendKwd("root") ;
    param->appendKwd("ifmove") ;
    param->appendKwd("forceOn") ;
    param->setObj(ctlBlk) ;
    param->setPushFunc(pushCbcGenCutParam) ;
    param->setLongHelp(
        "This switches on knapsack cuts (either at root or in entire tree). See branchAndCut for information on options."
    ) ;
    parameters.push_back(param) ;

    /*
      param = new CbcGenParam(CbcGenParam::LANDPCUTS,
    	"lift!AndProjectCuts","Whether to use lift-and-project cuts","off",
    	ctlBlk->landp_.action_) ;
      param->appendKwd("on") ;
      param->appendKwd("root") ;
      param->appendKwd("ifmove") ;
      param->appendKwd("forceOn") ;
      param->setObj(ctlBlk) ;
      param->setLongHelp(
    	"This switches on lift-and-project cuts (either at root or in entire tree). See branchAndCut for information on options."
    	) ;
      parameters.push_back(param) ;
    */

    param = new CbcGenParam(CbcGenParam::LOCALTREE,
                            "local!TreeSearch", "Whether to use local tree search", "off",
                            ctlBlk->localTree_.action_) ;
    param->appendKwd("on") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This switches on a local search algorithm when a solution is found.  This is from Fischetti and Lodi and is not really a heuristic although it can be used as one. When used from this program it has limited functionality.  It is not controlled by heuristicsOnOff."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::LOGLEVEL,
                            "log!Level", "Level of detail in cbc-generic output.",
                            -1, 999999, ctlBlk->logLvl_) ;
    param->setObj(ctlBlk) ;
    param->setPushFunc(pushCbcGenIntParam) ;
    param->setLongHelp(
        "If set to 0 then there should be no output in normal circumstances. A value of 1 is probably the best value for most uses, while 2 and 3 give more information."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::MIXEDCUTS,
                            "mixed!IntegerRoundingCuts",
                            "Whether to use Mixed Integer Rounding cuts", "off",
                            ctlBlk->mir_.action_) ;
    param->appendKwd("on") ;
    param->appendKwd("root") ;
    param->appendKwd("ifmove") ;
    param->appendKwd("forceOn") ;
    param->setObj(ctlBlk) ;
    param->setPushFunc(pushCbcGenCutParam) ;
    param->setLongHelp(
        "This switches on mixed integer rounding cuts (either at root or in entire tree).  See branchAndCut for information on options."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::USESOLUTION,
                            "force!Solution",
                            "Whether to use given solution as crash for BAB", "off", 0) ;
    param->appendKwd("on") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "If on then tries to branch to solution given by AMPL or priorities file."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::HELP, "help",
                            "Print out version, non-standard options and some help") ;
    param->setPushFunc(doHelpParam) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This prints out some help to get a user started. If you're seeing this message, you should be past that stage."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::IMPORT,
                            "import", "Import model from mps file",
                            ctlBlk->lastMpsIn_) ;
    param->setPushFunc(doImportParam) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This will read an MPS format file from the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to '', i.e., it must be set.  If you have libgz then it can read compressed files 'xxxxxxxx.gz'."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::SOLVECONTINUOUS, "initialS!olve",
                            "Solve to continuous optimum") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This just solves the problem to the continuous optimum, without adding any cuts."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::MESSAGES, "mess!ages",
                            "Controls whether standardised message prefix is printed", "off", 0) ;
    param->appendKwd("on") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "By default, messages have a standard prefix, such as:\n   Clp0005 2261  Objective 109.024 Primal infeas 944413 (758)\nbut this program turns this off to make it look more friendly.  It can be useful to turn them back on if you want to be able to 'grep' for particular messages or if you intend to override the behavior of a particular message."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::MIPLIB,
                            "miplib", "Do some of miplib test set") ;
    parameters.push_back(param) ;


    param = new CbcGenParam(CbcGenParam::OUTDUPROWS,
                            "outDup!licates",
                            "Takes duplicate rows, etc., out of the integer model", false) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::OUTPUTFORMAT,
                            "output!Format", "Which output format to use", 1, 6, 2) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "Normally export will be done using normal representation for numbers and two values per line.  You may want to do just one per line (for grep or suchlike) and you may wish to save with absolute accuracy using a coded version of the IEEE value. A value of 2 is normal. Otherwise, odd values give one value per line, even values two.  Values of 1 and 2 give normal format, 3 and 4 give greater precision, 5 and 6 give IEEE values.  When exporting a basis, 1 does not save values, 2 saves values, 3 saves with greater accuracy and 4 saves in IEEE format."
    ) ;
    parameters.push_back(param) ;

    /*
      In order for initialisation to work properly, the order of the keywords here
      must match the order of the enum IPPControl in CbcGenCtlBlk.hpp.
    */
    param = new CbcGenParam(CbcGenParam::PREPROCESS,
                            "preprocess", "Whether to use integer preprocessing", "off",
                            ctlBlk->preProcess_) ;
    param->appendKwd("on") ;
    param->appendKwd("save") ;
    param->appendKwd("equal") ;
    param->appendKwd("sos") ;
    param->appendKwd("trysos") ;
    param->appendKwd("equalall") ;
    param->appendKwd("strategy") ;
    param->setObj(ctlBlk) ;
    param->setPushFunc(pushCbcGenKwdParam) ;
    param->setLongHelp(
        "This tries to reduce size of the model in a similar way to presolve and it also tries to strengthen the model. This can be very useful and is worth trying.  save option saves on file presolved.mps.  equal will turn <= cliques into ==.  sos will create sos sets if all 0-1 in sets (well one extra is allowed) and no overlaps.  trysos is same but allows any number extra. equalall will turn all valid inequalities into equalities with integer slacks. strategy is as on but uses CbcStrategy."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::PRINTOPTIONS,
                            "pO!ptions", "Dubious print options", 0, COIN_INT_MAX, 0, false) ;
    param->setObj(ctlBlk) ;
    param->setPushFunc(pushCbcGenIntParam) ;
    param->setLongHelp(
        "If this is greater than 0 then presolve will give more information and branch and cut will give statistics"
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::PROBINGCUTS,
                            "probing!Cuts", "Whether to use Probing cuts", "off",
                            ctlBlk->probing_.action_) ;
    param->appendKwd("on") ;
    param->appendKwd("root") ;
    param->appendKwd("ifmove") ;
    param->appendKwd("forceOn") ;
    param->appendKwd("forceOnBut") ;
    param->setObj(ctlBlk) ;
    param->setPushFunc(pushCbcGenCutParam) ;
    param->setLongHelp(
        "This switches on probing cuts (either at root or in entire tree). See branchAndCut for information on options."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::INTPRINT,
                            "printi!ngOptions", "Print options", "normal", 0) ;
    param->appendKwd("integer") ;
    param->appendKwd("special") ;
    param->appendKwd("rows") ;
    param->appendKwd("all") ;
    param->setPushFunc(pushCbcGenKwdParam) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This changes the amount and format of printing a solution:\nnormal - nonzero column variables \ninteger - nonzero integer column variables\nspecial - in format suitable for OsiRowCutDebugger\nrows - nonzero column variables and row activities\nall - all column variables and row activities.\n\nFor non-integer problems 'integer' and 'special' act like 'normal'.  Also see printMask for controlling output."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::PRINTMASK,
                            "printM!ask",
                            "Control printing of solution with a regular expression", empty) ;
    param->setPushFunc(doPrintMaskParam) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "If set then only those names which match mask are printed in a solution. '?' matches any character and '*' matches any set of characters.  The default is '' (unset) so all variables are printed. This is only active if model has names."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::PRIORITYIN,
                            "prio!rityIn", "Import priorities etc from file", empty) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This will read a file with priorities from the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to '', i.e. it must be set.  This can not read from compressed files. File is in csv format with allowed headings - name, number, priority, direction, up, down, solution.  Exactly one of name and number must be given."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::REDSPLITCUTS,
                            "reduce!AndSplitCuts",
                            "Whether to use Reduce-and-Split cuts", "off",
                            ctlBlk->redSplit_.action_) ;
    param->appendKwd("on") ;
    param->appendKwd("root") ;
    param->appendKwd("ifmove") ;
    param->appendKwd("forceOn") ;
    param->setObj(ctlBlk) ;
    param->setPushFunc(pushCbcGenCutParam) ;
    param->setLongHelp(
        "This switches on reduce and split cuts (either at root or in entire tree). See branchAndCut for information on options."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::ROUNDING,
                            "round!ingHeuristic", "Whether to use Rounding heuristic", "off",
                            ctlBlk->rounding_.action_) ;
    param->appendKwd("on") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This switches on a simple (but effective) rounding heuristic at each node of tree."
    ) ;
    parameters.push_back(param) ;


    param = new CbcGenParam(CbcGenParam::EXIT, "quit", "Stops execution") ;
    param->setPushFunc(doExitParam) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This stops execution; end, exit, quit and stop are synonyms."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::DUMMY,
                            "sleep", "for debug", 0, 9999, 0, false) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "If passed to solver from ampl, then ampl will wait so that you can copy .nl file for debug."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::SOLUTION,
                            "solu!tion", "Prints solution to file",
                            std::string("stdout")) ;
    param->setPushFunc(doSolutionParam) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This will write a primitive solution file to the given file name.  It will use the default directory given by 'directory'.  A name of '$' will use the previous value for the name.  This is initialized to 'stdout'.  The amount of output can be varied using printi!ngOptions or printMask."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam ;
    CbcGenSolvers::setupSolverParam(*param) ;
    param->setObj(ctlBlk) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::SOS,
                            "sos!Options", "Whether to use SOS from AMPL", "off", 0) ;
    param->appendKwd("on") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "Normally if AMPL says there are SOS variables they should be used, but sometimes they should be turned off - this does so."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::EXIT, "stop", "Stops execution") ;
    param->setPushFunc(doExitParam) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This stops execution; end, exit, quit and stop are synonyms."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::STRENGTHEN,
                            "strengthen", "Create strengthened problem") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This creates a new problem by applying the root node cuts. All tight constraints will be in resulting problem."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::TIGHTENFACTOR, "tighten!Factor",
                            "Tighten bounds using value times largest activity at continuous solution",
                            1.0, 1.0e20) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This sleazy trick can help on some problems."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::TWOMIRCUTS,
                            "two!MirCuts",
                            "Whether to use Two phase Mixed Integer Rounding cuts", "off",
                            ctlBlk->twomir_.action_) ;
    param->appendKwd("on") ;
    param->appendKwd("root") ;
    param->appendKwd("ifmove") ;
    param->appendKwd("forceOn") ;
    param->setObj(ctlBlk) ;
    param->setPushFunc(pushCbcGenCutParam) ;
    param->setLongHelp(
        "This switches on two phase mixed integer rounding cuts (either at root or in entire tree). See branchAndCut for information on options."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::UNITTEST, "unitTest", "Do unit test") ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "This exercises the unit test."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::USERCBC,
                            "userCbc", "Hand coded Cbc stuff", 0, COIN_INT_MAX, 0, false) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "There are times (e.g., when using AMPL interface) when you may wish to do something unusual.  Look for USERCBC in main driver and modify sample code."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::VERBOSE,
                            "verbose", "Switches on longer help on single ?",
                            0, 15, ctlBlk->verbose_, false) ;
    param->setPushFunc(pushCbcGenIntParam) ;
    param->setObj(ctlBlk) ;
    param->setLongHelp(
        "Set to 1 to get short help with ? list, 2 to get long help."
    ) ;
    parameters.push_back(param) ;

    param = new CbcGenParam(CbcGenParam::SHOWUNIMP,
                            "unimp!lemented", "Report unimplemented commands.", false) ;
    param->setPushFunc(doUnimplementedParam) ;
    param->setObj(ctlBlk) ;
    parameters.push_back(param) ;

    numberParameters = parameters.size() ;
    assert (((unsigned) numberParameters) <= parameters.capacity()) ;

    return ;
}


void loadGenParamObj (const CoinParamVec paramVec, int first, int last,
                      CbcGenCtlBlk *ctlBlk)

{
    int i ;
    /*
      Load the cbc-generic object into the parameters
    */
    for (i = first ; i <= last ; i++) {
        CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(paramVec[i]) ;
        assert (genParam != 0) ;
        genParam->setObj(ctlBlk) ;
    }

    return ;
}


/* Functions to implement  cbc-generic (CbcGenParam) parameters */

/*
  Maintainer's utility, scan the parameters and report the ones that are
  unimplemented (i.e., have no pushFunc).
*/

int doUnimplementedParam (CoinParam *param)

{
    assert (param != 0) ;

    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;

    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;

    assert (ctlBlk->paramVec_ != 0) ;
    CoinParamVec &paramVec = *ctlBlk->paramVec_ ;

    int unimpCnt = 0 ;
    int maxAcross = 5 ;
    for (unsigned i = 0 ; i < paramVec.size() ; i++) {
        CoinParam *param = paramVec[i] ;
        if (param->pushFunc() == 0) {
            if (unimpCnt % maxAcross == 0) {
                std::cout << std::endl ;
            } else {
                std::cout << " " ;
            }
            std::cout << param->name() ;
            unimpCnt++ ;
        }
    }
    if (unimpCnt % maxAcross != 1) {
        std::cout << std::endl ;
    }
    std::cout << unimpCnt << " unimplemented parameters." << std::endl ;

    return (0) ;
}


/*
  Noop function. Mainly to eliminate commands from the list returned by
  doUnimplmentedParam.
*/

int doNothingParam (CoinParam *param)
{
    return (0) ;
}

/*
  Function to terminate command parsing by returning -1.
*/

int doExitParam (CoinParam *param)

{
    return (-1) ;
}

/*
  Function to print the current version.
*/

int doVersionParam (CoinParam *param)

{
    assert (param != 0) ;
    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;
    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;

    std::cout << "cbc-generic version " << ctlBlk->version_ << std::endl ;
    std::cout
        << "cbc-generic is experimental software. If you want a stable MIP "
        << "solver, please" << std::endl
        << "use cbc. If you discover bugs while using cbc-generic "
        << "please specify" << std::endl
        << "cbc-generic in the ticket description or email subject line."
        << std::endl ;

    return (0) ;
}

/*
  Function to handle help (HELP), `?' (GENERALQUERY), and `???'
  (FULLGENERALQUERY).
*/

int doHelpParam (CoinParam *param)

{
    assert (param != 0) ;
    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;
    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;

    CbcGenParam::CbcGenParamCode code = genParam->paramCode() ;

    int verbose = ctlBlk->verbose_ ;
    bool shortHelp = ((verbose & 0x01) ? true : false) ;
    bool longHelp = ((verbose & 0x02) ? true : false) ;
    bool hidden = ((verbose & 0x08) ? true : false) ;

    CoinParamVec *paramVec = ctlBlk->paramVec_ ;
    assert (paramVec != 0) ;
    /*
      Tune up the initial settings. FULLGENERALQUERY will print normally hidden
      params, and a request for long help overrules a request for short help.
    */
    if (code == CbcGenParam::FULLGENERALQUERY) {
        hidden = true ;
    }
    if (longHelp) {
        shortHelp = false ;
    }

    CoinParamUtils::printGenericHelp() ;

    std::cout << "\nAvailable commands are:" ;
    std::string pfx("  ") ;
    CoinParamUtils::printHelp(*paramVec, 0, paramVec->size() - 1, pfx,
                              shortHelp, longHelp, hidden) ;

    return (0) ;
}

/*
  Function to push a double-valued parameter.
*/

int pushCbcGenDblParam (CoinParam *param)

{
    assert (param != 0) ;
    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;
    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;

    double val = genParam->dblVal() ;
    CbcGenParam::CbcGenParamCode code = genParam->paramCode() ;

    int retval = 0 ;
    /*
      Figure out what we're doing and set the relevant field.
    */
    switch (code) {
    case CbcGenParam::DJFIX: {
        ctlBlk->djFix_.action_ = true ;
        ctlBlk->djFix_.threshold_ = val ;
        break ;
    }
    default: {
        std::cerr << "pushCbcGenDbl: no equivalent CbcGenCtlBlk field for "
                  << "parameter code `" << code << "'." << std::endl ;
        retval = -1 ;
        break ;
    }
    }

    return (retval) ;
}

/*
  Function to push an integer-valued parameter.
*/

int pushCbcGenIntParam (CoinParam *param)

{
    assert (param != 0) ;
    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;
    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;

    int val = genParam->intVal() ;
    CbcGenParam::CbcGenParamCode code = genParam->paramCode() ;

    int retval = 0 ;
    /*
      Figure out what we're doing and set the relevant field.
    */
    switch (code) {
    case CbcGenParam::CUTDEPTH: {
        ctlBlk->setCutDepth(val) ;
        break ;
    }
    case CbcGenParam::LOGLEVEL: {
        ctlBlk->setLogLevel(val) ;
        break ;
    }
    case CbcGenParam::PRINTOPTIONS: {
        ctlBlk->printOpt_ = val ;
        break ;
    }
    case CbcGenParam::VERBOSE: {
        ctlBlk->verbose_ = val ;
        break ;
    }
    default: {
        std::cerr << "pushCbcGenInt: no equivalent CbcGenCtlBlk field for "
                  << "parameter code `" << code << "'." << std::endl ;
        retval = -1 ;
        break ;
    }
    }

    return (retval) ;
}


/*
  Function to push a keyword-valued parameter. This is the catch-all function
  for keyword parameters that don't belong to any other useful grouping.
*/

int pushCbcGenKwdParam (CoinParam *param)

{
    assert (param != 0) ;
    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;
    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;

    std::string str = genParam->kwdVal() ;
    CbcGenParam::CbcGenParamCode code = genParam->paramCode() ;

    int retval = 0 ;
    /*
      Figure out what we're doing and set the relevant field.
    */
    switch (code) {
    case CbcGenParam::PREPROCESS: {
        if (str == "off") {
            ctlBlk->setIPPAction(CbcGenCtlBlk::IPPOff) ;
        } else if (str == "on") {
            ctlBlk->setIPPAction(CbcGenCtlBlk::IPPOn) ;
        } else if (str == "save") {
            ctlBlk->setIPPAction(CbcGenCtlBlk::IPPSave) ;
        } else if (str == "equal") {
            ctlBlk->setIPPAction(CbcGenCtlBlk::IPPEqual) ;
        } else if (str == "sos") {
            ctlBlk->setIPPAction(CbcGenCtlBlk::IPPSOS) ;
        } else if (str == "trysos") {
            ctlBlk->setIPPAction(CbcGenCtlBlk::IPPTrySOS) ;
        } else if (str == "equalall") {
            ctlBlk->setIPPAction(CbcGenCtlBlk::IPPEqualAll) ;
        } else if (str == "strategy") {
            ctlBlk->setIPPAction(CbcGenCtlBlk::IPPStrategy) ;
        } else {
            std::cerr
                << "pushCbcGenKwdParam(PREPROCESS): unrecognised keyword `"
                << str << "'." << std::endl ;
            retval = -1 ;
        }
        break ;
    }
    case CbcGenParam::COSTSTRATEGY: {
        if (str == "off") {
            ctlBlk->priorityAction_ = CbcGenCtlBlk::BPOff ;
        } else if (str == "priorities") {
            ctlBlk->priorityAction_ = CbcGenCtlBlk::BPCost ;
        }
        if (str == "columnOrder") {
            ctlBlk->priorityAction_ = CbcGenCtlBlk::BPOrder ;
        } else {
            std::cerr
                << "pushCbcGenKwdParam(COSTSTRATEGY): unrecognised keyword `"
                << str << "'." << std::endl ;
            retval = -1 ;
        }
        break ;
    }
    case CbcGenParam::INTPRINT: {
        if (str == "normal") {
            ctlBlk->printMode_ = 0 ;
        } else if (str == "integer") {
            ctlBlk->printMode_ = 1 ;
        } else if (str == "special") {
            ctlBlk->printMode_ = 2 ;
        } else if (str == "rows") {
            ctlBlk->printMode_ = 3 ;
        } else if (str == "all") {
            ctlBlk->printMode_ = 4 ;
        } else {
            std::cerr
                << "pushCbcGenKwdParam(INTPRINT): unrecognised keyword `"
                << str << "'." << std::endl ;
            retval = -1 ;
        }
        break ;
    }
    default: {
        std::cerr
            << "pushCbcGenKwdParam: unrecognised parameter code `"
            << code << "'." << std::endl ;
        retval = -1 ;
        break ;
    }
    }

    return (retval) ;
}

/*
  Function to push a string-valued parameter
*/

int pushCbcGenStrParam (CoinParam *param)

{
    assert (param != 0) ;
    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;
    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;

    std::string str = genParam->strVal() ;
    CbcGenParam::CbcGenParamCode code = genParam->paramCode() ;

    int retval = 0 ;
    /*
      Figure out what we're doing and set the relevant field.
    */
    switch (code) {
    case CbcGenParam::DIRECTORY: {
        char dirSep = CoinFindDirSeparator() ;
        if (str[str.length()-1] != dirSep) {
            str += dirSep ;
        }
        ctlBlk->dfltDirectory_ = str ;
        break ;
    }
    default: {
        std::cerr << "pushCbcGenStr: no equivalent CbcGenCtlBlk field for "
                  << "parameter code `" << code << "'." << std::endl ;
        retval = -1 ;
        break ;
    }
    }

    return (retval) ;
}

/*
  The various parameters to control cut generators can be
  grouped, as they all use the same set of keywords.
*/
int pushCbcGenCutParam (CoinParam *param)

{
    assert (param != 0) ;
    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;
    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;

    std::string str = genParam->kwdVal() ;
    CbcGenParam::CbcGenParamCode code = genParam->paramCode() ;
    /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we need
      to do is correct to 0 (no error) if we're successful.
    */
    int retval ;
    if (CoinParamUtils::isInteractive()) {
        retval = 1 ;
    } else {
        retval = -1 ;
    }
    /*
      First translate the keyword into the correct CGControl enum value.
    */
    CbcGenCtlBlk::CGControl action ;
    if (str == "off") {
        action = CbcGenCtlBlk::CGOff ;
    } else if (str == "on") {
        action = CbcGenCtlBlk::CGOn ;
    } else if (str == "root") {
        action = CbcGenCtlBlk::CGRoot ;
    } else if (str == "ifmove") {
        action = CbcGenCtlBlk::CGIfMove ;
    } else if (str == "forceOn") {
        action = CbcGenCtlBlk::CGForceOn ;
    } else if (str == "forceOnBut") {
        action = CbcGenCtlBlk::CGForceBut ;
    } else {
        std::cerr
            << "pushCbcGenCutParam: unrecognised keyword `"
            << str << "'." << std::endl ;
        return (retval) ;
    }
    /*
      Validate the parameter code and set a variable to separate cuts from
      heuristics.
    */
    bool isCut = false ;
    bool isHeuristic = false ;
    switch (code) {
    case CbcGenParam::CLIQUECUTS:
    case CbcGenParam::FLOWCUTS:
    case CbcGenParam::GOMORYCUTS:
    case CbcGenParam::KNAPSACKCUTS:
    case CbcGenParam::MIXEDCUTS:
    case CbcGenParam::PROBINGCUTS:
    case CbcGenParam::REDSPLITCUTS:
    case CbcGenParam::TWOMIRCUTS:
    case CbcGenParam::CUTSTRATEGY: {
        isCut = true ;
        break ;
    }
    case CbcGenParam::COMBINE:
    case CbcGenParam::FPUMP:
    case CbcGenParam::ROUNDING:
    case CbcGenParam::LOCALTREE:
    case CbcGenParam::HEURISTICSTRATEGY: {
        isHeuristic = true ;
        break ;
    }
    default: {
        std::cerr
            << "pushCbcGenCutParam: unrecognised parameter code `"
            << code << "'." << std::endl ;
        return (retval) ;
    }
    }
    /*
      See if the action is valid for the specified type. Heuristics are on or
      off; cuts can be any of the other codes. Only probing can use forceOnBut.
    */
    if (isHeuristic) {
        if (!(action == CbcGenCtlBlk::CGOff || action == CbcGenCtlBlk::CGOn)) {
            std::cerr
                << "pushCbcGenCutParam: only on or off is valid for a heuristic."
                << std::endl ;
            return (retval) ;
        }
    } else if (isCut) {
        if (action == CbcGenCtlBlk::CGForceBut &&
                code != CbcGenParam::PROBINGCUTS) {
            std::cerr
                << "pushCbcGenCutParam: forceOnBut is valid only for probing."
                << std::endl ;
            return (retval) ;
        }
    }
    /*
      We've done the basic checks; go ahead and set the relevant field in the
      control block. We shouldn't need the default case, but some compilers will
      complain if it's missing.
    */
    switch (code) {
    case CbcGenParam::CLIQUECUTS: {
        ctlBlk->setCliqueAction(action) ;
        break ;
    }
    case CbcGenParam::FLOWCUTS: {
        ctlBlk->setFlowAction(action) ;
        break ;
    }
    case CbcGenParam::GOMORYCUTS: {
        ctlBlk->setGomoryAction(action) ;
        break ;
    }
    case CbcGenParam::KNAPSACKCUTS: {
        ctlBlk->setKnapsackAction(action) ;
        break ;
    }
    case CbcGenParam::MIXEDCUTS: {
        ctlBlk->setMirAction(action) ;
        break ;
    }
    case CbcGenParam::PROBINGCUTS: {
        ctlBlk->setProbingAction(action) ;
        break ;
    }
    case CbcGenParam::REDSPLITCUTS: {
        ctlBlk->setRedSplitAction(action) ;
        break ;
    }
    case CbcGenParam::TWOMIRCUTS: {
        ctlBlk->setTwomirAction(action) ;
        break ;
    }
    case CbcGenParam::CUTSTRATEGY: {
        ctlBlk->setCliqueAction(action) ;
        ctlBlk->setFlowAction(action) ;
        ctlBlk->setGomoryAction(action) ;
        ctlBlk->setKnapsackAction(action) ;
        ctlBlk->setMirAction(action) ;
        ctlBlk->setProbingAction(action) ;
        ctlBlk->setRedSplitAction(action) ;
        ctlBlk->setTwomirAction(action) ;
        break ;
    }
    case CbcGenParam::COMBINE: {
        ctlBlk->setCombineAction(action) ;
        break ;
    }
    case CbcGenParam::FPUMP: {
        ctlBlk->setFPumpAction(action) ;
        break ;
    }
    case CbcGenParam::GREEDY: {
        ctlBlk->setGreedyCoverAction(action) ;
        ctlBlk->setGreedyEqualityAction(action) ;
        break ;
    }
    case CbcGenParam::LOCALTREE: {
        ctlBlk->setTreeLocalAction(action) ;
        break ;
    }
    case CbcGenParam::ROUNDING: {
        ctlBlk->setRoundingAction(action) ;
        break ;
    }
    case CbcGenParam::HEURISTICSTRATEGY: {
        ctlBlk->setCombineAction(action) ;
        ctlBlk->setFPumpAction(action) ;
        ctlBlk->setGreedyCoverAction(action) ;
        ctlBlk->setRoundingAction(action) ;
        ctlBlk->setTreeLocalAction(action) ;
        break ;
    }
    default: {
        std::cerr
            << "pushCbcGenCutParam: internal confusion!" << std::endl ;
        return (-1) ;
    }
    }

    return (0) ;
}

/*
  This routine imports a new constraint system into the solver.
*/

int doImportParam (CoinParam *param)

{
    assert (param != 0) ;
    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;
    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;
    /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we need
      to do is correct to 0 (no error) if we're successful.
    */
    int retval ;
    if (CoinParamUtils::isInteractive()) {
        retval = 1 ;
    } else {
        retval = -1 ;
    }
    /*
      Figure out where we're going to acquire this new model. As special cases,
      `$' says `use the previous input source' and `-' says `use stdin'.
    */
    std::string field = genParam->strVal() ;
    std::string fileName ;
    if (field == "$") {
        fileName = ctlBlk->lastMpsIn_ ;
        field = fileName ;
    } else if (field == "-") {
        fileName = "stdin" ;
        field = fileName ;
    } else {
        fileName = field ;
    }
    /*
      See if we can open a file. fileCoinReadable understands a fair bit about
      platforms and compressed files and will try variations of the file name
      (see the doxygen doc'n for details). The file name returned in field wil
      be the one that actually worked.
    */
    bool canOpen = fileCoinReadable(fileName, ctlBlk->dfltDirectory_) ;
    if (canOpen == false) {
        std::cout
            << "Unable to open file `" << fileName
            << "', original name '" << genParam->strVal() << "'." << std::endl ;
        return (retval) ;
    }
    /*
      We can find the file. Record the name. This requires a little finesse: what
      we want is the base file name (and extension(s), if present) but not the
      prefix, unless it's an absolute path.
    */
    if (!fileAbsPath(fileName)) {
        std::string::size_type pos = fileName.rfind(field) ;
        ctlBlk->lastMpsIn_ = fileName.substr(pos) ;
    } else {
        ctlBlk->lastMpsIn_ = fileName ;
    }
    /*
      Try to read the file. Standard OSI doesn't support the Clp extensions for
      keepImportNames and allowImportErrors. It should at least support
      keepImportNames. Status will be zero for a successful read.
    */
    OsiSolverInterface *lpSolver = ctlBlk->model_->solver() ;
    int status = lpSolver->readMps(fileName.c_str(), "") ;
    if (status) {
        std::cout
            << "There were " << status << " errors on input." << std::endl ;
        return (retval) ;
    }
    /*
      We have a model! Return success.
    */
    ctlBlk->goodModel_ = true ;

    return (0) ;
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

int doDebugParam (CoinParam *param)

{
    assert (param != 0) ;
    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;
    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;
    /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we need
      to do is correct to 0 (no error) if we're successful.
    */
    int retval ;
    if (CoinParamUtils::isInteractive()) {
        retval = 1 ;
    } else {
        retval = -1 ;
    }
    /*
       If the file name is `create' or `createAfterPre', we're just setting up to
       make a debug file the next time we do branch-and-cut.
    */
    std::string field = genParam->strVal() ;
    if (field == "create" || field == "createAfterPre") {
        ctlBlk->debugCreate_ = field ;
        return (0) ;
    }
    /*
      Figure out where we're going to acquire the debug file. As special cases,
      `$' says `use the previous input source' and `-' says `use stdin'.
    */
    std::string fileName ;
    if (field == "$") {
        fileName = ctlBlk->debugFile_ ;
        field = fileName ;
    } else if (field == "-") {
        fileName = "stdin" ;
        field = fileName ;
    } else {
        fileName = field ;
    }
    /*
      See if we can open a file. fileCoinReadable understands a fair bit about
      platforms and compressed files and will try variations of the file name
      (see the doxygen doc'n for details). The file name returned in field wil be
      the one that actually worked. No default prefix --- a debug file is assumed
      to always be in the current directory.
    */
    bool canOpen = fileCoinReadable(fileName, ctlBlk->dfltDirectory_) ;
    if (canOpen == false) {
        std::cout
            << "Unable to open file `" << fileName
            << "', original name '" << genParam->strVal() << "'." << std::endl ;
        return (retval) ;
    }
    /*
      We can find the file. Record the name. This requires a little finesse: what
      we want is the base file name (and extension(s), if present) but not the
      prefix, unless it's an absolute path.
    */
    if (!fileAbsPath(fileName)) {
        std::string::size_type pos = fileName.rfind(field) ;
        ctlBlk->lastMpsIn_ = fileName.substr(pos) ;
    } else {
        ctlBlk->lastMpsIn_ = fileName ;
    }
    /*
      Load the primal variable values into the debug solution vector.
    */
    int intUnused, numCols ;
    double dblUnused ;
    double *primals ;

    bool readOK = readSolution(fileName,
                               intUnused, numCols, dblUnused, 0, 0, &primals, 0) ;

    if (readOK) {
        if (ctlBlk->debugSol_.values_) {
            delete[] ctlBlk->debugSol_.values_ ;
        }
        ctlBlk->debugSol_.numCols_ = numCols ;
        ctlBlk->debugSol_.values_ = primals ;
        retval = 0 ;
    } else {
        if (primals) {
            delete[] primals ;
        }
    }

    return (retval) ;
}


/*
  Utility routine to save the current solution to a file. No formatting, and
  not intended to be portable in any way, shape, or form.
*/

void saveSolution (const OsiSolverInterface *osi, std::string fileName)

{
    FILE *fp = fopen(fileName.c_str(), "wb") ;

    if (fp) {
        int numberRows = osi->getNumRows() ;
        int numberColumns = osi->getNumCols() ;
        double objectiveValue = osi->getObjValue() ;

        fwrite(&numberRows, sizeof(int), 1, fp) ;
        fwrite(&numberColumns, sizeof(int), 1, fp) ;
        fwrite(&objectiveValue, sizeof(double), 1, fp) ;

        const double *primalRowSolution = osi->getRowActivity() ;
        const double *dualRowSolution = osi->getRowPrice() ;
        const double *primalColumnSolution = osi->getColSolution() ;
        const double *dualColumnSolution = osi->getReducedCost() ;

        fwrite(primalRowSolution, sizeof(double), numberRows, fp) ;
        fwrite(dualRowSolution, sizeof(double), numberRows, fp) ;
        fwrite(primalColumnSolution, sizeof(double), numberColumns, fp) ;
        fwrite(dualColumnSolution, sizeof(double), numberColumns, fp) ;

        fclose(fp) ;
    } else {
        std::cout
            << "saveSolution: Unable to open file `"
            << fileName << "'." << std::endl ;
    }

    return ;
}

/*
  Utility routine to read in a solution dump created by saveSolution. Generally
  we don't need all the info in this file, so the routine accepts a bunch of
  reference/pointer paramaters and fills in any that are non-null. It's the
  client's responsibility to dispose of space allocated for solution vectors.
  The parameters fileName, numRows, numCols, and objVal are mandatory. The rest
  can be null.
*/
bool readSolution (std::string fileName,
                   int &numRows, int &numCols, double &objVal,
                   double **rowActivity, double **dualVars,
                   double **primalVars, double **reducedCosts)

{
    FILE *fp = fopen(fileName.c_str(), "rb") ;
    bool retval = true ;

    numRows = -1 ;
    numCols = -1 ;
    objVal = 0 ;
    *rowActivity = 0 ;
    *dualVars = 0 ;
    *primalVars = 0 ;
    *reducedCosts = 0 ;

    if (fp) {
        fread(&numRows, sizeof(int), 1, fp) ;
        fread(&numCols, sizeof(int), 1, fp) ;
        fread(&objVal, sizeof(double), 1, fp) ;

        if (rowActivity != NULL) {
            *rowActivity = new double [numRows] ;
            fread(*rowActivity, sizeof(double), numRows, fp) ;
        } else {
            fseek(fp, numRows*sizeof(double), SEEK_CUR) ;
        }

        if (dualVars != NULL) {
            *dualVars = new double [numRows] ;
            fread(*dualVars, sizeof(double), numRows, fp) ;
        } else {
            fseek(fp, numRows*sizeof(double), SEEK_CUR) ;
        }

        if (primalVars != NULL) {
            *primalVars = new double [numCols] ;
            fread(*primalVars, sizeof(double), numCols, fp) ;
        } else {
            fseek(fp, numCols*sizeof(double), SEEK_CUR) ;
        }

        if (reducedCosts != NULL) {
            *reducedCosts = new double [numCols] ;
            fread(*reducedCosts, sizeof(double), numCols, fp) ;
        } else {
            fseek(fp, numCols*sizeof(double), SEEK_CUR) ;
        }

        fclose(fp) ;
    } else {
        std::cout
            << "readSolution: Unable to open file `"
            << fileName << "'." << std::endl ;
        retval = false ;
    }

    return (retval) ;
}


} // end namespace CbcGenParamUtils


