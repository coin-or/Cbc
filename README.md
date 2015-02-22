# CBC Version 2.9.2 README

Welcome to the README for the COIN Branch and Cut Solver (CBC). CBC is
distributed under the Eclipse Public License and is freely redistributable.
All source code and documentation is Copyright IBM and others. This README may
be redistributed freely.

## DOCUMENTATION

For a quick start guide, please see the INSTALL file in this distribution. A
(somehwat outdated) user's manual is available here:

http://www.coin-or.org/Cbc

More up-to-date automatically generated documentation of the source code can
be found here:

http://www.coin-or.org/Doxygen/Cbc/

Further information can be found here:

http://projects.coin-or.org/Cbc

## SUPPORT

### List Serve

CBC users should use the Cbc mailing list. To subscribe, go to 
http://list.coin-or.org/mailman/listinfo/cbc

### Bug Reports

Bug reports should be reported on the CBC development web site at

https://projects.coin-or.org/Cbc/newticket

## CHANGELOG

 * Release 5.9.2

   * Fix for proper installation with {{{DESTDIR}}}

 * Release 5.9.1

   * Fix for dependency linking
   * Minor bug fixes

 * Release 5.9.0 
   * Major algorithmic improvements

 * Release 5.8.13
   * Improved message handling
   * Miscellaneous bug fixes.

 * Release 5.8.12
   * Update for dependencies.

 * Release 5.8.11
   * Major overhaul of C interface
   * Fixes to SOS
   * Miscellaneous bug fixes

 * Release 5.8.10
   * More changes related to thread safety.
   * Fix bug in build system with Visual Studio compiler.
   * Miscellaneous bug fixes.

 * Release 5.8.9
   * Attempt to make Cbc thread safe.
   * Add parallel examples.
   * Add CbcSolverUsefulInfo.
   * Bug fixes.

 * Release 5.8.8

   * Added example to show how to use Cbc with installed libraries in MSVC++
   * Fixed inconsistency in addition of libCbcSolver to dependencies in
     {{{cbc_addlibs.txt}}}.

 * Release 5.8.7

   * Changed so that Doxygen builds LaTex
   * Fixes for build system

 * Release 5.8.6

   * Added option to explicitly link dependencies to comply with packaging
     requirements on Fedora and Debian, as well as allow building of MinGW
     DLLs.

 * Release 5.8.5
   * Minor fixes to build system

 * Release 5.8.4
   * Small bug fixes
   * Upgrades to build system

 * Release 2.8.3:

   * Fix for handling SOS.

 * Release 2.8.2:

   * Fixed recognition of Glpk source in main configure.
   * Minor bug fixes in CoinUtils, Clp, and Cbc.

 * Release 2.8.1:

   * Minor bug fixes

 * Release 2.8.0:

   * Introduced new secondaryStatus 8 to indicate that solving stopped due to
     an iteration limit.
   * Solution pool is now accessible via the command line and the CbcMain*
     interface.
   * New mipstart option to read an initial feasible solution from a file.
     Only values for discrete variables need to be provided.

   * Added Proximity Search heuristic by Fischetti and Monaci (off by
     default): The simplest way to switch it on using stand-alone version is
     ```-proximity on```.

     Proximity Search is the new "No-Neighborhood Search" 0-1 MIP refinement
     heuristic recently proposed by Fischetti and Monaci (2012). The idea is
     to define a sub-MIP without additional constraints but with a modified
     objective function intended to attract the search in the proximity of the
     incumbent. The approach works well for 0-1 MIPs whose solution landscape
     is not too irregular (meaning the there is reasonable probability of
     finding an improved solution by flipping a small number of binary
     variables), in particular when it is applied to the first heuristic
     solutions found at the root node.
   
   * An implementation of Zero-Half-Cuts by Alberto Caprara is now available.
     By default, these cuts are off. To use add to your command line
     -zerohalfCuts root (or other options) or just -zero. So far, they may
     help only on a small subset of problems and may need some tuning.

     The implementation of these cuts is described in G. Andreello, A.
     Caprara, and M. Fischetti "Embedding Cuts in a Branch and Cut Framework:
     a Computational Study with {0,1/2}-Cuts" INFORMS Journal on Computing
     19(2), 229-238, 2007 http://dx.doi.org/10.1287/ijoc.1050.0162

   * An alternative implementation of a reduce and split cut generator by
     Giacomo Nannicini is now available. By default, these cuts are off. To
     use add to your command line -reduce2AndSplitCuts root (or other
     options).

     The implementation of these cuts is described in G. Cornuejols and G.
     Nannicini "Practical strategies for generating rank-1 split cuts in
     mixed-integer linear programming" Mathematical Programming Computation
     3(4), 281-318, 2011 http://dx.doi.org/10.1007/s12532-011-0028-6

   * An alternative robust implementation of a Gomory cut generator by Giacomo
     Nannicini is now available. By default, these cuts are off. To use add to
     your command line -GMI root (or other options).

     The implementation of these cuts is described in G. Cornuejols, F.
     Margot, and G. Nannicini "On the safety of Gomory cut generators"
     http://faculty.sutd.edu.sg/~nannicini/index.php?page=publications

   * To encourage the use of some of the more exotic/expensive cut generators
     a parameter -slowcutpasses has been added. The idea is that the code does
     these cuts just a few times - less than the more usual cuts. The default
     is 10. The cut generators identified by "may be slow" at present are just
     Lift and project and ReduceAndSplit (both versions).
   
   * Allow initialization of random seed by user. Pseudo-random numbers are
     used in Cbc and Clp. In Clp they are used to break ties in degenerate
     problems, while in Cbc heuristics such as the Feasibility Pump use them
     to decide whether to round up or down. So if a different pseudo-random
     seed is given to Clp then you may get a different continuous optimum and
     so different cuts and heuristic solutions. This can be switched on by
     setting randomSeed for Clp and/or randomCbcSeed for Cbc. The special
     value of 0 tells code to use time of day for initial seed.
   
   * Building on this idea, Andrea Lodi, Matteo Fischetti, Michele Monaci,
     Domenico Salvagnin, Yuji Shinano, and Andrea Tramontani suggest that this
     idea be improved by running at the root node with multiple copies of
     solver, each with its own different seed and then passing in the
     solutions and cuts so that the main solver has a richer set of solutions
     and possibly stronger cuts. This is switched on by setting
     -multipleRootPasses. These can also be done in parallel.
 
   * Few changes to presolve for special variables and badly scaled problems
     (in CoinUtils).

   * New option -extraVariables <number> which switches on a trivial
    re-formulation that introduces extra integer variables to group together
    variables with same cost.

   * For some problems, cut generators and general branching work better if
     the problem would be infeasible if the cost is too high. If the new
     option -constraintFromCutoff is set, the objective function is added as a
     constraint which rhs is set to the current cutoff value (objective value
     of best known solution).

 * Release 2.7.8:

   * Change message when LP simplex iteration limit is hit from "Exiting on
     maximum nodes" to "Exiting on maximum number of iterations"
   * Fix for using overlapping SOS.
   * Fixes in buildsystem.

 * Release 2.7.7:

   * Fix to report interruption on user event if SIGINT is received by
     CbcSolver. model->status() should now be 5 if this event happened. Added
     method CbcModel::sayEventHappened() to make cbc stop due to an 'user
     event'.

   * Other minor fixes.

 * Release 2.7.6:

   * Fixes to build system.

   * Other minor fixes.

 * Release 2.7.5:

   * Fixes to get AMPL interface working again.

   * More fixes to MSVC++ files.

 * Release 2.7.4:

   * Minor bugfixes.

 * Release 2.7.3:

   * Minor bugfixes.

   * Fixes to MSVC++ files.

 * Release 2.7.2:

   * Allow row/column names for GMPL models.

   * Added CbcModel::haveMultiThreadSupport() to indicate whether Cbc library
     has been compiled with multithread support.

   * Added CbcModel::waitingForMiniBranchAndBound() to indicate whether
     sub-MIP heuristic is currently running.

   * Cbc shell should work with readline if configured with
     ```--enable-gnu-packages```.

   * Support for compressed input files (.gz, .bz2) is now enabled by default.

   * Fix problems with relative gap tolerance > 100% and further bugs.

   * Fixes for MSVC++ Version 9 files.

   * Minor fixes in buildsystem; update to BuildTools 0.7.1.

 * Release 2.7.1:

   * Fixes to MSVC++ files

 * Release 2.7.0:

   * License has been changed to the EPL.

   * Support for MSVC++ version 10 added.

   * Support for BuildTools version 0.7 to incorporate recent enhancements,
     including proper library versioning in Linux, prohibiting installation of
     private headers, etc.

   * Updated externals to new stable versions of dependent projects.

   * Improvements to heuristics. 

   * New options for cut generation. 

   * Improved reporting of results. 

   * Improvements to documentation. 

   * Minor bug fixes. 


