# CBC

Cbc (*C*oin-or *b*ranch and *c*ut) is an open-source mixed integer linear programming solver written in C++.
It can be used as a callable library or using a stand-alone executable.
It can be called through
AIMMS (through the [AIMMSlinks](https://github.com/coin-o/AIMMSlinks) project),
AMPL (natively),
[CMPL](https://github.com/coin-or/Cmpl),
GAMS (through the [GAMSlinks](https://github.com/coin-or/GAMSlinks) project),
[JuMP](https://github.com/JuliaOpt/JuMP.jl),
Mathematica,
[MiniZinc](http://www.minizinc.org/),
MPL (through the [CoinMP](https://github.com/coin-or/CoinMP) project),
[PuLP](https://github.com/coin-or/PuLP),
Python (e.g., [cbcpy](https://pypi.org/project/cbcpy)), and
[OpenSolver for Excel](http://opensolver.org), among others.

Cbc links to a number of other COIN-OR projects for additional functionality, including:
 * [Clp](https://github.com/coin-or/Clp) (the default solver for LP relaxations)
 * [Cgl](https://github.com/coin-or/Cgl) (for cut generation)
 * [CoinUtils](https://github.com/coin-or/CoinUtils) (for reading input files and various utilities)
For more information on supported platforms, links to dependent projects, current version, and more, click [here](http://www.coin-or.org/projects/Cbc.xml)

Cbc is developed by [John Forrest](http://www.fastercoin.com), now retired from IBM Research.
The project is currently managed by John Forrest, [Ted Ralphs](http://coral.ie.lehigh.edu/~ted/), Haroldo Gambini Santos, and the rest of the Cbc team (Dan Fylstra (Frontline), Lou Hafer (SFU), Bill Hart (Sandia), Bjarni Kristjannson (Maximal), Cindy Phillips (Sandia), Matthew Saltzman (Clemson), Edwin Straver (Frontline), Jean-Paul Watson (Sandia)).

Cbc is written in C++ and is released as open source code under the [Eclipse Public License (EPL)](http://www.opensource.org/licenses/eclipse-1.0) and is freely redistributable.
All source code and documentation is Copyright IBM and others. This README may be redistributed freely.

Cbc is available from the [COIN-OR initiative](http://www.coin-or.org/).
The Cbc website is https://github.com/coin-or/Cbc.

## CITE

[![DOI](https://zenodo.org/badge/30382416.svg)](https://zenodo.org/badge/latestdoi/30382416)

## CURRENT BUILD STATUS

[![Build Status](https://travis-ci.org/coin-or/Cbc.svg?branch=master)](https://travis-ci.org/coin-or/Cbc)

[![Build status](https://ci.appveyor.com/api/projects/status/l2hwifsxwhswng8y/branch/master?svg=true)](https://ci.appveyor.com/project/tkralphs/cbc/branch/master)


## DOWNLOAD

Binaries for most platforms are available for download from [Bintray](https://bintray.com/coin-or/download/Cbc)

[ ![Download](https://api.bintray.com/packages/coin-or/download/Cbc/images/download.svg?version=2.10) ](https://bintray.com/coin-or/download/Cbc/2.10/link)

 * *Linux*: On Debian/Ubuntu, Cbc is available in the package `coinor-cbc` and can be installed with apt. On Fedora, Cbc is available in the package `coin-or-Cbc`.
 * *Windows*: The easiest way to get Cbc on Windows is to download from *[Bintray](https://bintray.com/coin-or/download/Cbc)*, although an old interactive installer for the [COIN-OR Optimization Suite](http://www.coin-or.org/download/binary/CoinAll) is also still available.
 * *Mac OS X*: The easiest way to get Cbc on Mac OS X is through [Homebrew](https://brew.sh).
   * `brew tap coin-or-tools/coinor`
   * `brew install cbc`
 * AMPL also provides stand-alone [Cbc executables](http://ampl.com/products/solvers/open-source/#cbc) that can be used with (or without) AMPL.

Due to license incompatibilities, pre-compiled binaries lack some functionality.
If binaries are not available for your platform for the latest version and you would like to request them to be built and posted, feel free to let us know on the mailing list.

*Source code* can be obtained either by

 * Downloading a snapshot of the source code for the latest release version of Cbc from the [Cbc source code download page](http://www.coin-or.org/download/source/Cbc), or
 * Checking out the code from [Github](https://github.com/coin-or/Cbc) or using the `coinbrew` script (recommended). 

Below is a quick start guide for building on common platforms. More detailed
build instructions are
[https://coin-or.github.io/user_introduction.html](here) (this is a work in
progress).

## BUILDING from source

### Using CoinBrew

To build Cbc from source, obtain the `coinbrew` script from
https://coin-or.github.io/coinbrew/
and run

    /path/to/coinbrew fetch Cbc
    /path/to/coinbrew build =Cbc --prefix=/dir/to/install --test
    /path/to/coinbrew install Cbc

The `coinbrew` script will fetch [these](Dependencies) additional projects.

### Without CoinBrew (Expert users)

 0. Install [these Dependencies](Dependencies)
 1. Obtain the source code, e.g., from https://github.com/coin-or/Cbc
 2. Run `./configure -C` to generate makefiles
 3. Run `make` to build the CoinUtils library
 4. Run `make test` to build and run the CoinUtils unit test program
 5. Run `make install` to install library and header files.

### With Microsoft Visual Studio

For Microsoft Visual C++ users, there are project files for version 10
available in the `MSVisualStudio` directory. First, obtain the source code
using either a Windows git client or download a snapshot. In MSVC++ Version
10, open the solution file (this should be converted to whatever version of
MSVC+ you are using) and build the Cbc project. The code should build out of
the box with default settings.

It is also possible to build Cbc with the Visual Studio compiler from the
command line using the procedure for Unix-like environments, using the Msys2
shell or CYGWIN. This is the recommended and best-supported way of building
Cbc in Windows from source.

If you want to build a *parallel version* of CBC using Visual Studio you can
following instructions: (thanks to Tobias Stengel and Alexis Guigue).

Assumptions:

- A VS solution with all necessary projects (libCbc, libClp, libCbcSolver,
  libCgl, libCoinUtils, libOsi, libOsiCbc, libOsiClp). The project files can
  be found inside the "MSVisualStudio" folders. 

Steps (based on VS 2013):

1. for each of the lib* projects do:
   add `CBC_THREAD` under Properties -> Configuration Properties -> C/C++ ->
   Preprocessor -> Preprocessor Definitions (a subset of the lib* projects may
   be sufficient, but it doesn't hurt to do it for all) 

2. Link against a pthreads library.
   [PThreadsWin32](https://www.sourceware.org/pthreads-win32/) works (even in
   64 bits systems) and is distributed under the LGPL. If you decide to use
   the precompiled binaries: both pthreadVC2 and pthreadVS2 seem to work.
   Otherwise: third party VS project files for pthreads4win can be found on
   github.

   Note: If you use C++/Cli, make sure that no callback (eventHandlers, custom
   cut generators, custom heuristics, ...) contains managed code. Otherwise
   your code will crash at runtime with AssembyNotFoundExceptions inside the
   native threads created by Cbc. Even if not, problems with the GC are
   likely.

3. If you link statically against pthreads4win, you have to define
   PTW32_STATIC_LIB when building your program/Cbc (e.g. via Properties ->
   C/C++/Preprocessor -> Preprocessor Definitions) AND - only if you build
   pthreads yourself - when building pthreads. Linking pthreads dynamically
   works without additional preprocessor definitions.

4. pass "-threads" "yourNumber" to CbcMain1


## DOCUMENTATION

 * [INSTALL](INSTALL) file (partially outdated)
 * [User's Guide](https://coin-or.github.io/Cbc) (from 2005)
 * [Doxygen generated documentation](http://www.coin-or.org/Doxygen/Cbc/hierarchy.html)
 * Source code [examples](Cbc/examples)
 * [Cbc command-line guide](https://projects.coin-or.org/CoinBinary/export/1059/OptimizationSuite/trunk/Installer/files/doc/cbcCommandLine.pdf)


## SUPPORT

### List Serve

CBC users should use the Cbc mailing list. To subscribe, go to
http://list.coin-or.org/mailman/listinfo/cbc

### Bug Reports

Bug reports should be reported on the CBC development web site at

https://github.com/coin-or/Cbc/issues

## CHANGELOG

 * Release 2.10.4
   * Fix parsing of optional arguments to AMPL interface.

 * Release 2.10.3
   * Improve performance of some primal heuristics, incl. feasibility pump, by making integer slacks continuous
   * Added additional timelimit checks
   * Fixed initialization of Cbc_clone result
   * Additional bugfixes

 * Release 2.10.2
   * Bugfixes

 * Release 2.10.1
   * Fixed Cbc_clone in C interface
   * Fixed CbcMain1() call in examples/driver3.cpp
   * Fixed possible issue with MIPstart if presolve added variables
   * More minor bugfixes

 * Release 2.10.0
   * Improved handling of SOS, starting point, and symmetries
   * Improved performance of primal heuristics regarding the handling of
     implicit integer variables
   * Mini-B&B is now disabled when solving with multiple threads
   * Changed default value for zero half cuts parameter from off to ifmove
   * Added CbcModel::postProcessedSolver() to obtained LP after presolve
   * New option "PrepNames" to indicate whether column names should be
     kept in the pre-processed model
   * New option "sosPrioritize" to determine how to prioritize SOS
   * Added new event "generatedCuts"
   * CbcSolver can now read compressed .lp files (GZIP, BZIP2)
   * New functions in the C interface: Cbc_readLp, Cbc_writeLp,
     Cbc_addCol, Cbc_addRow, Cbc_getNumIntegers, Cbc_bestSolution,
     Cbc_getObjValue, Cbc_getRowNz, Cbc_getRowIndices, Cbc_getRowCoeffs,
     Cbc_getRowRHS, Cbc_getRowSense, Cbc_getColNz, Cbc_getColIndices,
     Cbc_getColCoeffs, Cbc_getReducedCost, Cbc_numberSavedSolutions,
     Cbc_savedSolution, Cbc_savedSolutionObj, Cbc_setMIPStart,
     Cbc_setMIPStartI, Cbc_addCutCallback, Osi_getNumCols, Osi_getColName,
     Osi_getColLower, Osi_getColUpper, Osi_isInteger, Osi_getNumRows,
     Osi_getRowNz, Osi_getRowIndices, Osi_getRowCoeffs, Osi_getRowRHS,
     Osi_getRowSense, Osi_getColSolution, OsiCuts_addRowCut,
     Cbc_getAllowableGap, Cbc_setAllowableGap, Cbc_getAllowableFractionGap,
     Cbc_setAllowableFractionGap, Cbc_getAllowablePercentageGap,
     Cbc_setAllowablePercentageGap, Cbc_getCutoff, Cbc_setCutoff,
     Cbc_getMaximumNodes, Cbc_setMaximumNodes, Cbc_getMaximumSolutions,
     Cbc_setMaximumSolutions, Cbc_getLogLevel, Cbc_setLogLevel,
     Cbc_getMaximumSeconds, Cbc_setMaximumSeconds
   * New action "guess" checks properties of the model to decide the best
     parameters for solving the LP relaxation.
   * New example inc.cpp to illustrate solution callback
   * New example driver5.cpp to illustrate user-defined branching rule
   * New example clpdriver.cpp to illustrate use of ClpEventHandler
   * Added support for using OsiHiGHS with CbcGeneric
   * Added MSVC 14 project files
   * Bugfixes

 * Release 2.9.10
   * Fix a numerical issue
   * Fix some memory leaks
   * Fix issue when root node is obviously infeasible
   * Performance improvements for mini-B&B
   * Fix name of bound in final message
   * Fix names in preprocessed problem

 * Release 2.9.9

   * Fixes for SOS2
   * Updates to mipstart
   * Switching to new build system
   * Updates for CI

 * Release 2.9.8

   * Update to most current releases of dependencies
   * Small bug fixes
   * Add support for automatic build and test with Travis and Appveyor

 * Release 2.9.7

   * Small bug fixes
   * Option to switch to line buffered output

 * Release 2.9.6

   * Small bug fixes

 * Release 2.9.5

   * Small bug fixes

 * Release 2.9.4

   * Small fixes for stability
   * Fixes for Doygen documentation generation

 * Release 2.9.3

   * Minor bug fixes

 * Release 2.9.2

   * Fix for proper installation with ```DESTDIR```

 * Release 2.9.1

   * Fix for dependency linking
   * Minor bug fixes

 * Release 2.9.0

   * Introduced specialized branching methods for dealing with "big Ms".
   * Introduced new methods for dealing with symmetry (requires installation of [nauty](http://pallini.di.uniroma1.it/))
   * Introduction of conflict cuts (off by default, turn on with `-constraint conflict`)

 * Release 2.8.13

   * Improved message handling
   * Miscellaneous bug fixes.

 * Release 2.8.12

   * Update for dependencies.

 * Release 2.8.11

   * Major overhaul of C interface
   * Fixes to SOS
   * Miscellaneous bug fixes

 * Release 2.8.10

   * More changes related to thread safety.
   * Fix bug in build system with Visual Studio compiler.
   * Miscellaneous bug fixes.

 * Release 2.8.9

   * Attempt to make Cbc thread safe.
   * Add parallel examples.
   * Add CbcSolverUsefulInfo.
   * Bug fixes.

 * Release 2.8.8

   * Added example to show how to use Cbc with installed libraries in MSVC++
   * Fixed inconsistency in addition of libCbcSolver to dependencies in
     {{{cbc_addlibs.txt}}}.

 * Release 2.8.7

   * Changed so that Doxygen builds LaTex
   * Fixes for build system

 * Release 2.8.6

   * Added option to explicitly link dependencies to comply with packaging
     requirements on Fedora and Debian, as well as allow building of MinGW
     DLLs.

 * Release 2.8.5

   * Minor fixes to build system

 * Release 2.8.4

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
