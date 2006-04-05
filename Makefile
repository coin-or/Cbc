# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.Cbc

# Makefile.location specifies the components of COIN that are locally
# available. Makefile.Cbc is where you specify compiler optimisation level and
# library type.

# From this makefile you can build three primary targets:
# * libCbc, the cbc branch-and-cut library
# * solve, a main program which is tailored to work with the clp solver
# * cbc, a main program which is capable of using solvers other than clp

# libCbc generally uses the OSI interface to consult the underlying solver. It
# contains some code that is specific to OsiClp, but it will work with any OSI
# solver that supports the required functionality. Using solvers other than clp
# should be considered a work in progress. If all you want to do is solve a few
# MIP problems, you probably want to leave CBC_SOLVERS = Clp (the default) and
# build the `solve' main program.

# If you want to experiment with other solvers, adjust CBC_SOLVERS below and
# build the `cbc' main program.  Be sure that the solvers are available and
# that the information in Makefile.location is correct for all specified
# solvers. As far as libCbc is concerned, the only thing that matters is
# presence or absence of clp. The `cbc' main program will create a set of OSI
# solver prototypes to match the list in CBC_SOLVERS. You can change solvers at
# run time.  The command `cbc -miplib' will exercise all specified solvers on
# the miplib examples. The -miplib option uses OsiCbc. If you want to
# completely exclude clp from the build, you should also examine and edit
# Osi/OsiCbc/Makefile.

# To specify a solver, add the name to the list here. Use the name given in
# Makefile.location, without the COIN_lib prefix, e.g., Clp, Cpx, Dylp, Glpk,
# etc. The default build uses Clp only. To use Clp and Dylp, you would say
# CBC_SOLVERS := Clp Dylp

CBC_SOLVERS := Clp

# Regardless of the number of solvers specified, it's a good idea to set the
# default solver. (Check the code in CbcMain if you must know what will happen
# if you don't.) This must match one of the solver names given above, but all
# in lower case letters (don't ask).

cbcDefaultSolver := clp

# You should not need to make changes below this line.
###############################################################################

# Bring on the boilerplate. Makefile.coin will bring in Makefile.<O/S> and
# Makefile.location.

export CoinDir := $(shell cd ..; pwd)
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin

# First check that all solvers specified by the user have been configured in
# Makefile.location.

cbcMissingSolvers := $(filter-out \
	$(patsubst COIN_libOsi%,%,$(CoinLibsDefined)),$(CBC_SOLVERS))

ifneq ($(cbcMissingSolvers),)
  $(foreach solver,$(cbcMissingSolvers), \
    $(warning $(solver) is not configured in Makefiles/Makefile.location.) \
    $(warning Probably the line 'CoinLibsDefined += COIN_lib$(solver)' is \
      commented out))
  $(error Please correct Makefile.location and try again.)
endif

$(info Building $(MAKECMDGOALS) with solvers $(CBC_SOLVERS))

# Generate appropriate defines for the compilation command from the value of
# CBC_SOLVERS and cbcDefaultSolver.

CBC_DEFINES := $(foreach solver,$(CBC_SOLVERS), \
		 $(patsubst COIN_HAS_%,CBC_USE_%,\
		   $(filter COIN_HAS_%,$($(solver)Define))))
CBC_DEFINES += CBC_DEFAULT_SOLVER="\"$(cbcDefaultSolver)\""

export CBC_DEFINES
export CBC_SOLVERS

# $(warning CBC_DEFINES is $(CBC_DEFINES))

# Pull together a dependency list for libCbc. You can't build libCbc without
# the Coin, Osi, and Cgl libraries.  The current (06.03.22) Osi makefile will
# use the set of solvers enabled in Makefile.location to determine which OSI
# interfaces to build. To make a long story short, since you can't add a solver
# to CBC_SOLVERS unless it's enabled in Makefile.location, we don't need to
# explicitly trigger a build for solver XXX, or for interface OsiXXX. Nor can
# we suppress the build of OsiCbc. (But note that disabling OsiCbc in
# Makefile.location will exclude it from the link of the cbc main program.) If
# this behaviour ever changes, the target list here will need to be upgraded to
# include XXX and OsiXXX for solvers XXX in CBC_SOLVERS.

libTgts := Coin Osi Cgl

# $(warning Dependency list is $(libTgts))

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install clean library unitTest cbc solver solve \
	libdepend libCbc doc

# These targets are for libCbc.

default: install

libCbc: library

install library: libdepend
	${MAKE} -f Makefile.Cbc $@

libdepend:
	$(foreach tgt,$(libTgts),(cd $(CoinDir)/$(tgt) && $(MAKE) install) ; )

# These targets are for the two main programs, cbc and solve. In each case,
# we'll do libCbc first, then go for the main program.

unitTest cbc: install
	(cd Test && ${MAKE} cbc)

solver solve: install
	(cd Test && ${MAKE} solver)

clean: 
	@rm -rf Junk
	@rm -rf $(UNAME)*
	@rm -rf dep
	@rm -rf Test/Junk
	@rm -rf Test/$(UNAME)*
	@rm -rf Test/dep
	@rm -f cbc
	@rm -f solve

doc:
	doxygen $(MakefileDir)/doxygen.conf
