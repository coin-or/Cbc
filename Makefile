# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.Cbc

# Define variables here to select the solvers that you want to incorporate into
# the cbc build. Be sure that the solvers are available and that the
# information in Makefile.location is correct for all selected solvers.
# Check Makefile.Cbc to specify the optimisation level of the build.

# Compile-time configuration for Cbc. There are basically two options:
# (1)
# Build cbc to use only clp as the underlying lp solver. This is the most
# efficient configuration if you wish to take full advantage of capabilities
# supported by clp but not available through the standard OSI interface, and
# you have no interest in experimenting with other solvers.
# (2)
# Build cbc so that it can use some mix of one or more OSI solvers, not
# necessarily including clp. Use this option if you want to experiment with
# other solvers.

# To specify a solver, add the name to the list here. Use the name given in
# Makefile.location, without the COIN_lib prefix, e.g., Clp, Cpx, Dylp, Glpk,
# etc. The default build uses Clp only. To use Clp and Dylp, you would say
# CBC_SOLVERS := Clp Dylp

CBC_SOLVERS := Clp Dylp

# Regardless of the number of solvers specified, it's a good idea to set the
# default solver. (All right, the real reason is it'll take too long to explain
# what happens if you don't. Check the code in CbcMain if you must know.)
# This must match one of the solver names given above, but all in lower case
# letters (don't ask).

cbcDefaultSolver := dylp

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

# Figure out the configuration based on the value of CBC_SOLVERS. We need to
# generate appropriate defines for the compilation command.

CBC_DEFINES :=

ifeq ($(CBC_SOLVERS),Clp)
  CBC_ONLY_CLP := 1
  CBC_DEFINES := CBC_ONLY_CLP CBC_USE_CLP
else
  CBC_ONLY_CLP := 0
  CBC_DEFINES := $(foreach solver,$(CBC_SOLVERS), \
      $(patsubst COIN_HAS_%,CBC_USE_%,\
	$(filter COIN_HAS_%,$($(solver)Define))))
endif
CBC_DEFINES += CBC_DEFAULT_SOLVER="\"$(cbcDefaultSolver)\""
export CBC_ONLY_CLP
export CBC_DEFINES
export CBC_SOLVERS

# $(warning CBC_DEFINES is $(CBC_DEFINES))

$(warning Building cbc with solvers $(CBC_SOLVERS))

# Pull together the full dependency list for cbc. You can't build cbc without
# the Coin, Osi, and Cgl libraries. Add Coin and Osi later, for technical
# reasons.

libTgts := Cgl

# This makefile fronts for two main programs down in the Test directory, cbc
# and solve. solve will always want OsiClp and Vol (note: not OsiVol). Add them
# here, if they're not already in CBC_SOLVERS, and add the list from
# CBC_SOLVERS.

libTgts += $(filter-out $(CBC_SOLVERS),Vol)
libTgts += $(patsubst %,Osi%,$(CBC_SOLVERS))
libTgts += $(filter-out $(libTgts),OsiClp)

# Relocate the OsiXXX targets to the Osi directory, then prepend Coin and Osi.

libTgts := Coin Osi $(patsubst Osi%,Osi/Osi%,$(libTgts))

$(warning Complete dependency list is $(libTgts))

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install clean library unitTest cbc solver solve \
	libdepend libCbc doc

default: install
libCbc: library

libCbc: library

ifneq ($(filter COIN_libOsiCbc,$(CoinLibsDefined)),)
	(cd $(CoinDir)/Osi/OsiCbc && $(MAKE) -f Makefile.lightweight install)
endif

install library: libdepend
	${MAKE} -f Makefile.Cbc $@

# Build the dependencies. OsiCbc is its own strange animal, pick it off
# separately.

libdepend:
	$(foreach tgt,$(libTgts),(cd $(CoinDir)/$(tgt) && $(MAKE) install) ; )
ifneq ($(filter COIN_libOsiCbc,$(CoinLibsDefined)),)
	(cd $(CoinDir)/Osi/OsiCbc && $(MAKE) -f Makefile.lightweight install)
endif

unitTest cbc: 
	(cd Test && ${MAKE} unitTest)

solver: 
	(cd Test && ${MAKE} solver)

solve: 
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
