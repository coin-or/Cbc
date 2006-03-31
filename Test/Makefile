# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.test
# - Makefile.solve

###############################################################################

export CoinDir := $(shell cd ../..; pwd)
export MakefileDir := $(CoinDir)/Makefiles

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default clean cbc solver

default: solver

# In spite of what you might think, OsiCbc will be built. See the explanation
# in ../Makefile. The exclusion here is just a hedge against change.

cbc:
ifneq ($(filter COIN_libOsiCbc,$(CoinLibsDefined)),)
	(cd $(CoinDir)/Osi/OsiCbc && $(MAKE) -f Makefile.lightweight install)
endif
	${MAKE} -f Makefile.test $@

solver:
	${MAKE} -f Makefile.solve $@

clean:
	@rm -rf Junk
	@rm -rf $(UNAME)
	@rm -rf dep
