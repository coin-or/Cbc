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


cbc:
	(cd $(CoinDir)/Osi/OsiCbc && $(MAKE) -f Makefile.lightweight install)
	${MAKE} -f Makefile.test $@

solver:
	${MAKE} -f Makefile.solve $@

clean:
	@rm -rf Junk
	@rm -rf $(UNAME)
	@rm -rf dep
