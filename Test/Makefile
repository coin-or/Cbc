# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.test

###############################################################################

export CoinDir := $(shell cd ../..; pwd)
export MakefileDir := $(CoinDir)/Makefiles

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default clean unitTest libdepend

default: unitTest

libdepend :
	(cd .. && ${MAKE} install)

unitTest: libdepend
	${MAKE} -f Makefile.test $*

clean:
	@rm -rf Junk
	@rm -rf $(UNAME)
	@rm -rf dep
