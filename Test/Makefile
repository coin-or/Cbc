# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.test
# - Makefile.solve

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

solver: libdepend
	${MAKE} -f Makefile.solve $*

clean:
	@rm -rf Junk
	@rm -rf $(UNAME)
	@rm -rf dep
