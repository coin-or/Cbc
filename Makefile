# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.Cbc

###############################################################################

export CoinDir := $(shell cd ..; pwd)
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install clean library unitTest solver solve libdepend libCbc doc

default: install
libCbc: library

libCbc: library

install library: libdepend
	${MAKE} -f Makefile.Cbc $@

# Uncomment below to build OsiDylp

libdepend:
	(cd $(CoinDir)/Coin && $(MAKE) install)
	(cd $(CoinDir)/Clp && $(MAKE) install)
ifeq ($(VolDefine),COIN_HAS_VOL)
	(cd $(CoinDir)/Vol && $(MAKE) install)
endif 
ifneq ($(filter COIN_libOsiClp,$(CoinLibsDefined)),)
	(cd $(CoinDir)/Osi/OsiClp && $(MAKE) install)
endif
ifneq ($(filter COIN_libOsiCbc,$(CoinLibsDefined)),)
	(cd $(CoinDir)/Osi/OsiCbc && $(MAKE) -f Makefile.lightweight install)
endif
#	(cd $(CoinDir)/Osi/OsiDylp && $(MAKE) install)
	(cd $(CoinDir)/Cgl && $(MAKE) install)

unitTest: 
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

doc:
	doxygen $(MakefileDir)/doxygen.conf
