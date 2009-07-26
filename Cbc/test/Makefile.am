# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                         unitTest for Cbc                             #
########################################################################

CBC_TEST_TGTS = gamstests

if COIN_HAS_CLP
  CBC_TEST_TGTS += test_cbc
endif
if CBC_BUILD_CBC_GENERIC
if COIN_HAS_DATA_SAMPLE
  CBC_TEST_TGTS += test_cbc_generic
endif
endif

TESTFLAGS =
if COIN_HAS_DATA_SAMPLE
if COIN_HAS_DATA_SAMPLE_PREINSTALLED
  TESTFLAGS += -dirSample $(DATA_SAMPLEDATADIR)
else
  TESTFLAGS += -dirSample $(DATA_SAMPLEOBJDIR)
endif
endif
TESTFLAGS += -unitTest
if COIN_HAS_DATA_MIPLIB3
if COIN_HAS_DATA_MIPLIB3_PREINSTALLED
  TESTFLAGS += -dirMiplib $(DATA_MIPLIB3DATADIR) -miplib
else
  TESTFLAGS += -dirMiplib $(DATA_MIPLIB3OBJDIR) -miplib
endif
endif

# We are using the solver executables to do the unit test

test: $(CBC_TEST_TGTS)

test_cbc: ../src/cbc$(EXEEXT)
	@echo "calling cbc with flags $(TESTFLAGS)"
	../src/cbc$(EXEEXT) $(TESTFLAGS)

#	if test x"$(DIRMIPLIB)" != x; then \
#		MIPLIBARG="-dirMiplib $(DIRMIPLIB) -miplib"; \
#	fi; \
#	if test x"$(DIRSAMPLE)" != x; then \
#		SAMPLEARG="-dirSample $(DIRSAMPLE)"; \
#	fi; 

if COIN_HAS_DATA_SAMPLE_PREINSTALLED
test_cbc_generic: ../src/cbc-generic$(EXEEXT)
	echo "A minimal test is all that cbc-generic can manage just yet."
	../src/cbc-generic$(EXEEXT) -directory $(DATA_SAMPLEDATADIR) -import p0033.mps -branch
else
test_cbc_generic: ../src/cbc-generic$(EXEEXT)
	echo "A minimal test is all that cbc-generic can manage just yet."
	../src/cbc-generic$(EXEEXT) -directory $(DATA_SAMPLEOBJDIR) -import p0033.mps -branch
endif

gamstests: gamsTest$(EXEEXT)
	./gamsTest$(EXEEXT)

.PHONY: test

bin_PROGRAMS = gamsTest
gamsTest_SOURCES = gamsTest.cpp
gamsTest_LDADD = ../src/libCbcSolver.la ../src/libCbc.la
AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/../src`

gamsTest_LDFLAGS = -rpath $(libdir)
LIBS = `cat $(CBCOBJDIR)/cbc_addlibs.txt` $(ADDLIBS)

if COIN_HAS_CGL_PREINSTALLED
  LIBS += -L`$(CYGPATH_W) $(CGLINSTDIR)/lib` -lCgl
else
  gamsTest_LDADD += $(CGLOBJDIR)/src/libCgl.la
endif

if COIN_HAS_CLP
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiClp
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp`
  gamsTest_LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la
endif
if COIN_HAS_CLP_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(CLPINSTDIR)/lib` -lClp \
	`cat $(CLPDOCDIR)/clp_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
  gamsTest_LDADD += $(CLPOBJDIR)/src/libClp.la
  LIBS += `cat $(CLPOBJDIR)/clp_addlibs.txt`
endif
endif

if COIN_HAS_CPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiCpx
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
  gamsTest_LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
endif
endif

if COIN_HAS_DYLP
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiDylp
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiDylp`
  gamsTest_LDADD += $(OSIOBJDIR)/src/OsiDylp/libOsiDylp.la
endif
if COIN_HAS_DYLP_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(DYLPINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(DYLPINSTDIR)/lib` -lDylp \
	`cat $(DYLPDOCDIR)/dylp_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(DYLPSRCDIR)/src/Dylp` \
	-I`$(CYGPATH_W) $(DYLPSRCDIR)/src/DylpStdLib` \
	-I`$(CYGPATH_W) $(DYLPOBJDIR)/inc`
  gamsTest_LDADD += $(DYLPOBJDIR)/src/Dylp/libDylp.la
  LIBS += `cat $(DYLPOBJDIR)/dylp_addlibs.txt`
endif
endif

if COIN_HAS_FMP
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiFmp
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiFmp`
  gamsTest_LDADD += $(OSIOBJDIR)/src/OsiFmp/libOsiFmp.la
endif
endif

if COIN_HAS_THIRDPARTY_GLPK
if COIN_HAS_THIRDPARTY_GLPK_PREINSTALLED
  LIBS += -L`$(CYGPATH_W) $(THIRDPARTY_GLPKINSTDIR)/lib/ThirdParty` -lcoinglpk
else
  gamsTest_LDADD += $(THIRDPARTY_GLPKOBJDIR)/libcoinglpk.la
endif
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiGlpk
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiGlpk`
  gamsTest_LDADD += $(OSIOBJDIR)/src/OsiGlpk/libOsiGlpk.la
endif
endif

if COIN_HAS_MSK
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(MSKINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiMsk
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiMsk`
  gamsTest_LDADD += $(OSIOBJDIR)/src/OsiMsk/libOsiMsk.la
endif
endif

if COIN_HAS_OSL
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSLINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiOsl
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiOsl`
  gamsTest_LDADD += $(OSIOBJDIR)/src/OsiOsl/libOsiOsl.la
endif
endif

if COIN_HAS_SPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SPXINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiSpx
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSpx`
  gamsTest_LDADD += $(OSIOBJDIR)/src/OsiSpx/libOsiSpx.la
endif
endif

if COIN_HAS_SYMPHONY
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiSym
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSym`
  gamsTest_LDADD += $(OSIOBJDIR)/src/OsiSym/libOsiSym.la
endif
if COIN_HAS_SYMPHONY_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMPHONYINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(SYMPHONYINSTDIR)/lib` -lSym \
	`cat $(SYMPHONYDOCDIR)/sym_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMPHONYSRCDIR)/src`
  gamsTest_LDADD += $(SYMPHONYOBJDIR)/src/libSym.la
  LIBS += `cat $(SYMPHONYOBJDIR)/sym_addlibs.txt`
endif
endif

if COIN_HAS_VOL
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiVol
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiVol`
  gamsTest_LDADD += $(OSIOBJDIR)/src/OsiVol/libOsiVol.la
endif
if COIN_HAS_VOL_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(VOLINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(VOLINSTDIR)/lib` -lVol \
	`cat $(VOLDOCDIR)/vol_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(VOLSRCDIR)/src`
  gamsTest_LDADD += $(VOLOBJDIR)/src/libVol.la
  LIBS += `cat $(VOLOBJDIR)/vol_addlibs.txt`
endif
endif

if COIN_HAS_XPR
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(XPRINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiXpr
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiXpr`
  gamsTest_LDADD += $(OSIOBJDIR)/src/OsiXpr/libOsiXpr.la
endif
endif

if COIN_HAS_ASL
  gamsTest_LDADD += $(ASLLIB)
  AM_CPPFLAGS += $(ASL_CPPFLAGS) -DCOIN_HAS_ASL
endif

if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsi \
	`cat $(OSIDOCDIR)/osi_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSIOBJDIR)/inc`
  gamsTest_LDADD += $(OSIOBJDIR)/src/libOsi.la
  LIBS += `cat $(OSIOBJDIR)/osi_addlibs.txt`
endif

if COIN_HAS_COINUTILS_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(COINUTILSINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(COINUTILSINSTDIR)/lib` -lCoinUtils \
	`cat $(COINUTILSDOCDIR)/coinutils_addlibs.txt`
else
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc`
  gamsTest_LDADD += $(COINUTILSOBJDIR)/src/libCoinUtils.la
  LIBS += `cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt`
endif


########################################################################
#                          Cleaning stuff                              #
########################################################################

# Here we list everything that is not generated by the compiler, e.g.,
# output files of a program

DISTCLEANFILES = yy.mps xx.mps
