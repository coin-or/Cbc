# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libCbc                                 #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libCbc.la

# List all source files for this library, including headers
libCbc_la_SOURCES = \
	CbcBranchActual.cpp CbcBranchActual.hpp \
	CbcBranchBase.cpp CbcBranchBase.hpp \
	CbcBranchCut.cpp CbcBranchCut.hpp \
	CbcBranchDynamic.cpp CbcBranchDynamic.hpp \
	CbcBranchLotsize.cpp CbcBranchLotsize.hpp \
	CbcCompareActual.cpp CbcCompareActual.hpp \
	CbcCompareBase.hpp \
	CbcCountRowCut.cpp CbcCountRowCut.hpp \
	CbcCutGenerator.cpp CbcCutGenerator.hpp \
	CbcEventHandler.cpp CbcEventHandler.hpp \
	CbcFathom.cpp CbcFathom.hpp \
	CbcFathomDynamicProgramming.cpp CbcFathomDynamicProgramming.hpp \
	CbcFeasibilityBase.hpp \
	CbcHeuristic.cpp CbcHeuristic.hpp \
	CbcHeuristicFPump.cpp CbcHeuristicFPump.hpp \
	CbcHeuristicGreedy.cpp CbcHeuristicGreedy.hpp \
	CbcHeuristicLocal.cpp CbcHeuristicLocal.hpp \
	CbcMessage.cpp CbcMessage.hpp \
	CbcModel.cpp CbcModel.hpp \
	CbcNode.cpp CbcNode.hpp \
	CbcStatistics.cpp CbcStatistics.hpp \
	CbcStrategy.cpp CbcStrategy.hpp \
	CbcTree.cpp CbcTree.hpp \
	CbcTreeLocal.cpp CbcTreeLocal.hpp

# This is for libtool (on Windows)
libCbc_la_LDFLAGS = #-no-undefined

########################################################################
#                            cbc program                               #
########################################################################

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory
bin_PROGRAMS = cbc

# List all source files for this executable, including headers
cbc_SOURCES = \
	CbcMain.cpp \
	CbcParam.cpp CbcParam.hpp \
	unitTest.cpp

# The following file is not to be distributed, since it is "borrowed" from Clp
nodist_cbc_SOURCES = CbcOrClpParam.cpp

cbc_LDADD = libCbc.la \
	$(CGLOBJDIR)/src/libCgl.la \
	$(OSIOBJDIR)/src/OsiCbc/libOsiCbc.la
CXXFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCbc`

# List all additionally required solver and Osi libraries
if COIN_HAS_CLP
  CXXFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
	$(CLPOBJDIR)/src/libClp.la
  LIBS  += `cat $(CLPOBJDIR)/clp_addlibs.txt`
endif

if COIN_HAS_CPX
  CXXFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
endif

if COIN_HAS_DYLP
  CXXFLAGS += -I`$(CYGPATH_W) $(DYLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiDylp`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiDylp/libOsiDylp.la \
	$(DYLPOBJDIR)/src/libDylp.la
  LIBS += `cat $(DYLPOBJDIR)/dylp_addlibs.txt`
endif

if COIN_HAS_FMP
  CXXFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiFmp`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiFmp/libOsiFmp.la
endif

if COIN_HAS_GLPK
  CXXFLAGS += -I`$(CYGPATH_W) $(GLPKINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiGlpk`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiGlpk/libOsiGlpk.la
endif

if COIN_HAS_MSK
  CXXFLAGS += -I`$(CYGPATH_W) $(MSKINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiMsk`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiMsk/libOsiMsk.la
endif

if COIN_HAS_OSL
  CXXFLAGS += -I`$(CYGPATH_W) $(OSLINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiOsl`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiOsl/libOsiOsl.la
endif

if COIN_HAS_SPX
  CXXFLAGS += -I`$(CYGPATH_W) $(SPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSpx`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiSpx/libOsiSpx.la
endif

if COIN_HAS_SYM
  CXXFLAGS += -I`$(CYGPATH_W) $(SYMSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSym`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiSym/libOsiSym.la \
	$(SYMOBJDIR)/src/libSym.la
  LIBS += `cat $(SYMOBJDIR)/sym_addlibs.txt`
endif

if COIN_HAS_VOL
  CXXFLAGS += -I`$(CYGPATH_W) $(VOLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiVol`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiVol/libOsiVol.la \
	$(VOLOBJDIR)/src/libVol.la
  LIBS += `cat $(VOLOBJDIR)/vol_addlibs.txt`
endif

if COIN_HAS_XPR
  CXXFLAGS += -I`$(CYGPATH_W) $(XPRINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiXpr`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiXpr/libOsiXpr.la
endif

# Finally the Osi Coin library (many of the above need it)
cbc_LDADD += $(OSIOBJDIR)/src/libOsi.la \
	$(COINUTILSOBJDIR)/src/libCoinUtils.la
LIBS += `cat $(OSIOBJDIR)/osi_addlibs.txt` \
	`cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt` $(ADDLIBS)

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
cbc_LDFLAGS = -rpath $(libdir)

########################################################################
#                           solve program                              #
########################################################################

#ToDo: I'm not sure if I'm doing the right thing here

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory
bin_PROGRAMS += solve

# List all source files for this executable, including headers
solve_SOURCES = \
	CoinSolve.cpp \
	Cbc_ampl.cpp Cbc_ampl.h \
	unitTestClp.cpp

# The following file is not to be distributed, since it is "borrowed" from Clp
nodist_solve_SOURCES = CbcOrClpParam.cpp

# Additional COIN libraries
solve_LDADD = libCbc.la \
	$(CGLOBJDIR)/src/libCgl.la \
	$(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
	$(OSIOBJDIR)/src/libOsi.la \
	$(CLPOBJDIR)/src/libClp.la \
	$(COINUTILSOBJDIR)/src/libCoinUtils.la

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
solve_LDFLAGS = -rpath $(libdir)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
CXXFLAGS += \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglClique` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglFlowCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglGomory` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglKnapsackCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding2` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglOddHole` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglPreProcess` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglProbing` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglRedSplit` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglTwomir` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src`

# List additional defines
DEFS += -DCOIN_NO_CLP_MESSAGE -DUSE_CBCCONFIG

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install'
include_HEADERS = \
	CbcBranchActual.hpp \
	CbcBranchBase.hpp \
	CbcCompareActual.hpp \
	CbcCompareBase.hpp \
	CbcCutGenerator.hpp \
	CbcEventHandler.hpp \
	CbcHeuristic.hpp \
	CbcHeuristicFPump.hpp \
	CbcHeuristicGreedy.hpp \
	CbcHeuristicLocal.hpp \
	CbcMessage.hpp \
	CbcModel.hpp \
	CbcNode.hpp \
	CbcStrategy.hpp \
	CbcTree.hpp \
	CbcTreeLocal.hpp

########################################################################
#     Generating the configure header that is to be distributed        #
########################################################################

CONFIG_DIST = CbcConfig.hpp
CONFIG_H = config_cbc.h
$(libCbc_la_OBJECTS) $(cbc_OBJECTS) $(solve_OBJECTS): ../inc/$(CONFIG_DIST)

nodist_include_HEADERS = ../inc/$(CONFIG_DIST)

../inc/$(CONFIG_DIST): ../inc/$(CONFIG_H)
	echo "/** Required defines from $(CONFIG_H) */" >bla.tmp
	echo "#ifndef __$(CONFIG_DIST)__" | sed -e 's/\./_/g' >>bla.tmp
	echo "#define __$(CONFIG_DIST)__" | sed -e 's/\./_/g' >>bla.tmp
	$(EGREP) 'COIN_|CBC_DEFAULT_SOLVER' $< >>bla.tmp
	echo "#endif" >>bla.tmp
	if test -r ../inc/$(CONFIG_DIST); then \
	  if diff ../inc/$(CONFIG_DIST) bla.tmp >/dev/null 2>&1 ; then \
	    rm -f bla.tmp; \
	  else \
	    rm -f ../inc/$(CONFIG_DIST); \
	    mv bla.tmp ../inc/$(CONFIG_DIST); \
	  fi; \
	else \
	  mv bla.tmp ../inc/$(CONFIG_DIST); \
	fi

CLEANFILES = bla.tmp

DISTCLEANFILES = ../inc/$(CONFIG_DIST)
