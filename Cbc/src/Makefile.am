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
libCbc_la_LDFLAGS = -no-undefined

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
CXXFLAGS += \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglProbing` \
	-I`$(CYGPATH_W) $(COINBASESRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src`

if COIN_HAS_CLP
  libCbc_la_SOURCES += Cbc_C_Interface.cpp Cbc_C_Interface.h
  CXXFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src`
endif

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
$(libCbc_la_OBJECTS): ../inc/$(CONFIG_DIST)

nodist_include_HEADERS = ../inc/$(CONFIG_DIST)

../inc/$(CONFIG_DIST): ../inc/$(CONFIG_H)
	echo "/** Required defines from $(CONFIG_H) */" >bla.tmp
	echo "#ifndef __$(CONFIG_DIST)__" | sed -e 's/\./_/g' >>bla.tmp
	echo "#define __$(CONFIG_DIST)__" | sed -e 's/\./_/g' >>bla.tmp
	$(EGREP) 'COIN_|CBC_DEFAULT_SOLVER' $< >>bla.tmp
	echo "#endif" >>bla.tmp
	if test -r ../inc/$(CONFIG_DIST); then \
	  if diff ./inc/$(CONFIG_DIST) bla.tmp >/dev/null 2>1 ; then \
	    rm -f bla.tmp; \
	  else \
	    rm -f ../inc/$(CONFIG_DIST); \
	    mv bla.tmp ../inc/$(CONFIG_DIST); \
	  fi; \
	else \
	  mv bla.tmp ../inc/$(CONFIG_DIST); \
	fi

CLEANFILES = bla.tmp

DISTCLEANFILES = ../inc/$(CONFIG_DIST) \
	CbcOrClpParam.cpp \
	CbcOrClpParam.hpp
