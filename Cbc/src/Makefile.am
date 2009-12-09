# Copyright (C) 2006, 2007 International Business Machines and others.
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
lib_LTLIBRARIES = libCbc.la libCbcSolver.la

# List all source files for this library, including headers
libCbc_la_SOURCES = \
	Cbc_C_Interface.cpp Cbc_C_Interface.h \
	CbcConfig.h \
	CbcBranchActual.hpp \
	CbcBranchAllDifferent.cpp CbcBranchAllDifferent.hpp \
	CbcBranchBase.hpp \
	CbcBranchCut.cpp CbcBranchCut.hpp \
	CbcBranchDecision.cpp CbcBranchDecision.hpp \
	CbcBranchDefaultDecision.cpp CbcBranchDefaultDecision.hpp \
	CbcBranchDynamic.cpp CbcBranchDynamic.hpp \
	CbcBranchingObject.cpp CbcBranchingObject.hpp \
	CbcBranchLotsize.cpp CbcBranchLotsize.hpp \
	CbcBranchToFixLots.cpp CbcBranchToFixLots.hpp \
	CbcCompareActual.hpp \
	CbcCompareBase.hpp \
	CbcCompareDefault.cpp CbcCompareDefault.hpp \
	CbcCompareDepth.cpp CbcCompareDepth.hpp \
	CbcCompareEstimate.cpp CbcCompareEstimate.hpp \
	CbcCompareObjective.cpp CbcCompareObjective.hpp \
	CbcConsequence.cpp CbcConsequence.hpp \
	CbcChooseVariable.hpp \
	CbcClique.cpp CbcClique.hpp \
	CbcCompare.hpp \
	CbcCountRowCut.cpp CbcCountRowCut.hpp \
	CbcCutGenerator.cpp CbcCutGenerator.hpp \
	CbcCutModifier.cpp CbcCutModifier.hpp \
	CbcCutSubsetModifier.cpp CbcCutSubsetModifier.hpp \
	CbcDummyBranchingObject.cpp CbcDummyBranchingObject.hpp \
	CbcEventHandler.cpp CbcEventHandler.hpp \
	CbcFathom.cpp CbcFathom.hpp \
	CbcFathomDynamicProgramming.cpp CbcFathomDynamicProgramming.hpp \
	CbcFeasibilityBase.hpp \
	CbcFixVariable.cpp CbcFixVariable.hpp \
	CbcFullNodeInfo.cpp CbcFullNodeInfo.hpp \
	CbcFollowOn.cpp CbcFollowOn.hpp \
	CbcGeneral.cpp CbcGeneral.hpp \
	CbcGeneralDepth.cpp CbcGeneralDepth.hpp \
	CbcHeuristic.cpp CbcHeuristic.hpp \
	CbcHeuristicDINS.cpp CbcHeuristicDINS.hpp \
	CbcHeuristicDive.cpp CbcHeuristicDive.hpp \
	CbcHeuristicDiveCoefficient.cpp CbcHeuristicDiveCoefficient.hpp \
	CbcHeuristicDiveFractional.cpp CbcHeuristicDiveFractional.hpp \
	CbcHeuristicDiveGuided.cpp CbcHeuristicDiveGuided.hpp \
	CbcHeuristicDiveLineSearch.cpp CbcHeuristicDiveLineSearch.hpp \
	CbcHeuristicDivePseudoCost.cpp CbcHeuristicDivePseudoCost.hpp \
	CbcHeuristicDiveVectorLength.cpp CbcHeuristicDiveLength.hpp \
	CbcHeuristicFPump.cpp CbcHeuristicFPump.hpp \
	CbcHeuristicGreedy.cpp CbcHeuristicGreedy.hpp \
	CbcHeuristicLocal.cpp CbcHeuristicLocal.hpp \
	CbcHeuristicPivotAndFix.cpp CbcHeuristicPivotAndFix.hpp \
	CbcHeuristicRandRound.cpp CbcHeuristicRandRound.hpp \
	CbcHeuristicRENS.cpp CbcHeuristicRENS.hpp \
	CbcHeuristicRINS.cpp CbcHeuristicRINS.hpp \
	CbcHeuristicVND.cpp CbcHeuristicVND.hpp \
	CbcMessage.cpp CbcMessage.hpp \
	CbcModel.cpp CbcModel.hpp \
	CbcNode.cpp CbcNode.hpp \
	CbcNodeInfo.cpp CbcNodeInfo.hpp \
	CbcNWay.cpp CbcNWay.hpp \
	CbcObject.cpp CbcObject.hpp \
	CbcObjectUpdateData.cpp CbcObjectUpdateData.hpp \
	CbcPartialNodeInfo.cpp CbcPartialNodeInfo.hpp \
	CbcSimpleInteger.cpp CbcSimpleInteger.hpp \
	CbcSimpleIntegerDynamicPseudoCost.cpp \
        CbcSimpleIntegerDynamicPseudoCost.hpp \
	CbcSimpleIntegerPseudoCost.cpp \
        CbcSimpleIntegerPseudoCost.hpp \
	CbcSOS.cpp CbcSOS.hpp \
	CbcStatistics.cpp CbcStatistics.hpp \
	CbcStrategy.cpp CbcStrategy.hpp \
	CbcSubProblem.cpp CbcSubProblem.hpp \
	CbcTree.cpp CbcTree.hpp \
	CbcTreeLocal.cpp CbcTreeLocal.hpp


libCbcSolver_la_SOURCES = \
	CbcCbcParam.cpp \
	Cbc_ampl.cpp Cbc_ampl.h \
	CbcLinked.cpp CbcLinked.hpp CbcLinkedUtils.cpp \
	unitTestClp.cpp CbcSolver.cpp \
	CbcSolverHeuristics.cpp CbcSolverHeuristics.hpp

#Stefan: added unitTestClp.cpp because CbcSolver defines reference to it

## If we have the Ampl solver library, we need to add additional things
## also see lower down
## new idea is to add in source code anyway
##libCbc_la_SOURCES += Cbc_ampl.cpp Cbc_ampl.h


# This is for libtool (on Windows)

libCbc_la_LDFLAGS = $(LT_LDFLAGS)
libCbcSolver_la_LDFLAGS = $(LT_LDFLAGS)

# Some common setup for cbc and cbc-generic.

bin_PROGRAMS =
AM_CPPFLAGS = 

LIBS += `cat $(CGLOBJDIR)/cgl_addlibs.txt` \
	`cat $(OSIOBJDIR)/osi_addlibs.txt` \
	`cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt` $(ADDLIBS)


########################################################################
#                            cbc program                               #
########################################################################

if COIN_HAS_CLP

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory
bin_PROGRAMS += cbc

# List all source files for this executable, including headers
cbc_SOURCES = \
	CoinSolve.cpp

# Additional COIN libraries
cbc_LDADD = libCbcSolver.la \
	libCbc.la \
	$(CGLOBJDIR)/src/libCgl.la \
	$(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
	$(OSIOBJDIR)/src/libOsi.la \
	$(CLPOBJDIR)/src/libClp.la \
	$(COINUTILSOBJDIR)/src/libCoinUtils.la

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
cbc_LDFLAGS = -rpath $(libdir)

# Add locations of include files.

AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
LIBS  += `cat $(CLPOBJDIR)/clp_addlibs.txt`

if COIN_HAS_CPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
endif

# If we have the Ampl solver library, we need to add additional things
if COIN_HAS_ASL
#cbc_SOURCES += Cbc_ampl.cpp Cbc_ampl.h
cbc_LDADD += $(ASLLIB)
AM_CPPFLAGS += $(ASL_CPPFLAGS) -DCOIN_HAS_ASL
endif

endif


########################################################################
#                        cbc-generic program                           #
########################################################################

if CBC_BUILD_CBC_GENERIC

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory

bin_PROGRAMS += cbc-generic

# List all source files for this executable, including headers

cbc_generic_SOURCES = \
	CbcGenBaB.cpp \
	CbcGenCbcParam.cpp CbcGenCbcParam.hpp CbcGenCbcParamUtils.cpp \
	CbcGenCtlBlk.cpp CbcGenCtlBlk.hpp \
	CbcGeneric.cpp \
	CbcGenMessages.cpp CbcGenMessages.hpp \
	CbcGenOsiParam.cpp CbcGenOsiParam.hpp CbcGenOsiParamUtils.cpp \
	CbcGenParam.cpp CbcGenParam.hpp CbcGenParamUtils.cpp \
	CbcGenSolution.cpp \
	CbcGenSolvers.cpp

cbc_generic_LDADD = libCbcSolver.la \
	libCbc.la \
	$(CGLOBJDIR)/src/libCgl.la

# List all additionally required solver and Osi libraries

if COIN_HAS_CLP
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
	$(CLPOBJDIR)/src/libClp.la
  LIBS  += `cat $(CLPOBJDIR)/clp_addlibs.txt`
endif

if COIN_HAS_CPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la 
endif

if COIN_HAS_DYLP
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(DYLPSRCDIR)/src/DylpStdLib` \
        -I`$(CYGPATH_W) $(DYLPSRCDIR)/src/Dylp` \
        -I`$(CYGPATH_W) $(DYLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiDylp`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiDylp/libOsiDylp.la \
	$(DYLPOBJDIR)/src/Dylp/libDylp.la
  LIBS += `cat $(DYLPOBJDIR)/dylp_addlibs.txt`
endif

if COIN_HAS_FMP
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiFmp`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiFmp/libOsiFmp.la
endif

if COIN_HAS_GLPK
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(GLPKINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiGlpk`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiGlpk/libOsiGlpk.la
endif

if COIN_HAS_MSK
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(MSKINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiMsk`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiMsk/libOsiMsk.la
endif

if COIN_HAS_OSL
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSLINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiOsl`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiOsl/libOsiOsl.la
endif

if COIN_HAS_SPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSpx`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiSpx/libOsiSpx.la
endif

if COIN_HAS_SYMPHONY
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMPHONYSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSym`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiSym/libOsiSym.la \
	$(SYMPHONYOBJDIR)/src/libSym.la
  LIBS += `cat $(SYMPHONYOBJDIR)/sym_addlibs.txt`
endif

if COIN_HAS_VOL
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(VOLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiVol`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiVol/libOsiVol.la \
	$(VOLOBJDIR)/src/libVol.la
  LIBS += `cat $(VOLOBJDIR)/vol_addlibs.txt`
endif

if COIN_HAS_XPR
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(XPRINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiXpr`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiXpr/libOsiXpr.la
endif

if COIN_HAS_ASL
  cbc_generic_LDADD += $(ASLLIB)
  AM_CPPFLAGS += $(ASL_CPPFLAGS) -DCOIN_HAS_ASL
endif

# Finally the Osi Coin library (many of the above need it)
cbc_generic_LDADD += $(OSIOBJDIR)/src/libOsi.la \
	$(COINUTILSOBJDIR)/src/libCoinUtils.la

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.

cbc_generic_LDFLAGS = -rpath $(libdir)

endif

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglClique` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglDuplicateRow` \
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
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglLandP` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglResidualCapacity` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglZeroHalf` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSIOBJDIR)/inc`

# List additional defines
AM_CPPFLAGS += -DCOIN_NO_CLP_MESSAGE -DUSE_CBCCONFIG

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	../inc/config_cbc.h \
	Cbc_C_Interface.h \
	CbcBranchActual.hpp \
	CbcBranchAllDifferent.hpp \
	CbcBranchBase.hpp \
	CbcBranchDynamic.hpp \
	CbcBranchCut.hpp \
        CbcBranchDecision.hpp \
	CbcBranchDefaultDecision.hpp \
	CbcBranchingObject.hpp \
	CbcBranchLotsize.hpp \
	CbcBranchToFixLots.hpp \
	CbcCompareActual.hpp \
	CbcCompareDefault.hpp \
	CbcCompareDepth.hpp \
	CbcCompareEstimate.hpp \
	CbcCompareObjective.hpp \
	CbcConsequence.hpp \
	CbcChooseVariable.hpp \
	CbcClique.hpp \
	CbcCompare.hpp \
	CbcCutGenerator.hpp \
	CbcCutModifier.hpp \
	CbcCutSubsetModifier.hpp \
	CbcDummyBranchingObject.hpp \
	CbcFathom.hpp \
	CbcEventHandler.hpp \
	CbcFeasibilityBase.hpp \
	CbcFixVariable.hpp \
	CbcFollowOn.hpp \
	CbcFullNodeInfo.hpp \
	CbcGeneral.hpp \
	CbcGeneralDepth.hpp \
	CbcHeuristic.hpp \
	CbcHeuristicDINS.hpp \
	CbcHeuristicDive.hpp \
	CbcHeuristicDiveCoefficient.hpp \
	CbcHeuristicDiveFractional.hpp \
	CbcHeuristicDiveGuided.hpp \
	CbcHeuristicDiveLineSearch.hpp \
	CbcHeuristicDivePseudoCost.hpp \
	CbcHeuristicDiveVectorLength.hpp \
	CbcHeuristicFPump.hpp \
	CbcHeuristicGreedy.hpp \
	CbcHeuristicLocal.hpp \
	CbcHeuristicPivotAndFix.hpp \
	CbcHeuristicRandRound.hpp \
	CbcHeuristicRENS.hpp \
	CbcHeuristicRINS.hpp \
	CbcHeuristicVND.hpp \
	CbcMessage.hpp \
	CbcModel.hpp \
	CbcNode.hpp \
	CbcNodeInfo.hpp \
	CbcNWay.hpp \
        CbcObject.hpp \
	CbcObjectUpdateData.hpp \
	CbcPartialNodeInfo.hpp \
	CbcSimpleInteger.hpp \
	CbcSimpleIntegerDynamicPseudoCost.hpp \
	CbcSimpleIntegerPseudoCost.hpp \
	CbcStrategy.hpp \
	CbcSolver.hpp \
	CbcSolverHeuristics.hpp \
	CbcSOS.hpp \
	CbcSubProblem.hpp \
	CbcTree.hpp \
	CbcLinked.hpp \
	CbcTreeLocal.hpp \
	ClpConstraintAmpl.hpp \
	ClpAmplObjective.hpp 

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = CbcConfig.h

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includecoindir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/$(ConfigHeader)
