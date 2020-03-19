# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

########################################################################
#                               libCbc                                 #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libCbc.la libCbcSolver.la

# List all source files for this library, including headers
libCbc_la_SOURCES = \
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
	CbcHeuristicDW.cpp CbcHeuristicDW.hpp \
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
	CbcSymmetry.cpp CbcSymmetry.hpp \
	CbcThread.cpp CbcThread.hpp \
	CbcTree.cpp CbcTree.hpp \
	CbcTreeLocal.cpp CbcTreeLocal.hpp

libCbcSolver_la_SOURCES = \
	Cbc_C_Interface.cpp Cbc_C_Interface.h \
	CbcCbcParam.cpp \
	CbcLinked.cpp CbcLinked.hpp CbcLinkedUtils.cpp \
	unitTestClp.cpp CbcSolver.cpp \
	CbcSolverHeuristics.cpp CbcSolverHeuristics.hpp \
	CbcSolverAnalyze.cpp CbcSolverAnalyze.hpp \
	CbcMipStartIO.cpp CbcMipStartIO.hpp \
	CbcSolverExpandKnapsack.cpp CbcSolverExpandKnapsack.hpp

# List all additionally required libraries

libCbc_la_LIBADD = $(CBCLIB_LFLAGS)
libCbcSolver_la_LIBADD = $(CBCSOLVERLIB_LFLAGS) libCbc.la

# Some common setup for cbc and cbc-generic.

libCbc_la_CPPFLAGS = -DCBC_BUILD $(CBCLIB_CFLAGS)
libCbcSolver_la_CPPFLAGS = -DCBCSOLVER_BUILD $(CBCSOLVERLIB_CFLAGS)

AM_CPPFLAGS =
if COIN_HAS_NAUTY
  AM_CPPFLAGS += -DCOIN_HAS_NTY $(NAUTY_CFLAGS)
endif

########################################################################
#                            cbc program                               #
########################################################################

bin_PROGRAMS =

if COIN_HAS_CLP

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory
bin_PROGRAMS += cbc

# List all source files for this executable, including headers
cbc_SOURCES = CoinSolve.cpp

# Additional libraries
cbc_LDADD = libCbcSolver.la libCbc.la

# not really accurate to use these CFLAGS, but they should have everything
cbc_CPPFLAGS = $(CBCSOLVERLIB_CFLAGS)

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

cbc_generic_LDADD = libCbcSolver.la libCbc.la $(CBCGENERIC_LFLAGS)

# List all additionally required solver and Osi libraries
# the linker flags for all available LP solvers should already be included in $CBCGENERIC_LFLAGS,
# so just add compiler flags here (all we can think of)

AM_CPPFLAGS += $(CBCGENERIC_CFLAGS)

endif

########################################################################
#                            Additional flags                          #
########################################################################

# List additional defines
AM_CPPFLAGS += -DCOIN_NO_CLP_MESSAGE -DUSE_CBCCONFIG

# Additional options for libtool
AM_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(pkgincludedir)
includecoin_HEADERS = \
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
	CbcCompareBase.hpp \
	CbcCompareDefault.hpp \
	CbcCompareDepth.hpp \
	CbcCompareEstimate.hpp \
	CbcCompareObjective.hpp \
	CbcConsequence.hpp \
	CbcClique.hpp \
	CbcCompare.hpp \
	CbcCountRowCut.hpp \
	CbcCutGenerator.hpp \
	CbcCutModifier.hpp \
	CbcCutSubsetModifier.hpp \
	CbcDummyBranchingObject.hpp \
	CbcFathom.hpp \
	CbcEventHandler.hpp \
	CbcFathomDynamicProgramming.hpp \
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
	CbcHeuristicDW.hpp \
	CbcMessage.hpp \
	CbcModel.hpp \
	CbcNode.hpp \
	CbcNodeInfo.hpp \
	CbcNWay.hpp \
	CbcObject.hpp \
	CbcObjectUpdateData.hpp \
	CbcParam.hpp \
	CbcPartialNodeInfo.hpp \
	CbcSimpleInteger.hpp \
	CbcSimpleIntegerDynamicPseudoCost.hpp \
	CbcSimpleIntegerPseudoCost.hpp \
	CbcStrategy.hpp \
	CbcSolver.hpp \
	CbcMipStartIO.hpp \
	CbcSOS.hpp \
	CbcSubProblem.hpp \
	CbcTree.hpp \
	CbcLinked.hpp \
	CbcTreeLocal.hpp \
	ClpConstraintAmpl.hpp \
	ClpAmplObjective.hpp 

#######################################################################
# Create the Config.h file that has all public defines and install it #
#######################################################################

install-exec-local:
	$(install_sh_DATA) config_cbc.h $(DESTDIR)$(includecoindir)/CbcConfig.h
	$(install_sh_DATA) config_cbcsolver.h $(DESTDIR)$(includecoindir)/CbcSolverConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/CbcConfig.h $(DESTDIR)$(includecoindir)/CbcSolverConfig.h
