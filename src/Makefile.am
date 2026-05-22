# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

########################################################################
#                             libmipster                               #
########################################################################

# Name of the library compiled in this directory.
# We want it to be installed in the 'lib' directory
lib_LTLIBRARIES = libmipster.la

# List all source files for this library, including headers
libmipster_la_SOURCES = \
	CbcConfig.h \
	CbcBab.cpp \
	CbcBranchActual.hpp \
	CbcBranchAllDifferent.cpp CbcBranchAllDifferent.hpp \
	CbcBranchBase.hpp \
	CbcBranchCut.cpp CbcBranchCut.hpp \
	CbcBranchDecision.cpp CbcBranchDecision.hpp \
	CbcBranchDefaultDecision.cpp CbcBranchDefaultDecision.hpp \
	CbcBranchingRanker.cpp CbcBranchingRanker.hpp \
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
	CbcHeuristicDiveConfigurable.cpp CbcHeuristicDiveConfigurable.hpp \
	CbcDiveScheduler.cpp CbcDiveScheduler.hpp \
	CbcRootHeuristicSchedule.cpp CbcRootHeuristicSchedule.hpp \
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
	CbcOutput.cpp CbcOutput.hpp \
	CbcMipStart.cpp CbcMipStart.hpp \
	CbcModel.cpp CbcModel.hpp \
	CbcNode.cpp CbcNode.hpp \
	CbcNodeInfo.cpp CbcNodeInfo.hpp \
	CbcNWay.cpp CbcNWay.hpp \
	CbcObject.cpp CbcObject.hpp \
	CbcObjectUpdateData.cpp CbcObjectUpdateData.hpp \
	CbcParam.cpp CbcParam.hpp \
	CbcParamUtils.cpp CbcParamUtils.hpp \
	CbcPartialNodeInfo.cpp CbcPartialNodeInfo.hpp \
	CbcPartialSolutionPool.cpp CbcPartialSolutionPool.hpp \
	CbcSimpleInteger.cpp CbcSimpleInteger.hpp \
	CbcSimpleIntegerDynamicPseudoCost.cpp \
	CbcSimpleIntegerDynamicPseudoCost.hpp \
	CbcSimpleIntegerPseudoCost.cpp \
	CbcSimpleIntegerPseudoCost.hpp \
	CbcSOS.cpp CbcSOS.hpp CbcSolution.cpp \
	CbcSolver.cpp CbcSolverPreprocess.cpp \
	CbcSolverAnalyze.cpp CbcSolverAnalyze.hpp \
	CbcSolverCutSetup.cpp CbcSolverCutSetup.hpp \
	CbcSolverExpandKnapsack.cpp CbcSolverExpandKnapsack.hpp \
	CbcSolverHeuristics.cpp CbcSolverHeuristics.hpp \
	CbcParameters.cpp CbcParameters.hpp \
	CbcStatistics.cpp CbcStatistics.hpp \
	CbcSolverStatistics.cpp CbcSolverStatistics.hpp \
	CbcInstanceFeatures.cpp CbcInstanceFeatures.hpp \
	CbcLpParamScorer.cpp CbcLpParamScorer.hpp \
	CbcBoundPropagation.cpp CbcBoundPropagation.hpp \
	CbcStrategy.cpp CbcStrategy.hpp \
	CbcSubProblem.cpp CbcSubProblem.hpp \
	CbcThread.cpp CbcThread.hpp \
	CbcTree.cpp CbcTree.hpp \
	CbcTreeLocal.cpp CbcTreeLocal.hpp \
	Cbc_C_Interface.cpp Cbc_C_Interface.h \
	CbcLinked.cpp CbcLinked.hpp CbcLinkedUtils.cpp \
	unitTestClp.cpp

# List all additionally required libraries
libmipster_la_LIBADD = CoinUtils/libCoinUtils.la Clp/libClp.la Cgl/libCgl.la $(CBCLIB_LFLAGS)

# Some common setup for cbc and cbc-generic.
libmipster_la_CPPFLAGS = $(CBCLIB_CFLAGS) -I$(srcdir)/CoinUtils -I$(srcdir)/Clp -I$(srcdir)/Cgl

########################################################################
#                            cbc program                               #
########################################################################

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory
bin_PROGRAMS = mipster

# List all source files for this executable, including headers
mipster_SOURCES = CoinSolve.cpp

# Additional libraries
mipster_LDADD = libmipster.la

# not really accurate to use these CFLAGS, but they should have everything
mipster_CPPFLAGS = $(CBCLIB_CFLAGS) -I$(srcdir)/CoinUtils -I$(srcdir)/Clp -I$(srcdir)/Cgl

########################################################################
#                            Additional flags                          #
########################################################################

# Additional options for libtool
AM_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/mipster'
includecoindir = $(includedir)/mipster
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
	CbcGenMessages.hpp \
	CbcHeuristic.hpp \
	CbcHeuristicDINS.hpp \
	CbcHeuristicDive.hpp \
	CbcHeuristicDiveCoefficient.hpp \
	CbcHeuristicDiveConfigurable.hpp \
	CbcDiveScheduler.hpp \
	CbcRootHeuristicSchedule.hpp \
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
	CbcOutput.hpp \
	CbcModel.hpp \
	CbcNode.hpp \
	CbcNodeInfo.hpp \
	CbcNWay.hpp \
	CbcObject.hpp \
	CbcObjectUpdateData.hpp \
	CbcParam.hpp \
	CbcPartialNodeInfo.hpp \
	CbcPartialSolutionPool.hpp \
	CbcParameters.hpp \
	CbcSimpleInteger.hpp \
	CbcSimpleIntegerDynamicPseudoCost.hpp \
	CbcSimpleIntegerPseudoCost.hpp \
	CbcStrategy.hpp \
	CbcSolver.hpp \
	CbcMipStart.hpp \
	CbcSOS.hpp \
	CbcSubProblem.hpp \
	CbcTree.hpp \
	CbcLinked.hpp \
	CbcTreeLocal.hpp \
	ClpConstraintAmpl.hpp \
	CbcSolverStatistics.hpp \
	CbcInstanceFeatures.hpp \
	CbcLpParamScorer.hpp \
	CbcBoundPropagation.hpp \
	ClpAmplObjective.hpp

#######################################################################
# Create the Config.h file that has all public defines and install it #
#######################################################################

install-exec-local:
	$(install_sh_DATA) config_cbc.h $(DESTDIR)$(includecoindir)/CbcConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/CbcConfig.h $(DESTDIR)$(includecoindir)/CbcSolverConfig.h
