# Makefile.am for CoinUtils — merged into Cbc build tree.
# All utility sources are compiled as a single non-installed library
# that is linked into libCbc.la.

########################################################################
#                           libCoinUtils                               #
########################################################################

noinst_LTLIBRARIES = libCoinUtils.la

libCoinUtils_la_SOURCES = \
	CoinUtilsConfig.h \
	CoinAdjacencyVector.cpp \
	CoinAlloc.cpp \
	CoinBoundPropagation.cpp \
	CoinBronKerbosch.cpp \
	CoinBuild.cpp \
	CoinCliqueExtender.cpp \
	CoinCliqueList.cpp \
	CoinCliqueSet.cpp \
	CoinConflictGraph.cpp \
	CoinCutPool.cpp \
	CoinDenseFactorization.cpp \
	CoinDenseVector.cpp \
	CoinDynamicConflictGraph.cpp \
	CoinError.cpp \
	CoinFactorization1.cpp \
	CoinFactorization2.cpp \
	CoinFactorization3.cpp \
	CoinFactorization4.cpp \
	CoinFileIO.cpp \
	CoinFinite.cpp \
	CoinIndexedVector.cpp \
	CoinKnapsackRow.cpp \
	CoinLpIO.cpp \
	CoinMessage.cpp \
	CoinMessageHandler.cpp \
	CoinModel.cpp \
	CoinModelUseful.cpp \
	CoinModelUseful2.cpp \
	CoinMpsIO.cpp \
	CoinNodeHeap.cpp \
	CoinOddWheelSeparator.cpp \
	CoinOslFactorization.cpp \
	CoinOslFactorization2.cpp \
	CoinOslFactorization3.cpp \
	CoinPackedMatrix.cpp \
	CoinPackedVector.cpp \
	CoinPackedVectorBase.cpp \
	CoinParam.cpp \
	CoinParamUtils.cpp \
	CoinPostsolveMatrix.cpp \
	CoinPrePostsolveMatrix.cpp \
	CoinPresolveDoubleton.cpp \
	CoinPresolveDual.cpp \
	CoinPresolveDupcol.cpp \
	CoinPresolveEmpty.cpp \
	CoinPresolveFixed.cpp \
	CoinPresolveForcing.cpp \
	CoinPresolveHelperFunctions.cpp \
	CoinPresolveImpliedFree.cpp \
	CoinPresolveIsolated.cpp \
	CoinPresolveMatrix.cpp \
	CoinPresolveMonitor.cpp \
	CoinPresolvePsdebug.cpp \
	CoinPresolveSingleton.cpp \
	CoinPresolveSubst.cpp \
	CoinPresolveTighten.cpp \
	CoinPresolveTripleton.cpp \
	CoinPresolveUseless.cpp \
	CoinPresolveZeros.cpp \
	CoinRational.cpp \
	CoinSearchTree.cpp \
	CoinShallowPackedVector.cpp \
	CoinShortestPath.cpp \
	CoinSimpFactorization.cpp \
	CoinSnapshot.cpp \
	CoinStaticConflictGraph.cpp \
	CoinStructuredModel.cpp \
	CoinTable.cpp \
	CoinTerm.cpp \
	CoinWarmStartBasis.cpp \
	CoinWarmStartDual.cpp \
	CoinWarmStartPrimalDual.cpp \
	CoinWarmStartVector.cpp \
	CoinAdjacencyVector.hpp \
	CoinAlloc.hpp \
	CoinBoundPropagation.hpp \
	CoinBronKerbosch.hpp \
	CoinBuild.hpp \
	CoinCliqueExtender.hpp \
	CoinCliqueList.hpp \
	CoinCliqueSet.hpp \
	CoinColumnType.hpp \
	CoinConflictGraph.hpp \
	CoinCutPool.hpp \
	CoinDenseFactorization.hpp \
	CoinDenseVector.hpp \
	CoinDistance.hpp \
	CoinDynamicConflictGraph.hpp \
	CoinError.hpp \
	CoinFactorization.hpp \
	CoinFileIO.hpp \
	CoinFinite.hpp \
	CoinFloatEqual.hpp \
	CoinHelperFunctions.hpp \
	CoinIndexedVector.hpp \
	CoinKnapsackRow.hpp \
	CoinLpIO.hpp \
	CoinMessage.hpp \
	CoinMessageHandler.hpp \
	CoinModel.hpp \
	CoinModelUseful.hpp \
	CoinMpsIO.hpp \
	CoinNodeHeap.hpp \
	CoinOddWheelSeparator.hpp \
	CoinOslFactorization.hpp \
	CoinPackedMatrix.hpp \
	CoinPackedVector.hpp \
	CoinPackedVectorBase.hpp \
	CoinParam.hpp \
	CoinPragma.hpp \
	CoinPresolveDoubleton.hpp \
	CoinPresolveDual.hpp \
	CoinPresolveDupcol.hpp \
	CoinPresolveEmpty.hpp \
	CoinPresolveFixed.hpp \
	CoinPresolveForcing.hpp \
	CoinPresolveImpliedFree.hpp \
	CoinPresolveIsolated.hpp \
	CoinPresolveMatrix.hpp \
	CoinPresolveMonitor.hpp \
	CoinPresolvePsdebug.hpp \
	CoinPresolveSingleton.hpp \
	CoinPresolveSubst.hpp \
	CoinPresolveTighten.hpp \
	CoinPresolveTripleton.hpp \
	CoinPresolveUseless.hpp \
	CoinPresolveZeros.hpp \
	CoinRational.hpp \
	CoinRowType.hpp \
	CoinSearchTree.hpp \
	CoinShallowPackedVector.hpp \
	CoinShortestPath.hpp \
	CoinSignal.hpp \
	CoinSimpFactorization.hpp \
	CoinSmartPtr.hpp \
	CoinSnapshot.hpp \
	CoinSort.hpp \
	CoinStaticConflictGraph.hpp \
	CoinStructuredModel.hpp \
	CoinTable.hpp \
	CoinTerm.hpp \
	CoinTime.hpp \
	CoinUtility.hpp \
	CoinWarmStart.hpp \
	CoinWarmStartBasis.hpp \
	CoinWarmStartDual.hpp \
	CoinWarmStartPrimalDual.hpp \
	CoinWarmStartVector.hpp \
	CoinOslC.h \
	CoinTypes.h \
	config_coinutils_default.h \
	config_default.h \
	configall_system.h \
	configall_system_msc.h

AM_CPPFLAGS = $(CBCLIB_CFLAGS) -I$(srcdir)

AM_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#              Headers installed into $(includedir)/coin-or            #
########################################################################

includecoindir = $(includedir)/mipster
includecoin_HEADERS = \
	CoinUtilsConfig.h \
	CoinAlloc.hpp \
	CoinBuild.hpp \
	CoinDenseVector.hpp \
	CoinDistance.hpp \
	CoinError.hpp \
	CoinFactorization.hpp \
	CoinDenseFactorization.hpp \
	CoinFileIO.hpp \
	CoinFinite.hpp \
	CoinFloatEqual.hpp \
	CoinHelperFunctions.hpp \
	CoinIndexedVector.hpp \
	CoinLpIO.hpp \
	CoinMessage.hpp \
	CoinMessageHandler.hpp \
	CoinModel.hpp \
	CoinStructuredModel.hpp \
	CoinModelUseful.hpp \
	CoinMpsIO.hpp \
	CoinPackedMatrix.hpp \
	CoinPackedVector.hpp \
	CoinPackedVectorBase.hpp \
	CoinParam.hpp \
	CoinPragma.hpp \
	CoinPresolveDoubleton.hpp \
	CoinPresolveDual.hpp \
	CoinPresolveDupcol.hpp \
	CoinPresolveEmpty.hpp \
	CoinPresolveFixed.hpp \
	CoinPresolveForcing.hpp \
	CoinPresolveImpliedFree.hpp \
	CoinPresolveIsolated.hpp \
	CoinPresolveMatrix.hpp \
	CoinPresolveMonitor.hpp \
	CoinPresolvePsdebug.hpp \
	CoinPresolveSingleton.hpp \
	CoinPresolveSubst.hpp \
	CoinPresolveTighten.hpp \
	CoinPresolveTripleton.hpp \
	CoinPresolveUseless.hpp \
	CoinPresolveZeros.hpp \
	CoinRational.hpp \
	CoinSearchTree.hpp \
	CoinShallowPackedVector.hpp \
	CoinSignal.hpp \
	CoinSmartPtr.hpp \
	CoinSnapshot.hpp \
	CoinSort.hpp \
	CoinTime.hpp \
	CoinTypes.h \
	CoinUtility.hpp \
	CoinWarmStart.hpp \
	CoinWarmStartBasis.hpp \
	CoinWarmStartVector.hpp \
	CoinWarmStartDual.hpp \
	CoinWarmStartPrimalDual.hpp \
	CoinAdjacencyVector.hpp \
	CoinBronKerbosch.hpp \
	CoinCliqueExtender.hpp \
	CoinCliqueList.hpp \
	CoinCliqueSet.hpp \
	CoinConflictGraph.hpp \
	CoinCutPool.hpp \
	CoinDynamicConflictGraph.hpp \
	CoinKnapsackRow.hpp \
	CoinBoundPropagation.hpp \
	CoinTerm.hpp \
	CoinNodeHeap.hpp \
	CoinOddWheelSeparator.hpp \
	CoinShortestPath.hpp \
	CoinStaticConflictGraph.hpp \
	CoinColumnType.hpp \
	CoinTable.hpp \
	CoinOslFactorization.hpp \
	CoinSimpFactorization.hpp

noinst_HEADERS = \
	CoinRowType.hpp \
	CoinOslC.h \
	config_coinutils_default.h \
	config_default.h \
	configall_system.h \
	configall_system_msc.h
