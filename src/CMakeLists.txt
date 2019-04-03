set(NAME "CBC")

set(COIN_${NAME}_CHECKLEVEL "0" CACHE INTERNAL
  "${NAME} check level")
set(COIN_${NAME}_VERBOSITY "0" CACHE INTERNAL
  "${NAME} verbosity level")
configure_file(config.h.cmake.in config.h)
configure_file(config_cbc.h.cmake.in config_cbc.h)
configure_file(config_cbc.h.cmake.in CbcConfig.h)

set(_SRCS
  #CbcBranchBase.cpp
  CbcBranchAllDifferent.cpp
  CbcBranchCut.cpp
  CbcBranchDecision.cpp
  CbcBranchDefaultDecision.cpp
  CbcBranchDynamic.cpp
  CbcBranchingObject.cpp
  CbcBranchLotsize.cpp
  CbcBranchToFixLots.cpp
  CbcCompareDefault.cpp
  CbcCompareDepth.cpp
  CbcCompareEstimate.cpp
  CbcCompareObjective.cpp
  CbcConsequence.cpp
  CbcClique.cpp
  CbcCountRowCut.cpp
  CbcCutGenerator.cpp
  CbcCutModifier.cpp
  CbcCutSubsetModifier.cpp
  CbcDummyBranchingObject.cpp
  CbcEventHandler.cpp
  CbcFathom.cpp
  CbcFathomDynamicProgramming.cpp
  CbcFixVariable.cpp
  CbcFullNodeInfo.cpp
  CbcFollowOn.cpp
  CbcGeneral.cpp
  CbcGeneralDepth.cpp
  CbcHeuristic.cpp
  CbcHeuristicDINS.cpp
  CbcHeuristicDive.cpp
  CbcHeuristicDiveCoefficient.cpp
  CbcHeuristicDiveFractional.cpp
  CbcHeuristicDiveGuided.cpp
  CbcHeuristicDiveLineSearch.cpp
  CbcHeuristicDivePseudoCost.cpp
  CbcHeuristicDiveVectorLength.cpp
  CbcHeuristicFPump.cpp
  CbcHeuristicGreedy.cpp
  CbcHeuristicLocal.cpp
  CbcHeuristicPivotAndFix.cpp
  CbcHeuristicRandRound.cpp
  CbcHeuristicRENS.cpp
  CbcHeuristicRINS.cpp
  CbcHeuristicVND.cpp
  CbcHeuristicDW.cpp
  CbcMessage.cpp
  CbcModel.cpp
  CbcNode.cpp
  CbcNodeInfo.cpp
  CbcNWay.cpp
  CbcObject.cpp
  CbcObjectUpdateData.cpp
  CbcPartialNodeInfo.cpp
  CbcSimpleInteger.cpp
  CbcSimpleIntegerDynamicPseudoCost.cpp
  CbcSimpleIntegerPseudoCost.cpp
  CbcSOS.cpp
  CbcStatistics.cpp
  CbcStrategy.cpp
  CbcSubProblem.cpp
  CbcSymmetry.cpp
  CbcThread.cpp
  CbcTree.cpp
  CbcTreeLocal.cpp
  )

set(_HDRS
  CbcConfig.h
  CbcBranchActual.hpp
  CbcBranchAllDifferent.hpp
  CbcBranchBase.hpp
  CbcBranchCut.hpp
  CbcBranchDecision.hpp
  CbcBranchDefaultDecision.hpp
  CbcBranchDynamic.hpp
  CbcBranchingObject.hpp
  CbcBranchLotsize.hpp
  CbcBranchToFixLots.hpp
  CbcCompareActual.hpp
  CbcCompareBase.hpp
  CbcCompareDefault.hpp
  CbcCompareDepth.hpp
  CbcCompareEstimate.hpp
  CbcCompareObjective.hpp
  CbcConsequence.hpp
  CbcClique.hpp
  CbcCompare.hpp
  CbcCountRowCut.hpp
  CbcCutGenerator.hpp
  CbcCutModifier.hpp
  CbcCutSubsetModifier.hpp
  CbcDummyBranchingObject.hpp
  CbcEventHandler.hpp
  CbcFathom.hpp
  CbcFathomDynamicProgramming.hpp
  CbcFeasibilityBase.hpp
  CbcFixVariable.hpp
  CbcFullNodeInfo.hpp
  CbcFollowOn.hpp
  CbcGeneral.hpp
  CbcGeneralDepth.hpp
  CbcHeuristic.hpp
  CbcHeuristicDINS.hpp
  CbcHeuristicDive.hpp
  CbcHeuristicDiveCoefficient.hpp
  CbcHeuristicDiveFractional.hpp
  CbcHeuristicDiveGuided.hpp
  CbcHeuristicDiveLineSearch.hpp
  CbcHeuristicDivePseudoCost.hpp
  CbcHeuristicDiveVectorLength.hpp
  CbcHeuristicFPump.hpp
  CbcHeuristicGreedy.hpp
  CbcHeuristicLocal.hpp
  CbcHeuristicPivotAndFix.hpp
  CbcHeuristicRandRound.hpp
  CbcHeuristicRENS.hpp
  CbcHeuristicRINS.hpp
  CbcHeuristicVND.hpp
  CbcHeuristicDW.hpp
  CbcMessage.hpp
  CbcModel.hpp
  CbcNode.hpp
  CbcNodeInfo.hpp
  CbcNWay.hpp
  CbcObject.hpp
  CbcObjectUpdateData.hpp
  CbcParam.hpp
  CbcPartialNodeInfo.hpp
  CbcSimpleInteger.hpp
  CbcSimpleIntegerDynamicPseudoCost.hpp
  CbcSimpleIntegerPseudoCost.hpp
  CbcSOS.hpp
  CbcStatistics.hpp
  CbcStrategy.hpp
  CbcSubProblem.hpp
  CbcSymmetry.hpp
  CbcThread.hpp
  CbcTree.hpp
	CbcLinked.hpp
  CbcTreeLocal.hpp
  ClpConstraintAmpl.hpp
  ClpAmplObjective.hpp
  )

add_library(Cbc ${_SRCS} ${_HDRS})
target_include_directories(Cbc PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include/coin>)
target_compile_definitions(Cbc
  PUBLIC HAVE_CONFIG_H
  PRIVATE CBC_BUILD)
if(CMAKE_VERSION VERSION_LESS "3.8.2")
  set_property(TARGET Cbc PROPERTY CXX_STANDARD 11)
  set_property(TARGET Cbc PROPERTY CXX_STANDARD_REQUIRED ON)
else()
  target_compile_features(Cbc PUBLIC cxx_std_11)
endif()
target_link_libraries(Cbc PUBLIC
  Coin::Cgl
  Coin::Clp
  Coin::Osi
  Coin::CoinUtils
  )
set_target_properties(Cbc PROPERTIES
  PUBLIC_HEADER
  "${_HDRS};${CMAKE_CURRENT_BINARY_DIR}/config_cbc.h;${CMAKE_CURRENT_BINARY_DIR}/CbcConfig.h"
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
add_library(Coin::Cbc ALIAS Cbc)

# Install
include(GNUInstallDirs)
install(TARGETS Cbc
  EXPORT ${PROJECT_NAME}Targets
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

#################
##  CbcSolver  ##
#################
configure_file(config_cbcsolver.h.cmake.in config_cbcsolver.h)
configure_file(config_cbcsolver.h.cmake.in CbcSolverConfig.h)

set(_SRCS
  Cbc_C_Interface.cpp
  CbcCbcParam.cpp
  CbcLinked.cpp
  CbcLinkedUtils.cpp
  unitTestClp.cpp
  CbcSolver.cpp
  CbcSolverHeuristics.cpp
  CbcSolverAnalyze.cpp
  CbcMipStartIO.cpp
  CbcSolverExpandKnapsack.cpp
  )

set(_HDRS
  Cbc_C_Interface.h
  CbcLinked.hpp
  CbcSolver.hpp
  CbcSolverHeuristics.hpp
  CbcSolverAnalyze.hpp
  CbcMipStartIO.hpp
  CbcSolverExpandKnapsack.hpp
  )

add_library(CbcSolver ${_SRCS} ${_HDRS})
target_include_directories(CbcSolver PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include/coin>)
target_compile_definitions(CbcSolver
  PUBLIC HAVE_CONFIG_H
  PRIVATE CBC_BUILD CBC_THREAD_SAFE CBC_NO_INTERRUPT)
if(CMAKE_VERSION VERSION_LESS "3.8.2")
  set_property(TARGET CbcSolver PROPERTY CXX_STANDARD 11)
  set_property(TARGET CbcSolver PROPERTY CXX_STANDARD_REQUIRED ON)
else()
  target_compile_features(CbcSolver PUBLIC cxx_std_11)
endif()
target_link_libraries(CbcSolver PUBLIC Coin::Cbc)
set_target_properties(CbcSolver PROPERTIES
  PUBLIC_HEADER "${_HDRS};${CMAKE_CURRENT_BINARY_DIR}/CbcSolverConfig.h"
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
add_library(Coin::CbcSolver ALIAS CbcSolver)

# Install
install(TARGETS CbcSolver
  EXPORT ${PROJECT_NAME}Targets
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

##############
##  OsiCbc  ##
##############
add_subdirectory(OsiCbc)
