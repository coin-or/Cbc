set(NAME "CBC")

set(COIN_${NAME}_CHECKLEVEL "0" CACHE INTERNAL
  "${NAME} check level")
set(COIN_${NAME}_VERBOSITY "0" CACHE INTERNAL
  "${NAME} verbosity level")
configure_file(config.h.cmake.in config.h)
configure_file(config_cbc.h.cmake.in config_cbc.h)

set(_SRCS
  #src/CbcBranchBase.cpp
  src/CbcBranchAllDifferent.cpp
  src/CbcBranchCut.cpp
  src/CbcBranchDecision.cpp
  src/CbcBranchDefaultDecision.cpp
  src/CbcBranchDynamic.cpp
  src/CbcBranchLotsize.cpp
  src/CbcBranchToFixLots.cpp
  src/CbcBranchingObject.cpp
  src/CbcClique.cpp
  src/CbcCompareDefault.cpp
  src/CbcCompareDepth.cpp
  src/CbcCompareEstimate.cpp
  src/CbcCompareObjective.cpp
  src/CbcConsequence.cpp
  src/CbcCountRowCut.cpp
  src/CbcCutGenerator.cpp
  src/CbcCutModifier.cpp
  src/CbcCutSubsetModifier.cpp
  src/CbcDummyBranchingObject.cpp
  src/CbcEventHandler.cpp
  src/CbcFathom.cpp
  src/CbcFathomDynamicProgramming.cpp
  src/CbcFixVariable.cpp
  src/CbcFollowOn.cpp
  src/CbcFullNodeInfo.cpp
  src/CbcGeneral.cpp
  src/CbcGeneralDepth.cpp
  src/CbcHeuristic.cpp
  src/CbcHeuristicDINS.cpp
  src/CbcHeuristicDW.cpp
  src/CbcHeuristicDive.cpp
  src/CbcHeuristicDiveCoefficient.cpp
  src/CbcHeuristicDiveFractional.cpp
  src/CbcHeuristicDiveGuided.cpp
  src/CbcHeuristicDiveLineSearch.cpp
  src/CbcHeuristicDivePseudoCost.cpp
  src/CbcHeuristicDiveVectorLength.cpp
  src/CbcHeuristicFPump.cpp
  src/CbcHeuristicGreedy.cpp
  src/CbcHeuristicLocal.cpp
  src/CbcHeuristicPivotAndFix.cpp
  src/CbcHeuristicRENS.cpp
  src/CbcHeuristicRINS.cpp
  src/CbcHeuristicRandRound.cpp
  src/CbcHeuristicVND.cpp
  src/CbcMessage.cpp
  src/CbcModel.cpp
  src/CbcNWay.cpp
  src/CbcNode.cpp
  src/CbcNodeInfo.cpp
  src/CbcObject.cpp
  src/CbcObjectUpdateData.cpp
  src/CbcPartialNodeInfo.cpp
  src/CbcSOS.cpp
  src/CbcSimpleInteger.cpp
  src/CbcSimpleIntegerDynamicPseudoCost.cpp
  src/CbcSimpleIntegerPseudoCost.cpp
  src/CbcStatistics.cpp
  src/CbcStrategy.cpp
  src/CbcSubProblem.cpp
  src/CbcThread.cpp
  src/CbcTree.cpp
  src/CbcTreeLocal.cpp
  )

set(_HDRS
  src/CbcConfig.h
  src/CbcBranchActual.hpp
  src/CbcBranchAllDifferent.hpp
  src/CbcBranchBase.hpp
  src/CbcBranchCut.hpp
  src/CbcBranchDecision.hpp
  src/CbcBranchDefaultDecision.hpp
  src/CbcBranchDynamic.hpp
  src/CbcBranchLotsize.hpp
  src/CbcBranchToFixLots.hpp
  src/CbcBranchingObject.hpp
  src/CbcClique.hpp
  src/CbcCompare.hpp
  src/CbcCompareActual.hpp
  src/CbcCompareBase.hpp
  src/CbcCompareDefault.hpp
  src/CbcCompareDepth.hpp
  src/CbcCompareEstimate.hpp
  src/CbcCompareObjective.hpp
  src/CbcConsequence.hpp
  src/CbcCountRowCut.hpp
  src/CbcCutGenerator.hpp
  src/CbcCutModifier.hpp
  src/CbcCutSubsetModifier.hpp
  src/CbcDummyBranchingObject.hpp
  src/CbcEventHandler.hpp
  src/CbcFathom.hpp
  src/CbcFathomDynamicProgramming.hpp
  src/CbcFeasibilityBase.hpp
  src/CbcFixVariable.hpp
  src/CbcFollowOn.hpp
  src/CbcFullNodeInfo.hpp
  src/CbcGeneral.hpp
  src/CbcGeneralDepth.hpp
  src/CbcHeuristic.hpp
  src/CbcHeuristicDINS.hpp
  src/CbcHeuristicDW.hpp
  src/CbcHeuristicDive.hpp
  src/CbcHeuristicDiveCoefficient.hpp
  src/CbcHeuristicDiveFractional.hpp
  src/CbcHeuristicDiveGuided.hpp
  src/CbcHeuristicDiveLineSearch.hpp
  src/CbcHeuristicDivePseudoCost.hpp
  src/CbcHeuristicDiveVectorLength.hpp
  src/CbcHeuristicFPump.hpp
  src/CbcHeuristicGreedy.hpp
  src/CbcHeuristicLocal.hpp
  src/CbcHeuristicPivotAndFix.hpp
  src/CbcHeuristicRENS.hpp
  src/CbcHeuristicRINS.hpp
  src/CbcHeuristicRandRound.hpp
  src/CbcHeuristicVND.hpp
  src/CbcMessage.hpp
  src/CbcModel.hpp
  src/CbcNWay.hpp
  src/CbcNode.hpp
  src/CbcNodeInfo.hpp
  src/CbcObject.hpp
  src/CbcObjectUpdateData.hpp
  src/CbcParam.hpp
  src/CbcPartialNodeInfo.hpp
  src/CbcSOS.hpp
  src/CbcSimpleInteger.hpp
  src/CbcSimpleIntegerDynamicPseudoCost.hpp
  src/CbcSimpleIntegerPseudoCost.hpp
  src/CbcStrategy.hpp
  src/CbcSubProblem.hpp
  src/CbcTree.hpp
  src/CbcTreeLocal.hpp
  src/ClpAmplObjective.hpp
  src/ClpConstraintAmpl.hpp
  )

add_library(Cbc ${_SRCS} ${_HDRS})
target_include_directories(Cbc PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
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
  PUBLIC_HEADER "${_HDRS};${CMAKE_CURRENT_BINARY_DIR}/config_cbc.h"
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

##############
##  OsiCbc  ##
##############
set(_SRCS src/OsiCbc/OsiCbcSolverInterface.cpp)
set(_HDRS	src/OsiCbc/OsiCbcSolverInterface.hpp)

add_library(OsiCbc ${_SRCS} ${_HDRS})
target_include_directories(OsiCbc PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/OsiCbc>
  $<INSTALL_INTERFACE:include/coin>)
target_compile_definitions(OsiCbc
  PUBLIC HAVE_CONFIG_H
  PRIVATE CBC_BUILD)
if(CMAKE_VERSION VERSION_LESS "3.8.2")
  set_property(TARGET OsiCbc PROPERTY CXX_STANDARD 11)
  set_property(TARGET OsiCbc PROPERTY CXX_STANDARD_REQUIRED ON)
else()
  target_compile_features(OsiCbc PUBLIC cxx_std_11)
endif()
target_link_libraries(OsiCbc PUBLIC
  Coin::CoinUtils
  Coin::Osi
  Coin::Cbc
  )
set_target_properties(OsiCbc PROPERTIES
  PUBLIC_HEADER "${_HDRS}"
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
add_library(Coin::OsiCbc ALIAS OsiCbc)

# Install
install(TARGETS OsiCbc
  EXPORT ${PROJECT_NAME}Targets
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

#################
##  CbcSolver  ##
#################
set(_SRCS
  src/Cbc_C_Interface.cpp
  src/CbcCbcParam.cpp
  src/Cbc_ampl.cpp
  src/CbcLinked.cpp
  src/CbcLinkedUtils.cpp
  src/unitTestClp.cpp
  src/CbcSolver.cpp
  src/CbcSolverHeuristics.cpp
  src/CbcSolverAnalyze.cpp
  src/CbcMipStartIO.cpp
  src/CbcSolverExpandKnapsack.cpp)

set(_HDRS
  src/Cbc_C_Interface.h
  src/Cbc_ampl.h
  src/CbcLinked.hpp
  src/CbcSolver.hpp
  src/CbcMipStartIO.hpp)

add_library(CbcSolver ${_SRCS} ${_HDRS})
target_include_directories(CbcSolver PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
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
  PUBLIC_HEADER "${_HDRS}"
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
