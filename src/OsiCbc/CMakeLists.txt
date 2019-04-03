set(NAME "OSICBC")

configure_file(config_osicbc.h.cmake.in config_osicbc.h)
configure_file(config_osicbc.h.cmake.in OsiCbcConfig.h)

set(_SRCS OsiCbcSolverInterface.cpp)
set(_HDRS	OsiCbcSolverInterface.hpp)

add_library(OsiCbc ${_SRCS} ${_HDRS})
target_include_directories(OsiCbc PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
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
  PUBLIC_HEADER
  "${_HDRS};${CMAKE_CURRENT_BINARY_DIR}/config_osicbc.h;${CMAKE_CURRENT_BINARY_DIR}/OsiCbcConfig.h"
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
