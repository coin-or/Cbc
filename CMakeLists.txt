cmake_minimum_required(VERSION 3.5)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(ParseAc)
parse_ac(VERSION MAJOR MINOR PATCH)

project(Cbc VERSION ${VERSION})

# config options
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif(NOT CMAKE_BUILD_TYPE)

option(BUILD_SHARED_LIBS "" ON)
if(BUILD_SHARED_LIBS AND MSVC)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif(BUILD_SHARED_LIBS AND MSVC)

# config options
if(MSVC)
  # Build with multiple processes
  add_definitions(/MP)
  add_definitions(/D_CRT_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_DEPRECATE)
  # MSVC warning suppressions
  add_definitions(
    /wd4018 # 'expression' : signed/unsigned mismatch
    /wd4065 # switch statement contains 'default' but no 'case' labels
    /wd4101 # 'identifier' : unreferenced local variable
    /wd4102 # 'label' : unreferenced label
    /wd4244 # 'conversion' conversion from 'type1' to 'type2', possible loss of data
    /wd4267 # 'var' : conversion from 'size_t' to 'type', possible loss of data
    /wd4309 # 'conversion' : truncation of constant value
    /wd4805 # 'operation' : unsafe mix of type 'type1' and type 'type2' in operation.
    /wd4996 # The compiler encountered a deprecated declaration.
    )
endif()
if(APPLE)
  set(
    CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wno-inconsistent-missing-override -Wno-unused-command-line-argument -Wno-unused-result -Wno-exceptions"
    )
  set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9" CACHE STRING "Minimum OS X deployment version")
endif(APPLE)

# ZLIB
find_package(ZLIB)
if(ZLIB_FOUND)
  set(HAVE_ZLIB_H "1" CACHE INTERNAL "Use zlib")
  set(COIN_HAS_ZLIB "1" CACHE INTERNAL "Use zlib")
endif()

# PThread
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)
if(CMAKE_USE_PTHREADS_INIT)
  set(PTHREADS_FOUND TRUE)
else()
  set(PTHREADS_FOUND FALSE)
endif()

find_package(CoinUtils REQUIRED CONFIG)
find_package(Osi REQUIRED CONFIG)
find_package(Clp REQUIRED CONFIG)
find_package(Cgl REQUIRED CONFIG)

include(CheckEnv)
include(CTest)

add_subdirectory(Cbc)

include(GNUInstallDirs)
install(EXPORT ${PROJECT_NAME}Targets
  NAMESPACE Coin::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/Config.cmake.in
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  COMPATIBILITY SameMajorVersion)
install(
  FILES
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  COMPONENT Devel)
