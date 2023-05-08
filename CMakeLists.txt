# set required cmake version
cmake_minimum_required(VERSION 3.19)

# This avoids googletest complaining that this (IPO) policy is not set
cmake_policy(SET CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)

project(
  qmap
  LANGUAGES CXX
  DESCRIPTION "MQT QMAP - A library for mapping of quantum circuits to quantum architectures")

# check whether the submodule ``modulename`` is correctly cloned in the ``/extern`` directory.
macro(CHECK_SUBMODULE_PRESENT modulename)
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/${modulename}/CMakeLists.txt")
    message(
      FATAL_ERROR
        "${modulename} submodule not cloned properly. \
        Please run `git submodule update --init --recursive` \
        from the main project directory")
  endif()
endmacro()

check_submodule_present(qfr)

check_submodule_present(LogicBlocks)

# Add path for custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# search for Z3
find_package(Z3 4.8.15)
if(NOT Z3_FOUND)
  message(WARNING "Did not find Z3. Exact library and other depending target will not be available")
endif()

# add main library code
add_subdirectory(src)

# add test code
option(BUILD_QMAP_TESTS "Also build tests for QMAP project")
if(BUILD_QMAP_TESTS)
  enable_testing()
  include(GoogleTest)
  add_subdirectory(test)
endif()

# add Python binding code
option(BINDINGS "Configure for building Python bindings")
if(BINDINGS)
  add_subdirectory(mqt/qmap)
endif()
