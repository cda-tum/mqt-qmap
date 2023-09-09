# set required cmake version
cmake_minimum_required(VERSION 3.19...3.27)

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

check_submodule_present(mqt-core)
check_submodule_present(LogicBlocks)

option(BUILD_MQT_QMAP_BINDINGS "Build the MQT QMAP Python bindings" OFF)
if(BUILD_MQT_QMAP_BINDINGS)
  # ensure that the BINDINGS option is set
  set(BINDINGS
      ON
      CACHE BOOL "Enable settings related to Python bindings" FORCE)
  # cmake-lint: disable=C0103
  set(Python_FIND_VIRTUALENV
      FIRST
      CACHE STRING "Give precedence to virtualenvs when searching for Python")
  # cmake-lint: disable=C0103
  set(Python_ARTIFACTS_INTERACTIVE
      ON
      CACHE BOOL "Prevent multiple searches for Python and instead cache the results.")
  # top-level call to find Python
  find_package(
    Python 3.8 REQUIRED
    COMPONENTS Interpreter Development.Module
    OPTIONAL_COMPONENTS Development.SABIModule)
endif()

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
option(BUILD_MQT_QMAP_TESTS "Also build tests for the MQT QMAP project" ON)
if(BUILD_MQT_QMAP_TESTS)
  enable_testing()
  include(GoogleTest)
  add_subdirectory(test)
endif()
