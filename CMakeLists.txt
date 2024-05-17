# set required cmake version
cmake_minimum_required(VERSION 3.19...3.28)

project(
  mqt-qmap
  LANGUAGES CXX
  DESCRIPTION "MQT QMAP - A library for mapping of quantum circuits to quantum architectures")

option(BUILD_MQT_QMAP_TESTS "Also build tests for the MQT QMAP project" ON)
option(BUILD_MQT_QMAP_BINDINGS "Build the MQT QMAP Python bindings" OFF)

if(BUILD_MQT_QMAP_BINDINGS)
  # ensure that the BINDINGS option is set
  set(BINDINGS
      ON
      CACHE INTERNAL "Enable settings related to Python bindings")
  # Some common settings for finding Python
  set(Python_FIND_VIRTUALENV
      FIRST
      CACHE STRING "Give precedence to virtualenvs when searching for Python")
  set(Python_FIND_FRAMEWORK
      LAST
      CACHE STRING "Prefer Brew/Conda to Apple framework Python")
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

include(cmake/ExternalDependencies.cmake)

# set the include directory for the build tree
set(MQT_QMAP_INCLUDE_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

# add main library code
add_subdirectory(src)

# add test code
if(BUILD_MQT_QMAP_TESTS)
  enable_testing()
  include(GoogleTest)
  add_subdirectory(test)
endif()

if(NOT TARGET mqt-qmap-uninstall)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in
                 ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)
  add_custom_target(mqt-qmap-uninstall COMMAND ${CMAKE_COMMAND} -P
                                               ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
