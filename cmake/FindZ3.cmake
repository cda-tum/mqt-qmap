# Code adapted from https://fossies.org/linux/llvm/cmake/modules/FindZ3.cmake
include(CheckCXXSourceRuns)

# if Z3_ROOT is provided, check there first
set(Z3_ROOT
    ""
    CACHE PATH "Root of Z3 distribution.")
if(DEFINED ENV{Z3_ROOT})
  set(Z3_ROOT $ENV{Z3_ROOT})
  message("Z3_ROOT: ${Z3_ROOT}")
endif()
if(NOT ${Z3_ROOT} STREQUAL "")
  find_path(
    Z3_CXX_INCLUDE_DIRS
    NAMES z3.h z3++.h
    NO_DEFAULT_PATH
    PATHS ${Z3_ROOT}/include
    PATH_SUFFIXES libz3 z3)

  find_library(
    Z3_LIBRARIES
    NAMES z3 libz3
    NO_DEFAULT_PATH
    PATHS ${Z3_ROOT}
    PATH_SUFFIXES lib bin)
endif()

# see if a config file is available
if(NOT Z3_CXX_INCLUDE_DIRS OR NOT Z3_LIBRARIES)
  find_package(Z3 CONFIG)
endif()

# try default paths as a last hope
if(NOT Z3_FOUND)
  find_path(
    Z3_CXX_INCLUDE_DIRS
    NAMES z3.h z3++.h
    PATH_SUFFIXES libz3 z3)
  find_library(
    Z3_LIBRARIES
    NAMES z3 libz3
    PATH_SUFFIXES lib bin)

  unset(Z3_VERSION_STRING)

  # Try to check version of Z3
  if(Z3_CXX_INCLUDE_DIRS AND (EXISTS "${Z3_CXX_INCLUDE_DIRS}/z3_version.h"))
    file(STRINGS "${Z3_CXX_INCLUDE_DIRS}/z3_version.h" z3_version_str
         REGEX "^#define[\t ]+Z3_MAJOR_VERSION[\t ]+.*")
    string(REGEX REPLACE "^.*Z3_MAJOR_VERSION[\t ]+([0-9]*).*$" "\\1" Z3_MAJOR "${z3_version_str}")

    file(STRINGS "${Z3_CXX_INCLUDE_DIRS}/z3_version.h" z3_version_str
         REGEX "^#define[\t ]+Z3_MINOR_VERSION[\t ]+.*")
    string(REGEX REPLACE "^.*Z3_MINOR_VERSION[\t ]+([0-9]*).*$" "\\1" Z3_MINOR "${z3_version_str}")

    file(STRINGS "${Z3_CXX_INCLUDE_DIRS}/z3_version.h" z3_version_str
         REGEX "^#define[\t ]+Z3_BUILD_NUMBER[\t ]+.*")
    string(REGEX REPLACE "^.*Z3_BUILD_NUMBER[\t ]+([0-9]*).*$" "\\1" Z3_BUILD "${z3_version_str}")

    set(Z3_VERSION_STRING ${Z3_MAJOR}.${Z3_MINOR}.${Z3_BUILD})
    unset(z3_version_str)
  endif()

  if(NOT Z3_VERSION_STRING)
    set(Z3_VERSION_STRING "0.0.0")
  endif()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(
    Z3
    REQUIRED_VARS Z3_LIBRARIES Z3_CXX_INCLUDE_DIRS
    VERSION_VAR Z3_VERSION_STRING)
  mark_as_advanced(Z3_CXX_INCLUDE_DIRS Z3_LIBRARIES)
endif()
