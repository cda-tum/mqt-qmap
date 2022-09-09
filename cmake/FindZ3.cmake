# Code adapted from https://fossies.org/linux/llvm/cmake/modules/FindZ3.cmake
include(CheckCXXSourceRuns)

# Function to check Z3's version
function(check_z3_version z3_include z3_lib)
  # The program that will be executed to print Z3's version.
  file(
    WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testz3.cpp
    "#include <assert.h>
        #include <z3.h>
        int main() {
          unsigned int major, minor, build, rev;
          Z3_get_version(&major, &minor, &build, &rev);
          printf(\"%u.%u.%u\", major, minor, build);
          return 0;
       }")

  # Get lib path
  get_filename_component(z3_lib_path ${z3_lib} PATH)

  try_run(
    Z3_RETURNCODE Z3_COMPILED ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testz3.cpp
    COMPILE_DEFINITIONS -I"${z3_include}" LINK_LIBRARIES -L${z3_lib_path} -lz3
    RUN_OUTPUT_VARIABLE SRC_OUTPUT)

  if(Z3_COMPILED)
    string(REGEX REPLACE "([0-9]*\\.[0-9]*\\.[0-9]*)" "\\1" z3_version "${SRC_OUTPUT}")
    set(Z3_VERSION_STRING
        ${z3_version}
        PARENT_SCOPE)
  endif()
endfunction(check_z3_version)

set(Z3_ROOT
    ""
    CACHE PATH "Root of Z3 distribution.")
if(DEFINED ENV{Z3_ROOT})
  set(Z3_ROOT $ENV{Z3_ROOT})
  message(STATUS "Z3_ROOT from environment: ${Z3_ROOT}")
endif()

# if Z3_ROOT is provided, check there first
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

  if(Z3_CXX_INCLUDE_DIRS AND Z3_LIBRARIES)
    message(STATUS "Z3_ROOT provided and includes and libraries found.")
    message(STATUS "Z3_CXX_INCLUDE_DIRS: ${Z3_CXX_INCLUDE_DIRS}")
    message(STATUS "Z3_LIBRARIES: ${Z3_LIBRARIES}")
  endif()
endif()

# see if a config file is available
if(NOT Z3_CXX_INCLUDE_DIRS OR NOT Z3_LIBRARIES)
  find_package(Z3 CONFIG)
  if(Z3_FOUND)
    message(STATUS "Found Z3 using config file")
  endif()
endif()

if(NOT Z3_FOUND)
  # if the include directories or libraries have not been found, look in the system paths
  if(NOT Z3_CXX_INCLUDE_DIRS OR NOT Z3_LIBRARIES)
    find_path(
      Z3_CXX_INCLUDE_DIRS
      NAMES z3.h z3++.h
      PATH_SUFFIXES libz3 z3)
    find_library(
      Z3_LIBRARIES
      NAMES z3 libz3
      PATH_SUFFIXES lib bin)

    if(Z3_CXX_INCLUDE_DIRS AND Z3_LIBRARIES)
      message(STATUS "Found Z3 includes and libraries in system paths.")
      message(STATUS "Z3_CXX_INCLUDE_DIRS: ${Z3_CXX_INCLUDE_DIRS}")
      message(STATUS "Z3_LIBRARIES: ${Z3_LIBRARIES}")
    endif()
  endif()

  # if they are still not found, try to find them with Python as a last resort
  if(NOT Z3_CXX_INCLUDE_DIRS OR NOT Z3_LIBRARIES)
    find_package(Python COMPONENTS Interpreter)
    if(Python_FOUND)
      execute_process(
        COMMAND ${Python_EXECUTABLE} -c "import os, z3; print(os.path.dirname(z3.__file__))"
        RESULT_VARIABLE Z3_PYTHON_ROOT)
      if(Z3_PYTHON_ROOT)
        find_path(
          Z3_CXX_INCLUDE_DIRS
          NAMES z3.h z3++.h
          NO_DEFAULT_PATH
          PATHS ${Z3_PYTHON_ROOT}/include
          PATH_SUFFIXES libz3 z3)

        find_library(
          Z3_LIBRARIES
          NAMES z3 libz3
          NO_DEFAULT_PATH
          PATHS ${Z3_PYTHON_ROOT}
          PATH_SUFFIXES lib bin)

        if(Z3_CXX_INCLUDE_DIRS AND Z3_LIBRARIES)
          message(STATUS "Found Z3 includes and libraries from Python installation.")
          message(STATUS "Z3_CXX_INCLUDE_DIRS: ${Z3_CXX_INCLUDE_DIRS}")
          message(STATUS "Z3_LIBRARIES: ${Z3_LIBRARIES}")
        endif()
      endif()
    endif()
  endif()

  # try to determine version of Z3
  unset(Z3_VERSION_STRING)

  # Try to check version by compiling a small program that prints Z3's version
  if(Z3_CXX_INCLUDE_DIRS AND Z3_LIBRARIES)
    check_z3_version(${Z3_CXX_INCLUDE_DIRS} ${Z3_LIBRARIES})
  endif()

  if(NOT Z3_VERSION_STRING AND (Z3_CXX_INCLUDE_DIRS AND EXISTS "${Z3_CXX_INCLUDE_DIRS}/z3_version.h"
                               ))
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
    message(STATUS "Could not determine Z3 version")
    set(Z3_VERSION_STRING "0.0.0")
  endif()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(
    Z3
    REQUIRED_VARS Z3_LIBRARIES Z3_CXX_INCLUDE_DIRS
    VERSION_VAR Z3_VERSION_STRING)
  mark_as_advanced(Z3_CXX_INCLUDE_DIRS Z3_LIBRARIES)
endif()
