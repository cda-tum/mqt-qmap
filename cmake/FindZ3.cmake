# Code adapted from https://fossies.org/linux/llvm/cmake/modules/FindZ3.cmake
INCLUDE(CheckCXXSourceRuns)

# Function to check Z3's version
function(check_z3_version z3_include z3_lib)
	# The program that will be executed to print Z3's version.
	file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testz3.cpp
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
			Z3_RETURNCODE
			Z3_COMPILED
			${CMAKE_BINARY_DIR}
			${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testz3.cpp
			COMPILE_DEFINITIONS -I"${z3_include}"
			LINK_LIBRARIES -L${z3_lib_path} -lz3
			RUN_OUTPUT_VARIABLE SRC_OUTPUT
	)

	if(Z3_COMPILED)
		string(REGEX REPLACE "([0-9]*\\.[0-9]*\\.[0-9]*)" "\\1"
		       z3_version "${SRC_OUTPUT}")
		set(Z3_VERSION_STRING ${z3_version} PARENT_SCOPE)
	endif()
endfunction(check_z3_version)

# if Z3_ROOT is provided, check there first
set(Z3_ROOT "" CACHE PATH "Root of Z3 distribution.")
if (DEFINED ENV{Z3_ROOT})
	set(Z3_ROOT $ENV{Z3_ROOT})
	message("Z3_ROOT: ${Z3_ROOT}")
endif ()
if (NOT ${Z3_ROOT} STREQUAL "")
	find_path(Z3_CXX_INCLUDE_DIRS NAMES z3.h z3++.h
	          NO_DEFAULT_PATH
	          PATHS ${Z3_ROOT}/include
	          PATH_SUFFIXES libz3 z3)

	find_library(Z3_LIBRARIES NAMES z3 libz3
	             NO_DEFAULT_PATH
	             PATHS ${Z3_ROOT}
	             PATH_SUFFIXES lib bin)
endif ()

# see if a config file is available
if (NOT Z3_CXX_INCLUDE_DIRS OR NOT Z3_LIBRARIES)
	find_package(Z3 CONFIG)
endif()

# try default paths as a last hope
if (NOT Z3_FOUND)
	find_path(Z3_CXX_INCLUDE_DIRS NAMES z3.h z3++.h
	          PATH_SUFFIXES libz3 z3)
	find_library(Z3_LIBRARIES NAMES z3 libz3
	             PATH_SUFFIXES lib bin)

	unset(Z3_VERSION_STRING)

	# Try to check version by compiling a small program that prints Z3's version
	if(Z3_CXX_INCLUDE_DIRS AND Z3_LIBRARIES)
		check_z3_version(${Z3_CXX_INCLUDE_DIRS} ${Z3_LIBRARIES})
	endif()

	if(NOT Z3_VERSION_STRING)
		set(Z3_VERSION_STRING "0.0.0")
	endif()

	include (FindPackageHandleStandardArgs)
	find_package_handle_standard_args(Z3
	                                  REQUIRED_VARS Z3_LIBRARIES Z3_CXX_INCLUDE_DIRS
	                                  VERSION_VAR Z3_VERSION_STRING)
	mark_as_advanced(Z3_CXX_INCLUDE_DIRS Z3_LIBRARIES)
endif ()
