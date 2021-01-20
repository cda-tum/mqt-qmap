cmake_minimum_required(VERSION 3.10...3.16)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(qmap
        VERSION 1.2.0
        DESCRIPTION "QMAP  - A JKQ library for mapping of quantum circuits to quantum architectures"
        LANGUAGES CXX)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
SET(COVERAGE OFF CACHE BOOL "Coverage")
option(BUILD_SHARED_LIBS "Build shared libraries (.dll/.so) instead of static ones (.lib/.a)" OFF)
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
	set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
	    STRING "Choose the type of build." FORCE)
	# Set the possible values of build type for cmake-gui
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
	             "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(GIT_SUBMODULE "Check submodules during build" ON)
macro(handle_submodule MODULENAME)
	find_package(Git QUIET)
	if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git" AND GIT_SUBMODULE)
		message(STATUS "${MODULENAME} update")
		execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --remote -- extern/${MODULENAME}
		                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		                RESULT_VARIABLE GIT_SUBMOD_RESULT)
		if(NOT GIT_SUBMOD_RESULT EQUAL "0")
			message(FATAL_ERROR "git submodule update --init --recursive --remote -- extern/${MODULENAME} failed with ${GIT_SUBMOD_RESULT}.")
		endif()
	endif()

	if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/${MODULENAME}/CMakeLists.txt")
		message(FATAL_ERROR "${MODULENAME} failed to download! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
	endif()
endmacro()

handle_submodule(qfr)

# Add path for custom modules
list (APPEND CMAKE_MODULE_PATH
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Z3)
if(Z3_FOUND AND NOT TARGET z3::z3lib)
	add_library(z3::z3lib IMPORTED INTERFACE)
	set_property(TARGET z3::z3lib PROPERTY
	             INTERFACE_INCLUDE_DIRECTORIES ${Z3_CXX_INCLUDE_DIRS})
	set_property(TARGET z3::z3lib PROPERTY
	             INTERFACE_LINK_LIBRARIES ${Z3_LIBRARIES})
endif()

set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package( Boost 1.50 COMPONENTS program_options)

add_subdirectory(src)

option(BUILD_QMAP_TESTS "Also build tests for QMAP project" OFF)
if (BUILD_QMAP_TESTS)
	enable_testing()
	include(GoogleTest)
	add_subdirectory(test)
endif ()

# only build apps if this is the main project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND Boost_FOUND)
	message(STATUS "Boost version: ${Boost_VERSION}")
	add_subdirectory(apps)
endif()

option(BUILD_QMAP_BINDINGS "Also build Python bindings for QMAP project" OFF)
if(BUILD_QMAP_BINDINGS AND Z3_FOUND)
	add_subdirectory(jkq/qmap)
endif()
