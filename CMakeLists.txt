cmake_minimum_required(VERSION 3.10)
project(QMAP
        VERSION 0.9.2
        DESCRIPTION "QMAP  - A library for mapping of quantum circuits to quantum architectures"
        LANGUAGES CXX)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

SET(COVERAGE OFF CACHE BOOL "Coverage")

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
	set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
	    STRING "Choose the type of build." FORCE)
	# Set the possible values of build type for cmake-gui
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
	             "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


FIND_PACKAGE( Boost 1.50 COMPONENTS program_options REQUIRED )

add_subdirectory(src)

#if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
#	enable_testing()
#	add_subdirectory(test)
#endif()