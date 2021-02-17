# add JKQ::qfr target
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/qfr" "extern/qfr")

macro(add_qmap_library LIBNAME SRCFILE)
	# main project library
	add_library(${PROJECT_NAME}_${LIBNAME}_lib
	            ${CMAKE_CURRENT_SOURCE_DIR}/utils.cpp
	            ${CMAKE_CURRENT_SOURCE_DIR}/Architecture.cpp
	            ${CMAKE_CURRENT_SOURCE_DIR}/Mapper.cpp
	            ${CMAKE_CURRENT_SOURCE_DIR}/${LIBNAME}/${SRCFILE}.cpp
	            ${${PROJECT_NAME}_SOURCE_DIR}/include/utils.hpp
	            ${${PROJECT_NAME}_SOURCE_DIR}/include/Architecture.hpp
	            ${${PROJECT_NAME}_SOURCE_DIR}/include/Mapper.hpp
	            ${${PROJECT_NAME}_SOURCE_DIR}/include/MappingResults.hpp
	            ${${PROJECT_NAME}_SOURCE_DIR}/include/MappingSettings.hpp
	            ${${PROJECT_NAME}_SOURCE_DIR}/include/${LIBNAME}/${SRCFILE}.hpp)

	# set include directories
	target_include_directories(${PROJECT_NAME}_${LIBNAME}_lib PUBLIC $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include/>)

	# set required C++ standard and disable compiler specific extensions
	target_compile_features(${PROJECT_NAME}_${LIBNAME}_lib PUBLIC cxx_std_17)
	set_target_properties(${PROJECT_NAME}_${LIBNAME}_lib PROPERTIES CMAKE_CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

	# add JKQ::qfr library
	target_link_libraries(${PROJECT_NAME}_${LIBNAME}_lib PUBLIC JKQ::qfr z3::z3lib)

	# enable interprocedural optimization if it is supported
	include(CheckIPOSupported)
	check_ipo_supported(RESULT ipo_supported)
	if(ipo_supported)
		set_target_properties(${PROJECT_NAME}_${LIBNAME}_lib PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
	endif()

	# set compiler flags depending on compiler
	if(MSVC)
		target_compile_options(${PROJECT_NAME}_${LIBNAME}_lib PUBLIC /utf-8)
	else()
		target_compile_options(${PROJECT_NAME}_${LIBNAME}_lib PUBLIC -Wall -Wextra $<$<CONFIG:DEBUG>:-Og>)
		if (BINDINGS)
			# adjust visibility settings for building Python bindings
			target_compile_options(${PROJECT_NAME}_${LIBNAME}_lib PUBLIC -fvisibility=hidden)
		endif()
		if(NOT DEPLOY)
			# only include machine-specific optimizations when building for the host machine
			target_compile_options(${PROJECT_NAME}_${LIBNAME}_lib PUBLIC -mtune=native -march=native)
		endif ()
	endif()

	if(GENERATE_POSITION_INDEPENDENT_CODE OR BINDINGS)
		include(CheckPIESupported)
		check_pie_supported()
		set_target_properties(${PROJECT_NAME}_${LIBNAME}_lib PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
	endif()

	# add coverage compiler and linker flag to the library and all targets that link against it, if COVERAGE is set
	if (COVERAGE)
		target_compile_options(${PROJECT_NAME}_${LIBNAME}_lib PUBLIC --coverage)
		target_link_libraries(${PROJECT_NAME}_${LIBNAME}_lib PUBLIC --coverage)
	endif()

	# add JKQ alias
	add_library(JKQ::${PROJECT_NAME}_${LIBNAME}_lib ALIAS ${PROJECT_NAME}_${LIBNAME}_lib)
endmacro()

if (Z3_FOUND)
	# exact mapper project library
	add_qmap_library(exact ExactMapper)
	# add Z3 library
	target_link_libraries(${PROJECT_NAME}_exact_lib PUBLIC z3::z3lib)
else()
	message(WARNING "Z3 library not found. Exact mapping library cannot be built without Z3 and will thus not be availbale as target")
endif ()

# heuristic mapper project library
add_qmap_library(heuristic HeuristicMapper)
