# allow loader to find lib
if(APPLE)
	set(basePoint @loader_path)
else()
	set(basePoint $ORIGIN)
endif()
include(GNUInstallDirs)
file(RELATIVE_PATH relDir
     ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
     ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
     )
set(CMAKE_INSTALL_RPATH ${basePoint} ${basePoint}/${relDir})

option(TRY_USING_INSTALLED_QFR_LIBRARY "Attempt searching for installed version of QFR library and use it if found" ON)
if (TRY_USING_INSTALLED_QFR_LIBRARY)
	find_package(qfr 1.0 QUIET)
	if (qfr_FOUND)
		message(STATUS "using installed QFR library")
	else()
		message(STATUS "using QFR submodule")
		add_subdirectory("${PROJECT_SOURCE_DIR}/extern/qfr" "extern/qfr")
	endif()
else()
	message(STATUS "using QFR submodule")
	add_subdirectory("${PROJECT_SOURCE_DIR}/extern/qfr" "extern/qfr")
endif ()

find_package(Z3 CONFIG)
if (Z3_FOUND)
	if(NOT TARGET z3::z3lib)
		add_library(z3::z3lib IMPORTED INTERFACE)
		set_property(TARGET z3::z3lib PROPERTY
		             INTERFACE_INCLUDE_DIRECTORIES ${Z3_CXX_INCLUDE_DIRS})
		set_property(TARGET z3::z3lib PROPERTY
		             INTERFACE_LINK_LIBRARIES ${Z3_LIBRARIES})
	endif()
	message(STATUS "Z3 version: ${Z3_VERSION_STRING}")

	add_library(${PROJECT_NAME}_exact_lib
	            ${CMAKE_CURRENT_SOURCE_DIR}/utils.cpp
	            ${CMAKE_CURRENT_SOURCE_DIR}/Architecture.cpp
	            ${CMAKE_CURRENT_SOURCE_DIR}/Mapper.cpp
	            ${CMAKE_CURRENT_SOURCE_DIR}/exact/ExactMapper.cpp
	            ${${PROJECT_NAME}_SOURCE_DIR}/include/utils.hpp
	            ${${PROJECT_NAME}_SOURCE_DIR}/include/Architecture.hpp
	            ${${PROJECT_NAME}_SOURCE_DIR}/include/Mapper.hpp
	            ${${PROJECT_NAME}_SOURCE_DIR}/include/MappingResults.hpp
	            ${${PROJECT_NAME}_SOURCE_DIR}/include/MappingSettings.hpp
	            ${${PROJECT_NAME}_SOURCE_DIR}/include/exact/ExactMapper.hpp)
	add_library(JKQ::${PROJECT_NAME}_exact_lib ALIAS ${PROJECT_NAME}_exact_lib)
	target_include_directories(${PROJECT_NAME}_exact_lib PUBLIC
	                           $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include/>)
	target_compile_features(${PROJECT_NAME}_exact_lib PUBLIC cxx_std_14)
	set_target_properties(${PROJECT_NAME}_exact_lib PROPERTIES CXX_EXTENSIONS OFF)
	target_link_libraries(${PROJECT_NAME}_exact_lib PUBLIC JKQ::qfr nlohmann_json::nlohmann_json PRIVATE z3::z3lib)

	# add coverage compiler and linker flag if COVERAGE is set
	if (COVERAGE)
		target_compile_options(${PROJECT_NAME}_exact_lib PRIVATE --coverage)
		target_link_libraries(${PROJECT_NAME}_exact_lib PRIVATE --coverage)
	endif()

	# check if interprocedural optimization is supported
	include(CheckIPOSupported)
	check_ipo_supported(RESULT ipo_supported)
	if(ipo_supported)
		set_target_properties(${PROJECT_NAME}_exact_lib PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
	endif()

	# set compiler flags depending on compiler
	if(MSVC)
		target_compile_options(${PROJECT_NAME}_exact_lib PRIVATE "/W4" "/GR" "$<$<CONFIG:RELEASE>:/O2 /GL>")
		target_link_options(${PROJECT_NAME}_exact_lib PRIVATE "$<$<CONFIG:Release>:/LTCG>")
	else()
		target_compile_options(${PROJECT_NAME}_exact_lib PRIVATE -Wall $<$<CONFIG:DEBUG>:-g3 -Og -Wextra -Wpedantic -pedantic> $<$<CONFIG:RELEASE>:-O3 -mtune=native -march=native -DNDEBUG -g0>)
	endif()
else()
	message(WARNING "Z3 library not found. Exact mapping library cannot be built without Z3")
endif ()

add_library(${PROJECT_NAME}_heuristic_lib
            ${CMAKE_CURRENT_SOURCE_DIR}/utils.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Architecture.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Mapper.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/heuristic/HeuristicMapper.cpp
            ${${PROJECT_NAME}_SOURCE_DIR}/include/utils.hpp
            ${${PROJECT_NAME}_SOURCE_DIR}/include/Architecture.hpp
            ${${PROJECT_NAME}_SOURCE_DIR}/include/Mapper.hpp
            ${${PROJECT_NAME}_SOURCE_DIR}/include/MappingResults.hpp
            ${${PROJECT_NAME}_SOURCE_DIR}/include/MappingSettings.hpp
            ${${PROJECT_NAME}_SOURCE_DIR}/include/heuristic/HeuristicMapper.hpp
            ${${PROJECT_NAME}_SOURCE_DIR}/include/heuristic/unique_priority_queue.hpp)
add_library(JKQ::${PROJECT_NAME}_heuristic_lib ALIAS ${PROJECT_NAME}_heuristic_lib)
target_include_directories(${PROJECT_NAME}_heuristic_lib PUBLIC
                           $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include/>)
target_compile_features(${PROJECT_NAME}_heuristic_lib PUBLIC cxx_std_14)
set_target_properties(${PROJECT_NAME}_heuristic_lib PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(${PROJECT_NAME}_heuristic_lib PUBLIC JKQ::qfr nlohmann_json::nlohmann_json)

# add coverage compiler and linker flag if COVERAGE is set
if (COVERAGE)
	target_compile_options(${PROJECT_NAME}_heuristic_lib PRIVATE --coverage)
	target_link_libraries(${PROJECT_NAME}_heuristic_lib PRIVATE --coverage)
endif()

# check if interprocedural optimization is supported
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
	set_target_properties(${PROJECT_NAME}_heuristic_lib PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# set compiler flags depending on compiler
if(MSVC)
	target_compile_options(${PROJECT_NAME}_heuristic_lib PRIVATE "/W4" "/GR" "$<$<CONFIG:RELEASE>:/O2 /GL>")
	target_link_options(${PROJECT_NAME}_heuristic_lib PRIVATE "$<$<CONFIG:Release>:/LTCG>")
else()
	target_compile_options(${PROJECT_NAME}_heuristic_lib PRIVATE -Wall $<$<CONFIG:DEBUG>:-g3 -Og -Wextra -Wpedantic -pedantic> $<$<CONFIG:RELEASE>:-O3 -mtune=native -march=native -DNDEBUG -g0>)
endif()

