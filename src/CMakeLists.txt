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

option(TRY_USING_INSTALLED_QFR_LIBRARY "Attempt searching for installed version of QFR library and use it if found" OFF)
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

if (Z3_FOUND)
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
	target_link_libraries(${PROJECT_NAME}_exact_lib PUBLIC JKQ::qfr z3::z3lib)

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
		target_compile_options(${PROJECT_NAME}_exact_lib PRIVATE -Wall $<$<CONFIG:DEBUG>:-g -Og -Wextra -Wpedantic -pedantic> $<$<CONFIG:RELEASE>:-O3 -mtune=native -march=native -g0>)
	endif()
else()
	message(WARNING "Z3 library not found. Exact mapping library cannot be built without Z3 and will thus not be availbale as target")
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
target_link_libraries(${PROJECT_NAME}_heuristic_lib PUBLIC JKQ::qfr)

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
	target_compile_options(${PROJECT_NAME}_heuristic_lib PRIVATE -Wall $<$<CONFIG:DEBUG>:-g -Og -Wextra -Wpedantic -pedantic> $<$<CONFIG:RELEASE>:-O3 -mtune=native -march=native -g0>)
endif()

# installation instructions

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
if (TARGET ${PROJECT_NAME}_exact_lib)
	install(TARGETS ${PROJECT_NAME}_exact_lib ${PROJECT_NAME}_heuristic_lib
	        EXPORT ${PROJECT_NAME}-targets
	        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
	set_target_properties(${PROJECT_NAME}_exact_lib PROPERTIES EXPORT_NAME ${PROJECT_NAME}_exact_lib)
else ()
	install(TARGETS ${PROJECT_NAME}_heuristic_lib
	        EXPORT ${PROJECT_NAME}-targets
	        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()
set_target_properties(${PROJECT_NAME}_heuristic_lib PROPERTIES EXPORT_NAME ${PROJECT_NAME}_heuristic_lib)

install(DIRECTORY ${${PROJECT_NAME}_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT ${PROJECT_NAME}-targets
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE JKQ::
        DESTINATION ${INSTALL_CONFIGDIR}
        )

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
		${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
		VERSION ${PROJECT_VERSION}
		COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${${PROJECT_NAME}_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
                              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
                              INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
                              )

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
        )

# exporting from the build tree
export(EXPORT ${PROJECT_NAME}-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake NAMESPACE JKQ::)

# Register package in user's package registry
export(PACKAGE ${PROJECT_NAME})
