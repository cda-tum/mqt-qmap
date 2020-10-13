# fix cmake warnings for googletest when building shared libs
if (APPLE AND BUILD_SHARED_LIBS)
	set(CMAKE_MACOSX_RPATH ON)
endif ()

if (NOT TARGET gtest OR NOT TARGET gmock)
	# Prevent overriding the parent project's compiler/linker settings on Windows
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
	add_subdirectory("${PROJECT_SOURCE_DIR}/extern/qfr/extern/dd_package/extern/googletest" "extern/qfr/extern/dd_package/extern/googletest" EXCLUDE_FROM_ALL)
	mark_as_advanced(
			BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
			gmock_build_tests gtest_build_samples gtest_build_tests
			gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
	)
	set_target_properties(gtest PROPERTIES FOLDER extern)
	set_target_properties(gtest_main PROPERTIES FOLDER extern)
	set_target_properties(gmock PROPERTIES FOLDER extern)
	set_target_properties(gmock_main PROPERTIES FOLDER extern)
	target_compile_features(gtest PUBLIC cxx_std_14)
	target_compile_features(gtest_main PUBLIC cxx_std_14)
	target_compile_features(gmock PUBLIC cxx_std_14)
	target_compile_features(gmock_main PUBLIC cxx_std_14)
endif()

macro(package_add_test TESTNAME LIBRARY)
	# create an exectuable in which the tests will be stored
	add_executable(${TESTNAME} ${ARGN})
	# link the Google test infrastructure, mocking library, and a default main fuction to
	# the test executable.  Remove g_test_main if writing your own main function.
	target_link_libraries(${TESTNAME} PRIVATE JKQ::${PROJECT_NAME}_${LIBRARY}_lib gmock gtest_main)
	set_property(TARGET ${TESTNAME} PROPERTY CXX_EXTENSIONS OFF)
	if (COVERAGE)
		target_compile_options(${TESTNAME} PUBLIC --coverage)
		target_link_libraries(${TESTNAME} PUBLIC --coverage)
	endif()
	if(MSVC)
		target_compile_options(${TESTNAME} PUBLIC "/W4" "$<$<CONFIG:RELEASE>:/O2 /GL>")
		target_link_options(${TESTNAME} PUBLIC "$<$<CONFIG:Release>:/LTCG>")
	else()
		target_compile_options(${TESTNAME} PUBLIC -Wall $<$<CONFIG:DEBUG>:-g3 -Og -Wextra -Wpedantic -pedantic> $<$<CONFIG:RELEASE>:-O3 -mtune=native -march=native -g>)
	endif()

	# gtest_discover_tests replaces gtest_add_tests,
	# see https://cmake.org/cmake/help/v3.10/module/GoogleTest.html for more options to pass to it
	gtest_discover_tests(${TESTNAME}
	                     # set a working directory so your project root so that you can find test data via paths relative to the project root
	                     WORKING_DIRECTORY ${PROJECT_DIR}
	                     PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
	                     )
	set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

package_add_test(${PROJECT_NAME}_heuristic_test heuristic
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_heuristic.cpp)

add_custom_command(TARGET ${PROJECT_NAME}_heuristic_test
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${PROJECT_NAME}_heuristic_test>/${PROJECT_NAME}_heuristic_test ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_heuristic_test
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../examples $<TARGET_FILE_DIR:${PROJECT_NAME}_heuristic_test>/examples
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${PROJECT_NAME}_heuristic_test>/examples ${CMAKE_BINARY_DIR}/examples
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../extern/architectures $<TARGET_FILE_DIR:${PROJECT_NAME}_heuristic_test>/architectures
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${PROJECT_NAME}_heuristic_test>/architectures ${CMAKE_BINARY_DIR}/architectures
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../extern/calibration $<TARGET_FILE_DIR:${PROJECT_NAME}_heuristic_test>/calibration
                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${PROJECT_NAME}_heuristic_test>/calibration ${CMAKE_BINARY_DIR}/calibration
                   COMMENT "Copying files and creating symlinks for ${PROJECT_NAME}_heuristic_test"
                   VERBATIM)

if (TARGET ${PROJECT_NAME}_exact_lib)
	package_add_test(${PROJECT_NAME}_exact_test exact
	                 ${CMAKE_CURRENT_SOURCE_DIR}/test_exact.cpp)

	add_custom_command(TARGET ${PROJECT_NAME}_exact_test
	                   POST_BUILD
	                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${PROJECT_NAME}_exact_test>/${PROJECT_NAME}_exact_test ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_exact_test
	                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../examples $<TARGET_FILE_DIR:${PROJECT_NAME}_exact_test>/examples
	                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${PROJECT_NAME}_exact_test>/examples ${CMAKE_BINARY_DIR}/examples
	                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../extern/architectures $<TARGET_FILE_DIR:${PROJECT_NAME}_exact_test>/architectures
	                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${PROJECT_NAME}_exact_test>/architectures ${CMAKE_BINARY_DIR}/architectures
	                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../extern/calibration $<TARGET_FILE_DIR:${PROJECT_NAME}_exact_test>/calibration
	                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${PROJECT_NAME}_exact_test>/calibration ${CMAKE_BINARY_DIR}/calibration
	                   COMMENT "Copying files and creating symlinks for ${PROJECT_NAME}_exact_test"
	                   VERBATIM)
endif ()


