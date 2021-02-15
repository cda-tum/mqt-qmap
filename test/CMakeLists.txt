if (NOT TARGET gtest OR NOT TARGET gmock)
	# Prevent overriding the parent project's compiler/linker settings on Windows
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
	add_subdirectory("${PROJECT_SOURCE_DIR}/extern/qfr/extern/dd_package/extern/googletest" "extern/qfr/extern/dd_package/extern/googletest" EXCLUDE_FROM_ALL)
	mark_as_advanced(
			BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
			gmock_build_tests gtest_build_samples gtest_build_tests
			gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
	)
	set_target_properties(gtest gtest_main gmock gmock_main PROPERTIES FOLDER extern)
	if (BINDINGS)
		# adjust visibility settings for building Python bindings
		target_compile_options(gtest PUBLIC -fvisibility=hidden)
		target_compile_options(gmock PUBLIC -fvisibility=hidden)
	endif()
endif()

macro(package_add_test TESTNAME LIBRARY)
	# create an exectuable in which the tests will be stored
	add_executable(${TESTNAME} ${ARGN})
	# link the Google test infrastructure and a default main fuction to the test executable.
	target_link_libraries(${TESTNAME} PRIVATE JKQ::${PROJECT_NAME}_${LIBRARY}_lib gmock gtest_main)
	# discover tests
	gtest_discover_tests(${TESTNAME} WORKING_DIRECTORY ${PROJECT_DIR} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}")
	set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)

	add_custom_command(TARGET ${TESTNAME}
	                   POST_BUILD
	                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${TESTNAME}>/${TESTNAME} ${CMAKE_BINARY_DIR}/${TESTNAME}
	                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../examples $<TARGET_FILE_DIR:${TESTNAME}>/examples
	                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${TESTNAME}>/examples ${CMAKE_BINARY_DIR}/examples
	                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../extern/architectures $<TARGET_FILE_DIR:${TESTNAME}>/architectures
	                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${TESTNAME}>/architectures ${CMAKE_BINARY_DIR}/architectures
	                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../extern/calibration $<TARGET_FILE_DIR:${TESTNAME}>/calibration
	                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${TESTNAME}>/calibration ${CMAKE_BINARY_DIR}/calibration
	                   COMMENT "Copying files and creating symlinks for ${TESTNAME}"
	                   VERBATIM)
endmacro()

package_add_test(${PROJECT_NAME}_heuristic_test heuristic
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_heuristic.cpp)

if (TARGET ${PROJECT_NAME}_exact_lib)
	package_add_test(${PROJECT_NAME}_exact_test exact
	                 ${CMAKE_CURRENT_SOURCE_DIR}/test_exact.cpp)
endif ()


