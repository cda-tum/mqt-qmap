if(NOT TARGET gtest OR NOT TARGET gmock)
  # Prevent overriding the parent project's compiler/linker settings on Windows
  set(gtest_force_shared_crt # cmake-lint: disable=C0103
      ON
      CACHE BOOL "" FORCE)
  add_subdirectory("${PROJECT_SOURCE_DIR}/extern/qfr/extern/dd_package/extern/googletest"
                   "extern/qfr/extern/dd_package/extern/googletest" EXCLUDE_FROM_ALL)
  mark_as_advanced(
    BUILD_GMOCK
    BUILD_GTEST
    BUILD_SHARED_LIBS
    gmock_build_tests
    gtest_build_samples
    gtest_build_tests
    gtest_disable_pthreads
    gtest_force_shared_crt
    gtest_hide_internal_symbols)
  set_target_properties(gtest gtest_main gmock gmock_main PROPERTIES FOLDER extern)
  if(BINDINGS)
    # adjust visibility settings for building Python bindings
    target_compile_options(gtest PUBLIC -fvisibility=hidden)
    target_compile_options(gmock PUBLIC -fvisibility=hidden)
  endif()
endif()

# macro to add a test executable for one of the project libraries
macro(PACKAGE_ADD_TEST testname library)
  # create an executable in which the tests will be stored
  add_executable(${testname} ${ARGN})
  # link the Google test infrastructure and a default main function to the test executable.
  target_link_libraries(${testname} PRIVATE MQT::${PROJECT_NAME}_${library}_lib gmock gtest_main)
  # discover tests
  gtest_discover_tests(
    ${testname}
    WORKING_DIRECTORY ${PROJECT_DIR}
    PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}")
  set_target_properties(${testname} PROPERTIES FOLDER tests)

  add_custom_command(
    TARGET ${testname}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${testname}>/${testname}
            ${CMAKE_BINARY_DIR}/${testname}
    COMMENT "Creating symlink for ${testname} in build directory"
    VERBATIM)
endmacro()

package_add_test(
  ${PROJECT_NAME}_heuristic_test heuristic ${CMAKE_CURRENT_SOURCE_DIR}/test_heuristic.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_general.cpp ${CMAKE_CURRENT_SOURCE_DIR}/test_architecture.cpp)

add_custom_command(
  TARGET ${PROJECT_NAME}_heuristic_test
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../examples
          $<TARGET_FILE_DIR:${PROJECT_NAME}_heuristic_test>/examples
  COMMAND ${CMAKE_COMMAND} -E create_symlink
          $<TARGET_FILE_DIR:${PROJECT_NAME}_heuristic_test>/examples ${CMAKE_BINARY_DIR}/examples
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../extern/architectures
          $<TARGET_FILE_DIR:${PROJECT_NAME}_heuristic_test>/architectures
  COMMAND
    ${CMAKE_COMMAND} -E create_symlink
    $<TARGET_FILE_DIR:${PROJECT_NAME}_heuristic_test>/architectures
    ${CMAKE_BINARY_DIR}/architectures
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../extern/calibration
          $<TARGET_FILE_DIR:${PROJECT_NAME}_heuristic_test>/calibration
  COMMAND
    ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${PROJECT_NAME}_heuristic_test>/calibration
    ${CMAKE_BINARY_DIR}/calibration
  COMMENT "Copying files and creating symlinks"
  VERBATIM)

if(TARGET ${PROJECT_NAME}_exact_lib)
  package_add_test(${PROJECT_NAME}_exact_test exact ${CMAKE_CURRENT_SOURCE_DIR}/test_exact.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_encodings.cpp)
endif()

if(TARGET ${PROJECT_NAME}_cliffordsynthesis_lib)
  configure_file(${CMAKE_SOURCE_DIR}/test/cliffordsynthesis/tableaus.json
                 ${CMAKE_CURRENT_BINARY_DIR}/cliffordsynthesis/tableaus.json COPYONLY)
  configure_file(${CMAKE_SOURCE_DIR}/test/cliffordsynthesis/circuits.json
                 ${CMAKE_CURRENT_BINARY_DIR}/cliffordsynthesis/circuits.json COPYONLY)
  package_add_test(
    ${PROJECT_NAME}_cliffordsynthesis_test cliffordsynthesis
    ${CMAKE_CURRENT_SOURCE_DIR}/test_tableau.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cliffordsynthesis/test_synthesis.cpp)
endif()
