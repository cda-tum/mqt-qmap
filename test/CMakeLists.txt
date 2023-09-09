if(NOT TARGET gtest OR NOT TARGET gmock)
  # Prevent overriding the parent project's compiler/linker settings on Windows
  set(gtest_force_shared_crt # cmake-lint: disable=C0103
      ON
      CACHE BOOL "" FORCE)
  add_subdirectory("${PROJECT_SOURCE_DIR}/extern/mqt-core/extern/googletest"
                   "extern/mqt-core/extern/googletest" EXCLUDE_FROM_ALL)
  set_target_properties(gtest gtest_main gmock gmock_main PROPERTIES FOLDER extern)
endif()

package_add_test(
  ${PROJECT_NAME}_heuristic_test ${PROJECT_NAME}_heuristic_lib
  ${CMAKE_CURRENT_SOURCE_DIR}/test_heuristic.cpp ${CMAKE_CURRENT_SOURCE_DIR}/test_general.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_architecture.cpp)

if(TARGET ${PROJECT_NAME}_exact_lib)
  package_add_test(
    ${PROJECT_NAME}_exact_test ${PROJECT_NAME}_exact_lib ${CMAKE_CURRENT_SOURCE_DIR}/test_exact.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_encodings.cpp)
endif()

if(TARGET ${PROJECT_NAME}_cliffordsynthesis_lib)
  configure_file(${CMAKE_SOURCE_DIR}/test/cliffordsynthesis/tableaus.json
                 ${CMAKE_CURRENT_BINARY_DIR}/cliffordsynthesis/tableaus.json COPYONLY)
  configure_file(${CMAKE_SOURCE_DIR}/test/cliffordsynthesis/circuits.json
                 ${CMAKE_CURRENT_BINARY_DIR}/cliffordsynthesis/circuits.json COPYONLY)
  package_add_test(
    ${PROJECT_NAME}_cliffordsynthesis_test ${PROJECT_NAME}_cliffordsynthesis_lib
    ${CMAKE_CURRENT_SOURCE_DIR}/test_tableau.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cliffordsynthesis/test_synthesis.cpp)
endif()
