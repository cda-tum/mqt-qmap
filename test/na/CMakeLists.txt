if(TARGET MQT::QMapNA)
  file(GLOB_RECURSE NA_TEST_SOURCES *.cpp)
  list(FILTER NA_TEST_SOURCES EXCLUDE REGEX "test_exe.cpp")
  package_add_test(mqt-qmap-na-test MQT::QMapNA ${NA_TEST_SOURCES})

  add_executable(mqt-qmap-na-test-exe test_exe.cpp)
  target_link_libraries(mqt-qmap-na-test-exe MQT::QMapNA MQT::ProjectOptions
    MQT::ProjectWarnings)
endif()