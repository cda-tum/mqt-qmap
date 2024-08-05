if(TARGET MQT::NALAC)
  file(GLOB_RECURSE NA_TEST_SOURCES *.cpp)
  package_add_test(mqt-qmap-nalac-test MQT::NALAC ${NA_TEST_SOURCES})
endif()
