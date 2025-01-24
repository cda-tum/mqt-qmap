if(TARGET MQT::NALAC)
  file(GLOB NA_TEST_SOURCES *.cpp)
  package_add_test(mqt-qmap-nalac-test MQT::NALAC ${NA_TEST_SOURCES})
endif()
