if(TARGET MQT::ZONAR)
  file(GLOB NA_TEST_SOURCES *.cpp)
  package_add_test(mqt-qmap-zonar-test MQT::ZONAR ${NA_TEST_SOURCES})
endif()
