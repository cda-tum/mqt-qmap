if(TARGET MQT::NASP)
  file(GLOB NA_TEST_SOURCES *.cpp)
  package_add_test(mqt-qmap-nasp-test MQT::NASP ${NA_TEST_SOURCES})
endif()
