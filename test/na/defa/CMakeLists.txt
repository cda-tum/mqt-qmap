if(TARGET MQT::DEFA)
  file(GLOB NA_TEST_SOURCES *.cpp)
  package_add_test(mqt-qmap-defa-test MQT::DEFA ${NA_TEST_SOURCES})
endif()
