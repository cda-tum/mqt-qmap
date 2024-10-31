if(TARGET MQT::QMapSCExact)
  file(GLOB_RECURSE SC_EXACT_TEST_SOURCES *.cpp)
  package_add_test(mqt-qmap-sc-exact-test MQT::QMapSCExact ${SC_EXACT_TEST_SOURCES})
endif()
