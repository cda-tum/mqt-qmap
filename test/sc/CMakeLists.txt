if(TARGET MQT::QMapSC)
  file(GLOB SC_TEST_SOURCES *.cpp)
  package_add_test(mqt-qmap-sc-test MQT::QMapSC ${SC_TEST_SOURCES})
endif()

add_subdirectory(heuristic)
add_subdirectory(exact)
