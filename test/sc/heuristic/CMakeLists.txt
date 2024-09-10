if(TARGET MQT::QMapSCHeuristic)
  file(GLOB_RECURSE SC_HEURISTIC_TEST_SOURCES *.cpp)
  package_add_test(mqt-qmap-sc-heuristic-test MQT::QMapSCHeuristic ${SC_HEURISTIC_TEST_SOURCES})
endif()
