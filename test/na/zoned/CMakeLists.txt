if(TARGET MQT::QMapNAZoned)
  file(GLOB SOURCES *.cpp)
  file(GLOB CIRCUITS circuits/*.qasm)
  # make list of files a comma separated list of strings
  string(REPLACE ";" "\",\"" CIRCUITS "\"${CIRCUITS}\"")
  package_add_test(mqt-qmap-na-zoned-test MQT::QMapNAZoned ${SOURCES})
  target_link_libraries(mqt-qmap-na-zoned-test PRIVATE MQT::CoreCircuitOptimizer MQT::CoreQASM)
  target_compile_definitions(mqt-qmap-na-zoned-test PRIVATE TEST_CIRCUITS=${CIRCUITS})
  # also provide individual targets for every file
  foreach(TEST_SOURCE ${SOURCES})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    string(TOLOWER ${TEST_NAME} TEST_NAME)
    string(REPLACE "_" "-" TEST_NAME ${TEST_NAME})
    package_add_test(mqt-qmap-na-zoned-${TEST_NAME} MQT::QMapNAZoned ${TEST_SOURCE})
    target_link_libraries(mqt-qmap-na-zoned-${TEST_NAME} PRIVATE MQT::CoreCircuitOptimizer
                                                                 MQT::CoreQASM)
  endforeach()
  target_compile_definitions(mqt-qmap-na-zoned-test-compiler PRIVATE TEST_CIRCUITS=${CIRCUITS})
endif()
