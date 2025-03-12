if(TARGET MQT::AZAC)
  file(GLOB AZAC_TEST_SOURCES *.cpp)
  package_add_test(mqt-qmap-azac-test MQT::AZAC ${AZAC_TEST_SOURCES})
  target_link_libraries(mqt-qmap-azac-test PRIVATE MQT::CoreCircuitOptimizer MQT::CoreQASM)
  # also provide individual targets for every file
  foreach(TEST_SOURCE ${AZAC_TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    string(TOLOWER ${TEST_NAME} TEST_NAME)
    string(REPLACE "_" "-" TEST_NAME ${TEST_NAME})
    package_add_test(mqt-qmap-azac-${TEST_NAME} MQT::AZAC ${TEST_SOURCE})
    target_link_libraries(mqt-qmap-azac-${TEST_NAME} PRIVATE MQT::CoreCircuitOptimizer
                                                             MQT::CoreQASM)
  endforeach()
endif()
