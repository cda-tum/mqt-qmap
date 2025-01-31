if(TARGET MQT::AZAC)
    file(GLOB AZAC_TEST_SOURCES *.cpp)
    package_add_test(mqt-qmap-azac-test MQT::AZAC ${AZAC_TEST_SOURCES})
endif()
