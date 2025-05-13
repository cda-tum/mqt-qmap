if(TARGET MQT::QMapNA)
  file(GLOB NA_TEST_SOURCES *.cpp)
  package_add_test(mqt-qmap-na-test MQT::QMapNA ${NA_TEST_SOURCES})
endif()

add_subdirectory(nalac)
add_subdirectory(nasp)
add_subdirectory(zoned)
