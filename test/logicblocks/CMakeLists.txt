if(TARGET MQT::LogicBlocks)
  file(GLOB_RECURSE LOGICBLOCKS_TEST_SOURCES *.cpp)
  package_add_test(mqt-logic-blocks-test MQT::LogicBlocks ${LOGICBLOCKS_TEST_SOURCES})
endif()
