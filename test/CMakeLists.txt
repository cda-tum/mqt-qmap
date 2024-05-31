package_add_test(
  mqt-qmap-heuristic-test MQT::QMapHeuristic ${CMAKE_CURRENT_SOURCE_DIR}/test_heuristic.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_general.cpp ${CMAKE_CURRENT_SOURCE_DIR}/test_architecture.cpp)

if(TARGET MQT::QMapExact)
  package_add_test(mqt-qmap-exact-test MQT::QMapExact ${CMAKE_CURRENT_SOURCE_DIR}/test_exact.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_encodings.cpp)
endif()

add_subdirectory(na)

if(TARGET MQT::QMapCliffordSynthesis)
  configure_file(${PROJECT_SOURCE_DIR}/test/cliffordsynthesis/tableaus.json
                 ${CMAKE_CURRENT_BINARY_DIR}/cliffordsynthesis/tableaus.json COPYONLY)
  configure_file(${PROJECT_SOURCE_DIR}/test/cliffordsynthesis/circuits.json
                 ${CMAKE_CURRENT_BINARY_DIR}/cliffordsynthesis/circuits.json COPYONLY)
  package_add_test(
    mqt-qmap-cliffordsynthesis-test MQT::QMapCliffordSynthesis
    ${CMAKE_CURRENT_SOURCE_DIR}/test_tableau.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cliffordsynthesis/test_synthesis.cpp)
endif()

if(TARGET MQT::LogicBlocks)
  package_add_test(mqt-logic-blocks-test MQT::LogicBlocks
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_logicblocks.cpp)
endif()

if(TARGET MQT::QMapHybrid)
  file(COPY ${PROJECT_SOURCE_DIR}/test/hybridmap DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  package_add_test(${PROJECT_NAME}_hybridmap_test MQT::QMapHybrid
                   ${CMAKE_CURRENT_SOURCE_DIR}/hybridmap/test_hybridmap.cpp)
endif()
