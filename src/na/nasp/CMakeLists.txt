set(MQT_QMAP_NASP_TARGET_NAME ${MQT_QMAP_TARGET_NAME}-nasp)

if(NOT TARGET ${MQT_QMAP_NASP_TARGET_NAME})
  file(GLOB_RECURSE NASP_HEADERS ${MQT_QMAP_INCLUDE_BUILD_DIR}/na/nasp/*.hpp)
  file(GLOB_RECURSE NASP_SOURCES *.cpp)

  add_library(${MQT_QMAP_NASP_TARGET_NAME} ${NASP_HEADERS} ${NASP_SOURCES})

  target_include_directories(${MQT_QMAP_NASP_TARGET_NAME} PUBLIC $<BUILD_INTERFACE:${MQT_QMAP_INCLUDE_BUILD_DIR}/na/nasp>)

  target_link_libraries(${MQT_QMAP_NASP_TARGET_NAME}
          PUBLIC MQT::QMapNA nlohmann_json::nlohmann_json z3::z3lib yaml-cpp
          PRIVATE MQT::ProjectOptions MQT::ProjectWarnings)

  add_library(MQT::NASP ALIAS ${MQT_QMAP_NASP_TARGET_NAME})
endif()
