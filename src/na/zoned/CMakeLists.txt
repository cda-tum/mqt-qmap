set(MQT_QMAP_NA_ZONED_TARGET_NAME ${MQT_QMAP_TARGET_NAME}-na-zoned)

if(NOT TARGET ${MQT_QMAP_NA_ZONED_TARGET_NAME})
  file(GLOB HEADERS ${MQT_QMAP_INCLUDE_BUILD_DIR}/na/zoned/*.hpp)
  file(GLOB SOURCES *.cpp)

  add_library(${MQT_QMAP_NA_ZONED_TARGET_NAME} ${HEADERS} ${SOURCES})

  target_include_directories(${MQT_QMAP_NA_ZONED_TARGET_NAME}
                             PUBLIC $<BUILD_INTERFACE:${MQT_QMAP_INCLUDE_BUILD_DIR}>)

  target_link_libraries(
    ${MQT_QMAP_NA_ZONED_TARGET_NAME}
    PUBLIC MQT::CoreIR MQT::CoreNA nlohmann_json::nlohmann_json
    PRIVATE MQT::ProjectOptions MQT::ProjectWarnings)

  add_library(MQT::QMapNAZoned ALIAS ${MQT_QMAP_NA_ZONED_TARGET_NAME})
endif()
