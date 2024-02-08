set(MQT_QMAP_NA_TARGET_NAME ${MQT_QMAP_TARGET_NAME}-na)

if(NOT TARGET ${MQT_QMAP_NA_TARGET_NAME})
  file(GLOB_RECURSE NA_HEADERS ${MQT_QMAP_INCLUDE_BUILD_DIR}/na/*.hpp)
  file(GLOB_RECURSE NA_SOURCES *.cpp)

  add_library(
    ${MQT_QMAP_NA_TARGET_NAME}
    ${NA_HEADERS}
    ${NA_SOURCES})

  target_include_directories(${MQT_QMAP_NA_TARGET_NAME} PUBLIC $<BUILD_INTERFACE:${MQT_QMAP_INCLUDE_BUILD_DIR}/na>)

  target_link_libraries(${MQT_QMAP_NA_TARGET_NAME} PUBLIC MQT::Core nlohmann_json::nlohmann_json)
  target_link_libraries(${MQT_QMAP_NA_TARGET_NAME} PRIVATE MQT::ProjectOptions MQT::ProjectWarnings)
endif()
