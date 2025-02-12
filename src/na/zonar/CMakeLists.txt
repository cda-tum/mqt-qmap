set(MQT_QMAP_ZONAR_TARGET_NAME ${MQT_QMAP_TARGET_NAME}-zonar)

if(NOT TARGET ${MQT_QMAP_ZONAR_TARGET_NAME})
  file(GLOB ZONAR_HEADERS ${MQT_QMAP_INCLUDE_BUILD_DIR}/na/zonar/*.hpp)
  file(GLOB ZONAR_SOURCES *.cpp)

  add_library(${MQT_QMAP_ZONAR_TARGET_NAME} ${ZONAR_HEADERS} ${ZONAR_SOURCES})

  target_include_directories(${MQT_QMAP_ZONAR_TARGET_NAME}
                             PUBLIC $<BUILD_INTERFACE:${MQT_QMAP_INCLUDE_BUILD_DIR}>)

  target_link_libraries(
    ${MQT_QMAP_ZONAR_TARGET_NAME}
    PUBLIC MQT::QMapNA
    PRIVATE MQT::ProjectOptions MQT::ProjectWarnings)

  add_library(MQT::ZONAR ALIAS ${MQT_QMAP_ZONAR_TARGET_NAME})
endif()
