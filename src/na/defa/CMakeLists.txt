set(MQT_QMAP_DEFA_TARGET_NAME ${MQT_QMAP_TARGET_NAME}-defa)

if(NOT TARGET ${MQT_QMAP_DEFA_TARGET_NAME})
  file(GLOB DEFA_HEADERS ${MQT_QMAP_INCLUDE_BUILD_DIR}/na/defa/*.hpp)
  file(GLOB DEFA_SOURCES *.cpp)

  add_library(${MQT_QMAP_DEFA_TARGET_NAME} ${DEFA_HEADERS} ${DEFA_SOURCES})

  target_include_directories(${MQT_QMAP_DEFA_TARGET_NAME}
                             PUBLIC $<BUILD_INTERFACE:${MQT_QMAP_INCLUDE_BUILD_DIR}>)

  target_link_libraries(
    ${MQT_QMAP_DEFA_TARGET_NAME}
    PUBLIC MQT::QMapNA
    PRIVATE MQT::ProjectOptions MQT::ProjectWarnings)

  add_library(MQT::DEFA ALIAS ${MQT_QMAP_DEFA_TARGET_NAME})
endif()
