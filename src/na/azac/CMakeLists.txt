set(MQT_QMAP_AZAC_TARGET_NAME ${MQT_QMAP_TARGET_NAME}-azac)

if(NOT TARGET ${MQT_QMAP_AZAC_TARGET_NAME})
  file(GLOB AZAC_HEADERS ${MQT_QMAP_INCLUDE_BUILD_DIR}/na/azac/*.hpp)
  file(GLOB AZAC_SOURCES *.cpp)

#  add_library(${MQT_QMAP_AZAC_TARGET_NAME} ${AZAC_HEADERS} ${AZAC_SOURCES})
  add_library(${MQT_QMAP_AZAC_TARGET_NAME}
          ${MQT_QMAP_INCLUDE_BUILD_DIR}/na/azac/Utils.hpp
          ${MQT_QMAP_INCLUDE_BUILD_DIR}/na/azac/Architecture.hpp
          Architecture.cpp
  )

  target_include_directories(${MQT_QMAP_AZAC_TARGET_NAME}
                             PUBLIC $<BUILD_INTERFACE:${MQT_QMAP_INCLUDE_BUILD_DIR}>)

  target_link_libraries(${MQT_QMAP_AZAC_TARGET_NAME}
          PUBLIC MQT::CoreIR nlohmann_json::nlohmann_json
          PRIVATE MQT::ProjectOptions MQT::ProjectWarnings)

  add_library(MQT::AZAC ALIAS ${MQT_QMAP_AZAC_TARGET_NAME})
endif()
