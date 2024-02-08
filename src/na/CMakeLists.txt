if(NOT TARGET ${MQT_CORE_TARGET_NAME}-na)
  file(GLOB_RECURSE NA_HEADERS ${MQT_CORE_INCLUDE_BUILD_DIR}/na/*.hpp)

  # add NA package library
  add_library(
    ${MQT_CORE_TARGET_NAME}-na
    ${NA_HEADERS}
    Architecture.cpp
    Configuration.cpp
    NeutralAtomMapper.cpp)

  target_link_libraries(${MQT_CORE_TARGET_NAME}-na PUBLIC MQT::Core)
  target_link_libraries(${MQT_CORE_TARGET_NAME}-na PRIVATE MQT::ProjectOptions
    MQT::ProjectWarnings)
endif()
