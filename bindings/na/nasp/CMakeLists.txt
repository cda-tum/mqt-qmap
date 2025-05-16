# Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
# Copyright (c) 2025 Munich Quantum Software Company GmbH
# All rights reserved.
#
# SPDX-License-Identifier: MIT
#
# Licensed under the MIT License

# set target name
set(TARGET_NAME ${MQT_QMAP_TARGET_NAME}-na-nasp-bindings)
# replace '-' with '_' because python imports cannot contain '-'
string(REPLACE "-" "_" TARGET_NAME ${TARGET_NAME})

# declare the Python module
pybind11_add_module(
  # Name of the extension
  ${TARGET_NAME}
  # Prefer thin LTO if available
  THIN_LTO
  # Optimize the bindings for size
  OPT_SIZE
  # Source code goes here
  nasp.cpp)
target_compile_definitions(${TARGET_NAME} PRIVATE MQT_QMAP_MODULE_NAME=${TARGET_NAME})
target_link_libraries(${TARGET_NAME} PRIVATE MQT::NASP MQT::CoreQASM MQT::ProjectOptions
                                             MQT::ProjectWarnings pybind11_json)

# Install directive for scikit-build-core
install(
  TARGETS ${TARGET_NAME}
  DESTINATION ./na/state_preparation
  COMPONENT ${MQT_QMAP_TARGET_NAME}_Python)

# Install the Python stub files in editable mode for better IDE support
if(SKBUILD_STATE STREQUAL "editable")
  file(GLOB PYI_FILES ${PROJECT_SOURCE_DIR}/python/mqt/qmap/na/state_preparation/*.pyi)
  install(
    FILES ${PYI_FILES}
    DESTINATION ./na/state_preparation
    COMPONENT ${MQT_QMAP_TARGET_NAME}_Python)
endif()
