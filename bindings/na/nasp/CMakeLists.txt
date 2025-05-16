# Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
# Copyright (c) 2025 Munich Quantum Software Company GmbH
# All rights reserved.
#
# SPDX-License-Identifier: MIT
#
# Licensed under the MIT License

# set target name
set(TARGET_NAME ${MQT_QMAP_TARGET_NAME}-na-nasp-bindings)

# set the module name, this must be a valid Python module name
set(MODULE_NAME "state_preparation")

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

# the library name must be the same as the module name
set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME ${MODULE_NAME})
target_compile_definitions(${TARGET_NAME} PRIVATE MQT_QMAP_MODULE_NAME=${MODULE_NAME})

target_link_libraries(${TARGET_NAME} PRIVATE MQT::NASP MQT::CoreQASM MQT::ProjectOptions
                                             MQT::ProjectWarnings pybind11_json)

# Install directive for scikit-build-core
install(
  TARGETS ${TARGET_NAME}
  DESTINATION ./na
  COMPONENT ${MQT_QMAP_TARGET_NAME}_Python)
