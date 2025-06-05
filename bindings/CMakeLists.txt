# Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
# Copyright (c) 2025 Munich Quantum Software Company GmbH
# All rights reserved.
#
# SPDX-License-Identifier: MIT
#
# Licensed under the MIT License

if(APPLE)
  set(BASEPOINT @loader_path)
else()
  set(BASEPOINT $ORIGIN)
endif()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
list(
  APPEND
  CMAKE_INSTALL_RPATH
  ${BASEPOINT}
  ${BASEPOINT}/${CMAKE_INSTALL_LIBDIR}
  ${BASEPOINT}/../core/${CMAKE_INSTALL_LIBDIR}
  ${BASEPOINT}/../core/lib
  ${BASEPOINT}/../core/lib64
  ${BASEPOINT}/../../core/${CMAKE_INSTALL_LIBDIR}
  ${BASEPOINT}/../../core/lib
  ${BASEPOINT}/../../core/lib64)

add_mqt_qmap_binding(
  pyqmap
  bindings.cpp
  LINK_LIBS
  MQT::QMapSCExact
  MQT::QMapSCHeuristic
  MQT::QMapCliffordSynthesis
  MQT::QMapHybrid
  MQT::NASP
  MQT::QMapNAZoned
  MQT::CoreQASM
  pybind11_json)
target_compile_definitions(pyqmap PRIVATE Z3_FOUND)

# add all modules
add_subdirectory(clifford_synthesis)
add_subdirectory(na)
add_subdirectory(sc)
