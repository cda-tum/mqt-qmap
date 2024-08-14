set(MQT_QMAP_TARGET_NAME mqt-qmap)

# macro to add core libraries
macro(ADD_INTERNAL_LIBRARY libname)
  set(lib ${libname})

  # set include directories
  target_include_directories(${lib} PUBLIC $<BUILD_INTERFACE:${MQT_QMAP_INCLUDE_BUILD_DIR}>)

  # add MQT::qfr library
  target_link_libraries(
    ${lib}
    PUBLIC MQT::CoreIR MQT::CoreCircuitOptimizer nlohmann_json::nlohmann_json
    PRIVATE MQT::ProjectOptions MQT::ProjectWarnings)
endmacro()

# macro to add mapping libraries
macro(ADD_QMAP_LIBRARY libname srcfile)
  set(lib ${MQT_QMAP_TARGET_NAME}-${libname})

  # main project library
  add_library(
    ${lib}
    ${libname}/${srcfile}.cpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/${libname}/${srcfile}.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/Architecture.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/configuration
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/DataLogger.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/heuristic/UniquePriorityQueue.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/Mapper.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/MappingResults.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/utils.hpp
    Architecture.cpp
    configuration/Configuration.cpp
    DataLogger.cpp
    Mapper.cpp
    utils.cpp)

  add_internal_library(${lib})
endmacro()

# macro to add synthesis libraries
macro(ADD_SYNTHESIS_LIBRARY libname srcfile)
  set(lib ${MQT_QMAP_TARGET_NAME}-${libname})

  # main project library
  add_library(
    ${lib}
    ${libname}/${srcfile}.cpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/${libname}/${srcfile}.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/cliffordsynthesis/Configuration.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/cliffordsynthesis/encoding/GateEncoder.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/cliffordsynthesis/encoding/MultiGateEncoder.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/cliffordsynthesis/encoding/ObjectiveEncoder.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/cliffordsynthesis/encoding/SATEncoder.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/cliffordsynthesis/encoding/SingleGateEncoder.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/cliffordsynthesis/encoding/TableauEncoder.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/cliffordsynthesis/Results.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/cliffordsynthesis/Tableau.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/cliffordsynthesis/TargetMetric.hpp
    ${MQT_QMAP_INCLUDE_BUILD_DIR}/utils.hpp
    cliffordsynthesis/encoding/GateEncoder.cpp
    cliffordsynthesis/encoding/MultiGateEncoder.cpp
    cliffordsynthesis/encoding/ObjectiveEncoder.cpp
    cliffordsynthesis/encoding/SATEncoder.cpp
    cliffordsynthesis/encoding/SingleGateEncoder.cpp
    cliffordsynthesis/encoding/TableauEncoder.cpp
    cliffordsynthesis/Tableau.cpp
    utils.cpp)

  add_internal_library(${lib})
  target_link_libraries(${lib} PUBLIC MQT::CoreNA)
endmacro()

# macro to add hybrid neutral atom mapping libraries
macro(ADD_HYBRIDMAP_LIBRARY libname srcfile)
  set(lib ${MQT_QMAP_TARGET_NAME}-${libname})

  # main project library
  add_library(
    ${lib}
    ${libname}/${srcfile}.cpp
    ${PROJECT_SOURCE_DIR}/include/${libname}/${srcfile}.hpp
    ${PROJECT_SOURCE_DIR}/include/${libname}/NeutralAtomUtils.hpp
    ${PROJECT_SOURCE_DIR}/include/${libname}/NeutralAtomScheduler.hpp
    ${PROJECT_SOURCE_DIR}/include/${libname}/NeutralAtomArchitecture.hpp
    ${PROJECT_SOURCE_DIR}/include/${libname}/HardwareQubits.hpp
    ${PROJECT_SOURCE_DIR}/include/${libname}/Mapping.hpp
    ${PROJECT_SOURCE_DIR}/include/${libname}/NeutralAtomDefinitions.hpp
    ${PROJECT_SOURCE_DIR}/include/${libname}/MoveToAodConverter.hpp
    ${PROJECT_SOURCE_DIR}/include/${libname}/NeutralAtomLayer.hpp
    ${PROJECT_SOURCE_DIR}/include/${libname}/HybridAnimation.hpp
    hybridmap/NeutralAtomUtils.cpp
    hybridmap/NeutralAtomScheduler.cpp
    hybridmap/NeutralAtomArchitecture.cpp
    hybridmap/HardwareQubits.cpp
    hybridmap/Mapping.cpp
    hybridmap/MoveToAodConverter.cpp
    hybridmap/HybridAnimation.cpp
    hybridmap/NeutralAtomLayer.cpp)

  add_internal_library(${lib})
endmacro()

# exact mapper project library
if(Z3_FOUND)
  # add LogicBlocks library
  add_subdirectory(logicblocks)

  add_qmap_library(exact ExactMapper)
  target_link_libraries(${MQT_QMAP_TARGET_NAME}-exact PRIVATE MQT::LogicBlocks)

  add_synthesis_library(cliffordsynthesis CliffordSynthesizer)
  target_link_libraries(
    ${MQT_QMAP_TARGET_NAME}-cliffordsynthesis
    PUBLIC plog::plog
    PRIVATE MQT::LogicBlocks)

  # add MQT alias targets
  add_library(MQT::QMapExact ALIAS ${MQT_QMAP_TARGET_NAME}-exact)
  add_library(MQT::QMapCliffordSynthesis ALIAS ${MQT_QMAP_TARGET_NAME}-cliffordsynthesis)
endif()

# heuristic mapper project library
add_qmap_library(heuristic HeuristicMapper)

# hybrid neutral atom mapper project library
add_hybridmap_library(hybridmap HybridNeutralAtomMapper)

# add MQT alias targets
add_library(MQT::QMapHeuristic ALIAS ${MQT_QMAP_TARGET_NAME}-heuristic)
add_library(MQT::QMapHybrid ALIAS ${MQT_QMAP_TARGET_NAME}-hybridmap)

add_subdirectory(na)

if(BUILD_MQT_QMAP_BINDINGS)
  add_subdirectory(python)
endif()
