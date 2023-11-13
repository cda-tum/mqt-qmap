# add LogicBlocks target
set(BUILD_LB_TESTS
    OFF
    CACHE BOOL "Build LogicBlocks tests")
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/LogicBlocks" "extern/LogicBlocks" EXCLUDE_FROM_ALL)
# add MQT::Core target
set(BUILD_MQT_CORE_TESTS
    OFF
    CACHE BOOL "Build MQT Core tests")
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/mqt-core" "extern/mqt-core" EXCLUDE_FROM_ALL)

# macro to add core libraries
macro(ADD_INTERNAL_LIBRARY libname)
  set(lib ${libname})

  # set include directories
  target_include_directories(${lib} PUBLIC ${PROJECT_SOURCE_DIR}/include
                                           ${PROJECT_BINARY_DIR}/include)

  # add MQT::qfr library
  target_link_libraries(${lib} PUBLIC MQT::Core)

  # add MQT alias
  add_library(MQT::${lib} ALIAS ${lib})

endmacro()

# macro to add mapping libraries
macro(ADD_QMAP_LIBRARY libname srcfile)
  set(lib ${PROJECT_NAME}_${libname}_lib)

  # main project library
  add_library(
    ${lib}
    ${libname}/${srcfile}.cpp
    ${PROJECT_SOURCE_DIR}/include/${libname}/${srcfile}.hpp
    ${PROJECT_SOURCE_DIR}/include/Architecture.hpp
    ${PROJECT_SOURCE_DIR}/include/configuration
    ${PROJECT_SOURCE_DIR}/include/heuristic/UniquePriorityQueue.hpp
    ${PROJECT_SOURCE_DIR}/include/Mapper.hpp
    ${PROJECT_SOURCE_DIR}/include/MappingResults.hpp
    ${PROJECT_SOURCE_DIR}/include/utils.hpp
    Architecture.cpp
    configuration/Configuration.cpp
    Mapper.cpp
    utils.cpp)

  add_internal_library(${lib} ${srcfile})
endmacro()

# macro to add synthesis libraries
macro(ADD_SYNTHESIS_LIBRARY libname srcfile)
  set(lib ${PROJECT_NAME}_${libname}_lib)

  # main project library
  add_library(
    ${lib}
    ${libname}/${srcfile}.cpp
    ${PROJECT_SOURCE_DIR}/include/${libname}/${srcfile}.hpp
    ${PROJECT_SOURCE_DIR}/include/cliffordsynthesis/Configuration.hpp
    ${PROJECT_SOURCE_DIR}/include/cliffordsynthesis/encoding/GateEncoder.hpp
    ${PROJECT_SOURCE_DIR}/include/cliffordsynthesis/encoding/MultiGateEncoder.hpp
    ${PROJECT_SOURCE_DIR}/include/cliffordsynthesis/encoding/ObjectiveEncoder.hpp
    ${PROJECT_SOURCE_DIR}/include/cliffordsynthesis/encoding/SATEncoder.hpp
    ${PROJECT_SOURCE_DIR}/include/cliffordsynthesis/encoding/SingleGateEncoder.hpp
    ${PROJECT_SOURCE_DIR}/include/cliffordsynthesis/encoding/TableauEncoder.hpp
    ${PROJECT_SOURCE_DIR}/include/cliffordsynthesis/Results.hpp
    ${PROJECT_SOURCE_DIR}/include/cliffordsynthesis/Tableau.hpp
    ${PROJECT_SOURCE_DIR}/include/cliffordsynthesis/TargetMetric.hpp
    ${PROJECT_SOURCE_DIR}/include/utils.hpp
    ${PROJECT_SOURCE_DIR}/include/cliffordsynthesis/GateSet.hpp
    ${PROJECT_SOURCE_DIR}/include/cliffordsynthesis/Utils.hpp
    cliffordsynthesis/encoding/GateEncoder.cpp
    cliffordsynthesis/encoding/MultiGateEncoder.cpp
    cliffordsynthesis/encoding/ObjectiveEncoder.cpp
    cliffordsynthesis/encoding/SATEncoder.cpp
    cliffordsynthesis/encoding/SingleGateEncoder.cpp
    cliffordsynthesis/encoding/TableauEncoder.cpp
    cliffordsynthesis/Tableau.cpp
    cliffordsynthesis/Utils.cpp
    utils.cpp)

  add_internal_library(${lib} ${srcfile})
endmacro()

# exact mapper project library
add_qmap_library(exact ExactMapper)

# add LogicBlocks library
target_link_libraries(${PROJECT_NAME}_exact_lib PUBLIC logicblocks::Logic logicblocks::Encodings)

add_synthesis_library(cliffordsynthesis CliffordSynthesizer)

# add LogicBlocks library
target_link_libraries(${PROJECT_NAME}_cliffordsynthesis_lib PUBLIC logicblocks::Logic
                                                                   logicblocks::Encodings)

# heuristic mapper project library
add_qmap_library(heuristic HeuristicMapper)

if(BUILD_MQT_QMAP_BINDINGS)
  add_subdirectory(python)
endif()
