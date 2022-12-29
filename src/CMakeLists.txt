# add MQT::qfr target
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/qfr" "extern/qfr")
# add LogicBlocks target
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/LogicBlocks" "extern/LogicBlocks")

# macro to add core libraries
macro(ADD_INTERNAL_LIBRARY libname srcfile)
  set(lib ${libname})

  # set include directories
  target_include_directories(${lib}
                             PUBLIC $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include/>)

  # set required C++ standard and disable compiler specific extensions
  target_compile_features(${lib} PUBLIC cxx_std_17)
  set_target_properties(${lib} PROPERTIES CMAKE_CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

  # add MQT::qfr library
  target_link_libraries(${lib} PUBLIC MQT::qfr)

  # enable interprocedural optimization if it is supported
  include(CheckIPOSupported)
  check_ipo_supported(RESULT ipo_supported)
  if(ipo_supported)
    set_target_properties(${lib} PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()

  # set compiler flags depending on compiler
  if(MSVC)
    target_compile_options(${lib} PUBLIC /utf-8)
  else()
    target_compile_options(${lib} PUBLIC -Wall -Wextra $<$<CONFIG:DEBUG>:-g3 -ggdb>)
    if(BINDINGS)
      # adjust visibility settings for building Python bindings
      target_compile_options(${lib} PUBLIC -fvisibility=hidden)
    endif()
    if(NOT DEPLOY)
      # only include machine-specific optimizations when building for the host machine
      target_compile_options(${lib} PUBLIC -mtune=native)
      include(CheckCXXCompilerFlag)
      check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)
      if(HAS_MARCH_NATIVE)
        target_compile_options(${lib} PUBLIC -march=native)
      endif()
    endif()
  endif()

  if(GENERATE_POSITION_INDEPENDENT_CODE OR BINDINGS)
    include(CheckPIESupported)
    check_pie_supported()
    set_target_properties(${lib} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  endif()

  # add coverage compiler and linker flag to the library and all targets that link against it, if
  # COVERAGE is set
  if(COVERAGE)
    target_compile_options(${lib} PUBLIC --coverage)
    target_link_libraries(${lib} PUBLIC --coverage)
  endif()

  # add MQT alias
  add_library(MQT::${lib} ALIAS ${lib})
endmacro()
# macro to add core libraries
macro(ADD_QMAP_LIBRARY libname srcfile)
  set(lib ${PROJECT_NAME}_${libname}_lib)

  # main project library
  add_library(
    ${lib}
    ${CMAKE_CURRENT_SOURCE_DIR}/utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Architecture.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Mapper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/${libname}/${srcfile}.cpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/utils.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/Architecture.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/Mapper.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/MappingResults.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/configuration
    ${${PROJECT_NAME}_SOURCE_DIR}/include/${libname}/${srcfile}.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/heuristic/UniquePriorityQueue.hpp)

  add_internal_library(${lib} ${srcfile})
endmacro()
# macro to add core libraries
macro(ADD_SYNTHESIS_LIBRARY libname srcfile)
  set(lib ${PROJECT_NAME}_${libname}_lib)

  # main project library
  add_library(
    ${lib}
    ${CMAKE_CURRENT_SOURCE_DIR}/utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cliffordsynthesis/encoding/SATEncoder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cliffordsynthesis/encoding/TableauEncoder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cliffordsynthesis/encoding/GateEncoder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cliffordsynthesis/encoding/SingleGateEncoder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cliffordsynthesis/encoding/MultiGateEncoder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cliffordsynthesis/encoding/ObjectiveEncoder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cliffordsynthesis/Tableau.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/${libname}/${srcfile}.cpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/utils.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/cliffordsynthesis/Configuration.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/cliffordsynthesis/Results.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/cliffordsynthesis/encoding/SATEncoder.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/cliffordsynthesis/encoding/TableauEncoder.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/cliffordsynthesis/encoding/GateEncoder.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/cliffordsynthesis/encoding/SingleGateEncoder.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/cliffordsynthesis/encoding/MultiGateEncoder.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/cliffordsynthesis/encoding/ObjectiveEncoder.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/cliffordsynthesis/Tableau.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/cliffordsynthesis/TargetMetric.hpp
    ${${PROJECT_NAME}_SOURCE_DIR}/include/${libname}/${srcfile}.hpp)

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
# add Pthreads
target_link_libraries(${PROJECT_NAME}_cliffordsynthesis_lib PRIVATE Threads::Threads)

# heuristic mapper project library
add_qmap_library(heuristic HeuristicMapper)
