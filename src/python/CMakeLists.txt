if(NOT SKBUILD)
  # Manually detect the installed pybind11 package and import it into CMake.
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE pybind11_DIR)
  list(APPEND CMAKE_PREFIX_PATH "${pybind11_DIR}")
endif()

# Import pybind11 through CMake's find_package mechanism
find_package(pybind11 CONFIG REQUIRED)

# We are now ready to compile the actual extension module
pybind11_add_module(
  pyqmap
  # Prefer thin LTO if available
  THIN_LTO
  # Optimize the bindings for size
  OPT_SIZE
  # Source code goes here
  bindings.cpp)
target_compile_definitions(pyqmap PRIVATE Z3_FOUND)
target_link_libraries(pyqmap PRIVATE MQT::CorePython)
target_link_libraries(pyqmap PRIVATE qmap_exact_lib)
target_link_libraries(pyqmap PRIVATE qmap_heuristic_lib)
target_link_libraries(pyqmap PRIVATE qmap_cliffordsynthesis_lib)

# Install directive for scikit-build-core
install(
  TARGETS pyqmap
  DESTINATION .
  COMPONENT mqt-qmap_PythonModule)
