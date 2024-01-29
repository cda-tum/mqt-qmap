pybind11_add_module(
  pyqmap
  # Prefer thin LTO if available
  THIN_LTO
  # Optimize the bindings for size
  OPT_SIZE
  # Source code goes here
  bindings.cpp)
target_compile_definitions(pyqmap PRIVATE Z3_FOUND)
target_link_libraries(
  pyqmap
  PRIVATE MQT::QMapExact
          MQT::QMapHeuristic
          MQT::QMapCliffordSynthesis
          MQT::CorePython
          MQT::ProjectOptions
          MQT::ProjectWarnings
          pybind11_json)

# Install directive for scikit-build-core
install(
  TARGETS pyqmap
  DESTINATION .
  COMPONENT mqt-qmap_PythonModule)
