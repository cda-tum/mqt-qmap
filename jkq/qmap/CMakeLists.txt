pybind11_add_module(py${PROJECT_NAME} bindings.cpp)
target_link_libraries(py${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_exact_lib ${PROJECT_NAME}_heuristic_lib JKQ::pyqfr pybind11_json)
set_target_properties(${PROJECT_NAME}_exact_lib PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
set_target_properties(${PROJECT_NAME}_heuristic_lib PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

# LTO is disabled since a bug in manylinux2014 prevents the project from linking with it
if (${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux")
	set_target_properties(${PROJECT_NAME}_exact_lib PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
	set_target_properties(${PROJECT_NAME}_heuristic_lib PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
	set_target_properties(qfr PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
	set_target_properties(dd_package PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
endif ()
