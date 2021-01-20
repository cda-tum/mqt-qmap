pybind11_add_module(py${PROJECT_NAME} bindings.cpp)
target_link_libraries(py${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_exact_lib ${PROJECT_NAME}_heuristic_lib JKQ::pyqfr pybind11_json)
set_target_properties(${PROJECT_NAME}_exact_lib PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
set_target_properties(${PROJECT_NAME}_heuristic_lib PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
