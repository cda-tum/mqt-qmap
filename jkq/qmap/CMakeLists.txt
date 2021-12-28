pybind11_add_module(py${PROJECT_NAME} bindings.cpp)
target_link_libraries(py${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_heuristic_lib JKQ::qfr_python pybind11_json)
if (Z3_FOUND)
	target_compile_definitions(py${PROJECT_NAME} PUBLIC Z3_FOUND)
	target_link_libraries(py${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_exact_lib)
endif()

# LTO is disabled since a bug in manylinux2014 prevents the project from linking with it
if (${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux")
	set_target_properties(${PROJECT_NAME}_exact_lib PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
	set_target_properties(${PROJECT_NAME}_heuristic_lib PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
	set_target_properties(qfr PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
	set_target_properties(DDPackage PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
endif ()
