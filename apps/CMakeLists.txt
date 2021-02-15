if(NOT TARGET Boost::program_options)
	add_library(Boost::program_options IMPORTED INTERFACE)
	set_property(TARGET Boost::program_options PROPERTY
	             INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
	set_property(TARGET Boost::program_options PROPERTY
	             INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()

macro(add_qmap_executable APPNAME)
	add_executable(${PROJECT_NAME}_${APPNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${APPNAME}_app.cpp)
	target_link_libraries(${PROJECT_NAME}_${APPNAME} PUBLIC JKQ::${PROJECT_NAME}_${APPNAME}_lib Boost::program_options)
	add_custom_command(TARGET ${PROJECT_NAME}_${APPNAME}
	                   POST_BUILD
	                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${PROJECT_NAME}_${APPNAME}>/${PROJECT_NAME}_${APPNAME} ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_${APPNAME}
	                   COMMENT "Creating symlink for ${PROJECT_NAME}_${APPNAME}"
	                   VERBATIM)
endmacro()

if (Boost_FOUND)
	if (TARGET ${PROJECT_NAME}_exact_lib)
		add_qmap_executable(exact)
	endif()

	add_qmap_executable(heuristic)
else()
	message(STATUS "Did not find Boost! Commandline applications will not be an available target!")
endif ()

