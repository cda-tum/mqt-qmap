if(NOT TARGET Boost::program_options)
	add_library(Boost::program_options IMPORTED INTERFACE)
	set_property(TARGET Boost::program_options PROPERTY
	             INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
	set_property(TARGET Boost::program_options PROPERTY
	             INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()

macro(add_qmap_executable APPNAME)
	add_executable(${PROJECT_NAME}_${APPNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${APPNAME}_app.cpp)
	target_compile_features(${PROJECT_NAME}_${APPNAME} PUBLIC cxx_std_14)
	set_target_properties(${PROJECT_NAME}_${APPNAME} PROPERTIES CXX_EXTENSIONS OFF)
	target_link_libraries(${PROJECT_NAME}_${APPNAME} PUBLIC JKQ::${PROJECT_NAME}_${APPNAME}_lib Boost::program_options)
	add_custom_command(TARGET ${PROJECT_NAME}_${APPNAME}
	                   POST_BUILD
	                   COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_DIR:${PROJECT_NAME}_${APPNAME}>/${PROJECT_NAME}_${APPNAME} ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_${APPNAME}
	                   COMMENT "Creating symlink for ${PROJECT_NAME}_${APPNAME}"
	                   VERBATIM)
	if(MSVC)
		target_compile_options(${PROJECT_NAME}_${APPNAME} PUBLIC "/W4" "$<$<CONFIG:RELEASE>:/O2 /GL>")
		target_link_options(${PROJECT_NAME}_${APPNAME} PUBLIC "$<$<CONFIG:Release>:/LTCG>")
	else()
		target_compile_options(${PROJECT_NAME}_${APPNAME} PUBLIC -Wall $<$<CONFIG:DEBUG>:-g3 -Og -Wextra -Wpedantic -pedantic> $<$<CONFIG:RELEASE>:-O3 -mtune=native -march=native -DNDEBUG -g0>)
	endif()

	# add coverage compiler and linker flag if COVERAGE is set
	if (COVERAGE)
		target_compile_options(${PROJECT_NAME}_${APPNAME} PRIVATE --coverage)
		target_link_libraries(${PROJECT_NAME}_${APPNAME} PRIVATE --coverage)
	endif()
endmacro()

include(GNUInstallDirs)
if (Boost_FOUND)
	if (TARGET ${PROJECT_NAME}_exact_lib)
		add_qmap_executable(exact)
		install(TARGETS ${PROJECT_NAME}_exact
		        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

		set_target_properties(${PROJECT_NAME}_exact PROPERTIES EXPORT_NAME ${PROJECT_NAME}_exact)
	endif()

	add_qmap_executable(heuristic)
	install(TARGETS ${PROJECT_NAME}_heuristic
	        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

	set_target_properties(${PROJECT_NAME}_heuristic PROPERTIES EXPORT_NAME ${PROJECT_NAME}_heuristic)
else()
	message(STATUS "Did not find Boost! Commandline applications will not be an available target!")
endif ()

