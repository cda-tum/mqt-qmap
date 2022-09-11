# search for Boost program_options
set(BOOST_USE_MULTITHREADED ON)
set(BOOST_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.50 COMPONENTS program_options)

if(NOT TARGET Boost::program_options)
  add_library(Boost::program_options IMPORTED INTERFACE)
  set_property(TARGET Boost::program_options PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                                      ${Boost_INCLUDE_DIR})
  set_property(TARGET Boost::program_options PROPERTY INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()

# macro to add a command line executable
macro(ADD_QMAP_EXECUTABLE appname)
  add_executable(${PROJECT_NAME}_${appname} ${CMAKE_CURRENT_SOURCE_DIR}/${appname}_app.cpp)
  target_link_libraries(${PROJECT_NAME}_${appname} PUBLIC MQT::${PROJECT_NAME}_${appname}_lib
                                                          Boost::program_options)
  add_custom_command(
    TARGET ${PROJECT_NAME}_${appname}
    POST_BUILD
    COMMAND
      ${CMAKE_COMMAND} -E create_symlink
      $<TARGET_FILE_DIR:${PROJECT_NAME}_${appname}>/${PROJECT_NAME}_${appname}
      ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_${appname}
    COMMENT "Creating symlink for ${PROJECT_NAME}_${appname}"
    VERBATIM)
endmacro()

if(Boost_FOUND)
  message(STATUS "Boost version: ${Boost_VERSION}")
  if(TARGET ${PROJECT_NAME}_exact_lib)
    add_qmap_executable(exact)
  endif()

	add_qmap_executable(heuristic)

	add_qmap_executable(cliffordsynthesis)
else()
  message(STATUS "Did not find Boost! Commandline applications will not be an available target!")
endif()
