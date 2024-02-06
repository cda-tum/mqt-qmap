# neutral atoms mapper library
add_library(
  ${PROJECT_NAME}_na_lib
  na/NeutralAtomMapper.cpp ${PROJECT_SOURCE_DIR}/include/na/NeutralAtomMapper.hpp
  ${PROJECT_SOURCE_DIR}/include/na/Architecture.hpp ${PROJECT_SOURCE_DIR}/include/Mapper.hpp
  na/Architecture.cpp Mapper.cpp)
add_internal_library(${PROJECT_NAME}_na_lib NeutralAtomMapper)