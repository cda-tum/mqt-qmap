# Declare all external dependencies and make sure that they are available.

include(FetchContent)
set(FETCH_PACKAGES "")

# search for Z3
find_package(Z3 4.8.15)
if(NOT Z3_FOUND)
  message(
    WARNING "Did not find Z3. Exact mapper and Clifford synthesis libraries will not be available")
endif()

if(BUILD_MQT_QMAP_BINDINGS)
  # Manually detect the installed mqt-core package.
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m mqt.core --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE mqt-core_DIR
    ERROR_QUIET)

  # Add the detected directory to the CMake prefix path.
  if(mqt-core_DIR)
    list(APPEND CMAKE_PREFIX_PATH "${mqt-core_DIR}")
    message(STATUS "Found mqt-core package: ${mqt-core_DIR}")
  endif()

  if(NOT SKBUILD)
    # Manually detect the installed pybind11 package.
    execute_process(
      COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
      OUTPUT_STRIP_TRAILING_WHITESPACE
      OUTPUT_VARIABLE pybind11_DIR)

    # Add the detected directory to the CMake prefix path.
    list(APPEND CMAKE_PREFIX_PATH "${pybind11_DIR}")
  endif()

  # add pybind11 library
  find_package(pybind11 2.13.5 CONFIG REQUIRED)
endif()

# cmake-format: off
set(MQT_CORE_VERSION 3.0.0
    CACHE STRING "MQT Core version")
set(MQT_CORE_REV "b7c3cfc9b94df6a7a5f9ba217ffcbf58e56dbc67"
    CACHE STRING "MQT Core identifier (tag, branch or commit hash)")
set(MQT_CORE_REPO_OWNER "cda-tum"
    CACHE STRING "MQT Core repository owner (change when using a fork)")
# cmake-format: on
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
  FetchContent_Declare(
    mqt-core
    GIT_REPOSITORY https://github.com/${MQT_CORE_REPO_OWNER}/mqt-core.git
    GIT_TAG ${MQT_CORE_REV}
    FIND_PACKAGE_ARGS ${MQT_CORE_VERSION})
  list(APPEND FETCH_PACKAGES mqt-core)
else()
  find_package(mqt-core ${MQT_CORE_VERSION} QUIET)
  if(NOT mqt-core_FOUND)
    FetchContent_Declare(
      mqt-core
      GIT_REPOSITORY https://github.com/${MQT_CORE_REPO_OWNER}/mqt-core.git
      GIT_TAG ${MQT_CORE_REV})
    list(APPEND FETCH_PACKAGES mqt-core)
  endif()
endif()

set(JSON_VERSION
    3.11.3
    CACHE STRING "nlohmann_json version")
set(JSON_URL https://github.com/nlohmann/json/releases/download/v${JSON_VERSION}/json.tar.xz)
set(JSON_SystemInclude
    ON
    CACHE INTERNAL "Treat the library headers like system headers")
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
  FetchContent_Declare(nlohmann_json URL ${JSON_URL} FIND_PACKAGE_ARGS ${JSON_VERSION})
  list(APPEND FETCH_PACKAGES nlohmann_json)
else()
  find_package(nlohmann_json ${JSON_VERSION} QUIET)
  if(NOT nlohmann_json_FOUND)
    FetchContent_Declare(nlohmann_json URL ${JSON_URL})
    list(APPEND FETCH_PACKAGES nlohmann_json)
  endif()
endif()

set(PLOG_VERSION
    1.1.10
    CACHE STRING "Plog version")
set(PLOG_URL https://github.com/SergiusTheBest/plog/archive/refs/tags/${PLOG_VERSION}.tar.gz)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
  FetchContent_Declare(plog URL ${PLOG_URL} FIND_PACKAGE_ARGS ${PLOG_VERSION})
  list(APPEND FETCH_PACKAGES plog)
else()
  find_package(plog ${PLOG_VERSION} QUIET)
  if(NOT plog_FOUND)
    FetchContent_Declare(plog URL ${PLOG_URL})
    list(APPEND FETCH_PACKAGES plog)
  endif()
endif()

if(BUILD_MQT_QMAP_TESTS)
  set(gtest_force_shared_crt
      ON
      CACHE BOOL "" FORCE)
  set(GTEST_VERSION
      1.14.0
      CACHE STRING "Google Test version")
  set(GTEST_URL https://github.com/google/googletest/archive/refs/tags/v${GTEST_VERSION}.tar.gz)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
    FetchContent_Declare(googletest URL ${GTEST_URL} FIND_PACKAGE_ARGS ${GTEST_VERSION} NAMES GTest)
    list(APPEND FETCH_PACKAGES googletest)
  else()
    find_package(googletest ${GTEST_VERSION} QUIET NAMES GTest)
    if(NOT googletest_FOUND)
      FetchContent_Declare(googletest URL ${GTEST_URL})
      list(APPEND FETCH_PACKAGES googletest)
    endif()
  endif()
endif()

if(BUILD_MQT_QMAP_BINDINGS)
  # add pybind11_json library
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
    FetchContent_Declare(
      pybind11_json
      GIT_REPOSITORY https://github.com/pybind/pybind11_json
      FIND_PACKAGE_ARGS)
    list(APPEND FETCH_PACKAGES pybind11_json)
  else()
    find_package(pybind11_json QUIET)
    if(NOT pybind11_json_FOUND)
      FetchContent_Declare(pybind11_json GIT_REPOSITORY https://github.com/pybind/pybind11_json)
      list(APPEND FETCH_PACKAGES pybind11_json)
    endif()
  endif()
endif()

# Make all declared dependencies available.
FetchContent_MakeAvailable(${FETCH_PACKAGES})
