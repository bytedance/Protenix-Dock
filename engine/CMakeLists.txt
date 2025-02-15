cmake_minimum_required(VERSION 3.12)
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.16.0")
  cmake_policy(SET CMP0095 NEW)
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

project(protenixdock LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

option(BDOCK_DOUBLE "Use double precision for parameters (maybe slower, better for debug)." ON)
if(BDOCK_DOUBLE)
  add_definitions(-DENABLE_DOUBLE_PRECISION=1)
endif()

option(BDOCK_PYSDK "If set, Build a Python module. Otherwise, build command-line binaries." OFF)
if(BDOCK_PYSDK)
  set(BDOCK_VETOS OFF)
else()
  set(BDOCK_VETOS ON)
  add_definitions(-DENABLE_VETOS_ACCESS=1)
endif()

# Add a custom command that produces version.cpp, plus a dummy output that's not actually
# produced, in order to force version.cmake to always be re-run before the build.
add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/_version.cpp
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/version.cpp           
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/version.cmake)
add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/_version.cpp)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/_version.cpp
                            PROPERTIES SYMBOLIC TRUE)

find_package(Python3 REQUIRED COMPONENTS Development.Module)
find_package(Boost 1.72.0 REQUIRED COMPONENTS log program_options python serialization)
include_directories(${Boost_INCLUDE_DIRS})

if(BDOCK_VETOS)
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
  find_package(VeTos)
  if(NOT VeTos_FOUND)
    message("VeTos not FOUND, so git clone it and build.")
    set(TOS_INSTALL_DIR "${CMAKE_BINARY_DIR}/ve-tos")
    include(ExternalProject)
    ExternalProject_Add(
      VeTos
      URL https://github.com/volcengine/ve-tos-cpp-sdk/archive/refs/tags/2.6.12.tar.gz
      URL_HASH MD5=a679a1e2b0f7fbf9f10fc2a6222edc0f
      PREFIX         ${CMAKE_BINARY_DIR}/external
      CMAKE_ARGS     -DCMAKE_INSTALL_PREFIX=${TOS_INSTALL_DIR} -DBUILD_SHARED_LIB=ON
    )
    set(VeTos_INCLUDE_DIRS ${TOS_INSTALL_DIR}/include)
    set(VeTos_LIBRARIES ${TOS_INSTALL_DIR}/lib/libve-tos-cpp-sdk-lib.so)
  endif()
endif()

include_directories("${CMAKE_SOURCE_DIR}/include")
set(Dock_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")
include_directories(${Dock_SOURCE_DIR})

file(GLOB BdockOpt_SOURCES "${Dock_SOURCE_DIR}/core/*.cpp"
                           "${Dock_SOURCE_DIR}/ext/*.cpp"
                           "${Dock_SOURCE_DIR}/lib/*.cpp")
list(APPEND BdockOpt_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

add_library(Dock_LIBRARY SHARED ${BdockOpt_SOURCES})
if(BDOCK_VETOS AND NOT VeTos_FOUND)
  add_dependencies(Dock_LIBRARY VeTos)
endif()
set_target_properties(Dock_LIBRARY PROPERTIES OUTPUT_NAME "dock")
target_link_libraries(Dock_LIBRARY ${Boost_LOG_LIBRARY}
                                   ${Boost_SERIALIZATION_LIBRARY})
if(BDOCK_VETOS)
  target_include_directories(Dock_LIBRARY PRIVATE ${VeTos_INCLUDE_DIRS})
  target_link_libraries(Dock_LIBRARY ${VeTos_LIBRARIES})
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # Suppress warnings from `TosClient.h`
  target_compile_options(Dock_LIBRARY PRIVATE "-Wno-deprecated-declarations")
endif()

if(NOT BDOCK_PYSDK)
add_executable(bdock ${Dock_SOURCE_DIR}/cli/bdock.cpp)
target_link_libraries(bdock Dock_LIBRARY ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(bdopt ${Dock_SOURCE_DIR}/cli/bdopt.cpp)
target_link_libraries(bdopt Dock_LIBRARY ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(bdgen ${Dock_SOURCE_DIR}/cli/bdgen.cpp)
target_link_libraries(bdgen Dock_LIBRARY ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(bdscr ${Dock_SOURCE_DIR}/cli/bdscr.cpp)
target_link_libraries(bdscr Dock_LIBRARY ${Boost_PROGRAM_OPTIONS_LIBRARY})
endif()

if(BDOCK_PYSDK)
  add_library(_protenix_dock MODULE ${Dock_SOURCE_DIR}/pyb/api.cpp)
  target_include_directories(_protenix_dock PRIVATE ${Python3_INCLUDE_DIRS})
  target_link_libraries(_protenix_dock Dock_LIBRARY
                                      Python3::Module ${Boost_PYTHON_LIBRARY})
  set_target_properties(_protenix_dock PROPERTIES PREFIX "")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Suppress warnings from `TosClient.h`
    target_compile_options(_protenix_dock PRIVATE "-Wno-deprecated-declarations")
  endif()
endif()

# An OS-dependent trick: https://stackoverflow.com/q/53428219/28276974
if(BDOCK_VETOS)
  set_target_properties(Dock_LIBRARY PROPERTIES INSTALL_RPATH "\${ORIGIN}")
endif()
if(BDOCK_PYSDK)
  set_target_properties(_protenix_dock PROPERTIES INSTALL_RPATH "\${ORIGIN}")
else()
  set_target_properties(bdock bdopt bdgen bdscr PROPERTIES
                        INSTALL_RPATH "\${ORIGIN}/../lib")
endif()

include(GNUInstallDirs)
install(TARGETS Dock_LIBRARY LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
if(BDOCK_PYSDK)
  install(TARGETS Dock_LIBRARY _protenix_dock
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
else()
  install(TARGETS bdock bdopt bdgen bdscr
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
if(BDOCK_VETOS)
  install(FILES ${VeTos_LIBRARIES} DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

option(BDOCK_UTESTS "Build unit tests." OFF)
if(BDOCK_UTESTS)
  if(NOT BDOCK_DOUBLE)
    message(FATAL_ERROR "To build unit tests, please enable option BDOCK_DOUBLE!")
  endif()
  if(BDOCK_PYSDK)
    message(FATAL_ERROR "To build unit tests, please disable option BDOCK_PYSDK!")
  endif()

  include(FetchContent)
  FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/refs/tags/release-1.12.0.tar.gz
      URL_HASH MD5=c483d6c24846cea536a9c2b44085bc14
  )
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  set(Dock_UnitTest_SOURCE_DIR "${CMAKE_SOURCE_DIR}/utest")
  file(GLOB Dock_UnitTest_SOURCES "${Dock_UnitTest_SOURCE_DIR}/*.cpp")

  enable_testing()
  add_executable(bdock_test ${Dock_UnitTest_SOURCES})
  target_link_libraries(bdock_test GTest::gtest_main Dock_LIBRARY)

  include(GoogleTest)
  gtest_discover_tests(bdock_test)
endif()
