cmake_minimum_required(VERSION 2.8)
project(stlbfgs)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(RELATIVE_BIN_DIR bin/)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${RELATIVE_BIN_DIR}/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if (MSVC)
    # warning level 4 (and all warnings as errors, /WX)
    add_compile_options(/W4)
else()
    # lots of warnings and all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic)
endif()

file(GLOB SOURCES *.cpp *.h)
add_library(stlbfgs ${SOURCES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(stlbfgs OpenMP::OpenMP_CXX)
endif()
set_target_properties(stlbfgs PROPERTIES SOVERSION 1.0)

option(STLBFGS_UNIT_TESTS "STLBFGS_UNIT_TESTS" OFF)
if (STLBFGS_UNIT_TESTS)
    add_subdirectory(tests)
endif()

add_executable(stlbfgs-helloworld helloworld.cpp)
target_link_libraries(stlbfgs-helloworld ${CMAKE_DL_LIBS} stlbfgs)
IF (NOT WIN32)
  target_link_libraries(stlbfgs-helloworld m)
ENDIF()

