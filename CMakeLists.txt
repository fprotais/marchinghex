cmake_minimum_required(VERSION 2.8)
project(marchinghex)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



add_subdirectory(ultimaille)
include_directories(ultimaille ultimaille/ext)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if (WIN32) 
	set(BUILD_SHARED_LIBS OFF)
endif()

add_library(libmarchinghex src/intersections.cpp src/marchingcube.cpp src/projection.cpp src/locale_iterative_smoother.cpp src/intersections.h src/marchingcube.h src/projection.h src/locale_iterative_smoother.h)
target_link_libraries(libmarchinghex ultimaille)

add_executable(marchinghex_bimaterial src/marchinghex_bimaterial.cpp)
target_link_libraries(marchinghex_bimaterial ultimaille libmarchinghex)

add_executable(marchinghex_hexmeshing src/marchinghex_hexmeshing.cpp)
target_link_libraries(marchinghex_hexmeshing ultimaille libmarchinghex)

add_executable(gridgenerator src/gridgenerator.cpp)
target_link_libraries(gridgenerator ultimaille libmarchinghex)

add_executable(make_examples src/mhexample.cpp)
target_link_libraries(make_examples ultimaille libmarchinghex)


if(OpenMP_CXX_FOUND)
    target_link_libraries(libmarchinghex OpenMP::OpenMP_CXX)
	target_link_libraries(marchinghex_bimaterial OpenMP::OpenMP_CXX)
    target_link_libraries(marchinghex_hexmeshing OpenMP::OpenMP_CXX)
    target_link_libraries(gridgenerator OpenMP::OpenMP_CXX)
    target_link_libraries(make_examples OpenMP::OpenMP_CXX)
endif()

