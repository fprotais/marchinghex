cmake_minimum_required(VERSION 2.8)
project(marchinghex)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_subdirectory(ultimaille)
include_directories(ultimaille)

#add_library(mcl intersections.cpp marchingcube.cpp)
#set_property(TARGET mcl PROPERTY PUBLIC_HEADER intersections.h marchingcube.h)
#target_link_libraries(mcl ultimaille)

add_executable(marchinghex main.cpp)
target_link_libraries(marchinghex ultimaille)


