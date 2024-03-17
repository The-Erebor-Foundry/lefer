
cmake_minimum_required(VERSION 3.22)

project(lefer CXX)
add_library(lefer STATIC src/main.cpp)


add_executable(examples1 examples/src/even_spaced_curves.cpp)
target_include_directories(examples1 PUBLIC src)
target_link_libraries(examples1 lefer)
