
cmake_minimum_required(VERSION 3.22)

project(lefer CXX)
include_directories(PUBLIC ./../include/)
link_libraries(m)
add_library(lefer src/main.cpp)

