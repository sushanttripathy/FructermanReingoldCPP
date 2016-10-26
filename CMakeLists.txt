cmake_minimum_required(VERSION 3.6)
project(FR)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp fructermanreingold.cpp fructermanreingold.h csv.h liboptions.cpp liboptions.hpp trim.hpp)
add_executable(FR ${SOURCE_FILES})
link_directories(/usr/local/lib /usr/lib)
target_link_libraries(FR pthread)