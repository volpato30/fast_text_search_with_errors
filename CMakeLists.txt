cmake_minimum_required(VERSION 3.12)
project(wild_type_match)

set(CMAKE_CXX_STANDARD 14)
include_directories(./include)

add_executable(wild_type_match main.cpp include/hamming1_search.h include/string_split.h)