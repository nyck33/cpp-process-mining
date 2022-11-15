cmake_minimum_required(VERSION 3.22)
project(cpp_process_mining)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(cpp_process_mining
        mim.h
        mimest.cpp
        utils.h)
