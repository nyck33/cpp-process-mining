cmake_minimum_required(VERSION 3.22)
project(cpp_process_mining)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(cpp_process_mining
        mim.h
        mimest.cpp
        utils.h)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


