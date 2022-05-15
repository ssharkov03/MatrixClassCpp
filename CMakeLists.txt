cmake_minimum_required(VERSION 3.20)
project(matrix_cpp)

set(CMAKE_CXX_STANDARD 14)

add_executable(matrix
        matrix.cpp
        matrix.h)

add_executable(hamming
        hamming.cpp
        hamming.h)

add_executable(test
        test.cpp)

