cmake_minimum_required(VERSION 2.8.4)
project(kalman)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(EIGEN_INC_DIR /Users/emanuele/eigen)
include_directories(${EIGEN_INC_DIR})

set(SOURCE_FILES kalman.cc kalman-test2D.cpp)
add_executable(kalman-test2D ${SOURCE_FILES})
