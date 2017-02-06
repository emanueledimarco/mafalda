cmake_minimum_required(VERSION 2.8.4)
project(kalman)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -g3 -Wall -O0")

set(EIGEN_INC_DIR /Users/emanuele/eigen)
include_directories(${EIGEN_INC_DIR})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

#---Create  a main program using the library
set(SOURCE_FILES kalman.cc LinearTrackFinder.cc TracksAnalyzer.cc btf-tracking.cpp)
add_executable(btf-tracking ${SOURCE_FILES})
target_link_libraries(btf-tracking ${ROOT_LIBRARIES})

