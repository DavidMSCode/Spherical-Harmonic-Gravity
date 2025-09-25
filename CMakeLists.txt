# CMAKE BUILD FILE for Spherical Harmonic Gravity Library
cmake_minimum_required(VERSION 3.10)
project(SphericalHarmonicGravity)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Main SHG library
include_directories(include)
file(GLOB SHG_SRC "src/*.cpp")
add_library(SHG src/SHG.cpp)

#main executable
add_executable(Main src/main.cpp)
target_link_libraries(Main SHG)
