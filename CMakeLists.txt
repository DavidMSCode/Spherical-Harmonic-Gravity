# CMAKE BUILD FILE for Spherical Harmonic Gravity Library
cmake_minimum_required(VERSION 3.10)
project(SphericalHarmonicGravity)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Main SHG library
include_directories(include)
file(GLOB SHG_SRC "src/*.cpp")
add_library(SHG src/SHG.cpp)

# Main user applications
add_executable(GravityCalc src/egm2008_gravity_calc.cpp)
target_link_libraries(GravityCalc SHG)

add_executable(PotentialCalc src/egm2008_potential_calc.cpp)
target_link_libraries(PotentialCalc SHG)

add_executable(Demo src/Demo.cpp)
target_link_libraries(Demo SHG)

