# CMAKE BUILD FILE for Spherical Harmonic Gravity Library
cmake_minimum_required(VERSION 3.10)
project(SphericalHarmonicGravity VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Add external libraries
add_subdirectory(extern/APCGrav)

#Main SHG library
include_directories(include)
file(GLOB SHG_SRC "src/*.cpp")
add_library(SHG SHARED src/SHG.cpp)
target_include_directories(SHG
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
		$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Main user applications
add_executable(GravityCalc src/egm2008_gravity_calc.cpp)
target_link_libraries(GravityCalc SHG)

add_executable(PotentialCalc src/egm2008_potential_calc.cpp)
target_link_libraries(PotentialCalc SHG)

add_executable(Demo src/Demo.cpp)
target_link_libraries(Demo SHG)

# Install targets and headers
install(TARGETS SHG
	EXPORT SHGTargets
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Package config
set(SHG_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/SHG)

write_basic_package_version_file(
	${CMAKE_CURRENT_BINARY_DIR}/SHGConfigVersion.cmake
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/SHGConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/SHGConfig.cmake
	INSTALL_DESTINATION ${SHG_INSTALL_CMAKEDIR}
)

install(FILES
	${CMAKE_CURRENT_BINARY_DIR}/SHGConfig.cmake
	${CMAKE_CURRENT_BINARY_DIR}/SHGConfigVersion.cmake
	DESTINATION ${SHG_INSTALL_CMAKEDIR}
)

install(EXPORT SHGTargets
	FILE SHGTargets.cmake
	NAMESPACE SHG::
	DESTINATION ${SHG_INSTALL_CMAKEDIR}
#Share include directories with higher level cmakelists
target_include_directories(SHG PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

