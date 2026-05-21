# CMAKE BUILD FILE for Spherical Harmonic Gravity Library
cmake_minimum_required(VERSION 3.10)
project(SphericalHarmonicGravity VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

enable_testing()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Geopotential artifact bins (shared with Geopotential.jl)
option(SHG_DOWNLOAD_GEOPOTENTIAL_ARTIFACTS "Download geopotential artifact bins during build" ON)
option(SHG_INSTALL_GEOPOTENTIAL_ARTIFACTS "Install geopotential artifact bins with SHG" ON)
set(SHG_GEOPOTENTIAL_ARTIFACT_URL "https://github.com/DavidMSCode/geopotential-data/releases/download/v1.0.0/geopotential_bins.tar.gz" CACHE STRING "URL for the geopotential artifact tarball")
set(SHG_GEOPOTENTIAL_ARTIFACT_ARCHIVE "" CACHE FILEPATH "Optional local geopotential artifact tarball to use instead of downloading")
set(SHG_GEOPOTENTIAL_ARTIFACTS_DIR "${CMAKE_BINARY_DIR}/geopotential_bins" CACHE PATH "Directory where geopotential artifact bins are stored")
set(SHG_GEOPOTENTIAL_ARTIFACTS_BIN_DIR "${SHG_GEOPOTENTIAL_ARTIFACTS_DIR}/bin" CACHE PATH "Directory containing extracted geopotential model bins")
set(SHG_GEOPOTENTIAL_ARTIFACTS_STAMP "${SHG_GEOPOTENTIAL_ARTIFACTS_DIR}/.download_complete")
set(SHG_GEOPOTENTIAL_INSTALL_ARTIFACTS_DIR "${CMAKE_INSTALL_DATADIR}/SHG/geopotential_bins")

add_custom_command(
	OUTPUT "${SHG_GEOPOTENTIAL_ARTIFACTS_STAMP}"
	COMMAND ${CMAKE_COMMAND}
			-DSHG_GEOPOTENTIAL_ARTIFACT_URL=${SHG_GEOPOTENTIAL_ARTIFACT_URL}
			-DSHG_GEOPOTENTIAL_ARTIFACT_ARCHIVE=${SHG_GEOPOTENTIAL_ARTIFACT_ARCHIVE}
			-DSHG_GEOPOTENTIAL_ARTIFACTS_DIR=${SHG_GEOPOTENTIAL_ARTIFACTS_DIR}
			-P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/DownloadGeopotentialArtifacts.cmake
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/DownloadGeopotentialArtifacts.cmake
	COMMENT "Preparing geopotential artifact bins"
	VERBATIM
)

add_custom_target(download_geopotential_artifacts DEPENDS "${SHG_GEOPOTENTIAL_ARTIFACTS_STAMP}")


#Main SHG library
include_directories(include)
file(GLOB SHG_SRC "src/*.cpp")
add_library(SHG SHARED src/SHG.cpp)
target_include_directories(SHG
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
		$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if(SHG_DOWNLOAD_GEOPOTENTIAL_ARTIFACTS)
	add_dependencies(SHG download_geopotential_artifacts)
endif()

# Main user applications
add_executable(GravityCalc src/egm2008_gravity_calc.cpp)
target_link_libraries(GravityCalc SHG)

add_executable(PotentialCalc src/egm2008_potential_calc.cpp)
target_link_libraries(PotentialCalc SHG)

add_executable(Demo src/Demo.cpp)
target_link_libraries(Demo SHG)

# Profiling utility
add_executable(EGM2008ProfileRunner tests/egm2008_profile_runner.cpp)
target_link_libraries(EGM2008ProfileRunner SHG)
add_executable(EGM2008ProfileAccel tests/egm2008_profile_accel.cpp)
target_link_libraries(EGM2008ProfileAccel SHG)
add_executable(EGM2008ProfilePotential tests/egm2008_profile_potential.cpp)
target_link_libraries(EGM2008ProfilePotential SHG)

# Regression tests (CTest)
add_executable(EGM2008RegressionTest tests/egm2008_regression_ctest.cpp)
target_link_libraries(EGM2008RegressionTest SHG)
add_executable(EGM2008BaselineGenerator tests/egm2008_generate_baseline.cpp)
target_link_libraries(EGM2008BaselineGenerator SHG)
add_executable(EGM2008WorkspaceToggleTest tests/egm2008_workspace_toggle_ctest.cpp)
target_link_libraries(EGM2008WorkspaceToggleTest SHG)
add_executable(EGM2008ThreadSafetyTest tests/egm2008_thread_safety_ctest.cpp)
target_link_libraries(EGM2008ThreadSafetyTest SHG)
add_executable(EGM2008ArtifactParityTest tests/egm2008_artifact_parity_ctest.cpp)
target_link_libraries(EGM2008ArtifactParityTest SHG)
add_executable(EGM2008MethodBenchmark tests/egm2008_method_benchmark.cpp)
target_link_libraries(EGM2008MethodBenchmark SHG)
if(SHG_DOWNLOAD_GEOPOTENTIAL_ARTIFACTS)
	add_dependencies(EGM2008ArtifactParityTest download_geopotential_artifacts)
	add_dependencies(EGM2008MethodBenchmark download_geopotential_artifacts)
endif()
add_test(
	NAME EGM2008MethodBenchmark
	COMMAND EGM2008MethodBenchmark 720 25 8 --human
)
set_tests_properties(EGM2008MethodBenchmark PROPERTIES WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(
	NAME EGM2008MethodBenchmarkCSV
	COMMAND EGM2008MethodBenchmark 720 25 8 --csv ${CMAKE_BINARY_DIR}/EGM2008MethodBenchmark.csv
)
set_tests_properties(EGM2008MethodBenchmarkCSV PROPERTIES WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(
	NAME EGM2008Regression
	COMMAND EGM2008RegressionTest ${CMAKE_SOURCE_DIR}/tests/egm2008_regression_baseline.csv
)
set_tests_properties(EGM2008Regression PROPERTIES WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(
	NAME EGM2008WorkspaceToggle
	COMMAND EGM2008WorkspaceToggleTest
)
set_tests_properties(EGM2008WorkspaceToggle PROPERTIES WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(
	NAME EGM2008ThreadSafety
	COMMAND EGM2008ThreadSafetyTest
)
set_tests_properties(EGM2008ThreadSafety PROPERTIES WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(
	NAME EGM2008ArtifactParity
	COMMAND EGM2008ArtifactParityTest
)
set_tests_properties(EGM2008ArtifactParity PROPERTIES WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

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

# Install geopotential artifact bins if they are available in the build tree.
# This keeps the repo source-only and copies only the downloaded build cache.
if(SHG_INSTALL_GEOPOTENTIAL_ARTIFACTS AND EXISTS "${SHG_GEOPOTENTIAL_ARTIFACTS_BIN_DIR}/models.json")
	install(DIRECTORY ${SHG_GEOPOTENTIAL_ARTIFACTS_BIN_DIR}/
		DESTINATION ${SHG_GEOPOTENTIAL_INSTALL_ARTIFACTS_DIR}/bin
	)
endif()

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
	PATH_VARS SHG_GEOPOTENTIAL_INSTALL_ARTIFACTS_DIR
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
)

#Share include directories with higher level cmakelists
target_include_directories(SHG PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

