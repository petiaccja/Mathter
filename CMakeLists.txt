cmake_minimum_required(VERSION 3.24.0)

# Project version
if (NOT DEFINED MATHTER_VERSION)
	set(MATHTER_VERSION 0.0.1)
endif()

# Project
enable_language(CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(Mathter VERSION ${MATHTER_VERSION})

# Project options
option(MATHTER_USE_XSIMD "Uses XSimd for vectorization of math routines. Uses scalar fallback if turned off." ON)
set(MATHTER_INIT_MODE "DEFAULT" CACHE STRING "Set default initialization of Mathter types.")
set_property(CACHE MATHTER_INIT_MODE PROPERTY STRINGS DEFAULT NULL INVALID UNINITIALIZED)
option(MATHTER_BUILD_TESTS "Include or exclude tests from the generated project." ON)
option(MATHTER_BUILD_BENCHMARKS "Include or exclude bechmarks from the generated project." ON)
option(ENABLE_LLVM_COV "Adds compiler flags to generate LLVM source-based code coverage. Only works with Clang." OFF)

# Global compiler flags
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
	if (ENABLE_LLVM_COV)
		message("Using source-based coverage")
		add_compile_options("-fprofile-instr-generate" "-fcoverage-mapping" "-mllvm" "-enable-name-compression=false")
		add_link_options("-fprofile-instr-generate" "-fcoverage-mapping")
	endif()
endif()

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Subdirectories
add_subdirectory(include/Mathter)
if (${MATHTER_BUILD_TESTS})
	add_subdirectory(test)
endif()
if (${MATHTER_BUILD_BENCHMARKS})
	add_subdirectory(benchmark)
endif()

# Installation
install(TARGETS Mathter EXPORT MathterTargets
	FILE_SET headers DESTINATION "include"
	FILE_SET natvis DESTINATION "include"
)

install(EXPORT MathterTargets
        FILE MathterTargets.cmake
        NAMESPACE Mathter::
        DESTINATION lib/cmake/${PROJECT_NAME})

include(CMakePackageConfigHelpers)

set(INCLUDE_INSTALL_DIR "include" CACHE PATH "Location of header files" )

configure_package_config_file(
	${PROJECT_NAME}Config.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
	INSTALL_DESTINATION lib/cmake/${PROJECT_NAME}
	PATH_VARS INCLUDE_INSTALL_DIR
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  COMPATIBILITY SameMajorVersion
)

install(
	FILES
		${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
		${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
	DESTINATION lib/cmake/${PROJECT_NAME}
)