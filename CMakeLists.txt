## Mathter CMake

cmake_minimum_required(VERSION 3.1.0)

# Cpp flags
enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)

option(ENABLE_LLVM_COV "Adds compiler flags to generate LLVM source-based code coverage. Only works with Clang." OFF)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
	if (ENABLE_LLVM_COV)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
	endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2 -lpthread")
	if (CMAKE_SYSTEM_NAME MATCHES "Windows")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
	endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")

elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /permissive-")
	set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS} /PROFILE")
endif()


# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Bin)


# Enumerate sources and create project
project(Mathter)

set(sources
	main.cpp
	Measure.cpp
	TestIntersection.cpp
	TestMatrix.cpp
	TestQuaternion.cpp
	TestVector.cpp
	TestGenerators.hpp
	TestGeometry.cpp
	Example.cpp)
	
set(headers
	Measure.hpp)
	
file(GLOB mathter_headers Mathter/*.hpp)
set (mathter_headers ${mathter_headers} Mathter/Mathter.natvis)
file(GLOB mathter_module_headers Mathter/MatrixModules/*.hpp)
file(GLOB mathter_swizzle_headers Mathter/Swizzle/*.hpp)

include_directories(${CMAKE_SOURCE_DIR}/Lib)
add_executable(Mathter ${sources} ${headers} ${mathter_headers} ${mathter_module_headers} ${mathter_swizzle_headers})

# Source groups
source_group("Test" FILES ${sources} ${headers})
source_group("Mathter" FILES ${mathter_headers})
source_group("Mathter\\Modules" FILES ${mathter_module_headers})
source_group("Mathter\\Swizzle" FILES ${mathter_swizzle_headers})

