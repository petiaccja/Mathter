## Mathter CMake

cmake_minimum_required(VERSION 3.1.0)

# Cpp flags
enable_language(CXX)
set(CMAKE_CXX_STANDARD 14)


if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")

elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2 -lpthread")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")

elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /permissive-")
	set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS} /PROFILE")
endif()


# Enumerate sources and create project
project(Mathter)

set(sources
	main.cpp
	Measure.cpp
	TestIntersection.cpp
	TestMatrix.cpp
	TestQuaternion.cpp
	TestVector.cpp
	TestGenerators.hpp)
	
set(headers
	Measure.hpp)
	
file(GLOB mathter_headers Mathter/*.hpp)
set (mathter_headers ${mathter_headers} Mathter/Mathter.natvis)
file(GLOB mathter_module_headers Mathter/MatrixModules/*.hpp)
file(GLOB mathter_swizzle_headers Mathter/Swizzle/*.hpp)

add_executable(Mathter ${sources} ${headers} ${mathter_headers} ${mathter_module_headers} ${mathter_swizzle_headers})
include_directories(${CMAKE_SOURCE_DIR}/Lib)

# Source groups
source_group("Test" FILES ${sources} ${headers})
source_group("Mathter" FILES ${mathter_headers})
source_group("Mathter\\Modules" FILES ${mathter_module_headers})
source_group("Mathter\\Swizzle" FILES ${mathter_swizzle_headers})

