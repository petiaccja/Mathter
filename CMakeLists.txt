## Mathter CMake

cmake_minimum_required(VERSION 3.0.0)

# Cpp flags
set(CMAKE_CXX_STANDARD 14)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	# using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -msse2 -lpthread)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	# using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	# using Visual Studio C++
endif()


# Enumerate sources and create project
project(Mathter)

set(sources
	main.cpp
	Measure.cpp
	TestIntersection.cpp
	TestMatrix.cpp
	TestQuaternion.cpp
	TestVector.cpp)
	
set(headers
	Measure.hpp)
	
file(GLOB mathter_headers Mathter/*.hpp)
file(GLOB mathter_module_headers Mathter/MatrixModules/*.hpp)
file(GLOB mathter_swizzle_headers Mathter/Swizzle/*.hpp)

add_executable(Mathter ${sources} ${headers} ${mathter_headers} ${mathter_swizzle_headers})
include_directories(${CMAKE_SOURCE_DIR}/Lib)

# Source groups
# TODO