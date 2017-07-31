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


project(Mathter)

set(sources
	main.cpp
	Example.cpp
	Measure.cpp
	TestIntersection.cpp
	TestMatrix.cpp
	TestQuaternion.cpp
	TestVector.cpp)
	
set(headers
	Measure.hpp
	Mathter/Geometry.hpp
	Mathter/Matrix.hpp
	Mathter/Quaternion.hpp
	Mathter/Simd.hpp
	Mathter/Simd_SSE2.hpp
	Mathter/Utility.hpp
	Mathter/Vector.hpp)
	
add_executable(Mathter ${sources} ${headers})
include_directories(${CMAKE_SOURCE_DIR}/Lib)



if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	# using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	target_link_libraries(Mathter optimized pthread "${CMAKE_SOURCE_DIR}/Lib/gcc/release64/libgtest.a")
	target_link_libraries(Mathter debug pthread "${CMAKE_SOURCE_DIR}/Lib/gcc/debug64/libgtest.a")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	# using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	# using Visual Studio C++
endif()
