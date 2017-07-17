cmake_minimum_required(VERSION 3.0.0)

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
include_directories(${CMAKE_SOURCE_DIR})
