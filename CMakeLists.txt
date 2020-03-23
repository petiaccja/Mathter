## Mathter CMake
project(Mathter)

cmake_minimum_required(VERSION 3.10.0)

# Cpp flags
enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)

option(ENABLE_LLVM_COV "Adds compiler flags to generate LLVM source-based code coverage. Only works with Clang." OFF)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
	if (ENABLE_LLVM_COV)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping -mllvm -enable-name-compression=false")
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

set(INIT_MODE "DEFAULT" CACHE STRING "Set default initialization of Mathter types.")
set_property(CACHE INIT_MODE PROPERTY STRINGS DEFAULT NULL INVALID UNINITIALIZED)

if (INIT_MODE STREQUAL "NULL")
	add_definitions("-DMATHTER_NULL_INITIALIZE=1")
elseif(INIT_MODE STREQUAL "INVALID")
	add_definitions("-DMATHTER_INVALID_INITIALIZE=1")
elseif(INIT_MODE STREQUAL "UNINITIALIZED")
	add_definitions("-DMATHTER_DONT_INITIALIZE=1")
endif()

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Bin)

# Include paths
include_directories(${CMAKE_SOURCE_DIR}/Lib)

# Subdirectories
add_subdirectory(Test)

