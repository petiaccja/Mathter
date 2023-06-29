cmake_minimum_required(VERSION 3.24.0)

## Mathter CMake
project(Mathter)

# Cpp flags
enable_language(CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
	add_compile_options("/MP")
	set(RELWITHDEBINFO_OPTIONS "/Ob2")
	set(RELWITHDEBINFO_REMOVE "/Ob1")
	foreach(option ${RELWITHDEBINFO_REMOVE})
		string(REPLACE ${option} "" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
	endforeach()
	add_compile_options("$<$<CONFIG:RELWITHDEBINFO>:${RELEASE_OPTIONS}>")
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

# Subdirectories
add_subdirectory(include/Mathter)
add_subdirectory(test)

