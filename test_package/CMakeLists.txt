cmake_minimum_required(VERSION 3.25)

project(TestMathterPackage)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(TestMathterPackage)

target_sources(TestMathterPackage
	PRIVATE
		main.cpp
)

find_package(Mathter 99.0.0)
target_link_libraries(TestMathterPackage Mathter::Mathter)