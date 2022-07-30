add_library(Mathter INTERFACE)

message("${CMAKE_CURRENT_SOURCE_DIR}/..")
target_include_directories(Mathter INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/..")

target_sources(Mathter
	INTERFACE
		"Mathter.natvis"
)