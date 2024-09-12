add_executable(Examples)


target_sources(Examples
	PRIVATE
		"main.cpp"
		"Math.hpp"
		"Camera.hpp"
		"Camera.cpp"
		"Object.hpp"
		"Object.cpp"
)


target_link_libraries(Examples Mathter)