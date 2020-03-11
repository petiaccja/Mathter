//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)

#include <iostream>

#define CATCH_CONFIG_RUNNER
#include <Catch2/catch.hpp>

using namespace std;

int main(int argc, char* argv[]) {
	int ret = Catch::Session().run(argc, argv);
	return ret;
}