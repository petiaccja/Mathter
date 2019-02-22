//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)

#include "Mathter/Vector.hpp"
#include "Mathter/Matrix.hpp"
#include "Measure.hpp"

#include <iostream>

#ifdef _MSC_VER
#include <conio.h>
#else
#include <iostream>
char _getch() { return std::cin.get(); }
#endif

#define CATCH_CONFIG_RUNNER
#include <Catch2/catch.hpp>


using namespace std;
using namespace mathter;



int main(int argc, char* argv[]) {
	int ret = Catch::Session().run(argc, argv);

#ifdef NDEBUG
	cout << endl;
	cout << "Performance measurements:" << endl << "-------------------------------------------------------" << endl;
	Measure();
	cout << endl;

	//cout << endl;
	//cout << "Press any key to exit..." << endl;
	//_getch();
#endif
	return ret;
}