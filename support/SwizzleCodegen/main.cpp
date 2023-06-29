//L=============================================================================
//L This software is distributed under the MIT license.
//L Copyright 2021 Péter Kardos
//L=============================================================================


#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <cassert>


using namespace std;


std::string GenSwizzleDefinitions(int vdim, int sdim) {
	std::string lines;
	
	// vdim specifies which members we can use:
	// 1: x
	// 2: x, y
	// 3: x, y, z
	// 4: x, y, z, w
	// sdim specifies how many of them to use:
	// 2: ??
	// 3: ???
	// 4: ????

	// get the last number we count until
	// number is written in base vdim, and has sdim biggest digits
	// i.e. vdim=3, we can use x(0),y(1),z(2)
	// sdim = 2, so we have 2 digits
	// maxValue is going to be zz or 22(base 3) = 2*3^1 + 2*3^0
	int maxValue = 0;
	for (int i = 0; i < sdim; ++i) {
		maxValue *= vdim;
		maxValue += vdim - 1;
	}

	std::array<int, 16> digits; // 16 just to make sure, 4 is enough
	std::array<char, 16> letters = {
		'x', 'y', 'z', 'w',
		'A', 'B', 'C', 'D',
		'E', 'F', 'G', 'H',
		'I', 'J', 'K', 'L',
	};
	for (int members = 0; members <= maxValue; ++members) {
		// extract base vdim digits from the number
		int currentMembers = members;
		for (int i = 0; i < sdim; ++i) {
			digits[sdim - i - 1] = currentMembers % vdim;
			currentMembers /= vdim;
		}

		std::stringstream line;
		line << "Swizzle<ST";
		for (int i = 0; i < sdim; ++i) {
			line << ", " << digits[i];
		}
		line << "> ";
		for (int i = 0; i < sdim; ++i) {
			line << letters[digits[i]];
		}
		line << ";" << endl;

		lines += line.str();
	}

	return lines;
}


int main() {
	std::string outputFileDir = GENDIR;


	for (int vdim = 1; vdim <= 4; ++vdim) {
		std::stringstream ss;
		std::ofstream outputFile;
		ss << "\\Swizzle\\Swizzle_" << vdim << ".inc.hpp";
		std::string outputFileName = outputFileDir + ss.str();
		outputFile.open(outputFileName, ios::trunc);
		assert(outputFile.is_open());

		for (int sdim = 2; sdim <= 4; ++sdim) {
			std::string lines = GenSwizzleDefinitions(vdim, sdim);
			outputFile << lines;
		}

		outputFile.close();
	}


	return 0;
}