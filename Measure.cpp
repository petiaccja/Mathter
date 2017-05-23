#include "Measure.hpp"

#include <random>
#include <chrono>
#include <iostream>
#include <string>

#include "Mathter/Vector.hpp"
#include "Mathter/Matrix.hpp"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h>
#endif


using namespace std;
using namespace mathter;


template <class Type, int Rows1, int Columns1, eMatrixLayout Layout1, int Rows2 = Columns1, int Columns2 = Rows1, eMatrixLayout Layout2 = Layout1, bool Packed = false>
double MeasureMatrixMultiplication() {
	using LeftT = Matrix<Type, Rows1, Columns1, eMatrixOrder::FOLLOW_VECTOR, Layout1, Packed>;
	using RightT = Matrix<Type, Rows2, Columns2, eMatrixOrder::FOLLOW_VECTOR, Layout2, Packed>;
	using ResultT = typename decltype(LeftT()*RightT());

	constexpr int repeatCount = 50;
	constexpr int iterationCount = 100'000;
	std::vector<LeftT> left(iterationCount);
	std::vector<RightT> right(iterationCount);
	std::vector<ResultT> result(iterationCount);


	memset(left.data(), 2, left.size() * sizeof(LeftT));
	memset(right.data(), 3, right.size() * sizeof(RightT));


	std::chrono::high_resolution_clock::time_point startTime;
	std::chrono::high_resolution_clock::time_point endTime;


	startTime = std::chrono::high_resolution_clock::now();
	for (int j = 0; j < repeatCount; ++j) {
		for (int i = 0; i < iterationCount; ++i) {
			result[i] = left[i] * right[i];
		}
	}
	endTime = std::chrono::high_resolution_clock::now();


	double totalTime = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count() * 1e-9;
	return totalTime / iterationCount / repeatCount;
}


template <class Type, int Rows, int Columns, eMatrixLayout Layout, bool Packed = false>
double MeasureMatrixAddition() {
	using MatrixT = Matrix<Type, Rows, Columns, eMatrixOrder::FOLLOW_VECTOR, Layout, Packed>;

	constexpr int repeatCount = 50;
	constexpr int iterationCount = 100'000;
	std::vector<MatrixT> left(iterationCount);
	std::vector<MatrixT> right(iterationCount);
	std::vector<MatrixT> result(iterationCount);


	memset(left.data(), 2, left.size() * sizeof(MatrixT));
	memset(right.data(), 3, right.size() * sizeof(MatrixT));


	std::chrono::high_resolution_clock::time_point startTime;
	std::chrono::high_resolution_clock::time_point endTime;


	startTime = std::chrono::high_resolution_clock::now();
	for (int j = 0; j < repeatCount; ++j) {
		for (int i = 0; i < iterationCount; ++i) {
			result[i] = left[i] + right[i];
		}
	}
	endTime = std::chrono::high_resolution_clock::now();


	double totalTime = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count() * 1e-9;
	return totalTime / iterationCount / repeatCount;
}



void Measure() {
	constexpr auto ROW = eMatrixLayout::ROW_MAJOR;
	constexpr auto COL = eMatrixLayout::COLUMN_MAJOR;

#ifdef _WIN32
	BOOL success;
	cout << "[Initialize]" << endl;
	success = SetPriorityClass(GetCurrentProcess(), REALTIME_PRIORITY_CLASS);
	cout << (success ? "Process set to highest priority." : "Failed to set process priority class.");
	cout << endl;
	success = SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL);
	cout << (success ? "Thread set to highest priority." : "Failed to set thread priority class.");
	cout << endl << endl;
#endif

	cout << "Running test, please wait..." << endl << endl;


	// Multiplication measures
	std::string mulLabels[7] = {
		"type", "dim1", "dim2", "layout1", "layout2", "packed", "time/op"
	};

	constexpr int mulNumCases = 10;
	std::string mulConfigs[mulNumCases][6] = {
		{ "float",	"2x2",	"2x2",	"R",	"R",	"F" },
		{ "float",	"3x3",	"3x3",	"R",	"R",	"F" },
		{ "float",	"4x4",	"4x4",	"R",	"R",	"F" },
		{ "float",	"4x4",	"4x4",	"R",	"C",	"F" },
		{ "float",	"4x4",	"4x4",	"C",	"R",	"F" },
		{ "float",	"4x4",	"4x4",	"C",	"C",	"F" },
		{ "double",	"4x4",	"4x4",	"R",	"R",	"F" },
		{ "float",	"4x3",	"3x4",	"R",	"R",	"F" },
		{ "float",	"3x4",	"4x3",	"R",	"R",	"F" },
		{ "float",	"3x4",	"4x3",	"R",	"R",	"T" },
	};
	double mulTimes[mulNumCases] = {
		MeasureMatrixMultiplication<float, 2, 2, ROW, 2, 2, ROW, false>(),
		MeasureMatrixMultiplication<float, 3, 3, ROW, 3, 3, ROW, false>(),
		MeasureMatrixMultiplication<float, 4, 4, ROW, 4, 4, ROW, false>(),
		MeasureMatrixMultiplication<float, 4, 4, ROW, 4, 4, COL, false>(),
		MeasureMatrixMultiplication<float, 4, 4, COL, 4, 4, ROW, false>(),
		MeasureMatrixMultiplication<float, 4, 4, COL, 4, 4, COL, false>(),
		MeasureMatrixMultiplication<double, 4, 4, ROW, 4, 4, ROW, false>(),
		MeasureMatrixMultiplication<float, 4, 3, ROW, 3, 4, ROW, false>(),
		MeasureMatrixMultiplication<float, 3, 4, ROW, 4, 3, ROW, false>(),
		MeasureMatrixMultiplication<float, 3, 4, ROW, 4, 3, ROW, true>(),
	};

	cout << "[Matrix multiplication]" << endl;
	for (auto label : mulLabels) {
		cout << label << "\t";
	}
	cout << endl;
	for (int i = 0; i < 8 * 7; ++i) {
		cout << "-";
	}
	cout << endl;

	for (int i = 0; i < mulNumCases; ++i) {
		for (auto config : mulConfigs[i]) {
			cout << config << "\t";
		}
		printf("%.02f ns", mulTimes[i] * 1e9);
		cout << endl;
	}
	cout << endl;



	// Addition measures
	std::string addLabels[5] = {
		"type", "dim", "layout", "packed", "time/op"
	};

	constexpr int addNumCases = 8;
	std::string addConfigs[addNumCases][4] = {
		{ "float",	"1x1",	"R",	"F" },
		{ "float",	"2x2",	"R",	"F" },
		{ "float",	"3x3",	"R",	"F" },
		{ "float",	"4x4",	"R",	"F" },
		{ "double",	"4x4",	"R",	"F" },
		{ "float",	"4x3",	"R",	"F" },
		{ "float",	"3x4",	"R",	"F" },
		{ "float",	"3x4",	"R",	"T" },
	};
	double addTimes[addNumCases] = {
		MeasureMatrixAddition<float, 1, 1, ROW, false>(),
		MeasureMatrixAddition<float, 2, 2, ROW, false>(),
		MeasureMatrixAddition<float, 3, 3, ROW, false>(),
		MeasureMatrixAddition<float, 4, 4, ROW, false>(),
		MeasureMatrixAddition<double, 4, 4, ROW, false>(),
		MeasureMatrixAddition<float, 4, 3, ROW, false>(),
		MeasureMatrixAddition<float, 3, 4, ROW, false>(),
		MeasureMatrixAddition<float, 3, 4, ROW, true>(),
	};

	cout << "[Matrix addition]" << endl;
	for (auto label : addLabels) {
		cout << label << "\t";
	}
	cout << endl;
	for (int i = 0; i < 8 * 7; ++i) {
		cout << "-";
	}
	cout << endl;

	for (int i = 0; i < addNumCases; ++i) {
		for (auto config : addConfigs[i]) {
			cout << config << "\t";
		}
		printf("%.02f ns", addTimes[i] * 1e9);
		cout << endl;
	}
	cout << endl;
}