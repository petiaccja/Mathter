#include "Measure.hpp"

#include <random>
#include <chrono>
#include <iostream>
#include <string>
#include <numeric>
#include <cstring>

#include "Mathter/Vector.hpp"
#include "Mathter/Matrix.hpp"
#include "Mathter/Quaternion.hpp"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#undef NOMINMAX
#include <Windows.h>
#elif __unix__
#include <pthread.h>
#include <sched.h>
#endif

#ifdef _MSC_VER
#include <intrin.h>
#define USE_RDTSC
#elif __GNUC__
#include <x86intrin.h>
#define USE_RDTSC
#endif


using namespace std;
using namespace mathter;



template <class Type, int Rows1, int Columns1, eMatrixLayout Layout1, int Rows2 = Columns1, int Columns2 = Rows1, eMatrixLayout Layout2 = Layout1, bool Packed = false>
double MeasureMatrixMultiplication(double* cycles = nullptr) {
	using LeftT = Matrix<Type, Rows1, Columns1, eMatrixOrder::FOLLOW_VECTOR, Layout1, Packed>;
	using RightT = Matrix<Type, Rows2, Columns2, eMatrixOrder::FOLLOW_VECTOR, Layout2, Packed>;
	using ResultT = decltype(LeftT()*RightT());

	constexpr int repeatCount = 50;
	constexpr int iterationCount = 100'000;
	std::vector<LeftT> left(iterationCount);
	std::vector<RightT> right(iterationCount);
	std::vector<ResultT> result(iterationCount);


	memset(left.data(), 2, left.size() * sizeof(LeftT));
	memset(right.data(), 3, right.size() * sizeof(RightT));


	std::chrono::high_resolution_clock::time_point startTime;
	std::chrono::high_resolution_clock::time_point endTime;
	int64_t startCycle = 0, endCycle = -1;


	startTime = std::chrono::high_resolution_clock::now();
#ifdef USE_RDTSC
	startCycle = (int64_t)__rdtsc();
#endif
	for (int j = 0; j < repeatCount; ++j) {
		for (int i = 0; i < iterationCount; ++i) {
			result[i] = left[i] * right[i];
		}
	}
#ifdef USE_RDTSC
	endCycle = (int64_t)__rdtsc();
#endif
	endTime = std::chrono::high_resolution_clock::now();


	if (cycles) {
		*cycles = double(endCycle - startCycle) / iterationCount / repeatCount;
	}
	double totalTime = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count() * 1e-9;
	return totalTime / iterationCount / repeatCount;
}



template <class Type, int Rows, int Columns, eMatrixLayout Layout, bool Packed = false>
double MeasureMatrixAddition(double* cycles = nullptr) {
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
	int64_t startCycle = 0, endCycle = -1;


	startTime = std::chrono::high_resolution_clock::now();
#ifdef USE_RDTSC
	startCycle = (int64_t)__rdtsc();
#endif
	for (int j = 0; j < repeatCount; ++j) {
		for (int i = 0; i < iterationCount; ++i) {
			result[i] = left[i] + right[i];
		}
	}
#ifdef USE_RDTSC
	endCycle = (int64_t)__rdtsc();
#endif
	endTime = std::chrono::high_resolution_clock::now();


	if (cycles) {
		*cycles = double(endCycle - startCycle) / iterationCount / repeatCount;
	}
	double totalTime = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count() * 1e-9;
	return totalTime / iterationCount / repeatCount;
}



template <class Type, int VDim, int MDim, eMatrixLayout Layout, bool Packed>
double MeasureVectorMatrixMultiplication(double* cycles) {
	using MatrixT = Matrix<Type, MDim, MDim, eMatrixOrder::FOLLOW_VECTOR, Layout, Packed>;
	using VectorT = Vector<Type, VDim, Packed>;

	constexpr int repeatCount = 50;
	constexpr int iterationCount = 100'000;
	std::vector<VectorT> left(iterationCount);
	std::vector<MatrixT> right(iterationCount);
	std::vector<VectorT> result(iterationCount);


	memset(left.data(), 2, left.size() * sizeof(VectorT));
	memset(right.data(), 3, right.size() * sizeof(MatrixT));


	std::chrono::high_resolution_clock::time_point startTime;
	std::chrono::high_resolution_clock::time_point endTime;
	int64_t startCycle = 0, endCycle = -1;


	startTime = std::chrono::high_resolution_clock::now();
#ifdef USE_RDTSC
	startCycle = (int64_t)__rdtsc();
#endif
	for (int j = 0; j < repeatCount; ++j) {
		for (int i = 0; i < iterationCount; ++i) {
			result[i] = left[i] * right[i];
		}
	}
#ifdef USE_RDTSC
	endCycle = (int64_t)__rdtsc();
#endif
	endTime = std::chrono::high_resolution_clock::now();


	if (cycles) {
		*cycles = double(endCycle - startCycle) / iterationCount / repeatCount;
	}
	double totalTime = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count() * 1e-9;
	return totalTime / iterationCount / repeatCount;
}



template <class Type, bool Packed>
double MeasureVectorQuatMultiplication(double* cycles) {
	using QuatT = Quaternion<Type, Packed>;
	using VectorT = Vector<Type, 3, Packed>;

	constexpr int repeatCount = 50;
	constexpr int iterationCount = 100'000;
	std::vector<VectorT> left(iterationCount);
	std::vector<QuatT> right(iterationCount);
	std::vector<VectorT> result(iterationCount);


	memset(left.data(), 2, left.size() * sizeof(VectorT));
	memset(right.data(), 3, right.size() * sizeof(QuatT));


	std::chrono::high_resolution_clock::time_point startTime;
	std::chrono::high_resolution_clock::time_point endTime;
	int64_t startCycle = 0, endCycle = -1;


	startTime = std::chrono::high_resolution_clock::now();
#ifdef USE_RDTSC
	startCycle = (int64_t)__rdtsc();
#endif
	for (int j = 0; j < repeatCount; ++j) {
		for (int i = 0; i < iterationCount; ++i) {
			result[i] = left[i] * right[i];
		}
	}
#ifdef USE_RDTSC
	endCycle = (int64_t)__rdtsc();
#endif
	endTime = std::chrono::high_resolution_clock::now();


	if (cycles) {
		*cycles = double(endCycle - startCycle) / iterationCount / repeatCount;
	}
	double totalTime = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count() * 1e-9;
	return totalTime / iterationCount / repeatCount;
}


template <class Type, bool Packed>
double MeasureQuatMultiplication(double* cycles) {
	using QuatT = Quaternion<Type, Packed>;

	constexpr int repeatCount = 50;
	constexpr int iterationCount = 100'000;
	std::vector<QuatT> left(iterationCount);
	std::vector<QuatT> right(iterationCount);
	std::vector<QuatT> result(iterationCount);


	memset(left.data(), 2, left.size() * sizeof(QuatT));
	memset(right.data(), 3, right.size() * sizeof(QuatT));


	std::chrono::high_resolution_clock::time_point startTime;
	std::chrono::high_resolution_clock::time_point endTime;
	int64_t startCycle = 0, endCycle = -1;


	startTime = std::chrono::high_resolution_clock::now();
#ifdef USE_RDTSC
	startCycle = (int64_t)__rdtsc();
#endif
	for (int j = 0; j < repeatCount; ++j) {
		for (int i = 0; i < iterationCount; ++i) {
			result[i] = left[i] * right[i];
		}
	}
#ifdef USE_RDTSC
	endCycle = (int64_t)__rdtsc();
#endif
	endTime = std::chrono::high_resolution_clock::now();


	if (cycles) {
		*cycles = double(endCycle - startCycle) / iterationCount / repeatCount;
	}
	double totalTime = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count() * 1e-9;
	return totalTime / iterationCount / repeatCount;
}


template <class Type>
double MeasureSvd2x2(double* cycles) {
	constexpr int repeatCount = 50;
	constexpr int iterationCount = 100'000;

	std::vector<Matrix<Type, 2, 2>> matrices(iterationCount);

	mt19937_64 rne(1000);
	std::uniform_real_distribution<Type> rng(-1, 1);

	for (auto& M : matrices) {
		M = {
			rng(rne), rng(rne),
			rng(rne), rng(rne)
		};
	}


	std::chrono::high_resolution_clock::time_point startTime;
	std::chrono::high_resolution_clock::time_point endTime;
	int64_t startCycle = 0, endCycle = -1;

	startTime = std::chrono::high_resolution_clock::now();
#ifdef USE_RDTSC
	startCycle = (int64_t)__rdtsc();
#endif
	Type c1, s1, c2, s2, d1, d2;
	for (int j = 0; j < repeatCount; ++j) {
		for (auto& M : matrices) {
			impl::Svd2x2Helper(M, c1, s1, c2, s2, d1, d2);
		}
	}
#ifdef USE_RDTSC
	endCycle = (int64_t)__rdtsc();
#endif
	endTime = std::chrono::high_resolution_clock::now();


	if (cycles) {
		*cycles = double(endCycle - startCycle) / iterationCount / repeatCount;
	}
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
	cout << (success ? "Process set to highest priority." : "Failed to set process priority class - times may have jitter.");
	cout << endl;
	success = SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL);
	cout << (success ? "Thread set to highest priority." : "Failed to set thread priority class - times may have jitter.");
	cout << endl;
	success = SetThreadAffinityMask(GetCurrentThread(), 1) != 0;
	cout << (success ? "Thread is limited to 1st CPU core." : "Failed to set thread affinity - cycle count may be incorrect.");
	cout << endl << endl;
#elif __unix__
	cout << "[Initialize]" << endl;
	pthread_t this_thread = pthread_self();
	sched_param params;
	params.sched_priority = sched_get_priority_max(SCHED_FIFO);
	int ret = pthread_setschedparam(this_thread, SCHED_FIFO, &params);
	cout << (ret == 0 ? "Thread set to highest priority." : "Failed to set thread priority class - times may have jitter.");
	cout << endl << endl;
#endif

	cout << "Running test, please wait..." << endl << endl;


	// Multiplication measures
	std::string mulLabels[8] = {
		"type", "dim1", "dim2", "layout1", "layout2", "packed", "time/op", "\tcycles"
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
	std::array<double, mulNumCases> mulCycles;
	std::array<double, mulNumCases> mulTimes = {
		MeasureMatrixMultiplication<float, 2, 2, ROW, 2, 2, ROW, false>(&mulCycles[0]),
		MeasureMatrixMultiplication<float, 3, 3, ROW, 3, 3, ROW, false>(&mulCycles[1]),
		MeasureMatrixMultiplication<float, 4, 4, ROW, 4, 4, ROW, false>(&mulCycles[2]),
		MeasureMatrixMultiplication<float, 4, 4, ROW, 4, 4, COL, false>(&mulCycles[3]),
		MeasureMatrixMultiplication<float, 4, 4, COL, 4, 4, ROW, false>(&mulCycles[4]),
		MeasureMatrixMultiplication<float, 4, 4, COL, 4, 4, COL, false>(&mulCycles[5]),
		MeasureMatrixMultiplication<double, 4, 4, ROW, 4, 4, ROW, false>(&mulCycles[6]),
		MeasureMatrixMultiplication<float, 4, 3, ROW, 3, 4, ROW, false>(&mulCycles[7]),
		MeasureMatrixMultiplication<float, 3, 4, ROW, 4, 3, ROW, false>(&mulCycles[8]),
		MeasureMatrixMultiplication<float, 3, 4, ROW, 4, 3, ROW, true>(&mulCycles[9]),
	};

	cout << "[Matrix multiplication]" << endl;
	for (auto label : mulLabels) {
		cout << label << "\t";
	}
	cout << endl;
	for (int i = 0; i < 8 * 9; ++i) {
		cout << "-";
	}
	cout << endl;

	for (int i = 0; i < mulNumCases; ++i) {
		for (auto config : mulConfigs[i]) {
			cout << config << "\t";
		}
		printf("%-6.02f ns\t", mulTimes[i] * 1e9);
		printf("%.02f\t", mulCycles[i]);
		cout << endl;
	}
	cout << endl;



	// Addition measures
	std::string addLabels[6] = {
		"type", "dim", "layout", "packed", "time/op", "\tcycles"
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
	std::array<double, addNumCases> addCycles;
	std::array<double, addNumCases> addTimes = {
		MeasureMatrixAddition<float, 1, 1, ROW, false>(&addCycles[0]),
		MeasureMatrixAddition<float, 2, 2, ROW, false>(&addCycles[1]),
		MeasureMatrixAddition<float, 3, 3, ROW, false>(&addCycles[2]),
		MeasureMatrixAddition<float, 4, 4, ROW, false>(&addCycles[3]),
		MeasureMatrixAddition<double, 4, 4, ROW, false>(&addCycles[4]),
		MeasureMatrixAddition<float, 4, 3, ROW, false>(&addCycles[5]),
		MeasureMatrixAddition<float, 3, 4, ROW, false>(&addCycles[6]),
		MeasureMatrixAddition<float, 3, 4, ROW, true>(&addCycles[7]),
	};

	cout << "[Matrix addition]" << endl;
	for (auto label : addLabels) {
		cout << label << "\t";
	}
	cout << endl;
	for (int i = 0; i < 8 * 9; ++i) {
		cout << "-";
	}
	cout << endl;

	for (int i = 0; i < addNumCases; ++i) {
		for (auto config : addConfigs[i]) {
			cout << config << "\t";
		}
		printf("%-6.02f ns\t", addTimes[i] * 1e9);
		printf("%.02f\t", addCycles[i]);
		cout << endl;
	}
	cout << endl;


	// Vector-Matrix measures
	std::string vecmatMulLabels[6] = {
		"type", "dim", "layout", "packed", "time/op", "\tcycles"
	};

	constexpr int vecmatMulNumCases = 4;
	std::string vecmatMulConfigs[vecmatMulNumCases][4] = {
		{ "float",	"3*3x3",	"R",	"F" },
		{ "float",	"4*4x4",	"R",	"F" },
		{ "float",	"3*4x4",	"R",	"F" },
		{ "float",	"4*4x4",	"R",	"T" },
	};
	std::array<double, vecmatMulNumCases> vecmatMulCycles;
	std::array<double, vecmatMulNumCases> vecmatMulTimes = {
		MeasureVectorMatrixMultiplication<float, 3, 3, ROW, false>(&vecmatMulCycles[0]),
		MeasureVectorMatrixMultiplication<float, 4, 4, ROW, false>(&vecmatMulCycles[1]),
		MeasureVectorMatrixMultiplication<float, 3, 4, ROW, false>(&vecmatMulCycles[2]),
		MeasureVectorMatrixMultiplication<float, 4, 4, ROW, true>(&vecmatMulCycles[3]),
	};

	cout << "[Vector-matrix multiplication]" << endl;
	for (auto label : vecmatMulLabels) {
		cout << label << "\t";
	}
	cout << endl;
	for (int i = 0; i < 8 * 9; ++i) {
		cout << "-";
	}
	cout << endl;

	for (int i = 0; i < vecmatMulNumCases; ++i) {
		for (auto config : vecmatMulConfigs[i]) {
			cout << config << "\t";
		}
		printf("%-6.02f ns\t", vecmatMulTimes[i] * 1e9);
		printf("%.02f\t", vecmatMulCycles[i]);
		cout << endl;
	}
	cout << endl;


	// Vector-Quaternion measures
	std::string vecquatMulLabels[6] = {
		"type", "packed", "time/op", "\tcycles"
	};

	constexpr int vecquatMulNumCases = 3;
	std::string vecquatMulConfigs[vecquatMulNumCases][2] = {
		{ "float",	"F" },
		{ "double",	"F" },
		{ "float",	"T" },
	};
	std::array<double, vecquatMulNumCases> vecquatMulCycles;
	std::array<double, vecquatMulNumCases> vecquatMulTimes = {
		MeasureVectorQuatMultiplication<float, false>(&vecquatMulCycles[0]),
		MeasureVectorQuatMultiplication<double, false>(&vecquatMulCycles[1]),
		MeasureVectorQuatMultiplication<float, true>(&vecquatMulCycles[2]),
	};

	cout << "[Vector-quaternion multiplication]" << endl;
	for (auto label : vecquatMulLabels) {
		cout << label << "\t";
	}
	cout << endl;
	for (int i = 0; i < 8 * 9; ++i) {
		cout << "-";
	}
	cout << endl;

	for (int i = 0; i < vecquatMulNumCases; ++i) {
		for (auto config : vecquatMulConfigs[i]) {
			cout << config << "\t";
		}
		printf("%-6.02f ns\t", vecquatMulTimes[i] * 1e9);
		printf("%.02f\t", vecquatMulCycles[i]);
		cout << endl;
	}
	cout << endl;


	// Quaternion-Quaternion measures
	std::string quatMulLabels[6] = {
		"type", "packed", "time/op", "\tcycles"
	};

	constexpr int quatMulNumCases = 3;
	std::string quatMulConfigs[quatMulNumCases][2] = {
		{ "float",	"F" },
		{ "double",	"F" },
		{ "float",	"T" },
	};
	std::array<double, quatMulNumCases> quatMulCycles;
	std::array<double, quatMulNumCases> quatMulTimes = {
		MeasureQuatMultiplication<float, false>(&quatMulCycles[0]),
		MeasureQuatMultiplication<double, false>(&quatMulCycles[1]),
		MeasureQuatMultiplication<float, true>(&quatMulCycles[2]),
	};

	cout << "[Quaternion product]" << endl;
	for (auto label : quatMulLabels) {
		cout << label << "\t";
	}
	cout << endl;
	for (int i = 0; i < 8 * 9; ++i) {
		cout << "-";
	}
	cout << endl;

	for (int i = 0; i < quatMulNumCases; ++i) {
		for (auto config : quatMulConfigs[i]) {
			cout << config << "\t";
		}
		printf("%-6.02f ns\t", quatMulTimes[i] * 1e9);
		printf("%.02f\t", quatMulCycles[i]);
		cout << endl;
	}
	cout << endl;


	// SVD 2x2 measures
	double svdCycles, svdCyclesDouble;
	double svdTime, svdTimeDouble;
	svdTime = MeasureSvd2x2<float>(&svdCycles);
	svdTimeDouble = MeasureSvd2x2<double>(&svdCyclesDouble);

	cout << "[SVD 2x2]" << endl;
	cout << "type\t" << "\ttime/op" << "\tcycles" << endl;
	cout << "float\t";
	printf("%-6.02f ns\t", svdTime * 1e9);
	printf("%.02f\t", svdCycles);
	cout << endl << "double\t";
	printf("%-6.02f ns\t", svdTimeDouble * 1e9);
	printf("%.02f\t", svdCyclesDouble);

	cout << endl << endl;


	// Estimate CPU frequency
	double sumTime = std::accumulate(mulTimes.begin(), mulTimes.end(), 0.0);
	double sumCycles = std::accumulate(mulCycles.begin(), mulCycles.end(), 0.0);
	cout << "Estimated CPU clock frequency: " << (sumCycles / sumTime) / 1e9 << " GHz" << endl;

}