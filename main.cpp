//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)

#include "Mathter/Vector.hpp"
#include "Mathter/Matrix.hpp"
#include "Measure.hpp"

#include <iostream>
#include <chrono>
#include <vector>
#include <random>
#include <iomanip>

#ifdef _MSC_VER
#include <conio.h>
#else
#include <iostream>
char _getch() { return std::cin.get(); }
#endif


#include <gtest/gtest.h>



using namespace std;
using namespace mathter;


struct PlainMat4 {
	__m128 stripes[4];

	PlainMat4() = default;
	PlainMat4(float _11, float _21, float _31, float _41,
			  float _12, float _22, float _32, float _42,
			  float _13, float _23, float _33, float _43,
			  float _14, float _24, float _34, float _44)
	{
		(*this)(0, 0) = _11; (*this)(1, 0) = _21; (*this)(2, 0) = _31; (*this)(3, 0) = _41;
		(*this)(0, 1) = _12; (*this)(1, 1) = _22; (*this)(2, 1) = _32; (*this)(3, 1) = _42;
		(*this)(0, 2) = _13; (*this)(1, 2) = _23; (*this)(2, 2) = _33; (*this)(3, 2) = _43;
		(*this)(0, 3) = _14; (*this)(1, 3) = _24; (*this)(2, 3) = _34; (*this)(3, 3) = _44;
	}

	float& operator()(int x, int y) {
		return reinterpret_cast<float*>(&stripes[y])[x];
	}
	float operator()(int x, int y) const {
		return reinterpret_cast<const float*>(&stripes[y])[x];
	}
};

PlainMat4 operator*(const PlainMat4& lhs, const PlainMat4& rhs) {
	PlainMat4 res;

	res.stripes[0] = _mm_mul_ps(rhs.stripes[0], _mm_set1_ps(lhs(0, 0)));
	res.stripes[1] = _mm_mul_ps(rhs.stripes[0], _mm_set1_ps(lhs(0, 1)));
	res.stripes[2] = _mm_mul_ps(rhs.stripes[0], _mm_set1_ps(lhs(0, 2)));
	res.stripes[3] = _mm_mul_ps(rhs.stripes[0], _mm_set1_ps(lhs(0, 3)));

	res.stripes[0] = _mm_add_ps(res.stripes[0], _mm_mul_ps(rhs.stripes[1], _mm_set1_ps(lhs(1, 0))));
	res.stripes[1] = _mm_add_ps(res.stripes[1], _mm_mul_ps(rhs.stripes[1], _mm_set1_ps(lhs(1, 1))));
	res.stripes[2] = _mm_add_ps(res.stripes[2], _mm_mul_ps(rhs.stripes[1], _mm_set1_ps(lhs(1, 2))));
	res.stripes[3] = _mm_add_ps(res.stripes[3], _mm_mul_ps(rhs.stripes[1], _mm_set1_ps(lhs(1, 3))));

	res.stripes[0] = _mm_add_ps(res.stripes[0], _mm_mul_ps(rhs.stripes[2], _mm_set1_ps(lhs(2, 0))));
	res.stripes[1] = _mm_add_ps(res.stripes[1], _mm_mul_ps(rhs.stripes[2], _mm_set1_ps(lhs(2, 1))));
	res.stripes[2] = _mm_add_ps(res.stripes[2], _mm_mul_ps(rhs.stripes[2], _mm_set1_ps(lhs(2, 2))));
	res.stripes[3] = _mm_add_ps(res.stripes[3], _mm_mul_ps(rhs.stripes[2], _mm_set1_ps(lhs(2, 3))));

	res.stripes[0] = _mm_add_ps(res.stripes[0], _mm_mul_ps(rhs.stripes[3], _mm_set1_ps(lhs(3, 0))));
	res.stripes[1] = _mm_add_ps(res.stripes[1], _mm_mul_ps(rhs.stripes[3], _mm_set1_ps(lhs(3, 1))));
	res.stripes[2] = _mm_add_ps(res.stripes[2], _mm_mul_ps(rhs.stripes[3], _mm_set1_ps(lhs(3, 2))));
	res.stripes[3] = _mm_add_ps(res.stripes[3], _mm_mul_ps(rhs.stripes[3], _mm_set1_ps(lhs(3, 3))));

	return res;
}

std::ostream& operator<<(std::ostream& os, const PlainMat4& mat) {
	for (int y = 0; y < 4; ++y) {
		os << "{";
		for (int x = 0; x < 4; ++x) {
			os << mat(x, y) << (x == 4 - 1 ? "" : "\t");
		}
		os << "}\n";
	}
	return os;
}


double MatMulSpeedTestPlain() {
	using LeftT = PlainMat4;
	using RightT = PlainMat4;
	using ResultT = decltype(LeftT()*RightT());

	constexpr int iterationCount = 100'000;
	std::vector<LeftT> left(iterationCount);
	std::vector<RightT> right(iterationCount);
	std::vector<ResultT> result(iterationCount);

	std::minstd_rand rne;
	std::uniform_real_distribution<float> rng;

	for (int i = 0; i < iterationCount; ++i) {
		LeftT& l = left[i];
		RightT& r = right[i];

		for (int x = 0; x < 4; ++x) {
			for (int y = 0; y < 4; ++y) {
				l(x, y) = 2;
			}
		}
		for (int x = 0; x < 4; ++x) {
			for (int y = 0; y < 4; ++y) {
				r(x, y) = 2;
			}
		}
	}

	std::chrono::high_resolution_clock::time_point startTime;
	std::chrono::high_resolution_clock::time_point endTime;

	startTime = std::chrono::high_resolution_clock::now();
	for (int j = 0; j < 100; ++j) {
		for (int i = 0; i < iterationCount; ++i) {
			result[i] = left[i] * right[i];
		}
	}
	endTime = std::chrono::high_resolution_clock::now();

	return chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count() * 1e-9;
}


enum class eDiag {
	MAIN,
	OFF,
	NONE,
};

std::ostream& operator<<(std::ostream& os, eDiag diag) {
	os << [diag] {
		switch (diag) {
			case eDiag::MAIN:
				return "main";
			case eDiag::OFF:
				return "off";
			case eDiag::NONE:
				return "no";
		}
	}();
	return os;
}

struct SvdMetrics {
	eDiag s_diag;
	eDiag uau_diag;
	eDiag vav_diag;
	bool uv_same;
	float det;
	float z1, z2, t1, t2;
};

eDiag IsDiag(const Matrix<float, 2, 2>& M) {
	float dmain = abs(M(0, 0)) + abs(M(1, 1));
	float doff = abs(M(1, 0)) + abs(M(0, 1));
	float tolerance = 10;
	if (dmain > tolerance*doff) {
		return eDiag::MAIN;
	}
	else if (doff > tolerance*dmain) {
		return eDiag::OFF;
	}
	else {
		return eDiag::NONE;
	}
};

void SVD2(Matrix<float, 2, 2> A, Matrix<float, 2, 2>& U, Matrix<float, 2, 2>& S, Matrix<float, 2, 2>& V) {
	float s1, c1, s2, c2, d1, d2;
	Svd2x2Helper(A, c1, s1, c2, s2, d1, d2);

	// write out matrices
	U = {
		c1, s1,
		-s1, c1
	};

	V = {
		c2, -s2,
		s2, c2
	};

	S = {
		d1,
		0,
		0,
		d2,
	};
}

void QR2(const Matrix<float, 2, 2>& A, Matrix<float, 2, 2>& R, Matrix<float, 2, 2>& Q) {
	float x, y, z, c, s;
	Rq2x2Helper(A, x, y, z, c, s);
	R = {
		x,y,
		0,z
	};
	Q = {
		c,s,
		-s,c
	};
}


float ErrorSvd(Matrix<float, 2, 2> A, Matrix<float, 2, 2>& U, Matrix<float, 2, 2>& S, Matrix<float, 2, 2>& V) {
	auto Aapprox = U*S*V;

	float norm = A.Norm();
	auto E = A-Aapprox;
	float enorm = E.Norm();
	if (norm == 0) {
		norm = 1;
	}
	return enorm / norm;
}


void TestSvd(int seed = 1000) {
	std::mt19937_64 rne(seed);
	std::uniform_real_distribution<float> rng(-1, 1);

	std::vector<Matrix<float, 2, 2>> matrices(200'000'000);
	for (auto& M : matrices) {
		M = {
			rng(rne), rng(rne),
			rng(rne), rng(rne)
		};
	}

	float maxError = 0;
	int maxErrorIdx = 0;
	double avgError = 0;

	Matrix<float, 2, 2> U, S, V;

	cout << "sdiag\tpred\tdet\tuvsame\tz1\tz2\tt1\tt2\tc1\ts1" << endl;
	for (size_t i = 0; i<matrices.size(); ++i) {
		auto& M = matrices[i];
		SVD2(M, U, S, V);

		float error = ErrorSvd(M, U, S, V);
		if (isnan(error)) {
			cout << M << endl;
		}
		if (error > maxError) {
			maxError = error;
			maxErrorIdx = i;
		}
		avgError += error;
	}
	avgError /= matrices.size();

	cout << "Error: " << endl;
	cout << "   max: " << maxError << endl;
	cout << "   avg: " << avgError << endl;

	cout << "Max err matrix: " << matrices[maxErrorIdx] << endl;
}


int main(int argc, char* argv[]) {
	//cout << setprecision(10);
	//for (int i = 0; i<10; ++i) {
	//	TestSvd(i*1000);
	//}

	//Matrix<float, 2, 2> M = {
	//	-0.8339210749, -0.7153375745, -0.6953243613, 0.8105902672
	//};
	//Matrix<float, 2, 2> U, S, V;
	//SVD2(M, U, S, V);


	//cout << U << endl;
	//cout << S << endl;
	//cout << V << endl;
	//cout << endl;
	//cout << M << endl;
	//cout << U*S*V << endl;


	//_getch();
	//return 0;

	//Vector<float, 4, false>::DumpLayout(cout);
	//Matrix<float, 4, 4>::DumpLayout(cout);

	//srand(clock());
	//Matrix<float, 4, 4> m1 = {
	//	rand(),2,3,4,
	//	5,6,7,8,
	//	9,10,11,12,
	//	13,14,15,16
	//};
	//Matrix<float, 4, 4> m2 = m1*m1;
	//cout << m2;

	//PlainMat4 mp1 = {
	//	m1(0,0),2,3,4,
	//	5,6,7,8,
	//	9,10,11,12,
	//	13,14,15,16
	//};
	//PlainMat4 mp2;
	//mp2 = mp1*mp1;
	//cout << mp2;


	::testing::InitGoogleTest(&argc, argv);
	auto ret = RUN_ALL_TESTS();
	cout << endl;

#ifdef NDEBUG
	cout << "Performance measurements:" << endl << "-------------------------------------------------------" << endl;
	Measure();
	cout << endl;
#endif

	cout << "Press any key to exit...";
	_getch();
	return ret;
}