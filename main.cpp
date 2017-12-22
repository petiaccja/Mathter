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

SvdMetrics SVD2(Matrix<float, 2, 2> A, Matrix<float, 2, 2>& U, Matrix<float, 2, 2>& S, Matrix<float, 2, 2>& V, bool flip = false) {
	float s1, c1, s2, c2, d1, d2;
	//auto internals = Svd2x2Helper(A, c1, s1, c2, s2, d1, d2, flip);
	Svd2x2HelperQR(A, c1, s1, c2, s2, d1, d2);
	//Svd2x2Positive(c1, s1, c2, s2, d1, d2);

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

	auto AAt = A*A.Transposed();
	auto uau = U.Transposed()*AAt*U;

	auto AtA = A.Transposed()*A;
	auto vav = V*AtA*V.Transposed();

	auto uav = U.Transposed()*A*V.Transposed();

	SvdMetrics metrics;
	metrics.s_diag = IsDiag(uav);
	metrics.uau_diag = IsDiag(uau);
	metrics.vav_diag = IsDiag(vav);
	metrics.uv_same = abs(uau(0, 0)) > abs(uau(1, 1)) && abs(vav(0, 0)) > abs(vav(1, 1));
	metrics.det = A(0, 0)*A(1, 1) - A(1, 0)*A(0, 1);
	metrics.z1 = 0;//internals.z1;
	metrics.z2 = 0;//internals.z2;
	metrics.t1 = 0;//internals.t1;
	metrics.t2 = 0;//internals.t2;

	return metrics;
}

void QR2(const Matrix<float, 2, 2>& A, Matrix<float, 2, 2>& R, Matrix<float, 2, 2>& Q) {
	float x, y, z, c, s;
	Qr2x2Helper(A, x, y, z, c, s);
	R = {
		x,y,
		0,z
	};
	Q = {
		c,s,
		-s,c
	};
}

bool Pred(const SvdMetrics& m, const Matrix<float, 2, 2>& M, const Matrix<float, 2, 2>& U, const Matrix<float, 2, 2>& S, const Matrix<float, 2, 2>& V) {
	float a = M(0,0);
	float b = M(0,1);
	float c = M(1,0);
	float d = M(1,1);

	float c1 = -U(0, 0);
	float s1 = U(0, 1);
	float c2 = V(0, 0);
	float s2 = V(0, 1);

	return c1 > c2;
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

	int uvsameButWrong = 0;
	int good = 0;
	int nondiag = 0;
	int uvsame = 0;
	int uvnotsameAndWrong = 0;
	float maxError = 0;
	int maxErrorIdx = 0;
	double avgError = 0;


	Matrix<float, 2, 2> U, S, V;

	cout << "sdiag\tpred\tdet\tuvsame\tz1\tz2\tt1\tt2\tc1\ts1" << endl;
	for (size_t i = 0; i<matrices.size(); ++i) {
		auto& M = matrices[i];
		SvdMetrics metrics = SVD2(M, U, S, V);

		if (false) {
			cout << metrics.s_diag << '\t';
			cout << Pred(metrics, M, U, S, V) << '\t';
			cout << sign(metrics.det) << '\t';
			cout << (metrics.uv_same ? "yes" : "no") << '\t';
			cout << sign(metrics.z1) << '\t';
			cout << sign(metrics.z2) << '\t';
			cout << metrics.t1 << '\t';
			cout << metrics.t2 << '\t';
			cout << U.Determinant() << '\t';
			cout << V.Determinant() << '\t';
			cout << M;
			cout << endl;
		}

		if (metrics.uv_same && metrics.s_diag != eDiag::MAIN) {
			++uvsameButWrong;
		}
		if (metrics.s_diag == eDiag::MAIN) {
			++good;
		}
		else {
			cout << "BAD: " << M << endl;
		}
		//if (metrics.s_diag == eDiag::NONE) {
		//	++nondiag;
		//}
		//if (metrics.uv_same) {
		//	++uvsame;
		//}
		//if (!metrics.uv_same && metrics.s_diag != eDiag::MAIN) {
		//	++uvnotsameAndWrong;
		//}
		//if (metrics.s_diag == eDiag::NONE) {
		//	//cout << M << "\t\t" << U.Transposed()*M*V.Transposed() << endl;
		//}

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

	cout << "Good, bad: " << good << ", " << matrices.size() - good << endl;
	//cout << "UV same: " << uvsame << endl;
	//cout << "UV same but wrong: " << uvsameButWrong << endl;
	//cout << "UV NOT same and wrong: " << uvnotsameAndWrong << endl;
	//cout << "Non diagonal: " << nondiag << endl;
	cout << "Error: " << endl;
	cout << "   max: " << maxError << endl;
	cout << "   avg: " << avgError << endl;

	cout << "Max err matrix: " << matrices[maxErrorIdx] << endl;
}


int main(int argc, char* argv[]) {
	cout << setprecision(10);


	for (int i = 1; i<=10; ++i) {
		TestSvd(i);
	}
	_getch();
	return 0;

	Matrix<float, 2, 2> M, Mc, U, S, V, R, Q;

	M = {
		//0.5224171877, -0.9644564986, 0.9644948244, 0.5225217342
		//0.5225211877, -0.9644944986, 0.9644948244, 0.5225217342
		//1, -1, 0.5, sqrt(2)
		//-1.2e-37, 2.3e-37, 1e-20, 0
		//0.307, -0.589, -0.77, 0.221
		//-0.8877753615, -0.7042343616, 0.00017619133, -5.531311035e-05
		//-0.94581604, -0.1507509947, 0.0005884170532, 0.0005111694336
		-0.7336711884, 0.7314484119, 0.3577359915, 0.306984067
	};
	//M.DecomposeSVD(U, S, V);

	QR2(M, R, Q);
	Mc = R*Q;
	cout << "R = " << R << endl;
	cout << "Q = " << Q << endl << endl;
	cout << "M  = " << M << endl << "M' = " << Mc << endl;
	//_getch();
	//return 0;
	cout << endl;

	SVD2(M, U, S, V);

	cout << "U = " << U << endl << "S = " << S << endl << "V = " << V << endl;
	cout << endl;

	cout << "A  = " << M << endl;
	Mc = U*S*V;
	cout << "A' = " << Mc << endl;
	cout << "err = " << ErrorSvd(M, U, S, V) << endl;

	_getch();
	return 0;


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