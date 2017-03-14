#include "Mathter\Vector.hpp"
#include "Mathter\Matrix.hpp"

#include <iostream>
#include <chrono>
#include <vector>
#include <random>

using namespace std;
using namespace mathter;


template <class T, int Columns, int Rows, eMatrixOrder Order>
std::ostream& operator<<(std::ostream& os, const Matrix<T, Columns, Rows, Order>& mat) {
	for (int y = 0; y < mat.Height(); ++y) {
		os << "{";
		for (int x = 0; x < mat.Width(); ++x) {
			os << mat(x, y) << (x == mat.Width() - 1 ? "" : "\t");
		}
		os << "}\n";
	}
	return os;
}

template <class T, int D>
std::ostream& operator<<(std::ostream& os, const Vector<T, D>& v) {
	os << "{";
	for (int x = 0; x < D; ++x) {
		os << v(x) << (x == D - 1 ? "" : "\t");
	}
	os << "}";
	return os;
}


struct PlainMat3 {
	__m128 stripes[3];

	float& operator()(int x, int y) {
		return stripes[y].m128_f32[x];
	}
	float operator()(int x, int y) const {
		return stripes[y].m128_f32[x];
	}
};

PlainMat3 operator*(const PlainMat3& lhs, const PlainMat3& rhs) {
	PlainMat3 res;
	//res.stripes[0] = _mm_mul_ps(rhs.stripes[0], _mm_set1_ps(lhs(0, 0)));
	//res.stripes[1] = _mm_mul_ps(rhs.stripes[0], _mm_set1_ps(lhs(0, 1)));
	//res.stripes[2] = _mm_mul_ps(rhs.stripes[0], _mm_set1_ps(lhs(0, 2)));

	//res.stripes[0] = _mm_add_ps(res.stripes[0], _mm_mul_ps(rhs.stripes[1], _mm_set1_ps(lhs(1, 0))));
	//res.stripes[1] = _mm_add_ps(res.stripes[1], _mm_mul_ps(rhs.stripes[1], _mm_set1_ps(lhs(1, 1))));
	//res.stripes[2] = _mm_add_ps(res.stripes[2], _mm_mul_ps(rhs.stripes[1], _mm_set1_ps(lhs(1, 2))));

	//res.stripes[0] = _mm_add_ps(res.stripes[0], _mm_mul_ps(rhs.stripes[2], _mm_set1_ps(lhs(2, 0))));
	//res.stripes[1] = _mm_add_ps(res.stripes[1], _mm_mul_ps(rhs.stripes[2], _mm_set1_ps(lhs(2, 1))));
	//res.stripes[2] = _mm_add_ps(res.stripes[2], _mm_mul_ps(rhs.stripes[2], _mm_set1_ps(lhs(2, 2))));

	__m128 scalarMultiplier;
	for (int y = 0; y < 3; ++y) {
		scalarMultiplier = _mm_set1_ps(lhs(0, y));
		res.stripes[y] = _mm_mul_ps(scalarMultiplier, rhs.stripes[0]);
	}
	for (int x = 1; x < 3; ++x) {
		for (int y = 0; y < 3; ++y) {
			scalarMultiplier = _mm_set1_ps(lhs(x, y));
			__m128 tmp = _mm_mul_ps(scalarMultiplier, rhs.stripes[x]);
			res.stripes[y] = _mm_add_ps(tmp, res.stripes[y]);
		}
	}

	return res;
}

std::ostream& operator<<(std::ostream& os, const PlainMat3& mat) {
	for (int y = 0; y < 3; ++y) {
		os << "{";
		for (int x = 0; x < 3; ++x) {
			os << mat(x, y) << (x == 3 - 1 ? "" : "\t");
		}
		os << "}\n";
	}
	return os;
}


template <class T, int Col1, int Row1, int Col2, int Row2>
double MatMulSpeedTest() {
	using LeftT = Matrix<T, Col1, Row1>;
	using RightT = Matrix<T, Col2, Row2>;
	using ResultT = typename decltype(LeftT()*RightT());

	constexpr int iterationCount = 100'000;
	std::vector<LeftT> left(iterationCount);
	std::vector<RightT> right(iterationCount);
	std::vector<ResultT> result(iterationCount);

	std::minstd_rand rne;
	std::uniform_real_distribution<T> rng;

	for (int i = 0; i < iterationCount; ++i) {
		LeftT& l = left[i];
		RightT& r = right[i];

		for (int x = 0; x < Col1; ++x) {
			for (int y = 0; y < Row1; ++y) {
				l(x, y) = 2;
			}
		}
		for (int x = 0; x < Col2; ++x) {
			for (int y = 0; y < Row2; ++y) {
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


int main() {
	Vector<float, 4> v1{ 1.f, 2.f, 3.f, 4.f };
	Vector<float, 4> v2{ 2.f, 3.f, 4.f, 5.f };
	Vector<float, 4> v3 = v1*v2;

	Vector<float, 8> v4(v2, v3);
	v4.Set(1, 2, 3, 4, 5, 6, 7, 8);
	v4.Set(v3, v2);

	cout << v3 << endl << endl;

	Vector<float, 3> v(1.0f, 2.0f, 3.0f);
	Vector<float, 3> u(1);
	Vector<float, 4> w(v);
	Vector<float, 3> c(v);
	v.Set(1, 2, 3);
	w.Set(2.5f, v);

	auto d = Vector<float, 3>::Cross(v, u);

	using Matrix3x3 = typename std::conditional<false, PlainMat3, Matrix<float, 3, 3>>::type;
	using Matrix4x4 = Matrix<float, 4, 4>;


	// correct test mat4x4
	Matrix4x4 mat1 = {
		1,2,3,3,
		4,5,6,6,
		7,8,9,9,
		7,8,9,9
	};
	Matrix4x4 mat2 = {
		5,4,3,3,
		6,5,4,4,
		7,6,5,5,
		7,6,5,5
	};

	cout << mat1 * mat2 << endl;
	cout << mat1 + mat2 << endl;


	// correct test mat3x3
	Matrix3x3 mat3;
	Matrix3x3 mat4;

	mat3(0, 0) = 1;		mat3(1, 0) = 2;		mat3(2, 0) = 3;
	mat3(0, 1) = 4;		mat3(1, 1) = 5;		mat3(2, 1) = 6;
	mat3(0, 2) = 7;		mat3(1, 2) = 8;		mat3(2, 2) = 9;

	mat4(0, 0) = 5;		mat4(1, 0) = 4;		mat4(2, 0) = 3;
	mat4(0, 1) = 6;		mat4(1, 1) = 5;		mat4(2, 1) = 4;
	mat4(0, 2) = 7;		mat4(1, 2) = 6;		mat4(2, 2) = 5;

	cout << mat3 * mat4 << endl;

	// correct test mat4x2 & 2x4
	Matrix<float, 4, 2> mat7;
	Matrix<float, 2, 4> mat8;

	mat7(0, 0) = 1;		mat7(1, 0) = 2;		mat7(2, 0) = 3;		mat7(3, 0) = 4;
	mat7(0, 1) = 5;		mat7(1, 1) = 6;		mat7(2, 1) = 7;		mat7(3, 1) = 8;

	mat8(0, 0) = 1;		mat8(1, 0) = 5;
	mat8(0, 1) = 2;		mat8(1, 1) = 6;
	mat8(0, 2) = 3;		mat8(1, 2) = 7;
	mat8(0, 3) = 4;		mat8(1, 3) = 8;

	cout << mat7 << " x " << endl;
	cout << mat8 << " = " << endl;
	cout << mat7 * mat8 << " & " << endl;
	cout << mat8 * mat7 << endl;


	double elapsed;
	elapsed = MatMulSpeedTest<float, 2, 2, 2, 2>();
	cout << "time 2x2 x 2x2:\t" << elapsed * 1000 << " ms" << endl;
	elapsed = MatMulSpeedTest<float, 3, 3, 3, 3>();
	cout << "time 3x3 x 3x3:\t" << elapsed * 1000 << " ms" << endl;
	elapsed = MatMulSpeedTest<float, 4, 4, 4, 4>();
	cout << "time 4x4 x 4x4:\t" << elapsed * 1000 << " ms" << endl;
	//elapsed = MatMulSpeedTest<float, 8, 8, 8, 8>();
	//cout << "time 8x8 x 8x8:\t" << elapsed * 1000 << " ms" << endl;
	elapsed = MatMulSpeedTest<float, 4, 2, 2, 4>();
	cout << "time 4x2 x 2x4:\t" << elapsed * 1000 << " ms" << endl;
	elapsed = MatMulSpeedTest<float, 2, 4, 4, 2>();
	cout << "time 2x4 x 4x2:\t" << elapsed * 1000 << " ms" << endl;
	elapsed = MatMulSpeedTest<float, 4, 3, 3, 4>();
	cout << "time 4x3 x 3x4:\t" << elapsed * 1000 << " ms" << endl;
	elapsed = MatMulSpeedTest<float, 3, 4, 4, 3>();
	cout << "time 3x4 x 4x3:\t" << elapsed * 1000 << " ms" << endl;
	

	return 0;
}