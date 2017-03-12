#include "Mathter\Vector.hpp"
#include "Mathter\Matrix.hpp"

#include <iostream>
#include <chrono>
#include <vector>
using namespace std;


template <class T, int Columns, int Rows, eMatrixLayout Layout, eMatrixOrder Order>
std::ostream& operator<<(std::ostream& os, const Matrix<T, Columns, Rows, Layout, Order>& mat) {
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


int main() {
	Vector<float, 4> v1{ 1.f, 2.f, 3.f, 4.f };
	Vector<float, 4> v2{ 2.f, 3.f, 4.f, 5.f };
	Vector<float, 4> v3 = v1*v2;

	cout << v3 << endl << endl;

	Vector<float, 3> v(1.0f, 2.0f, 3.0f);
	Vector<float, 3> u(1);
	Vector<float, 4> w(v, 1.2f);
	Vector<float, 3> c(v);
	v.Set(1, 2, 3);
	w.Set(2.5f, v);

	auto d = Vector<float, 3>::Cross(v, u);

	using Matrix3x3 = typename std::conditional<false, PlainMat3, Matrix<float, 3, 3, eMatrixLayout::ROW_MAJOR>>::type;
	using Matrix4x4 = Matrix<float, 4, 4, eMatrixLayout::ROW_MAJOR>;


	// correct test mat4x4
	Matrix4x4 mat1;
	Matrix4x4 mat2;

	mat1(0, 0) = 1;		mat1(1, 0) = 2;		mat1(2, 0) = 3;		mat1(3, 0) = 3;
	mat1(0, 1) = 4;		mat1(1, 1) = 5;		mat1(2, 1) = 6;		mat1(3, 1) = 6;
	mat1(0, 2) = 7;		mat1(1, 2) = 8;		mat1(2, 2) = 9;		mat1(3, 2) = 9;
	mat1(0, 3) = 7;		mat1(1, 3) = 8;		mat1(2, 3) = 9;		mat1(3, 3) = 9;

	mat2(0, 0) = 5;		mat2(1, 0) = 4;		mat2(2, 0) = 3;		mat2(3, 0) = 3;
	mat2(0, 1) = 6;		mat2(1, 1) = 5;		mat2(2, 1) = 4;		mat2(3, 1) = 4;
	mat2(0, 2) = 7;		mat2(1, 2) = 6;		mat2(2, 2) = 5;		mat2(3, 2) = 5;
	mat2(0, 3) = 7;		mat2(1, 3) = 6;		mat2(2, 3) = 5;		mat2(3, 3) = 5;

	cout << mat1 * mat2 << endl;


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


	// speed test mat4x4
	constexpr int count = 1000000;
	std::chrono::high_resolution_clock::time_point startTime;
	std::chrono::high_resolution_clock::time_point endTime;

	//std::vector<Matrix4x4> mats1(count);
	//std::vector<Matrix4x4> mats2(count);
	//std::vector<Matrix4x4> mats3(count);

	//srand(chrono::high_resolution_clock::now().time_since_epoch().count());
	//for (int k = 0; k < count; ++k) {
	//	for (int i = 0; i < 4; ++i) {
	//		for (int j = 0; j < 4; ++j) {
	//			mats1[k](i, j) = rand() % 100; 
	//			mats2[k](i, j) = rand() % 100;
	//		}
	//	}
	//}

	//startTime = std::chrono::high_resolution_clock::now();
	//for (int i = 0; i < 10; ++i) {
	//	for (int i = 0; i < count; ++i) {
	//		mats3[i] = mats1[1] * mats2[i];
	//	}
	//}
	//endTime = std::chrono::high_resolution_clock::now();

	//mats1 = mats2 = mats3 = {};

	//cout << "elapsed = " << chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count() * 1e-6 << " us" << endl;


	// speed test mat3x3
	std::vector<Matrix3x3> mats4(count);
	std::vector<Matrix3x3> mats5(count);
	std::vector<Matrix3x3> mats6(count);

	for (int k = 0; k < count; ++k) {
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				mats4[k](i, j) = rand() % 100;
				mats5[k](i, j) = rand() % 100;
			}
		}
	}

	startTime = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 10; ++i) {
		for (int i = 0; i < count; ++i) {
			mats6[i] = mats4[1] * mats5[i];
		}
	}
	endTime = std::chrono::high_resolution_clock::now();

	mats4 = mats5 = mats6 = {};

	cout << "elapsed = " << chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count() * 1e-6 << " us" << endl;

	return 0;
}