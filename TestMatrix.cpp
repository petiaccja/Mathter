//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)
#include <gtest/gtest.h>

#include "Mathter/Matrix.hpp"

#include <random>

using namespace mathter;


template <class T, int Rows, int Columns, eMatrixOrder Order = eMatrixOrder::FOLLOW_VECTOR, bool Packed = false>
using MatrixC = Matrix<T, Rows, Columns, Order, eMatrixLayout::COLUMN_MAJOR, Packed>;


int Ranint() {
	static std::mt19937 rne;
	static std::uniform_int_distribution<int> rng(0, 100);
	return rng(rne);
}

//------------------------------------------------------------------------------
// Helper macros
//------------------------------------------------------------------------------

// Test for both layouts
#define CALL_LAYOUT(Name, Case) \
TEST(Name, Case) { \
Name##_##Case<eMatrixLayout::ROW_MAJOR>(); \
Name##_##Case<eMatrixLayout::COLUMN_MAJOR>(); \
}

#define TEST_LAYOUT(Name, Case) \
template <eMatrixLayout Layout> \
void Name##_##Case(); \
CALL_LAYOUT(Name, Case) \
template <eMatrixLayout Layout> \
void Name##_##Case()


// Test for both orders
#define CALL_ORDER(Name, Case) \
TEST(Name, Case) { \
Name##_##Case<eMatrixOrder::FOLLOW_VECTOR>(); \
Name##_##Case<eMatrixOrder::PRECEDE_VECTOR>(); \
}

#define TEST_ORDER(Name, Case) \
template <eMatrixOrder Order> \
void Name##_##Case(); \
CALL_ORDER(Name, Case) \
template <eMatrixOrder Order> \
void Name##_##Case()


// Test with every combination of Order and Layout
#define CALL_LAYOUTxORDER(Name, Case) \
TEST(Name, Case) { \
Name##_##Case<eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR>(); \
Name##_##Case<eMatrixOrder::PRECEDE_VECTOR, eMatrixLayout::ROW_MAJOR>(); \
Name##_##Case<eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR>(); \
Name##_##Case<eMatrixOrder::PRECEDE_VECTOR, eMatrixLayout::COLUMN_MAJOR>(); \
}

#define TEST_LAYOUTxORDER(Name, Case) \
template <eMatrixOrder Order, eMatrixLayout Layout> \
void Name##_##Case(); \
CALL_LAYOUTxORDER(Name, Case) \
template <eMatrixOrder Order, eMatrixLayout Layout> \
void Name##_##Case()


// Test with every combination of Layout pairs
#define CALL_LAYOUT_SQUARED(Name, Case) \
TEST(Name, Case) { \
Name##_##Case<eMatrixLayout::ROW_MAJOR,		eMatrixLayout::ROW_MAJOR>(); \
Name##_##Case<eMatrixLayout::ROW_MAJOR,		eMatrixLayout::COLUMN_MAJOR>(); \
Name##_##Case<eMatrixLayout::COLUMN_MAJOR,	eMatrixLayout::ROW_MAJOR>(); \
Name##_##Case<eMatrixLayout::COLUMN_MAJOR,	eMatrixLayout::COLUMN_MAJOR>(); \
}

#define TEST_LAYOUT_SQUARED(Name, Case) \
template <eMatrixLayout Layout1, eMatrixLayout Layout2> \
void Name##_##Case(); \
CALL_LAYOUT_SQUARED(Name, Case) \
template <eMatrixLayout Layout1, eMatrixLayout Layout2> \
void Name##_##Case()

// Test with every combination of Order and Layout pairs
#define CALL_LAYOUTxORDER_SQUARED(Name, Case) \
TEST(Name, Case) { \
Name##_##Case<eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::ROW_MAJOR,		eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::ROW_MAJOR	>(); \
Name##_##Case<eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::ROW_MAJOR,		eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::COLUMN_MAJOR	>(); \
Name##_##Case<eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::ROW_MAJOR,		eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::ROW_MAJOR	>(); \
Name##_##Case<eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::ROW_MAJOR,		eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::COLUMN_MAJOR	>(); \
Name##_##Case<eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::COLUMN_MAJOR,	eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::ROW_MAJOR	>(); \
Name##_##Case<eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::COLUMN_MAJOR,	eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::COLUMN_MAJOR	>(); \
Name##_##Case<eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::COLUMN_MAJOR,	eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::ROW_MAJOR	>(); \
Name##_##Case<eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::COLUMN_MAJOR,	eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::COLUMN_MAJOR	>(); \
Name##_##Case<eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::ROW_MAJOR,		eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::ROW_MAJOR	>(); \
Name##_##Case<eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::ROW_MAJOR,		eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::COLUMN_MAJOR	>(); \
Name##_##Case<eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::ROW_MAJOR,		eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::ROW_MAJOR	>(); \
Name##_##Case<eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::ROW_MAJOR,		eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::COLUMN_MAJOR	>(); \
Name##_##Case<eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::COLUMN_MAJOR,	eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::ROW_MAJOR	>(); \
Name##_##Case<eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::COLUMN_MAJOR,	eMatrixOrder::FOLLOW_VECTOR,	eMatrixLayout::COLUMN_MAJOR	>(); \
Name##_##Case<eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::COLUMN_MAJOR,	eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::ROW_MAJOR	>(); \
Name##_##Case<eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::COLUMN_MAJOR,	eMatrixOrder::PRECEDE_VECTOR,	eMatrixLayout::COLUMN_MAJOR	>(); \
}

#define TEST_LAYOUTxORDER_SQUARED(Name, Case) \
template <eMatrixOrder Order1, eMatrixLayout Layout1, eMatrixOrder Order2, eMatrixLayout Layout2> \
void Name##_##Case(); \
CALL_LAYOUTxORDER_SQUARED(Name, Case) \
template <eMatrixOrder Order1, eMatrixLayout Layout1, eMatrixOrder Order2, eMatrixLayout Layout2> \
void Name##_##Case()



//------------------------------------------------------------------------------
// Matrix tests
//------------------------------------------------------------------------------

TEST_LAYOUTxORDER(Matrix, CtorIndexer) {
	Matrix<float, 3, 3, Order, Layout> m = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	Matrix<float, 3, 3, Order, Layout> n;
	n(0, 0) = 1;	n(0, 1) = 2;	n(0, 2) = 3;
	n(1, 0) = 4;	n(1, 1) = 5;	n(1, 2) = 6;
	n(2, 0) = 7;	n(2, 1) = 8;	n(2, 2) = 9;

	ASSERT_TRUE(m == n);
}



TEST_LAYOUT_SQUARED(Matrix, Add) {
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR, Layout1> m1 = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR, Layout2> m2 = {
		7,6,5,
		4,3,2,
		1,0,-1,
	};

	decltype(m1 + m2) rexp = {
		8,8,8,
		8,8,8,
		8,8,8,
	};

	ASSERT_TRUE(m1 + m2 == rexp);
}


TEST_LAYOUT_SQUARED(Matrix, Sub) {
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR, Layout1> m1 = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR, Layout2> m2 = {
		2,3,4,
		5,6,7,
		8,9,10,
	};

	decltype(m1 - m2) rexp = {
		-1, -1, -1,
		-1, -1, -1,
		-1, -1, -1,
	};

	ASSERT_TRUE(m1 - m2 == rexp);
}


TEST_LAYOUT_SQUARED(Matrix, MulSquare) {
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR, Layout1> m = {
		1,	2,	3,
		4,	5,	6,
		7,	8,	9
	};
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR, Layout2> n = {
		5,	6,	8,
		1,	3,	5,
		7,	8,	4
	};
	decltype(m*n) exp = {
		28,	36,	30,
		67,	87,	81,
		106,138,132
	};

	ASSERT_TRUE(m*n == exp);

	Matrix<float, 5, 5, eMatrixOrder::FOLLOW_VECTOR, Layout1> m5 = {
		1,	2,	3,	4,	5 ,
		6,	7,	8,	9,	10,
		11,	12,	13,	14,	15,
		16,	17,	18,	19,	20,
		21,	22,	23,	24,	25
	};
	Matrix<float, 5, 5, eMatrixOrder::FOLLOW_VECTOR, Layout2> n5 = {
		9,	8,	7,	6,	5,
		4,	2,	7,	3,	5,
		3,	6,	2,	7,	2,
		9,	4,	1,	4,	7,
		5,	7,	5,	5,	1
	};
	decltype(m5*n5) exp5 = {
		87,	81,	56,	74,	54,
		237,216,166,199,154,
		387,351,276,324,254,
		537,486,386,449,354,
		687,621,496,574,454
	};

	ASSERT_TRUE(m5*n5 == exp5);
}


TEST(Matrix, Identity) {
	Matrix<float, 3, 3> m = Matrix<float, 3, 3>::Identity();
	Matrix<float, 3, 3> mexp = {
		1,0,0,
		0,1,0,
		0,0,1,
	};
	
	ASSERT_TRUE(m == mexp);

	Matrix<float, 3, 5> m5 = Matrix<float, 3, 5>::Identity();
	Matrix<float, 3, 5> mexp5 = {
		1,0,0,0,0,
		0,1,0,0,0,
		0,0,1,0,0,
	};

	ASSERT_TRUE(m5 == mexp5);
}

TEST(Matrix, Zero) {
	Matrix<float, 3, 4> m = Matrix<float, 3, 4>::Zero();
	Matrix<float, 3, 4> mexp = {
		0,0,0,0,
		0,0,0,0,
		0,0,0,0,
	};

	ASSERT_TRUE(m == mexp);
}


TEST(Matrix, LU_Decomp) {
	Matrix<float, 3, 3> A = {
		3, -0.1f, -0.2f,
		0.1f, 7, -0.3f,
		0.3f, -0.2f, 10
	};

	Matrix<float, 3, 3> L, U;
	A.DecomposeLU(L, U);

	for (int i = 0; i < A.RowCount(); ++i) {
		for (int j = 0; j < i - 1; ++j) {
			ASSERT_FLOAT_EQ(U(i, j), 0.0f);
			ASSERT_FLOAT_EQ(L(j, i), 0.0f);
		}
	}

	auto Mprod = L*U;
	ASSERT_TRUE(A.AlmostEqual(Mprod));
}


TEST(Matrix, LU_Solve) {
	Matrix<float, 3, 3> A = {
		3, -0.1f, -0.2f,
		0.1f, 7, -0.3f,
		0.3f, -0.2f, 10
	};
	Vector<float, 3> b = { 7.85, -19.3, 71.4 };
	Vector<float, 3> x;
	Vector<float, 3> xexp = {3, -2.5, 7};

	bool solved = A.DecompositionLU().Solve(x, b);

	ASSERT_TRUE(solved);
	ASSERT_TRUE(x.AlmostEqual(xexp));
}


TEST(Matrix, QR_Decomp) {
	// example from wikipedia SVD article
	Matrix<float, 5, 4> A1 = Matrix<float, 4, 5>{
		1, 0, 0, 1, 2,
		0, 0, 3, 0, 0,
		0, 0, 0, 0, 0,
		0, 2, 0, 0, 0,
	}.Transposed();
	Matrix<float, 5, 4> R1;
	Matrix<float, 5, 5> Q1;

	A1.DecomposeQR(Q1, R1);
	Matrix<float, 5, 4> A1assembled = Q1*R1;
	ASSERT_TRUE(A1assembled.AlmostEqual(A1));


	// the same matrix as the LU
	Matrix<float, 3, 3> A2 = {
		3, -0.1f, -0.2f,
		0.1f, 7, -0.3f,
		0.3f, -0.2f, 10
		//12, -51, 4,
		//6, 167, -68,
		//-4, 24, -41
	};
	Matrix<float, 3, 3> R2;
	Matrix<float, 3, 3> Q2;

	A2.DecomposeQR(Q2, R2);

	Matrix<float, 3, 3> A2assembled = Q2*R2;
	ASSERT_TRUE(A2assembled.AlmostEqual(A2));
}


TEST(Matrix, SVD_Decomp) {
	// example from wikipedia SVD article
	Matrix<float, 5, 4> A1 = Matrix<float, 4, 5>{
		1, 0, 0, 1, 2,
		0, 0, 3, 0, 0,
		0, 0, 0, 0, 0,
		0, 2, 0, 0, 0,
	}.Transposed();
	Matrix<float, 5, 4> U1;
	Matrix<float, 4, 4> S1;
	Matrix<float, 4, 4> V1;

	A1.DecomposeSVD(U1, S1, V1);
	auto A1assembled = U1*S1*V1;
	ASSERT_TRUE(A1.AlmostEqual(A1assembled));


	// the same matrix as the LU
	Matrix<float, 3, 3> A2 = {
		3, -0.1f, -0.2f,
		0.1f, 7, -0.3f,
		0.3f, -0.2f, 10
	};

	Matrix<float, 3, 3> U2, S2, V2;
	A2.DecomposeSVD(U2, S2, V2);
	auto A2assembled = U2*S2*V2;
	ASSERT_TRUE(A2assembled.AlmostEqual(A2));

	
	1 == 1;
}


TEST(Matrix, Transpose) {
	Matrix<float, 4, 2> m = {
		1,2,
		3,4,
		5,6,
		7,8,
	};
	Matrix<float, 2, 4> mT = m.Transposed();
	Matrix<float, 2, 4> mexp = {
		1,3,5,7,
		2,4,6,8,
	};

	ASSERT_TRUE(mT == mexp);
}


TEST_LAYOUT(Matrix, Determinant) {
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR, Layout> m = {
		1,3,2,
		4,5,6,
		7,8,9,
	};
	float det = m.Determinant();

	ASSERT_FLOAT_EQ(det, 9.0f);

	m = {
		1,2,3,
		4,5,6,
		7,8,9
	};
	det = m.Determinant();
	ASSERT_FLOAT_EQ(det, 0.0f);

	Matrix<double, 5, 5, eMatrixOrder::FOLLOW_VECTOR, Layout> m5 = {
		5, 7, 3, 6, 4,
		4, 7, 4, 6, 3,
		6, 2, 8, 9, 7,
		1, 2, 7, 4, 8,
		5, 9, 7, 1, 5
	};
	det = m5.Determinant();
	ASSERT_FLOAT_EQ(det, 4134);
}


TEST_LAYOUT(Matrix, Trace) {
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR, Layout> m = {
		1,3,2,
		4,5,6,
		7,8,9,
	};
	float trace = m.Trace();

	ASSERT_FLOAT_EQ(trace, 15.f);

	Matrix<double, 5, 5, eMatrixOrder::FOLLOW_VECTOR, Layout> m5 = {
		5, 7, 3, 6, 4,
		4, 7, 4, 6, 3,
		6, 2, 8, 9, 7,
		1, 2, 7, 4, 8,
		5, 9, 7, 1, 5
	};
	trace = m5.Trace();
	ASSERT_FLOAT_EQ(trace, 29);
}


TEST(Matrix, Inverse) {
	Matrix<float, 3, 3> m = {
		1,3,2,
		4,5,6,
		7,8,9,
	};
	Matrix<float, 3, 3> mI = m.Inverse();
	Matrix<float, 3, 3> mexp = {
		-0.333333,	-1.222222,	0.888889,
		0.666667,	-0.555556,	0.222222,
		-0.333333,	1.444444,	-0.777778,
	};

	Matrix<float, 5, 5> n = {
		1,56,8,4,3,
		4,2,7,8,4,
		1,5,7,4,3,
		9,5,3,8,4,
		7,2,83,46,4,		
	};
	Matrix<float, 5, 5> nI = n.Inverse();
	Matrix<float, 5, 5> iden = n*nI;
	Matrix<float, 5, 5> idenexp;
	idenexp.SetIdentity();

	ASSERT_TRUE(mexp.AlmostEqual(mI));
	ASSERT_TRUE(idenexp.AlmostEqual(iden));
}


TEST(Matrix, Rotation2D) {
	auto m = Matrix<float, 2, 2>::Rotation(1.f);
	auto m3 = Matrix<float, 3, 3>::Rotation(1.f);
	Matrix<float, 2, 2> mexp = {
		0.54030, 0.84147,
		-0.84147, 0.54030
	};
	Matrix<float, 3, 3> m3exp = {
		0.54030, 0.84147, 0,
		-0.84147, 0.54030, 0,
		0,0,1,
	};
	
	ASSERT_TRUE(m.AlmostEqual(mexp));
	ASSERT_TRUE(m3.AlmostEqual(m3exp));
}


TEST(Matrix, RotationPrincipal) {
	auto m1 = Matrix<float, 3, 3>::RotationX(1.f);
	Matrix<float, 3, 3> mexp = {
		1.000000, 0.000000, 0.000000,
		0.000000, 0.540302, 0.841471,
		0.000000, -0.841471, 0.540302
	};
	ASSERT_TRUE(m1.AlmostEqual(mexp));


	auto m2 = Matrix<float, 4, 3>::RotationY(1.f);
	Matrix<float, 4, 3> m2exp = {
		0.540302, 0.000000, -0.841471,
		0.000000, 1.000000, 0.000000,
		0.841471, 0.000000, 0.540302,
		0,			 0,			0
	};
	ASSERT_TRUE(m2.AlmostEqual(m2exp));


	auto m3 = Matrix<float, 3, 4, eMatrixOrder::PRECEDE_VECTOR>::RotationZ(1.f);
	Matrix<float, 4, 3> m3exp = {
		0.540302, 0.841471, 0.000000,
		-0.841471, 0.540302, 0.000000,
		0.000000, 0.000000, 1.000000,
		0,			0,			0
	};
	ASSERT_TRUE(m3.AlmostEqual(m3exp.Transposed()));

	auto m4 = Matrix<float, 4, 4>::RotationZ(1.f);
	Matrix<float, 4, 4> m4exp = {
		0.540302, 0.841471, 0.000000, 0,
		-0.841471, 0.540302, 0.000000, 0,
		0.000000, 0.000000, 1.000000, 0,
		0,0,0,1
	};
	ASSERT_TRUE(m4.AlmostEqual(m4exp));
}


TEST(Matrix, RotationAxisAngle) {
	auto m = Matrix<float, 3, 3>::RotationAxisAngle(Vector<float, 3>(1,2,3).Normalized(), 1.0f);
	Matrix<float, 3, 3> mexp = {
		0.573138, 0.740349, -0.351279, -0.609007, 0.671645, 0.421906, 0.548292, -0.027879, 0.835822
	};
	ASSERT_TRUE(m.AlmostEqual(mexp));

	auto m4 = Matrix<float, 4, 4>::RotationAxisAngle(Vector<float, 3>(1, 2, 3).Normalized(), 1.0f);
	Matrix<float, 4, 4> m4exp = {
		0.573138, 0.740349, -0.351279, 0,
		-0.609007, 0.671645, 0.421906, 0,
		0.548292, -0.027879, 0.835822, 0,
		0,		0,			0,		1
	};
	ASSERT_TRUE(m4.AlmostEqual(m4exp));
}


TEST(Matrix, Scale) {
	auto m = Matrix<float, 5, 5>::Scale(1, 2, 3, 4, 5);
	Vector<float, 5> v(2, 6, 3, 7, 5);
	auto m3 = Matrix<float, 3, 3>::Scale(Vector<float, 3>{1, 2, 3});

	auto vt1 = v*Vector<float, 5>{ 1, 2, 3, 4, 5 };
	auto vt2 = v*m;

	ASSERT_TRUE(vt1 == vt2);
}


TEST(Matrix, Translation) {
	auto m33 = Matrix<float, 3, 3>::Translation(1, 2);
	auto m = Matrix<float, 6, 5>::Translation(Vector<float, 5>{ 1,2,3,4,5 });
	Vector<float, 5> v(1,2,3,4,5);
	v = v*m;
	Vector<float, 5> vexp(2,4,6,8,10);
	ASSERT_TRUE(v == vexp);
	
	auto m2 = Matrix<float, 3, 3>::Translation(1, 2);
	Matrix<float, 3, 3> m2exp = {
		1,0,0,
		0,1,0,
		1,2,1,
	};
	ASSERT_TRUE(m2 == m2exp);

	auto m3 = Matrix<float, 3, 2>::Translation(Vector<float, 2>(1, 2));
	Matrix<float, 3, 2> m3exp = {
		1,0,
		0,1,
		1,2,
	};
	ASSERT_TRUE(m3 == m3exp);

	auto m4 = Matrix<float, 2, 3, eMatrixOrder::PRECEDE_VECTOR>::Translation(Vector<float, 2>(1, 2));
	Matrix<float, 2, 3, eMatrixOrder::PRECEDE_VECTOR> m4exp = {
		1,0,1,
		0,1,2
	};
	ASSERT_TRUE(m4 == m4exp);
}


TEST(Matrix, Perspective) {
	Vector<float, 4> worldFrustum[2] = {
		{ -0.25f, -0.140625f, 0.5f, 1 },
		{ 5.0f, 2.8125f, 10.f, 1 }
	};
	Vector<float, 4> ndcFrustum[2];

	// Z forward
	auto m = Matrix<float, 4, 4>::Perspective(53.13010235f/180.f*3.1415926f, 16.f/9.f, 0.5f, 10, 0, 1);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	ASSERT_TRUE(ndcFrustum[0].AlmostEqual({ -1, -1, 0, 1 }));
	ASSERT_TRUE(ndcFrustum[1].AlmostEqual({ 1, 1, 1, 1 }));

	// Z backward in NDC
	m = Matrix<float, 4, 4>::Perspective(53.13010235f / 180.f*3.1415926f, 16.f / 9.f, 0.5f, 10, 1, -1);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	ASSERT_TRUE(ndcFrustum[0].AlmostEqual({ -1, -1, 1, 1 }));
	ASSERT_TRUE(ndcFrustum[1].AlmostEqual({ 1, 1, -1, 1 }));

	// Z backward in world
	m = Matrix<float, 4, 4>::Perspective(53.13010235f / 180.f*3.1415926f, 16.f / 9.f, -0.5f, -10, 0, 1);
	worldFrustum[0].z *= -1;
	worldFrustum[1].z *= -1;
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	ASSERT_TRUE(ndcFrustum[0].AlmostEqual({ -1, -1, 0, 1 }));
	ASSERT_TRUE(ndcFrustum[1].AlmostEqual({ 1, 1, 1, 1 }));

	// Z backward in world && NDC
	m = Matrix<float, 4, 4>::Perspective(53.13010235f / 180.f*3.1415926f, 16.f / 9.f, -0.5f, -10, 1, -1);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	ASSERT_TRUE(ndcFrustum[0].AlmostEqual({ -1, -1, 1, 1 }));
	ASSERT_TRUE(ndcFrustum[1].AlmostEqual({ 1, 1, -1, 1 }));
}

TEST(Matrix, Orthographic) {
	Vector<float, 3> worldFrustum[2] = {
		{ -0.25f, -0.44444444f, 0.5f },
		{ 5.0f, 8.8888888f, 10.f }
	};
	Vector<float, 3> ndcFrustum[2];

	// Z forward
	auto m = Matrix<float, 4, 4>::Orthographic(worldFrustum[0], worldFrustum[1], 0, 1);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;

	ASSERT_TRUE(ndcFrustum[0].AlmostEqual({ -1, -1, 0 }));
	ASSERT_TRUE(ndcFrustum[1].AlmostEqual({ 1, 1, 1 }));
}


TEST(Matrix, View) {
	auto m = Matrix<float, 4, 4>::LookAt({ -6,-5,-5 }, { -1,0,0 }, Vector<float, 3>{0, 0, 1});

	Vector<float, 3> p = { 0, -1, 0 };
	Vector<float, 3> pt = p*m;
	Vector<float, 3> pexp = { sqrt(2),0,8.66025403 };
	ASSERT_TRUE(pexp.AlmostEqual(pt));	

	m = Matrix<float, 4, 4>::LookAt({ 0,0,0 }, { 0,5,0 }, Vector<float, 3>{0, 0, 1}, true, false, false);
	p = { 1, 4, 1 };
	pt = p*m;
	pexp = { 1, 1, 4 };

	ASSERT_TRUE(pexp.AlmostEqual(pt));

	m = Matrix<float, 4, 4>::LookAt({ 0,0,0 }, { 5,0,0 }, Vector<float, 3>{0, 0, 1}, true, false, false);
	p = { 1, 4, 1 };
	pt = p*m;
	pexp = { -4, 1, 1 };

	ASSERT_TRUE(pexp.AlmostEqual(pt));
}


TEST(Matrix, Submatrix) {
	Matrix<char, 5, 5>  m1 = {
		'a','b','c','d','e',
		'f','g','h','i','j',
		'k','l','m','n','o',
		'p','q','r','s','t',
		'u','v','w','x','y',
	};

	Matrix<char, 5, 5>  m2 = {
		'z','z','z','z','z',
		'z','z','z','z','z',
		'z','z','z','z','z',
		'z','z','z','z','z',
		'z','z','z','z','z',
	};

	Matrix<char, 5, 5>  r = {
		'z','z','z','p','q',
		'z','z','z','u','v',
		'c','d','e','z','z',
		'h','i','j','z','z',
		'm','n','o','z','z',
	};

	Matrix<char, 2, 2> sm = m1.Submatrix<2, 2>(3, 0);
	m2.Submatrix<3, 3>(2, 0) = m1.Submatrix<3, 3>(0, 2);
	m2.Submatrix<2, 2>(0, 3) = sm;
	ASSERT_EQ(m2, r);

	m2.Column(4) = Vector<float, 5>('0');
	r(0, 4) = r(1, 4) = r(2, 4) = r(3, 4) = r(4, 4) = '0';
	ASSERT_EQ(m2, r);


	Vector<char, 3> v = m1.Submatrix<3, 1>(0, 0);
	Vector<char, 3> vr = { 'a', 'f', 'k' };
	ASSERT_EQ(v, vr);
	v = m1.Submatrix<1, 3>(0, 0);
	vr = {'a', 'b', 'c'};
	ASSERT_EQ(v, vr);




	// compile error as it should be
	// v = m1.Submatrix<2, 3>(0, 0);


	// compile error as it should be 
	//const Matrix<char, 5, 5>& m2c = m2;
	//m2c.Submatrix<3, 3>(2, 0) = m1.Submatrix<3, 3>(0, 2);

}


TEST(Matrix, IOParse) {
	Matrix<float, 2, 2> parsed;

	std::string successCases[] = {
		"[3.14, 2.718; 0.57, 6.63]",
		"[3.14  2.718; 0.57  6.63  ]",
		"[[3.14, 2.718];\n[0.57, 6.63]]",
		"3.14, 2.718; 0.57, 6.63   ]",
	};
	Matrix<float, 2, 2> expected{ 
		3.14f, 2.718f, 
		0.57f, 6.63f 
	};

	for (const auto& c : successCases) {
		const char* end;
		parsed = strtomat<decltype(parsed)>(c.c_str(), &end);
		ASSERT_TRUE(end != c.c_str());
		ASSERT_EQ(parsed, expected);
	}

	std::string failureCases[] = {
		"[3.14, 2.718\n 0.57, 6.63]", // missing row delimiter
		"[3.14  2.718h 0.57  6.63  ]", // invalid delimiter
		"[3.14, 2.718; 0.57, 6.63; 1.38, 6.02]", // too many rows
		"[3.14, 2.718 ]", // too few rows
	};

	for (const auto& c : failureCases) {
		const char* end;
		parsed = strtomat<decltype(parsed)>(c.c_str(), &end);
		ASSERT_TRUE(end == c.c_str());
	}
}