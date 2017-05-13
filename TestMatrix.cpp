#pragma warning(disable: 4244)
#include <gtest\gtest.h>

#include "Mathter\Matrix.hpp"

using namespace mathter;

template <class T, int Rows, int Columns, eMatrixOrder Order = eMatrixOrder::FOLLOW_VECTOR, bool Packed = false>
using MatrixC = Matrix<T, Rows, Columns, Order, eMatrixLayout::COLUMN_MAJOR, Packed>;



TEST(Matrix, Ctor_And_Indexer_R) {
	Matrix<float, 3, 3> m = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	Matrix<float, 3, 3> n;
	n(0, 0) = 1;	n(0, 1) = 2;	n(0, 2) = 3;
	n(1, 0) = 4;	n(1, 1) = 5;	n(1, 2) = 6;
	n(2, 0) = 7;	n(2, 1) = 8;	n(2, 2) = 9;

	ASSERT_EQ(m, n);
}


TEST(Matrix, Ctor_And_Indexer_C) {
	MatrixC<float, 3, 3> m = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	MatrixC<float, 3, 3> n;
	n(0, 0) = 1;	n(0, 1) = 2;	n(0, 2) = 3;
	n(1, 0) = 4;	n(1, 1) = 5;	n(1, 2) = 6;
	n(2, 0) = 7;	n(2, 1) = 8;	n(2, 2) = 9;

	ASSERT_EQ(m, n);
}


TEST(Matrix, Add_RR_CC_CR_RC) {
	Matrix<float, 3, 3> m1 = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	MatrixC<float, 3, 3> m2 = {
		7,6,5,
		4,3,2,
		1,0,-1,
	};

	Matrix<float, 3, 3> rexp;
	MatrixC<float, 3, 3> cexp;

	rexp = {
		8,8,8,
		8,8,8,
		8,8,8,
	};
	ASSERT_TRUE(rexp.AlmostEqual(m1 + m2));

	cexp = rexp;
	ASSERT_TRUE(cexp.AlmostEqual(m2 + m1));

	rexp = m1 * 2;
	ASSERT_TRUE(rexp.AlmostEqual(m1 + m1));
	
	cexp = 2 * m2;
	ASSERT_TRUE(cexp.AlmostEqual(m2 + m2));
}


TEST(Matrix, Sub_RR_CC_CR_RC) {
	Matrix<float, 3, 3> m1 = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	MatrixC<float, 3, 3> m2 = {
		2,3,4,
		5,6,7,
		8,9,10,
	};

	Matrix<float, 3, 3> rexp;
	MatrixC<float, 3, 3> cexp;

	rexp = {
		-1, -1, -1,
		-1, -1, -1,
		-1, -1, -1,
	};
	ASSERT_TRUE(rexp.AlmostEqual(m1 - m2));

	cexp = -rexp;
	ASSERT_TRUE(cexp.AlmostEqual(m2 - m1));

	rexp = m1 * 0;
	ASSERT_TRUE(rexp.AlmostEqual(m1 - m1));

	cexp = 0 * m2;
	ASSERT_TRUE(cexp.AlmostEqual(m2 - m2));
}


TEST(Matrix, Multiply_Square_RR) {
	Matrix<float, 3, 3> m = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	m = m*m;
	Matrix<float, 3, 3> mexp = {
		30, 36, 42,
		66, 81, 96,
		102,126,150,
	};


	Matrix<double, 5, 5> n = {
		1,2,3,4,5,
		6,7,8,9,10,
		11,12,13,14,15,
		16,17,18,19,20,
		21,22,23,24,25,
	};
	n = n*n;
	Matrix<double, 5, 5> nexp = {
		215,	230,	245,	260,	275,
		490,	530,	570,	610,	650,
		765,	830,	895,	960,	1025,
		1040,	1130,	1220,	1310,	1400,
		1315,	1430,	1545,	1660,	1775,
	};


	ASSERT_EQ(m, mexp);
	ASSERT_EQ(n, nexp);
}


TEST(Matrix, Multiply_NonSquare_RR) {
	Matrix<float, 4, 2> m = {
		1,2,
		3,4,
		5,6,
		7,8,
	};
	Matrix<float, 2, 4> n = {
		1,2,3,4,
		5,6,7,8,
	};
	Matrix<float, 2, 2> nm = n*m;
	Matrix<float, 4, 4> mn = m*n;

	Matrix<float, 2, 2> nmexp = {
		50,	60,
		114, 140,
	};
	Matrix<float, 4, 4> mnexp = {
		11,	14,	17,	20,
		23,	30,	37,	44,
		35,	46,	57,	68,
		47,	62,	77,	92,
	};

	ASSERT_EQ(mn, mnexp);
	ASSERT_EQ(mn, mnexp);
}


TEST(Matrix, Multiply_Square_CC) {
	MatrixC<float, 3, 3> m = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	m = m*m;
	MatrixC<float, 3, 3> mexp = {
		30, 36, 42,
		66, 81, 96,
		102,126,150,
	};


	MatrixC<double, 5, 5> n = {
		1,2,3,4,5,
		6,7,8,9,10,
		11,12,13,14,15,
		16,17,18,19,20,
		21,22,23,24,25,
	};
	n = n*n;
	MatrixC<double, 5, 5> nexp = {
		215,	230,	245,	260,	275,
		490,	530,	570,	610,	650,
		765,	830,	895,	960,	1025,
		1040,	1130,	1220,	1310,	1400,
		1315,	1430,	1545,	1660,	1775,
	};


	ASSERT_EQ(m, mexp);
	ASSERT_EQ(n, nexp);
}


TEST(Matrix, Multiply_NonSquare_CC) {
	MatrixC<float, 4, 2> m = {
		1,2,
		3,4,
		5,6,
		7,8,
	};
	MatrixC<float, 2, 4> n = {
		1,2,3,4,
		5,6,7,8,
	};
	MatrixC<float, 2, 2> nm = n*m;
	MatrixC<float, 4, 4> mn = m*n;

	MatrixC<float, 2, 2> nmexp = {
		50,	60,
		114, 140,
	};
	MatrixC<float, 4, 4> mnexp = {
		11,	14,	17,	20,
		23,	30,	37,	44,
		35,	46,	57,	68,
		47,	62,	77,	92,
	};

	ASSERT_EQ(mn, mnexp);
	ASSERT_EQ(mn, mnexp);
}


TEST(Matrix, Multiply_Square_RC) {
	Matrix<float, 3, 3> mr = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	MatrixC<float, 3, 3> mc = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	auto m = mr*mc;
	Matrix<float, 3, 3> mexp = {
		30, 36, 42,
		66, 81, 96,
		102,126,150,
	};


	Matrix<double, 5, 5> nr = {
		1,2,3,4,5,
		6,7,8,9,10,
		11,12,13,14,15,
		16,17,18,19,20,
		21,22,23,24,25,
	};
	MatrixC<double, 5, 5> nc = {
		1,2,3,4,5,
		6,7,8,9,10,
		11,12,13,14,15,
		16,17,18,19,20,
		21,22,23,24,25,
	};
	auto n = nr*nc;
	Matrix<double, 5, 5> nexp = {
		215,	230,	245,	260,	275,
		490,	530,	570,	610,	650,
		765,	830,	895,	960,	1025,
		1040,	1130,	1220,	1310,	1400,
		1315,	1430,	1545,	1660,	1775,
	};


	ASSERT_EQ(m, mexp);
	ASSERT_EQ(n, nexp);
}



TEST(Matrix, Identity) {
	Matrix<float, 3, 3> m = Matrix<float, 3, 3>::Identity();
	Matrix<float, 3, 3> mexp = {
		1,0,0,
		0,1,0,
		0,0,1,
	};

	ASSERT_EQ(m, mexp);
}

TEST(Matrix, Zero) {
	Matrix<float, 3, 4> m = Matrix<float, 3, 4>::Zero();
	Matrix<float, 3, 4> mexp = {
		0,0,0,0,
		0,0,0,0,
		0,0,0,0,
	};

	ASSERT_EQ(m, mexp);
}


TEST(Matrix, LU_Decomp) {
	Matrix<float, 3, 3> A = {
		3, -0.1f, -0.2f,
		0.1f, 7, -0.3f,
		0.3f, -0.2f, 10
	};

	Matrix<float, 3, 3> L, U;
	A.DecomposeLU(L, U);

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

	ASSERT_EQ(mT, mexp);
}


TEST(Matrix, Determinant) {
	Matrix<float, 3, 3> m = {
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
}


TEST(Matrix, Trace) {
	Matrix<float, 3, 3> m = {
		1,3,2,
		4,5,6,
		7,8,9,
	};
	float trace = m.Trace();

	ASSERT_FLOAT_EQ(trace, 15.f);
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
	auto m = Matrix<float, 3, 3>::RotationX(1.f);
	Matrix<float, 3, 3> mexp = {
		1.000000, 0.000000, 0.000000, 0.000000, 0.540302, 0.841471, 0.000000, -0.841471, 0.540302
	};
	ASSERT_TRUE(m.AlmostEqual(mexp));


	m = Matrix<float, 3, 3>::RotationY(1.f);
	mexp = {
		0.540302, 0.000000, -0.841471, 0.000000, 1.000000, 0.000000, 0.841471, 0.000000, 0.540302
	};
	ASSERT_TRUE(m.AlmostEqual(mexp));


	m = Matrix<float, 3, 3>::RotationZ(1.f);
	mexp = {
		0.540302, 0.841471, 0.000000, -0.841471, 0.540302, 0.000000, 0.000000, 0.000000, 1.000000
	};
	ASSERT_TRUE(m.AlmostEqual(mexp));
}


TEST(Matrix, RotationAxisAngle) {
	auto m = Matrix<float, 3, 3>::RotationAxisAngle(Vector<float, 3>(1,2,3).Normalized(), 1.0f);
	Matrix<float, 3, 3> mexp = {
		0.573138, 0.740349, -0.351279, -0.609007, 0.671645, 0.421906, 0.548292, -0.027879, 0.835822
	};
	ASSERT_TRUE(m.AlmostEqual(mexp));
}


TEST(Matrix, Scale) {
	auto m = Matrix<float, 5, 5>::Scale(1, 2, 3, 4, 5);
	Vector<float, 5> v(2, 6, 3, 7, 5);

	auto vt1 = v*Vector<float, 5>{ 1, 2, 3, 4, 5 };
	auto vt2 = v*m;

	ASSERT_EQ(vt1, vt2);
}


TEST(Matrix, Translation) {
	auto m = Matrix<float, 6, 5>::Translation(Vector<float, 5>{ 1,2,3,4,5 });
	auto m2 = Matrix<float, 3, 3>::Translation(1, 2);
	auto m3 = Matrix<float, 3, 3>::Translation(Vector<float, 2>(1, 2));
	Vector<float, 5> v(1,2,3,4,5);
	v = (v|1)*m;

	Vector<float, 5> vexp(2,4,6,8,10);

	ASSERT_EQ(v, vexp);
}


TEST(Matrix, Perspective) {
	Vector<float, 4> worldFrustum[2] = {
		{ -0.25f, -0.44444444f, 0.5f, 1 },
		{ 5.0f, 8.8888888f, 10.f, 1 }
	};
	Vector<float, 4> ndcFrustum[2];

	// Z forward
	auto m = Matrix<float, 4, 4>::Perspective(53.13010235/180.f*3.1415926f, 16.f/9.f, 0.5, 10, 0, 1);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	ASSERT_TRUE(ndcFrustum[0].AlmostEqual({ -1, -1, 0, 1 }));
	ASSERT_TRUE(ndcFrustum[1].AlmostEqual({ 1, 1, 1, 1 }));

	// Z backward in NDC
	m = Matrix<float, 4, 4>::Perspective(53.13010235 / 180.f*3.1415926f, 16.f / 9.f, 0.5, 10, 1, -1);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	ASSERT_TRUE(ndcFrustum[0].AlmostEqual({ -1, -1, 1, 1 }));
	ASSERT_TRUE(ndcFrustum[1].AlmostEqual({ 1, 1, -1, 1 }));

	// Z backward in world
	m = Matrix<float, 4, 4>::Perspective(53.13010235 / 180.f*3.1415926f, 16.f / 9.f, -0.5, -10, 0, 1);
	worldFrustum[0].z *= -1;
	worldFrustum[1].z *= -1;
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	ASSERT_TRUE(ndcFrustum[0].AlmostEqual({ -1, -1, 0, 1 }));
	ASSERT_TRUE(ndcFrustum[1].AlmostEqual({ 1, 1, 1, 1 }));

	// Z backward in world && NDC
	m = Matrix<float, 4, 4>::Perspective(53.13010235 / 180.f*3.1415926f, 16.f / 9.f, -0.5, -10, 1, -1);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	ASSERT_TRUE(ndcFrustum[0].AlmostEqual({ -1, -1, 1, 1 }));
	ASSERT_TRUE(ndcFrustum[1].AlmostEqual({ 1, 1, -1, 1 }));
}