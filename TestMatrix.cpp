#include <gtest\gtest.h>

#include "Mathter\Matrix.hpp"

using namespace mathter;



TEST(Matrix, CtorIndex) {
	Matrix<float, 3, 3> m = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	Matrix<float, 3, 3> n;
	n(0, 0) = 1;	n(1, 0) = 2;	n(2, 0) = 3;
	n(0, 1) = 4;	n(1, 1) = 5;	n(2, 1) = 6;
	n(0, 2) = 7;	n(1, 2) = 8;	n(2, 2) = 9;

	ASSERT_EQ(m, n);
}


TEST(Matrix, MulSquare) {
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


TEST(Matrix, MulNonsquare) {
	Matrix<float, 2, 4> m = {
		1,2,
		3,4,
		5,6,
		7,8,
	};
	Matrix<float, 4, 2> n = {
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


TEST(Matrix, MulNonsquareColmajor) {
	Matrix<float, 2, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR> m = {
		1,2,
		3,4,
		5,6,
		7,8,
	};
	Matrix<float, 4, 2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR> n = {
		1,2,3,4,
		5,6,7,8,
	};
	Matrix<float, 2, 2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR> nm = n*m;
	Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR> mn = m*n;

	Matrix<float, 2, 2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR> nmexp = {
		50,	60,
		114, 140,
	};
	Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR> mnexp = {
		11,	14,	17,	20,
		23,	30,	37,	44,
		35,	46,	57,	68,
		47,	62,	77,	92,
	};

	ASSERT_EQ(mn, mnexp);
	ASSERT_EQ(mn, mnexp);
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
	Matrix<float, 4, 3> m = Matrix<float, 4, 3>::Zero();
	Matrix<float, 4, 3> mexp = {
		0,0,0,0,
		0,0,0,0,
		0,0,0,0,
	};

	ASSERT_EQ(m, mexp);
}


TEST(Matrix, Transpose) {
	Matrix<float, 2, 4> m = {
		1,2,
		3,4,
		5,6,
		7,8,
	};
	Matrix<float, 4, 2> mT = m.Transposed();
	Matrix<float, 4, 2> mexp = {
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
	Matrix<float, 3, 3> mI = m.Inverted();
	Matrix<float, 3, 3> mexp = {
		-0.333333,	-1.222222,	0.888889,
		0.666667,	-0.555556,	0.222222,
		-0.333333,	1.444444,	-0.777778,
	};

	ASSERT_TRUE(mexp.AlmostEqual(mI));
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
	auto m = Matrix<float, 5, 6>::Translation(Vector<float, 5>{ 1,2,3,4,5 });
	auto m2 = Matrix<float, 3, 3>::Translation(1, 2);
	Vector<float, 5> v(1,2,3,4,5);
	v = (v|1)*m;

	Vector<float, 5> vexp(2,4,6,8,10);

	ASSERT_EQ(v, vexp);

}