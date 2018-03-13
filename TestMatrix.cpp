//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)

#include <Catch2/catch.hpp>

#include "Mathter/Matrix.hpp"
#include "TestGenerators.hpp"

#include <random>
#include <complex>

using namespace mathter;


int Ranint() {
	static std::mt19937 rne;
	static std::uniform_int_distribution<int> rng(0, 100);
	return rng(rne);
}


//------------------------------------------------------------------------------
// Matrix tests
//------------------------------------------------------------------------------


TEST_CASE_VARIANT("Matrix - Constructor & indexer", "[Matrix]", TypesAll, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> m = {
			Type(1), Type(2), Type(3),
			Type(4), Type(5), Type(6),
			Type(7), Type(8), Type(9),
		};
		MatrixT<3, 3> n;
		n(0, 0) = Type(1);	n(0, 1) = Type(2);	n(0, 2) = Type(3);
		n(1, 0) = Type(4);	n(1, 1) = Type(5);	n(1, 2) = Type(6);
		n(2, 0) = Type(7);	n(2, 1) = Type(8);	n(2, 2) = Type(9);

		REQUIRE(m == n);
	}
}


TEST_CASE_VARIANT_2("Matrix - Addition", "[Matrix]",
	TypesFloating, OrdersFollow, LayoutsAll, PackedFalse,
	TypesFloating, OrdersFollow, LayoutsAll, PackedFalse)
{
	MatrixT1<3, 3> m1 = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	MatrixT2<3, 3> m2 = {
		7,6,5,
		4,3,2,
		1,0,-1,
	};

	decltype(m1 + m2) rexp = {
		8,8,8,
		8,8,8,
		8,8,8,
	};

	REQUIRE(m1 + m2 == rexp);
}


TEST_CASE_VARIANT_2("Matrix - Subtraction", "[Matrix]",
	TypesFloating, OrdersFollow, LayoutsAll, PackedFalse,
	TypesFloating, OrdersFollow, LayoutsAll, PackedFalse)
{
	MatrixT1<3, 3> m1 = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	MatrixT2<3, 3> m2 = {
		2,3,4,
		5,6,7,
		8,9,10,
	};

	decltype(m1 - m2) rexp = {
		-1, -1, -1,
		-1, -1, -1,
		-1, -1, -1,
	};

	REQUIRE(m1 - m2 == rexp);
}


TEST_CASE_VARIANT_2("Matrix - Multiply square (unpacked)", "[Matrix]",
	TypesFloating, OrdersFollow, LayoutsAll, PackedFalse,
	TypesFloating, OrdersFollow, LayoutsAll, PackedFalse)
{
	SECTION(SECTIONNAME2) {
		MatrixT1<3, 3> m = {
			1,	2,	3,
			4,	5,	6,
			7,	8,	9
		};
		MatrixT2<3, 3> n = {
			5,	6,	8,
			1,	3,	5,
			7,	8,	4
		};
		decltype(m*n) exp = {
			28,	36,	30,
			67,	87,	81,
			106,138,132
		};

		REQUIRE(m*n == exp);

		MatrixT1<5, 5> m5 = {
			1,	2,	3,	4,	5 ,
			6,	7,	8,	9,	10,
			11,	12,	13,	14,	15,
			16,	17,	18,	19,	20,
			21,	22,	23,	24,	25
		};
		MatrixT2<5, 5> n5 = {
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

		REQUIRE(m5*n5 == exp5);
	}
}

TEST_CASE_VARIANT_2("Matrix - Multiply square (packed)", "[Matrix]",
	TypesFloating, OrdersFollow, LayoutsAll, PackedTrue,
	TypesFloating, OrdersFollow, LayoutsAll, PackedTrue)
{
	SECTION(SECTIONNAME2) {
		MatrixT1<3, 3> m = {
			1,	2,	3,
			4,	5,	6,
			7,	8,	9
		};
		MatrixT2<3, 3> n = {
			5,	6,	8,
			1,	3,	5,
			7,	8,	4
		};
		decltype(m*n) exp = {
			28,	36,	30,
			67,	87,	81,
			106,138,132
		};

		REQUIRE(m*n == exp);

		MatrixT1<5, 5> m5 = {
			1,	2,	3,	4,	5 ,
			6,	7,	8,	9,	10,
			11,	12,	13,	14,	15,
			16,	17,	18,	19,	20,
			21,	22,	23,	24,	25
		};
		MatrixT2<5, 5> n5 = {
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

		REQUIRE(m5*n5 == exp5);
	}
}


TEST_CASE("Matrix - Identity", "[Matrix]") {
	Matrix<float, 3, 3> m = Matrix<float, 3, 3>::Identity();
	Matrix<float, 3, 3> mexp = {
		1,0,0,
		0,1,0,
		0,0,1,
	};
	
	REQUIRE(m == mexp);

	Matrix<float, 3, 5> m5 = Matrix<float, 3, 5>::Identity();
	Matrix<float, 3, 5> mexp5 = {
		1,0,0,0,0,
		0,1,0,0,0,
		0,0,1,0,0,
	};

	REQUIRE(m5 == mexp5);
}

TEST_CASE("Matrix - Zero", "[Matrix]") {
	Matrix<float, 3, 4> m = Matrix<float, 3, 4>::Zero();
	Matrix<float, 3, 4> mexp = {
		0,0,0,0,
		0,0,0,0,
		0,0,0,0,
	};

	REQUIRE(m == mexp);
}


TEST_CASE_VARIANT("Matrix - LU decomposition", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> A = {
			3, -0.1, -0.2,
			0.1, 7, -0.3,
			0.3, -0.2, 10
		};

		MatrixT<3, 3> L, U;
		A.DecomposeLU(L, U);

		for (int i = 0; i < A.RowCount(); ++i) {
			for (int j = 0; j < i - 1; ++j) {
				REQUIRE(U(i, j) == Approx(0.0));
				REQUIRE(L(j, i) == Approx(0.0));
			}
		}

		auto Mprod = L * U;
		REQUIRE(A.Approx() == Mprod);
	}
}


TEST_CASE_VARIANT("Matrix - LU solve", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> A = {
			3, -0.1f, -0.2f,
			0.1f, 7, -0.3f,
			0.3f, -0.2f, 10
		};
		Vector<Type, 3, Packed> b = { 7.85, -19.3, 71.4 };
		Vector<Type, 3, Packed> x;
		Vector<Type, 3, Packed> xexp = { 3, -2.5, 7 };

		x = A.DecompositionLU().Solve(b);
		REQUIRE(x.Approx() == xexp);
	}
}


TEST_CASE_VARIANT("Matrix - LUP decomposition", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> A = {
			3, -0.1f, -0.2f,
			0.3f, -0.2f, 10,
			0.1f, 7, -0.3f,
		};

		MatrixT<3, 3> L, U;
		Vector<int, 3, false> P;
		A.DecomposeLUP(L, U, P);

		for (int i = 0; i < A.RowCount(); ++i) {
			for (int j = 0; j < i - 1; ++j) {
				REQUIRE(U(i, j) == Approx(0.0f));
				REQUIRE(L(j, i) == Approx(0.0f));
			}
		}

		MatrixT<3, 3> Pm = decltype(Pm)::Zero();
		for (int i : P) {
			Pm(i, P(i)) = 1.0f;
		}

		auto Mprod = Pm.Transposed()*L*U;
		REQUIRE(A.Approx() == Mprod);
	}
}


TEST_CASE_VARIANT("Matrix - LUP solve", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<4, 4> A = {
			1,3,4,6,
			3,6,2,6,
			9,2,6,7,
			6,2,7,5,
		};
		Vector<Type, 4, Packed> b = { 3,4,2,8 };
		Vector<Type, 4, Packed> x;
		Vector<Type, 4, Packed> xexp = { -94.f / 497, 895.f / 497, 1000.f / 497, -850.f / 497 };

		x = A.DecompositionLUP().Solve(b);

		REQUIRE(x.Approx() == xexp);
	}
}


TEST_CASE("Matrix - LUP decomposition singular", "[Matrix]") {
	Matrix<float, 3, 3> A = {
		1, 0, 0,
		0, 0, 1,
		0, -1, 0,
	};

	Matrix<float, 3, 3> L, U;
	Vector<int, 3, false> P;
	A.DecomposeLUP(L, U, P);

	for (int i = 0; i < A.RowCount(); ++i) {
		for (int j = 0; j < i - 1; ++j) {
			REQUIRE(U(i, j) == Approx(0.0f));
			REQUIRE(L(j, i) == Approx(0.0f));
		}
	}

	Matrix<float, 3, 3> Pm = decltype(Pm)::Zero();
	for (int i : P) {
		Pm(i, P(i)) = 1.0f;
	}

	auto Mprod = Pm.Transposed()*L*U;
	REQUIRE(A.Approx() == Mprod);
}


TEST_CASE_VARIANT("Matrix - QR decomposition", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		// example from wikipedia SVD article
		MatrixT<5, 4> A1 = MatrixT<4, 5>{
			1, 0, 0, 1, 2,
			0, 0, 3, 0, 0,
			0, 0, 0, 0, 0,
			0, 2, 0, 0, 0,
		}.Transposed();
		MatrixT<5, 4> R1;
		MatrixT<5, 5> Q1;

		A1.DecomposeQR(Q1, R1);
		MatrixT<5, 4> A1assembled = Q1 * R1;
		REQUIRE(A1assembled.Approx() == A1);


		// the same matrix as the LU
		MatrixT<3, 3> A2 = {
			3, -0.1f, -0.2f,
			0.1f, 7, -0.3f,
			0.3f, -0.2f, 10
			//12, -51, 4,
			//6, 167, -68,
			//-4, 24, -41
		};
		MatrixT<3, 3> R2;
		MatrixT<3, 3> Q2;

		A2.DecomposeQR(Q2, R2);

		MatrixT<3, 3> A2assembled = Q2 * R2;
		REQUIRE(A2assembled.Approx() == A2);
	}
}


TEST_CASE_VARIANT("Matrix - SVD", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		// example from wikipedia SVD article
		MatrixT<5, 4> A1 = MatrixT<4, 5>{
			1, 0, 0, 1, 2,
			0, 0, 3, 0, 0,
			0, 0, 0, 0, 0,
			0, 2, 0, 0, 0,
		}.Transposed();
		MatrixT<5, 4> U1;
		MatrixT<4, 4> S1;
		MatrixT<4, 4> V1;

		A1.DecomposeSVD(U1, S1, V1);
		auto A1assembled = U1 * S1*V1;
		REQUIRE(A1.Approx() == A1assembled);

		MatrixT<4, 4> U1T;
		MatrixT<4, 4> S1T;
		MatrixT<4, 5> V1T;
		A1.Transposed().DecomposeSVD(U1T, S1T, V1T);
		auto A1Tassembled = U1T * S1T*V1T;
		REQUIRE(A1Tassembled.Approx() == A1.Transposed());

		// the same matrix as the LU
		MatrixT<3, 3> A2 = {
			3, -0.1f, -0.2f,
			0.1f, 7, -0.3f,
			0.3f, -0.2f, 10
		};

		MatrixT<3, 3> U2, S2, V2;
		A2.DecomposeSVD(U2, S2, V2);
		auto A2assembled = U2 * S2*V2;
		REQUIRE(A2assembled.Approx() == A2);
	}
}


TEST_CASE_VARIANT("Matrix - Transpose", "[Matrix]", TypesAll, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<4, 2> m = {
			1,2,
			3,4,
			5,6,
			7,8,
		};
		MatrixT<2, 4> mT = m.Transposed();
		MatrixT<2, 4> mexp = {
			1,3,5,7,
			2,4,6,8,
		};

		REQUIRE(mT == mexp);
	}
}


TEST_CASE_VARIANT("Matrix - Determinant", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> m = {
			1,3,2,
			4,5,6,
			7,8,9,
		};
		auto det = m.Determinant();

		REQUIRE(Approx(det) == 9.0f);

		MatrixT<5, 5> m5 = {
			5, 7, 3, 6, 4,
			4, 7, 4, 6, 3,
			6, 2, 8, 9, 7,
			1, 2, 7, 4, 8,
			5, 9, 7, 1, 5
		};
		det = m5.Determinant();
		REQUIRE(Approx(det) == 4134);
	}
}


TEST_CASE_VARIANT("Matrix - Trace", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> m = {
			1,3,2,
			4,5,6,
			7,8,9,
		};
		auto trace = m.Trace();

		REQUIRE(Approx(trace) == 15.f);

		MatrixT<5, 5> m5 = {
			5, 7, 3, 6, 4,
			4, 7, 4, 6, 3,
			6, 2, 8, 9, 7,
			1, 2, 7, 4, 8,
			5, 9, 7, 1, 5
		};
		trace = m5.Trace();
		REQUIRE(Approx(trace) == 29);
	}
}


TEST_CASE_VARIANT("Matrix - Inverse", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> m = {
			1,3,2,
			4,5,6,
			7,8,9,
		};
		MatrixT<3, 3> mI = m.Inverse();
		MatrixT<3, 3> mexp = {
			-0.333333,	-1.222222,	0.888889,
			0.666667,	-0.555556,	0.222222,
			-0.333333,	1.444444,	-0.777778,
		};

		MatrixT<5, 5> n = {
			1,56,8,4,3,
			4,2,7,8,4,
			1,5,7,4,3,
			9,5,3,8,4,
			7,2,83,46,4,
		};
		MatrixT<5, 5> nI = n.Inverse();
		MatrixT<5, 5> iden = n * nI;
		MatrixT<5, 5> idenexp;
		idenexp.SetIdentity();

		REQUIRE(mexp.Approx() == mI);
		REQUIRE(idenexp.Approx() == iden);
	}
}


TEST_CASE("Matrix - Rotation2D", "[Matrix]") {
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
	
	REQUIRE(m.Approx() == mexp);
	REQUIRE(m3.Approx() == m3exp);
}


TEST_CASE("Matrix - RotationPrincipal", "[Matrix]") {
	auto m1 = Matrix<float, 3, 3>::RotationX(1.f);
	Matrix<float, 3, 3> mexp = {
		1.000000, 0.000000, 0.000000,
		0.000000, 0.540302, 0.841471,
		0.000000, -0.841471, 0.540302
	};
	REQUIRE(m1.Approx() == mexp);


	auto m2 = Matrix<float, 4, 3>::RotationY(1.f);
	Matrix<float, 4, 3> m2exp = {
		0.540302, 0.000000, -0.841471,
		0.000000, 1.000000, 0.000000,
		0.841471, 0.000000, 0.540302,
		0,			 0,			0
	};
	REQUIRE(m2.Approx() == m2exp);


	auto m3 = Matrix<float, 3, 4, eMatrixOrder::PRECEDE_VECTOR>::RotationZ(1.f);
	Matrix<float, 4, 3> m3exp = {
		0.540302, 0.841471, 0.000000,
		-0.841471, 0.540302, 0.000000,
		0.000000, 0.000000, 1.000000,
		0,			0,			0
	};
	REQUIRE(m3.Approx() == m3exp.Transposed());

	auto m4 = Matrix<float, 4, 4>::RotationZ(1.f);
	Matrix<float, 4, 4> m4exp = {
		0.540302, 0.841471, 0.000000, 0,
		-0.841471, 0.540302, 0.000000, 0,
		0.000000, 0.000000, 1.000000, 0,
		0,0,0,1
	};
	REQUIRE(m4.Approx() == m4exp);
}


TEST_CASE("Matrix - RotationAxisAngle", "[Matrix]") {
	auto m = Matrix<float, 3, 3>::RotationAxisAngle(Vector<float, 3>(1,2,3).Normalized(), 1.0f);
	Matrix<float, 3, 3> mexp = {
		0.573138, 0.740349, -0.351279, -0.609007, 0.671645, 0.421906, 0.548292, -0.027879, 0.835822
	};
	REQUIRE(m.Approx() == mexp);

	auto m4 = Matrix<float, 4, 4>::RotationAxisAngle(Vector<float, 3>(1, 2, 3).Normalized(), 1.0f);
	Matrix<float, 4, 4> m4exp = {
		0.573138, 0.740349, -0.351279, 0,
		-0.609007, 0.671645, 0.421906, 0,
		0.548292, -0.027879, 0.835822, 0,
		0,		0,			0,		1
	};
	REQUIRE(m4.Approx() == m4exp);
}


TEST_CASE("Matrix - Scale", "[Matrix]") {
	auto m = Matrix<float, 5, 5>::Scale(1, 2, 3, 4, 5);
	Vector<float, 5> v(2, 6, 3, 7, 5);
	auto m3 = Matrix<float, 3, 3>::Scale(Vector<float, 3>{1, 2, 3});

	auto vt1 = v*Vector<float, 5>{ 1, 2, 3, 4, 5 };
	auto vt2 = v*m;

	REQUIRE(vt1 == vt2);
}


TEST_CASE("Matrix - Translation", "[Matrix]") {
	auto m33 = Matrix<float, 3, 3>::Translation(1, 2);
	auto m = Matrix<float, 6, 5>::Translation(Vector<float, 5>{ 1,2,3,4,5 });
	Vector<float, 5> v(1,2,3,4,5);
	v = v*m;
	Vector<float, 5> vexp(2,4,6,8,10);
	REQUIRE(v == vexp);
	
	auto m2 = Matrix<float, 3, 3>::Translation(1, 2);
	Matrix<float, 3, 3> m2exp = {
		1,0,0,
		0,1,0,
		1,2,1,
	};
	REQUIRE(m2 == m2exp);

	auto m3 = Matrix<float, 3, 2>::Translation(Vector<float, 2>(1, 2));
	Matrix<float, 3, 2> m3exp = {
		1,0,
		0,1,
		1,2,
	};
	REQUIRE(m3 == m3exp);

	auto m4 = Matrix<float, 2, 3, eMatrixOrder::PRECEDE_VECTOR>::Translation(Vector<float, 2>(1, 2));
	Matrix<float, 2, 3, eMatrixOrder::PRECEDE_VECTOR> m4exp = {
		1,0,1,
		0,1,2
	};
	REQUIRE(m4 == m4exp);
}


TEST_CASE("Matrix - Perspective", "[Matrix]") {
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

	REQUIRE(ndcFrustum[0].Approx() == Vector<float, 4>{ -1, -1, 0, 1 });
	REQUIRE(ndcFrustum[1].Approx() == Vector<float, 4>{ 1, 1, 1, 1 });

	// Z backward in NDC
	m = Matrix<float, 4, 4>::Perspective(53.13010235f / 180.f*3.1415926f, 16.f / 9.f, 0.5f, 10, 1, -1);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	REQUIRE(ndcFrustum[0].Approx() == Vector<float, 4>{ -1, -1, 1, 1 });
	REQUIRE(ndcFrustum[1].Approx() == Vector<float, 4>{ 1, 1, -1, 1 });

	// Z backward in world
	m = Matrix<float, 4, 4>::Perspective(53.13010235f / 180.f*3.1415926f, 16.f / 9.f, -0.5f, -10, 0, 1);
	worldFrustum[0].z *= -1;
	worldFrustum[1].z *= -1;
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	REQUIRE(ndcFrustum[0].Approx() == Vector<float, 4>{ -1, -1, 0, 1 });
	REQUIRE(ndcFrustum[1].Approx() == Vector<float, 4>{ 1, 1, 1, 1 });

	// Z backward in world && NDC
	m = Matrix<float, 4, 4>::Perspective(53.13010235f / 180.f*3.1415926f, 16.f / 9.f, -0.5f, -10, 1, -1);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	REQUIRE(ndcFrustum[0].Approx() == Vector<float, 4>{ -1, -1, 1, 1 });
	REQUIRE(ndcFrustum[1].Approx() == Vector<float, 4>{ 1, 1, -1, 1 });
}

TEST_CASE("Matrix - Orthographic", "[Matrix]") {
	Vector<float, 3> worldFrustum[2] = {
		{ -0.25f, -0.44444444f, 0.5f },
		{ 5.0f, 8.8888888f, 10.f }
	};
	Vector<float, 3> ndcFrustum[2];

	// Z forward
	auto m = Matrix<float, 4, 4>::Orthographic(worldFrustum[0], worldFrustum[1], 0, 1);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;

	REQUIRE(ndcFrustum[0].Approx() == Vector<float, 3>{ -1, -1, 0 });
	REQUIRE(ndcFrustum[1].Approx() == Vector<float, 3>{ 1, 1, 1 });
}


TEST_CASE("Matrix - View", "[Matrix]") {
	auto m = Matrix<float, 4, 4>::LookAt({ -6,-5,-5 }, { -1,0,0 }, Vector<float, 3>{0, 0, 1});

	Vector<float, 3> p = { 0, -1, 0 };
	Vector<float, 3> pt = p*m;
	Vector<float, 3> pexp = { sqrt(2),0,8.66025403 };
	REQUIRE(pexp.Approx() == pt);	

	m = Matrix<float, 4, 4>::LookAt({ 0,0,0 }, { 0,5,0 }, Vector<float, 3>{0, 0, 1}, true, false, false);
	p = { 1, 4, 1 };
	pt = p*m;
	pexp = { 1, 1, 4 };

	REQUIRE(pexp.Approx() == pt);

	m = Matrix<float, 4, 4>::LookAt({ 0,0,0 }, { 5,0,0 }, Vector<float, 3>{0, 0, 1}, true, false, false);
	p = { 1, 4, 1 };
	pt = p*m;
	pexp = { -4, 1, 1 };

	REQUIRE(pexp.Approx() == pt);
}


TEST_CASE("Matrix - Submatrix", "[Matrix]") {
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
	REQUIRE(m2 == r);

	m2.Column(4) = Vector<float, 5>('0');
	r(0, 4) = r(1, 4) = r(2, 4) = r(3, 4) = r(4, 4) = '0';
	REQUIRE(m2 == r);


	Vector<char, 3> v = m1.Submatrix<3, 1>(0, 0);
	Vector<char, 3> vr = { 'a', 'f', 'k' };
	REQUIRE(v == vr);
	v = m1.Submatrix<1, 3>(0, 0);
	vr = {'a', 'b', 'c'};
	REQUIRE(v == vr);
	

	// compile error as it should be
	// v = m1.Submatrix<2, 3>(0, 0);


	// compile error as it should be 
	//const Matrix<char, 5, 5>& m2c = m2;
	//m2c.Submatrix<3, 3>(2, 0) = m1.Submatrix<3, 3>(0, 2);

}


TEST_CASE("Matrix - IOParse", "[Matrix]") {
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
		REQUIRE(end != c.c_str());
		REQUIRE(parsed == expected);
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
		REQUIRE(end == c.c_str());
	}
}