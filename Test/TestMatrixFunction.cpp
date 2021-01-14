//L=============================================================================
//L This software is distributed under the MIT license.
//L Copyright 2021 Péter Kardos
//L=============================================================================

#pragma warning(disable : 4244)

#include "../Mathter/Common/Approx.hpp"
#include "../Mathter/Matrix.hpp"
#include "TestGenerators.hpp"

#include <Catch2/catch.hpp>
#include <complex>


using namespace mathter;



TEST_CASE_VARIANT("Matrix - Trace", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> m = {
			1, 3, 2,
			4, 5, 6,
			7, 8, 9
		};
		auto trace = Trace(m);

		REQUIRE(Approx(trace) == 15.f);

		MatrixT<5, 5> m5 = {
			5, 7, 3, 6, 4,
			4, 7, 4, 6, 3,
			6, 2, 8, 9, 7,
			1, 2, 7, 4, 8,
			5, 9, 7, 1, 5
		};
		trace = Trace(m5);
		REQUIRE(Approx(trace) == 29);
	}
}


TEST_CASE_VARIANT("Matrix - Transpose", "[Matrix]", TypesAll, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<4, 2> m = {
			1, 2,
			3, 4,
			5, 6,
			7, 8
		};
		MatrixT<2, 4> mT = Transpose(m);
		MatrixT<2, 4> mexp = {
			1, 3, 5, 7,
			2, 4, 6, 8
		};

		REQUIRE(mT == mexp);
	}
}


TEST_CASE_VARIANT("Matrix - Determinant small matrix", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<2, 2> m2 = {
			1, 3,
			4, 5
		};
		REQUIRE(ApproxVec(Determinant(m2)) == -7);

		MatrixT<4, 4> m4 = {
			1, 3, 2, 1,
			4, 5, 6, 2,
			7, 8, 9, 3,
			1, 2, 3, 4
		};
		REQUIRE(ApproxVec(Determinant(m4)) == 27);



		MatrixT<3, 3> m3 = {
			1, 3, 2,
			4, 5, 6,
			7, 8, 9
		};
		REQUIRE(Approx(Determinant(m3)) == 9);
	}
}


TEST_CASE_VARIANT("Matrix - Determinant", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<5, 5> m5 = {
			5, 7, 3, 6, 4,
			4, 7, 4, 6, 3,
			6, 2, 8, 9, 7,
			1, 2, 7, 4, 8,
			5, 9, 7, 1, 5
		};
		REQUIRE(Approx(Determinant(m5)) == 4134);
	}
}


TEST_CASE_VARIANT("Matrix - Inverse Small Matrix", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<2, 2> m2 = {
			1, 3,
			4, 5
		};
		MatrixT<2, 2> mI2 = Inverse(m2);
		MatrixT<2, 2> mexp2 = {
			-0.714286, 0.428571,
			0.571429, -0.142857
		};

		REQUIRE(ApproxVec(mI2) == mexp2);

		MatrixT<3, 3> m3 = {
			1, 3, 2,
			4, 5, 6,
			7, 8, 9
		};
		MatrixT<3, 3> mI3 = Inverse(m3);
		MatrixT<3, 3> mexp3 = {
			-0.333333, -1.222222, 0.888889,
			0.666667, -0.555556, 0.222222,
			-0.333333, 1.444444, -0.777778
		};

		REQUIRE(ApproxVec(mI3) == mexp3);

		MatrixT<4, 4> m4 = {
			1, 3, 2, 1,
			4, 5, 6, 2,
			7, 8, 9, 3,
			1, 2, 3, 4
		};
		MatrixT<4, 4> mI4 = Inverse(m4);
		MatrixT<4, 4> mexp4 = {
			-0.333333, -1.296296, 0.925926, 0.037037,
			0.666667, -0.407407, 0.148148, -0.074074,
			-0.333333, 1.592593, -0.851852, -0.074074,
			0, -0.666667, 0.333333, 0.333333
		};

		REQUIRE(ApproxVec(mI4) == mexp4);
	}
}


TEST_CASE_VARIANT("Matrix - Inverse", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<5, 5> n = {
			1, 56, 8, 4, 3,
			4, 2, 7, 8, 4,
			1, 5, 7, 4, 3,
			9, 5, 3, 8, 4,
			7, 2, 83, 46, 4
		};
		MatrixT<5, 5> nI = Inverse(n);
		MatrixT<5, 5> iden = n * nI;
		MatrixT<5, 5> idenexp;
		idenexp = Identity();

		REQUIRE(ApproxVec(idenexp) == iden);
	}
}


TEST_CASE_VARIANT("Matrix - Norm", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	Vector<Type, 8> v = { 1, 2, 3, 4, 5, 6, 7, 8 };
	MatrixT<2, 4> m = { 1, 2, 3, 4, 5, 6, 7, 8 };
	REQUIRE(Approx(Length(v)) == Norm(m));
}