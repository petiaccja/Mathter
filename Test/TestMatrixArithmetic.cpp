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


//------------------------------------------------------------------------------
// Matrix-matrix
//------------------------------------------------------------------------------

TEST_CASE_VARIANT("Matrix - Addition", "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> m1 = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
		MatrixT<3, 3> m2 = { 7, 6, 5, 4, 3, 2, 1, 0, -1 };
		decltype(m1 + m2) rexp1 = { 8, 8, 8, 8, 8, 8, 8, 8, 8 };

		MatrixT<4, 5> m3 = { 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4 };
		MatrixT<4, 5> m4 = { 4, 3, 2, 1, 4, 3, 2, 1, 4, 3, 2, 1, 4, 3, 2, 1, 4, 3, 2, 1 };
		decltype(m3 + m4) rexp2 = { 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 };


		MatrixT<2, 2> m5 = { 1, 2, 3, 4 };
		MatrixT<2, 2> m6 = { 4, 3, 2, 1 };
		decltype(m5 + m6) rexp3 = { 5, 5, 5, 5 };

		REQUIRE(m1 + m2 == rexp1);
		REQUIRE(m3 + m4 == rexp2);
		REQUIRE(m5 + m6 == rexp3);
	}
}


TEST_CASE_VARIANT("Matrix - Subtraction", "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> m1 = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
		MatrixT<3, 3> m2 = { 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		decltype(m1 - m2) rexp1 = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };

		MatrixT<2, 2> m3 = { 1, 2, 3, 4 };
		MatrixT<2, 2> m4 = { 2, 3, 4, 5 };
		decltype(m3 - m4) rexp2 = { -1, -1, -1, -1 };

		REQUIRE(m1 - m2 == rexp1);
		REQUIRE(m3 - m4 == rexp2);
	}
}


TEST_CASE_VARIANT("Matrix - Multiply square", "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<2, 2> m2 = {
			1, 2,
			3, 4
		};
		MatrixT<2, 2> n2 = {
			5, 6,
			7, 8
		};
		decltype(m2 * n2) exp2 = {
			19, 22,
			43, 50
		};

		REQUIRE(m2 * n2 == exp2);

		MatrixT<3, 3> m = {
			1, 2, 3,
			4, 5, 6,
			7, 8, 9
		};
		MatrixT<3, 3> n = {
			5, 6, 8,
			1, 3, 5,
			7, 8, 4
		};
		decltype(m * n) exp = {
			28, 36, 30,
			67, 87, 81,
			106, 138, 132
		};

		REQUIRE(m * n == exp);

		MatrixT<5, 5> m5 = {
			1, 2, 3, 4, 5,
			6, 7, 8, 9, 10,
			11, 12, 13, 14, 15,
			16, 17, 18, 19, 20,
			21, 22, 23, 24, 25
		};
		MatrixT<5, 5> n5 = {
			9, 8, 7, 6, 5,
			4, 2, 7, 3, 5,
			3, 6, 2, 7, 2,
			9, 4, 1, 4, 7,
			5, 7, 5, 5, 1
		};
		decltype(m5 * n5) exp5 = {
			87, 81, 56, 74, 54,
			237, 216, 166, 199, 154,
			387, 351, 276, 324, 254,
			537, 486, 386, 449, 354,
			687, 621, 496, 574, 454
		};

		REQUIRE(m5 * n5 == exp5);
	}
}

TEST_CASE_VARIANT("Matrix - Multiply arbitrary", "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<2, 4> m2 = {
			1, 2, 3, 4,
			3, 4, 5, 6
		};
		MatrixT<4, 2> n2 = {
			5,
			6,
			7,
			8,
			6,
			4,
			4,
			9,
		};
		decltype(m2 * n2) exp21 = {
			53, 70,
			97, 124
		};
		decltype(n2 * m2) exp22 = {
			23, 34, 45, 56,
			31, 46, 61, 76,
			18, 28, 38, 48,
			31, 44, 57, 70
		};

		REQUIRE(m2 * n2 == exp21);
		REQUIRE(n2 * m2 == exp22);
	}
}



//------------------------------------------------------------------------------
// Matrix-matrix opposite layout
//------------------------------------------------------------------------------

TEST_CASE_VARIANT("Matrix - Addition opposite layout", "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> m1 = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
		MatrixTOL<3, 3> m2 = { 7, 6, 5, 4, 3, 2, 1, 0, -1 };
		decltype(m1 + m2) rexp1 = { 8, 8, 8, 8, 8, 8, 8, 8, 8 };

		MatrixT<4, 5> m3 = { 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4 };
		MatrixTOL<4, 5> m4 = { 4, 3, 2, 1, 4, 3, 2, 1, 4, 3, 2, 1, 4, 3, 2, 1, 4, 3, 2, 1 };
		decltype(m3 + m4) rexp2 = { 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 };


		MatrixT<2, 2> m5 = { 1, 2, 3, 4 };
		MatrixTOL<2, 2> m6 = { 4, 3, 2, 1 };
		decltype(m5 + m6) rexp3 = { 5, 5, 5, 5 };

		REQUIRE(m1 + m2 == rexp1);
		REQUIRE(m3 + m4 == rexp2);
		REQUIRE(m5 + m6 == rexp3);
	}
}


TEST_CASE_VARIANT("Matrix - Subtraction opposite layout", "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> m1 = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
		MatrixTOL<3, 3> m2 = { 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		decltype(m1 - m2) rexp1 = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };

		MatrixT<2, 2> m3 = { 1, 2, 3, 4 };
		MatrixTOL<2, 2> m4 = { 2, 3, 4, 5 };
		decltype(m3 - m4) rexp2 = { -1, -1, -1, -1 };

		REQUIRE(m1 - m2 == rexp1);
		REQUIRE(m3 - m4 == rexp2);
	}
}


TEST_CASE_VARIANT("Matrix - Multiply square opposite layout", "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<2, 2> m2 = {
			1, 2,
			3, 4
		};
		MatrixTOL<2, 2> n2 = {
			5, 6,
			7, 8
		};
		decltype(m2 * n2) exp2 = {
			19, 22,
			43, 50
		};

		REQUIRE(m2 * n2 == exp2);

		MatrixT<3, 3> m = {
			1, 2, 3,
			4, 5, 6,
			7, 8, 9
		};
		MatrixTOL<3, 3> n = {
			5, 6, 8,
			1, 3, 5,
			7, 8, 4
		};
		decltype(m * n) exp = {
			28, 36, 30,
			67, 87, 81,
			106, 138, 132
		};

		REQUIRE(m * n == exp);

		MatrixT<5, 5> m5 = {
			1, 2, 3, 4, 5,
			6, 7, 8, 9, 10,
			11, 12, 13, 14, 15,
			16, 17, 18, 19, 20,
			21, 22, 23, 24, 25
		};
		MatrixTOL<5, 5> n5 = {
			9, 8, 7, 6, 5,
			4, 2, 7, 3, 5,
			3, 6, 2, 7, 2,
			9, 4, 1, 4, 7,
			5, 7, 5, 5, 1
		};
		decltype(m5 * n5) exp5 = {
			87, 81, 56, 74, 54,
			237, 216, 166, 199, 154,
			387, 351, 276, 324, 254,
			537, 486, 386, 449, 354,
			687, 621, 496, 574, 454
		};

		REQUIRE(m5 * n5 == exp5);
	}
}

TEST_CASE_VARIANT("Matrix - Multiply arbitrary opposite layout", "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<2, 4> m2 = {
			1, 2, 3, 4,
			3, 4, 5, 6
		};
		MatrixTOL<4, 2> n2 = {
			5,
			6,
			7,
			8,
			6,
			4,
			4,
			9,
		};
		decltype(m2 * n2) exp21 = {
			53, 70,
			97, 124
		};
		decltype(n2 * m2) exp22 = {
			23, 34, 45, 56,
			31, 46, 61, 76,
			18, 28, 38, 48,
			31, 44, 57, 70
		};

		REQUIRE(m2 * n2 == exp21);
		REQUIRE(n2 * m2 == exp22);
	}
}



//------------------------------------------------------------------------------
// Matrix-scalar
//------------------------------------------------------------------------------

#define NAME "asd"
#define OPERATOR *

#define TEST_CASE_MATRIX_SCALAR_OP(NAME, OPERATOR)                                                                   \
	TEST_CASE_VARIANT("Matrix - Matrix-scalar " NAME, "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) { \
		const MatrixT<2, 2> sm = {                                                                                   \
			1, 2,                                                                                                    \
			3, 4                                                                                                     \
		};                                                                                                           \
                                                                                                                     \
		const MatrixT<2, 5> m = {                                                                                    \
			1, 2, 5, 6, 9,                                                                                           \
			3, 4, 7, 8, 10                                                                                           \
		};                                                                                                           \
                                                                                                                     \
		const Type b = 27;                                                                                           \
		auto smr = sm OPERATOR b;                                                                                    \
		auto mr = m OPERATOR b;                                                                                      \
                                                                                                                     \
		for (auto i = 0; i < sm.RowCount(); ++i) {                                                                   \
			for (auto j = 0; j < sm.ColumnCount(); ++j) {                                                            \
				REQUIRE(Approx(sm(i, j) OPERATOR b) == smr(i, j));                                                         \
			}                                                                                                        \
		}                                                                                                            \
		for (auto i = 0; i < m.RowCount(); ++i) {                                                                    \
			for (auto j = 0; j < m.ColumnCount(); ++j) {                                                             \
				REQUIRE(Approx(m(i, j) OPERATOR b) == mr(i, j));                                                           \
			}                                                                                                        \
		}                                                                                                            \
	}

#define TEST_CASE_SCALAR_MATRIX_OP(NAME, OPERATOR)                                                                   \
	TEST_CASE_VARIANT("Matrix - Scalar-matrix " NAME, "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) { \
		const MatrixT<2, 2> sm = {                                                                                   \
			1, 2,                                                                                                    \
			3, 4                                                                                                     \
		};                                                                                                           \
                                                                                                                     \
		const MatrixT<2, 5> m = {                                                                                    \
			1, 2, 5, 6, 9,                                                                                           \
			3, 4, 7, 8, 10                                                                                           \
		};                                                                                                           \
                                                                                                                     \
		const Type b = 27;                                                                                           \
		auto smr = b OPERATOR sm;                                                                                    \
		auto mr = b OPERATOR m;                                                                                      \
                                                                                                                     \
		for (auto i = 0; i < sm.RowCount(); ++i) {                                                                   \
			for (auto j = 0; j < sm.ColumnCount(); ++j) {                                                            \
				REQUIRE(Approx(b OPERATOR sm(i, j)) == smr(i, j));                                                          \
			}                                                                                                        \
		}                                                                                                            \
		for (auto i = 0; i < m.RowCount(); ++i) {                                                                    \
			for (auto j = 0; j < m.ColumnCount(); ++j) {                                                             \
				REQUIRE(Approx(b OPERATOR m(i, j)) == mr(i, j));                                                            \
			}                                                                                                        \
		}                                                                                                            \
	}


#define TEST_CASE_MATRIX_SCALAR_COMPOUND_OP(NAME, OPERATOR)                                                                   \
	TEST_CASE_VARIANT("Matrix - Matrix-scalar compound " NAME, "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) { \
		const MatrixT<2, 2> sm = {                                                                                            \
			1, 2,                                                                                                             \
			3, 4                                                                                                              \
		};                                                                                                                    \
                                                                                                                              \
		const MatrixT<2, 5> m = {                                                                                             \
			1, 2, 5, 6, 9,                                                                                                    \
			3, 4, 7, 8, 10                                                                                                    \
		};                                                                                                                    \
                                                                                                                              \
		const Type b = 27;                                                                                                    \
		auto smr = sm;                                                                                                        \
		auto mr = m;                                                                                                          \
		smr OPERATOR b;                                                                                                       \
		mr OPERATOR b;                                                                                                        \
                                                                                                                              \
		for (auto i = 0; i < sm.RowCount(); ++i) {                                                                            \
			for (auto j = 0; j < sm.ColumnCount(); ++j) {                                                                     \
				auto elem = sm(i, j);                                                                                         \
				elem OPERATOR b;                                                                                              \
				REQUIRE(Approx(elem) == smr(i, j));                                                                                   \
			}                                                                                                                 \
		}                                                                                                                     \
		for (auto i = 0; i < m.RowCount(); ++i) {                                                                             \
			for (auto j = 0; j < m.ColumnCount(); ++j) {                                                                      \
				auto elem = m(i, j);                                                                                         \
				elem OPERATOR b;                                                                                              \
				REQUIRE(Approx(elem) == mr(i, j));                                                                                   \
			}                                                                                                                 \
		}                                                                                                                     \
	}

TEST_CASE_MATRIX_SCALAR_OP("mul", *)
TEST_CASE_MATRIX_SCALAR_OP("div", /)

TEST_CASE_SCALAR_MATRIX_OP("mul", *)
TEST_CASE_SCALAR_MATRIX_OP("div", /)

TEST_CASE_MATRIX_SCALAR_COMPOUND_OP("mul", *=)
TEST_CASE_MATRIX_SCALAR_COMPOUND_OP("div", /=)