// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../TestGenerators.hpp"

#include <Mathter/Common/Approx.hpp>
#include <Mathter/Common/TypeTraits.hpp>
#include <Mathter/Matrix.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>



using namespace mathter;
using Catch::Approx;


//------------------------------------------------------------------------------
// Matrix-matrix
//------------------------------------------------------------------------------

using TypeListAll = TestTypeList<TypesAll, PackedAll, OrdersAll, LayoutsAll>;
using TypeListFloating = TestTypeList<TypesFloating, PackedAll, OrdersAll, LayoutsAll>;


TEMPLATE_LIST_TEST_CASE("Matrix - Add", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using M33I = invert_layout_t<M33>;
		M33 a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

		// Regular layout
		M33 b = { 5, 7, 3, 6, 2, 7, 3, 7, 3 };
		decltype(a + b) r = { 6, 9, 6, 10, 7, 13, 10, 15, 12 };
		REQUIRE(a + b == r);

		// Inverted layout
		M33 bi = { 5, 7, 3, 6, 2, 7, 3, 7, 3 };
		decltype(a + bi) ri = { 6, 9, 6, 10, 7, 13, 10, 15, 12 };
		REQUIRE(a + bi == ri);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Sub", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using M33I = invert_layout_t<M33>;
		M33 a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

		// Regular layout
		M33 b = { 5, 7, 3, 6, 2, 7, 3, 7, 3 };
		decltype(a - b) r = { -4, -5, 0, -2, 3, -1, 4, 1, 6 };
		REQUIRE(a - b == r);

		// Inverted layout
		M33 bi = { 5, 7, 3, 6, 2, 7, 3, 7, 3 };
		decltype(a - bi) ri = { -4, -5, 0, -2, 3, -1, 4, 1, 6 };
		REQUIRE(a - bi == ri);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Mul square", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using M33I = invert_layout_t<M33>;
		M33 a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

		// Regular layout
		M33 b = { 5, 7, 3, 6, 2, 7, 3, 7, 3 };
		decltype(a * b) r = { 26, 32, 26, 68, 80, 65, 110, 128, 104 };
		REQUIRE(a * b == r);

		// Inverted layout
		M33 bi = { 5, 7, 3, 6, 2, 7, 3, 7, 3 };
		decltype(a * bi) ri = { 26, 32, 26, 68, 80, 65, 110, 128, 104 };
		REQUIRE(a * bi == ri);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Mul non-square", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M35 = typename TestType::template Matrix<3, 5>;
		using M52 = typename TestType::template Matrix<5, 2>;
		using M52I = invert_layout_t<M52>;
		M35 a = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

		// Regular layout
		M52 b = { 5, 7, 3, 6, 2, 7, 3, 7, 3, 9 };
		decltype(a * b) r = { 44, 113, 124, 293, 204, 473 };
		REQUIRE(a * b == r);

		// Inverted layout
		M52I bi = { 5, 7, 3, 6, 2, 7, 3, 7, 3, 9 };
		decltype(a * bi) ri = { 44, 113, 124, 293, 204, 473 };
		REQUIRE(a * bi == ri);
	}
}


//------------------------------------------------------------------------------
// Matrix-scalar
//------------------------------------------------------------------------------

template <class MatrixT, class Fun>
auto VerifyScalarResults(const MatrixT& original, Fun fun, const MatrixT& result) {
	for (auto i = 0; i < original.RowCount(); ++i) {
		for (auto j = 0; j < original.ColumnCount(); ++j) {
			REQUIRE(Approx(fun(original(i, j))) == result(i, j));
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - scalar ops", "[Matrix]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using T = scalar_type_t<M33>;

		M33 a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
		T s = T(7);

		SECTION("mul") {
			VerifyScalarResults(
				a, [&s](auto v) { return v * s; }, a * s);
		}
		SECTION("div") {
			VerifyScalarResults(
				a, [&s](auto v) { return v / s; }, a / s);
		}
		SECTION("rmul") {
			VerifyScalarResults(
				a, [&s](auto v) { return s * v; }, s * a);
		}
		SECTION("rdiv") {
			VerifyScalarResults(
				a, [&s](auto v) { return s / v; }, s / a);
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - scalar compound ops", "[Matrix]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using T = scalar_type_t<M33>;

		M33 a = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
		T s = T(7);

		SECTION("mul") {
			auto c = a;
			c *= s;
			VerifyScalarResults(
				a, [&s](auto v) { return v * s; }, c);
		}
		SECTION("div") {
			auto c = a;
			c /= s;
			VerifyScalarResults(
				a, [&s](auto v) { return v / s; }, c);
		}
	}
}