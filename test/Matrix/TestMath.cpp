// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Matrix/Comparison.hpp>
#include <Mathter/Matrix/Math.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>


using namespace mathter;
using namespace test_util;
using namespace std::complex_literals;


TEMPLATE_LIST_TEST_CASE("Matrix - Min", "[Matrix]",
						decltype(MatrixCaseList<ScalarsFloatAndInt32, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M23 = typename TestType::template Matrix<2, 3>;

	const M23 m = { -1, 3, 2, 6, 3, 2 };
	REQUIRE(Min(m) == -1);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Max", "[Matrix]",
						decltype(MatrixCaseList<ScalarsFloatAndInt32, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M23 = typename TestType::template Matrix<2, 3>;

	const M23 m = { -1, 3, 2, 6, 3, 2 };
	REQUIRE(Max(m) == 6);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Sum", "[Matrix]",
						decltype(MatrixCaseList<ScalarsFloatAndInt32, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M23 = typename TestType::template Matrix<2, 3>;

	const M23 m = { -1, 3, 2, 6, 3, 2 };
	REQUIRE(Sum(m) == 15);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Abs", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M23 = typename TestType::template Matrix<2, 3>;

	const M23 m = { -1, 3, 2, 6, 3, 2 };
	const auto r = Abs(m);
	static_assert(!is_complex_v<scalar_type_t<std::decay_t<decltype(r)>>>);
	REQUIRE(r(0, 0) == 1);
	REQUIRE(r(0, 1) == 3);
	REQUIRE(r(0, 2) == 2);
	REQUIRE(r(1, 0) == 6);
	REQUIRE(r(1, 1) == 3);
	REQUIRE(r(1, 2) == 2);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Real (of real matrix)", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M23 = typename TestType::template Matrix<2, 3>;

	const M23 m = { -1, 3, 2, 6, 3, 2 };
	const auto r = Real(m);
	static_assert(!is_complex_v<scalar_type_t<std::decay_t<decltype(r)>>>);
	REQUIRE(r(0, 0) == -1);
	REQUIRE(r(0, 1) == 3);
	REQUIRE(r(0, 2) == 2);
	REQUIRE(r(1, 0) == 6);
	REQUIRE(r(1, 1) == 3);
	REQUIRE(r(1, 2) == 2);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Imag (of real matrix)", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M23 = typename TestType::template Matrix<2, 3>;

	const M23 m = { -1, 3, 2, 6, 3, 2 };
	const auto r = Imag(m);
	static_assert(!is_complex_v<scalar_type_t<std::decay_t<decltype(r)>>>);
	REQUIRE(r(0, 0) == 0);
	REQUIRE(r(0, 1) == 0);
	REQUIRE(r(0, 2) == 0);
	REQUIRE(r(1, 0) == 0);
	REQUIRE(r(1, 1) == 0);
	REQUIRE(r(1, 2) == 0);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Real (of complex matrix)", "[Matrix]",
						decltype(MatrixCaseList<ScalarsComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M22 = typename TestType::template Matrix<2, 2>;

	const M22 m = { 1.0f, -2.0if, 2.0f + 1if, 1.1f };
	const auto r = Real(m);
	static_assert(!is_complex_v<scalar_type_t<std::decay_t<decltype(r)>>>);
	REQUIRE(r(0, 0) == 1.0f);
	REQUIRE(r(0, 1) == 0.0f);
	REQUIRE(r(1, 0) == 2.0f);
	REQUIRE(r(1, 1) == 1.1f);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Imag (of complex matrix)", "[Matrix]",
						decltype(MatrixCaseList<ScalarsComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M22 = typename TestType::template Matrix<2, 2>;

	const M22 m = { 1.0f, -2.0if, 2.0f + 1if, 1.1f };
	const auto r = Imag(m);
	static_assert(!is_complex_v<scalar_type_t<std::decay_t<decltype(r)>>>);
	REQUIRE(r(0, 0) == 0.0f);
	REQUIRE(r(0, 1) == -2.0f);
	REQUIRE(r(1, 0) == 1.0f);
	REQUIRE(r(1, 1) == 0.0f);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Norm", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M24 = typename TestType::template Matrix<2, 4>;
	using V8 = Vector<scalar_type_t<M24>, 8>;

	const V8 v = { 1, 2, 3, 4, 5, 6, 7, 8 };
	const M24 m = { 1, 2, 3, 4, 5, 6, 7, 8 };
	REQUIRE(Catch::Approx(Length(v)) == Norm(m));
}


TEMPLATE_LIST_TEST_CASE("Matrix - NormPrecise", "[Matrix]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M22 = typename TestType::template Matrix<2, 2>;

	SECTION("Underflow") {
		const M22 m = { 1e-30f, 1e-30f, 1e-30f, 1e-30f };
		REQUIRE(Catch::Approx(2e-30f) == NormPrecise(m));
	}
	SECTION("Denormal") {
		const M22 m = { 1e-39f, 1e-39f, 1e-39f, 1e-39f };
		REQUIRE(Catch::Approx(2e-39f) == NormPrecise(m));
	}
	SECTION("Zero") {
		const M22 m = { 0, 0, 0, 0 };
		REQUIRE(0 == NormPrecise(m));
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Trace", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M33 = typename TestType::template Matrix<3, 3>;

	const M33 m = {
		1, 3, 2,
		4, 5, 6,
		7, 8, 9
	};

	REQUIRE(std::real(Trace(m)) == Catch::Approx(15.f));
}


TEMPLATE_LIST_TEST_CASE("Matrix - Transpose", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M42 = typename TestType::template Matrix<4, 2>;
	using M24 = typename TestType::template Matrix<2, 4>;

	const M42 value = {
		1, 2,
		3, 4,
		5, 6,
		7, 8
	};
	const M24 transpose = Transpose(value);
	const M24 expected = {
		1, 3, 5, 7,
		2, 4, 6, 8
	};

	REQUIRE(transpose == expected);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Conjugate transpose", "[Matrix]",
						decltype(MatrixCaseList<ScalarsComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M42 = typename TestType::template Matrix<4, 2>;
	using M24 = typename TestType::template Matrix<2, 4>;

	const M42 value = {
		1, 2,
		3, 4.f + 2if,
		5, 6,
		7, 8
	};
	const M24 conjugateTranspose = ConjTranspose(value);
	const M24 expected = {
		1, 3, 5, 7,
		2, 4.f - 2if, 6, 8
	};

	REQUIRE(conjugateTranspose == expected);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Bareiss algorithm", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	SECTION("Zero matrix") {
		using M33 = typename TestType::template Matrix<3, 3>;
		const M33 m2 = {
			0, 0, 0,
			0, 0, 0,
			0, 0, 0
		};
		REQUIRE(std::real(impl::BareissAlgorithm(m2)) == Catch::Approx(0));
	}
	SECTION("No swap") {
		using M33 = typename TestType::template Matrix<3, 3>;
		const M33 m2 = {
			2, 0, 0,
			0, 1, 0,
			0, 0, 3
		};
		REQUIRE(std::real(impl::BareissAlgorithm(m2)) == Catch::Approx(6));
	}
	SECTION("Swap") {
		using M33 = typename TestType::template Matrix<3, 3>;
		const M33 m2 = {
			2, 2, 0,
			1, 1, 1,
			0, 3, 3
		};
		REQUIRE(std::real(impl::BareissAlgorithm(m2)) == Catch::Approx(-6));
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Determinant", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	// Small matrices have hand-rolled explicit formulas that need to be tested separately.
	// Large matrices use general decompositions and solvers.
	SECTION("2x2") {
		using M22 = typename TestType::template Matrix<2, 2>;
		const M22 m2 = {
			1, 3,
			4, 5
		};
		REQUIRE(std::real(Determinant(m2)) == Catch::Approx(-7));
	}
	SECTION("3x3") {
		using M33 = typename TestType::template Matrix<3, 3>;
		const M33 m3 = {
			1, 3, 2,
			4, 5, 6,
			7, 8, 9
		};
		REQUIRE(std::real(Determinant(m3)) == Catch::Approx(9));
	}
	SECTION("4x4") {
		using M44 = typename TestType::template Matrix<4, 4>;
		const M44 m4 = {
			1, 3, 2, 1,
			4, 5, 6, 2,
			7, 8, 9, 3,
			1, 2, 3, 4
		};
		REQUIRE(std::real(Determinant(m4)) == Catch::Approx(27));
	}
	SECTION("Large") {
		using M55 = typename TestType::template Matrix<5, 5>;
		const M55 m5 = {
			5, 7, 3, 6, 4,
			4, 7, 4, 6, 3,
			6, 2, 8, 9, 7,
			1, 2, 7, 4, 8,
			5, 9, 7, 1, 5
		};
		REQUIRE(std::real(Determinant(m5)) == Catch::Approx(4134));
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Inverse", "[Matrix]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	// Small matrices have hand-rolled explicit formulas that need to be tested separately.
	// Large matrices use general decompositions and solvers.
	SECTION("2x2") {
		using M22 = typename TestType::template Matrix<2, 2>;

		const M22 value = {
			1, 3,
			4, 5
		};
		const M22 inverse = Inverse(value);
		const M22 expected = {
			-5 / 7.0, 3 / 7.0,
			4 / 7.0, -1 / 7.0
		};

		REQUIRE(inverse == test_util::Approx(expected));
	}
	SECTION("3x3") {
		using M33 = typename TestType::template Matrix<3, 3>;

		const M33 value = {
			1, 3, 2,
			4, 5, 6,
			7, 8, 9
		};
		const M33 inverse = Inverse(value);
		const M33 expected = {
			-1 / 3.0, -11 / 9.0, 8 / 9.0,
			2 / 3.0, -5 / 9.0, 2 / 9.0,
			-1 / 3.0, 13 / 9.0, -7 / 9.0
		};
		REQUIRE(inverse == test_util::Approx(expected));
	}
	SECTION("4x4") {
		using M44 = typename TestType::template Matrix<4, 4>;

		const M44 value = {
			1, 3, 2, 1,
			4, 5, 6, 2,
			7, 8, 9, 3,
			1, 2, 3, 4
		};
		const M44 inverse = Inverse(value);
		const M44 expected = {
			-1 / 3.0, -35 / 27.0, 25 / 27.0, 1 / 27.0,
			2 / 3.0, -11 / 27.0, 4 / 27.0, -2 / 27.0,
			-1 / 3.0, 43 / 27.0, -23 / 27.0, -2 / 27.0,
			0.0, -2 / 3.0, 1 / 3.0, 1 / 3.0
		};

		REQUIRE(inverse == test_util::Approx(expected));
	}
	SECTION("Large") {
		using M55 = typename TestType::template Matrix<5, 5>;

		const M55 value = {
			1, 56, 8, 4, 3,
			4, 2, 7, 8, 4,
			1, 5, 7, 4, 3,
			9, 5, 3, 8, 4,
			7, 2, 83, 46, 4
		};
		const M55 inverse = Inverse(value);
		const M55 expected = {
			-19 / 649.0, -34648 / 81125.0, 18107 / 81125.0, 22021 / 81125.0, 828 / 81125.0,
			13 / 649.0, 89 / 16225.0, -401 / 16225.0, -28 / 16225.0, -4 / 16225.0,
			-14 / 649.0, -4539 / 16225.0, 4226 / 16225.0, 1428 / 16225.0, 204 / 16225.0,
			1 / 22.0, 1559 / 2750.0, -1481 / 2750.0, -543 / 2750.0, 1 / 2750.0,
			-2 / 59.0, 222 / 7375.0, 3027 / 7375.0, -219 / 7375.0, -242 / 7375.0
		};

		REQUIRE(inverse == test_util::Approx(expected, 1e-6f));
	}
}