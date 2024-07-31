// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Cases.hpp"

#include <Mathter/Matrix/Matrix.hpp>

#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Matrix - Construct default", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;

	// This test seems dumb, but it makes sure the default constructor of the template compiles.
	// Initialization is handled by the underlying Vectors, so no need to test.
	const Mat m;
}


TEMPLATE_LIST_TEST_CASE("Matrix - Construct from elements", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<2, 3>;
	using Scalar = scalar_type_t<Mat>;

	Mat m = {
		1, 2, 3,
		4, 5, 6
	};

	SECTION("Mutable") {
		REQUIRE(m(0, 0) == static_cast<Scalar>(1));
		REQUIRE(m(0, 1) == static_cast<Scalar>(2));
		REQUIRE(m(0, 2) == static_cast<Scalar>(3));
		REQUIRE(m(1, 0) == static_cast<Scalar>(4));
		REQUIRE(m(1, 1) == static_cast<Scalar>(5));
		REQUIRE(m(1, 2) == static_cast<Scalar>(6));
	}
	SECTION("Const") {
		const auto& c = m;
		REQUIRE(c(0, 0) == static_cast<Scalar>(1));
		REQUIRE(c(0, 1) == static_cast<Scalar>(2));
		REQUIRE(c(0, 2) == static_cast<Scalar>(3));
		REQUIRE(c(1, 0) == static_cast<Scalar>(4));
		REQUIRE(c(1, 1) == static_cast<Scalar>(5));
		REQUIRE(c(1, 2) == static_cast<Scalar>(6));
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Construct row/column", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	SECTION("Row") {
		using Mat = typename TestType::template Matrix<1, 3>;
		using Scalar = scalar_type_t<Mat>;

		Mat m = Vector<Scalar, 3>{ 1, 2, 3 };
		const auto& c = m;

		REQUIRE(m(0, 0) == static_cast<Scalar>(1));
		REQUIRE(m(0, 1) == static_cast<Scalar>(2));
		REQUIRE(m(0, 2) == static_cast<Scalar>(3));
		REQUIRE(m(0) == static_cast<Scalar>(1));
		REQUIRE(m(1) == static_cast<Scalar>(2));
		REQUIRE(m(2) == static_cast<Scalar>(3));
		REQUIRE(c(0) == static_cast<Scalar>(1));
		REQUIRE(c(1) == static_cast<Scalar>(2));
		REQUIRE(c(2) == static_cast<Scalar>(3));
	}
	SECTION("Column") {
		using Mat = typename TestType::template Matrix<3, 1>;
		using Scalar = scalar_type_t<Mat>;

		Mat m = Vector<Scalar, 3>{ 1, 2, 3 };
		const auto& c = m;

		REQUIRE(m(0, 0) == static_cast<Scalar>(1));
		REQUIRE(m(1, 0) == static_cast<Scalar>(2));
		REQUIRE(m(2, 0) == static_cast<Scalar>(3));
		REQUIRE(m(0) == static_cast<Scalar>(1));
		REQUIRE(m(1) == static_cast<Scalar>(2));
		REQUIRE(m(2) == static_cast<Scalar>(3));
		REQUIRE(c(0) == static_cast<Scalar>(1));
		REQUIRE(c(1) == static_cast<Scalar>(2));
		REQUIRE(c(2) == static_cast<Scalar>(3));
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Access row/column", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;
	using Scalar = scalar_type_t<Mat>;

	const Mat m = {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9
	};

	SECTION("Row") {
		const auto column = m.Column(1);
		REQUIRE(column[0] == static_cast<Scalar>(2));
		REQUIRE(column[1] == static_cast<Scalar>(5));
		REQUIRE(column[2] == static_cast<Scalar>(8));
	}
	SECTION("Column") {
		const auto row = m.Row(1);
		REQUIRE(row[0] == static_cast<Scalar>(4));
		REQUIRE(row[1] == static_cast<Scalar>(5));
		REQUIRE(row[2] == static_cast<Scalar>(6));
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Modify element", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<2, 2>;
	using Scalar = scalar_type_t<Mat>;

	Mat m = {
		1, 2,
		4, 5
	};


	m(1, 0) = 10;

	REQUIRE(m(0, 0) == static_cast<Scalar>(1));
	REQUIRE(m(0, 1) == static_cast<Scalar>(2));
	REQUIRE(m(1, 0) == static_cast<Scalar>(10));
	REQUIRE(m(1, 1) == static_cast<Scalar>(5));
}


TEMPLATE_LIST_TEST_CASE("Matrix - Modify row/column", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;
	using Scalar = scalar_type_t<Mat>;

	const Mat m = {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9
	};

	SECTION("Row") {
		auto copy = m;
		copy.Column(1, Vector(10, 11, 12));
		REQUIRE(copy(0, 0) == static_cast<Scalar>(1));
		REQUIRE(copy(1, 0) == static_cast<Scalar>(4));
		REQUIRE(copy(2, 0) == static_cast<Scalar>(7));

		REQUIRE(copy(0, 1) == static_cast<Scalar>(10));
		REQUIRE(copy(1, 1) == static_cast<Scalar>(11));
		REQUIRE(copy(2, 1) == static_cast<Scalar>(12));

		REQUIRE(copy(0, 2) == static_cast<Scalar>(3));
		REQUIRE(copy(1, 2) == static_cast<Scalar>(6));
		REQUIRE(copy(2, 2) == static_cast<Scalar>(9));
	}
	SECTION("Column") {
		auto copy = m;
		copy.Row(1, Vector(10, 11, 12));
		REQUIRE(copy(0, 0) == static_cast<Scalar>(1));
		REQUIRE(copy(0, 1) == static_cast<Scalar>(2));
		REQUIRE(copy(0, 2) == static_cast<Scalar>(3));

		REQUIRE(copy(1, 0) == static_cast<Scalar>(10));
		REQUIRE(copy(1, 1) == static_cast<Scalar>(11));
		REQUIRE(copy(1, 2) == static_cast<Scalar>(12));

		REQUIRE(copy(2, 0) == static_cast<Scalar>(7));
		REQUIRE(copy(2, 1) == static_cast<Scalar>(8));
		REQUIRE(copy(2, 2) == static_cast<Scalar>(9));
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Extract submatrix", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<4, 4>;
	using Scalar = scalar_type_t<Mat>;

	const Mat m = {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12,
		13, 14, 15, 16
	};

	const auto submatrix = m.Extract<2, 3>(1, 1);

	REQUIRE(submatrix(0, 0) == static_cast<Scalar>(6));
	REQUIRE(submatrix(0, 1) == static_cast<Scalar>(7));
	REQUIRE(submatrix(0, 2) == static_cast<Scalar>(8));
	REQUIRE(submatrix(1, 0) == static_cast<Scalar>(10));
	REQUIRE(submatrix(1, 1) == static_cast<Scalar>(11));
	REQUIRE(submatrix(1, 2) == static_cast<Scalar>(12));
}


TEMPLATE_LIST_TEST_CASE("Matrix - Insert submatrix", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<4, 4>;
	using Sub = typename TestType::template Matrix<2, 3>;
	using Scalar = scalar_type_t<Mat>;

	Mat m = {
		1, 2, 3, 4,
		5, 0, 0, 0,
		9, 0, 0, 0,
		13, 14, 15, 16
	};

	const Sub submatrix = {
		6, 7, 8,
		10, 11, 12
	};

	m.Insert(1, 1, submatrix);

	REQUIRE(m(0, 0) == static_cast<Scalar>(1));
	REQUIRE(m(0, 1) == static_cast<Scalar>(2));
	REQUIRE(m(0, 2) == static_cast<Scalar>(3));
	REQUIRE(m(0, 3) == static_cast<Scalar>(4));

	REQUIRE(m(1, 0) == static_cast<Scalar>(5));
	REQUIRE(m(1, 1) == static_cast<Scalar>(6));
	REQUIRE(m(1, 2) == static_cast<Scalar>(7));
	REQUIRE(m(1, 3) == static_cast<Scalar>(8));

	REQUIRE(m(2, 0) == static_cast<Scalar>(9));
	REQUIRE(m(2, 1) == static_cast<Scalar>(10));
	REQUIRE(m(2, 2) == static_cast<Scalar>(11));
	REQUIRE(m(2, 3) == static_cast<Scalar>(12));

	REQUIRE(m(3, 0) == static_cast<Scalar>(13));
	REQUIRE(m(3, 1) == static_cast<Scalar>(14));
	REQUIRE(m(3, 2) == static_cast<Scalar>(15));
	REQUIRE(m(3, 3) == static_cast<Scalar>(16));
}


TEMPLATE_LIST_TEST_CASE("Matrix - Convert to vector", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	SECTION("Row") {
		using Mat = typename TestType::template Matrix<1, 3>;
		using Scalar = scalar_type_t<Mat>;

		const Mat m = { 1, 2, 3 };
		const Vector v = m;

		REQUIRE(v[0] == static_cast<Scalar>(1));
		REQUIRE(v[1] == static_cast<Scalar>(2));
		REQUIRE(v[2] == static_cast<Scalar>(3));
	}
	SECTION("Column") {
		using Mat = typename TestType::template Matrix<3, 1>;
		using Scalar = scalar_type_t<Mat>;

		const Mat m = { 1, 2, 3 };
		const Vector v = m;

		REQUIRE(v[0] == static_cast<Scalar>(1));
		REQUIRE(v[1] == static_cast<Scalar>(2));
		REQUIRE(v[2] == static_cast<Scalar>(3));
	}
}