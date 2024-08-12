// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Cases.hpp"

#include <Mathter/Matrix/Arithmetic.hpp>

#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Matrix - Algo: solve lower triangular (precede)", "[Matrix]",
						decltype(BinaryCaseList<MatrixCaseList<ScalarsFloating, OrdersPrecede, LayoutsAll, PackingsAll>,
												MatrixCaseList<ScalarsFloating, OrdersPrecede, LayoutsAll, PackingsAll>>{})) {
	using MatCoeff = typename TestType::Lhs::template Matrix<3, 3>;
	using MatRhs = typename TestType::Lhs::template Matrix<3, 1>;

	const MatCoeff L = {
		2, 0, 0,
		-2, -1, 0,
		6, 2, -2
	};
	const MatRhs b = {
		8,
		-14,
		20
	};

	const auto x = SolveLowerTriangular(L, b);
	REQUIRE(x(0, 0) == 4);
	REQUIRE(x(1, 0) == 6);
	REQUIRE(x(2, 0) == 8);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Algo: solve upper triangular (precede)", "[Matrix]",
						decltype(BinaryCaseList<MatrixCaseList<ScalarsFloating, OrdersPrecede, LayoutsAll, PackingsAll>,
												MatrixCaseList<ScalarsFloating, OrdersPrecede, LayoutsAll, PackingsAll>>{})) {
	using MatCoeff = typename TestType::Lhs::template Matrix<3, 3>;
	using MatRhs = typename TestType::Lhs::template Matrix<3, 1>;

	const MatCoeff U = {
		-2, 2, 0,
		0, 2, -2,
		0, 0, -2
	};
	const MatRhs b = {
		4,
		-4,
		-16,
	};

	const auto x = SolveUpperTriangular(U, b);
	REQUIRE(x(0, 0) == 4);
	REQUIRE(x(1, 0) == 6);
	REQUIRE(x(2, 0) == 8);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Algo: solve lower triangular (follow)", "[Matrix]",
						decltype(BinaryCaseList<MatrixCaseList<ScalarsFloating, OrdersFollow, LayoutsAll, PackingsAll>,
												MatrixCaseList<ScalarsFloating, OrdersFollow, LayoutsAll, PackingsAll>>{})) {
	using MatCoeff = typename TestType::Lhs::template Matrix<3, 3>;
	using MatLhs = typename TestType::Lhs::template Matrix<1, 3>;

	const MatCoeff L = {
		2, 0, 0,
		-4, 1, 0,
		0, -1, 2
	};
	const MatLhs b = {
		-16, -2, 16
	};

	const auto x = SolveLowerTriangular(L, b);
	REQUIRE(x(0, 0) == 4);
	REQUIRE(x(0, 1) == 6);
	REQUIRE(x(0, 2) == 8);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Algo: solve upper triangular (follow)", "[Matrix]",
						decltype(BinaryCaseList<MatrixCaseList<ScalarsFloating, OrdersFollow, LayoutsAll, PackingsAll>,
												MatrixCaseList<ScalarsFloating, OrdersFollow, LayoutsAll, PackingsAll>>{})) {
	using MatCoeff = typename TestType::Lhs::template Matrix<3, 3>;
	using MatLhs = typename TestType::Lhs::template Matrix<1, 3>;

	const MatCoeff U = {
		3, -1, 0,
		0, 1, -4,
		0, 0, 2
	};
	const MatLhs b = {
		12, 2, -8
	};

	const auto x = SolveUpperTriangular(U, b);
	REQUIRE(x(0, 0) == 4);
	REQUIRE(x(0, 1) == 6);
	REQUIRE(x(0, 2) == 8);
}