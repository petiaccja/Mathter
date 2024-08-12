// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Cases.hpp"

#include <Mathter/Matrix/Comparison.hpp>

#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Matrix - Compare (same order)", "[Matrix]",
						decltype(BinaryCaseList<MatrixCaseList<ScalarsFloatAndInt32, OrdersFollow, LayoutsAll, PackingsAll>,
												MatrixCaseList<ScalarsFloatAndInt32, OrdersFollow, LayoutsAll, PackingsAll>>{})) {
	using MatLhs = typename TestType::Lhs::template Matrix<2, 3>;
	using MatRhs = typename TestType::Rhs::template Matrix<2, 3>;

	const MatLhs lhs = { 1, 2, 3, 4, 5, 6 };
	const MatRhs same = { 1, 2, 3, 4, 5, 6 };
	const MatRhs different = { 1, 3, 3, 4, 5, 6 };

	REQUIRE(lhs == same);
	REQUIRE(lhs != different);
}


TEMPLATE_LIST_TEST_CASE("Matrix - Compare (different order)", "[Matrix]",
						decltype(BinaryCaseList<MatrixCaseList<ScalarsFloatAndInt32, OrdersPrecede, LayoutsAll, PackingsAll>,
												MatrixCaseList<ScalarsFloatAndInt32, OrdersFollow, LayoutsAll, PackingsAll>>{})) {
	using MatLhs = typename TestType::Lhs::template Matrix<3, 2>;
	using MatRhs = typename TestType::Rhs::template Matrix<2, 3>;

	const MatLhs lhs = { 1, 4, 2, 5, 3, 6 };
	const MatRhs same = { 1, 2, 3, 4, 5, 6 };
	const MatRhs different = { 1, 3, 3, 4, 5, 6 };

	REQUIRE(lhs == same);
	REQUIRE(lhs != different);
}