// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../ApplyTransform.hpp"
#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Transforms/TranslationBuilder.hpp>

#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;
using namespace std::complex_literals;


TEMPLATE_LIST_TEST_CASE("Transform: Translation", "[Transforms]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using TestMat = typename TestType::template Matrix<2, 2>;
	using TestVector = Vector<scalar_type_t<TestMat>, 2, is_packed_v<TestMat>>;
	const TestVector testVector = { 2, 2 };
	const TestVector expectedVector = { 5, 6 };

	SECTION("Homogeneous") {
		using M33 = typename TestType::template Matrix<3, 3>;

		if constexpr (order_v<M33> == eMatrixOrder::FOLLOW_VECTOR) {
			const M33 value = Translation(3, 4);
			const M33 expected = {
				1, 0, 0,
				0, 1, 0,
				3, 4, 1
			};
			REQUIRE(value == expected);
			REQUIRE(ApplyTransform(value, testVector) == expectedVector);
		}
		else {
			const M33 value = Translation(3, 4);
			const M33 expected = {
				1, 0, 3,
				0, 1, 4,
				0, 0, 1
			};
			REQUIRE(value == expected);
			REQUIRE(ApplyTransform(value, testVector) == expectedVector);
		}
	}
	SECTION("Slim -- follow") {
		using M32 = typename TestType::template Matrix<3, 2>;

		if constexpr (order_v<M32> == eMatrixOrder::FOLLOW_VECTOR) {
			const M32 value = Translation(3, 4);
			const M32 expected = {
				1, 0,
				0, 1,
				3, 4
			};
			REQUIRE(value == expected);
			REQUIRE(ApplyTransform(value, testVector) == expectedVector);
		}
	}
	SECTION("Slim -- precede") {
		using M23 = typename TestType::template Matrix<2, 3>;

		if constexpr (order_v<M23> == eMatrixOrder::PRECEDE_VECTOR) {
			const M23 value = Translation(3, 4);
			const M23 expected = {
				1, 0, 3,
				0, 1, 4
			};
			REQUIRE(value == expected);
			REQUIRE(ApplyTransform(value, testVector) == expectedVector);
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Translation APIs", "[Transforms]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;

	const Mat mVec = Translation(Vector(1, 2));
	const Mat mArr = Translation(std::array{ 1, 2 });
	const Mat mItems = Translation(1, 2);

	REQUIRE(mVec == mArr);
	REQUIRE(mItems == mArr);
}