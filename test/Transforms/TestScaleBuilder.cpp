// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../ApplyTransform.hpp"
#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Transforms/ScaleBuilder.hpp>

#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;
using namespace std::complex_literals;


TEMPLATE_LIST_TEST_CASE("Transform: Scale", "[Transforms]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using TestMat = typename TestType::template Matrix<2, 2>;
	using TestVector = Vector<scalar_type_t<TestMat>, 2, is_packed_v<TestMat>>;
	const TestVector testVector = { 2, 2 };
	const TestVector expectedVector = { 6, 8 };

	SECTION("Full") {
		using M22 = typename TestType::template Matrix<2, 2>;

		const M22 value = Scale(3, 4);
		const M22 expected = {
			3, 0,
			0, 4
		};
		REQUIRE(value == expected);
		REQUIRE(ApplyTransform(value, testVector) == expectedVector);
	}
	SECTION("Homogeneous") {
		using M33 = typename TestType::template Matrix<3, 3>;

		const M33 value = Scale(3, 4);
		const M33 expected = {
			3, 0, 0,
			0, 4, 0,
			0, 0, 1
		};
		REQUIRE(value == expected);
		REQUIRE(ApplyTransform(value, testVector) == expectedVector);
	}
	SECTION("With translation -- follow") {
		using M32 = typename TestType::template Matrix<3, 2>;

		if constexpr (order_v<M32> == eMatrixOrder::FOLLOW_VECTOR) {
			const M32 value = Scale(3, 4);
			const M32 expected = {
				3, 0,
				0, 4,
				0, 0
			};
			REQUIRE(value == expected);
			REQUIRE(ApplyTransform(value, testVector) == expectedVector);
		}
	}
	SECTION("With translation -- precede") {
		using M23 = typename TestType::template Matrix<2, 3>;

		if constexpr (order_v<M23> == eMatrixOrder::PRECEDE_VECTOR) {
			const M23 value = Scale(3, 4);
			const M23 expected = {
				3, 0, 0,
				0, 4, 0
			};
			REQUIRE(value == expected);
			REQUIRE(ApplyTransform(value, testVector) == expectedVector);
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Scale APIs", "[Transforms]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<2, 2>;

	const Mat mVec = Scale(Vector(1, 2));
	const Mat mArr = Scale(std::array{ 1, 2 });
	const Mat mItems = Scale(1, 2);

	REQUIRE(mVec == mArr);
	REQUIRE(mItems == mArr);
}