// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../ApplyTransform.hpp"
#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Transforms/Rotation2DBuilder.hpp>

#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;
using namespace std::complex_literals;


TEMPLATE_LIST_TEST_CASE("Transform: Rotation 2D", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using TestMat = typename TestType::template Matrix<2, 2>;
	using Scalar = scalar_type_t<TestMat>;
	using TestVector = Vector<Scalar, 2, is_packed_v<TestMat>>;
	const TestVector testVector = { 3, 4 };
	const TestVector expectedVector = { 0, 5 };
	const Scalar angle = 0.64350110879328;

	SECTION("Full") {
		using M22 = typename TestType::template Matrix<2, 2>;

		const M22 value = Rotation(angle);
		const bool v = ApplyTransform(value, testVector) == test_util::Approx(expectedVector, 2e-6f);
		REQUIRE(ApplyTransform(value, testVector) == test_util::Approx(expectedVector, 2e-6f));
	}
	SECTION("Homogeneous") {
		using M33 = typename TestType::template Matrix<3, 3>;

		const M33 value = Rotation(angle);
		REQUIRE(ApplyTransform(value, testVector) == test_util::Approx(expectedVector, 2e-6f));
	}
	SECTION("With translation -- follow") {
		using M32 = typename TestType::template Matrix<3, 2>;

		if constexpr (order_v<M32> == eMatrixOrder::FOLLOW_VECTOR) {
			const M32 value = Rotation(angle);
			REQUIRE(ApplyTransform(value, testVector) == test_util::Approx(expectedVector, 2e-6f));
		}
	}
	SECTION("With translation -- precede") {
		using M23 = typename TestType::template Matrix<2, 3>;

		if constexpr (order_v<M23> == eMatrixOrder::PRECEDE_VECTOR) {
			const M23 value = Rotation(angle);
			REQUIRE(ApplyTransform(value, testVector) == test_util::Approx(expectedVector, 2e-6f));
		}
	}
}
