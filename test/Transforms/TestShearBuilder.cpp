// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../ApplyTransform.hpp"
#include "../Approx.hpp"
#include "../Cases.hpp"
#include "../Rotation.hpp"

#include <Mathter/Transforms/ShearBuilder.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


template <class Mat, class Builder, class Vec>
void TestCaseShearMatrix(const Builder& builder, const Vec& testPoint, const Vec& expectedPoint) {
	const Mat m = builder;
	const auto result = ApplyTransform(m, testPoint);

	REQUIRE(result == test_util::Approx(expectedPoint));
	if constexpr (row_count_v<Mat> == column_count_v<Mat>) {
		const auto det = Determinant(m);
		REQUIRE(std::real(det) == Catch::Approx(1));
		REQUIRE(std::imag(det) == Catch::Approx(0));
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Shear -- aligned / behaviour", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M33 = typename TestType::template Matrix<3, 3>;
	using Scalar = scalar_type_t<M33>;
	using Vec = Vector<Scalar, 3, is_packed_v<M33>>;

	const Vec testPoint = { 1, 2, 3 };

	SECTION("Positive slope") {
		const Scalar slope = 0.5f;

		const auto builder = Shear(slope, 1, 2);
		const Vec expected = { 1, 3.5f, 3 };

		TestCaseShearMatrix<M33>(builder, testPoint, expected);
	}
	SECTION("Negative slope") {
		const Scalar slope = -0.5f;

		const auto builder = Shear(slope, 1, 2);
		const Vec expected = { 1, 0.5f, 3 };

		TestCaseShearMatrix<M33>(builder, testPoint, expected);
	}
	SECTION("Inverted axes") {
		const Scalar slope = 0.5f;

		const auto builder = Shear(slope, 2, 1);
		const Vec expected = { 1, 2, 4 };

		TestCaseShearMatrix<M33>(builder, testPoint, expected);
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Shear -- aligned / API", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M33 = typename TestType::template Matrix<3, 3>;
	using Scalar = scalar_type_t<M33>;
	using Vec = Vector<Scalar, 3, is_packed_v<M33>>;

	const Vec testPoint = { 1, 2, 3 };
	const Vec expectedPoint = { 1, 3.5f, 3 };
	const Scalar slope = 0.5f;
	const auto builder = Shear(slope, 1, 2);

	SECTION("3x3") {
		using Mat = typename TestType::template Matrix<3, 3>;
		TestCaseShearMatrix<Mat>(builder, testPoint, expectedPoint);
	}
	SECTION("4x4") {
		using Mat = typename TestType::template Matrix<4, 4>;
		TestCaseShearMatrix<Mat>(builder, testPoint, expectedPoint);
	}
	SECTION("3x4") {
		using Mat = typename TestType::template Matrix<3, 4>;
		if constexpr (order_v<Mat> == eMatrixOrder::PRECEDE_VECTOR) {
			TestCaseShearMatrix<Mat>(builder, testPoint, expectedPoint);
		}
	}
	SECTION("4x3") {
		using Mat = typename TestType::template Matrix<4, 3>;
		if constexpr (order_v<Mat> == eMatrixOrder::FOLLOW_VECTOR) {
			TestCaseShearMatrix<Mat>(builder, testPoint, expectedPoint);
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Shear -- general / behaviour", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M33 = typename TestType::template Matrix<3, 3>;
	using Scalar = scalar_type_t<M33>;
	using Vec = Vector<Scalar, 3, is_packed_v<M33>>;

	const auto rotationAxis = Normalize(Vec(4, 2, 5));
	const auto rotationAngle = Scalar(0.896345786489);

	const auto directionAxis = test_util::Rotate(Vec(1, 0, 0), rotationAxis, rotationAngle);
	const auto modulatorAxis = test_util::Rotate(Vec(0, 0, 1), rotationAxis, rotationAngle);
	const auto testPoint = test_util::Rotate(Vec(1, 2, 3), rotationAxis, rotationAngle);

	SECTION("Positive slope") {
		const Scalar slope = 0.5f;

		const auto builder = Shear(slope, directionAxis, modulatorAxis);
		const auto expected = test_util::Rotate(Vec(2.5f, 2, 3), rotationAxis, rotationAngle);

		TestCaseShearMatrix<M33>(builder, testPoint, expected);
	}
	SECTION("Negative slope") {
		const Scalar slope = -0.5f;

		const auto builder = Shear(slope, directionAxis, modulatorAxis);
		const auto expected = test_util::Rotate(Vec(-0.5f, 2, 3), rotationAxis, rotationAngle);

		TestCaseShearMatrix<M33>(builder, testPoint, expected);
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Shear -- general / API", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M33 = typename TestType::template Matrix<3, 3>;
	using Scalar = scalar_type_t<M33>;
	using Vec = Vector<Scalar, 3, is_packed_v<M33>>;

	const Scalar slope = 0.5f;
	const auto directionAxis = Vec(1, 0, 0);
	const auto modulatorAxis = Vec(0, 0, 1);

	const Vec testPoint = { 1, 2, 3 };
	const Vec expectedPoint = { 2.5f, 2, 3 };

	const auto builder = Shear(slope, directionAxis, modulatorAxis);

	SECTION("3x3") {
		using Mat = typename TestType::template Matrix<3, 3>;
		TestCaseShearMatrix<Mat>(builder, testPoint, expectedPoint);
	}
	SECTION("4x4") {
		using Mat = typename TestType::template Matrix<4, 4>;
		TestCaseShearMatrix<Mat>(builder, testPoint, expectedPoint);
	}
	SECTION("3x4") {
		using Mat = typename TestType::template Matrix<3, 4>;
		if constexpr (order_v<Mat> == eMatrixOrder::PRECEDE_VECTOR) {
			TestCaseShearMatrix<Mat>(builder, testPoint, expectedPoint);
		}
	}
	SECTION("4x3") {
		using Mat = typename TestType::template Matrix<4, 3>;
		if constexpr (order_v<Mat> == eMatrixOrder::FOLLOW_VECTOR) {
			TestCaseShearMatrix<Mat>(builder, testPoint, expectedPoint);
		}
	}
}