// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../ApplyTransform.hpp"
#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Transforms/OrthographicBuilder.hpp>

#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;
using namespace std::complex_literals;


template <class Mat, class Builder, class Vec>
void TestCaseOrthographicMatrix(const Builder& builder, const Vec& testPoint, const Vec& expectedPoint) {
	const Mat m = builder;
	const auto result = ApplyTransform(m, testPoint);
	REQUIRE(result == test_util::Approx(expectedPoint, 2e-6f));
}


TEMPLATE_LIST_TEST_CASE("Transform: Orthographic / behaviour", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using ExampleMat = typename TestType::template Matrix<1, 1>;
	using Mat = typename TestType::template Matrix<4, 4>;
	using Scalar = scalar_type_t<ExampleMat>;
	using Vec = Vector<Scalar, 3, is_packed_v<ExampleMat>>;

	const Vec minBounds = { -1, -2, -3 };
	const Vec maxBounds = { 3, 2, 5 };

	const std::array testPoints = {
		minBounds,
		maxBounds,
		Vec(0, 1, -1),
	};


	SECTION("Regular Z-planes") {
		const Scalar projNearPlane(1);
		const Scalar projFarPlane(2);

		const auto builder = Orthographic(minBounds, maxBounds, projNearPlane, projFarPlane);
		const std::array expectedPoints = {
			Vec(-1, -1, 1),
			Vec(1, 1, 2),
			Vec(-0.5f, 0.5f, 1.25f),
		};

		for (size_t i = 0; i < testPoints.size(); ++i) {
			TestCaseOrthographicMatrix<Mat>(builder, testPoints[i], expectedPoints[i]);
		}
	}
	SECTION("Inverted Z-planes") {
		const Scalar projNearPlane(2);
		const Scalar projFarPlane(1);

		const auto builder = Orthographic(minBounds, maxBounds, projNearPlane, projFarPlane);
		const std::array expectedPoints = {
			Vec(-1, -1, 2),
			Vec(1, 1, 1),
			Vec(-0.5f, 0.5f, 1.75f),
		};

		for (size_t i = 0; i < testPoints.size(); ++i) {
			TestCaseOrthographicMatrix<Mat>(builder, testPoints[i], expectedPoints[i]);
		}
	}
	SECTION("Negative regular Z-planes") {
		const Scalar projNearPlane(-1);
		const Scalar projFarPlane(-2);

		const auto builder = Orthographic(minBounds, maxBounds, projNearPlane, projFarPlane);
		const std::array expectedPoints = {
			Vec(-1, -1, -1),
			Vec(1, 1, -2),
			Vec(-0.5f, 0.5f, -1.25f),
		};

		for (size_t i = 0; i < testPoints.size(); ++i) {
			TestCaseOrthographicMatrix<Mat>(builder, testPoints[i], expectedPoints[i]);
		}
	}
	SECTION("Negative inverted Z-planes") {
		const Scalar projNearPlane(-2);
		const Scalar projFarPlane(-1);

		const auto builder = Orthographic(minBounds, maxBounds, projNearPlane, projFarPlane);
		const std::array expectedPoints = {
			Vec(-1, -1, -2),
			Vec(1, 1, -1),
			Vec(-0.5f, 0.5f, -1.75f),
		};

		for (size_t i = 0; i < testPoints.size(); ++i) {
			TestCaseOrthographicMatrix<Mat>(builder, testPoints[i], expectedPoints[i]);
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Orthographic / API", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using ExampleMat = typename TestType::template Matrix<1, 1>;
	using Scalar = scalar_type_t<ExampleMat>;
	using Vec = Vector<Scalar, 3, is_packed_v<ExampleMat>>;

	const Vec minBounds = { -1, -2, -3 };
	const Vec maxBounds = { 3, 2, 5 };
	const Scalar projNearPlane(0);
	const Scalar projFarPlane(1);

	const std::array testPoints = {
		minBounds,
		maxBounds,
		Vec(0, 1, -1),
	};
	const std::array expectedPoints = {
		Vec(-1, -1, 0),
		Vec(1, 1, 1),
		Vec(-0.5f, 0.5f, 0.25f),
	};

	const auto builder = Orthographic(minBounds, maxBounds, projNearPlane, projFarPlane);

	SECTION("Homogeneous") {
		using Mat = typename TestType::template Matrix<4, 4>;
		for (size_t i = 0; i < testPoints.size(); ++i) {
			TestCaseOrthographicMatrix<Mat>(builder, testPoints[i], expectedPoints[i]);
		}
	}
	SECTION("Slim -- follow") {
		using Mat = typename TestType::template Matrix<4, 3>;

		if constexpr (order_v<Mat> == eMatrixOrder::FOLLOW_VECTOR) {
			for (size_t i = 0; i < testPoints.size(); ++i) {
				TestCaseOrthographicMatrix<Mat>(builder, testPoints[i], expectedPoints[i]);
			}
		}
	}
	SECTION("Slim -- precede") {
		using Mat = typename TestType::template Matrix<3, 4>;

		if constexpr (order_v<Mat> == eMatrixOrder::PRECEDE_VECTOR) {
			for (size_t i = 0; i < testPoints.size(); ++i) {
				TestCaseOrthographicMatrix<Mat>(builder, testPoints[i], expectedPoints[i]);
			}
		}
	}
}