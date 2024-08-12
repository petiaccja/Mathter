// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../ApplyTransform.hpp"
#include "../Approx.hpp"
#include "../Cases.hpp"
#include "../Rotation.hpp"

#include <Mathter/Transforms/ViewBuilder.hpp>
#include <Mathter/Vector/Concat.hpp>

#include <array>
#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;
using namespace std::complex_literals;


template <class Mat, class Builder, class Vec, size_t N>
void TestCaseLookAtMatrix(const Builder& builder, const std::array<Vec, N>& testPoints, const std::array<Vec, N>& expectedPoints) {
	const Mat m = builder;

	for (size_t i = 0; i < N; ++i) {
		const auto result = ApplyTransform(m, testPoints[i]) | scalar_type_t<Vec>(1);
		const auto& expected = expectedPoints[i];
		REQUIRE(result == test_util::Approx(expected | scalar_type_t<Vec>(1), 1e-6f));
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: View 3D / generic", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using ExampleMat = typename TestType::template Matrix<1, 1>;
	using Scalar = scalar_type_t<ExampleMat>;
	using Vec = Vector<Scalar, 3, is_packed_v<ExampleMat>>;

	const auto map = [](const Vec& v) {
		const auto rotationAxis = Normalize(Vec(5, 2, 5));
		const auto rotationAngle = Scalar(1.89463572574);
		return Rotate(v, rotationAxis, rotationAngle);
	};

	const auto eye = map(Vec(0, 10, 0));
	const auto target = map(Vec(0, 3, 0));
	const auto up = map(Vec(0, 0.5f, 1));


	const std::array testPoints = {
		map(Vec(0, 10, 0)),
		map(Vec(0, 4, 1)),
		map(Vec(1, 4, 0)),
	};

	SECTION("API") {
		const bool zForward = true;
		const bool flipX = false;
		const bool flipY = false;
		const auto builder = LookAt(eye, target, up, zForward, flipX, flipY);

		const std::array expectedPoints = {
			Vec(0, 0, 0),
			Vec(0, 1, 6),
			Vec(-1, 0, 6),
		};

		SECTION("4x4") {
			using Mat = typename TestType::template Matrix<4, 4>;
			TestCaseLookAtMatrix<Mat>(builder, testPoints, expectedPoints);
		}
		SECTION("3x4") {
			using Mat = typename TestType::template Matrix<3, 4>;
			if constexpr (order_v<Mat> == eMatrixOrder::PRECEDE_VECTOR) {
				TestCaseLookAtMatrix<Mat>(builder, testPoints, expectedPoints);
			}
		}
		SECTION("4x3") {
			using Mat = typename TestType::template Matrix<4, 3>;
			if constexpr (order_v<Mat> == eMatrixOrder::FOLLOW_VECTOR) {
				TestCaseLookAtMatrix<Mat>(builder, testPoints, expectedPoints);
			}
		}
	}
	SECTION("Z backwards") {
		const bool zForward = false;
		const bool flipX = false;
		const bool flipY = false;
		const auto builder = LookAt(eye, target, up, zForward, flipX, flipY);

		const std::array expectedPoints = {
			Vec(0, 0, 0),
			Vec(0, 1, -6),
			Vec(-1, 0, -6),
		};

		using Mat = typename TestType::template Matrix<4, 4>;
		TestCaseLookAtMatrix<Mat>(builder, testPoints, expectedPoints);
	}
	SECTION("Flip X") {
		const bool zForward = true;
		const bool flipX = true;
		const bool flipY = false;
		const auto builder = LookAt(eye, target, up, zForward, flipX, flipY);

		const std::array expectedPoints = {
			Vec(0, 0, 0),
			Vec(0, 1, 6),
			Vec(1, 0, 6),
		};

		using Mat = typename TestType::template Matrix<4, 4>;
		TestCaseLookAtMatrix<Mat>(builder, testPoints, expectedPoints);
	}
	SECTION("Flip Y") {
		const bool zForward = true;
		const bool flipX = false;
		const bool flipY = true;
		const auto builder = LookAt(eye, target, up, zForward, flipX, flipY);

		const std::array expectedPoints = {
			Vec(0, 0, 0),
			Vec(0, -1, 6),
			Vec(-1, 0, 6),
		};

		using Mat = typename TestType::template Matrix<4, 4>;
		TestCaseLookAtMatrix<Mat>(builder, testPoints, expectedPoints);
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: View 2D", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using ExampleMat = typename TestType::template Matrix<1, 1>;
	using Scalar = scalar_type_t<ExampleMat>;
	using Vec = Vector<Scalar, 2, is_packed_v<ExampleMat>>;

	const auto map = [](const Vec& v) {
		const Scalar angle(1.89463572574);
		const auto [c, s] = std::tuple(std::cos(angle), std::sin(angle));
		return Vec(c * v.x - s * v.y, s * v.x + c * v.y);
	};

	const auto eye = map(Vec(0, 10));
	const auto target = map(Vec(0, 3));


	const std::array testPoints = {
		map(Vec(0, 10)),
		map(Vec(1, 4)),
	};

	SECTION("Z forward") {
		const bool zForward = true;
		const bool flipX = false;
		const auto builder = LookAt(eye, target, zForward, flipX);

		const std::array expectedPoints = {
			Vec(0, 0),
			Vec(-1, 6),
		};

		using Mat = typename TestType::template Matrix<3, 3>;
		TestCaseLookAtMatrix<Mat>(builder, testPoints, expectedPoints);
	}
	SECTION("Z backwards") {
		const bool zForward = false;
		const bool flipX = false;
		const auto builder = LookAt(eye, target, zForward, flipX);

		const std::array expectedPoints = {
			Vec(0, 0),
			Vec(-1, -6),
		};

		using Mat = typename TestType::template Matrix<3, 3>;
		TestCaseLookAtMatrix<Mat>(builder, testPoints, expectedPoints);
	}
	SECTION("Flip X") {
		const bool zForward = true;
		const bool flipX = true;
		const auto builder = LookAt(eye, target, zForward, flipX);

		const std::array expectedPoints = {
			Vec(0, 0),
			Vec(1, 6),
		};

		using Mat = typename TestType::template Matrix<3, 3>;
		TestCaseLookAtMatrix<Mat>(builder, testPoints, expectedPoints);
	}
}