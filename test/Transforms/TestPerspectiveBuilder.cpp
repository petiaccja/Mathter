// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../ApplyTransform.hpp"
#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Transforms/PerspectiveBuilder.hpp>

#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;
using namespace std::complex_literals;


TEMPLATE_LIST_TEST_CASE("Transform: Perspective 3D / generic", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using ExampleMat = typename TestType::template Matrix<1, 1>;
	using Mat = typename TestType::template Matrix<4, 4>;
	using Scalar = scalar_type_t<ExampleMat>;
	using Vec = Vector<Scalar, 3, is_packed_v<ExampleMat>>;

	const Scalar fovX(0.92729521800161); // A 1 unit wide piece at 2 units away fits the screen exactly.
	const Scalar aspectRatio(2.0);

	const std::array projPlanes = {
		std::pair{ Scalar(1), Scalar(2) },
		std::pair{ Scalar(2), Scalar(1) },
		std::pair{ Scalar(-1), Scalar(-2) },
		std::pair{ Scalar(-2), Scalar(-1) },
		std::pair{ Scalar(-1), Scalar(1) },
		std::pair{ Scalar(1), Scalar(-1) },
	};

	SECTION("Positive near & far plane") {
		const Scalar nearPlane(0.5);
		const Scalar farPlane(2.0);

		for (const auto& [projNearPlane, projFarPlane] : projPlanes) {
			const std::array<Vec, 3> testPoints = {
				Vec(-0.25f, -0.125f, nearPlane),
				Vec(1.0f, 0.5f, farPlane),
				Vec(0.4f, 0.1f, 1.6f),
			};

			const std::array<Vec, 2> expectedPoints = {
				Vec(-1.0f, -1.0f, projNearPlane),
				Vec(1.0f, 1.0f, projFarPlane),
			};

			const auto builder = Perspective(fovX, aspectRatio, nearPlane, farPlane, projNearPlane, projFarPlane);
			const Mat m = builder;
			REQUIRE(ApplyTransform(m, testPoints[0]) == test_util::Approx(expectedPoints[0], 1e-6f));
			REQUIRE(ApplyTransform(m, testPoints[1]) == test_util::Approx(expectedPoints[1], 1e-6f));
			const auto insidePoint = ApplyTransform(m, testPoints[2]);
			REQUIRE(std::abs(insidePoint[0]) <= Scalar(1.0));
			REQUIRE(std::abs(insidePoint[1]) <= Scalar(1.0));
			REQUIRE(std::min(projNearPlane, projFarPlane) <= insidePoint[2]);
			REQUIRE(insidePoint[2] <= std::max(projNearPlane, projFarPlane));
		}
	}
	SECTION("Negative near & far plane") {
		const Scalar nearPlane(-0.5);
		const Scalar farPlane(-2.0);

		for (const auto& [projNearPlane, projFarPlane] : projPlanes) {
			const std::array testPoints = {
				Vec(-0.25f, -0.125f, nearPlane),
				Vec(1.0f, 0.5f, farPlane),
				Vec(0.4f, 0.1f, -1.6f),
			};

			const std::array expectedPoints = {
				Vec(-1.0f, -1.0f, projNearPlane),
				Vec(1.0f, 1.0f, projFarPlane),
			};

			const auto builder = Perspective(fovX, aspectRatio, nearPlane, farPlane, projNearPlane, projFarPlane);
			const Mat m = builder;
			REQUIRE(ApplyTransform(m, testPoints[0]) == test_util::Approx(expectedPoints[0], 1e-6f));
			REQUIRE(ApplyTransform(m, testPoints[1]) == test_util::Approx(expectedPoints[1], 1e-6f));
			const auto insidePoint = ApplyTransform(m, testPoints[2]);
			REQUIRE(std::abs(insidePoint[0]) <= Scalar(1.0));
			REQUIRE(std::abs(insidePoint[1]) <= Scalar(1.0));
			REQUIRE(std::min(projNearPlane, projFarPlane) <= insidePoint[2]);
			REQUIRE(insidePoint[2] <= std::max(projNearPlane, projFarPlane));
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Perspective 2D", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	// This test does not need to test all combinations exhaustively, the generic function
	// is tested exhaustively via the 3D test.
	// This is just to ensure the 2D interface works as intended.

	using ExampleMat = typename TestType::template Matrix<1, 1>;
	using Mat = typename TestType::template Matrix<3, 3>;
	using Scalar = scalar_type_t<ExampleMat>;
	using Vec = Vector<Scalar, 2, is_packed_v<ExampleMat>>;

	const Scalar fovX(0.92729521800161); // A 1 unit wide piece at 2 units away fits the screen exactly.
	const Scalar nearPlane(0.5);
	const Scalar farPlane(2);
	const Scalar projNearPlane(0);
	const Scalar projFarPlane(1);
	const auto builder = Perspective(fovX, nearPlane, farPlane, projNearPlane, projFarPlane);

	const std::array testPoints = {
		Vec(-0.25f, nearPlane),
		Vec(1.0f, farPlane),
	};

	const std::array expectedPoints = {
		Vec(-1.0f, projNearPlane),
		Vec(1.0f, projFarPlane),
	};

	const Mat m = builder;
	REQUIRE(ApplyTransform(m, testPoints[0]) == test_util::Approx(expectedPoints[0], 1e-6f));
	REQUIRE(ApplyTransform(m, testPoints[1]) == test_util::Approx(expectedPoints[1], 1e-6f));
}