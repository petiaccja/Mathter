// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Geometry/Triangle.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>

using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Triangle: Constructor", "[Triangle]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 2>;

	const Vec a(1, 0);
	const Vec b(5, 0);
	const Vec c(1, 3);

	SECTION("Individual points") {
		const Triangle<Scalar, 2> triangle(a, b, c);
		REQUIRE(triangle.corners[0] == a);
		REQUIRE(triangle.corners[1] == b);
		REQUIRE(triangle.corners[2] == c);
	}
	SECTION("Array of points") {
		const Triangle<Scalar, 2> triangle(std::array{ a, b, c });
		REQUIRE(triangle.corners[0] == a);
		REQUIRE(triangle.corners[1] == b);
		REQUIRE(triangle.corners[2] == c);
	}
}


TEMPLATE_LIST_TEST_CASE("Triangle: converting ctor", "[Triangle]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using ScalarLhs = typename TestType::Lhs::Scalar;
	using ScalarRhs = typename TestType::Rhs::Scalar;
	using VecLhs = Vector<ScalarLhs, 2>;
	using TriangleLhs = Triangle<ScalarLhs, 2>;
	using TriangleRhs = Triangle<ScalarRhs, 2>;

	const VecLhs a(1, 0);
	const VecLhs b(5, 0);
	const VecLhs c(1, 3);

	const auto original = TriangleLhs(a, b, c);
	const auto converted = TriangleRhs(original);

	REQUIRE(original.corners[0] == test_util::Approx(converted.corners[0]));
	REQUIRE(original.corners[1] == test_util::Approx(converted.corners[1]));
	REQUIRE(original.corners[2] == test_util::Approx(converted.corners[2]));
}


TEMPLATE_LIST_TEST_CASE("Triangle: area", "[Triangle]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 2>;

	const Vec a(1, 0);
	const Vec b(5, 0);
	const Vec c(1, 3);

	const Triangle triangle(a, b, c);
	REQUIRE(triangle.Area() == Catch::Approx(6));
}


TEMPLATE_LIST_TEST_CASE("Triangle: centroid", "[Triangle]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 2>;

	const Vec a(1, 0);
	const Vec b(5, 0);
	const Vec c(1, 3);

	const Triangle triangle(a, b, c);
	const auto centroid = triangle.Centroid();

	// Distance of the two corners from the line the crosses the third corner and the centroid must be equal.
	const auto lineA = Normalize(Cross(centroid - a));
	const auto lineB = Normalize(Cross(centroid - b));
	const auto lineC = Normalize(Cross(centroid - c));
	REQUIRE(Dot(lineA, centroid - b) == Catch::Approx(Dot(lineA, c - centroid)));
	REQUIRE(Dot(lineB, centroid - a) == Catch::Approx(Dot(lineB, c - centroid)));
	REQUIRE(Dot(lineC, centroid - a) == Catch::Approx(Dot(lineC, b - centroid)));
}