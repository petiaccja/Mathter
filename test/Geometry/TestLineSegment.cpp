// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Geometry/LineSegment.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>

using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("LineSegment: Construct from base & dir & length", "[LineSegment]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 3>;

	const Vec point1 = { 1, 2, 3 };
	const auto direction = Normalize(Vec(4, 2, 6));
	const Vec point2 = point1 + direction * Scalar(2.5);

	const LineSegment<Scalar, 3> lineSegment(point1, direction, Scalar(2.5));

	REQUIRE(lineSegment.Length() == Catch::Approx(2.5));
	REQUIRE(lineSegment.Direction() == test_util::Approx(direction));
	REQUIRE(lineSegment.Start() == test_util::Approx(point1));
	REQUIRE(lineSegment.End() == test_util::Approx(point2));

	const auto line = lineSegment.Line();
	REQUIRE(line.base == test_util::Approx(point1));
	REQUIRE(line.direction == test_util::Approx(direction));
}


TEMPLATE_LIST_TEST_CASE("LineSegment: Construct from endpoints", "[LineSegment]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 3>;

	const Vec point1 = { 1, 2, 3 };
	const auto direction = Normalize(Vec(4, 2, 6));
	const Vec point2 = point1 + direction * Scalar(2.5);

	const LineSegment<Scalar, 3> lineSegment(point1, point2);

	REQUIRE(lineSegment.Start() == test_util::Approx(point1));
	REQUIRE(lineSegment.End() == test_util::Approx(point2));
}


TEMPLATE_LIST_TEST_CASE("LineSegment: Interpolation", "[LineSegment]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 3>;

	const Vec point1 = { 1, 2, 3 };
	const auto direction = Normalize(Vec(2, 3, 6));
	const Vec point2 = point1 + direction * Scalar(21);

	const LineSegment<Scalar, 3> lineSegment(point1, point2);

	REQUIRE(lineSegment.Interpol(Scalar(1.0 / 3.0)) == test_util::Approx(Vec(3, 5, 9)));
	REQUIRE(lineSegment.End() == test_util::Approx(point2));
}