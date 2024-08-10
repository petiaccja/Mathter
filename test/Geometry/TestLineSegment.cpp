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


TEMPLATE_LIST_TEST_CASE("LineSegment: converting ctor", "[LineSegment]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using ScalarLhs = typename TestType::Lhs::Scalar;
	using ScalarRhs = typename TestType::Rhs::Scalar;
	using VecLhs = Vector<ScalarLhs, 3>;
	using LineLhs = LineSegment<ScalarLhs, 3>;
	using LineRhs = LineSegment<ScalarRhs, 3>;

	const auto point1 = VecLhs(1, 2, 3);
	const auto point2 = VecLhs(4, 2, 6);

	const auto original = LineLhs(point1, point2);
	const auto converted = LineRhs(original);

	REQUIRE(original.point1 == test_util::Approx(converted.point1));
	REQUIRE(original.point2 == test_util::Approx(converted.point2));
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


TEMPLATE_LIST_TEST_CASE("LineSegment: InterpolOf", "[LineSegment]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 3>;

	const Vec point1 = { 1, 2, 3 };
	const auto direction = Normalize(Vec(2, 3, 6));
	const Vec point2 = point1 + direction * Scalar(21);
	const auto perpendicular = Normalize(Cross(direction, Vec(0, 0, 1)));

	const LineSegment<Scalar, 3> lineSegment(point1, point2);

	REQUIRE(lineSegment.InterpolOf(point1 + direction * 6) == Catch::Approx(6.0 / 21.0));
	REQUIRE(lineSegment.InterpolOf(point1 + direction * 6 + perpendicular * 2) == Catch::Approx(6.0 / 21.0));
}