// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Geometry/Line.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>

using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Line: construct from base & dir", "[Line]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 3>;

	const Vec base = { 1, 2, 3 };
	const auto direction = Normalize(Vec(4, 2, 6));

	const Line<Scalar, 3> line(base, direction);

	REQUIRE(line.base == base);
	REQUIRE(line.Base() == base);
	REQUIRE(line.direction == direction);
	REQUIRE(line.Direction() == direction);
}


TEMPLATE_LIST_TEST_CASE("Line: construct from points", "[Line]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 3>;

	const Vec point1 = { 1, 2, 3 };
	const auto direction = Normalize(Vec(4, 2, 6));
	const Vec point2 = point1 + direction * Scalar(2.5);

	const auto line = Line<Scalar, 3>::Through(point1, point2);

	REQUIRE(line.base == test_util::Approx(point1));
	REQUIRE(line.Base() == test_util::Approx(point1));
	REQUIRE(line.direction == test_util::Approx(direction));
	REQUIRE(line.Direction() == test_util::Approx(direction));
}


TEMPLATE_LIST_TEST_CASE("Line: converting ctor", "[Line]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using ScalarLhs = typename TestType::Lhs::Scalar;
	using ScalarRhs = typename TestType::Rhs::Scalar;
	using VecLhs = Vector<ScalarLhs, 3>;
	using LineLhs = Line<ScalarLhs, 3>;
	using LineRhs = Line<ScalarRhs, 3>;

	const VecLhs base = { 1, 2, 3 };
	const auto direction = Normalize(VecLhs(4, 2, 6));

	const auto original = LineLhs(base, direction);
	const auto converted = LineRhs(original);

	REQUIRE(original.base == test_util::Approx(converted.base));
	REQUIRE(original.direction == test_util::Approx(converted.direction));
}


TEMPLATE_LIST_TEST_CASE("Line: PointAt", "[Line]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 3>;

	const Vec base = { 1, 2, 3 };
	const auto direction = Normalize(Vec(2, 3, 6));

	const Line<Scalar, 3> line(base, direction);

	const auto result = line.PointAt(Scalar(7));
	REQUIRE(result == test_util::Approx(Vec(3, 5, 9)));
}