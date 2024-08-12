// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Geometry/BezierCurve.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>

using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Bezier: construct from array", "[Bezier]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 2>;

	const std::array controlPoints = {
		Vec(1, 2),
		Vec(3, 4),
		Vec(5, 6),
	};

	const BezierCurve curve(controlPoints);

	REQUIRE(curve.controlPoints[0] == controlPoints[0]);
	REQUIRE(curve.controlPoints[1] == controlPoints[1]);
	REQUIRE(curve.controlPoints[2] == controlPoints[2]);
}


TEMPLATE_LIST_TEST_CASE("Bezier: construct from points", "[Bezier]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 2>;

	const std::array controlPoints = {
		Vec(1, 2),
		Vec(3, 4),
		Vec(5, 6),
	};

	const BezierCurve curve(controlPoints[0], controlPoints[1], controlPoints[2]);

	REQUIRE(curve.controlPoints[0] == controlPoints[0]);
	REQUIRE(curve.controlPoints[1] == controlPoints[1]);
	REQUIRE(curve.controlPoints[2] == controlPoints[2]);
}


TEMPLATE_LIST_TEST_CASE("Bezier: covnerting constructor", "[Bezier]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using ScalarLhs = typename TestType::Lhs::Scalar;
	using ScalarRhs = typename TestType::Rhs::Scalar;
	using VecLhs = Vector<ScalarLhs, 2>;
	using BezierLhs = BezierCurve<ScalarLhs, 2, 2>;
	using BezierRhs = BezierCurve<ScalarRhs, 2, 2>;

	const std::array controlPoints = {
		VecLhs(1, 2),
		VecLhs(3, 4),
		VecLhs(5, 6),
	};

	const auto original = BezierLhs(controlPoints);
	const auto converted = BezierRhs(original);

	REQUIRE(original.controlPoints[0] == test_util::Approx(converted.controlPoints[0]));
	REQUIRE(original.controlPoints[1] == test_util::Approx(converted.controlPoints[1]));
	REQUIRE(original.controlPoints[2] == test_util::Approx(converted.controlPoints[2]));
}


TEMPLATE_LIST_TEST_CASE("Bezier: interpolation", "[Bezier]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 2>;

	const std::array controlPoints = {
		Vec(0, 0),
		Vec(1, 1),
		Vec(3, 0),
	};

	const BezierCurve curve(controlPoints[0], controlPoints[1], controlPoints[2]);

	const auto t25 = curve(Scalar(0.25));
	const auto t50 = curve(Scalar(0.5));
	const auto t75 = curve(Scalar(0.75));

	REQUIRE(t25 == test_util::Approx(Vec(0.5625f, 0.375f)));
	REQUIRE(t50 == test_util::Approx(Vec(1.25f, 0.5f)));
	REQUIRE(t75 == test_util::Approx(Vec(2.0625f, 0.375f)));
}