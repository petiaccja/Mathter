// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Geometry/Ray.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>

using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Ray: colinear line", "[Ray]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 3>;

	const Vec base = { 1, 2, 3 };
	const auto direction = Normalize(Vec(4, 2, 6));

	const Ray<Scalar, 3> ray(base, direction);
	const auto line = ray.Line();

	REQUIRE(ray.PointAt(Scalar(1)) == test_util::Approx(line.PointAt(Scalar(1))));
	REQUIRE(ray.PointAt(Scalar(2)) == test_util::Approx(line.PointAt(Scalar(2))));
}


TEMPLATE_LIST_TEST_CASE("Ray: converting ctor", "[Ray]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using ScalarLhs = typename TestType::Lhs::Scalar;
	using ScalarRhs = typename TestType::Rhs::Scalar;
	using VecLhs = Vector<ScalarLhs, 3>;
	using RayLhs = Ray<ScalarLhs, 3>;
	using RayRhs = Ray<ScalarRhs, 3>;

	const VecLhs base = { 1, 2, 3 };
	const auto direction = Normalize(VecLhs(4, 2, 6));

	const auto original = RayLhs(base, direction);
	const auto converted = RayRhs(original);

	REQUIRE(original.base == test_util::Approx(converted.base));
	REQUIRE(original.direction == test_util::Approx(converted.direction));
}