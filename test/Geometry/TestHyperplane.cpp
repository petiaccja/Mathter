// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Geometry/Hyperplane.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>

using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Hyperplane: construct from base & normal", "[Hyperplane]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 3>;

	const Vec base = { 1, 2, 3 };
	const auto normal = Normalize(Vec(4, 2, 6));
	const auto offset1 = Normalize(Cross(normal, Vec(1, 0, 0)));
	const auto offset2 = Normalize(Cross(normal, offset1));

	const Hyperplane<Scalar, 3> plane(base, normal);

	REQUIRE(plane.Normal() == test_util::Approx(normal));
	REQUIRE(std::abs(plane.Distance(base)) < 1e-6f);
	REQUIRE(std::abs(plane.Distance(plane.Base())) < 1e-6f);
	REQUIRE(std::abs(plane.Distance(base + offset1)) < 1e-6f);
	REQUIRE(std::abs(plane.Distance(base + offset2)) < 1e-6f);
}


TEMPLATE_LIST_TEST_CASE("Hyperplane: construct from normal & scalar", "[Hyperplane]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 3>;

	const Vec base = { 1, 2, 3 };
	const auto normal = Normalize(Vec(4, 2, 6));
	const auto offset1 = Normalize(Cross(normal, Vec(1, 0, 0)));
	const auto offset2 = Normalize(Cross(normal, offset1));

	const Hyperplane<Scalar, 3> plane(normal, Dot(normal, base));

	REQUIRE(plane.Normal() == test_util::Approx(normal));
	REQUIRE(std::abs(plane.Distance(base)) < 1e-6f);
	REQUIRE(std::abs(plane.Distance(plane.Base())) < 1e-6f);
	REQUIRE(std::abs(plane.Distance(base + offset1)) < 1e-6f);
	REQUIRE(std::abs(plane.Distance(base + offset2)) < 1e-6f);
}


TEMPLATE_LIST_TEST_CASE("Hyperplane: convert to line", "[Hyperplane]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 2>;

	const Vec base = { 1, 2 };
	const auto normal = Normalize(Vec(4, 2));

	const Hyperplane<Scalar, 2> plane(base, normal);
	const Line<Scalar, 2> line(plane);

	REQUIRE(std::abs(plane.Distance(line.PointAt(Scalar(1)))) < 1e-6f);
	REQUIRE(std::abs(plane.Distance(line.PointAt(Scalar(2)))) < 1e-6f);
	REQUIRE(std::abs(plane.Distance(line.PointAt(Scalar(3)))) < 1e-6f);
}


TEMPLATE_LIST_TEST_CASE("Hyperplane: construct from line", "[Hyperplane]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	using Vec = Vector<Scalar, 2>;

	const Vec base = { 1, 2 };
	const auto direction = Normalize(Vec(4, 2));

	const Line<Scalar, 2> line(base, direction);
	const Hyperplane<Scalar, 2> plane(line);

	REQUIRE(std::abs(plane.Distance(line.PointAt(Scalar(1)))) < 1e-6f);
	REQUIRE(std::abs(plane.Distance(line.PointAt(Scalar(2)))) < 1e-6f);
	REQUIRE(std::abs(plane.Distance(line.PointAt(Scalar(3)))) < 1e-6f);
}