// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"
#include "../Rotation.hpp"

#include <Mathter/Geometry/Intersection.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>

using namespace mathter;
using namespace test_util;


template <class Vec>
auto Misalign(const Vec& v) {
	const auto axis = Normalize(Vec(5, 2, 3));
	const auto angle = scalar_type_t<Vec>(0.345824657);
	return Rotate(v, axis, angle);
}


TEMPLATE_LIST_TEST_CASE("Intersection: hyperplane - line", "[Intersection]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using S1 = typename TestType::Lhs::Scalar;
	using S2 = typename TestType::Rhs::Scalar;
	using VecL = Vector<S1, 3>;
	using VecP = Vector<S2, 3>;

	SECTION("Intersecting") {
		const Line<S1, 3> line(Misalign(VecL(0, 1, 1)), Misalign(VecL(0, 1, 0)));
		const Hyperplane<S2, 3> plane(Misalign(VecP(0, 6, 1)), Misalign(Normalize(VecP(3, 2, 6))));

		const auto intersection = Intersect(line, plane);
		REQUIRE(intersection == Intersect(plane, line)); // Flipped operands.
		REQUIRE(intersection.has_value());
		REQUIRE(intersection.value() == test_util::Approx(Misalign(VecP(0, 6, 1)), 1e-6f));
	}
	SECTION("Almost parallel") {
		const Line<S1, 3> line(VecL(0, 1, 1), Normalize(VecL(0, 1, 0)));
		const Hyperplane<S2, 3> plane(VecP(0, 6, 1), Normalize(VecP(0, 1e-5f, 1)));
		const auto intersection = Intersect(line, plane);
		REQUIRE(intersection == Intersect(plane, line)); // Flipped operands.
		REQUIRE(intersection.has_value());

		// The floating point accuracy on this one is really poor.
		// To my understanding, this is because shifting the `base` of the plane from
		// the original to one that lies on the plane's normal axis shifts the plane
		// along the normal just a tiny bit, as the new base cannot be exactly represented
		// by floats. This tiny shift is magnified when a line is nearly parallel to the plane.
		// Once the original base is recalculated as the scalar in the plane's equation, there
		// is nothing that can be done. The scalar can be represented as two floats using a
		// compensated sum, but I'm not sure if it's worth the hassle. Even with the compensated
		// scalar, the new base cannot be more accurately represented, so intersections may work,
		// but conversions to lines are not improved.
		REQUIRE(intersection.value() == test_util::Approx(VecP(0, 6, 1), 1e-2f));
	}
	SECTION("Parallel") {
		const Line<S1, 3> line(VecL(0, 1, 1), Normalize(VecL(0, 1, 0)));
		const Hyperplane<S2, 3> plane(VecP(0, 6, 1), Normalize(VecP(0, 0, 1)));
		const auto intersection = Intersect(line, plane);
		REQUIRE(intersection == Intersect(plane, line)); // Flipped operands.
		REQUIRE(!intersection.has_value());
	}
}


TEMPLATE_LIST_TEST_CASE("Intersection: hyperplane - line segment", "[Intersection]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using S1 = typename TestType::Lhs::Scalar;
	using S2 = typename TestType::Rhs::Scalar;
	using VecL = Vector<S1, 3>;
	using VecP = Vector<S2, 3>;

	SECTION("Intersecting") {
		const LineSegment<S1, 3> line(VecL(2, 3, -1), VecL(2, 3, 1));
		const Hyperplane<S2, 3> plane(VecP(0, 0, 0), Normalize(VecP(0, 0, 1)));

		const auto intersection = Intersect(line, plane);
		REQUIRE(intersection == Intersect(plane, line)); // Flipped operands.
		REQUIRE(intersection.has_value());
		REQUIRE(intersection.value() == test_util::Approx(VecP(2, 3, 0), 1e-6f));
	}
	SECTION("Short start") {
		const LineSegment<S1, 3> line(VecL(2, 3, 0.001f), VecL(2, 3, 1));
		const Hyperplane<S2, 3> plane(VecP(0, 0, 0), Normalize(VecP(0, 0, 1)));

		const auto intersection = Intersect(line, plane);
		REQUIRE(intersection == Intersect(plane, line)); // Flipped operands.
		REQUIRE(!intersection.has_value());
	}
	SECTION("Short end") {
		const LineSegment<S1, 3> line(VecL(2, 3, -1), VecL(2, 3, -0.001f));
		const Hyperplane<S2, 3> plane(VecP(0, 0, 0), Normalize(VecP(0, 0, 1)));

		const auto intersection = Intersect(line, plane);
		REQUIRE(intersection == Intersect(plane, line)); // Flipped operands.
		REQUIRE(!intersection.has_value());
	}
}


TEMPLATE_LIST_TEST_CASE("Intersection: hyperplane - ray", "[Intersection]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using S1 = typename TestType::Lhs::Scalar;
	using S2 = typename TestType::Rhs::Scalar;
	using VecL = Vector<S1, 3>;
	using VecP = Vector<S2, 3>;

	SECTION("Intersecting") {
		const Ray<S1, 3> ray(VecL(2, 3, -1), Normalize(VecL(0, 0, 1)));
		const Hyperplane<S2, 3> plane(VecP(0, 0, 0), Normalize(VecP(0, 0, 1)));

		const auto intersection = Intersect(ray, plane);
		REQUIRE(intersection == Intersect(plane, ray)); // Flipped operands.
		REQUIRE(intersection.has_value());
		REQUIRE(intersection.value() == test_util::Approx(VecP(2, 3, 0), 1e-6f));
	}
	SECTION("Short start") {
		const LineSegment<S1, 3> ray(VecL(2, 3, 0.001f), VecL(0, 0, 1));
		const Hyperplane<S2, 3> plane(VecP(0, 0, 0), Normalize(VecP(0, 0, 1)));

		const auto intersection = Intersect(ray, plane);
		REQUIRE(intersection == Intersect(plane, ray)); // Flipped operands.
		REQUIRE(!intersection.has_value());
	}
}


TEMPLATE_LIST_TEST_CASE("Intersection: line<2> - line<2>", "[Intersection]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using S1 = typename TestType::Lhs::Scalar;
	using S2 = typename TestType::Rhs::Scalar;
	using VecL = Vector<S1, 2>;
	using VecP = Vector<S2, 2>;

	SECTION("Intersecting") {
		const Line<S1, 2> line1(VecL(0, 0), Normalize(VecL(1, 0)));
		const Line<S2, 2> line2(VecP(0, 1), Normalize(VecP(1, -1)));

		const auto intersection = Intersect(line1, line2);
		REQUIRE(intersection.has_value());
		REQUIRE(intersection.value() == test_util::Approx(VecP(1, 0), 1e-6f));
	}
	SECTION("Parallel") {
		const Line<S1, 2> line1(VecL(0, 0), Normalize(VecL(1, 0)));
		const Line<S2, 2> line2(VecP(0, 1), Normalize(VecP(1, 0)));

		const auto intersection = Intersect(line1, line2);
		REQUIRE(!intersection.has_value());
	}
}


TEMPLATE_LIST_TEST_CASE("Intersection: line<2> - line segment<2>", "[Intersection]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using S1 = typename TestType::Lhs::Scalar;
	using S2 = typename TestType::Rhs::Scalar;
	using VecL = Vector<S1, 2>;
	using VecP = Vector<S2, 2>;

	SECTION("Intersecting") {
		const Line<S1, 2> line1(VecL(0, 0), Normalize(VecL(1, 0)));
		const LineSegment<S2, 2> line2(VecP(0, 1), VecP(2, -1));

		const auto intersection = Intersect(line1, line2);
		REQUIRE(intersection == Intersect(line2, line1)); // Flipped operands.
		REQUIRE(intersection.has_value());
		REQUIRE(intersection.value() == test_util::Approx(VecP(1, 0), 1e-6f));
	}
	SECTION("Short start") {
		const Line<S1, 2> line1(VecL(0, 0), Normalize(VecL(1, 0)));
		const LineSegment<S2, 2> line2(VecP(0, 1), VecP(2, 0.001f));

		const auto intersection = Intersect(line1, line2);
		REQUIRE(intersection == Intersect(line2, line1)); // Flipped operands.
		REQUIRE(!intersection.has_value());
	}
	SECTION("Short end") {
		const Line<S1, 2> line1(VecL(0, 0), Normalize(VecL(1, 0)));
		const LineSegment<S2, 2> line2(VecP(0, -0.001f), VecP(2, -1));

		const auto intersection = Intersect(line1, line2);
		REQUIRE(intersection == Intersect(line2, line1)); // Flipped operands.
		REQUIRE(!intersection.has_value());
	}
}


TEMPLATE_LIST_TEST_CASE("Intersection: line segment<2> - line segment<2>", "[Intersection]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using S1 = typename TestType::Lhs::Scalar;
	using S2 = typename TestType::Rhs::Scalar;
	using VecL = Vector<S1, 2>;
	using VecP = Vector<S2, 2>;

	SECTION("Intersecting") {
		const LineSegment<S1, 2> line1(VecL(-2, 0), VecL(2, 0));
		const LineSegment<S2, 2> line2(VecP(0, 1), VecP(2, -1));

		const auto intersection = Intersect(line1, line2);
		REQUIRE(intersection.has_value());
		REQUIRE(intersection.value() == test_util::Approx(VecP(1, 0), 1e-6f));
	}
	SECTION("Short start") {
		const LineSegment<S1, 2> line1(VecL(-2, 0), VecL(2, 0));
		const LineSegment<S2, 2> line2(VecP(0, 1), VecP(2, 0.001f));

		const auto intersection = Intersect(line1, line2);
		REQUIRE(!intersection.has_value());
	}
	SECTION("Short end") {
		const LineSegment<S1, 2> line1(VecL(-2, 0), VecL(2, 0));
		const LineSegment<S2, 2> line2(VecP(0, -0.001f), VecP(2, -1));

		const auto intersection = Intersect(line1, line2);
		REQUIRE(!intersection.has_value());
	}
	SECTION("Outside") {
		const LineSegment<S1, 2> line1(VecL(-2, 0), VecL(2, 0));
		const LineSegment<S2, 2> line2(VecP(4, 2), VecP(4, 1));

		const auto intersection = Intersect(line1, line2);
		REQUIRE(!intersection.has_value());
	}
}


TEMPLATE_LIST_TEST_CASE("Intersection: ray<3> - triangle<3>", "[Intersection]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using S1 = typename TestType::Lhs::Scalar;
	using S2 = typename TestType::Rhs::Scalar;
	using VecL = Vector<S1, 3>;
	using VecP = Vector<S2, 3>;

	const Triangle<S1, 3> triangle(Misalign(VecL(0, 0, 0)), Misalign(VecL(1, 0, 0)), Misalign(VecL(0, 1, 0)));


	SECTION("Misses side 1") {
		Ray<S2, 3> ray(Misalign(VecP(0.5f, -0.001f, -1.0f)), Misalign(VecP(0, 0, 1)));

		const auto intersection = Intersect(ray, triangle);
		REQUIRE(intersection == Intersect(triangle, ray));
		REQUIRE(!intersection.has_value());
	}
	SECTION("Misses side 2") {
		Ray<S2, 3> ray(Misalign(VecP(0.5001f, 0.5001f, -1.0f)), Misalign(VecP(0, 0, 1)));

		const auto intersection = Intersect(ray, triangle);
		REQUIRE(intersection == Intersect(triangle, ray));
		REQUIRE(!intersection.has_value());
	}
	SECTION("Misses side 3") {
		Ray<S2, 3> ray(Misalign(VecP(-0.001f, 0.5f, -1.0f)), Misalign(VecP(0, 0, 1)));

		const auto intersection = Intersect(ray, triangle);
		REQUIRE(intersection == Intersect(triangle, ray));
		REQUIRE(!intersection.has_value());
	}
	SECTION("Hit") {
		Ray<S2, 3> ray(Misalign(VecP(0.001f, 0.5f, -1.0f)), Misalign(VecP(0, 0, 1)));

		const auto intersection = Intersect(ray, triangle);
		REQUIRE(intersection == Intersect(triangle, ray));
		REQUIRE(intersection.has_value());
		REQUIRE(intersection.value() == test_util::Approx(Misalign(VecP(0.001f, 0.5f, 0.0f)), 1e-6f));
	}
	SECTION("Short") {
		Ray<S2, 3> ray(Misalign(VecP(0.001f, 0.5f, 0.001f)), Misalign(VecP(0, 0, 1)));

		const auto intersection = Intersect(ray, triangle);
		REQUIRE(intersection == Intersect(triangle, ray));
		REQUIRE(!intersection.has_value());
	}
}