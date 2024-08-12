// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "Approx.hpp"
#include "Cases.hpp"

#include <Mathter/Utility.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Utility: Deg2Rad", "[Utility]",
						decltype(ScalarCaseList<ScalarsFloatingAndComplex>{})) {
	using Scalar = typename TestType::Scalar;
	const Scalar degree(60);
	const Scalar expected(1.047197551196597631);
	const auto result = Deg2Rad(degree);
	REQUIRE(std::real(result) == Catch::Approx(std::real(expected)));
}


TEMPLATE_LIST_TEST_CASE("Utility: Rag2Deg", "[Utility]",
						decltype(ScalarCaseList<ScalarsFloatingAndComplex>{})) {
	using Scalar = typename TestType::Scalar;
	const Scalar radians(1.047197551196597631);
	const Scalar expected(60);
	const auto result = Rad2Deg(radians);
	REQUIRE(std::real(result) == Catch::Approx(std::real(expected)));
}


TEMPLATE_LIST_TEST_CASE("Utility: Clamp / scalar", "[Utility]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	const Scalar lower(-3);
	const Scalar upper(3);
	SECTION("Below") {
		REQUIRE(Clamp(Scalar(-4), lower, upper) == lower);
	}
	SECTION("Above") {
		REQUIRE(Clamp(Scalar(4), lower, upper) == upper);
	}
	SECTION("Between") {
		REQUIRE(Clamp(Scalar(2), lower, upper) == Scalar(2));
	}
}


TEMPLATE_LIST_TEST_CASE("Utility: Saturate / scalar", "[Utility]",
						decltype(ScalarCaseList<ScalarsFloating>{})) {
	using Scalar = typename TestType::Scalar;
	SECTION("Below") {
		REQUIRE(Saturate(Scalar(-1)) == Scalar(0));
	}
	SECTION("Above") {
		REQUIRE(Saturate(Scalar(2)) == Scalar(1));
	}
	SECTION("Between") {
		REQUIRE(Saturate(Scalar(0.5)) == Scalar(0.5));
	}
}


TEMPLATE_LIST_TEST_CASE("Utility: Clamp / vector", "[Utility]",
						decltype(VectorCaseList<ScalarsFloating, PackingsAll>{})) {
	using Vec = typename TestType::template Vector<3>;
	using Scalar = scalar_type_t<Vec>;
	const Scalar lower(-3);
	const Scalar upper(3);
	const Vec value = { -4, 4, 2 };
	const auto result = Clamp(value, lower, upper);
	REQUIRE(result[0] == lower);
	REQUIRE(result[1] == upper);
	REQUIRE(result[2] == value[2]);
}


TEMPLATE_LIST_TEST_CASE("Utility: Saturate / vector", "[Utility]",
						decltype(VectorCaseList<ScalarsFloating, PackingsAll>{})) {
	using Vec = typename TestType::template Vector<3>;
	using Scalar = scalar_type_t<Vec>;
	const Vec value = { -1, 2, 0.5f };
	const auto result = Saturate(value);
	REQUIRE(result[0] == Scalar(0));
	REQUIRE(result[1] == Scalar(1));
	REQUIRE(result[2] == value[2]);
}