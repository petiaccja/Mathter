// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Transforms/IdentityBuilder.hpp>

#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;
using namespace std::complex_literals;


TEMPLATE_LIST_TEST_CASE("Transform: Identity - Matrix", "[Transforms]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using S = scalar_type_t<typename TestType::template Matrix<1, 1>>;

	SECTION("Rectangular") {
		using M23 = typename TestType::template Matrix<2, 3>;
		const M23 value = Identity();

		REQUIRE(value(0, 0) == static_cast<S>(1));
		REQUIRE(value(0, 1) == static_cast<S>(0));
		REQUIRE(value(0, 2) == static_cast<S>(0));
		REQUIRE(value(1, 0) == static_cast<S>(0));
		REQUIRE(value(1, 1) == static_cast<S>(1));
		REQUIRE(value(1, 2) == static_cast<S>(0));
	}
	SECTION("Square") {
		using M22 = typename TestType::template Matrix<2, 2>;
		const M22 value = Identity();

		REQUIRE(value(0, 0) == static_cast<S>(1));
		REQUIRE(value(0, 1) == static_cast<S>(0));
		REQUIRE(value(1, 0) == static_cast<S>(0));
		REQUIRE(value(1, 1) == static_cast<S>(1));
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Identity - Quaternion", "[Transforms]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	const Quat value = Identity();
	REQUIRE(value.s == 1);
	REQUIRE(value.i == 0);
	REQUIRE(value.j == 0);
	REQUIRE(value.k == 0);
}