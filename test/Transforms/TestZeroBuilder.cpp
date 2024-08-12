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


TEMPLATE_LIST_TEST_CASE("Transform: Zero - Matrix", "[Transforms]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {

	using M23 = typename TestType::template Matrix<2, 3>;
	using S = scalar_type_t<M23>;

	const M23 value = Zero();

	REQUIRE(value(0, 0) == static_cast<S>(0));
	REQUIRE(value(0, 1) == static_cast<S>(0));
	REQUIRE(value(0, 2) == static_cast<S>(0));
	REQUIRE(value(1, 0) == static_cast<S>(0));
	REQUIRE(value(1, 1) == static_cast<S>(0));
	REQUIRE(value(1, 2) == static_cast<S>(0));
}


TEMPLATE_LIST_TEST_CASE("Transform: Zero - Vector", "[Transforms]",
						decltype(VectorCaseList<ScalarsAll, PackingsAll>{})) {
	using Vec = typename TestType::template Vector<3>;
	using S = scalar_type_t<Vec>;

	const Vec value = Zero();
	REQUIRE(value[0] == static_cast<S>(0));
	REQUIRE(value[1] == static_cast<S>(0));
	REQUIRE(value[2] == static_cast<S>(0));
}


TEMPLATE_LIST_TEST_CASE("Transform: Zero - Quaternion", "[Transforms]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	const Quat value = Zero();
	REQUIRE(value.s == 0);
	REQUIRE(value.i == 0);
	REQUIRE(value.j == 0);
	REQUIRE(value.k == 0);
}