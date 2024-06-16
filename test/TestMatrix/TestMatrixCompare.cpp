// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../TestGenerators.hpp"

#include <Mathter/Common/Approx.hpp>
#include <Mathter/Common/Traits.hpp>
#include <Mathter/Matrix.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>


using namespace mathter;
using Catch::Approx;


using TypeListAll = TestTypeList<TypesAll, PackedAll, OrdersAll, LayoutsAll>;


TEMPLATE_LIST_TEST_CASE("Matrix - Compare", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M32 = typename TestType::template Matrix<3, 2>;
		using M32L = invert_layout_t<M32>;
		using M23O = invert_order_t<M32>;
		M32 value = { 1, 2, 3, 4, 5, 6 };

		SECTION("equal values") {
			M32L valueL = { 1, 2, 3, 4, 5, 6 };
			M23O valueO = { 1, 3, 5, 2, 4, 6 };

			REQUIRE(value == valueL);
			REQUIRE(value == valueO);
			REQUIRE(!(value != valueL));
			REQUIRE(!(value != valueO));
		}

		SECTION("unequal values") {
			M32L valueL = { 1, 2, 3, 4, 5, 7 };
			M23O valueO = { 1, 3, 5, 2, 4, 7 };

			REQUIRE(value != valueL);
			REQUIRE(value != valueO);
			REQUIRE(!(value == valueL));
			REQUIRE(!(value == valueO));
		}
	}
}