// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Quaternion/Comparison.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Quaternion - Comparison", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using QuatRhs = typename TestType::Rhs::Quat;

	const QuatLhs lhs(1, 2, 3, 4);
	const QuatRhs same(1, 2, 3, 4);
	const QuatRhs different(7, 2, 3, 4);

	REQUIRE(lhs == same);
	REQUIRE(!(lhs != same));
	REQUIRE(lhs != different);
	REQUIRE(!(lhs == different));
}