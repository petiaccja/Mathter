// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Cases.hpp"

#include <Mathter/Vector/Comparison.hpp>

#include <catch2/catch_template_test_macros.hpp>

using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Vector - Compare vec x vec", "[Vector]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloatAndInt32, PackingsAll>,
												VectorCaseList<ScalarsFloatAndInt32, PackingsAll>>{})) {

	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using VecRhs = typename TestType::Rhs::template Vector<3>;

	const VecLhs lhs = { 1, 2, 3 };
	const VecRhs equal = { 1, 2, 3 };
	const VecRhs notEqual = { 3, 2, 7 };

	REQUIRE(lhs == equal);
	REQUIRE(lhs != notEqual);
}


TEMPLATE_LIST_TEST_CASE("Vector - Compare vec x swizzle", "[Vector]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloatAndInt32, PackingsAll>,
												SwizzleCaseList<ScalarsFloatAndInt32, PackingsAll>>{})) {

	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using SwizRhs = typename TestType::Rhs::template Swizzle<3, 1, 2, 0>;

	const VecLhs lhs = { 1, 2, 3 };
	const SwizRhs equal = { 3, 1, 2 };
	const SwizRhs notEqual = { 7, 3, 0 };

	REQUIRE(lhs == equal);
	REQUIRE(lhs != notEqual);
}


TEMPLATE_LIST_TEST_CASE("Vector - Compare swizzle x vec", "[Vector]",
						decltype(BinaryCaseList<SwizzleCaseList<ScalarsFloatAndInt32, PackingsAll>,
												VectorCaseList<ScalarsFloatAndInt32, PackingsAll>>{})) {

	using SwizLhs = typename TestType::Lhs::template Swizzle<3, 1, 2, 0>;
	using VecRhs = typename TestType::Rhs::template Vector<3>;

	const VecRhs rhs = { 1, 2, 3 };
	const SwizLhs equal = { 3, 1, 2 };
	const SwizLhs notEqual = { 7, 3, 0 };

	REQUIRE(equal == rhs);
	REQUIRE(notEqual != rhs);
}


TEMPLATE_LIST_TEST_CASE("Vector - Compare swizzle x swizzle", "[Vector]",
						decltype(BinaryCaseList<SwizzleCaseList<ScalarsFloatAndInt32, PackingsAll>,
												SwizzleCaseList<ScalarsFloatAndInt32, PackingsAll>>{})) {

	using SwizLhs = typename TestType::Lhs::template Swizzle<3, 1, 2, 0>;
	using SwizRhs = typename TestType::Lhs::template Swizzle<3, 2, 0, 1>;

	const SwizLhs rhs = { 1, 2, 3 }; // [2, 3, 1]
	const SwizRhs equal = { 3, 1, 2 };
	const SwizRhs notEqual = { 7, 3, 0 };

	REQUIRE(equal == rhs);
	REQUIRE(notEqual != rhs);
}


TEMPLATE_LIST_TEST_CASE("Vector - Compare swizzle x scalar", "[Vector]",
						decltype(SwizzleCaseList<ScalarsFloatAndInt32, PackingsAll>{})) {

	using Swiz = typename TestType::template Swizzle<3, 1>;

	const Swiz lhs = { 1, 2, 3 };

	REQUIRE(lhs == 2);
	REQUIRE(lhs != 3);
}


TEMPLATE_LIST_TEST_CASE("Vector - Compare scalar x swiz", "[Vector]",
						decltype(SwizzleCaseList<ScalarsFloatAndInt32, PackingsAll>{})) {

	using Swiz = typename TestType::template Swizzle<3, 1>;

	const Swiz rhs = { 1, 2, 3 };

	REQUIRE(2 == rhs);
	REQUIRE(3 != rhs);
}