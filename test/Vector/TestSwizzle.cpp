// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Cases.hpp"

#include <Mathter/Vector/Comparison.hpp>
#include <Mathter/Vector/Vector.hpp>

#include <catch2/catch_template_test_macros.hpp>

using namespace mathter;
using namespace test_util;



TEMPLATE_LIST_TEST_CASE("Swizzle - operator[]", "[Swizzle]",
						decltype(SwizzleCaseList<ScalarsFloating, PackingsAll>{})) {

	using Swiz = typename TestType::template Swizzle<3, 2, 0>;

	Swiz swiz{ 1, 2, 3 };
	const auto& cswiz = swiz;

	SECTION("mutable") {
		REQUIRE(swiz[0] == 3);
		REQUIRE(swiz[1] == 1);
	}
	SECTION("const") {
		REQUIRE(cswiz[0] == 3);
		REQUIRE(cswiz[1] == 1);
	}
}


TEMPLATE_LIST_TEST_CASE("Swizzle - operator()", "[Swizzle]",
						decltype(SwizzleCaseList<ScalarsFloating, PackingsAll>{})) {

	using Swiz = typename TestType::template Swizzle<3, 2, 0>;

	Swiz swiz{ 1, 2, 3 };
	const auto& cswiz = swiz;

	SECTION("mutable") {
		REQUIRE(swiz(0) == 3);
		REQUIRE(swiz(1) == 1);
	}
	SECTION("const") {
		REQUIRE(cswiz(0) == 3);
		REQUIRE(cswiz(1) == 1);
	}
}


TEMPLATE_LIST_TEST_CASE("Swizzle - assign swizzle", "[Swizzle]",
						decltype(BinaryCaseList<SwizzleCaseList<ScalarsInt32, PackingsAll>,
												SwizzleCaseList<ScalarsAll, PackingsAll>>{})) {

	using SwizSrc = typename TestType::Lhs::template Swizzle<3, 2, 0>;
	using SwizDst = typename TestType::Rhs::template Swizzle<3, 1, 2>;

	SwizSrc src{ 1, 2, 3 }; // [3, 1]
	SwizDst dst{ 0, 0, 0 }; // [_, *0, *1]

	dst = src;

	REQUIRE(dst.array[0] == scalar_type_t<SwizDst>(0));
	REQUIRE(dst.array[1] == scalar_type_t<SwizDst>(3));
	REQUIRE(dst.array[2] == scalar_type_t<SwizDst>(1));
}


TEMPLATE_LIST_TEST_CASE("Swizzle - assign vector", "[Swizzle]",
						decltype(BinaryCaseList<SwizzleCaseList<ScalarsAll, PackingsAll>,
												VectorCaseList<ScalarsInt32, PackingsAll>>{})) {

	using SwizDst = typename TestType::Lhs::template Swizzle<3, 1, 2>;
	using VecSrc = typename TestType::Rhs::template Vector<2>;

	VecSrc src{ 3, 1 }; // [3, 1]
	SwizDst dst{ 0, 0, 0 }; // [_, *0, *1]

	dst = src;

	REQUIRE(dst.array[0] == scalar_type_t<SwizDst>(0));
	REQUIRE(dst.array[1] == scalar_type_t<SwizDst>(3));
	REQUIRE(dst.array[2] == scalar_type_t<SwizDst>(1));
}


TEMPLATE_LIST_TEST_CASE("Swizzle - <1>: assign scalar", "[Swizzle]",
						decltype(SwizzleCaseList<ScalarsAll, PackingsAll>{})) {

	using SwizDst = typename TestType::template Swizzle<3, 1>;

	SwizDst dst{ 0, 0, 0 }; // [_, *0, *1]

	dst = 5;

	REQUIRE(dst.array[0] == scalar_type_t<SwizDst>(0));
	REQUIRE(dst.array[1] == scalar_type_t<SwizDst>(5));
	REQUIRE(dst.array[2] == scalar_type_t<SwizDst>(0));
}


TEMPLATE_LIST_TEST_CASE("Swizzle - <1>: convert to scalar", "[Swizzle]",
						decltype(SwizzleCaseList<ScalarsAll, PackingsAll>{})) {

	using SwizDst = typename TestType::template Swizzle<3, 1>;

	SwizDst dst{ 1, 2, 3 }; // [_, *0, *1]

	const scalar_type_t<SwizDst> value = dst;

	REQUIRE(value == scalar_type_t<SwizDst>(2));
}