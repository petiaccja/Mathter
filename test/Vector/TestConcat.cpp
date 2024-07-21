// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================


#include <Mathter/Vector/Concat.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace mathter;


TEST_CASE("Vector - concat", "[Vector]") {
	using Vec = Vector<int, 3, false>;

	const Vec src = { 1, 2, 3 };

	SECTION("Vec | Vec") {
		const auto result = src | src;
		REQUIRE(result.Dimension() == 6);
		REQUIRE(result[0] == 1);
		REQUIRE(result[1] == 2);
		REQUIRE(result[2] == 3);
		REQUIRE(result[3] == 1);
		REQUIRE(result[4] == 2);
		REQUIRE(result[5] == 3);
	}
	SECTION("Vec | Swizzle") {
		const auto result = src | src.yx;
		REQUIRE(result.Dimension() == 5);
		REQUIRE(result[0] == 1);
		REQUIRE(result[1] == 2);
		REQUIRE(result[2] == 3);
		REQUIRE(result[3] == 2);
		REQUIRE(result[4] == 1);
	}
	SECTION("Vec | Scalar") {
		const auto result = src | 9;
		REQUIRE(result.Dimension() == 4);
		REQUIRE(result[0] == 1);
		REQUIRE(result[1] == 2);
		REQUIRE(result[2] == 3);
		REQUIRE(result[3] == 9);
	}
	SECTION("Swizzle | Swizzle") {
		const auto result = src.zx | src.yx;
		REQUIRE(result.Dimension() == 4);
		REQUIRE(result[0] == 3);
		REQUIRE(result[1] == 1);
		REQUIRE(result[2] == 2);
		REQUIRE(result[3] == 1);
	}
	SECTION("Swizzle | Scalar") {
		const auto result = src.yx | 9;
		REQUIRE(result.Dimension() == 3);
		REQUIRE(result[0] == 2);
		REQUIRE(result[1] == 1);
		REQUIRE(result[2] == 9);
	}
}