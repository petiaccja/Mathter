// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include <Mathter/Vector/SIMDUtil.hpp>

#include <catch2/catch_test_macros.hpp>


using namespace mathter;


TEST_CASE("SIMDUtil - GetBatchDim", "[SIMDUtil]") {
	REQUIRE(GetBatchSize(1, false) == 1);
	REQUIRE(GetBatchSize(2, false) == 2);
	REQUIRE(GetBatchSize(3, false) == 4);
	REQUIRE(GetBatchSize(4, false) == 4);
	REQUIRE(GetBatchSize(5, false) == 8);
	REQUIRE(GetBatchSize(6, false) == 8);
	REQUIRE(GetBatchSize(7, false) == 8);
	REQUIRE(GetBatchSize(8, false) == 8);

	REQUIRE(GetBatchSize(1, true) == 1);
	REQUIRE(GetBatchSize(2, true) == 2);
	REQUIRE(GetBatchSize(3, true) == 3);
	REQUIRE(GetBatchSize(4, true) == 4);
	REQUIRE(GetBatchSize(5, true) == 5);
	REQUIRE(GetBatchSize(6, true) == 6);
	REQUIRE(GetBatchSize(7, true) == 7);
	REQUIRE(GetBatchSize(8, true) == 8);
}


TEST_CASE("SIMDUtil - IsBatched", "[SIMDUtil]") {
#ifdef MATHTER_ENABLE_SIMD
	SECTION("SSE2") {
		if constexpr (xsimd::sse2::supported()) {
			REQUIRE(IsBatched<float, 3, false>());
			REQUIRE(IsBatched<float, 4, false>());
		}
	}
#endif
	SECTION("Packed") {
		REQUIRE(!IsBatched<float, 3, true>());
		REQUIRE(!IsBatched<float, 4, true>());
	}
}


TEST_CASE("SIMDUtil - GetStorageSize", "[SIMDUtil]") {
#ifdef MATHTER_ENABLE_SIMD
	SECTION("SSE2") {
		if constexpr (xsimd::sse2::supported()) {
			REQUIRE(GetStorageSize<float, 3, false>() == 4);
			REQUIRE(GetStorageSize<float, 4, false>() == 4);
			REQUIRE(GetStorageSize<float, 17, false>() == 17);
		}
	}
#endif
	SECTION("Packed") {
		REQUIRE(GetStorageSize<float, 3, true>() == 3);
		REQUIRE(GetStorageSize<float, 4, true>() == 4);
		REQUIRE(GetStorageSize<float, 17, true>() == 17);
	}
}


TEST_CASE("SIMDUtil - Alignment", "[SIMDUtil]") {
#ifdef MATHTER_ENABLE_SIMD
	SECTION("SSE2") {
		if constexpr (xsimd::sse2::supported()) {
			REQUIRE(GetStorageAlignment<float, 3, false>() == 16);
			REQUIRE(GetStorageAlignment<float, 4, false>() == 16);
			REQUIRE(GetStorageAlignment<float, 17, false>() == 4);
		}
	}
#endif
	SECTION("Packed") {
		REQUIRE(GetStorageAlignment<float, 3, true>() == 4);
		REQUIRE(GetStorageAlignment<float, 4, true>() == 4);
		REQUIRE(GetStorageAlignment<float, 17, true>() == 4);
	}
}