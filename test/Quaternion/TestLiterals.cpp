// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Quaternion/Literals.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


TEST_CASE("Quaternion - Literals (float)", "[Quaternion]") {
	using namespace quat_literals;

	SECTION("Real") {
		const Quaternion q = 1.0f;
		static_assert(std::is_same_v<scalar_type_t<std::decay_t<decltype(q)>>, float>);
		REQUIRE(q.s == 1);
		REQUIRE(q.i == 0);
		REQUIRE(q.j == 0);
		REQUIRE(q.k == 0);
	}
	SECTION("i") {
		const auto q = 1.0_if;
		static_assert(std::is_same_v<scalar_type_t<std::decay_t<decltype(q)>>, float>);
		REQUIRE(q.s == 0);
		REQUIRE(q.i == 1);
		REQUIRE(q.j == 0);
		REQUIRE(q.k == 0);
	}
	SECTION("j") {
		const auto q = 1.0_jf;
		static_assert(std::is_same_v<scalar_type_t<std::decay_t<decltype(q)>>, float>);
		REQUIRE(q.s == 0);
		REQUIRE(q.i == 0);
		REQUIRE(q.j == 1);
		REQUIRE(q.k == 0);
	}
	SECTION("k") {
		const auto q = 1.0_kf;
		static_assert(std::is_same_v<scalar_type_t<std::decay_t<decltype(q)>>, float>);
		REQUIRE(q.s == 0);
		REQUIRE(q.i == 0);
		REQUIRE(q.j == 0);
		REQUIRE(q.k == 1);
	}
}


TEST_CASE("Quaternion - Literals (double)", "[Quaternion]") {
	using namespace quat_literals;

	SECTION("Real") {
		const Quaternion q = 1.0;
		static_assert(std::is_same_v<scalar_type_t<std::decay_t<decltype(q)>>, double>);
		REQUIRE(q.s == 1);
		REQUIRE(q.i == 0);
		REQUIRE(q.j == 0);
		REQUIRE(q.k == 0);
	}
	SECTION("i") {
		const auto q = 1.0_i;
		static_assert(std::is_same_v<scalar_type_t<std::decay_t<decltype(q)>>, double>);
		REQUIRE(q.s == 0);
		REQUIRE(q.i == 1);
		REQUIRE(q.j == 0);
		REQUIRE(q.k == 0);
	}
	SECTION("j") {
		const auto q = 1.0_j;
		static_assert(std::is_same_v<scalar_type_t<std::decay_t<decltype(q)>>, double>);
		REQUIRE(q.s == 0);
		REQUIRE(q.i == 0);
		REQUIRE(q.j == 1);
		REQUIRE(q.k == 0);
	}
	SECTION("k") {
		const auto q = 1.0_k;
		static_assert(std::is_same_v<scalar_type_t<std::decay_t<decltype(q)>>, double>);
		REQUIRE(q.s == 0);
		REQUIRE(q.i == 0);
		REQUIRE(q.j == 0);
		REQUIRE(q.k == 1);
	}
}


TEST_CASE("Quaternion - Literals (long double)", "[Quaternion]") {
	using namespace quat_literals;

	SECTION("Real") {
		const Quaternion q = 1.0l;
		static_assert(std::is_same_v<scalar_type_t<std::decay_t<decltype(q)>>, long double>);
		REQUIRE(q.s == 1);
		REQUIRE(q.i == 0);
		REQUIRE(q.j == 0);
		REQUIRE(q.k == 0);
	}
	SECTION("i") {
		const auto q = 1.0_il;
		static_assert(std::is_same_v<scalar_type_t<std::decay_t<decltype(q)>>, long double>);
		REQUIRE(q.s == 0);
		REQUIRE(q.i == 1);
		REQUIRE(q.j == 0);
		REQUIRE(q.k == 0);
	}
	SECTION("j") {
		const auto q = 1.0_jl;
		static_assert(std::is_same_v<scalar_type_t<std::decay_t<decltype(q)>>, long double>);
		REQUIRE(q.s == 0);
		REQUIRE(q.i == 0);
		REQUIRE(q.j == 1);
		REQUIRE(q.k == 0);
	}
	SECTION("k") {
		const auto q = 1.0_kl;
		static_assert(std::is_same_v<scalar_type_t<std::decay_t<decltype(q)>>, long double>);
		REQUIRE(q.s == 0);
		REQUIRE(q.i == 0);
		REQUIRE(q.j == 0);
		REQUIRE(q.k == 1);
	}
}