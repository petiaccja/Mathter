// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Vector/Comparison.hpp>
#include <Mathter/Vector/Vector.hpp>

#include <catch2/catch_template_test_macros.hpp>
#include <complex>
#include <cstring>
#include <new>

using namespace mathter;
using namespace test_util;


#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-attributes"
#endif


[[gnu::noinline, msvc::noinline]] void TestDefaultInit_Init(void* data, size_t size) {
	std::memset(data, 0xCC, size);
}

template <class T>
[[gnu::noinline, msvc::noinline]] void TestDefaultInit_Construct(void* data) {
	new (data) T;
}

[[gnu::noinline, msvc::noinline]] void TestDefaultInit_Copy(void* actual, const void* data, size_t size) {
	std::memcpy(actual, data, size);
}


TEMPLATE_LIST_TEST_CASE("Vector - Default initializer", "[Vector]",
						decltype(VectorCaseList<ScalarsFloatingAndComplex, PackingsAll>{})) {

	using Vec = typename TestType::template Vector<3>;

#ifdef NDEBUG
	{
		constexpr auto totalSize = sizeof(Vec);
		constexpr auto validSize = sizeof(scalar_type_t<Vec>) * 3;

		alignas(alignof(Vec)) uint8_t bytes[totalSize];
		uint8_t expected[totalSize];
		uint8_t actual[totalSize];

		std::memset(expected, 0xCC, validSize);
		std::memset(expected + validSize, 0x00, totalSize - validSize);

		TestDefaultInit_Init(bytes, totalSize);
		TestDefaultInit_Construct<Vec>(bytes);
		TestDefaultInit_Copy(actual, bytes, totalSize);

		REQUIRE(std::memcmp(actual, expected, totalSize) == 0);
	}
#else
	{
		const Vec v;
		for (auto& element : v) {
			REQUIRE(std::isnan(std::real(element)));
		}
	}
#endif
}


TEMPLATE_LIST_TEST_CASE("Vector - Ctor: conversion", "[Vector]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloatAndInt32AndComplex, PackingsAll>,
												VectorCaseList<ScalarsFloatAndInt32AndComplex, PackingsAll>>{})) {

	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using VecRhs = typename TestType::Lhs::template Vector<3>;

	VecLhs lhs(1, 2, 3);
	VecRhs rhs(lhs);

	REQUIRE(rhs[0] == scalar_type_t<VecRhs>(1));
	REQUIRE(rhs[1] == scalar_type_t<VecRhs>(2));
	REQUIRE(rhs[2] == scalar_type_t<VecRhs>(3));
}


TEMPLATE_LIST_TEST_CASE("Vector - Ctor: homogeneous upcast", "[Vector]",
						decltype(VectorCaseList<ScalarsAll, PackingsAll>{})) {

	using VecSrc = typename TestType::template Vector<3>;
	using VecDst = typename TestType::template Vector<4>;

	VecSrc src(2, 4, 6);
	VecDst dst(src);

	REQUIRE(dst[0] == scalar_type_t<VecDst>(2));
	REQUIRE(dst[1] == scalar_type_t<VecDst>(4));
	REQUIRE(dst[2] == scalar_type_t<VecDst>(6));
	REQUIRE(dst[3] == scalar_type_t<VecDst>(1));
}


TEMPLATE_LIST_TEST_CASE("Vector - Ctor: homogeneous downcast", "[Vector]",
						decltype(VectorCaseList<ScalarsAll, PackingsAll>{})) {

	using VecSrc = typename TestType::template Vector<4>;
	using VecDst = typename TestType::template Vector<3>;

	VecSrc src(2, 4, 6, 2);
	VecDst dst(src);

	REQUIRE(dst[0] == scalar_type_t<VecSrc>(1));
	REQUIRE(dst[1] == scalar_type_t<VecSrc>(2));
	REQUIRE(dst[2] == scalar_type_t<VecSrc>(3));
}


TEMPLATE_LIST_TEST_CASE("Vector - Ctor: pointer to elements", "[Vector]",
						decltype(VectorCaseList<ScalarsAll, PackingsAll>{})) {

	using Vec = typename TestType::template Vector<3>;

	std::array<scalar_type_t<Vec>, 3> elements = { 1, 2, 3 };

	Vec vec(elements.data());

	REQUIRE(vec[0] == scalar_type_t<Vec>(1));
	REQUIRE(vec[1] == scalar_type_t<Vec>(2));
	REQUIRE(vec[2] == scalar_type_t<Vec>(3));
}


TEMPLATE_LIST_TEST_CASE("Vector - Ctor: same elements", "[Vector]",
						decltype(VectorCaseList<ScalarsAll, PackingsAll>{})) {

	using Vec = typename TestType::template Vector<3>;

	Vec vec(1);

	REQUIRE(vec[0] == scalar_type_t<Vec>(1));
	REQUIRE(vec[1] == scalar_type_t<Vec>(1));
	REQUIRE(vec[2] == scalar_type_t<Vec>(1));
}


TEMPLATE_LIST_TEST_CASE("Vector - Ctor: from swizzle", "[Vector]",
						decltype(VectorCaseList<ScalarsAll, PackingsAll>{})) {

	using Vec = typename TestType::template Vector<3>;
	using Swiz = Swizzle<int, 4, is_packed_v<Vec>, 3, 0, 1>;

	Swiz swiz{ 1, 2, 3, 4 };
	Vec vec(swiz);

	REQUIRE(vec[0] == scalar_type_t<Vec>(4));
	REQUIRE(vec[1] == scalar_type_t<Vec>(1));
	REQUIRE(vec[2] == scalar_type_t<Vec>(2));
}


TEST_CASE("Vector - concat helpers: GetConcatDim", "[Vector]") {
	static_assert(impl::GetConcatDim<int, Vector<float, 2, false>>() == 3);
	static_assert(impl::GetConcatDim<Swizzle<float, 3, false, 1, 2>, Vector<float, 2, false>>() == 4);
}


TEST_CASE("Vector - concat helpers: GetConcatPacking", "[Vector]") {
	static_assert(impl::GetConcatPacking<Vector<float, 2, false>, Vector<float, 2, false>>() == false);
	static_assert(impl::GetConcatPacking<Vector<float, 2, false>, Vector<float, 2, true>>() == false);
	static_assert(impl::GetConcatPacking<Vector<float, 2, true>, Vector<float, 2, true>>() == true);
}


TEMPLATE_LIST_TEST_CASE("Vector - Ctor: from parts", "[Vector]",
						decltype(VectorCaseList<ScalarsAll, PackingsAll>{})) {

	using Vec = typename TestType::template Vector<4>;
	using Swiz = Swizzle<int, 4, is_packed_v<Vec>, 3, 0>;

	const Swiz swiz{ 1, 2, 3, 4 };
	SECTION("Explicit vector type") {
		const Vec vec(swiz, 6, 4);

		REQUIRE(vec[0] == scalar_type_t<Vec>(4));
		REQUIRE(vec[1] == scalar_type_t<Vec>(1));
		REQUIRE(vec[2] == scalar_type_t<Vec>(6));
		REQUIRE(vec[3] == scalar_type_t<Vec>(4));
	}
	SECTION("CTAD") {
		const auto vec = Vector(swiz, 6, 4);
		using DeducedVec = std::decay_t<decltype(vec)>;

		REQUIRE(vec[0] == scalar_type_t<DeducedVec>(4));
		REQUIRE(vec[1] == scalar_type_t<DeducedVec>(1));
		REQUIRE(vec[2] == scalar_type_t<DeducedVec>(6));
		REQUIRE(vec[3] == scalar_type_t<DeducedVec>(4));

		REQUIRE(dimension_v<DeducedVec> == 4);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - operator[]", "[Vector]",
						decltype(VectorCaseList<ScalarsFloat32, PackingsAll>{})) {

	using Vec = typename TestType::template Vector<3>;

	Vec vec;
	vec[0] = 1;
	vec[1] = 2;
	vec[2] = 3;

	REQUIRE(vec[0] == 1);
	REQUIRE(vec[1] == 2);
	REQUIRE(vec[2] == 3);
}


TEMPLATE_LIST_TEST_CASE("Vector - operator()", "[Vector]",
						decltype(VectorCaseList<ScalarsFloat32, PackingsAll>{})) {

	using Vec = typename TestType::template Vector<3>;

	Vec vec;
	vec(0) = 1;
	vec(1) = 2;
	vec(2) = 3;

	REQUIRE(vec(0) == 1);
	REQUIRE(vec(1) == 2);
	REQUIRE(vec(2) == 3);
}


TEMPLATE_LIST_TEST_CASE("Vector - iterators", "[Vector]",
						decltype(VectorCaseList<ScalarsFloat32, PackingsAll>{})) {

	using Vec = typename TestType::template Vector<3>;

	Vec vec(1, 2, 3);
	const Vec& cvec = vec;

	SECTION("mutable") {
		auto it = vec.begin();

		REQUIRE(*it++ == 1);
		REQUIRE(*it++ == 2);
		REQUIRE(*it++ == 3);
		REQUIRE(it == vec.end());
	}
	SECTION("const by ref") {
		auto it = cvec.begin();

		REQUIRE(*it++ == 1);
		REQUIRE(*it++ == 2);
		REQUIRE(*it++ == 3);
		REQUIRE(it == cvec.end());
	}
	SECTION("const by request") {
		auto it = vec.cbegin();

		REQUIRE(*it++ == 1);
		REQUIRE(*it++ == 2);
		REQUIRE(*it++ == 3);
		REQUIRE(it == vec.cend());
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - data", "[Vector]",
						decltype(VectorCaseList<ScalarsFloat32, PackingsAll>{})) {

	using Vec = typename TestType::template Vector<3>;

	Vec vec(1, 2, 3);
	const Vec& cvec = vec;

	SECTION("mutable") {
		REQUIRE(vec.Data()[0] == 1);
		REQUIRE(vec.Data()[1] == 2);
		REQUIRE(vec.Data()[2] == 3);
	}
	SECTION("const") {
		REQUIRE(cvec.Data()[0] == 1);
		REQUIRE(cvec.Data()[1] == 2);
		REQUIRE(cvec.Data()[2] == 3);
	}
}