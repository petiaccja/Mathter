// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Vector/Comparison.hpp>
#include <Mathter/Vector/Vector.hpp>

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <new>

using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Vector - Default initializer", "[Vector]",
						decltype(VectorCaseList<ScalarsFloatingAndComplex, PackingsAll>{})) {

	using Vec = typename TestType::template Vector<3>;

	alignas(alignof(Vec)) std::array<uint8_t, sizeof(Vec)> memoryRegion;

	std::fill(memoryRegion.begin(), memoryRegion.end(), uint8_t(0xCC));


	const auto ptr = reinterpret_cast<Vec*>(memoryRegion.data());
	new (ptr) Vec;

#ifdef NDEBUG
	REQUIRE(std::all_of(memoryRegion.begin(), memoryRegion.end(), [](const auto& v) { return v == 0xCC; }));
#else
	REQUIRE(std::all_of(ptr->begin(), ptr->end(), [](const auto& v) { return std::isnan(std::real(v)); }));
#endif
}


TEMPLATE_LIST_TEST_CASE("Vector - Ctor: conversion", "[Vector]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloatAndInt32AndComplex, PackingsAll>,
												VectorCaseList<ScalarsFloatAndInt32AndComplex, PackingsAll>>{})) {

	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using VecRhs = typename TestType::Lhs::template Vector<3>;

	VecLhs lhs(1, 2, 3);
	VecRhs rhs(lhs);

	REQUIRE(lhs[0] == scalar_type_t<VecLhs>(1));
	REQUIRE(lhs[1] == scalar_type_t<VecLhs>(2));
	REQUIRE(lhs[2] == scalar_type_t<VecLhs>(3));
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


TEMPLATE_LIST_TEST_CASE("Vector - Ctor: from parts", "[Vector]",
						decltype(VectorCaseList<ScalarsAll, PackingsAll>{})) {

	using Vec = typename TestType::template Vector<4>;
	using Swiz = Swizzle<int, 4, is_packed_v<Vec>, 3, 0>;

	Swiz swiz{ 1, 2, 3, 4 };
	Vec vec(swiz, 6, 4);

	REQUIRE(vec[0] == scalar_type_t<Vec>(4));
	REQUIRE(vec[1] == scalar_type_t<Vec>(1));
	REQUIRE(vec[2] == scalar_type_t<Vec>(6));
	REQUIRE(vec[3] == scalar_type_t<Vec>(4));
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