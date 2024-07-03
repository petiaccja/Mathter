// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../TestGenerators.hpp"

#include <Mathter/Common/Approx.hpp>
#include <Mathter/Vector.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstring>
#include <new>

using namespace mathter;
using Catch::Approx;


using TypeListReal = TestTypeList<TypesReal, PackedAll>;
using TypeListAll = TestTypeList<TypesAll, PackedAll>;


TEST_CASE("Vector - GetBatchDim", "[Vector]") {
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


TEST_CASE("Vector deterministic default initializer", "[Init]") {
	using VecT = Vector<float, 3>;
	alignas(alignof(VecT)) std::array<uint8_t, sizeof(VecT)> rawData;
	std::memset(rawData.data(), 0xCC, rawData.size());

	for (auto& v : rawData) {
		REQUIRE(v == uint8_t(0xCC));
	}

	new (rawData.data()) VecT;

#ifdef NDEBUG
	for (auto& v : rawData) {
		REQUIRE(v == uint8_t(0xCC));
	}
#else
	for (auto& v : *reinterpret_cast<const VecT*>(rawData.data())) {
		REQUIRE(std::isnan(v));
	}
#endif
}



TEMPLATE_LIST_TEST_CASE("Vector - CtorAll", "[Vector]", TypeListReal) {
	SECTION(TestType::Name()) {
		using Vec1 = typename TestType::template Vector<1>;
		using Vec2 = typename TestType::template Vector<2>;
		using Vec3 = typename TestType::template Vector<3>;
		using Vec4 = typename TestType::template Vector<4>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec1 v1(10);

		REQUIRE(v1[0] == 10);

		Vec2 v2(10);

		REQUIRE(v2[0] == 10);
		REQUIRE(v2[1] == 10);

		Vec3 v3(10);

		REQUIRE(v3[0] == 10);
		REQUIRE(v3[1] == 10);
		REQUIRE(v3[2] == 10);

		Vec4 v4(10);

		REQUIRE(v4[0] == 10);
		REQUIRE(v4[1] == 10);
		REQUIRE(v4[2] == 10);
		REQUIRE(v4[3] == 10);

		Vec5 v5(10);

		REQUIRE(v5[0] == 10);
		REQUIRE(v5[1] == 10);
		REQUIRE(v5[2] == 10);
		REQUIRE(v5[3] == 10);
		REQUIRE(v5[4] == 10);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Ctor data pointer", "[Vector]", TypeListReal) {
	SECTION(TestType::Name()) {
		using Vec1 = typename TestType::template Vector<1>;
		using Vec2 = typename TestType::template Vector<2>;
		using Vec3 = typename TestType::template Vector<3>;
		using Vec4 = typename TestType::template Vector<4>;
		using Vec5 = typename TestType::template Vector<5>;

		const double data[5] = { 1, 2, 3, 4, 5 };

		Vec1 v1(data);

		REQUIRE(v1[0] == 1);

		Vec2 v2(data);

		REQUIRE(v2[0] == 1);
		REQUIRE(v2[1] == 2);

		Vec3 v3(data);

		REQUIRE(v3[0] == 1);
		REQUIRE(v3[1] == 2);
		REQUIRE(v3[2] == 3);

		Vec4 v4(data);

		REQUIRE(v4[0] == 1);
		REQUIRE(v4[1] == 2);
		REQUIRE(v4[2] == 3);
		REQUIRE(v4[3] == 4);

		Vec5 v5(data);

		REQUIRE(v5[0] == 1);
		REQUIRE(v5[1] == 2);
		REQUIRE(v5[2] == 3);
		REQUIRE(v5[3] == 4);
		REQUIRE(v5[4] == 5);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Ctor conversion", "[Vector]", TypeListReal) {
	SECTION(TestType::Name()) {
		using Vec1 = typename TestType::template Vector<1>;
		using Vec2 = typename TestType::template Vector<2>;
		using Vec3 = typename TestType::template Vector<3>;
		using Vec4 = typename TestType::template Vector<4>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec1 v1 = Vector<double, 1, false>(1);

		REQUIRE(v1[0] == 1);

		Vec2 v2 = Vector<double, 2, false>(1, 2);

		REQUIRE(v2[0] == 1);
		REQUIRE(v2[1] == 2);

		Vec3 v3 = Vector<double, 3, false>(1, 2, 3);

		REQUIRE(v3[0] == 1);
		REQUIRE(v3[1] == 2);
		REQUIRE(v3[2] == 3);

		Vec4 v4 = Vector<double, 4, false>(1, 2, 3, 4);

		REQUIRE(v4[0] == 1);
		REQUIRE(v4[1] == 2);
		REQUIRE(v4[2] == 3);
		REQUIRE(v4[3] == 4);

		Vec5 v5 = Vector<double, 5, false>(1, 2, 3, 4, 5);

		REQUIRE(v5[0] == 1);
		REQUIRE(v5[1] == 2);
		REQUIRE(v5[2] == 3);
		REQUIRE(v5[3] == 4);
		REQUIRE(v5[4] == 5);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Homogeneous cast", "[Vector]", TypeListReal) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec4 = typename TestType::template Vector<4>;
		using Vec5 = typename TestType::template Vector<5>;
		using Vec6 = typename TestType::template Vector<6>;

		Vec3 v3(0, 0, 0);
		Vec4 u3 = (Vec4)v3;
		Vec3 d3 = (Vec3)u3;

		REQUIRE(u3 == Vec4(0, 0, 0, 1));
		REQUIRE(v3 == d3);

		Vec5 v5(0, 0, 0, 0, 0);
		Vec6 u5 = (Vec6)v5;
		Vec5 d5 = (Vec5)u5;

		REQUIRE(u5 == Vec6(0, 0, 0, 0, 0, 1));
		REQUIRE(v5 == d5);
	}
}



TEMPLATE_LIST_TEST_CASE("Vector - Ctor scalar", "[Vector]", TypeListReal) {
	SECTION(TestType::Name()) {
		using Vec1 = typename TestType::template Vector<1>;
		using Vec2 = typename TestType::template Vector<2>;
		using Vec3 = typename TestType::template Vector<3>;
		using Vec4 = typename TestType::template Vector<4>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec1 v1(1);

		REQUIRE(v1[0] == 1);

		Vec2 v2(1, 2);

		REQUIRE(v2[0] == 1);
		REQUIRE(v2[1] == 2);

		Vec3 v3(1, 2, 3);

		REQUIRE(v3[0] == 1);
		REQUIRE(v3[1] == 2);
		REQUIRE(v3[2] == 3);

		Vec4 v4(1, 2, 3, 4);

		REQUIRE(v4[0] == 1);
		REQUIRE(v4[1] == 2);
		REQUIRE(v4[2] == 3);
		REQUIRE(v4[3] == 4);

		Vec5 v5(1, 2, 3, 4, 5);

		REQUIRE(v5[0] == 1);
		REQUIRE(v5[1] == 2);
		REQUIRE(v5[2] == 3);
		REQUIRE(v5[3] == 4);
		REQUIRE(v5[4] == 5);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Ctor mixed", "[Vector]", TypeListReal) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec4 = typename TestType::template Vector<4>;
		using Vec5 = typename TestType::template Vector<5>;

		Vector<double, 2> vd(1, 2);
		Vector<float, 2> vf(3, 4);

		Vec3 v3 = { vd, 3 };

		REQUIRE(v3[0] == 1);
		REQUIRE(v3[1] == 2);
		REQUIRE(v3[2] == 3);

		Vec4 v4 = { vd, vf };

		REQUIRE(v4[0] == 1);
		REQUIRE(v4[1] == 2);
		REQUIRE(v4[2] == 3);
		REQUIRE(v4[3] == 4);

		Vec5 v5 = { vd, 0, vf };

		REQUIRE(v5[0] == 1);
		REQUIRE(v5[1] == 2);
		REQUIRE(v5[2] == 0);
		REQUIRE(v5[3] == 3);
		REQUIRE(v5[4] == 4);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Ctor mixed swizzle", "[Vector]", TypeListReal) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vector<double, 4> source(1, 2, 3, 4);

		Vec3 v3 = { source.wxy };

		REQUIRE(v3[0] == 4);
		REQUIRE(v3[1] == 1);
		REQUIRE(v3[2] == 2);

		Vec5 v5 = { source.xy, 0, source.zw };

		REQUIRE(v5[0] == 1);
		REQUIRE(v5[1] == 2);
		REQUIRE(v5[2] == 0);
		REQUIRE(v5[3] == 3);
		REQUIRE(v5[4] == 4);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Operator []", "[Vector]", TypeListReal) {
	SECTION(TestType::Name()) {
		using Vec4 = typename TestType::template Vector<4>;

		Vec4 source(0, 1, 2, 3);

		REQUIRE(source[0] == 0);
		REQUIRE(source[1] == 1);
		REQUIRE(source[2] == 2);
		REQUIRE(source[3] == 3);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Operator ()", "[Vector]", TypeListReal) {
	SECTION(TestType::Name()) {
		using Vec4 = typename TestType::template Vector<4>;

		Vec4 source(0, 1, 2, 3);

		REQUIRE(source(0) == 0);
		REQUIRE(source(1) == 1);
		REQUIRE(source(2) == 2);
		REQUIRE(source(3) == 3);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Iterators", "[Vector]", TypeListReal) {
	SECTION(TestType::Name()) {
		using Vec4 = typename TestType::template Vector<4>;
		using Type = scalar_type_t<Vec4>;

		Vec4 source(5, 6, 7, 8);

		Type first = Type(5);
		for (const auto& v : source) {
			REQUIRE(v == first);
			first = first + Type(1);
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Swizzle", "[Vector]", TypeListReal) {
	SECTION(TestType::Name()) {
		using Vec2 = typename TestType::template Vector<2>;
		using Vec3 = typename TestType::template Vector<3>;
		using Vec4 = typename TestType::template Vector<4>;

		Vec2 v2 = { 1, 2 };
		Vec3 v3 = { 1, 2, 3 };
		Vec4 v4 = { 1, 2, 3, 4 };

		REQUIRE(Vec2(v2.yx) == Vec2(2, 1));
		REQUIRE(Vec3(v2.yxy) == Vec3(2, 1, 2));
		REQUIRE(Vec4(v2.yxyx) == Vec4(2, 1, 2, 1));

		REQUIRE(Vec2(v3.yz) == Vec2(2, 3));
		REQUIRE(Vec3(v3.yzy) == Vec3(2, 3, 2));
		REQUIRE(Vec4(v3.yzyx) == Vec4(2, 3, 2, 1));

		REQUIRE(Vec2(v4.wz) == Vec2(4, 3));
		REQUIRE(Vec3(v4.wzy) == Vec3(4, 3, 2));
		REQUIRE(Vec4(v4.wzyx) == Vec4(4, 3, 2, 1));
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Comparison", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 value3(1, 2, 3);
		Vec3 equal3(1, 2, 3);
		Vec3 unequal3(1, 2, 10);

		REQUIRE(value3 == equal3);
		REQUIRE(!(value3 == unequal3));
		REQUIRE(!(value3 != equal3));
		REQUIRE(value3 != unequal3);

		Vec5 value5(1, 2, 3, 4, 5);
		Vec5 equal5(1, 2, 3, 4, 5);
		Vec5 unequal5(1, 2, 3, 4, 10);

		REQUIRE(value5 == equal5);
		REQUIRE(!(value5 == unequal5));
		REQUIRE(value5 != unequal5);
		REQUIRE(!(value5 != equal5));
	}
}


TEST_CASE("Vector - IOParse", "[Vector]") {
	Vector<float, 3> parsed;

	std::string successCases[] = {
		"3.14, 2.718, 0.57",
		"  (3.14,2.718  ,  0.57\t  )  ",
		"[3.14\t2.718\t0.57]",
		"{    3.14  2.718, 0.57 } ",
	};
	Vector<float, 3> expected{ 3.14f, 2.718f, 0.57f };

	for (const auto& c : successCases) {
		const char* end;
		parsed = strtovec<decltype(parsed)>(c.c_str(), &end);
		REQUIRE(end != c.c_str());
		REQUIRE(parsed == expected);
	}

	std::string failureCases[] = {
		"[3. 14, 2.718, 0.57]", // more elements
		"(3.1q, 2.718, 0.57\t)", // invalid element
		"[3.14, 2.718, 0.57}", // unmatching brackets
		"{ 3.14, 0.57 } ", // not enough elements
		"[ 3.14, , 2.718, 0.57 ]", // empty elements
	};

	for (const auto& c : failureCases) {
		const char* end;
		parsed = strtovec<Vector<float, 3, false>>(c.c_str(), &end);
		REQUIRE(end == c.c_str());
	}
}