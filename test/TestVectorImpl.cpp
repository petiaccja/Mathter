//L=============================================================================
//L This software is distributed under the MIT license.
//L Copyright 2021 Péter Kardos
//L=============================================================================

#pragma warning(disable : 4244)

#include <Mathter/Common/Approx.hpp>
#include <Mathter/Vector.hpp>
#include "TestGenerators.hpp"

#include <catch2/catch.hpp>
#include <new>
#include <cstring>

using namespace mathter;



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



TEST_CASE_VEC_VARIANT("Vector - CtorAll", "[Vector]", TypesReal, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<1> v1(10);

		REQUIRE(v1.data[0] == 10);

		VectorT<2> v2(10);

		REQUIRE(v2.data[0] == 10);
		REQUIRE(v2.data[1] == 10);

		VectorT<3> v3(10);

		REQUIRE(v3.data[0] == 10);
		REQUIRE(v3.data[1] == 10);
		REQUIRE(v3.data[2] == 10);

		VectorT<4> v4(10);

		REQUIRE(v4.data[0] == 10);
		REQUIRE(v4.data[1] == 10);
		REQUIRE(v4.data[2] == 10);
		REQUIRE(v4.data[3] == 10);

		VectorT<5> v5(10);

		REQUIRE(v5.data[0] == 10);
		REQUIRE(v5.data[1] == 10);
		REQUIRE(v5.data[2] == 10);
		REQUIRE(v5.data[3] == 10);
		REQUIRE(v5.data[4] == 10);
	}
}


TEST_CASE_VEC_VARIANT("Vector - Ctor data pointer", "[Vector]", TypesReal, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		const double data[5] = { 1, 2, 3, 4, 5 };

		VectorT<1> v1(data);

		REQUIRE(v1.data[0] == 1);

		VectorT<2> v2(data);

		REQUIRE(v2.data[0] == 1);
		REQUIRE(v2.data[1] == 2);

		VectorT<3> v3(data);

		REQUIRE(v3.data[0] == 1);
		REQUIRE(v3.data[1] == 2);
		REQUIRE(v3.data[2] == 3);

		VectorT<4> v4(data);

		REQUIRE(v4.data[0] == 1);
		REQUIRE(v4.data[1] == 2);
		REQUIRE(v4.data[2] == 3);
		REQUIRE(v4.data[3] == 4);

		VectorT<5> v5(data);

		REQUIRE(v5.data[0] == 1);
		REQUIRE(v5.data[1] == 2);
		REQUIRE(v5.data[2] == 3);
		REQUIRE(v5.data[3] == 4);
		REQUIRE(v5.data[4] == 5);
	}
}


TEST_CASE_VEC_VARIANT("Vector - Ctor conversion", "[Vector]", TypesReal, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<1> v1 = Vector<double, 1, false>(1);

		REQUIRE(v1.data[0] == 1);

		VectorT<2> v2 = Vector<double, 2, false>(1, 2);

		REQUIRE(v2.data[0] == 1);
		REQUIRE(v2.data[1] == 2);

		VectorT<3> v3 = Vector<double, 3, false>(1, 2, 3);

		REQUIRE(v3.data[0] == 1);
		REQUIRE(v3.data[1] == 2);
		REQUIRE(v3.data[2] == 3);

		VectorT<4> v4 = Vector<double, 4, false>(1, 2, 3, 4);

		REQUIRE(v4.data[0] == 1);
		REQUIRE(v4.data[1] == 2);
		REQUIRE(v4.data[2] == 3);
		REQUIRE(v4.data[3] == 4);

		VectorT<5> v5 = Vector<double, 5, false>(1, 2, 3, 4, 5);

		REQUIRE(v5.data[0] == 1);
		REQUIRE(v5.data[1] == 2);
		REQUIRE(v5.data[2] == 3);
		REQUIRE(v5.data[3] == 4);
		REQUIRE(v5.data[4] == 5);
	}
}


TEST_CASE_VEC_VARIANT("Vector - Homogeneous cast", "[Vector]", TypesReal, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> v3(0, 0, 0);
		VectorT<4> u3 = (VectorT<4>)v3;
		VectorT<3> d3 = (VectorT<3>)u3;

		REQUIRE(u3 == VectorT<4>(0, 0, 0, 1));
		REQUIRE(v3 == d3);

		VectorT<5> v5(0, 0, 0, 0, 0);
		VectorT<6> u5 = (VectorT<6>)v5;
		VectorT<5> d5 = (VectorT<5>)u5;

		REQUIRE(u5 == VectorT<6>(0, 0, 0, 0, 0, 1));
		REQUIRE(v5 == d5);
	}
}



TEST_CASE_VEC_VARIANT("Vector - Ctor scalar", "[Vector]", TypesReal, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<1> v1(1);

		REQUIRE(v1.data[0] == 1);

		VectorT<2> v2(1, 2);

		REQUIRE(v2.data[0] == 1);
		REQUIRE(v2.data[1] == 2);

		VectorT<3> v3(1, 2, 3);

		REQUIRE(v3.data[0] == 1);
		REQUIRE(v3.data[1] == 2);
		REQUIRE(v3.data[2] == 3);

		VectorT<4> v4(1, 2, 3, 4);

		REQUIRE(v4.data[0] == 1);
		REQUIRE(v4.data[1] == 2);
		REQUIRE(v4.data[2] == 3);
		REQUIRE(v4.data[3] == 4);

		VectorT<5> v5(1, 2, 3, 4, 5);

		REQUIRE(v5.data[0] == 1);
		REQUIRE(v5.data[1] == 2);
		REQUIRE(v5.data[2] == 3);
		REQUIRE(v5.data[3] == 4);
		REQUIRE(v5.data[4] == 5);
	}
}


TEST_CASE_VEC_VARIANT("Vector - Ctor mixed", "[Vector]", TypesReal, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Vector<double, 2> vd(1, 2);
		Vector<float, 2> vf(3, 4);

		VectorT<3> v3 = { vd, 3 };

		REQUIRE(v3.data[0] == 1);
		REQUIRE(v3.data[1] == 2);
		REQUIRE(v3.data[2] == 3);

		VectorT<4> v4 = { vd, vf };

		REQUIRE(v4.data[0] == 1);
		REQUIRE(v4.data[1] == 2);
		REQUIRE(v4.data[2] == 3);
		REQUIRE(v4.data[3] == 4);

		VectorT<5> v5 = { vd, 0, vf };

		REQUIRE(v5.data[0] == 1);
		REQUIRE(v5.data[1] == 2);
		REQUIRE(v5.data[2] == 0);
		REQUIRE(v5.data[3] == 3);
		REQUIRE(v5.data[4] == 4);
	}
}


TEST_CASE_VEC_VARIANT("Vector - Ctor mixed swizzle", "[Vector]", TypesReal, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Vector<double, 4> source(1, 2, 3, 4);

		VectorT<3> v3 = { source.wxy };

		REQUIRE(v3.data[0] == 4);
		REQUIRE(v3.data[1] == 1);
		REQUIRE(v3.data[2] == 2);

		VectorT<5> v5 = { source.xy, 0, source.zw };

		REQUIRE(v5.data[0] == 1);
		REQUIRE(v5.data[1] == 2);
		REQUIRE(v5.data[2] == 0);
		REQUIRE(v5.data[3] == 3);
		REQUIRE(v5.data[4] == 4);
	}
}


TEST_CASE_VEC_VARIANT("Vector - Operator []", "[Vector]", TypesReal, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<4> source(0, 1, 2, 3);

		REQUIRE(source[0] == 0);
		REQUIRE(source[1] == 1);
		REQUIRE(source[2] == 2);
		REQUIRE(source[3] == 3);
	}
}


TEST_CASE_VEC_VARIANT("Vector - Operator ()", "[Vector]", TypesReal, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<4> source(0, 1, 2, 3);

		REQUIRE(source(0) == 0);
		REQUIRE(source(1) == 1);
		REQUIRE(source(2) == 2);
		REQUIRE(source(3) == 3);
	}
}


TEST_CASE_VEC_VARIANT("Vector - Iterators", "[Vector]", TypesReal, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<4> source(5, 6, 7, 8);

		Type first = Type(5);
		for (const auto& v : source) {
			REQUIRE(v == first);
			first = first + Type(1);
		}
	}
}


TEST_CASE_VEC_VARIANT("Vector - Swizzle", "[Vector]", TypesReal, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<2> v2 = { 1, 2 };
		VectorT<3> v3 = { 1, 2, 3 };
		VectorT<4> v4 = { 1, 2, 3, 4 };

		REQUIRE(VectorT<2>(v2.yx) == VectorT<2>(2, 1));
		REQUIRE(VectorT<3>(v2.yxy) == VectorT<3>(2, 1, 2));
		REQUIRE(VectorT<4>(v2.yxyx) == VectorT<4>(2, 1, 2, 1));

		REQUIRE(VectorT<2>(v3.yz) == VectorT<2>(2, 3));
		REQUIRE(VectorT<3>(v3.yzy) == VectorT<3>(2, 3, 2));
		REQUIRE(VectorT<4>(v3.yzyx) == VectorT<4>(2, 3, 2, 1));

		REQUIRE(VectorT<2>(v4.wz) == VectorT<2>(4, 3));
		REQUIRE(VectorT<3>(v4.wzy) == VectorT<3>(4, 3, 2));
		REQUIRE(VectorT<4>(v4.wzyx) == VectorT<4>(4, 3, 2, 1));
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
		parsed = strtovec<float, 3, false>(c.c_str(), &end);
		REQUIRE(end == c.c_str());
	}
}