// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../TestGenerators.hpp"

#include <Mathter/Common/Approx.hpp>
#include <Mathter/Matrix.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>
#include <cstring>
#include <new>



using namespace mathter;
using Catch::Approx;


using TypeListAll = TestTypeList<TypesFloating, PackedAll, OrdersAll, LayoutsAll>;
using TypeListFloating = TestTypeList<TypesFloating, PackedAll, OrdersAll, LayoutsAll>;


TEST_CASE("Matrix deterministic default initializer", "[Init]") {
	using MatT = Matrix<float, 3, 3>;
	alignas(alignof(MatT)) std::array<uint8_t, sizeof(MatT)> rawData;
	std::memset(rawData.data(), 0xCC, rawData.size());

	for (auto& v : rawData) {
		REQUIRE(v == uint8_t(0xCC));
	}

	new (rawData.data()) MatT;

#ifdef NDEBUG
	for (auto& v : rawData) {
		REQUIRE(v == uint8_t(0xCC));
	}
#else
	const MatT& m = *reinterpret_cast<const MatT*>(rawData.data());
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			REQUIRE(std::isnan(m(i, j)));
		}
	}
#endif
}



TEMPLATE_LIST_TEST_CASE("Matrix - Constructor & indexer", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M = typename TestType::template Matrix<3, 3>;
		using Type = scalar_type_t<M>;

		M m = {
			Type(1), Type(2), Type(3),
			Type(4), Type(5), Type(6),
			Type(7), Type(8), Type(9)
		};
		M n;
		n(0, 0) = Type(1);
		n(0, 1) = Type(2);
		n(0, 2) = Type(3);
		n(1, 0) = Type(4);
		n(1, 1) = Type(5);
		n(1, 2) = Type(6);
		n(2, 0) = Type(7);
		n(2, 1) = Type(8);
		n(2, 2) = Type(9);

		REQUIRE(m == n);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Converting constructor", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M = typename TestType::template Matrix<3, 3>;
		constexpr auto Order = order_v<M>;

		Matrix<float, 3, 3, Order> m = {
			1, 2, 3,
			4, 5, 6,
			7, 8, 9
		};
		const auto n = M(m);

		REQUIRE(m == ApproxVec(n));
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Thin matrix from vector", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M13 = typename TestType::template Matrix<1, 3>;
		using M31 = typename TestType::template Matrix<3, 1>;
		using V = typename TestType::template Vector<3>;

		V v = { 1, 2, 3 };

		M31 m1 = v;
		M13 m2 = v;

		REQUIRE(m1(0, 0) == 1);
		REQUIRE(m1(1, 0) == 2);
		REQUIRE(m1(2, 0) == 3);

		REQUIRE(m2(0, 0) == 1);
		REQUIRE(m2(0, 1) == 2);
		REQUIRE(m2(0, 2) == 3);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Thin matrix to vector", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M13 = typename TestType::template Matrix<1, 3>;
		using M31 = typename TestType::template Matrix<3, 1>;
		using V = typename TestType::template Vector<3>;

		V vexp = { 1, 2, 3 };

		M31 m1 = { 1, 2, 3 };
		M13 m2 = { 1, 2, 3 };

		V v1 = m1;
		V v2 = m2;

		REQUIRE(v1 == vexp);
		REQUIRE(v2 == vexp);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Thin matrix short indexing", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M13 = typename TestType::template Matrix<1, 3>;
		using M31 = typename TestType::template Matrix<3, 1>;

		M31 m1 = { 1, 2, 3 };
		M13 m2 = { 1, 2, 3 };

		REQUIRE(m1(0) == 1);
		REQUIRE(m1(1) == 2);
		REQUIRE(m1(2) == 3);

		REQUIRE(m2(0) == 1);
		REQUIRE(m2(1) == 2);
		REQUIRE(m2(2) == 3);
	}
}


TEST_CASE("Matrix - Extract", "[Matrix]") {
	REQUIRE(false && "not implemented");
}


TEST_CASE("Matrix - Insert", "[Matrix]") {
	REQUIRE(false && "not implemented");
}


TEST_CASE("Matrix - IOParse", "[Matrix]") {
	Matrix<float, 2, 2> parsed;

	std::string successCases[] = {
		"[3.14, 2.718; 0.57, 6.63]",
		"[3.14  2.718; 0.57  6.63  ]",
		"[[3.14, 2.718];\n[0.57, 6.63]]",
		"3.14, 2.718; 0.57, 6.63   ]",
	};
	Matrix<float, 2, 2> expected{
		3.14f, 2.718f,
		0.57f, 6.63f
	};

	for (const auto& c : successCases) {
		const char* end;
		parsed = strtomat<decltype(parsed)>(c.c_str(), &end);
		REQUIRE(end != c.c_str());
		REQUIRE(parsed == expected);
	}

	std::string failureCases[] = {
		"[3.14, 2.718\n 0.57, 6.63]", // missing row delimiter
		"[3.14  2.718h 0.57  6.63  ]", // invalid delimiter
		"[3.14, 2.718; 0.57, 6.63; 1.38, 6.02]", // too many rows
		"[3.14, 2.718 ]", // too few rows
	};

	for (const auto& c : failureCases) {
		const char* end;
		parsed = strtomat<decltype(parsed)>(c.c_str(), &end);
		REQUIRE(end == c.c_str());
	}
}