//L=============================================================================
//L This software is distributed under the MIT license.
//L Copyright 2021 Péter Kardos
//L=============================================================================

#pragma warning(disable : 4244)

#include <Mathter/Common/Approx.hpp>
#include <Mathter/Matrix.hpp>
#include "TestGenerators.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <complex>
#include <new>
#include <cstring>


using namespace mathter;
using Catch::Approx;



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



TEST_CASE_VARIANT("Matrix - Constructor & indexer", "[Matrix]", TypesAll, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> m = {
			Type(1), Type(2), Type(3),
			Type(4), Type(5), Type(6),
			Type(7), Type(8), Type(9)
		};
		MatrixT<3, 3> n;
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


TEST_CASE_VARIANT("Matrix - Converting constructor", "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		Matrix<float, 3, 3, Order> m = {
			1, 2, 3,
			4, 5, 6,
			7, 8, 9
		};
		MatrixT<3, 3> n = MatrixT<3, 3>(m);


		REQUIRE(m == ApproxVec(n));
	}
}


TEST_CASE_VARIANT("Matrix - Thin matrix from vector", "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		Vector<Type, 3> v = { 1, 2, 3 };

		MatrixT<3, 1> m1 = v;
		MatrixT<1, 3> m2 = v;

		REQUIRE(m1(0, 0) == 1);
		REQUIRE(m1(1, 0) == 2);
		REQUIRE(m1(2, 0) == 3);

		REQUIRE(m2(0, 0) == 1);
		REQUIRE(m2(0, 1) == 2);
		REQUIRE(m2(0, 2) == 3);
	}
}


TEST_CASE_VARIANT("Matrix - Thin matrix to vector", "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		Vector<Type, 3> vexp = { 1, 2, 3 };

		MatrixT<3, 1> m1 = { 1, 2, 3 };
		MatrixT<1, 3> m2 = { 1, 2, 3 };

		Vector<Type, 3> v1 = m1;
		Vector<Type, 3> v2 = m2;

		REQUIRE(v1 == vexp);
		REQUIRE(v2 == vexp);
	}
}


TEST_CASE_VARIANT("Matrix - Thin matrix short indexing", "[Matrix]", TypesFloating, OrdersAll, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 1> m1 = { 1, 2, 3 };
		MatrixT<1, 3> m2 = { 1, 2, 3 };

		REQUIRE(m1(0) == 1);
		REQUIRE(m1(1) == 2);
		REQUIRE(m1(2) == 3);

		REQUIRE(m2(0) == 1);
		REQUIRE(m2(1) == 2);
		REQUIRE(m2(2) == 3);
	}
}

TEST_CASE("Matrix - Submatrix", "[Matrix]") {
	Matrix<char, 5, 5> m1 = {
		'a', 'b', 'c', 'd', 'e',
		'f', 'g', 'h', 'i', 'j',
		'k', 'l', 'm', 'n', 'o',
		'p', 'q', 'r', 's', 't',
		'u', 'v', 'w', 'x', 'y'
	};

	Matrix<char, 5, 5> m2 = {
		'z', 'z', 'z', 'z', 'z',
		'z', 'z', 'z', 'z', 'z',
		'z', 'z', 'z', 'z', 'z',
		'z', 'z', 'z', 'z', 'z',
		'z', 'z', 'z', 'z', 'z'
	};

	Matrix<char, 5, 5> r = {
		'z', 'z', 'z', 'p', 'q',
		'z', 'z', 'z', 'u', 'v',
		'c', 'd', 'e', 'z', 'z',
		'h', 'i', 'j', 'z', 'z',
		'm', 'n', 'o', 'z', 'z'
	};

	Matrix<char, 2, 2> sm = m1.Submatrix<2, 2>(3, 0);
	m2.Submatrix<3, 3>(2, 0) = m1.Submatrix<3, 3>(0, 2);
	m2.Submatrix<2, 2>(0, 3) = sm;
	REQUIRE(m2 == r);

	m2.Column(4) = Vector<float, 5>('0');
	r(0, 4) = r(1, 4) = r(2, 4) = r(3, 4) = r(4, 4) = '0';
	REQUIRE(m2 == r);


	Vector<char, 3> v = m1.Submatrix<3, 1>(0, 0);
	Vector<char, 3> vr = { 'a', 'f', 'k' };
	REQUIRE(v == vr);
	v = m1.Submatrix<1, 3>(0, 0);
	vr = { 'a', 'b', 'c' };
	REQUIRE(v == vr);
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