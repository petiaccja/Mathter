//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)

#include <Catch2/catch.hpp>

#include "Mathter/Vector.hpp"

using namespace mathter;



TEST_CASE("Vector - CtorAll", "[Vector]") {
	Vector<float, 3> v(5.f);

	REQUIRE(v(0) == 5.0f);
	REQUIRE(v(1) == 5.0f);
	REQUIRE(v(2) == 5.0f);

	Vector<double, 5> w(5.0);

	REQUIRE(w(0) == 5.0);
	REQUIRE(w(1) == 5.0);
	REQUIRE(w(2) == 5.0);
	REQUIRE(w(3) == 5.0);
	REQUIRE(w(4) == 5.0);
}

TEST_CASE("Vector - CtorConcat", "[Vector]") {
	Vector<double, 2> arg1(1, 2);
	Vector<float, 2> arg3(4, 5);

	Vector<float, 3> a(arg1, 3);
	Vector<float, 3> b(1, 2, 3);
	REQUIRE(a.Approx() == b);

	Vector<double, 5> c(arg1, 3, arg3);
	Vector<double, 5> d(1, 2, 3, 4, 5);
	REQUIRE(c == d);
}


TEST_CASE("Vector - VectorAdd", "[Vector]") {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> c(5, 7, 9);

	REQUIRE(a + b == c);

	Vector<double, 5> d(1, 2, 3, 4, 5);
	Vector<double, 5> e(4, 5, 6, 7, 8);
	Vector<double, 5> f(5, 7, 9, 11, 13);

	REQUIRE(d + e == f);
}


TEST_CASE("Vector - VectorSub", "[Vector]") {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> c(-3, -3, -3);

	REQUIRE(a - b == c);

	Vector<double, 5> d(1, 2, 3, 4, 5);
	Vector<double, 5> e(4, 5, 6, 7, 8);
	Vector<double, 5> f(-3, -3, -3, -3, -3);

	REQUIRE(d - e == f);
}



TEST_CASE("Vector - VectorMultiply", "[Vector]") {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> c(4, 10, 18);

	REQUIRE(a*b == c);

	Vector<double, 5> d(1, 2, 3, 4, 5);
	Vector<double, 5> e(4, 5, 6, 7, 8);
	Vector<double, 5> f(4, 10, 18, 28, 40);

	REQUIRE(d*e == f);
}


TEST_CASE("Vector - VectorDiv", "[Vector]") {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> c(0.25f, 0.4f, 0.5f);

	REQUIRE(a / b == c);

	Vector<double, 5> d(2, 6, 6, 12, 10);
	Vector<double, 5> e(1, 2, 3, 4, 5);
	Vector<double, 5> f(2, 3, 2, 3, 2);

	REQUIRE(d / e == f);
}


TEST_CASE("Vector - Dot", "[Vector]") {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	float r1 = Vector<float, 3>::Dot(a, b);

	REQUIRE(r1 == 32);

	Vector<double, 5> c(1, 2, 3, 2, 1);
	Vector<double, 5> d(4, 5, 6, 5, 4);
	double r2 = Vector<double, 5>::Dot(c, d);
	REQUIRE(r2 == 46);
}


TEST_CASE("Vector - Cross", "[Vector]") {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> r = Vector<float, 3>::Cross(a, b);
	Vector<float, 3> rexp(-3, 6, -3);

	REQUIRE(r == rexp);
}


TEST_CASE("Vector - CrossND", "[Vector]") {
	// Simple 3D cross product
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> r = Vector<float, 3>::Cross(a, b);
	Vector<float, 3> rexp(-3, 6, -3);

	REQUIRE(r == rexp);

	// 2D cross product, that is, rotate by 90 degree
	Vector<float, 2> a2(1, 2);
	Vector<float, 2> r2 = a2.Cross(a2);
	Vector<float, 2> r2exp(-2, 1);

	REQUIRE(r2.Approx() == r2exp);

	// 4D cross product
	Vector<float, 4> a4(1, 2, 3, 4);
	Vector<float, 4> b4(4, 2, 6, 3);
	Vector<float, 4> c4(3, 6, 4, -9);
	Vector<float, 4> r4 = Cross(a4, b4, c4);

	float dotprod = std::abs(Dot(a4, r4)) + std::abs(Dot(b4, r4)) + std::abs(Dot(c4, r4));
	REQUIRE(dotprod < 1e-5f);
}


TEST_CASE("Vector - Swizzle", "[Vector]") {
	Vector<int, 3> v1 = { 1,2,3 };
	Vector<int, 6> v2 = { v1.zx, v1.yzyx };
	Vector<int, 6> v3 = v1.zyx | v1.zyx;
	Vector<int, 6> v2exp = { 3,1,2,3,2,1 };
	Vector<int, 6> v3exp = { 3,2,1,3,2,1 };

	REQUIRE(v2 == v2exp);
	REQUIRE(v3 == v3exp);

	Vector<int, 4> v4{ 1,2,3,4 };
	v4.yxwz = v4.wzyx; // wzxy=4321 -> v4=3412
	Vector<int, 4> v4exp = { 3, 4, 1, 2 };

	REQUIRE(v4 == v4exp);

	v4 = { 1,2,3,4 };
	v4 = v4.xxzz;
	v4exp = { 1,1,3,3 };
	REQUIRE(v4 == v4exp);

	v4 = v1.zyx | 1.0f;
	v4 = 1.0f | v1.zyx;
}


TEST_CASE("Vector - TruncateExtend", "[Vector]") {
	Vector<float, 3> v(1,2,3);
	Vector<float, 4> u = (Vector<float, 4>)v;
	Vector<float, 3> v2 = (Vector<float, 3>)u;

	REQUIRE(u == Vector<float, 4>(1, 2, 3, 1));
	REQUIRE(v == v2);
}


TEST_CASE("Vector - IOParse", "[Vector]") {
	Vector<float, 3> parsed;

	std::string successCases[] = {
		"3.14, 2.718, 0.57",
		"  (3.14,2.718  ,  0.57\t  )  ",
		"[3.14\t2.718\t0.57]",
		"{    3.14  2.718, 0.57 } ",
	};
	Vector<float, 3> expected{3.14f, 2.718f, 0.57f};

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