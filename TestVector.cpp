//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)
#include <gtest\gtest.h>

#include "Mathter\Vector.hpp"

using namespace mathter;



TEST(Vector, CtorAll) {
	Vector<float, 3> v(5.f);

	ASSERT_EQ(v(0), 5.0f);
	ASSERT_EQ(v(1), 5.0f);
	ASSERT_EQ(v(2), 5.0f);

	Vector<double, 5> w(5.0);

	ASSERT_EQ(w(0), 5.0);
	ASSERT_EQ(w(1), 5.0);
	ASSERT_EQ(w(2), 5.0);
	ASSERT_EQ(w(3), 5.0);
	ASSERT_EQ(w(4), 5.0);
}

TEST(Vector, CtorConcat) {
	Vector<double, 2> arg1(1, 2);
	Vector<float, 2> arg3(4, 5);

	Vector<float, 3> a(arg1, 3);
	Vector<float, 3> b(1, 2, 3);
	ASSERT_TRUE(a == b);

	Vector<double, 5> c(arg1, 3, arg3);
	Vector<double, 5> d(1, 2, 3, 4, 5);
	ASSERT_TRUE(c == d);
}


TEST(Vector, VectorAdd) {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> c(5, 7, 9);

	ASSERT_TRUE(a + b == c);

	Vector<double, 5> d(1, 2, 3, 4, 5);
	Vector<double, 5> e(4, 5, 6, 7, 8);
	Vector<double, 5> f(5, 7, 9, 11, 13);

	ASSERT_TRUE(d + e == f);
}


TEST(Vector, VectorSub) {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> c(-3, -3, -3);

	ASSERT_TRUE(a - b == c);

	Vector<double, 5> d(1, 2, 3, 4, 5);
	Vector<double, 5> e(4, 5, 6, 7, 8);
	Vector<double, 5> f(-3, -3, -3, -3, -3);

	ASSERT_TRUE(d - e == f);
}



TEST(Vector, VectorMultiply) {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> c(4, 10, 18);

	ASSERT_TRUE(a*b == c);

	Vector<double, 5> d(1, 2, 3, 4, 5);
	Vector<double, 5> e(4, 5, 6, 7, 8);
	Vector<double, 5> f(4, 10, 18, 28, 40);

	ASSERT_TRUE(d*e == f);
}


TEST(Vector, VectorDiv) {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> c(0.25f, 0.4f, 0.5f);

	ASSERT_TRUE(a / b == c);

	Vector<double, 5> d(2, 6, 6, 12, 10);
	Vector<double, 5> e(1, 2, 3, 4, 5);
	Vector<double, 5> f(2, 3, 2, 3, 2);

	ASSERT_TRUE(d / e == f);
}


TEST(Vector, Dot) {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	float r1 = Vector<float, 3>::Dot(a, b);

	ASSERT_EQ(r1, 32);

	Vector<double, 5> c(1, 2, 3, 2, 1);
	Vector<double, 5> d(4, 5, 6, 5, 4);
	double r2 = Vector<double, 5>::Dot(c, d);
	ASSERT_EQ(r2, 46);
}


TEST(Vector, Cross) {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> r = Vector<float, 3>::Cross(a, b);
	Vector<float, 3> rexp(-3, 6, -3);

	ASSERT_TRUE(r == rexp);
}


TEST(Vector, CrossND) {
	// Simple 3D cross product
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> r = Vector<float, 3>::Cross(a, b);
	Vector<float, 3> rexp(-3, 6, -3);

	ASSERT_TRUE(r == rexp);

	// 2D cross product, that is, rotate by 90 degree
	Vector<float, 2> a2(1, 2);
	Vector<float, 2> r2 = a2.Cross(a2);
	Vector<float, 2> r2exp(-2, 1);

	ASSERT_TRUE(r2.AlmostEqual(r2exp));

	// 4D cross product
	Vector<float, 4> a4(1, 2, 3, 4);
	Vector<float, 4> b4(4, 2, 6, 3);
	Vector<float, 4> c4(3, 6, 4, -9);
	Vector<float, 4> r4 = Cross(a4, b4, c4);

	float dotprod = abs(Dot(a4, r4)) + abs(Dot(b4, r4)) + abs(Dot(c4, r4));
	ASSERT_TRUE(dotprod < 1e-5f);
}


TEST(Vector, Swizzle) {
	Vector<int, 3> v1 = { 1,2,3 };
	Vector<int, 6> v2 = { v1.zx, v1.yzyx };
	Vector<int, 6> v3 = v1.zyx | v1.zyx;
	Vector<int, 6> v2exp = { 3,1,2,3,2,1 };
	Vector<int, 6> v3exp = { 3,2,1,3,2,1 };

	ASSERT_EQ(v2, v2exp);
	ASSERT_EQ(v3, v3exp);

	Vector<int, 4> v4{ 1,2,3,4 };
	v4.yxwz = v4.wzyx; // wzxy=4321 -> v4=3412
	Vector<int, 4> v4exp = { 3, 4, 1, 2 };

	ASSERT_TRUE(v4 == v4exp);

	v4 = { 1,2,3,4 };
	v4 = v4.xxzz;
	v4exp = { 1,1,3,3 };
	ASSERT_TRUE(v4 == v4exp);

	v4 = v1.zyx | 1.0f;
	v4 = 1.0f | v1.zyx;
}