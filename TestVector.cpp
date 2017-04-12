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
	ASSERT_EQ(a, b);

	Vector<double, 5> c(arg1, 3, arg3);
	Vector<double, 5> d(1, 2, 3, 4, 5);
	ASSERT_EQ(c, d);
}


TEST(Vector, VectorAdd) {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> c(5, 7, 9);

	ASSERT_EQ(a + b, c);

	Vector<double, 5> d(1, 2, 3, 4, 5);
	Vector<double, 5> e(4, 5, 6, 7, 8);
	Vector<double, 5> f(5, 7, 9, 11, 13);

	ASSERT_EQ(d + e, f);
}


TEST(Vector, VectorSub) {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> c(-3, -3, -3);

	ASSERT_EQ(a - b, c);

	Vector<double, 5> d(1, 2, 3, 4, 5);
	Vector<double, 5> e(4, 5, 6, 7, 8);
	Vector<double, 5> f(-3, -3, -3, -3, -3);

	ASSERT_EQ(d - e, f);
}



TEST(Vector, VectorMultiply) {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> c(4, 10, 18);

	ASSERT_EQ(a*b, c);

	Vector<double, 5> d(1, 2, 3, 4, 5);
	Vector<double, 5> e(4, 5, 6, 7, 8);
	Vector<double, 5> f(4, 10, 18, 28, 40);

	ASSERT_EQ(d*e, f);
}


TEST(Vector, VectorDiv) {
	Vector<float, 3> a(1, 2, 3);
	Vector<float, 3> b(4, 5, 6);
	Vector<float, 3> c(0.25f, 0.4f, 0.5f);

	ASSERT_EQ(a / b, c);

	Vector<double, 5> d(2, 6, 6, 12, 10);
	Vector<double, 5> e(1, 2, 3, 4, 5);
	Vector<double, 5> f(2, 3, 2, 3, 2);

	ASSERT_EQ(d / e, f);
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

	ASSERT_EQ(r, rexp);
}