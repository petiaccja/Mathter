//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)
#include <gtest/gtest.h>

#include "Mathter/Quaternion.hpp"


using namespace mathter;

// Expected results based on:
// http://www.andre-gaschler.com/rotationconverter/


TEST(Quaternion, Ctor) {
	Quaternion<float> q1(1, 2, 3, 4);
	ASSERT_EQ(q1.w, 1);
	ASSERT_EQ(q1.x, 2);
	ASSERT_EQ(q1.y, 3);
	ASSERT_EQ(q1.z, 4);

	Quaternion<float> q2(1.f, Vector<float, 3>{ 2,3,4 });
	ASSERT_EQ(q2.w, 1);
	ASSERT_EQ(q2.x, 2);
	ASSERT_EQ(q2.y, 3);
	ASSERT_EQ(q2.z, 4);

	Quaternion<float> q3(Vector<float, 3>{2, 3, 4});
	ASSERT_EQ(q3.w, 0);
	ASSERT_EQ(q3.x, 2);
	ASSERT_EQ(q3.y, 3);
	ASSERT_EQ(q3.z, 4);
}


TEST(Quaternion, AxisAngle) {
	Quaternion<float> q = Quaternion<float>::AxisAngle(Vector<float, 3>{ 1,2,3 }.Normalized(), 0.83f);
	Quaternion<float> qexp = { 0.9151163f, 0.107757f, 0.2155141f, 0.3232711f }; 
	ASSERT_TRUE(q.AlmostEqual(qexp));
}

TEST(Quaternion, QueryAxisAngle) {
	Vector<float, 3> axis = { 1,2,3 };
	axis.Normalize();
	float angle = 0.83f;

	Quaternion<float> q = Quaternion<float>::AxisAngle(axis, angle);

	ASSERT_TRUE(axis.AlmostEqual(q.Axis()));
	ASSERT_FLOAT_EQ(angle, q.Angle());

	q = { 1, 0, 0, 0 };
	axis = { 1, 0, 0 };
	auto qaxis = q.Axis();
	ASSERT_TRUE(axis.AlmostEqual(q.Axis()));
	ASSERT_FLOAT_EQ(0.0f, q.Angle());
}



TEST(Quaternion, ToMatrix) {
	Quaternion<float> q = { 0.9151163f, 0.107757f, 0.2155141f, 0.3232711f };
	Matrix<float, 3, 3, eMatrixOrder::PRECEDE_VECTOR> m331 = (decltype(m331))q;
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR> m332 = (decltype(m332))q;
	Matrix<float, 3, 4, eMatrixOrder::PRECEDE_VECTOR> m43 = (decltype(m43))q;
	Matrix<float, 4, 3, eMatrixOrder::FOLLOW_VECTOR> m34 = (decltype(m34))q;
	Matrix<float, 4, 4, eMatrixOrder::PRECEDE_VECTOR> m441 = (decltype(m441))q;
	Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR> m442 = (decltype(m442))q;

	Matrix<float, 3, 3, eMatrixOrder::PRECEDE_VECTOR> m331exp = {
		0.6980989, -0.5452151, 0.4641104,
		0.6381077, 0.7677684, -0.0578815,
		-0.3247714, 0.3365594, 0.8838842
	};
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR> m332exp = m331exp;
	Matrix<float, 3, 4, eMatrixOrder::PRECEDE_VECTOR> m43exp = {
		0.6980989, -0.5452151, 0.4641104,	0,
		0.6381077, 0.7677684, -0.0578815,	0,
		-0.3247714, 0.3365594, 0.8838842,	0,
	};
	Matrix<float, 4, 3, eMatrixOrder::FOLLOW_VECTOR> m34exp = m43exp;
	Matrix<float, 4, 4, eMatrixOrder::PRECEDE_VECTOR> m441exp = {
		0.6980989, -0.5452151, 0.4641104,	0,
		0.6381077, 0.7677684, -0.0578815,	0,
		-0.3247714, 0.3365594, 0.8838842,	0,
		0,			0,			0,			1
	};
	Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR> m442exp = m441exp;

	ASSERT_TRUE(m331.AlmostEqual(m331exp));
	ASSERT_TRUE(m332.AlmostEqual(m332exp));
	ASSERT_TRUE(m34.AlmostEqual(m34exp));
	ASSERT_TRUE(m43.AlmostEqual(m43exp));
	ASSERT_TRUE(m441.AlmostEqual(m441exp));
	ASSERT_TRUE(m441.AlmostEqual(m441exp));
}


// Only works if ToMatrix works.
TEST(Quaternion, FromMatrix) {
	Quaternion<float> q = { 0.9151163f, 0.107757f, 0.2155141f, 0.3232711f };
	Matrix<float, 3, 3, eMatrixOrder::PRECEDE_VECTOR> m331 = (decltype(m331))q;
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR> m332 = (decltype(m332))q;
	Matrix<float, 3, 4, eMatrixOrder::PRECEDE_VECTOR> m43 = (decltype(m43))q;
	Matrix<float, 4, 3, eMatrixOrder::FOLLOW_VECTOR> m34 = (decltype(m34))q;
	Matrix<float, 4, 4, eMatrixOrder::PRECEDE_VECTOR> m441 = (decltype(m441))q;
	Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR> m442 = (decltype(m442))q;

	ASSERT_TRUE(q.AlmostEqual(Quaternion<float>(m331)));
	ASSERT_TRUE(q.AlmostEqual(Quaternion<float>(m332)));
	ASSERT_TRUE(q.AlmostEqual(Quaternion<float>(m43)));
	ASSERT_TRUE(q.AlmostEqual(Quaternion<float>(m34)));
	ASSERT_TRUE(q.AlmostEqual(Quaternion<float>(m441)));
	ASSERT_TRUE(q.AlmostEqual(Quaternion<float>(m442)));
}


TEST(Quaternion, AddSub) {
	Quaternion<float> q1 = { 1,2,3,4 };
	Quaternion<float> q2 = { 4,5,6,3 };
	Quaternion<float> q3 = q1 + q2;
	Quaternion<float> q4 = q1 - q2;
	Quaternion<float> q3exp = { 5,7,9,7 };
	Quaternion<float> q4exp = { -3,-3,-3,1 };

	ASSERT_TRUE(q3exp.AlmostEqual(q3));
	ASSERT_TRUE(q4exp.AlmostEqual(q4));
}



TEST(Quaternion, Product) {
	using namespace quat_literals;

	Quaternion<float> q1 = { 1,2,3,4 };
	Quaternion<float> q2 = { 4,5,6,3 };
	Quaternion<float> q3 = q1 * q2;
	Quaternion<float> q3exp = -36 + -2_i + 32_j + 16_k; // dont use this notation, I wrote it just for fun

	ASSERT_TRUE(q3exp.AlmostEqual(q3));
}


TEST(Quaternion, VectorRotation) {
	Quaternion<float> q = Quaternion<float>::AxisAngle(Vector<float, 3>{ 1, 2, 3 }.Normalized(), 0.83f);
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR> M = Matrix<float, 3, 3>::RotationAxisAngle(Vector<float, 3>{ 1, 2, 3 }.Normalized(), 0.83f);

	Vector<float, 3> v = { 3,2,6 };

	auto v1 = v*q;
	auto v2 = v*M;

	ASSERT_TRUE(v1.AlmostEqual(v2));
}


TEST(Quaternion, Chaining) {
	Vector<float, 3> axis1 = { 1,2,3 };
	Vector<float, 3> axis2 = { 3,1,2 };
	axis1.Normalize();
	axis2.Normalize();
	float angle1 = 0.83f;
	float angle2 = 1.92f;

	Quaternion<float> q1 = Quaternion<float>::AxisAngle(axis1, angle1);
	Quaternion<float> q2 = Quaternion<float>::AxisAngle(axis2, angle2);
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR> M1 = Matrix<float, 3, 3>::RotationAxisAngle(axis1, angle1);
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR> M2 = Matrix<float, 3, 3>::RotationAxisAngle(axis2, angle2);

	Vector<float, 3> v = { 3,2,6 };

	auto v1 = v*(q2*q1);
	auto v2 = v*(M1*M2);

	ASSERT_TRUE(v1.AlmostEqual(v2));
}



TEST(Quaternion, ExpLog) {
	Quaternion<float> q(1.0f, 2.0f, 0.5f, -0.7f);

	Quaternion<float> p = Quaternion<float>::Exp(Quaternion<float>::Log(q));
	
	ASSERT_TRUE(q.AlmostEqual(p));
}


TEST(Quaternion, Pow) {
	Quaternion<float> q(1.0f, 2.0f, 0.5f, -0.7f);

	Quaternion<float> p = Quaternion<float>::Pow(q, 3);
	Quaternion<float> pexp = q*q*q;

	ASSERT_TRUE(p.AlmostEqual(pexp));
}