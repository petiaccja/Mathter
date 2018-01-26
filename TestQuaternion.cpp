//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)

#include <Catch2/catch.hpp>

#include "Mathter/Quaternion.hpp"


using namespace mathter;

// Expected results based on:
// http://www.andre-gaschler.com/rotationconverter/


TEST_CASE("Quaternion - Ctor", "[Quaternion]") {
	Quaternion<float> q1(1, 2, 3, 4);
	REQUIRE(q1.w == 1);
	REQUIRE(q1.x == 2);
	REQUIRE(q1.y == 3);
	REQUIRE(q1.z == 4);

	Quaternion<float> q2(1.f, Vector<float, 3>{ 2,3,4 });
	REQUIRE(q2.w == 1);
	REQUIRE(q2.x == 2);
	REQUIRE(q2.y == 3);
	REQUIRE(q2.z == 4);

	Quaternion<float> q3(Vector<float, 3>{2, 3, 4});
	REQUIRE(q3.w == 0);
	REQUIRE(q3.x == 2);
	REQUIRE(q3.y == 3);
	REQUIRE(q3.z == 4);
}


TEST_CASE("Quaternion - AxisAngle", "[Quaternion]") {
	Quaternion<float> q = Quaternion<float>::AxisAngle(Vector<float, 3>{ 1,2,3 }.Normalized(), 0.83f);
	Quaternion<float> qexp = { 0.9151163f, 0.107757f, 0.2155141f, 0.3232711f }; 
	REQUIRE(q.Approx() == qexp);
}

TEST_CASE("Quaternion - QueryAxisAngle", "[Quaternion]") {
	Vector<float, 3> axis = { 1,2,3 };
	axis.Normalize();
	float angle = 0.83f;

	Quaternion<float> q = Quaternion<float>::AxisAngle(axis, angle);

	REQUIRE(axis.Approx() == q.Axis());
	REQUIRE(Approx(angle) == q.Angle());

	q = { 1, 0, 0, 0 };
	axis = { 1, 0, 0 };
	auto qaxis = q.Axis();
	REQUIRE(axis.Approx() == q.Axis());
	REQUIRE(Approx(0.0f) == q.Angle());
}



TEST_CASE("Quaternion - ToMatrix", "[Quaternion]") {
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

	REQUIRE(m331.Approx() == m331exp);
	REQUIRE(m332.Approx() == m332exp);
	REQUIRE(m34. Approx() == m34exp);
	REQUIRE(m43. Approx() == m43exp);
	REQUIRE(m441.Approx() == m441exp);
	REQUIRE(m441.Approx() == m441exp);
}


// Only works if ToMatrix works.
TEST_CASE("Quaternion - FromMatrix", "[Quaternion]") {
	Quaternion<float> q = { 0.9151163f, 0.107757f, 0.2155141f, 0.3232711f };
	Matrix<float, 3, 3, eMatrixOrder::PRECEDE_VECTOR> m331 = (decltype(m331))q;
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR> m332 = (decltype(m332))q;
	Matrix<float, 3, 4, eMatrixOrder::PRECEDE_VECTOR> m43 = (decltype(m43))q;
	Matrix<float, 4, 3, eMatrixOrder::FOLLOW_VECTOR> m34 = (decltype(m34))q;
	Matrix<float, 4, 4, eMatrixOrder::PRECEDE_VECTOR> m441 = (decltype(m441))q;
	Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR> m442 = (decltype(m442))q;

	REQUIRE(q.Approx() == Quaternion<float>(m331));
	REQUIRE(q.Approx() == Quaternion<float>(m332));
	REQUIRE(q.Approx() == Quaternion<float>(m43));
	REQUIRE(q.Approx() == Quaternion<float>(m34));
	REQUIRE(q.Approx() == Quaternion<float>(m441));
	REQUIRE(q.Approx() == Quaternion<float>(m442));
}


TEST_CASE("Quaternion - AddSub", "[Quaternion]") {
	Quaternion<float> q1 = { 1,2,3,4 };
	Quaternion<float> q2 = { 4,5,6,3 };
	Quaternion<float> q3 = q1 + q2;
	Quaternion<float> q4 = q1 - q2;
	Quaternion<float> q3exp = { 5,7,9,7 };
	Quaternion<float> q4exp = { -3,-3,-3,1 };

	REQUIRE(q3exp.Approx() == q3);
	REQUIRE(q4exp.Approx() == q4);
}



TEST_CASE("Quaternion - Product", "[Quaternion]") {
	using namespace quat_literals;

	Quaternion<float> q1 = { 1,2,3,4 };
	Quaternion<float> q2 = { 4,5,6,3 };
	Quaternion<float> q3 = q1 * q2;
	Quaternion<float> q3exp = -36 + -2_i + 32_j + 16_k; // dont use this notation, I wrote it just for fun

	REQUIRE(q3exp.Approx() == q3);
}


TEST_CASE("Quaternion - VectorRotation", "[Quaternion]") {
	Quaternion<float> q = Quaternion<float>::AxisAngle(Vector<float, 3>{ 1, 2, 3 }.Normalized(), 0.83f);
	Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR> M = Matrix<float, 3, 3>::RotationAxisAngle(Vector<float, 3>{ 1, 2, 3 }.Normalized(), 0.83f);

	Vector<float, 3> v = { 3,2,6 };

	auto v1 = v*q;
	auto v2 = v*M;

	REQUIRE(v1.Approx() == v2);
}


TEST_CASE("Quaternion - Chaining", "[Quaternion]") {
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

	REQUIRE(v1.Approx() == v2);
}



TEST_CASE("Quaternion - ExpLog", "[Quaternion]") {
	Quaternion<float> q(1.0f, 2.0f, 0.5f, -0.7f);

	Quaternion<float> p = Quaternion<float>::Exp(Quaternion<float>::Log(q));
	
	REQUIRE(q.Approx() == p);
}


TEST_CASE("Quaternion - Pow", "[Quaternion]") {
	Quaternion<float> q(1.0f, 2.0f, 0.5f, -0.7f);

	Quaternion<float> p = Quaternion<float>::Pow(q, 3);
	Quaternion<float> pexp = q*q*q;

	REQUIRE(p.Approx() == pexp);
}