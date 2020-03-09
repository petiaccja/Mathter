//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)

#include <Catch2/catch.hpp>

#include "Mathter/Quaternion.hpp"
#include "Mathter/Approx.hpp"
#include "TestGenerators.hpp"


using namespace mathter;
using namespace quat_literals;

// Expected results based on:
// http://www.andre-gaschler.com/rotationconverter/


TEST_CASE_VEC_VARIANT("Quaternion - Ctor", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q1(1, 2, 3, 4);
		REQUIRE(q1.w == 1);
		REQUIRE(q1.x == 2);
		REQUIRE(q1.y == 3);
		REQUIRE(q1.z == 4);

		QuatT q2(1.f, VectorT<3>{ 2, 3, 4 });
		REQUIRE(q2.w == 1);
		REQUIRE(q2.x == 2);
		REQUIRE(q2.y == 3);
		REQUIRE(q2.z == 4);

		QuatT q3(VectorT<3>{2, 3, 4});
		REQUIRE(q3.w == 0);
		REQUIRE(q3.x == 2);
		REQUIRE(q3.y == 3);
		REQUIRE(q3.z == 4);
	}
}


TEST_CASE_VEC_VARIANT("Quaternion - AxisAngle", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q = QuatT::AxisAngle(VectorT<3>{ 1, 2, 3 }.Normalized(), 0.83f);
		QuatT qexp = { 0.9151163f, 0.107757f, 0.2155141f, 0.3232711f };
		REQUIRE(ApproxVec(q) == qexp);
	}
}

TEST_CASE_VEC_VARIANT("Quaternion - QueryAxisAngle", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> axis = { 1,2,3 };
		axis.Normalize();
		float angle = 0.83f;

		QuatT q = QuatT::AxisAngle(axis, angle);

		REQUIRE(ApproxVec(axis) == q.Axis());
		REQUIRE(Approx(angle) == q.Angle());

		q = { 1, 0, 0, 0 };
		axis = { 1, 0, 0 };
		auto qaxis = q.Axis();
		REQUIRE(ApproxVec(axis) == q.Axis());
		REQUIRE(Approx(0.0f) == q.Angle());
	}
}



TEST_CASE_VEC_VARIANT("Quaternion - ToMatrix", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q = { 0.9151163f, 0.107757f, 0.2155141f, 0.3232711f };
		Matrix<Type, 3, 3, eMatrixOrder::PRECEDE_VECTOR> m331 = (decltype(m331))q;
		Matrix<Type, 3, 3, eMatrixOrder::FOLLOW_VECTOR> m332 = (decltype(m332))q;
		Matrix<Type, 3, 4, eMatrixOrder::PRECEDE_VECTOR> m43 = (decltype(m43))q;
		Matrix<Type, 4, 3, eMatrixOrder::FOLLOW_VECTOR> m34 = (decltype(m34))q;
		Matrix<Type, 4, 4, eMatrixOrder::PRECEDE_VECTOR> m441 = (decltype(m441))q;
		Matrix<Type, 4, 4, eMatrixOrder::FOLLOW_VECTOR> m442 = (decltype(m442))q;

		Matrix<Type, 3, 3, eMatrixOrder::PRECEDE_VECTOR> m331exp = {
			0.6980989, -0.5452151, 0.4641104,
			0.6381077, 0.7677684, -0.0578815,
			-0.3247714, 0.3365594, 0.8838842
		};
		Matrix<Type, 3, 3, eMatrixOrder::FOLLOW_VECTOR> m332exp = m331exp;
		Matrix<Type, 3, 4, eMatrixOrder::PRECEDE_VECTOR> m43exp = {
			0.6980989, -0.5452151, 0.4641104,	0,
			0.6381077, 0.7677684, -0.0578815,	0,
			-0.3247714, 0.3365594, 0.8838842,	0,
		};
		Matrix<Type, 4, 3, eMatrixOrder::FOLLOW_VECTOR> m34exp = m43exp;
		Matrix<Type, 4, 4, eMatrixOrder::PRECEDE_VECTOR> m441exp = {
			0.6980989, -0.5452151, 0.4641104,	0,
			0.6381077, 0.7677684, -0.0578815,	0,
			-0.3247714, 0.3365594, 0.8838842,	0,
			0,			0,			0,			1
		};
		Matrix<Type, 4, 4, eMatrixOrder::FOLLOW_VECTOR> m442exp = m441exp;

		REQUIRE(ApproxVec(m331) == m331exp);
		REQUIRE(ApproxVec(m332) == m332exp);
		REQUIRE(ApproxVec(m34) == m34exp);
		REQUIRE(ApproxVec(m43) == m43exp);
		REQUIRE(ApproxVec(m441) == m441exp);
		REQUIRE(ApproxVec(m441) == m441exp);
	}
}


// Only works if ToMatrix works.
TEST_CASE_VEC_VARIANT("Quaternion - FromMatrix", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q = { 0.9151163f, 0.107757f, 0.2155141f, 0.3232711f };
		Matrix<Type, 3, 3, eMatrixOrder::PRECEDE_VECTOR> m331 = (decltype(m331))q;
		Matrix<Type, 3, 3, eMatrixOrder::FOLLOW_VECTOR> m332 = (decltype(m332))q;
		Matrix<Type, 3, 4, eMatrixOrder::PRECEDE_VECTOR> m43 = (decltype(m43))q;
		Matrix<Type, 4, 3, eMatrixOrder::FOLLOW_VECTOR> m34 = (decltype(m34))q;
		Matrix<Type, 4, 4, eMatrixOrder::PRECEDE_VECTOR> m441 = (decltype(m441))q;
		Matrix<Type, 4, 4, eMatrixOrder::FOLLOW_VECTOR> m442 = (decltype(m442))q;

		REQUIRE(ApproxVec(q) == QuatT(m331));
		REQUIRE(ApproxVec(q) == QuatT(m332));
		REQUIRE(ApproxVec(q) == QuatT(m43));
		REQUIRE(ApproxVec(q) == QuatT(m34));
		REQUIRE(ApproxVec(q) == QuatT(m441));
		REQUIRE(ApproxVec(q) == QuatT(m442));
	}
}


TEST_CASE_VEC_VARIANT("Quaternion - AddSub", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q1 = { 1,2,3,4 };
		QuatT q2 = { 4,5,6,3 };
		QuatT q3 = q1 + q2;
		QuatT q4 = q1 - q2;
		QuatT q3exp = { 5,7,9,7 };
		QuatT q4exp = { -3,-3,-3,1 };

		REQUIRE(ApproxVec(q3exp) == q3);
		REQUIRE(ApproxVec(q4exp) == q4);
	}
}



TEST_CASE_VEC_VARIANT("Quaternion - Product", "[Quaternion]", TypeCases<double>, PackingCases<false>) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q1 = { 1,2,3,4 };
		QuatT q2 = { 4,5,6,3 };
		QuatT q3 = q1 * q2;
		QuatT q3exp = -36 + -2_i + 32_j + 16_k; // dont use this notation, I wrote it just for fun

		REQUIRE(ApproxVec(q3exp) == q3);
	}
}


TEST_CASE_VEC_VARIANT("Quaternion - VectorRotation", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		auto q = QuatT::AxisAngle(VectorT<3>{ 1, 2, 3 }.Normalized(), 0.83f);
		auto M = Matrix<Type, 3, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, Packed>::RotationAxisAngle(VectorT<3>{ 1, 2, 3 }.Normalized(), 0.83f);

		VectorT<3> v = { 3,2,6 };

		auto v1 = v*q;
		auto v2 = v*M;

		REQUIRE(ApproxVec(v1) == v2);
	}
}


TEST_CASE_VEC_VARIANT("Quaternion - Chaining", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> axis1 = { 1,2,3 };
		VectorT<3> axis2 = { 3,1,2 };
		axis1.Normalize();
		axis2.Normalize();
		float angle1 = 0.83f;
		float angle2 = 1.92f;

		QuatT q1 = QuatT::AxisAngle(axis1, angle1);
		QuatT q2 = QuatT::AxisAngle(axis2, angle2);
		auto M1 = Matrix<Type, 3, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, Packed>::RotationAxisAngle(axis1, angle1);
		auto M2 = Matrix<Type, 3, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, Packed>::RotationAxisAngle(axis2, angle2);

		VectorT<3> v = { 3,2,6 };

		auto v1 = v*(q2*q1);
		auto v2 = v*(M1*M2);

		REQUIRE(ApproxVec(v1) == v2);
	}
}



TEST_CASE_VEC_VARIANT("Quaternion - ExpLog", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q(1.0f, 2.0f, 0.5f, -0.7f);

		QuatT p = QuatT::Exp(QuatT::Log(q));

		REQUIRE(ApproxVec(q) == p);
	}
}


TEST_CASE_VEC_VARIANT("Quaternion - Pow", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q(1.0f, 2.0f, 0.5f, -0.7f);

		QuatT p = QuatT::Pow(q, 3);
		QuatT pexp = q*q*q;

		REQUIRE(ApproxVec(p) == pexp);
	}
}