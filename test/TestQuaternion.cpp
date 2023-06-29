// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "TestGenerators.hpp"

#include <Mathter/Common/Approx.hpp>
#include <Mathter/Quaternion.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstring>
#include <new>


using namespace mathter;
using namespace quat_literals;
using Catch::Approx;


// Expected results based on:
// http://www.andre-gaschler.com/rotationconverter/


TEST_CASE("Quaternion deterministic default initializer", "[Init]") {
	using QuatT = Quaternion<float>;
	alignas(alignof(QuatT)) std::array<uint8_t, sizeof(QuatT)> rawData;
	std::memset(rawData.data(), 0xCC, rawData.size());

	for (auto& v : rawData) {
		REQUIRE(v == uint8_t(0xCC));
	}

	new (rawData.data()) QuatT;

#ifdef NDEBUG
	for (auto& v : rawData) {
		REQUIRE(v == uint8_t(0xCC));
	}
#else
	auto& q = *reinterpret_cast<const QuatT*>(rawData.data());
	REQUIRE(std::isnan(q.x));
	REQUIRE(std::isnan(q.y));
	REQUIRE(std::isnan(q.z));
	REQUIRE(std::isnan(q.w));
#endif
}


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

		QuatT q3(VectorT<3>{ 2, 3, 4 });
		REQUIRE(q3.w == 0);
		REQUIRE(q3.x == 2);
		REQUIRE(q3.y == 3);
		REQUIRE(q3.z == 4);
	}
}


TEST_CASE_VEC_VARIANT("Quaternion - AxisAngle", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q = RotationAxisAngle(Normalize(VectorT<3>{ 1, 2, 3 }), 0.83f);
		QuatT qexp = { 0.9151163f, 0.107757f, 0.2155141f, 0.3232711f };
		REQUIRE(ApproxVec(q) == qexp);
	}
}

TEST_CASE_VEC_VARIANT("Quaternion - TriAxis", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q = RotationAxis3<1, 2, 0>(1.0f, 0.8f, 1.2f);
		QuatT qexp = QuatT(RotationAxisAngle(VectorT<3>{ 1, 0, 0 }, 1.2f))
					 * QuatT(RotationAxisAngle(VectorT<3>{ 0, 0, 1 }, 0.8f))
					 * QuatT(RotationAxisAngle(VectorT<3>{ 0, 1, 0 }, 1.0f));
		REQUIRE(ApproxVec(q) == qexp);
	}
}

TEST_CASE_VEC_VARIANT("Quaternion - QueryAxisAngle", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> axis = { 1, 2, 3 };
		axis = Normalize(axis);
		float angle = 0.83f;

		QuatT q = RotationAxisAngle(axis, angle);

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
		auto m332exp = matrix_representation_cast<Matrix<Type, 3, 3, eMatrixOrder::FOLLOW_VECTOR>>(m331exp);
		Matrix<Type, 3, 4, eMatrixOrder::PRECEDE_VECTOR> m43exp = {
			0.6980989,
			-0.5452151,
			0.4641104,
			0,
			0.6381077,
			0.7677684,
			-0.0578815,
			0,
			-0.3247714,
			0.3365594,
			0.8838842,
			0,
		};
		auto m34exp = matrix_representation_cast<Matrix<Type, 4, 3, eMatrixOrder::FOLLOW_VECTOR>>(m43exp);
		Matrix<Type, 4, 4, eMatrixOrder::PRECEDE_VECTOR> m441exp = {
			0.6980989, -0.5452151, 0.4641104, 0,
			0.6381077, 0.7677684, -0.0578815, 0,
			-0.3247714, 0.3365594, 0.8838842, 0,
			0, 0, 0, 1
		};
		auto m442exp = matrix_representation_cast<Matrix<Type, 4, 4, eMatrixOrder::FOLLOW_VECTOR>>(m441exp);

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
		QuatT q1 = { 1, 2, 3, 4 };
		QuatT q2 = { 4, 5, 6, 3 };
		QuatT q3 = q1 + q2;
		QuatT q4 = q1 - q2;
		QuatT q3exp = { 5, 7, 9, 7 };
		QuatT q4exp = { -3, -3, -3, 1 };

		REQUIRE(ApproxVec(q3exp) == q3);
		REQUIRE(ApproxVec(q4exp) == q4);
	}
}



TEST_CASE_VEC_VARIANT("Quaternion - Product", "[Quaternion]", TypeCases<double>, PackingCases<false>) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q1 = { 1, 2, 3, 4 };
		QuatT q2 = { 4, 5, 6, 3 };
		QuatT q3 = q1 * q2;
		QuatT q3exp = -36 + -2_i + 32_j + 16_k; // dont use this notation, I wrote it just for fun

		REQUIRE(ApproxVec(q3exp) == q3);
	}
}


TEST_CASE_VEC_VARIANT("Quaternion - VectorRotation", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q = RotationAxisAngle(Normalize(VectorT<3>{ 1, 2, 3 }), 0.83f);
		Matrix<Type, 3, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, Packed> M = RotationAxisAngle(Normalize(VectorT<3>{ 1, 2, 3 }), 0.83f);

		VectorT<3> v = { 3, 2, 6 };

		auto v1 = v * q;
		auto v2 = v * M;

		REQUIRE(ApproxVec(v1) == v2);
	}
}


TEST_CASE_VEC_VARIANT("Quaternion - Chaining", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> axis1 = { 1, 2, 3 };
		VectorT<3> axis2 = { 3, 1, 2 };
		axis1 = Normalize(axis1);
		axis2 = Normalize(axis2);
		float angle1 = 0.83f;
		float angle2 = 1.92f;

		QuatT q1 = RotationAxisAngle(axis1, angle1);
		QuatT q2 = RotationAxisAngle(axis2, angle2);
		Matrix<Type, 3, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, Packed> M1 = RotationAxisAngle(axis1, angle1);
		Matrix<Type, 3, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, Packed> M2 = RotationAxisAngle(axis2, angle2);

		VectorT<3> v = { 3, 2, 6 };

		auto v1 = v * (q2 * q1);
		auto v2 = v * (M1 * M2);

		REQUIRE(ApproxVec(v1) == v2);
	}
}



TEST_CASE_VEC_VARIANT("Quaternion - ExpLog", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q(1.0f, 2.0f, 0.5f, -0.7f);

		QuatT p = Exp(Log(q));

		REQUIRE(ApproxVec(q) == p);
	}
}


TEST_CASE_VEC_VARIANT("Quaternion - Pow", "[Quaternion]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		QuatT q(1.0f, 2.0f, 0.5f, -0.7f);

		QuatT p = Pow(q, Type(3));
		QuatT pexp = q * q * q;

		REQUIRE(ApproxVec(p) == pexp);
	}
}