//==============================================================================
// This software is distributed under The Unlicense.
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable : 4244)

#include "../Mathter/Common/Approx.hpp"
#include "../Mathter/Matrix.hpp"
#include "TestGenerators.hpp"

#include <Catch2/catch.hpp>
#include <complex>

using namespace mathter;


TEST_CASE_VARIANT("Matrix - Identity", "[Transform]", TypesAll, OrdersAll, LayoutsAll, PackedAll) {
	Matrix<float, 3, 3> m = Identity();
	Matrix<float, 3, 3> mexp = {
		1, 0, 0,
		0, 1, 0,
		0, 0, 1
	};

	REQUIRE(m == mexp);

	Matrix<float, 3, 5> m5 = Identity();
	Matrix<float, 3, 5> mexp5 = {
		1, 0, 0, 0, 0,
		0, 1, 0, 0, 0,
		0, 0, 1, 0, 0
	};

	REQUIRE(m5 == mexp5);
}

TEST_CASE_VARIANT("Matrix - Zero", "[Transform]", TypesAll, OrdersAll, LayoutsAll, PackedAll) {
	Matrix<float, 3, 4> m = Zero();
	Matrix<float, 3, 4> mexp = {
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0
	};

	REQUIRE(m == mexp);
}



TEST_CASE("Matrix - Rotation2D", "[Matrix]") {
	Matrix<float, 2, 2> m = Rotation(1.f);
	Matrix<float, 3, 3> m3 = Rotation(1.f);
	Matrix<float, 2, 2> mexp = {
		0.54030, 0.84147,
		-0.84147, 0.54030
	};
	Matrix<float, 3, 3> m3exp = {
		0.54030, 0.84147, 0,
		-0.84147, 0.54030, 0,
		0, 0, 1
	};

	REQUIRE(ApproxVec(m) == mexp);
	REQUIRE(ApproxVec(m3) == m3exp);
}


TEST_CASE("Matrix - RotationPrincipal", "[Matrix]") {
	Matrix<float, 3, 3> m1 = RotationX(1.f);
	Matrix<float, 3, 3> mexp = {
		1.000000, 0.000000, 0.000000,
		0.000000, 0.540302, 0.841471,
		0.000000, -0.841471, 0.540302
	};
	REQUIRE(ApproxVec(m1) == mexp);


	Matrix<float, 4, 3> m2 = RotationY(1.f);
	Matrix<float, 4, 3> m2exp = {
		0.540302, 0.000000, -0.841471,
		0.000000, 1.000000, 0.000000,
		0.841471, 0.000000, 0.540302,
		0, 0, 0
	};
	REQUIRE(ApproxVec(m2) == m2exp);


	Matrix<float, 3, 4, eMatrixOrder::PRECEDE_VECTOR> m3 = RotationZ(1.f);
	Matrix<float, 4, 3> m3exp = {
		0.540302, 0.841471, 0.000000,
		-0.841471, 0.540302, 0.000000,
		0.000000, 0.000000, 1.000000,
		0, 0, 0
	};
	REQUIRE(ApproxVec(m3) == Transpose(m3exp));

	Matrix<float, 4, 4> m4 = RotationZ(1.f);
	Matrix<float, 4, 4> m4exp = {
		0.540302, 0.841471, 0.000000, 0,
		-0.841471, 0.540302, 0.000000, 0,
		0.000000, 0.000000, 1.000000, 0,
		0, 0, 0, 1
	};
	REQUIRE(ApproxVec(m4) == m4exp);
}


TEST_CASE("Matrix - RotationTriAxis", "[Matrix]") {
	Matrix<float, 3, 3> m1 = RotationAxis3<0, 1, 1>(1.f, 1.0f, -1.0f);
	Matrix<float, 3, 3> mexp = {
		1.000000, 0.000000, 0.000000,
		0.000000, 0.540302, 0.841471,
		0.000000, -0.841471, 0.540302
	};
	REQUIRE(ApproxVec(m1) == mexp);


	Matrix<float, 4, 3> m2 = RotationAxis3<0, 1, 2>(0.0f, 1.f, 0.0f);
	Matrix<float, 4, 3> m2exp = {
		0.540302, 0.000000, -0.841471,
		0.000000, 1.000000, 0.000000,
		0.841471, 0.000000, 0.540302,
		0, 0, 0
	};
	REQUIRE(ApproxVec(m2) == m2exp);


	Matrix<float, 3, 4, eMatrixOrder::PRECEDE_VECTOR> m3 = RotationAxis3<1, 1, 2>(1.0f, -1.0f, 1.f);
	Matrix<float, 4, 3> m3exp = {
		0.540302, 0.841471, 0.000000,
		-0.841471, 0.540302, 0.000000,
		0.000000, 0.000000, 1.000000,
		0, 0, 0
	};
	REQUIRE(ApproxVec(m3) == Transpose(m3exp));

	Matrix<float, 4, 4> m4 = RotationAxis3<0, 0, 2>(-1.0f, 1.0f, 1.f);
	Matrix<float, 4, 4> m4exp = {
		0.540302, 0.841471, 0.000000, 0,
		-0.841471, 0.540302, 0.000000, 0,
		0.000000, 0.000000, 1.000000, 0,
		0, 0, 0, 1
	};
	REQUIRE(ApproxVec(m4) == m4exp);
}


TEST_CASE("Matrix - RotationAxisAngle", "[Matrix]") {
	Matrix<float, 3, 3> m = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
	Matrix<float, 3, 3> mexp = {
		0.573138, 0.740349, -0.351279, -0.609007, 0.671645, 0.421906, 0.548292, -0.027879, 0.835822
	};
	REQUIRE(ApproxVec(m) == mexp);

	Matrix<float, 4, 4> m4 = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
	Matrix<float, 4, 4> m4exp = {
		0.573138, 0.740349, -0.351279, 0,
		-0.609007, 0.671645, 0.421906, 0,
		0.548292, -0.027879, 0.835822, 0,
		0, 0, 0, 1
	};
	REQUIRE(ApproxVec(m4) == m4exp);
}


TEST_CASE("Matrix - Scale", "[Matrix]") {
	Matrix<float, 5, 5> m = Scale(1, 2, 3, 4, 5);
	Vector<float, 5> v(2, 6, 3, 7, 5);
	Matrix<float, 3, 3> m3 = Scale(Vector<float, 3>{ 1, 2, 3 });

	auto vt1 = v * Vector<float, 5>{ 1, 2, 3, 4, 5 };
	auto vt2 = v * m;

	REQUIRE(vt1 == vt2);
}


TEST_CASE("Matrix - Translation", "[Matrix]") {
	Matrix<float, 3, 3> m33 = Translation(1, 2);
	Matrix<float, 6, 5> m = Translation(Vector<float, 5>{ 1, 2, 3, 4, 5 });
	Vector<float, 5> v(1, 2, 3, 4, 5);
	v = v * m;
	Vector<float, 5> vexp(2, 4, 6, 8, 10);
	REQUIRE(v == vexp);

	Matrix<float, 3, 3> m2 = Translation(1, 2);
	Matrix<float, 3, 3> m2exp = {
		1, 0, 0,
		0, 1, 0,
		1, 2, 1
	};
	REQUIRE(m2 == m2exp);

	Matrix<float, 3, 2> m3 = Translation(Vector<float, 2>(1, 2));
	Matrix<float, 3, 2> m3exp = {
		1, 0,
		0, 1,
		1, 2
	};
	REQUIRE(m3 == m3exp);

	Matrix<float, 2, 3, eMatrixOrder::PRECEDE_VECTOR> m4 = Translation(Vector<float, 2>(1, 2));
	Matrix<float, 2, 3, eMatrixOrder::PRECEDE_VECTOR> m4exp = {
		1, 0, 1,
		0, 1, 2
	};
	REQUIRE(m4 == m4exp);
}


TEST_CASE("Matrix - Perspective", "[Matrix]") {
	Vector<float, 4> worldFrustum[2] = {
		{ -0.25f, -0.140625f, 0.5f, 1 },
		{ 5.0f, 2.8125f, 10.f, 1 }
	};
	Vector<float, 4> ndcFrustum[2];

	// Z forward
	Matrix<float, 4, 4> m = Perspective(53.13010235f / 180.f * 3.1415926f, 16.f / 9.f, 0.5f, 10.f, 0.f, 1.f);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	REQUIRE(ApproxVec(ndcFrustum[0]) == Vector<float, 4>{ -1, -1, 0, 1 });
	REQUIRE(ApproxVec(ndcFrustum[1]) == Vector<float, 4>{ 1, 1, 1, 1 });

	// Z backward in NDC
	m = Perspective(53.13010235f / 180.f * 3.1415926f, 16.f / 9.f, 0.5f, 10.f, 1.f, -1.f);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	REQUIRE(ApproxVec(ndcFrustum[0]) == Vector<float, 4>{ -1, -1, 1, 1 });
	REQUIRE(ApproxVec(ndcFrustum[1]) == Vector<float, 4>{ 1, 1, -1, 1 });

	// Z backward in world
	m = Perspective(53.13010235f / 180.f * 3.1415926f, 16.f / 9.f, -0.5f, -10.f, 0.f, 1.f);
	worldFrustum[0].z *= -1;
	worldFrustum[1].z *= -1;
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	REQUIRE(ApproxVec(ndcFrustum[0]) == Vector<float, 4>{ -1, -1, 0, 1 });
	REQUIRE(ApproxVec(ndcFrustum[1]) == Vector<float, 4>{ 1, 1, 1, 1 });

	// Z backward in world && NDC
	m = Perspective(53.13010235f / 180.f * 3.1415926f, 16.f / 9.f, -0.5f, -10.f, 1.f, -1.f);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;
	ndcFrustum[0] /= ndcFrustum[0].w;
	ndcFrustum[1] /= ndcFrustum[1].w;

	REQUIRE(ApproxVec(ndcFrustum[0]) == Vector<float, 4>{ -1, -1, 1, 1 });
	REQUIRE(ApproxVec(ndcFrustum[1]) == Vector<float, 4>{ 1, 1, -1, 1 });
}

TEST_CASE("Matrix - Orthographic", "[Matrix]") {
	Vector<float, 3> worldFrustum[2] = {
		{ -0.25f, -0.44444444f, 0.5f },
		{ 5.0f, 8.8888888f, 10.f }
	};
	Vector<float, 3> ndcFrustum[2];

	// Z forward
	Matrix<float, 4, 4> m = Orthographic(worldFrustum[0], worldFrustum[1], 0.f, 1.f);
	ndcFrustum[0] = worldFrustum[0] * m;
	ndcFrustum[1] = worldFrustum[1] * m;

	REQUIRE(ApproxVec(ndcFrustum[0]) == Vector<float, 3>{ -1, -1, 0 });
	REQUIRE(ApproxVec(ndcFrustum[1]) == Vector<float, 3>{ 1, 1, 1 });
}


TEST_CASE("Matrix - View", "[Matrix]") {
	Matrix<float, 4, 4> m = LookAt({ -6, -5, -5 }, { -1, 0, 0 }, Vector<float, 3>{ 0, 0, 1 }, true, false, false);

	Vector<float, 3> p = { 0, -1, 0 };
	Vector<float, 3> pt = p * m;
	Vector<float, 3> pexp = { sqrt(2), 0, 8.66025403 };
	REQUIRE(ApproxVec(pexp) == pt);

	m = LookAt({ 0, 0, 0 }, { 0, 5, 0 }, Vector<float, 3>{ 0, 0, 1 }, true, false, false);
	p = { 1, 4, 1 };
	pt = p * m;
	pexp = { 1, 1, 4 };

	REQUIRE(ApproxVec(pexp) == pt);

	m = LookAt({ 0, 0, 0 }, { 5, 0, 0 }, Vector<float, 3>{ 0, 0, 1 }, true, false, false);
	p = { 1, 4, 1 };
	pt = p * m;
	pexp = { -4, 1, 1 };

	REQUIRE(ApproxVec(pexp) == pt);
}