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
	MatrixT<3, 3> m = Identity();
	MatrixT<3, 3> mexp = {
		1, 0, 0,
		0, 1, 0,
		0, 0, 1
	};

	REQUIRE(m == mexp);

	MatrixT<3, 5> m5 = Identity();
	MatrixT<3, 5> mexp5 = {
		1, 0, 0, 0, 0,
		0, 1, 0, 0, 0,
		0, 0, 1, 0, 0
	};

	REQUIRE(m5 == mexp5);
}

TEST_CASE_VARIANT("Matrix - Zero", "[Transform]", TypesAll, OrdersAll, LayoutsAll, PackedAll) {
	MatrixT<3, 4> m = Zero();
	MatrixT<3, 4> mexp = {
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0
	};

	REQUIRE(m == mexp);
}



TEST_CASE_VARIANT("Matrix - Rotation2D FOLLOW", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	MatrixT<2, 2> m22 = Rotation(1.f);
	MatrixT<3, 2> m32 = Rotation(1.f);
	MatrixT<3, 3> m33 = Rotation(1.f);

	MatrixT<2, 2> m22exp = {
		0.54030, 0.84147,
		-0.84147, 0.54030
	};
	MatrixT<3, 2> m32exp = {
		0.54030, 0.84147,
		-0.84147, 0.54030,
		0, 0
	};
	MatrixT<3, 3> m33exp = {
		0.54030, 0.84147, 0,
		-0.84147, 0.54030, 0,
		0, 0, 1
	};
	REQUIRE(ApproxVec(m22) == m22exp);
	REQUIRE(ApproxVec(m32) == m32exp);
	REQUIRE(ApproxVec(m33) == m33exp);


	MatrixTOO<2, 2> m22p = Rotation(1.f);
	MatrixTOO<2, 3> m32p = Rotation(1.f);
	MatrixTOO<3, 3> m33p = Rotation(1.f);

	MatrixTOO<2, 2> m22expp = {
		0.54030, -0.84147,
		0.84147, 0.54030
	};
	MatrixTOO<2, 3> m32expp = {
		0.54030, -0.84147, 0,
		0.84147, 0.54030, 0
	};
	MatrixTOO<3, 3> m33expp = {
		0.54030, -0.84147, 0,
		0.84147, 0.54030, 0,
		0, 0, 1
	};
	REQUIRE(ApproxVec(m22p) == m22expp);
	REQUIRE(ApproxVec(m32p) == m32expp);
	REQUIRE(ApproxVec(m33p) == m33expp);
}



TEST_CASE_VARIANT("Matrix - RotationPrincipal", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	MatrixT<3, 3> m33 = RotationX(1.f);
	MatrixT<3, 3> m33exp = {
		1.000000, 0.000000, 0.000000,
		0.000000, 0.540302, 0.841471,
		0.000000, -0.841471, 0.540302
	};
	REQUIRE(ApproxVec(m33) == m33exp);

	MatrixT<4, 3> m43 = RotationY(1.f);
	MatrixT<4, 3> m43exp = {
		0.540302, 0.000000, -0.841471,
		0.000000, 1.000000, 0.000000,
		0.841471, 0.000000, 0.540302,
		0, 0, 0
	};
	REQUIRE(ApproxVec(m43) == m43exp);

	MatrixT<4, 4> m44 = RotationZ(1.f);
	MatrixT<4, 4> m44exp = {
		0.540302, 0.841471, 0.000000, 0,
		-0.841471, 0.540302, 0.000000, 0,
		0.000000, 0.000000, 1.000000, 0,
		0, 0, 0, 1
	};
	REQUIRE(ApproxVec(m44) == m44exp);


	MatrixTOO<3, 3> m33p = RotationX(1.f);
	MatrixTOO<3, 3> m33expp = matrix_representation_cast<MatrixTOO<3, 3>>(m33exp);
	REQUIRE(ApproxVec(m33p) == m33expp);

	MatrixTOO<3, 4> m34p = RotationY(1.f);
	MatrixTOO<3, 4> m34expp = matrix_representation_cast<MatrixTOO<3, 4>>(m43exp);
	REQUIRE(ApproxVec(m34p) == m34expp);

	MatrixTOO<4, 4> m44p = RotationZ(1.f);
	MatrixTOO<4, 4> m44expp = matrix_representation_cast<MatrixTOO<4, 4>>(m44exp);
	REQUIRE(ApproxVec(m44p) == m44expp);
}


TEST_CASE_VARIANT("Matrix - RotationTriAxis", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	MatrixT<3, 3> m33 = RotationAxis3<0, 1, 1>(1.f, 1.0f, -1.0f);
	MatrixT<3, 3> m33exp = {
		1.000000, 0.000000, 0.000000,
		0.000000, 0.540302, 0.841471,
		0.000000, -0.841471, 0.540302
	};
	REQUIRE(ApproxVec(m33) == m33exp);

	MatrixT<4, 3> m43 = RotationAxis3<0, 1, 2>(0.0f, 1.f, 0.0f);
	MatrixT<4, 3> m43exp = {
		0.540302, 0.000000, -0.841471,
		0.000000, 1.000000, 0.000000,
		0.841471, 0.000000, 0.540302,
		0, 0, 0
	};
	REQUIRE(ApproxVec(m43) == m43exp);

	MatrixT<4, 4> m44 = RotationAxis3<0, 0, 2>(-1.0f, 1.0f, 1.f);
	MatrixT<4, 4> m44exp = {
		0.540302, 0.841471, 0.000000, 0,
		-0.841471, 0.540302, 0.000000, 0,
		0.000000, 0.000000, 1.000000, 0,
		0, 0, 0, 1
	};
	REQUIRE(ApproxVec(m44) == m44exp);



	MatrixTOO<3, 3> m33p = RotationAxis3<0, 1, 1>(1.f, 1.0f, -1.0f);
	MatrixTOO<3, 3> m33expp = matrix_representation_cast<MatrixTOO<3, 3>>(m33exp);
	REQUIRE(ApproxVec(m33p) == m33expp);

	MatrixTOO<3, 4> m34p = RotationAxis3<0, 1, 2>(0.0f, 1.f, 0.0f);
	MatrixTOO<3, 4> m34expp = matrix_representation_cast<MatrixTOO<3, 4>>(m43exp);
	REQUIRE(ApproxVec(m34p) == m34expp);

	MatrixTOO<4, 4> m44p = RotationAxis3<0, 0, 2>(-1.0f, 1.0f, 1.f);
	MatrixTOO<4, 4> m44expp = matrix_representation_cast<MatrixTOO<4, 4>>(m44exp);
	REQUIRE(ApproxVec(m44p) == m44expp);
}


TEST_CASE_VARIANT("Matrix - RotationAxisAngle", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	MatrixT<3, 3> m33 = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
	MatrixT<3, 3> m33exp = {
		0.573138, 0.740349, -0.351279,
		-0.609007, 0.671645, 0.421906,
		0.548292, -0.027879, 0.835822
	};
	REQUIRE(ApproxVec(m33) == m33exp);

	MatrixT<4, 3> m43 = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
	MatrixT<4, 3> m43exp = {
		0.573138, 0.740349, -0.351279,
		-0.609007, 0.671645, 0.421906,
		0.548292, -0.027879, 0.835822,
		0, 0, 0
	};
	REQUIRE(ApproxVec(m43) == m43exp);

	MatrixT<4, 4> m44 = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
	MatrixT<4, 4> m44exp = {
		0.573138, 0.740349, -0.351279, 0,
		-0.609007, 0.671645, 0.421906, 0,
		0.548292, -0.027879, 0.835822, 0,
		0, 0, 0, 1
	};
	REQUIRE(ApproxVec(m44) == m44exp);


	MatrixTOO<3, 3> m33p = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
	MatrixTOO<3, 3> m33expp = matrix_representation_cast<MatrixTOO<3, 3>>(m33exp);
	REQUIRE(ApproxVec(m33p) == m33expp);

	MatrixTOO<3, 4> m34p = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
	MatrixTOO<3, 4> m34expp = matrix_representation_cast<MatrixTOO<3, 4>>(m43exp);
	REQUIRE(ApproxVec(m34p) == m34expp);

	MatrixTOO<4, 4> m44p = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
	MatrixTOO<4, 4> m44expp = matrix_representation_cast<MatrixTOO<4, 4>>(m44exp);
	REQUIRE(ApproxVec(m44p) == m44expp);
}


TEST_CASE_VARIANT("Matrix - Scale", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	MatrixT<5, 5> m = Scale(1, 2, 3, 4, 5);
	Vector<Type, 5, Packed> v(2, 6, 3, 7, 5);
	MatrixT<3, 3> m3 = Scale(Vector<Type, 3>{ 1, 2, 3 });

	auto vt1 = v * Vector<Type, 5, Packed>{ 1, 2, 3, 4, 5 };
	auto vt2 = v * m;

	REQUIRE(vt1 == vt2);

	MatrixTOO<5, 5> mp = Scale(1, 2, 3, 4, 5);
	REQUIRE(ApproxVec(mp) == m);
}


TEST_CASE_VARIANT("Matrix - Translation", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	// FOLLOW
	MatrixT<3, 3> m2d_33a = Translation(1, 2);
	MatrixT<3, 3> m2d_33b = Translation(Vector<Type, 2, Packed>(1, 2));
	MatrixT<3, 3> m2d_33exp = {
		1, 0, 0,
		0, 1, 0,
		1, 2, 1
	};
	REQUIRE(ApproxVec(m2d_33a) == m2d_33exp);
	REQUIRE(ApproxVec(m2d_33b) == m2d_33exp);


	MatrixT<3, 2> m2d_32 = Translation(1, 2);
	MatrixT<3, 2> m2d_32exp = {
		1, 0,
		0, 1,
		1, 2
	};
	REQUIRE(ApproxVec(m2d_32) == m2d_32exp);

	// PRECEDE
	MatrixTOO<3, 3> m2d_33p = Translation(1, 2);
	MatrixTOO<3, 3> m2d_33expp = {
		1, 0, 1,
		0, 1, 2,
		0, 0, 1
	};
	REQUIRE(ApproxVec(m2d_33a) == m2d_33exp);


	MatrixTOO<2, 3> m2d_23p = Translation(1, 2);
	MatrixTOO<2, 3> m2d_23expp = {
		1, 0, 1,
		0, 1, 2
	};
	REQUIRE(ApproxVec(m2d_23p) == m2d_23expp);


	// LARGE
	Vector<Type, 5, Packed> t = { 1, 2, 3, 4, 5 };
	MatrixT<6, 5> m = Translation(t);
	Vector<Type, 5, Packed> v(1, 2, 3, 4, 5);
	v *= m;
	Vector<Type, 5, Packed> vexp(2, 4, 6, 8, 10);
	REQUIRE(v == vexp);
}


template <class T, bool Packed>
struct Frustum {
	Frustum(T nearPlane, T farPlane, T angle, T aspect) {
		T dxdz = std::tan(angle / 2);
		T dydz = dxdz / aspect;
		Vector<T, 3> dir = { dxdz, dydz, 1 };
		for (int i = 0; i < 4; ++i) {
			T xside = XSide(i);
			T yside = YSide(i);
			Vector<T, 3> side = { xside, yside, 1 };
			near[i] = side * dir * Vector<T, 3>(std::abs(nearPlane), std::abs(nearPlane), nearPlane);
			far[i] = side * dir * Vector<T, 3>(std::abs(farPlane), std::abs(farPlane), farPlane);
		}

		for (int i = 0; i < 16; ++i) {
			Vector<T, 3> scale = { 0.3f, 0.3f, 1 };
			float t = float(i + 1) / float(17);
			float distance = t * farPlane + (1 - t) * nearPlane;
			middle[i] = dir * scale * Vector<T, 3>(std::abs(distance), std::abs(distance), distance);
		}
	}

	static T XSide(int i) {
		return (i == 0 || i == 3) ? -1 : 1;
	}
	static T YSide(int i) {
		return (i == 0 || i == 1) ? -1 : 1;
	}

	Vector<T, 3, Packed> near[4];
	Vector<T, 3, Packed> far[4];
	Vector<T, 3, Packed> middle[16];
};


TEST_CASE_VARIANT("Matrix - Perspective", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	constexpr int numCases = 5;

	Type fov = Deg2Rad(Type(45));
	Type ar = 2.0f;

	Type nearSet[numCases] = { 0.5f, -0.5f, 0.5f, -0.5f, 0.5f };
	Type farSet[numCases] = { 10.0f, -10.0f, 10.0f, -10.0f, 10.f };
	Type pnearSet[numCases] = { 0.0f, 0.0f, 1.0f, 1.0f, 2.0f };
	Type pfarSet[numCases] = { 1.0f, 1.0f, -1.0f, -1.0f, 3.0f };
	const char* sections[numCases] = {
		"View Z+, NDC 0:1",
		"View Z-, NDC 0:1",
		"View Z+, NDC 1:-1",
		"View Z-, NDC 1:-1",
		"View Z+, NDC 3:2",
	};

	for (int i = 0; i < numCases; ++i) {
		SECTION(sections[i]) {
			Type near = nearSet[i];
			Type far = farSet[i];
			Type pnear = pnearSet[i];
			Type pfar = pfarSet[i];

			Frustum<Type, Packed> frustum(near, far, fov, ar);

			MatrixT<4, 4> m = Perspective(fov, ar, near, far, pnear, pfar);

			for (int i = 0; i < 4; ++i) {
				auto vt = frustum.near[i] * m;
				REQUIRE(vt.z == Approx(pnear));
				REQUIRE(vt.x == Approx(frustum.XSide(i)));
				REQUIRE(vt.y == Approx(frustum.YSide(i)));
			}
			for (int i = 0; i < 4; ++i) {
				auto vt = frustum.far[i] * m;
				REQUIRE(vt.z == Approx(pfar));
				REQUIRE(vt.x == Approx(frustum.XSide(i)));
				REQUIRE(vt.y == Approx(frustum.YSide(i)));
			}
			for (const auto& v : frustum.middle) {
				auto vt = v * m;
				REQUIRE(vt.z <= Approx(std::max(pnear, pfar)));
				REQUIRE(vt.z >= Approx(std::min(pfar, pnear)));
			}

			MatrixTOO<4, 4> mp = Perspective(fov, ar, near, far, pnear, pfar);
			REQUIRE(ApproxVec(mp) == matrix_representation_cast<decltype(mp)>(m));
		}
	}
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


template <class T, bool Packed>
struct Basis {
	using Vec = Vector<T, 3, Packed>;
	Basis() {
		basis1 = Normalize(Vec{ -1, 3, 0 });
		basis2 = Normalize(Vec{ 3, 1, 0 });
		basis3 = Normalize(Vec{ 0, 0, 1 });
		center = { 6, 5, 8 };
		assert(Approx(0) == Dot(basis1, basis2));
		assert(Approx(0) == Dot(basis1, basis3));
		assert(Approx(0) == Dot(basis3, basis2));
	}

	Vec Express(Vec v) {
		return v(0) * basis1 + v(1) * basis2 + v(2) * basis3 + center;
	}

	Vec basis1;
	Vec basis2;
	Vec basis3;
	Vec center;
};



TEST_CASE_VARIANT("Matrix - View", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	Basis<Type, Packed> basis;

	using Vec = Vector<Type, 3, Packed>;

	std::array<Vec, 6> viewVecs = {
		Vec{ 1, 2, 3 },
		Vec{ 5, -5, 3 },
		Vec{ 1, 7, -1 },
		Vec{ 9, 3, -2 },
		Vec{ 9, 3, 4 },
		Vec{ -4, -3, 4 },
	};
	std::array<Vec, 6> worldVecs;
	for (int i = 0; i < 6; ++i) {
		worldVecs[i] = basis.Express(viewVecs[i]);
	}

	Vec eye = basis.center;
	Vec target = basis.center + 2 * basis.basis1;
	Vec up = Normalize(basis.basis3 + Type(0.1) * basis.basis1);

	MatrixT<4, 4> m = LookAt(eye, target, up, true, false, false);
	MatrixT<4, 4> mfff = LookAt(eye, target, up, false, false, false);
	MatrixT<4, 4> mftf = LookAt(eye, target, up, false, true, false);
	MatrixT<4, 4> mftt = LookAt(eye, target, up, false, true, true);

	REQUIRE((basis.center + basis.basis1) * m == ApproxVec(Vec(0, 0, 1)));
	REQUIRE((basis.center + basis.basis2) * m == ApproxVec(Vec(1, 0, 0)));
	REQUIRE((basis.center + basis.basis3) * m == ApproxVec(Vec(0, 1, 0)));

	REQUIRE((basis.center + basis.basis1) * mfff == ApproxVec(Vec(0, 0, -1)));
	REQUIRE((basis.center + basis.basis2) * mfff == ApproxVec(Vec(1, 0, 0)));
	REQUIRE((basis.center + basis.basis3) * mfff == ApproxVec(Vec(0, 1, 0)));

	REQUIRE((basis.center + basis.basis1) * mftf == ApproxVec(Vec(0, 0, -1)));
	REQUIRE((basis.center + basis.basis2) * mftf == ApproxVec(Vec(-1, 0, 0)));
	REQUIRE((basis.center + basis.basis3) * mftf == ApproxVec(Vec(0, 1, 0)));

	REQUIRE((basis.center + basis.basis1) * mftt == ApproxVec(Vec(0, 0, -1)));
	REQUIRE((basis.center + basis.basis2) * mftt == ApproxVec(Vec(-1, 0, 0)));
	REQUIRE((basis.center + basis.basis3) * mftt == ApproxVec(Vec(0, -1, 0)));

	for (int i = 0; i < 6; ++i) {
		REQUIRE(worldVecs[i] * m == ApproxVec(Vec(viewVecs[i].yzx)));
	}

	MatrixTOO<4, 4> mp = LookAt(eye, target, up, true, false, false);
	REQUIRE(ApproxVec(mp) == matrix_representation_cast<decltype(mp)>(m));
}


TEST_CASE_VARIANT("Matrix - Perspective 2D", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	using Vec = Vector<Type, 2, Packed>;

	std::array<Vec, 4> frustum = {
		Vec(-2, 2), Vec(2, 2),
		Vec(-1, 1), Vec(1, 1)
	};

	MatrixT<3, 3> m = Perspective(Deg2Rad(90.f), 1.0f, 2.0f, 0.0f, 1.0f);

	std::array<Vec, 4> ndcexp = {
		Vec(-1, 1), Vec(1, 1),
		Vec(-1, 0), Vec(1, 0)
	};

	for (int i = 0; i < 4; ++i) {
		REQUIRE(frustum[i] * m == ApproxVec(ndcexp[i]));
	}
}


TEST_CASE_VARIANT("Matrix - View 2D", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	using Vec = Vector<Type, 2, Packed>;

	Vec eye = { 3, 4 };
	Vec target = { 6, 5 };
	Vec test = { 4, 4 };

	MatrixT<3, 3> m = LookAt(eye, target, true, false);

	REQUIRE(eye * m == ApproxVec(Vec(0, 0)));
	REQUIRE(Normalize(target * m) == ApproxVec(Vec(0, 1)));
	REQUIRE((test * m).x > 0);
	REQUIRE((test * m).y > 0);
}