// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "TestGenerators.hpp"

#include <Mathter/Common/Approx.hpp>
#include <Mathter/Matrix.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>

using namespace mathter;
using Catch::Approx;


using TypeListFloating = TestTypeList<TypesFloating, PackedAll, OrdersAll, LayoutsAll>;
using TypeListFloatingFollow = TestTypeList<TypesFloating, PackedAll, OrdersFollow, LayoutsAll>;



TEMPLATE_LIST_TEST_CASE("Matrix - Identity", "[Transform]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using M35 = typename TestType::template Matrix<3, 5>;

		M33 m = Identity();
		M33 mexp = {
			1, 0, 0,
			0, 1, 0,
			0, 0, 1
		};

		REQUIRE(m == mexp);

		M35 m5 = Identity();
		M35 mexp5 = {
			1, 0, 0, 0, 0,
			0, 1, 0, 0, 0,
			0, 0, 1, 0, 0
		};

		REQUIRE(m5 == mexp5);
	}
}

TEMPLATE_LIST_TEST_CASE("Matrix - Zero", "[Transform]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using M34 = typename TestType::template Matrix<3, 4>;

		M34 m = Zero();
		M34 mexp = {
			0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0
		};

		REQUIRE(m == mexp);
	}
}



TEMPLATE_LIST_TEST_CASE("Matrix - Rotation2D FOLLOW", "[Matrix]", TypeListFloatingFollow) {
	SECTION(TestType::Name()) {
		using M22 = typename TestType::template Matrix<2, 2>;
		using M32 = typename TestType::template Matrix<3, 2>;
		using M33 = typename TestType::template Matrix<3, 3>;
		using M22I = invert_order_t<M22>;
		using M23I = invert_order_t<M32>;
		using M33I = invert_order_t<M33>;

		M22 m22 = Rotation(1.f);
		M32 m32 = Rotation(1.f);
		M33 m33 = Rotation(1.f);

		M22 m22exp = {
			0.54030, 0.84147,
			-0.84147, 0.54030
		};
		M32 m32exp = {
			0.54030, 0.84147,
			-0.84147, 0.54030,
			0, 0
		};
		M33 m33exp = {
			0.54030, 0.84147, 0,
			-0.84147, 0.54030, 0,
			0, 0, 1
		};
		REQUIRE(ApproxVec(m22) == m22exp);
		REQUIRE(ApproxVec(m32) == m32exp);
		REQUIRE(ApproxVec(m33) == m33exp);


		M22I m22p = Rotation(1.f);
		M23I m32p = Rotation(1.f);
		M33I m33p = Rotation(1.f);

		M22I m22expp = {
			0.54030, -0.84147,
			0.84147, 0.54030
		};
		M23I m32expp = {
			0.54030, -0.84147, 0,
			0.84147, 0.54030, 0
		};
		M33I m33expp = {
			0.54030, -0.84147, 0,
			0.84147, 0.54030, 0,
			0, 0, 1
		};
		REQUIRE(ApproxVec(m22p) == m22expp);
		REQUIRE(ApproxVec(m32p) == m32expp);
		REQUIRE(ApproxVec(m33p) == m33expp);
	}
}



TEMPLATE_LIST_TEST_CASE("Matrix - RotationPrincipal", "[Matrix]", TypeListFloatingFollow) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using M43 = typename TestType::template Matrix<4, 3>;
		using M44 = typename TestType::template Matrix<4, 4>;
		using M33I = invert_order_t<M33>;
		using M34I = invert_order_t<M43>;
		using M44I = invert_order_t<M44>;

		M33 m33 = RotationX(1.f);
		M33 m33exp = {
			1.000000, 0.000000, 0.000000,
			0.000000, 0.540302, 0.841471,
			0.000000, -0.841471, 0.540302
		};
		REQUIRE(ApproxVec(m33) == m33exp);

		M43 m43 = RotationY(1.f);
		M43 m43exp = {
			0.540302, 0.000000, -0.841471,
			0.000000, 1.000000, 0.000000,
			0.841471, 0.000000, 0.540302,
			0, 0, 0
		};
		REQUIRE(ApproxVec(m43) == m43exp);

		M44 m44 = RotationZ(1.f);
		M44 m44exp = {
			0.540302, 0.841471, 0.000000, 0,
			-0.841471, 0.540302, 0.000000, 0,
			0.000000, 0.000000, 1.000000, 0,
			0, 0, 0, 1
		};
		REQUIRE(ApproxVec(m44) == m44exp);


		M33I m33p = RotationX(1.f);
		M33I m33expp = matrix_representation_cast<M33I>(m33exp);
		REQUIRE(ApproxVec(m33p) == m33expp);

		M34I m34p = RotationY(1.f);
		M34I m34expp = matrix_representation_cast<M34I>(m43exp);
		REQUIRE(ApproxVec(m34p) == m34expp);

		M44I m44p = RotationZ(1.f);
		M44I m44expp = matrix_representation_cast<M44I>(m44exp);
		REQUIRE(ApproxVec(m44p) == m44expp);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - RotationTriAxis", "[Matrix]", TypeListFloatingFollow) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using M43 = typename TestType::template Matrix<4, 3>;
		using M44 = typename TestType::template Matrix<4, 4>;
		using M33I = invert_order_t<M33>;
		using M34I = invert_order_t<M43>;
		using M44I = invert_order_t<M44>;

		M33 m33 = RotationAxis3<0, 1, 1>(1.f, 1.0f, -1.0f);
		M33 m33exp = {
			1.000000, 0.000000, 0.000000,
			0.000000, 0.540302, 0.841471,
			0.000000, -0.841471, 0.540302
		};
		REQUIRE(ApproxVec(m33) == m33exp);

		M43 m43 = RotationAxis3<0, 1, 2>(0.0f, 1.f, 0.0f);
		M43 m43exp = {
			0.540302, 0.000000, -0.841471,
			0.000000, 1.000000, 0.000000,
			0.841471, 0.000000, 0.540302,
			0, 0, 0
		};
		REQUIRE(ApproxVec(m43) == m43exp);

		M44 m44 = RotationAxis3<0, 0, 2>(-1.0f, 1.0f, 1.f);
		M44 m44exp = {
			0.540302, 0.841471, 0.000000, 0,
			-0.841471, 0.540302, 0.000000, 0,
			0.000000, 0.000000, 1.000000, 0,
			0, 0, 0, 1
		};
		REQUIRE(ApproxVec(m44) == m44exp);



		M33I m33p = RotationAxis3<0, 1, 1>(1.f, 1.0f, -1.0f);
		M33I m33expp = matrix_representation_cast<M33I>(m33exp);
		REQUIRE(ApproxVec(m33p) == m33expp);

		M34I m34p = RotationAxis3<0, 1, 2>(0.0f, 1.f, 0.0f);
		M34I m34expp = matrix_representation_cast<M34I>(m43exp);
		REQUIRE(ApproxVec(m34p) == m34expp);

		M44I m44p = RotationAxis3<0, 0, 2>(-1.0f, 1.0f, 1.f);
		M44I m44expp = matrix_representation_cast<M44I>(m44exp);
		REQUIRE(ApproxVec(m44p) == m44expp);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - RotationAxisAngle", "[Matrix]", TypeListFloatingFollow) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using M43 = typename TestType::template Matrix<4, 3>;
		using M44 = typename TestType::template Matrix<4, 4>;
		using M33I = invert_order_t<M33>;
		using M34I = invert_order_t<M43>;
		using M44I = invert_order_t<M44>;


		M33 m33 = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
		M33 m33exp = {
			0.573138, 0.740349, -0.351279,
			-0.609007, 0.671645, 0.421906,
			0.548292, -0.027879, 0.835822
		};
		REQUIRE(ApproxVec(m33) == m33exp);

		M43 m43 = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
		M43 m43exp = {
			0.573138, 0.740349, -0.351279,
			-0.609007, 0.671645, 0.421906,
			0.548292, -0.027879, 0.835822,
			0, 0, 0
		};
		REQUIRE(ApproxVec(m43) == m43exp);

		M44 m44 = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
		M44 m44exp = {
			0.573138, 0.740349, -0.351279, 0,
			-0.609007, 0.671645, 0.421906, 0,
			0.548292, -0.027879, 0.835822, 0,
			0, 0, 0, 1
		};
		REQUIRE(ApproxVec(m44) == m44exp);


		M33I m33p = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
		M33I m33expp = matrix_representation_cast<M33I>(m33exp);
		REQUIRE(ApproxVec(m33p) == m33expp);

		M34I m34p = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
		M34I m34expp = matrix_representation_cast<M34I>(m43exp);
		REQUIRE(ApproxVec(m34p) == m34expp);

		M44I m44p = RotationAxisAngle(Normalize(Vector<float, 3>(1, 2, 3)), 1.0f);
		M44I m44expp = matrix_representation_cast<M44I>(m44exp);
		REQUIRE(ApproxVec(m44p) == m44expp);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Scale", "[Matrix]", TypeListFloatingFollow) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using M55 = typename TestType::template Matrix<5, 5>;
		using V5 = typename TestType::template Vector<5>;
		using V3 = typename TestType::template Vector<3>;
		using M55I = invert_order_t<M55>;

		M55 m = Scale(1, 2, 3, 4, 5);
		V5 v(2, 6, 3, 7, 5);
		M33 m3 = Scale(V3{ 1, 2, 3 });

		auto vt1 = v * V5{ 1, 2, 3, 4, 5 };
		auto vt2 = v * m;

		REQUIRE(vt1 == vt2);

		M55 mp = Scale(1, 2, 3, 4, 5);
		REQUIRE(ApproxVec(mp) == m);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Translation", "[Matrix]", TypeListFloatingFollow) {
	SECTION(TestType::Name()) {
		using M32 = typename TestType::template Matrix<3, 2>;
		using M33 = typename TestType::template Matrix<3, 3>;
		using M65 = typename TestType::template Matrix<6, 5>;
		using V5 = typename TestType::template Vector<5>;
		using M23I = invert_order_t<M32>;
		using M33I = invert_order_t<M33>;

		// FOLLOW
		M33 m2d_33a = Translation(1, 2);
		M33 m2d_33b = Translation(Vector<float, 2, false>(1, 2));
		M33 m2d_33exp = {
			1, 0, 0,
			0, 1, 0,
			1, 2, 1
		};
		REQUIRE(ApproxVec(m2d_33a) == m2d_33exp);
		REQUIRE(ApproxVec(m2d_33b) == m2d_33exp);


		M32 m2d_32 = Translation(1, 2);
		M32 m2d_32exp = {
			1, 0,
			0, 1,
			1, 2
		};
		REQUIRE(ApproxVec(m2d_32) == m2d_32exp);

		// PRECEDE
		M33I m2d_33p = Translation(1, 2);
		M33I m2d_33expp = {
			1, 0, 1,
			0, 1, 2,
			0, 0, 1
		};
		REQUIRE(ApproxVec(m2d_33a) == m2d_33exp);


		M23I m2d_23p = Translation(1, 2);
		M23I m2d_23expp = {
			1, 0, 1,
			0, 1, 2
		};
		REQUIRE(ApproxVec(m2d_23p) == m2d_23expp);


		// LARGE
		V5 t = { 1, 2, 3, 4, 5 };
		M65 m = Translation(t);
		V5 v(1, 2, 3, 4, 5);
		v *= m;
		V5 vexp(2, 4, 6, 8, 10);
		REQUIRE(v == vexp);
	}
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


TEMPLATE_LIST_TEST_CASE("Matrix - Perspective", "[Matrix]", TypeListFloatingFollow) {
	SECTION(TestType::Name()) {
		using M44 = typename TestType::template Matrix<4, 4>;
		using M44I = invert_order_t<M44>;
		using Type = scalar_type_t<M44>;
		constexpr auto Packed = is_packed_v<M44>;

		constexpr int numCases = 5;

		float fov = Deg2Rad(Type(45));
		float ar = 2.0f;

		float nearSet[numCases] = { 0.5f, -0.5f, 0.5f, -0.5f, 0.5f };
		float farSet[numCases] = { 10.0f, -10.0f, 10.0f, -10.0f, 10.f };
		float pnearSet[numCases] = { 0.0f, 0.0f, 1.0f, 1.0f, 2.0f };
		float pfarSet[numCases] = { 1.0f, 1.0f, -1.0f, -1.0f, 3.0f };
		const char* sections[numCases] = {
			"View Z+, NDC 0:1",
			"View Z-, NDC 0:1",
			"View Z+, NDC 1:-1",
			"View Z-, NDC 1:-1",
			"View Z+, NDC 3:2",
		};

		for (int i = 0; i < numCases; ++i) {
			SECTION(sections[i]) {
				float near = nearSet[i];
				float far = farSet[i];
				float pnear = pnearSet[i];
				float pfar = pfarSet[i];

				Frustum<Type, Packed> frustum(near, far, fov, ar);

				M44 m = Perspective(fov, ar, near, far, pnear, pfar);

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

				M44I mp = Perspective(fov, ar, near, far, pnear, pfar);
				REQUIRE(ApproxVec(mp) == matrix_representation_cast<decltype(mp)>(m));
			}
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Orthographic", "[Matrix]", TypeListFloatingFollow) {
	SECTION(TestType::Name()) {
		using M44 = typename TestType::template Matrix<4, 4>;
		using Vec = typename TestType::template Vector<3>;
		using Type = scalar_type_t<M44>;

		Vec worldFrustum[2] = {
			{ -0.25f, -0.44444444f, 0.5f },
			{ 5.0f, 8.8888888f, 10.f }
		};
		Vec ndcFrustum[2];

		// Z forward
		M44 m = Orthographic(Vec(worldFrustum[0]), Vec(worldFrustum[1]), Type(0), Type(1));
		ndcFrustum[0] = worldFrustum[0] * m;
		ndcFrustum[1] = worldFrustum[1] * m;

		REQUIRE(ApproxVec(ndcFrustum[0]) == Vec{ -1, -1, 0 });
		REQUIRE(ApproxVec(ndcFrustum[1]) == Vec{ 1, 1, 1 });
	}
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



TEMPLATE_LIST_TEST_CASE("Matrix - View", "[Matrix]", TypeListFloatingFollow) {
	SECTION(TestType::Name()) {
		using M44 = typename TestType::template Matrix<4, 4>;
		using M44I = invert_order_t<M44>;
		using Vec = typename TestType::template Vector<3>;
		using Type = scalar_type_t<M44>;
		constexpr auto Packed = is_packed_v<M44>;

		Basis<Type, Packed> basis;

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

		Vector<float, 3, false> eye = basis.center;
		Vector<float, 3, false> target = basis.center + 2 * basis.basis1;
		Vector<float, 3, false> up = Normalize(basis.basis3 + Type(0.1) * basis.basis1);

		M44 m = LookAt(eye, target, up, true, false, false);
		M44 mfff = LookAt(eye, target, up, false, false, false);
		M44 mftf = LookAt(eye, target, up, false, true, false);
		M44 mftt = LookAt(eye, target, up, false, true, true);

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

		M44I mp = LookAt(eye, target, up, true, false, false);
		REQUIRE(ApproxVec(mp) == matrix_representation_cast<decltype(mp)>(m));
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Perspective 2D", "[Matrix]", TypeListFloatingFollow) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using Vec = typename TestType::template Vector<2>;

		std::array<Vec, 4> frustum = {
			Vec(-2, 2), Vec(2, 2),
			Vec(-1, 1), Vec(1, 1)
		};

		M33 m = Perspective(Deg2Rad(90.f), 1.0f, 2.0f, 0.0f, 1.0f);

		std::array<Vec, 4> ndcexp = {
			Vec(-1, 1), Vec(1, 1),
			Vec(-1, 0), Vec(1, 0)
		};

		for (int i = 0; i < 4; ++i) {
			REQUIRE(frustum[i] * m == ApproxVec(ndcexp[i]));
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - View 2D", "[Matrix]", TypeListFloatingFollow) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using Vec = typename TestType::template Vector<2>;

		Vector<float, 2, false> eye = { 3, 4 };
		Vector<float, 2, false> target = { 6, 5 };
		Vector<float, 2, false> test = { 4, 4 };

		M33 m = LookAt(eye, target, true, false);

		REQUIRE(Vec(eye) * m == ApproxVec(Vec(0, 0)));
		REQUIRE(Normalize(Vec(target) * m) == ApproxVec(Vec(0, 1)));
		REQUIRE((Vec(test) * m).x > 0);
		REQUIRE((Vec(test) * m).y > 0);
	}
}