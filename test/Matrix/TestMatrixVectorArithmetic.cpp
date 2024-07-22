// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../TestGenerators.hpp"

#include <Mathter/Common/Approx.hpp>
#include <Mathter/Matrix.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>



using namespace mathter;
using Catch::Approx;


using TypeListFollow = TestTypeList<TypesAll, PackedAll, OrdersFollow, LayoutsAll>;


TEMPLATE_LIST_TEST_CASE("Matrix - Matrix-vector square multiplication", "[Matrix]", TypeListFollow) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using M33I = invert_order_t<M33>;
		using V3 = typename TestType::template Vector<3>;

		M33 m = {
			1, 2, 3,
			4, 5, 6,
			7, 8, 9
		};
		M33I mT = matrix_representation_cast<M33I>(m);
		V3 v = {
			5,
			7,
			11
		};

		const auto p1 = mT * v;
		const auto p2 = v * m;

		V3 expected = {
			110,
			133,
			156
		};

		REQUIRE(p1 == ApproxVec(expected));
		REQUIRE(p2 == ApproxVec(expected));
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Matrix-vector non-square multiplication", "[Matrix]", TypeListFollow) {
	SECTION(TestType::Name()) {
		using M43 = typename TestType::template Matrix<4, 3>;
		using M34I = invert_order_t<M43>;
		using V3 = typename TestType::template Vector<3>;
		using V4 = typename TestType::template Vector<4>;

		M43 m = {
			1, 2, 3,
			4, 5, 6,
			7, 8, 9,
			6, 7, 8
		};
		M34I mT = matrix_representation_cast<M34I>(m);
		V4 v = {
			5,
			7,
			11,
			1
		};

		auto p1 = mT * v;
		auto p2 = v * m;

		V3 expected = {
			116,
			140,
			164
		};

		REQUIRE(p1 == ApproxVec(expected));
		REQUIRE(p2 == ApproxVec(expected));
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Matrix-vector implicit affine multiplication", "[Matrix]", TypeListFollow) {
	SECTION(TestType::Name()) {
		using M43 = typename TestType::template Matrix<4, 3>;
		using M34I = invert_order_t<M43>;
		using V3 = typename TestType::template Vector<3>;
		using V4 = typename TestType::template Vector<4>;

		M43 m = {
			1, 2, 3,
			4, 5, 6,
			7, 8, 9,
			6, 7, 8
		};
		M34I mT = matrix_representation_cast<M34I>(m);
		V3 v = {
			5,
			7,
			11
		};

		auto p1 = mT * v;
		auto p2 = v * m;

		V3 expected = {
			116,
			140,
			164
		};

		REQUIRE(p1 == ApproxVec(expected));
		REQUIRE(p2 == ApproxVec(expected));
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Matrix-vector implicit homogeneous multiplication", "[Matrix]", TypeListFollow) {
	SECTION(TestType::Name()) {
		using M44 = typename TestType::template Matrix<4, 4>;
		using M44I = invert_order_t<M44>;
		using V3 = typename TestType::template Vector<3>;
		using V4 = typename TestType::template Vector<4>;

		M44 m = {
			1, 2, 3, 3,
			4, 5, 6, 7,
			7, 8, 9, 2,
			6, 7, 8, 3
		};
		M44I mT = matrix_representation_cast<M44I>(m);
		V3 v = {
			5,
			7,
			11
		};

		auto p1 = mT * v;
		auto p2 = v * m;

		V3 expected = {
			116.f / 89.f,
			140.f / 89.f,
			164.f / 89.f
		};

		REQUIRE(p1 == ApproxVec(expected));
		REQUIRE(p2 == ApproxVec(expected));
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Matrix-vector compound multiplication", "[Matrix]", TypeListFollow) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using M43 = typename TestType::template Matrix<4, 3>;
		using M44 = typename TestType::template Matrix<4, 4>;
		using V3 = typename TestType::template Vector<3>;

		M33 m33 = {
			1, 2, 3,
			4, 5, 6,
			7, 8, 9
		};
		M43 m43 = {
			1, 2, 3,
			4, 5, 6,
			7, 8, 9,
			6, 7, 8
		};
		M44 m44 = {
			1, 2, 3, 3,
			4, 5, 6, 7,
			7, 8, 9, 2,
			6, 7, 8, 3
		};

		const V3 v = { 1, 2, 3 };

		auto v1 = v;
		auto v2 = v;
		auto v3 = v;

		v1 *= m33;
		v2 *= m43;
		v3 *= m44;

		REQUIRE(ApproxVec(v * m33) == v1);
		REQUIRE(ApproxVec(v * m43) == v2);
		REQUIRE(ApproxVec(v * m44) == v3);
	}
}