//L=============================================================================
//L This software is distributed under the MIT license.
//L Copyright 2021 Péter Kardos
//L=============================================================================

#pragma warning(disable : 4244)

#include <Mathter/Common/Approx.hpp>
#include <Mathter/Matrix.hpp>
#include "TestGenerators.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <complex>


using namespace mathter;
using Catch::Approx;


TEST_CASE_VARIANT("Matrix - Matrix-vector square multiplication", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	MatrixT<3, 3> m = {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9
	};
	MatrixTOO<3, 3> mT = matrix_representation_cast<MatrixTOO<3, 3>>(m);
	Vector<Type, 3, Packed> v = {
		5,
		7,
		11
	};

	auto p1 = mT * v;
	auto p2 = v * m;

	Vector<Type, 3, Packed> expected = {
		52,
		121,
		190
	};

	REQUIRE(p1 == ApproxVec(expected));
	REQUIRE(p2 == ApproxVec(expected));
}


TEST_CASE_VARIANT("Matrix - Matrix-vector non-square multiplication", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	MatrixT<4, 3> m = {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		6, 7, 8
	};
	MatrixTOO<3, 4> mT = matrix_representation_cast<MatrixTOO<3, 4>>(m);
	Vector<Type, 4, Packed> v = {
		5,
		7,
		11,
		1
	};

	auto p1 = mT * v;
	auto p2 = v * m;

	Vector<Type, 3, Packed> expected = {
		116,
		140,
		164
	};

	REQUIRE(p1 == ApproxVec(expected));
	REQUIRE(p2 == ApproxVec(expected));
}


TEST_CASE_VARIANT("Matrix - Matrix-vector implicit affine multiplication", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	MatrixT<4, 3> m = {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		6, 7, 8
	};
	MatrixTOO<3, 4> mT = matrix_representation_cast<MatrixTOO<3, 4>>(m);
	Vector<Type, 3, Packed> v = {
		5,
		7,
		11
	};

	auto p1 = mT * v;
	auto p2 = v * m;

	Vector<Type, 3, Packed> expected = {
		116,
		140,
		164
	};

	REQUIRE(p1 == ApproxVec(expected));
	REQUIRE(p2 == ApproxVec(expected));
}


TEST_CASE_VARIANT("Matrix - Matrix-vector implicit homogeneous multiplication", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	MatrixT<4, 4> m = {
		1, 2, 3, 3,
		4, 5, 6, 7, 
		7, 8, 9, 2,
		6, 7, 8, 3
	};
	MatrixTOO<4, 4> mT = matrix_representation_cast<MatrixTOO<4, 4>>(m);
	Vector<Type, 3, Packed> v = {
		5,
		7,
		11
	};

	auto p1 = mT * v;
	auto p2 = v * m;

	Vector<Type, 3, Packed> expected = {
		116.f/89.f,
		140.f/89.f,
		164.f/89.f
	};

	REQUIRE(p1 == ApproxVec(expected));
	REQUIRE(p2 == ApproxVec(expected));
}


TEST_CASE_VARIANT("Matrix - Matrix-vector compound multiplication", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	MatrixT<3, 3> m33 = {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9
	};
	MatrixT<4, 3> m43 = {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		6, 7, 8
	};
	MatrixT<4, 4> m44 = {
		1, 2, 3, 3,
		4, 5, 6, 7,
		7, 8, 9, 2,
		6, 7, 8, 3
	};

	const Vector<Type, 3, Packed> v = { 1, 2, 3 };

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