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


using TypeListAll = TestTypeList<TypesFloating, PackedAll, OrdersAll, LayoutsAll>;
using TypeListFloating = TestTypeList<TypesFloating, PackedAll, OrdersAll, LayoutsAll>;


TEST_CASE("Matrix - Min", "[Matrix]") {
	using M22 = Matrix<float, 2, 2>;

	SECTION("Underflow") {
		M22 m = { -1, 3, 2, 6 };
		REQUIRE(Min(m) == -1);
	}
}


TEST_CASE("Matrix - Max", "[Matrix]") {
	using M22 = Matrix<float, 2, 2>;

	SECTION("Underflow") {
		M22 m = { -1, 3, 2, 6 };
		REQUIRE(Max(m) == 6);
	}
}


TEST_CASE("Matrix - Abs", "[Matrix]") {
	using M22 = Matrix<float, 2, 2>;

	SECTION("Underflow") {
		M22 m = { -1, 3, 2, 6 };
		M22 e = { 1, 3, 2, 6 };
		REQUIRE(Abs(m) == e);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Trace", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M33 = typename TestType::template Matrix<3, 3>;
		using M55 = typename TestType::template Matrix<5, 5>;

		M33 m = {
			1, 3, 2,
			4, 5, 6,
			7, 8, 9
		};
		auto trace = Trace(m);

		REQUIRE(Approx(trace) == 15.f);

		M55 m5 = {
			5, 7, 3, 6, 4,
			4, 7, 4, 6, 3,
			6, 2, 8, 9, 7,
			1, 2, 7, 4, 8,
			5, 9, 7, 1, 5
		};
		trace = Trace(m5);
		REQUIRE(Approx(trace) == 29);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Transpose", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M42 = typename TestType::template Matrix<4, 2>;
		using M24 = typename TestType::template Matrix<2, 4>;

		M42 m = {
			1, 2,
			3, 4,
			5, 6,
			7, 8
		};
		M24 mT = Transpose(m);
		M24 mexp = {
			1, 3, 5, 7,
			2, 4, 6, 8
		};

		REQUIRE(mT == mexp);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Determinant small matrix", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M22 = typename TestType::template Matrix<2, 2>;
		using M33 = typename TestType::template Matrix<3, 3>;
		using M44 = typename TestType::template Matrix<4, 4>;

		M22 m2 = {
			1, 3,
			4, 5
		};
		REQUIRE(ApproxVec(Determinant(m2)) == -7);

		M44 m4 = {
			1, 3, 2, 1,
			4, 5, 6, 2,
			7, 8, 9, 3,
			1, 2, 3, 4
		};
		REQUIRE(ApproxVec(Determinant(m4)) == 27);

		M33 m3 = {
			1, 3, 2,
			4, 5, 6,
			7, 8, 9
		};
		REQUIRE(Approx(Determinant(m3)) == 9);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Determinant", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M55 = typename TestType::template Matrix<5, 5>;

		M55 m5 = {
			5, 7, 3, 6, 4,
			4, 7, 4, 6, 3,
			6, 2, 8, 9, 7,
			1, 2, 7, 4, 8,
			5, 9, 7, 1, 5
		};
		REQUIRE(Approx(Determinant(m5)) == 4134);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Inverse Small Matrix", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M22 = typename TestType::template Matrix<2, 2>;
		using M33 = typename TestType::template Matrix<3, 3>;
		using M44 = typename TestType::template Matrix<4, 4>;

		M22 m2 = {
			1, 3,
			4, 5
		};
		M22 mI2 = Inverse(m2);
		M22 mexp2 = {
			-0.714286, 0.428571,
			0.571429, -0.142857
		};

		REQUIRE(ApproxVec(mI2) == mexp2);

		M33 m3 = {
			1, 3, 2,
			4, 5, 6,
			7, 8, 9
		};
		M33 mI3 = Inverse(m3);
		M33 mexp3 = {
			-0.333333, -1.222222, 0.888889,
			0.666667, -0.555556, 0.222222,
			-0.333333, 1.444444, -0.777778
		};

		REQUIRE(ApproxVec(mI3) == mexp3);

		M44 m4 = {
			1, 3, 2, 1,
			4, 5, 6, 2,
			7, 8, 9, 3,
			1, 2, 3, 4
		};
		M44 mI4 = Inverse(m4);
		M44 mexp4 = {
			-0.333333, -1.296296, 0.925926, 0.037037,
			0.666667, -0.407407, 0.148148, -0.074074,
			-0.333333, 1.592593, -0.851852, -0.074074,
			0, -0.666667, 0.333333, 0.333333
		};

		REQUIRE(ApproxVec(mI4) == mexp4);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Inverse", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using M55 = typename TestType::template Matrix<5, 5>;

		M55 n = {
			1, 56, 8, 4, 3,
			4, 2, 7, 8, 4,
			1, 5, 7, 4, 3,
			9, 5, 3, 8, 4,
			7, 2, 83, 46, 4
		};
		M55 nI = Inverse(n);
		M55 iden = n * nI;
		M55 idenexp;
		idenexp = Identity();

		REQUIRE(ApproxVec(idenexp) == iden);
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - Norm", "[Matrix]", TypeListAll) {
	SECTION(TestType::Name()) {
		using V8 = typename TestType::template Vector<8>;
		using M24 = typename TestType::template Matrix<2, 4>;

		V8 v = { 1, 2, 3, 4, 5, 6, 7, 8 };
		M24 m = { 1, 2, 3, 4, 5, 6, 7, 8 };
		REQUIRE(Approx(Length(v)) == Norm(m));
	}
}


TEST_CASE("Matrix - Norm precise", "[Matrix]") {
	using M22 = Matrix<float, 2, 2>;

	SECTION("Underflow") {
		M22 m = { 1e-30f, 1e-30f, 1e-30f, 1e-30f };
		REQUIRE(Approx(2e-30f) == NormPrecise(m));
	}
	SECTION("Denormal") {
		M22 m = { 1e-39f, 1e-39f, 1e-39f, 1e-39f };
		REQUIRE(Approx(2e-39f) == NormPrecise(m));
	}
}