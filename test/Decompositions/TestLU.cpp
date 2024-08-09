// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../ApplyTransform.hpp"
#include "../Approx.hpp"
#include "../Cases.hpp"
#include "../MatrixUtil.hpp"

#include <Mathter/Decompositions/DecomposeLU.hpp>
#include <Mathter/Matrix/Arithmetic.hpp>
#include <Mathter/Transforms/Rotation3DBuilder.hpp>
#include <Mathter/Transforms/ScaleBuilder.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


template <class Mat>
static void VerifyDecomposition(const Mat& A, const Mat& L, const Mat& U) {
	using Scalar = scalar_type_t<Mat>;
	using Real = remove_complex_t<Scalar>;
	constexpr auto tolerance = Real(15) * std::numeric_limits<Real>::epsilon();

	REQUIRE(NormPrecise(ZeroLowerTriangle(L)) < tolerance * NormPrecise(L));
	REQUIRE(NormPrecise(ZeroUpperTriangle(U)) < tolerance * NormPrecise(U));

	REQUIRE(L * U == test_util::Approx(A));
}


template <class Mat, class Perm>
static void VerifyDecomposition(const Mat& A, const Mat& L, const Mat& U, const Perm& P) {
	using Decomp = std::decay_t<decltype(DecomposeLUP(A))>;
	using Scalar = scalar_type_t<Mat>;
	using Real = remove_complex_t<Scalar>;
	constexpr auto tolerance = Real(15) * std::numeric_limits<Real>::epsilon();


	const auto PM = Decomp::ExpandPermutation(P);

	REQUIRE(NormPrecise(ZeroLowerTriangle(L)) < tolerance * NormPrecise(L));
	REQUIRE(NormPrecise(ZeroUpperTriangle(U)) < tolerance * NormPrecise(U));

	REQUIRE(L * U == test_util::Approx(PM * A));
}


TEMPLATE_LIST_TEST_CASE("LU decomposition: real", "[LU]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;

	const Mat m = {
		1.92f, 1.17f, 0.85f,
		0.78f, 0.09f, -1.21f,
		3.98f, 0.07f, -2.92f
	};

	SECTION("LU") {
		const auto [L, U] = DecomposeLU(m);
		VerifyDecomposition(m, L, U);
	}
	SECTION("LUP") {
		const auto [L, U, P] = DecomposeLUP(m);
		REQUIRE(P == decltype(P){ 2, 0, 1 }); // Make sure permutation is not identity.
		VerifyDecomposition(m, L, U, P);
	}
}


TEMPLATE_LIST_TEST_CASE("LU decomposition: complex", "[LU]",
						decltype(MatrixCaseList<ScalarsComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;

	using namespace std::complex_literals;

	const Mat m = {
		1.f + 0.7if, 0.5f + 1.2if, -1.3f - 0.9if,
		2.f - 0.2if, -1.1f - 1.0if, -0.7f + 0.6if,
		-1.f + 0.3if, 0.3f + 1.2if, 0.3f - 0.1if
	};

	SECTION("LU") {
		const auto [L, U] = DecomposeLU(m);
		VerifyDecomposition(m, L, U);
	}
	SECTION("LUP") {
		const auto [L, U, P] = DecomposeLUP(m);
		REQUIRE(P == decltype(P){ 1, 0, 2 }); // Make sure permutation is not identity.
		VerifyDecomposition(m, L, U, P);
	}
}


TEMPLATE_LIST_TEST_CASE("LU decomposition: zero matrix", "[LU]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;
	using T = scalar_type_t<Mat>;

	const Mat m = Zero();

	SECTION("LUP") {
		const auto [L, U, P] = DecomposeLUP(m);
		REQUIRE(NormPrecise(L * U) < 1e-6f);
	}
}


TEMPLATE_LIST_TEST_CASE("LU decomposition: pivoting correctness", "[LU]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	// Small matrices, like 3x3, do not have enough chances for pivoting
	// for things to go awry. This 6x6 should uncover pivoting problems.
	using Mat = typename TestType::template Matrix<6, 6>;

	const Mat m = {
		1.92f, 1.17f, 0.85f, -0.78f, 0.94f, 0.11f,
		0.78f, 0.09f, -1.21f, 1.76f, -0.56f, 1.34f,
		3.98f, 0.07f, -2.92f, -0.01f, 0.44f, -0.23,
		1.45f, -0.67, -0.93f, 0.45f, 0.68f, -0.35f,
		0.62f, 0.93f, 0.65f, 0.28f, 0.34f, 0.62f,
		0.95f, 0.35f, 0.67f, 0.84f, 0.82f, 0.45f
	};

	SECTION("LU") {
		const auto [L, U] = DecomposeLU(m);
		VerifyDecomposition(m, L, U);
	}
	SECTION("LUP") {
		const auto [L, U, P] = DecomposeLUP(m);
		REQUIRE(P == decltype(P){ 2, 0, 3, 1, 4, 5 }); // Make sure permutation is not identity.
		VerifyDecomposition(m, L, U, P);
	}
}


TEMPLATE_LIST_TEST_CASE("LU decomposition: solve system of equations", "[LU]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;
	using Scalar = scalar_type_t<Mat>;
	using Vec = Vector<Scalar, 3, false>;

	const Mat m = {
		1.92f, 1.17f, 0.85f,
		0.78f, 0.09f, -1.21f,
		3.98f, 0.07f, -2.92f
	};
	const Vec b = { 4, 6, 8 };

	SECTION("LU") {
		const auto LU = DecomposeLU(m);
		const auto x = LU.Solve(b);
		REQUIRE(ApplyTransform(m, x) == test_util::Approx(b));
	}
	SECTION("LUP") {
		const auto LUP = DecomposeLUP(m);
		REQUIRE(LUP.P == decltype(LUP.P){ 2, 0, 1 }); // Make sure permutation is not identity.
		const auto x = LUP.Solve(b);
		REQUIRE(ApplyTransform(m, x) == test_util::Approx(b));
	}
}


TEMPLATE_LIST_TEST_CASE("LU decomposition: inverse", "[LU]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;

	const Mat m = {
		1.92f, 1.17f, 0.85f,
		0.78f, 0.09f, -1.21f,
		3.98f, 0.07f, -2.92f
	};

	SECTION("LU") {
		const auto LU = DecomposeLU(m);
		const auto inverse = LU.Inverse();
		REQUIRE(m * inverse == test_util::Approx(Mat(Identity()), 3e-6f));
	}
	SECTION("LUP") {
		const auto LUP = DecomposeLUP(m);
		const auto inverse = LUP.Inverse();
		REQUIRE(m * inverse == test_util::Approx(Mat(Identity())));
	}
}