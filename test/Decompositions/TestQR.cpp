// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../ApplyTransform.hpp"
#include "../Approx.hpp"
#include "../Cases.hpp"
#include "../MatrixUtil.hpp"

#include <Mathter/Decompositions/DecomposeQR.hpp>
#include <Mathter/Matrix/Arithmetic.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


template <class MatA, class MatQ, class MatR>
static void VerifyDecomposition(const MatA& A, const MatQ& Q, const MatR& R) {
	using Scalar = scalar_type_t<MatA>;
	using Real = remove_complex_t<Scalar>;

	constexpr auto tolerance = Real(15) * std::numeric_limits<Real>::epsilon();

	if constexpr (row_count_v<MatQ> == column_count_v<MatQ>) {
		const auto detQ = Determinant(Q);
		REQUIRE(std::abs(detQ) == Catch::Approx(1.0));
	}

	const auto QtQ = ConjTranspose(Q) * Q;
	REQUIRE(QtQ == test_util::Approx(decltype(QtQ)(Identity()), tolerance));
	REQUIRE(Q * R == test_util::Approx(A, tolerance));

	REQUIRE(NormPrecise(ZeroUpperTriangle(R)) < tolerance * NormPrecise(R));
}


TEMPLATE_LIST_TEST_CASE("QR decomposition: square matrix", "[QR]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;
	using T = scalar_type_t<Mat>;

	// This is the example from Wikipedia's page on QR decomposition.
	// https://en.wikipedia.org/wiki/QR_decomposition#Example_2
	const Mat m = {
		T(12), T(-51), T(4),
		T(6), T(167), T(-68),
		T(-4), T(24), T(-41)
	};
	const Mat qExpected = {
		T(-6.0 / 7), T(69.0 / 175), T(-58.0 / 175),
		T(-3.0 / 7), T(-158.0 / 175), T(6.0 / 175),
		T(2.0 / 7), T(-6.0 / 35), T(-33.0 / 35)
	};
	const Mat rExpected = {
		T(-14), T(-21), T(14),
		T(0), T(-175), T(70),
		T(0), T(0), T(35)
	};

	SECTION("QR") {
		const auto [Q, R] = DecomposeQR(m);
		REQUIRE(Q == test_util::Approx(qExpected));
		REQUIRE(R == test_util::Approx(rExpected));
		VerifyDecomposition(m, Q, R);
	}
	SECTION("LQ") {
		const auto mt = FlipLayoutAndOrder(m);
		const auto rtExpected = FlipLayoutAndOrder(rExpected);
		const auto qtExpected = FlipLayoutAndOrder(qExpected);
		const auto [L, Q] = DecomposeLQ(mt);
		REQUIRE(L == test_util::Approx(rtExpected));
		REQUIRE(Q == test_util::Approx(qtExpected));
		REQUIRE(L * Q == test_util::Approx(mt));
	}
}


TEMPLATE_LIST_TEST_CASE("QR decomposition: tall/wide matrix", "[QR]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<4, 3>;

	const Mat m = {
		1, 2, 3,
		0, 1, 2,
		-1, 0, 1,
		-2, -1, 0
	};

	SECTION("QR") {
		const auto [Q, R] = DecomposeQR(m);
		VerifyDecomposition(m, Q, R);
	}
	SECTION("LQ") {
		const auto mt = FlipLayoutAndOrder(m);
		const auto [L, Q] = DecomposeLQ(mt);
		REQUIRE(L * Q == test_util::Approx(mt));
	}
}


TEMPLATE_LIST_TEST_CASE("QR decomposition: zero matrix", "[QR]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;

	const Mat m = Zero();

	const auto [Q, R] = DecomposeQR(m);

	const auto QtQ = ConjTranspose(Q) * Q;
	REQUIRE(QtQ == test_util::Approx(decltype(QtQ)(Identity()), 1e-6f));
	REQUIRE(Max(Abs(R)) < 1e-6f);
}


TEMPLATE_LIST_TEST_CASE("QR decomposition: complex matrix", "[QR]",
						decltype(MatrixCaseList<ScalarsComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;
	using namespace std::complex_literals;

	const Mat m = {
		1.f + 0.7if, 0.5f + 1.2if, -1.3f - 0.9if,
		2.f - 0.2if, -1.1f - 1.0if, -0.7f + 0.6if,
		-1.f + 0.3if, 0.3f + 1.2if, 0.3f - 0.1if
	};

	SECTION("QR") {
		const auto [Q, R] = DecomposeQR(m);
		VerifyDecomposition(m, Q, R);
	}
	SECTION("LQ") {
		const auto mt = FlipLayoutAndOrder(m);
		const auto [L, Q] = DecomposeLQ(mt);
		REQUIRE(L * Q == test_util::Approx(mt));
	}
}


TEMPLATE_LIST_TEST_CASE("QR decomposition: compute real inverse", "[QR]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;
	using namespace std::complex_literals;

	const Mat m = {
		0.89f, 1.94f, -2.33f,
		1.74f, -0.89f, -1.11f,
		2.55f, 1.38f, 0.62f
	};

	SECTION("QR") {
		const auto QR = DecomposeQR(m);
		const auto inverse = QR.Inverse();
		REQUIRE(m * inverse == test_util::Approx(Mat(Identity())));
	}
	SECTION("LQ") {
		const auto mt = FlipLayoutAndOrder(m);
		const auto LQ = DecomposeLQ(mt);
		const auto inverse = LQ.Inverse();
		REQUIRE(mt * inverse == test_util::Approx(std::decay_t<decltype(mt)>(Identity())));
	}
}


TEMPLATE_LIST_TEST_CASE("QR decomposition: compute complex inverse", "[QR]",
						decltype(MatrixCaseList<ScalarsComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;
	using namespace std::complex_literals;

	const Mat m = {
		1.f + 0.7if, 0.5f + 1.2if, -1.3f - 0.9if,
		2.f - 0.2if, -1.1f - 1.0if, -0.7f + 0.6if,
		-1.f + 0.3if, 0.3f + 1.2if, 0.3f - 0.1if
	};

	SECTION("QR") {
		const auto QR = DecomposeQR(m);
		const auto inverse = QR.Inverse();
		REQUIRE(m * inverse == test_util::Approx(Mat(Identity())));
	}
	SECTION("LQ") {
		const auto mt = FlipLayoutAndOrder(m);
		const auto LQ = DecomposeLQ(mt);
		const auto inverse = LQ.Inverse();
		REQUIRE(mt * inverse == test_util::Approx(std::decay_t<decltype(mt)>(Identity())));
	}
}


TEMPLATE_LIST_TEST_CASE("QR decomposition: compute real pseudoinverse", "[QR]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<4, 3>;
	using namespace std::complex_literals;

	const Mat m = {
		0.89f, 1.94f, -2.33f,
		1.74f, -0.89f, -1.11f,
		2.55f, 1.38f, 0.62f,
		0.77f, -0.98f, 0.13f
	};

	SECTION("QR") {
		const auto QR = DecomposeQR(m);
		const auto pseudoinverse = QR.Inverse();
		static_assert(row_count_v<std::decay_t<decltype(pseudoinverse)>> == 3);
		static_assert(column_count_v<std::decay_t<decltype(pseudoinverse)>> == 4);
		// The Moore-Penrose conditions
		REQUIRE(m * pseudoinverse * m == test_util::Approx(m));
		REQUIRE(pseudoinverse * m * pseudoinverse == test_util::Approx(pseudoinverse));
		REQUIRE(ConjTranspose(m * pseudoinverse) == test_util::Approx(m * pseudoinverse));
		REQUIRE(ConjTranspose(pseudoinverse * m) == test_util::Approx(pseudoinverse * m));
	}
	SECTION("LQ") {
		const auto mt = FlipLayoutAndOrder(m);
		const auto LQ = DecomposeLQ(mt);
		const auto pseudoinverse = LQ.Inverse();
		REQUIRE(mt * pseudoinverse * mt == test_util::Approx(mt));
	}
}


TEMPLATE_LIST_TEST_CASE("QR decomposition: compute complex pseudoinverse", "[QR]",
						decltype(MatrixCaseList<ScalarsComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<4, 3>;
	using namespace std::complex_literals;

	const Mat m = {
		1.f + 0.7if, 0.5f + 1.2if, -1.3f - 0.9if,
		2.f - 0.2if, -1.1f - 1.0if, -0.7f + 0.6if,
		-1.f + 0.3if, 0.3f + 1.2if, 0.3f - 0.1if,
		-0.4f + 1.3f, 0.6f - 0.2f, -0.9f + 0.7f
	};

	SECTION("QR") {
		const auto QR = DecomposeQR(m);
		const auto pseudoinverse = QR.Inverse();
		static_assert(row_count_v<std::decay_t<decltype(pseudoinverse)>> == 3);
		static_assert(column_count_v<std::decay_t<decltype(pseudoinverse)>> == 4);
		// The Moore-Penrose conditions
		REQUIRE(m * pseudoinverse * m == test_util::Approx(m));
		REQUIRE(pseudoinverse * m * pseudoinverse == test_util::Approx(pseudoinverse));
		REQUIRE(ConjTranspose(m * pseudoinverse) == test_util::Approx(m * pseudoinverse));
		REQUIRE(ConjTranspose(pseudoinverse * m) == test_util::Approx(pseudoinverse * m));
	}
	SECTION("LQ") {
		const auto mt = FlipLayoutAndOrder(m);
		const auto LQ = DecomposeLQ(mt);
		const auto pseudoinverse = LQ.Inverse();
		REQUIRE(mt * pseudoinverse * mt == test_util::Approx(mt));
	}
}


TEMPLATE_LIST_TEST_CASE("QR decomposition: solve system of equations", "[QR]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersPrecede, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;
	using Vec = Vector<double, 3, false>;
	using namespace std::complex_literals;

	const Mat m = {
		0.89f, 1.94f, -2.33f,
		1.74f, -0.89f, -1.11f,
		2.55f, 1.38f, 0.62f
	};

	const Vec b = { 4, 6, 8 };

	SECTION("QR") {
		const auto QR = DecomposeQR(m);
		const auto x = QR.Solve(b);
		REQUIRE(m * x == test_util::Approx(b, 1e-6f));
	}
	SECTION("LQ") {
		const auto mt = FlipLayoutAndOrder(m);
		const auto LQ = DecomposeLQ(mt);
		const auto x = LQ.Solve(b);
		REQUIRE(x * mt == test_util::Approx(b, 1e-6f));
	}
}


TEMPLATE_LIST_TEST_CASE("QR decomposition: solve multiple systems of equations", "[QR]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersPrecede, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<3, 3>;
	using MatB = typename TestType::template Matrix<3, 2>;
	using namespace std::complex_literals;

	const Mat m = {
		0.89f, 1.94f, -2.33f,
		1.74f, -0.89f, -1.11f,
		2.55f, 1.38f, 0.62f
	};

	const MatB b = {
		4, 5,
		6, 7,
		8, 9
	};

	SECTION("QR") {
		const auto QR = DecomposeQR(m);
		const auto x = QR.Solve(b);
		REQUIRE(m * x.Column(0) == test_util::Approx(b.Column(0), 1e-6f));
		REQUIRE(m * x.Column(1) == test_util::Approx(b.Column(1), 1e-6f));
	}
	SECTION("LQ") {
		const auto mt = FlipLayoutAndOrder(m);
		const auto LQ = DecomposeLQ(mt);
		const auto x = LQ.Solve(FlipLayoutAndOrder(b));
		REQUIRE(x.Row(0) * mt == test_util::Approx(b.Column(0), 1e-6f));
		REQUIRE(x.Row(1) * mt == test_util::Approx(b.Column(1), 1e-6f));
	}
}


TEMPLATE_LIST_TEST_CASE("QR decomposition: solve least squares problem", "[QR]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersPrecede, LayoutsAll, PackingsAll>{})) {
	using Mat = typename TestType::template Matrix<4, 3>;
	using Vec = Vector<double, 4, false>;
	using namespace std::complex_literals;

	const Mat m = {
		0.89f, 1.94f, -2.33f,
		1.74f, -0.89f, -1.11f,
		2.55f, 1.38f, 0.62f,
		0.77f, -0.98f, 0.13f
	};

	const Vec b = { 4, 6, 8, 10 };

	SECTION("QR") {
		const auto QR = DecomposeQR(m);
		const auto x = QR.Solve(b);
		// A reference solution is calculated using the pseudinverse.
		const auto pseudoinverse = QR.Inverse();
		const auto xRef = pseudoinverse * b;
		REQUIRE(x == test_util::Approx(xRef, 1e-6f));
	}
	SECTION("LQ") {
		const auto mt = FlipLayoutAndOrder(m);
		const auto LQ = DecomposeLQ(mt);
		const auto x = LQ.Solve(b);
		const auto pseudoinverse = LQ.Inverse();
		const auto xRef = b * pseudoinverse;
		REQUIRE(x == test_util::Approx(xRef, 1e-6f));
	}
}


TEMPLATE_LIST_TEST_CASE("QR decomposition: QR/LQ selection", "[QR]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	SECTION("Rectangular") {
		// These only need to compile. Behaviour is tested separately.
		using Mat = typename TestType::template Matrix<4, 3>;

		const Mat m = {
			1, 2, 3,
			0, 1, 2,
			-1, 0, 1,
			-2, -1, 0
		};
		const auto mt = FlipLayoutAndOrder(m);

		REQUIRE_NOTHROW(DecomposeQRorLQ(m));
		REQUIRE_NOTHROW(DecomposeQRorLQ(mt));
	}
	SECTION("Square") {
		using Mat = typename TestType::template Matrix<3, 3>;
		using Vec = Vector<scalar_type_t<Mat>, 3, is_packed_v<Mat>>;

		const Mat m = {
			0.89f, 1.94f, -2.33f,
			1.74f, -0.89f, -1.11f,
			2.55f, 1.38f, 0.62f
		};
		const Vec b = { 4, 6, 8 };

		const auto dec = DecomposeQRorLQ(m);
		const auto x = dec.Solve(b);
		REQUIRE(ApplyTransform(m, x) == test_util::Approx(b));
	}
}