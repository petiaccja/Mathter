// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Decompositions/DecomposeSVD.hpp>
#include <Mathter/IoStream.hpp>
#include <Mathter/Transforms/Rotation2DBuilder.hpp>
#include <Mathter/Transforms/ScaleBuilder.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <cmath>
#include <complex>


using namespace mathter;
using namespace test_util;


template <class T, class Out = T>
auto ExpandRQ2x2(const mathter::impl::DecompositionRQ2x2<T>& rq) {
	using M22 = Matrix<Out, 2, 2>;
	return std::tuple{
		M22{ rq.r11, rq.r12, 0.0f, rq.r22 },
		M22{ rq.cq, -rq.sq, rq.sq, rq.cq },
	};
}


template <class T, class Out = T>
auto ExpandSVD2x2(const mathter::impl::DecompositionSVD2x2<T>& svd) {
	using M22 = Matrix<Out, 2, 2>;
	return std::tuple{
		M22{ svd.cu, -svd.su, svd.su, svd.cu },
		M22{ svd.s11, 0.0f, 0.0f, svd.s22 },
		M22{ svd.cv, -svd.sv, svd.sv, svd.cv },
	};
}


template <class T>
auto ExpandDiagSymm2x2(std::tuple<T, T> v) {
	using M22 = Matrix<T, 2, 2>;
	const auto [cv, sv] = v;
	return M22{ cv, -sv, sv, cv };
}


template <class T>
auto ExpandDiagSymm2x2(T a11, T aoff, T a22) {
	using M22 = Matrix<T, 2, 2>;
	return M22{ a11, aoff, aoff, a22 };
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed, class Out = T>
auto ExpandSVD(const DecompositionSVD<T, Rows, Columns, Order, Layout, Packed>& svd) {
	using MPP = Matrix<T, svd.PDim, svd.PDim, Order, Layout, Packed>;
	return std::tuple{ svd.U, MPP(Scale(svd.S)), svd.V };
}


static constexpr std::array hammerExponents = {
	-10.f,
	-4.0f,
	-2.0f,
	-1.0f,
	-0.5f,
	0.0f,
	0.5f,
	0.5f,
	1.0f,
	1.5f,
	2.0f,
	4.0f,
	8.0f,
	17.0f,
	18.f,
	19.f,
	24.f,
	32.f,
};


static constexpr std::array hammerAngleU = {
	0.0,
	0.0000001,
	0.000001,
	0.00001,
	0.0001,
	0.001,
	0.01,
	1.5707963267948966,
	3.1315926535897933,
	3.1405926535897932,
	3.141492653589793,
	3.141582653589793,
	3.141591653589793,
	3.1415925535897933,
	3.141592653589793,
};

static constexpr std::array hammerMagS22 = {
	0.0,
	1e-22,
	1e-16,
	0.0000001,
	0.000001,
	0.00001,
	0.0001,
	0.001,
	0.999,
	0.9999,
	0.99999,
	0.999999,
	0.9999999,
	1.0,
};

static constexpr std::array hammerAngleV = hammerAngleU;


template <class T>
auto GenFromExp(T exponent) {
	const auto value = std::pow(T(10), -std::abs(exponent));
	return value * (exponent < 0 ? -T(1) : T(1));
}


template <class MatA, class MatU, class MatS, class MatV>
void VerifySVD(MatA A, MatU U, MatS S, MatV V, scalar_type_t<MatA> tolerance) {
	const auto USV = U * S * V;
	const auto UTU = Transpose(U) * U;
	const auto VTV = Transpose(V) * V;

	INFO("A =   " << A);
	INFO("U =   " << U);
	INFO("S =   " << S);
	INFO("V =   " << V);

	REQUIRE_THAT(Determinant(U), Catch::Matchers::WithinAbs(1.0f, tolerance) || Catch::Matchers::WithinAbs(-1.0f, tolerance));
	REQUIRE_THAT(Determinant(V), Catch::Matchers::WithinAbs(1.0f, tolerance) || Catch::Matchers::WithinAbs(-1.0f, tolerance));
	REQUIRE(UTU == test_util::Approx(decltype(UTU)(Identity()), tolerance));
	REQUIRE(VTV == test_util::Approx(decltype(VTV)(Identity()), tolerance));
	REQUIRE(A == test_util::Approx(USV, tolerance));
}


TEST_CASE("Matrix - RQ 2x2", "[Matrix]") {
	using M22 = Matrix<float, 2, 2>;

	SECTION("Zero") {
		const M22 A = {
			0.0f, 0.0f,
			0.0f, 0.0f
		};
		const auto rq = mathter::impl::DecomposeRQ2x2(A);
		REQUIRE(rq.r11 == 0.0f);
		REQUIRE(rq.r12 == 0.0f);
		REQUIRE(rq.r22 == 0.0f);
		REQUIRE(rq.cq == 1.0f);
		REQUIRE(rq.sq == 0.0f);
	}
	SECTION("Upper triangular") {
		const M22 A = {
			1.0f, 2.0f,
			0.0f, 3.0f
		};
		const auto rq = mathter::impl::DecomposeRQ2x2(A);
		REQUIRE(rq.r11 == 1.0f);
		REQUIRE(rq.r12 == 2.0f);
		REQUIRE(rq.r22 == 3.0f);
		REQUIRE(rq.cq == 1.0f);
		REQUIRE(rq.sq == 0.0f);
	}
	SECTION("Lower single") {
		const M22 A = {
			0.0f, 0.0f,
			1.0f, 0.0f
		};
		const auto rq = mathter::impl::DecomposeRQ2x2(A);
		const auto [R, Q] = ExpandRQ2x2(rq);
		REQUIRE(Determinant(Q) == Catch::Approx(1));
		REQUIRE(NormPrecise(A - R * Q) < 1e-6f * NormPrecise(A));
	}
	SECTION("Upper single") {
		const M22 A = {
			0.0f, 1.0f,
			0.0f, 0.0f
		};
		const auto rq = mathter::impl::DecomposeRQ2x2(A);
		const auto [R, Q] = ExpandRQ2x2(rq);
		REQUIRE(Determinant(Q) == Catch::Approx(1));
		REQUIRE(NormPrecise(A - R * Q) < 1e-6f * NormPrecise(A));
	}
	SECTION("Tiny lower singular") {
		const M22 A = {
			1.0f, 0.0f,
			1e-25f, 0.0f
		};
		const auto rq = mathter::impl::DecomposeRQ2x2(A);
		const auto [R, Q] = ExpandRQ2x2(rq);
		REQUIRE(Determinant(Q) == Catch::Approx(1));
		REQUIRE(NormPrecise(A - R * Q) < 1e-6f * NormPrecise(A));
	}
	SECTION("Tiny lower diminished") {
		const M22 A = {
			1.0f, 0.0f,
			1e-25f, 1.0f
		};
		const auto rq = mathter::impl::DecomposeRQ2x2(A);
		const auto [R, Q] = ExpandRQ2x2(rq);
		REQUIRE(Determinant(Q) == Catch::Approx(1));
		REQUIRE(NormPrecise(A - R * Q) < 1e-6f * NormPrecise(A));
	}
	SECTION("r11 cancellation") {
		const M22 A = {
			0.999'999f, 0.999'999f,
			1.000'001f, 1.000'002f
		};
		const auto rq = mathter::impl::DecomposeRQ2x2(A);
		const auto [R, Q] = ExpandRQ2x2(rq);
		REQUIRE(Determinant(Q) == Catch::Approx(1));
		REQUIRE(NormPrecise(A - R * Q) < 1e-6f * NormPrecise(A));
	}
}


TEST_CASE("Matrix - RQ 2x2 hammer", "[Matrix]") {
	using M22 = Matrix<float, 2, 2>;

	for (auto a11Exp : hammerExponents) {
		for (auto a12Exp : hammerExponents) {
			for (auto a21Exp : hammerExponents) {
				for (auto a22Exp : hammerExponents) {
					const M22 A = {
						GenFromExp<float>(a11Exp), GenFromExp<float>(a12Exp),
						GenFromExp<float>(a21Exp), GenFromExp<float>(a22Exp)
					};
					const auto rq = mathter::impl::DecomposeRQ2x2(A);
					const auto [R, Q] = ExpandRQ2x2(rq);
					INFO("A = " << A);
					INFO("R = " << R);
					INFO("Q = " << Q);
					REQUIRE(Determinant(Q) == Catch::Approx(1));
					REQUIRE(NormPrecise(A - R * Q) < 1e-6f * NormPrecise(A));
				}
			}
		}
	}
}


TEST_CASE("Matrix - SVD 2x2 (2-sided core)", "[Matrix]") {
	using M22 = Matrix<float, 2, 2>;

	SECTION("Zero") {
		const M22 A = {
			0.0f, 0.0f,
			0.0f, 0.0f
		};
		const auto svd = mathter::impl::DecomposeSVD2x2(A);
		REQUIRE(svd.cu == 1.0f);
		REQUIRE(svd.su == 0.0f);
		REQUIRE(svd.s11 == 0.0f);
		REQUIRE(svd.s22 == 0.0f);
		REQUIRE(svd.cv == 1.0f);
		REQUIRE(svd.sv == 0.0f);
	}
	SECTION("Top left") {
		const M22 A = {
			1.0f, 0.0f,
			0.0f, 0.0f
		};
		const auto svd = mathter::impl::DecomposeSVD2x2(A);
		const auto [U, S, V] = ExpandSVD2x2(svd);
		VerifySVD(A, U, S, V, 1e-6f);
	}
	SECTION("Top right") {
		const M22 A = {
			0.0f, 1.0f,
			0.0f, 0.0f
		};
		const auto svd = mathter::impl::DecomposeSVD2x2(A);
		const auto [U, S, V] = ExpandSVD2x2(svd);
		VerifySVD(A, U, S, V, 1e-6f);
	}
	SECTION("Bottom right") {
		const M22 A = {
			0.0f, 0.0f,
			0.0f, 1.0f
		};
		const auto svd = mathter::impl::DecomposeSVD2x2(A);
		const auto [U, S, V] = ExpandSVD2x2(svd);
		VerifySVD(A, U, S, V, 1e-6f);
	}
	SECTION("Cancellation within z") {
		const M22 A = {
			1.0f, 0.7071067f,
			0.0f, 0.7071067f
		};
		const auto svd = mathter::impl::DecomposeSVD2x2(A);
		const auto [U, S, V] = ExpandSVD2x2(svd);
		VerifySVD(A, U, S, V, 1e-6f);
	}
	SECTION("Positive z") {
		const M22 A = {
			1.0f, 0.05f,
			0.0f, 0.05f
		};
		const auto svd = mathter::impl::DecomposeSVD2x2(A);
		const auto [U, S, V] = ExpandSVD2x2(svd);
		VerifySVD(A, U, S, V, 1e-6f);
	}
	SECTION("Negative z") {
		const M22 A = {
			0.05f, 1.0f,
			0.0f, 1.0f
		};
		const auto svd = mathter::impl::DecomposeSVD2x2(A);
		const auto [U, S, V] = ExpandSVD2x2(svd);
		VerifySVD(A, U, S, V, 1e-6f);
	}
	SECTION("RQ preconditioner") {
		const M22 A = {
			0.8f, 1.3f,
			0.5f, 0.9f
		};
		const auto svd = mathter::impl::DecomposeSVD2x2(A);
		const auto [U, S, V] = ExpandSVD2x2(svd);
		VerifySVD(A, U, S, V, 1e-6f);
	}
	SECTION("Diagonal degeneration") {
		const M22 A = {
			1.0f, 1e-22f,
			0.0f, 1.0f
		};
		const auto svd = mathter::impl::DecomposeSVD2x2(A);
		const auto [U, S, V] = ExpandSVD2x2(svd);
		VerifySVD(A, U, S, V, 1e-6f);
	}
	SECTION("Diagonal singularity") {
		const M22 A = {
			1.0f, 0.0f,
			0.0f, 1.0f
		};
		const auto svd = mathter::impl::DecomposeSVD2x2(A);
		const auto [U, S, V] = ExpandSVD2x2(svd);
		VerifySVD(A, U, S, V, 1e-6f);
	}
}


TEST_CASE("Matrix - SVD 2x2 general hammer (2-sided core)", "[Matrix]") {
	using M22 = Matrix<float, 2, 2>;

	for (auto a11Exp : hammerExponents) {
		for (auto a12Exp : hammerExponents) {
			for (auto a22Exp : hammerExponents) {
				const M22 A = {
					GenFromExp<float>(a11Exp), GenFromExp<float>(a12Exp),
					0.0f, GenFromExp<float>(a22Exp)
				};
				const auto svd = mathter::impl::DecomposeSVD2x2(A);
				const auto [U, S, V] = ExpandSVD2x2(svd);
				VerifySVD(A, U, S, V, 1e-6f);
			}
		}
	}
}


TEST_CASE("Matrix - SVD 2x2 edge case hammer (2-sided core)", "[Matrix]") {
	using M22 = Matrix<float, 2, 2>;
	using M22d = Matrix<float, 2, 2>;

	for (auto u : hammerAngleU) {
		for (auto s22 : hammerMagS22) {
			for (auto v : hammerAngleV) {
				const M22d Uref = Rotation(u);
				const M22d Sref = Scale(1.0f, s22);
				const M22d Vref = Rotation(v);
				const auto Aref = Uref * Sref * Vref;
				const auto A = M22(Aref);

				const auto svd = mathter::impl::DecomposeSVD2x2(A);
				const auto [U, S, V] = ExpandSVD2x2(svd);
				VerifySVD(A, U, S, V, 1e-6f);
			}
		}
	}
}


TEST_CASE("Matrix - Diagonalize symmetric 2x2", "[Matrix]") {
	SECTION("General") {
		const auto a11 = 1.0f;
		const auto aoff = 0.7f;
		const auto a22 = 0.9f;

		const auto A = ExpandDiagSymm2x2(a11, aoff, a22);
		const auto V = ExpandDiagSymm2x2(mathter::impl::DiagonalizeSymmetric2x2(a11, aoff, a22));
		const auto D = Transpose(V) * A * V;
		INFO("A = " << A);
		INFO("V = " << V);
		INFO("D = " << D);
		REQUIRE(D(0, 1) < 1e-6f * NormPrecise(A));
		REQUIRE(D(1, 0) < 1e-6f * NormPrecise(A));
		REQUIRE(Determinant(D) - Determinant(A) < 1e-6f * NormPrecise(A));
	}
	SECTION("Cancellation z") {
		const auto a11 = 1.0f;
		const auto aoff = 0.0f;
		const auto a22 = 0.999999f;

		const auto A = ExpandDiagSymm2x2(a11, aoff, a22);
		const auto V = ExpandDiagSymm2x2(mathter::impl::DiagonalizeSymmetric2x2(a11, aoff, a22));
		const auto D = Transpose(V) * A * V;
		INFO("A = " << A);
		INFO("V = " << V);
		INFO("D = " << D);
		REQUIRE(D(0, 1) < 1e-6f * NormPrecise(A));
		REQUIRE(D(1, 0) < 1e-6f * NormPrecise(A));
		REQUIRE(Determinant(D) - Determinant(A) < 1e-6f * NormPrecise(A));
	}
	SECTION("Singularity z") {
		const auto a11 = 1.0f;
		const auto aoff = 0.0f;
		const auto a22 = 1.0f;

		const auto A = ExpandDiagSymm2x2(a11, aoff, a22);
		const auto V = ExpandDiagSymm2x2(mathter::impl::DiagonalizeSymmetric2x2(a11, aoff, a22));
		const auto D = Transpose(V) * A * V;
		INFO("A = " << A);
		INFO("V = " << V);
		INFO("D = " << D);
		REQUIRE(D(0, 1) < 1e-6f * NormPrecise(A));
		REQUIRE(D(1, 0) < 1e-6f * NormPrecise(A));
		REQUIRE(Determinant(D) - Determinant(A) < 1e-6f * NormPrecise(A));
	}
	SECTION("Near singularity p1") {
		const auto a11 = 1.0f;
		const auto aoff = 1e-21f;
		const auto a22 = 1.0f;

		const auto A = ExpandDiagSymm2x2(a11, aoff, a22);
		const auto V = ExpandDiagSymm2x2(mathter::impl::DiagonalizeSymmetric2x2(a11, aoff, a22));
		const auto D = Transpose(V) * A * V;
		INFO("A = " << A);
		INFO("V = " << V);
		INFO("D = " << D);
		REQUIRE(D(0, 1) < 1e-6f * NormPrecise(A));
		REQUIRE(D(1, 0) < 1e-6f * NormPrecise(A));
		REQUIRE(Determinant(D) - Determinant(A) < 1e-6f * NormPrecise(A));
	}
	SECTION("Tiny") {
		const auto a11 = 1e-21f;
		const auto aoff = 1e-25f;
		const auto a22 = 1e-21f;

		const auto A = ExpandDiagSymm2x2(a11, aoff, a22);
		const auto V = ExpandDiagSymm2x2(mathter::impl::DiagonalizeSymmetric2x2(a11, aoff, a22));
		const auto D = Transpose(V) * A * V;
		INFO("A = " << A);
		INFO("V = " << V);
		INFO("D = " << D);
		REQUIRE(D(0, 1) < 1e-6f * NormPrecise(A));
		REQUIRE(D(1, 0) < 1e-6f * NormPrecise(A));
		REQUIRE(Determinant(D) - Determinant(A) < 1e-6f * NormPrecise(A));
	}
}


TEST_CASE("Matrix - SVD 2x2 DEBUG (1-sided)", "[Matrix]") {
	using M22 = Matrix<float, 2, 2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR>;

	const M22 A = {
		-1e-10, -0.0001, 0, -1e-10
	};
	const auto [Uref, Sref, Vref] = ExpandSVD2x2(mathter::impl::DecomposeSVD2x2(A));
	INFO("U_ref = " << Uref);
	INFO("S_ref = " << Sref);
	INFO("V_ref = " << Vref);
	const auto svd = mathter::impl::DecomposeSVDJacobiOneSided(A);
	const auto svdRef = mathter::impl::DecomposeSVDJacobiTwoSided(A);
	const auto [U, S, V] = ExpandSVD(svd);
	VerifySVD(A, U, S, V, 1e-6f);
}


TEST_CASE("Matrix - SVD 2x2 general hammer (1-sided)", "[Matrix]") {
	using M22 = Matrix<float, 2, 2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR>;

	for (auto a11Exp : hammerExponents) {
		for (auto a12Exp : hammerExponents) {
			for (auto a22Exp : hammerExponents) {
				const M22 A = {
					GenFromExp<float>(a11Exp), GenFromExp<float>(a12Exp),
					0.0f, GenFromExp<float>(a22Exp)
				};
				const auto [U, S, V] = mathter::impl::DecomposeSVDJacobiOneSided(A);
				VerifySVD(A, U, M22(Scale(S)), V, 1e-6f);
			}
		}
	}
}


TEST_CASE("Matrix - SVD 2x2 case hammer (1-sided)", "[Matrix]") {
	using M22 = Matrix<float, 2, 2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR>;
	using M22d = Matrix<float, 2, 2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR>;

	for (auto u : hammerAngleU) {
		for (auto s22 : hammerMagS22) {
			for (auto v : hammerAngleV) {
				const M22d Uref = Rotation(u);
				const M22d Sref = Scale(1.0f, s22);
				const M22d Vref = Rotation(v);
				const auto Aref = Uref * Sref * Vref;
				const auto A = M22(Aref);

				const auto [U, S, V] = mathter::impl::DecomposeSVDJacobiOneSided(A);
				VerifySVD(A, U, M22(Scale(S)), V, 1e-6f);
			}
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Matrix - SVD square matrix", "[Matrix]",
						decltype(MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using Mat = Matrix<float, 4, 4>;


	const Mat A = {
		1, 3, 4, -1,
		3, -1, -3, 2,
		4, -2, -1, 0,
		2, -4, 1, 3
	};

	SECTION("1-sided algorithm") {
		const auto [U, S, V] = mathter::DecomposeSVD(A, SVDAlgorithmOneSided);
		VerifySVD(A, U, Mat(Scale(S)), V, 1e-6f);
	}
	SECTION("2-sided algorithm") {
		const auto [U, S, V] = mathter::DecomposeSVD(A, SVDAlgorithmTwoSided);
		VerifySVD(A, U, Mat(Scale(S)), V, 1e-6f);
	}
}