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



TEST_CASE_VARIANT("Matrix - LU decomposition", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> A = {
			3, -0.1, -0.2,
			0.1, 7, -0.3,
			0.3, -0.2, 10
		};

		auto [L, U] = DecomposeLU(A);

		for (int i = 0; i < A.RowCount(); ++i) {
			for (int j = 0; j < i - 1; ++j) {
				REQUIRE(U(i, j) == Approx(0.0));
				REQUIRE(L(j, i) == Approx(0.0));
			}
		}

		auto Mprod = L * U;
		REQUIRE(ApproxVec(A) == Mprod);
	}
}


TEST_CASE_VARIANT("Matrix - LU solve", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> A = {
			3, -0.1f, -0.2f,
			0.1f, 7, -0.3f,
			0.3f, -0.2f, 10
		};
		Vector<Type, 3, Packed> b = { 7.85, -19.3, 71.4 };
		Vector<Type, 3, Packed> x;
		Vector<Type, 3, Packed> xexp = { 3, -2.5, 7 };

		x = DecomposeLU(A).Solve(b);
		REQUIRE(ApproxVec(x) == xexp);
	}
}


TEST_CASE_VARIANT("Matrix - LUP decomposition", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<3, 3> A = {
			3, -0.1f, -0.2f,
			0.3f, -0.2f, 10,
			0.1f, 7, -0.3f
		};

		auto [L, U, P] = DecomposeLUP(A);

		for (int i = 0; i < A.RowCount(); ++i) {
			for (int j = 0; j < i - 1; ++j) {
				REQUIRE(U(i, j) == Approx(0.0f));
				REQUIRE(L(j, i) == Approx(0.0f));
			}
		}

		MatrixT<3, 3> Pm = Zero();
		for (int i : P) {
			Pm(i, P(i)) = 1.0f;
		}

		auto Mprod = Transpose(Pm) * L * U;
		REQUIRE(ApproxVec(A) == Mprod);
	}
}


TEST_CASE_VARIANT("Matrix - LUP solve", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		MatrixT<4, 4> A = {
			1, 3, 4, 6,
			3, 6, 2, 6,
			9, 2, 6, 7,
			6, 2, 7, 5
		};
		Vector<Type, 4, Packed> b = { 3, 4, 2, 8 };
		Vector<Type, 4, Packed> x;
		Vector<Type, 4, Packed> xexp = { -94.f / 497, 895.f / 497, 1000.f / 497, -850.f / 497 };

		x = DecomposeLUP(A).Solve(b);

		REQUIRE(ApproxVec(x) == xexp);
	}
}


TEST_CASE("Matrix - LUP decomposition singular", "[Matrix]") {
	Matrix<float, 3, 3> A = {
		1, 0, 0,
		0, 0, 1,
		0, -1, 0
	};

	auto [L, U, P] = DecomposeLUP(A);

	for (int i = 0; i < A.RowCount(); ++i) {
		for (int j = 0; j < i - 1; ++j) {
			REQUIRE(U(i, j) == Approx(0.0f));
			REQUIRE(L(j, i) == Approx(0.0f));
		}
	}

	Matrix<float, 3, 3> Pm = Zero();
	for (int i : P) {
		Pm(i, P(i)) = 1.0f;
	}

	auto Mprod = Transpose(Pm) * L * U;
	REQUIRE(ApproxVec(A) == Mprod);
}


TEST_CASE_VARIANT("Matrix - QR decomposition", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		// example from wikipedia SVD article
		MatrixT<5, 4> A1 = Transpose(MatrixT<4, 5>{
			1, 0, 0, 1,
			2, 0, 0, 3,
			0, 0, 0, 0,
			0, 0, 0, 0,
			2, 0, 0, 0 });
		auto [Q1, R1] = DecomposeQR(A1);
		MatrixT<5, 4> A1assembled = Q1 * R1;
		REQUIRE(ApproxVec(A1assembled) == A1);


		// the same matrix as the LU
		MatrixT<3, 3> A2 = {
			3, -0.1f, -0.2f,
			0.1f, 7, -0.3f,
			0.3f, -0.2f, 10
		};

		auto [Q2, R2] = DecomposeQR(A2);

		MatrixT<3, 3> A2assembled = Q2 * R2;
		REQUIRE(ApproxVec(A2assembled) == A2);
	}
}


TEST_CASE_VARIANT("Matrix - SVD", "[Matrix]", TypesFloating, OrdersFollow, LayoutsAll, PackedAll) {
	SECTION(SECTIONNAME) {
		// example from wikipedia SVD article
		MatrixT<5, 4> A1 = Transpose(MatrixT<4, 5>{
			1, 0, 0, 1, 2,
			0, 0, 3, 0, 0,
			0, 0, 0, 0, 0,
			0, 2, 0, 0, 0 });

		auto [U1, S1, V1] = DecomposeSVD(A1);
		auto A1assembled = U1 * S1 * V1;
		REQUIRE(ApproxVec(A1) == A1assembled);

		auto [U1T, S1T, V1T] = DecomposeSVD(Transpose(A1));
		auto A1Tassembled = U1T * S1T * V1T;
		REQUIRE(ApproxVec(A1Tassembled) == Transpose(A1));

		// The same matrix as the LU
		MatrixT<3, 3> A2 = {
			3, -0.1f, -0.2f,
			0.1f, 7, -0.3f,
			0.3f, -0.2f, 10
		};

		auto [U2, S2, V2] = DecomposeSVD(A2);
		auto A2assembled = U2 * S2 * V2;
		REQUIRE(ApproxVec(A2assembled) == A2);
	}
}


TEST_CASE("Matrix - SVD Identity", "[Matrix]") {
	Matrix<float, 2, 2> m = Identity();

	auto svd = DecomposeSVD(m);
	REQUIRE(svd.U == ApproxVec(Matrix<float, 2, 2>(Identity())));
	REQUIRE(svd.S == ApproxVec(Matrix<float, 2, 2>(Identity())));
	REQUIRE(svd.V == ApproxVec(Matrix<float, 2, 2>(Identity())));

	Matrix<float, 4, 4> m4 = Identity();
	auto svd4 = DecomposeSVD(m4);
	REQUIRE(svd4.U == ApproxVec(Matrix<float, 4, 4>(Identity())));
	REQUIRE(svd4.S == ApproxVec(Matrix<float, 4, 4>(Identity())));
	REQUIRE(svd4.V == ApproxVec(Matrix<float, 4, 4>(Identity())));
}