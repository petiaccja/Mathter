#pragma once

#include "MatrixModule.hpp"

namespace mathter {


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
class MatrixLU {
	using MatrixT = Matrix<T, Rows, Columns, Order, Layout, Packed>;
public:
	friend MatrixT;
	using Inherit = impl::Empty<MatrixLU>;
};

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
class MatrixLU<T, Dim, Dim, Order, Layout, Packed> {
	using MatrixT = Matrix<T, Dim, Dim, Order, Layout, Packed>;
	MatrixT& self() { return *static_cast<MatrixT*>(this); }
	const MatrixT& self() const { return *static_cast<const MatrixT*>(this); }
public:
	void DecomposeLU(MatrixT& L, MatrixT& U) const;
	mathter::DecompositionLU<T, Dim, Order, Layout, Packed> DecompositionLU() const {
		return mathter::DecompositionLU<T, Dim, Order, Layout, Packed>(self());
	}
public:
	friend MatrixT;
	using Inherit = MatrixLU;
};



// From: https://www.gamedev.net/resources/_/technical/math-and-physics/matrix-inversion-using-lu-decomposition-r3637
template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
void MatrixLU<T, Dim, Dim, Order, Layout, Packed>::DecomposeLU(MatrixT& L, MatrixT& U) const {
	const auto& A = self();
	constexpr int n = Dim;

	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			L(i, j) = 0;
		}
		for (int j = 0; j <= i; ++j) {
			U(i, j) = i == j;
		}
	}

	// Crout's algorithm
	for (int i = 0; i < n; ++i) {
		L(i, 0) = A(i, 0);
	}
	for (int j = 1; j < n; ++j) {
		U(0, j) = A(0, j) / L(0, 0);
	}

	for (int j = 1; j < n-1; ++j) {
		for (int i = j; i < n; ++i) {
			float Lij;
			Lij = A(i, j);
			for (int k = 0; k <= j - 1; ++k) {
				Lij -= L(i, k)*U(k, j);
			}
			L(i, j) = Lij;
		}
		for (int k = j; k < n; ++k) {
			float Ujk;
			Ujk = A(j, k);
			for (int i = 0; i <= j - 1; ++i) {
				Ujk -= L(j, i)*U(i, k);
			}
			Ujk /= L(j, j);
			U(j, k) = Ujk;
		}
	}

	L(n - 1, n - 1) = A(n - 1, n - 1);
	for (int k = 0; k < n - 1; ++k) {
		L(n - 1, n - 1) -= L(n - 1, k)*U(k, n - 1);
	}
}




}